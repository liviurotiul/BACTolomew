"""AMRFinderPlus (--plus, and --organism when the species has one) on the Bakta outputs.

AMRFinderPlus reports acquired genes (AMR, STRESS and VIRULENCE elements) and, for the
organisms it knows, point mutations in chromosomal genes (gyrA_S81L, parC_S84L, ...). The
acquired genes become the run's gene evidence (pipeline.BaktaEvidence); the POINT rows (and
POINT_DISRUPT, a disrupted gene) are read from here onto the Bakta CDS that carries them
(ResistanceMutations / ResistanceMutationClasses). A species without an AMRFinderPlus
organism still gets the acquired genes, without point mutation calls.
"""

from __future__ import annotations

import csv
import json
import logging
import os
import re
import subprocess
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import pandas as pd

from pipeline.RunSteps import ProgressCallback, report_progress
from pipeline.ToolEnv import configured_amrfinder_organism, tool_cmd

_LOG = logging.getLogger(__name__)

AMRFINDER_DIR = "amrfinder"
AMRFINDER_SUMMARY = "summary.json"
AMRFINDER_DB_SUBDIR = Path("amrfinderplus-db") / "latest"
AMRFINDER_THREADS_PER_GENOME = 4
POINT_SUBTYPES = ("POINT", "POINT_DISRUPT")
# A POINT_DISRUPT call (a gene disrupted by a stop or frameshift) keeps its symbol with this
# suffix, so the GUI says the gene is disrupted rather than point mutated.
DISRUPTION_SUFFIX = " (gene disrupted)"
MUTATION_COLUMNS = ("ResistanceMutations", "ResistanceMutationClasses")

# mlst (PubMLST) scheme names -> AMRFinderPlus taxgroup. Schemes shared by several taxgroups
# (e.g. "neisseria") are left out: the organism must be unambiguous.
MLST_SCHEME_ORGANISMS: dict[str, str] = {
    "abaumannii": "Acinetobacter_baumannii",
    "abaumannii_2": "Acinetobacter_baumannii",
    "bcc": "Burkholderia_cepacia",
    "bpseudomallei": "Burkholderia_pseudomallei",
    "campylobacter": "Campylobacter",
    "cfreundii": "Citrobacter_freundii",
    "cdifficile": "Clostridioides_difficile",
    "cdiphtheriae": "Corynebacterium_diphtheriae",
    "ecloacae": "Enterobacter_cloacae",
    "efaecalis": "Enterococcus_faecalis",
    "efaecium": "Enterococcus_faecium",
    "ecoli": "Escherichia",
    "ecoli_achtman_4": "Escherichia",
    "hinfluenzae": "Haemophilus_influenzae",
    "koxytoca": "Klebsiella_oxytoca",
    "klebsiella": "Klebsiella_pneumoniae",
    "kpneumoniae": "Klebsiella_pneumoniae",
    "paeruginosa": "Pseudomonas_aeruginosa",
    "senterica": "Salmonella",
    "senterica_achtman_2": "Salmonella",
    "saureus": "Staphylococcus_aureus",
    "spseudintermedius": "Staphylococcus_pseudintermedius",
    "sagalactiae": "Streptococcus_agalactiae",
    "spneumoniae": "Streptococcus_pneumoniae",
    "spyogenes": "Streptococcus_pyogenes",
    "vcholerae": "Vibrio_cholerae",
    "vcholerae_2": "Vibrio_cholerae",
    "vparahaemolyticus": "Vibrio_parahaemolyticus",
    "vvulnificus": "Vibrio_vulnificus",
}


@dataclass(frozen=True)
class AmrFinderResult:
    output_dir: Path
    summary_json: Path
    organism: str | None
    point_rows: int
    failed: tuple[str, ...] = ()


def amrfinder_db_dir(bakta_db: str | os.PathLike | None) -> Path | None:
    """The AMRFinderPlus DB that ships inside the Bakta DB (resolved, so the version shows)."""
    if not bakta_db:
        return None
    path = Path(bakta_db).expanduser() / AMRFINDER_DB_SUBDIR
    return path.resolve() if path.is_dir() else None


def amrfinder_db_version(db_dir: Path) -> str:
    try:
        return (db_dir / "version.txt").read_text(encoding="utf-8").strip()
    except OSError:
        return db_dir.name


def load_taxgroups(db_dir: Path) -> list[str]:
    out: list[str] = []
    try:
        lines = (db_dir / "taxgroup.tsv").read_text(encoding="utf-8").splitlines()
    except OSError:
        return out
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        out.append(line.split("\t", 1)[0].strip())
    return out


def _match_taxgroup(name: str, taxgroups: Sequence[str]) -> str | None:
    """A species or genus name ('Acinetobacter baumannii', 'acinetobacter_baumannii') -> taxgroup."""
    key = re.sub(r"[\s_]+", "_", str(name or "").strip()).lower()
    if not key:
        return None
    by_key = {group.lower(): group for group in taxgroups}
    if key in by_key:
        return by_key[key]
    genus = key.split("_", 1)[0]
    return by_key.get(genus)


def resolve_amrfinder_organism(
    output_path: Path,
    db_dir: Path,
    *,
    schema_dir: str | None = None,
) -> tuple[str | None, str]:
    """The run's AMRFinderPlus --organism and how it was chosen.

    Order: the config override, the majority MLST scheme, then the cgMLST schema's species.
    Returns (None, reason) when nothing maps to a taxgroup the DB has mutations for.
    """
    taxgroups = load_taxgroups(db_dir)
    if not taxgroups:
        return None, f"AMRFinderPlus DB has no taxgroup.tsv: {db_dir}"

    configured = configured_amrfinder_organism()
    if configured:
        match = _match_taxgroup(configured, taxgroups)
        if match:
            return match, f"configured organism '{configured}'"
        return None, f"configured organism '{configured}' is not an AMRFinderPlus taxgroup"

    mlst_json = Path(output_path) / "mlst" / "results.json"
    try:
        rows = json.loads(mlst_json.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        rows = []
    schemes = Counter(
        str(row.get("scheme") or "").strip().lower()
        for row in rows
        if isinstance(row, Mapping) and str(row.get("scheme") or "").strip() not in {"", "-"}
    )
    if schemes:
        (top_scheme, top_count), *rest = schemes.most_common()
        if rest and rest[0][1] == top_count:
            return None, "MLST schemes are tied, so the organism is ambiguous: " + ", ".join(
                f"{name} ({count})" for name, count in schemes.most_common()
            )
        organism = MLST_SCHEME_ORGANISMS.get(top_scheme)
        if organism and organism in taxgroups:
            return organism, f"MLST scheme '{top_scheme}' ({top_count}/{sum(schemes.values())} genomes)"

    if schema_dir:
        match = _match_taxgroup(Path(schema_dir).name, taxgroups)
        if match:
            return match, f"cgMLST schema '{Path(schema_dir).name}'"

    label = ", ".join(schemes) if schemes else "none"
    return None, f"no AMRFinderPlus organism for this species (MLST schemes: {label})"


def _bakta_inputs(output_path: Path, stem: str) -> tuple[Path, Path, Path]:
    sample_dir = Path(output_path) / "bakta" / stem
    return (
        sample_dir / f"{stem}.faa",
        sample_dir / f"{stem}.fna",
        sample_dir / f"{stem}.gff3",
    )


def amrfinder_tsv_path(output_path: Path, stem: str) -> Path:
    return Path(output_path) / AMRFINDER_DIR / f"{stem}.tsv"


def _run_one(stem: str, output_path: Path, db_dir: Path, organism: str | None, threads: int) -> Path:
    faa, fna, gff = _bakta_inputs(output_path, stem)
    missing = [str(p) for p in (faa, fna, gff) if not p.is_file()]
    if missing:
        raise RuntimeError(f"{stem}: missing Bakta output(s) {', '.join(missing)}")
    out_tsv = amrfinder_tsv_path(output_path, stem)
    tmp_tsv = out_tsv.with_name(f".{out_tsv.name}.tmp")
    cmd = [
        *tool_cmd("amrfinder"),
        "--protein", str(faa),
        "--nucleotide", str(fna),
        "--gff", str(gff),
        "--annotation_format", "bakta",
        *(["--organism", organism] if organism else []),
        "--plus",
        "--database", str(db_dir),
        "--threads", str(max(1, int(threads))),
        "--name", stem,
        "--output", str(tmp_tsv),
        "--quiet",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if res.returncode != 0 or not tmp_tsv.is_file():
        detail = (res.stderr or res.stdout or f"exit code {res.returncode}").strip().splitlines()
        raise RuntimeError(f"{stem}: amrfinder failed: {' | '.join(detail[-3:])}")
    os.replace(tmp_tsv, out_tsv)
    return out_tsv


def _amrfinder_version() -> str:
    try:
        res = subprocess.run([*tool_cmd("amrfinder"), "--version"], capture_output=True, text=True, check=False)
    except FileNotFoundError:
        return "not found"
    text = (res.stdout or res.stderr or "").strip().splitlines()
    return text[-1].strip() if text else "unknown"


def read_amrfinder_summary(output_path: Path) -> dict[str, Any] | None:
    try:
        data = json.loads((Path(output_path) / AMRFINDER_DIR / AMRFINDER_SUMMARY).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def amrfinder_reuse_issue(output_path: Path, genomes: Sequence[Any], plan: "AmrFinderPlan") -> str:
    """Why the recorded calls do not fit this organism, DB or Bakta annotation ('' when they do)."""
    summary = read_amrfinder_summary(output_path)
    if summary is None:
        return "no AMRFinderPlus summary"
    if (summary.get("organism") or None) != plan.organism:
        return f"organism changed ({summary.get('organism')} -> {plan.organism})"
    if plan.db_dir is None:
        return ""
    if summary.get("skipped"):
        return "AMRFinderPlus did not run last time"
    if str(summary.get("db_version") or "") != plan.db_version:
        return f"AMRFinderPlus DB changed ({summary.get('db_version')} -> {plan.db_version})"
    if summary.get("failed"):
        return f"AMRFinderPlus failed for {len(summary['failed'])} genome(s) last time"
    for genome in genomes:
        stem = Path(str(genome)).stem
        report, proteins = amrfinder_tsv_path(output_path, stem), _bakta_inputs(output_path, stem)[0]
        try:
            if report.stat().st_mtime_ns < proteins.stat().st_mtime_ns:
                return f"Bakta reannotated {stem}"
        except OSError:
            return f"missing report or Bakta proteins for {stem}"
    return ""


@dataclass(frozen=True)
class AmrFinderPlan:
    organism: str | None
    note: str
    db_dir: Path | None
    db_version: str


def plan_amrfinder(output_path: Path, *, bakta_db: str | None, schema_dir: str | None = None) -> AmrFinderPlan:
    """The organism and DB the step would use now (no DB: the step is skipped; no organism:
    AMRFinderPlus runs without point mutation calls)."""
    db_dir = amrfinder_db_dir(bakta_db)
    if db_dir is None:
        return AmrFinderPlan(None, f"AMRFinderPlus DB not found under the Bakta DB ({bakta_db})", None, "")
    organism, note = resolve_amrfinder_organism(output_path, db_dir, schema_dir=schema_dir)
    return AmrFinderPlan(organism, note, db_dir, amrfinder_db_version(db_dir))


def run_amrfinder_mutations(
    genome_files: Sequence[Path],
    output_path: Path,
    logger: logging.Logger,
    *,
    plan: AmrFinderPlan,
    threads: int = 1,
    progress: ProgressCallback = None,
) -> AmrFinderResult:
    """Run AMRFinderPlus on every genome's Bakta outputs and write amrfinder/summary.json.

    Without an AMRFinderPlus DB the step writes a summary that says so and makes no calls.
    When the species has no AMRFinderPlus organism it runs without --organism (acquired
    genes only: the run simply has no point mutation calls). A genome
    that fails is listed under 'failed' in the summary (no report, so no calls for it) and
    does not stop the run; the next run retries it.
    """
    output_path = Path(output_path)
    out_dir = output_path / AMRFINDER_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    stems = [Path(g).stem for g in genome_files]
    organism, db_dir, db_version = plan.organism, plan.db_dir, plan.db_version
    for stale in out_dir.glob("*.tsv"):
        stale.unlink()

    summary: dict[str, Any] = {
        "organism": organism,
        "organism_source": plan.note,
        "amrfinder_version": _amrfinder_version() if db_dir else None,
        "db_dir": str(db_dir) if db_dir else None,
        "db_version": db_version or None,
        "samples": sorted(stems) if db_dir else [],
        "rows": 0,
        "point_rows": 0,
    }
    if db_dir is None:
        logger.warning("Skipping AMRFinderPlus: %s.", plan.note)
        summary["skipped"] = plan.note
        _write_summary(out_dir, summary)
        return AmrFinderResult(out_dir, out_dir / AMRFINDER_SUMMARY, None, 0)
    if organism is None:
        logger.warning("AMRFinderPlus runs without --organism, so without point mutation calls: %s.", plan.note)

    per_genome = min(AMRFINDER_THREADS_PER_GENOME, max(1, int(threads)))
    workers = max(1, min(len(stems), int(threads) // per_genome))
    logger.info(
        "Running AMRFinderPlus (--organism %s, DB %s; %s) on %d genome(s), %d at a time.",
        organism or "none", db_version, plan.note, len(stems), workers,
    )
    failed: dict[str, str] = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(_run_one, stem, output_path, db_dir, organism, per_genome): stem
            for stem in stems
        }
        report_progress(progress, 0, len(stems))
        for done, future in enumerate(as_completed(futures), start=1):
            try:
                future.result()
            except Exception as exc:
                failed[futures[future]] = str(exc)
            report_progress(progress, done, len(stems))
    if failed:
        summary["failed"] = dict(sorted(failed.items()))
        logger.warning(
            "AMRFinderPlus failed for %d of %d genome(s); they have no point mutation calls "
            "(the next run retries them): %s",
            len(failed), len(stems), "; ".join(list(summary["failed"].values())[:5]),
        )

    calls = load_amrfinder_calls(output_path, stems)
    summary["rows"] = int(len(calls))
    summary["point_rows"] = int(len(_amr_point_rows(calls)))
    _write_summary(out_dir, summary)
    logger.info("AMRFinderPlus: %d rows, %d AMR point mutation rows.", summary["rows"], summary["point_rows"])
    return AmrFinderResult(
        out_dir, out_dir / AMRFINDER_SUMMARY, organism, summary["point_rows"], tuple(sorted(failed))
    )


def _write_summary(out_dir: Path, summary: Mapping[str, Any]) -> None:
    path = out_dir / AMRFINDER_SUMMARY
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    os.replace(tmp, path)


_CALL_COLUMNS = {
    "Name": "Sample",
    "Protein id": "ProteinId",
    "Contig id": "ContigId",
    "Start": "Start",
    "Stop": "Stop",
    "Strand": "Strand",
    "Element symbol": "Symbol",
    "Element name": "ElementName",
    "Type": "Type",
    "Subtype": "Subtype",
    "Class": "Class",
    "Subclass": "Subclass",
}


def load_amrfinder_calls(output_path: Path, samples: Sequence[str] | None = None) -> pd.DataFrame:
    """Every row of the per-genome AMRFinderPlus reports (all types), with short column names."""
    out_dir = Path(output_path) / AMRFINDER_DIR
    if samples is None:
        paths = sorted(out_dir.glob("*.tsv"))
    else:
        paths = [amrfinder_tsv_path(output_path, s) for s in samples]
    frames: list[pd.DataFrame] = []
    for path in paths:
        if not path.is_file():
            continue
        with path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if not rows:
            continue
        frame = pd.DataFrame(rows)
        frame = frame.rename(columns=_CALL_COLUMNS)
        if "Sample" not in frame.columns:
            frame["Sample"] = path.stem
        frames.append(frame[[c for c in _CALL_COLUMNS.values() if c in frame.columns]])
    if not frames:
        return pd.DataFrame(columns=list(_CALL_COLUMNS.values()))
    return pd.concat(frames, ignore_index=True)


def _join_unique(values) -> str:
    seen: list[str] = []
    for value in values:
        text = str(value or "").strip()
        if text and text.upper() != "NA" and text not in seen:
            seen.append(text)
    return ";".join(seen)


def _mutation_classes(class_text: str, subclass_text: str) -> list[str]:
    """'QUINOLONE' / 'QUINOLONE' -> ['quinolone']; 'BETA-LACTAM' / 'CARBAPENEM' -> both."""
    out: list[str] = []
    for text in (subclass_text, class_text):
        for part in str(text or "").split("/"):
            label = part.strip().lower()
            if label and label != "na" and label not in out:
                out.append(label)
    return out


def _amr_point_rows(calls: pd.DataFrame) -> pd.DataFrame:
    """The AMR point mutation rows (Type AMR, Subtype POINT or POINT_DISRUPT).

    STRESS or VIRULENCE point mutations (reported with --plus) are left out, so the
    X tag always means resistance.
    """
    if calls.empty or "Subtype" not in calls.columns or "Type" not in calls.columns:
        return calls.iloc[0:0]
    keep = calls["Subtype"].astype(str).str.upper().isin(POINT_SUBTYPES)
    keep &= calls["Type"].astype(str).str.upper().eq("AMR")
    return calls[keep]


def point_mutation_maps(output_path: Path) -> tuple[dict[tuple[str, str], tuple[str, str]], dict[tuple[str, str, int, int], tuple[str, str]]]:
    """POINT calls keyed by (sample, locus tag) and by (sample, contig, start, stop).

    The values are (';'-joined mutation symbols, ';'-joined drug classes).
    """
    calls = _amr_point_rows(load_amrfinder_calls(output_path)).copy()
    if calls.empty:
        return {}, {}
    disrupted = calls["Subtype"].astype(str).str.upper().eq("POINT_DISRUPT")
    calls.loc[disrupted, "Symbol"] = calls.loc[disrupted, "Symbol"].astype(str) + DISRUPTION_SUFFIX
    calls["Classes"] = [
        ";".join(_mutation_classes(c, s)) for c, s in zip(calls.get("Class", ""), calls.get("Subclass", ""))
    ]
    by_locus: dict[tuple[str, str], tuple[str, str]] = {}
    by_coords: dict[tuple[str, str, int, int], tuple[str, str]] = {}
    for (sample, protein, contig, start, stop), group in calls.groupby(
        ["Sample", "ProteinId", "ContigId", "Start", "Stop"], sort=False, dropna=False
    ):
        value = (_join_unique(group["Symbol"]), _join_unique(";".join(group["Classes"]).split(";")))
        protein = str(protein or "").strip()
        if protein and protein.upper() != "NA":
            by_locus[(str(sample), protein)] = value
        try:
            by_coords[(str(sample), _contig_key(contig), int(start), int(stop))] = value
        except (TypeError, ValueError):
            continue
    return by_locus, by_coords


def _contig_key(value: Any) -> str:
    """'contig_12', 12 and '12' all key as '12' (Bakta rows keep only the contig number)."""
    text = str(value if value is not None else "").strip()
    if text.lower().startswith("contig_"):
        text = text.split("_", 1)[1]
    try:
        number = float(text)
    except ValueError:
        return text
    return str(int(number)) if number.is_integer() else text


def attach_resistance_mutations(df: pd.DataFrame, output_path: Path) -> pd.DataFrame:
    """Add ResistanceMutations / ResistanceMutationClasses ('' when none) to Bakta CDS rows.

    A CDS is matched by its Bakta locus tag (AMRFinderPlus 'Protein id'), falling back to
    the contig and coordinates.
    """
    by_locus, by_coords = point_mutation_maps(output_path)
    mutations = [""] * len(df)
    classes = [""] * len(df)
    used: set[tuple[str, str]] = set()
    if (by_locus or by_coords) and not df.empty:
        samples = df["Sample"].astype(str).tolist()
        loci = df["Locus"].astype(str).tolist() if "Locus" in df.columns else [""] * len(df)
        contigs = df["Contig"].tolist() if "Contig" in df.columns else [None] * len(df)
        starts = pd.to_numeric(df.get("Start"), errors="coerce").tolist()
        stops = pd.to_numeric(df.get("Stop"), errors="coerce").tolist()
        for i, (sample, locus, contig, start, stop) in enumerate(zip(samples, loci, contigs, starts, stops)):
            hit = by_locus.get((sample, locus))
            if hit is None and pd.notna(start) and pd.notna(stop):
                hit = by_coords.get((sample, _contig_key(contig), int(start), int(stop)))
            if hit is not None:
                mutations[i], classes[i] = hit
                used.add((sample, hit[0]))
    # A call outside every CDS (e.g. a promoter mutation from the AMR_DNA tables) has no
    # row to tag: say so in the log rather than lose it silently.
    calls = {(key[0], value[0]) for found in (by_locus, by_coords) for key, value in found.items()}
    unmatched = sorted(calls - used)
    if unmatched:
        _LOG.warning(
            "%d AMRFinderPlus point mutation call(s) match no Bakta CDS and are not shown: %s",
            len(unmatched), ", ".join(f"{s} {sym}" for s, sym in unmatched[:10]),
        )
    out = df.copy()
    out["ResistanceMutations"] = mutations
    out["ResistanceMutationClasses"] = classes
    return out


__all__ = [
    "AMRFINDER_DIR",
    "AMRFINDER_SUMMARY",
    "AmrFinderPlan",
    "AmrFinderResult",
    "DISRUPTION_SUFFIX",
    "MUTATION_COLUMNS",
    "POINT_SUBTYPES",
    "amrfinder_db_dir",
    "amrfinder_db_version",
    "amrfinder_reuse_issue",
    "amrfinder_tsv_path",
    "attach_resistance_mutations",
    "load_amrfinder_calls",
    "plan_amrfinder",
    "point_mutation_maps",
    "read_amrfinder_summary",
    "resolve_amrfinder_organism",
    "run_amrfinder_mutations",
]
