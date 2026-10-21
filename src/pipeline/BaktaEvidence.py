"""Per-gene evidence from the run's own outputs (step "bakta_evidence").

BACTolomew ships no marker databases of its own; the gene categories come from what the
pipeline already computes:

- AMR, point mutations (X), stress/metal/biocide (BMRG) and AMRFinderPlus "plus"
  virulence: the amrfinder/ step's reports (AMRFinderPlus --plus, its DB from the Bakta DB),
  not Bakta's own AMRFinderPlus hits (a subset of them). A call on no Bakta CDS (a gene cut
  by a contig end, a partial or frameshifted copy Bakta did not call) is kept on its own
  contig coordinates with no locus tag; the explorer shows it as a gene of its own.
- VIRULENCE: Bakta's VFDB expert hits and VFDB: db_xrefs.
- MOBILE: Bakta's IS expert hits and IS: db_xrefs, and Bakta products naming an IS element
  or a transposase (marked as a product name match).
- PHAGE: Bakta has no phage database; only Bakta products with clear phage terms, marked
  as a product name match.
- PLASMID: MOB-suite's replicon and relaxase (and MPF) biomarker hits, on the CDS they cover.

One row per (gene, source, reference) is written to evidence/<sample>.tsv, with
evidence/summary.json recording the inputs (sha256), the sources each tag had and the counts.
The explorer reads these rows as it reads the marker search hits (the same tags and views).
"""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import os
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import pandas as pd

from pipeline.AmrFinderWorkflow import (
    POINT_SUBTYPES,
    amrfinder_tsv_path,
    read_amrfinder_summary,
)

EVIDENCE_DIR = "evidence"
EVIDENCE_SUMMARY = "summary.json"
# Bump when the rules below change: a run's evidence is then rebuilt on resume.
EVIDENCE_VERSION = 3

EVIDENCE_COLUMNS = (
    "sample", "contig", "start", "stop", "strand", "locus_tag",
    "category", "tag", "source", "method", "id", "gene", "product",
    "identity", "coverage", "class", "subclass",
)

SOURCE_BAKTA_VFDB = "Bakta VFDB"
SOURCE_BAKTA_IS = "Bakta IS"
SOURCE_PRODUCT_NAME = "Bakta product name"
SOURCE_AMRFINDER = "AMRFinderPlus"
SOURCE_MOBSUITE = "MOB-suite"
# Every source name the step writes (a user database cannot take one, pipeline.UserDatabases).
EVIDENCE_SOURCES = (SOURCE_BAKTA_VFDB, SOURCE_BAKTA_IS, SOURCE_PRODUCT_NAME, SOURCE_AMRFINDER, SOURCE_MOBSUITE)
# A product name match is a keyword rule on Bakta's product, not a database search.
KEYWORD_METHOD = "product name match"
# AMRFinderPlus's method for a gene with a premature stop codon (a disrupted copy).
DISRUPTED_METHOD = "INTERNAL_STOP"

CATEGORY_TAGS: dict[str, str] = {
    "AMR": "A",
    "BMRG": "B",
    "VIRULENCE": "V",
    "MOBILE": "M",
    "PHAGE": "P",
    "PLASMID": "",
}
POINT_TAG = "X"

# VFDB's top level virulence factor classes (the VFC ids of Bakta's VFDB: xrefs).
VFDB_CLASSES: dict[str, str] = {
    "VFC0001": "Adherence",
    "VFC0083": "Invasion",
    "VFC0086": "Effector delivery system",
    "VFC0204": "Motility",
    "VFC0235": "Exotoxin",
    "VFC0251": "Exoenzyme",
    "VFC0258": "Immune modulation",
    "VFC0271": "Biofilm",
    "VFC0272": "Nutritional/Metabolic factor",
    "VFC0282": "Stress survival",
    "VFC0301": "Regulation",
    "VFC0315": "Post-translational modification",
    "VFC0325": "Antimicrobial activity/Competitive advantage",
    "VFC0346": "Others",
}

# IS elements and transposases by name. IS names are matched case sensitively (IS26,
# ISAba1, ISL3, IS3-like), so words such as "is" never count.
_MOBILE_PRODUCT_RE = re.compile(
    r"(?i:transposase|insertion sequence|insertion element)"
    r"|\bIS[A-Z]?[A-Za-z]{0,3}\d+|\bIS\d*-like"
)
# Clear phage terms. "holin" as a word (not choline), "portal protein" (not a transport portal).
_PHAGE_PRODUCT_RE = re.compile(
    r"(?i)phage|terminase|portal protein|capsid|tail fib(?:er|re)|tape[ -]measure|\bholins?\b|baseplate"
)
# Never phage evidence by name: generic integrases, transposases and recombinases (often
# "phage integrase family"), host genes named after phage (phage shock, anti-phage defense,
# abortive infection, FxsA phage exclusion), type VI secretion baseplate subunits and the
# Phage-Barnase toxin nuclease domain.
_PHAGE_EXCLUDE_RE = re.compile(
    r"(?i:integrase|transposase|recombinase|resolvase|phage shock|anti-?phage|phage (?:resistance|defen[cs]e)"
    r"|abortive infection|exclusion|FxsA|type VI|T6SS|Phage-Barnase)|\bTss[A-Z]"
)
_IS_FAMILY_RE = re.compile(r"\b(IS[A-Za-z]*\d+|Tn\d+)(?:/\w+)? family\b")

_LOG = logging.getLogger(__name__)


# Matching ------------------------------------------------------------------

def mobile_product_match(product: str) -> str:
    """The IS / transposase term of a product ('' when it names none)."""
    match = _MOBILE_PRODUCT_RE.search(str(product or ""))
    return match.group(0) if match else ""


def phage_product_match(product: str) -> str:
    """The phage term of a product ('' when it names none or is excluded)."""
    text = str(product or "")
    if _PHAGE_EXCLUDE_RE.search(text):
        return ""
    match = _PHAGE_PRODUCT_RE.search(text)
    return match.group(0) if match else ""


def contig_key(value: Any) -> str:
    """'contig_12', 12 and '12' all key as '12'."""
    text = str(value if value is not None else "").strip()
    if text.lower().startswith("contig_"):
        text = text.split("_", 1)[1]
    try:
        number = float(text)
    except ValueError:
        return text
    return str(int(number)) if number.is_integer() else text


def _pct(value: Any, *, fraction: bool = False) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return ""
    if number != number:  # NaN
        return ""
    return f"{number * 100 if fraction else number:.1f}"


def _text(value: Any) -> str:
    text = str(value if value is not None else "").strip()
    return "" if text.upper() in {"NA", "NAN", "NONE", "-"} else text


# Sources -------------------------------------------------------------------

def _cds_features(bakta: Mapping[str, Any]) -> list[dict[str, Any]]:
    return [f for f in bakta.get("features") or [] if isinstance(f, dict) and f.get("type") == "cds"]


def _row(sample: str, feature: Mapping[str, Any], **fields: Any) -> dict[str, str]:
    category = fields.pop("category")
    row = {
        "sample": sample,
        "contig": str(feature.get("sequence") or ""),
        "start": str(feature.get("start") or ""),
        "stop": str(feature.get("stop") or ""),
        "strand": str(feature.get("strand") or ""),
        "locus_tag": str(feature.get("locus") or ""),
        "category": category,
        "tag": fields.pop("tag", CATEGORY_TAGS.get(category, "")),
    }
    for col in EVIDENCE_COLUMNS:
        if col not in row:
            row[col] = _text(fields.get(col, ""))
    return row


def _best_expert(feature: Mapping[str, Any], source: str) -> dict[str, Any] | None:
    hits = [e for e in feature.get("expert") or [] if isinstance(e, dict) and e.get("source") == source]
    if not hits:
        return None
    return max(hits, key=lambda e: (float(e.get("identity") or 0), float(e.get("score") or 0)))


def _xrefs(values: Iterable[Any], prefix: str) -> list[str]:
    return [str(x).split(":", 1)[1] for x in values or [] if str(x).startswith(prefix + ":")]


def bakta_evidence_rows(sample: str, bakta: Mapping[str, Any]) -> list[dict[str, str]]:
    """VFDB and IS (expert hits, else db_xrefs) and the MOBILE / PHAGE product name matches."""
    rows: list[dict[str, str]] = []
    for feature in _cds_features(bakta):
        xrefs = feature.get("db_xrefs") or []
        product = str(feature.get("product") or "")

        vf = _best_expert(feature, "VFDB")
        vf_xrefs = _xrefs([*(vf.get("db_xrefs") or [] if vf else []), *xrefs], "VFDB")
        vfg = next((x for x in vf_xrefs if x.startswith("VFG")), "")
        vfc = next((x for x in vf_xrefs if x.startswith("VFC")), "")
        if vf or vfg:
            rows.append(_row(
                sample, feature, category="VIRULENCE", source=SOURCE_BAKTA_VFDB,
                method="Bakta expert" if vf else "Bakta xref",
                id=(vf or {}).get("id") or vfg,
                gene=(vf or {}).get("gene") or feature.get("gene"),
                product=(vf or {}).get("product") or product,
                identity=_pct((vf or {}).get("identity"), fraction=True),
                coverage=_pct((vf or {}).get("subject_cov"), fraction=True),
                **{"class": VFDB_CLASSES.get(vfc, vfc), "subclass": vfc},
            ))

        is_hit = _best_expert(feature, "IS")
        is_names = _xrefs([*(is_hit.get("db_xrefs") or [] if is_hit else []), *xrefs], "IS")
        if is_hit or is_names:
            is_product = (is_hit or {}).get("product") or product
            family = _IS_FAMILY_RE.search(is_product)
            rows.append(_row(
                sample, feature, category="MOBILE", source=SOURCE_BAKTA_IS,
                method="Bakta expert" if is_hit else "Bakta xref",
                id=is_names[0] if is_names else (is_hit or {}).get("id"),
                gene=(is_hit or {}).get("gene") or feature.get("gene"),
                product=is_product,
                identity=_pct((is_hit or {}).get("identity"), fraction=True),
                coverage=_pct((is_hit or {}).get("subject_cov"), fraction=True),
                **{"class": family.group(1) if family else ""},
            ))

        mobile_term = mobile_product_match(product)
        if mobile_term:
            rows.append(_row(
                sample, feature, category="MOBILE", source=SOURCE_PRODUCT_NAME,
                method=KEYWORD_METHOD, id=mobile_term, gene=feature.get("gene"), product=product,
            ))
        phage_term = phage_product_match(product)
        if phage_term:
            rows.append(_row(
                sample, feature, category="PHAGE", source=SOURCE_PRODUCT_NAME,
                method=KEYWORD_METHOD, id=phage_term, gene=feature.get("gene"), product=product,
            ))
    return rows


def _feature_index(features: Sequence[Mapping[str, Any]]) -> tuple[dict[str, Mapping[str, Any]], dict[str, list[Mapping[str, Any]]]]:
    by_locus = {str(f.get("locus")): f for f in features if f.get("locus")}
    by_contig: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for f in features:
        by_contig[contig_key(f.get("sequence"))].append(f)
    return by_locus, by_contig


def _overlapping_cds(
    by_contig: Mapping[str, list[Mapping[str, Any]]], contig: Any, start: Any, stop: Any
) -> Mapping[str, Any] | None:
    """The CDS that a hit covers best: at least half of the shorter of the two overlap."""
    try:
        lo, hi = sorted((int(float(start)), int(float(stop))))
    except (TypeError, ValueError):
        return None
    best, best_overlap = None, 0
    for f in by_contig.get(contig_key(contig), ()):
        f_lo, f_hi = sorted((int(f.get("start") or 0), int(f.get("stop") or 0)))
        overlap = min(hi, f_hi) - max(lo, f_lo) + 1
        if overlap <= 0 or overlap * 2 < min(hi - lo + 1, f_hi - f_lo + 1):
            continue
        if overlap > best_overlap:
            best, best_overlap = f, overlap
    return best


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _call_site(call: Mapping[str, Any]) -> dict[str, str]:
    """A call's own contig coordinates, as a feature without a locus tag."""
    return {
        "sequence": _text(call.get("Contig id")),
        "start": _text(call.get("Start")),
        "stop": _text(call.get("Stop")),
        "strand": str(call.get("Strand") or "").strip(),  # "-" is a strand here, not a missing value
        "locus": "",
    }


def amrfinder_evidence_rows(sample: str, report: Path, features: Sequence[Mapping[str, Any]]) -> list[dict[str, str]]:
    """AMR (A), AMR point mutations (X), STRESS (BMRG) and VIRULENCE rows of one
    AMRFinderPlus report, on the CDS of their Bakta locus tag (else the CDS they cover, at
    least half of the shorter of the two). A gene call on no CDS keeps its own contig
    coordinates and an empty locus tag; a point mutation outside every CDS (a promoter) is
    left out, it has no gene to tag."""
    by_locus, by_contig = _feature_index(features)
    rows: list[dict[str, str]] = []
    for call in _read_tsv_rows(report):
        kind = _text(call.get("Type")).upper()
        subtype = _text(call.get("Subtype")).upper()
        point = subtype in POINT_SUBTYPES
        if kind == "AMR":
            category, tag = "AMR", POINT_TAG if point else "A"
        elif point:
            continue  # STRESS / VIRULENCE point mutations: X always means resistance
        elif kind == "STRESS":
            category, tag = "BMRG", "B"
        elif kind == "VIRULENCE":
            category, tag = "VIRULENCE", "V"
        else:
            continue
        feature = by_locus.get(_text(call.get("Protein id")))
        if feature is None:
            feature = _overlapping_cds(by_contig, call.get("Contig id"), call.get("Start"), call.get("Stop"))
        if feature is None:
            if point:
                continue  # outside every CDS (a promoter mutation): no gene to tag
            feature = _call_site(call)
        rows.append(_row(
            sample, feature, category=category, tag=tag, source=SOURCE_AMRFINDER,
            method=call.get("Method"),
            id=call.get("Closest reference accession") or call.get("HMM accession"),
            gene=call.get("Element symbol"),
            product=call.get("Element name"),
            identity=_pct(call.get("% Identity to reference")),
            coverage=_pct(call.get("% Coverage of reference")),
            **{"class": call.get("Class"), "subclass": call.get("Subclass")},
        ))
    return rows


def mobsuite_biomarkers_path(output_dir: Path, sample: str) -> Path:
    return Path(output_dir) / "mobsuite" / "per_sample" / sample / "biomarkers.blast.txt"


def mobsuite_contig_report_path(output_dir: Path, sample: str) -> Path:
    return Path(output_dir) / "mobsuite" / "per_sample" / sample / "contig_report.txt"


def mobsuite_plasmid_contigs(contig_report: Path) -> set[str] | None:
    """Contig keys MOB-suite calls plasmid (None when the report is missing)."""
    if not contig_report.is_file():
        return None
    return {
        contig_key(row.get("contig_id"))
        for row in _read_tsv_rows(contig_report)
        if _text(row.get("molecule_type")).lower() == "plasmid"
    }


def mobsuite_evidence_rows(
    sample: str,
    biomarkers: Path,
    features: Sequence[Mapping[str, Any]],
    plasmid_contigs: set[str] | None = None,
) -> list[dict[str, str]]:
    """MOB-suite replicon, relaxase and MPF hits (PLASMID) on the CDS they cover, only on
    contigs MOB-suite calls plasmid (every contig when `plasmid_contigs` is None)."""
    _by_locus, by_contig = _feature_index(features)
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for hit in _read_tsv_rows(biomarkers):
        kind = _text(hit.get("biomarker")).lower()
        if kind not in {"replicon", "relaxase", "mpf"}:
            continue
        if plasmid_contigs is not None and contig_key(hit.get("sseqid")) not in plasmid_contigs:
            continue
        feature = _overlapping_cds(by_contig, hit.get("sseqid"), hit.get("sstart"), hit.get("send"))
        if feature is None:
            continue
        reference = _text(hit.get("qseqid"))
        name = reference.rsplit("|", 1)[-1] if reference else ""
        key = (str(feature.get("locus") or ""), name)
        if key in seen:
            continue
        seen.add(key)
        rows.append(_row(
            sample, feature, category="PLASMID", source=SOURCE_MOBSUITE, method="mob_recon",
            id=reference, gene=name, product=f"{kind} {name}".strip(),
            identity=_pct(hit.get("pident")), coverage=_pct(hit.get("qcovhsp")),
            **{"class": kind},
        ))
    return rows


# The step --------------------------------------------------------------------

def evidence_dir(output_dir: Path) -> Path:
    return Path(output_dir) / EVIDENCE_DIR


def evidence_tsv_path(output_dir: Path, sample: str) -> Path:
    return evidence_dir(output_dir) / f"{sample}.tsv"


def bakta_json_path(output_dir: Path, sample: str) -> Path:
    return Path(output_dir) / "bakta" / sample / f"{sample}.json"


def _sha256(path: Path) -> str | None:
    try:
        digest = hashlib.sha256()
        with path.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 22), b""):
                digest.update(chunk)
        return digest.hexdigest()
    except OSError:
        return None


def _amrfinder_ready(output_dir: Path) -> bool:
    """AMRFinderPlus ran (with or without an organism) and wrote reports."""
    summary = read_amrfinder_summary(output_dir) or {}
    return bool(summary.get("samples")) and not summary.get("skipped")


def _point_mutations_ready(output_dir: Path) -> bool:
    summary = read_amrfinder_summary(output_dir) or {}
    return _amrfinder_ready(output_dir) and bool(summary.get("organism"))


def _sample_inputs(output_dir: Path, sample: str) -> dict[str, str | None]:
    report = amrfinder_tsv_path(output_dir, sample)
    return {
        "bakta_json": _sha256(bakta_json_path(output_dir, sample)),
        "amrfinder": _sha256(report) if _amrfinder_ready(output_dir) and report.is_file() else None,
        "mobsuite_biomarkers": _sha256(mobsuite_biomarkers_path(output_dir, sample)),
        "mobsuite_contigs": _sha256(mobsuite_contig_report_path(output_dir, sample)),
    }


def read_evidence_summary(output_dir: Path) -> dict[str, Any] | None:
    try:
        data = json.loads((evidence_dir(output_dir) / EVIDENCE_SUMMARY).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def evidence_reuse_issue(output_dir: Path, genomes: Sequence[Any]) -> str:
    """Why the recorded evidence does not fit the current inputs ('' when it does)."""
    summary = read_evidence_summary(output_dir)
    if summary is None:
        return "no evidence summary"
    if summary.get("version") != EVIDENCE_VERSION:
        return f"evidence rules changed (version {summary.get('version')} -> {EVIDENCE_VERSION})"
    recorded = summary.get("inputs") if isinstance(summary.get("inputs"), dict) else {}
    samples = [Path(str(g)).stem for g in genomes]
    if sorted(recorded) != sorted(samples):
        return "the genome set changed"
    for sample in samples:
        if not evidence_tsv_path(output_dir, sample).is_file():
            return f"missing evidence table for {sample}"
        if recorded.get(sample) != _sample_inputs(output_dir, sample):
            return f"inputs changed for {sample}"
    return ""


def _sources_by_tag(output_dir: Path, *, bakta_ok: bool, mobsuite_ok: bool) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    """The evidence sources each tag had in this run, and the ones it lacked."""
    amr = _amrfinder_ready(output_dir)
    point = _point_mutations_ready(output_dir)
    planned: dict[str, list[tuple[str, bool]]] = {
        "A": [(SOURCE_AMRFINDER, amr)],
        "B": [(SOURCE_AMRFINDER, amr)],
        "V": [(SOURCE_BAKTA_VFDB, bakta_ok), (SOURCE_AMRFINDER, amr)],
        "M": [(SOURCE_BAKTA_IS, bakta_ok), (SOURCE_PRODUCT_NAME, bakta_ok)],
        "P": [(SOURCE_PRODUCT_NAME, bakta_ok)],
        "X": [(SOURCE_AMRFINDER, point)],
        "PLASMID": [(SOURCE_MOBSUITE, mobsuite_ok)],
    }
    tags = {tag: [name for name, ok in items if ok] for tag, items in planned.items()}
    missing = {tag: [name for name, ok in items if not ok] for tag, items in planned.items()}
    return tags, {tag: names for tag, names in missing.items() if names}


def build_evidence(
    output_dir: Path,
    genomes: Sequence[Any],
    logger: logging.Logger | None = None,
) -> dict[str, Any]:
    """Write evidence/<sample>.tsv for every genome and evidence/summary.json; returns the summary."""
    log = logger or _LOG
    output_dir = Path(output_dir)
    out_dir = evidence_dir(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    samples = [Path(str(g)).stem for g in genomes]
    amr_ready = _amrfinder_ready(output_dir)
    amr_summary = read_amrfinder_summary(output_dir) or {}
    amr_failed = set((amr_summary.get("failed") or {}).keys()) if isinstance(amr_summary.get("failed"), dict) else set()

    counts: Counter = Counter()
    unplaced: dict[str, int] = {}
    genes: dict[tuple[str, str], set[tuple[str, str]]] = defaultdict(set)
    bakta_db: dict[str, Any] = {}
    no_biomarkers: list[str] = []
    inputs: dict[str, dict[str, str | None]] = {}
    for sample in samples:
        json_path = bakta_json_path(output_dir, sample)
        bakta = json.loads(json_path.read_text(encoding="utf-8"))
        if not bakta_db and isinstance(bakta.get("version"), dict):
            bakta_db = {"bakta": bakta["version"].get("bakta"), **(bakta["version"].get("db") or {})}
        features = _cds_features(bakta)
        rows = bakta_evidence_rows(sample, bakta)
        report = amrfinder_tsv_path(output_dir, sample)
        if amr_ready and report.is_file():
            amr_rows = amrfinder_evidence_rows(sample, report, features)
            off_cds = [r for r in amr_rows if not r["locus_tag"]]
            if off_cds:
                unplaced[sample] = len(off_cds)
                log.info(
                    "AMRFinderPlus: %d call(s) of %s on no Bakta CDS, kept on their contig coordinates: %s.",
                    len(off_cds), sample,
                    ", ".join(f"{r['gene']} ({r['contig']}:{r['start']}-{r['stop']})" for r in off_cds),
                )
            rows += amr_rows
        biomarkers = mobsuite_biomarkers_path(output_dir, sample)
        if biomarkers.is_file():
            plasmid_contigs = mobsuite_plasmid_contigs(mobsuite_contig_report_path(output_dir, sample))
            rows += mobsuite_evidence_rows(sample, biomarkers, features, plasmid_contigs)
        else:
            no_biomarkers.append(sample)
        rows.sort(key=lambda r: (contig_key(r["contig"]).zfill(8), int(r["start"] or 0), r["category"], r["source"]))
        path = evidence_tsv_path(output_dir, sample)
        tmp = path.with_name(f".{path.name}.tmp")
        with tmp.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(EVIDENCE_COLUMNS), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp, path)
        for r in rows:
            counts[(r["category"], r["source"])] += 1
            genes[(r["category"], r["source"])].add((sample, r["locus_tag"] or unplaced_locus(r)))
        inputs[sample] = _sample_inputs(output_dir, sample)

    for stale in out_dir.glob("*.tsv"):
        if stale.stem not in samples:
            stale.unlink()

    tags, missing = _sources_by_tag(
        output_dir, bakta_ok=bool(samples), mobsuite_ok=len(no_biomarkers) < len(samples)
    )
    summary: dict[str, Any] = {
        "version": EVIDENCE_VERSION,
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "bakta_db": bakta_db,
        "amrfinder": {
            "used": amr_ready,
            "organism": amr_summary.get("organism"),
            "db_version": amr_summary.get("db_version"),
            "failed": sorted(amr_failed),
            # Gene calls on no Bakta CDS, per genome (kept on their contig coordinates).
            "unplaced": dict(sorted(unplaced.items())),
        },
        "tags": tags,
        "missing": missing,
        "rows": {f"{c}|{s}": n for (c, s), n in sorted(counts.items())},
        "genes": {f"{c}|{s}": len(v) for (c, s), v in sorted(genes.items())},
        "inputs": inputs,
    }
    path = out_dir / EVIDENCE_SUMMARY
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    os.replace(tmp, path)
    log.info(
        "Bakta evidence: %d rows for %d genome(s) (%s).",
        sum(counts.values()), len(samples),
        ", ".join(f"{c} {s}: {len(v)} genes" for (c, s), v in sorted(genes.items())) or "none",
    )
    if unplaced:
        log.info(
            "AMRFinderPlus: %d gene call(s) in %d genome(s) on no Bakta CDS.",
            sum(unplaced.values()), len(unplaced),
        )
    return summary


def evidence_record(output_dir: Path) -> dict[str, Any] | None:
    """The run_config / run_report record of the evidence sources (None before the step ran)."""
    summary = read_evidence_summary(output_dir)
    if summary is None or summary.get("version") != EVIDENCE_VERSION:
        return None
    return {
        "tags": summary.get("tags") or {},
        "missing": summary.get("missing") or {},
        "bakta_db": summary.get("bakta_db") or {},
        "amrfinder": summary.get("amrfinder") or {},
    }


def validate_evidence_outputs(output_dir: Path, genomes: Sequence[Any]) -> list[str]:
    summary = read_evidence_summary(output_dir)
    if summary is None:
        return [f"Missing or invalid evidence summary: {evidence_dir(output_dir) / EVIDENCE_SUMMARY}"]
    errors: list[str] = []
    for genome in genomes:
        path = evidence_tsv_path(output_dir, Path(str(genome)).stem)
        try:
            header = path.open("r", encoding="utf-8").readline().rstrip("\n").split("\t")
        except OSError:
            errors.append(f"Missing evidence table: {path}")
            continue
        if tuple(header) != EVIDENCE_COLUMNS:
            errors.append(f"Unexpected evidence table header: {path}")
    return errors


# Reading (the explorer) -------------------------------------------------------

def read_evidence(output_dir: Path) -> pd.DataFrame | None:
    """Every evidence row of the run (strings), or None for a run without the step."""
    if read_evidence_summary(output_dir) is None:
        return None
    frames = [
        pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
        for p in sorted(evidence_dir(output_dir).glob("*.tsv"))
    ]
    frames = [f for f in frames if not f.empty]
    if not frames:
        return pd.DataFrame(columns=list(EVIDENCE_COLUMNS))
    return pd.concat(frames, ignore_index=True)


def unplaced_locus(row: Mapping[str, Any]) -> str:
    """The key of a row without a locus tag (an AMRFinderPlus call on no Bakta CDS): its
    contig coordinates, '@contig_3:18-1217'."""
    return f"@{row.get('contig') or ''}:{row.get('start') or ''}-{row.get('stop') or ''}"


def unplaced_evidence(evidence: pd.DataFrame | None) -> pd.DataFrame:
    """The gene rows that lie on no Bakta CDS, one per call site (the explorer adds each
    as a gene of its own)."""
    if evidence is None or evidence.empty:
        return pd.DataFrame(columns=list(EVIDENCE_COLUMNS))
    rows = evidence[(evidence["locus_tag"] == "") & (evidence["tag"] != POINT_TAG)]
    return rows.drop_duplicates(subset=["sample", "contig", "start", "stop"]).reset_index(drop=True)


def _hit_product(row: Mapping[str, Any]) -> str:
    """The text the gene card shows for a row: what was found, its class and, for an
    AMRFinderPlus call, whether the gene is disrupted or lies on no Bakta CDS."""
    gene, product, ident = row.get("gene") or "", row.get("product") or "", row.get("id") or ""
    source = row.get("source")
    if source == SOURCE_PRODUCT_NAME:
        return f"{product} (matched: {ident})"
    if source == SOURCE_BAKTA_IS:
        return f"{ident}: {product}" if ident else product
    head = f"{gene}: {product}" if gene and gene not in product else product
    extra = [x for x in (ident, row.get("class") or "") if x]
    if source == SOURCE_AMRFINDER and row.get("subclass") and row.get("subclass") != row.get("class"):
        extra.append(row["subclass"])
    text = f"{head} ({', '.join(extra)})" if extra else head
    if source == SOURCE_AMRFINDER and row.get("method") == DISRUPTED_METHOD:
        text += "; disrupted: internal stop codon"
    if not row.get("locus_tag"):
        text += "; on no Bakta CDS"
    return text


def evidence_marker_hits(evidence: pd.DataFrame, annotation_ids: Mapping[tuple[str, str], str]) -> pd.DataFrame:
    """The evidence rows as marker hits (DB = the source) of the Bakta CDS they name, keyed
    by (sample, locus tag) -> AnnotationID (a row without a locus tag by (sample,
    `unplaced_locus(row)`)). Point mutation rows are left out: the X tag and the gene
    card's mutation section read them from the AMRFinderPlus reports."""
    cols = [
        "AnnotationID", "Product", "Identity", "Subject_coverage", "Method", "MarkerCategory",
        "Category", "DB", "SubjectIdentifier", "ResistanceMetadataSource", "ResistanceAlleles",
        "ResistanceDrugClasses", "ResistanceAntibiotics",
    ]
    if evidence is None or evidence.empty:
        return pd.DataFrame(columns=cols)
    ev = evidence[evidence["tag"] != POINT_TAG]
    records: list[dict[str, Any]] = []
    for row in ev.to_dict("records"):
        ann_id = annotation_ids.get((row["sample"], row["locus_tag"] or unplaced_locus(row)))
        if not ann_id:
            continue
        amr = row["source"] == SOURCE_AMRFINDER and row["category"] in {"AMR", "BMRG"}
        records.append({
            "AnnotationID": ann_id,
            "Product": _hit_product(row),
            "Identity": row.get("identity"),
            "Subject_coverage": row.get("coverage"),
            "Method": row.get("method") or "",
            "MarkerCategory": row["category"],
            "Category": row["category"],
            "DB": row["source"],
            "SubjectIdentifier": row.get("id") or pd.NA,
            "ResistanceMetadataSource": SOURCE_AMRFINDER if amr else pd.NA,
            "ResistanceAlleles": (row.get("gene") or pd.NA) if amr else pd.NA,
            "ResistanceDrugClasses": (row.get("class") or pd.NA) if amr else pd.NA,
            "ResistanceAntibiotics": (row.get("subclass") or pd.NA) if amr else pd.NA,
        })
    out = pd.DataFrame(records, columns=cols)
    for col in ("Identity", "Subject_coverage"):
        out[col] = pd.to_numeric(out[col], errors="coerce").astype("Float64")
    return out


__all__ = [
    "CATEGORY_TAGS",
    "DISRUPTED_METHOD",
    "EVIDENCE_COLUMNS",
    "EVIDENCE_DIR",
    "EVIDENCE_SOURCES",
    "EVIDENCE_SUMMARY",
    "EVIDENCE_VERSION",
    "KEYWORD_METHOD",
    "SOURCE_AMRFINDER",
    "SOURCE_BAKTA_IS",
    "SOURCE_BAKTA_VFDB",
    "SOURCE_MOBSUITE",
    "SOURCE_PRODUCT_NAME",
    "VFDB_CLASSES",
    "amrfinder_evidence_rows",
    "bakta_evidence_rows",
    "build_evidence",
    "evidence_marker_hits",
    "evidence_record",
    "evidence_reuse_issue",
    "evidence_tsv_path",
    "mobile_product_match",
    "mobsuite_evidence_rows",
    "phage_product_match",
    "read_evidence",
    "read_evidence_summary",
    "unplaced_evidence",
    "unplaced_locus",
    "validate_evidence_outputs",
]
