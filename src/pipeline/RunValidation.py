from __future__ import annotations

import hashlib
import json
import logging
import csv
import re
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

from pipeline.AmrFinderWorkflow import AMRFINDER_DIR, amrfinder_tsv_path, read_amrfinder_summary
from pipeline.BaktaEvidence import validate_evidence_outputs
from pipeline.ChewBBACAWorkflow import CHEWBBACA_OUT_FILES, ensure_chewbbaca_out
from pipeline.AssemblyQualityWorkflow import ASSEMBLY_QC_OUT_FILES
from pipeline.CoreAlignment import (
    CORE_ALIGNMENT_DIR,
    CORE_ALIGNMENT_OUT_FILES,
    ORDERED_ALIGNMENT_REL,
    file_sha256,
    read_ordering_summary,
)
from pipeline.CoreTreeWorkflow import (
    IQTREE_DIR,
    REFIT_DIR,
    REFIT_TREE_ID,
    current_refit_stamp,
    read_core_tree_stamp,
    tree_stamp_path,
)
from pipeline.FastGearWorkflow import validate_sample_clades_tsv
from pipeline.PopulationStructure import method_dir, method_label, run_population_method
from pipeline.GubbinsWorkflow import (
    MASKED_ALIGNMENT_REL,
    SUMMARY_SCHEMA_VERSION as GUBBINS_SCHEMA_VERSION,
    expected_gubbins_files,
    read_gubbins_summary,
)
from pipeline.TreesIndex import LEGACY_SNIPPY_TREE_REL, SNIPPY_RAW_REL, TREES_JSON
from pipeline.MlstWorkflow import MLST_OUT_FILES
from pipeline.MobSuiteWorkflow import MOBSUITE_OUT_FILES


STEP_BAKTA = "bakta"
STEP_ASSEMBLY_QC = "assembly_qc"
STEP_MLST = "mlst"
STEP_MOBSUITE = "mobsuite"
STEP_PANAROO = "panaroo"
STEP_CORE_ALIGNMENT = "core_alignment"
STEP_GUBBINS = "gubbins"
STEP_CORE_TREE = "core_tree"
STEP_PHYLOGENY = "phylogeny"
STEP_FASTGEAR = "fastgear_clades"
STEP_CHEWBBACA = "chewbbaca"
STEP_DIAMOND_MARKERS = "diamond_markers"
STEP_DIAMOND_REFERENCE = "diamond_reference_loci"
STEP_COLLECT_ANNOTATIONS = "collect_annotations"
STEP_AMRFINDER_MUTATIONS = "amrfinder_mutations"
STEP_BAKTA_EVIDENCE = "bakta_evidence"

RESULT_REQUIRED_STEPS: tuple[str, ...] = (
    STEP_ASSEMBLY_QC,
    STEP_BAKTA,
    STEP_MLST,
    STEP_MOBSUITE,
    STEP_PANAROO,
    STEP_GUBBINS,
    STEP_PHYLOGENY,
    STEP_FASTGEAR,
    STEP_CHEWBBACA,
    STEP_DIAMOND_MARKERS,
    STEP_DIAMOND_REFERENCE,
    STEP_COLLECT_ANNOTATIONS,
)
# Steps added after older runs were made: checked when a run recorded them, not required.
RESULT_OPTIONAL_STEPS: tuple[str, ...] = (
    STEP_AMRFINDER_MUTATIONS,
    STEP_BAKTA_EVIDENCE,
    STEP_CORE_ALIGNMENT,
    STEP_CORE_TREE,
)


def _normalize_path_values(values: Iterable[Any]) -> list[Path]:
    out: list[Path] = []
    for value in values:
        text = str(value or "").strip()
        if text:
            out.append(Path(text).expanduser())
    return out


def _stems(values: Iterable[Any]) -> list[str]:
    return [path.stem for path in _normalize_path_values(values)]


GENOME_SET_RECORD_FILES: tuple[str, ...] = ("run_config.json", "run_state.json", "run_report.json")


def recorded_genome_samples(output_dir: Path) -> tuple[set[str], str] | None:
    """Sample names of the genome set an existing run folder was built from.

    Reads inputs.genomes from run_config.json, then run_state.json, then
    run_report.json; a folder with none of them but a Panaroo table falls back to
    the table's sample columns. Returns (samples, source) or None for a fresh folder.
    """
    root = Path(output_dir)
    for name in GENOME_SET_RECORD_FILES:
        payload = read_json_dict(root / name)
        inputs = payload.get("inputs") if isinstance(payload, Mapping) else None
        genomes = inputs.get("genomes") if isinstance(inputs, Mapping) else None
        if isinstance(genomes, list) and genomes:
            return set(_stems(genomes)), name
    samples = panaroo_sample_names(root)
    if samples:
        return set(samples), "panaroo/out/gene_presence_absence.csv"
    return None


def genome_set_change_error(output_dir: Path, genomes: Sequence[Any]) -> str | None:
    """Message refusing to reuse a run folder built from a different genome set."""
    recorded = recorded_genome_samples(output_dir)
    if recorded is None:
        return None
    previous, source = recorded
    current = set(_stems(genomes))
    if previous == current:
        return None
    added = sorted(current - previous)
    removed = sorted(previous - current)

    def _names(values: list[str]) -> str:
        shown = ", ".join(values[:10])
        return shown + (f" and {len(values) - 10} more" if len(values) > 10 else "")

    parts = [
        f"Output folder {Path(output_dir)} already holds a run for a different genome set "
        f"({len(previous)} samples recorded in {source}, {len(current)} given)."
    ]
    if added:
        parts.append(f"Added: {_names(added)}.")
    if removed:
        parts.append(f"Removed: {_names(removed)}.")
    parts.append(
        "Its cohort-wide outputs (Panaroo, Gubbins, trees, fastGEAR, chewBBACA, DIAMOND "
        "marker hits and annotation IDs) belong to the recorded set and would be reused "
        "or attached to the wrong genes. Use a new output folder for this genome set."
    )
    return " ".join(parts)


def panaroo_sample_names(output_dir: Path) -> list[str]:
    path = Path(output_dir) / "panaroo" / "out" / "gene_presence_absence.csv"
    try:
        with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            header = next(csv.reader(handle), [])
    except OSError:
        return []
    return [str(x).strip() for x in header[3:] if str(x).strip()]


def chewbbaca_sample_names(alleles_path: Path) -> list[str]:
    with Path(alleles_path).open("r", encoding="utf-8", errors="ignore") as handle:
        rows = [line.split("\t", 1)[0].strip() for line in handle]
    return [x for x in rows[1:] if x]


def newick_leaf_names(text: str) -> list[str]:
    """Leaf labels of a Newick tree (quoted labels unquoted, comments dropped)."""
    text = re.sub(r"\[[^\]]*\]", "", text)
    return [
        m.group(1).strip().strip("'").replace("''", "'")
        for m in re.finditer(r"(?<=[(,])\s*('(?:[^']|'')*'|[^:,();\s]+)", text)
    ]


def _sample_set_errors(
    found: Iterable[str],
    expected: Iterable[str],
    label: str,
    *,
    allow_missing: bool = False,
) -> list[str]:
    """Errors when a cohort output's samples differ from the run's genomes.

    With allow_missing, samples the tool left out (Gubbins gap filtering, chewBBACA
    inputs without predicted CDSs) are accepted; samples outside the genome set
    always mean the output belongs to another set.
    """
    found_set = set(found)
    expected_set = set(expected)
    missing = [] if allow_missing else sorted(expected_set - found_set)
    extra = sorted(found_set - expected_set)
    if not missing and not extra:
        return []
    detail = []
    if missing:
        detail.append(f"missing {len(missing)} ({', '.join(missing[:5])})")
    if extra:
        detail.append(f"unexpected {len(extra)} ({', '.join(extra[:5])})")
    return [f"{label} samples do not match the run's genomes: {'; '.join(detail)}"]


def tree_sample_names(path: Path, *, reference_stem: str | None = None, genomes: Sequence[Any] | None = None) -> list[str]:
    """Leaf names of a tree file mapped to sample names.

    A Snippy tree names its reference "Reference". When the reference is one of the
    run's genomes that leaf stands for it; when it is an external genome it is not a
    sample and is dropped.
    """
    leaves = newick_leaf_names(Path(path).read_text(encoding="utf-8", errors="ignore"))
    if reference_stem and "Reference" in leaves and reference_stem not in leaves:
        expected = set(_stems(genomes or []))
        if reference_stem in expected:
            leaves = [reference_stem if x == "Reference" else x for x in leaves]
        else:
            leaves = [x for x in leaves if x != "Reference"]
    return leaves


def _tree_sample_errors(
    path: Path,
    genomes: Sequence[Any],
    label: str,
    *,
    reference_stem: str | None = None,
    allow_missing: bool = False,
) -> list[str]:
    try:
        leaves = tree_sample_names(path, reference_stem=reference_stem, genomes=genomes)
    except OSError as exc:
        return [f"Failed reading {label}: {path} ({exc})"]
    return _sample_set_errors(leaves, _stems(genomes), f"{label} {path.name}", allow_missing=allow_missing)


GENOME_SET_STAMP = "genome_set.json"


def genome_set_fingerprint(genomes: Iterable[Any]) -> str:
    """Order-free fingerprint of a genome set (its sorted, unique sample names)."""
    joined = "\n".join(sorted(set(_stems(genomes))))
    return hashlib.sha256(joined.encode("utf-8")).hexdigest()


def write_genome_set_stamp(directory: Path, genomes: Sequence[Any], **extra: Any) -> Path:
    """Record which genome set produced a cohort-wide output folder."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    payload = {
        "fingerprint": genome_set_fingerprint(genomes),
        "samples": sorted(set(_stems(genomes))),
        **extra,
    }
    path = directory / GENOME_SET_STAMP
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    tmp.replace(path)
    return path


def read_genome_set_stamp(directory: Path) -> dict[str, Any] | None:
    return read_json_dict(Path(directory) / GENOME_SET_STAMP)


def genome_set_stamp_matches(directory: Path, genomes: Sequence[Any]) -> bool | None:
    """True/False when the folder has a stamp for this/another genome set; None without one."""
    stamp = read_genome_set_stamp(directory)
    if stamp is None:
        return None
    return str(stamp.get("fingerprint") or "") == genome_set_fingerprint(genomes)


def _stamp_errors(
    directory: Path,
    genomes: Sequence[Any],
    label: str,
    *,
    tool_input_matches: Callable[[], bool] | None = None,
) -> tuple[list[str], bool]:
    """(errors, allow_missing) from a cohort folder's genome-set stamp.

    A stamp for another genome set is an error. A matching stamp means this run's
    tool produced the output, so samples it left out are accepted. A folder without
    a stamp (written before stamps existed) accepts left-out samples only when the
    tool's own input record shows it was given this genome set (tool_input_matches);
    otherwise every genome must be present.
    """
    match = genome_set_stamp_matches(directory, genomes)
    if match is False:
        return [f"{label} was built for a different genome set ({Path(directory) / GENOME_SET_STAMP})"], False
    if match is None and tool_input_matches is not None:
        return [], bool(tool_input_matches())
    return [], bool(match)


def _gubbins_input_matches(output_dir: Path, genomes: Sequence[Any]) -> bool:
    """The alignment staged for Gubbins holds exactly this genome set (runs made before
    the ordered core alignment staged a copy in gubbins/input)."""
    path = Path(output_dir) / "gubbins" / "input" / "core_gene_alignment.aln"
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            names = {line[1:].strip() for line in handle if line.startswith(">")}
    except OSError:
        return False
    return names == set(_stems(genomes))


def _chewbbaca_input_count_matches(output_dir: Path, genomes: Sequence[Any]) -> bool:
    """chewBBACA's logging_info.txt reports one input per genome of this set."""
    path = Path(output_dir) / "chewbbaca" / "out" / "logging_info.txt"
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return False
    found = re.search(r"^Number of inputs:\s*(\d+)\s*$", text, flags=re.MULTILINE)
    return bool(found) and int(found.group(1)) == len(set(_stems(genomes)))


MARKER_HITS_STAMP = "query_fasta.json"


def write_marker_hits_stamp(
    hits_dir: Path,
    *,
    query_sha256: str,
    genomes: Sequence[Any],
    marker_sequences: Sequence[Any],
) -> Path:
    """Record the annotation FASTA (AnnotationIDs and sequences) the marker hits were searched from.

    The hit TSVs are keyed by AnnotationID, which numbers every CDS of the cohort, so
    they are only valid for the exact query FASTA that produced them.
    """
    hits_dir = Path(hits_dir)
    hits_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "query_fasta_sha256": str(query_sha256),
        "genome_set_fingerprint": genome_set_fingerprint(genomes),
        "markers": sorted(p.stem for p in _normalize_path_values(marker_sequences)),
    }
    path = hits_dir / MARKER_HITS_STAMP
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    tmp.replace(path)
    return path


def marker_hits_stamp_matches(hits_dir: Path, query_sha256: str) -> bool:
    """True only when the hits were searched from this exact annotation FASTA.

    A folder without the stamp (written before it existed, possibly for another
    genome set) is treated as stale.
    """
    stamp = read_json_dict(Path(hits_dir) / MARKER_HITS_STAMP)
    return bool(stamp) and str(stamp.get("query_fasta_sha256") or "") == str(query_sha256)


def chewbbaca_samples_left_out(output_dir: Path, genomes: Sequence[Any]) -> list[str]:
    out_dir = ensure_chewbbaca_out(Path(output_dir), logger=logging.getLogger(__name__))
    if out_dir is None:
        return []
    try:
        found = set(chewbbaca_sample_names(out_dir / "results_alleles.tsv"))
    except OSError:
        return []
    return sorted(set(_stems(genomes)) - found)


def _snippy_reference_stem(output_dir: Path) -> str | None:
    try:
        text = (Path(output_dir) / "snippy" / "reference.txt").read_text(encoding="utf-8").strip()
    except OSError:
        return None
    return Path(text).stem if text else None


def file_has_bytes(path: Path) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def file_exists(path: Path) -> bool:
    try:
        return path.is_file()
    except OSError:
        return False


def read_json_dict(path: Path) -> dict[str, Any] | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    return payload if isinstance(payload, dict) else None


def expected_bakta_json_paths(output_dir: Path, genomes: Sequence[Any]) -> list[Path]:
    return [
        Path(output_dir) / "bakta" / stem / f"{stem}.json"
        for stem in _stems(genomes)
    ]


def expected_bakta_gff_paths(output_dir: Path, genomes: Sequence[Any]) -> list[Path]:
    return [
        Path(output_dir) / "bakta" / stem / f"{stem}.gff3"
        for stem in _stems(genomes)
    ]


def expected_marker_db_paths(output_dir: Path, marker_sequences: Sequence[Any]) -> list[Path]:
    return [
        Path(output_dir) / "dmnd" / "queries" / f"{path.stem}.dmnd"
        for path in _normalize_path_values(marker_sequences)
    ]


def expected_marker_hit_paths(output_dir: Path, marker_sequences: Sequence[Any]) -> list[Path]:
    return [
        Path(output_dir) / "annotations" / "bakta_vs_queries" / f"bakta_annotations_nt-vs-{path.stem}.tsv"
        for path in _normalize_path_values(marker_sequences)
    ]


def expected_reference_hit_paths(
    output_dir: Path,
    genomes: Sequence[Any],
    reference_locus_sequences: Sequence[Any],
) -> list[Path]:
    genome_paths = _normalize_path_values(genomes)
    locus_paths = _normalize_path_values(reference_locus_sequences)
    return [
        Path(output_dir) / "annotations" / "reference_loci" / f"{genome.stem}-vs-{locus.stem}.tsv"
        for genome in genome_paths
        for locus in locus_paths
    ]


def expected_mlst_paths(output_dir: Path) -> list[Path]:
    return [Path(output_dir) / "mlst" / name for name in MLST_OUT_FILES]


def expected_assembly_qc_paths(output_dir: Path) -> list[Path]:
    return [Path(output_dir) / "assembly_qc" / name for name in ASSEMBLY_QC_OUT_FILES]


def expected_mobsuite_paths(output_dir: Path) -> list[Path]:
    return [Path(output_dir) / "mobsuite" / name for name in MOBSUITE_OUT_FILES]


def expected_chewbbaca_paths(output_dir: Path) -> list[Path]:
    return [Path(output_dir) / "chewbbaca" / "out" / name for name in CHEWBBACA_OUT_FILES]


def expected_gubbins_paths(output_dir: Path) -> list[Path]:
    """Gubbins outputs from its summary: one cohort run (old runs and single mode) or
    one folder per lineage that ran, the masked core alignment and the block tables; a
    step where no group ran has only its summary."""
    return [Path(output_dir) / rel for rel in expected_gubbins_files(read_gubbins_summary(output_dir))]


def snippy_tree_path(output_dir: Path) -> Path:
    """The Snippy tree: trees/raw/snippy_core.raxml.nwk, or the old trees/snippy_raxml.tree."""
    root = Path(output_dir)
    legacy = root / LEGACY_SNIPPY_TREE_REL
    raw = root / SNIPPY_RAW_REL
    return legacy if legacy.is_file() and not raw.is_file() else raw


def expected_phylogeny_paths(output_dir: Path) -> list[Path]:
    return [snippy_tree_path(output_dir)]


def expected_core_alignment_paths(output_dir: Path) -> list[Path]:
    return [Path(output_dir) / CORE_ALIGNMENT_DIR / rel for rel in CORE_ALIGNMENT_OUT_FILES]


def expected_core_tree_paths(output_dir: Path) -> list[Path]:
    root = Path(output_dir)
    stamp = read_core_tree_stamp(root) or {}
    tree_id = str((stamp.get("params") or {}).get("tree_id") or "core_recomb")
    return [
        root / IQTREE_DIR / f"{tree_id}.treefile",
        root / IQTREE_DIR / "run.json",
        tree_stamp_path(root, tree_id),
        root / "trees" / f"{tree_id}.nwk",
        root / "trees" / TREES_JSON,
    ]


def validate_bakta_outputs(output_dir: Path, genomes: Sequence[Any]) -> list[str]:
    errors: list[str] = []
    for json_path, gff_path in zip(
        expected_bakta_json_paths(output_dir, genomes),
        expected_bakta_gff_paths(output_dir, genomes),
    ):
        if not file_has_bytes(json_path):
            errors.append(f"Missing or empty Bakta JSON: {json_path}")
            continue
        if read_json_dict(json_path) is None:
            errors.append(f"Invalid Bakta JSON: {json_path}")
        if not file_has_bytes(gff_path):
            errors.append(f"Missing or empty Bakta GFF3: {gff_path}")
    return errors


def validate_assembly_qc_outputs(output_dir: Path, genomes: Sequence[Any]) -> list[str]:
    errors: list[str] = []
    paths = expected_assembly_qc_paths(output_dir)
    for path in paths:
        if not file_has_bytes(path):
            errors.append(f"Missing or empty assembly QC output: {path}")
    if errors:
        return errors

    expected_samples = set(_stems(genomes))
    summary_tsv = Path(output_dir) / "assembly_qc" / "summary.tsv"
    summary_json = Path(output_dir) / "assembly_qc" / "summary.json"
    summary = read_json_dict(summary_json)
    if summary is None:
        errors.append(f"Invalid assembly QC summary JSON: {summary_json}")

    try:
        with summary_tsv.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except Exception as exc:
        return [f"Failed to parse assembly QC summary TSV: {exc}"]

    if not rows:
        errors.append(f"Assembly QC summary has no rows: {summary_tsv}")
        return errors

    observed_samples = {str(row.get("Sample") or "").strip() for row in rows if str(row.get("Sample") or "").strip()}
    missing = sorted(sample for sample in expected_samples if sample not in observed_samples)
    if missing:
        errors.append(f"Assembly QC results missing samples: {', '.join(missing)}")

    for row in rows:
        sample = str(row.get("Sample") or "").strip() or "<missing sample>"
        for col in ("Quality", "CheckM2Completeness", "CheckM2Contamination", "QuastN50", "QuastContigs"):
            text = str(row.get(col) or "").strip()
            try:
                float(text)
            except Exception:
                errors.append(f"Assembly QC {col} is missing or non-numeric for {sample}.")
    return errors


def validate_mlst_outputs(output_dir: Path, genomes: Sequence[Any]) -> list[str]:
    errors: list[str] = []
    paths = expected_mlst_paths(output_dir)
    for path in paths:
        if not file_has_bytes(path):
            errors.append(f"Missing or empty MLST output: {path}")
    if errors:
        return errors

    results_json = next((path for path in paths if path.name == "results.json"), None)
    if results_json is None:
        return ["MLST results.json is missing."]
    try:
        data = json.loads(results_json.read_text(encoding="utf-8"))
    except Exception as exc:
        return [f"Failed to parse MLST results.json: {exc}"]
    if not isinstance(data, list) or not data:
        return ["MLST results.json is empty or malformed."]

    expected_samples = set(_stems(genomes))
    observed = {
        str(row.get("sample") or "").strip()
        for row in data
        if isinstance(row, Mapping)
    }
    missing = sorted(sample for sample in expected_samples if sample not in observed)
    if missing:
        errors.append(f"MLST results missing samples: {', '.join(missing)}")
    return errors


def validate_mobsuite_outputs(output_dir: Path, genomes: Sequence[Any]) -> list[str]:
    errors: list[str] = []
    contigs_tsv, plasmids_tsv, summary_json = expected_mobsuite_paths(output_dir)
    if not file_has_bytes(contigs_tsv):
        errors.append(f"Missing or empty MOB-suite contig table: {contigs_tsv}")
    if not file_has_bytes(plasmids_tsv):
        errors.append(f"Missing or empty MOB-suite plasmid table: {plasmids_tsv}")
    if not file_has_bytes(summary_json):
        errors.append(f"Missing or empty MOB-suite summary: {summary_json}")
        return errors

    summary = read_json_dict(summary_json)
    if summary is None:
        errors.append(f"Invalid MOB-suite summary JSON: {summary_json}")
        return errors

    samples = summary.get("samples")
    if not isinstance(samples, list):
        errors.append(f"Invalid MOB-suite summary samples list: {summary_json}")
        return errors
    expected_samples = set(_stems(genomes))
    observed_samples = {
        str(item.get("sample") or "").strip()
        for item in samples
        if isinstance(item, Mapping)
    }
    missing = sorted(sample for sample in expected_samples if sample not in observed_samples)
    if missing:
        errors.append(f"MOB-suite results missing samples: {', '.join(missing)}")
    failed = [
        str(item.get("sample") or "").strip()
        for item in samples
        if isinstance(item, Mapping) and str(item.get("status") or "").strip().lower() != "success"
    ]
    if failed:
        errors.append(f"MOB-suite samples failed: {', '.join(sorted(x for x in failed if x))}")
    return errors


def validate_panaroo_outputs(output_dir: Path, genomes: Sequence[Any] | None = None) -> list[str]:
    errors: list[str] = []
    for name in ("gene_presence_absence.csv", "core_gene_alignment.aln"):
        path = Path(output_dir) / "panaroo" / "out" / name
        if not file_has_bytes(path):
            errors.append(f"Missing or empty Panaroo output: {path}")
    if not errors and genomes:
        errors.extend(
            _sample_set_errors(panaroo_sample_names(output_dir), _stems(genomes), "Panaroo gene_presence_absence.csv")
        )
    return errors


def validate_core_alignment_outputs(output_dir: Path, genomes: Sequence[Any] | None = None) -> list[str]:
    errors: list[str] = []
    for path in expected_core_alignment_paths(output_dir):
        if not file_has_bytes(path):
            errors.append(f"Missing or empty core alignment output: {path}")
    if errors:
        return errors
    summary = read_ordering_summary(output_dir)
    if summary is None:
        return [f"Invalid core alignment record: {Path(output_dir) / CORE_ALIGNMENT_DIR / 'ordering.json'}"]
    if genomes:
        alignment = Path(output_dir) / ORDERED_ALIGNMENT_REL
        with alignment.open("rb") as handle:
            names = [line[1:].strip().decode("utf-8", "replace") for line in handle if line.startswith(b">")]
        errors.extend(_sample_set_errors(names, _stems(genomes), f"Ordered core alignment {alignment.name}"))
    return errors


def validate_gubbins_outputs(output_dir: Path, genomes: Sequence[Any] | None = None) -> list[str]:
    errors: list[str] = []
    allow_missing = False
    summary = read_gubbins_summary(output_dir)
    if summary is not None and summary.get("ran") is False:
        # No group ran (too few genomes, lineages too small): the summary says so.
        return [] if summary.get("skip_reason") else [f"Gubbins skip record without a reason: {output_dir}"]
    per_group = summary is not None and summary.get("schema_version") == GUBBINS_SCHEMA_VERSION
    if genomes:
        errors, allow_missing = _stamp_errors(
            Path(output_dir) / "gubbins",
            genomes,
            "Gubbins output",
            tool_input_matches=lambda: _gubbins_input_matches(output_dir, genomes),
        )
    for path in expected_gubbins_paths(output_dir):
        if not file_has_bytes(path):
            errors.append(f"Missing or empty Gubbins output: {path}")
            continue
        if path.name == "summary.json" and read_json_dict(path) is None:
            errors.append(f"Invalid Gubbins summary JSON: {path}")
        if path.suffix in {".tre", ".nwk", ".tree"}:
            try:
                text = path.read_text(encoding="utf-8", errors="ignore").strip()
            except OSError as exc:
                errors.append(f"Failed reading Gubbins tree: {path} ({exc})")
                continue
            if not text.endswith(";"):
                errors.append(f"Invalid Gubbins Newick tree output: {path}")
            elif genomes:
                # A lineage's tree holds that lineage only; Gubbins may leave gappy genomes out.
                errors.extend(
                    _tree_sample_errors(path, genomes, "Gubbins tree", allow_missing=allow_missing or per_group)
                )
    if per_group and genomes and not errors:
        masked = Path(output_dir) / MASKED_ALIGNMENT_REL
        with masked.open("rb") as handle:
            names = [line[1:].strip().decode("utf-8", "replace") for line in handle if line.startswith(b">")]
        errors.extend(_sample_set_errors(names, _stems(genomes), f"Masked core alignment {masked.name}"))
    return errors


def validate_core_tree_outputs(output_dir: Path, genomes: Sequence[Any] | None = None) -> list[str]:
    errors: list[str] = []
    for path in expected_core_tree_paths(output_dir):
        if not file_has_bytes(path):
            errors.append(f"Missing or empty core tree output: {path}")
    if errors:
        return errors
    index = read_json_dict(Path(output_dir) / "trees" / TREES_JSON)
    if index is None or not isinstance(index.get("trees"), list):
        return [f"Invalid tree index: {Path(output_dir) / 'trees' / TREES_JSON}"]
    tree_path = expected_core_tree_paths(output_dir)[3]
    entry = next((e for e in index["trees"] if isinstance(e, Mapping) and e.get("file") == tree_path.name), None)
    if entry is None:
        errors.append(f"The core tree {tree_path.name} is not listed in {TREES_JSON}")
    if genomes:
        errors.extend(_tree_sample_errors(tree_path, genomes, "Core tree", allow_missing=True))
    # A current refit (branch lengths including recombination) is listed with its file; a
    # skipped or failed one is a record only, and one kept from an older tree is not listed.
    core_raw = Path(output_dir) / IQTREE_DIR / "core_recomb.treefile"
    if entry is not None and current_refit_stamp(output_dir, core_raw) is not None:
        refit_path = Path(output_dir) / "trees" / f"{REFIT_TREE_ID}.nwk"
        for path in (Path(output_dir) / REFIT_DIR / f"{REFIT_TREE_ID}.treefile", refit_path):
            if not file_has_bytes(path):
                errors.append(f"Missing or empty refit core tree output: {path}")
        sets = entry.get("branch_length_sets") or []
        if not any(isinstance(item, Mapping) and item.get("file") == refit_path.name for item in sets):
            errors.append(f"The refit core tree {refit_path.name} is not listed in {TREES_JSON}")
    return errors


def validate_phylogeny_outputs(output_dir: Path, genomes: Sequence[Any] | None = None) -> list[str]:
    """The Snippy tree (the core genome tree is the core_tree step since 2026-10)."""
    errors: list[str] = []
    for path in expected_phylogeny_paths(output_dir):
        if not file_has_bytes(path):
            errors.append(f"Missing or empty phylogeny output: {path}")
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="ignore").strip()
        except OSError as exc:
            errors.append(f"Failed reading phylogeny tree: {path} ({exc})")
            continue
        if not text.endswith(";"):
            errors.append(f"Invalid Newick tree output: {path}")
        elif genomes:
            errors.extend(
                _tree_sample_errors(
                    path,
                    genomes,
                    "Phylogeny tree",
                    reference_stem=_snippy_reference_stem(output_dir),
                )
            )
    return errors


def population_output_files(output_dir: Path, method: str) -> dict[str, Path]:
    """The outputs of the population structure step (fastgear_clades) for a method."""
    folder = method_dir(output_dir, method)
    if method == "fastbaps":
        return {"sample_clades": folder / "sample_clades.tsv", "fastbaps_clusters": folder / "fastbaps_clusters.tsv"}
    return {
        "sample_clades": folder / "sample_clades.tsv",
        "lineage_information": folder / "output" / "lineage_information.txt",
        "result_mat": folder / "fastgear_res.mat",
    }


def validate_fastgear_outputs(output_dir: Path, genomes: Sequence[Any], method: str | None = "") -> list[str]:
    """The population structure step's outputs, of the run's recorded method (or of
    ``method``): fastGEAR in fastgear/ (every run before FastBAPS), FastBAPS in fastbaps/."""
    method = run_population_method(output_dir) if method == "" else method
    if method is None:
        return ["The run has no population structure method (neither fastGEAR nor FastBAPS was available)."]
    label = method_label(method)
    files = population_output_files(output_dir, method)

    errors: list[str] = []
    for path in files.values():
        if not file_has_bytes(path):
            errors.append(f"Missing or empty {label} output: {path}")
    if errors:
        return errors

    try:
        validate_sample_clades_tsv(files["sample_clades"], expected_samples=_stems(genomes))
    except Exception as exc:
        errors.append(f"Invalid {label} sample_clades.tsv: {exc}")
    return errors


def validate_chewbbaca_outputs(output_dir: Path, genomes: Sequence[Any] | None = None) -> list[str]:
    errors: list[str] = []
    for path in expected_chewbbaca_paths(output_dir):
        if not file_has_bytes(path):
            errors.append(f"Missing or empty chewBBACA output: {path}")
    if errors:
        return errors

    out_dir = ensure_chewbbaca_out(Path(output_dir), logger=logging.getLogger(__name__))
    if out_dir is None:
        return ["chewBBACA output directory could not be resolved."]
    alleles_path = out_dir / "results_alleles.tsv"
    if not file_has_bytes(alleles_path):
        return [f"Missing or empty chewBBACA alleles table: {alleles_path}"]
    try:
        with alleles_path.open("r", encoding="utf-8", errors="ignore") as handle:
            line_count = sum(1 for _ in handle)
    except OSError as exc:
        return [f"Failed to read chewBBACA alleles table: {exc}"]
    if line_count <= 1:
        errors.append(f"chewBBACA alleles table has no data rows: {alleles_path}")
    elif genomes:
        try:
            found = chewbbaca_sample_names(alleles_path)
        except OSError as exc:
            return [f"Failed to read chewBBACA alleles table: {exc}"]
        stamp_errors, allow_missing = _stamp_errors(
            Path(output_dir) / "chewbbaca",
            genomes,
            "chewBBACA output",
            tool_input_matches=lambda: _chewbbaca_input_count_matches(output_dir, genomes),
        )
        errors.extend(stamp_errors)
        errors.extend(
            _sample_set_errors(
                found,
                _stems(genomes),
                "chewBBACA results_alleles.tsv",
                allow_missing=allow_missing,
            )
        )
    return errors


def validate_amrfinder_outputs(output_dir: Path, genomes: Sequence[Any]) -> list[str]:
    """summary.json, and one AMRFinderPlus report per genome unless the step was skipped.

    Genomes the summary lists as failed have no report (the step goes on without them).
    """
    summary = read_amrfinder_summary(output_dir)
    if summary is None:
        return [f"Missing or invalid AMRFinderPlus summary: {Path(output_dir) / AMRFINDER_DIR / 'summary.json'}"]
    if summary.get("skipped") or not summary.get("samples"):
        return []
    errors = _sample_set_errors(
        [str(s) for s in summary.get("samples") or []],
        _stems(genomes),
        "AMRFinderPlus summary",
        allow_missing=False,
    )
    failed = summary.get("failed") if isinstance(summary.get("failed"), Mapping) else {}
    for stem in _stems(genomes):
        if stem in failed:
            continue
        path = amrfinder_tsv_path(output_dir, stem)
        try:
            header = path.open("r", encoding="utf-8").readline()
        except OSError:
            errors.append(f"Missing AMRFinderPlus report: {path}")
            continue
        if "Element symbol" not in header or "Subtype" not in header:
            errors.append(f"Unexpected AMRFinderPlus report header: {path}")
    return errors


def validate_marker_db_outputs(output_dir: Path, marker_sequences: Sequence[Any]) -> list[str]:
    errors: list[str] = []
    for path in expected_marker_db_paths(output_dir, marker_sequences):
        if not file_has_bytes(path):
            errors.append(f"Missing or empty DIAMOND marker DB: {path}")
    return errors


def validate_marker_hit_outputs(output_dir: Path, marker_sequences: Sequence[Any]) -> list[str]:
    errors: list[str] = []
    for path in expected_marker_hit_paths(output_dir, marker_sequences):
        if not file_exists(path):
            errors.append(f"Missing DIAMOND marker hit TSV: {path}")
    return errors


def validate_reference_hit_outputs(
    output_dir: Path,
    genomes: Sequence[Any],
    reference_locus_sequences: Sequence[Any],
) -> list[str]:
    errors: list[str] = []
    for path in expected_reference_hit_paths(output_dir, genomes, reference_locus_sequences):
        if not file_exists(path):
            errors.append(f"Missing DIAMOND reference-locus hit TSV: {path}")
    return errors


def _validator_inputs(run_cfg: Mapping[str, Any] | None) -> tuple[list[Any], list[Any], list[Any]]:
    if not isinstance(run_cfg, Mapping):
        return [], [], []
    inputs = run_cfg.get("inputs")
    if not isinstance(inputs, Mapping):
        return [], [], []
    genomes = list(inputs.get("genomes") or [])
    marker_sequences = list(inputs.get("marker_sequences") or [])
    reference_locus_sequences = list(inputs.get("reference_locus_sequences") or [])
    return genomes, marker_sequences, reference_locus_sequences


def validate_step_outputs(
    step_name: str,
    output_dir: Path,
    *,
    run_cfg: Mapping[str, Any] | None,
) -> list[str]:
    genomes, marker_sequences, reference_locus_sequences = _validator_inputs(run_cfg)
    validators = {
        STEP_ASSEMBLY_QC: lambda: validate_assembly_qc_outputs(output_dir, genomes),
        STEP_BAKTA: lambda: validate_bakta_outputs(output_dir, genomes),
        STEP_MLST: lambda: validate_mlst_outputs(output_dir, genomes),
        STEP_MOBSUITE: lambda: validate_mobsuite_outputs(output_dir, genomes),
        STEP_PANAROO: lambda: validate_panaroo_outputs(output_dir, genomes),
        STEP_CORE_ALIGNMENT: lambda: validate_core_alignment_outputs(output_dir, genomes),
        STEP_GUBBINS: lambda: validate_gubbins_outputs(output_dir, genomes),
        STEP_CORE_TREE: lambda: validate_core_tree_outputs(output_dir, genomes),
        STEP_PHYLOGENY: lambda: validate_phylogeny_outputs(output_dir, genomes),
        STEP_FASTGEAR: lambda: validate_fastgear_outputs(output_dir, genomes),
        STEP_CHEWBBACA: lambda: validate_chewbbaca_outputs(output_dir, genomes),
        STEP_DIAMOND_MARKERS: lambda: validate_marker_db_outputs(output_dir, marker_sequences),
        STEP_DIAMOND_REFERENCE: lambda: validate_reference_hit_outputs(
            output_dir,
            genomes,
            reference_locus_sequences,
        ),
        STEP_COLLECT_ANNOTATIONS: lambda: validate_marker_hit_outputs(output_dir, marker_sequences),
        STEP_AMRFINDER_MUTATIONS: lambda: validate_amrfinder_outputs(output_dir, genomes),
        STEP_BAKTA_EVIDENCE: lambda: validate_evidence_outputs(output_dir, genomes),
    }
    validator = validators.get(str(step_name))
    if validator is None:
        return []
    return validator()


def skippable_steps(run_cfg: Mapping[str, Any] | None) -> frozenset[str]:
    """Required steps a run skips without failing: the marker search when no marker database
    was installed, the reference locus search when no compound database was, the population
    structure when the tools image had neither fastGEAR nor FastBAPS."""
    _genomes, marker_sequences, reference_locus_sequences = _validator_inputs(run_cfg)
    out: set[str] = set()
    inputs = run_cfg.get("inputs") if isinstance(run_cfg, Mapping) else None
    population = inputs.get("population_structure") if isinstance(inputs, Mapping) else None
    if isinstance(population, Mapping) and population.get("method") is None:
        out.add(STEP_FASTGEAR)
    if not marker_sequences:
        out.add(STEP_DIAMOND_MARKERS)
    if not reference_locus_sequences:
        out.add(STEP_DIAMOND_REFERENCE)
    return frozenset(out)


def validate_completed_run(
    output_dir: Path,
    *,
    run_cfg: Mapping[str, Any] | None,
    steps: Mapping[str, Any] | None,
) -> tuple[dict[str, dict[str, Any]], list[str]]:
    validations: dict[str, dict[str, Any]] = {}
    all_errors: list[str] = []
    steps = steps if isinstance(steps, Mapping) else {}

    # An optional step counts once the run recorded it, unless it was skipped (an error in
    # it does not stop the run; the outputs then simply lack its data).
    optional = tuple(
        name for name in RESULT_OPTIONAL_STEPS
        if isinstance(steps.get(name), Mapping) and str(steps[name].get("status") or "").lower() != "skipped"
    )
    skippable = skippable_steps(run_cfg)
    for step_name in (*RESULT_REQUIRED_STEPS, *optional):
        step_info = steps.get(step_name)
        status = str(step_info.get("status") or "").strip().lower() if isinstance(step_info, Mapping) else ""
        step_errors: list[str]
        if status == "skipped" and step_name in skippable:
            step_errors = []
        elif status != "success":
            step_errors = [f"Step not successful: {step_name} ({status or 'missing'})"]
        else:
            step_errors = validate_step_outputs(step_name, output_dir, run_cfg=run_cfg)
        validations[step_name] = {
            "status": status or "missing",
            "outputs_validated": not step_errors,
            "validation_errors": step_errors,
        }
        all_errors.extend(step_errors)

    return validations, all_errors


__all__ = [
    "RESULT_OPTIONAL_STEPS",
    "RESULT_REQUIRED_STEPS",
    "STEP_AMRFINDER_MUTATIONS",
    "STEP_BAKTA_EVIDENCE",
    "STEP_ASSEMBLY_QC",
    "STEP_BAKTA",
    "STEP_MLST",
    "STEP_MOBSUITE",
    "STEP_PANAROO",
    "STEP_CORE_ALIGNMENT",
    "STEP_GUBBINS",
    "STEP_CORE_TREE",
    "STEP_PHYLOGENY",
    "STEP_FASTGEAR",
    "STEP_CHEWBBACA",
    "STEP_DIAMOND_MARKERS",
    "STEP_DIAMOND_REFERENCE",
    "STEP_COLLECT_ANNOTATIONS",
    "expected_marker_db_paths",
    "skippable_steps",
    "expected_marker_hit_paths",
    "expected_assembly_qc_paths",
    "expected_reference_hit_paths",
    "expected_mlst_paths",
    "expected_mobsuite_paths",
    "expected_chewbbaca_paths",
    "expected_core_alignment_paths",
    "expected_core_tree_paths",
    "expected_gubbins_paths",
    "expected_phylogeny_paths",
    "snippy_tree_path",
    "file_has_bytes",
    "GENOME_SET_STAMP",
    "MARKER_HITS_STAMP",
    "chewbbaca_samples_left_out",
    "file_sha256",
    "genome_set_change_error",
    "genome_set_fingerprint",
    "genome_set_stamp_matches",
    "marker_hits_stamp_matches",
    "read_genome_set_stamp",
    "tree_sample_names",
    "write_genome_set_stamp",
    "write_marker_hits_stamp",
    "newick_leaf_names",
    "recorded_genome_samples",
    "read_json_dict",
    "validate_amrfinder_outputs",
    "validate_assembly_qc_outputs",
    "validate_bakta_outputs",
    "validate_mlst_outputs",
    "validate_mobsuite_outputs",
    "validate_panaroo_outputs",
    "validate_core_alignment_outputs",
    "validate_core_tree_outputs",
    "validate_gubbins_outputs",
    "validate_phylogeny_outputs",
    "validate_fastgear_outputs",
    "population_output_files",
    "validate_chewbbaca_outputs",
    "validate_marker_db_outputs",
    "validate_marker_hit_outputs",
    "validate_reference_hit_outputs",
    "validate_step_outputs",
    "validate_completed_run",
]
