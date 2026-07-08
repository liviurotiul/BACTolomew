from __future__ import annotations

import json
import logging
import csv
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from pipeline.ChewBBACAWorkflow import CHEWBBACA_OUT_FILES, ensure_chewbbaca_out
from pipeline.AssemblyQualityWorkflow import ASSEMBLY_QC_OUT_FILES
from pipeline.FastGearWorkflow import validate_sample_clades_tsv
from pipeline.GubbinsWorkflow import GUBBINS_OUT_FILES
from pipeline.MlstWorkflow import MLST_OUT_FILES
from pipeline.MobSuiteWorkflow import MOBSUITE_OUT_FILES


STEP_BAKTA = "bakta"
STEP_ASSEMBLY_QC = "assembly_qc"
STEP_MLST = "mlst"
STEP_MOBSUITE = "mobsuite"
STEP_PANAROO = "panaroo"
STEP_GUBBINS = "gubbins"
STEP_PHYLOGENY = "phylogeny"
STEP_FASTGEAR = "fastgear_clades"
STEP_CHEWBBACA = "chewbbaca"
STEP_DIAMOND_MARKERS = "diamond_markers"
STEP_DIAMOND_REFERENCE = "diamond_reference_loci"
STEP_COLLECT_ANNOTATIONS = "collect_annotations"

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


def _normalize_path_values(values: Iterable[Any]) -> list[Path]:
    out: list[Path] = []
    for value in values:
        text = str(value or "").strip()
        if text:
            out.append(Path(text).expanduser())
    return out


def _stems(values: Iterable[Any]) -> list[str]:
    return [path.stem for path in _normalize_path_values(values)]


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
    return [Path(output_dir) / "gubbins" / rel for rel in GUBBINS_OUT_FILES]


def expected_phylogeny_paths(output_dir: Path) -> list[Path]:
    root = Path(output_dir)
    return [
        root / "trees" / "core_gene_tree.nwk",
        root / "trees" / "snippy_raxml.tree",
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


def validate_panaroo_outputs(output_dir: Path) -> list[str]:
    errors: list[str] = []
    for name in ("gene_presence_absence.csv", "core_gene_alignment.aln"):
        path = Path(output_dir) / "panaroo" / "out" / name
        if not file_has_bytes(path):
            errors.append(f"Missing or empty Panaroo output: {path}")
    return errors


def validate_gubbins_outputs(output_dir: Path) -> list[str]:
    errors: list[str] = []
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
    return errors


def validate_phylogeny_outputs(output_dir: Path) -> list[str]:
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
    return errors


def validate_fastgear_outputs(output_dir: Path, genomes: Sequence[Any]) -> list[str]:
    fastgear_dir = Path(output_dir) / "fastgear"
    sample_clades = fastgear_dir / "sample_clades.tsv"
    lineage_information = fastgear_dir / "output" / "lineage_information.txt"
    result_mat = fastgear_dir / "fastgear_res.mat"

    errors: list[str] = []
    for path in (sample_clades, lineage_information, result_mat):
        if not file_has_bytes(path):
            errors.append(f"Missing or empty fastGEAR output: {path}")
    if errors:
        return errors

    try:
        validate_sample_clades_tsv(sample_clades, expected_samples=_stems(genomes))
    except Exception as exc:
        errors.append(f"Invalid fastGEAR sample_clades.tsv: {exc}")
    return errors


def validate_chewbbaca_outputs(output_dir: Path) -> list[str]:
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
        STEP_PANAROO: lambda: validate_panaroo_outputs(output_dir),
        STEP_GUBBINS: lambda: validate_gubbins_outputs(output_dir),
        STEP_PHYLOGENY: lambda: validate_phylogeny_outputs(output_dir),
        STEP_FASTGEAR: lambda: validate_fastgear_outputs(output_dir, genomes),
        STEP_CHEWBBACA: lambda: validate_chewbbaca_outputs(output_dir),
        STEP_DIAMOND_MARKERS: lambda: validate_marker_db_outputs(output_dir, marker_sequences),
        STEP_DIAMOND_REFERENCE: lambda: validate_reference_hit_outputs(
            output_dir,
            genomes,
            reference_locus_sequences,
        ),
        STEP_COLLECT_ANNOTATIONS: lambda: validate_marker_hit_outputs(output_dir, marker_sequences),
    }
    validator = validators.get(str(step_name))
    if validator is None:
        return []
    return validator()


def validate_completed_run(
    output_dir: Path,
    *,
    run_cfg: Mapping[str, Any] | None,
    steps: Mapping[str, Any] | None,
) -> tuple[dict[str, dict[str, Any]], list[str]]:
    validations: dict[str, dict[str, Any]] = {}
    all_errors: list[str] = []
    steps = steps if isinstance(steps, Mapping) else {}

    for step_name in RESULT_REQUIRED_STEPS:
        step_info = steps.get(step_name)
        status = str(step_info.get("status") or "").strip().lower() if isinstance(step_info, Mapping) else ""
        step_errors: list[str]
        if status != "success":
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
    "RESULT_REQUIRED_STEPS",
    "STEP_ASSEMBLY_QC",
    "STEP_BAKTA",
    "STEP_MLST",
    "STEP_MOBSUITE",
    "STEP_PANAROO",
    "STEP_GUBBINS",
    "STEP_PHYLOGENY",
    "STEP_FASTGEAR",
    "STEP_CHEWBBACA",
    "STEP_DIAMOND_MARKERS",
    "STEP_DIAMOND_REFERENCE",
    "STEP_COLLECT_ANNOTATIONS",
    "expected_marker_db_paths",
    "expected_marker_hit_paths",
    "expected_assembly_qc_paths",
    "expected_reference_hit_paths",
    "expected_mlst_paths",
    "expected_mobsuite_paths",
    "expected_chewbbaca_paths",
    "expected_gubbins_paths",
    "expected_phylogeny_paths",
    "file_has_bytes",
    "read_json_dict",
    "validate_assembly_qc_outputs",
    "validate_bakta_outputs",
    "validate_mlst_outputs",
    "validate_mobsuite_outputs",
    "validate_panaroo_outputs",
    "validate_gubbins_outputs",
    "validate_phylogeny_outputs",
    "validate_fastgear_outputs",
    "validate_chewbbaca_outputs",
    "validate_marker_db_outputs",
    "validate_marker_hit_outputs",
    "validate_reference_hit_outputs",
    "validate_step_outputs",
    "validate_completed_run",
]
