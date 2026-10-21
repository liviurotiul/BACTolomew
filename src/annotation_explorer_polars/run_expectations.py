from __future__ import annotations

import json
import os
import socket
from pathlib import Path
from typing import Optional, Sequence

from pipeline.RunValidation import skippable_steps

from .constants import RUN_LOCK_FILENAME, RUN_REPORT_FILENAME, RUN_STATE_FILENAME

_FASTGEAR_EXPECTED_RELATIVE = (
    Path("fastgear") / "fastgear_res.mat",
    Path("fastgear") / "output" / "lineage_information.txt",
    Path("fastgear") / "sample_clades.tsv",
)
# Runs whose expected_outputs.json predates the tree files (the phylogeny rework of
# 2026-10 lists trees/trees.json and the raw Snippy tree itself).
_PHYLOGENY_EXPECTED_RELATIVE = (
    Path("trees") / "core_gene_tree.nwk",
    Path("trees") / "snippy_raxml.tree",
)
_GUBBINS_EXPECTED_RELATIVE = (
    Path("gubbins") / "core_gene.final_tree.tre",
    Path("gubbins") / "core_gene.filtered_polymorphic_sites.fasta",
    Path("gubbins") / "core_gene.recombination_predictions.gff",
    Path("gubbins") / "summary.json",
)
_CHEWBBACA_EXPECTED_RELATIVE = (
    Path("chewbbaca") / "out" / "results_alleles.tsv",
    Path("chewbbaca") / "out" / "results_statistics.tsv",
    Path("chewbbaca") / "out" / "loci_summary_stats.tsv",
    Path("chewbbaca") / "out" / "logging_info.txt",
)
_MLST_EXPECTED_RELATIVE = (
    Path("mlst") / "results.csv",
    Path("mlst") / "results.json",
)
_MOBSUITE_EXPECTED_RELATIVE = (
    Path("mobsuite") / "contigs.tsv",
    Path("mobsuite") / "plasmids.tsv",
    Path("mobsuite") / "summary.json",
)
_ASSEMBLY_QC_EXPECTED_RELATIVE = (
    Path("assembly_qc") / "summary.tsv",
    Path("assembly_qc") / "summary.json",
    Path("assembly_qc") / "quast" / "report.tsv",
    Path("assembly_qc") / "checkm2" / "quality_report.tsv",
)
_IGNORED_LEGACY_SEGMENTS = (
    "/snippy/",
    "/variant_calling/",
    "/platon/",
)

_PENDING_LABELS = {
    "Run files": "run files",
    "Bakta": "Bakta",
    "Core alignment": "ordered core alignment",
    "Gubbins": "Gubbins",
    "Phylogeny": "phylogeny trees",
    "DIAMOND marker DBs": "marker DBs",
    "DIAMOND marker hits": "marker hits",
    "DIAMOND loci DBs": "reference-locus DBs",
    "DIAMOND loci hits": "reference-locus hits",
    "Panaroo": "Panaroo",
    "fastGEAR": "fastGEAR clades",
    "chewBBACA": "chewBBACA",
    "MLST": "MLST",
    "MOB-suite": "MOB-suite",
    "Assembly QC": "assembly QC",
    "AMRFinderPlus": "AMRFinderPlus",
    "Bakta evidence": "Bakta evidence",
}
_REQUIRED_STEP_LABELS = {
    "assembly_qc": "assembly QC",
    "bakta": "Bakta",
    "mlst": "MLST",
    "mobsuite": "MOB-suite",
    "panaroo": "Panaroo",
    "gubbins": "Gubbins",
    "phylogeny": "phylogeny",
    "fastgear_clades": "fastGEAR clades",
    "chewbbaca": "chewBBACA",
    "diamond_markers": "marker hits",
    "diamond_reference_loci": "reference-locus hits",
    "collect_annotations": "annotation collection",
}


def _dedupe_paths(paths: Sequence[Path]) -> list[Path]:
    seen: set[str] = set()
    out: list[Path] = []
    for path in paths:
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        out.append(path)
    return out


def _load_run_config(output_dir: Path) -> dict | None:
    run_cfg_path = output_dir / "run_config.json"
    if not run_cfg_path.exists():
        return None
    return _load_json_dict(run_cfg_path)


def _load_json_dict(path: Path) -> dict | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    return payload if isinstance(payload, dict) else None


def _pid_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    return True


def run_is_active(output_dir: Path) -> bool:
    lock_path = Path(output_dir) / RUN_LOCK_FILENAME
    if not lock_path.exists():
        return False
    data = _load_json_dict(lock_path)
    if data is None:
        return True
    host = str(data.get("host") or "").strip()
    pid = int(data.get("pid") or 0)
    if host and host != socket.gethostname():
        return True
    if _pid_alive(pid):
        return True
    try:
        lock_path.unlink()
    except Exception:
        pass
    return False


def load_run_state(output_dir: Path) -> dict | None:
    return _load_json_dict(Path(output_dir) / RUN_STATE_FILENAME)


def load_run_report(output_dir: Path) -> dict | None:
    return _load_json_dict(Path(output_dir) / RUN_REPORT_FILENAME)


def _missing_expected_output_message(output_dir: Path) -> str | None:
    missing_groups: list[str] = []
    seen: set[str] = set()
    for path in resolve_expected_files(output_dir):
        if expected_output_exists(output_dir, path):
            continue
        group = classify_expected_path(path)
        if group is None or group in seen:
            continue
        seen.add(group)
        missing_groups.append(group)

    if not missing_groups:
        return None
    labels = [_PENDING_LABELS.get(group, group) for group in missing_groups]
    return "Missing required outputs: " + ", ".join(labels)


def _required_step_validation_errors(report: dict, run_cfg: dict | None) -> list[str]:
    steps = report.get("steps") if isinstance(report.get("steps"), dict) else {}
    missing: list[str] = []
    incomplete: list[str] = []
    # A search skipped because none of its databases was installed is not incomplete.
    skippable = skippable_steps(run_cfg)

    for step_name, label in _REQUIRED_STEP_LABELS.items():
        step = steps.get(step_name) if isinstance(steps, dict) else None
        if not isinstance(step, dict):
            missing.append(label)
            continue
        status = str(step.get("status") or "").strip().lower()
        if status == "skipped" and step_name in skippable:
            continue
        if status != "success" or step.get("outputs_validated") is False:
            incomplete.append(label)

    errors: list[str] = []
    if missing:
        errors.append("Missing required steps: " + ", ".join(missing))
    if incomplete:
        errors.append("Required steps not successful: " + ", ".join(incomplete))
    return errors


def load_run_status(output_dir: Path) -> dict[str, object]:
    output_dir = Path(output_dir)
    active = run_is_active(output_dir)
    state = load_run_state(output_dir) or {}
    report = load_run_report(output_dir) or {}
    has_final_report = bool(report)

    state_status = str(state.get("status") or "").strip().lower()
    report_status = str(report.get("status") or "").strip().lower()
    if active and (state_status == "running" or not state):
        raw_status = "running"
    else:
        raw_status = state_status or report_status
    if not raw_status:
        raw_status = "running" if active else "unknown"
    if raw_status == "running" and not active:
        raw_status = "interrupted"

    outputs_validated_raw = state.get("outputs_validated")
    if outputs_validated_raw is None:
        outputs_validated_raw = report.get("outputs_validated")
    outputs_validated = bool(outputs_validated_raw)

    validation_errors_raw = state.get("validation_errors")
    if validation_errors_raw is None:
        validation_errors_raw = report.get("validation_errors")
    validation_errors = (
        [str(item).strip() for item in validation_errors_raw if str(item).strip()]
        if isinstance(validation_errors_raw, list)
        else []
    )

    if raw_status == "success" and not has_final_report:
        raw_status = "partial"
        validation_errors = validation_errors or ["Missing final run report."]
    if raw_status == "success" and not outputs_validated:
        raw_status = "partial"
        validation_errors = validation_errors or ["Run outputs were not validated successfully."]

    if raw_status == "success":
        success_errors = _required_step_validation_errors(report, _load_run_config(output_dir))
        missing_outputs = _missing_expected_output_message(output_dir)
        if missing_outputs:
            success_errors.insert(0, missing_outputs)
        if success_errors:
            raw_status = "partial"
            outputs_validated = False
            validation_errors = validation_errors or success_errors

    error_text = str(state.get("error") or report.get("error") or "").strip()
    if not error_text and validation_errors:
        error_text = validation_errors[0]

    return {
        "status": raw_status,
        "active": active,
        "current_step": str(state.get("current_step") or "").strip(),
        "error": error_text,
        "run_id": str(state.get("run_id") or report.get("run_id") or "").strip(),
        "sample_set_id": str(state.get("sample_set_id") or report.get("sample_set_id") or "").strip(),
        "updated_at": str(state.get("updated_at") or report.get("end_time") or "").strip(),
        "outputs_validated": outputs_validated,
        "validation_errors": validation_errors,
        "has_final_report": has_final_report,
        # The parsed run_state.json, so a caller that also shows the steps reads it once.
        "state": state,
    }


def run_allows_view(status_info: dict[str, object]) -> bool:
    return (
        str(status_info.get("status") or "").strip().lower() == "success"
        and bool(status_info.get("outputs_validated"))
        and bool(status_info.get("has_final_report"))
    )


def run_allows_resume(status_info: dict[str, object]) -> bool:
    if bool(status_info.get("active")):
        return False
    status = str(status_info.get("status") or "").strip().lower()
    return status != "success"


def describe_run_status(status_info: dict[str, object], pending_text: str | None = None) -> Optional[str]:
    status = str(status_info.get("status") or "").strip().lower()
    current_step = str(status_info.get("current_step") or "").strip()
    error = str(status_info.get("error") or "").strip()

    if status == "running":
        if current_step:
            return f"Running: {current_step}"
        return "Running"
    if status == "success":
        return "Completed"
    if status == "stopped":
        return "Stopped by user. Click Resume to continue."
    if status == "interrupted":
        return "Interrupted. Click Resume to continue."
    if status == "failed":
        return f"Failed: {error}" if error else "Failed. Click Resume to continue."
    if status == "partial":
        if error:
            return f"Partial: {error}"
        return pending_text or "Partial run. Click Resume to continue."
    return pending_text or "Status unknown. Click Resume to continue."


def _run_uses_panaroo(output_dir: Path, expected_files: Sequence[Path]) -> bool:
    if any("/panaroo/out/" in str(path) for path in expected_files):
        return True
    panaroo_out = output_dir / "panaroo" / "out"
    return panaroo_out.exists()


def _should_expect_fastgear(output_dir: Path, expected_files: Sequence[Path]) -> bool:
    if any("/fastgear/" in str(path) for path in expected_files):
        return True
    run_cfg = _load_run_config(output_dir)
    if run_cfg is None:
        return False
    inputs = run_cfg.get("inputs")
    if not isinstance(inputs, dict):
        return False
    return _run_uses_panaroo(output_dir, expected_files)


def _should_expect_phylogeny(output_dir: Path, expected_files: Sequence[Path]) -> bool:
    # A run that lists its tree files is taken as is (new runs no longer make
    # core_gene_tree.nwk); only an old list without them gets the old names.
    if any("/trees/" in str(path) for path in expected_files):
        return False
    return _run_uses_panaroo(output_dir, expected_files)


def _should_expect_gubbins(output_dir: Path, expected_files: Sequence[Path]) -> bool:
    # A run that lists its Gubbins files is taken as is (a skipped Gubbins has only
    # its summary).
    if any("/gubbins/" in str(path) for path in expected_files):
        return False
    return _run_uses_panaroo(output_dir, expected_files)


def resolve_expected_files(output_dir: Path) -> list[Path]:
    expected_path = output_dir / "expected_outputs.json"
    if not expected_path.exists():
        return []
    try:
        data = json.loads(expected_path.read_text(encoding="utf-8"))
    except Exception:
        return []
    files = data.get("files")
    if not isinstance(files, list) or not files:
        return []

    resolved: list[Path] = []
    for raw in files:
        path = Path(raw)
        if not path.is_absolute():
            path = output_dir / path
        path_str = str(path)
        if any(segment in path_str for segment in _IGNORED_LEGACY_SEGMENTS):
            continue
        resolved.append(path)

    if _should_expect_gubbins(output_dir, resolved):
        resolved.extend(output_dir / rel for rel in _GUBBINS_EXPECTED_RELATIVE)
    if _should_expect_phylogeny(output_dir, resolved):
        resolved.extend(output_dir / rel for rel in _PHYLOGENY_EXPECTED_RELATIVE)
    if _should_expect_fastgear(output_dir, resolved):
        resolved.extend(output_dir / rel for rel in _FASTGEAR_EXPECTED_RELATIVE)
    resolved.extend(output_dir / rel for rel in _CHEWBBACA_EXPECTED_RELATIVE)
    resolved.extend(output_dir / rel for rel in _MLST_EXPECTED_RELATIVE)
    resolved.extend(output_dir / rel for rel in _MOBSUITE_EXPECTED_RELATIVE)
    resolved.extend(output_dir / rel for rel in _ASSEMBLY_QC_EXPECTED_RELATIVE)
    return _dedupe_paths(resolved)


def expected_output_exists(output_dir: Path, expected_path: Path) -> bool:
    if expected_path.exists():
        return True

    s = str(expected_path)
    if any(segment in s for segment in _IGNORED_LEGACY_SEGMENTS):
        return False
    return False


def classify_expected_path(path: Path) -> Optional[str]:
    name = path.name
    s = str(path)
    if name in {"run_config.json", "run_report.json", RUN_STATE_FILENAME}:
        return "Run files"
    if "/panaroo/out/" in s:
        return "Panaroo"
    if "/core_alignment/" in s:
        return "Core alignment"
    if "/gubbins/" in s:
        return "Gubbins"
    if "/trees/" in s or name in {"core_gene_tree.nwk", "snippy_raxml.tree", "snippy_raxml.newick"}:
        return "Phylogeny"
    if "/fastgear/" in s:
        return "fastGEAR"
    if "/chewbbaca/out/" in s:
        return "chewBBACA"
    if "/mlst/" in s:
        return "MLST"
    if "/mobsuite/" in s:
        return "MOB-suite"
    if "/assembly_qc/" in s:
        return "Assembly QC"
    if "/amrfinder/" in s:
        return "AMRFinderPlus"
    if "/evidence/" in s:
        return "Bakta evidence"
    if "/annotations/bakta_vs_queries/" in s:
        return "DIAMOND marker hits"
    if "/annotations/reference_loci/" in s:
        return "DIAMOND loci hits"
    if "/dmnd/queries/" in s:
        return "DIAMOND marker DBs"
    if "/dmnd/reference_loci/" in s:
        return "DIAMOND loci DBs"
    if "/bakta/" in s:
        return "Bakta"
    return None


def summarize_pending_outputs(output_dir: Path, expected_files: Sequence[Path]) -> Optional[str]:
    missing_groups: list[str] = []
    seen: set[str] = set()
    for path in expected_files:
        if expected_output_exists(output_dir, path):
            continue
        group = classify_expected_path(path)
        if group is None or group in seen:
            continue
        seen.add(group)
        missing_groups.append(group)

    if not missing_groups:
        return None
    labels = [_PENDING_LABELS.get(group, group) for group in missing_groups]
    if len(labels) == 1 and labels[0] == "fastGEAR clades":
        return "fastGEAR not yet run. Click Resume to run it."
    if len(labels) == 1 and labels[0] == "chewBBACA":
        return "chewBBACA not yet run. Click Resume to run it."
    if len(labels) == 1 and labels[0] == "MLST":
        return "MLST not yet run. Click Resume to run it."
    if len(labels) == 1 and labels[0] == "MOB-suite":
        return "MOB-suite not yet run. Click Resume to run it."
    if len(labels) == 1 and labels[0] == "assembly QC":
        return "Assembly QC not yet run. Click Resume to run it."
    if len(labels) == 1 and labels[0] == "Gubbins":
        return "Gubbins not yet run. Click Resume to run it."
    return "Pending: " + ", ".join(labels)


def pending_outputs_text(output_dir: Path) -> Optional[str]:
    """What a partial run still lacks ("Gubbins not yet run. Click Resume to run it.",
    "Pending: ..."), from expected_outputs.json; None when nothing is missing."""
    output_dir = Path(output_dir)
    expected_files = resolve_expected_files(output_dir)
    if not expected_files:
        return None
    return summarize_pending_outputs(output_dir, expected_files)


__all__ = [
    "describe_run_status",
    "classify_expected_path",
    "expected_output_exists",
    "load_run_report",
    "load_run_state",
    "load_run_status",
    "pending_outputs_text",
    "resolve_expected_files",
    "run_allows_resume",
    "run_allows_view",
    "run_is_active",
    "summarize_pending_outputs",
]
