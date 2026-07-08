from __future__ import annotations

import json
import hashlib
import subprocess
import tempfile
import os
import socket
from datetime import datetime
from pathlib import Path

import pandas as pd

from pipeline.AnnotationUtils import (
    assign_annotation_ids,
    compute_category,
    infer_annotation_category_and_source,
    write_annotation_fasta,
)
from pipeline.AssemblyQualityWorkflow import (
    ASSEMBLY_QC_OUT_FILES,
    run_assembly_quality as run_assembly_quality_tool,
)
from pipeline.BaktaWorkflow import parse_multiple_bakta_json, run_bakta_on_fastas
from pipeline.ChewBBACAWorkflow import (
    CHEWBBACA_OUT_FILES,
    ensure_chewbbaca_out,
    run_chewbbaca as run_chewbbaca_tool,
    validate_chewbbaca_settings,
)
from pipeline.DataInput import DataInput
from pipeline.DbCatalog import MARKER_DB_CATALOG_PATH, load_db_catalog_index
from pipeline.DiamondWorkflow import (
    build_diamond_db_for_files,
    collapse_diamond_annotations,
    load_diamond_tsvs,
    read_annotations,
    run_diamond_annotations,
)
from pipeline.FastGearWorkflow import run_fastgear_clades, validate_sample_clades_tsv
from pipeline.GubbinsWorkflow import (
    GUBBINS_OUT_FILES,
    run_gubbins as run_gubbins_tool,
)
from pipeline.Logger import setup_logger
from pipeline.MlstWorkflow import (
    MLST_OUT_FILES,
    run_mlst as run_mlst_tool,
)
from pipeline.MobSuiteWorkflow import (
    MOBSUITE_OUT_FILES,
    load_mobsuite_contig_predictions,
    run_mobsuite as run_mobsuite_tool,
)
from pipeline.ResistanceMetadata import (
    RESISTANCE_MARKER_COLUMNS,
    build_annotation_resistance_maps,
    enrich_marker_hits_with_resistance,
)
from pipeline.PhylogenyWorkflow import (
    run_panaroo as run_panaroo_tool,
    ensure_panaroo_out,
    run_phylogeny as run_phylogeny_tool,
)
from pipeline.RunValidation import (
    file_has_bytes as run_file_has_bytes,
    read_json_dict as run_read_json_dict,
    validate_bakta_outputs,
    validate_assembly_qc_outputs,
    validate_chewbbaca_outputs,
    validate_completed_run,
    validate_fastgear_outputs,
    validate_gubbins_outputs,
    validate_marker_db_outputs,
    validate_marker_hit_outputs,
    validate_mlst_outputs,
    validate_panaroo_outputs,
    validate_mobsuite_outputs,
    validate_phylogeny_outputs,
    validate_reference_hit_outputs,
)
from pipeline.ToolEnv import (
    conda_tool_cmd,
    configured_assembly_qc_settings,
    configured_chewbbaca_settings,
    configured_fastgear_settings,
    configured_mobsuite_settings,
    make_container_context,
    reset_container_context,
    set_container_context,
)
from pipeline.Utils import assign_gene_ids

logger = setup_logger()

RUN_CONFIG_FILENAME = "run_config.json"
RUN_REPORT_FILENAME = "run_report.json"
RUN_STATE_FILENAME = "run_state.json"
RUN_LOCK_FILENAME = "run.lock.json"
RUN_STOP_FILENAME = "run.stop"
ANNOTATION_ID_STRATEGY = "sorted_bakta_json_v1"

ANNOTATION_COLUMNS = [
    "Sample", "Contig", "Start", "Stop", "Strand",
    "Gene", "Product",
    "Identity", "Evalue", "Bitscore",
    "Method", "MarkerCategory", "NCBI_NRP_ID", "Refseq_ID", "Id", "AnnotationCategory",
    "Locus", "StartType", "RBS", "Protein_MW", "pI", "Hypothetical",
    "Nucleotide_Sequence", "Amino_Acid_Sequence",
    *RESISTANCE_MARKER_COLUMNS,
    "MobSuitePrediction",
    "MobSuiteIsPlasmid",
    "MobSuiteSourceContig",
    "MobSuitePlasmidID",
    "MobSuiteMoleculeType",
    "MobSuitePrimaryClusterID",
    "MobSuiteSecondaryClusterID",
    "MobSuiteLength",
    "MobSuiteRepliconTypes",
    "MobSuiteRelaxaseTypes",
    "MobSuiteMPFType",
    "MobSuiteOriTTypes",
    "MobSuitePredictedMobility",
    "MobSuiteMashNearestNeighbor",
    "MobSuiteMashNeighborDistance",
    "AnnotationID",
]

BAKTA_OUTPUT_FILES = [
    ".embl",
    ".faa",
    ".ffn",
    ".fna",
    ".gbff",
    ".gff3",
    ".hypotheticals.faa",
    ".hypotheticals.tsv",
    ".inference.tsv",
    ".json",
    ".log",
    ".tsv",
    ".txt",
]

PANAROO_OUT_FILES = [
    "gene_presence_absence.csv",
    "core_gene_alignment.aln",
]

MOBSUITE_OUTPUT_FILES = list(MOBSUITE_OUT_FILES)

FASTGEAR_OUT_FILES = [
    Path("fastgear_res.mat"),
    Path("output") / "lineage_information.txt",
    Path("sample_clades.tsv"),
]

GUBBINS_OUTPUT_FILES = list(GUBBINS_OUT_FILES)

PHYLOGENY_OUT_FILES = [
    Path("trees") / "core_gene_tree.nwk",
    Path("trees") / "snippy_raxml.tree",
]

BAKTA_REUSE_SUFFIXES = [
    ".json",
    ".gff3",
]


class RunStoppedError(RuntimeError):
    pass


def _compute_thread_plan(
    total_cores: int,
    *,
    n_genomes: int,
    n_marker_dbs: int,
    n_locus_dbs: int,
) -> dict[str, int]:
    total = max(1, int(total_cores))
    panaroo_threads = min(total, max(4, min(32, 2 * max(1, int(n_genomes)))))
    return {
        "total_cores": total,
        "bakta": min(total, 32),
        "assembly_qc": min(total, 32),
        "mobsuite": min(total, 32),
        "chewbbaca": min(total, 32),
        "mlst": min(total, max(1, min(64, int(n_genomes)))),
        "diamond_makedb": min(total, 32),
        "diamond_blastx": min(total, 64),
        "panaroo": panaroo_threads,
        "gubbins": min(total, 64),
        "phylogeny_snippy": min(total, 16),
        "phylogeny_raxml": min(total, 64),
        "n_genomes": int(n_genomes),
        "n_marker_dbs": int(n_marker_dbs),
        "n_locus_dbs": int(n_locus_dbs),
    }


def _ensure_annotation_columns(df: pd.DataFrame) -> pd.DataFrame:
    for col in ANNOTATION_COLUMNS:
        if col not in df.columns:
            df[col] = pd.NA
    return df


def _now_iso() -> str:
    return datetime.now().astimezone().isoformat()


def _compute_sample_set_id(genomes: list[Path]) -> str:
    samples = sorted({p.stem for p in genomes if p is not None})
    joined = "|".join(samples)
    digest = hashlib.sha1(joined.encode("utf-8")).hexdigest()[:12]
    return f"sample_{digest}"


def _new_attempt_id(sample_set_id: str) -> str:
    stamp = datetime.now().astimezone().strftime("%Y%m%d_%H%M%S_%f")
    return f"{sample_set_id}__run_{stamp}"


def _compute_run_id(genomes: list[Path]) -> str:
    return _compute_sample_set_id(genomes)


def _collect_output_files(root: Path) -> list[str]:
    if not root.exists():
        return []
    files = [p for p in root.rglob("*") if p.is_file()]
    return sorted(str(p) for p in files)


def _dedupe_paths(paths: list[Path]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for p in paths:
        s = str(p)
        if s in seen:
            continue
        seen.add(s)
        out.append(s)
    return out


def _all_files_exist(paths: list[Path]) -> bool:
    return bool(paths) and all(path.is_file() for path in paths)


def _read_json_dict(path: Path) -> dict | None:
    return run_read_json_dict(path)


def _file_has_bytes(path: Path) -> bool:
    return run_file_has_bytes(path)


def _validate_bakta_reuse_outputs(CFG: DataInput) -> bool:
    return not validate_bakta_outputs(CFG.output, CFG.genomes)


def _validate_assembly_qc_reuse_outputs(CFG: DataInput) -> bool:
    return not validate_assembly_qc_outputs(CFG.output, CFG.genomes)


def _validate_mlst_reuse_outputs(CFG: DataInput) -> bool:
    return not validate_mlst_outputs(CFG.output, CFG.genomes)


def _validate_mobsuite_reuse_outputs(CFG: DataInput) -> bool:
    return not validate_mobsuite_outputs(CFG.output, CFG.genomes)


def _validate_chewbbaca_reuse_outputs(CFG: DataInput) -> bool:
    return not validate_chewbbaca_outputs(CFG.output)


def _validate_gubbins_reuse_outputs(CFG: DataInput) -> bool:
    return not validate_gubbins_outputs(CFG.output)


def _validate_dmnd_files(paths: list[Path]) -> bool:
    return bool(paths) and all(_file_has_bytes(path) for path in paths)


def _expected_bakta_reuse_paths(CFG: DataInput) -> list[Path]:
    expected: list[Path] = []
    for genome in CFG.genomes:
        stem = genome.stem
        sample_dir = CFG.output / "bakta" / stem
        for suffix in BAKTA_REUSE_SUFFIXES:
            expected.append(sample_dir / f"{stem}{suffix}")
    return expected


def _expected_marker_db_paths(CFG: DataInput) -> list[Path]:
    return [CFG.output / "dmnd" / "queries" / f"{seq.stem}.dmnd" for seq in CFG.marker_sequences]


def _expected_assembly_qc_paths(CFG: DataInput) -> list[Path]:
    return [CFG.output / "assembly_qc" / rel for rel in ASSEMBLY_QC_OUT_FILES]


def _expected_marker_hit_paths(CFG: DataInput) -> list[Path]:
    return [
        CFG.output / "annotations" / "bakta_vs_queries" / f"bakta_annotations_nt-vs-{seq.stem}.tsv"
        for seq in CFG.marker_sequences
    ]


def _expected_reference_hit_paths(CFG: DataInput) -> list[Path]:
    return [
        CFG.output / "annotations" / "reference_loci" / f"{genome.stem}-vs-{seq.stem}.tsv"
        for genome in CFG.genomes
        for seq in CFG.reference_locus_sequences
    ]


def _expected_chewbbaca_paths(CFG: DataInput) -> list[Path]:
    return [CFG.output / "chewbbaca" / "out" / fname for fname in CHEWBBACA_OUT_FILES]


def _expected_mlst_paths(CFG: DataInput) -> list[Path]:
    return [CFG.output / "mlst" / fname for fname in MLST_OUT_FILES]


def _expected_mobsuite_paths(CFG: DataInput) -> list[Path]:
    return [CFG.output / "mobsuite" / fname for fname in MOBSUITE_OUTPUT_FILES]


def _expected_gubbins_paths(CFG: DataInput) -> list[Path]:
    return [CFG.output / "gubbins" / rel for rel in GUBBINS_OUTPUT_FILES]


def _expected_phylogeny_paths(CFG: DataInput) -> list[Path]:
    return [CFG.output / rel for rel in PHYLOGENY_OUT_FILES]


def _build_expected_outputs(CFG: DataInput, *, run_id: str) -> list[str]:
    output = CFG.output
    genomes = [Path(p) for p in CFG.genomes]

    expected: list[Path] = []

    expected.append(output / "expected_outputs.json")
    expected.append(output / RUN_CONFIG_FILENAME)
    expected.append(output / RUN_STATE_FILENAME)
    expected.append(output / RUN_REPORT_FILENAME)

    # Bakta outputs per genome
    for g in genomes:
        stem = g.stem
        sample_dir = output / "bakta" / stem
        for suffix in BAKTA_OUTPUT_FILES:
            expected.append(sample_dir / f"{stem}{suffix}")

    for path in _expected_assembly_qc_paths(CFG):
        expected.append(path)

    # DIAMOND DBs
    for seq in CFG.marker_sequences:
        expected.append(output / "dmnd" / "queries" / f"{seq.stem}.dmnd")
    for seq in CFG.reference_locus_sequences:
        expected.append(output / "dmnd" / "reference_loci" / f"{seq.stem}.dmnd")

    # DIAMOND annotations
    for seq in CFG.marker_sequences:
        expected.append(
            output / "annotations" / "bakta_vs_queries" / f"bakta_annotations_nt-vs-{seq.stem}.tsv"
        )
    for g in genomes:
        for seq in CFG.reference_locus_sequences:
            expected.append(
                output / "annotations" / "reference_loci" / f"{g.stem}-vs-{seq.stem}.tsv"
            )

    # Panaroo (minimal expected outputs)
    for fname in PANAROO_OUT_FILES:
        expected.append(output / "panaroo" / "out" / fname)

    for path in _expected_gubbins_paths(CFG):
        expected.append(path)

    for path in _expected_phylogeny_paths(CFG):
        expected.append(path)

    # fastGEAR clades from Panaroo core alignment
    for relpath in FASTGEAR_OUT_FILES:
        expected.append(output / "fastgear" / relpath)

    for path in _expected_chewbbaca_paths(CFG):
        expected.append(path)
    for path in _expected_mlst_paths(CFG):
        expected.append(path)
    for path in _expected_mobsuite_paths(CFG):
        expected.append(path)

    return _dedupe_paths(expected)


def _run_version_cmds(cmds: list[list[str]]) -> str | None:
    for cmd in cmds:
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, check=False)
        except FileNotFoundError:
            continue
        out = (res.stdout or "") + (res.stderr or "")
        for line in out.splitlines():
            s = line.strip()
            if s:
                return s
    return None


def _build_marker_resistance_metadata_inputs(marker_paths: list[Path]) -> dict[str, dict]:
    try:
        marker_index = load_db_catalog_index(MARKER_DB_CATALOG_PATH, logger=logger)
    except Exception as exc:
        logger.warning("Failed to load marker DB catalog for resistance metadata: %s", exc)
        return {}

    return marker_index.resistance_metadata_map(
        marker_paths,
        logger=logger,
        label="marker database",
    )


def _build_run_inputs(CFG: DataInput, *, thread_plan: dict[str, int]) -> dict[str, object]:
    inputs: dict[str, object] = {
        "genomes": [str(p) for p in CFG.genomes],
        "marker_sequences": [str(p) for p in CFG.marker_sequences],
        "reference_locus_sequences": [str(p) for p in CFG.reference_locus_sequences],
        "snippy_reference": str(CFG.snippy_reference) if CFG.snippy_reference else None,
        "bakta_db": str(CFG.bakta_db) if CFG.bakta_db else None,
        "output": str(CFG.output),
        "cores": CFG.total_cores,
        "window_length": CFG.window_length,
        "diamond": CFG.diamond,
        "bakta": CFG.bakta,
        "dmnd_threads": CFG.dmnd_threads,
        "bakta_threads": CFG.bakta_threads,
        "assembly_qc_threads": thread_plan["assembly_qc"],
        "mobsuite_threads": thread_plan["mobsuite"],
        "chewbbaca_threads": thread_plan["chewbbaca"],
        "mlst_threads": thread_plan["mlst"],
        "panaroo_threads": thread_plan["panaroo"],
        "gubbins_threads": thread_plan["gubbins"],
        "phylogeny_snippy_threads": thread_plan["phylogeny_snippy"],
        "phylogeny_raxml_threads": thread_plan["phylogeny_raxml"],
        "marker_categories": CFG.marker_categories,
        "reference_locus_categories": CFG.reference_locus_categories,
        "chewbbaca": configured_chewbbaca_settings(CFG.chewbbaca_config),
        "mobsuite": configured_mobsuite_settings(CFG.mobsuite_config),
        "assembly_qc": configured_assembly_qc_settings(CFG.assembly_qc_config),
        "fastgear": CFG.fastgear_config or configured_fastgear_settings(),
    }
    marker_resistance_metadata = _build_marker_resistance_metadata_inputs(CFG.marker_sequences)
    if marker_resistance_metadata:
        inputs["marker_resistance_metadata"] = marker_resistance_metadata
    return inputs


def _build_run_config(
    CFG: DataInput,
    *,
    sample_set_id: str,
    run_id: str,
    annotation_id_strategy: str,
    thread_plan: dict[str, int],
) -> dict:
    inputs = _build_run_inputs(CFG, thread_plan=thread_plan)
    return {
        "created_at": _now_iso(),
        "sample_set_id": sample_set_id,
        "run_id": run_id,
        "annotation_id_strategy": annotation_id_strategy,
        "inputs": inputs,
        "thread_plan": thread_plan,
        "output_dir": str(CFG.output),
    }


def _build_run_report(
    CFG: DataInput,
    *,
    sample_set_id: str,
    run_id: str,
    thread_plan: dict[str, int],
    start_time: str,
    steps: dict,
    tools_used: set[str],
    status: str,
    error: str | None,
    outputs_validated: bool,
    validation_errors: list[str],
) -> dict:
    report: dict = {
        "start_time": start_time,
        "end_time": _now_iso(),
        "status": status,
        "outputs_validated": bool(outputs_validated),
        "sample_set_id": sample_set_id,
        "run_id": run_id,
        "inputs": _build_run_inputs(CFG, thread_plan=thread_plan),
        "steps": steps,
        "tools": {
            "used": sorted(tools_used),
            "versions": _collect_versions(tools_used),
        },
        "output_dir": str(CFG.output),
        "output_files": _collect_output_files(CFG.output),
    }
    if error:
        report["error"] = error
    if validation_errors:
        report["validation_errors"] = list(validation_errors)
    return report


def _build_run_state(
    CFG: DataInput,
    *,
    sample_set_id: str,
    run_id: str,
    thread_plan: dict[str, int],
    start_time: str,
    steps: dict,
    status: str,
    current_step: str | None,
    error: str | None,
    outputs_validated: bool,
    validation_errors: list[str],
    report_written: bool,
) -> dict:
    payload: dict[str, object] = {
        "start_time": start_time,
        "updated_at": _now_iso(),
        "status": status,
        "outputs_validated": bool(outputs_validated),
        "sample_set_id": sample_set_id,
        "run_id": run_id,
        "current_step": current_step,
        "inputs": _build_run_inputs(CFG, thread_plan=thread_plan),
        "steps": steps,
        "output_dir": str(CFG.output),
        "stop_requested": _stop_requested(CFG.output),
        "report_written": bool(report_written),
    }
    if status in {"success", "failed", "stopped", "interrupted", "partial"}:
        payload["end_time"] = _now_iso()
    if error:
        payload["error"] = error
    if validation_errors:
        payload["validation_errors"] = list(validation_errors)
    return payload


def _write_run_lock(output_dir: Path, *, run_id: str, start_time: str) -> None:
    data = {
        "run_id": run_id,
        "pid": os.getpid(),
        "host": socket.gethostname(),
        "start_time": start_time,
    }
    _write_json(output_dir / RUN_LOCK_FILENAME, data)


def _clear_run_lock(output_dir: Path) -> None:
    try:
        (output_dir / RUN_LOCK_FILENAME).unlink()
    except FileNotFoundError:
        return
    except Exception:
        return


def _stop_requested(output_dir: Path) -> bool:
    try:
        return (output_dir / RUN_STOP_FILENAME).exists()
    except Exception:
        return False


def _collect_versions(tools_used: set[str]) -> dict[str, str]:
    def _cmds(
        tool: str,
        *args_variants: list[str],
        executables: tuple[str, ...] | None = None,
    ) -> list[list[str]]:
        out: list[list[str]] = []
        executables = executables or (tool,)
        for executable in executables:
            env_cmd = conda_tool_cmd(tool, executable=executable)
            for args in args_variants:
                if env_cmd:
                    out.append([*env_cmd, *args])
                out.append([executable, *args])
        return out

    version_cmds = {
        "bakta": _cmds("bakta", ["--version"], ["-V"]),
        "chewbbaca": _cmds(
            "chewbbaca",
            ["--version"],
            ["--help"],
            executables=("chewBBACA.py", "chewbbaca", "chewBBACA"),
        ),
        "mlst": _cmds("mlst", ["--version"], ["--help"]),
        "quast": _cmds("quast", ["--version"], ["--help"], executables=("quast.py", "quast")),
        "checkm2": _cmds("checkm2", ["--version"], ["--help"]),
        "mob_recon": _cmds("mob_recon", ["--version"], ["--help"]),
        "mob_typer": _cmds("mob_typer", ["--version"], ["--help"]),
        "diamond": _cmds("diamond", ["--version"], ["version"], ["-v"]),
        "panaroo": _cmds("panaroo", ["--version"], ["-v"], ["--help"]),
        "run_gubbins.py": _cmds(
            "run_gubbins.py",
            ["--version"],
            ["--help"],
            executables=("run_gubbins.py", "gubbins"),
        ),
        "snippy": _cmds("snippy", ["--version"], ["--help"]),
        "raxml-ng": _cmds("raxml-ng", ["--version"], ["--help"]),
        "docker": _cmds("docker", ["--version"]),
    }
    versions: dict[str, str] = {}
    for tool in sorted(tools_used):
        cmds = version_cmds.get(tool)
        if not cmds:
            versions[tool] = "unknown"
            continue
        v = _run_version_cmds(cmds)
        versions[tool] = v if v is not None else "not found"
    return versions


def _write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as f:
            json.dump(data, f, indent=2, default=str)
            f.flush()
            os.fsync(f.fileno())
            tmp_path = Path(f.name)
        os.replace(tmp_path, path)
    finally:
        if tmp_path is not None and tmp_path.exists():
            try:
                tmp_path.unlink()
            except Exception:
                pass


def run_bakta(CFG: DataInput, *, bakta_threads: int) -> bool:
    try:
        logger.info("Running Bakta on genomes.")
        db = (Path(CFG.bakta_db).expanduser().resolve() if CFG.bakta_db else None)
        if db and not db.exists():
            raise FileNotFoundError(f"Bakta DB not found: {db}")
        run_bakta_on_fastas(
            fasta_paths=CFG.genomes,
            output_path=CFG.output,
            logger=logger,
            bakta_db=str(db) if db else "",
            bakta_threads=bakta_threads,
        )
        logger.info("Bakta completed.")
        return True
    except Exception as e:
        logger.error(f"Bakta workflow failed: {e}")
        raise


def run_assembly_quality(
    CFG: DataInput,
    *,
    assembly_qc_threads: int,
):
    try:
        logger.info("Running assembly quality verification with QUAST and CheckM2.")
        return run_assembly_quality_tool(
            genome_files=CFG.genomes,
            output_path=CFG.output,
            logger=logger,
            threads=assembly_qc_threads,
            config=CFG.assembly_qc_config,
        )
    except Exception as exc:
        logger.error("Assembly QC workflow failed: %s", exc)
        raise


def _collect_bakta_gffs_for_genomes(CFG: DataInput) -> list[Path]:
    gff_candidates = [p for p in (CFG.output / "bakta").rglob("*.gff*") if p.is_file()]
    expected_stems = [p.stem for p in CFG.genomes]
    expected_set = set(expected_stems)

    gff_map: dict[str, Path] = {}
    for p in gff_candidates:
        stem = p.stem
        if stem not in expected_set:
            continue
        existing = gff_map.get(stem)
        if existing is None:
            gff_map[stem] = p
        elif existing.suffix.lower() != ".gff3" and p.suffix.lower() == ".gff3":
            gff_map[stem] = p

    missing = [s for s in expected_stems if s not in gff_map]
    if missing:
        logger.error(
            "Missing Bakta GFFs for %d genome(s): %s",
            len(missing),
            ", ".join(missing),
        )
        return []

    return [gff_map[s] for s in expected_stems if s in gff_map]


def run_panaroo(
    CFG: DataInput,
    *,
    ran_bakta: bool,
    panaroo_threads: int,
) -> Path | None:
    if not ran_bakta:
        logger.info("Skipping Panaroo: Bakta did not run, so no GFFs available.")
        return None
    try:
        gff_files = _collect_bakta_gffs_for_genomes(CFG)
        if not gff_files:
            logger.warning("Skip Panaroo: no GFF files found under output/bakta.")
            return None

        logger.info(f"Running Panaroo on {len(gff_files)} Bakta GFFs.")
        return run_panaroo_tool(
            gff_files=gff_files,
            output_path=CFG.output,
            logger=logger,
            threads=panaroo_threads,
        )
    except Exception as e:
        logger.error(f"Panaroo workflow failed: {e}")
        raise


def run_gubbins(
    CFG: DataInput,
    *,
    gubbins_threads: int,
):
    try:
        logger.info("Running Gubbins on Panaroo core alignment.")
        return run_gubbins_tool(
            output_path=CFG.output,
            logger=logger,
            threads=gubbins_threads,
        )
    except Exception as exc:
        logger.error("Gubbins workflow failed: %s", exc)
        raise


def run_phylogeny(
    CFG: DataInput,
    *,
    snippy_threads: int,
    raxml_threads: int,
):
    try:
        logger.info("Running phylogeny tree inference from Panaroo and Snippy outputs.")
        return run_phylogeny_tool(
            genome_files=CFG.genomes,
            output_path=CFG.output,
            logger=logger,
            snippy_threads=snippy_threads,
            raxml_threads=raxml_threads,
            reference_genome=CFG.snippy_reference,
        )
    except Exception as exc:
        logger.error("Phylogeny workflow failed: %s", exc)
        raise


def run_chewbbaca(
    CFG: DataInput,
    *,
    chewbbaca_threads: int,
) -> Path | None:
    try:
        logger.info("Running chewBBACA AlleleCall on genomes.")
        return run_chewbbaca_tool(
            genome_files=CFG.genomes,
            output_path=CFG.output,
            logger=logger,
            threads=chewbbaca_threads,
            config=CFG.chewbbaca_config,
        )
    except Exception as exc:
        logger.error("chewBBACA workflow failed: %s", exc)
        raise


def run_mlst(
    CFG: DataInput,
    *,
    mlst_threads: int,
) -> Path:
    try:
        logger.info("Running mlst sequence typing on genomes.")
        return run_mlst_tool(
            genome_files=CFG.genomes,
            output_path=CFG.output,
            logger=logger,
            threads=mlst_threads,
        )
    except Exception as exc:
        logger.error("mlst workflow failed: %s", exc)
        raise


def run_mobsuite(
    CFG: DataInput,
    *,
    mobsuite_threads: int,
):
    try:
        logger.info("Running MOB-suite plasmid reconstruction and typing on genomes.")
        return run_mobsuite_tool(
            genome_files=CFG.genomes,
            output_path=CFG.output,
            logger=logger,
            threads=mobsuite_threads,
            config=CFG.mobsuite_config,
        )
    except Exception as exc:
        logger.error("MOB-suite workflow failed: %s", exc)
        raise


def run_diamond(CFG: DataInput, *, diamond_makedb_threads: int) -> list[Path]:
    try:
        logger.info("Running DIAMOND on MarkerSequence FASTAs.")
        dmnd_build_folder = CFG.output / "dmnd" / "queries"
        dmnd_build_folder.mkdir(parents=True, exist_ok=True)

        query_db_files = build_diamond_db_for_files(
            CFG.marker_sequences, dmnd_build_folder, logger, diamond_makedb_threads
        )

        annotations_folder = CFG.output / "annotations" / "queries"
        annotations_folder.mkdir(parents=True, exist_ok=True)
        return query_db_files
    except Exception as e:
        logger.error(f"DIAMOND (queries) workflow failed: {e}")
        raise


def run_reference_locus_search(
    CFG: DataInput,
    *,
    diamond_makedb_threads: int,
    diamond_blastx_threads: int,
) -> list[Path]:
    try:
        logger.info("Running DIAMOND for reference locus sequences against genomes.")
        db_folder = CFG.output / "dmnd" / "reference_loci"
        db_folder.mkdir(parents=True, exist_ok=True)
        locus_db_files = build_diamond_db_for_files(
            CFG.reference_locus_sequences, db_folder, logger, diamond_makedb_threads
        )
        if not locus_db_files:
            logger.warning("No reference locus DBs available; skipping reference locus search.")
            return []

        out_dir = CFG.output / "annotations" / "reference_loci"
        out_dir.mkdir(parents=True, exist_ok=True)
        run_diamond_annotations(
            genome_files=CFG.genomes,
            db_files=locus_db_files,
            output_path=out_dir,
            logger=logger,
            threads=diamond_blastx_threads,
            sanitize=True,
        )

        tsvs = sorted(out_dir.rglob("*.tsv"))
        hits_df = load_diamond_tsvs(tsvs, CFG.reference_locus_categories, logger)
        logger.info(
            "Reference locus hits loaded: %d rows from %d TSV(s).",
            len(hits_df),
            len(tsvs),
        )
        return tsvs
    except Exception as e:
        logger.error(f"Reference locus DIAMOND workflow failed: {e}")
        raise


def _normalize_contig_for_mobsuite(value) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    if text.lower().startswith("contig_"):
        text = text.split("_", 1)[1]
    try:
        fval = float(text)
    except Exception:
        return text
    if fval.is_integer():
        return str(int(fval))
    return text


def _attach_mobsuite_predictions(df: pd.DataFrame, output_path: Path) -> pd.DataFrame:
    predictions = load_mobsuite_contig_predictions(output_path)
    if not predictions:
        for col in (
            "MobSuitePrediction",
            "MobSuiteIsPlasmid",
            "MobSuiteSourceContig",
            "MobSuitePlasmidID",
            "MobSuiteMoleculeType",
            "MobSuitePrimaryClusterID",
            "MobSuiteSecondaryClusterID",
            "MobSuiteLength",
            "MobSuiteRepliconTypes",
            "MobSuiteRelaxaseTypes",
            "MobSuiteMPFType",
            "MobSuiteOriTTypes",
            "MobSuitePredictedMobility",
            "MobSuiteMashNearestNeighbor",
            "MobSuiteMashNeighborDistance",
        ):
            if col not in df.columns:
                df[col] = pd.NA
        if "MobSuiteIsPlasmid" in df.columns:
            df["MobSuiteIsPlasmid"] = df["MobSuiteIsPlasmid"].fillna("false")
        return df

    def _lookup(row: pd.Series) -> dict[str, str]:
        sample = str(row.get("Sample") or "").strip()
        contig = _normalize_contig_for_mobsuite(row.get("Contig"))
        return predictions.get((sample, contig), {})

    hits = df.apply(_lookup, axis=1)
    df["MobSuitePrediction"] = hits.apply(lambda item: item.get("Prediction") or pd.NA)
    df["MobSuiteIsPlasmid"] = hits.apply(lambda item: item.get("IsPlasmid") or "false")
    df["MobSuiteSourceContig"] = hits.apply(lambda item: item.get("SourceContig") or pd.NA)
    df["MobSuitePlasmidID"] = hits.apply(lambda item: item.get("PlasmidID") or pd.NA)
    df["MobSuiteMoleculeType"] = hits.apply(lambda item: item.get("MoleculeType") or pd.NA)
    df["MobSuitePrimaryClusterID"] = hits.apply(lambda item: item.get("PrimaryClusterID") or pd.NA)
    df["MobSuiteSecondaryClusterID"] = hits.apply(lambda item: item.get("SecondaryClusterID") or pd.NA)
    df["MobSuiteLength"] = hits.apply(lambda item: item.get("Size") or pd.NA)
    df["MobSuiteRepliconTypes"] = hits.apply(lambda item: item.get("RepliconTypes") or pd.NA)
    df["MobSuiteRelaxaseTypes"] = hits.apply(lambda item: item.get("RelaxaseTypes") or pd.NA)
    df["MobSuiteMPFType"] = hits.apply(lambda item: item.get("MPFType") or pd.NA)
    df["MobSuiteOriTTypes"] = hits.apply(lambda item: item.get("OriTTypes") or pd.NA)
    df["MobSuitePredictedMobility"] = hits.apply(lambda item: item.get("PredictedMobility") or pd.NA)
    df["MobSuiteMashNearestNeighbor"] = hits.apply(lambda item: item.get("MashNearestNeighbor") or pd.NA)
    df["MobSuiteMashNeighborDistance"] = hits.apply(lambda item: item.get("MashNeighborDistance") or pd.NA)
    return df


def collect_annotations(
    CFG: DataInput,
    *,
    ran_bakta: bool,
    query_db_files: list[Path],
    diamond_blastx_threads: int,
    reuse_marker_hits: bool = False,
    annotation_id_strategy: str = ANNOTATION_ID_STRATEGY,
) -> pd.DataFrame:
    marker_resistance_metadata = _build_marker_resistance_metadata_inputs(CFG.marker_sequences)

    if ran_bakta:
        logger.info("Getting Bakta annotations")
        bakta_json_files = sorted((CFG.output / "bakta").rglob("*.json"), key=lambda p: str(p))
        if bakta_json_files:
            logger.info(f"Found {len(bakta_json_files)} Bakta JSON files for parsing.")
            bakta_df = parse_multiple_bakta_json(bakta_json_files, logger)
        else:
            logger.warning("No Bakta JSON files found for parsing.")
            bakta_df = pd.DataFrame()
    else:
        bakta_df = pd.DataFrame()

    bakta_df = _ensure_annotation_columns(bakta_df)

    if not bakta_df.empty:
        bakta_df["Method"] = "bakta"
        bakta_df["AnnotationCategory"] = bakta_df.apply(compute_category, axis=1)
        bakta_df = assign_annotation_ids(bakta_df, prefix="BKT")
        bakta_df = _ensure_annotation_columns(bakta_df)
        bakta_df = bakta_df[ANNOTATION_COLUMNS]

        with tempfile.TemporaryDirectory(dir=CFG.output) as tdir:
            tdir = Path(tdir)
            fa_path = write_annotation_fasta(
                bakta_df,
                tdir,
                logger,
                id_col="AnnotationID",
                seq_col="Nucleotide_Sequence",
                filename="bakta_annotations_nt.fasta",
            )

            if fa_path is not None and query_db_files:
                annotations_folder = CFG.output / "annotations" / "bakta_vs_queries"
                annotations_folder.mkdir(parents=True, exist_ok=True)
                if reuse_marker_hits and _all_files_exist(_expected_marker_hit_paths(CFG)):
                    logger.info("Reusing existing DIAMOND marker hit TSVs under %s.", annotations_folder)
                else:
                    run_diamond_annotations(
                        genome_files=[fa_path],
                        db_files=query_db_files,
                        output_path=annotations_folder,
                        logger=logger,
                        threads=diamond_blastx_threads,
                        sanitize=False,
                    )
                    logger.info("DIAMOND on Bakta annotations vs MarkerSequence DBs completed.")

    else:
        logger.warning("No Bakta annotations parsed; nothing to export to FASTA.")
        bakta_df = bakta_df[ANNOTATION_COLUMNS]

    diamond_annotations = read_annotations(CFG.output / "annotations" / "bakta_vs_queries")
    if not diamond_annotations.empty:
        logger.info(f"Total DIAMOND annotations parsed (raw): {len(diamond_annotations)}")
        diamond_annotations = collapse_diamond_annotations(diamond_annotations)
        diamond_annotations = enrich_marker_hits_with_resistance(
            diamond_annotations,
            marker_resistance_metadata,
            logger=logger,
        )
        logger.info(
            f"DIAMOND annotations after collapsing per (AnnotationID, Method): "
            f"{len(diamond_annotations)}"
        )

        marker_cols = [c for c in diamond_annotations.columns if c != "AnnotationID"]
        marker_map = (
            diamond_annotations
            .groupby("AnnotationID", sort=False)[marker_cols]
            .apply(lambda g: g.to_dict(orient="records"))
            .to_dict()
        )
    else:
        logger.warning("No DIAMOND annotations parsed.")
        marker_map = {}

    if not bakta_df.empty:
        if "AnnotationID" not in bakta_df.columns:
            logger.warning("Bakta DataFrame has no AnnotationID; cannot attach markers.")
        else:
            bakta_df["Markers"] = bakta_df["AnnotationID"].map(marker_map)
            bakta_df["Markers"] = bakta_df["Markers"].apply(
                lambda v: v if isinstance(v, list) else []
            )

            consensus = bakta_df["Markers"].apply(infer_annotation_category_and_source)
            bakta_df["ConsensusCategory"] = consensus.str[0]
            bakta_df["ConsensusSource"] = consensus.str[1]
            bakta_df["ConsensusDonorProduct"] = consensus.str[2]

            resistance_maps = build_annotation_resistance_maps(diamond_annotations)
            for col, col_map in resistance_maps.items():
                bakta_df[col] = bakta_df["AnnotationID"].map(col_map).replace({"": pd.NA})

            logger.info(
                "Attached marker lists, consensus category/source/donor product, and resistance metadata "
                "to Bakta annotations."
            )
    else:
        logger.warning("Bakta annotations are empty; no markers attached.")

    if not bakta_df.empty:
        bakta_df = assign_gene_ids(bakta_df)
        bakta_df = _attach_mobsuite_predictions(bakta_df, CFG.output)

    # No CSVs are written; annotations are kept in-memory only.
    return bakta_df


def run_workflow(CFG: DataInput):
    # Establish the container-execution context for this run: every tool_cmd()
    # call below wraps `conda run -n <env> <exe>` into `docker run`/`singularity
    # exec` with the run's dirs bind-mounted at identical paths (or stays native
    # under BACTOLOMEW_NATIVE / runtime:native for local dev).
    _ctx = make_container_context(
        output_dir=CFG.output,
        extra_paths=[
            *CFG.genomes,
            *CFG.marker_sequences,
            *CFG.reference_locus_sequences,
            CFG.snippy_reference,
            CFG.bakta_db,
        ],
    )
    _token = set_container_context(_ctx)
    if _ctx.runtime != "native":
        logger.info(f"Tool runtime: {_ctx.runtime} (image {_ctx.image or _ctx.sif_path})")
    try:
        return _run_workflow_impl(CFG)
    finally:
        reset_container_context(_token)


def _run_workflow_impl(CFG: DataInput):
    validate_chewbbaca_settings(CFG.chewbbaca_config)
    try:
        CFG.output.mkdir(parents=True, exist_ok=True)
        logger.info(f"Output folder is ready: {CFG.output}")
    except Exception as e:
        msg = f"Failed to create output folder '{CFG.output}': {e}"
        logger.error(msg)
        raise RuntimeError(msg) from e

    thread_plan = _compute_thread_plan(
        CFG.total_cores,
        n_genomes=len(CFG.genomes),
        n_marker_dbs=len(CFG.marker_sequences),
        n_locus_dbs=len(CFG.reference_locus_sequences),
    )

    existing_run_cfg = _read_json_dict(CFG.output / RUN_CONFIG_FILENAME) or {}
    sample_set_id = str(existing_run_cfg.get("sample_set_id") or "").strip() or _compute_sample_set_id(CFG.genomes)
    configured_run_id = str(CFG.run_id).strip() if getattr(CFG, "run_id", None) else ""
    existing_run_id = str(existing_run_cfg.get("run_id") or "").strip()
    annotation_id_strategy = str(existing_run_cfg.get("annotation_id_strategy") or "").strip() or ANNOTATION_ID_STRATEGY
    run_id = configured_run_id or existing_run_id or _new_attempt_id(sample_set_id)
    start_time = _now_iso()

    run_config = _build_run_config(
        CFG,
        sample_set_id=sample_set_id,
        run_id=run_id,
        annotation_id_strategy=annotation_id_strategy,
        thread_plan=thread_plan,
    )
    _write_json(CFG.output / RUN_CONFIG_FILENAME, run_config)

    try:
        (CFG.output / RUN_STOP_FILENAME).unlink()
    except FileNotFoundError:
        pass
    except Exception:
        pass

    expected_outputs = _build_expected_outputs(CFG, run_id=run_id)
    _write_json(
        CFG.output / "expected_outputs.json",
        {
            "created_at": _now_iso(),
            "sample_set_id": sample_set_id,
            "run_id": run_id,
            "output_dir": str(CFG.output),
            "files": expected_outputs,
        },
    )
    _write_run_lock(CFG.output, run_id=run_id, start_time=start_time)
    steps: dict = {}
    tools_used: set[str] = set()
    error: Exception | None = None
    error_text: str | None = None
    current_step: str | None = None
    outputs_validated = False
    validation_errors: list[str] = []

    def _flush_run_state(
        status: str,
        *,
        error_override: str | None = None,
        report_written: bool = False,
    ) -> None:
        try:
            _write_json(
                CFG.output / RUN_STATE_FILENAME,
                _build_run_state(
                    CFG,
                    sample_set_id=sample_set_id,
                    run_id=run_id,
                    thread_plan=thread_plan,
                    start_time=start_time,
                    steps=steps,
                    status=status,
                    current_step=current_step,
                    error=error_override,
                    outputs_validated=outputs_validated,
                    validation_errors=validation_errors,
                    report_written=report_written,
                ),
            )
        except Exception as exc:
            logger.warning("Failed to persist %s: %s", RUN_STATE_FILENAME, exc)

    _flush_run_state("running")

    def _check_stop() -> None:
        if _stop_requested(CFG.output):
            _flush_run_state("stopped", error_override="Run stopped by user")
            raise RunStoppedError("Run stopped by user")

    def _step_start(name: str) -> None:
        nonlocal current_step
        current_step = name
        steps[name] = {
            "status": "running",
            "start_time": _now_iso(),
            "outputs_validated": False,
            "validation_errors": [],
        }
        _flush_run_state("running")

    def _step_end(name: str, status: str, **info) -> None:
        nonlocal current_step
        step = steps.get(name, {})
        step.update({"status": status, "end_time": _now_iso()})
        step.update(info)
        if "outputs_validated" not in step:
            step["outputs_validated"] = status == "success"
        validation_errors_step = step.get("validation_errors")
        if not isinstance(validation_errors_step, list):
            step["validation_errors"] = []
        steps[name] = step
        current_step = None
        live_status = "running"
        if status == "failed":
            live_status = "failed"
        elif status == "partial":
            live_status = "partial"
        error_override = str(step.get("error") or "").strip() or None
        _flush_run_state(live_status, error_override=error_override)

    def _finish_step_success(name: str, *, validator_errors: list[str] | None = None, **info) -> None:
        errors = list(validator_errors or [])
        if errors:
            _step_end(
                name,
                "failed",
                outputs_validated=False,
                validation_errors=errors,
                **info,
            )
            raise RuntimeError(f"{name} outputs failed validation: {errors[0]}")
        _step_end(
            name,
            "success",
            outputs_validated=True,
            validation_errors=[],
            **info,
        )

    all_marker_sequences = CFG.marker_sequences + CFG.reference_locus_sequences
    for f in CFG.genomes:
        logger.info(f"Genome: {f.name}")
    for f in all_marker_sequences:
        logger.info(f"MarkerSequence: {f.name}")

    try:
        # Assembly quality verification
        _check_stop()
        _step_start("assembly_qc")
        assembly_qc_paths_existing = _expected_assembly_qc_paths(CFG)
        if _validate_assembly_qc_reuse_outputs(CFG):
            _finish_step_success(
                "assembly_qc",
                reason="reused_existing_output",
                output=str(CFG.output / "assembly_qc"),
                threads=thread_plan["assembly_qc"],
            )
        else:
            if _all_files_exist(assembly_qc_paths_existing):
                logger.warning("Existing assembly QC outputs are incomplete or invalid; rerunning QUAST and CheckM2.")
            try:
                tools_used.add("quast")
                tools_used.add("checkm2")
                assembly_qc_result = run_assembly_quality(
                    CFG,
                    assembly_qc_threads=thread_plan["assembly_qc"],
                )
            except Exception as e:
                _step_end("assembly_qc", "failed", error=str(e))
                raise
            else:
                _finish_step_success(
                    "assembly_qc",
                    validator_errors=validate_assembly_qc_outputs(CFG.output, CFG.genomes),
                    output=str(assembly_qc_result.output_dir),
                    summary=str(assembly_qc_result.summary_tsv),
                    sample_count=int(assembly_qc_result.sample_count),
                    threads=thread_plan["assembly_qc"],
                )

        # Bakta
        _check_stop()
        _step_start("bakta")
        bakta_reuse_paths = _expected_bakta_reuse_paths(CFG)
        if _all_files_exist(bakta_reuse_paths) and _validate_bakta_reuse_outputs(CFG):
            ran_bakta = True
            _finish_step_success(
                "bakta",
                reason="reused_existing_output",
                existing=str(CFG.output / "bakta"),
                threads=thread_plan["bakta"],
            )
        else:
            if _all_files_exist(bakta_reuse_paths):
                logger.warning("Existing Bakta outputs are incomplete or invalid; rerunning Bakta.")
            try:
                tools_used.add("bakta")
                ran_bakta = run_bakta(CFG, bakta_threads=thread_plan["bakta"])
            except Exception as e:
                _step_end("bakta", "failed", error=str(e))
                raise
            else:
                bakta_json_files = list((CFG.output / "bakta").rglob("*.json"))
                validation = validate_bakta_outputs(CFG.output, CFG.genomes)
                if not bakta_json_files and not validation:
                    validation = ["Bakta completed without producing any JSON outputs."]
                _finish_step_success(
                    "bakta",
                    validator_errors=validation,
                    json_count=len(bakta_json_files),
                    threads=thread_plan["bakta"],
                )

        # Existing outputs to decide whether to re-run downstream tools.
        panaroo_out_existing = ensure_panaroo_out(CFG.output, logger)
        panaroo_gene_presence_existing = (
            (panaroo_out_existing / "gene_presence_absence.csv")
            if panaroo_out_existing and (panaroo_out_existing / "gene_presence_absence.csv").is_file()
            else None
        )
        fastgear_dir_existing = CFG.output / "fastgear"
        fastgear_sample_clades_existing = fastgear_dir_existing / "sample_clades.tsv"
        fastgear_lineage_existing = fastgear_dir_existing / "output" / "lineage_information.txt"
        fastgear_mat_existing = fastgear_dir_existing / "fastgear_res.mat"
        gubbins_paths_existing = _expected_gubbins_paths(CFG)
        phylogeny_paths_existing = _expected_phylogeny_paths(CFG)
        chewbbaca_out_existing = ensure_chewbbaca_out(CFG.output, logger)
        mlst_paths_existing = _expected_mlst_paths(CFG)
        mobsuite_paths_existing = _expected_mobsuite_paths(CFG)
        query_db_files_existing = _expected_marker_db_paths(CFG)
        marker_hit_paths_existing = _expected_marker_hit_paths(CFG)
        ref_hit_paths_existing = _expected_reference_hit_paths(CFG)

        # MLST
        _check_stop()
        _step_start("mlst")
        if _validate_mlst_reuse_outputs(CFG):
            _finish_step_success(
                "mlst",
                reason="reused_existing_output",
                output=str(CFG.output / "mlst"),
                threads=thread_plan["mlst"],
            )
        else:
            if _all_files_exist(mlst_paths_existing):
                logger.warning("Existing MLST outputs are incomplete or invalid; rerunning mlst.")
            try:
                tools_used.add("mlst")
                mlst_out = run_mlst(CFG, mlst_threads=thread_plan["mlst"])
            except Exception as e:
                _step_end("mlst", "failed", error=str(e))
                raise
            else:
                _finish_step_success(
                    "mlst",
                    validator_errors=validate_mlst_outputs(CFG.output, CFG.genomes),
                    output=str(mlst_out),
                    threads=thread_plan["mlst"],
                )

        # MOB-suite plasmid reconstruction and contig assignment
        _check_stop()
        _step_start("mobsuite")
        if _validate_mobsuite_reuse_outputs(CFG):
            _finish_step_success(
                "mobsuite",
                reason="reused_existing_output",
                output=str(CFG.output / "mobsuite"),
                threads=thread_plan["mobsuite"],
            )
        else:
            if _all_files_exist(mobsuite_paths_existing):
                logger.warning("Existing MOB-suite outputs are incomplete or invalid; rerunning MOB-suite.")
            try:
                tools_used.add("mob_recon")
                tools_used.add("mob_typer")
                mobsuite_result = run_mobsuite(CFG, mobsuite_threads=thread_plan["mobsuite"])
            except Exception as e:
                _step_end("mobsuite", "failed", error=str(e))
                raise
            else:
                _finish_step_success(
                    "mobsuite",
                    validator_errors=validate_mobsuite_outputs(CFG.output, CFG.genomes),
                    output=str(mobsuite_result.output_dir),
                    contigs=str(mobsuite_result.contigs_tsv),
                    plasmids=str(mobsuite_result.plasmids_tsv),
                    sample_count=int(mobsuite_result.sample_count),
                    plasmid_contig_count=int(mobsuite_result.plasmid_contig_count),
                    plasmid_count=int(mobsuite_result.plasmid_count),
                    threads=thread_plan["mobsuite"],
                )

        # Panaroo
        _check_stop()
        _step_start("panaroo")
        panaroo_reused_existing = False
        if not ran_bakta:
            _step_end("panaroo", "skipped", reason="bakta_not_run", outputs_validated=False, validation_errors=[])
            panaroo_out = None
        elif panaroo_gene_presence_existing:
            panaroo_out = panaroo_out_existing
            panaroo_reused_existing = True
            _finish_step_success(
                "panaroo",
                validator_errors=validate_panaroo_outputs(CFG.output),
                reason="reused_existing_output",
                existing=str(panaroo_gene_presence_existing),
                output=str(panaroo_out) if panaroo_out else None,
                threads=thread_plan["panaroo"],
            )
        else:
            try:
                tools_used.add("panaroo")
                panaroo_out = run_panaroo(
                    CFG,
                    ran_bakta=ran_bakta,
                    panaroo_threads=thread_plan["panaroo"],
                )
            except Exception as e:
                _step_end("panaroo", "failed", error=str(e))
                raise
            else:
                if not panaroo_out:
                    _step_end(
                        "panaroo",
                        "failed",
                        reason="panaroo_failed_or_no_gff",
                        threads=thread_plan["panaroo"],
                        outputs_validated=False,
                        validation_errors=["Panaroo produced no output directory."],
                    )
                    raise RuntimeError("Panaroo produced no output directory.")
                _finish_step_success(
                    "panaroo",
                    validator_errors=validate_panaroo_outputs(CFG.output),
                    output=str(panaroo_out),
                    threads=thread_plan["panaroo"],
                )

        # Gubbins recombination filtering on the Panaroo core alignment.
        _check_stop()
        _step_start("gubbins")
        gubbins_regenerated = False
        if not panaroo_out:
            _step_end(
                "gubbins",
                "failed",
                reason="no_panaroo_output",
                outputs_validated=False,
                validation_errors=["Panaroo output required before Gubbins can run."],
            )
            raise RuntimeError("Panaroo output required before Gubbins can run.")
        elif _validate_gubbins_reuse_outputs(CFG):
            _finish_step_success(
                "gubbins",
                reason="reused_existing_output",
                output=str(CFG.output / "gubbins"),
                final_tree=str(CFG.output / "gubbins" / "core_gene.final_tree.tre"),
                filtered_alignment=str(CFG.output / "gubbins" / "core_gene.filtered_polymorphic_sites.fasta"),
                recombination_gff=str(CFG.output / "gubbins" / "core_gene.recombination_predictions.gff"),
                threads=thread_plan["gubbins"],
            )
        else:
            if _all_files_exist(gubbins_paths_existing):
                logger.warning("Existing Gubbins outputs are incomplete or invalid; rerunning Gubbins.")
            try:
                tools_used.add("run_gubbins.py")
                gubbins_result = run_gubbins(CFG, gubbins_threads=thread_plan["gubbins"])
            except Exception as e:
                _step_end("gubbins", "failed", error=str(e))
                raise
            else:
                gubbins_regenerated = True
                _finish_step_success(
                    "gubbins",
                    validator_errors=validate_gubbins_outputs(CFG.output),
                    output=str(gubbins_result.output_dir),
                    input_alignment=str(gubbins_result.input_alignment),
                    final_tree=str(gubbins_result.final_tree),
                    filtered_alignment=str(gubbins_result.filtered_alignment),
                    recombination_gff=str(gubbins_result.recombination_gff),
                    threads=thread_plan["gubbins"],
                )

        # Phylogeny trees
        _check_stop()
        _step_start("phylogeny")
        if not panaroo_out:
            _step_end(
                "phylogeny",
                "failed",
                reason="no_panaroo_output",
                outputs_validated=False,
                validation_errors=["Panaroo output required before phylogeny tree inference."],
            )
            raise RuntimeError("Panaroo output required before phylogeny tree inference.")
        elif not gubbins_regenerated and not validate_phylogeny_outputs(CFG.output):
            _finish_step_success(
                "phylogeny",
                validator_errors=[],
                reason="reused_existing_output",
                reference=str(CFG.snippy_reference) if CFG.snippy_reference else None,
                core_tree=str(CFG.output / "trees" / "core_gene_tree.nwk"),
                snippy_tree=str(CFG.output / "trees" / "snippy_raxml.tree"),
                snippy_threads=thread_plan["phylogeny_snippy"],
                raxml_threads=thread_plan["phylogeny_raxml"],
            )
        else:
            if gubbins_regenerated:
                logger.info("Gubbins output was regenerated; rebuilding phylogeny from filtered alignment.")
            elif _all_files_exist(phylogeny_paths_existing):
                logger.warning("Existing phylogeny outputs are incomplete or invalid; rebuilding them.")
            try:
                tools_used.add("snippy")
                tools_used.add("raxml-ng")
                phylogeny_result = run_phylogeny(
                    CFG,
                    snippy_threads=thread_plan["phylogeny_snippy"],
                    raxml_threads=thread_plan["phylogeny_raxml"],
                )
            except Exception as e:
                _step_end("phylogeny", "failed", error=str(e))
                raise
            else:
                _finish_step_success(
                    "phylogeny",
                    validator_errors=validate_phylogeny_outputs(CFG.output),
                    reference=str(phylogeny_result.reference_genome),
                    core_tree=str(phylogeny_result.core_tree),
                    snippy_tree=str(phylogeny_result.snippy_tree),
                    snippy_threads=thread_plan["phylogeny_snippy"],
                    raxml_threads=thread_plan["phylogeny_raxml"],
                )

        # fastGEAR clades
        _check_stop()
        _step_start("fastgear_clades")
        reuse_fastgear_output = False
        if not panaroo_out:
            _step_end(
                "fastgear_clades",
                "failed",
                reason="no_panaroo_output",
                outputs_validated=False,
                validation_errors=["Panaroo output required for fastGEAR clade inference."],
            )
            raise RuntimeError(
                "Panaroo output required for fastGEAR clade inference. "
                "Upstream Bakta or Panaroo failed, or no Bakta GFFs were produced."
            )
        elif (
            panaroo_reused_existing
            and fastgear_sample_clades_existing.is_file()
            and fastgear_lineage_existing.is_file()
            and fastgear_mat_existing.is_file()
        ):
            try:
                validate_sample_clades_tsv(
                    fastgear_sample_clades_existing,
                    expected_samples=[p.stem for p in CFG.genomes],
                )
            except Exception as exc:
                logger.warning("Existing fastGEAR sample clades are invalid; re-running fastGEAR: %s", exc)
            else:
                _finish_step_success(
                    "fastgear_clades",
                    validator_errors=validate_fastgear_outputs(CFG.output, CFG.genomes),
                    reason="reused_existing_output",
                    sample_clades=str(fastgear_sample_clades_existing),
                    lineage_information=str(fastgear_lineage_existing),
                    result_mat=str(fastgear_mat_existing),
                )
                reuse_fastgear_output = True

        if not reuse_fastgear_output and not panaroo_reused_existing:
            logger.info("Panaroo output was regenerated; fastGEAR will be rerun against the new core alignment.")

        if not reuse_fastgear_output:
            try:
                tools_used.add("docker")
                fastgear_result = run_fastgear_clades(
                    CFG.output,
                    expected_samples=[p.stem for p in CFG.genomes],
                    logger=logger,
                    fastgear_config=CFG.fastgear_config or configured_fastgear_settings(),
                )
            except Exception as e:
                _step_end("fastgear_clades", "failed", error=str(e))
                raise
            else:
                _finish_step_success(
                    "fastgear_clades",
                    validator_errors=validate_fastgear_outputs(CFG.output, CFG.genomes),
                    sample_clades=str(fastgear_result.sample_clades_tsv),
                    lineage_information=str(fastgear_result.lineage_information),
                    result_mat=str(fastgear_result.result_mat),
                )

        # chewBBACA allele calling
        _check_stop()
        _step_start("chewbbaca")
        chewbbaca_expected_paths = _expected_chewbbaca_paths(CFG)
        if _validate_chewbbaca_reuse_outputs(CFG):
            _finish_step_success(
                "chewbbaca",
                reason="reused_existing_output",
                existing=str((chewbbaca_out_existing or (CFG.output / "chewbbaca" / "out"))),
                output=str((chewbbaca_out_existing or (CFG.output / "chewbbaca" / "out"))),
                threads=thread_plan["chewbbaca"],
            )
        else:
            if _all_files_exist(chewbbaca_expected_paths):
                logger.warning("Existing chewBBACA outputs are incomplete or invalid; rerunning chewBBACA.")
            try:
                tools_used.add("chewbbaca")
                chewbbaca_out = run_chewbbaca(
                    CFG,
                    chewbbaca_threads=thread_plan["chewbbaca"],
                )
            except Exception as e:
                _step_end("chewbbaca", "failed", error=str(e))
                raise
            else:
                _finish_step_success(
                    "chewbbaca",
                    validator_errors=validate_chewbbaca_outputs(CFG.output),
                    output=str(chewbbaca_out) if chewbbaca_out else str(CFG.output / "chewbbaca" / "out"),
                    threads=thread_plan["chewbbaca"],
                )

        # DIAMOND marker DBs
        _check_stop()
        _step_start("diamond_markers")
        if _validate_dmnd_files(query_db_files_existing):
            query_db_files = query_db_files_existing
            _finish_step_success(
                "diamond_markers",
                validator_errors=validate_marker_db_outputs(CFG.output, CFG.marker_sequences),
                reason="reused_existing_output",
                db_count=len(query_db_files),
                threads=thread_plan["diamond_makedb"],
            )
        else:
            if _all_files_exist(query_db_files_existing):
                logger.warning("Existing DIAMOND marker DBs are incomplete; rebuilding them.")
            try:
                tools_used.add("diamond")
                query_db_files = run_diamond(CFG, diamond_makedb_threads=thread_plan["diamond_makedb"])
            except Exception as e:
                _step_end("diamond_markers", "failed", error=str(e))
                raise
            else:
                _finish_step_success(
                    "diamond_markers",
                    validator_errors=validate_marker_db_outputs(CFG.output, CFG.marker_sequences),
                    db_count=len(query_db_files),
                    threads=thread_plan["diamond_makedb"],
                )

        # Reference locus search
        _check_stop()
        _step_start("diamond_reference_loci")
        if not validate_reference_hit_outputs(CFG.output, CFG.genomes, CFG.reference_locus_sequences):
            ref_hits_tsvs = ref_hit_paths_existing
            _finish_step_success(
                "diamond_reference_loci",
                reason="reused_existing_output",
                tsv_count=len(ref_hits_tsvs),
                threads=thread_plan["diamond_blastx"],
            )
        else:
            try:
                tools_used.add("diamond")
                ref_hits_tsvs = run_reference_locus_search(
                    CFG,
                    diamond_makedb_threads=thread_plan["diamond_makedb"],
                    diamond_blastx_threads=thread_plan["diamond_blastx"],
                )
            except Exception as e:
                _step_end("diamond_reference_loci", "failed", error=str(e))
                raise
            else:
                _finish_step_success(
                    "diamond_reference_loci",
                    validator_errors=validate_reference_hit_outputs(
                        CFG.output,
                        CFG.genomes,
                        CFG.reference_locus_sequences,
                    ),
                    tsv_count=len(ref_hits_tsvs),
                    threads=thread_plan["diamond_blastx"],
                )

        # Collect annotations
        _check_stop()
        _step_start("collect_annotations")
        try:
            ann_df = collect_annotations(
                CFG,
                ran_bakta=ran_bakta,
                query_db_files=query_db_files,
                diamond_blastx_threads=thread_plan["diamond_blastx"],
                reuse_marker_hits=not validate_marker_hit_outputs(CFG.output, CFG.marker_sequences),
                annotation_id_strategy=annotation_id_strategy,
            )
        except Exception as e:
            _step_end("collect_annotations", "failed", error=str(e))
            raise
        else:
            status = "success" if not ann_df.empty else "partial"
            if status == "success":
                _finish_step_success(
                    "collect_annotations",
                    validator_errors=validate_marker_hit_outputs(CFG.output, CFG.marker_sequences),
                    rows=int(len(ann_df)),
                    threads=thread_plan["diamond_blastx"],
                )
            else:
                _step_end(
                    "collect_annotations",
                    status,
                    rows=int(len(ann_df)),
                    threads=thread_plan["diamond_blastx"],
                    outputs_validated=False,
                    validation_errors=[],
                )

    except Exception as e:
        error = e
        error_text = str(e)
    finally:
        step_validations: dict[str, dict] = {}
        if isinstance(error, RunStoppedError):
            status = "stopped"
        else:
            status = "failed" if error is not None else "success"
        if error is None:
            step_validations, validation_errors = validate_completed_run(
                CFG.output,
                run_cfg=run_config,
                steps=steps,
            )
            for step_name, validation_info in step_validations.items():
                step = steps.get(step_name)
                if not isinstance(step, dict):
                    continue
                step["outputs_validated"] = bool(validation_info.get("outputs_validated"))
                step["validation_errors"] = list(validation_info.get("validation_errors") or [])
                if step.get("status") == "success" and not step["outputs_validated"]:
                    step["status"] = "partial"
            statuses = [v.get("status") for v in steps.values() if isinstance(v, dict)]
            if any(s == "failed" for s in statuses):
                status = "failed"
                error_text = error_text or "one or more steps failed"
            elif any(s == "partial" for s in statuses):
                status = "partial"
                if validation_errors and not error_text:
                    error_text = "run outputs failed final validation"
            outputs_validated = status == "success" and not validation_errors
        else:
            outputs_validated = False
            validation_errors = []

        run_report = _build_run_report(
            CFG,
            sample_set_id=sample_set_id,
            run_id=run_id,
            thread_plan=thread_plan,
            start_time=start_time,
            steps=steps,
            tools_used=tools_used,
            status=status,
            error=error_text,
            outputs_validated=outputs_validated,
            validation_errors=validation_errors,
        )
        _write_json(CFG.output / RUN_REPORT_FILENAME, run_report)
        _flush_run_state(status, error_override=error_text, report_written=True)
        _clear_run_lock(CFG.output)
        print("\nRUN REPORT\n" + json.dumps(run_report, indent=2, default=str))

    if error is not None:
        raise error


__all__ = ["run_workflow", "logger"]
