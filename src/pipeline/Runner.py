from __future__ import annotations

import json
import hashlib
import subprocess
import tempfile
import os
import socket
import threading
from datetime import datetime
from pathlib import Path

import pandas as pd

from pipeline import __version__
from pipeline.AmrFinderWorkflow import (
    MUTATION_COLUMNS,
    amrfinder_reuse_issue,
    attach_resistance_mutations,
    plan_amrfinder,
    read_amrfinder_summary,
    run_amrfinder_mutations,
)
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
from pipeline.BaktaEvidence import (
    EVIDENCE_DIR,
    EVIDENCE_SUMMARY,
    build_evidence,
    evidence_record,
    evidence_reuse_issue,
    evidence_tsv_path,
    validate_evidence_outputs,
)
from pipeline.BaktaWorkflow import parse_multiple_bakta_json, run_bakta_on_fastas
from pipeline.ChewBBACAWorkflow import (
    CHEWBBACA_OUT_FILES,
    CgmlstSchemaMissing,
    chewbbaca_schema_stamp_issue,
    ensure_chewbbaca_out,
    run_chewbbaca as run_chewbbaca_tool,
    validate_chewbbaca_settings,
)
from pipeline.DataInput import DataInput
from pipeline.DbCatalog import (
    COMPOUND_DB_CATALOG_PATH,
    MARKER_DB_CATALOG_PATH,
    load_db_catalog_index,
    run_database_record,
)
from pipeline.DiamondWorkflow import (
    build_diamond_db_for_files,
    collapse_diamond_annotations,
    dmnd_is_current,
    hit_tsv_is_current,
    load_diamond_tsvs,
    read_annotations,
    run_diamond_annotations,
)
from pipeline.CoreAlignment import (
    CORE_ALIGNMENT_DIR,
    CORE_ALIGNMENT_OUT_FILES,
    core_alignment_reuse_issue,
    ordered_alignment_sha256,
    read_ordering_summary,
    run_core_alignment,
)
from pipeline.CoreTreeWorkflow import (
    IQTREE_DIR,
    core_tree_reuse_issue,
    read_core_tree_stamp,
    refresh_gubbins_notes,
    run_core_tree,
    run_recombination_refit,
)
from pipeline.FastBapsWorkflow import fastbaps_reuse_issue, run_fastbaps_clades
from pipeline.FastGearWorkflow import fastgear_engine, fastgear_reuse_issue, run_fastgear_clades
from pipeline.PopulationStructure import (
    PopulationChoice,
    configured_population_choice,
    method_label,
    population_method_explicit,
)
from pipeline.GubbinsWorkflow import (
    GroupPlan,
    expected_gubbins_files,
    gubbins_reuse_issue,
    planned_gubbins_files,
    read_fastgear_lineages,
    read_gubbins_summary,
    resolve_gubbins_groups,
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
    drop_resistance_target_hits,
    enrich_marker_hits_with_resistance,
)
from pipeline.PhylogenyWorkflow import (
    run_panaroo as run_panaroo_tool,
    ensure_panaroo_out,
    run_phylogeny as run_phylogeny_tool,
)
from pipeline.RunSteps import ProgressCallback, StatusCallback, StepProgress, planned_steps
from pipeline.RunValidation import (
    chewbbaca_samples_left_out,
    file_has_bytes as run_file_has_bytes,
    file_sha256,
    genome_set_change_error,
    MARKER_HITS_STAMP,
    marker_hits_stamp_matches,
    snippy_tree_path,
    write_genome_set_stamp,
    write_marker_hits_stamp,
    read_json_dict as run_read_json_dict,
    validate_amrfinder_outputs,
    validate_bakta_outputs,
    validate_assembly_qc_outputs,
    validate_chewbbaca_outputs,
    validate_completed_run,
    validate_core_alignment_outputs,
    validate_core_tree_outputs,
    population_output_files,
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
from pipeline.ContainerRuntime import database_provenance, environment_provenance
from pipeline.ToolEnv import (
    conda_tool_cmd,
    configured_image_ref,
    configured_assembly_qc_settings,
    configured_chewbbaca_settings,
    configured_fastgear_settings,
    configured_mobsuite_settings,
    configured_phylo_settings,
    current_container_context,
    make_container_context,
    phylo_outgroup_names,
    reset_container_context,
    resource_budget_text,
    set_container_context,
)
from pipeline.TreesIndex import (
    SNIPPY_RAW_REL,
    TREES_JSON,
    adopt_legacy_snippy_tree,
    invalidate_trees,
    snippy_reference_stem,
    write_trees_index,
)
from pipeline.Utils import assign_gene_ids

logger = setup_logger()

RUN_CONFIG_FILENAME = "run_config.json"
RUN_REPORT_FILENAME = "run_report.json"
RUN_STATE_FILENAME = "run_state.json"
# A per-genome step's progress rewrites run_state.json at most this often (seconds).
PROGRESS_FLUSH_INTERVAL = 1.0
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
    *MUTATION_COLUMNS,
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

# The Snippy tree (raw, unrooted) and the rooted tree index written from every tree.
PHYLOGENY_OUT_FILES = [
    SNIPPY_RAW_REL,
    Path("trees") / TREES_JSON,
]

BAKTA_REUSE_SUFFIXES = [
    ".json",
    ".gff3",
]


class GenomeSetChangedError(RuntimeError):
    """The output folder holds a run built from a different genome set."""


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
        "amrfinder": min(total, 32),
        "assembly_qc": min(total, 32),
        "mobsuite": min(total, 32),
        "chewbbaca": min(total, 32),
        "mlst": min(total, max(1, min(64, int(n_genomes)))),
        "diamond_makedb": min(total, 32),
        "diamond_blastx": min(total, 64),
        "panaroo": panaroo_threads,
        "gubbins": min(total, 64),
        "core_tree": min(total, 16),
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
    if validate_chewbbaca_outputs(CFG.output, CFG.genomes):
        return False
    # Calls made with another build of the schema (e.g. before the P01 rebuild) are redone.
    # Only the pipeline checks this; the GUI still opens runs whose calls predate the stamp.
    schema_issue = chewbbaca_schema_stamp_issue(CFG.output)
    if schema_issue:
        logger.warning("%s; rerunning chewBBACA.", schema_issue)
        return False
    return True


def _validate_dmnd_files(paths: list[Path]) -> bool:
    return bool(paths) and all(_file_has_bytes(path) for path in paths)


def _marker_dmnds_current(CFG: DataInput) -> bool:
    """Every marker .dmnd was built from the current marker FASTA (a rebuilt DB makes it stale)."""
    return all(
        dmnd_is_current(seq, CFG.output / "dmnd" / "queries" / f"{seq.stem}.dmnd")
        for seq in CFG.marker_sequences
    )


def _marker_hits_current(CFG: DataInput) -> bool:
    """Every marker hit TSV is at least as new as the .dmnd it searched."""
    return all(
        hit_tsv_is_current(
            CFG.output / "annotations" / "bakta_vs_queries" / f"bakta_annotations_nt-vs-{seq.stem}.tsv",
            CFG.output / "dmnd" / "queries" / f"{seq.stem}.dmnd",
        )
        for seq in CFG.marker_sequences
    )


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


def _expected_core_alignment_paths(CFG: DataInput) -> list[Path]:
    return [CFG.output / CORE_ALIGNMENT_DIR / rel for rel in CORE_ALIGNMENT_OUT_FILES]


def _planned_gubbins_runs(
    CFG: DataInput, plan: GroupPlan | None, summary: dict | None = None
) -> tuple[list[str], bool]:
    """(ids of the Gubbins groups expected to run, single mode). Before fastGEAR the
    groups are unknown: one run on all genomes stands in until the gubbins step starts
    and rewrites expected_outputs.json with the resolved groups; after the step its
    summary gives the groups that ran (a group can still be skipped for too few
    variable columns once its alignment is staged)."""
    if summary is not None:
        ran = [str(g.get("id")) for g in summary.get("groups") or [] if isinstance(g, dict) and g.get("ran")]
        return ran, summary.get("resolved_mode") == "single"
    if plan is None:
        enough = len(CFG.genomes) >= int(configured_phylo_settings(CFG.phylo_config)["gubbins_min_genomes"])
        return (["all"] if enough else []), True
    return [g.id for g in plan.groups if g.skip_reason is None], plan.single


def _expected_gubbins_paths(CFG: DataInput, plan: GroupPlan | None = None, summary: dict | None = None) -> list[Path]:
    if summary is not None:
        return [CFG.output / rel for rel in expected_gubbins_files(summary)]
    ran, single = _planned_gubbins_runs(CFG, plan)
    return [CFG.output / rel for rel in planned_gubbins_files(ran, single)]


def _expected_core_tree_paths(CFG: DataInput, plan: GroupPlan | None = None, summary: dict | None = None) -> list[Path]:
    tree_id = "core_recomb" if _planned_gubbins_runs(CFG, plan, summary)[0] else "core"
    return [
        CFG.output / IQTREE_DIR / f"{tree_id}.treefile",
        CFG.output / "trees" / f"{tree_id}.nwk",
    ]


def _expected_phylogeny_paths(CFG: DataInput) -> list[Path]:
    return [CFG.output / rel for rel in PHYLOGENY_OUT_FILES]


def _build_expected_outputs(
    CFG: DataInput,
    *,
    run_id: str,
    gubbins_plan: GroupPlan | None = None,
    gubbins_summary: dict | None = None,
    population: PopulationChoice | None = None,
) -> list[str]:
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

    for path in _expected_core_alignment_paths(CFG):
        expected.append(path)

    for path in _expected_gubbins_paths(CFG, gubbins_plan, gubbins_summary):
        expected.append(path)

    for path in _expected_core_tree_paths(CFG, gubbins_plan, gubbins_summary):
        expected.append(path)

    for path in _expected_phylogeny_paths(CFG):
        expected.append(path)

    # Population structure (fastGEAR or FastBAPS) from the ordered core alignment; none
    # when the tools image has neither method.
    method = population.method if population is not None else "fastgear"
    if method is not None:
        expected.extend(population_output_files(output, method).values())

    for path in _expected_chewbbaca_paths(CFG):
        expected.append(path)
    for path in _expected_mlst_paths(CFG):
        expected.append(path)
    for path in _expected_mobsuite_paths(CFG):
        expected.append(path)
    expected.append(output / "amrfinder" / "summary.json")
    expected.append(output / EVIDENCE_DIR / EVIDENCE_SUMMARY)
    for g in genomes:
        expected.append(evidence_tsv_path(output, g.stem))

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


def _build_run_inputs(
    CFG: DataInput, *, thread_plan: dict[str, int], population: PopulationChoice | None = None
) -> dict[str, object]:
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
        "amrfinder_threads": thread_plan["amrfinder"],
        "assembly_qc_threads": thread_plan["assembly_qc"],
        "mobsuite_threads": thread_plan["mobsuite"],
        "chewbbaca_threads": thread_plan["chewbbaca"],
        "mlst_threads": thread_plan["mlst"],
        "panaroo_threads": thread_plan["panaroo"],
        "gubbins_threads": thread_plan["gubbins"],
        "core_tree_threads": thread_plan["core_tree"],
        "phylogeny_snippy_threads": thread_plan["phylogeny_snippy"],
        "phylogeny_raxml_threads": thread_plan["phylogeny_raxml"],
        "marker_categories": CFG.marker_categories,
        "reference_locus_categories": CFG.reference_locus_categories,
        "chewbbaca": configured_chewbbaca_settings(CFG.chewbbaca_config),
        "mobsuite": configured_mobsuite_settings(CFG.mobsuite_config),
        "assembly_qc": configured_assembly_qc_settings(CFG.assembly_qc_config),
        "fastgear": CFG.fastgear_config or configured_fastgear_settings(),
        "phylo": configured_phylo_settings(CFG.phylo_config),
    }
    if population is not None:
        # The population structure method the run uses (PopulationStructure: the explorer
        # and the validators find its folder from this record).
        inputs["population_structure"] = population.record()
    marker_resistance_metadata = _build_marker_resistance_metadata_inputs(CFG.marker_sequences)
    if marker_resistance_metadata:
        inputs["marker_resistance_metadata"] = marker_resistance_metadata
    inputs["databases"] = _run_databases(CFG)
    return inputs


# The missing source of the I tag in a run without a compound database.
NO_COMPOUND_SOURCE = "a reference locus (compound) database"


def _run_databases(CFG: DataInput) -> dict[str, dict]:
    """The marker and compound databases the run searches and the catalogs' ones that are
    not installed, and the evidence sources each tag had (Bakta, AMRFinderPlus, MOB-suite;
    once the bakta_evidence step ran). The GUI shows a tag without any as not installed."""
    out = {
        "marker": run_database_record(CFG.marker_sequences, MARKER_DB_CATALOG_PATH),
        "compound": run_database_record(CFG.reference_locus_sequences, COMPOUND_DB_CATALOG_PATH),
    }
    evidence = evidence_record(CFG.output)
    if evidence is not None:
        if not out["compound"].get("used"):
            # I (overlaps a reference locus) has no source but a compound database.
            evidence["missing"] = {**(evidence.get("missing") or {}), "I": [NO_COMPOUND_SOURCE]}
        out["evidence"] = evidence
    return out


def _remove_stale_database_outputs(CFG: DataInput) -> None:
    """Delete the DIAMOND databases and hit tables of marker and compound databases this run
    no longer searches (removed since the folder's last run), so they are never read."""
    marker_stems = {p.stem for p in CFG.marker_sequences}
    locus_stems = {p.stem for p in CFG.reference_locus_sequences}
    stale: list[Path] = []
    for folder, stems in (
        (CFG.output / "dmnd" / "queries", marker_stems),
        (CFG.output / "dmnd" / "reference_loci", locus_stems),
    ):
        for path in folder.glob("*.dmnd*") if folder.is_dir() else ():
            if path.name.split(".dmnd", 1)[0] not in stems:
                stale.append(path)
    for folder, stems in (
        (CFG.output / "annotations" / "bakta_vs_queries", marker_stems),
        (CFG.output / "annotations" / "reference_loci", locus_stems),
    ):
        for path in folder.glob("*-vs-*.tsv") if folder.is_dir() else ():
            if path.stem.split("-vs-", 1)[1] not in stems:
                stale.append(path)
    for path in stale:
        try:
            path.unlink()
        except OSError as exc:
            logger.warning("Could not remove %s of a database this run no longer searches: %s", path, exc)
        else:
            logger.info("Removed %s: its database is no longer in this run.", path)
    # The marker hits' stamp lists only the databases whose hit tables are left.
    stamp_path = CFG.output / "annotations" / "bakta_vs_queries" / MARKER_HITS_STAMP
    stamp = _read_json_dict(stamp_path)
    if stamp and isinstance(stamp.get("markers"), list):
        kept = [name for name in stamp["markers"] if name in marker_stems]
        if kept != stamp["markers"]:
            _write_json(stamp_path, {**stamp, "markers": kept})
            logger.info("Dropped %s from %s.", ", ".join(sorted(set(stamp["markers"]) - set(kept))), stamp_path)


def _database_warnings(databases: dict[str, dict]) -> list[dict[str, str]]:
    """Log each database that is not installed; the run-level warnings about the catalog
    entries whose file is missing, under the step that would have searched them."""
    warnings: list[dict[str, str]] = []
    for kind, step in (("marker", "diamond_markers"), ("compound", "diamond_reference_loci")):
        record = databases.get(kind) or {}
        missing: list[str] = []
        for db in record.get("not_installed") or []:
            category = f" ({db['category']})" if db.get("category") else ""
            logger.warning(f"{kind.capitalize()} database not installed, skipped: {db['name']}{category}: {db['path']}")
            missing.append(f"{db['name']}{category}")
        if missing:
            warnings.append({
                "step": step,
                "code": "databases_not_installed",
                "message": "Not installed, not searched: " + ", ".join(missing),
            })
        if not record.get("used"):
            # The normal case: the tags come from Bakta, AMRFinderPlus and MOB-suite.
            what = "marker search" if kind == "marker" else "reference locus search"
            logger.info(f"No {kind} database is installed: the {what} is skipped.")
    return warnings


def _build_run_config(
    CFG: DataInput,
    *,
    sample_set_id: str,
    run_id: str,
    annotation_id_strategy: str,
    thread_plan: dict[str, int],
    environment: dict | None = None,
    environment_history: list[dict] | None = None,
    population: PopulationChoice | None = None,
) -> dict:
    inputs = _build_run_inputs(CFG, thread_plan=thread_plan, population=population)
    config = {
        "created_at": _now_iso(),
        "bactolomew_version": __version__,
        "sample_set_id": sample_set_id,
        "run_id": run_id,
        "annotation_id_strategy": annotation_id_strategy,
        "inputs": inputs,
        "thread_plan": thread_plan,
        "output_dir": str(CFG.output),
    }
    if environment is not None:
        config["environment"] = environment
    if environment_history:
        config["environment_history"] = environment_history
    return config


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
    warnings: list[dict[str, str]],
    environment: dict | None = None,
    environment_history: list[dict] | None = None,
    population: PopulationChoice | None = None,
) -> dict:
    report: dict = {
        "start_time": start_time,
        "end_time": _now_iso(),
        "status": status,
        "bactolomew_version": __version__,
        "outputs_validated": bool(outputs_validated),
        "sample_set_id": sample_set_id,
        "run_id": run_id,
        "inputs": _build_run_inputs(CFG, thread_plan=thread_plan, population=population),
        "steps": steps,
        "warnings": list(warnings),
        "tools": {
            "used": sorted(tools_used),
            "versions": _collect_versions(tools_used),
        },
        "output_dir": str(CFG.output),
        "output_files": _collect_output_files(CFG.output),
    }
    if environment is not None:
        report["environment"] = environment
    if environment_history:
        report["environment_history"] = environment_history
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
    inputs: dict[str, object],
    start_time: str,
    steps: dict,
    status: str,
    current_step: str | None,
    error: str | None,
    outputs_validated: bool,
    validation_errors: list[str],
    report_written: bool,
    warnings: list[dict[str, str]],
) -> dict:
    payload: dict[str, object] = {
        "start_time": start_time,
        "updated_at": _now_iso(),
        "status": status,
        "outputs_validated": bool(outputs_validated),
        "sample_set_id": sample_set_id,
        "run_id": run_id,
        "current_step": current_step,
        "planned_steps": planned_steps(
            (inputs.get("population_structure") or {"method": "fastgear"}).get("method")
        ),
        "inputs": inputs,
        "steps": steps,
        "warnings": list(warnings),
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
        "amrfinder": _cmds("amrfinder", ["--version"]),
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
        "iqtree2": _cmds("iqtree2", ["--version"]),
        "fastbaps": _cmds("fastbaps", ["-e", 'cat(as.character(packageVersion("fastbaps")))'],
                          executables=("Rscript",)),
        "minimap2": _cmds("minimap2", ["--version"]),
        "docker": _cmds("docker", ["--version"]),
        "podman": [["podman", "--version"]],
        "singularity": [["singularity", "--version"], ["apptainer", "--version"]],
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


def run_bakta(CFG: DataInput, *, bakta_threads: int, progress: ProgressCallback = None) -> bool:
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
            progress=progress,
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


def run_phylogeny(
    CFG: DataInput,
    *,
    snippy_threads: int,
    raxml_threads: int,
    progress: ProgressCallback = None,
):
    try:
        logger.info("Running Snippy and RAxML-NG for the whole-genome SNP tree.")
        return run_phylogeny_tool(
            genome_files=CFG.genomes,
            output_path=CFG.output,
            logger=logger,
            snippy_threads=snippy_threads,
            raxml_threads=raxml_threads,
            reference_genome=CFG.snippy_reference,
            progress=progress,
        )
    except Exception as exc:
        logger.error("Phylogeny workflow failed: %s", exc)
        raise


def run_chewbbaca(
    CFG: DataInput,
    *,
    chewbbaca_threads: int,
    status: StatusCallback = None,
) -> Path | None:
    """The chewBBACA step; a Stop during its cgMLST install raises RunStoppedError."""
    from pipeline.CgmlstInstall import Cancelled

    try:
        logger.info("Running chewBBACA AlleleCall on genomes.")
        return run_chewbbaca_tool(
            genome_files=CFG.genomes,
            output_path=CFG.output,
            logger=logger,
            threads=chewbbaca_threads,
            config=CFG.chewbbaca_config,
            status=status,
            should_cancel=lambda: _stop_requested(CFG.output),
        )
    except Cancelled as exc:
        raise RunStoppedError("Run stopped by user") from exc
    except Exception as exc:
        logger.error("chewBBACA workflow failed: %s", exc)
        raise


def run_mlst(
    CFG: DataInput,
    *,
    mlst_threads: int,
    progress: ProgressCallback = None,
) -> Path:
    try:
        logger.info("Running mlst sequence typing on genomes.")
        return run_mlst_tool(
            genome_files=CFG.genomes,
            output_path=CFG.output,
            logger=logger,
            threads=mlst_threads,
            progress=progress,
        )
    except Exception as exc:
        logger.error("mlst workflow failed: %s", exc)
        raise


def run_mobsuite(
    CFG: DataInput,
    *,
    mobsuite_threads: int,
    progress: ProgressCallback = None,
):
    try:
        logger.info("Running MOB-suite plasmid reconstruction and typing on genomes.")
        return run_mobsuite_tool(
            genome_files=CFG.genomes,
            output_path=CFG.output,
            logger=logger,
            threads=mobsuite_threads,
            config=CFG.mobsuite_config,
            progress=progress,
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
    progress: ProgressCallback = None,
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
            progress=progress,
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
                # The hit TSVs are keyed by AnnotationID, so they are valid only for
                # the exact query FASTA (IDs and sequences) they were searched from.
                query_sha = file_sha256(fa_path)
                query_matches = marker_hits_stamp_matches(annotations_folder, query_sha)
                if (
                    reuse_marker_hits
                    and query_matches
                    and _all_files_exist(_expected_marker_hit_paths(CFG))
                ):
                    logger.info("Reusing existing DIAMOND marker hit TSVs under %s.", annotations_folder)
                else:
                    if not query_matches and any(p.exists() for p in _expected_marker_hit_paths(CFG)):
                        logger.warning(
                            "Existing DIAMOND marker hits under %s were not searched from the current "
                            "annotation set (AnnotationIDs or sequences differ, or no record of it); "
                            "recomputing them.",
                            annotations_folder,
                        )
                    run_diamond_annotations(
                        genome_files=[fa_path],
                        db_files=query_db_files,
                        output_path=annotations_folder,
                        logger=logger,
                        threads=diamond_blastx_threads,
                        sanitize=False,
                        # Only a changed query invalidates every TSV; a rebuilt marker DB makes
                        # only its own TSV stale, which run_diamond_annotations checks per DB.
                        force=not query_matches,
                    )
                    if _all_files_exist(_expected_marker_hit_paths(CFG)):
                        write_marker_hits_stamp(
                            annotations_folder,
                            query_sha256=query_sha,
                            genomes=CFG.genomes,
                            marker_sequences=CFG.marker_sequences,
                        )
                    logger.info("DIAMOND on Bakta annotations vs MarkerSequence DBs completed.")

    else:
        logger.warning("No Bakta annotations parsed; nothing to export to FASTA.")
        bakta_df = bakta_df[ANNOTATION_COLUMNS]

    diamond_annotations = read_annotations(CFG.output / "annotations" / "bakta_vs_queries")
    if not diamond_annotations.empty:
        logger.info(f"Total DIAMOND annotations parsed (raw): {len(diamond_annotations)}")
        diamond_annotations = drop_resistance_target_hits(
            diamond_annotations,
            marker_resistance_metadata,
            logger=logger,
        )
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
        bakta_df = attach_resistance_mutations(bakta_df, CFG.output)

    # No CSVs are written; annotations are kept in-memory only.
    return bakta_df


def _check_phylo_settings(CFG: DataInput, settings: dict) -> None:
    """Refuse, before anything is written, phylo settings that name unknown samples."""
    samples = {p.stem for p in CFG.genomes}
    unknown = [name for name in phylo_outgroup_names(settings) if name not in samples]
    if unknown:
        raise ValueError(f"phylo.outgroup: not samples of this run: {', '.join(unknown)}")
    order_reference = str(settings.get("order_reference") or "").strip()
    if order_reference and order_reference not in samples and not Path(order_reference).expanduser().is_file():
        raise ValueError(
            f"phylo.order_reference: {order_reference!r} is neither a sample of this run nor a genome FASTA file"
        )


def _core_alignment_info(summary: dict) -> dict:
    keys = (
        "method",
        "ordering_genome",
        "genes_total",
        "genes_placed",
        "genes_appended",
        "genes_reversed",
        "genes_multi_copy",
        "diversity",
        "predicted_gubbins_dropouts",
    )
    info = {("ordering_method" if k == "method" else k): summary.get(k) for k in keys}
    info["alignment"] = (summary.get("output") or {}).get("path")
    info["sha256"] = (summary.get("output") or {}).get("sha256")
    info["source_alignment"] = ((summary.get("inputs") or {}).get("source_alignment") or {}).get("path")
    return info


def _gubbins_info(summary: dict) -> dict:
    params = summary.get("params") or {}
    return {
        "mode": summary.get("mode"),
        "resolved_mode": summary.get("resolved_mode"),
        "groups": [
            {
                **{key: group.get(key) for key in ("id", "genomes", "ran", "skip_reason", "reused", "converged", "iterations_run")},
                "longest_branch_share": (group.get("longest_branch_share") or {}).get("share"),
            }
            for group in summary.get("groups") or []
            if isinstance(group, dict)
        ],
        "skipped": not summary.get("ran", True),
        "skip_reason": summary.get("skip_reason"),
        "converged": summary.get("converged"),
        "iterations_run": summary.get("iterations_run"),
        "max_iterations": summary.get("max_iterations"),
        "seed": params.get("seed"),
        "threads": params.get("threads"),
        "converge_method": params.get("converge_method"),
        "outgroup": params.get("outgroup"),
        "diversity": summary.get("diversity"),
        "samples_left_out": summary.get("samples_left_out") or [],
        "input_sha256": summary.get("input_sha256"),
        "warnings": [w.get("code") for w in summary.get("warnings") or []],
    }


def _refit_info(refit: dict) -> dict:
    """Step info of the core tree's second branch length set (including recombination)."""
    return {
        key: refit.get(key)
        for key in ("status", "skip_reason", "message", "error", "model", "threads", "wall_time_s", "tree_length_ratio")
        if refit.get(key) is not None
    }


def _core_tree_info(output: Path, run_stamp: dict) -> dict:
    params = run_stamp.get("params") or {}
    tree_id = params.get("tree_id")
    return {
        "tree_id": tree_id,
        "file": str(output / "trees" / f"{tree_id}.nwk") if tree_id else None,
        "model": run_stamp.get("model"),
        "support_type": params.get("support_type"),
        "replicates": params.get("bootstraps"),
        "fast_mode": params.get("fast_mode"),
        "seed": params.get("seed"),
        "threads": params.get("threads"),
        "sites": run_stamp.get("sites"),
        "variable_sites": run_stamp.get("variable_sites"),
        "ufboot_converged": run_stamp.get("ufboot_converged"),
        "composition_failed": run_stamp.get("composition_failed"),
        "identical_collapsed": run_stamp.get("identical_collapsed"),
        "rf_to_gubbins": run_stamp.get("rf_to_gubbins"),
    }


def _external_order_reference(CFG: DataInput) -> list[Path]:
    """phylo.order_reference when it names a genome FASTA (minimap2 must see it)."""
    value = str(configured_phylo_settings(CFG.phylo_config).get("order_reference") or "").strip()
    path = Path(value).expanduser() if value else None
    return [path] if path is not None and path.is_file() else []


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
            *_external_order_reference(CFG),
        ],
    )
    _token = set_container_context(_ctx)
    if _ctx.runtime != "native":
        logger.info(f"Tool runtime: {_ctx.runtime} (image {_ctx.image or _ctx.sif_path})")
    try:
        logger.info(f"Resource budget: {resource_budget_text()}; this run plans {CFG.total_cores} cores")
    except Exception as exc:  # noqa: BLE001 - the log line is informational
        logger.warning(f"Resource budget could not be read: {exc}")
    try:
        # Provenance: one image inspect and one small read of the image's tool versions.
        environment = environment_provenance(_ctx, configured_image=configured_image_ref())
        _add_database_provenance(environment, CFG)
        return _run_workflow_impl(CFG, environment=environment)
    finally:
        reset_container_context(_token)


def _add_database_provenance(environment: dict, CFG: DataInput) -> None:
    """Add the databases outside the image (Bakta, CheckM2, an external mob_suite DB) to
    `environment["databases"]`, next to the in-image ones environment_provenance read."""
    try:
        host = database_provenance(
            bakta_db=CFG.bakta_db,
            checkm2_db=configured_assembly_qc_settings(CFG.assembly_qc_config).get("checkm2_database_path"),
            mobsuite_db_dir=configured_mobsuite_settings(CFG.mobsuite_config).get("database_directory"),
        )
    except Exception as exc:
        host = {"error": f"database provenance failed: {exc}"}
    if host:
        environment["databases"] = {**(environment.get("databases") or {}), **host}


# The parts of an `environment` block that identify what computed a run's outputs; a resume
# under a different one keeps the earlier block in `environment_history`.
_ENVIRONMENT_IDENTITY_KEYS = ("runtime", "image", "image_id", "sif_path", "sif_sha256", "tool_versions", "databases")


def _environment_identity(environment: dict | None) -> str:
    if not isinstance(environment, dict):
        return ""
    return json.dumps({k: environment.get(k) for k in _ENVIRONMENT_IDENTITY_KEYS}, sort_keys=True, default=str)


def _environment_history(existing_run_cfg: dict, environment: dict | None) -> list[dict]:
    """Earlier environments of this output folder, oldest first. run_config.json and
    run_report.json describe the current invocation; when a run folder is resumed with
    another image or other databases, outputs reused from before were computed by an
    earlier environment, so each distinct earlier one is kept with the time it was
    recorded (`recorded_at`, the earlier run_config.json's created_at)."""
    history = [h for h in existing_run_cfg.get("environment_history") or [] if isinstance(h, dict)]
    previous = existing_run_cfg.get("environment")
    if existing_run_cfg.get("created_at") and previous is None and not history:
        # A folder first run before run provenance existed: what made it is unknown.
        history.append({"recorded_at": existing_run_cfg.get("created_at"), "environment": None,
                        "note": "this folder's earlier runs recorded no environment"})
    elif isinstance(previous, dict) and environment is not None:
        seen = {_environment_identity(h.get("environment")) for h in history}
        identity = _environment_identity(previous)
        if identity != _environment_identity(environment) and identity not in seen:
            history.append({"recorded_at": existing_run_cfg.get("created_at"), "environment": previous})
    return history


def _run_workflow_impl(CFG: DataInput, *, environment: dict | None = None):
    # Refuse an output folder built from another genome set before anything is
    # written: its cohort-wide outputs and annotation IDs belong to that set.
    genome_set_error = genome_set_change_error(CFG.output, CFG.genomes)
    if genome_set_error:
        logger.error(genome_set_error)
        raise GenomeSetChangedError(genome_set_error)
    validate_chewbbaca_settings(CFG.chewbbaca_config)
    phylo_settings = configured_phylo_settings(CFG.phylo_config)
    _check_phylo_settings(CFG, phylo_settings)
    # The population structure method, from what the tools image holds; a method set for
    # the run or in the config that the image lacks stops the run here, before anything is
    # written (the default, FastBAPS, falls back to fastGEAR on an image without it).
    population = configured_population_choice(
        phylo_settings, explicit=population_method_explicit(CFG.phylo_config)
    )
    logger.info(
        "Population structure: %s (phylo.population_method %s; in the tools image: %s).",
        method_label(population.method) if population.method else "none",
        population.requested,
        ", ".join(method_label(m) for m in population.available) if population.available is not None
        else "not checked",
    )
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
    environment_history = _environment_history(existing_run_cfg, environment)

    run_config = _build_run_config(
        CFG,
        sample_set_id=sample_set_id,
        run_id=run_id,
        annotation_id_strategy=annotation_id_strategy,
        thread_plan=thread_plan,
        environment=environment,
        environment_history=environment_history,
        population=population,
    )
    _write_json(CFG.output / RUN_CONFIG_FILENAME, run_config)

    try:
        (CFG.output / RUN_STOP_FILENAME).unlink()
    except FileNotFoundError:
        pass
    except Exception:
        pass

    # The arguments of the last expected_outputs.json write (a later write keeps them).
    expected_args: dict = {}

    def _write_expected_outputs(gubbins_plan: GroupPlan | None = None, gubbins_summary: dict | None = None) -> None:
        if gubbins_plan is not None or gubbins_summary is not None:
            expected_args.update(gubbins_plan=gubbins_plan, gubbins_summary=gubbins_summary)
        _write_json(
            CFG.output / "expected_outputs.json",
            {
                "created_at": _now_iso(),
                "sample_set_id": sample_set_id,
                "run_id": run_id,
                "output_dir": str(CFG.output),
                "files": _build_expected_outputs(CFG, run_id=run_id, population=population, **expected_args),
            },
        )

    _write_expected_outputs()
    _write_run_lock(CFG.output, run_id=run_id, start_time=start_time)
    steps: dict = {}
    tools_used: set[str] = set()
    error: Exception | None = None
    error_text: str | None = None
    current_step: str | None = None
    outputs_validated = False
    validation_errors: list[str] = []
    # Run-level warnings ({step, code, message}): conditions that do not stop the run but
    # must be seen (Gubbins not converged, lineage structure, a tree left out, ...).
    run_warnings: list[dict[str, str]] = []

    def _add_run_warnings(items) -> None:
        for item in items or []:
            if isinstance(item, dict) and item.get("code") and dict(item) not in run_warnings:
                run_warnings.append(dict(item))

    # Warnings of the last trees.json write: each write replaces them, so a tree that was
    # stale before the Snippy step and rebuilt by it leaves no warning behind.
    index_warnings: list[dict[str, str]] = []

    def _write_trees_index() -> dict:
        """Root the raw trees and rewrite trees/trees.json; returns the default tree's rooting."""
        for item in index_warnings:
            if item in run_warnings:
                run_warnings.remove(item)
        index_warnings.clear()
        others = list(run_warnings)
        payload, new_warnings = write_trees_index(
            CFG.output, CFG.genomes, phylo_settings, logger, run_warnings=others
        )
        for item in new_warnings:
            if dict(item) not in run_warnings:
                index_warnings.append(dict(item))
        _add_run_warnings(index_warnings)
        return payload.get("rooting") or {}

    # The run's inputs do not change while it runs: built once, not on every write.
    run_state_inputs = _build_run_inputs(CFG, thread_plan=thread_plan, population=population)
    # run_state.json is written whole (temp file, then rename) under this lock, so a
    # per-genome progress update never interleaves with a step boundary's write.
    state_lock = threading.RLock()
    # The running steps' "progress" (genome counts, the chewBBACA step's cgMLST install).
    step_progress = StepProgress(steps, state_lock, lambda: _flush_run_state("running"),
                                 interval=PROGRESS_FLUSH_INTERVAL)

    def _flush_run_state(
        status: str,
        *,
        error_override: str | None = None,
        report_written: bool = False,
    ) -> None:
        with state_lock:
            step_progress.written()
            try:
                _write_json(
                    CFG.output / RUN_STATE_FILENAME,
                    _build_run_state(
                        CFG,
                        sample_set_id=sample_set_id,
                        run_id=run_id,
                        inputs=run_state_inputs,
                        start_time=start_time,
                        steps=steps,
                        status=status,
                        current_step=current_step,
                        error=error_override,
                        outputs_validated=outputs_validated,
                        validation_errors=validation_errors,
                        report_written=report_written,
                        warnings=run_warnings,
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
        with state_lock:
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
        with state_lock:
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
    _remove_stale_database_outputs(CFG)
    for f in CFG.genomes:
        logger.info(f"Genome: {f.name}")
    for f in all_marker_sequences:
        logger.info(f"MarkerSequence: {f.name}")
    _add_run_warnings(_database_warnings(run_state_inputs["databases"]))

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
                ran_bakta = run_bakta(
                    CFG, bakta_threads=thread_plan["bakta"], progress=step_progress.genome_progress("bakta")
                )
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
        chewbbaca_out_existing = ensure_chewbbaca_out(CFG.output, logger)
        mlst_paths_existing = _expected_mlst_paths(CFG)
        mobsuite_paths_existing = _expected_mobsuite_paths(CFG)
        query_db_files_existing = _expected_marker_db_paths(CFG)
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
                mlst_out = run_mlst(CFG, mlst_threads=thread_plan["mlst"], progress=step_progress.genome_progress("mlst"))
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

        # AMRFinderPlus point mutations (needs the Bakta outputs and the MLST species)
        _check_stop()
        _step_start("amrfinder_mutations")
        amrfinder_plan = plan_amrfinder(
            CFG.output,
            bakta_db=CFG.bakta_db,
            schema_dir=configured_chewbbaca_settings(CFG.chewbbaca_config).get("schema_dir"),
        )
        amrfinder_issue = amrfinder_reuse_issue(CFG.output, CFG.genomes, amrfinder_plan)
        if not amrfinder_issue and not validate_amrfinder_outputs(CFG.output, CFG.genomes):
            amrfinder_summary = read_amrfinder_summary(CFG.output) or {}
            _finish_step_success(
                "amrfinder_mutations",
                reason="reused_existing_output",
                output=str(CFG.output / "amrfinder"),
                summary=str(CFG.output / "amrfinder" / "summary.json"),
                organism=amrfinder_plan.organism,
                organism_source=amrfinder_summary.get("organism_source"),
                db_version=amrfinder_plan.db_version or None,
                point_rows=amrfinder_summary.get("point_rows"),
                threads=thread_plan["amrfinder"],
            )
        else:
            if (CFG.output / "amrfinder").is_dir():
                logger.warning("Existing AMRFinderPlus outputs cannot be reused (%s); rerunning.",
                               amrfinder_issue or "incomplete or invalid")
            try:
                if amrfinder_plan.db_dir:
                    tools_used.add("amrfinder")
                amrfinder_result = run_amrfinder_mutations(
                    CFG.genomes,
                    CFG.output,
                    logger,
                    plan=amrfinder_plan,
                    threads=thread_plan["amrfinder"],
                    progress=step_progress.genome_progress("amrfinder_mutations"),
                )
            except Exception as e:
                # Optional step: the run goes on without point mutation calls.
                logger.warning("AMRFinderPlus point mutations skipped after an error: %s", e)
                _step_end("amrfinder_mutations", "skipped", reason=f"error: {e}", outputs_validated=False)
            else:
                _finish_step_success(
                    "amrfinder_mutations",
                    validator_errors=validate_amrfinder_outputs(CFG.output, CFG.genomes),
                    output=str(amrfinder_result.output_dir),
                    summary=str(amrfinder_result.summary_json),
                    organism=amrfinder_result.organism,
                    organism_source=amrfinder_plan.note,
                    db_version=amrfinder_plan.db_version or None,
                    point_rows=amrfinder_result.point_rows,
                    failed_genomes=list(amrfinder_result.failed) or None,
                    threads=thread_plan["amrfinder"],
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
                mobsuite_result = run_mobsuite(
                    CFG, mobsuite_threads=thread_plan["mobsuite"], progress=step_progress.genome_progress("mobsuite")
                )
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

        # Gene evidence from Bakta, AMRFinderPlus and MOB-suite (the tags of the annotations)
        _check_stop()
        _step_start("bakta_evidence")
        if not ran_bakta:
            _step_end("bakta_evidence", "skipped", reason="no Bakta annotation", outputs_validated=False)
        else:
            evidence_issue = evidence_reuse_issue(CFG.output, CFG.genomes)
            if not evidence_issue and not validate_evidence_outputs(CFG.output, CFG.genomes):
                _finish_step_success(
                    "bakta_evidence",
                    reason="reused_existing_output",
                    output=str(CFG.output / EVIDENCE_DIR),
                    sources=(evidence_record(CFG.output) or {}).get("tags"),
                )
            else:
                if (CFG.output / EVIDENCE_DIR).is_dir():
                    logger.info("Rebuilding the gene evidence (%s).", evidence_issue or "incomplete or invalid")
                try:
                    evidence_summary = build_evidence(CFG.output, CFG.genomes, logger)
                except Exception as e:
                    _step_end("bakta_evidence", "failed", error=str(e))
                    raise
                _finish_step_success(
                    "bakta_evidence",
                    validator_errors=validate_evidence_outputs(CFG.output, CFG.genomes),
                    output=str(CFG.output / EVIDENCE_DIR),
                    sources=evidence_summary.get("tags"),
                    genes=evidence_summary.get("genes"),
                    amrfinder_unplaced=(evidence_summary.get("amrfinder") or {}).get("unplaced"),
                )
            # run_config.json names the evidence sources the run had (the explorer reads them).
            run_config["inputs"]["databases"] = _run_databases(CFG)
            _write_json(CFG.output / RUN_CONFIG_FILENAME, run_config)

        # Panaroo
        _check_stop()
        _step_start("panaroo")
        if not ran_bakta:
            _step_end("panaroo", "skipped", reason="bakta_not_run", outputs_validated=False, validation_errors=[])
            panaroo_out = None
        elif panaroo_gene_presence_existing and not validate_panaroo_outputs(CFG.output, CFG.genomes):
            panaroo_out = panaroo_out_existing
            _finish_step_success(
                "panaroo",
                validator_errors=validate_panaroo_outputs(CFG.output, CFG.genomes),
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
                    validator_errors=validate_panaroo_outputs(CFG.output, CFG.genomes),
                    output=str(panaroo_out),
                    threads=thread_plan["panaroo"],
                )

        # Genome-ordered core alignment: the one input of Gubbins, the core tree and fastGEAR.
        # Each of these steps is reused only when its input stamp still matches (sha256 of
        # its input, its parameters, the tool version); otherwise it reruns with the reason.
        _check_stop()
        _step_start("core_alignment")
        if not panaroo_out:
            _step_end(
                "core_alignment",
                "failed",
                reason="no_panaroo_output",
                outputs_validated=False,
                validation_errors=["Panaroo output required before the core alignment can be ordered."],
            )
            raise RuntimeError("Panaroo output required before the core alignment can be ordered.")
        core_alignment_issue = core_alignment_reuse_issue(CFG.output, CFG.genomes, phylo_settings)
        if core_alignment_issue is None and not validate_core_alignment_outputs(CFG.output, CFG.genomes):
            ordering = read_ordering_summary(CFG.output) or {}
            _add_run_warnings(ordering.get("warnings"))
            _finish_step_success(
                "core_alignment",
                reason="reused_existing_output",
                **_core_alignment_info(ordering),
            )
        else:
            logger.info("Ordering the core alignment: %s.", core_alignment_issue or "existing output is invalid")
            previous_sha = ordered_alignment_sha256(CFG.output)
            try:
                core_alignment_result = run_core_alignment(CFG.output, CFG.genomes, phylo_settings, logger)
            except Exception as e:
                invalidate_trees(CFG.output, ("core_recomb", "core", "gubbins_final"))
                _step_end("core_alignment", "failed", error=str(e))
                raise
            else:
                # The trees downstream are deleted only when the ordered alignment changed;
                # an identical alignment keeps Gubbins and the core tree reusable.
                if core_alignment_result.summary.get("output", {}).get("sha256") != previous_sha:
                    invalidate_trees(CFG.output, ("core_recomb", "core", "gubbins_final"))
                else:
                    logger.info("The ordered core alignment is unchanged; the trees built on it stay valid.")
                _add_run_warnings(core_alignment_result.warnings)
                _finish_step_success(
                    "core_alignment",
                    validator_errors=validate_core_alignment_outputs(CFG.output, CFG.genomes),
                    **_core_alignment_info(core_alignment_result.summary),
                )

        # Population structure (fastGEAR or FastBAPS, pipeline/PopulationStructure.py) on the
        # ordered core alignment: the lineages Gubbins runs per and the explorer's clades and
        # clusters. The step keeps the id fastgear_clades for both methods. A failure does
        # not stop the phylogeny: Gubbins runs once on all genomes (phylo.gubbins_mode auto),
        # the core tree and Snippy follow, and the run then fails with the error. Without a
        # method (the tools image has neither) the step is skipped and Gubbins runs once.
        _check_stop()
        _step_start("fastgear_clades")
        population_method = population.method
        population_name = method_label(population_method) if population_method else "population structure"
        fastgear_config = CFG.fastgear_config or configured_fastgear_settings()
        fastgear_failure: Exception | None = None
        fastgear_error: str | None = None
        if population_method is None:
            fastgear_error = population.note
            logger.warning("Population structure skipped: %s; Gubbins runs once on all genomes.", fastgear_error)
            _step_end("fastgear_clades", "skipped", reason=fastgear_error, method=None, outputs_validated=False)
        else:
            population_files = {
                key: str(path) for key, path in population_output_files(CFG.output, population_method).items()
            }
            if population_method == "fastbaps":
                population_issue = fastbaps_reuse_issue(CFG.output, int(phylo_settings["seed"]))
            else:
                population_issue = fastgear_reuse_issue(CFG.output, fastgear_config)
            if population_issue is None and not validate_fastgear_outputs(CFG.output, CFG.genomes, population_method):
                _finish_step_success(
                    "fastgear_clades", reason="reused_existing_output", method=population_method, **population_files
                )
            else:
                logger.info("Running %s: %s.", population_name, population_issue or "existing output is invalid")
                try:
                    if population_method == "fastbaps":
                        tools_used.add("fastbaps")
                        run_fastbaps_clades(
                            CFG.output,
                            expected_samples=[p.stem for p in CFG.genomes],
                            logger=logger,
                            seed=int(phylo_settings["seed"]),
                        )
                    else:
                        tools_used.add(fastgear_engine(current_container_context() or make_container_context()))
                        run_fastgear_clades(
                            CFG.output,
                            expected_samples=[p.stem for p in CFG.genomes],
                            logger=logger,
                            fastgear_config=fastgear_config,
                        )
                except RunStoppedError:
                    raise
                except Exception as e:
                    _step_end("fastgear_clades", "failed", error=str(e), method=population_method)
                    fastgear_failure = e
                    logger.error("%s failed (%s); the phylogeny steps continue without lineages.", population_name, e)
                else:
                    try:
                        _finish_step_success(
                            "fastgear_clades",
                            validator_errors=validate_fastgear_outputs(CFG.output, CFG.genomes, population_method),
                            method=population_method,
                            **population_files,
                        )
                    except RuntimeError as e:
                        fastgear_failure = e
        fastgear_lineages = None
        if fastgear_failure is not None:
            fastgear_error = str(fastgear_failure)
        elif population_method is not None:
            try:
                fastgear_lineages = read_fastgear_lineages(CFG.output)
            except Exception as e:
                fastgear_error = f"{population_method}/sample_clades.tsv is unreadable: {e}"

        # Gubbins per fastGEAR lineage (or once on all genomes) and the masked core alignment.
        _check_stop()
        _step_start("gubbins")
        try:
            gubbins_plan = resolve_gubbins_groups(
                [p.stem for p in CFG.genomes], phylo_settings, fastgear_lineages, fastgear_error,
                lineage_source=population_name,
            )
        except RuntimeError as e:
            _step_end("gubbins", "failed", error=str(e))
            raise
        _write_expected_outputs(gubbins_plan)
        gubbins_issue = gubbins_reuse_issue(
            CFG.output,
            phylo_settings,
            threads=thread_plan["gubbins"],
            lineages=fastgear_lineages,
            lineage_error=fastgear_error,
            lineage_source=population_name,
        )
        if gubbins_issue is None and not validate_gubbins_outputs(CFG.output, CFG.genomes):
            gubbins_summary = read_gubbins_summary(CFG.output) or {}
            _write_expected_outputs(gubbins_summary=gubbins_summary)
            _add_run_warnings(gubbins_summary.get("warnings"))
            _finish_step_success(
                "gubbins",
                reason="reused_existing_output",
                **_gubbins_info(gubbins_summary),
            )
        else:
            logger.info("Running Gubbins: %s.", gubbins_issue or "existing output is invalid")
            # The trees built on the Gubbins output follow its sha256 (core_tree reuse check
            # and trees.json); they are deleted here only when Gubbins fails.
            try:
                tools_used.add("run_gubbins.py")
                gubbins_result = run_gubbins_tool(
                    CFG.output,
                    phylo_settings,
                    logger,
                    threads=thread_plan["gubbins"],
                    lineages=fastgear_lineages,
                    lineage_error=fastgear_error,
                    lineage_source=population_name,
                )
            except Exception as e:
                invalidate_trees(CFG.output, ("core_recomb", "gubbins_final"))
                _step_end("gubbins", "failed", error=str(e))
                raise
            else:
                write_genome_set_stamp(CFG.output / "gubbins", CFG.genomes, step="gubbins")
                _write_expected_outputs(gubbins_summary=gubbins_result.summary)
                _add_run_warnings(gubbins_result.warnings)
                _finish_step_success(
                    "gubbins",
                    validator_errors=validate_gubbins_outputs(CFG.output, CFG.genomes),
                    **_gubbins_info(gubbins_result.summary),
                )

        # Core genome tree: IQ-TREE on the recombination-masked core alignment (the default tree).
        _check_stop()
        _step_start("core_tree")
        core_tree_issue = core_tree_reuse_issue(CFG.output, phylo_settings, threads=thread_plan["core_tree"])
        if core_tree_issue is None and refresh_gubbins_notes(CFG.output):
            logger.info("Core tree reused; its Gubbins notes were updated from gubbins/summary.json.")
        core_tree_stamp = read_core_tree_stamp(CFG.output) or {}
        if core_tree_issue is None:
            _add_run_warnings(core_tree_stamp.get("warnings"))
        else:
            logger.info("Building the core tree: %s.", core_tree_issue)
            try:
                tools_used.add("iqtree2")
                run_core_tree(CFG.output, phylo_settings, logger, threads=thread_plan["core_tree"])
            except Exception as e:
                _step_end("core_tree", "failed", error=str(e))
                raise
            core_tree_stamp = read_core_tree_stamp(CFG.output) or {}
            _add_run_warnings(core_tree_stamp.get("warnings"))
        # Second branch length set of the same tree: refit on the unmasked alignment (own
        # reuse stamp; a failure is a warning, the core tree stands).
        refit = run_recombination_refit(CFG.output, phylo_settings, logger, threads=thread_plan["core_tree"])
        if refit.get("status") in ("built", "failed"):
            tools_used.add("iqtree2")
        _add_run_warnings(refit.get("warnings"))
        _write_trees_index()
        _finish_step_success(
            "core_tree",
            validator_errors=validate_core_tree_outputs(CFG.output, CFG.genomes),
            **({"reason": "reused_existing_output"} if core_tree_issue is None else {}),
            **_core_tree_info(CFG.output, core_tree_stamp),
            recombination_branch_lengths=_refit_info(refit),
        )

        # Snippy whole-genome SNP tree (secondary, not recombination-corrected).
        _check_stop()
        _step_start("phylogeny")
        if not validate_phylogeny_outputs(CFG.output, CFG.genomes):
            adopt_legacy_snippy_tree(CFG.output, logger)
            rooting = _write_trees_index()
            _finish_step_success(
                "phylogeny",
                validator_errors=validate_phylogeny_outputs(CFG.output, CFG.genomes),
                reason="reused_existing_output",
                reference=snippy_reference_stem(CFG.output),
                snippy_tree=str(snippy_tree_path(CFG.output)),
                rooting=rooting,
                snippy_threads=thread_plan["phylogeny_snippy"],
                raxml_threads=thread_plan["phylogeny_raxml"],
            )
        else:
            if (CFG.output / "snippy").is_dir():
                logger.warning("Existing Snippy tree is missing or invalid; rerunning Snippy.")
            try:
                tools_used.add("snippy")
                tools_used.add("raxml-ng")
                phylogeny_result = run_phylogeny(
                    CFG,
                    snippy_threads=thread_plan["phylogeny_snippy"],
                    raxml_threads=thread_plan["phylogeny_raxml"],
                    progress=step_progress.genome_progress("phylogeny"),
                )
            except Exception as e:
                _step_end("phylogeny", "failed", error=str(e))
                raise
            else:
                _add_run_warnings(phylogeny_result.warnings)
                rooting = _write_trees_index()
                _finish_step_success(
                    "phylogeny",
                    validator_errors=validate_phylogeny_outputs(CFG.output, CFG.genomes),
                    reference=str(phylogeny_result.reference_genome),
                    snippy_tree=str(phylogeny_result.snippy_tree),
                    rooting=rooting,
                    snippy_threads=thread_plan["phylogeny_snippy"],
                    raxml_threads=thread_plan["phylogeny_raxml"],
                )

        if fastgear_failure is not None:
            # Gubbins, the core tree and Snippy ran on what was available; the run still
            # fails on the population structure step (fastGEAR or FastBAPS), before
            # anything else runs.
            raise fastgear_failure

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
                    status=step_progress.step_status("chewbbaca"),
                )
            except RunStoppedError:
                _check_stop()
                raise
            except Exception as e:
                # A missing cgMLST scheme keeps its code (the run card says what to do).
                code = {"code": e.code} if isinstance(e, CgmlstSchemaMissing) else {}
                _step_end("chewbbaca", "failed", error=str(e), **code)
                raise
            else:
                write_genome_set_stamp(CFG.output / "chewbbaca", CFG.genomes, step="chewbbaca")
                chewbbaca_left_out = chewbbaca_samples_left_out(CFG.output, CFG.genomes)
                if chewbbaca_left_out:
                    logger.warning(
                        "chewBBACA left %d sample(s) out of results_alleles.tsv (no CDSs "
                        "predicted): %s",
                        len(chewbbaca_left_out),
                        ", ".join(chewbbaca_left_out),
                    )
                _finish_step_success(
                    "chewbbaca",
                    validator_errors=validate_chewbbaca_outputs(CFG.output, CFG.genomes),
                    samples_left_out=chewbbaca_left_out,
                    output=str(chewbbaca_out) if chewbbaca_out else str(CFG.output / "chewbbaca" / "out"),
                    threads=thread_plan["chewbbaca"],
                )

        # DIAMOND marker DBs
        _check_stop()
        _step_start("diamond_markers")
        if not CFG.marker_sequences:
            query_db_files = []
            _step_end("diamond_markers", "skipped", reason="no marker database installed",
                      outputs_validated=True, validation_errors=[])
        elif _validate_dmnd_files(query_db_files_existing) and _marker_dmnds_current(CFG):
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
                logger.warning(
                    "Existing DIAMOND marker DBs are incomplete or older than their marker FASTA; "
                    "rebuilding them."
                )
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
        if not CFG.reference_locus_sequences:
            _step_end("diamond_reference_loci", "skipped", reason="no compound database installed",
                      outputs_validated=True, validation_errors=[])
        elif not validate_reference_hit_outputs(CFG.output, CFG.genomes, CFG.reference_locus_sequences):
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
                    progress=step_progress.genome_progress("diamond_reference_loci"),
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
                reuse_marker_hits=(
                    not validate_marker_hit_outputs(CFG.output, CFG.marker_sequences)
                    and _marker_hits_current(CFG)
                ),
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
            warnings=run_warnings,
            environment=environment,
            environment_history=environment_history,
            population=population,
        )
        _write_json(CFG.output / RUN_REPORT_FILENAME, run_report)
        _flush_run_state(status, error_override=error_text, report_written=True)
        _clear_run_lock(CFG.output)
        print("\nRUN REPORT\n" + json.dumps(run_report, indent=2, default=str))

    if error is not None:
        raise error


__all__ = ["run_workflow", "logger"]
