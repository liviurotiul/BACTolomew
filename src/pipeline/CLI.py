from __future__ import annotations

import argparse
import os
from pathlib import Path

def _load_gui():
    from plotting import GUI  # Local import to avoid PyQt dependency for CLI-only runs
    return GUI.run_gui_and_get_config
from pipeline.ChewBBACAWorkflow import validate_chewbbaca_settings
from pipeline.DataInput import DataInput
from pipeline.DbCatalog import (
    MARKER_DB_CATALOG_PATH,
    COMPOUND_DB_CATALOG_PATH,
    load_db_catalog_index,
)
from pipeline.FastaUtils import collect_fastas, FASTA_EXTS
from pipeline.ToolEnv import (
    configured_assembly_qc_settings,
    configured_chewbbaca_settings,
    configured_fastgear_settings,
    configured_mobsuite_settings,
    effective_core_budget,
)


def read_arguments():
    parser = argparse.ArgumentParser(description="Prepare genome and MarkerSequence FASTAs for DIAMOND.")
    parser.add_argument("--GUI", action="store_true")
    parser.add_argument("--genomes", help="Folder with assembled genome FASTA files")
    parser.add_argument("--queries", required=False, help="FASTA file or folder with MarkerSequence gene/protein sequences")
    parser.add_argument("--compound", type=str, required=False, help="FASTA folder with compound elements in fasta format")
    parser.add_argument("--output", help="Output folder for downstream results")
    parser.add_argument("--bakta_db", default=None, help="Path to Bakta database (optional, overrides env var)")
    parser.add_argument("--mobsuite_db", default=None, help="Path to MOB-suite database directory (optional, overrides bactolomew_config.json)")
    parser.add_argument("--checkm2_db", default=None, help="Path to CheckM2 database DMND file (optional, overrides bactolomew_config.json)")
    parser.add_argument("--dmnd_threads", type=int, default=None, help="Threads for DIAMOND (default: all cores)")
    parser.add_argument("--bakta_threads", type=int, default=None, help="Threads for Bakta   (default: all cores)")
    parser.add_argument("--cores", type=int, default=None, help="Total cores available for thread planning")
    parser.add_argument("--window_length", type=int, default=500, help="Window length for clustering (default: 10000)")

    args = parser.parse_args()

    if not args.GUI and not args.genomes:
        parser.error("--genomes is required when not using --GUI")

    return args


def workflow_cfg(args, logger) -> DataInput:
    if not args.GUI:
        genome_path = Path(args.genomes).expanduser().resolve()
        output_path = Path(args.output).expanduser().resolve()
        if not args.queries:
            logger.error("--queries is required (marker sequences).")
            raise SystemExit(1)
        if not args.compound:
            logger.error("--compound is required (reference locus sequences).")
            raise SystemExit(1)
        if not args.bakta_db:
            logger.error("--bakta_db is required.")
            raise SystemExit(1)

        marker_dir = Path(args.queries).expanduser().resolve()
        ce_path = Path(args.compound).expanduser().resolve()
        bakta_db = args.bakta_db
        cpu_total = effective_core_budget()  # min(physical cores, configured max_threads)
        dmnd_threads = args.dmnd_threads or cpu_total
        bakta_threads = args.bakta_threads or cpu_total
        total_cores = args.cores or cpu_total
        window_length = args.window_length if args.window_length and args.window_length > 0 else 500

        genome_files = collect_fastas(genome_path)
        if not genome_files:
            logger.error(f"No FASTA files found in genome folder: {genome_path}")
            raise SystemExit(1)
        logger.info(f"Found {len(genome_files)} genome FASTA files.")

        marker_sequences = []
        if marker_dir.is_file() and marker_dir.suffix.lower() in FASTA_EXTS:
            marker_sequences = [marker_dir]
        elif marker_dir.is_dir():
            marker_sequences = collect_fastas(marker_dir)
        else:
            logger.error(f"MarkerSequence input is neither a valid FASTA file nor a folder: {marker_dir}")
        if marker_sequences:
            logger.info(f"Found {len(marker_sequences)} MarkerSequence FASTA file(s).")
        else:
            logger.error(f"No MarkerSequence FASTA files found in: {marker_dir}")
            raise SystemExit(1)

        ce_files = []
        if ce_path.is_file() and ce_path.suffix.lower() in FASTA_EXTS:
            ce_files = [ce_path]
        elif ce_path.is_dir():
            ce_files = collect_fastas(ce_path)
        else:
            logger.error(f"Compound elements input is neither a valid FASTA file nor a folder: {ce_path}")
        if ce_files:
            logger.info(f"Found {len(ce_files)} compound element FASTA file(s).")
        else:
            logger.error(f"No compound element FASTA files found in: {ce_path}")
            raise SystemExit(1)

        try:
            marker_index = load_db_catalog_index(MARKER_DB_CATALOG_PATH, logger=logger)
        except Exception as e:
            logger.error(f"Failed to load marker DB catalog: {e}")
            raise SystemExit(1)

        try:
            locus_index = load_db_catalog_index(COMPOUND_DB_CATALOG_PATH, logger=logger)
        except Exception as e:
            logger.error(f"Failed to load compound DB catalog: {e}")
            raise SystemExit(1)

        marker_cat_map = marker_index.category_map(
            marker_sequences, logger=logger, label="marker database"
        )
        loc_cat_map = locus_index.category_map(
            ce_files, logger=logger, label="compound database"
        )
        fastgear_config = configured_fastgear_settings()
        if not fastgear_config.get("image"):
            logger.error("fastGEAR Docker image not configured in bactolomew_config.json (fastgear.image).")
            raise SystemExit(1)
        mobsuite_config = configured_mobsuite_settings(
            {"database_directory": args.mobsuite_db} if args.mobsuite_db else None
        )
        assembly_qc_config = configured_assembly_qc_settings(
            {"checkm2_database_path": args.checkm2_db} if args.checkm2_db else None
        )
        chewbbaca_config = configured_chewbbaca_settings()
        try:
            validate_chewbbaca_settings(chewbbaca_config)
        except Exception as exc:
            logger.error(str(exc))
            raise SystemExit(1)

        return DataInput(
            genomes=genome_files,
            marker_sequences=marker_sequences,
            reference_locus_sequences=ce_files,
            output=output_path,
            snippy_reference=None,
            bakta_db=bakta_db,
            dmnd_threads=dmnd_threads,
            bakta_threads=bakta_threads,
            window_length=window_length,
            diamond=True,
            bakta=True,
            total_cores=total_cores,
            marker_categories=marker_cat_map,
            reference_locus_categories=loc_cat_map,
            fastgear_config=fastgear_config,
            mobsuite_config=mobsuite_config,
            chewbbaca_config=chewbbaca_config,
            assembly_qc_config=assembly_qc_config,
        )
    try:
        run_gui_and_get_config = _load_gui()
    except Exception as e:
        logger.error(f"GUI mode requested but PyQt/GUI components could not be loaded: {e}")
        raise SystemExit(1)

    cfg = run_gui_and_get_config(args, show=True)
    if cfg is None:
        raise SystemExit(1)
    return cfg


__all__ = ["read_arguments", "workflow_cfg"]
