from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional, Sequence

from pipeline.ChewBBACAWorkflow import validate_chewbbaca_settings
from pipeline.DataInput import DataInput
from pipeline.DbCatalog import (
    MARKER_DB_CATALOG_PATH,
    COMPOUND_DB_CATALOG_PATH,
    DbCatalogIndex,
    load_db_catalog_index,
)
from pipeline.FastaUtils import collect_fastas, FASTA_EXTS
from pipeline.PopulationStructure import resume_phylo_overrides
from pipeline.ToolEnv import (
    _load_app_config,
    configured_assembly_qc_settings,
    configured_chewbbaca_settings,
    configured_fastgear_settings,
    configured_mobsuite_settings,
    effective_core_budget,
)

RUN_CONFIG_FILENAME = "run_config.json"  # pipeline.Runner's; not imported to keep this import light


def add_run_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """The options of a headless pipeline run. The databases default to the ones in
    bactolomew_config.json (first-time setup writes them)."""
    parser.add_argument("--genomes", metavar="DIR", help="Folder with assembled genome FASTA files "
                                          "(with --resume: default the run's recorded genomes)")
    parser.add_argument("--out", "--output", dest="output", metavar="DIR", help="Output folder for the run")
    parser.add_argument("--resume", action="store_true",
                        help="continue the run already in the output folder (its finished steps are kept)")
    parser.add_argument("--threads", "--cores", dest="cores", metavar="N", type=int, default=None,
                        help="Total cores the run plans with (default: the cores this process may use, "
                             "capped by resources.max_threads)")
    parser.add_argument("--queries", default=None,
                        help="FASTA file or folder with marker sequences (default: the marker DB catalog)")
    parser.add_argument("--compound", default=None,
                        help="FASTA file or folder with reference locus sequences (default: the compound DB catalog)")
    parser.add_argument("--bakta_db", default=None, help="Bakta database folder (default: bakta_db in bactolomew_config.json)")
    parser.add_argument("--mobsuite_db", default=None, help="MOB-suite database directory (default: bactolomew_config.json)")
    parser.add_argument("--checkm2_db", default=None, help="CheckM2 database DMND file (default: bactolomew_config.json)")
    parser.add_argument("--dmnd_threads", type=int, default=None, help="Threads for DIAMOND (default: --threads)")
    parser.add_argument("--bakta_threads", type=int, default=None, help="Threads for Bakta (default: --threads)")
    parser.add_argument("--window_length", type=int, default=500, help="Window length for clustering (default: 500)")
    parser.add_argument("--population-method", dest="population_method", default=None,
                        choices=("auto", "fastgear", "fastbaps"),
                        help="Population structure method (default: phylo.population_method, fastbaps; "
                             "auto: FastBAPS when the tools image has it, else fastGEAR; with --resume: "
                             "the run's method)")
    return parser


def read_arguments(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = add_run_arguments(argparse.ArgumentParser(description="Run the BACTolomew pipeline on a genome folder."))
    args = parser.parse_args(argv)
    if not args.output:
        parser.error("--out is required")
    if not args.genomes and not args.resume:
        parser.error("--genomes is required (or --resume to continue the run in --out)")
    return args


def _catalog_files(index) -> list[Path]:
    """Every sequence file a DB catalog lists (installed or not)."""
    out: set[Path] = set()
    for info in index.catalog.values():
        raw = info.get("db_path") or info.get("path") or info.get("db-path")
        if raw:
            out.add(Path(str(raw)).expanduser().resolve())
    return sorted(out)


def _fasta_inputs(arg: Optional[str], catalog_path: Path, label: str, logger) -> tuple[list[Path], dict[str, str]]:
    """The FASTAs given on the command line (file or folder), else the installed ones the
    catalog lists (the run logs and records the others as not installed); with their
    categories from the catalog. A missing catalog is an empty one (no database installed)."""
    try:
        index = load_db_catalog_index(catalog_path, logger=logger)
    except FileNotFoundError:
        index = DbCatalogIndex.from_catalog({})
    except Exception as e:
        logger.error(f"Failed to load {label} catalog: {e}")
        raise SystemExit(1)
    if arg:
        path = Path(arg).expanduser().resolve()
        if path.is_file() and path.suffix.lower() in FASTA_EXTS:
            files = [path]
        elif path.is_dir():
            files = collect_fastas(path)
        else:
            logger.error(f"{label} input is neither a valid FASTA file nor a folder: {path}")
            raise SystemExit(1)
        if not files:
            logger.error(f"No {label} FASTA files found ({arg}).")
            raise SystemExit(1)
    else:
        files = [p for p in _catalog_files(index) if p.is_file()]
    logger.info(f"Found {len(files)} {label} FASTA file(s).")
    return files, index.category_map(files, logger=logger, label=f"{label} database")


def _recorded_genomes(output_path: Path) -> list[Path]:
    try:
        run_cfg = json.loads((output_path / RUN_CONFIG_FILENAME).read_text(encoding="utf-8"))
    except Exception:
        return []
    inputs = run_cfg.get("inputs") if isinstance(run_cfg, dict) else None
    return [Path(p) for p in (inputs or {}).get("genomes") or []]


def workflow_cfg(args, logger) -> DataInput:
    output_path = Path(args.output).expanduser().resolve()
    has_run = (output_path / RUN_CONFIG_FILENAME).is_file()
    resume = bool(getattr(args, "resume", False))
    if resume and not has_run:
        logger.error(f"--resume: no run to continue in {output_path} (no {RUN_CONFIG_FILENAME}).")
        raise SystemExit(1)
    if has_run and not resume:
        logger.error(f"{output_path} already holds a run; pass --resume to continue it or choose another --output.")
        raise SystemExit(1)

    if args.genomes:
        genome_path = Path(args.genomes).expanduser().resolve()
        genome_files = collect_fastas(genome_path)
        if not genome_files:
            logger.error(f"No FASTA files found in genome folder: {genome_path}")
            raise SystemExit(1)
    else:
        genome_files = _recorded_genomes(output_path)
        if not genome_files:
            logger.error(f"{RUN_CONFIG_FILENAME} in {output_path} lists no genomes; pass --genomes.")
            raise SystemExit(1)
    logger.info(f"Found {len(genome_files)} genome FASTA files.")

    bakta_db = args.bakta_db or str(_load_app_config().get("bakta_db") or "").strip()
    if not bakta_db:
        logger.error("No Bakta database: run `bactolomew setup` or pass --bakta_db.")
        raise SystemExit(1)

    cpu_total = effective_core_budget()  # usable cores (affinity, cgroup, SLURM), capped by max_threads
    total_cores = args.cores or cpu_total
    if total_cores > cpu_total:
        logger.warning(f"--threads {total_cores} is more than the {cpu_total} cores this process may use.")
    dmnd_threads = args.dmnd_threads or total_cores
    bakta_threads = args.bakta_threads or total_cores
    window_length = args.window_length if args.window_length and args.window_length > 0 else 500

    marker_sequences, marker_cat_map = _fasta_inputs(args.queries, MARKER_DB_CATALOG_PATH, "marker", logger)
    ce_files, loc_cat_map = _fasta_inputs(args.compound, COMPOUND_DB_CATALOG_PATH, "compound", logger)

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

    try:
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
            phylo_config=(
                resume_phylo_overrides(output_path, args.population_method) if resume
                else ({"population_method": args.population_method} if args.population_method else {})
            ),
        )
    except ValueError as exc:
        logger.error(str(exc))
        raise SystemExit(1)


__all__ = ["add_run_arguments", "read_arguments", "workflow_cfg"]
