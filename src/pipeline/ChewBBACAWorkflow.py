from __future__ import annotations

import json
import logging
import shutil
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any, Mapping

from pipeline.ToolEnv import (
    configured_chewbbaca_settings,
    current_container_context,
    get_tool_env_prefix,
    project_root,
    tool_cmd,
)

CHEWBBACA_OUT_FILES: tuple[str, ...] = (
    "results_alleles.tsv",
    "results_statistics.tsv",
    "loci_summary_stats.tsv",
    "logging_info.txt",
)

_CHEWBBACA_EXECUTABLES: tuple[str, ...] = (
    "chewBBACA.py",
    "chewbbaca",
    "chewBBACA",
)

CHEWBBACA_SCHEMA_CONFIG = ".schema_config"
CHEWBBACA_GENE_LIST = ".genes_list"
CHEWBBACA_SHORT_DIR = "short"
CHEWBBACA_SCHEME_CATALOG = "scheme_species_catalog.json"
MLST_RESULTS_JSON = "results.json"


def resolve_chewbbaca_settings(raw_config: Mapping[str, Any] | None = None) -> dict[str, Any]:
    settings = dict(configured_chewbbaca_settings(raw_config))
    for key in ("schema_dir", "training_file", "genes_list", "loci_list", "catalog_path"):
        value = str(settings.get(key) or "").strip()
        settings[key] = str(Path(value).expanduser().resolve()) if value else ""
    available_schemas = settings.get("available_schemas")
    if isinstance(available_schemas, Mapping):
        resolved_schemas: dict[str, str] = {}
        for label, value in available_schemas.items():
            label_text = str(label).strip()
            path_text = str(value).strip()
            if not label_text or not path_text:
                continue
            resolved_schemas[label_text] = str(Path(path_text).expanduser().resolve())
        settings["available_schemas"] = resolved_schemas
    else:
        settings["available_schemas"] = {}
    extra_args = settings.get("extra_args")
    if not isinstance(extra_args, list):
        settings["extra_args"] = []
    return settings


def _is_internal_schema_genes_list(path_value: str) -> bool:
    text = str(path_value or "").strip()
    return bool(text) and Path(text).name == CHEWBBACA_GENE_LIST


def _default_catalog_path(settings: Mapping[str, Any]) -> Path:
    catalog_path = str(settings.get("catalog_path") or "").strip()
    if catalog_path:
        return Path(catalog_path).expanduser().resolve()

    schema_dir = str(settings.get("schema_dir") or "").strip()
    if schema_dir:
        return Path(schema_dir).expanduser().resolve().parent / CHEWBBACA_SCHEME_CATALOG

    return (project_root() / "DB" / "cgmlst" / CHEWBBACA_SCHEME_CATALOG).resolve()


def _load_mlst_rows(output_path: Path) -> list[dict[str, Any]]:
    path = output_path / "mlst" / MLST_RESULTS_JSON
    if not path.is_file():
        return []
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise RuntimeError(f"Failed to read MLST results from {path}: {exc}") from exc
    if not isinstance(raw, list):
        raise RuntimeError(f"Invalid MLST results format in {path}: expected a JSON list.")
    rows: list[dict[str, Any]] = []
    for item in raw:
        if isinstance(item, Mapping):
            rows.append(dict(item))
    return rows


def _load_scheme_catalog(catalog_path: Path) -> list[dict[str, Any]]:
    try:
        raw = json.loads(catalog_path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise RuntimeError(f"Failed to read cgMLST scheme catalog {catalog_path}: {exc}") from exc

    if isinstance(raw, Mapping):
        entries_raw = raw.get("schemes")
    else:
        entries_raw = raw
    if not isinstance(entries_raw, list):
        raise RuntimeError(f"Invalid cgMLST scheme catalog format in {catalog_path}: expected a list of schemes.")

    entries: list[dict[str, Any]] = []
    seen_schema_ids: set[str] = set()
    seen_aliases: dict[str, str] = {}
    for item in entries_raw:
        if not isinstance(item, Mapping):
            continue
        schema_id = str(item.get("schema_id") or "").strip()
        schema_dir = str(item.get("schema_dir") or schema_id).strip()
        if not schema_id or not schema_dir:
            continue
        species_raw = item.get("species")
        species = [str(value).strip() for value in species_raw] if isinstance(species_raw, list) else []
        species = [value for value in species if value]
        mlst_schemes_raw = item.get("mlst_schemes")
        mlst_schemes = (
            [str(value).strip() for value in mlst_schemes_raw]
            if isinstance(mlst_schemes_raw, list)
            else []
        )
        mlst_schemes = [value for value in mlst_schemes if value]
        if schema_id in seen_schema_ids:
            raise RuntimeError(f"Duplicate schema_id in cgMLST scheme catalog: {schema_id}")
        seen_schema_ids.add(schema_id)
        for alias in mlst_schemes:
            normalized = alias.lower()
            previous = seen_aliases.get(normalized)
            if previous is not None and previous != schema_id:
                raise RuntimeError(
                    f"Ambiguous MLST scheme alias '{alias}' in cgMLST catalog: "
                    f"{previous} and {schema_id}"
                )
            seen_aliases[normalized] = schema_id
        entries.append(
            {
                "schema_id": schema_id,
                "schema_dir": schema_dir,
                "species": species,
                "mlst_schemes": mlst_schemes,
            }
        )
    return entries


def _resolve_catalog_schema_path(entry: Mapping[str, Any], catalog_path: Path, settings: Mapping[str, Any]) -> Path:
    schema_id = str(entry.get("schema_id") or "").strip()
    available_schemas = settings.get("available_schemas")
    if isinstance(available_schemas, Mapping):
        configured = str(available_schemas.get(schema_id) or "").strip()
        if configured:
            return Path(configured).expanduser().resolve()
    schema_dir = str(entry.get("schema_dir") or "").strip()
    path = Path(schema_dir).expanduser()
    if path.is_absolute():
        return path.resolve()
    return (catalog_path.parent / path).resolve()


def _select_runtime_schema(
    settings: Mapping[str, Any],
    output_path: Path,
    logger: logging.Logger,
) -> dict[str, Any]:
    resolved = dict(settings)
    if not bool(resolved.get("auto_select_from_mlst")):
        return resolved

    mlst_rows = _load_mlst_rows(output_path)
    if not mlst_rows:
        logger.info("chewBBACA schema auto-selection skipped: no MLST results found under %s.", output_path / "mlst")
        return resolved

    catalog_path = _default_catalog_path(resolved)
    if not catalog_path.is_file():
        logger.info("chewBBACA schema auto-selection skipped: cgMLST catalog not found at %s.", catalog_path)
        return resolved

    entries = _load_scheme_catalog(catalog_path)
    alias_to_entry = {
        alias.lower(): entry
        for entry in entries
        for alias in entry.get("mlst_schemes", [])
    }

    counts: Counter[str] = Counter()
    labels: dict[str, str] = {}
    raw_schemes: dict[str, Counter[str]] = {}
    mapped_entries: dict[str, dict[str, Any]] = {}

    for row in mlst_rows:
        scheme = str(row.get("scheme") or "").strip()
        if not scheme:
            continue
        entry = alias_to_entry.get(scheme.lower())
        if entry is None:
            key = f"mlst:{scheme.lower()}"
            labels[key] = scheme
        else:
            key = str(entry["schema_id"])
            mapped_entries[key] = entry
            species = entry.get("species") or []
            labels[key] = str(species[0]) if species else key
        counts[key] += 1
        raw = raw_schemes.setdefault(key, Counter())
        raw[scheme] += 1

    if not counts:
        raise RuntimeError("MLST completed but produced no scheme assignments that could be used for cgMLST selection.")

    top_count = max(counts.values())
    top_keys = [key for key, value in counts.items() if value == top_count]
    if len(top_keys) > 1:
        detail = ", ".join(
            f"{labels.get(key, key)} ({counts[key]})"
            for key in sorted(top_keys, key=lambda item: labels.get(item, item).lower())
        )
        raise RuntimeError(
            "Multiple species are tied in MLST results, so cgMLST cannot auto-select one schema: "
            f"{detail}"
        )

    selected_key = top_keys[0]
    selected_entry = mapped_entries.get(selected_key)
    if selected_entry is None:
        mlst_label = ", ".join(sorted(raw_schemes.get(selected_key, {}).keys()))
        raise RuntimeError(
            "The most frequent MLST species/scheme has no bundled cgMLST schema mapping: "
            f"{mlst_label}"
        )

    schema_path = _resolve_catalog_schema_path(selected_entry, catalog_path, resolved)
    resolved["schema_dir"] = str(schema_path)
    genes_list = schema_path / CHEWBBACA_GENE_LIST
    resolved["genes_list"] = str(genes_list.resolve()) if genes_list.is_file() else ""
    resolved["catalog_path"] = str(catalog_path.resolve())

    species = selected_entry.get("species") or []
    logger.info(
        "Selected cgMLST schema '%s' for dominant MLST species '%s' (%d/%d genomes; MLST schemes: %s).",
        selected_entry["schema_id"],
        species[0] if species else selected_entry["schema_id"],
        top_count,
        sum(counts.values()),
        ", ".join(sorted(raw_schemes.get(selected_key, {}).keys())),
    )
    return resolved


def chewbbaca_schema_issues(schema_path: Path) -> list[str]:
    path = Path(schema_path).expanduser()
    if not path.exists() or not path.is_dir():
        return [f"chewBBACA schema directory not found: {path}"]

    entries = list(path.iterdir())
    locus_fastas = sorted(
        p.name for p in entries
        if p.is_file() and p.name.lower().endswith(".fasta")
    )
    fastq_files = sorted(
        p.name for p in entries
        if p.is_file() and (
            p.name.lower().endswith(".fastq")
            or p.name.lower().endswith(".fastq.gz")
        )
    )

    issues: list[str] = []
    if not locus_fastas:
        if fastq_files:
            issues.append(
                "directory contains FASTQ files, but chewBBACA schemas require one locus FASTA "
                "(.fasta) per gene at the schema root"
            )
        else:
            issues.append("no locus FASTA files (.fasta) were found at the schema root")

    short_dir = path / CHEWBBACA_SHORT_DIR
    if not short_dir.exists() or not short_dir.is_dir():
        issues.append(
            "missing required short/ directory with representative allele FASTA files"
        )
    else:
        short_fastas = [
            p.name for p in short_dir.iterdir()
            if p.is_file() and p.name.lower().endswith(".fasta")
        ]
        if not short_fastas:
            issues.append("short/ exists but contains no representative allele FASTA files")

    config_path = path / CHEWBBACA_SCHEMA_CONFIG
    if not config_path.is_file():
        if locus_fastas:
            issues.append(
                "missing .schema_config; this looks like a legacy/external schema and must be "
                "adapted with chewBBACA PrepExternalSchema before AlleleCall"
            )
        else:
            issues.append("missing required .schema_config file")

    return issues


def validate_chewbbaca_settings(raw_config: Mapping[str, Any] | None = None) -> dict[str, Any]:
    settings = resolve_chewbbaca_settings(raw_config)
    schema_dir = str(settings.get("schema_dir") or "").strip()
    if not schema_dir:
        raise ValueError(
            "chewBBACA schema directory not configured. "
            "Set bactolomew_config.json chewbbaca.schema_dir before running the pipeline."
        )

    schema_path = Path(schema_dir)
    schema_issues = chewbbaca_schema_issues(schema_path)
    if schema_issues:
        details = "\n".join(f"- {issue}" for issue in schema_issues)
        raise ValueError(f"Invalid chewBBACA schema directory: {schema_path}\n{details}")

    training_file = str(settings.get("training_file") or "").strip()
    if training_file:
        training_path = Path(training_file)
        if not training_path.exists() or not training_path.is_file():
            raise ValueError(f"chewBBACA training file not found: {training_path}")

    genes_list = str(settings.get("genes_list") or "").strip()
    if genes_list:
        genes_path = Path(genes_list)
        if not genes_path.exists() or not genes_path.is_file():
            raise ValueError(f"chewBBACA genes list not found: {genes_path}")
        if not _is_internal_schema_genes_list(genes_list):
            raise ValueError(
                "chewBBACA genes_list must point to the schema's internal .genes_list file; "
                "use chewbbaca.loci_list or chewbbaca.gl for an AlleleCall subset list."
            )

    loci_list = str(settings.get("loci_list") or "").strip()
    if loci_list:
        loci_path = Path(loci_list)
        if not loci_path.exists() or not loci_path.is_file():
            raise ValueError(f"chewBBACA loci list not found: {loci_path}")

    return settings


def _candidate_cmds(*args: str) -> list[list[str]]:
    prefix = get_tool_env_prefix("chewbbaca")
    out: list[list[str]] = []
    for executable in _CHEWBBACA_EXECUTABLES:
        if prefix:
            out.append(["conda", "run", "-p", prefix, executable, *args])
        out.append([executable, *args])
    return out


def _resolve_chewbbaca_base_cmd() -> list[str] | None:
    # Inside a container runtime, chewBBACA runs in the tools image's `chewbbaca`
    # env exactly like every other tool (the run's dirs are bind-mounted at
    # identical paths), so no host env or PATH probe is needed.
    ctx = current_container_context()
    if ctx is not None and ctx.runtime in ("docker", "singularity"):
        return tool_cmd("chewbbaca", "chewBBACA.py")
    # Native / dev runtime: probe the configured host env or PATH.
    for probe_args in (("--help",), ("--version",)):
        for cmd in _candidate_cmds(*probe_args):
            try:
                res = subprocess.run(cmd, capture_output=True, text=True, check=False)
            except FileNotFoundError:
                continue
            if res.returncode == 0:
                return cmd[:-len(probe_args)]
    return None


def _find_chewbbaca_output_dir(search_root: Path, *, preferred_dir: Path | None = None) -> Path | None:
    candidates = {p.parent for p in search_root.rglob("results_alleles.tsv") if p.is_file()}
    if not candidates:
        return None

    def _rank(path: Path) -> tuple[int, int, str]:
        pref = 1
        if preferred_dir is not None and path == preferred_dir:
            pref = 0
        try:
            depth = len(path.relative_to(search_root).parts)
        except ValueError:
            depth = 999
        return (pref, depth, str(path))

    return sorted(candidates, key=_rank)[0]


def ensure_chewbbaca_out(output_path: Path, logger: logging.Logger) -> Path | None:
    chew_root = output_path / "chewbbaca"
    out_dir = chew_root / "out"
    if not chew_root.exists():
        return None

    source_dir = _find_chewbbaca_output_dir(chew_root, preferred_dir=out_dir)
    if source_dir is None:
        return None

    if source_dir != out_dir:
        if out_dir.exists():
            shutil.rmtree(out_dir, ignore_errors=True)
        out_dir.mkdir(parents=True, exist_ok=True)
        for item in source_dir.iterdir():
            target = out_dir / item.name
            if target.exists():
                if target.is_dir():
                    shutil.rmtree(target, ignore_errors=True)
                else:
                    target.unlink()
            shutil.move(str(item), str(target))
        logger.info("chewBBACA outputs consolidated: %s -> %s", source_dir, out_dir)
    else:
        out_dir.mkdir(parents=True, exist_ok=True)

    if not (out_dir / "results_alleles.tsv").is_file():
        return None
    return out_dir


def run_chewbbaca(
    genome_files: list[Path],
    output_path: Path,
    logger: logging.Logger,
    *,
    threads: int = 1,
    config: Mapping[str, Any] | None = None,
) -> Path | None:
    settings = resolve_chewbbaca_settings(config)
    settings = _select_runtime_schema(settings, output_path, logger)
    settings = validate_chewbbaca_settings(settings)

    genomes = [Path(p) for p in genome_files if Path(p).exists() and Path(p).is_file()]
    if not genomes:
        logger.error("No genome FASTA files provided to chewBBACA.")
        return None

    base_cmd = _resolve_chewbbaca_base_cmd()
    if base_cmd is None:
        logger.error("chewBBACA executable not available (configured env or PATH).")
        return None

    chew_root = output_path / "chewbbaca"
    if chew_root.exists() and chew_root.is_dir():
        for item in chew_root.iterdir():
            if item.is_dir():
                shutil.rmtree(item)
            else:
                item.unlink()
    chew_root.mkdir(parents=True, exist_ok=True)
    stdout_log = chew_root / "command.stdout.log"
    stderr_log = chew_root / "command.stderr.log"

    with tempfile.TemporaryDirectory(dir=chew_root) as tmpdir:
        tmp_root = Path(tmpdir)
        input_list = tmp_root / "genomes.txt"
        input_list.write_text(
            "\n".join(str(path.resolve()) for path in genomes) + "\n",
            encoding="utf-8",
        )

        tmp_out = tmp_root / "out"
        cmd = [
            *base_cmd,
            "AlleleCall",
            "-i",
            str(input_list),
            "-g",
            str(settings["schema_dir"]),
            "-o",
            str(tmp_out),
            "--cpu",
            str(max(1, int(threads))),
        ]
        training_file = str(settings.get("training_file") or "").strip()
        if training_file:
            cmd.extend(["--ptf", training_file])
        loci_list = str(settings.get("loci_list") or "").strip()
        if loci_list:
            cmd.extend(["--gl", loci_list])
        extra_args = settings.get("extra_args") or []
        if isinstance(extra_args, list):
            cmd.extend(str(arg) for arg in extra_args if str(arg).strip())

        logger.info("Running chewBBACA AlleleCall on %d genomes -> %s", len(genomes), tmp_out)
        logger.info("CMD: %s", " ".join(cmd))
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, check=False)
            stdout_log.write_text(res.stdout or "", encoding="utf-8")
            stderr_log.write_text(res.stderr or "", encoding="utf-8")
            if res.returncode != 0:
                stderr_text = (res.stderr or "").strip()
                stdout_text = (res.stdout or "").strip()
                detail = stderr_text or stdout_text or f"exit code {res.returncode}"
                raise RuntimeError(
                    "chewBBACA AlleleCall failed. "
                    f"{detail}\nstdout log: {stdout_log}\nstderr log: {stderr_log}"
                )
            logger.info("chewBBACA completed.")
        except Exception as exc:
            logger.error("chewBBACA failed: %s", exc)
            raise

        out_dir = ensure_chewbbaca_out(output_path, logger)
        if out_dir is None:
            detail = ""
            try:
                stderr_text = stderr_log.read_text(encoding="utf-8").strip()
            except Exception:
                stderr_text = ""
            if stderr_text:
                detail = f" Last stderr: {stderr_text.splitlines()[-1]}"
            raise RuntimeError(
                "chewBBACA completed but results_alleles.tsv was not found in the output."
                f"{detail}\nstdout log: {stdout_log}\nstderr log: {stderr_log}"
            )
        return out_dir


__all__ = [
    "CHEWBBACA_OUT_FILES",
    "CHEWBBACA_GENE_LIST",
    "CHEWBBACA_SCHEMA_CONFIG",
    "CHEWBBACA_SHORT_DIR",
    "chewbbaca_schema_issues",
    "ensure_chewbbaca_out",
    "resolve_chewbbaca_settings",
    "run_chewbbaca",
    "validate_chewbbaca_settings",
]
