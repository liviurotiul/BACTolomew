from __future__ import annotations

import hashlib
import json
import logging
import pickle
import re
import shutil
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any, Callable, Mapping

from pipeline.ChewbbacaSchemaBuild import schema_in_use
from pipeline.RunSteps import StatusCallback, report_status
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
CHEWBBACA_SELF_SCORES = "self_scores"
MLST_RESULTS_JSON = "results.json"
# chewBBACA prints this for every BLAST match whose representative has no self-score;
# the match is then dropped, so novel alleles silently become missing data.
_MISSING_SELF_SCORE_RE = re.compile(r"Could not get the self-score")
# Written next to a run's chewBBACA output: which schema build produced the calls.
CHEWBBACA_SCHEMA_STAMP = "schema_stamp.json"
_SCHEMA_BUILD_INFO = Path("schema_build_report") / "build_info.json"
# Why the run's species has no usable cgMLST schema (CgmlstSchemaMissing.code).
CGMLST_SCHEMA_NOT_INSTALLED = "cgmlst_schema_not_installed"
CGMLST_NO_SCHEMA_FOR_SPECIES = "cgmlst_no_schema_for_species"
CGMLST_TERMS_NOT_ACCEPTED = "cgmlst_terms_not_accepted"
CGMLST_INSTALL_FAILED = "cgmlst_install_failed"


class CgmlstSchemaMissing(RuntimeError):
    """The run's dominant MLST species has no installed cgMLST schema.

    ``code`` is CGMLST_SCHEMA_NOT_INSTALLED when the catalog knows a schema that is not
    on disk (the chewBBACA step installs it from cgMLST.org), CGMLST_NO_SCHEMA_FOR_SPECIES
    when there is none to install (no catalog entry matches the MLST scheme, or cgMLST.org
    has no scheme for it), CGMLST_TERMS_NOT_ACCEPTED and CGMLST_INSTALL_FAILED when the
    install could not run. Every one but the first fails the chewBBACA step (the message
    says what to do before Resume).
    """

    def __init__(self, message: str, *, code: str, schema_id: str = "", catalog_path: Path | None = None,
                 species: str = "") -> None:
        super().__init__(message)
        self.code = code
        self.schema_id = schema_id
        self.catalog_path = catalog_path
        self.species = species


def short_species(name: str) -> str:
    """"Acinetobacter baumannii" -> "A. baumannii" ("Brucella spp." stays as it is)."""
    words = str(name or "").split()
    if len(words) >= 2 and words[0][:1].isupper() and words[1].lower().rstrip(".") not in ("spp", "sp"):
        return f"{words[0][0]}. {' '.join(words[1:])}"
    return str(name or "")


def _no_scheme_message(species: str, mlst_scheme: str = "") -> str:
    """'No cgMLST scheme for <species>: import one in Settings > Add schema ...'.

    Without a known species (an MLST scheme no catalog entry maps) the MLST scheme names it.
    """
    if not species:
        return (
            f"No cgMLST scheme for the MLST scheme {mlst_scheme}: import one in Settings > Add schema "
            f"with its species and {mlst_scheme} as its MLST scheme (or python -m pipeline.CgmlstInstall "
            f"--import <schema files> --species \"<genus species>\" --mlst {mlst_scheme}), then Resume the run."
        )
    return (
        f"No cgMLST scheme for {species}: import one in Settings > Add schema (or python -m "
        f"pipeline.CgmlstInstall --import <schema files> --species \"{species}\"), then Resume the run."
    )


def schema_build_fingerprint(schema_path: Path) -> str | None:
    """Hash of the schema's build_info.json (None without one).

    Only the build record is hashed, not the schema: inferred-mode AlleleCalls add
    alleles to the schema, which must not invalidate the calls they produced.
    """
    try:
        return hashlib.sha256((Path(schema_path) / _SCHEMA_BUILD_INFO).read_bytes()).hexdigest()
    except OSError:
        return None


def write_chewbbaca_schema_stamp(output_path: Path, schema_path: Path) -> Path:
    path = Path(output_path) / "chewbbaca" / CHEWBBACA_SCHEMA_STAMP
    payload = {
        "schema_dir": str(Path(schema_path).resolve()),
        "build_info_sha256": schema_build_fingerprint(schema_path),
    }
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    tmp.replace(path)
    return path


def chewbbaca_schema_stamp_issue(output_path: Path) -> str:
    """Why a run's chewBBACA calls may come from another schema build ('' when current)."""
    path = Path(output_path) / "chewbbaca" / CHEWBBACA_SCHEMA_STAMP
    try:
        stamp = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return f"chewBBACA calls have no schema stamp ({path}); they predate the current schema checks"
    schema_dir = str(stamp.get("schema_dir") or "") if isinstance(stamp, dict) else ""
    if not schema_dir:
        return f"chewBBACA schema stamp names no schema ({path})"
    if stamp.get("build_info_sha256") != schema_build_fingerprint(Path(schema_dir)):
        return f"chewBBACA schema {schema_dir} was rebuilt after these calls were made"
    return ""


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


def schema_representative_ids(schema_path: Path) -> set[str]:
    """Sequence IDs of every representative allele in the schema's short/ directory."""
    short_dir = Path(schema_path) / CHEWBBACA_SHORT_DIR
    ids: set[str] = set()
    if not short_dir.is_dir():
        return ids
    for fasta in short_dir.glob("*.fasta"):
        try:
            with fasta.open("r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if line.startswith(">"):
                        token = line[1:].strip().split()
                        if token:
                            ids.add(token[0])
        except OSError:
            continue
    return ids


class _PlainDataUnpickler(pickle.Unpickler):
    """Loads only plain data (dict, list, tuple, str, numbers): any class or function
    reference is refused, so a planted pickle cannot run code on the host."""

    def find_class(self, module: str, name: str):  # noqa: D401 - pickle hook
        raise pickle.UnpicklingError(f"refused global {module}.{name}")


def load_self_scores(schema_path: Path) -> dict | None:
    """The schema's cached BLASTp self-scores (short/self_scores), or None when absent/unreadable.

    chewBBACA's AlleleCall writes this cache (a dict of plain tuples) and it is never
    copied on import; it is read with :class:`_PlainDataUnpickler` all the same.
    """
    path = Path(schema_path) / CHEWBBACA_SHORT_DIR / CHEWBBACA_SELF_SCORES
    if not path.is_file():
        return None
    try:
        with path.open("rb") as fh:
            data = _PlainDataUnpickler(fh).load()
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


def self_score_issue(schema_path: Path) -> str:
    """Why the cached self-scores do not cover the representatives ("" when they do or are absent).

    chewBBACA computes short/self_scores once (the first AlleleCall on a schema)
    and reuses it forever. Every representative must have an entry: a BLAST match
    against a representative without one is discarded, so the allele call is lost.
    Extra entries are harmless (AlleleCall keeps the scores of candidates that
    were not added as representatives).
    """
    scores = load_self_scores(schema_path)
    if scores is None:
        return ""
    reps = schema_representative_ids(schema_path)
    missing = reps - set(scores)
    if not missing:
        return ""
    example = ", ".join(sorted(missing)[:3])
    return (
        f"{len(scores)} self-score entries for {len(reps)} representatives "
        f"({len(missing)} representatives have none, e.g. {example})"
    )


def _drop_self_scores(schema_path: Path) -> None:
    try:
        (Path(schema_path) / CHEWBBACA_SHORT_DIR / CHEWBBACA_SELF_SCORES).unlink()
    except FileNotFoundError:
        pass


def check_self_scores_before_call(schema_path: Path, logger: logging.Logger) -> None:
    """Delete a stale self-score cache so AlleleCall recomputes it for every representative."""
    issue = self_score_issue(schema_path)
    if issue:
        _drop_self_scores(schema_path)
        logger.warning(
            "chewBBACA schema self-score cache was incomplete (%s); deleted %s so AlleleCall "
            "recomputes it.",
            issue,
            Path(schema_path) / CHEWBBACA_SHORT_DIR / CHEWBBACA_SELF_SCORES,
        )


def check_self_scores_after_call(schema_path: Path, stdout_text: str) -> None:
    """Fail when AlleleCall ran without a self-score for every representative."""
    issue = self_score_issue(schema_path)
    if load_self_scores(schema_path) is None:
        issue = "AlleleCall did not write short/self_scores"
    missing_lines = len(_MISSING_SELF_SCORE_RE.findall(stdout_text or ""))
    if not issue and not missing_lines:
        return
    _drop_self_scores(schema_path)
    parts = [issue] if issue else []
    if missing_lines:
        parts.append(
            f"chewBBACA logged {missing_lines} 'Could not get the self-score' lines, so those "
            "BLAST matches were dropped"
        )
    raise RuntimeError(
        "chewBBACA allele calls are unreliable: " + "; ".join(parts) + ". The schema's "
        "self-score cache was deleted; rebuild the schema with PrepExternalSchema if this "
        f"repeats. Schema: {schema_path}"
    )


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
        raise CgmlstSchemaMissing(
            _no_scheme_message("", mlst_label),
            code=CGMLST_NO_SCHEMA_FOR_SPECIES,
            catalog_path=catalog_path,
            species=mlst_label,
        )

    schema_path = _resolve_catalog_schema_path(selected_entry, catalog_path, resolved)
    if not schema_path.exists():
        species = selected_entry.get("species") or [selected_entry["schema_id"]]
        raise CgmlstSchemaMissing(
            f"The cgMLST schema for {species[0]} ('{selected_entry['schema_id']}') is not installed "
            f"({schema_path}).",
            code=CGMLST_SCHEMA_NOT_INSTALLED,
            schema_id=str(selected_entry["schema_id"]),
            catalog_path=catalog_path,
            species=str(species[0]),
        )
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


def schema_orientation_issue(schema_path: Path) -> str:
    """Non-empty when representatives are stored on the reverse strand (an unadapted schema).

    chewBBACA expects every allele in coding orientation. A schema assembled from
    cgMLST.org/Ridom FASTAs without PrepExternalSchema keeps reference-strand
    alleles, and those loci can then never be called (LNF in every genome).
    """
    from pipeline.ChewbbacaSchemaBuild import cds_orientation, read_fasta_file, schema_translation_table

    short_dir = Path(schema_path) / CHEWBBACA_SHORT_DIR
    if not short_dir.is_dir():
        return ""
    table = schema_translation_table(schema_path)
    reverse = total = 0
    for fasta in short_dir.glob("*.fasta"):
        for _header, seq in read_fasta_file(fasta):
            total += 1
            if cds_orientation(seq, table) == "reverse":
                reverse += 1
    if not reverse:
        return ""
    return (
        f"{reverse} of {total} representative alleles are coding only on the reverse strand; "
        "rebuild the schema with chewBBACA PrepExternalSchema (Settings, \"Add schema\", "
        "or python -m pipeline.CgmlstInstall --install <species> --replace)"
    )


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
    if not schema_issues:
        orientation_issue = schema_orientation_issue(schema_path)
        if orientation_issue:
            schema_issues.append(orientation_issue)
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
    if ctx is not None and ctx.runtime in ("docker", "podman", "singularity"):
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


def cgmlst_terms_accepted(catalog_path: Path | None = None, config: Mapping[str, Any] | None = None) -> bool:
    """Whether the cgMLST.org terms were accepted for the configured cgMLST folder or for
    the folder of the run's catalog (``catalog_path``, else the one of ``config``), where
    the chewBBACA step installs a missing scheme."""
    from pipeline.CgmlstInstall import default_dest, terms_accepted

    if catalog_path is None:
        try:
            catalog_path = _default_catalog_path(resolve_chewbbaca_settings(config))
        except Exception:  # noqa: BLE001 - the configured folder decides
            catalog_path = None
    folders = [default_dest()] + ([Path(catalog_path).parent] if catalog_path else [])
    return any(terms_accepted(folder) for folder in folders)


def _install_missing_schema(
    missing: CgmlstSchemaMissing,
    logger: logging.Logger,
    *,
    threads: int,
    status: StatusCallback = None,
    should_cancel: Callable[[], bool] | None = None,
) -> Path:
    """Install the run's missing schema from cgMLST.org into the folder of the run's
    catalog, where it is registered; returns the installed schema.

    ``status`` follows the install ("Installing cgMLST scheme · A. baumannii", share);
    ``should_cancel`` (the run's Stop) ends it with CgmlstInstall.Cancelled.
    The cgMLST.org terms count as accepted when they were accepted for that folder or for
    the configured cgMLST folder (the GUI asks for them before such a run starts; the
    command line records them with --accept-terms). Raises CgmlstSchemaMissing, with what
    to do before Resume, when they were not, when cgMLST.org has no scheme for the
    species, or when the install fails.
    """
    from pipeline.CgmlstInstall import (
        Cancelled,
        _resolve_org,
        default_dest,
        fetch_org_index,
        install,
        load_catalog,
    )

    species = missing.species or missing.schema_id
    fail = dict(schema_id=missing.schema_id, catalog_path=missing.catalog_path, species=species)
    dest = Path(missing.catalog_path).parent if missing.catalog_path else default_dest()
    # The terms first: nothing is fetched from cgMLST.org before they are accepted.
    if not cgmlst_terms_accepted(missing.catalog_path or default_dest() / CHEWBBACA_SCHEME_CATALOG):
        raise CgmlstSchemaMissing(
            f"The cgMLST scheme for {species} is not installed, and downloading it from cgMLST.org "
            "needs their terms of use accepted once: run python -m pipeline.CgmlstInstall --accept-terms "
            "(or accept them in Settings > Download from cgMLST.org), then Resume the run.",
            code=CGMLST_TERMS_NOT_ACCEPTED, **fail,
        )
    try:
        index = fetch_org_index()
    except Exception as exc:  # noqa: BLE001 - reported with what to do
        raise CgmlstSchemaMissing(
            f"The cgMLST scheme for {species} is not installed and cgMLST.org could not be reached "
            f"to install it ({exc}). Resume the run to try again.",
            code=CGMLST_INSTALL_FAILED, **fail,
        ) from exc
    scheme, _why = _resolve_org(missing.schema_id, index, load_catalog(Path(missing.catalog_path))
                                if missing.catalog_path else [])
    if scheme is None or not scheme.bacterial:
        raise CgmlstSchemaMissing(_no_scheme_message(species), code=CGMLST_NO_SCHEMA_FOR_SPECIES, **fail)
    text = f"Installing cgMLST scheme · {short_species(species)}"
    logger.info("Installing the missing cgMLST schema '%s' from cgMLST.org into %s.", missing.schema_id, dest)
    report_status(status, text, 0.0)
    try:
        results = install([scheme.slug], dest=dest, cpu=max(1, int(threads)), accept=True, index=index,
                          progress=logger.info, should_cancel=should_cancel,
                          on_fraction=lambda _sid, share: report_status(status, text, share))
    except Exception as exc:  # noqa: BLE001 - reported with what to do
        results, error = [], str(exc) or type(exc).__name__
    else:
        error = ""
    finally:
        report_status(status, "", None)
    done = [r for r in results if r.status in ("installed", "skipped") and r.path and Path(r.path).is_dir()]
    if not done and should_cancel is not None and should_cancel():
        raise Cancelled("cancelled")
    if not done:
        reason = error or next((r.reason for r in results if r.reason), "nothing was installed")
        raise CgmlstSchemaMissing(
            f"Installing the cgMLST scheme for {species} from cgMLST.org failed: {reason}. Resume the run "
            "to try again, or install it in Settings > Download from cgMLST.org.",
            code=CGMLST_INSTALL_FAILED, **fail,
        )
    logger.info("Installed the cgMLST schema '%s' in %s.", missing.schema_id, done[0].path)
    return Path(done[0].path)


def run_may_install_schema(output_path: Path | None, config: Mapping[str, Any] | None = None) -> bool:
    """Whether a run's chewBBACA step may install a scheme from cgMLST.org (no network):
    True when its MLST result selects a schema that is not installed, or when there is no
    MLST result yet (the species is not known before the run); False without MLST-based
    auto-selection or when the selected schema is installed or not installable."""
    quiet = logging.getLogger("bactolomew.chewbbaca.preflight")
    quiet.disabled = True
    try:
        settings = resolve_chewbbaca_settings(config)
        if not settings.get("auto_select_from_mlst"):
            return False
        if output_path is None or not _load_mlst_rows(Path(output_path)):
            return True
        _select_runtime_schema(settings, Path(output_path), quiet)
    except CgmlstSchemaMissing as missing:
        return missing.code == CGMLST_SCHEMA_NOT_INSTALLED
    except Exception:  # noqa: BLE001 - the step reports it
        return False
    return False


def prepare_chewbbaca_schema(
    output_path: Path,
    logger: logging.Logger,
    *,
    threads: int = 1,
    config: Mapping[str, Any] | None = None,
    status: StatusCallback = None,
    should_cancel: Callable[[], bool] | None = None,
) -> dict[str, Any]:
    """The validated chewBBACA settings of a run: the schema is auto-selected from its MLST
    result and, when the run's species has a cgMLST.org scheme that is not installed,
    installed first (``status`` follows the install; ``should_cancel``, the run's Stop, ends
    it with CgmlstInstall.Cancelled). Raises CgmlstSchemaMissing (the step fails) when
    there is no schema to use."""
    settings = resolve_chewbbaca_settings(config)
    try:
        settings = _select_runtime_schema(settings, output_path, logger)
    except CgmlstSchemaMissing as missing:
        if missing.code != CGMLST_SCHEMA_NOT_INSTALLED:
            raise
        installed_path = _install_missing_schema(missing, logger, threads=threads, status=status,
                                                 should_cancel=should_cancel)
        settings = resolve_chewbbaca_settings(config)
        # The installed folder wins over a stale available_schemas path of the run's settings.
        settings["available_schemas"] = {**settings["available_schemas"],
                                         missing.schema_id: str(installed_path.resolve())}
        settings = _select_runtime_schema(settings, output_path, logger)
    return validate_chewbbaca_settings(settings)


def run_chewbbaca(
    genome_files: list[Path],
    output_path: Path,
    logger: logging.Logger,
    *,
    threads: int = 1,
    config: Mapping[str, Any] | None = None,
    status: StatusCallback = None,
    should_cancel: Callable[[], bool] | None = None,
) -> Path | None:
    """chewBBACA AlleleCall on the run's genomes, after preparing its schema (see
    prepare_chewbbaca_schema); raises CgmlstInstall.Cancelled when ``should_cancel`` ends
    the work before the call."""
    settings = prepare_chewbbaca_schema(output_path, logger, threads=threads, config=config, status=status,
                                        should_cancel=should_cancel)
    if should_cancel is not None and should_cancel():
        from pipeline.CgmlstInstall import Cancelled

        raise Cancelled("cancelled")

    genomes = [Path(p) for p in genome_files if Path(p).exists() and Path(p).is_file()]
    if not genomes:
        logger.error("No genome FASTA files provided to chewBBACA.")
        return None

    base_cmd = _resolve_chewbbaca_base_cmd()
    if base_cmd is None:
        logger.error("chewBBACA executable not available (configured env or PATH).")
        return None

    # CgmlstInstall refuses to replace this schema while the call runs.
    with schema_in_use(Path(settings["schema_dir"])):
        check_self_scores_before_call(Path(settings["schema_dir"]), logger)

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
                check_self_scores_after_call(Path(settings["schema_dir"]), res.stdout or "")
                write_chewbbaca_schema_stamp(output_path, Path(settings["schema_dir"]))
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
    "CGMLST_INSTALL_FAILED",
    "CGMLST_NO_SCHEMA_FOR_SPECIES",
    "CGMLST_SCHEMA_NOT_INSTALLED",
    "CGMLST_TERMS_NOT_ACCEPTED",
    "CHEWBBACA_OUT_FILES",
    "CHEWBBACA_GENE_LIST",
    "CHEWBBACA_SCHEMA_CONFIG",
    "CHEWBBACA_SCHEMA_STAMP",
    "CHEWBBACA_SHORT_DIR",
    "check_self_scores_after_call",
    "CgmlstSchemaMissing",
    "check_self_scores_before_call",
    "chewbbaca_schema_issues",
    "chewbbaca_schema_stamp_issue",
    "cgmlst_terms_accepted",
    "load_self_scores",
    "prepare_chewbbaca_schema",
    "schema_orientation_issue",
    "schema_build_fingerprint",
    "schema_representative_ids",
    "self_score_issue",
    "short_species",
    "ensure_chewbbaca_out",
    "resolve_chewbbaca_settings",
    "run_chewbbaca",
    "run_may_install_schema",
    "validate_chewbbaca_settings",
    "write_chewbbaca_schema_stamp",
]
