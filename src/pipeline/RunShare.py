"""Shared runs: one archive with what the explorer needs to show a run, its import
into another install's runs folder, and the checks that say what may not match.

A share package is a zip (ZIP64) with one top folder, the run folder: the run's
small result files (``SHARE_PATTERNS``), its data tables in ``share/tables/`` and a
manifest (``share/manifest.json``: the BACTolomew version that exported it, the run,
the tools record, the format versions, every file with its size and sha256). It
never holds the files only needed to build the tables (Bakta JSON, DIAMOND hits),
nor the working data of the tools, nor the exporting machine's paths (its text files
name paths inside the run relative to it, others by their file name only). An
imported run gets ``share/imported.json``.

The compatibility check only reports: nothing it finds blocks opening a run. No Qt
and no ``annotation_explorer_polars`` imports.
"""

from __future__ import annotations

import fnmatch
import hashlib
import json
import os
import re
import shutil
import stat
import uuid
import zipfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path, PurePosixPath
from typing import Any, Callable, Iterable, Mapping, Optional, Sequence

SHARE_FORMAT_VERSION = 1
SHARE_DIR = "share"
MANIFEST = "share/manifest.json"
IMPORTED = "share/imported.json"
TABLES_DIR = "share/tables"
METADATA_CSV = "share/metadata.csv"
SCHEME_SPECIES = "share/scheme_species.json"
TABLES_META = "tables_meta.json"
ARCHIVE_SUFFIX = ".bactolomew.zip"
# Imports are unpacked here first, then moved into SHARED_RUNS_DIR.
STAGING_PREFIX = ".import-"
SHARED_RUNS_DIR = "shared"

# The run files the explorer reads (with the tables supplied): about 28 MB for 153
# genomes. Globs relative to the run folder; ``**`` spans folders.
SHARE_PATTERNS = (
    "run_config.json",
    "run_report.json",
    "run_state.json",
    "trees/**",
    "gubbins/*.json",
    "gubbins/core_gene.final_tree.tre",
    "fastgear/sample_clades.tsv",
    "fastgear/output/lineage_information.txt",
    "fastbaps/sample_clades.tsv",
    "fastbaps/*.json",
    "panaroo/out/gene_presence_absence.csv",
    "panaroo/out/summary_statistics.txt",
    "chewbbaca/out/results_alleles.tsv",
    "chewbbaca/out/results_statistics.tsv",
    "chewbbaca/out/loci_summary_stats.tsv",
    "chewbbaca/out/logging_info.txt",
    "chewbbaca/schema_stamp.json",
    "chewbbaca/genome_set.json",
    "mlst/results.csv",
    "mlst/results.json",
    "mobsuite/plasmids.tsv",
    "mobsuite/contigs.tsv",
    "mobsuite/summary.json",
    "assembly_qc/summary.tsv",
    "assembly_qc/summary.json",
    "assembly_qc/quast/report.tsv",
    "assembly_qc/checkm2/quality_report.tsv",
    "amrfinder/*.tsv",
    "amrfinder/summary.json",
    "evidence/*.tsv",
    "evidence/summary.json",
    "core_alignment/ordering.json",
    "core_alignment/gene_map.tsv",
    "snippy/reference.txt",
    "snippy/core/snippy_raxml.tree",
)
# The genome assemblies, optional (about 4 MB per genome).
ASSEMBLY_PATTERN = "bakta/*/*.fna"
# Never copied: the run lock, the stop request and the tools' home folder.
EXCLUDE = frozenset({"run.lock.json", "run.stop", ".bactolomew_home"})
RUN_FILES = ("run_config.json", "run_report.json", "run_state.json")

# Already compressed: stored, not deflated.
_STORED_SUFFIXES = (".parquet", ".fna", ".gz", ".zip", ".zst")
# Text files whose paths of the exporting machine are rewritten (trees hold none).
_TEXT_SUFFIXES = (".json", ".tsv", ".csv", ".txt", ".log", ".iqtree")
# An absolute path of at least two parts (or one under ~), not inside a word, a
# relative path, a ratio or a URL.
_PATH_PART = r"[^/\s\"'<>|;,()\[\]{}=]"
_ABSOLUTE_PATH = re.compile(rf"(?<![\w.~/)\]-])(?:~|/{_PATH_PART}+)(?:/{_PATH_PART}*)+")
# A home folder anywhere, also inside forms the path pattern skips ("-o/home/u/x",
# "file:///home/u/x"): replaced by "~" so the user name never leaves the machine.
_HOME_FOLDER = re.compile(rf"/(?:home|Users)/{_PATH_PART}+|/root(?=/|\b)")
# The tools image's own folders, the same on every machine: kept.
_IMAGE_PREFIXES = ("/opt/", "/usr/", "/bin/", "/etc/", "/work/")
_MANIFEST_MAX_BYTES = 16 * 1024 * 1024
_JSON_MAX_BYTES = 64 * 1024 * 1024
_CHUNK = 1024 * 1024


class ShareError(RuntimeError):
    """A share package that cannot be written or read (the message is for the user)."""


class ShareCancelled(ShareError):
    def __init__(self) -> None:
        super().__init__("Cancelled.")


@dataclass(frozen=True)
class Issue:
    """One compatibility finding: ``level`` "warning" (shown on the run card and
    before opening) or "info" (Run Details only). ``short`` fits the card line."""

    level: str
    code: str
    text: str
    short: str = ""


@dataclass(frozen=True)
class ImportResult:
    path: Path
    run_id: str
    damaged: tuple[str, ...] = ()
    missing: tuple[str, ...] = ()
    renamed_from: Optional[str] = None


ProgressFn = Optional[Callable[[int, int, str], None]]
CancelFn = Optional[Callable[[], bool]]


# ----------------------------------------------------------------- helpers


def parse_version(text: Any) -> Optional[tuple[int, ...]]:
    """'0.2.0' -> (0, 2, 0); a suffix after the numbers is ignored; None when unreadable."""
    match = re.match(r"\s*v?(\d+(?:\.\d+)*)", str(text or ""))
    if not match:
        return None
    return tuple(int(part) for part in match.group(1).split("."))


def _read_json(path: Path, limit: int = _JSON_MAX_BYTES) -> Any:
    try:
        if path.stat().st_size > limit:
            return None
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _json_dict(path: Path) -> dict:
    data = _read_json(path)
    return data if isinstance(data, dict) else {}


def _now_iso() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(_CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _check(cancelled: CancelFn) -> None:
    if cancelled is not None and cancelled():
        raise ShareCancelled()


def _safe_name(text: str, fallback: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text or "")).strip("._")
    return slug or fallback


def rebase_recorded_path(path: Any, recorded_root: Any, current_root: Any) -> Path:
    """A path a run recorded under its folder on the machine that made it
    (``recorded_root``), moved to the folder the run is in now; any other path,
    or a run that did not move, is returned as given."""
    p = Path(str(path)).expanduser()
    if not recorded_root or current_root is None:
        return p
    old, new = Path(str(recorded_root)).expanduser(), Path(current_root)
    if not p.is_absolute():
        return new / p
    try:
        rel = p.relative_to(old)
    except ValueError:
        return p
    return new / rel


def recorded_root(run_dir: Path) -> Optional[str]:
    """The folder the run was made in, as its run_config.json recorded it."""
    value = _json_dict(Path(run_dir) / "run_config.json").get("output_dir")
    return str(value) if value else None


def is_shared_run(run_dir: Path) -> bool:
    """A run unpacked from a share package (its tables in ``share/tables``)."""
    return (Path(run_dir) / MANIFEST).is_file()


def is_imported_run(run_dir: Path) -> bool:
    """A run brought in with Import Run... (a share package or a copied run folder)."""
    return (Path(run_dir) / IMPORTED).is_file()


def read_manifest(source: Path) -> dict:
    """The manifest of a share package (a ``.zip``) or of a run folder; {} when
    there is none. Data only: its size is capped and nothing in it is run."""
    source = Path(source)
    if source.is_dir():
        data = _read_json(source / MANIFEST, _MANIFEST_MAX_BYTES)
        return data if isinstance(data, dict) else {}
    try:
        with zipfile.ZipFile(source) as archive:
            info = _manifest_member(archive)
            if info is None or info.file_size > _MANIFEST_MAX_BYTES:
                return {}
            data = json.loads(archive.read(info).decode("utf-8"))
    except (OSError, ValueError, zipfile.BadZipFile, KeyError):
        return {}
    return data if isinstance(data, dict) else {}


def _manifest_member(archive: zipfile.ZipFile) -> Optional[zipfile.ZipInfo]:
    for info in archive.infolist():
        parts = PurePosixPath(info.filename).parts
        if len(parts) == 3 and "/".join(parts[1:]) == MANIFEST:
            return info
    return None


def _matches(rel: str, pattern: str) -> bool:
    if pattern.endswith("/**"):
        return rel.startswith(pattern[:-2])
    return fnmatch.fnmatchcase(rel, pattern) and rel.count("/") == pattern.count("/")


def share_files(run_dir: Path, *, include_assemblies: bool = False) -> list[tuple[str, int]]:
    """(relative path, size) of the run files a share package holds, tables aside."""
    run_dir = Path(run_dir)
    patterns = SHARE_PATTERNS + ((ASSEMBLY_PATTERN,) if include_assemblies else ())
    out: list[tuple[str, int]] = []
    seen: set[str] = set()
    for pattern in patterns:
        if "*" not in pattern:
            candidates: Iterable[Path] = [run_dir / pattern]
        else:
            # The folder before the first wildcard ("bakta" for "bakta/*/*.fna").
            root = run_dir / pattern.split("*", 1)[0].rpartition("/")[0]
            candidates = root.rglob("*") if root.is_dir() else []
        for path in candidates:
            try:
                if path.is_symlink() or not path.is_file():
                    continue
                rel = path.relative_to(run_dir).as_posix()
            except (OSError, ValueError):
                continue
            if rel in seen or path.name in EXCLUDE or rel.startswith(SHARE_DIR + "/"):
                continue
            if "*" in pattern and not _matches(rel, pattern):
                continue
            seen.add(rel)
            out.append((rel, path.stat().st_size))
    return sorted(out)


def _table_files(tables_dir: Optional[Path]) -> list[Path]:
    if tables_dir is None or not Path(tables_dir).is_dir():
        return []
    return sorted(p for p in Path(tables_dir).iterdir() if p.is_file() and p.suffix == ".parquet")


def estimate_share_size(run_dir: Path, tables_dir: Optional[Path] = None) -> dict[str, int]:
    """Bytes of the package's parts: ``run_files``, ``tables`` and ``assemblies``
    (zip overhead aside; parquet is stored, the text files shrink when deflated)."""
    run_files = sum(size for _rel, size in share_files(run_dir))
    with_fna = sum(size for _rel, size in share_files(run_dir, include_assemblies=True))
    tables = sum(p.stat().st_size for p in _table_files(tables_dir))
    return {"run_files": run_files, "tables": tables, "assemblies": with_fna - run_files}


def _run_status(run_dir: Path) -> str:
    for name in ("run_state.json", "run_report.json"):
        status = str(_json_dict(Path(run_dir) / name).get("status") or "").strip().lower()
        if status:
            return status
    return "unknown"


def _environment_copy(environment: Any) -> Optional[dict]:
    """The run's tools record without local paths (database folders, image files)."""
    if not isinstance(environment, Mapping):
        return None

    def _clean(value: Any) -> Any:
        if isinstance(value, Mapping):
            return {k: _clean(v) for k, v in value.items() if k not in {"path", "sif_path"}}
        if isinstance(value, list):
            return [_clean(v) for v in value]
        return value

    return _clean(environment)


def _schemas(run_dir: Path, tables_schema: Optional[int]) -> dict[str, Any]:
    return {
        "tables": tables_schema,
        "trees": _json_dict(run_dir / "trees" / "trees.json").get("schema_version"),
        "gubbins_summary": _json_dict(run_dir / "gubbins" / "summary.json").get("schema_version"),
        "core_alignment_ordering": _json_dict(run_dir / "core_alignment" / "ordering.json").get("schema_version"),
    }


def _path_roots(run_dir: Path) -> list[str]:
    """The run folder as this machine and the machine that made it name it, longest first."""
    roots = {str(run_dir), str(run_dir.resolve())}
    made_in = recorded_root(run_dir)
    if made_in:
        roots.add(made_in)
    return sorted((root.rstrip("/") for root in roots if root.rstrip("/")), key=len, reverse=True)


def _strip_path(path: str, roots: Sequence[str]) -> str:
    """A path of the exporting machine as the package keeps it: relative to the run
    folder when inside it ("." for the folder), the image's folders as they are,
    any other path (genomes, databases, schemes, home folders) only its name."""
    for root in roots:
        if path.rstrip("/") == root:
            return "."
        if path.startswith(root + "/"):
            return path[len(root) + 1:]
    if path.startswith(_IMAGE_PREFIXES):
        return path
    if _HOME_FOLDER.fullmatch(path.rstrip("/")):
        return "~"
    return PurePosixPath(path.rstrip("/")).name


def _strip_text(text: str, roots: Sequence[str]) -> str:
    text = _ABSOLUTE_PATH.sub(lambda match: _strip_path(match.group(0), roots), text)
    return _HOME_FOLDER.sub("~", text)


def _strip_json(value: Any, roots: Sequence[str]) -> Any:
    if isinstance(value, str):
        return _strip_text(value, roots)
    if isinstance(value, Mapping):
        return {_strip_text(str(k), roots): _strip_json(v, roots) for k, v in value.items()}
    if isinstance(value, list):
        return [_strip_json(v, roots) for v in value]
    return value


def _shared_copy(path: Path, roots: Sequence[str]) -> Optional[bytes]:
    """A text file's content without the exporting machine's paths, or None when
    nothing changes (the file is then copied as it is)."""
    if not path.name.endswith(_TEXT_SUFFIXES) or path.stat().st_size > _JSON_MAX_BYTES:
        return None
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return None
    if not _ABSOLUTE_PATH.search(text):
        return None
    if path.suffix == ".json":
        try:
            data = json.loads(text)
        except ValueError:
            data = None
        if data is not None:
            return json.dumps(_strip_json(data, roots), indent=2).encode("utf-8")
    return _strip_text(text, roots).encode("utf-8")


def _scheme_catalog(scheme_species: Mapping[str, str]) -> dict:
    """``share/scheme_species.json``: MLST scheme -> species in the cgMLST scheme
    catalog's own format, so RunMetadata reads it like the catalog."""
    return {
        "schemes": [
            {"species": [species], "mlst_schemes": [scheme]}
            for scheme, species in sorted(scheme_species.items())
            if scheme and species
        ]
    }


# ----------------------------------------------------------------- export


def export_run(
    run_dir: Path,
    tables_dir: Optional[Path],
    dest: Path,
    *,
    include_assemblies: bool = False,
    metadata_csv: Optional[Path] = None,
    scheme_species: Optional[Mapping[str, str]] = None,
    app_version: str = "",
    progress: ProgressFn = None,
    cancelled: CancelFn = None,
) -> Path:
    """Write the run's share package to ``dest`` (through ``dest.part``, renamed
    when complete). ``tables_dir`` holds the run's parquet tables and their
    ``tables_meta.json``. Only completed, idle runs are shared."""
    run_dir = Path(run_dir)
    dest = Path(dest)
    if (run_dir / "run.lock.json").exists():
        raise ShareError("This run is still running.")
    if _run_status(run_dir) != "success":
        raise ShareError("Only completed runs can be shared.")
    tables = _table_files(tables_dir)
    if not tables:
        raise ShareError("The run has no data tables yet. Open it once, then export it.")
    config = _json_dict(run_dir / "run_config.json")
    report = _json_dict(run_dir / "run_report.json")
    run_id = str(config.get("run_id") or run_dir.name)
    top = _safe_name(run_id, "run")
    tables_meta = _json_dict(Path(tables_dir) / TABLES_META)
    tables_schema = tables_meta.get("schema_version")

    # (archive path inside the run folder, source file or bytes). The run files keep
    # no paths of this machine: the importing side finds everything in the run folder.
    roots = _path_roots(run_dir)
    entries: list[tuple[str, Any]] = []
    for rel, _size in share_files(run_dir, include_assemblies=include_assemblies):
        _check(cancelled)
        stripped = _shared_copy(run_dir / rel, roots)
        entries.append((rel, stripped if stripped is not None else run_dir / rel))
    entries += [(f"{TABLES_DIR}/{p.name}", p) for p in tables]
    # The cache's meta names this machine's folder: only the format is kept.
    entries.append((f"{TABLES_DIR}/{TABLES_META}", json.dumps({"schema_version": tables_schema}, indent=2).encode()))
    if metadata_csv is not None:
        entries.append((METADATA_CSV, Path(metadata_csv)))
    if scheme_species:
        entries.append((SCHEME_SPECIES, json.dumps(_scheme_catalog(scheme_species), indent=2).encode()))

    total = sum(len(src) if isinstance(src, bytes) else src.stat().st_size for _rel, src in entries)
    files: list[dict[str, Any]] = []
    part = dest.with_name(dest.name + ".part")
    dest.parent.mkdir(parents=True, exist_ok=True)
    done = 0
    try:
        with zipfile.ZipFile(part, "w", allowZip64=True) as archive:
            for rel, src in entries:
                _check(cancelled)
                info = zipfile.ZipInfo(f"{top}/{rel}", date_time=datetime.now().timetuple()[:6])
                info.compress_type = (
                    zipfile.ZIP_STORED if rel.endswith(_STORED_SUFFIXES) else zipfile.ZIP_DEFLATED
                )
                info.external_attr = 0o644 << 16
                digest = hashlib.sha256()
                size = 0
                with archive.open(info, "w", force_zip64=True) as out:
                    if isinstance(src, bytes):
                        out.write(src)
                        digest.update(src)
                        size = len(src)
                    else:
                        with Path(src).open("rb") as handle:
                            for chunk in iter(lambda: handle.read(_CHUNK), b""):
                                _check(cancelled)
                                out.write(chunk)
                                digest.update(chunk)
                                size += len(chunk)
                                done += len(chunk)
                                if progress is not None:
                                    progress(done, total, rel)
                if isinstance(src, bytes):
                    done += size
                files.append({"path": rel, "size": size, "sha256": digest.hexdigest()})
            manifest = {
                "format_version": SHARE_FORMAT_VERSION,
                "app_version": app_version,
                "exported_at": _now_iso(),
                "run_id": run_id,
                "sample_set_id": str(config.get("sample_set_id") or ""),
                "run_created_at": config.get("created_at"),
                "status": _run_status(run_dir),
                "sample_count": len((config.get("inputs") or {}).get("genomes") or []),
                "made_with": config.get("bactolomew_version") or report.get("bactolomew_version"),
                "environment": _environment_copy(config.get("environment") or report.get("environment")),
                "tool_versions": (report.get("tools") or {}).get("versions") or {},
                "schemas": _schemas(run_dir, tables_schema),
                "includes_assemblies": bool(include_assemblies),
                "metadata": "csv" if metadata_csv is not None else "run",
                "files": files,
            }
            archive.writestr(f"{top}/{MANIFEST}", json.dumps(manifest, indent=2))
        _check(cancelled)
        os.replace(part, dest)
    except BaseException:
        part.unlink(missing_ok=True)
        raise
    if progress is not None:
        progress(total, total, "")
    return dest


# ----------------------------------------------------------------- import


def _run_ids(runs_root: Path) -> set[str]:
    out: set[str] = set()
    if not Path(runs_root).is_dir():
        return out
    for cfg in Path(runs_root).rglob("run_config.json"):
        if any(part.startswith(STAGING_PREFIX) for part in cfg.parts):
            continue
        run_id = str(_json_dict(cfg).get("run_id") or "")
        if run_id:
            out.add(run_id)
    return out


def _destination(runs_root: Path, sample_set_id: str, run_id: str) -> Path:
    base = Path(runs_root) / SHARED_RUNS_DIR / _safe_name(sample_set_id, "shared")
    name = _safe_name(run_id, "run")
    dest = base / name
    n = 2
    while dest.exists():
        dest = base / f"{name}_{n}"
        n += 1
    return dest


def _unique_run_id(runs_root: Path, run_id: str) -> Optional[str]:
    """A new run id when ``run_id`` is taken in this runs folder (its cache and
    remembered metadata are kept per run id), else None."""
    taken = _run_ids(runs_root)
    if run_id not in taken:
        return None
    n = 1
    while True:
        candidate = f"{run_id}_shared" + (f"{n}" if n > 1 else "")
        if candidate not in taken:
            return candidate
        n += 1


def run_id_taken(runs_root: Path, run_id: str) -> bool:
    """A run with this run id is already in the runs folder (an import becomes a copy)."""
    return bool(run_id) and str(run_id) in _run_ids(runs_root)


def import_run_id(source: Path) -> str:
    """The run id a package (``.zip``) or a run folder would be imported with."""
    source = Path(source)
    if source.is_dir():
        return str(_json_dict(source / "run_config.json").get("run_id") or "")
    return str(read_manifest(source).get("run_id") or "")


def _rename_run(run_dir: Path, new_id: str) -> list[str]:
    rewritten: list[str] = []
    for name in RUN_FILES:
        path = run_dir / name
        data = _read_json(path)
        if isinstance(data, dict) and "run_id" in data:
            data["run_id"] = new_id
            path.write_text(json.dumps(data, indent=2), encoding="utf-8")
            rewritten.append(name)
    return rewritten


def _finish_import(
    staging: Path,
    runs_root: Path,
    *,
    source: Path,
    manifest: Mapping[str, Any],
    damaged: list[str],
    missing: list[str],
) -> ImportResult:
    config = _json_dict(staging / "run_config.json")
    if not config:
        raise ShareError("This is not a BACTolomew run: run_config.json is missing.")
    run_id = str(config.get("run_id") or manifest.get("run_id") or "run")
    new_id = _unique_run_id(runs_root, run_id)
    rewritten = _rename_run(staging, new_id) if new_id else []
    record = {
        "imported_at": _now_iso(),
        "source": source.name,
        "kind": "package" if manifest else "folder",
        "format_version": manifest.get("format_version"),
        "damaged": damaged,
        "missing": missing,
    }
    if new_id:
        record.update(run_id_from=run_id, rewritten=rewritten)
    (staging / SHARE_DIR).mkdir(exist_ok=True)
    (staging / IMPORTED).write_text(json.dumps(record, indent=2), encoding="utf-8")
    dest = _destination(runs_root, str(config.get("sample_set_id") or manifest.get("sample_set_id") or ""), new_id or run_id)
    dest.parent.mkdir(parents=True, exist_ok=True)
    os.replace(staging, dest)
    return ImportResult(
        path=dest, run_id=new_id or run_id, damaged=tuple(damaged), missing=tuple(missing), renamed_from=run_id if new_id else None
    )


def _member_rel(info: zipfile.ZipInfo, top: Optional[str]) -> tuple[str, str]:
    """(top folder, path inside it) of a zip member; ShareError for an unsafe name."""
    name = info.filename
    mode = info.external_attr >> 16
    if stat.S_ISLNK(mode):
        raise ShareError(f"The package holds a link, which is not allowed: {name}")
    if name.startswith(("/", "\\")) or "\\" in name or re.match(r"^[A-Za-z]:", name):
        raise ShareError(f"The package holds an absolute path: {name}")
    parts = PurePosixPath(name).parts
    if any(part in ("..", ".") for part in parts) or len(parts) < 2:
        raise ShareError(f"The package holds an unsafe path: {name}")
    if top is not None and parts[0] != top:
        raise ShareError("The package must hold one run folder.")
    return parts[0], "/".join(parts[1:])


def import_archive(
    archive_path: Path,
    runs_root: Path,
    *,
    max_bytes: int = 256 * 1024**3,
    progress: ProgressFn = None,
    cancelled: CancelFn = None,
) -> ImportResult:
    """Unpack a share package into ``runs_root/shared/<sample set>/<run>``: names are
    checked (no absolute paths, ``..`` or links), sizes capped, every file checked
    against its manifest sha256 (a damaged or missing file is reported, not fatal).
    A run id already in this runs folder gets a ``_shared`` suffix."""
    archive_path = Path(archive_path)
    runs_root = Path(runs_root)
    staging = runs_root / f"{STAGING_PREFIX}{uuid.uuid4().hex[:12]}"
    try:
        archive = zipfile.ZipFile(archive_path)
    except (OSError, zipfile.BadZipFile) as exc:
        raise ShareError(f"Not a share package: {exc}") from exc
    with archive:
        members = [info for info in archive.infolist() if not info.is_dir()]
        top: Optional[str] = None
        rels: list[tuple[zipfile.ZipInfo, str]] = []
        for info in members:
            top, rel = _member_rel(info, top)
            rels.append((info, rel))
        manifest = read_manifest(archive_path)
        if not manifest:
            raise ShareError("Not a share package: share/manifest.json is missing.")
        total = sum(info.file_size for info, _rel in rels)
        if total > max_bytes:
            raise ShareError("The package is larger than allowed.")
        runs_root.mkdir(parents=True, exist_ok=True)
        free = shutil.disk_usage(runs_root).free
        if total > free:
            raise ShareError(f"Not enough disk space: {total / 1e9:.1f} GB needed, {free / 1e9:.1f} GB free.")
        expected = {
            str(entry.get("path")): entry
            for entry in manifest.get("files") or []
            if isinstance(entry, Mapping) and entry.get("path")
        }
        damaged: list[str] = []
        done = 0
        try:
            staging.mkdir(parents=True)
            for info, rel in rels:
                _check(cancelled)
                target = staging / rel
                target.parent.mkdir(parents=True, exist_ok=True)
                digest = hashlib.sha256()
                with archive.open(info) as src, target.open("wb") as out:
                    for chunk in iter(lambda: src.read(_CHUNK), b""):
                        _check(cancelled)
                        out.write(chunk)
                        digest.update(chunk)
                        done += len(chunk)
                        if progress is not None:
                            progress(done, total, rel)
                entry = expected.get(rel)
                if entry is not None and entry.get("sha256") and entry.get("sha256") != digest.hexdigest():
                    damaged.append(rel)
            missing = sorted(rel for rel in expected if not (staging / rel).is_file())
            return _finish_import(
                staging, runs_root, source=archive_path, manifest=manifest, damaged=sorted(damaged), missing=missing
            )
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise


def import_folder(
    src: Path,
    runs_root: Path,
    *,
    progress: ProgressFn = None,
    cancelled: CancelFn = None,
) -> ImportResult:
    """Copy a run folder (a full run, or an unpacked share package) into
    ``runs_root/shared/...``, without its lock, stop request or tools home."""
    src = Path(src).resolve()
    runs_root = Path(runs_root)
    if not (src / "run_config.json").is_file():
        raise ShareError("This folder is not a BACTolomew run: run_config.json is missing.")
    try:
        src.relative_to(runs_root.resolve())
    except ValueError:
        pass
    else:
        raise ShareError("This run is already in the run list.")
    manifest = read_manifest(src)
    expected = {
        str(entry.get("path")): entry
        for entry in manifest.get("files") or []
        if isinstance(entry, Mapping) and entry.get("path")
    }
    files: list[tuple[Path, str]] = []
    for root, dirs, names in os.walk(src):
        dirs[:] = [d for d in dirs if d not in EXCLUDE and not os.path.islink(os.path.join(root, d))]
        for name in names:
            path = Path(root) / name
            if name in EXCLUDE or path.is_symlink():
                continue
            rel = path.relative_to(src).as_posix()
            if rel == IMPORTED:
                continue
            files.append((path, rel))
    total = sum(p.stat().st_size for p, _rel in files)
    runs_root.mkdir(parents=True, exist_ok=True)
    if total > shutil.disk_usage(runs_root).free:
        raise ShareError("Not enough disk space to copy this run.")
    staging = runs_root / f"{STAGING_PREFIX}{uuid.uuid4().hex[:12]}"
    damaged: list[str] = []
    done = 0
    try:
        staging.mkdir(parents=True)
        for path, rel in files:
            _check(cancelled)
            target = staging / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)
            done += path.stat().st_size
            if progress is not None:
                progress(done, total, rel)
            entry = expected.get(rel)
            if entry is not None and entry.get("sha256") and entry.get("sha256") != _sha256(target):
                damaged.append(rel)
        missing = sorted(rel for rel in expected if not (staging / rel).is_file())
        return _finish_import(staging, runs_root, source=src, manifest=manifest, damaged=damaged, missing=missing)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def remove_imported_run(run_dir: Path) -> None:
    """Delete an imported run's folder (only imported runs: never a run made here)."""
    run_dir = Path(run_dir)
    if not is_imported_run(run_dir):
        raise ShareError("Only imported runs can be removed here.")
    shutil.rmtree(run_dir)


# ----------------------------------------------------------- compatibility


def _flatten_tool_versions(record: Any) -> dict[str, str]:
    """tool -> version from a tool_versions.json record (``envs.<env>.tools``) or a
    flat mapping."""
    out: dict[str, str] = {}
    if not isinstance(record, Mapping):
        return out
    envs = record.get("envs")
    if isinstance(envs, Mapping):
        for env in envs.values():
            tools = env.get("tools") if isinstance(env, Mapping) else None
            for name, info in (tools or {}).items() if isinstance(tools, Mapping) else ():
                version = info.get("version") if isinstance(info, Mapping) else info
                if version:
                    out[str(name)] = str(version)
        return out
    for name, version in record.items():
        if isinstance(version, (str, int, float)) and str(version).strip():
            out[str(name)] = str(version)
    return out


def _facts(source: Any) -> dict[str, Any]:
    """What the checks look at, from a run folder or from a manifest (an archive
    not unpacked yet)."""
    if isinstance(source, Mapping):
        manifest = dict(source)
        schemas = manifest.get("schemas") or {}
        return {
            "shared": True,
            "format_version": manifest.get("format_version"),
            "made_with": manifest.get("made_with"),
            "exported_with": manifest.get("app_version"),
            "status": manifest.get("status"),
            "environment": manifest.get("environment"),
            "tables_schema": schemas.get("tables"),
            "trees_schema": schemas.get("trees"),
            "gubbins_schema": schemas.get("gubbins_summary"),
            "has_tables": any(str(f.get("path", "")).endswith(".parquet") for f in manifest.get("files") or [] if isinstance(f, Mapping)),
            "missing": [],
            "damaged": [],
        }
    run_dir = Path(source)
    manifest = read_manifest(run_dir) if is_shared_run(run_dir) else {}
    imported = _json_dict(run_dir / IMPORTED)
    config = _json_dict(run_dir / "run_config.json")
    report = _json_dict(run_dir / "run_report.json")
    shared = bool(manifest)
    missing: list[str] = list(imported.get("missing") or [])
    if shared:
        listed = [str(f.get("path")) for f in manifest.get("files") or [] if isinstance(f, Mapping) and f.get("path")]
        missing = sorted(set(missing) | {rel for rel in listed if not (run_dir / rel).is_file()})
    tables_dir = run_dir / TABLES_DIR
    return {
        "shared": shared,
        "format_version": manifest.get("format_version") if shared else None,
        "made_with": config.get("bactolomew_version") or report.get("bactolomew_version") or manifest.get("made_with"),
        "exported_with": manifest.get("app_version") if shared else None,
        "status": _run_status(run_dir),
        "environment": config.get("environment") or report.get("environment"),
        "tables_schema": _json_dict(tables_dir / TABLES_META).get("schema_version") if shared else None,
        "trees_schema": _json_dict(run_dir / "trees" / "trees.json").get("schema_version"),
        "gubbins_schema": _json_dict(run_dir / "gubbins" / "summary.json").get("schema_version"),
        "has_tables": bool(_table_files(tables_dir)) if shared else None,
        "missing": missing,
        "damaged": list(imported.get("damaged") or []),
    }


def _names(paths: Sequence[str], limit: int = 3) -> str:
    shown = ", ".join(paths[:limit])
    return shown + (f" and {len(paths) - limit} more" if len(paths) > limit else "")


def _files_word(paths: Sequence[str]) -> str:
    return "file" if len(paths) == 1 else "files"


def _count(paths: Sequence[str]) -> str:
    """'1 file' / '3 files'."""
    return f"{len(paths)} {_files_word(paths)}"


def _newer(a: Any, b: Any) -> bool:
    try:
        return a is not None and b is not None and int(a) > int(b)
    except (TypeError, ValueError):
        return False


def check_compatibility(
    source: Any,
    *,
    app_version: str,
    table_schema: Optional[int] = None,
    trees_schema: Optional[int] = None,
    gubbins_schema: Optional[int] = None,
    image_ref: Optional[str] = None,
    image_tool_versions: Optional[Mapping[str, Any]] = None,
) -> list[Issue]:
    """What may not match between a run (a folder, or a package's manifest) and this
    install. Warnings first; none of them blocks opening the run. The versions of
    this install come in as arguments (``current_versions()`` has the pipeline's)."""
    facts = _facts(source)
    issues: list[Issue] = []

    def warn(code: str, text: str, short: str) -> None:
        issues.append(Issue("warning", code, text, short))

    def info(code: str, text: str) -> None:
        issues.append(Issue("info", code, text, text))

    if _newer(facts["format_version"], SHARE_FORMAT_VERSION):
        warn(
            "format_newer",
            f"The package format (v{facts['format_version']}) is newer than this app reads (v{SHARE_FORMAT_VERSION}).",
            "newer package format",
        )
    app = parse_version(app_version)
    # The exporter's version is not the run's: it is checked on its own below.
    made = facts["made_with"]
    made_v = parse_version(made)
    if made_v is None:
        info("version_unknown", "The BACTolomew version that made this run is not recorded (made before 0.2.0).")
    elif app is not None and made_v[:2] > app[:2]:
        warn("app_newer", f"Made with BACTolomew {made}; this app is {app_version}.", f"made with BACTolomew {made}")
    elif app is not None and made_v < app:
        info("app_older", f"Made with an older BACTolomew ({made}); this app is {app_version}.")
    exported = parse_version(facts["exported_with"])
    if exported is not None and app is not None and exported[:2] > app[:2] and exported != made_v:
        warn(
            "exporter_newer",
            f"Shared with BACTolomew {facts['exported_with']}; this app is {app_version}.",
            f"shared with BACTolomew {facts['exported_with']}",
        )
    if facts["shared"] and table_schema is not None:
        if facts["has_tables"] is False:
            warn("no_tables", "The package has no data tables: the run cannot be shown.", "no data tables")
        elif facts["tables_schema"] is not None and facts["tables_schema"] != table_schema:
            warn(
                "table_schema",
                f"The data tables are format v{facts['tables_schema']}; this app reads v{table_schema} and cannot"
                " rebuild them without the full run.",
                "other data table format",
            )
    if _newer(facts["trees_schema"], trees_schema):
        warn(
            "trees_schema",
            f"The tree index is format v{facts['trees_schema']}; this app reads up to v{trees_schema}.",
            "newer tree index",
        )
    if _newer(facts["gubbins_schema"], gubbins_schema):
        warn(
            "gubbins_schema",
            f"The Gubbins summary is format v{facts['gubbins_schema']}; this app reads up to v{gubbins_schema}.",
            "newer Gubbins summary",
        )
    if facts["damaged"]:
        damaged = list(facts["damaged"])
        warn(
            "damaged",
            f"{_count(damaged)} damaged (checksum differs): {_names(damaged)}.",
            f"{len(damaged)} damaged {_files_word(damaged)}",
        )
    if facts["missing"]:
        missing = list(facts["missing"])
        warn("missing", f"{_count(missing)} missing: {_names(missing)}.", f"{len(missing)} missing {_files_word(missing)}")
    status = str(facts["status"] or "unknown")
    if status != "success":
        warn("status", f"The run did not complete (status: {status}).", "run not completed")

    environment = facts["environment"] if isinstance(facts["environment"], Mapping) else None
    if environment is None:
        info("no_environment", "No record of the tools image (made before 0.2.0).")
    else:
        image = str(environment.get("image") or "")
        names = {image, str(environment.get("configured_image") or "")}
        if image and image_ref and image_ref not in names:
            info("image_differs", f"Tools image {image}; this install uses {image_ref}.")
        theirs = _flatten_tool_versions(environment.get("tool_versions"))
        ours = _flatten_tool_versions(image_tool_versions)
        differ = sorted(name for name in theirs if name in ours and theirs[name] != ours[name])
        if differ:
            detail = ", ".join(f"{n} {theirs[n]} (here {ours[n]})" for n in differ[:4])
            more = f" and {len(differ) - 4} more" if len(differ) > 4 else ""
            info("tool_versions", f"{len(differ)} tool versions differ from this install's image: {detail}{more}.")
    return sorted(issues, key=lambda issue: issue.level != "warning")


def warning_summary(issues: Sequence[Issue]) -> str:
    """The run card line: '2 warnings: made with BACTolomew 0.3.0'; '' without warnings."""
    warnings = [issue for issue in issues if issue.level == "warning"]
    if not warnings:
        return ""
    noun = "warning" if len(warnings) == 1 else "warnings"
    return f"{len(warnings)} {noun}: {warnings[0].short or warnings[0].text}"


def current_versions() -> dict[str, Any]:
    """This install's versions the pipeline knows: the app, the tree index and
    Gubbins summary formats, the tools image its config names and that image's tool
    versions (the table format is the explorer's, ``cache.TABLE_CACHE_SCHEMA_VERSION``).
    A viewer install or a native runtime uses no image: both are None then."""
    from pipeline import __version__
    from pipeline.ToolEnv import (
        APP_MODE_VIEWER,
        _load_app_config,
        _resolve_runtime,
        configured_app_mode,
        configured_image_ref,
        project_root,
    )

    config = _load_app_config()
    uses_image = configured_app_mode(config) != APP_MODE_VIEWER and _resolve_runtime(config) != "native"
    out: dict[str, Any] = {"app_version": __version__, "image_ref": configured_image_ref(config) if uses_image else None}
    try:
        from pipeline.TreesIndex import SCHEMA_VERSION as trees_schema

        out["trees_schema"] = trees_schema
    except Exception:
        pass
    try:
        from pipeline.GubbinsWorkflow import SUMMARY_SCHEMA_VERSION

        out["gubbins_schema"] = SUMMARY_SCHEMA_VERSION
    except Exception:
        pass
    lock = _read_json(project_root() / "docker" / "locks" / "tool_versions.json") if uses_image else None
    out["image_tool_versions"] = lock if isinstance(lock, Mapping) else None
    return out


# ----------------------------------------------------------------- run info


def _fmt_bytes(size: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024 or unit == "GB":
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} GB"


def format_size(size: int) -> str:
    return _fmt_bytes(int(size))


def run_info(run_dir: Path) -> dict[str, list[tuple[str, str]]]:
    """Run Details: (label, value) rows per section ("Run", "Made with",
    "Analysis", "Shared"; a section without rows is left out)."""
    run_dir = Path(run_dir)
    config = _json_dict(run_dir / "run_config.json")
    report = _json_dict(run_dir / "run_report.json")
    inputs = config.get("inputs") if isinstance(config.get("inputs"), Mapping) else {}
    manifest = read_manifest(run_dir) if is_shared_run(run_dir) else {}
    imported = _json_dict(run_dir / IMPORTED)
    environment = config.get("environment") or report.get("environment") or manifest.get("environment") or {}
    sections: dict[str, list[tuple[str, str]]] = {}

    def add(section: str, label: str, value: Any) -> None:
        text = str(value).strip() if value is not None else ""
        if text:
            sections.setdefault(section, []).append((label, text))

    genomes = inputs.get("genomes") if isinstance(inputs.get("genomes"), list) else []
    add("Run", "Run", config.get("run_id") or run_dir.name)
    add("Run", "Sample set", config.get("sample_set_id"))
    add("Run", "Genomes", len(genomes) if genomes else manifest.get("sample_count"))
    status = _run_status(run_dir)
    add("Run", "Status", "completed" if status == "success" else status)
    add("Run", "Started", report.get("start_time") or config.get("created_at"))
    add("Run", "Finished", report.get("end_time"))
    add("Run", "Folder", str(run_dir))

    add("Made with", "BACTolomew", config.get("bactolomew_version") or report.get("bactolomew_version") or manifest.get("made_with") or "not recorded")
    if isinstance(environment, Mapping) and environment:
        add("Made with", "Runtime", environment.get("runtime"))
        add("Made with", "Tools image", environment.get("image"))
        image_id = str(environment.get("image_id") or "")
        add("Made with", "Image id", image_id[:19] if image_id.startswith("sha256:") else image_id)
        add("Made with", "Image built", environment.get("image_created"))
        dbs = environment.get("databases") if isinstance(environment.get("databases"), Mapping) else {}
        bakta = dbs.get("bakta") if isinstance(dbs.get("bakta"), Mapping) else {}
        if bakta:
            add("Made with", "Bakta database", " ".join(
                bit for bit in (f"v{bakta.get('major')}.{bakta.get('minor')}" if bakta.get("major") is not None else "",
                                str(bakta.get("type") or ""), str(bakta.get("date") or "")) if bit
            ))
            add("Made with", "AMRFinderPlus database", bakta.get("amrfinderplus_db"))
        mob = dbs.get("mobsuite") if isinstance(dbs.get("mobsuite"), Mapping) else {}
        add("Made with", "MOB-suite database", mob.get("download_date"))
        tools = _flatten_tool_versions(environment.get("tool_versions"))
        if tools:
            add("Made with", "Tool versions", ", ".join(f"{k} {v}" for k, v in sorted(tools.items())))
    else:
        add("Made with", "Tools image", "not recorded (made before 0.2.0)")
    versions = (report.get("tools") or {}).get("versions") if isinstance(report.get("tools"), Mapping) else None
    if isinstance(versions, Mapping) and versions:
        add("Made with", "Reported tools", "; ".join(f"{k}: {v}" for k, v in sorted(versions.items())))

    population = inputs.get("population_structure")
    method = population.get("method") if isinstance(population, Mapping) else None
    if not method:
        method = "fastGEAR" if (run_dir / "fastgear" / "sample_clades.tsv").is_file() else (
            "FastBAPS" if (run_dir / "fastbaps" / "sample_clades.tsv").is_file() else None
        )
    add("Analysis", "Population structure", {"fastgear": "fastGEAR", "fastbaps": "FastBAPS", "none": "none"}.get(str(method), method))
    gubbins = _json_dict(run_dir / "gubbins" / "summary.json")
    if gubbins:
        mode = str(gubbins.get("resolved_mode") or gubbins.get("mode") or "")
        ran = "ran" if gubbins.get("ran") else "not run"
        add("Analysis", "Gubbins", f"{mode.replace('_', ' ')}, {ran}" if mode else ran)
    trees = _json_dict(run_dir / "trees" / "trees.json")
    tree_list = trees.get("trees") if isinstance(trees.get("trees"), list) else []
    labels = [str(t.get("label") or t.get("id")) for t in tree_list if isinstance(t, Mapping)]
    add("Analysis", "Trees", "; ".join(labels))
    chew = inputs.get("chewbbaca") if isinstance(inputs.get("chewbbaca"), Mapping) else {}
    schema_dir = str(chew.get("schema_dir") or "")
    add("Analysis", "cgMLST scheme", Path(schema_dir).name.replace("_", " ") if schema_dir else None)
    markers = [Path(str(p)).stem for p in inputs.get("marker_sequences") or []]
    loci = [Path(str(p)).stem for p in inputs.get("reference_locus_sequences") or []]
    add("Analysis", "Marker databases", ", ".join(markers) if markers else "none")
    add("Analysis", "Reference loci", ", ".join(loci) if loci else "none")

    if manifest:
        add("Shared", "Exported", manifest.get("exported_at"))
        add("Shared", "Exported with", f"BACTolomew {manifest.get('app_version')}" if manifest.get("app_version") else None)
        add("Shared", "Package format", f"v{manifest.get('format_version')}")
        add("Shared", "Assemblies", "included" if manifest.get("includes_assemblies") else "not included")
        add("Shared", "Metadata", "CSV included" if manifest.get("metadata") == "csv" else "run data only")
    if imported:
        add("Shared", "Imported", imported.get("imported_at"))
        add("Shared", "Imported from", imported.get("source"))
        if imported.get("run_id_from"):
            add("Shared", "Original run id", imported.get("run_id_from"))
    return sections


__all__ = [
    "ARCHIVE_SUFFIX",
    "ASSEMBLY_PATTERN",
    "EXCLUDE",
    "IMPORTED",
    "ImportResult",
    "Issue",
    "MANIFEST",
    "METADATA_CSV",
    "SCHEME_SPECIES",
    "SHARE_DIR",
    "SHARE_FORMAT_VERSION",
    "SHARE_PATTERNS",
    "STAGING_PREFIX",
    "ShareCancelled",
    "ShareError",
    "TABLES_DIR",
    "check_compatibility",
    "current_versions",
    "estimate_share_size",
    "export_run",
    "format_size",
    "import_archive",
    "import_folder",
    "import_run_id",
    "is_imported_run",
    "is_shared_run",
    "parse_version",
    "read_manifest",
    "rebase_recorded_path",
    "recorded_root",
    "remove_imported_run",
    "run_id_taken",
    "run_info",
    "share_files",
    "warning_summary",
]
