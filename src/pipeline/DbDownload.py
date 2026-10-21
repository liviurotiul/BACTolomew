"""Install the databases whose tools download them themselves.

One function per database, each running the tool's own download command in the tools
image through ToolEnv (identity bind mounts, the configured runtime):

* ``bakta-full``: ``bakta_db download --type full`` (the light database is no longer
  installed; one installed earlier still works). bakta_db then
  runs ``amrfinder_update --force_update`` on the DB's ``amrfinderplus-db`` with the
  image's amrfinder, so the bundled AMRFinderPlus DB matches the image's software (a
  v5.1 DB downloaded earlier on the host carries a format 3.12 AMRF DB that
  amrfinder 4.0.23 rejects); the pipeline's amrfinder step uses this database. With
  ``amrfinder_db=`` a pinned release installed by ``install_amrfinderplus`` is copied in
  instead.
* ``checkm2``: ``checkm2 database --download --path``.
* ``amrfinderplus``: ``amrfinder_update --database`` (the newest release of the image's
  database format), or with ``release=`` NCBI's versioned database folder
  (``<format>/<release>/``) fetched over HTTPS and indexed with ``amrfinder_index``.
  Its ReferenceGeneCatalog.txt and AMR_CDS feed the AMRF marker files (MarkerDbBuild).
* ``mobsuite`` and ``mlst``: shipped in the image (mob_init ran at build time; the mlst
  package carries its PubMLST data); only verified.

Every install downloads into ``.<target>.partial`` next to the target and renames it into
place only when the download finished and the result was verified. An existing database
is never deleted or changed: by default it is kept and nothing is downloaded
(``if_exists="keep"``), or the new one goes beside it (``if_exists="beside"``, target
suffixed with the date). What was installed (tool, database version, date, path, size,
checksum or the tool's own version file) is written to ``bactolomew_install.json`` in the
target. No Qt; progress lines go to a callback.

Usage:
    PYTHONPATH=src python -m pipeline.DbDownload --list
    PYTHONPATH=src python -m pipeline.DbDownload --install bakta-full checkm2 --dest /data/tooldbs [--dry-run]
    PYTHONPATH=src python -m pipeline.DbDownload --install bakta-full --set-default
    PYTHONPATH=src python -m pipeline.DbDownload --install amrfinderplus --amrf-release 2025-07-16.1 --dry-run
    PYTHONPATH=src python -m pipeline.DbDownload --verify mobsuite mlst
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
import threading
import urllib.request
from dataclasses import dataclass
from datetime import datetime
from functools import wraps
from pathlib import Path
from typing import Any, Callable, Mapping, Optional, Sequence

from pipeline.ContainerRuntime import _stream, read_image_metadata
from pipeline.ToolEnv import (
    ContainerContext,
    _container_prefix,
    _load_app_config,
    _tool_env_key,
    get_tool_env_prefix,
    make_container_context,
    project_root,
    save_app_config,
)

LogFn = Callable[[str], None]

INSTALL_RECORD = "bactolomew_install.json"
PARTIAL_MARKER = ".bactolomew_partial"
AMRF_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database"
# Database format of the image's amrfinder (ncbi-amrfinderplus 4.0.23, docker/locks).
AMRF_DEFAULT_FORMAT = "4.0"
CHECKM2_DMND = Path("CheckM2_database") / "uniref100.KO.1.dmnd"
_GB = 1_000_000_000  # sizes are shown in GB (10^9 bytes) everywhere

# What every install asks for on top of its peak: the AMRFinderPlus update (about 0.3 GB a
# release), the install record, logs and some room left on the disk.
SPACE_MARGIN_BYTES = 5 * _GB
# The full Bakta DB that the image's bakta (1.10.4, DB schema v5) downloads: db.tar.gz of
# the v5.1 record (zenodo 10522951, 39.4 GB), unpacked into db/ with its AMRFinderPlus DB
# (76.3 GB, measured on a v5.1 install). bakta_db keeps the archive until it is unpacked,
# so the peak is 115.7 GB; with the margin the check asks for 120.7 GB.
BAKTA_FULL_ARCHIVE_BYTES = 39_401_553_410
BAKTA_FULL_UNPACKED_BYTES = 76_302_905_767
BAKTA_FULL_PEAK_BYTES = BAKTA_FULL_ARCHIVE_BYTES + BAKTA_FULL_UNPACKED_BYTES
BAKTA_FULL_NEED_BYTES = BAKTA_FULL_PEAK_BYTES + SPACE_MARGIN_BYTES
# CheckM2: a 1.71 GB archive unpacked to 3.08 GB (measured), kept side by side: 4.8 GB.
CHECKM2_ARCHIVE_BYTES = 1_710_000_000
CHECKM2_PEAK_BYTES = CHECKM2_ARCHIVE_BYTES + 3_082_508_797
# The tools image as a Singularity/Apptainer .sif: the docker image is 22 to 29 GB
# unpacked; its squashfs is estimated (not measured) at up to 15 GB.
SIF_NEED_BYTES = 15 * _GB


@dataclass(frozen=True)
class DbSpec:
    db_id: str
    label: str
    tool: str               # ToolEnv tool name (picks the conda env in the image)
    default_dirname: str    # target folder name under --dest
    need_bytes: int         # peak free space the install needs (archive + unpacked)
    installable: bool = True
    download_bytes: int = 0  # about what the tool downloads (shown before an install)


SPECS: dict[str, DbSpec] = {s.db_id: s for s in (
    DbSpec("bakta-full", "Bakta DB (full)", "bakta", "bakta-full", BAKTA_FULL_NEED_BYTES,
           download_bytes=BAKTA_FULL_ARCHIVE_BYTES),
    DbSpec("checkm2", "CheckM2 DIAMOND DB", "checkm2", "checkm2", CHECKM2_PEAK_BYTES + SPACE_MARGIN_BYTES,
           download_bytes=CHECKM2_ARCHIVE_BYTES),
    DbSpec("amrfinderplus", "AMRFinderPlus DB", "amrfinder", "amrfinderplus", 1 * _GB),
    DbSpec("mobsuite", "MOB-suite DB (in the image)", "mob_init", "", 0, installable=False),
    DbSpec("mlst", "MLST PubMLST data (in the mlst package)", "mlst", "", 0, installable=False),
)}


class DbInstallError(RuntimeError):
    pass


class NotEnoughSpace(DbInstallError):
    def __init__(self, folder: Path, free_bytes: int, need_bytes: int, detail: str = "") -> None:
        super().__init__(
            f"Not enough free space in {folder}: {free_bytes / _GB:.1f} GB free, about "
            f"{need_bytes / _GB:.0f} GB needed{f' ({detail})' if detail else ''}"
        )
        self.folder = folder
        self.free_bytes = free_bytes
        self.need_bytes = need_bytes


def _log(log: Optional[LogFn], msg: str) -> None:
    (log or print)(msg)


# --------------------------------------------------------------------------- #
# Commands through ToolEnv
# --------------------------------------------------------------------------- #
# The temporary container HOMEs (make_container_context without output_dir) made during
# one public call on this thread; removed when that call returns (_drops_temp_homes).
_TEMP_HOMES = threading.local()


def _drops_temp_homes(fn: Callable[..., Any]) -> Callable[..., Any]:
    @wraps(fn)
    def call(*args: Any, **kwargs: Any) -> Any:
        depth = getattr(_TEMP_HOMES, "depth", 0)
        if depth == 0:
            _TEMP_HOMES.paths = []
        _TEMP_HOMES.depth = depth + 1
        try:
            return fn(*args, **kwargs)
        finally:
            _TEMP_HOMES.depth = depth
            if depth == 0:
                for home in _TEMP_HOMES.paths:
                    shutil.rmtree(home, ignore_errors=True)
                _TEMP_HOMES.paths = []
    return call


def _context(config: Optional[Mapping[str, Any]], paths: Sequence[Path], workdir: Path) -> ContainerContext:
    ctx = make_container_context(config, extra_paths=[str(p) for p in paths], workdir=str(workdir))
    home = Path(ctx.home_dir)
    if (getattr(_TEMP_HOMES, "depth", 0) and home.name.startswith("bactolomew-home-")
            and home.parent == Path(tempfile.gettempdir())):
        _TEMP_HOMES.paths.append(home)
    return ctx


def tool_argv(ctx: ContainerContext, tool: str, exe: str, args: Sequence[Any], *, workdir: Path) -> list[str]:
    """`exe args` in the tool's env: inside the image for docker/singularity, through the
    configured conda env for native. Docker runs get `--init`, so a stop (SIGTERM, which
    `docker run` passes on) ends the container's tool too."""
    if ctx.runtime == "native":
        prefix = get_tool_env_prefix(tool)
        base = ["conda", "run", "--no-capture-output", "-p", prefix, exe] if prefix else [exe]
    else:
        base = _container_prefix(ctx, _tool_env_key(tool), exe, str(workdir))
        if ctx.runtime in ("docker", "podman") and base[1:2] == ["run"]:
            base.insert(2, "--init")
    return [*base, *(str(a) for a in args)]


def _run(cmd: list[str], log: Optional[LogFn], what: str) -> None:
    rc = _stream(cmd, log)
    if rc != 0:
        raise DbInstallError(f"{what} failed (exit {rc})")


# --------------------------------------------------------------------------- #
# Partial folder, final placement, install record
# --------------------------------------------------------------------------- #
def _non_empty_dir(path: Path) -> bool:
    return path.is_dir() and any(path.iterdir())


def resolve_target(target: Path, if_exists: str) -> Optional[Path]:
    """The folder to install into, or None when `target` holds data and is kept."""
    target = Path(target).expanduser().absolute()
    if not target.exists() or (target.is_dir() and not any(target.iterdir())):
        return target
    if if_exists == "beside":
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        return target.with_name(f"{target.name}-{stamp}")
    if if_exists == "keep":
        return None
    raise DbInstallError(f"{target} exists; choose to keep it or install beside it")


def partial_dir(target: Path) -> Path:
    return target.with_name(f".{target.name}.partial")


def _free_bytes(path: Path) -> int:
    return shutil.disk_usage(path).free


def _open_partial(target: Path, need_bytes: int, log: Optional[LogFn]) -> Path:
    part = partial_dir(target)
    if part.exists():
        if not (part / PARTIAL_MARKER).is_file():
            raise DbInstallError(f"{part} exists and was not made by this installer; move it away first")
        _log(log, f"Removing the leftover partial download {part}")
        shutil.rmtree(part)
    part.parent.mkdir(parents=True, exist_ok=True)
    free = _free_bytes(part.parent)
    if need_bytes and free < need_bytes:
        raise NotEnoughSpace(part.parent, free, need_bytes)
    part.mkdir()
    (part / PARTIAL_MARKER).write_text(datetime.now().astimezone().isoformat(), encoding="utf-8")
    return part


def _discard_partial(part: Path) -> None:
    if (part / PARTIAL_MARKER).is_file():
        shutil.rmtree(part, ignore_errors=True)


def _tree_size(path: Path) -> int:
    total = 0
    for root, _dirs, files in os.walk(path):
        for name in files:
            try:
                total += (Path(root) / name).lstat().st_size
            except OSError:
                pass
    return total


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def _finish(part: Path, target: Path, record: dict[str, Any], log: Optional[LogFn]) -> dict[str, Any]:
    """Write the record into the verified partial folder and rename it into place."""
    (part / PARTIAL_MARKER).unlink(missing_ok=True)
    record = {
        **record,
        "installed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "size_bytes": _tree_size(part),
    }
    final = dict(record, path=_rebase(record.get("path"), part, target))
    (part / INSTALL_RECORD).write_text(json.dumps(final, indent=2), encoding="utf-8")
    if target.exists():  # an empty folder only: resolve_target let it through
        target.rmdir()
    part.rename(target)
    _log(log, f"Installed {record.get('database')} at {final['path']}")
    return final


def _rebase(path: Any, part: Path, target: Path) -> str:
    p = Path(str(path or part))
    try:
        return str(target / p.relative_to(part))
    except ValueError:
        return str(p)


def read_install_record(target: Path) -> Optional[dict]:
    try:
        raw = json.loads((Path(target) / INSTALL_RECORD).read_text(encoding="utf-8"))
        return raw if isinstance(raw, dict) else None
    except Exception:
        return None


def _kept(target: Path, path: Path, database: str, log: Optional[LogFn]) -> dict[str, Any]:
    _log(log, f"{database}: keeping the existing database at {target} (nothing downloaded)")
    return read_install_record(target) or {"database": database, "path": str(path), "kept": True}


def _dry(log: Optional[LogFn], cmds: Sequence[Sequence[str]], target: Path, extra: str = "") -> dict[str, Any]:
    _log(log, f"[dry-run] target {target} (via {partial_dir(target)})")
    for cmd in cmds:
        _log(log, "[dry-run] $ " + " ".join(cmd))
    if extra:
        _log(log, f"[dry-run] {extra}")
    return {"dry_run": True, "target": str(target), "commands": [list(c) for c in cmds]}


def _image_ref(ctx: ContainerContext) -> str:
    return ctx.sif_path if ctx.runtime == "singularity" else (ctx.image if ctx.runtime in ("docker", "podman") else "native")


# --------------------------------------------------------------------------- #
# AMRFinderPlus database checks (shared by the Bakta and the standalone install)
# --------------------------------------------------------------------------- #
def _amrf_check_cmd(ctx: ContainerContext, db_dir: Path, workdir: Path) -> list[str]:
    # Exit 0 only when the image's amrfinder can load this database (format and indexes).
    return tool_argv(ctx, "amrfinder", "amrfinder", ["--database", db_dir, "--database_version"], workdir=workdir)


def _amrf_release(db_dir: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for key, name in (("release", "version.txt"), ("format", "database_format_version.txt")):
        f = db_dir / name
        if f.is_file():
            out[key] = f.read_text(encoding="utf-8", errors="replace").strip()
    return out


def _copy_amrf_release(src: Path, amrf_root: Path, log: Optional[LogFn]) -> str:
    """Put an installed AMRF release (the folder holding version.txt) into `amrf_root` as
    `<release>/` plus a `latest` link, the layout amrfinder_update writes."""
    src = Path(src).expanduser().resolve()
    release = _amrf_release(src).get("release")
    if not release:
        raise DbInstallError(f"{src} is not an AMRFinderPlus database (no version.txt)")
    amrf_root.mkdir(parents=True, exist_ok=True)
    dest = amrf_root / release
    if dest.exists():  # inside our partial folder: bakta_db's own copy of this release
        shutil.rmtree(dest)
    _log(log, f"Copying AMRFinderPlus release {release} from {src}")
    shutil.copytree(src, dest, symlinks=False)
    latest = amrf_root / "latest"
    if latest.is_symlink() or latest.exists():
        latest.unlink()
    latest.symlink_to(release)
    return release


# --------------------------------------------------------------------------- #
# Installers
# --------------------------------------------------------------------------- #
def _bakta_db_dir(base: Path) -> Optional[Path]:
    # db-light: a light database installed before the full one became the only choice.
    for cand in (base / "db", base / "db-light", base):
        if (cand / "version.json").is_file():
            return cand
    return None


@_drops_temp_homes
def install_bakta(
    target: Path,
    *,
    amrfinder_db: Optional[Path] = None,
    if_exists: str = "keep",
    config: Optional[Mapping[str, Any]] = None,
    log: Optional[LogFn] = None,
    dry_run: bool = False,
) -> dict[str, Any]:
    """`bakta_db download --type full` into `target`; record["path"] is the Bakta DB folder
    (`target/db`) to put in the config's `bakta_db`."""
    spec = SPECS["bakta-full"]
    requested = Path(target).expanduser().absolute()
    final = resolve_target(requested, if_exists)
    if final is None:
        return _kept(requested, _bakta_db_dir(requested) or requested, spec.label, log)
    part = partial_dir(final)
    ctx = _context(config, [final.parent], part)
    inner = part / "db"
    download = tool_argv(ctx, "bakta", "bakta_db", ["download", "--output", part, "--type", "full"], workdir=part)
    check = _amrf_check_cmd(ctx, inner / "amrfinderplus-db" / "latest", part)
    if dry_run:
        extra = f"then copy the AMRFinderPlus release {amrfinder_db} into {inner}/amrfinderplus-db" if amrfinder_db else ""
        return _dry(log, [download, check], final, extra)
    part = _open_partial(final, spec.need_bytes, log)
    try:
        _log(log, f"Downloading the {spec.label} into {part} ...")
        _run(download, log, "bakta_db download")
        db_dir = _bakta_db_dir(part)
        if db_dir is None or not (db_dir / "bakta.db").is_file():
            raise DbInstallError(f"bakta_db download finished but no Bakta DB was found under {part}")
        amrf_root = db_dir / "amrfinderplus-db"
        if amrfinder_db:
            _copy_amrf_release(Path(amrfinder_db), amrf_root, log)
        _log(log, "Checking the bundled AMRFinderPlus DB with the image's amrfinder ...")
        _run(_amrf_check_cmd(ctx, amrf_root / "latest", part), log, "AMRFinderPlus DB check")
        version = json.loads((db_dir / "version.json").read_text(encoding="utf-8"))
        record = {
            "database": spec.db_id,
            "tool": "bakta_db",
            "image": _image_ref(ctx),
            "path": str(db_dir),
            "version": {k: version.get(k) for k in ("major", "minor", "date", "type", "doi")},
            "version_file": "version.json",
            "amrfinderplus_db": _amrf_release(amrf_root / "latest"),
        }
        return _finish(part, final, record, log)
    except BaseException:
        _discard_partial(part)
        raise


@_drops_temp_homes
def install_checkm2(
    target: Path,
    *,
    if_exists: str = "keep",
    config: Optional[Mapping[str, Any]] = None,
    log: Optional[LogFn] = None,
    dry_run: bool = False,
    checksum: bool = True,
) -> dict[str, Any]:
    """`checkm2 database --download` into `target`; record["path"] is the DIAMOND file for
    the config's `assembly_qc.checkm2_database_path`."""
    spec = SPECS["checkm2"]
    requested = Path(target).expanduser().absolute()
    final = resolve_target(requested, if_exists)
    if final is None:
        return _kept(requested, requested / CHECKM2_DMND, spec.label, log)
    part = partial_dir(final)
    ctx = _context(config, [final.parent], part)
    # --no_write_json_db: the image's checkm2 install folder is read-only; runs pass the path.
    download = tool_argv(ctx, "checkm2", "checkm2",
                         ["database", "--download", "--path", part, "--no_write_json_db"], workdir=part)
    if dry_run:
        return _dry(log, [download], final)
    part = _open_partial(final, spec.need_bytes, log)
    try:
        _log(log, f"Downloading the {spec.label} into {part} ...")
        _run(download, log, "checkm2 database --download")
        dmnd = part / CHECKM2_DMND
        if not dmnd.is_file():
            found = sorted(part.rglob("*.dmnd"))
            if not found:
                raise DbInstallError(f"checkm2 finished but no DIAMOND database was found under {part}")
            dmnd = found[0]
        if dmnd.stat().st_size == 0:
            raise DbInstallError(f"{dmnd} is empty")
        record = {
            "database": spec.db_id,
            "tool": "checkm2 database --download",
            "image": _image_ref(ctx),
            "path": str(dmnd),
            "file_size": dmnd.stat().st_size,
        }
        if checksum:
            _log(log, f"Hashing {dmnd.name} ...")
            record["sha256"] = _sha256(dmnd)
        return _finish(part, final, record, log)
    except BaseException:
        _discard_partial(part)
        raise


_HREF_RE = re.compile(r'href="([^"?/][^"]*)"')


def _amrf_release_files(url: str) -> list[str]:
    """File names in NCBI's HTML listing of one release folder."""
    with urllib.request.urlopen(url, timeout=60) as resp:
        html = resp.read().decode("utf-8", errors="replace")
    names = [n for n in _HREF_RE.findall(html) if "/" not in n and not n.startswith(("http", "#"))]
    return sorted(set(names))


def _fetch(url: str, dest: Path) -> None:
    tmp = dest.with_name(dest.name + ".part")
    with urllib.request.urlopen(url, timeout=120) as resp, open(tmp, "wb") as fh:
        shutil.copyfileobj(resp, fh, 1 << 20)
    tmp.rename(dest)


@_drops_temp_homes
def install_amrfinderplus(
    target: Path,
    *,
    release: str = "",
    db_format: str = AMRF_DEFAULT_FORMAT,
    if_exists: str = "keep",
    config: Optional[Mapping[str, Any]] = None,
    log: Optional[LogFn] = None,
    dry_run: bool = False,
) -> dict[str, Any]:
    """The AMRFinderPlus DB in amrfinder_update's layout (`target/<release>/`, `latest`).
    Without `release`: `amrfinder_update` (newest release of the image's format). With it:
    NCBI's `<db_format>/<release>/` folder over HTTPS, then `amrfinder_index`."""
    spec = SPECS["amrfinderplus"]
    requested = Path(target).expanduser().absolute()
    final = resolve_target(requested, if_exists)
    if final is None:
        return _kept(requested, requested / "latest", spec.label, log)
    part = partial_dir(final)
    ctx = _context(config, [final.parent], part)
    check = _amrf_check_cmd(ctx, part / "latest", part)
    if release:
        url = f"{AMRF_BASE_URL}/{db_format}/{release}/"
        index = tool_argv(ctx, "amrfinder", "amrfinder_index", [part / release], workdir=part)
        cmds = [index, check]
        extra = f"files listed at {url} are downloaded into {part / release} before amrfinder_index"
    else:
        url = ""
        cmds = [tool_argv(ctx, "amrfinder", "amrfinder_update", ["--database", part, "--force_update"], workdir=part), check]
        extra = ""
    if dry_run:
        return _dry(log, cmds, final, extra)
    part = _open_partial(final, spec.need_bytes, log)
    try:
        if release:
            rel_dir = part / release
            rel_dir.mkdir()
            names = _amrf_release_files(url)
            if "version.txt" not in names:
                raise DbInstallError(f"{url} does not list an AMRFinderPlus release (no version.txt)")
            for i, name in enumerate(names, 1):
                _log(log, f"[{i}/{len(names)}] {url}{name}")
                _fetch(url + name, rel_dir / name)
            got = _amrf_release(rel_dir).get("release")
            if got != release:
                raise DbInstallError(f"{url}version.txt says {got!r}, expected {release!r}")
            (part / "latest").symlink_to(release)
            _run(cmds[0], log, "amrfinder_index")
        else:
            _log(log, f"Downloading the newest {spec.label} release into {part} ...")
            _run(cmds[0], log, "amrfinder_update")
        _log(log, "Checking the database with the image's amrfinder ...")
        _run(check, log, "AMRFinderPlus DB check")
        rel_dir = (part / "latest").resolve()
        info = _amrf_release(rel_dir)
        record = {
            "database": spec.db_id,
            "tool": "amrfinder_index" if release else "amrfinder_update",
            "image": _image_ref(ctx),
            "path": str(part / info.get("release", "latest")),
            "version": info,
            "version_file": "version.txt",
            "pinned": bool(release),
            "source": url or f"{AMRF_BASE_URL}/ (amrfinder_update)",
            # Inputs of the AMRF marker FASTAs (MarkerDbBuild, DB/metadata/AMRF).
            "marker_inputs": {
                name: str(part / info.get("release", "latest") / name)
                for name in ("ReferenceGeneCatalog.txt", "AMR_CDS.fa", "AMRProt.fa")
                if (rel_dir / name).is_file()
            },
        }
        record["marker_inputs"] = {k: _rebase(v, part, final) for k, v in record["marker_inputs"].items()}
        return _finish(part, final, record, log)
    except BaseException:
        _discard_partial(part)
        raise


@_drops_temp_homes
def verify_mobsuite(*, config: Optional[Mapping[str, Any]] = None, log: Optional[LogFn] = None,
                    dry_run: bool = False) -> dict[str, Any]:
    """The MOB-suite DB is baked into the image by mob_init (or an external
    mobsuite.database_directory): report it, download nothing."""
    cfg = dict(config) if isinstance(config, Mapping) else _load_app_config()
    mob = cfg.get("mobsuite") if isinstance(cfg.get("mobsuite"), Mapping) else {}
    external = str(mob.get("database_directory") or "").strip()
    if external:
        status = Path(external).expanduser() / "status.txt"
        ok = status.is_file()
        _log(log, f"MOB-suite: external database {external}: {'ok' if ok else 'status.txt missing'}")
        return {"database": "mobsuite", "location": "external", "path": external, "ok": ok}
    ctx = _context(cfg, [], Path("/tmp"))
    if dry_run or ctx.runtime == "native":
        _log(log, "[dry-run] read the image's MOB-suite databases/status.txt" if dry_run
             else "MOB-suite: native runtime, the conda env's mob_init database is used")
        return {"database": "mobsuite", "location": "image", "dry_run": dry_run}
    found = read_image_metadata(ctx)["databases"].get("mobsuite")
    _log(log, f"MOB-suite: {'database in the image, ' + json.dumps(found) if found else 'no database found in the image'}")
    return {"database": "mobsuite", "location": "image", "ok": bool(found), **(found or {})}


@_drops_temp_homes
def verify_mlst(*, config: Optional[Mapping[str, Any]] = None, log: Optional[LogFn] = None,
                dry_run: bool = False) -> dict[str, Any]:
    """The PubMLST snapshot ships in the mlst package: report its version."""
    workdir = Path("/tmp")
    ctx = _context(config, [], workdir)
    cmd = tool_argv(ctx, "mlst", "mlst", ["--version"], workdir=workdir)
    if dry_run:
        _log(log, "[dry-run] $ " + " ".join(cmd))
        return {"database": "mlst", "dry_run": True, "commands": [cmd]}
    res = subprocess.run(cmd, capture_output=True, text=True)
    text = (res.stdout or res.stderr or "").strip().splitlines()
    version = text[-1] if text else ""
    _log(log, f"MLST: {version or 'mlst did not run'} (PubMLST data inside the package)")
    return {"database": "mlst", "location": "image", "ok": res.returncode == 0, "version": version}


# `<exe> --version` of the tools whose databases first-time setup installs.
TOOL_CHECKS = (("bakta", "bakta"), ("checkm2", "checkm2"))


@_drops_temp_homes
def verify_tools(*, config: Optional[Mapping[str, Any]] = None, log: Optional[LogFn] = None) -> None:
    """Start each tool of TOOL_CHECKS once (`--version`) in the tools image, or its conda env
    for native; raises DbInstallError naming the first one that does not start."""
    workdir = Path("/tmp")
    ctx = _context(config, [], workdir)
    for tool, exe in TOOL_CHECKS:
        _run(tool_argv(ctx, tool, exe, ["--version"], workdir=workdir), log, f"{exe} --version")


def install(db_id: str, dest_root: Path, **kwargs: Any) -> dict[str, Any]:
    """Install one database by id under `dest_root/<default folder>`."""
    spec = SPECS.get(db_id)
    if spec is None:
        raise DbInstallError(f"unknown database {db_id!r}; known: {', '.join(SPECS)}")
    if not spec.installable:
        verify = verify_mobsuite if db_id == "mobsuite" else verify_mlst
        return verify(config=kwargs.get("config"), log=kwargs.get("log"), dry_run=kwargs.get("dry_run", False))
    target = Path(dest_root) / spec.default_dirname
    if db_id.startswith("bakta-"):
        return install_bakta(target, **kwargs)
    if db_id == "checkm2":
        return install_checkm2(target, **kwargs)
    return install_amrfinderplus(target, **kwargs)


def config_update(db_id: str, record: Mapping[str, Any]) -> Optional[dict[str, Any]]:
    """The bactolomew_config.json keys that make an installed (or, for a dry run, planned)
    database the one runs use: `bakta_db` or `assembly_qc.checkm2_database_path`; None
    for a database the config does not name."""
    if record.get("dry_run"):
        target = Path(str(record.get("target") or ""))
        path = target / "db" if db_id.startswith("bakta-") else target / CHECKM2_DMND
    else:
        path = Path(str(record.get("path") or ""))
    if db_id.startswith("bakta-"):
        return {"bakta_db": str(path)}
    if db_id == "checkm2":
        return {"assembly_qc": {"checkm2_database_path": str(path)}}
    return None


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Install the databases whose tools download them themselves.")
    parser.add_argument("--list", action="store_true", help="list the databases and exit")
    parser.add_argument("--install", nargs="+", metavar="ID", default=[], help="database ids to install")
    parser.add_argument("--verify", nargs="+", metavar="ID", default=[], help="mobsuite / mlst checks")
    parser.add_argument("--dest", type=Path, default=project_root() / "tool_dbs",
                        help="folder that receives one subfolder per database (default: <repo>/tool_dbs)")
    parser.add_argument("--beside", action="store_true",
                        help="when a target already holds a database, install the new one beside it")
    parser.add_argument("--amrf-release", default="", help="pinned AMRFinderPlus release, e.g. 2025-07-16.1")
    parser.add_argument("--amrf-format", default=AMRF_DEFAULT_FORMAT, help="AMRFinderPlus database format folder")
    parser.add_argument("--bakta-amrf", type=Path, default=None,
                        help="installed AMRFinderPlus release folder to copy into the Bakta DB")
    parser.add_argument("--set-default", action="store_true",
                        help="make the installed Bakta or CheckM2 database the one runs use "
                             "(writes bakta_db / assembly_qc.checkm2_database_path in bactolomew_config.json)")
    parser.add_argument("--dry-run", action="store_true", help="print the exact commands, change nothing")
    args = parser.parse_args(argv)

    if args.list or not (args.install or args.verify):
        for spec in SPECS.values():
            how = f"-> {args.dest / spec.default_dirname}" if spec.installable else "(verify only)"
            print(f"{spec.db_id:14s} {spec.label:42s} {how}")
        return 0
    if_exists = "beside" if args.beside else "keep"
    for db_id in [*args.install, *args.verify]:
        kwargs: dict[str, Any] = {"log": print, "dry_run": args.dry_run}
        spec = SPECS.get(db_id)
        if spec is not None and spec.installable:
            kwargs["if_exists"] = if_exists
            if db_id == "amrfinderplus":
                kwargs.update(release=args.amrf_release, db_format=args.amrf_format)
            elif db_id.startswith("bakta-") and args.bakta_amrf:
                kwargs["amrfinder_db"] = args.bakta_amrf
        try:
            record = install(db_id, args.dest, **kwargs)
        except DbInstallError as exc:
            print(f"ERROR {db_id}: {exc}")
            return 1
        if not args.dry_run:
            print(json.dumps(record, indent=2))
        update = config_update(db_id, record) if args.set_default else None
        if update is not None:
            if args.dry_run:
                print(f"[dry-run] would set in bactolomew_config.json: {json.dumps(update)}")
            else:
                save_app_config(update)
                print(f"Set in bactolomew_config.json: {json.dumps(update)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
