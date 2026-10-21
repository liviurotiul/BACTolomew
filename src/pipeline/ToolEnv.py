from __future__ import annotations

import contextvars
import json
import logging
import os
import re
import shlex
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, List, Mapping, Optional, Tuple

from pipeline.Resources import detect_budget

REQUIRED_TOOL_ENVS: tuple[str, ...] = (
    "bakta",
    "diamond",
    "mobsuite",
    "phylogeny",
    "gubbins",
    "panaroo",
    "chewbbaca",
    "mlst",
    "assembly_qc",
)

DEV_BAKTA_ENV_PREFIX = "/Shared/bakta"
DEV_BAKTA_DB_PATH = "/Shared/bakta_db_full/db"
DEFAULT_FASTGEAR_DOCKER_BIN = "docker"
# The third-party fastGEAR image (last pushed 2022-12-01) that docker/Dockerfile copies
# into the tools image, pinned by its Docker Hub digest. A local `docker save` archive
# (config container.fastgear_archive) keeps it buildable if Docker Hub loses it.
DEFAULT_FASTGEAR_IMAGE = (
    "pgcbioinfo/fastgear:latest@sha256:afb299e7fb327129ebf2acc1bcae00e0a54cd6fbe71d114d0a9f1642a2216b79"
)
# The unpinned reference older configs carry (fastgear.image): read as the pinned default.
LEGACY_FASTGEAR_IMAGE = "pgcbioinfo/fastgear:latest"
DEFAULT_FASTGEAR_LAUNCHER = "/opt/run_fastGEAR.sh"
DEFAULT_FASTGEAR_MCR_ROOT = "/opt/v901"
DEFAULT_FASTGEAR_SPECS_FILE = "/opt/fG_input_specs.txt"
DEFAULT_FASTGEAR_WORKDIR = "/work"
CHECKM2_DB_RELATIVE_PATH = Path("checkm2_db") / "CheckM2_database" / "uniref100.KO.1.dmnd"
CHECKM2_LOWMEM_BELOW_GB = 32  # see configured_assembly_qc_settings


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_app_config() -> dict:
    path = project_root() / "bactolomew_config.json"
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


def _deep_merge(base: dict, updates: Mapping[str, Any]) -> dict:
    for key, value in updates.items():
        if isinstance(value, Mapping) and isinstance(base.get(key), Mapping):
            base[key] = _deep_merge(dict(base[key]), value)
        else:
            base[key] = value
    return base


def save_app_config(updates: Mapping[str, Any]) -> dict:
    """Merge `updates` into bactolomew_config.json and write it back atomically.

    Unlike first-time setup (which rewrites the whole file), this preserves every
    existing key: nested dicts are merged, so e.g. updating resources.max_threads
    leaves resources.max_ram_gb and all other blocks untouched.
    """
    path = project_root() / "bactolomew_config.json"
    try:
        current = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
        if not isinstance(current, dict):
            current = {}
    except Exception:
        current = {}
    merged = _deep_merge(current, dict(updates))
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(merged, indent=2, sort_keys=True), encoding="utf-8")
    os.replace(tmp, path)  # atomic on POSIX
    return merged


def configured_resources() -> dict[str, int]:
    """Global compute limits from config `resources` block, with sane defaults.

    max_threads defaults to every core this process may use (no cap); max_ram_gb
    defaults to 0 (no limit).
    """
    config = _load_app_config()
    raw = config.get("resources") if isinstance(config.get("resources"), Mapping) else {}
    cpu = detect_budget().cores

    def _pos_int(value: Any, default: int) -> int:
        try:
            n = int(value)
        except (TypeError, ValueError):
            return default
        return n if n > 0 else default

    return {
        "max_threads": _pos_int(raw.get("max_threads"), cpu),
        "max_ram_gb": _pos_int(raw.get("max_ram_gb"), 0),
    }


def effective_core_budget() -> int:
    """The core budget a new run should plan against: min(usable cores, config cap). Usable
    cores are the machine's, capped by the CPU affinity, the cgroup CPU quota and the
    SLURM allocation (pipeline.Resources)."""
    return max(1, min(detect_budget().cores, configured_resources()["max_threads"]))


def configured_max_ram_gb() -> int:
    """User-configured RAM cap in GB, or 0 when unset (no enforcement)."""
    return configured_resources()["max_ram_gb"]


def effective_ram_gb() -> Optional[float]:
    """The RAM a run may use in GB: the detected memory (machine, cgroup limit, SLURM
    allocation) capped by resources.max_ram_gb; None when neither is known."""
    detected = detect_budget().memory_gb
    cap = configured_max_ram_gb()
    known = [v for v in (detected, float(cap) if cap else None) if v]
    return min(known) if known else None


def resource_budget_text() -> str:
    """One line for a run's log: the core and RAM budget and what limits each."""
    budget = detect_budget()
    cores, cap = effective_core_budget(), configured_resources()["max_threads"]
    text = budget.describe()
    if cores < budget.cores and cap == cores:
        text += f"; resources.max_threads caps it at {cores} cores"
    ram_cap = configured_max_ram_gb()
    if ram_cap and (budget.memory_gb is None or ram_cap < budget.memory_gb):
        text += f"; resources.max_ram_gb caps it at {ram_cap} GB"
    return text


def _normalize_shell_tokens(value: Any) -> list[str]:
    if isinstance(value, str):
        return [token for token in shlex.split(value) if token]
    if isinstance(value, (list, tuple)):
        out: list[str] = []
        for item in value:
            text = str(item).strip()
            if text:
                out.append(text)
        return out
    return []


def _tool_env_key(tool_name: str) -> str:
    key_map = {
        "bakta": "bakta",
        "amrfinder": "bakta",
        "diamond": "diamond",
        "mob_recon": "mobsuite",
        "mob_typer": "mobsuite",
        "mob_init": "mobsuite",
        "mobsuite": "mobsuite",
        "snippy": "phylogeny",
        "snippy-core": "phylogeny",
        "raxml-ng": "phylogeny",
        "run_gubbins.py": "gubbins",
        "gubbins": "gubbins",
        "iqtree2": "gubbins",
        "mask_gubbins_aln.py": "gubbins",
        "minimap2": "phylogeny",
        "panaroo": "panaroo",
        "chewbbaca": "chewbbaca",
        "mlst": "mlst",
        "quast": "assembly_qc",
        "quast.py": "assembly_qc",
        "checkm2": "assembly_qc",
    }
    return key_map.get(tool_name, tool_name)


def get_tool_env_prefix(tool_name: str) -> Optional[str]:
    config = _load_app_config()
    envs = config.get("envs")
    env_key = _tool_env_key(tool_name)
    dev_mode = bool(config.get("dev_mode"))

    if isinstance(envs, dict):
        raw = envs.get(env_key)
        if isinstance(raw, str) and raw.strip():
            return raw.strip()

    project_env = project_root() / "envs" / env_key
    if (project_env / "conda-meta").is_dir():
        return str(project_env)

    # Dev fallback used by current setup.
    if dev_mode and env_key in {"bakta", "diamond"}:
        return DEV_BAKTA_ENV_PREFIX
    return None


# ---------------------------------------------------------------------------
# Containerized execution (Docker / Podman / Singularity) of the pipeline tools.
#
# Every tool call funnels through tool_cmd()/conda_tool_cmd(). Instead of
# `conda run -p <env>` on the host, we wrap `conda run -n <env> <exe>` in a
# `docker run` / `singularity exec` against one image that carries every tool,
# bind-mounting the host dirs a run touches at IDENTICAL paths so the existing
# absolute-path command args resolve unchanged inside the container.
# ---------------------------------------------------------------------------

TOOLS_IMAGE = "bactolomew-tools:0.2.0"
DEFAULT_RUNTIME = "native"  # backward-compatible default when config has no `runtime`
CONTAINER_RUNTIMES = ("docker", "podman", "singularity")  # config `runtime` values besides native
OCI_RUNTIMES = ("docker", "podman")  # engines that run `<exe> run ... <image> <argv>`


@dataclass(frozen=True)
class ContainerContext:
    runtime: str                        # "docker" | "podman" | "singularity" | "native"
    image: str = ""                     # docker image ref
    sif_path: str = ""                  # singularity .sif path
    bind_roots: Tuple[str, ...] = ()    # identity mount roots
    uid: int = 0
    gid: int = 0
    home_dir: str = "/tmp"              # writable HOME inside the container
    default_workdir: str = "/tmp"
    mem_limit_gb: int = 0              # per-container RAM cap in GB (0 = no limit)
    shm_size: str = ""                 # docker --shm-size ("" = docker's default)
    engine_exe: str = ""               # docker, podman, singularity or apptainer ("" = the runtime's name)
    podman: bool = False               # Podman runs the containers (runtime podman, or the podman-docker shim)
    rootless: bool = False             # rootless Docker or Podman (see container_user_args)
    selinux: bool = False              # SELinux is enforcing on the host


_CURRENT_CTX: "contextvars.ContextVar[Optional[ContainerContext]]" = contextvars.ContextVar(
    "bactolomew_container_ctx", default=None
)
# Module-global mirror of the active context. contextvars do NOT propagate into
# ThreadPoolExecutor worker threads (several workflows parallelize per-genome), so
# the global is what those threads read. Safe because only one run executes at a time.
_ACTIVE_CTX: Optional[ContainerContext] = None


def set_container_context(ctx: Optional[ContainerContext]) -> "contextvars.Token":
    global _ACTIVE_CTX
    _ACTIVE_CTX = ctx
    return _CURRENT_CTX.set(ctx)


def reset_container_context(token: "contextvars.Token") -> None:
    global _ACTIVE_CTX
    _ACTIVE_CTX = None
    try:
        _CURRENT_CTX.reset(token)
    except Exception:
        pass


def current_container_context() -> Optional[ContainerContext]:
    return _CURRENT_CTX.get() or _ACTIVE_CTX


def _is_within(path: Any, root: Any) -> bool:
    try:
        Path(str(path)).resolve().relative_to(Path(str(root)).resolve())
        return True
    except Exception:
        return False


def compute_bind_roots(paths: Iterable[Any]) -> List[str]:
    """Minimal set of top-level host dirs to identity-bind so every path resolves."""
    dirs: list[Path] = []
    for p in paths:
        if not p:
            continue
        pp = Path(str(p)).expanduser()
        try:
            pp = pp.resolve()
        except Exception:
            pp = pp.absolute()
        dirs.append(pp if pp.is_dir() else pp.parent)
    dirs.append(Path("/tmp"))
    roots: list[Path] = []
    for d in sorted(set(dirs), key=lambda x: len(str(x))):
        if not any(d == r or _is_within(d, r) for r in roots):
            roots.append(d)
    return [str(r) for r in roots]


def _resolve_runtime(config: Mapping[str, Any]) -> str:
    if os.environ.get("BACTOLOMEW_NATIVE"):
        return "native"
    rt = str(config.get("runtime") or "").strip().lower()
    rt = "singularity" if rt == "apptainer" else rt
    return rt if rt in {*CONTAINER_RUNTIMES, "native"} else DEFAULT_RUNTIME


def configured_image_ref(config: Optional[Mapping[str, Any]] = None) -> str:
    """The image reference the config names (registry prefix applied, digest kept)."""
    config = config if isinstance(config, Mapping) else _load_app_config()
    container = config.get("container") if isinstance(config.get("container"), Mapping) else {}
    return _resolve_image_ref(container)


def image_has_repository_path(image: str) -> bool:
    """True when `image` names a repository path (`ghcr.io/x/y:t`, `localhost:5000/y`,
    `owner/y`): a reference a registry can serve as it is. A bare `name[:tag][@digest]`
    is a local name that `container.registry` prefixes."""
    return "/" in str(image).strip().partition("@")[0]  # a tag or digest never holds "/"


def podman_image_ref(image: str) -> str:
    """`image` as Podman should name it: a Docker Hub path without a registry host
    (`owner/name:tag`) gets `docker.io/`, since Podman's short-name policy can refuse it;
    bare local names (`bactolomew-tools:0.2.0`) and qualified references are unchanged."""
    ref = str(image).strip()
    first, sep, _rest = ref.partition("/")
    if not sep or "." in first or ":" in first or first == "localhost":
        return ref
    return f"docker.io/{ref}"


def _resolve_image_ref(container: Mapping[str, Any]) -> str:
    image = str(container.get("image") or TOOLS_IMAGE).strip()
    registry = str(container.get("registry") or "").strip().rstrip("/")
    if registry and not image_has_repository_path(image):
        return f"{registry}/{image}"
    return image


def make_container_context(
    config: Optional[Mapping[str, Any]] = None,
    *,
    output_dir: Optional[Any] = None,
    extra_paths: Iterable[Any] = (),
    workdir: Optional[Any] = None,
) -> ContainerContext:
    """Build a run/call-scoped context from config + the paths the work will touch."""
    config = dict(config) if isinstance(config, Mapping) else _load_app_config()
    runtime = _resolve_runtime(config)
    container = config.get("container") if isinstance(config.get("container"), Mapping) else {}

    paths: list[Any] = list(extra_paths)
    if output_dir:
        paths.append(output_dir)
    dbs = config.get("databases") if isinstance(config.get("databases"), Mapping) else {}
    paths.append(dbs.get("root_dir"))
    paths.append(config.get("bakta_db"))
    paths.append(configured_assembly_qc_settings().get("checkm2_database_path"))
    cb = configured_chewbbaca_settings()
    paths.extend([cb.get("schema_dir"), cb.get("catalog_path")])
    paths.append(configured_mobsuite_settings().get("database_directory"))

    # Writable HOME inside the container (under output_dir when we have one, else tmp).
    if output_dir:
        home = str(Path(str(output_dir)) / ".bactolomew_home")
    else:
        home = tempfile.mkdtemp(prefix="bactolomew-home-")
    try:
        Path(home).mkdir(parents=True, exist_ok=True)
    except Exception:
        home = tempfile.mkdtemp(prefix="bactolomew-home-")
    paths.append(home)

    bind_roots = compute_bind_roots([p for p in paths if p])
    wd = str(workdir) if workdir else (str(output_dir) if output_dir else "/tmp")
    image = _resolve_image_ref(container)
    engine = None
    if runtime in CONTAINER_RUNTIMES:
        # Which engine runs the containers (rootless, Podman, SELinux), probed once per process.
        from pipeline.ContainerRuntime import engine_info, resolve_pinned_image

        engine = engine_info(runtime)
        if runtime in OCI_RUNTIMES and "@sha256:" in image:
            # A pinned digest that the local image lacks: warn and run the local name[:tag].
            image = resolve_pinned_image(image, log=logging.getLogger(__name__).warning, engine=engine.exe)
    return ContainerContext(
        runtime=runtime,
        image=image,
        sif_path=str(container.get("sif_path") or ""),
        bind_roots=tuple(bind_roots),
        uid=os.getuid(),
        gid=os.getgid(),
        home_dir=home,
        default_workdir=wd,
        mem_limit_gb=configured_max_ram_gb(),
        shm_size=docker_shm_size(container.get("shm_size"), configured_max_ram_gb()),
        engine_exe=engine.exe if engine else "",
        podman=bool(engine and engine.podman),
        rootless=bool(engine and engine.rootless),
        selinux=bool(engine and engine.selinux),
    )


def _effective_context() -> ContainerContext:
    ctx = _CURRENT_CTX.get() or _ACTIVE_CTX
    return ctx if ctx is not None else make_container_context()


# Docker gives a container a private 64 MB /dev/shm, and Gubbins' multithreaded
# ancestral reconstruction (pyjar shared memory) dies with a bus error on a 150-genome
# core alignment there. Every tool container gets a bounded /dev/shm instead
# (config container.shm_size, docker syntax). A tmpfs size is a ceiling, not an
# allocation: pages count only when written. Shared memory counts against the
# container's --memory cgroup, so with resources.max_ram_gb set the size is capped at
# half of it. Singularity uses the host's /dev/shm.
DEFAULT_SHM_SIZE = "8g"
_SHM_UNITS = {"": 1, "b": 1, "k": 1 << 10, "m": 1 << 20, "g": 1 << 30}


def _shm_bytes(value: str) -> int:
    found = re.fullmatch(r"\s*(\d+)\s*([bkmg]?)b?\s*", value.lower())
    if not found:
        raise ValueError(f"container.shm_size: {value!r} is not a docker size such as 8g or 512m")
    return int(found.group(1)) * _SHM_UNITS[found.group(2)]


def docker_shm_size(value: Any, mem_limit_gb: int = 0) -> str:
    """The docker --shm-size value: container.shm_size (default 8g), at most half of
    the container memory cap when one is set."""
    text = str(value).strip() if value not in (None, "") else DEFAULT_SHM_SIZE
    size = _shm_bytes(text)
    if size <= 0:
        raise ValueError("container.shm_size: must be positive")
    if mem_limit_gb and mem_limit_gb > 0:
        cap = int(mem_limit_gb) * (1 << 30) // 2
        if size > cap:
            return f"{cap // (1 << 20)}m"
    return text.lower().replace(" ", "")


def _docker_shm_args(ctx: ContainerContext) -> list[str]:
    return ["--shm-size", ctx.shm_size] if ctx.shm_size else []


def singularity_bin() -> Optional[str]:
    return shutil.which("singularity") or shutil.which("apptainer")


def native_container_runtime() -> str:
    """The engine a native context runs its containers on (fastGEAR still runs in the
    tools image): docker when installed, else podman."""
    return "podman" if not shutil.which("docker") and shutil.which("podman") else "docker"


def _engine_traits(ctx: ContainerContext) -> tuple[bool, bool, bool]:
    """(podman, rootless, selinux) of the engine running `ctx`'s containers; a native
    context probes the one it falls back to (native_container_runtime)."""
    if ctx.runtime != "native":
        return ctx.podman, ctx.rootless, ctx.selinux
    from pipeline.ContainerRuntime import engine_info

    engine = engine_info(native_container_runtime())
    return engine.podman, engine.rootless, engine.selinux


def oci_exe(ctx: ContainerContext) -> str:
    """The executable that runs `ctx`'s containers; for a native context the probed
    docker or podman (native_container_runtime)."""
    if ctx.runtime in OCI_RUNTIMES:
        return ctx.engine_exe or ctx.runtime
    from pipeline.ContainerRuntime import engine_info

    return engine_info(native_container_runtime()).exe


def container_user_args(ctx: ContainerContext) -> list[str]:
    """The flags that make files written to the bind mounts belong to the host user.
    Rootful Docker or Podman: --user uid:gid. Rootless Podman: --userns=keep-id (the user
    keeps its uid inside; --user would name an unmapped subordinate uid). Rootless Docker:
    none (the container's root is the host user). With SELinux enforcing, label=disable
    keeps the identity mounts usable without relabeling the user's folders. A native
    context reads these from the engine it runs on (native_container_runtime)."""
    podman, rootless, selinux = _engine_traits(ctx)
    if podman and rootless:
        args = ["--userns=keep-id"]
    elif rootless:
        args = []
    else:
        args = ["--user", f"{ctx.uid}:{ctx.gid}"]
    if selinux:
        args += ["--security-opt", "label=disable"]
    return args


def _oci_run_cmd(ctx: ContainerContext, image: str, argv: list[str], *, env: Mapping[str, Any],
                 binds: Iterable[str], workdir: str) -> list[str]:
    cmd = [oci_exe(ctx), "run", "--rm", *container_user_args(ctx), *_docker_shm_args(ctx)]
    if ctx.mem_limit_gb and ctx.mem_limit_gb > 0:
        cmd += ["--memory", f"{ctx.mem_limit_gb}g"]
    for k, v in env.items():
        cmd += ["-e", f"{k}={v}"]
    for root in binds:
        cmd += ["-v", f"{root}:{root}"]
    if workdir:
        cmd += ["-w", workdir]
    return [*cmd, podman_image_ref(image) if _engine_traits(ctx)[0] else image, *argv]


def _singularity_exec_cmd(ctx: ContainerContext, argv: list[str], *, env: Mapping[str, Any],
                          binds: Iterable[str], workdir: str) -> list[str]:
    """`singularity exec` (or apptainer). HOME goes through --home: Apptainer refuses
    --env HOME=..."""
    env = dict(env)
    home = str(env.pop("HOME", "") or ctx.home_dir)
    cmd = [ctx.engine_exe or singularity_bin() or "singularity", "exec", "--cleanenv", "--home", home]
    for k, v in env.items():
        cmd += ["--env", f"{k}={v}"]
    for root in binds:
        cmd += ["--bind", root]
    if workdir:
        cmd += ["--pwd", workdir]
    return [*cmd, ctx.sif_path, *argv]


def _container_prefix(ctx: ContainerContext, env: str, exe: str, workdir: Optional[str]) -> list[str]:
    wd = str(workdir) if workdir else ctx.default_workdir
    inner = ["conda", "run", "--no-capture-output", "-n", env, exe]
    if ctx.runtime == "singularity":
        return _singularity_exec_cmd(ctx, inner, env={}, binds=ctx.bind_roots, workdir=wd)
    return _oci_run_cmd(ctx, ctx.image, inner, env={"HOME": ctx.home_dir}, binds=ctx.bind_roots, workdir=wd)


def container_image_cmd(
    ctx: ContainerContext,
    argv: Iterable[Any],
    *,
    workdir: Optional[Any] = None,
    extra_env: Optional[Mapping[str, Any]] = None,
    extra_binds: Iterable[Any] = (),
) -> list[str]:
    """Run an arbitrary argv directly inside the tools image (no conda-run wrapper).

    Used for tools that are not a `conda run -n env exe` (e.g. fastGEAR's MCR launcher).
    Container-only: a `native` context falls back to docker with the tools image.
    """
    wd = str(workdir) if workdir else ctx.default_workdir
    binds = list(ctx.bind_roots) + [str(b) for b in extra_binds if b]
    env = dict(extra_env or {})
    argv = [str(a) for a in argv]
    if ctx.runtime == "singularity":
        return _singularity_exec_cmd(ctx, argv, env=env, binds=binds, workdir=wd)
    return _oci_run_cmd(ctx, ctx.image or TOOLS_IMAGE, argv, env=env, binds=binds, workdir=wd)


def conda_tool_cmd(tool_name: str, executable: str | None = None) -> Optional[list[str]]:
    ctx = _effective_context()
    exe = executable or tool_name
    if ctx.runtime == "native":
        prefix = get_tool_env_prefix(tool_name)
        if not prefix:
            return None
        return ["conda", "run", "-p", prefix, exe]
    return _container_prefix(ctx, _tool_env_key(tool_name), exe, None)


def tool_cmd(tool_name: str, executable: str | None = None, *, workdir: str | None = None) -> list[str]:
    ctx = _effective_context()
    exe = executable or tool_name
    if ctx.runtime == "native":
        cmd = conda_tool_cmd(tool_name, executable=executable)
        return cmd if cmd else [exe]
    return _container_prefix(ctx, _tool_env_key(tool_name), exe, workdir)


def runtime_is_container(config: Optional[Mapping[str, Any]] = None) -> bool:
    config = config if isinstance(config, Mapping) else _load_app_config()
    return _resolve_runtime(config) in CONTAINER_RUNTIMES


def gui_tool_cmd(
    tool_name: str,
    executable: str | None = None,
    *,
    mount_paths: Iterable[Any] = (),
    workdir: str | None = None,
) -> list[str]:
    """Build a tool command for an ad-hoc GUI call (no run context).

    Constructs a fresh, one-off container context from config with `mount_paths`
    bound, without touching the module-global run context. Falls back to native.
    """
    ctx = make_container_context(extra_paths=mount_paths, workdir=workdir)
    exe = executable or tool_name
    if ctx.runtime == "native":
        prefix = get_tool_env_prefix(tool_name)
        return ["conda", "run", "-p", prefix, exe] if prefix else [exe]
    return _container_prefix(ctx, _tool_env_key(tool_name), exe, workdir)


def configured_bakta_db_path() -> str:
    config = _load_app_config()
    raw = config.get("bakta_db")
    if isinstance(raw, str) and raw.strip():
        return raw.strip()
    return ""


def configured_amrfinder_organism() -> str:
    """Optional AMRFinderPlus --organism override (config `amrfinder.organism`); '' when unset."""
    raw = _load_app_config().get("amrfinder")
    return str(raw.get("organism") or "").strip() if isinstance(raw, Mapping) else ""


def configured_database_dirs() -> dict[str, str]:
    """Current on-disk locations of the reference-DB sub-directories, read from
    bactolomew_config.json's `databases` block. Used to re-resolve absolute paths
    baked into an older run_config.json after the project or DB has moved/renamed."""
    config = _load_app_config()
    dbs = config.get("databases") if isinstance(config.get("databases"), Mapping) else {}
    return {
        "root_dir": str(dbs.get("root_dir") or "").strip(),
        "marker_sequences_dir": str(dbs.get("marker_sequences_dir") or "").strip(),
        "compound_sequences_dir": str(dbs.get("compound_sequences_dir") or "").strip(),
        "cgmlst_dir": str(dbs.get("cgmlst_dir") or "").strip(),
    }


def relocate_db_path(stored: Any, current_dir: Optional[str]) -> str:
    """Re-resolve a DB input path recorded in a run against the current config.

    If the stored path still exists, keep it. Otherwise, if `current_dir` holds an
    entry of the same basename, use that (handles a moved/renamed project or DB).
    Fall back to the stored value so validation still reports the original path when
    no better candidate exists.
    """
    text = str(stored).strip()
    if not text:
        return text
    p = Path(text).expanduser()
    if p.exists():
        return str(p)
    if current_dir:
        candidate = Path(current_dir) / p.name
        if candidate.exists():
            return str(candidate)
    return text


def configured_mobsuite_settings(raw_config: Mapping[str, Any] | None = None) -> dict[str, Any]:
    config = _load_app_config()
    app_raw = config.get("mobsuite")
    app_settings_raw = app_raw if isinstance(app_raw, Mapping) else {}
    settings_sources: list[Mapping[str, Any]] = [app_settings_raw]
    if isinstance(raw_config, Mapping):
        settings_sources.append(raw_config)

    settings: dict[str, Any] = {
        "database_directory": "",
        "extra_args": [],
    }

    for settings_raw in settings_sources:
        db_path = str(
            settings_raw.get("database_directory")
            or settings_raw.get("db_path")
            or settings_raw.get("database")
            or settings_raw.get("db")
            or ""
        ).strip()
        if db_path:
            settings["database_directory"] = db_path

        extra_args = _normalize_shell_tokens(settings_raw.get("extra_args"))
        if extra_args:
            settings["extra_args"] = extra_args

    return settings


def configured_fastgear_settings(raw_config: Mapping[str, Any] | None = None) -> dict[str, Any]:
    if raw_config is None:
        config = _load_app_config()
        raw = config.get("fastgear")
    else:
        raw = raw_config

    settings_raw = raw if isinstance(raw, Mapping) else {}
    settings: dict[str, Any] = {
        "docker_bin": DEFAULT_FASTGEAR_DOCKER_BIN,
        "image": DEFAULT_FASTGEAR_IMAGE,
        "launcher": DEFAULT_FASTGEAR_LAUNCHER,
        "mcr_root": DEFAULT_FASTGEAR_MCR_ROOT,
        "specs_file": DEFAULT_FASTGEAR_SPECS_FILE,
        "workdir": DEFAULT_FASTGEAR_WORKDIR,
    }

    image = str(settings_raw.get("image") or "").strip()
    if image and image != LEGACY_FASTGEAR_IMAGE:
        settings["image"] = image

    docker_bin = str(settings_raw.get("docker_bin") or DEFAULT_FASTGEAR_DOCKER_BIN).strip()
    if docker_bin:
        settings["docker_bin"] = docker_bin

    workdir = str(settings_raw.get("workdir") or DEFAULT_FASTGEAR_WORKDIR).strip()
    if workdir.startswith("/"):
        settings["workdir"] = workdir

    launcher = str(settings_raw.get("launcher") or "").strip()
    if launcher.startswith("/"):
        settings["launcher"] = launcher

    mcr_root = str(settings_raw.get("mcr_root") or "").strip()
    if mcr_root.startswith("/"):
        settings["mcr_root"] = mcr_root

    specs_file = str(settings_raw.get("specs_file") or "").strip()
    if specs_file.startswith("/"):
        settings["specs_file"] = specs_file

    extra_args = _normalize_shell_tokens(settings_raw.get("extra_args"))
    if extra_args:
        settings["extra_args"] = extra_args

    return settings


def configured_fastgear_image(config: Optional[Mapping[str, Any]] = None) -> tuple[str, str]:
    """(image, archive) of the fastGEAR base image a local build of the tools image copies:
    container.fastgear_image (else fastgear.image, else the pinned default) and
    container.fastgear_archive, a `docker save` archive loaded when the image is missing."""
    config = config if isinstance(config, Mapping) else _load_app_config()
    container = config.get("container") if isinstance(config.get("container"), Mapping) else {}
    fastgear = config.get("fastgear") if isinstance(config.get("fastgear"), Mapping) else None
    image = str(container.get("fastgear_image") or "").strip() or configured_fastgear_settings(fastgear or {})["image"]
    archive = str(container.get("fastgear_archive") or "").strip()
    return image, str(Path(archive).expanduser()) if archive else ""


# Phylogeny settings (config section `phylo`): one set of automatic rules for any
# cohort, each overridable. See PHYLOGENY_REWORK_2026-10-10.md section 3.8.
PHYLO_DEFAULTS: dict[str, Any] = {
    "core_alignment": "full",            # Panaroo alignment for every core consumer: full | filtered
    "order_reference": "",               # empty = auto (most contiguous genome); sample name or FASTA path
    "ordering_min_placed": 0.95,         # auto ordering genome must place this share of core genes
    "outgroup": "",                      # sample name or comma list; empty = midpoint rooting
    "seed": 1,                           # Gubbins and IQ-TREE seed
    "gubbins_iterations": 20,
    "gubbins_mode": "auto",              # auto (per fastGEAR lineage when >= 2) | single | per_lineage
    "gubbins_min_lineage_size": 5,       # lineages with fewer genomes are not run (left unmasked)
    "gubbins_converge_method": "recombination",  # stop when the recombination calls repeat
    "gubbins_filter_percentage": 25.0,
    "gubbins_tree_builder": "raxml",
    "gubbins_min_genomes": 4,
    "gubbins_diversity_warn": 0.01,      # mean pairwise core p-distance
    "gubbins_max_diversity_warn": 0.02,  # max pairwise core p-distance
    "gubbins_longest_branch_warn": 0.5,  # share of clonal-frame SNPs on one branch
    "iqtree_model": "MFP",
    "iqtree_bootstraps": 1000,
    "iqtree_bnni": True,
    "iqtree_fast_above_genomes": 1000,
    # Second branch length set of core_recomb: its topology refit on the unmasked alignment.
    "recombination_branch_lengths": True,
}
_PHYLO_CHOICES: dict[str, tuple[str, ...]] = {
    "core_alignment": ("full", "filtered"),
    "gubbins_converge_method": ("weighted_robinson_foulds", "robinson_foulds", "recombination"),
    "gubbins_tree_builder": ("raxml", "raxmlng", "iqtree", "iqtree-fast", "fasttree", "hybrid", "rapidnj", "veryfasttree"),
    "gubbins_mode": ("auto", "single", "per_lineage"),
}


def configured_phylo_settings(raw_config: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Effective phylogeny settings: defaults, then config `phylo`, then raw_config.

    Values of the wrong type or outside their choices raise ValueError, so a typo in
    the config stops the run at its start instead of silently using a default.
    """
    config = _load_app_config()
    sources: list[Mapping[str, Any]] = []
    if isinstance(config.get("phylo"), Mapping):
        sources.append(config["phylo"])
    if isinstance(raw_config, Mapping):
        sources.append(raw_config)
    settings: dict[str, Any] = dict(PHYLO_DEFAULTS)
    for source in sources:
        for key, value in source.items():
            if key not in PHYLO_DEFAULTS or value is None:
                continue
            default = PHYLO_DEFAULTS[key]
            try:
                if isinstance(default, bool):
                    if isinstance(value, str):
                        value = value.strip().lower() in {"1", "true", "yes", "on"}
                    else:
                        value = bool(value)
                elif isinstance(default, int):
                    value = int(value)
                elif isinstance(default, float):
                    value = float(value)
                elif isinstance(value, (list, tuple)):
                    # A list of names (phylo.outgroup) is kept as a comma-joined string.
                    value = ",".join(str(item).strip() for item in value if str(item).strip())
                else:
                    value = str(value).strip()
            except (TypeError, ValueError):
                raise ValueError(f"phylo.{key}: invalid value {value!r}") from None
            choices = _PHYLO_CHOICES.get(key)
            if choices and value not in choices:
                raise ValueError(f"phylo.{key}: {value!r} is not one of {', '.join(choices)}")
            settings[key] = value
    if not 0.0 < settings["ordering_min_placed"] <= 1.0:
        raise ValueError("phylo.ordering_min_placed: must be in (0, 1]")
    if settings["gubbins_iterations"] < 1:
        raise ValueError("phylo.gubbins_iterations: must be at least 1")
    if settings["gubbins_min_genomes"] < 3:
        # Gubbins itself refuses fewer than three sequences.
        raise ValueError("phylo.gubbins_min_genomes: must be at least 3")
    if settings["gubbins_min_lineage_size"] < 3:
        raise ValueError("phylo.gubbins_min_lineage_size: must be at least 3")
    if settings["iqtree_bootstraps"] < 1000:
        # IQ-TREE refuses UFBoot with fewer than 1000 replicates.
        raise ValueError("phylo.iqtree_bootstraps: IQ-TREE ultrafast bootstrap needs at least 1000")
    return settings


def phylo_outgroup_names(settings: Mapping[str, Any]) -> list[str]:
    return [name for name in dict.fromkeys(
        part.strip() for part in str(settings.get("outgroup") or "").split(",")
    ) if name]


def configured_chewbbaca_settings(raw_config: Mapping[str, Any] | None = None) -> dict[str, Any]:
    config = _load_app_config()
    app_raw = config.get("chewbbaca")
    app_settings_raw = app_raw if isinstance(app_raw, Mapping) else {}
    settings_sources: list[Mapping[str, Any]] = [app_settings_raw]
    if isinstance(raw_config, Mapping):
        settings_sources.append(raw_config)

    settings: dict[str, Any] = {
        "schema_dir": "",
        "training_file": "",
        "genes_list": "",
        "loci_list": "",
        "catalog_path": "",
        "available_schemas": {},
        "auto_select_from_mlst": False,
    }
    auto_select_set = False

    for settings_raw in settings_sources:
        schema_dir = str(
            settings_raw.get("schema_dir")
            or settings_raw.get("schema_directory")
            or ""
        ).strip()
        if schema_dir:
            settings["schema_dir"] = schema_dir

        training_file = str(
            settings_raw.get("training_file")
            or settings_raw.get("ptf")
            or ""
        ).strip()
        if training_file:
            settings["training_file"] = training_file

        genes_list_raw = str(settings_raw.get("genes_list") or "").strip()
        if genes_list_raw:
            if Path(genes_list_raw).name == ".genes_list":
                settings["genes_list"] = genes_list_raw
            else:
                settings["loci_list"] = genes_list_raw

        loci_list = str(
            settings_raw.get("loci_list")
            or settings_raw.get("subset_loci_list")
            or settings_raw.get("gl")
            or ""
        ).strip()
        if loci_list:
            settings["loci_list"] = loci_list

        catalog_path = str(
            settings_raw.get("catalog_path")
            or settings_raw.get("scheme_catalog_path")
            or settings_raw.get("cgmlst_scheme_catalog_path")
            or ""
        ).strip()
        if catalog_path:
            settings["catalog_path"] = catalog_path

        available_schemas_raw = settings_raw.get("available_schemas")
        if isinstance(available_schemas_raw, Mapping):
            merged_schemas = dict(settings.get("available_schemas") or {})
            for label, path_value in available_schemas_raw.items():
                label_text = str(label).strip()
                path_text = str(path_value).strip()
                if label_text and path_text:
                    merged_schemas[label_text] = path_text
            settings["available_schemas"] = merged_schemas

        if "auto_select_from_mlst" in settings_raw:
            settings["auto_select_from_mlst"] = bool(settings_raw.get("auto_select_from_mlst"))
            auto_select_set = True

        extra_args = _normalize_shell_tokens(settings_raw.get("extra_args"))
        if extra_args:
            settings["extra_args"] = extra_args

    if not auto_select_set:
        settings["auto_select_from_mlst"] = bool(
            settings["catalog_path"] or settings["available_schemas"] or settings["schema_dir"]
        )

    return settings


def configured_assembly_qc_settings(raw_config: Mapping[str, Any] | None = None) -> dict[str, Any]:
    config = _load_app_config()
    app_raw = config.get("assembly_qc")
    app_settings_raw = app_raw if isinstance(app_raw, Mapping) else {}
    settings_sources: list[Mapping[str, Any]] = [app_settings_raw]
    if isinstance(raw_config, Mapping):
        settings_sources.append(raw_config)

    settings: dict[str, Any] = {
        "quast_min_contig": 500,
        "quast_extra_args": [],
        "checkm2_database_path": "",
        "checkm2_lowmem": False,
        "checkm2_extra_args": [],
    }

    for settings_raw in settings_sources:
        min_contig = settings_raw.get("quast_min_contig", settings_raw.get("min_contig"))
        if min_contig not in (None, ""):
            try:
                settings["quast_min_contig"] = max(0, int(min_contig))
            except Exception:
                pass

        quast_extra_args = _normalize_shell_tokens(settings_raw.get("quast_extra_args"))
        if quast_extra_args:
            settings["quast_extra_args"] = quast_extra_args

        db_path = str(
            settings_raw.get("checkm2_database_path")
            or settings_raw.get("checkm2_db")
            or settings_raw.get("database_path")
            or ""
        ).strip()
        if db_path:
            settings["checkm2_database_path"] = db_path

        if "checkm2_lowmem" in settings_raw:
            settings["checkm2_lowmem"] = bool(settings_raw.get("checkm2_lowmem"))
        elif "lowmem" in settings_raw:
            settings["checkm2_lowmem"] = bool(settings_raw.get("lowmem"))

        checkm2_extra_args = _normalize_shell_tokens(settings_raw.get("checkm2_extra_args"))
        if checkm2_extra_args:
            settings["checkm2_extra_args"] = checkm2_extra_args

    if not str(settings.get("checkm2_database_path") or "").strip():
        env_path = get_tool_env_prefix("checkm2")
        if env_path:
            default_db = Path(env_path) / CHECKM2_DB_RELATIVE_PATH
            if default_db.is_file():
                settings["checkm2_database_path"] = str(default_db)

    # Auto-enable CheckM2 low-memory mode when the run's RAM (the configured cap, or the
    # machine / cgroup / SLURM memory detected) is tight. CheckM2 runs DIAMOND over the
    # uniref100.KO database with a large block size, which takes a big share of a 32 GB
    # node with many threads, and other steps of the run share that memory; --lowmem
    # (smaller DIAMOND blocks: slower, a few GB) avoids the out-of-memory kill that a
    # cluster slot or a container gets instead of swapping. The bound keeps the earlier
    # max_ram_gb rule (below 32 GB); a 32 GB machine reports a little less and gets it.
    ram_gb = effective_ram_gb()
    if ram_gb and ram_gb < CHECKM2_LOWMEM_BELOW_GB:
        settings["checkm2_lowmem"] = True

    return settings


__all__ = [
    "REQUIRED_TOOL_ENVS",
    "TOOLS_IMAGE",
    "CONTAINER_RUNTIMES",
    "OCI_RUNTIMES",
    "configured_image_ref",
    "podman_image_ref",
    "container_user_args",
    "oci_exe",
    "native_container_runtime",
    "singularity_bin",
    "configured_fastgear_image",
    "DEFAULT_RUNTIME",
    "ContainerContext",
    "save_app_config",
    "configured_resources",
    "effective_core_budget",
    "configured_max_ram_gb",
    "effective_ram_gb",
    "resource_budget_text",
    "docker_shm_size",
    "compute_bind_roots",
    "make_container_context",
    "set_container_context",
    "reset_container_context",
    "current_container_context",
    "DEV_BAKTA_ENV_PREFIX",
    "DEV_BAKTA_DB_PATH",
    "DEFAULT_FASTGEAR_DOCKER_BIN",
    "DEFAULT_FASTGEAR_IMAGE",
    "DEFAULT_FASTGEAR_LAUNCHER",
    "DEFAULT_FASTGEAR_MCR_ROOT",
    "DEFAULT_FASTGEAR_SPECS_FILE",
    "DEFAULT_FASTGEAR_WORKDIR",
    "project_root",
    "get_tool_env_prefix",
    "conda_tool_cmd",
    "tool_cmd",
    "configured_amrfinder_organism",
    "configured_bakta_db_path",
    "configured_database_dirs",
    "relocate_db_path",
    "configured_mobsuite_settings",
    "configured_fastgear_settings",
    "configured_chewbbaca_settings",
    "configured_assembly_qc_settings",
    "configured_phylo_settings",
    "phylo_outgroup_names",
    "PHYLO_DEFAULTS",
]
