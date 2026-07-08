from __future__ import annotations

import contextvars
import json
import os
import shlex
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, List, Mapping, Optional, Tuple

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
DEFAULT_FASTGEAR_IMAGE = "pgcbioinfo/fastgear:latest"
DEFAULT_FASTGEAR_LAUNCHER = "/opt/run_fastGEAR.sh"
DEFAULT_FASTGEAR_MCR_ROOT = "/opt/v901"
DEFAULT_FASTGEAR_SPECS_FILE = "/opt/fG_input_specs.txt"
DEFAULT_FASTGEAR_WORKDIR = "/work"
CHECKM2_DB_RELATIVE_PATH = Path("checkm2_db") / "CheckM2_database" / "uniref100.KO.1.dmnd"


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
# Containerized execution (Docker / Singularity) of the pipeline tools.
#
# Every tool call funnels through tool_cmd()/conda_tool_cmd(). Instead of
# `conda run -p <env>` on the host, we wrap `conda run -n <env> <exe>` in a
# `docker run` / `singularity exec` against one image that carries every tool,
# bind-mounting the host dirs a run touches at IDENTICAL paths so the existing
# absolute-path command args resolve unchanged inside the container.
# ---------------------------------------------------------------------------

TOOLS_IMAGE = "bactolomew-tools:0.1.0"
DEFAULT_RUNTIME = "native"  # backward-compatible default when config has no `runtime`


@dataclass(frozen=True)
class ContainerContext:
    runtime: str                        # "docker" | "singularity" | "native"
    image: str = ""                     # docker image ref
    sif_path: str = ""                  # singularity .sif path
    bind_roots: Tuple[str, ...] = ()    # identity mount roots
    uid: int = 0
    gid: int = 0
    home_dir: str = "/tmp"              # writable HOME inside the container
    default_workdir: str = "/tmp"


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
    return rt if rt in {"docker", "singularity", "native"} else DEFAULT_RUNTIME


def _resolve_image_ref(container: Mapping[str, Any]) -> str:
    image = str(container.get("image") or TOOLS_IMAGE).strip()
    registry = str(container.get("registry") or "").strip().rstrip("/")
    if registry and "/" not in image.split(":", 1)[0]:
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
    return ContainerContext(
        runtime=runtime,
        image=_resolve_image_ref(container),
        sif_path=str(container.get("sif_path") or ""),
        bind_roots=tuple(bind_roots),
        uid=os.getuid(),
        gid=os.getgid(),
        home_dir=home,
        default_workdir=wd,
    )


def _effective_context() -> ContainerContext:
    ctx = _CURRENT_CTX.get() or _ACTIVE_CTX
    return ctx if ctx is not None else make_container_context()


def _container_prefix(ctx: ContainerContext, env: str, exe: str, workdir: Optional[str]) -> list[str]:
    wd = str(workdir) if workdir else ctx.default_workdir
    inner = ["conda", "run", "--no-capture-output", "-n", env, exe]
    if ctx.runtime == "docker":
        cmd = [
            "docker", "run", "--rm",
            "--user", f"{ctx.uid}:{ctx.gid}",
            "-e", f"HOME={ctx.home_dir}",
        ]
        for root in ctx.bind_roots:
            cmd += ["-v", f"{root}:{root}"]
        if wd:
            cmd += ["-w", wd]
        cmd += [ctx.image, *inner]
        return cmd
    # singularity (wired but untested this iteration)
    cmd = ["singularity", "exec", "--cleanenv", "--home", ctx.home_dir]
    for root in ctx.bind_roots:
        cmd += ["--bind", root]
    if wd:
        cmd += ["--pwd", wd]
    cmd += [ctx.sif_path, *inner]
    return cmd


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
    runtime = ctx.runtime if ctx.runtime in ("docker", "singularity") else "docker"
    wd = str(workdir) if workdir else ctx.default_workdir
    binds = list(ctx.bind_roots) + [str(b) for b in extra_binds if b]
    env = dict(extra_env or {})
    argv = [str(a) for a in argv]
    if runtime == "docker":
        cmd = ["docker", "run", "--rm", "--user", f"{ctx.uid}:{ctx.gid}"]
        for k, v in env.items():
            cmd += ["-e", f"{k}={v}"]
        for root in binds:
            cmd += ["-v", f"{root}:{root}"]
        if wd:
            cmd += ["-w", wd]
        cmd += [ctx.image or TOOLS_IMAGE, *argv]
        return cmd
    cmd = ["singularity", "exec", "--cleanenv"]
    for k, v in env.items():
        cmd += ["--env", f"{k}={v}"]
    for root in binds:
        cmd += ["--bind", root]
    if wd:
        cmd += ["--pwd", wd]
    cmd += [ctx.sif_path, *argv]
    return cmd


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
    return _resolve_runtime(config) in ("docker", "singularity")


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
    return DEV_BAKTA_DB_PATH


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
    if image:
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

    return settings


__all__ = [
    "REQUIRED_TOOL_ENVS",
    "TOOLS_IMAGE",
    "DEFAULT_RUNTIME",
    "ContainerContext",
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
    "configured_bakta_db_path",
    "configured_mobsuite_settings",
    "configured_fastgear_settings",
    "configured_chewbbaca_settings",
    "configured_assembly_qc_settings",
]
