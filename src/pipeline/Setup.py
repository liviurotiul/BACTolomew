"""First-time setup without a window (Qt free): the tool image or environments, the Bakta
(full) and CheckM2 databases, the DB/ folder and the config, run as `SetupJob` steps.
The GUI dialog (`annotation_explorer_polars.first_time_setup`) and `bactolomew_main.py setup`
both drive it."""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import threading
from dataclasses import dataclass, replace
from datetime import datetime
from pathlib import Path
from typing import Callable, Iterable, Optional

from pipeline.ChewBBACAWorkflow import chewbbaca_schema_issues
from pipeline.DbCatalog import (
    COMPOUND_DB_CATALOG_PATH,
    MARKER_DB_CATALOG_PATH,
    load_db_catalog,
)
from pipeline.ToolEnv import (
    DEV_BAKTA_DB_PATH,
    REQUIRED_TOOL_ENVS,
    TOOLS_IMAGE,
    configured_assembly_qc_settings,
    configured_chewbbaca_settings,
    configured_fastgear_image,
    configured_fastgear_settings,
    project_root,
)


LogFn = Callable[[str], None]


DB_DIRNAME = "DB"
DB_CATALOGS_DIRNAME = "catalogs"
CHEWBBACA_GENE_LIST = ".genes_list"
CGMLST_SCHEME_CATALOG_FILENAME = "scheme_species_catalog.json"
# Our species catalog of the cgMLST schemes (species and MLST aliases per scheme), copied
# into DB/cgmlst by setup when that folder has none.
CGMLST_SCHEME_CATALOG_SOURCE = Path(__file__).resolve().parent / "data" / CGMLST_SCHEME_CATALOG_FILENAME
YAML_DIRNAME = "yamls"
ENVS_DIRNAME = "envs"
BAKTA_DB_DIRNAME = "bakta_db"
CHECKM2_DB_DIRNAME = "checkm2_db"
CHECKM2_DB_RELATIVE_PATH = Path("CheckM2_database") / "uniref100.KO.1.dmnd"
DOCKER_CLIENT_SENTINEL = ".docker_client"
FASTGEAR_IMAGE_SENTINEL = ".fastgear_docker_image"
TOOLS_IMAGE_SENTINEL = ".bactolomew_tools_image"
MANAGED_TOOL_ENVS = REQUIRED_TOOL_ENVS

def _normalized(path: Path) -> Path:
    path = Path(path).expanduser()
    try:
        return path.resolve()
    except Exception:
        return path


def _same_path(left: Path, right: Path) -> bool:
    return _normalized(left) == _normalized(right)


def _is_within(path: Path, root: Path) -> bool:
    try:
        _normalized(path).relative_to(_normalized(root))
        return True
    except Exception:
        return False


def _project_db_root(root: Path) -> Path:
    return root / DB_DIRNAME


def _project_marker_dir(root: Path) -> Path:
    return _project_db_root(root) / "seq" / "markers"


def _project_compound_dir(root: Path) -> Path:
    return _project_db_root(root) / "seq" / "compound"


def _project_metadata_dir(root: Path) -> Path:
    return _project_db_root(root) / "metadata"


def _project_cgmlst_dir(root: Path) -> Path:
    return _project_db_root(root) / "cgmlst"


def _project_cgmlst_catalog_path(root: Path) -> Path:
    return _project_cgmlst_dir(root) / CGMLST_SCHEME_CATALOG_FILENAME


def _project_catalog_dir(root: Path) -> Path:
    return _project_db_root(root) / DB_CATALOGS_DIRNAME


def _project_marker_catalog_path(root: Path) -> Path:
    return _project_catalog_dir(root) / MARKER_DB_CATALOG_PATH.name


def _project_compound_catalog_path(root: Path) -> Path:
    return _project_catalog_dir(root) / COMPOUND_DB_CATALOG_PATH.name


def _project_db_layout_ready(root: Path) -> bool:
    required = [
        _project_marker_dir(root),
        _project_compound_dir(root),
        _project_metadata_dir(root),
        _project_catalog_dir(root),
        _project_marker_catalog_path(root),
        _project_compound_catalog_path(root),
    ]
    return all(path.exists() for path in required)


def ensure_db_layout(root: Path, log: Optional[LogFn] = None) -> Path:
    """Make `DB/` as setup leaves it: empty marker and compound catalogs (no database is
    shipped; users add their own in Settings), their sequence folders, the metadata
    folder and the cgMLST folder with our scheme species catalog. What exists is kept."""
    for path in (_project_marker_dir(root), _project_compound_dir(root), _project_metadata_dir(root),
                 _project_catalog_dir(root), _project_cgmlst_dir(root)):
        _ensure_dir_preserve(path, path.name)
    for catalog in (_project_marker_catalog_path(root), _project_compound_catalog_path(root)):
        if not catalog.exists():
            catalog.write_text("{}\n", encoding="utf-8")
            if log is not None:
                log(f"Created the empty catalog {catalog}")
    scheme_catalog = _project_cgmlst_catalog_path(root)
    if not scheme_catalog.exists() and CGMLST_SCHEME_CATALOG_SOURCE.is_file():
        shutil.copyfile(CGMLST_SCHEME_CATALOG_SOURCE, scheme_catalog)
        if log is not None:
            log(f"Copied the cgMLST scheme catalog to {scheme_catalog}")
    return _project_db_root(root)


def _is_conda_env(path: Path) -> bool:
    return path.exists() and path.is_dir() and (path / "conda-meta").is_dir()


def _yaml_pinned_versions(yaml_path: Path) -> dict[str, str]:
    pins: dict[str, str] = {}
    if not yaml_path.exists() or not yaml_path.is_file():
        return pins

    in_dependencies = False
    for raw_line in yaml_path.read_text(encoding="utf-8").splitlines():
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped == "dependencies:":
            in_dependencies = True
            continue
        if not in_dependencies or not stripped.startswith("- "):
            continue

        spec = stripped[2:].strip()
        if not spec or spec.endswith(":"):
            continue
        if "=" not in spec:
            continue

        name, remainder = spec.split("=", 1)
        version = remainder.split("=", 1)[0].strip()
        name = name.strip()
        if name and version:
            pins[name] = version
    return pins


def _conda_meta_versions(env_path: Path) -> dict[str, str]:
    conda_meta = env_path / "conda-meta"
    versions: dict[str, str] = {}
    if not conda_meta.exists() or not conda_meta.is_dir():
        return versions

    for meta_path in sorted(conda_meta.glob("*.json")):
        try:
            raw = json.loads(meta_path.read_text(encoding="utf-8"))
        except Exception:
            continue
        if not isinstance(raw, dict):
            continue
        name = str(raw.get("name") or "").strip()
        version = str(raw.get("version") or "").strip()
        if name and version:
            versions[name] = version
    return versions


def _conda_env_matches_yaml(env_path: Path, yaml_path: Path) -> bool:
    expected = _yaml_pinned_versions(yaml_path)
    if not expected:
        return True
    observed = _conda_meta_versions(env_path)
    if not observed:
        return False
    for name, version in expected.items():
        if observed.get(name) != version:
            return False
    return True


def expected_tool_env_paths(root: Path) -> dict[str, Path]:
    env_root = root / ENVS_DIRNAME
    return {name: env_root / name for name in MANAGED_TOOL_ENVS}


def _is_bakta_db_dir(path: Path) -> bool:
    return path.exists() and path.is_dir() and (path / "version.json").is_file()


def resolve_bakta_db_path(base: Path) -> Optional[Path]:
    base = _normalized(base)
    candidates = [base / "db", base]
    if base.exists() and base.is_dir():
        candidates.extend(p for p in sorted(base.iterdir()) if p.is_dir())
    for cand in candidates:
        if _is_bakta_db_dir(cand):
            return cand
    return None


def _docker_client_available(docker_bin: str) -> bool:
    return shutil.which(docker_bin) is not None


def _docker_image_available(docker_bin: str, image: str) -> bool:
    try:
        res = subprocess.run(
            [docker_bin, "image", "inspect", image],
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        return False
    return res.returncode == 0


def _ensure_fastgear_image(settings: dict[str, object], *, log: LogFn, archive: str = "") -> dict[str, object]:
    """Make the fastGEAR image present: loaded from `archive` (container.fastgear_archive)
    when missing, else pulled."""
    docker_bin = str(settings.get("docker_bin") or "docker").strip() or "docker"
    image = str(settings.get("image") or "").strip()
    if not image:
        raise RuntimeError("fastGEAR Docker image is not configured.")
    if not _docker_client_available(docker_bin):
        raise RuntimeError(f"Docker client not found in PATH: {docker_bin}")
    if _docker_image_available(docker_bin, image):
        return settings

    from pipeline import ContainerRuntime as CR

    CR.ensure_fastgear_image(image, archive, engine=docker_bin, log=log)
    return settings


def _catalog_targets_dir(catalog_path: Path, expected_dir: Path) -> bool:
    """The catalog exists and every database it lists lies in `expected_dir`. The
    databases themselves are optional: a run searches the installed ones."""
    if not catalog_path.exists() or not catalog_path.is_file():
        return False
    if not expected_dir.exists() or not expected_dir.is_dir():
        return False

    try:
        raw = load_db_catalog(catalog_path)
    except Exception:
        return False
    if not isinstance(raw, dict):
        return False

    expected_dir = _normalized(expected_dir)
    for info in raw.values():
        if not isinstance(info, dict):
            return False
        db_path_raw = info.get("db_path") or info.get("path") or info.get("db-path")
        if not db_path_raw:
            return False
        if _normalized(Path(str(db_path_raw))).parent != expected_dir:
            return False
    return True


def _db_config_outdated(root: Path, config: dict) -> bool:
    if not _project_db_layout_ready(root):
        return False

    db_root = _project_db_root(root)
    marker_dir = _project_marker_dir(root)
    compound_dir = _project_compound_dir(root)
    metadata_dir = _project_metadata_dir(root)
    cgmlst_dir = _project_cgmlst_dir(root)
    marker_catalog = _project_marker_catalog_path(root)
    compound_catalog = _project_compound_catalog_path(root)

    db_cfg = config.get("databases") if isinstance(config.get("databases"), dict) else {}
    expected_db_paths = {
        "root_dir": db_root,
        "marker_sequences_dir": marker_dir,
        "compound_sequences_dir": compound_dir,
        "marker_catalog_path": marker_catalog,
        "compound_catalog_path": compound_catalog,
        "metadata_dir": metadata_dir,
        "cgmlst_dir": cgmlst_dir,
    }
    for key, expected in expected_db_paths.items():
        raw = str(db_cfg.get(key) or "").strip()
        if not raw:
            return True
        if not _same_path(Path(raw), expected):
            return True

    available_schemas = _discover_chewbbaca_schemas(cgmlst_dir)
    existing_chewbbaca = config.get("chewbbaca") if isinstance(config.get("chewbbaca"), dict) else {}
    chewbbaca_settings = configured_chewbbaca_settings(existing_chewbbaca)
    schema_dir = str(chewbbaca_settings.get("schema_dir") or "").strip()
    genes_list = str(chewbbaca_settings.get("genes_list") or "").strip()
    preferred_schema = _preferred_schema_path(available_schemas, existing_chewbbaca)
    if preferred_schema is not None:
        if not schema_dir or not _same_path(Path(schema_dir), preferred_schema):
            return True
        preferred_genes = preferred_schema / CHEWBBACA_GENE_LIST
        if preferred_genes.is_file():
            if not genes_list or not _same_path(Path(genes_list), preferred_genes):
                return True

    return False


def collect_setup_missing_paths(root: Path, config: dict, *, dev_mode: bool = False) -> list[Path]:
    catalogs_dir = _project_catalog_dir(root)
    marker_catalog = _project_marker_catalog_path(root)
    compound_catalog = _project_compound_catalog_path(root)
    db_root = _project_db_root(root)
    marker_dir = _project_marker_dir(root)
    compound_dir = _project_compound_dir(root)
    metadata_dir = _project_metadata_dir(root)

    missing: list[Path] = []
    if not _catalog_targets_dir(marker_catalog, marker_dir):
        missing.append(marker_catalog)
    if not _catalog_targets_dir(compound_catalog, compound_dir):
        missing.append(compound_catalog)

    required_dirs = [
        catalogs_dir,
        marker_dir,
        compound_dir,
        metadata_dir,
        root / "runs",
        root / "cache",
    ]
    missing.extend(p for p in required_dirs if not (p.exists() and p.is_dir()))

    config_refresh_needed = _db_config_outdated(root, config)
    if config_refresh_needed:
        missing.append(root / "bactolomew_config.json")

    db_cfg = config.get("databases") if isinstance(config.get("databases"), dict) else {}
    configured_root_raw = str(db_cfg.get("root_dir") or "").strip()
    if not config_refresh_needed and configured_root_raw and not _same_path(Path(configured_root_raw), db_root):
        missing.append(db_root)

    expected_db_paths = {
        "marker_sequences_dir": marker_dir,
        "compound_sequences_dir": compound_dir,
        "marker_catalog_path": marker_catalog,
        "compound_catalog_path": compound_catalog,
    }
    legacy_db_root = root / "databases"
    if not config_refresh_needed:
        for key, expected in expected_db_paths.items():
            raw = str(db_cfg.get(key) or "").strip()
            if not raw:
                continue
            configured = Path(raw)
            if not _same_path(configured, expected):
                missing.append(expected)
            if _is_within(configured, legacy_db_root):
                missing.append(db_root)

    if not dev_mode:
        bakta_db_raw = str(config.get("bakta_db") or "").strip()
        if not bakta_db_raw:
            missing.append(root / BAKTA_DB_DIRNAME)
        else:
            bakta_db = Path(bakta_db_raw).expanduser()
            if resolve_bakta_db_path(bakta_db) is None:
                missing.append(bakta_db)

    from pipeline.ToolEnv import runtime_is_container

    if runtime_is_container(config):
        # Container mode: one tools image replaces the per-tool conda envs + the
        # separate fastGEAR image. Flag the image (or the runtime) as missing.
        from pipeline import ContainerRuntime as _CR

        container = config.get("container") if isinstance(config.get("container"), dict) else {}
        image = str(container.get("image") or "").strip() or TOOLS_IMAGE
        runtime = str(config.get("runtime") or "docker").strip().lower()
        runtime = "singularity" if runtime == "apptainer" else runtime
        if runtime == "singularity":
            if not _CR.sif_present(str(container.get("sif_path") or "")):
                missing.append(root / TOOLS_IMAGE_SENTINEL)
        elif not _CR.runtime_status(runtime).ok:
            missing.append(root / DOCKER_CLIENT_SENTINEL)
        elif not _CR.docker_image_present(image, _CR.engine_exe(runtime)):
            missing.append(root / TOOLS_IMAGE_SENTINEL)
    else:
        # Native (dev) mode: per-tool conda envs under envs/.
        env_cfg = config.get("envs") if isinstance(config.get("envs"), dict) else {}
        expected = expected_tool_env_paths(root)
        for name in expected.keys():
            default_path = expected[name]
            raw = env_cfg.get(name)
            if not isinstance(raw, str) or not raw.strip():
                missing.append(default_path)
                continue
            env_path = Path(raw).expanduser()
            if not _is_conda_env(env_path):
                missing.append(env_path)
                continue
            if name == "bakta":
                yaml_path = root / YAML_DIRNAME / "bakta.yaml"
                if not _conda_env_matches_yaml(env_path, yaml_path):
                    missing.append(env_path)

        fastgear_settings = configured_fastgear_settings(config.get("fastgear"))
        docker_bin = str(fastgear_settings.get("docker_bin") or "docker")
        image = str(fastgear_settings.get("image") or "").strip()
        if not _docker_client_available(docker_bin):
            missing.append(root / DOCKER_CLIENT_SENTINEL)
        elif image and not _docker_image_available(docker_bin, image):
            missing.append(root / FASTGEAR_IMAGE_SENTINEL)

    chewbbaca_settings = configured_chewbbaca_settings(config.get("chewbbaca"))
    for key in ("schema_dir", "training_file", "genes_list", "loci_list"):
        value = str(chewbbaca_settings.get(key) or "").strip()
        if not value:
            continue
        path = Path(value).expanduser()
        if key == "schema_dir":
            if config_refresh_needed:
                continue
            if chewbbaca_schema_issues(path):
                missing.append(path)
        elif not path.exists() or not path.is_file():
            missing.append(path)

    dedup: list[Path] = []
    seen: set[str] = set()
    for p in missing:
        key = str(p)
        if key in seen:
            continue
        seen.add(key)
        dedup.append(p)
    return dedup


def _discover_chewbbaca_schemas(cgmlst_root: Path) -> dict[str, Path]:
    if not cgmlst_root.exists() or not cgmlst_root.is_dir():
        return {}

    schemas: dict[str, Path] = {}
    for child in sorted(cgmlst_root.iterdir(), key=lambda path: path.name.lower()):
        if not child.is_dir():
            continue
        if chewbbaca_schema_issues(child):
            continue
        schemas[child.name] = _normalized(child)
    return schemas


def _preferred_schema_path(
    available_schemas: dict[str, Path],
    existing_chewbbaca: Optional[dict],
) -> Optional[Path]:
    if not available_schemas:
        return None

    existing_settings = configured_chewbbaca_settings(existing_chewbbaca if isinstance(existing_chewbbaca, dict) else None)
    existing_schema = str(existing_settings.get("schema_dir") or "").strip()
    if existing_schema:
        existing_name = Path(existing_schema).name.lower()
        for schema_path in available_schemas.values():
            if schema_path.name.lower() == existing_name:
                return schema_path

    return next(iter(available_schemas.values()))


def _ensure_dir_preserve(path: Path, label: str) -> None:
    if path.exists():
        if not path.is_dir():
            raise RuntimeError(f"{label} path exists but is not a directory: {path}")
        return
    path.mkdir(parents=True, exist_ok=True)


def _run_cmd(cmd: list[str], *, log: LogFn, env: Optional[dict[str, str]] = None) -> None:
    from pipeline import ContainerRuntime as CR

    rc = CR._stream(cmd, log, env=env)
    if rc != 0:
        raise RuntimeError(f"{' '.join(cmd[:4])} failed (exit {rc})")


def _create_tool_envs(root: Path, *, log: LogFn, tools: Optional[Iterable[str]] = None) -> dict[str, Path]:
    if shutil.which("conda") is None:
        raise RuntimeError("conda was not found in PATH. Install conda to create tool environments.")

    yaml_root = root / YAML_DIRNAME
    env_root = root / ENVS_DIRNAME
    env_root.mkdir(parents=True, exist_ok=True)

    tool_names = list(tools) if tools is not None else list(MANAGED_TOOL_ENVS)
    env_paths: dict[str, Path] = {}
    for tool in tool_names:
        yaml_path = yaml_root / f"{tool}.yaml"
        if not yaml_path.exists() or not yaml_path.is_file():
            raise RuntimeError(f"Missing environment YAML: {yaml_path}")

        env_path = env_root / tool
        if _is_conda_env(env_path):
            if tool == "bakta" and not _conda_env_matches_yaml(env_path, yaml_path):
                cmd = ["conda", "env", "update", "-p", str(env_path), "-f", str(yaml_path), "--prune"]
                _run_cmd(cmd, log=log)
                if not _conda_env_matches_yaml(env_path, yaml_path):
                    raise RuntimeError(f"Failed syncing conda env to YAML: {env_path}")
            env_paths[tool] = env_path.resolve()
            continue
        if env_path.exists() and not _is_conda_env(env_path):
            raise RuntimeError(f"Env path exists but is not a conda env: {env_path}")

        cmd = ["conda", "env", "create", "-p", str(env_path), "-f", str(yaml_path)]
        _run_cmd(cmd, log=log)
        if not _is_conda_env(env_path):
            raise RuntimeError(f"Failed creating conda env: {env_path}")
        env_paths[tool] = env_path.resolve()

    return env_paths


def _tool_envs_ready(root: Path) -> bool:
    for tool, env_path in expected_tool_env_paths(root).items():
        if not _is_conda_env(env_path):
            return False
        if tool == "bakta" and not _conda_env_matches_yaml(env_path, root / YAML_DIRNAME / "bakta.yaml"):
            return False
    return True


def _ensure_bakta_db(bakta_env: Optional[Path], db_root: Path, *, log: LogFn) -> Path:
    db_root = db_root.resolve()
    db_root.mkdir(parents=True, exist_ok=True)

    existing = resolve_bakta_db_path(db_root)
    if existing is not None:
        return existing

    download = ["bakta_db", "download", "--output", str(db_root), "--type", "full"]
    if bakta_env is not None:
        download = ["conda", "run", "--no-capture-output", "-p", str(bakta_env), *download]
    _run_cmd(download, log=log)

    resolved = resolve_bakta_db_path(db_root)
    if resolved is None:
        raise RuntimeError(f"Bakta DB download finished but no DB was found under {db_root}")
    return resolved


def _native_checkm2_db_path(assembly_qc_env: Path) -> Path:
    return Path(assembly_qc_env).expanduser() / CHECKM2_DB_DIRNAME / CHECKM2_DB_RELATIVE_PATH


def _ensure_checkm2_db(assembly_qc_env: Optional[Path], *, log: LogFn) -> Path:
    if assembly_qc_env is None:
        raise RuntimeError("assembly_qc conda environment was not created.")
    env_path = Path(assembly_qc_env).expanduser().resolve()
    if not _is_conda_env(env_path):
        raise RuntimeError(f"assembly_qc path is not a conda env: {env_path}")

    db_root = env_path / CHECKM2_DB_DIRNAME
    db_path = _native_checkm2_db_path(env_path)
    if db_path.is_file():
        return db_path.resolve()

    db_root.mkdir(parents=True, exist_ok=True)
    mpl_cache = env_path / "var" / "matplotlib"
    mpl_cache.mkdir(parents=True, exist_ok=True)
    run_env = dict(os.environ)
    run_env["MPLCONFIGDIR"] = str(mpl_cache)

    _run_cmd(
        [
            "conda", "run", "--no-capture-output", "-p", str(env_path),
            "checkm2", "database", "--download", "--path", str(db_root),
        ],
        log=log,
        env=run_env,
    )
    if db_path.is_file():
        return db_path.resolve()

    candidates = sorted(db_root.rglob("*.dmnd"))
    if candidates:
        return candidates[0].resolve()
    raise RuntimeError(f"CheckM2 database download completed but no DIAMOND DB was found under: {db_root}")


def _read_config(root: Path) -> dict:
    config_path = root / "bactolomew_config.json"
    if not config_path.exists():
        return {}
    try:
        raw = json.loads(config_path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return raw if isinstance(raw, dict) else {}


def _write_setup_config(
    root: Path,
    existing_config: dict,
    *,
    bakta_db: Path,
    checkm2_db: Path,
    envs: dict[str, str],
    dev_mode: bool,
    log: LogFn,
    container: bool = False,
) -> None:
    """bactolomew_config.json after setup: the DB/ paths, the cgMLST schemas, the tool
    databases and envs; every other key of the existing config is kept."""
    db_root = _project_db_root(root)
    metadata_dir = _project_metadata_dir(root)
    cgmlst_dir = _project_cgmlst_dir(root)

    available_schemas = _discover_chewbbaca_schemas(cgmlst_dir)
    existing_chewbbaca = existing_config.get("chewbbaca") if isinstance(existing_config.get("chewbbaca"), dict) else {}
    selected_schema = _preferred_schema_path(available_schemas, existing_chewbbaca)

    db_config: dict[str, object] = {
        "root_dir": str(db_root),
        "marker_sequences_dir": str(_project_marker_dir(root)),
        "compound_sequences_dir": str(_project_compound_dir(root)),
        "marker_catalog_path": str(_project_marker_catalog_path(root)),
        "compound_catalog_path": str(_project_compound_catalog_path(root)),
        "metadata_dir": str(metadata_dir),
        "cgmlst_dir": str(cgmlst_dir),
    }
    cgmlst_catalog_path = _project_cgmlst_catalog_path(root)
    if cgmlst_catalog_path.is_file():
        db_config["cgmlst_scheme_catalog_path"] = str(cgmlst_catalog_path)

    chewbbaca_config = configured_chewbbaca_settings(existing_chewbbaca)
    chewbbaca_config["available_schemas"] = {label: str(path) for label, path in available_schemas.items()}
    if cgmlst_catalog_path.is_file():
        chewbbaca_config["catalog_path"] = str(cgmlst_catalog_path)
        chewbbaca_config["auto_select_from_mlst"] = True
    else:
        chewbbaca_config["catalog_path"] = ""
        chewbbaca_config["auto_select_from_mlst"] = False
    if selected_schema is not None:
        chewbbaca_config["schema_dir"] = str(selected_schema)
        genes_list = selected_schema / CHEWBBACA_GENE_LIST
        chewbbaca_config["genes_list"] = str(genes_list) if genes_list.is_file() else ""
    else:
        chewbbaca_config["schema_dir"] = ""
        chewbbaca_config["genes_list"] = ""
    chewbbaca_config["training_file"] = ""

    config = {
        "bakta_db": str(bakta_db),
        "envs_root": str((root / ENVS_DIRNAME).resolve()),
        "envs": envs,
        "dev_mode": bool(dev_mode),
        "databases": db_config,
        "created_at": datetime.now().astimezone().isoformat(),
    }

    existing_fastgear = existing_config.get("fastgear")
    fastgear_settings = configured_fastgear_settings(
        existing_fastgear if isinstance(existing_fastgear, dict) else None
    )
    existing_assembly_qc = existing_config.get("assembly_qc")
    # In a container runtime fastGEAR runs from the tools image: no separate image to fetch.
    config["fastgear"] = fastgear_settings if container else _ensure_fastgear_image(
        fastgear_settings, log=log, archive=configured_fastgear_image(existing_config)[1])
    config["chewbbaca"] = chewbbaca_config
    config["mobsuite"] = {}
    assembly_qc_config = configured_assembly_qc_settings(
        existing_assembly_qc if isinstance(existing_assembly_qc, dict) else None
    )
    assembly_qc_config["checkm2_database_path"] = str(checkm2_db)
    config["assembly_qc"] = assembly_qc_config

    for key, value in existing_config.items():  # runtime, container, resources, ...
        config.setdefault(key, value)

    (root / "bactolomew_config.json").write_text(json.dumps(config, indent=2, sort_keys=True), encoding="utf-8")


# --------------------------------------------------------------------------- #
# First-time setup as steps (no Qt): what is installed already is kept, the rest is
# installed in order. The setup window (annotation_explorer_polars.first_time_setup)
# shows the steps and runs them on a thread; `bactolomew_main.py setup` prints them.
# --------------------------------------------------------------------------- #
# Under the Bakta step's name until it is installed: the full database's free space
# (DbDownload.BAKTA_FULL_NEED_BYTES: 115.7 GB while it unpacks plus a 5 GB margin) and time
# (a 39.4 GB download, then unpacking 76 GB).
BAKTA_HINT = "Needs about 121 GB free, takes 1 to 2 hours"

STEP_WAITING = "waiting"
STEP_RUNNING = "running"
STEP_DONE = "done"
STEP_FAILED = "failed"


class SetupCancelled(Exception):
    pass


class SetupStepError(RuntimeError):
    """A step failed; `short` is the one line its row shows (the log has the rest)."""

    def __init__(self, short: str, detail: str = "") -> None:
        super().__init__(detail or short)
        self.short = short


@dataclass
class SetupStep:
    key: str
    label: str
    size: str = ""        # what it downloads, e.g. "1.7 GB"
    hint: str = ""        # under its name until it is done, e.g. the free space it needs
    weight: float = 1.0   # its share of the window's progress bar
    state: str = STEP_WAITING
    text: str = "waiting"
    note: str = ""        # the failure, on its row
    detail: str = ""      # tooltip: where it is installed
    fraction: float = 0.0


def _size_text(n: int) -> str:
    if n <= 0:
        return ""
    gb = n / 1e9
    return f"{gb:.1f} GB" if gb >= 1 else f"{n / 1e6:.0f} MB"


def _short_error(exc: BaseException) -> str:
    from pipeline.DbDownload import NotEnoughSpace

    if isinstance(exc, SetupStepError):
        return exc.short
    if isinstance(exc, NotEnoughSpace):
        return (f"Needs about {exc.need_bytes / 1e9:.0f} GB free to unpack "
                f"({exc.free_bytes / 1e9:.1f} GB free)")
    text = (str(exc).strip().splitlines() or [type(exc).__name__])[0]
    failed = re.match(r"^\S.*download.* failed \(exit (-?\d+)\)$", text)
    if failed:
        if int(failed.group(1)) in (125, 126, 127):  # docker run itself, not the tool
            return f"The download could not start in the tool image (exit {failed.group(1)})"
        return "Download failed. Check the connection and retry."
    if len(text) > 110:
        cut = text.find(". ")
        text = text[:cut + 1] if 0 < cut <= 110 else text[:107].rstrip() + "..."
    return text


_PERCENT_RE = re.compile(r"(\d{1,3})(?:\.\d+)?%\s*\|")
_BUILD_STEP_RE = re.compile(r"^#\d+ \[[^\]]*?(\d+)/(\d+)\]")
_LAYER_RE = re.compile(r"^([0-9a-f]{6,64}): (Pulling fs layer|Waiting|Downloading|Verifying Checksum|"
                       r"Download complete|Extracting|Pull complete|Already exists)")
# Phases a download step goes through after the transfer: (text in a line, status, progress).
_DB_PHASES = (
    ("checksum", "verifying", 0.78),
    ("extract", "unpacking", 0.8),
    ("amrfinderplus database", "updating", 0.86),
    ("checking the bundled", "verifying", 0.94),
    ("verifying", "verifying", 0.9),
    ("hashing", "verifying", 0.95),
)
_DOWNLOAD_SHARE = 0.75  # the transfer's share of a database step


def setup_fraction(steps: Iterable[SetupStep]) -> float:
    """The whole setup's progress (0..1): each step by its weight. A failed step counts
    nothing (Retry starts it over)."""
    steps = list(steps)
    total = sum(s.weight for s in steps) or 1.0
    return sum(s.weight * (1.0 if s.state == STEP_DONE else 0.0 if s.state == STEP_FAILED else s.fraction)
               for s in steps) / total


class _LineProgress:
    """A step's output lines read as (fraction 0..1, short status): tqdm percentages of
    a download, docker pull layers, docker build steps and the phases after a download."""

    def __init__(self, key: str) -> None:
        self.key = key
        self.fraction = 0.0
        self.text = ""
        self._layers: set[str] = set()
        self._pulled: set[str] = set()

    def _set(self, fraction: float, text: str) -> None:
        self.fraction = max(self.fraction, min(0.99, fraction))
        self.text = text

    def feed(self, line: str) -> None:
        text = line.strip()
        if self.key == "image":
            layer = _LAYER_RE.match(text)
            if layer:
                self._layers.add(layer.group(1))
                if layer.group(2) in ("Pull complete", "Already exists"):
                    self._pulled.add(layer.group(1))
                done = len(self._pulled) / max(1, len(self._layers))
                self._set(0.95 * done, f"{int(100 * done)}%")
                return
            build = _BUILD_STEP_RE.match(text)
            if build and int(build.group(2)):
                done = int(build.group(1)) / int(build.group(2))
                self._set(0.95 * done, f"{int(100 * done)}%")
                return
            if text.startswith("Pulling "):
                self._set(0.0, "downloading")
            elif text.startswith("Building image") or text.startswith("container.build_on_missing"):
                self._set(0.0, "building")
            return
        percent = _PERCENT_RE.search(text)
        if percent and self.fraction < _DOWNLOAD_SHARE:
            pct = min(100, int(percent.group(1)))
            self._set(_DOWNLOAD_SHARE * pct / 100, f"{pct}%")
            return
        low = text.lower()
        for needle, status, fraction in _DB_PHASES:
            if needle in low and fraction >= self.fraction:
                self._set(fraction, status)
                return


class SetupJob:
    """The steps of first-time setup, run in order on any thread.

    Tool image (or the tool environments in native mode), the full Bakta database (with
    its AMRFinderPlus database; its free space is checked before anything starts), the
    CheckM2 database and a quick check that makes the `DB/` folder and writes the config.
    Whatever is installed already is found first and shown done; nothing is deleted or
    replaced (a light Bakta database installed earlier is kept and used).
    `on_step(step)` gets a copy of a step at each change, `log(line)` every output line.
    """

    def __init__(
        self,
        root: Optional[Path] = None,
        *,
        dev_mode: bool = False,
        on_step: Optional[Callable[[SetupStep], None]] = None,
        log: Optional[LogFn] = None,
        config: Optional[dict] = None,
    ) -> None:
        from pipeline.DbDownload import SPECS
        from pipeline.ToolEnv import runtime_is_container

        self.root = Path(root) if root is not None else project_root()
        self.dev_mode = bool(dev_mode)
        self.config = dict(config) if config is not None else _read_config(self.root)
        self.container = runtime_is_container(self.config) and not self.dev_mode
        self._on_step = on_step or (lambda _step: None)
        self._log_fn = log or print
        self._cancel = threading.Event()
        self.thread_id: Optional[int] = None
        self._found: dict[str, object] = {}  # step key -> what it installed or found

        self.steps = [
            SetupStep("image", "Tool image" if self.container else "Tool environments", weight=2),
            SetupStep("bakta", "Bakta database", weight=4,
                      size="" if self.dev_mode else _size_text(SPECS["bakta-full"].download_bytes),
                      hint="" if self.dev_mode else BAKTA_HINT),
            SetupStep("checkm2", "CheckM2 database", size=_size_text(SPECS["checkm2"].download_bytes), weight=4),
            SetupStep("check", "Quick check", weight=0.5),
        ]

    # ---- state ----
    def cancel(self) -> None:
        from pipeline import ContainerRuntime as CR

        self._cancel.set()
        if self.thread_id is not None:
            CR.stop_stream(self.thread_id)

    def rearm(self) -> None:
        """Clear a stop before running again (on the starting thread, so a stop that
        comes while the worker starts is not lost)."""
        self._cancel.clear()

    @property
    def cancelled(self) -> bool:
        return self._cancel.is_set()

    def _update(self, step: SetupStep, **changes: object) -> None:
        for name, value in changes.items():
            setattr(step, name, value)
        self._on_step(replace(step))

    def _log(self, line: str) -> None:
        if self._cancel.is_set():
            raise SetupCancelled()
        self._log_fn(line)

    # ---- run ----
    def run(self) -> str:
        """Run every step not done yet; "done", "failed" (the failed step says why) or
        "cancelled". Run again to retry from the failed step (after `rearm()` when it
        was stopped)."""
        for step in self.steps:
            if step.state != STEP_DONE:
                self._update(step, state=STEP_WAITING, text="waiting", note="", fraction=0.0)
        for step in self.steps:
            if step.state != STEP_DONE and step.key != "check":
                try:
                    found = getattr(self, f"_find_{step.key}")()
                except Exception:  # noqa: BLE001 - a failed lookup means "install it"
                    found = None
                if found is not None:
                    self._found[step.key] = found
                    self._update(step, state=STEP_DONE, text="done", detail=f"Already installed: {found}")
        short = self._space_problem()
        if short is not None:  # before anything downloads (the full Bakta database needs about 121 GB)
            step, row_text, log_text = short
            self._log_fn(f"ERROR: {log_text}")
            self._update(step, state=STEP_FAILED, text="failed", note=row_text)
            return "failed"
        for step in self.steps:
            if step.state == STEP_DONE:
                continue
            if self._cancel.is_set():
                return "cancelled"
            self._update(step, state=STEP_RUNNING, text="checking" if step.key == "check" else "starting")
            progress = _LineProgress(step.key)

            def step_log(line: str, step: SetupStep = step, progress: _LineProgress = progress) -> None:
                self._log(line)
                before = (progress.fraction, progress.text)
                progress.feed(line)
                if progress.text and (progress.fraction, progress.text) != before:
                    self._update(step, fraction=progress.fraction, text=progress.text)

            self._log(f"== {step.label} ==")
            try:
                found = getattr(self, f"_run_{step.key}")(step_log)
            except BaseException as exc:  # noqa: BLE001 - shown on the step's row
                if self._cancel.is_set() or isinstance(exc, SetupCancelled):
                    self._update(step, state=STEP_WAITING, text="stopped", fraction=0.0)
                    return "cancelled"
                if not isinstance(exc, Exception):
                    raise
                self._log_fn(f"ERROR: {exc}")
                self._update(step, state=STEP_FAILED, text="failed", note=_short_error(exc))
                return "failed"
            self._found[step.key] = found
            self._update(step, state=STEP_DONE, text="done", fraction=1.0,
                         detail=str(found) if found else step.detail)
        return "done"

    def status(self) -> list[tuple[SetupStep, object]]:
        """Each step with what is installed for it (None: to install; the check always runs)."""
        out: list[tuple[SetupStep, object]] = []
        for step in self.steps:
            try:
                found = None if step.key == "check" else getattr(self, f"_find_{step.key}")()
            except Exception:  # noqa: BLE001 - a failed lookup means "install it"
                found = None
            out.append((step, found))
        return out

    def dry_run(self) -> list[str]:
        """What `run()` would do, as lines, changing nothing: each step found installed
        (and where) or to be installed, with the download commands of the databases."""
        from pipeline import DbDownload

        lines: list[str] = []
        for step, found in self.status():
            if found is not None:
                lines.append(f"{step.label}: already installed ({found})")
                continue
            size = f", {step.size}" if step.size else ""
            hint = f" ({step.hint[:1].lower()}{step.hint[1:]})" if step.hint else ""
            lines.append(f"{step.label}: {'run' if step.key == 'check' else 'install'}{size}{hint}")
            if step.key == "bakta":
                folder = self._space_folder()
                lines.append(f"  free space in {folder}: {DbDownload._free_bytes(folder) / 1e9:.0f} GB")
            elif step.key == "check" and not _project_db_layout_ready(self.root):
                lines.append(f"  make {_project_db_root(self.root)} (empty marker and compound catalogs, "
                             "the cgMLST scheme catalog)")
            if not self.container:
                continue
            if step.key == "image":
                runtime, image, registry, sif, build = self._container_settings()
                target = sif if runtime == "singularity" else image
                how = "pull" + (" (or build when the pull fails)" if build else "")
                lines.append(f"  {how} {target} with {runtime}" + (f" from {registry}" if registry else ""))
            elif step.key in ("bakta", "checkm2"):
                quiet = lambda _line: None  # noqa: E731 - the commands are listed below
                record = (DbDownload.install_bakta(self._bakta_target(), config=self.config, log=quiet, dry_run=True)
                          if step.key == "bakta" else
                          DbDownload.install_checkm2(self.root / CHECKM2_DB_DIRNAME, config=self.config, log=quiet,
                                                     dry_run=True))
                lines.extend("  $ " + " ".join(map(str, cmd)) for cmd in record.get("commands") or [])
        return lines

    # ---- tool image / tool environments ----
    def _container_settings(self) -> tuple[str, str, str, str, bool]:
        container = self.config.get("container") if isinstance(self.config.get("container"), dict) else {}
        image = str(container.get("image") or "").strip() or TOOLS_IMAGE
        runtime = str(self.config.get("runtime") or "docker").strip().lower()
        runtime = "singularity" if runtime == "apptainer" else runtime
        registry = str(container.get("registry") or "").strip()
        sif = str(container.get("sif_path") or "").strip()
        return runtime, image, registry, sif, container.get("build_on_missing", True) is not False

    def _find_image(self) -> Optional[str]:
        from pipeline import ContainerRuntime as CR

        if not self.container:
            return str(self.root / ENVS_DIRNAME) if _tool_envs_ready(self.root) else None
        runtime, image, _registry, sif, _build = self._container_settings()
        if runtime == "singularity":
            return sif if CR.sif_present(sif) else None
        if not CR.runtime_status(runtime).ok:
            return None
        return image if CR.docker_image_present(image, CR.engine_exe(runtime)) else None

    def _run_image(self, log: LogFn) -> str:
        from pipeline import ContainerRuntime as CR

        if not self.container:
            _create_tool_envs(self.root, log=log)
            return str(self.root / ENVS_DIRNAME)
        runtime, image, registry, sif, build_on_missing = self._container_settings()
        status = CR.runtime_status(runtime)
        if not status.ok:
            raise SetupStepError(status.message, "\n".join(filter(None, (status.message, status.details))))
        fastgear_image, fastgear_archive = configured_fastgear_image(self.config)
        CR.ensure_tools_image(runtime, image=image, registry=registry, sif_path=sif, repo_root=self.root,
                              log=log, build_on_missing=build_on_missing,
                              fastgear_image=fastgear_image, fastgear_archive=fastgear_archive)
        return sif if runtime == "singularity" else image

    def _env_path(self, tool: str) -> Path:
        return expected_tool_env_paths(self.root)[tool]

    # ---- Bakta ----
    def _bakta_target(self) -> Path:
        return self.root / BAKTA_DB_DIRNAME

    @staticmethod
    def _existing(folder: Path) -> Path:
        """The nearest existing folder (the one whose disk receives ``folder``)."""
        while not folder.exists() and folder != folder.parent:
            folder = folder.parent
        return folder

    def _space_folder(self) -> Path:
        """The folder whose disk receives the Bakta database (its nearest existing one)."""
        return self._existing(self._bakta_target().parent)

    def _space_problem(self) -> Optional[tuple[SetupStep, str, str]]:
        """(step, row text, log text) when a disk has less free space than the steps still
        to install there need together: the full Bakta database while it unpacks
        (115.7 GB), CheckM2 (4.8 GB) and, with Singularity, the tool image's .sif (about
        15 GB), plus one 5 GB margin per disk. The step is the largest one on that disk.
        None when there is room (or in dev mode)."""
        from pipeline import DbDownload

        if self.dev_mode:
            return None
        by_key = {step.key: step for step in self.steps}
        pending: list[tuple[str, Path, int, str]] = []
        if "bakta" in by_key and by_key["bakta"].state != STEP_DONE:
            pending.append(("bakta", self._bakta_target().parent, DbDownload.BAKTA_FULL_PEAK_BYTES, "Bakta database"))
        if self.container and "checkm2" in by_key and by_key["checkm2"].state != STEP_DONE:
            pending.append(("checkm2", self.root, DbDownload.CHECKM2_PEAK_BYTES, "CheckM2"))
        runtime, _image, _registry, sif, _build = self._container_settings()
        if self.container and runtime == "singularity" and sif and "image" in by_key and by_key["image"].state != STEP_DONE:
            pending.append(("image", Path(sif).expanduser().parent, DbDownload.SIF_NEED_BYTES, "tool image .sif"))
        disks: dict[int, list[tuple[str, Path, int, str]]] = {}
        for item in pending:
            folder = self._existing(item[1])
            disks.setdefault(folder.stat().st_dev, []).append((item[0], folder, item[2], item[3]))
        for items in disks.values():
            need = sum(item[2] for item in items) + DbDownload.SPACE_MARGIN_BYTES
            folder = items[0][1]
            free = DbDownload._free_bytes(folder)
            if free >= need:
                continue
            detail = ", ".join(f"{label} {size / 1e9:.1f} GB" for _key, _f, size, label in items)
            detail += f", {DbDownload.SPACE_MARGIN_BYTES / 1e9:.0f} GB margin"
            exc = DbDownload.NotEnoughSpace(folder, free, need, detail)
            step = by_key[max(items, key=lambda item: item[2])[0]]
            return step, _short_error(exc), str(exc)
        return None

    def _find_bakta(self) -> Optional[Path]:
        """The configured Bakta database (full, or a light one installed earlier), else
        one in the default folder."""
        if self.dev_mode:
            return Path(DEV_BAKTA_DB_PATH)
        configured = str(self.config.get("bakta_db") or "").strip()
        if configured:
            found = resolve_bakta_db_path(Path(configured).expanduser())
            if found is not None:
                return found
        return resolve_bakta_db_path(self._bakta_target())

    def _run_bakta(self, log: LogFn) -> Path:
        from pipeline import DbDownload

        target = self._bakta_target()
        if not self.container:
            return _ensure_bakta_db(self._env_path("bakta"), target, log=log)
        record = DbDownload.install_bakta(target, config=self.config, log=log)
        found = resolve_bakta_db_path(Path(record["path"]))
        if found is None:
            raise SetupStepError(f"No Bakta database found in {target}")
        return found

    # ---- CheckM2 ----
    def _find_checkm2(self) -> Optional[Path]:
        if not self.container:
            path = _native_checkm2_db_path(self._env_path("assembly_qc"))
            return path if path.is_file() else None
        qc = self.config.get("assembly_qc") if isinstance(self.config.get("assembly_qc"), dict) else {}
        configured = str(qc.get("checkm2_database_path") or "").strip()
        if configured and Path(configured).expanduser().is_file():
            return Path(configured).expanduser()
        default = self.root / CHECKM2_DB_DIRNAME / CHECKM2_DB_RELATIVE_PATH
        return default if default.is_file() else None

    def _run_checkm2(self, log: LogFn) -> Path:
        from pipeline import DbDownload

        if not self.container:
            return _ensure_checkm2_db(self._env_path("assembly_qc"), log=log)
        record = DbDownload.install_checkm2(self.root / CHECKM2_DB_DIRNAME, config=self.config, log=log)
        path = Path(record["path"])
        if not path.is_file():
            raise SetupStepError("No CheckM2 database found", f"No CheckM2 DIAMOND database at {path}")
        return path

    # ---- quick check, then the config ----
    def _run_check(self, log: LogFn) -> str:
        from pipeline import DbDownload

        root = self.root
        _ensure_dir_preserve(root / "runs", "runs")
        _ensure_dir_preserve(root / "cache", "cache")
        ensure_db_layout(root, log)
        bakta_db = Path(self._found["bakta"])  # type: ignore[arg-type]
        checkm2_db = Path(self._found["checkm2"])  # type: ignore[arg-type]
        if not self.dev_mode and not (bakta_db / "version.json").is_file():
            raise SetupStepError("The Bakta database is incomplete", f"No version.json in {bakta_db}")
        if not checkm2_db.is_file() or checkm2_db.stat().st_size == 0:
            raise SetupStepError("The CheckM2 database is incomplete", f"{checkm2_db} is missing or empty")
        if self.container:
            try:
                DbDownload.verify_tools(config=self.config, log=log)
            except DbDownload.DbInstallError as exc:
                tool = str(exc).split(" ", 1)[0]
                raise SetupStepError(f"{tool} does not start in the tool image", str(exc)) from exc
            envs: dict[str, str] = {}
        else:
            envs = {name: str(path.resolve()) for name, path in expected_tool_env_paths(root).items()}
        log("Writing bactolomew_config.json")
        _write_setup_config(
            root,
            _read_config(root),
            bakta_db=bakta_db,
            checkm2_db=checkm2_db,
            envs=envs,
            dev_mode=self.dev_mode,
            log=log,
            container=self.container,
        )
        return str(root / "bactolomew_config.json")
