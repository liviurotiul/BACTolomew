from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable, List, Optional

from PyQt6 import QtCore, QtWidgets

from pipeline.ChewBBACAWorkflow import chewbbaca_schema_issues
from pipeline.DbCatalog import (
    COMPOUND_DB_CATALOG_PATH,
    MARKER_DB_CATALOG_PATH,
    known_marker_resistance_metadata,
    load_db_catalog,
    normalize_resistance_metadata,
)
from pipeline.ToolEnv import (
    DEV_BAKTA_DB_PATH,
    DEV_BAKTA_ENV_PREFIX,
    REQUIRED_TOOL_ENVS,
    configured_assembly_qc_settings,
    configured_chewbbaca_settings,
    configured_fastgear_settings,
)
from pipeline.Utils import cat_specs_no_other

from .cache import project_root


FASTA_EXTS = (
    ".fa", ".fna", ".fasta", ".faa",
    ".fa.gz", ".fna.gz", ".fasta.gz", ".faa.gz",
)

DB_DIRNAME = "DB"
DB_LATEST_DIRNAME = "DB_latest"
DB_CATALOGS_DIRNAME = "catalogs"
DB_DEFAULT_SOURCE_ROOT = Path("/Data/BACTolomewDB")
CHEWBBACA_GENE_LIST = ".genes_list"
CGMLST_SCHEME_CATALOG_FILENAME = "scheme_species_catalog.json"
ARCHIVE_SUFFIXES = (
    ".tar.xz",
    ".txz",
    ".tar.gz",
    ".tgz",
    ".tar.bz2",
    ".tbz2",
    ".tar",
    ".zip",
)

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


def _ensure_fastgear_image(settings: dict[str, object]) -> dict[str, object]:
    docker_bin = str(settings.get("docker_bin") or "docker").strip() or "docker"
    image = str(settings.get("image") or "").strip()
    if not image:
        raise RuntimeError("fastGEAR Docker image is not configured.")
    if not _docker_client_available(docker_bin):
        raise RuntimeError(f"Docker client not found in PATH: {docker_bin}")
    if _docker_image_available(docker_bin, image):
        return settings

    subprocess.run([docker_bin, "pull", image], check=True)
    if not _docker_image_available(docker_bin, image):
        raise RuntimeError(f"Docker image pull finished but image is still unavailable: {image}")
    return settings


def _is_fasta(path: Path) -> bool:
    name = path.name.lower()
    return any(name.endswith(ext) for ext in FASTA_EXTS)


def _strip_fasta_suffix(name: str) -> str:
    lname = name.lower()
    for ext in sorted(FASTA_EXTS, key=len, reverse=True):
        if lname.endswith(ext):
            return name[: -len(ext)]
    return Path(name).stem


def _archive_suffix(path: Path) -> str:
    lower_name = path.name.lower()
    for suffix in ARCHIVE_SUFFIXES:
        if lower_name.endswith(suffix):
            return suffix
    return ""


def _is_supported_archive(path: Path) -> bool:
    return path.exists() and path.is_file() and bool(_archive_suffix(path))


def _default_db_archive_path() -> Optional[Path]:
    candidates = [
        DB_DEFAULT_SOURCE_ROOT / "DBs" / "DB_v100.tar.xz",
        DB_DEFAULT_SOURCE_ROOT / "DB_v100.tar.xz",
    ]
    for candidate in candidates:
        if _is_supported_archive(candidate):
            return _normalized(candidate)
    return None


def _catalog_targets_dir(catalog_path: Path, expected_dir: Path) -> bool:
    if not catalog_path.exists() or not catalog_path.is_file():
        return False
    if not expected_dir.exists() or not expected_dir.is_dir():
        return False

    try:
        raw = load_db_catalog(catalog_path)
    except Exception:
        return False
    if not isinstance(raw, dict) or not raw:
        return False

    expected_dir = _normalized(expected_dir)
    found = 0
    for info in raw.values():
        if not isinstance(info, dict):
            return False
        db_path_raw = info.get("db_path") or info.get("path") or info.get("db-path")
        if not db_path_raw:
            return False
        db_path = _normalized(Path(str(db_path_raw)))
        if db_path.parent != expected_dir:
            return False
        if not db_path.exists() or not db_path.is_file():
            return False
        found += 1
    return found > 0


def _existing_marker_catalog_ready(root: Path) -> bool:
    return _catalog_targets_dir(_project_marker_catalog_path(root), _project_marker_dir(root))


def _existing_compound_catalog_ready(root: Path) -> bool:
    return _catalog_targets_dir(_project_compound_catalog_path(root), _project_compound_dir(root))


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
        image = str(container.get("image") or "").strip() or "bactolomew-tools:0.1.0"
        runtime = str(config.get("runtime") or "docker").strip().lower()
        if runtime == "singularity":
            if not _CR.sif_present(str(container.get("sif_path") or "")):
                missing.append(root / TOOLS_IMAGE_SENTINEL)
        elif not _CR.docker_available():
            missing.append(root / DOCKER_CLIENT_SENTINEL)
        elif not _CR.docker_image_present(image):
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


def _format_missing_item(path: Path) -> str:
    name = path.name
    if name == DOCKER_CLIENT_SENTINEL:
        return "Docker client in PATH"
    if name == FASTGEAR_IMAGE_SENTINEL:
        settings = configured_fastgear_settings()
        return f"fastGEAR Docker image: {settings.get('image')}"
    return str(path)


def _category_names() -> list[str]:
    return [name for name, _ in cat_specs_no_other]


def _normalize_category(value: str) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    categories = _category_names()
    if text in categories:
        return text
    lower_map = {cat.lower(): cat for cat in categories}
    return lower_map.get(text.lower(), text)


def _infer_category_from_name(name: str) -> str:
    categories = _category_names()
    base = _strip_fasta_suffix(name)
    candidates: list[str] = [base]
    if "_" in base:
        candidates.insert(0, base.split("_")[-1])
    if "-" in base:
        candidates.insert(0, base.split("-")[-1])

    lower_map = {cat.lower(): cat for cat in categories}
    for candidate in candidates:
        key = candidate.strip()
        if not key:
            continue
        if key in categories:
            return key
        mapped = lower_map.get(key.lower())
        if mapped:
            return mapped

    low_base = base.lower()
    for category in sorted(categories, key=len, reverse=True):
        if category.lower() in low_base:
            return category
    return ""


def _load_catalog_template(path: Path) -> dict[str, object]:
    if not path.exists() or not path.is_file():
        return {}
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return raw if isinstance(raw, dict) else {}


def _find_template_entry(template: dict[str, object], fasta_path: Path) -> Optional[dict[str, object]]:
    stem = _strip_fasta_suffix(fasta_path.name)
    for key in (stem, fasta_path.stem, fasta_path.name):
        info = template.get(key)
        if isinstance(info, dict):
            return info

    for info in template.values():
        if not isinstance(info, dict):
            continue
        db_path_raw = info.get("db_path") or info.get("path") or info.get("db-path")
        if not db_path_raw:
            continue
        db_name = _strip_fasta_suffix(Path(str(db_path_raw)).name)
        if db_name == stem:
            return info
    return None


def _find_metadata_file(metadata_root: Path, source: str, basename: str) -> Optional[Path]:
    if not basename:
        return None

    search_roots: list[Path] = []
    source = str(source or "").strip()
    if source:
        search_roots.append(metadata_root / source)
    search_roots.append(metadata_root)

    seen: set[Path] = set()
    for root in search_roots:
        root = _normalized(root)
        if root in seen or not root.exists() or not root.is_dir():
            continue
        seen.add(root)

        direct = root / basename
        if direct.is_file():
            return direct
        raw = root / "raw" / basename
        if raw.is_file():
            return raw
        for candidate in root.rglob(basename):
            if candidate.is_file():
                return candidate
    return None


def _remap_resistance_metadata(value: object, metadata_root: Path) -> Optional[dict[str, object]]:
    normalized = normalize_resistance_metadata(value)
    if normalized is None:
        return None

    source = str(normalized.get("source") or "").strip()
    remapped: list[str] = []
    for raw_file in normalized.get("files") or []:
        basename = Path(str(raw_file)).name
        resolved = _find_metadata_file(metadata_root, source, basename)
        if resolved is not None:
            remapped.append(str(resolved))

    if not remapped:
        return None

    payload: dict[str, object] = {"files": remapped}
    if source:
        payload["source"] = source
    return payload


def _build_catalog_payload(
    sequence_dir: Path,
    template_path: Path,
    *,
    metadata_root: Optional[Path] = None,
) -> dict[str, dict[str, object]]:
    if not sequence_dir.exists() or not sequence_dir.is_dir():
        raise RuntimeError(f"Database directory not found: {sequence_dir}")

    template = _load_catalog_template(template_path)
    fasta_paths = sorted(
        (path for path in sequence_dir.iterdir() if path.is_file() and _is_fasta(path)),
        key=lambda path: path.name.lower(),
    )
    if not fasta_paths:
        raise RuntimeError(f"No FASTA databases were found in {sequence_dir}")

    payload: dict[str, dict[str, object]] = {}
    for fasta_path in fasta_paths:
        stem = _strip_fasta_suffix(fasta_path.name)
        template_entry = _find_template_entry(template, fasta_path)
        category = ""
        if isinstance(template_entry, dict):
            category = _normalize_category(str(template_entry.get("category") or ""))
        if not category:
            category = _infer_category_from_name(fasta_path.name)

        entry: dict[str, object] = {
            "category": category,
            "db_path": str(fasta_path),
        }

        if metadata_root is not None:
            template_meta = None
            if isinstance(template_entry, dict):
                template_meta = _remap_resistance_metadata(template_entry.get("resistance_metadata"), metadata_root)
            if template_meta is None:
                template_meta = _remap_resistance_metadata(
                    known_marker_resistance_metadata(stem, fasta_path),
                    metadata_root,
                )
            if template_meta is not None:
                entry["resistance_metadata"] = template_meta

        key = stem
        if key in payload:
            suffix = 2
            while f"{stem}_{suffix}" in payload:
                suffix += 1
            key = f"{stem}_{suffix}"
        payload[key] = entry

    return payload


def _write_catalog(path: Path, payload: dict[str, dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _is_db_latest_root(path: Path) -> bool:
    path = Path(path).expanduser()
    return (
        path.exists()
        and path.is_dir()
        and (path / "seq" / "markers").is_dir()
        and (path / "seq" / "compound").is_dir()
    )


def _candidate_dirs(path: Path) -> list[Path]:
    path = _normalized(path)
    candidates: list[Path] = []
    if path.is_file():
        candidates.append(path.parent)
    else:
        candidates.append(path)

    for base in list(candidates):
        candidates.append(base / DB_LATEST_DIRNAME)
        if base.name == DB_LATEST_DIRNAME:
            candidates.append(base.parent)

    dedup: list[Path] = []
    seen: set[Path] = set()
    for candidate in candidates:
        candidate = _normalized(candidate)
        if candidate in seen:
            continue
        seen.add(candidate)
        dedup.append(candidate)
    return dedup


def _find_db_latest_root(path: Path) -> Optional[Path]:
    for candidate in _candidate_dirs(path):
        if _is_db_latest_root(candidate):
            return candidate
        if not candidate.exists() or not candidate.is_dir():
            continue
        for child in sorted(candidate.iterdir()):
            if child.is_dir() and _is_db_latest_root(child):
                return child
    return None


@dataclass
class InstalledDbBundle:
    source_archive: Path


def _install_db_bundle(source_archive: Path, project_db_root: Path) -> InstalledDbBundle:
    source_archive = _normalized(source_archive)
    if not _is_supported_archive(source_archive):
        raise RuntimeError(f"Selected database archive is missing or unsupported: {source_archive}")

    project_db_root.parent.mkdir(parents=True, exist_ok=True)
    temp_root = Path(tempfile.mkdtemp(prefix="bactolomew_db_", dir=str(project_db_root.parent)))
    try:
        unpack_dir = temp_root / "unpacked"
        unpack_dir.mkdir(parents=True, exist_ok=True)
        shutil.unpack_archive(str(source_archive), str(unpack_dir))

        latest_source_dir = _find_db_latest_root(unpack_dir)
        if latest_source_dir is None:
            raise RuntimeError(
                "Selected database archive did not unpack to a usable BACTolomew bundle. "
                f"Expected {DB_LATEST_DIRNAME}/ with seq/, metadata/, cgmlst/, and catalogs/."
            )

        prepared_tree = temp_root / DB_DIRNAME
        shutil.copytree(latest_source_dir, prepared_tree)

        if project_db_root.exists():
            if not project_db_root.is_dir():
                raise RuntimeError(f"Project DB path exists but is not a directory: {project_db_root}")
            shutil.rmtree(project_db_root)
        project_db_root.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(prepared_tree), str(project_db_root))
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)

    return InstalledDbBundle(source_archive=source_archive)


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


class PathSelector(QtWidgets.QWidget):
    def __init__(
        self,
        *,
        caption: str,
        initial_dir: Optional[Path] = None,
        file_mode: bool = False,
        file_filter: str = "All files (*)",
    ) -> None:
        super().__init__()
        self._caption = caption
        self._initial_dir = initial_dir or Path.home()
        self._file_mode = bool(file_mode)
        self._file_filter = file_filter

        self.edit = QtWidgets.QLineEdit()
        self.button = QtWidgets.QPushButton("Browse...")
        self.button.clicked.connect(self._browse)

        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.edit, 1)
        layout.addWidget(self.button, 0)

    def _browse(self) -> None:
        current = self.path()
        start_dir = Path(current) if current else self._initial_dir
        if self._file_mode:
            path, _ = QtWidgets.QFileDialog.getOpenFileName(
                self,
                self._caption,
                str(start_dir),
                self._file_filter,
            )
            if path:
                self.set_path(Path(path))
            return
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, self._caption, str(start_dir))
        if folder:
            self.set_path(Path(folder))

    def set_path(self, path: Path) -> None:
        self.edit.setText(str(_normalized(path)))

    def path(self) -> str:
        return self.edit.text().strip()


@dataclass
class FirstTimeSetupResult:
    dev_mode: bool
    bakta_db_type: str
    bakta_db_path: Optional[Path]
    db_archive_path: Optional[Path]


class FirstTimeSetupDialog(QtWidgets.QDialog):
    def __init__(
        self,
        missing: Optional[list[Path]] = None,
        parent: Optional[QtWidgets.QWidget] = None,
        *,
        dev_mode: bool = False,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("First-Time Setup")
        self.resize(760, 440)
        self._result: Optional[FirstTimeSetupResult] = None
        self._dev_mode = bool(dev_mode)

        root = project_root()
        local_db_ready = _project_db_layout_ready(root)
        local_bakta_db = resolve_bakta_db_path(root / BAKTA_DB_DIRNAME)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(14, 14, 14, 14)
        layout.setSpacing(10)

        if local_db_ready:
            intro_text = f"Local DB found at {root / DB_DIRNAME}. Setup can reuse it."
        else:
            intro_text = "Choose the BACTolomew DB archive."
        intro = QtWidgets.QLabel(intro_text)
        intro.setWordWrap(True)
        layout.addWidget(intro)

        fastgear_note = QtWidgets.QLabel(
            "fastGEAR uses Docker. The image will be pulled if needed."
        )
        fastgear_note.setWordWrap(True)
        layout.addWidget(fastgear_note)

        if missing:
            missing_text = "\n".join(f"- {_format_missing_item(p)}" for p in missing)
            missing_label = QtWidgets.QLabel("Needs setup:\n" + missing_text)
            missing_label.setWordWrap(True)
            layout.addWidget(missing_label)

        bakta_group = QtWidgets.QGroupBox("Tool Environments And Tool DBs")
        bakta_layout = QtWidgets.QGridLayout(bakta_group)
        bakta_layout.addWidget(
            QtWidgets.QLabel(
                f"Envs: {root / ENVS_DIRNAME}\n"
                f"YAMLs: {root / YAML_DIRNAME}"
            ),
            0,
            0,
            1,
            2,
        )
        self.bakta_db_type_combo = QtWidgets.QComboBox()
        self.bakta_db_type_combo.addItems(["full", "light"])
        if self._dev_mode:
            self.bakta_db_type_combo.setEnabled(False)
            bakta_layout.addWidget(
                QtWidgets.QLabel(
                    f"Dev mode: Bakta is fixed to:\n{DEV_BAKTA_DB_PATH}"
                ),
                1,
                0,
                1,
                2,
            )
        else:
            bakta_layout.addWidget(QtWidgets.QLabel("Bakta DB folder"), 1, 0)
            self.bakta_db_selector = PathSelector(
                caption="Select Bakta DB folder",
                initial_dir=(local_bakta_db.parent if local_bakta_db is not None else root / BAKTA_DB_DIRNAME),
                file_mode=False,
            )
            if local_bakta_db is not None:
                self.bakta_db_selector.set_path(local_bakta_db)
            bakta_layout.addWidget(self.bakta_db_selector, 1, 1)
            bakta_layout.addWidget(QtWidgets.QLabel("If empty, download a local Bakta DB."), 2, 0, 1, 2)
            bakta_layout.addWidget(QtWidgets.QLabel("DB type if downloading"), 3, 0)
            bakta_layout.addWidget(self.bakta_db_type_combo, 3, 1)
            bakta_layout.addWidget(
                QtWidgets.QLabel(f"Local Bakta DB path: {root / BAKTA_DB_DIRNAME}"),
                4,
                0,
                1,
                2,
            )
        bakta_layout.addWidget(
            QtWidgets.QLabel("MOB-suite initializes its own database when the MOB-suite env is installed."),
            5,
            0,
            1,
            2,
        )
        layout.addWidget(bakta_group)

        db_group = QtWidgets.QGroupBox("Database Bundle")
        db_layout = QtWidgets.QVBoxLayout(db_group)
        db_layout.addWidget(
            QtWidgets.QLabel(
                "Choose a DB archive to install or replace the local DB. Leave empty to reuse `DB/`."
            )
        )
        default_archive = _default_db_archive_path()
        initial_dir = default_archive.parent if default_archive is not None else (
            DB_DEFAULT_SOURCE_ROOT if DB_DEFAULT_SOURCE_ROOT.exists() else Path.home()
        )
        self.db_source_selector = PathSelector(
            caption="Select BACTolomew database archive",
            initial_dir=initial_dir,
            file_mode=True,
            file_filter="Archives (*.tar.xz *.txz *.tar.gz *.tgz *.tar.bz2 *.tbz2 *.tar *.zip);;All files (*)",
        )
        if not _project_db_layout_ready(root) and default_archive is not None:
            self.db_source_selector.set_path(default_archive)
        db_layout.addWidget(self.db_source_selector)
        db_layout.addWidget(
            QtWidgets.QLabel(
                f"Local DB ready: {'yes' if local_db_ready else 'no'}\n"
                f"Path: {root / DB_DIRNAME}"
            )
        )
        layout.addWidget(db_group)

        btns = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel
        )
        btns.accepted.connect(self._on_accept)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)

    def result(self) -> Optional[FirstTimeSetupResult]:
        return self._result

    def _on_accept(self) -> None:
        result = self._validate()
        if result is None:
            return
        self._result = result
        self.accept()

    def _validate(self) -> Optional[FirstTimeSetupResult]:
        errors: List[str] = []
        root = project_root()

        bakta_db_path: Optional[Path] = None
        db_archive_path: Optional[Path] = None
        source_text = self.db_source_selector.path()
        if source_text:
            candidate = Path(source_text).expanduser()
            if not candidate.exists() or not candidate.is_file():
                errors.append(f"Database archive not found: {candidate}")
            else:
                db_archive_path = _normalized(candidate)
                if not _is_supported_archive(db_archive_path):
                    errors.append(
                        "Selected database file is not a supported archive "
                        f"({', '.join(ARCHIVE_SUFFIXES)})."
                    )
        elif not _project_db_layout_ready(root):
            errors.append("Select the BACTolomew database archive.")

        if not self._dev_mode:
            bakta_text = self.bakta_db_selector.path().strip()
            if bakta_text:
                candidate = Path(bakta_text).expanduser()
                resolved_bakta = resolve_bakta_db_path(candidate)
                if resolved_bakta is None:
                    errors.append(f"Bakta DB not found under: {candidate}")
                else:
                    bakta_db_path = _normalized(resolved_bakta)

        if errors:
            QtWidgets.QMessageBox.critical(self, "Setup error", "\n".join(errors))
            return None

        bakta_db_type = self.bakta_db_type_combo.currentText().strip().lower()
        if bakta_db_type not in {"full", "light"}:
            bakta_db_type = "full"

        return FirstTimeSetupResult(
            dev_mode=self._dev_mode,
            bakta_db_type=bakta_db_type,
            bakta_db_path=bakta_db_path,
            db_archive_path=db_archive_path,
        )


def _ensure_dir_preserve(path: Path, label: str) -> None:
    if path.exists():
        if not path.is_dir():
            raise RuntimeError(f"{label} path exists but is not a directory: {path}")
        return
    path.mkdir(parents=True, exist_ok=True)


def _run_cmd(cmd: list[str], *, cwd: Optional[Path] = None, env: Optional[dict[str, str]] = None) -> None:
    subprocess.run(cmd, check=True, cwd=str(cwd) if cwd else None, env=env)


def _create_tool_envs(root: Path, tools: Optional[Iterable[str]] = None) -> dict[str, Path]:
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
                _run_cmd(cmd)
                if not _conda_env_matches_yaml(env_path, yaml_path):
                    raise RuntimeError(f"Failed syncing conda env to YAML: {env_path}")
            env_paths[tool] = env_path.resolve()
            continue
        if env_path.exists() and not _is_conda_env(env_path):
            raise RuntimeError(f"Env path exists but is not a conda env: {env_path}")

        cmd = ["conda", "env", "create", "-p", str(env_path), "-f", str(yaml_path)]
        _run_cmd(cmd)
        if not _is_conda_env(env_path):
            raise RuntimeError(f"Failed creating conda env: {env_path}")
        env_paths[tool] = env_path.resolve()

    return env_paths


def _ensure_bakta_db(bakta_env: Optional[Path], db_root: Path, db_type: str) -> Path:
    db_root = db_root.resolve()
    db_root.mkdir(parents=True, exist_ok=True)

    existing = resolve_bakta_db_path(db_root)
    if existing is not None:
        return existing

    if bakta_env is not None:
        cmd = [
            "conda",
            "run",
            "-p",
            str(bakta_env),
            "bakta_db",
            "download",
            "--output",
            str(db_root),
            "--type",
            db_type,
        ]
    else:
        cmd = [
            "bakta_db",
            "download",
            "--output",
            str(db_root),
            "--type",
            db_type,
        ]
    _run_cmd(cmd)

    resolved = resolve_bakta_db_path(db_root)
    if resolved is None:
        raise RuntimeError(f"Bakta DB download finished but no DB was found under {db_root}")
    return resolved


def _ensure_checkm2_db(assembly_qc_env: Optional[Path]) -> Path:
    if assembly_qc_env is None:
        raise RuntimeError("assembly_qc conda environment was not created.")
    env_path = Path(assembly_qc_env).expanduser().resolve()
    if not _is_conda_env(env_path):
        raise RuntimeError(f"assembly_qc path is not a conda env: {env_path}")

    db_root = env_path / CHECKM2_DB_DIRNAME
    db_path = db_root / CHECKM2_DB_RELATIVE_PATH
    if db_path.is_file():
        return db_path.resolve()

    db_root.mkdir(parents=True, exist_ok=True)
    mpl_cache = env_path / "var" / "matplotlib"
    mpl_cache.mkdir(parents=True, exist_ok=True)
    run_env = dict(os.environ)
    run_env["MPLCONFIGDIR"] = str(mpl_cache)

    _run_cmd(
        [
            "conda",
            "run",
            "-p",
            str(env_path),
            "checkm2",
            "database",
            "--download",
            "--path",
            str(db_root),
        ],
        env=run_env,
    )
    if db_path.is_file():
        return db_path.resolve()

    candidates = sorted(db_root.rglob("*.dmnd"))
    if candidates:
        return candidates[0].resolve()
    raise RuntimeError(f"CheckM2 database download completed but no DIAMOND DB was found under: {db_root}")


def apply_first_time_setup(result: FirstTimeSetupResult) -> None:
    root = project_root()
    runs_dir = root / "runs"
    cache_dir = root / "cache"
    db_root = _project_db_root(root)
    catalogs_dir = _project_catalog_dir(root)
    marker_catalog_path = _project_marker_catalog_path(root)
    compound_catalog_path = _project_compound_catalog_path(root)
    marker_dir = _project_marker_dir(root)
    compound_dir = _project_compound_dir(root)
    metadata_dir = _project_metadata_dir(root)
    cgmlst_dir = _project_cgmlst_dir(root)
    config_path = root / "bactolomew_config.json"

    existing_config: dict = {}
    if config_path.exists():
        try:
            existing_raw = json.loads(config_path.read_text(encoding="utf-8"))
        except Exception:
            existing_raw = {}
        if isinstance(existing_raw, dict):
            existing_config = existing_raw

    _ensure_dir_preserve(runs_dir, "runs")
    _ensure_dir_preserve(cache_dir, "cache")

    existing_db_cfg = existing_config.get("databases") if isinstance(existing_config.get("databases"), dict) else {}
    installed_bundle: Optional[InstalledDbBundle] = None
    if result.db_archive_path is not None:
        installed_bundle = _install_db_bundle(result.db_archive_path, db_root)
    elif not _project_db_layout_ready(root):
        raise RuntimeError("No project-local DB is installed and no database archive was provided.")

    if not catalogs_dir.exists() or not catalogs_dir.is_dir():
        raise RuntimeError(f"Installed DB is missing catalog directory: {catalogs_dir}")
    if not _catalog_targets_dir(marker_catalog_path, marker_dir):
        raise RuntimeError(
            "Installed marker DB catalog is missing or invalid. "
            f"Expected a bundle-owned catalog at {marker_catalog_path} targeting {marker_dir}."
        )
    if not _catalog_targets_dir(compound_catalog_path, compound_dir):
        raise RuntimeError(
            "Installed compound DB catalog is missing or invalid. "
            f"Expected a bundle-owned catalog at {compound_catalog_path} targeting {compound_dir}."
        )

    available_schemas = _discover_chewbbaca_schemas(cgmlst_dir)
    existing_chewbbaca = existing_config.get("chewbbaca") if isinstance(existing_config.get("chewbbaca"), dict) else {}
    selected_schema = _preferred_schema_path(available_schemas, existing_chewbbaca)

    if result.dev_mode:
        env_paths = _create_tool_envs(root)
        envs = {name: str(path) for name, path in env_paths.items()}
        bakta_db = Path(DEV_BAKTA_DB_PATH)
    else:
        env_paths = _create_tool_envs(root)
        if result.bakta_db_path is not None:
            resolved_bakta = resolve_bakta_db_path(result.bakta_db_path)
            if resolved_bakta is None:
                raise RuntimeError(f"Bakta DB not found under: {result.bakta_db_path}")
            bakta_db = resolved_bakta
        else:
            bakta_db = _ensure_bakta_db(env_paths.get("bakta"), root / BAKTA_DB_DIRNAME, result.bakta_db_type)
        envs = {name: str(path) for name, path in env_paths.items()}

    checkm2_db = _ensure_checkm2_db(env_paths.get("assembly_qc"))

    db_config: dict[str, object] = {
        "root_dir": str(db_root),
        "marker_sequences_dir": str(marker_dir),
        "compound_sequences_dir": str(compound_dir),
        "marker_catalog_path": str(marker_catalog_path),
        "compound_catalog_path": str(compound_catalog_path),
        "metadata_dir": str(metadata_dir),
        "cgmlst_dir": str(cgmlst_dir),
    }
    cgmlst_catalog_path = _project_cgmlst_catalog_path(root)
    if cgmlst_catalog_path.is_file():
        db_config["cgmlst_scheme_catalog_path"] = str(cgmlst_catalog_path)
    manifest_path = metadata_dir / "manifest.tsv"
    if manifest_path.is_file():
        db_config["manifest_path"] = str(manifest_path)

    if installed_bundle is not None:
        db_config["source_archive_path"] = str(installed_bundle.source_archive)
    else:
        for key in ("source_archive_path", "manifest_path"):
            value = str(existing_db_cfg.get(key) or "").strip()
            if value and key not in db_config:
                db_config[key] = value

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
        "dev_mode": bool(result.dev_mode),
        "databases": db_config,
        "created_at": datetime.now().astimezone().isoformat(),
    }

    existing_fastgear = existing_config.get("fastgear")
    fastgear_settings = configured_fastgear_settings(
        existing_fastgear if isinstance(existing_fastgear, dict) else None
    )
    existing_assembly_qc = existing_config.get("assembly_qc")
    config["fastgear"] = _ensure_fastgear_image(fastgear_settings)
    config["chewbbaca"] = chewbbaca_config
    config["mobsuite"] = {}
    assembly_qc_config = configured_assembly_qc_settings(
        existing_assembly_qc if isinstance(existing_assembly_qc, dict) else None
    )
    assembly_qc_config["checkm2_database_path"] = str(checkm2_db)
    config["assembly_qc"] = assembly_qc_config

    config_path.write_text(json.dumps(config, indent=2, sort_keys=True), encoding="utf-8")
