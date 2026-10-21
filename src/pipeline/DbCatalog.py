from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional, Sequence

from pipeline.Utils import CATEGORIES


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DB_ROOT = PROJECT_ROOT / "DB"
CATALOG_DIR = DB_ROOT / "catalogs"
MARKER_DB_CATALOG_PATH = CATALOG_DIR / "marker_db_catalog.json"
COMPOUND_DB_CATALOG_PATH = CATALOG_DIR / "compound_db_catalog.json"
RESISTANCE_METADATA_ROOT = DB_ROOT / "metadata"
AMRF_METADATA_DIR = RESISTANCE_METADATA_ROOT / "AMRF"
CARD_METADATA_DIR = RESISTANCE_METADATA_ROOT / "CARD"

KNOWN_MARKER_RESISTANCE_METADATA: dict[str, dict[str, Any]] = {
    "AMRF_core-AMR": {
        "source": "AMRF",
        "files": [
            str(AMRF_METADATA_DIR / "ReferenceGeneCatalog.txt"),
        ],
    },
    "AMRF_plus-AMR": {
        "source": "AMRF",
        "files": [
            str(AMRF_METADATA_DIR / "ReferenceGeneCatalog.txt"),
        ],
    },
    "AMRF_core-BMRG": {
        "source": "AMRF",
        "files": [
            str(AMRF_METADATA_DIR / "ReferenceGeneCatalog.txt"),
        ],
    },
    "AMRF_plus-BMRG": {
        "source": "AMRF",
        "files": [
            str(AMRF_METADATA_DIR / "ReferenceGeneCatalog.txt"),
        ],
    },
    "CARD_AMR": {
        "source": "CARD",
        "files": [
            str(CARD_METADATA_DIR / "aro_index.tsv"),
        ],
    },
}


def _read_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Catalog must be a JSON object: {path}")
    return data


def _resolve_catalog_path(catalog_path: Path, value: Any) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = catalog_path.parent / path
    try:
        return str(path.resolve())
    except Exception:
        return str(path)


def normalize_resistance_metadata(value: Any) -> dict[str, Any] | None:
    if not isinstance(value, dict):
        return None

    source = str(value.get("source") or value.get("resolver") or "").strip()
    files_raw = value.get("files") or value.get("paths") or value.get("metadata_files") or []
    if isinstance(files_raw, (str, Path)):
        files = [str(Path(files_raw).expanduser())]
    elif isinstance(files_raw, (list, tuple, set)):
        files = [str(Path(v).expanduser()) for v in files_raw if str(v).strip()]
    else:
        files = []

    if not source and not files:
        return None

    out: dict[str, Any] = {
        "files": files,
    }
    if source:
        out["source"] = source
    return out


def known_marker_resistance_metadata(name: str, db_path: str | Path | None = None) -> dict[str, Any] | None:
    candidates: list[str] = [str(name).strip()]
    if db_path is not None:
        candidates.append(Path(db_path).expanduser().stem)
    for candidate in candidates:
        if not candidate:
            continue
        meta = KNOWN_MARKER_RESISTANCE_METADATA.get(candidate)
        if meta:
            return normalize_resistance_metadata(meta)
    return None


def load_db_catalog(path: Path, *, logger: Optional[logging.Logger] = None) -> dict[str, dict[str, Any]]:
    if not path.exists():
        raise FileNotFoundError(f"Catalog not found: {path}")
    path = path.expanduser()
    raw = _read_json(path)
    out: dict[str, dict[str, Any]] = {}
    for name, info in raw.items():
        if not isinstance(info, dict):
            if logger:
                logger.warning("Catalog entry '%s' is not an object; skipping.", name)
            continue
        db_path = info.get("db_path") or info.get("path") or info.get("db-path")
        if not db_path:
            if logger:
                logger.warning("Catalog entry '%s' missing db_path; skipping.", name)
            continue
        category = info.get("category")
        normalized = dict(info)
        normalized["db_path"] = _resolve_catalog_path(path, db_path)
        normalized["category"] = "" if category is None else str(category)

        resistance_metadata = normalize_resistance_metadata(info.get("resistance_metadata"))
        if resistance_metadata is not None:
            resistance_metadata = {
                **resistance_metadata,
                "files": [_resolve_catalog_path(path, item) for item in resistance_metadata.get("files") or []],
            }
        else:
            resistance_metadata = known_marker_resistance_metadata(str(name), normalized["db_path"])
        if resistance_metadata is not None:
            normalized["resistance_metadata"] = resistance_metadata

        out[str(name)] = normalized
    return out


@dataclass
class DbCatalogIndex:
    catalog: dict[str, dict[str, Any]]
    by_path: dict[Path, str]
    by_stem: dict[str, str]

    @classmethod
    def from_catalog(cls, catalog: dict[str, dict[str, Any]]) -> "DbCatalogIndex":
        by_path: dict[Path, str] = {}
        by_stem: dict[str, str] = {}
        for name, info in (catalog or {}).items():
            p_raw = info.get("db_path") or info.get("path")
            if not p_raw:
                continue
            p = Path(p_raw).expanduser()
            try:
                p_res = p.resolve()
            except Exception:
                p_res = p
            by_path[p_res] = name
            by_stem.setdefault(p.stem, name)
        return cls(catalog=catalog or {}, by_path=by_path, by_stem=by_stem)

    def _match(self, path: Path) -> dict[str, Any] | None:
        try:
            if path.exists():
                p_res = path.resolve()
                name = self.by_path.get(p_res)
                if name:
                    return self.catalog.get(name)
        except Exception:
            pass

        if path.stem in self.catalog:
            return self.catalog.get(path.stem)
        if path.name in self.catalog:
            return self.catalog.get(path.name)

        name = self.by_stem.get(path.stem)
        if name:
            return self.catalog.get(name)
        return None

    def entry_for_path(self, path: Path) -> dict[str, Any] | None:
        return self._match(path)

    def category_for_path(
        self,
        path: Path,
        *,
        logger: Optional[logging.Logger] = None,
        label: str = "database",
    ) -> str:
        entry = self._match(path)
        if entry is None:
            if logger:
                logger.warning("No catalog entry for %s: %s", label, path)
            return ""

        cat = str(entry.get("category") or "").strip()
        if cat and cat not in CATEGORIES:
            if logger:
                logger.warning(
                    "Catalog category '%s' for %s is not in CATEGORIES; keeping as-is.",
                    cat, path.name
                )
        return cat

    def category_map(
        self,
        paths: list[Path],
        *,
        logger: Optional[logging.Logger] = None,
        label: str = "database",
    ) -> dict[str, str]:
        out: dict[str, str] = {}
        for p in paths:
            out[p.stem] = self.category_for_path(p, logger=logger, label=label)
        return out

    def resistance_metadata_for_path(
        self,
        path: Path,
        *,
        logger: Optional[logging.Logger] = None,
        label: str = "database",
    ) -> dict[str, Any] | None:
        entry = self._match(path)
        if entry is None:
            if logger:
                logger.warning("No catalog entry for %s: %s", label, path)
            return None
        meta = normalize_resistance_metadata(entry.get("resistance_metadata"))
        if meta is None:
            meta = known_marker_resistance_metadata(path.stem, entry.get("db_path"))
        return meta

    def resistance_metadata_map(
        self,
        paths: list[Path],
        *,
        logger: Optional[logging.Logger] = None,
        label: str = "database",
    ) -> dict[str, dict[str, Any]]:
        out: dict[str, dict[str, Any]] = {}
        for p in paths:
            meta = self.resistance_metadata_for_path(p, logger=logger, label=label)
            if meta is not None:
                out[p.stem] = meta
        return out


def load_db_catalog_index(path: Path, *, logger: Optional[logging.Logger] = None) -> DbCatalogIndex:
    catalog = load_db_catalog(path, logger=logger)
    return DbCatalogIndex.from_catalog(catalog)


@dataclass(frozen=True)
class CatalogDb:
    """One database a catalog lists; it is installed when its FASTA file exists. The
    marker and compound databases are optional: a run searches the installed ones."""
    name: str
    category: str
    path: Path

    @property
    def installed(self) -> bool:
        return self.path.is_file()


def catalog_databases(path: Path, *, logger: Optional[logging.Logger] = None) -> list[CatalogDb]:
    """Every database of a catalog, installed or not (raises when the catalog cannot be read)."""
    return [
        CatalogDb(name, str(info.get("category") or ""), Path(info["db_path"]))
        for name, info in load_db_catalog(path, logger=logger).items()
    ]


# sha256 of the database files, per (path, size, mtime_ns): a run records it more than once.
_SHA256_CACHE: dict[tuple[str, int, int], str] = {}


def _database_file_record(path: Path) -> dict[str, Any]:
    out: dict[str, Any] = {"path": str(path)}
    try:
        st = path.stat()
    except OSError:
        out["missing"] = True
        return out
    out["size"] = st.st_size
    out["mtime"] = datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat(timespec="seconds")
    key = (str(path), st.st_size, st.st_mtime_ns)
    if key not in _SHA256_CACHE:
        digest = hashlib.sha256()
        with path.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                digest.update(chunk)
        _SHA256_CACHE[key] = digest.hexdigest()
    out["sha256"] = _SHA256_CACHE[key]
    return out


def run_database_record(used: Sequence[Path], catalog_path: Path) -> dict[str, list[dict[str, Any]]]:
    """The databases of one kind a run searches ("used": name, category, path, size, mtime,
    sha256) and the catalog's databases it could not ("not_installed": their file is
    missing). Never raises: an unreadable catalog leaves the used ones unnamed."""
    try:
        listed = catalog_databases(catalog_path)
    except Exception:
        listed = []
    by_stem = {db.path.stem: db for db in listed}
    used_stems = {Path(p).stem for p in used}
    used_rows: list[dict[str, Any]] = []
    for raw in used:
        path = Path(raw)
        db = by_stem.get(path.stem)
        used_rows.append({
            "name": db.name if db else path.stem,
            "category": db.category if db else "",
            **_database_file_record(path),
        })
    missing = [
        {"name": db.name, "category": db.category, "path": str(db.path)}
        for db in listed
        if not db.installed and db.path.stem not in used_stems
    ]
    return {"used": used_rows, "not_installed": missing}


__all__ = [
    "CATALOG_DIR",
    "MARKER_DB_CATALOG_PATH",
    "COMPOUND_DB_CATALOG_PATH",
    "RESISTANCE_METADATA_ROOT",
    "AMRF_METADATA_DIR",
    "CARD_METADATA_DIR",
    "KNOWN_MARKER_RESISTANCE_METADATA",
    "normalize_resistance_metadata",
    "known_marker_resistance_metadata",
    "load_db_catalog",
    "DbCatalogIndex",
    "load_db_catalog_index",
    "CatalogDb",
    "catalog_databases",
    "run_database_record",
]
