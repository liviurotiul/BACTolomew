from __future__ import annotations

import hashlib
import json
import re
import time
from pathlib import Path
from typing import Any, Dict, Optional

from .polars_io import read_parquet_df
from pipeline.DbCatalog import MARKER_DB_CATALOG_PATH, load_db_catalog_index

CACHE_FILENAME = "ui_cache.json"
CACHE_DIRNAME = "cache"
RUNS_DIRNAME = "runs"
ANNOTATIONS_CACHE_FILENAME = "annotations.parquet"
REFERENCE_LOCUS_CACHE_FILENAME = "reference_locus_hits.parquet"
MARKER_HITS_CACHE_FILENAME = "marker_hits.parquet"
TABLE_CACHE_META_FILENAME = "tables_meta.json"
# Only the current cache schema is supported. Older metadata is invalidated
# and rebuilt in place instead of maintaining multiple legacy cache formats.
# v13: marker hits now collapse per (AnnotationID, DB) instead of (AnnotationID, Method),
#      so AMRF category databases are no longer merged. Bumping invalidates stale caches
#      built before the fix, forcing a rebuild from the existing DIAMOND TSVs.
TABLE_CACHE_SCHEMA_VERSION = 13


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def cache_root() -> Path:
    # Project root (same level as src/, runs/, catalogs/)
    return project_root() / CACHE_DIRNAME


def runs_root() -> Path:
    return project_root() / RUNS_DIRNAME


def ensure_cache_dir() -> Path:
    root = cache_root()
    root.mkdir(parents=True, exist_ok=True)
    return root


def ensure_runs_dir() -> Path:
    root = runs_root()
    root.mkdir(parents=True, exist_ok=True)
    return root


def _cache_path() -> Path:
    return ensure_cache_dir() / CACHE_FILENAME


def load_cache() -> Dict[str, str]:
    path = _cache_path()
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text())
    except Exception:
        return {}
    if not isinstance(data, dict):
        return {}
    out: Dict[str, str] = {}
    for key in ("output_dir", "metadata_csv"):
        val = data.get(key)
        if isinstance(val, str) and val.strip():
            out[key] = val
    return out


def save_cache(output_dir: Optional[Path], metadata_csv: Optional[Path]) -> None:
    data: Dict[str, Any] = {
        "output_dir": str(output_dir) if output_dir else "",
        "metadata_csv": str(metadata_csv) if metadata_csv else "",
    }
    path = _cache_path()
    path.write_text(json.dumps(data, indent=2, sort_keys=True))


def cache_run_dir(run_id: str) -> Path:
    root = ensure_cache_dir()
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", run_id or "run").strip("_") or "run"
    path = root / slug
    path.mkdir(parents=True, exist_ok=True)
    return path


def cached_table_paths(run_id: str) -> tuple[Path, Path, Path]:
    run_dir = cache_run_dir(run_id)
    return (
        run_dir / ANNOTATIONS_CACHE_FILENAME,
        run_dir / REFERENCE_LOCUS_CACHE_FILENAME,
        run_dir / MARKER_HITS_CACHE_FILENAME,
    )


def _cached_tables_meta_path(run_id: str) -> Path:
    return cache_run_dir(run_id) / TABLE_CACHE_META_FILENAME


def _read_cached_tables_meta(run_id: str) -> Dict[str, Any] | None:
    meta_path = _cached_tables_meta_path(run_id)
    if not meta_path.exists():
        return None
    try:
        payload = json.loads(meta_path.read_text())
    except Exception:
        return None
    return payload if isinstance(payload, dict) else None


def _normalize_output_dir(output_dir: Path) -> str:
    return str(Path(output_dir).expanduser().resolve())


def _run_config_metadata_paths(run_cfg_path: Path) -> list[Path]:
    try:
        payload = json.loads(run_cfg_path.read_text())
    except Exception:
        return []
    if not isinstance(payload, dict):
        return []
    inputs = payload.get("inputs")
    if not isinstance(inputs, dict):
        return []
    marker_meta = inputs.get("marker_resistance_metadata")
    resolved_meta = marker_meta if isinstance(marker_meta, dict) else None
    if resolved_meta is None:
        try:
            index = load_db_catalog_index(MARKER_DB_CATALOG_PATH)
        except Exception:
            resolved_meta = {}
        else:
            marker_paths = [
                Path(str(p)).expanduser()
                for p in (inputs.get("marker_sequences") or [])
                if str(p).strip()
            ]
            resolved_meta = index.resistance_metadata_map(marker_paths)

    out: list[Path] = []
    for raw in resolved_meta.values():
        meta = raw.get("resistance_metadata") if isinstance(raw, dict) and "resistance_metadata" in raw else raw
        if not isinstance(meta, dict):
            continue
        files = meta.get("files") or []
        if isinstance(files, (str, Path)):
            files = [files]
        if not isinstance(files, (list, tuple, set)):
            continue
        for file_path in files:
            text = str(file_path).strip()
            if text:
                out.append(Path(text).expanduser())
    return out


def _table_source_paths(output_dir: Path) -> list[Path]:
    root = Path(output_dir).expanduser().resolve()
    sources: list[Path] = []

    run_cfg_path = root / "run_config.json"
    if run_cfg_path.exists() and run_cfg_path.is_file():
        sources.append(run_cfg_path)
        for metadata_path in _run_config_metadata_paths(run_cfg_path):
            if metadata_path.exists() and metadata_path.is_file():
                sources.append(metadata_path.resolve())

    scan_specs = [
        (root / "bakta", "*.json"),
        (root / "chewbbaca" / "out", "*.tsv"),
        (root / "chewbbaca" / "out", "logging_info.txt"),
        (root / "mlst", "results.csv"),
        (root / "mlst", "results.json"),
        (root / "mobsuite", "contigs.tsv"),
        (root / "mobsuite", "plasmids.tsv"),
        (root / "mobsuite", "summary.json"),
        (root / "annotations" / "bakta_vs_queries", "*.tsv"),
        (root / "annotations" / "reference_loci", "*.tsv"),
    ]
    for base_dir, pattern in scan_specs:
        if not base_dir.exists():
            continue
        sources.extend(sorted(p for p in base_dir.rglob(pattern) if p.is_file()))
    return sources


def build_tables_source_fingerprint(output_dir: Path) -> str:
    root = Path(output_dir).expanduser().resolve()
    payload: Dict[str, Any] = {
        "output_dir": str(root),
        "files": [],
    }
    for path in _table_source_paths(root):
        try:
            stat = path.stat()
        except OSError:
            continue
        try:
            rel_path = str(path.relative_to(root))
        except ValueError:
            rel_path = str(path)
        payload["files"].append(
            {
                "path": rel_path,
                "size": int(stat.st_size),
                "mtime_ns": int(stat.st_mtime_ns),
            }
        )
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def cached_tables_exist(run_id: str, output_dir: Path) -> bool:
    ann_path, ref_path, marker_hits_path = cached_table_paths(run_id)
    if not (ann_path.exists() and ref_path.exists() and marker_hits_path.exists()):
        return False
    meta = _read_cached_tables_meta(run_id)
    if not isinstance(meta, dict):
        return False
    if meta.get("schema_version") != TABLE_CACHE_SCHEMA_VERSION:
        return False
    if meta.get("output_dir") != _normalize_output_dir(output_dir):
        return False
    expected_fingerprint = build_tables_source_fingerprint(output_dir)
    return meta.get("source_fingerprint") == expected_fingerprint


def write_cached_tables(run_id: str, output_dir: Path, ann_df, ref_df, marker_hits_df) -> tuple[Path, Path, Path]:
    ann_path, ref_path, marker_hits_path = cached_table_paths(run_id)
    try:
        ann_df.to_parquet(ann_path, index=False)
        ref_df.to_parquet(ref_path, index=False)
        marker_hits_df.to_parquet(marker_hits_path, index=False)
        meta_path = _cached_tables_meta_path(run_id)
        meta_path.write_text(
            json.dumps(
                {
                    "schema_version": TABLE_CACHE_SCHEMA_VERSION,
                    "output_dir": _normalize_output_dir(output_dir),
                    "source_fingerprint": build_tables_source_fingerprint(output_dir),
                },
                indent=2,
                sort_keys=True,
            )
        )
    except Exception as exc:
        raise RuntimeError(
            "Failed to write parquet cache tables. Install pyarrow to enable parquet support."
        ) from exc
    return ann_path, ref_path, marker_hits_path


def read_cached_tables(run_id: str):
    ann_path, ref_path, marker_hits_path = cached_table_paths(run_id)
    if not ann_path.exists() or not ref_path.exists() or not marker_hits_path.exists():
        raise FileNotFoundError("Cached parquet tables not found.")
    try:
        t0 = time.perf_counter()
        print(f"[file] reading {ann_path}")
        ann_df = read_parquet_df(ann_path)
        print(f"[file] read {ann_path} in {time.perf_counter() - t0:.2f}s")
        t1 = time.perf_counter()
        print(f"[file] reading {ref_path}")
        ref_df = read_parquet_df(ref_path)
        print(f"[file] read {ref_path} in {time.perf_counter() - t1:.2f}s")
        t2 = time.perf_counter()
        print(f"[file] reading {marker_hits_path}")
        marker_hits_df = read_parquet_df(marker_hits_path)
        print(f"[file] read {marker_hits_path} in {time.perf_counter() - t2:.2f}s")
    except Exception as exc:
        raise RuntimeError(
            "Failed to read parquet cache tables. Install pyarrow to enable parquet support."
        ) from exc
    return ann_df, ref_df, marker_hits_df


def run_cache_path(run_id: str) -> Path:
    return cache_run_dir(run_id) / CACHE_FILENAME


def load_run_cache(run_id: str) -> Dict[str, str]:
    path = run_cache_path(run_id)
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text())
    except Exception:
        return {}
    if not isinstance(data, dict):
        return {}
    out: Dict[str, str] = {}
    meta = data.get("metadata_csv")
    if isinstance(meta, str) and meta.strip():
        out["metadata_csv"] = meta
    return out


def save_run_cache(run_id: str, metadata_csv: Optional[Path]) -> None:
    data: Dict[str, Any] = {
        "metadata_csv": str(metadata_csv) if metadata_csv else "",
    }
    path = run_cache_path(run_id)
    path.write_text(json.dumps(data, indent=2, sort_keys=True))
