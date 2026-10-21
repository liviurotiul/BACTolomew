"""Dependency smoke test for tools/check_latest_deps.sh (no Qt window).

Imports every host dependency, then loads a small synthetic run the way the
app does: the annotation tables go through the parquet cache (written to a
temporary folder, never the project's cache/), are read back, built into
RunData and filtered by QueryEngine. Prints one line per check and exits 1 on
the first failure.

Usage (from the repository root):
  QT_QPA_PLATFORM=offscreen PYTHONPATH=src:tests VENV/bin/python tools/deps_smoke.py
"""
from __future__ import annotations

import importlib
import os
import sys
import tempfile
import traceback
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

MODULES = (
    "numpy", "pandas", "polars", "pyarrow", "scipy", "Bio",
    "matplotlib", "pyqtgraph", "PyQt6.QtWidgets", "PyQt6.QtWebEngineWidgets", "PyQt6Ads",
)


def _check(label, fn):
    try:
        detail = fn()
    except Exception:
        print(f"FAIL {label}")
        traceback.print_exc()
        sys.exit(1)
    print(f"ok   {label}" + (f": {detail}" if detail else ""))


def _imports():
    versions = []
    for name in MODULES:
        mod = importlib.import_module(name)
        version = getattr(mod, "__version__", "")
        if version:
            versions.append(f"{name}={version}")
    return " ".join(versions)


def _load_synthetic_run():
    from _qt_test_env import isolate_qt_settings

    isolate_qt_settings()
    import pandas as pd
    from test_refresh_parity import build_synthetic_run

    from annotation_explorer_polars import cache
    from annotation_explorer_polars.model import QueryEngine, RunData
    from annotation_explorer_polars.model.filter_state import FilterState

    ann_df, meta_df = build_synthetic_run()
    with tempfile.TemporaryDirectory(prefix="bactolomew-deps-") as tmp:
        cache.cache_root = lambda: Path(tmp) / "cache"  # never the project's cache/
        out_dir = Path(tmp) / "run"
        out_dir.mkdir()
        cache.write_cached_tables("smoke", out_dir, ann_df, pd.DataFrame(), pd.DataFrame())
        ann_back, ref_back, hits_back = cache.read_cached_tables("smoke")
        if len(ann_back) != len(ann_df):
            raise AssertionError(f"parquet round trip: {len(ann_back)} rows, expected {len(ann_df)}")
        run = RunData.build(
            ann_back, meta_df, ref_df=ref_back, marker_hits_df=hits_back, output_dir=out_dir, run_id="smoke"
        )
        engine = QueryEngine(run)
        everything = engine.compute(FilterState())
        searched = engine.compute(FilterState(gene="blaOXA"))
        if not searched.samples or len(searched.samples) >= len(everything.samples):
            raise AssertionError(
                f"gene search kept {len(searched.samples)} of {len(everything.samples)} samples"
            )
    return f"{len(ann_back)} annotations, {len(everything.samples)} samples, blaOXA in {len(searched.samples)}"


if __name__ == "__main__":
    _check("imports", _imports)
    _check("synthetic run through the parquet cache, RunData and QueryEngine", _load_synthetic_run)
