"""Run data and queries (plan 4.3). Shared read-only by every explore window of a run.

No module of this package loads Qt. The names are imported lazily, so
``model.run_data`` can be imported without polars-heavy query code.
"""
from __future__ import annotations

from importlib import import_module
from typing import Any

_LAZY_EXPORTS = {
    "AppliedResult": ".filter_state",
    "FilterState": ".filter_state",
    "PanData": ".pan_data",
    "QueryEngine": ".query_engine",
    "RunData": ".run_data",
    "SampleSubset": ".filter_state",
    "Selection": ".filter_state",
}

__all__ = sorted(_LAZY_EXPORTS)


def __getattr__(name: str) -> Any:
    module_name = _LAZY_EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(module_name, __name__), name)
    globals()[name] = value
    return value
