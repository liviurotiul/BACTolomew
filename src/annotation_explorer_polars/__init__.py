from __future__ import annotations

from importlib import import_module


_LAZY_EXPORTS = {
    "ExploreWindow": ("annotation_explorer_polars.windows.explore_window", "ExploreWindow"),
    "ANNOTATION_COLUMNS_TO_SHOW": ("annotation_explorer_polars.constants", "ANNOTATION_COLUMNS_TO_SHOW"),
    "ANNOTATIONS_FILENAME": ("annotation_explorer_polars.constants", "ANNOTATIONS_FILENAME"),
    "EXPECTED_ANNOTATION_COLUMNS": ("annotation_explorer_polars.constants", "EXPECTED_ANNOTATION_COLUMNS"),
    "EXPECTED_REFERENCE_LOCUS_COLUMNS": (
        "annotation_explorer_polars.constants",
        "EXPECTED_REFERENCE_LOCUS_COLUMNS",
    ),
    "REFERENCE_LOCUS_FILENAME": ("annotation_explorer_polars.constants", "REFERENCE_LOCUS_FILENAME"),
    "StartupDialog": ("annotation_explorer_polars.startup_dialog", "StartupDialog"),
    "apply_annotation_filters": ("annotation_explorer_polars.data", "apply_annotation_filters"),
    "apply_metadata_filters": ("annotation_explorer_polars.data", "apply_metadata_filters"),
    "load_annotations": ("annotation_explorer_polars.data", "load_annotations"),
    "load_cache": ("annotation_explorer_polars.cache", "load_cache"),
    "load_metadata": ("annotation_explorer_polars.data", "load_metadata"),
    "load_reference_locus_samples": ("annotation_explorer_polars.data", "load_reference_locus_samples"),
    "main": ("annotation_explorer_polars.app", "main"),
    "prepare_metadata": ("annotation_explorer_polars.data", "prepare_metadata"),
    "save_cache": ("annotation_explorer_polars.cache", "save_cache"),
}

__all__ = sorted(_LAZY_EXPORTS)


def __getattr__(name: str) -> object:
    target = _LAZY_EXPORTS.get(name)
    if target is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr_name = target
    module = import_module(module_name)
    value = getattr(module, attr_name)
    globals()[name] = value
    return value
