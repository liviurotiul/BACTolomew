"""The filter query and its applied result (window plan 4.4, Phase 3).

``FilterState`` is what one explore window asks the ``QueryEngine`` for and
``AppliedResult`` is what it gets back; ``Selection`` is the window's shared
selection. They replace ``RefreshSnapshot`` / ``RefreshPayload`` (formerly in
refresh_worker.py) and live here, next to the engine, so that the model package
loads no Qt module.

Decision D3: the annotation table's view mode (Contig, Around matches, All
matches) is not part of the query. The sample search always filters samples by
name and metadata; genes and products are searched only through ``gene`` and
``product``.
"""
from __future__ import annotations

import dataclasses
import json
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, Dict, FrozenSet, Iterable, List, Mapping, Optional, Sequence, Tuple

import polars as pl

MOLECULE_VALUES = ("plasmid", "chromosome")


def _text(value: Any) -> str:
    return "" if value is None else str(value).strip()


def _unique_texts(values: Iterable[Any]) -> Tuple[str, ...]:
    out: List[str] = []
    seen: set[str] = set()
    for raw in values or ():
        text = _text(raw)
        if text and text not in seen:
            seen.add(text)
            out.append(text)
    return tuple(out)


@dataclass(frozen=True)
class SampleSubset:
    """One sample-subset filter (the chip in the rail): its samples, label and kind."""

    samples: Tuple[str, ...]
    label: str = ""
    kind: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "samples", _unique_texts(self.samples))
        object.__setattr__(self, "label", _text(self.label))
        object.__setattr__(self, "kind", _text(self.kind))


def _normalize_metadata(value: Any) -> Tuple[Tuple[str, Tuple[str, ...]], ...]:
    items = value.items() if isinstance(value, Mapping) else (value or ())
    out: Dict[str, Tuple[str, ...]] = {}
    for column, values in items:
        column = _text(column)
        # The engine ignores empty value lists; so does the query key.
        normalized = tuple(sorted({str(v) for v in (values or ()) if v is not None}))
        if column and normalized:
            out[column] = normalized
    return tuple(sorted(out.items()))


@dataclass(frozen=True)
class FilterState:
    """The filters of one explore window. Hashable and order insensitive where
    the query is (predictors, tags, metadata values)."""

    sample_query: str = ""  # always sample name and metadata (D3)
    gene: str = ""  # "Gene/annotation contains" (also marker, predictor and method text)
    product: str = ""  # "Product contains"
    molecule: str = ""  # "", "plasmid" or "chromosome"
    predictors: FrozenSet[str] = frozenset()
    tags: FrozenSet[str] = frozenset()
    hide_hypothetical: bool = False
    hide_empty_contigs: bool = False
    metadata: Tuple[Tuple[str, Tuple[str, ...]], ...] = ()
    subset: Optional[SampleSubset] = None
    # A locked subset of a window opened from a clade and the like (Phase 7).
    window_subset: Optional[SampleSubset] = None
    # No view mode: it is a display option of each Annotations panel (D3).

    def __post_init__(self) -> None:
        set_ = object.__setattr__
        set_(self, "sample_query", _text(self.sample_query))
        set_(self, "gene", _text(self.gene))
        set_(self, "product", _text(self.product))
        molecule = _text(self.molecule).lower()
        set_(self, "molecule", molecule if molecule in MOLECULE_VALUES else "")
        set_(self, "predictors", frozenset(_unique_texts(self.predictors)))
        set_(self, "tags", frozenset(_unique_texts(self.tags)))
        set_(self, "hide_hypothetical", bool(self.hide_hypothetical))
        set_(self, "hide_empty_contigs", bool(self.hide_empty_contigs))
        set_(self, "metadata", _normalize_metadata(self.metadata))
        for name in ("subset", "window_subset"):
            subset = getattr(self, name)
            if subset is not None and not subset.samples:
                set_(self, name, None)

    # Construction

    @classmethod
    def from_query(
        cls,
        *,
        sample_query: Optional[str] = None,
        meta_filters: Optional[Mapping[str, Sequence[str]]] = None,
        subset_samples: Sequence[str] = (),
        ann_filters: Optional[Mapping[str, Any]] = None,
        hide_empty_contigs: bool = False,
    ) -> "FilterState":
        """From the pieces of the former refresh snapshot (tests, scripts)."""
        ann = dict(ann_filters or {})
        return cls(
            sample_query=_text(sample_query),
            gene=_text(ann.get("gene_text")),
            product=_text(ann.get("product_text")),
            molecule=_text(ann.get("molecule")),
            predictors=frozenset(ann.get("predictors") or ()),
            tags=frozenset(ann.get("tags") or ()),
            hide_hypothetical=bool(ann.get("hide_hypothetical")),
            hide_empty_contigs=bool(hide_empty_contigs),
            metadata=_normalize_metadata(meta_filters or {}),
            subset=SampleSubset(tuple(subset_samples)) if subset_samples else None,
        )

    @classmethod
    def from_widgets(cls, host: Any) -> "FilterState":
        """What the filter rail widgets of ``host`` show now.

        Duck typed (no Qt import): every widget is optional, so the small test
        harnesses that build only a few rail widgets can use it too.
        """

        def line_text(name: str) -> str:
            widget = getattr(host, name, None)
            return _text(widget.text()) if hasattr(widget, "text") else ""

        def selected(name: str) -> Tuple[str, ...]:
            method = getattr(host, name, None)
            if not callable(method):
                return ()
            try:
                return tuple(method())
            except Exception:
                return ()

        molecule = ""
        box = getattr(host, "molecule_filter_box", None)
        if hasattr(box, "currentData"):
            data = box.currentData()
            molecule = _text(data if data is not None else "").lower()
        # The subset chip is not a widget: it is the context's draft subset
        # (set through ctx.set_subset).
        draft = getattr(getattr(host, "ctx", None), "draft", None)
        subset = getattr(draft, "subset", None)
        return cls(
            sample_query=line_text("sample_query"),
            gene=line_text("gene_query"),
            product=line_text("product_query"),
            molecule=molecule,
            predictors=frozenset(selected("_selected_predictors")),
            tags=frozenset(selected("_selected_tags")),
            # Window filters shown in the Annotations panel, not in the rail:
            # the draft holds them (a panel may be closed).
            hide_hypothetical=bool(getattr(draft, "hide_hypothetical", False)),
            hide_empty_contigs=bool(getattr(draft, "hide_empty_contigs", False)),
            metadata=_normalize_metadata(getattr(host, "metadata_filters", None) or {}),
            subset=subset,
            window_subset=getattr(host, "_window_subset", None),
        )

    def replace(self, **changes: Any) -> "FilterState":
        return dataclasses.replace(self, **changes)

    # Reading

    def has_annotation_filter(self) -> bool:
        """Gene, product, molecule, predictors, tags or hide hypothetical (not the sample search)."""
        return bool(
            self.gene or self.product or self.molecule or self.predictors or self.tags or self.hide_hypothetical
        )

    def ann_filters(self) -> Dict[str, Any]:
        """The annotation filters in the engine's dict form."""
        return {
            "gene_text": self.gene or None,
            "product_text": self.product or None,
            "molecule": self.molecule,
            "predictors": sorted(self.predictors),
            "tags": sorted(self.tags),
            "hide_hypothetical": self.hide_hypothetical,
        }

    def metadata_filters(self) -> Dict[str, List[str]]:
        return {column: list(values) for column, values in self.metadata}

    def key(self) -> str:
        """Stable text form of everything that changes the result (subset labels do not)."""
        return json.dumps(
            {
                "sample_query": self.sample_query.lower(),
                "gene": self.gene.lower(),
                "product": self.product.lower(),
                "molecule": self.molecule,
                "predictors": sorted(p.lower() for p in self.predictors),
                "tags": sorted(t.upper() for t in self.tags),
                "hide_hypothetical": self.hide_hypothetical,
                "hide_empty_contigs": self.hide_empty_contigs,
                "metadata": [[column, list(values)] for column, values in self.metadata],
                "subset": list(self.subset.samples) if self.subset else None,
                "window_subset": list(self.window_subset.samples) if self.window_subset else None,
            },
            sort_keys=True,
            ensure_ascii=False,
        )


@dataclass(frozen=True, eq=False)
class AppliedResult:
    """The result of one FilterState. Immutable, so the engine's cache and every
    window share one object; the contig rows are read-only views."""

    # THE sample set of the window, in metadata order: the samples passing the
    # sample search, metadata filters and subsets, narrowed to the samples with
    # hits when an annotation filter is active. Every panel uses it: the sample
    # tree's rows, the map, the pan panels and the status bar count.
    samples: Tuple[str, ...]
    # The samples passing the sample search, metadata filters and subsets,
    # before the annotation narrowing (the annotations of their contigs are
    # the "base" rows the Contig table reads).
    allowed_samples: Tuple[str, ...]
    filtered_ann: pl.DataFrame
    metadata: pl.DataFrame  # metadata rows of samples (the tree's rows)
    contig_records: Mapping[str, Tuple[Mapping[str, Any], ...]]
    filtered_counts: Mapping[str, int]
    stats_text: str
    # The sample search matched no sample, but the same text is in gene or
    # product fields: the window shows the one-time D3 hint.
    sample_query_matches_annotations: bool = False

    def __post_init__(self) -> None:
        set_ = object.__setattr__
        set_(self, "samples", tuple(self.samples))
        set_(self, "allowed_samples", tuple(self.allowed_samples))
        set_(
            self,
            "contig_records",
            MappingProxyType({sample: tuple(rows) for sample, rows in dict(self.contig_records).items()}),
        )
        set_(self, "filtered_counts", MappingProxyType(dict(self.filtered_counts)))


@dataclass(frozen=True)
class HoverState:
    """What the pointer is over in one explore window: one or more samples.

    ``origin`` is the uid of the panel (or ``rail``, or ``info:<kind>:<key>``
    for a card) that produced it; the bus never delivers it back there.
    ``kind`` says what was hovered ("sample", "clade", "plasmid", ...).
    """

    samples: Tuple[str, ...] = ()
    origin: Optional[str] = None
    kind: str = "none"

    def __post_init__(self) -> None:
        samples = _unique_texts(self.samples)
        object.__setattr__(self, "samples", samples)
        object.__setattr__(self, "kind", _text(self.kind) or ("sample" if samples else "none"))
        if not samples:
            object.__setattr__(self, "kind", "none")

    @property
    def sample(self) -> Optional[str]:
        """The one hovered sample, or None when none or several are hovered."""
        return self.samples[0] if len(self.samples) == 1 else None

    def same_samples(self, other: "HoverState") -> bool:
        return self.samples == other.samples


@dataclass(frozen=True)
class Selection:
    """The window's shared selection: a sample's contig (or a plasmid's member contigs)."""

    sample: Optional[str] = None
    contig: Optional[str] = None
    contigs: Tuple[str, ...] = ()
    contig_label: Optional[str] = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "contigs", tuple(self.contigs or ()))
