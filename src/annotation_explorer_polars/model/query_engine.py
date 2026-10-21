"""QueryEngine: the filter queries over one RunData (window plan 4.3, Phase 2).

Pure and thread safe: ``compute`` runs on a pool thread and touches nothing
but the RunData and the engine's own caches (guarded by one lock). One engine
per RunData; several windows with the same filters share its caches.

``compute`` takes a ``FilterState`` and returns an immutable ``AppliedResult``
(Phase 3). The view mode of the annotation table is not part of the query
(decision D3): the sample search always filters samples, and the window's one
sample set narrows to the samples with hits only when an annotation filter is
active. Recorded by tests/test_refresh_parity.py. No Qt module is loaded.
"""
from __future__ import annotations

import math
import re
import threading
from collections import OrderedDict
from types import MappingProxyType
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import pandas as pd
import polars as pl

from ..constants import MISSING_FILTER_TOKEN, PLASMID_TRUE_VALUES
from .filter_state import AppliedResult, FilterState
from .run_data import RunData


COUNTRY_COLUMN_CANDIDATES: Tuple[str, ...] = ("Country", "_country")
LATITUDE_COLUMN_CANDIDATES: Tuple[str, ...] = ("Latitude", "latitude", "Lat", "lat")
LONGITUDE_COLUMN_CANDIDATES: Tuple[str, ...] = ("Longitude", "longitude", "Lon", "lon", "Long", "long")

_METADATA_CACHE_SIZE = 12
_ANNOTATION_CACHE_SIZE = 8
_RESULT_CACHE_SIZE = 8
_MAP_CACHE_SIZE = 8
# The per-country sample names, kept in the cached map payload but not
# returned by map_payload (map_country_samples reads them).
_COUNTRY_SAMPLES_KEY = "_countrySamples"


def resolve_metadata_column(columns: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    normalized = {str(column).strip().casefold(): str(column) for column in columns if str(column).strip()}
    for candidate in candidates:
        key = str(candidate).strip().casefold()
        if key in normalized:
            return normalized[key]
    return None


def normalize_sample_values(samples: Sequence[str]) -> List[str]:
    out: List[str] = []
    seen: set[str] = set()
    for raw in samples:
        val = str(raw).strip()
        if not val or val in seen:
            continue
        seen.add(val)
        out.append(val)
    return out


def collect_sample_list(df: pl.DataFrame) -> List[str]:
    """The distinct samples of ``df`` in row order (first appearance)."""
    if df.is_empty() or "Sample" not in df.columns:
        return []
    vals = (
        df.select(pl.col("Sample").cast(pl.Utf8).drop_nulls().unique(maintain_order=True))
        .to_series(0)
        .to_list()
    )
    return normalize_sample_values(vals)


def _filter_contigs_by_molecule(contigs_pl: pl.DataFrame, filters: Dict[str, Any]) -> pl.DataFrame:
    molecule = str(filters.get("molecule") or "").strip().lower()
    if molecule not in {"plasmid", "chromosome"} or contigs_pl.is_empty():
        return contigs_pl
    if "MobSuiteIsPlasmid" not in contigs_pl.columns:
        return contigs_pl.head(0) if molecule == "plasmid" else contigs_pl
    is_plasmid = (
        pl.col("MobSuiteIsPlasmid")
        .cast(pl.Utf8, strict=False)
        .fill_null("")
        .str.strip_chars()
        .str.to_lowercase()
        .is_in(list(PLASMID_TRUE_VALUES))
    )
    if molecule == "plasmid":
        return contigs_pl.lazy().filter(is_plasmid).collect()
    return contigs_pl.lazy().filter(~is_plasmid).collect()


def _build_contig_record_map(contigs_pl: pl.DataFrame) -> Dict[str, List[Mapping[str, Any]]]:
    """Contig rows per sample. Each row is a read-only view: the rows are shared
    by the result cache and every window the result is handed to."""
    contig_records: Dict[str, List[Mapping[str, Any]]] = {}
    for row in contigs_pl.iter_rows(named=True):
        sample_name = str(row.get("Sample") or "").strip()
        if not sample_name:
            continue
        row.pop("_contig_num", None)
        contig_records.setdefault(sample_name, []).append(MappingProxyType(row))
    return contig_records


def _build_filtered_counts(filtered_ann_pl: pl.DataFrame) -> Dict[str, int]:
    if filtered_ann_pl.is_empty():
        return {}
    counts: Dict[str, int] = {}
    for row in (
        filtered_ann_pl.lazy()
        .group_by("Sample")
        .agg(pl.len().alias("n_hits"))
        .collect()
        .iter_rows(named=True)
    ):
        sample_name = str(row.get("Sample") or "").strip()
        if sample_name:
            counts[sample_name] = int(row.get("n_hits") or 0)
    return counts


class QueryEngine:
    def __init__(self, run: RunData) -> None:
        self.run = run
        self._lock = threading.Lock()
        self._metadata_cache: "OrderedDict[object, pl.DataFrame]" = OrderedDict()
        self._annotation_cache: "OrderedDict[object, pl.DataFrame]" = OrderedDict()
        self._result_cache: "OrderedDict[str, AppliedResult]" = OrderedDict()
        self._map_cache: "OrderedDict[Tuple[str, ...], dict]" = OrderedDict()

    # Caches

    def _cache_get(self, cache: "OrderedDict[object, Any]", key: object) -> Optional[Any]:
        with self._lock:
            value = cache.get(key)
            if value is None:
                return None
            cache.move_to_end(key)
            return value

    def _cache_put(self, cache: "OrderedDict[object, Any]", key: object, value: Any, *, max_entries: int) -> Any:
        with self._lock:
            cache[key] = value
            cache.move_to_end(key)
            while len(cache) > max_entries:
                cache.popitem(last=False)
        return value

    # Empty frames with the public columns

    def empty_annotations(self) -> pl.DataFrame:
        return self.run.ann_pl.select(list(self.run.ann_public_columns)).head(0)

    def empty_metadata(self) -> pl.DataFrame:
        return self.run.meta_pl.select(list(self.run.meta_public_columns)).head(0)

    # Queries

    @staticmethod
    def _normalized_filter_key(column_filters: Dict[str, Sequence[str]]) -> Tuple[Tuple[str, Tuple[str, ...]], ...]:
        items: List[Tuple[str, Tuple[str, ...]]] = []
        for key in sorted(column_filters):
            values = tuple(sorted(str(v) for v in column_filters.get(key, []) if v is not None))
            items.append((str(key), values))
        return tuple(items)

    def metadata_subset(
        self,
        sample_query: Optional[str],
        column_filters: Dict[str, Sequence[str]],
    ) -> pl.DataFrame:
        """Metadata rows (public columns) matching the sample search and the metadata filters."""
        cache_key = (str(sample_query or "").strip().lower() or None, self._normalized_filter_key(column_filters))
        cached = self._cache_get(self._metadata_cache, cache_key)
        if isinstance(cached, pl.DataFrame):
            return cached

        meta_pl = self.run.meta_pl
        lf = meta_pl.lazy()
        q = str(sample_query or "").strip().lower()
        if q:
            lf = lf.filter(pl.col("__search_text").str.contains(re.escape(q)))

        for col, raw_values in column_filters.items():
            if col not in meta_pl.columns or not raw_values:
                continue
            values = [str(v) for v in raw_values if v is not None]
            include_missing = MISSING_FILTER_TOKEN in values
            values = [v for v in values if v != MISSING_FILTER_TOKEN]
            text_col = pl.col(col).cast(pl.Utf8).fill_null("").str.strip_chars()
            if include_missing and values:
                lf = lf.filter(pl.col(col).is_null() | (text_col == "") | text_col.is_in(values))
            elif include_missing:
                lf = lf.filter(pl.col(col).is_null() | (text_col == ""))
            elif values:
                lf = lf.filter(text_col.is_in(values))

        result = lf.select(list(self.run.meta_public_columns)).collect()
        return self._cache_put(self._metadata_cache, cache_key, result, max_entries=_METADATA_CACHE_SIZE)

    def filter_annotations(
        self,
        filters: Dict[str, Any],
        allowed_samples: Sequence[str],
    ) -> pl.DataFrame:
        """Annotation rows (public columns) of the allowed samples that pass the annotation filters."""
        allowed = normalize_sample_values(allowed_samples)
        if len(allowed_samples) == 0 or not allowed:
            return self.empty_annotations()

        cache_key = (
            tuple(allowed),
            str(filters.get("gene_text") or "").strip().lower() or None,
            str(filters.get("product_text") or "").strip().lower() or None,
            str(filters.get("molecule") or "").strip().lower() or None,
            tuple(sorted(str(v).strip().lower() for v in (filters.get("predictors") or []) if str(v).strip())),
            tuple(sorted(str(v).strip().upper() for v in (filters.get("tags") or []) if str(v).strip())),
            bool(filters.get("hide_hypothetical")),
        )
        cached = self._cache_get(self._annotation_cache, cache_key)
        if isinstance(cached, pl.DataFrame):
            return cached

        lf = self.run.ann_pl.lazy().filter(pl.col("__sample_norm").is_in(allowed))

        predictor_vals = [str(p).strip() for p in (filters.get("predictors") or []) if str(p).strip()]
        if predictor_vals:
            pred_expr = pl.lit(False)
            pred_col = pl.col("__predictors_norm")
            for predictor in predictor_vals:
                pat = r"(^|,\s*)" + re.escape(predictor.lower()) + r"(,|$)"
                pred_expr = pred_expr | pred_col.str.contains(pat)
            lf = lf.filter(pred_expr)

        tag_vals = [str(t).upper().strip() for t in (filters.get("tags") or []) if str(t).strip()]
        if tag_vals:
            tag_pat = "[" + "".join(re.escape(t) for t in tag_vals) + "]"
            lf = lf.filter(pl.col("__tags_upper").str.contains(tag_pat))

        gene_text = str(filters.get("gene_text") or "").strip()
        if gene_text:
            gene_pat = re.escape(gene_text.lower())
            lf = lf.filter(
                pl.any_horizontal(
                    [
                        pl.col("__gene_norm").str.contains(gene_pat),
                        pl.col("__marker_search_norm").str.contains(gene_pat),
                        pl.col("__marker_display_norm").str.contains(gene_pat),
                        pl.col("__predictors_norm").str.contains(gene_pat),
                        pl.col("__method_norm").str.contains(gene_pat),
                    ]
                )
            )

        product_text = str(filters.get("product_text") or "").strip()
        if product_text:
            lf = lf.filter(pl.col("__product_norm").str.contains(re.escape(product_text.lower())))

        if filters.get("hide_hypothetical"):
            lf = lf.filter(~pl.col("__is_hypothetical"))

        molecule = str(filters.get("molecule") or "").strip().lower()
        if molecule == "plasmid":
            lf = lf.filter(pl.col("__is_mobsuite_plasmid"))
        elif molecule == "chromosome":
            lf = lf.filter(~pl.col("__is_mobsuite_plasmid"))

        result = lf.select(list(self.run.ann_public_columns)).collect()
        return self._cache_put(self._annotation_cache, cache_key, result, max_entries=_ANNOTATION_CACHE_SIZE)

    def compute(self, filters: FilterState) -> AppliedResult:
        """The applied result of one FilterState (cached; the result is immutable,
        so every caller with the same filters gets the same object)."""
        cache_key = filters.key()
        cached = self._cache_get(self._result_cache, cache_key)
        if isinstance(cached, AppliedResult):
            return cached
        result = self._compute(filters)
        return self._cache_put(self._result_cache, cache_key, result, max_entries=_RESULT_CACHE_SIZE)

    def _compute(self, filters: FilterState) -> AppliedResult:
        # D3: the sample search always filters samples by name and metadata.
        meta_subset_pl = self.metadata_subset(filters.sample_query or None, filters.metadata_filters())
        allowed_samples = collect_sample_list(meta_subset_pl)
        for subset in (filters.subset, filters.window_subset):
            subset_samples = normalize_sample_values(subset.samples) if subset is not None else []
            if not subset_samples:
                continue
            subset_allowed = set(subset_samples)
            allowed_samples = [sample for sample in allowed_samples if sample in subset_allowed]
            if not meta_subset_pl.is_empty():
                meta_subset_pl = (
                    meta_subset_pl.lazy()
                    .filter(pl.col("Sample").cast(pl.Utf8).is_in(allowed_samples))
                    .collect()
                )

        ann_filters = filters.ann_filters()
        filtered_ann_pl = self.filter_annotations(ann_filters, allowed_samples)

        if not filtered_ann_pl.is_empty():
            contig_hits_pl = (
                filtered_ann_pl.lazy()
                .group_by(["Sample", "Contig"])
                .agg(pl.len().alias("n_hits"))
                .collect()
            )
        else:
            contig_hits_pl = pl.DataFrame(
                schema={"Sample": pl.Utf8, "Contig": pl.Utf8, "n_hits": pl.Int64}
            )

        contig_totals_pl = self.run.contig_totals_pl
        if allowed_samples:
            contigs_pl = (
                contig_totals_pl.lazy()
                .filter(pl.col("Sample").cast(pl.Utf8).is_in(allowed_samples))
                .join(contig_hits_pl.lazy(), on=["Sample", "Contig"], how="left")
                .with_columns(
                    [
                        pl.col("n_hits").fill_null(0).cast(pl.Int64),
                        pl.col("n_total").fill_null(0).cast(pl.Int64),
                        pl.col("Contig").cast(pl.Float64, strict=False).alias("_contig_num"),
                    ]
                )
                .sort(["Sample", "_contig_num", "Contig"], nulls_last=True)
                .collect()
            )
        else:
            contigs_pl = contig_totals_pl.head(0).with_columns(
                [
                    pl.lit(0).cast(pl.Int64).alias("n_hits"),
                    pl.lit(None).alias("_contig_num"),
                ]
            )

        contigs_pl = _filter_contigs_by_molecule(contigs_pl, ann_filters)
        # Only the user's "Hide contigs without hits" drops contigs (D3).
        if filters.hide_empty_contigs:
            contigs_pl = contigs_pl.filter(pl.col("n_hits") > 0)

        contig_records = _build_contig_record_map(contigs_pl)
        filtered_counts = _build_filtered_counts(filtered_ann_pl)
        # The window's one sample set: narrowed to the samples with hits only
        # by an annotation filter (the sample search narrows through the
        # metadata filter instead, D3). Metadata order.
        samples = allowed_samples
        meta_for_samples_pl = meta_subset_pl
        if filters.has_annotation_filter():
            with_hits = set(collect_sample_list(filtered_ann_pl))
            samples = [sample for sample in allowed_samples if sample in with_hits]
            if len(samples) != len(allowed_samples) and not meta_subset_pl.is_empty():
                meta_for_samples_pl = (
                    meta_subset_pl.lazy()
                    .filter(pl.col("Sample").cast(pl.Utf8).is_in(samples))
                    .collect()
                )

        return AppliedResult(
            samples=tuple(samples),
            allowed_samples=tuple(allowed_samples),
            filtered_ann=filtered_ann_pl,
            metadata=meta_for_samples_pl,
            contig_records=contig_records,
            filtered_counts=filtered_counts,
            # The status bar counts the window's one sample set (plan 4.4).
            stats_text=f"Samples: {len(samples)} | Annotations after filters: {filtered_ann_pl.height}",
            sample_query_matches_annotations=self._sample_query_is_an_annotation_search(filters),
        )

    def _sample_query_is_an_annotation_search(self, filters: FilterState) -> bool:
        """The sample search alone matches no sample, but the same text is in the
        fields that "Gene/annotation contains" or "Product contains" search."""
        query = filters.sample_query.lower()
        if not query or not self.metadata_subset(query, {}).is_empty():
            return False
        pattern = re.escape(query)
        columns = [
            column
            for column in (
                "__gene_norm",
                "__marker_search_norm",
                "__marker_display_norm",
                "__predictors_norm",
                "__method_norm",
                "__product_norm",
            )
            if column in self.run.ann_pl.columns
        ]
        if not columns:
            return False
        matches = (
            self.run.ann_pl.lazy()
            .filter(pl.any_horizontal([pl.col(column).str.contains(pattern) for column in columns]))
            .head(1)
            .collect()
        )
        return not matches.is_empty()

    def contig_annotations(
        self,
        *,
        allowed_samples: Optional[Sequence[str]] = None,
        filtered: Optional[pl.DataFrame] = None,
        sample: Optional[str] = None,
        contig: Optional[str] = None,
        contigs: Optional[Sequence[str]] = None,
    ) -> pd.DataFrame:
        """Annotation rows of one sample's contig(s) as pandas (public columns).

        The rows come from ``filtered`` (an applied result's filtered
        annotations) when it is given, else from every annotation of
        ``allowed_samples``. Nothing to read gives an empty frame.
        """
        public = list(self.run.ann_public_columns)
        if filtered is not None:
            if filtered.is_empty():
                return pd.DataFrame(columns=public)
            lf = filtered.lazy()
        else:
            if not allowed_samples:
                return pd.DataFrame(columns=public)
            lf = self.run.ann_pl.lazy().filter(pl.col("__sample_norm").is_in(list(allowed_samples)))
        if sample is not None:
            lf = lf.filter(pl.col("Sample").cast(pl.Utf8) == str(sample).strip())
        contig_values = [str(value).strip() for value in (contigs or []) if str(value).strip()]
        if contig_values:
            lf = lf.filter(pl.col("Contig").cast(pl.Utf8).is_in(contig_values))
        elif contig is not None:
            lf = lf.filter(pl.col("Contig").cast(pl.Utf8) == str(contig).strip())
        return lf.select(public).collect().to_pandas()

    def result_annotations(
        self,
        result: Optional[AppliedResult],
        *,
        hits: bool,
        sample: Optional[str] = None,
        contig: Optional[str] = None,
        contigs: Optional[Sequence[str]] = None,
    ) -> pd.DataFrame:
        """``contig_annotations`` of an applied result: its hits (the filtered
        annotations) or every annotation of its allowed samples. Before the
        first result (None) the frame is empty."""
        if result is None:
            return pd.DataFrame(columns=list(self.run.ann_public_columns))
        if hits:
            return self.contig_annotations(filtered=result.filtered_ann, sample=sample, contig=contig, contigs=contigs)
        return self.contig_annotations(
            allowed_samples=list(result.allowed_samples), sample=sample, contig=contig, contigs=contigs
        )

    def map_payload(self, samples: Sequence[str]) -> dict:
        """Map points (samples with coordinates) and per-country counts (the rest) for these samples.

        Cached by the sample tuple, so several maps of one window, or of
        windows with the same filters, compute it once. Each call returns a new dict; its lists are shared, read only."""
        key = tuple(str(s) for s in samples)
        cached = self._cache_get(self._map_cache, key)
        if cached is None:
            cached = self._cache_put(self._map_cache, key, self._compute_map_payload(key), max_entries=_MAP_CACHE_SIZE)
        payload = dict(cached)
        payload.pop(_COUNTRY_SAMPLES_KEY, None)
        return payload

    def map_country_samples(self, samples: Sequence[str]) -> Dict[str, List[str]]:
        """The samples of each country the map counts (country name as in the
        metadata, samples in metadata order), from the same cached computation
        as ``map_payload``: samples with coordinates and a country, and samples
        with a country only. Read only."""
        key = tuple(str(s) for s in samples)
        cached = self._cache_get(self._map_cache, key)
        if cached is None:
            cached = self._cache_put(self._map_cache, key, self._compute_map_payload(key), max_entries=_MAP_CACHE_SIZE)
        return cached.get(_COUNTRY_SAMPLES_KEY) or {}

    def _compute_map_payload(self, samples: Sequence[str]) -> dict:
        meta_pl = self.run.meta_pl
        samples = list(samples)
        if not samples or "Sample" not in meta_pl.columns:
            return {"points": [], "countries": {}, "maxCountryCount": 1}
        meta_columns = [str(column) for column in meta_pl.columns]
        lat_col = resolve_metadata_column(meta_columns, LATITUDE_COLUMN_CANDIDATES)
        lon_col = resolve_metadata_column(meta_columns, LONGITUDE_COLUMN_CANDIDATES)
        country_col = resolve_metadata_column(meta_columns, COUNTRY_COLUMN_CANDIDATES)
        meta = (
            meta_pl.lazy()
            .with_columns(
                [
                    pl.col("Sample").cast(pl.Utf8).str.strip_chars().alias("Sample"),
                    (
                        pl.col(lat_col).cast(pl.Float64, strict=False)
                        if lat_col
                        else pl.lit(None, dtype=pl.Float64)
                    ).alias("_lat"),
                    (
                        pl.col(lon_col).cast(pl.Float64, strict=False)
                        if lon_col
                        else pl.lit(None, dtype=pl.Float64)
                    ).alias("_lon"),
                    (
                        pl.col(country_col).cast(pl.Utf8).fill_null("").str.strip_chars()
                        if country_col
                        else pl.lit("", dtype=pl.Utf8)
                    ).alias("_country"),
                ]
            )
            .filter(pl.col("Sample").is_in(samples))
            .unique(subset=["Sample"], keep="first", maintain_order=True)
            .select(["Sample", "_lat", "_lon", "_country"])
            .collect()
        )
        if meta.is_empty():
            return {"points": [], "countries": {}, "maxCountryCount": 1, "total": len(samples)}

        points: List[Dict[str, Any]] = []
        countries: Dict[str, int] = {}
        country_samples: Dict[str, List[str]] = {}
        for row in meta.iter_rows(named=True):
            sample = str(row.get("Sample") or "").strip()
            if not sample:
                continue
            lat = row.get("_lat")
            lon = row.get("_lon")
            has_coords = (
                isinstance(lat, (float, int))
                and isinstance(lon, (float, int))
                and not (isinstance(lat, float) and math.isnan(lat))
                and not (isinstance(lon, float) and math.isnan(lon))
            )
            country_name = str(row.get("_country") or "").strip()
            if country_name:
                country_samples.setdefault(country_name, []).append(sample)
            if has_coords:
                points.append(
                    {
                        "sample": sample,
                        "lat": float(lat),
                        "lon": float(lon),
                        "count": 1,
                        "country": str(row.get("_country") or ""),
                    }
                )
            else:
                country = str(row.get("_country") or "").strip()
                if country:
                    countries[country] = countries.get(country, 0) + 1

        max_count = max(countries.values(), default=1)
        # total also counts samples with no location (in neither points nor countries)
        return {"points": points, "countries": countries, "maxCountryCount": max_count, "total": len(samples),
                _COUNTRY_SAMPLES_KEY: country_samples}

    # Info cards (the cards filter these shared frames themselves)

    def sample_bundle(self, sample: str) -> Dict[str, Any]:
        """Keyword arguments of a sample card for ``sample``."""
        run = self.run
        return {
            "sample": sample,
            "ann_df": run.ann_df,
            "meta_df": run.meta_df,
            "marker_hits_df": run.marker_hits_df,
            "reference_hits_df": run.reference_hits_df,
            "output_dir": run.output_dir,
        }

    def annotation_bundle(self) -> Dict[str, Any]:
        """Run data keyword arguments of a gene card (its row is passed separately).
        The neighborhood viewer it opens takes the rest from its owner window."""
        run = self.run
        return {
            "ann_df": run.ann_df,
            "meta_df": run.meta_df,
            "marker_hits_df": run.marker_hits_df,
        }
