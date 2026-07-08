from __future__ import annotations

from collections import OrderedDict
import re
import threading
from typing import Any, Dict, List, Optional, Sequence, Tuple

import pandas as pd
import polars as pl

from .constants import MARKER_SEARCH_TEXT_COLUMN, MISSING_FILTER_TOKEN, PLASMID_TRUE_VALUES
from .refresh_worker import RefreshPayload, RefreshSnapshot


class PolarsQueriesMixin:
    @staticmethod
    def _norm_text_expr(column: str) -> pl.Expr:
        return pl.col(column).cast(pl.Utf8, strict=False).fill_null("").str.strip_chars()


    @staticmethod
    def _cache_get(
        cache: "OrderedDict[object, object]",
        lock: threading.Lock,
        key: object,
    ) -> Optional[object]:
        with lock:
            value = cache.get(key)
            if value is None:
                return None
            cache.move_to_end(key)
            return value


    @staticmethod
    def _cache_put(
        cache: "OrderedDict[object, object]",
        lock: threading.Lock,
        key: object,
        value: object,
        *,
        max_entries: int,
    ) -> object:
        with lock:
            cache[key] = value
            cache.move_to_end(key)
            while len(cache) > max_entries:
                cache.popitem(last=False)
        return value


    def _prepare_metadata_polars(self, df: pl.DataFrame) -> pl.DataFrame:
        search_parts = [self._norm_text_expr(str(col)).str.to_lowercase() for col in df.columns]
        return (
            df.lazy()
            .with_columns(
                [
                    self._norm_text_expr("Sample").alias("__sample_norm"),
                    (
                        pl.concat_str(search_parts, separator=" ")
                        if search_parts
                        else pl.lit("", dtype=pl.Utf8)
                    ).alias("__search_text"),
                ]
            )
            .collect()
        )


    def _prepare_annotation_polars(self, df: pl.DataFrame) -> pl.DataFrame:
        def _expr_or_empty(column: str) -> pl.Expr:
            if column in df.columns:
                return self._norm_text_expr(column)
            return pl.lit("", dtype=pl.Utf8)

        return (
            df.lazy()
            .with_columns(
                [
                    _expr_or_empty("Sample").alias("__sample_norm"),
                    _expr_or_empty("Contig").alias("__contig_norm"),
                    _expr_or_empty("Gene").str.to_lowercase().alias("__gene_norm"),
                    _expr_or_empty("Product").str.to_lowercase().alias("__product_norm"),
                    _expr_or_empty("Predictors").str.to_lowercase().alias("__predictors_norm"),
                    _expr_or_empty("Method").str.to_lowercase().alias("__method_norm"),
                    _expr_or_empty("MarkerDisplay").str.to_lowercase().alias("__marker_display_norm"),
                    _expr_or_empty(MARKER_SEARCH_TEXT_COLUMN).str.to_lowercase().alias("__marker_search_norm"),
                    _expr_or_empty("Hypothetical").str.to_lowercase().alias("__hypothetical_norm"),
                    _expr_or_empty("MobSuiteIsPlasmid").str.to_lowercase().alias("__mobsuite_plasmid_norm"),
                    _expr_or_empty("Tags").str.to_uppercase().alias("__tags_upper"),
                ]
            )
            .with_columns(
                [
                    pl.concat_str(
                        [
                            pl.col("__sample_norm").str.to_lowercase(),
                            pl.col("__contig_norm").str.to_lowercase(),
                            pl.col("__gene_norm"),
                            pl.col("__product_norm"),
                            pl.col("__predictors_norm"),
                            pl.col("__method_norm"),
                            pl.col("__marker_display_norm"),
                            pl.col("__marker_search_norm"),
                            pl.col("__tags_upper").str.to_lowercase(),
                        ],
                        separator=" ",
                    ).alias("__search_text"),
                    pl.col("__hypothetical_norm").is_in(["true", "1", "yes"]).alias("__is_hypothetical"),
                    pl.col("__mobsuite_plasmid_norm").is_in(list(PLASMID_TRUE_VALUES)).alias("__is_mobsuite_plasmid"),
                ]
            )
            .collect()
        )


    def _empty_ann_public_pl(self) -> pl.DataFrame:
        return self._ann_pl.select(self._ann_public_columns).head(0)


    def _empty_meta_public_pl(self) -> pl.DataFrame:
        return self._meta_pl.select(self._meta_public_columns).head(0)


    def _build_annotation_hit_counts_polars(self) -> pl.DataFrame:
        cols = ["Sample", "Contig", "DB", "Hits"]
        frames: List[pl.DataFrame] = []
        if not self.marker_hits_df.empty:
            frames.append(pl.from_pandas(self.marker_hits_df[["Sample", "Contig", "DB"]]))
        if not self.reference_hits_df.empty:
            frames.append(pl.from_pandas(self.reference_hits_df[["Sample", "Contig", "DB"]]))
        if not frames:
            return pl.DataFrame(schema={"Sample": pl.Utf8, "Contig": pl.Utf8, "DB": pl.Utf8, "Hits": pl.Int64})

        combined = pl.concat(frames, how="vertical_relaxed")
        counts = (
            combined.lazy()
            .with_columns(
                [
                    pl.col("Sample").cast(pl.Utf8).str.strip_chars().alias("Sample"),
                    pl.col("Contig").cast(pl.Utf8).str.strip_chars().alias("Contig"),
                    pl.col("DB").cast(pl.Utf8).str.strip_chars().alias("DB"),
                ]
            )
            .filter(
                (pl.col("Sample").is_not_null())
                & (pl.col("Sample") != "")
                & (pl.col("Contig").is_not_null())
                & (pl.col("Contig") != "")
                & (pl.col("DB").is_not_null())
                & (pl.col("DB") != "")
            )
            .group_by(["Sample", "Contig", "DB"])
            .agg(pl.len().alias("Hits"))
            .sort(["Sample", "Contig", "DB"], nulls_last=True)
            .collect()
        )
        if counts.is_empty():
            return pl.DataFrame(schema={"Sample": pl.Utf8, "Contig": pl.Utf8, "DB": pl.Utf8, "Hits": pl.Int64})
        return counts.select(cols)


    def _normalize_sample_values(self, samples: Sequence[str]) -> List[str]:
        out: List[str] = []
        seen: set[str] = set()
        for raw in samples:
            val = str(raw).strip()
            if not val or val in seen:
                continue
            seen.add(val)
            out.append(val)
        return out


    def _collect_sample_list(self, df: pl.DataFrame) -> List[str]:
        if df.is_empty() or "Sample" not in df.columns:
            return []
        vals = (
            df.select(pl.col("Sample").cast(pl.Utf8).drop_nulls().unique())
            .to_series(0)
            .to_list()
        )
        return self._normalize_sample_values(vals)


    def _subset_annotations_pd(
        self,
        source: str,
        *,
        sample: Optional[str] = None,
        contig: Optional[str] = None,
        contigs: Optional[Sequence[str]] = None,
    ) -> pd.DataFrame:
        if source == "base":
            if not getattr(self, "_base_allowed_samples", None):
                return pd.DataFrame(columns=self._ann_public_columns)
            lf = self._ann_pl.lazy().filter(pl.col("__sample_norm").is_in(self._base_allowed_samples))
        else:
            df_pl = self._filtered_ann_pl
            if df_pl is None or df_pl.is_empty():
                return pd.DataFrame(columns=self._ann_public_columns)
            lf = df_pl.lazy()
        if sample is not None:
            lf = lf.filter(pl.col("Sample").cast(pl.Utf8) == str(sample).strip())
        contig_values = [
            str(value).strip()
            for value in (contigs or [])
            if str(value).strip()
        ]
        if contig_values:
            lf = lf.filter(pl.col("Contig").cast(pl.Utf8).is_in(contig_values))
        elif contig is not None:
            lf = lf.filter(pl.col("Contig").cast(pl.Utf8) == str(contig).strip())
        return lf.select(self._ann_public_columns).collect().to_pandas()


    @staticmethod
    def _normalized_filter_key(column_filters: Dict[str, Sequence[str]]) -> Tuple[Tuple[str, Tuple[str, ...]], ...]:
        items: List[Tuple[str, Tuple[str, ...]]] = []
        for key in sorted(column_filters):
            values = tuple(sorted(str(v) for v in column_filters.get(key, []) if v is not None))
            items.append((str(key), values))
        return tuple(items)


    def _metadata_cache_key(
        self,
        sample_query: Optional[str],
        column_filters: Dict[str, Sequence[str]],
    ) -> Tuple[Optional[str], Tuple[Tuple[str, Tuple[str, ...]], ...]]:
        query = str(sample_query or "").strip().lower() or None
        return query, self._normalized_filter_key(column_filters)


    def _annotation_cache_key(
        self,
        filters: Dict[str, Any],
        allowed_samples: Sequence[str],
        *,
        apply_global_query: bool,
        global_query: Optional[str],
    ) -> Tuple[Tuple[str, ...], Optional[str], Optional[str], Optional[str], Tuple[str, ...], Tuple[str, ...], bool, bool, Optional[str]]:
        return (
            tuple(self._normalize_sample_values(allowed_samples)),
            str(filters.get("gene_text") or "").strip().lower() or None,
            str(filters.get("product_text") or "").strip().lower() or None,
            str(filters.get("molecule") or "").strip().lower() or None,
            tuple(sorted(str(v).strip().lower() for v in (filters.get("predictors") or []) if str(v).strip())),
            tuple(sorted(str(v).strip().upper() for v in (filters.get("tags") or []) if str(v).strip())),
            bool(filters.get("hide_hypothetical")),
            bool(apply_global_query),
            str(global_query or "").strip().lower() or None,
        )


    def _metadata_subset_polars(
        self,
        sample_query: Optional[str],
        column_filters: Dict[str, Sequence[str]],
    ) -> pl.DataFrame:
        cache_key = self._metadata_cache_key(sample_query, column_filters)
        cached = self._cache_get(self._metadata_subset_cache, self._polars_cache_lock, cache_key)
        if isinstance(cached, pl.DataFrame):
            return cached

        lf = self._meta_pl.lazy()
        q = str(sample_query or "").strip().lower()
        if q:
            lf = lf.filter(pl.col("__search_text").str.contains(re.escape(q)))

        for col, raw_values in column_filters.items():
            if col not in self._meta_pl.columns or not raw_values:
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

        result = lf.select(self._meta_public_columns).collect()
        return self._cache_put(
            self._metadata_subset_cache,
            self._polars_cache_lock,
            cache_key,
            result,
            max_entries=12,
        )  # type: ignore[return-value]


    def _filter_annotations_polars(
        self,
        filters: Dict[str, Any],
        allowed_samples: Sequence[str],
        *,
        apply_global_query: bool = False,
        global_query: Optional[str] = None,
    ) -> pl.DataFrame:
        allowed = self._normalize_sample_values(allowed_samples)
        if len(allowed_samples) == 0:
            return self._empty_ann_public_pl()

        cache_key = self._annotation_cache_key(
            filters,
            allowed,
            apply_global_query=apply_global_query,
            global_query=global_query,
        )
        cached = self._cache_get(self._annotation_filter_cache, self._polars_cache_lock, cache_key)
        if isinstance(cached, pl.DataFrame):
            return cached

        lf = self._ann_pl.lazy()
        if allowed:
            lf = lf.filter(pl.col("__sample_norm").is_in(allowed))
        else:
            return self._empty_ann_public_pl()

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

        if apply_global_query:
            query = str(global_query or "").strip().lower()
            if query:
                lf = lf.filter(pl.col("__search_text").str.contains(re.escape(query)))

        result = lf.select(self._ann_public_columns).collect()
        return self._cache_put(
            self._annotation_filter_cache,
            self._polars_cache_lock,
            cache_key,
            result,
            max_entries=8,
        )  # type: ignore[return-value]


    @staticmethod
    def _snapshot_has_active_annotation_filter(snapshot: RefreshSnapshot) -> bool:
        filters = snapshot.ann_filters
        return any(
            [
                bool(filters.get("gene_text")),
                bool(filters.get("product_text")),
                bool(filters.get("molecule")),
                bool(filters.get("predictors")),
                bool(filters.get("tags")),
                bool(filters.get("hide_hypothetical")),
            ]
        )


    def _filter_contigs_by_molecule(
        self,
        contigs_pl: pl.DataFrame,
        filters: Dict[str, Any],
    ) -> pl.DataFrame:
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


    def _snapshot_should_show_global(self, snapshot: RefreshSnapshot) -> bool:
        search_text = str(snapshot.sample_query or "").strip()
        if search_text:
            return True
        return self._snapshot_has_active_annotation_filter(snapshot)


    @staticmethod
    def _build_contig_record_map(contigs_pl: pl.DataFrame) -> Dict[str, List[Dict[str, Any]]]:
        contig_records: Dict[str, List[Dict[str, Any]]] = {}
        for row in contigs_pl.iter_rows(named=True):
            sample_name = str(row.get("Sample") or "").strip()
            if not sample_name:
                continue
            row.pop("_contig_num", None)
            contig_records.setdefault(sample_name, []).append(row)
        return contig_records


    @staticmethod
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


    def _compute_refresh_payload(self, snapshot: RefreshSnapshot) -> RefreshPayload:
        meta_sample_q = None if snapshot.view_mode == "Global" else snapshot.sample_query
        meta_subset_pl = self._metadata_subset_polars(meta_sample_q, snapshot.meta_filters)
        allowed_samples = self._collect_sample_list(meta_subset_pl)
        subset_samples = self._normalize_sample_values(snapshot.sample_subset_samples)
        if subset_samples:
            subset_allowed = {sample for sample in subset_samples}
            allowed_samples = [sample for sample in allowed_samples if sample in subset_allowed]
            if not meta_subset_pl.is_empty():
                meta_subset_pl = (
                    meta_subset_pl.lazy()
                    .filter(pl.col("Sample").cast(pl.Utf8).is_in(allowed_samples))
                    .collect()
                )

        filtered_ann_pl = self._filter_annotations_polars(
            snapshot.ann_filters,
            allowed_samples,
            apply_global_query=snapshot.view_mode == "Global",
            global_query=snapshot.sample_query if snapshot.view_mode == "Global" else None,
        )

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

        if allowed_samples:
            contigs_pl = (
                self._contig_totals_pl.lazy()
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
            contigs_pl = self._contig_totals_pl.head(0).with_columns(
                [
                    pl.lit(0).cast(pl.Int64).alias("n_hits"),
                    pl.lit(None).alias("_contig_num"),
                ]
            )

        contigs_pl = self._filter_contigs_by_molecule(contigs_pl, snapshot.ann_filters)

        show_hits_only = snapshot.view_mode == "Global" and self._snapshot_should_show_global(snapshot)
        if show_hits_only:
            contigs_pl = contigs_pl.filter(pl.col("n_hits") > 0)
        elif snapshot.hide_empty_contigs:
            contigs_pl = contigs_pl.filter(pl.col("n_hits") > 0)

        contig_records = self._build_contig_record_map(contigs_pl)
        filtered_counts = self._build_filtered_counts(filtered_ann_pl)
        panel_allowed_samples = allowed_samples
        if str(snapshot.sample_query or "").strip() or self._snapshot_has_active_annotation_filter(snapshot):
            panel_allowed_samples = self._collect_sample_list(filtered_ann_pl)

        if show_hits_only:
            visible_samples = self._collect_sample_list(contigs_pl)
            if visible_samples:
                meta_for_tree_pl = (
                    meta_subset_pl.lazy()
                    .filter(pl.col("Sample").cast(pl.Utf8).is_in(visible_samples))
                    .collect()
                )
            else:
                meta_for_tree_pl = self._empty_meta_public_pl()
            tree_allowed_samples = visible_samples
        else:
            meta_for_tree_pl = meta_subset_pl
            tree_allowed_samples = allowed_samples

        return RefreshPayload(
            token=snapshot.token,
            ann_filters=dict(snapshot.ann_filters),
            allowed_samples=allowed_samples,
            panel_allowed_samples=panel_allowed_samples,
            filtered_ann_pl=filtered_ann_pl,
            meta_for_tree_pl=meta_for_tree_pl,
            contig_records=contig_records,
            filtered_counts=filtered_counts,
            tree_allowed_samples=tree_allowed_samples,
            stats_text=f"Samples: {len(allowed_samples)} | Annotations after filters: {filtered_ann_pl.height}",
        )


    def _collect_predictor_options(self) -> List[str]:
        if "Predictors" not in self._ann_pl.columns:
            return []
        values = (
            self._ann_pl.lazy()
            .select(pl.col("Predictors").cast(pl.Utf8).fill_null("").str.split(",").alias("__pred"))
            .explode("__pred")
            .with_columns(pl.col("__pred").cast(pl.Utf8).str.strip_chars().alias("__pred"))
            .filter(pl.col("__pred").is_not_null() & (pl.col("__pred") != ""))
            .select(pl.col("__pred").unique().sort())
            .collect()
            .get_column("__pred")
            .to_list()
        )
        return [str(v) for v in values]
