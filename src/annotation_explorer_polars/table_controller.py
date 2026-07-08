from __future__ import annotations

import bisect
import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import polars as pl
from PyQt6 import QtCore, QtGui, QtWidgets

from .alignment_dialog import AlignmentDialog
from .annotation_detail_dialog import AnnotationDetailDialog
from .annotation_table import AnnotationTableModel, DetailsDelegate, TagsDelegate
from .constants import ANNOTATION_COLUMNS_TO_SHOW, MISSING_FILTER_TOKEN
from .geo_view import GeoMapView
from .igv_view import GeneSegment, IgvTrackView
from .metadata_filter_dialog import MetadataValueDialog
from .panaroo_views import (
    EntropyTrackTab,
    PanCoreAccessoryTab,
    SimilarityMapTab,
    StructPresenceHeatmapTab,
)
from .sample_card_dialog import SampleCardDialog


class TableControllerMixin:
    def _annotation_table_columns(
        self,
        *,
        include_sample: bool = False,
        include_contig: bool = False,
    ) -> List[Tuple[str, str]]:
        columns: List[Tuple[str, str]] = []
        if include_sample:
            columns.append(("Sample", "Sample"))
        if include_contig:
            columns.append(("Contig", "Contig"))
        columns.extend(ANNOTATION_COLUMNS_TO_SHOW)
        return columns

    def _update_table_search_hint(self, filtered_ann_pl: pl.DataFrame) -> None:
        hint = getattr(self, "table_search_hint", None)
        if hint is None:
            return
        gene_query = self.gene_query.text().strip() if hasattr(self, "gene_query") else ""
        if not gene_query:
            hint.clear()
            hint.setVisible(False)
            return
        hits = int(filtered_ann_pl.height) if filtered_ann_pl is not None else 0
        samples = (
            len(self._collect_sample_list(filtered_ann_pl))
            if filtered_ann_pl is not None and not filtered_ann_pl.is_empty()
            else 0
        )
        hint.setText(
            f"Gene search: {samples} sample{'s' if samples != 1 else ''} found, "
            f"{hits} total hit{'s' if hits != 1 else ''} found"
        )
        hint.setVisible(True)

    def _load_global(self) -> None:
        if not self._should_show_global():
            self._set_table_model([], self._annotation_table_columns(include_sample=True, include_contig=True))
            self.status.setText("Global view | enter a search or apply an annotation filter to show results")
            if self.igv_view is not None:
                self.igv_view.clear()
            return
        df_pl = (
            self._filtered_ann_pl.lazy()
            .sort(["Sample", "Contig", "StartInt", "StopInt", "Gene", "Product"], nulls_last=True)
            .collect()
        )
        columns = self._annotation_table_columns(include_sample=True, include_contig=True)

        if df_pl.is_empty():
            self._set_table_model([], columns)
            self.status.setText("Global view | 0 annotations")
            if self.igv_view is not None:
                self.igv_view.clear()
            return

        row_items: List[Dict[str, Any]] = []
        prev_sample: Optional[str] = None
        for row_data in df_pl.iter_rows(named=True):
            sample = str(row_data.get("Sample") or "").strip()
            if prev_sample is not None and sample != prev_sample:
                row_items.append({"type": "separator"})
            row_items.append({"type": "data", "data": dict(row_data), "is_hit": False})
            prev_sample = sample

        self._set_table_model(row_items, columns)
        sample_count = len(self._collect_sample_list(df_pl))
        self.status.setText(
            f"Global view | {df_pl.height} annotations across {sample_count} samples"
        )

    def _show_annotation_details(self, index: QtCore.QModelIndex) -> None:
        self._open_annotation_details_for_row(index.row())

    def _open_sample_card_for_sample(self, sample: str) -> None:
        sample = str(sample or "").strip()
        if not sample:
            return
        dlg = SampleCardDialog(
            None,
            sample=sample,
            ann_df=self.ann_df,
            meta_df=self.meta_df,
            marker_hits_df=self.marker_hits_df,
            reference_hits_df=self.reference_hits_df,
            output_dir=self.output_dir,
        )
        self._show_floating_dialog(dlg)

    def _open_sample_card_for_row(self, row: int) -> None:
        model = self.table.model()
        if not isinstance(model, AnnotationTableModel):
            return
        row_data = model.row_data(row)
        if not row_data:
            return
        sample = str(row_data.get("Sample") or "").strip()
        self._open_sample_card_for_sample(sample)

    def _open_annotation_details_for_row(self, row: int) -> None:
        model = self.table.model()
        if not isinstance(model, AnnotationTableModel):
            return
        row_data = model.row_data(row)
        if not row_data:
            return
        dlg = AnnotationDetailDialog(
            None,
            row_data,
            row_data.keys(),
            ann_df=self.ann_df,
            meta_df=self.meta_df,
            marker_hits_df=self.marker_hits_df,
            output_dir=self.output_dir,
            pan_reference_sample=self._snippy_reference_sample,
            pan_fallback_tree=self._resolve_raxml_tree_path(),
            sample_hover_callback=self._set_pan_highlight,
        )
        self._show_floating_dialog(dlg)


    def _selected_contig_values(
        self,
        contig: str,
        contigs: Optional[Sequence[str]] = None,
    ) -> List[str]:
        out: List[str] = []
        seen: set[str] = set()
        for raw in (contigs or [contig]):
            text = str(raw or "").strip()
            if not text or text in seen:
                continue
            seen.add(text)
            out.append(text)
        if not out:
            text = str(contig or "").strip()
            if text:
                out.append(text)
        return out

    def _with_synthetic_igv_coords(
        self,
        rows_df: pd.DataFrame,
        contigs: Sequence[str],
    ) -> pd.DataFrame:
        if rows_df.empty or len(contigs) <= 1 or "Contig" not in rows_df.columns:
            return rows_df
        out = rows_df.copy()
        order = {str(contig): idx for idx, contig in enumerate(contigs)}
        out["__ContigOrder"] = out["Contig"].astype(str).map(order).fillna(len(order)).astype(int)
        out = out.sort_values(["__ContigOrder", "StartInt", "StopInt", "Gene", "Product"], kind="stable")

        offsets: Dict[str, int] = {}
        cursor = 0
        gap_bp = 1000
        for contig in contigs:
            contig_rows = out[out["Contig"].astype(str) == str(contig)]
            if contig_rows.empty:
                offsets[str(contig)] = cursor
                continue
            stop_series = pd.to_numeric(contig_rows.get("StopInt"), errors="coerce")
            if stop_series.isna().all():
                stop_series = pd.to_numeric(contig_rows.get("Stop"), errors="coerce")
            contig_len = int(stop_series.max()) if stop_series.notna().any() else 1
            contig_len = max(1, contig_len)
            offsets[str(contig)] = cursor
            cursor += contig_len + gap_bp

        starts = (
            pd.to_numeric(out["StartInt"], errors="coerce")
            if "StartInt" in out.columns
            else pd.Series(np.nan, index=out.index)
        )
        if "Start" in out.columns:
            starts = starts.fillna(pd.to_numeric(out["Start"], errors="coerce"))
        stops = (
            pd.to_numeric(out["StopInt"], errors="coerce")
            if "StopInt" in out.columns
            else pd.Series(np.nan, index=out.index)
        )
        if "Stop" in out.columns:
            stops = stops.fillna(pd.to_numeric(out["Stop"], errors="coerce"))
        offset_series = out["Contig"].astype(str).map(offsets).fillna(0).astype(int)
        out["__IgvStartInt"] = (starts.fillna(0).astype(int) + offset_series).astype(int)
        out["__IgvStopInt"] = (stops.fillna(0).astype(int) + offset_series).astype(int)
        return out


    def _on_igv_gene_double_clicked(self, row_data: Dict[str, Any]) -> None:
        if not row_data:
            return
        dlg = AnnotationDetailDialog(
            None,
            row_data,
            row_data.keys(),
            ann_df=self.ann_df,
            meta_df=self.meta_df,
            marker_hits_df=self.marker_hits_df,
            output_dir=self.output_dir,
            pan_reference_sample=self._snippy_reference_sample,
            pan_fallback_tree=self._resolve_raxml_tree_path(),
            sample_hover_callback=self._set_pan_highlight,
        )
        self._show_floating_dialog(dlg)


    def _set_table_model(self, row_items: List[Dict[str, Any]], columns: Sequence[Tuple[str, str]]) -> None:
        model = AnnotationTableModel(row_items, columns, self.table)
        self._table_model = model
        self.table.setModel(model)
        # Delegates are bound by numeric column index and persist across models.
        # Clear them first to avoid shifted delegates (e.g. details icon in Marker col).
        reset_upto = max(self._delegate_reset_column_count, len(columns))
        for idx in range(reset_upto):
            self.table.setItemDelegateForColumn(idx, None)
        self._delegate_reset_column_count = len(columns)
        self._table_row_bp_index = []
        self._table_bp_values = []
        for row_idx in range(model.rowCount()):
            row_data = model.row_data(row_idx)
            mid = self._row_midpoint_bp(row_data)
            if mid is None:
                continue
            self._table_row_bp_index.append((row_idx, int(mid)))
            self._table_bp_values.append(int(mid))
        if self.table.selectionModel() is not None:
            try:
                self.table.selectionModel().currentRowChanged.disconnect(self._on_table_row_changed)
            except Exception:
                pass
            self.table.selectionModel().currentRowChanged.connect(self._on_table_row_changed)
        self.table.clearSpans()
        for r in model.group_rows:
            self.table.setSpan(r, 0, 1, len(columns))
        for r in model.separator_rows:
            self.table.setSpan(r, 0, 1, len(columns))
            self.table.setRowHeight(r, 2)

        header = self.table.horizontalHeader()
        header.setStretchLastSection(False)
        header.setMinimumSectionSize(30)

        tags_col = next((i for i, (c, _) in enumerate(columns) if c == "Tags"), None)
        if tags_col is not None:
            if self._tags_delegate is None:
                self._tags_delegate = TagsDelegate(self.table)
            self.table.setItemDelegateForColumn(tags_col, self._tags_delegate)
            header.setSectionResizeMode(tags_col, QtWidgets.QHeaderView.ResizeMode.Fixed)
            self.table.setColumnWidth(tags_col, 108)

        details_col = next((i for i, (c, _) in enumerate(columns) if c == "Details"), None)
        if details_col is not None:
            if self._details_delegate is None:
                self._details_delegate = DetailsDelegate(self.table)
                self._details_delegate.clicked.connect(self._open_annotation_details_for_row)
            self.table.setItemDelegateForColumn(details_col, self._details_delegate)
            header.setSectionResizeMode(details_col, QtWidgets.QHeaderView.ResizeMode.Fixed)
            self.table.setColumnWidth(details_col, 30)

        metrics = QtGui.QFontMetrics(self.table.font())
        compact_numeric_width = max(72, metrics.horizontalAdvance("0" * 8) + 20)
        for col_name in ("Start", "Stop"):
            col_idx = next((i for i, (c, _) in enumerate(columns) if c == col_name), None)
            if col_idx is None:
                continue
            header.setSectionResizeMode(col_idx, QtWidgets.QHeaderView.ResizeMode.Fixed)
            self.table.setColumnWidth(col_idx, compact_numeric_width)

        strand_col = next((i for i, (c, _) in enumerate(columns) if c == "Strand"), None)
        if strand_col is not None:
            header.setSectionResizeMode(strand_col, QtWidgets.QHeaderView.ResizeMode.Fixed)
            self.table.setColumnWidth(strand_col, 74)

        for idx in range(len(columns)):
            if idx in {
                details_col,
                tags_col,
                strand_col,
                next((i for i, (c, _) in enumerate(columns) if c == "Start"), None),
                next((i for i, (c, _) in enumerate(columns) if c == "Stop"), None),
            }:
                continue
            header.setSectionResizeMode(idx, QtWidgets.QHeaderView.ResizeMode.Stretch)


    def _load_annotations(
        self,
        sample: str,
        contig: str,
        *,
        contigs: Optional[Sequence[str]] = None,
        display_label: Optional[str] = None,
    ) -> None:
        contig_str = str(contig)
        selected_contigs = self._selected_contig_values(contig_str, contigs)
        is_grouped_plasmid = len(selected_contigs) > 1
        table_columns = self._annotation_table_columns(include_contig=is_grouped_plasmid)
        display_name = str(display_label or f"Contig {contig_str}").strip()
        if self.view_mode == "Context":
            base_contig = self._subset_annotations_pd(
                "base", sample=sample, contig=contig_str, contigs=selected_contigs
            ).sort_values(["StartInt", "StopInt", "Gene", "Product"])
            hits = self._subset_annotations_pd(
                "filtered", sample=sample, contig=contig_str, contigs=selected_contigs
            ).sort_values(["StartInt", "StopInt", "Gene", "Product"])
            base_contig = self._with_synthetic_igv_coords(base_contig, selected_contigs)
            hits = self._with_synthetic_igv_coords(hits, selected_contigs)
            base_contig = base_contig.reset_index(drop=True)
            hit_uids = set(hits.get("AnnUID", pd.Series()).dropna().astype(str).tolist())
            base_contig["IsHit"] = base_contig["AnnUID"].astype(str).isin(hit_uids)
            hit_indices = base_contig.index[base_contig["IsHit"]].tolist()

            intervals = self._build_context_intervals(hit_indices, len(base_contig))
            focus_windows = self._context_windows_from_intervals(base_contig, intervals)
            self._update_igv_view(
                sample,
                contig_str,
                mode="Context",
                contigs=selected_contigs,
                display_label=display_name,
                focus_windows=focus_windows,
            )

            row_items: List[Dict[str, Any]] = []
            shown_rows = 0
            for block_idx, (start_idx, end_idx) in enumerate(intervals):
                if block_idx > 0:
                    row_items.append({"type": "separator"})
                block = base_contig.iloc[start_idx : end_idx + 1]
                shown_rows += len(block)
                for _, rec in block.iterrows():
                    row_data = rec.to_dict()
                    row_data.pop("IsHit", None)
                    row_items.append(
                        {
                            "type": "data",
                            "data": row_data,
                            "is_hit": bool(rec.get("IsHit")),
                        }
                    )
            self._set_table_model(row_items, table_columns)
            self.status.setText(
                f"Context view | Sample {sample} | {display_name} | "
                f"{len(hit_indices)} hits, {shown_rows} rows shown"
            )
            return

        subset = self._subset_annotations_pd(
            "filtered", sample=sample, contig=contig_str, contigs=selected_contigs
        ).sort_values(["StartInt", "StopInt", "Gene", "Product"])
        base_contig = self._subset_annotations_pd(
            "base", sample=sample, contig=contig_str, contigs=selected_contigs
        ).sort_values(["StartInt", "StopInt", "Gene", "Product"])
        base_contig = self._with_synthetic_igv_coords(base_contig, selected_contigs)
        subset = self._with_synthetic_igv_coords(subset, selected_contigs)
        base_contig = base_contig.reset_index(drop=True)

        hit_uids = set(subset.get("AnnUID", pd.Series()).dropna().astype(str).tolist())
        if hit_uids:
            base_contig["IsHit"] = base_contig["AnnUID"].astype(str).isin(hit_uids)
        else:
            base_contig["IsHit"] = False

        highlight_hits = self._has_active_annotation_filter()
        self._update_igv_view(
            sample,
            contig_str,
            mode="Contig",
            contigs=selected_contigs,
            display_label=display_name,
        )
        self._render_table(base_contig, highlight_hits=highlight_hits, columns=table_columns)
        hit_count = int(base_contig["IsHit"].sum()) if "IsHit" in base_contig.columns else 0
        self.status.setText(
            f"Contig view | Sample {sample} | {display_name} | "
            f"{hit_count} hits, {len(base_contig)} annotation rows"
        )


    def _render_table(
        self,
        rows_df: pd.DataFrame,
        highlight_hits: bool = False,
        columns: Optional[Sequence[Tuple[str, str]]] = None,
    ) -> None:
        rows_df = rows_df.reset_index(drop=True)
        rows = rows_df.to_dict("records")
        resolved_columns = list(columns) if columns is not None else self._annotation_table_columns()
        row_items: List[Dict[str, Any]] = []
        for rec in rows:
            row_data = dict(rec)
            row_data.pop("IsHit", None)
            is_hit = bool(rec.get("IsHit")) if highlight_hits else False
            row_items.append({"type": "data", "data": row_data, "is_hit": is_hit})
        self._set_table_model(row_items, resolved_columns)


    def _apply_view_mode(self) -> None:
        if self.view_mode == "Global":
            self._load_global()
            return
        if self.current_sample is None or self.current_contig is None:
            self._set_table_model([], self._annotation_table_columns())
            if self.igv_view is not None:
                self.igv_view.clear()
            return
        self._load_annotations(
            self.current_sample,
            self.current_contig,
            contigs=self.current_contigs,
            display_label=self.current_contig_label,
        )
