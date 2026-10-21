from __future__ import annotations

import contextlib
from typing import Any, Dict, Iterator, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import polars as pl
from PyQt6 import QtCore, QtGui, QtWidgets

from .annotation_table import GAP_ROW_HEIGHT, AnnotationTableModel, DetailsDelegate, TagsDelegate
from .constants import ANNOTATION_COLUMNS_TO_SHOW
from .filters_controller import TAG_OPTIONS
from .model.filter_state import Selection
from .model.run_data import not_installed_key_note
from .model.query_engine import collect_sample_list
from .session.requests import OpenGeneCard
from .shell.commands import INFO_OPEN_ANNOTATION, INFO_OPEN_SAMPLE


ALL_MATCHES_EMPTY_TEXT = "Enter a gene or product search, or apply an annotation filter, to list matching annotations."


def _plural(count: int, noun: str) -> str:
    return f"{int(count)} {noun}{'' if int(count) == 1 else 's'}"


# How spare table width is shared between text columns (Product gets the most).
_TEXT_COLUMN_WEIGHTS = {"Product": 3.0, "Sample": 1.2, "Marker": 1.2, "Gene": 1.0, "Predictors": 1.0, "Contig": 0.6}


class _ResizeWatcher(QtCore.QObject):
    """Calls back whenever the watched widget is resized."""

    def __init__(self, widget: QtWidgets.QWidget, callback) -> None:
        super().__init__(widget)
        self._callback = callback
        widget.installEventFilter(self)

    def eventFilter(self, obj: QtCore.QObject, event: QtCore.QEvent) -> bool:
        if event.type() == QtCore.QEvent.Type.Resize:
            self._callback()
        return False


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
        # The applied search, not the box (which may be ahead of the result).
        if not self._applied_filters().gene:
            hint.clear()
            hint.setVisible(False)
            return
        hits = int(filtered_ann_pl.height) if filtered_ann_pl is not None else 0
        samples = (
            len(collect_sample_list(filtered_ann_pl))
            if filtered_ann_pl is not None and not filtered_ann_pl.is_empty()
            else 0
        )
        hint.setText(
            f"Gene search: {_plural(samples, 'sample')}, {_plural(hits, 'hit')}"
        )
        hint.setVisible(True)

    # The Annotations panels (plan 4.8, D3). Every function below draws into
    # one panel (default: the first one, which is paired with the genome
    # track); only the paired panel touches the track.

    def _annotations_panels(self) -> List[Any]:
        """Every built Annotations panel of the window, the paired one first."""
        primary = getattr(self, "annotations_panel", None)
        panels = [primary] if primary is not None else []
        panel_set = getattr(self, "_panels", None)
        if panel_set is not None:
            for panel in panel_set.panels("annotations"):
                if panel is not primary:
                    panels.append(panel)
        return panels

    def _annotations_target(self, panel: Any = None) -> Any:
        return panel if panel is not None else getattr(self, "annotations_panel", None)

    def _is_paired_panel(self, panel: Any) -> bool:
        return panel is None or panel is getattr(self, "annotations_panel", None)

    def _panel_view_mode(self, panel: Any = None) -> str:
        target = self._annotations_target(panel)
        return target.view_mode() if target is not None else self.view_mode

    def _panel_table(self, panel: Any = None) -> QtWidgets.QTableView:
        target = self._annotations_target(panel)
        table = getattr(target, "table", None) if target is not None else None
        return table if table is not None else self.table

    def _render_annotations(self, panel: Any = None) -> None:
        """Redraw one Annotations panel (and, for the paired one, the track)
        from the applied result, its own view mode and the selection."""
        self._apply_view_mode(panel)

    def _render_annotations_for_selection(self) -> None:
        """The selection was set again (the same contig clicked): redraw the
        panels that show the selected contig, as a fresh click always did."""
        for panel in self._annotations_panels():
            if panel.view_mode() in {"Contig", "Context"}:
                self._render_annotations(panel)

    def _select_contig(
        self,
        sample: Optional[str],
        contig: Optional[str],
        contigs: Optional[Sequence[str]] = None,
        label: Optional[str] = None,
        *,
        origin: str = "",
    ) -> bool:
        """Set the shared selection in one step (one SelectionChanged event);
        returns whether it changed. None clears it."""
        if not sample or contig is None:
            selection = Selection()
        else:
            members = tuple(str(c) for c in (contigs or [contig]) if str(c).strip())
            selection = Selection(
                sample=str(sample),
                contig=str(contig),
                contigs=members or (str(contig),),
                contig_label=str(label or contig),
            )
        if selection == self.ctx.selection:
            return False
        self.ctx.select(selection, origin=origin)
        return True

    def _set_annotations_status(self, text: str, panel: Any = None) -> None:
        """The Annotations panel's own status line (D3: the table's mode text
        never goes to the window's shared status bar)."""
        target = self._annotations_target(panel)
        if target is not None:
            target.set_status_text(text)

    def _clear_paired_track(self, panel: Any) -> None:
        if self._is_paired_panel(panel) and self.igv_view is not None:
            self.igv_view.clear()

    def _load_global(self, panel: Any = None) -> None:
        """All matches: the applied result's filtered annotations across the sample set."""
        if not self._should_show_global():
            self._set_table_model(
                [], self._annotation_table_columns(include_sample=True, include_contig=True), panel=panel
            )
            self._set_annotations_status(ALL_MATCHES_EMPTY_TEXT, panel)
            self._clear_paired_track(panel)
            return
        df_pl = (
            self.ctx.result.filtered_ann.lazy()
            .sort(["Sample", "Contig", "StartInt", "StopInt", "Gene", "Product"], nulls_last=True)
            .collect()
        )
        columns = self._annotation_table_columns(include_sample=True, include_contig=True)

        if df_pl.is_empty():
            self._set_table_model([], columns, panel=panel)
            self._set_annotations_status("No annotations match the filters.", panel)
            self._clear_paired_track(panel)
            return

        # Samples are told apart by an alternating band plus a hairline at each new
        # sample (painted by the delegates), not by separator rows.
        row_items: List[Dict[str, Any]] = []
        prev_sample: Optional[str] = None
        band = -1
        for row_data in df_pl.iter_rows(named=True):
            sample = str(row_data.get("Sample") or "").strip()
            new_sample = sample != prev_sample
            if new_sample:
                band += 1
            row_items.append(
                {
                    "type": "data",
                    "data": dict(row_data),
                    "is_hit": False,
                    "band": band % 2,
                    "group_start": new_sample and prev_sample is not None,
                }
            )
            prev_sample = sample

        self._set_table_model(row_items, columns, panel=panel)
        sample_count = len(collect_sample_list(df_pl))
        self._set_annotations_status(
            f"{df_pl.height} annotation{'s' if df_pl.height != 1 else ''} across "
            f"{sample_count} sample{'s' if sample_count != 1 else ''}",
            panel,
        )

    def _open_sample_card_for_sample(self, sample: str) -> Any:
        """The sample card (an info window owned by this window; one per sample)."""
        sample = str(sample or "").strip()
        if not sample:
            return None
        return self.app_controller.commands.run(INFO_OPEN_SAMPLE, self, sample=sample)

    def _open_annotation_details_for_row(self, row: int, panel: Any = None) -> None:
        model = self._panel_table(panel).model()
        if not isinstance(model, AnnotationTableModel):
            return
        row_data = model.row_data(row)
        if not row_data:
            return
        self._open_gene_card(row_data)

    def _open_gene_card(self, row_data: Dict[str, Any]) -> Any:
        """The gene card for one annotation row (an info window owned by this
        window; one per gene). Its actions run as commands on this window."""
        if not row_data:
            return None
        return self.app_controller.commands.run(INFO_OPEN_ANNOTATION, self, row=row_data)


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
        if row_data:
            self.ctx.request(OpenGeneCard.from_row(row_data), origin="track")


    def _on_annotation_row_changed(
        self,
        panel: Any,
        current: QtCore.QModelIndex,
        previous: QtCore.QModelIndex,
    ) -> None:
        """A row of an Annotations panel became current. The paired panel
        drives the track; another panel in All matches selects the row's
        contig (shown by the panels in Contig or Around matches mode)."""
        if self._is_paired_panel(panel):
            self._on_table_row_changed(current, previous)
            return
        if not current.isValid() or panel.view_mode() != "Global":
            return
        model = panel.table.model() if panel.table is not None else None
        row_data = model.row_data(current.row()) if isinstance(model, AnnotationTableModel) else None
        sample = str((row_data or {}).get("Sample") or "").strip()
        contig = str((row_data or {}).get("Contig") or "").strip()
        if sample and contig:
            self._select_contig(sample, contig, [contig], f"Contig {contig}", origin=panel.uid)


    def _set_table_model(
        self,
        row_items: List[Dict[str, Any]],
        columns: Sequence[Tuple[str, str]],
        *,
        panel: Any = None,
    ) -> None:
        target = self._annotations_target(panel)
        table = self._panel_table(panel)
        paired = self._is_paired_panel(panel)
        model = AnnotationTableModel(row_items, columns, table)
        model.tags_header_tooltip = not_installed_key_note(self.run.tag_databases, TAG_OPTIONS)[1]
        target.model = model
        table.setModel(model)
        target.rendered_key = target.render_key()
        hover = getattr(self.ctx, "hover_state", None)
        if hover is not None and hover.samples and hover.origin != target.uid:
            target.tint_samples(hover.samples)
        if not table.isVisible():
            # A hidden view keeps its old scroll position against the new
            # rows until it is laid out again; start at the top when the
            # Annotations panel is shown (AnnotationsPanel.on_shown).
            table.scrollToTop()
            target.scroll_top_on_show = True
        # Delegates are bound by numeric column index and persist across models.
        # Clear them first to avoid shifted delegates (e.g. details icon in Marker col).
        reset_upto = max(target.delegate_reset_column_count, len(columns))
        for idx in range(reset_upto):
            table.setItemDelegateForColumn(idx, None)
        target.delegate_reset_column_count = len(columns)
        if paired:
            # New rows have no current row yet: the track drops its outline.
            if self.igv_view is not None:
                self.igv_view.set_selected_row(None)
            self._table_row_bp_index = []
            self._table_bp_values = []
            for row_idx in range(model.rowCount()):
                row_data = model.row_data(row_idx)
                mid = self._row_midpoint_bp(row_data)
                if mid is None:
                    continue
                self._table_row_bp_index.append((row_idx, int(mid)))
                self._table_bp_values.append(int(mid))
        if table.selectionModel() is not None:
            # A new model brings a new selection model; connect it once.
            table.selectionModel().currentRowChanged.connect(
                lambda current, previous, _panel=target: self._on_annotation_row_changed(_panel, current, previous)
            )
        table.clearSpans()
        for r in model.group_rows:
            table.setSpan(r, 0, 1, len(columns))
        for r in model.separator_rows:
            table.setSpan(r, 0, 1, len(columns))
            table.setRowHeight(r, GAP_ROW_HEIGHT)

        header = table.horizontalHeader()
        header.setStretchLastSection(False)
        header.setMinimumSectionSize(30)

        tags_col = next((i for i, (c, _) in enumerate(columns) if c == "Tags"), None)
        if tags_col is not None:
            if target.tags_delegate is None:
                target.tags_delegate = TagsDelegate(table)
            table.setItemDelegateForColumn(tags_col, target.tags_delegate)
            header.setSectionResizeMode(tags_col, QtWidgets.QHeaderView.ResizeMode.Fixed)
            table.setColumnWidth(tags_col, 108)

        details_col = next((i for i, (c, _) in enumerate(columns) if c == "Details"), None)
        if details_col is not None:
            if target.details_delegate is None:
                target.details_delegate = DetailsDelegate(table)
                target.details_delegate.clicked.connect(
                    lambda row, _panel=target: self._open_annotation_details_for_row(row, _panel)
                )
            table.setItemDelegateForColumn(details_col, target.details_delegate)
            header.setSectionResizeMode(details_col, QtWidgets.QHeaderView.ResizeMode.Fixed)
            table.setColumnWidth(details_col, 30)

        metrics = QtGui.QFontMetrics(table.font())
        compact_numeric_width = max(72, metrics.horizontalAdvance("0" * 8) + 20)
        for col_name in ("Start", "Stop"):
            col_idx = next((i for i, (c, _) in enumerate(columns) if c == col_name), None)
            if col_idx is None:
                continue
            header.setSectionResizeMode(col_idx, QtWidgets.QHeaderView.ResizeMode.Fixed)
            table.setColumnWidth(col_idx, compact_numeric_width)

        strand_col = next((i for i, (c, _) in enumerate(columns) if c == "Strand"), None)
        if strand_col is not None:
            header.setSectionResizeMode(strand_col, QtWidgets.QHeaderView.ResizeMode.Fixed)
            table.setColumnWidth(strand_col, 74)

        preferred_widths = {
            "Sample": 170,
            "Contig": 170,
            "Predictors": 130,
            "Gene": 150,
            "Product": 280,
            "Marker": 170,
        }
        target.text_columns = {}
        for idx in range(len(columns)):
            if idx in {
                details_col,
                tags_col,
                strand_col,
                next((i for i, (c, _) in enumerate(columns) if c == "Start"), None),
                next((i for i, (c, _) in enumerate(columns) if c == "Stop"), None),
            }:
                continue
            # Stretching every text column made headers and values collapse to
            # a few unreadable characters in a split workspace. Keep useful
            # widths as a floor (the horizontal scrollbar covers narrow panes)
            # and hand any spare width to the text columns in _fit_table_columns.
            header.setSectionResizeMode(idx, QtWidgets.QHeaderView.ResizeMode.Interactive)
            column_name = str(columns[idx][0])
            target.text_columns[idx] = (preferred_widths.get(column_name, 140), column_name)
        if target.fit_filter is None:
            target.fit_filter = _ResizeWatcher(
                table.viewport(), lambda _panel=target: self._fit_table_columns(_panel)
            )
        self._fit_table_columns(target)

    def _fit_table_columns(self, panel: Any = None) -> None:
        """Widen text columns (Product most) so the table fills the viewport."""
        target = self._annotations_target(panel)
        table = self._panel_table(panel)
        text_cols = getattr(target, "text_columns", None) or {}
        model = table.model()
        if model is None or not text_cols:
            return
        other = sum(
            table.columnWidth(i) for i in range(model.columnCount()) if i not in text_cols
        )
        base = sum(width for width, _ in text_cols.values())
        extra = max(0, table.viewport().width() - other - base)
        weights = {idx: _TEXT_COLUMN_WEIGHTS.get(name, 1.0) for idx, (_, name) in text_cols.items()}
        total_weight = sum(weights.values()) or 1.0
        given = 0
        last = max(text_cols)
        for idx, (width, _) in sorted(text_cols.items()):
            share = extra - given if idx == last else int(extra * weights[idx] / total_weight)
            given += share
            table.setColumnWidth(idx, width + share)


    def _load_annotations(
        self,
        sample: str,
        contig: str,
        *,
        contigs: Optional[Sequence[str]] = None,
        display_label: Optional[str] = None,
        panel: Any = None,
    ) -> None:
        contig_str = str(contig)
        selected_contigs = self._selected_contig_values(contig_str, contigs)
        is_grouped_plasmid = len(selected_contigs) > 1
        table_columns = self._annotation_table_columns(include_contig=is_grouped_plasmid)
        display_name = str(display_label or f"Contig {contig_str}").strip()
        paired = self._is_paired_panel(panel)
        if self._panel_view_mode(panel) == "Context":
            base_contig = self.engine.result_annotations(
                self.ctx.result, hits=False, sample=sample, contig=contig_str, contigs=selected_contigs
            ).sort_values(["StartInt", "StopInt", "Gene", "Product"])
            hits = self.engine.result_annotations(
                self.ctx.result, hits=True, sample=sample, contig=contig_str, contigs=selected_contigs
            ).sort_values(["StartInt", "StopInt", "Gene", "Product"])
            base_contig = self._with_synthetic_igv_coords(base_contig, selected_contigs)
            hits = self._with_synthetic_igv_coords(hits, selected_contigs)
            base_contig = base_contig.reset_index(drop=True)
            hit_uids = set(hits.get("AnnUID", pd.Series()).dropna().astype(str).tolist())
            base_contig["IsHit"] = base_contig["AnnUID"].astype(str).isin(hit_uids)
            hit_indices = base_contig.index[base_contig["IsHit"]].tolist()

            intervals = self._build_context_intervals(hit_indices, len(base_contig))
            if paired:
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
            prev_end: Optional[int] = None
            for start_idx, end_idx in intervals:
                if prev_end is not None:
                    hidden = max(0, start_idx - prev_end - 1)
                    row_items.append(
                        {"type": "separator", "label": f"{hidden} gene{'s' if hidden != 1 else ''} not shown"}
                    )
                prev_end = end_idx
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
            self._set_table_model(row_items, table_columns, panel=panel)
            counts = _plural(shown_rows, "row") + " shown"
            if self._has_active_annotation_filter():
                counts = f"{_plural(len(hit_indices), 'hit')}, {counts}"
            self._set_annotations_status(f"Sample {sample} | {display_name} | {counts}", panel)
            return

        subset = self.engine.result_annotations(
            self.ctx.result, hits=True, sample=sample, contig=contig_str, contigs=selected_contigs
        ).sort_values(["StartInt", "StopInt", "Gene", "Product"])
        base_contig = self.engine.result_annotations(
            self.ctx.result, hits=False, sample=sample, contig=contig_str, contigs=selected_contigs
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
        if paired:
            self._update_igv_view(
                sample,
                contig_str,
                mode="Contig",
                contigs=selected_contigs,
                display_label=display_name,
            )
        self._render_table(base_contig, highlight_hits=highlight_hits, columns=table_columns, panel=panel)
        hit_count = int(base_contig["IsHit"].sum()) if "IsHit" in base_contig.columns else 0
        counts = _plural(len(base_contig), "annotation row")
        if highlight_hits:
            counts = f"{_plural(hit_count, 'hit')}, {counts}"
        self._set_annotations_status(f"Sample {sample} | {display_name} | {counts}", panel)


    def _render_table(
        self,
        rows_df: pd.DataFrame,
        highlight_hits: bool = False,
        columns: Optional[Sequence[Tuple[str, str]]] = None,
        panel: Any = None,
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
        self._set_table_model(row_items, resolved_columns, panel=panel)


    @contextlib.contextmanager
    def _batched_annotations_render(self) -> Iterator[None]:
        """Draw each Annotations panel at most once for the steps inside
        (a selection plus a mode change: Open in Context view), at the end."""
        if getattr(self, "_annotations_render_batch", None) is not None:
            yield
            return
        batch: Dict[int, Any] = {}
        self._annotations_render_batch = batch
        try:
            yield
        finally:
            self._annotations_render_batch = None
            for target in batch.values():
                self._apply_view_mode(target)

    def _apply_view_mode(self, panel: Any = None) -> None:
        """Draw one Annotations panel (default: the paired one) for its mode."""
        target = self._annotations_target(panel)
        if target is None:
            return
        batch = getattr(self, "_annotations_render_batch", None)
        if batch is not None:
            batch[id(target)] = target
            return
        if self._panel_view_mode(panel) == "Global":
            self._load_global(panel)
            return
        selection = self.ctx.selection
        if selection.sample is None or selection.contig is None:
            self._set_table_model([], self._annotation_table_columns(), panel=panel)
            self._set_annotations_status("", panel)
            self._clear_paired_track(panel)
            return
        self._load_annotations(
            selection.sample,
            selection.contig,
            contigs=list(selection.contigs) or None,
            display_label=selection.contig_label,
            panel=panel,
        )
