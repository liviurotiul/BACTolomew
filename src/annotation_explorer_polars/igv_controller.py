from __future__ import annotations

import bisect
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from PyQt6 import QtCore, QtGui

from .annotation_table import AnnotationTableModel
from .context_model import gene_type_fill
from .igv_view import BREAK, GAP, LOCUS, GeneSegment
from .model.run_data import clean_text


class IgvControllerMixin:
    def _update_igv_view(
        self,
        sample: str,
        contig: str,
        *,
        mode: str,
        contigs: Optional[Sequence[str]] = None,
        display_label: Optional[str] = None,
        focus_windows: Optional[Sequence[Tuple[int, int]]] = None,
        highlight_bp: Optional[int] = None,
    ) -> None:
        def _clear_igv_state() -> None:
            self._igv_bp_segments = []
            self._igv_raw_contig_len = 0
            self._igv_display_contig_len = 0

        if self.igv_view is None:
            _clear_igv_state()
            return
        if mode not in {"Contig", "Context"}:
            self.igv_view.clear()
            _clear_igv_state()
            return
        if self.ctx.result is None or not self.ctx.result.allowed_samples:
            self.igv_view.clear()
            _clear_igv_state()
            return

        sample_name = str(sample).strip()
        contig_name = str(contig).strip()
        contig_values = [
            str(value).strip()
            for value in (contigs or [contig_name])
            if str(value).strip()
        ]
        if not contig_values and contig_name:
            contig_values = [contig_name]
        is_grouped_plasmid = len(contig_values) > 1
        contig_df = self.engine.result_annotations(self.ctx.result, hits=False, sample=sample_name, contig=contig_name, contigs=contig_values)
        ref_contig_df = pd.DataFrame()
        if (
            self.run.ref_df is not None
            and not self.run.ref_df.empty
            and {"Sample", "Contig"}.issubset(self.run.ref_df.columns)
        ):
            sample_series = self.run.ref_df["Sample"].astype(str).str.strip()
            contig_series = self.run.ref_df["Contig"].astype(str).str.strip()
            if is_grouped_plasmid:
                ref_contig_df = self.run.ref_df[(sample_series == sample_name) & (contig_series.isin(contig_values))]
            else:
                ref_contig_df = self.run.ref_df[(sample_series == sample_name) & (contig_series == contig_name)]

        if contig_df.empty and ref_contig_df.empty:
            self.igv_view.clear()
            _clear_igv_state()
            return

        separator_spans: List[Tuple[int, int, str]] = []

        def _add_synthetic_coords(rows: pd.DataFrame, offsets: Dict[str, int]) -> pd.DataFrame:
            if rows.empty or "Contig" not in rows.columns:
                return rows
            out = rows.copy()
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
            order = {str(name): idx for idx, name in enumerate(contig_values)}
            out["__ContigOrder"] = out["Contig"].astype(str).map(order).fillna(len(order)).astype(int)
            return out.sort_values(["__ContigOrder", "__IgvStartInt", "__IgvStopInt"], kind="stable")

        if is_grouped_plasmid:
            gap_bp = 1000
            offsets: Dict[str, int] = {}
            cursor = 0
            for member in contig_values:
                member_text = str(member)
                max_stop = 0
                if not contig_df.empty and "Contig" in contig_df.columns:
                    member_rows = contig_df[contig_df["Contig"].astype(str) == member_text]
                    if not member_rows.empty:
                        stop_series = pd.to_numeric(member_rows.get("StopInt"), errors="coerce")
                        if stop_series.isna().all():
                            stop_series = pd.to_numeric(member_rows.get("Stop"), errors="coerce")
                        if stop_series.notna().any():
                            max_stop = max(max_stop, int(stop_series.max()))
                if not ref_contig_df.empty and "Contig" in ref_contig_df.columns:
                    member_ref = ref_contig_df[ref_contig_df["Contig"].astype(str) == member_text]
                    if not member_ref.empty and "Stop" in member_ref.columns:
                        ref_stop = pd.to_numeric(member_ref["Stop"], errors="coerce")
                        if ref_stop.notna().any():
                            max_stop = max(max_stop, int(ref_stop.max()))
                member_len = max(1, max_stop)
                offsets[member_text] = cursor
                cursor += member_len
                if member_text != str(contig_values[-1]):
                    separator_spans.append((cursor + 1, cursor + gap_bp, member_text))
                    cursor += gap_bp
            contig_df = _add_synthetic_coords(contig_df, offsets)
            ref_contig_df = _add_synthetic_coords(ref_contig_df, offsets)

        max_stops: List[int] = []
        if not contig_df.empty:
            ann_stop = (
                pd.to_numeric(contig_df["__IgvStopInt"], errors="coerce").max()
                if "__IgvStopInt" in contig_df.columns
                else pd.to_numeric(contig_df["StopInt"], errors="coerce").max()
                if "StopInt" in contig_df.columns
                else pd.to_numeric(contig_df["Stop"], errors="coerce").max()
            )
            if pd.notna(ann_stop):
                max_stops.append(int(ann_stop))
        if not ref_contig_df.empty and "Stop" in ref_contig_df.columns:
            ref_stop = (
                pd.to_numeric(ref_contig_df["__IgvStopInt"], errors="coerce").max()
                if "__IgvStopInt" in ref_contig_df.columns
                else pd.to_numeric(ref_contig_df["Stop"], errors="coerce").max()
            )
            if pd.notna(ref_stop):
                max_stops.append(int(ref_stop))
        if not max_stops:
            self.igv_view.clear()
            _clear_igv_state()
            return
        contig_len = max(max_stops)

        normalized_focus: List[Tuple[int, int]] = []
        if mode == "Context" and focus_windows:
            for start, end in focus_windows:
                try:
                    s = int(start)
                    e = int(end)
                except Exception:
                    continue
                if s > e:
                    s, e = e, s
                s = max(1, s)
                e = min(contig_len, e)
                if s <= e:
                    normalized_focus.append((s, e))
            normalized_focus.sort(key=lambda x: (x[0], x[1]))
            merged_focus: List[Tuple[int, int]] = []
            for s, e in normalized_focus:
                if not merged_focus or s > merged_focus[-1][1] + 1:
                    merged_focus.append((s, e))
                else:
                    prev_s, prev_e = merged_focus[-1]
                    merged_focus[-1] = (prev_s, max(prev_e, e))
            normalized_focus = merged_focus

        # The applied annotation filters decide the hit arrows (D3: never the
        # sample search), so they key the drawn track too.
        applied = self._applied_filters()
        filter_signature = (
            applied.gene.lower(),
            applied.product.lower(),
            applied.molecule,
            tuple(sorted(applied.predictors)),
            tuple(sorted(applied.tags)),
            int(applied.hide_hypothetical),
            int(highlight_bp) if highlight_bp is not None else 0,
        )
        key = (
            sample_name,
            contig_name,
            tuple(contig_values),
            str(display_label or ""),
            mode,
            tuple((int(s), int(e)) for s, e in normalized_focus),
            filter_signature,
        )
        if not self._igv_dirty and self._igv_key == key:
            return

        mapped_contig_len = contig_len
        placeholder_segments: List[Tuple[int, int, int]] = []
        mapping_segments: List[Tuple[int, int, float, float]] = [
            (1, contig_len, 1.0, 1.0)
        ]

        def _identity_map(pos: int) -> float:
            return float(max(1, min(int(pos), contig_len)))

        map_pos = _identity_map
        if mode == "Context" and normalized_focus:
            # Make collapsed context blocks more visible in context mode.
            placeholder_visual_span = 650.0
            pieces: List[Tuple[int, int, float, float, bool]] = []
            src_cursor = 1
            dst_cursor = 1.0
            for s, e in normalized_focus:
                if s > src_cursor:
                    out_len = s - src_cursor
                    if out_len > 0:
                        out_scale = placeholder_visual_span / float(out_len)
                        out_start = dst_cursor
                        out_end = dst_cursor + out_len * out_scale
                        pieces.append((src_cursor, s - 1, out_start, out_scale, False))
                        placeholder_segments.append(
                            (
                                int(round(out_start)),
                                int(round(out_end)),
                                int(out_len),
                            )
                        )
                        dst_cursor = out_end
                in_len = e - s + 1
                pieces.append((s, e, dst_cursor, 1.0, True))
                dst_cursor += in_len
                src_cursor = e + 1
            if src_cursor <= contig_len:
                out_len = contig_len - src_cursor + 1
                if out_len > 0:
                    out_scale = placeholder_visual_span / float(out_len)
                    out_start = dst_cursor
                    out_end = dst_cursor + out_len * out_scale
                    pieces.append((src_cursor, contig_len, out_start, out_scale, False))
                    placeholder_segments.append(
                        (
                            int(round(out_start)),
                            int(round(out_end)),
                            int(out_len),
                        )
                    )
                    dst_cursor = out_end
            mapped_contig_len = max(1, int(round(dst_cursor)))

            def _map_pos(pos: int) -> float:
                p = max(1, min(int(pos), contig_len))
                for src_start, src_end, dst_start, scale, _is_focus in pieces:
                    if p <= src_end:
                        return dst_start + (p - src_start) * scale
                return float(mapped_contig_len)

            map_pos = _map_pos
            mapping_segments = [(a, b, c, d) for a, b, c, d, _f in pieces]

        def _overlaps_focus(start_pos: int, end_pos: int) -> bool:
            if not normalized_focus:
                return True
            s = min(int(start_pos), int(end_pos))
            e = max(int(start_pos), int(end_pos))
            for fs, fe in normalized_focus:
                if s <= fe and e >= fs:
                    return True
            return False

        segments: List[GeneSegment] = []
        for sep_start, sep_end, member_name in separator_spans:
            if mode == "Context" and normalized_focus and not _overlaps_focus(sep_start, sep_end):
                continue
            mapped_start = int(round(map_pos(sep_start)))
            mapped_end = int(round(map_pos(sep_end)))
            if mapped_end <= mapped_start:
                mapped_end = mapped_start + 1
            segments.append(
                GeneSegment(
                    start=mapped_start,
                    end=mapped_end,
                    strand="+",
                    label=f"break after {member_name}",
                    lane=0,
                    kind=BREAK,
                )
            )
        for _, row in contig_df.iterrows():
            start = (
                row.get("__IgvStartInt")
                if pd.notna(row.get("__IgvStartInt"))
                else row.get("StartInt") if pd.notna(row.get("StartInt")) else row.get("Start")
            )
            end = (
                row.get("__IgvStopInt")
                if pd.notna(row.get("__IgvStopInt"))
                else row.get("StopInt") if pd.notna(row.get("StopInt")) else row.get("Stop")
            )
            try:
                start_i = int(start)
                end_i = int(end)
            except Exception:
                continue
            if mode == "Context" and normalized_focus and not _overlaps_focus(start_i, end_i):
                continue
            gene_val = row.get("Gene")
            gene = "" if pd.isna(gene_val) else str(gene_val).strip()
            prod_val = row.get("Product")
            product = "" if pd.isna(prod_val) else str(prod_val).strip()
            label = gene or product or "gene"
            strand = str(row.get("Strand") or "+").strip()
            row_data = row.to_dict()
            # The neighborhood plot's gene type colors (a resistance mutation
            # reads as AMR).
            color = QtGui.QColor(gene_type_fill(row_data))
            lane = 2 if strand == "-" else 1
            mapped_start = int(round(map_pos(start_i)))
            mapped_end = int(round(map_pos(end_i)))
            segments.append(
                GeneSegment(
                    start=mapped_start,
                    end=mapped_end,
                    strand=strand,
                    label=label,
                    color=color,
                    lane=lane,
                    row_data=row_data,
                )
            )

        # Render reference loci as merged blocks (one block per locus label),
        # instead of showing each DIAMOND hit as a separate gene-like segment.
        ref_blocks_by_label: Dict[str, List[Tuple[int, int]]] = {}
        for _, row in ref_contig_df.iterrows():
            start = row.get("__IgvStartInt") if pd.notna(row.get("__IgvStartInt")) else row.get("Start")
            end = row.get("__IgvStopInt") if pd.notna(row.get("__IgvStopInt")) else row.get("Stop")
            try:
                start_i = int(start)
                end_i = int(end)
            except Exception:
                continue
            block_start = min(start_i, end_i)
            block_end = max(start_i, end_i)
            if mode == "Context" and normalized_focus and not _overlaps_focus(block_start, block_end):
                continue

            method_val = row.get("Method")
            method = "" if pd.isna(method_val) else str(method_val).strip()
            product_val = row.get("Product")
            product = "" if pd.isna(product_val) else str(product_val).strip()
            label = method or product or "reference_locus"
            ref_blocks_by_label.setdefault(label, []).append((block_start, block_end))

        for label, spans in ref_blocks_by_label.items():
            if not spans:
                continue
            merged_start = min(s for s, _ in spans)
            merged_end = max(e for _, e in spans)
            mapped_start = int(round(map_pos(merged_start)))
            mapped_end = int(round(map_pos(merged_end)))
            segments.append(
                GeneSegment(
                    start=mapped_start,
                    end=mapped_end,
                    strand="+",
                    label=label,
                    lane=0,
                    kind=LOCUS,
                )
            )

        if mode == "Context" and normalized_focus and placeholder_segments:
            context_edge_gap_px = 14.0
            for mapped_start, mapped_end, bp_count in placeholder_segments:
                if bp_count <= 0:
                    continue
                span = float(mapped_end - mapped_start)
                if span <= 1:
                    continue
                inset = min(context_edge_gap_px, max(0.0, span / 3.0))
                vis_start = int(round(mapped_start + inset))
                vis_end = int(round(mapped_end - inset))
                if vis_end <= vis_start:
                    vis_start = int(round(mapped_start))
                    vis_end = int(round(mapped_end))
                vis_end = max(vis_start + 1, vis_end)

                # Show context signifier on all 3 lanes for visibility.
                segments.append(
                    GeneSegment(
                        start=vis_start,
                        end=vis_end,
                        strand="+",
                        label=f"({bp_count}) bp context",
                            lane=1,
                        kind=GAP,
                    )
                )
                segments.append(
                    GeneSegment(
                        start=vis_start,
                        end=vis_end,
                        strand="+",
                        label="",
                            lane=0,
                        kind=GAP,
                    )
                )
                segments.append(
                    GeneSegment(
                        start=vis_start,
                        end=vis_end,
                        strand="-",
                        label="",
                            lane=2,
                        kind=GAP,
                    )
                )

        hit_bp_positions: List[int] = []
        if highlight_bp is not None:
            try:
                bp = max(1, min(int(highlight_bp), contig_len))
            except Exception:
                bp = None  # type: ignore[assignment]
            if bp is not None:
                if not normalized_focus or any(fs <= bp <= fe for fs, fe in normalized_focus):
                    hit_bp_positions.append(int(round(map_pos(bp))))
        else:
            show_hit_arrows = self._has_active_annotation_filter()
            if show_hit_arrows and self.ctx.result is not None and not self.ctx.result.filtered_ann.is_empty():
                hits_df = self.engine.result_annotations(self.ctx.result, hits=True, sample=sample_name, contig=contig_name, contigs=contig_values)
                if is_grouped_plasmid:
                    hits_df = _add_synthetic_coords(hits_df, offsets)
                for _, row in hits_df.iterrows():
                    start = (
                        row.get("__IgvStartInt")
                        if pd.notna(row.get("__IgvStartInt"))
                        else row.get("StartInt") if pd.notna(row.get("StartInt")) else row.get("Start")
                    )
                    end = (
                        row.get("__IgvStopInt")
                        if pd.notna(row.get("__IgvStopInt"))
                        else row.get("StopInt") if pd.notna(row.get("StopInt")) else row.get("Stop")
                    )
                    try:
                        start_i = int(start)
                        end_i = int(end)
                    except Exception:
                        continue
                    if normalized_focus and not _overlaps_focus(start_i, end_i):
                        continue
                    mid = int((start_i + end_i) / 2)
                    hit_bp_positions.append(int(round(map_pos(mid))))

        if hit_bp_positions:
            hit_bp_positions = sorted(
                set(max(1, min(int(bp), int(mapped_contig_len))) for bp in hit_bp_positions)
            )

        self._igv_bp_segments = mapping_segments
        self._igv_raw_contig_len = int(contig_len)
        self._igv_display_contig_len = int(mapped_contig_len)
        self.igv_view.set_track(
            segments,
            mapped_contig_len,
            bp_per_px=10.0,
            hit_bp_positions=hit_bp_positions,
        )
        self._igv_key = key
        self._igv_dirty = False
        if highlight_bp is None:
            # A new contig or mode opens on the selected gene, else on the
            # first searched gene (All matches centers on its row itself).
            target = self.igv_view.selected_gene_bp()
            if target is None and hit_bp_positions:
                target = hit_bp_positions[0]
            if target is not None:
                self.igv_view.center_on_bp(target)


    def _on_table_row_changed(
        self,
        current: QtCore.QModelIndex,
        _previous: QtCore.QModelIndex,
    ) -> None:
        model = self.table.model()
        row_data = (
            model.row_data(current.row())
            if current.isValid() and isinstance(model, AnnotationTableModel)
            else None
        )
        if self.igv_view is not None:
            # The table's current row is the track's selected gene, also when
            # scrolling the track moved it.
            self.igv_view.set_selected_row(row_data)
        if self._syncing_igv_to_table or not row_data:
            return
        if self.view_mode == "Global":
            self._show_global_context_in_igv(row_data)
            return
        if self.view_mode not in {"Contig", "Context"}:
            return
        if self.igv_view is None:
            return
        mid = self._row_midpoint_bp(row_data)
        if mid is None:
            return
        mapped_bp = self._map_raw_bp_to_igv_bp(mid)
        if mapped_bp is None:
            return
        # The table scrolled to show the row: keep that scroll from
        # recentering the track on the row in the middle of the table.
        self._table_scroll_timer.stop()
        self._animate_igv_center_to_bp(mapped_bp)


    def _animate_igv_center_to_bp(self, bp: int) -> None:
        if self.igv_view is None:
            return
        self._igv_sync_token += 1
        token = self._igv_sync_token
        self._syncing_table_to_igv = True
        started = self.igv_view.animate_center_on_bp(bp, duration_ms=self._sync_animation_ms)
        if not started:
            if token == self._igv_sync_token:
                self._syncing_table_to_igv = False
            return

        def _release() -> None:
            if token == self._igv_sync_token:
                self._syncing_table_to_igv = False

        QtCore.QTimer.singleShot(self._sync_animation_ms + 40, _release)


    def _animate_table_to_row_center(
        self,
        row: int,
        index: QtCore.QModelIndex,
    ) -> None:
        vbar = self.table.verticalScrollBar()
        viewport = self.table.viewport()
        if vbar is None or viewport is None:
            return
        row_height = max(1, self.table.rowHeight(row))
        row_pos = self.table.rowViewportPosition(row)
        desired_pos = int((viewport.height() - row_height) / 2)
        target = int(vbar.value() + row_pos - desired_pos)
        target = max(vbar.minimum(), min(vbar.maximum(), target))

        self._table_sync_token += 1
        token = self._table_sync_token
        self._syncing_igv_to_table = True
        self.table.setCurrentIndex(index)

        if abs(vbar.value() - target) <= 1:
            vbar.setValue(target)
            if token == self._table_sync_token:
                self._syncing_igv_to_table = False
            return

        if self._table_scroll_animation is not None:
            self._table_scroll_animation.stop()
            self._table_scroll_animation.deleteLater()
            self._table_scroll_animation = None

        anim = QtCore.QPropertyAnimation(vbar, b"value", self)
        anim.setDuration(self._sync_animation_ms)
        anim.setStartValue(int(vbar.value()))
        anim.setEndValue(int(target))
        anim.setEasingCurve(QtCore.QEasingCurve.Type.InOutCubic)

        def _finish() -> None:
            if self._table_scroll_animation is anim:
                self._table_scroll_animation = None
            if token == self._table_sync_token:
                self._syncing_igv_to_table = False

        anim.finished.connect(_finish)
        self._table_scroll_animation = anim
        anim.start()


    def _row_midpoint_bp(self, row_data: Optional[Dict[str, Any]]) -> Optional[int]:
        if not row_data:
            return None
        start = (
            row_data.get("__IgvStartInt")
            if pd.notna(row_data.get("__IgvStartInt"))
            else row_data.get("StartInt") if pd.notna(row_data.get("StartInt")) else row_data.get("Start")
        )
        end = (
            row_data.get("__IgvStopInt")
            if pd.notna(row_data.get("__IgvStopInt"))
            else row_data.get("StopInt") if pd.notna(row_data.get("StopInt")) else row_data.get("Stop")
        )
        try:
            start_i = int(start)
            end_i = int(end)
        except Exception:
            return None
        return int((start_i + end_i) / 2)


    def _map_raw_bp_to_igv_bp(self, raw_bp: int) -> Optional[int]:
        if self._igv_raw_contig_len <= 0:
            return None
        p = max(1, min(int(raw_bp), int(self._igv_raw_contig_len)))
        segments = self._igv_bp_segments
        if not segments:
            return p
        for src_start, src_end, dst_start, scale in segments:
            if p <= src_end:
                mapped = dst_start + (p - src_start) * scale
                return int(round(mapped))
        last_start, last_end, last_dst_start, last_scale = segments[-1]
        mapped = last_dst_start + (last_end - last_start) * last_scale
        return int(round(mapped))


    def _map_igv_bp_to_raw_bp(self, igv_bp: int) -> Optional[int]:
        if self._igv_raw_contig_len <= 0:
            return None
        p = max(1.0, min(float(igv_bp), float(max(1, self._igv_display_contig_len))))
        segments = self._igv_bp_segments
        if not segments:
            return int(round(p))
        for src_start, src_end, dst_start, scale in segments:
            dst_end = dst_start + (src_end - src_start) * scale
            if p <= dst_end:
                if scale == 0:
                    return src_start
                raw = src_start + (p - dst_start) / scale
                raw_i = int(round(raw))
                return max(1, min(raw_i, int(self._igv_raw_contig_len)))
        return int(self._igv_raw_contig_len)


    def _on_table_scroll_range_changed(self, _minimum: int, _maximum: int) -> None:
        # New rows or a resized table (a splitter, the window): Qt clamps the
        # scroll to the new range right after. That is not the user
        # scrolling, so the track must not follow it.
        self._table_scroll_range_changed = True
        QtCore.QTimer.singleShot(0, lambda: setattr(self, "_table_scroll_range_changed", False))


    def _on_table_scroll_changed(self, _value: int) -> None:
        if self._table_scroll_range_changed:
            self._table_scroll_range_changed = False
            return
        if self._syncing_igv_to_table or self.view_mode not in {"Contig", "Context"}:
            return
        self._table_scroll_timer.start()


    def _nearest_annotation_row(
        self,
        model: AnnotationTableModel,
        row: int,
    ) -> Optional[int]:
        row_count = model.rowCount()
        if row_count <= 0:
            return None
        base = max(0, min(int(row), row_count - 1))
        if model.row_data(base) is not None:
            return base
        for offset in range(1, row_count):
            up = base - offset
            down = base + offset
            if up >= 0 and model.row_data(up) is not None:
                return up
            if down < row_count and model.row_data(down) is not None:
                return down
        return None


    def _sync_igv_from_table_scroll(self) -> None:
        if self._syncing_igv_to_table or self.igv_view is None:
            return
        model = self.table.model()
        if not isinstance(model, AnnotationTableModel):
            return
        viewport = self.table.viewport()
        if viewport is None:
            return
        mid_y = max(0, viewport.height() // 2)
        row = self.table.rowAt(mid_y)
        if row < 0:
            row = self.table.rowAt(0)
        if row < 0:
            return
        row = self._nearest_annotation_row(model, row)
        if row is None:
            return
        row_data = model.row_data(row)
        mid = self._row_midpoint_bp(row_data)
        if mid is None:
            return
        mapped = self._map_raw_bp_to_igv_bp(mid)
        if mapped is None:
            return
        self._animate_igv_center_to_bp(mapped)


    def _on_igv_center_bp_changed(self, bp: int) -> None:
        if self._syncing_table_to_igv or self.view_mode not in {"Contig", "Context"}:
            return
        self._pending_igv_center_bp = int(bp)
        self._igv_scroll_timer.start()


    def _find_nearest_table_row_by_bp(self, raw_bp: int) -> Optional[int]:
        if not self._table_row_bp_index or not self._table_bp_values:
            return None
        bp = int(raw_bp)
        idx = bisect.bisect_left(self._table_bp_values, bp)
        if idx <= 0:
            return self._table_row_bp_index[0][0]
        if idx >= len(self._table_bp_values):
            return self._table_row_bp_index[-1][0]
        prev_bp = self._table_bp_values[idx - 1]
        next_bp = self._table_bp_values[idx]
        if abs(bp - prev_bp) <= abs(next_bp - bp):
            return self._table_row_bp_index[idx - 1][0]
        return self._table_row_bp_index[idx][0]


    def _sync_table_from_igv_center(self) -> None:
        if self._syncing_table_to_igv:
            return
        if self._pending_igv_center_bp is None:
            return
        model = self.table.model()
        if not isinstance(model, AnnotationTableModel):
            return
        raw_bp = self._map_igv_bp_to_raw_bp(self._pending_igv_center_bp)
        if raw_bp is None:
            return
        target_row = self._find_nearest_table_row_by_bp(raw_bp)
        if target_row is None:
            return
        index = model.index(target_row, 0)
        if not index.isValid():
            return
        self._animate_table_to_row_center(target_row, index)


    def _build_context_intervals(
        self,
        hit_indices: Sequence[int],
        total_rows: int,
    ) -> List[Tuple[int, int]]:
        if total_rows <= 0 or not hit_indices:
            return []
        window = max(0, int(self.context_window))
        sorted_hits = sorted({int(i) for i in hit_indices if 0 <= int(i) < total_rows})
        if not sorted_hits:
            return []

        intervals: List[Tuple[int, int]] = []
        for idx in sorted_hits:
            start = max(0, idx - window)
            end = min(total_rows - 1, idx + window)
            if not intervals or start > intervals[-1][1] + 1:
                intervals.append((start, end))
            else:
                prev_start, prev_end = intervals[-1]
                intervals[-1] = (prev_start, max(prev_end, end))
        return intervals


    def _context_windows_from_intervals(
        self,
        contig_df: pd.DataFrame,
        intervals: Sequence[Tuple[int, int]],
    ) -> List[Tuple[int, int]]:
        windows: List[Tuple[int, int]] = []
        if contig_df.empty or not intervals:
            return windows

        for start_idx, end_idx in intervals:
            block = contig_df.iloc[start_idx : end_idx + 1]
            if block.empty:
                continue
            start_series = (
                pd.to_numeric(block["__IgvStartInt"], errors="coerce")
                if "__IgvStartInt" in block.columns
                else pd.Series(np.nan, index=block.index)
            )
            stop_series = (
                pd.to_numeric(block["__IgvStopInt"], errors="coerce")
                if "__IgvStopInt" in block.columns
                else pd.Series(np.nan, index=block.index)
            )
            if start_series.isna().all():
                start_series = (
                    pd.to_numeric(block["StartInt"], errors="coerce")
                    if "StartInt" in block.columns
                    else pd.Series(np.nan, index=block.index)
                )
            if stop_series.isna().all():
                stop_series = (
                    pd.to_numeric(block["StopInt"], errors="coerce")
                    if "StopInt" in block.columns
                    else pd.Series(np.nan, index=block.index)
                )
            if start_series.isna().all():
                start_series = (
                    pd.to_numeric(block["Start"], errors="coerce")
                    if "Start" in block.columns
                    else pd.Series(np.nan, index=block.index)
                )
            if stop_series.isna().all():
                stop_series = (
                    pd.to_numeric(block["Stop"], errors="coerce")
                    if "Stop" in block.columns
                    else pd.Series(np.nan, index=block.index)
                )
            if start_series.isna().all() or stop_series.isna().all():
                continue
            win_start = int(start_series.min())
            win_end = int(stop_series.max())
            if win_start > win_end:
                win_start, win_end = win_end, win_start
            windows.append((win_start, win_end))

        if not windows:
            return []
        windows.sort(key=lambda x: (x[0], x[1]))
        merged: List[Tuple[int, int]] = [windows[0]]
        for start, end in windows[1:]:
            prev_start, prev_end = merged[-1]
            if start <= prev_end + 1:
                merged[-1] = (prev_start, max(prev_end, end))
            else:
                merged.append((start, end))
        return merged


    def _find_context_window_for_row(
        self,
        row_data: Dict[str, Any],
        flank_genes: int = 15,
    ) -> Optional[Tuple[str, str, Tuple[int, int]]]:
        if self.ctx.result is None or not self.ctx.result.allowed_samples:
            return None
        sample = clean_text(row_data.get("Sample"))
        contig = clean_text(row_data.get("Contig"))
        if not sample or not contig:
            return None

        contig_df = self.engine.result_annotations(self.ctx.result, hits=False, sample=sample, contig=contig).sort_values(
            ["StartInt", "StopInt", "Gene", "Product"]
        )
        if contig_df.empty:
            return None
        contig_df = contig_df.reset_index(drop=True)

        target_idx: Optional[int] = None
        ann_uid = clean_text(row_data.get("AnnUID"))
        if ann_uid and "AnnUID" in contig_df.columns:
            matches = contig_df.index[contig_df["AnnUID"].astype(str) == ann_uid].tolist()
            if matches:
                target_idx = int(matches[0])

        if target_idx is None:
            ann_id = clean_text(row_data.get("AnnotationID"))
            if ann_id and "AnnotationID" in contig_df.columns:
                matches = contig_df.index[contig_df["AnnotationID"].astype(str) == ann_id].tolist()
                if matches:
                    target_idx = int(matches[0])

        if target_idx is None:
            target_bp = self._row_midpoint_bp(row_data)
            if target_bp is not None:
                mids = (
                    pd.to_numeric(contig_df.get("StartInt"), errors="coerce").fillna(
                        pd.to_numeric(contig_df.get("Start"), errors="coerce")
                    )
                    + pd.to_numeric(contig_df.get("StopInt"), errors="coerce").fillna(
                        pd.to_numeric(contig_df.get("Stop"), errors="coerce")
                    )
                ) / 2.0
                if mids.notna().any():
                    target_idx = int((mids - float(target_bp)).abs().idxmin())

        if target_idx is None:
            return None

        flank = max(0, int(flank_genes))
        start_idx = max(0, target_idx - flank)
        end_idx = min(len(contig_df) - 1, target_idx + flank)
        windows = self._context_windows_from_intervals(contig_df, [(start_idx, end_idx)])
        if not windows:
            return None
        return sample, contig, windows[0]


    def _show_global_context_in_igv(self, row_data: Dict[str, Any]) -> None:
        if self.igv_view is None:
            return
        focus = self._find_context_window_for_row(row_data, flank_genes=15)
        if focus is None:
            self.igv_view.clear()
            self._igv_bp_segments = []
            self._igv_raw_contig_len = 0
            self._igv_display_contig_len = 0
            return
        sample, contig, window = focus
        target_bp = self._row_midpoint_bp(row_data)
        # The paired table's row selects its contig (one SelectionChanged; the
        # Annotations panels in Contig or Around matches mode load it).
        self._select_contig(sample, contig, [contig], f"Contig {contig}", origin="annotations")
        self._igv_dirty = True
        self._update_igv_view(
            sample,
            contig,
            mode="Context",
            focus_windows=[window],
            highlight_bp=target_bp,
        )
        mapped = self._map_raw_bp_to_igv_bp(target_bp) if target_bp is not None else None
        if mapped is not None and self.igv_view is not None:
            # Bring the row's gene, under the search pin, into view.
            self.igv_view.center_on_bp(mapped)
