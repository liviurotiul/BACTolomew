from __future__ import annotations

from collections import Counter
import html
import math
import re
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets

from .panaroo_common import _METADATA_STRIP_PALETTE, _file_signature, _natural_sort_key
from .panaroo_data import (
    _find_fastgear_sample_clades,
    _find_panaroo_out,
    _load_fastgear_clusters,
    _load_fastgear_lineages,
    _load_presence_matrix,
    _select_samples_from_presence,
)
from .panaroo_widgets import _MetaValueChip, _PgUnavailableWidget

try:
    from scipy.spatial.distance import pdist, squareform
except Exception:
    pdist = None
    squareform = None

try:
    import pyqtgraph as pg
except Exception:
    pg = None


class SimilarityMapTab(QtWidgets.QWidget):
    _PREVALENCE_MIN = 0.05
    _PREVALENCE_MAX = 0.95
    _MIN_GENES_AFTER_PREVALENCE_FILTER = 64
    _CLUSTER_PALETTE = [
        "#1d4ed8",
        "#dc2626",
        "#059669",
        "#7c3aed",
        "#ea580c",
        "#0f766e",
        "#be185d",
        "#4f46e5",
        "#65a30d",
        "#b45309",
        "#0891b2",
        "#4338ca",
        "#b91c1c",
        "#15803d",
        "#a21caf",
        "#0369a1",
    ]

    _MISSING_VALUE_LABEL = "Missing"

    def __init__(
        self,
        output_dir: Path,
        metadata_df: Optional[pd.DataFrame] = None,
        show_metadata_controls: bool = True,
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self._output_dir = Path(output_dir)
        self._sample_filter: Optional[set[str]] = None
        self._highlight_sample: Optional[str] = None
        self._heatmap_row_highlight: Optional[QtWidgets.QGraphicsRectItem] = None
        self._heatmap_col_highlight: Optional[QtWidgets.QGraphicsRectItem] = None
        self._heatmap_rows_multi_highlight: Optional[QtWidgets.QGraphicsPathItem] = None
        self._heatmap_cols_multi_highlight: Optional[QtWidgets.QGraphicsPathItem] = None
        self._matrix_sample_count = 0
        self._sample_index_by_name: Dict[str, int] = {}
        self._matrix_display_index_by_sample: Dict[str, int] = {}
        self._metadata_by_sample: Dict[str, Dict[str, object]] = {}
        self._metadata_columns: List[str] = []
        self._visible_similarity_samples: List[str] = []
        self._meta_hover_value: Optional[str] = None
        self._meta_hover_indices: Set[int] = set()
        self._meta_controls: Optional[QtWidgets.QWidget] = None
        self._meta_column_combo: Optional[QtWidgets.QComboBox] = None
        self._meta_chip_container: Optional[QtWidgets.QWidget] = None
        self._meta_chip_layout: Optional[QtWidgets.QGridLayout] = None

        self._set_metadata(metadata_df)
        self._show_metadata_controls = bool(show_metadata_controls)
        layout = QtWidgets.QVBoxLayout(self)
        self._info = QtWidgets.QLabel()
        self._info.setWordWrap(True)
        layout.addWidget(self._info)
        self._build_meta_controls(layout)
        if pg is None:
            layout.addWidget(_PgUnavailableWidget("PyQtGraph is required for the similarity map view."))
            return

        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal, self)
        split.setChildrenCollapsible(False)
        layout.addWidget(split, stretch=1)

        self._heatmap_plot = pg.PlotWidget()
        self._heatmap_plot.setBackground("#ffffff")
        self._heatmap_plot.setLabel("left", "Samples")
        self._heatmap_plot.setLabel("bottom", "Samples")
        self._heatmap_plot.getPlotItem().setMenuEnabled(False)
        self._heatmap_plot.invertY(True)
        self._heatmap_plot.getViewBox().setAspectLocked(True, ratio=1.0)
        self._heatmap_img = pg.ImageItem()
        self._heatmap_plot.addItem(self._heatmap_img)
        split.addWidget(self._heatmap_plot)

        self._scatter_plot = pg.PlotWidget()
        self._scatter_plot.setBackground("#ffffff")
        self._scatter_plot.setLabel("left", "MDS2")
        self._scatter_plot.setLabel("bottom", "MDS1")
        self._scatter_plot.getPlotItem().setMenuEnabled(False)
        split.addWidget(self._scatter_plot)
        self._scatter_hover_proxy = pg.SignalProxy(
            self._scatter_plot.scene().sigMouseMoved,
            rateLimit=60,
            slot=self._on_scatter_mouse_moved,
        )
        split.setSizes([740, 520])
        self.refresh()

    def _set_metadata(self, metadata_df: Optional[pd.DataFrame]) -> None:
        self._metadata_by_sample = {}
        self._metadata_columns = []
        if metadata_df is None or metadata_df.empty or "Sample" not in metadata_df.columns:
            return
        cols = [str(c) for c in metadata_df.columns if str(c) != "Sample"]
        if not cols:
            return
        self._metadata_columns = cols
        for _, row in metadata_df.iterrows():
            sample = str(row.get("Sample", "")).strip()
            if not sample or sample in self._metadata_by_sample:
                continue
            self._metadata_by_sample[sample] = {col: row.get(col, np.nan) for col in cols}

    def _build_meta_controls(self, parent_layout: QtWidgets.QVBoxLayout) -> None:
        box = QtWidgets.QWidget(self)
        v = QtWidgets.QVBoxLayout(box)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(4)

        row = QtWidgets.QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(6)
        lbl = QtWidgets.QLabel("Metadata")
        lbl.setStyleSheet("font-weight: 600; color: #2b3c4d;")
        row.addWidget(lbl)
        combo = QtWidgets.QComboBox()
        combo.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToContents)
        combo.currentIndexChanged.connect(self._on_meta_column_changed)
        row.addWidget(combo, stretch=1)
        v.addLayout(row)

        chips_wrap = QtWidgets.QWidget()
        chips_grid = QtWidgets.QGridLayout(chips_wrap)
        chips_grid.setContentsMargins(0, 0, 0, 0)
        chips_grid.setHorizontalSpacing(6)
        chips_grid.setVerticalSpacing(6)
        v.addWidget(chips_wrap)

        self._meta_controls = box
        self._meta_column_combo = combo
        self._meta_chip_container = chips_wrap
        self._meta_chip_layout = chips_grid
        parent_layout.addWidget(box)

        self._refresh_meta_columns()

    def _refresh_meta_columns(self) -> None:
        if self._meta_controls is None or self._meta_column_combo is None:
            return
        combo = self._meta_column_combo
        combo.blockSignals(True)
        current = combo.currentText().strip()
        combo.clear()
        combo.addItems(self._metadata_columns)
        if current and current in self._metadata_columns:
            combo.setCurrentText(current)
        combo.blockSignals(False)
        has_metadata = bool(self._metadata_columns) and self._show_metadata_controls
        self._meta_controls.setVisible(has_metadata)
        if has_metadata:
            self._rebuild_meta_value_chips()

    @classmethod
    def _norm_meta_value(cls, value: object) -> str:
        if pd.isna(value):
            return cls._MISSING_VALUE_LABEL
        txt = str(value).strip()
        if not txt:
            return cls._MISSING_VALUE_LABEL
        low = txt.lower()
        if low in {"none", "nan", "na", "n/a", "null", "nil", "missing", "unknown", "not available"}:
            return cls._MISSING_VALUE_LABEL
        return txt

    def _meta_values_for_visible_samples(self, column: str) -> Dict[str, List[str]]:
        global_values = {self._MISSING_VALUE_LABEL}
        for _sample, row in self._metadata_by_sample.items():
            raw = row.get(column, np.nan)
            global_values.add(self._norm_meta_value(raw))

        by_value: Dict[str, List[str]] = {v: [] for v in global_values}
        for sample in self._visible_similarity_samples:
            row = self._metadata_by_sample.get(sample)
            raw = row.get(column, np.nan) if row is not None else np.nan
            val = self._norm_meta_value(raw)
            by_value.setdefault(val, []).append(sample)
        return by_value

    def _clear_meta_chip_layout(self) -> None:
        if self._meta_chip_layout is None:
            return
        while self._meta_chip_layout.count():
            item = self._meta_chip_layout.takeAt(0)
            w = item.widget()
            if w is not None:
                w.deleteLater()

    def _rebuild_meta_value_chips(self) -> None:
        if self._meta_column_combo is None or self._meta_chip_layout is None:
            return
        column = self._meta_column_combo.currentText().strip()
        self._clear_meta_chip_layout()
        if not column:
            return
        grouped = self._meta_values_for_visible_samples(column)
        values = sorted(grouped.keys(), key=lambda v: (v != self._MISSING_VALUE_LABEL, v.lower()))
        max_cols = 5
        for i, value in enumerate(values):
            samples = grouped.get(value, [])
            chip = _MetaValueChip(value, len(samples))
            chip.hovered.connect(self._on_meta_value_hovered)
            chip.unhovered.connect(self._on_meta_value_unhovered)
            r = i // max_cols
            c = i % max_cols
            self._meta_chip_layout.addWidget(chip, r, c)
        self._meta_chip_layout.setColumnStretch(max_cols, 1)
        self._meta_chip_layout.setRowStretch((len(values) // max_cols) + 1, 1)

    def _on_meta_column_changed(self, _idx: int) -> None:
        self._meta_hover_value = None
        self._meta_hover_indices = set()
        self._rebuild_meta_value_chips()
        self._set_scatter_hover(None)

    def _on_meta_value_hovered(self, value: str) -> None:
        if self._meta_column_combo is None:
            return
        col = self._meta_column_combo.currentText().strip()
        if not col:
            return
        grouped = self._meta_values_for_visible_samples(col)
        samples = grouped.get(value, [])
        indices = {self._sample_index_by_name[s] for s in samples if s in self._sample_index_by_name}
        self._meta_hover_value = value
        self._meta_hover_indices = indices
        self._set_scatter_hover(None)
        self._apply_metadata_scatter_override()

    def _on_meta_value_unhovered(self, value: str) -> None:
        if self._meta_hover_value != value:
            return
        self._meta_hover_value = None
        self._meta_hover_indices = set()
        self._set_scatter_hover(None)

    def set_sample_filter(self, samples: Optional[Sequence[str]]) -> None:
        next_filter = {str(s).strip() for s in samples if str(s).strip()} if samples is not None else None
        if next_filter == self._sample_filter:
            return
        self._sample_filter = next_filter
        self.refresh()

    def highlight_sample(self, sample: Optional[str]) -> None:
        value = str(sample).strip() if sample is not None else ""
        self._highlight_sample = value or None
        if not self._show_metadata_controls:
            self._meta_hover_value = None
            self._meta_hover_indices = set()
        scatter_idx: Optional[int] = None
        if self._highlight_sample:
            scatter_idx = self._sample_index_by_name.get(self._highlight_sample)
        self._set_scatter_hover(scatter_idx)

    def highlight_samples(self, samples: Optional[Sequence[str]]) -> None:
        if not samples:
            self._meta_hover_value = None
            self._meta_hover_indices = set()
            self._set_scatter_hover(None)
            self._apply_metadata_scatter_override()
            return
        wanted = {str(s).strip() for s in samples if str(s).strip()}
        idxs = {idx for name, idx in self._sample_index_by_name.items() if name in wanted}
        self._meta_hover_value = "__external__"
        self._meta_hover_indices = idxs
        self._set_scatter_hover(None)
        self._apply_metadata_scatter_override()

    def _set_ticks(self, plot: "pg.PlotWidget", sample_labels: Sequence[str]) -> None:
        n = len(sample_labels)
        if n <= 0:
            plot.getAxis("bottom").setTicks([[]])
            plot.getAxis("left").setTicks([[]])
            return
        tick_step = max(1, n // 24)
        ticks_x = [(i, sample_labels[i]) for i in range(0, n, tick_step)]
        ticks_y = [(i, sample_labels[::-1][i]) for i in range(0, n, tick_step)]
        plot.getAxis("bottom").setTicks([ticks_x])
        plot.getAxis("left").setTicks([ticks_y])

    def _apply_square_heatmap_view(self, n_samples: int) -> None:
        vb = self._heatmap_plot.getViewBox()
        if vb is None:
            return
        size = float(max(1, n_samples))
        vb.setAspectLocked(True, ratio=1.0)
        vb.setRange(
            xRange=(0.0, size),
            yRange=(0.0, size),
            padding=0.0,
            disableAutoRange=True,
        )

    @staticmethod
    def _safe_levels(arr: np.ndarray, fallback: Tuple[float, float]) -> Tuple[float, float]:
        lo = float(np.nanmin(arr))
        hi = float(np.nanmax(arr))
        if not np.isfinite(lo):
            lo = fallback[0]
        if not np.isfinite(hi):
            hi = fallback[1]
        if hi <= lo:
            hi = lo + 1e-6
        return lo, hi

    @staticmethod
    def _classical_mds(dist: np.ndarray) -> np.ndarray:
        n = int(dist.shape[0])
        if n <= 1:
            return np.zeros((n, 2), dtype=float)
        d2 = np.square(np.asarray(dist, dtype=float))
        h = np.eye(n) - (np.ones((n, n), dtype=float) / float(n))
        b = -0.5 * h @ d2 @ h
        try:
            eigvals, eigvecs = np.linalg.eigh(b)
        except np.linalg.LinAlgError:
            return np.zeros((n, 2), dtype=float)
        order = np.argsort(eigvals)[::-1]
        eigvals = eigvals[order]
        eigvecs = eigvecs[:, order]
        pos = np.maximum(eigvals[:2], 0.0)
        coords = eigvecs[:, :2] * np.sqrt(pos[np.newaxis, :])
        if coords.shape[1] < 2:
            coords = np.hstack([coords, np.zeros((n, 2 - coords.shape[1]), dtype=float)])
        return coords[:, :2]

    @staticmethod
    def _kmedoids(dist: np.ndarray, k: int, max_iter: int = 50) -> np.ndarray:
        n = int(dist.shape[0])
        if n == 0:
            return np.zeros(0, dtype=int)
        if k <= 1 or n == 1:
            return np.zeros(n, dtype=int)
        k = max(2, min(int(k), n))

        mean_d = np.nanmean(dist, axis=1)
        first = int(np.nanargmin(mean_d))
        medoids = [first]
        while len(medoids) < k:
            to_m = dist[:, medoids]
            if to_m.ndim == 1:
                to_m = to_m[:, np.newaxis]
            nearest = np.nanmin(to_m, axis=1)
            for m in medoids:
                nearest[m] = -np.inf
            nxt = int(np.nanargmax(nearest))
            if nxt in medoids:
                remaining = [i for i in range(n) if i not in medoids]
                if not remaining:
                    break
                nxt = int(remaining[0])
            medoids.append(nxt)
        medoids = np.array(medoids[:k], dtype=int)

        labels = np.zeros(n, dtype=int)
        for _ in range(max_iter):
            d_to_m = dist[:, medoids]
            labels = np.argmin(d_to_m, axis=1).astype(int)
            new_medoids = medoids.copy()
            for cid in range(len(medoids)):
                members = np.flatnonzero(labels == cid)
                if members.size == 0:
                    dcur = np.min(d_to_m, axis=1)
                    for m in medoids:
                        dcur[m] = -np.inf
                    repl = int(np.argmax(dcur))
                    new_medoids[cid] = repl
                    continue
                sub = dist[np.ix_(members, members)]
                costs = np.sum(sub, axis=1)
                best_local = int(np.argmin(costs))
                new_medoids[cid] = int(members[best_local])
            new_medoids = np.array(sorted(set(new_medoids.tolist())), dtype=int)
            if new_medoids.size < k:
                remaining = [i for i in range(n) if i not in new_medoids]
                for r in remaining:
                    if new_medoids.size >= k:
                        break
                    new_medoids = np.append(new_medoids, int(r))
            new_medoids = new_medoids[:k]
            if np.array_equal(np.sort(medoids), np.sort(new_medoids)):
                break
            medoids = new_medoids
        d_to_m = dist[:, medoids]
        labels = np.argmin(d_to_m, axis=1).astype(int)
        return labels

    @staticmethod
    def _silhouette_score(dist: np.ndarray, labels: np.ndarray) -> float:
        n = int(dist.shape[0])
        if n <= 2:
            return 0.0
        uniq = np.unique(labels)
        if uniq.size <= 1:
            return 0.0
        svals = np.zeros(n, dtype=float)
        for i in range(n):
            own = labels[i]
            own_idx = np.flatnonzero(labels == own)
            if own_idx.size <= 1:
                svals[i] = 0.0
                continue
            a = float(np.mean(dist[i, own_idx[own_idx != i]]))
            b = np.inf
            for other in uniq:
                if other == own:
                    continue
                other_idx = np.flatnonzero(labels == other)
                if other_idx.size == 0:
                    continue
                b = min(b, float(np.mean(dist[i, other_idx])))
            if not np.isfinite(b):
                svals[i] = 0.0
                continue
            den = max(a, b)
            svals[i] = 0.0 if den <= 0.0 else ((b - a) / den)
        return float(np.mean(svals))

    def _cluster_from_dist(self, dist: np.ndarray, sample_labels: Sequence[str]) -> Tuple[np.ndarray, float]:
        n = int(dist.shape[0])
        if n <= 2:
            return np.zeros(n, dtype=int), 0.0
        k_max = min(8, n - 1)
        best_labels = np.zeros(n, dtype=int)
        best_score = -1.0
        for k in range(2, k_max + 1):
            labels = self._kmedoids(dist, k)
            score = self._silhouette_score(dist, labels)
            if score > best_score + 1e-9:
                best_score = score
                best_labels = labels
        if best_score < 0:
            best_score = 0.0
        return best_labels.astype(int), float(best_score)

    @staticmethod
    def _ordered_cluster_groups(labels: np.ndarray, sample_labels: Sequence[str]) -> List[Tuple[int, List[int]]]:
        groups: Dict[int, List[int]] = {}
        for idx, cid in enumerate(labels.tolist()):
            groups.setdefault(int(cid), []).append(int(idx))
        ordered = sorted(
            groups.items(),
            key=lambda kv: (
                -len(kv[1]),
                tuple(sorted(str(sample_labels[i]).lower() for i in kv[1])),
            ),
        )
        out: List[Tuple[int, List[int]]] = []
        for cid, members in ordered:
            out.append((int(cid), sorted(members, key=lambda i: str(sample_labels[i]).lower())))
        return out

    def _ordered_cluster_colors(
        self,
        labels: np.ndarray,
        sample_labels: Sequence[str],
    ) -> Dict[int, QtGui.QColor]:
        n = len(sample_labels)
        if n == 0:
            return {}
        ordered = self._ordered_cluster_groups(labels, sample_labels)
        out: Dict[int, QtGui.QColor] = {}
        for color_idx, (_cid, members) in enumerate(ordered):
            col = QtGui.QColor(self._CLUSTER_PALETTE[color_idx % len(self._CLUSTER_PALETTE)])
            for m in members:
                out[m] = col
        for i in range(n):
            out.setdefault(i, QtGui.QColor("#64748b"))
        return out

    @staticmethod
    def _spot_index(point: object) -> Optional[int]:
        try:
            data = point.data()
        except Exception:
            return None
        if isinstance(data, dict) and "index" in data:
            try:
                return int(data.get("index"))
            except Exception:
                return None
        return None

    def _ensure_heatmap_highlight_items(self) -> None:
        if (
            self._heatmap_row_highlight is not None
            and self._heatmap_col_highlight is not None
            and self._heatmap_rows_multi_highlight is not None
            and self._heatmap_cols_multi_highlight is not None
        ):
            return
        vb = self._heatmap_plot.getViewBox()
        if vb is None:
            return
        rows_multi = QtWidgets.QGraphicsPathItem()
        rows_multi.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
        rows_multi.setBrush(QtGui.QBrush(QtGui.QColor(249, 115, 22, 64)))
        rows_multi.setZValue(19)
        rows_multi.setVisible(False)
        vb.addItem(rows_multi)
        cols_multi = QtWidgets.QGraphicsPathItem()
        cols_multi.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
        cols_multi.setBrush(QtGui.QBrush(QtGui.QColor(249, 115, 22, 64)))
        cols_multi.setZValue(19)
        cols_multi.setVisible(False)
        vb.addItem(cols_multi)
        row = QtWidgets.QGraphicsRectItem()
        row.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
        row.setBrush(QtGui.QBrush(QtGui.QColor(15, 23, 42, 72)))
        row.setZValue(20)
        row.setVisible(False)
        vb.addItem(row)
        col = QtWidgets.QGraphicsRectItem()
        col.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
        col.setBrush(QtGui.QBrush(QtGui.QColor(15, 23, 42, 72)))
        col.setZValue(20)
        col.setVisible(False)
        vb.addItem(col)
        self._heatmap_row_highlight = row
        self._heatmap_col_highlight = col
        self._heatmap_rows_multi_highlight = rows_multi
        self._heatmap_cols_multi_highlight = cols_multi

    def _set_heatmap_multi_highlight(self, scatter_indices: Optional[Set[int]]) -> None:
        self._ensure_heatmap_highlight_items()
        if self._heatmap_rows_multi_highlight is None or self._heatmap_cols_multi_highlight is None:
            return
        if not scatter_indices or self._matrix_sample_count <= 0:
            self._heatmap_rows_multi_highlight.setVisible(False)
            self._heatmap_cols_multi_highlight.setVisible(False)
            return
        display_indices: Set[int] = set()
        for scatter_idx in scatter_indices:
            sample = self._scatter_sample_by_index.get(int(scatter_idx), "")
            if not sample:
                continue
            display_idx = self._matrix_display_index_by_sample.get(sample)
            if display_idx is None:
                continue
            display_indices.add(int(display_idx))
        if not display_indices:
            self._heatmap_rows_multi_highlight.setVisible(False)
            self._heatmap_cols_multi_highlight.setVisible(False)
            return
        n = float(self._matrix_sample_count)
        row_path = QtGui.QPainterPath()
        col_path = QtGui.QPainterPath()
        for display_idx in sorted(display_indices):
            x = float(display_idx)
            y = float((self._matrix_sample_count - 1) - display_idx)
            row_path.addRect(QtCore.QRectF(0.0, y, n, 1.0))
            col_path.addRect(QtCore.QRectF(x, 0.0, 1.0, n))
        self._heatmap_rows_multi_highlight.setPath(row_path)
        self._heatmap_cols_multi_highlight.setPath(col_path)
        self._heatmap_rows_multi_highlight.setVisible(True)
        self._heatmap_cols_multi_highlight.setVisible(True)

    def _set_heatmap_cross_highlight(self, scatter_idx: Optional[int]) -> None:
        self._ensure_heatmap_highlight_items()
        if self._heatmap_row_highlight is None or self._heatmap_col_highlight is None:
            return
        if self._meta_hover_indices:
            self._heatmap_row_highlight.setVisible(False)
            self._heatmap_col_highlight.setVisible(False)
            self._set_heatmap_multi_highlight(self._meta_hover_indices)
            return
        self._set_heatmap_multi_highlight(None)
        if scatter_idx is None or self._matrix_sample_count <= 0:
            self._heatmap_row_highlight.setVisible(False)
            self._heatmap_col_highlight.setVisible(False)
            return
        sample = self._scatter_sample_by_index.get(scatter_idx, "")
        if not sample:
            self._heatmap_row_highlight.setVisible(False)
            self._heatmap_col_highlight.setVisible(False)
            return
        display_idx = self._matrix_display_index_by_sample.get(sample)
        if display_idx is None:
            self._heatmap_row_highlight.setVisible(False)
            self._heatmap_col_highlight.setVisible(False)
            return
        n = float(self._matrix_sample_count)
        x = float(display_idx)
        y = float((self._matrix_sample_count - 1) - display_idx)
        self._heatmap_row_highlight.setRect(0.0, y, n, 1.0)
        self._heatmap_col_highlight.setRect(x, 0.0, 1.0, n)
        self._heatmap_row_highlight.setVisible(True)
        self._heatmap_col_highlight.setVisible(True)

    def _apply_metadata_scatter_override(self) -> None:
        if not self._scatter_points_by_index:
            self._set_heatmap_multi_highlight(None)
            return
        if not self._meta_hover_indices:
            self._set_heatmap_multi_highlight(None)
            return
        for i, point in self._scatter_points_by_index.items():
            col = self._scatter_color_by_index.get(i, QtGui.QColor("#64748b"))
            selected = i in self._meta_hover_indices
            if selected:
                point.setBrush(pg.mkBrush(col.red(), col.green(), col.blue(), 255))
                point.setPen(pg.mkPen(QtGui.QColor("#0f172a"), width=2.1))
                point.setSize(20)
            else:
                point.setBrush(pg.mkBrush(col.red(), col.green(), col.blue(), 35))
                point.setPen(pg.mkPen(QtGui.QColor("#cbd5e1"), width=0.8))
                point.setSize(8)
        if self._scatter_hover_text is not None:
            self._scatter_hover_text.setVisible(False)
        self._set_heatmap_cross_highlight(None)

    def _set_scatter_hover(self, idx: Optional[int]) -> None:
        self._hover_scatter_index = idx
        if not self._scatter_points_by_index:
            self._set_heatmap_cross_highlight(None)
            return

        if idx is None:
            for i, point in self._scatter_points_by_index.items():
                col = self._scatter_color_by_index.get(i, QtGui.QColor("#64748b"))
                point.setBrush(pg.mkBrush(col.red(), col.green(), col.blue(), 195))
                point.setPen(pg.mkPen(col.darker(130), width=1.0))
                point.setSize(8)
            if self._scatter_hover_text is not None:
                self._scatter_hover_text.setVisible(False)
            self._set_heatmap_cross_highlight(None)
            self._apply_metadata_scatter_override()
            return

        cluster_id = self._scatter_cluster_by_index.get(idx)
        for i, point in self._scatter_points_by_index.items():
            col = self._scatter_color_by_index.get(i, QtGui.QColor("#64748b"))
            same_cluster = cluster_id is not None and self._scatter_cluster_by_index.get(i) == cluster_id
            if same_cluster:
                alpha = 242
                size = 13
                pen = pg.mkPen(col.darker(145), width=1.4)
            else:
                alpha = 45
                size = 8
                pen = pg.mkPen(QtGui.QColor("#94a3b8"), width=0.8)
            if i == idx:
                alpha = 255
                size = 18
                pen = pg.mkPen(QtGui.QColor("#0f172a"), width=2.0)
            point.setBrush(pg.mkBrush(col.red(), col.green(), col.blue(), alpha))
            point.setPen(pen)
            point.setSize(size)
        self._set_heatmap_cross_highlight(idx)

        if self._scatter_hover_text is not None:
            sample = self._scatter_sample_by_index.get(idx, "")
            xy = self._scatter_coords_by_index.get(idx)
            if sample and xy is not None:
                text = html.escape(sample)
                self._scatter_hover_text.setHtml(
                    "<div style='background:rgba(255,255,255,220); border:1px solid #0f172a; "
                    "padding:2px 4px; border-radius:4px; color:#0f172a; font-size:11px;'>"
                    f"{text}</div>"
                )
                self._scatter_hover_text.setPos(float(xy[0]), float(xy[1]))
                self._scatter_hover_text.setVisible(True)
            else:
                self._scatter_hover_text.setVisible(False)
        self._apply_metadata_scatter_override()

    def _on_scatter_mouse_moved(self, evt: object) -> None:
        if pg is None or self._scatter_item is None:
            return
        if self._meta_hover_indices:
            return
        try:
            pos = evt[0] if isinstance(evt, (tuple, list)) else evt
        except Exception:
            pos = evt
        if pos is None:
            self._set_scatter_hover(None)
            return
        if not self._scatter_plot.sceneBoundingRect().contains(pos):
            self._set_scatter_hover(None)
            return
        vb = self._scatter_plot.getViewBox()
        if vb is None:
            self._set_scatter_hover(None)
            return
        mouse_point = vb.mapSceneToView(pos)
        points = self._scatter_item.pointsAt(mouse_point)
        if points is None:
            self._set_scatter_hover(None)
            return
        try:
            points_count = len(points)
        except Exception:
            points_count = int(np.asarray(points).size)
        if points_count == 0:
            self._set_scatter_hover(None)
            return
        best_idx: Optional[int] = None
        best_d2 = float("inf")
        mx = float(mouse_point.x())
        my = float(mouse_point.y())
        for point in points:
            idx = self._spot_index(point)
            if idx is None:
                continue
            xy = self._scatter_coords_by_index.get(idx)
            if xy is None:
                continue
            dx = float(xy[0]) - mx
            dy = float(xy[1]) - my
            d2 = (dx * dx) + (dy * dy)
            if d2 < best_d2:
                best_d2 = d2
                best_idx = idx
        self._set_scatter_hover(best_idx)

    def refresh(self) -> None:
        if pg is None:
            return
        self._heatmap_plot.clear()
        self._scatter_plot.clear()
        self._scatter_item = None
        self._scatter_hover_text = None
        self._hover_scatter_index = None
        self._heatmap_row_highlight = None
        self._heatmap_col_highlight = None
        self._heatmap_rows_multi_highlight = None
        self._heatmap_cols_multi_highlight = None
        self._matrix_sample_count = 0
        self._sample_index_by_name = {}
        self._matrix_display_index_by_sample = {}
        self._scatter_points_by_index = {}
        self._scatter_sample_by_index = {}
        self._scatter_cluster_by_index = {}
        self._scatter_color_by_index = {}
        self._scatter_coords_by_index = {}
        self._heatmap_img = pg.ImageItem()
        self._heatmap_plot.addItem(self._heatmap_img)

        out = _find_panaroo_out(self._output_dir)
        if out is None:
            self._info.setText("Panaroo output not found.")
            return
        try:
            presence, sample_labels, _, _ = _load_presence_matrix(out)
        except Exception as exc:
            self._info.setText(f"Failed to load gene presence matrix: {exc}")
            return
        presence, sample_labels = _select_samples_from_presence(presence, sample_labels, self._sample_filter)
        self._visible_similarity_samples = list(sample_labels)
        self._refresh_meta_columns()
        if presence.size == 0:
            self._info.setText("No selected samples available for similarity view.")
            return

        # Filter ultra-rare and ultra-common loci for a more informative distance map.
        # Keep a safe fallback to all loci if filtering leaves too few genes.
        total_genes = int(presence.shape[0])
        prevalence = presence.mean(axis=1, dtype=float)
        keep = (prevalence >= self._PREVALENCE_MIN) & (prevalence <= self._PREVALENCE_MAX)
        filtered_gene_count = int(np.count_nonzero(keep))
        used_prevalence_filter = filtered_gene_count >= self._MIN_GENES_AFTER_PREVALENCE_FILTER
        used_presence = presence[keep] if used_prevalence_filter else presence
        used_gene_count = int(used_presence.shape[0])

        # Use int32 to prevent uint8 overflow in intersection counts.
        sample_gene = used_presence.T.astype(np.int32, copy=False)
        inter = sample_gene @ sample_gene.T
        sums = sample_gene.sum(axis=1).astype(float)
        unions = sums[:, None] + sums[None, :] - inter
        sim = np.divide(inter, unions, out=np.ones_like(unions, dtype=float), where=unions > 0)
        dist = 1.0 - sim

        cmap = pg.ColorMap(
            [0.0, 0.5, 1.0],
            [(30, 64, 175), (56, 189, 248), (250, 204, 21)],
        )

        coords = self._classical_mds(dist)
        labels, silhouette = self._cluster_from_dist(dist, sample_labels)
        color_map = self._ordered_cluster_colors(labels, sample_labels)
        ordered_groups = self._ordered_cluster_groups(labels, sample_labels)
        matrix_order = [idx for _cid, members in ordered_groups for idx in members]
        if len(matrix_order) != len(sample_labels):
            matrix_order = list(range(len(sample_labels)))
        sim_display = sim[np.ix_(matrix_order, matrix_order)]
        labels_display = [sample_labels[i] for i in matrix_order]
        self._matrix_sample_count = len(labels_display)
        self._sample_index_by_name = {str(name): int(i) for i, name in enumerate(sample_labels)}
        self._matrix_display_index_by_sample = {str(name): int(i) for i, name in enumerate(labels_display)}

        sim_levels = self._safe_levels(sim_display, (0.0, 1.0))
        self._heatmap_img.setLookupTable(cmap.getLookupTable())
        self._heatmap_img.setImage(np.rot90(sim_display, k=1), levels=sim_levels, autoLevels=False)
        self._set_ticks(self._heatmap_plot, labels_display)
        self._apply_square_heatmap_view(len(labels_display))
        self._ensure_heatmap_highlight_items()
        spots = []
        for i in range(len(sample_labels)):
            col = color_map.get(i, QtGui.QColor("#64748b"))
            spots.append(
                {
                    "pos": (float(coords[i, 0]), float(coords[i, 1])),
                    "data": {"index": int(i), "sample": sample_labels[i], "cluster": int(labels[i])},
                    "brush": pg.mkBrush(col.red(), col.green(), col.blue(), 195),
                    "pen": pg.mkPen(col.darker(130)),
                    "size": 11,
                }
            )
        scatter = pg.ScatterPlotItem()
        scatter.addPoints(spots)
        self._scatter_plot.addItem(scatter)
        self._scatter_item = scatter
        for point in scatter.points():
            idx = self._spot_index(point)
            if idx is None:
                continue
            self._scatter_points_by_index[idx] = point
            self._scatter_sample_by_index[idx] = str(sample_labels[idx])
            self._scatter_cluster_by_index[idx] = int(labels[idx])
            self._scatter_color_by_index[idx] = color_map.get(idx, QtGui.QColor("#64748b"))
            self._scatter_coords_by_index[idx] = (float(coords[idx, 0]), float(coords[idx, 1]))

        self._scatter_hover_text = pg.TextItem(anchor=(0.0, 1.4))
        self._scatter_hover_text.setZValue(30)
        self._scatter_hover_text.setVisible(False)
        self._scatter_plot.addItem(self._scatter_hover_text)
        self.highlight_sample(self._highlight_sample)
        self._apply_metadata_scatter_override()

        tri = np.triu_indices_from(sim, k=1)
        svals = sim[tri]
        dvals = dist[tri]
        mean_s = float(np.nanmean(svals)) if svals.size else 1.0
        mean_d = float(np.nanmean(dvals)) if dvals.size else 0.0
        n_clusters = int(len(np.unique(labels))) if labels.size else 0
        filter_note = (
            f"prevalence filter: {int(self._PREVALENCE_MIN * 100)}-{int(self._PREVALENCE_MAX * 100)}% "
            f"(kept {used_gene_count}/{total_genes})"
            if used_prevalence_filter
            else (
                f"prevalence filter skipped (kept {filtered_gene_count}/{total_genes} < "
                f"{self._MIN_GENES_AFTER_PREVALENCE_FILTER})"
            )
        )
        self._info.setText(
            f"{out / 'gene_presence_absence.csv'} | samples: {len(sample_labels)} | genes: {used_gene_count}/{total_genes} | "
            f"mean pairwise Jaccard similarity: {mean_s:.4f} | mean pairwise Jaccard distance: {mean_d:.4f} | "
            f"clusters: {n_clusters} | silhouette: {silhouette:.3f} | {filter_note}"
        )



class StructPresenceHeatmapTab(QtWidgets.QWidget):
    def __init__(self, output_dir: Path, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self._output_dir = Path(output_dir)
        self._sample_filter: Optional[set[str]] = None
        self._visible_sample_cols: List[str] = []
        self._visible_variant_count: int = 0
        self._highlight_sample: Optional[str] = None
        self._highlight_samples: List[str] = []
        self._highlight_item: Optional[QtWidgets.QGraphicsRectItem] = None
        self._multi_highlight_item: Optional[QtWidgets.QGraphicsPathItem] = None
        self._struct_df_cache_sig: Optional[Tuple[str, int, int]] = None
        self._struct_df_cache: Optional[pd.DataFrame] = None
        layout = QtWidgets.QVBoxLayout(self)
        self._info = QtWidgets.QLabel()
        self._info.setWordWrap(True)
        layout.addWidget(self._info)
        if pg is None:
            layout.addWidget(_PgUnavailableWidget("PyQtGraph is required for the structural-variant heatmap view."))
            return
        self._plot = pg.PlotWidget()
        self._plot.setBackground("#ffffff")
        self._plot.getPlotItem().setMenuEnabled(False)
        self._plot.setLabel("left", "Structural variants")
        self._plot.setLabel("bottom", "Samples")
        self._plot.invertY(True)
        self._plot.setMouseEnabled(x=True, y=False)
        self._plot.getAxis("bottom").setStyle(showValues=False, hideOverlappingLabels=True, autoExpandTextSpace=True)
        self._plot.getAxis("bottom").setHeight(12)
        self._image = pg.ImageItem(axisOrder="row-major")
        self._plot.addItem(self._image)
        self._ensure_highlight_item()
        layout.addWidget(self._plot, stretch=1)
        self.refresh()

    def set_sample_filter(self, samples: Optional[Sequence[str]]) -> None:
        next_filter = {str(s).strip() for s in samples if str(s).strip()} if samples is not None else None
        if next_filter == self._sample_filter:
            return
        self._sample_filter = next_filter
        self.refresh()

    def highlight_sample(self, sample: Optional[str]) -> None:
        value = str(sample).strip() if sample is not None else ""
        next_sample = value or None
        if self._highlight_sample == next_sample:
            return
        self._highlight_sample = next_sample
        if next_sample is not None:
            self._highlight_samples = []
        self._apply_highlight()

    def highlight_samples(self, samples: Optional[Sequence[str]]) -> None:
        if not samples:
            self._highlight_samples = []
            self._apply_highlight()
            return
        seen = set()
        ordered: List[str] = []
        for raw in samples:
            s = str(raw).strip()
            if not s or s in seen:
                continue
            seen.add(s)
            ordered.append(s)
        self._highlight_samples = ordered
        if ordered:
            self._highlight_sample = None
        self._apply_highlight()

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        super().resizeEvent(event)
        if pg is not None:
            self._update_bottom_ticks(self._visible_sample_cols)

    @staticmethod
    def _order_samples_by_similarity(mat: np.ndarray, sample_cols: Sequence[str]) -> Tuple[np.ndarray, List[str]]:
        n = len(sample_cols)
        if n <= 2 or mat.size == 0:
            return mat, [str(s) for s in sample_cols]

        # Use int32 to avoid overflow when counting shared features.
        sample_feature = (mat.T > 0).astype(np.int32)
        inter = sample_feature @ sample_feature.T
        sums = sample_feature.sum(axis=1).astype(float)
        unions = sums[:, None] + sums[None, :] - inter
        sim = np.divide(inter, unions, out=np.zeros_like(unions, dtype=float), where=unions > 0)
        np.fill_diagonal(sim, 0.0)

        degree = sim.sum(axis=1)
        laplacian = np.diag(degree) - sim
        try:
            eigvals, eigvecs = np.linalg.eigh(laplacian)
            if eigvecs.shape[1] >= 2:
                key = eigvecs[:, 1]
            else:
                key = eigvecs[:, 0]
        except np.linalg.LinAlgError:
            key = np.arange(n, dtype=float)

        order = np.argsort(key, kind="mergesort")
        return mat[:, order], [str(sample_cols[i]) for i in order]

    def _update_bottom_ticks(self, sample_cols: Sequence[str]) -> None:
        axis = self._plot.getAxis("bottom")
        axis.setTicks([[]])

    def _ensure_highlight_item(self) -> None:
        if self._highlight_item is not None and self._multi_highlight_item is not None:
            return
        if self._multi_highlight_item is None:
            multi = QtWidgets.QGraphicsPathItem()
            multi.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
            multi.setBrush(QtGui.QBrush(QtGui.QColor(251, 146, 60, 70)))
            multi.setZValue(19)
            multi.setVisible(False)
            self._plot.getViewBox().addItem(multi)
            self._multi_highlight_item = multi
        item = QtWidgets.QGraphicsRectItem()
        item.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
        item.setBrush(QtGui.QBrush(QtGui.QColor(15, 23, 42, 64)))
        item.setZValue(20)
        item.setVisible(False)
        self._plot.getViewBox().addItem(item)
        self._highlight_item = item

    def _apply_highlight(self) -> None:
        self._ensure_highlight_item()
        if self._highlight_item is None:
            return
        if self._multi_highlight_item is not None:
            if not self._highlight_samples or not self._visible_sample_cols or self._visible_variant_count <= 0:
                self._multi_highlight_item.setVisible(False)
            else:
                img_rect = self._image.boundingRect()
                n_samples = max(1, len(self._visible_sample_cols))
                cell_w = float(img_rect.width()) / float(n_samples)
                rows_h = float(img_rect.height())
                indices = sorted(
                    {
                        self._visible_sample_cols.index(name)
                        for name in self._highlight_samples
                        if name in self._visible_sample_cols
                    }
                )
                if not indices:
                    self._multi_highlight_item.setVisible(False)
                else:
                    p = QtGui.QPainterPath()
                    for col in indices:
                        x0 = float(img_rect.left()) + float(col) * cell_w
                        p.addRect(QtCore.QRectF(x0, float(img_rect.top()), cell_w, rows_h))
                    self._multi_highlight_item.setPath(p)
                    self._multi_highlight_item.setVisible(True)
        if not self._highlight_sample or not self._visible_sample_cols or self._visible_variant_count <= 0:
            self._highlight_item.setVisible(False)
            return
        try:
            col = self._visible_sample_cols.index(self._highlight_sample)
        except ValueError:
            self._highlight_item.setVisible(False)
            return
        img_rect = self._image.boundingRect()
        n_samples = max(1, len(self._visible_sample_cols))
        cell_w = float(img_rect.width()) / float(n_samples)
        x0 = float(img_rect.left()) + float(col) * cell_w
        self._highlight_item.setRect(
            x0,
            float(img_rect.top()),
            cell_w,
            float(img_rect.height()),
        )
        self._highlight_item.setVisible(True)

    def _load_struct_df(self, path: Path) -> pd.DataFrame:
        sig = _file_signature(path)
        if self._struct_df_cache is None or self._struct_df_cache_sig != sig:
            self._struct_df_cache = read_csv_df(path, sep="\t")
            self._struct_df_cache_sig = sig
        return self._struct_df_cache

    def refresh(self) -> None:
        if pg is None:
            return
        self._plot.clear()
        self._image = pg.ImageItem(axisOrder="row-major")
        self._plot.addItem(self._image)
        self._highlight_item = None
        self._multi_highlight_item = None
        self._ensure_highlight_item()

        out = _find_panaroo_out(self._output_dir)
        if out is None:
            self._visible_sample_cols = []
            self._visible_variant_count = 0
            self._info.setText("Panaroo output not found.")
            self._apply_highlight()
            return
        path = out / "struct_presence_absence.Rtab"
        if not path.is_file():
            self._visible_sample_cols = []
            self._visible_variant_count = 0
            self._info.setText("Panaroo struct_presence_absence.Rtab not found.")
            self._apply_highlight()
            return
        try:
            df = self._load_struct_df(path)
        except Exception as exc:
            self._visible_sample_cols = []
            self._visible_variant_count = 0
            self._info.setText(f"Failed to read {path.name}: {exc}")
            self._apply_highlight()
            return
        if df.empty or len(df.columns) < 2:
            self._visible_sample_cols = []
            self._visible_variant_count = 0
            self._info.setText(f"{path.name} has no matrix content.")
            self._apply_highlight()
            return

        sample_cols = [c for c in df.columns if c != "Gene"]
        if self._sample_filter is not None:
            sample_cols = [c for c in sample_cols if str(c).strip() in self._sample_filter]
        if not sample_cols:
            self._visible_sample_cols = []
            self._visible_variant_count = 0
            self._info.setText("No selected samples available for structural-variant view.")
            self._apply_highlight()
            return
        mat = df[sample_cols].fillna(0).astype(float).to_numpy()
        mat, sample_cols = self._order_samples_by_similarity(mat, sample_cols)
        if mat.shape[0] > 4000:
            step = int(math.ceil(mat.shape[0] / 4000.0))
            mat = mat[::step, :]
            note = f"rows downsampled x{step}"
        else:
            note = ""

        cmap = pg.ColorMap([0.0, 1.0], [(245, 245, 245), (2, 132, 199)])
        self._image.setLookupTable(cmap.getLookupTable())
        self._image.setImage(mat, levels=(0, 1), autoLevels=False)
        self._visible_sample_cols = list(sample_cols)
        self._visible_variant_count = int(mat.shape[0])
        self._update_bottom_ticks(self._visible_sample_cols)
        n_samples = mat.shape[1]
        n_variants = mat.shape[0]
        self._plot.setXRange(-0.5, max(0.5, n_samples - 0.5), padding=0.0)
        self._plot.setYRange(-0.5, max(0.5, n_variants - 0.5), padding=0.0)
        self._plot.setLimits(
            xMin=-0.5,
            xMax=max(0.5, n_samples - 0.5),
            yMin=-0.5,
            yMax=max(0.5, n_variants - 0.5),
            minYRange=max(1.0, float(n_variants)),
            maxYRange=max(1.0, float(n_variants)),
        )
        self._info.setText(
            f"{path} | structural features: {df.shape[0]} | samples: {len(sample_cols)}"
            + (f" | {note}" if note else "")
            + " | ordered by SV similarity"
        )
        self._apply_highlight()


class EntropyTrackTab(QtWidgets.QWidget):
    def __init__(self, output_dir: Path, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self._output_dir = Path(output_dir)
        self._sample_filter: Optional[set[str]] = None
        self._gene_query: str = ""
        self._entropy_df_cache_sig: Optional[Tuple[str, int, int]] = None
        self._entropy_df_cache: Optional[pd.DataFrame] = None
        self._gpa_df_cache_sig: Optional[Tuple[str, int, int]] = None
        self._gpa_df_cache: Optional[pd.DataFrame] = None
        layout = QtWidgets.QVBoxLayout(self)
        self._info = QtWidgets.QLabel()
        self._info.setWordWrap(True)
        layout.addWidget(self._info)
        if pg is None:
            layout.addWidget(_PgUnavailableWidget("PyQtGraph is required for the entropy track view."))
            return

        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal, self)
        split.setChildrenCollapsible(False)
        layout.addWidget(split, stretch=1)

        self._plot = pg.PlotWidget()
        self._plot.setBackground("#ffffff")
        self._plot.setLabel("left", "Entropy")
        self._plot.setLabel("bottom", "Displayed locus rank")
        self._plot.showGrid(x=False, y=True, alpha=0.25)
        split.addWidget(self._plot)

        self._top = QtWidgets.QTableWidget(0, 6)
        self._top.setHorizontalHeaderLabels(
            ["Cluster", "Gene symbol(s)", "Entropy", "Rank", "Prevalence", "Annotation"]
        )
        self._top.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self._top.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self._top.horizontalHeader().setSectionResizeMode(2, QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self._top.horizontalHeader().setSectionResizeMode(3, QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self._top.horizontalHeader().setSectionResizeMode(4, QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self._top.horizontalHeader().setSectionResizeMode(5, QtWidgets.QHeaderView.ResizeMode.Stretch)
        self._top.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self._top.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self._top.verticalHeader().setVisible(False)
        split.addWidget(self._top)
        split.setSizes([760, 560])
        self.refresh()

    def set_sample_filter(self, samples: Optional[Sequence[str]]) -> None:
        next_filter = {str(s).strip() for s in samples if str(s).strip()} if samples is not None else None
        if next_filter == self._sample_filter:
            return
        self._sample_filter = next_filter
        self.refresh()

    def set_gene_query(self, query: Optional[str]) -> None:
        next_query = str(query or "").strip()
        if next_query == self._gene_query:
            return
        self._gene_query = next_query
        self.refresh()

    def _load_entropy_df(self, path: Path) -> pd.DataFrame:
        sig = _file_signature(path)
        if self._entropy_df_cache is None or self._entropy_df_cache_sig != sig:
            self._entropy_df_cache = read_csv_df(path, header=None, names=["Locus", "Entropy"], usecols=[0, 1])
            self._entropy_df_cache_sig = sig
        return self._entropy_df_cache

    def _load_gpa_df(self, path: Path) -> pd.DataFrame:
        sig = _file_signature(path)
        if self._gpa_df_cache is None or self._gpa_df_cache_sig != sig:
            self._gpa_df_cache = read_csv_df(path, low_memory=False)
            self._gpa_df_cache_sig = sig
        return self._gpa_df_cache

    def refresh(self) -> None:
        if pg is None:
            return
        self._plot.clear()
        self._top.setRowCount(0)

        out = _find_panaroo_out(self._output_dir)
        if out is None:
            self._info.setText("Panaroo output not found.")
            return
        path = out / "alignment_entropy.csv"
        if not path.is_file():
            self._info.setText("Panaroo alignment_entropy.csv not found.")
            return
        try:
            df = self._load_entropy_df(path).copy()
        except Exception as exc:
            self._info.setText(f"Failed to read {path.name}: {exc}")
            return
        df["Locus"] = df["Locus"].astype(str)
        df["Cluster"] = df["Locus"].str.replace(r"\.aln(?:\.fas)?$", "", regex=True)
        df["Entropy"] = pd.to_numeric(df["Entropy"], errors="coerce")
        df = df.dropna(subset=["Entropy"]).reset_index(drop=True)
        if df.empty:
            self._info.setText("No numeric entropy values found.")
            return

        gpa_path = out / "gene_presence_absence.csv"
        if not gpa_path.is_file():
            self._info.setText("Panaroo gene_presence_absence.csv not found.")
            return
        try:
            gpa = self._load_gpa_df(gpa_path)
        except Exception as exc:
            self._info.setText(f"Failed to read {gpa_path.name}: {exc}")
            return

        sample_cols = [c for c in gpa.columns if c not in _META_COLUMNS]
        if self._sample_filter is not None:
            sample_cols = [c for c in sample_cols if str(c).strip() in self._sample_filter]
        prev_denom = len(sample_cols)
        if prev_denom > 0:
            presence = gpa[sample_cols].fillna("").astype(str)
            present = (presence.ne("")) & (presence.ne("0"))
            prevalence_count = present.sum(axis=1).astype(int)
        else:
            prevalence_count = pd.Series(0, index=gpa.index, dtype=int)

        gene_col = gpa["Gene"].astype(str) if "Gene" in gpa.columns else pd.Series("", index=gpa.index, dtype=str)
        symbols_col = (
            gpa["Non-unique Gene name"].fillna("").astype(str)
            if "Non-unique Gene name" in gpa.columns
            else pd.Series("", index=gpa.index, dtype=str)
        )
        annotation_col = (
            gpa["Annotation"].fillna("").astype(str)
            if "Annotation" in gpa.columns
            else pd.Series("", index=gpa.index, dtype=str)
        )
        meta = pd.DataFrame(
            {
                "Cluster": gene_col,
                "GeneSymbols": symbols_col,
                "Annotation": annotation_col,
                "PrevalenceCount": prevalence_count.astype(int),
            }
        )
        meta = meta.drop_duplicates(subset=["Cluster"], keep="first").reset_index(drop=True)
        merged = df.merge(meta, on="Cluster", how="left")
        merged["GeneSymbols"] = merged["GeneSymbols"].fillna("").astype(str)
        merged["Annotation"] = merged["Annotation"].fillna("").astype(str)
        merged["Rank"] = (
            merged["Entropy"].rank(method="first", ascending=False).astype(int)
        )
        merged = merged.sort_values("Entropy", ascending=False).reset_index(drop=True)

        all_y = merged["Entropy"].to_numpy(dtype=float)
        all_x = np.arange(1, len(all_y) + 1, dtype=float)
        curve = pg.PlotCurveItem(x=all_x, y=all_y, pen=pg.mkPen("#7c3aed", width=1.6))
        self._plot.addItem(curve)
        dots = pg.ScatterPlotItem(
            x=all_x,
            y=all_y,
            pen=pg.mkPen(91, 33, 182, 70),
            brush=pg.mkBrush(124, 58, 237, 50),
            size=3,
        )
        self._plot.addItem(dots)
        self._plot.setXRange(1, float(max(2, len(all_y))), padding=0.02)
        self._plot.setYRange(0, float(max(0.0001, all_y.max())) * 1.05, padding=0.02)

        query = self._gene_query.strip()
        if query:
            q_norm = re.sub(r"[^a-z0-9]+", "", query.lower())
            meta_cluster_norm = meta["Cluster"].astype(str).str.lower().str.replace(r"[^a-z0-9]+", "", regex=True)
            meta_symbols_norm = meta["GeneSymbols"].astype(str).str.lower().str.replace(r"[^a-z0-9]+", "", regex=True)
            meta_annotation_norm = meta["Annotation"].astype(str).str.lower().str.replace(r"[^a-z0-9]+", "", regex=True)
            meta_hit_mask = (
                meta["Cluster"].str.contains(query, case=False, na=False, regex=False)
                | meta["GeneSymbols"].str.contains(query, case=False, na=False, regex=False)
                | meta["Annotation"].str.contains(query, case=False, na=False, regex=False)
            )
            if q_norm:
                meta_hit_mask = (
                    meta_hit_mask
                    | meta_cluster_norm.str.contains(q_norm, na=False, regex=False)
                    | meta_symbols_norm.str.contains(q_norm, na=False, regex=False)
                    | meta_annotation_norm.str.contains(q_norm, na=False, regex=False)
                )
            subset = meta[meta_hit_mask].copy()
            subset = subset.merge(
                merged[["Cluster", "Entropy", "Rank"]].copy(),
                on="Cluster",
                how="left",
            )
            subset = subset.sort_values(["Entropy", "Cluster"], ascending=[False, True], na_position="last").reset_index(drop=True)

            hit_subset = subset.dropna(subset=["Entropy"])
            if not hit_subset.empty:
                hit_x = hit_subset["Rank"].to_numpy(dtype=float)
                hit_y = hit_subset["Entropy"].to_numpy(dtype=float)
                hit_dots = pg.ScatterPlotItem(
                    x=hit_x,
                    y=hit_y,
                    pen=pg.mkPen("#9a3412"),
                    brush=pg.mkBrush(249, 115, 22, 210),
                    size=8,
                )
                self._plot.addItem(hit_dots)
            mode_text = f"query='{query}'"
        else:
            subset = merged.head(60).copy()
            mode_text = "top high-entropy loci"

        self._top.setRowCount(len(subset))
        for i, row in subset.iterrows():
            prevalence = (
                f"{int(row['PrevalenceCount'])}/{prev_denom} ({(100.0 * float(row['PrevalenceCount']) / float(prev_denom)):.1f}%)"
                if prev_denom > 0
                else "0/0"
            )
            entropy_v = row.get("Entropy")
            rank_v = row.get("Rank")
            self._top.setItem(i, 0, QtWidgets.QTableWidgetItem(str(row["Cluster"])))
            self._top.setItem(i, 1, QtWidgets.QTableWidgetItem(str(row["GeneSymbols"])))
            self._top.setItem(
                i,
                2,
                QtWidgets.QTableWidgetItem("N/A" if pd.isna(entropy_v) else f"{float(entropy_v):.6f}"),
            )
            self._top.setItem(
                i,
                3,
                QtWidgets.QTableWidgetItem("N/A" if pd.isna(rank_v) else str(int(rank_v))),
            )
            self._top.setItem(i, 4, QtWidgets.QTableWidgetItem(prevalence))
            self._top.setItem(i, 5, QtWidgets.QTableWidgetItem(str(row["Annotation"])))

        hits_info = ""
        if query:
            entropy_hits = int(subset["Entropy"].notna().sum()) if "Entropy" in subset.columns else 0
            hits_info = f" | matches: {len(subset)} | with entropy: {entropy_hits}"
        self._info.setText(
            f"{path} | loci: {len(merged)} | showing: {len(subset)} ({mode_text}) | "
            f"mean entropy: {float(merged['Entropy'].mean()):.6f} | max entropy: {float(merged['Entropy'].max()):.6f}"
            + hits_info
        )
