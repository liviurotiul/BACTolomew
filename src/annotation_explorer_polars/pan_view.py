from __future__ import annotations

import math
import time
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets

from .panaroo_style import _apply_analysis_panel_style, _log_view_event
from .polars_io import read_csv_df

try:
    import pyqtgraph as pg
except Exception:
    pg = None


_ROARY_META_COLUMNS = {
    "Gene",
    "Non-unique Gene name",
    "Annotation",
    "No. isolates",
    "No. sequences",
    "Avg sequences per isolate",
    "Min group size nuc",
    "Max group size nuc",
    "Avg group size nuc",
    "Genome Fragment",
    "Order within Fragment",
    "Accessory Fragment",
    "Accessory Order with Fragment",
    "QC",
}


class PanHeatmapWidget(QtWidgets.QWidget):
    viewChanged = QtCore.pyqtSignal()
    sampleHovered = QtCore.pyqtSignal(object)

    def __init__(
        self,
        output_dir: Path,
        parent: Optional[QtWidgets.QWidget] = None,
        *,
        auto_refresh: bool = True,
    ) -> None:
        super().__init__(parent)
        _apply_analysis_panel_style(self)
        self._output_dir = Path(output_dir)
        _log_view_event(f"created Panaroo gene presence view for {self._output_dir}")
        self._plot: Optional[QtWidgets.QWidget] = None
        self._image_item = None
        self._info = QtWidgets.QLabel()
        self._info.setWordWrap(True)
        self._data_loaded = False
        self._data_path: Optional[Path] = None
        self._presence: Optional[np.ndarray] = None
        self._gene_labels: List[str] = []
        self._annotation_labels: List[str] = []
        self._sample_labels: List[str] = []
        self._gene_series: Optional[pd.Series] = None
        self._annotation_series: Optional[pd.Series] = None
        self._filter_samples: Optional[List[str]] = None
        self._filter_gene_text: Optional[str] = None
        self._filter_product_text: Optional[str] = None
        self._sample_order: Optional[List[str]] = None
        self._highlight_sample: Optional[str] = None
        self._highlight_samples: List[str] = []
        self._hover_sample: Optional[str] = None
        self._highlight_item: Optional[QtWidgets.QGraphicsRectItem] = None
        self._multi_highlight_item = None
        self._display_sample_labels: List[str] = []
        self._display_shape: Tuple[int, int] = (0, 0)
        self._build_ui()
        if auto_refresh:
            self.refresh()

    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._info)
        self._info.setVisible(False)

        if pg is None:
            msg = QtWidgets.QLabel("PyQtGraph is not available. Install it to enable the pan heatmap.")
            msg.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
            msg.setWordWrap(True)
            layout.addWidget(msg, stretch=1)
            return

        pg.setConfigOptions(imageAxisOrder="row-major")
        plot = pg.PlotWidget()
        plot.setBackground("#ffffff")
        plot.showGrid(x=False, y=False)
        plot.setMouseEnabled(x=True, y=False)
        plot.invertY(True)
        plot.setMouseTracking(True)
        plot.viewport().setMouseTracking(True)
        plot.viewport().installEventFilter(self)
        plot.getPlotItem().setMenuEnabled(False)
        plot.getPlotItem().setLabel("bottom", "Genes")
        plot.getPlotItem().setLabel("left", "Samples")
        plot.getPlotItem().vb.sigRangeChanged.connect(lambda *_: self.viewChanged.emit())

        self._image_item = pg.ImageItem()
        plot.addItem(self._image_item)
        multi_highlight = QtWidgets.QGraphicsPathItem()
        multi_highlight.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
        multi_highlight.setBrush(QtGui.QBrush(QtGui.QColor(251, 146, 60, 72)))
        multi_highlight.setZValue(18)
        multi_highlight.setVisible(False)
        plot.getViewBox().addItem(multi_highlight)
        self._multi_highlight_item = multi_highlight

        highlight = QtWidgets.QGraphicsRectItem()
        highlight.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
        highlight.setBrush(QtGui.QBrush(QtGui.QColor(100, 181, 246, 120)))
        highlight.setZValue(20)
        highlight.setVisible(False)
        plot.getViewBox().addItem(highlight)
        self._highlight_item = highlight

        self._highlight_opacity = QtWidgets.QGraphicsOpacityEffect()
        highlight.setGraphicsEffect(self._highlight_opacity)
        self._highlight_opacity.setOpacity(0.0)
        self._highlight_anim = QtCore.QPropertyAnimation(self._highlight_opacity, b"opacity", self)
        self._highlight_anim.setDuration(180)
        self._highlight_anim.setEasingCurve(QtCore.QEasingCurve.Type.InOutQuad)
        layout.addWidget(plot, stretch=1)
        self._plot = plot

    def refresh(self) -> None:
        if pg is None or self._image_item is None:
            return
        path = self._find_gene_presence_absence()
        if path is None:
            print(f"[file] Panaroo gene_presence_absence.csv not found under {self._output_dir}")
            self._info.setText("Panaroo gene presence/absence data not found.")
            self._info.setVisible(True)
            self._image_item.setImage(np.zeros((1, 1)), levels=(0, 1), autoLevels=False)
            return

        if not self._data_loaded or self._data_path != path:
            self._load_data(path)

        if not self._data_loaded or self._presence is None:
            return

        self._apply_filters()

    def highlight_sample(self, sample: Optional[str]) -> None:
        self._highlight_sample = sample
        self._apply_highlight()

    def highlight_samples(self, samples: Optional[Sequence[str]]) -> None:
        if not samples:
            self._highlight_samples = []
        else:
            seen = set()
            ordered: List[str] = []
            for raw in samples:
                name = str(raw).strip()
                if not name or name in seen:
                    continue
                seen.add(name)
                ordered.append(name)
            self._highlight_samples = ordered
        self._apply_highlight()

    def eventFilter(self, obj: QtCore.QObject, event: QtCore.QEvent) -> bool:
        if self._plot is not None and obj is self._plot.viewport():
            if event.type() == QtCore.QEvent.Type.MouseMove:
                pos = event.position().toPoint() if hasattr(event, "position") else event.pos()
                scene_pos = self._plot.mapToScene(pos)
                self._update_hover_from_scene(scene_pos)
            elif event.type() == QtCore.QEvent.Type.Leave:
                self._set_hover_sample(None)
        return super().eventFilter(obj, event)

    def set_filters(
        self,
        allowed_samples: Optional[List[str]],
        gene_text: Optional[str],
        product_text: Optional[str],
    ) -> None:
        self._filter_samples = allowed_samples
        self._filter_gene_text = gene_text.strip() if gene_text else None
        self._filter_product_text = product_text.strip() if product_text else None
        self._apply_filters()

    def set_sample_order(self, sample_order: Optional[List[str]]) -> None:
        self._sample_order = list(sample_order) if sample_order else None
        self._apply_filters()

    def display_samples(self) -> List[str]:
        return list(self._display_sample_labels)

    def sample_viewport_y(self, sample: str) -> Optional[float]:
        if pg is None or self._plot is None:
            return None
        if not sample or sample not in self._display_sample_labels:
            return None
        row = self._display_sample_labels.index(sample)
        vb = self._plot.getPlotItem().vb
        scene_pt = vb.mapViewToScene(QtCore.QPointF(0.0, float(row) + 0.5))
        view_pt = self._plot.mapFromScene(scene_pt)
        return float(view_pt.y())

    def _find_gene_presence_absence(self) -> Optional[Path]:
        if not self._output_dir.exists():
            return None

        # Annotation explorer should display the Panaroo pangenome matrix.
        panaroo_primary = self._output_dir / "panaroo" / "out" / "gene_presence_absence.csv"
        if panaroo_primary.is_file():
            return panaroo_primary

        panaroo_root = self._output_dir / "panaroo"
        if not panaroo_root.exists():
            return None

        candidates = list(panaroo_root.rglob("gene_presence_absence.csv"))
        if not candidates:
            return None

        candidates.sort(key=lambda p: p.stat().st_mtime, reverse=True)
        return candidates[0]

    def _load_data(self, path: Path) -> None:
        self._data_loaded = False
        self._data_path = path
        t0 = time.perf_counter()
        print(f"[file] reading {path}")
        df = read_csv_df(path, low_memory=False)
        print(f"[file] read {path} in {time.perf_counter() - t0:.2f}s")
        if "Gene" not in df.columns:
            self._info.setText("Panaroo gene presence/absence data is missing the Gene column.")
            self._info.setVisible(True)
            return

        sample_cols = [c for c in df.columns if c not in _ROARY_META_COLUMNS]
        if not sample_cols:
            self._info.setText("Panaroo gene presence/absence data has no sample columns.")
            self._info.setVisible(True)
            return

        raw = df[sample_cols].fillna("").astype(str)
        presence = raw.ne("") & raw.ne("0")
        self._presence = presence.to_numpy(dtype=np.uint8)
        self._gene_labels = df["Gene"].astype(str).tolist()
        if "Annotation" in df.columns:
            self._annotation_labels = df["Annotation"].fillna("").astype(str).tolist()
        else:
            self._annotation_labels = [""] * len(self._gene_labels)
        self._sample_labels = [str(c) for c in sample_cols]
        self._gene_series = pd.Series(self._gene_labels)
        self._annotation_series = pd.Series(self._annotation_labels)
        self._data_loaded = True

    def _apply_filters(self) -> None:
        if self._presence is None:
            return
        genes_total = len(self._gene_labels)
        samples_total = len(self._sample_labels)

        gene_mask = pd.Series(True, index=range(genes_total))
        if self._filter_gene_text and self._gene_series is not None:
            gene_mask &= self._gene_series.str.contains(self._filter_gene_text, case=False, na=False)
        if self._filter_product_text and self._annotation_series is not None:
            gene_mask &= self._annotation_series.str.contains(self._filter_product_text, case=False, na=False)
        gene_idx = gene_mask[gene_mask].index.to_numpy()

        if self._filter_samples is None:
            sample_idx = np.arange(samples_total, dtype=int)
        else:
            allowed = set(self._filter_samples)
            sample_idx = np.array(
                [i for i, s in enumerate(self._sample_labels) if s in allowed],
                dtype=int,
            )

        if self._sample_order:
            rank = {name: i for i, name in enumerate(self._sample_order)}
            ordered = sorted(
                sample_idx.tolist(),
                key=lambda idx: (rank.get(self._sample_labels[idx], len(rank)), idx),
            )
            sample_idx = np.array(ordered, dtype=int)

        if gene_idx.size == 0 or sample_idx.size == 0:
            self._image_item.setImage(np.zeros((1, 1)), levels=(0, 1), autoLevels=False)
            self._lock_vertical_range(rows=1)
            self._apply_ticks([], [])
            self._info.setText("No genes or samples match the current filters.")
            self._info.setVisible(True)
            self.viewChanged.emit()
            return

        matrix = self._presence[np.ix_(gene_idx, sample_idx)]
        gene_labels = [self._gene_labels[i] for i in gene_idx]
        sample_labels = [self._sample_labels[i] for i in sample_idx]

        matrix, gene_labels, sample_labels, note = self._downsample(matrix, gene_labels, sample_labels)
        # Rotate 90° counter-clockwise: x=genes, y=samples (reversed order)
        matrix = np.rot90(matrix, k=1)
        x_labels = gene_labels
        y_labels = list(reversed(sample_labels))
        cmap = pg.ColorMap([0.0, 1.0], [(240, 240, 240), (30, 136, 229)])
        self._image_item.setLookupTable(cmap.getLookupTable())
        self._image_item.setImage(matrix, levels=(0, 1), autoLevels=False)
        self._lock_vertical_range(rows=matrix.shape[0])

        self._display_sample_labels = y_labels
        self._display_shape = matrix.shape
        self._apply_highlight()
        self._apply_ticks(x_labels, y_labels)
        if note:
            _log_view_event(f"Panaroo gene presence view downsampled: {note}")
        self._info.clear()
        self._info.setVisible(False)
        self.viewChanged.emit()

    def _sample_from_scene_pos(self, scene_pos: QtCore.QPointF) -> Optional[str]:
        if pg is None or self._plot is None:
            return None
        if not self._display_sample_labels:
            return None
        vb = self._plot.getPlotItem().vb
        if not vb.sceneBoundingRect().contains(scene_pos):
            return None
        view_pt = vb.mapSceneToView(scene_pos)
        row = int(math.floor(float(view_pt.y())))
        if row < 0 or row >= len(self._display_sample_labels):
            return None
        return self._display_sample_labels[row]

    def _set_hover_sample(self, sample: Optional[str]) -> None:
        if sample == self._hover_sample:
            return
        self._hover_sample = sample
        self.sampleHovered.emit(sample)

    def _update_hover_from_scene(self, scene_pos: QtCore.QPointF) -> None:
        sample = self._sample_from_scene_pos(scene_pos)
        self._set_hover_sample(sample)

    def _lock_vertical_range(self, rows: int) -> None:
        if pg is None or self._plot is None:
            return
        rows_f = float(max(1, rows))
        vb = self._plot.getPlotItem().vb
        vb.enableAutoRange(axis=pg.ViewBox.YAxis, enable=False)
        vb.setLimits(
            yMin=0.0,
            yMax=rows_f,
            minYRange=rows_f,
            maxYRange=rows_f,
        )
        vb.setYRange(0.0, rows_f, padding=0.0)

    def _downsample(
        self,
        matrix: np.ndarray,
        genes: List[str],
        samples: List[str],
        max_genes: int = 2000,
        max_samples: int = 2000,
    ) -> Tuple[np.ndarray, List[str], List[str], str]:
        n_genes, n_samples = matrix.shape
        gene_step = max(1, math.ceil(n_genes / max_genes))
        sample_step = max(1, math.ceil(n_samples / max_samples))
        note = ""
        if gene_step > 1 or sample_step > 1:
            matrix = matrix[::gene_step, ::sample_step]
            genes = genes[::gene_step]
            samples = samples[::sample_step]
            note = f"downsampled (genes x{gene_step}, samples x{sample_step})"
        return matrix, genes, samples, note

    def _apply_ticks(self, x_labels: List[str], y_labels: List[str]) -> None:
        if pg is None or self._plot is None:
            return
        plot = self._plot
        if len(x_labels) <= 50:
            x_ticks = [(i + 0.5, s) for i, s in enumerate(x_labels)]
            plot.getAxis("bottom").setTicks([x_ticks])
        else:
            plot.getAxis("bottom").setTicks([])

        if len(y_labels) <= 50:
            y_ticks = [(i + 0.5, g) for i, g in enumerate(y_labels)]
            plot.getAxis("left").setTicks([y_ticks])
        else:
            plot.getAxis("left").setTicks([])

    def _apply_highlight(self) -> None:
        if self._multi_highlight_item is not None:
            if not self._highlight_samples or not self._display_sample_labels:
                self._multi_highlight_item.setVisible(False)
            else:
                rows, cols = self._display_shape
                if rows == 0 or cols == 0:
                    self._multi_highlight_item.setVisible(False)
                else:
                    row_idx = sorted(
                        {
                            self._display_sample_labels.index(name)
                            for name in self._highlight_samples
                            if name in self._display_sample_labels
                        }
                    )
                    if not row_idx:
                        self._multi_highlight_item.setVisible(False)
                    else:
                        path = QtGui.QPainterPath()
                        for row in row_idx:
                            path.addRect(QtCore.QRectF(0.0, float(row), float(cols), 1.0))
                        self._multi_highlight_item.setPath(path)
                        self._multi_highlight_item.setVisible(True)

        if self._highlight_item is None:
            return
        if not self._highlight_sample or not self._display_sample_labels:
            self._highlight_item.setVisible(False)
            if hasattr(self, "_highlight_opacity"):
                self._highlight_opacity.setOpacity(0.0)
            return
        try:
            row = self._display_sample_labels.index(self._highlight_sample)
        except ValueError:
            self._highlight_item.setVisible(False)
            return
        rows, cols = self._display_shape
        if rows == 0 or cols == 0:
            self._highlight_item.setVisible(False)
            if hasattr(self, "_highlight_opacity"):
                self._highlight_opacity.setOpacity(0.0)
            return
        self._highlight_item.setRect(0, row, cols, 1)
        self._highlight_item.setVisible(True)
        if hasattr(self, "_highlight_anim"):
            self._highlight_anim.stop()
            self._highlight_anim.setStartValue(self._highlight_opacity.opacity())
            self._highlight_anim.setEndValue(1.0)
            self._highlight_anim.start()
