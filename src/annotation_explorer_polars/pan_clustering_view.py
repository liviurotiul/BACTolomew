"""The Pangenome clustering view: how a part of the pangenome makes the
isolates cluster.

The band at the bottom (``prevalence_band.py``) picks the genes by
prevalence; the plot shows the classical MDS of the Jaccard distances
between the isolates' profiles over those genes, in 2D (MDS1 against MDS2)
or in 3D (drawn here with plain pyqtgraph: ``scatter3d.py``), colored by a
metadata column or the k-means cluster and shaped by the k-means cluster.
The compute (``model/pan_clustering.py``) runs on the thread pool, cached by
sample set, range and k; colors, shapes and the 2D/3D switch never
recompute.
"""
from __future__ import annotations

import html
import threading
import time
from collections import OrderedDict
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from PyQt6 import QtCore, QtGui, QtWidgets

import pyqtgraph as pg

from .cgmlst_force_plot import _FORCE_SYMBOLS
from .model import pan_clustering as pc
from .panaroo_common import CLUSTER_HUES, _natural_sort_key, cluster_colors
from .panaroo_data import _load_presence_matrix
from .panaroo_force import GraphKeyEntry, _GraphKey
from .panaroo_style import _fastgear_lineage_colors
from .panaroo_widgets import _high_contrast_color
from .prevalence_band import PrevalenceRangeBand, image_from_shares, percent_text
from .refresh_worker import RefreshWorker
from .scatter3d import (
    CAMERA_DISTANCE,
    DEPTH_FADES,
    DEPTH_SIZES,
    VIEW_EXTENT,
    DepthScatterItem,
    OrbitViewBox,
    Projection3D,
    ReferenceBox,
    depth_levels,
    fade_toward,
    sphere_brush,
)
from .settings_column import (
    ElidedHint,
    SettingsColumn,
    WrapLayout,
    SettingsStrip,
    auto_fold,
    select_combo_data,
    settings_card,
    settings_combo,
    make_chip,
    settings_control,
    stacked_field,
)
from .theme import COLORS, out_of_set_color

KMEANS = "__kmeans__"
KMEANS_LABEL = "MDS cluster"
ISOLATES_WINDOW = "window"
ISOLATES_ALL = "all"
MODE_2D = "2d"
MODE_3D = "3d"
_POINT_PX = 10.0
# 3D points are shaded spheres, a bit bigger than the flat 2D markers.
_SPHERE_PX = 15.0
# Hovered isolates (here or in another panel) are redrawn in this color.
_HIGHLIGHT_COLOR = "#e6ff00"
_HOVER_PX = 8.0
_RANGE_DEBOUNCE_MS = 250
_FRAME_MS = 16
# A numeric column with more values than this is colored by quantile bins.
_NUMERIC_MAX_CATEGORIES = 12
_NUMERIC_BINS = 5
# Count labels on stacked points (isolates with identical profiles), at most.
_STACK_LABELS_MAX = 80
# The exported plots: this wide, as tall as the points need (within bounds).
_EXPORT_WIDTH = 1600
_EXPORT_HEIGHT = (900, 1500)
# The 3D reference box: tick marks and the axis titles' gap out of its
# outline, in pixels.
_BOX_TICK_PX = 5.0
_BOX_TITLE_GAP_PX = 7.0


def _trait_label(name: str) -> str:
    """The fastGEAR cluster is not the MDS cluster: say which one it is."""
    return "fastGEAR cluster" if name == "Cluster" else name


def _categorical_color(index: int, count: int) -> QtGui.QColor:
    """The app's cluster hues while they last (well apart in their first
    entries), then the golden angle colors."""
    if count <= len(CLUSTER_HUES):
        return QtGui.QColor(CLUSTER_HUES[index])
    return QtGui.QColor(_high_contrast_color(index))


def _number(value: str) -> Optional[float]:
    try:
        number = float(str(value).replace(",", ""))
    except ValueError:
        return None
    return number if np.isfinite(number) else None


def numeric_bins(values: Sequence[str]) -> Optional[Tuple[List[str], List[str]]]:
    """A continuous numeric column (every value a number, more than
    ``_NUMERIC_MAX_CATEGORIES`` of them) as quantile bins: each value's bin
    label ("" for no value) and the labels from low to high; None for a
    categorical column."""
    present = [v for v in values if v]
    if len(set(present)) <= _NUMERIC_MAX_CATEGORIES:
        return None
    numbers = [_number(v) for v in present]
    if any(n is None for n in numbers):
        return None
    data = np.asarray(numbers, dtype=float)
    edges = np.unique(np.quantile(data, np.linspace(0.0, 1.0, _NUMERIC_BINS + 1)))
    if edges.size < 2:
        return None

    def fmt(x: float) -> str:
        return f"{x:,.0f}" if abs(x) >= 1000 else f"{x:.3g}"

    labels = [f"{fmt(edges[i])} to {fmt(edges[i + 1])}" for i in range(edges.size - 1)]
    out = []
    for value in values:
        number = _number(value) if value else None
        if number is None:
            out.append("")
            continue
        index = int(np.clip(np.searchsorted(edges, number, side="right") - 1, 0, len(labels) - 1))
        out.append(labels[index])
    return out, labels


def ramp_colors(count: int) -> List[QtGui.QColor]:
    """Sequential colors, low to high (viridis without its palest end)."""
    cmap = pg.colormap.get("viridis")
    stops = np.linspace(0.0, 0.85, max(1, count))
    return [cmap.map(float(t), mode="qcolor") for t in stops]


class _Shared:
    """The run's presence matrix and the clustering caches, shared by every
    Pangenome clustering panel of the run (``PanData.memo``)."""

    def __init__(self, pan: Any) -> None:
        self._pan = pan
        self._lock = threading.Lock()
        self._matrix: Optional[Tuple[np.ndarray, List[str], List[str], List[str]]] = None
        self._slices: "OrderedDict[Any, pc.PresenceSlice]" = OrderedDict()
        self.cache = pc.ClusteringCache()

    def matrix(self) -> Tuple[np.ndarray, List[str], List[str], List[str]]:
        with self._lock:
            if self._matrix is None:
                out = self._pan.panaroo_out()
                if out is None:
                    raise RuntimeError("This run has no Panaroo output.")
                self._matrix = _load_presence_matrix(Path(out))
            return self._matrix

    def slice_for(self, key: Any, samples: Optional[Sequence[str]]) -> pc.PresenceSlice:
        with self._lock:
            cached = self._slices.get(key)
            if cached is not None:
                self._slices.move_to_end(key)
                return cached
        presence, names, genes, annotations = self.matrix()
        slice_ = pc.presence_slice(presence, names, genes, annotations, samples)
        with self._lock:
            self._slices[key] = slice_
            while len(self._slices) > 4:
                self._slices.popitem(last=False)
        return slice_


def shared_state(pan: Any) -> _Shared:
    return pan.memo(("pan_mds", "shared"), lambda: _Shared(pan))


class PanClusteringView(QtWidgets.QWidget):
    """Settings column, plot (2D or 3D), key, hover card and band."""

    samplesHovered = QtCore.pyqtSignal(object)
    sampleOpenRequested = QtCore.pyqtSignal(str)
    # ({"name", "samples", "label"}, in a new window)
    subsetRequested = QtCore.pyqtSignal(object, bool)
    computed = QtCore.pyqtSignal()

    _GRAPH_MIN_WIDTH = 420

    def __init__(
        self,
        pan: Any,
        *,
        traits: Optional[Mapping[str, Mapping[str, str]]] = None,
        hover_metadata: Optional[Mapping[str, Mapping[str, str]]] = None,
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self.setObjectName("PanClusteringView")
        self._pan = pan
        self._shared = shared_state(pan)
        self._traits: Dict[str, Dict[str, str]] = {str(k): dict(v) for k, v in (traits or {}).items()}
        self._hover_meta: Dict[str, Dict[str, str]] = {str(k): dict(v) for k, v in (hover_metadata or {}).items()}
        self._window_samples: Optional[Tuple[str, ...]] = None
        self._isolates = ISOLATES_WINDOW
        self._k = pc.DEFAULT_K
        self._seed = pc.DEFAULT_SEED
        self._color_by = "Clade" if "Clade" in self._traits else KMEANS
        self._shape_by = KMEANS
        self._mode = MODE_2D
        self.result: Optional[pc.ClusteringResult] = None
        self._missing = 0
        self._disposed = False
        # Compute: the latest request wins; one runs, one waits.
        self._token = 0
        self._running = False
        self._pending: Optional[Tuple[Any, ...]] = None
        self._shown_params: Optional[Tuple[Any, ...]] = None
        self._compute_started = 0.0
        self._last_seconds = 0.0
        # Drawn points (indices into result.samples) in drawing order.
        self._draw_order = np.zeros(0, dtype=int)
        self._draw_xy = np.zeros((0, 2))
        self._draw_z = np.zeros(0)
        self._hover_idx: Optional[int] = None
        # Isolates with identical profiles: each one's group at its point.
        self._stack_of: Dict[int, np.ndarray] = {}
        self._color_binned = False
        self._external: set = set()
        self._projection = Projection3D()
        self._box = ReferenceBox()
        self._box_corners = np.zeros((0, 2))
        self._rotate_pending = [0.0, 0.0, 0.0]

        self._build_ui()

        self._range_timer = QtCore.QTimer(self)
        self._range_timer.setSingleShot(True)
        self._range_timer.setInterval(_RANGE_DEBOUNCE_MS)
        self._range_timer.timeout.connect(self.request_compute)
        self._frame_timer = QtCore.QTimer(self)
        self._frame_timer.setSingleShot(True)
        self._frame_timer.setInterval(_FRAME_MS)
        self._frame_timer.timeout.connect(self._apply_orbit)

    # ================================================================ UI
    def _build_ui(self) -> None:
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)
        top = QtWidgets.QHBoxLayout()
        top.setContentsMargins(0, 0, 0, 0)
        top.setSpacing(0)
        outer.addLayout(top, stretch=1)

        strip = SettingsStrip(self)
        strip.expandRequested.connect(lambda: self.set_settings_open(True))
        top.addWidget(strip)
        self._settings_strip = strip
        column = self._build_settings_column()
        top.addWidget(column)
        self._settings_column = column

        top.addWidget(self._build_plot_area(), stretch=1)

        self.band = PrevalenceRangeBand(self)
        self.band.set_range(pc.DEFAULT_LOW, pc.DEFAULT_HIGH)
        self.band.rangeEdited.connect(lambda *_a: self._range_timer.start())
        self.band.rangeCommitted.connect(self._on_range_committed)
        outer.addWidget(self.band)

        column.hide()
        self._settings_auto_folded = True
        self._settings_user_opened = False
        self._settings_fit_timer = QtCore.QTimer(self)
        self._settings_fit_timer.setSingleShot(True)
        self._settings_fit_timer.setInterval(0)
        self._settings_fit_timer.timeout.connect(self._fit_settings)

    def _build_settings_column(self) -> SettingsColumn:
        column = SettingsColumn(self, object_name="FilterRail", max_width=300)
        column.collapseRequested.connect(lambda: self.set_settings_open(False))

        card, layout = settings_card(column.content, "Points")
        self._color_combo = settings_combo(card, self._color_options(), self._color_by, self._set_color_by)
        stacked_field(card, layout, "Color by", self._color_combo)
        self._shape_combo = settings_combo(
            card,
            [(KMEANS, KMEANS_LABEL, "One shape per MDS cluster (k-means)"), ("", "None", "Every point a disc")],
            self._shape_by,
            self._set_shape_by,
        )
        stacked_field(card, layout, "Shape by", self._shape_combo)
        column.add_card(card)

        card, layout = settings_card(column.content, "Clustering")
        self._k_spin = QtWidgets.QSpinBox(card)
        self._k_spin.setRange(2, 12)
        self._k_spin.setValue(self._k)
        self._k_spin.setKeyboardTracking(False)
        self._k_spin.setToolTip("Number of MDS clusters: k-means on the three MDS axes (fixed seed)")
        self._k_spin.valueChanged.connect(self._set_k)
        stacked_field(card, layout, "MDS clusters (k)", settings_control(self._k_spin))
        self._isolates_combo = settings_combo(
            card,
            [
                (ISOLATES_WINDOW, "Window samples", "The window's sample set; prevalence over these isolates"),
                (ISOLATES_ALL, "All isolates", "Every isolate of the run, as the script; the window's set stands out"),
            ],
            self._isolates,
            self._set_isolates,
        )
        stacked_field(card, layout, "Isolates", self._isolates_combo)
        column.add_card(card)

        card, layout = settings_card(column.content, "Axes")
        chips = QtWidgets.QWidget(card)
        self._axes_chips = WrapLayout(chips, spacing=6, line_spacing=6)
        self._axes_chips.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(chips)
        self._stats_label = QtWidgets.QLabel("", card)
        self._stats_label.setWordWrap(True)
        self._stats_label.setProperty("role", "muted")
        layout.addWidget(self._stats_label)
        column.add_card(card)

        column.finish(min_width=248)
        column.setFixedWidth(column.maximumWidth())
        return column

    def _build_plot_area(self) -> QtWidgets.QWidget:
        area = QtWidgets.QWidget(self)
        area.setObjectName("PanClusteringPlotArea")
        area.setStyleSheet(f"QWidget#PanClusteringPlotArea {{ background: {COLORS['surface']}; }}")
        layout = QtWidgets.QVBoxLayout(area)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        bar = QtWidgets.QWidget(area)
        bar_layout = QtWidgets.QHBoxLayout(bar)
        bar_layout.setContentsMargins(10, 6, 10, 2)
        bar_layout.setSpacing(0)
        self._mode_buttons: Dict[str, QtWidgets.QToolButton] = {}
        group = QtWidgets.QButtonGroup(bar)
        group.setExclusive(True)
        for mode, text, tip, side in (
            (MODE_2D, "2D", "MDS1 against MDS2", "left"),
            (MODE_3D, "3D", "The first three MDS axes: drag to turn, wheel to zoom, double-click to reset", "right"),
        ):
            button = QtWidgets.QToolButton(bar)
            button.setText(text)
            button.setToolTip(tip)
            button.setCheckable(True)
            button.setChecked(mode == self._mode)
            button.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
            button.setFixedSize(38, 24)
            radius = "border-top-left-radius: 6px; border-bottom-left-radius: 6px;" if side == "left" else (
                "border-top-right-radius: 6px; border-bottom-right-radius: 6px; border-left: none;"
            )
            button.setStyleSheet(
                f"QToolButton {{ background: {COLORS['surface']}; color: {COLORS['text_muted']};"
                f" border: 1px solid {COLORS['border_strong']}; border-radius: 0px; {radius} font-weight: 600; }}"
                f" QToolButton:hover {{ background: {COLORS['surface_hover']}; }}"
                f" QToolButton:checked {{ background: {COLORS['primary_soft']}; color: {COLORS['primary_ink']}; }}"
            )
            button.clicked.connect(lambda _checked=False, m=mode: self.set_mode(m))
            group.addButton(button)
            bar_layout.addWidget(button)
            self._mode_buttons[mode] = button
        bar_layout.addSpacing(12)
        self._status = QtWidgets.QLabel("", bar)
        self._status.setStyleSheet(f"color: {COLORS['text_muted']};")
        bar_layout.addWidget(self._status)
        bar_layout.addSpacing(12)
        self._hint = ElidedHint("Drag to turn · wheel to zoom · double-click to reset", bar)
        self._hint.setToolTip(self._hint.text())
        self._hint.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
        self._hint.setStyleSheet(f"color: {COLORS['text_muted']}; font-size: 11px;")
        self._hint.setVisible(False)
        bar_layout.addWidget(self._hint, stretch=1)
        bar_layout.addSpacing(8)
        self._export_button = QtWidgets.QToolButton(bar)
        self._export_button.setText("Export…")
        self._export_button.setToolTip(
            "Save the cluster assignments, the samples and unique genes of each MDS cluster,"
            " the parameters and the 2D and 3D plots to a folder"
        )
        self._export_button.clicked.connect(self._on_export_clicked)
        bar_layout.addWidget(self._export_button)
        layout.addWidget(bar)

        self._vb = OrbitViewBox()
        self._vb.setMenuEnabled(False)
        self._plot = pg.PlotWidget(viewBox=self._vb)
        # Repaint the whole view on every change: with the default minimal
        # updates, moving the hover highlight left the old one painted.
        self._plot.setViewportUpdateMode(QtWidgets.QGraphicsView.ViewportUpdateMode.FullViewportUpdate)
        self._plot.setBackground(COLORS["surface"])
        item = self._plot.getPlotItem()
        item.setMenuEnabled(False)
        item.hideButtons()
        for name in ("left", "bottom"):
            axis = item.getAxis(name)
            axis.setPen(pg.mkPen(COLORS["border_strong"]))
            axis.setTextPen(pg.mkPen(COLORS["text_muted"]))
            axis.enableAutoSIPrefix(False)
        self._plot.setAspectLocked(True, ratio=1.0)
        layout.addWidget(self._plot, stretch=1)

        # The 3D reference box, behind every point: the far walls' fill and
        # grid, the wireframe with its tick marks, the axis titles.
        pane_fill = QtGui.QColor(COLORS["window"])
        pane_fill.setAlpha(90)
        self._box_panes = QtWidgets.QGraphicsPathItem()
        self._box_panes.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
        self._box_panes.setBrush(QtGui.QBrush(pane_fill))
        self._box_panes.setZValue(-8)
        self._plot.addItem(self._box_panes)
        self._box_grid = pg.PlotDataItem(pen=pg.mkPen(COLORS["border"], width=1.0), connect="pairs")
        self._box_grid.setZValue(-7)
        self._plot.addItem(self._box_grid)
        self._box_edges = pg.PlotDataItem(pen=pg.mkPen(COLORS["border_strong"], width=1.0), connect="pairs")
        self._box_edges.setZValue(-6)
        self._plot.addItem(self._box_edges)
        self._axis_labels: List[pg.TextItem] = []
        for i in range(3):
            label = pg.TextItem(f"MDS{i + 1}", color=COLORS["text_muted"], anchor=(0.5, 0.5))
            label.setZValue(-4)
            self._plot.addItem(label)
            self._axis_labels.append(label)
        self._scatter = DepthScatterItem(pxMode=True)
        self._plot.addItem(self._scatter)
        self._stack_labels: List[pg.TextItem] = []
        self._rings = pg.ScatterPlotItem(pxMode=True)
        self._rings.setZValue(20)
        self._plot.addItem(self._rings)
        card = QtGui.QColor(COLORS["surface"])
        card.setAlpha(235)
        self._hover_card = pg.TextItem(
            color=COLORS["text_strong"],
            anchor=(0.0, 1.4),
            border=pg.mkPen(COLORS["border_strong"]),
            fill=pg.mkBrush(card),
        )
        self._hover_card.setZValue(30)
        self._hover_card.setVisible(False)
        self._plot.addItem(self._hover_card)
        self._key = _GraphKey(self._vb)

        self._message = QtWidgets.QLabel("", self._plot)
        self._message.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self._message.setWordWrap(True)
        self._message.setStyleSheet(f"color: {COLORS['text_muted']}; background: transparent;")
        self._message.setAttribute(QtCore.Qt.WidgetAttribute.WA_TransparentForMouseEvents, True)

        self._vb.orbitDragged.connect(self._on_orbit_dragged)
        # The box's tick marks and title gaps are in pixels: redrawn when the
        # view's range or size changes.
        self._vb.sigRangeChanged.connect(self._on_view_scaled)
        self._vb.sigResized.connect(self._on_view_scaled)
        self._vb.orbitZoomed.connect(self._on_orbit_zoomed)
        self._move_proxy = pg.SignalProxy(self._plot.scene().sigMouseMoved, rateLimit=60, slot=self._on_mouse_moved)
        self._plot.scene().sigMouseClicked.connect(self._on_mouse_clicked)
        self._plot.viewport().installEventFilter(self)
        self._plot.installEventFilter(self)
        self._show_axes_for_mode()
        return area

    # ================================================================ Settings column fold
    def settings_open(self) -> bool:
        return not self._settings_column.isHidden()

    def set_settings_open(self, visible: bool) -> None:
        self._settings_auto_folded = False
        self._settings_user_opened = bool(visible)
        self._show_settings(bool(visible))

    def _show_settings(self, visible: bool) -> None:
        self._settings_column.setVisible(visible)
        self._settings_strip.setVisible(not visible)

    def _fit_settings(self) -> None:
        if not self.isVisible():
            return
        small = self.width() < self._settings_column.width() + self._GRAPH_MIN_WIDTH
        if not small:
            self._settings_user_opened = False
        elif self._settings_user_opened:
            return
        self._settings_auto_folded = auto_fold(
            small, self.settings_open(), self._settings_auto_folded, self._show_settings
        )

    def showEvent(self, event: QtGui.QShowEvent) -> None:  # type: ignore[override]
        super().showEvent(event)
        self._settings_fit_timer.start()

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:  # type: ignore[override]
        super().resizeEvent(event)
        self._settings_fit_timer.start()

    def eventFilter(self, watched: QtCore.QObject, event: QtCore.QEvent) -> bool:
        if watched is self._plot and event.type() == QtCore.QEvent.Type.Resize:
            self._message.setGeometry(self._plot.rect().adjusted(24, 24, -24, -24))
        elif watched is self._plot.viewport() and event.type() == QtCore.QEvent.Type.Leave:
            self._set_hover(None)
        return super().eventFilter(watched, event)

    # ================================================================ options
    def _color_options(self) -> List[Tuple[str, str, str]]:
        options = [(KMEANS, KMEANS_LABEL, "The MDS cluster (k-means) of each isolate")]
        names = [name for name in ("Clade", "Cluster", "ST") if name in self._traits]
        names += sorted((n for n in self._traits if n not in names), key=str.casefold)
        options += [(name, _trait_label(name), f"Color by {_trait_label(name)}") for name in names]
        return options

    def set_traits(self, traits: Mapping[str, Mapping[str, str]], hover_metadata: Mapping[str, Mapping[str, str]]) -> None:
        self._traits = {str(k): dict(v) for k, v in traits.items()}
        self._hover_meta = {str(k): dict(v) for k, v in hover_metadata.items()}
        with QtCore.QSignalBlocker(self._color_combo):
            self._color_combo.clear()
            for value, text, tip in self._color_options():
                self._color_combo.addItem(text, value)
                self._color_combo.setItemData(self._color_combo.count() - 1, tip, QtCore.Qt.ItemDataRole.ToolTipRole)
        if self._color_by != KMEANS and self._color_by not in self._traits:
            self._color_by = KMEANS
        select_combo_data(self._color_combo, self._color_by)
        self._restyle()

    def _set_color_by(self, value: Any) -> None:
        self._color_by = str(value or KMEANS)
        self._restyle()

    def _set_shape_by(self, value: Any) -> None:
        self._shape_by = str(value or "")
        self._restyle()

    def _set_k(self, value: int) -> None:
        self._k = int(value)
        self.request_compute()

    def _set_isolates(self, value: Any) -> None:
        self._isolates = ISOLATES_ALL if value == ISOLATES_ALL else ISOLATES_WINDOW
        self.request_compute()

    def set_mode(self, mode: str) -> None:
        mode = MODE_3D if mode == MODE_3D else MODE_2D
        self._mode_buttons[mode].setChecked(True)
        if mode == self._mode:
            return
        self._mode = mode
        self._set_hover(None)
        self._show_axes_for_mode()
        self._restyle(fit=True)

    def mode(self) -> str:
        return self._mode

    def _show_axes_for_mode(self) -> None:
        item = self._plot.getPlotItem()
        is_3d = self._mode == MODE_3D
        ok = self.result is not None and self.result.ok
        self._vb.set_orbit(is_3d)
        self._hint.setVisible(is_3d and ok)
        for name in ("left", "bottom"):
            # Nothing drawn: no axes, so no stale ticks or labels behind the message.
            item.showAxis(name, not is_3d and ok)
        for box_item in (self._box_panes, self._box_grid, self._box_edges):
            box_item.setVisible(is_3d and ok)
        for label in self._axis_labels:
            label.setVisible(is_3d and ok)

    # ================================================================ window link
    def set_window_samples(self, samples: Sequence[str]) -> None:
        """The window's sample set (recomputes with window samples; with all
        isolates only the out of set styling changes)."""
        self._window_samples = tuple(str(s) for s in samples)
        if self._isolates == ISOLATES_WINDOW:
            self.request_compute()
        else:
            if self.result is None:
                self.request_compute()
            self._restyle()

    def highlight_samples(self, samples: Optional[Sequence[str]]) -> None:
        self._external = {str(s) for s in (samples or ()) if str(s)}
        if self._hover_idx is not None and not self._plot.viewport().underMouse():
            # The pointer left the plot without a Leave event (e.g. over an
            # overlay): its own hover must not stay highlighted.
            self._set_hover(None)
            return
        self._update_rings()

    # ================================================================ compute
    def _params(self) -> Optional[Tuple[Any, ...]]:
        if self._window_samples is None:
            return None
        if self._isolates == ISOLATES_ALL:
            key: Any = ("all",)
            samples: Optional[Tuple[str, ...]] = None
        else:
            samples = tuple(sorted(set(self._window_samples)))
            key = ("set", samples)
        return (key, samples, self.band.low(), self.band.high(), self._k, self._seed)

    def request_compute(self) -> None:
        self._range_timer.stop()
        params = self._params()
        if params is None or self._disposed:
            return
        if self._running:
            # Even when it is the shown one: the running compute is another.
            self._pending = params
            return
        if params == self._shown_params and self.result is not None:
            return
        self._start(params)

    def _start(self, params: Tuple[Any, ...]) -> None:
        self._token += 1
        token = self._token
        self._running = True
        self._pending = None
        self._compute_started = time.perf_counter()
        key, samples, low, high, k, seed = params
        shared = self._shared
        self._status.setText("Computing...")

        def compute() -> Tuple[int, Tuple[Any, ...], Any, Any, int]:
            slice_ = shared.slice_for(key, samples)
            result = shared.cache.compute(slice_, key, low, high, k, seed)
            image = pc.band_image(slice_, result.row_order())
            missing = len(samples) - slice_.n_isolates if samples is not None else 0
            return token, params, result, image, missing

        worker = RefreshWorker(token, compute)
        worker.signals.finished.connect(self._on_computed)
        worker.signals.failed.connect(self._on_failed)
        QtCore.QThreadPool.globalInstance().start(worker)

    def _on_computed(self, payload: object) -> None:
        if self._disposed:
            return
        token, params, result, image, missing = payload  # type: ignore[misc]
        self._running = False
        if token != self._token:
            return
        self._last_seconds = time.perf_counter() - self._compute_started
        self._shown_params = params
        self._missing = int(missing)
        self._set_result(result, image)
        if self._pending is not None and self._pending != params:
            self._start(self._pending)
        else:
            self._pending = None
            self._status.setText("")
        self.computed.emit()

    def _on_failed(self, _token: int, error: str) -> None:
        if self._disposed:
            return
        self._running = False
        self._status.setText("")
        last = error.strip().splitlines()[-1] if error.strip() else "Unknown error"
        self.result = None
        self._shown_params = None
        self._set_hover(None)
        self._find_stacks()
        self._restyle()
        self._show_axes_for_mode()
        self._update_axes_text()
        self._show_message(f"The clustering failed: {last}")
        if self._pending is not None:
            self._start(self._pending)

    def is_computing(self) -> bool:
        return self._running or self._range_timer.isActive()

    def _on_range_committed(self, _low: float, _high: float) -> None:
        self.request_compute()

    # ================================================================ result
    def _set_result(self, result: pc.ClusteringResult, image: np.ndarray) -> None:
        self.result = result
        self.band.set_prevalence(result.slice.prevalence)
        self.band.set_image(image_from_shares(image))
        self.band.set_dropped(result.dropped)
        self._projection.set_points(result.coords if result.ok else np.zeros((0, 3)))
        self._box.fit(self._projection)
        self._set_hover(None)
        self._find_stacks()
        self._update_axes_text()
        self._update_axis_titles()
        self._show_axes_for_mode()
        if not result.ok:
            self._show_message(result.message or "Nothing to draw.")
        else:
            self._show_message("")
        self._restyle(fit=True)

    def _find_stacks(self) -> None:
        """Isolates with identical profiles sit at one point: the points
        shared by several (``_stack_of``: the others at each point)."""
        result = self.result
        self._stack_of = {}
        if result is None or not result.ok:
            return
        _unique, inverse, counts = np.unique(
            np.round(result.coords, 12), axis=0, return_inverse=True, return_counts=True
        )
        inverse = np.asarray(inverse).ravel()
        for group in np.flatnonzero(counts > 1):
            members = np.flatnonzero(inverse == group)
            for idx in members:
                self._stack_of[int(idx)] = members

    def _stacked(self, idx: int) -> np.ndarray:
        """``idx`` and the isolates at the same point, ``idx`` first."""
        members = self._stack_of.get(int(idx))
        if members is None:
            return np.asarray([idx])
        return np.concatenate([[idx], members[members != idx]])

    def _show_message(self, text: str) -> None:
        self._message.setText(text)
        self._message.setVisible(bool(text))
        self._message.setGeometry(self._plot.rect().adjusted(24, 24, -24, -24))
        self._message.raise_()

    def _update_axes_text(self) -> None:
        result = self.result
        while self._axes_chips.count():
            item = self._axes_chips.takeAt(0)
            if item is not None and item.widget() is not None:
                item.widget().deleteLater()
        if result is not None and result.ok:
            tip = (
                "This axis' share of the total inertia: its eigenvalue over the trace of the double-centered"
                " matrix, as R's eig / sum(eig). Jaccard is not Euclidean, so the trace also counts negative"
                " eigenvalues and the shares can be larger than over the positive ones alone."
            )
            for i, value in enumerate(result.explained):
                chip = make_chip(f"MDS{i + 1} {percent_text(value)}", tip)
                chip.setParent(self._axes_chips.parentWidget())
                self._axes_chips.addWidget(chip)
        if result is None:
            self._stats_label.setText("")
            return
        lines = [f"{len(result.samples):,} isolates, {result.n_genes:,} genes, {self._last_seconds:.2f} s"]
        if result.dropped:
            lines.append(f"{len(result.dropped):,} without genes in range")
        if self._missing:
            lines.append(f"{self._missing:,} window samples have no Panaroo profile")
        if self._stack_of:
            points = len({id(m) for m in self._stack_of.values()})
            lines.append(f"{len(self._stack_of):,} isolates share {points:,} points (identical profiles)")
        self._stats_label.setText("\n".join(lines))
        seconds = ", ".join(f"{name} {value * 1000:.0f} ms" for name, value in result.seconds.items())
        self._stats_label.setToolTip(f"Compute: {seconds}" if seconds else "")

    # ================================================================ styling
    def _value_of(self, trait: str, sample: str, index: int) -> str:
        if trait == KMEANS:
            result = self.result
            return str(int(result.labels[index])) if result is not None and result.ok else ""
        return self._traits.get(trait, {}).get(sample, "")

    def _values(self, trait: str) -> List[str]:
        result = self.result
        if result is None or not result.ok:
            return []
        return [self._value_of(trait, s, i) for i, s in enumerate(result.samples)]

    @staticmethod
    def _sorted_values(values: Sequence[str], numeric: bool) -> List[str]:
        unique = {v for v in values if v}
        if numeric:
            return sorted(unique, key=lambda v: (int(v) if v.isdigit() else 10**9, v))
        return sorted(unique, key=str.casefold)

    def _fixed_colors(self, trait: str) -> Dict[str, QtGui.QColor]:
        """The colors the Phylogeny gives the fastGEAR clades and clusters."""
        values = self._traits.get(trait, {})
        if trait == "Clade":
            return dict(_fastgear_lineage_colors(values))
        if trait == "Cluster":
            colors = cluster_colors(list(values.values()), sort_key=_natural_sort_key)
            return {value: QtGui.QColor(color) for value, color in colors.items()}
        return {}

    def _in_set(self) -> np.ndarray:
        result = self.result
        if result is None or not result.ok:
            return np.zeros(0, dtype=bool)
        if self._isolates == ISOLATES_WINDOW or self._window_samples is None:
            return np.ones(len(result.samples), dtype=bool)
        window = set(self._window_samples)
        return np.asarray([s in window for s in result.samples], dtype=bool)

    def _restyle(self, fit: bool = False) -> None:
        """Colors, shapes and the key from the current options."""
        result = self.result
        if result is None or not result.ok:
            self._scatter.clear()
            self._rings.clear()
            self._key.set_sections([])
            self._draw_order = np.zeros(0, dtype=int)
            self._draw_xy = np.zeros((0, 2))
            self._place_stack_labels()
            return
        in_set = self._in_set()
        self._in_set_mask = in_set
        colors = self._values(self._color_by)
        bins = numeric_bins(colors) if self._color_by != KMEANS else None
        self._color_binned = bins is not None
        if bins is not None:
            # A continuous column: quantile bins on a sequential ramp.
            colors, palette_values = bins
            palette = ramp_colors(len(palette_values))
        else:
            palette_values = self._sorted_values([c for c, s in zip(colors, in_set) if s], self._color_by == KMEANS)
            fixed = self._fixed_colors(self._color_by)
            palette = [
                QtGui.QColor(fixed[v]) if v in fixed else _categorical_color(i, len(palette_values))
                for i, v in enumerate(palette_values)
            ]
        self._color_values = palette_values
        missing_color = QtGui.QColor(COLORS["tree_inactive"])
        out_color = out_of_set_color()
        palette += [missing_color, out_color]
        index_of = {v: i for i, v in enumerate(palette_values)}
        missing_idx, out_idx = len(palette_values), len(palette_values) + 1
        color_idx = np.asarray(
            [out_idx if not s else index_of.get(c, missing_idx) for c, s in zip(colors, in_set)], dtype=int
        )
        self._color_idx = color_idx
        self._palette = palette
        background = QtGui.QColor(COLORS["surface"])
        bg = (background.red(), background.green(), background.blue())
        # Flat per color for 2D; shaded spheres for 3D, per color and depth level.
        table = np.empty(len(palette), dtype=object)
        shaded = np.empty((len(palette), len(DEPTH_SIZES)), dtype=object)
        for ci, color in enumerate(palette):
            table[ci] = QtGui.QBrush(color)
            rgba = color.getRgb()
            for level, fade in enumerate(DEPTH_FADES):
                faded = QtGui.QColor(*fade_toward(rgba, bg, fade))
                shaded[ci, level] = sphere_brush(faded)
        self._brush_table = table
        self._shaded_table = shaded
        shapes = self._values(KMEANS) if self._shape_by == KMEANS else [""] * len(result.samples)
        shape_values = self._sorted_values(shapes, True)
        symbol_of = {v: _FORCE_SYMBOLS[i % len(_FORCE_SYMBOLS)] for i, v in enumerate(shape_values)}
        self._symbol_values = shape_values
        self._symbol_of = symbol_of
        self._symbols = np.asarray([symbol_of.get(s, "o") for s in shapes], dtype=object)
        # Out of set points first, so the set is drawn over them.
        self._base_order = np.concatenate([np.flatnonzero(~in_set), np.flatnonzero(in_set)])
        outline = pg.mkPen(QtGui.QColor(31, 44, 59, 150), width=0.8)
        faint = pg.mkPen(out_of_set_color(QtGui.QColor(31, 44, 59, 150)), width=0.8)
        self._pens = np.asarray([outline if s else faint for s in in_set], dtype=object)
        self._key_sections_update(colors, shapes, in_set)
        self._scatter.drop_level_styles()
        self._redraw(fit=fit)

    def _key_sections_update(self, colors: Sequence[str], shapes: Sequence[str], in_set: np.ndarray) -> None:
        sections: List[Tuple[str, List[GraphKeyEntry]]] = []
        color_counts: Dict[str, int] = {}
        missing = 0
        for value, keep in zip(colors, in_set):
            if not keep:
                continue
            if value:
                color_counts[value] = color_counts.get(value, 0) + 1
            else:
                missing += 1
        color_label = KMEANS_LABEL if self._color_by == KMEANS else _trait_label(self._color_by)
        if self._color_binned:
            color_label += " (quantiles)"
        entries = [
            GraphKeyEntry(value, color_counts.get(value, 0), color=self._palette[i])
            for i, value in enumerate(self._color_values)
        ]
        if missing:
            entries.append(GraphKeyEntry("No value", missing, color=QtGui.QColor(COLORS["tree_inactive"])))
        sections.append((f"Color: {color_label}", entries))
        if self._shape_by == KMEANS and self._color_by != KMEANS and self._mode == MODE_2D:
            shape_counts: Dict[str, int] = {}
            for value, keep in zip(shapes, in_set):
                if keep and value:
                    shape_counts[value] = shape_counts.get(value, 0) + 1
            sections.append((
                f"Shape: {KMEANS_LABEL}",
                [GraphKeyEntry(v, shape_counts.get(v, 0), symbol=self._symbol_of[v]) for v in self._symbol_values],
            ))
        self._key.set_sections(sections)

    # ================================================================ drawing
    def _redraw(self, fit: bool = False) -> None:
        result = self.result
        if result is None or not result.ok or not hasattr(self, "_brush_table"):
            return
        if self._mode == MODE_2D:
            order = self._base_order
            xy = result.coords[order, :2]
            z = np.zeros(order.size)
        else:
            all_xy, all_z = self._projection.project()
            # Back to front; out of set points under the set at any depth.
            depth_order = np.argsort(all_z, kind="stable")
            in_set = self._in_set_mask
            order = np.concatenate([depth_order[~in_set[depth_order]], depth_order[in_set[depth_order]]])
            xy, z = all_xy[order], all_z[order]
            levels = depth_levels(z)
            self._draw_box()
        self._draw_order, self._draw_xy, self._draw_z = order, xy, z
        if self._mode == MODE_2D:
            self._scatter.drop_level_styles()
            self._scatter.setData(
                x=xy[:, 0],
                y=xy[:, 1],
                brush=list(self._brush_table[self._color_idx[order]]),
                symbol=list(self._symbols[order]),
                size=_POINT_PX,
                pen=list(self._pens[order]),
            )
        else:
            if not self._scatter.has_level_styles():
                n_levels = len(DEPTH_SIZES)
                self._scatter.set_level_styles(
                    self._shaded_table[self._color_idx[None, :], np.arange(n_levels)[:, None]],
                    np.repeat(_SPHERE_PX * np.asarray(DEPTH_SIZES)[:, None], self._color_idx.size, axis=1),
                    ["o"] * self._color_idx.size,
                    list(self._pens),
                )
            self._scatter.show_frame(xy, order, levels)
        if fit:
            self._fit_view()
        self._key.set_nodes(xy, 0.5 * (_SPHERE_PX if self._mode == MODE_3D else _POINT_PX))
        if fit:
            self._key.place()
        self._place_stack_labels()
        self._update_rings()
        if self._hover_idx is not None:
            self._place_hover_card()

    def _place_stack_labels(self) -> None:
        """A count beside each point that several isolates share."""
        groups: Dict[int, np.ndarray] = {}
        for members in self._stack_of.values():
            groups.setdefault(id(members), members)
        stacks = sorted(groups.values(), key=len, reverse=True)[:_STACK_LABELS_MAX]
        while len(self._stack_labels) < len(stacks):
            label = pg.TextItem("", color=COLORS["text_muted"], anchor=(-0.15, 0.85))
            label.setZValue(15)
            self._plot.addItem(label)
            self._stack_labels.append(label)
        position = {int(p): i for i, p in enumerate(self._draw_order)}
        for label, members in zip(self._stack_labels, stacks):
            row = position.get(int(members[0]))
            if row is None:
                label.setVisible(False)
                continue
            label.setText(f"×{len(members)}")
            label.setPos(float(self._draw_xy[row, 0]), float(self._draw_xy[row, 1]))
            label.setVisible(True)
        for label in self._stack_labels[len(stacks):]:
            label.setVisible(False)

    def _on_view_scaled(self, *_args: object) -> None:
        result = self.result
        if self._mode == MODE_3D and result is not None and result.ok and hasattr(self, "_brush_table"):
            self._draw_box()

    def _box_px(self) -> float:
        """View units per pixel (the 3D view is aspect locked)."""
        px = self._vb.viewPixelSize()[0]
        return float(px) if px and np.isfinite(px) else 2.0 * VIEW_EXTENT / 800.0

    def _draw_box(self) -> None:
        """The reference box as the cloud is turned now."""
        px = self._box_px()
        frame = self._box.frame(self._projection, _BOX_TICK_PX * px, _BOX_TITLE_GAP_PX * px)
        path = QtGui.QPainterPath()
        for pane in frame.panes:
            path.addPolygon(QtGui.QPolygonF([QtCore.QPointF(float(x), float(y)) for x, y in pane]))
            path.closeSubpath()
        self._box_panes.setPath(path)
        self._box_grid.setData(frame.grid[:, 0], frame.grid[:, 1])
        self._box_edges.setData(frame.edges[:, 0], frame.edges[:, 1])
        (x0, x1), (y0, y1) = self._vb.viewRange()
        for label, pos, anchor in zip(self._axis_labels, frame.labels, frame.anchors):
            if tuple(label.anchor) != anchor:
                label.setAnchor(anchor)
            # Kept in the view (a title beside an edge seen end on sits out
            # past the box's corners).
            rect = label.boundingRect()
            w, h = rect.width() * px, rect.height() * px
            ax, ay = anchor
            x = min(max(float(pos[0]), x0 + ax * w), max(x0 + ax * w, x1 - (1.0 - ax) * w))
            y = min(max(float(pos[1]), y0 + (1.0 - ay) * h), max(y0 + (1.0 - ay) * h, y1 - ay * h))
            label.setPos(x, y)
        self._box_corners = frame.corners

    def _content_bounds(self) -> Tuple[np.ndarray, np.ndarray]:
        """The drawn points (and in 3D the box and its titles), low and high corners."""
        xy = self._draw_xy if self._draw_xy.size else np.zeros((1, 2))
        if self._mode == MODE_3D:
            px = self._box_px()
            parts = [xy, self._box_corners]
            for label in self._axis_labels:
                rect = label.boundingRect()
                w, h = rect.width() * px, rect.height() * px
                ax, ay = label.anchor
                x0, top = label.pos().x() - ax * w, label.pos().y() + ay * h
                parts.append(np.array([[x0, top - h], [x0 + w, top]]))
            xy = np.vstack(parts)
        return xy.min(axis=0), xy.max(axis=0)

    def _view_extent_3d(self) -> float:
        """Half the fixed 3D view: the box at any turn (a sphere of the box's
        radius r is drawn at most r D / sqrt(D^2 - r^2) from the center)."""
        r = min(self._box.radius(), 0.9 * CAMERA_DISTANCE)
        if r <= 0:
            return VIEW_EXTENT
        return 1.04 * r * CAMERA_DISTANCE / float(np.sqrt(CAMERA_DISTANCE**2 - r * r))

    def _fit_view(self, tight: bool = False) -> None:
        """2D: the points; 3D: the fixed view, or the cloud as it is turned
        (``tight``, for the export)."""
        vb = self._vb
        if self._mode == MODE_3D:
            vb.enableAutoRange(enable=False)
            if tight:
                lo, hi = self._content_bounds()
                pad = 0.04 * float(np.max(np.maximum(hi - lo, 1e-6)))
                lo, hi = lo - pad, hi + pad
                # The box reaches the corners: room for the open key beside it.
                area = vb.boundingRect()
                key = self._key.reserved_size(area.width(), area.height())
                if key is not None and area.width() > key.width():
                    span = np.maximum(hi - lo, 1e-6)
                    scale = min((area.width() - key.width()) / span[0], area.height() / span[1])
                    hi[0] += key.width() / scale
                vb.setRange(xRange=(lo[0], hi[0]), yRange=(lo[1], hi[1]), padding=0.0)
            else:
                extent = self._view_extent_3d()
                vb.setRange(xRange=(-extent, extent), yRange=(-extent, extent), padding=0.0)
            return
        result = self.result
        xy = result.coords[:, :2] if result is not None and result.ok else np.zeros((1, 2))
        lo, hi = xy.min(axis=0), xy.max(axis=0)
        span = np.maximum(hi - lo, 1e-6)
        pad = 0.08 * span.max()
        vb.setRange(
            xRange=(lo[0] - pad, hi[0] + pad), yRange=(lo[1] - pad, hi[1] + pad), padding=0.0
        )
        self._update_axis_titles()

    def _update_axis_titles(self) -> None:
        item = self._plot.getPlotItem()
        result = self.result
        if result is not None and result.ok:
            item.setLabel("bottom", f"MDS1 ({percent_text(result.explained[0])})", color=COLORS["text_muted"])
            item.setLabel("left", f"MDS2 ({percent_text(result.explained[1])})", color=COLORS["text_muted"])
            for i, label in enumerate(self._axis_labels):
                label.setText(f"MDS{i + 1} ({percent_text(result.explained[i])})", color=COLORS["text_muted"])
        else:
            item.setLabel("bottom", "")
            item.setLabel("left", "")

    # ---------------------------------------------------------------- 3D orbit
    def _on_orbit_dragged(self, dx: float, dy: float) -> None:
        self._rotate_pending[0] += dx
        self._rotate_pending[1] += dy
        if not self._frame_timer.isActive():
            self._frame_timer.start()

    def _on_orbit_zoomed(self, steps: float) -> None:
        self._rotate_pending[2] += steps
        if not self._frame_timer.isActive():
            self._frame_timer.start()

    def _apply_orbit(self) -> None:
        dx, dy, steps = self._rotate_pending
        self._rotate_pending = [0.0, 0.0, 0.0]
        if dx or dy:
            self._projection.rotate(dx, dy)
        if steps:
            self._projection.zoom_by(steps)
        self._redraw()

    def reset_view(self) -> None:
        if self._mode == MODE_3D:
            self._projection.reset()
        self._redraw(fit=True)

    # ================================================================ hover
    def _point_at(self, scene_pos: QtCore.QPointF) -> Optional[int]:
        """The drawn point under ``scene_pos`` (index into result.samples):
        the nearest within reach, in 3D the frontmost."""
        if self._draw_xy.size == 0 or not self._vb.sceneBoundingRect().contains(scene_pos):
            return None
        local = self._vb.mapSceneToView(scene_pos)
        px = self._vb.viewPixelSize()
        if not px[0] or not px[1]:
            return None
        dx = (self._draw_xy[:, 0] - local.x()) / px[0]
        dy = (self._draw_xy[:, 1] - local.y()) / px[1]
        dist = np.hypot(dx, dy)
        reach = max(_HOVER_PX, 0.5 * _SPHERE_PX * max(DEPTH_SIZES) + 1.0) if self._mode == MODE_3D else _HOVER_PX
        near = np.flatnonzero(dist <= reach)
        if near.size == 0:
            return None
        if self._mode == MODE_3D:
            pick = near[int(np.argmax(self._draw_z[near]))]
        else:
            # The nearest; among equals the one drawn last (on top).
            pick = near[np.lexsort((-near, dist[near]))[0]]
        return int(self._draw_order[pick])

    def _on_mouse_moved(self, event: Tuple[QtCore.QPointF]) -> None:
        if QtWidgets.QApplication.mouseButtons() != QtCore.Qt.MouseButton.NoButton:
            return
        self._set_hover(self._point_at(event[0]))

    def _set_hover(self, idx: Optional[int]) -> None:
        if idx == self._hover_idx:
            return
        self._hover_idx = idx
        result = self.result
        if idx is None or result is None or not result.ok:
            self._hover_idx = None
            self._hover_card.setVisible(False)
            self._key.set_faded(None)
            self.samplesHovered.emit([])
        else:
            self._hover_card.setHtml(self._hover_html(idx))
            self._hover_card.setVisible(True)
            self._place_hover_card()
            self.samplesHovered.emit([result.samples[int(i)] for i in self._stacked(idx)])
        self._update_rings()

    def hovered_sample(self) -> Optional[str]:
        if self._hover_idx is None or self.result is None:
            return None
        return self.result.samples[self._hover_idx]

    def _hover_html(self, idx: int) -> str:
        result = self.result
        assert result is not None
        sample = result.samples[idx]
        muted = COLORS["text_muted"]
        lines = [f"<b>{html.escape(sample)}</b>"]
        shown = set()
        if self._color_by != KMEANS:
            value = self._value_of(self._color_by, sample, idx) or "no value"
            lines.append(f"{html.escape(_trait_label(self._color_by))}: {html.escape(value)}")
            shown.add(self._color_by.casefold())
        lines.append(f"{KMEANS_LABEL}: {int(result.labels[idx])}")
        for name in ("Clade", "Cluster"):
            value = self._traits.get(name, {}).get(sample, "")
            if value and name.casefold() not in shown:
                lines.append(f"<span style='color:{muted}'>{_trait_label(name)}: {html.escape(value)}</span>")
                shown.add(name.casefold())
        for label, value in self._hover_meta.get(sample, {}).items():
            if label.casefold() not in shown:
                lines.append(f"<span style='color:{muted}'>{html.escape(label)}: {html.escape(value)}</span>")
        if self._isolates == ISOLATES_ALL and not self._in_set_mask[idx]:
            lines.append(f"<span style='color:{muted}'>Outside the window's sample set</span>")
        others = self._stacked(idx)[1:]
        if others.size:
            shown_names = [html.escape(result.samples[int(i)]) for i in others[:8]]
            more = f"<br>and {others.size - len(shown_names):,} more" if others.size > len(shown_names) else ""
            lines.append(
                f"<b>{others.size + 1} isolates at this point</b> (identical profiles in this range):"
                f"<br>{'<br>'.join(shown_names)}{more}"
            )
        return "<br>".join(lines)

    def _xy_of(self, idx: int) -> Optional[Tuple[float, float]]:
        where = np.flatnonzero(self._draw_order == idx)
        if where.size == 0:
            return None
        x, y = self._draw_xy[int(where[0])]
        return float(x), float(y)

    def _place_hover_card(self) -> None:
        """Above and right of its point, flipped left or below when that side
        would leave the plot."""
        if self._hover_idx is None:
            return
        coords = self._xy_of(self._hover_idx)
        if coords is None:
            return
        card = self._hover_card
        x, y = coords
        node = self._vb.mapViewToScene(QtCore.QPointF(x, y))
        area = self._vb.sceneBoundingRect()
        size = card.boundingRect()
        w, h = float(size.width()), float(size.height())
        right_room = area.right() - node.x()
        left_room = node.x() - area.left()
        ax = 0.0 if right_room >= w or right_room >= left_room else 1.0
        above_room = node.y() - area.top()
        below_room = area.bottom() - node.y()
        ay = 1.4 if above_room >= 1.4 * h or above_room >= below_room else -0.4
        if tuple(card.anchor) != (ax, ay):
            card.setAnchor((ax, ay))
        card.setPos(x, y)
        self._key.set_faded(QtCore.QRectF(node.x() - ax * w, node.y() - ay * h, w, h))

    def _update_rings(self) -> None:
        result = self.result
        if result is None or not result.ok or self._draw_order.size == 0:
            self._rings.clear()
            return
        wanted = set()
        if self._external:
            wanted = {i for i, s in enumerate(result.samples) if s in self._external}
        if self._hover_idx is not None:
            wanted.add(self._hover_idx)
        if not wanted:
            self._rings.clear()
            self._plot.viewport().update()
            return
        mask = np.isin(self._draw_order, np.fromiter(wanted, dtype=int))
        xy = self._draw_xy[mask]
        color = QtGui.QColor(_HIGHLIGHT_COLOR)
        if self._mode == MODE_3D:
            symbols: Any = "o"
            # Each highlight as big as its sphere at its depth, plus a little.
            size: Any = _SPHERE_PX * np.asarray(DEPTH_SIZES)[depth_levels(self._draw_z[mask])] + 3.0
            brush = sphere_brush(color)
        else:
            symbols = list(self._symbols[self._draw_order[mask]])
            size = _POINT_PX + 3.0
            brush = QtGui.QBrush(color)
        self._rings.setData(
            x=xy[:, 0],
            y=xy[:, 1],
            symbol=symbols,
            size=size,
            brush=brush,
            pen=pg.mkPen(QtGui.QColor(31, 44, 59, 200), width=0.9),
        )
        self._plot.viewport().update()

    # ================================================================ clicks
    def _on_mouse_clicked(self, event: Any) -> None:
        idx = self._point_at(event.scenePos())
        result = self.result
        if event.button() == QtCore.Qt.MouseButton.LeftButton and event.double():
            if idx is not None and result is not None:
                self.sampleOpenRequested.emit(result.samples[idx])
            else:
                self.reset_view()
            event.accept()
        elif event.button() == QtCore.Qt.MouseButton.RightButton and idx is not None and result is not None:
            event.accept()
            self._cluster_menu(int(result.labels[idx]), event.screenPos())

    def cluster_payload(self, cluster: int) -> Optional[Dict[str, Any]]:
        result = self.result
        if result is None or not result.ok:
            return None
        members = [s for s, label in zip(result.samples, result.labels) if int(label) == int(cluster)]
        if not members:
            return None
        return {"name": str(cluster), "samples": members, "label": f"{KMEANS_LABEL} {cluster}"}

    def _cluster_menu(self, cluster: int, screen_pos: QtCore.QPointF) -> None:
        payload = self.cluster_payload(cluster)
        if payload is None:
            return
        menu = QtWidgets.QMenu(self)
        count = len(payload["samples"])
        only = menu.addAction(f"Show only {KMEANS_LABEL} {cluster} ({count}) in this window")
        new = menu.addAction(f"Open {KMEANS_LABEL} {cluster} in a new window")
        chosen = menu.exec(QtCore.QPoint(int(screen_pos.x()), int(screen_pos.y())))
        if chosen is only:
            self.subsetRequested.emit(payload, False)
        elif chosen is new:
            self.subsetRequested.emit(payload, True)

    # ================================================================ export
    def _on_export_clicked(self) -> None:
        if self.result is None or not self.result.ok:
            QtWidgets.QMessageBox.information(self, "Export", "There is no clustering to export yet.")
            return
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "Export the clustering to a folder")
        if not folder:
            return
        try:
            paths = self.export_to(Path(folder))
        except Exception as exc:  # disk full, no permission
            QtWidgets.QMessageBox.warning(self, "Export", f"The export failed: {exc}")
            return
        QtWidgets.QMessageBox.information(
            self, "Export", f"Saved {len(paths)} files to {folder}:\n" + "\n".join(p.name for p in paths)
        )

    def export_prefix(self) -> str:
        result = self.result
        assert result is not None

        def pct(v: float) -> str:
            return percent_text(v).rstrip("%")

        return f"pan_mds_{pct(result.low)}-{pct(result.high)}_k{result.k}_"

    def export_to(self, folder: Path) -> List[Path]:
        """Write the tables, the parameters and the plots; returns the paths."""
        result = self.result
        if result is None or not result.ok:
            return []
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        prefix = self.export_prefix()
        color_label = "" if self._color_by == KMEANS else self._color_by
        clusters, members, unique = pc.export_tables(
            result, color_label, self._traits.get(self._color_by, {}) if color_label else None
        )
        paths: List[Path] = []
        for name, frame in (("clusters.csv", clusters), ("cluster_samples.csv", members), ("unique_genes.csv", unique)):
            path = folder / f"{prefix}{name}"
            frame.write_csv(path)
            paths.append(path)
        path = folder / f"{prefix}parameters.txt"
        path.write_text(self._parameters_text(), encoding="utf-8")
        paths.append(path)
        paths += self._export_plots(folder, prefix)
        return paths

    def _parameters_text(self) -> str:
        result = self.result
        assert result is not None
        lines = [
            "BACTolomew pangenome clustering",
            f"Isolates: {'all isolates of the run' if self._isolates == ISOLATES_ALL else 'the window sample set'}",
            f"Prevalence range: more than {percent_text(result.low)} and less than {percent_text(result.high)}"
            + (" (100% keeps genes present in every isolate)" if result.high >= 1.0 else ""),
            f"Analyzed isolates: {result.slice.n_isolates}",
            f"Genes in range: {result.n_genes} of {len(result.slice.genes)}",
            "Distance: Jaccard on binary presence/absence",
            "Ordination: classical MDS (Torgerson, 3 axes)",
            "Explained (eigenvalue / trace): "
            + ", ".join(f"MDS{i + 1} {percent_text(v)}" for i, v in enumerate(result.explained)),
            f"MDS clusters (Cluster_Number): k-means on the 3 MDS axes, k = {result.k}, seed = {result.seed},"
            f" k-means++ starts, best of {pc.DEFAULT_N_INIT}",
            f"Clustered isolates: {len(result.samples)}",
            f"Isolates with no genes in range (not clustered): {len(result.dropped)}",
        ]
        lines += [f"  {name}" for name in result.dropped]
        lines.append(f"Color: {KMEANS_LABEL if self._color_by == KMEANS else _trait_label(self._color_by)}")
        return "\n".join(lines) + "\n"

    def _export_size(self) -> Tuple[int, int]:
        """A fixed width, and the height the drawn content needs (the plot is
        aspect locked), within bounds: the same files whatever the dock."""
        lo, hi = self._content_bounds()
        span = np.maximum(hi - lo, 1e-6)
        width = float(_EXPORT_WIDTH)
        if self._mode == MODE_3D:
            # The key goes beside the box (see ``_fit_view``).
            key = self._key.reserved_size(width, float(_EXPORT_HEIGHT[1]))
            if key is not None:
                width = max(0.5 * width, width - key.width())
        height = int(np.clip(width * span[1] / span[0], *_EXPORT_HEIGHT))
        return _EXPORT_WIDTH, height + (60 if self._mode == MODE_2D else 0)

    def _export_plots(self, folder: Path, prefix: str) -> List[Path]:
        """The 2D (PNG and SVG) and 3D (PNG) plots, drawn off screen at a
        fixed size; the 3D one as it is turned now, fitted to the cloud. The
        on-screen view (mode, zoom, pan) is left as it was."""
        from pyqtgraph import exporters

        paths: List[Path] = []
        mode = self._mode
        x_range, y_range = self._vb.viewRange()
        self._set_hover(None)
        area = self._plot.parentWidget()
        layout = area.layout()
        index = layout.indexOf(self._plot)
        host = QtWidgets.QWidget()
        host.setAttribute(QtCore.Qt.WidgetAttribute.WA_DontShowOnScreen, True)
        self._plot.setParent(host)
        self._plot.show()
        host.show()
        try:
            for target in (MODE_2D, MODE_3D):
                self._mode = target
                self._show_axes_for_mode()
                self._redraw()
                width, height = self._export_size()
                host.resize(width, height)
                self._plot.setGeometry(0, 0, width, height)
                self._plot.getPlotItem().layout.activate()
                self._redraw()
                self._fit_view(tight=True)
                self._key.place()
                image = exporters.ImageExporter(self._plot.getPlotItem())
                image.parameters()["width"] = width
                path = folder / f"{prefix}plot_{target}.png"
                image.export(str(path))
                paths.append(path)
                if target == MODE_2D:
                    path = folder / f"{prefix}plot_2d.svg"
                    self._render_svg(path)
                    paths.append(path)
        finally:
            self._plot.setParent(area)
            layout.insertWidget(index, self._plot, stretch=1)
            self._plot.show()
            host.deleteLater()
            self._mode = mode
            self._mode_buttons[mode].setChecked(True)
            self._show_axes_for_mode()
            self._redraw()
            self._vb.setRange(xRange=x_range, yRange=y_range, padding=0.0)
            self._key.place()
        return paths

    def _render_svg(self, path: Path) -> None:
        """The plot as it is drawn, as vectors (pyqtgraph's SVG exporter
        fails on this Qt's path strings)."""
        from PyQt6.QtSvg import QSvgGenerator

        size = self._plot.size()
        generator = QSvgGenerator()
        generator.setFileName(str(path))
        generator.setSize(size)
        generator.setViewBox(QtCore.QRect(QtCore.QPoint(0, 0), size))
        generator.setTitle("Pangenome clustering")
        painter = QtGui.QPainter(generator)
        # Points as shapes, not as the cached bitmaps drawn on screen.
        items = (self._scatter, self._rings)
        for item in items:
            item.setExportMode(True, {"antialias": True, "painter": painter})
        try:
            self._plot.render(painter)
        finally:
            painter.end()
            for item in items:
                item.setExportMode(False)

    # ================================================================ state
    def save_state(self) -> Dict[str, Any]:
        state = {
            "low": self.band.low(),
            "high": self.band.high(),
            "k": self._k,
            "isolates": self._isolates,
            "color_by": self._color_by,
            "shape_by": self._shape_by,
            "mode": self._mode,
            "view3d": self._projection.state(),
        }
        if not self._settings_auto_folded:
            state["settings_open"] = self.settings_open()
        return state

    def restore_state(self, state: Mapping[str, Any]) -> None:
        try:
            low, high = float(state.get("low", self.band.low())), float(state.get("high", self.band.high()))
            if 0.0 <= low < high <= 1.0:
                self.band.set_range(low, high)
        except (TypeError, ValueError):
            pass
        try:
            k = int(state.get("k", self._k))
            self._k = max(2, min(12, k))
        except (TypeError, ValueError):
            pass
        with QtCore.QSignalBlocker(self._k_spin):
            self._k_spin.setValue(self._k)
        if state.get("isolates") in (ISOLATES_WINDOW, ISOLATES_ALL):
            self._isolates = str(state["isolates"])
            select_combo_data(self._isolates_combo, self._isolates)
        color = str(state.get("color_by") or "")
        if color == KMEANS or color in self._traits:
            self._color_by = color
            select_combo_data(self._color_combo, color)
        if state.get("shape_by") in (KMEANS, ""):
            self._shape_by = str(state["shape_by"])
            select_combo_data(self._shape_combo, self._shape_by)
        self._projection.set_state(state.get("view3d"))
        if "settings_open" in state:
            self.set_settings_open(bool(state["settings_open"]))
        self.set_mode(str(state.get("mode") or MODE_2D))
        self._restyle()
        self.request_compute()

    def dispose(self) -> None:
        self._disposed = True
        self._token += 1
        for timer in (self._range_timer, self._frame_timer, self._settings_fit_timer):
            timer.stop()
        self._move_proxy.disconnect()
