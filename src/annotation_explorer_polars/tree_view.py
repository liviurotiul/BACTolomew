from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np
from PyQt6 import QtCore, QtGui, QtWidgets


@dataclass
class _TreeNode:
    name: str = ""
    children: List["_TreeNode"] = field(default_factory=list)
    branch_length: Optional[float] = None

    @property
    def is_leaf(self) -> bool:
        return len(self.children) == 0


@dataclass
class _LayoutNode:
    node: _TreeNode
    parent: Optional["_LayoutNode"] = None
    depth: int = 0
    number: int = 1
    children: List["_LayoutNode"] = field(default_factory=list)
    prelim: float = 0.0
    mod: float = 0.0
    shift: float = 0.0
    change: float = 0.0
    ancestor: Optional["_LayoutNode"] = None
    thread: Optional["_LayoutNode"] = None
    left_sibling: Optional["_LayoutNode"] = None
    leftmost_sibling: Optional["_LayoutNode"] = None
    x: float = 0.0
    y: float = 0.0

    @property
    def is_leaf(self) -> bool:
        return len(self.children) == 0

    def left(self) -> Optional["_LayoutNode"]:
        if self.children:
            return self.children[0]
        return self.thread

    def right(self) -> Optional["_LayoutNode"]:
        if self.children:
            return self.children[-1]
        return self.thread


@dataclass
class _GroupSummary:
    counts: Dict[str, int]
    total: int
    dominant: Optional[str]
    purity: float
    leaf_order: List[str]
    group_order: List[str]
    structural_score: float


class _NewickParser:
    def __init__(self, text: str) -> None:
        self.text = text.strip()
        self.n = len(self.text)
        self.i = 0

    def parse(self) -> _TreeNode:
        if not self.text:
            raise ValueError("Empty Newick string.")
        node = self._parse_subtree()
        self._skip_ws()
        if self._peek() == ";":
            self.i += 1
        return node

    def _peek(self) -> str:
        if self.i >= self.n:
            return ""
        return self.text[self.i]

    def _skip_ws(self) -> None:
        while self.i < self.n and self.text[self.i].isspace():
            self.i += 1

    def _parse_subtree(self) -> _TreeNode:
        self._skip_ws()
        if self._peek() == "(":
            self.i += 1
            children: List[_TreeNode] = []
            while True:
                children.append(self._parse_subtree())
                self._skip_ws()
                ch = self._peek()
                if ch == ",":
                    self.i += 1
                    continue
                if ch == ")":
                    self.i += 1
                    break
                raise ValueError(f"Invalid Newick near index {self.i}.")
            name = self._parse_label()
            branch = self._parse_branch_length()
            return _TreeNode(name=name, children=children, branch_length=branch)

        name = self._parse_label()
        branch = self._parse_branch_length()
        return _TreeNode(name=name, branch_length=branch)

    def _parse_label(self) -> str:
        self._skip_ws()
        if self.i >= self.n:
            return ""
        if self._peek() == "'":
            self.i += 1
            parts: List[str] = []
            while self.i < self.n:
                ch = self.text[self.i]
                if ch == "'":
                    self.i += 1
                    if self.i < self.n and self.text[self.i] == "'":
                        parts.append("'")
                        self.i += 1
                        continue
                    break
                parts.append(ch)
                self.i += 1
            return "".join(parts).strip()

        start = self.i
        while self.i < self.n and self.text[self.i] not in ",():;":
            self.i += 1
        return self.text[start:self.i].strip()

    def _parse_branch_length(self) -> Optional[float]:
        self._skip_ws()
        if self._peek() != ":":
            return None
        self.i += 1
        start = self.i
        while self.i < self.n and self.text[self.i] not in ",() ;":
            self.i += 1
        raw = self.text[start:self.i].strip()
        if not raw:
            return None
        try:
            return float(raw)
        except Exception:
            return None


class _LeafLabelItem(QtWidgets.QGraphicsSimpleTextItem):
    def __init__(self, label_text: str, leaf_key: str, owner: "RaxmlTreeView") -> None:
        super().__init__(label_text)
        self._leaf_name = leaf_key
        self._owner = owner
        hover_enabled = owner.hover_interactions_enabled()
        self.setAcceptHoverEvents(hover_enabled)
        if hover_enabled:
            self.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)

    def hoverEnterEvent(self, event: Optional[QtWidgets.QGraphicsSceneHoverEvent]) -> None:
        self._owner._on_leaf_hover_enter(self._leaf_name)
        super().hoverEnterEvent(event)

    def hoverLeaveEvent(self, event: Optional[QtWidgets.QGraphicsSceneHoverEvent]) -> None:
        self._owner._on_leaf_hover_leave(self._leaf_name)
        super().hoverLeaveEvent(event)

    def mouseDoubleClickEvent(self, event: Optional[QtWidgets.QGraphicsSceneMouseEvent]) -> None:
        self._owner._on_leaf_double_click(self._leaf_name)
        if event is not None:
            event.accept()
        super().mouseDoubleClickEvent(event)


class _ActionTextItem(QtWidgets.QGraphicsSimpleTextItem):
    def __init__(
        self,
        text: str,
        callback: Callable[[], None],
        *,
        enable_hover: bool = True,
        hover_enter_callback: Optional[Callable[[], None]] = None,
        hover_leave_callback: Optional[Callable[[], None]] = None,
        parent: Optional[QtWidgets.QGraphicsItem] = None,
    ) -> None:
        super().__init__(text, parent)
        self._callback = callback
        self._hover_enter_callback = hover_enter_callback
        self._hover_leave_callback = hover_leave_callback
        self._hover_scale = 1.05
        self._base_font: Optional[QtGui.QFont] = None
        self._base_brush: Optional[QtGui.QBrush] = None
        self.setAcceptHoverEvents(bool(enable_hover))
        if enable_hover:
            self.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)

    def _capture_base_style(self) -> None:
        if self._base_font is None:
            self._base_font = QtGui.QFont(self.font())
        if self._base_brush is None:
            self._base_brush = QtGui.QBrush(self.brush())

    def hoverEnterEvent(self, event: Optional[QtWidgets.QGraphicsSceneHoverEvent]) -> None:
        self._capture_base_style()
        if self._hover_enter_callback is not None:
            self._hover_enter_callback()
        hover_font = QtGui.QFont(self._base_font) if self._base_font is not None else QtGui.QFont(self.font())
        self.setFont(hover_font)
        hover_brush = QtGui.QBrush(QtGui.QColor("#ea580c"))
        self.setBrush(hover_brush)
        self.setOpacity(1.0)
        self.setScale(self._hover_scale)
        super().hoverEnterEvent(event)

    def hoverLeaveEvent(self, event: Optional[QtWidgets.QGraphicsSceneHoverEvent]) -> None:
        if self._hover_leave_callback is not None:
            self._hover_leave_callback()
        if self._base_font is not None:
            self.setFont(QtGui.QFont(self._base_font))
        if self._base_brush is not None:
            self.setBrush(QtGui.QBrush(self._base_brush))
        self.setOpacity(1.0)
        self.setScale(1.0)
        super().hoverLeaveEvent(event)

    def mousePressEvent(self, event: Optional[QtWidgets.QGraphicsSceneMouseEvent]) -> None:
        if event is not None and event.button() == QtCore.Qt.MouseButton.LeftButton:
            self._callback()
            event.accept()
        super().mousePressEvent(event)


class _ActionPathItem(QtWidgets.QGraphicsPathItem):
    def __init__(
        self,
        path: QtGui.QPainterPath,
        callback: Callable[[], None],
        parent: Optional[QtWidgets.QGraphicsItem] = None,
    ) -> None:
        super().__init__(path, parent)
        self._callback = callback
        self.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)

    def mousePressEvent(self, event: Optional[QtWidgets.QGraphicsSceneMouseEvent]) -> None:
        if event is not None and event.button() == QtCore.Qt.MouseButton.LeftButton:
            self._callback()
            event.accept()
        super().mousePressEvent(event)


class _TreeEdgeItem(QtWidgets.QGraphicsPathItem):
    def __init__(
        self,
        path: QtGui.QPainterPath,
        owner: "RaxmlTreeView",
        subtree_leaves: Sequence[str],
        base_pen: QtGui.QPen,
        hover_pen: QtGui.QPen,
    ) -> None:
        super().__init__(path)
        self._owner = owner
        self._subtree_leaves = tuple(str(s) for s in subtree_leaves if str(s).strip())
        self._base_pen = QtGui.QPen(base_pen)
        self._hover_pen = QtGui.QPen(hover_pen)
        self.setPen(self._base_pen)
        self.setAcceptHoverEvents(owner.hover_interactions_enabled())

    def subtree_leaves(self) -> Tuple[str, ...]:
        return self._subtree_leaves

    def set_hovered(self, active: bool) -> None:
        try:
            self.setPen(self._hover_pen if active else self._base_pen)
        except RuntimeError:
            # Underlying C++ item can be deleted after scene.clear() while stale
            # Python references still exist briefly.
            return

    def hoverEnterEvent(self, event: Optional[QtWidgets.QGraphicsSceneHoverEvent]) -> None:
        self._owner._on_subtree_hover_enter(self._subtree_leaves)
        super().hoverEnterEvent(event)

    def hoverLeaveEvent(self, event: Optional[QtWidgets.QGraphicsSceneHoverEvent]) -> None:
        self._owner._on_subtree_hover_leave(self._subtree_leaves)
        super().hoverLeaveEvent(event)


class _NodeMarkerItem(QtWidgets.QGraphicsEllipseItem):
    def __init__(
        self,
        rect: QtCore.QRectF,
        base_pen: QtGui.QPen,
        base_brush: QtGui.QBrush,
        hover_pen: QtGui.QPen,
        hover_brush: QtGui.QBrush,
    ) -> None:
        super().__init__(rect)
        self._base_pen = QtGui.QPen(base_pen)
        self._base_brush = QtGui.QBrush(base_brush)
        self._hover_pen = QtGui.QPen(hover_pen)
        self._hover_brush = QtGui.QBrush(hover_brush)
        self.setPen(self._base_pen)
        self.setBrush(self._base_brush)

    def set_hovered(self, active: bool) -> None:
        try:
            self.setPen(self._hover_pen if active else self._base_pen)
            self.setBrush(self._hover_brush if active else self._base_brush)
        except RuntimeError:
            return


class _MatrixStripItem(QtWidgets.QGraphicsObject):
    def __init__(
        self,
        rect: QtCore.QRectF,
        *,
        title: str,
        rows: Sequence[Tuple[str, float, float]],
        matrix: np.ndarray,
        column_labels: Sequence[str],
        column_sizes: Sequence[int],
        present_color: QtGui.QColor,
        absent_color: QtGui.QColor,
        binned: bool,
        active_samples: Optional[Set[str]] = None,
        parent: Optional[QtWidgets.QGraphicsItem] = None,
    ) -> None:
        super().__init__(parent)
        self._rect = QtCore.QRectF(rect)
        self._title = str(title or "").strip() or "Panaroo genes"
        self._rows = [
            (str(sample).strip(), float(top), float(bottom))
            for sample, top, bottom in rows
            if str(sample).strip() and float(bottom) > float(top)
        ]
        self._row_index = {sample: idx for idx, (sample, _top, _bottom) in enumerate(self._rows)}
        self._matrix = np.asarray(matrix, dtype=np.float32)
        self._column_labels = tuple(str(label or "").strip() for label in column_labels)
        self._column_sizes = tuple(max(1, int(size)) for size in column_sizes)
        self._present_color = QtGui.QColor(present_color)
        self._absent_color = QtGui.QColor(absent_color)
        self._binned = bool(binned)
        self._active_samples = set(active_samples) if active_samples is not None else None
        self._highlight_samples: Set[str] = set()
        self._base_image = self._build_image()
        self.setAcceptHoverEvents(True)

    def boundingRect(self) -> QtCore.QRectF:
        return QtCore.QRectF(self._rect).adjusted(-2.0, -2.0, 2.0, 2.0)

    def set_highlighted_samples(self, samples: Optional[Set[str]]) -> None:
        next_set = set(samples or [])
        if next_set == self._highlight_samples:
            return
        self._highlight_samples = next_set
        self.update()

    def paint(
        self,
        painter: QtGui.QPainter,
        option: Optional[QtWidgets.QStyleOptionGraphicsItem],
        widget: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        del option, widget
        if self._base_image.isNull() or not self._rows:
            return
        painter.save()
        painter.setRenderHint(QtGui.QPainter.RenderHint.SmoothPixmapTransform, False)
        n_cols = max(1, self._base_image.width())
        image_left = self._rect.left() + 1.0
        image_width = max(1.0, self._rect.width() - 2.0)
        for row_idx, (sample, top, bottom) in enumerate(self._rows):
            row_top = top + (1.0 if row_idx == 0 else 0.0)
            row_bottom = bottom - (1.0 if row_idx == len(self._rows) - 1 else 0.0)
            height = row_bottom - row_top
            if height <= 0.0:
                continue
            target = QtCore.QRectF(image_left, row_top, image_width, height)
            source = QtCore.QRectF(0.0, float(row_idx), float(n_cols), 1.0)
            painter.drawImage(target, self._base_image, source)
            if sample in self._highlight_samples:
                pen = QtGui.QPen(QtGui.QColor("#f97316"))
                pen.setWidthF(2.0)
                painter.setPen(pen)
                painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
                painter.drawRect(target.adjusted(0.5, 0.5, -0.5, -0.5))
        border_pen = QtGui.QPen(QtGui.QColor("#cbd5e1"))
        border_pen.setWidthF(1.0)
        painter.setPen(border_pen)
        painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
        painter.drawRect(self._rect.adjusted(0.5, 0.5, -0.5, -0.5))
        painter.restore()

    def hoverMoveEvent(self, event: Optional[QtWidgets.QGraphicsSceneHoverEvent]) -> None:
        if event is None:
            return
        point = event.pos()
        row_idx = self._row_index_at(point.y())
        col_idx = self._column_index_at(point.x())
        if row_idx is None or col_idx is None:
            QtWidgets.QToolTip.hideText()
            super().hoverMoveEvent(event)
            return
        sample, _top, _bottom = self._rows[row_idx]
        label = self._column_labels[col_idx] if col_idx < len(self._column_labels) else f"Bin {col_idx + 1}"
        column_size = self._column_sizes[col_idx] if col_idx < len(self._column_sizes) else 1
        value = float(self._matrix[row_idx, col_idx]) if row_idx < self._matrix.shape[0] and col_idx < self._matrix.shape[1] else 0.0
        if self._binned or column_size > 1:
            present_count = int(round(value * column_size))
            percent = value * 100.0
            status = f"Presence in bin: {percent:.1f}% ({present_count}/{column_size})"
        else:
            status = "Present" if value >= 0.5 else "Absent"
        tooltip = f"{sample}\n{self._title}: {label}\n{status}"
        screen_pos = event.screenPos()
        if isinstance(screen_pos, QtCore.QPointF):
            screen_pos = screen_pos.toPoint()
        QtWidgets.QToolTip.showText(screen_pos, tooltip)
        super().hoverMoveEvent(event)

    def hoverLeaveEvent(self, event: Optional[QtWidgets.QGraphicsSceneHoverEvent]) -> None:
        QtWidgets.QToolTip.hideText()
        super().hoverLeaveEvent(event)

    def _row_index_at(self, y: float) -> Optional[int]:
        for idx, (_sample, top, bottom) in enumerate(self._rows):
            if top <= y <= bottom:
                return idx
        return None

    def _column_index_at(self, x: float) -> Optional[int]:
        if self._rect.width() <= 0.0 or x < self._rect.left() or x > self._rect.right():
            return None
        rel = (x - self._rect.left()) / self._rect.width()
        idx = int(math.floor(rel * self._matrix.shape[1]))
        idx = max(0, min(idx, self._matrix.shape[1] - 1))
        return idx

    def _build_image(self) -> QtGui.QImage:
        rows = len(self._rows)
        cols = int(self._matrix.shape[1]) if self._matrix.ndim == 2 else 0
        if rows <= 0 or cols <= 0:
            return QtGui.QImage()
        image = QtGui.QImage(cols, rows, QtGui.QImage.Format.Format_ARGB32_Premultiplied)
        image.fill(QtCore.Qt.GlobalColor.transparent)
        for row_idx, (sample, _top, _bottom) in enumerate(self._rows):
            row_active = self._active_samples is None or sample in self._active_samples
            for col_idx in range(cols):
                image.setPixelColor(col_idx, row_idx, self._cell_color(float(self._matrix[row_idx, col_idx]), active=row_active))
        return image

    def _cell_color(self, value: float, *, active: bool) -> QtGui.QColor:
        clamped = max(0.0, min(1.0, float(value)))
        if clamped <= 0.0:
            color = QtGui.QColor(self._absent_color)
        elif clamped >= 1.0:
            color = QtGui.QColor(self._present_color)
        else:
            color = self._blend_colors(self._absent_color, self._present_color, clamped)
        if active:
            return color
        faded = self._blend_colors(color, QtGui.QColor("#f8fafc"), 0.72)
        faded.setAlpha(190)
        return faded

    @staticmethod
    def _blend_colors(left: QtGui.QColor, right: QtGui.QColor, weight: float) -> QtGui.QColor:
        w = max(0.0, min(1.0, float(weight)))
        return QtGui.QColor(
            int(round(left.red() + (right.red() - left.red()) * w)),
            int(round(left.green() + (right.green() - left.green()) * w)),
            int(round(left.blue() + (right.blue() - left.blue()) * w)),
            int(round(left.alpha() + (right.alpha() - left.alpha()) * w)),
        )


class RaxmlTreeView(QtWidgets.QGraphicsView):
    viewChanged = QtCore.pyqtSignal()
    leafHovered = QtCore.pyqtSignal(object)
    subtreeHovered = QtCore.pyqtSignal(object)
    leafDoubleClicked = QtCore.pyqtSignal(object)
    cladeClicked = QtCore.pyqtSignal(object)
    clusterClicked = QtCore.pyqtSignal(object)

    def __init__(
        self,
        parent: Optional[QtWidgets.QWidget] = None,
        *,
        layout_mode: str = "circular",
        preserve_branch_lengths: bool = False,
        optimize_groups: bool = False,
        use_group_backgrounds: bool = False,
        show_group_labels: bool = True,
        show_cluster_annotations: bool = True,
        enable_hover_interactions: bool = True,
    ) -> None:
        super().__init__(parent)
        self._scene = QtWidgets.QGraphicsScene(self)
        self.setScene(self._scene)
        self.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        self.setDragMode(QtWidgets.QGraphicsView.DragMode.ScrollHandDrag)
        self.setViewportUpdateMode(QtWidgets.QGraphicsView.ViewportUpdateMode.BoundingRectViewportUpdate)
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.setTransformationAnchor(QtWidgets.QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setResizeAnchor(QtWidgets.QGraphicsView.ViewportAnchor.AnchorViewCenter)
        self.setMinimumWidth(280)
        mode = (layout_mode or "circular").strip().lower()
        self._layout_mode = mode if mode in {"circular", "rectangular"} else "circular"
        self._preserve_branch_lengths = bool(preserve_branch_lengths)
        self._optimize_groups = bool(optimize_groups)
        self._use_group_backgrounds = bool(use_group_backgrounds)
        self._show_group_labels = bool(show_group_labels)
        self._show_cluster_annotations = bool(show_cluster_annotations)
        self._enable_hover_interactions = bool(enable_hover_interactions)

        self._root: Optional[_TreeNode] = None
        self._root_version = 0
        # When True, the whole tree is scaled to fit the viewport on render/resize.
        # Manual wheel-zoom turns this off until the next tree/layout load.
        self._auto_fit = True
        self._active_samples: Optional[set[str]] = None
        self._reference_sample: Optional[str] = None
        self._leaf_group_map: Dict[str, str] = {}
        self._group_color_map: Dict[str, QtGui.QColor] = {}
        self._leaf_display_group_map: Dict[str, str] = {}
        self._display_group_color_map: Dict[str, QtGui.QColor] = {}
        self._leaf_cluster_map: Dict[str, str] = {}
        self._leaf_hover_metadata: Dict[str, Dict[str, str]] = {}
        self._leaf_side_metadata_map: Dict[str, str] = {}
        self._side_metadata_color_map: Dict[str, QtGui.QColor] = {}
        self._side_metadata_title: Optional[str] = None
        self._side_metadata_missing_label = "Missing"
        self._side_metadata_tracks: List[Dict[str, object]] = []
        self._bottom_legend_sections: List[Dict[str, object]] = []
        self._matrix_track_items: List[_MatrixStripItem] = []
        self._render_root_cache_key: Optional[Tuple[int, bool, Tuple[Tuple[str, str], ...]]] = None
        self._render_root_cache: Optional[_TreeNode] = None
        self._leaf_order: List[str] = []
        self._leaf_pos_scene: Dict[str, QtCore.QPointF] = {}
        self._leaf_label_items: Dict[str, _LeafLabelItem] = {}
        self._edge_items_by_leaf: Dict[str, List[_TreeEdgeItem]] = {}
        self._marker_items_by_leaf: Dict[str, List[_NodeMarkerItem]] = {}
        self._hover_leaf: Optional[str] = None
        self._active_leaf_set: Set[str] = set()
        self._active_edge_items: Set[_TreeEdgeItem] = set()
        self._active_marker_items: Set[_NodeMarkerItem] = set()
        self._hover_subtree_key: Optional[Tuple[str, ...]] = None
        self._action_hover_key: Optional[Tuple[str, ...]] = None
        self._label_brush_cache: Dict[Tuple[int, int, int, int], QtGui.QBrush] = {}

        self._base_edge_pen = QtGui.QPen(QtGui.QColor("#000000"))
        self._base_edge_pen.setWidthF(1.2)
        self._base_edge_pen.setCapStyle(QtCore.Qt.PenCapStyle.RoundCap)
        self._base_edge_pen.setJoinStyle(QtCore.Qt.PenJoinStyle.RoundJoin)
        self._hover_edge_pen = QtGui.QPen(QtGui.QColor("#f97316"))
        self._hover_edge_pen.setWidthF(2.4)
        self._hover_edge_pen.setCapStyle(QtCore.Qt.PenCapStyle.RoundCap)
        self._hover_edge_pen.setJoinStyle(QtCore.Qt.PenJoinStyle.RoundJoin)

        self._draw_error("No tree loaded.")

    def hover_interactions_enabled(self) -> bool:
        return self._enable_hover_interactions

    def set_tree_file(self, path: Path) -> None:
        try:
            text = Path(path).read_text(encoding="utf-8")
            self.set_tree_text(text)
        except Exception as exc:
            self._root = None
            self._draw_error(f"Failed to load tree: {exc}")
        
    def set_tree_text(self, text: str) -> None:
        parser = _NewickParser(str(text or ""))
        self._root = parser.parse()
        self._root_version += 1
        self._auto_fit = True  # a freshly loaded tree starts fitted to the window
        self._invalidate_render_root_cache()
        self._render_tree()

    def set_active_samples(self, samples: Optional[Sequence[str]]) -> None:
        if samples is None:
            self._active_samples = None
        else:
            self._active_samples = {self._normalize_label(s) for s in samples if str(s).strip()}
        if self._root is not None:
            self._render_tree()

    def set_reference_sample(self, sample: Optional[str]) -> None:
        value = self._normalize_label(sample) if sample else ""
        self._reference_sample = value or None
        if self._root is not None:
            self._render_tree()

    def set_preserve_branch_lengths(self, enabled: bool) -> None:
        enabled = bool(enabled)
        if enabled == self._preserve_branch_lengths:
            return
        self._preserve_branch_lengths = enabled
        if self._root is not None:
            self._render_tree()

    def set_optimize_groups(self, enabled: bool) -> None:
        enabled = bool(enabled)
        if enabled == self._optimize_groups:
            return
        self._optimize_groups = enabled
        self._invalidate_render_root_cache()
        if self._root is not None:
            self._render_tree()

    def set_use_group_backgrounds(self, enabled: bool) -> None:
        enabled = bool(enabled)
        if enabled == self._use_group_backgrounds:
            return
        self._use_group_backgrounds = enabled
        if self._root is not None:
            self._render_tree()

    def set_leaf_groups(
        self,
        groups: Optional[Mapping[str, object]],
        *,
        group_colors: Optional[Mapping[str, object]] = None,
    ) -> None:
        normalized_groups: Dict[str, str] = {}
        if groups:
            for raw_leaf, raw_group in groups.items():
                leaf_name = self._normalize_label(str(raw_leaf))
                group_name = str(raw_group or "").strip()
                if not leaf_name or not group_name:
                    continue
                normalized_groups[leaf_name] = group_name

        normalized_colors: Dict[str, QtGui.QColor] = {}
        if group_colors:
            for raw_group, raw_color in group_colors.items():
                group_name = str(raw_group or "").strip()
                color = self._coerce_color(raw_color)
                if not group_name or color is None:
                    continue
                normalized_colors[group_name] = color

        if (
            normalized_groups == self._leaf_group_map
            and self._same_group_colors(normalized_colors, self._group_color_map)
        ):
            return
        self._leaf_group_map = normalized_groups
        self._group_color_map = normalized_colors
        self._invalidate_render_root_cache()
        if self._root is not None:
            self._render_tree()

    def set_leaf_display_groups(
        self,
        groups: Optional[Mapping[str, object]],
        *,
        group_colors: Optional[Mapping[str, object]] = None,
    ) -> None:
        normalized_groups: Dict[str, str] = {}
        if groups:
            for raw_leaf, raw_group in groups.items():
                leaf_name = self._normalize_label(str(raw_leaf))
                group_name = str(raw_group or "").strip()
                if not leaf_name or not group_name:
                    continue
                normalized_groups[leaf_name] = group_name

        normalized_colors: Dict[str, QtGui.QColor] = {}
        if group_colors:
            for raw_group, raw_color in group_colors.items():
                group_name = str(raw_group or "").strip()
                color = self._coerce_color(raw_color)
                if not group_name or color is None:
                    continue
                normalized_colors[group_name] = color

        if (
            normalized_groups == self._leaf_display_group_map
            and self._same_group_colors(normalized_colors, self._display_group_color_map)
        ):
            return
        self._leaf_display_group_map = normalized_groups
        self._display_group_color_map = normalized_colors
        if self._root is not None:
            self._render_tree()

    def set_leaf_clusters(self, clusters: Optional[Mapping[str, object]]) -> None:
        normalized_clusters: Dict[str, str] = {}
        if clusters:
            for raw_leaf, raw_cluster in clusters.items():
                leaf_name = self._normalize_label(str(raw_leaf))
                cluster_name = str(raw_cluster or "").strip()
                if not leaf_name or not cluster_name:
                    continue
                normalized_clusters[leaf_name] = cluster_name

        if normalized_clusters == self._leaf_cluster_map:
            return
        self._leaf_cluster_map = normalized_clusters
        if self._root is not None:
            self._render_tree()

    def set_leaf_hover_metadata(
        self,
        metadata: Optional[Mapping[str, Mapping[str, object]]],
    ) -> None:
        normalized: Dict[str, Dict[str, str]] = {}
        for raw_leaf, raw_fields in (metadata or {}).items():
            leaf_name = self._normalize_label(str(raw_leaf))
            if not leaf_name or not isinstance(raw_fields, Mapping):
                continue
            fields: Dict[str, str] = {}
            for raw_label, raw_value in raw_fields.items():
                label = str(raw_label).strip()
                value = str(raw_value).strip() if raw_value is not None else ""
                if not label or not value:
                    continue
                fields[label] = value
            if fields:
                normalized[leaf_name] = fields
        if normalized == self._leaf_hover_metadata:
            return
        self._leaf_hover_metadata = normalized
        self._refresh_leaf_tooltips()

    def set_leaf_side_metadata(
        self,
        values: Optional[Mapping[str, object]],
        *,
        title: Optional[str] = None,
        value_colors: Optional[Mapping[str, object]] = None,
        missing_label: str = "Missing",
    ) -> None:
        if values is None:
            self.set_leaf_side_metadata_tracks(None)
            return
        normalized_values: Dict[str, str] = {}
        if values:
            for raw_leaf, raw_value in values.items():
                leaf_name = self._normalize_label(str(raw_leaf))
                if not leaf_name:
                    continue
                value_text = str(raw_value or "").strip()
                normalized_values[leaf_name] = value_text or str(missing_label or "Missing").strip() or "Missing"

        normalized_colors: Dict[str, QtGui.QColor] = {}
        if value_colors:
            for raw_value, raw_color in value_colors.items():
                value_name = str(raw_value or "").strip()
                color = self._coerce_color(raw_color)
                if not value_name or color is None:
                    continue
                normalized_colors[value_name] = color

        next_title = str(title or "").strip() or None
        next_missing = str(missing_label or "Missing").strip() or "Missing"
        if (
            normalized_values == self._leaf_side_metadata_map
            and self._same_group_colors(normalized_colors, self._side_metadata_color_map)
            and next_title == self._side_metadata_title
            and next_missing == self._side_metadata_missing_label
        ):
            return
        self._leaf_side_metadata_map = normalized_values
        self._side_metadata_color_map = normalized_colors
        self._side_metadata_title = next_title
        self._side_metadata_missing_label = next_missing
        self.set_leaf_side_metadata_tracks(
            [
                {
                    "title": next_title,
                    "values": normalized_values,
                    "value_colors": normalized_colors,
                    "missing_label": next_missing,
                }
            ]
        )

    def set_leaf_side_metadata_tracks(
        self,
        tracks: Optional[Sequence[Mapping[str, object]]],
    ) -> None:
        normalized_tracks: List[Dict[str, object]] = []
        for track in tracks or []:
            mode = str(track.get("mode") or "categorical").strip().lower() or "categorical"
            title = str(track.get("title") or "").strip() or None
            missing_label = str(track.get("missing_label") or "Missing").strip() or "Missing"
            if mode == "heatmap_matrix":
                raw_matrix = track.get("matrix")
                matrix = np.asarray(raw_matrix, dtype=np.float32) if raw_matrix is not None else np.zeros((0, 0), dtype=np.float32)
                if matrix.ndim != 2 or matrix.shape[0] <= 0 or matrix.shape[1] <= 0:
                    continue
                raw_sample_order = track.get("sample_order") or []
                sample_order = tuple(
                    self._normalize_label(str(sample))
                    for sample in raw_sample_order
                    if self._normalize_label(str(sample))
                )
                if len(sample_order) != matrix.shape[0]:
                    continue
                raw_column_labels = track.get("column_labels") or []
                column_labels = tuple(str(label or "").strip() or f"Column {idx + 1}" for idx, label in enumerate(raw_column_labels))
                if len(column_labels) != matrix.shape[1]:
                    continue
                raw_column_sizes = track.get("column_sizes") or []
                column_sizes = tuple(max(1, int(size)) for size in raw_column_sizes)
                if len(column_sizes) != matrix.shape[1]:
                    column_sizes = tuple(1 for _ in range(matrix.shape[1]))
                present_color = self._coerce_color(track.get("present_color")) or QtGui.QColor("#2563eb")
                absent_color = self._coerce_color(track.get("absent_color")) or QtGui.QColor("#e2e8f0")
                normalized_tracks.append(
                    {
                        "mode": mode,
                        "title": title,
                        "matrix": np.array(matrix, dtype=np.float32, copy=True),
                        "sample_order": sample_order,
                        "column_labels": column_labels,
                        "column_sizes": column_sizes,
                        "present_color": QtGui.QColor(present_color),
                        "absent_color": QtGui.QColor(absent_color),
                        "binned": bool(track.get("binned")),
                    }
                )
                continue
            raw_values = track.get("values")
            raw_colors = track.get("value_colors")

            normalized_values: Dict[str, str] = {}
            if isinstance(raw_values, Mapping):
                for raw_leaf, raw_value in raw_values.items():
                    leaf_name = self._normalize_label(str(raw_leaf))
                    if not leaf_name:
                        continue
                    value_text = str(raw_value or "").strip()
                    normalized_values[leaf_name] = value_text or missing_label

            normalized_colors: Dict[str, QtGui.QColor] = {}
            if isinstance(raw_colors, Mapping):
                for raw_value, raw_color in raw_colors.items():
                    value_name = str(raw_value or "").strip()
                    color = self._coerce_color(raw_color)
                    if not value_name or color is None:
                        continue
                    normalized_colors[value_name] = color

            if not normalized_values:
                continue
            normalized_tracks.append(
                {
                    "mode": mode,
                    "title": title,
                    "values": normalized_values,
                    "value_colors": normalized_colors,
                    "missing_label": missing_label,
                }
            )

        if self._same_side_metadata_tracks(normalized_tracks, self._side_metadata_tracks):
            return
        self._side_metadata_tracks = normalized_tracks
        if normalized_tracks:
            first_track = normalized_tracks[0]
            if str(first_track.get("mode") or "categorical") == "heatmap_matrix":
                self._leaf_side_metadata_map = {}
                self._side_metadata_color_map = {}
                self._side_metadata_title = str(first_track.get("title") or "").strip() or None
                self._side_metadata_missing_label = "Missing"
            else:
                self._leaf_side_metadata_map = dict(first_track.get("values", {}))
                self._side_metadata_color_map = dict(first_track.get("value_colors", {}))
                self._side_metadata_title = str(first_track.get("title") or "").strip() or None
                self._side_metadata_missing_label = str(first_track.get("missing_label") or "Missing").strip() or "Missing"
        else:
            self._leaf_side_metadata_map = {}
            self._side_metadata_color_map = {}
            self._side_metadata_title = None
            self._side_metadata_missing_label = "Missing"
        if self._root is not None:
            self._render_tree()

    def set_bottom_legend_sections(
        self,
        sections: Optional[Sequence[Mapping[str, object]]],
    ) -> None:
        normalized_sections: List[Dict[str, object]] = []
        for section in sections or []:
            title = str(section.get("title") or "").strip()
            raw_items = section.get("items")
            if not title or not isinstance(raw_items, Sequence):
                continue
            items: List[Tuple[str, QtGui.QColor]] = []
            for raw_item in raw_items:
                if not isinstance(raw_item, Sequence) or len(raw_item) < 2:
                    continue
                label = str(raw_item[0] or "").strip()
                color = self._coerce_color(raw_item[1])
                if not label or color is None:
                    continue
                items.append((label, QtGui.QColor(color)))
            if items:
                normalized_sections.append({"title": title, "items": items})
        if self._same_bottom_legend_sections(normalized_sections, self._bottom_legend_sections):
            return
        self._bottom_legend_sections = normalized_sections
        if self._root is not None:
            self._render_tree()

    def leaf_order(self) -> List[str]:
        return list(self._leaf_order)

    def export_figure(
        self,
        path: Path,
        *,
        image_scale: float = 4.0,
        image_dpi: int = 300,
    ) -> None:
        target_path = Path(path).expanduser()
        suffix = target_path.suffix.lower()
        if suffix not in {".pdf", ".png", ".jpg", ".jpeg"}:
            raise ValueError(f"Unsupported export format: {target_path.suffix or '<none>'}")

        scene_rect = QtCore.QRectF(self.sceneRect())
        if scene_rect.width() <= 1.0 or scene_rect.height() <= 1.0:
            raise ValueError("Tree scene is empty.")

        if suffix == ".pdf":
            self._export_pdf(target_path, scene_rect, dpi=image_dpi)
            return
        self._export_image(target_path, scene_rect, scale=image_scale, dpi=image_dpi)

    def _export_image(
        self,
        path: Path,
        scene_rect: QtCore.QRectF,
        *,
        scale: float,
        dpi: int,
    ) -> None:
        export_scale = max(1.0, float(scale))
        width = max(1, int(math.ceil(scene_rect.width() * export_scale)))
        height = max(1, int(math.ceil(scene_rect.height() * export_scale)))

        image = QtGui.QImage(width, height, QtGui.QImage.Format.Format_ARGB32_Premultiplied)
        dots_per_meter = int(round(float(dpi) / 0.0254))
        image.setDotsPerMeterX(dots_per_meter)
        image.setDotsPerMeterY(dots_per_meter)
        image.fill(QtGui.QColor("#ffffff"))

        painter = QtGui.QPainter(image)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        painter.setRenderHint(QtGui.QPainter.RenderHint.TextAntialiasing, True)
        painter.setRenderHint(QtGui.QPainter.RenderHint.SmoothPixmapTransform, True)
        try:
            self._scene.render(
                painter,
                QtCore.QRectF(0.0, 0.0, float(width), float(height)),
                scene_rect,
            )
        finally:
            painter.end()

        if not image.save(str(path)):
            raise IOError(f"Failed to save image to {path}")

    def _export_pdf(
        self,
        path: Path,
        scene_rect: QtCore.QRectF,
        *,
        dpi: int,
    ) -> None:
        writer = QtGui.QPdfWriter(str(path))
        writer.setResolution(max(72, int(dpi)))
        writer.setPageMargins(QtCore.QMarginsF(0.0, 0.0, 0.0, 0.0))

        width_inches = max(1.0, float(scene_rect.width())) / 96.0
        height_inches = max(1.0, float(scene_rect.height())) / 96.0
        page_size = QtGui.QPageSize(
            QtCore.QSizeF(width_inches, height_inches),
            QtGui.QPageSize.Unit.Inch,
            "neoTreeExport",
            QtGui.QPageSize.SizeMatchPolicy.ExactMatch,
        )
        writer.setPageSize(page_size)

        painter = QtGui.QPainter(writer)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        painter.setRenderHint(QtGui.QPainter.RenderHint.TextAntialiasing, True)
        painter.setRenderHint(QtGui.QPainter.RenderHint.SmoothPixmapTransform, True)
        try:
            target_rect = QtCore.QRectF(0.0, 0.0, float(writer.width()), float(writer.height()))
            painter.fillRect(target_rect, QtGui.QColor("#ffffff"))
            self._scene.render(painter, target_rect, scene_rect)
        finally:
            painter.end()

    def leaf_viewport_y(self, sample: str) -> Optional[float]:
        key = self._normalize_label(sample)
        pt = self._leaf_pos_scene.get(key)
        if pt is None:
            return None
        view_pt = self.mapFromScene(pt)
        return float(view_pt.y())

    def highlight_sample(self, sample: Optional[str]) -> None:
        key = self._normalize_label(sample) if sample else None
        if key and key not in self._leaf_label_items:
            key = None
        self._hover_subtree_key = None
        self._set_active_leaves([key] if key else None, emit_signals=False)

    def highlight_samples(self, samples: Optional[Sequence[str]]) -> None:
        if not samples:
            self._hover_subtree_key = None
            self._set_active_leaves(None, emit_signals=False)
            return
        normalized = [self._normalize_label(s) for s in samples if str(s).strip()]
        self._hover_subtree_key = None
        self._set_active_leaves(normalized, emit_signals=False)

    def clear_tree(self, message: str = "No tree loaded.") -> None:
        self._root = None
        self._invalidate_render_root_cache()
        self._draw_error(message)

    def _invalidate_render_root_cache(self) -> None:
        self._render_root_cache_key = None
        self._render_root_cache = None

    def _draw_error(self, message: str) -> None:
        self._scene.clear()
        self._leaf_order = []
        self._leaf_pos_scene = {}
        self._leaf_label_items = {}
        self._edge_items_by_leaf = {}
        self._marker_items_by_leaf = {}
        self._hover_leaf = None
        self._active_leaf_set = set()
        self._active_edge_items = set()
        self._active_marker_items = set()
        self._hover_subtree_key = None
        self._action_hover_key = None

        txt = QtWidgets.QGraphicsSimpleTextItem(message)
        txt.setBrush(QtGui.QBrush(QtGui.QColor("#64748b")))
        txt.setPos(8, 8)
        self._scene.addItem(txt)
        rect = QtCore.QRectF(0, 0, max(300.0, txt.boundingRect().width() + 24.0), 120.0)
        self._scene.setSceneRect(rect)
        self.resetTransform()
        self.centerOn(rect.center())
        self.viewChanged.emit()

    def _normalize_label(self, raw: str) -> str:
        s = (raw or "").strip()
        if s.startswith("'") and s.endswith("'") and len(s) >= 2:
            s = s[1:-1]
        return s.strip()

    def _canonical_leaf_name(self, raw: str) -> str:
        name = self._normalize_label(raw)
        if name.lower() == "reference" and self._reference_sample:
            return self._reference_sample
        return name

    def _display_leaf_name(self, canonical_name: str) -> str:
        if self._reference_sample and canonical_name == self._reference_sample:
            return f"{canonical_name} (reference)"
        return canonical_name

    @staticmethod
    def _coerce_color(value: object) -> Optional[QtGui.QColor]:
        if isinstance(value, QtGui.QColor):
            return QtGui.QColor(value)
        color = QtGui.QColor(str(value))
        if not color.isValid():
            return None
        return color

    @staticmethod
    def _same_group_colors(
        left: Mapping[str, QtGui.QColor],
        right: Mapping[str, QtGui.QColor],
    ) -> bool:
        if set(left.keys()) != set(right.keys()):
            return False
        for key, color in left.items():
            other = right.get(key)
            if other is None or color.rgba() != other.rgba():
                return False
        return True

    def _leaf_group(self, leaf_name: str) -> Optional[str]:
        group = self._leaf_group_map.get(leaf_name)
        return group if group else None

    def _group_color(self, group_name: Optional[str]) -> QtGui.QColor:
        if not group_name:
            return QtGui.QColor("#000000")
        color = self._group_color_map.get(group_name)
        return QtGui.QColor(color) if color is not None else QtGui.QColor("#000000")

    def _display_group(self, leaf_name: str) -> Optional[str]:
        group = self._leaf_display_group_map.get(leaf_name)
        return group if group else None

    def _display_group_color(self, group_name: Optional[str]) -> QtGui.QColor:
        if not group_name:
            return QtGui.QColor("#000000")
        color = self._display_group_color_map.get(group_name)
        return QtGui.QColor(color) if color is not None else QtGui.QColor("#000000")

    @staticmethod
    def _blend_colors(base: QtGui.QColor, target: QtGui.QColor, ratio: float) -> QtGui.QColor:
        weight = max(0.0, min(1.0, float(ratio)))
        keep = 1.0 - weight
        return QtGui.QColor(
            int(round(base.red() * keep + target.red() * weight)),
            int(round(base.green() * keep + target.green() * weight)),
            int(round(base.blue() * keep + target.blue() * weight)),
            int(round(base.alpha() * keep + target.alpha() * weight)),
        )

    def _subtree_active_state(self, subtree_leaves: Sequence[str]) -> str:
        if self._active_samples is None or not subtree_leaves:
            return "all"
        active_count = sum(1 for leaf in subtree_leaves if leaf in self._active_samples)
        if active_count <= 0:
            return "inactive"
        if active_count < len(subtree_leaves):
            return "partial"
        return "active"

    def _visual_color_for_subtree(
        self,
        subtree_leaves: Sequence[str],
        active_color: QtGui.QColor,
        *,
        inactive_color: str = "#cbd5e1",
        partial_color: str = "#94a3b8",
        partial_blend: float = 0.45,
    ) -> QtGui.QColor:
        state = self._subtree_active_state(subtree_leaves)
        if state == "inactive":
            return QtGui.QColor(inactive_color)
        if state == "partial":
            partial = QtGui.QColor(partial_color)
            return self._blend_colors(active_color, partial, partial_blend)
        return QtGui.QColor(active_color)

    def _leaf_base_color(self, leaf_name: str) -> QtGui.QColor:
        display_group = self._display_group(leaf_name)
        if display_group:
            return self._display_group_color(display_group)
        if self._use_group_backgrounds:
            return QtGui.QColor("#000000")
        return self._group_color(self._leaf_group(leaf_name))

    def _leaf_cluster(self, leaf_name: str) -> Optional[str]:
        cluster = self._leaf_cluster_map.get(leaf_name)
        return cluster if cluster else None

    def _leaf_side_metadata(self, leaf_name: str) -> str:
        value = self._leaf_side_metadata_map.get(leaf_name, self._side_metadata_missing_label)
        return value or self._side_metadata_missing_label

    def _side_metadata_color(self, value: str) -> QtGui.QColor:
        color = self._side_metadata_color_map.get(value)
        if color is None:
            return QtGui.QColor("#e2e8f0")
        return QtGui.QColor(color)

    def _same_side_metadata_tracks(
        self,
        left: Sequence[Mapping[str, object]],
        right: Sequence[Mapping[str, object]],
    ) -> bool:
        if len(left) != len(right):
            return False
        for left_track, right_track in zip(left, right):
            left_mode = str(left_track.get("mode") or "categorical")
            right_mode = str(right_track.get("mode") or "categorical")
            if left_mode != right_mode:
                return False
            if str(left_track.get("title") or "") != str(right_track.get("title") or ""):
                return False
            if left_mode == "heatmap_matrix":
                if tuple(left_track.get("sample_order", ())) != tuple(right_track.get("sample_order", ())):
                    return False
                if tuple(left_track.get("column_labels", ())) != tuple(right_track.get("column_labels", ())):
                    return False
                if tuple(left_track.get("column_sizes", ())) != tuple(right_track.get("column_sizes", ())):
                    return False
                if bool(left_track.get("binned")) != bool(right_track.get("binned")):
                    return False
                if self._coerce_color(left_track.get("present_color")) != self._coerce_color(right_track.get("present_color")):
                    return False
                if self._coerce_color(left_track.get("absent_color")) != self._coerce_color(right_track.get("absent_color")):
                    return False
                left_matrix = np.asarray(left_track.get("matrix"), dtype=np.float32)
                right_matrix = np.asarray(right_track.get("matrix"), dtype=np.float32)
                if left_matrix.shape != right_matrix.shape or not np.array_equal(left_matrix, right_matrix):
                    return False
                continue
            if str(left_track.get("missing_label") or "") != str(right_track.get("missing_label") or ""):
                return False
            if dict(left_track.get("values", {})) != dict(right_track.get("values", {})):
                return False
            if not self._same_group_colors(
                dict(left_track.get("value_colors", {})),
                dict(right_track.get("value_colors", {})),
            ):
                return False
        return True

    def _same_bottom_legend_sections(
        self,
        left: Sequence[Mapping[str, object]],
        right: Sequence[Mapping[str, object]],
    ) -> bool:
        if len(left) != len(right):
            return False
        for left_section, right_section in zip(left, right):
            if str(left_section.get("title") or "") != str(right_section.get("title") or ""):
                return False
            left_items = list(left_section.get("items", []))
            right_items = list(right_section.get("items", []))
            if len(left_items) != len(right_items):
                return False
            for (left_label, left_color), (right_label, right_color) in zip(left_items, right_items):
                if str(left_label or "") != str(right_label or ""):
                    return False
                if self._coerce_color(left_color) != self._coerce_color(right_color):
                    return False
        return True

    def _subtree_group(self, subtree_leaves: Sequence[str]) -> Optional[str]:
        current: Optional[str] = None
        for leaf_name in subtree_leaves:
            group = self._leaf_group(leaf_name)
            if not group:
                return None
            if current is None:
                current = group
                continue
            if group != current:
                return None
        return current

    def _subtree_display_group(self, subtree_leaves: Sequence[str]) -> Optional[str]:
        current: Optional[str] = None
        for leaf_name in subtree_leaves:
            group = self._display_group(leaf_name)
            if not group:
                return None
            if current is None:
                current = group
                continue
            if group != current:
                return None
        return current

    def _edge_color_for_subtree(self, subtree_leaves: Sequence[str]) -> QtGui.QColor:
        display_group = self._subtree_display_group(subtree_leaves)
        if display_group:
            return self._visual_color_for_subtree(subtree_leaves, self._display_group_color(display_group))
        if self._use_group_backgrounds:
            return self._visual_color_for_subtree(subtree_leaves, QtGui.QColor("#000000"))
        group = self._subtree_group(subtree_leaves)
        if group:
            return self._visual_color_for_subtree(subtree_leaves, self._group_color(group))
        return self._visual_color_for_subtree(subtree_leaves, QtGui.QColor("#000000"))

    def _leaf_tooltip(self, leaf_name: str, label_text: str) -> str:
        fields = self._leaf_hover_metadata.get(leaf_name, {})
        if not fields:
            return label_text
        lines = [label_text]
        for label, value in fields.items():
            lines.append(f"{label}: {value}")
        return "\n".join(lines)

    def _refresh_leaf_tooltips(self) -> None:
        for leaf_name, label_item in self._leaf_label_items.items():
            label_text = self._display_leaf_name(leaf_name)
            label_item.setToolTip(self._leaf_tooltip(leaf_name, label_text))

    def _group_background_brush(
        self,
        group_name: Optional[str],
        group_leaves: Sequence[str],
    ) -> QtGui.QBrush:
        color = self._group_color(group_name)
        state = self._subtree_active_state(group_leaves)
        if state == "inactive":
            color = QtGui.QColor("#e2e8f0")
            color.setAlpha(34)
        elif state == "partial":
            color.setAlpha(24)
        else:
            color.setAlpha(56)
        return QtGui.QBrush(color)

    def _edge_pen_for_color(self, color: QtGui.QColor, *, width: float = 1.2) -> QtGui.QPen:
        pen = QtGui.QPen(color)
        pen.setWidthF(width)
        pen.setCapStyle(QtCore.Qt.PenCapStyle.RoundCap)
        pen.setJoinStyle(QtCore.Qt.PenJoinStyle.RoundJoin)
        return pen

    def _group_signature(self) -> Tuple[Tuple[str, str], ...]:
        return tuple(sorted(self._leaf_group_map.items()))

    def _should_optimize_groups(self) -> bool:
        if not self._optimize_groups:
            return False
        unique = {group for group in self._leaf_group_map.values() if str(group).strip()}
        return len(unique) >= 2

    def _clone_tree(self, node: _TreeNode) -> _TreeNode:
        cloned = _TreeNode(
            name=str(node.name or ""),
            branch_length=node.branch_length,
        )
        cloned.children = [self._clone_tree(child) for child in node.children]
        return cloned

    def _group_for_leaf_sequence(self, leaf_name: str) -> str:
        group = self._leaf_group(leaf_name)
        return group if group is not None else ""

    def _score_group_sequence(self, groups: Sequence[str]) -> float:
        if not groups:
            return 0.0
        score = 0.0
        transitions = 0
        missing = 0
        prev: Optional[str] = None
        run_len = 0
        for raw_group in groups:
            group = str(raw_group or "")
            if not group:
                missing += 1
            if prev is None:
                prev = group
                run_len = 1
                continue
            if group == prev:
                run_len += 1
                continue
            if prev:
                score += float(run_len * run_len)
            transitions += 1
            prev = group
            run_len = 1
        if prev:
            score += float(run_len * run_len)
        score -= float(transitions) * 1.35
        score -= float(missing) * 0.75
        return score

    def _child_sort_key(self, summary: _GroupSummary) -> Tuple[int, str, float, int, str]:
        dominant = summary.dominant or ""
        first_leaf = summary.leaf_order[0] if summary.leaf_order else ""
        return (
            1 if not dominant else 0,
            dominant,
            -summary.purity,
            -summary.total,
            first_leaf,
        )

    def _node_structural_bonus(self, child_summaries: Sequence[_GroupSummary]) -> float:
        if len(child_summaries) <= 1:
            return 0.0
        child_dominants = [summary.dominant for summary in child_summaries if summary.dominant]
        distinct = len(set(child_dominants))
        weighted_purity = sum(summary.purity * float(max(1, summary.total)) for summary in child_summaries)
        mixed_penalty = sum((1.0 - summary.purity) * float(max(1, summary.total)) for summary in child_summaries)
        bonus = 0.0
        if distinct > 1:
            bonus += weighted_purity * (1.0 + 0.35 * float(distinct - 1))
        else:
            bonus -= mixed_penalty * 0.6
        if (
            len(child_summaries) == 2
            and child_summaries[0].dominant
            and child_summaries[1].dominant
            and child_summaries[0].dominant != child_summaries[1].dominant
        ):
            bonus += 0.5 * float(child_summaries[0].total + child_summaries[1].total)
        return bonus

    def _merge_group_summaries(
        self,
        child_pairs: Sequence[Tuple[_TreeNode, _GroupSummary]],
    ) -> _GroupSummary:
        counts: Dict[str, int] = {}
        total = 0
        leaf_order: List[str] = []
        group_order: List[str] = []
        structural_score = 0.0
        child_summaries = [summary for _, summary in child_pairs]
        for summary in child_summaries:
            total += summary.total
            leaf_order.extend(summary.leaf_order)
            group_order.extend(summary.group_order)
            structural_score += summary.structural_score
            for group, count in summary.counts.items():
                counts[group] = counts.get(group, 0) + int(count)

        dominant: Optional[str] = None
        purity = 0.0
        if counts:
            dominant = max(sorted(counts), key=lambda key: (counts[key], key))
            purity = float(counts[dominant]) / float(max(1, total))

        structural_score += self._node_structural_bonus(child_summaries)
        return _GroupSummary(
            counts=counts,
            total=total,
            dominant=dominant,
            purity=purity,
            leaf_order=leaf_order,
            group_order=group_order,
            structural_score=structural_score,
        )

    def _arrangement_score(
        self,
        child_pairs: Sequence[Tuple[_TreeNode, _GroupSummary]],
    ) -> Tuple[float, int, float, str]:
        merged = self._merge_group_summaries(child_pairs)
        sequence_score = self._score_group_sequence(merged.group_order)
        total_score = sequence_score + merged.structural_score
        distinct = len({summary.dominant for _, summary in child_pairs if summary.dominant})
        first_leaf = merged.leaf_order[0] if merged.leaf_order else ""
        return total_score, distinct, merged.purity, first_leaf

    def _optimize_child_pairs(
        self,
        child_pairs: Sequence[Tuple[_TreeNode, _GroupSummary]],
    ) -> List[Tuple[_TreeNode, _GroupSummary]]:
        if len(child_pairs) <= 1:
            return list(child_pairs)
        if len(child_pairs) <= 6:
            best_pairs: Optional[List[Tuple[_TreeNode, _GroupSummary]]] = None
            best_score: Optional[Tuple[float, int, float, str]] = None
            for perm in itertools.permutations(child_pairs):
                perm_list = list(perm)
                score = self._arrangement_score(perm_list)
                if best_score is None or score > best_score:
                    best_score = score
                    best_pairs = perm_list
            return best_pairs if best_pairs is not None else list(child_pairs)

        ordered = sorted(child_pairs, key=lambda pair: self._child_sort_key(pair[1]))
        best: List[Tuple[_TreeNode, _GroupSummary]] = [ordered[0]]
        for pair in ordered[1:]:
            chosen: Optional[List[Tuple[_TreeNode, _GroupSummary]]] = None
            chosen_score: Optional[Tuple[float, int, float, str]] = None
            for idx in range(len(best) + 1):
                candidate = best[:idx] + [pair] + best[idx:]
                score = self._arrangement_score(candidate)
                if chosen_score is None or score > chosen_score:
                    chosen_score = score
                    chosen = candidate
            best = chosen if chosen is not None else best
        return best

    def _optimize_tree_for_groups(self, node: _TreeNode) -> _GroupSummary:
        if node.is_leaf:
            leaf_name = self._canonical_leaf_name(node.name or "")
            group_name = self._group_for_leaf_sequence(leaf_name)
            counts = {group_name: 1} if group_name else {}
            return _GroupSummary(
                counts=counts,
                total=1,
                dominant=group_name or None,
                purity=1.0 if group_name else 0.0,
                leaf_order=[leaf_name],
                group_order=[group_name],
                structural_score=0.0,
            )

        child_pairs: List[Tuple[_TreeNode, _GroupSummary]] = []
        for child in node.children:
            child_pairs.append((child, self._optimize_tree_for_groups(child)))
        child_pairs = self._optimize_child_pairs(child_pairs)
        node.children = [child for child, _ in child_pairs]
        return self._merge_group_summaries(child_pairs)

    def _tree_to_graph(
        self,
        node: _TreeNode,
    ) -> Tuple[Dict[int, List[Tuple[int, float]]], Dict[int, str], List[Tuple[int, int, float]]]:
        adjacency: Dict[int, List[Tuple[int, float]]] = {}
        names: Dict[int, str] = {}
        edges: List[Tuple[int, int, float]] = []

        def walk(current: _TreeNode) -> None:
            node_id = id(current)
            names[node_id] = (
                self._canonical_leaf_name(current.name or "")
                if current.is_leaf
                else self._normalize_label(current.name or "")
            )
            adjacency.setdefault(node_id, [])
            for child in current.children:
                child_id = id(child)
                length = child.branch_length
                if length is None or not math.isfinite(length) or length < 0.0:
                    length = 0.0
                adjacency.setdefault(child_id, [])
                adjacency[node_id].append((child_id, float(length)))
                adjacency[child_id].append((node_id, float(length)))
                edges.append((node_id, child_id, float(length)))
                walk(child)

        walk(node)
        return adjacency, names, edges

    def _build_rooted_from_graph(
        self,
        node_id: int,
        parent_id: Optional[int],
        adjacency: Mapping[int, Sequence[Tuple[int, float]]],
        names: Mapping[int, str],
        branch_length: Optional[float],
    ) -> _TreeNode:
        node = _TreeNode(
            name=str(names.get(node_id, "")),
            branch_length=branch_length,
        )
        for neighbor_id, length in adjacency.get(node_id, ()):
            if parent_id is not None and neighbor_id == parent_id:
                continue
            child = self._build_rooted_from_graph(
                neighbor_id,
                node_id,
                adjacency,
                names,
                float(length),
            )
            node.children.append(child)
        return node

    def _root_graph_on_edge(
        self,
        left_id: int,
        right_id: int,
        edge_length: float,
        adjacency: Mapping[int, Sequence[Tuple[int, float]]],
        names: Mapping[int, str],
    ) -> _TreeNode:
        left_branch = max(0.0, float(edge_length) * 0.5)
        right_branch = max(0.0, float(edge_length) - left_branch)
        root = _TreeNode(name="", branch_length=None)
        root.children = [
            self._build_rooted_from_graph(left_id, right_id, adjacency, names, left_branch),
            self._build_rooted_from_graph(right_id, left_id, adjacency, names, right_branch),
        ]
        return root

    def _candidate_score(self, summary: _GroupSummary) -> Tuple[float, int, float, int]:
        distinct_groups = len({group for group in summary.group_order if group})
        missing_groups = sum(1 for group in summary.group_order if not group)
        total_score = self._score_group_sequence(summary.group_order) + summary.structural_score
        return total_score, distinct_groups, summary.purity, -missing_groups

    def _resolved_render_root(self) -> _TreeNode:
        if self._root is None:
            raise RuntimeError("No tree loaded.")
        if not self._should_optimize_groups():
            return self._root

        cache_key = (
            self._root_version,
            self._optimize_groups,
            self._group_signature(),
        )
        if self._render_root_cache_key == cache_key and self._render_root_cache is not None:
            return self._render_root_cache

        adjacency, names, edges = self._tree_to_graph(self._root)
        best_root = self._clone_tree(self._root)
        best_summary = self._optimize_tree_for_groups(best_root)
        best_score = self._candidate_score(best_summary)

        for left_id, right_id, edge_length in edges:
            candidate_root = self._root_graph_on_edge(
                left_id,
                right_id,
                edge_length,
                adjacency,
                names,
            )
            candidate_summary = self._optimize_tree_for_groups(candidate_root)
            candidate_score = self._candidate_score(candidate_summary)
            if candidate_score > best_score:
                best_root = candidate_root
                best_score = candidate_score

        self._render_root_cache_key = cache_key
        self._render_root_cache = best_root
        return best_root

    def _collect_depths(
        self,
        node: _TreeNode,
        depth: int,
        depth_map: Dict[int, int],
        leaves: List[_TreeNode],
    ) -> None:
        depth_map[id(node)] = depth
        if node.is_leaf:
            leaves.append(node)
            return
        for ch in node.children:
            self._collect_depths(ch, depth + 1, depth_map, leaves)

    def _collect_branch_distances(
        self,
        node: _TreeNode,
        cumulative: float,
        distance_map: Dict[int, float],
    ) -> float:
        distance_map[id(node)] = cumulative
        max_distance = cumulative
        for child in node.children:
            branch_length = child.branch_length
            if branch_length is None or not math.isfinite(branch_length) or branch_length < 0.0:
                branch_length = 0.0
            max_distance = max(
                max_distance,
                self._collect_branch_distances(child, cumulative + float(branch_length), distance_map),
            )
        return max_distance

    def _collect_subtree_leaf_names(
        self,
        node: _TreeNode,
        leaf_name_map: Dict[int, str],
        out: Dict[int, List[str]],
    ) -> List[str]:
        if node.is_leaf:
            name = leaf_name_map.get(id(node), "")
            leaves = [name] if name else []
            out[id(node)] = leaves
            return leaves
        acc: List[str] = []
        for child in node.children:
            acc.extend(self._collect_subtree_leaf_names(child, leaf_name_map, out))
        out[id(node)] = acc
        return acc

    def _build_layout_tree(
        self,
        node: _TreeNode,
        *,
        parent: Optional[_LayoutNode] = None,
        depth: int = 0,
        number: int = 1,
        left_sibling: Optional[_LayoutNode] = None,
        leftmost_sibling: Optional[_LayoutNode] = None,
        node_map: Optional[Dict[int, _LayoutNode]] = None,
    ) -> _LayoutNode:
        layout = _LayoutNode(
            node=node,
            parent=parent,
            depth=depth,
            number=number,
            left_sibling=left_sibling,
            leftmost_sibling=leftmost_sibling,
        )
        layout.ancestor = layout
        if node_map is not None:
            node_map[id(node)] = layout

        first_child: Optional[_LayoutNode] = None
        prev_child: Optional[_LayoutNode] = None
        for idx, child in enumerate(node.children, start=1):
            child_layout = self._build_layout_tree(
                child,
                parent=layout,
                depth=depth + 1,
                number=idx,
                left_sibling=prev_child,
                leftmost_sibling=first_child,
                node_map=node_map,
            )
            layout.children.append(child_layout)
            if first_child is None:
                first_child = child_layout
            prev_child = child_layout
        return layout

    def _move_subtree(self, wl: _LayoutNode, wr: _LayoutNode, shift: float) -> None:
        subtrees = max(1, wr.number - wl.number)
        delta = shift / float(subtrees)
        wr.change -= delta
        wr.shift += shift
        wl.change += delta
        wr.prelim += shift
        wr.mod += shift

    def _execute_shifts(self, v: _LayoutNode) -> None:
        shift = 0.0
        change = 0.0
        for w in reversed(v.children):
            w.prelim += shift
            w.mod += shift
            change += w.change
            shift += w.shift + change

    def _ancestor(self, vil: _LayoutNode, v: _LayoutNode, default_ancestor: _LayoutNode) -> _LayoutNode:
        parent = v.parent
        if parent is None:
            return default_ancestor
        if vil.ancestor is not None and vil.ancestor in parent.children:
            return vil.ancestor
        return default_ancestor

    def _apportion(self, v: _LayoutNode, default_ancestor: _LayoutNode, distance: float) -> _LayoutNode:
        w = v.left_sibling
        if w is None:
            return default_ancestor

        vir = v
        vor = v
        vil = w
        vol = v.leftmost_sibling
        if vol is None:
            return default_ancestor

        sir = v.mod
        sor = v.mod
        sil = vil.mod
        sol = vol.mod

        while vil.right() is not None and vir.left() is not None:
            vil = vil.right()  # type: ignore[assignment]
            vir = vir.left()  # type: ignore[assignment]
            vol_next = vol.left()
            vor_next = vor.right()
            if vol_next is None or vor_next is None:
                break
            vol = vol_next
            vor = vor_next
            vor.ancestor = v
            shift = (vil.prelim + sil) - (vir.prelim + sir) + distance
            if shift > 0.0:
                move_ancestor = self._ancestor(vil, v, default_ancestor)
                self._move_subtree(move_ancestor, v, shift)
                sir += shift
                sor += shift
            sil += vil.mod
            sir += vir.mod
            sol += vol.mod
            sor += vor.mod

        if vil.right() is not None and vor.right() is None:
            vor.thread = vil.right()
            vor.mod += sil - sor
        elif vir.left() is not None and vol.left() is None:
            vol.thread = vir.left()
            vol.mod += sir - sol
            default_ancestor = v
        return default_ancestor

    def _first_walk(self, v: _LayoutNode, distance: float = 1.0) -> _LayoutNode:
        if v.is_leaf:
            if v.left_sibling is not None:
                v.prelim = v.left_sibling.prelim + distance
            else:
                v.prelim = 0.0
            return v

        default_ancestor = v.children[0]
        for child in v.children:
            self._first_walk(child, distance=distance)
            default_ancestor = self._apportion(child, default_ancestor, distance)

        self._execute_shifts(v)
        midpoint = 0.5 * (v.children[0].prelim + v.children[-1].prelim)
        if v.left_sibling is not None:
            v.prelim = v.left_sibling.prelim + distance
            v.mod = v.prelim - midpoint
        else:
            v.prelim = midpoint
        return v

    def _second_walk(self, v: _LayoutNode, m: float = 0.0, depth: int = 0, min_x: Optional[float] = None) -> float:
        v.x = v.prelim + m
        v.y = float(depth)
        current_min = v.x if min_x is None else min(min_x, v.x)
        for child in v.children:
            current_min = self._second_walk(child, m + v.mod, depth + 1, current_min)
        return current_min

    def _third_walk(self, v: _LayoutNode, shift: float) -> None:
        v.x += shift
        for child in v.children:
            self._third_walk(child, shift)

    def _iter_layout_nodes(self, root: _LayoutNode) -> List[_LayoutNode]:
        out: List[_LayoutNode] = []
        stack: List[_LayoutNode] = [root]
        while stack:
            cur = stack.pop()
            out.append(cur)
            for child in reversed(cur.children):
                stack.append(child)
        return out

    def _collect_layout_leaves(self, node: _LayoutNode, leaves: List[_LayoutNode]) -> None:
        if node.is_leaf:
            leaves.append(node)
            return
        for child in node.children:
            self._collect_layout_leaves(child, leaves)

    def _render_tree(self) -> None:
        if self._root is None:
            self._draw_error("No tree loaded.")
            return

        # Preserve logical highlight state across redraw, but never reuse old
        # QGraphicsItem pointers because scene.clear() deletes them.
        pending_active = list(self._active_leaf_set)
        pending_subtree = self._hover_subtree_key
        self._scene.clear()
        self._leaf_order = []
        self._leaf_pos_scene = {}
        self._leaf_label_items = {}
        self._edge_items_by_leaf = {}
        self._marker_items_by_leaf = {}
        self._matrix_track_items = []
        self._hover_leaf = None
        self._active_leaf_set = set()
        self._active_edge_items = set()
        self._active_marker_items = set()
        self._hover_subtree_key = None
        self._action_hover_key = None

        try:
            render_root = self._resolved_render_root()
        except Exception as exc:
            self._draw_error(f"Failed to render tree: {exc}")
            return

        depth_map: Dict[int, int] = {}
        leaves: List[_TreeNode] = []
        self._collect_depths(render_root, 0, depth_map, leaves)
        if not leaves:
            self._draw_error("Tree has no leaves.")
            return

        leaf_name_map: Dict[int, str] = {}
        for i, leaf in enumerate(leaves, start=1):
            raw = leaf.name or f"leaf_{i}"
            leaf_name_map[id(leaf)] = self._canonical_leaf_name(raw)

        if self._layout_mode == "rectangular":
            self._render_tree_rectangular(render_root, depth_map, leaves, leaf_name_map)
        else:
            self._render_tree_circular(render_root, depth_map, leaves, leaf_name_map)

        if pending_subtree:
            self._set_active_leaves(list(pending_subtree), emit_signals=False)
            self._hover_subtree_key = tuple(self._ordered_leaf_subset(pending_subtree))
        elif pending_active:
            self._set_active_leaves(pending_active, emit_signals=False)

        if self._auto_fit:
            self._fit_to_viewport()

        self.viewChanged.emit()

    def _render_tree_circular(
        self,
        root: _TreeNode,
        depth_map: Dict[int, int],
        leaves: List[_TreeNode],
        leaf_name_map: Dict[int, str],
    ) -> None:
        n_leaves = len(leaves)
        max_depth = max(depth_map.values()) if depth_map else 0
        if n_leaves == 1:
            leaf_angles = {id(leaves[0]): -math.pi / 2.0}
        else:
            start_angle = -math.pi / 2.0
            step = (2.0 * math.pi) / float(n_leaves)
            leaf_angles = {id(leaf): start_angle + i * step for i, leaf in enumerate(leaves)}

        angle_map: Dict[int, float] = {}

        def assign_angle(node: _TreeNode) -> float:
            if node.is_leaf:
                ang = leaf_angles[id(node)]
                angle_map[id(node)] = ang
                return ang
            child_angles = [assign_angle(ch) for ch in node.children]
            sx = sum(math.cos(a) for a in child_angles)
            sy = sum(math.sin(a) for a in child_angles)
            ang = math.atan2(sy, sx)
            angle_map[id(node)] = ang
            return ang

        assign_angle(root)

        inner_r = 40.0
        depth_step = 62.0
        if max_depth > 14:
            depth_step = 44.0

        radius_map: Dict[int, float] = {
            node_id: inner_r + float(depth) * depth_step
            for node_id, depth in depth_map.items()
        }

        def point(radius: float, angle: float) -> QtCore.QPointF:
            return QtCore.QPointF(math.cos(angle) * radius, math.sin(angle) * radius)

        def to_qt_degrees(angle_radians: float) -> float:
            # Our angle convention is clockwise (because +Y points down in scene coords).
            # QPainterPath.arcTo expects Qt's counter-clockwise angle convention.
            return -math.degrees(angle_radians)

        def normalized_sweep(start_deg: float, end_deg: float) -> float:
            sweep = end_deg - start_deg
            while sweep <= -180.0:
                sweep += 360.0
            while sweep > 180.0:
                sweep -= 360.0
            return sweep

        subtree_leaf_map: Dict[int, List[str]] = {}
        self._collect_subtree_leaf_names(root, leaf_name_map, subtree_leaf_map)

        def draw_edges(node: _TreeNode, path_stack: List[_TreeEdgeItem]) -> None:
            parent_angle = angle_map[id(node)]
            parent_r = radius_map[id(node)]
            parent_point = point(parent_r, parent_angle)
            for child in node.children:
                child_angle = angle_map[id(child)]
                child_r = radius_map[id(child)]
                child_point = point(child_r, child_angle)

                path = QtGui.QPainterPath(parent_point)
                if abs(child_angle - parent_angle) > 1e-8:
                    rect = QtCore.QRectF(-parent_r, -parent_r, 2.0 * parent_r, 2.0 * parent_r)
                    start_deg = to_qt_degrees(parent_angle)
                    end_deg = to_qt_degrees(child_angle)
                    sweep_deg = normalized_sweep(start_deg, end_deg)
                    path.arcTo(rect, start_deg, sweep_deg)
                path.lineTo(child_point)
                subtree_leaves = subtree_leaf_map.get(id(child), [])
                item = _TreeEdgeItem(
                    path,
                    self,
                    subtree_leaves,
                    self._edge_pen_for_color(self._edge_color_for_subtree(subtree_leaves)),
                    self._hover_edge_pen,
                )
                self._scene.addItem(item)

                stack = path_stack + [item]
                if child.is_leaf:
                    leaf_name = leaf_name_map[id(child)]
                    self._edge_items_by_leaf[leaf_name] = stack
                draw_edges(child, stack)

        draw_edges(root, [])

        outer_r = inner_r + float(max_depth) * depth_step
        label_r = outer_r + 28.0

        for leaf in leaves:
            leaf_name = leaf_name_map[id(leaf)]
            leaf_angle = angle_map[id(leaf)]
            leaf_r = radius_map[id(leaf)]
            leaf_point = point(leaf_r, leaf_angle)
            self._leaf_order.append(leaf_name)
            self._leaf_pos_scene[leaf_name] = leaf_point

            label_text = self._display_leaf_name(leaf_name)
            label_item = _LeafLabelItem(label_text, leaf_name, self)
            label_item.setBrush(QtGui.QBrush(self._label_color(leaf_name)))
            label_item.setToolTip(self._leaf_tooltip(leaf_name, label_text))
            br = label_item.boundingRect()
            label_point = point(label_r, leaf_angle)
            if math.cos(leaf_angle) >= 0.0:
                label_item.setPos(label_point.x() + 5.0, label_point.y() - br.height() / 2.0)
            else:
                label_item.setPos(label_point.x() - br.width() - 5.0, label_point.y() - br.height() / 2.0)
            self._scene.addItem(label_item)
            self._leaf_label_items[leaf_name] = label_item

        bounds = self._scene.itemsBoundingRect().adjusted(-20.0, -20.0, 20.0, 20.0)
        self._scene.setSceneRect(bounds)
        self.resetTransform()
        self.centerOn(bounds.center())

    def _render_tree_rectangular(
        self,
        root: _TreeNode,
        depth_map: Dict[int, int],
        leaves: List[_TreeNode],
        leaf_name_map: Dict[int, str],
    ) -> None:
        layout_map: Dict[int, _LayoutNode] = {}
        layout_root = self._build_layout_tree(root, node_map=layout_map)
        self._first_walk(layout_root, distance=1.0)
        min_x = self._second_walk(layout_root)
        if min_x < 0.0:
            self._third_walk(layout_root, -min_x)

        layout_leaves: List[_LayoutNode] = []
        self._collect_layout_leaves(layout_root, layout_leaves)
        if not layout_leaves:
            self._draw_error("Tree has no leaves.")
            return

        n_leaves = len(layout_leaves)
        max_depth = max((node.depth for node in layout_map.values()), default=0)

        left_margin = 28.0
        top_margin = 24.0

        # Keep labels readable and vertically ordered.
        leaf_step = 34.0
        if n_leaves > 100:
            leaf_step = 30.0
        if n_leaves > 180:
            leaf_step = 26.0
        if n_leaves > 300:
            leaf_step = 22.0
        if n_leaves > 500:
            leaf_step = 18.0
        label_min_step = float(QtGui.QFontMetrics(self.font()).height()) + 4.0
        if leaf_step < label_min_step:
            leaf_step = label_min_step

        tree_height = max(220.0, float(max(1, n_leaves - 1)) * leaf_step)
        target_width = max(260.0, tree_height * 0.38)
        depth_step = target_width / float(max(1, max_depth))
        depth_step = max(18.0, min(56.0, depth_step))

        distance_map: Dict[int, float] = {}
        max_distance = self._collect_branch_distances(root, 0.0, distance_map)
        use_branch_lengths = self._preserve_branch_lengths and max_distance > 0.0
        branch_scale = target_width / float(max_distance) if use_branch_lengths else 0.0

        # Enforce strict vertical leaf ranking:
        # leaf_0 is above leaf_1 above leaf_2, etc.
        y_by_layout: Dict[int, float] = {}
        for idx, leaf_layout in enumerate(layout_leaves):
            y_by_layout[id(leaf_layout)] = top_margin + float(idx) * leaf_step

        def assign_internal_y(layout_node: _LayoutNode) -> float:
            node_id = id(layout_node)
            if node_id in y_by_layout:
                return y_by_layout[node_id]
            if not layout_node.children:
                y_by_layout[node_id] = top_margin
                return top_margin
            child_ys = [assign_internal_y(ch) for ch in layout_node.children]
            y = sum(child_ys) / float(len(child_ys)) if child_ys else top_margin
            y_by_layout[node_id] = y
            return y

        assign_internal_y(layout_root)

        pos_map: Dict[int, QtCore.QPointF] = {}
        for layout_node in self._iter_layout_nodes(layout_root):
            if use_branch_lengths:
                px = left_margin + float(distance_map.get(id(layout_node.node), 0.0)) * branch_scale
            else:
                px = left_margin + float(layout_node.depth) * depth_step
            py = y_by_layout[id(layout_node)]
            pos_map[id(layout_node.node)] = QtCore.QPointF(px, py)

        subtree_leaf_map: Dict[int, List[str]] = {}
        self._collect_subtree_leaf_names(root, leaf_name_map, subtree_leaf_map)

        def draw_edges(layout_node: _LayoutNode, path_stack: List[_TreeEdgeItem]) -> None:
            parent_point = pos_map[id(layout_node.node)]
            if not layout_node.children:
                return

            ordered_children = sorted(layout_node.children, key=lambda ch: pos_map[id(ch.node)].y())
            for child_layout in ordered_children:
                child_point = pos_map[id(child_layout.node)]
                curve = QtGui.QPainterPath(parent_point)
                curve.lineTo(parent_point.x(), child_point.y())
                curve.lineTo(child_point)
                subtree_leaves = subtree_leaf_map.get(id(child_layout.node), [])
                edge_item = _TreeEdgeItem(
                    curve,
                    self,
                    subtree_leaves,
                    self._edge_pen_for_color(self._edge_color_for_subtree(subtree_leaves)),
                    self._hover_edge_pen,
                )
                self._scene.addItem(edge_item)

                stack = path_stack + [edge_item]
                if child_layout.is_leaf:
                    leaf_name = leaf_name_map[id(child_layout.node)]
                    self._edge_items_by_leaf[leaf_name] = stack
                draw_edges(child_layout, stack)

        draw_edges(layout_root, [])

        label_offset = max(10.0, min(22.0, depth_step * 0.35))
        connector_gap = 3.0
        label_right_by_leaf: Dict[str, float] = {}

        for layout_leaf in layout_leaves:
            leaf = layout_leaf.node
            leaf_name = leaf_name_map[id(leaf)]
            leaf_point = pos_map[id(leaf)]
            connector_color = self._edge_color_for_subtree([leaf_name])
            self._leaf_order.append(leaf_name)
            self._leaf_pos_scene[leaf_name] = leaf_point

            label_text = self._display_leaf_name(leaf_name)
            label_item = _LeafLabelItem(label_text, leaf_name, self)
            label_item.setBrush(QtGui.QBrush(self._label_color(leaf_name)))
            label_item.setToolTip(self._leaf_tooltip(leaf_name, label_text))
            br = label_item.boundingRect()
            label_x = leaf_point.x() + label_offset
            label_item.setPos(label_x, leaf_point.y() - br.height() / 2.0)
            self._scene.addItem(label_item)
            self._leaf_label_items[leaf_name] = label_item
            label_right_by_leaf[leaf_name] = label_x + br.width()

            connector = QtGui.QPainterPath(leaf_point)
            connector.lineTo(label_x - connector_gap, leaf_point.y())
            connector_item = _TreeEdgeItem(
                connector,
                self,
                [leaf_name],
                self._edge_pen_for_color(connector_color, width=1.0),
                self._hover_edge_pen,
            )
            self._scene.addItem(connector_item)
            if leaf_name in self._edge_items_by_leaf:
                self._edge_items_by_leaf[leaf_name] = self._edge_items_by_leaf[leaf_name] + [connector_item]
            else:
                self._edge_items_by_leaf[leaf_name] = [connector_item]

        track_anchor_by_leaf = dict(label_right_by_leaf)
        content_right_by_leaf = self._draw_side_metadata_tracks(
            layout_leaves,
            pos_map,
            leaf_name_map,
            track_anchor_by_leaf,
        )

        if self._use_group_backgrounds:
            self._draw_group_backgrounds(
                layout_leaves,
                pos_map,
                leaf_name_map,
                label_right_by_leaf,
                left_margin=left_margin,
            )
        if self._show_cluster_annotations:
            self._draw_cluster_annotations(layout_leaves, pos_map, leaf_name_map, content_right_by_leaf)
        if self._bottom_legend_sections:
            content_bounds = self._scene.itemsBoundingRect()
            self._draw_bottom_legend_sections(
                left_margin,
                content_bounds.bottom() + 18.0,
                max(content_bounds.right(), max(content_right_by_leaf.values(), default=content_bounds.right())),
            )

        bounds = self._scene.itemsBoundingRect().adjusted(-20.0, -20.0, 20.0, 20.0)
        self._scene.setSceneRect(bounds)
        self.resetTransform()
        self.centerOn(bounds.center())

    def _draw_side_metadata_tracks(
        self,
        layout_leaves: Sequence[_LayoutNode],
        pos_map: Mapping[int, QtCore.QPointF],
        leaf_name_map: Mapping[int, str],
        label_right_by_leaf: Mapping[str, float],
    ) -> Dict[str, float]:
        if not self._side_metadata_tracks or not layout_leaves or not label_right_by_leaf:
            return dict(label_right_by_leaf)

        track_right_by_leaf = dict(label_right_by_leaf)
        track_gap = 16.0
        inter_track_gap = 10.0
        cell_width = 34.0
        cell_radius = 4.0
        current_left = max(label_right_by_leaf.values()) + track_gap
        layout_count = len(layout_leaves)
        total_tracks = len(self._side_metadata_tracks)
        for track in self._side_metadata_tracks:
            mode = str(track.get("mode") or "categorical").strip().lower() or "categorical"
            if mode == "heatmap_matrix":
                track_right = self._draw_matrix_side_metadata_track(
                    track,
                    layout_leaves,
                    pos_map,
                    leaf_name_map,
                    current_left,
                    total_tracks=total_tracks,
                )
                for layout_leaf in layout_leaves:
                    leaf_name = leaf_name_map.get(id(layout_leaf.node), "")
                    if leaf_name:
                        track_right_by_leaf[leaf_name] = track_right
                current_left = track_right + inter_track_gap
                continue
            track_left = current_left
            track_right = track_left + cell_width
            title_text = str(track.get("title") or "").strip()
            missing_label = str(track.get("missing_label") or "Missing").strip() or "Missing"
            value_map = dict(track.get("values", {}))
            color_map = dict(track.get("value_colors", {}))

            if title_text:
                title_font = QtGui.QFont(self.font())
                base_size = title_font.pointSizeF()
                if base_size <= 0:
                    base_size = 9.0
                if total_tracks <= 1:
                    title_size = max(34.0, base_size + 20.0)
                elif total_tracks <= 4:
                    title_size = max(25.0, base_size + 9.0)
                else:
                    title_size = max(22.0, base_size + 5.0)
                title_font.setPointSizeF(title_size)
                title_font.setBold(True)
                title_item = QtWidgets.QGraphicsSimpleTextItem(title_text)
                title_item.setFont(title_font)
                title_item.setBrush(QtGui.QBrush(QtGui.QColor("#334155")))
                title_item.setTransformOriginPoint(title_item.boundingRect().left(), title_item.boundingRect().bottom())
                title_item.setRotation(-45.0)
                title_item.setPos(track_left + 18.0, -32.0)
                self._scene.addItem(title_item)

            for idx, layout_leaf in enumerate(layout_leaves):
                leaf_name = leaf_name_map.get(id(layout_leaf.node), "")
                if not leaf_name:
                    continue
                y_center = pos_map[id(layout_leaf.node)].y()
                prev_y = pos_map[id(layout_leaves[idx - 1].node)].y() if idx > 0 else y_center
                next_y = pos_map[id(layout_leaves[idx + 1].node)].y() if idx + 1 < layout_count else y_center
                top_y = y_center - ((y_center - prev_y) * 0.5 if idx > 0 else 9.0)
                bottom_y = y_center + ((next_y - y_center) * 0.5 if idx + 1 < layout_count else 9.0)
                rect_top = top_y + 1.5
                rect_bottom = bottom_y - 1.5
                rect_height = rect_bottom - rect_top
                if rect_height <= 1.0:
                    track_right_by_leaf[leaf_name] = track_right
                    continue

                value = str(value_map.get(leaf_name, missing_label) or missing_label)
                color = QtGui.QColor(color_map.get(value, QtGui.QColor("#e2e8f0")))
                if self._active_samples is not None and leaf_name not in self._active_samples:
                    color = QtGui.QColor("#e2e8f0")
                elif value == missing_label:
                    color.setAlpha(210)

                rect = QtCore.QRectF(track_left, rect_top, cell_width, rect_height)
                path = QtGui.QPainterPath()
                path.addRoundedRect(rect, cell_radius, cell_radius)
                cell_item = QtWidgets.QGraphicsPathItem(path)
                if leaf_name in self._active_leaf_set:
                    outline = QtGui.QPen(QtGui.QColor("#f97316"))
                    outline.setWidthF(1.8)
                    cell_item.setPen(outline)
                else:
                    cell_item.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
                cell_item.setBrush(QtGui.QBrush(color))
                label = title_text or "Metadata"
                cell_item.setToolTip(f"{leaf_name}\n{label}: {value}")
                self._scene.addItem(cell_item)
                track_right_by_leaf[leaf_name] = track_right
            current_left = track_right + inter_track_gap

        return track_right_by_leaf

    def _draw_matrix_side_metadata_track(
        self,
        track: Mapping[str, object],
        layout_leaves: Sequence[_LayoutNode],
        pos_map: Mapping[int, QtCore.QPointF],
        leaf_name_map: Mapping[int, str],
        track_left: float,
        *,
        total_tracks: int,
    ) -> float:
        title_text = str(track.get("title") or "").strip() or "Panaroo genes"
        matrix = np.asarray(track.get("matrix"), dtype=np.float32)
        if matrix.ndim != 2 or matrix.shape[0] <= 0 or matrix.shape[1] <= 0:
            return track_left
        sample_order = tuple(str(sample).strip() for sample in track.get("sample_order", ()) if str(sample).strip())
        row_lookup: Dict[str, int] = {sample: idx for idx, sample in enumerate(sample_order)}
        row_spans: List[Tuple[str, float, float]] = []
        matrix_row_indices: List[int] = []
        layout_count = len(layout_leaves)
        for idx, layout_leaf in enumerate(layout_leaves):
            leaf_name = leaf_name_map.get(id(layout_leaf.node), "")
            row_idx = row_lookup.get(leaf_name)
            if row_idx is None:
                continue
            y_center = pos_map[id(layout_leaf.node)].y()
            prev_y = pos_map[id(layout_leaves[idx - 1].node)].y() if idx > 0 else y_center
            next_y = pos_map[id(layout_leaves[idx + 1].node)].y() if idx + 1 < layout_count else y_center
            top_y = y_center - ((y_center - prev_y) * 0.5 if idx > 0 else 9.0) + 1.5
            bottom_y = y_center + ((next_y - y_center) * 0.5 if idx + 1 < layout_count else 9.0) - 1.5
            if bottom_y <= top_y:
                continue
            row_spans.append((leaf_name, top_y, bottom_y))
            matrix_row_indices.append(row_idx)
        if not row_spans:
            return track_left
        ordered_matrix = matrix[np.asarray(matrix_row_indices, dtype=int), :]

        if title_text:
            title_font = QtGui.QFont(self.font())
            base_size = title_font.pointSizeF()
            if base_size <= 0:
                base_size = 9.0
            if total_tracks <= 1:
                title_size = max(34.0, base_size + 20.0)
            elif total_tracks <= 4:
                title_size = max(25.0, base_size + 9.0)
            else:
                title_size = max(22.0, base_size + 5.0)
            title_font.setPointSizeF(title_size)
            title_font.setBold(True)
            title_item = QtWidgets.QGraphicsSimpleTextItem(title_text)
            title_item.setFont(title_font)
            title_item.setBrush(QtGui.QBrush(QtGui.QColor("#334155")))
            title_item.setTransformOriginPoint(title_item.boundingRect().left(), title_item.boundingRect().bottom())
            title_item.setRotation(-45.0)
            title_item.setPos(track_left + 18.0, -32.0)
            self._scene.addItem(title_item)

        n_cols = matrix.shape[1]
        strip_width = min(1600.0, max(320.0, float(n_cols) * 2.0))
        top = min(top for _sample, top, _bottom in row_spans)
        bottom = max(bottom for _sample, _top, bottom in row_spans)
        item = _MatrixStripItem(
            QtCore.QRectF(track_left, top, strip_width, bottom - top),
            title=title_text,
            rows=row_spans,
            matrix=ordered_matrix,
            column_labels=tuple(track.get("column_labels", ())),
            column_sizes=tuple(track.get("column_sizes", ())),
            present_color=self._coerce_color(track.get("present_color")) or QtGui.QColor("#2563eb"),
            absent_color=self._coerce_color(track.get("absent_color")) or QtGui.QColor("#e2e8f0"),
            binned=bool(track.get("binned")),
            active_samples=self._active_samples,
        )
        item.set_highlighted_samples(self._active_leaf_set)
        self._scene.addItem(item)
        self._matrix_track_items.append(item)
        return track_left + strip_width

    def _draw_bottom_legend_sections(
        self,
        left_x: float,
        top_y: float,
        right_x: float,
    ) -> None:
        if not self._bottom_legend_sections:
            return
        items_per_column = 3
        swatch_size = 36.0
        label_gap = 14.0
        column_gap = 36.0
        row_gap = 18.0
        section_gap = 64.0
        x_cursor = float(left_x)
        top = float(top_y)

        title_font = QtGui.QFont(self.font())
        title_font.setBold(True)
        title_font.setPointSizeF(max(28.0, title_font.pointSizeF() * 2.6))
        text_font = QtGui.QFont(self.font())
        text_font.setPointSizeF(max(20.0, text_font.pointSizeF() * 2.2))
        title_metrics = QtGui.QFontMetricsF(title_font)
        text_metrics = QtGui.QFontMetricsF(text_font)

        for section in self._bottom_legend_sections:
            title = str(section.get("title") or "").strip()
            items = list(section.get("items", []))
            if not title or not items:
                continue

            title_item = QtWidgets.QGraphicsSimpleTextItem(title)
            title_item.setFont(title_font)
            title_item.setBrush(QtGui.QBrush(QtGui.QColor("#334155")))
            title_item.setPos(x_cursor, top)
            self._scene.addItem(title_item)
            section_content_y = top + title_metrics.height() + 14.0

            legend_items: List[Tuple[str, QtGui.QColor, float, float]] = []
            for label, color in items:
                label_text = str(label or "").strip()
                if not label_text:
                    continue
                text_width = text_metrics.horizontalAdvance(label_text)
                item_width = swatch_size + label_gap + text_width
                legend_items.append((label_text, QtGui.QColor(color), text_width, item_width))
            if not legend_items:
                continue

            column_count = max(1, int(math.ceil(len(legend_items) / float(items_per_column))))
            column_widths: List[float] = [0.0] * column_count
            for idx, (_label_text, _color, _text_width, item_width) in enumerate(legend_items):
                col = idx // items_per_column
                column_widths[col] = max(column_widths[col], item_width)

            column_x: List[float] = []
            section_x = float(x_cursor)
            for width in column_widths:
                column_x.append(section_x)
                section_x += width + column_gap

            row_height = max(text_metrics.height(), swatch_size) + row_gap

            for idx, (label_text, color, _text_width, _item_width) in enumerate(legend_items):
                row = idx % items_per_column
                col = idx // items_per_column
                x = column_x[col]
                item_y = section_content_y + row * row_height

                swatch_rect = QtCore.QRectF(x, item_y + 4.0, swatch_size, swatch_size)
                swatch_item = QtWidgets.QGraphicsRectItem(swatch_rect)
                swatch_pen = QtGui.QPen(QtGui.QColor("#94a3b8"))
                swatch_pen.setWidthF(2.0)
                swatch_item.setPen(swatch_pen)
                swatch_item.setBrush(QtGui.QBrush(QtGui.QColor(color)))
                self._scene.addItem(swatch_item)

                text_item = QtWidgets.QGraphicsSimpleTextItem(label_text)
                text_item.setFont(text_font)
                text_item.setBrush(QtGui.QBrush(QtGui.QColor("#334155")))
                text_item.setPos(x + swatch_size + label_gap, item_y)
                self._scene.addItem(text_item)

            section_width = max(
                title_metrics.horizontalAdvance(title),
                sum(column_widths) + max(0, column_count - 1) * column_gap,
            )
            x_cursor += section_width + section_gap

    def _draw_cluster_annotations(
        self,
        layout_leaves: Sequence[_LayoutNode],
        pos_map: Mapping[int, QtCore.QPointF],
        leaf_name_map: Mapping[int, str],
        label_right_by_leaf: Mapping[str, float],
    ) -> None:
        if not self._leaf_cluster_map or not layout_leaves:
            return

        cluster_runs: List[Tuple[str, List[str], float, float]] = []
        cluster_leaves_by_name: Dict[str, List[str]] = {}
        current_cluster: Optional[str] = None
        current_leaves: List[str] = []
        start_y = 0.0
        end_y = 0.0

        def flush_run() -> None:
            nonlocal current_cluster, current_leaves, start_y, end_y
            if current_cluster and current_leaves:
                cluster_runs.append((current_cluster, list(current_leaves), start_y, end_y))
            current_cluster = None
            current_leaves = []
            start_y = 0.0
            end_y = 0.0

        for layout_leaf in layout_leaves:
            leaf_name = leaf_name_map.get(id(layout_leaf.node), "")
            cluster_name = self._leaf_cluster(leaf_name)
            if not cluster_name:
                flush_run()
                continue
            leaf_y = pos_map[id(layout_leaf.node)].y()
            if cluster_name != current_cluster:
                flush_run()
                current_cluster = cluster_name
                current_leaves = [leaf_name]
                cluster_leaves_by_name.setdefault(cluster_name, []).append(leaf_name)
                start_y = leaf_y
                end_y = leaf_y
                continue
            current_leaves.append(leaf_name)
            cluster_leaves_by_name.setdefault(cluster_name, []).append(leaf_name)
            end_y = leaf_y
        flush_run()

        if not cluster_runs:
            return

        font = QtGui.QFont(self.font())
        base_size = font.pointSizeF()
        if base_size <= 0:
            base_size = 9.0
        font.setPointSizeF(max(16.0, min(24.0, base_size + 5.0)))
        font.setBold(True)
        metrics = QtGui.QFontMetrics(font)
        bracket_gap = 7.0
        bracket_width = 10.0
        label_gap = 8.0

        for cluster_name, cluster_leaves, raw_top, raw_bottom in cluster_runs:
            linked_cluster_leaves = cluster_leaves_by_name.get(cluster_name, cluster_leaves)
            top_y = raw_top
            bottom_y = raw_bottom
            if math.isclose(top_y, bottom_y):
                top_y -= 6.0
                bottom_y += 6.0

            dominant_group = self._subtree_group(cluster_leaves)
            if self._use_group_backgrounds:
                base_color = QtGui.QColor("#000000")
            else:
                base_color = self._group_color(dominant_group) if dominant_group else QtGui.QColor("#475569")
            color = self._visual_color_for_subtree(cluster_leaves, base_color)
            pen = self._edge_pen_for_color(color, width=1.6)

            cluster_right = max(label_right_by_leaf.get(leaf_name, 0.0) for leaf_name in cluster_leaves)
            inner_x = cluster_right + bracket_gap
            outer_x = inner_x + bracket_width

            bracket = QtGui.QPainterPath(QtCore.QPointF(inner_x, top_y))
            bracket.lineTo(outer_x, top_y)
            bracket.lineTo(outer_x, bottom_y)
            bracket.lineTo(inner_x, bottom_y)
            cluster_payload = {"name": cluster_name, "samples": list(linked_cluster_leaves)}
            bracket_item = QtWidgets.QGraphicsPathItem(bracket)
            bracket_item.setPen(pen)
            self._scene.addItem(bracket_item)

            label_text = f"Cluster {cluster_name}"
            label_item = _ActionTextItem(
                label_text,
                callback=lambda payload=cluster_payload: self._on_cluster_click(payload),
                enable_hover=self._enable_hover_interactions,
                hover_enter_callback=lambda leaves=tuple(linked_cluster_leaves): self._on_action_hover_enter(leaves),
                hover_leave_callback=lambda leaves=tuple(linked_cluster_leaves): self._on_action_hover_leave(leaves),
            )
            label_item.setFont(font)
            label_item.setBrush(QtGui.QBrush(color))
            label_item.setToolTip(label_text)
            label_y = ((top_y + bottom_y) * 0.5) - (metrics.height() * 0.5)
            label_item.setPos(outer_x + label_gap, label_y)
            self._scene.addItem(label_item)

    def _draw_group_backgrounds(
        self,
        layout_leaves: Sequence[_LayoutNode],
        pos_map: Mapping[int, QtCore.QPointF],
        leaf_name_map: Mapping[int, str],
        label_right_by_leaf: Mapping[str, float],
        *,
        left_margin: float,
    ) -> None:
        if not self._leaf_group_map or not layout_leaves:
            return

        group_runs: List[Tuple[str, List[str], float, float]] = []
        current_group: Optional[str] = None
        current_leaves: List[str] = []
        start_y = 0.0
        end_y = 0.0

        def flush_run() -> None:
            nonlocal current_group, current_leaves, start_y, end_y
            if current_group and current_leaves:
                group_runs.append((current_group, list(current_leaves), start_y, end_y))
            current_group = None
            current_leaves = []
            start_y = 0.0
            end_y = 0.0

        for layout_leaf in layout_leaves:
            leaf_name = leaf_name_map.get(id(layout_leaf.node), "")
            group_name = self._leaf_group(leaf_name)
            if not group_name:
                flush_run()
                continue
            leaf_y = pos_map[id(layout_leaf.node)].y()
            if group_name != current_group:
                flush_run()
                current_group = group_name
                current_leaves = [leaf_name]
                start_y = leaf_y
                end_y = leaf_y
                continue
            current_leaves.append(leaf_name)
            end_y = leaf_y
        flush_run()

        if not group_runs:
            return

        clade_left_shift = 18.0
        left_x = max(0.0, left_margin - (16.0 + clade_left_shift))
        row_pad = 12.0
        min_height = 18.0
        radius = 10.0
        text_right_pad = 40.0
        text_left_pad = 40.0

        for group_name, group_leaves, raw_top, raw_bottom in group_runs:
            top_y = raw_top - row_pad
            bottom_y = raw_bottom + row_pad
            if (bottom_y - top_y) < min_height:
                center_y = (top_y + bottom_y) * 0.5
                top_y = center_y - (min_height * 0.5)
                bottom_y = center_y + (min_height * 0.5)

            rect_height = max(min_height, bottom_y - top_y)
            label_text = f"Clade {group_name}"
            text_font = QtGui.QFont(self.font())
            base_size = text_font.pointSizeF()
            if base_size <= 0:
                base_size = 9.0
            text_font.setPointSizeF(max(20.0, min(40.0, base_size + max(10.0, rect_height * 0.16))))
            text_font.setBold(True)
            text_metrics = QtGui.QFontMetrics(text_font)

            base_right_x = max(label_right_by_leaf.get(leaf_name, left_x + 40.0) for leaf_name in group_leaves)
            if self._show_group_labels:
                extra_right = max(
                    285.0,
                    float(text_metrics.horizontalAdvance(label_text)) + text_left_pad + text_right_pad + 110.0,
                )
                right_x = base_right_x + extra_right
            else:
                right_x = base_right_x + 10.0
            rect = QtCore.QRectF(
                left_x,
                top_y,
                max(24.0, right_x - left_x),
                rect_height,
            )
            path = QtGui.QPainterPath()
            path.addRoundedRect(rect, radius, radius)
            bg_item = QtWidgets.QGraphicsPathItem(path)
            bg_item.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
            bg_item.setBrush(self._group_background_brush(group_name, group_leaves))
            bg_item.setZValue(-10.0)
            self._scene.addItem(bg_item)

            if self._show_group_labels:
                clade_payload = {"name": group_name, "samples": list(group_leaves)}
                text_item = _ActionTextItem(
                    label_text,
                    callback=lambda payload=clade_payload: self._on_clade_click(payload),
                    enable_hover=self._enable_hover_interactions,
                    hover_enter_callback=lambda leaves=tuple(group_leaves): self._on_action_hover_enter(leaves),
                    hover_leave_callback=lambda leaves=tuple(group_leaves): self._on_action_hover_leave(leaves),
                )
                text_item.setFont(text_font)
                text_color = self._visual_color_for_subtree(group_leaves, QtGui.QColor("#000000"))
                text_item.setBrush(QtGui.QBrush(text_color))
                text_item.setToolTip(label_text)
                text_width = float(text_metrics.horizontalAdvance(label_text))
                text_height = float(text_metrics.height())
                text_x = rect.right() - text_right_pad - text_width - clade_left_shift
                text_y = rect.center().y() - (text_height * 0.5)
                text_item.setPos(text_x, text_y)
                self._scene.addItem(text_item)

    def _iter_internal_nodes(self, root: _TreeNode) -> List[_TreeNode]:
        nodes: List[_TreeNode] = []
        stack: List[Tuple[_TreeNode, bool]] = [(root, True)]
        while stack:
            node, is_root = stack.pop()
            if not node.is_leaf and not is_root:
                nodes.append(node)
            for child in reversed(node.children):
                stack.append((child, False))
        return nodes

    def _label_color(self, leaf_name: str) -> QtGui.QColor:
        base = self._leaf_base_color(leaf_name)
        if self._action_hover_key is not None and leaf_name in self._active_leaf_set:
            return QtGui.QColor("#ea580c")
        if self._active_samples is None or leaf_name in self._active_samples:
            return base
        return QtGui.QColor("#94a3b8")

    def _label_brush(self, color: QtGui.QColor) -> QtGui.QBrush:
        key = (int(color.red()), int(color.green()), int(color.blue()), int(color.alpha()))
        brush = self._label_brush_cache.get(key)
        if brush is None:
            brush = QtGui.QBrush(QtGui.QColor(key[0], key[1], key[2], key[3]))
            self._label_brush_cache[key] = brush
        return brush

    def _ordered_leaf_subset(self, leaves: Sequence[str]) -> List[str]:
        wanted = {self._normalize_label(x) for x in leaves if str(x).strip()}
        if not wanted:
            return []
        return [name for name in self._leaf_order if name in wanted]

    def _set_active_leaves(self, leaves: Optional[Sequence[str]], emit_signals: bool = True) -> None:
        ordered = self._ordered_leaf_subset(leaves or [])
        next_set = set(ordered)
        if next_set == self._active_leaf_set:
            return

        if self._active_edge_items:
            for item in self._active_edge_items:
                item.set_hovered(False)
        if self._active_marker_items:
            for item in self._active_marker_items:
                item.set_hovered(False)
        for old_leaf in self._active_leaf_set:
            old_label = self._leaf_label_items.get(old_leaf)
            if old_label is not None:
                old_label.setBrush(self._label_brush(self._label_color(old_leaf)))

        active_items: Set[_TreeEdgeItem] = set()
        active_markers: Set[_NodeMarkerItem] = set()
        for leaf in ordered:
            for item in self._edge_items_by_leaf.get(leaf, []):
                active_items.add(item)
            for marker in self._marker_items_by_leaf.get(leaf, []):
                active_markers.add(marker)
            label_item = self._leaf_label_items.get(leaf)
            if label_item is not None:
                label_item.setBrush(self._label_brush(QtGui.QColor("#ea580c")))

        for item in active_items:
            item.set_hovered(True)
        for item in active_markers:
            item.set_hovered(True)

        self._active_edge_items = active_items
        self._active_marker_items = active_markers
        self._active_leaf_set = next_set
        self._hover_leaf = ordered[0] if len(ordered) == 1 else None
        for item in self._matrix_track_items:
            item.set_highlighted_samples(self._active_leaf_set)

        if emit_signals:
            if len(ordered) == 1:
                self.leafHovered.emit(ordered[0])
            else:
                self.leafHovered.emit(None)
            if len(ordered) >= 2:
                self.subtreeHovered.emit(ordered)
            elif not ordered:
                self.subtreeHovered.emit(None)

    def _on_leaf_hover_enter(self, leaf_name: str) -> None:
        key = (leaf_name,)
        self._hover_subtree_key = key
        self._set_active_leaves([leaf_name], emit_signals=True)

    def _on_leaf_double_click(self, leaf_name: str) -> None:
        self.leafDoubleClicked.emit(leaf_name)

    def _on_clade_click(self, payload: Mapping[str, object]) -> None:
        self.cladeClicked.emit(dict(payload))

    def _on_cluster_click(self, payload: Mapping[str, object]) -> None:
        self.clusterClicked.emit(dict(payload))

    def _on_action_hover_enter(self, leaves: Sequence[str]) -> None:
        ordered = self._ordered_leaf_subset(leaves)
        key = tuple(ordered)
        self._action_hover_key = key if key else None
        self._hover_subtree_key = key if key else self._hover_subtree_key
        self._set_active_leaves(ordered, emit_signals=True)

    def _on_action_hover_leave(self, leaves: Sequence[str]) -> None:
        ordered = self._ordered_leaf_subset(leaves)
        key = tuple(ordered)
        if self._action_hover_key != key:
            return
        self._action_hover_key = None
        if self._hover_subtree_key == key:
            self._hover_subtree_key = None
        self._set_active_leaves(None, emit_signals=True)

    def _on_leaf_hover_leave(self, leaf_name: str) -> None:
        key = (leaf_name,)
        if self._hover_subtree_key == key:
            self._hover_subtree_key = None
            self._set_active_leaves(None, emit_signals=True)

    def _on_subtree_hover_enter(self, leaves: Sequence[str]) -> None:
        ordered = self._ordered_leaf_subset(leaves)
        key = tuple(ordered)
        self._hover_subtree_key = key
        self._set_active_leaves(ordered, emit_signals=True)

    def _on_subtree_hover_leave(self, leaves: Sequence[str]) -> None:
        ordered = self._ordered_leaf_subset(leaves)
        key = tuple(ordered)
        if self._hover_subtree_key == key:
            self._hover_subtree_key = None
            self._set_active_leaves(None, emit_signals=True)

    def leaveEvent(self, event: QtCore.QEvent) -> None:
        self._hover_subtree_key = None
        self._set_active_leaves(None, emit_signals=True)
        super().leaveEvent(event)

    def _fit_to_viewport(self) -> None:
        """Scale the whole tree so it fits the current viewport (keeping aspect)."""
        if self._root is None:
            return
        viewport = self.viewport()
        if viewport is None or viewport.width() <= 0 or viewport.height() <= 0:
            return
        rect = self._scene.itemsBoundingRect()
        if rect.isNull() or rect.isEmpty():
            rect = self.sceneRect()
        if rect.isNull() or rect.isEmpty():
            return
        # Stretch to fill the window in both dimensions (no aspect lock).
        self.fitInView(rect, QtCore.Qt.AspectRatioMode.IgnoreAspectRatio)

    def fit_to_window(self) -> None:
        """Public: re-enable auto-fit and scale the tree to the window."""
        self._auto_fit = True
        self._fit_to_viewport()

    def wheelEvent(self, event: QtGui.QWheelEvent) -> None:
        delta = event.angleDelta().y()
        if delta == 0:
            super().wheelEvent(event)
            return
        factor = 1.15 if delta > 0 else (1.0 / 1.15)
        current = self.transform().m11()
        new_scale = current * factor
        if new_scale < 0.03 or new_scale > 120.0:
            event.accept()
            return
        self._auto_fit = False  # manual zoom takes over until the tree reloads
        self.scale(factor, factor)
        self.viewChanged.emit()
        event.accept()

    def mouseDoubleClickEvent(self, event: QtGui.QMouseEvent) -> None:
        # Double-click on empty background re-fits the tree to the window.
        if self.itemAt(event.pos()) is None:
            self.fit_to_window()
            event.accept()
            return
        super().mouseDoubleClickEvent(event)

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        super().resizeEvent(event)
        if self._auto_fit:
            self._fit_to_viewport()
        self.viewChanged.emit()

    def scrollContentsBy(self, dx: int, dy: int) -> None:
        super().scrollContentsBy(dx, dy)
        self.viewChanged.emit()
