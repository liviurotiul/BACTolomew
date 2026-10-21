from __future__ import annotations

from contextlib import contextmanager
import math
from dataclasses import dataclass, field
from decimal import Decimal
from pathlib import Path
from typing import Callable, Dict, Iterator, List, Mapping, Optional, Sequence, Set, Tuple

from PyQt6 import QtCore, QtGui, QtWidgets

from .phylo import (
    ORDER_AS_COMPUTED,
    ORDER_MODES,
    ROOT_AS_IN_FILE,
    ROOT_GROUPS,
    ROOT_MIDPOINT,
    ROOT_MODES,
    BranchBreaks,
    GroupOptimizer,
    GroupSummary,
    NewickParser as _NewickParser,
    TreeNode as _TreeNode,
    analyze_branch_lengths,
    clone_tree as _clone_tree,
    collapse_root_grade,
    edge_length,
    follow_order,
    full_subtree,
    group_root,
    group_starts,
    groups_contiguous,
    ladderize,
    leaf_names as _leaf_names,
    long_branch_breaks,
    midpoint_root as _midpoint_root,
    prune_to as _prune_to,
    resolve_collapsed,
    root_grades,
)
from .row_range_note import RowRangeNote, row_range_text, visible_rows
from .theme import COLORS, out_of_set_color


_TREE_CANVAS_COLOR = COLORS["surface"]
_TREE_LABEL_COLOR = COLORS["text"]
_TREE_BRANCH_COLOR = COLORS["tree_branch"]
_TREE_INACTIVE_COLOR = COLORS["tree_inactive"]
# Hover emphasis (the shared amber of every panel's hover): the hovered
# root-to-tip path, the tip dot and the name chip; the hovered labels use
# the darker ink so they stay readable.
_TREE_HOVER_COLOR = COLORS["hover"]
_TREE_HOVER_INK = COLORS["hover_ink"]
# The hovered path's width in device pixels (a cosmetic pen: the same at
# every zoom) and its z above every other branch.
_HOVER_PATH_WIDTH_PX = 3.0
_HOVER_Z = 50.0
# While a hover is active the rest of the tree is drawn at this opacity.
_HOVER_DIM_OPACITY = 0.35
# The tip dot's diameter and the name chips (device pixels), and how many
# chips a multi-sample hover shows before a single "+N more" chip.
_HOVER_DOT_DIAMETER_PX = 7.0
_HOVER_CHIP_CAP = 3
_HOVER_CHIP_GAP_PX = 8.0
_HOVER_CHIP_MARGIN_PX = 4.0
_TREE_CLUSTER_COLOR = "#7a8796"
_TREE_REFERENCE_COLOR = "#0f766e"
_CLADE_TINT_ALPHA = 16
_CLADE_PARTIAL_TINT_ALPHA = 12
_CLADE_INACTIVE_TINT_ALPHA = 10
_CLADE_RAIL_WIDTH = 4.0
_LEAF_FONT_SIZE = 9.5
_CLUSTER_FONT_SIZE = 10.5
_CLADE_FONT_SIZE = 11.5
# The long branch rule of the Distance drawing (phylo.long_branch_breaks): an
# edge longer than this many times the 95th percentile root to tip depth of
# the rest of the tree is drawn short, with a break mark and its true length.
LONG_BRANCH_BREAK_FACTOR = 10.0
# The break mark (two slanted strokes across the edge) and the true length
# label above it, in device pixels.
_BREAK_MARK_HALF_HEIGHT_PX = 5.0
_BREAK_MARK_GAP_PX = 3.0
_BREAK_LABEL_FONT_SIZE = 8.0
# A collapsed subtree takes _collapsed_row_count rows (2 for 2 genomes, 3 up
# to 10, then a little more for bigger ones, at most _COLLAPSED_MAX_ROWS, so
# collapsing really shortens the tree and size still shows); its triangle
# leaves _COLLAPSED_ROW_GAP rows free and is at least _COLLAPSED_MIN_PX tall
# on screen (within its rows); it is filled with its color at
# _COLLAPSED_FILL_ALPHA of its alpha, and labeled in _COLLAPSED_LABEL_FONT_SIZE
# (device pixels: readable in a tree fitted to a short panel).
_COLLAPSED_MAX_ROWS = 5.0
_COLLAPSED_ROW_GAP = 0.6
_COLLAPSED_MIN_PX = 16.0
_COLLAPSED_LABEL_FONT_SIZE = 8.5
_COLLAPSED_FILL_ALPHA = 0.28
# The collapse handle (device pixels): its radius at rest, collapsed and
# hovered, and the radius that takes the clicks.
_HANDLE_RADIUS_PX = 3.5
_HANDLE_COLLAPSED_RADIUS_PX = 5.0
_HANDLE_HOVER_RADIUS_PX = 6.5
_HANDLE_HIT_RADIUS_PX = 7.5
# A clade or cluster in several subtrees (not monophyletic) gets a handle on
# its _SPLIT_GROUP_MAX_HANDLES biggest subtrees of _SPLIT_GROUP_MIN_HANDLE_SIZE
# genomes or more (a right click collapses any other).
_SPLIT_GROUP_MAX_HANDLES = 6
_SPLIT_GROUP_MIN_HANDLE_SIZE = 3
# A right click this close to a node (viewport pixels) offers to collapse it.
_NODE_MENU_RADIUS_PX = 12.0
# Navigation. One zoom step (a wheel notch, a toolbar button, +/-) and the
# trackpad pixels that make one; the camera zoom's lower limit (of the fit)
# and the leaf labels' largest size on screen (times their font size); the
# screen stable annotations grow with the camera up to _STABLE_ZOOM_MAX;
# a stretch stays within _STRETCH_MIN to _STRETCH_MAX of the fit on each
# axis; panning always keeps _PAN_KEEP_PX of the drawing in view.
_ZOOM_STEP = 1.2
_PIXELS_PER_ZOOM_STEP = 60.0
_CAMERA_ZOOM_MIN = 0.5
_LABEL_MAX_SCALE = 4.0
_STABLE_ZOOM_MAX = 2.5
_STRETCH_MIN = 0.25
_STRETCH_MAX = 40.0
_PAN_KEEP_PX = 40.0
# The scene units of air around the drawing in the scene rect and the fit.
_SCENE_RECT_PAD = 20.0

_SCENE_DATA_KIND = 0
_SCENE_DATA_NAME = 1
_SCENE_DATA_RUN_INDEX = 2
_SCENE_DATA_RUN_COUNT = 3
_SCENE_DATA_VALUE = 4
# The leaves of a clade or cluster run (its chip's priority when chips overlap).
_SCENE_DATA_RUN_SIZE = 5

# Root modes (ROOT_*, from phylo): ROOT_GROUPS reroots on the edge
# that keeps the groups together, as ``reroot_for_groups`` does.


_IGNORES_PARENT_OPACITY = QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIgnoresParentOpacity


def _node_support(node: _TreeNode) -> Optional[float]:
    """A node's branch support: its ``support`` when the parser gives one,
    else its label when that is a number (IQ-TREE writes the support of an
    edge as the label of the node below it). None for a leaf or a name."""
    if node.is_leaf:
        return None
    value = getattr(node, "support", None)
    if value is None:
        try:
            value = float(str(node.name or "").strip())
        except ValueError:
            return None
    value = float(value)
    return value if math.isfinite(value) else None


def _support_by_split(root: _TreeNode, leaf_name: Callable[[str], str]) -> Tuple[Dict[frozenset, float], frozenset]:
    """The support of each edge of the file tree by the split it makes (the
    side without the first leaf, so any rooting of the tree finds it), and
    every leaf. A rerooted or rotated drawing looks its edges up here, so the
    values stay with their edges whatever root is drawn."""
    leaves_of: Dict[int, frozenset] = {}
    order: List[_TreeNode] = []
    stack = [root]
    while stack:
        node = stack.pop()
        order.append(node)
        stack.extend(node.children)
    for node in reversed(order):
        if node.is_leaf:
            leaves_of[id(node)] = frozenset((leaf_name(node.name or ""),))
        else:
            leaves_of[id(node)] = frozenset().union(*(leaves_of[id(child)] for child in node.children))
    every = leaves_of[id(root)]
    if not every:
        return {}, every
    anchor = min(every)
    supports: Dict[frozenset, float] = {}
    for node in order:
        if node is root:
            continue
        value = _node_support(node)
        if value is None:
            continue
        side = leaves_of[id(node)]
        key = every - side if anchor in side else side
        if 1 < len(key) < len(every) - 1:
            supports.setdefault(key, value)
    return supports, every


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
        menu_callback: Optional[Callable[[QtCore.QPoint], None]] = None,
        parent: Optional[QtWidgets.QGraphicsItem] = None,
    ) -> None:
        super().__init__(text, parent)
        self._callback = callback
        self._hover_enter_callback = hover_enter_callback
        self._hover_leave_callback = hover_leave_callback
        # Right click: a small menu of what the label can do (screen position).
        self._menu_callback = menu_callback
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
        # The shared hover ink, not primary blue (cluster 1's color).
        hover_brush = QtGui.QBrush(QtGui.QColor(_TREE_HOVER_INK))
        self.setBrush(hover_brush)
        self.setOpacity(1.0)
        # The label under the pointer is not dimmed with the rest of the tree.
        self.setFlag(_IGNORES_PARENT_OPACITY, True)
        super().hoverEnterEvent(event)

    def hoverLeaveEvent(self, event: Optional[QtWidgets.QGraphicsSceneHoverEvent]) -> None:
        if self._hover_leave_callback is not None:
            self._hover_leave_callback()
        if self._base_font is not None:
            self.setFont(QtGui.QFont(self._base_font))
        if self._base_brush is not None:
            self.setBrush(QtGui.QBrush(self._base_brush))
        self.setOpacity(1.0)
        self.setFlag(_IGNORES_PARENT_OPACITY, False)
        super().hoverLeaveEvent(event)

    def mousePressEvent(self, event: Optional[QtWidgets.QGraphicsSceneMouseEvent]) -> None:
        if event is not None and event.button() == QtCore.Qt.MouseButton.LeftButton:
            self._callback()
            event.accept()
            return
        super().mousePressEvent(event)

    def contextMenuEvent(self, event: Optional[QtWidgets.QGraphicsSceneContextMenuEvent]) -> None:
        if event is not None and self._menu_callback is not None:
            self._menu_callback(event.screenPos())
            event.accept()
            return
        super().contextMenuEvent(event)


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
        # Branch support only in the hover text, never drawn; a broken
        # edge says its true length first.
        support = owner.support_text(self._subtree_leaves) if len(self._subtree_leaves) > 1 else ""
        note = owner.branch_break_note(self._subtree_leaves)
        text = " ".join(part for part in (note, support) if part)
        if text:
            self.setToolTip(text)

    def subtree_leaves(self) -> Tuple[str, ...]:
        return self._subtree_leaves

    def paint(self, painter, option, widget=None) -> None:
        # On screen, with the rows squeezed more than the width (a tall tree
        # fitted to its panel), a scene width pen would draw the horizontal
        # branches far thinner than the vertical ones: draw a thin pen of
        # the same pixels both ways. Exports (no widget) keep the scene pen.
        pen = self.pen()
        if widget is not None and not pen.isCosmetic():
            transform = painter.worldTransform()
            scale_x = math.hypot(transform.m11(), transform.m12())
            scale_y = math.hypot(transform.m21(), transform.m22())
            if abs(scale_x - scale_y) > 1e-3 * max(scale_x, scale_y):
                device_pen = QtGui.QPen(pen)
                device_pen.setCosmetic(True)
                device_pen.setWidthF(max(0.8, pen.widthF() * min(1.0, max(scale_x, scale_y))))
                painter.setPen(device_pen)
                painter.setBrush(self.brush())
                painter.drawPath(self.path())
                return
        super().paint(painter, option, widget)

    def set_hovered(self, active: bool) -> None:
        try:
            self.setPen(self._hover_pen if active else self._base_pen)
            # Above every other branch, and not dimmed with the rest.
            self.setZValue(_HOVER_Z if active else 0.0)
            self.setFlag(_IGNORES_PARENT_OPACITY, active)
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


def _collapsed_row_count(count: int) -> float:
    """The rows a collapsed subtree of ``count`` drawn genomes takes."""
    if count <= 10:
        return float(max(1, min(count, 3)))
    return min(_COLLAPSED_MAX_ROWS, 3.0 + 1.5 * math.log10(count / 10.0))


class _CollapsedCladeItem(QtWidgets.QGraphicsPathItem):
    """A collapsed subtree: a triangle from its node to the depth of its
    deepest leaf, colored like its branches (its label, "Cluster 4 · 83
    genomes", is a _CollapsedLabelItem beside it). Hover highlights its
    genomes (in other panels too); double-click expands it; right click
    offers Expand.

    ``band`` (rectangular layout): (apex, far x, half height, half of its
    rows): the triangle is drawn at least _COLLAPSED_MIN_PX tall on screen,
    never beyond its rows; else ``path`` as it is (a wedge of the circle)."""

    def __init__(
        self,
        path: QtGui.QPainterPath,
        owner: "RaxmlTreeView",
        members: Sequence[str],
        color: QtGui.QColor,
        band: Optional[Tuple[QtCore.QPointF, float, float, float]] = None,
    ) -> None:
        super().__init__(path)
        self._owner = owner
        self._members = tuple(members)
        self._band = band
        self._drawn_half = band[2] if band is not None else 0.0
        fill = QtGui.QColor(color)
        fill.setAlphaF(fill.alphaF() * _COLLAPSED_FILL_ALPHA)
        self._base_pen = owner._edge_pen_for_color(color)
        self._base_brush = QtGui.QBrush(fill)
        hover_fill = QtGui.QColor(_TREE_HOVER_COLOR)
        hover_fill.setAlphaF(_COLLAPSED_FILL_ALPHA)
        self._hover_brush = QtGui.QBrush(hover_fill)
        self._hover_pen = QtGui.QPen(owner._hover_edge_pen)
        self.setPen(self._base_pen)
        self.setBrush(self._base_brush)
        self.setData(_SCENE_DATA_KIND, "collapsed-clade")
        self.setAcceptHoverEvents(owner.hover_interactions_enabled())

    def _triangle(self, half: float) -> QtGui.QPainterPath:
        apex, far_x, _half, _band = self._band  # type: ignore[misc]
        path = QtGui.QPainterPath(apex)
        path.lineTo(far_x, apex.y() - half)
        path.lineTo(far_x, apex.y() + half)
        path.closeSubpath()
        return path

    def boundingRect(self) -> QtCore.QRectF:
        if self._band is None:
            return super().boundingRect()
        apex, far_x, _half, band_half = self._band
        pad = self.pen().widthF() + 1.0
        return QtCore.QRectF(apex.x(), apex.y() - band_half, far_x - apex.x(), 2.0 * band_half).adjusted(
            -pad, -pad, pad, pad
        )

    def shape(self) -> QtGui.QPainterPath:
        if self._band is None:
            return super().shape()
        return self._triangle(self._drawn_half)

    def paint(self, painter, option, widget=None) -> None:
        if self._band is None:
            super().paint(painter, option, widget)
            return
        _apex, _far_x, half, band_half = self._band
        scale_y = abs(float(painter.worldTransform().m22())) or 1.0
        # On screen the least height grows with the camera like the labels.
        least = _COLLAPSED_MIN_PX * (self._owner._stable_scale() if widget is not None else 1.0)
        drawn = min(band_half, max(half, least / 2.0 / scale_y))
        # The hover shape follows what is drawn (within the bounding rect).
        self._drawn_half = drawn
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        painter.setPen(self.pen())
        painter.setBrush(self.brush())
        painter.drawPath(self._triangle(drawn))

    def members(self) -> Tuple[str, ...]:
        return self._members

    def set_hovered(self, active: bool) -> None:
        try:
            self.setPen(self._hover_pen if active else self._base_pen)
            self.setBrush(self._hover_brush if active else self._base_brush)
            self.setZValue(_HOVER_Z if active else 0.0)
            self.setFlag(_IGNORES_PARENT_OPACITY, active)
        except RuntimeError:
            return

    def hoverEnterEvent(self, event: Optional[QtWidgets.QGraphicsSceneHoverEvent]) -> None:
        self._owner._on_subtree_hover_enter(self._members)
        super().hoverEnterEvent(event)

    def hoverLeaveEvent(self, event: Optional[QtWidgets.QGraphicsSceneHoverEvent]) -> None:
        self._owner._on_subtree_hover_leave(self._members)
        super().hoverLeaveEvent(event)

    def mouseDoubleClickEvent(self, event: Optional[QtWidgets.QGraphicsSceneMouseEvent]) -> None:
        # Redrawn later: the scene must not delete this item in its own event.
        members = self._members
        QtCore.QTimer.singleShot(0, lambda: self._owner.expand_subtree(members))
        if event is not None:
            event.accept()

    def contextMenuEvent(self, event: Optional[QtWidgets.QGraphicsSceneContextMenuEvent]) -> None:
        if event is not None:
            self._owner._show_node_menu(self._members, event.screenPos())
            event.accept()


class _CollapsedLabelItem(QtWidgets.QGraphicsSimpleTextItem):
    """A collapsed subtree's label, in device pixels (screen stable; scaled
    with the scene in an export), vertically centered on its position and
    left of it when ``align_right`` (the left half of a circular tree)."""

    def __init__(self, text: str, font: QtGui.QFont, *, align_right: bool = False) -> None:
        super().__init__(text)
        self.setFont(font)
        self._align_right = bool(align_right)
        # The label on one line ("Cluster 4 · 83 genomes").
        self.full_text = text
        self.setData(_SCENE_DATA_KIND, "collapsed-label")

    def _offset(self) -> QtCore.QPointF:
        rect = super().boundingRect()
        return QtCore.QPointF(-rect.width() if self._align_right else 0.0, -rect.height() / 2.0)

    def boundingRect(self) -> QtCore.QRectF:
        return super().boundingRect().translated(self._offset())

    def shape(self) -> QtGui.QPainterPath:
        path = QtGui.QPainterPath()
        path.addRect(self.boundingRect())
        return path

    def paint(self, painter, option, widget=None) -> None:
        painter.translate(self._offset())
        super().paint(painter, option, widget)


class _CollapseHandleItem(QtWidgets.QGraphicsItem):
    """The small round handle on the node where a clade or cluster starts
    (or on a collapsed node): click collapses or expands that subtree. It
    is a faint ring until the pointer is over it, then a filled dot with a
    minus (or plus, when collapsed) in the group's color. Device pixels
    (screen stable); never exported."""

    def __init__(
        self,
        owner: "RaxmlTreeView",
        members: Sequence[str],
        color: QtGui.QColor,
        *,
        collapsed: bool,
        label: str,
    ) -> None:
        super().__init__()
        self._owner = owner
        self._members = tuple(members)
        self._color = QtGui.QColor(color)
        self._collapsed = bool(collapsed)
        self._hovered = False
        self.label = label
        self.setData(_SCENE_DATA_KIND, "collapse-handle")
        self.setData(_SCENE_DATA_NAME, label)
        self.setAcceptHoverEvents(True)
        self.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.setZValue(_HOVER_Z + 4.0)
        count = len(self._members)
        noun = "genome" if count == 1 else "genomes"
        action = "Expand" if collapsed else "Collapse"
        self.setToolTip(f"{action} {label} ({count} {noun})" if label else f"{action} these {count} {noun}")

    def members(self) -> Tuple[str, ...]:
        return self._members

    def is_collapsed(self) -> bool:
        return self._collapsed

    def boundingRect(self) -> QtCore.QRectF:
        r = _HANDLE_HIT_RADIUS_PX
        return QtCore.QRectF(-r, -r, 2.0 * r, 2.0 * r)

    def shape(self) -> QtGui.QPainterPath:
        path = QtGui.QPainterPath()
        path.addEllipse(self.boundingRect())
        return path

    def paint(self, painter, option, widget=None) -> None:
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        strong = self._hovered or self._collapsed
        radius = _HANDLE_HOVER_RADIUS_PX if self._hovered else (
            _HANDLE_COLLAPSED_RADIUS_PX if self._collapsed else _HANDLE_RADIUS_PX
        )
        if strong:
            ring = QtGui.QPen(QtGui.QColor(_TREE_CANVAS_COLOR))
            ring.setWidthF(1.2)
            painter.setPen(ring)
            painter.setBrush(QtGui.QBrush(self._color))
        else:
            # Unobtrusive: a thin ring in the group's color over the canvas.
            pen_color = QtGui.QColor(self._color)
            pen_color.setAlphaF(pen_color.alphaF() * 0.75)
            pen = QtGui.QPen(pen_color)
            pen.setWidthF(1.2)
            painter.setPen(pen)
            painter.setBrush(QtGui.QBrush(QtGui.QColor(_TREE_CANVAS_COLOR)))
        painter.drawEllipse(QtCore.QPointF(0.0, 0.0), radius, radius)
        if not strong:
            return
        glyph = QtGui.QPen(QtGui.QColor(_TREE_CANVAS_COLOR))
        glyph.setWidthF(1.4)
        glyph.setCapStyle(QtCore.Qt.PenCapStyle.RoundCap)
        painter.setPen(glyph)
        arm = radius * 0.5
        painter.drawLine(QtCore.QPointF(-arm, 0.0), QtCore.QPointF(arm, 0.0))
        if self._collapsed:
            painter.drawLine(QtCore.QPointF(0.0, -arm), QtCore.QPointF(0.0, arm))

    def hoverEnterEvent(self, event: Optional[QtWidgets.QGraphicsSceneHoverEvent]) -> None:
        self._hovered = True
        # Shown at full strength over a tree dimmed for the hover.
        self.setFlag(_IGNORES_PARENT_OPACITY, True)
        self.update()
        if self._owner.hover_interactions_enabled():
            self._owner._on_subtree_hover_enter(self._members)
        super().hoverEnterEvent(event)

    def hoverLeaveEvent(self, event: Optional[QtWidgets.QGraphicsSceneHoverEvent]) -> None:
        self._hovered = False
        self.setFlag(_IGNORES_PARENT_OPACITY, False)
        self.update()
        if self._owner.hover_interactions_enabled():
            self._owner._on_subtree_hover_leave(self._members)
        super().hoverLeaveEvent(event)

    def mousePressEvent(self, event: Optional[QtWidgets.QGraphicsSceneMouseEvent]) -> None:
        if event is not None and event.button() == QtCore.Qt.MouseButton.LeftButton:
            members = self._members
            # Redrawn later: the scene must not delete this item in its own event.
            QtCore.QTimer.singleShot(0, lambda: self._owner.toggle_collapsed(members))
            event.accept()
            return
        super().mousePressEvent(event)

    def mouseDoubleClickEvent(self, event: Optional[QtWidgets.QGraphicsSceneMouseEvent]) -> None:
        # A quick second click is another toggle, not the view's refit.
        self.mousePressEvent(event)

    def contextMenuEvent(self, event: Optional[QtWidgets.QGraphicsSceneContextMenuEvent]) -> None:
        if event is not None:
            self._owner._show_node_menu(self._members, event.screenPos())
            event.accept()


class _BranchBreakItem(QtWidgets.QGraphicsItem):
    """The break mark of a shortened edge: a gap in the branch crossed by
    two short parallel slanted strokes, in device pixels (screen stable;
    scaled with the scene in an export), with the edge's true length as a
    child label above it, or below it (``label_below``) for every second
    broken edge, so the labels of nearby broken edges never overlap."""

    def __init__(self, true_length: float, label_text: str, *, label_below: bool = False) -> None:
        super().__init__()
        self._half = _BREAK_MARK_HALF_HEIGHT_PX
        self._gap = _BREAK_MARK_GAP_PX
        self.setData(_SCENE_DATA_KIND, "branch-break")
        self.setData(_SCENE_DATA_VALUE, float(true_length))
        # Above the branches (a hovered path is _HOVER_Z).
        self.setZValue(_HOVER_Z + 5.0)
        label = QtWidgets.QGraphicsSimpleTextItem(label_text, self)
        font = QtGui.QFont()
        font.setPointSizeF(_BREAK_LABEL_FONT_SIZE)
        label.setFont(font)
        label.setBrush(QtGui.QBrush(QtGui.QColor(COLORS["text_muted"])))
        bounds = label.boundingRect()
        label_y = self._half + 2.0 if label_below else -self._half - 2.0 - bounds.height()
        label.setPos(-bounds.width() / 2.0, label_y)
        label.setData(_SCENE_DATA_KIND, "branch-break-label")
        label.setData(_SCENE_DATA_VALUE, float(true_length))
        self.label = label

    def boundingRect(self) -> QtCore.QRectF:
        width = self._gap + self._half + 2.0
        return QtCore.QRectF(-width, -self._half - 1.0, 2.0 * width, 2.0 * self._half + 2.0)

    def paint(self, painter, option, widget=None) -> None:
        half, gap = self._half, self._gap
        # The gap: the branch under it is covered with the canvas color.
        painter.fillRect(QtCore.QRectF(-gap, -2.0, 2.0 * gap, 4.0), QtGui.QColor(_TREE_CANVAS_COLOR))
        pen = QtGui.QPen(QtGui.QColor(COLORS["text"]))
        pen.setWidthF(1.3)
        pen.setCapStyle(QtCore.Qt.PenCapStyle.RoundCap)
        painter.setPen(pen)
        slant = half * 0.6
        for x in (-gap, gap):
            painter.drawLine(QtCore.QPointF(x - slant, half), QtCore.QPointF(x + slant, -half))


class _TreeContentItem(QtWidgets.QGraphicsItem):
    """The one parent of every drawn tree item: a hover dims the whole tree
    with a single ``setOpacity`` on it (the hovered items ignore it). It
    draws nothing and has no bounds of its own, so fits and exports see
    only its children."""

    def __init__(self) -> None:
        super().__init__()
        self.setFlag(QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemHasNoContents, True)

    def boundingRect(self) -> QtCore.QRectF:
        return QtCore.QRectF()

    def paint(self, painter, option, widget=None) -> None:  # pragma: no cover - no contents
        return


def _postorder_nodes(root: _TreeNode) -> List[_TreeNode]:
    """Every node of ``root``'s tree after its children."""
    order: List[_TreeNode] = []
    stack = [root]
    while stack:
        node = stack.pop()
        order.append(node)
        stack.extend(node.children)
    order.reverse()
    return order


@dataclass
class HoverOverlay:
    """Where the hover emphasis is drawn over the tree, in viewport pixels:
    one dot per hovered tip on screen, at most ``_HOVER_CHIP_CAP`` name
    chips and a "+N more" chip for the rest."""

    dots: List[Tuple[str, QtCore.QRectF]] = field(default_factory=list)
    chips: List[Tuple[str, str, QtCore.QRectF]] = field(default_factory=list)
    more: Optional[Tuple[str, QtCore.QRectF]] = None
    # Leader lines (dot center, chip edge) for chips the stack moved off their tip.
    leaders: List[Tuple[QtCore.QPointF, QtCore.QPointF]] = field(default_factory=list)
    # Hovered tips scrolled out of view: one marker per side they are on
    # ("above", "below", "left", "right"), its text and rect at that edge, so
    # a dimmed tree always says where the hovered samples are.
    edges: List[Tuple[str, str, QtCore.QRectF]] = field(default_factory=list)


class RaxmlTreeView(QtWidgets.QGraphicsView):
    viewChanged = QtCore.pyqtSignal()
    leafHovered = QtCore.pyqtSignal(object)
    subtreeHovered = QtCore.pyqtSignal(object)
    leafDoubleClicked = QtCore.pyqtSignal(object)
    cladeClicked = QtCore.pyqtSignal(object)
    clusterClicked = QtCore.pyqtSignal(object)
    # A clade or cluster label's menu: open it in a new explore window (plan 3.10).
    cladeWindowRequested = QtCore.pyqtSignal(object)
    clusterWindowRequested = QtCore.pyqtSignal(object)
    branchDistanceAvailabilityChanged = QtCore.pyqtSignal(bool, str)
    branchScaleModeChanged = QtCore.pyqtSignal(bool)
    # Whether the Distance drawing has long branches to break (the rule
    # finds some, whether or not true lengths are shown).
    branchBreaksAvailabilityChanged = QtCore.pyqtSignal(bool)

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
        reroot_for_groups: bool = True,
        prune_inactive: bool = False,
        root_mode: str = ROOT_AS_IN_FILE,
        order_mode: str = ORDER_AS_COMPUTED,
        scale_bar_position: str = "top",
    ) -> None:
        super().__init__(parent)
        self._scene = QtWidgets.QGraphicsScene(self)
        self.setScene(self._scene)
        self.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        self.setBackgroundBrush(QtGui.QColor(_TREE_CANVAS_COLOR))
        self.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        # Navigation is the view's own (wheelEvent, the mouse and key
        # events): no scroll bars, no hand drag mode, zooms anchored by
        # _scale_view.
        self.setDragMode(QtWidgets.QGraphicsView.DragMode.NoDrag)
        self.setViewportUpdateMode(QtWidgets.QGraphicsView.ViewportUpdateMode.BoundingRectViewportUpdate)
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.setTransformationAnchor(QtWidgets.QGraphicsView.ViewportAnchor.NoAnchor)
        self.setResizeAnchor(QtWidgets.QGraphicsView.ViewportAnchor.AnchorViewCenter)
        self.setFocusPolicy(QtCore.Qt.FocusPolicy.StrongFocus)
        self.setMinimumWidth(280)
        mode = (layout_mode or "circular").strip().lower()
        self._layout_mode = mode if mode in {"circular", "rectangular"} else "circular"
        self._preserve_branch_lengths = bool(preserve_branch_lengths)
        self._optimize_groups = bool(optimize_groups)
        self._use_group_backgrounds = bool(use_group_backgrounds)
        self._show_group_labels = bool(show_group_labels)
        self._show_cluster_annotations = bool(show_cluster_annotations)
        self._enable_hover_interactions = bool(enable_hover_interactions)
        # The group optimizer may try every edge as the root to make groups
        # contiguous; without this it only rotates children (the file root,
        # or the midpoint root, stays).
        self._reroot_for_groups = bool(reroot_for_groups)
        # Inactive samples (set_active_samples) are pruned from the drawn
        # tree instead of greyed.
        self._prune_inactive = bool(prune_inactive)
        self._root_mode = root_mode if root_mode in ROOT_MODES else ROOT_AS_IN_FILE
        # The leaf order (phylo ORDER_*): ladderized by clade size or depth
        # with the clade, cluster and branch color bands kept together, or as
        # the root and clade rotation leave it.
        self._order_mode = order_mode if order_mode in ORDER_MODES else ORDER_AS_COMPUTED
        self._scale_bar_position = "bottom" if scale_bar_position == "bottom" else "top"
        # The scale bar's unit (None: generic branch-length units).
        self._branch_length_unit: Optional[str] = None
        # What the file's numeric internal labels are ("UFBoot", "SH-aLRT";
        # None: no support) and the value from which an edge is well
        # supported; shown in a branch's hover text only.
        self._support_kind: Optional[str] = None
        self._support_threshold: Optional[float] = None
        self._support_key: Optional[Tuple[object, ...]] = None
        self._support_map: Dict[frozenset, float] = {}
        self._support_leaves: frozenset = frozenset()
        # Another fit's branch lengths drawn on the loaded topology
        # (set_branch_lengths; None: the file's own), the scale bar's name for
        # the lengths shown, why the last ones asked for are not drawn, and
        # the drawn tree with them (the resolved root, its key, the copy).
        self._lengths_root: Optional[_TreeNode] = None
        self._lengths_version = 0
        self._lengths_label: Optional[str] = None
        self._own_lengths_label: Optional[str] = None
        self._lengths_error = ""
        self._lengths_display: Optional[Tuple[_TreeNode, Tuple[object, ...], _TreeNode]] = None
        # The long branch rule (Distance only): the edges the drawing breaks,
        # and True to draw every edge at its true length instead.
        self._show_true_branch_lengths = False
        self._branch_breaks = BranchBreaks()
        self._branch_breaks_drawn = False
        # The tooltip of each broken edge (subtree leaves -> text).
        self._branch_break_notes: Dict[Tuple[str, ...], str] = {}
        self._branch_breaks_announced = False

        self._root: Optional[_TreeNode] = None
        self._root_version = 0
        self._branch_distance_available = False
        self._branch_distance_status = "Load a tree to inspect branch distances."
        self._device_stable_items: List[QtWidgets.QGraphicsItem] = []
        # When True, the whole tree is scaled to fit the viewport on render/resize.
        # Zooming, stretching or panning turns this off until Fit (or the
        # next tree/layout load).
        self._auto_fit = True
        self._fitting = False  # re-entrancy guard for fitInView/resize
        # The camera zoom over the layout scales (the view transform is
        # their product; 1 at every fit) and the leaf labels' last squeeze
        # (_apply_camera_to_items; None: not applied to the drawn labels).
        self._camera_zoom = 1.0
        self._label_squeeze: Optional[Tuple[float, float]] = None
        # A left or middle drag on the background: where it last was, its
        # button, and whether it moved far enough to pan.
        self._pan_from: Optional[QtCore.QPointF] = None
        self._pan_button = QtCore.Qt.MouseButton.NoButton
        self._panning = False
        # The sub-pixel part of the scrolls so far (small wheel deltas add up).
        self._scroll_rest = (0.0, 0.0)
        # The leaves a redraw keeps at their place on screen (a collapsed or
        # expanded subtree; None: the row nearest the middle).
        self._redraw_anchor: Optional[frozenset] = None
        # The fit's content rect for the drawn tree (None: not measured).
        self._fit_rect_cache: Optional[QtCore.QRectF] = None
        # Batches viewChanged while several scroll bars move as one.
        self._quiet_scroll = False
        # The scales the chips were last decluttered at.
        self._declutter_key: Optional[Tuple[float, float, float, int]] = None
        # A branch squeeze grows the collapsed labels in scene units: the
        # tree is drawn again when the squeezing stops.
        self._restretch_timer = QtCore.QTimer(self)
        self._restretch_timer.setSingleShot(True)
        self._restretch_timer.setInterval(250)
        self._restretch_timer.timeout.connect(self._render_tree_if_drawn)
        # The last zoom's anchor: (viewport point, scene point, camera after it).
        self._zoom_anchor: Optional[Tuple[QtCore.QPointF, QtCore.QPointF, Tuple[int, int, float, float]]] = None
        self._active_samples: Optional[set[str]] = None
        self._reference_sample: Optional[str] = None
        self._leaf_group_map: Dict[str, str] = {}
        self._group_color_map: Dict[str, QtGui.QColor] = {}
        self._leaf_display_group_map: Dict[str, str] = {}
        self._display_group_color_map: Dict[str, QtGui.QColor] = {}
        # Display groups color the leaf labels too (False: branches only,
        # the RaxmlTreeView branch_colors_only option).
        self._display_groups_color_labels = True
        self._leaf_cluster_map: Dict[str, str] = {}
        self._leaf_hover_metadata: Dict[str, Dict[str, str]] = {}
        self._bottom_legend_sections: List[Dict[str, object]] = []
        self._render_root_cache_key: Optional[Tuple[object, ...]] = None
        self._render_root_cache: Optional[_TreeNode] = None
        # The ladderized drawing of the grouped root (_resolved_render_root);
        # _band_version counts the changes of the group, cluster and display
        # group maps it keeps together. With pruning on, the full tree's
        # ladder order (_full_ladder_order) is followed, so a filter only
        # removes rows.
        self._band_version = 0
        self._ordered_root_key: Optional[Tuple[object, ...]] = None
        self._ordered_root: Optional[_TreeNode] = None
        self._full_order_key: Optional[Tuple[object, ...]] = None
        self._full_order: Optional[List[str]] = None
        self._base_root_cache_key: Optional[Tuple[object, ...]] = None
        self._base_root_cache: Optional[_TreeNode] = None
        # The full tree's group root (``root_mode="groups"``), kept so a
        # filter only prunes it: the root never moves with the filter.
        self._group_full_root_key: Optional[Tuple[object, ...]] = None
        self._group_full_root: Optional[_TreeNode] = None
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
        # The parent of every tree item (dimmed while a hover is active).
        self._content_item: Optional[_TreeContentItem] = None
        self._hover_dimmed = False
        # Collapsed subtrees, each kept as the set of its leaf names (in the
        # order they were collapsed): found again in every drawing of the
        # tree (another order, a filter), not drawn while the tree has no
        # such subtree (another tree source or root).
        self._collapsed: List[frozenset] = []
        # Of the last drawing: the collapsed sets it has no subtree for, and
        # each other set's drawn leaves (fewer in a tree pruned to the
        # sample set).
        self._collapsed_unmatched: Set[frozenset] = set()
        self._collapsed_drawn_of: Dict[frozenset, frozenset] = {}
        self._reset_collapse_drawing()

        self._base_edge_pen = QtGui.QPen(QtGui.QColor(_TREE_BRANCH_COLOR))
        self._base_edge_pen.setWidthF(1.0)
        self._base_edge_pen.setCapStyle(QtCore.Qt.PenCapStyle.RoundCap)
        self._base_edge_pen.setJoinStyle(QtCore.Qt.PenJoinStyle.RoundJoin)
        # Cosmetic: the hovered path is _HOVER_PATH_WIDTH_PX on screen at
        # every zoom (a scene-width pen was under a pixel at the fit zoom).
        self._hover_edge_pen = QtGui.QPen(QtGui.QColor(_TREE_HOVER_COLOR))
        self._hover_edge_pen.setWidthF(_HOVER_PATH_WIDTH_PX)
        self._hover_edge_pen.setCosmetic(True)
        self._hover_edge_pen.setCapStyle(QtCore.Qt.PenCapStyle.RoundCap)
        self._hover_edge_pen.setJoinStyle(QtCore.Qt.PenJoinStyle.RoundJoin)

        # "Rows 1 to 60 of 153" while the rectangular tree shows only part
        # of its rows (set_row_note_enabled; off by default).
        self._row_note: Optional[RowRangeNote] = None
        # The (x, y) scales of the last fit (None: none yet), the stretch limits.
        self._fit_scales: Optional[Tuple[float, float]] = None
        self.viewChanged.connect(self._update_row_note)
        self.viewChanged.connect(self._declutter_group_chips)

        self._draw_error("No tree loaded.")

    def hover_interactions_enabled(self) -> bool:
        return self._enable_hover_interactions

    def layout_mode(self) -> str:
        return self._layout_mode

    def set_layout_mode(self, mode: str) -> None:
        normalized = str(mode or "").strip().lower()
        if normalized not in {"circular", "rectangular"}:
            raise ValueError(f"Unsupported tree layout: {mode!r}")
        if normalized == self._layout_mode:
            return
        self._layout_mode = normalized
        self._auto_fit = True
        if self._root is not None:
            self._render_tree()

    def preserve_branch_lengths_enabled(self) -> bool:
        return self._preserve_branch_lengths

    def branch_distance_available(self) -> bool:
        return self._branch_distance_available

    def branch_distance_status_message(self) -> str:
        return self._branch_distance_status

    @staticmethod
    def _analyze_branch_lengths(root: _TreeNode) -> Tuple[bool, str]:
        return analyze_branch_lengths(root)

    def _set_branch_distance_status(self, available: bool, message: str) -> None:
        self._branch_distance_available = bool(available)
        self._branch_distance_status = str(message)
        self.branchDistanceAvailabilityChanged.emit(
            self._branch_distance_available,
            self._branch_distance_status,
        )

    def load_tree_file(
        self,
        path: Path,
        *,
        root_mode: str,
        branch_length_unit: Optional[str] = None,
        support_kind: Optional[str] = None,
        support_threshold: Optional[float] = None,
        lengths_file: Optional[Path] = None,
        lengths_label: Optional[str] = None,
        own_lengths_label: Optional[str] = None,
    ) -> None:
        """Load a tree source in one render: its file, the root to draw, the
        unit the scale bar names (None: generic branch lengths), its
        support (``support_kind`` such as "UFBoot" or "SH-aLRT", well
        supported from ``support_threshold``; None: no support hover) and
        its branch lengths (``set_branch_lengths``: ``lengths_file`` draws
        another fit's lengths on this topology, ``lengths_label`` and
        ``own_lengths_label`` name its lengths and the file's own on the
        scale bar)."""
        mode = str(root_mode or "").strip().lower()
        if mode not in ROOT_MODES:
            raise ValueError(f"Unsupported tree root: {root_mode!r}")
        self._root_mode = mode
        self._branch_length_unit = str(branch_length_unit or "").strip() or None
        self._support_kind = str(support_kind or "").strip() or None
        self._support_threshold = float(support_threshold) if support_threshold is not None else None
        self._lengths_label = str(lengths_label or "").strip() or None
        self._own_lengths_label = str(own_lengths_label or "").strip() or None
        self.set_tree_file(path, lengths_file=lengths_file)

    def set_tree_file(self, path: Path, *, lengths_file: Optional[Path] = None) -> None:
        try:
            text = Path(path).read_text(encoding="utf-8")
        except Exception as exc:
            self.clear_tree(f"Failed to load tree: {exc}")
            return
        # An unreadable lengths file draws the tree with its own lengths.
        lengths_text: Optional[str] = None
        lengths_error = ""
        if lengths_file is not None:
            try:
                lengths_text = Path(lengths_file).read_text(encoding="utf-8")
            except OSError as exc:
                lengths_error = f"Could not read the branch lengths: {exc}"
        try:
            self.set_tree_text(text, lengths_text=lengths_text)
        except Exception as exc:
            self.clear_tree(f"Failed to load tree: {exc}")
            return
        if lengths_error:
            self._lengths_error = lengths_error

    def set_tree_text(self, text: str, *, lengths_text: Optional[str] = None) -> None:
        parser = _NewickParser(str(text or ""))
        root = parser.parse()
        self._root = root
        self._lengths_error = ""
        self._lengths_root = self._matching_lengths_tree(lengths_text) if lengths_text is not None else None
        self._lengths_version += 1
        available, message = self._analyze_branch_lengths(root)
        self._set_branch_distance_status(available, message)
        if self._preserve_branch_lengths and not available:
            self._preserve_branch_lengths = False
            self.branchScaleModeChanged.emit(False)
        self._root_version += 1
        self._auto_fit = True  # a freshly loaded tree starts fitted to the window
        self._invalidate_render_root_cache()
        self._render_tree()

    def set_branch_lengths(
        self,
        lengths_file: Optional[Path],
        *,
        label: Optional[str] = None,
        own_label: Optional[str] = None,
        unit: Optional[str] = None,
    ) -> bool:
        """Draw the loaded tree with the branch lengths of ``lengths_file``
        (another fit of the same topology, such as the core tree refit with
        recombination), or its own (None). Only lengths change: the leaf
        order, root, clusters and support stay the loaded tree's, and the
        long branch rule and the scale bar follow the new lengths. ``label``
        names those lengths on the scale bar and ``own_label`` the tree's
        own (drawn for None or a refused file), ``unit`` is their unit. False
        (and the tree's own lengths drawn) when the file cannot be read or
        is not the same tree (``branch_lengths_error``)."""
        self._lengths_label = str(label or "").strip() or None
        self._own_lengths_label = str(own_label or "").strip() or None
        if unit is not None:
            self._branch_length_unit = str(unit).strip() or None
        lengths_root: Optional[_TreeNode] = None
        if lengths_file is not None and self._root is not None:
            try:
                text = Path(lengths_file).read_text(encoding="utf-8")
            except OSError as exc:
                self._lengths_error = f"Could not read the branch lengths: {exc}"
            else:
                lengths_root = self._matching_lengths_tree(text)
        elif lengths_file is None:
            self._lengths_error = ""
        self._lengths_root = lengths_root
        self._lengths_version += 1
        if self._root is not None:
            self._render_tree()
        return lengths_file is None or lengths_root is not None

    def branch_lengths_error(self) -> str:
        """Why the last branch lengths asked for are not drawn; empty when
        they are."""
        return self._lengths_error

    def branch_lengths_replaced(self) -> bool:
        """True while another fit's branch lengths are drawn."""
        return self._lengths_root is not None

    def _split_lengths(
        self, root: _TreeNode, keep: Optional[frozenset] = None
    ) -> Tuple[Dict[frozenset, float], frozenset]:
        """Each edge's length by the split it makes of the leaves ``keep``
        (every leaf by default; the side without the first leaf, so any
        rooting finds it), and the leaves split. A split several edges make
        (a pruned tree's merged chain, or the two root edges) gets their
        lengths summed; an edge that does not split ``keep`` is left out."""
        leaves_of: Dict[int, frozenset] = {}
        order = _postorder_nodes(root)
        for node in order:
            if node.is_leaf:
                leaves_of[id(node)] = frozenset((self._canonical_leaf_name(node.name or ""),))
            else:
                leaves_of[id(node)] = frozenset().union(*(leaves_of[id(child)] for child in node.children))
        every = leaves_of[id(root)] if keep is None else keep
        if not every:
            return {}, every
        anchor = min(every)
        lengths: Dict[frozenset, float] = {}
        for node in order:
            if node is root:
                continue
            side = leaves_of[id(node)] if keep is None else leaves_of[id(node)] & keep
            if not side or len(side) >= len(every):
                continue
            key = every - side if anchor in side else side
            lengths[key] = lengths.get(key, 0.0) + edge_length(node)
        return lengths, every

    def _matching_lengths_tree(self, text: str) -> Optional[_TreeNode]:
        """``text`` parsed when it is the loaded tree with other branch
        lengths (the same leaves and splits, usable lengths); else None and
        ``_lengths_error`` says why."""
        if self._root is None:
            return None
        try:
            lengths_root = _NewickParser(str(text or "")).parse()
        except Exception as exc:
            self._lengths_error = f"Could not read the branch lengths: {exc}"
            return None
        own, own_leaves = self._split_lengths(self._root)
        other, other_leaves = self._split_lengths(lengths_root)
        if own_leaves != other_leaves or own.keys() != other.keys():
            self._lengths_error = "These branch lengths belong to a different tree (other genomes or topology)."
            return None
        available, message = analyze_branch_lengths(lengths_root)
        if not available:
            self._lengths_error = message
            return None
        self._lengths_error = ""
        return lengths_root

    def _with_branch_lengths(self, root: _TreeNode) -> _TreeNode:
        """The drawn tree ``root`` with the lengths of ``set_branch_lengths``
        (a copy, same node order; ``root`` itself without them). Every drawn
        edge takes the length of the same split in the other fit; edges that
        share a split (the two edges of the drawn root, a node left with one
        child) share its length in the proportion of their own, so the root
        stays on the same edge at the same relative place."""
        lengths_root = self._lengths_root
        if lengths_root is None:
            return root
        cached = self._lengths_display
        key = (self._lengths_version, self._reference_sample)
        if cached is not None and cached[0] is root and cached[1] == key:
            return cached[2]
        display = _clone_tree(root)
        order = _postorder_nodes(display)
        leaves_of: Dict[int, frozenset] = {}
        for node in order:
            if node.is_leaf:
                leaves_of[id(node)] = frozenset((self._canonical_leaf_name(node.name or ""),))
            else:
                leaves_of[id(node)] = frozenset().union(*(leaves_of[id(child)] for child in node.children))
        keep = leaves_of[id(display)]
        if keep:
            lengths, _every = self._split_lengths(lengths_root, keep)
            anchor = min(keep)
            edges: Dict[frozenset, List[_TreeNode]] = {}
            for node in order:
                side = leaves_of[id(node)]
                if node is display or len(side) >= len(keep):
                    continue
                edges.setdefault(keep - side if anchor in side else side, []).append(node)
            for split, nodes in edges.items():
                total = lengths.get(split)
                if total is None:
                    continue
                own = [edge_length(node) for node in nodes]
                own_total = sum(own)
                for node, value in zip(nodes, own):
                    share = value / own_total if own_total > 0.0 else 1.0 / len(nodes)
                    node.branch_length = total * share
        self._lengths_display = (root, key, display)
        return display

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

    def set_preserve_branch_lengths(self, enabled: bool) -> bool:
        enabled = bool(enabled)
        if enabled and self._root is not None and not self._branch_distance_available:
            enabled = False
        if enabled == self._preserve_branch_lengths:
            return self._preserve_branch_lengths
        self._preserve_branch_lengths = enabled
        if self._root is not None:
            self._render_tree()
        self.branchScaleModeChanged.emit(self._preserve_branch_lengths)
        return self._preserve_branch_lengths

    def show_true_branch_lengths(self) -> bool:
        """True when long branches are drawn at their true length (the
        "True branch lengths" toggle), False when the rule breaks them."""
        return self._show_true_branch_lengths

    def set_show_true_branch_lengths(self, enabled: bool) -> None:
        enabled = bool(enabled)
        if enabled == self._show_true_branch_lengths:
            return
        self._show_true_branch_lengths = enabled
        self._auto_fit = True
        if self._root is not None:
            self._render_tree()

    def branch_breaks_available(self) -> bool:
        """Whether the drawn tree has long branches the Distance drawing
        breaks (True also while true lengths are shown)."""
        return bool(self._branch_breaks)

    def branch_breaks_drawn(self) -> bool:
        """Whether the drawing shows shortened, marked branches now."""
        return self._branch_breaks_drawn

    def broken_branches(self) -> List[Tuple[Tuple[str, ...], float]]:
        """The broken edges of the drawn tree: the leaves under each and
        its true length (empty when none, or while true lengths are shown)."""
        if not self._branch_breaks_drawn:
            return []
        names = {id(leaf): self._canonical_leaf_name(leaf.name or "") for leaf in self._iter_leaves(self._branch_breaks.nodes)}
        out: List[Tuple[Tuple[str, ...], float]] = []
        for node, length in zip(self._branch_breaks.nodes, self._branch_breaks.true_lengths):
            out.append((tuple(names[id(leaf)] for leaf in self._iter_leaves((node,))), length))
        return out

    @staticmethod
    def _iter_leaves(nodes: Sequence[_TreeNode]) -> List[_TreeNode]:
        out: List[_TreeNode] = []
        stack = list(reversed(list(nodes)))
        while stack:
            node = stack.pop()
            if node.is_leaf:
                out.append(node)
            else:
                stack.extend(reversed(node.children))
        return out

    def _announce_branch_breaks(self) -> None:
        available = bool(self._branch_breaks)
        if available != self._branch_breaks_announced:
            self._branch_breaks_announced = available
            self.branchBreaksAvailabilityChanged.emit(available)

    def _short_unit(self) -> str:
        """The scale bar unit, short enough for a label on a branch."""
        unit = self._branch_length_unit or ""
        if unit.startswith("substitutions/site"):
            return "subst/site"
        if unit.startswith("SNPs"):
            return "SNPs"
        return unit

    def camera_zoom(self) -> float:
        """The camera zoom over the drawing (1: as fitted, or as stretched)."""
        return self._camera_zoom

    def is_fitted(self) -> bool:
        """Whether the tree is fitted to the panel (refitted on resize)."""
        return self._auto_fit

    def _layout_scales(self) -> Tuple[float, float]:
        """The (x, y) scales of the drawing without the camera zoom: the
        fit's, or as Ctrl+Shift and Alt with the wheel stretched them."""
        transform = self.transform()
        zoom = self._camera_zoom
        return abs(float(transform.m11())) / zoom, abs(float(transform.m22())) / zoom

    def _stable_scale(self) -> float:
        """How much the screen stable annotations (chips, handles, break
        marks, the scale bar) grow with the camera zoom."""
        return min(self._camera_zoom, _STABLE_ZOOM_MAX)

    def zoom_by(self, factor: float, anchor: Optional[QtCore.QPointF] = None) -> None:
        """Camera zoom by ``factor`` around ``anchor`` (viewport pixels;
        default the viewport's center): the whole canvas scales like an
        image, the branches, labels, chips and marks together. Leaves the
        fitted mode. Limits: half the fit's size, and leaf labels at
        _LABEL_MAX_SCALE times their size."""
        factor = float(factor)
        if not math.isfinite(factor) or factor <= 0.0:
            raise ValueError("Zoom factor must be a positive finite number.")
        if not self._leaf_order:
            return
        zoom = self._camera_zoom
        scale_x, scale_y = self._layout_scales()
        top = self._camera_zoom_top(scale_x, scale_y)
        # A limit never moves the zoom against the step (a stretch may have
        # left it past one).
        floor = min(zoom, self._camera_zoom_floor(scale_x, scale_y))
        target = min(max(zoom, top), max(floor, zoom * factor))
        if abs(target - zoom) <= 1e-9 * zoom:
            return
        self._camera_zoom = target
        step = target / zoom
        self._scale_view(step, step, anchor)

    @staticmethod
    def _camera_zoom_top(scale_x: float, scale_y: float) -> float:
        """The most camera zoom: leaf labels (drawn at the smaller layout
        scale) at _LABEL_MAX_SCALE times their size, and never less than 4."""
        label_scale = min(scale_x, scale_y)
        return max(4.0, _LABEL_MAX_SCALE / label_scale) if label_scale > 0.0 else 4.0

    def _camera_zoom_floor(self, scale_x: float, scale_y: float) -> float:
        """The least camera zoom: the drawing _CAMERA_ZOOM_MIN of the size
        that fits the viewport as it is now (the panel may have changed
        size since the last fit)."""
        viewport = self.viewport()
        if self._fit_rect_cache is None:
            self._fit_rect_cache = self._fit_rect()
        rect = self._fit_rect_cache
        if viewport is None or rect is None or scale_x <= 0.0 or scale_y <= 0.0:
            return _CAMERA_ZOOM_MIN
        fits = min(
            float(viewport.width()) / max(1.0, float(rect.width()) * scale_x),
            float(viewport.height()) / max(1.0, float(rect.height()) * scale_y),
        )
        return _CAMERA_ZOOM_MIN * fits

    def stretch_by(
        self,
        factor_x: float = 1.0,
        factor_y: float = 1.0,
        anchor: Optional[QtCore.QPointF] = None,
    ) -> None:
        """Stretch or squeeze the rectangular drawing along x (the branch
        lengths) and y (the row spacing) around ``anchor``; the labels and
        annotations keep their size. Leaves the fitted mode. Limits: a
        quarter of the fit to _STRETCH_MAX times it on each axis, and no
        stretch past the camera zoom's top (the leaf labels' largest size)."""
        if self._layout_mode != "rectangular" or not self._leaf_order:
            return
        scale_x, scale_y = self._layout_scales()
        if scale_x <= 0.0 or scale_y <= 0.0:
            return
        fit_x, fit_y = self._fit_scales or (scale_x, scale_y)

        def limited(scale: float, factor: float, fit: float) -> float:
            # Never against the step (the limits may have moved since).
            wanted = scale * float(factor)
            return min(max(scale, fit * _STRETCH_MAX), max(min(scale, fit * _STRETCH_MIN), wanted))

        new_x = limited(scale_x, factor_x, fit_x)
        new_y = limited(scale_y, factor_y, fit_y)
        zoom = self._camera_zoom
        # The leaf labels (at the smaller scale) stay within the zoom's top.
        most = max(min(scale_x, scale_y), _LABEL_MAX_SCALE / zoom)
        if zoom > self._camera_zoom_top(new_x, new_y) and min(new_x, new_y) > most:
            new_x = max(scale_x, most) if new_x > scale_x else new_x
            new_y = max(scale_y, most) if new_y > scale_y else new_y
        if abs(new_x - scale_x) <= 1e-9 * scale_x and abs(new_y - scale_y) <= 1e-9 * scale_y:
            return
        self._scale_view(new_x / scale_x, new_y / scale_y, anchor)
        if new_x < scale_x and self._collapsed_label_items:
            self._restretch_timer.start()

    def _scale_view(self, factor_x: float, factor_y: float, anchor: Optional[QtCore.QPointF]) -> None:
        """Scale the view and keep the scene point under ``anchor`` there;
        leaves the fitted mode."""
        viewport = self.viewport()
        if viewport is None:
            return
        if anchor is None:
            anchor = QtCore.QPointF(viewport.width() / 2.0, viewport.height() / 2.0)
        anchor = QtCore.QPointF(anchor)
        last = self._zoom_anchor
        if last is not None and last[0] == anchor and last[2] == self._camera_state():
            # Steps at one pointer position keep the first step's scene
            # point (whole pixel scrolls would let it drift).
            scene_point = last[1]
        else:
            scene_point = self._scene_point_at(anchor)
        quiet = self._quiet_scroll
        self._quiet_scroll = True
        try:
            self._leave_fit()
            self.scale(factor_x, factor_y)
            self._apply_camera_to_items()
            self._update_camera_scene_rect()
            self._put_scene_point_at(scene_point, anchor)
        finally:
            self._quiet_scroll = quiet
        self._zoom_anchor = (anchor, scene_point, self._camera_state())
        self.viewChanged.emit()

    def _camera_state(self) -> Tuple[int, int, float, float]:
        transform = self.transform()
        return (
            self.horizontalScrollBar().value(),
            self.verticalScrollBar().value(),
            float(transform.m11()),
            float(transform.m22()),
        )

    def _scene_point_at(self, viewport_point: QtCore.QPointF) -> QtCore.QPointF:
        """The scene point under a viewport point (sub-pixel)."""
        inverse, invertible = self.viewportTransform().inverted()
        return inverse.map(viewport_point) if invertible else QtCore.QPointF(self.mapToScene(viewport_point.toPoint()))

    def _put_scene_point_at(self, scene_point: QtCore.QPointF, viewport_point: QtCore.QPointF) -> None:
        """Scroll so ``scene_point`` is under ``viewport_point`` (no
        viewChanged: the caller sends it)."""
        moved = self.viewportTransform().map(scene_point) - viewport_point
        self._move_scroll_bars(int(round(moved.x())), int(round(moved.y())))

    def _move_scroll_bars(self, dx: int, dy: int) -> bool:
        """Move both scroll bars as one step without a viewChanged; whether
        the view moved."""
        bar_x, bar_y = self.horizontalScrollBar(), self.verticalScrollBar()
        before = (bar_x.value(), bar_y.value())
        quiet = self._quiet_scroll
        self._quiet_scroll = True
        try:
            if dx:
                bar_x.setValue(bar_x.value() + dx)
            if dy:
                bar_y.setValue(bar_y.value() + dy)
        finally:
            self._quiet_scroll = quiet
        return (bar_x.value(), bar_y.value()) != before

    def _scroll_by(self, dx: float, dy: float, *, keep_fit: bool = False) -> bool:
        """Pan the camera by (dx, dy) viewport pixels (the sub-pixel rest
        is kept for the next scroll); leaves the fitted mode unless
        ``keep_fit`` (then it only moves within the fitted drawing).
        Whether the view moved."""
        if not self._leaf_order:
            return False
        rest_x, rest_y = self._scroll_rest
        dx += rest_x
        dy += rest_y
        step_x, step_y = int(round(dx)), int(round(dy))
        self._scroll_rest = (dx - step_x, dy - step_y)
        if not step_x and not step_y:
            return False
        if not keep_fit:
            self._leave_fit()
        moved = self._move_scroll_bars(step_x, step_y)
        if moved:
            self.viewChanged.emit()
        return moved

    def _leave_fit(self) -> None:
        """From the fitted mode to a free camera, the view where it is."""
        if not self._auto_fit:
            return
        viewport = self.viewport()
        middle = QtCore.QPointF(viewport.width() / 2.0, viewport.height() / 2.0)
        center = self._scene_point_at(middle)
        self._auto_fit = False
        self._update_camera_scene_rect()
        self._put_scene_point_at(center, middle)

    def _update_camera_scene_rect(self) -> None:
        """The scene rect over everything drawn; out of the fitted mode it
        reaches a viewport beyond it on every side (less _PAN_KEEP_PX), so
        the camera can go anywhere that keeps some of the tree in view
        and a zoom anchored under the pointer is never clamped."""
        bounds = self._update_scene_rect_for_device_stable_items()
        viewport = self.viewport()
        if self._auto_fit or viewport is None or bounds.isNull() or bounds.isEmpty():
            return
        transform = self.transform()
        scale_x = abs(float(transform.m11())) or 1.0
        scale_y = abs(float(transform.m22())) or 1.0
        margin_x = max(0.0, float(viewport.width()) - _PAN_KEEP_PX) / scale_x
        margin_y = max(0.0, float(viewport.height()) - _PAN_KEEP_PX) / scale_y
        # From the drawing itself, not its padding (scene units: many
        # pixels at a high zoom).
        pad = _SCENE_RECT_PAD
        content = bounds.adjusted(pad, pad, -pad, -pad)
        self._scene.setSceneRect(content.adjusted(-margin_x, -margin_y, margin_x, margin_y))

    def _apply_camera_to_items(self) -> None:
        """Size the screen stable items for the camera zoom, and draw the
        leaf labels undistorted (a uniform scale, the smaller of the
        layout's two) so a stretch does not deform them."""
        scale = self._stable_scale()
        for item in self._device_stable_items:
            if item is not None and item.scene() is self._scene and abs(item.scale() - scale) > 1e-9:
                item.setScale(scale)
        scale_x, scale_y = self._layout_scales()
        uniform = min(scale_x, scale_y)
        squeeze = (1.0, 1.0)
        if self._layout_mode == "rectangular" and uniform > 0.0:
            squeeze = (uniform / scale_x, uniform / scale_y)
        if self._label_squeeze is not None and all(
            abs(a - b) <= 1e-9 for a, b in zip(squeeze, self._label_squeeze)
        ):
            return
        self._label_squeeze = squeeze
        self._squeeze_leaf_labels(*squeeze)

    def _squeeze_leaf_labels(self, factor_x: float, factor_y: float) -> None:
        """Scale every leaf label about its left middle (its row)."""
        for label in self._leaf_label_items.values():
            if label is None or label.scene() is not self._scene:
                continue
            if abs(factor_x - 1.0) <= 1e-9 and abs(factor_y - 1.0) <= 1e-9:
                label.setTransform(QtGui.QTransform())
                continue
            middle = label.boundingRect().height() / 2.0
            label.setTransform(QtGui.QTransform().translate(0.0, middle).scale(factor_x, factor_y).translate(0.0, -middle))

    def _cover_viewport_at(self, center: QtCore.QPointF, scale: float, scale_y: Optional[float] = None) -> None:
        """Grow the scene rect to the viewport's span around ``center``:
        vertically always, horizontally where the scene is narrower than the
        viewport. The view centers a scene smaller than itself, so the fit
        would not land where it is put (a fit by the height leaves the tree
        narrower than the viewport)."""
        viewport = self.viewport()
        scale_y = float(scale if scale_y is None else scale_y)
        if viewport is None or scale <= 0.0 or scale_y <= 0.0:
            return
        rect = QtCore.QRectF(self._scene.sceneRect())
        width = float(viewport.width()) / scale
        height = float(viewport.height()) / scale_y
        if rect.width() < width:
            rect.setLeft(min(rect.left(), center.x() - width / 2.0))
            rect.setRight(max(rect.right(), center.x() + width / 2.0))
        rect.setTop(min(rect.top(), center.y() - height / 2.0))
        rect.setBottom(max(rect.bottom(), center.y() + height / 2.0))
        self._scene.setSceneRect(rect)

    def zoom_in(self) -> None:
        self.zoom_by(_ZOOM_STEP)

    def zoom_out(self) -> None:
        self.zoom_by(1.0 / _ZOOM_STEP)

    def reset_view(self) -> None:
        self.fit_to_window()

    def set_optimize_groups(self, enabled: bool) -> None:
        enabled = bool(enabled)
        if enabled == self._optimize_groups:
            return
        self._optimize_groups = enabled
        self._invalidate_render_root_cache()
        if self._root is not None:
            self._render_tree()

    def root_mode(self) -> str:
        return self._root_mode

    def set_root_mode(self, mode: str) -> None:
        """Draw the tree rooted as in the file, at its midpoint, or on the
        edge that keeps the groups together (ROOT_GROUPS)."""
        normalized = str(mode or "").strip().lower()
        if normalized not in ROOT_MODES:
            raise ValueError(f"Unsupported tree root: {mode!r}")
        if normalized == self._root_mode:
            return
        self._root_mode = normalized
        self._invalidate_render_root_cache()
        if self._root is not None:
            self._render_tree()

    def order_mode(self) -> str:
        return self._order_mode

    def set_order_mode(self, mode: str) -> None:
        """Draw the leaves ladderized by clade size or by depth, or as
        computed (phylo ORDER_*); rotations only, the bands stay together."""
        normalized = str(mode or "").strip().lower()
        if normalized not in ORDER_MODES:
            raise ValueError(f"Unsupported tree order: {mode!r}")
        if normalized == self._order_mode:
            return
        self._order_mode = normalized
        if self._root is not None:
            self._render_tree()

    def set_prune_inactive(self, enabled: bool) -> None:
        """Prune the inactive samples from the drawn tree (True) or grey them."""
        enabled = bool(enabled)
        if enabled == self._prune_inactive:
            return
        self._prune_inactive = enabled
        self._invalidate_render_root_cache()
        if self._root is not None:
            self._render_tree()

    def branch_break_note(self, leaves: Sequence[str]) -> str:
        """The hover text of a broken edge ("Branch shortened for display;
        true length 0.0123 substitutions/site."); empty for any other edge."""
        return self._branch_break_notes.get(tuple(str(leaf) for leaf in leaves), "")

    def support_text(self, leaves: Sequence[str]) -> str:
        """The hover text of the edge above ``leaves`` ("UFBoot 97, well
        supported at ≥ 95"); empty without support or for an edge the file
        tree does not have (a pruned drawing)."""
        if self._support_kind is None or self._root is None:
            return ""
        key = (self._root_version, self._reference_sample)
        if self._support_key != key:
            self._support_map, self._support_leaves = _support_by_split(self._root, self._canonical_leaf_name)
            self._support_key = key
        every = self._support_leaves
        side = frozenset(str(leaf) for leaf in leaves)
        if not side or not every or not side < every:
            return ""
        split = every - side if min(every) in side else side
        value = self._support_map.get(split)
        if value is None:
            return ""
        text = f"{self._support_kind} {value:g}"
        limit = self._support_threshold
        if limit is None:
            return text
        if value >= limit:
            return f"{text}, well supported at \u2265 {limit:g}"
        return f"{text}, not well supported (< {limit:g})"

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
        self._band_version += 1
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
        if normalized_groups != self._leaf_display_group_map:
            self._band_version += 1
        self._leaf_display_group_map = normalized_groups
        self._display_group_color_map = normalized_colors
        if self._root is not None:
            self._render_tree()

    def set_display_groups_color_labels(self, enabled: bool) -> None:
        """Whether the display groups color the leaf labels as well as the
        branches (True, the default) or the branches only."""
        enabled = bool(enabled)
        if enabled == self._display_groups_color_labels:
            return
        self._display_groups_color_labels = enabled
        if self._root is not None:
            self._render_tree()

    def tree_leaf_names(self) -> List[str]:
        """The leaf names of the loaded tree file (before any pruning or
        rerooting), in file order; empty without a tree."""
        if self._root is None:
            return []
        names: List[str] = []
        stack = [self._root]
        while stack:
            node = stack.pop()
            if node.is_leaf:
                names.append(self._canonical_leaf_name(node.name or ""))
            else:
                stack.extend(reversed(node.children))
        return names

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
        self._band_version += 1
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
        """Every drawn leaf, top to bottom; a collapsed subtree's leaves are
        listed in its row, in their order (a linked heatmap keeps them all)."""
        return list(self._leaf_order)

    # ---- collapsed subtrees: a triangle in place of a clade or cluster ----

    def _leaf_set(self, leaves: Sequence[str]) -> frozenset:
        return frozenset(self._normalize_label(str(leaf)) for leaf in leaves or () if str(leaf).strip())

    def collapsed_groups(self) -> List[Tuple[str, ...]]:
        """The collapsed subtrees, each as its sorted leaf names, in the
        order they were collapsed (the panel state keeps them). A set the
        drawn tree has no subtree for (another tree source or root) is left
        out, though kept: it is drawn again on the way back."""
        return [tuple(sorted(members)) for members in self._collapsed if members not in self._collapsed_unmatched]

    def set_collapsed_groups(self, groups: Optional[Sequence[Sequence[str]]]) -> None:
        """Collapse these subtrees (each given by its leaf names) and expand
        every other. A set the drawn tree has no subtree for is not drawn
        and left out of ``collapsed_groups``."""
        wanted: List[frozenset] = []
        for group in groups or ():
            if isinstance(group, str):
                continue
            members = self._leaf_set(group)
            if len(members) >= 2 and members not in wanted:
                wanted.append(members)
        if wanted == self._collapsed:
            return
        self._collapsed = wanted
        self._redraw_for_collapse()

    def _stored_sets(self, members: frozenset) -> List[frozenset]:
        """The collapsed sets drawn as ``members`` (a pruned drawing draws
        a set as its leaves in the sample set), or ``members`` itself."""
        found = [stored for stored, drawn in self._collapsed_drawn_of.items() if drawn == members]
        if members in self._collapsed and members not in found:
            found.append(members)
        return found

    def is_collapsed(self, leaves: Sequence[str]) -> bool:
        return bool(self._stored_sets(self._leaf_set(leaves)))

    def collapse_subtree(self, leaves: Sequence[str]) -> None:
        """Draw the subtree of exactly these leaves as a triangle. In a
        tree pruned to the sample set, the whole subtree of the unpruned
        tree is kept collapsed (unless it holds another clade or cluster),
        so another filter still finds it."""
        members = self._leaf_set(leaves)
        if len(members) < 2 or self._stored_sets(members):
            return
        if self._prune_inactive and self._active_samples is not None and self._root is not None:
            full = full_subtree(self._unpruned_root(), members, frozenset(self._leaf_order), self._canonical_leaf_name)
            # Never one that takes in another clade or cluster.
            if all(
                {groups.get(name) for name in full} == {groups.get(name) for name in members}
                for groups in (self._leaf_group_map, self._leaf_cluster_map)
            ):
                members = full
        if members in self._collapsed:
            self._collapsed.remove(members)
        self._collapsed.append(members)
        self._redraw_anchor = members
        self._redraw_for_collapse()

    def expand_subtree(self, leaves: Sequence[str]) -> None:
        stored = self._stored_sets(self._leaf_set(leaves))
        if not stored:
            return
        self._collapsed = [members for members in self._collapsed if members not in stored]
        self._redraw_anchor = frozenset().union(*stored)
        self._redraw_for_collapse()

    def toggle_collapsed(self, leaves: Sequence[str]) -> None:
        if self.is_collapsed(leaves):
            self.expand_subtree(leaves)
        else:
            self.collapse_subtree(leaves)

    def expand_all(self) -> None:
        if not self._collapsed:
            return
        self._collapsed = []
        self._redraw_for_collapse()

    def collapse_handles(self) -> List[Tuple[str, Tuple[str, ...], bool]]:
        """The drawn collapse handles: their group's name ("Cluster 4"; ""
        for an unnamed node), leaves in row order, and whether collapsed."""
        out: List[Tuple[str, Tuple[str, ...], bool]] = []
        for item in self._collapse_handle_items:
            try:
                out.append((item.label, item.members(), item.is_collapsed()))
            except RuntimeError:
                continue
        return out

    def collapsed_rows(self) -> List[Tuple[str, Tuple[str, ...]]]:
        """The drawn collapsed subtrees, top to bottom: their label
        ("Cluster 4 · 83 genomes") and their leaves in row order."""
        position = {name: index for index, name in enumerate(self._leaf_order)}
        rows = sorted(self._collapsed_drawn.values(), key=lambda entry: position.get(entry[0][0], 0))
        return [(label, members) for members, label, _anchor in rows]

    def _redraw_for_collapse(self) -> None:
        # The hovered items are about to be replaced: no hover outlives them.
        if self._active_leaf_set:
            self._hover_subtree_key = None
            self._action_hover_key = None
            self._set_active_leaves(None, emit_signals=True)
        if self._root is not None:
            self._render_tree()

    def _unpruned_root(self) -> _TreeNode:
        """The tree a pruned drawing is cut from (its root and rotation)."""
        if (
            self._root_mode == ROOT_GROUPS
            and not self._reroot_for_groups
            and self._group_full_root is not None
            and self._should_optimize_groups()
        ):
            return self._group_full_root
        return _midpoint_root(self._root) if self._root_mode == ROOT_MIDPOINT else self._root

    def _with_collapsed_grade(self, root: _TreeNode) -> _TreeNode:
        """``root`` with a collapsed clade or cluster that is a grade at the
        root (phylo.root_grades: a root placed inside it) gathered under one
        node, so it is drawn as a triangle."""
        if not self._collapsed:
            return root
        names = frozenset(_leaf_names(root, self._canonical_leaf_name))
        for groups in (self._leaf_group_map, self._leaf_cluster_map):
            for _group, grade in root_grades(root, groups, self._canonical_leaf_name):
                if any((members & names) == grade for members in self._collapsed):
                    gathered = collapse_root_grade(root, grade, self._canonical_leaf_name)
                    if gathered is not None:
                        return gathered
        return root

    def _collapsed_nodes(self, root: _TreeNode) -> Dict[int, frozenset]:
        """The nodes of the drawn tree to collapse (by id). The sets it has
        no subtree for are kept, not drawn (``_collapsed_unmatched``)."""
        self._collapsed_unmatched = set()
        self._collapsed_drawn_of = {}
        if not self._collapsed:
            return {}
        pruned = self._prune_inactive and self._active_samples is not None
        nodes, unmatched, drawn_of = resolve_collapsed(
            root, self._collapsed, self._canonical_leaf_name, pruned=pruned
        )
        self._collapsed_unmatched = set(unmatched)
        self._collapsed_drawn_of = drawn_of
        return nodes

    def _handle_node_ids(self, root: _TreeNode) -> Set[int]:
        """The nodes where a clade or cluster starts (phylo.group_starts);
        a group in several subtrees gets its biggest ones only."""
        ids: Set[int] = set()
        for groups in (self._leaf_group_map, self._leaf_cluster_map):
            by_group: Dict[str, List[Tuple[int, _TreeNode]]] = {}
            for group, node in group_starts(root, groups, self._canonical_leaf_name):
                by_group.setdefault(group, []).append((len(_leaf_names(node, self._canonical_leaf_name)), node))
            for starts in by_group.values():
                if len(starts) > 1:
                    starts = [start for start in starts if start[0] >= _SPLIT_GROUP_MIN_HANDLE_SIZE]
                    starts = sorted(starts, key=lambda start: -start[0])[:_SPLIT_GROUP_MAX_HANDLES]
                ids.update(id(node) for _size, node in starts)
        return ids

    def _grade_handles(self, root: _TreeNode, order: Sequence[str]) -> List[List[str]]:
        """The clades and clusters that are a grade at the root (see
        ``_with_collapsed_grade``), each as its leaves in ``order``: their
        handle sits on the root."""
        out: List[List[str]] = []
        for groups in (self._leaf_group_map, self._leaf_cluster_map):
            for _group, grade in root_grades(root, groups, self._canonical_leaf_name):
                members = [name for name in order if name in grade]
                if members not in out:
                    out.append(members)
        return out

    @staticmethod
    def _group_text(kind: str, name: str) -> str:
        return name if name.lower().startswith(kind.lower()) else f"{kind} {name}"

    def _subtree_name(self, members: Sequence[str], sizes: Mapping[Tuple[str, str], int]) -> Tuple[str, int]:
        """The group a subtree is ("Cluster 4"; "Clade 2"; "Clade 1 ·
        Cluster 1" when it is a whole clade that is one cluster; "" when it
        is neither one clade nor one cluster) and the drawn genomes of the
        last group named (0 for none). ``sizes``: _drawn_group_sizes."""
        clades = {self._leaf_group_map.get(name) for name in members}
        clusters = {self._leaf_cluster_map.get(name) for name in members}
        clade = next(iter(clades)) if len(clades) == 1 else None
        cluster = next(iter(clusters)) if len(clusters) == 1 else None
        parts: List[str] = []
        total = 0
        if clade and (not cluster or sizes.get(("clade", clade)) == len(members)):
            parts.append(self._group_text("Clade", clade))
            total = sizes.get(("clade", clade), 0)
        if cluster:
            parts.append(self._group_text("Cluster", cluster))
            total = sizes.get(("cluster", cluster), 0)
        return " \u00b7 ".join(parts), total

    def _drawn_group_sizes(self, leaf_names: Sequence[str]) -> Dict[Tuple[str, str], int]:
        """How many drawn genomes each clade and cluster has."""
        sizes: Dict[Tuple[str, str], int] = {}
        for name in leaf_names:
            for kind, groups in (("clade", self._leaf_group_map), ("cluster", self._leaf_cluster_map)):
                group = groups.get(name)
                if group:
                    sizes[(kind, group)] = sizes.get((kind, group), 0) + 1
        return sizes

    def _collapsed_text(self, name: str, total: int, members: Sequence[str]) -> str:
        """"Cluster 4 · 83 genomes"; "Cluster 3 · 46 of 60 genomes" for one
        of the subtrees of a group that is not monophyletic."""
        count = len(members)
        if total > count:
            genomes = f"{count} of {total} genomes"
        else:
            genomes = f"{count} genome" if count == 1 else f"{count} genomes"
        return f"{name} \u00b7 {genomes}" if name else genomes

    def _collapsed_tooltip(self, text: str, members: Sequence[str]) -> str:
        lines = [text]
        if self._active_samples is not None:
            inside = sum(1 for name in members if name in self._active_samples)
            lines.append(f"{inside} of {len(members)} in the sample set")
        lines.append("Double-click, or click its handle, to expand")
        return "\n".join(lines)

    def _add_collapsed_clade(
        self,
        path: QtGui.QPainterPath,
        members: Sequence[str],
        sizes: Mapping[Tuple[str, str], int],
        *,
        band: Optional[Tuple[QtCore.QPointF, float, float, float]] = None,
        align_right: bool = False,
    ) -> Tuple[_CollapsedCladeItem, _CollapsedLabelItem]:
        """Add a collapsed subtree's triangle (``path``, or ``band``: see
        _CollapsedCladeItem) and its label (screen stable; the caller
        places it)."""
        color = self._edge_color_for_subtree(members)
        item = _CollapsedCladeItem(path, self, members, color, band)
        name, total = self._subtree_name(members, sizes)
        text = self._collapsed_text(name, total, members)
        tooltip = self._collapsed_tooltip(text, members)
        item.setToolTip(tooltip)
        self._scene.addItem(item)
        font = QtGui.QFont(self.font())
        font.setPointSizeF(_COLLAPSED_LABEL_FONT_SIZE)
        font.setWeight(QtGui.QFont.Weight.DemiBold)
        # Two lines (the group, then its genomes): narrow beside the rows.
        shown = f"{name}\n{text[len(name) + 3:]}" if name else text
        label = _CollapsedLabelItem(shown, font, align_right=align_right)
        label.full_text = text
        label_color = QtGui.QColor(color)
        if self._subtree_active_state(members) != "inactive":
            # Readable: the group's color, never a faded one.
            label_color.setAlpha(255)
        label.setBrush(QtGui.QBrush(label_color))
        label.setToolTip(tooltip)
        self._add_device_stable_item(label)
        self._collapsed_label_items.append(label)
        return item, label

    def _add_collapse_handle(
        self,
        point: QtCore.QPointF,
        members: Sequence[str],
        name: str,
        *,
        collapsed: bool,
    ) -> None:
        handle = _CollapseHandleItem(
            self,
            members,
            self._edge_color_for_subtree(members),
            collapsed=collapsed,
            label=name,
        )
        handle.setPos(point)
        self._add_device_stable_item(handle)
        self._collapse_handle_items.append(handle)

    def node_menu(self, leaves: Sequence[str]) -> Optional[QtWidgets.QMenu]:
        """The right click menu of a node (``leaves``: the leaves under it):
        Collapse or Expand it, and Expand all while anything is collapsed;
        None when it offers nothing (a leaf, the root, no collapse)."""
        members = self._leaf_set(leaves)
        menu = QtWidgets.QMenu(self)
        if self.is_collapsed(members):
            action = menu.addAction("Expand")
            action.setObjectName("TreeExpandAction")
            action.triggered.connect(lambda _c=False, m=tuple(members): self._later(lambda: self.expand_subtree(m)))
        elif len(members) >= 2 and len(members) < len(self._leaf_order):
            action = menu.addAction("Collapse")
            action.setObjectName("TreeCollapseAction")
            action.triggered.connect(lambda _c=False, m=tuple(members): self._later(lambda: self.collapse_subtree(m)))
        if self._collapsed_drawn:
            action = menu.addAction("Expand all")
            action.setObjectName("TreeExpandAllAction")
            action.triggered.connect(lambda _c=False: self._later(self.expand_all))
        if menu.isEmpty():
            menu.deleteLater()
            return None
        return menu

    @staticmethod
    def _later(callback: Callable[[], None]) -> None:
        # A redraw deletes the scene's items: never inside one's own event.
        QtCore.QTimer.singleShot(0, callback)

    def _show_node_menu(self, leaves: Sequence[str], screen_pos: QtCore.QPoint) -> None:
        menu = self.node_menu(leaves)
        if menu is None:
            return
        menu.exec(screen_pos)
        menu.deleteLater()

    def _node_near(self, pos: QtCore.QPoint) -> Tuple[str, ...]:
        """The leaves under the drawn internal node nearest to ``pos``
        (viewport pixels), within _NODE_MENU_RADIUS_PX; else under the
        branch there; empty when neither."""
        transform = self.viewportTransform()
        best: Tuple[str, ...] = ()
        best_distance = _NODE_MENU_RADIUS_PX
        for scene_point, members in self._node_points:
            point = transform.map(scene_point)
            distance = math.hypot(point.x() - pos.x(), point.y() - pos.y())
            if distance <= best_distance:
                best, best_distance = tuple(members), distance
        if best:
            return best
        item = self.itemAt(pos)
        if isinstance(item, _TreeEdgeItem) and len(item.subtree_leaves()) >= 2:
            return item.subtree_leaves()
        return ()

    def contextMenuEvent(self, event: QtGui.QContextMenuEvent) -> None:
        item = self.itemAt(event.pos())
        while item is not None and not isinstance(item, (_ActionTextItem, _CollapsedCladeItem, _CollapseHandleItem)):
            item = item.parentItem()
        if item is None and self._leaf_order:
            # Not a label, triangle or handle (they have their own menus).
            menu = self.node_menu(self._node_near(event.pos()))
            if menu is not None:
                menu.exec(event.globalPos())
                menu.deleteLater()
                event.accept()
                return
        super().contextMenuEvent(event)

    @contextmanager
    def _without_collapse_handles(self) -> Iterator[None]:
        """Exports draw the triangles but not the handles."""
        shown = [item for item in self._collapse_handle_items if item.scene() is self._scene and item.isVisible()]
        for item in shown:
            item.setVisible(False)
        try:
            yield
        finally:
            for item in shown:
                item.setVisible(True)

    def export_figure(
        self,
        path: Path,
        *,
        image_scale: float = 4.0,
        image_dpi: int = 300,
        caption: Optional[str] = None,
    ) -> None:
        """Save the drawn tree as PDF, PNG or JPEG; ``caption`` (for example
        the tree source and its root) is printed above it."""
        target_path = Path(path).expanduser()
        suffix = target_path.suffix.lower()
        if suffix not in {".pdf", ".png", ".jpg", ".jpeg"}:
            raise ValueError(f"Unsupported export format: {target_path.suffix or '<none>'}")

        # Screen annotations ignore the view transform so they remain legible
        # while zooming. Temporarily restore normal scene scaling for export so
        # PNG resolution and PDF vectors scale as expected. A hover's emphasis
        # (thick path, dimmed tree) is not part of the figure.
        with self._scene_scaled_annotations_for_export(), self._without_hover_emphasis(), \
                self._without_collapse_handles(), self._export_caption(caption):
            scene_rect = self._scene.itemsBoundingRect().adjusted(-20.0, -20.0, 20.0, 20.0)
            if scene_rect.width() <= 1.0 or scene_rect.height() <= 1.0:
                raise ValueError("Tree scene is empty.")

            if suffix == ".pdf":
                self._export_pdf(target_path, scene_rect, dpi=image_dpi)
                return
            self._export_image(target_path, scene_rect, scale=image_scale, dpi=image_dpi)

    @contextmanager
    def _export_caption(self, caption: Optional[str]) -> Iterator[None]:
        """A caption line above the tree while it is exported."""
        text = str(caption or "").strip()
        bounds = self._scene.itemsBoundingRect()
        if not text or bounds.isEmpty():
            yield
            return
        item = QtWidgets.QGraphicsSimpleTextItem(text)
        item.setBrush(QtGui.QColor(_TREE_LABEL_COLOR))
        item.setPos(bounds.left(), bounds.top() - item.boundingRect().height() - 8.0)
        self._scene.addItem(item)
        try:
            yield
        finally:
            self._scene.removeItem(item)

    @contextmanager
    def _without_hover_emphasis(self) -> Iterator[None]:
        active = list(self._active_leaf_set)
        if active:
            self._set_active_leaves(None, emit_signals=False)
        try:
            yield
        finally:
            if active:
                self._set_active_leaves(active, emit_signals=False)

    @contextmanager
    def _scene_scaled_annotations_for_export(self) -> Iterator[None]:
        ignore_flag = QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations
        registered = [
            item
            for item in self._device_stable_items
            if item is not None and item.scene() is self._scene
        ]
        previous = [(item, bool(item.flags() & ignore_flag), float(item.scale())) for item in registered]
        # Chips hidden on screen for overlapping (_declutter_group_chips) are
        # all exported: scaled with the scene they do not overlap.
        hidden = [item for item in registered if not item.isVisible()]
        # The figure is the whole tree at its own size: no camera zoom on the
        # annotations, no squeeze of the leaf labels.
        squeeze = self._label_squeeze
        try:
            for item, ignored, _scale in previous:
                if ignored:
                    item.setFlag(ignore_flag, False)
                item.setScale(1.0)
            for item in hidden:
                item.setVisible(True)
            if squeeze is not None:
                self._squeeze_leaf_labels(1.0, 1.0)
            yield
        finally:
            for item, ignored, scale in previous:
                item.setFlag(ignore_flag, ignored)
                item.setScale(scale)
            for item in hidden:
                item.setVisible(False)
            if squeeze is not None:
                self._squeeze_leaf_labels(*squeeze)

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
            "treeExport",
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

    def visible_row_range(self) -> Optional[Tuple[int, int]]:
        """The drawn rows in view (first, last; zero based, inclusive, a row
        counts when its center is in the viewport), rectangular layout only;
        None otherwise or when no row is in view."""
        if self._layout_mode != "rectangular" or not self._leaf_order:
            return None
        viewport = self.viewport()
        if viewport is None or viewport.height() <= 0:
            return None
        total = len(self._leaf_order)
        first = self.leaf_viewport_y(self._leaf_order[0])
        if first is None:
            return None
        if total > 1:
            last = self.leaf_viewport_y(self._leaf_order[-1])
            if last is None:
                return None
            pitch = (float(last) - float(first)) / float(total - 1)
        else:
            pitch = float(viewport.height())
        if not math.isfinite(pitch) or pitch <= 0.0:
            return None
        if self._collapsed_drawn:
            # Rows of unequal height: the leaves placed in the viewport
            # (a collapsed subtree's along its triangle).
            transform = self.viewportTransform()
            height = float(viewport.height())
            inside = [
                index
                for index, name in enumerate(self._leaf_order)
                if 0.0 <= transform.map(self._leaf_pos_scene[name]).y() <= height
            ]
            return (inside[0], inside[-1]) if inside else None
        # Row i spans i to i + 1 with its center (the leaf) at i + 0.5.
        top_row = 0.5 - float(first) / pitch
        return visible_rows(top_row, float(viewport.height()) / pitch, total)

    def row_range_text(self) -> str:
        """'Rows 1 to 60 of 153' while part of the rows is scrolled away ('':
        every row in view, or no rows)."""
        if self._layout_mode != "rectangular" or not self._leaf_order:
            return ""
        return row_range_text(self.visible_row_range(), len(self._leaf_order))

    def set_row_note_enabled(self, enabled: bool) -> None:
        """Show the row range note over the tree while it shows only part of
        its rows."""
        if enabled and self._row_note is None:
            self._row_note = RowRangeNote(self, self.viewport())
        elif not enabled and self._row_note is not None:
            self._row_note.deleteLater()
            self._row_note = None
        self._update_row_note()

    def row_note(self) -> Optional[RowRangeNote]:
        return self._row_note

    def _update_row_note(self) -> None:
        if self._row_note is not None:
            self._row_note.set_note(self.row_range_text())

    def leaf_viewport_y(self, sample: str) -> Optional[float]:
        key = self._normalize_label(sample)
        pt = self._leaf_pos_scene.get(key)
        if pt is None:
            return None
        # Sub-pixel (mapFromScene rounds to whole pixels).
        return float(self.viewportTransform().map(pt).y())

    def highlight_sample(self, sample: Optional[str]) -> None:
        key = self._normalize_label(sample) if sample else None
        if key and key not in self._leaf_pos_scene:
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
        self._lengths_root = None
        self._lengths_display = None
        self._set_branch_distance_status(False, "Distance unavailable: no tree is loaded.")
        self._invalidate_render_root_cache()
        self._draw_error(message)

    def _invalidate_render_root_cache(self) -> None:
        self._render_root_cache_key = None
        self._render_root_cache = None
        self._base_root_cache_key = None
        self._base_root_cache = None
        self._group_full_root_key = None
        self._group_full_root = None
        self._ordered_root_key = None
        self._ordered_root = None
        self._full_order_key = None
        self._full_order = None

    def _reset_collapse_drawing(self) -> None:
        """Forget the drawn collapse items (before a redraw)."""
        # The drawn collapsed subtrees: members in row order, label and the
        # scene point a hover marks (the middle of the triangle's far side).
        self._collapsed_drawn: Dict[frozenset, Tuple[Tuple[str, ...], str, QtCore.QPointF]] = {}
        self._collapsed_of_leaf: Dict[str, frozenset] = {}
        self._collapse_handle_items: List[_CollapseHandleItem] = []
        self._collapsed_label_items: List[QtWidgets.QGraphicsSimpleTextItem] = []
        # The drawn internal nodes (not the root) a right click can collapse.
        self._node_points: List[Tuple[QtCore.QPointF, Sequence[str]]] = []

    def _draw_error(self, message: str) -> None:
        self._reset_collapse_drawing()
        self._content_item = None
        self._hover_dimmed = False
        self._scene.clear()
        self._device_stable_items = []
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
        self._branch_breaks = BranchBreaks()
        self._branch_breaks_drawn = False
        self._branch_break_notes = {}
        self._announce_branch_breaks()

        txt = QtWidgets.QGraphicsSimpleTextItem(message)
        txt.setBrush(QtGui.QBrush(QtGui.QColor("#64748b")))
        txt.setPos(8, 8)
        self._scene.addItem(txt)
        rect = QtCore.QRectF(0, 0, max(300.0, txt.boundingRect().width() + 24.0), 120.0)
        self._scene.setSceneRect(rect)
        self.resetTransform()
        self._camera_zoom = 1.0
        # The next tree is fitted (no camera to keep over a message).
        self._auto_fit = True
        self._fit_rect_cache = None
        self._declutter_key = None
        self.centerOn(rect.center())
        self.viewChanged.emit()

    def _empty_prune_message(self) -> str:
        count = len(self._active_samples or ())
        if count == 0:
            return "No genomes match the current filters."
        if count == 1:
            return "The filtered genome is not in the tree."
        return f"None of the {count} filtered genomes are in the tree."

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

    def _leaf_font(self, leaf_name: str) -> QtGui.QFont:
        font = QtGui.QFont(self.font())
        font.setPointSizeF(_LEAF_FONT_SIZE)
        if self._reference_sample and leaf_name == self._reference_sample:
            font.setWeight(QtGui.QFont.Weight.DemiBold)
        else:
            font.setWeight(QtGui.QFont.Weight.Normal)
        return font

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
            return QtGui.QColor(_TREE_BRANCH_COLOR)
        color = self._group_color_map.get(group_name)
        return QtGui.QColor(color) if color is not None else QtGui.QColor(_TREE_BRANCH_COLOR)

    def _display_group(self, leaf_name: str) -> Optional[str]:
        group = self._leaf_display_group_map.get(leaf_name)
        return group if group else None

    def _display_group_color(self, group_name: Optional[str]) -> QtGui.QColor:
        if not group_name:
            return QtGui.QColor(_TREE_BRANCH_COLOR)
        color = self._display_group_color_map.get(group_name)
        return QtGui.QColor(color) if color is not None else QtGui.QColor(_TREE_BRANCH_COLOR)

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
    ) -> QtGui.QColor:
        """The color of what is drawn for ``subtree_leaves``: its own while
        any of them is in the sample set (a branch on the way to an in set
        genome keeps its color), the faint out of set grey when none is."""
        if self._subtree_active_state(subtree_leaves) == "inactive":
            return out_of_set_color()
        return QtGui.QColor(active_color)

    def _leaf_base_color(self, leaf_name: str) -> QtGui.QColor:
        display_group = self._display_group(leaf_name) if self._display_groups_color_labels else None
        if display_group:
            return self._display_group_color(display_group)
        if self._use_group_backgrounds:
            return QtGui.QColor(_TREE_LABEL_COLOR)
        return self._group_color(self._leaf_group(leaf_name))

    def _leaf_cluster(self, leaf_name: str) -> Optional[str]:
        cluster = self._leaf_cluster_map.get(leaf_name)
        return cluster if cluster else None

    @staticmethod
    def _compact_cluster_label(cluster_name: str) -> str:
        value = str(cluster_name or "").strip()
        lowered = value.lower()
        if lowered.startswith("cluster"):
            value = value[len("cluster"):].strip(" :-_")
        return f"C{value}" if value else "Cluster"

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
            return self._visual_color_for_subtree(subtree_leaves, QtGui.QColor(_TREE_BRANCH_COLOR))
        group = self._subtree_group(subtree_leaves)
        if group:
            return self._visual_color_for_subtree(subtree_leaves, self._group_color(group))
        return self._visual_color_for_subtree(subtree_leaves, QtGui.QColor(_TREE_BRANCH_COLOR))

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
            color = QtGui.QColor(_TREE_INACTIVE_COLOR)
            color.setAlpha(_CLADE_INACTIVE_TINT_ALPHA)
        elif state == "partial":
            color.setAlpha(_CLADE_PARTIAL_TINT_ALPHA)
        else:
            color.setAlpha(_CLADE_TINT_ALPHA)
        return QtGui.QBrush(color)

    def _edge_pen_for_color(self, color: QtGui.QColor, *, width: float = 1.0) -> QtGui.QPen:
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

    def _group_for_leaf_sequence(self, leaf_name: str) -> str:
        group = self._leaf_group(leaf_name)
        return group if group is not None else ""

    def _optimize_tree_for_groups(self, node: _TreeNode) -> GroupSummary:
        """Rotate children so leaf groups run contiguously (phylo)."""
        return GroupOptimizer(self._group_for_leaf_sequence, self._canonical_leaf_name).optimize(node)

    def _graph_node_name(self, node: _TreeNode) -> str:
        """A rerooted node's name: canonical for a leaf, unquoted otherwise."""
        if node.is_leaf:
            return self._canonical_leaf_name(node.name or "")
        return self._normalize_label(node.name or "")

    def _base_render_root(self) -> Optional[_TreeNode]:
        """The tree to draw before any group rotation: the file tree, rooted
        at its midpoint when asked, then pruned to the active samples when
        inactive ones are hidden (the root stays where the full tree has it,
        so a filter never moves it). None when pruning keeps no leaf."""
        if self._root is None:
            raise RuntimeError("No tree loaded.")
        prune_key = (
            frozenset(self._active_samples)
            if self._prune_inactive and self._active_samples is not None
            else None
        )
        if self._root_mode in (ROOT_AS_IN_FILE, ROOT_GROUPS) and prune_key is None:
            # The render cache keys on the base identity: reset it, or a
            # cleared filter would match the last pruned tree's cache.
            self._base_root_cache_key = None
            self._base_root_cache = None
            return self._root
        cache_key = (self._root_version, self._root_mode, prune_key)
        if self._base_root_cache_key == cache_key:
            return self._base_root_cache
        root = _midpoint_root(self._root) if self._root_mode == ROOT_MIDPOINT else self._root
        if prune_key is not None:
            root = _prune_to(root, set(prune_key), self._canonical_leaf_name)
        self._base_root_cache_key = cache_key
        self._base_root_cache = root
        return root

    def _resolved_render_root(self) -> Optional[_TreeNode]:
        """The tree drawn: the grouped root (``_grouped_render_root``) in
        the order mode's leaf order. Ladderizing keeps every clade, cluster
        and branch color band as contiguous as the grouped root has it; a
        pruned tree follows the full tree's ladder order when it can, so a
        filter only removes rows and never reorders them."""
        grouped = self._grouped_render_root()
        if grouped is None or self._order_mode == ORDER_AS_COMPUTED:
            return grouped
        key = (self._order_mode, self._reference_sample, self._band_version)
        if self._ordered_root_key is not None and self._ordered_root_key[0] is grouped and self._ordered_root_key[1:] == key:
            return self._ordered_root
        bands = [band for band in (self._leaf_group_map, self._leaf_cluster_map, self._leaf_display_group_map) if band]
        ordered = None
        full_order = self._full_ladder_order(bands) if self._prune_inactive and self._active_samples is not None else None
        if full_order is not None:
            ordered = follow_order(grouped, full_order, leaf_name=self._canonical_leaf_name, keep_together=bands)
        if ordered is None:
            ordered = ladderize(grouped, leaf_name=self._canonical_leaf_name, mode=self._order_mode, keep_together=bands)
        self._ordered_root = ordered
        self._ordered_root_key = (grouped, *key)
        return ordered

    def _full_ladder_order(self, bands: Sequence[Mapping[str, str]]) -> Optional[List[str]]:
        """The leaf order of the unpruned tree, grouped and ladderized as it
        is drawn without a filter; None when the pruned tree is rerooted on
        its own (``reroot_for_groups``), as its root moves with the filter."""
        optimize = self._should_optimize_groups()
        if optimize and (self._reroot_for_groups or (self._root_mode == ROOT_GROUPS and self._group_full_root is None)):
            return None
        key = (self._root_version, self._root_mode, optimize, self._order_mode, self._reference_sample, self._band_version)
        if self._full_order_key == key:
            return self._full_order
        if optimize and self._root_mode == ROOT_GROUPS:
            full = self._group_full_root
        else:
            full = _midpoint_root(self._root) if self._root_mode == ROOT_MIDPOINT else self._root
            if optimize:
                full = _clone_tree(full)
                self._optimize_tree_for_groups(full)
        ladder = ladderize(full, leaf_name=self._canonical_leaf_name, mode=self._order_mode, keep_together=bands)
        self._full_order = _leaf_names(ladder, self._canonical_leaf_name)
        self._full_order_key = key
        return self._full_order

    def _grouped_render_root(self) -> Optional[_TreeNode]:
        base_root = self._base_render_root()
        if base_root is None or not self._should_optimize_groups():
            return base_root

        cache_key = (
            self._root_version,
            self._base_root_cache_key,
            self._root_mode,
            self._optimize_groups,
            self._reroot_for_groups,
            self._group_signature(),
        )
        if self._render_root_cache_key == cache_key and self._render_root_cache is not None:
            return self._render_root_cache

        pruned = self._prune_inactive and self._active_samples is not None
        if self._root_mode == ROOT_GROUPS and not self._reroot_for_groups and pruned:
            # The group root and rotation are chosen on the full tree, then pruned, so a filter
            # never moves the root; rotated again only if the pruned clades
            # are no longer contiguous.
            full_key = (self._root_version, self._group_signature())
            if self._group_full_root_key != full_key or self._group_full_root is None:
                self._group_full_root = group_root(
                    self._root,
                    self._group_for_leaf_sequence,
                    self._canonical_leaf_name,
                    node_name=self._graph_node_name,
                )
                self._group_full_root_key = full_key
            best_root = _prune_to(self._group_full_root, set(self._active_samples or ()), self._canonical_leaf_name)
            if best_root is not None:
                names = _leaf_names(best_root, self._canonical_leaf_name)
                group_map = {name: self._group_for_leaf_sequence(name) for name in names}
                if not groups_contiguous(best_root, group_map, self._canonical_leaf_name):
                    self._optimize_tree_for_groups(best_root)
        elif self._reroot_for_groups or self._root_mode == ROOT_GROUPS:
            # The reroot search over every edge (phylo.group_root).
            best_root = group_root(
                base_root,
                self._group_for_leaf_sequence,
                self._canonical_leaf_name,
                node_name=self._graph_node_name,
            )
        else:
            best_root = _clone_tree(base_root)
            self._optimize_tree_for_groups(best_root)

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
        drawn_lengths: Optional[Mapping[int, float]] = None,
    ) -> float:
        """Each node's x in branch units (``drawn_lengths``: the length an
        edge is drawn with instead of its own, by the child's id)."""
        distance_map[id(node)] = cumulative
        max_distance = cumulative
        for child in node.children:
            branch_length = child.branch_length
            if branch_length is None or not math.isfinite(branch_length) or branch_length < 0.0:
                branch_length = 0.0
            if drawn_lengths is not None and id(child) in drawn_lengths:
                branch_length = drawn_lengths[id(child)]
            max_distance = max(
                max_distance,
                self._collect_branch_distances(child, cumulative + float(branch_length), distance_map, drawn_lengths),
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
        collapsed: Optional[Mapping[int, object]] = None,
    ) -> _LayoutNode:
        """The layout tree of ``node``; a ``collapsed`` node (by id) is a
        layout leaf (its subtree is drawn as one triangle)."""
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
        children = node.children if collapsed is None or id(node) not in collapsed else []
        for idx, child in enumerate(children, start=1):
            child_layout = self._build_layout_tree(
                child,
                parent=layout,
                depth=depth + 1,
                number=idx,
                left_sibling=prev_child,
                leftmost_sibling=first_child,
                node_map=node_map,
                collapsed=collapsed,
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

        # Preserve logical highlight and manual viewport state across redraw, but
        # never reuse old QGraphicsItem pointers because scene.clear() deletes them.
        pending_transform: Optional[QtGui.QTransform] = None
        pending_center: Optional[QtCore.QPointF] = None
        # The rows kept at their height on screen: (leaves, scene y).
        pending_rows: Optional[Tuple[Sequence[str], float]] = None
        viewport = self.viewport()
        if not self._auto_fit and viewport is not None:
            pending_transform = QtGui.QTransform(self.transform())
            pending_center = self._scene_point_at(QtCore.QPointF(viewport.width() / 2.0, viewport.height() / 2.0))
            pending_rows = self._anchor_rows(pending_center)
        self._redraw_anchor = None
        self._restretch_timer.stop()
        self._fit_rect_cache = None
        self._declutter_key = None

        pending_active = list(self._active_leaf_set)
        pending_subtree = self._hover_subtree_key
        self._reset_collapse_drawing()
        self._device_stable_items = []
        self._label_squeeze = None
        self._zoom_anchor = None
        self._content_item = None
        self._hover_dimmed = False
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
        self._branch_breaks = BranchBreaks()
        self._branch_breaks_drawn = False
        self._branch_break_notes = {}

        try:
            render_root = self._resolved_render_root()
        except Exception as exc:
            self._draw_error(f"Failed to render tree: {exc}")
            return
        if render_root is None:
            self._draw_error(self._empty_prune_message())
            return
        # The order, root and clades are chosen on the loaded tree; only the
        # drawn lengths come from another fit.
        render_root = self._with_collapsed_grade(self._with_branch_lengths(render_root))

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
        self._announce_branch_breaks()
        if self._leaf_order:
            self._adopt_tree_items()

        if pending_subtree:
            self._set_active_leaves(list(pending_subtree), emit_signals=False)
            self._hover_subtree_key = tuple(self._ordered_leaf_subset(pending_subtree))
        elif pending_active:
            self._set_active_leaves(pending_active, emit_signals=False)

        if self._auto_fit:
            self._fit_to_viewport()
        elif pending_transform is not None and pending_center is not None:
            # The camera stays where it was (the same zoom, stretch and
            # center); rows added or taken out above move the center so the
            # anchor rows stay where they were on screen.
            self.setTransform(pending_transform)
            self._apply_camera_to_items()
            self._update_camera_scene_rect()
            if pending_rows is not None:
                now = self._rows_scene_y(pending_rows[0])
                if now is not None:
                    pending_center = QtCore.QPointF(pending_center.x(), pending_center.y() + now - pending_rows[1])
            middle = QtCore.QPointF(viewport.width() / 2.0, viewport.height() / 2.0)
            self._put_scene_point_at(pending_center, middle)

        self.viewChanged.emit()

    def _render_tree_if_drawn(self) -> None:
        if self._root is not None and self._leaf_order:
            self._render_tree()

    def _rows_scene_y(self, leaves: Sequence[str]) -> Optional[float]:
        """The scene y of the middle of these leaves' rows (a collapsed one
        at its triangle); None: none is drawn, or the layout is circular."""
        if self._layout_mode != "rectangular":
            return None
        ys: List[float] = []
        for name in leaves:
            point = self._leaf_pos_scene.get(name)
            if point is None:
                key = self._collapsed_of_leaf.get(name)
                drawn = self._collapsed_drawn.get(key) if key is not None else None
                point = drawn[2] if drawn is not None else None
            if point is not None:
                ys.append(float(point.y()))
        if not ys:
            return None
        return (min(ys) + max(ys)) / 2.0

    def _anchor_rows(self, center: QtCore.QPointF) -> Optional[Tuple[Sequence[str], float]]:
        """The rows a redraw keeps on screen: the subtree being collapsed or
        expanded, else the row nearest ``center``; with their scene y."""
        if self._layout_mode != "rectangular" or not self._leaf_pos_scene:
            return None
        leaves: Sequence[str]
        if self._redraw_anchor:
            leaves = tuple(self._redraw_anchor)
        else:
            nearest = min(self._leaf_pos_scene.items(), key=lambda entry: abs(float(entry[1].y()) - center.y()))
            leaves = (nearest[0],)
        y = self._rows_scene_y(leaves)
        return (leaves, y) if y is not None else None

    def _render_tree_circular(
        self,
        root: _TreeNode,
        depth_map: Dict[int, int],
        leaves: List[_TreeNode],
        leaf_name_map: Dict[int, str],
    ) -> None:
        collapsed = self._collapsed_nodes(root)
        handle_ids = self._handle_node_ids(root)
        # The rows around the circle: a leaf, or a collapsed subtree
        # (_collapsed_row_count rows wide).
        rows: List[_TreeNode] = []
        row_width: Dict[int, float] = {}
        stack = [root]
        while stack:
            node = stack.pop()
            if node.is_leaf or id(node) in collapsed:
                rows.append(node)
                row_width[id(node)] = (
                    _collapsed_row_count(sum(1 for _leaf in self._iter_leaves((node,))))
                    if id(node) in collapsed
                    else 1.0
                )
            else:
                stack.extend(reversed(node.children))
        max_depth = max(depth_map.values()) if depth_map else 0
        slots = sum(row_width.values())
        step = (2.0 * math.pi) / float(max(1.0, slots))
        start_angle = -math.pi / 2.0
        row_angles: Dict[int, float] = {}
        if len(rows) == 1 and id(rows[0]) not in collapsed:
            row_angles[id(rows[0])] = start_angle
        else:
            cursor = 0.0
            for row in rows:
                width = row_width[id(row)]
                row_angles[id(row)] = start_angle + (cursor + (width - 1.0) / 2.0) * step
                cursor += width

        angle_map: Dict[int, float] = {}

        def assign_angle(node: _TreeNode) -> float:
            if id(node) in row_angles:
                ang = row_angles[id(node)]
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
        group_sizes = self._drawn_group_sizes(list(leaf_name_map.values()))
        outer_r = inner_r + float(max_depth) * depth_step
        label_r = outer_r + 28.0
        half_angle_of = {node_id: (row_width[node_id] - _COLLAPSED_ROW_GAP) * step / 2.0 for node_id in collapsed}

        # Each collapsed subtree's wedge: from its node to the radius of its
        # deepest leaf; its leaves sit along the far side.
        triangle_by_node: Dict[int, _CollapsedCladeItem] = {}
        label_by_node: Dict[int, _CollapsedLabelItem] = {}
        far_r_by_node: Dict[int, float] = {}
        for node_id, _members_set in collapsed.items():
            node = next(row for row in rows if id(row) == node_id)
            ang = angle_map[node_id]
            half_angle = half_angle_of[node_id]
            apex = point(radius_map[node_id], ang)
            far_r = max((radius_map[id(leaf)] for leaf in self._iter_leaves((node,))), default=radius_map[node_id])
            far_r = max(far_r, radius_map[node_id] + depth_step * 0.5)
            far_r_by_node[node_id] = far_r
            members = [name for name in subtree_leaf_map.get(node_id, []) if name]
            path = QtGui.QPainterPath(apex)
            path.lineTo(point(far_r, ang - half_angle))
            path.lineTo(point(far_r, ang + half_angle))
            path.closeSubpath()
            triangle, label = self._add_collapsed_clade(
                path, members, group_sizes, align_right=math.cos(ang) < 0.0
            )
            triangle_by_node[node_id] = triangle
            label_by_node[node_id] = label
            key = frozenset(members)
            self._collapsed_drawn[key] = (tuple(members), label.full_text, point(far_r, ang))
            for name_ in members:
                self._collapsed_of_leaf[name_] = key

        def draw_edges(node: _TreeNode, path_stack: List[QtWidgets.QGraphicsItem]) -> None:
            if id(node) in collapsed:
                return
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
                triangle = triangle_by_node.get(id(child))
                if triangle is not None:
                    for name in subtree_leaves:
                        self._edge_items_by_leaf[name] = stack + [triangle]
                elif child.is_leaf:
                    leaf_name = leaf_name_map[id(child)]
                    self._edge_items_by_leaf[leaf_name] = stack
                draw_edges(child, stack)

        draw_edges(root, [])

        # The internal nodes a right click finds, and the collapse handles
        # (a clade or cluster that is a grade at the root has its on the root).
        root_point = point(radius_map[id(root)], angle_map[id(root)])
        for members in self._grade_handles(root, subtree_leaf_map.get(id(root), [])):
            self._add_collapse_handle(
                root_point, members, self._subtree_name(members, group_sizes)[0], collapsed=False
            )
        stack = list(root.children)
        while stack:
            node = stack.pop()
            if node.is_leaf:
                continue
            members = [name for name in subtree_leaf_map.get(id(node), []) if name]
            node_point = point(radius_map[id(node)], angle_map[id(node)])
            self._node_points.append((node_point, members))
            if id(node) in handle_ids or id(node) in collapsed:
                self._add_collapse_handle(
                    node_point,
                    members,
                    self._subtree_name(members, group_sizes)[0],
                    collapsed=id(node) in collapsed,
                )
            if id(node) not in collapsed:
                stack.extend(node.children)

        for row in rows:
            row_angle = angle_map[id(row)]
            label = label_by_node.get(id(row))
            if label is not None:
                members = [name for name in subtree_leaf_map.get(id(row), []) if name]
                count = len(members)
                half_angle = half_angle_of[id(row)]
                for index, name in enumerate(members):
                    share = 0.5 if count == 1 else index / float(count - 1)
                    self._leaf_order.append(name)
                    self._leaf_pos_scene[name] = point(far_r_by_node[id(row)], row_angle - half_angle + share * 2.0 * half_angle)
                label_point = point(label_r, row_angle)
                label.setPos(label_point.x() + (5.0 if math.cos(row_angle) >= 0.0 else -5.0), label_point.y())
                continue
            leaf_name = leaf_name_map[id(row)]
            leaf_point = point(radius_map[id(row)], row_angle)
            self._leaf_order.append(leaf_name)
            self._leaf_pos_scene[leaf_name] = leaf_point

            label_text = self._display_leaf_name(leaf_name)
            label_item = _LeafLabelItem(label_text, leaf_name, self)
            label_item.setFont(self._leaf_font(leaf_name))
            label_item.setBrush(QtGui.QBrush(self._label_color(leaf_name)))
            label_item.setToolTip(self._leaf_tooltip(leaf_name, label_text))
            br = label_item.boundingRect()
            label_point = point(label_r, row_angle)
            if math.cos(row_angle) >= 0.0:
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
        collapsed = self._collapsed_nodes(root)
        handle_ids = self._handle_node_ids(root)
        layout_map: Dict[int, _LayoutNode] = {}
        layout_root = self._build_layout_tree(root, node_map=layout_map, collapsed=collapsed)
        self._first_walk(layout_root, distance=1.0)
        min_x = self._second_walk(layout_root)
        if min_x < 0.0:
            self._third_walk(layout_root, -min_x)

        # The rows: a leaf, or a collapsed subtree (_collapsed_row_count rows).
        layout_leaves: List[_LayoutNode] = []
        self._collect_layout_leaves(layout_root, layout_leaves)
        if not layout_leaves:
            self._draw_error("Tree has no leaves.")
            return

        row_width: Dict[int, float] = {
            id(row.node): (
                _collapsed_row_count(sum(1 for _leaf in self._iter_leaves((row.node,))))
                if id(row.node) in collapsed
                else 1.0
            )
            for row in layout_leaves
        }
        n_rows = sum(row_width.values())
        # The whole tree's depth, collapsed subtrees too (their triangles
        # reach their deepest leaf).
        max_depth = max(depth_map.values(), default=0)

        left_margin = 28.0
        top_margin = 24.0

        # Keep labels readable and vertically ordered.
        leaf_step = 34.0
        if n_rows > 100:
            leaf_step = 30.0
        if n_rows > 180:
            leaf_step = 26.0
        if n_rows > 300:
            leaf_step = 22.0
        if n_rows > 500:
            leaf_step = 18.0
        label_min_step = float(QtGui.QFontMetrics(self._leaf_font("")).height()) + 4.0
        if leaf_step < label_min_step:
            leaf_step = label_min_step

        tree_height = max(220.0, float(max(1.0, n_rows - 1.0)) * leaf_step)
        target_width = max(260.0, tree_height * 0.38)
        depth_step = target_width / float(max(1, max_depth))
        depth_step = max(18.0, min(56.0, depth_step))

        # The long branch rule (Distance only): the edges it breaks are drawn
        # with a short length unless true lengths are asked for.
        drawn_lengths: Optional[Dict[int, float]] = None
        if self._preserve_branch_lengths and self._branch_distance_available:
            self._branch_breaks = long_branch_breaks(root, LONG_BRANCH_BREAK_FACTOR)
            if self._branch_breaks and not self._show_true_branch_lengths:
                drawn_lengths = {id(node): self._branch_breaks.drawn_length for node in self._branch_breaks.nodes}
        distance_map: Dict[int, float] = {}
        max_distance = self._collect_branch_distances(root, 0.0, distance_map, drawn_lengths)
        use_branch_lengths = (
            self._preserve_branch_lengths
            and self._branch_distance_available
            and max_distance > 0.0
        )
        self._branch_breaks_drawn = use_branch_lengths and drawn_lengths is not None
        branch_scale = target_width / float(max_distance) if use_branch_lengths else 0.0

        # Enforce strict vertical row ranking: row 0 above row 1, etc.; a
        # collapsed subtree is centered in its rows.
        y_by_layout: Dict[int, float] = {}
        cursor = top_margin
        for leaf_layout in layout_leaves:
            if id(leaf_layout.node) in collapsed:
                width = row_width[id(leaf_layout.node)]
                y_by_layout[id(leaf_layout)] = cursor + (width - 1.0) * leaf_step / 2.0
                cursor += width * leaf_step
            else:
                y_by_layout[id(leaf_layout)] = cursor
                cursor += leaf_step
        last_row_y = cursor - leaf_step

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

        def node_x(node: _TreeNode) -> float:
            if use_branch_lengths:
                return left_margin + float(distance_map.get(id(node), 0.0)) * branch_scale
            return left_margin + float(depth_map.get(id(node), 0)) * depth_step

        pos_map: Dict[int, QtCore.QPointF] = {}
        for layout_node in self._iter_layout_nodes(layout_root):
            pos_map[id(layout_node.node)] = QtCore.QPointF(node_x(layout_node.node), y_by_layout[id(layout_node)])

        subtree_leaf_map: Dict[int, List[str]] = {}
        self._collect_subtree_leaf_names(root, leaf_name_map, subtree_leaf_map)
        group_sizes = self._drawn_group_sizes(list(leaf_name_map.values()))
        if self._branch_breaks_drawn:
            # Set before the edges are made: a broken edge's hover says it.
            unit = self._branch_length_unit or "branch-length units"
            for node, length in zip(self._branch_breaks.nodes, self._branch_breaks.true_lengths):
                key = tuple(name for name in subtree_leaf_map.get(id(node), []) if str(name).strip())
                self._branch_break_notes[key] = (
                    f"Branch shortened for display; true length {self._format_break_length(length)} {unit}."
                )

        # Each collapsed subtree's triangle: from its node to its deepest
        # leaf, its rows tall less _COLLAPSED_ROW_GAP; its leaves sit along
        # the far side, top to bottom (where a hover from another panel
        # marks them).
        triangle_by_node: Dict[int, _CollapsedCladeItem] = {}
        label_by_node: Dict[int, _CollapsedLabelItem] = {}
        far_x_by_node: Dict[int, float] = {}
        half_height_of: Dict[int, float] = {}
        for node_id, members_set in collapsed.items():
            layout_node = layout_map[node_id]
            node = layout_node.node
            apex = pos_map[node_id]
            far_x = max((node_x(leaf) for leaf in self._iter_leaves((node,))), default=apex.x())
            far_x = max(far_x, apex.x() + depth_step)
            far_x_by_node[node_id] = far_x
            width = row_width[node_id]
            half_height = (width - _COLLAPSED_ROW_GAP) * leaf_step / 2.0
            half_height_of[node_id] = half_height
            members = [name for name in subtree_leaf_map.get(node_id, []) if name]
            path = QtGui.QPainterPath(apex)
            path.lineTo(far_x, apex.y() - half_height)
            path.lineTo(far_x, apex.y() + half_height)
            path.closeSubpath()
            band = (QtCore.QPointF(apex), far_x, half_height, (width - 0.15) * leaf_step / 2.0)
            triangle, label = self._add_collapsed_clade(path, members, group_sizes, band=band)
            triangle_by_node[node_id] = triangle
            label_by_node[node_id] = label
            key = frozenset(members)
            anchor = QtCore.QPointF(far_x, apex.y())
            self._collapsed_drawn[key] = (tuple(members), label.full_text, anchor)
            for name_ in members:
                self._collapsed_of_leaf[name_] = key

        def draw_edges(layout_node: _LayoutNode, path_stack: List[QtWidgets.QGraphicsItem]) -> None:
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
                triangle = triangle_by_node.get(id(child_layout.node))
                if triangle is not None:
                    for name in subtree_leaves:
                        self._edge_items_by_leaf[name] = stack + [triangle]
                elif child_layout.is_leaf:
                    leaf_name = leaf_name_map[id(child_layout.node)]
                    self._edge_items_by_leaf[leaf_name] = stack
                draw_edges(child_layout, stack)

        draw_edges(layout_root, [])
        # The internal nodes a right click finds, and the collapse handles
        # where a clade or cluster starts (and on every collapsed node; a
        # clade or cluster that is a grade at the root has its on the root).
        for members in self._grade_handles(root, subtree_leaf_map.get(id(root), [])):
            self._add_collapse_handle(
                pos_map[id(root)], members, self._subtree_name(members, group_sizes)[0], collapsed=False
            )
        for layout_node in self._iter_layout_nodes(layout_root):
            node = layout_node.node
            if layout_node.parent is None or node.is_leaf:
                continue
            members = [name for name in subtree_leaf_map.get(id(node), []) if name]
            self._node_points.append((QtCore.QPointF(pos_map[id(node)]), members))
            if id(node) in handle_ids or id(node) in collapsed:
                self._add_collapse_handle(
                    pos_map[id(node)],
                    members,
                    self._subtree_name(members, group_sizes)[0],
                    collapsed=id(node) in collapsed,
                )
        if self._branch_breaks_drawn:
            short_unit = self._short_unit()
            marks = []
            for node, length in zip(self._branch_breaks.nodes, self._branch_breaks.true_lengths):
                layout_node = layout_map.get(id(node))
                if layout_node is None or layout_node.parent is None:
                    continue
                child_point = pos_map[id(node)]
                parent_x = pos_map[id(layout_node.parent.node)].x()
                marks.append((child_point.y(), parent_x, child_point.x(), length))
            # Broken edges can lie a few pixels apart (the two root edges of a
            # root on the long branch, in a short panel) or on one row (two
            # broken edges on one path). Top to bottom, then left to right,
            # their labels go above and below in turn and their marks sit at
            # 35% and 65% of the edge in turn, so no two marks or labels meet.
            marks.sort(key=lambda mark: (mark[0], mark[1]))
            for index, (y, parent_x, child_x, length) in enumerate(marks):
                text = f"{self._format_break_length(length)} {short_unit}".strip()
                share = 0.5 if len(marks) == 1 else (0.35, 0.65)[index % 2]
                mark = _BranchBreakItem(length, text, label_below=index % 2 == 1)
                mark.setPos(parent_x + share * (child_x - parent_x), y)
                self._add_device_stable_item(mark)
        if use_branch_lengths:
            if self._scale_bar_position == "bottom":
                # Under the tree, below the last row.
                bar_y = last_row_y + leaf_step * 0.5 + 16.0
            else:
                bar_y = top_margin - 14.0
            self._draw_distance_scale_bar(
                left_margin=left_margin,
                bar_y=bar_y,
                branch_scale=branch_scale,
                max_distance=max_distance,
            )

        label_offset = max(10.0, min(22.0, depth_step * 0.35))
        connector_gap = 3.0
        label_right_by_leaf: Dict[str, float] = {}
        # Each row's leaves (one, or a collapsed subtree's) and whether collapsed.
        rows: List[Tuple[List[str], bool]] = []
        # The collapsed subtrees' labels (screen stable), their x and leaves.
        collapsed_labels: List[Tuple[_CollapsedLabelItem, float, List[str]]] = []

        for layout_leaf in layout_leaves:
            leaf = layout_leaf.node
            leaf_point = pos_map[id(leaf)]
            label = label_by_node.get(id(leaf))
            if label is not None:
                members = [name for name in subtree_leaf_map.get(id(leaf), []) if name]
                far_x = far_x_by_node[id(leaf)]
                half_height = half_height_of[id(leaf)]
                label_x = far_x + label_offset * 0.6
                label.setPos(label_x, leaf_point.y())
                collapsed_labels.append((label, label_x, members))
                count = len(members)
                for index, name in enumerate(members):
                    share = 0.5 if count == 1 else index / float(count - 1)
                    y = leaf_point.y() - half_height + share * 2.0 * half_height
                    self._leaf_order.append(name)
                    self._leaf_pos_scene[name] = QtCore.QPointF(far_x, y)
                rows.append((members, True))
                continue
            leaf_name = leaf_name_map[id(leaf)]
            connector_color = self._edge_color_for_subtree([leaf_name])
            self._leaf_order.append(leaf_name)
            self._leaf_pos_scene[leaf_name] = leaf_point
            rows.append(([leaf_name], False))

            label_text = self._display_leaf_name(leaf_name)
            label_item = _LeafLabelItem(label_text, leaf_name, self)
            label_item.setFont(self._leaf_font(leaf_name))
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

        if collapsed_labels:
            # A collapsed subtree's label is in device pixels: the clade bands
            # and cluster brackets beside it leave its width at the x scale
            # it will be shown at (the view's; a fit's: the viewport over the
            # content, estimated twice, the second time with the labels).
            scale_x = abs(float(self.transform().m11())) or 1.0
            if self._auto_fit and self.viewport() is not None:
                left = left_margin - 20.0
                right = max(label_right_by_leaf.values(), default=left_margin) + 80.0
                width = float(max(1, self.viewport().width()))
                for _pass in range(2):
                    scale_x = min(1.0, width / max(1.0, right - left))
                    right = max(
                        [right]
                        + [x + item.boundingRect().width() / scale_x + 80.0 for item, x, _m in collapsed_labels]
                    )
            # Out of the fit, the labels are sized for the camera zoom.
            grown = 1.0 if self._auto_fit else self._stable_scale()
            for item, x, names in collapsed_labels:
                for name in names:
                    label_right_by_leaf[name] = x + (item.boundingRect().width() * grown + 8.0) / max(0.05, scale_x)
        track_anchor_by_leaf = dict(label_right_by_leaf)
        if self._use_group_backgrounds:
            track_anchor_by_leaf = self._draw_group_backgrounds(
                self._annotation_entries(rows, self._leaf_group),
                min((pos_map[id(row.node)].x() for row in layout_leaves), default=0.0),
                label_right_by_leaf,
            )
        content_right_by_leaf = track_anchor_by_leaf

        if self._show_cluster_annotations:
            self._draw_cluster_annotations(self._annotation_entries(rows, self._leaf_cluster), content_right_by_leaf)
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

    def _annotation_entries(
        self,
        rows: Sequence[Tuple[List[str], bool]],
        group_of: Callable[[str], Optional[str]],
    ) -> List[Tuple[str, float]]:
        """The rows as (leaf, y) for the clade bands or cluster brackets: a
        collapsed subtree of one group gives its leaves (along its triangle,
        so the band or bracket spans it); one of several groups gives a
        break ("")."""
        entries: List[Tuple[str, float]] = []
        for names, collapsed in rows:
            if collapsed and len({group_of(name) for name in names}) > 1:
                entries.append(("", self._leaf_pos_scene[names[0]].y()))
                continue
            entries.extend((name, self._leaf_pos_scene[name].y()) for name in names)
        return entries

    @staticmethod
    def _format_branch_scale_value(value: float) -> str:
        """The scale bar value in plain decimals (0.000002, not 2e-06)."""
        return format(Decimal(f"{float(value):.12g}").normalize(), "f")

    @classmethod
    def _format_break_length(cls, value: float) -> str:
        """A broken edge's true length in plain decimals: 3 significant
        digits below 100 (substitutions/site), else to one decimal with
        thousands separators (Gubbins SNP counts: 36,991.5, not 37000)."""
        value = float(value)
        if abs(value) >= 100.0:
            return f"{value:,.1f}".removesuffix(".0")
        return cls._format_branch_scale_value(float(f"{value:.3g}"))

    @staticmethod
    def _nice_branch_scale_value(max_distance: float) -> float:
        target = float(max_distance) * 0.2
        if not math.isfinite(target) or target <= 0.0:
            return 0.0
        magnitude = 10.0 ** math.floor(math.log10(target))
        for multiplier in (5.0, 2.0, 1.0):
            candidate = multiplier * magnitude
            if candidate <= target * (1.0 + 1e-12):
                return candidate
        return magnitude

    def _draw_distance_scale_bar(
        self,
        *,
        left_margin: float,
        bar_y: float,
        branch_scale: float,
        max_distance: float,
    ) -> None:
        scale_value = self._nice_branch_scale_value(max_distance)
        scale_width = scale_value * float(branch_scale)
        if scale_value <= 0.0 or scale_width <= 0.0:
            return

        tick_height = 4.0
        path = QtGui.QPainterPath(QtCore.QPointF(left_margin, bar_y - tick_height))
        path.lineTo(left_margin, bar_y + tick_height)
        path.moveTo(left_margin, bar_y)
        path.lineTo(left_margin + scale_width, bar_y)
        path.moveTo(left_margin + scale_width, bar_y - tick_height)
        path.lineTo(left_margin + scale_width, bar_y + tick_height)

        bar_item = QtWidgets.QGraphicsPathItem(path)
        pen = QtGui.QPen(QtGui.QColor(_TREE_BRANCH_COLOR))
        pen.setWidthF(1.2)
        pen.setCapStyle(QtCore.Qt.PenCapStyle.RoundCap)
        bar_item.setPen(pen)
        bar_item.setData(_SCENE_DATA_KIND, "distance-scale-bar")
        bar_item.setData(_SCENE_DATA_VALUE, scale_value)
        self._scene.addItem(bar_item)

        label_text = (
            f"{self._format_branch_scale_value(scale_value)} "
            f"{self._branch_length_unit or 'branch-length units'}"
        )
        # Name the lengths actually drawn: a refused fit draws the tree's own.
        lengths_name = self._lengths_label if self._lengths_root is not None else self._own_lengths_label
        if lengths_name:
            label_text += f", {lengths_name[:1].lower()}{lengths_name[1:]}"
        if self._branch_breaks_drawn:
            label_text += " (broken branches not to scale)"
        label_item = QtWidgets.QGraphicsSimpleTextItem(label_text)
        font = QtGui.QFont(self.font())
        font.setPointSizeF(8.5)
        font.setWeight(QtGui.QFont.Weight.Medium)
        label_item.setFont(font)
        label_item.setBrush(QtGui.QBrush(QtGui.QColor(_TREE_BRANCH_COLOR)))
        label_bounds = label_item.boundingRect()
        label_item.setData(_SCENE_DATA_KIND, "distance-scale-label")
        label_item.setData(_SCENE_DATA_VALUE, scale_value)
        # Readable at any fit or zoom (a tree stretched to a short panel
        # squashes scene text): the label keeps its size on screen, anchored
        # beside the bar's end and centered on it.
        anchor = QtWidgets.QGraphicsRectItem(QtCore.QRectF())
        anchor.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
        label_item.setParentItem(anchor)
        label_item.setPos(0.0, -label_bounds.height() * 0.5)
        anchor.setPos(left_margin + scale_width + 7.0, bar_y)
        self._add_device_stable_item(anchor)

    def _draw_bottom_legend_sections(
        self,
        left_x: float,
        top_y: float,
        right_x: float,
    ) -> None:
        if not self._bottom_legend_sections:
            return
        items_per_column = 5
        swatch_size = 12.0
        label_gap = 7.0
        column_gap = 24.0
        row_gap = 5.0
        section_gap = 36.0
        x_cursor = float(left_x)
        top = float(top_y)

        title_font = QtGui.QFont(self.font())
        title_font.setWeight(QtGui.QFont.Weight.DemiBold)
        title_font.setPointSizeF(_CLADE_FONT_SIZE)
        text_font = QtGui.QFont(self.font())
        text_font.setPointSizeF(_LEAF_FONT_SIZE)
        title_metrics = QtGui.QFontMetricsF(title_font)
        text_metrics = QtGui.QFontMetricsF(text_font)

        for section in self._bottom_legend_sections:
            title = str(section.get("title") or "").strip()
            items = list(section.get("items", []))
            if not title or not items:
                continue

            title_item = QtWidgets.QGraphicsSimpleTextItem(title)
            title_item.setFont(title_font)
            title_item.setBrush(QtGui.QBrush(QtGui.QColor(_TREE_LABEL_COLOR)))
            title_item.setPos(x_cursor, top)
            self._scene.addItem(title_item)
            section_content_y = top + title_metrics.height() + 7.0

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

                swatch_y = item_y + max(0.0, (text_metrics.height() - swatch_size) * 0.5)
                swatch_rect = QtCore.QRectF(x, swatch_y, swatch_size, swatch_size)
                swatch_item = QtWidgets.QGraphicsRectItem(swatch_rect)
                swatch_pen = QtGui.QPen(QtGui.QColor(COLORS["border_strong"]))
                swatch_pen.setWidthF(0.8)
                swatch_item.setPen(swatch_pen)
                swatch_item.setBrush(QtGui.QBrush(QtGui.QColor(color)))
                self._scene.addItem(swatch_item)

                text_item = QtWidgets.QGraphicsSimpleTextItem(label_text)
                text_item.setFont(text_font)
                text_item.setBrush(QtGui.QBrush(QtGui.QColor(_TREE_LABEL_COLOR)))
                text_item.setPos(x + swatch_size + label_gap, item_y)
                self._scene.addItem(text_item)

            section_width = max(
                title_metrics.horizontalAdvance(title),
                sum(column_widths) + max(0, column_count - 1) * column_gap,
            )
            x_cursor += section_width + section_gap

    @staticmethod
    def _device_stable_local_bounds(
        item: QtWidgets.QGraphicsItem,
    ) -> QtCore.QRectF:
        """The item's extent in device pixels from its anchor (with its
        camera scale, _apply_camera_to_items)."""
        bounds = item.boundingRect().united(item.childrenBoundingRect())
        scale = float(item.scale())
        if abs(scale - 1.0) <= 1e-9:
            return bounds
        return QtCore.QRectF(bounds.x() * scale, bounds.y() * scale, bounds.width() * scale, bounds.height() * scale)

    def _declutter_group_chips(self) -> None:
        """Hide the clade and cluster chips that would overlap a chip of a
        bigger run on screen (a tall tree fitted to a short panel puts
        them closer than their height); the others show again when the
        zoom leaves them room. Exports show all of them."""
        chips = [
            item
            for item in self._device_stable_items
            if item is not None
            and item.scene() is self._scene
            and item.data(_SCENE_DATA_KIND) in ("clade-label", "cluster-label")
        ]
        if len(chips) < 2:
            for item in chips:
                item.setVisible(True)
            return
        transform = self.viewportTransform()
        # Which chips clash depends on the scales only, not on a pan.
        key = (float(transform.m11()), float(transform.m22()), self._stable_scale(), len(chips))
        if key == self._declutter_key:
            return
        self._declutter_key = key
        placed: List[Tuple[int, float, QtWidgets.QGraphicsItem, QtCore.QRectF]] = []
        for item in chips:
            anchor = transform.map(item.scenePos())
            rect = self._device_stable_local_bounds(item).translated(anchor)
            size = item.data(_SCENE_DATA_RUN_SIZE)
            placed.append((-int(size or 0), float(rect.top()), item, rect))
        placed.sort(key=lambda entry: (entry[0], entry[1]))
        kept: List[QtCore.QRectF] = []
        for _size, _top, item, rect in placed:
            # A pixel of air between chips.
            clash = any(rect.adjusted(-1.0, -1.0, 1.0, 1.0).intersects(other) for other in kept)
            if item.isVisible() == clash:
                item.setVisible(not clash)
            if not clash:
                kept.append(rect)

    def _add_device_stable_item(
        self,
        item: QtWidgets.QGraphicsItem,
    ) -> None:
        item.setFlag(
            QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations,
            True,
        )
        self._scene.addItem(item)
        self._device_stable_items.append(item)

    def _draw_cluster_annotations(
        self,
        entries: Sequence[Tuple[str, float]],
        label_right_by_leaf: Mapping[str, float],
    ) -> None:
        """The cluster brackets and chips right of the rows (``entries``:
        each drawn leaf and its y, top to bottom; "" breaks a run)."""
        if not self._leaf_cluster_map or not entries:
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

        for leaf_name, leaf_y in entries:
            cluster_name = self._leaf_cluster(leaf_name) if leaf_name else None
            if not cluster_name:
                flush_run()
                continue
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

        cluster_run_counts: Dict[str, int] = {}
        for cluster_name, _cluster_leaves, _top, _bottom in cluster_runs:
            cluster_run_counts[cluster_name] = cluster_run_counts.get(cluster_name, 0) + 1
        cluster_run_indices: Dict[str, int] = {}

        font = QtGui.QFont(self.font())
        font.setPointSizeF(_CLUSTER_FONT_SIZE)
        font.setWeight(QtGui.QFont.Weight.DemiBold)
        bracket_gap = 8.0
        bracket_width = 8.0
        label_gap = 7.0
        chip_pad_x = 6.0
        chip_pad_y = 3.0

        for cluster_name, cluster_leaves, raw_top, raw_bottom in cluster_runs:
            linked_cluster_leaves = list(
                dict.fromkeys(cluster_leaves_by_name.get(cluster_name, cluster_leaves))
            )
            run_index = cluster_run_indices.get(cluster_name, 0) + 1
            cluster_run_indices[cluster_name] = run_index
            run_count = cluster_run_counts.get(cluster_name, 1)
            top_y = raw_top
            bottom_y = raw_bottom
            if math.isclose(top_y, bottom_y):
                top_y -= 6.0
                bottom_y += 6.0

            color = self._visual_color_for_subtree(
                cluster_leaves,
                QtGui.QColor(_TREE_CLUSTER_COLOR),
            )
            pen = self._edge_pen_for_color(color, width=1.5)

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
            bracket_item.setData(_SCENE_DATA_KIND, "cluster-bracket")
            bracket_item.setData(_SCENE_DATA_NAME, cluster_name)
            bracket_item.setData(_SCENE_DATA_RUN_INDEX, run_index)
            bracket_item.setData(_SCENE_DATA_RUN_COUNT, run_count)
            self._scene.addItem(bracket_item)

            short_label = self._compact_cluster_label(cluster_name)
            genome_count = len(linked_cluster_leaves)
            genome_word = "genome" if genome_count == 1 else "genomes"
            tooltip = f"Cluster {cluster_name}"
            if run_count > 1:
                tooltip += f" · block {run_index} of {run_count}"
            tooltip += f" · {genome_count} {genome_word}"
            bracket_item.setToolTip(tooltip)

            label_item = _ActionTextItem(
                short_label,
                callback=lambda payload=cluster_payload: self._on_cluster_click(payload),
                menu_callback=lambda pos, payload=cluster_payload: self._show_group_menu("cluster", payload, pos),
                enable_hover=self._enable_hover_interactions,
                hover_enter_callback=lambda leaves=tuple(linked_cluster_leaves): self._on_action_hover_enter(leaves),
                hover_leave_callback=lambda leaves=tuple(linked_cluster_leaves): self._on_action_hover_leave(leaves),
            )
            label_item.setFont(font)
            label_item.setBrush(QtGui.QBrush(color))
            label_item.setToolTip(tooltip)
            label_item.setData(_SCENE_DATA_KIND, "cluster-label")
            label_item.setData(_SCENE_DATA_NAME, cluster_name)
            label_item.setData(_SCENE_DATA_RUN_INDEX, run_index)
            label_item.setData(_SCENE_DATA_RUN_COUNT, run_count)
            label_bounds = label_item.boundingRect()

            chip_rect = label_bounds.adjusted(
                -chip_pad_x,
                -chip_pad_y,
                chip_pad_x,
                chip_pad_y,
            )
            chip_path = QtGui.QPainterPath()
            chip_path.addRoundedRect(chip_rect, 5.0, 5.0)
            chip_item = QtWidgets.QGraphicsPathItem(chip_path)
            chip_pen_color = QtGui.QColor(color)
            chip_pen_color.setAlpha(120)
            chip_pen = QtGui.QPen(chip_pen_color)
            chip_pen.setWidthF(0.9)
            chip_item.setPen(chip_pen)
            chip_fill = QtGui.QColor(_TREE_CANVAS_COLOR)
            chip_fill.setAlpha(244)
            chip_item.setBrush(QtGui.QBrush(chip_fill))
            chip_item.setToolTip(tooltip)
            chip_item.setData(_SCENE_DATA_KIND, "cluster-chip")
            chip_item.setData(_SCENE_DATA_NAME, cluster_name)
            chip_item.setData(_SCENE_DATA_RUN_INDEX, run_index)
            chip_item.setData(_SCENE_DATA_RUN_COUNT, run_count)
            chip_item.setAcceptedMouseButtons(QtCore.Qt.MouseButton.NoButton)
            chip_item.setParentItem(label_item)
            chip_item.setFlag(
                QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemStacksBehindParent,
                True,
            )
            chip_item.setZValue(-1.0)
            chip_left = outer_x + label_gap
            chip_top = ((top_y + bottom_y) * 0.5) - (chip_rect.height() * 0.5)
            label_item.setPos(
                chip_left - chip_rect.left(),
                chip_top - chip_rect.top(),
            )
            label_item.setData(_SCENE_DATA_RUN_SIZE, len(cluster_leaves))
            self._add_device_stable_item(label_item)

    def _draw_group_backgrounds(
        self,
        entries: Sequence[Tuple[str, float]],
        fallback_left: float,
        label_right_by_leaf: Mapping[str, float],
    ) -> Dict[str, float]:
        """The clade bands right of the rows (``entries``: each drawn leaf
        and its y, top to bottom; "" breaks a run)."""
        content_right_by_leaf = dict(label_right_by_leaf)
        if not self._leaf_group_map or not entries:
            return content_right_by_leaf

        group_runs: List[Tuple[str, List[str], float, float]] = []
        group_leaves_by_name: Dict[str, List[str]] = {}
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

        for leaf_name, leaf_y in entries:
            group_name = self._leaf_group(leaf_name) if leaf_name else None
            if not group_name:
                flush_run()
                continue
            if group_name != current_group:
                flush_run()
                current_group = group_name
                current_leaves = [leaf_name]
                group_leaves_by_name.setdefault(group_name, []).append(leaf_name)
                start_y = leaf_y
                end_y = leaf_y
                continue
            current_leaves.append(leaf_name)
            group_leaves_by_name.setdefault(group_name, []).append(leaf_name)
            end_y = leaf_y
        flush_run()

        if not group_runs:
            return content_right_by_leaf

        run_counts: Dict[str, int] = {}
        for group_name, _group_leaves, _top, _bottom in group_runs:
            run_counts[group_name] = run_counts.get(group_name, 0) + 1
        run_indices: Dict[str, int] = {}

        row_pad = 8.0
        min_height = 18.0
        lane_pad = 4.0
        rail_gap = 8.0
        chip_gap = 8.0
        chip_pad_x = 7.0
        chip_pad_y = 3.0

        text_font = QtGui.QFont(self.font())
        text_font.setPointSizeF(_CLADE_FONT_SIZE)
        text_font.setWeight(QtGui.QFont.Weight.DemiBold)
        text_metrics = QtGui.QFontMetricsF(text_font)
        for group_name, group_leaves, raw_top, raw_bottom in group_runs:
            linked_group_leaves = list(
                dict.fromkeys(group_leaves_by_name.get(group_name, group_leaves))
            )
            run_index = run_indices.get(group_name, 0) + 1
            run_indices[group_name] = run_index
            run_count = run_counts.get(group_name, 1)
            top_y = raw_top - row_pad
            bottom_y = raw_bottom + row_pad
            if (bottom_y - top_y) < min_height:
                center_y = (top_y + bottom_y) * 0.5
                top_y = center_y - (min_height * 0.5)
                bottom_y = center_y + (min_height * 0.5)

            rect_height = max(min_height, bottom_y - top_y)
            genome_count = len(linked_group_leaves)
            genome_word = "genome" if genome_count == 1 else "genomes"
            label_text = f"Clade {group_name} · {genome_count} {genome_word}"
            text_width = float(text_metrics.horizontalAdvance(label_text))
            text_height = float(text_metrics.height())

            base_right_x = max(
                (
                    label_right_by_leaf.get(leaf_name, fallback_left + 40.0)
                    for leaf_name in group_leaves
                ),
                default=fallback_left + 40.0,
            )
            rail_x = base_right_x + rail_gap
            if self._show_group_labels:
                chip_left = rail_x + _CLADE_RAIL_WIDTH + chip_gap - chip_pad_x
                if rect_height >= text_height + 16.0:
                    text_top = top_y + 6.0
                else:
                    text_top = ((top_y + bottom_y) * 0.5) - (text_height * 0.5)
                chip_top = text_top - chip_pad_y
                chip_rect = QtCore.QRectF(
                    chip_left,
                    chip_top,
                    text_width + (chip_pad_x * 2.0),
                    text_height + (chip_pad_y * 2.0),
                )
                right_x = chip_rect.right() + lane_pad
            else:
                chip_left = rail_x
                chip_top = top_y
                chip_rect = QtCore.QRectF()
                right_x = rail_x + _CLADE_RAIL_WIDTH + lane_pad

            band_rect = QtCore.QRectF(
                rail_x,
                top_y,
                max(24.0, right_x - rail_x),
                rect_height,
            )
            path = QtGui.QPainterPath()
            path.addRoundedRect(band_rect, 6.0, 6.0)
            bg_item = QtWidgets.QGraphicsPathItem(path)
            bg_item.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
            bg_item.setBrush(self._group_background_brush(group_name, group_leaves))
            bg_item.setZValue(-10.0)
            bg_item.setData(_SCENE_DATA_KIND, "clade-band")
            bg_item.setData(_SCENE_DATA_NAME, group_name)
            bg_item.setData(_SCENE_DATA_RUN_INDEX, run_index)
            bg_item.setData(_SCENE_DATA_RUN_COUNT, run_count)
            self._scene.addItem(bg_item)

            accent_color = self._visual_color_for_subtree(
                group_leaves,
                self._group_color(group_name),
            )
            rail_item = QtWidgets.QGraphicsRectItem(
                QtCore.QRectF(rail_x, top_y, _CLADE_RAIL_WIDTH, rect_height)
            )
            rail_item.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
            rail_item.setBrush(QtGui.QBrush(accent_color))
            rail_item.setZValue(-9.0)
            rail_item.setData(_SCENE_DATA_KIND, "clade-rail")
            rail_item.setData(_SCENE_DATA_NAME, group_name)
            rail_item.setData(_SCENE_DATA_RUN_INDEX, run_index)
            rail_item.setData(_SCENE_DATA_RUN_COUNT, run_count)
            self._scene.addItem(rail_item)

            if self._show_group_labels:
                chip_path = QtGui.QPainterPath()
                local_chip_rect = QtCore.QRectF(
                    -chip_pad_x,
                    -chip_pad_y,
                    chip_rect.width(),
                    chip_rect.height(),
                )
                chip_path.addRoundedRect(local_chip_rect, 6.0, 6.0)
                chip_item = QtWidgets.QGraphicsPathItem(chip_path)
                chip_pen_color = QtGui.QColor(accent_color)
                chip_pen_color.setAlpha(110)
                chip_pen = QtGui.QPen(chip_pen_color)
                chip_pen.setWidthF(0.8)
                chip_item.setPen(chip_pen)
                chip_fill = QtGui.QColor(_TREE_CANVAS_COLOR)
                chip_fill.setAlpha(235)
                chip_item.setBrush(QtGui.QBrush(chip_fill))
                chip_item.setZValue(-1.0)
                chip_item.setData(_SCENE_DATA_KIND, "clade-chip")
                chip_item.setData(_SCENE_DATA_NAME, group_name)
                chip_item.setData(_SCENE_DATA_RUN_INDEX, run_index)
                chip_item.setData(_SCENE_DATA_RUN_COUNT, run_count)

                clade_payload = {"name": group_name, "samples": list(linked_group_leaves)}
                text_item = _ActionTextItem(
                    label_text,
                    callback=lambda payload=clade_payload: self._on_clade_click(payload),
                    menu_callback=lambda pos, payload=clade_payload: self._show_group_menu("clade", payload, pos),
                    enable_hover=self._enable_hover_interactions,
                    hover_enter_callback=lambda leaves=tuple(linked_group_leaves): self._on_action_hover_enter(leaves),
                    hover_leave_callback=lambda leaves=tuple(linked_group_leaves): self._on_action_hover_leave(leaves),
                )
                text_item.setFont(text_font)
                text_item.setBrush(QtGui.QBrush(QtGui.QColor(_TREE_LABEL_COLOR)))
                tooltip = label_text
                if run_count > 1:
                    tooltip += f" · block {run_index} of {run_count}"
                text_item.setToolTip(tooltip)
                text_item.setData(_SCENE_DATA_KIND, "clade-label")
                text_item.setData(_SCENE_DATA_NAME, group_name)
                text_item.setData(_SCENE_DATA_RUN_INDEX, run_index)
                text_item.setData(_SCENE_DATA_RUN_COUNT, run_count)
                chip_item.setAcceptedMouseButtons(QtCore.Qt.MouseButton.NoButton)
                chip_item.setParentItem(text_item)
                chip_item.setFlag(
                    QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemStacksBehindParent,
                    True,
                )
                chip_item.setZValue(-1.0)
                text_item.setPos(
                    chip_left + chip_pad_x,
                    chip_top + chip_pad_y,
                )
                text_item.setData(_SCENE_DATA_RUN_SIZE, len(group_leaves))
                self._add_device_stable_item(text_item)

            for leaf_name in group_leaves:
                content_right_by_leaf[leaf_name] = max(
                    content_right_by_leaf.get(leaf_name, 0.0),
                    right_x,
                )

        return content_right_by_leaf

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
            return QtGui.QColor(_TREE_HOVER_INK)
        if self._active_samples is not None and leaf_name not in self._active_samples:
            return out_of_set_color()
        if self._reference_sample and leaf_name == self._reference_sample:
            return QtGui.QColor(_TREE_REFERENCE_COLOR)
        return base

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
                self._set_label_hovered(old_label, False)

        active_items: Set[_TreeEdgeItem] = set()
        active_markers: Set[_NodeMarkerItem] = set()
        for leaf in ordered:
            for item in self._edge_items_by_leaf.get(leaf, []):
                active_items.add(item)
            for marker in self._marker_items_by_leaf.get(leaf, []):
                active_markers.add(marker)
            label_item = self._leaf_label_items.get(leaf)
            if label_item is not None:
                label_item.setBrush(self._label_brush(QtGui.QColor(_TREE_HOVER_INK)))
                self._set_label_hovered(label_item, True)

        for item in active_items:
            item.set_hovered(True)
        for item in active_markers:
            item.set_hovered(True)

        self._active_edge_items = active_items
        self._active_marker_items = active_markers
        self._active_leaf_set = next_set
        self._hover_leaf = ordered[0] if len(ordered) == 1 else None
        self._set_hover_dimmed(bool(next_set))
        # The dot and chips are drawn over the viewport (drawForeground).
        viewport = self.viewport()
        if viewport is not None:
            viewport.update()

        if emit_signals:
            if len(ordered) == 1:
                self.leafHovered.emit(ordered[0])
            else:
                self.leafHovered.emit(None)
            if len(ordered) >= 2:
                self.subtreeHovered.emit(ordered)
            elif not ordered:
                self.subtreeHovered.emit(None)

    # ---- hover emphasis: thick amber path, tip dot and chip, dim the rest ----

    def _adopt_tree_items(self) -> None:
        """Put every drawn item under one parent so a hover dims the tree
        with a single opacity change (no per-item walk per hover)."""
        content = _TreeContentItem()
        self._scene.addItem(content)
        # Ascending stacking order keeps the drawing order of equal z items.
        for item in self._scene.items(QtCore.Qt.SortOrder.AscendingOrder):
            if item is content or item.parentItem() is not None:
                continue
            item.setParentItem(content)
        self._content_item = content

    @staticmethod
    def _set_label_hovered(label: QtWidgets.QGraphicsItem, active: bool) -> None:
        try:
            label.setZValue(_HOVER_Z if active else 0.0)
            label.setFlag(_IGNORES_PARENT_OPACITY, active)
        except RuntimeError:
            return

    def _set_hover_dimmed(self, dimmed: bool) -> None:
        if dimmed == self._hover_dimmed:
            return
        self._hover_dimmed = dimmed
        content = self._content_item
        if content is not None:
            try:
                content.setOpacity(_HOVER_DIM_OPACITY if dimmed else 1.0)
            except RuntimeError:
                return

    def hover_dimmed(self) -> bool:
        """Whether the tree is dimmed for an active hover."""
        return self._hover_dimmed

    def _hover_chip_texts(self, leaf_name: str) -> Tuple[str, str]:
        """The chip's sample name and its clade and cluster ("Clade 2 ·
        Cluster 1"), each when known."""
        details: List[str] = []
        clade = self._leaf_group_map.get(leaf_name)
        if clade:
            details.append(clade if clade.lower().startswith("clade") else f"Clade {clade}")
        # Only the cluster map: display groups color the branches and can be
        # clades, so they would be misread as a cluster here.
        cluster = self._leaf_cluster_map.get(leaf_name)
        if cluster:
            details.append(cluster if cluster.lower().startswith("cluster") else f"Cluster {cluster}")
        return self._display_leaf_name(leaf_name), " \u00b7 ".join(details)

    def _hover_chip_fonts(self) -> Tuple[QtGui.QFont, QtGui.QFont]:
        detail_font = QtGui.QFont(self.font())
        detail_font.setPointSizeF(9.0)
        name_font = QtGui.QFont(detail_font)
        name_font.setWeight(QtGui.QFont.Weight.DemiBold)
        return name_font, detail_font

    def hover_overlay(self) -> HoverOverlay:
        """The dots and chips of the active hover in viewport pixels (empty
        without a hover). They are sized in pixels, so they do not scale with
        the zoom, and every chip is kept inside the viewport: it sits right of
        its tip, flips to the left near the right edge and is clamped."""
        overlay = HoverOverlay()
        viewport = self.viewport()
        if not self._active_leaf_set or viewport is None:
            return overlay
        view_rect = QtCore.QRectF(viewport.rect())
        transform = self.viewportTransform()
        radius = _HOVER_DOT_DIAMETER_PX / 2.0
        visible: List[Tuple[str, QtCore.QPointF]] = []
        offscreen: Dict[str, int] = {}
        # A collapsed subtree's hovered leaves get one dot and chip (its
        # label, and how many are hovered when not all).
        group_texts: Dict[str, Tuple[str, str]] = {}
        seen_groups: Set[frozenset] = set()
        for leaf in self._leaf_order:
            if leaf not in self._active_leaf_set:
                continue
            count = 1
            group = self._collapsed_of_leaf.get(leaf)
            if group is not None:
                if group in seen_groups:
                    continue
                seen_groups.add(group)
                members, label, scene_pt = self._collapsed_drawn[group]
                count = sum(1 for name in members if name in self._active_leaf_set)
                group_texts[label] = (label, "" if count == len(members) else f"{count} hovered")
                leaf = label
            else:
                scene_pt = self._leaf_pos_scene.get(leaf)
            if scene_pt is None:
                continue
            pt = transform.map(scene_pt)
            if not view_rect.contains(pt):
                if pt.y() < view_rect.top():
                    side = "above"
                elif pt.y() > view_rect.bottom():
                    side = "below"
                else:
                    side = "right" if pt.x() > view_rect.right() else "left"
                offscreen[side] = offscreen.get(side, 0) + count
                continue
            visible.append((leaf, pt))
            overlay.dots.append(
                (leaf, QtCore.QRectF(pt.x() - radius, pt.y() - radius, 2.0 * radius, 2.0 * radius))
            )

        name_font, detail_font = self._hover_chip_fonts()
        name_fm = QtGui.QFontMetricsF(name_font)
        detail_fm = QtGui.QFontMetricsF(detail_font)
        pad_x, pad_y = 7.0, 3.0
        chip_h = max(name_fm.height(), detail_fm.height()) + 2.0 * pad_y
        margin = _HOVER_CHIP_MARGIN_PX
        max_w = max(40.0, view_rect.width() - 2.0 * margin)
        stack_gap = 2.0

        # Off-screen markers first: the name chips stack between them.
        top_limit = view_rect.top() + margin
        bottom_limit = view_rect.bottom() - margin
        if offscreen:
            overlay.edges = self._hover_edge_markers(offscreen, view_rect, chip_h, detail_fm, pad_x)
            for side, _text, rect in overlay.edges:
                if side == "above":
                    top_limit = max(top_limit, rect.bottom() + stack_gap)
                elif side == "below":
                    bottom_limit = min(bottom_limit, rect.top() - stack_gap)
        if not visible:
            return overlay

        def place(width: float, pt: QtCore.QPointF) -> QtCore.QRectF:
            width = min(width, max_w)
            x = pt.x() + _HOVER_CHIP_GAP_PX
            if x + width > view_rect.right() - margin:
                x = pt.x() - _HOVER_CHIP_GAP_PX - width  # flip to the left
            x = max(view_rect.left() + margin, min(x, view_rect.right() - margin - width))
            y = pt.y() - chip_h / 2.0
            return QtCore.QRectF(x, y, width, chip_h)

        # Only as many chips as the viewport has rows for: the rest fold into
        # "+N more", and a viewport too short for even one chip keeps the dots.
        # The off-screen markers count the tips out of view, so "+N more"
        # counts only the tips on screen without a chip.
        rows = int((bottom_limit - top_limit + stack_gap) // (chip_h + stack_gap))
        if rows <= 0:
            return overlay
        total = len(visible)
        shown_count = min(_HOVER_CHIP_CAP, total)
        if shown_count < total or shown_count > rows:
            shown_count = min(shown_count, rows - 1)  # one row for "+N more"

        chips: List[Tuple[str, str, QtCore.QRectF]] = []
        tips: List[QtCore.QPointF] = []
        for leaf, pt in visible[:shown_count]:
            name, detail = group_texts.get(leaf) or self._hover_chip_texts(leaf)
            width = name_fm.horizontalAdvance(name) + 2.0 * pad_x
            if detail:
                width += 8.0 + detail_fm.horizontalAdvance(detail)
            chips.append((name, detail, place(width, pt)))
            tips.append(pt)
        rest = total - len(chips)
        more: Optional[Tuple[str, QtCore.QRectF]] = None
        if rest > 0:
            text = f"+{rest} more" if chips else f"{rest} samples"
            width = min(detail_fm.horizontalAdvance(text) + 2.0 * pad_x, max_w)
            if chips:
                last = chips[-1][2]
                more = (text, QtCore.QRectF(last.left(), last.bottom() + stack_gap, width, chip_h))
            else:
                more = (text, place(width, visible[0][1]))

        # Stack the chips (and "+N more") without overlaps, then keep the
        # stack inside the viewport.
        rects = [rect for _n, _d, rect in chips] + ([more[1]] if more else [])
        for prev, rect in zip(rects, rects[1:]):
            if rect.top() < prev.bottom() + stack_gap:
                rect.moveTop(prev.bottom() + stack_gap)
        shift = 0.0
        if rects[-1].bottom() > bottom_limit:
            shift = bottom_limit - rects[-1].bottom()
        if rects[0].top() + shift < top_limit:
            shift = top_limit - rects[0].top()
        if shift:
            for rect in rects:
                rect.translate(0.0, shift)
        # A chip pushed away from its tip by the stack is tied back to its
        # dot with a thin leader line, so it is clear which chip is whose.
        for (_n, _d, rect), pt in zip(chips, tips):
            if abs(rect.center().y() - pt.y()) <= 2.0:
                continue
            edge_x = rect.left() if rect.center().x() >= pt.x() else rect.right()
            overlay.leaders.append((QtCore.QPointF(pt), QtCore.QPointF(edge_x, rect.center().y())))
        overlay.chips = chips
        overlay.more = more
        return overlay

    @staticmethod
    def _hover_edge_markers(
        offscreen: Mapping[str, int],
        view_rect: QtCore.QRectF,
        chip_h: float,
        detail_fm: QtGui.QFontMetricsF,
        pad_x: float,
    ) -> List[Tuple[str, str, QtCore.QRectF]]:
        """One marker per side with hovered tips out of view ("2 hovered
        samples below"), centered on that edge (the side markers sit at the
        middle of their edge). An arrow is drawn in the marker, so the text
        leaves room for it. A viewport too short for a top and a bottom
        marker gets one marker with the total."""
        margin = _HOVER_CHIP_MARGIN_PX
        arrow_w = 12.0
        max_w = max(40.0, view_rect.width() - 2.0 * margin)
        sides = [side for side in ("above", "below", "left", "right") if offscreen.get(side)]
        if view_rect.height() < 3.0 * chip_h + 2.0 * margin and len(sides) > 1:
            count = sum(offscreen.values())
            noun = "sample" if count == 1 else "samples"
            sides = ["off"]
            texts = {"off": f"{count} hovered {noun} off screen"}
        else:
            texts = {}
            for side in sides:
                count = offscreen[side]
                noun = "sample" if count == 1 else "samples"
                where = side if side in ("above", "below") else f"to the {side}"
                texts[side] = f"{count} hovered {noun} {where}"
        markers: List[Tuple[str, str, QtCore.QRectF]] = []
        for side in sides:
            text = texts[side]
            width = min(detail_fm.horizontalAdvance(text) + 2.0 * pad_x + arrow_w, max_w)
            if side == "left":
                x = view_rect.left() + margin
            elif side == "right":
                x = view_rect.right() - margin - width
            else:
                x = view_rect.center().x() - width / 2.0
            if side == "below":
                y = view_rect.bottom() - margin - chip_h
            elif side in ("left", "right"):
                y = view_rect.center().y() - chip_h / 2.0
            else:
                y = view_rect.top() + margin
            markers.append((side, text, QtCore.QRectF(x, y, width, chip_h)))
        return markers

    def drawForeground(self, painter: QtGui.QPainter, rect: QtCore.QRectF) -> None:
        super().drawForeground(painter, rect)
        if not self._active_leaf_set:
            return
        overlay = self.hover_overlay()
        if not overlay.dots and not overlay.chips and not overlay.edges:
            return
        painter.save()
        try:
            # Viewport pixels: the dot and chips never scale with the zoom.
            painter.resetTransform()
            painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
            hover = QtGui.QColor(_TREE_HOVER_COLOR)
            ring = QtGui.QPen(QtGui.QColor(COLORS["surface"]))
            ring.setWidthF(1.5)
            if overlay.leaders:
                leader = QtGui.QPen(hover)
                leader.setWidthF(1.0)
                painter.setPen(leader)
                for start, end in overlay.leaders:
                    painter.drawLine(start, end)
            painter.setPen(ring)
            painter.setBrush(hover)
            for _leaf, dot in overlay.dots:
                painter.drawEllipse(dot)
            name_font, detail_font = self._hover_chip_fonts()
            border = QtGui.QPen(hover)
            border.setWidthF(1.5)
            fill = QtGui.QColor(COLORS["surface"])
            text_pen = QtGui.QPen(QtGui.QColor(COLORS["text_strong"]))
            muted_pen = QtGui.QPen(QtGui.QColor(COLORS["text_muted"]))
            flags = int(QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignVCenter)
            for name, detail, chip in overlay.chips:
                painter.setPen(border)
                painter.setBrush(fill)
                painter.drawRoundedRect(chip, 4.0, 4.0)
                inner = chip.adjusted(7.0, 0.0, -7.0, 0.0)
                # A chip clamped to a narrow viewport elides its text.
                name_fm = QtGui.QFontMetricsF(name_font)
                shown = name_fm.elidedText(name, QtCore.Qt.TextElideMode.ElideRight, inner.width())
                painter.setFont(name_font)
                painter.setPen(text_pen)
                painter.drawText(inner, flags, shown)
                room = inner.width() - name_fm.horizontalAdvance(shown) - 8.0
                if detail and room > 12.0:
                    detail_fm = QtGui.QFontMetricsF(detail_font)
                    painter.setFont(detail_font)
                    painter.setPen(muted_pen)
                    painter.drawText(
                        inner.adjusted(inner.width() - room, 0.0, 0.0, 0.0),
                        flags,
                        detail_fm.elidedText(detail, QtCore.Qt.TextElideMode.ElideRight, room),
                    )
            if overlay.more is not None:
                text, chip = overlay.more
                painter.setPen(border)
                painter.setBrush(fill)
                painter.drawRoundedRect(chip, 4.0, 4.0)
                painter.setFont(detail_font)
                painter.setPen(muted_pen)
                painter.drawText(chip.adjusted(7.0, 0.0, -7.0, 0.0), flags, text)
            for side, text, chip in overlay.edges:
                painter.setPen(border)
                painter.setBrush(fill)
                painter.drawRoundedRect(chip, 4.0, 4.0)
                # A small amber arrow toward the hidden tips, then the count.
                cx = chip.left() + 7.0 + 4.5
                cy = chip.center().y()
                if side == "above":
                    points = [QtCore.QPointF(cx - 4.5, cy + 2.5), QtCore.QPointF(cx + 4.5, cy + 2.5), QtCore.QPointF(cx, cy - 3.0)]
                elif side == "below":
                    points = [QtCore.QPointF(cx - 4.5, cy - 2.5), QtCore.QPointF(cx + 4.5, cy - 2.5), QtCore.QPointF(cx, cy + 3.0)]
                elif side == "left":
                    points = [QtCore.QPointF(cx + 2.5, cy - 4.5), QtCore.QPointF(cx + 2.5, cy + 4.5), QtCore.QPointF(cx - 3.0, cy)]
                elif side == "right":
                    points = [QtCore.QPointF(cx - 2.5, cy - 4.5), QtCore.QPointF(cx - 2.5, cy + 4.5), QtCore.QPointF(cx + 3.0, cy)]
                else:
                    points = []
                if points:
                    painter.setPen(QtCore.Qt.PenStyle.NoPen)
                    painter.setBrush(hover)
                    painter.drawPolygon(QtGui.QPolygonF(points))
                detail_fm = QtGui.QFontMetricsF(detail_font)
                inner = chip.adjusted(7.0 + 12.0, 0.0, -7.0, 0.0)
                painter.setFont(detail_font)
                painter.setPen(text_pen)
                painter.drawText(
                    inner, flags, detail_fm.elidedText(text, QtCore.Qt.TextElideMode.ElideRight, inner.width())
                )
        finally:
            painter.restore()

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

    def group_menu(self, kind: str, payload: Mapping[str, object]) -> QtWidgets.QMenu:
        """The menu of a clade or cluster label: show only its samples in this
        window (what a left click does) or open them in a new explore window."""
        menu = QtWidgets.QMenu(self)
        noun = "clade" if kind == "clade" else "cluster"
        only = menu.addAction(f"Show only this {noun}")
        only.setObjectName("TreeGroupFilterAction")
        new_window = menu.addAction(f"Open {noun} in new explore window")
        new_window.setObjectName("TreeGroupWindowAction")
        new_window.setToolTip(
            f"A new explore window that shows only the samples of this {noun}; this window stays as it is"
        )
        data = dict(payload)
        if kind == "clade":
            only.triggered.connect(lambda _c=False: self._on_clade_click(data))
            new_window.triggered.connect(lambda _c=False: self.cladeWindowRequested.emit(dict(data)))
        else:
            only.triggered.connect(lambda _c=False: self._on_cluster_click(data))
            new_window.triggered.connect(lambda _c=False: self.clusterWindowRequested.emit(dict(data)))
        return menu

    def _show_group_menu(self, kind: str, payload: Mapping[str, object], screen_pos: QtCore.QPoint) -> None:
        menu = self.group_menu(kind, payload)
        menu.exec(screen_pos)
        menu.deleteLater()

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

    def _update_scene_rect_for_device_stable_items(
        self,
        scale: Optional[float] = None,
        scale_y: Optional[float] = None,
    ) -> QtCore.QRectF:
        """Grow the scene rect over the screen stable items at the given
        scales (x, and y when it differs; default: the view's)."""
        bounds = self._scene.itemsBoundingRect()
        if bounds.isNull() or bounds.isEmpty():
            bounds = QtCore.QRectF(self.sceneRect())
        if bounds.isNull() or bounds.isEmpty():
            return bounds

        if scale is None:
            scale_x_raw = float(self.transform().m11())
            scale_y_raw = float(self.transform().m22())
        else:
            scale_x_raw = float(scale)
            scale_y_raw = float(scale_y if scale_y is not None else scale)
        effective_scale_x = max(0.03, abs(scale_x_raw))
        effective_scale_y = max(0.03, abs(scale_y_raw))
        effective_bounds = QtCore.QRectF(bounds)
        for item in self._device_stable_items:
            if item is None or item.scene() is not self._scene:
                continue
            local_bounds = self._device_stable_local_bounds(item)
            anchor = item.scenePos()
            stable_bounds = QtCore.QRectF(
                anchor.x() + (local_bounds.left() / effective_scale_x),
                anchor.y() + (local_bounds.top() / effective_scale_y),
                local_bounds.width() / effective_scale_x,
                local_bounds.height() / effective_scale_y,
            )
            effective_bounds = effective_bounds.united(stable_bounds)

        pad = _SCENE_RECT_PAD
        effective_bounds = effective_bounds.adjusted(-pad, -pad, pad, pad)
        self._scene.setSceneRect(effective_bounds)
        return effective_bounds

    def _fit_content_rect(self) -> QtCore.QRectF:
        """The scene rect of the tree items that scale with the zoom (the
        screen stable labels and chips left out: their scene size depends on
        the zoom)."""
        stable = {id(item) for item in self._device_stable_items if item is not None}
        rect = QtCore.QRectF()
        for item in self._scene.items():
            if item.parentItem() is not None or id(item) in stable:
                continue
            local = item.boundingRect().united(item.childrenBoundingRect())
            bounds = item.mapRectToScene(local)
            if bounds.isNull() or bounds.isEmpty():
                continue
            rect = bounds if rect.isNull() else rect.united(bounds)
        return rect

    def _stable_vertical_extent(self, top: float, scale: float) -> Tuple[float, float]:
        """Viewport pixels above ``top`` and below it that the screen stable
        items take at ``scale`` (the extent of everything, from the scene y
        ``top``)."""
        above = 0.0
        below = 0.0
        for item in self._device_stable_items:
            if item is None or item.scene() is not self._scene:
                continue
            local = self._device_stable_local_bounds(item)
            anchor = (float(item.scenePos().y()) - top) * scale
            above = max(above, -(anchor + float(local.top())))
            below = max(below, anchor + float(local.bottom()))
        return above, below

    def _fit_rect(self) -> Optional[QtCore.QRectF]:
        """The scene rect the fit puts in the viewport (None: nothing to
        fit)."""
        if self._layout_mode == "rectangular":
            rect = self._fit_content_rect()
        else:
            rect = self._scene.itemsBoundingRect()
        rect = rect.adjusted(-20.0, -20.0, 20.0, 20.0)
        if rect.isNull() or rect.isEmpty():
            rect = self.sceneRect()
        if rect.isNull() or rect.isEmpty():
            return None
        return rect

    def _rectangular_fit_scales(self, rect: QtCore.QRectF) -> Optional[Tuple[float, float]]:
        """The (x, y) scales of the rectangular fit: y puts every row (and the
        screen stable chips around the tree) in the viewport height; x fits
        the tree's width, and is never below y. A tree that is tall for its
        panel keeps the width fit across (x above y) instead of shrinking to
        a sliver as wide as a uniform scale would leave it."""
        viewport = self.viewport()
        if viewport is None or viewport.width() <= 0 or viewport.height() <= 0:
            return None
        available_width = max(1.0, float(viewport.width()) - 16.0)
        available_height = max(1.0, float(viewport.height()) - 16.0)
        scale = available_width / max(1.0, float(rect.width()))
        for item in self._device_stable_items:
            if item is None or item.scene() is not self._scene:
                continue
            local_bounds = self._device_stable_local_bounds(item)
            anchor_offset = item.scenePos().x() - rect.left()
            fixed_right = max(0.0, float(local_bounds.right()))
            scalable_span = max(0.0, float(anchor_offset))
            if scalable_span <= 0.0 or fixed_right >= available_width:
                continue
            scale = min(
                scale,
                (available_width - fixed_right) / scalable_span,
            )
        # Every row in the height: the scaled tree plus the screen stable
        # chips around it fit the viewport (a few passes, as the chips'
        # share depends on the scale). In a viewport too short for the
        # chips, the tree alone fits (chips may stick out, the rows do not).
        tree_only = available_height / max(1.0, float(rect.height()))
        width_scale = scale
        for _ in range(4):
            above, below = self._stable_vertical_extent(float(rect.top()), scale)
            extent = above + max(float(rect.height()) * scale, below)
            if extent <= available_height + 0.5:
                break
            if float(rect.height()) * scale >= below:
                # The tree sets the extent: scale it into what the chips
                # above leave.
                scale = max(1e-3, (available_height - above) / max(1.0, float(rect.height())))
            else:
                scale *= max(0.05, available_height / extent)
        scale = max(scale, min(tree_only * 0.5, width_scale))
        scale_y = max(0.002, min(120.0, scale))
        scale_x = max(scale_y, min(120.0, width_scale))
        return scale_x, scale_y

    def _fit_to_viewport(self) -> None:
        """Fit every row in. Rectangular: the y scale fits the rows in the
        height, the x scale the tree's width (never below y; see
        _rectangular_fit_scales).
        Circular: the whole tree, undistorted. The camera zoom is 1 again."""
        if self._root is None or self._fitting:
            return
        viewport = self.viewport()
        if viewport is None or viewport.width() <= 0 or viewport.height() <= 0:
            return
        if not self._leaf_order:
            # A message (empty filter, no leaves) stays at its own size.
            self.resetTransform()
            self._camera_zoom = 1.0
            self.centerOn(self.sceneRect().center())
            return
        self._fitting = True
        try:
            # The fit measures the annotations and labels at their own size.
            self.resetTransform()
            self._camera_zoom = 1.0
            self._apply_camera_to_items()
            rect = self._fit_rect()
            if rect is None:
                return
            if self._layout_mode == "rectangular":
                scales = self._rectangular_fit_scales(rect)
                if scales is None:
                    return
                scale_x, scale = scales
                self._fit_scales = (scale_x, scale)
                self.scale(scale_x, scale)
                self._apply_camera_to_items()
                self._update_scene_rect_for_device_stable_items(scale_x, scale)
                available_width = max(1.0, float(viewport.width()) - 16.0)
                above, below = self._stable_vertical_extent(float(rect.top()), scale)
                extent_bottom = max(float(rect.height()) * scale, below)
                center_y = float(rect.top()) + ((extent_bottom - above) * 0.5) / scale
                center_x = rect.left() + (available_width * 0.5 / scale_x)
                # A scene smaller than the viewport would be centered by the
                # view whatever centerOn says.
                self._cover_viewport_at(QtCore.QPointF(center_x, center_y), scale_x, scale)
                self.centerOn(center_x, center_y)
            else:
                self.fitInView(rect, QtCore.Qt.AspectRatioMode.KeepAspectRatio)
                scale = abs(float(self.transform().m11()))
                self._fit_scales = (scale, scale)
        finally:
            self._fitting = False

    def fit_to_window(self) -> None:
        """Public: back to the fitted mode, the whole tree in the panel."""
        self._auto_fit = True
        self._fit_to_viewport()
        self.viewChanged.emit()

    @staticmethod
    def _wheel_delta(point: QtCore.QPoint) -> float:
        """The larger of a wheel delta's axes (a platform may turn the wheel
        sideways with Shift or Alt held)."""
        return float(point.y() if abs(point.y()) >= abs(point.x()) else point.x())

    def wheelEvent(self, event: QtGui.QWheelEvent) -> None:
        """Ctrl: camera zoom under the pointer; Ctrl+Shift: stretch the
        branch lengths; Alt: stretch the rows; Shift: scroll sideways;
        alone: scroll up and down (while fitted only within the fitted
        drawing; a wheel with nothing to scroll goes to the parent).
        Trackpads (pixel deltas) and wheels (angle deltas) both move smoothly."""
        event.accept()
        modifiers = event.modifiers()
        ctrl = bool(modifiers & QtCore.Qt.KeyboardModifier.ControlModifier)
        shift = bool(modifiers & QtCore.Qt.KeyboardModifier.ShiftModifier)
        alt = bool(modifiers & QtCore.Qt.KeyboardModifier.AltModifier)
        pixels = event.pixelDelta()
        if ctrl or alt:
            if not pixels.isNull():
                steps = self._wheel_delta(pixels) / _PIXELS_PER_ZOOM_STEP
            else:
                steps = self._wheel_delta(event.angleDelta()) / 120.0
            if steps == 0.0:
                return
            factor = _ZOOM_STEP ** steps
            anchor = event.position()
            if ctrl and shift:
                self.stretch_by(factor, 1.0, anchor)
            elif ctrl:
                self.zoom_by(factor, anchor)
            else:
                self.stretch_by(1.0, factor, anchor)
            return
        if not pixels.isNull():
            dx, dy = float(pixels.x()), float(pixels.y())
        else:
            # Lines per notch as the platform sets them, 20 pixels a line.
            line = 20.0 * max(1, QtWidgets.QApplication.wheelScrollLines()) / 120.0
            dx, dy = event.angleDelta().x() * line, event.angleDelta().y() * line
        if shift:
            dx, dy = (dy if abs(dy) >= abs(dx) else dx), 0.0
        fitted = self._auto_fit
        if not self._scroll_by(-dx, -dy, keep_fit=fitted) and fitted:
            event.ignore()

    def mousePressEvent(self, event: QtGui.QMouseEvent) -> None:
        # The items first (handles, chips, labels); a left or middle press
        # nothing took starts a pan.
        super().mousePressEvent(event)
        button = event.button()
        if button not in (QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.MouseButton.MiddleButton):
            return
        if event.isAccepted() and button == QtCore.Qt.MouseButton.LeftButton:
            return
        self._pan_from = event.position()
        self._pan_button = button
        self._panning = False
        event.accept()

    def _end_pan(self) -> None:
        if self._panning:
            self.viewport().unsetCursor()
        self._pan_from = None
        self._pan_button = QtCore.Qt.MouseButton.NoButton
        self._panning = False

    def mouseMoveEvent(self, event: QtGui.QMouseEvent) -> None:
        if self._pan_from is not None and not (event.buttons() & self._pan_button):
            # Its release went elsewhere (a menu, another window).
            self._end_pan()
        if self._pan_from is None:
            super().mouseMoveEvent(event)
            return
        moved = event.position() - self._pan_from
        if not self._panning:
            if moved.manhattanLength() < QtWidgets.QApplication.startDragDistance():
                event.accept()
                return
            self._panning = True
            self.viewport().setCursor(QtCore.Qt.CursorShape.ClosedHandCursor)
        self._pan_from = event.position()
        self._scroll_by(-moved.x(), -moved.y())
        event.accept()

    def mouseReleaseEvent(self, event: QtGui.QMouseEvent) -> None:
        if self._pan_from is not None and event.button() == self._pan_button:
            self._end_pan()
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def focusOutEvent(self, event: QtGui.QFocusEvent) -> None:
        self._end_pan()
        super().focusOutEvent(event)

    def mouseDoubleClickEvent(self, event: QtGui.QMouseEvent) -> None:
        # A left double-click no item takes (a triangle, a label) re-fits
        # the tree to the window.
        super().mouseDoubleClickEvent(event)
        if event.button() == QtCore.Qt.MouseButton.LeftButton and not event.isAccepted():
            self.fit_to_window()
            event.accept()

    def keyPressEvent(self, event: QtGui.QKeyEvent) -> None:
        """+ (or =) and - zoom, the arrows pan, 0 or F fits."""
        modifiers = event.modifiers() & ~(
            QtCore.Qt.KeyboardModifier.ShiftModifier | QtCore.Qt.KeyboardModifier.KeypadModifier
        )
        if modifiers != QtCore.Qt.KeyboardModifier.NoModifier:
            super().keyPressEvent(event)
            return
        key = QtCore.Qt.Key
        viewport = self.viewport()
        step_x = max(20.0, viewport.width() * 0.1)
        step_y = max(20.0, viewport.height() * 0.1)
        pans = {
            key.Key_Left: (-step_x, 0.0),
            key.Key_Right: (step_x, 0.0),
            key.Key_Up: (0.0, -step_y),
            key.Key_Down: (0.0, step_y),
        }
        pressed = event.key()
        if pressed in (key.Key_Plus, key.Key_Equal):
            self.zoom_in()
        elif pressed in (key.Key_Minus, key.Key_Underscore):
            self.zoom_out()
        elif pressed in (key.Key_0, key.Key_F):
            self.fit_to_window()
        elif pressed in pans:
            self._scroll_by(*pans[pressed])
        else:
            super().keyPressEvent(event)
            return
        event.accept()

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        super().resizeEvent(event)
        if self._auto_fit:
            self._fit_to_viewport()
        elif self._leaf_order:
            # The camera stays: the same zoom about the same center (the
            # view kept it; the new scene rect must not move it).
            viewport = self.viewport()
            middle = QtCore.QPointF(viewport.width() / 2.0, viewport.height() / 2.0)
            center = self._scene_point_at(middle)
            self._update_camera_scene_rect()
            self._put_scene_point_at(center, middle)
        self.viewChanged.emit()

    def scrollContentsBy(self, dx: int, dy: int) -> None:
        super().scrollContentsBy(dx, dy)
        if self._active_leaf_set:
            # The chips are placed in viewport pixels: a scroll blit would
            # carry them along with the tree.
            viewport = self.viewport()
            if viewport is not None:
                viewport.update()
        if not self._quiet_scroll:
            self.viewChanged.emit()
