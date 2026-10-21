"""The genome track: one contig's genes in three lanes (reference loci,
plus strand, minus strand), with an overview of the whole contig at the
bottom of the panel that works as its scroll bar.

It looks like the neighborhood plot (``context_canvas``): the same gene
arrows (``gene_arrow_path``), gene type colors (``context_model.gene_type_fill``,
set by ``IgvControllerMixin``), theme chrome, hover halo and selection outline,
and the same resistance mutation badge. Searched genes are marked by a pin
above the lanes with a faint guide down through them.

The lanes are centered vertically in the space above the overview; when the
panel is shorter than they need, the lanes and the overview get more compact
(down to ``MIN_LANE_PITCH`` / ``MIN_OVERVIEW_HEIGHT``) before the lanes scroll.
The overview (``_ContigOverview``) marks the visible range with a window:
drag it to pan, click elsewhere to center there, drag its edges to zoom.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Tuple

from PyQt6 import QtCore, QtGui, QtWidgets

from .context_canvas import (
    BASELINE_COLOR,
    CANVAS_BACKGROUND,
    GAP_BACKGROUND,
    GAP_TEXT_COLOR,
    HOVER_COLOR,
    HOVER_HALO_COLOR,
    MIN_TEXT_PX,
    MUTED_TEXT_COLOR,
    ROW_ALT_BACKGROUND,
    SELECTED_COLOR,
    TEXT_COLOR,
    badge_inset_px,
    gene_arrow_path,
)
from .resistance_mutation import mutation_note, paint_mutation_badge
from .theme import COLORS


# Segment kinds: a CDS (an arrow), a merged reference locus block, the joint
# between two members of a grouped plasmid, and a folded stretch of contig
# between two shown windows ("N bp context").
GENE = "gene"
LOCUS = "locus"
BREAK = "break"
GAP = "gap"

LANE_LABELS = ("Reference loci", "Plus strand", "Minus strand")
LANE_COUNT = len(LANE_LABELS)
RAIL_WIDTH = 96
# Natural lane geometry, and the compact limits used when the panel is short.
# The top margin holds the search pins (their halo reaches 26 px up).
TOP_MARGIN = 30.0
MIN_TOP_MARGIN = 26.0
LANE_PITCH = 36.0
MIN_LANE_PITCH = 21.0
LANE_PITCH_STEP = 3.0
GENE_HEIGHT = 18.0
BOTTOM_MARGIN = 6.0
MIN_BOTTOM_MARGIN = 2.0
# Height of the overview strip at the bottom of the panel.
OVERVIEW_HEIGHT = 40
MIN_OVERVIEW_HEIGHT = 26
# The least height of the whole track (the panel asks no more of its dock).
TRACK_MIN_HEIGHT = MIN_OVERVIEW_HEIGHT + 10
# Grab distance of the overview window's edges (zoom), in pixels.
OVERVIEW_EDGE_GRAB = 5.0
MIN_WINDOW_PX = 8.0
MIN_BP_PER_PX = 0.5
# A lane pitch change in a short panel redraws the contig; when that is slow,
# it waits until the panel height has settled (a splitter drag).
RELAYOUT_DELAY_MS = 150
BADGE_PX = 12.0
SEGMENT_OUTLINE = QtGui.QColor(COLORS["surface"])
SEGMENT_OUTLINE.setAlpha(215)
LOCUS_FILL = QtGui.QColor(COLORS["primary_soft"])
LOCUS_OUTLINE = QtGui.QColor(COLORS["primary"])
LOCUS_OUTLINE.setAlpha(150)
BREAK_FILL = QtGui.QColor(COLORS["border_strong"])
PIN_COLOR = SELECTED_COLOR
PIN_HALO = QtGui.QColor(COLORS["primary_soft"])
PIN_GUIDE = QtGui.QColor(COLORS["primary"])
PIN_GUIDE.setAlpha(70)


@dataclass
class GeneSegment:
    start: int
    end: int
    strand: str
    label: str
    # The fill of a GENE (the other kinds take theirs from the theme).
    color: QtGui.QColor = field(default_factory=QtGui.QColor)
    lane: int = 0
    row_data: Optional[Dict[str, Any]] = None
    kind: str = GENE


def _row_key(row: Optional[Dict[str, Any]]) -> Optional[Tuple[str, ...]]:
    """Identity of an annotation row, to find its gene again after a redraw."""
    if not row:
        return None
    uid = str(row.get("AnnUID") or "").strip()
    if uid:
        return ("uid", uid)
    return tuple(str(row.get(column) or "") for column in ("Sample", "Contig", "Start", "Stop"))


def _pin_path(x: float, tip_y: float) -> QtGui.QPainterPath:
    """A slim map pin whose point touches ``(x, tip_y)``."""
    radius = 6.5
    head_y = tip_y - 12.5
    head = QtGui.QPainterPath()
    head.addEllipse(QtCore.QPointF(x, head_y), radius, radius)
    point = QtGui.QPainterPath()
    point.moveTo(x - radius * 0.72, head_y + radius * 0.62)
    point.quadTo(x - 1.2, tip_y - 3.0, x, tip_y)
    point.quadTo(x + 1.2, tip_y - 3.0, x + radius * 0.72, head_y + radius * 0.62)
    point.closeSubpath()
    return head.united(point)


class _MutationBadgeItem(QtWidgets.QGraphicsItem):
    """The resistance mutation badge on a gene's lower right corner."""

    def __init__(self, center: QtCore.QPointF, note: str, row_index: Optional[int]) -> None:
        super().__init__()
        self._rect = QtCore.QRectF(-BADGE_PX / 2.0, -BADGE_PX / 2.0, BADGE_PX, BADGE_PX)
        self.setPos(center)
        self.setToolTip(note)
        if row_index is not None:
            self.setData(0, row_index)
        self.setZValue(25)

    def boundingRect(self) -> QtCore.QRectF:
        return self._rect

    def paint(
        self,
        painter: QtGui.QPainter,
        option: QtWidgets.QStyleOptionGraphicsItem,
        widget: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        del option, widget
        paint_mutation_badge(painter, self._rect)


class _LaneRail(QtWidgets.QWidget):
    """Lane names left of the track; they stay put while the track scrolls."""

    def __init__(self, view: "IgvTrackView") -> None:
        super().__init__(view)
        self.view = view
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_TransparentForMouseEvents, True)

    def paintEvent(self, event: QtGui.QPaintEvent) -> None:
        del event
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        painter.fillRect(self.rect(), CANVAS_BACKGROUND)
        font = QtGui.QFont()
        font.setPixelSize(11)
        painter.setFont(font)
        width = float(self.width())
        if self.view.has_track():
            pitch = self.view.lane_pitch()
            for lane, text in enumerate(LANE_LABELS):
                center = self.view.mapFromScene(QtCore.QPointF(0.0, self.view.lane_center(lane))).y()
                band = QtCore.QRectF(0.0, center - pitch / 2.0, width, pitch)
                if lane % 2:
                    painter.fillRect(band, ROW_ALT_BACKGROUND)
                painter.setPen(MUTED_TEXT_COLOR)
                painter.drawText(band.adjusted(12.0, 0.0, -6.0, 0.0), QtCore.Qt.AlignmentFlag.AlignVCenter, text)
        painter.setPen(QtGui.QPen(BASELINE_COLOR, 1.0))
        painter.drawLine(QtCore.QPointF(width - 0.5, 0.0), QtCore.QPointF(width - 0.5, float(self.height())))


class _ContigOverview(QtWidgets.QWidget):
    """The whole contig in miniature at the bottom of the track, used as its
    scroll bar: the window marks the visible range; drag it to pan, click
    elsewhere to center there, drag its left or right edge to zoom. The
    contig axis spans the same x range as the track's lanes."""

    PAN = "pan"
    LEFT = "left"
    RIGHT = "right"
    JUMP = "jump"

    def __init__(self, view: "IgvTrackView") -> None:
        super().__init__(view)
        self.view = view
        self.setObjectName("TrackOverview")
        self.setMouseTracking(True)
        self.setToolTip("Drag the window to move along the contig, click to jump, drag its edges to zoom")
        self._drag: Optional[str] = None
        # PAN: the bp under the cursor minus the window's left bp. LEFT and
        # RIGHT: the bp of the edge that stays put.
        self._drag_anchor = 0.0
        # An edge drag on a contig too slow to redraw live: the bp range the
        # window will show on release.
        self._pending: Optional[Tuple[float, float]] = None
        self._font = QtGui.QFont()
        self._font.setPixelSize(MIN_TEXT_PX)
        # The genes in miniature, drawn once per track and axis (a
        # chromosome has thousands): panning only moves the window over it.
        self._genes_pixmap: Optional[QtGui.QPixmap] = None
        self._genes_pixmap_key: Optional[Tuple[float, ...]] = None
        # The selected gene, looked up once per selection or track.
        self._selected_gene: Optional[GeneSegment] = None
        self._selected_known = False

    def invalidate(self) -> None:
        """The track's genes changed: redraw the miniature."""
        self._genes_pixmap = None
        self.selection_changed()

    def selection_changed(self) -> None:
        self._selected_known = False
        self.update()

    # Geometry --------------------------------------------------------------
    def axis(self) -> QtCore.QRectF:
        """The contig's x range and the mini lanes' height. It lies under the
        track's lanes: a contig shorter than the view spans only the width
        it is drawn at there."""
        left = float(RAIL_WIDTH) + 6.0
        right = float(RAIL_WIDTH + self.view.viewport().width()) - 6.0
        drawn = self.view.contig_width_px()
        if 0.0 < drawn < right - left:
            right = left + max(drawn, 2.0 * MIN_WINDOW_PX)
        return QtCore.QRectF(left, 6.0, max(1.0, right - left), max(4.0, self.height() - 11.0))

    def x_for_bp(self, bp: float) -> float:
        axis = self.axis()
        length = max(1, self.view.contig_length())
        return axis.left() + max(0.0, min(1.0, float(bp) / length)) * axis.width()

    def bp_at(self, x: float) -> float:
        axis = self.axis()
        return (float(x) - axis.left()) / axis.width() * max(1, self.view.contig_length())

    def window_rect(self) -> Optional[QtCore.QRectF]:
        visible = self._pending or self.view.visible_bp_range()
        if visible is None:
            return None
        axis = self.axis()
        x0 = self.x_for_bp(visible[0])
        x1 = self.x_for_bp(visible[1])
        if x1 - x0 < MIN_WINDOW_PX:
            # A long contig's window is widened around its middle to stay grabbable.
            middle = (x0 + x1) / 2.0
            x0, x1 = middle - MIN_WINDOW_PX / 2.0, middle + MIN_WINDOW_PX / 2.0
        return QtCore.QRectF(x0, axis.top() - 3.0, x1 - x0, axis.height() + 6.0)

    def _hit(self, x: float) -> Optional[str]:
        """What a press at ``x`` grabs: an edge of the window, its body or nothing."""
        window = self.window_rect()
        if window is None:
            return None
        x0, x1 = window.left(), window.right()
        grab = OVERVIEW_EDGE_GRAB
        # Edges are grabbed up to ``grab`` outside the window and as far
        # inside as leaves the middle third of a narrow window for panning.
        inside = min(grab, (x1 - x0) / 3.0) if x1 - x0 >= 12.0 else -0.5
        if x0 - grab <= x <= x0 + inside:
            return self.LEFT
        if x1 - inside <= x <= x1 + grab:
            return self.RIGHT
        if x0 <= x <= x1:
            return self.PAN
        axis = self.axis()
        if axis.left() <= x <= axis.right():
            return self.JUMP
        return None

    # Mouse -----------------------------------------------------------------
    def _cursor_for(self, part: Optional[str]) -> QtCore.Qt.CursorShape:
        if part in (self.LEFT, self.RIGHT):
            return QtCore.Qt.CursorShape.SizeHorCursor
        if part == self.PAN:
            return QtCore.Qt.CursorShape.ClosedHandCursor if self._drag else QtCore.Qt.CursorShape.OpenHandCursor
        if part == self.JUMP:
            return QtCore.Qt.CursorShape.PointingHandCursor
        return QtCore.Qt.CursorShape.ArrowCursor

    def mousePressEvent(self, event: QtGui.QMouseEvent) -> None:
        if event.button() != QtCore.Qt.MouseButton.LeftButton or not self.view.has_track():
            return super().mousePressEvent(event)
        x = event.position().x()
        part = self._hit(x)
        visible = self.view.visible_bp_range()
        if visible is None or part is None:
            # The rail column and the margins beside the axis do nothing.
            return super().mousePressEvent(event)
        if part == self.LEFT:
            self._drag_anchor = visible[1]
        elif part == self.RIGHT:
            self._drag_anchor = visible[0]
        else:
            if part == self.JUMP:
                # Jump: center the view here, then keep dragging the window.
                self.view.center_on_bp(int(round(self.bp_at(x))))
                visible = self.view.visible_bp_range() or visible
            part = self.PAN
            self._drag_anchor = self.bp_at(x) - visible[0]
        self._drag = part
        self.view.stop_center_animation()
        self.setCursor(self._cursor_for(part))
        event.accept()

    def mouseMoveEvent(self, event: QtGui.QMouseEvent) -> None:
        x = event.position().x()
        if self._drag is None:
            self.setCursor(self._cursor_for(self._hit(x)))
            return super().mouseMoveEvent(event)
        bp = max(0.0, min(float(self.view.contig_length()), self.bp_at(x)))
        if self._drag == self.PAN:
            self.view.scroll_to_left_bp(bp - self._drag_anchor)
            event.accept()
            return
        # An edge dragged past the other one stops at the closest zoom.
        length = float(self.view.contig_length())
        if self._drag == self.LEFT:
            span = self.view.clamp_bp_span(max(0.0, self._drag_anchor - bp))
            start, end = self._drag_anchor - span, self._drag_anchor
        else:
            span = self.view.clamp_bp_span(max(0.0, bp - self._drag_anchor))
            start, end = self._drag_anchor, self._drag_anchor + span
        # The zoom limit can push the range past a contig end: shift it back.
        shift = max(0.0, -start) - max(0.0, end - length)
        start, end = max(0.0, start + shift), min(length, end + shift)
        if self.view.redraws_fast():
            self.view.show_bp_range(start, end, keep=self.RIGHT if self._drag == self.LEFT else self.LEFT)
        else:
            self._pending = (start, end)
            self.update()
        event.accept()

    def mouseReleaseEvent(self, event: QtGui.QMouseEvent) -> None:
        if self._drag is None:
            return super().mouseReleaseEvent(event)
        pending, self._pending = self._pending, None
        if pending is not None:
            keep = self.RIGHT if self._drag == self.LEFT else self.LEFT
            self.view.show_bp_range(pending[0], pending[1], keep=keep)
        self._drag = None
        self.setCursor(self._cursor_for(self._hit(event.position().x())))
        event.accept()

    def wheelEvent(self, event: QtGui.QWheelEvent) -> None:
        delta = event.angleDelta().y()
        if delta == 0 or not self.view.has_track():
            return super().wheelEvent(event)
        self.view.zoom_by(delta / 120.0)
        event.accept()

    # Painting --------------------------------------------------------------
    def _gene_box(self, gene: GeneSegment) -> QtCore.QRectF:
        axis = self.axis()
        lane_gap = 1.5
        lane_h = max(1.5, (axis.height() - 2 * lane_gap) / LANE_COUNT)
        lane = max(0, min(int(gene.lane), LANE_COUNT - 1))
        x0 = self.x_for_bp(float(min(gene.start, gene.end)))
        x1 = max(self.x_for_bp(float(max(gene.start, gene.end))), x0 + 1.0)
        return QtCore.QRectF(x0, axis.top() + lane * (lane_h + lane_gap), x1 - x0, lane_h)

    def _draw_genes(self) -> QtGui.QPixmap:
        ratio = self.devicePixelRatioF()
        pixmap = QtGui.QPixmap(QtCore.QSize(max(1, int(self.width() * ratio)), max(1, int(self.height() * ratio))))
        pixmap.setDevicePixelRatio(ratio)
        pixmap.fill(QtCore.Qt.GlobalColor.transparent)
        painter = QtGui.QPainter(pixmap)
        painter.setPen(QtCore.Qt.PenStyle.NoPen)
        for gene in self.view.genes():
            if gene.kind == GENE:
                color = QtGui.QColor(gene.color)
            elif gene.kind == LOCUS:
                color = QtGui.QColor(LOCUS_OUTLINE)
            else:
                color = QtGui.QColor(BREAK_FILL)
            painter.setBrush(QtGui.QBrush(color))
            painter.drawRect(self._gene_box(gene))
        painter.end()
        return pixmap

    def _selected_box(self) -> Optional[QtCore.QRectF]:
        if not self._selected_known:
            key = self.view.selected_key()
            self._selected_gene = None
            if key is not None:
                self._selected_gene = next(
                    (g for g in self.view.genes() if g.kind == GENE and _row_key(g.row_data) == key), None
                )
            self._selected_known = True
        return None if self._selected_gene is None else self._gene_box(self._selected_gene)

    def paintEvent(self, event: QtGui.QPaintEvent) -> None:
        del event
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        painter.fillRect(self.rect(), CANVAS_BACKGROUND)
        painter.setPen(QtGui.QPen(BASELINE_COLOR, 1.0))
        painter.drawLine(QtCore.QPointF(0.0, 0.5), QtCore.QPointF(float(self.width()), 0.5))
        painter.drawLine(
            QtCore.QPointF(RAIL_WIDTH - 0.5, 0.0), QtCore.QPointF(RAIL_WIDTH - 0.5, float(self.height()))
        )
        if not self.view.has_track():
            return
        painter.setFont(self._font)
        painter.setPen(MUTED_TEXT_COLOR)
        painter.drawText(
            QtCore.QRectF(12.0, 0.0, RAIL_WIDTH - 18.0, float(self.height())),
            QtCore.Qt.AlignmentFlag.AlignVCenter,
            "Overview",
        )

        axis = self.axis()
        frame = axis.adjusted(-4.0, -3.0, 4.0, 3.0)
        painter.setPen(QtGui.QPen(QtGui.QColor(COLORS["border"]), 1.0))
        painter.setBrush(QtGui.QBrush(QtGui.QColor(COLORS["surface_subtle"])))
        painter.drawRoundedRect(frame, 4.0, 4.0)

        key = (self.width(), self.height(), axis.left(), axis.width(), axis.height())
        if self._genes_pixmap is None or self._genes_pixmap_key != key:
            self._genes_pixmap = self._draw_genes()
            self._genes_pixmap_key = key
        painter.drawPixmap(0, 0, self._genes_pixmap)
        selected = self._selected_box()

        # Searched genes: a pin head on the frame's top edge and a thin line down.
        for bp in self.view.hit_bp_positions():
            x = self.x_for_bp(float(bp))
            painter.setPen(QtGui.QPen(PIN_COLOR, 1.2))
            painter.drawLine(QtCore.QPointF(x, frame.top() + 2.0), QtCore.QPointF(x, frame.bottom() - 1.0))
            painter.setPen(QtGui.QPen(QtGui.QColor(COLORS["surface"]), 1.0))
            painter.setBrush(QtGui.QBrush(PIN_COLOR))
            painter.drawEllipse(QtCore.QPointF(x, frame.top() + 1.5), 3.2, 3.2)

        if selected is not None:
            painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
            painter.setPen(QtGui.QPen(SELECTED_COLOR, 1.6))
            painter.drawRoundedRect(selected.adjusted(-2.0, -2.0, 2.0, 2.0), 2.0, 2.0)

        window = self.window_rect()
        if window is not None:
            window = window.intersected(frame.adjusted(0.0, 1.0, 0.0, -1.0))
            soft = QtGui.QColor(COLORS["primary_soft"])
            soft.setAlpha(120)
            painter.setBrush(QtGui.QBrush(soft))
            painter.setPen(QtGui.QPen(SELECTED_COLOR, 1.4))
            painter.drawRoundedRect(window, 3.0, 3.0)
            if window.width() >= 14.0:
                # Grips on the edges say they can be dragged (zoom).
                grip = QtGui.QPen(SELECTED_COLOR, 1.6)
                grip.setCapStyle(QtCore.Qt.PenCapStyle.RoundCap)
                painter.setPen(grip)
                mid = window.center().y()
                for x in (window.left() + 3.0, window.right() - 3.0):
                    painter.drawLine(QtCore.QPointF(x, mid - 4.0), QtCore.QPointF(x, mid + 4.0))


class IgvTrackView(QtWidgets.QGraphicsView):
    geneDoubleClicked = QtCore.pyqtSignal(dict)
    viewCenterBpChanged = QtCore.pyqtSignal(int)

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self._scene = QtWidgets.QGraphicsScene(self)
        self.setScene(self._scene)
        self.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.setBackgroundBrush(QtGui.QBrush(CANVAS_BACKGROUND))
        self.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        self.setRenderHint(QtGui.QPainter.RenderHint.TextAntialiasing, True)
        self.setDragMode(QtWidgets.QGraphicsView.DragMode.ScrollHandDrag)
        self.setViewportUpdateMode(QtWidgets.QGraphicsView.ViewportUpdateMode.BoundingRectViewportUpdate)
        self.setTransformationAnchor(QtWidgets.QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        # resizeEvent keeps the middle of the view itself.
        self.setResizeAnchor(QtWidgets.QGraphicsView.ViewportAnchor.NoAnchor)
        self.setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignTop)
        # The overview is the horizontal scroll bar; the lanes scroll
        # vertically only when the panel is shorter than their compact size.
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.setMouseTracking(True)
        self.viewport().setMouseTracking(True)
        self._top_margin = TOP_MARGIN
        self._lane_pitch = LANE_PITCH
        self._bottom_margin = BOTTOM_MARGIN
        self._overview_height = OVERVIEW_HEIGHT
        self._pending_lanes = (self._top_margin, self._lane_pitch, self._bottom_margin)
        self._relayout_timer = QtCore.QTimer(self)
        self._relayout_timer.setSingleShot(True)
        self._relayout_timer.setInterval(RELAYOUT_DELAY_MS)
        self._relayout_timer.timeout.connect(self._apply_lane_layout)
        self._resizing = False
        self.setViewportMargins(RAIL_WIDTH, 0, 0, self._overview_height)
        self._rail = _LaneRail(self)
        self._overview = _ContigOverview(self)
        # The dock can be dragged very short: the overview stays whole at the
        # bottom and the lanes scroll above it.
        self.setMinimumHeight(TRACK_MIN_HEIGHT)

        self._bp_per_px = 10.0
        self._last_render_ms = 0.0
        self._contig_len = 0
        self._genes: List[GeneSegment] = []
        self._hit_bp_positions: List[int] = []
        self._center_animation: Optional[QtCore.QPropertyAnimation] = None
        # Gene arrows of the current scene with their segments (hover and
        # selection outlines follow their paths).
        self._gene_items: List[Tuple[QtWidgets.QGraphicsPathItem, GeneSegment]] = []
        # Annotation rows of the scene's items, by the index in their data(0).
        self._rows: List[Dict[str, Any]] = []
        self._hover_halo: Optional[QtWidgets.QGraphicsPathItem] = None
        self._hover_outline: Optional[QtWidgets.QGraphicsPathItem] = None
        self._selection_outline: Optional[QtWidgets.QGraphicsPathItem] = None
        self._hovered_index: Optional[int] = None
        self._selected_key: Optional[Tuple[str, ...]] = None
        self._label_font = QtGui.QFont()
        self._label_font.setPixelSize(MIN_TEXT_PX)
        self._label_metrics = QtGui.QFontMetricsF(self._label_font)
        self._gap_font = QtGui.QFont(self._label_font)
        self._gap_font.setItalic(True)
        self._gap_metrics = QtGui.QFontMetricsF(self._gap_font)
        self.horizontalScrollBar().valueChanged.connect(self._on_hscroll_changed)

    # Lane geometry ---------------------------------------------------------
    def lane_pitch(self) -> float:
        return self._lane_pitch

    def lane_center(self, lane: int) -> float:
        lane_idx = max(0, min(int(lane), LANE_COUNT - 1))
        return self._top_margin + lane_idx * self._lane_pitch + self._lane_pitch / 2.0

    def _lanes_bottom(self) -> float:
        return self._top_margin + LANE_COUNT * self._lane_pitch

    def _fit_to_height(self) -> None:
        """Size the overview strip and the lanes for the panel's height:
        natural sizes when they fit, compact ones (down to the minimums)
        when the panel is shorter; below that the lanes scroll."""
        height = self.contentsRect().height()
        natural = TOP_MARGIN + LANE_COUNT * LANE_PITCH + BOTTOM_MARGIN
        overview = int(max(MIN_OVERVIEW_HEIGHT, min(OVERVIEW_HEIGHT, height - natural)))
        lanes_height = float(height - overview)
        if lanes_height >= natural:
            top, pitch, bottom = TOP_MARGIN, LANE_PITCH, BOTTOM_MARGIN
        else:
            # The pitch goes down in steps (the spare pixels go to centering);
            # a slow contig is redrawn once the height settles.
            room = (lanes_height - MIN_TOP_MARGIN - MIN_BOTTOM_MARGIN) / LANE_COUNT
            pitch = max(MIN_LANE_PITCH, MIN_LANE_PITCH + LANE_PITCH_STEP * ((room - MIN_LANE_PITCH) // LANE_PITCH_STEP))
            top, bottom = MIN_TOP_MARGIN, MIN_BOTTOM_MARGIN
        if overview != self._overview_height:
            self._overview_height = overview
            self.setViewportMargins(RAIL_WIDTH, 0, 0, overview)
        self._pending_lanes = (top, pitch, bottom)
        if self._pending_lanes == (self._top_margin, self._lane_pitch, self._bottom_margin):
            self._relayout_timer.stop()
        elif self.has_track() and not self.redraws_fast():
            self._relayout_timer.start()
        else:
            self._apply_lane_layout()
        self._update_scene_rect()

    def _apply_lane_layout(self) -> None:
        """Redraw the lanes at the pitch ``_fit_to_height`` chose, keeping
        the scroll position."""
        self._relayout_timer.stop()
        if self._pending_lanes == (self._top_margin, self._lane_pitch, self._bottom_margin):
            return
        if not self.has_track():
            self._top_margin, self._lane_pitch, self._bottom_margin = self._pending_lanes
            return
        hbar = self.horizontalScrollBar()
        scroll = hbar.value()
        resizing, self._resizing = self._resizing, True
        try:
            self._render_scene()
            hbar.setValue(scroll)
            self._reset_vertical_scroll()
        finally:
            self._resizing = resizing

    def _reset_vertical_scroll(self) -> None:
        """When the lanes scroll: show the pins and the first two lanes
        if they fit, else the middle of the lanes."""
        vbar = self.verticalScrollBar()
        view_height = float(self.viewport().height())
        if view_height >= self._top_margin + 2 * self._lane_pitch:
            vbar.setValue(vbar.minimum())
            return
        top = self.lane_center(1) - view_height / 2.0 - self._scene.sceneRect().top()
        vbar.setValue(vbar.minimum() + int(round(top)))

    def _update_scene_rect(self) -> None:
        """Center the lanes vertically in the viewport when it is taller than
        they need (the pins above them stay in view); otherwise the scene is
        just the lanes and scrolls."""
        if self._contig_len <= 0:
            return
        content = self._lanes_bottom() + self._bottom_margin
        view_height = float(self.viewport().height())
        middle = self._top_margin + LANE_COUNT * self._lane_pitch / 2.0
        top = min(0.0, middle - view_height / 2.0)
        bottom = max(content, top + view_height)
        self._scene.setSceneRect(0.0, top, self._contig_len / self._bp_per_px, bottom - top)

    # Track data ------------------------------------------------------------
    def has_track(self) -> bool:
        return self._contig_len > 0

    def contig_length(self) -> int:
        return self._contig_len

    def contig_width_px(self) -> float:
        """The contig's drawn width in the lanes."""
        return self._contig_len / self._bp_per_px if self._contig_len > 0 else 0.0

    def genes(self) -> List[GeneSegment]:
        return self._genes

    def hit_bp_positions(self) -> List[int]:
        return self._hit_bp_positions

    def selected_key(self) -> Optional[Tuple[str, ...]]:
        return self._selected_key

    def clear(self) -> None:
        if self._center_animation is not None:
            self._center_animation.stop()
            self._center_animation = None
        self._reset_scene_items()
        self._scene.clear()
        self._contig_len = 0
        self._genes = []
        self._hit_bp_positions = []
        self.viewport().update()
        self._rail.update()
        self._overview.invalidate()

    def set_track(
        self,
        genes: Iterable[GeneSegment],
        contig_length: int,
        bp_per_px: float = 10.0,
        hit_bp_positions: Optional[Iterable[int]] = None,
    ) -> None:
        if self._center_animation is not None:
            self._center_animation.stop()
            self._center_animation = None
        self._genes = list(genes)
        self._bp_per_px = max(bp_per_px, 0.1)
        self._contig_len = max(int(contig_length), 0)
        self._hit_bp_positions = []
        if hit_bp_positions is not None:
            for bp in hit_bp_positions:
                try:
                    b = int(bp)
                except Exception:
                    continue
                if b > 0:
                    self._hit_bp_positions.append(b)
            if self._hit_bp_positions:
                self._hit_bp_positions = sorted(set(self._hit_bp_positions))
        self._overview.invalidate()
        self._render_scene()
        QtCore.QTimer.singleShot(0, self._emit_view_center_bp_changed)

    def set_bp_per_px(self, bp_per_px: float) -> None:
        self._bp_per_px = max(bp_per_px, 0.1)
        self._render_scene()

    def set_selected_row(self, row: Optional[Dict[str, Any]]) -> None:
        """Outline the gene of this annotation row (the table's current row)."""
        self._selected_key = _row_key(row)
        self._apply_selection()
        self._overview.selection_changed()

    def _reset_scene_items(self) -> None:
        self._gene_items = []
        self._rows = []
        self._hover_halo = None
        self._hover_outline = None
        self._selection_outline = None
        self._hovered_index = None

    def _render_scene(self) -> None:
        # Any redraw takes the lane pitch a deferred re-layout is waiting for.
        self._relayout_timer.stop()
        self._top_margin, self._lane_pitch, self._bottom_margin = self._pending_lanes
        started = time.perf_counter()
        self._render_scene_items()
        self._last_render_ms = (time.perf_counter() - started) * 1000.0

    def redraws_fast(self) -> bool:
        """Whether a zoom can redraw the track live (while an overview edge
        is dragged) instead of when the drag ends."""
        return self._last_render_ms <= 40.0

    def _render_scene_items(self) -> None:
        self._scene.clear()
        self._reset_scene_items()
        self._rail.update()
        self._overview.update()
        if self._contig_len <= 0:
            return
        self._update_scene_rect()

        for gene in self._genes:
            self._add_segment(gene)

        for bp in self._hit_bp_positions:
            self._add_hit_pin(bp)

        self._hover_halo = self._overlay_item(QtGui.QPen(HOVER_HALO_COLOR, 6.0), 20)
        self._hover_outline = self._overlay_item(QtGui.QPen(HOVER_COLOR, 1.8), 21)
        self._selection_outline = self._overlay_item(QtGui.QPen(SELECTED_COLOR, 2.2), 22)
        self._apply_selection()
        self.viewport().update()

    def drawBackground(self, painter: QtGui.QPainter, rect: QtCore.QRectF) -> None:
        """Lane bands and baselines across the whole viewport, so they run on
        past the end of a contig shorter than the view, as the rail's do."""
        painter.fillRect(rect, CANVAS_BACKGROUND)
        if self._contig_len <= 0:
            return
        painter.save()
        painter.setPen(QtCore.Qt.PenStyle.NoPen)
        pitch = self._lane_pitch
        for lane in range(1, LANE_COUNT, 2):
            center = self.lane_center(lane)
            painter.fillRect(
                QtCore.QRectF(rect.left(), center - pitch / 2.0, rect.width(), pitch),
                ROW_ALT_BACKGROUND,
            )
        pen = QtGui.QPen(BASELINE_COLOR, 1.0)
        pen.setCosmetic(True)
        painter.setPen(pen)
        for lane in range(LANE_COUNT):
            center = self.lane_center(lane)
            painter.drawLine(QtCore.QPointF(rect.left(), center), QtCore.QPointF(rect.right(), center))
        painter.restore()

    def selected_gene_bp(self) -> Optional[int]:
        """Middle of the selected gene, when it is drawn on this track."""
        if self._selected_key is None:
            return None
        for _item, gene in self._gene_items:
            if _row_key(gene.row_data) == self._selected_key:
                return int((gene.start + gene.end) / 2)
        return None

    def _overlay_item(self, pen: QtGui.QPen, z: float) -> QtWidgets.QGraphicsPathItem:
        pen.setCosmetic(True)
        pen.setJoinStyle(QtCore.Qt.PenJoinStyle.RoundJoin)
        item = QtWidgets.QGraphicsPathItem()
        item.setPen(pen)
        item.setBrush(QtGui.QBrush(QtCore.Qt.BrushStyle.NoBrush))
        item.setAcceptedMouseButtons(QtCore.Qt.MouseButton.NoButton)
        item.setAcceptHoverEvents(False)
        item.setZValue(z)
        item.hide()
        self._scene.addItem(item)
        return item

    def center_on_bp(self, bp: int) -> None:
        if self._contig_len <= 0:
            return
        target = self._target_horizontal_scroll_for_bp(bp)
        if target is not None:
            self.horizontalScrollBar().setValue(target)
        self._emit_view_center_bp_changed()

    def animate_center_on_bp(self, bp: int, duration_ms: int = 200) -> bool:
        target = self._target_horizontal_scroll_for_bp(bp)
        if target is None:
            return False
        hbar = self.horizontalScrollBar()
        current = hbar.value()
        target_i = int(target)
        if abs(current - target_i) <= 1:
            hbar.setValue(target_i)
            self._emit_view_center_bp_changed()
            return False

        if self._center_animation is not None:
            self._center_animation.stop()
            self._center_animation = None

        anim = QtCore.QPropertyAnimation(hbar, b"value", self)
        anim.setDuration(max(1, int(duration_ms)))
        anim.setStartValue(int(current))
        anim.setEndValue(target_i)
        anim.setEasingCurve(QtCore.QEasingCurve.Type.InOutCubic)

        def _cleanup() -> None:
            if self._center_animation is anim:
                self._center_animation = None
            self._emit_view_center_bp_changed()

        anim.finished.connect(_cleanup)
        self._center_animation = anim
        anim.start()
        return True

    def _target_horizontal_scroll_for_bp(self, bp: int) -> Optional[int]:
        if self._contig_len <= 0 or self._bp_per_px <= 0:
            return None
        bp_clamped = max(1, min(int(bp), self._contig_len))
        x = bp_clamped / self._bp_per_px
        viewport_rect = self.viewport().rect()
        half_width = viewport_rect.width() / 2.0 if not viewport_rect.isEmpty() else 0.0
        target = int(round(x - half_width))
        hbar = self.horizontalScrollBar()
        return max(hbar.minimum(), min(hbar.maximum(), target))

    def _emit_view_center_bp_changed(self) -> None:
        if self._contig_len <= 0:
            return
        viewport_rect = self.viewport().rect()
        if viewport_rect.isEmpty():
            return
        scene_center = self.mapToScene(viewport_rect.center())
        # Scene X is measured in rendered pixels, while consumers of this
        # signal expect biological coordinates.  Genes are placed at
        # ``bp / _bp_per_px`` in _add_segment(), so apply the inverse
        # conversion here before synchronising the annotation table.
        bp = int(round(scene_center.x() * self._bp_per_px))
        bp = max(1, min(bp, self._contig_len))
        self.viewCenterBpChanged.emit(bp)

    def _on_hscroll_changed(self, _value: int) -> None:
        if not self._resizing:
            self._emit_view_center_bp_changed()
        self.viewport().update()
        self._overview.update()

    def stop_center_animation(self) -> None:
        if self._center_animation is not None:
            self._center_animation.stop()
            self._center_animation = None

    def visible_bp_range(self) -> Optional[Tuple[float, float]]:
        """The contig range in view, in bp (the whole contig when it is
        shorter than the view)."""
        if self._contig_len <= 0:
            return None
        left = float(self.horizontalScrollBar().value())
        right = left + float(self.viewport().width())
        if self._contig_len / self._bp_per_px <= self.viewport().width():
            return 0.0, float(self._contig_len)
        return max(0.0, left * self._bp_per_px), min(float(self._contig_len), right * self._bp_per_px)

    def scroll_to_left_bp(self, bp: float) -> None:
        """Put ``bp`` at the left edge of the view (the overview's pan)."""
        if self._contig_len <= 0:
            return
        self.horizontalScrollBar().setValue(int(round(bp / self._bp_per_px)))

    def _max_bp_per_px(self) -> float:
        # Zooming out stops when the whole contig fits the view (a short
        # contig already drawn smaller than the view stays as it is).
        return max(self._bp_per_px, self._contig_len / max(1.0, float(self.viewport().width())))

    def clamp_bp_span(self, span: float) -> float:
        """The bp span the view can show, nearest to ``span`` (zoom limits)."""
        width = max(1.0, float(self.viewport().width()))
        return max(MIN_BP_PER_PX, min(float(span) / width, self._max_bp_per_px())) * width

    def show_bp_range(self, start_bp: float, end_bp: float, keep: str = "") -> None:
        """Zoom so ``start_bp`` to ``end_bp`` fills the view (the overview's
        edge drag). The scale is clamped; ``keep`` ("left" or "right") names
        the end that stays where it was when it is."""
        if self._contig_len <= 0:
            return
        width = max(1.0, float(self.viewport().width()))
        span = max(1.0, float(end_bp) - float(start_bp))
        bp_per_px = max(MIN_BP_PER_PX, min(span / width, self._max_bp_per_px()))
        left = float(end_bp) - bp_per_px * width if keep == _ContigOverview.RIGHT else float(start_bp)
        if abs(bp_per_px - self._bp_per_px) > 1e-6:
            self._bp_per_px = bp_per_px
            self._render_scene()
        self.scroll_to_left_bp(left)
        self.viewport().update()
        self._overview.update()

    def zoom_by(self, steps: float, anchor_x: Optional[float] = None) -> None:
        """Zoom in (positive ``steps``) or out, keeping the bp under the
        viewport x ``anchor_x`` (the view's middle by default) in place."""
        if self._contig_len <= 0 or self._bp_per_px <= 0:
            return
        if anchor_x is None:
            anchor_x = self.viewport().width() / 2.0
        bp_at_anchor = (self.horizontalScrollBar().value() + anchor_x) * self._bp_per_px
        new_bp = self._bp_per_px * (0.75 ** steps)
        new_bp = max(MIN_BP_PER_PX, min(new_bp, self._max_bp_per_px()))
        if abs(new_bp - self._bp_per_px) < 1e-3:
            return
        self._bp_per_px = new_bp
        self._render_scene()
        self.horizontalScrollBar().setValue(int(bp_at_anchor / self._bp_per_px - anchor_x))
        self.viewport().update()

    def _segment_rect(self, gene: GeneSegment) -> Optional[QtCore.QRectF]:
        start = min(gene.start, gene.end)
        end = max(gene.start, gene.end)
        if end <= 0:
            return None
        x0 = max(0.0, (start - 1) / self._bp_per_px)
        x1 = max(x0 + 1.0, end / self._bp_per_px)
        top = self.lane_center(gene.lane) - GENE_HEIGHT / 2.0
        return QtCore.QRectF(x0, top, x1 - x0, GENE_HEIGHT)

    def _add_segment(self, gene: GeneSegment) -> None:
        rect = self._segment_rect(gene)
        if rect is None:
            return
        tooltip = gene.label
        note = ""
        if gene.kind == GENE:
            path = gene_arrow_path(rect, "-" if str(gene.strand).strip() == "-" else "+", 1.0)
            fill = QtGui.QColor(gene.color)
            fill.setAlpha(245)
            pen = QtGui.QPen(SEGMENT_OUTLINE, 0.8)
            note = mutation_note(gene.row_data)
            tooltip = self._gene_tooltip(gene, note)
        else:
            path = QtGui.QPainterPath()
            path.addRoundedRect(rect, 3.0, 3.0)
            if gene.kind == LOCUS:
                fill, pen = QtGui.QColor(LOCUS_FILL), QtGui.QPen(LOCUS_OUTLINE, 1.0)
            elif gene.kind == GAP:
                fill = QtGui.QColor(GAP_BACKGROUND)
                pen = QtGui.QPen(BREAK_FILL, 1.0, QtCore.Qt.PenStyle.DashLine)
            else:
                fill, pen = QtGui.QColor(BREAK_FILL), QtGui.QPen(QtCore.Qt.PenStyle.NoPen)
        pen.setCosmetic(True)
        item = QtWidgets.QGraphicsPathItem(path)
        item.setBrush(QtGui.QBrush(fill))
        item.setPen(pen)
        item.setZValue(10)
        if tooltip:
            item.setToolTip(tooltip)
        row_index: Optional[int] = None
        if gene.row_data is not None:
            # An index into ``_rows``: a dict in item data is converted to a
            # QVariantMap, which made redrawing a chromosome much slower.
            row_index = len(self._rows)
            self._rows.append(gene.row_data)
            item.setData(0, row_index)
        self._scene.addItem(item)
        if gene.kind == GENE:
            item.setData(1, len(self._gene_items))
            self._gene_items.append((item, gene))
        self._add_segment_label(item, gene, rect, fill)
        if note:
            x = rect.right() - badge_inset_px(rect.width(), str(gene.strand).strip() != "-")
            self._scene.addItem(_MutationBadgeItem(QtCore.QPointF(x, rect.bottom()), note, row_index))

    @staticmethod
    def _gene_tooltip(gene: GeneSegment, note: str) -> str:
        row = gene.row_data or {}
        product = str(row.get("Product") or "").strip()
        lines = [gene.label or "Unnamed gene"]
        if product and product.lower() not in {"nan", gene.label.lower()}:
            lines.append(product)
        if note:
            lines.append(note)
        start, stop = row.get("Start"), row.get("Stop")
        try:
            lines.append(f"{int(start):,} to {int(stop):,} ({'minus' if str(gene.strand).strip() == '-' else 'plus'} strand)")
        except (TypeError, ValueError):
            pass
        lines.append("Double-click for the full gene card")
        return "\n".join(lines)

    def _add_segment_label(
        self,
        item: QtWidgets.QGraphicsPathItem,
        gene: GeneSegment,
        rect: QtCore.QRectF,
        fill: QtGui.QColor,
    ) -> None:
        label = gene.label or ""
        if not label or gene.kind == BREAK:
            return
        font = self._gap_font if gene.kind == GAP else self._label_font
        metrics = self._gap_metrics if gene.kind == GAP else self._label_metrics
        body = QtCore.QRectF(rect)
        if gene.kind == GENE:
            # The arrow head takes up to 8 px: center the name on the body.
            head = min(rect.width() * 0.42, 8.0)
            if str(gene.strand).strip() == "-":
                body.setLeft(rect.left() + head)
            else:
                body.setRight(rect.right() - head)
        available = body.width() - 6.0
        if available < 14.0:
            return
        text = metrics.elidedText(label, QtCore.Qt.TextElideMode.ElideRight, available)
        if not text or text == "\u2026":
            return
        if gene.kind == GAP:
            color = GAP_TEXT_COLOR
        elif gene.kind == GENE and fill.lightness() < 145:
            color = QtGui.QColor(COLORS["surface"])
        else:
            color = TEXT_COLOR
        text_item = QtWidgets.QGraphicsSimpleTextItem(text, item)
        text_item.setFont(font)
        text_item.setBrush(QtGui.QBrush(color))
        text_width = metrics.horizontalAdvance(text)
        text_item.setPos(
            body.center().x() - text_width / 2.0,
            body.center().y() - metrics.height() / 2.0,
        )

    def _add_hit_pin(self, bp: int) -> None:
        if self._contig_len <= 0:
            return
        x = max(1, min(int(bp), self._contig_len)) / self._bp_per_px
        tip_y = self._top_margin - 3.0
        # A faint guide behind the genes leads the eye from the pin down to
        # the gene, without covering its name.
        guide = self._scene.addLine(x, tip_y, x, self._lanes_bottom(), QtGui.QPen(PIN_GUIDE, 1.0))
        guide.setZValue(5)
        halo = self._scene.addEllipse(
            QtCore.QRectF(x - 10.5, tip_y - 23.0, 21.0, 21.0),
            QtGui.QPen(QtCore.Qt.PenStyle.NoPen),
            QtGui.QBrush(PIN_HALO),
        )
        halo.setZValue(29)
        pin = QtWidgets.QGraphicsPathItem(_pin_path(x, tip_y))
        pen = QtGui.QPen(QtGui.QColor(COLORS["surface"]), 1.2)
        pen.setCosmetic(True)
        pin.setPen(pen)
        pin.setBrush(QtGui.QBrush(PIN_COLOR))
        pin.setZValue(30)
        pin.setToolTip("Gene matching the search")
        self._scene.addItem(pin)
        dot = self._scene.addEllipse(
            QtCore.QRectF(x - 2.2, tip_y - 14.7, 4.4, 4.4),
            QtGui.QPen(QtCore.Qt.PenStyle.NoPen),
            QtGui.QBrush(QtGui.QColor(COLORS["surface"])),
        )
        dot.setZValue(31)

    def _apply_selection(self) -> None:
        outline = self._selection_outline
        if outline is None:
            return
        for item, gene in self._gene_items:
            if self._selected_key is not None and _row_key(gene.row_data) == self._selected_key:
                outline.setPath(item.path())
                outline.show()
                return
        outline.hide()

    def _set_hovered(self, index: Optional[int]) -> None:
        if index == self._hovered_index or self._hover_halo is None or self._hover_outline is None:
            return
        self._hovered_index = index
        if index is None or not 0 <= index < len(self._gene_items):
            self._hover_halo.hide()
            self._hover_outline.hide()
            return
        path = self._gene_items[index][0].path()
        for overlay in (self._hover_halo, self._hover_outline):
            overlay.setPath(path)
            overlay.show()

    def _payload_from_item(self, item: Optional[QtWidgets.QGraphicsItem]) -> Optional[Dict[str, Any]]:
        cur = item
        while cur is not None:
            index = cur.data(0)
            if isinstance(index, int) and 0 <= index < len(self._rows):
                return self._rows[index]
            cur = cur.parentItem()
        return None

    def _gene_index_at(self, pos: QtCore.QPoint) -> Optional[int]:
        for item in self.items(pos):
            cur: Optional[QtWidgets.QGraphicsItem] = item
            while cur is not None:
                index = cur.data(1)
                if isinstance(index, int):
                    return index
                cur = cur.parentItem()
        return None

    def mouseMoveEvent(self, event: QtGui.QMouseEvent) -> None:
        super().mouseMoveEvent(event)
        if event.buttons() == QtCore.Qt.MouseButton.NoButton:
            self._set_hovered(self._gene_index_at(event.position().toPoint()))

    def leaveEvent(self, event: QtCore.QEvent) -> None:
        self._set_hovered(None)
        super().leaveEvent(event)

    def mouseDoubleClickEvent(self, event: QtGui.QMouseEvent) -> None:
        pos = event.position().toPoint()
        payload = None
        for item in self.items(pos):
            payload = self._payload_from_item(item)
            if payload is not None:
                break
        if payload is not None:
            self.geneDoubleClicked.emit(payload)
            event.accept()
            return
        super().mouseDoubleClickEvent(event)

    def wheelEvent(self, event: QtGui.QWheelEvent) -> None:
        delta = event.angleDelta().y()
        if self._contig_len <= 0 or delta == 0:
            return super().wheelEvent(event)
        self.zoom_by(delta / 120.0, event.position().x())
        event.accept()

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        # QGraphicsView's own resize anchoring (and a re-layout for a new
        # lane pitch) can move the scroll position: keep the bp in the middle
        # of the view, and the lanes' top (the search pins) in view.
        # A resize is not a move along the contig, so it does not report
        # one (the paired table would follow the transient scroll).
        # The event is the viewport's (QAbstractScrollArea forwards it): the
        # viewport already has its new size, the scroll value is the old one.
        center_bp = None
        if self.has_track():
            old_width = event.oldSize().width()
            if old_width <= 0:
                old_width = self.viewport().width()
            center_bp = (self.horizontalScrollBar().value() + old_width / 2.0) * self._bp_per_px
        self._resizing = True
        try:
            super().resizeEvent(event)
            self._fit_to_height()
            if center_bp is not None:
                self.horizontalScrollBar().setValue(
                    int(round(center_bp / self._bp_per_px - self.viewport().width() / 2.0))
                )
            self._reset_vertical_scroll()
        finally:
            self._resizing = False
        area = self.contentsRect()
        self._rail.setGeometry(area.left(), area.top(), RAIL_WIDTH, self.viewport().height())
        self._overview.setGeometry(
            area.left(), area.bottom() + 1 - self._overview_height, area.width(), self._overview_height
        )
        self._rail.update()
        self._overview.update()
        self.viewport().update()
