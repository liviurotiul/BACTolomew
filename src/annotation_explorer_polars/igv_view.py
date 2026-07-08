from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional

from PyQt6 import QtCore, QtGui, QtWidgets


@dataclass
class GeneSegment:
    start: int
    end: int
    strand: str
    label: str
    color: QtGui.QColor
    lane: int = 0
    row_data: Optional[Dict[str, Any]] = None


class IgvTrackView(QtWidgets.QGraphicsView):
    geneDoubleClicked = QtCore.pyqtSignal(dict)
    viewCenterBpChanged = QtCore.pyqtSignal(int)

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self._scene = QtWidgets.QGraphicsScene(self)
        self.setScene(self._scene)
        self.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        self.setDragMode(QtWidgets.QGraphicsView.DragMode.ScrollHandDrag)
        self.setViewportUpdateMode(QtWidgets.QGraphicsView.ViewportUpdateMode.BoundingRectViewportUpdate)
        self.setTransformationAnchor(QtWidgets.QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setResizeAnchor(QtWidgets.QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.setMinimumHeight(210)

        self._lane_count = 3
        self._top_margin = 28
        self._bottom_margin = 20
        self._track_height = 24
        self._lane_gap = 12
        self._arrow_size = 8
        self._lane_labels = {
            0: "Ref loci",
            1: "Plus (+)",
            2: "Minus (-)",
        }
        self._mini_top_gap = 10
        self._mini_height = 46
        self._mini_inner_margin = 4
        self._mini_lane_gap = 2
        self._bp_per_px = 10.0
        self._contig_len = 0
        self._genes: List[GeneSegment] = []
        self._hit_bp_positions: List[int] = []
        self._center_animation: Optional[QtCore.QPropertyAnimation] = None
        self.horizontalScrollBar().valueChanged.connect(self._on_hscroll_changed)

    def clear(self) -> None:
        if self._center_animation is not None:
            self._center_animation.stop()
            self._center_animation = None
        self._scene.clear()
        self._contig_len = 0
        self._genes = []
        self._hit_bp_positions = []
        self.viewport().update()

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
        self._render_scene()
        QtCore.QTimer.singleShot(0, self._emit_view_center_bp_changed)

    def set_bp_per_px(self, bp_per_px: float) -> None:
        self._bp_per_px = max(bp_per_px, 0.1)
        self._render_scene()

    def _render_scene(self) -> None:
        self._scene.clear()
        if self._contig_len <= 0:
            return

        width = self._contig_len / self._bp_per_px
        height = (
            self._top_margin
            + self._lane_count * self._track_height
            + (self._lane_count - 1) * self._lane_gap
            + self._bottom_margin
        )
        self._scene.setSceneRect(0, 0, width, height)

        for lane in range(self._lane_count):
            lane_top = self._lane_top(lane)
            lane_center = lane_top + self._track_height / 2
            baseline = QtWidgets.QGraphicsLineItem(0, lane_center, width, lane_center)
            baseline.setPen(QtGui.QPen(QtGui.QColor("#c9cdd3"), 1))
            self._scene.addItem(baseline)

            label_text = self._lane_labels.get(lane, f"Lane {lane + 1}")
            label = QtWidgets.QGraphicsSimpleTextItem(label_text)
            label_font = QtGui.QFont()
            label_font.setPointSize(8)
            label_font.setBold(True)
            label.setFont(label_font)
            label.setBrush(QtGui.QBrush(QtGui.QColor("#5b6570")))
            label.setPos(4, max(0.0, lane_top - 14))
            label.setZValue(100)
            self._scene.addItem(label)

        for gene in self._genes:
            self._add_gene_item(gene)

        for bp in self._hit_bp_positions:
            self._add_hit_arrow(bp)
        self.viewport().update()

    def center_on_bp(self, bp: int) -> None:
        if self._contig_len <= 0:
            return
        bp_clamped = max(1, min(int(bp), self._contig_len))
        x = bp_clamped / self._bp_per_px
        center_y = self._lane_top(1) + self._track_height / 2
        self.centerOn(x, center_y)
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
        bp = int(round(scene_center.x()))
        bp = max(1, min(bp, self._contig_len))
        self.viewCenterBpChanged.emit(bp)

    def _on_hscroll_changed(self, _value: int) -> None:
        self._emit_view_center_bp_changed()
        self.viewport().update()

    def _lane_top(self, lane: int) -> float:
        lane_idx = max(0, min(int(lane), self._lane_count - 1))
        return self._top_margin + lane_idx * (self._track_height + self._lane_gap)

    def _add_gene_item(self, gene: GeneSegment) -> None:
        start = min(gene.start, gene.end)
        end = max(gene.start, gene.end)
        if end <= 0:
            return

        x0 = max(0.0, (start - 1) / self._bp_per_px)
        x1 = max(x0 + 1.0, end / self._bp_per_px)

        y = self._lane_top(getattr(gene, "lane", 0))
        total_w = x1 - x0
        arrow_w = min(float(self._arrow_size), total_w)
        body_rect: Optional[QtCore.QRectF] = None

        is_reverse = str(gene.strand).strip() == "-"
        if is_reverse:
            head_base_x = x0 + arrow_w
            body_w = max(0.0, x1 - head_base_x)
            if body_w > 0:
                body_rect = QtCore.QRectF(head_base_x, y, body_w, self._track_height)
            arrow = QtGui.QPolygonF(
                [
                    QtCore.QPointF(head_base_x, y),
                    QtCore.QPointF(x0, y + self._track_height / 2),
                    QtCore.QPointF(head_base_x, y + self._track_height),
                ]
            )
        else:
            head_base_x = x1 - arrow_w
            body_w = max(0.0, head_base_x - x0)
            if body_w > 0:
                body_rect = QtCore.QRectF(x0, y, body_w, self._track_height)
            arrow = QtGui.QPolygonF(
                [
                    QtCore.QPointF(head_base_x, y),
                    QtCore.QPointF(x1, y + self._track_height / 2),
                    QtCore.QPointF(head_base_x, y + self._track_height),
                ]
            )

        body_item: Optional[QtWidgets.QGraphicsRectItem] = None
        if body_rect is not None:
            body_item = QtWidgets.QGraphicsRectItem(body_rect)
            body_item.setBrush(QtGui.QBrush(gene.color))
            body_item.setPen(QtGui.QPen(QtGui.QColor("#333333"), 0.6))
            body_item.setToolTip(gene.label)
            if gene.row_data is not None:
                body_item.setData(0, gene.row_data)
            self._scene.addItem(body_item)

        label = gene.label or ""
        if label:
            font = QtGui.QFont()
            font.setPointSize(8)
            metrics = QtGui.QFontMetrics(font)
            label_rect = body_rect if body_rect is not None else QtCore.QRectF(x0, y, total_w, self._track_height)
            max_w = label_rect.width() - 4
            if max_w >= 6:
                text = metrics.elidedText(label, QtCore.Qt.TextElideMode.ElideRight, int(max_w))
                if text:
                    text_item = QtWidgets.QGraphicsSimpleTextItem(text, body_item)
                    text_item.setFont(font)
                    text_item.setBrush(QtGui.QBrush(QtGui.QColor("#1a1a1a")))
                    text_item.setPos(
                        label_rect.x() + 2,
                        label_rect.y() + (label_rect.height() - metrics.height()) / 2,
                    )
                    text_item.setZValue(25)

        # Arrowhead is always inside [x0, x1], so body + head equals full gene length.
        arrow_item = QtWidgets.QGraphicsPolygonItem(arrow)
        arrow_item.setBrush(QtGui.QBrush(gene.color))
        arrow_item.setPen(QtGui.QPen(QtGui.QColor("#333333"), 0.6))
        arrow_item.setToolTip(gene.label)
        if gene.row_data is not None:
            arrow_item.setData(0, gene.row_data)
        self._scene.addItem(arrow_item)

    def _add_hit_arrow(self, bp: int) -> None:
        if self._contig_len <= 0:
            return
        bp_clamped = max(1, min(int(bp), self._contig_len))
        x = bp_clamped / self._bp_per_px

        lane0_top = self._lane_top(0)
        gap_from_lane = 10.0
        shaft_len = 18.0
        head_half_w = 6.5
        head_h = 4.5
        tip_y = max(6.0, lane0_top - gap_from_lane)
        shaft_top = max(1.0, tip_y - shaft_len)

        line = QtWidgets.QGraphicsLineItem(x, shaft_top, x, tip_y - head_h)
        line.setPen(QtGui.QPen(QtGui.QColor("#000000"), 1.8))
        line.setZValue(300)
        self._scene.addItem(line)

        arrow = QtGui.QPolygonF(
            [
                QtCore.QPointF(x - head_half_w, tip_y - head_h),
                QtCore.QPointF(x + head_half_w, tip_y - head_h),
                QtCore.QPointF(x, tip_y + head_h),
            ]
        )
        arrow_item = QtWidgets.QGraphicsPolygonItem(arrow)
        arrow_item.setBrush(QtGui.QBrush(QtGui.QColor("#000000")))
        arrow_item.setPen(QtGui.QPen(QtGui.QColor("#000000"), 1.0))
        arrow_item.setZValue(300)
        self._scene.addItem(arrow_item)

    def drawForeground(self, painter: QtGui.QPainter, rect: QtCore.QRectF) -> None:
        super().drawForeground(painter, rect)
        if self._contig_len <= 0:
            return
        vrect = self.viewport().rect()
        if vrect.isEmpty():
            return

        mini_h = min(float(self._mini_height), max(20.0, float(vrect.height()) - 2.0))
        mini_rect = QtCore.QRectF(
            float(vrect.left()),
            float(vrect.bottom()) - mini_h + 1.0,
            float(vrect.width()),
            mini_h,
        )
        if mini_rect.width() <= 2 or mini_rect.height() <= 2:
            return

        painter.save()
        painter.resetTransform()
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)

        painter.setPen(QtGui.QPen(QtGui.QColor("#8f99a3"), 1.0))
        painter.setBrush(QtGui.QBrush(QtGui.QColor("#f5f7fa")))
        painter.drawRect(mini_rect)

        inner = mini_rect.adjusted(
            0.0,
            float(self._mini_inner_margin),
            0.0,
            -float(self._mini_inner_margin),
        )
        lane_gap = float(self._mini_lane_gap)
        lane_h = max(2.5, (inner.height() - 2 * lane_gap) / 3.0)
        lane_y = [inner.top() + i * (lane_h + lane_gap) for i in range(3)]

        painter.setPen(QtGui.QPen(QtGui.QColor("#d0d5dc"), 0.8))
        for y in lane_y:
            cy = y + lane_h / 2.0
            painter.drawLine(QtCore.QPointF(inner.left(), cy), QtCore.QPointF(inner.right(), cy))

        def _x_for_bp(bp: float) -> float:
            ratio = 0.0 if self._contig_len <= 0 else max(0.0, min(1.0, float(bp) / float(self._contig_len)))
            return inner.left() + ratio * inner.width()

        for gene in self._genes:
            lane = max(0, min(int(getattr(gene, "lane", 1)), 2))
            start_bp = float(min(gene.start, gene.end))
            end_bp = float(max(gene.start, gene.end))
            x0 = _x_for_bp(start_bp)
            x1 = _x_for_bp(end_bp)
            if x1 <= x0:
                x1 = min(inner.right(), x0 + 1.0)
            g_rect = QtCore.QRectF(x0, lane_y[lane] + 0.2, max(1.0, x1 - x0), max(1.0, lane_h - 0.4))
            painter.setPen(QtCore.Qt.PenStyle.NoPen)
            painter.setBrush(QtGui.QBrush(gene.color))
            painter.drawRect(g_rect)

        painter.setPen(QtGui.QPen(QtGui.QColor("#000000"), 1.0))
        for bp in self._hit_bp_positions:
            x = _x_for_bp(float(bp))
            painter.drawLine(QtCore.QPointF(x, mini_rect.top() + 1.0), QtCore.QPointF(x, mini_rect.bottom() - 1.0))

        scene_rect = self.sceneRect()
        if scene_rect.width() > 0:
            visible = self.mapToScene(vrect).boundingRect()
            left_ratio = max(0.0, min(1.0, (visible.left() - scene_rect.left()) / scene_rect.width()))
            right_ratio = max(0.0, min(1.0, (visible.right() - scene_rect.left()) / scene_rect.width()))
            if right_ratio < left_ratio:
                left_ratio, right_ratio = right_ratio, left_ratio
            vx0 = inner.left() + left_ratio * inner.width()
            vx1 = inner.left() + right_ratio * inner.width()
            if vx1 <= vx0:
                vx1 = min(inner.right(), vx0 + 1.0)
            viewport_rect = QtCore.QRectF(vx0, mini_rect.top(), max(1.0, vx1 - vx0), mini_rect.height())
            painter.setBrush(QtGui.QBrush(QtGui.QColor(0, 0, 0, 0)))
            painter.setPen(QtGui.QPen(QtGui.QColor("#20252a"), 1.6))
            painter.drawRect(viewport_rect)

        painter.restore()

    def _payload_from_item(self, item: Optional[QtWidgets.QGraphicsItem]) -> Optional[Dict[str, Any]]:
        cur = item
        while cur is not None:
            try:
                payload = cur.data(0)
            except Exception:
                payload = None
            if isinstance(payload, dict):
                return payload
            cur = cur.parentItem()
        return None

    def mouseDoubleClickEvent(self, event: QtGui.QMouseEvent) -> None:
        if self._scene is None:
            return super().mouseDoubleClickEvent(event)
        pos = event.position().toPoint() if hasattr(event, "position") else event.pos()
        scene_pos = self.mapToScene(pos)
        item = self._scene.itemAt(scene_pos, QtGui.QTransform())
        payload = self._payload_from_item(item)
        if payload is not None:
            self.geneDoubleClicked.emit(payload)
            event.accept()
            return
        super().mouseDoubleClickEvent(event)

    def wheelEvent(self, event: QtGui.QWheelEvent) -> None:
        if self._contig_len <= 0 or self._bp_per_px <= 0:
            return super().wheelEvent(event)
        delta = event.angleDelta().y()
        if delta == 0:
            return super().wheelEvent(event)

        view_pos = event.position().toPoint() if hasattr(event, "position") else event.pos()
        scene_pos = self.mapToScene(view_pos)
        bp_at_cursor = scene_pos.x() * self._bp_per_px

        steps = delta / 120.0
        # Faster + smoother zoom curve
        factor = 0.75 ** steps
        new_bp = self._bp_per_px * factor
        new_bp = max(0.5, min(new_bp, 500.0))
        if abs(new_bp - self._bp_per_px) < 1e-3:
            return
        self._bp_per_px = new_bp
        self._render_scene()

        new_scene_x = bp_at_cursor / self._bp_per_px
        hbar = self.horizontalScrollBar()
        hbar.setValue(int(new_scene_x - view_pos.x()))
        self.viewport().update()
        event.accept()

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        super().resizeEvent(event)
        self.viewport().update()
