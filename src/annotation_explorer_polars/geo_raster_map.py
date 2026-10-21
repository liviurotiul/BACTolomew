"""Dependency-free, OpenGL-free world map for the sample geographic view.

A QGraphicsView (raster viewport - never QOpenGLWidget, so it cannot hit the GL
crash) that draws countries as a choropleth by sample count plus per-sample
points, with a fixed title and colour-scale legend. Same public API as the other
map views: update_data(points, countries, max_count) + highlight_sample(s), and
countryHovered(key) with the normalized name of the country with samples under
the pointer ("" for none).
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from PyQt6 import QtCore, QtGui, QtWidgets

from .theme import COLORS

# ---- theme ----
NO_DATA_FILL = "#eef1f5"
BORDER = "#c2ccd6"
OCEAN = "#f4f8fc"
POINT_FILL = "#b11242"
POINT_EDGE = "#ffffff"
HILITE_FILL = COLORS["hover"]
INK = "#1f2d3d"
_RAMP_LO = (222, 235, 247)
_RAMP_HI = (8, 69, 148)

_ASSETS = Path(__file__).resolve().parent / "assets"
# scene coordinate space: equirectangular 2:1
_W = 2000.0
_H = 1000.0

_ALIASES = {
    "usa": "united states of america",
    "united states": "united states of america",
    "uk (scotland)": "united kingdom",
    "uk": "united kingdom",
    "scotland": "united kingdom",
    "great britain": "united kingdom",
}


def _norm(name: str) -> str:
    key = re.sub(r"[^a-z0-9 ]+", "", str(name or "").strip().lower()).strip()
    return _ALIASES.get(key, key)


def _project(lon: float, lat: float) -> Tuple[float, float]:
    return (float(lon) + 180.0) / 360.0 * _W, (90.0 - float(lat)) / 180.0 * _H


def _count_color(count: int, max_count: int) -> QtGui.QColor:
    if count <= 0 or max_count <= 0:
        return QtGui.QColor(NO_DATA_FILL)
    t = math.log1p(count) / math.log1p(max_count)
    t = max(0.0, min(1.0, t))
    r = round(_RAMP_LO[0] + (_RAMP_HI[0] - _RAMP_LO[0]) * t)
    g = round(_RAMP_LO[1] + (_RAMP_HI[1] - _RAMP_LO[1]) * t)
    b = round(_RAMP_LO[2] + (_RAMP_HI[2] - _RAMP_LO[2]) * t)
    return QtGui.QColor(r, g, b)


# A country (or island) drawn smaller than SMALL_RING_PX on screen is hovered
# from up to HOVER_SLOP_PX away (Singapore, Fiji at world zoom are a pixel or
# two); larger countries hover only under the pointer, so a border is exact.
HOVER_SLOP_PX = 4.0
SMALL_RING_PX = 8.0


def ring_bounds(entries: Sequence[Tuple[str, QtGui.QPainterPath]]) -> List[Tuple[str, QtCore.QRectF]]:
    """(country, bounds) of every ring of these outlines."""
    out: List[Tuple[str, QtCore.QRectF]] = []
    for name, path in entries:
        for polygon in path.toSubpathPolygons():
            out.append((name, polygon.boundingRect()))
    return out


def small_country_near(
    rings: Sequence[Tuple[str, QtCore.QRectF]], scene_pos: QtCore.QPointF, scale: float
) -> str:
    """The country with a ring drawn under SMALL_RING_PX at ``scale`` (view
    pixels per scene unit) within HOVER_SLOP_PX of the point, the nearest
    one; "" for none."""
    if scale <= 0.0:
        return ""
    slop = HOVER_SLOP_PX / scale
    small = SMALL_RING_PX / scale
    x, y = float(scene_pos.x()), float(scene_pos.y())
    best = ""
    best_d = math.inf
    for name, r in rings:
        if r.width() > small or r.height() > small:
            continue
        dx = max(r.left() - x, 0.0, x - r.right())
        dy = max(r.top() - y, 0.0, y - r.bottom())
        if dx > slop or dy > slop:
            continue
        d = math.hypot(dx, dy)
        if d < best_d:
            best, best_d = name, d
    return best

_WORLD_CACHE: Optional[List[Tuple[str, QtGui.QPainterPath]]] = None


def _load_world() -> List[Tuple[str, QtGui.QPainterPath]]:
    """[(normalized_name, QPainterPath in scene coords), ...], cached."""
    global _WORLD_CACHE
    if _WORLD_CACHE is not None:
        return _WORLD_CACHE
    path_file = None
    for name in ("world_countries.geojson", "world_countries_50m.geojson", "world_countries_10m.geojson"):
        if (_ASSETS / name).exists():
            path_file = _ASSETS / name
            break
    world: List[Tuple[str, QtGui.QPainterPath]] = []
    if path_file is not None:
        try:
            data = json.loads(path_file.read_text())
        except Exception:
            data = {"features": []}
        for feat in data.get("features", []):
            props = feat.get("properties", {})
            name = props.get("NAME") or props.get("ADMIN") or props.get("SOVEREIGNT") or ""
            geom = feat.get("geometry") or {}
            gtype = geom.get("type")
            coords = geom.get("coordinates") or []
            polys = coords if gtype == "MultiPolygon" else ([coords] if gtype == "Polygon" else [])
            qpath = QtGui.QPainterPath()
            qpath.setFillRule(QtCore.Qt.FillRule.OddEvenFill)
            for poly in polys:
                for ring in poly:
                    if len(ring) < 2:
                        continue
                    x0, y0 = _project(ring[0][0], ring[0][1])
                    qpath.moveTo(x0, y0)
                    for lon, lat in ring[1:]:
                        x, y = _project(lon, lat)
                        qpath.lineTo(x, y)
                    qpath.closeSubpath()
            if not qpath.isEmpty():
                world.append((_norm(name), qpath))
    _WORLD_CACHE = world
    return world


class RasterGeoMapView(QtWidgets.QGraphicsView):
    # The country under the pointer (normalized like the counts), "" when the
    # pointer is over no country with samples or leaves the map. Emitted only
    # when it changes.
    countryHovered = QtCore.pyqtSignal(str)

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self._scene = QtWidgets.QGraphicsScene(self)
        self._scene.setSceneRect(0, 0, _W, _H)
        self.setScene(self._scene)
        # raster viewport only (QGraphicsView default) - never QOpenGLWidget.
        self.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        self.setRenderHint(QtGui.QPainter.RenderHint.TextAntialiasing, True)
        self.setDragMode(QtWidgets.QGraphicsView.DragMode.ScrollHandDrag)
        self.setTransformationAnchor(QtWidgets.QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setResizeAnchor(QtWidgets.QGraphicsView.ViewportAnchor.AnchorViewCenter)
        # Scrollbars must be AsNeeded (not off) or ScrollHandDrag can't pan;
        # keep them thin/subtle rather than hidden.
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.setInteractive(True)
        self.setMouseTracking(True)
        self.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.setBackgroundBrush(QtGui.QBrush(QtGui.QColor(OCEAN)))
        self.setStyleSheet(
            "QScrollBar:horizontal, QScrollBar:vertical { background: transparent; width: 8px; height: 8px; }"
            "QScrollBar::handle { background: rgba(31,45,61,0.28); border-radius: 4px; }"
            "QScrollBar::add-line, QScrollBar::sub-line { width: 0; height: 0; }"
        )
        self._fit_scale = 1.0

        self._points: List[dict] = []
        self._countries: Dict[str, int] = {}
        self._max_count: int = 1
        self._highlight: set[str] = set()
        self._point_items: List[Tuple[str, QtWidgets.QGraphicsEllipseItem]] = []
        self._did_fit = False
        self._hover_key = ""
        self.viewport().setMouseTracking(True)  # pointer moves without a button report the hovered country
        # (country, outline, bounds) of the countries with samples (_draw_world).
        self._hit_paths: List[Tuple[str, QtGui.QPainterPath, QtCore.QRectF]] = []
        # Ring bounds of the countries with samples (hover slop), built on
        # the first pointer move after a rebuild.
        self._hit_rings: Optional[List[Tuple[str, QtCore.QRectF]]] = None

        self._border_pen = QtGui.QPen(QtGui.QColor(BORDER))
        self._border_pen.setWidthF(0.6)
        self._border_pen.setCosmetic(True)
        self._border_pen.setJoinStyle(QtCore.Qt.PenJoinStyle.RoundJoin)

        self._draw_world(empty=True)

    # ---- public API (matches the other map views) ----
    def update_data(self, points: List[dict], countries: Dict[str, int], max_count: int) -> None:
        self._points = list(points or [])
        self._countries = {str(k): int(v) for k, v in (countries or {}).items()}
        # scale the choropleth to the biggest of country-counts / passed max
        self._max_count = max([int(max_count or 1), *self._countries.values()] or [1])
        self._draw_world(empty=False)

    def highlight_sample(self, sample: Optional[str]) -> None:
        self.highlight_samples([sample] if sample else [])

    def highlight_samples(self, samples: Optional[Sequence[str]]) -> None:
        self._highlight = {str(s).strip() for s in (samples or []) if str(s).strip()}
        self._restyle_points()

    # ---- rendering ----
    def _draw_world(self, empty: bool) -> None:
        self._scene.clear()
        self._point_items = []
        self._hit_paths = []
        self._hit_rings = None
        counts_norm: Dict[str, int] = {}
        for name, count in self._countries.items():
            counts_norm[_norm(name)] = counts_norm.get(_norm(name), 0) + int(count)
        # A country with sample points is hoverable too (its points, not its fill, show them).
        hoverable = set() if empty else {key for key, count in counts_norm.items() if count > 0}
        if not empty:
            hoverable.update(_norm(pt.get("country", "")) for pt in self._points)

        for norm_name, qpath in _load_world():
            count = 0 if empty else counts_norm.get(norm_name, 0)
            item = self._scene.addPath(qpath, self._border_pen, QtGui.QBrush(_count_color(count, self._max_count)))
            item.setZValue(0)
            if norm_name in hoverable:
                self._hit_paths.append((norm_name, qpath, qpath.boundingRect()))

        if not empty:
            for pt in self._points:
                try:
                    x, y = _project(pt["lon"], pt["lat"])
                except (KeyError, TypeError, ValueError):
                    continue
                dot = self._scene.addEllipse(-4.0, -4.0, 8.0, 8.0, self._point_pen(), QtGui.QBrush(QtGui.QColor(POINT_FILL)))
                dot.setPos(x, y)
                dot.setZValue(5)
                dot.setFlag(QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations, True)
                sample = str(pt.get("sample") or "").strip()
                dot.setToolTip(f"{sample}  ·  {pt.get('country', '')}".strip(" ·"))
                self._point_items.append((sample, dot))
            self._restyle_points()

        self._did_fit = False
        self._fit()

    # ---- pointer hover ----
    def country_at(self, scene_pos: QtCore.QPointF) -> str:
        """The country with samples at this scene point, "" for none."""
        for name, path, bounds in self._hit_paths:
            if bounds.contains(scene_pos) and path.contains(scene_pos):
                return name
        return ""

    def country_near(self, scene_pos: QtCore.QPointF) -> str:
        """The country under the point, else a small one a few pixels away
        (small_country_near)."""
        name = self.country_at(scene_pos)
        if name:
            return name
        if self._hit_rings is None:
            self._hit_rings = ring_bounds([(n, path) for n, path, _bounds in self._hit_paths])
        return small_country_near(self._hit_rings, scene_pos, abs(float(self.transform().m11())))

    def reset_hover(self) -> None:
        """Forget the hovered country without emitting (the panel already
        cleared its hover), so the next pointer move over it reports it again."""
        self._hover_key = ""

    def _set_hover_key(self, key: str) -> None:
        if key != self._hover_key:
            self._hover_key = key
            self.countryHovered.emit(key)

    def mouseMoveEvent(self, event: QtGui.QMouseEvent) -> None:
        super().mouseMoveEvent(event)
        if event.buttons() == QtCore.Qt.MouseButton.NoButton:  # not while dragging the map
            self._set_hover_key(self.country_near(self.mapToScene(event.position().toPoint())))

    def viewportEvent(self, event: QtCore.QEvent) -> bool:
        if event.type() == QtCore.QEvent.Type.Leave:  # the pointer left the map
            self._set_hover_key("")
        return super().viewportEvent(event)

    @staticmethod
    def _point_pen() -> QtGui.QPen:
        pen = QtGui.QPen(QtGui.QColor(POINT_EDGE))
        pen.setWidthF(1.1)
        return pen

    def _restyle_points(self) -> None:
        for sample, item in self._point_items:
            hot = sample in self._highlight
            if self._highlight and hot:
                item.setRect(-6.0, -6.0, 12.0, 12.0)
                item.setBrush(QtGui.QBrush(QtGui.QColor(HILITE_FILL)))
                item.setZValue(9)
            else:
                item.setRect(-4.0, -4.0, 8.0, 8.0)
                item.setBrush(QtGui.QBrush(QtGui.QColor(POINT_FILL)))
                item.setZValue(5)
            if self._highlight and not hot:
                item.setOpacity(0.45)
            else:
                item.setOpacity(1.0)

    def _fit(self) -> None:
        self.resetTransform()
        self.fitInView(self._scene.sceneRect(), QtCore.Qt.AspectRatioMode.KeepAspectRatio)
        self._fit_scale = self.transform().m11() or 1.0
        self._did_fit = True

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        super().resizeEvent(event)
        if not self._did_fit:
            self._fit()

    def showEvent(self, event: QtGui.QShowEvent) -> None:
        super().showEvent(event)
        self._fit()

    def wheelEvent(self, event: QtGui.QWheelEvent) -> None:
        delta = event.angleDelta().y() or event.pixelDelta().y()
        if delta == 0:
            return
        factor = 1.2 if delta > 0 else 1 / 1.2
        cur = self.transform().m11()
        fit = self._fit_scale or cur or 1.0
        target = cur * factor
        # allow zooming from the fitted view up to ~50x; never further out than fit
        if target < fit * 0.98 or target > fit * 50:
            event.accept()
            return
        self.scale(factor, factor)
        event.accept()

    def mouseDoubleClickEvent(self, event: QtGui.QMouseEvent) -> None:
        self._fit()

    # fixed title + legend overlay (viewport pixels, unaffected by pan/zoom)
    def drawForeground(self, painter: QtGui.QPainter, rect: QtCore.QRectF) -> None:
        painter.save()
        painter.resetTransform()
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        vw = self.viewport().width()
        vh = self.viewport().height()
        ink = QtGui.QColor(INK)

        title_font = QtGui.QFont("Sans Serif", 12)
        title_font.setBold(True)
        painter.setFont(title_font)
        painter.setPen(ink)
        painter.drawText(14, 24, "Sample geographic distribution")

        if self._max_count > 0 and self._countries:
            self._draw_legend(painter, vw, vh, ink)
        painter.restore()

    def _draw_legend(self, painter: QtGui.QPainter, vw: int, vh: int, ink: QtGui.QColor) -> None:
        stops = sorted({1, max(2, self._max_count // 8), max(4, self._max_count // 4),
                        max(8, self._max_count // 2), self._max_count})
        sw, sh, gap = 30, 13, 2
        n = len(stops)
        total = n * sw + (n - 1) * gap
        x0 = vw - total - 16
        y0 = vh - 34
        font = QtGui.QFont("Sans Serif", 8)
        font.setBold(True)
        painter.setFont(font)
        fm = QtGui.QFontMetrics(font)
        pad = 6
        panel = QtCore.QRectF(x0 - pad, y0 - 8 - fm.height(), total + 2 * pad, sh + fm.height() + 8 + 2 * pad)
        painter.setPen(QtCore.Qt.PenStyle.NoPen)
        painter.setBrush(QtGui.QColor(255, 255, 255, 210))
        painter.drawRoundedRect(panel, 5, 5)
        painter.setPen(ink)
        painter.drawText(int(x0), int(y0 - 8), "samples per country")
        swatch_pen = QtGui.QPen(QtGui.QColor(BORDER))
        swatch_pen.setWidthF(0.6)
        for i, s in enumerate(stops):
            x = x0 + i * (sw + gap)
            painter.setPen(swatch_pen)
            painter.setBrush(_count_color(s, self._max_count))
            painter.drawRect(QtCore.QRectF(x, y0, sw, sh))
            painter.setPen(ink)
            label = str(s)
            tw = fm.horizontalAdvance(label)
            painter.drawText(int(x + (sw - tw) / 2), int(y0 + sh + fm.ascent() + 2), label)
