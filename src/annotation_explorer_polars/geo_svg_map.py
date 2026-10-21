"""SVG choropleth map: the offline/no-GPU fallback for the Geo panel.

A vector world choropleth (countries coloured by sample count, no points) rendered
as an SVG document inside a QGraphicsView, so it pans/zooms crisply with NO OpenGL
and NO network. Same API as the other maps: update_data(points, countries,
max_count) + highlight_sample(s), and countryHovered(key) with the
normalized name of the country with samples under the pointer ("" for none).
"""
from __future__ import annotations

import json
from typing import Dict, List, Optional, Sequence, Tuple

from PyQt6 import QtCore, QtGui, QtSvg, QtSvgWidgets, QtWidgets

from .geo_raster_map import (
    _ASSETS, _count_color, _norm, BORDER, INK, NO_DATA_FILL, OCEAN, ring_bounds, small_country_near,
)
from .theme import COLORS

_W, _H = 2000.0, 1000.0
HILITE = COLORS["hover"]

_WORLD: Optional[List[Tuple[str, List[List[Tuple[float, float]]]]]] = None
# Projected outlines, built once per process as plain Python data (no Qt
# objects at module level): each country's rings in map coordinates and its
# SVG path data.
_PROJECTED: Optional[Dict[str, List[List[Tuple[float, float]]]]] = None
_PATH_DATA: Optional[List[Tuple[str, str]]] = None


def _load_world() -> List[Tuple[str, List[List[Tuple[float, float]]]]]:
    global _WORLD
    if _WORLD is not None:
        return _WORLD
    out: List[Tuple[str, List[List[Tuple[float, float]]]]] = []
    # 50m gives a clean vector map without the weight of the 10m file.
    for fname in ("world_countries_50m.geojson", "world_countries.geojson", "world_countries_10m.geojson"):
        path = _ASSETS / fname
        if not path.exists():
            continue
        try:
            data = json.loads(path.read_text())
        except Exception:
            break
        for feat in data.get("features", []):
            props = feat.get("properties", {})
            name = props.get("NAME") or props.get("ADMIN") or props.get("SOVEREIGNT") or ""
            geom = feat.get("geometry") or {}
            gtype = geom.get("type")
            coords = geom.get("coordinates") or []
            polys = coords if gtype == "MultiPolygon" else ([coords] if gtype == "Polygon" else [])
            rings: List[List[Tuple[float, float]]] = []
            for poly in polys:
                for ring in poly:
                    if len(ring) >= 2:
                        rings.append([(float(x), float(y)) for x, y in ring])
            if rings:
                out.append((_norm(name), rings))
        break
    _WORLD = out
    return out


def _project(lon: float, lat: float) -> Tuple[float, float]:
    return (lon + 180.0) / 360.0 * _W, (90.0 - lat) / 180.0 * _H


def _projected() -> Dict[str, List[List[Tuple[float, float]]]]:
    """Each country's rings in map coordinates."""
    global _PROJECTED
    if _PROJECTED is not None:
        return _PROJECTED
    out: Dict[str, List[List[Tuple[float, float]]]] = {}
    for norm_name, rings in _load_world():
        out.setdefault(norm_name, []).extend([_project(lon, lat) for lon, lat in ring] for ring in rings)
    _PROJECTED = out
    return out


def _path_data() -> List[Tuple[str, str]]:
    """(country, SVG path data) for every country, in drawing order."""
    global _PATH_DATA
    if _PATH_DATA is not None:
        return _PATH_DATA
    out: List[Tuple[str, str]] = []
    for norm_name, rings in _load_world():
        d = []
        for ring in rings:
            x0, y0 = _project(ring[0][0], ring[0][1])
            d.append(f"M{x0:.1f} {y0:.1f}")
            for lon, lat in ring[1:]:
                x, y = _project(lon, lat)
                d.append(f"L{x:.1f} {y:.1f}")
            d.append("Z")
        out.append((norm_name, "".join(d)))
    _PATH_DATA = out
    return out


class SvgGeoView(QtWidgets.QGraphicsView):
    # The country under the pointer (normalized like the counts), "" when the
    # pointer is over no country with samples or leaves the map. Emitted only
    # when it changes.
    countryHovered = QtCore.pyqtSignal(str)

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self._scene = QtWidgets.QGraphicsScene(self)
        self.setScene(self._scene)
        self.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        self.setRenderHint(QtGui.QPainter.RenderHint.SmoothPixmapTransform, True)
        self.setDragMode(QtWidgets.QGraphicsView.DragMode.ScrollHandDrag)
        self.setTransformationAnchor(QtWidgets.QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.setBackgroundBrush(QtGui.QBrush(QtGui.QColor(OCEAN)))
        self.setStyleSheet(
            "QScrollBar:horizontal,QScrollBar:vertical{background:transparent;width:8px;height:8px;}"
            "QScrollBar::handle{background:rgba(31,45,61,0.28);border-radius:4px;}"
            "QScrollBar::add-line,QScrollBar::sub-line{width:0;height:0;}"
        )
        self._counts: Dict[str, int] = {}
        self._sample_country: Dict[str, str] = {}
        self._max = 1
        self._highlight: set[str] = set()
        # Hovered countries are one overlay path above the SVG, so a hover
        # never re-renders the map or resets the zoom. The scene owns the
        # item; a new one is made with each document (_rebuild).
        self._highlight_item: Optional[QtWidgets.QGraphicsPathItem] = None
        # Country outlines as QPainterPaths, built when first hovered.
        self._outlines: Dict[str, QtGui.QPainterPath] = {}
        self._item: Optional[QtSvgWidgets.QGraphicsSvgItem] = None
        self._renderer: Optional[QtSvg.QSvgRenderer] = None
        self._fit_scale = 1.0
        self._did_fit = False
        self._hover_key = ""
        self.viewport().setMouseTracking(True)  # pointer moves without a button report the hovered country
        # (country, outline, bounds) of the countries with samples, built on
        # the first pointer move after a data change.
        self._hit_paths: Optional[List[Tuple[str, QtGui.QPainterPath, QtCore.QRectF]]] = None
        # Their ring bounds, for the hover slop of small countries.
        self._hit_rings: Optional[List[Tuple[str, QtCore.QRectF]]] = None

    # ---- API ----
    def update_data(self, points: List[dict], countries: Dict[str, int], max_count: int,
                    total: Optional[int] = None) -> None:
        # total is accepted for parity with the web map; only the web map draws the on-map pie.
        counts: Dict[str, int] = {}
        for name, c in (countries or {}).items():
            counts[_norm(name)] = counts.get(_norm(name), 0) + int(c)
        self._sample_country = {}
        for p in points or []:
            country = _norm(p.get("country", ""))
            sample = str(p.get("sample") or "").strip()
            if country:
                counts[country] = counts.get(country, 0) + 1
                if sample:
                    self._sample_country[sample] = country
        self._counts = counts
        self._max = max([int(max_count or 1), *counts.values()] or [1])
        self._hit_paths = None
        self._hit_rings = None
        self._rebuild()

    def highlight_sample(self, sample: Optional[str]) -> None:
        self.highlight_samples([sample] if sample else [])

    def highlight_samples(self, samples: Optional[Sequence[str]]) -> None:
        highlight = {
            self._sample_country[str(s).strip()]
            for s in (samples or []) if str(s).strip() in self._sample_country
        }
        if highlight == self._highlight:
            return
        self._highlight = highlight
        self._draw_highlight()

    def highlighted_countries(self) -> set:
        return set(self._highlight)

    def _draw_highlight(self) -> None:
        item = self._highlight_item
        if item is None:
            return
        outline = QtGui.QPainterPath()
        for name in sorted(self._highlight):
            path = self._outline(name)
            if path is not None:
                outline.addPath(path)
        item.setPath(outline)

    def _outline(self, name: str) -> Optional[QtGui.QPainterPath]:
        path = self._outlines.get(name)
        if path is None:
            rings = _projected().get(name)
            if not rings:
                return None
            path = QtGui.QPainterPath()
            for ring in rings:
                path.addPolygon(QtGui.QPolygonF([QtCore.QPointF(x, y) for x, y in ring]))
                path.closeSubpath()
            self._outlines[name] = path
        return path

    # ---- pointer hover ----
    def country_at(self, scene_pos: QtCore.QPointF) -> str:
        """The country with samples at this scene point, "" for none."""
        if self._hit_paths is None:
            hits = []
            for name, count in self._counts.items():
                path = self._outline(name) if count > 0 else None
                if path is not None:
                    hits.append((name, path, path.boundingRect()))
            self._hit_paths = hits
        for name, path, bounds in self._hit_paths:
            if bounds.contains(scene_pos) and path.contains(scene_pos):
                return name
        return ""

    def country_near(self, scene_pos: QtCore.QPointF) -> str:
        """The country under the point, else a small one a few pixels away
        (Singapore, Fiji: small_country_near)."""
        name = self.country_at(scene_pos)
        if name:
            return name
        if self._hit_rings is None:
            self._hit_rings = ring_bounds([(n, path) for n, path, _bounds in self._hit_paths or []])
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

    # ---- rendering ----
    def _svg(self) -> str:
        parts = [
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {int(_W)} {int(_H)}">',
            f'<rect width="{int(_W)}" height="{int(_H)}" fill="{OCEAN}"/>',
        ]
        # Hovered countries are not part of the document (overlay path);
        # their outlines are projected here, once, so a hover is cheap.
        _projected()
        for norm_name, d in _path_data():
            count = self._counts.get(norm_name, 0)
            fill = _count_color(count, self._max).name() if count > 0 else NO_DATA_FILL
            parts.append(f'<path d="{d}" fill="{fill}" stroke="{BORDER}" stroke-width="0.6"/>')
        parts.append(self._legend_svg())
        parts.append("</svg>")
        return "".join(parts)

    def _legend_svg(self) -> str:
        stops = sorted({1, max(2, self._max // 8), max(4, self._max // 4), max(8, self._max // 2), self._max})
        sw, x0, y0 = 60, 40, _H - 70
        out = [f'<text x="{x0}" y="{y0 - 12}" font-family="sans-serif" font-size="26" font-weight="bold" fill="{INK}">samples per country</text>']
        for i, s in enumerate(stops):
            x = x0 + i * (sw + 4)
            out.append(f'<rect x="{x}" y="{y0}" width="{sw}" height="26" fill="{_count_color(s, self._max).name()}" stroke="{BORDER}"/>')
            out.append(f'<text x="{x + sw/2:.0f}" y="{y0 + 48}" font-family="sans-serif" font-size="22" fill="{INK}" text-anchor="middle">{s}</text>')
        return "".join(out)

    def _rebuild(self) -> None:
        self._scene.clear()
        # The scene owns the renderer (its items are deleted before it), so a
        # pending paint never meets a renderer that Python already freed.
        old, self._renderer = self._renderer, QtSvg.QSvgRenderer(bytearray(self._svg(), "utf-8"), self._scene)
        if old is not None:
            old.deleteLater()
        self._item = QtSvgWidgets.QGraphicsSvgItem()
        self._item.setSharedRenderer(self._renderer)
        self._scene.addItem(self._item)
        overlay = QtWidgets.QGraphicsPathItem()
        overlay.setPen(QtGui.QPen(QtGui.QColor(BORDER), 0.6))
        overlay.setBrush(QtGui.QBrush(QtGui.QColor(HILITE)))
        overlay.setZValue(1)
        self._scene.addItem(overlay)
        self._highlight_item = overlay
        self._scene.setSceneRect(0, 0, _W, _H)
        self._did_fit = False
        self._fit()
        self._draw_highlight()

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
        fit = self._fit_scale or 1.0
        if self.transform().m11() * factor < fit * 0.98 or self.transform().m11() * factor > fit * 60:
            return
        self.scale(factor, factor)

    def mouseDoubleClickEvent(self, event: QtGui.QMouseEvent) -> None:
        self._fit()
