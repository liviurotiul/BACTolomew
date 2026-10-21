"""Native Qt slippy tile map - the "web map" look (real OSM/Carto tiles, smooth
pan + deep zoom) with NO WebEngine and NO OpenGL, so it cannot crash or go black.

Countries are drawn as a choropleth (fill + coloured outline) by sample count,
overlaid on the tiles. Public API matches the other map views:
update_data(points, countries, max_count) + highlight_sample(s), and
countryHovered(key) with the normalized name of the country with samples under
the pointer ("" for none).
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from PyQt6 import QtCore, QtGui, QtNetwork, QtWidgets

from .geo_raster_map import _ASSETS, _count_color, _norm, ring_bounds, small_country_near
from .theme import COLORS

TILE = 256
BORDER_DARKEN = 135
HILITE = COLORS["hover"]
INK = "#334155"
TILE_URL = "https://a.basemaps.cartocdn.com/light_all/{z}/{x}/{y}.png"
ATTRIBUTION = "© OpenStreetMap  © CARTO"
MIN_Z, MAX_Z = 2, 18
_CACHE = Path.home() / ".cache" / "bactolomew_tiles"


def _lonlat_to_world(lon: float, lat: float, z: int) -> Tuple[float, float]:
    n = TILE * (2 ** z)
    lat = max(-85.0511, min(85.0511, float(lat)))
    x = (float(lon) + 180.0) / 360.0 * n
    s = math.sin(math.radians(lat))
    y = (0.5 - math.log((1 + s) / (1 - s)) / (4 * math.pi)) * n
    return x, y


def _world_to_lonlat(x: float, y: float, z: int) -> Tuple[float, float]:
    n = TILE * (2 ** z)
    lon = x / n * 360.0 - 180.0
    ty = math.pi * (1 - 2 * y / n)
    lat = math.degrees(math.atan(math.sinh(ty)))
    return lon, lat


_WORLD_LONLAT: Optional[List[Tuple[str, List[List[Tuple[float, float]]]]]] = None


def _load_world_lonlat() -> List[Tuple[str, List[List[Tuple[float, float]]]]]:
    """[(normalized_name, [ring_of_(lon,lat), ...]), ...], cached."""
    global _WORLD_LONLAT
    if _WORLD_LONLAT is not None:
        return _WORLD_LONLAT
    out: List[Tuple[str, List[List[Tuple[float, float]]]]] = []
    # Prefer high-resolution outlines so the choropleth edges match the crisp
    # borders in the map tiles (10m > 50m > 110m).
    for fname in ("world_countries_10m.geojson", "world_countries_50m.geojson", "world_countries.geojson"):
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
    _WORLD_LONLAT = out
    return out


class TileMapView(QtWidgets.QGraphicsView):
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
        self.setTransformationAnchor(QtWidgets.QGraphicsView.ViewportAnchor.NoAnchor)
        self.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.setBackgroundBrush(QtGui.QBrush(QtGui.QColor("#e7edf2")))
        self.setStyleSheet(
            "QScrollBar:horizontal,QScrollBar:vertical{background:transparent;width:9px;height:9px;}"
            "QScrollBar::handle{background:rgba(51,65,85,0.30);border-radius:4px;}"
            "QScrollBar::add-line,QScrollBar::sub-line{width:0;height:0;}"
        )

        self._z = MIN_Z
        self._points: List[dict] = []
        self._country_counts: Dict[str, int] = {}
        # Ring bounds of the drawn countries, for the hover slop of small
        # ones (built on the first pointer move after a rebuild).
        self._hit_rings: Optional[List[Tuple[str, QtCore.QRectF]]] = None
        self._sample_country: Dict[str, str] = {}
        self._max_count = 1
        self._highlight_countries: set[str] = set()
        self._tile_items: Dict[Tuple[int, int, int], QtWidgets.QGraphicsPixmapItem] = {}
        self._pending: set[Tuple[int, int, int]] = set()
        self._country_items: Dict[str, QtWidgets.QGraphicsPathItem] = {}
        self._did_initial = False
        self._hover_key = ""
        self.viewport().setMouseTracking(True)  # pointer moves without a button report the hovered country

        self._net = QtNetwork.QNetworkAccessManager(self)
        _CACHE.mkdir(parents=True, exist_ok=True)

        self._load_timer = QtCore.QTimer(self)
        self._load_timer.setSingleShot(True)
        self._load_timer.setInterval(60)
        self._load_timer.timeout.connect(self._load_visible)
        self.horizontalScrollBar().valueChanged.connect(self._schedule_load)
        self.verticalScrollBar().valueChanged.connect(self._schedule_load)

    # ---- public API ----
    def update_data(self, points: List[dict], countries: Dict[str, int], max_count: int) -> None:
        self._points = [p for p in (points or []) if _has_coords(p)]
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
        self._country_counts = counts
        self._max_count = max([int(max_count or 1), *counts.values()] or [1])
        self._did_initial = False
        self._rebuild()
        self._fit_view()

    def highlight_sample(self, sample: Optional[str]) -> None:
        self.highlight_samples([sample] if sample else [])

    def highlight_samples(self, samples: Optional[Sequence[str]]) -> None:
        self._highlight_countries = {
            self._sample_country[str(s).strip()]
            for s in (samples or [])
            if str(s).strip() in self._sample_country
        }
        self._restyle_countries()

    # ---- scene build ----
    def _rebuild(self) -> None:
        self._scene.clear()
        self._tile_items.clear()
        self._pending.clear()
        self._country_items = {}
        self._hit_rings = None
        n = TILE * (2 ** self._z)
        self._scene.setSceneRect(0, 0, n, n)
        self._place_countries()
        self._schedule_load()

    def _place_countries(self) -> None:
        for norm_name, rings in _load_world_lonlat():
            count = self._country_counts.get(norm_name, 0)
            if count <= 0:
                continue
            path = QtGui.QPainterPath()
            path.setFillRule(QtCore.Qt.FillRule.OddEvenFill)
            for ring in rings:
                x0, y0 = _lonlat_to_world(ring[0][0], ring[0][1], self._z)
                path.moveTo(x0, y0)
                for lon, lat in ring[1:]:
                    x, y = _lonlat_to_world(lon, lat, self._z)
                    path.lineTo(x, y)
                path.closeSubpath()
            item = self._scene.addPath(path)
            item.setZValue(3)
            self._country_items[norm_name] = item
        self._restyle_countries()

    def _restyle_countries(self) -> None:
        for norm_name, item in self._country_items.items():
            count = self._country_counts.get(norm_name, 0)
            base = _count_color(count, self._max_count)
            hot = norm_name in self._highlight_countries
            if hot:
                fill = QtGui.QColor(HILITE)
                fill.setAlpha(150)
                pen = QtGui.QPen(QtGui.QColor(HILITE).darker(120))
                pen.setWidthF(2.2)
            else:
                fill = QtGui.QColor(base)
                fill.setAlpha(70 if self._highlight_countries else 120)
                pen = QtGui.QPen(base.darker(BORDER_DARKEN))
                pen.setWidthF(1.5)
            pen.setCosmetic(True)
            pen.setJoinStyle(QtCore.Qt.PenJoinStyle.RoundJoin)
            item.setBrush(QtGui.QBrush(fill))
            item.setPen(pen)
            item.setZValue(6 if hot else 3)

    # ---- pointer hover ----
    def country_at(self, scene_pos: QtCore.QPointF) -> str:
        """The country with samples at this scene point, "" for none (only
        countries with samples have a path item)."""
        for name, item in self._country_items.items():
            path = item.path()
            if path.controlPointRect().contains(scene_pos) and path.contains(scene_pos):
                return name
        return ""

    def country_near(self, scene_pos: QtCore.QPointF) -> str:
        """The country under the point, else a small one a few pixels away
        (small_country_near)."""
        name = self.country_at(scene_pos)
        if name:
            return name
        if self._hit_rings is None:
            self._hit_rings = ring_bounds([(n, item.path()) for n, item in self._country_items.items()])
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

    # ---- tiles ----
    def _schedule_load(self) -> None:
        self._load_timer.start()

    def _load_visible(self) -> None:
        rect = self.mapToScene(self.viewport().rect()).boundingRect()
        n_tiles = 2 ** self._z
        x0 = max(0, int(rect.left() // TILE) - 1)
        x1 = min(n_tiles - 1, int(rect.right() // TILE) + 1)
        y0 = max(0, int(rect.top() // TILE) - 1)
        y1 = min(n_tiles - 1, int(rect.bottom() // TILE) + 1)
        for tx in range(x0, x1 + 1):
            for ty in range(y0, y1 + 1):
                key = (self._z, tx, ty)
                if key in self._tile_items or key in self._pending:
                    continue
                self._request_tile(key)

    def _request_tile(self, key: Tuple[int, int, int]) -> None:
        z, x, y = key
        cache_file = _CACHE / f"{z}_{x}_{y}.png"
        if cache_file.exists():
            pix = QtGui.QPixmap(str(cache_file))
            if not pix.isNull():
                self._place_tile(key, pix)
                return
        self._pending.add(key)
        req = QtNetwork.QNetworkRequest(QtCore.QUrl(TILE_URL.format(z=z, x=x, y=y)))
        req.setHeader(QtNetwork.QNetworkRequest.KnownHeaders.UserAgentHeader, "BACTolomewMapLab/1.0")
        reply = self._net.get(req)
        reply.finished.connect(lambda r=reply, k=key, f=cache_file: self._on_reply(r, k, f))

    def _on_reply(self, reply: QtNetwork.QNetworkReply, key: Tuple[int, int, int], cache_file: Path) -> None:
        self._pending.discard(key)
        try:
            if reply.error() == QtNetwork.QNetworkReply.NetworkError.NoError:
                data = reply.readAll()
                pix = QtGui.QPixmap()
                if pix.loadFromData(data) and not pix.isNull():
                    try:
                        cache_file.write_bytes(bytes(data))
                    except Exception:
                        pass
                    if key[0] == self._z:
                        self._place_tile(key, pix)
        finally:
            reply.deleteLater()

    def _place_tile(self, key: Tuple[int, int, int], pix: QtGui.QPixmap) -> None:
        if key in self._tile_items or key[0] != self._z:
            return
        item = self._scene.addPixmap(pix)
        item.setPos(key[1] * TILE, key[2] * TILE)
        item.setZValue(0)
        item.setTransformationMode(QtCore.Qt.TransformationMode.SmoothTransformation)
        self._tile_items[key] = item

    # ---- view / interaction ----
    def _set_zoom(self, z: int, anchor_lonlat: Optional[Tuple[float, float]] = None,
                  anchor_view: Optional[QtCore.QPoint] = None) -> None:
        z = max(MIN_Z, min(MAX_Z, z))
        if z == self._z:
            return
        self._z = z
        self._rebuild()
        if anchor_lonlat is not None:
            wx, wy = _lonlat_to_world(anchor_lonlat[0], anchor_lonlat[1], z)
            if anchor_view is not None:
                delta = anchor_view - QtCore.QPoint(self.viewport().width() // 2, self.viewport().height() // 2)
                self.centerOn(wx - delta.x(), wy - delta.y())
            else:
                self.centerOn(wx, wy)
        self._schedule_load()

    def wheelEvent(self, event: QtGui.QWheelEvent) -> None:
        delta = event.angleDelta().y() or event.pixelDelta().y()
        if delta == 0:
            return
        vp = event.position().toPoint()
        scene_pt = self.mapToScene(vp)
        lon, lat = _world_to_lonlat(scene_pt.x(), scene_pt.y(), self._z)
        self._set_zoom(self._z + (1 if delta > 0 else -1), (lon, lat), vp)
        event.accept()

    def _fit_view(self) -> None:
        pts = self._points
        if not pts:
            self._z = MIN_Z
            self._rebuild()
            self.centerOn(TILE * 2, TILE * 1.6)
            return
        lons = [p["lon"] for p in pts]
        lats = [p["lat"] for p in pts]
        min_lon, max_lon, min_lat, max_lat = min(lons), max(lons), min(lats), max(lats)
        vw = max(320, self.viewport().width())
        vh = max(240, self.viewport().height())
        best = MIN_Z
        for z in range(MIN_Z, 12):
            x0, y0 = _lonlat_to_world(min_lon, max_lat, z)
            x1, y1 = _lonlat_to_world(max_lon, min_lat, z)
            if (x1 - x0) < vw * 0.9 and (y1 - y0) < vh * 0.9:
                best = z
            else:
                break
        self._z = best
        self._rebuild()
        wx, wy = _lonlat_to_world((min_lon + max_lon) / 2, (min_lat + max_lat) / 2, best)
        self.centerOn(wx, wy)
        self._did_initial = True

    def showEvent(self, event: QtGui.QShowEvent) -> None:
        super().showEvent(event)
        if not self._did_initial and self._points:
            self._fit_view()
        else:
            self._schedule_load()

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        super().resizeEvent(event)
        self._schedule_load()

    def drawForeground(self, painter: QtGui.QPainter, rect: QtCore.QRectF) -> None:
        painter.save()
        painter.resetTransform()
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        vw, vh = self.viewport().width(), self.viewport().height()
        font = QtGui.QFont("Sans Serif", 8)
        painter.setFont(font)
        fm = QtGui.QFontMetrics(font)
        tw = fm.horizontalAdvance(ATTRIBUTION) + 12
        painter.setPen(QtCore.Qt.PenStyle.NoPen)
        painter.setBrush(QtGui.QColor(255, 255, 255, 190))
        painter.drawRect(QtCore.QRectF(vw - tw, vh - fm.height() - 4, tw, fm.height() + 4))
        painter.setPen(QtGui.QColor(INK))
        painter.drawText(int(vw - tw + 6), int(vh - 5), ATTRIBUTION)
        self._draw_legend(painter, vw, vh)
        painter.restore()

    def _draw_legend(self, painter: QtGui.QPainter, vw: int, vh: int) -> None:
        if not self._country_counts:
            return
        stops = sorted({1, max(2, self._max_count // 8), max(4, self._max_count // 4),
                        max(8, self._max_count // 2), self._max_count})
        sw, sh, gap = 26, 12, 2
        x0, y0 = 12, vh - 30
        font = QtGui.QFont("Sans Serif", 8)
        font.setBold(True)
        painter.setFont(font)
        fm = QtGui.QFontMetrics(font)
        total = len(stops) * (sw + gap)
        painter.setPen(QtCore.Qt.PenStyle.NoPen)
        painter.setBrush(QtGui.QColor(255, 255, 255, 205))
        painter.drawRoundedRect(QtCore.QRectF(x0 - 6, y0 - 8 - fm.height(), total + 10, sh + fm.height() + 12), 5, 5)
        painter.setPen(QtGui.QColor(INK))
        painter.drawText(x0, y0 - 8, "samples per country")
        for i, s in enumerate(stops):
            x = x0 + i * (sw + gap)
            painter.setPen(QtCore.Qt.PenStyle.NoPen)
            painter.setBrush(_count_color(s, self._max_count))
            painter.drawRect(QtCore.QRectF(x, y0, sw, sh))
            painter.setPen(QtGui.QColor(INK))
            label = str(s)
            painter.drawText(int(x + (sw - fm.horizontalAdvance(label)) / 2), int(y0 + sh + fm.ascent() + 1), label)


def _has_coords(p: dict) -> bool:
    try:
        lat = float(p["lat"])
        lon = float(p["lon"])
        return not (math.isnan(lat) or math.isnan(lon))
    except (KeyError, TypeError, ValueError):
        return False
