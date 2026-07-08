"""SVG choropleth map — the offline/no-GPU fallback for the Geo panel.

A vector world choropleth (countries coloured by sample count, no points) rendered
as an SVG document inside a QGraphicsView, so it pans/zooms crisply with NO OpenGL
and NO network. Same API as the other maps: update_data(points, countries,
max_count) + highlight_sample(s).
"""
from __future__ import annotations

import json
from typing import Dict, List, Optional, Sequence, Tuple

from PyQt6 import QtCore, QtGui, QtSvg, QtSvgWidgets, QtWidgets

from .geo_raster_map import _ASSETS, _count_color, _norm, BORDER, INK, NO_DATA_FILL, OCEAN

_W, _H = 2000.0, 1000.0
HILITE = "#f59e0b"

_WORLD: Optional[List[Tuple[str, List[List[Tuple[float, float]]]]]] = None


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


class SvgGeoView(QtWidgets.QGraphicsView):
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
        self._item: Optional[QtSvgWidgets.QGraphicsSvgItem] = None
        self._renderer: Optional[QtSvg.QSvgRenderer] = None
        self._fit_scale = 1.0
        self._did_fit = False

    # ---- API ----
    def update_data(self, points: List[dict], countries: Dict[str, int], max_count: int) -> None:
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
        self._rebuild()

    def highlight_sample(self, sample: Optional[str]) -> None:
        self.highlight_samples([sample] if sample else [])

    def highlight_samples(self, samples: Optional[Sequence[str]]) -> None:
        self._highlight = {
            self._sample_country[str(s).strip()]
            for s in (samples or []) if str(s).strip() in self._sample_country
        }
        self._rebuild()

    # ---- rendering ----
    def _svg(self) -> str:
        parts = [
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {int(_W)} {int(_H)}">',
            f'<rect width="{int(_W)}" height="{int(_H)}" fill="{OCEAN}"/>',
        ]
        for norm_name, rings in _load_world():
            count = self._counts.get(norm_name, 0)
            if norm_name in self._highlight:
                fill = HILITE
            elif count > 0:
                fill = _count_color(count, self._max).name()
            else:
                fill = NO_DATA_FILL
            d = []
            for ring in rings:
                x0, y0 = _project(ring[0][0], ring[0][1])
                d.append(f"M{x0:.1f} {y0:.1f}")
                for lon, lat in ring[1:]:
                    x, y = _project(lon, lat)
                    d.append(f"L{x:.1f} {y:.1f}")
                d.append("Z")
            parts.append(f'<path d="{"".join(d)}" fill="{fill}" stroke="{BORDER}" stroke-width="0.6"/>')
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
        self._renderer = QtSvg.QSvgRenderer(bytearray(self._svg(), "utf-8"))
        self._item = QtSvgWidgets.QGraphicsSvgItem()
        self._item.setSharedRenderer(self._renderer)
        self._scene.addItem(self._item)
        self._scene.setSceneRect(0, 0, _W, _H)
        self._did_fit = False
        self._fit()

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
