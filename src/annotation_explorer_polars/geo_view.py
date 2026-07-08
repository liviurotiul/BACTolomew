from __future__ import annotations

import json
import time
import math
import re
import unicodedata
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from PyQt6 import QtCore, QtGui, QtWidgets
from .constants import COUNTRY_ALIASES
try:
    from PyQt6.QtOpenGLWidgets import QOpenGLWidget
except Exception:
    QOpenGLWidget = None


_DEFAULT_PEN = QtGui.QPen(QtGui.QColor("#000000"))
_DEFAULT_PEN.setWidthF(0.25)
_LINE_WIDTH_STEPS = [
    (0.0, 1.0, 0.42),
    (1.0, 1.8, 0.32),
    (1.8, 3.0, 0.2),
    (3.0, float("inf"), 0.12),
]
_SIMPLIFY_TOLERANCES = {"low": 0.5, "high": 0.2, "ultra": 0.05}


def _normalize_country_name(name: str) -> str:
    if not name:
        return ""
    folded = unicodedata.normalize("NFKD", str(name))
    folded = "".join(ch for ch in folded if not unicodedata.combining(ch))
    key = re.sub(r"[^a-z0-9]+", " ", folded.lower()).strip()
    alias = COUNTRY_ALIASES.get(key)
    if alias:
        folded_alias = unicodedata.normalize("NFKD", str(alias))
        folded_alias = "".join(ch for ch in folded_alias if not unicodedata.combining(ch))
        key = re.sub(r"[^a-z0-9]+", " ", folded_alias.lower()).strip()
    return key


def _geo_cache_dir() -> Path:
    from .cache import ensure_cache_dir

    root = ensure_cache_dir()
    geo_dir = root / "geo"
    geo_dir.mkdir(parents=True, exist_ok=True)
    return geo_dir


def _simplified_geojson_path(src_path: Path, label: str, tolerance: float) -> Path:
    tol_tag = f"{tolerance:.3f}".replace(".", "p")
    return _geo_cache_dir() / f"{src_path.stem}_{label}_simp_{tol_tag}.geojson"


def _dist_to_segment_sq(
    px: float, py: float, ax: float, ay: float, bx: float, by: float
) -> float:
    dx = bx - ax
    dy = by - ay
    if dx == 0 and dy == 0:
        return (px - ax) ** 2 + (py - ay) ** 2
    t = ((px - ax) * dx + (py - ay) * dy) / (dx * dx + dy * dy)
    if t <= 0.0:
        cx, cy = ax, ay
    elif t >= 1.0:
        cx, cy = bx, by
    else:
        cx = ax + t * dx
        cy = ay + t * dy
    return (px - cx) ** 2 + (py - cy) ** 2


def _douglas_peucker(points: List[List[float]], tolerance: float) -> List[List[float]]:
    if len(points) <= 2:
        return points
    tol_sq = tolerance * tolerance
    keep = [False] * len(points)
    keep[0] = True
    keep[-1] = True
    stack = [(0, len(points) - 1)]
    while stack:
        start, end = stack.pop()
        ax, ay = points[start]
        bx, by = points[end]
        max_dist = -1.0
        index = None
        for i in range(start + 1, end):
            px, py = points[i]
            dist = _dist_to_segment_sq(px, py, ax, ay, bx, by)
            if dist > max_dist:
                max_dist = dist
                index = i
        if max_dist > tol_sq and index is not None:
            keep[index] = True
            stack.append((start, index))
            stack.append((index, end))
    return [pt for i, pt in enumerate(points) if keep[i]]


def _simplify_ring(ring: List[List[float]], tolerance: float) -> List[List[float]]:
    if not ring or len(ring) <= 4:
        return ring
    closed = ring[0] == ring[-1]
    pts = ring[:-1] if closed else ring
    simplified = _douglas_peucker(pts, tolerance)
    if len(simplified) < 3:
        return ring
    if closed and simplified[0] != simplified[-1]:
        simplified.append(simplified[0])
    return simplified


def _simplify_geojson_data(data: dict, tolerance: float) -> dict:
    features = data.get("features") if isinstance(data, dict) else None
    if not features:
        return data
    for feature in features:
        if not isinstance(feature, dict):
            continue
        geom = feature.get("geometry")
        if not isinstance(geom, dict):
            continue
        gtype = geom.get("type")
        coords = geom.get("coordinates")
        if not coords or not gtype:
            continue
        if gtype == "Polygon":
            new_polys = []
            for ring in coords:
                if isinstance(ring, list):
                    new_polys.append(_simplify_ring(ring, tolerance))
            geom["coordinates"] = new_polys
        elif gtype == "MultiPolygon":
            new_multi = []
            for poly in coords:
                if not isinstance(poly, list):
                    continue
                new_poly = []
                for ring in poly:
                    if isinstance(ring, list):
                        new_poly.append(_simplify_ring(ring, tolerance))
                if new_poly:
                    new_multi.append(new_poly)
            geom["coordinates"] = new_multi
    return data


def _ensure_simplified_geojson(src_path: Path, label: str, tolerance: float) -> Path:
    if tolerance <= 0 or "minimal" in src_path.name.lower():
        return src_path
    out_path = _simplified_geojson_path(src_path, label, tolerance)
    try:
        if out_path.exists() and out_path.stat().st_mtime >= src_path.stat().st_mtime:
            return out_path
    except OSError:
        pass
    try:
        t0 = time.perf_counter()
        print(f"[geo] simplifying {src_path} (tol={tolerance})")
        with src_path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
        data = _simplify_geojson_data(data, tolerance)
        out_path.write_text(json.dumps(data))
        print(f"[geo] wrote {out_path} in {time.perf_counter() - t0:.2f}s")
        return out_path
    except Exception as exc:
        print(f"[geo] simplification failed for {src_path}: {exc}")
        return src_path


def _best_name(props: dict) -> str:
    for key in ("name", "NAME", "ADMIN", "name_long", "NAME_LONG"):
        value = props.get(key)
        if value:
            return str(value)
    return "Unknown"


def _count_color(count: int, max_count: int) -> QtGui.QColor:
    if count <= 0:
        return QtGui.QColor("#eeeeee")
    denom = max(max_count, 1)
    t = count / denom
    if t > 0.8:
        return QtGui.QColor("#bd0026")
    if t > 0.6:
        return QtGui.QColor("#f03b20")
    if t > 0.4:
        return QtGui.QColor("#fd8d3c")
    if t > 0.2:
        return QtGui.QColor("#feb24c")
    return QtGui.QColor("#fed976")


def _lonlat_to_scene(lon: float, lat: float) -> QtCore.QPointF:
    x = lon + 180.0
    y = 90.0 - lat
    return QtCore.QPointF(x, y)


@dataclass
class CountryShape:
    name_key: str
    display_name: str
    item: QtWidgets.QGraphicsPathItem
    bounds: QtCore.QRectF


class GeoMapView(QtWidgets.QGraphicsView):
    def __init__(
        self,
        geojson_path: Path,
        high_res_path: Optional[Path] = None,
        ultra_res_path: Optional[Path] = None,
        parent: Optional[QtWidgets.QWidget] = None,
        use_opengl: bool = True,
    ) -> None:
        super().__init__(parent)
        self._scene = QtWidgets.QGraphicsScene(self)
        self._scene.setItemIndexMethod(QtWidgets.QGraphicsScene.ItemIndexMethod.BspTreeIndex)
        self.setScene(self._scene)
        # An OpenGL viewport is smoother, but on machines where GL context
        # creation fails it hard-crashes the whole app. Callers on such setups
        # pass use_opengl=False to fall back to the raster viewport.
        if use_opengl and QOpenGLWidget is not None:
            self.setViewport(QOpenGLWidget())
        self.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, False)
        self.setDragMode(QtWidgets.QGraphicsView.DragMode.ScrollHandDrag)
        self.setViewportUpdateMode(QtWidgets.QGraphicsView.ViewportUpdateMode.BoundingRectViewportUpdate)
        self.setTransformationAnchor(QtWidgets.QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setResizeAnchor(QtWidgets.QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setOptimizationFlag(QtWidgets.QGraphicsView.OptimizationFlag.DontSavePainterState, True)
        self.setOptimizationFlag(QtWidgets.QGraphicsView.OptimizationFlag.DontAdjustForAntialiasing, True)
        self.setCacheMode(QtWidgets.QGraphicsView.CacheModeFlag.CacheBackground)

        self._country_shapes: Dict[str, List[CountryShape]] = {}
        self._all_country_shapes: List[CountryShape] = []
        self._layers: Dict[str, Tuple[QtWidgets.QGraphicsItemGroup, Dict[str, List[CountryShape]], List[CountryShape]]] = {}
        self._current_lod = "low"
        self._current_group: Optional[QtWidgets.QGraphicsItemGroup] = None
        self._points_group: Optional[QtWidgets.QGraphicsItemGroup] = None
        self._point_items_by_sample: Dict[str, List[QtWidgets.QGraphicsEllipseItem]] = {}
        self._highlight_sample: Optional[str] = None
        self._highlight_samples: set[str] = set()
        self._did_initial_fit = False
        self._country_counts: Dict[str, int] = {}
        self._max_country_count: int = 1
        self._lod_norm_threshold_hi_in = 1.0
        self._lod_norm_threshold_hi_out = 0.9
        self._lod_norm_threshold_ultra_in = 2.0
        self._lod_norm_threshold_ultra_out = 1.8
        self._debug_lod_text = ""
        self._lod_timer = QtCore.QTimer(self)
        self._lod_timer.setSingleShot(True)
        self._lod_timer.setInterval(120)
        self._lod_timer.timeout.connect(self._update_lod)
        self._zoom_anim_timer = QtCore.QTimer(self)
        self._zoom_anim_timer.setInterval(16)
        self._zoom_anim_timer.timeout.connect(self._on_zoom_anim)
        self._zoom_anim_active = False
        self._zoom_anim_start_scale = 1.0
        self._zoom_anim_target_scale = 1.0
        self._zoom_anim_start_time = 0.0
        self._zoom_anim_duration_ms = 140

        low_tol = _SIMPLIFY_TOLERANCES.get("low", 0.0)
        high_tol = _SIMPLIFY_TOLERANCES.get("high", 0.0)
        self._lod_paths: Dict[str, Path] = {
            "low": _ensure_simplified_geojson(geojson_path, "low", low_tol)
        }
        if high_res_path and high_res_path.exists():
            self._lod_paths["high"] = _ensure_simplified_geojson(high_res_path, "high", high_tol)
        if ultra_res_path and ultra_res_path.exists():
            ultra_tol = _SIMPLIFY_TOLERANCES.get("ultra", 0.05)
            self._lod_paths["ultra"] = _ensure_simplified_geojson(ultra_res_path, "ultra", ultra_tol)

        self._load_layer("low")
        self._activate_layer("low")
        self._scene.setSceneRect(0.0, 0.0, 360.0, 180.0)
        self._fit_world()

    def highlight_sample(self, sample: Optional[str]) -> None:
        value = str(sample).strip() if sample is not None else ""
        self._highlight_sample = value or None
        if self._highlight_sample:
            self._highlight_samples = set()
        self._apply_points_highlight()

    def highlight_samples(self, samples: Optional[Sequence[str]]) -> None:
        if not samples:
            self._highlight_samples = set()
            self._apply_points_highlight()
            return
        self._highlight_samples = {str(s).strip() for s in samples if str(s).strip()}
        if self._highlight_samples:
            self._highlight_sample = None
        self._apply_points_highlight()

    def _apply_points_highlight(self) -> None:
        if not self._point_items_by_sample:
            return
        wanted = set(self._highlight_samples)
        if self._highlight_sample:
            wanted = {self._highlight_sample}
        has_wanted = bool(wanted)
        for sample, items in self._point_items_by_sample.items():
            selected = sample in wanted if has_wanted else False
            if selected:
                pen = QtGui.QPen(QtGui.QColor("#c2410c"))
                brush = QtGui.QBrush(QtGui.QColor("#f97316"))
                z = 12
            else:
                pen = QtGui.QPen(QtGui.QColor("#1f77b4"))
                base = QtGui.QColor("#1f77b4")
                if has_wanted:
                    base.setAlpha(70)
                brush = QtGui.QBrush(base)
                z = 10
            for item in items:
                item.setPen(pen)
                item.setBrush(brush)
                item.setZValue(z)

    def _fit_world(self) -> None:
        self.fitInView(self._scene.sceneRect(), QtCore.Qt.AspectRatioMode.KeepAspectRatio)
        self._update_lod()

    def _load_layer(self, lod: str) -> None:
        if lod in self._layers:
            return
        path = self._lod_paths.get(lod)
        if not path:
            return
        group, shapes_map, shapes_list = self._load_geojson(path)
        self._layers[lod] = (group, shapes_map, shapes_list)

    def _activate_layer(self, lod: str) -> None:
        if lod not in self._layers:
            return
        if self._current_group is not None:
            self._current_group.setVisible(False)
        group, shapes_map, shapes_list = self._layers[lod]
        group.setVisible(True)
        self._current_group = group
        self._country_shapes = shapes_map
        self._all_country_shapes = shapes_list
        self._current_lod = lod
        self._apply_pen_width()
        self._apply_country_counts()

    def _load_geojson(
        self, path: Path
    ) -> Tuple[QtWidgets.QGraphicsItemGroup, Dict[str, List[CountryShape]], List[CountryShape]]:
        if not path.exists():
            raise FileNotFoundError(str(path))
        t0 = time.perf_counter()
        print(f"[file] reading {path}")
        with path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
        print(f"[file] read {path} in {time.perf_counter() - t0:.2f}s")

        features = data.get("features") if isinstance(data, dict) else None
        if not features:
            raise ValueError("GeoJSON does not contain any features.")

        group = QtWidgets.QGraphicsItemGroup()
        shapes_map: Dict[str, List[CountryShape]] = {}
        shapes_list: List[CountryShape] = []

        for feature in features:
            if not isinstance(feature, dict):
                continue
            geom = feature.get("geometry") or {}
            props = feature.get("properties") or {}
            name = _best_name(props)
            name_key = _normalize_country_name(name)
            if not name_key:
                continue

            path_item = self._geometry_to_path(geom)
            if path_item is None:
                continue
            item = QtWidgets.QGraphicsPathItem(path_item)
            item.setCacheMode(QtWidgets.QGraphicsItem.CacheMode.DeviceCoordinateCache)
            item.setPen(_DEFAULT_PEN)
            item.setBrush(QtGui.QBrush(QtGui.QColor("#eeeeee")))
            item.setToolTip(name)
            item.setZValue(0)
            group.addToGroup(item)
            bounds = item.boundingRect()
            shape = CountryShape(name_key=name_key, display_name=name, item=item, bounds=bounds)
            shapes_map.setdefault(name_key, []).append(shape)
            shapes_list.append(shape)

        self._scene.addItem(group)
        group.setVisible(False)
        return group, shapes_map, shapes_list

    def _geometry_to_path(self, geom: dict) -> Optional[QtGui.QPainterPath]:
        gtype = geom.get("type")
        coords = geom.get("coordinates")
        if not gtype or not coords:
            return None

        polygons: Iterable = []
        if gtype == "Polygon":
            polygons = [coords]
        elif gtype == "MultiPolygon":
            polygons = coords
        else:
            return None

        path = QtGui.QPainterPath()
        path.setFillRule(QtCore.Qt.FillRule.OddEvenFill)

        for polygon in polygons:
            if not polygon:
                continue
            for ring in polygon:
                if not ring:
                    continue
                ring_points: List[QtCore.QPointF] = []
                for coord in ring:
                    if coord is None or len(coord) < 2:
                        continue
                    lon = float(coord[0])
                    lat = float(coord[1])
                    ring_points.append(_lonlat_to_scene(lon, lat))
                if len(ring_points) < 3:
                    continue
                path.moveTo(ring_points[0])
                for pt in ring_points[1:]:
                    path.lineTo(pt)
                path.closeSubpath()

        return path if not path.isEmpty() else None

    def update_data(self, points: List[dict], countries: Dict[str, int], max_count: int) -> None:
        normalized_counts: Dict[str, int] = {}
        for raw, count in countries.items():
            key = _normalize_country_name(str(raw))
            if not key:
                continue
            normalized_counts[key] = normalized_counts.get(key, 0) + int(count)

        self._country_counts = normalized_counts
        self._max_country_count = max_count
        self._apply_country_counts()

        self._render_points(points)
        self._auto_fit(points, normalized_counts)

    def _apply_country_counts(self) -> None:
        for shape in self._all_country_shapes:
            count = self._country_counts.get(shape.name_key, 0)
            shape.item.setBrush(QtGui.QBrush(_count_color(count, self._max_country_count)))
            if count:
                shape.item.setToolTip(f"{shape.display_name}: {count}")
            else:
                shape.item.setToolTip(shape.display_name)

    def _apply_pen_width(self) -> None:
        width = self._current_line_width()
        pen = QtGui.QPen(QtGui.QColor("#000000"))
        pen.setWidthF(width)
        for shape in self._all_country_shapes:
            shape.item.setPen(pen)
        if self._points_group is not None:
            self._points_group.setZValue(10)

    def _current_line_width(self) -> float:
        ratio = self._normalized_scale()
        for _, hi, width in _LINE_WIDTH_STEPS:
            if ratio <= hi:
                return width
        return _LINE_WIDTH_STEPS[-1][2]

    def _fit_scale_for_viewport(self) -> float:
        scene_rect = self.sceneRect()
        vp = self.viewport().rect()
        if scene_rect.width() <= 0 or scene_rect.height() <= 0 or vp.width() <= 0 or vp.height() <= 0:
            return 1.0
        sx = float(vp.width()) / float(scene_rect.width())
        sy = float(vp.height()) / float(scene_rect.height())
        fit_scale = min(sx, sy)
        return max(fit_scale, 1e-9)

    def _normalized_scale(self) -> float:
        current = self.transform().m11()
        if current <= 0:
            return 1.0
        return float(current) / self._fit_scale_for_viewport()

    def _update_debug_overlay(self) -> None:
        norm = self._normalized_scale()
        if "ultra" in self._lod_paths:
            if norm >= self._lod_norm_threshold_ultra_in:
                zone = "ultra in"
            elif norm >= self._lod_norm_threshold_hi_in:
                zone = "high in"
            elif norm <= self._lod_norm_threshold_hi_out:
                zone = "low in"
            else:
                zone = "hold"
            self._debug_lod_text = (
                f"LOD={self._current_lod} | norm={norm:.3f} | "
                f"high:{self._lod_norm_threshold_hi_in:.2f}/{self._lod_norm_threshold_hi_out:.2f} "
                f"ultra:{self._lod_norm_threshold_ultra_in:.2f}/{self._lod_norm_threshold_ultra_out:.2f} "
                f"| {zone}"
            )
        else:
            if norm >= self._lod_norm_threshold_hi_in:
                zone = "above hi"
            elif norm <= self._lod_norm_threshold_hi_out:
                zone = "below lo"
            else:
                zone = "between"
            self._debug_lod_text = (
                f"LOD={self._current_lod} | norm={norm:.3f} | "
                f"hi:{self._lod_norm_threshold_hi_in:.2f}/{self._lod_norm_threshold_hi_out:.2f} | {zone}"
            )
        self.viewport().update()

    def _render_points(self, points: List[dict]) -> None:
        if self._points_group is not None:
            self._scene.removeItem(self._points_group)
            self._points_group = None
        self._point_items_by_sample = {}

        if not points:
            return

        group = QtWidgets.QGraphicsItemGroup()
        for p in points:
            lat = p.get("lat")
            lon = p.get("lon")
            if lat is None or lon is None or not math.isfinite(lat) or not math.isfinite(lon):
                continue
            pt = _lonlat_to_scene(float(lon), float(lat))
            radius = 3.8
            item = QtWidgets.QGraphicsEllipseItem(
                pt.x() - radius, pt.y() - radius, radius * 2, radius * 2
            )
            item.setCacheMode(QtWidgets.QGraphicsItem.CacheMode.DeviceCoordinateCache)
            item.setPen(QtGui.QPen(QtGui.QColor("#1f77b4")))
            item.setBrush(QtGui.QBrush(QtGui.QColor("#1f77b4")))
            item.setZValue(10)
            label = p.get("sample") or "Sample"
            count = p.get("count") or 0
            item.setToolTip(f"{label} ({count})")
            group.addToGroup(item)
            sample_key = str(label).strip()
            if sample_key:
                self._point_items_by_sample.setdefault(sample_key, []).append(item)
        self._scene.addItem(group)
        self._points_group = group
        self._apply_points_highlight()

    def _auto_fit(self, points: List[dict], countries: Dict[str, int]) -> None:
        bounds = self._points_bounds(points)
        if bounds is None and countries:
            bounds = self._countries_bounds(countries.keys())
        if bounds is None:
            self._fit_world()
            return

        padded = QtCore.QRectF(bounds)
        pad = 8.0
        padded.adjust(-pad, -pad, pad, pad)
        self.fitInView(padded, QtCore.Qt.AspectRatioMode.KeepAspectRatio)
        self._update_lod()

    def _points_bounds(self, points: List[dict]) -> Optional[QtCore.QRectF]:
        xs: List[float] = []
        ys: List[float] = []
        for p in points:
            lat = p.get("lat")
            lon = p.get("lon")
            if lat is None or lon is None or not math.isfinite(lat) or not math.isfinite(lon):
                continue
            pt = _lonlat_to_scene(float(lon), float(lat))
            xs.append(pt.x())
            ys.append(pt.y())
        if not xs:
            return None
        return QtCore.QRectF(min(xs), min(ys), max(xs) - min(xs), max(ys) - min(ys))

    def _countries_bounds(self, country_keys: Iterable[str]) -> Optional[QtCore.QRectF]:
        rect: Optional[QtCore.QRectF] = None
        for key in country_keys:
            shapes = self._country_shapes.get(key)
            if not shapes:
                continue
            for shape in shapes:
                rect = shape.bounds if rect is None else rect.united(shape.bounds)
        return rect

    def wheelEvent(self, event: QtGui.QWheelEvent) -> None:
        delta = event.angleDelta().y()
        if delta == 0:
            delta = event.pixelDelta().y()
        if delta == 0:
            return
        # Smooth continuous zoom with short animation between steps.
        factor = math.pow(1.0015, delta)
        base_scale = (
            self._zoom_anim_target_scale if self._zoom_anim_active else self.transform().m11()
        )
        target_scale = base_scale * factor
        self._start_zoom_animation(target_scale)
        event.accept()

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        super().resizeEvent(event)
        if not self._did_initial_fit:
            self._fit_world()
            self._did_initial_fit = True
        self._update_debug_overlay()
        self._schedule_lod_update()

    def _schedule_lod_update(self) -> None:
        if "high" not in self._lod_paths and "ultra" not in self._lod_paths:
            self._update_debug_overlay()
            return
        if self._lod_timer.isActive():
            self._lod_timer.stop()
        self._lod_timer.start()

    def _update_lod(self) -> None:
        has_high = "high" in self._lod_paths
        has_ultra = "ultra" in self._lod_paths
        if not has_high and not has_ultra:
            self._update_debug_overlay()
            return
        norm = self._normalized_scale()
        if self._current_lod == "low":
            if has_high and norm > self._lod_norm_threshold_hi_in:
                self._load_layer("high")
                self._activate_layer("high")
            elif has_ultra and not has_high and norm > self._lod_norm_threshold_ultra_in:
                self._load_layer("ultra")
                self._activate_layer("ultra")
        elif self._current_lod == "high":
            if has_ultra and norm > self._lod_norm_threshold_ultra_in:
                self._load_layer("ultra")
                self._activate_layer("ultra")
            elif norm < self._lod_norm_threshold_hi_out:
                self._activate_layer("low")
        elif self._current_lod == "ultra":
            if has_high and norm < self._lod_norm_threshold_ultra_out:
                self._activate_layer("high")
            elif (not has_high) and norm < self._lod_norm_threshold_ultra_out:
                self._activate_layer("low")
        self._apply_pen_width()
        self._update_debug_overlay()

    def _start_zoom_animation(self, target_scale: float) -> None:
        current_scale = self.transform().m11()
        if current_scale <= 0:
            return
        self._zoom_anim_start_scale = current_scale
        self._zoom_anim_target_scale = max(0.001, target_scale)
        self._zoom_anim_start_time = time.perf_counter()
        self._zoom_anim_active = True
        if not self._zoom_anim_timer.isActive():
            self._zoom_anim_timer.start()

    def _on_zoom_anim(self) -> None:
        if not self._zoom_anim_active:
            return
        elapsed = (time.perf_counter() - self._zoom_anim_start_time) * 1000.0
        t = min(1.0, elapsed / float(self._zoom_anim_duration_ms))
        # Smoothstep easing for gentle acceleration/deceleration.
        eased = t * t * (3.0 - 2.0 * t)
        desired = self._zoom_anim_start_scale + (
            (self._zoom_anim_target_scale - self._zoom_anim_start_scale) * eased
        )
        current = self.transform().m11()
        if current > 0:
            ratio = desired / current
            if ratio > 0:
                self.scale(ratio, ratio)
        self._update_debug_overlay()
        if t >= 1.0:
            self._zoom_anim_active = False
            self._zoom_anim_timer.stop()
            self._schedule_lod_update()

    def drawForeground(self, painter: QtGui.QPainter, rect: QtCore.QRectF) -> None:  # type: ignore[override]
        super().drawForeground(painter, rect)
        if not self._debug_lod_text:
            return
        painter.save()
        painter.resetTransform()
        font = painter.font()
        font.setPointSize(max(9, font.pointSize()))
        painter.setFont(font)
        metrics = QtGui.QFontMetrics(font)
        text_rect = metrics.boundingRect(self._debug_lod_text)
        margin = 8
        pad_x = 8
        pad_y = 5
        box = QtCore.QRectF(
            float(margin),
            float(margin),
            float(text_rect.width() + pad_x * 2),
            float(text_rect.height() + pad_y * 2),
        )
        painter.setPen(QtCore.Qt.PenStyle.NoPen)
        painter.setBrush(QtGui.QColor(0, 0, 0, 170))
        painter.drawRoundedRect(box, 6.0, 6.0)
        painter.setPen(QtGui.QColor("#ffffff"))
        painter.drawText(
            box.adjusted(float(pad_x), float(pad_y), -float(pad_x), -float(pad_y)),
            int(QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignVCenter),
            self._debug_lod_text,
        )
        painter.restore()
