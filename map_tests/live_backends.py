"""Registry + dispatcher for the Map Lab's live map widgets.

Each renderer lives in live_<name>.py exposing build_widget(parent) -> QWidget.
"""
from __future__ import annotations
import importlib
from typing import Optional

# name, label, tag ("safe" = no OpenGL, cannot crash | "gl" = needs GPU GL)
CATALOG = [
    ("tilemap", "Tile map — real basemap (native)", "safe", "OSM/Carto tiles + country choropleth, drag-pan + wheel-zoom. Web-map look, NO WebEngine, NO GL. Needs internet for tiles."),
    ("raster", "Raster (Qt) — in the app now", "safe", "QGraphicsView choropleth. No GL, cannot crash. This is what ships."),
    ("svg", "SVG vector", "safe", "QtSvg vector map. No GL."),
    ("matplotlib", "matplotlib", "safe", "Agg canvas + pan/zoom toolbar. No GL."),
    ("cartopy", "Cartopy", "safe", "Projected map + toolbar, offline. No GL."),
    ("plotly", "plotly (WebEngine)", "gl", "Interactive plotly in Chromium — needs GPU OpenGL."),
    ("webengine", "Inline SVG (WebEngine)", "gl", "Self-contained SVG map in Chromium — needs GPU OpenGL."),
    ("maplibre", "MapLibre GL — real web map", "gl", "Vector web map: Carto basemap + choropleth + points, smooth pan/zoom. WebGL — needs GPU (post-reboot) + internet."),
]

NAMES = [c[0] for c in CATALOG]


def meta(name: str) -> Optional[dict]:
    for n, label, tag, desc in CATALOG:
        if n == name:
            return {"name": n, "label": label, "tag": tag, "desc": desc}
    return None


def build(name: str, parent=None):
    """Import live_<name> and return its live widget."""
    module = importlib.import_module(f"live_{name}")
    return module.build_widget(parent)
