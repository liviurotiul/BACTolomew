"""LIVE native slippy tile map for the Map Lab — real OSM/Carto tiles + points,
pan/zoom, NO WebEngine and NO OpenGL. Uses the shipping-candidate widget."""
from __future__ import annotations
import sys

_SRC = "/home/bals/Documents/dev/BACTolomew/src"
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

import fixture  # noqa: E402
from annotation_explorer_polars.geo_tile_map import TileMapView  # noqa: E402


def build_widget(parent=None):
    view = TileMapView(parent)
    points = [
        {"sample": p["sample"], "lat": p["lat"], "lon": p["lon"], "country": p["country"]}
        for p in fixture.POINTS
    ]
    view.update_data(points=points, countries=fixture.COUNTRIES, max_count=fixture.MAX_COUNT)
    return view
