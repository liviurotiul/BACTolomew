"""Shared geo helpers: world polygons, equirectangular projection, colour ramp,
and matching fixture country names to the Natural Earth geojson."""
from __future__ import annotations
import json
import re
from pathlib import Path
from typing import Dict, List, Tuple

_ASSETS = Path(__file__).resolve().parents[1] / "src/annotation_explorer_polars/assets"
GEOJSON_PATH = _ASSETS / "world_countries.geojson"

# theme
NO_DATA_FILL = "#eef1f5"
BORDER = "#c2ccd6"
OCEAN = "#f4f8fc"
POINT_FILL = "#b11242"
POINT_EDGE = "#ffffff"
INK = "#1f2d3d"
_RAMP_LO = (222, 235, 247)   # light blue
_RAMP_HI = (8, 69, 148)      # dark blue

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


def load_world() -> List[Tuple[str, List[List[Tuple[float, float]]]]]:
    """Return [(country_name, [ring_of_(lon,lat), ...]), ...]."""
    data = json.loads(GEOJSON_PATH.read_text())
    out: List[Tuple[str, List[List[Tuple[float, float]]]]] = []
    for feat in data["features"]:
        props = feat.get("properties", {})
        name = props.get("NAME") or props.get("ADMIN") or props.get("SOVEREIGNT") or ""
        geom = feat.get("geometry") or {}
        gtype = geom.get("type")
        coords = geom.get("coordinates") or []
        rings: List[List[Tuple[float, float]]] = []
        if gtype == "Polygon":
            polys = [coords]
        elif gtype == "MultiPolygon":
            polys = coords
        else:
            polys = []
        for poly in polys:
            for ring in poly:
                rings.append([(float(x), float(y)) for x, y in ring])
        if rings:
            out.append((name, rings))
    return out


def counts_by_geojson_name(countries: Dict[str, int]) -> Dict[str, int]:
    """Map fixture country counts onto normalized names used by the geojson."""
    out: Dict[str, int] = {}
    for name, count in countries.items():
        out[_norm(name)] = out.get(_norm(name), 0) + int(count)
    return out


def match_count(geojson_name: str, counts_norm: Dict[str, int]) -> int:
    return int(counts_norm.get(_norm(geojson_name), 0))


def project(lon: float, lat: float, w: float, h: float) -> Tuple[float, float]:
    """Equirectangular (plate carrée) lon/lat -> pixel x/y (y down)."""
    x = (float(lon) + 180.0) / 360.0 * w
    y = (90.0 - float(lat)) / 180.0 * h
    return x, y


def count_color(count: int, max_count: int) -> str:
    if count <= 0 or max_count <= 0:
        return NO_DATA_FILL
    import math
    # log scaling so small counts stay visible against big ones
    t = math.log1p(count) / math.log1p(max_count)
    t = max(0.0, min(1.0, t))
    r = round(_RAMP_LO[0] + (_RAMP_HI[0] - _RAMP_LO[0]) * t)
    g = round(_RAMP_LO[1] + (_RAMP_HI[1] - _RAMP_LO[1]) * t)
    b = round(_RAMP_LO[2] + (_RAMP_HI[2] - _RAMP_LO[2]) * t)
    return f"#{r:02x}{g:02x}{b:02x}"


def text_on(count: int, max_count: int) -> str:
    """Readable label colour for a choropleth cell."""
    import math
    t = 0.0 if max_count <= 0 else math.log1p(count) / math.log1p(max_count)
    return "#ffffff" if t >= 0.55 else INK
