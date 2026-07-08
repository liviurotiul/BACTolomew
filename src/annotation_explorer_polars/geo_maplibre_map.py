"""MapLibre GL choropleth map — the primary Geo-panel map when a GPU/WebEngine is
available: a real vector web map (Carto basemap + country choropleth coloured by
sample count, NO points). Runs in a QWebEngineView.

Emits renderUnavailable(str) if WebGL/network/page init fails, so the panel can
swap to the offline SVG map. (Creating the WebView itself is gated by the panel's
WebEngine probe, so we never realize it on a machine with broken GL.)
"""
from __future__ import annotations

import json
import tempfile
from typing import Dict, List, Optional, Sequence

from PyQt6 import QtCore, QtWebEngineCore, QtWidgets
from PyQt6.QtWebEngineWidgets import QWebEngineView

from .geo_raster_map import _ASSETS, _count_color, _norm

TILE_URL = "https://a.basemaps.cartocdn.com/light_all/{z}/{x}/{y}.png"
MAPLIBRE_JS = "https://unpkg.com/maplibre-gl@4.7.1/dist/maplibre-gl.js"
MAPLIBRE_CSS = "https://unpkg.com/maplibre-gl@4.7.1/dist/maplibre-gl.css"
HILITE = "#f59e0b"


class MapLibreGeoView(QWebEngineView):
    renderUnavailable = QtCore.pyqtSignal(str)

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        settings = self.settings()
        for attr in ("LocalContentCanAccessRemoteUrls", "LocalContentCanAccessFileUrls"):
            settings.setAttribute(getattr(QtWebEngineCore.QWebEngineSettings.WebAttribute, attr), True)
        self._tmp: List[str] = []
        self._counts: Dict[str, int] = {}
        self._sample_country: Dict[str, str] = {}
        self._max = 1
        self._points: List[dict] = []
        self._ready = False
        self._pending_hot: List[str] = []
        self.loadFinished.connect(self._on_load)
        self.titleChanged.connect(self._on_title)

    # ---- API ----
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
        self._counts = counts
        self._max = max([int(max_count or 1), *counts.values()] or [1])
        self._reload()

    def highlight_sample(self, sample: Optional[str]) -> None:
        self.highlight_samples([sample] if sample else [])

    def highlight_samples(self, samples: Optional[Sequence[str]]) -> None:
        hot = sorted({
            self._sample_country[str(s).strip()]
            for s in (samples or []) if str(s).strip() in self._sample_country
        })
        if self._ready:
            self.page().runJavaScript(f"window.__setHot && window.__setHot({json.dumps(hot)});")
        else:
            self._pending_hot = hot

    # ---- html ----
    def _choropleth_geojson(self) -> dict:
        src = next((_ASSETS / n for n in
                    ("world_countries_10m.geojson", "world_countries_50m.geojson", "world_countries.geojson")
                    if (_ASSETS / n).exists()), None)
        feats = []
        if src is not None:
            data = json.loads(src.read_text())
            for feat in data.get("features", []):
                props = feat.get("properties", {})
                name = props.get("NAME") or props.get("ADMIN") or props.get("SOVEREIGNT") or ""
                key = _norm(name)
                count = self._counts.get(key, 0)
                if count <= 0:
                    continue
                feats.append({
                    "type": "Feature",
                    "properties": {"id": key, "name": name, "count": count,
                                   "color": _count_color(count, self._max).name()},
                    "geometry": feat.get("geometry"),
                })
        return {"type": "FeatureCollection", "features": feats}

    def _bounds(self) -> Optional[list]:
        if not self._points:
            return None
        lons = [p["lon"] for p in self._points]
        lats = [p["lat"] for p in self._points]
        return [[min(lons), min(lats)], [max(lons), max(lats)]]

    def _legend_html(self) -> str:
        stops = sorted({1, max(2, self._max // 8), max(4, self._max // 4), max(8, self._max // 2), self._max})
        cells = "".join(
            f'<i style="background:{_count_color(s, self._max).name()}"></i>{s} '
            for s in stops
        )
        return f'<b>samples / country</b><br/>{cells}'

    def _html(self) -> str:
        countries = json.dumps(self._choropleth_geojson())
        bounds = json.dumps(self._bounds())
        legend = self._legend_html()
        return _TEMPLATE.format(
            css=MAPLIBRE_CSS, js=MAPLIBRE_JS, tile=TILE_URL, hilite=HILITE,
            countries=countries, bounds=bounds, legend=legend,
        )

    def _reload(self) -> None:
        self._ready = False
        tmp = tempfile.NamedTemporaryFile("w", suffix=".html", delete=False, encoding="utf-8")
        tmp.write(self._html())
        tmp.flush()
        tmp.close()
        self._tmp.append(tmp.name)
        self.load(QtCore.QUrl.fromLocalFile(tmp.name))

    def _on_load(self, ok: bool) -> None:
        if not ok:
            self.renderUnavailable.emit("page failed to load")

    def _on_title(self, title: str) -> None:
        if title == "MAPFAIL":
            self.renderUnavailable.emit("WebGL/map initialisation failed")
        elif title == "MAPOK":
            self._ready = True
            if self._pending_hot:
                self.highlight_samples([])  # cleared below; reapply via stored countries
                self.page().runJavaScript(f"window.__setHot && window.__setHot({json.dumps(self._pending_hot)});")
                self._pending_hot = []


def _has_coords(p: dict) -> bool:
    try:
        return not (p["lat"] is None or p["lon"] is None) and float(p["lat"]) == float(p["lat"])
    except (KeyError, TypeError, ValueError):
        return False


_TEMPLATE = """<!DOCTYPE html><html><head><meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<link href="{css}" rel="stylesheet"/><script src="{js}"></script>
<style>html,body,#map{{margin:0;height:100%;width:100%}}
.legend{{position:absolute;left:10px;bottom:22px;background:rgba(255,255,255,.92);padding:8px 11px;border-radius:6px;
 font:12px sans-serif;color:#1f2d3d;box-shadow:0 1px 5px rgba(0,0,0,.2)}}
.legend i{{display:inline-block;width:22px;height:11px;margin-right:2px;vertical-align:middle}}</style></head>
<body><div id="map"></div><div class="legend">{legend}</div>
<script>
var COUNTRIES={countries}, BOUNDS={bounds};
function fail(){{document.title='MAPFAIL';}}
setTimeout(function(){{if(document.title!=='MAPOK')fail();}},9000);
try{{
 if(typeof maplibregl==='undefined'){{fail();}}   // CDN failed to load
 else{{
  var style={{version:8,sources:{{
    carto:{{type:'raster',tiles:['{tile}'],tileSize:256,attribution:'© OpenStreetMap © CARTO'}},
    countries:{{type:'geojson',data:COUNTRIES,promoteId:'id'}}}},
   layers:[{{id:'bg',type:'raster',source:'carto'}},
    {{id:'choro',type:'fill',source:'countries',paint:{{
      'fill-color':['case',['boolean',['feature-state','hot'],false],'{hilite}',['get','color']],
      'fill-opacity':0.62}}}},
    {{id:'edge',type:'line',source:'countries',paint:{{'line-color':'#40608f','line-width':0.7}}}}]}};
  var map=new maplibregl.Map({{container:'map',style:style,center:[20,30],zoom:1.6,attributionControl:true}});
  map.addControl(new maplibregl.NavigationControl());
  map.on('error',function(){{fail();}});
  map.on('load',function(){{
    if(BOUNDS){{try{{map.fitBounds(BOUNDS,{{padding:55,maxZoom:5,duration:0}});}}catch(e){{}}}}
    document.title='MAPOK';
  }});
  window.__hot=[];
  window.__setHot=function(list){{
    window.__hot.forEach(function(id){{map.setFeatureState({{source:'countries',id:id}},{{hot:false}});}});
    window.__hot=list||[];
    window.__hot.forEach(function(id){{map.setFeatureState({{source:'countries',id:id}},{{hot:true}});}});
  }};
 }}
}}catch(e){{fail();}}
</script></body></html>"""
