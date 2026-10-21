"""MapLibre GL choropleth map: the primary Geo-panel map when a GPU/WebEngine is
available (country choropleth coloured by sample count, NO points). Runs in a
QWebEngineView.

Fully offline: the MapLibre library (assets/maplibre/) and the country shapes
(assets/world_countries_*.geojson) ship with the app, so the map needs no
internet, tile service or API key. The page is loaded once per view; later
update_data() calls recolour the live map instead of reloading it.

Emits countryHovered(str) with the normalized name of the country with samples
under the pointer ("" when the pointer is over no such country or leaves the
map). The page reports it through its title ("MAPHOVER:<id>"), the same local
channel as MAPOK/MAPFAIL, so it needs no QWebChannel and works offline.

Emits renderUnavailable(str) only if the map cannot start at all (WebGL or the
page fails), so the panel can swap to the offline SVG map. (Creating the WebView
itself is gated by the panel's WebEngine probe, so we never realize it on a
machine with broken GL.)
"""
from __future__ import annotations

import json
import math
import os
import tempfile
from typing import Dict, List, Optional, Sequence

from PyQt6 import QtCore, QtWebEngineCore, QtWidgets
from PyQt6.QtWebEngineWidgets import QWebEngineView

from .geo_raster_map import _ASSETS, _count_color, _norm
from .theme import COLORS

MAPLIBRE_JS = (_ASSETS / "maplibre" / "maplibre-gl.js").as_uri()
MAPLIBRE_CSS = (_ASSETS / "maplibre" / "maplibre-gl.css").as_uri()
HILITE = COLORS["hover"]
WATER = "#dfe6eb"
LAND = "#f4f4f1"
BORDER = "#c3ccd5"
PIE_ON = "#386baa"
PIE_OFF = "#b8c0ca"

_HOVER_TITLE = "MAPHOVER:"

_WORLD_JSON: Optional[str] = None  # every country shape, built once and shared by all map views
_WORLD_IDS: frozenset = frozenset()


def _world_json() -> str:
    """All countries from the bundled geojson, with ids keyed like the sample counts (_norm)."""
    global _WORLD_JSON, _WORLD_IDS
    if _WORLD_JSON is None:
        src = next((_ASSETS / n for n in
                    ("world_countries_10m.geojson", "world_countries_50m.geojson", "world_countries.geojson")
                    if (_ASSETS / n).exists()), None)
        feats = []
        if src is not None:
            for feat in json.loads(src.read_text()).get("features", []):
                props = feat.get("properties", {})
                name = props.get("NAME") or props.get("ADMIN") or props.get("SOVEREIGNT") or ""
                feats.append({"type": "Feature", "properties": {"id": _norm(name), "name": name},
                              "geometry": feat.get("geometry")})
        _WORLD_IDS = frozenset(f["properties"]["id"] for f in feats if f["properties"]["id"])
        _WORLD_JSON = json.dumps({"type": "FeatureCollection", "features": feats})
    return _WORLD_JSON


def _world_ids() -> frozenset:
    _world_json()
    return _WORLD_IDS


def _hover_from_title(title: str) -> Optional[str]:
    """The country id a MAPHOVER title reports ("" for none), None for any other title."""
    if title.startswith(_HOVER_TITLE):
        return title[len(_HOVER_TITLE):].strip()
    return None


class MapLibreGeoView(QWebEngineView):
    renderUnavailable = QtCore.pyqtSignal(str)
    # The country under the pointer (normalized like the counts), "" for none.
    # Emitted only when it changes.
    countryHovered = QtCore.pyqtSignal(str)

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        settings = self.settings()
        for attr in ("LocalContentCanAccessRemoteUrls", "LocalContentCanAccessFileUrls"):
            settings.setAttribute(getattr(QtWebEngineCore.QWebEngineSettings.WebAttribute, attr), True)
        self._tmp: Optional[str] = None
        self._counts: Dict[str, int] = {}
        self._sample_country: Dict[str, str] = {}
        self._max = 1
        self._points: List[dict] = []
        self._total = 0
        self._on_map = 0
        self._page_requested = False
        self._ready = False
        self._data_dirty = False
        self._pending_hot: List[str] = []
        self._hover_key = ""
        self.loadFinished.connect(self._on_load)
        self.titleChanged.connect(self._on_title)

    # ---- API ----
    def update_data(self, points: List[dict], countries: Dict[str, int], max_count: int,
                    total: Optional[int] = None) -> None:
        # total: samples in the current selection, including those with no location at all
        # (they are not in points/countries); feeds the "on map / not on map" pie.
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
        located = len(points or []) + sum(int(c) for c in (countries or {}).values())
        self._total = located if total is None else int(total)
        ids = _world_ids()
        self._on_map = sum(c for key, c in counts.items() if key in ids)
        if not self._page_requested:
            self._load_page()
        elif self._ready:
            self._push_data()
        else:
            self._data_dirty = True  # pushed once the page reports MAPOK

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
    def _colors(self) -> Dict[str, str]:
        return {key: _count_color(count, self._max).name() for key, count in self._counts.items() if count > 0}

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

    def _pie_html(self) -> str:
        total = self._total
        if total <= 0:
            return ""
        on = min(self._on_map, total)
        circ = 2 * math.pi * 8  # an r=8 circle with a 16-wide stroke fills a radius-16 pie
        rows = ""
        for label, n, color in (("on map", on, PIE_ON), ("not on map", total - on, PIE_OFF)):
            pct = 100 * n / total
            rows += (f'<tr><td><i class="sw" style="background:{color}"></i>{label}</td>'
                     f'<td class="n">{n}</td><td class="n">{"&lt;1" if 0 < pct < 1 else round(pct)}%</td></tr>')
        tip = ("On map: the sample's country is drawn on this map. Not on map: no country in the "
               "metadata, or a country name the map does not recognise.")
        return (f'<svg width="40" height="40" viewBox="0 0 32 32"><title>{tip}</title>'
                f'<circle cx="16" cy="16" r="16" fill="{PIE_OFF}"/>'
                f'<circle cx="16" cy="16" r="8" fill="none" stroke="{PIE_ON}" stroke-width="16" '
                f'stroke-dasharray="{circ * on / total:.3f} {circ:.3f}" transform="rotate(-90 16 16)"/></svg>'
                f'<div title="{tip}"><b>geolocation</b><table>{rows}</table></div>')

    def _html(self) -> str:
        return _TEMPLATE.format(
            css=MAPLIBRE_CSS, js=MAPLIBRE_JS, hilite=HILITE, water=WATER, land=LAND, border=BORDER,
            world=_world_json(), colors=json.dumps(self._colors()),
            bounds=json.dumps(self._bounds()), legend=self._legend_html(), pie=self._pie_html(),
        )

    def _load_page(self) -> None:
        # Loaded once per view; every later update goes through _push_data.
        self._page_requested = True
        self._ready = False
        tmp = tempfile.NamedTemporaryFile("w", suffix=".html", delete=False, encoding="utf-8")
        tmp.write(self._html())
        tmp.flush()
        tmp.close()
        self._tmp = tmp.name
        self.load(QtCore.QUrl.fromLocalFile(tmp.name))

    def _push_data(self) -> None:
        self._data_dirty = False
        payload = {"colors": self._colors(), "bounds": self._bounds(), "legend": self._legend_html(),
                   "pie": self._pie_html()}
        self.page().runJavaScript(f"window.__setData && window.__setData({json.dumps(payload)});")

    def _on_load(self, ok: bool) -> None:
        # The page is in memory now, so the temp HTML is no longer needed.
        if self._tmp is not None:
            try:
                os.remove(self._tmp)
            except OSError:
                pass
            self._tmp = None
        if not ok and not self._ready:
            self.renderUnavailable.emit("page failed to load")

    def _on_title(self, title: str) -> None:
        hovered = _hover_from_title(title)
        if hovered is not None:
            self._set_hover_key(hovered if self._counts.get(hovered, 0) > 0 else "")
            return
        if title == "MAPFAIL":
            self.renderUnavailable.emit("map failed to start (WebGL or page init)")
        elif title == "MAPOK":
            self._ready = True
            if self._data_dirty:
                self._push_data()
            if self._pending_hot:
                self.page().runJavaScript(f"window.__setHot && window.__setHot({json.dumps(self._pending_hot)});")
                self._pending_hot = []


    def _set_hover_key(self, key: str) -> None:
        if key != self._hover_key:
            self._hover_key = key
            self.countryHovered.emit(key)

    def reset_hover(self) -> None:
        """Forget the hovered country without emitting (the panel already
        cleared its hover), here and in the page, so the next pointer move over
        the same country reports it again."""
        self._hover_key = ""
        if self._ready:
            self.page().runJavaScript("window.__resetHover && window.__resetHover();")

    def leaveEvent(self, event: QtCore.QEvent) -> None:
        super().leaveEvent(event)
        self._set_hover_key("")


def _has_coords(p: dict) -> bool:
    try:
        return not (p["lat"] is None or p["lon"] is None) and float(p["lat"]) == float(p["lat"])
    except (KeyError, TypeError, ValueError):
        return False


_TEMPLATE = """<!DOCTYPE html><html><head><meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<link href="{css}" rel="stylesheet"/><script src="{js}"></script>
<style>html,body,#map{{margin:0;height:100%;width:100%;background:{water}}}
.overlay{{position:absolute;left:10px;bottom:22px;display:flex;flex-direction:column;align-items:flex-start;gap:8px}}
.card{{background:rgba(255,255,255,.92);padding:8px 11px;border-radius:6px;
 font:12px sans-serif;color:#1f2d3d;box-shadow:0 1px 5px rgba(0,0,0,.2)}}
.card:empty{{display:none}}
.legend i{{display:inline-block;width:22px;height:11px;margin-right:2px;vertical-align:middle}}
.pie{{display:flex;align-items:center;gap:10px}}
.pie table{{border-collapse:collapse;margin-top:2px}}
.pie td{{padding:1px 0}}
.pie td.n{{text-align:right;padding-left:10px}}
.pie .sw{{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:5px;vertical-align:-1px}}</style></head>
<body><div id="map"></div><div class="overlay"><div class="card pie">{pie}</div><div class="card legend">{legend}</div></div>
<script>
var WORLD={world}, COLORS={colors}, BOUNDS={bounds};
// Everything is local, so the only failure is the map not starting at all (no WebGL, broken page).
var loaded=false;
function fail(){{if(!loaded)document.title='MAPFAIL';}}
setTimeout(function(){{if(!loaded)fail();}},30000);
try{{
 if(typeof maplibregl==='undefined'){{fail();}}   // bundled script missing
 else{{
  var colored=['to-boolean',['feature-state','color']];
  var hot=['boolean',['feature-state','hot'],false];
  var style={{version:8,sources:{{countries:{{type:'geojson',data:WORLD,promoteId:'id'}}}},
   layers:[{{id:'water',type:'background',paint:{{'background-color':'{water}'}}}},
    {{id:'land',type:'fill',source:'countries',paint:{{'fill-color':'{land}'}}}},
    {{id:'choro',type:'fill',source:'countries',paint:{{
      'fill-color':['case',hot,'{hilite}',['coalesce',['feature-state','color'],'{land}']],
      'fill-opacity':['case',hot,0.62,colored,0.62,0]}}}},
    {{id:'edge',type:'line',source:'countries',paint:{{
      'line-color':['case',colored,'#40608f','{border}'],
      'line-width':['case',colored,0.7,0.4]}}}}]}};
  var map=new maplibregl.Map({{container:'map',style:style,center:[20,30],zoom:1.6,attributionControl:false}});
  map.addControl(new maplibregl.NavigationControl());
  map.on('error',function(){{fail();}});  // no-op once the map has started
  var fitted=null;
  function fit(b){{  // refit only when the bounds change, so filtering keeps the user's view
    var key=JSON.stringify(b);
    if(key===fitted)return;
    fitted=key;
    try{{if(b)map.fitBounds(b,{{padding:55,maxZoom:5,duration:0}});else map.jumpTo({{center:[20,30],zoom:1.6}});}}catch(e){{}}
  }}
  var current={{}};
  function paint(colors){{  // recolour through feature-state; no shapes are resent
    Object.keys(current).forEach(function(id){{if(!(id in colors))map.setFeatureState({{source:'countries',id:id}},{{color:null}});}});
    Object.keys(colors).forEach(function(id){{map.setFeatureState({{source:'countries',id:id}},{{color:colors[id]}});}});
    current=colors;
  }}
  map.on('load',function(){{
    loaded=true;
    paint(COLORS);
    fit(BOUNDS);
    document.title='MAPOK';
  }});
  window.__setData=function(d){{  // in-place update, no page reload
    paint(d.colors);
    document.querySelector('.legend').innerHTML=d.legend;
    document.querySelector('.pie').innerHTML=d.pie;
    fit(d.bounds);
  }};
  // The country with samples under the pointer goes to Python through the
  // title (MAPHOVER:<id>, empty id for none), once per change and at most
  // once per frame; a country without samples counts as none.
  var hoverId='', hoverPoint=null, hoverFrame=0;
  function sendHover(id){{
    if(id===hoverId)return;
    hoverId=id;
    document.title='MAPHOVER:'+id;
  }}
  function pickHover(){{
    hoverFrame=0;
    if(!hoverPoint){{sendHover('');return;}}
    var id='';
    try{{
      var fs=map.queryRenderedFeatures(hoverPoint,{{layers:['land']}});
      if(fs.length){{var key=String(fs[0].properties.id||'');if(key in current)id=key;}}
    }}catch(e){{}}
    sendHover(id);
  }}
  // MapLibre fires mousemove during a drag too: no hover while a button is
  // down or the map is dragged; the country under the pointer is picked again
  // when the drag ends.
  var dragging=false, lastPoint=null;
  function schedule(){{if(!hoverFrame)hoverFrame=requestAnimationFrame(pickHover);}}
  map.on('mousemove',function(e){{
    lastPoint=e.point;
    if(dragging||(e.originalEvent&&e.originalEvent.buttons))return;
    hoverPoint=e.point;
    schedule();
  }});
  function leave(){{hoverPoint=null;if(hoverFrame){{cancelAnimationFrame(hoverFrame);hoverFrame=0;}}sendHover('');}}
  map.on('mouseout',function(){{lastPoint=null;leave();}});
  map.on('dragstart',function(){{dragging=true;leave();}});
  map.on('dragend',function(e){{
    dragging=false;
    var o=e&&e.originalEvent;  // the release; its position is where the pointer rests
    if(o&&typeof o.clientX==='number'){{
      var r=map.getCanvasContainer().getBoundingClientRect();
      lastPoint=[o.clientX-r.left,o.clientY-r.top];
    }}
    if(lastPoint){{hoverPoint=lastPoint;schedule();}}
  }});
  // The panel cleared its hover (hidden, left): forget the last country so the
  // next pointer move over it reports it again.
  window.__resetHover=function(){{
    hoverPoint=null;if(hoverFrame){{cancelAnimationFrame(hoverFrame);hoverFrame=0;}}
    hoverId='';document.title='MAPHOVER:';
  }};
  window.__hot=[];
  window.__setHot=function(list){{
    window.__hot.forEach(function(id){{map.setFeatureState({{source:'countries',id:id}},{{hot:false}});}});
    window.__hot=list||[];
    window.__hot.forEach(function(id){{map.setFeatureState({{source:'countries',id:id}},{{hot:true}});}});
  }};
 }}
}}catch(e){{fail();}}
</script></body></html>"""
