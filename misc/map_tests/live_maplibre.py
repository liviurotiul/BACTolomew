"""LIVE real MapLibre GL web map for the Map Lab - a proper interactive web map:
Carto raster basemap + a data-driven country choropleth + sample points, smooth
vector pan/zoom in a QWebEngineView. Needs a working GL stack (post-reboot) and
internet (maplibre-gl from CDN + Carto tiles)."""
from __future__ import annotations
import json
import tempfile

import fixture
import geo_data
from PyQt6 import QtCore
from PyQt6.QtWebEngineWidgets import QWebEngineView

_TMP: list = []


def _countries_geojson() -> dict:
    # Use high-resolution outlines so the choropleth edges match the map tiles,
    # but keep the payload small by including only countries that have samples.
    counts = geo_data.counts_by_geojson_name(fixture.COUNTRIES)
    src = geo_data._ASSETS / "world_countries_10m.geojson"
    if not src.exists():
        src = geo_data._ASSETS / "world_countries_50m.geojson"
    if not src.exists():
        src = geo_data.GEOJSON_PATH
    data = json.loads(src.read_text())
    feats = []
    for feat in data.get("features", []):
        props = feat.get("properties", {})
        name = props.get("NAME") or props.get("ADMIN") or props.get("SOVEREIGNT") or ""
        count = geo_data.match_count(name, counts)
        if count <= 0:
            continue
        feats.append({
            "type": "Feature",
            "properties": {"name": name, "count": count},
            "geometry": feat.get("geometry"),
        })
    return {"type": "FeatureCollection", "features": feats}


def _points_geojson() -> dict:
    feats = []
    for p in fixture.POINTS:
        feats.append({
            "type": "Feature",
            "properties": {"sample": p["sample"], "country": p["country"]},
            "geometry": {"type": "Point", "coordinates": [p["lon"], p["lat"]]},
        })
    return {"type": "FeatureCollection", "features": feats}


def _html() -> str:
    mx = max(1, fixture.MAX_COUNT)
    countries = json.dumps(_countries_geojson())
    points = json.dumps(_points_geojson())
    return """<!DOCTYPE html><html><head><meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<link href="https://unpkg.com/maplibre-gl@4.7.1/dist/maplibre-gl.css" rel="stylesheet"/>
<script src="https://unpkg.com/maplibre-gl@4.7.1/dist/maplibre-gl.js"></script>
<style>html,body,#map{margin:0;height:100%;width:100%}
.legend{position:absolute;left:10px;bottom:20px;background:rgba(255,255,255,.9);padding:8px 10px;border-radius:6px;
 font:12px sans-serif;color:#1f2d3d;box-shadow:0 1px 4px rgba(0,0,0,.2)}
.legend i{display:inline-block;width:22px;height:11px;margin-right:2px;vertical-align:middle}</style></head>
<body><div id="map"></div>
<div class="legend"><b>samples / country</b><br/><span id="lg"></span></div>
<script>
const COUNTRIES=__COUNTRIES__, POINTS=__POINTS__, MX=__MX__;
const ramp=[0,'#eef1f5',1,'#deebf7',Math.max(2,MX/8|0),'#9ecae1',Math.max(4,MX/3|0),'#4292c6',MX,'#084594'];
const style={version:8,sources:{
  carto:{type:'raster',tiles:['https://a.basemaps.cartocdn.com/light_all/{z}/{x}/{y}.png'],tileSize:256,
    attribution:'© OpenStreetMap © CARTO'},
  countries:{type:'geojson',data:COUNTRIES}, pts:{type:'geojson',data:POINTS}},
 layers:[{id:'bg',type:'raster',source:'carto'},
  {id:'choro',type:'fill',source:'countries',filter:['>',['get','count'],0],
   paint:{'fill-color':['interpolate',['linear'],['get','count'],...ramp],'fill-opacity':0.6}},
  {id:'edge',type:'line',source:'countries',filter:['>',['get','count'],0],
   paint:{'line-color':'#4a6fa5','line-width':0.6}},
  {id:'pts',type:'circle',source:'pts',
   paint:{'circle-radius':4,'circle-color':'#b11242','circle-stroke-color':'#fff','circle-stroke-width':1}}]};
const map=new maplibregl.Map({container:'map',style,center:[25,42],zoom:2.4,attributionControl:true});
map.addControl(new maplibregl.NavigationControl());
map.on('load',()=>{
  const b=new maplibregl.LngLatBounds();
  POINTS.features.forEach(f=>b.extend(f.geometry.coordinates));
  if(!b.isEmpty()) map.fitBounds(b,{padding:60,maxZoom:5,duration:0});
});
let lg=''; for(let i=0;i<ramp.length;i+=2){lg+='<i style="background:'+ramp[i+1]+'"></i>'+ramp[i]+' ';}
document.getElementById('lg').innerHTML=lg;
</script></body></html>""".replace("__COUNTRIES__", countries).replace("__POINTS__", points).replace("__MX__", str(mx))


def build_widget(parent=None):
    from PyQt6 import QtWebEngineCore

    view = QWebEngineView(parent)
    settings = view.settings()
    for attr in ("LocalContentCanAccessRemoteUrls", "LocalContentCanAccessFileUrls"):
        settings.setAttribute(getattr(QtWebEngineCore.QWebEngineSettings.WebAttribute, attr), True)
    tmp = tempfile.NamedTemporaryFile("w", suffix=".html", delete=False, encoding="utf-8")
    tmp.write(_html())
    tmp.flush()
    tmp.close()
    _TMP.append(tmp.name)
    view.load(QtCore.QUrl.fromLocalFile(tmp.name))
    return view
