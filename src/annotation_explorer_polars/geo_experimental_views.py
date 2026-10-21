from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from typing import Dict, List, Optional, Sequence
import unicodedata

from PyQt6 import QtCore, QtWidgets
from .constants import COUNTRY_ALIASES


def _ordered_unique_samples(samples: Optional[Sequence[str]]) -> List[str]:
    if not samples:
        return []
    seen: set[str] = set()
    ordered: List[str] = []
    for raw in samples:
        sample = str(raw).strip()
        if not sample or sample in seen:
            continue
        seen.add(sample)
        ordered.append(sample)
    return ordered


def _map_status_text(points: List[dict], countries: Dict[str, int]) -> str:
    point_count = len(points)
    country_only = int(sum(int(v) for v in countries.values()))
    if point_count and country_only:
        return f"{point_count} mapped sample(s); {country_only} sample(s) only have country metadata."
    if point_count:
        return f"{point_count} mapped sample(s)."
    if country_only:
        return f"No coordinate-level points available. {country_only} sample(s) only have country metadata."
    return "No geospatial metadata available for the current filter."


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


_WEBENGINE_PROBE_CODE = """
try:
    import resource
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))  # no core file if we crash
except Exception:
    pass
from PyQt6 import QtCore, QtWidgets
QtCore.QCoreApplication.setAttribute(
    QtCore.Qt.ApplicationAttribute.AA_ShareOpenGLContexts,
    True,
)
app = QtWidgets.QApplication(["bactolomew-webengine-probe"])
from PyQt6.QtWebEngineWidgets import QWebEngineView
view = QWebEngineView()
# Realize the view (creates the native surface + OpenGL context) WITHOUT a
# visible window. On a display with broken GL this fails/crashes right here --
# but safely, because it is this throwaway subprocess, not the real app. That is
# the whole point: it must exercise the SAME real-display GL path the actual view
# will, so we never report OK and then crash the main process.
view.setAttribute(QtCore.Qt.WidgetAttribute.WA_DontShowOnScreen, True)
view.resize(48, 48)
view.show()
view.setHtml("<html><body>probe</body></html>")
for _ in range(25):
    app.processEvents(QtCore.QEventLoop.ProcessEventsFlag.AllEvents, 20)
print("WEBENGINE_OK")
""".strip()


_WEBENGINE_RUNTIME_CACHE: Optional[tuple[bool, str]] = None


def _cached_webengine_runtime_status() -> Optional[tuple[bool, str]]:
    return _WEBENGINE_RUNTIME_CACHE


def _webengine_runtime_status() -> tuple[bool, str]:
    global _WEBENGINE_RUNTIME_CACHE
    if _WEBENGINE_RUNTIME_CACHE is not None:
        return _WEBENGINE_RUNTIME_CACHE
    env = os.environ.copy()
    # Do NOT force offscreen: the probe must run on the same real display/GL as
    # the app, or WebEngine passes here and then crashes when the real view is
    # created on a machine whose hardware GL is broken.
    env.pop("QT_QPA_PLATFORM", None)
    try:
        proc = subprocess.run(
            [sys.executable, "-c", _WEBENGINE_PROBE_CODE],
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=15,
            check=False,
        )
    except Exception as exc:
        _WEBENGINE_RUNTIME_CACHE = (False, str(exc))
        return _WEBENGINE_RUNTIME_CACHE

    stdout = (proc.stdout or "").strip()
    stderr = (proc.stderr or "").strip()
    if proc.returncode == 0 and "WEBENGINE_OK" in stdout:
        _WEBENGINE_RUNTIME_CACHE = (True, "")
        return _WEBENGINE_RUNTIME_CACHE

    details = stderr or stdout or f"Qt WebEngine probe failed with exit code {proc.returncode}."
    details = details[-1200:]
    _WEBENGINE_RUNTIME_CACHE = (False, details)
    return _WEBENGINE_RUNTIME_CACHE


class _WebEngineProbeSignals(QtCore.QObject):
    finished = QtCore.pyqtSignal(bool, str)


class _WebEngineProbeWorker(QtCore.QRunnable):
    def __init__(self) -> None:
        super().__init__()
        self.signals = _WebEngineProbeSignals()
        self.setAutoDelete(True)

    def run(self) -> None:
        ok, details = _webengine_runtime_status()
        self.signals.finished.emit(bool(ok), str(details or ""))


def _maplibre_html(country_geojson_name: str, country_aliases_json: str) -> str:
    html = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>BACTolomew Geo Web</title>
  <link
    rel="stylesheet"
    href="https://unpkg.com/maplibre-gl@4.7.1/dist/maplibre-gl.css"
  />
  <style>
    html, body, #app, #map {
      margin: 0;
      padding: 0;
      width: 100%;
      height: 100%;
      overflow: hidden;
      background: #dbeafe;
      font-family: "Segoe UI", sans-serif;
    }
    #app {
      position: relative;
    }
    #map {
      position: absolute;
      inset: 0;
    }
    #banner {
      position: absolute;
      top: 12px;
      left: 12px;
      z-index: 10;
      max-width: min(520px, calc(100vw - 24px));
      padding: 10px 12px;
      border-radius: 10px;
      background: rgba(255, 255, 255, 0.92);
      color: #0f172a;
      box-shadow: 0 8px 28px rgba(15, 23, 42, 0.18);
      backdrop-filter: blur(6px);
      font-size: 13px;
      line-height: 1.35;
    }
    #banner.error {
      background: rgba(127, 29, 29, 0.92);
      color: #fff7ed;
    }
    #banner strong {
      display: block;
      margin-bottom: 4px;
      font-size: 13px;
      letter-spacing: 0.02em;
      text-transform: uppercase;
    }
    .maplibregl-popup-content {
      font: 13px/1.35 "Segoe UI", sans-serif;
      padding: 10px 12px;
    }
  </style>
</head>
<body>
  <div id="app">
    <div id="map"></div>
    <div id="banner"><strong>Geo Web</strong><span id="banner-text">Loading map…</span></div>
  </div>
  <script src="https://unpkg.com/maplibre-gl@4.7.1/dist/maplibre-gl.js"></script>
  <script>
    (function () {
      const banner = document.getElementById("banner");
      const bannerText = document.getElementById("banner-text");
      const state = {
        points: [],
        countries: {},
        countriesNormalized: {},
        maxCount: 1,
        highlighted: []
      };
      let map = null;
      let popup = null;
      let loaded = false;
      let countriesLoaded = false;
      let countriesGeojsonRaw = null;
      const sourceId = "bactolomew-samples";
      const layerId = "bactolomew-sample-circles";
      const strokeLayerId = "bactolomew-sample-strokes";
      const countrySourceId = "bactolomew-countries";
      const countryFillLayerId = "bactolomew-country-fill";
      const countryLineLayerId = "bactolomew-country-lines";
      const countryLabelLayerId = "bactolomew-country-labels";
      const countryGeojsonUrl = "__COUNTRY_GEOJSON__";
      const aliasMap = __COUNTRY_ALIASES__;

      function setBanner(text, isError) {
        bannerText.textContent = text || "";
        banner.classList.toggle("error", !!isError);
      }

      function normalizeCountryName(name) {
        if (!name) {
          return "";
        }
        let key = String(name)
          .normalize("NFKD")
          .replace(/[\\u0300-\\u036f]/g, "")
          .toLowerCase()
          .replace(/[^a-z0-9]+/g, " ")
          .trim();
        const alias = aliasMap[key];
        if (alias) {
          key = String(alias)
            .normalize("NFKD")
            .replace(/[\\u0300-\\u036f]/g, "")
            .toLowerCase()
            .replace(/[^a-z0-9]+/g, " ")
            .trim();
        }
        return key;
      }

      function bestCountryName(props) {
        if (!props) {
          return "Unknown";
        }
        return props.name || props.NAME || props.ADMIN || props.name_long || props.NAME_LONG || "Unknown";
      }

      function countryColor(count, maxCount) {
        if (!count || count <= 0) {
          return "#e5e7eb";
        }
        const denom = Math.max(1, Number(maxCount || 1));
        const t = count / denom;
        if (t > 0.8) return "#bd0026";
        if (t > 0.6) return "#f03b20";
        if (t > 0.4) return "#fd8d3c";
        if (t > 0.2) return "#feb24c";
        return "#fed976";
      }

      function geojsonFromState() {
        const highlighted = new Set((state.highlighted || []).map(String));
        return {
          type: "FeatureCollection",
          features: (state.points || []).map((point) => {
            const sample = String(point.sample || "");
            return {
              type: "Feature",
              properties: {
                sample: sample,
                count: Number(point.count || 0),
                selected: highlighted.has(sample),
                country: String(point.country || "")
              },
              geometry: {
                type: "Point",
                coordinates: [Number(point.lon), Number(point.lat)]
              }
            };
          })
        };
      }

      function countriesGeojsonFromState() {
        if (!countriesGeojsonRaw || !countriesGeojsonRaw.features) {
          return { type: "FeatureCollection", features: [] };
        }
        return {
          type: "FeatureCollection",
          features: countriesGeojsonRaw.features.map((feature) => {
            const props = Object.assign({}, feature.properties || {});
            const key = normalizeCountryName(bestCountryName(props));
            const count = Number((state.countriesNormalized || {})[key] || 0);
            props.__countryKey = key;
            props.__count = count;
            props.__color = countryColor(count, state.maxCount || 1);
            props.__label = bestCountryName(props);
            return {
              type: "Feature",
              properties: props,
              geometry: feature.geometry
            };
          })
        };
      }

      function computeBoundsFromGeojson(data) {
        if (!data || !data.features || !data.features.length) {
          return null;
        }
        let minLon = null;
        let minLat = null;
        let maxLon = null;
        let maxLat = null;

        function visitCoords(coords) {
          if (!Array.isArray(coords)) {
            return;
          }
          if (
            coords.length >= 2 &&
            typeof coords[0] === "number" &&
            typeof coords[1] === "number"
          ) {
            const lon = Number(coords[0]);
            const lat = Number(coords[1]);
            if (!Number.isFinite(lon) || !Number.isFinite(lat)) {
              return;
            }
            minLon = minLon === null ? lon : Math.min(minLon, lon);
            minLat = minLat === null ? lat : Math.min(minLat, lat);
            maxLon = maxLon === null ? lon : Math.max(maxLon, lon);
            maxLat = maxLat === null ? lat : Math.max(maxLat, lat);
            return;
          }
          coords.forEach(visitCoords);
        }

        data.features.forEach((feature) => {
          if (feature && feature.geometry && feature.geometry.coordinates) {
            visitCoords(feature.geometry.coordinates);
          }
        });
        if (minLon === null || minLat === null || maxLon === null || maxLat === null) {
          return null;
        }
        return [[minLon, minLat], [maxLon, maxLat]];
      }

      function fitToPoints() {
        if (!map) {
          return;
        }
        if (!state.points || state.points.length === 0) {
          if (countriesLoaded && state.countriesNormalized && Object.keys(state.countriesNormalized).length > 0) {
            const activeCountries = countriesGeojsonFromState();
            const features = activeCountries.features.filter(
              (feature) => Number((feature.properties || {}).__count || 0) > 0
            );
            const bounds = computeBoundsFromGeojson({ type: "FeatureCollection", features: features });
            if (bounds) {
              map.fitBounds(bounds, { padding: 50, duration: 650, maxZoom: 4.6 });
              return;
            }
          }
          map.easeTo({ center: [8, 18], zoom: 1.15, duration: 600 });
          return;
        }
        const bounds = new maplibregl.LngLatBounds();
        state.points.forEach((point) => bounds.extend([Number(point.lon), Number(point.lat)]));
        if (state.points.length === 1) {
          map.easeTo({ center: bounds.getCenter(), zoom: 4.8, duration: 600 });
          return;
        }
        map.fitBounds(bounds, { padding: 50, duration: 650, maxZoom: 5.8 });
      }

      function renderData() {
        if (!map || !loaded) {
          return;
        }
        const data = geojsonFromState();
        const source = map.getSource(sourceId);
        if (source) {
          source.setData(data);
        }
        if (countriesLoaded) {
          const countrySource = map.getSource(countrySourceId);
          if (countrySource) {
            countrySource.setData(countriesGeojsonFromState());
          }
        }
        const mappedCount = (state.points || []).length;
        if (mappedCount > 0) {
          setBanner(
            mappedCount + " mapped sample(s)." +
              (Object.keys(state.countries || {}).length ? " Some samples only have country metadata." : ""),
            false
          );
        } else if (Object.keys(state.countries || {}).length) {
          setBanner("No coordinate-level points available for this filter.", false);
        } else {
          setBanner("No geospatial metadata available for this filter.", false);
        }
      }

      function loadCountries() {
        return fetch(countryGeojsonUrl)
          .then((response) => {
            if (!response.ok) {
              throw new Error("HTTP " + response.status);
            }
            return response.json();
          })
          .then((data) => {
            countriesGeojsonRaw = data;
            countriesLoaded = true;
          })
          .catch((error) => {
            countriesLoaded = false;
            setBanner("Country overlay failed to load: " + String(error), true);
          });
      }

      function ensureMap() {
        if (typeof maplibregl === "undefined") {
          setBanner(
            "MapLibre assets failed to load. This tab needs network access to the demo style/CDN.",
            true
          );
          return;
        }
        map = new maplibregl.Map({
          container: "map",
          style: {
            version: 8,
            glyphs: "https://demotiles.maplibre.org/font/{fontstack}/{range}.pbf",
            sources: {},
            layers: [
              {
                id: "background",
                type: "background",
                paint: {
                  "background-color": "#dbeafe"
                }
              }
            ]
          },
          center: [8, 18],
          zoom: 1.15,
          attributionControl: false
        });
        map.addControl(new maplibregl.NavigationControl({ visualizePitch: false }), "top-right");
        map.on("load", function () {
          loaded = true;
          loadCountries().finally(function () {
            if (countriesLoaded) {
              map.addSource(countrySourceId, {
                type: "geojson",
                data: countriesGeojsonFromState()
              });
              map.addLayer({
                id: countryFillLayerId,
                type: "fill",
                source: countrySourceId,
                paint: {
                  "fill-color": ["get", "__color"],
                  "fill-opacity": 0.82
                }
              });
              map.addLayer({
                id: countryLineLayerId,
                type: "line",
                source: countrySourceId,
                paint: {
                  "line-color": "#334155",
                  "line-width": [
                    "interpolate",
                    ["linear"],
                    ["zoom"],
                    0, 0.55,
                    2, 0.75,
                    4, 1.1
                  ],
                  "line-opacity": 0.92
                }
              });
              map.addLayer({
                id: countryLabelLayerId,
                type: "symbol",
                source: countrySourceId,
                layout: {
                  "text-field": ["get", "__label"],
                  "text-font": ["Open Sans Semibold"],
                  "text-size": [
                    "interpolate",
                    ["linear"],
                    ["zoom"],
                    0, 10,
                    2, 11.5,
                    4, 13
                  ],
                  "text-max-width": 9,
                  "text-allow-overlap": false,
                  "text-ignore-placement": false,
                  "symbol-placement": "point"
                },
                paint: {
                  "text-color": "#0f172a",
                  "text-halo-color": "#ffffff",
                  "text-halo-width": 1.8,
                  "text-halo-blur": 0.6,
                  "text-opacity": 0.97
                }
              });
              map.on("click", countryFillLayerId, function (e) {
                const feature = e.features && e.features[0];
                if (!feature) {
                  return;
                }
                const props = feature.properties || {};
                if (popup) {
                  popup.remove();
                }
                popup = new maplibregl.Popup({ closeButton: false, offset: 12 })
                  .setLngLat(e.lngLat)
                  .setHTML(
                    "<strong>" + String(bestCountryName(props)) + "</strong><br/>" +
                    "count: " + String(props.__count || 0)
                  )
                  .addTo(map);
              });
              map.on("mouseenter", countryFillLayerId, function () {
                map.getCanvas().style.cursor = "pointer";
              });
              map.on("mouseleave", countryFillLayerId, function () {
                map.getCanvas().style.cursor = "";
              });
            }
            map.addSource(sourceId, {
              type: "geojson",
              data: geojsonFromState()
            });
            map.addLayer({
              id: layerId,
              type: "circle",
              source: sourceId,
              paint: {
                "circle-radius": [
                  "case",
                  ["boolean", ["get", "selected"], false],
                  9,
                  6
                ],
                "circle-color": [
                  "case",
                  ["boolean", ["get", "selected"], false],
                  "#f97316",
                  "#2563eb"
                ],
                "circle-opacity": 0.88
              }
            });
            map.addLayer({
              id: strokeLayerId,
              type: "circle",
              source: sourceId,
              paint: {
                "circle-radius": [
                  "case",
                  ["boolean", ["get", "selected"], false],
                  11.5,
                  7.8
                ],
                "circle-stroke-width": [
                  "case",
                  ["boolean", ["get", "selected"], false],
                  2.8,
                  1.8
                ],
                "circle-stroke-color": [
                  "case",
                  ["boolean", ["get", "selected"], false],
                  "#7c2d12",
                  "#ffffff"
                ],
                "circle-color": "rgba(0,0,0,0)"
              }
            });
            map.on("click", layerId, function (e) {
              const feature = e.features && e.features[0];
              if (!feature) {
                return;
              }
              const props = feature.properties || {};
              const coordinates = feature.geometry && feature.geometry.coordinates;
              if (!coordinates) {
                return;
              }
              if (popup) {
                popup.remove();
              }
              popup = new maplibregl.Popup({ closeButton: false, offset: 12 })
                .setLngLat(coordinates)
                .setHTML(
                  "<strong>" + String(props.sample || "Sample") + "</strong><br/>" +
                  (props.country ? String(props.country) + "<br/>" : "") +
                  "count: " + String(props.count || 0)
                )
                .addTo(map);
            });
            map.on("mouseenter", layerId, function () {
              map.getCanvas().style.cursor = "pointer";
            });
            map.on("mouseleave", layerId, function () {
              map.getCanvas().style.cursor = "";
            });
            renderData();
            fitToPoints();
          });
        });
        map.on("error", function (event) {
          if (event && event.error) {
            setBanner("Map rendering error: " + String(event.error), true);
          }
        });
      }

      window.__setBACTolomewPayload = function (payload) {
        payload = payload || {};
        state.points = Array.isArray(payload.points) ? payload.points : [];
        state.countries = payload.countries || {};
        state.countriesNormalized = payload.countriesNormalized || {};
        state.maxCount = Number(payload.maxCount || 1);
        renderData();
      };

      window.__setBACTolomewHighlightedSamples = function (samples) {
        state.highlighted = Array.isArray(samples) ? samples.map(String) : [];
        renderData();
      };

      ensureMap();
      window.setTimeout(function () {
        if (!map) {
          setBanner(
            "Map view is unavailable. If the page stays blank, the embedded browser or network resources are missing.",
            true
          );
        }
      }, 5000);
      window.setTimeout(function () {
        if (map && loaded) {
          fitToPoints();
        }
      }, 350);
      window.__fitBACTolomewPoints = fitToPoints;
    })();
  </script>
</body>
</html>
"""
    html = html.replace('"__COUNTRY_GEOJSON__"', json.dumps(country_geojson_name))
    html = html.replace("__COUNTRY_ALIASES__", country_aliases_json)
    return html


class GeoMapLibreView(QtWidgets.QWidget):
    """Geo-panel map container: the MapLibre GL web map when a GPU/WebEngine is
    available, otherwise the offline SVG choropleth. Same API the panel expects:
    update_data(points, countries, max_count) + highlight_sample(s), and
    countryHovered(key): the normalized name of the country with samples under
    the pointer in the active map ("" for none), relayed from the backend.

    Force the offline SVG map with BACTOLOMEW_SVG_MAP=1.
    """

    countryHovered = QtCore.pyqtSignal(str)

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self._active = None
        self._is_web = False
        self._last_update: Optional[tuple] = None
        self._last_highlight: List[str] = []
        self._probe_worker: Optional[_WebEngineProbeWorker] = None
        self._message_label: Optional[QtWidgets.QLabel] = None

        self._layout = QtWidgets.QVBoxLayout(self)
        self._layout.setContentsMargins(0, 0, 0, 0)

        if os.environ.get("BACTOLOMEW_SVG_MAP", "") in ("1", "true", "True"):
            self._activate_svg()
            return

        cached = _cached_webengine_runtime_status()
        if cached is None:
            self._message_label = self._build_message_label("Loading map…")
            self._layout.addWidget(self._message_label)
            worker = _WebEngineProbeWorker()
            worker.signals.finished.connect(self._on_webengine_probe_finished)
            self._probe_worker = worker
            QtCore.QThreadPool.globalInstance().start(worker)
            return
        self._on_webengine_probe_finished(bool(cached[0]), str(cached[1] or ""))

    # ---- backend selection ----
    def _on_webengine_probe_finished(self, ok: bool, details: str) -> None:
        self._probe_worker = None
        if self._active is not None:
            return
        if ok:
            self._activate_maplibre()
        else:
            self._activate_svg()

    def _activate_maplibre(self) -> None:
        try:
            from .geo_maplibre_map import MapLibreGeoView

            view = MapLibreGeoView(self)
        except Exception:
            self._activate_svg()
            return
        view.renderUnavailable.connect(self._on_maplibre_failed)
        self._set_active(view, is_web=True)

    def _on_maplibre_failed(self, _reason: str) -> None:
        # MapLibre couldn't render (WebGL/network) - drop to the offline SVG map.
        if self._is_web:
            self._activate_svg()

    def _activate_svg(self) -> None:
        try:
            from .geo_svg_map import SvgGeoView

            view = SvgGeoView(self)
        except Exception as exc:  # noqa: BLE001
            self._set_message(f"The map failed to load:\n\n{exc}")
            return
        self._set_active(view, is_web=False)

    def _set_active(self, view: QtWidgets.QWidget, *, is_web: bool) -> None:
        if self._message_label is not None:
            self._layout.removeWidget(self._message_label)
            self._message_label.deleteLater()
            self._message_label = None
        if self._active is not None:
            self._layout.removeWidget(self._active)
            try:
                self._active.countryHovered.disconnect(self.countryHovered)
            except (AttributeError, TypeError):
                pass
            self._active.deleteLater()
            self.countryHovered.emit("")  # the old map's hover ends with it
        self._active = view
        hovered = getattr(view, "countryHovered", None)
        if hovered is not None:
            hovered.connect(self.countryHovered)
        self._is_web = is_web
        self._layout.addWidget(view)
        if self._last_update is not None:
            points, countries, max_count, total = self._last_update
            view.update_data(points, countries, max_count, total=total)
        if self._last_highlight:
            view.highlight_samples(self._last_highlight)

    def _set_message(self, text: str) -> None:
        if self._message_label is None:
            self._message_label = self._build_message_label(text)
            self._layout.addWidget(self._message_label)
        else:
            self._message_label.setText(text)

    @staticmethod
    def _build_message_label(text: str) -> QtWidgets.QLabel:
        label = QtWidgets.QLabel(text)
        label.setWordWrap(True)
        label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        label.setStyleSheet("padding: 18px; color: #64748b;")
        return label

    # ---- API forwarded to the active backend ----
    def update_data(self, points: List[dict], countries: Dict[str, int], max_count: int,
                    total: Optional[int] = None) -> None:
        self._last_update = (list(points or []), dict(countries or {}), int(max_count or 1), total)
        if self._active is not None:
            self._active.update_data(points, countries, max_count, total=total)

    def highlight_sample(self, sample: Optional[str]) -> None:
        self.highlight_samples([sample] if sample is not None and str(sample).strip() else [])

    def highlight_samples(self, samples: Optional[Sequence[str]]) -> None:
        self._last_highlight = _ordered_unique_samples(samples)
        if self._active is not None:
            self._active.highlight_samples(self._last_highlight)

    def reset_hover(self) -> None:
        """Forget the active map's hovered country (see countryHovered)."""
        reset = getattr(self._active, "reset_hover", None)
        if reset is not None:
            reset()


class _GeoQmlBridge(QtCore.QObject):
    pointsChanged = QtCore.pyqtSignal()
    highlightSamplesChanged = QtCore.pyqtSignal()
    statusTextChanged = QtCore.pyqtSignal()

    def __init__(self, parent: Optional[QtCore.QObject] = None) -> None:
        super().__init__(parent)
        self._points: List[dict] = []
        self._highlight_samples: List[str] = []
        self._status_text = ""

    @QtCore.pyqtProperty("QVariantList", notify=pointsChanged)
    def points(self) -> List[dict]:
        return list(self._points)

    @QtCore.pyqtProperty("QVariantList", notify=highlightSamplesChanged)
    def highlightSamples(self) -> List[str]:
        return list(self._highlight_samples)

    @QtCore.pyqtProperty(str, notify=statusTextChanged)
    def statusText(self) -> str:
        return self._status_text

    def set_payload(self, points: List[dict], countries: Dict[str, int], _max_count: int) -> None:
        self._points = list(points)
        self._status_text = _map_status_text(points, countries)
        self.pointsChanged.emit()
        self.statusTextChanged.emit()

    def set_highlight_samples(self, samples: Optional[Sequence[str]]) -> None:
        self._highlight_samples = _ordered_unique_samples(samples)
        self.highlightSamplesChanged.emit()


class GeoQtLocationView(QtWidgets.QWidget):
    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self._bridge = _GeoQmlBridge(self)
        self._message_label: Optional[QtWidgets.QLabel] = None
        self._quick = None

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        try:
            from PyQt6.QtQuickWidgets import QQuickWidget
        except Exception as exc:
            self._message_label = self._build_message_label(
                f"Qt Quick/QML is not available in this environment: {exc}"
            )
            layout.addWidget(self._message_label)
            return
        missing_runtime: List[str] = []
        for mod_name in ("PyQt6.QtQml", "PyQt6.QtPositioning", "PyQt6.QtLocation"):
            if importlib.util.find_spec(mod_name) is None:
                missing_runtime.append(mod_name.rsplit(".", 1)[-1])
        if missing_runtime:
            self._message_label = self._build_message_label(
                "GeoQML is unavailable because these Qt modules are not installed: "
                + ", ".join(missing_runtime)
            )
            layout.addWidget(self._message_label)
            return

        quick = QQuickWidget(self)
        quick.setResizeMode(QQuickWidget.ResizeMode.SizeRootObjectToView)
        quick.statusChanged.connect(self._on_quick_status_changed)
        quick.rootContext().setContextProperty("bridge", self._bridge)
        qml_path = Path(__file__).resolve().parent / "assets" / "geo_qtlocation.qml"
        quick.setSource(QtCore.QUrl.fromLocalFile(str(qml_path.resolve())))
        layout.addWidget(quick)
        self._quick = quick

    @staticmethod
    def _build_message_label(text: str) -> QtWidgets.QLabel:
        label = QtWidgets.QLabel(text)
        label.setWordWrap(True)
        label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        label.setStyleSheet("padding: 18px; color: #7f1d1d;")
        return label

    def _on_quick_status_changed(self, status: object) -> None:
        if self._quick is None:
            return
        try:
            from PyQt6.QtQuickWidgets import QQuickWidget
            is_error = status == QQuickWidget.Status.Error
        except Exception:
            is_error = False
        if not is_error:
            return
        errors = self._quick.errors()
        details = "\n".join(str(err) for err in errors) if errors else "Unknown QML error."
        self._quick.hide()
        if self._message_label is None:
            self._message_label = self._build_message_label(details)
            self.layout().addWidget(self._message_label)  # type: ignore[union-attr]
        else:
            self._message_label.setText(details)

    def update_data(self, points: List[dict], countries: Dict[str, int], max_count: int) -> None:
        self._bridge.set_payload(points, countries, max_count)

    def highlight_sample(self, sample: Optional[str]) -> None:
        samples = [str(sample).strip()] if sample is not None and str(sample).strip() else []
        self.highlight_samples(samples)

    def highlight_samples(self, samples: Optional[Sequence[str]]) -> None:
        self._bridge.set_highlight_samples(samples)
