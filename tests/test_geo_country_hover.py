"""Hovering a country on the Map hovers that country's samples in the window.

Every map backend reports the country with samples under the pointer
(``countryHovered``); the Geo panel turns it into a window hover of the
samples the map counts for that country (the window's sample set), with the
panel as origin. Leaving the country or the map clears it; a country without
samples sends nothing.

Synthetic run (tests/test_refresh_parity.py): S1 and S4 in Romania and S6 in
Italy have coordinates, S2 (Italy) and S3 (Spain) have a country only, S5 has
no location.
"""
from __future__ import annotations

import os
import shutil
import tempfile
import types
import unittest
from pathlib import Path
from typing import Any, List
from unittest import mock

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # run from tests/
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

from PyQt6 import QtCore, QtGui, QtWidgets

# The web map's module imports QtWebEngineWidgets, which must happen before the
# QApplication exists (no web view is created here, nothing is started).
try:
    from PyQt6 import QtWebEngineWidgets  # noqa: F401

    _WEB_IMPORT_ERROR = ""
except ImportError as exc:
    _WEB_IMPORT_ERROR = str(exc)

try:
    from tests.test_dock_host import _app, _spin
    from tests.test_panel_contract import _WindowCase
    from tests.test_run_data import build_parity_run
except ImportError:
    from test_dock_host import _app, _spin
    from test_panel_contract import _WindowCase
    from test_run_data import build_parity_run

from annotation_explorer_polars.docking import docking_available
from annotation_explorer_polars.model.filter_state import SampleSubset
from annotation_explorer_polars.model.query_engine import QueryEngine

# A point inside each country (lon, lat).
BUCHAREST = (26.1, 44.43)
MADRID = (-3.7, 40.4)
PARIS = (2.35, 48.85)
ATLANTIC = (-35.0, 35.0)
SINGAPORE = (103.82, 1.35)
LUXEMBOURG = (6.13, 49.61)
BISCAY = (-5.0, 44.3)  # just off Spain's north coast

POINTS = [
    {"sample": "S1", "lat": 44.43, "lon": 26.1, "country": "Romania"},
    {"sample": "S4", "lat": 45.75, "lon": 21.23, "country": "Romania"},
    {"sample": "S6", "lat": 41.9, "lon": 12.5, "country": "Italy"},
]
COUNTRIES = {"Italy": 1, "Spain": 1}


def _move(view: QtWidgets.QGraphicsView, scene_pos: QtCore.QPointF, buttons: Any = None) -> None:
    """A real MouseMove on the map's viewport at this scene point."""
    viewport = view.viewport()
    pos = view.mapFromScene(scene_pos)
    pressed = QtCore.Qt.MouseButton.NoButton if buttons is None else buttons
    event = QtGui.QMouseEvent(
        QtCore.QEvent.Type.MouseMove,
        QtCore.QPointF(pos),
        QtCore.QPointF(viewport.mapToGlobal(pos)),
        QtCore.Qt.MouseButton.NoButton,
        pressed,
        QtCore.Qt.KeyboardModifier.NoModifier,
    )
    QtWidgets.QApplication.sendEvent(viewport, event)


def _leave(view: QtWidgets.QGraphicsView) -> None:
    QtWidgets.QApplication.sendEvent(view.viewport(), QtCore.QEvent(QtCore.QEvent.Type.Leave))


class EngineCountrySamplesTests(unittest.TestCase):
    def setUp(self) -> None:
        self.output_dir = Path(tempfile.mkdtemp(prefix="bactolomew_geo_hover_"))
        self.addCleanup(shutil.rmtree, self.output_dir, True)
        self.engine = QueryEngine(build_parity_run(self.output_dir))

    def test_every_counted_sample_is_listed_under_its_country(self) -> None:
        samples = ["S1", "S2", "S3", "S4", "S5", "S6"]
        self.assertEqual(
            self.engine.map_country_samples(samples),
            {"Romania": ["S1", "S4"], "Italy": ["S2", "S6"], "Spain": ["S3"]},
        )
        # Only the given sample set.
        self.assertEqual(self.engine.map_country_samples(["S1", "S2", "S5"]), {"Romania": ["S1"], "Italy": ["S2"]})
        self.assertEqual(self.engine.map_country_samples([]), {})

    def test_the_payload_keeps_its_keys_and_is_computed_once(self) -> None:
        samples = ["S1", "S2", "S3"]
        with mock.patch.object(self.engine, "_compute_map_payload", wraps=self.engine._compute_map_payload) as compute:
            payload = self.engine.map_payload(samples)
            self.engine.map_country_samples(samples)
        self.assertEqual(compute.call_count, 1)
        self.assertEqual(set(payload), {"points", "countries", "maxCountryCount", "total"})


class _BackendCase(unittest.TestCase):
    """One offscreen map view, fed the synthetic payload."""

    view_class: Any = None
    # A country drawn a pixel or two wide at the default zoom.
    small_country = ("Singapore", "singapore", SINGAPORE)

    @classmethod
    def setUpClass(cls) -> None:
        cls.app = _app()

    def _view(self) -> Any:
        view = self._make_view()
        self.addCleanup(view.deleteLater)
        view.resize(800, 400)
        view.show()
        _spin(self.app, 0.02)
        view.update_data(POINTS, COUNTRIES, 1)
        _spin(self.app, 0.02)
        self.hovered: List[str] = []
        view.countryHovered.connect(self.hovered.append)
        return view

    def _make_view(self) -> Any:
        return self.view_class()

    def _scene(self, view: Any, lonlat: tuple) -> QtCore.QPointF:
        raise NotImplementedError

    def _zoom_in(self, view: Any) -> None:
        view.scale(4.0, 4.0)
        view.centerOn(self._scene(view, (12.0, 45.0)))

    def _check_backend(self) -> None:
        view = self._view()
        self._zoom_in(view)
        _move(view, self._scene(view, BUCHAREST))
        self.assertEqual(self.hovered, ["romania"])
        _move(view, self._scene(view, (26.2, 44.5)))  # same country: nothing new
        self.assertEqual(self.hovered, ["romania"])
        _move(view, self._scene(view, MADRID))  # a country with a country-only sample
        self.assertEqual(self.hovered, ["romania", "spain"])
        _move(view, self._scene(view, PARIS))  # no samples: the hover ends
        self.assertEqual(self.hovered, ["romania", "spain", ""])
        _move(view, self._scene(view, ATLANTIC))  # still none: nothing sent
        self.assertEqual(self.hovered, ["romania", "spain", ""])
        _move(view, self._scene(view, BUCHAREST))
        _move(view, self._scene(view, MADRID), buttons=QtCore.Qt.MouseButton.LeftButton)  # dragging the map
        self.assertEqual(self.hovered[-1], "romania")
        _leave(view)
        self.assertEqual(self.hovered[-1], "")
        self.assertEqual(self.hovered.count(""), 2)


    def _check_small_country(self) -> None:
        """Final review 2026-10-08: Singapore (one sample) is a pixel at the
        default zoom and could not be hovered. A country drawn that small
        hovers from a few pixels away; a large one only under the pointer."""
        name, key, lonlat = self.small_country
        view = self._view()
        view.update_data([], {name: 1, "Spain": 1}, 1)
        _spin(self.app, 0.02)
        scale = abs(float(view.transform().m11()))
        center = self._scene(view, lonlat)
        _move(view, QtCore.QPointF(center.x(), center.y() + 3.0 / scale))
        self.assertEqual(self.hovered, [key])
        _move(view, QtCore.QPointF(center.x(), center.y() + 12.0 / scale))
        self.assertEqual(self.hovered, [key, ""])
        _move(view, self._scene(view, BISCAY))
        self.assertEqual(self.hovered, [key, ""])
        _move(view, self._scene(view, MADRID))
        self.assertEqual(self.hovered, [key, "", "spain"])


class SmallCountryNearTests(unittest.TestCase):
    def test_only_small_rings_within_the_slop_nearest_first(self) -> None:
        from annotation_explorer_polars.geo_raster_map import small_country_near

        rings = [
            ("big", QtCore.QRectF(0.0, 0.0, 100.0, 100.0)),
            ("tiny", QtCore.QRectF(110.0, 50.0, 2.0, 2.0)),
            ("other", QtCore.QRectF(118.0, 50.0, 2.0, 2.0)),
        ]
        # Scale 1: 4 px slop, under 8 px counts as small.
        self.assertEqual(small_country_near(rings, QtCore.QPointF(103.0, 51.0), 1.0), "")
        self.assertEqual(small_country_near(rings, QtCore.QPointF(108.5, 51.0), 1.0), "tiny")
        self.assertEqual(small_country_near(rings, QtCore.QPointF(116.5, 51.0), 1.0), "other")
        self.assertEqual(small_country_near(rings, QtCore.QPointF(130.0, 51.0), 1.0), "")
        # Zoomed in 5x the rings are 10 px: hovered only under the pointer.
        self.assertEqual(small_country_near(rings, QtCore.QPointF(109.5, 51.0), 5.0), "")


class SvgMapCountryHoverTests(_BackendCase):
    def _make_view(self) -> Any:
        from annotation_explorer_polars.geo_svg_map import SvgGeoView

        return SvgGeoView()

    def _scene(self, view: Any, lonlat: tuple) -> QtCore.QPointF:
        from annotation_explorer_polars.geo_svg_map import _project

        return QtCore.QPointF(*_project(*lonlat))

    def test_hover_reports_the_country_with_samples(self) -> None:
        self._check_backend()

    def test_a_small_country_hovers_from_a_few_pixels_away(self) -> None:
        self._check_small_country()

    def test_a_hover_does_not_highlight_or_redraw_the_map(self) -> None:
        view = self._view()
        with mock.patch.object(view, "_rebuild", wraps=view._rebuild) as rebuild:
            _move(view, self._scene(view, BUCHAREST))
        rebuild.assert_not_called()
        self.assertEqual(self.hovered, ["romania"])
        self.assertEqual(view.highlighted_countries(), set())

    def test_new_data_drops_countries_without_samples(self) -> None:
        view = self._view()
        view.update_data([], {"Spain": 2}, 2)
        _move(view, self._scene(view, BUCHAREST))
        self.assertEqual(self.hovered, [])
        _move(view, self._scene(view, MADRID))
        self.assertEqual(self.hovered, ["spain"])


class RasterMapCountryHoverTests(_BackendCase):
    # The offline raster map's coarse outlines have no Singapore.
    small_country = ("Luxembourg", "luxembourg", LUXEMBOURG)

    def _make_view(self) -> Any:
        from annotation_explorer_polars.geo_raster_map import RasterGeoMapView

        return RasterGeoMapView()

    def _scene(self, view: Any, lonlat: tuple) -> QtCore.QPointF:
        from annotation_explorer_polars.geo_raster_map import _project

        return QtCore.QPointF(*_project(*lonlat))

    def test_hover_reports_the_country_with_samples(self) -> None:
        self._check_backend()

    def test_a_small_country_hovers_from_a_few_pixels_away(self) -> None:
        self._check_small_country()


class TileMapCountryHoverTests(_BackendCase):
    def _make_view(self) -> Any:
        from annotation_explorer_polars import geo_tile_map

        # No tile requests (no network in tests); the country shapes are local.
        patcher = mock.patch.object(geo_tile_map.TileMapView, "_load_visible", lambda self: None)
        patcher.start()
        self.addCleanup(patcher.stop)
        return geo_tile_map.TileMapView()

    def _scene(self, view: Any, lonlat: tuple) -> QtCore.QPointF:
        from annotation_explorer_polars.geo_tile_map import _lonlat_to_world

        return QtCore.QPointF(*_lonlat_to_world(lonlat[0], lonlat[1], view._z))

    def _zoom_in(self, view: Any) -> None:
        view.centerOn(self._scene(view, (12.0, 45.0)))

    def test_hover_reports_the_country_with_samples(self) -> None:
        self._check_backend()

    def test_a_small_country_hovers_from_a_few_pixels_away(self) -> None:
        self._check_small_country()


@unittest.skipIf(_WEB_IMPORT_ERROR, f"QtWebEngineWidgets not importable here: {_WEB_IMPORT_ERROR}")
class MapLibreBridgeTests(unittest.TestCase):
    """The web map reports the hovered country through its title; the Python
    side of that channel (a WebEngine view is not created offscreen)."""

    def test_hover_titles(self) -> None:
        from annotation_explorer_polars.geo_maplibre_map import _hover_from_title

        self.assertEqual(_hover_from_title("MAPHOVER:romania"), "romania")
        self.assertEqual(_hover_from_title("MAPHOVER:"), "")
        self.assertIsNone(_hover_from_title("MAPOK"))
        self.assertIsNone(_hover_from_title("MAPFAIL"))

    def test_titles_become_country_hovers(self) -> None:
        from annotation_explorer_polars.geo_maplibre_map import MapLibreGeoView

        hovered: List[str] = []
        fake = types.SimpleNamespace(
            _counts={"romania": 2, "spain": 1, "france": 0},
            _hover_key="",
            countryHovered=types.SimpleNamespace(emit=hovered.append),
            _ready=True, _data_dirty=False, _pending_hot=[],
        )
        fake._set_hover_key = types.MethodType(MapLibreGeoView._set_hover_key, fake)
        on_title = types.MethodType(MapLibreGeoView._on_title, fake)
        for title in ("MAPHOVER:romania", "MAPHOVER:romania", "MAPHOVER:spain", "MAPHOVER:france",
                      "MAPHOVER:", "MAPHOVER:nowhere", "MAPHOVER:spain", "MAPHOVER:"):
            on_title(title)
        self.assertEqual(hovered, ["romania", "spain", "", "spain", ""])

    def test_the_page_reports_hovers_offline(self) -> None:
        from annotation_explorer_polars.geo_maplibre_map import _TEMPLATE

        self.assertIn("document.title='MAPHOVER:'+id", _TEMPLATE)
        self.assertIn("map.on('mousemove'", _TEMPLATE)
        self.assertIn("map.on('mouseout',", _TEMPLATE)
        self.assertIn("requestAnimationFrame(pickHover)", _TEMPLATE)
        self.assertIn("if(key in current)id=key", _TEMPLATE)  # only countries with samples
        self.assertNotIn("qwebchannel", _TEMPLATE.lower())
        self.assertNotIn("http", _TEMPLATE)

    def test_the_page_does_not_hover_while_dragging(self) -> None:
        from annotation_explorer_polars.geo_maplibre_map import _TEMPLATE

        # MapLibre fires mousemove during a drag; the page skips those moves
        # and picks the country again when the drag ends.
        self.assertIn("if(dragging||(e.originalEvent&&e.originalEvent.buttons))return;", _TEMPLATE)
        self.assertIn("map.on('dragstart',function(){{dragging=true;leave();}})", _TEMPLATE)
        self.assertIn("map.on('dragend',", _TEMPLATE)
        self.assertIn("window.__resetHover=function()", _TEMPLATE)

    def test_reset_hover_forgets_the_country_here_and_in_the_page(self) -> None:
        from annotation_explorer_polars.geo_maplibre_map import MapLibreGeoView

        scripts: List[str] = []
        hovered: List[str] = []
        page = types.SimpleNamespace(runJavaScript=scripts.append)
        fake = types.SimpleNamespace(
            _counts={"romania": 2}, _hover_key="", _ready=True, _data_dirty=False, _pending_hot=[],
            countryHovered=types.SimpleNamespace(emit=hovered.append), page=lambda: page,
        )
        fake._set_hover_key = types.MethodType(MapLibreGeoView._set_hover_key, fake)
        on_title = types.MethodType(MapLibreGeoView._on_title, fake)
        on_title("MAPHOVER:romania")
        MapLibreGeoView.reset_hover(fake)
        self.assertEqual(fake._hover_key, "")
        self.assertEqual(hovered, ["romania"])  # a reset does not emit
        self.assertEqual(scripts, ["window.__resetHover && window.__resetHover();"])
        on_title("MAPHOVER:")  # the page's reset title
        on_title("MAPHOVER:romania")  # the same country again is reported again
        self.assertEqual(hovered, ["romania", "romania"])


class ContainerRelayTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = _app()

    def test_the_container_relays_the_active_backend(self) -> None:
        from annotation_explorer_polars import geo_experimental_views as gev

        with mock.patch.dict(os.environ, {"BACTOLOMEW_SVG_MAP": "1"}):
            container = gev.GeoMapLibreView()
        self.addCleanup(container.deleteLater)
        hovered: List[str] = []
        container.countryHovered.connect(hovered.append)
        container._active.countryHovered.emit("romania")
        self.assertEqual(hovered, ["romania"])
        old = container._active
        container._activate_svg()  # a backend switch ends the old map's hover
        self.assertEqual(hovered, ["romania", ""])
        old.countryHovered.emit("spain")  # the old map is disconnected
        container._active.countryHovered.emit("italy")
        self.assertEqual(hovered, ["romania", "", "italy"])


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class GeoPanelCountryHoverTests(_WindowCase):
    """The real window (WebEngine off, so the Map is the offline SVG map)."""

    def _map(self, win: Any) -> Any:
        from annotation_explorer_polars.geo_svg_map import SvgGeoView

        manager = win._panel_manager
        self.assertEqual(manager.focus_panel("geo"), "geo")
        self.assertTrue(self._wait_for(lambda: manager.panel("geo") is not None and not manager.frame("geo").is_dirty(), 10.0))
        self._idle(win)
        panel = manager.panel("geo")
        view = panel.view._active
        self.assertIsInstance(view, SvgGeoView)
        view.scale(4.0, 4.0)
        view.centerOn(self._scene((12.0, 45.0)))
        return panel, view

    @staticmethod
    def _scene(lonlat: tuple) -> QtCore.QPointF:
        from annotation_explorer_polars.geo_svg_map import _project

        return QtCore.QPointF(*_project(*lonlat))

    def _hover_now(self, win: Any) -> tuple:
        win.ctx.flush_hover()
        return win.ctx.hover_state

    def test_hovering_a_country_hovers_its_samples(self) -> None:
        win = self._window()
        panel, view = self._map(win)
        published: List[Any] = []
        win.ctx.hoverChanged.connect(published.append)

        _move(view, self._scene(BUCHAREST))
        state = self._hover_now(win)
        self.assertEqual(state.samples, ("S1", "S4"))
        self.assertEqual(state.origin, "geo")
        self.assertEqual(state.kind, "samples")
        # No echo into the map's own highlight.
        _spin(self.app, 0.1)
        self.assertEqual(view.highlighted_countries(), set())

        _move(view, self._scene(MADRID))  # a sample with a country but no coordinates
        state = self._hover_now(win)
        self.assertEqual(state.samples, ("S3",))
        self.assertEqual(state.kind, "sample")

        _move(view, self._scene(PARIS))  # no samples: cleared
        self.assertEqual(self._hover_now(win).samples, ())
        count = len(published)
        _move(view, self._scene(ATLANTIC))
        _move(view, self._scene(PARIS))
        self._hover_now(win)
        _spin(self.app, 0.05)
        self.assertEqual(len(published), count)  # a country without samples sends nothing

        _move(view, self._scene(BUCHAREST))
        self.assertEqual(self._hover_now(win).samples, ("S1", "S4"))
        _leave(view)
        self.assertEqual(self._hover_now(win).samples, ())

    def test_leaving_the_panel_clears_the_hover(self) -> None:
        win = self._window()
        panel, view = self._map(win)
        _move(view, self._scene(BUCHAREST))
        self.assertEqual(self._hover_now(win).samples, ("S1", "S4"))
        QtWidgets.QApplication.sendEvent(panel, QtCore.QEvent(QtCore.QEvent.Type.Leave))
        self.assertEqual(self._hover_now(win).samples, ())

    def test_a_hidden_panel_reports_the_same_country_again(self) -> None:
        win = self._window()
        panel, view = self._map(win)
        _move(view, self._scene(BUCHAREST))
        self.assertEqual(self._hover_now(win).samples, ("S1", "S4"))
        panel.on_hidden()  # e.g. another panel maximized while the pointer rests here
        self.assertEqual(self._hover_now(win).samples, ())
        _move(view, self._scene(BUCHAREST))  # no Leave reached the map in between
        self.assertEqual(self._hover_now(win).samples, ("S1", "S4"))

    def test_only_the_samples_of_the_sample_set(self) -> None:
        win = self._window()
        panel, view = self._map(win)
        win.ctx.set_subset(SampleSubset(("S1", "S2", "S3"), label="three"))
        self._idle(win)
        _move(view, self._scene(BUCHAREST))
        self.assertEqual(self._hover_now(win).samples, ("S1",))
        _move(view, self._scene((12.5, 41.9)))  # Italy: S2 only (S6 is not in the set)
        self.assertEqual(self._hover_now(win).samples, ("S2",))
        # The set changes under a hovered country: the hover follows it.
        win.ctx.set_subset(None)
        self._idle(win)
        self.assertEqual(self._hover_now(win).samples, ("S2", "S6"))
        win.ctx.set_subset(SampleSubset(("S1",), label="one"))
        self._idle(win)
        self.assertEqual(self._hover_now(win).samples, ())

    def test_the_sample_tree_shows_the_country_hover(self) -> None:
        win = self._window()
        panel, view = self._map(win)
        tree = win._panel_manager.panel("sample_tree")
        self.assertIsNotNone(tree)
        got: List[Any] = []
        real = tree.show_hover

        def record(hover: Any) -> None:
            got.append(tuple(hover.samples))
            real(hover)

        tree.show_hover = record
        _move(view, self._scene(BUCHAREST))
        self.assertTrue(self._wait_for(lambda: ("S1", "S4") in got, 2.0), got)
        _leave(view)
        self.assertTrue(self._wait_for(lambda: got and got[-1] == (), 2.0), got)


if __name__ == "__main__":
    unittest.main()
