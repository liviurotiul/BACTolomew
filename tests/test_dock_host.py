"""Docking area tests (window architecture plan, Phase 1).

Offscreen. Skipped when PyQt6Ads is not installed. The first class tests
``DockHost`` on its own (a rail plus the docking area in a plain window); the
second builds the real ``ExploreWindow`` on the synthetic run of the
golden tests and checks the panels, menus, lazy build, visibility gating and
layout persistence. Real-input drags need a real X server and are covered by
the Phase 0 spike and checklist C; here drags are scripted with QTest.
"""
from __future__ import annotations

import gc
import os
import shutil
import tempfile
import time
import unittest
from unittest import mock
from pathlib import Path
from typing import Any, List, Optional

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets, sip
from PyQt6.QtTest import QTest

from annotation_explorer_polars.docking import docking_available

try:
    from tests._pan_fixture import write_pan_outputs
    from tests.test_refresh_parity import build_synthetic_run
except ImportError:
    from _pan_fixture import write_pan_outputs
    from test_refresh_parity import build_synthetic_run

if docking_available():
    from annotation_explorer_polars.docking import ads_config
    from annotation_explorer_polars.docking.dock_host import DockHost, Placement, ads
    from annotation_explorer_polars.panels import manager as manager_module
    from annotation_explorer_polars.panels.base import PanelFrame, PanelSpec
    from annotation_explorer_polars.panels.manager import CLOSED_STACK_LIMIT, PanelManager

CLASSIC_GROUPS = [["sample_tree"], ["annotations"], ["track"], ["geo"], ["phylogeny"], ["gene_presence_absence"]]


def _app() -> QtWidgets.QApplication:
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from annotation_explorer_polars.theme import apply_application_theme

    apply_application_theme(app)
    return app


def _spin(app: QtWidgets.QApplication, seconds: float = 0.05) -> None:
    end = time.monotonic() + seconds
    while True:
        app.processEvents(QtCore.QEventLoop.ProcessEventsFlag.AllEvents, 20)
        if time.monotonic() >= end:
            break
        time.sleep(0.005)


def _sorted_groups(groups: List[List[str]]) -> List[List[str]]:
    return sorted(sorted(group) for group in groups)


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class DockHostTests(unittest.TestCase):
    UIDS = ("tree", "table", "track", "map", "pan")

    @classmethod
    def setUpClass(cls) -> None:
        cls.app = _app()

    def setUp(self) -> None:
        self.window = QtWidgets.QMainWindow()
        self.splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        self.rail = QtWidgets.QLabel("rail")
        self.rail.setObjectName("FilterRail")
        self.rail.setMinimumWidth(200)
        self.splitter.addWidget(self.rail)
        self.host = DockHost()
        self.splitter.addWidget(self.host)
        self.window.setCentralWidget(self.splitter)
        self.frames = {}
        for uid in self.UIDS:
            self._register(self.host, uid)
        self.host.place_view("tree", Placement("left", outer=True))
        self.host.place_view("table", Placement("right", relative_to="tree"))
        self.host.place_view("track", Placement("bottom", relative_to="table"))
        self.host.place_view("map", Placement("bottom", outer=True))
        self.host.place_view("pan", Placement("tab", relative_to="map"))
        self.window.resize(1200, 800)
        self.window.show()
        _spin(self.app, 0.1)

    def tearDown(self) -> None:
        self.window.close()
        self.window.deleteLater()
        _spin(self.app, 0.02)

    def _register(self, host: "DockHost", uid: str) -> "PanelFrame":
        frame = PanelFrame(PanelSpec(uid, uid.title(), "Test", placeholder=uid))
        host.register_view(uid, uid.title(), frame)
        self.frames.setdefault(id(host), {})[uid] = frame
        return frame

    def _groups(self, host: Optional["DockHost"] = None) -> List[List[str]]:
        return _sorted_groups((host or self.host).groups())

    def test_add_tab_and_split(self) -> None:
        self.assertEqual(self._groups(), [["map", "pan"], ["table"], ["track"], ["tree"]])
        self.assertEqual(self.host.group_uids("pan"), ["map", "pan"])
        self.assertTrue(self.host.move_view("pan", Placement("right", relative_to="map")))
        _spin(self.app)
        self.assertEqual(self._groups(), [["map"], ["pan"], ["table"], ["track"], ["tree"]])
        self.assertTrue(self.host.move_view("track", Placement("tab", relative_to="tree")))
        self.assertEqual(self.host.group_uids("tree"), ["tree", "track"])

    def test_lone_panel_cannot_split_its_own_group(self) -> None:
        self.assertFalse(self.host.move_view("tree", Placement("right", relative_to="tree")))
        self.assertEqual(self._groups(), [["map", "pan"], ["table"], ["track"], ["tree"]])

    def test_close_and_reopen_in_place(self) -> None:
        closed: List[str] = []
        self.host.viewClosed.connect(closed.append)
        neighbours = self.host.group_uids("map")
        self.host.close_view("map")
        _spin(self.app)
        self.assertEqual(closed, ["map"])
        self.assertFalse(self.host.is_open("map"))
        self.assertTrue(self.host.is_placed("map"))
        self.assertTrue(self.host.open_view("map"))
        _spin(self.app)
        self.assertTrue(self.host.is_open("map"))
        self.assertEqual(self.host.group_uids("map"), neighbours)

    def test_save_restore_round_trip_keeps_closed_panels_closed(self) -> None:
        self.host.move_view("pan", Placement("top", relative_to="tree"))
        self.host.close_view("track")
        _spin(self.app)
        state = self.host.save_layout()
        before = self._groups()
        self.host.move_view("pan", Placement("tab", relative_to="map"))
        self.host.open_view("track")
        _spin(self.app)
        self.assertTrue(self.host.restore_layout(state))
        _spin(self.app)
        self.assertEqual(self._groups(), before)
        self.assertFalse(self.host.is_open("track"))
        self.assertTrue(self.host.open_view("track"))
        self.assertEqual(self.host.floating_count(), 0)

    def test_restore_into_a_fresh_host(self) -> None:
        # A newly built window has no docks in the layout yet; restore must
        # still bring every saved panel back (regression found in Phase 1).
        state = self.host.save_layout()
        other = DockHost()
        for uid in self.UIDS:
            self._register(other, uid)
        other.resize(1000, 700)
        other.show()
        try:
            self.assertTrue(other.restore_layout(state))
            _spin(self.app)
            self.assertEqual(self._groups(other), self._groups())
        finally:
            other.close()
            other.deleteLater()

    def test_unknown_uid_in_state_is_dropped(self) -> None:
        self._register(self.host, "ghost")
        self.host.place_view("ghost", Placement("tab", relative_to="tree"))
        state = self.host.save_layout()
        other = DockHost()
        for uid in self.UIDS:
            self._register(other, uid)
        other.show()
        try:
            self.assertTrue(other.restore_layout(state))
            _spin(self.app)
            self.assertNotIn("ghost", other.view_uids())
            self.assertEqual(self._groups(other), [["map", "pan"], ["table"], ["track"], ["tree"]])
        finally:
            other.close()
            other.deleteLater()

    def test_corrupt_state_returns_false_and_keeps_layout(self) -> None:
        before = self._groups()
        self.assertFalse(self.host.restore_layout(b"not a dock state"))
        self.assertFalse(self.host.restore_layout(b""))
        self.assertEqual(self._groups(), before)

    def test_every_panel_is_not_floatable(self) -> None:
        for uid in self.host.view_uids():
            self.assertFalse(self.host.view_is_floatable(uid), uid)

    def test_panel_missing_from_a_restored_state_opens_docked(self) -> None:
        # ADS floats a dock it does not know about when toggled open; DockHost
        # must place it instead.
        self.host.close_view("pan")
        state = self.host.save_layout()
        other = DockHost()
        for uid in ("tree", "table", "track", "map"):
            self._register(other, uid)
        other.show()
        try:
            self.assertTrue(other.restore_layout(state))
            self._register(other, "late")
            self.assertTrue(other.open_view("late"))
            _spin(self.app)
            self.assertTrue(other.is_open("late"))
            self.assertEqual(other.floating_count(), 0)
        finally:
            other.close()
            other.deleteLater()

    def test_floating_container_is_redocked(self) -> None:
        # Force ADS to create a floating container through the one path found
        # in Phase 1 (toggling open a dock that a restored state left
        # unassigned) and check the guard puts it back into the window.
        self.host.close_view("pan")
        state = self.host.save_layout()
        self._register(self.host, "stray")
        self.host.place_view("stray", Placement("tab", relative_to="tree"))
        self.assertTrue(self.host.restore_layout(state))
        redocked: List[str] = []
        self.host.floatingRedocked.connect(redocked.append)
        self.host._docks["stray"].toggleView(True)  # bypasses DockHost on purpose
        _spin(self.app, 0.2)
        self.assertEqual(redocked, ["stray"])
        self.assertEqual(self.host.floating_count(), 0)
        self.assertTrue(self.host.is_open("stray"))
        self.assertTrue(self.window.isAncestorOf(self.host._docks["stray"].dockAreaWidget()))

    def test_scripted_drags_never_create_a_floating_container(self) -> None:
        created: List[object] = []
        self.host._manager.floatingWidgetCreated.connect(created.append)
        before = self._groups()
        for uid, target in (
            ("map", QtCore.QPoint(3000, 3000)),  # far outside the window
            ("tree", self.rail.mapToGlobal(self.rail.rect().center())),  # over the rail
        ):
            tab = self.host._docks[uid].tabWidget()
            start = QtCore.QPoint(8, tab.height() // 2)
            QTest.mousePress(tab, QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.KeyboardModifier.NoModifier, start)
            end = tab.mapFromGlobal(target)
            for step in range(1, 21):
                point = start + (end - start) * step / 20
                QTest.mouseMove(tab, QtCore.QPoint(int(point.x()), int(point.y())))
                _spin(self.app, 0.005)
            QTest.mouseRelease(tab, QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.KeyboardModifier.NoModifier, end)
            _spin(self.app, 0.05)
            self.assertEqual(self.host.floating_count(), 0, uid)
        self.assertEqual(created, [])
        self.assertEqual(self._groups(), before)
        for uid in self.UIDS:
            # ADS parks background tabs without a parent, so check the group.
            area = self.host._docks[uid].dockAreaWidget()
            self.assertTrue(self.window.isAncestorOf(area), uid)

    def test_double_click_on_a_tab_maximizes_and_never_floats(self) -> None:
        tab = self.host._docks["table"].tabWidget()
        QTest.mouseDClick(tab, QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.KeyboardModifier.NoModifier, QtCore.QPoint(8, 5))
        _spin(self.app)
        self.assertTrue(self.host.is_maximized())
        self.assertEqual(self.host.floating_count(), 0)
        frames = self.frames[id(self.host)]
        bottom = self.host.current_in_group("map")
        self.assertTrue(frames["table"].isVisible())
        self.assertFalse(frames["tree"].isVisible())
        self.assertFalse(frames[bottom].isVisible())
        self.assertTrue(self.rail.isVisible())
        QTest.mouseDClick(tab, QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.KeyboardModifier.NoModifier, QtCore.QPoint(8, 5))
        _spin(self.app)
        self.assertFalse(self.host.is_maximized())
        self.assertTrue(frames["tree"].isVisible())
        self.assertTrue(frames[bottom].isVisible())

    def test_save_while_maximized_saves_the_full_layout(self) -> None:
        before = self._groups()
        self.assertTrue(self.host.maximize_group("table"))
        state = self.host.save_layout()
        self.assertFalse(self.host.is_maximized())
        self.assertTrue(self.host.restore_layout(state))
        _spin(self.app)
        self.assertEqual(self._groups(), before)
        self.assertTrue(self.frames[id(self.host)]["tree"].isVisible())

    def test_ads_context_menus_are_replaced(self) -> None:
        requested: List[str] = []
        self.host.tabMenuRequested.connect(lambda uid, _pos: requested.append(f"tab:{uid}"))
        self.host.groupMenuRequested.connect(lambda uid, _pos: requested.append(f"group:{uid}"))
        tab = self.host._docks["map"].tabWidget()
        event = QtGui.QContextMenuEvent(QtGui.QContextMenuEvent.Reason.Mouse, QtCore.QPoint(5, 5), tab.mapToGlobal(QtCore.QPoint(5, 5)))
        QtWidgets.QApplication.sendEvent(tab, event)
        title_bar = self.host._docks["map"].dockAreaWidget().titleBar()
        event = QtGui.QContextMenuEvent(QtGui.QContextMenuEvent.Reason.Mouse, QtCore.QPoint(2, 2), title_bar.mapToGlobal(QtCore.QPoint(2, 2)))
        QtWidgets.QApplication.sendEvent(title_bar, event)
        _spin(self.app)
        self.assertEqual(requested, ["tab:map", f"group:{self.host.current_in_group('map')}"])
        self.assertIsNone(QtWidgets.QApplication.activePopupWidget())

    def test_group_strip_has_add_and_menu_buttons(self) -> None:
        # The buttons are Python-made children of ADS-made title bars; their
        # connections must survive garbage collection of the title bar's
        # temporary wrapper (regression found in Phase 1).
        gc.collect()
        _spin(self.app)
        title_bar = self.host._docks["map"].dockAreaWidget().titleBar()
        add = title_bar.findChild(QtWidgets.QToolButton, "DockGroupAddButton")
        menu = title_bar.findChild(QtWidgets.QToolButton, "DockGroupMenuButton")
        self.assertIsNotNone(add)
        self.assertIsNotNone(menu)
        self.assertTrue(sip.ispycreated(add))
        asked: List[str] = []
        menus: List[str] = []
        self.host.addToGroupRequested.connect(asked.append)
        self.host.groupMenuRequested.connect(lambda uid, _pos: menus.append(uid))
        add.click()
        menu.click()
        current = self.host.current_in_group("map")
        self.assertEqual(asked, [current])
        self.assertEqual(menus, [current])

    def test_filter_handles_generically_typed_wrappers(self) -> None:
        # A wrapper first made through a generic path (a Python style's
        # polish, an application event filter) is typed as a plain Qt class;
        # the tab and title bar must still be recognised.
        requested: List[str] = []
        self.host.tabMenuRequested.connect(lambda uid, _pos: requested.append(f"tab:{uid}"))
        self.host.groupMenuRequested.connect(lambda uid, _pos: requested.append(f"group:{uid}"))
        tab = sip.cast(self.host._docks["track"].tabWidget(), QtWidgets.QFrame)
        event = QtGui.QContextMenuEvent(QtGui.QContextMenuEvent.Reason.Mouse, QtCore.QPoint(2, 2), QtCore.QPoint(2, 2))
        self.assertTrue(self.host.eventFilter(tab, event))
        title_bar = sip.cast(self.host._docks["track"].dockAreaWidget().titleBar(), QtWidgets.QFrame)
        self.assertTrue(self.host.eventFilter(title_bar, event))
        self.assertEqual(requested, ["tab:track", "group:track"])

    def test_neighbor_group(self) -> None:
        self.assertEqual(self.host.neighbor_group("tree", "right"), "table")
        self.assertIn(self.host.neighbor_group("tree", "bottom"), ("map", "pan"))
        self.assertEqual(self.host.neighbor_group("track", "top"), "table")
        self.assertIsNone(self.host.neighbor_group("tree", "left"))

    def test_empty_state_page(self) -> None:
        page = QtWidgets.QLabel("empty")
        self.host.set_empty_page(page)
        changes: List[bool] = []
        self.host.emptyChanged.connect(changes.append)
        for uid in self.UIDS:
            self.host.close_view(uid)
        _spin(self.app)
        self.assertTrue(self.host.is_empty())
        self.assertTrue(page.isVisible())
        self.host.open_view("tree")
        _spin(self.app)
        self.assertFalse(self.host.is_empty())
        self.assertFalse(page.isVisible())
        self.assertEqual(changes, [True, False])

    def test_rail_is_not_part_of_the_docking_area(self) -> None:
        self.assertFalse(self.host.isAncestorOf(self.rail))
        self.host.clear_layout()
        self.host.place_view("tree", Placement("left", outer=True))
        _spin(self.app)
        self.assertTrue(self.rail.isVisible())
        self.assertIs(self.rail.parentWidget(), self.splitter)

    def test_global_ads_config(self) -> None:
        # Cheap guard for the provisional flags in ads_config.py.
        manager = ads.CDockManager
        flags = manager.eConfigFlag
        for name in ads_config.CONFIG_FLAGS_ON:
            self.assertTrue(manager.testConfigFlag(getattr(flags, name)), name)
        for name in ads_config.CONFIG_FLAGS_OFF:
            self.assertFalse(manager.testConfigFlag(getattr(flags, name)), name)
        for name in ("DoubleClickUndocksWidget", "DockAreaHasUndockButton"):
            self.assertIn(name, ads_config.CONFIG_FLAGS_OFF)

    def test_close_while_maximized_leaves_no_empty_group(self) -> None:
        # Close a panel that is alone in a group hidden by Maximize: ending
        # the maximize must not show that now empty group again.
        self.assertTrue(self.host.maximize_group("table"))
        self.host.close_view("tree")
        _spin(self.app)
        self.assertFalse(self.host.is_maximized())
        area = self.host._docks["tree"].dockAreaWidget()
        self.assertEqual(area.openDockWidgetsCount(), 0)
        self.assertFalse(area.isVisible())
        frames = self.frames[id(self.host)]
        for uid in ("table", "track"):
            self.assertTrue(frames[uid].isVisible(), uid)
        self.assertTrue(frames[self.host.current_in_group("map")].isVisible())
        # A whole column emptied while hidden stays hidden too.
        self.assertTrue(self.host.maximize_group("tree") is False)
        self.assertTrue(self.host.maximize_group("map"))
        self.host.close_view("table")
        self.host.close_view("track")
        _spin(self.app)
        for uid in ("table", "track"):
            self.assertFalse(self.host._docks[uid].dockAreaWidget().isVisible(), uid)
        self.assertTrue(frames[self.host.current_in_group("map")].isVisible())

    def test_panel_closed_before_a_rebuild_gets_its_tab_back(self) -> None:
        # Reset Layout after a close: ADS kept the tab of the closed dock
        # hidden when the dock was taken out and placed again.
        self.host.close_view("track")
        self.host.close_view("tree")
        _spin(self.app)
        self.host.clear_layout()
        self.host.place_view("tree", Placement("left", outer=True))
        self.host.place_view("track", Placement("right", relative_to="tree"))
        self.host.place_view("table", Placement("bottom", relative_to="track"))
        _spin(self.app)
        for uid in ("tree", "track", "table"):
            self.assertTrue(self.host.is_open(uid), uid)
            self.assertTrue(self.host._docks[uid].tabWidget().isVisible(), uid)

    def test_set_current_tab_moves_the_focus_highlight(self) -> None:
        self.host.focus_view("pan")
        _spin(self.app)
        self.host.set_current_tab("map")
        _spin(self.app)
        self.assertEqual(self.host._manager.focusedDockWidget().objectName(), "map")
        self.assertEqual(self.host.active_view(), "map")
        # Another group's focus is left alone.
        self.host.focus_view("tree")
        self.host.set_current_tab("pan")
        _spin(self.app)
        self.assertEqual(self.host._manager.focusedDockWidget().objectName(), "tree")

    def test_double_click_on_a_divider_evens_out_the_two_sides(self) -> None:
        splitter = self.host._docks["tree"].dockAreaWidget().parentWidget()
        self.assertIsInstance(splitter, QtWidgets.QSplitter)
        total = sum(splitter.sizes())
        splitter.setSizes([total // 4, total - total // 4])
        _spin(self.app)
        handle = splitter.handle(1)
        QTest.mouseDClick(handle, QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.KeyboardModifier.NoModifier, handle.rect().center())
        _spin(self.app)
        left, right = splitter.sizes()
        self.assertLessEqual(abs(left - right), 2)
        self.assertEqual(self.host.floating_count(), 0)

    def test_double_click_on_a_group_title_bar_maximizes(self) -> None:
        title_bar = self.host._docks["track"].dockAreaWidget().titleBar()
        point = QtCore.QPoint(title_bar.width() - 120, title_bar.height() // 2)
        QTest.mouseDClick(title_bar, QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.KeyboardModifier.NoModifier, point)
        _spin(self.app)
        self.assertTrue(self.host.is_maximized())
        self.assertEqual(self.host.maximized_view(), "track")
        self.assertEqual(self.host.floating_count(), 0)
        QTest.mouseDClick(title_bar, QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.KeyboardModifier.NoModifier, point)
        _spin(self.app)
        self.assertFalse(self.host.is_maximized())

    def test_middle_click_closes_a_tab(self) -> None:
        tab = self.host._docks["table"].tabWidget()
        point = QtCore.QPoint(8, tab.height() // 2)
        # ADS closes only while the cursor is still over the tab. Put the
        # cursor back afterwards: ADS drops use the cursor position too.
        old_pos = QtGui.QCursor.pos()
        self.addCleanup(QtGui.QCursor.setPos, old_pos)
        QtGui.QCursor.setPos(tab.mapToGlobal(point))
        QTest.mouseClick(tab, QtCore.Qt.MouseButton.MiddleButton, QtCore.Qt.KeyboardModifier.NoModifier, point)
        _spin(self.app)
        self.assertFalse(self.host.is_open("table"))
        self.assertTrue(self.host.open_view("table"))
        self.assertTrue(self.host.is_open("table"))


def _fake_specs(
    extra: int = 0,
    unavailable: tuple = (),
    requires: Optional[dict] = None,
    multi: tuple = (),
) -> List["PanelSpec"]:
    """Fake panel kinds. As in the real window the track requires the first
    Annotations panel (``requires`` overrides that map); ``multi`` lists the
    multi-instance kinds."""
    kinds = ["sample_tree", "annotations", "track", "geo", "phylogeny", "gene_presence_absence", "pan_metadata", "cgmlst_force"]
    kinds += [f"extra{i}" for i in range(extra)]
    requires = {"track": "annotations"} if requires is None else requires
    specs = []
    for kind in kinds:
        available = (lambda: (False, "Not available: test reason")) if kind in unavailable else (lambda: (True, ""))
        specs.append(
            PanelSpec(
                kind,
                kind.replace("_", " ").title(),
                "Test",
                placeholder=kind,
                available=available,
                single_instance=kind not in multi,
                requires=requires.get(kind, ""),
            )
        )
    return specs


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class PanelManagerTests(unittest.TestCase):
    """PanelManager on fake panels: placement rules, menus, palette, closed stack."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.app = _app()

    def _make(self, **kwargs: Any) -> "PanelManager":
        window = QtWidgets.QMainWindow()
        manager = PanelManager(window, _fake_specs(**kwargs))
        splitter = QtWidgets.QSplitter()
        splitter.addWidget(QtWidgets.QLabel("rail"))
        splitter.addWidget(manager.host)
        window.setCentralWidget(splitter)
        window.resize(1300, 850)
        window.show()
        manager.build_default()
        _spin(self.app, 0.1)
        self.addCleanup(window.deleteLater)
        self.addCleanup(window.close)
        return manager

    @staticmethod
    def _action(menu: QtWidgets.QMenu, path: str) -> QtGui.QAction:
        parts = path.split(" > ")
        for part in parts[:-1]:
            menu = next(a.menu() for a in menu.actions() if a.text() == part)
        return next(a for a in menu.actions() if a.text() == parts[-1])

    def test_new_panel_opens_as_a_tab_of_the_active_group(self) -> None:
        manager = self._make()
        manager.host.focus_view("annotations")
        self.assertTrue(manager.open_panel("pan_metadata"))
        _spin(self.app)
        self.assertEqual(manager.host.group_uids("annotations"), ["annotations", "pan_metadata"])
        self.assertTrue(manager.host.is_current_tab("pan_metadata"))

    def test_single_instance_kind_focuses_the_existing_panel(self) -> None:
        manager = self._make()
        groups = _sorted_groups(manager.host.groups())
        manager.host.set_current_tab("geo")
        manager.host.focus_view("sample_tree")
        self.assertTrue(manager.open_panel("phylogeny"))
        _spin(self.app)
        self.assertEqual(_sorted_groups(manager.host.groups()), groups)
        self.assertTrue(manager.host.is_current_tab("phylogeny"))
        self.assertEqual(manager.host.active_view(), "phylogeny")

    def test_unavailable_kinds_are_disabled_and_left_out_of_classic(self) -> None:
        from annotation_explorer_polars.panels.palette import PanelPalette

        manager = self._make(unavailable=("phylogeny", "cgmlst_force"))
        open_uids = set(manager.host.open_uids())
        self.assertNotIn("phylogeny", open_uids)
        self.assertEqual(open_uids, {"sample_tree", "annotations", "track", "geo", "gene_presence_absence"})
        self.assertFalse(manager.open_panel("phylogeny"))
        self.assertFalse(manager.open_panel("cgmlst_force", to_side=True))
        self.assertFalse(manager.is_open("phylogeny"))
        palette = PanelPalette(manager.specs(), is_open=manager.is_open)
        try:
            by_kind = {item.data(0, QtCore.Qt.ItemDataRole.UserRole + 1): item for item in palette.items()}
            for kind in ("phylogeny", "cgmlst_force"):
                self.assertFalse(by_kind[kind].flags() & QtCore.Qt.ItemFlag.ItemIsEnabled, kind)
                self.assertEqual(by_kind[kind].toolTip(0), "Not available: test reason")
            self.assertTrue(by_kind["geo"].flags() & QtCore.Qt.ItemFlag.ItemIsEnabled)
        finally:
            palette.deleteLater()
        menu_bar = QtWidgets.QMenuBar()
        self.addCleanup(menu_bar.deleteLater)
        manager.install_menus(menu_bar)
        manager._fill_panels_menu()
        action = self._action(manager._panels_menu, "Phylogeny")
        self.assertFalse(action.isEnabled())
        self.assertEqual(action.toolTip(), "Not available: test reason")

    def test_closed_stack_keeps_the_last_ten_in_order(self) -> None:
        manager = self._make(extra=6)
        kinds = [spec.kind for spec in manager.specs()]
        for kind in kinds:
            manager.open_panel(kind)
        _spin(self.app)
        for kind in kinds:
            manager.close_panel(kind)
        _spin(self.app)
        self.assertEqual(manager.closed_stack(), kinds[-CLOSED_STACK_LIMIT:])
        # Reopened panels leave the stack; one reopened by other means is skipped.
        manager.open_panel(kinds[-2])
        self.assertEqual(manager.reopen_closed(), kinds[-1])
        self.assertEqual(manager.reopen_closed(), kinds[-3])
        _spin(self.app)
        self.assertTrue(manager.is_open(kinds[-1]))
        self.assertTrue(manager.is_open(kinds[-3]))
        self.assertEqual(manager.closed_stack(), kinds[-CLOSED_STACK_LIMIT:-3])

    def test_palette_from_a_group_opens_to_the_side_of_that_group(self) -> None:
        manager = self._make()
        manager.host.focus_view("annotations")

        class FakePalette(QtCore.QObject):
            opened = QtCore.pyqtSignal(str, bool)

            def __init__(self, *_args: Any, **_kwargs: Any) -> None:
                super().__init__()

            def exec(self) -> int:
                self.opened.emit("cgmlst_force", True)
                return 1

        original = manager_module.PanelPalette
        manager_module.PanelPalette = FakePalette
        try:
            manager.show_palette(target="sample_tree")
        finally:
            manager_module.PanelPalette = original
        _spin(self.app)
        self.assertEqual(manager.host.group_uids("cgmlst_force"), ["cgmlst_force"])
        self.assertEqual(manager.host.neighbor_group("sample_tree", "right"), "cgmlst_force")
        self.assertNotEqual(manager.host.neighbor_group("annotations", "right"), "cgmlst_force")

    def test_reset_layout_restores_default_titles(self) -> None:
        manager = self._make()
        manager.rename_panel("geo", "World map")
        self.assertEqual(manager.host.view_title("geo"), "World map")
        manager.reset_layout()
        _spin(self.app)
        self.assertEqual(manager.host.view_title("geo"), "Geo")
        self.assertNotIn("title", next(p for p in manager.panel_list() if p["kind"] == "geo")["params"])

    def test_reset_after_closing_panels_shows_every_tab(self) -> None:
        manager = self._make()
        for uid in list(manager.host.open_uids()):
            manager.close_panel(uid)
        _spin(self.app)
        manager.reset_layout()
        _spin(self.app)
        for uid in manager.host.open_uids():
            self.assertTrue(manager.host._docks[uid].tabWidget().isVisible(), uid)
        manager.close_panel("track")
        manager.reset_layout()
        _spin(self.app)
        self.assertTrue(manager.host._docks["track"].tabWidget().isVisible())

    def test_reset_moves_the_focus_highlight_to_a_visible_tab(self) -> None:
        manager = self._make()
        host = manager.host
        host.focus_view("phylogeny")
        _spin(self.app)
        for uid in list(host.open_uids()):
            manager.close_panel(uid)
        _spin(self.app)
        manager.reset_layout()
        _spin(self.app)
        focused = host._manager.focusedDockWidget()
        if focused is not None:
            self.assertTrue(host.is_current_tab(focused.objectName()), focused.objectName())
        self.assertTrue(host.is_current_tab("geo"))

    def test_tab_menu_actions(self) -> None:
        manager = self._make()
        host = manager.host
        # Alone in its group: splitting, New Group and Close Others are off.
        menu = manager.build_tab_menu("annotations")
        for path in ("Split Right", "Split Down", "Move to Group > New Group", "Close Others"):
            self.assertFalse(self._action(menu, path).isEnabled(), path)
        self.assertTrue(self._action(menu, "Close All in Group").isEnabled())
        # Move to Group > Below: into the nearest group below (the track).
        self._action(menu, "Move to Group > Below").trigger()
        _spin(self.app)
        self.assertEqual(host.group_uids("annotations"), ["track", "annotations"])
        # Now not alone: Split Right makes a new group to the right.
        menu = manager.build_tab_menu("annotations")
        self.assertTrue(self._action(menu, "Split Right").isEnabled())
        self._action(menu, "Split Right").trigger()
        _spin(self.app)
        self.assertEqual(host.group_uids("annotations"), ["annotations"])
        self.assertEqual(host.neighbor_group("track", "right"), "annotations")
        # Split Down from a two-tab group.
        host.move_view("phylogeny", Placement("tab", relative_to="geo"))
        _spin(self.app)
        self._action(manager.build_tab_menu("phylogeny"), "Split Down").trigger()
        _spin(self.app)
        self.assertEqual(host.group_uids("phylogeny"), ["phylogeny"])
        self.assertEqual(host.neighbor_group("geo", "bottom"), "phylogeny")
        # New Group splits right.
        self._action(manager.build_tab_menu("phylogeny"), "Move to Group > Above").trigger()
        _spin(self.app)
        self.assertEqual(host.group_uids("phylogeny"), ["geo", "phylogeny"])
        self._action(manager.build_tab_menu("phylogeny"), "Move to Group > New Group").trigger()
        _spin(self.app)
        self.assertEqual(host.neighbor_group("geo", "right"), "phylogeny")
        # Close Others and Close All in Group.
        host.move_view("phylogeny", Placement("tab", relative_to="geo"))
        manager.open_panel("pan_metadata", placement=Placement("tab", relative_to="geo"))
        _spin(self.app)
        self._action(manager.build_tab_menu("geo"), "Close Others").trigger()
        _spin(self.app)
        self.assertEqual(host.group_uids("geo"), ["geo"])
        self.assertFalse(manager.is_open("phylogeny") or manager.is_open("pan_metadata"))
        self._action(manager.build_tab_menu("geo"), "Close All in Group").trigger()
        _spin(self.app)
        self.assertFalse(manager.is_open("geo"))
        self._action(manager.build_tab_menu("track"), "Close").trigger()
        _spin(self.app)
        self.assertFalse(manager.is_open("track"))
        self.assertEqual(host.floating_count(), 0)

    def test_group_menu_actions(self) -> None:
        manager = self._make()
        host = manager.host
        menu = manager.build_group_menu("sample_tree")
        self.assertFalse(self._action(menu, "Split Group Right").isEnabled())
        self.assertFalse(self._action(menu, "Reopen Closed Panel").isEnabled())
        # Split Group Right / Down move the visible panel out of a two-tab group.
        host.move_view("phylogeny", Placement("tab", relative_to="geo"))
        host.set_current_tab("phylogeny")
        self._action(manager.build_group_menu("phylogeny"), "Split Group Right").trigger()
        _spin(self.app)
        self.assertEqual(host.neighbor_group("geo", "right"), "phylogeny")
        host.move_view("phylogeny", Placement("tab", relative_to="geo"))
        self._action(manager.build_group_menu("phylogeny"), "Split Group Down").trigger()
        _spin(self.app)
        self.assertEqual(host.neighbor_group("geo", "bottom"), "phylogeny")
        # Maximize Group, then Restore Group.
        self._action(manager.build_group_menu("annotations"), "Maximize Group").trigger()
        self.assertTrue(host.is_maximized())
        self._action(manager.build_group_menu("annotations"), "Restore Group").trigger()
        self.assertFalse(host.is_maximized())
        # Close Group closes every tab of the group; Reopen Closed Panel brings
        # the last one back in place.
        host.move_view("phylogeny", Placement("tab", relative_to="geo"))
        _spin(self.app)
        self._action(manager.build_group_menu("geo"), "Close Group").trigger()
        _spin(self.app)
        self.assertFalse(manager.is_open("geo") or manager.is_open("phylogeny"))
        reopen = self._action(manager.build_group_menu("annotations"), "Reopen Closed Panel")
        self.assertTrue(reopen.isEnabled())
        reopen.trigger()
        _spin(self.app)
        reopened = [k for k in ("geo", "phylogeny") if manager.is_open(k)]
        self.assertEqual(len(reopened), 1)
        self.assertEqual(host.neighbor_group("track", "bottom"), reopened[0])

    # ------------------------------------------- the track requires Annotations
    def _close_quietly(self, manager: "PanelManager", uid: str) -> None:
        """Close a dock without the manager's rules (to build a state that
        has the track without Annotations, as an older saved layout may)."""
        manager.host._quiet += 1
        try:
            manager.host._docks[uid].toggleView(False)
        finally:
            manager.host._quiet -= 1

    def test_closing_annotations_closes_the_track_by_every_route(self) -> None:
        manager = self._make()
        menu_bar = QtWidgets.QMenuBar()
        self.addCleanup(menu_bar.deleteLater)
        manager.install_menus(menu_bar)

        def close_by_panels_menu() -> None:
            manager._fill_panels_menu()
            self._action(manager._panels_menu, "Annotations").trigger()

        routes = {
            "close_panel": lambda: manager.close_panel("annotations"),
            "tab menu Close": lambda: self._action(manager.build_tab_menu("annotations"), "Close").trigger(),
            "tab menu Close All in Group": lambda: self._action(
                manager.build_tab_menu("annotations"), "Close All in Group"
            ).trigger(),
            "group menu Close Group": lambda: self._action(
                manager.build_group_menu("annotations"), "Close Group"
            ).trigger(),
            "Panels menu": close_by_panels_menu,
            "dock close": lambda: manager.host._docks["annotations"].toggleView(False),
        }
        for name, close in routes.items():
            manager.reset_layout()
            _spin(self.app)
            self.assertTrue(manager.is_open("track"), name)
            close()
            _spin(self.app)
            self.assertFalse(manager.is_open("annotations"), name)
            self.assertFalse(manager.is_open("track"), name)
            self.assertTrue(manager.is_open("sample_tree"), name)
            self.assertEqual(manager.closed_stack()[-2:], ["annotations", "track"], name)
        self.assertEqual(manager.host.floating_count(), 0)

    def test_closing_the_track_leaves_annotations_open(self) -> None:
        manager = self._make()
        manager.close_panel("track")
        _spin(self.app)
        self.assertFalse(manager.is_open("track"))
        self.assertTrue(manager.is_open("annotations"))

    def test_close_others_from_the_track_keeps_annotations(self) -> None:
        manager = self._make()
        host = manager.host
        host.move_view("track", Placement("tab", relative_to="annotations"))
        manager.open_panel("pan_metadata", placement=Placement("tab", relative_to="annotations"))
        _spin(self.app)
        self._action(manager.build_tab_menu("track"), "Close Others").trigger()
        _spin(self.app)
        self.assertTrue(manager.is_open("track") and manager.is_open("annotations"))
        self.assertFalse(manager.is_open("pan_metadata"))
        # From Annotations' tab Close Others closes the track.
        self._action(manager.build_tab_menu("annotations"), "Close Others").trigger()
        _spin(self.app)
        self.assertTrue(manager.is_open("annotations"))
        self.assertFalse(manager.is_open("track"))

    def test_opening_the_track_opens_annotations_too(self) -> None:
        manager = self._make()
        classic = _sorted_groups(manager.host.groups())
        manager.close_panel("annotations")
        _spin(self.app)
        # Both were closed together: both reopen where they were.
        self.assertTrue(manager.open_panel("track"))
        _spin(self.app)
        self.assertTrue(manager.is_open("annotations") and manager.is_open("track"))
        self.assertEqual(_sorted_groups(manager.host.groups()), classic)
        self.assertEqual(manager.closed_stack(), [])
        # A track never placed (a preset without it): Annotations reopens and
        # the track goes below it.
        manager.apply_preset("geo_table")
        _spin(self.app)
        self.assertFalse(manager.host.is_placed("track"))
        manager.close_panel("annotations")
        _spin(self.app)
        self.assertTrue(manager.open_panel("track"))
        _spin(self.app)
        self.assertTrue(manager.is_open("annotations"))
        self.assertEqual(manager.host.neighbor_group("track", "top"), "annotations")
        self.assertEqual(manager.host.group_uids("track"), ["track"])
        self.assertEqual(manager.host.floating_count(), 0)

    def test_opening_the_track_beside_open_annotations_goes_below_it(self) -> None:
        manager = self._make()
        manager.apply_preset("geo_table")
        _spin(self.app)
        manager.host.focus_view("geo")
        self.assertTrue(manager.open_panel("track"))
        _spin(self.app)
        self.assertEqual(manager.host.neighbor_group("track", "top"), "annotations")
        # An explicit placement is honored while Annotations is open.
        manager.close_panel("track")
        self.assertTrue(manager.open_panel("track", placement=Placement("tab", relative_to="geo")))
        _spin(self.app)
        self.assertEqual(manager.host.group_uids("track"), ["geo", "track"])

    def test_palette_opens_the_track_with_annotations(self) -> None:
        manager = self._make()
        manager.close_panel("annotations")
        manager.apply_preset("geo_table")  # neither is placed with the track now
        manager.close_panel("annotations")
        _spin(self.app)
        self.assertFalse(manager.is_open("annotations") or manager.is_open("track"))

        class FakePalette(QtCore.QObject):
            opened = QtCore.pyqtSignal(str, bool)

            def __init__(self, *_args: Any, **_kwargs: Any) -> None:
                super().__init__()

            def exec(self) -> int:
                self.opened.emit("track", True)
                return 1

        original = manager_module.PanelPalette
        manager_module.PanelPalette = FakePalette
        try:
            manager.show_palette(target="geo")
        finally:
            manager_module.PanelPalette = original
        _spin(self.app)
        self.assertTrue(manager.is_open("annotations") and manager.is_open("track"))
        self.assertEqual(manager.host.neighbor_group("geo", "right"), "annotations")
        self.assertEqual(manager.host.neighbor_group("track", "top"), "annotations")
        self.assertEqual(manager.host.floating_count(), 0)

    def test_palette_and_panels_menu_say_annotations_opens_too(self) -> None:
        from annotation_explorer_polars.panels.palette import PanelPalette

        manager = self._make()
        menu_bar = QtWidgets.QMenuBar()
        self.addCleanup(menu_bar.deleteLater)
        manager.install_menus(menu_bar)
        manager._fill_panels_menu()
        self.assertEqual(
            self._action(manager._panels_menu, "Track").toolTip(),
            "Opens with Annotations; closing Annotations closes it",
        )
        manager.close_panel("annotations")
        palette = PanelPalette(manager.specs(), is_open=manager.is_open, available=manager.is_available)
        try:
            item = palette.item_for("track")
            self.assertEqual(item.toolTip(0), "Open Track (opens Annotations too)")
        finally:
            palette.deleteLater()

    def test_reopen_closed_brings_annotations_and_the_track_back(self) -> None:
        manager = self._make()
        manager.close_panel("annotations")
        _spin(self.app)
        self.assertEqual(manager.reopen_closed(), "track")
        _spin(self.app)
        self.assertTrue(manager.is_open("annotations") and manager.is_open("track"))
        self.assertEqual(manager._reopenable(), [])

    def test_a_restored_track_without_annotations_is_closed(self) -> None:
        manager = self._make()
        self._close_quietly(manager, "annotations")
        self.assertTrue(manager.is_open("track"))
        state, panels = manager.host.save_layout(), manager.panel_list()
        manager.reset_layout()
        _spin(self.app)
        self.assertTrue(manager.restore_layout(state, panels))
        _spin(self.app)
        self.assertFalse(manager.is_open("track"))
        self.assertTrue(manager.is_open("sample_tree"))
        self.assertNotIn("track", manager.closed_stack())
        manager.reset_layout()
        _spin(self.app)
        self.assertTrue(manager.apply_saved_layout(state, panels))
        _spin(self.app)
        self.assertFalse(manager.is_open("track"))
        self.assertEqual(manager.closed_stack(), [])

    def test_a_saved_layout_with_only_the_track_changes_nothing(self) -> None:
        manager = self._make()
        for uid in list(manager.host.open_uids()):
            if uid != "track":
                self._close_quietly(manager, uid)
        self.assertEqual(manager.host.open_uids(), ["track"])
        state, panels = manager.host.save_layout(), manager.panel_list()
        manager.reset_layout()
        _spin(self.app)
        before = _sorted_groups(manager.host.groups())
        self.assertFalse(manager.apply_saved_layout(state, panels))
        _spin(self.app)
        self.assertEqual(_sorted_groups(manager.host.groups()), before)
        self.assertTrue(manager.is_open("track"))

    def test_restoring_a_layout_with_only_the_track_reports_failure(self) -> None:
        # Startup and session restores call restore_layout directly: a layout
        # left empty by the requirement rule must say so, so the window builds
        # Classic instead of opening an empty docking area.
        manager = self._make()
        for uid in list(manager.host.open_uids()):
            if uid != "track":
                self._close_quietly(manager, uid)
        state, panels = manager.host.save_layout(), manager.panel_list()
        manager.host.clear_layout()
        _spin(self.app)
        self.assertFalse(manager.restore_layout(state, panels))
        self.assertEqual(manager.host.open_uids(), [])
        manager.build_default()
        _spin(self.app)
        self.assertTrue(manager.is_open("annotations") and manager.is_open("track"))

    def test_close_others_is_disabled_when_only_the_required_panel_is_beside(self) -> None:
        manager = self._make()
        manager.host.move_view("track", Placement("tab", relative_to="annotations"))
        _spin(self.app)
        self.assertEqual(sorted(manager.host.group_uids("track")), ["annotations", "track"])
        self.assertFalse(self._action(manager.build_tab_menu("track"), "Close Others").isEnabled())
        # From Annotations the track is another panel, so it can be closed.
        self.assertTrue(self._action(manager.build_tab_menu("annotations"), "Close Others").isEnabled())

    def test_reopening_a_squashed_pair_restores_the_classic_split(self) -> None:
        manager = self._make()
        host = manager.host
        _spin(self.app)
        host.set_split_sizes("track", [5, 95])
        _spin(self.app)
        sizes = host.split_sizes("track")
        self.assertLess(sizes[0], 0.2 * sum(sizes))
        manager.close_panel("annotations")
        _spin(self.app)
        self.assertTrue(manager.open_panel("track"))
        _spin(self.app)
        sizes = host.split_sizes("track")
        self.assertEqual(len(sizes), 2)
        self.assertGreater(min(sizes), 0.2 * sum(sizes))

    def test_without_annotations_the_track_is_unavailable(self) -> None:
        manager = self._make(unavailable=("annotations",))
        self.assertFalse(manager.is_open("track"))
        self.assertFalse(manager.is_open("annotations"))
        available, reason = manager.is_available("track")
        self.assertFalse(available)
        self.assertEqual(reason, "Needs Annotations, which is not available for this run")
        self.assertFalse(manager.open_panel("track"))
        self.assertEqual(manager.host.floating_count(), 0)

    def test_a_duplicated_annotations_panel_does_not_keep_the_track(self) -> None:
        manager = self._make(multi=("annotations",))
        extra = manager.duplicate_panel("annotations", "right")
        self.assertIsNotNone(extra)
        _spin(self.app)
        # Closing the extra instance leaves the track (fed by the first one).
        manager.close_panel(extra)
        _spin(self.app)
        self.assertTrue(manager.is_open("track"))
        manager.reopen_closed()
        _spin(self.app)
        self.assertTrue(manager.is_open(extra))
        # Closing the first instance closes the track; the extra stays open.
        manager.close_panel("annotations")
        _spin(self.app)
        self.assertFalse(manager.is_open("track"))
        self.assertTrue(manager.is_open(extra))
        # Opening the track brings the first instance back with it.
        self.assertTrue(manager.open_panel("track"))
        _spin(self.app)
        self.assertTrue(manager.is_open("annotations") and manager.is_open("track"))
        self.assertTrue(manager.is_open(extra))


class _Explorer:
    """The real window on the synthetic run of the golden tests."""

    @staticmethod
    def build(output_dir: str) -> Any:
        from annotation_explorer_polars.windows.explore_window import ExploreWindow

        ann_df, meta_df = build_synthetic_run()
        return ExploreWindow(
            ann_df,
            meta_df,
            ref_df=pd.DataFrame(),
            marker_hits_df=pd.DataFrame(),
            output_dir=Path(output_dir),
        )

    @staticmethod
    def show_and_settle(app: QtWidgets.QApplication, win: Any, seconds: float = 1.2) -> None:
        done = {"v": win.ctx.first_result_done}
        win.firstRefreshFinished.connect(lambda: done.__setitem__("v", True))
        win.resize(1500, 950)
        win.show()
        end = time.monotonic() + 20
        while not done["v"] and time.monotonic() < end:
            _spin(app, 0.02)
        # Past the 700 ms delay before heavy panels may build.
        _spin(app, seconds)


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class ExploreWindowDockingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = _app()
        import annotation_explorer_polars.geo_experimental_views as gev

        cls._saved_webengine_cache = gev._WEBENGINE_RUNTIME_CACHE
        gev._WEBENGINE_RUNTIME_CACHE = (False, "disabled in tests")

    @classmethod
    def tearDownClass(cls) -> None:
        import annotation_explorer_polars.geo_experimental_views as gev

        gev._WEBENGINE_RUNTIME_CACHE = cls._saved_webengine_cache

    def setUp(self) -> None:
        self.output_dir = tempfile.mkdtemp(prefix="bactolomew_dock_run_")
        # A tree and Panaroo output, so Classic has its Phylogeny and Gene
        # presence/absence panels.
        write_pan_outputs(self.output_dir)
        self.windows: List[Any] = []
        self._clear_last_window()

    def tearDown(self) -> None:
        for win in self.windows:
            win.close()
            win.deleteLater()
        _spin(self.app, 0.05)
        shutil.rmtree(self.output_dir, ignore_errors=True)
        self._clear_last_window()

    @staticmethod
    def _clear_last_window() -> None:
        # Global fallback for runs without a saved layout; keep tests independent.
        from annotation_explorer_polars.shell.layout_store import LAST_WINDOW_GROUP

        settings = QtCore.QSettings()
        settings.remove(LAST_WINDOW_GROUP)
        settings.sync()

    def _wait_refresh(self, win: Any) -> None:
        end = time.monotonic() + 10
        while time.monotonic() < end and (
            win.ctx.scheduler.is_running() or win.ctx.scheduler.debounce_timer.isActive() or win._panels.has_pending_delivery()
        ):
            _spin(self.app, 0.02)
        _spin(self.app, 0.3)

    def _window(self, show: bool = True, output_dir: Optional[str] = "") -> Any:
        if output_dir == "":
            output_dir = self.output_dir
        if output_dir is None:
            from annotation_explorer_polars.windows.explore_window import ExploreWindow

            ann_df, meta_df = build_synthetic_run()
            win = ExploreWindow(
                ann_df, meta_df, ref_df=pd.DataFrame(), marker_hits_df=pd.DataFrame(), output_dir=None
            )
        else:
            win = _Explorer.build(output_dir)
        self.windows.append(win)
        if show:
            _Explorer.show_and_settle(self.app, win)
        return win

    def test_classic_layout_and_rail(self) -> None:
        win = self._window()
        manager = win._panel_manager
        self.assertIsNotNone(manager)
        self.assertFalse(hasattr(win, "_legacy_layout"))
        self.assertEqual(_sorted_groups(manager.host.groups()), _sorted_groups(CLASSIC_GROUPS))
        self.assertTrue(manager.host.is_current_tab("geo"))
        rail = win._filters_panel
        self.assertEqual(rail.objectName(), "FilterRail")
        self.assertIs(rail.parentWidget(), win._workspace_splitter)
        self.assertFalse(manager.host.isAncestorOf(rail))
        self.assertIs(win._workspace_splitter.widget(1), manager.host)

    def test_reset_layout_leaves_rail_and_filters_alone(self) -> None:
        win = self._window()
        manager = win._panel_manager
        win.gene_query.setText("lysR")
        width = win._workspace_splitter.sizes()[0]
        manager.close_panel("annotations")
        manager.host.move_view("phylogeny", Placement("right", relative_to="sample_tree"))
        manager.reset_layout()
        _spin(self.app, 0.2)
        self.assertEqual(_sorted_groups(manager.host.groups()), _sorted_groups(CLASSIC_GROUPS))
        self.assertTrue(win._filters_panel.isVisible())
        self.assertEqual(win._workspace_splitter.sizes()[0], width)
        self.assertEqual(win.gene_query.text(), "lysR")
        self.assertEqual(manager.closed_stack(), [])

    def test_heavy_panels_build_lazily_on_first_show(self) -> None:
        win = self._window()
        manager = win._panel_manager
        # The Classic bottom row is shown and built; the compact cgMLST MST is not open.
        self.assertIsNotNone(manager.panel("geo").view)
        self.assertIsNotNone(manager.panel("phylogeny").tab)
        self.assertIsNone(manager.panel("pan_force"))
        self.assertFalse(manager.frame("pan_force").is_built())
        self.assertIsNone(manager.panel("cgmlst_force"))
        manager.open_panel("pan_force")
        _spin(self.app, 0.3)
        self.assertIsNotNone(manager.panel("pan_force").tab)
        self.assertTrue(manager.frame("pan_force").is_built())
        manager.close_panel("pan_force")
        _spin(self.app, 0.2)
        # No Plot menu: the cgMLST MST is a panel like any other.
        self.assertNotIn("Plot", [a.text() for a in win.menuBar().actions()])
        self.assertFalse(hasattr(win, "_show_phandango_style_plot"))
        manager.open_panel("cgmlst_force")
        _spin(self.app, 0.3)
        self.assertIsNotNone(manager.panel("cgmlst_force").plot)
        self.assertTrue(manager.is_visible("cgmlst_force"))
        # The legacy plot dialogs are gone (Phase 10).
        self.assertFalse(hasattr(win, "_phandango_plot_dialog"))
        self.assertFalse(hasattr(win, "_cgmlst_force_plot_dialog"))

    def test_refresh_reaches_every_visible_heavy_panel(self) -> None:
        win = self._window()
        manager = win._panel_manager
        manager.open_panel("phylogeny")
        _spin(self.app, 0.3)
        # Side by side in Classic: both heavy panels visible at once.
        self.assertEqual(manager.host.neighbor_group("geo", "right"), "phylogeny")
        self.assertTrue(manager.is_visible("geo") and manager.is_visible("phylogeny"))
        win.gene_query.setText("blaOXA")
        end = time.monotonic() + 10
        while time.monotonic() < end and (win.ctx.scheduler.is_running() or win.ctx.scheduler.debounce_timer.isActive() or win._panels.has_pending_delivery()):
            _spin(self.app, 0.02)
        _spin(self.app, 0.4)
        self.assertFalse(manager.frame("phylogeny").is_dirty())
        self.assertFalse(manager.frame("geo").is_dirty())
        # Never-opened panels stay dirty until they are shown.
        self.assertTrue(manager.frame("cgmlst_force").is_dirty())
        # A hidden panel (Phylogeny as a background tab) stays dirty and
        # catches up when shown.
        manager.host.move_view("phylogeny", Placement("tab", relative_to="geo"))
        manager.host.set_current_tab("geo")
        _spin(self.app, 0.1)
        self.assertFalse(manager.is_visible("phylogeny"))
        win.gene_query.setText("lysR")
        end = time.monotonic() + 10
        while time.monotonic() < end and (win.ctx.scheduler.is_running() or win.ctx.scheduler.debounce_timer.isActive() or win._panels.has_pending_delivery()):
            _spin(self.app, 0.02)
        _spin(self.app, 0.4)
        self.assertTrue(manager.frame("phylogeny").is_dirty())
        manager.open_panel("phylogeny")
        _spin(self.app, 0.3)
        self.assertFalse(manager.frame("phylogeny").is_dirty())

    def test_the_genome_track_lives_and_dies_with_annotations(self) -> None:
        # Feedback 2026-10-07: the track cannot exist without Annotations.
        win = self._window()
        manager = win._panel_manager
        manager.close_panel("annotations")
        _spin(self.app, 0.2)
        self.assertFalse(manager.is_open("track"))
        # A filter changed while both are closed.
        win.gene_query.setText("blaOXA")
        self._wait_refresh(win)
        # View > Panels > Genome track opens Annotations too, below which the
        # track sits, and the track draws the selected contig.
        manager._fill_panels_menu()
        action = next(a for a in manager._panels_menu.actions() if a.text() == "Genome track")
        self.assertEqual(action.toolTip(), "Opens with Annotations; closing Annotations closes it")
        action.trigger()
        _spin(self.app, 0.3)
        self.assertTrue(manager.is_open("annotations") and manager.is_open("track"))
        self.assertEqual(manager.host.neighbor_group("track", "top"), "annotations")
        win._open_tree_contig_payload({"sample": "S1", "contig": "3", "contig_label": "Contig 3"}, mode="Contig")
        self._wait_refresh(win)
        self.assertTrue(win.table.model().rowCount() > 0)
        self.assertTrue(win.igv_view.scene().items())
        self.assertEqual(manager.host.floating_count(), 0)

    def test_menus_have_no_detach_entries_and_no_shortcuts(self) -> None:
        win = self._window()
        manager = win._panel_manager
        menus = [manager.build_tab_menu("geo"), manager.build_group_menu("geo")]
        view_menu = next(a.menu() for a in win.menuBar().actions() if a.text() == "View")
        manager._fill_panels_menu()
        menus.append(view_menu)

        def all_actions(menu: QtWidgets.QMenu) -> List[QtGui.QAction]:
            out = []
            for action in menu.actions():
                out.append(action)
                if action.menu() is not None:
                    out.extend(all_actions(action.menu()))
            return out

        texts = []
        for menu in menus:
            for action in all_actions(menu):
                texts.append(action.text())
                self.assertTrue(action.shortcut().isEmpty(), action.text())
        # The Window menu's only shortcut is New Panel... (Ctrl+Shift+P).
        window_menu = win._window_menu
        window_menu.aboutToShow.emit()
        shortcuts = {a.text(): a.shortcut().toString() for a in all_actions(window_menu) if not a.shortcut().isEmpty()}
        self.assertEqual(shortcuts, {"New Panel...": "Ctrl+Shift+P"})
        joined = " ".join(texts).lower()
        for banned in ("detach", "float", "undock", "move to explore window", "new explore window"):
            self.assertNotIn(banned, joined)
        for expected in ("Split Right", "Split Down", "Move to Group", "Rename...", "Close", "Close Others",
                         "Close All in Group", "Open Panel Here...", "Reopen Closed Panel",
                         "Split Group Right", "Split Group Down", "Maximize Group", "Close Group",
                         "Reset Layout"):
            self.assertIn(expected, texts)
        # The panel picker moved to Window > New Panel... (one entry point).
        self.assertNotIn("Open Panel...", texts)

    def test_reopen_closed_panel_and_empty_state(self) -> None:
        win = self._window()
        manager = win._panel_manager
        for uid in list(manager.host.open_uids()):
            manager.close_panel(uid)
        _spin(self.app, 0.1)
        self.assertTrue(manager.host.is_empty())
        page = win.findChild(QtWidgets.QWidget, "DockEmptyState")
        self.assertTrue(page.isVisible())
        self.assertTrue(win._filters_panel.isVisible())
        # The last panel closed (the last of the bottom row) comes back.
        self.assertEqual(manager.reopen_closed(), "gene_presence_absence")
        _spin(self.app, 0.1)
        self.assertFalse(manager.host.is_empty())
        page.findChild(QtWidgets.QPushButton, "EmptyStateRestoreLayout")  # exists
        manager.reset_layout()
        _spin(self.app, 0.1)
        self.assertEqual(_sorted_groups(manager.host.groups()), _sorted_groups(CLASSIC_GROUPS))

    def test_palette_lists_every_kind_and_opens_beside(self) -> None:
        from annotation_explorer_polars.panels.palette import PanelPalette

        win = self._window()
        manager = win._panel_manager
        palette = PanelPalette(manager.specs(), is_open=manager.is_open)
        try:
            kinds = {item.data(0, QtCore.Qt.ItemDataRole.UserRole + 1) for item in palette.items()}
            self.assertEqual(kinds, {s.kind for s in manager.specs()})
            palette.search.setText("cgmlst")
            visible = [item for item in palette.items() if not item.isHidden()]
            self.assertEqual([item.data(0, QtCore.Qt.ItemDataRole.UserRole + 1) for item in visible], ["cgmlst_force", "pan_force"])
        finally:
            palette.deleteLater()
        manager.host.focus_view("annotations")
        self.assertTrue(manager.open_panel("cgmlst_force", to_side=True))
        _spin(self.app, 0.2)
        self.assertEqual(manager.host.neighbor_group("annotations", "right"), "cgmlst_force")

    def test_layout_persists_per_run_without_filters(self) -> None:
        win = self._window()
        manager = win._panel_manager
        manager.host.move_view("phylogeny", Placement("left", relative_to="sample_tree"))
        manager.rename_panel("geo", "World map")
        win.annotations_panel.set_view_mode("Context")
        win.gene_query.setText("lysR")
        win._apply_rail_state(333, False)
        _spin(self.app, 0.3)
        saved = _sorted_groups(manager.host.groups())
        win.close()
        _spin(self.app, 0.1)

        again = self._window()
        manager2 = again._panel_manager
        self.assertEqual(_sorted_groups(manager2.host.groups()), saved)
        self.assertEqual(again.view_mode, "Context")
        self.assertEqual(manager2.host.view_title("geo"), "World map")
        self.assertEqual(again.gene_query.text(), "")  # filters start fresh (D5)
        self.assertEqual(again._rail_state()[0], 333)
        self.assertTrue(again._geometry_restored)

    def test_corrupt_saved_layout_falls_back_to_classic(self) -> None:
        from annotation_explorer_polars.shell.layout_store import LayoutStore, WindowLayout, run_key

        LayoutStore().save_window(
            run_key(Path(self.output_dir)),
            WindowLayout(dock_state=b"garbage", panels=[{"uid": "x", "kind": "gone", "params": {}}]),
        )
        win = self._window()
        self.assertEqual(_sorted_groups(win._panel_manager.host.groups()), _sorted_groups(CLASSIC_GROUPS))

    def test_old_format_version_builds_classic(self) -> None:
        from annotation_explorer_polars.shell.layout_store import LayoutStore, WindowLayout, run_key

        win = self._window()
        manager = win._panel_manager
        manager.host.move_view("phylogeny", Placement("left", relative_to="sample_tree"))
        state = manager.host.save_layout()
        win.close()
        LayoutStore().save_window(run_key(Path(self.output_dir)), WindowLayout(dock_state=state, version=0))
        again = self._window()
        self.assertEqual(_sorted_groups(again._panel_manager.host.groups()), _sorted_groups(CLASSIC_GROUPS))

    def test_the_legacy_layout_flag_is_gone(self) -> None:
        # Phase 10: BACTOLOMEW_LEGACY_LAYOUT no longer exists; the docking
        # area is the only layout.
        import annotation_explorer_polars.docking as docking

        self.assertFalse(hasattr(docking, "LEGACY_LAYOUT_ENV"))
        self.assertFalse(hasattr(docking, "legacy_layout_requested"))
        with mock.patch.dict(os.environ, {"BACTOLOMEW_LEGACY_LAYOUT": "1"}):
            win = self._window(show=False)
        self.assertIsNotNone(win._panel_manager)
        self.assertIsNotNone(win.findChild(QtWidgets.QWidget, "DockHost"))
        for name in ("bottom_tabs", "_vertical_splitter", "_top_splitter", "_table_splitter", "view_mode_box"):
            self.assertFalse(hasattr(win, name), name)

    def test_show_filters_collapses_and_restores_the_rail(self) -> None:
        # Window > Show Filters does what the removed top bar toggle did.
        win = self._window()
        host = win._panel_manager.host
        win._apply_rail_state(320, False)
        _spin(self.app, 0.3)
        self.assertEqual(win._workspace_splitter.sizes()[0], 320)
        host_width = host.width()
        self.assertTrue(win.filters_action.isChecked())
        self.assertTrue(win.filters_action.shortcut().isEmpty())
        win.filters_action.trigger()
        _spin(self.app, 0.2)
        self.assertTrue(win._filters_panel.isHidden())
        # Hidden, the rail folds to its strip, which brings it back; the
        # header chevron hides it again, all one state with the action.
        self.assertTrue(win._filters_strip.isVisible())
        self.assertGreater(host.width(), host_width + 250)
        self.assertEqual(_sorted_groups(host.groups()), _sorted_groups(CLASSIC_GROUPS))
        win._filters_strip.expand_button.click()
        _spin(self.app, 0.2)
        self.assertFalse(win._filters_panel.isHidden())
        self.assertTrue(win._filters_strip.isHidden())
        self.assertTrue(win.filters_action.isChecked())
        win._filters_panel.collapse_button.click()
        _spin(self.app, 0.2)
        self.assertTrue(win._filters_panel.isHidden())
        self.assertTrue(win._filters_strip.isVisible())
        self.assertFalse(win.filters_action.isChecked())
        # Saved collapsed and restored collapsed (v2 store, not migration).
        win.close()
        _spin(self.app, 0.1)
        again = self._window()
        self.assertTrue(again._filters_panel.isHidden())
        self.assertFalse(again.filters_action.isChecked())
        again._window_menu.aboutToShow.emit()
        self.assertIn(again.filters_action, again._window_menu.actions())
        again.filters_action.trigger()
        _spin(self.app, 0.3)
        self.assertFalse(again._filters_panel.isHidden())
        self.assertTrue(again._filters_strip.isHidden())
        self.assertEqual(again._workspace_splitter.sizes()[0], 320)

    def test_no_top_bar_and_new_panel_in_the_window_menu(self) -> None:
        win = self._window()
        # The Filters toggle and the Panels button are gone, and their bar.
        for name in ("ExplorerCommandBar", "FilterRailToggle", "PanelPaletteButton"):
            self.assertIsNone(win.findChild(QtWidgets.QWidget, name), name)
        for name in ("filters_toggle", "panel_palette_button", "_build_command_bar"):
            self.assertFalse(hasattr(win, name), name)
        # The window stripe sits right above the rail and the docking area.
        central = win.centralWidget().layout()
        self.assertEqual(central.count(), 2)
        self.assertIs(central.itemAt(0).widget(), win._window_stripe)
        # Window > New Panel... with its shortcut, right after the window commands.
        menu = win._window_menu
        menu.aboutToShow.emit()
        texts = [a.text() for a in menu.actions() if not a.isSeparator()]
        self.assertEqual(texts[:5], ["New Explore Window", "Duplicate Window", "New Panel...", "Layouts", "Show Filters"])
        action = win._new_panel_action
        self.assertIn(action, menu.actions())
        self.assertEqual(action.shortcut(), QtGui.QKeySequence("Ctrl+Shift+P"))
        self.assertEqual(action.shortcutContext(), QtCore.Qt.ShortcutContext.WindowShortcut)
        self.assertIn(action, win.actions())  # works with the menu closed
        # Filled again, the actions are the window's and survive clear().
        menu.aboutToShow.emit()
        self.assertFalse(sip.isdeleted(action))
        self.assertEqual(sum(1 for a in menu.actions() if a is action), 1)
        with mock.patch.object(win._panel_manager, "show_palette") as palette:
            action.trigger()
        palette.assert_called_once_with()

    def test_ctrl_shift_p_opens_the_picker_over_its_own_window(self) -> None:
        from annotation_explorer_polars.panels import manager as manager_module

        first = self._window()
        second = self._window()
        parents: List[Any] = []

        class FakePalette(QtCore.QObject):
            opened = QtCore.pyqtSignal(str, bool)

            def __init__(self, *_args: Any, parent: Any = None, hint: Any = None, **_kwargs: Any) -> None:
                super().__init__()
                parents.append(parent)
                hints.append(hint)

            def exec(self) -> int:
                return 0

        hints: List[Any] = []
        # No other shortcut of the window uses the same keys.
        keys = QtGui.QKeySequence("Ctrl+Shift+P")
        for win in (first, second):
            bound = [a for a in win.findChildren(QtGui.QAction) if keys in a.shortcuts()]
            bound += [s for s in win.findChildren(QtGui.QShortcut) if s.key() == keys]
            self.assertEqual(bound, [win._new_panel_action])
        with mock.patch.object(manager_module, "PanelPalette", FakePalette):
            for win in (second, first):
                win.activateWindow()
                win.raise_()
                QTest.qWaitForWindowActive(win, 2000)
                target = win.table if win.table.isVisible() else win
                target.setFocus()
                _spin(self.app, 0.05)
                QTest.keyClick(
                    target,
                    QtCore.Qt.Key.Key_P,
                    QtCore.Qt.KeyboardModifier.ControlModifier | QtCore.Qt.KeyboardModifier.ShiftModifier,
                )
                _spin(self.app, 0.05)
        self.assertEqual(parents, [second, first])
        # A main window adds a tab of the active group by default.
        self.assertTrue(all("tab of the active group" in hint for hint in hints), hints)

    def test_autosave_saves_changes_but_not_while_maximized(self) -> None:
        from annotation_explorer_polars.shell.layout_store import LayoutStore, run_key

        win = self._window()
        manager = win._panel_manager
        key = run_key(Path(self.output_dir))
        store = LayoutStore()
        self.assertEqual(win._layout_autosave_timer.interval(), 60_000)
        self.assertTrue(win._layout_autosave_timer.isActive())
        win._autosave_explore_layout()  # the window was resized after the restore
        self.assertIsNotNone(store.load_window(key))
        QtCore.QSettings().remove(LayoutStore.window_group(key))
        win._autosave_explore_layout()
        self.assertIsNone(store.load_window(key))  # nothing changed since
        manager.host.move_view("phylogeny", Placement("right", relative_to="geo"))
        _spin(self.app, 0.2)
        self.assertTrue(manager.host.maximize_group("annotations"))
        win._autosave_explore_layout()
        self.assertIsNone(store.load_window(key))  # skipped while maximized
        self.assertTrue(manager.host.is_maximized())
        manager.host.end_maximize()
        win._autosave_explore_layout()
        saved = store.load_window(key)
        self.assertIsNotNone(saved)
        self.assertEqual(saved.dock_state, manager.host.save_layout())

    def test_new_run_starts_with_the_last_window_size_rail_and_options(self) -> None:
        win = self._window()
        win._apply_rail_state(301, False)
        win.annotations_panel.set_view_mode("Context")
        win._panel_manager.rename_panel("geo", "World map")
        win._panel_manager.host.move_view("phylogeny", Placement("left", relative_to="sample_tree"))
        _spin(self.app, 0.3)
        win.close()
        _spin(self.app, 0.1)
        other_run = tempfile.mkdtemp(prefix="bactolomew_dock_run2_")
        self.addCleanup(shutil.rmtree, other_run, True)
        write_pan_outputs(other_run)
        again = self._window(show=False, output_dir=other_run)
        self.assertTrue(again._geometry_restored)
        # (Offscreen clamps restored sizes to its 800x600 screen, so compare
        # what was stored instead of the size.)
        from annotation_explorer_polars.shell.layout_store import LayoutStore, run_key

        saved = LayoutStore().load_window(run_key(Path(self.output_dir)))
        self.assertEqual(LayoutStore().load_last_window().geometry, saved.geometry)
        self.assertEqual(again._rail_state(), (301, False))
        self.assertEqual(again.view_mode, "Context")
        manager = again._panel_manager
        # The dock layout and the tab titles belong to their run: Classic.
        self.assertEqual(_sorted_groups(manager.host.groups()), _sorted_groups(CLASSIC_GROUPS))
        self.assertEqual(manager.host.view_title("geo"), "Map")
        self.assertEqual(again.gene_query.text(), "")

    def test_table_starts_at_the_top_when_shown_after_a_hidden_refresh(self) -> None:
        win = self._window()
        manager = win._panel_manager
        manager.host.move_view("annotations", Placement("tab", relative_to="sample_tree"))
        manager.host.set_current_tab("sample_tree")
        _spin(self.app, 0.2)
        self.assertFalse(manager.is_visible("annotations"))
        win.gene_query.setText("blaOXA")
        self._wait_refresh(win)
        self.assertTrue(win.annotations_panel.scroll_top_on_show)
        win.table.verticalScrollBar().setValue(win.table.verticalScrollBar().maximum())
        manager.host.set_current_tab("annotations")
        _spin(self.app, 0.3)
        self.assertTrue(manager.is_visible("annotations"))
        self.assertFalse(win.annotations_panel.scroll_top_on_show)
        self.assertEqual(win.table.verticalScrollBar().value(), 0)

    def test_window_without_run_output_leaves_out_its_panels(self) -> None:
        win = self._window(output_dir=None)
        manager = win._panel_manager
        groups = _sorted_groups(manager.host.groups())
        self.assertEqual(groups, [["annotations"], ["geo"], ["sample_tree"], ["track"]])
        for kind in ("phylogeny", "cgmlst_force"):
            available, reason = manager.is_available(kind)
            self.assertFalse(available, kind)
            self.assertEqual(reason, "Not available: this run has no output folder.")
            self.assertFalse(manager.open_panel(kind), kind)
        self.assertFalse(manager.is_open("cgmlst_force"))
        self.assertEqual(manager.host.floating_count(), 0)


if __name__ == "__main__":
    unittest.main()
