"""Layout presets, user layouts and the default for new windows (window plan
Phase 8, 3.8, 3.11).

- Built-in presets (Window > Layouts): Classic, Genome focus, Pangenome and
  Geo + table arrange the docking area only; the rail, the filters and the
  sample set never change, and no compute starts.
- "Save Layout As..." stores the dock state and the panel list (no filters,
  no rail state) under ``explore/v2/layouts/<name>``; applying it in any
  window of the run, or another run, gives the same arrangement.
- "Set as Default for New Windows" (``explore/v2/default``): New Explore
  Window, a window opened from a subset and a run with no saved layout start
  with it instead of Classic; "Use Classic for New Windows" removes it.
- No new keyboard shortcut (D8), no floating panel (D7), nothing crosses
  windows (D2): applying a layout in one window leaves the others alone.
"""
from __future__ import annotations

import os
import shutil
import tempfile
import unittest
from pathlib import Path
from typing import Any, List
from unittest import mock

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets

try:
    from tests.test_dock_host import CLASSIC_GROUPS, _app, _fake_specs, _sorted_groups, _spin
    from tests.test_refresh_parity import build_synthetic_run
    from tests.test_window_manager import _document_panels, _RealWindowCase, _wait
except ImportError:
    from test_dock_host import CLASSIC_GROUPS, _app, _fake_specs, _sorted_groups, _spin
    from test_refresh_parity import build_synthetic_run
    from test_window_manager import _document_panels, _RealWindowCase, _wait

from annotation_explorer_polars.docking import docking_available
from annotation_explorer_polars.model import RunData
from annotation_explorer_polars.shell import commands as cmd
from annotation_explorer_polars.shell.app_controller import AppController
from annotation_explorer_polars.shell.layout_store import (
    DEFAULT_GROUP,
    LAYOUT_NAME_MAX,
    USER_LAYOUTS_GROUP,
    LayoutStore,
    WindowLayout,
)

if docking_available():
    from annotation_explorer_polars.docking import layouts
    from annotation_explorer_polars.panels.manager import PanelManager

CLASSIC_SORTED = _sorted_groups(CLASSIC_GROUPS)
# Classic on the synthetic run, which has no tree and no Panaroo output: the
# bottom row is the map alone.
CLASSIC_RUN_SORTED = [["annotations"], ["geo"], ["sample_tree"], ["track"]]
GENOME_FOCUS_GROUPS = [["annotations"], ["sample_tree"], ["track"]]
# Fake panels have no Plasmidome: Phylogeny and Gene presence/absence above
# the map.
PANGENOME_GROUPS = [["gene_presence_absence"], ["geo"], ["phylogeny"]]
# The synthetic run (no Panaroo output and no tree, plasmid clusters from its
# annotations): the map beside Plasmidome (Phase 9); Phylogeny needs a tree
# and Gene presence/absence Panaroo output.
PANGENOME_WINDOW_GROUPS = [["geo"], ["plasmidome"]]
GEO_TABLE_GROUPS = [["annotations"], ["geo"]]


def _clear_settings() -> None:
    settings = QtCore.QSettings()
    settings.remove("explore")
    settings.sync()


class LayoutStoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        # Keep the application alive (a collected QApplication takes the
        # theme's module-level helpers with it).
        cls.app = _app()

    def setUp(self) -> None:
        _clear_settings()
        self.addCleanup(_clear_settings)
        self.store = LayoutStore()

    @staticmethod
    def _layout(tag: bytes = b"state") -> WindowLayout:
        return WindowLayout(
            geometry=b"geometry",
            rail_width=333,
            rail_collapsed=True,
            dock_state=tag,
            panels=[{"uid": "annotations", "kind": "annotations", "params": {"view_mode": "Global", "title": "Hits"}}],
        )

    def test_user_layouts_round_trip_sorted_by_name(self) -> None:
        self.assertEqual(self.store.layout_names(), [])
        for name in ("beta", "Alpha", "gamma"):
            self.store.save_user_layout(name, self._layout(name.encode()))
        self.assertEqual(self.store.layout_names(), ["Alpha", "beta", "gamma"])
        loaded = self.store.load_user_layout("beta")
        self.assertEqual(loaded.dock_state, b"beta")
        self.assertEqual(loaded.panels, self._layout().panels)  # display options and titles kept

    def test_a_user_layout_keeps_no_filters_rail_or_geometry(self) -> None:
        self.store.save_user_layout("Mine", self._layout())
        settings = QtCore.QSettings()
        settings.beginGroup(f"{USER_LAYOUTS_GROUP}/Mine")
        try:
            self.assertEqual(sorted(settings.childKeys()), ["dockState", "name", "panels", "version"])
        finally:
            settings.endGroup()
        loaded = self.store.load_user_layout("Mine")
        self.assertEqual(loaded.geometry, b"")
        self.assertFalse(loaded.rail_collapsed)

    def test_names_with_any_characters_round_trip(self) -> None:
        names = ["a/b", "x\\y", "50% done", "Map + tree: Europe", "Ünïcode"]
        for name in names:
            self.store.save_user_layout(name, self._layout(name.encode("utf-8")))
        self.assertEqual(sorted(self.store.layout_names()), sorted(names))
        for name in names:
            self.assertTrue(self.store.has_user_layout(name), name)
            self.assertEqual(self.store.load_user_layout(name).dock_state, name.encode("utf-8"))

    def test_names_are_cleaned_and_saving_again_replaces(self) -> None:
        self.assertEqual(LayoutStore.clean_layout_name("  two\n  lines  "), "two lines")
        self.assertEqual(len(LayoutStore.clean_layout_name("x" * 200)), LAYOUT_NAME_MAX)
        with self.assertRaises(ValueError):
            self.store.save_user_layout("   ", self._layout())
        self.store.save_user_layout(" Mine ", self._layout(b"one"))
        self.store.save_user_layout("Mine", self._layout(b"two"))
        self.assertEqual(self.store.layout_names(), ["Mine"])
        self.assertEqual(self.store.load_user_layout("Mine").dock_state, b"two")

    def test_names_that_differ_only_in_case_are_one_layout(self) -> None:
        store = LayoutStore()
        store.save_user_layout("Map", self._layout(b"one"))
        self.assertTrue(store.has_user_layout("map"))
        self.assertEqual(store.existing_layout_name("MAP"), "Map")
        self.assertEqual(store.load_user_layout("map").dock_state, b"one")
        self.assertEqual(store.save_user_layout("map", self._layout(b"two")), "map")
        self.assertEqual(store.layout_names(), ["map"])  # the new spelling is kept
        self.assertEqual(store.load_user_layout("Map").dock_state, b"two")
        self.assertTrue(store.delete_user_layout("MAP"))
        self.assertEqual(store.layout_names(), [])

    def test_delete(self) -> None:
        self.store.save_user_layout("Mine", self._layout())
        self.assertFalse(self.store.delete_user_layout("Other"))
        self.assertTrue(self.store.delete_user_layout("Mine"))
        self.assertEqual(self.store.layout_names(), [])
        self.assertIsNone(self.store.load_user_layout("Mine"))

    def test_default_for_new_windows(self) -> None:
        self.assertFalse(self.store.has_default())
        self.assertIsNone(self.store.load_default())
        self.store.save_default(self._layout(b"default"))
        self.assertTrue(self.store.has_default())
        self.assertEqual(self.store.load_default().dock_state, b"default")
        settings = QtCore.QSettings()
        settings.beginGroup(DEFAULT_GROUP)
        try:
            self.assertEqual(sorted(settings.childKeys()), ["dockState", "panels", "version"])
        finally:
            settings.endGroup()
        self.store.clear_default()
        self.assertFalse(self.store.has_default())
        self.assertIsNone(self.store.load_default())


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class PresetTests(unittest.TestCase):
    """PanelManager presets on fake panels."""

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

    def _groups(self, manager: "PanelManager") -> List[List[str]]:
        _spin(self.app, 0.1)
        return _sorted_groups(manager.host.groups())

    def test_every_preset_builds_its_arrangement(self) -> None:
        manager = self._make()
        expected = {
            layouts.PRESET_GENOME_FOCUS: GENOME_FOCUS_GROUPS,
            layouts.PRESET_PANGENOME: PANGENOME_GROUPS,
            layouts.PRESET_GEO_TABLE: GEO_TABLE_GROUPS,
            layouts.PRESET_CLASSIC: CLASSIC_SORTED,
        }
        self.assertEqual({p.key for p in layouts.PRESETS}, set(expected))
        for key, groups in expected.items():
            self.assertTrue(manager.apply_preset(key), key)
            self.assertEqual(self._groups(manager), groups, key)
            self.assertEqual(manager.host.floating_count(), 0)
            self.assertFalse(manager.host.is_empty())

    def test_preset_proportions(self) -> None:
        manager = self._make()
        manager.apply_preset(layouts.PRESET_GENOME_FOCUS)
        _spin(self.app, 0.3)
        tree = manager.frame("sample_tree").width()
        table = manager.frame("annotations").width()
        self.assertLess(tree, table * 0.5)  # a narrow tree (26:74)
        manager.apply_preset(layouts.PRESET_PANGENOME)
        _spin(self.app, 0.3)
        tree, heatmap = manager.frame("phylogeny"), manager.frame("gene_presence_absence")
        self.assertGreater(tree.height(), manager.frame("geo").height())
        # Side by side at the same height (their rows line up), 35:65.
        self.assertEqual(tree.height(), heatmap.height())
        self.assertLess(tree.width(), heatmap.width() * 0.75)
        manager.apply_preset(layouts.PRESET_CLASSIC)
        _spin(self.app, 0.3)
        # Classic: the bottom row is about 30 percent of the height, map,
        # phylogeny and gene presence/absence side by side (40:25:35).
        bottom = manager.frame("geo").height()
        self.assertEqual(manager.frame("phylogeny").height(), bottom)
        self.assertEqual(manager.frame("gene_presence_absence").height(), bottom)
        self.assertLess(bottom, manager.frame("sample_tree").height() * 0.6)
        widths = [manager.frame(uid).width() for uid in ("geo", "phylogeny", "gene_presence_absence")]
        self.assertGreater(widths[0], widths[1])
        self.assertGreater(widths[2], widths[1])
        manager.apply_preset(layouts.PRESET_GEO_TABLE)
        _spin(self.app, 0.3)
        geo, table = manager.frame("geo").width(), manager.frame("annotations").width()
        self.assertLess(abs(geo - table), 0.1 * (geo + table))  # side by side, evenly

    def test_a_preset_leaves_out_panels_the_run_lacks(self) -> None:
        manager = self._make(unavailable=("phylogeny",))
        self.assertTrue(manager.apply_preset(layouts.PRESET_PANGENOME))
        self.assertEqual(self._groups(manager), [["gene_presence_absence"], ["geo"]])

    def test_presets_never_place_the_track_without_annotations(self) -> None:
        # Feedback 2026-10-07: the Genome track cannot exist without Annotations.
        manager = self._make(unavailable=("annotations",))
        self.assertEqual(
            self._groups(manager), [["gene_presence_absence"], ["geo"], ["phylogeny"], ["sample_tree"]]
        )
        self.assertTrue(manager.apply_preset(layouts.PRESET_GENOME_FOCUS))
        self.assertEqual(self._groups(manager), [["sample_tree"]])
        self.assertFalse(manager.is_open("track"))
        self.assertEqual(layouts.preset(layouts.PRESET_GENOME_FOCUS).needs_any, ("sample_tree", "annotations"))
        manager = self._make(unavailable=("annotations", "sample_tree"))
        available, reason = manager.preset_available(layouts.PRESET_GENOME_FOCUS)
        self.assertFalse(available)
        self.assertNotIn("Track", reason)
        # The builders on their own leave the track out too.
        manager.host.clear_layout()
        layouts.build_classic(manager.host, ["sample_tree", "track", "geo"])
        layouts.build_genome_focus(manager.host, ["sample_tree", "track"])
        self.assertFalse(manager.host.is_placed("track"))

    def test_a_preset_without_any_of_its_main_panels_is_refused(self) -> None:
        manager = self._make(unavailable=("phylogeny", "gene_presence_absence"))
        before = self._groups(manager)
        available, reason = manager.preset_available(layouts.PRESET_PANGENOME)
        self.assertFalse(available)
        self.assertIn("Phylogeny", reason)
        self.assertFalse(manager.apply_preset(layouts.PRESET_PANGENOME))
        self.assertEqual(self._groups(manager), before)
        self.assertFalse(manager.apply_preset("no_such_preset"))

    def test_a_preset_resets_titles_ends_maximize_and_closes_the_rest(self) -> None:
        manager = self._make()
        manager.open_panel("pan_metadata")
        manager.rename_panel("geo", "Europe")
        manager.host.maximize_group("sample_tree")
        self.assertTrue(manager.apply_preset(layouts.PRESET_GEO_TABLE))
        self.assertFalse(manager.host.is_maximized())
        self.assertEqual(manager.title("geo"), "Geo")
        self.assertFalse(manager.is_open("pan_metadata"))
        self.assertFalse(manager.is_open("sample_tree"))
        self.assertEqual(manager.closed_stack(), [])
        # A panel the preset closed opens again from the palette or menu.
        self.assertTrue(manager.open_panel("sample_tree"))

    def test_a_saved_arrangement_applies_and_a_broken_one_changes_nothing(self) -> None:
        manager = self._make()
        manager.host.move_view("phylogeny", layouts.Placement("right", relative_to="sample_tree"))
        manager.rename_panel("phylogeny", "Tree")
        saved = self._groups(manager)
        state, panels = manager.host.save_layout(), manager.panel_list()
        manager.apply_preset(layouts.PRESET_GEO_TABLE)
        self.assertEqual(manager.title("phylogeny"), "Phylogeny")
        self.assertTrue(manager.apply_saved_layout(state, panels))
        self.assertEqual(self._groups(manager), saved)
        self.assertEqual(manager.title("phylogeny"), "Tree")
        self.assertFalse(manager.apply_saved_layout(b"not a dock state", panels))
        self.assertEqual(self._groups(manager), saved)
        self.assertFalse(manager.apply_saved_layout(b"", panels))
        self.assertEqual(self._groups(manager), saved)
        self.assertEqual(manager.host.floating_count(), 0)

    def test_a_layout_saved_while_a_group_is_maximized_keeps_the_full_sizes(self) -> None:
        manager = self._make()
        manager.apply_preset(layouts.PRESET_GEO_TABLE)
        _spin(self.app, 0.3)
        geo = manager.frame("geo").width()
        normal = manager.host.save_layout()
        self.assertTrue(manager.host.maximize_group("annotations"))
        _spin(self.app, 0.2)
        maximized = manager.host.save_layout()  # ends the maximize first
        self.assertFalse(manager.host.is_maximized())
        panels = manager.panel_list()
        for state in (normal, maximized):
            manager.apply_preset(layouts.PRESET_CLASSIC)
            _spin(self.app, 0.2)
            self.assertTrue(manager.apply_saved_layout(state, panels))
            _spin(self.app, 0.3)
            self.assertLess(abs(manager.frame("geo").width() - geo), 0.1 * geo)

    def test_has_open_panels(self) -> None:
        manager = self._make()
        self.assertTrue(manager.has_open_panels())
        for uid in list(manager.host.open_uids()):
            manager.close_panel(uid)
        self.assertFalse(manager.has_open_panels())

    def test_a_saved_arrangement_closes_kinds_this_run_lacks(self) -> None:
        source = self._make()
        state, panels = source.host.save_layout(), source.panel_list()
        target = self._make(unavailable=("phylogeny",))
        target.apply_preset(layouts.PRESET_GEO_TABLE)
        self.assertTrue(target.apply_saved_layout(state, panels))
        self.assertFalse(target.is_open("phylogeny"))
        self.assertEqual(
            self._groups(target), [["annotations"], ["gene_presence_absence"], ["geo"], ["sample_tree"], ["track"]]
        )


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class WindowLayoutsTests(_RealWindowCase):
    """The Window > Layouts menu and the default for new windows on real
    windows of the synthetic run."""

    def setUp(self) -> None:
        super().setUp()
        self.extra_dirs: List[str] = []

    def tearDown(self) -> None:
        super().tearDown()
        for path in self.extra_dirs:
            shutil.rmtree(path, ignore_errors=True)

    def _controller(self, output_dir: Any = "default", **kwargs: Any) -> AppController:
        if output_dir == "default":
            return AppController(self._run_data(), **kwargs)
        ann_df, meta_df = build_synthetic_run()
        run = RunData.build(
            ann_df,
            meta_df,
            ref_df=pd.DataFrame(),
            marker_hits_df=pd.DataFrame(),
            output_dir=Path(output_dir) if output_dir is not None else None,
        )
        return AppController(run, **kwargs)

    def _first(self, controller: AppController) -> Any:
        win = controller.open_first_window()
        self.windows.append(win)
        self._settle(win)
        return win

    def _open(self, controller: AppController, **kwargs: Any) -> Any:
        win = controller.open_explore(**kwargs)
        self.windows.append(win)
        self._settle(win)
        return win

    def _groups(self, win: Any) -> List[List[str]]:
        _spin(self.app, 0.1)
        return _sorted_groups(win._panel_manager.host.groups())

    @staticmethod
    def _layouts_menu(win: Any) -> QtWidgets.QMenu:
        win._window_menu.aboutToShow.emit()
        action = next(a for a in win._window_menu.actions() if a.menu() is not None and a.text() == "Layouts")
        menu = action.menu()
        menu.aboutToShow.emit()
        return menu

    @staticmethod
    def _named(menu: QtWidgets.QMenu, name: str) -> QtGui.QAction:
        return next(a for a in menu.actions() if a.objectName() == name)

    def _computes(self, win: Any) -> List[int]:
        counter = [0]
        original = win.engine.compute

        def counting(*args: Any, **kwargs: Any) -> Any:
            counter[0] += 1
            return original(*args, **kwargs)

        patcher = mock.patch.object(win.engine, "compute", side_effect=counting)
        patcher.start()
        self.addCleanup(patcher.stop)
        return counter

    # ----------------------------------------------------------- the menu

    def test_the_layouts_menu_sits_in_the_window_menu_without_shortcuts(self) -> None:
        win = self._first(self._controller())
        menu = self._layouts_menu(win)
        texts = [a.text() for a in menu.actions() if not a.isSeparator()]
        self.assertEqual(
            texts,
            [
                "Built-in layouts", "Classic", "Genome focus", "Pangenome", "Geo + table",
                "Saved layouts", "No saved layouts",
                "Save Layout As...", "Set as Default for New Windows", "Use Classic for New Windows",
                "Delete Saved Layout",
            ],
        )
        for action in menu.actions():
            self.assertTrue(action.shortcut().isEmpty(), action.text())
        self.assertFalse(self._named(menu, "ClearDefaultLayoutAction").isEnabled())
        self.assertFalse(self._named(menu, "NoSavedLayouts").isEnabled())
        self.assertFalse(menu.findChild(QtWidgets.QMenu, "DeleteLayoutMenu").isEnabled())
        for key in ("classic", "genome_focus", "pangenome", "geo_table"):
            self.assertTrue(self._named(menu, f"LayoutPreset_{key}").isEnabled(), key)
        # Filling the Window menu again keeps the submenu alive.
        win._window_menu.aboutToShow.emit()
        win._window_menu.aboutToShow.emit()
        self.assertIs(self._layouts_menu(win), win._layouts_menu)
        self._layouts_menu(win)
        # Refilling never leaves submenus behind.
        self.assertEqual(len(win._layouts_menu.findChildren(QtWidgets.QMenu)), 1)
        self.assertEqual(len(win._window_menu.findChildren(QtWidgets.QMenu)), 2)
        for name in (cmd.LAYOUT_SAVE, cmd.LAYOUT_APPLY, cmd.LAYOUT_SET_DEFAULT, cmd.LAYOUT_CLEAR_DEFAULT, cmd.LAYOUT_DELETE):
            self.assertTrue(win.app_controller.commands.has(name), name)
        self.assertEqual(win.findChildren(QtGui.QShortcut), [])

    def test_pangenome_is_disabled_without_a_run_folder(self) -> None:
        win = self._first(self._controller(output_dir=None))
        menu = self._layouts_menu(win)
        action = self._named(menu, "LayoutPreset_pangenome")
        self.assertFalse(action.isEnabled())
        self.assertIn("Not available", action.toolTip())
        self.assertFalse(win.app_controller.commands.run(cmd.LAYOUT_APPLY, win, preset="pangenome"))

    # -------------------------------------------------------- presets

    def test_a_preset_never_touches_the_rail_filters_or_sample_set(self) -> None:
        win = self._first(self._controller())
        win.gene_query.setText("blaOXA")
        win._set_sample_subset_filter(label="Three", samples=["S1", "S2", "S3"])
        self.assertTrue(_wait(self.app, lambda: not win.ctx.is_busy(), 20))
        _spin(self.app, 0.2)
        rail_width = win._rail_state()
        applied, result, samples = win.ctx.applied, win.ctx.result, win.ctx.sample_set()
        computes = self._computes(win)
        menu = self._layouts_menu(win)
        for key, groups in (
            ("genome_focus", GENOME_FOCUS_GROUPS),
            ("pangenome", PANGENOME_WINDOW_GROUPS),
            ("geo_table", GEO_TABLE_GROUPS),
            ("classic", CLASSIC_RUN_SORTED),
        ):
            self._named(menu, f"LayoutPreset_{key}").trigger()
            self.assertEqual(self._groups(win), groups, key)
            _spin(self.app, 0.3)
            self.assertIs(win.ctx.applied, applied)
            self.assertIs(win.ctx.result, result)
            self.assertEqual(win.ctx.sample_set(), samples)
            self.assertEqual(win.gene_query.text(), "blaOXA")
            self.assertEqual(win._rail_state(), rail_width)
            self.assertFalse(win._filters_panel.isHidden())
        self.assertEqual(computes[0], 0)
        self.assertEqual(win._panel_manager.host.floating_count(), 0)

    def test_panels_shown_by_a_preset_catch_up_with_the_window_result(self) -> None:
        win = self._first(self._controller())
        manager = win._panel_manager
        self.assertFalse(manager.is_open("plasmidome"))
        win._set_sample_subset_filter(label="Two", samples=["S1", "S2"])
        self.assertTrue(_wait(self.app, lambda: not win.ctx.is_busy(), 20))
        samples = win.ctx.sample_set()
        win.app_controller.commands.run(cmd.LAYOUT_APPLY, win, preset="pangenome")
        for kind in ("plasmidome", "geo"):
            frame = manager.frame(kind)
            self.assertTrue(_wait(self.app, lambda f=frame: f.isVisible() and not f.is_dirty(), 10), kind)
        self.assertEqual(win.ctx.sample_set(), samples)

    def test_documents_stay_open_when_a_layout_is_applied(self) -> None:
        win = self._first(self._controller())
        win.app_controller.commands.run(
            cmd.PANEL_OPEN_NEIGHBOURHOOD,
            win,
            rows=[],
            title="Gene neighborhoods: geneA",
            key=("test", "geneA"),
        )
        _spin(self.app, 0.3)
        documents = _document_panels(win, "neighbourhood")
        self.assertEqual(len(documents), 1)
        uid = documents[0].uid
        manager = win._panel_manager
        self.assertTrue(manager.apply_preset("geo_table"))
        self.assertTrue(manager.is_open(uid))
        win._save_layout_as("Mine")
        manager.apply_preset("genome_focus")
        self.assertTrue(win._apply_layout(name="Mine"))
        self.assertTrue(manager.is_open(uid))
        self.assertEqual(manager.host.floating_count(), 0)

    # ------------------------------------------------------- user layouts

    def test_save_layout_as_and_apply_it_in_another_window(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        manager = first._panel_manager
        extra = manager.duplicate_panel("annotations", "right")
        manager.panel(extra).set_view_mode("Global")
        manager.rename_panel(extra, "All hits")
        manager.close_panel("geo")
        saved_groups = self._groups(first)
        menu = self._layouts_menu(first)
        with mock.patch.object(QtWidgets.QInputDialog, "getText", return_value=("  My layout ", True)):
            self._named(menu, "SaveLayoutAsAction").trigger()
        self.assertEqual(controller.layouts.layout_names(), ["My layout"])
        second = self._open(controller, source=first)
        second.gene_query.setText("geneB")
        self.assertTrue(_wait(self.app, lambda: not second.ctx.is_busy(), 20))
        applied = second.ctx.applied
        self.assertEqual(self._groups(second), CLASSIC_RUN_SORTED)
        menu = self._layouts_menu(second)
        action = next(a for a in menu.actions() if a.objectName() == "UserLayoutAction" and a.data() == "My layout")
        action.trigger()
        _spin(self.app, 0.3)
        groups = self._groups(second)
        self.assertEqual(len(groups), len(saved_groups))
        copies = [uid for uid in second._panel_manager.uids("annotations") if uid != "annotations"]
        self.assertEqual(len(copies), 1)
        self.assertEqual(second._panel_manager.title(copies[0]), "All hits")
        self.assertEqual(second._panel_manager.panel(copies[0]).view_mode(), "Global")
        self.assertFalse(second._panel_manager.is_open("geo"))
        # Filters stay; the first window is untouched.
        self.assertIs(second.ctx.applied, applied)
        self.assertEqual(second.gene_query.text(), "geneB")
        self.assertEqual(self._groups(first), saved_groups)

    def test_save_layout_as_asks_before_replacing(self) -> None:
        win = self._first(self._controller())
        store = win.app_controller.layouts
        self.assertEqual(win._save_layout_as("Mine"), "Mine")
        original = store.load_user_layout("Mine").dock_state
        win._panel_manager.apply_preset("geo_table")
        with mock.patch.object(QtWidgets.QInputDialog, "getText", return_value=("Mine", True)), mock.patch.object(
            QtWidgets.QMessageBox, "question", return_value=QtWidgets.QMessageBox.StandardButton.No
        ) as question:
            self.assertIsNone(win.app_controller.commands.run(cmd.LAYOUT_SAVE, win))
        question.assert_called_once()
        self.assertEqual(store.load_user_layout("Mine").dock_state, original)
        with mock.patch.object(QtWidgets.QInputDialog, "getText", return_value=("Mine", True)), mock.patch.object(
            QtWidgets.QMessageBox, "question", return_value=QtWidgets.QMessageBox.StandardButton.Yes
        ):
            self.assertEqual(win.app_controller.commands.run(cmd.LAYOUT_SAVE, win), "Mine")
        self.assertNotEqual(store.load_user_layout("Mine").dock_state, original)
        with mock.patch.object(QtWidgets.QInputDialog, "getText", return_value=("Other", False)):
            self.assertIsNone(win.app_controller.commands.run(cmd.LAYOUT_SAVE, win))
        with mock.patch.object(QtWidgets.QInputDialog, "getText", return_value=("   ", True)):
            self.assertIsNone(win.app_controller.commands.run(cmd.LAYOUT_SAVE, win))
        self.assertEqual(store.layout_names(), ["Mine"])

    def test_a_broken_saved_layout_leaves_the_panels_alone(self) -> None:
        win = self._first(self._controller())
        win._panel_manager.apply_preset("genome_focus")
        before = self._groups(win)
        store = win.app_controller.layouts
        store.save_user_layout("Broken", WindowLayout(dock_state=b"garbage", panels=[]))
        self.assertFalse(win.app_controller.commands.run(cmd.LAYOUT_APPLY, win, layout="Broken"))
        self.assertEqual(self._groups(win), before)
        self.assertIn("could not be applied", win.layout_message.text())
        self.assertTrue(win.layout_message.isVisibleTo(win.statusBar()))
        # The sample count stays in view next to the message.
        self.assertTrue(win.status.isVisibleTo(win.statusBar()))
        self.assertFalse(win.app_controller.commands.run(cmd.LAYOUT_APPLY, win, layout="Missing"))
        self.assertEqual(self._groups(win), before)

    def test_delete_a_saved_layout_from_the_menu(self) -> None:
        win = self._first(self._controller())
        win._save_layout_as("One")
        win._save_layout_as("Two")
        menu = self._layouts_menu(win)
        delete_menu = menu.findChild(QtWidgets.QMenu, "DeleteLayoutMenu")
        self.assertTrue(delete_menu.isEnabled())
        self.assertEqual([a.text() for a in delete_menu.actions()], ["One", "Two"])
        with mock.patch.object(
            QtWidgets.QMessageBox, "question", return_value=QtWidgets.QMessageBox.StandardButton.No
        ):
            delete_menu.actions()[0].trigger()
        self.assertEqual(win.app_controller.layouts.layout_names(), ["One", "Two"])
        with mock.patch.object(
            QtWidgets.QMessageBox, "question", return_value=QtWidgets.QMessageBox.StandardButton.Yes
        ):
            delete_menu.actions()[0].trigger()
        self.assertEqual(win.app_controller.layouts.layout_names(), ["Two"])

    # ------------------------------------------------- default for new windows

    def test_set_as_default_feeds_new_windows_and_subset_windows(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        first._panel_manager.apply_preset("geo_table")
        first.gene_query.setText("geneA")
        self.assertTrue(_wait(self.app, lambda: not first.ctx.is_busy(), 20))
        menu = self._layouts_menu(first)
        self._named(menu, "SetDefaultLayoutAction").trigger()
        self.assertTrue(controller.layouts.has_default())
        new_window = self._open(controller, source=first)
        self.assertEqual(self._groups(new_window), GEO_TABLE_GROUPS)
        self.assertEqual(new_window.gene_query.text(), "")  # a clean rail
        self.assertFalse(new_window._filters_panel.isHidden())
        from annotation_explorer_polars.model.filter_state import SampleSubset

        subset = SampleSubset(label="Clade 1", samples=("S1", "S2"), kind="clade")
        subset_window = self._open(controller, source=first, window_subset=subset)
        self.assertEqual(self._groups(subset_window), GEO_TABLE_GROUPS)
        # Window 1 keeps its own layout and filters.
        self.assertEqual(self._groups(first), GEO_TABLE_GROUPS)
        self.assertEqual(first.gene_query.text(), "geneA")
        # Duplicate Window still copies the source's own layout.
        first._panel_manager.apply_preset("genome_focus")
        duplicate = controller.duplicate_explore(first)
        self.windows.append(duplicate)
        self._settle(duplicate)
        self.assertEqual(self._groups(duplicate), GENOME_FOCUS_GROUPS)
        # Use Classic for New Windows.
        menu = self._layouts_menu(first)
        clear = self._named(menu, "ClearDefaultLayoutAction")
        self.assertTrue(clear.isEnabled())
        clear.trigger()
        self.assertFalse(controller.layouts.has_default())
        classic = self._open(controller, source=first)
        self.assertEqual(self._groups(classic), CLASSIC_RUN_SORTED)

    def test_a_run_without_a_saved_layout_starts_with_the_default(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        manager = first._panel_manager
        extra = manager.duplicate_panel("annotations", "right")
        manager.panel(extra).set_view_mode("Context")
        manager.apply_preset("genome_focus")
        manager.duplicate_panel("annotations", "tab")
        expected = self._groups(first)
        self.assertTrue(first.app_controller.commands.run(cmd.LAYOUT_SET_DEFAULT, first))
        first.close()
        other_dir = tempfile.mkdtemp(prefix="bactolomew_layout_run_")
        self.extra_dirs.append(other_dir)
        other = self._first(self._controller(output_dir=other_dir))
        self.assertEqual(len(self._groups(other)), len(expected))
        self.assertEqual(
            sorted(len(group) for group in self._groups(other)), sorted(len(group) for group in expected)
        )
        self.assertEqual(len(other._panel_manager.uids("annotations")), 2)
        # The run's own saved layout wins over the default.
        controller = self._controller()
        reopened = self._first(controller)
        self.assertEqual(self._groups(reopened), expected)

    def test_a_default_naming_panels_a_run_lacks_closes_them(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        first._panel_manager.apply_preset("pangenome")
        first.app_controller.commands.run(cmd.LAYOUT_SET_DEFAULT, first)
        first.close()
        no_folder = self._first(self._controller(output_dir=None))
        manager = no_folder._panel_manager
        for kind in ("phylogeny", "gene_presence_absence", "plasmidome"):
            self.assertFalse(manager.is_open(kind), kind)
        self.assertEqual(self._groups(no_folder), [["geo"]])

    def test_an_unusable_default_falls_back_to_classic(self) -> None:
        controller = self._controller()
        controller.layouts.save_default(WindowLayout(dock_state=b"garbage", panels=[]))
        first = self._first(controller)
        self.assertEqual(self._groups(first), CLASSIC_RUN_SORTED)
        second = self._open(controller, source=first)
        self.assertEqual(self._groups(second), CLASSIC_RUN_SORTED)

    # ------------------------------------------------- verification follow-ups

    def _open_document(self, win: Any, gene: str = "geneA") -> str:
        win.app_controller.commands.run(
            cmd.PANEL_OPEN_NEIGHBOURHOOD,
            win,
            rows=[],
            title=f"Gene neighborhoods: {gene}",
            key=("test", gene),
        )
        _spin(self.app, 0.3)
        documents = _document_panels(win, "neighbourhood")
        self.assertEqual(len(documents), 1)
        return documents[0].uid

    def test_messages_never_hide_the_sample_count(self) -> None:
        win = self._first(self._controller())
        bar = win.statusBar()
        self.assertTrue(win.status.isVisibleTo(bar))
        self.assertEqual(win._save_layout_as("Mine"), "Mine")
        self.assertEqual(win.layout_message.text(), 'Layout "Mine" saved')
        self.assertTrue(win.layout_message.isVisibleTo(bar))
        self.assertTrue(win.status.isVisibleTo(bar))
        self.assertEqual(bar.currentMessage(), "")
        win._layout_message_timer.timeout.emit()
        self.assertFalse(win.layout_message.isVisibleTo(bar))

    def test_menu_tooltips_show_and_ampersands_stay_literal(self) -> None:
        win = self._first(self._controller(output_dir=None))
        win._save_layout_as("Map & table")
        menu = self._layouts_menu(win)
        self.assertTrue(menu.toolTipsVisible())
        pangenome = self._named(menu, "LayoutPreset_pangenome")
        self.assertFalse(pangenome.isEnabled())
        self.assertIn("Not available for this run", pangenome.toolTip())
        action = next(a for a in menu.actions() if a.objectName() == "UserLayoutAction")
        self.assertEqual(action.data(), "Map & table")
        self.assertEqual(action.text(), "Map && table")
        self.assertEqual(action.iconText(), "Map & table")
        delete_menu = menu.findChild(QtWidgets.QMenu, "DeleteLayoutMenu")
        self.assertEqual([a.iconText() for a in delete_menu.actions()], ["Map & table"])
        action.trigger()
        self.assertTrue(win._panel_manager.has_open_panels())

    def test_an_empty_arrangement_is_not_saved(self) -> None:
        win = self._first(self._controller())
        manager = win._panel_manager
        self._open_document(win)
        for uid in list(manager.host.open_uids()):
            if not manager.frame(uid).spec.document:
                manager.close_panel(uid)
        self.assertFalse(manager.has_open_panels())
        store = win.app_controller.layouts
        with mock.patch.object(QtWidgets.QInputDialog, "getText") as ask:
            self.assertIsNone(win.app_controller.commands.run(cmd.LAYOUT_SAVE, win))
        ask.assert_not_called()
        self.assertIsNone(win._save_layout_as("Empty"))
        self.assertFalse(win.app_controller.commands.run(cmd.LAYOUT_SET_DEFAULT, win))
        self.assertEqual(store.layout_names(), [])
        self.assertFalse(store.has_default())
        self.assertIn("empty layout", win.layout_message.text())

    def test_document_panels_keep_a_usable_column(self) -> None:
        win = self._first(self._controller())
        manager = win._panel_manager
        uid = self._open_document(win)
        host_width = manager.host.width()

        def share() -> float:
            _spin(self.app, 0.4)
            return manager.frame(uid).width() / float(host_width)

        for key in ("classic", "geo_table", "genome_focus", "pangenome"):
            self.assertTrue(manager.apply_preset(key), key)
            self.assertTrue(manager.is_open(uid), key)
            self.assertTrue(0.25 < share() < 0.42, (key, share()))
        # The other groups keep the preset's own proportions.
        manager.apply_preset("geo_table")
        _spin(self.app, 0.4)
        geo, table = manager.frame("geo").width(), manager.frame("annotations").width()
        self.assertLess(abs(geo - table), 0.1 * (geo + table))
        win._save_layout_as("With document")
        manager.apply_preset("genome_focus")
        self.assertTrue(win._apply_layout(name="With document"))
        self.assertTrue(0.25 < share() < 0.42, share())

    def test_a_layout_saved_with_a_document_applies_where_that_document_is_not_open(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        self._open_document(first)
        first._panel_manager.apply_preset("geo_table")
        self.assertEqual(first.app_controller.commands.run(cmd.LAYOUT_SAVE, first, layout="Doc"), "Doc")
        first.app_controller.commands.run(cmd.LAYOUT_SET_DEFAULT, first)
        second = self._open(controller, source=first)
        # The default carried the document's dock name: ignored, no empty group.
        self.assertEqual(self._groups(second), GEO_TABLE_GROUPS)
        self.assertEqual(_document_panels(second, "neighbourhood"), [])
        second._panel_manager.apply_preset("classic")
        self.assertTrue(second.app_controller.commands.run(cmd.LAYOUT_APPLY, second, layout="Doc"))
        self.assertEqual(self._groups(second), GEO_TABLE_GROUPS)
        self.assertEqual(second._panel_manager.host.floating_count(), 0)
        self.assertEqual(len(_document_panels(first, "neighbourhood")), 1)

    def test_layouts_never_carry_or_change_the_rail(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        first._apply_rail_state(400, True)
        _spin(self.app, 0.2)
        self.assertTrue(first._rail_state()[1])
        first._panel_manager.apply_preset("geo_table")
        first.app_controller.commands.run(cmd.LAYOUT_SAVE, first, layout="Collapsed")
        first.app_controller.commands.run(cmd.LAYOUT_SET_DEFAULT, first)
        second = self._open(controller, source=first)
        self.assertFalse(second._rail_state()[1])  # a clean rail, not the source's
        second._apply_rail_state(300, False)
        _spin(self.app, 0.2)
        before = second._rail_state()
        self.assertTrue(second.app_controller.commands.run(cmd.LAYOUT_APPLY, second, layout="Collapsed"))
        _spin(self.app, 0.2)
        self.assertEqual(second._rail_state(), before)
        self.assertFalse(second._filters_panel.isHidden())

    def test_reset_layout_builds_classic_even_with_a_default_set(self) -> None:
        win = self._first(self._controller())
        win._panel_manager.apply_preset("geo_table")
        self.assertTrue(win.app_controller.commands.run(cmd.LAYOUT_SET_DEFAULT, win))
        win._panel_manager.reset_layout()
        self.assertEqual(self._groups(win), CLASSIC_RUN_SORTED)
        win._panel_manager.apply_preset("geo_table")
        win.app_controller.commands.run(cmd.LAYOUT_RESET, win)
        self.assertEqual(self._groups(win), CLASSIC_RUN_SORTED)
        self.assertTrue(win.app_controller.layouts.has_default())

    def test_classic_when_neither_the_layout_nor_the_previous_one_applies(self) -> None:
        win = self._first(self._controller())
        manager = win._panel_manager
        uid = self._open_document(win)
        manager.apply_preset("geo_table")
        state = manager.host.save_layout()
        with mock.patch.object(manager.host, "restore_layout", return_value=False):
            self.assertFalse(manager.apply_saved_layout(state, manager.panel_list()))
        _spin(self.app, 0.2)
        groups = self._groups(win)
        self.assertIn([uid], groups)
        groups.remove([uid])
        self.assertEqual(groups, CLASSIC_RUN_SORTED)
        self.assertEqual(manager.host.floating_count(), 0)

    def test_a_layout_saved_by_command_reaches_a_later_launch(self) -> None:
        win = self._first(self._controller())
        win._panel_manager.apply_preset("genome_focus")
        self.assertEqual(win.app_controller.commands.run(cmd.LAYOUT_SAVE, win, layout="Focus"), "Focus")
        win.close()
        later = self._first(self._controller(output_dir=None))
        self.assertEqual(later.app_controller.layouts.layout_names(), ["Focus"])
        self.assertTrue(later.app_controller.commands.run(cmd.LAYOUT_APPLY, later, layout="Focus"))
        self.assertEqual(self._groups(later), GENOME_FOCUS_GROUPS)

    def test_an_empty_session_slot_starts_with_the_default(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        first._panel_manager.apply_preset("geo_table")
        first.app_controller.commands.run(cmd.LAYOUT_SET_DEFAULT, first)
        first._panel_manager.apply_preset("genome_focus")
        first.close()
        controller = self._controller()
        key = controller.run_key()
        self.assertIsNotNone(controller.layouts.load_window(key, 1))
        self.assertIsNone(controller.layouts.load_window(key, 2))
        controller.layouts.save_session(key, [1, 2])
        reopened = self._first(controller)
        windows = controller.windows.explore_windows()
        self.windows += [w for w in windows if w not in self.windows]
        self.assertEqual(len(windows), 2)
        self._settle(windows[1])
        self.assertEqual(self._groups(reopened), GENOME_FOCUS_GROUPS)  # its own slot
        self.assertEqual(self._groups(windows[1]), GEO_TABLE_GROUPS)  # nothing saved: the default


if __name__ == "__main__":
    unittest.main()
