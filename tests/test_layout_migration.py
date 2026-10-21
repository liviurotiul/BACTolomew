"""Layout store (window architecture plan 3.11), the retired pre-docking keys
and the migration of format version 1 layouts.

Format version 1 named neoPan, the Pan tree, neoTree and Gene presence; they
became Phylogeny and Gene presence/absence (feedback 2026-10-07). A version 1
layout (window slots, so sessions; the last window; the default; user
layouts) is migrated in memory when loaded and restores without error.

The old build kept its window state in the QSettings group
``annotationExplorer``. Phases 1 to 9 read it once for migration; Phase 10
retired that read with the legacy layout, so the docking build neither reads
nor writes the old group: a window starts from its ``explore/v2/...`` layout,
else Classic.
"""
from __future__ import annotations

import os
import shutil
import tempfile
import time
import unittest
from pathlib import Path
from typing import Any, List

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

import pandas as pd
from PyQt6 import QtCore, QtWidgets

from annotation_explorer_polars.docking import docking_available
from annotation_explorer_polars.shell.layout_store import (
    FORMAT_VERSION,
    LayoutStore,
    WindowLayout,
    clamp_rail_width,
    migrate_layout,
    run_key,
)

try:
    from tests._pan_fixture import write_pan_outputs
    from tests.test_refresh_parity import build_synthetic_run
except ImportError:
    from _pan_fixture import write_pan_outputs
    from test_refresh_parity import build_synthetic_run

# The pre-docking QSettings group (no longer read or written).
OLD_GROUP = "annotationExplorer"
# Classic on the synthetic run (no tree, no Panaroo output: the map alone below).
CLASSIC_GROUPS = [["annotations"], ["geo"], ["sample_tree"], ["track"]]


def _spin(app: QtWidgets.QApplication, seconds: float) -> None:
    end = time.monotonic() + seconds
    while True:
        app.processEvents(QtCore.QEventLoop.ProcessEventsFlag.AllEvents, 20)
        if time.monotonic() >= end:
            break
        time.sleep(0.005)


def _clear_settings() -> None:
    settings = QtCore.QSettings()
    settings.remove(OLD_GROUP)
    settings.remove("explore")
    settings.sync()


class LayoutStoreTests(unittest.TestCase):
    def setUp(self) -> None:
        self.dir = tempfile.mkdtemp(prefix="bactolomew_layout_store_")
        self.settings = QtCore.QSettings(os.path.join(self.dir, "s.ini"), QtCore.QSettings.Format.IniFormat)
        self.store = LayoutStore(self.settings)

    def tearDown(self) -> None:
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_run_key_is_stable_and_per_run(self) -> None:
        self.assertEqual(run_key(Path("/a/run1")), run_key(Path("/a/run1/")))
        self.assertNotEqual(run_key(Path("/a/run1")), run_key(Path("/a/run2")))
        self.assertEqual(len(run_key(None)), 40)

    def test_window_round_trip(self) -> None:
        key = run_key(Path("/some/run"))
        self.assertIsNone(self.store.load_window(key))
        layout = WindowLayout(
            geometry=b"geo",
            rail_width=300,
            rail_collapsed=True,
            dock_state=b"\x00\x01state",
            panels=[{"uid": "annotations", "kind": "annotations", "params": {"view_mode": "Global"}}],
        )
        self.store.save_window(key, layout)
        loaded = self.store.load_window(key)
        self.assertEqual(loaded, layout)
        self.assertEqual(loaded.version, FORMAT_VERSION)
        self.assertTrue(self.settings.contains(f"explore/v2/run/{key}/windows/1/dockState"))

    def test_last_window_keeps_no_dock_state_and_no_titles(self) -> None:
        self.assertIsNone(self.store.load_last_window())
        layout = WindowLayout(
            geometry=b"geo",
            rail_width=301,
            rail_collapsed=True,
            dock_state=b"\x00state",
            panels=[
                {"uid": "annotations", "kind": "annotations", "params": {"view_mode": "Global", "title": "Genes"}},
                {"uid": "geo", "kind": "geo", "params": {"title": "World map"}},
            ],
        )
        self.store.save_last_window(layout)
        last = self.store.load_last_window()
        self.assertEqual((last.geometry, last.rail_width, last.rail_collapsed), (b"geo", 301, True))
        self.assertEqual(last.dock_state, b"")
        self.assertEqual([p["params"] for p in last.panels], [{"view_mode": "Global"}, {}])
        self.assertTrue(self.settings.contains("explore/v2/lastWindow/version"))
        # The run's own entry is untouched.
        self.assertIsNone(self.store.load_window(run_key(Path("/some/run"))))

    def test_rail_width_is_clamped(self) -> None:
        self.assertEqual(clamp_rail_width(10), 248)
        self.assertEqual(clamp_rail_width(9999), 420)
        self.assertEqual(clamp_rail_width("junk"), 288)

    def test_the_store_has_no_old_key_reader(self) -> None:
        self.assertFalse(hasattr(self.store, "read_legacy_once"))
        self.assertFalse(hasattr(self.store, "mark_legacy_migrated"))


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class OldKeysTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        from annotation_explorer_polars.theme import apply_application_theme

        apply_application_theme(cls.app)
        import annotation_explorer_polars.geo_experimental_views as gev

        cls._saved_webengine_cache = gev._WEBENGINE_RUNTIME_CACHE
        gev._WEBENGINE_RUNTIME_CACHE = (False, "disabled in tests")

    @classmethod
    def tearDownClass(cls) -> None:
        import annotation_explorer_polars.geo_experimental_views as gev

        gev._WEBENGINE_RUNTIME_CACHE = cls._saved_webengine_cache

    def setUp(self) -> None:
        _clear_settings()
        self.output_dir = tempfile.mkdtemp(prefix="bactolomew_migration_run_")
        self.windows: List[Any] = []

    def tearDown(self) -> None:
        for win in self.windows:
            win.close()
            win.deleteLater()
        _spin(self.app, 0.05)
        _clear_settings()
        shutil.rmtree(self.output_dir, ignore_errors=True)

    def _window(self) -> Any:
        from annotation_explorer_polars.windows.explore_window import ExploreWindow

        ann_df, meta_df = build_synthetic_run()
        win = ExploreWindow(
            ann_df,
            meta_df,
            ref_df=pd.DataFrame(),
            marker_hits_df=pd.DataFrame(),
            output_dir=Path(self.output_dir),
        )
        self.windows.append(win)
        return win

    def _write_old_keys(self, *, visible: bool, width: int, mode: str) -> bytes:
        geometry_source = QtWidgets.QWidget()
        geometry_source.setGeometry(20, 30, 711, 533)
        geometry = geometry_source.saveGeometry()
        geometry_source.deleteLater()
        settings = QtCore.QSettings()
        settings.beginGroup(OLD_GROUP)
        settings.setValue("geometry", geometry)
        settings.setValue("filtersVisible", visible)
        settings.setValue("filtersWidth", width)
        settings.setValue("viewMode", mode)
        settings.endGroup()
        settings.sync()
        return bytes(geometry)

    @staticmethod
    def _groups(win: Any) -> List[List[str]]:
        return sorted(sorted(g) for g in win._panel_manager.host.groups())

    def test_old_keys_are_ignored(self) -> None:
        # Phase 10: no migration any more; the window starts as a first run
        # (Classic, open rail, Contig) and never touches the old group.
        self._write_old_keys(visible=False, width=333, mode="Global")
        win = self._window()
        self.assertFalse(win._geometry_restored)
        self.assertFalse(win._filters_panel.isHidden())
        self.assertTrue(win.filters_action.isChecked())
        self.assertEqual(win.view_mode, "Contig")
        self.assertEqual(self._groups(win), CLASSIC_GROUPS)
        settings = QtCore.QSettings()
        self.assertFalse(settings.contains(f"{OLD_GROUP}/migrated"))
        win.close()
        _spin(self.app, 0.05)
        settings.sync()
        # Closing saves v2 only; the old keys stay as they were.
        self.assertEqual(settings.value(f"{OLD_GROUP}/viewMode"), "Global")
        self.assertEqual(int(settings.value(f"{OLD_GROUP}/filtersWidth")), 333)
        self.assertFalse(settings.contains(f"{OLD_GROUP}/workspaceSplitter"))
        self.assertFalse(settings.contains(f"{OLD_GROUP}/migrated"))
        saved = LayoutStore().load_window(run_key(Path(self.output_dir)))
        self.assertIsNotNone(saved)
        self.assertFalse(saved.rail_collapsed)

    def test_saved_v2_layout_wins_over_old_keys(self) -> None:
        self._write_old_keys(visible=False, width=333, mode="Global")
        LayoutStore().save_window(
            run_key(Path(self.output_dir)),
            WindowLayout(rail_width=260, rail_collapsed=False, dock_state=b"", panels=[]),
        )
        win = self._window()
        self.assertFalse(win._filters_panel.isHidden())
        self.assertEqual(win._filters_last_width, 260)
        self.assertEqual(win.view_mode, "Contig")
        # Empty dock state: Classic.
        self.assertEqual(self._groups(win), CLASSIC_GROUPS)
        # The old keys were not touched or marked.
        self.assertFalse(QtCore.QSettings().contains(f"{OLD_GROUP}/migrated"))



def _xml(state: bytes) -> str:
    return bytes(QtCore.qUncompress(QtCore.QByteArray(state))).decode("utf-8")


def _v1_xml(areas: str) -> bytes:
    """An uncompressed version 1 dock state with these areas."""
    return (
        '<?xml version="1.0" encoding="UTF-8"?><QtAdvancedDockingSystem Version="1" UserVersion="0" '
        f'Containers="1"><Container Floating="0">{areas}</Container></QtAdvancedDockingSystem>'
    ).encode("utf-8")


class FormatMigrationTests(unittest.TestCase):
    """migrate_layout on hand-made version 1 layouts (no window)."""

    @staticmethod
    def _widgets(state: bytes) -> List[tuple]:
        import xml.etree.ElementTree as ET

        data = state if state.lstrip().startswith(b"<") else bytes(QtCore.qUncompress(QtCore.QByteArray(state)))
        root = ET.fromstring(data)
        return [(w.get("Name"), w.get("Closed")) for w in root.iter("Widget")]

    def test_each_retired_kind_maps_to_its_new_kind(self) -> None:
        for old, new in (
            ("neopan", "phylogeny"),
            ("pan_tree", "phylogeny"),
            ("neotree", "phylogeny"),
            ("gene_presence", "gene_presence_absence"),
        ):
            with self.subTest(old=old):
                state = _v1_xml(f'<Area Tabs="1" Current="{old}"><Widget Name="{old}" Closed="0"/></Area>')
                layout = migrate_layout(
                    WindowLayout(dock_state=state, panels=[{"uid": old, "kind": old, "params": {}}], version=1)
                )
                self.assertEqual(layout.version, FORMAT_VERSION)
                self.assertEqual(self._widgets(layout.dock_state), [(new, "0")])
                self.assertIn(f'Current="{new}"', layout.dock_state.decode("utf-8"))
                self.assertEqual(layout.panels[0], {"uid": new, "kind": new, "params": {}})

    def test_a_compressed_state_stays_compressed_and_suffixed_uids_keep_their_suffix(self) -> None:
        xml = _v1_xml(
            '<Area Tabs="2" Current="pan_tree-91c0"><Widget Name="geo" Closed="0"/>'
            '<Widget Name="pan_tree-91c0" Closed="0"/></Area>'
            '<Area Tabs="1" Current="gene_presence-ab12"><Widget Name="gene_presence-ab12" Closed="0"/></Area>'
        )
        state = bytes(QtCore.qCompress(QtCore.QByteArray(xml), 9))
        panels = [
            {"uid": "pan_tree-91c0", "kind": "pan_tree", "params": {"layout": "rectangular", "branch_scale": "topology", "title": "Mine"}},
            {"uid": "gene_presence-ab12", "kind": "gene_presence", "params": {}},
            {"uid": "geo", "kind": "geo", "params": {"title": "World"}},
        ]
        layout = migrate_layout(WindowLayout(dock_state=state, panels=panels, version=1))
        self.assertFalse(layout.dock_state.startswith(b"<"))
        text = _xml(layout.dock_state)
        self.assertIn('Name="phylogeny-91c0"', text)
        self.assertIn('Name="gene_presence_absence-ab12"', text)
        self.assertIn('Current="phylogeny-91c0"', text)
        self.assertNotIn("pan_tree", text)
        by_uid = {p["uid"]: p for p in layout.panels}
        # The Pan tree's display options and the custom title carry over.
        self.assertEqual(
            by_uid["phylogeny-91c0"],
            {"uid": "phylogeny-91c0", "kind": "phylogeny", "params": {"layout": "rectangular", "branch_scale": "topology", "title": "Mine"}},
        )
        self.assertEqual(by_uid["gene_presence_absence-ab12"]["kind"], "gene_presence_absence")
        self.assertEqual(by_uid["geo"], panels[2])

    def test_neopan_and_neotree_give_one_open_phylogeny_and_neopan_brings_its_heatmap(self) -> None:
        xml = _v1_xml(
            '<Splitter Orientation="-" Count="2"><Area Tabs="2" Current="neotree">'
            '<Widget Name="sample_tree" Closed="0"/><Widget Name="neotree" Closed="0"/></Area>'
            '<Area Tabs="2" Current="neopan"><Widget Name="geo" Closed="0"/><Widget Name="neopan" Closed="0"/></Area>'
            "<Sizes>300 300 </Sizes></Splitter>"
        )
        panels = [
            {"uid": "neotree", "kind": "neotree", "params": {"pane_share": 0.4}},
            {"uid": "neopan", "kind": "neopan", "params": {"title": "Pan views"}},
        ]
        layout = migrate_layout(WindowLayout(dock_state=xml, panels=panels, version=1))
        widgets = self._widgets(layout.dock_state)
        names = [name for name, _closed in widgets]
        self.assertEqual(names.count("phylogeny"), 1)
        # neoTree came first and was open: it keeps the uid; neoPan's dock is
        # a closed extra instance (dropped by the restore).
        self.assertIn(("phylogeny", "0"), widgets)
        extra = [name for name, closed in widgets if name.startswith("phylogeny-")]
        self.assertEqual(len(extra), 1)
        self.assertIn((extra[0], "1"), widgets)
        text = layout.dock_state.decode("utf-8")
        self.assertIn('Current="phylogeny"', text)
        self.assertIn('Current="geo"', text)  # the area's current dock was closed
        by_uid = {p["uid"]: p for p in layout.panels}
        self.assertEqual(by_uid["phylogeny"]["params"], {})  # neoTree's options are dropped
        self.assertEqual(by_uid[extra[0]], {"uid": extra[0], "kind": "phylogeny", "params": {"title": "Pan views"}})
        # The open neoPan showed a heatmap: Gene presence/absence opens right
        # of the open Phylogeny.
        self.assertEqual(
            by_uid["gene_presence_absence"],
            {
                "uid": "gene_presence_absence",
                "kind": "gene_presence_absence",
                "params": {},
                "place": {"where": "right", "relative_to": "phylogeny"},
            },
        )

    def test_a_closed_neopan_or_an_open_gene_presence_adds_no_heatmap(self) -> None:
        closed = _v1_xml('<Area Tabs="2" Current="geo"><Widget Name="geo" Closed="0"/><Widget Name="neopan" Closed="1"/></Area>')
        layout = migrate_layout(WindowLayout(dock_state=closed, panels=[], version=1))
        self.assertNotIn("gene_presence_absence", [p["uid"] for p in layout.panels])
        both = _v1_xml(
            '<Area Tabs="1" Current="neopan"><Widget Name="neopan" Closed="0"/></Area>'
            '<Area Tabs="1" Current="gene_presence"><Widget Name="gene_presence" Closed="0"/></Area>'
        )
        layout = migrate_layout(WindowLayout(dock_state=both, panels=[], version=1))
        self.assertFalse(any("place" in p for p in layout.panels))

    def test_other_versions_and_layouts_without_retired_kinds_are_untouched(self) -> None:
        current = WindowLayout(dock_state=b"\x00state", panels=[{"uid": "neopan", "kind": "neopan", "params": {}}], version=FORMAT_VERSION)
        self.assertIs(migrate_layout(current), current)
        newer = WindowLayout(dock_state=b"x", panels=[], version=FORMAT_VERSION + 1)
        self.assertIs(migrate_layout(newer), newer)
        plain = _v1_xml('<Area Tabs="1" Current="geo"><Widget Name="geo" Closed="0"/></Area>')
        layout = migrate_layout(WindowLayout(dock_state=plain, panels=[{"uid": "geo", "kind": "geo", "params": {}}], version=1))
        self.assertEqual(layout.dock_state, plain)
        self.assertEqual(layout.version, FORMAT_VERSION)
        # A state that cannot be read is kept as it is (the restore decides).
        broken = migrate_layout(WindowLayout(dock_state=b"not a dock state", panels=[], version=1))
        self.assertEqual(broken.dock_state, b"not a dock state")

    def test_the_last_window_panel_list_is_renamed_without_duplicates(self) -> None:
        panels = [
            {"uid": "pan_tree", "kind": "pan_tree", "params": {"branch_scale": "topology"}},
            {"uid": "neopan", "kind": "neopan", "params": {}},
            {"uid": "gene_presence", "kind": "gene_presence", "params": {}},
            {"uid": "annotations", "kind": "annotations", "params": {"view_mode": "Global"}},
        ]
        layout = migrate_layout(WindowLayout(dock_state=b"", panels=panels, version=1))
        self.assertEqual(
            layout.panels,
            [
                {"uid": "phylogeny", "kind": "phylogeny", "params": {"branch_scale": "topology"}},
                {"uid": "gene_presence_absence", "kind": "gene_presence_absence", "params": {}},
                {"uid": "annotations", "kind": "annotations", "params": {"view_mode": "Global"}},
            ],
        )
        self.assertEqual(layout.dock_state, b"")


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class SavedV1LayoutTests(unittest.TestCase):
    """Real windows restore version 1 layouts, sessions and user layouts that
    name neoPan, the Pan tree, neoTree and Gene presence."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        from annotation_explorer_polars.theme import apply_application_theme

        apply_application_theme(cls.app)
        import annotation_explorer_polars.geo_experimental_views as gev

        cls._saved_webengine_cache = gev._WEBENGINE_RUNTIME_CACHE
        gev._WEBENGINE_RUNTIME_CACHE = (False, "disabled in tests")

    @classmethod
    def tearDownClass(cls) -> None:
        import annotation_explorer_polars.geo_experimental_views as gev

        gev._WEBENGINE_RUNTIME_CACHE = cls._saved_webengine_cache

    def setUp(self) -> None:
        _clear_settings()
        self.output_dir = tempfile.mkdtemp(prefix="bactolomew_v1_run_")
        write_pan_outputs(self.output_dir)
        self.windows: List[Any] = []

    def tearDown(self) -> None:
        for win in self.windows:
            win.close()
            win.deleteLater()
        _spin(self.app, 0.05)
        _clear_settings()
        shutil.rmtree(self.output_dir, ignore_errors=True)

    @staticmethod
    def _old_state(arrange: Any) -> bytes:
        """A real version 1 dock state: ADS docks with the old names."""
        from annotation_explorer_polars.docking.dock_host import DockHost, Placement

        host = DockHost()
        for uid in ("sample_tree", "annotations", "track", "geo", "neopan", "pan_tree", "neotree", "gene_presence"):
            host.register_view(uid, uid, QtWidgets.QLabel(uid))
        arrange(host, Placement)
        host.resize(1200, 800)
        host.show()
        QtWidgets.QApplication.processEvents()
        state = host.save_layout()
        host.deleteLater()
        return state

    def _window(self, output_dir: Any = None) -> Any:
        from annotation_explorer_polars.windows.explore_window import ExploreWindow

        ann_df, meta_df = build_synthetic_run()
        win = ExploreWindow(
            ann_df,
            meta_df,
            ref_df=pd.DataFrame(),
            marker_hits_df=pd.DataFrame(),
            output_dir=Path(output_dir or self.output_dir),
        )
        self.windows.append(win)
        win.show()
        _spin(self.app, 0.5)
        return win

    @staticmethod
    def _groups(win: Any) -> List[List[str]]:
        return sorted(sorted(g) for g in win._panel_manager.host.groups())

    def _classic_v1(self) -> WindowLayout:
        def arrange(host: Any, Placement: Any) -> None:
            host.place_view("sample_tree", Placement("left", outer=True))
            host.place_view("annotations", Placement("right", relative_to="sample_tree"))
            host.place_view("track", Placement("bottom", relative_to="annotations"))
            host.place_view("geo", Placement("bottom", outer=True))
            host.place_view("neopan", Placement("tab", relative_to="geo"))

        panels = [
            {"uid": uid, "kind": uid, "params": {}}
            for uid in ("sample_tree", "annotations", "track", "geo", "neopan", "pan_tree", "neotree", "gene_presence")
        ]
        return WindowLayout(dock_state=self._old_state(arrange), panels=panels, version=1)

    def test_a_v1_classic_window_slot_restores_with_phylogeny_and_its_heatmap(self) -> None:
        key = run_key(Path(self.output_dir))
        LayoutStore().save_window(key, self._classic_v1())
        settings = QtCore.QSettings()
        self.assertEqual(int(settings.value(f"explore/v2/run/{key}/windows/1/version")), 1)
        win = self._window()
        manager = win._panel_manager
        self.assertTrue(manager.is_open("phylogeny"))
        self.assertTrue(manager.is_open("gene_presence_absence"))
        self.assertEqual(manager.host.neighbor_group("phylogeny", "right"), "gene_presence_absence")
        self.assertIn(["geo", "phylogeny"], self._groups(win))
        for uid in ("sample_tree", "annotations", "track", "geo"):
            self.assertTrue(manager.is_open(uid), uid)
        # Nothing of the old kinds is left, and no closed duplicate either.
        self.assertEqual(manager.uids("phylogeny"), ["phylogeny"])
        self.assertEqual(manager.host.floating_count(), 0)
        # Loading wrote nothing: the slot is still version 1 until the next save.
        settings.sync()
        self.assertEqual(int(settings.value(f"explore/v2/run/{key}/windows/1/version")), 1)
        win._save_explore_layout()
        self.assertEqual(LayoutStore().load_window(key).version, FORMAT_VERSION)
        saved = LayoutStore().load_window(key)
        self.assertNotIn(b"neopan", bytes(QtCore.qUncompress(QtCore.QByteArray(saved.dock_state))))

    def test_a_v1_session_with_pan_tree_neotree_and_gene_presence_reopens(self) -> None:
        from annotation_explorer_polars.model import RunData
        from annotation_explorer_polars.shell.app_controller import AppController

        def arrange(host: Any, Placement: Any) -> None:
            host.place_view("pan_tree", Placement("left", outer=True))
            host.place_view("neotree", Placement("tab", relative_to="pan_tree"))
            host.place_view("gene_presence", Placement("right", relative_to="pan_tree"))
            host.place_view("geo", Placement("bottom", outer=True))

        key = run_key(Path(self.output_dir))
        store = LayoutStore()
        store.save_window(key, self._classic_v1(), 1)
        store.save_window(
            key,
            WindowLayout(
                dock_state=self._old_state(arrange),
                panels=[
                    {"uid": "pan_tree", "kind": "pan_tree", "params": {"branch_scale": "topology"}},
                    {"uid": "neotree", "kind": "neotree", "params": {}},
                    {"uid": "gene_presence", "kind": "gene_presence", "params": {}},
                ],
                version=1,
            ),
            2,
        )
        store.save_session(key, [1, 2])
        ann_df, meta_df = build_synthetic_run()
        run = RunData.build(
            ann_df, meta_df, ref_df=pd.DataFrame(), marker_hits_df=pd.DataFrame(), output_dir=Path(self.output_dir)
        )
        controller = AppController(run)
        controller.open_first_window()
        windows = controller.windows.explore_windows()
        self.windows += windows
        _spin(self.app, 0.8)
        self.assertEqual(len(windows), 2)
        second = windows[1]
        manager = second._panel_manager
        self.assertTrue(manager.is_open("phylogeny"))
        self.assertTrue(manager.is_open("gene_presence_absence"))
        self.assertTrue(manager.is_open("geo"))
        self.assertEqual(manager.uids("phylogeny"), ["phylogeny"])
        self.assertEqual(manager.host.neighbor_group("phylogeny", "right"), "gene_presence_absence")
        self.assertEqual(manager.frame("phylogeny").save_state().get("branch_scale"), "topology")
        # The first window (Classic with neoPan) has its tree and heatmap too.
        first = windows[0]._panel_manager
        self.assertTrue(first.is_open("phylogeny") and first.is_open("gene_presence_absence"))

    def test_v1_user_and_default_layouts_apply(self) -> None:
        from annotation_explorer_polars.shell.layout_store import DEFAULT_GROUP

        store = LayoutStore()
        store.save_user_layout("Old pan", self._classic_v1())
        store.save_default(self._classic_v1())
        settings = QtCore.QSettings()
        settings.setValue(f"{DEFAULT_GROUP}/version", 1)
        settings.sync()
        self.assertEqual(store.load_user_layout("Old pan").version, FORMAT_VERSION)
        # A run with no saved layout starts with the (version 1) default.
        win = self._window()
        manager = win._panel_manager
        self.assertTrue(manager.is_open("phylogeny") and manager.is_open("gene_presence_absence"))
        manager.apply_preset("geo_table")
        _spin(self.app, 0.2)
        self.assertFalse(manager.is_open("phylogeny"))
        self.assertTrue(win._apply_layout(name="Old pan"))
        _spin(self.app, 0.3)
        self.assertTrue(manager.is_open("phylogeny") and manager.is_open("gene_presence_absence"))
        self.assertEqual(manager.host.floating_count(), 0)


if __name__ == "__main__":
    unittest.main()
