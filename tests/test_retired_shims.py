"""Phase 10 of the window plan (8.11): the shims are retired.

- The explore window is ``windows.explore_window.ExploreWindow``; the
  ``AnnotationExplorer`` alias, ``explorer_window.py``, the
  ``PolarsQueriesMixin`` shim and the window's delegating properties are gone.
- The legacy layout flag and the old QSettings reads are gone, so PyQt6Ads is
  required: without it the app says so at startup, with the install command,
  instead of crashing later.
- Hover reaches every visible heavy panel within the 30 ms of plan 8.12 while
  a mouse sweep still redraws a heavy panel at most once per coalescing
  interval.
"""
from __future__ import annotations

import importlib
import os
import time
import unittest
from typing import List
from unittest import mock

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
    from tests._timing import within_budget
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings
    from _timing import within_budget

isolate_qt_settings()

from PyQt6 import QtCore, QtWidgets

import annotation_explorer_polars
from annotation_explorer_polars import docking
from annotation_explorer_polars.docking import dock_host
from annotation_explorer_polars.session import bus as bus_events
from annotation_explorer_polars.session import explore_context
from annotation_explorer_polars.session.bus import Delivery
from annotation_explorer_polars.session.explore_context import ExploreContext
from annotation_explorer_polars.windows.explore_window import ExploreWindow


def _app() -> QtWidgets.QApplication:
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


def _spin(seconds: float) -> None:
    app = _app()
    end = time.monotonic() + seconds
    while True:
        app.processEvents(QtCore.QEventLoop.ProcessEventsFlag.AllEvents, 5)
        if time.monotonic() >= end:
            break
        time.sleep(0.001)


class RetiredNamesTests(unittest.TestCase):
    def test_the_window_has_its_new_name_only(self) -> None:
        self.assertIs(annotation_explorer_polars.ExploreWindow, ExploreWindow)
        with self.assertRaises(AttributeError):
            annotation_explorer_polars.AnnotationExplorer  # noqa: B018
        for module in ("explorer_window", "polars_queries", "hover_sync"):
            with self.assertRaises(ModuleNotFoundError, msg=module):
                importlib.import_module(f"annotation_explorer_polars.{module}")

    def test_the_window_has_no_delegating_properties(self) -> None:
        # RunData fields, selection and result views, the refresh views and
        # the legacy layout's widgets: read them from run / ctx / panels.
        retired = (
            "ann_df", "meta_df", "ref_df", "marker_hits_df", "_ann_pl", "_meta_pl", "annotation_db_columns",
            "meta_columns", "predictor_options", "output_dir", "run_id", "_snippy_reference_sample",
            "current_sample", "current_contig", "current_contigs", "current_contig_label",
            "ann_filters", "_base_allowed_samples", "_filtered_ann_pl", "_applied_result",
            "_sample_subset_filter_samples", "_sample_subset_filter_label", "_sample_subset_filter_kind",
            "_table_scroll_top_on_show", "_refresh_timer", "_refresh_running", "_refresh_pending_snapshot",
            "_first_refresh_done", "_clean_text", "_collect_sample_list", "_subset_annotations_pd",
            "_empty_ann_public_pl", "_build_gene_workspace", "_save_ui_state", "_restore_ui_state",
            "_ensure_phandango_plot_dialog", "_ensure_cgmlst_force_plot_dialog", "_open_document_window",
            "_on_bottom_tab_changed", "_apply_pan_subset_filter_request", "_on_pan_clade_filter_requested",
            "_on_pan_cluster_filter_requested", "_on_pan_plasmid_filter_requested", "_set_view_mode_silently",
            "_build_table_view_switch", "_sync_table_view_box", "_on_view_mode_changed",
        )
        present = [name for name in retired if hasattr(ExploreWindow, name)]
        self.assertEqual(present, [])

    def test_neotree_is_gone(self) -> None:
        # Feedback 2026-10-07: the neoTree panel and its code were removed;
        # the tree model and the clade reroot live in ``phylo``.
        for module in (
            "neotree",
            "neotree.canvas",
            "panels.neotree",
            "model.gene_families",
            "model.recombination",
        ):
            with self.assertRaises(ModuleNotFoundError, msg=module):
                importlib.import_module(f"annotation_explorer_polars.{module}")
        from annotation_explorer_polars.panels import base

        self.assertFalse(hasattr(base, "NEOTREE"))
        for name in (
            "_create_neotree_panel", "_create_neo_tree_tab", "_show_phandango_style_plot",
            "_collect_neo_tree_trait_catalog", "_neo_tree_default_selected_traits",
            "_set_neo_tree_trait_present_color", "_pan_tabs",
        ):
            self.assertFalse(hasattr(ExploreWindow, name), name)
        # neoTree's leftovers: the tree's side columns and the tree source
        # switch had no caller left.
        from annotation_explorer_polars.panaroo_core_accessory import PanCoreAccessoryTab
        from annotation_explorer_polars.tree_view import RaxmlTreeView

        for name in ("set_leaf_side_metadata", "set_leaf_side_metadata_tracks", "_draw_side_metadata_tracks"):
            self.assertFalse(hasattr(RaxmlTreeView, name), name)
        for name in ("set_tree_source", "tree_source_key", "tree_source_label", "tree_source_path"):
            self.assertFalse(hasattr(PanCoreAccessoryTab, name), name)

    def test_neopan_pan_tree_and_gene_presence_are_gone(self) -> None:
        # Feedback 2026-10-07 (d): Phylogeny and Gene presence/absence replace
        # neoPan, the Pan tree and Gene presence.
        with self.assertRaises(ModuleNotFoundError):
            importlib.import_module("annotation_explorer_polars.panels.neopan")
        from annotation_explorer_polars.panels import base, pan
        from annotation_explorer_polars.panaroo_core_accessory import PanCoreAccessoryTab

        for name in ("NEOPAN", "PAN_TREE", "GENE_PRESENCE"):
            self.assertFalse(hasattr(base, name), name)
        for name in ("NeoPanPanel", "PanTreePanel", "GenePresencePanel"):
            self.assertFalse(hasattr(pan, name), name)
        for name in ("_create_neopan_panel", "_create_neopan_tab"):
            self.assertFalse(hasattr(ExploreWindow, name), name)
        self.assertEqual(base.PAN_VIEW_KINDS[:2], ("phylogeny", "gene_presence_absence"))
        for name in ("_create_selector_panel", "all_geo_views", "view_for_key", "viewSamplesHovered"):
            self.assertFalse(hasattr(PanCoreAccessoryTab, name), name)

    def test_no_legacy_layout_flag(self) -> None:
        self.assertFalse(hasattr(docking, "LEGACY_LAYOUT_ENV"))
        self.assertFalse(hasattr(docking, "legacy_layout_requested"))


class MissingDockingLibraryTests(unittest.TestCase):
    """Plan 6.3: without the legacy fallback a missing PyQt6Ads must give a
    clear message with the install hint, not a crash."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.app = _app()

    def _without_ads(self):
        return mock.patch.multiple(dock_host, ads=None, _ADS_IMPORT_ERROR="ModuleNotFoundError: No module named 'PyQt6Ads'")

    def test_the_message_names_the_package_and_the_install_command(self) -> None:
        with self._without_ads():
            self.assertFalse(docking.docking_available())
            message = docking.docking_unavailable_message()
            with self.assertRaises(docking.DockingUnavailableError) as raised:
                docking.require_docking()
        self.assertIn("PyQt6Ads", message)
        self.assertIn("pip install 'PyQt6Ads>=5.1,<6'", message)
        self.assertIn("requirements.txt", message)
        self.assertIn("No module named 'PyQt6Ads'", message)
        self.assertEqual(str(raised.exception), message)
        self.assertNotIn(chr(0x2014), message)  # em dash
        self.assertNotIn(chr(0x2013), message)  # en dash

    def test_the_install_command_survives_a_path_with_spaces(self) -> None:
        import shlex

        with self._without_ads(), mock.patch("sys.executable", "/opt/my envs/py 3/bin/python"):
            hint = docking.docking_install_hint()
        command = next(line.strip() for line in hint.splitlines() if "pip install" in line)
        self.assertEqual(
            shlex.split(command), ["/opt/my envs/py 3/bin/python", "-m", "pip", "install", "PyQt6Ads>=5.1,<6"]
        )

    def test_an_explore_window_refuses_to_build_with_the_message(self) -> None:
        with self._without_ads():
            with self.assertRaises(docking.DockingUnavailableError) as raised:
                ExploreWindow()
            with self.assertRaises(docking.DockingUnavailableError):
                dock_host.DockHost()
        self.assertIn("pip install", str(raised.exception))
        self.assertEqual([w for w in QtWidgets.QApplication.topLevelWidgets() if isinstance(w, ExploreWindow)], [])

    def test_the_app_stops_at_startup_with_a_message_box(self) -> None:
        from annotation_explorer_polars import app

        shown: List[str] = []

        def fake_exec(box: QtWidgets.QMessageBox) -> int:
            shown.append(box.text() + "\n" + box.informativeText())
            return 0

        with self._without_ads(), mock.patch.object(QtWidgets.QMessageBox, "exec", fake_exec), \
                mock.patch("sys.stderr") as stderr:
            self.assertFalse(app._docking_ready())
        self.assertEqual(len(shown), 1)
        self.assertIn("PyQt6Ads", shown[0])
        self.assertIn("pip install", shown[0])
        # The headline is not repeated in the details.
        self.assertEqual(shown[0].count(docking.DOCKING_UNAVAILABLE_HEADLINE), 1)
        self.assertTrue(stderr.write.called)
        # With the library the check passes silently.
        with mock.patch.object(QtWidgets.QMessageBox, "exec", fake_exec):
            self.assertEqual(app._docking_ready(), docking.docking_available())
        self.assertEqual(len(shown), 1 if docking.docking_available() else 2)

    def test_main_returns_before_any_dialog_without_the_library(self) -> None:
        from annotation_explorer_polars import app

        with mock.patch.object(app, "_docking_ready", return_value=False), \
                mock.patch.object(app, "_ensure_startup_ready") as ready, \
                mock.patch.object(app, "StartupDialog") as startup, \
                mock.patch.object(app.QtWidgets, "QApplication", return_value=self.app), \
                mock.patch.object(app, "apply_application_theme"):
            status = app.main()
        # A launcher can tell a start that failed from a normal quit.
        self.assertEqual(status, 1)
        ready.assert_not_called()
        startup.assert_not_called()


_NO_ADS_LAUNCH = r"""
import runpy, sys, tempfile
sys.modules["PyQt6Ads"] = None  # as if PyQt6Ads were not installed
from PyQt6 import QtCore, QtWidgets
fmt, scope = QtCore.QSettings.Format, QtCore.QSettings.Scope
QtCore.QSettings.setDefaultFormat(fmt.IniFormat)
settings_dir = tempfile.mkdtemp(prefix="bactolomew_noads_settings_")
for each_format in (fmt.IniFormat, fmt.NativeFormat):
    for each_scope in (scope.UserScope, scope.SystemScope):
        QtCore.QSettings.setPath(each_format, each_scope, settings_dir)

def fake_exec(box):
    print("BOX", box.windowTitle(), "|", box.text(), "|", box.informativeText().replace(chr(10), " "), flush=True)
    return 0

QtWidgets.QMessageBox.exec = fake_exec
QtWidgets.QDialog.exec = lambda self: sys.exit("a dialog opened: " + type(self).__name__)
sys.argv = [sys.argv[1]]
runpy.run_path(sys.argv[0], run_name="__main__")
"""


class MissingDockingLibraryLaunchTests(unittest.TestCase):
    """The real entry point with PyQt6Ads blocked, in a fresh interpreter: no
    module imports PyQt6Ads at import time, and the start fails cleanly."""

    def test_the_launcher_shows_the_box_and_exits_with_status_1(self) -> None:
        import subprocess
        import sys
        from pathlib import Path

        src = Path(annotation_explorer_polars.__file__).resolve().parent.parent
        env = dict(os.environ, QT_QPA_PLATFORM="offscreen", PYTHONPATH=str(src), PYTHONDONTWRITEBYTECODE="1")
        done = subprocess.run(
            [sys.executable, "-c", _NO_ADS_LAUNCH, str(src / "annotation_explorer_qt.py")],
            env=env,
            capture_output=True,
            text=True,
            timeout=120,
        )
        boxes = [line for line in done.stdout.splitlines() if line.startswith("BOX ")]
        self.assertEqual(done.returncode, 1, done.stderr[-2000:])
        self.assertEqual(len(boxes), 1, done.stdout)
        title, text, details = (part.strip() for part in boxes[0][4:].split("|"))
        self.assertEqual(title, "PyQt6Ads is missing")
        self.assertIn("PyQt6Ads", text)
        self.assertIn("-m pip install 'PyQt6Ads>=5.1,<6'", details)
        self.assertNotIn(text, details)  # the box says it once
        self.assertIn("pip install", done.stderr)
        self.assertNotIn("Traceback", done.stderr)


class HoverBudgetTests(unittest.TestCase):
    """Plan 8.12: hover reaches every visible heavy panel within 30 ms, and a
    sweep still redraws a heavy panel at most once per interval."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.app = _app()

    def setUp(self) -> None:
        # Windows an earlier test module left to deleteLater are deleted in
        # the first event loop that runs: delete them now, outside the timing.
        QtCore.QCoreApplication.sendPostedEvents(None, QtCore.QEvent.Type.DeferredDelete.value)
        self.app.processEvents()

    def test_the_coalescing_intervals_fit_the_budget(self) -> None:
        total = explore_context.HOVER_COALESCE_MS + bus_events.HEAVY_HOVER_COALESCE_MS
        # Room is left for the handlers and the event loop (Phase 10
        # measured about 2 ms on top on the reference run).
        self.assertLessEqual(total, 24)
        self.assertGreaterEqual(explore_context.HOVER_COALESCE_MS, 1)
        self.assertGreaterEqual(bus_events.HEAVY_HOVER_COALESCE_MS, 8)

    def _heavy_listeners(self, ctx: ExploreContext, count: int) -> List[List[float]]:
        got: List[List[float]] = []
        for i in range(count):
            times: List[float] = []
            got.append(times)
            ctx.bus.subscribe(
                f"heavy{i}",
                bus_events.HOVER_CHANGED,
                lambda _hover, _t=times: _t.append(time.perf_counter()),
                delivery=Delivery.COALESCED,
                coalesce_ms=bus_events.HEAVY_HOVER_COALESCE_MS,
            )
        return got

    def test_a_hover_reaches_four_heavy_panels_within_30_ms(self) -> None:
        # The bus alone (the two timers); tests/test_pan_panels.py times the
        # real Map, neoPan, Pan tree and Gene presence panels.
        ctx = ExploreContext(engine=None)
        self.addCleanup(ctx.deleteLater)
        got = self._heavy_listeners(ctx, 4)

        def median_ms() -> float:
            latencies = []
            for i in range(5):
                for times in got:
                    times.clear()
                start = time.perf_counter()
                ctx.hover([f"S{i}"], origin="sample_tree")
                loop = QtCore.QEventLoop()
                QtCore.QTimer.singleShot(200, loop.quit)
                timer = QtCore.QTimer()
                timer.setInterval(1)
                timer.timeout.connect(lambda: loop.quit() if all(got) else None)
                timer.start()
                loop.exec()
                timer.stop()
                self.assertTrue(all(got))
                latencies.append(max(t[0] for t in got) - start)
            # The median, so one slow pass of a busy test machine does not fail it.
            return sorted(latencies)[len(latencies) // 2] * 1000

        ok, seen = within_budget(median_ms, 30)
        self.assertTrue(ok, f"median hover latency {seen} ms")

    def test_a_sweep_redraws_a_heavy_panel_at_most_once_per_interval(self) -> None:
        ctx = ExploreContext(engine=None)
        self.addCleanup(ctx.deleteLater)
        got = self._heavy_listeners(ctx, 1)
        start = time.perf_counter()
        for i in range(60):  # one hover per ms for 60 ms
            ctx.hover([f"S{i}"], origin="sample_tree")
            _spin(0.001)
        elapsed_ms = (time.perf_counter() - start) * 1000
        _spin(0.1)
        deliveries = len(got[0])
        self.assertGreaterEqual(deliveries, 1)
        self.assertLessEqual(deliveries, int(elapsed_ms / bus_events.HEAVY_HOVER_COALESCE_MS) + 2)
        # The last hover always arrives.
        self.assertEqual(ctx.hover_state.samples, ("S59",))


    def test_the_offline_map_shows_hover_without_re_rendering(self) -> None:
        # The SVG fallback map rebuilt its whole document on every hover
        # (about 130 ms on the reference run); hover is an overlay now.
        from annotation_explorer_polars.geo_svg_map import SvgGeoView

        view = SvgGeoView()
        self.addCleanup(view.deleteLater)
        view.resize(600, 300)
        view.show()
        _spin(0.02)
        view.update_data([{"country": "Romania", "sample": "A"}, {"country": "China", "sample": "B"}], {}, 1)
        renderer = view._renderer
        view.scale(2.0, 2.0)
        zoom = view.transform().m11()
        with mock.patch.object(view, "_rebuild", wraps=view._rebuild) as rerender:
            view.highlight_samples(["B"])
        rerender.assert_not_called()  # structure, not wall-clock time
        self.assertEqual(view.highlighted_countries(), {"china"})
        self.assertFalse(view._highlight_item.path().isEmpty())
        self.assertIs(view._renderer, renderer)  # no new document
        self.assertEqual(view.transform().m11(), zoom)  # the user's zoom stays
        view.highlight_samples(["A", "B", "nobody"])
        self.assertEqual(view.highlighted_countries(), {"china", "romania"})
        two = view._highlight_item.path().boundingRect()
        self.assertGreater(two.width(), 0)
        view.highlight_samples([])
        self.assertTrue(view._highlight_item.path().isEmpty())
        # A new result keeps the hovered countries highlighted.
        view.highlight_samples(["A"])
        view.update_data([{"country": "Romania", "sample": "A"}], {}, 1)
        self.assertFalse(view._highlight_item.path().isEmpty())

if __name__ == "__main__":
    unittest.main()
