"""Neighbourhood and alignment panels in the real window (window plan Phase 6, 8.7).

- The neighbourhood comparison is a document panel of its explore window:
  opened beside its source (plan 3.6 rule 3), one per comparison, with an
  inline busy state while a gene card's rows are collected.
- Its rows follow the window's one sample set (D4): no metadata filters of
  its own; the rail's active filters are shown read only. Display options
  (molecule, sort, flank, color, labels, mirror, merge) stay and are its
  panel params.
- Hover goes both ways through the context; "Open sample(s) in table"
  narrows the window through ``set_subset``.
- Ctrl+W and Ctrl+F act only inside the panel (plan 3.7); a second Esc
  closes nothing.
- Document panels are not saved with the layout (3.11), survive Reset
  Layout, come back with Reopen Closed Panel and are listed in the palette
  only while open.
- The alignment result is a document panel too, its rows limited to the
  sample set.

- Without panel mode the body (``NeighbourhoodView``) is a standalone viewer
  with its own metadata filters (the ``IgvCollectionDialog`` window shim
  went in Phase 10; the tests build the body directly).

Offscreen on the synthetic run of the golden tests; the window tests are
skipped when PyQt6Ads is not installed.
"""
from __future__ import annotations

import os
import shutil
import tempfile
import threading
import time
import unittest
from typing import Any, List
from unittest import mock

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

from PyQt6 import QtCore, QtGui, QtTest, QtWidgets

from annotation_explorer_polars.docking import docking_available
from annotation_explorer_polars import igv_collection_view
from annotation_explorer_polars.igv_collection_view import NeighbourhoodView
from annotation_explorer_polars.model.filter_state import SampleSubset
from annotation_explorer_polars.panels.alignment import AlignmentPanel
from annotation_explorer_polars.panels.neighbourhood import NeighbourhoodPanel, filter_notes

try:
    from tests.test_dock_host import _app, _Explorer, _spin
except ImportError:
    from test_dock_host import _app, _Explorer, _spin

ALIGNMENT_TEXT = ">a|S1|1|g\nMKA\n>b|S2|1|g\nMKA\n>c|S5|1|g\nMRA\n"


def _rows(win: Any, annotation_id: str = "S1:1:1") -> List[dict]:
    from annotation_explorer_polars.annotation_detail_dialog import AnnotationDetailDialog

    row = win.run.ann_df.loc[win.run.ann_df["AnnotationID"] == annotation_id].iloc[0].to_dict()
    card = AnnotationDetailDialog(None, row, list(row), ann_df=win.run.ann_df)
    try:
        return card._build_context_collection_rows(window_genes=30)
    finally:
        card.deleteLater()


def _row_samples(panel: Any) -> List[str]:
    return sorted({str(r.get("sample")) for r in panel.view._visible_rows})


class _Case(unittest.TestCase):
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
        self.output_dir = tempfile.mkdtemp(prefix="bactolomew_neighbourhood_run_")
        self.windows: List[Any] = []
        self._clear()
        NeighbourhoodPanel.settings_panel_open = None

    def tearDown(self) -> None:
        for win in self.windows:
            win.close()
            win.deleteLater()
        _spin(self.app, 0.05)
        shutil.rmtree(self.output_dir, ignore_errors=True)
        self._clear()
        NeighbourhoodPanel.settings_panel_open = None

    @staticmethod
    def _clear() -> None:
        settings = QtCore.QSettings()
        settings.remove("explore")
        settings.remove("annotationExplorer")
        settings.sync()

    def _window(self) -> Any:
        win = _Explorer.build(self.output_dir)
        self.windows.append(win)
        _Explorer.show_and_settle(self.app, win, seconds=0.3)
        return win

    def _idle(self, win: Any, extra: float = 0.3) -> None:
        end = time.monotonic() + 10
        while time.monotonic() < end and (
            win.ctx.scheduler.is_running() or win.ctx.scheduler.debounce_timer.isActive() or win._panels.has_pending_delivery()
        ):
            _spin(self.app, 0.02)
        _spin(self.app, extra)

    def _wait(self, predicate: Any, timeout: float = 5.0) -> bool:
        end = time.monotonic() + timeout
        while time.monotonic() < end:
            if predicate():
                return True
            _spin(self.app, 0.02)
        return bool(predicate())

    def _open(self, win: Any, **kwargs: Any) -> NeighbourhoodPanel:
        kwargs.setdefault("rows", _rows(win))
        kwargs.setdefault("title", "Gene neighborhoods: lysR")
        kwargs.setdefault("source", "gene_card")
        uid = win._open_neighbourhood_view(**kwargs)
        self.assertIsNotNone(uid)
        _spin(self.app, 0.4)
        panel = win._panels.panel(uid)
        self.assertIsInstance(panel, NeighbourhoodPanel)
        return panel


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class NeighbourhoodPanelTests(_Case):
    def test_rows_follow_the_window_sample_set(self) -> None:
        win = self._window()
        panel = self._open(win)
        everyone = _row_samples(panel)
        self.assertGreater(len(everyone), 3)
        self.assertEqual(panel.view.scope_label.text(), "Every row is in this window's samples.")
        # An annotation filter narrows the window's sample set: the rows follow.
        win.gene_query.setText("blaOXA")
        self._idle(win)
        samples = set(win.ctx.sample_set())
        self.assertTrue(samples)
        self.assertLess(len(samples), 6)
        self.assertEqual(set(_row_samples(panel)), set(everyone) & samples)
        self.assertIn("are in this window's samples", panel.view.scope_label.text())
        self.assertEqual(panel.view._scope_notes, ['Gene: "blaOXA"'])
        # (In the Settings column, which a narrow panel starts folded.)
        self.assertFalse(panel.view.scope_chips.isHidden())
        # A subset from another panel narrows it too.
        win.gene_query.setText("")
        self._idle(win)
        self.assertEqual(_row_samples(panel), everyone)
        win.ctx.set_subset(SampleSubset(("S2",), label="Clade 2", kind="clade"), origin="phylogeny")
        self._idle(win)
        self.assertEqual(_row_samples(panel), ["S2"])
        self.assertEqual(panel.view._scope_notes, ["Clade 2"])
        # The reference scores stay anchored to the comparison's own reference.
        self.assertTrue(all("context_similarity" in row for row in panel.view._visible_rows))
        win.ctx.set_subset(None)
        self._idle(win)
        self.assertEqual(_row_samples(panel), everyone)

    def test_no_metadata_filters_of_its_own_but_every_display_option(self) -> None:
        win = self._window()
        panel = self._open(win, flank_sizes=(5, 10, 15, 30), initial_flank=15)
        view = panel.view
        self.assertTrue(view.is_panel_mode())
        self.assertFalse(view.add_filter_btn.isVisibleTo(view))
        self.assertFalse(view.add_filter_btn.isEnabled())
        # No search of its own: the window's search (the rail) sets its samples.
        self.assertFalse(view.sample_query.isVisibleTo(view))
        # Display options stay: molecule, sort, flank, color, labels, mirror, merge.
        for name in ("molecule_box", "sort_column_box", "sort_order_box", "flank_box", "color_by_box",
                     "labels_box", "normalize_orientation_check", "collapse_btn", "group_molecule_check"):
            self.assertTrue(getattr(view, name).isEnabled(), name)
        # They are panel params, restored as saved.
        view.flank_box.setCurrentIndex(view.flank_box.findData(5))
        view.color_by_box.setCurrentIndex(view.color_by_box.findData("family"))
        view.collapse_btn.setChecked(True)
        state = panel.save_state()
        self.assertEqual((state["flank"], state["color_by"], state["merge"]), (5, "family", True))
        other = self._open(win, key=("other",), flank_sizes=(5, 10, 15, 30), initial_flank=15)
        other.restore_state(state)
        self.assertEqual(other.view._context_flank, 5)
        self.assertEqual(other.view.view.color_mode, "family")
        self.assertTrue(other.view._collapsed)
        self.assertEqual(other.save_state(), state)

    def test_filter_notes_name_what_narrows_the_sample_set(self) -> None:
        from annotation_explorer_polars.model.filter_state import FilterState

        filters = FilterState(
            sample_query="ro", gene="bla", metadata={"Country": ["Romania", "Italy"]},
            subset=SampleSubset(("S1",), label="Clade 3"), hide_empty_contigs=True,
        )
        self.assertEqual(filter_notes(filters), ['Search: "ro"', "Country: Italy, Romania", "Clade 3", 'Gene: "bla"'])
        self.assertEqual(filter_notes(FilterState()), [])

    def test_hover_goes_both_ways_through_the_context(self) -> None:
        win = self._window()
        panel = self._open(win)
        canvas = panel.view.view
        # Out: a hovered row hovers its samples in the window, with the panel's uid.
        row = canvas.layout_result.rows[1]
        canvas._set_hover(None, row)
        win.ctx.flush_hover()
        self.assertEqual(win.ctx.hover_state.samples, tuple(row.row.samples))
        self.assertEqual(win.ctx.hover_state.origin, panel.uid)
        self.assertEqual(canvas.linked_samples, frozenset())  # no echo
        canvas._set_hover(None, None)
        win.ctx.flush_hover()
        self.assertEqual(win.ctx.hover_state.samples, ())
        # In: samples hovered elsewhere tint their rows (rail and plot).
        win.ctx.hover(["S2"], origin="sample_tree")
        win.ctx.flush_hover()
        self.assertTrue(self._wait(lambda: canvas.linked_samples == frozenset({"S2"})))
        linked = [r for r in canvas.layout_result.rows if canvas.row_is_linked(r)]
        self.assertTrue(linked)
        self.assertTrue(all("S2" in r.row.samples for r in linked))
        self.assertEqual(len(canvas._linked_row_items), len(linked))
        # Kept across a redraw (a flank change rebuilds the scene).
        panel.view._set_context_flank(5) if panel.view._flank_sizes else panel.view._apply_filters()
        self.assertEqual(len(canvas._linked_row_items), len(linked))
        win.ctx.clear_hover("sample_tree")
        win.ctx.flush_hover()
        _spin(self.app, 0.1)
        self.assertEqual(canvas._linked_row_items, [])
        # A row hovered when the panel closes is not left hovered.
        canvas._set_hover(None, canvas.layout_result.rows[0])
        win.ctx.flush_hover()
        self.assertTrue(win.ctx.hover_state.samples)
        win._panel_manager.close_panel(panel.uid)
        win.ctx.flush_hover()
        self.assertEqual(win.ctx.hover_state.samples, ())

    def test_no_tree_copy_and_a_phylogeny_panel_follows_the_rows(self) -> None:
        # The panel has no Phylogeny button or tree of its own; a Phylogeny
        # panel of the same window gets the hovered rows through the window.
        try:
            from tests._pan_fixture import write_pan_outputs
        except ImportError:
            from _pan_fixture import write_pan_outputs

        write_pan_outputs(self.output_dir)
        win = self._window()
        panel = self._open(win)
        view = panel.view
        self.assertFalse(hasattr(view, "phylogeny_btn"))
        self.assertFalse(hasattr(view, "panCopyCreated"))
        self.assertNotIn("pan_reference_sample", panel._view_options)
        manager = win._panel_manager
        self.assertTrue(manager.open_panel("phylogeny"))
        win._panels.set_ready()
        self.assertTrue(self._wait(lambda: manager.panel("phylogeny") is not None, timeout=15.0))
        phylogeny = manager.panel("phylogeny")
        hovers: List[Any] = []
        with mock.patch.object(phylogeny, "show_hover", side_effect=lambda hover: hovers.append(hover.samples)):
            canvas = view.view
            row = canvas.layout_result.rows[1]
            canvas._set_hover(None, row)
            win.ctx.flush_hover()
            self.assertTrue(self._wait(lambda: tuple(row.row.samples) in hovers))
            canvas._set_hover(None, None)
            win.ctx.flush_hover()
            self.assertTrue(self._wait(lambda: bool(hovers) and hovers[-1] == ()))

    def test_open_samples_in_table_narrows_the_window(self) -> None:
        win = self._window()
        panel = self._open(win)
        subsets: List[Any] = []
        win.ctx.subsetChanged.connect(lambda subset, origin: subsets.append((subset, origin)))
        texts: List[str] = []
        original = QtWidgets.QMenu.exec
        QtWidgets.QMenu.exec = lambda menu, *args: texts.extend(a.text() for a in menu.actions())
        try:
            panel.view._show_row_menu(dict(panel.view.view.layout_result.rows[0].row.source), QtCore.QPoint())
        finally:
            QtWidgets.QMenu.exec = original
        self.assertIn("Show only this sample in the window", texts)
        self.assertNotIn("Open sample in table", texts)
        panel.view._open_samples(["S3"])
        self._idle(win)
        self.assertEqual(subsets[-1], (SampleSubset(("S3",), label="Sample S3", kind="sample"), panel.uid))
        self.assertEqual(win.ctx.sample_set(), ("S3",))
        self.assertEqual(_row_samples(panel), ["S3"])

    def test_gene_double_click_opens_the_card_in_this_window(self) -> None:
        win = self._window()
        panel = self._open(win)
        row = win.run.ann_df.loc[win.run.ann_df["AnnotationID"] == "S2:1:1"].iloc[0].to_dict()
        panel.view.geneDoubleClicked.emit(row)
        cards = win.app_controller.windows.info_windows(win)
        self.assertEqual(len(cards), 1)
        self.assertEqual(cards[0]._row_data["AnnotationID"], "S2:1:1")

    def test_shortcuts_act_only_inside_the_panel(self) -> None:
        win = self._window()
        panel = self._open(win)
        other = self._open(win, key=("other",), title="Gene neighborhoods: other")
        shortcuts = {s.key().toString(): s for s in panel.findChildren(QtGui.QShortcut)}
        self.assertEqual(sorted(shortcuts), ["Ctrl+F", "Ctrl+W"])
        for shortcut in shortcuts.values():
            self.assertEqual(shortcut.context(), QtCore.Qt.ShortcutContext.WidgetWithChildrenShortcut)
        # The window itself has no shortcut (D8): only the panels have theirs.
        outside = [s for s in win.findChildren(QtGui.QShortcut) if not panel.isAncestorOf(s.parent()) and not other.isAncestorOf(s.parent())]
        self.assertEqual(outside, [])
        # Ctrl+F opens the settings and focuses this panel's search.
        shortcuts["Ctrl+F"].activated.emit()
        self.assertTrue(panel.view.settings_open())
        self.assertFalse(other.view.settings_open())
        # Ctrl+W closes only this panel (through the window's docking area).
        manager = win._panel_manager
        shortcuts["Ctrl+W"].activated.emit()
        _spin(self.app, 0.1)
        self.assertIsNone(manager.frame(panel.uid))
        self.assertTrue(manager.is_open(other.uid))
        self.assertTrue(win.isVisible())

    def test_escape_clears_the_pin_and_closes_nothing(self) -> None:
        win = self._window()
        panel = self._open(win)
        view = panel.view
        canvas = view.view
        canvas._select_gene(canvas.layout_result.rows[0].genes[0])
        self.assertTrue(view.pin_chip.isVisibleTo(view))
        QtTest.QTest.keyClick(canvas.graphics_view, QtCore.Qt.Key.Key_Escape)
        self.assertFalse(canvas.has_selection())
        self.assertFalse(view.synteny_detail_panel.isVisibleTo(view))
        QtTest.QTest.keyClick(canvas.graphics_view, QtCore.Qt.Key.Key_Escape)
        _spin(self.app, 0.05)
        self.assertTrue(win._panel_manager.is_open(panel.uid))
        self.assertTrue(win.isVisible())

    def test_one_panel_per_comparison_beside_its_source(self) -> None:
        win = self._window()
        manager = win._panel_manager
        manager.host.focus_view("sample_tree")
        first = self._open(win, key=("gene", "a"))
        # Rule 3: a tab of the group right of the source (Annotations).
        self.assertIn(first.uid, manager.host.group_uids("annotations"))
        again = win._open_neighbourhood_view(rows=_rows(win), title="x", key=("gene", "a"))
        self.assertEqual(again, first.uid)
        self.assertEqual(len(manager.frames("neighbourhood")), 1)
        # Further comparisons join the open one's group as tabs instead of
        # splitting the window again, whatever their source.
        manager.host.focus_view(first.uid)
        second = self._open(win, key=("gene", "b"))
        third = self._open(win, key=("gene", "c"), beside="annotations")
        group = manager.host.group_uids("annotations")
        self.assertIn(second.uid, group)
        self.assertIn(third.uid, group)
        for panel in (first, second, third):
            manager.close_panel(panel.uid)
        _spin(self.app, 0.1)
        # Beside a source with no group to its right: a new group on its right.
        fourth = self._open(win, key=("gene", "d"), beside="annotations")
        right = manager.host.neighbor_group("annotations", "right")
        self.assertIsNotNone(right)
        self.assertIn(fourth.uid, manager.host.group_uids(right))
        # The new group gets most of the split (it was a large window).
        _spin(self.app, 0.3)
        source_area = manager.host._docks["annotations"].dockAreaWidget()
        new_area = manager.host._docks[fourth.uid].dockAreaWidget()
        self.assertGreater(new_area.width(), source_area.width())
        self.assertEqual(manager.host.floating_count(), 0)

    def test_inline_busy_state_while_the_rows_load(self) -> None:
        win = self._window()
        release = threading.Event()
        rows = _rows(win)

        def slow() -> List[dict]:
            release.wait(5)
            return rows

        panel = self._open(win, rows=None, loader=slow, key=("slow",))
        self.assertTrue(panel.is_busy())
        self.assertIsNone(panel.view)
        self.assertEqual(panel.message.text(), "Loading gene neighborhoods...")
        self.assertTrue(panel.busy_bar.isVisibleTo(panel))
        self.assertEqual(win._panel_manager.title(panel.uid), "lysR neighborhoods")
        release.set()
        self.assertTrue(self._wait(lambda: panel.view is not None))
        self.assertFalse(panel.is_busy())
        self.assertIsInstance(panel.view, NeighbourhoodView)
        self.assertEqual(panel.document.rows, rows)
        # Nothing found, or a failure: said inside the panel.
        empty = self._open(win, rows=None, loader=lambda: [], key=("empty",))
        self.assertTrue(self._wait(lambda: not empty.is_busy()))
        self.assertEqual(empty.message.text(), "No neighborhoods were found for this gene.")
        self.assertFalse(empty.busy_bar.isVisibleTo(empty))

        def broken() -> List[dict]:
            raise RuntimeError("no annotations")

        failed = self._open(win, rows=None, loader=broken, key=("broken",))
        self.assertTrue(self._wait(lambda: not failed.is_busy()))
        self.assertIn("RuntimeError: no annotations", failed.message.text())

    def test_a_panel_closed_while_loading_ignores_its_rows(self) -> None:
        win = self._window()
        release = threading.Event()
        rows = _rows(win)
        panel = self._open(win, rows=None, loader=lambda: (release.wait(5), rows)[1], key=("closing",))
        win._panel_manager.close_panel(panel.uid)
        release.set()
        _spin(self.app, 0.3)
        self.assertIsNone(panel.view)

    def test_not_saved_with_the_layout_kept_by_reset_and_reopened_from_memory(self) -> None:
        win = self._window()
        manager = win._panel_manager
        panel = self._open(win, key=("keep",), flank_sizes=(5, 10, 15, 30), initial_flank=15)
        uid = panel.uid
        panel.view.flank_box.setCurrentIndex(panel.view.flank_box.findData(10))
        self.assertNotIn("neighbourhood", [entry["kind"] for entry in manager.panel_list()])
        # Reset Layout keeps the comparison open (in a column on the right).
        manager.reset_layout()
        _spin(self.app, 0.2)
        self.assertTrue(manager.is_open(uid))
        self.assertEqual(manager.host.floating_count(), 0)
        # Close: disposed; Reopen Closed Panel rebuilds it with its options.
        manager.close_panel(uid)
        _spin(self.app, 0.1)
        self.assertIsNone(manager.frame(uid))
        self.assertEqual(manager.reopen_closed(), uid)
        _spin(self.app, 0.3)
        reopened = manager.panel(uid)
        self.assertIsNot(reopened, panel)
        self.assertEqual(reopened.view._context_flank, 10)
        self.assertEqual(reopened.title(), "lysR neighborhoods")
        # A rebuilt window does not bring the comparison back (3.11), and its
        # saved layout still restores.
        win.close()
        _spin(self.app, 0.1)
        again = self._window()
        self.assertEqual(again._panel_manager.frames("neighbourhood"), [])
        self.assertTrue(again._panel_manager.is_open("sample_tree"))
        self.assertEqual(again._panel_manager.host.floating_count(), 0)

    def test_the_palette_lists_only_open_comparisons(self) -> None:
        win = self._window()
        manager = win._panel_manager
        seen: List[List[str]] = []

        def fake_exec(palette: Any) -> int:
            seen.append([str(item.data(0, QtCore.Qt.ItemDataRole.UserRole + 1)) for item in palette.items()])
            return 0

        with mock.patch("annotation_explorer_polars.panels.palette.PanelPalette.exec", fake_exec):
            manager.show_palette()
            panel = self._open(win)
            manager.show_palette()
        self.assertNotIn("neighbourhood", seen[0])
        self.assertNotIn("alignment", seen[0])
        self.assertIn(panel.uid, seen[1])
        # Choosing it shows the panel.
        manager.host.set_current_tab("annotations")
        with mock.patch(
            "annotation_explorer_polars.panels.palette.PanelPalette.exec",
            lambda palette: palette._choose(palette.item_for(panel.uid), False),
        ):
            manager.show_palette()
        self.assertTrue(manager.is_visible(panel.uid))

    def test_no_duplicate_panel_entry_for_a_comparison(self) -> None:
        win = self._window()
        panel = self._open(win)
        menu = win._panel_manager.build_tab_menu(panel.uid)
        self.assertNotIn("Duplicate Panel", [a.text() for a in menu.actions()])
        self.assertIsNone(win._panel_manager.duplicate_panel(panel.uid))


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class AlignmentPanelTests(_Case):
    def test_rows_follow_the_sample_set_and_hover_goes_both_ways(self) -> None:
        win = self._window()
        uid = win._open_alignment_view(alignment_text=ALIGNMENT_TEXT, status_text="done", title="Alignment: g")
        _spin(self.app, 0.3)
        panel = win._panel_manager.panel(uid)
        self.assertIsInstance(panel, AlignmentPanel)
        self.assertEqual(win._panel_manager.title(uid), "Alignment: g")
        view = panel.view
        self.assertEqual(view.row_samples(), ["S1", "S2", "S5", None])
        self.assertEqual([view.is_row_shown(r) for r in range(4)], [True, True, True, True])
        win.ctx.set_subset(SampleSubset(("S1", "S5"), label="two"), origin="phylogeny")
        self._idle(win)
        # The consensus row stays.
        self.assertEqual([view.is_row_shown(r) for r in range(4)], [True, False, True, True])
        view._on_hover_row(2)
        win.ctx.flush_hover()
        self.assertEqual((win.ctx.hover_state.samples, win.ctx.hover_state.origin), (("S5",), uid))
        win.ctx.hover(["S1"], origin="sample_tree")
        win.ctx.flush_hover()
        # (Delivered on the bus's next turn; a loaded machine may take longer.)
        self.assertTrue(self._wait(lambda: view.linked_rows() == [0]), view.linked_rows())
        # No Close button inside a panel (its tab closes it).
        self.assertNotIn("Close", [b.text() for b in view.findChildren(QtWidgets.QPushButton)])


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class NeighbourhoodPanelReviewTests(_Case):
    """Phase 6 review fixes: retry, one collection per document, grouping,
    catching up, both openers, real key dispatch, scoped header."""

    def test_a_failed_or_empty_load_loads_again_when_asked_again(self) -> None:
        win = self._window()
        rows = _rows(win)

        def broken() -> List[dict]:
            raise ValueError("boom")

        failed = self._open(win, rows=None, loader=broken, key=("gene", "x"))
        self.assertTrue(self._wait(lambda: not failed.is_busy()))
        self.assertIn("ValueError: boom", failed.message.text())
        self.assertTrue(failed.retry_button.isVisibleTo(failed))
        # Compare neighborhoods again: the same panel, loading again.
        uid = win._open_neighbourhood_view(rows=None, loader=lambda: rows, title="Gene neighborhoods: lysR",
                                           source="gene_card", key=("gene", "x"))
        self.assertEqual(uid, failed.uid)
        self.assertTrue(self._wait(lambda: failed.view is not None))
        # The panel's own Try again button.
        calls: List[int] = []

        def flaky() -> List[dict]:
            calls.append(1)
            if len(calls) == 1:
                raise RuntimeError("first")
            return rows

        second = self._open(win, rows=None, loader=flaky, key=("flaky",))
        self.assertTrue(self._wait(lambda: not second.is_busy()))
        self.assertIn("RuntimeError: first", second.message.text())
        second.retry_button.click()
        self.assertTrue(self._wait(lambda: second.view is not None))
        self.assertEqual(len(calls), 2)
        # A panel that shows its rows is only focused: nothing loads again.
        win._open_neighbourhood_view(rows=None, loader=flaky, title="x", key=("flaky",))
        _spin(self.app, 0.2)
        self.assertEqual(len(calls), 2)
        # Nothing found: no Try again, but asking again loads again.
        empty = self._open(win, rows=None, loader=lambda: [], key=("empty",))
        self.assertTrue(self._wait(lambda: not empty.is_busy()))
        self.assertFalse(empty.retry_button.isVisibleTo(empty))
        win._open_neighbourhood_view(rows=None, loader=lambda: rows, title="x", key=("empty",))
        self.assertTrue(self._wait(lambda: empty.view is not None))

    def test_a_reopened_panel_joins_the_running_collection_and_follows_the_sample_set(self) -> None:
        win = self._window()
        manager = win._panel_manager
        rows = _rows(win)
        release = threading.Event()
        calls: List[int] = []

        def slow() -> List[dict]:
            calls.append(1)
            release.wait(5)
            return rows

        panel = self._open(win, rows=None, loader=slow, key=("closing",))
        manager.close_panel(panel.uid)
        _spin(self.app, 0.1)
        # The sample set changes while nothing shows the comparison.
        win.ctx.set_subset(SampleSubset(("S2",), label="S2"), origin="test")
        self._idle(win)
        uid = manager.reopen_closed()
        reopened = manager.panel(uid)
        self.assertTrue(reopened.is_busy())
        release.set()
        self.assertTrue(self._wait(lambda: reopened.view is not None))
        self.assertEqual(len(calls), 1)  # one collection, not two
        self.assertIsNone(panel.view)
        # Rows that arrive after the sample set changed follow it.
        self.assertEqual(_row_samples(reopened), ["S2"])
        # Closed, narrowed elsewhere, reopened: still following.
        manager.close_panel(uid)
        _spin(self.app, 0.1)
        win.ctx.set_subset(SampleSubset(("S3",), label="S3"), origin="test")
        self._idle(win)
        again = manager.panel(manager.reopen_closed())
        _spin(self.app, 0.3)
        self.assertEqual(_row_samples(again), ["S3"])

    def test_rows_that_arrive_late_follow_a_sample_set_set_while_loading(self) -> None:
        win = self._window()
        rows = _rows(win)
        release = threading.Event()
        panel = self._open(win, rows=None, loader=lambda: (release.wait(5), rows)[1], key=("late",))
        win.ctx.set_subset(SampleSubset(("S3",), label="S3"), origin="test")
        self._idle(win)
        release.set()
        self.assertTrue(self._wait(lambda: panel.view is not None))
        self.assertEqual(_row_samples(panel), ["S3"])
        self.assertTrue(panel.view.title_chips["genomes"].text().startswith("1 of "))

    def test_comparisons_collect_in_one_group(self) -> None:
        win = self._window()
        manager = win._panel_manager
        manager.host.focus_view("sample_tree")
        first = self._open(win, key=("gene", "a"))
        groups_before = len(manager.host.group_uids("sample_tree"))
        for key in ("b", "c"):
            # Each from the active view (the newest comparison, which has no
            # group on its right): a tab, not another split.
            other = self._open(win, key=("gene", key))
            self.assertEqual(manager.host.group_uids(other.uid), manager.host.group_uids(first.uid))
        self.assertEqual(len(manager.host.group_uids("sample_tree")), groups_before)
        self.assertIsNone(manager.host.neighbor_group(first.uid, "right"))
        # A second card of the same gene, in another sample: told apart.
        self.assertEqual(first.title(), "lysR neighborhoods")
        twin = self._open(win, key=("gene", "twin"), reference_sample="S2")
        self.assertEqual(twin.title(), "lysR neighborhoods (S2)")
        self.assertEqual(manager.title(twin.uid), "lysR neighborhoods (S2)")
        self.assertEqual(twin.view.title_text(), "Gene neighborhoods: lysR")

    def test_a_background_tab_catches_up_when_shown(self) -> None:
        win = self._window()
        manager = win._panel_manager
        panel = self._open(win)
        everyone = _row_samples(panel)
        canvas = panel.view.view
        manager.host.set_current_tab("annotations")
        _spin(self.app, 0.1)
        self.assertFalse(manager.is_visible(panel.uid))
        win.ctx.set_subset(SampleSubset(("S2",), label="S2"), origin="test")
        self._idle(win)
        win.ctx.hover(["S2"], origin="sample_tree")
        win.ctx.flush_hover()
        _spin(self.app, 0.1)
        self.assertEqual(_row_samples(panel), everyone)  # results only while visible
        win.ctx.clear_hover("sample_tree")
        win.ctx.flush_hover()
        manager.host.set_current_tab(panel.uid)
        _spin(self.app, 0.3)
        self.assertEqual(_row_samples(panel), ["S2"])
        # A hover is momentary: once shown, the panel takes the next one.
        win.ctx.hover(["S2"], origin="sample_tree")
        win.ctx.flush_hover()
        self.assertTrue(self._wait(lambda: canvas.linked_samples == frozenset({"S2"})))

    def _new_windows(self, win: Any) -> List[Any]:
        others = [w for w in win.app_controller.windows.explore_windows() if w is not win]
        self.windows += [w for w in others if w not in self.windows]
        return others

    def test_open_sample_from_a_plasmid_context_narrows_the_window(self) -> None:
        from annotation_explorer_polars.session.requests import OpenPlasmidContext

        win = self._window()
        win.ctx.request(OpenPlasmidContext.from_payload({"rows": _rows(win), "name": "pX"}), origin="sample_tree")
        _spin(self.app, 0.4)
        # In a new neighbourhood window: the comparison alone in its docking area.
        others = self._new_windows(win)
        self.assertEqual(len(others), 1)
        new = others[0]
        self.assertEqual(win._panel_manager.frames("neighbourhood"), [])
        manager = new._panel_manager
        panels = manager.frames("neighbourhood")
        self.assertEqual(len(panels), 1)
        panel = panels[0].panel()
        self.assertEqual(manager.host.open_uids(), [panel.uid])
        self.assertFalse(manager.is_open("annotations"))
        self._idle(new)
        subsets: List[Any] = []
        new.ctx.subsetChanged.connect(lambda subset, origin: subsets.append((subset, origin)))
        before = win.ctx.sample_set()
        panel.view._open_samples(["S3"])
        self._idle(new)
        self.assertEqual(subsets[-1], (SampleSubset(("S3",), label="Sample S3", kind="sample"), panel.uid))
        self.assertEqual(_row_samples(panel), ["S3"])
        # A window without a table opens the Sample tree (left) and
        # Annotations (below the comparison), showing that window's samples.
        self.assertTrue(manager.is_open("sample_tree"))
        self.assertTrue(manager.is_open("annotations"))
        below = manager.host.neighbor_group(panel.uid, "bottom")
        self.assertIsNotNone(below)
        self.assertIn("annotations", manager.host.group_uids(below))
        self.assertIn("sample_tree", manager.host.group_uids(manager.host.neighbor_group(panel.uid, "left")))
        self.assertEqual(list(new.ctx.sample_set()), ["S3"])
        # The origin window does not change.
        self.assertEqual(win.ctx.sample_set(), before)
        self.assertIsNone(win.ctx.draft.subset)
        # Asked again with a table open: nothing more opens.
        uids = sorted(manager.host.open_uids())
        panel.view._open_samples(["S3", "S4"])
        self._idle(new)
        self.assertEqual(sorted(manager.host.open_uids()), uids)

    def test_compare_neighborhoods_opens_in_a_new_window(self) -> None:
        win = self._window()
        row = win.run.ann_df.loc[win.run.ann_df["AnnotationID"] == "S1:1:1"].iloc[0].to_dict()
        card = win._open_gene_card(row)
        card._context_button.click()
        _spin(self.app, 0.3)
        self.assertEqual(win._panel_manager.frames("neighbourhood"), [])
        others = self._new_windows(win)
        self.assertEqual(len(others), 1)
        manager = others[0]._panel_manager
        panels = manager.frames("neighbourhood")
        self.assertEqual(len(panels), 1)
        panel = panels[0].panel()
        self.assertTrue(self._wait(lambda: panel.view is not None, timeout=10))
        gene = str(row.get("Gene") or "gene")
        self.assertEqual(panel.title(), f"{gene} neighborhoods")

    def test_ctrl_shift_p_in_a_neighbourhood_window_opens_its_picker(self) -> None:
        # Window > New Panel... and its shortcut work in a neighbourhood
        # window too, also with the focus inside the comparison panel, and
        # the picker opens over that window, not the one it came from.
        from annotation_explorer_polars.panels import manager as manager_module

        win = self._window()
        row = win.run.ann_df.loc[win.run.ann_df["AnnotationID"] == "S1:1:1"].iloc[0].to_dict()
        card = win._open_gene_card(row)
        card._context_button.click()
        _spin(self.app, 0.3)
        new = self._new_windows(win)[0]
        panel = new._panel_manager.frames("neighbourhood")[0].panel()
        self.assertTrue(self._wait(lambda: panel.view is not None, timeout=10))
        new._window_menu.aboutToShow.emit()
        self.assertIn(new._new_panel_action, new._window_menu.actions())
        self.assertEqual(new._new_panel_action.shortcut().toString(), "Ctrl+Shift+P")
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
        with mock.patch.object(manager_module, "PanelPalette", FakePalette):
            new.activateWindow()
            new.raise_()
            QtTest.QTest.qWaitForWindowActive(new, 2000)
            target = panel.view.view
            target.setFocus()
            _spin(self.app, 0.05)
            QtTest.QTest.keyClick(
                target,
                QtCore.Qt.Key.Key_P,
                QtCore.Qt.KeyboardModifier.ControlModifier | QtCore.Qt.KeyboardModifier.ShiftModifier,
            )
            _spin(self.app, 0.05)
        self.assertEqual(parents, [new])
        # The footer says what a default Open does here: a split beside the
        # comparison, not a tab of the active group (final check of round 3).
        self.assertEqual(len(hints), 1)
        self.assertIn("left of the comparison", hints[0])
        self.assertNotIn("tab", hints[0])

    def test_the_in_window_command_opens_beside_the_active_group(self) -> None:
        from annotation_explorer_polars.shell.commands import PANEL_OPEN_NEIGHBOURHOOD

        win = self._window()
        manager = win._panel_manager
        manager.host.focus_view("sample_tree")
        uid = win.app_controller.commands.run(
            PANEL_OPEN_NEIGHBOURHOOD, win, rows=_rows(win), title="Gene neighborhoods: lysR", source="gene_card"
        )
        _spin(self.app, 0.3)
        self.assertIn(uid, manager.host.group_uids(manager.host.neighbor_group("sample_tree", "right")))
        self.assertEqual(self._new_windows(win), [])

    def test_shortcuts_follow_the_keyboard_focus(self) -> None:
        win = self._window()
        panel = self._open(win)
        manager = win._panel_manager
        win.activateWindow()
        _spin(self.app, 0.1)
        table = next(
            w for w in manager.host._docks["sample_tree"].widget().findChildren(QtWidgets.QAbstractItemView)
            if w.isVisible()
        )
        table.setFocus()
        _spin(self.app, 0.05)
        panel.view.set_settings_open(False)
        QtTest.QTest.keyClick(table, QtCore.Qt.Key.Key_F, QtCore.Qt.KeyboardModifier.ControlModifier)
        QtTest.QTest.keyClick(table, QtCore.Qt.Key.Key_W, QtCore.Qt.KeyboardModifier.ControlModifier)
        _spin(self.app, 0.1)
        self.assertFalse(panel.view.settings_open())
        self.assertTrue(manager.is_open(panel.uid))
        graphics = panel.view.view.graphics_view
        graphics.setFocus()
        _spin(self.app, 0.05)
        QtTest.QTest.keyClick(graphics, QtCore.Qt.Key.Key_F, QtCore.Qt.KeyboardModifier.ControlModifier)
        _spin(self.app, 0.05)
        self.assertTrue(panel.view.settings_open())
        # No search box in a panel: Ctrl+F only opens the settings.
        self.assertIs(QtWidgets.QApplication.focusWidget(), panel.view.molecule_box)
        QtTest.QTest.keyClick(panel.view.molecule_box, QtCore.Qt.Key.Key_W, QtCore.Qt.KeyboardModifier.ControlModifier)
        _spin(self.app, 0.1)
        self.assertIsNone(manager.frame(panel.uid))
        self.assertTrue(manager.is_open("sample_tree"))

    def test_no_rows_in_the_sample_set_points_to_the_rail(self) -> None:
        win = self._window()
        panel = self._open(win)
        win.ctx.set_subset(SampleSubset(("no-such-sample",), label="none"), origin="test")
        self._idle(win)
        self.assertEqual(_row_samples(panel), [])
        self.assertEqual(panel.view.view.empty_text, NeighbourhoodView.OUT_OF_SCOPE_TEXT)
        self.assertTrue(panel.view.title_chips["genomes"].text().startswith("0 of "))
        win.ctx.set_subset(None)
        self._idle(win)
        self.assertEqual(panel.view.view.empty_text, "")
        self.assertFalse(panel.view.title_chips["genomes"].text().startswith("0 of "))


class NeighbourhoodViewPanelModeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = _app()

    def test_a_small_panel_folds_the_settings_column_until_it_grows(self) -> None:
        class Memory:
            settings_panel_open: Any = True

        rows = [
            {"sample": f"S{i}", "contig": "1",
             "genes": [{"start": 0, "end": 900, "strand": "+", "label": "g", "is_center": True,
                        "row_data": {"Gene": "g"}}]}
            for i in range(4)
        ]
        view = NeighbourhoodView(None, rows=rows, title="Gene neighborhoods: g", panel_mode=True,
                                 settings_memory=Memory)
        self.addCleanup(view.deleteLater)
        view.resize(1600, 900)
        view.show()
        _spin(self.app, 0.1)
        self.assertTrue(view.settings_open())
        view.resize(500, 360)
        _spin(self.app, 0.1)
        self.assertFalse(view.settings_open())
        self.assertTrue(view.settings_wanted())  # what the panel saves
        self.assertTrue(Memory.settings_panel_open)  # not the session's choice
        view.resize(1600, 900)
        _spin(self.app, 0.1)
        self.assertTrue(view.settings_open())
        # Folded by the user: it stays folded when the panel grows again.
        view.resize(500, 360)
        _spin(self.app, 0.1)
        view.settings_expand_btn.click()
        self.assertTrue(view.settings_open())
        view.settings_collapse_btn.click()
        view.resize(1600, 900)
        _spin(self.app, 0.1)
        self.assertFalse(view.settings_open())
        self.assertFalse(view.settings_wanted())
        # The footer hint ends in '...' instead of being cut; text stays whole.
        self.assertEqual(view.footer_hint.minimumSizeHint().width(), 0)
        self.assertEqual(view.footer_hint.toolTip(), view.footer_hint.text())


class AlignmentReviewTests(_Case):
    def test_rows_match_samples_whose_names_have_special_characters(self) -> None:
        from annotation_explorer_polars.annotation_detail_dialog import AlignmentView

        hovered: List[Any] = []
        view = AlignmentView(None, alignment_text=">a|ST492_isolate_1|1|g\nMKA\n>b|S2|1|g\nMKA\n",
                             status_text="", panel_mode=True, sample_hover_callback=hovered.append)
        self.addCleanup(view.deleteLater)
        view.set_sample_scope(["ST492 isolate#1", "S9"])
        self.assertEqual([view.is_row_shown(r) for r in range(3)], [True, False, True])
        view.set_linked_samples(["ST492 isolate#1"])
        self.assertEqual(view.linked_rows(), [0])
        view._on_hover_row(0)
        self.assertEqual(hovered, ["ST492 isolate#1"])  # the window's spelling

    @unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
    def test_a_hidden_or_reopened_alignment_follows_the_sample_set(self) -> None:
        win = self._window()
        manager = win._panel_manager
        uid = win._open_alignment_view(alignment_text=ALIGNMENT_TEXT, status_text="done", title="Alignment: g")
        _spin(self.app, 0.3)
        view = manager.panel(uid).view
        manager.host.set_current_tab("annotations")
        win.ctx.set_subset(SampleSubset(("S1", "S5"), label="two"), origin="phylogeny")
        self._idle(win)
        manager.host.set_current_tab(uid)
        _spin(self.app, 0.3)
        self.assertEqual([view.is_row_shown(r) for r in range(4)], [True, False, True, True])
        manager.close_panel(uid)
        _spin(self.app, 0.1)
        win.ctx.set_subset(SampleSubset(("S2",), label="one"), origin="phylogeny")
        self._idle(win)
        reopened = manager.panel(manager.reopen_closed())
        _spin(self.app, 0.3)
        self.assertEqual([reopened.view.is_row_shown(r) for r in range(4)], [False, True, False, True])
        # Each alignment run is its own panel, in the same group.
        second = win._open_alignment_view(alignment_text=ALIGNMENT_TEXT, status_text="done", title="Alignment: g")
        _spin(self.app, 0.2)
        self.assertNotEqual(second, reopened.uid)
        self.assertEqual(manager.host.group_uids(second), manager.host.group_uids(reopened.uid))


class StandaloneViewTests(unittest.TestCase):
    """The body without panel mode: a standalone viewer (Phase 10 retired the
    ``IgvCollectionDialog`` window shim that wrapped it)."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.app = _app()

    def test_the_window_shims_are_gone(self) -> None:
        from annotation_explorer_polars import annotation_detail_dialog

        self.assertFalse(hasattr(igv_collection_view, "IgvCollectionDialog"))
        self.assertFalse(hasattr(annotation_detail_dialog, "AlignmentResultDialog"))

    def test_a_standalone_viewer_keeps_its_own_filters_and_window_shortcuts(self) -> None:
        rows = [
            {"sample": f"S{i}", "contig": "1", "sample_meta": {"Country": "Romania" if i % 2 else "Italy"},
             "genes": [{"start": 0, "end": 900, "strand": "+", "label": "g", "is_center": True,
                        "row_data": {"Gene": "g"}}]}
            for i in range(4)
        ]
        NeighbourhoodView.settings_panel_open = None
        view = NeighbourhoodView(None, rows=rows, title="Gene neighborhoods: g")
        self.addCleanup(view.deleteLater)
        self.assertFalse(view.is_panel_mode())
        self.assertEqual(view.title_text(), "Gene neighborhoods: g")
        self.assertTrue(view.add_filter_btn.isEnabled())
        self.assertEqual(view.sample_query.placeholderText(), "Search samples and metadata...")
        view._metadata_filters = {"Country": ["Italy"]}
        view._apply_filters()
        self.assertEqual(sorted(r["sample"] for r in view._visible_rows), ["S0", "S2"])
        shortcuts = view.findChildren(QtGui.QShortcut)
        self.assertTrue(shortcuts)
        self.assertTrue(all(s.context() == QtCore.Qt.ShortcutContext.WindowShortcut for s in shortcuts))

    def test_the_second_escape_passes_to_the_viewer_s_container(self) -> None:
        rows = [
            {"sample": f"S{i}", "contig": "1",
             "genes": [{"start": 0, "end": 900, "strand": "+", "label": "g", "is_center": True,
                        "row_data": {"Gene": "g"}}]}
            for i in range(3)
        ]
        host = QtWidgets.QDialog()
        self.addCleanup(host.deleteLater)
        layout = QtWidgets.QVBoxLayout(host)
        view = NeighbourhoodView(host, rows=rows, title="Gene neighborhoods: g")
        layout.addWidget(view)
        host.resize(1200, 700)
        host.show()
        _spin(self.app, 0.1)
        canvas = view.view
        canvas._select_gene(canvas.layout_result.rows[0].genes[0])
        QtTest.QTest.keyClick(canvas.graphics_view, QtCore.Qt.Key.Key_Escape)
        self.assertFalse(canvas.has_selection())
        self.assertTrue(host.isVisible())  # the first Esc clears the pin
        QtTest.QTest.keyClick(canvas.graphics_view, QtCore.Qt.Key.Key_Escape)
        _spin(self.app, 0.05)
        self.assertFalse(host.isVisible())  # the second is the container's (a dialog closes)


if __name__ == "__main__":
    unittest.main()
