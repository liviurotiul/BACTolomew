"""Multiple independent explore windows (window plan Phase 7, 8.8, 3.10, 3.11).

- New Explore Window: a clean rail and the default layout.
- Duplicate Window: a copy of the window's filters and layout, once; the two
  windows are independent afterwards.
- "Open clade / cluster / plasmid in new explore window": a new window whose
  rail carries the subset as a locked window subset.
- Window numbers and colors: ``[n]`` in the title and a color stripe while
  several windows are open; cards show their owner window's number and color.
- One rail and one context per window; nothing crosses windows (D7): no
  filter, hover, selection or result is ever delivered to another window's
  panels, and closing one window leaves the other untouched.
- No move-panel command and no cross-window drag (D2).
- Multi-window persistence: each window saves into the slot of its number;
  the windows open when the run is left are reopened with it (layout only, D5).
- Neighbourhood windows: Compare neighborhoods and plasmid contexts open in a
  new explore window that copies the requesting window's filters, sort and
  rail, holds only the comparison at first, and is independent afterwards.
"""
from __future__ import annotations

import os
import unittest
from types import SimpleNamespace
from typing import Any, List
from unittest import mock

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

from PyQt6 import QtCore, QtGui, QtWidgets

try:
    from tests.test_dock_host import _spin
    from tests.test_window_manager import _alive, _document_panels, _RealWindowCase, _wait
except ImportError:
    from test_dock_host import _spin
    from test_window_manager import _alive, _document_panels, _RealWindowCase, _wait

from annotation_explorer_polars.docking import docking_available
from annotation_explorer_polars.docking.dock_host import Placement
from annotation_explorer_polars.model.filter_state import FilterState, SampleSubset
from annotation_explorer_polars.panels.pan import request_pan_subset_window
from annotation_explorer_polars.session import bus as bus_events
from annotation_explorer_polars.session.requests import OpenSubsetInNewWindow
from annotation_explorer_polars.shell import commands as cmd
from annotation_explorer_polars.shell.app_controller import CASCADE_OFFSET, AppController
from annotation_explorer_polars.shell.layout_store import START_DEFAULT_LAYOUT, START_EMPTY_LAYOUT
from annotation_explorer_polars.shell.window_manager import WindowManager
from annotation_explorer_polars.theme import window_color
from annotation_explorer_polars.tree_view import RaxmlTreeView


def _menu_action(menu: QtWidgets.QMenu, name: str) -> QtGui.QAction:
    menu.aboutToShow.emit()
    return next(a for a in menu.actions() if a.objectName() == name)


class _BusRecorder:
    """Every event one window's bus delivers to a subscriber of each type."""

    def __init__(self, ctx: Any, owner: str = "recorder") -> None:
        self.events: List[tuple] = []
        for event in (
            bus_events.RESULT_APPLIED,
            bus_events.HOVER_CHANGED,
            bus_events.SELECTION_CHANGED,
            bus_events.BUSY_CHANGED,
        ):
            ctx.bus.subscribe(owner, event, lambda payload, _e=event: self.events.append((_e, payload)))


def _large_screen() -> Any:
    """The offscreen platform's screen is 800 x 800, smaller than a window:
    pretend it is large so the cascade fits."""
    return mock.patch.object(QtGui.QScreen, "availableGeometry", return_value=QtCore.QRect(0, 0, 4000, 3000))


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class MultiWindowTests(_RealWindowCase):
    def _controller(self, **kwargs: Any) -> AppController:
        return AppController(self._run_data(), **kwargs)

    def _open(self, controller: AppController, **kwargs: Any) -> Any:
        win = controller.open_explore(**kwargs)
        self.windows.append(win)
        self._settle(win)
        return win

    def _first(self, controller: AppController) -> Any:
        win = controller.open_first_window()
        self.windows.append(win)
        self._settle(win)
        return win

    @staticmethod
    def _clean() -> FilterState:
        """What a clean rail applies ("Hide empty contigs" starts checked)."""
        return FilterState(hide_empty_contigs=True)

    def _apply(self, win: Any) -> None:
        self.assertTrue(_wait(self.app, lambda: not win.ctx.is_busy(), 20))
        _spin(self.app, 0.1)

    # ------------------------------------------------------------ opening

    def test_new_explore_window_has_a_clean_rail_and_the_default_layout(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        first.move(100, 100)
        first.gene_query.setText("blaOXA")
        first._set_sample_subset_filter(label="Three", samples=["S1", "S2", "S3"])
        manager = first._panel_manager
        manager.close_panel("geo")
        self._apply(first)
        action = _menu_action(first._window_menu, "NewExploreWindowAction")
        self.assertTrue(action.shortcut().isEmpty())
        before = len(controller.windows.explore_windows())
        with _large_screen():
            action.trigger()
        _spin(self.app, 0.1)
        windows = controller.windows.explore_windows()
        self.assertEqual(len(windows), before + 1)
        second = windows[-1]
        self.windows.append(second)
        self._settle(second)
        self.assertEqual(controller.windows.number(second), 2)
        self.assertEqual(second.ctx.applied, self._clean())
        # The Sample tree and Annotations panels show the window filter defaults.
        self.assertTrue(second._panel_manager.panel("sample_tree").hide_empty_contigs_chk.isChecked())
        self.assertFalse(second.annotations_panel.hide_hypo_chk.isChecked())
        self.assertEqual(second.gene_query.text(), "")
        self.assertIsNone(second.ctx.draft.subset)
        self.assertEqual(len(second.ctx.sample_set()), len(first.run.meta_df))
        self.assertTrue(second._panel_manager.is_open("geo"))  # Classic, not window 1's layout
        self.assertFalse(manager.is_open("geo"))
        # Window 1 keeps its filters.
        self.assertEqual(first.ctx.applied.gene, "blaOXA")
        self.assertIsNotNone(first.ctx.applied.subset)
        # Opened from window 1: its size, offset down and right.
        self.assertEqual(second.size(), first.size())
        self.assertEqual(second.pos(), first.pos() + QtCore.QPoint(CASCADE_OFFSET, CASCADE_OFFSET))

    def test_titles_numbers_and_colors_while_several_windows_are_open(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        run_name = first.run.output_dir.name
        self.assertEqual(first.windowTitle(), f"BACTolomew - {run_name}")
        self.assertFalse(first._window_stripe.isVisibleTo(first))
        second = self._open(controller, source=first)
        self.assertEqual(first.windowTitle(), f"BACTolomew - {run_name} [1]")
        self.assertEqual(second.windowTitle(), f"BACTolomew - {run_name} [2]")
        self.assertTrue(first._window_stripe.isVisibleTo(first))
        self.assertTrue(second._window_stripe.isVisibleTo(second))
        self.assertIn(window_color(1), first._window_stripe.styleSheet())
        self.assertIn(window_color(2), second._window_stripe.styleSheet())
        self.assertNotEqual(window_color(1), window_color(2))
        self.assertEqual((first.window_number, second.window_number), (1, 2))
        # Cards: owner number in the title, owner color as a stripe.
        card = second._open_sample_card_for_sample("S1")
        self.assertEqual(card.windowTitle(), "Sample card: S1 [2]")
        self.assertEqual(card.owner_color(), window_color(2))
        self.assertTrue(card._owner_stripe.isVisibleTo(card))
        # The Window menu: each window once (its title has the number), with
        # its color swatch, and the card.
        menu = first._window_menu
        menu.aboutToShow.emit()
        texts = [a.text() for a in menu.actions()]
        self.assertIn(first.windowTitle(), texts)
        self.assertIn(second.windowTitle(), texts)
        self.assertNotIn(f"{first.windowTitle()} [1]", texts)
        entry = next(a for a in menu.actions() if a.text() == second.windowTitle())
        self.assertFalse(entry.icon().isNull())
        # The only shortcut is New Panel... (Ctrl+Shift+P).
        self.assertEqual(
            {a.text(): a.shortcut().toString() for a in menu.actions() if not a.shortcut().isEmpty()},
            {"New Panel...": "Ctrl+Shift+P"},
        )
        # Window 2 closes: number, stripe and badge go away in window 1.
        second.close()
        _spin(self.app, 0.1)
        self.assertEqual(first.windowTitle(), f"BACTolomew - {run_name}")
        self.assertFalse(first._window_stripe.isVisibleTo(first))
        self.assertIs(card.owner(), first)
        self.assertEqual(card.windowTitle(), "Sample card: S1")
        self.assertEqual(card.owner_color(), "")
        self.assertFalse(card._owner_stripe.isVisibleTo(card))
        # A new window takes the lowest free number again.
        third = self._open(controller)
        self.assertEqual(third.window_number, 2)
        self.assertEqual(third.windowTitle(), f"BACTolomew - {run_name} [2]")

    def test_duplicate_window_copies_filters_and_layout_once(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        predictor = first.run.predictor_options[0]
        first.sample_query.setText("s")
        first.gene_query.setText("blaOXA")
        first.product_query.setText("lactamase")
        first.molecule_filter_box.setCurrentIndex(first.molecule_filter_box.findData("chromosome"))
        first.annotations_panel.hide_hypo_chk.setChecked(True)
        first._panel_manager.panel("sample_tree").hide_empty_contigs_chk.setChecked(False)
        model = first._predictor_model
        for row in range(model.rowCount()):
            if model.item(row).data(QtCore.Qt.ItemDataRole.UserRole) == predictor:
                model.item(row).setCheckState(QtCore.Qt.CheckState.Checked)
        first.metadata_filters = {"Country": ["Italy", "Romania"]}
        first._render_metadata_filters()
        first._apply_rail_filters()
        first._set_sample_subset_filter(label="Five", samples=["S1", "S2", "S3", "S4", "S5"])
        self._apply(first)
        manager = first._panel_manager
        extra = manager.open_new_instance("annotations")
        manager.panel(extra).set_view_mode("Context")
        manager.close_panel("geo")
        manager.rename_panel("phylogeny", "Pan views")
        _spin(self.app, 0.2)
        source_filters = first.ctx.draft
        source_layout = sorted(manager.uids("annotations"))
        computes_before = first.ctx.scheduler.compute_count

        action = _menu_action(first._window_menu, "DuplicateWindowAction")
        action.trigger()
        _spin(self.app, 0.1)
        dup = controller.windows.explore_windows()[-1]
        self.windows.append(dup)
        self._settle(dup)
        self._apply(dup)
        # Filters: the rail shows them and they are applied as they were.
        self.assertEqual(dup.ctx.applied, source_filters)
        self.assertEqual(FilterState.from_widgets(dup), source_filters)
        self.assertEqual(dup.sample_query.text(), "s")
        self.assertEqual(dup.gene_query.text(), "blaOXA")
        self.assertEqual(dup.product_query.text(), "lactamase")
        self.assertEqual(dup.molecule_filter_box.currentData(), "chromosome")
        self.assertTrue(dup.annotations_panel.hide_hypo_chk.isChecked())
        self.assertTrue(dup.ctx.applied.hide_hypothetical)
        dup_tree = dup._panel_manager.panel("sample_tree")
        self.assertFalse(dup_tree.hide_empty_contigs_chk.isChecked())
        self.assertFalse(dup.ctx.applied.hide_empty_contigs)
        self.assertEqual(dup._selected_predictors(), [predictor])
        self.assertEqual(dup.metadata_filters, {"Country": ["Italy", "Romania"]})
        self.assertEqual(dup.ctx.draft.subset.label, "Five")
        self.assertEqual(dup.ctx.sample_set(), first.ctx.sample_set())
        self.assertEqual(dup.ctx.scheduler.compute_count, 1)  # one compute for the copy
        # Layout: the same panels, options and titles.
        dup_manager = dup._panel_manager
        self.assertEqual(sorted(dup_manager.uids("annotations")), source_layout)
        self.assertEqual(dup_manager.panel(extra).view_mode(), "Context")
        self.assertFalse(dup_manager.is_open("geo"))
        self.assertEqual(dup_manager.title("phylogeny"), "Pan views")
        self.assertEqual(
            [sorted(g) for g in dup_manager.host.groups()],
            [sorted(g) for g in manager.host.groups()],
        )
        # Independent from now on.
        dup.gene_query.setText("")
        dup._clear_sample_subset_filter()
        self._apply(dup)
        self._apply(first)
        self.assertEqual(first.ctx.draft, source_filters)
        self.assertEqual(first.gene_query.text(), "blaOXA")
        self.assertEqual(first.ctx.scheduler.compute_count, computes_before)
        dup_manager.close_panel(extra)
        self.assertTrue(manager.is_open(extra))

    def test_window_filters_stay_in_effect_while_their_panels_are_closed(self) -> None:
        # "Hide hypothetical" (Annotations) and "Hide contigs without hits"
        # (Sample tree) are window filters shown in panels: closing a panel
        # changes nothing, and the reopened panel shows what is in effect.
        controller = self._controller()
        win = self._first(controller)
        manager = win._panel_manager
        win.annotations_panel.hide_hypo_chk.setChecked(True)
        manager.panel("sample_tree").hide_empty_contigs_chk.setChecked(False)
        self._apply(win)
        manager.close_panel("annotations")
        manager.close_panel("sample_tree")
        _spin(self.app, 0.1)
        self.assertFalse(manager.is_open("annotations"))
        self.assertFalse(manager.is_open("sample_tree"))
        computes = win.ctx.scheduler.compute_count
        win.gene_query.setText("blaOXA")
        self._apply(win)
        self.assertEqual(win.ctx.scheduler.compute_count, computes + 1)
        self.assertTrue(win.ctx.applied.hide_hypothetical)
        self.assertFalse(win.ctx.applied.hide_empty_contigs)
        manager.open_panel("annotations")
        manager.open_panel("sample_tree")
        _spin(self.app, 0.1)
        self.assertTrue(win.annotations_panel.hide_hypo_chk.isChecked())
        self.assertFalse(manager.panel("sample_tree").hide_empty_contigs_chk.isChecked())
        self.assertEqual(win.ctx.scheduler.compute_count, computes + 1)

    def test_rail_note_shows_hide_hypothetical_while_annotations_is_closed(self) -> None:
        # Review fix: with no Annotations panel open, "Hide hypothetical" on
        # is named in the rail, with a link that opens the panel again.
        # "Hide contigs without hits" sits on the Sample tree and needs no note.
        controller = self._controller()
        win = self._first(controller)
        note = win.window_options_note
        manager = win._panel_manager
        self.assertFalse(note.isVisibleTo(win))  # defaults: nothing to say
        win.annotations_panel.hide_hypo_chk.setChecked(True)
        self._apply(win)
        self.assertFalse(note.isVisibleTo(win))  # the panel shows it
        # A preset without an Annotations panel (restored without per-panel
        # signals) shows the note too.
        if win._panel_manager.apply_preset("pangenome") and not manager.is_open("annotations"):
            _spin(self.app, 0.1)
            self.assertTrue(note.isVisibleTo(win))
            manager.apply_preset("classic")
            _spin(self.app, 0.1)
            self.assertFalse(note.isVisibleTo(win))
        manager.close_panel("annotations")
        _spin(self.app, 0.1)
        self.assertTrue(note.isVisibleTo(win))
        self.assertIn("Hide hypothetical is on", note.text())
        self.assertNotIn("contigs", note.text())
        win.ctx.update_draft(win.ctx.draft.replace(hide_empty_contigs=False))
        self._apply(win)
        self.assertNotIn("contigs", note.text())
        note.linkActivated.emit("annotations")
        _spin(self.app, 0.1)
        self.assertTrue(manager.is_open("annotations"))
        self.assertFalse(note.isVisibleTo(win))
        self.assertTrue(win.annotations_panel.hide_hypo_chk.isChecked())
        # Back to the default while the panel is closed (for example from a
        # copied draft): the note goes away.
        manager.close_panel("annotations")
        _spin(self.app, 0.1)
        self.assertTrue(note.isVisibleTo(win))
        win.ctx.update_draft(win.ctx.draft.replace(hide_hypothetical=False))
        self._apply(win)
        self.assertFalse(note.isVisibleTo(win))

    def test_hide_empty_contigs_lives_on_the_sample_tree(self) -> None:
        # Feedback 2026-10-07: "Hide contigs without hits" moved from the
        # Annotations strip to the Sample tree panel, the only panel it
        # changes. Annotations instances never have it.
        controller = self._controller()
        win = self._first(controller)
        manager = win._panel_manager
        panel = win.annotations_panel
        tree = manager.panel("sample_tree")
        self.assertIsNone(panel.findChild(QtWidgets.QCheckBox, "HideEmptyContigsCheck"))
        self.assertTrue(panel.hide_hypo_chk.isVisibleTo(panel))
        check = tree.findChild(QtWidgets.QCheckBox, "HideEmptyContigsCheck")
        self.assertIs(check, tree.hide_empty_contigs_chk)
        self.assertTrue(check.isVisibleTo(tree))
        self.assertTrue(check.isChecked())
        computes = win.ctx.scheduler.compute_count
        check.click()
        self._apply(win)
        self.assertFalse(win.ctx.draft.hide_empty_contigs)
        self.assertFalse(win.ctx.applied.hide_empty_contigs)
        self.assertEqual(win.ctx.scheduler.compute_count, computes + 1)
        # A duplicated Annotations panel has no such checkbox either.
        extra = manager.duplicate_panel("annotations")
        _spin(self.app, 0.1)
        self.assertIsNone(manager.panel(extra).findChild(QtWidgets.QCheckBox, "HideEmptyContigsCheck"))
        # A layout without a Sample tree, then back: the new tree shows the
        # value in effect.
        manager.apply_preset("geo_table")
        _spin(self.app, 0.1)
        self.assertFalse(win.ctx.applied.hide_empty_contigs)
        manager.apply_preset("classic")
        _spin(self.app, 0.1)
        self.assertTrue(manager.is_open("sample_tree"))
        tree = manager.panel("sample_tree")
        self.assertFalse(tree.hide_empty_contigs_chk.isChecked())
        self.assertFalse(win.ctx.applied.hide_empty_contigs)

    def test_duplicate_does_not_copy_document_panels(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        first._open_neighbourhood_view(rows=[], title="Gene neighborhoods: x", key=("gene", "x"))
        _spin(self.app, 0.2)
        self.assertEqual(len(_document_panels(first, "neighbourhood")), 1)
        dup = controller.duplicate_explore(first)
        self.windows.append(dup)
        self._settle(dup)
        self.assertEqual(_document_panels(dup, "neighbourhood"), [])
        self.assertEqual(len(_document_panels(first, "neighbourhood")), 1)

    def test_open_a_clade_in_a_new_window_locks_its_subset(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        all_samples = tuple(first.ctx.sample_set())
        panel = SimpleNamespace(ctx=first.ctx, uid="phylogeny")
        request_pan_subset_window(panel, {"name": "7", "samples": ["S2", "S4", "S9"]}, "clade")
        _spin(self.app, 0.1)
        windows = controller.windows.explore_windows()
        self.assertEqual(len(windows), 2)
        sub = windows[-1]
        self.windows.append(sub)
        self._settle(sub)
        self._apply(sub)
        subset = sub.window_subset
        self.assertEqual((subset.label, subset.kind, subset.samples), ("Clade 7", "clade", ("S2", "S4", "S9")))
        self.assertEqual(sub.ctx.applied.window_subset, subset)
        self.assertEqual(set(sub.ctx.sample_set()), {"S2", "S4"})
        self.assertIsNone(sub.ctx.applied.subset)
        self.assertTrue(sub.windowTitle().endswith("[2] - Clade 7"))
        # The rail shows it as a locked chip: no Remove button.
        chip = sub.findChild(QtWidgets.QWidget, "WindowSubsetChip")
        self.assertIsNotNone(chip)
        texts = [w.text() for w in chip.findChildren(QtWidgets.QLabel)]
        self.assertIn("Clade 7 | n=3", texts)
        self.assertIn("Locked", texts)
        self.assertEqual(chip.findChildren(QtWidgets.QPushButton), [])
        # Default layout and clean rail otherwise; window 1 is untouched.
        self.assertTrue(sub._panel_manager.is_open("geo"))
        self.assertEqual(first.ctx.applied, self._clean())
        self.assertEqual(tuple(first.ctx.sample_set()), all_samples)
        # The window subset survives every rail edit: another subset narrows
        # inside it, clearing that subset brings back the clade, not all samples.
        sub._set_sample_subset_filter(label="Two", samples=["S2", "S3"])
        self._apply(sub)
        self.assertEqual(set(sub.ctx.sample_set()), {"S2"})
        sub._clear_sample_subset_filter()
        sub.gene_query.setText("")
        sub.metadata_filters = {"Country": ["Italy"]}
        sub._apply_rail_filters()
        self._apply(sub)
        self.assertEqual(sub.ctx.applied.window_subset, subset)
        sub.metadata_filters = {}
        sub._render_metadata_filters()
        sub._apply_rail_filters()
        self._apply(sub)
        self.assertEqual(set(sub.ctx.sample_set()), {"S2", "S4"})
        self.assertIsNotNone(sub.findChild(QtWidgets.QWidget, "WindowSubsetChip"))
        # Duplicating the clade window keeps its lock.
        dup = controller.duplicate_explore(sub)
        self.windows.append(dup)
        self._settle(dup)
        self._apply(dup)
        self.assertEqual(dup.window_subset, subset)
        self.assertEqual(set(dup.ctx.sample_set()), {"S2", "S4"})

    def test_the_request_and_the_command_open_subset_windows(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        cluster = SampleSubset(("S1", "S3"), label="Cluster A", kind="cluster")
        first.ctx.request(OpenSubsetInNewWindow(cluster), origin="phylogeny")
        _spin(self.app, 0.1)
        sub = controller.windows.explore_windows()[-1]
        self.windows.append(sub)
        self._settle(sub)
        self.assertEqual(sub.window_subset, cluster)
        plasmid = SampleSubset(("S5",), label="Plasmid pX", kind="plasmid")
        third = controller.commands.run(cmd.WINDOW_OPEN_SUBSET, first, subset=plasmid)
        self.windows.append(third)
        self._settle(third)
        self.assertEqual(tuple(third.ctx.sample_set()), ("S5",))
        self.assertEqual(controller.windows.number(third), 3)

    def test_the_pan_views_offer_open_in_new_explore_window(self) -> None:
        # Clade and cluster labels of the pan trees: a right-click menu.
        tree = RaxmlTreeView()
        self.addCleanup(tree.deleteLater)
        got: List[Any] = []
        tree.cladeWindowRequested.connect(lambda payload: got.append(("clade", payload)))
        tree.clusterWindowRequested.connect(lambda payload: got.append(("cluster", payload)))
        tree.cladeClicked.connect(lambda payload: got.append(("filter", payload)))
        payload = {"name": "7", "samples": ["S1", "S2"]}
        menu = tree.group_menu("clade", payload)
        self.assertEqual([a.text() for a in menu.actions()], ["Show only this clade", "Open clade in new explore window"])
        self.assertTrue(all(a.shortcut().isEmpty() for a in menu.actions()))
        menu.actions()[1].trigger()
        menu.actions()[0].trigger()
        menu = tree.group_menu("cluster", payload)
        self.assertEqual(menu.actions()[1].text(), "Open cluster in new explore window")
        menu.actions()[1].trigger()
        self.assertEqual(got, [("clade", payload), ("filter", payload), ("cluster", payload)])
        # The plasmid table's row menu.
        from annotation_explorer_polars.plasmidome_view import PlasmidomePanel

        plasmidome = PlasmidomePanel.__new__(PlasmidomePanel)
        QtWidgets.QWidget.__init__(plasmidome)
        plasmidome.table = QtWidgets.QTableView(plasmidome)
        self.addCleanup(plasmidome.deleteLater)
        requested: List[Any] = []
        plasmidome.sampleSubsetWindowRequested.connect(requested.append)
        menu = PlasmidomePanel.row_menu(plasmidome, {"name": "pX", "samples": ["S5"], "cluster_key": "c"})
        texts = [a.text() for a in menu.actions()]
        self.assertEqual(texts, ["Filter by this plasmid", "Open plasmid in new explore window", "Open context"])
        menu.actions()[1].trigger()
        self.assertEqual(requested, [{"name": "pX", "samples": ["S5"]}])

    def test_every_pan_panel_forwards_the_new_window_request(self) -> None:
        from annotation_explorer_polars.panels.neighbourhood import NeighbourhoodPanel
        from annotation_explorer_polars.panels.pan import PAN_VIEW_PANELS

        class _Tab(QtWidgets.QWidget):
            samplesHovered = QtCore.pyqtSignal(object)
            cladeFilterRequested = QtCore.pyqtSignal(object)
            clusterFilterRequested = QtCore.pyqtSignal(object)
            plasmidFilterRequested = QtCore.pyqtSignal(object)
            subsetWindowRequested = QtCore.pyqtSignal(object, str)
            plasmidContextRequested = QtCore.pyqtSignal(object)
            sampleOpenRequested = QtCore.pyqtSignal(object)

            def single_view_widget(self) -> None:
                return None

        ctx = SimpleNamespace(requests=[])
        ctx.request = lambda request, origin="": ctx.requests.append((request, origin))
        payload = {"name": "7", "samples": ["S1", "S2"]}
        # The neighbourhood panel has no tree copy to wire any more.
        self.assertFalse(hasattr(NeighbourhoodPanel, "_wire_pan_copy"))
        # The pan panels (Phase 9).
        for cls, wire in [(cls, lambda p, t: p.set_tab(t)) for cls in PAN_VIEW_PANELS]:
            panel = cls.__new__(cls)
            QtWidgets.QWidget.__init__(panel)
            QtWidgets.QVBoxLayout(panel)
            panel.ctx, panel.uid = ctx, f"{cls.__name__}-uid"
            tab = _Tab()
            wire(panel, tab)
            tab.subsetWindowRequested.emit(payload, "clade")
            request, origin = ctx.requests[-1]
            self.assertIsInstance(request, OpenSubsetInNewWindow, cls.__name__)
            self.assertEqual(request.subset, SampleSubset(("S1", "S2"), label="Clade 7", kind="clade"))
            self.assertEqual(origin, panel.uid)
            panel.deleteLater()
        self.assertEqual(len(ctx.requests), len(PAN_VIEW_PANELS))

    def _assert_same_dock_state(self, got: bytes, expected: bytes) -> None:
        import re

        def parse(state: bytes) -> Any:
            xml = bytes(QtCore.qUncompress(QtCore.QByteArray(state))).decode()
            sizes = [[int(v) for v in group.split()] for group in re.findall(r"<Sizes>([^<]*)</Sizes>", xml)]
            return re.sub(r"<Sizes>[^<]*</Sizes>", "<Sizes/>", xml), sizes

        got_xml, got_sizes = parse(got)
        expected_xml, expected_sizes = parse(expected)
        self.assertEqual(got_xml, expected_xml)
        self.assertEqual([len(s) for s in got_sizes], [len(s) for s in expected_sizes])
        self.assertEqual(got_sizes, expected_sizes)

    def test_duplicate_copies_the_rail_width_and_the_dock_split(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        manager = first._panel_manager
        first._apply_rail_state(360, False)
        _spin(self.app, 0.3)
        width, collapsed = first._rail_state()
        self.assertEqual((width, collapsed), (360, False))
        dup = controller.duplicate_explore(first)
        self.windows.append(dup)
        self._settle(dup)
        _spin(self.app, 0.3)
        self.assertEqual(dup._rail_state(), (360, False))
        # The same arrangement, splitter sizes to the pixel: the copy is
        # restored again once the new window is laid out at its size.
        self._assert_same_dock_state(dup._layout_copy().dock_state, first._layout_copy().dock_state)
        # A copy of the copy does not drift either.
        dup2 = controller.duplicate_explore(dup)
        self.windows.append(dup2)
        self._settle(dup2)
        _spin(self.app, 0.3)
        self._assert_same_dock_state(dup2._layout_copy().dock_state, first._layout_copy().dock_state)
        self.assertEqual(
            [sorted(g) for g in dup._panel_manager.host.groups()],
            [sorted(g) for g in manager.host.groups()],
        )
        # A collapsed rail is copied collapsed.
        first._apply_rail_state(300, True)
        _spin(self.app, 0.2)
        self.assertTrue(first._rail_state()[1])
        collapsed_dup = controller.duplicate_explore(first)
        self.windows.append(collapsed_dup)
        self._settle(collapsed_dup)
        _spin(self.app, 0.2)
        self.assertTrue(collapsed_dup._rail_state()[1])
        self.assertTrue(collapsed_dup._filters_panel.isHidden())

    def test_duplicate_while_a_compute_is_pending_copies_what_the_rail_shows(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        first.gene_query.setText("blaOXA")
        first._apply_rail_filters()
        # No waiting: window 1 has not applied the edit yet.
        dup = controller.duplicate_explore(first)
        self.windows.append(dup)
        self._settle(dup)
        self._apply(dup)
        self._apply(first)
        self.assertEqual(dup.gene_query.text(), "blaOXA")
        self.assertEqual(dup.ctx.applied, first.ctx.applied)
        self.assertEqual(dup.ctx.applied.gene, "blaOXA")
        self.assertEqual(dup.ctx.scheduler.compute_count, 1)

    def test_window_menu_entry_raises_its_window(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        second = self._open(controller, source=first)
        menu = first._window_menu
        menu.aboutToShow.emit()
        entry = next(a for a in menu.actions() if a.text() == second.windowTitle())
        with mock.patch.object(second, "activateWindow") as activate, mock.patch.object(second, "raise_") as raise_:
            entry.trigger()
        activate.assert_called_once_with()
        raise_.assert_called_once_with()

    def test_new_windows_never_land_on_an_open_window(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        area = QtCore.QRect(0, 0, 4000, 3000)
        first.move(area.topLeft())
        with _large_screen():
            second = self._open(controller, source=first)
            self.assertEqual(second.pos(), QtCore.QPoint(CASCADE_OFFSET, CASCADE_OFFSET))
            # Window 2 in the bottom-right corner: the next one cannot cascade
            # from it, and restarts at the top left without covering window 1.
            second.move(area.right() - second.width() + 1, area.bottom() - second.height() + 1)
            _spin(self.app, 0.1)
            third = self._open(controller, source=second)
        self.assertNotEqual(third.pos(), first.pos())
        self.assertNotEqual(third.pos(), second.pos())
        self.assertTrue(area.contains(third.pos()))

    def test_session_windows_without_saved_geometry_cascade(self) -> None:
        controller = self._controller()
        controller.layouts.save_session(controller.run_key(), [1, 2, 3])
        first = self._first(controller)
        windows = controller.windows.explore_windows()
        self.windows += windows
        self.assertEqual(len(windows), 3)
        for win in windows[1:]:
            self._settle(win)
        positions = [w.pos() for w in windows]
        self.assertIs(windows[0], first)
        self.assertEqual(len(set((p.x(), p.y()) for p in positions)), 3)
        self.assertEqual(positions[1], positions[0] + QtCore.QPoint(CASCADE_OFFSET, CASCADE_OFFSET))

    # ------------------------------------------------------------ independence (D7)

    def test_two_contexts_never_deliver_to_each_others_panels(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        second = self._open(controller, source=first)
        self.assertIsNot(first.ctx, second.ctx)
        self.assertIsNot(first.ctx.bus, second.ctx.bus)
        self.assertIsNot(first.filter_rail, second.filter_rail)
        self.assertIs(first.engine, second.engine)  # shared read only, not linked
        recorder = _BusRecorder(second.ctx)
        computes = second.ctx.scheduler.compute_count
        frames_two = {f.uid: f for f in second._panel_manager.frames()}
        rows_before = second.tree.topLevelItemCount()
        # Filters, a subset, hover, selection and an open request in window 1.
        first.gene_query.setText("blaOXA")
        first._set_sample_subset_filter(label="Three", samples=["S1", "S2", "S3"])
        self._apply(first)
        first.ctx.hover(["S1"], origin="sample_tree")
        first.ctx.flush_hover()
        first.ctx.select(first.ctx.selection.__class__(sample="S1", contig="1"), origin="sample_tree")
        first._panel_manager.open_new_instance("geo")
        _spin(self.app, 0.3)
        self.assertEqual(recorder.events, [])
        self.assertEqual(second.ctx.scheduler.compute_count, computes)
        self.assertEqual(second.ctx.applied, self._clean())
        self.assertEqual(second.ctx.hover_state.samples, ())
        self.assertIsNone(second.ctx.selection.sample)
        self.assertEqual(second.hover_status.text(), "")
        self.assertEqual(second.tree.topLevelItemCount(), rows_before)
        self.assertEqual({f.uid for f in second._panel_manager.frames()}, set(frames_two))
        # And the other way round.
        # The sample search box of window 1 too.
        first.sample_query.setText("S1")
        self._apply(first)
        _spin(self.app, 0.2)
        self.assertEqual(recorder.events, [])
        self.assertEqual(second.sample_query.text(), "")
        self.assertEqual(second.ctx.applied, self._clean())
        # And the other way round.
        recorder_one = _BusRecorder(first.ctx, owner="recorder-1")
        hover_one = first.ctx.hover_state
        selection_one = first.ctx.selection
        second.ctx.hover(["S4"], origin="sample_tree")
        second.ctx.flush_hover()
        second.ctx.select(second.ctx.selection.__class__(sample="S4", contig="1"), origin="sample_tree")
        second.sample_query.setText("S4")
        second.product_query.setText("transposase")
        self._apply(second)
        _spin(self.app, 0.2)
        self.assertEqual(recorder_one.events, [])
        self.assertIs(first.ctx.hover_state, hover_one)
        self.assertIs(first.ctx.selection, selection_one)
        self.assertEqual(first.ctx.selection.sample, "S1")
        self.assertEqual(first.sample_query.text(), "S1")
        self.assertEqual(first.ctx.applied.product, "")
        # Each window's own bus delivered its own events.
        self.assertTrue(any(e == bus_events.RESULT_APPLIED for e, _p in recorder.events))
        self.assertTrue(any(e == bus_events.HOVER_CHANGED for e, _p in recorder.events))

    def test_closing_one_window_leaves_the_other_untouched(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        second = self._open(controller, source=first)
        first.gene_query.setText("blaOXA")
        self._apply(first)
        applied, result = first.ctx.applied, first.ctx.result
        frames = [(f.uid, first._panel_manager.is_open(f.uid)) for f in first._panel_manager.frames()]
        groups = first._panel_manager.host.groups()
        second.ctx.hover(["S2"], origin="sample_tree")
        second.ctx.flush_hover()
        computes = first.ctx.scheduler.compute_count
        second.close()
        _spin(self.app, 0.2)
        self.assertFalse(_alive(second))
        self.assertIs(first.ctx.applied, applied)
        self.assertIs(first.ctx.result, result)
        self.assertEqual(first.ctx.scheduler.compute_count, computes)
        self.assertEqual([(f.uid, first._panel_manager.is_open(f.uid)) for f in first._panel_manager.frames()], frames)
        self.assertEqual(first._panel_manager.host.groups(), groups)
        self.assertEqual(first.gene_query.text(), "blaOXA")
        # Window 1 still works.
        recorder = _BusRecorder(first.ctx)
        first.ctx.hover(["S3"], origin="annotations")
        first.ctx.flush_hover()
        self.assertEqual(first.hover_status.text(), "Hovered: S3")
        self.assertTrue(any(e == bus_events.HOVER_CHANGED for e, _p in recorder.events))
        self.assertEqual(controller.windows.explore_windows(), [first])

    def test_panels_never_move_between_windows(self) -> None:
        # D2: no API moves a panel to another window, and no menu offers it.
        qt_names = set(dir(QtCore.QObject))
        for owner in (WindowManager, AppController):
            names = [n for n in dir(owner) if ("move" in n.lower() or "take" in n.lower()) and n not in qt_names]
            self.assertEqual(names, [], owner.__name__)
        controller = self._controller()
        first = self._first(controller)
        manager = first._panel_manager
        tab_menu = manager.build_tab_menu("geo")
        group_menu = manager.build_group_menu("geo")

        def texts(menu: QtWidgets.QMenu) -> List[str]:
            out: List[str] = []
            for action in menu.actions():
                out.append(action.text())
                if action.menu() is not None:
                    out += texts(action.menu())
            return out

        for text in texts(tab_menu) + texts(group_menu):
            self.assertNotIn("window", text.lower())
            self.assertNotIn("float", text.lower())
        self.assertFalse(hasattr(manager.host, "move_to_window"))

    # ------------------------------------------------------------ persistence

    def test_each_window_saves_into_the_slot_of_its_number(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        second = self._open(controller, source=first)
        first._panel_manager.close_panel("geo")
        second._panel_manager.close_panel("phylogeny")
        _spin(self.app, 0.1)
        first._save_explore_layout()
        second._save_explore_layout()
        key = controller.run_key()
        one = controller.layouts.load_window(key, 1)
        two = controller.layouts.load_window(key, 2)
        self.assertIsNotNone(one)
        self.assertIsNotNone(two)
        self.assertNotEqual(one.dock_state, two.dock_state)
        # A window restored from a slot gets that slot's layout.
        restored = self._open(controller, layout=2)
        self.assertFalse(restored._panel_manager.is_open("phylogeny"))
        self.assertTrue(restored._panel_manager.is_open("geo"))

    def test_the_autosave_timer_writes_the_window_s_own_slot(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        second = self._open(controller, source=first)
        key = controller.run_key()
        first._save_explore_layout()
        slot_one = controller.layouts.load_window(key, 1)
        self.assertIsNone(controller.layouts.load_window(key, 2))
        second._panel_manager.close_panel("phylogeny")
        _spin(self.app, 0.1)
        second._layout_autosave_timer.timeout.emit()
        two = controller.layouts.load_window(key, 2)
        self.assertIsNotNone(two)
        self.assertEqual(two.dock_state, second._layout_copy().dock_state)
        self.assertEqual(controller.layouts.load_window(key, 1), slot_one)

    def test_quit_reopens_every_window_with_its_layout_but_fresh_filters(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        second = self._open(controller, source=first)
        clade = self._open(controller, window_subset=SampleSubset(("S1", "S2"), label="Clade 1", kind="clade"))
        first.gene_query.setText("blaOXA")
        self._apply(first)
        second._panel_manager.close_panel("geo")
        second._panel_manager.panel("annotations").set_view_mode("Global")
        _spin(self.app, 0.1)
        quit_action = next(a for a in first._file_menu.actions() if a.objectName() == "QuitAction")
        self.assertTrue(quit_action.shortcut().isEmpty())
        quit_action.trigger()
        _spin(self.app, 0.3)
        for win in (first, second, clade):
            self.assertFalse(_alive(win))
        self.assertEqual(controller.windows.explore_windows(), [])
        self.assertEqual(controller.layouts.load_session(controller.run_key()), [1, 2])
        # The run opens again: two windows (the clade window is not reopened),
        # with their layouts and a clean rail each (D5).
        again = self._controller()
        one = self._first(again)
        windows = again.windows.explore_windows()
        self.windows += windows
        self.assertEqual(len(windows), 2)
        two = windows[1]
        self._settle(two)
        self.assertEqual((one.window_number, two.window_number), (1, 2))
        self.assertEqual(one.ctx.applied, self._clean())
        self.assertEqual(one.gene_query.text(), "")
        self.assertTrue(one._panel_manager.is_open("geo"))
        self.assertFalse(two._panel_manager.is_open("geo"))
        self.assertEqual(two._panel_manager.panel("annotations").view_mode(), "Global")
        self.assertIsNone(two.window_subset)

    def test_closing_windows_one_by_one_reopens_the_last_one(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        second = self._open(controller, source=first)
        second._panel_manager.close_panel("phylogeny")
        _spin(self.app, 0.1)
        first.close()
        _spin(self.app, 0.1)
        second.close()
        _spin(self.app, 0.1)
        self.assertEqual(controller.layouts.load_session(controller.run_key()), [2])
        # Reopened as window 1 with window 2's layout.
        again = self._controller()
        one = self._first(again)
        self.assertEqual(len(again.windows.explore_windows()), 1)
        self.assertEqual(one.window_number, 1)
        self.assertFalse(one._panel_manager.is_open("phylogeny"))

    def test_a_clade_window_closed_last_reopens_window_one(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        first._panel_manager.close_panel("geo")
        clade = self._open(controller, window_subset=SampleSubset(("S1",), label="Clade 1", kind="clade"))
        first.close()
        _spin(self.app, 0.1)
        clade.close()
        _spin(self.app, 0.1)
        self.assertEqual(controller.layouts.load_session(controller.run_key()), [])
        again = self._controller()
        one = self._first(again)
        self.assertIsNone(one.window_subset)
        self.assertFalse(one._panel_manager.is_open("geo"))  # slot 1

    def test_the_app_quitting_with_windows_open_saves_them(self) -> None:
        app = self.app
        before = app.quitOnLastWindowClosed()
        self.addCleanup(app.setQuitOnLastWindowClosed, before)
        with mock.patch.object(QtCore.QCoreApplication, "quit"):
            controller = self._controller(manage_quit=True)
            first = self._first(controller)
            second = self._open(controller, source=first)
            second._panel_manager.close_panel("geo")
            _spin(self.app, 0.1)
            controller._on_about_to_quit()
        key = controller.run_key()
        self.assertEqual(controller.layouts.load_session(key), [1, 2])
        self.assertIsNotNone(controller.layouts.load_window(key, 2))
        restored = self._open(self._controller(), layout=2)
        self.assertFalse(restored._panel_manager.is_open("geo"))

    def test_new_window_and_close_commands(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        second = controller.commands.run(cmd.WINDOW_NEW, first)
        self.windows.append(second)
        self._settle(second)
        self.assertEqual(second._start_layout, START_DEFAULT_LAYOUT)
        close_action = next(a for a in second._file_menu.actions() if a.objectName() == "CloseWindowAction")
        self.assertTrue(close_action.shortcut().isEmpty())
        close_action.trigger()
        _spin(self.app, 0.1)
        self.assertFalse(_alive(second))
        self.assertEqual(controller.windows.explore_windows(), [first])
        for name in (cmd.WINDOW_NEW, cmd.WINDOW_DUPLICATE, cmd.WINDOW_OPEN_SUBSET, cmd.WINDOW_CLOSE, cmd.APP_QUIT):
            self.assertTrue(controller.commands.has(name), name)
        # No QShortcut anywhere in the window (D8) beyond the neighbourhood
        # panel's own, which this window does not have open.
        self.assertEqual(first.findChildren(QtGui.QShortcut), [])


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class NeighbourhoodWindowTests(_RealWindowCase):
    """``window.open_neighbourhood``: the comparison in a window of its own."""

    def _controller(self) -> AppController:
        return AppController(self._run_data())

    def _first(self, controller: AppController) -> Any:
        win = controller.open_first_window()
        self.windows.append(win)
        self._settle(win)
        return win

    def _apply(self, win: Any) -> None:
        self.assertTrue(_wait(self.app, lambda: not win.ctx.is_busy(), 20))
        _spin(self.app, 0.1)

    def _compare(self, controller: AppController, origin: Any, name: str = "lysR", **kwargs: Any) -> Any:
        before = list(controller.windows.explore_windows())
        kwargs.setdefault("rows", [])
        result = controller.commands.run(
            cmd.WINDOW_OPEN_NEIGHBOURHOOD,
            origin,
            title=f"Gene neighborhoods: {name}",
            source="gene_card",
            key=("gene", name),
            **kwargs,
        )
        if result not in before:
            self.windows.append(result)
            self._settle(result)
            self._apply(result)
        return result

    def test_it_copies_every_filter_the_sort_and_the_rail(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        predictor = first.run.predictor_options[0]
        first.sample_query.setText("s")
        first.gene_query.setText("blaOXA")
        first.product_query.setText("lactamase")
        first.molecule_filter_box.setCurrentIndex(first.molecule_filter_box.findData("chromosome"))
        first.annotations_panel.hide_hypo_chk.setChecked(True)
        first._panel_manager.panel("sample_tree").hide_empty_contigs_chk.setChecked(False)
        model = first._predictor_model
        for row in range(model.rowCount()):
            if model.item(row).data(QtCore.Qt.ItemDataRole.UserRole) == predictor:
                model.item(row).setCheckState(QtCore.Qt.CheckState.Checked)
        first.metadata_filters = {"Country": ["Italy", "Romania"]}
        first._render_metadata_filters()
        first._apply_rail_filters()
        first._set_sample_subset_filter(label="Five", samples=["S1", "S2", "S3", "S4", "S5"])
        sort_box, order_box = first.sample_sort_column_box, first.sample_sort_order_box
        sort_box.setCurrentIndex(sort_box.findData("Country"))
        order_box.setCurrentIndex(order_box.findData("desc"))
        first._apply_rail_state(360, False)
        self._apply(first)
        _spin(self.app, 0.3)
        source_filters = first.ctx.draft
        new = self._compare(controller, first)
        self.assertIsNot(new, first)
        self.assertEqual(new.ctx.applied, source_filters)
        self.assertEqual(FilterState.from_widgets(new), source_filters)
        self.assertEqual(new.sample_query.text(), "s")
        self.assertEqual(new.gene_query.text(), "blaOXA")
        self.assertEqual(new.product_query.text(), "lactamase")
        self.assertEqual(new.molecule_filter_box.currentData(), "chromosome")
        self.assertTrue(new.ctx.applied.hide_hypothetical)
        self.assertFalse(new.ctx.applied.hide_empty_contigs)
        self.assertEqual(new._selected_predictors(), [predictor])
        self.assertEqual(new.metadata_filters, {"Country": ["Italy", "Romania"]})
        self.assertEqual(new.ctx.draft.subset.label, "Five")
        self.assertEqual(new.ctx.sample_set(), first.ctx.sample_set())
        self.assertEqual(new.ctx.scheduler.compute_count, 1)  # one compute for the copy
        self.assertEqual(new.sample_sort_column_box.currentData(), "Country")
        self.assertEqual(new.sample_sort_order_box.currentData(), "desc")
        _spin(self.app, 0.3)
        self.assertEqual(new._rail_state(), (360, False))
        # Only the comparison in its docking area; nothing of window 1's layout.
        manager = new._panel_manager
        viewers = _document_panels(new, "neighbourhood")
        self.assertEqual(len(viewers), 1)
        self.assertEqual(manager.host.open_uids(), [viewers[0].uid])
        self.assertEqual(manager.closed_stack(), [])
        self.assertTrue(first._panel_manager.is_open("annotations"))
        self.assertEqual(_document_panels(first, "neighbourhood"), [])
        # The panel follows its window's sample set.
        self.assertEqual(viewers[0]._scope[0], tuple(new.ctx.sample_set()))

    def test_the_two_windows_are_independent_afterwards(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        first.gene_query.setText("blaOXA")
        self._apply(first)
        narrowed = first.ctx.sample_set()
        self.assertLess(len(narrowed), len(first.run.meta_df))
        new = self._compare(controller, first)
        self.assertEqual(new.ctx.sample_set(), narrowed)
        first_computes = first.ctx.scheduler.compute_count
        new.gene_query.setText("")
        self._apply(new)
        self.assertEqual(len(new.ctx.sample_set()), len(first.run.meta_df))
        self.assertEqual(first.ctx.sample_set(), narrowed)
        self.assertEqual(first.gene_query.text(), "blaOXA")
        self.assertEqual(first.ctx.scheduler.compute_count, first_computes)
        viewer = _document_panels(new, "neighbourhood")[0]
        self.assertTrue(_wait(self.app, lambda: len(viewer._scope[0]) == len(first.run.meta_df)))
        # And the other way round.
        new_computes = new.ctx.scheduler.compute_count
        first._set_sample_subset_filter(label="Two", samples=["S1", "S2"])
        self._apply(first)
        self.assertIsNone(new.ctx.draft.subset)
        self.assertEqual(new.ctx.scheduler.compute_count, new_computes)

    def test_it_is_an_ordinary_explore_window(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        new = self._compare(controller, first)
        manager = new._panel_manager
        viewer = _document_panels(new, "neighbourhood")[0]
        self.assertEqual(new.windowTitle(), f"BACTolomew - {first.run.output_dir.name} [2] - lysR neighborhoods")
        # Panels open into it like into any window (no search box of their own).
        controller.commands.run(cmd.PANEL_OPEN, new, kind="annotations")
        _spin(self.app, 0.2)
        self.assertTrue(manager.is_open("annotations"))
        self.assertTrue(manager.is_open(viewer.uid))
        self.assertEqual(manager.panel("annotations").findChildren(QtWidgets.QLineEdit), [])
        # Reset Layout builds Classic and keeps the comparison.
        controller.commands.run(cmd.LAYOUT_RESET, new)
        _spin(self.app, 0.2)
        self.assertTrue(manager.is_open("sample_tree"))
        self.assertTrue(manager.is_open(viewer.uid))
        # A second comparison from it opens another window, cloned from it.
        new.gene_query.setText("blaOXA")
        third = self._compare(controller, new, name="dnaA")
        self.assertNotIn(third, (first, new))
        self.assertEqual(third.ctx.applied.gene, "blaOXA")
        self.assertEqual(len(_document_panels(new, "neighbourhood")), 1)
        # Asking for an open comparison again shows its window: no new one.
        self.assertIs(self._compare(controller, first), new)
        self.assertIs(self._compare(controller, third), new)
        self.assertEqual(len(controller.windows.explore_windows()), 3)
        # The title names the comparison while it is open.
        manager.close_panel(viewer.uid)
        _spin(self.app, 0.1)
        self.assertEqual(new.windowTitle(), f"BACTolomew - {first.run.output_dir.name} [2]")
        self.assertTrue(manager.reopen_closed())
        _spin(self.app, 0.1)
        self.assertTrue(new.windowTitle().endswith(" - lysR neighborhoods"))

    def test_two_cards_of_one_gene_get_told_apart_across_windows(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        one = controller.commands.run(
            cmd.WINDOW_OPEN_NEIGHBOURHOOD, first, rows=[], title="Gene neighborhoods: lysR",
            source="gene_card", key=("gene", "a"), reference_sample="S1",
        )
        two = controller.commands.run(
            cmd.WINDOW_OPEN_NEIGHBOURHOOD, first, rows=[], title="Gene neighborhoods: lysR",
            source="gene_card", key=("gene", "b"), reference_sample="S2",
        )
        self.windows += [one, two]
        self.assertIsNot(one, two)
        self.assertEqual(_document_panels(one, "neighbourhood")[0].title(), "lysR neighborhoods")
        self.assertEqual(_document_panels(two, "neighbourhood")[0].title(), "lysR neighborhoods (S2)")
        self.assertTrue(two.windowTitle().endswith(" - lysR neighborhoods (S2)"))

    def test_a_locked_window_subset_stays_locked(self) -> None:
        controller = self._controller()
        clade = controller.open_explore(window_subset=SampleSubset(("S1", "S2"), label="Clade 1", kind="clade"))
        self.windows.append(clade)
        self._settle(clade)
        new = self._compare(controller, clade)
        self.assertEqual(new.window_subset, clade.window_subset)
        self.assertEqual(new.ctx.applied.window_subset, clade.window_subset)
        self.assertTrue(new.windowTitle().endswith(" - Clade 1 - lysR neighborhoods"))

    def test_it_is_not_reopened_with_the_run(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        first._apply_rail_state(300, False)
        new = self._compare(controller, first)
        self.assertFalse(controller._reopens(new))
        new._save_explore_layout()
        last = controller.layouts.load_last_window()
        self.assertTrue(last is None or "neighbourhood" not in str(last.panels))
        first._save_explore_layout()
        new._save_explore_layout()  # does not replace the last saved window
        self.assertEqual(controller.layouts.load_last_window().panels, first._last_saved_layout.panels)
        controller.quit_all()
        _spin(self.app, 0.2)
        self.assertEqual(controller.layouts.load_session(controller.run_key()), [1])
        # Closed last, it leaves an empty session (window 1 reopens).
        again = self._controller()
        one = self._first(again)
        self.assertEqual(len(again.windows.explore_windows()), 1)
        other = self._compare(again, one)
        one.close()
        _spin(self.app, 0.1)
        other.close()
        _spin(self.app, 0.1)
        self.assertEqual(again.layouts.load_session(again.run_key()), [])

    def test_it_never_writes_a_window_slot(self) -> None:
        # Window numbers are the lowest free one, so a neighbourhood window
        # can be [1]: it must not replace window 1's saved layout.
        controller = self._controller()
        first = self._first(controller)
        first._save_explore_layout()
        key = first._layout_run_key
        saved = controller.layouts.load_window(key, 1)
        self.assertIsNotNone(saved)
        second = self._compare(controller, first)
        second._save_explore_layout()
        self.assertIsNone(controller.layouts.load_window(key, second.window_number))
        first.close()
        _spin(self.app, 0.1)
        third = self._compare(controller, second, name="dnaA")
        self.assertEqual(third.window_number, 1)
        third._save_explore_layout()
        third.close()
        _spin(self.app, 0.1)
        kept = controller.layouts.load_window(key, 1)
        self.assertIsNotNone(kept)
        self.assertEqual(kept.panels, saved.panels)
        self.assertEqual(kept.rail_width, saved.rail_width)

    def test_a_panel_opened_from_the_picker_splits_beside_the_comparison(self) -> None:
        # The picker's default Open (no target, not to the side) puts a
        # panel beside the comparison: as a tab it would hide the rows whose
        # hover it shows (a Phylogeny in the real app; this run has no tree,
        # so the Map, which takes the hover too, stands in).
        controller = self._controller()
        first = self._first(controller)
        new = self._compare(controller, first)
        manager = new._panel_manager
        viewer = _document_panels(new, "neighbourhood")[0]
        manager.host.focus_view(viewer.uid)
        _spin(self.app, 0.1)
        self.assertEqual(manager.open_new_instance("geo"), "geo")
        _spin(self.app, 0.2)
        self.assertEqual(sorted(manager.host.groups()), sorted([["geo"], [viewer.uid]]))
        self.assertTrue(manager.host.is_current_tab("geo"))
        self.assertTrue(manager.host.is_current_tab(viewer.uid))
        # Heavy, so built once the window's heavy panels start (700 ms).
        self.assertTrue(_wait(self.app, lambda: manager.panel("geo") is not None))
        geo = manager.panel("geo")
        self.assertLess(geo.mapTo(new, QtCore.QPoint(0, 0)).x(), viewer.mapTo(new, QtCore.QPoint(0, 0)).x())
        # With an ordinary panel active, the rule is unchanged: a tab.
        manager.host.focus_view("geo")
        _spin(self.app, 0.1)
        self.assertEqual(manager.open_new_instance("annotations"), "annotations")
        _spin(self.app, 0.2)
        self.assertIn("geo", manager.host.group_uids("annotations"))
        self.assertNotIn(viewer.uid, manager.host.group_uids("annotations"))

    def test_show_only_in_window_adds_a_missing_sample_tree(self) -> None:
        controller = self._controller()
        first = self._first(controller)
        new = self._compare(controller, first)
        manager = new._panel_manager
        viewer = _document_panels(new, "neighbourhood")[0]
        controller.commands.run(cmd.PANEL_OPEN, new, kind="annotations")
        _spin(self.app, 0.2)
        self.assertTrue(manager.is_open("annotations"))
        self.assertFalse(manager.is_open("sample_tree"))
        # A tab of the comparison's group (a default open splits beside it).
        manager.host.move_view("annotations", Placement("tab", relative_to=viewer.uid))
        _spin(self.app, 0.1)
        manager.host.focus_view(viewer.uid)
        _spin(self.app, 0.1)
        self.assertTrue(manager.host.is_current_tab(viewer.uid))
        self.assertFalse(manager.host.is_current_tab("annotations"))
        viewer._open_samples_in_window(["S1", "S2"])
        _spin(self.app, 0.2)
        # Annotations needs the Sample tree; the comparison stays in front.
        self.assertTrue(manager.is_open("sample_tree"))
        self.assertTrue(manager.host.is_current_tab(viewer.uid))
        self.assertEqual(len(manager.uids("annotations")), 1)
        self.assertEqual(new.ctx.draft.subset.samples, ("S1", "S2"))

    def test_a_failed_open_falls_back_to_classic(self) -> None:
        from annotation_explorer_polars.windows.explore_window import ExploreWindow

        controller = self._controller()
        first = self._first(controller)
        with mock.patch.object(ExploreWindow, "_open_neighbourhood_view", return_value=None):
            new = self._compare(controller, first)
        self.assertIsNot(new, first)
        self.assertTrue(new._panel_manager.is_open("annotations"))
        self.assertTrue(new.isVisible())

    def test_the_empty_start_layout_builds_nothing(self) -> None:
        controller = self._controller()
        win = controller.open_explore(layout=START_EMPTY_LAYOUT)
        self.windows.append(win)
        self._settle(win)
        self.assertEqual(win._panel_manager.host.open_uids(), [])
        self.assertFalse(win._filters_panel.isHidden())
        self.assertIsNotNone(win.ctx.result)


class PanTabForwardingTests(unittest.TestCase):
    """The pan views' own widgets up to PanCoreAccessoryTab.subsetWindowRequested."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def test_the_pan_tab_forwards_tree_and_plasmidome_requests_with_their_kind(self) -> None:
        import tempfile
        from pathlib import Path

        import annotation_explorer_polars.panaroo_core_accessory as pca

        try:
            from tests.test_tree_view import TREE_TEXT
        except ImportError:
            from test_tree_view import TREE_TEXT

        class _Plasmidome(QtWidgets.QWidget):
            samplesHovered = QtCore.pyqtSignal(object)
            sampleSubsetRequested = QtCore.pyqtSignal(object)
            sampleSubsetWindowRequested = QtCore.pyqtSignal(object)
            contextRequested = QtCore.pyqtSignal(object)

            def __init__(self, *_args: Any, parent: Any = None, **_kwargs: Any) -> None:
                super().__init__(parent)

            def set_sample_filter(self, _selected: Any) -> None:
                pass

        with tempfile.TemporaryDirectory() as directory:
            tree_path = Path(directory) / "tree.nwk"
            tree_path.write_text(TREE_TEXT, encoding="utf-8")
            tab = pca.PanCoreAccessoryTab(
                Path(directory),
                None,
                reference_sample=None,
                tree_file=tree_path,
                preserve_branch_lengths=True,
                show_info_label=False,
            )
            self.addCleanup(tab.close)
            got: List[Any] = []
            tab.subsetWindowRequested.connect(lambda payload, kind: got.append((kind, payload)))
            clade = {"name": "7", "samples": ["A", "B"]}
            cluster = {"name": "3", "samples": ["C"]}
            tab.tree.cladeWindowRequested.emit(clade)
            tab.tree.clusterWindowRequested.emit(cluster)
            # A plasmidome view built by the tab (as a Plasmidome panel's is).
            with mock.patch.object(pca, "PlasmidomePanel", _Plasmidome):
                view = tab._build_view("plasmidome", tab)
            self.assertIsInstance(view, _Plasmidome)
            plasmid = {"name": "pX", "samples": ["E", "F"]}
            view.sampleSubsetWindowRequested.emit(plasmid)
            self.assertEqual(got, [("clade", clade), ("cluster", cluster), ("plasmid", plasmid)])

    def test_right_click_on_a_clade_label_opens_the_menu_and_a_left_click_filters(self) -> None:
        try:
            from tests.test_tree_view import CLUSTERS, GROUP_COLORS, GROUPS, TREE_TEXT
        except ImportError:
            from test_tree_view import CLUSTERS, GROUP_COLORS, GROUPS, TREE_TEXT
        from PyQt6 import QtTest

        from annotation_explorer_polars.tree_view import _SCENE_DATA_KIND, _SCENE_DATA_NAME

        view = RaxmlTreeView(
            layout_mode="rectangular",
            optimize_groups=False,
            use_group_backgrounds=True,
            show_group_labels=True,
            show_cluster_annotations=True,
        )
        self.addCleanup(view.close)
        view.resize(1000, 700)
        view.set_leaf_groups(GROUPS, group_colors=GROUP_COLORS)
        view.set_leaf_clusters(CLUSTERS)
        view.set_tree_text(TREE_TEXT)
        view.show()
        self.app.processEvents()
        got: List[Any] = []
        view.cladeWindowRequested.connect(lambda payload: got.append(("clade-window", payload["name"])))
        view.clusterWindowRequested.connect(lambda payload: got.append(("cluster-window", payload["name"])))
        view.cladeClicked.connect(lambda payload: got.append(("clade-filter", payload["name"])))

        def label_point(kind: str) -> tuple:
            item = next(i for i in view.scene().items() if i.data(_SCENE_DATA_KIND) == kind)
            # Labels ignore the view's zoom: map through their device transform.
            point = item.deviceTransform(view.viewportTransform()).map(item.boundingRect().center())
            return item.data(_SCENE_DATA_NAME), point.toPoint()

        shown: List[List[str]] = []

        def fake_exec(menu: QtWidgets.QMenu, *_args: Any) -> None:
            shown.append([a.text() for a in menu.actions()])
            menu.actions()[1].trigger()

        for kind in ("clade-label", "cluster-label"):
            name, point = label_point(kind)
            event = QtGui.QContextMenuEvent(
                QtGui.QContextMenuEvent.Reason.Mouse, point, view.viewport().mapToGlobal(point)
            )
            with mock.patch.object(QtWidgets.QMenu, "exec", new=fake_exec):
                QtWidgets.QApplication.sendEvent(view.viewport(), event)
            self.app.processEvents()
            noun = "clade" if kind == "clade-label" else "cluster"
            self.assertEqual(got[-1], (f"{noun}-window", name))
            self.assertEqual(shown[-1][1], f"Open {noun} in new explore window")
        # A left click on a clade label still filters this window, no new window.
        name, point = label_point("clade-label")
        count = len(got)
        QtTest.QTest.mouseClick(view.viewport(), QtCore.Qt.MouseButton.LeftButton, pos=point)
        self.app.processEvents()
        self.assertEqual(got[count:], [("clade-filter", name)])


if __name__ == "__main__":
    unittest.main()
