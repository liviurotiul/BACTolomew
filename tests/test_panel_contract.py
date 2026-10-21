"""Panel protocol and multi-instance panels in the real window (window plan, Phase 4).

Plan 8.5: ``tests/test_panel_contract.py``, parametrised over the registry
(create, bind, show, hide, dispose, JSON state round trip), plus the Phase 4
behavior of the window: multi-instance panels (map, neoPan, a second
Annotations panel with its own view mode, D3) with Duplicate Panel, Split
Right and Move to Group, hover through the context in both directions (the
Pan tree included), open requests, ``set_subset`` and the status bar's active-panel
line (3.9). Offscreen on the synthetic run of the golden tests; skipped when
PyQt6Ads is not installed.
"""
from __future__ import annotations

import gc
import json
import os
import shutil
import tempfile
import time
import unittest
import weakref
from typing import Any, List, Optional
from unittest import mock

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

from PyQt6 import QtCore, QtGui, QtWidgets

from annotation_explorer_polars.docking import docking_available

try:
    from tests.test_dock_host import _app, _Explorer, _spin
    from tests._pan_fixture import write_pan_outputs
except ImportError:
    from test_dock_host import _app, _Explorer, _spin
    from _pan_fixture import write_pan_outputs

if docking_available():
    from annotation_explorer_polars.docking.dock_host import Placement
    from annotation_explorer_polars.panels.annotations import AnnotationsPanel
    from annotation_explorer_polars.panels.alignment import AlignmentDocument, AlignmentPanel
    from annotation_explorer_polars.panels.base import PAN_MDS, PAN_VIEW_KINDS, ExplorePanel
    from annotation_explorer_polars.panels.neighbourhood import NeighbourhoodDocument
    from annotation_explorer_polars.session.requests import OpenPanel, params_tuple


def _neighbourhood_rows(win: Any, annotation_id: str = "S1:1:1") -> List[dict]:
    """Gene neighborhood rows of the synthetic run, as a gene card collects them."""
    from annotation_explorer_polars.annotation_detail_dialog import AnnotationDetailDialog

    row = win.run.ann_df.loc[win.run.ann_df["AnnotationID"] == annotation_id].iloc[0].to_dict()
    card = AnnotationDetailDialog(None, row, list(row), ann_df=win.run.ann_df)
    try:
        return card._build_context_collection_rows(window_genes=30)
    finally:
        card.deleteLater()


def _open_document(win: Any, kind: str) -> str:
    """Open a document panel of ``kind`` (neighbourhood or alignment) with data."""
    manager = win._panel_manager
    if kind == "neighbourhood":
        document: Any = NeighbourhoodDocument(
            title="Gene neighborhoods: lysR", rows=_neighbourhood_rows(win), source="gene_card", key=("test", kind)
        )
    else:
        document = AlignmentDocument(
            alignment_text=">a|S1|1|g\nMKA\n>b|S2|1|g\nMKA\n>c|S5|1|g\nMRA\n", title="Alignment: g"
        )
    uid = manager.open_document(kind, document)
    assert uid is not None
    return uid


def _menu_action(menu: QtWidgets.QMenu, path: str) -> Any:
    parts = path.split(" > ")
    for part in parts[:-1]:
        menu = next(a.menu() for a in menu.actions() if a.text() == part)
    return next(a for a in menu.actions() if a.text() == parts[-1])


class _WindowCase(unittest.TestCase):
    """A real window per test (QSettings isolated; WebEngine off)."""

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
        self.output_dir = tempfile.mkdtemp(prefix="bactolomew_panels_run_")
        self.windows: List[Any] = []
        self._clear_saved_layouts()

    def tearDown(self) -> None:
        for win in self.windows:
            win.close()
            win.deleteLater()
        _spin(self.app, 0.05)
        shutil.rmtree(self.output_dir, ignore_errors=True)
        self._clear_saved_layouts()

    @staticmethod
    def _clear_saved_layouts() -> None:
        settings = QtCore.QSettings()
        settings.remove("explore")
        settings.remove("annotationExplorer")
        settings.sync()

    def _window(self, show: bool = True) -> Any:
        win = _Explorer.build(self.output_dir)
        self.windows.append(win)
        if show:
            _Explorer.show_and_settle(self.app, win)
        return win

    def _idle(self, win: Any, extra: float = 0.3) -> None:
        end = time.monotonic() + 10
        while time.monotonic() < end and (
            win.ctx.scheduler.is_running() or win.ctx.scheduler.debounce_timer.isActive() or win._panels.has_pending_delivery()
        ):
            _spin(self.app, 0.02)
        _spin(self.app, extra)

    def _wait_for(self, predicate: Any, timeout: float = 3.0) -> bool:
        end = time.monotonic() + timeout
        while time.monotonic() < end:
            if predicate():
                return True
            _spin(self.app, 0.02)
        return bool(predicate())

    def _search_gene(self, win: Any, text: str) -> None:
        win.gene_query.setText(text)
        self._idle(win)

    def _duplicate_annotations(self, win: Any, where: str = "right") -> Any:
        uid = win._panel_manager.duplicate_panel("annotations", where)
        self.assertIsNotNone(uid)
        _spin(self.app, 0.2)
        panel = win._panel_manager.panel(uid)
        self.assertIsInstance(panel, AnnotationsPanel)
        return panel

    @staticmethod
    def _set_mode(panel: Any, mode: str) -> None:
        box = panel.mode_box()
        box.setCurrentIndex(box.findData(mode))

    def _collect(self) -> None:
        # PyQt drops the slots of deleted senders on a later deferred delete.
        for _ in range(3):
            QtCore.QCoreApplication.sendPostedEvents(None, QtCore.QEvent.Type.DeferredDelete.value)
            _spin(self.app, 0.02)
        gc.collect()

    @staticmethod
    def _mouse(widget: QtWidgets.QWidget, kind: Any, pos: Optional[QtCore.QPoint] = None) -> None:
        """Send a real MouseMove (or a Leave) to a viewport, as the pointer would."""
        if kind == QtCore.QEvent.Type.Leave:
            QtWidgets.QApplication.sendEvent(widget, QtCore.QEvent(QtCore.QEvent.Type.Leave))
            return
        local = QtCore.QPointF(pos)
        event = QtGui.QMouseEvent(
            QtCore.QEvent.Type.MouseMove,
            local,
            QtCore.QPointF(widget.mapToGlobal(pos)),
            QtCore.Qt.MouseButton.NoButton,
            QtCore.Qt.MouseButton.NoButton,
            QtCore.Qt.KeyboardModifier.NoModifier,
        )
        QtWidgets.QApplication.sendEvent(widget, event)


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class PanelContractTests(_WindowCase):
    """Every kind of the registry, through the panel protocol (plan 4.8).

    The run folder has pangenome outputs (tests/_pan_fixture.py), so the
    standalone pan panels of Phase 9 are available and go through the
    contract with the rest."""

    def setUp(self) -> None:
        super().setUp()
        write_pan_outputs(self.output_dir)

    def test_every_kind_follows_the_panel_contract(self) -> None:
        win = self._window()
        manager = win._panel_manager
        self._search_gene(win, "blaOXA")
        for spec in manager.registry.specs():
            kind = spec.kind
            with self.subTest(kind=kind):
                # Create and bind (lazily for heavy kinds), show. A document
                # kind is opened with its data and has a suffixed uid.
                if spec.document:
                    uid = _open_document(win, kind)
                    self.assertNotEqual(uid, kind)
                else:
                    self.assertTrue(manager.open_panel(kind))
                    uid = kind
                _spin(self.app, 0.4)
                frame = manager.frame(uid)
                self.assertTrue(manager.is_visible(uid))
                self.assertTrue(frame.is_built())
                panel = frame.panel()
                self.assertIsInstance(panel, ExplorePanel)
                self.assertEqual((panel.KIND, panel.uid), (kind, uid))
                self.assertIs(panel.ctx, win.ctx)
                self.assertEqual(panel.SINGLE_INSTANCE, spec.single_instance)
                self.assertEqual(panel.HEAVY, spec.heavy)
                self.assertEqual(panel.NEEDS_PARAMS, spec.document)
                self.assertFalse(frame.is_dirty())
                # JSON state round trip (panel params of the layout store).
                state = frame.save_state()
                frame.restore_state(json.loads(json.dumps(state)))
                self.assertEqual(frame.save_state(), state)
                # Hidden: a heavy panel goes dirty and catches up when shown.
                # (A document panel is disposed on close: hide it behind
                # another tab of its group instead.)
                if spec.document:
                    manager.host.move_view(uid, Placement("tab", relative_to="sample_tree"))
                    manager.host.set_current_tab("sample_tree")
                else:
                    manager.close_panel(kind)
                _spin(self.app, 0.1)
                self.assertFalse(manager.is_visible(uid))
                self._search_gene(win, "ompA" if win.gene_query.text() == "blaOXA" else "blaOXA")
                self.assertEqual(frame.is_dirty(), spec.heavy)
                manager.open_panel(uid)
                _spin(self.app, 0.4)
                self.assertFalse(frame.is_dirty())
        self.assertEqual(manager.host.floating_count(), 0)

    def test_disposed_panels_get_nothing_more(self) -> None:
        win = self._window()
        manager = win._panel_manager
        for kind in ("cgmlst_force", *PAN_VIEW_KINDS, PAN_MDS):
            manager.open_panel(kind)
            _spin(self.app, 0.3)
        frames = {frame.uid: frame for frame in manager.frames()}
        panels = {uid: frame.panel() for uid, frame in frames.items()}
        self.assertTrue(all(panels.values()))
        for uid, frame in frames.items():
            frame.dispose()
            self.assertEqual(win.ctx.bus.subscriptions(uid), [], uid)
            self.assertIsNone(panels[uid].ctx, uid)
        self.assertEqual(sorted(win.ctx.bus.owners()), ["rail", "window"])
        failures: List[str] = []
        win.ctx.bus.handlerFailed.connect(lambda owner, event, text: failures.append(owner))
        self._search_gene(win, "blaOXA")
        self.assertEqual(failures, [])

    def test_every_kind_unbinds_and_disposes(self) -> None:
        win = self._window()
        manager = win._panel_manager
        self._search_gene(win, "blaOXA")
        uids = {}
        for spec in manager.registry.specs():
            if spec.document:
                uids[spec.kind] = _open_document(win, spec.kind)
            else:
                manager.open_panel(spec.kind)
                uids[spec.kind] = spec.kind
            _spin(self.app, 0.3)
        for spec in manager.registry.specs():
            kind = spec.kind
            uid = uids[kind]
            with self.subTest(kind=kind):
                frame = manager.frame(uid)
                panel = frame.panel()
                self.assertIsNotNone(panel)
                self.assertTrue(win.ctx.bus.subscriptions(uid))
                frame.dispose()
                self.assertTrue(frame.is_disposed())
                self.assertEqual(win.ctx.bus.subscriptions(uid), [])
                self.assertIsNone(panel.ctx)
                frame.dispose()  # twice is harmless
        self.assertEqual(sorted(win.ctx.bus.owners()), ["rail", "window"])
        failures: List[str] = []
        win.ctx.bus.handlerFailed.connect(lambda owner, event, text: failures.append(owner))
        win.ctx.hover(["S1"], origin="rail")
        win.ctx.flush_hover()
        _spin(self.app, 0.1)
        self.assertEqual(failures, [])

    def test_closed_extra_instances_of_every_multi_instance_kind_do_not_leak(self) -> None:
        win = self._window()
        manager = win._panel_manager
        manager.open_panel("phylogeny")
        _spin(self.app, 0.4)
        for kind in [spec.kind for spec in manager.registry.specs() if not spec.single_instance and not spec.document]:
            with self.subTest(kind=kind):
                if not manager.is_open(kind):
                    # The pan panels are not in Classic: open the first one.
                    self.assertTrue(manager.open_panel(kind))
                    _spin(self.app, 0.4)
                uid = manager.duplicate_panel(kind, "right")
                self.assertIsNotNone(uid)
                _spin(self.app, 0.6)
                panel = manager.panel(uid)
                self.assertIsNotNone(panel)
                self.assertTrue(win.ctx.bus.subscriptions(uid))
                ref_panel = weakref.ref(panel)
                ref_frame = weakref.ref(manager.frame(uid))
                del panel
                manager.close_panel(uid)
                _spin(self.app, 0.1)
                self.assertIsNone(manager.frame(uid))
                self.assertEqual(win.ctx.bus.subscriptions(uid), [])
                self._collect()
                self.assertIsNone(ref_panel(), uid)
                self.assertIsNone(ref_frame(), uid)
        self.assertEqual(manager.host.floating_count(), 0)


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class MultiInstancePanelTests(_WindowCase):
    def test_a_second_annotations_panel_has_its_own_view_mode(self) -> None:
        win = self._window()
        manager = win._panel_manager
        self._search_gene(win, "blaOXA")
        dup = self._duplicate_annotations(win)
        self.assertNotEqual(dup.uid, "annotations")
        self.assertTrue(dup.uid.startswith("annotations-"))
        # Two open: each title says what it shows, the same mode numbered.
        self.assertEqual(manager.title("annotations"), "Annotations · Contig")
        self.assertEqual(manager.title(dup.uid), "Annotations · Contig 2")
        self.assertEqual(manager.host.neighbor_group("annotations", "right"), dup.uid)
        self.assertIs(dup.table.model(), dup.model)  # drew the current result at once
        computes = win.ctx.scheduler.compute_count
        result = win.ctx.result
        status = win.status.text()
        self._set_mode(dup, "Global")
        # D3 per instance: no compute, the other panel and the shared bar unchanged.
        self.assertEqual(dup.view_mode(), "Global")
        self.assertEqual(manager.host.view_title(dup.uid), "Annotations · All matches")
        self.assertEqual(manager.title("annotations"), "Annotations · Contig")
        self.assertEqual(win.view_mode, "Contig")
        self.assertEqual(dup.table.model().rowCount(), 4)
        self.assertEqual(dup.status_text(), "4 annotations across 3 samples")
        self.assertEqual(win.annotations_panel.status_text(), "")
        self.assertEqual(win.ctx.scheduler.compute_count, computes)
        self.assertIs(win.ctx.result, result)
        self.assertEqual(win.status.text(), status)
        # A selection loads the contig in the panels showing a contig only
        # (Open annotations switches the Annotations panel used last).
        manager.host.focus_view("annotations")
        win._open_tree_contig_payload({"sample": "S1", "contig": "3", "contig_label": "Contig 3"}, mode="Contig")
        self.assertEqual(win.annotations_panel.status_text(), "Sample S1 | Contig 3 | 1 hit, 2 annotation rows")
        self.assertEqual(dup.status_text(), "4 annotations across 3 samples")
        self._set_mode(dup, "Contig")
        self.assertEqual(dup.status_text(), "Sample S1 | Contig 3 | 1 hit, 2 annotation rows")
        # A row of an All matches duplicate selects its contig for the others.
        self._set_mode(dup, "Global")
        model = dup.table.model()
        row = next(r for r in range(model.rowCount()) if model.row_data(r).get("Sample") == "S4")
        dup.table.setCurrentIndex(model.index(row, 0))
        self.assertEqual(win.ctx.selection.sample, "S4")
        self.assertTrue(win.annotations_panel.status_text().startswith("Sample S4 |"))
        self.assertEqual(win.ctx.scheduler.compute_count, computes)

    def test_instance_titles_follow_the_detail_and_keep_custom_titles(self) -> None:
        win = self._window()
        manager = win._panel_manager
        self.assertEqual(manager.title("annotations"), "Annotations")
        dup = self._duplicate_annotations(win)
        self._set_mode(dup, "Context")
        self.assertEqual(manager.title(dup.uid), "Annotations · Around matches")
        # A title set with Rename is never overwritten, and saved alone.
        manager.rename_panel("annotations", "Main table")
        self._set_mode(win.annotations_panel, "Global")
        self.assertEqual(manager.host.view_title("annotations"), "Main table")
        saved = {entry["uid"]: entry["params"] for entry in manager.panel_list()}
        self.assertEqual(saved["annotations"].get("title"), "Main table")
        self.assertNotIn("title", saved[dup.uid])
        # Back to one open instance: the plain title, also when the one left
        # open is the extra instance (the closed first one stays registered).
        manager.rename_panel("annotations", "")
        manager.close_panel("annotations")
        _spin(self.app, 0.1)
        self.assertEqual(manager.host.view_title(dup.uid), "Annotations")
        manager.open_panel("annotations")
        _spin(self.app, 0.1)
        manager.close_panel(dup.uid)
        _spin(self.app, 0.1)
        self.assertEqual(manager.title("annotations"), "Annotations")
        self.assertEqual(manager.host.view_title("annotations"), "Annotations")

    def test_tab_menu_split_right_duplicate_and_move_to_group(self) -> None:
        win = self._window()
        manager = win._panel_manager
        menu = manager.build_tab_menu("geo")
        self.assertTrue(_menu_action(menu, "Duplicate Panel").isEnabled())
        _menu_action(menu, "Split Right").trigger()
        _spin(self.app, 0.4)
        maps = manager.uids("geo")
        self.assertEqual(len(maps), 2)
        second = maps[1]
        self.assertEqual(manager.host.group_uids(second), [second])
        self.assertTrue(manager.is_visible("geo") and manager.is_visible(second))
        self.assertFalse(manager.frame(second).is_dirty())
        # Both maps show the one sample set; the payload is computed once.
        self._search_gene(win, "blaOXA")
        self.assertFalse(manager.frame("geo").is_dirty() or manager.frame(second).is_dirty())
        self.assertIn(win.ctx.result.samples, win.engine._map_cache)
        # Duplicate Panel adds a tab to the same group; Move to Group moves it.
        _menu_action(manager.build_tab_menu(second), "Duplicate Panel").trigger()
        _spin(self.app, 0.3)
        third = manager.uids("geo")[2]
        self.assertEqual(manager.host.group_uids(second), [second, third])
        _menu_action(manager.build_tab_menu(third), "Move to Group > Left").trigger()
        _spin(self.app, 0.2)
        self.assertIn(third, manager.host.group_uids("geo"))
        # Single-instance kinds have no Duplicate Panel; their Split Right moves.
        texts = [a.text() for a in manager.build_tab_menu("sample_tree").actions()]
        self.assertNotIn("Duplicate Panel", texts)
        self.assertEqual(manager.host.floating_count(), 0)

    def test_two_maps_compute_the_payload_once(self) -> None:
        win = self._window()
        manager = win._panel_manager
        _menu_action(manager.build_tab_menu("geo"), "Split Right").trigger()
        _spin(self.app, 0.4)
        second = manager.uids("geo")[1]
        engine = win.engine
        with mock.patch.object(engine, "_compute_map_payload", wraps=engine._compute_map_payload) as compute:
            self._search_gene(win, "blaOXA")
            self.assertFalse(manager.frame("geo").is_dirty() or manager.frame(second).is_dirty())
            keys = [tuple(c.args[0]) for c in compute.call_args_list]
        self.assertEqual(keys.count(tuple(win.ctx.result.samples)), 1)

    def test_split_down_opens_a_copy_below(self) -> None:
        write_pan_outputs(self.output_dir)
        win = self._window()
        manager = win._panel_manager
        manager.open_panel("phylogeny")
        _spin(self.app, 0.4)
        _menu_action(manager.build_tab_menu("phylogeny"), "Split Down").trigger()
        _spin(self.app, 0.6)
        uids = manager.uids("phylogeny")
        self.assertEqual(len(uids), 2)
        copy = uids[1]
        self.assertTrue(manager.is_open("phylogeny") and manager.is_open(copy))
        self.assertEqual(manager.host.group_uids(copy), [copy])
        self.assertEqual(manager.host.neighbor_group(copy, "top"), "phylogeny")
        self.assertEqual(manager.host.floating_count(), 0)

    def test_a_reopened_extra_panel_gets_its_own_split_back(self) -> None:
        win = self._window()
        manager = win._panel_manager
        dup = self._duplicate_annotations(win, "right")
        uid = dup.uid
        self._set_mode(dup, "Global")
        del dup
        self.assertEqual(manager.host.group_uids(uid), [uid])
        manager.host._docks[uid].requestCloseDockWidget()
        _spin(self.app, 0.1)
        self.assertIsNone(manager.frame(uid))
        self.assertEqual(manager.reopen_closed(), uid)
        _spin(self.app, 0.2)
        self.assertEqual(manager.host.group_uids(uid), [uid])
        self.assertEqual(manager.host.neighbor_group("annotations", "right"), uid)
        self.assertEqual(manager.panel(uid).view_mode(), "Global")

    def test_open_panel_requests_open_a_new_instance_with_its_options(self) -> None:
        win = self._window()
        manager = win._panel_manager
        self._search_gene(win, "blaOXA")
        with mock.patch.object(win, "_load_global", wraps=win._load_global) as load_global:
            win.ctx.request(
                OpenPanel("annotations", params_tuple({"view_mode": "Global"}), beside="annotations"),
                origin="sample_tree",
            )
            _spin(self.app, 0.2)
        uids = manager.uids("annotations")
        self.assertEqual(len(uids), 2)
        panel = manager.panel(uids[1])
        self.assertEqual(panel.view_mode(), "Global")
        self.assertEqual(win.view_mode, "Contig")
        self.assertEqual(manager.host.group_uids(uids[1]), ["annotations", uids[1]])
        self.assertEqual(panel.status_text(), "4 annotations across 3 samples")
        # Drawn once, in its own mode (not first in Contig).
        self.assertEqual(load_global.call_count, 1)
        # A single-instance kind is opened (or focused).
        win.ctx.request(OpenPanel("cgmlst_force"), origin="sample_tree")
        _spin(self.app, 0.4)
        self.assertTrue(manager.is_visible("cgmlst_force"))

    def test_open_in_context_view_draws_the_table_once(self) -> None:
        win = self._window()
        self._search_gene(win, "blaOXA")
        self.assertEqual(win.view_mode, "Contig")
        with mock.patch.object(win, "_load_annotations", wraps=win._load_annotations) as load:
            win._open_tree_contig_payload({"sample": "S3", "contig": "2"}, mode="Context")
        self.assertEqual(win.view_mode, "Context")
        self.assertEqual(load.call_count, 1)
        self.assertTrue(win.annotations_panel.status_text().startswith("Sample S3 | 2 |"))

    def test_a_second_phylogeny_builds_its_own_view(self) -> None:
        write_pan_outputs(self.output_dir)
        win = self._window()
        manager = win._panel_manager
        uid = manager.duplicate_panel("phylogeny", "right")
        _spin(self.app, 1.0)
        panel = manager.panel(uid)
        self.assertIsNotNone(panel)
        self.assertIsNotNone(panel.tab)
        first = manager.panel("phylogeny")
        self.assertTrue(first is None or first.tab is not panel.tab)
        self.assertFalse(manager.frame(uid).is_dirty())

    def test_closing_an_extra_instance_disposes_it_and_reopen_brings_it_back(self) -> None:
        win = self._window()
        manager = win._panel_manager
        dup = self._duplicate_annotations(win, "tab")
        uid = dup.uid
        self._set_mode(dup, "Global")
        manager.rename_panel(uid, "Hits")
        ref = weakref.ref(dup)
        del dup
        _menu_action(manager.build_tab_menu(uid), "Close").trigger()
        _spin(self.app, 0.1)
        self.assertIsNone(manager.frame(uid))
        self.assertFalse(manager.host.has_view(uid))
        self.assertEqual(win.ctx.bus.subscriptions(uid), [])
        self.assertEqual(manager.closed_stack()[-1], uid)
        # PyQt drops the slots of deleted senders on a later deferred delete.
        for _ in range(3):
            QtCore.QCoreApplication.sendPostedEvents(None, QtCore.QEvent.Type.DeferredDelete.value)
            _spin(self.app, 0.02)
        gc.collect()
        self.assertIsNone(ref())
        # Reopen Closed Panel rebuilds it beside where it was, with its options.
        self.assertEqual(manager.reopen_closed(), uid)
        _spin(self.app, 0.2)
        again = manager.panel(uid)
        self.assertEqual(again.view_mode(), "Global")
        self.assertEqual(manager.title(uid), "Hits")
        self.assertEqual(manager.host.group_uids(uid), ["annotations", uid])
        self.assertTrue(again.status_text().endswith("samples") or again.status_text().startswith("Enter"))
        # The tab's own close button (ADS close request) disposes it as well.
        manager.host._docks[uid].requestCloseDockWidget()
        _spin(self.app, 0.1)
        self.assertIsNone(manager.frame(uid))
        self.assertEqual(manager.host.floating_count(), 0)

    def test_extra_instances_persist_with_the_layout_and_reset_closes_them(self) -> None:
        win = self._window()
        dup = self._duplicate_annotations(win)
        self._set_mode(dup, "Global")
        uid = dup.uid
        groups = sorted(sorted(g) for g in win._panel_manager.host.groups())
        win.close()
        _spin(self.app, 0.1)
        again = self._window()
        manager = again._panel_manager
        self.assertEqual(sorted(sorted(g) for g in manager.host.groups()), groups)
        restored = manager.panel(uid)
        self.assertIsInstance(restored, AnnotationsPanel)
        self.assertEqual(restored.view_mode(), "Global")
        self.assertEqual(restored.mode_box().currentData(), "Global")
        self.assertEqual(again.view_mode, "Contig")
        manager.reset_layout()
        _spin(self.app, 0.2)
        self.assertEqual(manager.uids("annotations"), ["annotations"])
        self.assertEqual(again.ctx.bus.subscriptions(uid), [])

    def test_the_active_annotations_panel_line_is_appended_to_the_status_bar(self) -> None:
        win = self._window()
        manager = win._panel_manager
        self._search_gene(win, "blaOXA")
        dup = self._duplicate_annotations(win)
        self._set_mode(dup, "Global")
        manager.host.focus_view(dup.uid)
        _spin(self.app, 0.1)
        self.assertEqual(win.active_panel_status.text(), "| 4 annotations across 3 samples")
        self.assertFalse(win.active_panel_status.isHidden())
        stats = win.status.text()
        self._set_mode(dup, "Contig")
        self.assertEqual(win.active_panel_status.text(), "")
        self.assertTrue(win.active_panel_status.isHidden())
        self.assertEqual(win.status.text(), stats)  # the applied stats never change with a view mode
        manager.host.focus_view("geo")
        _spin(self.app, 0.1)
        self.assertTrue(win.active_panel_status.isHidden())

    def test_open_in_context_view_shows_the_annotations_panel_used_last(self) -> None:
        win = self._window()
        manager = win._panel_manager
        manager.host.focus_view("annotations")
        manager.host.move_view("annotations", Placement("tab", relative_to="sample_tree"))
        manager.host.set_current_tab("sample_tree")
        _spin(self.app, 0.2)
        self.assertFalse(manager.is_visible("annotations"))
        computes = win.ctx.scheduler.compute_count
        win._open_tree_contig_payload({"sample": "S3", "contig": "2"}, mode="Context")
        _spin(self.app, 0.2)
        self.assertTrue(manager.is_visible("annotations"))
        self.assertEqual(win.view_mode, "Context")
        self.assertEqual(win.annotations_panel.status_text(), "Sample S3 | 2 | 2 rows shown")
        self.assertEqual(win.ctx.scheduler.compute_count, computes)
        # With a second panel used last, that one switches and the first stays.
        dup = self._duplicate_annotations(win)
        manager.host.focus_view(dup.uid)
        win._open_tree_contig_payload({"sample": "S1", "contig": "3"}, mode="Context")
        self.assertEqual(dup.view_mode(), "Context")
        self.assertTrue(dup.status_text().startswith("Sample S1 | 3 |"))


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class LinkingThroughTheContextTests(_WindowCase):
    """Hover, subsets and open requests go through the context (plan 4.5, 4.6)."""

    def test_tree_hover_reaches_the_other_visible_panels_and_not_back(self) -> None:
        win = self._window()
        manager = win._panel_manager
        self._search_gene(win, "blaOXA")
        dup = self._duplicate_annotations(win)
        self._set_mode(dup, "Global")
        tree_panel = manager.panel("sample_tree")
        with mock.patch.object(win._panel_manager.panel("geo").view, "highlight_samples") as map_highlight, \
                mock.patch.object(tree_panel, "show_hover") as tree_inbound:
            tree_panel._hover_from_tree("S1")
            self.assertEqual(win._tree_hover_shown, ("S1",))  # the tree's own feedback, at once
            win.ctx.flush_hover()
            _spin(self.app, 0.1)
            self.assertEqual(win.ctx.hover_state.samples, ("S1",))
            self.assertEqual(win.ctx.hover_state.origin, "sample_tree")
            self.assertEqual(win.hover_status.text(), "Hovered: S1")
            map_highlight.assert_called_with(["S1"])
            tree_inbound.assert_not_called()  # no echo
        self.assertEqual(dup.model.hover_samples(), frozenset({"S1"}))
        # The table's own hover goes out the same way, with its uid as origin.
        dup._hover_from_table("S3")
        win.ctx.flush_hover()
        self.assertEqual(win.ctx.hover_state.origin, dup.uid)
        self.assertEqual(win._tree_hover_shown, ("S3",))
        dup._hover_from_table(None)
        win.ctx.flush_hover()
        self.assertEqual(win._tree_hover_shown, ())

    def test_hover_handed_from_the_tree_to_a_table_row_of_the_same_sample_keeps_the_tree_bold(self) -> None:
        win = self._window()
        manager = win._panel_manager
        self._search_gene(win, "blaOXA")
        dup = self._duplicate_annotations(win)
        self._set_mode(dup, "Global")
        tree_panel = manager.panel("sample_tree")
        tree_panel._hover_from_tree("S1")
        win.ctx.flush_hover()
        # Off the tree and onto a table row of S1 within one coalescing window.
        tree_panel._hover_from_tree(None)
        self.assertEqual(win._tree_hover_shown, ())
        dup._hover_from_table("S1")
        win.ctx.flush_hover()
        self.assertEqual(win.ctx.hover_state.origin, dup.uid)
        self.assertEqual(win._tree_hover_shown, ("S1",))

    def test_hover_comes_from_the_viewport_event_filters_and_survives_a_move(self) -> None:
        win = self._window()
        manager = win._panel_manager
        self._search_gene(win, "blaOXA")
        tree = win.tree
        item = next(
            tree.topLevelItem(i)
            for i in range(tree.topLevelItemCount())
            if win._sample_from_tree_item(tree.topLevelItem(i))
        )
        sample = win._sample_from_tree_item(item)

        def hover_tree() -> None:
            _spin(self.app, 0.1)
            self._mouse(tree.viewport(), QtCore.QEvent.Type.MouseMove, tree.visualItemRect(item).center())
            win.ctx.flush_hover()
            self.assertEqual(win.ctx.hover_state.samples, (sample,))
            self.assertEqual(win.ctx.hover_state.origin, "sample_tree")
            self._mouse(tree.viewport(), QtCore.QEvent.Type.Leave)
            win.ctx.flush_hover()
            self.assertEqual(win.ctx.hover_state.samples, ())

        hover_tree()
        # Moved to another group: the panel's own filter still works.
        self.assertTrue(manager.host.move_view("sample_tree", Placement("right", relative_to="geo")))
        hover_tree()
        # The table's viewport filter (All matches lists several samples).
        dup = self._duplicate_annotations(win)
        self._set_mode(dup, "Global")
        _spin(self.app, 0.1)
        table = dup.table
        model = table.model()
        row = next(r for r in range(model.rowCount()) if model.row_data(r).get("Sample") == "S3")
        self._mouse(table.viewport(), QtCore.QEvent.Type.MouseMove, table.visualRect(model.index(row, 0)).center())
        win.ctx.flush_hover()
        self.assertEqual((win.ctx.hover_state.samples, win.ctx.hover_state.origin), (("S3",), dup.uid))
        self.assertEqual(win._tree_hover_shown, ("S3",))
        manager.host.move_view(dup.uid, Placement("tab", relative_to="geo"))
        _spin(self.app, 0.2)
        self._mouse(table.viewport(), QtCore.QEvent.Type.Leave)
        win.ctx.flush_hover()
        self.assertEqual(win.ctx.hover_state.samples, ())

    def test_cgmlst_force_hover_works_both_ways(self) -> None:
        win = self._window()
        manager = win._panel_manager
        # Beside the map, so the sample tree (its hover partner here) stays shown.
        self.assertTrue(manager.open_panel("cgmlst_force", placement=Placement("right", relative_to="geo")))
        _spin(self.app, 0.5)
        self.assertTrue(manager.is_visible("sample_tree"))
        self.assertTrue(manager.is_visible("cgmlst_force"))
        panel = manager.panel("cgmlst_force")
        self.assertIsNotNone(panel)
        plot = panel.plot
        # Out.
        plot.sampleHovered.emit("S2")
        win.ctx.flush_hover()
        self.assertEqual((win.ctx.hover_state.samples, win.ctx.hover_state.origin), (("S2",), "cgmlst_force"))
        self.assertEqual(win._tree_hover_shown, ("S2",))
        # In (coalesced for a heavy panel), one sample and several.
        with mock.patch.object(plot, "highlight_sample") as one, mock.patch.object(plot, "highlight_samples") as several:
            win.ctx.hover(["S4"], origin="sample_tree")
            win.ctx.flush_hover()
            self.assertTrue(self._wait_for(lambda: one.called))
            one.assert_called_with("S4")
            win.ctx.hover(["S1", "S5"], origin="sample_tree")
            win.ctx.flush_hover()
            self.assertTrue(self._wait_for(lambda: several.called))
            several.assert_called_with(["S1", "S5"])

    def test_a_maximized_group_hides_the_others_for_delivery(self) -> None:
        win = self._window()
        manager = win._panel_manager
        self.assertTrue(manager.is_visible("geo"))
        self.assertTrue(manager.host.maximize_group("annotations"))
        _spin(self.app, 0.2)
        self.assertFalse(manager.is_visible("geo"))
        with mock.patch.object(manager.panel("geo"), "show_hover") as geo_hover:
            win.ctx.hover(["S2"], origin="sample_tree")
            win.ctx.flush_hover()
            _spin(self.app, 0.1)
            geo_hover.assert_not_called()
        self._search_gene(win, "blaOXA")
        self.assertTrue(manager.frame("geo").is_dirty())
        manager.host.end_maximize()
        _spin(self.app, 0.5)
        self.assertTrue(manager.is_visible("geo"))
        self.assertFalse(manager.frame("geo").is_dirty())

    def test_a_move_redraws_a_panel_at_most_once(self) -> None:
        win = self._window()
        manager = win._panel_manager
        self._search_gene(win, "blaOXA")
        panel = manager.panel("geo")
        with mock.patch.object(panel, "apply_result", wraps=panel.apply_result) as redraw:
            self.assertTrue(manager.host.move_view("geo", Placement("right", relative_to="annotations")))
            _spin(self.app, 0.5)
            self.assertLessEqual(redraw.call_count, 1)
            self._search_gene(win, "ompA")
            self.assertEqual(redraw.call_count, 1)  # the new result, once
        self.assertFalse(manager.frame("geo").is_dirty())

    def test_phylogeny_hover_works_both_ways(self) -> None:
        write_pan_outputs(self.output_dir)
        win = self._window()
        manager = win._panel_manager
        # Classic has it in the bottom row, beside the sample tree.
        self.assertTrue(manager.open_panel("phylogeny"))
        _spin(self.app, 0.5)
        tab = win._panel_manager.panel("phylogeny").tab
        self.assertIsNotNone(tab)
        self.assertTrue(tab._enable_tree_hover_interactions)
        # Out: a subtree hovers several samples at once.
        tab.samplesHovered.emit(["S1", "S3"])
        win.ctx.flush_hover()
        self.assertEqual(win.ctx.hover_state.samples, ("S1", "S3"))
        self.assertEqual(win.ctx.hover_state.origin, "phylogeny")
        self.assertEqual(win._tree_hover_shown, ("S1", "S3"))
        self.assertEqual(win.hover_status.text(), "Hovered: 2 samples")
        # In: another panel's hover reaches the visible Phylogeny.
        self.assertTrue(manager.is_visible("phylogeny"))
        with mock.patch.object(tab, "highlight_sample") as one, mock.patch.object(tab, "highlight_samples") as several:
            win.ctx.hover(["S4"], origin="sample_tree")
            win.ctx.flush_hover()
            # Heavy panels get hover coalesced (30 ms); the tree may be busy
            # redrawing its pruned tree, so wait for the delivery.
            self.assertTrue(self._wait_for(lambda: one.called))
            one.assert_called_with("S4")
            win.ctx.hover(["S2", "S5"], origin="sample_tree")
            win.ctx.flush_hover()
            self.assertTrue(self._wait_for(lambda: several.called))
            several.assert_called_with(["S2", "S5"])

    def test_phylogeny_shares_every_hovered_sample(self) -> None:
        write_pan_outputs(self.output_dir)
        win = self._window()
        manager = win._panel_manager
        manager.open_panel("phylogeny")
        _spin(self.app, 0.5)
        tab = win._panel_manager.panel("phylogeny").tab
        # A subtree of the tree hovers three samples (no hover hub in the
        # tab since Phase 9: it goes straight to the panel).
        tab.tree.subtreeHovered.emit(["S1", "S2", "S6"])
        win.ctx.flush_hover()
        self.assertEqual(win.ctx.hover_state.samples, ("S1", "S2", "S6"))
        self.assertEqual(win.ctx.hover_state.origin, "phylogeny")
        self.assertEqual(win.ctx.hover_state.kind, "samples")

    def test_hidden_panels_drop_hover_and_catch_up_on_show(self) -> None:
        write_pan_outputs(self.output_dir)
        win = self._window()
        manager = win._panel_manager
        manager.host.move_view("phylogeny", Placement("tab", relative_to="geo"))
        _spin(self.app, 0.5)
        manager.host.set_current_tab("geo")
        _spin(self.app, 0.2)
        self.assertFalse(manager.is_visible("phylogeny"))
        with mock.patch.object(win._panel_manager.panel("phylogeny").tab, "highlight_sample") as tree_highlight:
            win.ctx.hover(["S2"], origin="sample_tree")
            win.ctx.flush_hover()
            _spin(self.app, 0.1)
            tree_highlight.assert_not_called()
            manager.host.set_current_tab("phylogeny")
            _spin(self.app, 0.3)
            # Replayed on show (the tab switch then moves the focus, which
            # clears the hover, as a tab change always did).
            tree_highlight.assert_any_call("S2")

    def test_pan_subset_requests_go_through_set_subset(self) -> None:
        write_pan_outputs(self.output_dir)
        win = self._window()
        changes: List[Any] = []
        origins: List[str] = []
        win.ctx.subsetChanged.connect(lambda subset, origin: (changes.append(subset), origins.append(origin)))
        win._panel_manager.open_panel("phylogeny")
        _spin(self.app, 0.5)
        # A clade click in the tree: the panel narrows the set, as its own uid.
        win._panel_manager.panel("phylogeny").tab.cladeFilterRequested.emit({"name": "A", "samples": ["S1", "S2"]})
        self._idle(win)
        self.assertEqual(changes[-1].samples, ("S1", "S2"))
        self.assertEqual(origins[-1], "phylogeny")
        self.assertEqual(win.ctx.applied.subset.label, "Clade A")
        self.assertEqual(win.ctx.result.samples, ("S1", "S2"))
        self.assertEqual(list(win.ctx.draft.subset.samples), ["S1", "S2"])
        labels = [w.text() for w in win.filters_container.findChildren(QtWidgets.QLabel)]
        self.assertIn("Clade A | n=2", labels)
        win._clear_sample_subset_filter()
        self._idle(win)
        self.assertIsNone(win.ctx.applied.subset)
        self.assertEqual(len(win.ctx.result.samples), 6)

    def test_open_requests_are_carried_out_by_the_window(self) -> None:
        from annotation_explorer_polars.session.requests import OpenGeneCard, OpenSampleCard

        win = self._window()
        with mock.patch.object(win, "_open_sample_card_for_sample") as sample_card, \
                mock.patch.object(win, "_open_gene_card") as gene_card:
            win.ctx.request(OpenSampleCard("S2"), origin="phylogeny")
            win.ctx.request(OpenGeneCard.from_row({"Sample": "S1", "Gene": "blaOXA"}), origin="track")
        sample_card.assert_called_once_with("S2")
        gene_card.assert_called_once_with({"Sample": "S1", "Gene": "blaOXA"})
        # A gene card's alignment result is a panel of this window; its row
        # hover goes into this window's context with the panel's uid.
        row = {"Sample": "S1", "Contig": "3", "Start": "1", "Stop": "9", "Strand": "+", "Gene": "blaOXA"}
        card = win._open_gene_card(row)  # closed with the window in tearDown
        with mock.patch.object(win, "_open_alignment_view") as alignment:
            card._show_alignment_result(">a|S1|3|blaOXA\nMA\n", "done")
        self.assertEqual(alignment.call_args.kwargs["title"], "Alignment: blaOXA")
        uid = win._open_alignment_view(alignment_text=">a|S5|3|blaOXA\nMA\n", status_text="done")
        panel = win._panel_manager.panel(uid)
        self.assertIsInstance(panel, AlignmentPanel)
        panel.view._set_hover_sample("S5")
        win.ctx.flush_hover()
        self.assertEqual(win.ctx.hover_state.samples, ("S5",))
        self.assertEqual(win.ctx.hover_state.origin, uid)

    def test_a_slow_compute_shows_the_busy_indicator(self) -> None:
        win = self._window()
        real_compute = win.engine.compute

        def slow(filters: Any) -> Any:
            time.sleep(0.6)
            return real_compute(filters)

        with mock.patch.object(win.engine, "compute", side_effect=slow):
            win.gene_query.setText("lysR")
            _spin(self.app, 0.6)
            self.assertFalse(win.busy_status.isHidden())
            self.assertEqual(win.busy_status.text(), "Updating...")
            self._idle(win)
        self.assertTrue(win.busy_status.isHidden())


if __name__ == "__main__":
    unittest.main()
