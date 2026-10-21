"""ExploreContext, its apply scheduler and the D3 view mode change (window plan, Phase 3).

The scheduler tests use a fake engine on a private thread pool, so they can
hold a compute and check single flight, latest wins and stale drops. The
window tests use the golden synthetic run (tests/test_refresh_parity.py).
"""
from __future__ import annotations

import os
import threading
import time
import unittest
from types import SimpleNamespace
from typing import Any, List, Optional
from unittest import mock

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
    from tests.test_refresh_parity import _ExplorerFixture, _wait_idle
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings
    from test_refresh_parity import _ExplorerFixture, _wait_idle

isolate_qt_settings()

import polars as pl
from PyQt6 import QtCore, QtWidgets

from annotation_explorer_polars.model.filter_state import AppliedResult, FilterState, SampleSubset, Selection
from annotation_explorer_polars.session.explore_context import ExploreContext
from annotation_explorer_polars.shell.layout_store import HINT_SAMPLE_SEARCH_D3, HINTS_GROUP


def _result(tag: str) -> AppliedResult:
    return AppliedResult(
        samples=(tag,),
        allowed_samples=(tag,),
        filtered_ann=pl.DataFrame(),
        metadata=pl.DataFrame(),
        contig_records={},
        filtered_counts={},
        stats_text=tag,
    )


class _FakeEngine:
    """Records every compute; ``hold`` makes the next computes wait for ``release``."""

    def __init__(self) -> None:
        self.calls: List[FilterState] = []
        self.gate = threading.Event()
        self.gate.set()
        self.fail_on: Optional[str] = None
        self._lock = threading.Lock()

    def hold(self) -> None:
        self.gate.clear()

    def release(self) -> None:
        self.gate.set()

    def compute(self, filters: FilterState) -> AppliedResult:
        with self._lock:
            self.calls.append(filters)
        self.gate.wait(10)
        if self.fail_on is not None and filters.gene == self.fail_on:
            raise RuntimeError("compute failed on purpose")
        return _result(filters.gene or "none")


class ExploreContextSchedulerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def setUp(self) -> None:
        self.engine = _FakeEngine()
        self.pool = QtCore.QThreadPool()
        self.pool.setMaxThreadCount(2)
        self.ctx = ExploreContext(self.engine, pool=self.pool)
        self.applied: List[tuple] = []
        self.failures: List[str] = []
        self.firsts: List[bool] = []
        self.ctx.resultApplied.connect(lambda f, r: self.applied.append((f, r)))
        self.ctx.applyFailed.connect(self.failures.append)
        self.ctx.firstResultApplied.connect(lambda: self.firsts.append(True))

    def tearDown(self) -> None:
        self.engine.release()
        self.pool.waitForDone(5000)
        self.app.processEvents()

    def _wait(self, condition, timeout: float = 5.0) -> None:
        end = time.monotonic() + timeout
        while time.monotonic() < end:
            self.app.processEvents(QtCore.QEventLoop.ProcessEventsFlag.AllEvents, 20)
            if condition():
                return
            time.sleep(0.005)
        self.fail("timed out")

    def test_three_text_edits_make_one_compute(self) -> None:
        for text in ("b", "bl", "bla"):
            self.ctx.update_draft(FilterState(gene=text), debounce=True)
        self.assertEqual(self.engine.calls, [])
        self._wait(lambda: bool(self.applied) and self.ctx.scheduler.is_idle())
        self.assertEqual([f.gene for f in self.engine.calls], ["bla"])
        self.assertEqual(self.ctx.applied, FilterState(gene="bla"))
        self.assertEqual(self.ctx.result.samples, ("bla",))
        self.assertEqual(self.ctx.scheduler.compute_count, 1)

    def test_a_burst_of_immediate_edits_makes_one_compute(self) -> None:
        # One action that sets several controls (and text typed just before)
        # is one request: immediate edits apply on the next event-loop pass.
        self.ctx.update_draft(FilterState(gene="oxa"), debounce=True)
        self.ctx.update_draft(FilterState(gene="oxa", molecule="plasmid"))
        self.ctx.update_draft(FilterState(gene="oxa", molecule="plasmid", hide_hypothetical=True))
        self.ctx.update_draft(FilterState(gene="oxa", molecule="plasmid", hide_hypothetical=True, product="beta"), debounce=True)
        self.assertEqual(self.engine.calls, [])
        self._wait(lambda: bool(self.applied) and self.ctx.scheduler.is_idle())
        self.assertEqual(self.engine.calls, [FilterState(gene="oxa", molecule="plasmid", hide_hypothetical=True, product="beta")])

    def test_latest_wins_and_stale_results_are_dropped(self) -> None:
        self.engine.hold()
        self.ctx.update_draft(FilterState(gene="first"))
        self._wait(lambda: len(self.engine.calls) == 1)
        # While "first" computes, two more requests: only the last one runs next.
        self.ctx.update_draft(FilterState(gene="second"))
        self.ctx.update_draft(FilterState(gene="third"))
        self.assertTrue(self.ctx.is_busy())
        self.engine.release()
        self._wait(lambda: self.ctx.scheduler.is_idle() and bool(self.applied))
        self.assertEqual([f.gene for f in self.engine.calls], ["first", "third"])
        # "first" finished after newer requests were made: dropped, never applied.
        self.assertEqual([f.gene for f, _r in self.applied], ["third"])
        self.assertEqual(self.ctx.result.stats_text, "third")

    def test_an_unchanged_draft_is_not_computed_again_unless_forced(self) -> None:
        self.ctx.update_draft(FilterState(gene="x"))
        self._wait(lambda: self.ctx.scheduler.is_idle() and bool(self.applied))
        self.ctx.update_draft(FilterState(gene="x"))
        self.ctx.update_draft(FilterState(gene="x"), debounce=True)
        self._wait(lambda: self.ctx.scheduler.is_idle())
        self.assertEqual(len(self.engine.calls), 1)
        self.ctx.apply()
        self._wait(lambda: self.ctx.scheduler.is_idle() and len(self.applied) == 2)
        self.assertEqual(len(self.engine.calls), 2)

    def test_first_result_signal_fires_once_also_on_failure(self) -> None:
        self.engine.fail_on = "broken"
        self.ctx.update_draft(FilterState(gene="broken"))
        self._wait(lambda: bool(self.failures))
        self.assertEqual(self.firsts, [True])
        self.assertIn("compute failed on purpose", self.failures[0])
        self.assertIsNone(self.ctx.result)
        self.assertTrue(self.ctx.first_result_done)
        self.ctx.update_draft(FilterState(gene="fine"))
        self._wait(lambda: bool(self.applied))
        self.assertEqual(self.firsts, [True])
        # A failed request can be asked for again.
        self.engine.fail_on = None
        self.ctx.update_draft(FilterState(gene="broken"))
        self._wait(lambda: len(self.applied) == 2)

    def test_the_pool_thread_runs_only_the_engine(self) -> None:
        threads: List[str] = []
        compute = self.engine.compute

        def recording(filters: FilterState) -> AppliedResult:
            threads.append(threading.current_thread().name)
            return compute(filters)

        self.engine.compute = recording  # type: ignore[method-assign]
        self.ctx.update_draft(FilterState(gene="x"))
        self._wait(lambda: bool(self.applied))
        self.assertNotEqual(threads[0], threading.main_thread().name)

    def test_select_emits_once_per_change(self) -> None:
        seen: List[tuple] = []
        self.ctx.selectionChanged.connect(lambda sel, origin: seen.append((sel, origin)))
        sel = Selection(sample="S1", contig="2", contigs=["2"], contig_label="Contig 2")
        self.ctx.select(sel, origin="tree")
        self.ctx.select(Selection(sample="S1", contig="2", contigs=("2",), contig_label="Contig 2"), origin="tree")
        self.assertEqual(seen, [(sel, "tree")])
        self.assertEqual(self.ctx.selection.contigs, ("2",))


class FilterStateTests(unittest.TestCase):
    def test_normalized_and_order_insensitive(self) -> None:
        a = FilterState(
            sample_query=" Romania ",
            predictors=frozenset({"B", "A"}),
            metadata=(("Country", ("Italy", "Spain")), ("Year", ())),
            molecule="Plasmid",
        )
        b = FilterState(sample_query="Romania", predictors=frozenset({"A", "B"}), metadata={"Country": ["Spain", "Italy"]}, molecule="plasmid")
        self.assertEqual(a, b)
        self.assertEqual(a.key(), b.key())
        self.assertEqual(a.metadata_filters(), {"Country": ["Italy", "Spain"]})
        self.assertEqual(FilterState(molecule="other").molecule, "")
        # An empty subset is no subset; a subset label does not change the query key.
        self.assertIsNone(FilterState(subset=SampleSubset(())).subset)
        self.assertEqual(
            FilterState(subset=SampleSubset(("S1",), label="A")).key(),
            FilterState(subset=SampleSubset(("S1",), label="B")).key(),
        )

    def test_the_sample_search_is_not_an_annotation_filter(self) -> None:
        self.assertFalse(FilterState(sample_query="blaOXA").has_annotation_filter())
        for changes in ({"gene": "x"}, {"product": "x"}, {"molecule": "plasmid"}, {"predictors": frozenset({"p"})},
                        {"tags": frozenset({"A"})}, {"hide_hypothetical": True}):
            with self.subTest(changes=changes):
                self.assertTrue(FilterState(**changes).has_annotation_filter())
        self.assertFalse(FilterState(hide_empty_contigs=True).has_annotation_filter())

    def test_ann_filters_dict(self) -> None:
        state = FilterState(gene="blaOXA", tags=frozenset({"V", "A"}))
        self.assertEqual(
            state.ann_filters(),
            {"gene_text": "blaOXA", "product_text": None, "molecule": "", "predictors": [], "tags": ["A", "V"], "hide_hypothetical": False},
        )
        self.assertEqual(FilterState.from_query(ann_filters=state.ann_filters()), state)

    def test_from_widgets_tolerates_missing_widgets(self) -> None:
        class Host:
            pass

        self.assertEqual(FilterState.from_widgets(Host()), FilterState())

    def test_from_widgets_takes_the_window_filter_flags_from_the_draft(self) -> None:
        # "Hide hypothetical" and "Hide contigs without hits" are not rail
        # widgets (they sit in the Annotations and Sample tree panels): the
        # draft holds them.
        host = SimpleNamespace(ctx=SimpleNamespace(draft=FilterState(hide_hypothetical=True, hide_empty_contigs=True)))
        state = FilterState.from_widgets(host)
        self.assertTrue(state.hide_hypothetical)
        self.assertTrue(state.hide_empty_contigs)


class ExploreWindowContextTests(unittest.TestCase):
    """The window over its context: one sample set, D3 view modes, the hint, sorting."""

    win: Any = None

    @classmethod
    def setUpClass(cls) -> None:
        cls.win = _ExplorerFixture.build(cls)
        cls.win.resize(1400, 900)
        _wait_idle(cls.app, cls.win)

    @classmethod
    def tearDownClass(cls) -> None:
        _ExplorerFixture.teardown(cls, cls.win)

    def setUp(self) -> None:
        win = self.win
        QtCore.QSettings().remove(f"{HINTS_GROUP}/{HINT_SAMPLE_SEARCH_D3}")
        win.sample_query.clear()
        win.gene_query.clear()
        win.product_query.clear()
        win._panel_manager.panel("sample_tree").hide_empty_contigs_chk.setChecked(False)
        win.annotations_panel.hide_hypo_chk.setChecked(False)
        win._clear_sample_subset_filter()
        win.sample_sort_column_box.setCurrentIndex(0)
        win.sample_sort_order_box.setCurrentIndex(0)
        if win.view_mode != "Contig":
            win.annotations_panel.set_view_mode("Contig")
        win.refresh()
        _wait_idle(self.app, win)
        win.sample_search_hint.setVisible(False)
        win._sample_search_hint_query = None

    def _search_gene(self, text: str) -> None:
        self.win.gene_query.setText(text)
        _wait_idle(self.app, self.win)

    def _tree_order(self) -> List[str]:
        tree = self.win.tree
        return [tree.topLevelItem(i).text(0) for i in range(tree.topLevelItemCount())]

    def _tree_samples_with_hits(self) -> List[str]:
        result = self.win.ctx.result
        return [s for s in self._tree_order() if int(result.filtered_counts.get(s, 0)) > 0]

    def _heavy_panel_samples(self) -> dict:
        """What every heavy panel receives for the applied result. The window
        is hidden, so its heavy panel frames hold the result until shown
        (plan 3.9); deliver it to each now and record what the panel gets."""
        from annotation_explorer_polars.panels.cgmlst_force import CgmlstForcePanel
        from annotation_explorer_polars.panels.geo import GeoPanel
        from annotation_explorer_polars.panels.pan import PanMetadataPanel

        win = self.win
        seen: dict = {}

        def record(panel: Any, _filters: Any, result: Any) -> None:
            seen[panel.KIND] = tuple(result.samples)
            if panel.KIND == "geo":
                seen["map"] = win.engine.map_payload(result.samples)

        win._panels.set_ready()
        with mock.patch.object(GeoPanel, "apply_result", autospec=True, side_effect=record), \
                mock.patch.object(PanMetadataPanel, "apply_result", autospec=True, side_effect=record), \
                mock.patch.object(CgmlstForcePanel, "apply_result", autospec=True, side_effect=record):
            for kind in ("geo", "pan_metadata", "cgmlst_force"):
                frame = win._panels.frame(kind)
                sub = frame._result_sub
                payload = sub.latest() if sub.dirty else (win.ctx.applied, win.ctx.result)
                self.assertIs(payload[1], win.ctx.result, kind)
                frame._deliver_result(payload)
        return seen

    def test_every_panel_gets_the_same_sample_set(self) -> None:
        self._search_gene("blaOXA")
        win = self.win
        samples = win.ctx.result.samples
        self.assertEqual(samples, ("S1", "S3", "S4"))
        seen = self._heavy_panel_samples()
        for name in ("geo", "pan_metadata", "cgmlst_force"):
            self.assertEqual(seen[name], samples, name)
        self.assertEqual(seen["map"]["total"], len(samples))
        self.assertEqual(sorted(p["sample"] for p in seen["map"]["points"]), ["S1", "S4"])
        # The tree's rows are the sample set itself (plan 4.4), and the status
        # bar counts it (checklist A: "the status bar count matches").
        self.assertEqual(tuple(self._tree_order()), samples)
        self.assertEqual(tuple(self._tree_samples_with_hits()), samples)
        self.assertEqual(win.status.text(), "Samples: 3 | Annotations after filters: 4")
        # Without an annotation filter the sample set is every allowed sample.
        self._search_gene("")
        self.assertEqual(self._tree_order(), ["S1", "S2", "S3", "S4", "S5", "S6"])
        self.assertEqual(win.status.text(), f"Samples: 6 | Annotations after filters: {win.ctx.result.filtered_ann.height}")
        # A gene with no hit leaves an empty sample set and an empty tree.
        self._search_gene("no-such-gene")
        self.assertEqual(win.ctx.result.samples, ())
        self.assertEqual(self._tree_order(), [])
        self.assertTrue(win.status.text().startswith("Samples: 0 |"))

    def test_d3_switching_the_table_mode_starts_no_compute(self) -> None:
        win = self.win
        self._search_gene("blaOXA")
        result = win.ctx.result
        status = win.status.text()
        computes = win.ctx.scheduler.compute_count
        tree = self._tree_order()
        panels = self._heavy_panel_samples()
        win._open_tree_contig_payload({"sample": "S1", "contig": "3", "contig_label": "Contig 3"}, mode="Contig")
        for mode in ("Context", "Global", "Contig", "Global", "Context"):
            with self.subTest(mode=mode):
                win.annotations_panel.set_view_mode(mode)
                self.assertTrue(win.ctx.scheduler.is_idle())
                self.assertEqual(win.view_mode, mode)
                self.assertIs(win.ctx.result, result)
                self.assertEqual(win.ctx.scheduler.compute_count, computes)
                self.assertEqual(self._tree_order(), tree)
                self.assertEqual(self._heavy_panel_samples(), panels)
                self.assertEqual(win.status.text(), status)
        # The table redrew for each mode; its own line says what it shows.
        win.annotations_panel.set_view_mode("Global")
        self.assertEqual(win.table.model().rowCount(), 4)
        self.assertEqual(win.annotations_panel.status_text(), "4 annotations across 3 samples")
        win.annotations_panel.set_view_mode("Contig")
        self.assertEqual(win.annotations_panel.status_text(), "Sample S1 | Contig 3 | 1 hit, 2 annotation rows")
        # The panel's own API does the same (D3: set_view_mode redraws, no compute).
        win.annotations_panel.set_view_mode("Global")
        self.assertEqual(win.view_mode, "Global")
        self.assertEqual(win.annotations_panel.mode_box().currentData(), "Global")
        self.assertEqual(win.table.model().rowCount(), 4)
        self.assertEqual(win.ctx.scheduler.compute_count, computes)

    def test_open_in_context_view_starts_no_compute(self) -> None:
        win = self.win
        computes = win.ctx.scheduler.compute_count
        seen: List[Any] = []
        win.ctx.selectionChanged.connect(lambda sel, _origin: seen.append(sel))
        try:
            win._open_tree_contig_payload({"sample": "S3", "contig": "2"}, mode="Context")
        finally:
            win.ctx.selectionChanged.disconnect()
        self.assertEqual(win.view_mode, "Context")
        self.assertEqual(win.ctx.selection, Selection(sample="S3", contig="2", contigs=("2",), contig_label="2"))
        self.assertTrue(seen)
        self.assertEqual(win.ctx.scheduler.compute_count, computes)
        self.assertTrue(win.ctx.scheduler.is_idle())
        # No annotation filter: no hit count (every row would count as a hit).
        self.assertEqual(win.annotations_panel.status_text(), "Sample S3 | 2 | 2 rows shown")

    def test_gene_name_in_the_sample_search_filters_samples_and_hints_once(self) -> None:
        win = self.win
        win.annotations_panel.set_view_mode("Global")
        win.sample_query.setText("blaoxa")
        _wait_idle(self.app, win)
        self.assertEqual(win.ctx.result.samples, ())
        self.assertEqual(win.ctx.result.allowed_samples, ())
        self.assertFalse(win.sample_search_hint.isHidden())
        self.assertEqual(
            win.sample_search_hint.text(),
            "The sample search filters samples only. To find genes, use Gene/annotation contains or Product contains.",
        )
        self.assertTrue(QtCore.QSettings().value(f"{HINTS_GROUP}/{HINT_SAMPLE_SEARCH_D3}", False, type=bool))
        # The hint goes with its text and never comes back.
        win.sample_query.setText("blaoxa-23")
        _wait_idle(self.app, win)
        self.assertTrue(win.sample_search_hint.isHidden())
        win.sample_query.setText("ompa")
        _wait_idle(self.app, win)
        self.assertTrue(win.ctx.result.sample_query_matches_annotations)
        self.assertTrue(win.sample_search_hint.isHidden())

    def test_no_hint_when_the_text_names_samples_or_nothing(self) -> None:
        win = self.win
        for text in ("romania", "nothing-matches-this"):
            with self.subTest(text=text):
                win.sample_query.setText(text)
                _wait_idle(self.app, win)
                self.assertTrue(win.sample_search_hint.isHidden())
        self.assertFalse(QtCore.QSettings().value(f"{HINTS_GROUP}/{HINT_SAMPLE_SEARCH_D3}", False, type=bool))

    def test_sorting_is_tree_local(self) -> None:
        win = self.win
        computes = win.ctx.scheduler.compute_count
        result = win.ctx.result
        self.assertEqual(self._tree_order(), ["S1", "S2", "S3", "S4", "S5", "S6"])
        win.sample_sort_order_box.setCurrentIndex(win.sample_sort_order_box.findData("desc"))
        self.assertEqual(self._tree_order(), ["S6", "S5", "S4", "S3", "S2", "S1"])
        win.sample_sort_column_box.setCurrentIndex(win.sample_sort_column_box.findData("Country"))
        # Countries descending; S5 has none and comes last.
        self.assertEqual(self._tree_order(), ["S3", "S4", "S1", "S6", "S2", "S5"])
        self.assertTrue(win.ctx.scheduler.is_idle())
        self.assertEqual(win.ctx.scheduler.compute_count, computes)
        self.assertIs(win.ctx.result, result)

    def test_the_selection_lives_in_the_context_only(self) -> None:
        # Phase 10: the window's current_* views over the selection are gone;
        # panels and the window read and change ctx.selection.
        win = self.win
        for name in ("current_sample", "current_contig", "current_contigs", "current_contig_label"):
            self.assertFalse(hasattr(win, name), name)
        win._open_tree_contig_payload({"sample": "S2", "contig": "1", "contig_label": "Contig 1"}, mode="Contig")
        self.assertEqual(win.ctx.selection.sample, "S2")
        self.assertEqual(win.ctx.selection.contig, "1")

    def test_shared_state_reads_the_applied_result_not_the_widgets(self) -> None:
        win = self.win
        self._search_gene("blaOXA")
        # Typed but not applied yet (still debouncing): everything on screen
        # still belongs to the applied search.
        win.gene_query.setText("ompA")
        self.assertTrue(win.ctx.scheduler.debounce_timer.isActive())
        self.assertEqual(win._applied_filters().gene, "blaOXA")
        self.assertEqual(win.ctx.applied.ann_filters()["gene_text"], "blaOXA")
        self.assertEqual(win._sample_tree_search_query(), "blaOXA")
        self.assertEqual(win.ctx.sample_set(), ("S1", "S3", "S4"))
        self.assertTrue(win._has_active_annotation_filter())
        _wait_idle(self.app, win)
        self.assertEqual(win._applied_filters().gene, "ompA")
        self.assertEqual(win.ctx.sample_set(), ("S1", "S2", "S4", "S5"))

    def test_the_hint_stays_hidden_in_a_new_window_once_shown(self) -> None:
        # The flag lives in QSettings, so a fresh window (no per-window hint
        # state) never shows it again.
        win = self.win
        QtCore.QSettings().setValue(f"{HINTS_GROUP}/{HINT_SAMPLE_SEARCH_D3}", True)
        win._sample_search_hint_query = None
        win.sample_search_hint.setVisible(False)
        win.sample_query.setText("blaoxa")
        _wait_idle(self.app, win)
        self.assertTrue(win.ctx.result.sample_query_matches_annotations)
        self.assertTrue(win.sample_search_hint.isHidden())
        holder = type("SecondWindow", (), {})
        second = _ExplorerFixture.build(holder)
        try:
            _wait_idle(self.app, second)
            second.sample_query.setText("blaoxa")
            _wait_idle(self.app, second)
            self.assertTrue(second.ctx.result.sample_query_matches_annotations)
            self.assertTrue(second.sample_search_hint.isHidden())
        finally:
            _ExplorerFixture.teardown(holder, second)

    def test_a_selection_outside_the_new_sample_set_is_dropped(self) -> None:
        win = self.win
        win._open_tree_contig_payload({"sample": "S2", "contig": "1", "contig_label": "Contig 1"}, mode="Contig")
        self.assertEqual(win.ctx.selection.sample, "S2")
        self.assertGreater(win.table.model().rowCount(), 0)
        self._search_gene("blaOXA")  # S2 has no blaOXA hit
        self.assertNotIn("S2", win.ctx.result.samples)
        self.assertEqual(win.ctx.selection, Selection())
        self.assertEqual(win.table.model().rowCount(), 0)
        self.assertEqual(win.annotations_panel.status_text(), "")
        # A selection inside the sample set stays.
        win._open_tree_contig_payload({"sample": "S1", "contig": "3", "contig_label": "Contig 3"}, mode="Contig")
        self._search_gene("blaOXA-23")
        self.assertIn("S1", win.ctx.result.samples)
        self.assertEqual(win.ctx.selection.sample, "S1")

    def test_the_track_follows_the_applied_annotation_filters_only(self) -> None:
        # D3 step 2: the track's hit arrows and its cache key read the applied
        # annotation filters, never the sample search or the boxes ahead of
        # the result.
        win = self.win
        drawn: List[dict] = []
        real_set_track = win.igv_view.set_track

        def record(*args: Any, **kwargs: Any) -> Any:
            drawn.append(dict(kwargs))
            return real_set_track(*args, **kwargs)

        with mock.patch.object(win.igv_view, "set_track", side_effect=record):
            self._search_gene("blaOXA")
            win._open_tree_contig_payload({"sample": "S1", "contig": "3", "contig_label": "Contig 3"}, mode="Contig")
            key = win._igv_key
            arrows = drawn[-1]["hit_bp_positions"]
            self.assertTrue(arrows)
            # A gene typed but not applied yet: a redraw keeps the applied key.
            win.gene_query.setText("ompA")
            self.assertTrue(win.ctx.scheduler.debounce_timer.isActive())
            win._igv_dirty = True
            win._apply_view_mode()
            self.assertEqual(win._igv_key, key)
            self.assertEqual(drawn[-1]["hit_bp_positions"], arrows)
            win.gene_query.setText("blaOXA")
            _wait_idle(self.app, win)
            # The sample search never makes hit arrows nor keys the track.
            win.sample_query.setText("S1")
            _wait_idle(self.app, win)
            self.assertEqual(win._igv_key, key)
            self._search_gene("")
            self.assertIn("S1", win.ctx.result.samples)
            win._igv_dirty = True
            win._apply_view_mode()
            self.assertNotEqual(win._igv_key, key)
            self.assertEqual(drawn[-1]["hit_bp_positions"], [])

    def test_typing_during_a_slow_compute_ends_on_the_last_search(self) -> None:
        # Checklist A: hold a real window compute, type a second gene, release;
        # everything on screen must belong to the last search.
        win = self.win
        gate = threading.Event()
        started = threading.Event()
        real_compute = win.engine.compute

        def slow_compute(filters: FilterState) -> AppliedResult:
            if filters.gene == "blaOXA":
                started.set()
                gate.wait(10)
            return real_compute(filters)

        with mock.patch.object(win.engine, "compute", side_effect=slow_compute):
            self.addCleanup(gate.set)
            win.gene_query.setText("blaOXA")
            end = time.monotonic() + 5
            while not started.is_set() and time.monotonic() < end:
                self.app.processEvents(QtCore.QEventLoop.ProcessEventsFlag.AllEvents, 20)
            self.assertTrue(started.is_set())
            win.gene_query.setText("ompA")
            for _ in range(30):  # let the debounce fire while the first compute is held
                self.app.processEvents(QtCore.QEventLoop.ProcessEventsFlag.AllEvents, 20)
                time.sleep(0.01)
            gate.set()
            _wait_idle(self.app, win)
        expected = win.engine.compute(FilterState(gene="ompA"))
        self.assertEqual(win._applied_filters().gene, "ompA")
        self.assertEqual(win.ctx.applied.ann_filters()["gene_text"], "ompA")
        self.assertEqual(win.ctx.result.samples, expected.samples)
        self.assertEqual(tuple(self._tree_order()), expected.samples)
        self.assertEqual(win.ctx.sample_set(), expected.samples)
        self.assertEqual(win.status.text(), expected.stats_text)
        self.assertEqual(win._sample_tree_search_query(), "ompA")

    def test_the_hovered_sample_shows_in_the_status_bar(self) -> None:
        win = self.win
        win.ctx.hover(["S4"], origin="test")
        win.ctx.flush_hover()
        self.assertEqual(win.hover_status.text(), "Hovered: S4")
        win.ctx.hover(["S1", "S3"], origin="test", kind="samples")
        win.ctx.flush_hover()
        self.assertEqual(win.hover_status.text(), "Hovered: 2 samples")
        win.ctx.clear_hover("test")
        win.ctx.flush_hover()
        self.assertEqual(win.hover_status.text(), "")

    def test_window_filter_checkboxes_live_in_their_panels(self) -> None:
        # Feedback 2026-10-07: the rail lost "Hide hypothetical" and "Hide
        # contigs without hits"; every Annotations instance shows the first,
        # the Sample tree the second, and both write the window's draft.
        win = self.win
        self.assertFalse(hasattr(win, "hide_hypo_chk"))
        self.assertFalse(hasattr(win, "hide_empty_contigs_chk"))
        rail_texts = {c.text() for c in win._filters_panel.findChildren(QtWidgets.QCheckBox)} if win._filters_panel else set()
        self.assertNotIn("Hide hypothetical", rail_texts)
        self.assertNotIn("Hide contigs without hits", rail_texts)
        paired = win.annotations_panel
        manager = win._panel_manager
        tree_panel = manager.panel("sample_tree")
        self.assertIsNone(paired.findChild(QtWidgets.QCheckBox, "HideEmptyContigsCheck"))
        self.assertIs(tree_panel.findChild(QtWidgets.QCheckBox, "HideEmptyContigsCheck"), tree_panel.hide_empty_contigs_chk)
        self.assertIsNone(tree_panel.findChild(QtWidgets.QCheckBox, "HideHypotheticalCheck"))
        extra_uid = manager.open_new_instance("annotations")
        try:
            _wait_idle(self.app, win)
            extra = manager.panel(extra_uid)
            # A new instance shows the current window state.
            self.assertFalse(extra.hide_hypo_chk.isChecked())
            self.assertIsNone(extra.findChild(QtWidgets.QCheckBox, "HideEmptyContigsCheck"))
            self.assertFalse(tree_panel.hide_empty_contigs_chk.isChecked())
            computes = win.ctx.scheduler.compute_count
            extra.hide_hypo_chk.setChecked(True)
            _wait_idle(self.app, win)
            self.assertTrue(win.ctx.applied.hide_hypothetical)
            self.assertTrue(paired.hide_hypo_chk.isChecked())
            self.assertEqual(win.ctx.scheduler.compute_count, computes + 1)
            tree_panel.hide_empty_contigs_chk.setChecked(True)
            _wait_idle(self.app, win)
            self.assertTrue(win.ctx.applied.hide_empty_contigs)
            self.assertTrue(win.ctx.draft.hide_empty_contigs)
            # The draft change does not echo back into a second toggle.
            self.assertTrue(tree_panel.hide_empty_contigs_chk.isChecked())
            self.assertEqual(win.ctx.scheduler.compute_count, computes + 2)
            # The rail clears the flag (draft edited elsewhere): the
            # checkbox follows without writing the draft again.
            win.ctx.update_draft(win.ctx.draft.replace(hide_empty_contigs=False))
            _wait_idle(self.app, win)
            self.assertFalse(tree_panel.hide_empty_contigs_chk.isChecked())
            self.assertEqual(win.ctx.scheduler.compute_count, computes + 3)
            tree_panel.hide_empty_contigs_chk.setChecked(True)
            _wait_idle(self.app, win)
            # Rail edits and refresh() keep the flags (they come from the draft).
            self._search_gene("blaOXA")
            win.refresh()
            _wait_idle(self.app, win)
            self.assertTrue(win.ctx.applied.hide_hypothetical)
            self.assertTrue(win.ctx.applied.hide_empty_contigs)
            self.assertEqual(win.ctx.applied.gene.lower(), "blaoxa")
        finally:
            manager.close_panel(extra_uid)
            _wait_idle(self.app, win)


if __name__ == "__main__":
    unittest.main()
