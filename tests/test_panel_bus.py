"""The event bus, the context's hover and the panel frames' delivery (window plan, Phase 4).

Plan 8.5: no echo; hidden panels go dirty and replay on show; unbind stops
delivery; multi-sample hover; a raising handler is isolated; no leaks after
dispose (checked with weakref). Offscreen, no window: fake panels in plain
frames on a real ``ExploreContext`` (its engine is never asked to compute;
results are handed to the context as the scheduler would).
"""
from __future__ import annotations

import gc
import os
import time
import unittest
import weakref
from typing import Any, List
from unittest import mock

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

import polars as pl
from PyQt6 import QtCore, QtWidgets

from annotation_explorer_polars.model.filter_state import (
    AppliedResult,
    FilterState,
    HoverState,
    SampleSubset,
    Selection,
)
from annotation_explorer_polars.panels.base import ExplorePanel, PanelFrame, PanelRegistry, PanelSpec
from annotation_explorer_polars.session import bus as bus_events
from annotation_explorer_polars.session.bus import Delivery, PanelBus
from annotation_explorer_polars.session.explore_context import ExploreContext
from annotation_explorer_polars.session.requests import FocusPanel, OpenSampleCard, params_tuple


def _app() -> QtWidgets.QApplication:
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


def _spin(seconds: float = 0.05) -> None:
    app = _app()
    end = time.monotonic() + seconds
    while True:
        app.processEvents(QtCore.QEventLoop.ProcessEventsFlag.AllEvents, 10)
        if time.monotonic() >= end:
            break
        time.sleep(0.003)


def _result(*samples: str) -> AppliedResult:
    return AppliedResult(
        samples=samples,
        allowed_samples=samples,
        filtered_ann=pl.DataFrame(),
        metadata=pl.DataFrame(),
        contig_records={},
        filtered_counts={},
        stats_text=f"Samples: {len(samples)}",
    )


class FakePanel(ExplorePanel):
    KIND = "fake"
    TITLE = "Fake"

    def __init__(self) -> None:
        super().__init__()
        self.results: List[Any] = []
        self.hovers: List[Any] = []
        self.selections: List[Any] = []
        self.shown = 0
        self.hidden = 0
        self.disposed = False
        self.unbound = False
        self.state = {"option": "a"}

    def apply_result(self, filters: Any, result: Any) -> None:
        self.results.append(result)

    def show_hover(self, hover: Any) -> None:
        self.hovers.append(hover)

    def show_selection(self, selection: Any) -> None:
        self.selections.append(selection)

    def on_shown(self) -> None:
        self.shown += 1

    def on_hidden(self) -> None:
        self.hidden += 1

    def unbind(self) -> None:
        self.unbound = True
        super().unbind()

    def dispose(self) -> None:
        self.disposed = True

    def save_state(self) -> dict:
        return dict(self.state)

    def restore_state(self, state: dict) -> None:
        self.state = dict(state)


def _spec(kind: str = "fake", *, heavy: bool = False, raising: bool = False) -> PanelSpec:
    def factory(_uid: str, _parent: QtWidgets.QWidget) -> ExplorePanel:
        panel = FakePanel()
        if raising:
            def boom(*_args: Any) -> None:
                raise RuntimeError("broken panel")

            panel.apply_result = boom  # type: ignore[assignment]
            panel.show_hover = boom  # type: ignore[assignment]
        return panel

    return PanelSpec(kind, kind.title(), "Test", heavy=heavy, factory=factory)


class PanelBusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = _app()

    def setUp(self) -> None:
        self.bus = PanelBus()
        self.addCleanup(self.bus.deleteLater)

    def test_no_echo_to_the_origin(self) -> None:
        got: dict = {"a": [], "b": []}
        self.bus.subscribe("a", bus_events.HOVER_CHANGED, got["a"].append)
        self.bus.subscribe("b", bus_events.HOVER_CHANGED, got["b"].append)
        self.bus.publish(bus_events.HOVER_CHANGED, "x", origin="a")
        self.assertEqual(got, {"a": [], "b": ["x"]})
        self.bus.publish(bus_events.HOVER_CHANGED, "y", origin=None)
        self.assertEqual(got, {"a": ["y"], "b": ["x", "y"]})

    def test_publishing_while_dispatching_is_queued_not_recursive(self) -> None:
        order: List[str] = []

        def first(payload: str) -> None:
            order.append(f"first:{payload}")
            if payload == "one":
                self.bus.publish(bus_events.SELECTION_CHANGED, "two")
                order.append("first:after-publish")

        self.bus.subscribe("p1", bus_events.SELECTION_CHANGED, first)
        self.bus.subscribe("p2", bus_events.SELECTION_CHANGED, lambda p: order.append(f"second:{p}"))
        self.bus.publish(bus_events.SELECTION_CHANGED, "one")
        self.assertEqual(
            order,
            ["first:one", "first:after-publish", "second:one", "first:two", "second:two"],
        )

    def test_a_raising_handler_is_isolated(self) -> None:
        got: List[str] = []
        failures: List[tuple] = []
        self.bus.handlerFailed.connect(lambda owner, event, text: failures.append((owner, event, text)))

        def broken(_payload: Any) -> None:
            raise ValueError("boom")

        self.bus.subscribe("broken", bus_events.RESULT_APPLIED, broken)
        self.bus.subscribe("fine", bus_events.RESULT_APPLIED, got.append)
        with self.assertLogs("annotation_explorer_polars.session.bus", level="ERROR") as logs:
            self.bus.publish(bus_events.RESULT_APPLIED, "r")
        self.assertEqual(got, ["r"])
        self.assertEqual(failures, [("broken", bus_events.RESULT_APPLIED, "ValueError: boom")])
        self.assertIn("panel broken", logs.output[0])

    def test_unsubscribe_owner_stops_delivery(self) -> None:
        got: List[str] = []
        self.bus.subscribe("p", bus_events.HOVER_CHANGED, got.append)
        self.bus.subscribe("p", bus_events.RESULT_APPLIED, got.append)
        self.bus.unsubscribe_owner("p")
        self.bus.publish(bus_events.HOVER_CHANGED, "h")
        self.bus.publish(bus_events.RESULT_APPLIED, "r")
        self.assertEqual(got, [])
        self.assertEqual(self.bus.subscriptions("p"), [])

    def test_coalesced_delivers_only_the_latest(self) -> None:
        got: List[str] = []
        self.bus.subscribe("p", bus_events.HOVER_CHANGED, got.append, delivery=Delivery.COALESCED, coalesce_ms=30)
        for value in ("a", "b", "c"):
            self.bus.publish(bus_events.HOVER_CHANGED, value)
        self.assertEqual(got, [])
        _spin(0.1)
        self.assertEqual(got, ["c"])

    def test_when_visible_waits_for_the_owner_to_be_shown(self) -> None:
        got: List[str] = []
        visible = {"v": False}
        sub = self.bus.subscribe(
            "p",
            bus_events.RESULT_APPLIED,
            got.append,
            delivery=Delivery.WHEN_VISIBLE,
            delay_ms=40,
            visible=lambda: visible["v"],
        )
        self.bus.publish(bus_events.RESULT_APPLIED, "r1")
        self.bus.publish(bus_events.RESULT_APPLIED, "r2")
        _spin(0.1)
        self.assertEqual(got, [])
        self.assertTrue(sub.dirty)
        visible["v"] = True
        self.bus.notify_shown("p")
        _spin(0.05)
        self.assertEqual(got, ["r2"])  # only the latest
        self.assertFalse(sub.dirty)
        # Visible: delivered after the delay, not at once.
        self.bus.publish(bus_events.RESULT_APPLIED, "r3")
        self.assertEqual(got, ["r2"])
        self.assertTrue(self.bus.has_pending())
        _spin(0.1)
        self.assertEqual(got, ["r2", "r3"])


class ContextHoverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = _app()

    def setUp(self) -> None:
        self.ctx = ExploreContext(engine=None)
        self.addCleanup(self.ctx.deleteLater)
        self.seen: List[HoverState] = []
        self.ctx.bus.subscribe("listener", bus_events.HOVER_CHANGED, self.seen.append)

    def test_hover_is_coalesced_and_carries_every_sample(self) -> None:
        self.ctx.hover(["S1"], origin="tree")
        self.ctx.hover(["S2", "S3", "S2"], origin="phylogeny", kind="samples")
        self.assertEqual(self.seen, [])
        _spin(0.06)
        self.assertEqual(len(self.seen), 1)
        state = self.seen[0]
        self.assertEqual(state.samples, ("S2", "S3"))
        self.assertEqual(state.origin, "phylogeny")
        self.assertIsNone(state.sample)
        self.assertEqual(self.ctx.hover_state, state)

    def test_hover_is_idempotent(self) -> None:
        self.ctx.hover(["S1"], origin="tree")
        self.ctx.flush_hover()
        self.ctx.hover(["S1"], origin="tree")
        self.ctx.flush_hover()
        self.assertEqual([s.samples for s in self.seen], [("S1",)])
        self.ctx.clear_hover("tree")
        self.ctx.flush_hover()
        self.assertEqual([s.samples for s in self.seen], [("S1",), ()])
        self.assertEqual(self.seen[-1].kind, "none")

    def test_the_same_samples_from_another_panel_reach_the_first_one(self) -> None:
        # Tree row X, then (within one coalescing window) a table row of X:
        # the tree never received its own hover, so it must get the table's.
        tree_seen: List[HoverState] = []
        self.ctx.bus.subscribe("sample_tree", bus_events.HOVER_CHANGED, tree_seen.append)
        self.ctx.hover(["X"], origin="sample_tree")
        self.ctx.flush_hover()
        self.assertEqual(tree_seen, [])
        self.ctx.hover([], origin="sample_tree")
        self.ctx.hover(["X"], origin="annotations")
        self.ctx.flush_hover()
        self.assertEqual([(h.samples, h.origin) for h in tree_seen], [(("X",), "annotations")])
        self.assertEqual(self.ctx.hover_state.origin, "annotations")
        # The same samples from the same producer again: nothing.
        self.ctx.hover(["X"], origin="annotations")
        self.ctx.flush_hover()
        self.assertEqual(len(tree_seen), 1)

    def test_no_cross_window_delivery(self) -> None:
        # Two windows, two contexts: each bus delivers only its own events.
        other = ExploreContext(engine=None)
        self.addCleanup(other.deleteLater)
        got: List[tuple] = []
        for name in (bus_events.HOVER_CHANGED, bus_events.SELECTION_CHANGED, bus_events.RESULT_APPLIED, bus_events.BUSY_CHANGED):
            other.bus.subscribe("other", name, lambda payload, _name=name: got.append((_name, payload)))
        self.ctx.hover(["S1"], origin="tree")
        self.ctx.flush_hover()
        self.ctx.select(Selection("S1", "1", ("1",), "1"), origin="tree")
        self.ctx._on_applied(FilterState(), _result("S1"))
        self.ctx._on_busy_changed(True)
        self.ctx._on_busy_changed(False)
        _spin(0.4)
        self.assertEqual(got, [])
        self.assertIsNot(other.bus, self.ctx.bus)
        self.assertEqual(other.hover_state, HoverState())
        self.assertEqual(other.selection, Selection())
        self.assertEqual([h.samples for h in self.seen], [("S1",)])

    def test_set_subset_edits_the_draft_and_tells_the_rail(self) -> None:
        changes: List[tuple] = []
        self.ctx.subsetChanged.connect(lambda subset, origin: changes.append((subset, origin)))
        with mock.patch.object(self.ctx.scheduler, "request") as request:
            subset = SampleSubset(("S1", "S2"), label="Clade 3", kind="clade")
            self.ctx.set_subset(subset, origin="phylogeny")
            self.ctx.set_subset(subset, origin="phylogeny")  # unchanged: nothing
            self.ctx.set_subset(None, origin="rail")
        self.assertEqual(changes, [(subset, "phylogeny"), (None, "rail")])
        self.assertIsNone(self.ctx.draft.subset)
        self.assertEqual(request.call_count, 2)
        self.assertEqual(request.call_args_list[0].args[0].subset, subset)

    def test_requests_go_to_the_window(self) -> None:
        got: List[tuple] = []
        self.ctx.openRequested.connect(lambda request, origin: got.append((request, origin)))
        self.ctx.request(OpenSampleCard("S1"), origin="phylogeny")
        self.ctx.request(FocusPanel("annotations", params_tuple({"view_mode": "Context"})), origin="sample_tree")
        self.assertEqual(got[0], (OpenSampleCard("S1"), "phylogeny"))
        self.assertEqual(got[1][0].params_dict(), {"view_mode": "Context"})

    def test_a_selection_outside_the_new_sample_set_is_dropped_before_panels_draw(self) -> None:
        seen_by_panel: List[Selection] = []
        self.ctx.bus.subscribe(
            "panel", bus_events.RESULT_APPLIED, lambda _payload: seen_by_panel.append(self.ctx.selection)
        )
        self.ctx.select(Selection("S9", "1", ("1",), "1"), origin="tree")
        self.ctx._on_applied(FilterState(), _result("S1"))
        self.assertEqual(seen_by_panel, [Selection()])
        self.assertEqual(self.ctx.selection, Selection())


class PanelFrameDeliveryTests(unittest.TestCase):
    """Frames gate delivery by visibility (plan 3.9)."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.app = _app()

    def setUp(self) -> None:
        self.ctx = ExploreContext(engine=None)
        self.addCleanup(self.ctx.deleteLater)
        self.host = QtWidgets.QWidget()
        self.host.resize(400, 300)
        QtWidgets.QVBoxLayout(self.host)
        self.addCleanup(self.host.deleteLater)
        self.addCleanup(self.host.hide)

    def _frame(self, uid: str = "fake", **kwargs: Any) -> PanelFrame:
        frame = PanelFrame(_spec(**kwargs), uid=uid)
        self.host.layout().addWidget(frame)
        frame.bind(self.ctx, ready=lambda: True)
        return frame

    def _apply(self, *samples: str) -> AppliedResult:
        result = _result(*samples)
        self.ctx._on_applied(FilterState(), result)
        return result

    def test_light_panels_draw_at_once_even_when_hidden(self) -> None:
        frame = self._frame()
        panel = frame.panel()
        self.assertIsInstance(panel, FakePanel)  # light: built when bound
        self.assertEqual(panel.uid, "fake")
        self.assertIs(panel.ctx, self.ctx)
        result = self._apply("S1")
        self.assertEqual(panel.results, [result])
        self.assertFalse(frame.is_dirty())

    def test_hidden_heavy_panels_go_dirty_and_replay_on_show(self) -> None:
        frame = self._frame(heavy=True)
        self.assertIsNone(frame.panel())  # heavy: built on first show
        self._apply("S1")
        latest = self._apply("S1", "S2")
        _spin(0.3)
        self.assertTrue(frame.is_dirty())
        self.assertIsNone(frame.panel())
        self.host.show()
        _spin(0.1)
        panel = frame.panel()
        self.assertIsInstance(panel, FakePanel)
        self.assertEqual(panel.results, [latest])  # the latest only, once
        self.assertFalse(frame.is_dirty())
        self.assertEqual(panel.shown, 1)
        # Visible: a new result arrives after the delay of the deferred refresh.
        newer = self._apply("S3")
        _spin(0.05)
        self.assertEqual(panel.results, [latest])
        _spin(bus_events.WHEN_VISIBLE_DELAY_MS / 1000 + 0.15)
        self.assertEqual(panel.results, [latest, newer])
        # Hidden again: dirty, and shown again: caught up.
        self.host.hide()
        self.assertEqual(panel.hidden, 1)
        newest = self._apply("S4")
        _spin(0.4)
        self.assertEqual(panel.results[-1], newer)
        self.assertTrue(frame.is_dirty())
        self.host.show()
        _spin(0.1)
        self.assertEqual(panel.results[-1], newest)

    def test_hover_is_dropped_while_hidden_and_replayed_on_show(self) -> None:
        frame = self._frame("light")
        heavy = self._frame("heavy", heavy=True)
        self.host.show()
        _spin(0.1)
        light_panel, heavy_panel = frame.panel(), heavy.panel()
        # On show each panel got the current hover (none) to clear stale marks.
        self.assertEqual([h.samples for h in light_panel.hovers], [()])
        light_panel.hovers.clear()
        heavy_panel.hovers.clear()
        self.ctx.hover(["S1", "S2"], origin="elsewhere", kind="samples")
        self.ctx.flush_hover()
        self.assertEqual([h.samples for h in light_panel.hovers], [("S1", "S2")])
        _spin(0.1)  # heavy panels get hover coalesced
        self.assertEqual([h.samples for h in heavy_panel.hovers], [("S1", "S2")])
        self.host.hide()
        self.ctx.hover(["S3"], origin="elsewhere")
        self.ctx.flush_hover()
        _spin(0.1)
        self.assertEqual([h.samples for h in light_panel.hovers], [("S1", "S2")])
        self.host.show()
        _spin(0.1)
        self.assertEqual(light_panel.hovers[-1].samples, ("S3",))
        self.assertEqual(heavy_panel.hovers[-1].samples, ("S3",))

    def test_no_echo_between_panels(self) -> None:
        first = self._frame("first")
        second = self._frame("second")
        self.host.show()
        _spin(0.05)
        first.panel().hovers.clear()
        second.panel().hovers.clear()
        self.ctx.hover(["S1"], origin="first")
        self.ctx.flush_hover()
        self.assertEqual(first.panel().hovers, [])
        self.assertEqual([h.samples for h in second.panel().hovers], [("S1",)])

    def test_selection_is_delivered_at_once_and_not_echoed(self) -> None:
        first = self._frame("first")
        second = self._frame("second", heavy=True)
        self.host.show()
        _spin(0.1)
        selection = Selection("S1", "1", ("1",), "1")
        self.ctx.select(selection, origin="first")
        # Immediate (plan 4.5): no event loop turn needed, heavy or not.
        self.assertEqual(second.panel().selections, [selection])
        self.assertEqual(first.panel().selections, [])

    def test_a_raising_panel_does_not_stop_the_others(self) -> None:
        broken = self._frame("broken", raising=True)
        fine = self._frame("fine")
        self.assertIsNotNone(broken.panel())
        with self.assertLogs("annotation_explorer_polars.session.bus", level="ERROR"):
            result = self._apply("S1")
        self.assertEqual(fine.panel().results, [result])

    def test_unbind_stops_delivery(self) -> None:
        frame = self._frame()
        panel = frame.panel()
        frame.unbind()
        self.assertTrue(panel.unbound)
        self._apply("S1")
        self.ctx.select(Selection("S1", "1", ("1",), "1"))
        self.ctx.hover(["S1"])
        self.ctx.flush_hover()
        self.assertEqual((panel.results, panel.selections, panel.hovers), ([], [], []))
        self.assertEqual(self.ctx.bus.subscriptions("fake"), [])

    def test_a_panel_opened_later_catches_up(self) -> None:
        result = self._apply("S1", "S2")
        light = self._frame("light")
        self.assertEqual(light.panel().results, [result])
        heavy = self._frame("heavy", heavy=True)
        self.assertTrue(heavy.is_dirty())
        self.host.show()
        _spin(0.1)
        self.assertEqual(heavy.panel().results, [result])

    def test_state_waits_for_a_lazy_build(self) -> None:
        frame = self._frame(heavy=True)
        frame.restore_state({"option": "b"})
        self.assertEqual(frame.save_state(), {"option": "b"})
        self.host.show()
        _spin(0.1)
        self.assertEqual(frame.panel().state, {"option": "b"})
        self.assertEqual(frame.save_state(), {"option": "b"})

    def test_lazy_builds_wait_until_the_window_is_ready(self) -> None:
        ready = {"v": False}
        frame = PanelFrame(_spec(heavy=True), uid="fake")
        self.host.layout().addWidget(frame)
        frame.bind(self.ctx, ready=lambda: ready["v"])
        self._apply("S1")
        self.host.show()
        _spin(0.4)
        self.assertIsNone(frame.panel())
        self.assertTrue(frame.is_dirty())
        ready["v"] = True
        frame.settle_shown()
        _spin(0.05)
        self.assertEqual(len(frame.panel().results), 1)

    def test_no_leaks_after_dispose(self) -> None:
        frame = self._frame("gone", heavy=True)
        self.host.show()
        _spin(0.1)
        panel = frame.panel()
        self.assertIsNotNone(panel)
        ref_panel = weakref.ref(panel)
        ref_frame = weakref.ref(frame)
        frame.dispose()
        self.assertTrue(panel.disposed and panel.unbound)
        self.assertEqual(self.ctx.bus.subscriptions("gone"), [])
        frame.setParent(None)
        frame.deleteLater()
        del frame, panel
        _spin(0.05)
        QtCore.QCoreApplication.sendPostedEvents(None, QtCore.QEvent.Type.DeferredDelete.value)
        gc.collect()
        self.assertIsNone(ref_panel())
        self.assertIsNone(ref_frame())
        # The bus keeps working for the rest.
        self._apply("S1")


class PanelRegistryTests(unittest.TestCase):
    def test_uids_and_kinds(self) -> None:
        registry = PanelRegistry([_spec("geo"), _spec("track")])
        self.assertEqual(registry.kinds(), ["geo", "track"])
        self.assertEqual(PanelRegistry.new_uid("geo", []), "geo")
        uid = PanelRegistry.new_uid("geo", ["geo"])
        self.assertRegex(uid, r"^geo-[0-9a-f]{4}$")
        self.assertEqual(PanelRegistry.kind_of(uid), "geo")
        with self.assertRaises(ValueError):
            PanelRegistry([_spec("geo"), _spec("geo")])


if __name__ == "__main__":
    unittest.main()
