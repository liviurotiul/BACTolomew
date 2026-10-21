"""ExploreContext: the shared session of one explore window (plan 4.4).

One per explore window, on the GUI thread. It is the only owner of the
window's linked state: the draft filters the rail shows while the user types,
the applied filters and their result (with the window's one sample set,
``result.samples``), the shared selection and the hover. Widgets display this
state; nobody reads the rail widgets to learn what is on screen. Contexts of
different windows never talk to each other.

Panels learn about changes through the context's ``bus`` (plan 4.5): results,
hover, selection and busy state, each delivered by the panel's frame with its
own delivery class. Panels change shared state only through the mutation API
below (``update_draft``, ``set_subset``, ``select``, ``hover``,
``clear_hover``, ``request``).

The annotation table's view mode is not here (D3): it is a display option of
each Annotations panel and never starts a compute.
"""
from __future__ import annotations

from typing import Any, Iterable, Optional

from PyQt6 import QtCore

from ..model.filter_state import AppliedResult, FilterState, HoverState, SampleSubset, Selection
from . import bus as bus_events
from .bus import PanelBus
from .scheduler import ApplyScheduler

__all__ = [
    "AppliedResult",
    "ExploreContext",
    "FilterState",
    "HoverState",
    "SampleSubset",
    "Selection",
]

# Hover is coalesced at the context (plan 4.9). With the heavy panels' own
# coalescing (session/bus.py) a hover reaches every visible heavy panel within
# the 30 ms of plan 8.12 (6 + 16 ms plus the handlers; Phase 10 measurement).
HOVER_COALESCE_MS = 6


class ExploreContext(QtCore.QObject):
    # Every panel redraws from (filters, result); emitted once per applied result.
    resultApplied = QtCore.pyqtSignal(object, object)
    # The current request failed (traceback text).
    applyFailed = QtCore.pyqtSignal(str)
    # Once, after the first result (or the first failure): dismisses the splash.
    firstResultApplied = QtCore.pyqtSignal()
    busyChanged = QtCore.pyqtSignal(bool)
    draftChanged = QtCore.pyqtSignal(object)
    selectionChanged = QtCore.pyqtSignal(object, str)  # Selection, origin
    hoverChanged = QtCore.pyqtSignal(object)  # HoverState
    subsetChanged = QtCore.pyqtSignal(object, str)  # SampleSubset or None, origin
    openRequested = QtCore.pyqtSignal(object, str)  # an open request, origin

    def __init__(
        self,
        engine: Any,
        parent: Optional[QtCore.QObject] = None,
        *,
        pool: Optional[QtCore.QThreadPool] = None,
    ) -> None:
        super().__init__(parent)
        self.engine = engine
        self.draft = FilterState()
        self.applied = FilterState()
        self.result: Optional[AppliedResult] = None
        self.selection = Selection()
        self.hover_state = HoverState()
        self._pending_hover: Optional[HoverState] = None
        self._first_done = False
        self.bus = PanelBus(self)
        self.scheduler = ApplyScheduler(engine, self, pool=pool)
        self.scheduler.applied.connect(self._on_applied)
        self.scheduler.failed.connect(self._on_failed)
        self.scheduler.busyChanged.connect(self._on_busy_changed)
        self._hover_timer = QtCore.QTimer(self)
        self._hover_timer.setSingleShot(True)
        self._hover_timer.setInterval(HOVER_COALESCE_MS)
        self._hover_timer.timeout.connect(self.flush_hover)

    # Reading

    @property
    def first_result_done(self) -> bool:
        return self._first_done

    def is_busy(self) -> bool:
        return not self.scheduler.is_idle()

    def sample_set(self) -> tuple:
        """The window's one sample set (empty before the first result)."""
        return self.result.samples if self.result is not None else ()

    # Mutation API (plan 4.4)

    def update_draft(
        self,
        draft: FilterState,
        *,
        debounce: bool = False,
        force: bool = False,
        origin: str = "rail",
    ) -> None:
        """The rail changed: text edits apply after the debounce delay, the
        other controls at once (on the next pass of the event loop, so one
        action that sets several controls makes one compute; it also takes any
        typed text still waiting in the debounce). ``force`` computes even if
        the filters did not change."""
        del origin  # the scheduler does not care who edited the draft
        if draft != self.draft:
            self.draft = draft
            self.draftChanged.emit(draft)
        self.scheduler.request(draft, debounce=debounce, force=force)

    def set_subset(self, subset: Optional[SampleSubset], origin: str = "") -> None:
        """Replace the one sample-subset filter (the chip in the rail), or clear
        it with None, and apply at once. This is the only way a panel narrows
        the window's sample set, and it narrows it for every panel."""
        if subset is not None and not subset.samples:
            subset = None
        if subset == self.draft.subset:
            return
        self.update_draft(self.draft.replace(subset=subset), origin=origin)
        self.subsetChanged.emit(subset, str(origin or ""))

    def apply(self, *, force: bool = True) -> None:
        """Apply the current draft at once (forced: computes even if nothing changed)."""
        self.scheduler.request(self.draft, force=force)

    def select(self, selection: Selection, origin: str = "") -> None:
        if selection == self.selection:
            return
        self.selection = selection
        self.selectionChanged.emit(selection, origin)
        self.bus.publish(bus_events.SELECTION_CHANGED, selection, origin or None)

    def hover(self, samples: Iterable[str], origin: Optional[str] = None, kind: str = "sample") -> None:
        """One or several samples are hovered in ``origin``. Coalesced: only the
        latest hover inside HOVER_COALESCE_MS is published; publishing the samples already
        shown does nothing."""
        self._pending_hover = HoverState(tuple(samples or ()), origin=origin, kind=kind)
        if not self._hover_timer.isActive():
            self._hover_timer.start()

    def clear_hover(self, origin: Optional[str] = None) -> None:
        self.hover((), origin=origin)

    def flush_hover(self) -> None:
        """Publish the pending hover now (the coalescing timer calls this)."""
        self._hover_timer.stop()
        state, self._pending_hover = self._pending_hover, None
        if state is None:
            return
        current = self.hover_state
        # Same samples from the same producer: nothing to redraw. A change of
        # origin is published even for the same samples, because the new
        # producer's previous state was never delivered to the old one (the
        # bus does not echo), e.g. tree row X -> table row of X.
        if state == current:
            return
        self.hover_state = state
        self.hoverChanged.emit(state)
        self.bus.publish(bus_events.HOVER_CHANGED, state, state.origin)

    def request(self, open_request: Any, origin: str = "") -> None:
        """A typed open request (session/requests.py); the window carries it out."""
        self.openRequested.emit(open_request, str(origin or ""))

    # Scheduler results

    def _on_applied(self, filters: FilterState, result: AppliedResult) -> None:
        self.applied = filters
        self.result = result
        # A selection outside the new sample set is gone from the tree; drop
        # it before anyone draws, so the table and track do not keep it.
        dropped = False
        selected = str(self.selection.sample or "").strip()
        if selected and selected not in set(result.samples):
            self.selection = Selection()
            dropped = True
        self.resultApplied.emit(filters, result)
        self.bus.publish(bus_events.RESULT_APPLIED, (filters, result), None)
        if dropped:
            self.selectionChanged.emit(self.selection, "apply")
            self.bus.publish(bus_events.SELECTION_CHANGED, self.selection, None)
        self._mark_first_done()

    def _on_failed(self, _filters: object, error_text: str) -> None:
        self.applyFailed.emit(error_text)
        self._mark_first_done()

    def _on_busy_changed(self, busy: bool) -> None:
        self.busyChanged.emit(bool(busy))
        self.bus.publish(bus_events.BUSY_CHANGED, bool(busy), None)

    def _mark_first_done(self) -> None:
        if not self._first_done:
            self._first_done = True
            self.firstResultApplied.emit()
