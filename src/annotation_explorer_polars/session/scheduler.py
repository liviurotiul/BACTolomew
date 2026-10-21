"""Apply pipeline of one explore window (plan 4.4, Phase 3).

Single flight and latest wins: at most one compute runs at a time; while it
runs, newer requests replace each other and only the last one is computed
next. Every request gets a generation token and a finished compute whose token
is no longer the latest is dropped. Text edits are debounced (180 ms); other
edits apply on the next pass of the event loop, so a burst of them (several
controls set by one action, or a "Clear" that unchecks a list) becomes one
request. The pool thread runs only ``engine.compute``, never window code.
"""
from __future__ import annotations

from typing import Any, Optional, Tuple

from PyQt6 import QtCore

from ..model.filter_state import AppliedResult, FilterState
from ..refresh_worker import RefreshWorker

DEBOUNCE_MS = 180


class ApplyScheduler(QtCore.QObject):
    applied = QtCore.pyqtSignal(object, object)  # FilterState, AppliedResult
    failed = QtCore.pyqtSignal(object, str)  # FilterState, traceback text
    busyChanged = QtCore.pyqtSignal(bool)

    def __init__(
        self,
        engine: Any,
        parent: Optional[QtCore.QObject] = None,
        *,
        pool: Optional[QtCore.QThreadPool] = None,
    ) -> None:
        super().__init__(parent)
        self._engine = engine
        self._pool = pool or QtCore.QThreadPool.globalInstance()
        self._generation = 0
        self._running = False
        self._busy = False
        self._pending: Optional[Tuple[int, FilterState]] = None
        self._running_filters: Optional[FilterState] = None
        self._wanted: Optional[FilterState] = None
        self._force_wanted = False
        self._last_requested: Optional[FilterState] = None
        self._worker: Optional[RefreshWorker] = None
        self.compute_count = 0  # computes started (tests and diagnostics)
        self.debounce_timer = QtCore.QTimer(self)
        self.debounce_timer.setSingleShot(True)
        self.debounce_timer.setInterval(DEBOUNCE_MS)
        self.debounce_timer.timeout.connect(self._flush_debounced)

    # State

    @property
    def generation(self) -> int:
        return self._generation

    def is_running(self) -> bool:
        return self._running

    def pending_filters(self) -> Optional[FilterState]:
        return self._pending[1] if self._pending is not None else None

    def is_idle(self) -> bool:
        return not (self._running or self._pending is not None or self.debounce_timer.isActive())

    # Requests

    def request(self, filters: FilterState, *, debounce: bool = False, force: bool = False) -> None:
        """Ask for ``filters`` to be applied: after the debounce delay, or on
        the next pass of the event loop. The latest filters asked for win.

        Without ``force`` a request equal to the last one sent is ignored (for
        example a text edited and edited back within the debounce delay).
        """
        self._wanted = filters
        self._force_wanted = self._force_wanted or force
        timer = self.debounce_timer
        if debounce and not (timer.isActive() and timer.interval() == 0):
            timer.setInterval(DEBOUNCE_MS)
        else:
            # Not debounced, or a prompt flush is already due: it takes the
            # latest filters (typed text included).
            timer.setInterval(0)
        timer.start()

    def _flush_debounced(self) -> None:
        force = self._force_wanted
        self._force_wanted = False
        if self._wanted is not None:
            self._submit(self._wanted, force=force)

    def _submit(self, filters: FilterState, *, force: bool) -> None:
        if not force and filters == self._last_requested:
            return
        self._last_requested = filters
        self._generation += 1
        self._pending = (self._generation, filters)
        if not self._running:
            self._start_next()

    def _start_next(self) -> None:
        pending = self._pending
        if pending is None:
            self._set_busy(self._running)
            return
        self._pending = None
        token, filters = pending
        self._running = True
        self._running_filters = filters
        self.compute_count += 1
        self._set_busy(True)
        engine = self._engine
        worker = RefreshWorker(token, lambda t=token, f=filters, e=engine: (t, e.compute(f)))
        worker.signals.finished.connect(self._on_finished)
        worker.signals.failed.connect(self._on_failed)
        self._worker = worker
        self._pool.start(worker)

    def _on_finished(self, payload: object) -> None:
        self._running = False
        self._worker = None
        filters = self._running_filters
        self._running_filters = None
        token, result = payload if isinstance(payload, tuple) else (None, None)
        if token == self._generation and isinstance(result, AppliedResult) and filters is not None:
            self.applied.emit(filters, result)
        self._start_next()

    def _on_failed(self, token: int, error_text: str) -> None:
        self._running = False
        self._worker = None
        filters = self._running_filters
        self._running_filters = None
        # A failed request may be asked for again.
        if filters == self._last_requested:
            self._last_requested = None
        if token == self._generation:
            self.failed.emit(filters, error_text)
        self._start_next()

    def _set_busy(self, busy: bool) -> None:
        busy = bool(busy)
        if busy != self._busy:
            self._busy = busy
            self.busyChanged.emit(busy)
