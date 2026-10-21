"""The event bus of one explore window (plan 4.5, Phase 4).

Each ``ExploreContext`` owns one ``PanelBus``. There are no cross-window
events: a bus only ever delivers to the subscribers of its own window.

Events (the ``payload`` each one carries):

- ``RESULT_APPLIED``: ``(filters, result)``, once per applied result;
- ``HOVER_CHANGED``: a ``HoverState`` (one or several samples, or none);
- ``SELECTION_CHANGED``: a ``Selection``;
- ``BUSY_CHANGED``: ``bool``, a compute is running or waiting.

Delivery classes:

- ``IMMEDIATE``: the handler runs now;
- ``COALESCED``: only the latest payload inside ``coalesce_ms`` is delivered;
- ``WHEN_VISIBLE``: the latest payload is stored; it is delivered after
  ``delay_ms`` while the subscriber is visible, otherwise when it becomes
  visible again (``notify_shown``).

Rules: an event is never delivered back to its origin (no echo); publishing
while a dispatch runs is queued, never recursive; every handler is wrapped so
that one broken panel cannot break the fan-out (the error is logged with the
subscriber's uid); subscriptions belong to their owner (a panel frame or the
window's chrome) and ``unsubscribe_owner`` drops all of them at once.
"""
from __future__ import annotations

import enum
import logging
from collections import deque
from typing import Any, Callable, Deque, Dict, List, Optional, Tuple

from PyQt6 import QtCore

_LOG = logging.getLogger(__name__)

RESULT_APPLIED = "result_applied"
HOVER_CHANGED = "hover_changed"
SELECTION_CHANGED = "selection_changed"
BUSY_CHANGED = "busy_changed"
EVENTS = (RESULT_APPLIED, HOVER_CHANGED, SELECTION_CHANGED, BUSY_CHANGED)

# Today's deferred panel refresh waited 250 ms after a result (plan 4.5).
WHEN_VISIBLE_DELAY_MS = 250
# Hover to heavy panels (plan 4.9): at most one redraw per 16 ms per panel
# during a mouse sweep; with the context's 6 ms the whole path stays under the
# 30 ms of plan 8.12 (was 30 ms here and 16 ms in the context until Phase 10).
HEAVY_HOVER_COALESCE_MS = 16


class Delivery(enum.Enum):
    IMMEDIATE = "immediate"
    COALESCED = "coalesced"
    WHEN_VISIBLE = "when_visible"


_NOTHING = object()


class Subscription:
    """One handler of one owner for one event, with its delivery state."""

    def __init__(
        self,
        bus: "PanelBus",
        owner: str,
        event: str,
        handler: Callable[[Any], None],
        delivery: Delivery,
        *,
        coalesce_ms: int = 0,
        delay_ms: int = 0,
        visible: Optional[Callable[[], bool]] = None,
        drop_when_hidden: bool = False,
    ) -> None:
        self.bus = bus
        self.owner = owner
        self.event = event
        self.handler = handler
        self.delivery = delivery
        self.coalesce_ms = max(0, int(coalesce_ms))
        self.delay_ms = max(0, int(delay_ms))
        self.visible = visible
        self.drop_when_hidden = bool(drop_when_hidden)
        self.active = True
        self._latest: Any = _NOTHING
        self._timer: Optional[QtCore.QTimer] = None

    # State (tests and diagnostics)

    @property
    def dirty(self) -> bool:
        """A payload is stored and not delivered yet."""
        return self._latest is not _NOTHING

    def is_pending(self) -> bool:
        """A timed delivery is scheduled."""
        return self._timer is not None and self._timer.isActive()

    def latest(self) -> Any:
        return None if self._latest is _NOTHING else self._latest

    # Delivery

    def _is_visible(self) -> bool:
        if self.visible is None:
            return True
        try:
            return bool(self.visible())
        except RuntimeError:  # the owner's widget is gone
            return False

    def _ensure_timer(self) -> QtCore.QTimer:
        if self._timer is None:
            timer = QtCore.QTimer(self.bus)
            timer.setSingleShot(True)
            timer.timeout.connect(self._on_timer)
            self._timer = timer
        return self._timer

    def offer(self, payload: Any) -> None:
        if not self.active:
            return
        if self.drop_when_hidden and not self._is_visible():
            # Hidden subscribers drop this event (hover); the owner replays
            # the current state when it is shown again.
            self._latest = _NOTHING
            if self._timer is not None:
                self._timer.stop()
            return
        if self.delivery is Delivery.IMMEDIATE:
            self._call(payload)
            return
        self._latest = payload
        if self.delivery is Delivery.COALESCED:
            timer = self._ensure_timer()
            if not timer.isActive():
                timer.start(self.coalesce_ms)
            return
        # WHEN_VISIBLE
        if self._is_visible():
            timer = self._ensure_timer()
            timer.start(self.delay_ms)
        elif self._timer is not None:
            self._timer.stop()

    def _on_timer(self) -> None:
        if not self.active or self._latest is _NOTHING:
            return
        if self.delivery is Delivery.WHEN_VISIBLE and not self._is_visible():
            return  # stays dirty; replayed on show
        if self.drop_when_hidden and not self._is_visible():
            self._latest = _NOTHING
            return
        payload, self._latest = self._latest, _NOTHING
        self._call(payload)

    def shown(self, delay_ms: int = 0) -> None:
        """The owner became visible: deliver a stored payload (WHEN_VISIBLE)."""
        if not self.active or self.delivery is not Delivery.WHEN_VISIBLE or self._latest is _NOTHING:
            return
        self._ensure_timer().start(max(0, int(delay_ms)))

    def flush(self) -> None:
        """Deliver a stored payload now if the owner can take it (tests)."""
        if self._timer is not None:
            self._timer.stop()
        self._on_timer()

    def _call(self, payload: Any) -> None:
        self.bus._run_handler(self, payload)

    def cancel(self) -> None:
        self.active = False
        self._latest = _NOTHING
        if self._timer is not None:
            self._timer.stop()
            self._timer.deleteLater()
            self._timer = None
        self.handler = _dropped_handler
        self.visible = None


def _dropped_handler(_payload: Any) -> None:
    return None


class PanelBus(QtCore.QObject):
    """Typed events of one explore window; see the module docstring."""

    # A handler raised (owner uid, event, error text); after it was logged.
    handlerFailed = QtCore.pyqtSignal(str, str, str)

    def __init__(self, parent: Optional[QtCore.QObject] = None) -> None:
        super().__init__(parent)
        self._subs: Dict[str, List[Subscription]] = {event: [] for event in EVENTS}
        self._queue: Deque[Tuple[str, Any, Optional[str]]] = deque()
        self._dispatching = False
        self.delivered = 0  # handler calls (tests and diagnostics)

    # Subscriptions

    def subscribe(
        self,
        owner: str,
        event: str,
        handler: Callable[[Any], None],
        *,
        delivery: Delivery = Delivery.IMMEDIATE,
        coalesce_ms: int = 0,
        delay_ms: int = 0,
        visible: Optional[Callable[[], bool]] = None,
        drop_when_hidden: bool = False,
    ) -> Subscription:
        if event not in self._subs:
            raise ValueError(f"unknown event {event!r}")
        sub = Subscription(
            self,
            str(owner),
            event,
            handler,
            delivery,
            coalesce_ms=coalesce_ms,
            delay_ms=delay_ms,
            visible=visible,
            drop_when_hidden=drop_when_hidden,
        )
        self._subs[event].append(sub)
        return sub

    def unsubscribe(self, sub: Subscription) -> None:
        subs = self._subs.get(sub.event, [])
        if sub in subs:
            subs.remove(sub)
        sub.cancel()

    def unsubscribe_owner(self, owner: str) -> None:
        for subs in self._subs.values():
            for sub in [s for s in subs if s.owner == owner]:
                subs.remove(sub)
                sub.cancel()

    def subscriptions(self, owner: Optional[str] = None, event: Optional[str] = None) -> List[Subscription]:
        out: List[Subscription] = []
        for name, subs in self._subs.items():
            if event is not None and name != event:
                continue
            out.extend(s for s in subs if owner is None or s.owner == owner)
        return out

    def owners(self) -> List[str]:
        seen: List[str] = []
        for subs in self._subs.values():
            for sub in subs:
                if sub.owner not in seen:
                    seen.append(sub.owner)
        return seen

    # Publishing

    def publish(self, event: str, payload: Any, origin: Optional[str] = None) -> None:
        if event not in self._subs:
            raise ValueError(f"unknown event {event!r}")
        self._queue.append((event, payload, origin))
        if self._dispatching:
            return  # queued: delivered after the current dispatch, never recursively
        self._dispatching = True
        try:
            while self._queue:
                name, data, source = self._queue.popleft()
                for sub in list(self._subs[name]):
                    if not sub.active:
                        continue
                    if source is not None and sub.owner == source:
                        continue  # no echo
                    sub.offer(data)
        finally:
            self._dispatching = False

    def notify_shown(self, owner: str, delay_ms: int = 0) -> None:
        """``owner`` became visible: replay its stored WHEN_VISIBLE payloads."""
        for sub in self.subscriptions(owner):
            sub.shown(delay_ms)

    def has_pending(self) -> bool:
        """A timed delivery is scheduled for a visible subscriber."""
        return any(sub.is_pending() for subs in self._subs.values() for sub in subs)

    def flush_pending(self) -> None:
        """Deliver every scheduled delivery now (tests)."""
        for subs in self._subs.values():
            for sub in list(subs):
                if sub.is_pending():
                    sub.flush()

    def _run_handler(self, sub: Subscription, payload: Any) -> None:
        try:
            sub.handler(payload)
        except Exception as exc:  # isolate one broken panel from the fan-out
            _LOG.exception("panel %s: %s handler failed", sub.owner, sub.event)
            self.handlerFailed.emit(sub.owner, sub.event, f"{type(exc).__name__}: {exc}")
        else:
            self.delivered += 1
