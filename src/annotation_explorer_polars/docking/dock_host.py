"""The docking area of an explore window (plan 3, 6.3).

``DockHost`` is the only module that imports PyQt6Ads. Everything else talks to
it through panel uids and ``Placement`` values, so the library can be swapped
for a hand-written splitter-and-tab tree (fallback 2 in plan 6.3) without
touching the window or the panels.

Panels never leave their explore window: every dock has DockWidgetFloatable
off, DockHost exposes no float, take or move-to-window call, and it redocks
any floating container ADS might still create through a path the flags miss.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional

from PyQt6 import QtCore, QtWidgets
from PyQt6 import sip

from . import ads_config

try:  # Imported here only; a missing PyQt6Ads is reported at startup (docking.require_docking).
    import PyQt6Ads as ads
except Exception as _ads_error:  # pragma: no cover - depends on the environment
    ads = None
    _ADS_IMPORT_ERROR = f"{type(_ads_error).__name__}: {_ads_error}"
else:
    _ADS_IMPORT_ERROR = ""

_LOG = logging.getLogger(__name__)

SIDES = ("left", "right", "top", "bottom")
DIRECTIONS = SIDES


def ads_available() -> bool:
    return ads is not None


def ads_import_error() -> str:
    """Why PyQt6Ads could not be imported ("" when it was)."""
    return _ADS_IMPORT_ERROR


@dataclass(frozen=True)
class Placement:
    """Where a panel goes inside the docking area.

    where: "tab" adds the panel as a tab of the target group; a side splits the
        target group (or, with ``outer``, the whole docking area) and puts the
        panel on that side.
    relative_to: a panel uid; its group is the target. None means the active
        group.
    outer: split the whole docking area instead of one group.
    """

    where: str = "tab"
    relative_to: Optional[str] = None
    outer: bool = False


def _alive(obj: object) -> bool:
    return obj is not None and not sip.isdeleted(obj)


def _same(a: object, b: object) -> bool:
    if not _alive(a) or not _alive(b):
        return False
    return sip.unwrapinstance(a) == sip.unwrapinstance(b)


class DockHost(QtWidgets.QWidget):
    """Groups, splits and tabs of one explore window; floating is disabled."""

    activeViewChanged = QtCore.pyqtSignal(str)
    viewOpened = QtCore.pyqtSignal(str)
    viewClosed = QtCore.pyqtSignal(str)
    emptyChanged = QtCore.pyqtSignal(bool)
    maximizedChanged = QtCore.pyqtSignal(bool)
    # Context menus and group buttons; the owner fills the menus (plan 3.5).
    tabMenuRequested = QtCore.pyqtSignal(str, QtCore.QPoint)
    groupMenuRequested = QtCore.pyqtSignal(str, QtCore.QPoint)
    addToGroupRequested = QtCore.pyqtSignal(str)
    floatingRedocked = QtCore.pyqtSignal(str)
    # A close (tab button, middle click) of a view registered with
    # ``custom_close``: the owner disposes the panel, then calls remove_view.
    viewCloseRequested = QtCore.pyqtSignal(str)

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        if ads is None:
            from . import require_docking

            require_docking()
        super().__init__(parent)
        self.setObjectName("DockHost")
        ads_config.apply_global_config(ads)

        self._stack = QtWidgets.QStackedLayout(self)
        self._stack.setContentsMargins(0, 0, 0, 0)
        # Parented to this plain widget (which sits in the window's splitter
        # next to the rail), never to the QMainWindow: the manager fills only
        # the docking area and never claims the rail.
        self._manager = ads.CDockManager(self)
        self._manager.setObjectName("DockManager")
        self._stack.addWidget(self._manager)
        self._empty_page: Optional[QtWidgets.QWidget] = None

        self._docks: Dict[str, "ads.CDockWidget"] = {}
        self._restoring = 0
        self._quiet = 0
        self._last_active_uid: Optional[str] = None
        self._maximized_uid: Optional[str] = None
        self._maximized_hidden: List[QtWidgets.QWidget] = []
        # The splitter sizes from before a maximize, put back when it ends.
        self._maximized_sizes: List[tuple] = []
        self._empty = True
        self._shown_once = False
        self._pending_sizes: List[tuple] = []
        self._pending_retries = 0

        self._manager.focusedDockWidgetChanged.connect(self._on_focused_dock_changed)
        # dockAreaCreated fires inside the area's constructor, before the
        # wrapper is complete; decorate once control returns to the loop.
        self._manager.dockAreaCreated.connect(self._schedule_decorate)
        self._decorate_pending = False
        self._group_buttons: set = set()
        self._manager.dockAreasAdded.connect(self._on_areas_added)
        self._manager.floatingWidgetCreated.connect(self._on_floating_widget_created)

    # ------------------------------------------------------------------ views
    def register_view(
        self,
        uid: str,
        title: str,
        frame: QtWidgets.QWidget,
        *,
        custom_close: bool = False,
    ) -> None:
        """Wrap ``frame`` in a dock; it is not placed until ``place_view``.

        ``custom_close`` (the extra instances of multi-instance panels, plan
        3.12): a close asks the owner through ``viewCloseRequested`` instead
        of hiding the dock, so the panel can be disposed and the dock deleted
        (``remove_view``). Other views are hidden on close and reopen in place.
        """
        if uid in self._docks:
            raise ValueError(f"duplicate panel uid {uid!r}")
        dock = ads.CDockWidget(self._manager, title)
        dock.setObjectName(uid)
        dock.setWidget(frame, ads.CDockWidget.eInsertMode.ForceNoScrollArea)
        ads_config.apply_dock_features(ads, dock, custom_close=custom_close)
        dock.tabWidget().installEventFilter(self)
        dock.viewToggled.connect(lambda opened, _uid=uid: self._on_view_toggled(_uid, opened))
        if custom_close:
            dock.closeRequested.connect(lambda _uid=uid: self.viewCloseRequested.emit(_uid))
        self._docks[uid] = dock

    def is_custom_close(self, uid: str) -> bool:
        dock = self._docks.get(uid)
        if dock is None:
            return False
        flag = ads.CDockWidget.DockWidgetFeature.CustomCloseHandling.value
        return bool(int(dock.features().value) & int(flag))

    def remove_view(self, uid: str) -> None:
        """Take a view out of the docking area for good and delete its dock
        (and the frame inside it). Used for closed extra instances."""
        dock = self._docks.pop(uid, None)
        if dock is None:
            return
        # As for any close: a maximized group is restored first.
        self.end_maximize()
        if self._last_active_uid == uid:
            self._last_active_uid = None
        self._quiet += 1
        try:
            if _alive(dock.dockAreaWidget()):
                self._manager.removeDockWidget(dock)
            dock.deleteLater()
        finally:
            self._quiet -= 1
        self._decorate_all_areas()
        self._update_empty()
        self.viewClosed.emit(uid)

    def view_uids(self) -> List[str]:
        return list(self._docks)

    def has_view(self, uid: str) -> bool:
        return uid in self._docks

    def is_placed(self, uid: str) -> bool:
        dock = self._docks.get(uid)
        return dock is not None and _alive(dock.dockAreaWidget())

    def is_open(self, uid: str) -> bool:
        dock = self._docks.get(uid)
        return dock is not None and self.is_placed(uid) and not dock.isClosed()

    def open_uids(self) -> List[str]:
        return [uid for uid in self._docks if self.is_open(uid)]

    def view_title(self, uid: str) -> str:
        dock = self._docks.get(uid)
        return dock.windowTitle() if dock is not None else ""

    def set_view_title(self, uid: str, title: str) -> None:
        dock = self._docks.get(uid)
        if dock is not None:
            dock.setWindowTitle(title)

    def view_is_floatable(self, uid: str) -> bool:
        dock = self._docks[uid]
        # PyQt6Ads feature values do not combine with & either; compare the ints.
        floatable = ads.CDockWidget.DockWidgetFeature.DockWidgetFloatable.value
        return bool(int(dock.features().value) & int(floatable))

    def floating_count(self) -> int:
        return len(self._manager.floatingWidgets())

    def place_view(self, uid: str, placement: Placement) -> bool:
        """Put a panel at ``placement`` (opening it if closed). False if impossible."""
        self.end_maximize()
        dock = self._docks[uid]
        target = self._target_area(placement, dock)
        where = placement.where
        if where == "tab":
            if target is None:
                self._manager.addDockWidget(ads.DockWidgetArea.LeftDockWidgetArea, dock)
            elif not _same(dock.dockAreaWidget(), target):
                self._manager.addDockWidgetTabToArea(dock, target)
        elif where in SIDES:
            own = dock.dockAreaWidget()
            if (
                target is not None
                and _same(own, target)
                and not dock.isClosed()
                and target.openDockWidgetsCount() <= 1
            ):
                # Splitting a group by its only panel would leave nothing behind.
                return False
            self._manager.addDockWidget(self._ads_area(where), dock, target)
        else:
            raise ValueError(f"unknown placement {where!r}")
        if dock.isClosed():
            # Safe now: the dock has a group, so ADS cannot float it.
            dock.toggleView(True)
        self._show_tab(dock)
        dock.setAsCurrentTab()
        self._decorate_all_areas()
        self._update_empty()
        return True

    def open_view(self, uid: str, placement: Optional[Placement] = None) -> bool:
        """Show a panel: focus it if open, reopen it in place if it was closed,
        otherwise place it (default: as a tab of the active group)."""
        dock = self._docks[uid]
        if self.is_open(uid) and placement is None:
            self.focus_view(uid)
            return True
        if placement is None and self.is_placed(uid):
            self.end_maximize()
            dock.toggleView(True)
            self._show_tab(dock)
            ok = True
        else:
            # Never toggleView(True) on a dock without a group: ADS would put it
            # in a floating container (found in Phase 1, PyQt6Ads 5.1.1).
            ok = self.place_view(uid, placement or Placement())
        if ok:
            self.focus_view(uid)
        self._update_empty()
        return ok

    def close_view(self, uid: str) -> None:
        dock = self._docks.get(uid)
        if dock is None or not self.is_open(uid):
            return
        if self.is_custom_close(uid):
            self.viewCloseRequested.emit(uid)
            return
        dock.toggleView(False)

    def focus_view(self, uid: str) -> None:
        dock = self._docks.get(uid)
        if dock is None or not self.is_open(uid):
            return
        dock.setAsCurrentTab()
        self._manager.setDockWidgetFocused(dock)
        self._last_active_uid = uid

    def set_current_tab(self, uid: str) -> None:
        """Make ``uid`` the visible tab of its group without focusing it.

        If the group's focused panel was another tab, the focus moves to
        ``uid``: ADS would otherwise keep drawing the now hidden tab as the
        focused one."""
        dock = self._docks.get(uid)
        if dock is None or not self.is_open(uid):
            return
        focused = self._manager.focusedDockWidget()
        dock.setAsCurrentTab()
        if (
            _alive(focused)
            and not _same(focused, dock)
            and _same(focused.dockAreaWidget(), dock.dockAreaWidget())
        ):
            self._manager.setDockWidgetFocused(dock)
            self._last_active_uid = uid

    def move_view(self, uid: str, placement: Placement) -> bool:
        if not self.is_open(uid):
            return False
        return self.place_view(uid, placement)

    def active_view(self) -> Optional[str]:
        dock = self._manager.focusedDockWidget()
        if _alive(dock) and self.is_open(dock.objectName()):
            return dock.objectName()
        if self._last_active_uid and self.is_open(self._last_active_uid):
            return self._last_active_uid
        for uid in self._docks:
            if self.is_open(uid):
                return uid
        return None

    # ----------------------------------------------------------------- groups
    def groups(self) -> List[List[str]]:
        out: List[List[str]] = []
        for index in range(self._manager.dockAreaCount()):
            area = self._manager.dockArea(index)
            uids = [d.objectName() for d in area.openedDockWidgets()]
            if uids:
                out.append(uids)
        return out

    def group_uids(self, uid: str) -> List[str]:
        dock = self._docks.get(uid)
        area = dock.dockAreaWidget() if dock is not None else None
        if not _alive(area):
            return []
        return [d.objectName() for d in area.openedDockWidgets()]

    def current_in_group(self, uid: str) -> Optional[str]:
        dock = self._docks.get(uid)
        area = dock.dockAreaWidget() if dock is not None else None
        current = area.currentDockWidget() if _alive(area) else None
        return current.objectName() if _alive(current) else None

    def is_current_tab(self, uid: str) -> bool:
        dock = self._docks.get(uid)
        return dock is not None and self.is_open(uid) and dock.isCurrentTab()

    def neighbor_group(self, uid: str, direction: str) -> Optional[str]:
        """A panel uid in the nearest visible group on ``direction`` of ``uid``'s group."""
        dock = self._docks.get(uid)
        area = dock.dockAreaWidget() if dock is not None else None
        if not _alive(area) or direction not in DIRECTIONS:
            return None
        src = self._area_rect(area)
        best: Optional[str] = None
        best_distance: Optional[int] = None
        for index in range(self._manager.dockAreaCount()):
            other = self._manager.dockArea(index)
            if _same(other, area) or not other.isVisible() or other.openDockWidgetsCount() == 0:
                continue
            rect = self._area_rect(other)
            if direction == "left":
                ok, distance, overlap = rect.right() < src.left(), src.left() - rect.right(), self._overlap(src.top(), src.bottom(), rect.top(), rect.bottom())
            elif direction == "right":
                ok, distance, overlap = rect.left() > src.right(), rect.left() - src.right(), self._overlap(src.top(), src.bottom(), rect.top(), rect.bottom())
            elif direction == "top":
                ok, distance, overlap = rect.bottom() < src.top(), src.top() - rect.bottom(), self._overlap(src.left(), src.right(), rect.left(), rect.right())
            else:
                ok, distance, overlap = rect.top() > src.bottom(), rect.top() - src.bottom(), self._overlap(src.left(), src.right(), rect.left(), rect.right())
            if not ok or overlap <= 0:
                continue
            if best_distance is None or distance < best_distance:
                current = other.currentDockWidget()
                if _alive(current):
                    best, best_distance = current.objectName(), distance
        return best

    def split_sizes(self, uid: str) -> List[int]:
        """Current sizes of the splitter that holds ``uid``'s group, else []."""
        dock = self._docks.get(uid)
        area = dock.dockAreaWidget() if dock is not None else None
        splitter = area.parentWidget() if _alive(area) else None
        if not isinstance(splitter, QtWidgets.QSplitter):
            return []
        return [int(s) for s in splitter.sizes()]

    def _outer_splitter(self, uid: str, outer: int) -> Optional[QtWidgets.QSplitter]:
        """The splitter ``outer`` levels above the one that holds ``uid``'s
        group (0: that one; 1: the splitter holding that splitter)."""
        dock = self._docks.get(uid)
        area = dock.dockAreaWidget() if dock is not None else None
        splitter = area.parentWidget() if _alive(area) else None
        for _level in range(max(0, int(outer))):
            if not isinstance(splitter, QtWidgets.QSplitter):
                return None
            splitter = splitter.parentWidget()
        return splitter if isinstance(splitter, QtWidgets.QSplitter) else None

    def set_split_sizes(self, uid: str, sizes: List[int], *, outer: int = 0) -> None:
        """Sizes of the splitter that holds ``uid``'s group (proportions); with
        ``outer`` the splitter that many levels above it (the rows around a
        row of groups side by side).

        Before the docking area is first shown the splitters have no real size
        yet, so the proportions are applied again right after the first show."""
        if outer:
            splitter = self._outer_splitter(uid, outer)
            if splitter is None or splitter.count() != len(sizes):
                return
            self._set_splitter_sizes(splitter, sizes)
            length = (
                splitter.width()
                if splitter.orientation() == QtCore.Qt.Orientation.Horizontal
                else splitter.height()
            )
            if not self._shown_once or length < 120 * len(sizes):
                self._pending_sizes.append((uid, [int(s) for s in sizes], int(outer)))
                if self._shown_once:
                    self._pending_retries = 0
                    QtCore.QTimer.singleShot(0, self._apply_pending_sizes)
            return
        dock = self._docks.get(uid)
        area = dock.dockAreaWidget() if dock is not None else None
        splitter = area.parentWidget() if _alive(area) else None
        if not isinstance(splitter, QtWidgets.QSplitter) or splitter.count() != len(sizes):
            return
        self._set_area_splitter_sizes(area, sizes)
        length = (
            splitter.width()
            if splitter.orientation() == QtCore.Qt.Orientation.Horizontal
            else splitter.height()
        )
        if not self._shown_once or length < 120 * len(sizes):
            # Not laid out yet (a new split gets its size a few passes later).
            self._pending_sizes.append((uid, [int(s) for s in sizes]))
            if self._shown_once:
                self._pending_retries = 0
                QtCore.QTimer.singleShot(0, self._apply_pending_sizes)

    def set_split_share(self, uid: str, share: int) -> None:
        """Give ``uid``'s group ``share`` percent of the splitter that holds
        it; the other groups of that splitter share the rest in their current
        proportions (or in the proportions still waiting to be applied)."""
        dock = self._docks.get(uid)
        area = dock.dockAreaWidget() if dock is not None else None
        splitter = area.parentWidget() if _alive(area) else None
        if not isinstance(splitter, QtWidgets.QSplitter) or splitter.count() < 2:
            return
        index = splitter.indexOf(area)
        count = splitter.count()
        weights: Optional[List[int]] = None
        keep = []
        for entry in self._pending_sizes:
            if len(entry) > 2:
                keep.append(entry)  # an outer splitter's proportions
                continue
            other = self._docks.get(entry[0])
            other_area = other.dockAreaWidget() if other is not None else None
            if _alive(other_area) and _same(other_area.parentWidget(), splitter):
                # Proportions set for this splitter before ``uid`` joined it.
                if len(entry[1]) == count - 1:
                    weights = [int(s) for s in entry[1]]
                continue
            keep.append(entry)
        self._pending_sizes = keep
        if weights is None:
            current = splitter.sizes()
            weights = [current[i] for i in range(count) if i != index]
            if not all(w > 0 for w in weights):
                weights = [1] * (count - 1)
        share = max(5, min(int(share), 95))
        total = max(1, sum(weights))
        others = iter(weights)
        sizes = [
            share * 100 if i == index else max(1, int(next(others) * (100 - share) * 100 / total))
            for i in range(count)
        ]
        self.set_split_sizes(uid, sizes)

    def _set_area_splitter_sizes(self, area, sizes: List[int]) -> None:
        # QSplitter hands surplus space out equally, not by weight, so the
        # proportions are scaled to the splitter's real length first.
        splitter = area.parentWidget()
        length = 0
        if isinstance(splitter, QtWidgets.QSplitter):
            length = (
                splitter.width()
                if splitter.orientation() == QtCore.Qt.Orientation.Horizontal
                else splitter.height()
            )
        total = max(1, sum(int(s) for s in sizes))
        scale_to = max(length, 10000)
        self._manager.setSplitterSizes(area, [max(1, int(int(s) * scale_to / total)) for s in sizes])

    @staticmethod
    def _set_splitter_sizes(splitter: QtWidgets.QSplitter, sizes: List[int]) -> None:
        # As _set_area_splitter_sizes, for a splitter that holds splitters.
        length = (
            splitter.width()
            if splitter.orientation() == QtCore.Qt.Orientation.Horizontal
            else splitter.height()
        )
        total = max(1, sum(int(s) for s in sizes))
        scale_to = max(length, 10000)
        splitter.setSizes([max(1, int(int(s) * scale_to / total)) for s in sizes])

    def showEvent(self, event) -> None:  # type: ignore[override]
        super().showEvent(event)
        if not self._shown_once:
            self._shown_once = True
            QtCore.QTimer.singleShot(0, self._apply_pending_sizes)

    def resizeEvent(self, event) -> None:  # type: ignore[override]
        super().resizeEvent(event)
        if self._shown_once and self._pending_sizes:
            QtCore.QTimer.singleShot(0, self._apply_pending_sizes)

    def _apply_pending_sizes(self) -> None:
        """Apply the proportions once each splitter has its real size (nested
        splitters get theirs a few layout passes after the window shows)."""
        if not self._pending_sizes:
            return
        waiting = []
        for entry in self._pending_sizes:
            uid, sizes = entry[0], entry[1]
            if len(entry) > 2:
                outer_splitter = self._outer_splitter(uid, entry[2])
                if outer_splitter is None or outer_splitter.count() != len(sizes):
                    continue
                length = (
                    outer_splitter.width()
                    if outer_splitter.orientation() == QtCore.Qt.Orientation.Horizontal
                    else outer_splitter.height()
                )
                if length < 120 * len(sizes):
                    waiting.append(entry)
                    continue
                self._set_splitter_sizes(outer_splitter, sizes)
                continue
            dock = self._docks.get(uid)
            area = dock.dockAreaWidget() if dock is not None else None
            splitter = area.parentWidget() if _alive(area) else None
            if not isinstance(splitter, QtWidgets.QSplitter):
                continue
            length = (
                splitter.width()
                if splitter.orientation() == QtCore.Qt.Orientation.Horizontal
                else splitter.height()
            )
            if length < 120 * len(sizes):
                waiting.append((uid, sizes))
                continue
            self._set_area_splitter_sizes(area, sizes)
        self._pending_sizes = waiting
        if waiting and self._pending_retries < 40:
            self._pending_retries += 1
            QtCore.QTimer.singleShot(25, self._apply_pending_sizes)

    def clear_pending_sizes(self) -> None:
        self._pending_sizes = []

    # ------------------------------------------------------------- maximize
    def is_maximized(self) -> bool:
        return self._maximized_uid is not None

    def maximized_view(self) -> Optional[str]:
        return self._maximized_uid

    def maximize_group(self, uid: str) -> bool:
        """Hide every other group; the rail is outside the docking area and stays."""
        self.end_maximize()
        dock = self._docks.get(uid)
        area = dock.dockAreaWidget() if dock is not None else None
        if not self.is_open(uid) or not _alive(area):
            return False
        root = self._manager.rootSplitter()
        hidden: List[QtWidgets.QWidget] = []
        sizes: List[tuple] = []
        widget: Optional[QtWidgets.QWidget] = area
        while _alive(widget) and not _same(widget, root):
            parent = widget.parentWidget()
            if isinstance(parent, QtWidgets.QSplitter):
                sizes.append((parent, parent.sizes()))
                for index in range(parent.count()):
                    sibling = parent.widget(index)
                    if not _same(sibling, widget) and not sibling.isHidden():
                        sibling.hide()
                        hidden.append(sibling)
            widget = parent
        if not hidden:
            return False
        self._maximized_uid = uid
        self._maximized_hidden = hidden
        self._maximized_sizes = sizes
        self.maximizedChanged.emit(True)
        return True

    def end_maximize(self) -> None:
        if self._maximized_uid is None:
            return
        hidden, self._maximized_hidden = self._maximized_hidden, []
        sizes, self._maximized_sizes = self._maximized_sizes, []
        self._maximized_uid = None
        for widget in hidden:
            # A panel closed while its group was hidden can leave that group
            # (or a whole splitter) empty; ADS hid it, so it stays hidden.
            if _alive(widget) and self._has_open_content(widget):
                widget.show()
        # Put the sizes from before the maximize back at once: the splitters
        # lay the re-shown groups out only later, and a layout saved right
        # now (Save Layout As, Set as Default, window close) would otherwise
        # keep the maximized group's sizes.
        for splitter, before in sizes:
            if not _alive(splitter) or splitter.count() != len(before) or not any(before):
                continue
            hidden_now = [splitter.widget(i).isHidden() for i in range(splitter.count())]
            splitter.setSizes([0 if hidden_now[i] else size for i, size in enumerate(before)])
        self.maximizedChanged.emit(False)

    def _has_open_content(self, widget: QtWidgets.QWidget) -> bool:
        area = self._as_ads(widget, ads.CDockAreaWidget)
        if area is not None:
            return area.openDockWidgetsCount() > 0
        if isinstance(widget, QtWidgets.QSplitter) or widget.inherits("QSplitter"):
            splitter = widget if isinstance(widget, QtWidgets.QSplitter) else sip.cast(widget, QtWidgets.QSplitter)
            return any(
                self._has_open_content(splitter.widget(index))
                for index in range(splitter.count())
                if splitter.widget(index) is not None
            )
        return True

    def toggle_maximize(self, uid: str) -> None:
        if self.is_maximized():
            self.end_maximize()
        else:
            self.maximize_group(uid)

    # ------------------------------------------------------------- layouts
    def clear_layout(self) -> None:
        """Take every panel out of the layout (used before rebuilding a default)."""
        self.end_maximize()
        self._quiet += 1
        try:
            for dock in self._docks.values():
                if _alive(dock.dockAreaWidget()):
                    self._manager.removeDockWidget(dock)
        finally:
            self._quiet -= 1
        self._update_empty()

    def save_layout(self) -> bytes:
        self.end_maximize()
        return bytes(self._manager.saveState())

    def restore_layout(self, state: bytes) -> bool:
        if not state:
            return False
        self.end_maximize()
        self._restoring += 1
        try:
            self.clear_pending_sizes()
            # ADS restores only docks it already knows: put every unplaced
            # dock into the layout first (as tabs of one group). Docks the
            # state does not mention end up unassigned and closed.
            self._quiet += 1
            try:
                anchor = None
                for dock in self._docks.values():
                    if _alive(dock.dockAreaWidget()):
                        anchor = dock.dockAreaWidget()
                        continue
                    if anchor is None:
                        anchor = self._manager.addDockWidget(ads.DockWidgetArea.LeftDockWidgetArea, dock)
                    else:
                        self._manager.addDockWidgetTabToArea(dock, anchor)
            finally:
                self._quiet -= 1
            ok = bool(self._manager.restoreState(QtCore.QByteArray(bytes(state))))
        except Exception:  # pragma: no cover - defensive, ADS returns False
            ok = False
        finally:
            self._restoring -= 1
        self._decorate_all_areas()
        self._update_empty()
        return ok

    def is_restoring(self) -> bool:
        return bool(self._restoring) or self._manager.isRestoringState()

    # ------------------------------------------------------------ empty state
    def set_empty_page(self, page: QtWidgets.QWidget) -> None:
        if self._empty_page is not None:
            self._stack.removeWidget(self._empty_page)
        self._empty_page = page
        self._stack.addWidget(page)
        self._update_empty()

    def is_empty(self) -> bool:
        return self._empty

    def _update_empty(self) -> None:
        empty = not self.open_uids()
        page = self._empty_page if empty and self._empty_page is not None else self._manager
        if self._stack.currentWidget() is not page:
            self._stack.setCurrentWidget(page)
        if empty != self._empty:
            self._empty = empty
            self.emptyChanged.emit(empty)

    # --------------------------------------------------------------- internals
    @staticmethod
    def _ads_area(where: str):
        areas = ads.DockWidgetArea
        return {
            "left": areas.LeftDockWidgetArea,
            "right": areas.RightDockWidgetArea,
            "top": areas.TopDockWidgetArea,
            "bottom": areas.BottomDockWidgetArea,
        }[where]

    def _active_area(self):
        uid = self.active_view()
        if uid is None:
            return None
        area = self._docks[uid].dockAreaWidget()
        return area if _alive(area) else None

    def _target_area(self, placement: Placement, dock) -> Optional[object]:
        if placement.outer:
            return None
        if placement.relative_to is not None:
            if self.is_open(placement.relative_to):
                return self._docks[placement.relative_to].dockAreaWidget()
        return self._active_area()

    def _area_rect(self, area) -> QtCore.QRect:
        top_left = area.mapTo(self._manager, QtCore.QPoint(0, 0))
        return QtCore.QRect(top_left, area.size())

    @staticmethod
    def _overlap(a0: int, a1: int, b0: int, b1: int) -> int:
        return min(a1, b1) - max(a0, b0)

    @staticmethod
    def _show_tab(dock) -> None:
        """ADS keeps a closed dock's tab hidden when the dock is taken out of
        the layout and added again (Reset Layout after a close): show it."""
        tab = dock.tabWidget()
        if not dock.isClosed() and _alive(tab) and tab.isHidden():
            tab.show()

    def _on_view_toggled(self, uid: str, opened: bool) -> None:
        if self._quiet or self.is_restoring():
            return
        if not opened and self.is_maximized():
            self.end_maximize()
        self._update_empty()
        if opened:
            self.viewOpened.emit(uid)
        else:
            self.viewClosed.emit(uid)

    def _on_focused_dock_changed(self, _old, now) -> None:
        if not _alive(now):
            return
        uid = now.objectName()
        if uid in self._docks:
            self._last_active_uid = uid
            self.activeViewChanged.emit(uid)

    def _on_areas_added(self) -> None:
        # A drop on the docking area's outer edge while a group is maximized.
        if self.is_maximized() and not self.is_restoring():
            self.end_maximize()
        # A drop can make new splitters (and dividers) without a new group.
        self._schedule_decorate()

    def _schedule_decorate(self, _area=None) -> None:
        if not self._decorate_pending:
            self._decorate_pending = True
            QtCore.QTimer.singleShot(0, self._decorate_all_areas)

    def _decorate_all_areas(self) -> None:
        self._decorate_pending = False
        for index in range(self._manager.dockAreaCount()):
            area = self._manager.dockArea(index)
            if _alive(area):
                self._decorate_area(area)
        self._hook_splitter_handles()

    def _hook_splitter_handles(self) -> None:
        """Double-clicking a divider evens out the two sides (plan 3.3)."""
        for handle in self._manager.findChildren(QtWidgets.QSplitterHandle):
            if _alive(handle) and not handle.property("bactolomewHandle"):
                handle.setProperty("bactolomewHandle", True)
                handle.installEventFilter(self)

    def even_out_divider(self, handle: QtWidgets.QSplitterHandle) -> bool:
        """Give the two sides of ``handle`` the same size (their total is kept)."""
        splitter = handle.splitter()
        if not _alive(splitter):
            return False
        for index in range(1, splitter.count()):
            if _same(splitter.handle(index), handle):
                sizes = splitter.sizes()
                total = sizes[index - 1] + sizes[index]
                sizes[index - 1] = total // 2
                sizes[index] = total - total // 2
                splitter.setSizes(sizes)
                return True
        return False

    def _decorate_area(self, area) -> None:
        """Add the "+" (open a panel here) and group-menu buttons to a group strip
        and take over its context menu and double-click."""
        title_bar = area.titleBar()
        if not _alive(title_bar) or title_bar.property("bactolomewDecorated"):
            return
        title_bar.setProperty("bactolomewDecorated", True)
        title_bar.installEventFilter(self)
        add_button = QtWidgets.QToolButton(title_bar)
        add_button.setObjectName("DockGroupAddButton")
        add_button.setText("+")
        add_button.setAutoRaise(True)
        add_button.setToolTip("Open a panel in this group")
        add_button.setAccessibleName("Open a panel in this group")
        add_button.setProperty("groupAction", "add")
        menu_button = QtWidgets.QToolButton(title_bar)
        menu_button.setObjectName("DockGroupMenuButton")
        menu_button.setText("⋯")
        menu_button.setAutoRaise(True)
        menu_button.setToolTip("Group actions")
        menu_button.setAccessibleName("Group actions")
        menu_button.setProperty("groupAction", "menu")
        for button in (add_button, menu_button):
            # The buttons are Python-made children of an ADS-made title bar:
            # SIP ties their wrappers to the title bar's short-lived wrapper,
            # so keep them here (dropped when Qt deletes the button) or their
            # connections vanish with the wrapper. The strip is resolved from
            # the button at click time; no closure holds a title bar.
            self._group_buttons.add(button)
            button.destroyed.connect(self._forget_group_button)
            button.clicked.connect(self._on_group_button_clicked)
        tabs_menu = title_bar.button(ads.TitleBarButton.TitleBarButtonTabsMenu)
        index = title_bar.indexOf(tabs_menu) if _alive(tabs_menu) else -1
        if index < 0:
            index = 2
        title_bar.insertWidget(index, add_button)
        title_bar.insertWidget(index + 1, menu_button)

    @staticmethod
    def _as_ads(obj: QtCore.QObject, cls):
        """``obj`` as the ADS class ``cls``, or None. Identified by the C++ class
        name: a Python wrapper first created through a generic path (a Python
        QProxyStyle.polish, an application event filter) is typed as a plain
        Qt class, so isinstance() alone misses ADS widgets."""
        if not _alive(obj):
            return None
        if isinstance(obj, cls):
            return obj
        if obj.metaObject().className() == f"ads::{cls.__name__}":
            return sip.cast(obj, cls)
        return None

    @QtCore.pyqtSlot()
    def _on_group_button_clicked(self) -> None:
        button = self.sender()
        if isinstance(button, QtWidgets.QToolButton) and _alive(button):
            self._on_group_button(button, str(button.property("groupAction") or ""))

    @QtCore.pyqtSlot(QtCore.QObject)
    def _forget_group_button(self, obj: QtCore.QObject) -> None:
        for button in list(self._group_buttons):
            if sip.isdeleted(button) or button is obj:
                self._group_buttons.discard(button)

    def _on_group_button(self, button: QtWidgets.QToolButton, what: str) -> None:
        title_bar = self._as_ads(button.parentWidget(), ads.CDockAreaTitleBar)
        if title_bar is None:
            return
        pos = button.mapToGlobal(QtCore.QPoint(0, button.height()))
        self._emit_for_title_bar(title_bar, what, pos)

    def _emit_for_title_bar(self, title_bar, what: str, pos: Optional[QtCore.QPoint] = None) -> None:
        if not _alive(title_bar):
            return
        area = title_bar.dockAreaWidget()
        current = area.currentDockWidget() if _alive(area) else None
        if not _alive(current):
            return
        uid = current.objectName()
        if what == "add":
            self.addToGroupRequested.emit(uid)
        else:
            self.groupMenuRequested.emit(uid, pos or QtCore.QPoint())

    def eventFilter(self, obj: QtCore.QObject, event: QtCore.QEvent) -> bool:  # type: ignore[override]
        etype = event.type()
        if etype not in (QtCore.QEvent.Type.ContextMenu, QtCore.QEvent.Type.MouseButtonDblClick):
            return super().eventFilter(obj, event)
        if (
            etype == QtCore.QEvent.Type.MouseButtonDblClick
            and _alive(obj)
            and obj.inherits("QSplitterHandle")
            and event.button() == QtCore.Qt.MouseButton.LeftButton
        ):
            handle = obj if isinstance(obj, QtWidgets.QSplitterHandle) else sip.cast(obj, QtWidgets.QSplitterHandle)
            return self.even_out_divider(handle)
        tab = self._as_ads(obj, ads.CDockWidgetTab)
        if tab is not None:
            dock = tab.dockWidget()
            if not _alive(dock):
                return False
            uid = dock.objectName()
            if etype == QtCore.QEvent.Type.ContextMenu:
                # Replaces ADS's own tab menu (Detach, Close, Close Others).
                self.tabMenuRequested.emit(uid, event.globalPos())
            elif event.button() == QtCore.Qt.MouseButton.LeftButton:
                self.toggle_maximize(uid)
            return True
        title_bar = self._as_ads(obj, ads.CDockAreaTitleBar)
        if title_bar is not None:
            if etype == QtCore.QEvent.Type.ContextMenu:
                # Replaces ADS's group menu (Detach Group, Close Group, ...).
                self._emit_for_title_bar(title_bar, "menu", event.globalPos())
            elif event.button() == QtCore.Qt.MouseButton.LeftButton:
                area = title_bar.dockAreaWidget()
                current = area.currentDockWidget() if _alive(area) else None
                if _alive(current):
                    self.toggle_maximize(current.objectName())
            return True
        return super().eventFilter(obj, event)

    def _on_floating_widget_created(self, _container) -> None:
        _LOG.warning("ADS created a floating container; redocking it into the window")
        QtCore.QTimer.singleShot(0, self._redock_floating_widgets)

    def _redock_floating_widgets(self) -> None:
        for container in list(self._manager.floatingWidgets()):
            if not _alive(container):
                continue
            for dock in list(container.dockWidgets()):
                target = None
                active = self.active_view()
                if active is not None and active != dock.objectName():
                    target = self._docks[active].dockAreaWidget()
                if _alive(target) and not container.isAncestorOf(target):
                    self._manager.addDockWidgetTabToArea(dock, target)
                else:
                    self._manager.addDockWidget(ads.DockWidgetArea.RightDockWidgetArea, dock)
                if dock.isClosed():
                    dock.toggleView(True)
                self.floatingRedocked.emit(dock.objectName())
            if _alive(container):
                container.hide()
                container.deleteLater()
        self._update_empty()
