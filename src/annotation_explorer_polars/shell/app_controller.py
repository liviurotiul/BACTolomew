"""The application object of one loaded run (plan 4.2, ``AppController``).

It replaces the app's former ``state['win']``: it owns the run's ``RunData``
and ``QueryEngine`` (shared read only by every explore window), the
``WindowManager``, the ``CommandRegistry`` and the ``LayoutStore``.

Explore windows (plan 3.10): any number of them, each fully independent (its
own rail, context, sample set and panels; D7). New Explore Window starts with
a clean rail and the default layout, Duplicate Window with a copy of another
window's filters and layout, a window opened from a clade, cluster or
plasmid with that subset locked in its rail, and a neighbourhood window
(Compare neighborhoods, a plasmid context): a copy of the requesting window's
filters and rail with only the comparison panel in its docking area. Sharing the read-only RunData
and engine is not linking: no state flows from one window to another.

Persistence (plan 3.11): each window saves its layout into the slot of its
number. The windows open when the run is left (by closing the last window, by
Quit, or by the app quitting) are recorded as the run's session; reopening
the run opens one window per recorded slot (layout only, D5).

``app.main`` makes one with ``manage_quit=True`` (the app quits when the last
explore window closes). A window built on its own (tests, scripts, the
app-check harness) makes a private controller without the quit policy.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, List, Optional, Tuple

from PyQt6 import QtCore, QtWidgets

from ..model import QueryEngine, RunData
from ..model.filter_state import FilterState, SampleSubset
from .commands import CommandRegistry, register_explore_commands
from .layout_store import START_DEFAULT_LAYOUT, START_EMPTY_LAYOUT, LayoutStore, run_key
from .window_manager import WindowManager, _bring_to_front

if TYPE_CHECKING:  # pragma: no cover
    from ..windows.explore_window import ExploreWindow

DEFAULT_WINDOW_SIZE = (1400, 900)
# A window opened from another one sits this far down and right of it.
CASCADE_OFFSET = 32


class AppController(QtCore.QObject):
    def __init__(
        self,
        run_data: RunData,
        engine: Optional[QueryEngine] = None,
        *,
        manage_quit: bool = False,
        open_run_list: Optional[Callable[["AppController", QtWidgets.QWidget], None]] = None,
        parent: Optional[QtCore.QObject] = None,
    ) -> None:
        super().__init__(parent)
        if engine is not None and engine.run is not run_data:
            raise ValueError("The engine was built for another run's RunData.")
        self.run: RunData = run_data
        self.engine: QueryEngine = engine if engine is not None else QueryEngine(run_data)
        self.layouts = LayoutStore()
        self.commands = CommandRegistry()
        self._quitting = False
        # File > Run List... (app.main only): closes the run's windows and shows
        # the run list; a controller without it has no such menu entry.
        self._open_run_list = open_run_list
        self._reloading = False
        quit_app = None
        if manage_quit:
            app = QtWidgets.QApplication.instance()
            if app is not None:
                app.setQuitOnLastWindowClosed(False)
                # The app may also quit with windows still open (for example
                # when the desktop session ends): save them first.
                app.aboutToQuit.connect(self._on_about_to_quit)
            self._quit_app = QtCore.QCoreApplication.quit
            quit_app = self._quit_unless_reloading
        self.windows = WindowManager(
            self.commands,
            window_factory=self._build_explore_window,
            quit_app=quit_app,
            parent=self,
        )
        self.windows.exploreClosed.connect(self._on_explore_closed)
        register_explore_commands(self.commands, self)

    # ------------------------------------------------------------ windows

    def _build_explore_window(self, **kwargs: Any) -> "ExploreWindow":
        from ..windows.explore_window import ExploreWindow

        window = ExploreWindow(run_data=self.run, engine=self.engine, app_controller=self, **kwargs)
        # A closed explore window is gone for good (its panels, map and
        # timers with it); the WindowManager has handled its cards by then.
        window.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose, True)
        return window

    def open_explore(
        self,
        *,
        show: bool = True,
        layout: Any = START_DEFAULT_LAYOUT,
        filters: Optional[FilterState] = None,
        window_subset: Optional[SampleSubset] = None,
        source: Optional[QtWidgets.QWidget] = None,
    ) -> "ExploreWindow":
        """A new explore window on this run. By default (New Explore Window)
        it has a clean rail and the default layout; ``layout`` may instead be
        a saved slot number or a ``WindowLayout`` to copy, ``filters`` the
        rail's starting filters and ``window_subset`` its locked subset. A
        window opened from ``source`` takes its size, offset from it."""
        window = self.windows.open_explore(layout=layout, initial_filters=filters, window_subset=window_subset)
        if show:
            self._show_new_window(window, source)
        return window

    def _show_new_window(self, window: QtWidgets.QWidget, source: Optional[QtWidgets.QWidget]) -> None:
        """Place, show and (when opened from ``source``) raise a new window."""
        # A window opened from another one always cascades from it (even
        # when it restored a saved geometry).
        if source is not None or not getattr(window, "_geometry_restored", False):
            self._place_new_window(window, source)
        window.show()
        if source is not None:
            window.raise_()
            window.activateWindow()

    def open_neighbourhood_window(self, origin: Any, /, *, beside: Optional[str] = None, **kwargs: Any) -> Any:
        """Compare neighborhoods (a gene card) or a plasmid context (a pan
        view): a new explore window holding the comparison panel. It starts
        with a copy of ``origin``'s filters (what its rail shows: searches,
        rail filters, the subset chip, the locked window subset and the panel
        window options such as Hide hypothetical), its sort and its rail
        width; its docking area holds only the comparison. It is an ordinary
        explore window: its filters change independently of ``origin`` from
        then on, and the comparison follows its own sample set.

        A comparison that is open already (same ``key``) in any window of the
        run is shown there instead (loading again when its load failed).
        ``beside`` (a panel of ``origin``) does not apply to the new window.
        ``kwargs`` are those of ``panel.open_neighbourhood`` (its ``source`` is
        the kind of comparison). Returns the
        window that shows the comparison."""
        del beside
        from ..panels.base import NEIGHBOURHOOD

        found = self.find_document(NEIGHBOURHOOD, kwargs.get("key"), prefer=origin)
        if found is not None:
            window, _uid = found
            window._open_neighbourhood_view(**kwargs)  # focuses it (and retries a failed load)
            _bring_to_front(window)
            return window
        window = self.open_explore(layout=START_EMPTY_LAYOUT, filters=origin.ctx.draft, show=False)
        clone = getattr(window, "_clone_rail_from", None)
        if callable(clone):
            clone(origin)
        uid = window._open_neighbourhood_view(**kwargs)
        if uid is None:
            # Never leave the new window empty.
            window._panel_manager.build_default()
        mark = getattr(window, "_mark_neighbourhood_window", None)
        if callable(mark):
            mark()
        else:
            window._neighbourhood_window = True
        self._show_new_window(window, origin)
        return window

    def find_document(self, kind: str, key: Any, *, prefer: Any = None) -> Optional[Tuple[Any, str]]:
        """(window, uid) of the open document panel of ``kind`` with ``key`` in
        any explore window of the run (``prefer`` first), else None."""
        if key is None:
            return None
        windows = self.windows.explore_windows()
        if prefer in windows:
            windows = [prefer] + [w for w in windows if w is not prefer]
        for window in windows:
            manager = getattr(window, "_panel_manager", None)
            uid = manager.find_document(kind, key) if manager is not None else None
            if uid is not None and manager.is_open(uid):
                return window, uid
        return None

    def duplicate_explore(self, source: Any) -> "ExploreWindow":
        """Window > Duplicate Window: a new window that starts with a copy of
        ``source``'s filters (what its rail shows, its subset chip and window
        subset included) and layout. Nothing is shared afterwards."""
        layout = source._layout_copy() if hasattr(source, "_layout_copy") else None
        window = self.open_explore(
            layout=layout if layout is not None else START_DEFAULT_LAYOUT,
            filters=source.ctx.draft,
            source=source,
        )
        if getattr(window, "_pending_dock_state", None) is not None:
            # Its dock state is copied again from ``source`` once laid out.
            window._pending_dock_source = source
        return window

    def _place_new_window(self, window: QtWidgets.QWidget, source: Optional[QtWidgets.QWidget]) -> None:
        if source is None or not source.isVisible():
            window.resize(*DEFAULT_WINDOW_SIZE)
            return
        window.resize(source.size())
        target = source.pos() + QtCore.QPoint(CASCADE_OFFSET, CASCADE_OFFSET)
        screen = source.screen()
        if screen is not None:
            area = screen.availableGeometry()
            if target.x() + window.width() > area.right() or target.y() + window.height() > area.bottom():
                target = self._free_cascade_position(window, area)
        window.move(target)

    def _free_cascade_position(self, window: QtWidgets.QWidget, area: QtCore.QRect) -> QtCore.QPoint:
        """The new window ran off the screen: restart the cascade at the
        screen's top left, on the first step no other open explore window
        sits on (so it never covers one exactly)."""
        taken = [
            other.pos()
            for other in self.windows.explore_windows()
            if other is not window and other.isVisible()
        ]
        step = QtCore.QPoint(CASCADE_OFFSET, CASCADE_OFFSET)
        for index in range(64):
            target = area.topLeft() + step * index
            if not area.contains(target):
                break
            if not any((pos - target).manhattanLength() < CASCADE_OFFSET for pos in taken):
                return target
        return area.topLeft() + step

    def open_first_window(self) -> "ExploreWindow":
        """The windows of the run's last session (one per recorded slot, the
        first numbered 1), or one window with slot 1's layout."""
        slots = self.layouts.load_session(self.run_key()) or [1]
        windows: List[Any] = []
        for slot in slots:
            window = self.open_explore(layout=slot, show=False)
            # A reopened window without saved geometry cascades from the one
            # before it instead of landing on top of it.
            if not getattr(window, "_geometry_restored", False):
                self._place_new_window(window, windows[-1] if windows else None)
            window.show()
            windows.append(window)
        first = windows[0]
        if len(windows) > 1:
            first.raise_()
            first.activateWindow()
        return first

    # ------------------------------------------------------------ session

    def run_key(self) -> str:
        return run_key(self.run.output_dir)

    @staticmethod
    def _reopens(window: Any) -> bool:
        """A window with a locked subset is not reopened with the run, nor a
        neighbourhood window (its comparison cannot be restored)."""
        if getattr(window, "_neighbourhood_window", False):
            return False
        return getattr(window, "window_subset", None) is None

    def _record_session(self, windows: List[Any]) -> None:
        slots = [self.windows.number(w) for w in windows if self._reopens(w)]
        self.layouts.save_session(self.run_key(), [slot for slot in slots if slot])

    def _on_explore_closed(self, window: Any, number: int, last: bool) -> None:
        """The last window closed: the session is that one window."""
        if self._quitting or not last:
            return
        self.layouts.save_session(self.run_key(), [number] if number and self._reopens(window) else [])

    def quit_all(self) -> None:
        """File > Quit: record every open window as the run's session, then
        close them all (each saves its layout); the last close quits the app
        (with the quit policy)."""
        if self._quitting:
            return
        self._quitting = True
        try:
            windows = self.windows.explore_windows()
            self._record_session(windows)
            for window in windows:
                window.close()
        finally:
            self._quitting = False

    # ------------------------------------------------------------ run list

    def can_open_run_list(self) -> bool:
        return self._open_run_list is not None

    def open_run_list(self, window: QtWidgets.QWidget) -> None:
        """File > Run List...: close this run's windows (like Quit, the session
        is recorded) and choose a run again; everything about runs lives there."""
        if self._open_run_list is not None:
            self._open_run_list(self, window)

    def close_for_reload(self) -> bool:
        """Close every explore window like Quit (the session is recorded)
        without quitting the app, so the run list can show. False when a
        window stayed open."""
        self._reloading = True
        try:
            self.quit_all()
        finally:
            self._reloading = False
        if self.windows.explore_windows():
            return False
        app = QtWidgets.QApplication.instance()
        if app is not None:
            try:
                app.aboutToQuit.disconnect(self._on_about_to_quit)
            except TypeError:
                pass
        return True

    def _quit_unless_reloading(self) -> None:
        if not self._reloading:
            self._quit_app()

    def _on_about_to_quit(self) -> None:
        windows = self.windows.explore_windows()
        if not windows:
            return  # left by closing the last window: recorded already
        self._record_session(windows)
        for window in windows:
            window._save_explore_layout()
