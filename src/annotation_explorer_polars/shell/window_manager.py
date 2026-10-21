"""The windows of one run (plan 4.2 and section 5, ``WindowManager``).

- Explore windows: numbered (1, 2, ..., the lowest free number), each with
  its color (``theme.window_color``), tracked in most-recently-active order
  through ``QApplication.focusChanged``. While several are open, each window
  and each card shows its (owner) window's number and color. Each window is
  independent (D7): its own rail, context, sample set and panels; there is no
  API to move a panel between windows and nothing is linked across them.
- Info windows (cards): one per ``InfoKey``; reopening raises the existing
  card. Each card has an owner explore window, the target of its commands.
  When the owner closes, its cards rebind to the most recently active explore
  window; with none left they close.
- Quit policy: with ``quit_app`` set, the app quits when the last explore
  window closes (each window saves its layout in its own closeEvent first);
  info windows never keep the app alive.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, Dict, List, Optional

from PyQt6 import QtCore, QtGui, QtWidgets

from ..theme import window_color
from .commands import WINDOW_DUPLICATE, WINDOW_NEW
from ..windows.info_window import InfoKey, InfoWindow

if TYPE_CHECKING:  # pragma: no cover
    from .commands import CommandRegistry


def _alive(widget: Optional[QtCore.QObject]) -> bool:
    if widget is None:
        return False
    try:
        from PyQt6 import sip

        return not sip.isdeleted(widget)
    except Exception:
        return True


def _bring_to_front(window: QtWidgets.QWidget) -> None:
    if window.isMinimized():
        window.showNormal()
    window.show()
    window.raise_()
    window.activateWindow()


def color_swatch(color: str, size: int = 12) -> QtGui.QIcon:
    """A small rounded square of ``color`` (Window menu entries)."""
    pixmap = QtGui.QPixmap(size, size)
    pixmap.fill(QtCore.Qt.GlobalColor.transparent)
    painter = QtGui.QPainter(pixmap)
    painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
    painter.setPen(QtCore.Qt.PenStyle.NoPen)
    painter.setBrush(QtGui.QColor(color))
    painter.drawRoundedRect(QtCore.QRectF(1, 1, size - 2, size - 2), 3, 3)
    painter.end()
    return QtGui.QIcon(pixmap)


def add_menu_header(menu: QtWidgets.QMenu, text: str) -> QtGui.QAction:
    """A disabled bold entry that labels a menu section (QMenu.addSection
    draws only a line, without its text, in the app's style)."""
    action = menu.addAction(text)
    action.setEnabled(False)
    font = action.font()
    font.setBold(True)
    action.setFont(font)
    return action


class WindowManager(QtCore.QObject):
    # An explore window closed: the window, the number it had and whether it
    # was the last one (emitted before the quit callback runs).
    exploreClosed = QtCore.pyqtSignal(object, int, bool)

    def __init__(
        self,
        commands: Optional["CommandRegistry"] = None,
        *,
        window_factory: Optional[Callable[..., QtWidgets.QWidget]] = None,
        quit_app: Optional[Callable[[], None]] = None,
        parent: Optional[QtCore.QObject] = None,
    ) -> None:
        super().__init__(parent)
        self._commands = commands
        self._window_factory = window_factory
        self._quit_app = quit_app
        self._explore: List[QtWidgets.QWidget] = []  # most recently active first
        self._numbers: Dict[int, int] = {}  # id(window) -> number
        self._infos: Dict[InfoKey, InfoWindow] = {}
        app = QtWidgets.QApplication.instance()
        if isinstance(app, QtWidgets.QApplication):
            app.focusChanged.connect(self._on_focus_changed)

    # ------------------------------------------------------------ explore

    def open_explore(self, **kwargs: Any) -> QtWidgets.QWidget:
        """Build a new explore window with the factory (it registers itself);
        ``kwargs`` say how it starts (layout, filters, window subset)."""
        if self._window_factory is None:
            raise RuntimeError("This WindowManager has no explore window factory")
        window = self._window_factory(**kwargs)
        if window not in self._explore:
            self.add_explore(window)
        return window

    def next_number(self) -> int:
        """The number the next explore window gets: the lowest free one."""
        used = set(self._numbers.values())
        number = 1
        while number in used:
            number += 1
        return number

    def add_explore(self, window: QtWidgets.QWidget, number: Optional[int] = None) -> int:
        """Register an explore window; returns its number (``number`` when it
        is free, else the lowest free one)."""
        if window in self._explore:
            return self.number(window)
        if number is None or int(number) < 1 or int(number) in set(self._numbers.values()):
            number = self.next_number()
        number = int(number)
        self._numbers[id(window)] = number
        self._explore.insert(0, window)
        closed = getattr(window, "windowClosed", None)
        if closed is not None:
            closed.connect(lambda _w=window: self._on_explore_closed(_w))
        window.destroyed.connect(lambda *_a, _w=window: self._on_explore_closed(_w))
        self._refresh_badges()
        return number

    def explore_windows(self) -> List[QtWidgets.QWidget]:
        """Open explore windows in number order."""
        return sorted(self._explore, key=self.number)

    def number(self, window: QtWidgets.QWidget) -> int:
        return self._numbers.get(id(window), 0)

    def color(self, window: Optional[QtWidgets.QWidget]) -> str:
        """The window's color, or "" for a window that is not registered."""
        number = self.number(window) if window is not None else 0
        return window_color(number) if number else ""

    def several_open(self) -> bool:
        return len(self._explore) > 1

    def active_explore(self) -> Optional[QtWidgets.QWidget]:
        """The explore window that had the focus last."""
        return self._explore[0] if self._explore else None

    def note_active(self, window: QtWidgets.QWidget) -> None:
        if window in self._explore and self._explore[0] is not window:
            self._explore.remove(window)
            self._explore.insert(0, window)

    def _on_focus_changed(self, _old: Optional[QtWidgets.QWidget], new: Optional[QtWidgets.QWidget]) -> None:
        if new is None or not _alive(new):
            return
        top = new.window()
        if top in self._explore:
            self.note_active(top)

    def _on_explore_closed(self, window: QtWidgets.QWidget) -> None:
        if not _alive(self) or window not in self._explore:
            # At app quit the manager can go before a window it tracks.
            return
        self._explore.remove(window)
        number = self._numbers.pop(id(window), 0)
        successor = self.active_explore()
        for card in list(self._infos.values()):
            if card.owner() is not window or not _alive(card):
                continue
            if successor is not None:
                card.set_owner(successor)  # orphaned: rebind (plan section 5)
            else:
                card.close()
        self._refresh_badges()
        self.exploreClosed.emit(window, number, successor is None)
        if successor is not None:
            return
        # The last explore window closed: nothing may keep the app alive.
        for card in list(self._infos.values()):
            if _alive(card):
                card.close()
        if self._quit_app is not None:
            self._quit_app()

    # ------------------------------------------------------------ info cards

    def open_info(
        self,
        key: InfoKey,
        factory: Callable[[], InfoWindow],
        *,
        owner: Optional[QtWidgets.QWidget] = None,
    ) -> InfoWindow:
        """Show the card for ``key``: raise the open one, or build it. A card
        reopened from another explore window is rebound to that window."""
        owner = owner if owner is not None else self.active_explore()
        card = self._infos.get(key)
        if card is not None and _alive(card) and card.isVisible():
            if owner is not None and card.owner() is not owner:
                card.set_owner(owner, self._badge_for(owner))
            _bring_to_front(card)
            return card
        if card is not None:
            self._forget_info(key, card)
        card = factory()
        card.bind(key, owner, self._commands)
        self._show_owner_identity(card)
        card.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self._infos[key] = card
        card.finished.connect(lambda _r=0, _k=key, _c=card: self._forget_info(_k, _c))
        card.destroyed.connect(lambda *_a, _k=key, _c=card: self._forget_info(_k, _c))
        card.ownerChanged.connect(lambda _o, _c=card: self._on_card_owner_changed(_c))
        _bring_to_front(card)
        return card

    def info_window(self, key: InfoKey) -> Optional[InfoWindow]:
        card = self._infos.get(key)
        return card if card is not None and _alive(card) else None

    def info_windows(self, owner: Optional[QtWidgets.QWidget] = None) -> List[InfoWindow]:
        cards = [card for card in self._infos.values() if _alive(card)]
        if owner is not None:
            cards = [card for card in cards if card.owner() is owner]
        return cards

    def close_info(self, key: InfoKey) -> bool:
        card = self.info_window(key)
        if card is None:
            return False
        card.close()
        return True

    def _forget_info(self, key: InfoKey, card: InfoWindow) -> None:
        if self._infos.get(key) is card:
            del self._infos[key]

    def _on_card_owner_changed(self, card: InfoWindow) -> None:
        if _alive(card):
            self._show_owner_identity(card)

    def _badge_for(self, owner: Optional[QtWidgets.QWidget]) -> str:
        """The owner's number, shown only while several explore windows are open."""
        if owner is None or len(self._explore) < 2:
            return ""
        number = self.number(owner)
        return str(number) if number else ""

    def _show_owner_identity(self, card: InfoWindow) -> None:
        """The card's owner badge: number in the title, color stripe (both only
        while several explore windows are open)."""
        badge = self._badge_for(card.owner())
        card.set_owner_badge(badge)
        card.set_owner_color(self.color(card.owner()) if badge else "")

    def _refresh_badges(self) -> None:
        several = self.several_open()
        for window in list(self._explore):
            identity = getattr(window, "set_window_identity", None)
            if callable(identity) and _alive(window):
                number = self.number(window)
                identity(number, several, window_color(number))
        for card in self.info_windows():
            self._show_owner_identity(card)

    # ------------------------------------------------------------ Window menu

    def fill_window_menu(
        self,
        menu: QtWidgets.QMenu,
        current: Optional[QtWidgets.QWidget] = None,
        layouts_menu: Optional[QtWidgets.QMenu] = None,
        *,
        new_panel_action: Optional[QtGui.QAction] = None,
        filters_action: Optional[QtGui.QAction] = None,
    ) -> None:
        """The Window menu: New Explore Window and Duplicate Window (run on
        ``current``), then ``current``'s New Panel..., Layouts submenu and
        Show Filters, then the open explore windows and the open cards; each
        entry raises its window. ``current`` owns its two actions (clear()
        leaves them alive). Only New Panel... has a shortcut (Ctrl+Shift+P,
        set on the action by the window; D8 holds for the rest)."""
        menu.clear()
        commands = self._commands
        if current is not None and commands is not None and commands.has(WINDOW_NEW):
            new_action = menu.addAction("New Explore Window")
            new_action.setObjectName("NewExploreWindowAction")
            new_action.setToolTip("A new explore window on this run, with all samples and the default layout")
            new_action.triggered.connect(lambda _c=False, _w=current: commands.run(WINDOW_NEW, _w))
            duplicate = menu.addAction("Duplicate Window")
            duplicate.setObjectName("DuplicateWindowAction")
            duplicate.setToolTip(
                "A new explore window that starts with this window's filters and layout; "
                "from then on the two are independent"
            )
            duplicate.triggered.connect(lambda _c=False, _w=current: commands.run(WINDOW_DUPLICATE, _w))
            if new_panel_action is not None or layouts_menu is not None or filters_action is not None:
                menu.addSeparator()
            if new_panel_action is not None:
                menu.addAction(new_panel_action)
            if layouts_menu is not None:
                # The window owns the submenu (clear() keeps it alive).
                menu.addMenu(layouts_menu)
            if filters_action is not None:
                menu.addAction(filters_action)
            menu.addSeparator()
        add_menu_header(menu, "Explore windows")
        several = self.several_open()
        for window in self.explore_windows():
            title = window.windowTitle() or "Explore window"
            if several and not callable(getattr(window, "set_window_identity", None)):
                # Explore windows carry their number in their own title.
                title = f"{title} [{self.number(window)}]"
            action = menu.addAction(title)
            action.setCheckable(True)
            action.setChecked(window is current)
            if several:
                action.setIcon(color_swatch(self.color(window)))
            action.triggered.connect(lambda _c=False, _w=window: _bring_to_front(_w))
        menu.addSeparator()
        add_menu_header(menu, "Info cards")
        cards = self.info_windows()
        if not cards:
            empty = menu.addAction("No open cards")
            empty.setEnabled(False)
        for card in cards:
            action = menu.addAction(card.menu_label())
            if several and card.owner() is not None:
                action.setIcon(color_swatch(self.color(card.owner())))
            action.triggered.connect(lambda _c=False, _w=card: _bring_to_front(_w))
