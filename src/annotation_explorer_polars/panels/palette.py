"""The panel palette (plan 3.6): a searchable list of panel kinds.

It always opens panels in the explore window it was invoked from. Each entry
opens the panel (default placement) or, with its "Open to the Side" button,
splits it to the right of the active group. Kinds whose data is missing for the
run are listed disabled with the reason as tooltip. The Context group lists
only the comparisons that are open (they need a gene or a plasmid to open).
"""
from __future__ import annotations

from typing import Callable, List, Optional, Sequence, Tuple

from PyQt6 import QtCore, QtWidgets

from .base import PanelSpec

# Palette groups in display order (plan 3.6).
GROUP_ORDER = ("Samples", "Annotations", "Pangenome", "Phylogeny", "Typing", "Geography", "Context")

_KIND_ROLE = QtCore.Qt.ItemDataRole.UserRole + 1
_ROW_HEIGHT = 30

# The footer: where Open and Open to the Side put a panel.
DEFAULT_HINT = "Open adds the panel as a tab of the active group; Open to the Side splits it to the right."


class PanelPalette(QtWidgets.QDialog):
    """Popup list of panels. ``opened(kind, to_side)`` is emitted on choice."""

    opened = QtCore.pyqtSignal(str, bool)

    def __init__(
        self,
        specs: Sequence[PanelSpec],
        *,
        is_open: Callable[[str], bool],
        title: str = "New panel",
        parent: Optional[QtWidgets.QWidget] = None,
        documents: Sequence[Tuple[str, str, str]] = (),
        available: Optional[Callable[[str], Tuple[bool, str]]] = None,
        hint: Optional[str] = None,
    ) -> None:
        super().__init__(parent)
        self.setObjectName("PanelPalette")
        self.setWindowTitle(title)
        self.setModal(True)
        self._specs = list(specs)
        self._is_open = is_open
        # Availability by kind (it also counts a required panel); defaults to
        # the spec's own check.
        self._available = available
        # Open document panels (uid, title, group): listed to show them; a
        # comparison is opened from a gene card or a pan view, never here.
        self._documents = list(documents)

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)

        self.search = QtWidgets.QLineEdit(self)
        self.search.setObjectName("PanelPaletteSearch")
        self.search.setPlaceholderText("Search panels...")
        self.search.setClearButtonEnabled(True)
        self.search.textChanged.connect(self._apply_filter)
        self.search.installEventFilter(self)
        layout.addWidget(self.search)

        self.tree = QtWidgets.QTreeWidget(self)
        self.tree.setObjectName("PanelPaletteList")
        self.tree.setHeaderHidden(True)
        self.tree.setColumnCount(2)
        self.tree.setRootIsDecorated(False)
        self.tree.setIndentation(10)
        self.tree.setUniformRowHeights(False)
        self.tree.itemActivated.connect(lambda item, _col: self._choose(item, False))
        layout.addWidget(self.tree, stretch=1)

        # The window says where a default Open goes (a neighbourhood window
        # splits beside its comparison instead of adding a tab).
        footer = QtWidgets.QLabel(hint or DEFAULT_HINT, self)
        footer.setObjectName("PanelPaletteHint")
        footer.setWordWrap(True)
        footer.setProperty("role", "muted")
        layout.addWidget(footer)

        self._items: List[QtWidgets.QTreeWidgetItem] = []
        self._populate()
        self.tree.header().setStretchLastSection(False)
        self.tree.header().setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.tree.header().setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.resize(460, 540)
        self._select_first_enabled()

    def _populate(self) -> None:
        by_group: dict = {}
        for spec in self._specs:
            by_group.setdefault(spec.group, []).append(spec)
        for uid, title, group in self._documents:
            by_group.setdefault(group, []).append((uid, title))
        groups = [g for g in GROUP_ORDER if g in by_group] + [g for g in by_group if g not in GROUP_ORDER]
        for group in groups:
            header = QtWidgets.QTreeWidgetItem(self.tree, [group])
            header.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled)
            font = header.font(0)
            font.setBold(True)
            header.setFont(0, font)
            header.setFirstColumnSpanned(True)
            header.setExpanded(True)
            header.setSizeHint(0, QtCore.QSize(0, _ROW_HEIGHT))
            for spec in by_group[group]:
                if isinstance(spec, tuple):
                    self._add_document_item(header, *spec)
                    continue
                available, reason = self._available(spec.kind) if self._available is not None else spec.available()
                label = spec.title + ("  (open)" if self._is_open(spec.uid) else "")
                required = self._required_title(spec)
                item = QtWidgets.QTreeWidgetItem(header, [label])
                item.setData(0, _KIND_ROLE, spec.kind)
                item.setSizeHint(0, QtCore.QSize(0, _ROW_HEIGHT))
                if available:
                    item.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsSelectable)
                    item.setToolTip(0, f"Open {spec.title}" + (f" (opens {required} too)" if required else ""))
                else:
                    item.setFlags(QtCore.Qt.ItemFlag.NoItemFlags)
                    item.setToolTip(0, reason or "Not available for this run")
                self._items.append(item)
                side = QtWidgets.QToolButton(self.tree)
                side.setObjectName("PanelPaletteSideButton")
                side.setText("Open to the Side")
                side.setAutoRaise(True)
                side.setEnabled(available)
                side.setToolTip(
                    f"Open {spec.title} to the right of the active group" + (f" (opens {required} too)" if required else "")
                    if available
                    else (reason or "Not available for this run")
                )
                side.clicked.connect(lambda _c=False, it=item: self._choose(it, True))
                self.tree.setItemWidget(item, 1, side)

    def _required_title(self, spec: PanelSpec) -> str:
        """The title of the closed panel that opening ``spec`` opens too."""
        if not spec.requires or self._is_open(spec.requires):
            return ""
        required = next((other for other in self._specs if other.uid == spec.requires), None)
        return required.title if required is not None else ""

    def _add_document_item(self, header: QtWidgets.QTreeWidgetItem, uid: str, title: str) -> None:
        item = QtWidgets.QTreeWidgetItem(header, [f"{title}  (open)"])
        item.setData(0, _KIND_ROLE, uid)
        item.setSizeHint(0, QtCore.QSize(0, _ROW_HEIGHT))
        item.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsSelectable)
        item.setToolTip(0, f"Show {title}")
        self._items.append(item)

    def items(self) -> List[QtWidgets.QTreeWidgetItem]:
        return list(self._items)

    def item_for(self, kind: str) -> Optional[QtWidgets.QTreeWidgetItem]:
        for item in self._items:
            if item.data(0, _KIND_ROLE) == kind:
                return item
        return None

    def _apply_filter(self, text: str) -> None:
        needle = text.strip().casefold()
        for index in range(self.tree.topLevelItemCount()):
            header = self.tree.topLevelItem(index)
            visible_children = 0
            group_match = bool(needle) and needle in header.text(0).casefold()
            for child_index in range(header.childCount()):
                child = header.child(child_index)
                match = not needle or group_match or needle in child.text(0).casefold()
                child.setHidden(not match)
                visible_children += int(match)
            header.setHidden(visible_children == 0)
        self._select_first_enabled()

    def _select_first_enabled(self) -> None:
        for item in self._items:
            if not item.isHidden() and item.flags() & QtCore.Qt.ItemFlag.ItemIsEnabled:
                self.tree.setCurrentItem(item)
                return
        self.tree.setCurrentItem(None)

    def _choose(self, item: Optional[QtWidgets.QTreeWidgetItem], to_side: bool) -> None:
        if item is None:
            return
        kind = item.data(0, _KIND_ROLE)
        if not kind or not (item.flags() & QtCore.Qt.ItemFlag.ItemIsEnabled):
            return
        self.opened.emit(str(kind), bool(to_side))
        self.accept()

    def eventFilter(self, obj: QtCore.QObject, event: QtCore.QEvent) -> bool:  # type: ignore[override]
        if obj is self.search and event.type() == QtCore.QEvent.Type.KeyPress:
            key = event.key()
            if key in (QtCore.Qt.Key.Key_Return, QtCore.Qt.Key.Key_Enter):
                self._choose(self.tree.currentItem(), False)
                return True
            if key in (QtCore.Qt.Key.Key_Down, QtCore.Qt.Key.Key_Up):
                QtWidgets.QApplication.sendEvent(self.tree, event)
                return True
        return super().eventFilter(obj, event)
