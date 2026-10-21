"""Info windows: cards about one entity (plan section 5, D6, D10).

An info window is always a separate top-level window, never docked, tabbed or
stacked into an explore window. It is a snapshot taken when it opens and never
follows filters. The ``WindowManager`` owns every card: it opens one card per
``InfoKey`` (reopening raises the existing card), binds it to the explore
window that asked for it (its owner) and rebinds it when that window closes.

A card's actions are named commands run against its owner window
(``run_command``); ``COMMANDS`` lists the ones a card kind may run. The sample
card has none: it is only a card (D10).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, ClassVar, Mapping, Optional, Tuple

from PyQt6 import QtCore, QtGui, QtWidgets

if TYPE_CHECKING:  # pragma: no cover
    from ..shell.commands import CommandRegistry

SAMPLE = "sample"
GENE = "gene"

# Height of the owner window's color stripe along a card's top edge.
OWNER_STRIPE_HEIGHT = 4


@dataclass(frozen=True)
class InfoKey:
    """Identity of a card: its kind and the entity it describes."""

    kind: str
    entity_id: Tuple[str, ...]


def _text(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none", "<na>"} else text


def _coordinate(value: Any) -> str:
    """A position as text, the same whether it came as 100, "100" or 100.0."""
    text = _text(value)
    try:
        number = float(text)
    except ValueError:
        return text
    return str(int(number)) if number.is_integer() else text


def sample_key(sample: str) -> InfoKey:
    return InfoKey(SAMPLE, (_text(sample),))


def gene_key(row: Mapping[str, Any]) -> InfoKey:
    """``("gene", sample, contig, locus or start-end-strand)`` (plan section 5).

    A row with neither a locus nor coordinates falls back to its AnnotationID,
    and without one to all of its fields, so two different genes never share
    a card."""
    sample = _text(row.get("Sample"))
    contig = _text(row.get("Contig"))
    locus = _text(row.get("Locus"))
    if not locus:
        start = _coordinate(row.get("Start"))
        stop = _coordinate(row.get("Stop"))
        if start or stop:
            locus = f"{start}-{stop}{_text(row.get('Strand'))}"
        elif _text(row.get("AnnotationID")):
            locus = f"id:{_text(row.get('AnnotationID'))}"
        else:
            fields = sorted((str(name), _text(value)) for name, value in row.items())
            locus = "row:" + "|".join(f"{name}={value}" for name, value in fields if value)
    return InfoKey(GENE, (sample, contig, locus))


class InfoWindow(QtWidgets.QDialog):
    """Base class of the cards: a top-level window (``Qt.Window``) bound to an
    owner explore window. It stays a QDialog so that Esc closes a card and the
    dialog styling applies, as before."""

    COMMANDS: ClassVar[Tuple[str, ...]] = ()

    ownerChanged = QtCore.pyqtSignal(object)

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setWindowFlags(QtCore.Qt.WindowType.Window)
        self.setWindowModality(QtCore.Qt.WindowModality.NonModal)
        self._info_key: Optional[InfoKey] = None
        self._owner: Optional[QtWidgets.QWidget] = None
        self._commands: Optional["CommandRegistry"] = None
        self._card_title = ""
        self._owner_badge = ""
        self._owner_color = ""
        # The owner window's color as a thin stripe along the top edge (only
        # while several explore windows are open; plan 3.10 and section 5).
        self._owner_stripe = QtWidgets.QFrame(self)
        self._owner_stripe.setObjectName("InfoWindowOwnerStripe")
        self._owner_stripe.setFixedHeight(OWNER_STRIPE_HEIGHT)
        self._owner_stripe.setAttribute(QtCore.Qt.WidgetAttribute.WA_TransparentForMouseEvents, True)
        self._owner_stripe.hide()

    # Identity and owner (set by the WindowManager)

    def info_key(self) -> Optional[InfoKey]:
        return self._info_key

    def owner(self) -> Optional[QtWidgets.QWidget]:
        return self._owner

    def commands(self) -> Tuple[str, ...]:
        """The commands this card can run on its owner window."""
        return tuple(self.COMMANDS)

    def bind(
        self,
        key: InfoKey,
        owner: Optional[QtWidgets.QWidget],
        commands: Optional["CommandRegistry"],
    ) -> None:
        self._info_key = key
        self._commands = commands
        self.set_owner(owner)

    def set_owner(self, owner: Optional[QtWidgets.QWidget], badge: str = "") -> None:
        changed = owner is not self._owner
        self._owner = owner
        self.set_owner_badge(badge)
        if changed:
            self._owner_changed()
            self.ownerChanged.emit(owner)

    def set_owner_badge(self, badge: str) -> None:
        """The owner window's number, shown in the title when several explore
        windows are open (plan section 5)."""
        self._owner_badge = str(badge or "")
        self._apply_title()

    def owner_badge(self) -> str:
        return self._owner_badge

    def set_owner_color(self, color: str) -> None:
        """The owner window's color (a stripe along the top), or "" for none."""
        self._owner_color = str(color or "")
        stripe = self._owner_stripe
        if self._owner_color:
            stripe.setStyleSheet(f"QFrame#InfoWindowOwnerStripe {{ background: {self._owner_color}; border: none; }}")
            self._place_owner_stripe()
            stripe.show()
            stripe.raise_()
        else:
            stripe.hide()

    def owner_color(self) -> str:
        return self._owner_color

    def _place_owner_stripe(self) -> None:
        self._owner_stripe.setGeometry(0, 0, self.width(), OWNER_STRIPE_HEIGHT)

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        super().resizeEvent(event)
        self._place_owner_stripe()
        if self._owner_stripe.isVisible():
            self._owner_stripe.raise_()

    def run_command(self, name: str, **kwargs: Any) -> bool:
        """Run one of this card's commands on its owner window."""
        if name not in self.COMMANDS:
            raise ValueError(f"{type(self).__name__} has no command {name!r}")
        if self._owner is None or self._commands is None:
            return False
        self._commands.run(name, self._owner, **kwargs)
        return True

    def can_run_commands(self) -> bool:
        return self._owner is not None and self._commands is not None

    def _owner_changed(self) -> None:
        """Subclass hook: the owner window changed (or the first one was set)."""

    # Title

    def set_card_title(self, title: str) -> None:
        self._card_title = str(title or "")
        self._apply_title()

    def card_title(self) -> str:
        return self._card_title

    def menu_label(self) -> str:
        """The card's entry in the Window menu (its title, owner badge last)."""
        return self.windowTitle() or self._card_title

    def _apply_title(self) -> None:
        badge = f" [{self._owner_badge}]" if self._owner_badge else ""
        self.setWindowTitle(f"{self._card_title}{badge}")

