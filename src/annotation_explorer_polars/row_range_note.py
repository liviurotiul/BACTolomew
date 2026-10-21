"""A small note over a row view that shows only part of its rows.

The Phylogeny tree shows every row at fit; zoomed in, it shows a window of
rows (the Gene presence/absence heatmap always shows every row). The note
("Rows 1 to 60 of 153") says the other rows are scrolled away, not missing.
It sits in the bottom left corner of the rows' viewport, takes no mouse
events and is hidden while every row is in view.
"""
from __future__ import annotations

import math
from typing import Optional, Tuple

from PyQt6 import QtCore, QtWidgets

from .theme import COLORS

# Pixels between the note and the viewport's edges.
NOTE_MARGIN_PX = 6


def visible_rows(top: float, rows: float, total: int) -> Optional[Tuple[int, int]]:
    """The rows (first, last; zero based, inclusive) whose centers lie in the
    span of ``rows`` rows from ``top`` (row i spans i to i + 1), or None when
    no row center is in it."""
    if total <= 0 or not (math.isfinite(top) and math.isfinite(rows)) or rows <= 0:
        return None
    first = max(0, int(math.ceil(top - 0.5 - 1e-6)))
    last = min(total - 1, int(math.floor(top + rows - 0.5 + 1e-6)))
    if last < first:
        return None
    return first, last


def row_range_text(span: Optional[Tuple[int, int]], total: int) -> str:
    """'Rows 1 to 60 of 153' for a part of the rows; '' when every row is in
    view (or there are none)."""
    if total <= 0:
        return ""
    if span is None:
        return f"No rows in view of {total}"
    first, last = span
    if first <= 0 and last >= total - 1:
        return ""
    if first == last:
        return f"Row {first + 1} of {total}"
    return f"Rows {first + 1} to {last + 1} of {total}"


class RowRangeNote(QtWidgets.QLabel):
    """The note, a child of ``host`` drawn over ``target`` (a widget inside
    ``host``, usually a scroll area's viewport). It is not a child of the
    viewport itself: a scrolled viewport moves its children along."""

    def __init__(
        self,
        host: QtWidgets.QWidget,
        target: QtWidgets.QWidget,
    ) -> None:
        super().__init__(host)
        self.setObjectName("RowRangeNote")
        self._target = target
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_TransparentForMouseEvents, True)
        self.setStyleSheet(
            "QLabel#RowRangeNote {"
            f" color: {COLORS['text_muted']};"
            " background: rgba(255, 255, 255, 220);"
            f" border: 1px solid {COLORS['border']};"
            " border-radius: 4px;"
            " padding: 1px 6px;"
            " font-size: 8pt;"
            "}"
        )
        self.hide()
        target.installEventFilter(self)
        if target is not host:
            host.installEventFilter(self)

    def set_note(self, text: str) -> None:
        text = str(text or "")
        if text != self.text():
            self.setText(text)
            self.adjustSize()
        self.setVisible(bool(text))
        if text:
            self._place()
            self.raise_()

    def eventFilter(self, obj: QtCore.QObject, event: QtCore.QEvent) -> bool:
        if event.type() in (QtCore.QEvent.Type.Resize, QtCore.QEvent.Type.Move, QtCore.QEvent.Type.Show):
            if self.isVisible():
                self._place()
        return super().eventFilter(obj, event)

    def _place(self) -> None:
        host = self.parentWidget()
        target = self._target
        if host is None or target is None:
            return
        area = QtCore.QRect(QtCore.QPoint(0, 0), target.size())
        offset = target.mapTo(host, QtCore.QPoint(0, 0)) if target is not host else QtCore.QPoint(0, 0)
        rect = area.translated(offset)
        x = rect.left() + NOTE_MARGIN_PX
        y = rect.bottom() - NOTE_MARGIN_PX - self.height() + 1
        self.move(x, max(rect.top(), y))
