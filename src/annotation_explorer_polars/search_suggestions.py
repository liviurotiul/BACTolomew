"""Suggestion list for the left-rail "Gene/annotation contains" and "Product contains" fields.

Each suggestion shows the name on the left and, on the right, its copies across all genomes
with the number of genomes in parentheses, e.g. "blaOXA-23   66 (43)". Names are grouped the
way the filter matches them (lowercased), shown with their most common spelling. The counts
are computed once, the first time the field gets focus.
"""
from __future__ import annotations

from typing import Callable, List, Tuple

import polars as pl
from PyQt6 import QtCore, QtGui, QtWidgets

from .theme import COLORS

Suggestion = Tuple[str, int, int]  # name, copies in all genomes, genomes carrying it
_COUNT_COLOR = "#607183"
_POPUP_MIN_WIDTH = 480
_POPUP_QSS = (
    f"QTreeView {{ border: 1px solid {COLORS['border']}; border-radius: 6px; padding: 3px 0; }}"
    " QTreeView::item { padding: 2px 6px; }"
)


def count_suggestions(ann_pl: pl.DataFrame, column: str, norm_column: str) -> List[Suggestion]:
    if column not in ann_pl.columns or norm_column not in ann_pl.columns:
        return []
    sample_col = "__sample_norm" if "__sample_norm" in ann_pl.columns else "Sample"
    counts = (
        ann_pl.lazy()
        .filter(pl.col(norm_column) != "")
        .group_by(norm_column)
        .agg(
            pl.col(column).cast(pl.Utf8).str.strip_chars().mode().sort().first().alias("name"),
            pl.len().alias("copies"),
            pl.col(sample_col).n_unique().alias("genomes"),
        )
        .sort(["copies", "name"], descending=[True, False])
        .collect()
    )
    return [(str(n), int(c), int(g)) for n, c, g in counts.select(["name", "copies", "genomes"]).iter_rows()]


def attach_count_suggestions(line_edit: QtWidgets.QLineEdit, load: Callable[[], List[Suggestion]]) -> None:
    _LazySuggestions(line_edit, load)


class _LazySuggestions(QtCore.QObject):
    """Installs the completer on the field's first focus-in (before QLineEdit wires it up)."""

    def __init__(self, line_edit: QtWidgets.QLineEdit, load: Callable[[], List[Suggestion]]) -> None:
        super().__init__(line_edit)
        self._edit = line_edit
        self._load = load
        line_edit.installEventFilter(self)

    def eventFilter(self, obj: QtCore.QObject, event: QtCore.QEvent) -> bool:
        if obj is self._edit and event.type() == QtCore.QEvent.Type.FocusIn:
            self._edit.removeEventFilter(self)
            self._edit.setCompleter(_build_completer(self._edit, self._load()))
        return False


def _build_completer(line_edit: QtWidgets.QLineEdit, rows: List[Suggestion]) -> QtWidgets.QCompleter:
    model = QtGui.QStandardItemModel(0, 2, line_edit)
    muted = QtGui.QBrush(QtGui.QColor(_COUNT_COLOR))
    right = QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter
    for name, copies, genomes in rows:
        tip = f"{name}\n{copies:,} copies in {genomes} genomes"
        name_item = QtGui.QStandardItem(name)
        name_item.setToolTip(tip)
        count_item = QtGui.QStandardItem(f"{copies:,} ({genomes})")
        count_item.setTextAlignment(right)
        count_item.setForeground(muted)
        count_item.setToolTip(tip)
        model.appendRow([name_item, count_item])

    completer = QtWidgets.QCompleter(model, line_edit)
    completer.setCompletionColumn(0)
    completer.setCaseSensitivity(QtCore.Qt.CaseSensitivity.CaseInsensitive)
    completer.setFilterMode(QtCore.Qt.MatchFlag.MatchContains)
    completer.setMaxVisibleItems(12)

    popup = QtWidgets.QTreeView()
    popup.setHeaderHidden(True)
    popup.setRootIsDecorated(False)
    popup.setUniformRowHeights(True)
    popup.setAllColumnsShowFocus(True)
    popup.setMinimumWidth(_POPUP_MIN_WIDTH)  # product names are long; the rail can be narrow
    popup.setTextElideMode(QtCore.Qt.TextElideMode.ElideMiddle)  # keep the distinguishing end (e.g. "OXA-66")
    popup.setStyleSheet(_POPUP_QSS)
    completer.setPopup(popup)
    header = popup.header()
    header.setStretchLastSection(False)
    header.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.Stretch)
    header.setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
    return completer
