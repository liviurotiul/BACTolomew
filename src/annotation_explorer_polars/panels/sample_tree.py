"""The Sample tree panel (plan 4.6): the rows of the window's one sample set.

Single instance. The tree widget and its rendering are still the window's
(``TreeControllerMixin``); this panel binds them to the context: it redraws
from each applied result, bolds the rows of hovered samples (several at once
for a clade or a plasmid) and produces hover from its own viewport.

A small strip above the tree holds the window's "Hide contigs without hits"
checkbox (feedback 2026-10-07: moved here from the Annotations panel, which
had it from the rail). It is a window filter, not an option of the panel: it
lives in the window's draft (``FilterState.hide_empty_contigs``), the panel
shows the draft's value (it follows ``ctx.draftChanged`` with signals
blocked) and a user toggle emits ``windowOptionToggled``; the window writes
its draft and applies at once.
"""
from __future__ import annotations

from typing import Any, Optional, Tuple

from PyQt6 import QtCore, QtWidgets

from .base import SAMPLE_TREE, ExplorePanel


# (FilterState field, label, objectName, tooltip) of the window filter
# checkbox in the strip.
HIDE_EMPTY_CONTIGS_ITEM: Tuple[str, str, str, str] = (
    "hide_empty_contigs",
    "Hide contigs without hits",
    "HideEmptyContigsCheck",
    "Window filter: this tree lists only the contigs that have hits under the"
    " current filters. Untick to list every contig of each sample",
)


def build_sample_tree_options(parent: Optional[QtWidgets.QWidget] = None) -> Tuple[QtWidgets.QWidget, QtWidgets.QCheckBox]:
    """The strip above the sample tree with its window filter checkbox
    (laid out like the Annotations panel's "Table shows" strip)."""
    strip = QtWidgets.QWidget(parent)
    strip.setObjectName("SampleTreeOptions")
    layout = QtWidgets.QHBoxLayout(strip)
    layout.setContentsMargins(6, 4, 4, 5)
    layout.setSpacing(8)
    _field, label, name, tooltip = HIDE_EMPTY_CONTIGS_ITEM
    check = QtWidgets.QCheckBox(label, strip)
    check.setObjectName(name)
    check.setToolTip(tooltip)
    layout.addWidget(check)
    layout.addStretch(1)
    return strip, check


class SampleTreePanel(ExplorePanel):
    KIND = SAMPLE_TREE
    TITLE = "Sample tree"
    SINGLE_INSTANCE = True

    # The window filter checkbox was toggled by the user: (FilterState field, value).
    windowOptionToggled = QtCore.pyqtSignal(str, bool)

    def __init__(self, window: Any, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setObjectName("SampleTreePanel")
        self._window = window
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.options_strip, self.hide_empty_contigs_chk = build_sample_tree_options(self)
        field = HIDE_EMPTY_CONTIGS_ITEM[0]
        self.hide_empty_contigs_chk.toggled.connect(
            lambda value, _field=field: self.windowOptionToggled.emit(_field, bool(value))
        )
        layout.addWidget(self.options_strip)
        layout.addWidget(window._build_sample_tree_panel(), stretch=1)
        self.tree: QtWidgets.QTreeWidget = window.tree
        self._hovered: Optional[str] = None
        # Hover producer: this panel's own viewport (plan 3.9; event filters
        # survive re-parenting, so hover keeps working after a move).
        self.tree.setMouseTracking(True)
        self.tree.viewport().setMouseTracking(True)
        self.tree.viewport().installEventFilter(self)

    # ------------------------------------------------------------ protocol
    def bind(self, ctx: Any, uid: str) -> None:
        super().bind(ctx, uid)
        # The checkbox follows the draft (whoever edits it: the rail, Duplicate
        # Window, a window opened from another).
        signal = getattr(ctx, "draftChanged", None)
        if signal is not None:
            signal.connect(self.show_filters)
        self.show_filters(getattr(ctx, "draft", None))

    def unbind(self) -> None:
        signal = getattr(self.ctx, "draftChanged", None)
        if signal is not None:
            try:
                signal.disconnect(self.show_filters)
            except (TypeError, RuntimeError):
                pass
        super().unbind()

    def show_filters(self, filters: Any) -> None:
        """Show the window's filter flag in the checkbox (no signals)."""
        if filters is None:
            return
        value = bool(getattr(filters, HIDE_EMPTY_CONTIGS_ITEM[0], False))
        check = self.hide_empty_contigs_chk
        if check.isChecked() != value:
            with QtCore.QSignalBlocker(check):
                check.setChecked(value)

    def apply_result(self, filters: Any, result: Any) -> None:
        del filters, result
        self._window._render_sample_tree()

    def show_hover(self, hover: Any) -> None:
        self._hovered = None
        self._window._apply_tree_sample_hover(tuple(getattr(hover, "samples", ()) or ()))

    def eventFilter(self, obj: QtCore.QObject, event: QtCore.QEvent) -> bool:  # type: ignore[override]
        if obj is self.tree.viewport():
            kind = event.type()
            if kind == QtCore.QEvent.Type.MouseMove:
                pos = event.position().toPoint() if hasattr(event, "position") else event.pos()
                self._hover_from_tree(self._window._sample_from_tree_item(self.tree.itemAt(pos)))
            elif kind == QtCore.QEvent.Type.Leave:
                self._hover_from_tree(None)
        return super().eventFilter(obj, event)

    def _hover_from_tree(self, sample: Optional[str]) -> None:
        if sample == self._hovered:
            return
        self._hovered = sample
        # The tree's own feedback (bold row) is immediate; the other panels
        # follow through the context.
        self._window._apply_tree_sample_hover((sample,) if sample else ())
        if self.ctx is not None:
            self.ctx.hover([sample] if sample else [], origin=self.uid or None)
