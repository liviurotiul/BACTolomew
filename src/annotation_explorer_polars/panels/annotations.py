"""The Annotations panel (plan 4.8, D3; multi-instance from Phase 4).

Each instance owns its table, its view mode ("Contig", "Context", "Global";
shown as Contig, Around matches and All matches) and its status line. The mode
is a display option of that instance only: changing it redraws the table (and,
for the first instance, its paired genome track) from the window's applied
result and never starts a compute. The status line shows what the table lists,
so the view mode never writes into the window's shared status bar.

The first instance (uid ``annotations``) is paired with the genome track: the
table and track scroll sync, and the window's ``table`` and
``annotations_panel`` are this instance. Further instances (Duplicate Panel)
are tables only. Every instance has its own "Table shows" box.

The strip also holds the window's "Hide hypothetical" checkbox (feedback
2026-10-07: moved here from the rail; "Hide contigs without hits" moved on to
the Sample tree panel). It is a window filter, not an option of the instance:
it lives in the window's draft (``FilterState.hide_hypothetical``), every
instance shows the same value (it follows ``ctx.draftChanged``), and closing
every Annotations panel leaves it in effect. A toggle emits
``windowOptionToggled`` and the window writes its draft.

The rendering itself is still the window's table code (``TableControllerMixin``,
which takes the panel to draw into); this class holds the per-instance state
and the panel protocol.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple

from PyQt6 import QtCore, QtWidgets

from .base import ANNOTATIONS, ExplorePanel

VIEW_MODES = ("Contig", "Context", "Global")

# (internal mode, label, tooltip) of the "Table shows" drop-down.
VIEW_MODE_ITEMS: Tuple[Tuple[str, str, str], ...] = (
    ("Contig", "Contig", "Table: annotations of the selected contig"),
    ("Context", "Around matches", "Table: genes around each match on the selected contig"),
    ("Global", "All matches", "Table: matching annotations across all samples"),
)


# (FilterState field, label, objectName, tooltip) of the window filter
# checkboxes in the strip.
WINDOW_OPTION_ITEMS: Tuple[Tuple[str, str, str, str], ...] = (
    (
        "hide_hypothetical",
        "Hide hypothetical",
        "HideHypotheticalCheck",
        "Window filter, not only this table: hypothetical proteins do not count as hits"
        " in any panel of this window. Contig and Around matches still list them",
    ),
)


class _FlowRow(QtWidgets.QLayout):
    """Left-aligned row that wraps onto more lines when the panel is narrow,
    so the strip's controls never force a wide Annotations panel. Items in a
    line are centered vertically."""

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None, *, hspacing: int = 8, vspacing: int = 4) -> None:
        super().__init__(parent)
        self._items: List[QtWidgets.QLayoutItem] = []
        self._hspacing = hspacing
        self._vspacing = vspacing

    def addItem(self, item: QtWidgets.QLayoutItem) -> None:  # type: ignore[override]
        self._items.append(item)

    def count(self) -> int:  # type: ignore[override]
        return len(self._items)

    def itemAt(self, index: int) -> Optional[QtWidgets.QLayoutItem]:  # type: ignore[override]
        return self._items[index] if 0 <= index < len(self._items) else None

    def takeAt(self, index: int) -> Optional[QtWidgets.QLayoutItem]:  # type: ignore[override]
        return self._items.pop(index) if 0 <= index < len(self._items) else None

    def expandingDirections(self) -> QtCore.Qt.Orientation:  # type: ignore[override]
        return QtCore.Qt.Orientation(0)

    def hasHeightForWidth(self) -> bool:  # type: ignore[override]
        return True

    def heightForWidth(self, width: int) -> int:  # type: ignore[override]
        return self._arrange(QtCore.QRect(0, 0, width, 0), apply=False)

    def setGeometry(self, rect: QtCore.QRect) -> None:  # type: ignore[override]
        super().setGeometry(rect)
        self._arrange(rect, apply=True)

    def sizeHint(self) -> QtCore.QSize:  # type: ignore[override]
        # Preferred: everything on one line.
        visible = [item for item in self._items if not item.isEmpty()]
        width = sum(item.sizeHint().width() for item in visible) + self._hspacing * max(0, len(visible) - 1)
        height = max((item.sizeHint().height() for item in visible), default=0)
        m = self.contentsMargins()
        return QtCore.QSize(width + m.left() + m.right(), height + m.top() + m.bottom())

    def minimumSize(self) -> QtCore.QSize:  # type: ignore[override]
        size = QtCore.QSize()
        for item in self._items:
            if not item.isEmpty():
                size = size.expandedTo(item.minimumSize())
        m = self.contentsMargins()
        return size + QtCore.QSize(m.left() + m.right(), m.top() + m.bottom())

    def _arrange(self, rect: QtCore.QRect, *, apply: bool) -> int:
        m = self.contentsMargins()
        area = rect.adjusted(m.left(), m.top(), -m.right(), -m.bottom())
        lines: List[List[QtWidgets.QLayoutItem]] = []
        line: List[QtWidgets.QLayoutItem] = []
        x = 0
        for item in self._items:
            if item.isEmpty():
                continue
            width = item.sizeHint().width()
            if line and x + self._hspacing + width > area.width():
                lines.append(line)
                line, x = [], 0
            x = width if not line else x + self._hspacing + width
            line.append(item)
        if line:
            lines.append(line)
        y = area.y()
        for index, line in enumerate(lines):
            height = max(item.sizeHint().height() for item in line)
            if apply:
                x = area.x()
                for item in line:
                    hint = item.sizeHint()
                    width = min(hint.width(), max(0, area.right() + 1 - x))
                    item.setGeometry(QtCore.QRect(x, y + (height - hint.height()) // 2, width, hint.height()))
                    x += hint.width() + self._hspacing
            y += height + (self._vspacing if index < len(lines) - 1 else 0)
        return y - area.y() + m.top() + m.bottom()


def build_view_mode_switch(parent: Optional[QtWidgets.QWidget] = None) -> Tuple[QtWidgets.QWidget, QtWidgets.QComboBox]:
    """The small "Table shows" strip that sits on top of an annotation table
    (the window filter checkboxes are added by ``build_window_options``)."""
    strip = QtWidgets.QWidget(parent)
    strip.setObjectName("TableViewSwitch")
    # A wrapping row: the window filter checkboxes (build_window_options)
    # go onto a second line when the panel is narrow.
    layout = _FlowRow(strip)
    layout.setContentsMargins(2, 0, 0, 5)

    caption = QtWidgets.QLabel("Table shows", strip)
    caption.setProperty("role", "muted")
    layout.addWidget(caption)

    box = QtWidgets.QComboBox(strip)
    box.setObjectName("TableViewBox")
    box.setAccessibleName("Table shows")
    caption.setBuddy(box)
    for mode, label, tooltip in VIEW_MODE_ITEMS:
        box.addItem(label, mode)
        box.setItemData(box.count() - 1, tooltip, QtCore.Qt.ItemDataRole.ToolTipRole)
    box.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToContents)
    layout.addWidget(box)
    return strip, box


def build_window_options(strip: QtWidgets.QWidget) -> Dict[str, QtWidgets.QCheckBox]:
    """Add the window filter checkboxes to a "Table shows" strip, after its
    box. Returns them by FilterState field."""
    layout = strip.layout()
    checks: Dict[str, QtWidgets.QCheckBox] = {}
    for field, label, name, tooltip in WINDOW_OPTION_ITEMS:
        check = QtWidgets.QCheckBox(label, strip)
        check.setObjectName(name)
        check.setToolTip(tooltip)
        if layout is not None:
            layout.addWidget(check)
        checks[field] = check
    return checks


def show_mode_in_box(box: Optional[QtWidgets.QComboBox], mode: str) -> None:
    """Show ``mode`` in a "Table shows" box without emitting its signals."""
    if box is None:
        return
    index = box.findData(str(mode))
    if index >= 0 and index != box.currentIndex():
        with QtCore.QSignalBlocker(box):
            box.setCurrentIndex(index)
    box.setToolTip(box.itemData(box.currentIndex(), QtCore.Qt.ItemDataRole.ToolTipRole) or "")


class AnnotationsPanel(ExplorePanel):
    KIND = ANNOTATIONS
    TITLE = "Annotations"
    SINGLE_INSTANCE = False

    # Local signals (not on the bus): the window redraws the table and the track.
    viewModeChanged = QtCore.pyqtSignal(str)
    statusChanged = QtCore.pyqtSignal(str)
    # A window filter checkbox was toggled by the user: (FilterState field, value).
    windowOptionToggled = QtCore.pyqtSignal(str, bool)

    def __init__(
        self,
        switch: QtWidgets.QWidget,
        body: QtWidgets.QWidget,
        parent: Optional[QtWidgets.QWidget] = None,
        *,
        margins: tuple = (4, 4, 4, 0),
        window: Any = None,
        table: Optional[QtWidgets.QTableView] = None,
        mode_box: Optional[QtWidgets.QComboBox] = None,
        window_options: Optional[Dict[str, QtWidgets.QCheckBox]] = None,
    ) -> None:
        super().__init__(parent)
        self.setObjectName("AnnotationsPanel")
        self._view_mode = "Contig"
        self._window = window
        # Per-instance table state (the window's table code fills it).
        self.table: Optional[QtWidgets.QTableView] = table
        self.model: Any = None
        self.tags_delegate: Any = None
        self.details_delegate: Any = None
        self.delegate_reset_column_count = 0
        self.text_columns: Dict[int, Tuple[int, str]] = {}
        self.fit_filter: Any = None
        self.scroll_top_on_show = False
        # (result, selection, mode) last drawn; a selection event that would
        # draw the same again is skipped.
        self.rendered_key: Optional[Tuple[Any, ...]] = None
        self._hovered_sample: Optional[str] = None
        self._mode_box = mode_box
        self._window_options: Dict[str, QtWidgets.QCheckBox] = dict(window_options or {})
        self.hide_hypo_chk: Optional[QtWidgets.QCheckBox] = self._window_options.get("hide_hypothetical")
        for field, check in self._window_options.items():
            check.toggled.connect(lambda value, _field=field: self.windowOptionToggled.emit(_field, bool(value)))
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(*margins)
        layout.setSpacing(0)
        layout.addWidget(switch)
        layout.addWidget(body, stretch=1)
        self.status_label = QtWidgets.QLabel("", self)
        self.status_label.setObjectName("AnnotationsStatus")
        self.status_label.setProperty("role", "muted")
        self.status_label.setTextFormat(QtCore.Qt.TextFormat.PlainText)
        self.status_label.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        self.status_label.setContentsMargins(2, 4, 2, 2)
        self.status_label.setSizePolicy(QtWidgets.QSizePolicy.Policy.Ignored, QtWidgets.QSizePolicy.Policy.Fixed)
        self.status_label.setVisible(False)
        layout.addWidget(self.status_label)
        if mode_box is not None:
            show_mode_in_box(mode_box, self._view_mode)
            mode_box.currentIndexChanged.connect(self._on_mode_box_changed)
        if table is not None:
            # Hover producer: this panel's own viewport (plan 3.9).
            table.setMouseTracking(True)
            table.viewport().setMouseTracking(True)
            table.viewport().installEventFilter(self)

    # ------------------------------------------------------------ view mode
    def view_mode(self) -> str:
        return self._view_mode

    def set_view_mode(self, mode: str, *, notify: bool = True) -> bool:
        """Switch the mode; returns whether it changed. ``notify`` emits
        viewModeChanged (the window redraws this table, and its track)."""
        mode = str(mode or "")
        if mode not in VIEW_MODES or mode == self._view_mode:
            return False
        self._view_mode = mode
        show_mode_in_box(self._mode_box, mode)
        self.titleDetailChanged.emit()
        if notify:
            self.viewModeChanged.emit(mode)
        return True

    def _on_mode_box_changed(self, _index: int) -> None:
        box = self._mode_box
        if box is None:
            return
        box.setToolTip(box.itemData(box.currentIndex(), QtCore.Qt.ItemDataRole.ToolTipRole) or "")
        mode = box.currentData()
        if mode:
            self.set_view_mode(str(mode))

    def title_detail(self) -> str:
        """The view mode's label ("Contig", "Around matches", "All matches")."""
        return next((label for mode, label, _tip in VIEW_MODE_ITEMS if mode == self._view_mode), "")

    def mode_box(self) -> Optional[QtWidgets.QComboBox]:
        return self._mode_box

    # ------------------------------------------------------------ window filters
    def show_filters(self, filters: Any) -> None:
        """Show the window's filter flags in the checkboxes (no signals)."""
        if filters is None:
            return
        for field, check in self._window_options.items():
            value = bool(getattr(filters, field, False))
            if check.isChecked() != value:
                with QtCore.QSignalBlocker(check):
                    check.setChecked(value)

    # ------------------------------------------------------------ status line
    def status_text(self) -> str:
        return self.status_label.text()

    def set_status_text(self, text: str) -> None:
        text = str(text or "")
        changed = text != self.status_label.text()
        self.status_label.setText(text)
        self.status_label.setToolTip(text)
        self.status_label.setVisible(bool(text))
        if changed:
            self.statusChanged.emit(text)

    # ------------------------------------------------------------ protocol
    def bind(self, ctx: Any, uid: str) -> None:
        super().bind(ctx, uid)
        # The window filter checkboxes follow the draft (whoever edits it:
        # another instance, Duplicate Window, a window opened from another).
        signal = getattr(ctx, "draftChanged", None)
        if signal is not None and self._window_options:
            signal.connect(self.show_filters)
        self.show_filters(getattr(ctx, "draft", None))

    def unbind(self) -> None:
        signal = getattr(self.ctx, "draftChanged", None)
        if signal is not None and self._window_options:
            try:
                signal.disconnect(self.show_filters)
            except (TypeError, RuntimeError):
                pass
        super().unbind()

    def is_paired(self) -> bool:
        """The first instance: drives the genome track."""
        return self._window is not None and getattr(self._window, "annotations_panel", None) is self

    def apply_result(self, filters: Any, result: Any) -> None:
        del filters, result
        if self._window is not None:
            self._window._render_annotations(self)

    def show_selection(self, selection: Any) -> None:
        # Contig and Around matches show the selected contig; All matches does
        # not depend on the selection.
        if self._window is None or self._view_mode not in {"Contig", "Context"}:
            return
        if self.rendered_key == self.render_key():
            return
        self._window._render_annotations(self)

    def show_hover(self, hover: Any) -> None:
        self.tint_samples(tuple(getattr(hover, "samples", ()) or ()))

    def tint_samples(self, samples: Sequence[str]) -> None:
        """Row tint for samples hovered elsewhere (plan 4.6, optional). Only
        All matches lists several samples; in Contig and Around matches every
        row belongs to the selected sample, so a tint would cover the table."""
        if self._view_mode != "Global":
            samples = ()
        model = self.model
        if model is not None and hasattr(model, "set_hover_samples"):
            model.set_hover_samples(samples)

    def render_key(self) -> Tuple[Any, ...]:
        ctx = self.ctx
        result = getattr(ctx, "result", None)
        selection = getattr(ctx, "selection", None)
        return (id(result), selection, self._view_mode)

    def on_shown(self) -> None:
        if self.scroll_top_on_show and self.table is not None:
            # A table redrawn while hidden starts at the top when shown.
            self.scroll_top_on_show = False
            self.table.scrollToTop()

    def save_state(self) -> Dict[str, Any]:
        return {"view_mode": self._view_mode}

    def restore_state(self, state: Dict[str, Any]) -> None:
        mode = str((state or {}).get("view_mode") or "")
        if mode not in VIEW_MODES:
            return
        self.set_view_mode(mode, notify=False)
        # A panel that has not drawn yet (a new instance restoring its
        # options before its first result) is drawn by apply_result, once.
        if self._window is not None and self.rendered_key is not None and getattr(self.ctx, "result", None) is not None:
            self._window._render_annotations(self)

    # ------------------------------------------------------------ hover out
    def eventFilter(self, obj: QtCore.QObject, event: QtCore.QEvent) -> bool:  # type: ignore[override]
        table = self.table
        if table is not None and obj is table.viewport():
            kind = event.type()
            if kind == QtCore.QEvent.Type.MouseMove:
                pos = event.position().toPoint() if hasattr(event, "position") else event.pos()
                self._hover_from_table(self._sample_at(table.indexAt(pos)))
            elif kind == QtCore.QEvent.Type.Leave:
                self._hover_from_table(None)
        return super().eventFilter(obj, event)

    def _sample_at(self, index: QtCore.QModelIndex) -> Optional[str]:
        if not index.isValid():
            return None
        model = self.table.model() if self.table is not None else None
        row_data = model.row_data(index.row()) if hasattr(model, "row_data") else None
        if not row_data:
            return None
        selection = getattr(self.ctx, "selection", None)
        raw = row_data.get("Sample", getattr(selection, "sample", None))
        sample = str(raw).strip() if raw is not None else ""
        return sample or None

    def _hover_from_table(self, sample: Optional[str]) -> None:
        if sample == self._hovered_sample:
            return
        self._hovered_sample = sample
        if self.ctx is not None:
            self.ctx.hover([sample] if sample else [], origin=self.uid or None)
