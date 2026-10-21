"""The Settings column shared by the side panels (the gene neighborhood
comparison, the Phylogeny, the cgMLST MST and the explore window's
filter rail): captions above full width controls, the 30 px drop-downs, the
cards, the chips under a panel's title, the column with its header chevron,
scrolling cards and a pinned footer, the folded strip with its vertical
caption, and the rule that folds the column while its panel is small.
"""
from __future__ import annotations

from typing import Any, Callable, List, Optional, Sequence, Tuple

from PyQt6 import QtCore, QtGui, QtWidgets

from .theme import COLORS

# One height for every control in the Settings column, one gap between its
# cards and one margin around them; the folded strip's width.
SETTINGS_CONTROL_HEIGHT = 30
SETTINGS_MARGIN = 10
SETTINGS_STRIP_WIDTH = 28


class ElidedHint(QtWidgets.QLabel):
    """A one-line hint that ends in '...' when its panel is narrow instead of
    being cut mid-word (``text()`` stays the whole hint; the tooltip shows it)."""

    def minimumSizeHint(self) -> QtCore.QSize:
        return QtCore.QSize(0, super().minimumSizeHint().height())

    def paintEvent(self, event: QtGui.QPaintEvent) -> None:
        del event
        painter = QtGui.QPainter(self)
        rect = self.contentsRect()
        text = self.fontMetrics().elidedText(self.text(), QtCore.Qt.TextElideMode.ElideRight, rect.width())
        self.style().drawItemText(
            painter, rect, int(self.alignment().value), self.palette(), self.isEnabled(), text, self.foregroundRole()
        )


class WrapLayout(QtWidgets.QLayout):
    """Lays widgets left to right and wraps whole widgets onto a new line.

    Used for the title chips, so a narrow window gets a second line
    instead of a huge minimum width.  With ``align_last_right``
    the last widget sits at the right edge of its line, like a stretch.
    """

    def __init__(
        self,
        parent: Optional[QtWidgets.QWidget] = None,
        *,
        spacing: int = 8,
        line_spacing: int = 6,
        align_last_right: bool = False,
    ) -> None:
        super().__init__(parent)
        self._items: List[QtWidgets.QLayoutItem] = []
        self._spacing = int(spacing)
        self._line_spacing = int(line_spacing)
        self._align_last_right = bool(align_last_right)

    def addItem(self, item: QtWidgets.QLayoutItem) -> None:
        self._items.append(item)

    def count(self) -> int:
        return len(self._items)

    def itemAt(self, index: int) -> Optional[QtWidgets.QLayoutItem]:
        return self._items[index] if 0 <= index < len(self._items) else None

    def takeAt(self, index: int) -> Optional[QtWidgets.QLayoutItem]:
        return self._items.pop(index) if 0 <= index < len(self._items) else None

    def expandingDirections(self) -> QtCore.Qt.Orientation:
        return QtCore.Qt.Orientation(0)

    def hasHeightForWidth(self) -> bool:
        return True

    def heightForWidth(self, width: int) -> int:
        return self._arrange(QtCore.QRect(0, 0, width, 0), apply=False)

    def setGeometry(self, rect: QtCore.QRect) -> None:
        super().setGeometry(rect)
        self._arrange(rect, apply=True)

    def sizeHint(self) -> QtCore.QSize:
        return self.minimumSize()

    def minimumSize(self) -> QtCore.QSize:
        size = QtCore.QSize()
        for item in self._visible_items():
            size = size.expandedTo(item.minimumSize().expandedTo(item.sizeHint()))
        margins = self.contentsMargins()
        return size + QtCore.QSize(margins.left() + margins.right(), margins.top() + margins.bottom())

    def _visible_items(self) -> List[QtWidgets.QLayoutItem]:
        return [item for item in self._items if not item.isEmpty()]

    def _arrange(self, rect: QtCore.QRect, *, apply: bool) -> int:
        margins = self.contentsMargins()
        area = rect.adjusted(margins.left(), margins.top(), -margins.right(), -margins.bottom())
        items = self._visible_items()
        lines: List[List[Tuple[QtWidgets.QLayoutItem, int]]] = [[]]
        x = area.x()
        for item in items:
            width = item.sizeHint().width()
            if lines[-1] and x + width > area.right() + 1:
                lines.append([])
                x = area.x()
            lines[-1].append((item, x))
            x += width + self._spacing
        y = area.y()
        for line_index, line in enumerate(lines):
            height = max((item.sizeHint().height() for item, _x in line), default=0)
            if apply:
                for position, (item, item_x) in enumerate(line):
                    hint = item.sizeHint()
                    is_last = line_index == len(lines) - 1 and position == len(line) - 1
                    if is_last and self._align_last_right:
                        item_x = max(item_x, area.right() + 1 - hint.width())
                    item.setGeometry(
                        QtCore.QRect(item_x, y + (height - hint.height()) // 2, hint.width(), hint.height())
                    )
            y += height + (self._line_spacing if line_index < len(lines) - 1 else 0)
        return y - rect.y() + margins.bottom()


class VerticalCaption(QtWidgets.QWidget):
    """A muted caption read bottom to top, for the folded Settings strip."""

    clicked = QtCore.pyqtSignal()

    def __init__(self, text: str, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self._text = text
        self.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        font = self.font()
        font.setPixelSize(11)
        font.setWeight(QtGui.QFont.Weight.DemiBold)
        self.setFont(font)
        metrics = QtGui.QFontMetrics(font)
        self.setFixedSize(metrics.height() + 4, metrics.horizontalAdvance(text) + 8)

    def text(self) -> str:
        return self._text

    def mouseReleaseEvent(self, event: QtGui.QMouseEvent) -> None:
        if event.button() == QtCore.Qt.MouseButton.LeftButton and self.rect().contains(event.position().toPoint()):
            self.clicked.emit()
        super().mouseReleaseEvent(event)

    def paintEvent(self, event: QtGui.QPaintEvent) -> None:
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.TextAntialiasing, True)
        painter.setPen(QtGui.QColor(COLORS["text_muted"]))
        painter.translate(0.0, float(self.height()))
        painter.rotate(-90.0)
        painter.drawText(
            QtCore.QRectF(0.0, 0.0, float(self.height()), float(self.width())),
            QtCore.Qt.AlignmentFlag.AlignCenter,
            self._text,
        )
        painter.end()


# ---------------------------------------------------------------- builders


def settings_caption(text: str, parent: QtWidgets.QWidget) -> QtWidgets.QLabel:
    caption = QtWidgets.QLabel(text, parent)
    caption.setProperty("role", "settingsCaption")
    return caption


def settings_control(widget: QtWidgets.QWidget, height: int = SETTINGS_CONTROL_HEIGHT) -> QtWidgets.QWidget:
    widget.setFixedHeight(int(height))
    return widget


def stacked_field(
    card: QtWidgets.QWidget, layout: QtWidgets.QVBoxLayout, caption: str, control: QtWidgets.QWidget
) -> QtWidgets.QWidget:
    """A muted caption above a control that fills the card's width."""
    group = QtWidgets.QWidget(card)
    group_layout = QtWidgets.QVBoxLayout(group)
    group_layout.setContentsMargins(0, 0, 0, 0)
    group_layout.setSpacing(4)
    control.setParent(group)
    label = settings_caption(caption, group)
    label.setBuddy(control)
    group_layout.addWidget(label)
    group_layout.addWidget(control)
    layout.addWidget(group)
    return group


def settings_chevron(parent: QtWidgets.QWidget, text: str, tooltip: str) -> QtWidgets.QToolButton:
    button = QtWidgets.QToolButton(parent)
    button.setObjectName("SettingsChevron")
    button.setText(text)
    button.setToolTip(tooltip)
    button.setAccessibleName(tooltip)
    button.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
    button.setFixedSize(24, 24)
    return button


def settings_card(parent: QtWidgets.QWidget, title: str) -> Tuple[QtWidgets.QFrame, QtWidgets.QVBoxLayout]:
    """A white card with a section title, for the Settings column."""
    card = QtWidgets.QFrame(parent)
    card.setObjectName("ContextRowsCard")
    layout = QtWidgets.QVBoxLayout(card)
    layout.setContentsMargins(12, 12, 12, 12)
    layout.setSpacing(10)
    heading = QtWidgets.QLabel(title, card)
    heading.setProperty("role", "sectionTitle")
    layout.addWidget(heading)
    return card, layout


def make_chip(text: str, tooltip: str = "", variant: str = "") -> QtWidgets.QLabel:
    chip = QtWidgets.QLabel(text)
    chip.setProperty("role", "annotationChip")
    if variant:
        chip.setProperty("variant", variant)
    if tooltip:
        chip.setToolTip(tooltip)
    return chip


def link_button(parent: QtWidgets.QWidget, text: str, tooltip: str = "") -> QtWidgets.QToolButton:
    """A text button in the primary color (Clear, Defaults, Reset view)."""
    button = QtWidgets.QToolButton(parent)
    button.setObjectName("ContextLinkButton")
    button.setText(text)
    if tooltip:
        button.setToolTip(tooltip)
    button.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
    return button


def select_combo_data(box: QtWidgets.QComboBox, value: Any) -> None:
    """Show ``value`` in ``box`` without calling its handler."""
    index = box.findData(value)
    if index >= 0 and index != box.currentIndex():
        with QtCore.QSignalBlocker(box):
            box.setCurrentIndex(index)
        box.setToolTip(box.itemData(index, QtCore.Qt.ItemDataRole.ToolTipRole) or "")


def settings_combo(
    parent: QtWidgets.QWidget,
    options: Sequence[Tuple[Any, str, str]],
    current: Any,
    on_select: Callable[[Any], None],
) -> QtWidgets.QComboBox:
    """Drop-down of the Settings column: one item per (value, text, tooltip)
    option, calling ``on_select`` with the chosen value."""
    box = QtWidgets.QComboBox(parent)
    for value, text, tooltip in options:
        box.addItem(text, value)
        box.setItemData(box.count() - 1, tooltip, QtCore.Qt.ItemDataRole.ToolTipRole)
    select_combo_data(box, current)
    box.setToolTip(box.itemData(box.currentIndex(), QtCore.Qt.ItemDataRole.ToolTipRole) or "")
    box.currentIndexChanged.connect(
        lambda index, b=box: (
            b.setToolTip(b.itemData(index, QtCore.Qt.ItemDataRole.ToolTipRole) or ""),
            on_select(b.itemData(index)),
        )
    )
    settings_control(box)
    return box


def auto_fold(small: bool, is_open: bool, auto_folded: bool, set_open: Callable[[bool], None]) -> bool:
    """Fold the column while its panel is small and open it again when the
    panel grows, if it was folded for that reason (a user's fold stands).
    Returns whether the column is now folded only because of the size."""
    if small and is_open:
        set_open(False)
        return True
    if not small and auto_folded:
        set_open(True)
        return False
    return auto_folded


# ---------------------------------------------------------------- widgets


class SettingsStrip(QtWidgets.QFrame):
    """The folded Settings column: an expand chevron, a badge (hidden until
    it has text) and the vertical caption; both open the column again."""

    expandRequested = QtCore.pyqtSignal()

    def __init__(
        self, parent: Optional[QtWidgets.QWidget], tooltip: str = "Show settings", caption: str = "Settings"
    ) -> None:
        super().__init__(parent)
        self.setObjectName("ContextSettingsStrip")
        self.setFixedWidth(SETTINGS_STRIP_WIDTH)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(2, SETTINGS_MARGIN, 2, SETTINGS_MARGIN)
        layout.setSpacing(6)
        self.expand_button = settings_chevron(self, "›", tooltip)
        self.expand_button.clicked.connect(self.expandRequested.emit)
        layout.addWidget(self.expand_button, alignment=QtCore.Qt.AlignmentFlag.AlignHCenter)
        self.badge = QtWidgets.QLabel("", self)
        self.badge.setObjectName("SettingsStripBadge")
        self.badge.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.badge.setFixedHeight(16)
        self.badge.setMinimumWidth(16)
        self.badge.hide()
        layout.addWidget(self.badge, alignment=QtCore.Qt.AlignmentFlag.AlignHCenter)
        self.caption = VerticalCaption(caption, self)
        self.caption.setToolTip(tooltip)
        self.caption.clicked.connect(self.expandRequested.emit)
        layout.addWidget(self.caption, alignment=QtCore.Qt.AlignmentFlag.AlignHCenter)
        layout.addStretch(1)
        # Its minimum is the chevron's (minimumSizeHint), not the caption's.
        layout.setSizeConstraint(QtWidgets.QLayout.SizeConstraint.SetNoConstraint)
        self._layout = layout

    def minimumSizeHint(self) -> QtCore.QSize:
        # The chevron alone: the caption steps aside in a very short panel.
        margins = self._layout.contentsMargins()
        return QtCore.QSize(SETTINGS_STRIP_WIDTH, self.expand_button.height() + margins.top() + margins.bottom())

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        super().resizeEvent(event)
        # The chevron keeps its place; the caption hides when the strip is too
        # short to show it whole (the badge and the chevron still open it).
        layout = self._layout
        margins = layout.contentsMargins()
        needed = margins.top() + margins.bottom() + self.expand_button.height()
        # The badge counted always: it shows and hides without a resize.
        needed += 2 * layout.spacing() + self.caption.height() + self.badge.height()
        self.caption.setVisible(self.height() >= needed)


class SettingsColumn(QtWidgets.QWidget):
    """The open Settings column: a header (title, chips, fold chevron), the
    cards in a scroll area no taller than they are, and a pinned footer that
    stays in view and is ruled off while the cards scroll under it."""

    collapseRequested = QtCore.pyqtSignal()

    def __init__(
        self,
        parent: Optional[QtWidgets.QWidget],
        *,
        object_name: str = "FilterRail",
        max_width: int = 420,
        collapse_tooltip: str = "Hide settings",
        title: str = "Settings",
    ) -> None:
        super().__init__(parent)
        self.setObjectName(object_name)
        self.setMaximumWidth(int(max_width))
        margin = SETTINGS_MARGIN
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        self.header = QtWidgets.QWidget(self)
        self.header_layout = QtWidgets.QHBoxLayout(self.header)
        self.header_layout.setContentsMargins(margin + 2, margin, margin - 2, 0)
        self.header_layout.setSpacing(8)
        self.title = QtWidgets.QLabel(title, self.header)
        self.title.setProperty("role", "panelTitle")
        self.collapse_button = settings_chevron(self.header, "‹", collapse_tooltip)
        self.collapse_button.clicked.connect(self.collapseRequested.emit)
        self.header_layout.addWidget(self.title)
        self.header_layout.addStretch(1)
        self.header_layout.addWidget(self.collapse_button)
        outer.addWidget(self.header)

        # The cards scroll in a short window instead of squeezing together.
        scroll = QtWidgets.QScrollArea(self)
        scroll.setObjectName("SettingsScroll")
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.viewport().setAutoFillBackground(False)
        self.scroll = scroll
        outer.addWidget(scroll, stretch=1)
        content = QtWidgets.QWidget(scroll)
        content.setObjectName("SettingsContent")
        scroll.setWidget(content)
        self.content = content
        self.content_layout = QtWidgets.QVBoxLayout(content)
        # No bottom margin: the pinned footer below adds the card gap.
        self.content_layout.setContentsMargins(margin, margin, margin, 0)
        self.content_layout.setSpacing(margin)
        self.content_layout.setAlignment(QtCore.Qt.AlignmentFlag.AlignTop)

        footer = QtWidgets.QWidget(self)
        footer.setObjectName("SettingsFooter")
        # Styled, so its top rule can show while the cards above scroll.
        footer.setAttribute(QtCore.Qt.WidgetAttribute.WA_StyledBackground, True)
        self.footer = footer
        self.footer_layout = QtWidgets.QVBoxLayout(footer)
        outer.addWidget(footer)
        # In a tall panel the footer follows the last card and the column's
        # background continues below it.  No stretch factor: it only takes
        # what the capped scroll area leaves.
        outer.addStretch(0)
        content.installEventFilter(self)
        scroll.viewport().installEventFilter(self)
        self._fit_scroll()
        # Its right edge follows the cards above when their scroll bar shows.
        self._scroll_bar_holder().installEventFilter(self)
        self._sync_footer_inset()

    def insert_header_widget(self, widget: QtWidgets.QWidget) -> None:
        """A chip or badge after the title, before the stretch."""
        widget.setParent(self.header)
        self.header_layout.insertWidget(1, widget)

    def add_card(self, card: QtWidgets.QWidget) -> None:
        self.content_layout.addWidget(card)

    def finish(self, min_width: int = 248) -> None:
        """Once every card is in: end the cards with a stretch and make the
        column wide enough for its widest card plus a scroll bar, so no card
        is ever clipped and the column never scrolls sideways."""
        self.content_layout.addStretch(1)
        self._fit_scroll()
        self._sync_footer_inset()
        needed = max(self.content.minimumSizeHint().width(), self.footer.minimumSizeHint().width())
        needed += self.scroll.verticalScrollBar().sizeHint().width()
        self.setMinimumWidth(max(int(min_width), min(self.maximumWidth(), needed)))

    def _sync_footer_inset(self) -> None:
        """Line the footer up with the scrolling cards, and rule it off from
        them while they scroll under it."""
        holder = self._scroll_bar_holder()
        # Its own shown state: the column itself may still be hidden.
        scrolling = not holder.isHidden()
        inset = holder.width() if scrolling else 0
        margin = SETTINGS_MARGIN
        # The rule is drawn in the top margin, so the card gap stays as usual.
        margins = QtCore.QMargins(margin, margin, margin + inset, margin)
        if self.footer_layout.contentsMargins() != margins:
            self.footer_layout.setContentsMargins(margins)
            # Now, not on the next pass, so the footer never shows misaligned.
            self.footer_layout.activate()
        footer = self.footer
        if footer.property("divider") != scrolling:
            footer.setProperty("divider", scrolling)
            footer.style().unpolish(footer)
            footer.style().polish(footer)

    def _fit_scroll(self) -> None:
        """No taller than its cards, so no blank band opens above the footer
        (at the column's width: a wrapped hint is taller in a narrow one)."""
        content = self.scroll.widget()
        height = content.sizeHint().height()
        if content.hasHeightForWidth():
            height = max(height, content.heightForWidth(self.scroll.viewport().width()))
        self.scroll.setMaximumHeight(max(1, height))

    def _scroll_bar_holder(self) -> QtWidgets.QWidget:
        """What the scroll area shows and hides: Qt keeps the bar in a
        container widget and toggles that, not the bar."""
        bar = self.scroll.verticalScrollBar()
        holder = bar.parentWidget()
        return bar if holder is None or holder is self.scroll else holder

    def eventFilter(self, watched: QtCore.QObject, event: QtCore.QEvent) -> bool:
        if (watched is self.content and event.type() == QtCore.QEvent.Type.LayoutRequest) or (
            watched is self.scroll.viewport() and event.type() == QtCore.QEvent.Type.Resize
        ):
            self._fit_scroll()
        elif watched is self._scroll_bar_holder() and event.type() in (
            QtCore.QEvent.Type.ShowToParent,
            QtCore.QEvent.Type.HideToParent,
            QtCore.QEvent.Type.Resize,
        ):
            self._sync_footer_inset()
        return super().eventFilter(watched, event)
