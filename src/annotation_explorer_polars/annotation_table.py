from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple

import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets

from .resistance_mutation import MUTATION_TAG, mutation_note, paint_mutation_badge
from .theme import COLORS


def _tag_color(tag: str) -> str:
    colors = {
        "A": "#c62828",  # AMR
        "B": "#1b5e20",  # BMRG
        "V": "#6a1b9a",  # Virulence
        "M": "#ff8f00",  # Mobile
        "P": "#1565c0",  # Phage
        "I": "#f9a825",  # Compound overlap
        MUTATION_TAG: "#c62828",  # Resistance mutation: the AMR red letter, then the badge
    }
    return colors.get(tag.upper(), "#000000")


# Row grouping without separator rows: in the all-matches view each sample's rows share
# a band colour that alternates per sample, and a hairline marks where a new sample starts.
GROUP_START_ROLE = QtCore.Qt.ItemDataRole.UserRole + 21
_BAND_COLORS = (QtGui.QColor("#ffffff"), QtGui.QColor("#f3f6fa"))
_GROUP_RULE_COLOR = QtGui.QColor("#b9c6d4")
# Gap rows in the around-matches view ("N genes not shown"): quiet, not a black bar.
_GAP_BACKGROUND = QtGui.QColor("#f7f9fb")
_GAP_TEXT = QtGui.QColor("#8593a3")
GAP_ROW_HEIGHT = 22
# Rows of samples hovered in another panel of the window (plan 4.6, optional
# table tint): the theme's soft primary.
_HOVER_TINT = QtGui.QColor(COLORS["primary_soft"])


def paint_group_rule(painter: QtGui.QPainter, option: QtWidgets.QStyleOptionViewItem, index: QtCore.QModelIndex) -> None:
    if not index.data(GROUP_START_ROLE):
        return
    painter.save()
    painter.setPen(QtGui.QPen(_GROUP_RULE_COLOR, 1))
    painter.drawLine(option.rect.left(), option.rect.top(), option.rect.right(), option.rect.top())
    painter.restore()


class RowRuleDelegate(QtWidgets.QStyledItemDelegate):
    """Default delegate for the annotation table: normal cells plus the sample-start rule."""

    def paint(
        self,
        painter: QtGui.QPainter,
        option: QtWidgets.QStyleOptionViewItem,
        index: QtCore.QModelIndex,
    ) -> None:
        super().paint(painter, option, index)
        paint_group_rule(painter, option, index)


class TagsDelegate(QtWidgets.QStyledItemDelegate):
    def paint(
        self,
        painter: QtGui.QPainter,
        option: QtWidgets.QStyleOptionViewItem,
        index: QtCore.QModelIndex,
    ) -> None:
        opt = QtWidgets.QStyleOptionViewItem(option)
        self.initStyleOption(opt, index)
        tags = index.data(QtCore.Qt.ItemDataRole.UserRole) or ""
        opt.text = ""
        style = QtWidgets.QApplication.style() if opt.widget is None else opt.widget.style()
        style.drawControl(QtWidgets.QStyle.ControlElement.CE_ItemViewItem, opt, painter, opt.widget)
        paint_group_rule(painter, option, index)

        letters = [ch for ch in str(tags) if ch.strip()]
        if not letters:
            return

        painter.save()
        font = opt.font
        font.setBold(True)
        painter.setFont(font)
        metrics = QtGui.QFontMetrics(font)
        spacing = 6
        badge = max(10, min(14, metrics.height() - 2))

        badge_gap = 2

        def _width(ch: str) -> int:
            width = metrics.horizontalAdvance(ch.upper())
            return width + badge_gap + badge if ch.upper() == MUTATION_TAG else width

        total_width = sum(_width(ch) for ch in letters)
        total_width += spacing * max(0, len(letters) - 1)
        x = opt.rect.x() + max(0, (opt.rect.width() - total_width) // 2)
        y = opt.rect.y() + (opt.rect.height() + metrics.ascent() - metrics.descent()) // 2
        for ch in letters:
            painter.setPen(QtGui.QColor(_tag_color(ch)))
            painter.drawText(QtCore.QPoint(x, y), ch.upper())
            if ch.upper() == MUTATION_TAG:
                # The X is red like the AMR A, and the mutation badge follows it.
                left = x + metrics.horizontalAdvance(ch.upper()) + badge_gap
                top = opt.rect.y() + (opt.rect.height() - badge) / 2.0
                paint_mutation_badge(painter, QtCore.QRectF(left, top, badge, badge))
            x += _width(ch) + spacing
        painter.restore()


class DetailsDelegate(QtWidgets.QStyledItemDelegate):
    clicked = QtCore.pyqtSignal(int)

    def paint(
        self,
        painter: QtGui.QPainter,
        option: QtWidgets.QStyleOptionViewItem,
        index: QtCore.QModelIndex,
    ) -> None:
        if index.model().data(index, QtCore.Qt.ItemDataRole.DisplayRole) is None:
            return super().paint(painter, option, index)

        opt = QtWidgets.QStyleOptionViewItem(option)
        self.initStyleOption(opt, index)  # picks up the row's band / hit background
        opt.text = ""
        style = QtWidgets.QApplication.style() if opt.widget is None else opt.widget.style()
        style.drawControl(QtWidgets.QStyle.ControlElement.CE_ItemViewItem, opt, painter, opt.widget)
        paint_group_rule(painter, option, index)

        rect = option.rect.adjusted(6, 6, -6, -6)
        if rect.width() <= 6 or rect.height() <= 6:
            return
        fold = min(rect.width(), rect.height()) * 0.35
        pen = QtGui.QPen(QtGui.QColor("#4d4d4d"))
        pen.setWidthF(1.0)
        painter.save()
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        painter.setPen(pen)
        painter.setBrush(QtGui.QBrush(QtGui.QColor("#ffffff")))

        path = QtGui.QPainterPath()
        path.moveTo(rect.left(), rect.top())
        path.lineTo(rect.right() - fold, rect.top())
        path.lineTo(rect.right(), rect.top() + fold)
        path.lineTo(rect.right(), rect.bottom())
        path.lineTo(rect.left(), rect.bottom())
        path.closeSubpath()
        painter.drawPath(path)

        painter.drawLine(
            QtCore.QPointF(rect.right() - fold, rect.top()),
            QtCore.QPointF(rect.right() - fold, rect.top() + fold),
        )
        painter.drawLine(
            QtCore.QPointF(rect.right() - fold, rect.top() + fold),
            QtCore.QPointF(rect.right(), rect.top() + fold),
        )
        painter.restore()

    def editorEvent(
        self,
        event: QtCore.QEvent,
        model: QtCore.QAbstractItemModel,
        option: QtWidgets.QStyleOptionViewItem,
        index: QtCore.QModelIndex,
    ) -> bool:
        if event.type() == QtCore.QEvent.Type.MouseButtonRelease:
            if hasattr(model, "row_data"):
                row_data = model.row_data(index.row())
                if row_data is None:
                    return False
            mouse_event = event  # type: ignore[assignment]
            if hasattr(mouse_event, "button") and mouse_event.button() == QtCore.Qt.MouseButton.LeftButton:
                pos = mouse_event.position().toPoint() if hasattr(mouse_event, "position") else mouse_event.pos()
                if option.rect.contains(pos):
                    self.clicked.emit(index.row())
                    return True
        return super().editorEvent(event, model, option, index)


class AnnotationTableModel(QtCore.QAbstractTableModel):
    def __init__(
        self,
        rows: List[Dict[str, Any]],
        columns: Sequence[Tuple[str, str]],
        parent: Optional[QtCore.QObject] = None,
    ) -> None:
        super().__init__(parent)
        self._rows = rows
        self._columns = list(columns)
        self.group_rows = [i for i, r in enumerate(rows) if r.get("type") == "group"]
        self.separator_rows = [i for i, r in enumerate(rows) if r.get("type") == "separator"]
        self._hover_samples: frozenset = frozenset()
        # The Tags header's tooltip: the tags' databases not installed for the run.
        self.tags_header_tooltip = ""

    def set_hover_samples(self, samples: Sequence[str]) -> None:
        """Tint the rows of these samples (hovered in another panel)."""
        wanted = frozenset(str(s).strip() for s in (samples or ()) if str(s).strip())
        if wanted == self._hover_samples:
            return
        self._hover_samples = wanted
        if self._rows and self._columns:
            self.dataChanged.emit(
                self.index(0, 0),
                self.index(len(self._rows) - 1, len(self._columns) - 1),
                [QtCore.Qt.ItemDataRole.BackgroundRole],
            )

    def hover_samples(self) -> frozenset:
        return self._hover_samples

    def rowCount(self, _parent: QtCore.QModelIndex = QtCore.QModelIndex()) -> int:
        return len(self._rows)

    def columnCount(self, _parent: QtCore.QModelIndex = QtCore.QModelIndex()) -> int:
        return len(self._columns)

    def headerData(self, section: int, orientation: QtCore.Qt.Orientation, role: int = QtCore.Qt.ItemDataRole.DisplayRole):
        if (
            role == QtCore.Qt.ItemDataRole.ToolTipRole
            and orientation == QtCore.Qt.Orientation.Horizontal
            and 0 <= section < len(self._columns)
            and self._columns[section][0] == "Tags"
        ):
            return self.tags_header_tooltip or None
        if role != QtCore.Qt.ItemDataRole.DisplayRole:
            return None
        if orientation == QtCore.Qt.Orientation.Horizontal:
            if 0 <= section < len(self._columns):
                return self._columns[section][1]
        return None

    def flags(self, index: QtCore.QModelIndex) -> QtCore.Qt.ItemFlags:
        if not index.isValid():
            return QtCore.Qt.ItemFlag.NoItemFlags
        row = self._rows[index.row()]
        if row.get("type") == "separator":
            return QtCore.Qt.ItemFlag.NoItemFlags
        if row.get("type") == "group":
            return QtCore.Qt.ItemFlag.ItemIsEnabled
        return QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsSelectable

    def data(self, index: QtCore.QModelIndex, role: int = QtCore.Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        row_info = self._rows[index.row()]
        col_name = self._columns[index.column()][0] if 0 <= index.column() < len(self._columns) else None

        if row_info.get("type") == "separator":
            if role == QtCore.Qt.ItemDataRole.DisplayRole:
                return row_info.get("label", "") if index.column() == 0 else ""
            if role == QtCore.Qt.ItemDataRole.BackgroundRole:
                return _GAP_BACKGROUND
            if role == QtCore.Qt.ItemDataRole.ForegroundRole:
                return _GAP_TEXT
            if role == QtCore.Qt.ItemDataRole.TextAlignmentRole:
                return QtCore.Qt.AlignmentFlag.AlignCenter
            if role == QtCore.Qt.ItemDataRole.FontRole:
                font = QtGui.QFont()
                font.setItalic(True)
                font.setPointSizeF(max(7.0, font.pointSizeF() - 1))
                return font
            if role == QtCore.Qt.ItemDataRole.SizeHintRole:
                return QtCore.QSize(0, GAP_ROW_HEIGHT)
            return None

        if row_info.get("type") == "group":
            if role == QtCore.Qt.ItemDataRole.DisplayRole:
                return row_info.get("label") if index.column() == 0 else ""
            if role == QtCore.Qt.ItemDataRole.FontRole:
                font = QtGui.QFont()
                font.setBold(True)
                return font
            if role == QtCore.Qt.ItemDataRole.BackgroundRole:
                return QtGui.QColor("#f0f0f0")
            return None

        row_data = row_info.get("data") or {}
        is_hit = bool(row_info.get("is_hit", False))

        if role == QtCore.Qt.ItemDataRole.BackgroundRole:
            if self._hover_samples and str(row_data.get("Sample") or "").strip() in self._hover_samples:
                return _HOVER_TINT
            if is_hit:
                return QtGui.QColor("#f5e7c6")
            band = row_info.get("band")
            if band is not None:
                return _BAND_COLORS[int(band) % 2]
        if role == GROUP_START_ROLE:
            return bool(row_info.get("group_start"))

        if col_name == "Tags":
            if role == QtCore.Qt.ItemDataRole.UserRole:
                val = row_data.get(col_name, "")
                return "" if pd.isna(val) else str(val)
            if role == QtCore.Qt.ItemDataRole.DisplayRole:
                return ""
            if role == QtCore.Qt.ItemDataRole.TextAlignmentRole:
                return QtCore.Qt.AlignmentFlag.AlignCenter
            if role == QtCore.Qt.ItemDataRole.ToolTipRole:
                return mutation_note(row_data) or None
            return None

        if col_name == "Strand":
            if role == QtCore.Qt.ItemDataRole.TextAlignmentRole:
                return QtCore.Qt.AlignmentFlag.AlignCenter
            if role == QtCore.Qt.ItemDataRole.DisplayRole:
                val = row_data.get(col_name, "")
                return "" if pd.isna(val) else str(val)
            return None

        if col_name == "Details":
            if role == QtCore.Qt.ItemDataRole.DisplayRole:
                return ""
            if role == QtCore.Qt.ItemDataRole.TextAlignmentRole:
                return QtCore.Qt.AlignmentFlag.AlignCenter
            return None

        if role == QtCore.Qt.ItemDataRole.DisplayRole:
            val = row_data.get(col_name, "")
            return "" if pd.isna(val) else str(val)
        return None

    def row_data(self, row: int) -> Optional[Dict[str, Any]]:
        if row < 0 or row >= len(self._rows):
            return None
        row_info = self._rows[row]
        if row_info.get("type") == "group":
            return None
        return row_info.get("data")

    def column_key(self, column: int) -> Optional[str]:
        if column < 0 or column >= len(self._columns):
            return None
        return self._columns[column][0]
