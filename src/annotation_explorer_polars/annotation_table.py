from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple

import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets


def _tag_color(tag: str) -> str:
    colors = {
        "A": "#c62828",  # AMR
        "B": "#1b5e20",  # BMRG
        "V": "#6a1b9a",  # Virulence
        "M": "#ff8f00",  # Mobile
        "P": "#1565c0",  # Phage
        "I": "#f9a825",  # Compound overlap
    }
    return colors.get(tag.upper(), "#000000")


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

        letters = [ch for ch in str(tags) if ch.strip()]
        if not letters:
            return

        painter.save()
        font = opt.font
        font.setBold(True)
        painter.setFont(font)
        metrics = QtGui.QFontMetrics(font)
        spacing = 6
        total_width = sum(metrics.horizontalAdvance(ch) for ch in letters)
        total_width += spacing * max(0, len(letters) - 1)
        x = opt.rect.x() + max(0, (opt.rect.width() - total_width) // 2)
        y = opt.rect.y() + (opt.rect.height() + metrics.ascent() - metrics.descent()) // 2
        for ch in letters:
            painter.setPen(QtGui.QColor(_tag_color(ch)))
            painter.drawText(QtCore.QPoint(x, y), ch.upper())
            x += metrics.horizontalAdvance(ch) + spacing
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
        opt.text = ""
        style = QtWidgets.QApplication.style() if opt.widget is None else opt.widget.style()
        style.drawControl(QtWidgets.QStyle.ControlElement.CE_ItemViewItem, opt, painter, opt.widget)

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

    def rowCount(self, _parent: QtCore.QModelIndex = QtCore.QModelIndex()) -> int:
        return len(self._rows)

    def columnCount(self, _parent: QtCore.QModelIndex = QtCore.QModelIndex()) -> int:
        return len(self._columns)

    def headerData(self, section: int, orientation: QtCore.Qt.Orientation, role: int = QtCore.Qt.ItemDataRole.DisplayRole):
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
                return ""
            if role == QtCore.Qt.ItemDataRole.BackgroundRole:
                return QtGui.QColor("#000000")
            if role == QtCore.Qt.ItemDataRole.SizeHintRole:
                return QtCore.QSize(0, 2)
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

        if role == QtCore.Qt.ItemDataRole.BackgroundRole and is_hit:
            return QtGui.QColor("#f5e7c6")

        if col_name == "Tags":
            if role == QtCore.Qt.ItemDataRole.UserRole:
                val = row_data.get(col_name, "")
                return "" if pd.isna(val) else str(val)
            if role == QtCore.Qt.ItemDataRole.DisplayRole:
                return ""
            if role == QtCore.Qt.ItemDataRole.TextAlignmentRole:
                return QtCore.Qt.AlignmentFlag.AlignCenter
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
