from __future__ import annotations

import hashlib
import html
import json
import re
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import polars as pl
from PyQt6 import QtCore, QtGui, QtWidgets

from .constants import PLASMID_TRUE_VALUES
from .panels.base import ANNOTATIONS, SAMPLE_TREE
from .session.requests import FocusPanel, params_tuple
from .model.query_engine import collect_sample_list
from .model.run_data import (
    clean_text,
    mobsuite_plasmid_group_key,
)
from .plasmid_model import plasmid_cluster_key_from_row


_TREE_SORT_ROLE = QtCore.Qt.ItemDataRole.UserRole.value + 1


def _make_mouse_transparent(widget: QtWidgets.QWidget) -> None:
    widget.setAttribute(QtCore.Qt.WidgetAttribute.WA_TransparentForMouseEvents, True)
    for child in widget.findChildren(QtWidgets.QWidget):
        child.setAttribute(QtCore.Qt.WidgetAttribute.WA_TransparentForMouseEvents, True)


class _TreeInlineLabel(QtWidgets.QLabel):
    def mousePressEvent(self, event: Optional[QtGui.QMouseEvent]) -> None:
        event.ignore() if event is not None else None

    def mouseDoubleClickEvent(self, event: Optional[QtGui.QMouseEvent]) -> None:
        event.ignore() if event is not None else None


class _SampleExplorerRow(QtWidgets.QFrame):
    def __init__(
        self,
        *,
        sample: str,
        country: str,
        flag_icon: QtGui.QIcon,
        quality_text: str,
        quality_color: str,
        quality_tooltip: str,
        completeness_text: str,
        completeness_color: str,
        completeness_tooltip: str,
        contamination_text: str,
        contamination_color: str,
        contamination_tooltip: str,
        row_tooltip: str,
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self._highlighted = False
        self.setObjectName("sampleExplorerRow")
        self.setToolTip(row_tooltip)
        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(3, 4, 8, 4)
        layout.setSpacing(4)

        sample_icon = _TreeInlineLabel()
        sample_icon.setPixmap(self._sample_pixmap())
        sample_icon.setFixedSize(18, 24)
        sample_icon.setToolTip(row_tooltip)
        layout.addWidget(sample_icon)

        name = _TreeInlineLabel(sample)
        name.setObjectName("sampleName")
        name.setToolTip(row_tooltip)
        name.setSizePolicy(QtWidgets.QSizePolicy.Policy.Maximum, QtWidgets.QSizePolicy.Policy.Preferred)
        name.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.NoTextInteraction)
        layout.addWidget(name, stretch=0)

        flag = _TreeInlineLabel()
        flag.setPixmap(flag_icon.pixmap(QtCore.QSize(30, 22)))
        flag.setToolTip(f"Country: {country}" if country else "Country not available")
        layout.addWidget(flag, stretch=0)

        layout.addWidget(self._quality_label(quality_text, quality_color, quality_tooltip), stretch=0)
        layout.addWidget(self._status_badge(completeness_text, completeness_color, completeness_tooltip, min_width=42), stretch=0)
        layout.addWidget(self._status_badge(contamination_text, contamination_color, contamination_tooltip, min_width=42), stretch=0)
        layout.addStretch(1)
        self._apply_style()
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_TransparentForMouseEvents, False)

    def _sample_pixmap(self) -> QtGui.QPixmap:
        pix = QtGui.QPixmap(18, 24)
        pix.fill(QtCore.Qt.GlobalColor.transparent)
        painter = QtGui.QPainter(pix)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        painter.setPen(QtGui.QPen(QtGui.QColor("#334155"), 1.2))
        painter.setBrush(QtGui.QBrush(QtGui.QColor("#e2e8f0")))
        painter.drawRoundedRect(QtCore.QRectF(3.0, 5.0, 12.0, 14.0), 2.8, 2.8)
        painter.setPen(QtGui.QPen(QtGui.QColor("#64748b"), 1.0))
        painter.drawLine(QtCore.QPointF(6.0, 9.0), QtCore.QPointF(12.0, 9.0))
        painter.drawLine(QtCore.QPointF(6.0, 12.0), QtCore.QPointF(12.0, 12.0))
        painter.drawLine(QtCore.QPointF(6.0, 15.0), QtCore.QPointF(10.5, 15.0))
        painter.setPen(QtGui.QPen(QtGui.QColor("#0f766e"), 1.3))
        painter.setBrush(QtGui.QBrush(QtGui.QColor("#99f6e4")))
        painter.drawEllipse(QtCore.QRectF(9.0, 2.5, 6.5, 6.5))
        painter.end()
        return pix

    def _quality_label(self, text: str, color: str, tooltip: str) -> QtWidgets.QLabel:
        label = _TreeInlineLabel(text)
        label.setToolTip(tooltip)
        label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        label.setMinimumWidth(34)
        label.setFixedHeight(26)
        label.setStyleSheet(
            f"color: {color}; background: transparent; border: none; padding: 0 2px; "
            "font-weight: 900; font-size: 16px;"
        )
        return label

    def _status_badge(self, text: str, color: str, tooltip: str, *, min_width: int) -> QtWidgets.QLabel:
        label = _TreeInlineLabel(text)
        label.setToolTip(tooltip)
        label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        label.setMinimumWidth(min_width)
        label.setFixedHeight(22)
        label.setObjectName("metricBadge")
        label.setProperty("badgeColor", color)
        label.setStyleSheet(
            "QLabel#metricBadge {"
            f"  background: {color};"
            "  color: white;"
            "  border-radius: 11px;"
            "  padding: 1px 7px;"
            "  font-weight: 800;"
            "  font-size: 11px;"
            "}"
        )
        return label

    def set_highlighted(self, active: bool) -> None:
        if self._highlighted == active:
            return
        self._highlighted = active
        self._apply_style()

    def _apply_style(self) -> None:
        bg = "transparent"
        border = "transparent"
        self.setStyleSheet(
            "QFrame#sampleExplorerRow {"
            f"  background: {bg};"
            f"  border: 1px solid {border};"
            "  border-radius: 12px;"
            "}"
            "QLabel#sampleName {"
            "  color: #162434;"
            "  font-size: 13px;"
            "  font-weight: 800;"
            "}"
        )


class _ContigExplorerRow(QtWidgets.QFrame):
    def __init__(
        self,
        *,
        label: str,
        icon: QtGui.QIcon,
        molecule_text: str,
        molecule_color: str,
        n_total: int,
        tooltip: str,
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self.setObjectName("contigExplorerRow")
        self.setToolTip(tooltip)
        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(8, 2, 8, 2)
        layout.setSpacing(8)

        icon_label = QtWidgets.QLabel()
        icon_label.setPixmap(icon.pixmap(QtCore.QSize(22, 22)))
        icon_label.setFixedSize(24, 24)
        layout.addWidget(icon_label)

        name = QtWidgets.QLabel(label)
        name.setObjectName("contigName")
        name.setMinimumWidth(0)
        name.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Preferred)
        layout.addWidget(name, stretch=1)

        molecule = QtWidgets.QLabel(molecule_text)
        molecule.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        molecule.setFixedSize(24, 22)
        molecule.setStyleSheet(
            f"background: {molecule_color}; color: white; border-radius: 11px; font-weight: 900; font-size: 11px;"
        )
        layout.addWidget(molecule)

        total = QtWidgets.QLabel(str(n_total))
        total.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
        total.setMinimumWidth(36)
        total.setStyleSheet("color: #64748b; font-weight: 700; font-size: 11px;")
        layout.addWidget(total)
        self.setStyleSheet(
            "QFrame#contigExplorerRow { background: transparent; border: none; }"
            "QLabel#contigName { color: #52677c; font-size: 12px; font-weight: 700; }"
        )
        _make_mouse_transparent(self)


class _SortableTreeItem(QtWidgets.QTreeWidgetItem):
    """Tree item that sorts numeric columns by a stored numeric key
    (``Qt.UserRole``) and otherwise by case-insensitive text."""

    def __lt__(self, other: "QtWidgets.QTreeWidgetItem") -> bool:
        tree = self.treeWidget()
        col = tree.sortColumn() if tree is not None else 0
        a = self.data(col, _TREE_SORT_ROLE)
        b = other.data(col, _TREE_SORT_ROLE)
        if a is not None and b is not None:
            if isinstance(a, (int, float)) and isinstance(b, (int, float)):
                return float(a) < float(b)
            return str(a).casefold() < str(b).casefold()
        return self.text(col).casefold() < other.text(col).casefold()


class TreeControllerMixin:
    def _tree_icon_cache(self) -> Dict[Tuple[str, str, str], QtGui.QIcon]:
        cache = getattr(self, "_tree_explorer_icon_cache", None)
        if not isinstance(cache, dict):
            cache = {}
            self._tree_explorer_icon_cache = cache
        return cache

    def _tree_vector_icon(
        self,
        kind: str,
        *,
        color: str = "#2563eb",
        text: str = "",
        tooltip_key: str = "",
    ) -> QtGui.QIcon:
        key = (str(kind), str(color), str(text or tooltip_key))
        cache = self._tree_icon_cache()
        cached = cache.get(key)
        if cached is not None:
            return cached

        size = 24
        dpr = QtWidgets.QApplication.primaryScreen().devicePixelRatio() if QtWidgets.QApplication.primaryScreen() else 1.0
        pix = QtGui.QPixmap(int(size * dpr), int(size * dpr))
        pix.setDevicePixelRatio(dpr)
        pix.fill(QtCore.Qt.GlobalColor.transparent)
        painter = QtGui.QPainter(pix)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)

        base = QtGui.QColor(color)
        dark = QtGui.QColor(base).darker(130)
        light = QtGui.QColor(base).lighter(170)

        if kind == "sample":
            painter.setPen(QtGui.QPen(QtGui.QColor("#1e3a8a"), 1.4))
            painter.setBrush(QtGui.QBrush(QtGui.QColor("#dbeafe")))
            painter.drawRoundedRect(QtCore.QRectF(7.0, 4.0, 10.0, 16.0), 4.0, 4.0)
            painter.setPen(QtGui.QPen(QtGui.QColor("#1d4ed8"), 1.2))
            painter.setBrush(QtGui.QBrush(QtGui.QColor("#60a5fa")))
            painter.drawRoundedRect(QtCore.QRectF(8.3, 3.0, 7.4, 4.0), 1.6, 1.6)
            painter.setPen(QtGui.QPen(QtGui.QColor("#93c5fd"), 1.0))
            painter.drawLine(QtCore.QPointF(10.0, 9.0), QtCore.QPointF(10.0, 16.0))
        elif kind == "card":
            painter.setPen(QtGui.QPen(QtGui.QColor("#64748b"), 1.2))
            painter.setBrush(QtGui.QBrush(QtGui.QColor("#ffffff")))
            painter.drawRoundedRect(QtCore.QRectF(5.0, 4.0, 14.0, 16.0), 2.4, 2.4)
            painter.setPen(QtGui.QPen(QtGui.QColor("#93a4b7"), 1.1))
            painter.drawLine(QtCore.QPointF(8.0, 9.0), QtCore.QPointF(16.0, 9.0))
            painter.drawLine(QtCore.QPointF(8.0, 13.0), QtCore.QPointF(15.0, 13.0))
            painter.drawLine(QtCore.QPointF(8.0, 17.0), QtCore.QPointF(13.0, 17.0))
        elif kind == "plasmid":
            painter.setPen(QtGui.QPen(dark, 2.2))
            painter.setBrush(QtGui.QBrush(QtCore.Qt.BrushStyle.NoBrush))
            painter.drawEllipse(QtCore.QRectF(5.0, 5.0, 14.0, 14.0))
            painter.setPen(QtGui.QPen(base.lighter(115), 1.6))
            painter.drawArc(QtCore.QRectF(5.0, 5.0, 14.0, 14.0), 20 * 16, 110 * 16)
        elif kind == "chromosome":
            painter.setPen(QtGui.QPen(dark, 1.3))
            painter.setBrush(QtGui.QBrush(light))
            painter.drawRoundedRect(QtCore.QRectF(4.0, 8.0, 16.0, 8.0), 4.0, 4.0)
            painter.setPen(QtGui.QPen(base.darker(120), 1.1))
            painter.drawLine(QtCore.QPointF(8.0, 10.0), QtCore.QPointF(16.0, 14.0))
            painter.drawLine(QtCore.QPointF(8.0, 14.0), QtCore.QPointF(16.0, 10.0))
        elif kind == "country":
            painter.setPen(QtGui.QPen(dark, 1.0))
            painter.setBrush(QtGui.QBrush(base))
            painter.drawRoundedRect(QtCore.QRectF(2.5, 4.0, 19.0, 16.0), 5.5, 5.5)
            font = QtGui.QFont()
            font.setPointSize(7)
            font.setBold(True)
            painter.setFont(font)
            painter.setPen(QtGui.QColor("#ffffff"))
            painter.drawText(QtCore.QRectF(2.5, 4.0, 19.0, 16.0), int(QtCore.Qt.AlignmentFlag.AlignCenter), text[:3].upper() or "?")
        elif kind == "genome":
            painter.setPen(QtGui.QPen(QtGui.QColor("#0f766e"), 1.6))
            painter.setBrush(QtGui.QBrush(QtGui.QColor("#5eead4")))
            painter.drawEllipse(QtCore.QRectF(6.0, 6.0, 12.0, 12.0))
            painter.setPen(QtGui.QPen(QtGui.QColor("#0f766e"), 1.2))
            painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
            painter.drawArc(QtCore.QRectF(8.4, 8.4, 7.2, 7.2), 40 * 16, 150 * 16)
        else:
            painter.setPen(QtGui.QPen(dark, 1.0))
            painter.setBrush(QtGui.QBrush(light))
            painter.drawRoundedRect(QtCore.QRectF(4.0, 4.0, 16.0, 16.0), 4.0, 4.0)

        painter.end()
        icon = QtGui.QIcon(pix)
        cache[key] = icon
        return icon

    def _country_flag_icon(self, country: str) -> QtGui.QIcon:
        code = self._country_badge_text(country)
        if not code:
            return self._tree_vector_icon("country", color="#94a3b8", text="?")

        key = ("flag", code, str(country or ""))
        cache = self._tree_icon_cache()
        cached = cache.get(key)
        if cached is not None:
            return cached

        width, height = 30, 22
        dpr = QtWidgets.QApplication.primaryScreen().devicePixelRatio() if QtWidgets.QApplication.primaryScreen() else 1.0
        pix = QtGui.QPixmap(int(width * dpr), int(height * dpr))
        pix.setDevicePixelRatio(dpr)
        pix.fill(QtCore.Qt.GlobalColor.transparent)

        painter = QtGui.QPainter(pix)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        rect = QtCore.QRectF(3.0, 4.0, 24.0, 15.0)
        clip = QtGui.QPainterPath()
        clip.addRoundedRect(rect, 2.4, 2.4)
        painter.save()
        painter.setClipPath(clip)
        self._draw_country_flag(painter, rect, code)
        painter.restore()
        painter.setPen(QtGui.QPen(QtGui.QColor(15, 23, 42, 72), 1.0))
        painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
        painter.drawRoundedRect(rect, 2.4, 2.4)
        painter.end()

        icon = QtGui.QIcon(pix)
        cache[key] = icon
        return icon

    def _country_badge_icon(self, country: str) -> QtGui.QIcon:
        return self._country_flag_icon(country)

    @staticmethod
    def _fill_flag_rect(painter: QtGui.QPainter, rect: QtCore.QRectF, color: str) -> None:
        painter.fillRect(rect, QtGui.QColor(color))

    @staticmethod
    def _draw_horizontal_stripes(painter: QtGui.QPainter, rect: QtCore.QRectF, colors: Sequence[str], ratios: Optional[Sequence[float]] = None) -> None:
        if not colors:
            return
        if ratios is None or len(ratios) != len(colors):
            ratios = [1.0] * len(colors)
        total = float(sum(max(0.0, r) for r in ratios)) or 1.0
        y = rect.top()
        for idx, (color, ratio) in enumerate(zip(colors, ratios)):
            h = rect.height() * (float(ratio) / total)
            if idx == len(colors) - 1:
                h = rect.bottom() - y
            painter.fillRect(QtCore.QRectF(rect.left(), y, rect.width(), h), QtGui.QColor(color))
            y += h

    @staticmethod
    def _draw_vertical_stripes(painter: QtGui.QPainter, rect: QtCore.QRectF, colors: Sequence[str], ratios: Optional[Sequence[float]] = None) -> None:
        if not colors:
            return
        if ratios is None or len(ratios) != len(colors):
            ratios = [1.0] * len(colors)
        total = float(sum(max(0.0, r) for r in ratios)) or 1.0
        x = rect.left()
        for idx, (color, ratio) in enumerate(zip(colors, ratios)):
            w = rect.width() * (float(ratio) / total)
            if idx == len(colors) - 1:
                w = rect.right() - x
            painter.fillRect(QtCore.QRectF(x, rect.top(), w, rect.height()), QtGui.QColor(color))
            x += w

    def _draw_country_flag(self, painter: QtGui.QPainter, rect: QtCore.QRectF, code: str) -> None:
        code = str(code or "").strip().upper()

        vertical_specs = {
            "BE": ("#000000", "#ffd90c", "#ef3340"),
            "FR": ("#0055a4", "#ffffff", "#ef4135"),
            "IE": ("#169b62", "#ffffff", "#ff883e"),
            "IT": ("#009246", "#ffffff", "#ce2b37"),
            "MX": ("#006847", "#ffffff", "#ce1126"),
            "RO": ("#002b7f", "#fcd116", "#ce1126"),
        }
        horizontal_specs = {
            "AR": ("#75aadb", "#ffffff", "#75aadb"),
            "AT": ("#ed2939", "#ffffff", "#ed2939"),
            "BG": ("#ffffff", "#00966e", "#d62612"),
            "DE": ("#000000", "#dd0000", "#ffce00"),
            "EE": ("#4891d9", "#000000", "#ffffff"),
            "ES": ("#aa151b", "#f1bf00", "#aa151b"),
            "HU": ("#ce2939", "#ffffff", "#477050"),
            "LT": ("#fdb913", "#006a44", "#c1272d"),
            "LU": ("#ef3340", "#ffffff", "#00a3e0"),
            "NL": ("#ae1c28", "#ffffff", "#21468b"),
            "PL": ("#ffffff", "#dc143c"),
            "RU": ("#ffffff", "#0039a6", "#d52b1e"),
            "UA": ("#005bbb", "#ffd500"),
        }

        if code in vertical_specs:
            self._draw_vertical_stripes(painter, rect, vertical_specs[code])
            return
        if code in horizontal_specs:
            ratios = (1.0, 2.0, 1.0) if code == "ES" else None
            self._draw_horizontal_stripes(painter, rect, horizontal_specs[code], ratios=ratios)
            return

        if code == "CH":
            self._fill_flag_rect(painter, rect, "#da291c")
            painter.fillRect(QtCore.QRectF(rect.center().x() - 2.0, rect.top() + 3.0, 4.0, rect.height() - 6.0), QtGui.QColor("#ffffff"))
            painter.fillRect(QtCore.QRectF(rect.left() + 6.0, rect.center().y() - 2.0, rect.width() - 12.0, 4.0), QtGui.QColor("#ffffff"))
            return
        if code == "CN":
            self._fill_flag_rect(painter, rect, "#de2910")
            painter.setPen(QtCore.Qt.PenStyle.NoPen)
            painter.setBrush(QtGui.QColor("#ffde00"))
            painter.drawEllipse(QtCore.QRectF(rect.left() + 3.2, rect.top() + 3.0, 3.2, 3.2))
            return
        if code == "JP":
            self._fill_flag_rect(painter, rect, "#ffffff")
            painter.setPen(QtCore.Qt.PenStyle.NoPen)
            painter.setBrush(QtGui.QColor("#bc002d"))
            painter.drawEllipse(QtCore.QRectF(rect.center().x() - 4.0, rect.center().y() - 4.0, 8.0, 8.0))
            return
        if code == "IN":
            self._draw_horizontal_stripes(painter, rect, ("#ff9933", "#ffffff", "#138808"))
            painter.setPen(QtGui.QPen(QtGui.QColor("#000080"), 0.9))
            painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
            painter.drawEllipse(QtCore.QRectF(rect.center().x() - 2.1, rect.center().y() - 2.1, 4.2, 4.2))
            return
        if code == "US":
            stripe_h = rect.height() / 13.0
            for idx in range(13):
                color = "#b22234" if idx % 2 == 0 else "#ffffff"
                painter.fillRect(QtCore.QRectF(rect.left(), rect.top() + idx * stripe_h, rect.width(), stripe_h + 0.2), QtGui.QColor(color))
            painter.fillRect(QtCore.QRectF(rect.left(), rect.top(), rect.width() * 0.42, stripe_h * 7.0), QtGui.QColor("#3c3b6e"))
            painter.setPen(QtCore.Qt.PenStyle.NoPen)
            painter.setBrush(QtGui.QColor("#ffffff"))
            for row in range(3):
                for col in range(3):
                    painter.drawEllipse(QtCore.QRectF(rect.left() + 2.0 + col * 3.0, rect.top() + 1.6 + row * 2.1, 0.8, 0.8))
            return
        if code == "GB" or code == "UK":
            self._fill_flag_rect(painter, rect, "#012169")
            painter.setPen(QtGui.QPen(QtGui.QColor("#ffffff"), 4.0))
            painter.drawLine(rect.topLeft(), rect.bottomRight())
            painter.drawLine(rect.bottomLeft(), rect.topRight())
            painter.setPen(QtGui.QPen(QtGui.QColor("#c8102e"), 1.8))
            painter.drawLine(rect.topLeft(), rect.bottomRight())
            painter.drawLine(rect.bottomLeft(), rect.topRight())
            painter.setPen(QtGui.QPen(QtGui.QColor("#ffffff"), 5.0))
            painter.drawLine(QtCore.QPointF(rect.center().x(), rect.top()), QtCore.QPointF(rect.center().x(), rect.bottom()))
            painter.drawLine(QtCore.QPointF(rect.left(), rect.center().y()), QtCore.QPointF(rect.right(), rect.center().y()))
            painter.setPen(QtGui.QPen(QtGui.QColor("#c8102e"), 2.4))
            painter.drawLine(QtCore.QPointF(rect.center().x(), rect.top()), QtCore.QPointF(rect.center().x(), rect.bottom()))
            painter.drawLine(QtCore.QPointF(rect.left(), rect.center().y()), QtCore.QPointF(rect.right(), rect.center().y()))
            return
        if code == "CA":
            self._draw_vertical_stripes(painter, rect, ("#ff0000", "#ffffff", "#ff0000"), ratios=(1.0, 2.0, 1.0))
            painter.setPen(QtCore.Qt.PenStyle.NoPen)
            painter.setBrush(QtGui.QColor("#ff0000"))
            painter.drawPolygon(
                QtGui.QPolygonF(
                    [
                        QtCore.QPointF(rect.center().x(), rect.top() + 3.0),
                        QtCore.QPointF(rect.center().x() + 3.0, rect.center().y()),
                        QtCore.QPointF(rect.center().x() + 1.0, rect.bottom() - 3.0),
                        QtCore.QPointF(rect.center().x(), rect.bottom() - 5.0),
                        QtCore.QPointF(rect.center().x() - 1.0, rect.bottom() - 3.0),
                        QtCore.QPointF(rect.center().x() - 3.0, rect.center().y()),
                    ]
                )
            )
            return
        if code == "BR":
            self._fill_flag_rect(painter, rect, "#009b3a")
            painter.setPen(QtCore.Qt.PenStyle.NoPen)
            painter.setBrush(QtGui.QColor("#ffdf00"))
            painter.drawPolygon(
                QtGui.QPolygonF(
                    [
                        QtCore.QPointF(rect.center().x(), rect.top() + 2.0),
                        QtCore.QPointF(rect.right() - 3.0, rect.center().y()),
                        QtCore.QPointF(rect.center().x(), rect.bottom() - 2.0),
                        QtCore.QPointF(rect.left() + 3.0, rect.center().y()),
                    ]
                )
            )
            painter.setBrush(QtGui.QColor("#002776"))
            painter.drawEllipse(QtCore.QRectF(rect.center().x() - 3.4, rect.center().y() - 3.4, 6.8, 6.8))
            return
        if code == "TR":
            self._fill_flag_rect(painter, rect, "#e30a17")
            painter.setPen(QtCore.Qt.PenStyle.NoPen)
            painter.setBrush(QtGui.QColor("#ffffff"))
            painter.drawEllipse(QtCore.QRectF(rect.left() + 5.0, rect.center().y() - 4.0, 8.0, 8.0))
            painter.setBrush(QtGui.QColor("#e30a17"))
            painter.drawEllipse(QtCore.QRectF(rect.left() + 7.0, rect.center().y() - 3.2, 6.4, 6.4))
            painter.setBrush(QtGui.QColor("#ffffff"))
            painter.drawEllipse(QtCore.QRectF(rect.left() + 15.0, rect.center().y() - 1.6, 3.2, 3.2))
            return

        seed = int(hashlib.sha1(code.encode("utf-8")).hexdigest()[:8], 16)
        colors = [
            QtGui.QColor.fromHsv(seed % 360, 150, 190).name(),
            QtGui.QColor.fromHsv((seed >> 8) % 360, 88, 238).name(),
            QtGui.QColor.fromHsv((seed >> 16) % 360, 150, 190).name(),
        ]
        self._draw_horizontal_stripes(painter, rect, colors)

    @staticmethod
    def _country_badge_text(country: str) -> str:
        text = str(country or "").strip()
        if not text:
            return ""
        if re.fullmatch(r"[A-Za-z]{2,3}", text):
            return text.upper()
        normalized = re.sub(r"[^a-z0-9]+", " ", text.lower()).strip()
        aliases = {
            "usa": "US",
            "u s a": "US",
            "united states": "US",
            "united states of america": "US",
            "uk": "UK",
            "u k": "UK",
            "united kingdom": "UK",
            "great britain": "GB",
            "england": "GB",
            "scotland": "GB",
            "romania": "RO",
            "france": "FR",
            "germany": "DE",
            "italy": "IT",
            "spain": "ES",
            "portugal": "PT",
            "netherlands": "NL",
            "belgium": "BE",
            "switzerland": "CH",
            "austria": "AT",
            "poland": "PL",
            "greece": "GR",
            "turkey": "TR",
            "china": "CN",
            "india": "IN",
            "japan": "JP",
            "south korea": "KR",
            "korea": "KR",
            "australia": "AU",
            "new zealand": "NZ",
            "canada": "CA",
            "brazil": "BR",
            "argentina": "AR",
            "mexico": "MX",
            "south africa": "ZA",
        }
        if normalized in aliases:
            return aliases[normalized]
        words = [w for w in normalized.split() if w not in {"the", "of", "and", "republic"}]
        if len(words) >= 2:
            return (words[0][0] + words[1][0]).upper()
        compact = re.sub(r"[^A-Za-z0-9]+", "", text)
        return compact[:2].upper() if compact else ""

    def _country_from_metadata(self, meta_vals: Dict[str, Any]) -> str:
        if not isinstance(meta_vals, dict):
            return ""
        preferred: List[Tuple[int, str, Any]] = []
        for col, value in meta_vals.items():
            col_text = str(col or "").strip()
            norm_col = re.sub(r"[^a-z0-9]+", "", col_text.lower())
            if not norm_col or norm_col == "sample":
                continue
            priority = 99
            if norm_col in {"country", "origincountry", "isolationcountry", "countryregion", "nation"}:
                priority = 0
            elif "country" in norm_col:
                priority = 1
            elif norm_col in {"geolocname", "geolocation", "location", "origin"}:
                priority = 2
            elif "geoloc" in norm_col or "location" in norm_col:
                priority = 3
            if priority < 99:
                preferred.append((priority, col_text, value))
        for _priority, _col, value in sorted(preferred, key=lambda item: item[0]):
            text = clean_text(value)
            if not text:
                continue
            text = re.split(r"[:;,/|]", text, maxsplit=1)[0].strip()
            if text and text.lower() not in {"missing", "unknown", "not provided", "na", "n/a"}:
                return text
        return ""

    @staticmethod
    def _tree_metric_float(value: Any) -> Optional[float]:
        if value is None:
            return None
        try:
            if pd.isna(value):
                return None
        except Exception:
            pass
        text = str(value).strip().replace(",", "")
        if not text or text.lower() in {"nan", "none", "null", "<na>", "na", "n/a", "-", "."}:
            return None
        try:
            val = float(text)
        except Exception:
            return None
        return val if np.isfinite(val) else None

    @staticmethod
    def _tree_metric_int_text(value: Any) -> str:
        val = TreeControllerMixin._tree_metric_float(value)
        if val is None:
            return "--"
        return f"{int(round(val)):,}"

    @staticmethod
    def _tree_metric_float_text(value: Any, digits: int = 1, suffix: str = "") -> str:
        val = TreeControllerMixin._tree_metric_float(value)
        if val is None:
            return "--"
        return f"{val:.{digits}f}{suffix}"

    @staticmethod
    def _quality_badge_style(value: Any) -> Tuple[str, str]:
        score = TreeControllerMixin._tree_metric_float(value)
        if score is None:
            return "--", "#94a3b8"
        if score >= 90:
            color = "#059669"
        elif score >= 75:
            color = "#d97706"
        else:
            color = "#dc2626"
        return f"{int(round(score))}", color

    @staticmethod
    def _completeness_badge_style(value: Any) -> Tuple[str, str]:
        completeness = TreeControllerMixin._tree_metric_float(value)
        if completeness is None:
            return "CMP ?", "#94a3b8"
        if completeness >= 95.0:
            return "CMP", "#059669"
        if completeness >= 90.0:
            return "CMP", "#d97706"
        return "INC", "#dc2626"

    @staticmethod
    def _contamination_badge_style(value: Any) -> Tuple[str, str]:
        contamination = TreeControllerMixin._tree_metric_float(value)
        if contamination is None:
            return "CNT ?", "#94a3b8"
        if contamination <= 2.0:
            return "CNT", "#059669"
        if contamination <= 5.0:
            return "CNT", "#d97706"
        return "CONT", "#dc2626"

    @staticmethod
    def _completeness_cell(value: Any) -> Tuple[str, str]:
        """Completeness as a coloured percentage for the table column."""
        completeness = TreeControllerMixin._tree_metric_float(value)
        if completeness is None:
            return "--", "#94a3b8"
        if completeness >= 95.0:
            return f"{completeness:.1f}%", "#059669"
        if completeness >= 90.0:
            return f"{completeness:.1f}%", "#d97706"
        return f"{completeness:.1f}%", "#dc2626"

    @staticmethod
    def _contamination_cell(value: Any) -> Tuple[str, str]:
        contamination = TreeControllerMixin._tree_metric_float(value)
        if contamination is None:
            return "--", "#94a3b8"
        if contamination <= 2.0:
            return f"{contamination:.1f}%", "#059669"
        if contamination <= 5.0:
            return f"{contamination:.1f}%", "#d97706"
        return f"{contamination:.1f}%", "#dc2626"

    @staticmethod
    def _st_label(value: Any) -> str:
        subtype = TreeControllerMixin._tree_metric_float(value)
        if subtype is None:
            return "--"
        return f"ST{int(round(subtype))}"

    def _sample_qc_tooltip(self, sample: str, country: str, meta_vals: Dict[str, Any]) -> str:
        lines = [f"Sample: {sample}"]
        if country:
            lines.append(f"Country: {country}")
        return "\n".join(lines)

    @staticmethod
    def _format_qc_tooltip(
        title: str,
        rows: Sequence[Tuple[str, str]],
        *,
        detail: str = "",
        notes: str = "",
    ) -> str:
        body_rows = []
        for label, value in rows:
            body_rows.append(
                "<tr>"
                f"<td style='padding:2px 12px 2px 0; color:#e5edf7;'>{html.escape(str(label))}</td>"
                f"<td style='padding:2px 0; color:#ffffff; font-weight:700; text-align:right;'>{html.escape(str(value))}</td>"
                "</tr>"
            )
        detail_html = ""
        if detail:
            detail_html = (
                "<div style='margin-top:8px; padding-top:7px; border-top:1px solid rgba(255,255,255,0.18); "
                "color:#e5edf7; line-height:1.35;'>"
                f"{html.escape(detail)}"
                "</div>"
            )
        notes_html = ""
        if notes:
            notes_html = (
                "<div style='margin-top:8px; padding-top:7px; border-top:1px solid rgba(255,255,255,0.18); "
                "color:#e5edf7; line-height:1.35;'>"
                f"<b style='color:#ffffff;'>Notes:</b> {html.escape(notes)}"
                "</div>"
            )
        return (
            "<html><body>"
            "<div style='width:320px; color:#ffffff; font-size:12px;'>"
            f"<div style='font-weight:800; font-size:13px; margin-bottom:6px; color:#ffffff;'>{html.escape(title)}</div>"
            "<table cellspacing='0' cellpadding='0' style='width:100%;'>"
            f"{''.join(body_rows)}"
            "</table>"
            f"{detail_html}"
            f"{notes_html}"
            "</div>"
            "</body></html>"
        )

    def _quast_tooltip(self, meta_vals: Dict[str, Any]) -> str:
        quality = self._tree_metric_int_text(meta_vals.get("Quality"))
        notes = clean_text(meta_vals.get("QualityNotes", ""))
        rows = [
            ("Quality score", f"{quality}/100" if quality != "--" else "--"),
            ("Total length", f"{self._tree_metric_int_text(meta_vals.get('QuastTotalLength'))} bp"),
            ("Contigs", self._tree_metric_int_text(meta_vals.get("QuastContigs"))),
            ("Largest contig", f"{self._tree_metric_int_text(meta_vals.get('QuastLargestContig'))} bp"),
            ("N50 / L50", f"{self._tree_metric_int_text(meta_vals.get('QuastN50'))} / {self._tree_metric_int_text(meta_vals.get('QuastL50'))}"),
            ("GC", self._tree_metric_float_text(meta_vals.get("QuastGC"), 2, "%")),
            ("N bases", self._tree_metric_int_text(meta_vals.get("QuastNs"))),
            ("N per 100 kbp", self._tree_metric_float_text(meta_vals.get("QuastNsPer100kbp"), 2)),
        ]
        quast_notes = [note for note in notes.split("; ") if "CheckM2" not in note]
        return self._format_qc_tooltip("QUAST assembly metrics", rows, notes="; ".join(quast_notes))

    def _checkm2_tooltip(self, meta_vals: Dict[str, Any], *, focus: str = "") -> str:
        completeness = self._tree_metric_float_text(meta_vals.get("CheckM2Completeness"), 2, "%")
        contamination = self._tree_metric_float_text(meta_vals.get("CheckM2Contamination"), 2, "%")
        rows = [
            ("Completeness", completeness),
            ("Contamination", contamination),
        ]
        model = clean_text(meta_vals.get("CheckM2CompletenessModel", ""))
        if model:
            rows.append(("Completeness model", model))
        detail = ""
        if focus == "completeness":
            detail = "Badge: incomplete if completeness is below 90%."
        elif focus == "contamination":
            detail = "Badge: contaminated if contamination is above 5%."
        notes = clean_text(meta_vals.get("QualityNotes", ""))
        return self._format_qc_tooltip("CheckM2 genome quality", rows, detail=detail, notes=notes)

    def _contig_tooltip(
        self,
        *,
        sample: str,
        display_label: str,
        molecule_label: str,
        n_total: int,
        member_contigs: Sequence[str],
        contig: Dict[str, Any],
    ) -> str:
        lines = [
            f"{display_label}",
            f"Sample: {sample}",
            f"Type: {'Plasmid' if molecule_label == 'P' else 'Chromosome'}",
            f"Annotated genes: {n_total:,}",
        ]
        members = [str(member).strip() for member in member_contigs if str(member).strip()]
        if len(members) > 1:
            lines.append(f"Member contigs: {', '.join(members)}")
        for label, field in (
            ("Replicon", "MobSuiteRepliconTypes"),
            ("Relaxase", "MobSuiteRelaxaseTypes"),
            ("MPF", "MobSuiteMPFType"),
            ("oriT", "MobSuiteOriTTypes"),
            ("Mobility", "MobSuitePredictedMobility"),
        ):
            text = clean_text(contig.get(field, ""))
            if text:
                lines.append(f"{label}: {text}")
        return "\n".join(lines)

    def _rerender_visible_table_columns(self) -> None:
        tree = getattr(self, "tree", None)
        if tree is None:
            return
        fixed = set(getattr(self, "_fixed_tree_columns", (0,)))
        # Only the fixed columns show; the compact view also drops the header.
        compact = bool(getattr(self, "_sample_tree_compact_view", False))
        for col_idx in range(1 if compact else 0, tree.columnCount()):
            tree.setColumnHidden(col_idx, col_idx not in fixed)
        if compact:
            tree.setHeaderHidden(True)

    @staticmethod
    def _truthy_text(value: Any) -> bool:
        text = str(value or "").strip().lower()
        return text in PLASMID_TRUE_VALUES

    @staticmethod
    def _nonempty_text(value: Any) -> str:
        text = str(value or "").strip()
        return "" if text.lower() in {"", "nan", "none", "null", "<na>"} else text

    def _mobsuite_plasmid_group_key(self, row: Dict[str, Any]) -> str:
        return mobsuite_plasmid_group_key(row)

    def _display_contig_records(
        self,
        records: Sequence[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        grouped: Dict[str, Dict[str, Any]] = {}
        ordered: List[Dict[str, Any]] = []
        for raw in records:
            row = dict(raw)
            sample = self._nonempty_text(row.get("Sample"))
            contig = self._nonempty_text(row.get("Contig"))
            group_key = self._mobsuite_plasmid_group_key(row)
            if not sample or not contig or not group_key:
                row["__member_contigs"] = [contig] if contig else []
                row["__display_label"] = f"Contig {contig}"
                row["__is_plasmid_group"] = False
                ordered.append(row)
                continue

            group = grouped.get(group_key)
            if group is None:
                members = list(self.run.mobsuite_plasmid_members.get((sample, group_key), []))
                if not members:
                    members = [contig]
                row["__member_contigs"] = members
                row["__visible_member_contigs"] = [contig]
                plasmid_index = self.run.plasmid_cluster_index
                if plasmid_index is not None and hasattr(plasmid_index, "display_name_for_row"):
                    row["__display_label"] = plasmid_index.display_name_for_row(row)
                else:
                    row["__display_label"] = f"Plasmid {group_key}"
                row["__is_plasmid_group"] = True
                row["__plasmid_group_key"] = group_key
                grouped[group_key] = row
                ordered.append(row)
            else:
                visible_members = group.setdefault("__visible_member_contigs", [])
                if contig not in visible_members:
                    visible_members.append(contig)
                for field in (
                    "MobSuiteRepliconTypes",
                    "MobSuiteRelaxaseTypes",
                    "MobSuiteMPFType",
                    "MobSuiteOriTTypes",
                    "MobSuitePredictedMobility",
                    "MobSuiteMashNearestNeighbor",
                    "MobSuiteMashNeighborDistance",
                ):
                    if not self._nonempty_text(group.get(field)):
                        group[field] = row.get(field)
        return ordered

    @staticmethod
    def _tree_item_payload(item: Optional[QtWidgets.QTreeWidgetItem]) -> Dict[str, Any]:
        if item is None:
            return {}
        payload_raw = item.data(0, QtCore.Qt.ItemDataRole.UserRole)
        if not payload_raw:
            return {}
        try:
            payload = json.loads(payload_raw)
        except (TypeError, ValueError):
            return {}
        return payload if isinstance(payload, dict) else {}

    def _tree_item_state_key(
        self,
        item: Optional[QtWidgets.QTreeWidgetItem],
    ) -> Optional[Tuple[object, ...]]:
        payload = self._tree_item_payload(item)
        sample = str(payload.get("sample") or "").strip()
        if not sample:
            return None
        contig = str(payload.get("contig") or "").strip()
        if not contig:
            return ("sample", sample)
        members = tuple(
            str(member).strip()
            for member in (payload.get("contigs") or [])
            if str(member).strip()
        )
        return ("contig", sample, contig, members)

    def _snapshot_metadata_tree_state(self) -> Dict[str, object]:
        expanded_samples: set[str] = set()
        for index in range(self.tree.topLevelItemCount()):
            item = self.tree.topLevelItem(index)
            if item is None or not item.isExpanded():
                continue
            payload = self._tree_item_payload(item)
            sample = str(payload.get("sample") or "").strip()
            if sample:
                expanded_samples.add(sample)
        return {
            "expanded_samples": expanded_samples,
            "current_key": self._tree_item_state_key(self.tree.currentItem()),
            "vertical_scroll": self.tree.verticalScrollBar().value(),
            "horizontal_scroll": self.tree.horizontalScrollBar().value(),
        }

    def _restore_metadata_tree_state(self, state: Dict[str, object]) -> None:
        expanded_samples = {
            str(sample)
            for sample in state.get("expanded_samples", set())
            if str(sample).strip()
        }
        current_key = state.get("current_key")
        current_item: Optional[QtWidgets.QTreeWidgetItem] = None
        blocker = QtCore.QSignalBlocker(self.tree)
        try:
            for sample, sample_item in self._sample_tree_items.items():
                sample_item.setExpanded(sample in expanded_samples)
                if current_key == ("sample", sample):
                    current_item = sample_item
                if not isinstance(current_key, tuple) or len(current_key) < 3:
                    continue
                if current_key[0] != "contig" or current_key[1] != sample:
                    continue
                for child_index in range(sample_item.childCount()):
                    child = sample_item.child(child_index)
                    if self._tree_item_state_key(child) == current_key:
                        current_item = child
                        break
            if current_item is not None:
                self.tree.setCurrentItem(current_item)
        finally:
            del blocker

        self.tree.doItemsLayout()
        self.tree.verticalScrollBar().setValue(int(state.get("vertical_scroll", 0)))
        self.tree.horizontalScrollBar().setValue(int(state.get("horizontal_scroll", 0)))

    def _ensure_sample_tree_sort_sync(self) -> None:
        if bool(getattr(self, "_sample_tree_sort_sync_connected", False)):
            return
        self._sample_tree_sort_sync_connected = True
        self._sample_tree_header_sort_override = None
        self._sample_tree_programmatic_sort_change = False
        self.tree.header().sortIndicatorChanged.connect(self._on_sample_tree_header_sort_changed)
        for combo_name in ("sample_sort_column_box", "sample_sort_order_box"):
            combo = getattr(self, combo_name, None)
            if isinstance(combo, QtWidgets.QComboBox):
                combo.currentIndexChanged.connect(self._on_sample_tree_sidebar_sort_changed)

    def _on_sample_tree_header_sort_changed(
        self,
        column: int,
        order: QtCore.Qt.SortOrder,
    ) -> None:
        if bool(getattr(self, "_sample_tree_programmatic_sort_change", False)):
            return
        self._sample_tree_header_sort_override = (int(column), order)

    def _on_sample_tree_sidebar_sort_changed(self, _index: int) -> None:
        self._sample_tree_header_sort_override = None

    def _sample_tree_sidebar_sort_column(self, sort_column: str) -> int:
        if sort_column == "__sample__":
            return 0
        try:
            return 5 + self.meta_columns_for_tree.index(sort_column)
        except ValueError:
            return 0

    def _apply_sample_tree_sort_indicator(self) -> None:
        if getattr(self, "_sample_tree_header_sort_override", None) is not None:
            return
        sort_column, ascending = self._sample_sort_settings()
        tree_column = self._sample_tree_sidebar_sort_column(sort_column)
        order = (
            QtCore.Qt.SortOrder.AscendingOrder
            if ascending
            else QtCore.Qt.SortOrder.DescendingOrder
        )
        self._sample_tree_programmatic_sort_change = True
        try:
            self.tree.header().setSortIndicator(tree_column, order)
        finally:
            self._sample_tree_programmatic_sort_change = False

    def _sample_tree_search_query(self) -> str:
        """Search text the tree's counts were computed for, or "" when none.

        This is the applied gene and/or product search (the boxes may already
        be ahead of the result on screen while the user types). While it is
        non-empty, matching samples open to show their contigs.
        """
        applied = self._applied_filters()
        return " / ".join(text for text in (applied.gene, applied.product) if text)

    def _render_sample_tree(self) -> None:
        """Rebuild the sample tree from the applied result (no compute). Its
        rows are the window's one sample set (plan 4.4, D4)."""
        result = self.ctx.result
        if result is None:
            return
        self._populate_metadata_tree(
            result.metadata,
            result.contig_records,
            result.filtered_counts,
            result.samples,
        )
        self._rerender_visible_table_columns()

    def _resort_sample_tree(self, *_args: object) -> None:
        """The rail's sort boxes changed: sorting is a tree option, not a
        filter, so the tree re-sorts from the result it has (plan 4.4)."""
        self._sample_tree_header_sort_override = None
        self._render_sample_tree()

    @staticmethod
    def _copies_tooltip_line(count: int, query: str, scope: str = "sample") -> str:
        noun = "copy" if int(count) == 1 else "copies"
        return f"{int(count)} {noun} of {query} in this {scope}"

    def _search_expanded_samples(
        self,
        search_query: str,
        previous_state: Dict[str, object],
        sample_names: Sequence[str],
        filtered_counts: Dict[str, int],
    ) -> set[str]:
        """Samples to show expanded after a rebuild, given the search text.

        - search starts: remember the user's expansion, open every match;
        - search text changes: open every match again (zero-match samples closed);
        - same search text: every match is open except the ones the user
          collapsed by hand since the text last changed; those stay closed even
          if another filter briefly hides them and they match again;
        - search cleared: give back the expansion from before the search.
        Without a search this is just the previous expansion.
        Matching is case-insensitive, so a case-only edit is not a new search.
        """
        previous_expanded = set(previous_state.get("expanded_samples", set()))
        previous_query = str(getattr(self, "_sample_tree_expand_query", "") or "")
        previous_matching = set(getattr(self, "_sample_tree_expand_matching", set()) or set())
        collapsed = set(getattr(self, "_sample_tree_search_collapsed", set()) or set())
        matching = {s for s in sample_names if int(filtered_counts.get(s, 0)) > 0}
        if search_query:
            if not previous_query:
                self._sample_tree_pre_search_expanded = previous_expanded
                collapsed = set()
            elif search_query.casefold() != previous_query.casefold():
                collapsed = set()
            else:
                # Matches the user closed since the last rebuild join the set;
                # ones they reopened leave it.
                collapsed = (collapsed | (previous_matching - previous_expanded)) - previous_expanded
            expanded = matching - collapsed
        elif previous_query:
            expanded = set(getattr(self, "_sample_tree_pre_search_expanded", set()) or set())
            self._sample_tree_pre_search_expanded = set()
            collapsed = set()
        else:
            expanded = previous_expanded
            collapsed = set()
        self._sample_tree_expand_query = search_query
        self._sample_tree_expand_matching = matching if search_query else set()
        self._sample_tree_search_collapsed = collapsed
        return expanded

    def _populate_metadata_tree(
        self,
        meta_subset: pl.DataFrame,
        contig_records: Dict[str, List[Dict[str, Any]]],
        filtered_counts: Dict[str, int],
        allowed_samples: Sequence[str],
    ) -> None:
        self._ensure_sample_tree_sort_sync()
        tree_state = self._snapshot_metadata_tree_state()
        self.tree.setUpdatesEnabled(False)
        sorting_was_enabled = self.tree.isSortingEnabled()
        self.tree.setSortingEnabled(False)
        self._apply_sample_tree_sort_indicator()
        try:
            self.tree.clear()
            self._sample_tree_items = {}
            self._sample_tree_widgets = {}
            self._tree_hover_shown = ()
            self._tree_hover_base_fonts = {}
            allowed_set = {str(s).strip() for s in allowed_samples if str(s).strip()}
            sample_names = sorted(allowed_set) if allowed_set else collect_sample_list(meta_subset)
            search_query = self._sample_tree_search_query()
            search_active = bool(search_query)
            tree_state["expanded_samples"] = self._search_expanded_samples(
                search_query, tree_state, sample_names, filtered_counts
            )

            meta_cols_available = [col for col in self.meta_columns_for_tree if col in meta_subset.columns]
            qc_cols_available = [
                col
                for col in self.run.sample_qc_columns
                if col in meta_subset.columns and col not in meta_cols_available
            ]
            meta_by_sample: Dict[str, Dict[str, Any]] = {}
            if not meta_subset.is_empty() and "Sample" in meta_subset.columns:
                select_cols = ["Sample", *meta_cols_available, *qc_cols_available]
                meta_rows = (
                    meta_subset.lazy()
                    .with_columns(pl.col("Sample").cast(pl.Utf8).str.strip_chars().alias("Sample"))
                    .select(select_cols)
                    .unique(subset=["Sample"], keep="first")
                    .collect()
                )
                for row in meta_rows.iter_rows(named=True):
                    sample_name = str(row.get("Sample") or "").strip()
                    if sample_name and sample_name not in meta_by_sample:
                        meta_by_sample[sample_name] = row

            sample_names = self._sorted_sample_names(sample_names, meta_by_sample)

            sample_db_map: Dict[str, Dict[str, int]] = {}
            sample_total_map: Dict[str, int] = {}
            contig_db_map: Dict[Tuple[str, str], Dict[str, int]] = {}
            contig_total_map: Dict[Tuple[str, str], int] = {}

            if not self.run.annotation_hit_counts_pl.is_empty() and allowed_set:
                counts = (
                    self.run.annotation_hit_counts_pl.lazy()
                    .filter(pl.col("Sample").cast(pl.Utf8).is_in(list(allowed_set)))
                    .collect()
                )
                for row in counts.iter_rows(named=True):
                    sample_name = str(row.get("Sample") or "").strip()
                    contig_name = str(row.get("Contig") or "").strip()
                    db_name = str(row.get("DB") or "").strip()
                    if not sample_name or not db_name:
                        continue
                    try:
                        hit_count = int(row.get("Hits") or 0)
                    except Exception:
                        continue

                    sample_bucket = sample_db_map.setdefault(sample_name, {})
                    sample_bucket[db_name] = sample_bucket.get(db_name, 0) + hit_count
                    sample_total_map[sample_name] = sample_total_map.get(sample_name, 0) + hit_count

                    if contig_name:
                        key = (sample_name, contig_name)
                        contig_bucket = contig_db_map.setdefault(key, {})
                        contig_bucket[db_name] = contig_bucket.get(db_name, 0) + hit_count
                        contig_total_map[key] = contig_total_map.get(key, 0) + hit_count

            max_total_hits = max((sample_total_map.get(s, 0) for s in sample_names), default=0)
            max_db_hits = {
                db: max((sample_db_map.get(s, {}).get(db, 0) for s in sample_names), default=0)
                for db in self.run.annotation_db_columns
            }

            contig_copies_map: Dict[Tuple[str, str], int] = {}
            if search_active:
                for rec_sample, records in contig_records.items():
                    for rec in records:
                        rec_contig = str(rec.get("Contig") or "").strip()
                        if not rec_contig:
                            continue
                        try:
                            rec_hits = int(rec.get("n_hits") or 0)
                        except (TypeError, ValueError):
                            rec_hits = 0
                        key = (str(rec_sample).strip(), rec_contig)
                        contig_copies_map[key] = contig_copies_map.get(key, 0) + rec_hits

            db_total_col_idx = 5 + len(self.meta_columns_for_tree)
            db_first_col_idx = db_total_col_idx + 1

            def _heat_color(value: int, max_value: int) -> Tuple[QtGui.QColor, QtGui.QColor]:
                if value <= 0 or max_value <= 0:
                    return QtGui.QColor("#ffffff"), QtGui.QColor("#1f2933")
                ratio = max(0.0, min(1.0, float(value) / float(max_value)))
                start = (247, 251, 255)  # very light blue
                end = (8, 81, 156)       # dark blue
                r = int(round(start[0] + (end[0] - start[0]) * ratio))
                g = int(round(start[1] + (end[1] - start[1]) * ratio))
                b = int(round(start[2] + (end[2] - start[2]) * ratio))
                fg = QtGui.QColor("#ffffff") if ratio >= 0.6 else QtGui.QColor("#1f2933")
                return QtGui.QColor(r, g, b), fg

            def _set_heat_cell(item: QtWidgets.QTreeWidgetItem, col_idx: int, value: int, max_value: int) -> None:
                bg, fg = _heat_color(value, max_value)
                item.setBackground(col_idx, QtGui.QBrush(bg))
                item.setForeground(col_idx, QtGui.QBrush(fg))
                item.setData(col_idx, _TREE_SORT_ROLE, int(value))

            for sample in sample_names:
                meta_vals = meta_by_sample.get(sample, {})
                country = self._country_from_metadata(meta_vals)
                meta_fields: List[str] = []
                for col in self.meta_columns_for_tree:
                    val = meta_vals.get(col, "")
                    text = clean_text(val)
                    meta_fields.append(text if text.lower() != "nan" else "")

                sample_hits = int(filtered_counts.get(sample, 0))
                sample_db_counts = sample_db_map.get(sample, {})
                sample_db_total = int(sample_total_map.get(sample, 0))
                sample_atb_total = int(self.run.sample_atb_counts.get(sample, 0))
                sample_db_values = [str(int(sample_db_counts.get(db, 0))) for db in self.run.annotation_db_columns]
                sample_item = _SortableTreeItem(
                    [
                        sample,
                        "",
                        country or "--",
                        str(int(self.run.sample_total_counts.get(sample, 0))),
                        str(sample_atb_total),
                        *meta_fields,
                        str(sample_db_total),
                        *sample_db_values,
                    ]
                )
                sample_item.setIcon(0, self._tree_vector_icon("genome", color="#0f766e"))
                sample_item.setIcon(1, self._tree_vector_icon("card"))
                sample_item.setIcon(2, self._country_badge_icon(country))
                sample_tooltip = self._sample_qc_tooltip(sample, country, meta_vals)
                if search_active:
                    copies_line = self._copies_tooltip_line(sample_hits, search_query)
                    sample_tooltip = f"{sample_tooltip}\n{copies_line}"
                sample_item.setToolTip(0, sample_tooltip)
                sample_item.setToolTip(1, "Open sample card")
                sample_item.setToolTip(2, f"Country: {country}" if country else "Country not available")
                sample_item.setSizeHint(0, QtCore.QSize(0, 26))
                sample_item.setForeground(0, QtGui.QBrush(QtGui.QColor("#1f2d3d")))
                sample_item.setForeground(1, QtGui.QBrush(QtGui.QColor("#2f6fed")))
                sample_item.setForeground(2, QtGui.QBrush(QtGui.QColor("#3f5163")))
                sample_item.setForeground(3, QtGui.QBrush(QtGui.QColor("#4f6275")))
                sample_item.setForeground(4, QtGui.QBrush(QtGui.QColor("#4f6275")))
                sample_item.setData(0, QtCore.Qt.ItemDataRole.UserRole, json.dumps({"sample": sample}))
                sample_item.setData(0, _TREE_SORT_ROLE, sample.casefold())
                sample_item.setData(2, _TREE_SORT_ROLE, (country or "").casefold())
                sample_item.setData(
                    3,
                    _TREE_SORT_ROLE,
                    int(self.run.sample_total_counts.get(sample, 0)),
                )
                sample_item.setData(4, _TREE_SORT_ROLE, sample_atb_total)
                for offset, text in enumerate(meta_fields):
                    sample_item.setData(5 + offset, _TREE_SORT_ROLE, text.casefold())
                self._sample_tree_items[sample] = sample_item
                if sample_hits == 0:
                    brush = QtGui.QBrush(QtGui.QColor("#888888"))
                    for col_idx in range(5 + len(meta_fields)):
                        sample_item.setForeground(col_idx, brush)
                _set_heat_cell(sample_item, db_total_col_idx, sample_db_total, max_total_hits)
                for offset, db in enumerate(self.run.annotation_db_columns):
                    _set_heat_cell(
                        sample_item,
                        db_first_col_idx + offset,
                        int(sample_db_counts.get(db, 0)),
                        int(max_db_hits.get(db, 0)),
                    )
                self.tree.addTopLevelItem(sample_item)
                quality_text, quality_color = self._quality_badge_style(meta_vals.get("Quality"))
                completeness_text, completeness_color = self._completeness_cell(
                    meta_vals.get("CheckM2Completeness")
                )
                contamination_text, contamination_color = self._contamination_cell(
                    meta_vals.get("CheckM2Contamination")
                )
                quast_tooltip = self._quast_tooltip(meta_vals)
                completeness_tooltip = self._checkm2_tooltip(meta_vals, focus="completeness")
                contamination_tooltip = self._checkm2_tooltip(meta_vals, focus="contamination")

                # ST column: the CSV's subtype, else the run's own MLST ST (run data only)
                st_value = meta_vals.get("subtype") if "subtype" in meta_vals else meta_vals.get("ST")
                subtype_val = self._tree_metric_float(st_value)
                sample_item.setText(self._col_st, self._st_label(st_value))
                sample_item.setForeground(self._col_st, QtGui.QBrush(QtGui.QColor("#54657a")))
                sample_item.setToolTip(self._col_st, "Sequence type (MLST)")
                if subtype_val is not None:
                    sample_item.setData(self._col_st, _TREE_SORT_ROLE, float(subtype_val))

                # Quality / Completeness / Contamination as coloured numbers
                qc_font = sample_item.font(self._col_quality)
                qc_font.setWeight(QtGui.QFont.Weight.DemiBold)
                for col_idx, text, color, tooltip, value in (
                    (self._col_quality, quality_text, quality_color, quast_tooltip,
                     self._tree_metric_float(meta_vals.get("Quality"))),
                    (self._col_completeness, completeness_text, completeness_color, completeness_tooltip,
                     self._tree_metric_float(meta_vals.get("CheckM2Completeness"))),
                    (self._col_contamination, contamination_text, contamination_color, contamination_tooltip,
                     self._tree_metric_float(meta_vals.get("CheckM2Contamination"))),
                ):
                    sample_item.setText(col_idx, text)
                    sample_item.setForeground(col_idx, QtGui.QBrush(QtGui.QColor(color)))
                    sample_item.setFont(col_idx, qc_font)
                    sample_item.setToolTip(col_idx, tooltip)
                    if value is not None:
                        sample_item.setData(col_idx, _TREE_SORT_ROLE, float(value))
                for contig in self._display_contig_records(contig_records.get(sample, [])):
                    contig_name = str(contig["Contig"]).strip()
                    member_contigs = [
                        str(member).strip()
                        for member in (contig.get("__member_contigs") or [contig_name])
                        if str(member).strip()
                    ]
                    if not member_contigs and contig_name:
                        member_contigs = [contig_name]
                    contig_keys = [(sample, member) for member in member_contigs]
                    contig_total_hits = int(sum(contig_total_map.get(key, 0) for key in contig_keys))
                    contig_atb_total = int(sum(int(self.run.contig_atb_counts.get(key, 0)) for key in contig_keys))
                    contig_db_counts: Dict[str, int] = {}
                    for key in contig_keys:
                        for db, count in contig_db_map.get(key, {}).items():
                            contig_db_counts[db] = contig_db_counts.get(db, 0) + int(count)
                    is_plasmid = self._truthy_text(contig.get("MobSuiteIsPlasmid"))
                    plasmid_label = "P" if is_plasmid else "C"
                    is_plasmid_group = bool(contig.get("__is_plasmid_group"))
                    display_label = str(contig.get("__display_label") or f"Contig {contig_name}")
                    n_total = int(
                        sum(
                            int(self.run.contig_total_count_map.get((sample, member), 0))
                            for member in member_contigs
                        )
                    )
                    if n_total <= 0:
                        n_total = int(contig.get("n_total", 0))
                    contig_item = _SortableTreeItem(
                        [
                            display_label,
                            "",
                            plasmid_label,
                            str(n_total),
                            str(contig_atb_total),
                            *["" for _ in self.meta_columns_for_tree],
                            str(contig_total_hits),
                            *[
                                str(
                                    int(
                                        contig_db_counts.get(db, 0)
                                    )
                                )
                                for db in self.run.annotation_db_columns
                            ],
                        ]
                    )
                    molecule_icon = (
                        self._tree_vector_icon("plasmid", color="#16a34a")
                        if is_plasmid
                        else self._tree_vector_icon("chromosome", color="#2563eb")
                    )
                    contig_item.setIcon(0, molecule_icon)
                    contig_item.setIcon(2, molecule_icon)
                    contig_item.setSizeHint(0, QtCore.QSize(0, 31))
                    contig_item.setTextAlignment(1, int(QtCore.Qt.AlignmentFlag.AlignCenter | QtCore.Qt.AlignmentFlag.AlignVCenter))
                    contig_item.setTextAlignment(2, int(QtCore.Qt.AlignmentFlag.AlignCenter | QtCore.Qt.AlignmentFlag.AlignVCenter))
                    contig_item.setForeground(0, QtGui.QBrush(QtGui.QColor("#5f7184")))
                    contig_item.setForeground(2, QtGui.QBrush(QtGui.QColor("#16a34a" if is_plasmid else "#2563eb")))
                    contig_item.setForeground(3, QtGui.QBrush(QtGui.QColor("#718194")))
                    contig_item.setForeground(4, QtGui.QBrush(QtGui.QColor("#718194")))
                    contig_item.setData(0, _TREE_SORT_ROLE, display_label.casefold())
                    contig_item.setData(2, _TREE_SORT_ROLE, plasmid_label.casefold())
                    contig_item.setData(3, _TREE_SORT_ROLE, n_total)
                    contig_item.setData(4, _TREE_SORT_ROLE, contig_atb_total)
                    contig_font = contig_item.font(2)
                    contig_font.setBold(True)
                    contig_item.setFont(2, contig_font)
                    if is_plasmid_group:
                        member_label = ", ".join(member_contigs)
                        plasmid_index = self.run.plasmid_cluster_index
                        cluster_tooltip = ""
                        if plasmid_index is not None and hasattr(plasmid_index, "tooltip_for_cluster"):
                            try:
                                cluster_tooltip = plasmid_index.tooltip_for_cluster(
                                    plasmid_cluster_key_from_row(contig)
                                )
                            except Exception:
                                cluster_tooltip = ""
                        tooltip = f"MOB-suite plasmid group from contigs: {member_label}"
                        if cluster_tooltip:
                            tooltip = f"{cluster_tooltip}\n\n{tooltip}"
                        contig_item.setToolTip(0, tooltip)
                        contig_item.setToolTip(2, f"MOB-suite: plasmid group ({len(member_contigs)} contigs)")
                    else:
                        tooltip = self._contig_tooltip(
                            sample=sample,
                            display_label=display_label,
                            molecule_label=plasmid_label,
                            n_total=n_total,
                            member_contigs=member_contigs,
                            contig=contig,
                        )
                        contig_item.setToolTip(0, tooltip)
                        contig_item.setToolTip(2, "MOB-suite: plasmid contig" if is_plasmid else "MOB-suite: chromosomal contig")
                    contig_item.setData(
                        0,
                        QtCore.Qt.ItemDataRole.UserRole,
                        json.dumps(
                            {
                                "sample": sample,
                                "contig": contig_name,
                                "contigs": member_contigs,
                                "contig_label": display_label,
                                "is_plasmid_group": is_plasmid_group,
                            }
                        ),
                    )
                    if search_active:
                        contig_copies = int(sum(contig_copies_map.get(key, 0) for key in contig_keys))
                        copies_line = self._copies_tooltip_line(
                            contig_copies,
                            search_query,
                            "plasmid group" if is_plasmid_group else "contig",
                        )
                        contig_item.setToolTip(0, f"{contig_item.toolTip(0)}\n{copies_line}")
                    _set_heat_cell(contig_item, db_total_col_idx, contig_total_hits, max_total_hits)
                    for offset, db in enumerate(self.run.annotation_db_columns):
                        _set_heat_cell(
                            contig_item,
                            db_first_col_idx + offset,
                            int(contig_db_counts.get(db, 0)),
                            int(max_db_hits.get(db, 0)),
                        )
                    sample_item.addChild(contig_item)
            self._apply_tree_sample_hover(self._tree_hover_samples())
        finally:
            if sorting_was_enabled:
                self.tree.setSortingEnabled(True)
            self.tree.setUpdatesEnabled(True)
        self._restore_metadata_tree_state(tree_state)

    def _sample_sort_settings(self) -> Tuple[str, bool]:
        sort_col = "__sample__"
        ascending = True
        col_box = getattr(self, "sample_sort_column_box", None)
        if isinstance(col_box, QtWidgets.QComboBox):
            data = col_box.currentData()
            sort_col = str(data if data is not None else col_box.currentText() or "__sample__")
        order_box = getattr(self, "sample_sort_order_box", None)
        if isinstance(order_box, QtWidgets.QComboBox):
            data = order_box.currentData()
            order = str(data if data is not None else order_box.currentText()).strip().lower()
            ascending = order.startswith("asc")
        return sort_col, ascending

    def _sorted_sample_names(
        self,
        sample_names: Sequence[str],
        meta_by_sample: Dict[str, Dict[str, Any]],
    ) -> List[str]:
        names = [str(s).strip() for s in sample_names if str(s).strip()]
        if not names:
            return []
        sort_col, ascending = self._sample_sort_settings()
        if sort_col == "__sample__":
            return sorted(names, key=lambda s: s.lower(), reverse=not ascending)

        present: List[Tuple[str, str]] = []
        missing: List[str] = []
        for sample in names:
            row = meta_by_sample.get(sample, {})
            value = clean_text(row.get(sort_col, ""))
            if value:
                present.append((sample, value))
            else:
                missing.append(sample)

        present_sorted = sorted(
            present,
            key=lambda sv: (str(sv[1]).lower(), str(sv[0]).lower()),
            reverse=not ascending,
        )
        missing_sorted = sorted(missing, key=lambda s: s.lower(), reverse=not ascending)
        return [sample for sample, _value in present_sorted] + missing_sorted


    def _on_item_expanded(self, _item: QtWidgets.QTreeWidgetItem) -> None:
        if getattr(self, "_tree_bulk_toggle_in_progress", False):
            return
        item = _item
        if item is None or item.parent() is not None:
            return
        mods = QtWidgets.QApplication.keyboardModifiers()
        if not (mods & QtCore.Qt.KeyboardModifier.ControlModifier):
            return
        self._tree_bulk_toggle_in_progress = True
        try:
            self.tree.expandAll()
        finally:
            self._tree_bulk_toggle_in_progress = False


    def _on_item_collapsed(self, _item: QtWidgets.QTreeWidgetItem) -> None:
        if getattr(self, "_tree_bulk_toggle_in_progress", False):
            return
        item = _item
        if item is None or item.parent() is not None:
            return
        mods = QtWidgets.QApplication.keyboardModifiers()
        if not (mods & QtCore.Qt.KeyboardModifier.ControlModifier):
            return
        self._tree_bulk_toggle_in_progress = True
        try:
            self.tree.collapseAll()
        finally:
            self._tree_bulk_toggle_in_progress = False


    def _on_tree_clicked(self, item: QtWidgets.QTreeWidgetItem, col: int) -> None:
        payload_raw = item.data(0, QtCore.Qt.ItemDataRole.UserRole)
        if not payload_raw:
            return
        payload = json.loads(payload_raw)
        sample = payload.get("sample")
        contig = payload.get("contig")
        contigs = [
            str(member).strip()
            for member in (payload.get("contigs") or [])
            if str(member).strip()
        ]
        contig_label = str(payload.get("contig_label") or contig or "").strip()
        if col == 1 and sample and contig is None:
            open_sample_card = getattr(self, "_open_sample_card_for_sample", None)
            if callable(open_sample_card):
                open_sample_card(str(sample))
            return
        # One selection change (SelectionChanged on the bus): the Annotations
        # panels in Contig or Around matches mode load the contig, the paired
        # one with its track (plan 4.5).
        if sample and contig is not None:
            self._igv_dirty = True
            changed = self._select_contig(
                sample,
                contig,
                contigs or [str(contig)],
                contig_label or str(contig),
                origin=SAMPLE_TREE,
            )
            if not changed:
                self._render_annotations_for_selection()
        else:
            self._select_contig(None, None, origin=SAMPLE_TREE)

    def _on_tree_double_clicked(self, item: QtWidgets.QTreeWidgetItem, col: int) -> None:
        # Card column was folded into the row; open the sample card on a
        # double-click of a sample (top-level) row instead.
        if item is None or item.parent() is not None:
            return
        payload_raw = item.data(0, QtCore.Qt.ItemDataRole.UserRole)
        if not payload_raw:
            return
        try:
            payload = json.loads(payload_raw)
        except (TypeError, ValueError):
            return
        sample = payload.get("sample")
        if not sample or payload.get("contig") is not None:
            return
        open_sample_card = getattr(self, "_open_sample_card_for_sample", None)
        if callable(open_sample_card):
            open_sample_card(str(sample))


    def _show_sample_tree_context_menu(self, position: QtCore.QPoint) -> None:
        item = self.tree.itemAt(position)
        if item is None:
            return
        # Select the clicked row so it is obvious which one the menu acts on.
        self.tree.setCurrentItem(item)
        sample_item = item
        while sample_item.parent() is not None:
            sample_item = sample_item.parent()
        sample = self._sample_from_tree_item(sample_item)
        if not sample:
            return

        payload = self._tree_item_payload(item)
        contig = str(payload.get("contig") or "").strip()
        if contig:
            label = str(payload.get("contig_label") or contig).strip()
            menu = QtWidgets.QMenu(self.tree)
            header_action = menu.addAction(f"{sample} · {label}")
            header_action.setEnabled(False)
            menu.addSeparator()
            annotations_action = menu.addAction("Open annotations")
            context_action = menu.addAction("Open in Context view")
            copy_action = menu.addAction(
                "Copy contig IDs" if len(payload.get("contigs") or []) > 1 else "Copy contig ID"
            )
            menu.addSeparator()
            card_action = menu.addAction("Open sample card")

            chosen = menu.exec(self.tree.viewport().mapToGlobal(position))
            if chosen is annotations_action:
                self._open_tree_contig_payload(payload, mode="Contig")
            elif chosen is context_action:
                self._open_tree_contig_payload(payload, mode="Context")
            elif chosen is copy_action:
                members = [
                    str(member).strip()
                    for member in (payload.get("contigs") or [contig])
                    if str(member).strip()
                ]
                QtWidgets.QApplication.clipboard().setText("\n".join(members or [contig]))
                status = getattr(self, "status", None)
                if isinstance(status, QtWidgets.QLabel):
                    status.setText("Copied contig ID" if len(members) <= 1 else "Copied contig IDs")
            elif chosen is card_action:
                open_sample_card = getattr(self, "_open_sample_card_for_sample", None)
                if callable(open_sample_card):
                    open_sample_card(str(sample))
            return

        menu = QtWidgets.QMenu(self.tree)
        header_action = menu.addAction(sample)
        header_action.setEnabled(False)
        menu.addSeparator()
        expand_action = None
        collapse_action = None
        if sample_item.childCount() > 0:
            if sample_item.isExpanded():
                collapse_action = menu.addAction("Collapse")
            else:
                expand_action = menu.addAction("Expand")
            menu.addSeparator()
        expand_all_action = menu.addAction("Expand all samples")
        collapse_all_action = menu.addAction("Collapse all samples")
        menu.addSeparator()
        card_action = menu.addAction("Open sample card")

        chosen = menu.exec(self.tree.viewport().mapToGlobal(position))
        if chosen is None:
            return
        if chosen is expand_action:
            sample_item.setExpanded(True)
        elif chosen is collapse_action:
            sample_item.setExpanded(False)
        elif chosen is expand_all_action:
            self._expand_all_samples()
        elif chosen is collapse_all_action:
            self._collapse_all_samples()
        elif chosen is card_action:
            open_sample_card = getattr(self, "_open_sample_card_for_sample", None)
            if callable(open_sample_card):
                open_sample_card(str(sample))

    def _open_tree_contig_payload(self, payload: Dict[str, Any], *, mode: str) -> None:
        sample = str(payload.get("sample") or "").strip()
        contig = str(payload.get("contig") or "").strip()
        if not sample or not contig:
            return
        contigs = [
            str(member).strip()
            for member in (payload.get("contigs") or [])
            if str(member).strip()
        ]
        display_label = str(payload.get("contig_label") or contig).strip()
        self._igv_dirty = True
        # One draw per panel for the selection and the mode switch together.
        with self._batched_annotations_render():
            changed = self._select_contig(sample, contig, contigs or [contig], display_label or contig, origin=SAMPLE_TREE)
            if not changed:
                self._render_annotations_for_selection()
            # Switch the Annotations panel used last to the mode and show it
            # (D3: no compute; plan 4.5: an open request plus the selection).
            self.ctx.request(FocusPanel(ANNOTATIONS, params_tuple({"view_mode": mode})), origin=SAMPLE_TREE)

    def _expand_all_samples(self) -> None:
        self._tree_bulk_toggle_in_progress = True
        try:
            self.tree.expandAll()
        finally:
            self._tree_bulk_toggle_in_progress = False

    def _collapse_all_samples(self) -> None:
        self._tree_bulk_toggle_in_progress = True
        try:
            self.tree.collapseAll()
        finally:
            self._tree_bulk_toggle_in_progress = False


    def _tree_hover_samples(self) -> Tuple[str, ...]:
        """The samples hovered in this window (bold in the tree)."""
        ctx = getattr(self, "ctx", None)
        hover = getattr(ctx, "hover_state", None)
        return tuple(hover.samples) if hover is not None else ()

    def _apply_tree_sample_hover(self, samples: Any) -> None:
        """Bold the rows of the hovered samples (one, several or none)."""
        if samples is None:
            samples = ()
        elif isinstance(samples, str):
            samples = (samples,)
        wanted = tuple(dict.fromkeys(str(s) for s in samples if s))
        shown = tuple(getattr(self, "_tree_hover_shown", ()) or ())
        if wanted == shown:
            return
        base_fonts = getattr(self, "_tree_hover_base_fonts", None)
        if not isinstance(base_fonts, dict):
            base_fonts = {}
            self._tree_hover_base_fonts = base_fonts
        widgets = getattr(self, "_sample_tree_widgets", {})
        keep = set(wanted)
        for sample in shown:
            if sample in keep:
                continue
            old_item = self._sample_tree_items.get(sample)
            old_font = base_fonts.pop(sample, None)
            if old_item is not None and isinstance(old_font, QtGui.QFont):
                old_item.setFont(0, old_font)
            old_widget = widgets.get(sample)
            if old_widget is not None and hasattr(old_widget, "set_highlighted"):
                old_widget.set_highlighted(False)
        already = set(shown)
        for sample in wanted:
            if sample in already:
                continue
            item = self._sample_tree_items.get(sample)
            if item is not None:
                base_font = QtGui.QFont(item.font(0))
                base_fonts[sample] = base_font
                font = QtGui.QFont(base_font)
                font.setBold(True)
                item.setFont(0, font)
            widget = widgets.get(sample)
            if widget is not None and hasattr(widget, "set_highlighted"):
                widget.set_highlighted(True)
        self._tree_hover_shown = wanted


    def _sample_from_tree_item(self, item: Optional[QtWidgets.QTreeWidgetItem]) -> Optional[str]:
        if item is None:
            return None
        while item.parent() is not None:
            item = item.parent()
        payload = self._tree_item_payload(item)
        sample = payload.get("sample")
        return str(sample) if sample else None
