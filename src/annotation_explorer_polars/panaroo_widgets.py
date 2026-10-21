from __future__ import annotations

from collections import Counter
from typing import List, Optional, Sequence, Tuple

import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets

try:
    import pyqtgraph as pg
except Exception:
    pg = None


class _PgUnavailableWidget(QtWidgets.QWidget):
    def __init__(self, text: str, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        layout = QtWidgets.QVBoxLayout(self)
        label = QtWidgets.QLabel(text)
        label.setWordWrap(True)
        label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(label, stretch=1)


class _MetadataPlaceholderView(QtWidgets.QWidget):
    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(24, 24, 24, 24)
        layout.setSpacing(10)
        title = QtWidgets.QLabel("Metadata")
        title.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        title.setStyleSheet("font-size: 22px; font-weight: 700; color: #1f2c3b;")
        body = QtWidgets.QLabel("Metadata view placeholder")
        body.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        body.setStyleSheet("font-size: 14px; color: #607183;")
        layout.addStretch(1)
        layout.addWidget(title)
        layout.addWidget(body)
        layout.addStretch(1)


class _DonutChartWidget(QtWidgets.QWidget):
    _PALETTE = [
        "#1d4ed8",
        "#dc2626",
        "#059669",
        "#7c3aed",
        "#ea580c",
        "#0f766e",
        "#be185d",
        "#4f46e5",
        "#65a30d",
        "#b45309",
    ]

    def __init__(self, title: str, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self._title = str(title)
        self._items: List[Tuple[str, int]] = []
        self._total = 0
        self._center_text = ""
        self.setMinimumHeight(240)

    def set_items(self, items: Sequence[Tuple[str, int]], total: int, center_text: str) -> None:
        self._items = [(str(label), int(count)) for label, count in items if int(count) > 0]
        self._total = max(0, int(total))
        self._center_text = str(center_text)
        self.setToolTip("\n".join(f"{label}: {count}" for label, count in self._items))
        self.update()

    def paintEvent(self, event: QtGui.QPaintEvent) -> None:
        super().paintEvent(event)
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)

        rect = self.rect().adjusted(10, 10, -10, -10)
        title_font = QtGui.QFont(self.font())
        title_font.setBold(True)
        title_font.setPointSizeF(max(11.0, title_font.pointSizeF() + 1.5))
        painter.setFont(title_font)
        painter.setPen(QtGui.QColor("#1f2c3b"))
        painter.drawText(
            QtCore.QRectF(rect.left(), rect.top(), rect.width(), 28),
            int(QtCore.Qt.AlignmentFlag.AlignHCenter | QtCore.Qt.AlignmentFlag.AlignVCenter),
            self._title,
        )

        donut_box = QtCore.QRectF(rect.left(), rect.top() + 36, rect.width(), min(170.0, rect.height() - 82.0))
        diameter = min(donut_box.width(), donut_box.height())
        donut_rect = QtCore.QRectF(
            donut_box.center().x() - diameter * 0.5,
            donut_box.top(),
            diameter,
            diameter,
        )

        if self._total <= 0 or not self._items:
            painter.setPen(QtGui.QColor("#94a3b8"))
            painter.drawText(
                donut_rect,
                int(QtCore.Qt.AlignmentFlag.AlignCenter),
                "No data",
            )
            return

        start_angle = 90.0 * 16.0
        for idx, (_label, count) in enumerate(self._items):
            span_angle = -360.0 * 16.0 * (float(count) / float(self._total))
            color = QtGui.QColor(self._PALETTE[idx % len(self._PALETTE)])
            painter.setBrush(QtGui.QBrush(color))
            painter.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
            painter.drawPie(donut_rect, int(round(start_angle)), int(round(span_angle)))
            start_angle += span_angle

        hole_ratio = 0.56
        hole_rect = donut_rect.adjusted(
            donut_rect.width() * hole_ratio * 0.5,
            donut_rect.height() * hole_ratio * 0.5,
            -donut_rect.width() * hole_ratio * 0.5,
            -donut_rect.height() * hole_ratio * 0.5,
        )
        painter.setBrush(QtGui.QBrush(QtGui.QColor("#ffffff")))
        painter.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
        painter.drawEllipse(hole_rect)

        center_font = QtGui.QFont(self.font())
        center_font.setBold(True)
        center_font.setPointSizeF(max(12.0, center_font.pointSizeF() + 1.0))
        painter.setFont(center_font)
        painter.setPen(QtGui.QColor("#1f2c3b"))
        painter.drawText(
            hole_rect,
            int(QtCore.Qt.AlignmentFlag.AlignCenter),
            self._center_text or str(self._total),
        )

        legend_top = donut_rect.bottom() + 10.0
        legend_left = rect.left() + 6.0
        legend_width = rect.width() - 12.0
        legend_font = QtGui.QFont(self.font())
        legend_font.setPointSizeF(max(9.0, legend_font.pointSizeF()))
        painter.setFont(legend_font)
        line_h = 18.0
        for idx, (label, count) in enumerate(self._items[:5]):
            y = legend_top + float(idx) * line_h
            swatch_rect = QtCore.QRectF(legend_left, y + 3.0, 10.0, 10.0)
            painter.setBrush(QtGui.QBrush(QtGui.QColor(self._PALETTE[idx % len(self._PALETTE)])))
            painter.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
            painter.drawRoundedRect(swatch_rect, 2.0, 2.0)
            painter.setPen(QtGui.QColor("#475569"))
            text = f"{label} ({count})"
            painter.drawText(
                QtCore.QRectF(legend_left + 16.0, y, legend_width - 16.0, line_h),
                int(QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignVCenter),
                text,
            )


class _MetadataDonutPanel(QtWidgets.QWidget):
    _COUNTRY_CANDIDATES = ("Country",)
    _K_LOCUS_CANDIDATES = (
        "K locus",
        "K_locus",
        "KLocus",
        "KL",
        "KL type",
        "K type",
    )
    _SUBTYPE_CANDIDATES = (
        "Subtype",
        "subtype",
        "SubType",
    )

    def __init__(self, meta_df: pd.DataFrame, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self._meta_df = meta_df.copy()
        self._sample_filter: Optional[set[str]] = None

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(10)

        self._summary = QtWidgets.QLabel("")
        self._summary.setStyleSheet("color: #607183; font-weight: 600;")
        layout.addWidget(self._summary)

        charts_row = QtWidgets.QHBoxLayout()
        charts_row.setContentsMargins(0, 0, 0, 0)
        charts_row.setSpacing(12)
        layout.addLayout(charts_row, stretch=1)

        self._country_chart = _DonutChartWidget("Country", self)
        self._k_locus_chart = _DonutChartWidget("K locus", self)
        self._subtype_chart = _DonutChartWidget("Subtype", self)
        for chart in (self._country_chart, self._k_locus_chart, self._subtype_chart):
            charts_row.addWidget(chart, stretch=1)

        self.refresh()

    @staticmethod
    def _resolve_column(columns: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
        return _resolve_column(columns, candidates)

    def set_sample_filter(self, samples: Optional[Sequence[str]]) -> None:
        next_filter = {str(s).strip() for s in samples if str(s).strip()} if samples is not None else None
        if next_filter == self._sample_filter:
            return
        self._sample_filter = next_filter
        self.refresh()

    def _filtered_meta(self) -> pd.DataFrame:
        if self._meta_df.empty or "Sample" not in self._meta_df.columns:
            return self._meta_df.iloc[0:0].copy()
        if self._sample_filter is None:
            return self._meta_df.copy()
        return self._meta_df[self._meta_df["Sample"].astype(str).str.strip().isin(self._sample_filter)].copy()

    @staticmethod
    def _value_counts(df: pd.DataFrame, column: Optional[str]) -> Tuple[List[Tuple[str, int]], int]:
        if not column or column not in df.columns:
            return [], 0
        values = (
            df[column]
            .fillna("")
            .astype(str)
            .str.strip()
        )
        values = values[values != ""]
        if values.empty:
            return [], 0
        counts = Counter(values.tolist())
        ordered = sorted(counts.items(), key=lambda item: (-int(item[1]), str(item[0]).casefold()))
        if len(ordered) > 5:
            top = ordered[:4]
            other_count = sum(count for _label, count in ordered[4:])
            ordered = top + [("Other", other_count)]
        total = int(sum(count for _label, count in ordered))
        return ordered, total

    def refresh(self) -> None:
        df = self._filtered_meta()
        sample_count = int(len(df.index))
        self._summary.setText(f"Samples in view: {sample_count}")

        columns = list(df.columns)
        country_col = self._resolve_column(columns, self._COUNTRY_CANDIDATES)
        klocus_col = self._resolve_column(columns, self._K_LOCUS_CANDIDATES)
        subtype_col = self._resolve_column(columns, self._SUBTYPE_CANDIDATES)

        for chart, column in (
            (self._country_chart, country_col),
            (self._k_locus_chart, klocus_col),
            (self._subtype_chart, subtype_col),
        ):
            items, total = self._value_counts(df, column)
            center_text = str(total) if total > 0 else "0"
            chart.set_items(items, total, center_text)

def _resolve_column(columns: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    """The first candidate present in ``columns`` (exact, then case blind)."""
    available = {str(col): str(col) for col in columns}
    lowered = {str(col).casefold(): str(col) for col in columns}
    for candidate in candidates:
        if candidate in available:
            return available[candidate]
        match = lowered.get(str(candidate).casefold())
        if match:
            return match
    return None


def _high_contrast_color(index: int) -> QtGui.QColor:
    # Golden-angle hue distribution plus alternating saturation/value gives
    # stronger separation than a short fixed palette when there are many
    # categorical values, which matters for many-valued columns such as K locus.
    hue = (0.19 + (float(index) * 0.61803398875)) % 1.0
    saturation = 0.92 if (index % 2 == 0) else 0.72
    value = 0.90 if ((index // 2) % 2 == 0) else 0.68
    return QtGui.QColor.fromHsvF(hue, saturation, value, 1.0)


class _MetaValueChip(QtWidgets.QFrame):
    hovered = QtCore.pyqtSignal(str)
    unhovered = QtCore.pyqtSignal(str)

    def __init__(self, value: str, count: int, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self._value = str(value)
        self.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.setObjectName("MetaValueChip")
        self.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed,
        )
        self.setFixedHeight(30)
        self.setStyleSheet(
            "#MetaValueChip {"
            "  border: 1px solid #d8e2ed;"
            "  border-radius: 10px;"
            "  background: #f8fbff;"
            "  padding: 2px 8px;"
            "}"
            "#MetaValueChip:hover {"
            "  border: 1px solid #97b6f6;"
            "  background: #eef4ff;"
            "}"
        )
        row = QtWidgets.QHBoxLayout(self)
        row.setContentsMargins(6, 3, 6, 3)
        row.setSpacing(6)
        text = QtWidgets.QLabel(self._value)
        text.setStyleSheet("color: #1f2c3b;")
        row.addWidget(text)
        cnt = QtWidgets.QLabel(str(int(count)))
        cnt.setStyleSheet(
            "background:#e9eef6; color:#4f6275; border-radius:8px; padding:0 6px; min-width:16px;"
        )
        cnt.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        row.addWidget(cnt)

    def enterEvent(self, event: QtCore.QEvent) -> None:
        self.hovered.emit(self._value)
        super().enterEvent(event)

    def leaveEvent(self, event: QtCore.QEvent) -> None:
        self.unhovered.emit(self._value)
        super().leaveEvent(event)
