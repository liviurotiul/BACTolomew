from __future__ import annotations

from collections import Counter
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets

from .tree_view import RaxmlTreeView
from .panaroo_common import _METADATA_STRIP_PALETTE, _NEO_TREE_RESISTANCE_TOKEN, _natural_sort_key
from .panaroo_style import _apply_analysis_panel_style

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


class _DeferredViewPlaceholder(QtWidgets.QWidget):
    def __init__(self, title: str, message: str, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(24, 24, 24, 24)
        layout.setSpacing(10)
        title_label = QtWidgets.QLabel(str(title))
        title_label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        title_label.setStyleSheet("font-size: 22px; font-weight: 700; color: #1f2c3b;")
        body = QtWidgets.QLabel(str(message))
        body.setWordWrap(True)
        body.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        body.setStyleSheet("font-size: 14px; color: #607183;")
        layout.addStretch(1)
        layout.addWidget(title_label)
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
        available = {str(col): str(col) for col in columns}
        lowered = {str(col).casefold(): str(col) for col in columns}
        for candidate in candidates:
            if candidate in available:
                return available[candidate]
            match = lowered.get(str(candidate).casefold())
            if match:
                return match
        return None

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

def _high_contrast_color(index: int) -> QtGui.QColor:
    # Golden-angle hue distribution plus alternating saturation/value gives
    # stronger separation than a short fixed palette when there are many
    # categorical values, which is important for the neoTree K-locus strip.
    hue = (0.19 + (float(index) * 0.61803398875)) % 1.0
    saturation = 0.92 if (index % 2 == 0) else 0.72
    value = 0.90 if ((index // 2) % 2 == 0) else 0.68
    return QtGui.QColor.fromHsvF(hue, saturation, value, 1.0)


class _LegendSwatchRow(QtWidgets.QWidget):
    def __init__(
        self,
        label: str,
        color: QtGui.QColor,
        count: Optional[int],
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Maximum,
            QtWidgets.QSizePolicy.Policy.Preferred,
        )
        row = QtWidgets.QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(5)

        swatch = QtWidgets.QFrame(self)
        swatch.setFixedSize(12, 12)
        swatch.setStyleSheet(
            "border: 1px solid rgba(15, 23, 42, 0.18);"
            f"border-radius: 3px; background: {QtGui.QColor(color).name()};"
        )
        row.addWidget(swatch, 0, QtCore.Qt.AlignmentFlag.AlignTop)

        text_label = f"{label} ({int(count)})" if count is not None else str(label)
        text = QtWidgets.QLabel(text_label, self)
        text.setWordWrap(True)
        text.setStyleSheet("color: #334155;")
        text.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Maximum,
            QtWidgets.QSizePolicy.Policy.Preferred,
        )
        row.addWidget(text, 0, QtCore.Qt.AlignmentFlag.AlignTop)


class _CategoricalLegendPanel(QtWidgets.QWidget):
    def __init__(
        self,
        title: str,
        *,
        value_prefix: str = "",
        scrollable: bool = True,
        grid_columns: int = 1,
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self._title = str(title).strip() or "Legend"
        self._value_prefix = str(value_prefix)
        self._scrollable = bool(scrollable)
        self._grid_columns = max(1, int(grid_columns))
        self._sample_filter: Optional[set[str]] = None
        self._value_by_sample: Dict[str, str] = {}
        self._color_by_value: Dict[str, QtGui.QColor] = {}
        self._missing_label = "Missing"
        self._show_missing = False

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(2, 0, 0, 0)
        layout.setSpacing(4)

        title_label = QtWidgets.QLabel(self._title, self)
        title_label.setStyleSheet("color: #334155; font-weight: 700; padding-left: 4px;")
        layout.addWidget(title_label)

        self._legend_scroll: Optional[QtWidgets.QScrollArea] = None
        self._legend_wrap = QtWidgets.QWidget(self)
        if self._scrollable:
            self._legend_scroll = QtWidgets.QScrollArea(self)
            self._legend_scroll.setWidgetResizable(True)
            self._legend_scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
            self._legend_scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
            self._legend_scroll.setMinimumHeight(110)
            self._legend_scroll.setMaximumHeight(140)
            self._legend_layout = QtWidgets.QVBoxLayout(self._legend_wrap)
            self._legend_layout.setContentsMargins(0, 0, 0, 0)
            self._legend_layout.setSpacing(4)
            self._legend_scroll.setWidget(self._legend_wrap)
            layout.addWidget(self._legend_scroll, stretch=0)
        else:
            self._legend_layout = QtWidgets.QGridLayout(self._legend_wrap)
            self._legend_layout.setContentsMargins(0, 0, 0, 0)
            self._legend_layout.setHorizontalSpacing(8)
            self._legend_layout.setVerticalSpacing(2)
            layout.addWidget(self._legend_wrap, stretch=0, alignment=QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignTop)

        self._rebuild_legend()

    def set_sample_filter(self, samples: Optional[Sequence[str]]) -> None:
        next_filter = {str(s).strip() for s in samples if str(s).strip()} if samples is not None else None
        if next_filter == self._sample_filter:
            return
        self._sample_filter = next_filter
        self._rebuild_legend()

    def set_assignments(
        self,
        value_by_sample: Mapping[str, str],
        color_by_value: Mapping[str, object],
        *,
        missing_label: str = "Missing",
        show_missing: bool = False,
    ) -> None:
        normalized_values: Dict[str, str] = {}
        for raw_sample, raw_value in value_by_sample.items():
            sample = str(raw_sample or "").strip()
            if not sample:
                continue
            value = str(raw_value or "").strip() or str(missing_label or "Missing").strip() or "Missing"
            normalized_values[sample] = value

        normalized_colors: Dict[str, QtGui.QColor] = {}
        for raw_value, raw_color in color_by_value.items():
            value = str(raw_value or "").strip()
            color = QtGui.QColor(raw_color) if not isinstance(raw_color, QtGui.QColor) else QtGui.QColor(raw_color)
            if not value or not color.isValid():
                continue
            normalized_colors[value] = color

        self._value_by_sample = normalized_values
        self._color_by_value = normalized_colors
        self._missing_label = str(missing_label or "Missing").strip() or "Missing"
        self._show_missing = bool(show_missing)
        self._rebuild_legend()

    def _legend_counts(self) -> List[Tuple[str, int]]:
        samples = sorted(self._value_by_sample.keys(), key=_natural_sort_key)
        if self._sample_filter is not None:
            samples = [sample for sample in samples if sample in self._sample_filter]
        counts = Counter(self._value_by_sample.get(sample, self._missing_label) for sample in samples)
        ordered_values = sorted(
            counts.keys(),
            key=lambda value: (value == self._missing_label, _natural_sort_key(value)),
        )
        if not self._show_missing:
            ordered_values = [value for value in ordered_values if value != self._missing_label]
        return [(value, int(counts[value])) for value in ordered_values if int(counts[value]) > 0]

    def _display_label(self, value: str) -> str:
        if value == self._missing_label:
            return value
        if self._value_prefix:
            return f"{self._value_prefix}{value}"
        return value

    def _rebuild_legend(self) -> None:
        while self._legend_layout.count():
            item = self._legend_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
        counts = self._legend_counts()
        columns = self._grid_columns if not self._scrollable else 1
        if not counts:
            empty = QtWidgets.QLabel(f"No {self._title.lower()} data", self._legend_wrap)
            empty.setStyleSheet("color: #94a3b8;")
            if self._scrollable:
                self._legend_layout.addWidget(empty)
                self._legend_layout.addStretch(1)
            else:
                self._legend_layout.addWidget(empty, 0, 0, 1, columns)
                cast_layout = self._legend_layout
                cast_layout.setRowStretch(1, 1)
            return
        for idx, (value, count) in enumerate(counts):
            color = self._color_by_value.get(value, QtGui.QColor("#e2e8f0"))
            row_widget = _LegendSwatchRow(self._display_label(value), QtGui.QColor(color), count, self._legend_wrap)
            if self._scrollable:
                self._legend_layout.addWidget(row_widget)
            else:
                row = idx // columns
                col = idx % columns
                self._legend_layout.addWidget(
                    row_widget,
                    row,
                    col,
                    alignment=QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignTop,
                )
        if self._scrollable:
            self._legend_layout.addStretch(1)
        else:
            cast_layout = self._legend_layout
            cast_layout.setRowStretch((len(counts) + columns - 1) // columns, 1)


class _StaticLegendPanel(QtWidgets.QWidget):
    def __init__(
        self,
        title: str,
        items: Sequence[Tuple[str, object]],
        *,
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self._title = str(title).strip() or "Legend"
        self._items = [(str(label), QtGui.QColor(color)) for label, color in items if QtGui.QColor(color).isValid()]

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(2, 0, 0, 0)
        layout.setSpacing(4)

        title_label = QtWidgets.QLabel(self._title, self)
        title_label.setStyleSheet("color: #334155; font-weight: 700; padding-left: 4px;")
        layout.addWidget(title_label)

        wrap = QtWidgets.QWidget(self)
        grid = QtWidgets.QGridLayout(wrap)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(8)
        grid.setVerticalSpacing(2)
        for idx, (label, color) in enumerate(self._items):
            row_widget = _LegendSwatchRow(label, color, None, wrap)
            grid.addWidget(
                row_widget,
                idx // 2,
                idx % 2,
                alignment=QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignTop,
            )
        grid.setRowStretch((len(self._items) + 1) // 2, 1)
        layout.addWidget(wrap, 0, QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignTop)

class _CategoricalMetadataStripCanvas(QtWidgets.QWidget):
    def __init__(self, panel: "_CategoricalMetadataStripPanel", parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self._panel = panel
        self.setMinimumWidth(84)
        self.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding,
        )

    def paintEvent(self, event: Optional[QtGui.QPaintEvent]) -> None:
        _ = event
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        painter.fillRect(self.rect(), QtGui.QColor("#ffffff"))

        tree = self._panel.tree_view()
        order = tree.leaf_order() if tree is not None else []
        if tree is None or not order:
            painter.setPen(QtGui.QColor("#94a3b8"))
            painter.drawText(self.rect().adjusted(8, 8, -8, -8), QtCore.Qt.AlignmentFlag.AlignCenter, "No tree")
            return
        if not self._panel.has_metadata():
            painter.setPen(QtGui.QColor("#94a3b8"))
            painter.drawText(self.rect().adjusted(8, 8, -8, -8), QtCore.Qt.AlignmentFlag.AlignCenter, "No K locus")
            return

        y_positions: List[Tuple[str, float]] = []
        for sample in order:
            y = tree.leaf_viewport_y(sample)
            if y is None:
                continue
            y_positions.append((sample, float(y)))
        if not y_positions:
            painter.setPen(QtGui.QColor("#94a3b8"))
            painter.drawText(self.rect().adjusted(8, 8, -8, -8), QtCore.Qt.AlignmentFlag.AlignCenter, "No view")
            return

        left = 12.0
        width = max(24.0, float(self.width()) - 24.0)
        height = float(self.height())
        active_samples = self._panel.active_samples()
        highlight_sample = self._panel.highlight_sample()

        for idx, (sample, y_center) in enumerate(y_positions):
            prev_y = y_positions[idx - 1][1] if idx > 0 else y_center
            next_y = y_positions[idx + 1][1] if idx + 1 < len(y_positions) else y_center
            top = y_center - ((y_center - prev_y) * 0.5 if idx > 0 else 9.0)
            bottom = y_center + ((next_y - y_center) * 0.5 if idx + 1 < len(y_positions) else 9.0)
            rect_top = max(2.0, top + 1.5)
            rect_bottom = min(height - 2.0, bottom - 1.5)
            rect_height = rect_bottom - rect_top
            if rect_height <= 1.0:
                continue

            value = self._panel.sample_value(sample)
            color = QtGui.QColor(self._panel.color_for_value(value))
            if active_samples is not None and sample not in active_samples:
                color = QtGui.QColor("#e2e8f0")
            elif value == self._panel.missing_label():
                color.setAlpha(210)

            rect = QtCore.QRectF(left, rect_top, width, rect_height)
            painter.setPen(QtCore.Qt.PenStyle.NoPen)
            painter.setBrush(QtGui.QBrush(color))
            painter.drawRoundedRect(rect, 4.0, 4.0)

            if highlight_sample and sample == highlight_sample:
                outline = QtGui.QPen(QtGui.QColor("#f97316"))
                outline.setWidthF(2.0)
                painter.setPen(outline)
                painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
                painter.drawRoundedRect(rect.adjusted(0.5, 0.5, -0.5, -0.5), 4.0, 4.0)


class _CategoricalMetadataStripPanel(QtWidgets.QWidget):
    def __init__(
        self,
        meta_df: pd.DataFrame,
        *,
        title: str,
        column_candidates: Sequence[str],
        show_canvas: bool = True,
        legend_scrollable: bool = True,
        legend_columns: int = 1,
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self._meta_df = meta_df.copy()
        self._title = str(title).strip() or "Metadata"
        self._column_candidates = tuple(str(c) for c in column_candidates if str(c).strip())
        self._show_canvas = bool(show_canvas)
        self._legend_scrollable = bool(legend_scrollable)
        self._legend_columns = max(1, int(legend_columns))
        self._sample_filter: Optional[set[str]] = None
        self._highlight_sample: Optional[str] = None
        self._tree_view: Optional[RaxmlTreeView] = None
        self._value_by_sample: Dict[str, str] = {}
        self._color_by_value: Dict[str, QtGui.QColor] = {}
        self._resolved_column: Optional[str] = None
        self._missing_label = "Missing"

        if self._show_canvas:
            self.setMinimumWidth(210)
            self.setMaximumWidth(260)

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(2, 0, 0, 0)
        layout.setSpacing(4)

        title_label = QtWidgets.QLabel(self._title, self)
        title_label.setStyleSheet("color: #334155; font-weight: 700; padding-left: 4px;")
        layout.addWidget(title_label)

        self._canvas: Optional[_CategoricalMetadataStripCanvas] = None
        if self._show_canvas:
            self._canvas = _CategoricalMetadataStripCanvas(self, self)
            layout.addWidget(self._canvas, stretch=1)

        self._legend_scroll: Optional[QtWidgets.QScrollArea] = None
        self._legend_wrap = QtWidgets.QWidget(self)
        if self._legend_scrollable:
            self._legend_scroll = QtWidgets.QScrollArea(self)
            self._legend_scroll.setWidgetResizable(True)
            self._legend_scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
            self._legend_scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
            self._legend_scroll.setMinimumHeight(110)
            self._legend_scroll.setMaximumHeight(140)
            self._legend_layout = QtWidgets.QVBoxLayout(self._legend_wrap)
            self._legend_layout.setContentsMargins(0, 0, 0, 0)
            self._legend_layout.setSpacing(4)
            self._legend_scroll.setWidget(self._legend_wrap)
            layout.addWidget(self._legend_scroll, stretch=0)
        else:
            self._legend_layout = QtWidgets.QGridLayout(self._legend_wrap)
            self._legend_layout.setContentsMargins(0, 0, 0, 0)
            self._legend_layout.setHorizontalSpacing(8)
            self._legend_layout.setVerticalSpacing(2)
            layout.addWidget(self._legend_wrap, stretch=0, alignment=QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignTop)

        self.refresh()

    @staticmethod
    def _resolve_column(columns: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
        available = {str(col): str(col) for col in columns}
        lowered = {str(col).casefold(): str(col) for col in columns}
        for candidate in candidates:
            if candidate in available:
                return available[candidate]
            match = lowered.get(str(candidate).casefold())
            if match:
                return match
        return None

    def tree_view(self) -> Optional[RaxmlTreeView]:
        return self._tree_view

    def has_metadata(self) -> bool:
        return bool(self._resolved_column and self._value_by_sample)

    def missing_label(self) -> str:
        return self._missing_label

    def active_samples(self) -> Optional[set[str]]:
        return self._sample_filter

    def highlight_sample(self) -> Optional[str]:
        return self._highlight_sample

    def sample_value(self, sample: str) -> str:
        value = self._value_by_sample.get(str(sample).strip(), self._missing_label)
        return value or self._missing_label

    def color_for_value(self, value: str) -> QtGui.QColor:
        return QtGui.QColor(self._color_by_value.get(value, QtGui.QColor("#e2e8f0")))

    def title(self) -> str:
        return self._title

    def metadata_values(self) -> Dict[str, str]:
        return dict(self._value_by_sample)

    def metadata_colors(self) -> Dict[str, QtGui.QColor]:
        return {key: QtGui.QColor(color) for key, color in self._color_by_value.items()}

    def set_tree_view(self, tree: Optional[RaxmlTreeView]) -> None:
        if self._tree_view is tree:
            return
        if self._tree_view is not None and self._canvas is not None:
            try:
                self._tree_view.viewChanged.disconnect(self._canvas.update)
            except Exception:
                pass
        self._tree_view = tree
        if self._tree_view is not None and self._canvas is not None:
            self._tree_view.viewChanged.connect(self._canvas.update)
        if self._canvas is not None:
            self._canvas.update()
        self._rebuild_legend()

    def set_sample_filter(self, samples: Optional[Sequence[str]]) -> None:
        next_filter = {str(s).strip() for s in samples if str(s).strip()} if samples is not None else None
        if next_filter == self._sample_filter:
            return
        self._sample_filter = next_filter
        if self._canvas is not None:
            self._canvas.update()
        self._rebuild_legend()

    def set_highlight_sample(self, sample: Optional[str]) -> None:
        value = str(sample).strip() if sample is not None else ""
        next_sample = value or None
        if next_sample == self._highlight_sample:
            return
        self._highlight_sample = next_sample
        if self._canvas is not None:
            self._canvas.update()

    def _normalize_value(self, value: object) -> str:
        try:
            if pd.isna(value):
                return self._missing_label
        except Exception:
            pass
        text = str(value).strip()
        if not text:
            return self._missing_label
        if text.casefold() in {"nan", "none", "null", "na", "n/a", "<na>", "missing", "unknown"}:
            return self._missing_label
        return text

    def _rebuild_metadata(self) -> None:
        self._value_by_sample = {}
        self._color_by_value = {}
        self._resolved_column = None
        if self._meta_df.empty or "Sample" not in self._meta_df.columns:
            return
        self._resolved_column = self._resolve_column(list(self._meta_df.columns), self._column_candidates)
        if not self._resolved_column:
            return
        seen_samples: set[str] = set()
        for _, row in self._meta_df.iterrows():
            sample = str(row.get("Sample") or "").strip()
            if not sample or sample in seen_samples:
                continue
            seen_samples.add(sample)
            self._value_by_sample[sample] = self._normalize_value(row.get(self._resolved_column))

        values = sorted(
            {value for value in self._value_by_sample.values() if value != self._missing_label},
            key=lambda v: str(v).casefold(),
        )
        for idx, value in enumerate(values):
            self._color_by_value[value] = _high_contrast_color(idx)
        self._color_by_value[self._missing_label] = QtGui.QColor("#e2e8f0")

    def _legend_counts(self) -> List[Tuple[str, int]]:
        samples: List[str]
        if self._tree_view is not None and self._tree_view.leaf_order():
            samples = list(self._tree_view.leaf_order())
        else:
            samples = sorted(self._value_by_sample.keys(), key=str.casefold)
        if self._sample_filter is not None:
            samples = [sample for sample in samples if sample in self._sample_filter]
        counts = Counter(self.sample_value(sample) for sample in samples)
        ordered_values = sorted(
            counts.keys(),
            key=lambda value: (value == self._missing_label, str(value).casefold()),
        )
        return [(value, int(counts[value])) for value in ordered_values]

    def _rebuild_legend(self) -> None:
        while self._legend_layout.count():
            item = self._legend_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
        columns = self._legend_columns if not self._legend_scrollable else 1
        if not self.has_metadata():
            empty = QtWidgets.QLabel("No K locus metadata", self._legend_wrap)
            empty.setStyleSheet("color: #94a3b8;")
            if self._legend_scrollable:
                self._legend_layout.addWidget(empty)
                self._legend_layout.addStretch(1)
            else:
                self._legend_layout.addWidget(empty, 0, 0, 1, columns)
                cast_layout = self._legend_layout
                cast_layout.setRowStretch(1, 1)
            return
        for idx, (value, count) in enumerate(self._legend_counts()):
            row_widget = _LegendSwatchRow(value, self.color_for_value(value), count, self._legend_wrap)
            if self._legend_scrollable:
                self._legend_layout.addWidget(row_widget)
            else:
                row = idx // columns
                col = idx % columns
                self._legend_layout.addWidget(
                    row_widget,
                    row,
                    col,
                    alignment=QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignTop,
                )
        if self._legend_scrollable:
            self._legend_layout.addStretch(1)
        else:
            cast_layout = self._legend_layout
            cast_layout.setRowStretch((len(self._legend_counts()) + columns - 1) // columns, 1)

    def refresh(self) -> None:
        self._rebuild_metadata()
        self._rebuild_legend()
        if self._canvas is not None:
            self._canvas.update()


class _OverlayContentPanel(QtWidgets.QWidget):
    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        _apply_analysis_panel_style(self)
        self._overlay_widget: Optional[QtWidgets.QWidget] = None
        self._overlay_margin = 12
        self._overlay_right_offset = 15
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self._content_layout = layout

    def set_content_widget(self, widget: QtWidgets.QWidget) -> None:
        self._content_layout.addWidget(widget, stretch=1)

    def set_overlay_widget(self, widget: QtWidgets.QWidget) -> None:
        self._overlay_widget = widget
        widget.setParent(self)
        widget.raise_()
        widget.show()
        self._reposition_overlay()

    def _reposition_overlay(self) -> None:
        if self._overlay_widget is None:
            return
        hint = self._overlay_widget.sizeHint()
        width = max(hint.width(), self._overlay_widget.minimumSizeHint().width())
        height = max(hint.height(), self._overlay_widget.minimumSizeHint().height())
        x = max(
            self._overlay_margin,
            self.width() - width - self._overlay_margin - self._overlay_right_offset,
        )
        y = self._overlay_margin
        self._overlay_widget.setGeometry(x, y, width, height)
        self._overlay_widget.raise_()

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        super().resizeEvent(event)
        self._reposition_overlay()


class _OverlaySelectorFrame(QtWidgets.QFrame):
    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setObjectName("NeoPanOverlaySelector")
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_StyledBackground, True)
        self.setStyleSheet(
            "#NeoPanOverlaySelector {"
            "  background: rgba(255, 255, 255, 224);"
            "  border: 1px solid rgba(215, 224, 234, 230);"
            "  border-radius: 8px;"
            "}"
        )
        self._idle_opacity = 0.58
        self._active_opacity = 0.97
        self._opacity_effect = QtWidgets.QGraphicsOpacityEffect(self)
        self.setGraphicsEffect(self._opacity_effect)
        self._set_active(False)

    def register_hover_widget(self, widget: QtWidgets.QWidget) -> None:
        widget.installEventFilter(self)
        for child in widget.findChildren(QtWidgets.QWidget):
            child.installEventFilter(self)

    def _set_active(self, active: bool) -> None:
        self._opacity_effect.setOpacity(self._active_opacity if active else self._idle_opacity)

    def _sync_from_cursor(self) -> None:
        local_pos = self.mapFromGlobal(QtGui.QCursor.pos())
        self._set_active(self.rect().contains(local_pos))

    def enterEvent(self, event: QtCore.QEvent) -> None:
        self._set_active(True)
        super().enterEvent(event)

    def leaveEvent(self, event: QtCore.QEvent) -> None:
        super().leaveEvent(event)
        QtCore.QTimer.singleShot(0, self._sync_from_cursor)

    def eventFilter(self, obj: QtCore.QObject, event: QtCore.QEvent) -> bool:
        event_type = event.type()
        if event_type in (QtCore.QEvent.Type.Enter, QtCore.QEvent.Type.HoverEnter):
            self._set_active(True)
        elif event_type in (QtCore.QEvent.Type.Leave, QtCore.QEvent.Type.HoverLeave):
            QtCore.QTimer.singleShot(0, self._sync_from_cursor)
        return super().eventFilter(obj, event)


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
