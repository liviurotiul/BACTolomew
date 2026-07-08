from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence

from PyQt6 import QtCore, QtGui, QtWidgets

from .panaroo_force import CgmlstForceGraphView
from .panaroo_style import (
    _apply_analysis_group_style,
    _apply_analysis_panel_style,
    _log_view_event,
)
from .panaroo_widgets import (
    _LegendSwatchRow,
    _high_contrast_color,
)

_FORCE_SYMBOLS = (
    "o",
    "s",
    "t",
    "d",
    "star",
    "+",
    "x",
    "p",
    "h",
    "t1",
    "t2",
    "t3",
)


class _ShapeLegendGlyph(QtWidgets.QWidget):
    def __init__(
        self,
        symbol: str,
        *,
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self._symbol = str(symbol or "").strip() or "o"
        self.setFixedSize(18, 18)

    def paintEvent(self, event: QtGui.QPaintEvent) -> None:
        del event
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        pen = QtGui.QPen(QtGui.QColor("#0f172a"))
        pen.setWidthF(1.5)
        painter.setPen(pen)
        painter.setBrush(QtGui.QBrush(QtGui.QColor("#e2e8f0")))

        rect = QtCore.QRectF(2.0, 2.0, float(self.width()) - 4.0, float(self.height()) - 4.0)
        cx = rect.center().x()
        cy = rect.center().y()
        symbol = self._symbol

        if symbol == "o":
            painter.drawEllipse(rect)
            return
        if symbol == "s":
            painter.drawRect(rect)
            return
        if symbol in {"t", "t1", "t2", "t3"}:
            if symbol == "t3":
                points = [
                    QtCore.QPointF(rect.left(), cy),
                    QtCore.QPointF(rect.right(), rect.top()),
                    QtCore.QPointF(rect.right(), rect.bottom()),
                ]
            elif symbol == "t2":
                points = [
                    QtCore.QPointF(rect.left(), rect.top()),
                    QtCore.QPointF(rect.left(), rect.bottom()),
                    QtCore.QPointF(rect.right(), cy),
                ]
            elif symbol == "t1":
                points = [
                    QtCore.QPointF(cx, rect.bottom()),
                    QtCore.QPointF(rect.left(), rect.top()),
                    QtCore.QPointF(rect.right(), rect.top()),
                ]
            else:
                points = [
                    QtCore.QPointF(cx, rect.top()),
                    QtCore.QPointF(rect.left(), rect.bottom()),
                    QtCore.QPointF(rect.right(), rect.bottom()),
                ]
            painter.drawPolygon(QtGui.QPolygonF(points))
            return
        if symbol == "d":
            painter.drawPolygon(
                QtGui.QPolygonF(
                    [
                        QtCore.QPointF(cx, rect.top()),
                        QtCore.QPointF(rect.left(), cy),
                        QtCore.QPointF(cx, rect.bottom()),
                        QtCore.QPointF(rect.right(), cy),
                    ]
                )
            )
            return
        if symbol == "p":
            painter.drawPolygon(self._regular_polygon(rect, 5, rotation=-90.0))
            return
        if symbol == "h":
            painter.drawPolygon(self._regular_polygon(rect, 6, rotation=-90.0))
            return
        if symbol == "star":
            painter.drawPolygon(self._star_polygon(rect))
            return
        if symbol == "+":
            painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
            painter.drawLine(QtCore.QPointF(cx, rect.top()), QtCore.QPointF(cx, rect.bottom()))
            painter.drawLine(QtCore.QPointF(rect.left(), cy), QtCore.QPointF(rect.right(), cy))
            return
        if symbol == "x":
            painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
            painter.drawLine(QtCore.QPointF(rect.left(), rect.top()), QtCore.QPointF(rect.right(), rect.bottom()))
            painter.drawLine(QtCore.QPointF(rect.left(), rect.bottom()), QtCore.QPointF(rect.right(), rect.top()))
            return
        painter.drawEllipse(rect)

    @staticmethod
    def _regular_polygon(rect: QtCore.QRectF, sides: int, *, rotation: float) -> QtGui.QPolygonF:
        points = []
        radius = min(rect.width(), rect.height()) * 0.5
        center = rect.center()
        for idx in range(sides):
            angle = math.radians(rotation + (360.0 * idx / sides))
            points.append(
                QtCore.QPointF(
                    center.x() + radius * math.cos(angle),
                    center.y() + radius * math.sin(angle),
                )
            )
        return QtGui.QPolygonF(points)

    @staticmethod
    def _star_polygon(rect: QtCore.QRectF) -> QtGui.QPolygonF:
        points = []
        outer = min(rect.width(), rect.height()) * 0.5
        inner = outer * 0.45
        center = rect.center()
        for idx in range(10):
            angle = math.radians(-90.0 + idx * 36.0)
            radius = outer if idx % 2 == 0 else inner
            points.append(
                QtCore.QPointF(
                    center.x() + radius * math.cos(angle),
                    center.y() + radius * math.sin(angle),
                )
            )
        return QtGui.QPolygonF(points)


class _ShapeLegendRow(QtWidgets.QWidget):
    def __init__(
        self,
        label: str,
        symbol: str,
        count: int,
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Maximum,
            QtWidgets.QSizePolicy.Policy.Preferred,
        )
        row = QtWidgets.QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(6)
        row.addWidget(_ShapeLegendGlyph(symbol, parent=self), 0, QtCore.Qt.AlignmentFlag.AlignTop)

        text = QtWidgets.QLabel(f"{label} ({int(count)})", self)
        text.setWordWrap(True)
        text.setStyleSheet("color: #334155;")
        text.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Maximum,
            QtWidgets.QSizePolicy.Policy.Preferred,
        )
        row.addWidget(text, 0, QtCore.Qt.AlignmentFlag.AlignTop)


class _ForceLegendSection(QtWidgets.QGroupBox):
    def __init__(self, title: str, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(title, parent)
        _apply_analysis_group_style(self)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(8, 10, 8, 8)
        layout.setSpacing(6)
        self._layout = layout

    def clear(self) -> None:
        while self._layout.count():
            item = self._layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()

    def add_widget(self, widget: QtWidgets.QWidget) -> None:
        self._layout.addWidget(widget)

    def add_stretch(self) -> None:
        self._layout.addStretch(1)


class CgmlstForcePlotPanel(QtWidgets.QWidget):
    sampleHovered = QtCore.pyqtSignal(object)

    def __init__(
        self,
        output_dir: Path,
        *,
        reference_sample: Optional[str] = None,
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        _log_view_event(f"created standalone cgMLST force plot for {Path(output_dir)}")
        self._trait_values_by_name: Dict[str, Dict[str, str]] = {}
        self._color_by_value: Dict[str, QtGui.QColor] = {}
        self._symbol_by_value: Dict[str, str] = {}
        self._force_view = CgmlstForceGraphView(
            output_dir,
            reference_sample=reference_sample,
            parent=self,
        )
        self._force_view.sampleHovered.connect(self.sampleHovered.emit)

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal, self)
        split.setChildrenCollapsible(False)
        split.setHandleWidth(8)
        layout.addWidget(split, stretch=1)

        controls_panel = QtWidgets.QWidget(split)
        _apply_analysis_panel_style(controls_panel)
        controls_panel.setMinimumWidth(260)
        controls_panel.setMaximumWidth(340)
        controls_panel_layout = QtWidgets.QVBoxLayout(controls_panel)
        controls_panel_layout.setContentsMargins(10, 10, 10, 10)
        controls_panel_layout.setSpacing(8)

        controls = QtWidgets.QGroupBox("Settings", controls_panel)
        _apply_analysis_group_style(controls)
        controls_layout = QtWidgets.QFormLayout(controls)
        controls_layout.setContentsMargins(10, 14, 10, 10)
        controls_layout.setSpacing(6)

        self._color_combo = QtWidgets.QComboBox(controls)
        self._color_combo.currentIndexChanged.connect(self._apply_style_controls)
        controls_layout.addRow("Color nodes by", self._color_combo)

        self._shape_combo = QtWidgets.QComboBox(controls)
        self._shape_combo.currentIndexChanged.connect(self._apply_style_controls)
        controls_layout.addRow("Shape nodes by", self._shape_combo)

        controls_panel_layout.addWidget(controls, stretch=0)
        controls_panel_layout.addStretch(1)

        legend_panel = QtWidgets.QWidget(split)
        _apply_analysis_panel_style(legend_panel)
        legend_panel.setMinimumWidth(280)
        legend_panel.setMaximumWidth(380)
        legend_panel_layout = QtWidgets.QVBoxLayout(legend_panel)
        legend_panel_layout.setContentsMargins(10, 10, 10, 10)
        legend_panel_layout.setSpacing(8)

        legend_group = QtWidgets.QGroupBox("Legend", legend_panel)
        _apply_analysis_group_style(legend_group)
        legend_group_layout = QtWidgets.QVBoxLayout(legend_group)
        legend_group_layout.setContentsMargins(10, 14, 10, 10)
        legend_group_layout.setSpacing(8)
        self._legend_hint = QtWidgets.QLabel("Choose a color or shape grouping to see the legend.")
        self._legend_hint.setWordWrap(True)
        self._legend_hint.setStyleSheet("color: #607183;")
        legend_group_layout.addWidget(self._legend_hint)
        self._legend_scroll = QtWidgets.QScrollArea(legend_group)
        self._legend_scroll.setWidgetResizable(True)
        self._legend_scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self._legend_scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self._legend_wrap = QtWidgets.QWidget(self._legend_scroll)
        self._legend_layout = QtWidgets.QVBoxLayout(self._legend_wrap)
        self._legend_layout.setContentsMargins(0, 0, 0, 0)
        self._legend_layout.setSpacing(8)
        self._legend_scroll.setWidget(self._legend_wrap)
        legend_group_layout.addWidget(self._legend_scroll, stretch=1)
        legend_panel_layout.addWidget(legend_group, stretch=1)

        split.addWidget(controls_panel)
        split.addWidget(self._force_view)
        split.addWidget(legend_panel)
        split.setStretchFactor(0, 0)
        split.setStretchFactor(1, 1)
        split.setStretchFactor(2, 0)
        split.setSizes([280, 1040, 320])

        self._rebuild_trait_combos()

    def force_view(self) -> CgmlstForceGraphView:
        return self._force_view

    def set_reference_sample(self, sample: Optional[str], *, refresh: bool = True) -> None:
        self._force_view.set_reference_sample(sample, refresh=refresh)

    def set_sample_filter(self, samples: Optional[Sequence[str]]) -> None:
        self._force_view.set_sample_filter(samples, refresh=False)

    def set_fastgear_assignments(
        self,
        lineage_by_sample: Mapping[str, str],
        lineage_colors: Mapping[str, QtGui.QColor],
        cluster_by_sample: Mapping[str, str],
    ) -> None:
        self._force_view.set_fastgear_assignments(
            lineage_by_sample,
            lineage_colors,
            cluster_by_sample,
        )

    def set_sample_hover_metadata(self, metadata: Optional[Mapping[str, Mapping[str, object]]]) -> None:
        self._force_view.set_sample_hover_metadata(metadata)

    def set_trait_values(self, trait_values_by_name: Optional[Mapping[str, Mapping[str, object]]]) -> None:
        normalized: Dict[str, Dict[str, str]] = {}
        for raw_name, raw_values in (trait_values_by_name or {}).items():
            name = str(raw_name).strip()
            if not name or not isinstance(raw_values, Mapping):
                continue
            values: Dict[str, str] = {}
            for raw_sample, raw_value in raw_values.items():
                sample = str(raw_sample).strip()
                value = str(raw_value).strip() if raw_value is not None else ""
                if not sample or not value:
                    continue
                values[sample] = value
            if values:
                normalized[name] = values
        self._trait_values_by_name = normalized
        self._rebuild_trait_combos()
        self._apply_style_controls()

    def refresh(self) -> None:
        self._force_view.refresh()

    def ensure_fresh(self) -> None:
        self._force_view.ensure_fresh()

    def highlight_sample(self, sample: Optional[str]) -> None:
        self._force_view.highlight_sample(sample)

    def highlight_samples(self, samples: Optional[Sequence[str]]) -> None:
        self._force_view.highlight_samples(samples)

    def _trait_names(self) -> list[str]:
        special = ["Clade", "Cluster"]
        ordered: list[str] = []
        seen: set[str] = set()
        for name in special:
            if name in self._trait_values_by_name and name not in seen:
                seen.add(name)
                ordered.append(name)
        for name in sorted(self._trait_values_by_name, key=str.casefold):
            if name in seen:
                continue
            seen.add(name)
            ordered.append(name)
        return ordered

    def _rebuild_trait_combos(self) -> None:
        color_value = str(self._color_combo.currentData() or "").strip()
        shape_value = str(self._shape_combo.currentData() or "").strip()
        names = self._trait_names()

        color_default = "Clade" if "Clade" in self._trait_values_by_name else (names[0] if names else "")
        shape_default = ""

        with QtCore.QSignalBlocker(self._color_combo), QtCore.QSignalBlocker(self._shape_combo):
            self._color_combo.clear()
            self._shape_combo.clear()
            for name in names:
                self._color_combo.addItem(name, name)
            self._shape_combo.addItem("None", "")
            for name in names:
                self._shape_combo.addItem(name, name)

            target_color = color_value if color_value in names else color_default
            target_shape = shape_value if shape_value in names else shape_default
            self._set_combo_value(self._color_combo, target_color)
            self._set_combo_value(self._shape_combo, target_shape)

    @staticmethod
    def _set_combo_value(combo: QtWidgets.QComboBox, value: str) -> None:
        for idx in range(combo.count()):
            if str(combo.itemData(idx) or "").strip() == value:
                combo.setCurrentIndex(idx)
                return
        if combo.count() > 0:
            combo.setCurrentIndex(0)

    def _apply_style_controls(self) -> None:
        color_trait = str(self._color_combo.currentData() or "").strip()
        shape_trait = str(self._shape_combo.currentData() or "").strip()
        color_map = self._color_map_for_trait(color_trait) if color_trait else None
        symbol_map = self._symbol_map_for_trait(shape_trait) if shape_trait else None
        self._force_view.set_node_style_maps(
            color_by_sample=color_map,
            symbol_by_sample=symbol_map,
        )
        self._rebuild_legend_panel(color_trait, shape_trait)

    def _color_map_for_trait(self, trait_name: str) -> Dict[str, QtGui.QColor]:
        values_by_sample = self._trait_values_by_name.get(trait_name, {})
        unique_values = sorted({value for value in values_by_sample.values() if value}, key=str.casefold)
        palette = {
            value: _high_contrast_color(idx)
            for idx, value in enumerate(unique_values)
        }
        self._color_by_value = {value: QtGui.QColor(color) for value, color in palette.items()}
        return {
            sample: QtGui.QColor(palette[value])
            for sample, value in values_by_sample.items()
            if value in palette
        }

    def _symbol_map_for_trait(self, trait_name: str) -> Dict[str, str]:
        values_by_sample = self._trait_values_by_name.get(trait_name, {})
        unique_values = sorted({value for value in values_by_sample.values() if value}, key=str.casefold)
        symbol_by_value = {
            value: _FORCE_SYMBOLS[idx % len(_FORCE_SYMBOLS)]
            for idx, value in enumerate(unique_values)
        }
        self._symbol_by_value = dict(symbol_by_value)
        return {
            sample: symbol_by_value[value]
            for sample, value in values_by_sample.items()
            if value in symbol_by_value
        }

    def _trait_counts(self, trait_name: str) -> Dict[str, int]:
        counts: Dict[str, int] = {}
        for value in self._trait_values_by_name.get(trait_name, {}).values():
            text = str(value).strip()
            if not text:
                continue
            counts[text] = counts.get(text, 0) + 1
        return counts

    def _clear_legend_layout(self) -> None:
        while self._legend_layout.count():
            item = self._legend_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()

    def _rebuild_legend_panel(self, color_trait: str, shape_trait: str) -> None:
        self._clear_legend_layout()
        has_sections = False

        if color_trait:
            color_section = _ForceLegendSection(f"Colors: {color_trait}", self._legend_wrap)
            counts = self._trait_counts(color_trait)
            for value in sorted(self._color_by_value, key=str.casefold):
                color_section.add_widget(
                    _LegendSwatchRow(
                        value,
                        self._color_by_value[value],
                        counts.get(value, 0),
                        parent=color_section,
                    )
                )
            color_section.add_stretch()
            self._legend_layout.addWidget(color_section)
            has_sections = True

        if shape_trait:
            shape_section = _ForceLegendSection(f"Shapes: {shape_trait}", self._legend_wrap)
            counts = self._trait_counts(shape_trait)
            for value in sorted(self._symbol_by_value, key=str.casefold):
                shape_section.add_widget(
                    _ShapeLegendRow(
                        value,
                        self._symbol_by_value[value],
                        counts.get(value, 0),
                        parent=shape_section,
                    )
                )
            shape_section.add_stretch()
            self._legend_layout.addWidget(shape_section)
            has_sections = True

        self._legend_hint.setVisible(not has_sections)
        self._legend_scroll.setVisible(has_sections)
        self._legend_layout.addStretch(1)
