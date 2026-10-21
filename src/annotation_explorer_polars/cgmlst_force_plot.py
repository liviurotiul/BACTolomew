from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from PyQt6 import QtCore, QtGui, QtWidgets

from .panaroo_force import CgmlstForceGraphView, GraphKeyEntry
from .panaroo_style import _log_view_event
from .panaroo_widgets import _high_contrast_color
from .settings_column import (
    SettingsColumn,
    SettingsStrip,
    auto_fold,
    settings_card,
    settings_control,
    stacked_field,
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


class CgmlstForcePlotPanel(QtWidgets.QWidget):
    sampleHovered = QtCore.pyqtSignal(object)
    # The "Color nodes by" or "Shape nodes by" choice changed.
    styleChanged = QtCore.pyqtSignal()

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
        self._sample_filter: Optional[set[str]] = None
        # The (color, shape) traits styleChanged last reported.
        self._style_traits: Tuple[str, str] = ("", "")
        self._force_view = CgmlstForceGraphView(
            output_dir,
            reference_sample=reference_sample,
            parent=self,
        )
        self._force_view.sampleHovered.connect(self.sampleHovered.emit)

        # The Settings column (or its folded strip) beside the graph, as the
        # Phylogeny's: it folds while the panel is too narrow for both and
        # opens again when the panel grows, unless the user folded it.
        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        strip = SettingsStrip(self)
        strip.expandRequested.connect(lambda: self.set_settings_open(True))
        layout.addWidget(strip)
        self._settings_strip = strip

        column = SettingsColumn(self, object_name="FilterRail", max_width=300)
        column.collapseRequested.connect(lambda: self.set_settings_open(False))
        card, card_layout = settings_card(column.content, "Nodes")
        self._color_combo = settings_control(QtWidgets.QComboBox(card))
        self._color_combo.currentIndexChanged.connect(self._apply_style_controls)
        stacked_field(card, card_layout, "Color nodes by", self._color_combo)
        self._shape_combo = settings_control(QtWidgets.QComboBox(card))
        self._shape_combo.currentIndexChanged.connect(self._apply_style_controls)
        stacked_field(card, card_layout, "Shape nodes by", self._shape_combo)
        column.add_card(card)
        column.finish(min_width=248)
        column.setFixedWidth(column.maximumWidth())
        layout.addWidget(column)
        self._settings_column = column
        layout.addWidget(self._force_view, stretch=1)

        # Folded until the panel is shown at its size (``_fit_settings``).
        column.hide()
        self._settings_auto_folded = True
        # Opened by the user in a narrow panel: it stays open there until the
        # panel is wide enough for both (as a user's fold stands).
        self._settings_user_opened = False
        self._settings_fit_timer = QtCore.QTimer(self)
        self._settings_fit_timer.setSingleShot(True)
        self._settings_fit_timer.setInterval(0)
        self._settings_fit_timer.timeout.connect(self._fit_settings)

        self._rebuild_trait_combos()

    _GRAPH_MIN_WIDTH = 420

    def settings_open(self) -> bool:
        """True while the Settings column is open, False while folded."""
        return not self._settings_column.isHidden()

    def set_settings_open(self, visible: bool) -> None:
        """Open or fold the Settings column (the header chevron or the
        strip); the choice stands while the panel keeps its size."""
        self._settings_auto_folded = False
        self._settings_user_opened = bool(visible)
        self._show_settings(bool(visible))

    def _show_settings(self, visible: bool) -> None:
        self._settings_column.setVisible(visible)
        self._settings_strip.setVisible(not visible)

    def _fit_settings(self) -> None:
        if not self.isVisible():
            return
        small = self.width() < self._settings_column.width() + self._GRAPH_MIN_WIDTH
        if not small:
            self._settings_user_opened = False
        elif self._settings_user_opened:
            return
        self._settings_auto_folded = auto_fold(
            small, self.settings_open(), self._settings_auto_folded, self._show_settings
        )

    def showEvent(self, event: QtGui.QShowEvent) -> None:  # type: ignore[override]
        super().showEvent(event)
        self._settings_fit_timer.start()

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:  # type: ignore[override]
        super().resizeEvent(event)
        # Once the layout has settled at the new size.
        self._settings_fit_timer.start()

    def force_view(self) -> CgmlstForceGraphView:
        return self._force_view

    def set_reference_sample(self, sample: Optional[str], *, refresh: bool = True) -> None:
        self._force_view.set_reference_sample(sample, refresh=refresh)

    def set_sample_filter(self, samples: Optional[Sequence[str]]) -> None:
        self._force_view.set_sample_filter(samples, refresh=False)
        self._sample_filter = {str(s).strip() for s in samples} if samples is not None else None
        color_trait = str(self._color_combo.currentData() or "").strip()
        shape_trait = str(self._shape_combo.currentData() or "").strip()
        self._force_view.set_key_sections(self._key_sections(color_trait, shape_trait))

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

    def color_trait(self) -> str:
        """The trait the nodes are colored by ("" for none)."""
        return str(self._color_combo.currentData() or "").strip()

    def shape_trait(self) -> str:
        """The trait the node shapes follow ("" for none)."""
        return str(self._shape_combo.currentData() or "").strip()

    def set_style_traits(self, color: Optional[str] = None, shape: Optional[str] = None) -> None:
        """Color and shape the nodes by these traits (None keeps the current
        one; a trait the run does not have is ignored)."""
        names = self._trait_names()
        with QtCore.QSignalBlocker(self._color_combo), QtCore.QSignalBlocker(self._shape_combo):
            if color is not None and color in names:
                self._set_combo_value(self._color_combo, color)
            if shape is not None and (shape in names or shape == ""):
                self._set_combo_value(self._shape_combo, shape)
        self._apply_style_controls()

    def _apply_style_controls(self) -> None:
        color_trait = str(self._color_combo.currentData() or "").strip()
        shape_trait = str(self._shape_combo.currentData() or "").strip()
        color_map = self._color_map_for_trait(color_trait) if color_trait else None
        symbol_map = self._symbol_map_for_trait(shape_trait) if shape_trait else None
        self._force_view.set_node_style_maps(
            color_by_sample=color_map,
            symbol_by_sample=symbol_map,
        )
        self._force_view.set_key_sections(self._key_sections(color_trait, shape_trait))
        if (color_trait, shape_trait) != self._style_traits:
            self._style_traits = (color_trait, shape_trait)
            self.styleChanged.emit()

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
        """Samples per value, among the window's samples."""
        counts: Dict[str, int] = {}
        for sample, value in self._trait_values_by_name.get(trait_name, {}).items():
            if self._sample_filter is not None and sample not in self._sample_filter:
                continue
            text = str(value).strip()
            if not text:
                continue
            counts[text] = counts.get(text, 0) + 1
        return counts

    def _key_sections(self, color_trait: str, shape_trait: str) -> List[Tuple[str, List[GraphKeyEntry]]]:
        """The graph key: the color and the shape of each value, with counts."""
        sections: List[Tuple[str, List[GraphKeyEntry]]] = []
        if color_trait:
            counts = self._trait_counts(color_trait)
            sections.append((
                f"Color: {color_trait}",
                [
                    GraphKeyEntry(value, counts.get(value, 0), color=self._color_by_value[value])
                    for value in sorted(self._color_by_value, key=str.casefold)
                ],
            ))
        if shape_trait:
            counts = self._trait_counts(shape_trait)
            sections.append((
                f"Shape: {shape_trait}",
                [
                    GraphKeyEntry(value, counts.get(value, 0), symbol=self._symbol_by_value[value])
                    for value in sorted(self._symbol_by_value, key=str.casefold)
                ],
            ))
        return sections
