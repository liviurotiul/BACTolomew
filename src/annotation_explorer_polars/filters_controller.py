from __future__ import annotations

import bisect
import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import polars as pl
from PyQt6 import QtCore, QtGui, QtWidgets

from .alignment_dialog import AlignmentDialog
from .annotation_detail_dialog import AnnotationDetailDialog
from .annotation_table import AnnotationTableModel, DetailsDelegate, TagsDelegate
from .constants import ANNOTATION_COLUMNS_TO_SHOW, MISSING_FILTER_TOKEN
from .geo_view import GeoMapView
from .igv_view import GeneSegment, IgvTrackView
from .metadata_filter_dialog import MetadataValueDialog
from .panaroo_common import (
    _NEO_TREE_MGE_TOKEN,
    _NEO_TREE_PANAROO_TOKEN,
    _NEO_TREE_PLASMID_TOKEN,
    _NEO_TREE_VIRULENCE_TOKEN,
)
from .panaroo_data import _extract_mge_labels, _extract_virulence_labels, _find_panaroo_out
from .panaroo_style import _apply_analysis_group_style
from .panaroo_views import (
    EntropyTrackTab,
    PanCoreAccessoryTab,
    SimilarityMapTab,
    StructPresenceHeatmapTab,
)

MISSING_FILTER_LABEL = "(missing)"

_NEO_TREE_EXCLUDED_TRAIT_KEYS = {
    "sample",
    "latitude",
    "longitude",
    "lat",
    "lon",
    "long",
    "fastgearlineage",
    "fastgearcluster",
    "clade",
    "cluster",
}

NEO_TREE_RESISTANCE_TOKEN = "__neo_tree_resistance__"
NEO_TREE_RESISTANCE_LABEL = "AMR"
NEO_TREE_MGE_TOKEN = _NEO_TREE_MGE_TOKEN
NEO_TREE_MGE_LABEL = "MGEs"
NEO_TREE_VIRULENCE_TOKEN = _NEO_TREE_VIRULENCE_TOKEN
NEO_TREE_VIRULENCE_LABEL = "Virulence"
NEO_TREE_PANAROO_TOKEN = _NEO_TREE_PANAROO_TOKEN
NEO_TREE_PANAROO_LABEL = "Panaroo gene presence/absence"
NEO_TREE_PLASMID_TOKEN = _NEO_TREE_PLASMID_TOKEN
NEO_TREE_PLASMID_LABEL = "Plasmids / MOB clusters"

_NEO_TREE_PREFERRED_TRAIT_CANDIDATES: Tuple[Tuple[str, ...], ...] = (
    ("K locus", "K_locus", "KLocus", "KL", "KL type", "K type"),
    ("Country",),
)

TAG_OPTIONS: List[Tuple[str, str]] = [
    ("A", "AMR"),
    ("B", "BMRG"),
    ("V", "Virulence"),
    ("M", "Mobile"),
    ("P", "Phage"),
    ("I", "Compound overlap"),
]


class FiltersControllerMixin:
    @staticmethod
    def _neo_tree_trait_present_color_text(color: QtGui.QColor) -> str:
        return str(color.name()).upper()

    @staticmethod
    def _neo_tree_trait_button_fg(color: QtGui.QColor) -> str:
        luminance = (0.299 * color.red()) + (0.587 * color.green()) + (0.114 * color.blue())
        return "#0f172a" if luminance >= 170 else "#f8fafc"

    @staticmethod
    def _neo_tree_trait_color_value(color: object, fallback: str) -> QtGui.QColor:
        if isinstance(color, QtGui.QColor):
            resolved = QtGui.QColor(color)
        else:
            resolved = QtGui.QColor(str(color or ""))
        if not resolved.isValid():
            resolved = QtGui.QColor(fallback)
        return resolved

    def _neo_tree_trait_present_color(self, token: str) -> QtGui.QColor:
        if str(token) == NEO_TREE_RESISTANCE_TOKEN:
            return self._neo_tree_trait_color_value(getattr(self, "_neo_tree_amr_color", None), "#f97316")
        if str(token) == NEO_TREE_MGE_TOKEN:
            return self._neo_tree_trait_color_value(getattr(self, "_neo_tree_mge_color", None), "#0f766e")
        if str(token) == NEO_TREE_VIRULENCE_TOKEN:
            return self._neo_tree_trait_color_value(getattr(self, "_neo_tree_virulence_color", None), "#be185d")
        return QtGui.QColor("#475569")

    def _set_neo_tree_trait_present_color(self, token: str, color: QtGui.QColor) -> None:
        resolved = self._neo_tree_trait_color_value(color, "#475569")
        if str(token) == NEO_TREE_RESISTANCE_TOKEN:
            self._neo_tree_amr_color = resolved
        elif str(token) == NEO_TREE_MGE_TOKEN:
            self._neo_tree_mge_color = resolved
        elif str(token) == NEO_TREE_VIRULENCE_TOKEN:
            self._neo_tree_virulence_color = resolved
        else:
            return
        self._update_neo_tree_trait_color_controls()
        self._apply_neo_tree_trait_color_selection()

    def _apply_neo_tree_trait_color_button_style(
        self,
        button: Optional[QtWidgets.QPushButton],
        color: QtGui.QColor,
        tooltip: str,
    ) -> None:
        if button is None:
            return
        text_color = self._neo_tree_trait_button_fg(color)
        button.setText(self._neo_tree_trait_present_color_text(color))
        button.setToolTip(tooltip)
        button.setStyleSheet(
            "QPushButton {"
            f"background: {color.name()};"
            f"color: {text_color};"
            "border: 1px solid #cbd5e1;"
            "border-radius: 8px;"
            "padding: 4px 10px;"
            "font-weight: 700;"
            "}"
        )

    def _update_neo_tree_trait_color_controls(self) -> None:
        available = set(getattr(self, "_neo_tree_trait_available", []) or [])
        group = getattr(self, "neo_tree_trait_colors_group", None)
        amr_row = getattr(self, "neo_tree_amr_color_row", None)
        mge_row = getattr(self, "neo_tree_mge_color_row", None)
        virulence_row = getattr(self, "neo_tree_virulence_color_row", None)
        show_amr = NEO_TREE_RESISTANCE_TOKEN in available
        show_mge = NEO_TREE_MGE_TOKEN in available
        show_virulence = NEO_TREE_VIRULENCE_TOKEN in available
        if amr_row is not None:
            amr_row.setVisible(show_amr)
        if mge_row is not None:
            mge_row.setVisible(show_mge)
        if virulence_row is not None:
            virulence_row.setVisible(show_virulence)
        if group is not None:
            group.setVisible(show_amr or show_mge or show_virulence)
            if isinstance(group, QtWidgets.QGroupBox):
                _apply_analysis_group_style(group)
        self._apply_neo_tree_trait_color_button_style(
            getattr(self, "neo_tree_amr_color_button", None),
            self._neo_tree_trait_present_color(NEO_TREE_RESISTANCE_TOKEN),
            "Set the AMR present color used in the Phandango-style plot.",
        )
        self._apply_neo_tree_trait_color_button_style(
            getattr(self, "neo_tree_mge_color_button", None),
            self._neo_tree_trait_present_color(NEO_TREE_MGE_TOKEN),
            "Set the MGE present color used in the Phandango-style plot.",
        )
        self._apply_neo_tree_trait_color_button_style(
            getattr(self, "neo_tree_virulence_color_button", None),
            self._neo_tree_trait_present_color(NEO_TREE_VIRULENCE_TOKEN),
            "Set the virulence present color used in the Phandango-style plot.",
        )

    def _connect_neo_tree_trait_color_controls(self) -> None:
        if getattr(self, "_neo_tree_trait_color_signals_connected", False):
            return
        amr_button = getattr(self, "neo_tree_amr_color_button", None)
        mge_button = getattr(self, "neo_tree_mge_color_button", None)
        virulence_button = getattr(self, "neo_tree_virulence_color_button", None)
        if amr_button is not None:
            amr_button.clicked.connect(
                lambda: self._choose_neo_tree_trait_present_color(
                    NEO_TREE_RESISTANCE_TOKEN,
                    "Choose AMR present color",
                )
            )
        if mge_button is not None:
            mge_button.clicked.connect(
                lambda: self._choose_neo_tree_trait_present_color(
                    NEO_TREE_MGE_TOKEN,
                    "Choose MGE present color",
                )
            )
        if virulence_button is not None:
            virulence_button.clicked.connect(
                lambda: self._choose_neo_tree_trait_present_color(
                    NEO_TREE_VIRULENCE_TOKEN,
                    "Choose virulence present color",
                )
            )
        self._neo_tree_trait_color_signals_connected = True

    def _choose_neo_tree_trait_present_color(self, token: str, title: str) -> None:
        initial = self._neo_tree_trait_present_color(token)
        chosen = QtWidgets.QColorDialog.getColor(initial, self, title)
        if not chosen.isValid():
            return
        self._set_neo_tree_trait_present_color(token, chosen)

    def _apply_neo_tree_trait_color_selection(self) -> None:
        amr_color = self._neo_tree_trait_present_color(NEO_TREE_RESISTANCE_TOKEN)
        mge_color = self._neo_tree_trait_present_color(NEO_TREE_MGE_TOKEN)
        virulence_color = self._neo_tree_trait_present_color(NEO_TREE_VIRULENCE_TOKEN)
        for attr_name in ("pan_neo_tree_tab", "pan_neo_tab"):
            panel = getattr(self, attr_name, None)
            if panel is not None and hasattr(panel, "set_special_track_colors"):
                panel.set_special_track_colors(
                    amr_color=amr_color,
                    mge_color=mge_color,
                    virulence_color=virulence_color,
                )

    @staticmethod
    def _neo_tree_trait_key(value: object) -> str:
        return re.sub(r"[^a-z0-9]+", "", str(value or "").casefold())

    def _collect_neo_tree_trait_available_columns(self) -> List[str]:
        available: List[str] = []
        if not self.meta_df.empty:
            seen: set[str] = set()
            for raw_col in self.meta_columns:
                column = str(raw_col).strip()
                if not column or column in seen:
                    continue
                seen.add(column)
                key = self._neo_tree_trait_key(column)
                if not key or key in _NEO_TREE_EXCLUDED_TRAIT_KEYS:
                    continue
                if "latitude" in key or "longitude" in key:
                    continue
                series = self.meta_df[column].fillna("").astype(str).str.strip()
                cleaned = [
                    value
                    for value in series.tolist()
                    if value
                    and value.casefold() not in {"nan", "none", "null", "na", "n/a", "<na>", "missing", "unknown"}
                ]
                if not cleaned:
                    continue
                unique_values = sorted({value for value in cleaned})
                unique_count = len(unique_values)
                if unique_count < 2:
                    continue
                non_missing_count = len(cleaned)
                numeric_like = pd.to_numeric(pd.Series(cleaned, dtype="object"), errors="coerce").notna().sum()
                if numeric_like >= max(3, int(non_missing_count * 0.9)) and unique_count > 12:
                    continue
                if unique_count >= max(25, int(non_missing_count * 0.85)):
                    continue
                available.append(column)
        if getattr(self, "_neo_tree_has_resistance_traits", False):
            available.append(NEO_TREE_RESISTANCE_TOKEN)
        if getattr(self, "_neo_tree_has_mge_traits", False):
            available.append(NEO_TREE_MGE_TOKEN)
        if getattr(self, "_neo_tree_has_virulence_traits", False):
            available.append(NEO_TREE_VIRULENCE_TOKEN)
        plasmid_index = getattr(self, "_plasmid_cluster_index", None)
        if plasmid_index:
            available.append(NEO_TREE_PLASMID_TOKEN)
        if self.output_dir is not None and _find_panaroo_out(self.output_dir) is not None:
            available.append(NEO_TREE_PANAROO_TOKEN)
        return available

    @staticmethod
    def _neo_tree_trait_label(value: str) -> str:
        if str(value) == NEO_TREE_RESISTANCE_TOKEN:
            return NEO_TREE_RESISTANCE_LABEL
        if str(value) == NEO_TREE_MGE_TOKEN:
            return NEO_TREE_MGE_LABEL
        if str(value) == NEO_TREE_VIRULENCE_TOKEN:
            return NEO_TREE_VIRULENCE_LABEL
        if str(value) == NEO_TREE_PANAROO_TOKEN:
            return NEO_TREE_PANAROO_LABEL
        if str(value) == NEO_TREE_PLASMID_TOKEN:
            return NEO_TREE_PLASMID_LABEL
        return str(value)

    def _detect_neo_tree_resistance_traits(self) -> bool:
        df = getattr(self, "marker_hits_df", pd.DataFrame())
        if df is None or df.empty:
            return False
        if "DB" in df.columns:
            resistance_dbs = [
                str(db).strip()
                for db in sorted({str(db).strip() for db in df["DB"].dropna().astype(str).tolist() if str(db).strip()})
                if hasattr(self, "_tree_db_is_resistance") and self._tree_db_is_resistance(str(db).strip())
            ]
            if resistance_dbs:
                return True
        for col in ("ResistanceDrugClasses", "ResistanceAlleles", "ResistanceGeneFamilies"):
            if col in df.columns and df[col].fillna("").astype(str).str.strip().ne("").any():
                return True
        return False

    def _detect_neo_tree_mge_traits(self) -> bool:
        df = getattr(self, "marker_hits_df", pd.DataFrame())
        if df is None or df.empty:
            return False
        if "DB" in df.columns:
            mge_dbs = [
                str(db).strip()
                for db in sorted({str(db).strip() for db in df["DB"].dropna().astype(str).tolist() if str(db).strip()})
                if hasattr(self, "_tree_db_is_mge_source") and self._tree_db_is_mge_source(str(db).strip())
            ]
            if not mge_dbs:
                return False
            db_values = df["DB"].fillna("").astype(str).str.strip()
            subset = df[db_values.isin(mge_dbs)]
        else:
            subset = df
        if subset.empty:
            return False
        for row in subset.to_dict(orient="records"):
            if _extract_mge_labels(row):
                return True
        return False

    def _detect_neo_tree_virulence_traits(self) -> bool:
        df = getattr(self, "marker_hits_df", pd.DataFrame())
        if df is None or df.empty:
            return False
        if "DB" in df.columns:
            virulence_dbs = [
                str(db).strip()
                for db in sorted({str(db).strip() for db in df["DB"].dropna().astype(str).tolist() if str(db).strip()})
                if hasattr(self, "_tree_db_is_virulence_source") and self._tree_db_is_virulence_source(str(db).strip())
            ]
            if not virulence_dbs:
                return False
            db_values = df["DB"].fillna("").astype(str).str.strip()
            subset = df[db_values.isin(virulence_dbs)]
        else:
            subset = df
        if subset.empty:
            return False
        for row in subset.to_dict(orient="records"):
            if _extract_virulence_labels(row):
                return True
        return False

    def _neo_tree_default_selected_traits(self) -> List[str]:
        available = list(getattr(self, "_neo_tree_trait_available", []) or [])
        if not available:
            return []
        by_key = {self._neo_tree_trait_key(column): column for column in available}
        for candidate_group in _NEO_TREE_PREFERRED_TRAIT_CANDIDATES:
            for candidate in candidate_group:
                match = by_key.get(self._neo_tree_trait_key(candidate))
                if match:
                    return [match]
        return [available[0]]

    def _selected_neo_tree_traits(self) -> List[str]:
        selected = getattr(self, "_neo_tree_trait_selected", None)
        available = list(getattr(self, "_neo_tree_trait_available", []) or [])
        if not isinstance(selected, list):
            return list(available[:1]) if available else []
        selected_set = {str(value).strip() for value in selected if str(value).strip()}
        return [column for column in available if column in selected_set]

    def _init_neo_tree_trait_controls(self) -> None:
        group = getattr(self, "neo_tree_traits_group", None)
        layout = getattr(self, "neo_tree_traits_layout", None)
        if group is None or layout is None:
            return
        if isinstance(group, QtWidgets.QGroupBox):
            _apply_analysis_group_style(group)

        available = list(getattr(self, "_neo_tree_trait_available", []) or [])
        self._neo_tree_trait_ui_updating = True
        self._neo_tree_trait_checkboxes = {}
        self._neo_tree_trait_select_all_checkbox = None
        if not getattr(self, "_neo_tree_trait_selected", None):
            self._neo_tree_trait_selected = self._neo_tree_default_selected_traits()
        else:
            selected_set = {str(value).strip() for value in self._neo_tree_trait_selected}
            self._neo_tree_trait_selected = [column for column in available if column in selected_set]

        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()

        if available:
            select_all = QtWidgets.QCheckBox("Select all")
            select_all.setTristate(True)
            select_all.stateChanged.connect(self._on_neo_tree_trait_select_all_state_changed)
            layout.addWidget(select_all)
            self._neo_tree_trait_select_all_checkbox = select_all

            for column in available:
                cb = QtWidgets.QCheckBox(self._neo_tree_trait_label(column))
                cb.setChecked(column in self._neo_tree_trait_selected)
                cb.setStyleSheet("margin-left: 14px;")
                cb.stateChanged.connect(lambda _state, key=column: self._on_neo_tree_trait_item_state_changed(key))
                layout.addWidget(cb)
                self._neo_tree_trait_checkboxes[column] = cb

        layout.addStretch(1)
        self._connect_neo_tree_trait_color_controls()
        self._neo_tree_trait_ui_updating = False
        self._update_neo_tree_trait_group_state()
        self._update_neo_tree_trait_summary_label()
        self._update_neo_tree_trait_color_controls()

    def _update_neo_tree_trait_summary_label(self) -> None:
        label = getattr(self, "neo_tree_traits_summary_label", None)
        if label is None:
            return
        total = len(getattr(self, "_neo_tree_trait_available", []) or [])
        selected = len(self._selected_neo_tree_traits())
        label.setText(f"Visible neoTree traits: {selected}/{total}. Clades and clusters stay pinned.")

    def _update_neo_tree_trait_group_state(self) -> None:
        checkbox = getattr(self, "_neo_tree_trait_select_all_checkbox", None)
        if checkbox is None:
            return
        self._neo_tree_trait_ui_updating = True
        total = len(getattr(self, "_neo_tree_trait_available", []) or [])
        selected = len(self._selected_neo_tree_traits())
        if selected <= 0:
            checkbox.setCheckState(QtCore.Qt.CheckState.Unchecked)
        elif selected >= total:
            checkbox.setCheckState(QtCore.Qt.CheckState.Checked)
        else:
            checkbox.setCheckState(QtCore.Qt.CheckState.PartiallyChecked)
        self._neo_tree_trait_ui_updating = False

    def _sync_selected_neo_tree_traits_from_ui(self) -> None:
        available = list(getattr(self, "_neo_tree_trait_available", []) or [])
        self._neo_tree_trait_selected = [
            column
            for column in available
            if column in self._neo_tree_trait_checkboxes and self._neo_tree_trait_checkboxes[column].isChecked()
        ]

    def _apply_neo_tree_trait_selection(self) -> None:
        if getattr(self, "pan_neo_tree_tab", None) is not None:
            self.pan_neo_tree_tab.set_side_metadata_columns(self._selected_neo_tree_traits())

    def _on_neo_tree_trait_select_all_state_changed(self, state: int) -> None:
        if self._neo_tree_trait_ui_updating:
            return
        state_value = state.value if hasattr(state, "value") else int(state)
        if state_value == QtCore.Qt.CheckState.PartiallyChecked.value:
            return
        target_checked = state_value == QtCore.Qt.CheckState.Checked.value
        self._neo_tree_trait_ui_updating = True
        for checkbox in self._neo_tree_trait_checkboxes.values():
            checkbox.setChecked(target_checked)
        self._neo_tree_trait_ui_updating = False
        self._sync_selected_neo_tree_traits_from_ui()
        self._update_neo_tree_trait_group_state()
        self._update_neo_tree_trait_summary_label()
        self._apply_neo_tree_trait_selection()

    def _on_neo_tree_trait_item_state_changed(self, _column: str) -> None:
        if self._neo_tree_trait_ui_updating:
            return
        self._sync_selected_neo_tree_traits_from_ui()
        self._update_neo_tree_trait_group_state()
        self._update_neo_tree_trait_summary_label()
        self._apply_neo_tree_trait_selection()

    def _normalize_sample_subset_values(self, samples: Sequence[str]) -> List[str]:
        seen: set[str] = set()
        ordered: List[str] = []
        for raw in samples:
            sample = str(raw).strip()
            if not sample or sample in seen:
                continue
            seen.add(sample)
            ordered.append(sample)
        return ordered

    def _set_sample_subset_filter(
        self,
        *,
        label: str,
        samples: Sequence[str],
        kind: Optional[str] = None,
    ) -> None:
        normalized = self._normalize_sample_subset_values(samples)
        if not normalized:
            return
        self._sample_subset_filter_samples = normalized
        self._sample_subset_filter_label = str(label).strip() or None
        self._sample_subset_filter_kind = str(kind).strip() or None
        self._render_metadata_filters()
        self._schedule_refresh()

    def _clear_sample_subset_filter(self) -> None:
        if not getattr(self, "_sample_subset_filter_samples", None):
            return
        self._sample_subset_filter_samples = []
        self._sample_subset_filter_label = None
        self._sample_subset_filter_kind = None
        self._render_metadata_filters()
        self._schedule_refresh()

    def _show_add_filter_menu(self) -> None:
        menu = QtWidgets.QMenu(self)
        for col in sorted(self.meta_columns):
            label = f"{col} (filtered)" if col in self.metadata_filters else col
            action = menu.addAction(label)
            action.triggered.connect(lambda _checked=False, c=col: self._edit_metadata_filter(c))
        anchor = self.add_filter_btn.mapToGlobal(QtCore.QPoint(0, self.add_filter_btn.height()))
        menu.exec(anchor)


    def _init_table_column_controls(self) -> None:
        if not hasattr(self, "table_columns_layout") or self.table_columns_layout is None:
            return
        self._table_column_ui_updating = True
        self._table_column_group_checkboxes = {}
        self._table_column_checkboxes = {}
        self._table_column_group_members = {}
        if not self._table_column_selected:
            self._table_column_selected = self._table_default_selected_columns()
        else:
            selected_set = {str(col) for col in self._table_column_selected}
            self._table_column_selected = [
                key for key, _label in self._table_column_items() if key in selected_set
            ]

        while self.table_columns_layout.count():
            item = self.table_columns_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()

        for group_key, group_label, group_items in self._table_column_groups():
            if not group_items:
                continue
            group_box = QtWidgets.QGroupBox(group_label)
            group_box.setStyleSheet("QGroupBox { font-weight: 600; }")
            group_layout = QtWidgets.QVBoxLayout(group_box)
            group_layout.setContentsMargins(8, 10, 8, 6)
            group_layout.setSpacing(4)

            group_toggle = QtWidgets.QCheckBox("Select all")
            group_toggle.setTristate(True)
            group_toggle.stateChanged.connect(
                lambda state, g=group_key: self._on_table_column_group_state_changed(g, state)
            )
            group_layout.addWidget(group_toggle)
            self._table_column_group_checkboxes[group_key] = group_toggle

            members: List[str] = []
            for col_key, col_label in group_items:
                cb = QtWidgets.QCheckBox(col_label)
                cb.setChecked(col_key in self._table_column_selected)
                cb.stateChanged.connect(
                    lambda _state, key=col_key: self._on_table_column_item_state_changed(key)
                )
                cb.setStyleSheet("margin-left: 14px;")
                group_layout.addWidget(cb)
                self._table_column_checkboxes[col_key] = cb
                members.append(col_key)
            self._table_column_group_members[group_key] = members
            self.table_columns_layout.addWidget(group_box)

        self.table_columns_layout.addStretch(1)
        self._table_column_ui_updating = False
        self._update_table_column_group_states()
        self._update_table_column_summary_label()


    def _update_table_column_summary_label(self) -> None:
        label = getattr(self, "table_columns_summary_label", None)
        if label is None:
            return
        total = len(self._table_column_items())
        selected = len(self._selected_table_columns())
        label.setText(
            f"Visible sample columns: {selected}/{total}. "
            "Sample / Contig, Card, and Plasmid stay pinned."
        )


    def _selected_table_column_set(self) -> set[str]:
        return {str(col) for col in self._selected_table_columns()}


    def _update_table_column_group_states(self) -> None:
        self._table_column_ui_updating = True
        selected = self._selected_table_column_set()
        for group_key, checkbox in self._table_column_group_checkboxes.items():
            members = self._table_column_group_members.get(group_key, [])
            checked_count = sum(1 for key in members if key in selected)
            if checked_count == 0:
                checkbox.setCheckState(QtCore.Qt.CheckState.Unchecked)
            elif checked_count == len(members):
                checkbox.setCheckState(QtCore.Qt.CheckState.Checked)
            else:
                checkbox.setCheckState(QtCore.Qt.CheckState.PartiallyChecked)
        self._table_column_ui_updating = False


    def _sync_selected_table_columns_from_ui(self) -> None:
        available_order = [key for key, _label in self._table_column_items()]
        self._table_column_selected = [
            key
            for key in available_order
            if key in self._table_column_checkboxes
            and self._table_column_checkboxes[key].isChecked()
        ]


    def _on_table_column_group_state_changed(self, group_key: str, state: int) -> None:
        if self._table_column_ui_updating:
            return
        state_value = state.value if hasattr(state, "value") else int(state)
        if state_value == QtCore.Qt.CheckState.PartiallyChecked.value:
            return
        target_checked = state_value == QtCore.Qt.CheckState.Checked.value
        self._table_column_ui_updating = True
        for key in self._table_column_group_members.get(group_key, []):
            checkbox = self._table_column_checkboxes.get(key)
            if checkbox is not None:
                checkbox.setChecked(target_checked)
        self._table_column_ui_updating = False
        self._sync_selected_table_columns_from_ui()
        self._update_table_column_group_states()
        self._update_table_column_summary_label()
        self._rerender_visible_table_columns()


    def _on_table_column_item_state_changed(self, _column_key: str) -> None:
        if self._table_column_ui_updating:
            return
        self._sync_selected_table_columns_from_ui()
        self._update_table_column_group_states()
        self._update_table_column_summary_label()
        self._rerender_visible_table_columns()


    def _init_predictor_model(self) -> None:
        model = QtGui.QStandardItemModel(self.predictor_box)
        clear_item = QtGui.QStandardItem("Any (clear)")
        clear_item.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled)
        clear_item.setData("__clear__", QtCore.Qt.ItemDataRole.UserRole)
        model.appendRow(clear_item)
        for pred in self.predictor_options:
            item = QtGui.QStandardItem(pred)
            item.setFlags(
                QtCore.Qt.ItemFlag.ItemIsEnabled
                | QtCore.Qt.ItemFlag.ItemIsUserCheckable
            )
            item.setData(QtCore.Qt.CheckState.Unchecked, QtCore.Qt.ItemDataRole.CheckStateRole)
            item.setData(pred, QtCore.Qt.ItemDataRole.UserRole)
            model.appendRow(item)
        model.itemChanged.connect(self._on_predictor_model_changed)
        self._predictor_model = model
        self.predictor_box.setModel(model)
        self._update_predictor_box_text()


    def _selected_predictors(self) -> List[str]:
        preds: List[str] = []
        if not self._predictor_model:
            return preds
        for row in range(self._predictor_model.rowCount()):
            item = self._predictor_model.item(row)
            if not item:
                continue
            token = item.data(QtCore.Qt.ItemDataRole.UserRole)
            if not token or token == "__clear__":
                continue
            if item.checkState() == QtCore.Qt.CheckState.Checked:
                preds.append(str(token))
        return preds


    def _update_predictor_box_text(self) -> None:
        selected = self._selected_predictors()
        if selected:
            txt = f"Predictor ({','.join(sorted(selected))})"
        else:
            txt = "Predictor (any)"
        if self.predictor_box.isEditable() and self.predictor_box.lineEdit() is not None:
            self.predictor_box.lineEdit().setText(txt)
        else:
            self.predictor_box.setCurrentText(txt)


    def _toggle_predictor_item(self, index: QtCore.QModelIndex) -> None:
        if not self._predictor_model:
            return
        item = self._predictor_model.itemFromIndex(index)
        if not item:
            return
        token = item.data(QtCore.Qt.ItemDataRole.UserRole)
        if token == "__clear__":
            for row in range(self._predictor_model.rowCount()):
                it = self._predictor_model.item(row)
                if not it:
                    continue
                if it.data(QtCore.Qt.ItemDataRole.UserRole) == "__clear__":
                    continue
                it.setCheckState(QtCore.Qt.CheckState.Unchecked)
            self._update_predictor_box_text()
            self._schedule_refresh()
            return
        state = item.checkState()
        item.setCheckState(
            QtCore.Qt.CheckState.Unchecked
            if state == QtCore.Qt.CheckState.Checked
            else QtCore.Qt.CheckState.Checked
        )


    def _on_predictor_model_changed(self, _item: QtGui.QStandardItem) -> None:
        self._update_predictor_box_text()
        self._schedule_refresh()


    def _init_tag_model(self) -> None:
        model = QtGui.QStandardItemModel(self.tag_box)
        clear_item = QtGui.QStandardItem("Any (clear)")
        clear_item.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled)
        clear_item.setData("__clear__", QtCore.Qt.ItemDataRole.UserRole)
        model.appendRow(clear_item)
        for tag, label in TAG_OPTIONS:
            item = QtGui.QStandardItem(f"{tag} — {label}")
            item.setFlags(
                QtCore.Qt.ItemFlag.ItemIsEnabled
                | QtCore.Qt.ItemFlag.ItemIsUserCheckable
            )
            item.setData(QtCore.Qt.CheckState.Unchecked, QtCore.Qt.ItemDataRole.CheckStateRole)
            item.setData(tag, QtCore.Qt.ItemDataRole.UserRole)
            model.appendRow(item)
        model.itemChanged.connect(self._on_tag_model_changed)
        self._tag_model = model
        self.tag_box.setModel(model)
        self._update_tag_box_text()


    def _selected_tags(self) -> List[str]:
        tags: List[str] = []
        if not self._tag_model:
            return tags
        for row in range(self._tag_model.rowCount()):
            item = self._tag_model.item(row)
            if not item:
                continue
            tag = item.data(QtCore.Qt.ItemDataRole.UserRole)
            if not tag or tag == "__clear__":
                continue
            if item.checkState() == QtCore.Qt.CheckState.Checked:
                tags.append(str(tag))
        return tags


    def _update_tag_box_text(self) -> None:
        selected = self._selected_tags()
        if selected:
            txt = f"Tags ({','.join(sorted(selected))})"
        else:
            txt = "Tags (any)"
        if self.tag_box.isEditable() and self.tag_box.lineEdit() is not None:
            self.tag_box.lineEdit().setText(txt)
        else:
            self.tag_box.setCurrentText(txt)


    def _toggle_tag_item(self, index: QtCore.QModelIndex) -> None:
        if not self._tag_model:
            return
        item = self._tag_model.itemFromIndex(index)
        if not item:
            return
        tag = item.data(QtCore.Qt.ItemDataRole.UserRole)
        if tag == "__clear__":
            for row in range(self._tag_model.rowCount()):
                it = self._tag_model.item(row)
                if not it:
                    continue
                if it.data(QtCore.Qt.ItemDataRole.UserRole) == "__clear__":
                    continue
                it.setCheckState(QtCore.Qt.CheckState.Unchecked)
            self._update_tag_box_text()
            self._schedule_refresh()
            return
        state = item.checkState()
        item.setCheckState(
            QtCore.Qt.CheckState.Unchecked
            if state == QtCore.Qt.CheckState.Checked
            else QtCore.Qt.CheckState.Checked
        )


    def _on_tag_model_changed(self, _item: QtGui.QStandardItem) -> None:
        self._update_tag_box_text()
        self._schedule_refresh()


    def _edit_metadata_filter(self, column: str) -> None:
        values = self._metadata_value_options(column)
        current = set(self.metadata_filters.get(column, []))
        dlg = MetadataValueDialog(self, column, values, current)
        if dlg.exec() != QtWidgets.QDialog.DialogCode.Accepted:
            return
        selected = dlg.selected_values()
        if not selected:
            if column in self.metadata_filters:
                del self.metadata_filters[column]
        else:
            self.metadata_filters[column] = selected
        self._render_metadata_filters()
        self._schedule_refresh()


    def _metadata_value_options(self, column: str) -> List[Tuple[str, str]]:
        if column not in self._meta_pl.columns:
            return []
        values_lf = self._meta_pl.lazy().select(pl.col(column).cast(pl.Utf8).str.strip_chars().alias("v"))
        has_missing = (
            values_lf
            .filter(pl.col("v").is_null() | (pl.col("v") == ""))
            .select(pl.len().alias("n"))
            .collect()
            .item(0, "n")
            > 0
        )
        values = (
            values_lf
            .filter(pl.col("v").is_not_null() & (pl.col("v") != ""))
            .unique()
            .sort("v")
            .collect()
            .get_column("v")
            .to_list()
        )
        values = sorted({str(v) for v in values if str(v).strip()})
        options: List[Tuple[str, str]] = [(v, v) for v in values]
        if has_missing:
            options.insert(0, (MISSING_FILTER_LABEL, MISSING_FILTER_TOKEN))
        return options


    def _render_metadata_filters(self) -> None:
        while self.filters_layout.count():
            item = self.filters_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
        subset_samples = list(getattr(self, "_sample_subset_filter_samples", []) or [])
        subset_label = str(getattr(self, "_sample_subset_filter_label", "") or "").strip()
        if subset_samples:
            row = QtWidgets.QWidget()
            row_layout = QtWidgets.QHBoxLayout(row)
            row_layout.setContentsMargins(0, 0, 0, 0)
            label = QtWidgets.QLabel("Sample set")
            summary_text = subset_label or f"{len(subset_samples)} samples"
            summary = QtWidgets.QLabel(f"{summary_text} | n={len(subset_samples)}")
            summary.setToolTip(", ".join(subset_samples))
            remove_btn = QtWidgets.QPushButton("Remove")
            remove_btn.clicked.connect(self._clear_sample_subset_filter)
            row_layout.addWidget(label)
            row_layout.addWidget(summary, stretch=1)
            row_layout.addWidget(remove_btn)
            self.filters_layout.addWidget(row)
        for column, values in self.metadata_filters.items():
            row = QtWidgets.QWidget()
            row_layout = QtWidgets.QHBoxLayout(row)
            row_layout.setContentsMargins(0, 0, 0, 0)
            label = QtWidgets.QLabel(column)
            summary = QtWidgets.QLabel(self._format_filter_values(values))
            summary.setToolTip(self._format_filter_values(values, full=True))
            remove_btn = QtWidgets.QPushButton("Remove")
            remove_btn.clicked.connect(lambda _checked=False, c=column: self._remove_metadata_filter(c))
            row_layout.addWidget(label)
            row_layout.addWidget(summary, stretch=1)
            row_layout.addWidget(remove_btn)
            self.filters_layout.addWidget(row)
        self.filters_layout.addStretch(1)


    def _remove_metadata_filter(self, column: str) -> None:
        if column in self.metadata_filters:
            del self.metadata_filters[column]
            self._render_metadata_filters()
            self._schedule_refresh()


    def _format_filter_values(self, values: Sequence[str], full: bool = False) -> str:
        display: List[str] = []
        for v in values:
            if v == MISSING_FILTER_TOKEN:
                display.append(MISSING_FILTER_LABEL)
            else:
                display.append(str(v))
        if full or len(display) <= 5:
            return ", ".join(display)
        return ", ".join(display[:5]) + f", … (+{len(display) - 5} more)"


    def _schedule_refresh(self) -> None:
        if self._refresh_timer.isActive():
            self._refresh_timer.stop()
        self._refresh_timer.start()


    def _on_view_mode_changed(self, mode: str) -> None:
        self.view_mode = mode
        if mode == "Contig":
            self._igv_dirty = True
        self.refresh()


    def _current_metadata_filters(self) -> Tuple[Optional[str], Dict[str, List[str]]]:
        sample_q = self.sample_query.text().strip() or None
        return sample_q, {k: list(v) for k, v in self.metadata_filters.items()}


    def _current_annotation_filters(self) -> Dict[str, Any]:
        predictors = self._selected_predictors()
        tags = self._selected_tags()
        molecule = ""
        molecule_box = getattr(self, "molecule_filter_box", None)
        if isinstance(molecule_box, QtWidgets.QComboBox):
            data = molecule_box.currentData()
            molecule = str(data if data is not None else "").strip().lower()
        return {
            "gene_text": self.gene_query.text().strip() or None,
            "product_text": self.product_query.text().strip() or None,
            "molecule": molecule if molecule in {"plasmid", "chromosome"} else "",
            "predictors": predictors,
            "tags": tags,
            "hide_hypothetical": self.hide_hypo_chk.isChecked(),
        }


    def _has_active_annotation_filter(self) -> bool:
        f = self._current_annotation_filters()
        return any(
            [
                bool(f.get("gene_text")),
                bool(f.get("product_text")),
                bool(f.get("molecule")),
                bool(f.get("predictors")),
                bool(f.get("tags")),
                bool(f.get("hide_hypothetical")),
            ]
        )


    def _should_show_global(self) -> bool:
        search_text = self.sample_query.text().strip()
        if search_text:
            return True
        if self._has_active_annotation_filter():
            return True
        return False
