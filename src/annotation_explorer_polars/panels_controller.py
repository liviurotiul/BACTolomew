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
from .cgmlst_force_plot import CgmlstForcePlotPanel
from .geo_experimental_views import GeoMapLibreView
from .geo_view import GeoMapView
from .igv_collection_view import IgvCollectionDialog
from .igv_view import GeneSegment, IgvTrackView
from .metadata_filter_dialog import MetadataValueDialog
from .pan_view import PanHeatmapWidget
from .panaroo_data import (
    _find_fastgear_sample_clades,
    _load_fastgear_clusters,
    _load_fastgear_lineages,
)
from .panaroo_style import _apply_analysis_panel_style
from .panaroo_views import (
    EntropyTrackTab,
    PanCoreAccessoryTab,
    StructPresenceHeatmapTab,
)


class PanelsControllerMixin:
    _PANAROO_TAB_NAMES: Tuple[str, ...] = ("neoPan",)
    _PLOT_WIDGET_NAMES: Tuple[str, ...] = ("neoTree", "cgMLSTForce")
    _GEO_TAB_NAMES: Tuple[str, ...] = ("Geo",)
    _COUNTRY_TOOLTIP_CANDIDATES: Tuple[str, ...] = ("Country", "_country")
    _LATITUDE_TOOLTIP_CANDIDATES: Tuple[str, ...] = (
        "Latitude",
        "latitude",
        "Lat",
        "lat",
    )
    _LONGITUDE_TOOLTIP_CANDIDATES: Tuple[str, ...] = (
        "Longitude",
        "longitude",
        "Lon",
        "lon",
        "Long",
        "long",
    )
    _ST_TOOLTIP_CANDIDATES: Tuple[str, ...] = (
        "ST",
        "Sequence Type",
        "SequenceType",
        "sequence_type",
        "MLST ST",
        "MLST_ST",
    )
    _COLLECTION_DATE_TOOLTIP_CANDIDATES: Tuple[str, ...] = (
        "Collection date",
        "Collection Date",
        "Collection_date",
        "collection_date",
        "CollectionDate",
        "Date",
        "date",
    )

    def _bottom_tab_by_name(self, name: str) -> Optional[QtWidgets.QWidget]:
        for i in range(self.bottom_tabs.count()):
            if self.bottom_tabs.tabText(i) == name:
                return self.bottom_tabs.widget(i)
        return None

    def _prepare_tab_layout(self, host: Optional[QtWidgets.QWidget]) -> Optional[QtWidgets.QVBoxLayout]:
        if host is None:
            return None
        layout = host.layout()
        if layout is None:
            layout = QtWidgets.QVBoxLayout(host)
            layout.setContentsMargins(0, 0, 0, 0)
        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.setParent(None)
                widget.deleteLater()
        return layout

    def _set_bottom_tab_placeholder(self, name: str, text: str) -> None:
        tab = self._bottom_tab_by_name(name)
        layout = self._prepare_tab_layout(tab)
        if layout is None:
            return
        label = QtWidgets.QLabel(text)
        label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        label.setWordWrap(True)
        label.setStyleSheet("color: #64748b; font-size: 13px; padding: 18px;")
        layout.addWidget(label)


    def _on_bottom_tab_changed(self, _index: int) -> None:
        # Prevent stale hover highlights when switching between bottom tabs.
        self._clear_pan_highlight()
        self._ensure_active_bottom_tab_initialized()
        self._flush_deferred_panel_refresh()


    def _active_bottom_tab_name(self) -> Optional[str]:
        if not hasattr(self, "bottom_tabs"):
            return None
        index = self.bottom_tabs.currentIndex()
        if index < 0:
            return None
        return self.bottom_tabs.tabText(index)

    def _ensure_active_bottom_tab_initialized(self) -> None:
        self._ensure_bottom_tab_initialized(self._active_bottom_tab_name())

    def _ensure_bottom_tab_initialized(self, name: Optional[str]) -> None:
        if name == "Geo":
            if getattr(self, "map_web_view", None) is None:
                self._init_geo_map()
                self._update_map()
            return
        if name == "neoPan":
            if getattr(self, "pan_neo_tab", None) is None:
                self._init_panaroo_analytics_tabs()
                self._refresh_panaroo_tab(
                    "neoPan",
                    self._current_panel_refresh_samples(),
                    getattr(self, "ann_filters", {}),
                )
            return

    @staticmethod
    def _ordered_sample_names(samples: Sequence[str]) -> List[str]:
        seen: set[str] = set()
        ordered: List[str] = []
        for raw in samples:
            sample = str(raw).strip()
            if not sample or sample in seen:
                continue
            seen.add(sample)
            ordered.append(sample)
        return ordered


    def _all_geo_views(self) -> List[QtWidgets.QWidget]:
        views: List[QtWidgets.QWidget] = []
        for attr_name in ("map_web_view",):
            view = getattr(self, attr_name, None)
            if view is not None:
                views.append(view)
        if self.pan_neo_tab is not None:
            try:
                views.extend(self.pan_neo_tab.all_geo_views())
            except Exception:
                if self.pan_neo_tab.geo_view is not None:
                    views.append(self.pan_neo_tab.geo_view)
        return views

    def _init_geo_map_view(
        self,
        host: Optional[QtWidgets.QWidget],
        *,
        view_attr: str,
        title: str,
    ) -> None:
        if host is None:
            return
        if getattr(self, view_attr, None) is not None:
            return
        layout = self._prepare_tab_layout(host)
        if layout is None:
            return
        try:
            view = GeoMapLibreView(parent=host)
        except Exception as exc:
            msg = QtWidgets.QLabel(f"Failed to initialize {title}: {exc}")
            msg.setWordWrap(True)
            msg.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
            layout.addWidget(msg)
            return
        setattr(self, view_attr, view)
        layout.addWidget(view)

    def _init_geo_map(self) -> None:
        self._init_geo_map_view(
            self._bottom_tab_by_name("Geo"),
            view_attr="map_web_view",
            title="Geo map",
        )
        self._apply_pan_visual_highlight()

    @staticmethod
    def _normalize_plot_trait_value(value: object) -> str:
        try:
            if pd.isna(value):
                return ""
        except Exception:
            pass
        text = str(value).strip()
        if not text:
            return ""
        if text.casefold() in {"nan", "none", "null", "na", "n/a", "<na>", "missing", "unknown"}:
            return ""
        return text

    @staticmethod
    def _resolve_plot_metadata_column(columns: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
        normalized = {str(column).strip().casefold(): str(column) for column in columns if str(column).strip()}
        for candidate in candidates:
            key = str(candidate).strip().casefold()
            if key in normalized:
                return normalized[key]
        return None

    def _load_plot_mlst_sequence_types(self) -> Dict[str, str]:
        if self.output_dir is None:
            return {}
        path = self.output_dir / "mlst" / "results.csv"
        if not path.is_file():
            return {}
        try:
            df = pd.read_csv(path, dtype=str, low_memory=False)
        except Exception:
            return {}
        columns = {str(column).strip(): column for column in df.columns}
        sample_col = columns.get("sample")
        st_col = columns.get("sequence_type")
        if sample_col is None or st_col is None:
            return {}
        sequence_types: Dict[str, str] = {}
        for sample_value, st_value in df[[sample_col, st_col]].itertuples(index=False):
            sample = self._normalize_plot_trait_value(sample_value)
            sequence_type = self._normalize_plot_trait_value(st_value)
            if sample and sequence_type and sample not in sequence_types:
                sequence_types[sample] = sequence_type
        return sequence_types

    def _build_plot_sample_hover_metadata(self) -> Dict[str, Dict[str, str]]:
        metadata: Dict[str, Dict[str, str]] = {}
        if not self.meta_df.empty and "Sample" in self.meta_df.columns:
            columns = [str(column) for column in self.meta_df.columns]
            country_col = self._resolve_plot_metadata_column(columns, self._COUNTRY_TOOLTIP_CANDIDATES)
            st_col = self._resolve_plot_metadata_column(columns, self._ST_TOOLTIP_CANDIDATES)
            date_col = self._resolve_plot_metadata_column(columns, self._COLLECTION_DATE_TOOLTIP_CANDIDATES)
            for _idx, row in self.meta_df.iterrows():
                sample = self._normalize_plot_trait_value(row.get("Sample"))
                if not sample:
                    continue
                fields = metadata.setdefault(sample, {})
                if country_col and "Country" not in fields:
                    value = self._normalize_plot_trait_value(row.get(country_col))
                    if value:
                        fields["Country"] = value
                if st_col and "ST" not in fields:
                    value = self._normalize_plot_trait_value(row.get(st_col))
                    if value:
                        fields["ST"] = value
                if date_col and "Collection date" not in fields:
                    value = self._normalize_plot_trait_value(row.get(date_col))
                    if value:
                        fields["Collection date"] = value

        for sample, sequence_type in self._load_plot_mlst_sequence_types().items():
            if not sequence_type:
                continue
            metadata.setdefault(sample, {}).setdefault("ST", sequence_type)
        return {sample: fields for sample, fields in metadata.items() if fields}

    def _load_cgmlst_fastgear_maps(self) -> Tuple[Dict[str, str], Dict[str, str]]:
        if self.output_dir is None:
            return {}, {}
        sample_clades_path = _find_fastgear_sample_clades(self.output_dir)
        if sample_clades_path is None:
            return {}, {}
        try:
            lineage_by_sample = _load_fastgear_lineages(sample_clades_path)
        except Exception:
            lineage_by_sample = {}
        try:
            cluster_by_sample = _load_fastgear_clusters(sample_clades_path)
        except Exception:
            cluster_by_sample = {}
        return lineage_by_sample, cluster_by_sample

    def _build_cgmlst_trait_values(self) -> Dict[str, Dict[str, str]]:
        traits: Dict[str, Dict[str, str]] = {}
        if not self.meta_df.empty and "Sample" in self.meta_df.columns:
            for column in self.meta_df.columns:
                if str(column) == "Sample":
                    continue
                values: Dict[str, str] = {}
                for sample_value, raw_value in self.meta_df[["Sample", column]].itertuples(index=False):
                    sample = self._normalize_plot_trait_value(sample_value)
                    value = self._normalize_plot_trait_value(raw_value)
                    if sample and value:
                        values.setdefault(sample, value)
                if values:
                    traits[str(column)] = values

        st_map = self._load_plot_mlst_sequence_types()
        if st_map and "ST" not in traits:
            traits["ST"] = st_map

        lineage_by_sample, cluster_by_sample = self._load_cgmlst_fastgear_maps()
        if lineage_by_sample:
            traits["Clade"] = lineage_by_sample
        if cluster_by_sample:
            traits["Cluster"] = cluster_by_sample
        return traits

    def _init_plot_menu(self) -> None:
        menu_bar = self.menuBar()
        plot_menu = menu_bar.addMenu("Plot")
        phandango_action = plot_menu.addAction("Phandango style plot")
        phandango_action.triggered.connect(self._show_phandango_style_plot)
        phandango_action.setEnabled(self.output_dir is not None)
        cgmlst_force_action = plot_menu.addAction("cgMLST force graph")
        cgmlst_force_action.triggered.connect(self._show_cgmlst_force_plot)
        cgmlst_force_action.setEnabled(self.output_dir is not None)
        self._plot_menu = plot_menu
        self._phandango_plot_action = phandango_action
        self._cgmlst_force_plot_action = cgmlst_force_action

    def _current_panel_refresh_samples(self) -> List[str]:
        if getattr(self, "_deferred_panel_allowed_samples", None):
            return self._ordered_sample_names(self._deferred_panel_allowed_samples)
        if getattr(self, "_base_allowed_samples", None):
            return self._ordered_sample_names(self._base_allowed_samples)
        return []

    def _ensure_pan_neo_tree_tab(self, parent: Optional[QtWidgets.QWidget] = None) -> Optional[PanCoreAccessoryTab]:
        if self.output_dir is None:
            return None
        if self.pan_neo_tree_tab is not None:
            return self.pan_neo_tree_tab

        tree_tab = PanCoreAccessoryTab(
            self.output_dir,
            self.meta_df,
            reference_sample=self._snippy_reference_sample,
            fallback_tree=self._resolve_raxml_tree_path(),
            preserve_branch_lengths=True,
            color_tree_by_fastgear_clade=True,
            color_tree_by_fastgear_cluster=True,
            optimize_tree_by_fastgear_clade=True,
            annotate_fastgear_clusters=False,
            use_fastgear_clade_backgrounds=True,
            show_clade_legend=True,
            show_cluster_legend=True,
            show_heatmap=False,
            show_info_label=False,
            enable_tree_hover_interactions=False,
            side_metadata_title="K locus",
            side_metadata_candidates=("K locus", "K_locus", "KLocus", "KL", "KL type", "K type"),
            side_metadata_columns=self._selected_neo_tree_traits(),
            marker_hits_df=self.marker_hits_df,
            resistance_db_names=[db for db in self.annotation_db_columns if self._tree_db_is_resistance(db)],
            mge_db_names=[db for db in self.annotation_db_columns if self._tree_db_is_mge_source(db)],
            virulence_db_names=[db for db in self.annotation_db_columns if self._tree_db_is_virulence_source(db)],
            amr_present_color=self._neo_tree_amr_color,
            mge_present_color=self._neo_tree_mge_color,
            virulence_present_color=self._neo_tree_virulence_color,
            recompute_tree_from_sample_filter=True,
            parent=parent,
        )
        tree_tab.sampleHovered.connect(self._on_pan_sample_hovered)
        tree_tab.cladeFilterRequested.connect(self._on_pan_clade_filter_requested)
        tree_tab.clusterFilterRequested.connect(self._on_pan_cluster_filter_requested)
        tree_tab.sampleOpenRequested.connect(self._on_pan_sample_open_requested)
        self.pan_neo_tree_tab = tree_tab
        return tree_tab

    def _ensure_phandango_plot_dialog(self) -> Optional[QtWidgets.QDialog]:
        dlg = self._phandango_plot_dialog
        if dlg is not None:
            return dlg
        dlg = QtWidgets.QDialog(self)
        dlg.setWindowTitle("Phandango style plot")
        dlg.setModal(False)
        dlg.resize(1720, 940)
        tree_tab = self._ensure_pan_neo_tree_tab(parent=dlg)
        if tree_tab is None:
            dlg.deleteLater()
            return None
        layout = QtWidgets.QVBoxLayout(dlg)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(8)

        main_split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal, dlg)
        main_split.setChildrenCollapsible(False)
        main_split.setHandleWidth(8)
        layout.addWidget(main_split, stretch=1)

        controls_panel = QtWidgets.QWidget(main_split)
        _apply_analysis_panel_style(controls_panel)
        controls_layout = QtWidgets.QVBoxLayout(controls_panel)
        controls_layout.setContentsMargins(10, 10, 10, 10)
        controls_layout.setSpacing(8)

        if getattr(self, "neo_tree_traits_group", None) is not None:
            self.neo_tree_traits_group.setMinimumWidth(250)
            self.neo_tree_traits_group.setMaximumWidth(340)
            self.neo_tree_traits_group.setVisible(bool(getattr(self, "_neo_tree_trait_available", []) or []))
            controls_layout.addWidget(self.neo_tree_traits_group, stretch=1)
        else:
            controls_layout.addStretch(1)

        plot_split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal, main_split)
        plot_split.setChildrenCollapsible(False)
        plot_split.setHandleWidth(8)
        plot_split.addWidget(tree_tab)

        self._phandango_plot_splitter = plot_split
        self._phandango_pan_heatmap = None

        main_split.addWidget(controls_panel)
        main_split.addWidget(plot_split)
        main_split.setStretchFactor(0, 0)
        main_split.setStretchFactor(1, 1)
        main_split.setSizes([300, 1420])
        self._phandango_plot_dialog = dlg
        return dlg

    def _show_phandango_style_plot(self) -> None:
        dlg = self._ensure_phandango_plot_dialog()
        if dlg is None or self.pan_neo_tree_tab is None:
            return
        if "neoTree" in self._deferred_pan_tab_dirty:
            self._refresh_panaroo_tab(
                "neoTree",
                self._current_panel_refresh_samples(),
                self.ann_filters,
            )
            self._deferred_pan_tab_dirty.discard("neoTree")
        dlg.show()
        dlg.raise_()
        dlg.activateWindow()
        self._apply_pan_visual_highlight()

    def _ensure_cgmlst_force_plot_dialog(self) -> Optional[QtWidgets.QDialog]:
        if self.output_dir is None:
            return None
        dlg = self._cgmlst_force_plot_dialog
        if dlg is not None:
            return dlg

        dlg = QtWidgets.QDialog(self)
        dlg.setWindowTitle("cgMLST force graph")
        dlg.setModal(False)
        dlg.resize(1540, 920)

        layout = QtWidgets.QVBoxLayout(dlg)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(8)

        panel = CgmlstForcePlotPanel(
            self.output_dir,
            reference_sample=self._snippy_reference_sample,
            parent=dlg,
        )
        panel.sampleHovered.connect(self._on_pan_sample_hovered)
        panel.set_trait_values(self._build_cgmlst_trait_values())
        layout.addWidget(panel, stretch=1)

        self._cgmlst_force_plot_panel = panel
        self._cgmlst_force_plot_dialog = dlg
        return dlg

    def _show_cgmlst_force_plot(self) -> None:
        dlg = self._ensure_cgmlst_force_plot_dialog()
        panel = self._cgmlst_force_plot_panel
        if dlg is None or panel is None:
            return
        if "cgMLSTForce" in self._deferred_pan_tab_dirty:
            self._refresh_panaroo_tab(
                "cgMLSTForce",
                self._current_panel_refresh_samples(),
                self.ann_filters,
            )
            self._deferred_pan_tab_dirty.discard("cgMLSTForce")
        elif not panel.isVisible():
            self._refresh_panaroo_tab(
                "cgMLSTForce",
                self._current_panel_refresh_samples(),
                self.ann_filters,
            )
        dlg.show()
        dlg.raise_()
        dlg.activateWindow()
        self._apply_pan_visual_highlight()


    def _init_panaroo_analytics_tabs(self) -> None:
        if self.output_dir is None:
            return

        self.pan_neo_tab = PanCoreAccessoryTab(
            self.output_dir,
            self.meta_df,
            reference_sample=self._snippy_reference_sample,
            fallback_tree=self._resolve_raxml_tree_path(),
            ann_df=self.ann_df,
            preserve_branch_lengths=True,
            color_tree_by_fastgear_clade=True,
            optimize_tree_by_fastgear_clade=True,
            annotate_fastgear_clusters=True,
            use_fastgear_clade_backgrounds=True,
            enable_view_selector=True,
            virulence_db_names=[db for db in self.annotation_db_columns if self._tree_db_is_virulence_source(db)],
            amr_present_color=self._neo_tree_amr_color,
            mge_present_color=self._neo_tree_mge_color,
            virulence_present_color=self._neo_tree_virulence_color,
        )
        self.pan_neo_tab.sampleHovered.connect(self._on_pan_sample_hovered)
        self.pan_neo_tab.cladeFilterRequested.connect(self._on_pan_clade_filter_requested)
        self.pan_neo_tab.clusterFilterRequested.connect(self._on_pan_cluster_filter_requested)
        self.pan_neo_tab.plasmidFilterRequested.connect(self._on_pan_plasmid_filter_requested)
        self.pan_neo_tab.plasmidContextRequested.connect(self._on_pan_plasmid_context_requested)
        self.pan_neo_tab.sampleOpenRequested.connect(self._on_pan_sample_open_requested)
        tab_specs: List[Tuple[str, QtWidgets.QWidget]] = [
            ("neoPan", self.pan_neo_tab),
        ]

        for tab_name, widget in tab_specs:
            tab = self._bottom_tab_by_name(tab_name)
            if tab is None:
                continue
            layout = self._prepare_tab_layout(tab)
            if layout is None:
                continue
            layout.addWidget(widget)
        self._apply_pan_visual_highlight()


    def _entropy_gene_query(self, ann_filters: Optional[Dict[str, Any]] = None) -> str:
        gene_query = ""
        if ann_filters:
            gene_query = (
                str(ann_filters.get("gene_text") or "").strip()
                or str(ann_filters.get("product_text") or "").strip()
            )
        if not gene_query and self.view_mode == "Global":
            gene_query = self.sample_query.text().strip()
        return gene_query


    def _panaroo_widgets_for_tab(self, tab_name: str) -> List[QtWidgets.QWidget]:
        widgets: List[QtWidgets.QWidget] = []
        if tab_name == "neoPan":
            if self.pan_neo_tab is not None:
                widgets.append(self.pan_neo_tab)
            return widgets
        if tab_name == "neoTree":
            if self.pan_neo_tree_tab is not None:
                widgets.append(self.pan_neo_tree_tab)
            return widgets
        if tab_name == "cgMLSTForce":
            if self._cgmlst_force_plot_panel is not None:
                widgets.append(self._cgmlst_force_plot_panel)
            return widgets
        if tab_name == "StructVar":
            if self.pan_structvar_tab is not None:
                widgets.append(self.pan_structvar_tab)
            return widgets
        if tab_name == "Entropy":
            if self.pan_entropy_tab is not None:
                widgets.append(self.pan_entropy_tab)
            return widgets
        return widgets


    def _refresh_panaroo_tab(
        self,
        tab_name: str,
        allowed_samples: Sequence[str],
        ann_filters: Optional[Dict[str, Any]] = None,
    ) -> bool:
        widgets = self._panaroo_widgets_for_tab(tab_name)
        if not widgets:
            return False
        if tab_name == "cgMLSTForce":
            lineage_by_sample, cluster_by_sample = self._load_cgmlst_fastgear_maps()
            hover_metadata = self._build_plot_sample_hover_metadata()
            for widget in widgets:
                if not isinstance(widget, CgmlstForcePlotPanel):
                    continue
                widget.set_reference_sample(self._snippy_reference_sample, refresh=False)
                widget.set_fastgear_assignments(lineage_by_sample, {}, cluster_by_sample)
                widget.set_sample_hover_metadata(hover_metadata)
                widget.set_trait_values(self._build_cgmlst_trait_values())
                widget.set_sample_filter(list(allowed_samples))
                widget.ensure_fresh()
            return True
        for widget in widgets:
            widget.set_sample_filter(allowed_samples)
            if tab_name == "Entropy" and hasattr(widget, "set_gene_query"):
                widget.set_gene_query(self._entropy_gene_query(ann_filters))
        return True


    def _update_panaroo_sample_scoped_views(
        self,
        allowed_samples: Sequence[str],
        ann_filters: Optional[Dict[str, Any]] = None,
    ) -> None:
        for tab_name in self._PANAROO_TAB_NAMES:
            self._refresh_panaroo_tab(tab_name, allowed_samples, ann_filters)
        for plot_name in self._PLOT_WIDGET_NAMES:
            self._refresh_panaroo_tab(plot_name, allowed_samples, ann_filters)


    def _schedule_deferred_panel_refresh(
        self,
        allowed_samples: Sequence[str],
        ann_filters: Optional[Dict[str, Any]] = None,
        *,
        immediate: bool = False,
    ) -> None:
        self._deferred_panel_allowed_samples = self._ordered_sample_names(allowed_samples)
        self._deferred_panel_ann_filters = dict(ann_filters or {})
        self._deferred_geo_dirty = True
        self._deferred_pan_tab_dirty.update(self._PANAROO_TAB_NAMES)
        self._deferred_pan_tab_dirty.update(self._PLOT_WIDGET_NAMES)
        if self._deferred_panel_timer.isActive():
            self._deferred_panel_timer.stop()
        if immediate:
            self._flush_deferred_panel_refresh()
            return
        self._deferred_panel_timer.start()


    def _flush_deferred_panel_refresh(self) -> None:
        if self._deferred_panel_timer.isActive():
            self._deferred_panel_timer.stop()

        updated = False
        allowed_samples = list(self._deferred_panel_allowed_samples)
        ann_filters = dict(self._deferred_panel_ann_filters)

        active_tab = self._active_bottom_tab_name()
        if active_tab in self._GEO_TAB_NAMES or active_tab == "neoPan":
            if self._deferred_geo_dirty:
                self._update_map()
                self._deferred_geo_dirty = False
                updated = True
        if active_tab in self._PANAROO_TAB_NAMES and active_tab in self._deferred_pan_tab_dirty:
            if self._refresh_panaroo_tab(active_tab, allowed_samples, ann_filters):
                updated = True
            self._deferred_pan_tab_dirty.discard(active_tab)
        for plot_name in self._PLOT_WIDGET_NAMES:
            if plot_name not in self._deferred_pan_tab_dirty:
                continue
            if plot_name == "neoTree":
                dlg = self._phandango_plot_dialog
                if dlg is None or not dlg.isVisible():
                    continue
            if plot_name == "cgMLSTForce":
                dlg = self._cgmlst_force_plot_dialog
                if dlg is None or not dlg.isVisible():
                    continue
            if self._refresh_panaroo_tab(plot_name, allowed_samples, ann_filters):
                updated = True
            self._deferred_pan_tab_dirty.discard(plot_name)

        if updated:
            self._apply_pan_visual_highlight()

    def _apply_pan_subset_filter_request(self, payload: object, *, kind: str) -> None:
        if not isinstance(payload, dict):
            return
        name = str(payload.get("name") or "").strip()
        samples = payload.get("samples")
        if not name or not isinstance(samples, (list, tuple)):
            return
        label = f"{kind.title()} {name}"
        self._set_sample_subset_filter(label=label, samples=samples, kind=kind)

    def _on_pan_clade_filter_requested(self, payload: object) -> None:
        self._apply_pan_subset_filter_request(payload, kind="clade")

    def _on_pan_cluster_filter_requested(self, payload: object) -> None:
        self._apply_pan_subset_filter_request(payload, kind="cluster")

    def _on_pan_plasmid_filter_requested(self, payload: object) -> None:
        self._apply_pan_subset_filter_request(payload, kind="plasmid")

    def _on_pan_plasmid_context_requested(self, payload: object) -> None:
        if not isinstance(payload, dict):
            return
        rows = payload.get("rows")
        if not isinstance(rows, list) or not rows:
            QtWidgets.QMessageBox.information(
                self,
                "Plasmid context",
                "No plasmid instances were available for this cluster.",
            )
            return
        name = str(payload.get("name") or "Plasmid").strip() or "Plasmid"
        dlg = IgvCollectionDialog(
            self,
            rows=rows,
            title=f"Plasmid context | {name}",
            output_dir=self.output_dir,
            pan_reference_sample=self._snippy_reference_sample,
            pan_fallback_tree=self._resolve_raxml_tree_path(),
        )
        dlg.geneDoubleClicked.connect(self._on_igv_gene_double_clicked)
        self._show_floating_dialog(dlg)

    def _on_pan_sample_open_requested(self, sample: object) -> None:
        value = str(sample).strip() if sample is not None else ""
        if value:
            self._open_sample_card_for_sample(value)



    def _resolve_raxml_tree_path(self) -> Optional[Path]:
        if self.output_dir is None:
            return None
        if not self.output_dir.exists():
            return None
        direct_candidates = [
            self.output_dir / "trees" / "snippy_raxml.tree",
            self.output_dir / "trees" / "snippy_raxml.newick",
            self.output_dir / "snippy" / "core" / "snippy_raxml.tree",
            self.output_dir / "trees" / "core_gene_tree.nwk",
        ]
        for path in direct_candidates:
            if path.exists() and path.is_file():
                return path

        patterns = [
            "*raxml*.tree",
            "*raxml*.newick",
            "*RAxML*.tree",
            "*RAxML*.newick",
            "*RAxML_bestTree*",
            "*raxml_besttree*",
        ]
        matches: List[Path] = []
        for pattern in patterns:
            matches.extend(p for p in self.output_dir.rglob(pattern) if p.is_file())
        if not matches:
            return None

        def _score(path: Path) -> Tuple[int, float]:
            name = path.name.lower()
            parent = str(path.parent).lower()
            bonus = 0
            if "raxml" in name:
                bonus += 8
            if "snippy" in name or "snippy" in parent:
                bonus += 4
            if name == "core_gene_tree.nwk":
                bonus += 1
            if "tree" in name:
                bonus += 1
            try:
                mtime = path.stat().st_mtime
            except Exception:
                mtime = 0.0
            return bonus, mtime

        matches.sort(key=_score, reverse=True)
        return matches[0]


    def _build_map_payload(self) -> dict:
        if self._filtered_ann_pl is None or self._filtered_ann_pl.is_empty():
            return {"points": [], "countries": {}, "maxCountryCount": 1}

        if "Sample" not in self._filtered_ann_pl.columns or "Sample" not in self._meta_pl.columns:
            return {"points": [], "countries": {}, "maxCountryCount": 1}

        samples = self._collect_sample_list(self._filtered_ann_pl)
        if not samples:
            return {"points": [], "countries": {}, "maxCountryCount": 1}
        meta_columns = [str(column) for column in self._meta_pl.columns]
        lat_col = self._resolve_plot_metadata_column(meta_columns, self._LATITUDE_TOOLTIP_CANDIDATES)
        lon_col = self._resolve_plot_metadata_column(meta_columns, self._LONGITUDE_TOOLTIP_CANDIDATES)
        country_col = self._resolve_plot_metadata_column(meta_columns, self._COUNTRY_TOOLTIP_CANDIDATES)
        meta = (
            self._meta_pl.lazy()
            .with_columns(
                [
                    pl.col("Sample").cast(pl.Utf8).str.strip_chars().alias("Sample"),
                    (
                        pl.col(lat_col).cast(pl.Float64, strict=False)
                        if lat_col
                        else pl.lit(None, dtype=pl.Float64)
                    ).alias("_lat"),
                    (
                        pl.col(lon_col).cast(pl.Float64, strict=False)
                        if lon_col
                        else pl.lit(None, dtype=pl.Float64)
                    ).alias("_lon"),
                    (
                        pl.col(country_col).cast(pl.Utf8).fill_null("").str.strip_chars()
                        if country_col
                        else pl.lit("", dtype=pl.Utf8)
                    ).alias("_country"),
                ]
            )
            .filter(pl.col("Sample").is_in(samples))
            .unique(subset=["Sample"], keep="first")
            .select(["Sample", "_lat", "_lon", "_country"])
            .collect()
        )
        if meta.is_empty():
            return {"points": [], "countries": {}, "maxCountryCount": 1}

        points: List[Dict[str, Any]] = []
        countries: Dict[str, int] = {}

        for row in meta.iter_rows(named=True):
            sample = str(row.get("Sample") or "").strip()
            if not sample:
                continue
            count = 1
            lat = row.get("_lat")
            lon = row.get("_lon")
            has_coords = (
                isinstance(lat, (float, int))
                and isinstance(lon, (float, int))
                and not (isinstance(lat, float) and np.isnan(lat))
                and not (isinstance(lon, float) and np.isnan(lon))
            )
            if has_coords:
                points.append(
                    {
                        "sample": sample,
                        "lat": float(lat),
                        "lon": float(lon),
                        "count": count,
                        "country": str(row.get("_country") or ""),
                    }
                )
            else:
                country = str(row.get("_country") or "").strip()
                if country:
                    countries[country] = countries.get(country, 0) + 1

        max_count = max(countries.values(), default=1)
        return {"points": points, "countries": countries, "maxCountryCount": max_count}


    def _update_map(self) -> None:
        payload = self._build_map_payload()
        geo_views = self._all_geo_views()
        if not geo_views:
            return
        for view in geo_views:
            view.update_data(
                points=payload.get("points", []),
                countries=payload.get("countries", {}),
                max_count=payload.get("maxCountryCount", 1),
            )
