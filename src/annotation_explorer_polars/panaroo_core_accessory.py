from __future__ import annotations

import html
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets

from .geo_experimental_views import GeoMapLibreView
from .pan_view import PanHeatmapWidget
from .panaroo_force import CgmlstForceGraphView
from .tree_view import RaxmlTreeView
from .panaroo_common import (
    _NEO_TREE_MGE_TOKEN,
    _NEO_TREE_PANAROO_TOKEN,
    _NEO_TREE_PLASMID_TOKEN,
    _NEO_TREE_RESISTANCE_TOKEN,
    _NEO_TREE_VIRULENCE_TOKEN,
    _file_signature,
    _natural_sort_key,
)
from .plasmid_model import PlasmidClusterIndex, build_plasmid_cluster_index
from .plasmidome_view import PlasmidomePanel
from .panaroo_data import (
    _extract_mge_labels,
    _extract_virulence_labels,
    _load_core_alignment_cached,
    _load_presence_matrix,
    _build_subset_tree_newick_cached,
    _extract_resistance_gene_labels,
    _find_fastgear_sample_clades,
    _find_panaroo_out,
    _find_panaroo_alignment,
    _load_fastgear_clusters,
    _load_fastgear_lineages,
    _open_text_auto,
    _replace_newick_leaf_tokens,
)
from .panaroo_widgets import (
    _CategoricalLegendPanel,
    _CategoricalMetadataStripPanel,
    _DeferredViewPlaceholder,
    _MetadataDonutPanel,
    _StaticLegendPanel,
    _high_contrast_color,
)
from .panaroo_style import (
    _apply_analysis_panel_style,
    _fastgear_cluster_colors,
    _fastgear_lineage_colors,
    _log_view_event,
)


class PanCoreAccessoryTab(QtWidgets.QWidget):
    sampleHovered = QtCore.pyqtSignal(object)
    cladeFilterRequested = QtCore.pyqtSignal(object)
    clusterFilterRequested = QtCore.pyqtSignal(object)
    plasmidFilterRequested = QtCore.pyqtSignal(object)
    plasmidContextRequested = QtCore.pyqtSignal(object)
    sampleOpenRequested = QtCore.pyqtSignal(object)

    _COUNTRY_TOOLTIP_CANDIDATES = ("Country", "_country")
    _ST_TOOLTIP_CANDIDATES = (
        "ST",
        "Sequence Type",
        "SequenceType",
        "sequence_type",
        "MLST ST",
        "MLST_ST",
    )
    _COLLECTION_DATE_TOOLTIP_CANDIDATES = (
        "Collection date",
        "Collection Date",
        "Collection_date",
        "collection_date",
        "CollectionDate",
        "Date",
        "date",
    )

    def __init__(
        self,
        output_dir: Path,
        metadata_df: Optional[pd.DataFrame],
        reference_sample: Optional[str],
        fallback_tree: Optional[Path],
        *,
        preserve_branch_lengths: bool = False,
        color_tree_by_fastgear_clade: bool = False,
        optimize_tree_by_fastgear_clade: bool = False,
        annotate_fastgear_clusters: bool = False,
        use_fastgear_clade_backgrounds: bool = False,
        color_tree_by_fastgear_cluster: bool = False,
        show_clade_legend: bool = False,
        show_cluster_legend: bool = False,
        show_heatmap: bool = True,
        enable_view_selector: bool = False,
        show_info_label: bool = True,
        enable_tree_hover_interactions: bool = True,
        side_metadata_title: Optional[str] = None,
        side_metadata_candidates: Optional[Sequence[str]] = None,
        side_metadata_columns: Optional[Sequence[str]] = None,
        marker_hits_df: Optional[pd.DataFrame] = None,
        ann_df: Optional[pd.DataFrame] = None,
        resistance_db_names: Optional[Sequence[str]] = None,
        mge_db_names: Optional[Sequence[str]] = None,
        virulence_db_names: Optional[Sequence[str]] = None,
        amr_present_color: Optional[object] = None,
        mge_present_color: Optional[object] = None,
        virulence_present_color: Optional[object] = None,
        recompute_tree_from_sample_filter: bool = False,
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self._output_dir = Path(output_dir)
        self._metadata_df = metadata_df.copy() if metadata_df is not None else pd.DataFrame()
        self._sample_hover_metadata = self._build_sample_hover_metadata()
        self._reference_sample = reference_sample
        self._fallback_tree = fallback_tree
        self._last_geo_points: List[dict] = []
        self._last_geo_countries: Dict[str, int] = {}
        self._last_geo_max_count = 1
        self._preserve_branch_lengths = bool(preserve_branch_lengths)
        self._color_tree_by_fastgear_clade = bool(color_tree_by_fastgear_clade)
        self._optimize_tree_by_fastgear_clade = bool(optimize_tree_by_fastgear_clade)
        self._annotate_fastgear_clusters = bool(annotate_fastgear_clusters)
        self._use_fastgear_clade_backgrounds = bool(use_fastgear_clade_backgrounds)
        self._color_tree_by_fastgear_cluster = bool(color_tree_by_fastgear_cluster)
        self._show_clade_legend = bool(show_clade_legend)
        self._show_cluster_legend = bool(show_cluster_legend)
        self._show_heatmap = bool(show_heatmap)
        self._enable_view_selector = bool(enable_view_selector)
        self._show_info_label = bool(show_info_label)
        self._enable_tree_hover_interactions = bool(enable_tree_hover_interactions)
        self._side_metadata_title = str(side_metadata_title or "").strip()
        self._side_metadata_candidates = tuple(
            str(candidate).strip()
            for candidate in (side_metadata_candidates or [])
            if str(candidate).strip()
        )
        self._marker_hits_df = marker_hits_df.copy() if marker_hits_df is not None else pd.DataFrame()
        self._ann_df = ann_df.copy() if ann_df is not None else pd.DataFrame()
        self._resistance_db_names = {
            str(name).strip()
            for name in (resistance_db_names or [])
            if str(name).strip()
        }
        self._mge_db_names = {
            str(name).strip()
            for name in (mge_db_names or [])
            if str(name).strip()
        }
        self._virulence_db_names = {
            str(name).strip()
            for name in (virulence_db_names or [])
            if str(name).strip()
        }
        self._amr_present_color = self._coerce_present_color(amr_present_color, "#f97316")
        self._mge_present_color = self._coerce_present_color(mge_present_color, "#0f766e")
        self._virulence_present_color = self._coerce_present_color(virulence_present_color, "#be185d")
        self._plasmid_cluster_index: PlasmidClusterIndex = build_plasmid_cluster_index(
            output_dir=self._output_dir,
            ann_df=self._ann_df,
        )
        self._side_metadata_columns = self._normalize_side_metadata_columns(side_metadata_columns)
        if not self._side_metadata_columns and self._side_metadata_candidates:
            resolved = _CategoricalMetadataStripPanel._resolve_column(
                list(self._metadata_df.columns),
                self._side_metadata_candidates,
            )
            if resolved:
                self._side_metadata_columns = [resolved]
        self._recompute_tree_from_sample_filter = bool(recompute_tree_from_sample_filter)
        self._sample_filter: Optional[set[str]] = None
        self._loaded_tree_key: Optional[Tuple[object, ...]] = None
        self._loaded_sample_clades_sig: Optional[Tuple[str, int, int]] = None
        self._lineage_by_sample: Dict[str, str] = {}
        self._lineage_color_map: Dict[str, QtGui.QColor] = {}
        self._cluster_by_sample: Dict[str, str] = {}
        self.tree: Optional[RaxmlTreeView] = None
        self.heatmap: Optional[PanHeatmapWidget] = None
        self.geo_view: Optional[GeoMapLibreView] = None
        self.metadata_view: Optional[_MetadataDonutPanel] = None
        self.force_view: Optional[CgmlstForceGraphView] = None
        self.plasmidome_view: Optional[PlasmidomePanel] = None
        self._tree_views: List[RaxmlTreeView] = []
        self._heatmap_views: List[PanHeatmapWidget] = []
        self._geo_views: List[GeoMapLibreView] = []
        self._metadata_views: List[_MetadataDonutPanel] = []
        self._force_views: List[CgmlstForceGraphView] = []
        self._plasmidome_views: List[PlasmidomePanel] = []
        self.side_metadata_views: List[_CategoricalMetadataStripPanel] = []
        self._side_metadata_view_by_column: Dict[str, _CategoricalMetadataStripPanel] = {}
        self.clade_legend_view: Optional[_CategoricalLegendPanel] = None
        self.cluster_legend_view: Optional[_CategoricalLegendPanel] = None
        self.resistance_legend_view: Optional[_StaticLegendPanel] = None
        self.mge_legend_view: Optional[_StaticLegendPanel] = None
        self.virulence_legend_view: Optional[_StaticLegendPanel] = None
        self._legend_layout: Optional[QtWidgets.QHBoxLayout] = None
        self._legend_bar: Optional[QtWidgets.QWidget] = None
        self._detail_view_combo: Optional[QtWidgets.QComboBox] = None
        self._detail_stack: Optional[QtWidgets.QStackedWidget] = None
        self._info: Optional[QtWidgets.QLabel] = None
        self._export_button: Optional[QtWidgets.QPushButton] = None
        self._resistance_gene_order_cache: Optional[List[str]] = None
        self._resistance_presence_cache: Optional[Dict[str, set[str]]] = None
        self._mge_presence_cache: Optional[Dict[str, set[str]]] = None
        self._virulence_presence_cache: Optional[Dict[str, set[str]]] = None
        self._panaroo_presence_matrix: Optional[object] = None
        self._panaroo_presence_genes: List[str] = []
        self._panaroo_presence_products: List[str] = []
        self._panaroo_presence_samples: List[str] = []
        self._panaroo_sample_index: Dict[str, int] = {}
        self._pending_hover_samples: Tuple[str, ...] = tuple()
        self._pending_hover_emit_sample: Optional[str] = None
        self._pending_hover_source_tree: Optional[RaxmlTreeView] = None
        self._last_applied_hover_samples: Optional[Tuple[str, ...]] = None
        self._hover_update_timer = QtCore.QTimer(self)
        self._hover_update_timer.setInterval(24)
        self._hover_update_timer.setSingleShot(True)
        self._hover_update_timer.timeout.connect(self._flush_hover_highlight)

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        if self._show_info_label:
            self._info = QtWidgets.QLabel()
            self._info.setWordWrap(True)
            self._info.setStyleSheet("color: #607183; font-weight: 600;")
            self._info.setVisible(False)
            layout.addWidget(self._info)
        _log_view_event(f"created Panaroo accessory tab for {self._output_dir}")

        if self._show_heatmap:
            split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal, self)
            split.setChildrenCollapsible(False)
            split.setHandleWidth(8)
            layout.addWidget(split, stretch=1)
            if self._enable_view_selector:
                left_panel, left_tree, left_heatmap, left_geo, left_metadata, left_force, _left_combo, _left_stack = self._create_selector_panel(
                    split,
                    default_view="tree",
                )
                right_panel, right_tree, right_heatmap, right_geo, right_metadata, right_force, right_combo, right_stack = self._create_selector_panel(
                    split,
                    default_view="presence",
                )
                self.tree = left_tree
                self.heatmap = right_heatmap
                self.geo_view = right_geo
                self.metadata_view = right_metadata
                self.force_view = right_force
                self._detail_view_combo = right_combo
                self._detail_stack = right_stack
                split.addWidget(left_panel)
                split.addWidget(right_panel)
            else:
                tree = self._create_tree_view(split)
                self.tree = tree
                split.addWidget(tree)
                heatmap = PanHeatmapWidget(self._output_dir, parent=split)
                split.addWidget(heatmap)
                self.heatmap = heatmap
                self._heatmap_views.append(heatmap)
                heatmap.sampleHovered.connect(self._on_heatmap_sample_hover)
            split.setStretchFactor(0, 0)
            split.setStretchFactor(1, 1)
            split.setSizes([520, 860])
        else:
            tree = self._create_tree_view(self)
            self.tree = tree
            tree_panel = QtWidgets.QWidget(self)
            tree_panel_layout = QtWidgets.QVBoxLayout(tree_panel)
            tree_panel_layout.setContentsMargins(0, 0, 0, 0)
            tree_panel_layout.setSpacing(6)

            controls_row = QtWidgets.QHBoxLayout()
            controls_row.setContentsMargins(0, 0, 0, 0)
            controls_row.setSpacing(6)
            controls_row.addStretch(1)

            export_button = QtWidgets.QPushButton("Export Figure…", tree_panel)
            export_button.setToolTip("Save the current neoTree figure as a high-quality PDF, PNG, or JPEG.")
            export_button.clicked.connect(self._export_tree_figure)
            controls_row.addWidget(export_button, 0, QtCore.Qt.AlignmentFlag.AlignRight)
            tree_panel_layout.addLayout(controls_row)
            self._export_button = export_button

            tree_panel_layout.addWidget(tree, stretch=1)

            layout.addWidget(tree_panel, stretch=1)

        self.refresh()

    @staticmethod
    def _coerce_present_color(value: Optional[object], fallback: str) -> QtGui.QColor:
        if isinstance(value, QtGui.QColor):
            color = QtGui.QColor(value)
        else:
            color = QtGui.QColor(str(value or ""))
        if not color.isValid():
            color = QtGui.QColor(fallback)
        return color

    def set_special_track_colors(
        self,
        *,
        amr_color: Optional[object] = None,
        mge_color: Optional[object] = None,
        virulence_color: Optional[object] = None,
    ) -> None:
        changed = False
        if amr_color is not None:
            next_amr = self._coerce_present_color(amr_color, "#f97316")
            if next_amr != self._amr_present_color:
                self._amr_present_color = next_amr
                changed = True
        if mge_color is not None:
            next_mge = self._coerce_present_color(mge_color, "#0f766e")
            if next_mge != self._mge_present_color:
                self._mge_present_color = next_mge
                changed = True
        if virulence_color is not None:
            next_virulence = self._coerce_present_color(virulence_color, "#be185d")
            if next_virulence != self._virulence_present_color:
                self._virulence_present_color = next_virulence
                changed = True
        if not changed:
            return
        self._rebuild_side_metadata_legend_views()
        self._apply_side_metadata_state(refresh_views=True)

    def _all_tree_views(self) -> List[RaxmlTreeView]:
        views: List[RaxmlTreeView] = []
        seen: set[int] = set()
        for view in [self.tree, *self._tree_views]:
            if view is None:
                continue
            ident = id(view)
            if ident in seen:
                continue
            seen.add(ident)
            views.append(view)
        return views

    def _all_heatmaps(self) -> List[PanHeatmapWidget]:
        views: List[PanHeatmapWidget] = []
        seen: set[int] = set()
        for view in [self.heatmap, *self._heatmap_views]:
            if view is None:
                continue
            ident = id(view)
            if ident in seen:
                continue
            seen.add(ident)
            views.append(view)
        return views

    def all_geo_views(self) -> List[GeoMapLibreView]:
        views: List[GeoMapLibreView] = []
        seen: set[int] = set()
        for view in [self.geo_view, *self._geo_views]:
            if view is None:
                continue
            ident = id(view)
            if ident in seen:
                continue
            seen.add(ident)
            views.append(view)
        return views

    def sample_hover_metadata(self) -> Dict[str, Dict[str, str]]:
        return {
            sample: dict(fields)
            for sample, fields in self._sample_hover_metadata.items()
        }

    def _all_metadata_views(self) -> List[_MetadataDonutPanel]:
        views: List[_MetadataDonutPanel] = []
        seen: set[int] = set()
        for view in [self.metadata_view, *self._metadata_views]:
            if view is None:
                continue
            ident = id(view)
            if ident in seen:
                continue
            seen.add(ident)
            views.append(view)
        return views

    def _all_force_views(self) -> List[CgmlstForceGraphView]:
        views: List[CgmlstForceGraphView] = []
        seen: set[int] = set()
        for view in [self.force_view, *self._force_views]:
            if view is None:
                continue
            ident = id(view)
            if ident in seen:
                continue
            seen.add(ident)
            views.append(view)
        return views

    def _all_plasmidome_views(self) -> List[PlasmidomePanel]:
        views: List[PlasmidomePanel] = []
        seen: set[int] = set()
        for view in [self.plasmidome_view, *self._plasmidome_views]:
            if view is None:
                continue
            ident = id(view)
            if ident in seen:
                continue
            seen.add(ident)
            views.append(view)
        return views

    @staticmethod
    def _visible_views(views: Sequence[QtWidgets.QWidget]) -> List[QtWidgets.QWidget]:
        return [view for view in views if view is not None and view.isVisible()]

    @staticmethod
    def _normalize_metadata_value(value: object) -> str:
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

    def _load_mlst_sequence_types(self) -> Dict[str, str]:
        path = self._output_dir / "mlst" / "results.csv"
        if not path.is_file():
            return {}
        try:
            df = pd.read_csv(path, dtype=str, low_memory=False)
        except Exception:
            return {}
        required = {"sample", "sequence_type"}
        if not required.issubset({str(column).strip() for column in df.columns}):
            return {}
        sample_col = _CategoricalMetadataStripPanel._resolve_column(list(df.columns), ("sample",))
        st_col = _CategoricalMetadataStripPanel._resolve_column(list(df.columns), ("sequence_type",))
        if not sample_col or not st_col:
            return {}
        sequence_types: Dict[str, str] = {}
        for row in df[[sample_col, st_col]].itertuples(index=False):
            sample = self._normalize_metadata_value(getattr(row, sample_col))
            st_value = self._normalize_metadata_value(getattr(row, st_col))
            if sample and st_value and sample not in sequence_types:
                sequence_types[sample] = st_value
        return sequence_types

    def _build_sample_hover_metadata(self) -> Dict[str, Dict[str, str]]:
        metadata: Dict[str, Dict[str, str]] = {}
        if not self._metadata_df.empty and "Sample" in self._metadata_df.columns:
            columns = list(self._metadata_df.columns)
            country_col = _CategoricalMetadataStripPanel._resolve_column(columns, self._COUNTRY_TOOLTIP_CANDIDATES)
            st_col = _CategoricalMetadataStripPanel._resolve_column(columns, self._ST_TOOLTIP_CANDIDATES)
            date_col = _CategoricalMetadataStripPanel._resolve_column(columns, self._COLLECTION_DATE_TOOLTIP_CANDIDATES)
            for _idx, row in self._metadata_df.iterrows():
                sample = self._normalize_metadata_value(row.get("Sample"))
                if not sample:
                    continue
                fields = metadata.setdefault(sample, {})
                if country_col and "Country" not in fields:
                    value = self._normalize_metadata_value(row.get(country_col))
                    if value:
                        fields["Country"] = value
                if st_col and "ST" not in fields:
                    value = self._normalize_metadata_value(row.get(st_col))
                    if value:
                        fields["ST"] = value
                if date_col and "Collection date" not in fields:
                    value = self._normalize_metadata_value(row.get(date_col))
                    if value:
                        fields["Collection date"] = value

        for sample, sequence_type in self._load_mlst_sequence_types().items():
            if not sequence_type:
                continue
            metadata.setdefault(sample, {}).setdefault("ST", sequence_type)
        return {sample: fields for sample, fields in metadata.items() if fields}

    def _create_tree_view(self, parent: Optional[QtWidgets.QWidget]) -> RaxmlTreeView:
        tree = RaxmlTreeView(
            parent,
            layout_mode="rectangular",
            preserve_branch_lengths=self._preserve_branch_lengths,
            optimize_groups=self._optimize_tree_by_fastgear_clade,
            use_group_backgrounds=self._use_fastgear_clade_backgrounds,
            show_group_labels=not self._show_clade_legend,
            show_cluster_annotations=self._annotate_fastgear_clusters,
            enable_hover_interactions=self._enable_tree_hover_interactions,
        )
        tree.setMinimumWidth(280)
        tree.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Preferred,
            QtWidgets.QSizePolicy.Policy.Expanding,
        )
        tree.set_reference_sample(self._reference_sample)
        tree.set_leaf_hover_metadata(self._sample_hover_metadata)
        tree.cladeClicked.connect(self._on_clade_clicked)
        tree.clusterClicked.connect(self._on_cluster_clicked)
        tree.leafDoubleClicked.connect(self._on_tree_leaf_double_clicked)
        if self._enable_tree_hover_interactions:
            tree.leafHovered.connect(self._on_tree_leaf_hover)
            tree.subtreeHovered.connect(self._on_tree_subtree_hover)
        self._tree_views.append(tree)
        return tree

    @staticmethod
    def _first_stack_view(
        stack: QtWidgets.QStackedWidget,
        view_type: type,
    ) -> Optional[QtWidgets.QWidget]:
        for idx in range(stack.count()):
            widget = stack.widget(idx)
            if isinstance(widget, view_type):
                return widget
        return None

    @staticmethod
    def _lazy_placeholder(view_type: str, parent: QtWidgets.QWidget) -> _DeferredViewPlaceholder:
        labels = {
            "tree": ("Tree", "Select this view to build the tree panel."),
            "presence": ("Gene absence/presence", "Select this view to load the Panaroo presence heatmap."),
            "geo": ("Geo", "Select this view to initialize the geospatial map."),
            "metadata": ("Metadata", "Select this view to build the metadata summaries."),
            "force": ("cgMLST Force", "Select this view to build the animated cgMLST force graph."),
            "plasmidome": ("Plasmidome", "Select this view to list MOB-suite plasmid clusters."),
        }
        title, message = labels.get(view_type, ("View", "Select this view to initialize it."))
        placeholder = _DeferredViewPlaceholder(title, message, parent=parent)
        placeholder.setProperty("lazy_view_type", view_type)
        return placeholder

    def _create_selector_panel(
        self,
        parent: QtWidgets.QWidget,
        *,
        default_view: str,
    ) -> Tuple[
        QtWidgets.QWidget,
        Optional[RaxmlTreeView],
        Optional[PanHeatmapWidget],
        Optional[GeoMapLibreView],
        Optional[_MetadataDonutPanel],
        Optional[CgmlstForceGraphView],
        QtWidgets.QComboBox,
        QtWidgets.QStackedWidget,
    ]:
        detail_panel = QtWidgets.QFrame(parent)
        _apply_analysis_panel_style(detail_panel)
        detail_layout = QtWidgets.QVBoxLayout(detail_panel)
        detail_layout.setContentsMargins(0, 0, 0, 0)
        detail_layout.setSpacing(0)

        selector_bar = QtWidgets.QFrame(detail_panel)
        selector_bar.setObjectName("NeoPanViewSelectorBar")
        selector_bar.setAttribute(QtCore.Qt.WidgetAttribute.WA_StyledBackground, True)
        selector_bar.setStyleSheet(
            "#NeoPanViewSelectorBar {"
            "  background: #f8fafc;"
            "  border-bottom: 1px solid #d8e2ed;"
            "  border-top-left-radius: 8px;"
            "  border-top-right-radius: 8px;"
            "}"
        )
        selector_layout = QtWidgets.QHBoxLayout(selector_bar)
        selector_layout.setContentsMargins(12, 8, 12, 8)
        selector_layout.setSpacing(8)
        selector_label = QtWidgets.QLabel("View", selector_bar)
        selector_label.setStyleSheet("color: #314154; font-weight: 700;")
        selector_layout.addWidget(selector_label)

        detail_stack = QtWidgets.QStackedWidget(detail_panel)
        detail_combo = QtWidgets.QComboBox(selector_bar)
        detail_combo.addItem("Tree", "tree")
        detail_combo.addItem("Gene absence/presence", "presence")
        detail_combo.addItem("Geo", "geo")
        detail_combo.addItem("Metadata", "metadata")
        detail_combo.addItem("cgMLST Force", "force")
        detail_combo.addItem("Plasmidome", "plasmidome")
        detail_combo.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToContents)
        detail_combo.setMinimumWidth(180)
        detail_combo.setMaximumWidth(240)
        detail_combo.currentIndexChanged.connect(
            lambda index, stack=detail_stack: self._on_detail_view_changed(stack, index)
        )
        selector_layout.addWidget(detail_combo)
        selector_layout.addStretch(1)
        detail_layout.addWidget(selector_bar)
        detail_layout.addWidget(detail_stack, stretch=1)

        for view_type in ("tree", "presence", "geo", "metadata", "force", "plasmidome"):
            detail_stack.addWidget(self._lazy_placeholder(view_type, detail_stack))

        default_key = str(default_view or "tree").strip().lower()
        selected_index = 0
        for idx in range(detail_combo.count()):
            if str(detail_combo.itemData(idx) or "").strip().lower() == default_key:
                selected_index = idx
                break
        detail_combo.setCurrentIndex(selected_index)
        detail_stack.setCurrentIndex(selected_index)
        created = self._ensure_detail_view_ready(detail_stack, selected_index)
        tree = created if isinstance(created, RaxmlTreeView) else self._first_stack_view(detail_stack, RaxmlTreeView)
        heatmap = created if isinstance(created, PanHeatmapWidget) else self._first_stack_view(detail_stack, PanHeatmapWidget)
        geo_view = created if isinstance(created, GeoMapLibreView) else self._first_stack_view(detail_stack, GeoMapLibreView)
        metadata_view = created if isinstance(created, _MetadataDonutPanel) else self._first_stack_view(detail_stack, _MetadataDonutPanel)
        force_view = created if isinstance(created, CgmlstForceGraphView) else self._first_stack_view(detail_stack, CgmlstForceGraphView)
        return detail_panel, tree, heatmap, geo_view, metadata_view, force_view, detail_combo, detail_stack

    def _normalize_side_metadata_columns(self, columns: Optional[Sequence[str]]) -> List[str]:
        available = {str(column).strip() for column in self._metadata_df.columns if str(column).strip() and str(column) != "Sample"}
        normalized: List[str] = []
        seen: set[str] = set()
        for raw_column in columns or []:
            column = str(raw_column).strip()
            if not column or column in seen:
                continue
            if column not in {
                _NEO_TREE_RESISTANCE_TOKEN,
                _NEO_TREE_MGE_TOKEN,
                _NEO_TREE_VIRULENCE_TOKEN,
                _NEO_TREE_PANAROO_TOKEN,
                _NEO_TREE_PLASMID_TOKEN,
            } and column not in available:
                continue
            if column == _NEO_TREE_RESISTANCE_TOKEN and not self._resistance_gene_order():
                continue
            if column == _NEO_TREE_MGE_TOKEN and not self._mge_element_order():
                continue
            if column == _NEO_TREE_VIRULENCE_TOKEN and not self._virulence_gene_order():
                continue
            if column == _NEO_TREE_PANAROO_TOKEN and not self._panaroo_track_available():
                continue
            if column == _NEO_TREE_PLASMID_TOKEN and not self._plasmid_cluster_index:
                continue
            seen.add(column)
            normalized.append(column)
        return normalized

    def _panaroo_track_available(self) -> bool:
        out = _find_panaroo_out(self._output_dir)
        return out is not None

    def _ensure_panaroo_presence_cache(self) -> bool:
        if self._panaroo_presence_matrix is not None and self._panaroo_presence_samples:
            return True
        out = _find_panaroo_out(self._output_dir)
        if out is None:
            return False
        try:
            presence, sample_cols, genes, products = _load_presence_matrix(out)
        except Exception:
            return False
        self._panaroo_presence_matrix = presence
        self._panaroo_presence_samples = [str(sample).strip() for sample in sample_cols]
        self._panaroo_presence_genes = [str(gene).strip() for gene in genes]
        self._panaroo_presence_products = [str(product).strip() for product in products]
        self._panaroo_sample_index = {
            sample: idx
            for idx, sample in enumerate(self._panaroo_presence_samples)
            if sample
        }
        _log_view_event(
            f"loaded Panaroo gene presence matrix: {out / 'gene_presence_absence.csv'} "
            f"({len(self._panaroo_presence_genes)} genes, {len(self._panaroo_presence_samples)} samples)"
        )
        return True

    def _resistance_rows(self) -> pd.DataFrame:
        if self._marker_hits_df.empty or "Sample" not in self._marker_hits_df.columns:
            return self._marker_hits_df.iloc[0:0].copy()
        df = self._marker_hits_df
        if self._resistance_db_names and "DB" in df.columns:
            db_values = df["DB"].fillna("").astype(str).str.strip()
            filtered = df[db_values.isin(self._resistance_db_names)].copy()
            if not filtered.empty:
                return filtered
        return df.iloc[0:0].copy()

    def _resistance_presence_by_sample(self) -> Dict[str, set[str]]:
        if self._resistance_presence_cache is not None:
            return self._resistance_presence_cache

        presence: Dict[str, set[str]] = {}
        df = self._resistance_rows()
        if df.empty:
            self._resistance_presence_cache = presence
            return presence

        for row in df.to_dict(orient="records"):
            sample = str(row.get("Sample") or "").strip()
            if not sample:
                continue
            genes = _extract_resistance_gene_labels(row)
            if not genes:
                continue
            bucket = presence.setdefault(sample, set())
            bucket.update(genes)

        self._resistance_presence_cache = presence
        return presence

    def _resistance_gene_order(self) -> List[str]:
        presence = self._resistance_presence_by_sample()
        if not presence:
            return []
        if self._sample_filter is None:
            genes = {gene for values in presence.values() for gene in values if str(gene).strip()}
        else:
            genes = {
                gene
                for sample, values in presence.items()
                if sample in self._sample_filter
                for gene in values
                if str(gene).strip()
            }
        return sorted(genes, key=_natural_sort_key)

    def _resistance_track_samples(self) -> List[str]:
        return self._track_samples_from_keys(self._resistance_presence_by_sample().keys())

    def _track_samples_from_keys(self, fallback_keys: Sequence[str]) -> List[str]:
        if self.tree is not None:
            tree_order = [str(sample).strip() for sample in self.tree.leaf_order() if str(sample).strip()]
            if tree_order:
                return tree_order
        if self._sample_filter is not None:
            return sorted(self._sample_filter)
        if not self._metadata_df.empty and "Sample" in self._metadata_df.columns:
            seen: set[str] = set()
            ordered: List[str] = []
            for raw_sample in self._metadata_df["Sample"].fillna("").astype(str).tolist():
                sample = str(raw_sample).strip()
                if not sample or sample in seen:
                    continue
                seen.add(sample)
                ordered.append(sample)
            if ordered:
                return ordered
        return sorted({str(sample).strip() for sample in fallback_keys if str(sample).strip()})

    @staticmethod
    def _build_binary_presence_tracks(
        labels: Sequence[str],
        presence: Dict[str, set[str]],
        samples: Sequence[str],
        *,
        present_color: QtGui.QColor,
    ) -> List[Dict[str, object]]:
        if not labels or not samples:
            return []
        colors = {
            "Present": QtGui.QColor(present_color),
            "Absent": QtGui.QColor("#e2e8f0"),
        }
        tracks: List[Dict[str, object]] = []
        for label in labels:
            values = {
                sample: ("Present" if label in presence.get(sample, set()) else "Absent")
                for sample in samples
            }
            tracks.append(
                {
                    "title": label,
                    "values": values,
                    "value_colors": colors,
                    "missing_label": "Absent",
                }
            )
        return tracks

    def _resistance_tracks(self) -> List[Dict[str, object]]:
        gene_order = self._resistance_gene_order()
        if not gene_order:
            return []
        presence = self._resistance_presence_by_sample()
        samples = self._resistance_track_samples()
        return self._build_binary_presence_tracks(
            gene_order,
            presence,
            samples,
            present_color=self._amr_present_color,
        )

    def _mge_rows(self) -> pd.DataFrame:
        if self._marker_hits_df.empty or "Sample" not in self._marker_hits_df.columns:
            return self._marker_hits_df.iloc[0:0].copy()
        df = self._marker_hits_df
        if self._mge_db_names and "DB" in df.columns:
            db_values = df["DB"].fillna("").astype(str).str.strip()
            filtered = df[db_values.isin(self._mge_db_names)].copy()
            if not filtered.empty:
                return filtered
        return df.iloc[0:0].copy()

    def _mge_presence_by_sample(self) -> Dict[str, set[str]]:
        if self._mge_presence_cache is not None:
            return self._mge_presence_cache

        presence: Dict[str, set[str]] = {}
        df = self._mge_rows()
        if df.empty:
            self._mge_presence_cache = presence
            return presence

        for row in df.to_dict(orient="records"):
            sample = str(row.get("Sample") or "").strip()
            if not sample:
                continue
            labels = _extract_mge_labels(row)
            if not labels:
                continue
            bucket = presence.setdefault(sample, set())
            bucket.update(labels)

        self._mge_presence_cache = presence
        return presence

    def _mge_element_order(self) -> List[str]:
        presence = self._mge_presence_by_sample()
        if not presence:
            return []
        if self._sample_filter is None:
            labels = {label for values in presence.values() for label in values if str(label).strip()}
        else:
            labels = {
                label
                for sample, values in presence.items()
                if sample in self._sample_filter
                for label in values
                if str(label).strip()
            }
        return sorted(labels, key=_natural_sort_key)

    def _mge_track_samples(self) -> List[str]:
        return self._track_samples_from_keys(self._mge_presence_by_sample().keys())

    def _mge_tracks(self) -> List[Dict[str, object]]:
        element_order = self._mge_element_order()
        if not element_order:
            return []
        presence = self._mge_presence_by_sample()
        samples = self._mge_track_samples()
        return self._build_binary_presence_tracks(
            element_order,
            presence,
            samples,
            present_color=self._mge_present_color,
        )

    def _virulence_rows(self) -> pd.DataFrame:
        if self._marker_hits_df.empty or "Sample" not in self._marker_hits_df.columns:
            return self._marker_hits_df.iloc[0:0].copy()
        df = self._marker_hits_df
        if self._virulence_db_names and "DB" in df.columns:
            db_values = df["DB"].fillna("").astype(str).str.strip()
            filtered = df[db_values.isin(self._virulence_db_names)].copy()
            if not filtered.empty:
                return filtered
        return df.iloc[0:0].copy()

    def _virulence_presence_by_sample(self) -> Dict[str, set[str]]:
        if self._virulence_presence_cache is not None:
            return self._virulence_presence_cache

        presence: Dict[str, set[str]] = {}
        df = self._virulence_rows()
        if df.empty:
            self._virulence_presence_cache = presence
            return presence

        for row in df.to_dict(orient="records"):
            sample = str(row.get("Sample") or "").strip()
            if not sample:
                continue
            labels = _extract_virulence_labels(row)
            if not labels:
                continue
            bucket = presence.setdefault(sample, set())
            bucket.update(labels)

        self._virulence_presence_cache = presence
        return presence

    def _virulence_gene_order(self) -> List[str]:
        presence = self._virulence_presence_by_sample()
        if not presence:
            return []
        if self._sample_filter is None:
            labels = {label for values in presence.values() for label in values if str(label).strip()}
        else:
            labels = {
                label
                for sample, values in presence.items()
                if sample in self._sample_filter
                for label in values
                if str(label).strip()
            }
        return sorted(labels, key=_natural_sort_key)

    def _virulence_track_samples(self) -> List[str]:
        return self._track_samples_from_keys(self._virulence_presence_by_sample().keys())

    def _virulence_tracks(self) -> List[Dict[str, object]]:
        gene_order = self._virulence_gene_order()
        if not gene_order:
            return []
        presence = self._virulence_presence_by_sample()
        samples = self._virulence_track_samples()
        return self._build_binary_presence_tracks(
            gene_order,
            presence,
            samples,
            present_color=self._virulence_present_color,
        )

    def _plasmid_tracks(self) -> List[Dict[str, object]]:
        index = self._plasmid_cluster_index
        if not index:
            return []
        samples = self._track_samples_from_keys(tuple(index.presence_by_sample.keys()))
        cluster_keys = index.cluster_order(self._sample_filter)
        if not samples or not cluster_keys:
            return []
        labels = [index.display_name_for_cluster(key) for key in cluster_keys]
        presence = {
            sample: {
                index.display_name_for_cluster(key)
                for key in index.presence_by_sample.get(sample, set())
            }
            for sample in samples
        }
        return self._build_binary_presence_tracks(
            labels,
            presence,
            samples,
            present_color=QtGui.QColor("#16a34a"),
        )

    def _panaroo_track_samples(self) -> List[str]:
        if self.tree is not None:
            tree_order = [str(sample).strip() for sample in self.tree.leaf_order() if str(sample).strip()]
            if tree_order:
                return tree_order
        if self._sample_filter is not None:
            return sorted(self._sample_filter)
        if self._panaroo_presence_samples:
            return list(self._panaroo_presence_samples)
        return []

    def _panaroo_tracks(self) -> List[Dict[str, object]]:
        if not self._ensure_panaroo_presence_cache():
            return []
        if self._panaroo_presence_matrix is None:
            return []
        samples = [sample for sample in self._panaroo_track_samples() if sample in self._panaroo_sample_index]
        if not samples:
            return []
        sample_indices = [self._panaroo_sample_index[sample] for sample in samples]
        subset = self._panaroo_presence_matrix[:, sample_indices]
        counts = subset.sum(axis=1)
        if len(samples) > 1:
            keep_mask = (counts > 0) & (counts < len(samples))
            if not bool(keep_mask.any()):
                keep_mask = counts > 0
        else:
            keep_mask = counts > 0
        kept_gene_indices = np.flatnonzero(keep_mask).tolist()
        if not kept_gene_indices:
            return []

        max_display_columns = 1200
        gene_labels = [
            str(self._panaroo_presence_genes[idx]).strip() or f"gene_{idx + 1}"
            for idx in kept_gene_indices
        ]
        kept_subset = subset[kept_gene_indices, :]

        if len(kept_gene_indices) <= max_display_columns:
            display_matrix = kept_subset.T.astype(np.float32, copy=True)
            column_labels = tuple(gene_labels)
            column_sizes = tuple(1 for _ in gene_labels)
            binned = False
        else:
            edges = np.linspace(0, len(kept_gene_indices), num=max_display_columns + 1, dtype=int)
            display_matrix = np.empty((len(samples), max_display_columns), dtype=np.float32)
            column_labels_list: List[str] = []
            column_sizes_list: List[int] = []
            out_col = 0
            for start, end in zip(edges[:-1], edges[1:]):
                if end <= start:
                    continue
                block = kept_subset[start:end, :]
                display_matrix[:, out_col] = np.asarray(block.mean(axis=0), dtype=np.float32)
                first_label = gene_labels[start]
                last_label = gene_labels[end - 1]
                column_labels_list.append(first_label if end - start == 1 else f"{first_label} .. {last_label}")
                column_sizes_list.append(end - start)
                out_col += 1
            if out_col <= 0:
                return []
            display_matrix = display_matrix[:, :out_col]
            column_labels = tuple(column_labels_list)
            column_sizes = tuple(column_sizes_list)
            binned = True

        return [
            {
                "mode": "heatmap_matrix",
                "title": f"Panaroo genes ({len(kept_gene_indices)})",
                "sample_order": tuple(samples),
                "matrix": display_matrix,
                "column_labels": column_labels,
                "column_sizes": column_sizes,
                "present_color": QtGui.QColor("#2563eb"),
                "absent_color": QtGui.QColor("#e2e8f0"),
                "binned": binned,
            }
        ]

    def _categorical_metadata_track(self, column: str) -> Optional[Dict[str, object]]:
        resolved = _CategoricalMetadataStripPanel._resolve_column(list(self._metadata_df.columns), (column,))
        if not resolved or self._metadata_df.empty or "Sample" not in self._metadata_df.columns:
            return None
        missing_label = "Missing"
        values: Dict[str, str] = {}
        seen_samples: set[str] = set()
        for _, row in self._metadata_df.iterrows():
            sample = self._normalize_metadata_value(row.get("Sample"))
            if not sample or sample in seen_samples:
                continue
            seen_samples.add(sample)
            value = self._normalize_metadata_value(row.get(resolved)) or missing_label
            values[sample] = value
        if not values:
            return None
        unique_values = sorted(
            {value for value in values.values() if value != missing_label},
            key=lambda value: str(value).casefold(),
        )
        value_colors: Dict[str, QtGui.QColor] = {
            value: _high_contrast_color(idx)
            for idx, value in enumerate(unique_values)
        }
        value_colors[missing_label] = QtGui.QColor("#e2e8f0")
        return {
            "title": column,
            "values": values,
            "value_colors": value_colors,
            "missing_label": missing_label,
        }

    def _bottom_legend_sections(self) -> List[Dict[str, object]]:
        sections: List[Dict[str, object]] = []
        if self._show_clade_legend and self._lineage_color_map:
            clade_items = [
                (f"Clade {value}", QtGui.QColor(color))
                for value, color in sorted(self._lineage_color_map.items(), key=lambda item: _natural_sort_key(item[0]))
            ]
            if clade_items:
                sections.append({"title": "Clades", "items": clade_items})
        if self._show_cluster_legend and self._cluster_by_sample:
            cluster_colors = _fastgear_cluster_colors(self._cluster_by_sample)
            cluster_items = [
                (f"Cluster {value}", QtGui.QColor(color))
                for value, color in sorted(cluster_colors.items(), key=lambda item: _natural_sort_key(item[0]))
            ]
            if cluster_items:
                sections.append({"title": "Clusters", "items": cluster_items})
        for column in self._side_metadata_columns:
            if column == _NEO_TREE_RESISTANCE_TOKEN:
                sections.append(
                    {
                        "title": "AMR",
                        "items": [
                            ("Present", QtGui.QColor(self._amr_present_color)),
                            ("Absent", QtGui.QColor("#e2e8f0")),
                        ],
                    }
                )
                continue
            if column == _NEO_TREE_MGE_TOKEN:
                sections.append(
                    {
                        "title": "MGEs",
                        "items": [
                            ("Present", QtGui.QColor(self._mge_present_color)),
                            ("Absent", QtGui.QColor("#e2e8f0")),
                        ],
                    }
                )
                continue
            if column == _NEO_TREE_VIRULENCE_TOKEN:
                sections.append(
                    {
                        "title": "Virulence",
                        "items": [
                            ("Present", QtGui.QColor(self._virulence_present_color)),
                            ("Absent", QtGui.QColor("#e2e8f0")),
                        ],
                    }
                )
                continue
            if column == _NEO_TREE_PANAROO_TOKEN:
                sections.append(
                    {
                        "title": "Panaroo genes",
                        "items": [
                            ("Present", QtGui.QColor("#2563eb")),
                            ("Absent", QtGui.QColor("#e2e8f0")),
                        ],
                    }
                )
                continue
            if column == _NEO_TREE_PLASMID_TOKEN:
                sections.append(
                    {
                        "title": "Plasmids",
                        "items": [
                            ("Present", QtGui.QColor("#16a34a")),
                            ("Absent", QtGui.QColor("#e2e8f0")),
                        ],
                    }
                )
                continue
            track = self._categorical_metadata_track(column)
            if track is None:
                continue
            missing_label = str(track.get("missing_label") or "Missing").strip() or "Missing"
            color_map = dict(track.get("value_colors", {}))
            items = [
                (value, QtGui.QColor(color_map[value]))
                for value in sorted(color_map.keys(), key=lambda value: (value == missing_label, str(value).casefold()))
                if value in color_map
            ]
            if items:
                sections.append({"title": column, "items": items})
        return sections

    def _rebuild_side_metadata_legend_views(self) -> None:
        if self._legend_layout is None or self._legend_bar is None:
            return
        for view in self.side_metadata_views:
            self._legend_layout.removeWidget(view)
            view.setParent(None)
            view.deleteLater()
        self.side_metadata_views = []
        self._side_metadata_view_by_column = {}
        if self.resistance_legend_view is not None:
            self._legend_layout.removeWidget(self.resistance_legend_view)
            self.resistance_legend_view.setParent(None)
            self.resistance_legend_view.deleteLater()
            self.resistance_legend_view = None
        if self.mge_legend_view is not None:
            self._legend_layout.removeWidget(self.mge_legend_view)
            self.mge_legend_view.setParent(None)
            self.mge_legend_view.deleteLater()
            self.mge_legend_view = None
        if self.virulence_legend_view is not None:
            self._legend_layout.removeWidget(self.virulence_legend_view)
            self.virulence_legend_view.setParent(None)
            self.virulence_legend_view.deleteLater()
            self.virulence_legend_view = None

        insert_index = self._legend_layout.count()
        if insert_index > 0:
            last_item = self._legend_layout.itemAt(insert_index - 1)
            if last_item is not None and last_item.spacerItem() is not None:
                insert_index -= 1
        selected = sorted(self._sample_filter) if self._sample_filter is not None else None
        for column in self._side_metadata_columns:
            if column == _NEO_TREE_RESISTANCE_TOKEN:
                resistance_legend = _StaticLegendPanel(
                    "AMR",
                    [
                        ("Present", QtGui.QColor(self._amr_present_color)),
                        ("Absent", QtGui.QColor("#e2e8f0")),
                    ],
                    parent=self._legend_bar,
                )
                self._legend_layout.insertWidget(
                    insert_index,
                    resistance_legend,
                    0,
                    QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignTop,
                )
                insert_index += 1
                self.resistance_legend_view = resistance_legend
                continue
            if column == _NEO_TREE_MGE_TOKEN:
                mge_legend = _StaticLegendPanel(
                    "MGEs",
                    [
                        ("Present", QtGui.QColor(self._mge_present_color)),
                        ("Absent", QtGui.QColor("#e2e8f0")),
                    ],
                    parent=self._legend_bar,
                )
                self._legend_layout.insertWidget(
                    insert_index,
                    mge_legend,
                    0,
                    QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignTop,
                )
                insert_index += 1
                self.mge_legend_view = mge_legend
                continue
            if column == _NEO_TREE_VIRULENCE_TOKEN:
                virulence_legend = _StaticLegendPanel(
                    "Virulence",
                    [
                        ("Present", QtGui.QColor(self._virulence_present_color)),
                        ("Absent", QtGui.QColor("#e2e8f0")),
                    ],
                    parent=self._legend_bar,
                )
                self._legend_layout.insertWidget(
                    insert_index,
                    virulence_legend,
                    0,
                    QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignTop,
                )
                insert_index += 1
                self.virulence_legend_view = virulence_legend
                continue
            if column in {_NEO_TREE_PANAROO_TOKEN, _NEO_TREE_PLASMID_TOKEN}:
                # These render as compact matrix tracks. The bottom legend
                # already carries their present/absent colors.
                continue
            side_metadata = _CategoricalMetadataStripPanel(
                self._metadata_df,
                title=column,
                column_candidates=(column,),
                show_canvas=False,
                legend_scrollable=False,
                legend_columns=2,
                parent=self._legend_bar,
            )
            side_metadata.set_tree_view(self.tree)
            side_metadata.set_sample_filter(selected)
            self._legend_layout.insertWidget(
                insert_index,
                side_metadata,
                0,
                QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignTop,
            )
            insert_index += 1
            self.side_metadata_views.append(side_metadata)
            self._side_metadata_view_by_column[column] = side_metadata

    def set_side_metadata_columns(self, columns: Optional[Sequence[str]]) -> None:
        normalized = self._normalize_side_metadata_columns(columns)
        if normalized == self._side_metadata_columns:
            return
        self._side_metadata_columns = normalized
        self._rebuild_side_metadata_legend_views()
        self._apply_side_metadata_state(refresh_views=True)

    def _ensure_detail_view_ready(
        self,
        stack: Optional[QtWidgets.QStackedWidget],
        index: int,
    ) -> Optional[QtWidgets.QWidget]:
        if stack is None or not (0 <= index < stack.count()):
            return None
        widget = stack.widget(index)
        if widget is None:
            return None
        lazy_view_type = str(widget.property("lazy_view_type") or "").strip().lower()
        if not lazy_view_type:
            if isinstance(widget, CgmlstForceGraphView):
                widget.ensure_fresh()
            elif isinstance(widget, PanHeatmapWidget):
                widget.refresh()
            return widget

        placeholder = widget
        selected = sorted(self._sample_filter) if self._sample_filter is not None else None
        if lazy_view_type == "tree":
            new_view: QtWidgets.QWidget = self._create_tree_view(stack)
            if self.tree is None:
                self.tree = new_view  # type: ignore[assignment]
        elif lazy_view_type == "presence":
            heatmap = PanHeatmapWidget(self._output_dir, parent=stack, auto_refresh=False)
            self._heatmap_views.append(heatmap)
            heatmap.sampleHovered.connect(self._on_heatmap_sample_hover)
            self.heatmap = self.heatmap or heatmap
            tree_order = self.tree.leaf_order() if self.tree is not None else []
            heatmap.set_sample_order(list(reversed(tree_order)) if tree_order else None)
            heatmap.set_filters(selected, None, None)
            heatmap.refresh()
            new_view = heatmap
        elif lazy_view_type == "geo":
            geo_view = GeoMapLibreView(parent=stack)
            self._geo_views.append(geo_view)
            self.geo_view = self.geo_view or geo_view
            geo_view.update_data(self._last_geo_points, self._last_geo_countries, self._last_geo_max_count)
            new_view = geo_view
        elif lazy_view_type == "metadata":
            metadata_view = _MetadataDonutPanel(self._metadata_df, parent=stack)
            self._metadata_views.append(metadata_view)
            self.metadata_view = self.metadata_view or metadata_view
            metadata_view.set_sample_filter(selected)
            new_view = metadata_view
        elif lazy_view_type == "force":
            force_view = CgmlstForceGraphView(
                self._output_dir,
                reference_sample=self._reference_sample,
                parent=stack,
            )
            self._force_views.append(force_view)
            force_view.sampleHovered.connect(self._on_force_sample_hover)
            force_view.set_fastgear_assignments(
                self._lineage_by_sample,
                self._lineage_color_map,
                self._cluster_by_sample,
            )
            force_view.set_sample_hover_metadata(self._sample_hover_metadata)
            force_view.set_reference_sample(self._reference_sample, refresh=False)
            force_view.set_sample_filter(selected, refresh=False)
            force_view.refresh()
            if self.force_view is None:
                self.force_view = force_view
            new_view = force_view
        elif lazy_view_type == "plasmidome":
            plasmidome = PlasmidomePanel(
                self._ann_df,
                self._output_dir,
                self._plasmid_cluster_index,
                parent=stack,
            )
            self._plasmidome_views.append(plasmidome)
            self.plasmidome_view = self.plasmidome_view or plasmidome
            plasmidome.samplesHovered.connect(self._on_plasmidome_samples_hover)
            plasmidome.sampleSubsetRequested.connect(self._on_plasmid_clicked)
            plasmidome.contextRequested.connect(self._on_plasmid_context_requested)
            plasmidome.set_sample_filter(selected)
            new_view = plasmidome
        else:
            return widget

        stack.insertWidget(index, new_view)
        stack.removeWidget(placeholder)
        placeholder.deleteLater()
        stack.setCurrentWidget(new_view)
        if isinstance(new_view, RaxmlTreeView):
            self._apply_side_metadata_state(refresh_views=True)
            self._apply_tree_source()
            new_view.set_active_samples(selected)
        return new_view

    def _on_detail_view_changed(self, stack: Optional[QtWidgets.QStackedWidget], index: int) -> None:
        if stack is None:
            return
        if 0 <= index < stack.count():
            stack.setCurrentIndex(index)
            self._ensure_detail_view_ready(stack, index)

    def _apply_side_metadata_state(self, *, refresh_views: bool = False) -> None:
        if self.tree is None:
            return
        selected = sorted(self._sample_filter) if self._sample_filter is not None else None
        tracks: List[Dict[str, object]] = []
        for column in self._side_metadata_columns:
            if column == _NEO_TREE_RESISTANCE_TOKEN:
                tracks.extend(self._resistance_tracks())
                continue
            if column == _NEO_TREE_MGE_TOKEN:
                tracks.extend(self._mge_tracks())
                continue
            if column == _NEO_TREE_VIRULENCE_TOKEN:
                tracks.extend(self._virulence_tracks())
                continue
            if column == _NEO_TREE_PANAROO_TOKEN:
                tracks.extend(self._panaroo_tracks())
                continue
            if column == _NEO_TREE_PLASMID_TOKEN:
                tracks.extend(self._plasmid_tracks())
                continue
            track = self._categorical_metadata_track(column)
            if track is not None:
                tracks.append(track)
        for tree_view in self._all_tree_views():
            tree_view.set_leaf_side_metadata_tracks(tracks or None)
            tree_view.set_bottom_legend_sections(self._bottom_legend_sections())

    def set_sample_filter(self, samples: Optional[Sequence[str]]) -> None:
        next_filter = {str(s).strip() for s in samples if str(s).strip()} if samples is not None else None
        if next_filter == self._sample_filter:
            return
        self._sample_filter = next_filter
        selected = sorted(self._sample_filter) if self._sample_filter is not None else None
        if self._all_tree_views():
            self._apply_tree_source()
            for tree_view in self._all_tree_views():
                tree_view.set_active_samples(selected)
        for heatmap_view in self._all_heatmaps():
            heatmap_view.set_filters(selected, None, None)
        for metadata_view in self._all_metadata_views():
            metadata_view.set_sample_filter(selected)
        for force_view in self._all_force_views():
            force_view.set_sample_filter(selected, refresh=force_view.isVisible())
        for plasmidome_view in self._all_plasmidome_views():
            plasmidome_view.set_sample_filter(selected)
        if self.clade_legend_view is not None:
            self.clade_legend_view.set_sample_filter(selected)
        if self.cluster_legend_view is not None:
            self.cluster_legend_view.set_sample_filter(selected)
        self._apply_side_metadata_state(refresh_views=False)

    def update_geo_data(self, points: List[dict], countries: Dict[str, int], max_count: int) -> None:
        self._last_geo_points = list(points)
        self._last_geo_countries = dict(countries)
        self._last_geo_max_count = int(max_count)
        for geo_view in self.all_geo_views():
            geo_view.update_data(points, countries, max_count)

    def _apply_detail_highlight(
        self,
        sample: Optional[str] = None,
        samples: Optional[Sequence[str]] = None,
        *,
        visible_only: bool = True,
    ) -> None:
        normalized = [str(s).strip() for s in (samples or []) if str(s).strip()]
        if not normalized and sample is not None and str(sample).strip():
            normalized = [str(sample).strip()]

        heatmap_views = self._visible_views(self._all_heatmaps()) if visible_only else self._all_heatmaps()
        for heatmap_view in heatmap_views:
            if len(normalized) >= 2:
                heatmap_view.highlight_sample(None)
                heatmap_view.highlight_samples(normalized)
            elif normalized:
                heatmap_view.highlight_samples(None)
                heatmap_view.highlight_sample(normalized[0])
            else:
                heatmap_view.highlight_samples(None)
                heatmap_view.highlight_sample(None)
        geo_views = self._visible_views(self.all_geo_views()) if visible_only else self.all_geo_views()
        for geo_view in geo_views:
            if len(normalized) >= 2:
                geo_view.highlight_samples(normalized)
            elif normalized:
                geo_view.highlight_sample(normalized[0])
            else:
                geo_view.highlight_samples(None)

    def _schedule_hover_highlight(
        self,
        samples: Optional[Sequence[str]],
        emit_sample: Optional[str],
        *,
        source_tree: Optional[RaxmlTreeView] = None,
    ) -> None:
        normalized = tuple(str(s).strip() for s in (samples or []) if str(s).strip())
        self._pending_hover_samples = normalized
        self._pending_hover_emit_sample = str(emit_sample).strip() if emit_sample else None
        self._pending_hover_source_tree = source_tree
        self._hover_update_timer.start()

    def _flush_hover_highlight(self) -> None:
        samples = self._pending_hover_samples
        source_tree = self._pending_hover_source_tree
        if samples != self._last_applied_hover_samples:
            if len(samples) == 1:
                self.highlight_sample(samples[0], visible_only=True, exclude_tree=source_tree)
            elif len(samples) > 1:
                self.highlight_samples(samples, visible_only=True, exclude_tree=source_tree)
            else:
                self.highlight_sample(None, visible_only=True, exclude_tree=source_tree)
            self._last_applied_hover_samples = samples
        self.sampleHovered.emit(self._pending_hover_emit_sample or None)

    def _on_clade_clicked(self, payload: object) -> None:
        if not isinstance(payload, dict):
            return
        self.cladeFilterRequested.emit(dict(payload))

    def _on_cluster_clicked(self, payload: object) -> None:
        if not isinstance(payload, dict):
            return
        self.clusterFilterRequested.emit(dict(payload))

    def _on_plasmid_clicked(self, payload: object) -> None:
        if not isinstance(payload, dict):
            return
        self.plasmidFilterRequested.emit(dict(payload))

    def _on_plasmid_context_requested(self, payload: object) -> None:
        if not isinstance(payload, dict):
            return
        self.plasmidContextRequested.emit(dict(payload))

    def _on_plasmidome_samples_hover(self, samples: object) -> None:
        if not isinstance(samples, (list, tuple)):
            self._schedule_hover_highlight(None, None)
            return
        normalized = [str(sample).strip() for sample in samples if str(sample).strip()]
        if len(normalized) == 1:
            self._schedule_hover_highlight(normalized, normalized[0])
        elif len(normalized) > 1:
            self._schedule_hover_highlight(normalized, None)
        else:
            self._schedule_hover_highlight(None, None)

    def _on_tree_leaf_double_clicked(self, sample: object) -> None:
        value = str(sample).strip() if sample is not None else ""
        if value:
            self.sampleOpenRequested.emit(value)

    def _export_tree_figure(self) -> None:
        if self.tree is None:
            return

        default_dir = self._output_dir if self._output_dir.exists() else Path.home()
        default_path = default_dir / "neotree_export.pdf"
        selected_path, selected_filter = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Export neoTree figure",
            str(default_path),
            "PDF Document (*.pdf);;PNG Image (*.png);;JPEG Image (*.jpg *.jpeg)",
        )
        if not selected_path:
            return

        export_path = Path(selected_path)
        if not export_path.suffix:
            filter_lower = str(selected_filter or "").lower()
            if "png" in filter_lower:
                export_path = export_path.with_suffix(".png")
            elif "jpeg" in filter_lower or "jpg" in filter_lower:
                export_path = export_path.with_suffix(".jpg")
            else:
                export_path = export_path.with_suffix(".pdf")

        try:
            self.tree.export_figure(export_path)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(
                self,
                "Export failed",
                f"Could not export the neoTree figure:\n{exc}",
            )
            return

        QtWidgets.QMessageBox.information(
            self,
            "Export complete",
            f"Saved neoTree figure to:\n{export_path}",
        )

    def _load_fastgear_clades(
        self,
    ) -> Tuple[Dict[str, str], Dict[str, QtGui.QColor], Dict[str, str], Optional[Path]]:
        path = _find_fastgear_sample_clades(self._output_dir)
        if path is None:
            self._loaded_sample_clades_sig = None
            self._lineage_by_sample = {}
            self._lineage_color_map = {}
            self._cluster_by_sample = {}
            return {}, {}, {}, None

        sig = _file_signature(path)
        if self._loaded_sample_clades_sig != sig:
            try:
                lineage_by_sample = _load_fastgear_lineages(path)
                cluster_by_sample = _load_fastgear_clusters(path)
            except Exception:
                lineage_by_sample = {}
                cluster_by_sample = {}
            self._lineage_by_sample = lineage_by_sample
            self._lineage_color_map = _fastgear_lineage_colors(lineage_by_sample)
            self._cluster_by_sample = cluster_by_sample
            self._loaded_sample_clades_sig = sig
            _log_view_event(
                f"loaded fastGEAR assignments: {path} "
                f"({len(set(lineage_by_sample.values()))} clades, "
                f"{len(set(cluster_by_sample.values()))} clusters)"
            )
        return (
            dict(self._lineage_by_sample),
            dict(self._lineage_color_map),
            dict(self._cluster_by_sample),
            path,
        )

    def _find_panaroo_tree(self) -> Optional[Path]:
        out = _find_panaroo_out(self._output_dir)
        if out is None:
            return None
        patterns = [
            "*.nwk",
            "*.newick",
            "*tree*.nwk",
            "*tree*.newick",
            "*bestTree*",
        ]
        matches: List[Path] = []
        for pattern in patterns:
            matches.extend([p for p in out.rglob(pattern) if p.is_file()])
        if not matches:
            return None
        matches.sort(key=lambda p: p.stat().st_mtime, reverse=True)
        return matches[0]

    def _alignment_selected_samples(self, alignment_labels: Sequence[str]) -> List[str]:
        if self._sample_filter is None:
            return list(alignment_labels)
        selected_set = self._sample_filter
        return [label for label in alignment_labels if label in selected_set]

    def _build_subset_tree_from_alignment(
        self,
        alignment_path: Path,
        selected_samples: Sequence[str],
    ) -> Optional[str]:
        selected = tuple(str(sample).strip() for sample in selected_samples if str(sample).strip())
        if not selected:
            return None
        newick = _build_subset_tree_newick_cached(*_file_signature(alignment_path), selected)
        return newick or None

    def _apply_tree_source(self) -> None:
        tree_views = self._all_tree_views()
        if not tree_views:
            return

        raxml_tree = self._fallback_tree if self._fallback_tree is not None and self._fallback_tree.is_file() else None
        panaroo_tree = self._find_panaroo_tree()
        canonical_tree: Optional[Path] = None
        source_note = ""
        if raxml_tree is not None:
            canonical_tree = raxml_tree
            if "snippy_raxml" in raxml_tree.name:
                source_note = "RAxML SNP"
            elif raxml_tree.name == "core_gene_tree.nwk":
                source_note = "RAxML core gene"
            else:
                source_note = "RAxML"
        elif panaroo_tree is not None:
            canonical_tree = panaroo_tree
            source_note = "Panaroo"

        alignment_path = _find_panaroo_alignment(self._output_dir) if self._recompute_tree_from_sample_filter else None
        alignment_labels: Tuple[str, ...] = tuple()
        alignment_sig: Optional[Tuple[str, int, int]] = None
        if alignment_path is not None and alignment_path.is_file():
            try:
                alignment_sig = _file_signature(alignment_path)
                alignment_labels, _alignment_matrix = _load_core_alignment_cached(*alignment_sig)
                _log_view_event(f"loaded Panaroo core alignment: {alignment_path}")
            except Exception as exc:
                for tree_view in tree_views:
                    tree_view.clear_tree(f"Failed to load Panaroo alignment: {exc}")
                self._loaded_tree_key = None
                if self._show_info_label and self._info is not None:
                    self._info.setText(f"Failed to load Panaroo alignment: {exc}")
                    self._info.setVisible(True)
                return

        selected_alignment_samples = self._alignment_selected_samples(alignment_labels) if alignment_labels else []
        full_alignment_visible = (
            alignment_labels
            and len(selected_alignment_samples) == len(alignment_labels)
        )
        should_recompute = bool(
            self._recompute_tree_from_sample_filter
            and alignment_sig is not None
            and (
                canonical_tree is None
                or (self._sample_filter is not None and not full_alignment_visible)
            )
        )

        if should_recompute and alignment_path is not None:
            if not selected_alignment_samples:
                for tree_view in tree_views:
                    tree_view.clear_tree("No filtered samples are present in the Panaroo core alignment.")
                self._loaded_tree_key = None
                if self._show_info_label and self._info is not None:
                    self._info.setText("No filtered samples are present in the Panaroo core alignment.")
                    self._info.setVisible(True)
                return
            try:
                newick = self._build_subset_tree_from_alignment(alignment_path, selected_alignment_samples)
            except Exception as exc:
                for tree_view in tree_views:
                    tree_view.clear_tree(f"Failed to recompute subset tree: {exc}")
                self._loaded_tree_key = None
                if self._show_info_label and self._info is not None:
                    self._info.setText(f"Failed to recompute subset tree: {exc}")
                    self._info.setVisible(True)
                return
            if not newick:
                for tree_view in tree_views:
                    tree_view.clear_tree("Subset tree could not be computed from the Panaroo core alignment.")
                self._loaded_tree_key = None
                if self._show_info_label and self._info is not None:
                    self._info.setText("Subset tree could not be computed from the Panaroo core alignment.")
                    self._info.setVisible(True)
                return
            tree_key: Tuple[object, ...] = ("subset", *(alignment_sig or ("", 0, 0)), tuple(selected_alignment_samples))
            if self._loaded_tree_key != tree_key:
                for tree_view in tree_views:
                    tree_view.set_tree_text(newick)
                self._loaded_tree_key = tree_key
                _log_view_event(
                    f"loaded recomputed subset tree from {alignment_path} "
                    f"({len(selected_alignment_samples)} samples)"
                )
            if self._show_info_label and self._info is not None:
                self._info.clear()
                self._info.setVisible(False)
            return

        if canonical_tree is not None:
            tree_sig = _file_signature(canonical_tree)
            tree_key = ("file", *tree_sig)
            if self._loaded_tree_key != tree_key:
                for tree_view in tree_views:
                    tree_view.set_tree_file(canonical_tree)
                self._loaded_tree_key = tree_key
                _log_view_event(f"loaded {source_note} tree: {canonical_tree}")
            if self._show_info_label and self._info is not None:
                self._info.clear()
                self._info.setVisible(False)
            return

        for tree_view in tree_views:
            tree_view.clear_tree("No tree file found for this tab.")
        self._loaded_tree_key = None
        if self._show_info_label and self._info is not None:
            self._info.setText("No Panaroo or fallback tree file found.")
            self._info.setVisible(True)

    def _on_tree_leaf_hover(self, sample: object) -> None:
        name = str(sample).strip() if sample is not None else ""
        source_tree = self.sender()
        self._schedule_hover_highlight(
            [name] if name else None,
            name or None,
            source_tree=source_tree if isinstance(source_tree, RaxmlTreeView) else None,
        )

    def _on_tree_subtree_hover(self, samples: object) -> None:
        source_tree = self.sender()
        tree_source = source_tree if isinstance(source_tree, RaxmlTreeView) else None
        if not isinstance(samples, (list, tuple)):
            self._schedule_hover_highlight(None, None, source_tree=tree_source)
            return
        names = [str(s).strip() for s in samples if str(s).strip()]
        if len(names) >= 2:
            self._schedule_hover_highlight(names, None, source_tree=tree_source)
        elif len(names) == 1:
            self._schedule_hover_highlight(names, names[0], source_tree=tree_source)
        else:
            self._schedule_hover_highlight(None, None, source_tree=tree_source)

    def _on_heatmap_sample_hover(self, sample: object) -> None:
        name = str(sample).strip() if sample is not None else ""
        self._schedule_hover_highlight([name] if name else None, name or None)

    def _on_force_sample_hover(self, sample: object) -> None:
        name = str(sample).strip() if sample is not None else ""
        self._schedule_hover_highlight([name] if name else None, name or None)

    def highlight_sample(
        self,
        sample: Optional[str],
        *,
        visible_only: bool = True,
        exclude_tree: Optional[RaxmlTreeView] = None,
    ) -> None:
        value = str(sample).strip() if sample is not None else ""
        tree_views = self._visible_views(self._all_tree_views()) if visible_only else self._all_tree_views()
        for tree_view in tree_views:
            if exclude_tree is not None and tree_view is exclude_tree:
                continue
            tree_view.highlight_sample(value or None)
        self._apply_detail_highlight(value or None, visible_only=visible_only)
        force_views = self._visible_views(self._all_force_views()) if visible_only else self._all_force_views()
        for force_view in force_views:
            force_view.highlight_sample(value or None)

    def highlight_samples(
        self,
        samples: Optional[Sequence[str]],
        *,
        visible_only: bool = True,
        exclude_tree: Optional[RaxmlTreeView] = None,
    ) -> None:
        tree_views = self._visible_views(self._all_tree_views()) if visible_only else self._all_tree_views()
        for tree_view in tree_views:
            if exclude_tree is not None and tree_view is exclude_tree:
                continue
            tree_view.highlight_samples(samples)
        self._apply_detail_highlight(samples=samples, visible_only=visible_only)
        force_views = self._visible_views(self._all_force_views()) if visible_only else self._all_force_views()
        for force_view in force_views:
            force_view.highlight_samples(samples)

    def refresh(self) -> None:
        tree_views = self._all_tree_views()
        if not tree_views:
            return
        lineage_by_sample: Dict[str, str] = {}
        lineage_colors: Dict[str, QtGui.QColor] = {}
        cluster_by_sample: Dict[str, str] = {}
        cluster_colors: Dict[str, QtGui.QColor] = {}
        if self._color_tree_by_fastgear_clade or self._color_tree_by_fastgear_cluster or self._annotate_fastgear_clusters:
            lineage_by_sample, lineage_colors, cluster_by_sample, _ = self._load_fastgear_clades()
            cluster_colors = _fastgear_cluster_colors(cluster_by_sample)
        for force_view in self._all_force_views():
            force_view.set_fastgear_assignments(lineage_by_sample, lineage_colors, cluster_by_sample)
            force_view.set_sample_hover_metadata(self._sample_hover_metadata)
            force_view.set_reference_sample(self._reference_sample, refresh=False)
        for tree_view in tree_views:
            tree_view.set_leaf_hover_metadata(self._sample_hover_metadata)
            if self._color_tree_by_fastgear_cluster:
                tree_view.set_leaf_display_groups(cluster_by_sample, group_colors=cluster_colors)
            else:
                tree_view.set_leaf_display_groups(None, group_colors=None)
            if self._color_tree_by_fastgear_clade:
                tree_view.set_leaf_groups(lineage_by_sample, group_colors=lineage_colors)
            else:
                tree_view.set_leaf_groups(None, group_colors=None)
            if self._annotate_fastgear_clusters:
                tree_view.set_leaf_clusters(cluster_by_sample)
            else:
                tree_view.set_leaf_clusters(None)
        if self.clade_legend_view is not None:
            self.clade_legend_view.set_assignments(lineage_by_sample, lineage_colors, show_missing=False)
        if self.cluster_legend_view is not None:
            self.cluster_legend_view.set_assignments(cluster_by_sample, cluster_colors, show_missing=False)
        self._apply_side_metadata_state(refresh_views=True)
        self._apply_tree_source()

        tree_order = self.tree.leaf_order() if self.tree is not None else []
        for heatmap_view in self._all_heatmaps():
            if tree_order:
                heatmap_view.set_sample_order(list(reversed(tree_order)))
            else:
                heatmap_view.set_sample_order(None)
        selected = sorted(self._sample_filter) if self._sample_filter is not None else None
        for tree_view in tree_views:
            tree_view.set_active_samples(selected)
        for heatmap_view in self._all_heatmaps():
            heatmap_view.set_filters(selected, None, None)
        for metadata_view in self._all_metadata_views():
            metadata_view.set_sample_filter(selected)
        for force_view in self._all_force_views():
            force_view.set_sample_filter(selected, refresh=False)
            if force_view.isVisible():
                force_view.ensure_fresh()
        if self.clade_legend_view is not None:
            self.clade_legend_view.set_sample_filter(selected)
        if self.cluster_legend_view is not None:
            self.cluster_legend_view.set_sample_filter(selected)
        self._apply_side_metadata_state(refresh_views=True)
