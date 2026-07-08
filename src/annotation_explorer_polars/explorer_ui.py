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

from .annotation_detail_dialog import AnnotationDetailDialog
from .annotation_table import AnnotationTableModel, DetailsDelegate, TagsDelegate
from .constants import ANNOTATION_COLUMNS_TO_SHOW, MISSING_FILTER_TOKEN
from .geo_view import GeoMapView
from .igv_view import GeneSegment, IgvTrackView
from .metadata_filter_dialog import MetadataValueDialog

# Shared frame + column-header look for the two document tables (the left sample
# tree and the right annotation table) so they read as one consistent,
# "sample-table" style. Crucially this carries NO `::item` rule: a QTableView
# with an `::item` stylesheet stops honouring the model's BackgroundRole, which
# would wipe the annotation table's group / separator / hit-row shading. Vertical
# column separators on the table come from its grid (gridline-color) instead.
# No explicit header font-size, so the titles inherit the app font rather than
# the old cramped 11px.
_FRAME_HEADER_QSS = (
    "QTreeView, QTableView {"
    "  border: 1px solid #cdd6e0;"
    "  background: #ffffff;"
    "  alternate-background-color: #f6f8fb;"
    "  gridline-color: #e7ebf1;"
    "  outline: none;"
    "  color: #23303f;"
    "}"
    "QHeaderView::section {"
    "  border: none;"
    "  border-right: 1px solid #dbe2ea;"
    "  border-bottom: 1px solid #cdd6e0;"
    "  background: #eef2f7;"
    "  padding: 5px 8px;"
    "  color: #55677a;"
    "  font-weight: 600;"
    "}"
    "QHeaderView::section:hover { background: #e6ecf3; }"
    "QTreeView::branch { background: transparent; }"
)

# The tree additionally styles individual items (vertical separators, hover,
# selection). QTreeWidgetItem.setBackground still paints under these rules, so
# the tree keeps its heat cells; this block is the tree's original look minus the
# old 11px header font-size.
SHARED_ITEM_VIEW_QSS = _FRAME_HEADER_QSS + (
    "QTreeView::item {"
    "  border-right: 1px solid #e7ebf1;"
    "  padding: 3px 4px;"
    "}"
    "QTreeView::item:hover { background: #eef4fc; }"
    "QTreeView::item:selected {"
    "  background: #d6e6fb;"
    "  color: #16324f;"
    "  border-right: 1px solid #c3d8f2;"
    "}"
    "QTreeView::item:selected:!active {"
    "  background: #dde9f8;"
    "  color: #16324f;"
    "}"
)


class ExplorerUIMixin:
    def _build_ui(self) -> None:
        self.setWindowTitle("Annotation Explorer (PyQt6)")

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)

        main_layout = QtWidgets.QVBoxLayout(central)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(6)

        self.status = QtWidgets.QLabel()
        self.status.setVisible(False)
        main_layout.addWidget(
            self.status,
            0,
            QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter,
        )

        content_row = QtWidgets.QHBoxLayout()
        content_row.setContentsMargins(8, 8, 8, 8)
        content_row.setSpacing(8)
        main_layout.addLayout(content_row, stretch=1)

        filters_panel = QtWidgets.QWidget()
        self._filters_panel = filters_panel
        filters_layout = QtWidgets.QVBoxLayout(filters_panel)
        filters_layout.setContentsMargins(0, 0, 0, 0)
        filters_layout.setSpacing(8)

        meta_group = QtWidgets.QGroupBox("Sample")
        meta_layout = QtWidgets.QVBoxLayout(meta_group)
        self.sample_query = QtWidgets.QLineEdit()
        self.sample_query.setPlaceholderText("Search samples and metadata...")
        self.sample_query.setClearButtonEnabled(True)
        meta_layout.addWidget(self.sample_query)

        self.add_filter_btn = QtWidgets.QPushButton("Add filter")
        self.add_filter_btn.clicked.connect(self._show_add_filter_menu)
        self.add_filter_btn.setEnabled(bool(self.meta_columns))
        meta_layout.addWidget(self.add_filter_btn)

        self.filters_container = QtWidgets.QWidget()
        self.filters_layout = QtWidgets.QVBoxLayout(self.filters_container)
        self.filters_layout.setContentsMargins(0, 0, 0, 0)
        self.filters_layout.setSpacing(6)
        meta_layout.addWidget(self.filters_container)

        sort_group = QtWidgets.QGroupBox("Sorting")
        sort_layout = QtWidgets.QFormLayout(sort_group)
        self.sample_sort_column_box = QtWidgets.QComboBox()
        self.sample_sort_column_box.addItem("Sample", "__sample__")
        for col in sorted(self.meta_columns):
            self.sample_sort_column_box.addItem(col, col)
        self.sample_sort_order_box = QtWidgets.QComboBox()
        self.sample_sort_order_box.addItem("Ascending", "asc")
        self.sample_sort_order_box.addItem("Descending", "desc")
        sort_layout.addRow("Sort by", self.sample_sort_column_box)
        sort_layout.addRow("Order", self.sample_sort_order_box)
        meta_layout.addWidget(sort_group)

        self.hide_empty_contigs_chk = QtWidgets.QCheckBox("Hide contigs without hits")
        self.hide_empty_contigs_chk.setChecked(True)

        ann_group = QtWidgets.QGroupBox("Annotations")
        ann_layout = QtWidgets.QFormLayout(ann_group)

        self.table_search_hint = QtWidgets.QLabel("")
        self.table_search_hint.setVisible(False)
        self.table_search_hint.setStyleSheet(
            "color: #607183; padding: 2px 0 6px 0; font-size: 12px;"
        )
        ann_layout.addRow(self.table_search_hint)

        self.gene_query = QtWidgets.QLineEdit()
        ann_layout.addRow("Gene/annotation contains", self.gene_query)

        self.product_query = QtWidgets.QLineEdit()
        ann_layout.addRow("Product contains", self.product_query)

        self.molecule_filter_box = QtWidgets.QComboBox()
        self.molecule_filter_box.addItem("All", "")
        self.molecule_filter_box.addItem("Plasmid only", "plasmid")
        self.molecule_filter_box.addItem("Chromosome only", "chromosome")
        ann_layout.addRow("Molecule", self.molecule_filter_box)

        self.predictor_box = QtWidgets.QComboBox()
        self.predictor_box.setEditable(True)
        self.predictor_box.lineEdit().setReadOnly(True)
        self.predictor_box.lineEdit().setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft)
        self.predictor_box.setInsertPolicy(QtWidgets.QComboBox.InsertPolicy.NoInsert)
        self.predictor_box.setMaximumWidth(220)
        self.predictor_box.setMinimumContentsLength(12)
        self.predictor_box.view().pressed.connect(self._toggle_predictor_item)
        self._init_predictor_model()
        ann_layout.addRow("Predictor", self.predictor_box)

        self.tag_box = QtWidgets.QComboBox()
        self.tag_box.setEditable(True)
        self.tag_box.lineEdit().setReadOnly(True)
        self.tag_box.lineEdit().setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft)
        self.tag_box.setInsertPolicy(QtWidgets.QComboBox.InsertPolicy.NoInsert)
        self.tag_box.setMaximumWidth(220)
        self.tag_box.setMinimumContentsLength(12)
        self.tag_box.view().pressed.connect(self._toggle_tag_item)
        self._init_tag_model()
        ann_layout.addRow("Tags", self.tag_box)

        self.hide_hypo_chk = QtWidgets.QCheckBox()
        ann_layout.addRow("Hide hypothetical", self.hide_hypo_chk)

        self.hide_empty_contigs_chk.setText("")
        ann_layout.addRow("Hide contigs without hits", self.hide_empty_contigs_chk)

        view_label = QtWidgets.QLabel("View")
        self.view_mode_box = QtWidgets.QComboBox()
        self.view_mode_box.addItems(["Contig", "Context", "Global"])
        self.view_mode_box.currentTextChanged.connect(self._on_view_mode_changed)
        ann_layout.addRow(view_label, self.view_mode_box)

        table_group = QtWidgets.QGroupBox("Table")
        self.sample_table_columns_group = table_group
        table_layout = QtWidgets.QVBoxLayout(table_group)
        table_layout.setContentsMargins(8, 8, 8, 8)
        table_layout.setSpacing(6)

        self.table_columns_summary_label = QtWidgets.QLabel("")
        self.table_columns_summary_label.setStyleSheet("color: #607183; padding: 0 0 2px 0;")
        self.table_columns_summary_label.setWordWrap(True)
        table_layout.addWidget(self.table_columns_summary_label)

        self.table_columns_scroll = QtWidgets.QScrollArea()
        self.table_columns_scroll.setWidgetResizable(True)
        self.table_columns_scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.table_columns_scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.table_columns_scroll.setMinimumHeight(220)

        self.table_columns_container = QtWidgets.QWidget()
        self.table_columns_layout = QtWidgets.QVBoxLayout(self.table_columns_container)
        self.table_columns_layout.setContentsMargins(0, 0, 0, 0)
        self.table_columns_layout.setSpacing(8)
        self.table_columns_scroll.setWidget(self.table_columns_container)
        table_layout.addWidget(self.table_columns_scroll)
        self._init_table_column_controls()
        table_group.setVisible(False)

        self.neo_tree_traits_group = QtWidgets.QGroupBox("neoTree traits")
        neo_tree_layout = QtWidgets.QVBoxLayout(self.neo_tree_traits_group)
        neo_tree_layout.setContentsMargins(8, 8, 8, 8)
        neo_tree_layout.setSpacing(6)

        self.neo_tree_traits_summary_label = QtWidgets.QLabel("")
        self.neo_tree_traits_summary_label.setStyleSheet("color: #607183; padding: 0 0 2px 0;")
        self.neo_tree_traits_summary_label.setWordWrap(True)
        neo_tree_layout.addWidget(self.neo_tree_traits_summary_label)

        self.neo_tree_trait_colors_group = QtWidgets.QGroupBox("Track colors")
        neo_tree_trait_colors_layout = QtWidgets.QVBoxLayout(self.neo_tree_trait_colors_group)
        neo_tree_trait_colors_layout.setContentsMargins(8, 8, 8, 8)
        neo_tree_trait_colors_layout.setSpacing(6)

        self.neo_tree_amr_color_row = QtWidgets.QWidget()
        neo_tree_amr_color_row_layout = QtWidgets.QHBoxLayout(self.neo_tree_amr_color_row)
        neo_tree_amr_color_row_layout.setContentsMargins(0, 0, 0, 0)
        neo_tree_amr_color_row_layout.setSpacing(8)
        self.neo_tree_amr_color_label = QtWidgets.QLabel("AMR")
        self.neo_tree_amr_color_button = QtWidgets.QPushButton()
        self.neo_tree_amr_color_button.setMinimumWidth(110)
        neo_tree_amr_color_row_layout.addWidget(self.neo_tree_amr_color_label)
        neo_tree_amr_color_row_layout.addStretch(1)
        neo_tree_amr_color_row_layout.addWidget(self.neo_tree_amr_color_button)
        neo_tree_trait_colors_layout.addWidget(self.neo_tree_amr_color_row)

        self.neo_tree_mge_color_row = QtWidgets.QWidget()
        neo_tree_mge_color_row_layout = QtWidgets.QHBoxLayout(self.neo_tree_mge_color_row)
        neo_tree_mge_color_row_layout.setContentsMargins(0, 0, 0, 0)
        neo_tree_mge_color_row_layout.setSpacing(8)
        self.neo_tree_mge_color_label = QtWidgets.QLabel("MGEs")
        self.neo_tree_mge_color_button = QtWidgets.QPushButton()
        self.neo_tree_mge_color_button.setMinimumWidth(110)
        neo_tree_mge_color_row_layout.addWidget(self.neo_tree_mge_color_label)
        neo_tree_mge_color_row_layout.addStretch(1)
        neo_tree_mge_color_row_layout.addWidget(self.neo_tree_mge_color_button)
        neo_tree_trait_colors_layout.addWidget(self.neo_tree_mge_color_row)

        self.neo_tree_virulence_color_row = QtWidgets.QWidget()
        neo_tree_virulence_color_row_layout = QtWidgets.QHBoxLayout(self.neo_tree_virulence_color_row)
        neo_tree_virulence_color_row_layout.setContentsMargins(0, 0, 0, 0)
        neo_tree_virulence_color_row_layout.setSpacing(8)
        self.neo_tree_virulence_color_label = QtWidgets.QLabel("Virulence")
        self.neo_tree_virulence_color_button = QtWidgets.QPushButton()
        self.neo_tree_virulence_color_button.setMinimumWidth(110)
        neo_tree_virulence_color_row_layout.addWidget(self.neo_tree_virulence_color_label)
        neo_tree_virulence_color_row_layout.addStretch(1)
        neo_tree_virulence_color_row_layout.addWidget(self.neo_tree_virulence_color_button)
        neo_tree_trait_colors_layout.addWidget(self.neo_tree_virulence_color_row)
        neo_tree_layout.addWidget(self.neo_tree_trait_colors_group)

        self.neo_tree_traits_scroll = QtWidgets.QScrollArea()
        self.neo_tree_traits_scroll.setWidgetResizable(True)
        self.neo_tree_traits_scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.neo_tree_traits_scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.neo_tree_traits_scroll.setMinimumHeight(150)

        self.neo_tree_traits_container = QtWidgets.QWidget()
        self.neo_tree_traits_layout = QtWidgets.QVBoxLayout(self.neo_tree_traits_container)
        self.neo_tree_traits_layout.setContentsMargins(0, 0, 0, 0)
        self.neo_tree_traits_layout.setSpacing(8)
        self.neo_tree_traits_scroll.setWidget(self.neo_tree_traits_container)
        neo_tree_layout.addWidget(self.neo_tree_traits_scroll)
        self._init_neo_tree_trait_controls()

        filters_layout.addWidget(meta_group)
        filters_layout.addWidget(ann_group)
        filters_layout.addWidget(table_group)
        filters_layout.addStretch(1)
        content_row.addWidget(filters_panel, stretch=0)

        content_row.addWidget(self._build_gene_workspace(), stretch=1)

        self._set_bottom_tab_placeholder("Geo", "Geo map will load after the main view is ready.")
        self._set_bottom_tab_placeholder("neoPan", "neoPan will load when this tab is opened.")
        self._init_plot_menu()

        self._refresh_timer = QtCore.QTimer(self)
        self._refresh_timer.setSingleShot(True)
        self._refresh_timer.setInterval(180)
        self._refresh_timer.timeout.connect(self.refresh)
        for w in [self.sample_query, self.gene_query, self.product_query]:
            w.textChanged.connect(self._schedule_refresh)
        self.sample_sort_column_box.currentIndexChanged.connect(self._schedule_refresh)
        self.sample_sort_order_box.currentIndexChanged.connect(self._schedule_refresh)
        self.molecule_filter_box.currentIndexChanged.connect(self._schedule_refresh)
        self.hide_hypo_chk.stateChanged.connect(self._schedule_refresh)
        self.hide_empty_contigs_chk.stateChanged.connect(self._schedule_refresh)
        QtCore.QTimer.singleShot(0, self._update_side_panel_widths)
        QtCore.QTimer.singleShot(700, self._ensure_active_bottom_tab_initialized)

    def _build_gene_workspace(self) -> QtWidgets.QWidget:
        page = QtWidgets.QWidget()
        page_layout = QtWidgets.QVBoxLayout(page)
        page_layout.setContentsMargins(0, 0, 0, 0)
        page_layout.setSpacing(0)

        v_split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        v_split.setChildrenCollapsible(False)
        page_layout.addWidget(v_split, stretch=1)

        top_panel = QtWidgets.QWidget()
        top_layout = QtWidgets.QVBoxLayout(top_panel)
        top_layout.setContentsMargins(0, 0, 0, 0)
        top_layout.setSpacing(6)

        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        split.setChildrenCollapsible(False)
        self._top_splitter = split
        top_layout.addWidget(split, stretch=1)

        tree_panel = QtWidgets.QWidget()
        tree_layout = QtWidgets.QVBoxLayout(tree_panel)
        tree_layout.setContentsMargins(0, 0, 0, 0)
        tree_layout.setSpacing(4)

        self.tree = QtWidgets.QTreeWidget()
        # Geneious-Prime-style document table: real sortable columns, not a
        # single crammed column-0 widget.
        self._sample_tree_compact_view = False
        tree_headers = (
            ["Name", "Card", "Origin", "n_total", "ATB"]
            + self.meta_columns_for_tree
            + ["db_total_hits"]
            + self.annotation_db_columns
            + ["ST", "Quality", "Complete", "Contam."]
        )
        self.tree.setHeaderLabels(tree_headers)
        self.tree.setHeaderHidden(False)
        # The four fixed QC columns are appended after the existing (optional)
        # columns so the db-hit heat-cell indices below are unchanged.
        _qc_base = 5 + len(self.meta_columns_for_tree) + 1 + len(self.annotation_db_columns)
        self._col_st = _qc_base
        self._col_quality = _qc_base + 1
        self._col_completeness = _qc_base + 2
        self._col_contamination = _qc_base + 3
        self._fixed_tree_columns = (
            0, 2, self._col_st, self._col_quality, self._col_completeness, self._col_contamination
        )
        self.tree.setColumnWidth(2, 122)
        self.tree.setColumnWidth(3, 68)
        self.tree.setColumnWidth(4, 56)
        self.tree.setColumnWidth(self._col_st, 62)
        self.tree.setColumnWidth(self._col_quality, 82)
        self.tree.setColumnWidth(self._col_completeness, 96)
        self.tree.setColumnWidth(self._col_contamination, 88)
        self.tree.setAlternatingRowColors(True)
        self.tree.setUniformRowHeights(False)
        self.tree.setIndentation(16)
        self.tree.setIconSize(QtCore.QSize(16, 16))
        self.tree.setAnimated(True)
        self.tree.setAllColumnsShowFocus(True)
        self.tree.setRootIsDecorated(True)
        self.tree.setExpandsOnDoubleClick(False)
        self.tree.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.tree.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.tree.setSortingEnabled(True)
        header = self.tree.header()
        header.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.Stretch)
        header.setDefaultAlignment(
            QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignVCenter
        )
        header.setSortIndicator(0, QtCore.Qt.SortOrder.AscendingOrder)
        header.setSortIndicatorShown(True)
        self.tree.itemExpanded.connect(self._on_item_expanded)
        self.tree.itemCollapsed.connect(self._on_item_collapsed)
        self.tree.itemClicked.connect(self._on_tree_clicked)
        self.tree.itemDoubleClicked.connect(self._on_tree_double_clicked)
        self.tree.setContextMenuPolicy(QtCore.Qt.ContextMenuPolicy.CustomContextMenu)
        self.tree.customContextMenuRequested.connect(self._show_sample_tree_context_menu)
        self.tree.setMouseTracking(True)
        self.tree.viewport().setMouseTracking(True)
        self.tree.viewport().installEventFilter(self)
        self.tree.setStyleSheet(SHARED_ITEM_VIEW_QSS)
        for col_idx in range(1, len(tree_headers)):
            self.tree.setColumnHidden(col_idx, col_idx not in self._fixed_tree_columns)
        tree_layout.addWidget(self.tree, stretch=1)
        split.addWidget(tree_panel)

        table_panel = QtWidgets.QWidget()
        right_layout = QtWidgets.QVBoxLayout(table_panel)
        right_layout.setContentsMargins(0, 0, 0, 0)

        table_split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        table_split.setChildrenCollapsible(False)
        right_layout.addWidget(table_split)

        self.table = QtWidgets.QTableView()
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setAlternatingRowColors(True)
        # Same frame + header as the sample tree. Keep the grid on (not `::item`
        # borders) so the model's group / separator / hit-row backgrounds still
        # paint; the light gridline-color supplies the column separators.
        self.table.setStyleSheet(_FRAME_HEADER_QSS)
        self.table.setShowGrid(True)
        self.table.setMouseTracking(True)
        self.table.viewport().setMouseTracking(True)
        self.table.viewport().installEventFilter(self)
        header = self.table.horizontalHeader()
        header.setStretchLastSection(True)
        header.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.table.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.table.verticalScrollBar().valueChanged.connect(self._on_table_scroll_changed)
        table_split.addWidget(self.table)

        table_placeholder = QtWidgets.QWidget()
        placeholder_layout = QtWidgets.QVBoxLayout(table_placeholder)
        placeholder_layout.setContentsMargins(0, 0, 0, 0)
        placeholder_layout.setSpacing(0)

        self.igv_view = IgvTrackView(table_placeholder)
        self.igv_view.geneDoubleClicked.connect(self._on_igv_gene_double_clicked)
        self.igv_view.viewCenterBpChanged.connect(self._on_igv_center_bp_changed)
        placeholder_layout.addWidget(self.igv_view, stretch=1)
        table_split.addWidget(table_placeholder)
        table_split.setSizes([600, 200])

        split.addWidget(table_panel)
        split.setStretchFactor(0, 0)
        split.setStretchFactor(1, 1)
        split.setSizes([420, 820])

        v_split.addWidget(top_panel)

        self.bottom_panel = QtWidgets.QWidget()
        self.bottom_panel.setMinimumHeight(120)
        bottom_layout = QtWidgets.QHBoxLayout(self.bottom_panel)
        bottom_layout.setContentsMargins(0, 0, 0, 0)
        bottom_layout.setSpacing(8)

        self.bottom_tabs = QtWidgets.QTabWidget()
        for name in ["Geo", "neoPan"]:
            tab = QtWidgets.QWidget()
            self.bottom_tabs.addTab(tab, name)
        self.bottom_tabs.currentChanged.connect(self._on_bottom_tab_changed)
        bottom_layout.addWidget(self.bottom_tabs, stretch=1)
        v_split.addWidget(self.bottom_panel)
        v_split.setSizes([800, 220])
        return page

    def _update_side_panel_widths(self) -> None:
        target = max(1, int(round(float(max(1, self.width())) * 0.10)))
        if self._filters_panel is not None:
            self._filters_panel.setMinimumWidth(target)
            self._filters_panel.setMaximumWidth(target)


    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        super().resizeEvent(event)
        self._update_side_panel_widths()
