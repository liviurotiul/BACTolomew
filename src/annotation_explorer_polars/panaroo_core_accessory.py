from __future__ import annotations

from html import escape as html_escape
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets

from pipeline.PopulationStructure import run_population_label

from .pan_view import PanHeatmapWidget
from .panaroo_force import CgmlstForceGraphView
from .phylo import ORDER_AS_COMPUTED, ORDER_LADDER_DEPTH, ORDER_LADDER_SIZE, ORDER_MODES
from .settings_column import (
    SettingsColumn,
    auto_fold,
    SettingsStrip,
    select_combo_data,
    settings_card,
    settings_combo,
    settings_control,
    stacked_field,
)
from .theme import icon as theme_icon
from .tree_view import LONG_BRANCH_BREAK_FACTOR, ROOT_AS_IN_FILE, ROOT_GROUPS, ROOT_MIDPOINT, ROOT_MODES, RaxmlTreeView
from .panaroo_common import (
    _file_signature,
    _natural_sort_key,
    cluster_colors as fastgear_cluster_colors,
)
from .model.pan_data import (
    SINGLE_VIEWS,
    VIEW_FORCE,
    VIEW_METADATA,
    VIEW_PLASMIDOME,
    VIEW_PRESENCE,
    VIEW_TREE,
    ROOTING_MIDPOINT,
    ROOTING_OUTGROUP,
    BranchLengthSet,
    PanData,
    TreeIndex,
    TreeSource,
    build_sample_hover_metadata,
    find_panaroo_tree,
    normalize_metadata_value,
    tree_sources as run_tree_sources,
)
from .plasmid_model import PlasmidClusterIndex, build_plasmid_cluster_index
from .plasmidome_view import PlasmidomePanel, plasmid_gene_rows
from .panaroo_data import (
    _find_fastgear_sample_clades,
    _load_fastgear_clusters,
    _load_fastgear_lineages,
)
from .panaroo_widgets import _MetadataDonutPanel
from .panaroo_style import _fastgear_lineage_colors, _log_view_event


# "Color tree by": branches and, unless ``branch_colors_only``, tip labels.
BRANCH_COLOR_NONE = "none"
BRANCH_COLOR_CLADE = "clade"
BRANCH_COLOR_CLUSTER = "cluster"
BRANCH_COLOR_MODES = (BRANCH_COLOR_NONE, BRANCH_COLOR_CLADE, BRANCH_COLOR_CLUSTER)


# The Root dropdown: the root a tree file has, its midpoint, or the edge that
# keeps the fastGEAR clades together.
ROOT_CHOICES = (
    (ROOT_AS_IN_FILE, "As in file"),
    (ROOT_MIDPOINT, "Midpoint"),
    (ROOT_GROUPS, "fastGEAR clades"),
)
# The Order dropdown: the leaves ladderized (larger clades at the top) by
# clade size or by depth, with the clade and cluster bands kept together, or
# in the order the root and the clade rotation give.
ORDER_CHOICES = (
    (ORDER_LADDER_SIZE, "By clade size"),
    (ORDER_LADDER_DEPTH, "By depth"),
    (ORDER_AS_COMPUTED, "As computed"),
)
# The id of a tree file handed to the tab (``tree_file``) instead of the run's.
TREE_FILE_SOURCE_ID = "file"
# The Lengths dropdown's values besides the branch length set ids of the
# drawn tree: a tree with fewer than two sets draws "Distance"; "Topology"
# draws every branch the same length.
LENGTHS_DISTANCE = ":distance"
LENGTHS_TOPOLOGY = ":topology"
# The Settings column beside the tree opens by itself in a panel at least
# this wide (Classic's Phylogeny is narrower and starts folded), and folds
# while the panel has no room for it beside a tree this wide.
TREE_SETTINGS_OPEN_MIN_WIDTH = 700
TREE_MIN_WIDTH = 280


class _MiddleElideComboBox(QtWidgets.QComboBox):
    """A combo whose closed text is elided in the middle (Qt elides on the
    right), so both ends of a long item stay readable."""

    def paintEvent(self, _event: QtGui.QPaintEvent) -> None:
        painter = QtWidgets.QStylePainter(self)
        option = QtWidgets.QStyleOptionComboBox()
        self.initStyleOption(option)
        painter.drawComplexControl(QtWidgets.QStyle.ComplexControl.CC_ComboBox, option)
        field = self.style().subControlRect(
            QtWidgets.QStyle.ComplexControl.CC_ComboBox,
            option,
            QtWidgets.QStyle.SubControl.SC_ComboBoxEditField,
            self,
        )
        if self.fontMetrics().horizontalAdvance(option.currentText) > field.width():
            option.currentText = self.fontMetrics().elidedText(
                option.currentText, QtCore.Qt.TextElideMode.ElideMiddle, max(0, field.width() - 4)
            )
        painter.drawControl(QtWidgets.QStyle.ControlElement.CE_ComboBoxLabel, option)


class _NotesRow(QtWidgets.QWidget):
    """The notes under the tree source (genomes missing, genomes unmasked,
    the Gubbins chip): side by side when they fit, else one under the other,
    so a narrow panel never cuts a note off."""

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.row = QtWidgets.QBoxLayout(QtWidgets.QBoxLayout.Direction.LeftToRight, self)
        self.row.setContentsMargins(0, 0, 0, 0)
        self.row.setSpacing(10)
        self._notes: List[QtWidgets.QWidget] = []

    def add_note(self, note: QtWidgets.QWidget) -> None:
        self._notes.append(note)
        self.row.addWidget(note, 0, QtCore.Qt.AlignmentFlag.AlignLeft)

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        super().resizeEvent(event)
        self.fit()

    def fit(self) -> None:
        shown = [note for note in self._notes if not note.isHidden()]
        needed = sum(note.sizeHint().width() for note in shown) + self.row.spacing() * max(0, len(shown) - 1)
        direction = (
            QtWidgets.QBoxLayout.Direction.LeftToRight
            if self.width() >= needed
            else QtWidgets.QBoxLayout.Direction.TopToBottom
        )
        if self.row.direction() != direction:
            self.row.setDirection(direction)


class PanCoreAccessoryTab(QtWidgets.QWidget):
    """The pan views of a run: (``single_view``) the one view of a pan panel
    (tree, gene presence/absence heatmap, metadata, force graph or
    plasmidome), or else a tree alone.

    Hover (window plan 8.10, Phase 9): the tab has no hover hub of its own.
    A view shows its own hover (a tree marks the hovered leaves itself, a
    heatmap or force graph highlights its hovered sample) and the tab only
    reports it: ``samplesHovered`` with every hovered sample. The panel puts
    it into its window's context, and the window's hover comes back through
    ``highlight_sample`` / ``highlight_samples``.

    A tree has a row above it (zoom and Export) and a foldable Settings
    column beside it (``settings_column``, as the gene neighborhood
    comparison): the 'Tree' card (source, root, leaf order) and the
    'Branches' card (one Lengths dropdown: the tree's branch length sets or
    Distance, and Topology; long branches at true length). The column opens
    by itself in a wide panel and folds in a narrow one until the user picks.

    With ``show_clade_legend`` a tree view gets a legend row above it: the
    fastGEAR clades (their background colors, with the number of genomes
    drawn) and clusters (their branch colors).

    Tree sources: the run's trees (``PanData.tree_index``: ``trees/trees.json``,
    or the old tree files). The tree views draw the default one (the
    recombination-filtered core tree) unless the source dropdown picks
    another; the Root dropdown's tooltip says how the drawn tree is rooted,
    and the notes above the tree which genomes it lacks, and warn when
    Gubbins did not converge.
    The root follows the file for a trees.json tree (rooted by the pipeline)
    and the midpoint for an old, unrooted file, unless the Root dropdown
    picks another.
    """

    # Every hovered sample (a list; empty for none): a clade, a cluster or a
    # plasmid hovers several at once (plan 4.6).
    samplesHovered = QtCore.pyqtSignal(object)
    cladeFilterRequested = QtCore.pyqtSignal(object)
    clusterFilterRequested = QtCore.pyqtSignal(object)
    plasmidFilterRequested = QtCore.pyqtSignal(object)
    # A clade, cluster or plasmid to open in a new explore window: the
    # payload ({"name", "samples"}) and its kind ("clade", "cluster", "plasmid").
    subsetWindowRequested = QtCore.pyqtSignal(object, str)
    plasmidContextRequested = QtCore.pyqtSignal(object)
    sampleOpenRequested = QtCore.pyqtSignal(object)
    # Whether the drawn tree has usable branch lengths, and if not why.
    treeDistanceAvailabilityChanged = QtCore.pyqtSignal(bool, str)
    # The tree drawn is another of the run's trees (the source dropdown or
    # restored display options).
    treeSourceChanged = QtCore.pyqtSignal()

    def __init__(
        self,
        output_dir: Path,
        metadata_df: Optional[pd.DataFrame],
        reference_sample: Optional[str],
        tree_file: Optional[Path] = None,
        *,
        preserve_branch_lengths: bool = False,
        color_tree_by_fastgear_clade: bool = False,
        optimize_tree_by_fastgear_clade: bool = False,
        annotate_fastgear_clusters: bool = False,
        use_fastgear_clade_backgrounds: bool = False,
        color_tree_by_fastgear_cluster: bool = False,
        show_clade_legend: bool = False,
        show_cluster_legend: bool = False,
        show_info_label: bool = True,
        enable_tree_hover_interactions: bool = True,
        ann_df: Optional[pd.DataFrame] = None,
        reroot_tree_for_clades: bool = True,
        tree_root_mode: Optional[str] = None,
        tree_order_mode: str = ORDER_LADDER_SIZE,
        branch_color_mode: str = BRANCH_COLOR_NONE,
        prune_tree_to_sample_filter: bool = False,
        scale_bar_position: str = "top",
        branch_colors_only: bool = False,
        single_view: Optional[str] = None,
        pan_data: Optional[PanData] = None,
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self._output_dir = Path(output_dir)
        # The run's shared pan data (window plan 8.10): with it, the tab takes
        # the run's plasmid index, hover lines, tree and fastGEAR lookups
        # instead of building its own, and shares the frames it is given
        # (read only) instead of copying them.
        self._pan_data = pan_data
        # The method that found the clades and clusters (fastGEAR, or FastBAPS for a run
        # made with it): the name the Root choice, the legend and their tooltips use.
        self._population_name = (
            pan_data.population_label() if pan_data is not None else run_population_label(self._output_dir)
        )
        if single_view is not None and single_view not in SINGLE_VIEWS:
            raise ValueError(f"unknown pan view {single_view!r}")
        self._single_view = single_view
        if metadata_df is None:
            self._metadata_df = pd.DataFrame()
        else:
            self._metadata_df = metadata_df if pan_data is not None else metadata_df.copy()
        self._sample_hover_metadata = (
            pan_data.sample_hover_metadata()
            if pan_data is not None
            else build_sample_hover_metadata(self._metadata_df, self._output_dir)
        )
        self._reference_sample = reference_sample
        # A tree file to draw instead of the run's trees (standalone use).
        self._tree_file = Path(tree_file) if tree_file is not None else None
        self._tree_index_cache: Optional[TreeIndex] = None
        # The tree source picked in the dropdown (None: the run's default).
        self._tree_source_id: Optional[str] = None
        self._preserve_branch_lengths = bool(preserve_branch_lengths)
        # What the user asked for: Distance stays asked for while a tree
        # without branch lengths draws Topology, and comes back with a tree
        # that has them (``_preserve_branch_lengths`` is what is drawn).
        self._requested_branch_lengths = bool(preserve_branch_lengths)
        # "Long branches at true length": drawn at their true length instead
        # of shortened with a break mark (only while drawn to scale).
        self._show_true_branch_lengths = False
        # The branch length set asked for (trees.json ``branch_length_sets``,
        # such as "with_recombination"; None: the tree's default set). Kept
        # while a tree without sets is drawn.
        self._branch_length_set_id: Optional[str] = None
        self._tree_layout_mode = "rectangular"
        self._color_tree_by_fastgear_clade = bool(color_tree_by_fastgear_clade)
        self._optimize_tree_by_fastgear_clade = bool(optimize_tree_by_fastgear_clade)
        # Without it the clade optimizer only rotates children: the drawn root
        # is the file's (or ``tree_root_mode``'s).
        self._reroot_tree_for_clades = bool(reroot_tree_for_clades)
        # A sample filter prunes the tree to the filtered genomes instead of
        # greying the others.
        self._prune_tree_to_sample_filter = bool(prune_tree_to_sample_filter)
        self._scale_bar_position = str(scale_bar_position or "top")
        # The root picked (None: the source's own, ``_tree_root_mode``):
        # "groups" roots on the edge that keeps the clades together.
        self._tree_root_choice: Optional[str] = tree_root_mode if tree_root_mode in ROOT_MODES else None
        # The leaf order (ORDER_*): ladderized by clade size by default.
        self._tree_order_mode = tree_order_mode if tree_order_mode in ORDER_MODES else ORDER_LADDER_SIZE
        # "Color tree by": fastGEAR clade or cluster colored branches.
        self._branch_color_mode = branch_color_mode if branch_color_mode in BRANCH_COLOR_MODES else BRANCH_COLOR_NONE
        # With ``branch_colors_only`` display groups (cluster colored trees,
        # "Color tree by") color the branches only, never the label text.
        self._branch_colors_only = bool(branch_colors_only)
        self._annotate_fastgear_clusters = bool(annotate_fastgear_clusters)
        self._use_fastgear_clade_backgrounds = bool(use_fastgear_clade_backgrounds)
        self._color_tree_by_fastgear_cluster = bool(color_tree_by_fastgear_cluster)
        self._show_clade_legend = bool(show_clade_legend)
        self._show_cluster_legend = bool(show_cluster_legend)
        self._show_info_label = bool(show_info_label)
        self._enable_tree_hover_interactions = bool(enable_tree_hover_interactions)
        if ann_df is None:
            self._ann_df = pd.DataFrame()
        else:
            self._ann_df = ann_df if pan_data is not None else ann_df.copy()
        self._plasmid_cluster_index: PlasmidClusterIndex = (
            pan_data.plasmid_cluster_index
            if pan_data is not None
            else build_plasmid_cluster_index(output_dir=self._output_dir, ann_df=self._ann_df)
        )
        self._sample_filter: Optional[set[str]] = None
        self._loaded_tree_key: Optional[Tuple[object, ...]] = None
        self._loaded_sample_clades_sig: Optional[Tuple[str, int, int]] = None
        self._lineage_by_sample: Dict[str, str] = {}
        self._lineage_color_map: Dict[str, QtGui.QColor] = {}
        self._cluster_by_sample: Dict[str, str] = {}
        self.tree: Optional[RaxmlTreeView] = None
        self.heatmap: Optional[PanHeatmapWidget] = None
        self.metadata_view: Optional[_MetadataDonutPanel] = None
        self.force_view: Optional[CgmlstForceGraphView] = None
        self.plasmidome_view: Optional[PlasmidomePanel] = None
        self._tree_views: List[RaxmlTreeView] = []
        self._heatmap_views: List[PanHeatmapWidget] = []
        self._metadata_views: List[_MetadataDonutPanel] = []
        self._force_views: List[CgmlstForceGraphView] = []
        self._plasmidome_views: List[PlasmidomePanel] = []
        self._info: Optional[QtWidgets.QLabel] = None
        self._tree_toolbar: Optional[QtWidgets.QFrame] = None
        self._tree_settings: Optional[SettingsColumn] = None
        self._tree_settings_strip: Optional[SettingsStrip] = None
        # Folded only because the panel is small for now (``auto_fold``):
        # it opens again when the panel grows. True until the user or a saved
        # state picks, so a wide panel opens it by default.
        self._tree_settings_auto_folded = True
        self._tree_settings_picked = False
        self._tree_body_row: Optional[QtWidgets.QHBoxLayout] = None
        self._tree_layout_combo: Optional[QtWidgets.QComboBox] = None
        self._tree_lengths_combo: Optional[QtWidgets.QComboBox] = None
        self._tree_true_lengths_check: Optional[QtWidgets.QCheckBox] = None
        self._tree_source_combo: Optional[QtWidgets.QComboBox] = None
        self._tree_source_name: Optional[QtWidgets.QLabel] = None
        self._tree_root_combo: Optional[QtWidgets.QComboBox] = None
        self._tree_order_combo: Optional[QtWidgets.QComboBox] = None
        self._tree_notes: Optional[_NotesRow] = None
        self._tree_missing_label: Optional[QtWidgets.QLabel] = None
        self._tree_unmasked_label: Optional[QtWidgets.QLabel] = None
        self._gubbins_chip: Optional[QtWidgets.QLabel] = None
        self._tree_distance_available = False
        self._tree_distance_status = "Load a tree to inspect branch distances."
        self._export_button: Optional[QtWidgets.QAbstractButton] = None
        self._single_view_widget: Optional[QtWidgets.QWidget] = None

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        if self._show_info_label:
            self._info = QtWidgets.QLabel()
            self._info.setWordWrap(True)
            self._info.setProperty("role", "muted")
            self._info.setVisible(False)
            layout.addWidget(self._info)
        self._legend: Optional[QtWidgets.QLabel] = None
        # Where the view goes: beside the tree's Settings column for a tree.
        view_layout: QtWidgets.QVBoxLayout = layout
        view_parent: QtWidgets.QWidget = self
        if self._single_view in (None, VIEW_TREE):
            body, view_parent, view_layout = self._build_tree_body()
            layout.addWidget(body, stretch=1)
        _log_view_event(f"created Panaroo accessory tab for {self._output_dir}")

        if self._single_view is not None:
            # A pan panel: one view.
            view = self._build_view(self._single_view, view_parent)
            self._single_view_widget = view
            view_layout.addWidget(view, stretch=1)
            self._after_view_placed(view)
        else:
            tree = self._create_tree_view(view_parent)
            self.tree = tree
            view_layout.addWidget(tree, stretch=1)

        self.refresh()

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

    def sample_hover_metadata(self) -> Dict[str, Dict[str, str]]:
        return {
            sample: dict(fields)
            for sample, fields in self._sample_hover_metadata.items()
        }

    def _all_views(self) -> List[QtWidgets.QWidget]:
        return [
            *self._all_tree_views(),
            *self._all_heatmaps(),
            *self._all_metadata_views(),
            *self._all_force_views(),
            *self._all_plasmidome_views(),
        ]

    def single_view(self) -> Optional[str]:
        return self._single_view

    def single_view_widget(self) -> Optional[QtWidgets.QWidget]:
        """The view of a single view tab (a standalone pan panel's view)."""
        return self._single_view_widget

    def ensure_views_fresh(self) -> None:
        """Build the force graphs whose data changed while they were hidden."""
        for force_view in self._visible_views(self._all_force_views()):
            force_view.ensure_fresh()

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
        return normalize_metadata_value(value)

    def _create_tree_view(self, parent: Optional[QtWidgets.QWidget]) -> RaxmlTreeView:
        tree = RaxmlTreeView(
            parent,
            layout_mode=self._tree_layout_mode,
            preserve_branch_lengths=self._preserve_branch_lengths,
            optimize_groups=self._optimize_tree_by_fastgear_clade,
            use_group_backgrounds=self._use_fastgear_clade_backgrounds,
            show_group_labels=not self._show_clade_legend,
            show_cluster_annotations=self._annotate_fastgear_clusters,
            enable_hover_interactions=self._enable_tree_hover_interactions,
            reroot_for_groups=self._reroot_tree_for_clades,
            prune_inactive=self._prune_tree_to_sample_filter,
            root_mode=self._tree_root_mode(),
            order_mode=self._tree_order_mode,
            scale_bar_position=self._scale_bar_position,
        )
        # "Color tree by" colors the label text unless branch_colors_only.
        tree.set_display_groups_color_labels(not self._branch_colors_only)
        tree.set_show_true_branch_lengths(self._show_true_branch_lengths)
        tree.setMinimumWidth(TREE_MIN_WIDTH)
        tree.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Preferred,
            QtWidgets.QSizePolicy.Policy.Expanding,
        )
        tree.set_reference_sample(self._reference_sample)
        tree.set_leaf_hover_metadata(self._sample_hover_metadata)
        tree.cladeClicked.connect(self._on_clade_clicked)
        tree.clusterClicked.connect(self._on_cluster_clicked)
        tree.cladeWindowRequested.connect(lambda payload: self._on_subset_window_requested(payload, "clade"))
        tree.clusterWindowRequested.connect(lambda payload: self._on_subset_window_requested(payload, "cluster"))
        tree.leafDoubleClicked.connect(self._on_tree_leaf_double_clicked)
        tree.branchDistanceAvailabilityChanged.connect(
            self._on_tree_distance_availability_changed
        )
        tree.branchScaleModeChanged.connect(self._on_tree_scale_mode_changed)
        tree.branchBreaksAvailabilityChanged.connect(lambda _available: self._update_true_lengths_control())
        if self._enable_tree_hover_interactions:
            tree.leafHovered.connect(self._on_tree_leaf_hover)
            tree.subtreeHovered.connect(self._on_tree_subtree_hover)
        self._tree_views.append(tree)
        return tree

    def _build_tree_body(self) -> Tuple[QtWidgets.QWidget, QtWidgets.QWidget, QtWidgets.QVBoxLayout]:
        """The tree's Settings column (or its folded strip) beside the tree
        area: the zoom and Export row, the notes, the legend, then the tree
        (added by the caller to the returned area and layout)."""
        body = QtWidgets.QWidget(self)
        row = QtWidgets.QHBoxLayout(body)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(0)
        strip = SettingsStrip(body)
        strip.expandRequested.connect(lambda: self.set_tree_settings_open(True))
        row.addWidget(strip)
        column = self._build_tree_settings(body)
        row.addWidget(column)
        area = QtWidgets.QWidget(body)
        area_layout = QtWidgets.QVBoxLayout(area)
        area_layout.setContentsMargins(8, 0, 0, 0)
        area_layout.setSpacing(6)
        self._tree_toolbar = self._build_tree_toolbar(area)
        area_layout.addWidget(self._tree_toolbar)
        self._tree_notes = self._build_tree_notes(area)
        area_layout.addWidget(self._tree_notes)
        if self._show_clade_legend:
            self._legend = self._build_tree_legend(area)
            area_layout.addWidget(self._legend)
        row.addWidget(area, stretch=1)
        self._tree_settings = column
        self._tree_settings_strip = strip
        self._tree_body_row = row
        # Folded until the tab is shown at its size (``_fit_tree_settings``).
        column.hide()
        self._tree_settings_fit_timer = QtCore.QTimer(self)
        self._tree_settings_fit_timer.setSingleShot(True)
        self._tree_settings_fit_timer.setInterval(0)
        self._tree_settings_fit_timer.timeout.connect(self._fit_tree_settings)
        return body, area, area_layout

    def _build_tree_toolbar(self, parent: QtWidgets.QWidget) -> QtWidgets.QFrame:
        """The row above the tree: Fit, Zoom out, Zoom in and Export."""
        toolbar = QtWidgets.QFrame(parent)
        toolbar.setObjectName("TreeToolbar")
        toolbar.setAttribute(QtCore.Qt.WidgetAttribute.WA_StyledBackground, True)
        row = QtWidgets.QHBoxLayout(toolbar)
        row.setContentsMargins(8, 5, 8, 5)
        row.setSpacing(6)
        row.addWidget(self._build_zoom_pill(toolbar))
        row.addStretch(1)
        export_button = QtWidgets.QToolButton(toolbar)
        export_button.setText("Export…")
        export_button.setToolTip("Save the current tree as a high-quality PDF, PNG, or JPEG")
        export_button.clicked.connect(self.export_tree_figure)
        row.addWidget(export_button)
        self._export_button = export_button
        return toolbar

    def _build_tree_settings(self, parent: QtWidgets.QWidget) -> SettingsColumn:
        """The Settings column beside the tree: the 'Tree' card (source,
        root, leaf order) and the 'Branches' card (lengths, long branches)."""
        column = SettingsColumn(parent, object_name="FilterRail", max_width=300)
        column.collapseRequested.connect(lambda: self.set_tree_settings_open(False))
        column.add_card(self._build_tree_card(column.content))
        column.add_card(self._build_branches_card(column.content))
        column.finish(min_width=248)
        # The neighbourhood view's default column width (its widest card if
        # wider); painted while it lies over the tree in a narrow panel.
        column.setFixedWidth(column.maximumWidth())
        column.setAttribute(QtCore.Qt.WidgetAttribute.WA_StyledBackground, True)
        return column

    def _build_tree_card(self, parent: QtWidgets.QWidget) -> QtWidgets.QFrame:
        card, layout = settings_card(parent, "Tree")
        # The tree source: a dropdown of the run's trees, or the one tree's name.
        source = QtWidgets.QWidget(card)
        source_layout = QtWidgets.QVBoxLayout(source)
        source_layout.setContentsMargins(0, 0, 0, 0)
        source_layout.setSpacing(0)
        source_name = QtWidgets.QLabel(source)
        source_name.setObjectName("TreeSourceName")
        source_name.setWordWrap(True)
        source_layout.addWidget(source_name)
        source_combo = _MiddleElideComboBox(source)
        source_combo.setObjectName("TreeSourceCombo")
        source_combo.setAccessibleName("Tree source")
        # The column's width: the tree names are elided in the middle, keeping
        # the end that tells the trees apart ("recombination-filtered", "not
        # recombination-corrected"), and the full name is in the tooltip.
        source_combo.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        source_combo.setMinimumContentsLength(6)
        source_combo.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed)
        source_combo.setVisible(False)
        source_combo.currentIndexChanged.connect(self._on_tree_source_control_changed)
        source_layout.addWidget(settings_control(source_combo))
        stacked_field(card, layout, "Tree", source)

        if not self._annotation_rich_tree():
            layout_tooltip = "Choose the phylogenetic tree layout"
            layout_combo = settings_combo(
                card,
                [("rectangular", "Rectangular", layout_tooltip), ("circular", "Circular", layout_tooltip)],
                self._tree_layout_mode,
                self._on_tree_layout_control_changed,
            )
            stacked_field(card, layout, "Layout", layout_combo)
            self._tree_layout_combo = layout_combo

        root_combo = _MiddleElideComboBox(card)
        root_combo.setObjectName("TreeRootCombo")
        root_combo.setAccessibleName("Tree root")
        for mode, text in ROOT_CHOICES:
            root_combo.addItem(f"{self._population_name} clades" if mode == ROOT_GROUPS else text, mode)
        # A long "As in file (outgroup X)" is elided in the middle.
        root_combo.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        root_combo.setMinimumContentsLength(6)
        root_combo.currentIndexChanged.connect(self._on_tree_root_control_changed)
        stacked_field(card, layout, "Root", settings_control(root_combo))

        order_tooltip = (
            "Leaf order. By clade size, By depth: the tree ladderized, the larger clade drawn above the "
            "smaller one at every branching, by number of genomes (clade size) or by the longest branch path "
            f"to a tip (depth). The {self._population_name} clades and "
            "clusters stay together; only the drawing order changes, never the tree. As computed: the order "
            "the root and the clade grouping give."
        )
        order_combo = settings_combo(
            card,
            [(mode, text, order_tooltip) for mode, text in ORDER_CHOICES],
            self._tree_order_mode,
            self.set_tree_order_mode,
        )
        order_combo.setObjectName("TreeOrderCombo")
        order_combo.setAccessibleName("Tree leaf order")
        stacked_field(card, layout, "Order", order_combo)
        # The drop-down lists exist from the start: the app style makes a
        # combo's list window see-through when it polishes the combo, and
        # creating that window during an app-wide restyle (setStyle polishes
        # every widget) crashes Qt.
        source_combo.view()
        root_combo.view()
        order_combo.view()
        self._tree_source_name = source_name
        self._tree_source_combo = source_combo
        self._tree_root_combo = root_combo
        self._tree_order_combo = order_combo
        return card

    def _build_branches_card(self, parent: QtWidgets.QWidget) -> QtWidgets.QFrame:
        card, layout = settings_card(parent, "Branches")
        # Lengths: the drawn tree's branch length sets (the core tree's
        # mutations only or including recombination; Distance for a tree
        # with one), or Topology.
        lengths_combo = QtWidgets.QComboBox(card)
        lengths_combo.setObjectName("TreeLengthsCombo")
        lengths_combo.setAccessibleName("Branch lengths")
        lengths_combo.currentIndexChanged.connect(self._on_lengths_control_changed)
        lengths_combo.view()
        stacked_field(card, layout, "Lengths", settings_control(lengths_combo))
        self._tree_lengths_combo = lengths_combo
        true_lengths = QtWidgets.QCheckBox("Long branches at true length", card)
        true_lengths.setObjectName("TreeTrueLengthsCheck")
        true_lengths.setChecked(self._show_true_branch_lengths)
        true_lengths.toggled.connect(self.set_true_branch_lengths_shown)
        layout.addWidget(true_lengths)
        self._tree_true_lengths_check = true_lengths
        self._update_lengths_control()
        return card

    def _build_tree_notes(self, parent: QtWidgets.QWidget) -> _NotesRow:
        """The notes above the tree, when there is something to say: the
        genomes the tree lacks or has not recombination-masked, and a warning
        chip when Gubbins did not converge (a recombination-filtered tree only)."""
        notes = _NotesRow(parent)
        notes.setObjectName("TreeNotes")
        notes.setContentsMargins(4, 0, 4, 0)
        missing = QtWidgets.QLabel(notes)
        missing.setObjectName("TreeMissingNote")
        missing.setProperty("role", "muted")
        missing.setVisible(False)
        notes.add_note(missing)
        unmasked = QtWidgets.QLabel(notes)
        unmasked.setObjectName("TreeUnmaskedNote")
        unmasked.setProperty("role", "muted")
        unmasked.setVisible(False)
        notes.add_note(unmasked)
        chip = QtWidgets.QLabel("Gubbins did not converge", notes)
        chip.setObjectName("GubbinsConvergenceChip")
        chip.setProperty("role", "warningChip")
        chip.setVisible(False)
        notes.add_note(chip)
        notes.row.addStretch(1)
        notes.setVisible(False)
        self._tree_missing_label = missing
        self._tree_unmasked_label = unmasked
        self._gubbins_chip = chip
        return notes

    def _annotation_rich_tree(self) -> bool:
        """The circular renderer draws topology and leaf labels only: a tree
        with clade or cluster annotations is not offered that layout, which
        would silently hide them."""
        return bool(
            self._use_fastgear_clade_backgrounds
            or self._annotate_fastgear_clusters
            or self._show_clade_legend
            or self._show_cluster_legend
        )

    # ------------------------------------------------------------ settings column
    def tree_settings_open(self) -> bool:
        """True while the tree's Settings column is open, False while folded."""
        column = self._tree_settings
        return column is not None and not column.isHidden()

    def set_tree_settings_open(self, visible: bool) -> None:
        """Open or fold the Settings column (the header chevron, the strip, a
        saved panel state); the choice stands while the panel has room. Not
        shown yet, it is applied when the panel is shown at its size."""
        self._tree_settings_picked = True
        if not self.isVisible():
            self._tree_settings_auto_folded = bool(visible)
            self._show_tree_settings(False)
            return
        self._tree_settings_auto_folded = False
        self._show_tree_settings(bool(visible))

    def _tree_settings_room(self) -> bool:
        """Room for the column beside a usable tree."""
        column = self._tree_settings
        return column is not None and self.width() >= column.width() + TREE_MIN_WIDTH + 8

    def _show_tree_settings(self, visible: bool) -> None:
        column, strip, row = self._tree_settings, self._tree_settings_strip, self._tree_body_row
        if column is None or strip is None or row is None:
            return
        # Beside the tree when there is room; opened in a narrower panel it
        # lies over the tree instead of squeezing it (the strip stays below,
        # so the tree does not move) until it is folded or the panel resized.
        overlay = visible and not self._tree_settings_room()
        docked = row.indexOf(column) >= 0
        if overlay and docked:
            row.removeWidget(column)
        elif not overlay and not docked:
            row.insertWidget(1, column)
        if bool(column.property("overlay")) != overlay:
            column.setProperty("overlay", overlay)
            column.style().unpolish(column)
            column.style().polish(column)
        strip.setVisible(not visible or overlay)
        column.setVisible(visible)
        if overlay:
            body = column.parentWidget()
            column.setGeometry(0, 0, column.width(), body.height())
            column.raise_()

    def _fit_tree_settings(self) -> None:
        """Fold the column while the panel is small and open it again when
        it grows (``auto_fold``, as the neighbourhood view): small is narrower
        than ``TREE_SETTINGS_OPEN_MIN_WIDTH`` until the user picks, then too
        narrow for the column beside the tree."""
        if self._tree_settings is None or not self.isVisible():
            return
        if self._tree_settings_picked:
            small = not self._tree_settings_room()
        else:
            small = self.width() < TREE_SETTINGS_OPEN_MIN_WIDTH
        self._tree_settings_auto_folded = auto_fold(
            small,
            self.tree_settings_open(),
            self._tree_settings_auto_folded,
            self._show_tree_settings,
        )
        if self.tree_settings_open():
            # Docked again (or laid over the tree) at the new size.
            self._show_tree_settings(True)

    def showEvent(self, event: QtGui.QShowEvent) -> None:
        super().showEvent(event)
        if self._tree_settings is not None:
            self._tree_settings_fit_timer.start()

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        super().resizeEvent(event)
        if self._tree_settings is not None:
            # Once the layout has settled at the new size.
            self._tree_settings_fit_timer.start()

    def _build_tree_legend(self, parent: QtWidgets.QWidget) -> QtWidgets.QLabel:
        """The clade and cluster legend above the tree (display only; it
        wraps instead of scrolling, so it never moves the tree's rows)."""
        legend = QtWidgets.QLabel(parent)
        legend.setObjectName("TreeLegend")
        legend.setTextFormat(QtCore.Qt.TextFormat.RichText)
        legend.setWordWrap(True)
        legend.setContentsMargins(10, 3, 10, 3)
        legend.setProperty("role", "muted")
        legend.setVisible(False)
        return legend

    @staticmethod
    def _legend_chip(text: str, color: QtGui.QColor) -> str:
        luminance = (0.299 * color.red() + 0.587 * color.green() + 0.114 * color.blue()) / 255.0
        ink = "#1f2933" if luminance > 0.62 else "#ffffff"
        return (
            f'<span style="background-color:{color.name()}; color:{ink}; font-weight:600;">'
            f"&nbsp;{html_escape(text)}&nbsp;</span>"
        )

    def tree_legend_entries(self) -> Tuple[List[Tuple[str, int, str]], List[Tuple[str, str]]]:
        """What the legend lists: (clade, genomes, color) and (cluster,
        color), for the leaves the tree draws in the sample set (greyed
        leaves are not counted), in natural order."""
        tree = self.tree
        leaves = [str(s) for s in tree.leaf_order()] if tree is not None else []
        if self._sample_filter is not None:
            leaves = [leaf for leaf in leaves if leaf in self._sample_filter]
        clade_counts: Dict[str, int] = {}
        clusters: set = set()
        for leaf in leaves:
            clade = str(self._lineage_by_sample.get(leaf, "") or "").strip()
            if clade:
                clade_counts[clade] = clade_counts.get(clade, 0) + 1
            cluster = str(self._cluster_by_sample.get(leaf, "") or "").strip()
            if cluster:
                clusters.add(cluster)
        cluster_colors = self._cluster_colors()
        clades = [
            (clade, clade_counts[clade], QtGui.QColor(self._lineage_color_map.get(clade, QtGui.QColor("#94a3b8"))).name())
            for clade in sorted(clade_counts, key=_natural_sort_key)
        ]
        cluster_rows = [
            (cluster, QtGui.QColor(cluster_colors.get(cluster, QtGui.QColor("#94a3b8"))).name())
            for cluster in sorted(clusters, key=_natural_sort_key)
        ]
        return clades, cluster_rows

    def _update_tree_legend(self) -> None:
        legend = self._legend
        if legend is None:
            return
        clades, clusters = self.tree_legend_entries()
        parts: List[str] = []
        if clades and self._color_tree_by_fastgear_clade:
            chips = "&nbsp;&nbsp;".join(
                f"{self._legend_chip(name, QtGui.QColor(color))}&nbsp;{count}" for name, count, color in clades
            )
            parts.append(f"<b>Clade</b>&nbsp;&nbsp;{chips}")
        if clusters and (self._annotate_fastgear_clusters or self._branch_color_mode == BRANCH_COLOR_CLUSTER):
            chips = "&nbsp;".join(self._legend_chip(name, QtGui.QColor(color)) for name, color in clusters)
            parts.append(f"<b>Cluster</b>&nbsp;&nbsp;{chips}")
        legend.setText("&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;".join(parts))
        legend.setToolTip(
            f"Clade: {self._population_name} lineages, the tree's background bands (number of genomes in the "
            f"sample set). Cluster: {self._population_name} clusters, the branch colors and the C brackets."
            if parts
            else ""
        )
        legend.setVisible(bool(parts))

    def _build_zoom_pill(self, parent: QtWidgets.QWidget) -> QtWidgets.QFrame:
        """Fit, Zoom out and Zoom in as one segmented pill of icon buttons."""
        pill = QtWidgets.QFrame(parent)
        pill.setObjectName("TreeZoomPill")
        pill.setAttribute(QtCore.Qt.WidgetAttribute.WA_StyledBackground, True)
        # Never squeezed to an empty box in a narrow panel.
        pill.setSizePolicy(QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed)
        pill_row = QtWidgets.QHBoxLayout(pill)
        pill_row.setContentsMargins(1, 1, 1, 1)
        pill_row.setSpacing(0)
        buttons = [
            (
                "zoom_fit",
                "Fit the tree to the panel (double-click the background, 0 or F)",
                "Fit tree",
                self._fit_tree_views,
            ),
            ("zoom_out", "Zoom out (Ctrl+mouse wheel or - also zooms)", "Zoom out tree", self._zoom_out_tree_views),
            ("zoom_in", "Zoom in (Ctrl+mouse wheel or + also zooms)", "Zoom in tree", self._zoom_in_tree_views),
        ]
        for icon_name, tooltip, name, slot in buttons:
            button = QtWidgets.QToolButton(pill)
            button.setIcon(theme_icon(icon_name))
            button.setIconSize(QtCore.QSize(16, 16))
            button.setToolTip(tooltip)
            button.setAccessibleName(name)
            button.setAutoRaise(True)
            button.clicked.connect(slot)
            pill_row.addWidget(button)
        return pill

    def _interactive_tree_views(self) -> List[RaxmlTreeView]:
        tree_views = self._all_tree_views()
        visible = [tree for tree in tree_views if tree.isVisible()]
        return visible or tree_views

    def _active_tree_view(self) -> Optional[RaxmlTreeView]:
        tree_views = self._interactive_tree_views()
        if self.tree in tree_views:
            return self.tree
        return tree_views[0] if tree_views else self.tree

    def _fit_tree_views(self) -> None:
        for tree in self._interactive_tree_views():
            tree.reset_view()

    def fit_tree_views(self) -> None:
        """Fit the visible tree views to their panel (the Fit button)."""
        self._fit_tree_views()

    def _zoom_in_tree_views(self) -> None:
        for tree in self._interactive_tree_views():
            tree.zoom_in()

    def _zoom_out_tree_views(self) -> None:
        for tree in self._interactive_tree_views():
            tree.zoom_out()

    def _on_tree_layout_control_changed(self, value: object) -> None:
        mode = str(value or "rectangular").strip().lower()
        if mode not in {"rectangular", "circular"}:
            return
        self._tree_layout_mode = mode
        self._update_true_lengths_control()
        self._update_lengths_control()
        for tree in self._all_tree_views():
            tree.set_layout_mode(mode)

    def true_branch_lengths_shown(self) -> bool:
        """True when long branches are drawn at their true length (the "Long
        branches at true length" checkbox), False when they are shortened
        and marked."""
        return self._show_true_branch_lengths

    def set_true_branch_lengths_shown(self, enabled: bool) -> None:
        enabled = bool(enabled)
        check = self._tree_true_lengths_check
        if check is not None and check.isChecked() != enabled:
            with QtCore.QSignalBlocker(check):
                check.setChecked(enabled)
        if enabled == self._show_true_branch_lengths:
            return
        self._show_true_branch_lengths = enabled
        for tree in self._all_tree_views():
            tree.set_show_true_branch_lengths(enabled)

    def _update_true_lengths_control(self) -> None:
        """The toggle is enabled while the Distance drawing has long branches
        to break."""
        check = self._tree_true_lengths_check
        if check is None:
            return
        distance = self._tree_layout_mode == "rectangular" and self._preserve_branch_lengths
        breaks = any(tree.branch_breaks_available() for tree in self._all_tree_views())
        check.setEnabled(distance and breaks)
        if not distance:
            check.setToolTip("Long branches are shortened on the Distance scale only.")
        elif not breaks:
            check.setToolTip("No branch in this tree is long enough to be shortened.")
        else:
            check.setToolTip(
                f"True branch lengths. Off: a branch over {LONG_BRANCH_BREAK_FACTOR:g} times the typical root to tip depth "
                "of the rest of the tree is drawn shortened, with a break mark and its true length. On: every branch at "
                "its true length."
            )

    def _on_tree_distance_availability_changed(
        self,
        available: bool,
        message: str,
    ) -> None:
        self._tree_distance_available = bool(available)
        self._tree_distance_status = str(message)
        self._update_true_lengths_control()
        self._update_lengths_control()
        self.treeDistanceAvailabilityChanged.emit(self._tree_distance_available, self._tree_distance_status)

    def tree_distance_status(self) -> Tuple[bool, str]:
        """Whether the drawn tree has usable branch lengths, and the reason
        shown when it has not."""
        return self._tree_distance_available, self._tree_distance_status

    def branch_lengths_shown(self) -> bool:
        """True while branches are drawn to scale (Distance), False for
        Topology."""
        return self._preserve_branch_lengths

    def branch_lengths_requested(self) -> bool:
        """True when the user asked for Distance, even while a tree without
        branch lengths draws Topology (what a panel saves)."""
        return self._requested_branch_lengths

    def set_branch_lengths_shown(self, distance: bool) -> bool:
        """Draw branches to scale (Distance) or as a cladogram (Topology).
        Distance falls back to Topology when the tree has no usable lengths
        and is drawn again once a tree with lengths is loaded; returns what
        is drawn now."""
        requested = bool(distance)
        self._requested_branch_lengths = requested
        effective = requested
        for tree in self._all_tree_views():
            effective = tree.set_preserve_branch_lengths(requested) and effective
        self._on_tree_scale_mode_changed(effective)
        return self._preserve_branch_lengths

    def _on_tree_scale_mode_changed(self, enabled: bool) -> None:
        self._preserve_branch_lengths = bool(enabled)
        self._update_true_lengths_control()
        self._update_lengths_control()

    def _shown_length_set(self, source: Optional[TreeSource] = None) -> Optional[BranchLengthSet]:
        """The branch length set drawn for ``source`` (the drawn tree by
        default): the one asked for, else its default; None for a tree
        without sets."""
        source = source or self.current_tree_source()
        if source is None or not source.branch_length_sets:
            return None
        return source.branch_length_set(self._branch_length_set_id)

    def _lengths_error(self) -> str:
        """Why the other fit's lengths asked for are not drawn (empty: they
        are, or none were asked for)."""
        return next((tree.branch_lengths_error() for tree in self._all_tree_views() if tree.branch_lengths_error()), "")

    def _drawn_length_set(self, source: Optional[TreeSource] = None) -> Optional[BranchLengthSet]:
        """The branch length set actually drawn: the one asked for, or the
        default when its file was refused (another tree, unreadable)."""
        source = source or self.current_tree_source()
        shown = self._shown_length_set(source)
        if shown is not None and source is not None and self._lengths_error():
            return source.branch_length_sets[0]
        return shown

    def branch_length_set_id(self) -> Optional[str]:
        """The id of the branch length set drawn (None: the tree has none)."""
        drawn = self._drawn_length_set()
        return drawn.id if drawn is not None else None

    def branch_length_text(self) -> str:
        """Which lengths are drawn ("Branch lengths: including
        recombination; root and leaf order from the mutations only
        lengths"); empty for a tree without branch length sets."""
        source = self.current_tree_source()
        drawn = self._drawn_length_set(source)
        if drawn is None or source is None:
            return ""
        text = f"Branch lengths: {drawn.label[:1].lower()}{drawn.label[1:]}"
        default = source.branch_length_sets[0]
        if drawn != default:
            # The root (midpoint too) and the leaf order are the default
            # lengths', so rows never move when the lengths change.
            text += f"; root and leaf order from the {default.label[:1].lower()}{default.label[1:]} lengths"
        return text

    def set_branch_length_set(self, set_id: Optional[str]) -> None:
        """Draw the branch lengths of the set ``set_id`` of the drawn tree
        (None or its default id: the default). Only the lengths change: the
        leaf order, root, clusters and support stay; a tree without that set
        draws its default and the choice is kept for one with it."""
        source = self.current_tree_source()
        value = str(set_id).strip() if set_id else None
        if source is not None and source.branch_length_sets and value == source.branch_length_sets[0].id:
            value = None
        before = self._shown_length_set()
        self._branch_length_set_id = value
        shown = self._shown_length_set()
        if source is not None and shown is not None and shown != before and self._loaded_tree_key is not None:
            lengths_path = self._other_lengths_path(source, shown)
            for tree_view in self._all_tree_views():
                tree_view.set_branch_lengths(
                    lengths_path,
                    label=shown.label,
                    own_label=source.branch_length_sets[0].label,
                    unit=shown.branch_length_unit or source.branch_length_unit,
                )
            self._loaded_tree_key = self._tree_load_key(source, lengths_path)
        self._update_lengths_control()

    @staticmethod
    def _other_lengths_path(source: TreeSource, lengths: Optional[BranchLengthSet]) -> Optional[Path]:
        """The file of another fit's lengths to draw on ``source`` (None:
        the tree file's own)."""
        return lengths.path if lengths is not None and lengths.path != source.path else None

    @staticmethod
    def _tree_load_key(source: TreeSource, lengths_path: Optional[Path]) -> Tuple[object, ...]:
        """What the drawn tree was loaded from: the tree file and the other
        fit's lengths file, so a refit computed or redone for an unchanged
        tree is drawn on the next load."""
        lengths_signature: Tuple[object, ...] = ()
        if lengths_path is not None:
            try:
                lengths_signature = _file_signature(lengths_path)
            except OSError:
                lengths_signature = (str(lengths_path),)
        return (source.id, *_file_signature(source.path), *lengths_signature)

    def _update_lengths_control(self) -> None:
        """The Lengths dropdown: the drawn tree's branch length sets when it
        has two or more (else Distance), then Topology; the drawn one shown."""
        combo = self._tree_lengths_combo
        if combo is None:
            return
        source = self.current_tree_source()
        sets = source.branch_length_sets if source is not None else ()
        # The other fit's file could not be drawn (another tree): its own
        # lengths are, and its item says why.
        error = self._lengths_error()
        options: List[Tuple[str, str, str]] = []
        if len(sets) >= 2:
            for entry in sets:
                details = [f"Branch lengths: {entry.label}."]
                if entry.alignment:
                    details.append(f"Fitted on {entry.alignment}.")
                if entry.model:
                    details.append(f"Model: {entry.model}.")
                details.append(
                    "Same tree (topology, root, support and leaf order); only the branch lengths change."
                )
                if error and entry is not sets[0]:
                    details.append(error)
                options.append((entry.id, entry.label, " ".join(details)))
        else:
            options.append((LENGTHS_DISTANCE, "Distance", "Distance: cumulative Newick branch lengths."))
        options.append((LENGTHS_TOPOLOGY, "Topology", "Topology: equal spacing by branching depth."))
        if not self._preserve_branch_lengths:
            current = LENGTHS_TOPOLOGY
        elif len(sets) >= 2:
            drawn = self._drawn_length_set(source)
            current = drawn.id if drawn is not None else sets[0].id
        else:
            current = LENGTHS_DISTANCE
        with QtCore.QSignalBlocker(combo):
            if [(combo.itemData(i), combo.itemText(i)) for i in range(combo.count())] != [
                (value, text) for value, text, _tip in options
            ]:
                combo.clear()
                for value, text, _tip in options:
                    combo.addItem(text, value)
            model = combo.model()
            for row, (value, _text, tooltip) in enumerate(options):
                combo.setItemData(row, tooltip, QtCore.Qt.ItemDataRole.ToolTipRole)
                item = model.item(row) if isinstance(model, QtGui.QStandardItemModel) else None
                if item is not None:
                    item.setEnabled(not (error and len(sets) >= 2 and value not in (sets[0].id, LENGTHS_TOPOLOGY)))
            combo.setCurrentIndex(max(0, combo.findData(current)))
        rectangular = self._tree_layout_mode == "rectangular"
        combo.setEnabled(rectangular and self._tree_distance_available)
        if not rectangular:
            combo.setToolTip("Branch scaling is available in the rectangular tree layout.")
        elif not self._tree_distance_available:
            combo.setToolTip(self._tree_distance_status)
        else:
            combo.setToolTip(combo.itemData(combo.currentIndex(), QtCore.Qt.ItemDataRole.ToolTipRole) or "")

    def _on_lengths_control_changed(self, _index: int) -> None:
        combo = self._tree_lengths_combo
        if combo is None:
            return
        value = combo.currentData()
        if value == LENGTHS_TOPOLOGY:
            self.set_branch_lengths_shown(False)
            return
        if value != LENGTHS_DISTANCE:
            self.set_branch_length_set(value)
        self.set_branch_lengths_shown(True)

    def tree_display_options(self) -> Dict[str, Any]:
        """The tree's display options (a Phylogeny panel's params): layout,
        branch scale, and the tree source, root, order, branch lengths, long
        branches, the Settings column's fold when picked, and the collapsed
        subtrees (each as its sorted leaf names) when any."""
        options: Dict[str, Any] = {
            "layout": self._tree_layout_mode,
            "branch_scale": "distance" if self._requested_branch_lengths else "topology",
        }
        if self._tree_source_id is not None:
            options["source"] = self._tree_source_id
        if self._tree_root_choice is not None:
            options["root"] = self._tree_root_choice
        if self._tree_order_mode != ORDER_LADDER_SIZE:
            options["order"] = self._tree_order_mode
        if self._show_true_branch_lengths:
            options["true_lengths"] = "on"
        if self._branch_length_set_id is not None:
            options["branch_lengths"] = self._branch_length_set_id
        if self._tree_settings_picked:
            wanted = self.tree_settings_open() or self._tree_settings_auto_folded
            options["settings"] = "open" if wanted else "folded"
        tree = self._active_tree_view()
        collapsed = tree.collapsed_groups() if tree is not None else []
        if collapsed:
            options["collapsed"] = [list(group) for group in collapsed]
        return options

    def set_tree_display_options(self, options: Dict[str, Any]) -> None:
        """Apply saved tree display options (unknown or unavailable values
        are ignored). The Lengths dropdown shows ``branch_scale`` (Topology)
        or else ``branch_lengths`` (a set id), as saved before it replaced
        the Branch scale dropdown and the Branch lengths switch."""
        layout = str(options.get("layout") or "").strip().lower()
        combo = self._tree_layout_combo
        if combo is not None and layout and layout != self._tree_layout_mode:
            index = combo.findData(layout)
            if index >= 0:
                combo.setCurrentIndex(index)
        root = str(options.get("root") or "").strip().lower()
        root = root if root in ROOT_MODES else None
        # An id this run does not list draws the default.
        source_id = self._listed_tree_source_id(str(options.get("source") or "") or None)
        if "source" in options and source_id != self._tree_source_id:
            if root is not None:
                # Picked before the switch, so the new tree loads once, rooted.
                self._tree_root_choice = root
            self.select_tree_source(source_id)
        elif root is not None:
            self.set_tree_root_mode(root)
        self.set_true_branch_lengths_shown(str(options.get("true_lengths") or "").strip().lower() == "on")
        self.set_branch_length_set(str(options.get("branch_lengths") or "").strip() or None)
        order = str(options.get("order") or "").strip().lower()
        self.set_tree_order_mode(order if order in ORDER_MODES else ORDER_LADDER_SIZE)
        scale = str(options.get("branch_scale") or "").strip().lower()
        if scale in {"distance", "topology"}:
            self.set_branch_lengths_shown(scale == "distance")
        settings = str(options.get("settings") or "").strip().lower()
        if settings in {"open", "folded"}:
            self.set_tree_settings_open(settings == "open")
        # Last, on the tree drawn with the options above: sets it has no
        # subtree for are not drawn. A saved layout has no "collapsed" (it is
        # run data): applying one leaves the collapsed subtrees as they are.
        collapsed = options.get("collapsed")
        if isinstance(collapsed, list):
            groups = [group for group in collapsed if isinstance(group, (list, tuple))]
            for tree in self._all_tree_views():
                tree.set_collapsed_groups(groups)

    def _cluster_colors(self) -> Dict[str, QtGui.QColor]:
        """fastGEAR cluster colors: the app's cluster palette, in cluster order."""
        return {
            cluster: QtGui.QColor(color)
            for cluster, color in fastgear_cluster_colors(list(self._cluster_by_sample.values()), sort_key=_natural_sort_key).items()
        }

    def _build_view(self, view_type: str, parent: QtWidgets.QWidget) -> QtWidgets.QWidget:
        """Build the pan view of a single view tab, filtered to the current
        sample set."""
        selected = sorted(self._sample_filter) if self._sample_filter is not None else None
        if view_type == VIEW_TREE:
            new_view: QtWidgets.QWidget = self._create_tree_view(parent)
            if self.tree is None:
                self.tree = new_view  # type: ignore[assignment]
        elif view_type == VIEW_PRESENCE:
            heatmap = PanHeatmapWidget(self._output_dir, parent=parent, auto_refresh=False)
            self._heatmap_views.append(heatmap)
            heatmap.sampleHovered.connect(self._on_heatmap_sample_hover)
            self.heatmap = self.heatmap or heatmap
            tree_order = self._heatmap_leaf_order()
            heatmap.set_sample_order(list(reversed(tree_order)) if tree_order else None)
            heatmap.set_filters(selected, None, None)
            heatmap.refresh()
            new_view = heatmap
        elif view_type == VIEW_METADATA:
            metadata_view = _MetadataDonutPanel(self._metadata_df, parent=parent)
            self._metadata_views.append(metadata_view)
            self.metadata_view = self.metadata_view or metadata_view
            metadata_view.set_sample_filter(selected)
            new_view = metadata_view
        elif view_type == VIEW_FORCE:
            force_view = CgmlstForceGraphView(
                self._output_dir,
                reference_sample=self._reference_sample,
                parent=parent,
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
        elif view_type == VIEW_PLASMIDOME:
            plasmidome = PlasmidomePanel(
                self._ann_df,
                self._output_dir,
                self._plasmid_cluster_index,
                parent=parent,
                # Built once per run and shared by every Plasmidome view.
                gene_rows=(
                    self._pan_data.memo(("plasmid_gene_rows",), lambda: plasmid_gene_rows(self._ann_df))
                    if self._pan_data is not None
                    else None
                ),
            )
            self._plasmidome_views.append(plasmidome)
            self.plasmidome_view = self.plasmidome_view or plasmidome
            plasmidome.samplesHovered.connect(self._on_plasmidome_samples_hover)
            plasmidome.sampleSubsetRequested.connect(self._on_plasmid_clicked)
            plasmidome.sampleSubsetWindowRequested.connect(
                lambda payload: self._on_subset_window_requested(payload, "plasmid")
            )
            plasmidome.contextRequested.connect(self._on_plasmid_context_requested)
            plasmidome.set_sample_filter(selected)
            new_view = plasmidome
        else:
            raise ValueError(f"unknown pan view {view_type!r}")
        return new_view

    def _after_view_placed(self, view: QtWidgets.QWidget) -> None:
        if isinstance(view, RaxmlTreeView):
            selected = sorted(self._sample_filter) if self._sample_filter is not None else None
            self._apply_tree_source()
            view.set_active_samples(selected)
            self._update_tree_legend()

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
            self._update_tree_legend()
        for heatmap_view in self._all_heatmaps():
            heatmap_view.set_filters(selected, None, None)
        for metadata_view in self._all_metadata_views():
            metadata_view.set_sample_filter(selected)
        for force_view in self._all_force_views():
            force_view.set_sample_filter(selected, refresh=force_view.isVisible())
        for plasmidome_view in self._all_plasmidome_views():
            plasmidome_view.set_sample_filter(selected)

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
            self._show_samples_in(heatmap_view, normalized)
        plasmidome_views = (
            self._visible_views(self._all_plasmidome_views()) if visible_only else self._all_plasmidome_views()
        )
        for plasmidome_view in plasmidome_views:
            plasmidome_view.set_linked_samples(normalized)

    @staticmethod
    def _show_samples_in(view: QtWidgets.QWidget, samples: Sequence[str]) -> None:
        """Highlight ``samples`` in a heatmap or force graph (none clears)."""
        if isinstance(view, PanHeatmapWidget):
            if len(samples) >= 2:
                view.highlight_sample(None)
                view.highlight_samples(list(samples))
            elif samples:
                view.highlight_samples(None)
                view.highlight_sample(samples[0])
            else:
                view.highlight_samples(None)
                view.highlight_sample(None)
        elif len(samples) >= 2:
            view.highlight_samples(list(samples))
        else:
            view.highlight_sample(samples[0] if samples else None)

    def _report_hover(self, samples: Optional[Sequence[str]]) -> None:
        """A view of this tab was hovered: tell the panel (no hover hub here;
        the window's bus links the panels, see the class docstring)."""
        normalized = [str(s).strip() for s in (samples or []) if str(s).strip()]
        self.samplesHovered.emit(list(normalized))

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

    def _on_subset_window_requested(self, payload: object, kind: str) -> None:
        if not isinstance(payload, dict):
            return
        self.subsetWindowRequested.emit(dict(payload), str(kind))

    def _on_plasmid_context_requested(self, payload: object) -> None:
        if not isinstance(payload, dict):
            return
        self.plasmidContextRequested.emit(dict(payload))

    def _on_plasmidome_samples_hover(self, samples: object) -> None:
        values = samples if isinstance(samples, (list, tuple)) else None
        self._report_hover(values)

    def _on_tree_leaf_double_clicked(self, sample: object) -> None:
        value = str(sample).strip() if sample is not None else ""
        if value:
            self.sampleOpenRequested.emit(value)

    def export_tree_figure(self) -> None:
        tree = self._active_tree_view()
        if tree is None:
            return

        default_dir = self._output_dir if self._output_dir.exists() else Path.home()
        source = self.current_tree_source()
        # The file name and a caption say which tree and root the figure shows.
        lengths = self._drawn_length_set(source) if self._preserve_branch_lengths else None
        stem = f"tree_{source.id}" if source is not None else "tree_export"
        if lengths is not None and source is not None and lengths != source.branch_length_sets[0]:
            stem += f"_{lengths.id}"
        default_path = default_dir / f"{stem}.pdf"
        caption = f"{source.label}. {self.tree_root_text()}. {self.tree_order_text()}" if source is not None else None
        if caption and lengths is not None:
            caption = f"{caption}. {self.branch_length_text()}"
        if tree.branch_breaks_drawn():
            caption = f"{caption}. Long branches shortened (marked)." if caption else "Long branches shortened (marked)."
        selected_path, selected_filter = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Export tree figure",
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
            tree.export_figure(export_path, caption=caption)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(
                self,
                "Export failed",
                f"Could not export the tree figure:\n{exc}",
            )
            return

        QtWidgets.QMessageBox.information(
            self,
            "Export complete",
            f"Saved the tree figure to:\n{export_path}",
        )

    def _load_fastgear_clades(
        self,
    ) -> Tuple[Dict[str, str], Dict[str, QtGui.QColor], Dict[str, str], Optional[Path]]:
        path = (
            self._pan_data.fastgear_clades_path()
            if self._pan_data is not None
            else _find_fastgear_sample_clades(self._output_dir)
        )
        if path is None:
            self._loaded_sample_clades_sig = None
            self._lineage_by_sample = {}
            self._lineage_color_map = {}
            self._cluster_by_sample = {}
            return {}, {}, {}, None

        sig = _file_signature(path)
        if self._loaded_sample_clades_sig != sig:
            if self._pan_data is not None:
                lineages, clusters = self._pan_data.fastgear_assignments()
                lineage_by_sample, cluster_by_sample = dict(lineages), dict(clusters)
            else:
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
                f"loaded {self._population_name} assignments: {path} "
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
        """Panaroo's tree (only looked for when the run has no trees.json and
        none of the old tree files; the tab keeps the result, ``tree_index``)."""
        if self._pan_data is not None:
            return self._pan_data.panaroo_tree_path()
        return find_panaroo_tree(self._output_dir)

    # ------------------------------------------------------------ tree sources
    def tree_index(self) -> TreeIndex:
        """The trees the tab can draw and the default one: the run's
        (``PanData.tree_index``), or the one ``tree_file`` it was given."""
        if self._tree_file is not None:
            if not self._tree_file.is_file():
                return TreeIndex()
            source = TreeSource(id=TREE_FILE_SOURCE_ID, label=self._tree_file.name, path=self._tree_file, rooted=True)
            return TreeIndex(sources=(source,), default_id=source.id)
        if self._pan_data is not None:
            return self._pan_data.tree_index()
        if self._tree_index_cache is None:
            self._tree_index_cache = run_tree_sources(self._output_dir, self._find_panaroo_tree)
        return self._tree_index_cache

    def current_tree_source(self) -> Optional[TreeSource]:
        """The tree the views draw: the one picked in the source dropdown,
        else the run's default (a picked id the run does not list falls back
        to the default)."""
        index = self.tree_index()
        return index.source(self._tree_source_id) or index.default()

    def tree_source_id(self) -> Optional[str]:
        source = self.current_tree_source()
        return source.id if source is not None else None

    def select_tree_source(self, source_id: Optional[str]) -> None:
        """Draw the run's tree ``source_id`` (None or an id the run does not
        list: the default)."""
        value = self._listed_tree_source_id(source_id)
        if value == self._tree_source_id:
            return
        self._tree_source_id = value
        if self._all_tree_views():
            self._apply_tree_source()
        self._apply_heatmap_order()
        self.treeSourceChanged.emit()

    def _listed_tree_source_id(self, source_id: Optional[str]) -> Optional[str]:
        """``source_id`` as stored: None for the default or an id the run
        does not list."""
        index = self.tree_index()
        value = str(source_id) if source_id and index.source(str(source_id)) is not None else None
        return None if value == index.default_id else value

    def _on_tree_source_control_changed(self, _index: int) -> None:
        combo = self._tree_source_combo
        if combo is not None:
            self.select_tree_source(combo.currentData())

    def _tree_root_mode(self, source: Optional[TreeSource] = None) -> str:
        """The root drawn: the one picked in the Root dropdown, else the
        source's own (its file root for a trees.json tree, the midpoint for
        an old unrooted file). "As in file" is never used for an unrooted file."""
        source = source or self.current_tree_source()
        rooted = source is None or source.rooted
        choice = self._tree_root_choice
        if choice is None or (choice == ROOT_AS_IN_FILE and not rooted):
            return ROOT_AS_IN_FILE if rooted else ROOT_MIDPOINT
        return choice

    def tree_root_text(self) -> str:
        """How the drawn tree is rooted ("Root: midpoint", "Root: outgroup
        X", "Root: fastGEAR clades"); empty without a tree."""
        source = self.current_tree_source()
        if source is None:
            return ""
        mode = self._tree_root_mode(source)
        if mode == ROOT_GROUPS:
            return f"Root: {self._population_name} clades"
        if mode == ROOT_MIDPOINT:
            return "Root: midpoint"
        return f"Root: {self._file_root_name(source)}"

    def tree_order_text(self) -> str:
        """How the leaves are ordered ("Ladderized by clade size")."""
        if self._tree_order_mode == ORDER_LADDER_SIZE:
            return "Ladderized by clade size"
        if self._tree_order_mode == ORDER_LADDER_DEPTH:
            return "Ladderized by depth"
        return "Leaf order as computed"

    @staticmethod
    def _file_root_name(source: TreeSource) -> str:
        """The root a tree file has: "midpoint", "outgroup X", or "unrooted"
        for an old file."""
        if not source.rooted:
            return "unrooted"
        if source.rooting == ROOTING_OUTGROUP and source.outgroup:
            return f"outgroup {source.outgroup}"
        if source.rooting == ROOTING_MIDPOINT:
            return "midpoint"
        return "as in file"

    def _update_tree_controls(self) -> None:
        """The source and Root dropdowns and the notes above the tree, for
        the drawn tree."""
        index = self.tree_index()
        source = self.current_tree_source()
        several = len(index.sources) > 1
        combo = self._tree_source_combo
        if combo is not None:
            with QtCore.QSignalBlocker(combo):
                combo.clear()
                for entry in index.sources:
                    combo.addItem(entry.label, entry.id)
                    details = (entry.method, f"Model: {entry.model}." if entry.model else None, *entry.notes)
                    tooltip = " ".join(part for part in details if part)
                    combo.setItemData(combo.count() - 1, tooltip or entry.label, QtCore.Qt.ItemDataRole.ToolTipRole)
                combo.setCurrentIndex(combo.findData(source.id) if source is not None else -1)
            combo.setVisible(several)
            combo.setToolTip(f"The tree to draw: {source.label}" if source is not None else "")
        name = self._tree_source_name
        if name is not None:
            name.setText(source.label if source is not None else "")
            name.setToolTip((source.method or "") if source is not None else "")
            name.setVisible(not several)
        root_combo = self._tree_root_combo
        if root_combo is not None and source is not None:
            mode = self._tree_root_mode(source)
            model = root_combo.model()
            has_clades = len({c for c in self._lineage_by_sample.values() if str(c).strip()}) >= 2
            with QtCore.QSignalBlocker(root_combo):
                for row in range(root_combo.count()):
                    data = root_combo.itemData(row)
                    if data == ROOT_AS_IN_FILE:
                        root_combo.setItemText(row, f"As in file ({self._file_root_name(source)})")
                    enabled = (data != ROOT_AS_IN_FILE or source.rooted) and (data != ROOT_GROUPS or has_clades)
                    item = model.item(row) if isinstance(model, QtGui.QStandardItemModel) else None
                    if item is not None:
                        item.setEnabled(enabled)
                root_combo.setCurrentIndex(max(0, root_combo.findData(mode)))
            note = f" {source.rooting_note}" if source.rooting_note else ""
            root_combo.setToolTip(
                f"{self.tree_root_text()}.{note} As in file: the root the pipeline chose (midpoint, or the "
                "outgroup set for the run). Midpoint: the middle of the longest tip to tip path. "
                f"{self._population_name} clades: the edge that keeps the {self._population_name} clades together."
                if source.rooted
                else f"{self.tree_root_text()}. This tree file is unrooted: it is drawn rooted at its "
                f"midpoint unless you pick the {self._population_name} clades."
            )
        if self._tree_notes is None:
            return
        if source is None:
            self._tree_notes.setVisible(False)
            return
        missing = self._tree_missing_label
        count = len(source.samples_missing)
        if missing is not None:
            noun = "genome" if count == 1 else "genomes"
            missing.setText(f"{count} {noun} not in this tree")
            missing.setToolTip(", ".join(source.samples_missing))
            missing.setVisible(count > 0)
        unmasked = self._tree_unmasked_label
        unmasked_count = len(source.samples_unmasked)
        if unmasked is not None:
            noun = "genome" if unmasked_count == 1 else "genomes"
            unmasked.setText(f"{unmasked_count} {noun} not recombination-masked")
            # The pipeline's notes give the reason per genome group; without
            # them only the names are known.
            reasons = [note for note in source.notes if "not recombination-masked" in note]
            unmasked.setToolTip(
                "\n".join(reasons)
                if reasons
                else "Recombination was not masked in these genomes: " + ", ".join(source.samples_unmasked)
            )
            unmasked.setVisible(unmasked_count > 0)
        chip = self._gubbins_chip
        show_chip = source.recombination_filtered and index.gubbins_not_converged()
        if chip is not None:
            chip.setVisible(show_chip)
            text, tooltip = self._gubbins_chip_text(index)
            chip.setText(text)
            chip.setToolTip(tooltip if show_chip else "")
        if self._tree_notes is not None:
            self._tree_notes.setVisible(count > 0 or unmasked_count > 0 or show_chip)
            self._tree_notes.fit()

    def _gubbins_chip_text(self, index: TreeIndex) -> Tuple[str, str]:
        """The non-convergence chip's text and tooltip: per lineage, the
        lineages that did not converge and each one's iterations; else the
        one Gubbins run's."""
        def stopped(entry: Mapping[str, Any]) -> str:
            done, cap = entry.get("iterations_run"), entry.get("max_iterations")
            return f"after {done} of {cap} iterations" if done is not None and cap is not None else ""

        def lineage_line(entry: Mapping[str, Any]) -> str:
            name, after = f"Lineage {entry.get('id')}", stopped(entry)
            if entry.get("converged") is None:
                # Gubbins' log said nothing about convergence.
                ran = f" ({after.removeprefix('after ')})" if after else ""
                return f"{name}: convergence not reported{ran}."
            return f"{name} stopped {after}." if after else f"{name} did not converge."

        lineages = index.gubbins_lineages_not_converged()
        if lineages:
            ids = [str(entry.get("id")) for entry in lineages]
            noun = "lineage" if len(ids) == 1 else "lineages"
            lines = [lineage_line(entry) for entry in lineages]
            tooltip = (
                f"Gubbins ran once per {self._population_name} lineage; in "
                + ("this lineage" if len(ids) == 1 else "these lineages")
                + " it did not report convergence, so the recombination it removed may be incomplete. "
                + " ".join(lines)
            )
            return f"Gubbins did not converge ({noun} {', '.join(ids)})", tooltip
        after = stopped(index.gubbins or {})
        runs = f" It stopped {after}." if after else ""
        return "Gubbins did not converge", "Gubbins did not converge: the recombination it removed may be incomplete." + runs

    def _heatmap_leaf_order(self) -> List[str]:
        """The sample order of the presence heatmaps: the tab's tree, or (a
        single Gene presence/absence view) the order the same tree would draw."""
        if self.tree is not None:
            return [str(sample) for sample in self.tree.leaf_order()]
        if self._single_view == VIEW_PRESENCE:
            return self._reference_leaf_order()
        return []

    def _reference_leaf_order(self) -> List[str]:
        """Leaf order of the tab's tree source (the run's default for a Gene
        presence/absence panel) as the Phylogeny panel draws it with no
        filter (fastGEAR clades grouped, in the tab's leaf order), computed with a hidden tree that is
        dropped again; shared through the run's PanData (one computation per
        run and source)."""
        source = self.current_tree_source()
        if source is None:
            return []
        path = source.path
        root_mode = self._tree_root_mode(source)
        lineage_by_sample, lineage_colors, cluster_by_sample, _ = (
            self._load_fastgear_clades()
            if (self._color_tree_by_fastgear_clade or self._annotate_fastgear_clusters)
            else ({}, {}, {}, None)
        )
        key = (
            "leaf_order",
            source.id,
            _file_signature(path),
            tuple(sorted(lineage_by_sample.items())) if self._color_tree_by_fastgear_clade else (),
            tuple(sorted(cluster_by_sample.items())) if self._annotate_fastgear_clusters else (),
            self._optimize_tree_by_fastgear_clade,
            self._reroot_tree_for_clades,
            root_mode,
            self._tree_order_mode,
            str(self._reference_sample or ""),
        )

        def compute() -> Tuple[str, ...]:
            tree = RaxmlTreeView(
                self,
                layout_mode="rectangular",
                preserve_branch_lengths=self._preserve_branch_lengths,
                optimize_groups=self._optimize_tree_by_fastgear_clade,
                use_group_backgrounds=self._use_fastgear_clade_backgrounds,
                show_group_labels=not self._show_clade_legend,
                show_cluster_annotations=self._annotate_fastgear_clusters,
                enable_hover_interactions=False,
                reroot_for_groups=self._reroot_tree_for_clades,
                root_mode=root_mode,
                order_mode=self._tree_order_mode,
            )
            tree.hide()
            # Same as _create_tree_view: a snippy tree's "Reference" leaf
            # stands for the reference sample.
            tree.set_reference_sample(self._reference_sample)
            try:
                if self._color_tree_by_fastgear_clade:
                    tree.set_leaf_groups(lineage_by_sample, group_colors=lineage_colors)
                if self._annotate_fastgear_clusters:
                    tree.set_leaf_clusters(cluster_by_sample)
                tree.set_tree_file(path)
                return tuple(str(sample) for sample in tree.leaf_order())
            finally:
                tree.setParent(None)
                tree.deleteLater()

        if self._pan_data is not None:
            return list(self._pan_data.memo(key, compute))
        return list(compute())

    def _apply_tree_source(self) -> None:
        """Load the drawn tree source into the tree views, rooted as its
        root mode says. A sample filter never changes the tree file: the
        views grey or prune the filtered-out leaves (``set_active_samples``)."""
        tree_views = self._all_tree_views()
        if not tree_views:
            return
        source = self.current_tree_source()
        self._update_tree_controls()
        if source is not None:
            root_mode = self._tree_root_mode(source)
            lengths = self._shown_length_set(source)
            lengths_path = self._other_lengths_path(source, lengths)
            tree_key = self._tree_load_key(source, lengths_path)
            if self._loaded_tree_key != tree_key:
                for tree_view in tree_views:
                    tree_view.load_tree_file(
                        source.path,
                        root_mode=root_mode,
                        branch_length_unit=(lengths.branch_length_unit if lengths is not None else None)
                        or source.branch_length_unit,
                        support_kind=source.support_name(),
                        support_threshold=source.support_threshold,
                        lengths_file=lengths_path,
                        lengths_label=lengths.label if lengths is not None else None,
                        own_lengths_label=source.branch_length_sets[0].label if lengths is not None else None,
                    )
                self._loaded_tree_key = tree_key
                self._update_lengths_control()
                # A tree without lengths fell back to Topology: a tree with
                # them draws the Distance the user asked for again.
                if self._requested_branch_lengths and not self._preserve_branch_lengths:
                    self.set_branch_lengths_shown(True)
                _log_view_event(f"loaded tree {source.id} ({source.label}): {source.path}")
            else:
                for tree_view in tree_views:
                    tree_view.set_root_mode(root_mode)
            if self._show_info_label and self._info is not None:
                self._info.clear()
                self._info.setVisible(False)
            return

        for tree_view in tree_views:
            tree_view.clear_tree("No tree file found for this tab.")
        self._loaded_tree_key = None
        if self._show_info_label and self._info is not None:
            self._info.setText("This run has no phylogenetic tree.")
            self._info.setVisible(True)

    def tree_root_mode(self) -> str:
        """The root drawn (``ROOT_*``)."""
        return self._tree_root_mode()

    def set_tree_root_mode(self, mode: Optional[str]) -> None:
        """Root the tree views as in the file, at the midpoint, or on the
        edge that keeps the clades together ("groups"); None: the source's
        own root (the file's for a trees.json tree, else the midpoint)."""
        normalized = str(mode or "").strip().lower() or None
        if normalized is not None and normalized not in ROOT_MODES:
            raise ValueError(f"Unsupported tree root: {mode!r}")
        if normalized == self._tree_root_choice:
            return
        before = self._tree_root_mode()
        self._tree_root_choice = normalized
        self._update_tree_controls()
        if self._tree_root_mode() == before:
            return
        for tree_view in self._all_tree_views():
            tree_view.set_root_mode(self._tree_root_mode())
        # The leaf order changes with the root: the tracks and this tab's
        # presence heatmaps follow it.
        self._apply_heatmap_order()

    def tree_order_mode(self) -> str:
        """The leaf order drawn (``ORDER_*``)."""
        return self._tree_order_mode

    def set_tree_order_mode(self, mode: str) -> None:
        """Draw the leaves ladderized by clade size or depth, or as
        computed; the presence heatmaps follow the new order."""
        normalized = str(mode or "").strip().lower()
        if normalized not in ORDER_MODES:
            raise ValueError(f"Unsupported tree order: {mode!r}")
        if self._tree_order_combo is not None:
            select_combo_data(self._tree_order_combo, normalized)
        if normalized == self._tree_order_mode:
            return
        self._tree_order_mode = normalized
        for tree_view in self._all_tree_views():
            tree_view.set_order_mode(normalized)
        self._apply_heatmap_order()

    def _on_tree_root_control_changed(self, _index: int) -> None:
        combo = self._tree_root_combo
        if combo is not None:
            self.set_tree_root_mode(combo.currentData())

    def prune_tree_to_sample_filter(self) -> bool:
        return self._prune_tree_to_sample_filter

    def set_prune_tree_to_sample_filter(self, enabled: bool) -> None:
        """Hide the genomes outside the sample filter (prune the tree) or
        show them greyed."""
        enabled = bool(enabled)
        if enabled == self._prune_tree_to_sample_filter:
            return
        self._prune_tree_to_sample_filter = enabled
        for tree_view in self._all_tree_views():
            tree_view.set_prune_inactive(enabled)

    def branch_color_mode(self) -> str:
        return self._branch_color_mode

    def set_branch_color_mode(self, mode: str) -> None:
        """Color the branches by fastGEAR clade or cluster, or not at all
        ("none"). With ``branch_colors_only`` the labels stay neutral."""
        normalized = str(mode or "").strip().lower()
        if normalized not in BRANCH_COLOR_MODES:
            raise ValueError(f"Unsupported branch coloring: {mode!r}")
        if normalized == self._branch_color_mode:
            return
        self._branch_color_mode = normalized
        self._apply_branch_colors()
        self._update_tree_legend()

    def _apply_branch_colors(self) -> None:
        if self._branch_color_mode == BRANCH_COLOR_NONE and not self._color_tree_by_fastgear_cluster:
            for tree_view in self._all_tree_views():
                tree_view.set_leaf_display_groups(None, group_colors=None)
            return
        lineage_by_sample, lineage_colors, cluster_by_sample, _ = self._load_fastgear_clades()
        if self._branch_color_mode == BRANCH_COLOR_CLADE:
            groups, colors = lineage_by_sample, lineage_colors
        else:
            groups, colors = cluster_by_sample, self._cluster_colors()
        for tree_view in self._all_tree_views():
            tree_view.set_leaf_display_groups(groups, group_colors=colors)

    def _on_tree_leaf_hover(self, sample: object) -> None:
        # The tree marks its hovered leaf itself.
        name = str(sample).strip() if sample is not None else ""
        self._report_hover([name] if name else None)

    def _on_tree_subtree_hover(self, samples: object) -> None:
        values = samples if isinstance(samples, (list, tuple)) else None
        self._report_hover(values)

    def _on_heatmap_sample_hover(self, sample: object) -> None:
        name = str(sample).strip() if sample is not None else ""
        source = self.sender()
        if isinstance(source, PanHeatmapWidget):
            # A heatmap marks its own hovered row through its highlight.
            self._show_samples_in(source, [name] if name else [])
        self._report_hover([name] if name else None)

    def _on_force_sample_hover(self, sample: object) -> None:
        name = str(sample).strip() if sample is not None else ""
        source = self.sender()
        if isinstance(source, CgmlstForceGraphView):
            self._show_samples_in(source, [name] if name else [])
        self._report_hover([name] if name else None)

    def highlight_sample(self, sample: Optional[str], *, visible_only: bool = True) -> None:
        """Show one hovered sample (None clears) in the views."""
        value = str(sample).strip() if sample is not None else ""
        tree_views = self._visible_views(self._all_tree_views()) if visible_only else self._all_tree_views()
        for tree_view in tree_views:
            tree_view.highlight_sample(value or None)
        self._apply_detail_highlight(value or None, visible_only=visible_only)
        force_views = self._visible_views(self._all_force_views()) if visible_only else self._all_force_views()
        for force_view in force_views:
            force_view.highlight_sample(value or None)

    def highlight_samples(self, samples: Optional[Sequence[str]], *, visible_only: bool = True) -> None:
        """Show several hovered samples in the views."""
        tree_views = self._visible_views(self._all_tree_views()) if visible_only else self._all_tree_views()
        for tree_view in tree_views:
            tree_view.highlight_samples(samples)
        self._apply_detail_highlight(samples=samples, visible_only=visible_only)
        force_views = self._visible_views(self._all_force_views()) if visible_only else self._all_force_views()
        for force_view in force_views:
            force_view.highlight_samples(samples)

    def _apply_heatmap_order(self) -> None:
        """Order this tab's presence heatmaps like its tree (or, without a
        tree, like the tree it would draw)."""
        heatmaps = self._all_heatmaps()
        tree_order = self._heatmap_leaf_order() if heatmaps else []
        for heatmap_view in heatmaps:
            heatmap_view.set_sample_order(list(reversed(tree_order)) if tree_order else None)

    def refresh(self) -> None:
        tree_views = self._all_tree_views()
        lineage_by_sample: Dict[str, str] = {}
        lineage_colors: Dict[str, QtGui.QColor] = {}
        cluster_by_sample: Dict[str, str] = {}
        cluster_colors: Dict[str, QtGui.QColor] = {}
        if (
            self._color_tree_by_fastgear_clade
            or self._color_tree_by_fastgear_cluster
            or self._annotate_fastgear_clusters
        ):
            lineage_by_sample, lineage_colors, cluster_by_sample, _ = self._load_fastgear_clades()
            cluster_colors = self._cluster_colors()
        for force_view in self._all_force_views():
            force_view.set_fastgear_assignments(lineage_by_sample, lineage_colors, cluster_by_sample)
            force_view.set_sample_hover_metadata(self._sample_hover_metadata)
            force_view.set_reference_sample(self._reference_sample, refresh=False)
        for tree_view in tree_views:
            tree_view.set_leaf_hover_metadata(self._sample_hover_metadata)
            if self._color_tree_by_fastgear_cluster:
                tree_view.set_leaf_display_groups(cluster_by_sample, group_colors=cluster_colors)
            elif self._branch_color_mode == BRANCH_COLOR_NONE:
                tree_view.set_leaf_display_groups(None, group_colors=None)
            if self._color_tree_by_fastgear_clade:
                tree_view.set_leaf_groups(lineage_by_sample, group_colors=lineage_colors)
            else:
                tree_view.set_leaf_groups(None, group_colors=None)
            if self._annotate_fastgear_clusters:
                tree_view.set_leaf_clusters(cluster_by_sample)
            else:
                tree_view.set_leaf_clusters(None)
        if self._branch_color_mode != BRANCH_COLOR_NONE and not self._color_tree_by_fastgear_cluster:
            self._apply_branch_colors()
        self._apply_tree_source()

        self._apply_heatmap_order()
        selected = sorted(self._sample_filter) if self._sample_filter is not None else None
        for tree_view in tree_views:
            tree_view.set_active_samples(selected)
        self._update_tree_legend()
        for heatmap_view in self._all_heatmaps():
            heatmap_view.set_filters(selected, None, None)
        for metadata_view in self._all_metadata_views():
            metadata_view.set_sample_filter(selected)
        for force_view in self._all_force_views():
            force_view.set_sample_filter(selected, refresh=False)
            if force_view.isVisible():
                force_view.ensure_fresh()
