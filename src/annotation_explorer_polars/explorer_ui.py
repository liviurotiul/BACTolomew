from __future__ import annotations

import functools
from typing import Any, Dict, List, Optional, Tuple

from PyQt6 import QtCore, QtGui, QtWidgets

from .annotation_table import GAP_ROW_HEIGHT, RowRuleDelegate
from .igv_view import TRACK_MIN_HEIGHT, IgvTrackView
from .model.filter_state import FilterState
from .panels.annotations import AnnotationsPanel, build_view_mode_switch, build_window_options
from .panels.base import (
    ALIGNMENT,
    ANNOTATIONS,
    CGMLST_FORCE,
    GEO,
    NEIGHBOURHOOD,
    SAMPLE_TREE,
    TRACK,
    ExplorePanel,
    PanelSpec,
)
from .model.pan_data import VIEW_PRESENCE
from .panels.cgmlst_force import CgmlstForcePanel
from .panels.pan import PAN_VIEW_PANELS, PanForcePanel, PhylogenyPanel
from .panels.pan_clustering import PanClusteringPanel
from .panels.row_link import TreeRowLink
from .panels.sample_tree import SampleTreePanel
from .panels.track import TrackPanel
from .search_suggestions import attach_count_suggestions, count_suggestions
from .settings_column import SettingsColumn, SettingsStrip, settings_card, settings_control, stacked_field
from .docking import layouts as dock_layouts
from .shell.commands import (
    APP_QUIT,
    APP_RUN_LIST,
    LAYOUT_APPLY,
    LAYOUT_CLEAR_DEFAULT,
    LAYOUT_DELETE,
    LAYOUT_SAVE,
    LAYOUT_SET_DEFAULT,
    PANEL_OPEN,
    WINDOW_CLOSE,
)
from .shell.window_manager import add_menu_header
from .shell.layout_store import (
    FORMAT_VERSION,
    START_DEFAULT_LAYOUT,
    START_EMPTY_LAYOUT,
    RAIL_DEFAULT_WIDTH,
    RAIL_MAX_WIDTH,
    RAIL_MIN_WIDTH,
    WindowLayout,
    clamp_rail_width,
    run_key,
)

# Autosave of the window layout while it changes (plan 3.11).
LAYOUT_AUTOSAVE_MS = 60_000
# The first heavy panel builds wait until shortly after the main view is ready.
HEAVY_PANELS_START_MS = 700
# The busy indicator shows only for computes that take longer than this.
BUSY_INDICATOR_DELAY_MS = 300
# How long a layout confirmation stays in the status bar.
LAYOUT_MESSAGE_MS = 5000
# Height of the window color stripe under the menu bar (several windows).
WINDOW_STRIPE_HEIGHT = 4

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
# selection). Note: with a QTreeView::item rule in place the style does not
# paint QTreeWidgetItem.setBackground, so per-item background tints stay
# invisible; this block is the tree's original look minus the old 11px header
# font-size.
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


class _FilterMultiSelectComboBox(QtWidgets.QComboBox):
    """Read-only summary field whose full surface opens its checkable popup."""

    _POPUP_KEYS = {
        QtCore.Qt.Key.Key_Space,
        QtCore.Qt.Key.Key_Return,
        QtCore.Qt.Key.Key_Enter,
        QtCore.Qt.Key.Key_Down,
    }

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setEditable(True)
        self.setInsertPolicy(QtWidgets.QComboBox.InsertPolicy.NoInsert)
        self.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)

        editor = self.lineEdit()
        if editor is not None:
            editor.setReadOnly(True)
            editor.setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft)
            editor.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
            editor.installEventFilter(self)

    def eventFilter(self, watched: QtCore.QObject, event: QtCore.QEvent) -> bool:
        editor = self.lineEdit()
        if watched is editor and self.isEnabled():
            if event.type() == QtCore.QEvent.Type.MouseButtonRelease:
                if event.button() == QtCore.Qt.MouseButton.LeftButton:
                    self.showPopup()
                    event.accept()
                    return True
            elif event.type() == QtCore.QEvent.Type.KeyPress:
                if event.key() in self._POPUP_KEYS:
                    self.showPopup()
                    event.accept()
                    return True
        return super().eventFilter(watched, event)

    def showPopup(self) -> None:
        # Long predictor names would be elided in a list only as wide as the
        # field, so let the list grow to its widest row.
        view = self.view()
        scrollbar = view.verticalScrollBar().sizeHint().width()
        indicator = self.style().pixelMetric(QtWidgets.QStyle.PixelMetric.PM_IndicatorWidth)
        view.setMinimumWidth(view.sizeHintForColumn(0) + indicator + scrollbar + 40)
        super().showPopup()


class ExplorerUIMixin:
    def _build_ui(self) -> None:
        self._apply_window_title()

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)

        main_layout = QtWidgets.QVBoxLayout(central)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # The window's color, shown while several explore windows are open
        # (plan 3.10); its cards carry the same stripe.
        self._window_stripe = QtWidgets.QFrame(central)
        self._window_stripe.setObjectName("ExploreWindowStripe")
        self._window_stripe.setFixedHeight(WINDOW_STRIPE_HEIGHT)
        self._window_stripe.hide()
        main_layout.addWidget(self._window_stripe)

        content = QtWidgets.QWidget(central)
        content_layout = QtWidgets.QHBoxLayout(content)
        content_layout.setContentsMargins(8, 8, 8, 8)
        content_layout.setSpacing(0)
        main_layout.addWidget(content, stretch=1)

        # The hidden rail folds to a strip (as the Settings column of the
        # Phylogeny and neighbourhood panels); it shows exactly while Window >
        # Show Filters is off.
        filters_strip = SettingsStrip(content, "Show filters", caption="Filters")
        filters_strip.expandRequested.connect(lambda: self._set_filters_panel_visible(True))
        filters_strip.hide()
        self._filters_strip = filters_strip
        content_layout.addWidget(filters_strip)

        workspace_splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal, content)
        workspace_splitter.setObjectName("ExplorerWorkspaceSplitter")
        workspace_splitter.setChildrenCollapsible(True)
        workspace_splitter.setHandleWidth(6)
        self._workspace_splitter = workspace_splitter
        self._filters_last_width = RAIL_DEFAULT_WIDTH
        content_layout.addWidget(workspace_splitter, stretch=1)

        # The filter rail: a Settings column with the Sample, Sorting and
        # Annotations cards; its header chevron hides it like Show Filters.
        filters_panel = SettingsColumn(
            None,
            object_name="FilterRail",
            max_width=RAIL_MAX_WIDTH,
            collapse_tooltip="Hide filters",
            title="Filters",
        )
        filters_panel.collapseRequested.connect(lambda: self._set_filters_panel_visible(False))
        filters_panel.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Preferred,
            QtWidgets.QSizePolicy.Policy.Expanding,
        )
        self._filters_panel = filters_panel
        cards_parent = filters_panel.content

        meta_card, meta_layout = settings_card(cards_parent, "Sample")
        self.sample_query = settings_control(QtWidgets.QLineEdit(meta_card))
        self.sample_query.setPlaceholderText("Search samples and metadata...")
        self.sample_query.setClearButtonEnabled(True)
        meta_layout.addWidget(self.sample_query)
        # One-time hint (D3) for gene names typed into the sample search.
        self.sample_search_hint = QtWidgets.QLabel("", meta_card)
        self.sample_search_hint.setObjectName("SampleSearchHint")
        self.sample_search_hint.setWordWrap(True)
        self.sample_search_hint.setVisible(False)
        self.sample_search_hint.setProperty("role", "muted")
        self.sample_search_hint.setStyleSheet("padding: 0 0 2px 0; font-size: 12px;")
        meta_layout.addWidget(self.sample_search_hint)

        self.add_filter_btn = QtWidgets.QPushButton("Add filter", meta_card)
        self.add_filter_btn.clicked.connect(self._show_add_filter_menu)
        self.add_filter_btn.setEnabled(bool(self.run.meta_columns))
        meta_layout.addWidget(self.add_filter_btn)

        self.filters_container = QtWidgets.QWidget(meta_card)
        self.filters_layout = QtWidgets.QVBoxLayout(self.filters_container)
        self.filters_layout.setContentsMargins(0, 0, 0, 0)
        self.filters_layout.setSpacing(6)
        meta_layout.addWidget(self.filters_container)
        filters_panel.add_card(meta_card)

        sort_card, sort_layout = settings_card(cards_parent, "Sorting")
        self.sample_sort_column_box = settings_control(QtWidgets.QComboBox(sort_card))
        self.sample_sort_column_box.addItem("Sample", "__sample__")
        for col in sorted(self.run.meta_columns):
            self.sample_sort_column_box.addItem(col, col)
        self.sample_sort_order_box = settings_control(QtWidgets.QComboBox(sort_card))
        self.sample_sort_order_box.addItem("Ascending", "asc")
        self.sample_sort_order_box.addItem("Descending", "desc")
        stacked_field(sort_card, sort_layout, "Sort by", self.sample_sort_column_box)
        stacked_field(sort_card, sort_layout, "Order", self.sample_sort_order_box)
        filters_panel.add_card(sort_card)

        ann_card, ann_layout = settings_card(cards_parent, "Annotations")

        self.table_search_hint = QtWidgets.QLabel("", ann_card)
        self.table_search_hint.setVisible(False)
        self.table_search_hint.setWordWrap(True)
        self.table_search_hint.setProperty("role", "contextHint")
        ann_layout.addWidget(self.table_search_hint)

        self.gene_query = settings_control(QtWidgets.QLineEdit(ann_card))
        attach_count_suggestions(self.gene_query, lambda: count_suggestions(self.run.ann_pl, "Gene", "__gene_norm"))
        self.gene_query.setClearButtonEnabled(True)
        stacked_field(ann_card, ann_layout, "Gene/annotation contains", self.gene_query)

        self.product_query = settings_control(QtWidgets.QLineEdit(ann_card))
        attach_count_suggestions(
            self.product_query, lambda: count_suggestions(self.run.ann_pl, "Product", "__product_norm")
        )
        self.product_query.setClearButtonEnabled(True)
        stacked_field(ann_card, ann_layout, "Product contains", self.product_query)

        self.molecule_filter_box = settings_control(QtWidgets.QComboBox(ann_card))
        self.molecule_filter_box.addItem("All", "")
        self.molecule_filter_box.addItem("Plasmid only", "plasmid")
        self.molecule_filter_box.addItem("Chromosome only", "chromosome")
        stacked_field(ann_card, ann_layout, "Molecule", self.molecule_filter_box)

        self.predictor_box = settings_control(_FilterMultiSelectComboBox(ann_card))
        self.predictor_box.setAccessibleName("Predictor filter")
        self.predictor_box.setAccessibleDescription("Select one or more annotation predictors")
        self.predictor_box.setToolTip("Select one or more annotation predictors")
        self.predictor_box.setMinimumContentsLength(12)
        self.predictor_box.view().pressed.connect(self._toggle_predictor_item)
        self._init_predictor_model()
        stacked_field(ann_card, ann_layout, "Predictor", self.predictor_box)

        self.tag_box = settings_control(_FilterMultiSelectComboBox(ann_card))
        self.tag_box.setAccessibleName("Annotation tag filter")
        self.tag_box.setAccessibleDescription("Select one or more annotation tags")
        self.tag_box.setToolTip("Select one or more annotation tags")
        self.tag_box.setMinimumContentsLength(12)
        self.tag_box.view().pressed.connect(self._toggle_tag_item)
        self._init_tag_model()
        stacked_field(ann_card, ann_layout, "Tags", self.tag_box)

        # "Hide hypothetical" is in the Annotations panel's strip and "Hide
        # contigs without hits" in the Sample tree panel's (feedback
        # 2026-10-07); they stay window filters in the draft (see
        # _set_window_option). While no Annotations panel is open, this note
        # says when "Hide hypothetical" is on, so a window-wide filter is never
        # in effect unseen. "Hide contigs without hits" needs no note: it only
        # changes the Sample tree, which shows it.
        self.window_options_note = QtWidgets.QLabel("", ann_card)
        self.window_options_note.setObjectName("WindowOptionsNote")
        self.window_options_note.setWordWrap(True)
        self.window_options_note.setProperty("role", "muted")
        self.window_options_note.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.window_options_note.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.LinksAccessibleByMouse)
        self.window_options_note.setStyleSheet("padding: 2px 0 2px 0; font-size: 12px;")
        self.window_options_note.linkActivated.connect(self._open_annotations_for_options)
        self.window_options_note.setVisible(False)
        ann_layout.addWidget(self.window_options_note)
        filters_panel.add_card(ann_card)
        filters_panel.finish(min_width=RAIL_MIN_WIDTH)
        self._filters_min_width = filters_panel.minimumWidth()

        workspace_splitter.addWidget(filters_panel)
        # The rail on the left, the docking area on the right (plan 3.2).
        workspace_splitter.addWidget(self._build_docking_area())
        workspace_splitter.setCollapsible(0, True)
        workspace_splitter.setCollapsible(1, False)
        workspace_splitter.setStretchFactor(0, 0)
        workspace_splitter.setStretchFactor(1, 1)
        workspace_splitter.setSizes([RAIL_DEFAULT_WIDTH, 1112])
        workspace_splitter.splitterMoved.connect(self._on_workspace_splitter_moved)

        status_bar = self.statusBar()
        status_bar.setSizeGripEnabled(False)
        self.status = QtWidgets.QLabel("Ready", status_bar)
        self.status.setObjectName("ExplorerStatusText")
        self.status.setTextFormat(QtCore.Qt.TextFormat.PlainText)
        self.status.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        status_bar.addWidget(self.status)
        # The shared status bar: the applied stats (their "Samples: N" is the
        # size of the window's one sample set) and the hovered sample. The
        # annotation table's own lines are in the Annotations panel (D3).
        self.status.setToolTip(
            "Samples: the samples every panel of this window shows (sample search, metadata "
            "filters and subset; only samples with hits while an annotation filter is on)"
        )
        # While an Annotations panel is the active panel, its own line is
        # appended here (plan 3.9); the applied stats above never change with
        # a table's view mode.
        self.active_panel_status = QtWidgets.QLabel("", status_bar)
        self.active_panel_status.setObjectName("ExplorerActivePanelStatus")
        self.active_panel_status.setTextFormat(QtCore.Qt.TextFormat.PlainText)
        self.active_panel_status.setProperty("role", "muted")
        self.active_panel_status.setVisible(False)
        # It fills the rest of the bar, so the line follows the stats text.
        self.active_panel_status.setSizePolicy(QtWidgets.QSizePolicy.Policy.Ignored, QtWidgets.QSizePolicy.Policy.Preferred)
        status_bar.addWidget(self.active_panel_status, stretch=1)
        # Short confirmations (a layout saved, a new default, a layout that
        # could not be applied). A permanent widget, unlike showMessage(),
        # which would hide the sample count while it shows.
        self.layout_message = QtWidgets.QLabel("", status_bar)
        self.layout_message.setObjectName("ExplorerLayoutMessage")
        self.layout_message.setTextFormat(QtCore.Qt.TextFormat.PlainText)
        self.layout_message.setProperty("role", "muted")
        self.layout_message.setVisible(False)
        status_bar.addPermanentWidget(self.layout_message)
        self._layout_message_timer = QtCore.QTimer(self)
        self._layout_message_timer.setSingleShot(True)
        self._layout_message_timer.setInterval(LAYOUT_MESSAGE_MS)
        # Bound methods, not lambdas, on the window's timers: PyQt holds only a
        # weak reference to the window then, so a timer that fires while a
        # closed window is torn down cannot call a callable the garbage
        # collector has already freed (a segfault, seen in the test suite).
        self._layout_message_timer.timeout.connect(self._hide_layout_message)
        self.busy_status = QtWidgets.QLabel("Updating...", status_bar)
        self.busy_status.setObjectName("ExplorerBusyStatus")
        self.busy_status.setProperty("role", "muted")
        self.busy_status.setVisible(False)
        status_bar.addPermanentWidget(self.busy_status)
        self.hover_status = QtWidgets.QLabel("", status_bar)
        self.hover_status.setObjectName("ExplorerHoverStatus")
        self.hover_status.setTextFormat(QtCore.Qt.TextFormat.PlainText)
        status_bar.addPermanentWidget(self.hover_status)
        self._busy_timer = QtCore.QTimer(self)
        self._busy_timer.setSingleShot(True)
        self._busy_timer.setInterval(BUSY_INDICATOR_DELAY_MS)
        self._busy_timer.timeout.connect(self._show_busy_status)

        self._init_file_menu()
        self._panel_manager.install_menus(self.menuBar())
        self._init_window_menu()

        # The rail writes the context's draft: text debounces (180 ms in the
        # scheduler), the other controls apply immediately. Sorting is a tree
        # option, not a filter: it re-sorts the tree without a compute.
        for w in [self.sample_query, self.gene_query, self.product_query]:
            w.textChanged.connect(self._schedule_refresh)
        self.sample_sort_column_box.currentIndexChanged.connect(self._resort_sample_tree)
        self.sample_sort_order_box.currentIndexChanged.connect(self._resort_sample_tree)
        self.molecule_filter_box.currentIndexChanged.connect(self._apply_rail_filters)
        self._restore_explore_layout()
        QtCore.QTimer.singleShot(HEAVY_PANELS_START_MS, self._start_heavy_panels)

    # Shell menus (plan 4.2): File (Run List..., close this window, quit) and the
    # Window menu (new and duplicate windows, New Panel..., the layouts, Show
    # Filters, the open windows and cards). Everything about runs (metadata,
    # sharing, run details, settings) lives in the run list. Only New Panel...
    # has a shortcut (Ctrl+Shift+P, asked for by the user; D8 holds for
    # everything else).

    def _init_file_menu(self) -> None:
        file_menu = self.menuBar().addMenu("File")
        file_menu.setObjectName("FileMenu")
        if self.app_controller.can_open_run_list():
            run_list_action = file_menu.addAction("Run List...")
            run_list_action.setObjectName("RunListAction")
            run_list_action.setToolTip("Close this run's windows and choose a run")
            run_list_action.triggered.connect(lambda _c=False: self.app_controller.commands.run(APP_RUN_LIST, self))
            run_list_action.setMenuRole(QtGui.QAction.MenuRole.NoRole)
            file_menu.addSeparator()
        close_action = file_menu.addAction("Close Window")
        close_action.setObjectName("CloseWindowAction")
        close_action.setToolTip("Close this explore window; other explore windows stay open")
        close_action.triggered.connect(lambda _c=False: self.app_controller.commands.run(WINDOW_CLOSE, self))
        file_menu.addSeparator()
        quit_action = file_menu.addAction("Quit")
        quit_action.setObjectName("QuitAction")
        quit_action.setToolTip("Close every explore window of this run; they open again with the run")
        quit_action.triggered.connect(lambda _c=False: self.app_controller.commands.run(APP_QUIT, self))
        # Keep both out of the macOS application menu (and its native Cmd+Q, D8).
        for action in (close_action, quit_action):
            action.setMenuRole(QtGui.QAction.MenuRole.NoRole)
        self._file_menu = file_menu

    def _init_window_menu(self) -> None:
        menu = self.menuBar().addMenu("Window")
        menu.setObjectName("WindowMenu")
        # Window > New Panel... (Ctrl+Shift+P): the panel picker of this
        # window. The action belongs to the window, not the menu (the menu
        # is cleared and filled each time it opens), so its shortcut works
        # whenever this window is active, menu opened or not.
        new_panel = QtGui.QAction("New Panel...", self)
        new_panel.setObjectName("NewPanelAction")
        new_panel.setToolTip("Pick a panel to open in this window")
        new_panel.setShortcut(QtGui.QKeySequence("Ctrl+Shift+P"))
        new_panel.setShortcutContext(QtCore.Qt.ShortcutContext.WindowShortcut)
        new_panel.setMenuRole(QtGui.QAction.MenuRole.NoRole)
        new_panel.triggered.connect(lambda _c=False: self.app_controller.commands.run(PANEL_OPEN, self))
        self.addAction(new_panel)
        self._new_panel_action = new_panel
        # Window > Show Filters: the rail is fixed and shown by default; this
        # hides it for a wider docking area and brings it back (no shortcut).
        show_filters = QtGui.QAction("Show Filters", self)
        show_filters.setObjectName("ShowFiltersAction")
        show_filters.setToolTip("Show or hide the sample and annotation filters on the left")
        show_filters.setCheckable(True)
        show_filters.setChecked(True)
        show_filters.toggled.connect(self._set_filters_panel_visible)
        self.filters_action = show_filters
        # Window > Layouts (plan 3.8): presets and user layouts of the
        # docking area.
        layouts_menu = QtWidgets.QMenu("Layouts", menu)
        layouts_menu.setObjectName("LayoutsMenu")
        # The descriptions, and why a preset is disabled, are tooltips.
        layouts_menu.setToolTipsVisible(True)
        layouts_menu.aboutToShow.connect(self._fill_layouts_menu)
        self._layouts_menu: QtWidgets.QMenu = layouts_menu
        # Kept for the window's life (clear() would leave a new submenu
        # behind on every fill).
        delete_menu = QtWidgets.QMenu("Delete Saved Layout", layouts_menu)
        delete_menu.setObjectName("DeleteLayoutMenu")
        self._delete_layout_menu = delete_menu
        # Filled each time it opens (and once when the window has registered
        # with the WindowManager, in ExploreWindow.__init__).
        menu.aboutToShow.connect(
            lambda: self.app_controller.windows.fill_window_menu(
                menu, self, self._layouts_menu, new_panel_action=new_panel, filters_action=show_filters
            )
        )
        self._window_menu = menu

    def _fill_layouts_menu(self) -> None:
        """Window > Layouts: the built-in presets, the saved user layouts,
        then Save Layout As..., the default for new windows and deleting a
        saved layout. Every entry runs a layout command on this window; none
        has a shortcut (D8)."""
        menu = self._layouts_menu
        manager = getattr(self, "_panel_manager", None)
        if menu is None or manager is None:
            return
        commands = self.app_controller.commands
        store = self.app_controller.layouts
        menu.clear()
        add_menu_header(menu, "Built-in layouts")
        for preset in dock_layouts.PRESETS:
            action = menu.addAction(preset.title)
            action.setObjectName(f"LayoutPreset_{preset.key}")
            available, reason = manager.preset_available(preset.key)
            action.setEnabled(available)
            action.setToolTip(preset.description if available else reason)
            action.triggered.connect(
                lambda _c=False, _k=preset.key: commands.run(LAYOUT_APPLY, self, preset=_k)
            )
        names = store.layout_names()
        menu.addSeparator()
        add_menu_header(menu, "Saved layouts")
        if not names:
            empty = menu.addAction("No saved layouts")
            empty.setObjectName("NoSavedLayouts")
            empty.setEnabled(False)
        for name in names:
            # "&" would make a mnemonic; the name itself is in the data.
            action = menu.addAction(name.replace("&", "&&"))
            action.setObjectName("UserLayoutAction")
            action.setData(name)
            action.setToolTip("Arrange the panels as saved; filters and the sample set stay as they are")
            action.triggered.connect(lambda _c=False, _n=name: commands.run(LAYOUT_APPLY, self, layout=_n))
        menu.addSeparator()
        save = menu.addAction("Save Layout As...")
        save.setObjectName("SaveLayoutAsAction")
        save.setToolTip("Save the arrangement of panels under a name (no filters, no rail state)")
        save.triggered.connect(lambda _c=False: commands.run(LAYOUT_SAVE, self))
        set_default = menu.addAction("Set as Default for New Windows")
        set_default.setObjectName("SetDefaultLayoutAction")
        set_default.setToolTip(
            "New explore windows, and runs opened without a saved layout, start with this arrangement of panels"
        )
        set_default.triggered.connect(lambda _c=False: commands.run(LAYOUT_SET_DEFAULT, self))
        clear_default = menu.addAction("Use Classic for New Windows")
        clear_default.setObjectName("ClearDefaultLayoutAction")
        clear_default.setToolTip("New explore windows start with the Classic layout again")
        clear_default.setEnabled(store.has_default())
        clear_default.triggered.connect(lambda _c=False: self._clear_default_layout())
        delete_menu = self._delete_layout_menu
        delete_menu.clear()
        menu.addMenu(delete_menu)
        delete_menu.setEnabled(bool(names))
        for name in names:
            action = delete_menu.addAction(name.replace("&", "&&"))
            action.setData(name)
            action.triggered.connect(lambda _c=False, _n=name: self._delete_layout_interactive(_n))

    def _start_heavy_panels(self) -> None:
        """Heavy panels build lazily on first show; the first ones wait until
        shortly after the main view is ready."""
        panels = getattr(self, "_panels", None)
        if panels is not None:
            panels.set_ready()

    def _show_panel_palette(self) -> None:
        """Window > New Panel... (Ctrl+Shift+P): the panel picker, over this
        window."""
        manager = getattr(self, "_panel_manager", None)
        if manager is not None:
            manager.show_palette()

    def _build_sample_tree_panel(self) -> QtWidgets.QWidget:
        """The sample tree view (the Sample tree panel's content)."""
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
            + list(self.run.annotation_db_columns)
            + ["ST", "Quality", "Complete", "Contam."]
        )
        self.tree.setHeaderLabels(tree_headers)
        self.tree.setHeaderHidden(False)
        # The four fixed QC columns are appended after the existing (optional)
        # columns so the db-hit heat-cell indices below are unchanged.
        _qc_base = 5 + len(self.meta_columns_for_tree) + 1 + len(self.run.annotation_db_columns)
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
        self.tree.setHorizontalScrollMode(
            QtWidgets.QAbstractItemView.ScrollMode.ScrollPerPixel
        )
        self.tree.setTextElideMode(QtCore.Qt.TextElideMode.ElideMiddle)
        self.tree.setSortingEnabled(True)
        header = self.tree.header()
        # Keep the biological identifier readable even when QC columns require
        # horizontal scrolling. A stretching Name section could collapse to a
        # few pixels in a narrow tree pane.
        header.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.Interactive)
        self.tree.setColumnWidth(0, 220)
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
        self.tree.setStyleSheet(SHARED_ITEM_VIEW_QSS)
        for col_idx in range(1, len(tree_headers)):
            self.tree.setColumnHidden(col_idx, col_idx not in self._fixed_tree_columns)
        tree_layout.addWidget(self.tree, stretch=1)
        return tree_panel

    def _new_annotation_table(self, *, paired: bool) -> QtWidgets.QTableView:
        """An annotation table view (without its "Table shows" switch). The
        paired one (the first Annotations panel's) is the window's ``table``
        and drives the genome track when scrolled."""
        table = QtWidgets.QTableView()
        table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        table.setAlternatingRowColors(True)
        # Same frame + header as the sample tree. Keep the grid on (not `::item`
        # borders) so the model's group / separator / hit-row backgrounds still
        # paint; the light gridline-color supplies the column separators.
        table.setStyleSheet(_FRAME_HEADER_QSS)
        table.setShowGrid(True)
        table.setMouseTracking(True)
        table.setItemDelegate(RowRuleDelegate(table))
        # Gap rows ("N genes not shown") are shorter than data rows.
        table.verticalHeader().setMinimumSectionSize(GAP_ROW_HEIGHT)
        table.viewport().setMouseTracking(True)
        header = table.horizontalHeader()
        header.setStretchLastSection(True)
        header.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Stretch)
        table.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        table.setHorizontalScrollMode(
            QtWidgets.QAbstractItemView.ScrollMode.ScrollPerPixel
        )
        if paired:
            self.table = table
            table.verticalScrollBar().rangeChanged.connect(self._on_table_scroll_range_changed)
            table.verticalScrollBar().valueChanged.connect(self._on_table_scroll_changed)
        return table

    def _build_track_widget(self) -> QtWidgets.QWidget:
        """The genome track (the Genome track panel's content)."""
        table_placeholder = QtWidgets.QWidget()
        placeholder_layout = QtWidgets.QVBoxLayout(table_placeholder)
        placeholder_layout.setContentsMargins(0, 0, 0, 0)
        placeholder_layout.setSpacing(0)

        self.igv_view = IgvTrackView(table_placeholder)
        self.igv_view.geneDoubleClicked.connect(self._on_igv_gene_double_clicked)
        self.igv_view.viewCenterBpChanged.connect(self._on_igv_center_bp_changed)
        placeholder_layout.addWidget(self.igv_view, stretch=1)
        return table_placeholder

    def _create_annotations_panel(self, uid: str, parent: Optional[QtWidgets.QWidget] = None) -> ExplorePanel:
        """Factory of the Annotations kind: a table with its own "Table shows"
        box. The first instance is paired with the genome track (the
        window's ``table`` and ``annotations_panel``, whose view mode the
        track follows)."""
        strip, box = build_view_mode_switch()
        options = build_window_options(strip)
        paired = uid == ANNOTATIONS
        table = self._new_annotation_table(paired=paired)
        panel = AnnotationsPanel(
            strip, table, parent, window=self, table=table, mode_box=box, window_options=options
        )
        panel.statusChanged.connect(self._update_active_panel_status)
        panel.windowOptionToggled.connect(self._set_window_option)
        if paired:
            self.annotations_panel = panel
            panel.viewModeChanged.connect(self._on_annotations_view_mode_changed)
        else:
            panel.viewModeChanged.connect(lambda _mode, _panel=panel: self._on_extra_annotations_mode_changed(_panel))
        return panel

    def _create_sample_tree_panel(self, uid: str, parent: Optional[QtWidgets.QWidget] = None) -> ExplorePanel:
        """Factory of the Sample tree kind; its "Hide contigs without hits"
        checkbox writes the window's draft (see _set_window_option)."""
        del uid
        panel = SampleTreePanel(self, parent)
        panel.windowOptionToggled.connect(self._set_window_option)
        return panel

    def _set_window_option(self, field: str, value: bool) -> None:
        """A window filter checkbox was toggled ("Hide hypothetical" in an
        Annotations panel, "Hide contigs without hits" in the Sample tree):
        write it into the draft, with the rest of the rail as it is now, and
        apply at once. Every panel showing it follows through draftChanged."""
        if field not in {"hide_hypothetical", "hide_empty_contigs"}:
            return
        self.ctx.update_draft(FilterState.from_widgets(self).replace(**{field: bool(value)}))

    def _annotations_panel_open(self) -> bool:
        manager = getattr(self, "_panel_manager", None)
        if manager is None:
            return False
        return any(manager.is_open(uid) for uid in manager.uids(ANNOTATIONS))

    def _update_window_option_widgets(self, *_args: object) -> None:
        """Show a rail note while no Annotations panel is open and "Hide
        hypothetical" (set in that panel) is on. "Hide contigs without hits"
        sits on the Sample tree, the only panel it changes, so it needs none."""
        manager = getattr(self, "_panel_manager", None)
        if manager is None:
            return
        draft = getattr(self.ctx, "draft", None)
        hide_hypo = bool(getattr(draft, "hide_hypothetical", False))
        note = getattr(self, "window_options_note", None)
        if note is None:
            return
        changed: List[str] = []
        if hide_hypo:
            changed.append("Hide hypothetical is on")
        show = bool(changed) and not self._annotations_panel_open()
        if show:
            note.setText(
                "; ".join(changed)
                + " (an Annotations panel option for the whole window). "
                + '<a href="annotations">Open Annotations</a>'
            )
            note.setToolTip("This window filter is set in the Annotations panel, which is closed")
        note.setVisible(show)

    def _open_annotations_for_options(self, _link: str = "") -> None:
        manager = getattr(self, "_panel_manager", None)
        if manager is not None:
            manager.focus_panel(ANNOTATIONS)

    def _on_extra_annotations_mode_changed(self, panel: AnnotationsPanel) -> None:
        self._render_annotations(panel)
        self._update_active_panel_status()

    # ------------------------------------------------------------ docking area
    def _explore_panel_specs(self) -> List[PanelSpec]:
        """Panel kinds of this window, in palette order."""

        def needs_run_output() -> Tuple[bool, str]:
            if self.run.output_dir is None:
                return False, "Not available: this run has no output folder."
            return True, ""

        def phylogeny_detail(state: Dict[str, Any]) -> str:
            # The tree a Phylogeny panel not built yet will draw: its saved
            # source, else the run's default tree.
            if self.run.output_dir is None:
                return ""
            index = self.run.pan.tree_index()
            source = index.source(state.get("source")) or index.default()
            return source.short_label if source is not None else ""

        def pan_spec(panel_class: type) -> PanelSpec:
            return PanelSpec(
                panel_class.KIND,
                panel_class.TITLE,
                panel_class.GROUP,
                heavy=True,
                placeholder=f"{panel_class.TITLE} will load when this panel is shown.",
                available=functools.partial(self._pan_view_available, panel_class.VIEW),
                single_instance=False,
                factory=self._create_pan_view_panel,
                state_detail=phylogeny_detail if panel_class is PhylogenyPanel else None,
            )

        return [
            PanelSpec(SAMPLE_TREE, "Sample tree", "Samples", factory=self._create_sample_tree_panel),
            PanelSpec(
                ANNOTATIONS,
                "Annotations",
                "Annotations",
                single_instance=False,
                factory=self._create_annotations_panel,
            ),
            # The track draws what the first Annotations panel shows: it opens
            # and closes with it (PanelSpec.requires).
            PanelSpec(
                TRACK,
                "Genome track",
                "Annotations",
                factory=lambda _uid, parent: TrackPanel(self, parent),
                requires=ANNOTATIONS,
                min_height=TRACK_MIN_HEIGHT,
            ),
            PanelSpec(
                GEO,
                "Map",
                "Geography",
                heavy=True,
                placeholder="Geo map will load after the main view is ready.",
                single_instance=False,
                factory=self._create_geo_panel,
            ),
            # Pan panels (plan 8.10, Phase 9): one pan view each, disabled
            # with the reason when the run lacks its data. The Gene
            # presence/absence rows of a window follow its Phylogeny's leaf
            # order (panels/row_link.py). The compact cgMLST MST follows the
            # cgMLST MST below.
            *[pan_spec(panel_class) for panel_class in PAN_VIEW_PANELS if panel_class is not PanForcePanel],
            PanelSpec(
                PanClusteringPanel.KIND,
                PanClusteringPanel.TITLE,
                PanClusteringPanel.GROUP,
                heavy=True,
                placeholder=f"{PanClusteringPanel.TITLE} will load when this panel is shown.",
                available=functools.partial(self._pan_view_available, VIEW_PRESENCE),
                single_instance=False,
                factory=self._create_pan_mds_panel,
                state_detail=PanClusteringPanel.state_detail,
            ),
            PanelSpec(
                CGMLST_FORCE,
                CgmlstForcePanel.TITLE,
                "Typing",
                heavy=True,
                placeholder=f"The {CgmlstForcePanel.TITLE} will load when this panel is shown.",
                available=needs_run_output,
                single_instance=False,
                factory=self._create_cgmlst_force_panel,
                state_detail=CgmlstForcePanel.state_detail,
            ),
            pan_spec(PanForcePanel),
            # Document panels (plan 8.7): opened with their data from a gene
            # card or a pan view, one per comparison or alignment.
            PanelSpec(
                NEIGHBOURHOOD,
                "Gene neighborhoods",
                "Context",
                heavy=True,
                single_instance=False,
                document=True,
                factory=self._create_neighbourhood_panel,
            ),
            PanelSpec(
                ALIGNMENT,
                "Sequence alignment",
                "Context",
                heavy=True,
                single_instance=False,
                document=True,
                factory=self._create_alignment_panel,
            ),
        ]

    def _build_docking_area(self) -> QtWidgets.QWidget:
        """The docking area right of the rail: every view is a panel."""
        from .panels.manager import PanelManager  # imports dock_host (PyQt6Ads)

        manager = PanelManager(self, self._explore_panel_specs(), ctx=self.ctx)
        self._panel_manager = manager
        self._panels = manager
        # Phylogeny and Gene presence/absence rows of this window, linked.
        self._row_link = TreeRowLink(manager, self)
        manager.activePanelChanged.connect(self._on_active_panel_changed)
        # The rail note and the strip follow the draft and which panels are open.
        manager.host.viewOpened.connect(self._update_window_option_widgets)
        manager.host.viewClosed.connect(self._update_window_option_widgets)
        manager.instancesChanged.connect(self._update_window_option_widgets)
        manager.layoutApplied.connect(self._update_window_option_widgets)
        # Placing a panel that a layout left out emits no viewOpened, but
        # opening it focuses it.
        manager.activePanelChanged.connect(self._update_window_option_widgets)
        self.ctx.draftChanged.connect(self._update_window_option_widgets)
        return manager.host

    def _on_active_panel_changed(self, _uid: str) -> None:
        # Changing the active panel clears hover, as a bottom-tab change did,
        # and changes which panel line the status bar appends (plan 3.9).
        self.ctx.clear_hover()
        self._update_active_panel_status()

    def _active_annotations_panel(self) -> Optional[AnnotationsPanel]:
        manager = getattr(self, "_panel_manager", None)
        if manager is None:
            return None
        uid = manager.host.active_view()
        panel = manager.panel(uid) if uid else None
        return panel if isinstance(panel, AnnotationsPanel) else None

    def _update_active_panel_status(self, *_args: object) -> None:
        """Append the active Annotations panel's own line to the status bar."""
        label = getattr(self, "active_panel_status", None)
        if label is None:
            return
        panel = self._active_annotations_panel()
        text = panel.status_text() if panel is not None else ""
        label.setText(f"| {text}" if text else "")
        label.setToolTip(text)
        label.setVisible(bool(text))

    def _set_filters_panel_visible(self, visible: bool) -> None:
        panel = getattr(self, "_filters_panel", None)
        splitter = getattr(self, "_workspace_splitter", None)
        toggle = getattr(self, "filters_action", None)
        if panel is None or splitter is None:
            return

        if toggle is not None and toggle.isChecked() != bool(visible):
            blocker = QtCore.QSignalBlocker(toggle)
            toggle.setChecked(bool(visible))
            del blocker

        sizes = splitter.sizes()
        total = max(splitter.width(), sum(sizes), 800)
        self._sync_filters_strip(bool(visible))
        if visible:
            panel.setMinimumWidth(getattr(self, "_filters_min_width", RAIL_MIN_WIDTH))
            panel.setMaximumWidth(RAIL_MAX_WIDTH)
            panel.show()
            if not sizes or sizes[0] <= 1:
                rail_width = clamp_rail_width(getattr(self, "_filters_last_width", RAIL_DEFAULT_WIDTH))
                splitter.setSizes([rail_width, max(1, total - rail_width)])
        else:
            if sizes and sizes[0] > 1:
                self._filters_last_width = clamp_rail_width(sizes[0])
            panel.setMinimumWidth(0)
            panel.hide()

    def _on_workspace_splitter_moved(self, _position: int, _index: int) -> None:
        splitter = getattr(self, "_workspace_splitter", None)
        toggle = getattr(self, "filters_action", None)
        if splitter is None or toggle is None:
            return
        sizes = splitter.sizes()
        expanded = bool(sizes and sizes[0] > 1)
        if expanded:
            self._filters_last_width = clamp_rail_width(sizes[0])
        self._sync_filters_strip(expanded)
        if toggle.isChecked() != expanded:
            blocker = QtCore.QSignalBlocker(toggle)
            toggle.setChecked(expanded)
            del blocker

    def _sync_filters_strip(self, rail_shown: bool) -> None:
        """The folded rail's strip shows exactly while Show Filters is off
        (hidden from the menu or the chevron, or dragged shut)."""
        strip = getattr(self, "_filters_strip", None)
        if strip is not None and strip.isHidden() != rail_shown:
            strip.setVisible(not rail_shown)

    # ------------------------------------------------- layout persistence v2
    def _rail_state(self) -> Tuple[int, bool]:
        """(width, collapsed) of the filter rail, as saved per window."""
        splitter = getattr(self, "_workspace_splitter", None)
        panel = getattr(self, "_filters_panel", None)
        sizes = splitter.sizes() if splitter is not None else []
        if sum(sizes) <= 0:
            # Not laid out yet (never shown): only the hidden flag is known.
            visible = bool(panel is not None and not panel.isHidden())
        else:
            visible = bool(panel is not None and not panel.isHidden() and sizes[0] > 1)
        pending = getattr(self, "_pending_rail_width", None)
        if pending is not None:
            width = pending
        elif visible and self.isVisible() and sizes and sizes[0] > 1:
            width = sizes[0]
        else:
            width = getattr(self, "_filters_last_width", RAIL_DEFAULT_WIDTH)
        return clamp_rail_width(width), not visible

    def _apply_rail_state(self, width: int, collapsed: bool) -> None:
        width = clamp_rail_width(width)
        self._set_filters_panel_visible(not collapsed)
        # After the visibility change, which records the splitter's current
        # (not yet laid out) width as the last width.
        self._filters_last_width = width
        if not collapsed:
            # The splitter has no real width until the window is laid out.
            self._pending_rail_width = width
            self._apply_pending_rail_width()

    def _apply_pending_rail_width(self, attempt: int = 0) -> None:
        width = getattr(self, "_pending_rail_width", None)
        splitter = getattr(self, "_workspace_splitter", None)
        if width is None or splitter is None:
            return
        total = sum(splitter.sizes())
        if not self.isVisible() or total < width + 200:
            if self.isVisible() and attempt < 40:
                QtCore.QTimer.singleShot(25, lambda: self._apply_pending_rail_width(attempt + 1))
            return
        self._pending_rail_width = None
        splitter.setSizes([width, max(1, total - width)])
        if getattr(self, "_pending_dock_state", None) is not None:
            # The docking area has its final width only after this layout pass.
            QtCore.QTimer.singleShot(0, self._apply_pending_dock_state)

    def _apply_pending_dock_state(self) -> None:
        """Duplicate Window: put the copied dock state back once the window
        is laid out at its size. Restored while the window was still unsized,
        QSplitter rescaled its sizes and a splitter could end up 1 px off the
        source window's (and drift further with each copy of a copy).
        ``_pending_dock_source``: the window copied, set by the caller."""
        state = getattr(self, "_pending_dock_state", None)
        self._pending_dock_state = None
        manager = getattr(self, "_panel_manager", None)
        if not state or manager is None or not self.isVisible():
            return
        if getattr(self, "_pending_rail_width", None) is not None:
            # Still waiting for the rail: it calls back here when done.
            self._pending_dock_state = state
            return
        source, self._pending_dock_source = getattr(self, "_pending_dock_source", None), None
        uids = getattr(self, "_pending_dock_uids", None)
        if sorted(manager.host.open_uids()) != uids:
            return  # the user changed the arrangement meanwhile
        # The source window lost the height of the window color stripe when
        # this window opened (its splitters rescaled): copy its current state
        # when it still holds the same panels, so both round the same way.
        source_manager = getattr(source, "_panel_manager", None) if source is not None else None
        try:
            if source_manager is not None and sorted(source_manager.host.open_uids()) == uids:
                state = source_manager.host.save_layout()
        except RuntimeError:  # the source window is gone
            pass
        manager.host.restore_layout(state)

    def showEvent(self, event: QtGui.QShowEvent) -> None:
        super().showEvent(event)
        if getattr(self, "_pending_rail_width", None) is not None:
            QtCore.QTimer.singleShot(0, self._apply_pending_rail_width)
        elif getattr(self, "_pending_dock_state", None) is not None:
            QtCore.QTimer.singleShot(0, self._apply_pending_dock_state)

    def _restore_explore_layout(self) -> None:
        """The layout this window starts with (``self._start_layout``):

        - a slot number (the default: the window's own number): the layout
          saved in that slot of this run; for slot 1 with nothing saved, the
          last saved window's options; else Classic;
        - "default": the default layout with a clean rail (New Explore
          Window, a window opened from a clade and the like);
        - "empty": an empty docking area (the neighbourhood window, whose
          comparison panel opens into it once the window is built);
        - a ``WindowLayout``: a copy of another window's layout (Duplicate
          Window).

        The default layout is the arrangement set with "Set as Default for
        New Windows" (plan 3.8), else Classic. A run with no saved layout
        also starts with it (plan 3.11), taking the window size and rail
        state of the window saved last.

        Layout only: filters and searches start fresh (D5)."""
        manager = self._panel_manager
        store = self.app_controller.layouts
        self._layout_store = store
        self._layout_run_key = run_key(self.run.output_dir)
        start = getattr(self, "_start_layout", None)
        restored = False
        if isinstance(start, WindowLayout):
            if start.geometry:
                self._geometry_restored = bool(self.restoreGeometry(QtCore.QByteArray(start.geometry)))
            self._apply_rail_state(start.rail_width, start.rail_collapsed)
            if start.version == FORMAT_VERSION:
                restored = manager.restore_layout(start.dock_state, start.panels)
            if restored and not start.geometry:
                # A copy opens at the source window's size only once shown:
                # restore the dock state again then (see below).
                self._pending_dock_state = start.dock_state
                self._pending_dock_uids = sorted(manager.host.open_uids())
            self._finish_layout_restore(restored)
            return
        if start == START_DEFAULT_LAYOUT:
            self._set_filters_panel_visible(True)
            self._finish_layout_restore(self._start_with_default_arrangement())
            return
        if start == START_EMPTY_LAYOUT:
            # Nothing is built: the caller opens the window's first panel.
            self._set_filters_panel_visible(True)
            self._finish_layout_restore(True)
            return
        slot = int(start) if isinstance(start, int) and start > 0 else self._window_number
        saved = store.load_window(self._layout_run_key, slot)
        if saved is not None:
            if saved.geometry:
                self._geometry_restored = bool(self.restoreGeometry(QtCore.QByteArray(saved.geometry)))
            self._apply_rail_state(saved.rail_width, saved.rail_collapsed)
            if saved.version == FORMAT_VERSION:
                restored = manager.restore_layout(saved.dock_state, saved.panels)
        elif slot != 1:
            # Only the first window inherits the older fallbacks below.
            self._set_filters_panel_visible(True)
            restored = self._start_with_default_arrangement()
        else:
            # A run opened for the first time keeps the window size, the rail
            # and the panel options of the window saved last; the dock layout
            # itself is the default arrangement (Classic unless one is set).
            last = store.load_last_window()
            if last is not None and last.version == FORMAT_VERSION:
                if last.geometry:
                    self._geometry_restored = bool(self.restoreGeometry(QtCore.QByteArray(last.geometry)))
                self._apply_rail_state(last.rail_width, last.rail_collapsed)
                restored = self._start_with_default_arrangement()
                if not restored:
                    manager.apply_panel_list(last.panels, create_instances=False)
            else:
                self._set_filters_panel_visible(True)
                restored = self._start_with_default_arrangement()
        self._finish_layout_restore(restored)

    def _start_with_default_arrangement(self) -> bool:
        """Arrange the empty docking area as the default for new windows,
        if one is set and usable here; False leaves it empty (Classic is
        built next)."""
        default = self._layout_store.load_default()
        if default is None or default.version != FORMAT_VERSION or not default.dock_state:
            return False
        return self._panel_manager.apply_saved_layout(default.dock_state, default.panels)

    def _finish_layout_restore(self, restored: bool) -> None:
        manager = self._panel_manager
        if not restored:
            manager.build_default()
        self._last_saved_layout = self._explore_layout_snapshot()
        timer = QtCore.QTimer(self)
        timer.setInterval(LAYOUT_AUTOSAVE_MS)
        timer.timeout.connect(self._autosave_explore_layout)
        timer.start()
        self._layout_autosave_timer = timer

    def _explore_layout_snapshot(self) -> WindowLayout:
        manager = self._panel_manager
        width, collapsed = self._rail_state()
        return WindowLayout(
            geometry=bytes(self.saveGeometry()),
            rail_width=width,
            rail_collapsed=collapsed,
            dock_state=manager.host.save_layout(),
            panels=manager.panel_list(),
        )

    def _save_explore_layout(self) -> None:
        """Save into the slot of this window's number (plan 3.11)."""
        manager = getattr(self, "_panel_manager", None)
        store = getattr(self, "_layout_store", None)
        if manager is None or store is None:
            return
        if getattr(self, "_neighbourhood_window", False):
            # A neighbourhood window is not reopened and leaves no saved
            # state: its number may be 1, and its size and panels (a
            # comparison that cannot be restored) must not replace window 1's
            # layout or seed a run opened for the first time.
            return
        snapshot = self._explore_layout_snapshot()
        store.save_window(self._layout_run_key, snapshot, self.window_number)
        store.save_last_window(snapshot)
        self._last_saved_layout = snapshot

    def _layout_copy(self) -> Optional[WindowLayout]:
        """This window's layout for Duplicate Window (docking area, panel list
        with display options, rail state; no geometry: the copy opens offset
        from this window). None before the docking area exists."""
        if getattr(self, "_panel_manager", None) is None:
            return None
        snapshot = self._explore_layout_snapshot()
        snapshot.geometry = b""
        return snapshot

    # ------------------------------------------------- presets and user layouts
    def _arrangement_snapshot(self) -> Optional[WindowLayout]:
        """The docking area alone (dock state and panel list): what a user
        layout and the default for new windows keep (plan 3.8)."""
        manager = getattr(self, "_panel_manager", None)
        if manager is None:
            return None
        return WindowLayout(dock_state=manager.host.save_layout(), panels=manager.panel_list())

    def _layout_status(self, text: str) -> None:
        self.layout_message.setText(text)
        self.layout_message.setVisible(True)
        self._layout_message_timer.start()

    def _has_arrangement(self) -> bool:
        """False while no panel is open (the empty state): such an
        arrangement cannot be applied, so it is not saved."""
        manager = getattr(self, "_panel_manager", None)
        if manager is None:
            return False
        if manager.has_open_panels():
            return True
        self._layout_status("Open a panel first; an empty layout cannot be saved")
        return False

    def _save_layout_as(self, name: Optional[str] = None) -> Optional[str]:
        """Save this window's arrangement of panels as a user layout. Without
        ``name`` the user is asked for one (and whether to replace a layout
        of that name); with it, a layout of that name is replaced."""
        store = self.app_controller.layouts
        if not self._has_arrangement():
            return None
        if name is None:
            text, ok = QtWidgets.QInputDialog.getText(
                self,
                "Save Layout As",
                "Name of the layout (filters and the rail are not part of it):",
            )
            if not ok:
                return None
            name = store.clean_layout_name(text)
            if not name:
                return None
            existing = store.existing_layout_name(name)
            if existing is not None:
                answer = QtWidgets.QMessageBox.question(
                    self,
                    "Replace layout",
                    f'A layout named "{existing}" exists. Replace it?',
                )
                if answer != QtWidgets.QMessageBox.StandardButton.Yes:
                    return None
        name = store.clean_layout_name(name)
        snapshot = self._arrangement_snapshot()
        if not name or snapshot is None:
            return None
        saved = store.save_user_layout(name, snapshot)
        self._layout_status(f'Layout "{saved}" saved')
        return saved

    def _apply_layout(self, *, preset: Optional[str] = None, name: Optional[str] = None) -> bool:
        """Arrange the docking area as a built-in preset or a saved user
        layout. Never touches the rail, the filters or the sample set."""
        manager = getattr(self, "_panel_manager", None)
        if manager is None:
            return False
        if preset is not None:
            return manager.apply_preset(preset)
        store = self.app_controller.layouts
        saved = store.load_user_layout(name or "")
        if saved is None:
            self._layout_status(f'No saved layout is named "{name}"')
            return False
        if saved.version != FORMAT_VERSION or not manager.apply_saved_layout(saved.dock_state, saved.panels):
            self._layout_status(f'Layout "{name}" could not be applied; the panels were left as they were')
            return False
        return True

    def _set_default_layout(self) -> bool:
        if not self._has_arrangement():
            return False
        snapshot = self._arrangement_snapshot()
        if snapshot is None:
            return False
        self.app_controller.layouts.save_default(snapshot)
        self._layout_status("New explore windows will start with this layout")
        return True

    def _clear_default_layout(self) -> None:
        self.app_controller.commands.run(LAYOUT_CLEAR_DEFAULT, self)
        self._layout_status("New explore windows will start with the Classic layout")

    def _delete_layout_interactive(self, name: str) -> bool:
        answer = QtWidgets.QMessageBox.question(
            self,
            "Delete layout",
            f'Delete the saved layout "{name}"? Windows using it keep their panels.',
        )
        if answer != QtWidgets.QMessageBox.StandardButton.Yes:
            return False
        return bool(self.app_controller.commands.run(LAYOUT_DELETE, self, layout=name))

    def _autosave_explore_layout(self) -> None:
        """Every 60 s, save if the layout changed (never while a group is maximized)."""
        manager = getattr(self, "_panel_manager", None)
        if manager is None or manager.host.is_maximized():
            return
        snapshot = self._explore_layout_snapshot()
        if snapshot != getattr(self, "_last_saved_layout", None):
            self._save_explore_layout()

    def _hide_layout_message(self) -> None:
        self.layout_message.setVisible(False)

    def _show_busy_status(self) -> None:
        self.busy_status.setVisible(self.ctx.is_busy())

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        self._save_explore_layout()
        super().closeEvent(event)
        if event.isAccepted():
            # A closed window must not autosave over the run's layout slot,
            # nor show "Updating..." for a compute still finishing.
            for name in ("_layout_autosave_timer", "_busy_timer"):
                timer = getattr(self, name, None)
                if timer is not None:
                    timer.stop()
            closed = getattr(self, "windowClosed", None)
            if closed is not None:
                closed.emit()
