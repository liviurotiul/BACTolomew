from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets

from ..igv_view import IgvTrackView
from ..explorer_ui import ExplorerUIMixin
from ..filters_controller import FiltersControllerMixin
from ..tree_controller import TreeControllerMixin
from ..table_controller import TableControllerMixin
from ..igv_controller import IgvControllerMixin
from ..panels_controller import PanelsControllerMixin
from ..docking import require_docking
from ..model import QueryEngine, RunData
from ..model.filter_state import AppliedResult, FilterState, SampleSubset
from ..panels.annotations import AnnotationsPanel
from ..panels.manager import PanelManager
from ..session import bus as bus_events
from ..session.explore_context import ExploreContext
from ..shell.app_controller import AppController
from .filter_rail import FilterRail

# "Hide contigs without hits" in a new window (Sample tree panel checkbox).
DEFAULT_HIDE_EMPTY_CONTIGS = True


class ExploreWindow(
    ExplorerUIMixin,
    FiltersControllerMixin,
    TreeControllerMixin,
    TableControllerMixin,
    IgvControllerMixin,
    PanelsControllerMixin,
    QtWidgets.QMainWindow,
):
    # Emitted once, when the first refresh has populated the sample tree; used
    # to dismiss the loading splash exactly when the view is actually ready.
    firstRefreshFinished = QtCore.pyqtSignal()
    # Emitted when the window has closed (its layout is saved by then); the
    # WindowManager rebinds or closes its cards and applies the quit policy.
    windowClosed = QtCore.pyqtSignal()

    def __init__(
        self,
        ann_df: Optional[pd.DataFrame] = None,
        meta_df: Optional[pd.DataFrame] = None,
        ref_df: Optional[pd.DataFrame] = None,
        marker_hits_df: Optional[pd.DataFrame] = None,
        output_dir: Optional[Path] = None,
        *,
        run_data: Optional[RunData] = None,
        engine: Optional[QueryEngine] = None,
        app_controller: Optional[AppController] = None,
        layout: Any = None,
        initial_filters: Optional[FilterState] = None,
        window_subset: Optional[SampleSubset] = None,
    ) -> None:
        """``layout``: what the docking area starts with (None: the layout
        saved for this window's number; ``START_EMPTY_LAYOUT``: an empty
        docking area (the neighbourhood window); ``START_DEFAULT_LAYOUT``: Classic with
        a clean rail; a slot number; or a ``WindowLayout`` to copy).
        ``initial_filters``: the rail starts with these filters instead of a
        clean rail (Duplicate Window). ``window_subset``: the locked subset of
        a window opened from a clade, cluster or plasmid (plan 3.10).

        Raises ``DockingUnavailableError`` (with the install hint) when
        PyQt6Ads is missing: the docking area is the only layout."""
        require_docking()
        super().__init__()
        # The run's data is either prepared already (the app builds it in its
        # loader thread) or built here from the tables (tests, scripts).
        if run_data is None:
            if ann_df is None or meta_df is None:
                raise ValueError("ann_df and meta_df are required when no run_data is given.")
            run_data = RunData.build(
                ann_df,
                meta_df,
                ref_df=ref_df,
                marker_hits_df=marker_hits_df,
                output_dir=output_dir,
            )
        elif any(value is not None for value in (ann_df, meta_df, ref_df, marker_hits_df, output_dir)):
            raise ValueError("Pass either the run's tables or run_data, not both.")
        if engine is not None and engine.run is not run_data:
            raise ValueError("The engine was built for another run's RunData.")
        if app_controller is not None:
            if app_controller.run is not run_data:
                raise ValueError("The app controller belongs to another run's RunData.")
            if engine is not None and engine is not app_controller.engine:
                raise ValueError("Pass the app controller's own engine (or none).")
            engine = app_controller.engine
        self.run: RunData = run_data
        self.engine: QueryEngine = engine if engine is not None else QueryEngine(run_data)
        # The application object (plan 4.2): the window manager that owns this
        # window's info cards, the named commands and the layout store. A
        # window built on its own (tests, scripts) gets a private one.
        self.app_controller: AppController = (
            app_controller if app_controller is not None else AppController(run_data, self.engine, parent=self)
        )
        # This window's number (it registers under it once built) and how it
        # starts (plan 3.10, 3.11).
        self._window_number = self.app_controller.windows.next_number()
        self._several_windows = False
        self._start_layout = layout
        # A neighbourhood window (opened by Compare neighborhoods or a
        # plasmid context): its title names its comparison, and it is not
        # reopened with the run.
        self._neighbourhood_window = False
        if window_subset is None and initial_filters is not None:
            window_subset = initial_filters.window_subset
        # The locked subset of this window (FilterState.window_subset): set
        # once, never removed from the rail.
        self._window_subset: Optional[SampleSubset] = (
            window_subset if window_subset is not None and window_subset.samples else None
        )
        # show all metadata columns in the tree
        self.meta_columns_for_tree = list(self.run.meta_columns)
        self.metadata_filters: Dict[str, List[str]] = {}
        # The shared session of this window: draft and applied filters, the
        # applied result with the one sample set, the selection and the
        # hover. The window's chrome (status bar, rail) and every panel
        # follow it through its bus (plan 4.5).
        self.ctx = ExploreContext(self.engine, self)
        # A new window hides contigs without hits (the default the rail
        # checkbox used to give; the checkbox is now in the Sample tree
        # panel and the draft holds the value). initial_filters, when given,
        # replace it.
        self.ctx.draft = self.ctx.draft.replace(hide_empty_contigs=DEFAULT_HIDE_EMPTY_CONTIGS)
        self.ctx.applyFailed.connect(self._on_apply_failed)
        self.ctx.firstResultApplied.connect(self.firstRefreshFinished)
        self.ctx.openRequested.connect(self._on_open_request)
        bus = self.ctx.bus
        bus.subscribe("window", bus_events.RESULT_APPLIED, self._on_result_applied)
        bus.subscribe("window", bus_events.HOVER_CHANGED, self._show_hover_in_status)
        bus.subscribe("window", bus_events.BUSY_CHANGED, self._show_busy)
        self.filter_rail = FilterRail(self, self.ctx)
        self.annotations_panel: Optional[AnnotationsPanel] = None
        self._annotations_render_batch: Optional[Dict[int, Any]] = None
        self._sample_search_hint_query: Optional[str] = None
        self._tag_model: Optional[QtGui.QStandardItemModel] = None
        self._predictor_model: Optional[QtGui.QStandardItemModel] = None
        self.context_window = 30
        self.igv_view: Optional[IgvTrackView] = None
        self._igv_key: Optional[Tuple[Any, ...]] = None
        self._igv_dirty = False
        self._igv_bp_segments: List[Tuple[int, int, float, float]] = []
        self._igv_raw_contig_len = 0
        self._igv_display_contig_len = 0
        self._syncing_table_to_igv = False
        self._syncing_igv_to_table = False
        self._pending_igv_center_bp: Optional[int] = None
        self._table_row_bp_index: List[Tuple[int, int]] = []
        self._table_bp_values: List[int] = []
        self._workspace_splitter: Optional[QtWidgets.QSplitter] = None
        self._filters_panel: Optional[QtWidgets.QWidget] = None
        self._filters_strip: Optional[QtWidgets.QWidget] = None
        # The docking area and its panel instances (PanelManager); set up by
        # _build_ui.
        self._panel_manager: Optional[PanelManager] = None
        self._panels: Optional[PanelManager] = None
        self._geometry_restored = False
        self._sample_tree_items: Dict[str, QtWidgets.QTreeWidgetItem] = {}
        self._tree_hover_shown: Tuple[str, ...] = ()
        self._tree_bulk_toggle_in_progress = False
        # The paired table's scroll range just changed (see
        # _on_table_scroll_range_changed).
        self._table_scroll_range_changed = False
        self._table_scroll_timer = QtCore.QTimer(self)
        self._table_scroll_timer.setSingleShot(True)
        self._table_scroll_timer.setInterval(40)
        self._table_scroll_timer.timeout.connect(self._sync_igv_from_table_scroll)
        self._igv_scroll_timer = QtCore.QTimer(self)
        self._igv_scroll_timer.setSingleShot(True)
        self._igv_scroll_timer.setInterval(40)
        self._igv_scroll_timer.timeout.connect(self._sync_table_from_igv_center)
        self._sync_animation_ms = 200
        self._table_scroll_animation: Optional[QtCore.QPropertyAnimation] = None
        self._table_sync_token = 0
        self._igv_sync_token = 0

        self._build_ui()
        if initial_filters is not None:
            self._show_filters_in_rail(initial_filters.replace(window_subset=self._window_subset))
        self._render_metadata_filters()
        self.refresh()
        self.app_controller.windows.add_explore(self, self._window_number)
        window_menu = getattr(self, "_window_menu", None)
        if window_menu is not None:
            self.app_controller.windows.fill_window_menu(
                window_menu,
                self,
                getattr(self, "_layouts_menu", None),
                new_panel_action=getattr(self, "_new_panel_action", None),
                filters_action=getattr(self, "filters_action", None),
            )

    # Identity among the explore windows of the run (plan 3.10)

    @property
    def window_number(self) -> int:
        """This window's number (1, 2, ...)."""
        return self.app_controller.windows.number(self) or self._window_number

    @property
    def window_subset(self) -> Optional[SampleSubset]:
        """The locked subset of a window opened from a clade and the like."""
        return self._window_subset

    def set_window_identity(self, number: int, several: bool, color: str) -> None:
        """The WindowManager's numbering: ``[n]`` in the title and the color
        stripe, both only while several explore windows are open."""
        self._window_number = int(number) or self._window_number
        self._several_windows = bool(several)
        self._apply_window_title()
        stripe = getattr(self, "_window_stripe", None)
        if stripe is not None:
            stripe.setStyleSheet(f"QFrame#ExploreWindowStripe {{ background: {color}; border: none; }}")
            stripe.setVisible(bool(several and color))

    def _apply_window_title(self, *_args: object) -> None:
        """``BACTolomew - <run> [n] - <subset label> - <comparison>`` (plan
        3.10); the comparison only in a neighbourhood window, while one is open."""
        output_dir = self.run.output_dir
        run_name = Path(output_dir).name.strip() if output_dir is not None else ""
        title = f"BACTolomew - {run_name or 'Annotation Explorer'}"
        if getattr(self, "_several_windows", False):
            title += f" [{self.window_number}]"
        subset = getattr(self, "_window_subset", None)
        if subset is not None:
            title += f" - {subset.label or f'{len(subset.samples)} samples'}"
        document = self._title_document()
        if document:
            title += f" - {document}"
        self.setWindowTitle(title)

    def _title_document(self) -> str:
        """A neighbourhood window's comparison: the tab title of the first
        comparison opened in it that is still open, else ""."""
        manager = getattr(self, "_panel_manager", None)
        if not getattr(self, "_neighbourhood_window", False) or manager is None:
            return ""
        from ..panels.base import NEIGHBOURHOOD

        for uid in manager.uids(NEIGHBOURHOOD):
            if manager.is_open(uid):
                return manager.title(uid)
        return ""

    def _mark_neighbourhood_window(self) -> None:
        """This window was opened for a comparison (``open_neighbourhood_window``):
        its title follows the comparison, and it is not reopened with the run."""
        self._neighbourhood_window = True
        manager = self._panel_manager
        if manager is not None:
            manager.instancesChanged.connect(self._apply_window_title)
            manager.layoutApplied.connect(self._apply_window_title)
            manager.host.viewOpened.connect(self._apply_window_title)
            manager.host.viewClosed.connect(self._apply_window_title)
        self._apply_window_title()

    def _clone_rail_from(self, source: Any) -> None:
        """The rail's settings that are not filters, copied from ``source``
        (the neighbourhood window): the sample sort and the rail's width and
        collapsed state. The filters come with ``initial_filters``."""
        for name in ("sample_sort_column_box", "sample_sort_order_box"):
            mine = getattr(self, name, None)
            theirs = getattr(source, name, None)
            if not isinstance(mine, QtWidgets.QComboBox) or not isinstance(theirs, QtWidgets.QComboBox):
                continue
            index = mine.findData(theirs.currentData())
            if index >= 0 and index != mine.currentIndex():
                mine.setCurrentIndex(index)
        rail_state = getattr(source, "_rail_state", None)
        if callable(rail_state):
            width, collapsed = rail_state()
            self._apply_rail_state(width, collapsed)

    def _show_filters_in_rail(self, filters: FilterState) -> None:
        """Start the rail with ``filters`` (Duplicate Window): the rail widgets
        show them and the context's draft is them, so the first compute
        applies exactly these filters (with the window subset)."""
        widgets = [
            self.sample_query, self.gene_query, self.product_query, self.molecule_filter_box,
        ]
        models = [m for m in (self._predictor_model, self._tag_model) if m is not None]
        blockers = [QtCore.QSignalBlocker(w) for w in widgets + models]
        try:
            self.sample_query.setText(filters.sample_query)
            self.gene_query.setText(filters.gene)
            self.product_query.setText(filters.product)
            index = self.molecule_filter_box.findData(filters.molecule)
            self.molecule_filter_box.setCurrentIndex(max(0, index))
            for model, chosen in ((self._predictor_model, filters.predictors), (self._tag_model, filters.tags)):
                if model is None:
                    continue
                for row in range(model.rowCount()):
                    item = model.item(row)
                    token = item.data(QtCore.Qt.ItemDataRole.UserRole) if item is not None else None
                    if not token or token == "__clear__":
                        continue
                    checked = str(token) in chosen
                    item.setCheckState(QtCore.Qt.CheckState.Checked if checked else QtCore.Qt.CheckState.Unchecked)
        finally:
            for blocker in blockers:
                blocker.unblock()
        self._update_predictor_box_text()
        self._update_tag_box_text()
        self.metadata_filters = filters.metadata_filters()
        self.ctx.update_draft(filters, origin="rail")


    # The paired Annotations panel's view mode (D3: a display option of that
    # panel, not a filter); the genome track follows it.

    @property
    def view_mode(self) -> str:
        panel = getattr(self, "annotations_panel", None)
        return panel.view_mode() if panel is not None else "Contig"

    # Actions
    def refresh(self) -> None:
        """Apply what the rail shows now (computes even if nothing changed)."""
        self.ctx.update_draft(FilterState.from_widgets(self), force=True)


    def _on_apply_failed(self, error_text: str) -> None:
        print(error_text)
        self.status.setText("Filter refresh failed. See terminal for traceback.")


    def _on_result_applied(self, payload: Tuple[FilterState, AppliedResult]) -> None:
        """The window's chrome for a new result: the shared status bar (the
        rail and every panel have their own subscriptions)."""
        _filters, result = payload
        self.status.setText(result.stats_text)
        self._update_active_panel_status()

    def _show_hover_in_status(self, hover: Any) -> None:
        samples = tuple(getattr(hover, "samples", ()) or ())
        label = getattr(self, "hover_status", None)
        if label is None:
            return
        if len(samples) == 1:
            label.setText(f"Hovered: {samples[0]}")
        elif samples:
            label.setText(f"Hovered: {len(samples)} samples")
        else:
            label.setText("")

    def _show_busy(self, busy: bool) -> None:
        """A compute runs: show "Updating..." when it takes a while."""
        label = getattr(self, "busy_status", None)
        timer = getattr(self, "_busy_timer", None)
        if label is None or timer is None:
            return
        if busy:
            timer.start()
        else:
            timer.stop()
            label.setVisible(False)
