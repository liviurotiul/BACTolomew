"""The neighbourhood panel (plan 4.6, 8.7): a gene neighborhood or plasmid
context comparison docked in its explore window.

A document panel: every comparison is its own instance, opened with its rows
(or a loader that computes them, shown as an inline busy state meanwhile)
from a gene card ("Compare neighborhoods") or a pan view (plasmid context).
Opening the same comparison again focuses the open panel. Its rows are not
saved with the layout (3.11); its display options are its panel params.

Linking (4.6): the rows follow the window's one sample set (D4: the viewer's
own metadata filters are gone; the rail's active filters are shown read only
in the Show card); a hovered row hovers its samples in the window, and
samples hovered elsewhere tint their rows; "Open sample(s) in table" narrows
the window through ``set_subset``. The panel has no tree of its own: a
Phylogeny panel of the window highlights the hovered rows' genomes through
the window's hover. Ctrl+W and Ctrl+F act only inside the panel (3.7), and a
second Esc closes nothing.
"""
from __future__ import annotations

import traceback
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Hashable, List, Optional, Sequence, Tuple

from PyQt6 import QtCore, QtWidgets

from ..constants import MISSING_FILTER_TOKEN
from ..model.filter_state import SampleSubset
from ..session.requests import ClosePanel, OpenGeneCard, ShowAnnotationsTable
from .base import NEIGHBOURHOOD, ExplorePanel
from .pan import hover_kind

# The display options a neighbourhood panel saves (plan 3.11 panel params).
_STATE_KEYS = (
    "flank",
    "color_by",
    "labels",
    "mirror",
    "molecule",
    "sort_by",
    "sort_order",
    "split_molecules",
    "merge",
    "settings_open",
)


@dataclass
class NeighbourhoodDocument:
    """What one neighbourhood panel compares (kept in memory, never saved).

    ``rows`` are the comparison rows; until they exist ``loader`` computes
    them on a worker thread. ``key`` identifies the comparison: opening the
    same key again focuses the open panel.
    """

    title: str
    rows: Optional[List[Dict[str, Any]]] = None
    loader: Optional[Callable[[], List[Dict[str, Any]]]] = None
    source: str = "plasmid"  # "plasmid" (a pan view) or "gene_card"
    key: Optional[Hashable] = None
    flank_sizes: Optional[Tuple[int, ...]] = None
    initial_flank: Optional[int] = None
    # Shown when the loader finds no rows.
    empty_text: str = "No neighborhoods were found for this gene."
    # Added to the tab title when another open comparison has the same one
    # (two gene cards of one gene in different samples).
    tab_suffix: str = ""
    extra: Dict[str, Any] = field(default_factory=dict)
    # The collection running for these rows, shared by every panel showing
    # this document (a panel closed while loading and reopened joins it).
    job: Optional["_LoadJob"] = field(default=None, repr=False, compare=False)

    def tab_title(self) -> str:
        """Short, subject first, so narrow tabs keep what tells comparisons
        apart: 'blaOXA-23 neighborhoods', 'pX1 plasmid context'."""
        heading, _sep, subject = (self.title or "").partition(":")
        subject = subject.strip()
        if not subject:
            text = (self.title or "").strip()
        elif heading.strip().lower() == "gene neighborhoods":
            text = f"{subject} neighborhoods"
        elif heading.strip().lower() == "plasmid context":
            text = f"{subject} plasmid context"
        else:
            text = (self.title or "").strip()
        return f"{text} ({self.tab_suffix})" if text and self.tab_suffix else text


def filter_notes(filters: Any) -> List[str]:
    """Short read-only texts for the filters that narrow a window's sample
    set (the rail's sample search, metadata filters, subsets and annotation
    filters)."""
    if filters is None:
        return []
    notes: List[str] = []

    def values_text(values: Sequence[str]) -> str:
        shown = ["(missing)" if v == MISSING_FILTER_TOKEN else str(v) for v in values]
        if len(shown) > 3:
            return ", ".join(shown[:3]) + f" (+{len(shown) - 3})"
        return ", ".join(shown)

    if getattr(filters, "sample_query", ""):
        notes.append(f'Search: "{filters.sample_query}"')
    for column, values in getattr(filters, "metadata", ()) or ():
        notes.append(f"{column}: {values_text(values)}")
    for subset in (getattr(filters, "window_subset", None), getattr(filters, "subset", None)):
        if subset is not None:
            notes.append(subset.label or f"{len(subset.samples)} samples")
    if getattr(filters, "gene", ""):
        notes.append(f'Gene: "{filters.gene}"')
    if getattr(filters, "product", ""):
        notes.append(f'Product: "{filters.product}"')
    if getattr(filters, "molecule", ""):
        notes.append(f"Molecule: {filters.molecule}")
    if getattr(filters, "predictors", None):
        notes.append(f"Predictors: {values_text(sorted(filters.predictors))}")
    if getattr(filters, "tags", None):
        notes.append(f"Tags: {values_text(sorted(filters.tags))}")
    if getattr(filters, "hide_hypothetical", False):
        notes.append("No hypothetical proteins")
    return notes


class _LoaderSignals(QtCore.QObject):
    finished = QtCore.pyqtSignal(int, object)
    failed = QtCore.pyqtSignal(int, str)


class _LoaderWorker(QtCore.QRunnable):
    def __init__(self, loader: Callable[[], List[Dict[str, Any]]]) -> None:
        super().__init__()
        self._loader = loader
        self.signals = _LoaderSignals()
        self.setAutoDelete(True)

    def run(self) -> None:
        try:
            rows = self._loader()
        except Exception:
            self.signals.failed.emit(0, traceback.format_exc())
            return
        self.signals.finished.emit(0, rows)


class _LoadJob(QtCore.QObject):
    """One collection of a document's rows on the thread pool. It belongs to
    the document, not to a panel: closing the panel while it runs and
    reopening it joins the running collection instead of starting a second
    one over the same gene card's caches. When it ends, the rows (or the
    error) are stored on the document."""

    finished = QtCore.pyqtSignal(object)
    failed = QtCore.pyqtSignal(str)

    def __init__(self, document: "NeighbourhoodDocument") -> None:
        super().__init__()
        self.document = document
        worker = _LoaderWorker(document.loader)
        worker.signals.finished.connect(self._on_finished)
        worker.signals.failed.connect(self._on_failed)
        self._worker: Optional[_LoaderWorker] = worker
        QtCore.QThreadPool.globalInstance().start(worker)

    def is_running(self) -> bool:
        return self._worker is not None

    def _on_finished(self, _token: int, rows_obj: object) -> None:
        self._worker = None
        rows = [dict(row) for row in rows_obj if isinstance(row, dict)] if isinstance(rows_obj, list) else []
        document = self.document
        document.rows = rows
        document.loader = None
        if document.job is self:
            document.job = None
        self.finished.emit(rows)

    def _on_failed(self, _token: int, error_text: str) -> None:
        self._worker = None
        if self.document.job is self:
            self.document.job = None
        self.failed.emit(error_text)


class NeighbourhoodPanel(ExplorePanel):
    KIND = NEIGHBOURHOOD
    TITLE = "Gene neighborhoods"
    HEAVY = True
    NEEDS_PARAMS = True
    # The Settings column choice of this session for neighbourhood panels
    # (see NeighbourhoodView); each panel also saves its own as a param.
    settings_panel_open: Optional[bool] = None

    def __init__(
        self,
        document: NeighbourhoodDocument,
        *,
        view_options: Optional[Dict[str, Any]] = None,
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self.setObjectName("NeighbourhoodPanel")
        self.document = document
        self._view_options = dict(view_options or {})
        self.view: Any = None  # the NeighbourhoodView, once the rows exist
        self._pending_state: Dict[str, Any] = {}
        self._scope: Optional[Tuple[Tuple[str, ...], List[str]]] = None
        self._hover_samples: Tuple[str, ...] = ()
        self._hovering = False
        self._job: Optional[_LoadJob] = None
        self._disposed = False
        self._stack = QtWidgets.QStackedLayout(self)
        self._stack.setContentsMargins(0, 0, 0, 0)
        self.message_page = QtWidgets.QWidget(self)
        message_layout = QtWidgets.QVBoxLayout(self.message_page)
        message_layout.addStretch(1)
        self.message = QtWidgets.QLabel("", self.message_page)
        self.message.setObjectName("NeighbourhoodPanelMessage")
        self.message.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.message.setWordWrap(True)
        self.message.setProperty("role", "muted")
        message_layout.addWidget(self.message)
        self.busy_bar = QtWidgets.QProgressBar(self.message_page)
        self.busy_bar.setObjectName("NeighbourhoodPanelBusy")
        self.busy_bar.setRange(0, 0)  # indeterminate
        self.busy_bar.setTextVisible(False)
        self.busy_bar.setFixedWidth(220)
        self.busy_bar.setFixedHeight(6)
        message_layout.addWidget(self.busy_bar, alignment=QtCore.Qt.AlignmentFlag.AlignHCenter)
        # After a failed load: collect the rows again.
        self.retry_button = QtWidgets.QPushButton("Try again", self.message_page)
        self.retry_button.setObjectName("NeighbourhoodPanelRetry")
        self.retry_button.clicked.connect(self._load)
        message_layout.addWidget(self.retry_button, alignment=QtCore.Qt.AlignmentFlag.AlignHCenter)
        message_layout.addStretch(2)
        self._stack.addWidget(self.message_page)
        self._load()

    # ------------------------------------------------------------ content
    def is_busy(self) -> bool:
        return self._job is not None

    def title(self) -> str:
        return self.document.tab_title() or self.TITLE

    def _show_message(self, text: str, *, busy: bool = False, retry: bool = False) -> None:
        self.message.setText(text)
        self.busy_bar.setVisible(busy)
        self.retry_button.setVisible(retry)
        self._stack.setCurrentWidget(self.message_page)

    def _load(self) -> None:
        """Show the rows, or collect them first (the inline busy state that
        replaces the gene card's progress dialog)."""
        if self._disposed or self.view is not None or self._job is not None:
            return
        document = self.document
        if document.rows:
            self._build_view(document.rows)
            return
        if document.job is None and document.loader is None:
            self._show_message(document.empty_text)
            return
        self._show_message("Loading gene neighborhoods...", busy=True)
        job = document.job
        if job is None:
            job = document.job = _LoadJob(document)
        self._job = job
        job.finished.connect(self._on_rows_loaded)
        job.failed.connect(self._on_load_failed)

    def _detach_job(self) -> None:
        job, self._job = self._job, None
        if job is None:
            return
        for signal, slot in ((job.finished, self._on_rows_loaded), (job.failed, self._on_load_failed)):
            try:
                signal.disconnect(slot)
            except (TypeError, RuntimeError):
                pass

    def retry_with(self, document: NeighbourhoodDocument) -> None:
        """The same comparison asked for again (a second Compare
        neighborhoods click): a panel whose load failed or found nothing
        loads again, with the newer loader or rows. A panel showing rows or
        still loading is only focused."""
        if self._disposed or self.view is not None or self._job is not None:
            return
        if document is not self.document:
            if document.rows:
                self.document.rows = [dict(row) for row in document.rows]
                self.document.loader = None
            elif document.loader is not None:
                self.document.rows = None
                self.document.loader = document.loader
        self._load()

    def _on_rows_loaded(self, rows: object) -> None:
        self._detach_job()
        if self._disposed:
            return
        if not rows:
            self._show_message(self.document.empty_text)
            return
        self._build_view(list(rows))

    def _on_load_failed(self, error_text: str) -> None:
        self._detach_job()
        if self._disposed:
            return
        lines = error_text.strip().splitlines()
        self._show_message(
            "Could not build the gene neighborhoods.\n\n" + (lines[-1] if lines else "Unknown error."),
            retry=self.document.loader is not None,
        )

    def _build_view(self, rows: Sequence[Dict[str, Any]]) -> None:
        from ..igv_collection_view import NeighbourhoodView

        options = self._view_options
        document = self.document
        view = NeighbourhoodView(
            self,
            rows=rows,
            title=document.title,
            output_dir=options.get("output_dir"),
            flank_sizes=document.flank_sizes,
            initial_flank=document.initial_flank,
            open_samples_callback=self._open_samples_in_window,
            panel_mode=True,
            settings_memory=type(self),
        )
        self.view = view
        self._stack.addWidget(view)
        self._stack.setCurrentWidget(view)
        view.geneDoubleClicked.connect(self._on_gene_double_clicked)
        view.closeRequested.connect(self._request_close)
        view.view.contextSamplesHovered.connect(self._on_rows_hovered)
        if self._pending_state:
            state, self._pending_state = self._pending_state, {}
            self._apply_state(state)
        if self._scope is not None:
            view.set_sample_scope(*self._scope)
        if self._hover_samples:
            view.set_linked_samples(self._hover_samples)

    # ------------------------------------------------------------ shared state in
    def apply_result(self, filters: Any, result: Any) -> None:
        self._scope = (tuple(result.samples), filter_notes(filters))
        if self.view is not None:
            self.view.set_sample_scope(*self._scope)

    def show_hover(self, hover: Any) -> None:
        self._hover_samples = tuple(getattr(hover, "samples", ()) or ())
        if self.view is not None:
            self.view.set_linked_samples(self._hover_samples)

    # ------------------------------------------------------------ out
    def _on_rows_hovered(self, samples: object) -> None:
        if self.ctx is None:
            return
        values = [str(s).strip() for s in samples if str(s).strip()] if isinstance(samples, (list, tuple)) else []
        self._hovering = bool(values)
        self.ctx.hover(values, origin=self.uid or None, kind=hover_kind(values))

    def _on_gene_double_clicked(self, payload: Dict[str, Any]) -> None:
        from ..annotation_detail_dialog import normalize_detail_row

        if self.ctx is not None and payload:
            self.ctx.request(OpenGeneCard.from_row(normalize_detail_row(payload)), origin=self.uid)

    def _open_samples_in_window(self, samples: Sequence[str]) -> None:
        """Row menu "Open sample(s) in table": narrow this window's sample set
        to them (a chip in the rail), as every opener does (plan 8.7). A
        window without an Annotations panel (a neighbourhood window) also
        opens the Sample tree and Annotations, so the table shows."""
        values = [str(sample).strip() for sample in samples if str(sample).strip()]
        if not values or self.ctx is None:
            return
        label = f"Sample {values[0]}" if len(values) == 1 else f"{len(values)} samples"
        self.ctx.set_subset(SampleSubset(tuple(values), label=label, kind="sample"), origin=self.uid)
        self.ctx.request(ShowAnnotationsTable(self.uid or None), origin=self.uid)

    def _request_close(self) -> None:
        if self.ctx is not None:
            self.ctx.request(ClosePanel(self.uid), origin=self.uid)

    # ------------------------------------------------------------ life cycle
    def unbind(self) -> None:
        ctx = self.ctx
        if ctx is not None and self._hovering:
            ctx.clear_hover(origin=self.uid)
        super().unbind()

    def dispose(self) -> None:
        self._disposed = True
        # A collection still running goes on for the document (Reopen
        # Closed Panel joins it); this panel no longer listens.
        self._detach_job()

    # ------------------------------------------------------------ display options
    def save_state(self) -> Dict[str, Any]:
        view = self.view
        if view is None:
            return dict(self._pending_state)
        return {
            "flank": view._context_flank,
            "color_by": view._context_color_mode,
            "labels": view._context_label_mode,
            "mirror": bool(view.normalize_orientation_check.isChecked()),
            "molecule": view._molecule_filter,
            "sort_by": view.sort_column_box.currentData(),
            "sort_order": view.sort_order_box.currentData(),
            "split_molecules": bool(view.group_molecule_check.isChecked()),
            "merge": bool(view.collapse_btn.isChecked()),
            "settings_open": bool(view.settings_wanted()),
        }

    def restore_state(self, state: Dict[str, Any]) -> None:
        state = {k: v for k, v in dict(state or {}).items() if k in _STATE_KEYS}
        if self.view is None:
            self._pending_state.update(state)
            return
        self._apply_state(state)

    def _apply_state(self, state: Dict[str, Any]) -> None:
        view = self.view

        def choose(box: QtWidgets.QComboBox, value: Any) -> None:
            index = box.findData(value)
            if index >= 0 and index != box.currentIndex():
                box.setCurrentIndex(index)

        if "flank" in state and state["flank"] is not None:
            choose(view.flank_box, state["flank"])
        if "color_by" in state:
            choose(view.color_by_box, state["color_by"])
        if "labels" in state:
            choose(view.labels_box, state["labels"])
        if "mirror" in state:
            view.normalize_orientation_check.setChecked(bool(state["mirror"]))
        if "molecule" in state:
            choose(view.molecule_box, state["molecule"] or "")
        if "sort_by" in state:
            choose(view.sort_column_box, state["sort_by"])
        if "sort_order" in state:
            choose(view.sort_order_box, state["sort_order"])
        if "split_molecules" in state:
            view.group_molecule_check.setChecked(bool(state["split_molecules"]))
        if "merge" in state:
            view.collapse_btn.setChecked(bool(state["merge"]))
        if "settings_open" in state:
            view.set_settings_open(bool(state["settings_open"]))
