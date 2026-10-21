"""Panel protocol, registry and frames of an explore window (plan 4.8, Phase 4).

- ``ExplorePanel`` is the protocol every panel follows: it is bound to its
  window's ``ExploreContext`` once, at creation, and unbound once, on close; it
  draws from the applied result, shows the hover and the selection, and saves
  and restores its own display options.
- ``PanelSpec`` describes one panel kind (title, palette group, heavy or
  light, single or multi instance, the factory that builds its panel, and
  whether the run has the data it needs). ``PanelRegistry`` holds the specs of
  a window; it drives the palette, the menus and the "+" buttons. The filter
  rail is not in it (it is window chrome, D11).
- ``PanelFrame`` is the plain widget that holds one panel instance. ``DockHost``
  wraps it in a dock. The frame builds its panel (light panels at once, heavy panels
  lazily on first show), owns the panel's subscriptions on the context's bus
  and gates delivery by visibility: a heavy panel in a background tab, a
  closed panel or a group hidden by Maximize keeps the latest result, marks
  itself dirty and drops hover; when it becomes visible it replays the latest
  result and the current hover.

The first instance of a kind has the kind as uid (``geo``); further instances
of multi-instance kinds get a suffix (``geo-91c0``). Document kinds
(``PanelSpec.document``: the neighbourhood and alignment panels) have no
first instance: each is opened with its own data, always gets a suffixed uid
and is disposed on close.
"""
from __future__ import annotations

import logging
import secrets
from dataclasses import dataclass, field
from typing import Any, Callable, ClassVar, Dict, Iterable, List, Optional, Sequence, Tuple

from PyQt6 import QtCore, QtWidgets, sip

from ..session import bus as bus_events
from ..session.bus import Delivery

_LOG = logging.getLogger(__name__)

SAMPLE_TREE = "sample_tree"
ANNOTATIONS = "annotations"
TRACK = "track"
GEO = "geo"
CGMLST_FORCE = "cgmlst_force"
# Pan panels (plan 8.10, Phase 9): one pan view each. The Gene
# presence/absence rows of a window follow its Phylogeny's leaf order
# (panels/row_link.py); the heatmap always fits every row.
PHYLOGENY = "phylogeny"
GENE_PRESENCE_ABSENCE = "gene_presence_absence"
PAN_METADATA = "pan_metadata"
PAN_FORCE = "pan_force"
PLASMIDOME = "plasmidome"
PAN_VIEW_KINDS = (PHYLOGENY, GENE_PRESENCE_ABSENCE, PAN_METADATA, PAN_FORCE, PLASMIDOME)
# Pangenome clustering (panels/pan_clustering.py): MDS and k-means of the
# isolates over a prevalence range of the Panaroo genes; not a pan view tab.
PAN_MDS = "pan_mds"
# neoPan, the Pan tree and neoTree of older builds became Phylogeny, Gene
# presence became Gene presence/absence (saved layouts are migrated on load,
# shell/layout_store.py RETIRED_KINDS).
# Document panels (plan 8.7): opened with their data, many per window.
NEIGHBOURHOOD = "neighbourhood"
ALIGNMENT = "alignment"


def _always_available() -> Tuple[bool, str]:
    return True, ""


class ExplorePanel(QtWidgets.QWidget):
    """One view of an explore window (plan 4.8).

    A panel is bound to exactly one context for its whole life; there is no
    move or rebind call. Panels never call each other: they change shared
    state through the context's mutation API (``hover``, ``select``,
    ``set_subset``, ``request``) and learn about changes through the methods
    below, which their frame calls.
    """

    KIND: ClassVar[str] = ""
    TITLE: ClassVar[str] = ""
    SINGLE_INSTANCE: ClassVar[bool] = False
    HEAVY: ClassVar[bool] = False
    NEEDS_PARAMS: ClassVar[bool] = False

    # What ``title_detail`` returns changed (the manager retitles the tab).
    titleDetailChanged = QtCore.pyqtSignal()

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.ctx: Any = None
        self.uid = ""

    # Life cycle
    def bind(self, ctx: Any, uid: str) -> None:
        """Once, at creation."""
        self.ctx = ctx
        self.uid = str(uid)

    def unbind(self) -> None:
        """Once, on close."""
        self.ctx = None

    def dispose(self) -> None:
        """The panel is closed for good: drop timers, views and caches."""

    # Shared state in
    def apply_result(self, filters: Any, result: Any) -> None:
        """Draw from the applied filters and their result."""

    def show_hover(self, hover: Any) -> None:
        """Show the hovered samples (never called with this panel's own hover)."""

    def show_selection(self, selection: Any) -> None:
        """The shared selection changed."""

    def on_shown(self) -> None:
        """The panel became visible (after the replay of result and hover)."""

    def on_hidden(self) -> None:
        """The panel was hidden (a background tab, closed, or Maximize)."""

    # Display options (panel params, plan 3.11)
    def save_state(self) -> Dict[str, Any]:
        return {}

    def restore_state(self, state: Dict[str, Any]) -> None:
        del state

    def title(self) -> str:
        return self.TITLE

    def title_detail(self) -> str:
        """A few words on what this instance shows ("Core genome", "All
        matches"), or "". With several instances of a kind open in a window,
        each default title is "<Title> · <detail>" (``PanelManager``); emit
        ``titleDetailChanged`` when it changes."""
        return ""


# Builds the panel of one instance (its uid) as a child of its frame, so a
# web or GL view inside it is created in the window (never as a top-level
# widget that is re-parented later); None when the run lacks what it needs.
PanelFactory = Callable[[str, QtWidgets.QWidget], Optional[ExplorePanel]]


@dataclass(frozen=True)
class PanelSpec:
    """One panel kind: drives the palette, the menus and the frames."""

    kind: str
    title: str
    group: str  # palette group (plan 3.6)
    heavy: bool = False  # built lazily on first show; results only while visible
    placeholder: str = ""
    available: Callable[[], Tuple[bool, str]] = field(default=_always_available)
    single_instance: bool = True
    # Builds the panel of one instance; None leaves the frame a placeholder.
    factory: Optional[PanelFactory] = None
    # A document kind (plan 3.11, 8.7): instances are opened with their data
    # (``PanelSet.open_document``), never built empty, never part of the
    # palette's kinds or of a saved layout.
    document: bool = False
    # The uid of a panel this one cannot exist without (the Genome track
    # needs the first Annotations panel, which feeds it). PanelManager opens
    # that panel with this one, closes this one when that one closes, and
    # closes this one after a restore or preset that leaves it alone.
    requires: str = ""
    # The frame's minimum height. The dock area does not respect it (it can
    # be dragged down to about 70 px with its title bar) and clips a taller
    # frame from the bottom; a panel whose bottom matters (the Genome track's
    # overview) asks for less and handles short heights itself.
    min_height: int = 90
    # The ``title_detail`` of an instance not built yet (a background tab),
    # from the state it will be built with ({} for none); None: "".
    state_detail: Optional[Callable[[Dict[str, Any]], str]] = None

    @property
    def uid(self) -> str:
        """The uid of the first instance (the kind itself, stable across runs)."""
        return self.kind


class PanelRegistry:
    """The panel kinds of one explore window, in palette order."""

    def __init__(self, specs: Iterable[PanelSpec]) -> None:
        self._specs: Dict[str, PanelSpec] = {}
        for spec in specs:
            if spec.kind in self._specs:
                raise ValueError(f"duplicate panel kind {spec.kind!r}")
            self._specs[spec.kind] = spec

    def specs(self) -> List[PanelSpec]:
        return list(self._specs.values())

    def kinds(self) -> List[str]:
        return list(self._specs)

    def spec(self, kind: str) -> Optional[PanelSpec]:
        return self._specs.get(kind)

    def __contains__(self, kind: object) -> bool:
        return kind in self._specs

    def is_multi_instance(self, kind: str) -> bool:
        spec = self._specs.get(kind)
        # Document panels are many per window too, but each holds its own
        # data: they are never duplicated or opened empty.
        return spec is not None and not spec.single_instance and not spec.document

    def available(self, kind: str) -> Tuple[bool, str]:
        spec = self._specs.get(kind)
        if spec is None:
            return False, "Unknown panel"
        available = spec.available()
        required = self._specs.get(PanelRegistry.kind_of(spec.requires)) if spec.requires else None
        if available[0] and required is not None and not required.available()[0]:
            return False, f"Needs {required.title}, which is not available for this run"
        return available

    def dependents(self, uid: str) -> List[str]:
        """The kinds that require the panel ``uid`` (``PanelSpec.requires``)."""
        return [spec.kind for spec in self._specs.values() if spec.requires and spec.requires == uid]

    @staticmethod
    def new_uid(kind: str, taken: Iterable[str]) -> str:
        """The kind for the first instance, a short random suffix after that."""
        used = set(taken)
        if kind not in used:
            return kind
        while True:
            uid = f"{kind}-{secrets.token_hex(2)}"
            if uid not in used:
                return uid

    @staticmethod
    def kind_of(uid: str) -> str:
        return str(uid).split("-", 1)[0]


def _alive(obj: Optional[QtCore.QObject]) -> bool:
    return obj is not None and not sip.isdeleted(obj)


class PanelFrame(QtWidgets.QWidget):
    """The content widget of one dock: a placeholder until its panel is built.

    ``shownChanged`` reports real visibility (a background tab, a closed panel
    or a group hidden by Maximize counts as hidden).
    """

    shownChanged = QtCore.pyqtSignal(str, bool)
    panelBuilt = QtCore.pyqtSignal(str)
    # The panel's ``title_detail`` changed (uid).
    titleDetailChanged = QtCore.pyqtSignal(str)

    def __init__(
        self,
        spec: PanelSpec,
        parent: Optional[QtWidgets.QWidget] = None,
        *,
        uid: Optional[str] = None,
    ) -> None:
        super().__init__(parent)
        self.spec = spec
        self.uid = str(uid or spec.uid)
        self.kind = spec.kind
        self.setObjectName("PanelFrame")
        self.setProperty("panelKind", spec.kind)
        # Groups have minimum sizes, so a split never squashes a panel to zero.
        self.setMinimumSize(160, spec.min_height)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self._content: Optional[QtWidgets.QWidget] = None
        self._panel: Optional[ExplorePanel] = None
        self._ctx: Any = None
        self._ready: Callable[[], bool] = lambda: True
        self._subs: List[Any] = []
        self._result_sub: Any = None
        self._pending_state: Optional[Dict[str, Any]] = None
        self._building = False
        self._build_failed = False
        self._disposed = False
        self._settle_pending = False
        if spec.placeholder:
            self.set_placeholder(spec.placeholder)

    # ------------------------------------------------------------ content
    def content(self) -> Optional[QtWidgets.QWidget]:
        return self._content

    def panel(self) -> Optional[ExplorePanel]:
        return self._panel

    def title_detail(self) -> str:
        """The panel's ``title_detail``; before the build, the spec's
        ``state_detail`` of the state the panel will be built with."""
        if self._disposed:
            return ""
        panel = self._panel
        try:
            if panel is not None:
                return str(panel.title_detail() or "").strip()
            if self.spec.state_detail is not None:
                return str(self.spec.state_detail(dict(self._pending_state or {})) or "").strip()
            return ""
        except Exception:
            _LOG.exception("panel %s: title_detail failed", self.uid)
            return ""

    def _on_panel_title_detail_changed(self) -> None:
        if not self._disposed:
            self.titleDetailChanged.emit(self.uid)

    def is_built(self) -> bool:
        return self._content is not None

    def is_bound(self) -> bool:
        return self._ctx is not None

    def _clear(self) -> None:
        layout = self.layout()
        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget is not None and widget is not self._content:
                widget.setParent(None)
                widget.deleteLater()

    def set_placeholder(self, text: str) -> None:
        self._clear()
        label = QtWidgets.QLabel(text, self)
        label.setObjectName("PanelPlaceholder")
        label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        label.setWordWrap(True)
        label.setStyleSheet("color: #64748b; font-size: 13px; padding: 18px;")
        self.layout().addWidget(label)

    def set_content(self, widget: QtWidgets.QWidget) -> None:
        self._clear()
        self._content = widget
        self.layout().addWidget(widget)

    # ------------------------------------------------------------ binding
    def bind(self, ctx: Any, *, ready: Optional[Callable[[], bool]] = None) -> None:
        """Subscribe this frame's panel to ``ctx``'s bus (once, at creation)."""
        if self._ctx is not None or self._disposed:
            return
        self._ctx = ctx
        if ready is not None:
            self._ready = ready
        bus = ctx.bus
        heavy = self.spec.heavy
        self._result_sub = bus.subscribe(
            self.uid,
            bus_events.RESULT_APPLIED,
            self._deliver_result,
            delivery=Delivery.WHEN_VISIBLE if heavy else Delivery.IMMEDIATE,
            delay_ms=bus_events.WHEN_VISIBLE_DELAY_MS,
            visible=self._can_take_result,
        )
        hover_sub = bus.subscribe(
            self.uid,
            bus_events.HOVER_CHANGED,
            self._deliver_hover,
            delivery=Delivery.COALESCED if heavy else Delivery.IMMEDIATE,
            coalesce_ms=bus_events.HEAVY_HOVER_COALESCE_MS,
            visible=self._is_shown,
            drop_when_hidden=True,
        )
        selection_sub = bus.subscribe(self.uid, bus_events.SELECTION_CHANGED, self._deliver_selection)
        self._subs = [self._result_sub, hover_sub, selection_sub]
        if not heavy:
            self.ensure_built()
        result = getattr(ctx, "result", None)
        if result is not None:
            # A panel opened after the first result catches up at once (light)
            # or when it is shown (heavy).
            self._result_sub.offer((ctx.applied, result))

    def unbind(self) -> None:
        """Stop every delivery to this frame and unbind its panel."""
        ctx, self._ctx = self._ctx, None
        if ctx is not None:
            ctx.bus.unsubscribe_owner(self.uid)
        self._subs = []
        self._result_sub = None
        if self._panel is not None:
            try:
                self._panel.unbind()
            except Exception:
                _LOG.exception("panel %s: unbind failed", self.uid)

    def dispose(self) -> None:
        """Closed for good (a duplicate's close): unbind, then dispose the panel."""
        if self._disposed:
            return
        self._disposed = True
        panel = self._panel
        self.unbind()
        if panel is not None:
            try:
                panel.dispose()
            except Exception:
                _LOG.exception("panel %s: dispose failed", self.uid)

    def is_disposed(self) -> bool:
        return self._disposed

    def is_dirty(self) -> bool:
        """A result is waiting for this (hidden or not yet built) panel."""
        return bool(self._result_sub is not None and self._result_sub.dirty)

    # ------------------------------------------------------------ building
    def ensure_built(self) -> bool:
        """Build the panel if it can be built now (heavy panels wait until the
        window says the first heavy builds may start)."""
        if self._panel is not None:
            return True
        if self._disposed or self._building or self._build_failed or self.spec.factory is None:
            return False
        if self.spec.heavy and not self.spec.document and not self._ready():
            # (A document panel is opened on demand and shown at once.)
            return False
        self._building = True
        try:
            panel = self.spec.factory(self.uid, self)
        except Exception as exc:
            _LOG.exception("panel %s: build failed", self.uid)
            self._build_failed = True
            self.set_placeholder(f"{self.spec.title} could not be loaded: {exc}")
            return False
        finally:
            self._building = False
        if panel is None:
            # Not buildable for this run (for example no output folder); do
            # not try again.
            self._build_failed = True
            available, reason = self.spec.available()
            self.set_placeholder(
                reason if (not available and reason) else f"{self.spec.title} is not available for this run."
            )
            return False
        self._panel = panel
        if self._content is not panel:
            self.set_content(panel)
        detail_changed = getattr(panel, "titleDetailChanged", None)
        if detail_changed is not None:
            detail_changed.connect(self._on_panel_title_detail_changed)
        if self._ctx is not None:
            panel.bind(self._ctx, self.uid)
        if self._pending_state is not None:
            state, self._pending_state = self._pending_state, None
            self._restore_panel_state(panel, state)
        self.panelBuilt.emit(self.uid)
        return True

    # ------------------------------------------------------------ state
    def save_state(self) -> Dict[str, Any]:
        if self._panel is not None:
            try:
                return dict(self._panel.save_state() or {})
            except Exception:
                _LOG.exception("panel %s: save_state failed", self.uid)
                return {}
        return dict(self._pending_state or {})

    def restore_state(self, state: Dict[str, Any]) -> None:
        state = dict(state or {})
        if self._panel is None:
            self._pending_state = state
            self._on_panel_title_detail_changed()
            return
        self._restore_panel_state(self._panel, state)

    def _restore_panel_state(self, panel: ExplorePanel, state: Dict[str, Any]) -> None:
        try:
            panel.restore_state(state)
        except Exception:
            _LOG.exception("panel %s: restore_state failed", self.uid)

    # ------------------------------------------------------------ delivery
    def _is_shown(self) -> bool:
        return _alive(self) and self.isVisible()

    def _can_take_result(self) -> bool:
        if not self._is_shown():
            return False
        return self._panel is not None or (self.spec.factory is not None and not self._build_failed and self._ready())

    def _deliver_result(self, payload: Any) -> None:
        if not self.ensure_built():
            return
        filters, result = payload
        self._panel.apply_result(filters, result)
        # A redraw can drop highlights: show the current hover again.
        self._replay_hover()

    def _deliver_hover(self, hover: Any) -> None:
        if self._panel is not None:
            self._panel.show_hover(hover)

    def _deliver_selection(self, selection: Any) -> None:
        if self._panel is not None:
            self._panel.show_selection(selection)

    def _replay_hover(self) -> None:
        ctx = self._ctx
        if ctx is None or self._panel is None or not self._is_shown():
            return
        hover = getattr(ctx, "hover_state", None)
        if hover is None or (hover.origin is not None and hover.origin == self.uid):
            return
        self._panel.show_hover(hover)

    def settle_shown(self) -> None:
        """The frame is visible and settled: build (lazy), replay the stored
        result, show the current hover. Also called when the window allows the
        first heavy builds."""
        self._settle_pending = False
        if not self._is_shown() or self._ctx is None:
            return
        if not self.ensure_built():
            return
        self._ctx.bus.notify_shown(self.uid)
        self._replay_hover()
        try:
            self._panel.on_shown()
        except Exception:
            _LOG.exception("panel %s: on_shown failed", self.uid)

    def showEvent(self, event) -> None:  # type: ignore[override]
        super().showEvent(event)
        self.shownChanged.emit(self.uid, True)
        # Deliver after the show event (and any re-parenting on a drop) settles.
        if not self._settle_pending:
            self._settle_pending = True
            QtCore.QTimer.singleShot(0, self._settle_if_alive)

    def _settle_if_alive(self) -> None:
        if _alive(self):
            self.settle_shown()

    def hideEvent(self, event) -> None:  # type: ignore[override]
        super().hideEvent(event)
        self.shownChanged.emit(self.uid, False)
        if self._panel is not None:
            try:
                self._panel.on_hidden()
            except Exception:
                _LOG.exception("panel %s: on_hidden failed", self.uid)


class PanelSet(QtCore.QObject):
    """The panel instances of one explore window and their binding.

    The base of ``PanelManager``, which adds the docking area, the menus and
    the palette; frames are created when panels open.
    """

    panelBuilt = QtCore.pyqtSignal(str)

    def __init__(
        self,
        parent: Optional[QtCore.QObject],
        specs: "Sequence[PanelSpec] | PanelRegistry",
        *,
        ctx: Any = None,
    ) -> None:
        super().__init__(parent)
        self.registry = specs if isinstance(specs, PanelRegistry) else PanelRegistry(specs)
        self.ctx = ctx
        self._frames: Dict[str, PanelFrame] = {}
        # The data of open document panels, by uid (never saved).
        self._documents: Dict[str, Any] = {}
        self._ready = False

    # ------------------------------------------------------------- frames
    def _create_frame(self, spec: PanelSpec, uid: str, state: Optional[Dict[str, Any]] = None) -> PanelFrame:
        frame = PanelFrame(spec, uid=uid)
        if state:
            # Restored right after the build, before the first result is drawn.
            frame.restore_state(state)
        self._frames[uid] = frame
        frame.panelBuilt.connect(self.panelBuilt)
        self._on_frame_created(frame)
        if self.ctx is not None:
            frame.bind(self.ctx, ready=self.is_ready)
        return frame

    def _on_frame_created(self, frame: PanelFrame) -> None:
        """Hook for the docking manager (registers the frame's dock)."""

    def _drop_frame(self, uid: str) -> Optional[PanelFrame]:
        frame = self._frames.pop(uid, None)
        self._documents.pop(uid, None)
        if frame is not None:
            frame.dispose()
        return frame

    # ------------------------------------------------------------ documents
    def is_document(self, uid: str) -> bool:
        frame = self._frames.get(uid)
        return frame is not None and frame.spec.document

    def document(self, uid: str) -> Any:
        """The data a document panel was opened with (its factory reads it)."""
        return self._documents.get(uid)

    def find_document(self, kind: str, key: Any) -> Optional[str]:
        """The uid of the open document panel of ``kind`` with ``key``."""
        if key is None:
            return None
        for uid, document in self._documents.items():
            frame = self._frames.get(uid)
            if frame is not None and frame.kind == kind and getattr(document, "key", None) == key:
                return uid
        return None

    def _create_document_frame(
        self,
        kind: str,
        document: Any,
        *,
        uid: Optional[str] = None,
        state: Optional[Dict[str, Any]] = None,
    ) -> Optional[PanelFrame]:
        spec = self.registry.spec(kind)
        if spec is None or not spec.document:
            return None
        taken = set(self._frames) | {kind}
        if uid is None or uid in taken or PanelRegistry.kind_of(uid) != kind:
            uid = PanelRegistry.new_uid(kind, taken)
        self._documents[uid] = document
        frame = self._create_frame(spec, uid, state)
        # Built at once (it is shown right away), so its title is known.
        frame.ensure_built()
        return frame

    def document_title(self, uid: str) -> str:
        frame = self._frames.get(uid)
        panel = frame.panel() if frame is not None else None
        if panel is not None:
            return panel.title()
        document = self._documents.get(uid)
        tab_title = getattr(document, "tab_title", None)
        if callable(tab_title) and tab_title():
            return str(tab_title())
        return str(getattr(document, "title", "") or (frame.spec.title if frame is not None else uid))

    def frame(self, uid: str) -> Optional[PanelFrame]:
        return self._frames.get(uid)

    def frames(self, kind: Optional[str] = None) -> List[PanelFrame]:
        return [f for f in self._frames.values() if kind is None or f.kind == kind]

    def uids(self, kind: Optional[str] = None) -> List[str]:
        return [f.uid for f in self.frames(kind)]

    def panel(self, uid: str) -> Optional[ExplorePanel]:
        frame = self._frames.get(uid)
        return frame.panel() if frame is not None else None

    def panels(self, kind: Optional[str] = None) -> List[ExplorePanel]:
        return [f.panel() for f in self.frames(kind) if f.panel() is not None]

    def is_visible(self, uid: str) -> bool:
        """Really on screen: open, the current tab of its group and not hidden
        by Maximize (this is what panel delivery gates on)."""
        frame = self._frames.get(uid)
        return frame is not None and _alive(frame) and frame.isVisible()

    # -------------------------------------------------------- lazy builds
    def is_ready(self) -> bool:
        return self._ready

    def set_ready(self) -> None:
        """The main view is ready: heavy panels may build when shown; the
        visible ones build now."""
        self._ready = True
        for frame in list(self._frames.values()):
            if frame.isVisible():
                frame.settle_shown()

    def has_pending_delivery(self) -> bool:
        return self.ctx is not None and self.ctx.bus.has_pending()
