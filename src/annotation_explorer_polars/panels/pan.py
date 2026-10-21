"""Pan panels (plan 2.2, 3.6, 4.6 and 8.10, Phase 9).

Phylogeny, Gene presence/absence, Metadata summary, cgMLST MST (compact) and
Plasmidome:
each shows one pan view over the window's one sample set, as a panel that can
be docked anywhere in its explore window. They are heavy and multi-instance
(plan 3.6). Each wraps a ``PanCoreAccessoryTab`` built in single view mode by
the window's factory (``_create_pan_view_panel``) over the run's shared
``PanData`` (``run.pan``), so all pan panels of all windows share the run's
tree, fastGEAR, plasmid and hover data.

Phylogeny draws one of the run's trees (``PanData.tree_index``: the
recombination-filtered core tree by default, another from its source
dropdown) with the fastGEAR clades rotated together (clade backgrounds,
cluster brackets and cluster colored branches) and a clade and cluster
legend, rooted as its file says (an old unrooted file at its midpoint; the
Root dropdown picks another root); it greys the genomes outside the sample
set, as neoPan did. Gene presence/absence draws the Panaroo matrix, always with
every row fitted top to bottom. In one window the heatmap's rows are the
tree's leaves in the tree's order (``panels/row_link.py``, owned by the
window), but the tree's scroll and zoom never move them, and the heatmap
fills its own panel from the top (it does not leave the height of the tree's
toolbar and legend blank to line its rows up with the leaves). Alone, the
heatmap shows the tree order of the run.

Hover goes through the window's bus only: what is hovered in the panel goes
into the context with the panel's uid as origin, and other panels' hover comes
back through ``show_hover``. Clade, cluster and plasmid clicks narrow the
window's sample set (``set_subset``); a sample, a plasmid context or a subset
in a new explore window are open requests (``ctx.request``).
"""
from __future__ import annotations

from typing import Any, ClassVar, Dict, Optional, Sequence, Tuple

from PyQt6 import QtCore, QtWidgets

from ..model.filter_state import SampleSubset
from ..model.pan_data import VIEW_FORCE, VIEW_METADATA, VIEW_PLASMIDOME, VIEW_PRESENCE, VIEW_TREE
from ..session.requests import OpenPlasmidContext, OpenSampleCard, OpenSubsetInNewWindow
from .base import GENE_PRESENCE_ABSENCE, PAN_FORCE, PAN_METADATA, PHYLOGENY, PLASMIDOME, ExplorePanel


def pan_subset(payload: object, kind: str) -> Optional[SampleSubset]:
    """The sample subset a pan view's clade, cluster or plasmid click asks
    for (``{"name": ..., "samples": [...]}``), or None when it names none.
    A payload with a ``label`` is labeled with it as it is ("K locus = KL2");
    otherwise the label is the kind and the name ("Clade 2")."""
    if not isinstance(payload, dict):
        return None
    name = str(payload.get("name") or "").strip()
    samples = payload.get("samples")
    if not name or not isinstance(samples, (list, tuple)):
        return None
    label = str(payload.get("label") or "").strip() or f"{kind.title()} {name}"
    subset = SampleSubset(tuple(str(s) for s in samples), label=label, kind=kind)
    return subset if subset.samples else None


def request_pan_subset(panel: ExplorePanel, payload: object, kind: str) -> None:
    """Narrow the panel's window to the subset a pan view asked for (plan
    4.4 ``set_subset``, origin the panel's uid)."""
    subset = pan_subset(payload, kind)
    if subset is not None and panel.ctx is not None:
        panel.ctx.set_subset(subset, origin=panel.uid or "")


def request_pan_subset_window(panel: ExplorePanel, payload: object, kind: str) -> None:
    """Open the clade, cluster or plasmid a pan view asked for in a new explore
    window (plan 3.10); the panel's own window does not change."""
    subset = pan_subset(payload, kind)
    if subset is not None and panel.ctx is not None:
        panel.ctx.request(OpenSubsetInNewWindow(subset), origin=panel.uid or "")


def show_pan_hover(tab: Any, samples: Sequence[str]) -> None:
    """Highlight hovered samples in a PanCoreAccessoryTab (visible views only)."""
    samples = [str(s) for s in samples if str(s).strip()]
    if len(samples) > 1:
        tab.highlight_samples(samples)
    else:
        tab.highlight_sample(samples[0] if samples else None)


def hover_kind(samples: Sequence[str]) -> str:
    return "sample" if len(samples) == 1 else "samples"


def hovered_samples(samples: object) -> list:
    """A pan view's hover payload as a clean list (empty for none)."""
    if not isinstance(samples, (list, tuple)):
        return []
    return [str(s).strip() for s in samples if str(s).strip()]


class PanTabRequests:
    """The requests of a PanCoreAccessoryTab's views, for the panels that
    hold one: subsets of the window, subsets in a new window, plasmid context
    and sample cards."""

    def connect_pan_requests(self, tab: Any) -> None:
        tab.cladeFilterRequested.connect(self._on_clade_filter_requested)
        tab.clusterFilterRequested.connect(self._on_cluster_filter_requested)
        tab.plasmidFilterRequested.connect(self._on_plasmid_filter_requested)
        tab.subsetWindowRequested.connect(self._on_subset_window_requested)
        tab.plasmidContextRequested.connect(self._on_plasmid_context_requested)
        tab.sampleOpenRequested.connect(self._on_sample_open_requested)

    def _on_clade_filter_requested(self, payload: object) -> None:
        request_pan_subset(self, payload, "clade")  # type: ignore[arg-type]

    def _on_cluster_filter_requested(self, payload: object) -> None:
        request_pan_subset(self, payload, "cluster")  # type: ignore[arg-type]

    def _on_plasmid_filter_requested(self, payload: object) -> None:
        request_pan_subset(self, payload, "plasmid")  # type: ignore[arg-type]

    def _on_subset_window_requested(self, payload: object, kind: str) -> None:
        request_pan_subset_window(self, payload, kind)  # type: ignore[arg-type]

    def _on_plasmid_context_requested(self, payload: object) -> None:
        ctx = getattr(self, "ctx", None)
        if ctx is not None and isinstance(payload, dict):
            ctx.request(OpenPlasmidContext.from_payload(payload), origin=getattr(self, "uid", ""))

    def _on_sample_open_requested(self, sample: object) -> None:
        value = str(sample).strip() if sample is not None else ""
        ctx = getattr(self, "ctx", None)
        if value and ctx is not None:
            ctx.request(OpenSampleCard(value), origin=getattr(self, "uid", ""))


class PanViewPanel(PanTabRequests, ExplorePanel):
    """One pan view of the window's sample set (a single view tab)."""

    HEAVY = True
    # The PanCoreAccessoryTab single view this panel shows.
    VIEW: ClassVar[str] = ""
    # Its palette group.
    GROUP: ClassVar[str] = "Pangenome"

    # Shown, hidden, built or closed for good: the window's row link
    # (panels/row_link.py) checks which panels it links again.
    linkStateChanged = QtCore.pyqtSignal()

    def __init__(self, window: Any, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setObjectName(type(self).__name__)
        self._window = window
        self.tab: Any = None
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

    def set_tab(self, tab: Any) -> None:
        """The panel's single view tab (built by the window's factory as this
        panel's child)."""
        self.tab = tab
        self.layout().addWidget(tab)
        tab.samplesHovered.connect(self._on_samples_hovered)
        self.connect_pan_requests(tab)

    def view(self) -> Optional[QtWidgets.QWidget]:
        """The view widget this panel shows (tree, heatmap, donuts, graph or plasmid table)."""
        return self.tab.single_view_widget() if self.tab is not None else None

    def apply_result(self, filters: Any, result: Any) -> None:
        del filters
        self.tab.set_sample_filter(list(result.samples))

    def show_hover(self, hover: Any) -> None:
        show_pan_hover(self.tab, tuple(getattr(hover, "samples", ()) or ()))

    def on_shown(self) -> None:
        # A force graph whose data changed while hidden is built when shown.
        if self.tab is not None:
            self.tab.ensure_views_fresh()
        self.linkStateChanged.emit()

    def on_hidden(self) -> None:
        self.linkStateChanged.emit()

    def dispose(self) -> None:
        self.linkStateChanged.emit()

    def _on_samples_hovered(self, samples: object) -> None:
        values = hovered_samples(samples)
        if self.ctx is not None:
            self.ctx.hover(values, origin=self.uid or None, kind=hover_kind(values))


class PhylogenyPanel(PanViewPanel):
    """The run's tree, grouped by fastGEAR clade, with a clade and cluster
    legend. It leads the rows of the window's Gene presence/absence panels."""

    KIND = PHYLOGENY
    TITLE = "Phylogeny"
    VIEW = VIEW_TREE
    GROUP = "Phylogeny"

    # The tree's drawn rows (leaf order, layout) changed; its scroll and
    # zoom alone do not count (the heatmap never follows them).
    rowsViewChanged = QtCore.pyqtSignal()

    def set_tab(self, tab: Any) -> None:
        super().set_tab(tab)
        self._rows_view_key: Optional[Tuple[Any, ...]] = None
        tree = self.tree()
        if tree is not None:
            tree.viewChanged.connect(self._on_tree_view_changed)
            # Zoomed in: "Rows 1 to 60 of 153" says the rest is scrolled away.
            tree.set_row_note_enabled(True)
        source_changed = getattr(tab, "treeSourceChanged", None)
        if source_changed is not None:
            source_changed.connect(self.titleDetailChanged)

    def title_detail(self) -> str:
        """The tree drawn ("Core genome", "Snippy SNPs")."""
        current = getattr(self.tab, "current_tree_source", None)
        source = current() if callable(current) else None
        return source.short_label if source is not None else ""

    def tree(self) -> Any:
        return self.view()

    def _on_tree_view_changed(self) -> None:
        tree = self.tree()
        if tree is None:
            return
        key = (tree.isVisible(), tree.layout_mode(), tuple(tree.leaf_order()))
        if key != self._rows_view_key:
            self._rows_view_key = key
            self.rowsViewChanged.emit()

    # Display options (panel params, plan 3.11): tree layout, branch scale,
    # the Settings column's fold once picked, the tree source and root when
    # picked (an id the run does not list draws the run's default tree) and
    # the collapsed subtrees (by leaf names; the ones the drawn tree has no
    # subtree for are dropped). The source, root and collapsed subtrees stay
    # in the run's window slots; layouts leave them out
    # (layout_store.RUN_DATA_PARAMS).
    def save_state(self) -> Dict[str, Any]:
        return self.tab.tree_display_options() if self.tab is not None else {}

    def restore_state(self, state: Dict[str, Any]) -> None:
        if self.tab is not None and state:
            self.tab.set_tree_display_options(dict(state))

    # ------------------------------------------------------------ row link
    def row_order(self) -> Optional[Tuple[str, ...]]:
        """The tree's drawn leaves, top to bottom; None while it draws no
        rows (no tree, every genome filtered out, the circular layout)."""
        tree = self.tree()
        if tree is None or not tree.isVisible() or tree.layout_mode() != "rectangular":
            return None
        order = tuple(str(s) for s in tree.leaf_order())
        return order or None


class GenePresenceAbsencePanel(PanViewPanel):
    """The Panaroo gene presence/absence matrix, every row always fitted top
    to bottom. In a window with a Phylogeny panel its rows are the tree's
    leaves in the tree's order (panels/row_link.py); the tree's scroll and
    zoom do not move them."""

    KIND = GENE_PRESENCE_ABSENCE
    TITLE = "Gene presence/absence"
    VIEW = VIEW_PRESENCE

    def heatmap(self) -> Any:
        return self.view()

    def is_linked(self) -> bool:
        heatmap = self.heatmap()
        return heatmap is not None and heatmap.linked_rows() is not None

    def follow_rows(self, order: Optional[Sequence[str]]) -> None:
        """Show these rows (a tree's leaves, top to bottom), all fitted to
        the height (None: the heatmap's own rows)."""
        heatmap = self.heatmap()
        if heatmap is not None:
            heatmap.set_linked_rows(list(order) if order is not None else None)


class PanMetadataPanel(PanViewPanel):
    KIND = PAN_METADATA
    TITLE = "Metadata summary"
    VIEW = VIEW_METADATA


class PanForcePanel(PanViewPanel):
    """The cgMLST minimum spanning tree without its Settings column (the
    fastGEAR coloring of the pan views)."""

    KIND = PAN_FORCE
    TITLE = "cgMLST MST (compact)"
    # Beside the cgMLST MST in the palette.
    GROUP = "Typing"
    VIEW = VIEW_FORCE


class PlasmidomePanel(PanViewPanel):
    KIND = PLASMIDOME
    TITLE = "Plasmidome"
    VIEW = VIEW_PLASMIDOME


# Panel class of each pan kind (in palette order).
PAN_VIEW_PANELS = (PhylogenyPanel, GenePresenceAbsencePanel, PanMetadataPanel, PanForcePanel, PlasmidomePanel)
PAN_VIEW_PANEL_BY_KIND = {cls.KIND: cls for cls in PAN_VIEW_PANELS}
