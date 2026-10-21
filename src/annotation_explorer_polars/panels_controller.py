from __future__ import annotations

from typing import Callable, Dict, List, Optional, Sequence, Tuple

import pandas as pd
from PyQt6 import QtWidgets

from .docking.dock_host import Placement
from .cgmlst_force_plot import CgmlstForcePlotPanel
from .model.query_engine import COUNTRY_COLUMN_CANDIDATES, resolve_metadata_column
from .plasmid_model import plain_plasmid_name
from .panaroo_core_accessory import BRANCH_COLOR_CLUSTER
from .model.pan_data import VIEW_PRESENCE
from .pan_clustering_view import PanClusteringView
from .panaroo_views import PanCoreAccessoryTab
from .panels.alignment import AlignmentDocument, AlignmentPanel
from .panels.base import (
    ALIGNMENT,
    ANNOTATIONS,
    NEIGHBOURHOOD,
    SAMPLE_TREE,
    ExplorePanel,
    PanelRegistry,
)
from .panels.cgmlst_force import CgmlstForcePanel
from .panels.geo import GeoPanel
from .panels.neighbourhood import NeighbourhoodDocument, NeighbourhoodPanel
from .panels.pan import PAN_VIEW_PANEL_BY_KIND
from .panels.pan_clustering import PanClusteringPanel
from .session.requests import (
    ClosePanel,
    FocusPanel,
    OpenGeneCard,
    OpenPanel,
    OpenPlasmidContext,
    OpenSampleCard,
    OpenSubsetInNewWindow,
    ShowAnnotationsTable,
)
from .shell.commands import WINDOW_OPEN_NEIGHBOURHOOD, WINDOW_OPEN_SUBSET


class PanelsControllerMixin:
    """The window side of the panels: their factories, open
    requests and the pan views' subset requests (plan 4.8, Phase 4).

    Delivery (results only to visible heavy panels, replay on show, lazy
    builds) is the panel frames' job (panels/base.py); nothing here decides
    by tab names any more."""

    _COUNTRY_TOOLTIP_CANDIDATES: Tuple[str, ...] = COUNTRY_COLUMN_CANDIDATES
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

    # ------------------------------------------------------------ factories
    # One factory per kind (plan 4.8). ``uid`` is the instance's uid; the
    # first instance of a kind has the kind as uid.

    def _create_geo_panel(self, uid: str, parent: Optional[QtWidgets.QWidget] = None) -> Optional[ExplorePanel]:
        del uid
        return GeoPanel(self, parent)

    def _create_pan_tab(self, parent: QtWidgets.QWidget, *, single_view: str) -> Optional[PanCoreAccessoryTab]:
        """The one factory of the pan views: the one view of a pan panel.
        Every tab shares the run's PanData (plan 8.10) and draws the tree the
        same way, so a Gene presence/absence panel alone orders its rows like
        the Phylogeny panel: the run's default tree (the recombination-filtered
        core tree of ``trees/trees.json``; PanData.tree_index), rooted as its
        file says (an old unrooted file at its midpoint; the tree toolbar's
        Root dropdown picks another root, its source dropdown another tree)
        with the fastGEAR clades rotated together, genomes outside the sample set greyed (as neoPan drew them; a
        linked heatmap dims their rows), cluster colored branches with neutral
        labels, cluster brackets and the clade and cluster legend."""
        if self.run.output_dir is None:
            return None
        return PanCoreAccessoryTab(
            self.run.output_dir,
            self.run.meta_df,
            reference_sample=self.run.snippy_reference_sample,
            ann_df=self.run.ann_df,
            preserve_branch_lengths=True,
            color_tree_by_fastgear_clade=True,
            optimize_tree_by_fastgear_clade=True,
            annotate_fastgear_clusters=True,
            use_fastgear_clade_backgrounds=True,
            reroot_tree_for_clades=False,
            branch_color_mode=BRANCH_COLOR_CLUSTER,
            branch_colors_only=True,
            show_clade_legend=True,
            single_view=single_view,
            pan_data=self.run.pan,
            parent=parent,
        )

    def _pan_view_available(self, view: str) -> Tuple[bool, str]:
        """Whether this run has the data of a standalone pan panel (plan 3.6)."""
        if self.run.output_dir is None:
            return False, "Not available: this run has no output folder."
        return self.run.pan.availability(view)

    def _create_pan_view_panel(self, uid: str, parent: Optional[QtWidgets.QWidget] = None) -> Optional[ExplorePanel]:
        """A pan panel (Phylogeny, Gene presence/absence, Metadata summary,
        cgMLST MST (compact) or Plasmidome; plan 8.10): one pan view over the run's PanData.
        The window's row link (panels/row_link.py) finds the Phylogeny and
        Gene presence/absence panels when they are built."""
        panel_class = PAN_VIEW_PANEL_BY_KIND.get(PanelRegistry.kind_of(uid))
        if panel_class is None or not self._pan_view_available(panel_class.VIEW)[0]:
            return None
        panel = panel_class(self, parent)
        tab = self._create_pan_tab(panel, single_view=panel_class.VIEW)
        if tab is None:
            panel.deleteLater()
            return None
        panel.set_tab(tab)
        return panel

    def _create_pan_mds_panel(self, uid: str, parent: Optional[QtWidgets.QWidget] = None) -> Optional[ExplorePanel]:
        """A Pangenome clustering panel over the run's Panaroo matrix (PanData),
        colored by the cgMLST MST's traits (metadata, ST, Clade, Cluster)."""
        del uid
        if not self._pan_view_available(VIEW_PRESENCE)[0]:
            return None
        panel = PanClusteringPanel(self, parent)
        panel.set_view(
            PanClusteringView(
                self.run.pan,
                traits=self._build_cgmlst_trait_values(),
                hover_metadata=self.run.pan.sample_hover_metadata(),
                parent=panel,
            )
        )
        return panel

    def _create_cgmlst_force_panel(self, uid: str, parent: Optional[QtWidgets.QWidget] = None) -> Optional[ExplorePanel]:
        del uid
        if self.run.output_dir is None:
            return None
        panel = CgmlstForcePanel(self, parent)
        plot = self._build_cgmlst_force_content(panel)
        if plot is None:
            panel.deleteLater()
            return None
        panel.set_plot(plot)
        return panel

    def _create_neighbourhood_panel(self, uid: str, parent: Optional[QtWidgets.QWidget] = None) -> Optional[ExplorePanel]:
        """A document panel: built from the comparison it was opened with."""
        document = self._panels.document(uid)
        if not isinstance(document, NeighbourhoodDocument):
            return None
        # The run folder is where figures are exported by default.
        view_options = {"output_dir": self.run.output_dir}
        return NeighbourhoodPanel(document, view_options=view_options, parent=parent)

    def _create_alignment_panel(self, uid: str, parent: Optional[QtWidgets.QWidget] = None) -> Optional[ExplorePanel]:
        document = self._panels.document(uid)
        if not isinstance(document, AlignmentDocument):
            return None
        return AlignmentPanel(document, parent)

    # ------------------------------------------------------- open requests
    def _on_open_request(self, request: object, origin: str = "") -> None:
        """Carry out a panel's open request (``ctx.request``) in this window."""
        if isinstance(request, OpenSampleCard):
            self._open_sample_card_for_sample(request.sample)
        elif isinstance(request, OpenGeneCard):
            self._open_gene_card(request.row_dict())
        elif isinstance(request, OpenPlasmidContext):
            self._on_pan_plasmid_context_requested(request.payload_dict(), origin=origin)
        elif isinstance(request, FocusPanel):
            self._focus_panel(request.kind, request.params_dict())
        elif isinstance(request, OpenPanel):
            self._open_panel_request(request)
        elif isinstance(request, ClosePanel):
            self._close_panel(request.uid)
        elif isinstance(request, OpenSubsetInNewWindow):
            # A new, independent explore window; this one does not change.
            self.app_controller.commands.run(WINDOW_OPEN_SUBSET, self, subset=request.subset)
        elif isinstance(request, ShowAnnotationsTable):
            self._show_annotations_table(request.beside)

    def _show_annotations_table(self, beside: Optional[str] = None) -> None:
        """A neighbourhood panel's "Show only this sample in the window":
        the window is narrowed already. A window with no Annotations panel
        open (a neighbourhood window) also gets the Sample tree as the left
        column and Annotations below the panel, brought forward. A window
        whose Annotations panel is open keeps its current tabs (the menu
        item promises a narrowed window, not a switch away from the
        comparison) but gets the Sample tree when it is missing, since
        Annotations shows nothing without it."""
        manager = self._panel_manager
        if not manager.is_open(SAMPLE_TREE):
            manager.open_panel(SAMPLE_TREE, placement=Placement("left", outer=True))
        if any(manager.is_open(uid) for uid in manager.uids(ANNOTATIONS)):
            return
        anchor = beside if beside and manager.is_open(beside) else None
        manager.open_panel(ANNOTATIONS, placement=Placement("bottom", relative_to=anchor))
        if manager.is_open(ANNOTATIONS):
            manager.host.focus_view(ANNOTATIONS)

    def _open_panel_request(self, request: OpenPanel) -> Optional[str]:
        """Open a panel of ``request.kind`` (a new instance of a multi-instance
        kind that is open already) as a tab beside ``request.beside``, with
        ``request.params`` as its display options. Returns the uid opened."""
        manager = self._panel_manager
        kind = request.kind
        params = request.params_dict()
        beside = request.beside if request.beside and manager.is_open(request.beside) else None
        if manager.is_open(kind) and manager.registry.is_multi_instance(kind):
            # A new instance, drawn once with its options.
            placement = Placement("tab", relative_to=beside) if beside else None
            return manager.new_instance(kind, placement=placement, state=params or None)
        uid = manager.open_new_instance(kind, target=beside)
        frame = manager.frame(uid) if uid else None
        if frame is not None and params:
            frame.restore_state(params)
        return uid

    def _focus_panel(self, kind: str, params: Dict[str, object]) -> None:
        """Show the instance of ``kind`` used last and apply ``params`` to it
        (the tree's "Open in Context view": the Annotations panel's mode)."""
        manager = self._panel_manager
        uid = manager.focus_panel(kind)
        panel = manager.panel(uid) if uid else None
        mode = str(params.get("view_mode") or "")
        if kind == ANNOTATIONS and mode:
            if panel is None:
                panel = self.annotations_panel
            if panel is not None:
                # Emits viewModeChanged: the panel (and, for the paired one,
                # the track) redraws from the applied result.
                panel.set_view_mode(mode)

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

    _resolve_plot_metadata_column = staticmethod(resolve_metadata_column)

    def _load_plot_mlst_sequence_types(self) -> Dict[str, str]:
        if self.run.output_dir is None:
            return {}
        path = self.run.output_dir / "mlst" / "results.csv"
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
        if not self.run.meta_df.empty and "Sample" in self.run.meta_df.columns:
            columns = [str(column) for column in self.run.meta_df.columns]
            country_col = self._resolve_plot_metadata_column(columns, self._COUNTRY_TOOLTIP_CANDIDATES)
            st_col = self._resolve_plot_metadata_column(columns, self._ST_TOOLTIP_CANDIDATES)
            date_col = self._resolve_plot_metadata_column(columns, self._COLLECTION_DATE_TOOLTIP_CANDIDATES)
            for _idx, row in self.run.meta_df.iterrows():
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
        # The run's fastGEAR clades and clusters, read once per run (PanData).
        lineages, clusters = self.run.pan.fastgear_assignments()
        return dict(lineages), dict(clusters)

    def _build_cgmlst_trait_values(self) -> Dict[str, Dict[str, str]]:
        traits: Dict[str, Dict[str, str]] = {}
        if not self.run.meta_df.empty and "Sample" in self.run.meta_df.columns:
            for column in self.run.meta_df.columns:
                if str(column) == "Sample":
                    continue
                values: Dict[str, str] = {}
                for sample_value, raw_value in self.run.meta_df[["Sample", column]].itertuples(index=False):
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

    def _build_cgmlst_force_content(self, parent: QtWidgets.QWidget) -> Optional[CgmlstForcePlotPanel]:
        if self.run.output_dir is None:
            return None
        panel = CgmlstForcePlotPanel(
            self.run.output_dir,
            reference_sample=self.run.snippy_reference_sample,
            parent=parent,
        )
        panel.set_trait_values(self._build_cgmlst_trait_values())
        return panel

    def _refresh_cgmlst_force_plot(self, plot: CgmlstForcePlotPanel, samples: Sequence[str]) -> None:
        """Bring the cgMLST MST up to date for the sample set."""
        lineage_by_sample, cluster_by_sample = self._load_cgmlst_fastgear_maps()
        plot.set_reference_sample(self.run.snippy_reference_sample, refresh=False)
        plot.set_fastgear_assignments(lineage_by_sample, {}, cluster_by_sample)
        plot.set_sample_hover_metadata(self._build_plot_sample_hover_metadata())
        plot.set_trait_values(self._build_cgmlst_trait_values())
        plot.set_sample_filter(list(samples))
        plot.ensure_fresh()

    def _on_pan_plasmid_context_requested(self, payload: object, origin: str = "") -> None:
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
        name = plain_plasmid_name(payload.get("name")) or "Plasmid"
        # One panel per plasmid cluster and set of instances.
        members = tuple(sorted({(str(row.get("sample") or ""), str(row.get("contig") or "")) for row in rows if isinstance(row, dict)}))
        # In a new neighbourhood window with this window's filters (``origin``,
        # a panel of this window, does not apply there).
        del origin
        self.app_controller.commands.run(
            WINDOW_OPEN_NEIGHBOURHOOD,
            self,
            rows=rows,
            title=f"Plasmid context: {name}",
            source="plasmid",
            key=("plasmid", str(payload.get("name") or ""), members),
        )

    # The neighbourhood comparisons and alignment results are document panels
    # (plan 8.7). A comparison opens by default in a neighbourhood window of
    # its own (``window.open_neighbourhood``), which places it with
    # ``panel.open_neighbourhood``; an alignment opens beside its source.

    def _open_neighbourhood_view(
        self,
        *,
        title: str,
        rows: Optional[List[Dict[str, object]]] = None,
        loader: Optional[Callable[[], List[Dict[str, object]]]] = None,
        source: str = "plasmid",
        key: object = None,
        flank_sizes: Optional[Sequence[int]] = None,
        initial_flank: Optional[int] = None,
        beside: Optional[str] = None,
        reference_sample: str = "",
    ) -> Optional[str]:
        """Command ``panel.open_neighbourhood`` (in this window): the plasmid
        context of a pan view (``source="plasmid"``, with its ``rows``) or a
        gene card's neighborhoods (``source="gene_card"``, with a ``loader``
        that collects them while the panel shows its busy state). Opening a
        comparison that is open in this window (same ``key``) shows it.
        Returns the panel's uid."""
        document = NeighbourhoodDocument(
            title=title,
            rows=[dict(row) for row in rows] if rows is not None else None,
            loader=loader,
            source=source,
            key=key,
            flank_sizes=tuple(flank_sizes) if flank_sizes is not None else None,
            initial_flank=initial_flank,
        )
        if reference_sample and self._open_neighbourhood_tab_titles(document.key).count(document.tab_title()):
            # Two cards of one gene in different samples: tell the tabs apart.
            document.tab_suffix = str(reference_sample)
        return self._open_document_panel(NEIGHBOURHOOD, document, beside=beside)

    def _open_neighbourhood_tab_titles(self, key: object) -> List[str]:
        """Tab titles of the other open neighbourhood comparisons, in every
        explore window of the run (each comparison usually has its own
        neighbourhood window, so their titles must differ too)."""
        windows = [w for w in self.app_controller.windows.explore_windows() if w is not self]
        titles: List[str] = []
        for window in [self] + windows:
            panels = getattr(window, "_panels", None)
            if panels is None:
                continue
            for uid in panels.uids(NEIGHBOURHOOD):
                document = panels.document(uid)
                if document is not None and (key is None or getattr(document, "key", None) != key):
                    titles.append(document.tab_title())
        return titles

    def _open_alignment_view(self, *, alignment_text: str, status_text: str, title: str = "Sequence alignment") -> Optional[str]:
        """Command ``panel.open_alignment``: a gene card's alignment result, as
        a panel of this window beside its active group."""
        document = AlignmentDocument(alignment_text=alignment_text, status_text=status_text, title=title)
        return self._open_document_panel(ALIGNMENT, document)

    def _open_document_panel(self, kind: str, document: object, *, beside: Optional[str] = None) -> Optional[str]:
        return self._panel_manager.open_document(kind, document, beside=beside)

    def _close_panel(self, uid: str) -> None:
        """A panel's own Ctrl+W (``ClosePanel``): close only that panel."""
        self._panel_manager.close_panel(uid)

