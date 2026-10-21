"""Standalone pan panels over a shared PanData (window plan 8.10, Phase 9).

- ``PanData`` (``run.pan``): which pan outputs a run has, the folder lookups
  done once, the shared fastGEAR, hover and plasmid data, ``memo``; no Qt.
- The five panels (Phylogeny, Gene presence/absence, Metadata summary, cgMLST MST (compact),
  Plasmidome) in the registry and the palette, disabled with a reason when the
  run lacks their data; each shows its one view over the window's sample set.
  neoPan, the Pan tree and Gene presence are gone (feedback 2026-10-07).
- Hover only through the window's bus: the panels' hover goes out with their
  uid and comes back from other panels, never echoed; PanCoreAccessoryTab has
  no hover hub of its own, so the Phylogeny tree and the Gene presence/absence
  heatmap meet on the bus like any two panels.
- Clade, cluster and plasmid clicks narrow the window's sample set, the other
  requests go to the window.
- One PanData per run, shared by every pan panel of every window.
- The Pangenome preset: Phylogeny beside Gene presence/absence, the map beside
  Plasmidome; a layout saved before Phase 9 restores with them closed.

The linked rows of Phylogeny and Gene presence/absence are tested in
``tests/test_row_link.py``.

Offscreen on the synthetic run of the golden tests with the pan outputs of
``tests/_pan_fixture.py``; the window tests are skipped without PyQt6Ads.
"""
from __future__ import annotations

import os
import shutil
import tempfile
import time
import unittest
from pathlib import Path
from types import SimpleNamespace
from typing import Any, List
from unittest import mock

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

import pandas as pd
from PyQt6 import QtCore, QtWidgets

try:
    from tests._pan_fixture import CLADES, write_pan_outputs
    from tests.test_dock_host import _app, _Explorer, _fake_specs, _sorted_groups, _spin
    from tests.test_panel_contract import _WindowCase
    from tests.test_refresh_parity import build_synthetic_run
    from tests._timing import within_budget
except ImportError:
    from _pan_fixture import CLADES, write_pan_outputs
    from test_dock_host import _app, _Explorer, _fake_specs, _sorted_groups, _spin
    from test_panel_contract import _WindowCase
    from test_refresh_parity import build_synthetic_run
    from _timing import within_budget

from annotation_explorer_polars.docking import docking_available
from annotation_explorer_polars.model import RunData
from annotation_explorer_polars.model import pan_data as pan_data_mod
from annotation_explorer_polars.model.pan_data import PanData
from annotation_explorer_polars.session import bus as bus_events
from annotation_explorer_polars.session.bus import Delivery

PAN_KINDS = ("phylogeny", "gene_presence_absence", "pan_metadata", "pan_force", "plasmidome")

if docking_available():
    from annotation_explorer_polars.docking import layouts
    from annotation_explorer_polars.docking.dock_host import Placement
    from annotation_explorer_polars.panaroo_core_accessory import LENGTHS_TOPOLOGY, PanCoreAccessoryTab
    from annotation_explorer_polars.panels.base import PanelSpec
    from annotation_explorer_polars.panels.manager import PanelManager
    from annotation_explorer_polars.panels.pan import PanViewPanel
    from annotation_explorer_polars.session.requests import (
        OpenPlasmidContext,
        OpenSampleCard,
        OpenSubsetInNewWindow,
    )
    from annotation_explorer_polars.shell.app_controller import AppController


def _run(output_dir: Any) -> RunData:
    ann_df, meta_df = build_synthetic_run()
    return RunData.build(
        ann_df,
        meta_df,
        ref_df=pd.DataFrame(),
        marker_hits_df=pd.DataFrame(),
        output_dir=Path(output_dir) if output_dir is not None else None,
    )


class PanDataTests(unittest.TestCase):
    """PanData alone (no window, no Qt needed)."""

    def setUp(self) -> None:
        self.pan_dir = tempfile.mkdtemp(prefix="bactolomew_pan_run_")
        self.empty_dir = tempfile.mkdtemp(prefix="bactolomew_pan_empty_")
        write_pan_outputs(self.pan_dir)

    def tearDown(self) -> None:
        shutil.rmtree(self.pan_dir, ignore_errors=True)
        shutil.rmtree(self.empty_dir, ignore_errors=True)

    def test_the_run_carries_one_pan_data(self) -> None:
        run = _run(self.pan_dir)
        self.assertIsInstance(run.pan, PanData)
        self.assertEqual(run.pan.output_dir, Path(self.pan_dir))
        # The run's own frozen plasmid index and metadata, not copies.
        self.assertIs(run.pan.plasmid_cluster_index, run.plasmid_cluster_index)
        self.assertIs(run.pan.meta_df, run.meta_df)
        self.assertEqual(run.pan.reference_sample, run.snippy_reference_sample)

    def test_availability_follows_the_run_outputs(self) -> None:
        full = _run(self.pan_dir).pan
        for view in pan_data_mod.SINGLE_VIEWS:
            self.assertEqual(full.availability(view), (True, ""), view)
        empty = _run(self.empty_dir).pan
        reasons = {view: empty.availability(view) for view in pan_data_mod.SINGLE_VIEWS}
        self.assertEqual(reasons["tree"], (False, "Not available: this run has no phylogenetic tree."))
        self.assertEqual(reasons["presence"], (False, "Not available: this run has no Panaroo output."))
        self.assertEqual(reasons["force"], (False, "Not available: this run has no chewBBACA (cgMLST) output."))
        # The synthetic run's metadata and plasmid annotations are there anyway.
        self.assertEqual(reasons["metadata"], (True, ""))
        self.assertEqual(reasons["plasmidome"], (True, ""))
        no_folder = _run(None).pan
        for view in pan_data_mod.SINGLE_VIEWS:
            self.assertEqual(no_folder.availability(view), (False, "Not available: this run has no output folder."))
        no_plasmids = PanData(Path(self.empty_dir), pd.DataFrame({"Sample": ["S1"]}))
        self.assertEqual(
            no_plasmids.availability("plasmidome"),
            (False, "Not available: this run has no MOB-suite plasmid clusters."),
        )
        self.assertEqual(
            no_plasmids.availability("metadata"),
            (False, "Not available: this run has no sample metadata."),
        )
        self.assertFalse(full.availability("nonsense")[0])

    def test_files_and_tables(self) -> None:
        pan = _run(self.pan_dir).pan
        self.assertEqual(pan.tree_path(), Path(self.pan_dir) / "trees" / "core_gene_tree.nwk")
        self.assertEqual(pan.panaroo_out(), Path(self.pan_dir) / "panaroo" / "out")
        self.assertEqual(pan.chewbbaca_out(), Path(self.pan_dir) / "chewbbaca" / "out")
        lineages, clusters = pan.fastgear_assignments()
        self.assertEqual(dict(lineages), {s: lineage for s, (lineage, _c) in CLADES.items()})
        self.assertEqual(dict(clusters), {s: cluster for s, (_l, cluster) in CLADES.items()})
        hover = pan.sample_hover_metadata()
        self.assertEqual(hover["S1"]["Country"], "Romania")
        self.assertNotIn("S5", hover)  # no country, ST or date for S5
        hover["S1"]["Country"] = "changed"
        self.assertEqual(pan.sample_hover_metadata()["S1"]["Country"], "Romania")  # a copy each time

    def test_folder_lookups_happen_once_per_run(self) -> None:
        pan = _run(self.pan_dir).pan
        with mock.patch.object(pan_data_mod, "tree_sources", wraps=pan_data_mod.tree_sources) as spy:
            for _ in range(3):
                pan.tree_path()
                pan.availability("tree")
        self.assertEqual(spy.call_count, 1)

    def test_resolved_lookups_leave_no_folder_search_for_availability(self) -> None:
        # The loader resolves the lookups on its worker thread; the palette's
        # availability checks afterwards search nothing.
        pan = _run(self.pan_dir).pan
        pan.resolve_lookups()
        with mock.patch.object(pan_data_mod, "tree_sources") as raxml, mock.patch.object(
            pan_data_mod, "find_panaroo_tree"
        ) as panaroo, mock.patch.object(pan_data_mod, "_find_panaroo_out") as panaroo_out, mock.patch.object(
            pan_data_mod, "_find_chewbbaca_out"
        ) as chewbbaca:
            for view in pan_data_mod.SINGLE_VIEWS:
                self.assertEqual(pan.availability(view), (True, ""), view)
        for spy in (raxml, panaroo, panaroo_out, chewbbaca):
            spy.assert_not_called()

    def test_fastgear_tables_are_read_once_per_run(self) -> None:
        pan = _run(self.pan_dir).pan
        with mock.patch.object(pan_data_mod, "_load_fastgear_lineages", wraps=pan_data_mod._load_fastgear_lineages) as spy:
            for _ in range(3):
                pan.fastgear_assignments()
        self.assertEqual(spy.call_count, 1)

    def test_memo_computes_once_per_key(self) -> None:
        pan = _run(self.pan_dir).pan
        calls: List[str] = []
        for _ in range(3):
            self.assertEqual(pan.memo(("k", 1), lambda: calls.append("a") or "A"), "A")
        self.assertEqual(pan.memo(("k", 2), lambda: calls.append("b") or "B"), "B")
        self.assertEqual(calls, ["a", "b"])


def _reference_cluster_gene_stats(ann_df: pd.DataFrame, cluster_key: str, samples: set) -> dict:
    """The Plasmidome counts as computed before Phase 9 (straight from the
    annotation table, per cluster), the reference for ``plasmid_gene_rows``."""
    from annotation_explorer_polars.plasmid_model import plasmid_cluster_key_from_row
    from annotation_explorer_polars.plasmidome_view import _gene_label, _truthy

    zero = {"unique": 0, "arg": 0, "mobile": 0}
    if ann_df.empty or "MobSuiteIsPlasmid" not in ann_df.columns:
        return zero
    df = ann_df.loc[ann_df["MobSuiteIsPlasmid"].map(_truthy)].copy()
    if samples:
        df = df.loc[df["Sample"].fillna("").astype(str).str.strip().isin(samples)].copy()
    if df.empty:
        return zero
    df["__cluster_key"] = [plasmid_cluster_key_from_row(r) for r in df.to_dict(orient="records")]
    df = df.loc[df["__cluster_key"].astype(str).str.strip() == cluster_key]
    if df.empty:
        return zero
    labels = {label for label in (_gene_label(r) for r in df.to_dict(orient="records")) if label}
    tags = df.get("Tags", pd.Series("", index=df.index)).fillna("").astype(str)
    return {
        "unique": len(labels),
        "arg": int(tags.str.contains("A", regex=False).sum()),
        "mobile": int(tags.str.contains("M", regex=False).sum()),
    }


class PlasmidomeRowsTests(unittest.TestCase):
    """The Plasmidome table from the shared ``plasmid_gene_rows`` equals the
    per-cluster computation it replaced, unfiltered and filtered."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.app = _app()

    def setUp(self) -> None:
        from annotation_explorer_polars.plasmid_model import build_plasmid_cluster_index

        self.pan_dir = tempfile.mkdtemp(prefix="bactolomew_plasmidome_")
        write_pan_outputs(self.pan_dir)
        ann_df, _meta = build_synthetic_run()
        extra = ann_df.iloc[:4].copy()
        # A secondary-only cluster, a padded sample name, no tags and an
        # unknown gene name, so the keys and labels take every path.
        extra["Sample"] = [" S5 ", "S5", "S6", "S6"]
        extra["MobSuiteIsPlasmid"] = ["yes", "True", "1", "True"]
        extra["MobSuiteSecondaryClusterID"] = ["BB009", "BB009", "BB009", ""]
        extra["MobSuiteMashNearestNeighbor"] = ["", "", "", "NN1"]
        extra["Gene"] = ["geneX", "UNK:1", "geneX", ""]
        extra["Tags"] = ["A", None, "AM", "M"]
        self.ann_df = pd.concat([ann_df, extra], ignore_index=True)
        self.index = build_plasmid_cluster_index(output_dir=Path(self.pan_dir), ann_df=self.ann_df).freeze()

    def tearDown(self) -> None:
        shutil.rmtree(self.pan_dir, ignore_errors=True)

    def _panel(self, **kwargs: Any) -> Any:
        from annotation_explorer_polars.plasmidome_view import PlasmidomePanel

        panel = PlasmidomePanel(self.ann_df, Path(self.pan_dir), self.index, **kwargs)
        self.addCleanup(panel.deleteLater)
        return panel

    def test_counts_equal_the_per_cluster_computation(self) -> None:
        from annotation_explorer_polars.plasmidome_view import plasmid_gene_rows

        rows = plasmid_gene_rows(self.ann_df)
        self.assertEqual(list(rows.columns), ["__sample", "__cluster_key", "__gene_label", "__tags"])
        keys = set(rows["__cluster_key"])
        self.assertGreaterEqual(len(keys), 3)
        panel = self._panel(gene_rows=rows)
        for sample_set in (set(), {"S1", "S2"}, {"S5", "S6"}, {"S4"}):
            for key in keys:
                self.assertEqual(
                    panel._cluster_gene_stats(key, sample_set),
                    _reference_cluster_gene_stats(self.ann_df, key, sample_set),
                    (key, sample_set),
                )

    def test_open_context_sends_every_instance_whatever_the_sample_filter(self) -> None:
        # The comparison opens in a neighbourhood window that narrows the rows
        # to its own sample set, which may be widened later.
        panel = self._panel()
        sent: List[Any] = []
        panel.contextRequested.connect(sent.append)
        found = {}
        for key in self.index.clusters:
            panel._emit_context_requested(key)
            samples = sorted({row["sample"] for row in sent[-1]["rows"]})
            if len(samples) > 1:
                found[key] = samples
        self.assertTrue(found)
        key, everyone = next(iter(found.items()))
        panel.set_sample_filter([everyone[0]])
        panel._emit_context_requested(key)
        self.assertEqual(sorted({row["sample"] for row in sent[-1]["rows"]}), everyone)

    def test_the_table_follows_the_sample_filter(self) -> None:
        def table(panel: Any) -> list:
            return [
                (row["cluster_key"], tuple(row["samples"]), row["unique_genes"], row["arg_genes"], row["mobile_genes"])
                for row in panel._rows
            ]

        shared = self._panel()
        own = self._panel(gene_rows=None)
        self.assertEqual(table(shared), table(own))
        everything = table(shared)
        self.assertTrue(everything)
        shared.set_sample_filter(["S1", "S4"])
        self.assertTrue(table(shared))
        for key, samples, *_counts in table(shared):
            self.assertTrue(set(samples) <= {"S1", "S4"}, key)
            stats = _reference_cluster_gene_stats(self.ann_df, key, set(samples))
            self.assertEqual(tuple(_counts), (stats["unique"], stats["arg"], stats["mobile"]))
        # An empty sample set lists no cluster; no filter lists them all again.
        shared.set_sample_filter([])
        self.assertEqual(table(shared), [])
        self.assertEqual(shared.summary_label.text(), "0 plasmid cluster(s)")
        shared.set_sample_filter(None)
        self.assertEqual(table(shared), everything)


class _PanWindowCase(_WindowCase):
    """A real window on the synthetic run with pangenome outputs."""

    def setUp(self) -> None:
        super().setUp()
        write_pan_outputs(self.output_dir)

    def _open(self, win: Any, kind: str, **kwargs: Any) -> Any:
        manager = win._panel_manager
        self.assertTrue(manager.open_panel(kind, **kwargs), kind)
        self.assertTrue(self._wait_for(lambda: manager.panel(kind) is not None and manager.is_visible(kind)), kind)
        _spin(self.app, 0.3)
        return manager.panel(kind)


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class PanPanelTests(_PanWindowCase):
    def test_the_registry_and_palette_list_the_five_panels(self) -> None:
        win = self._window()
        manager = win._panel_manager
        kinds = manager.registry.kinds()
        # The compact cgMLST MST is listed after the cgMLST MST (Typing).
        pan_kinds = [kind for kind in PAN_KINDS if kind != "pan_force"]
        self.assertEqual(kinds[kinds.index("geo") + 1 : kinds.index("geo") + 5], pan_kinds)
        self.assertEqual(kinds[kinds.index("cgmlst_force") + 1], "pan_force")
        for retired in ("neopan", "pan_tree", "gene_presence", "neotree"):
            self.assertNotIn(retired, kinds)
        titles = {}
        for kind in PAN_KINDS:
            spec = manager.registry.spec(kind)
            self.assertEqual(spec.group, {"phylogeny": "Phylogeny", "pan_force": "Typing"}.get(kind, "Pangenome"))
            self.assertTrue(spec.heavy)
            self.assertTrue(manager.registry.is_multi_instance(kind))
            self.assertEqual(manager.is_available(kind), (True, ""))
            titles[kind] = spec.title
            # Classic has the tree and its heatmap in the bottom row.
            self.assertEqual(manager.is_open(kind), kind in ("phylogeny", "gene_presence_absence"))
        self.assertEqual(
            titles,
            {
                "phylogeny": "Phylogeny",
                "gene_presence_absence": "Gene presence/absence",
                "pan_metadata": "Metadata summary",
                "pan_force": "cgMLST MST (compact)",
                "plasmidome": "Plasmidome",
            },
        )
        from annotation_explorer_polars.panels.palette import PanelPalette

        palette = PanelPalette(manager.specs(), is_open=manager.is_open)
        try:
            entries = {
                item.data(0, QtCore.Qt.ItemDataRole.UserRole + 1): item for item in palette._items
            }
            for kind in PAN_KINDS:
                self.assertIn(kind, entries)
                self.assertFalse(entries[kind].isDisabled(), kind)
        finally:
            palette.deleteLater()

    def test_each_panel_shows_its_one_view_of_the_sample_set(self) -> None:
        win = self._window()
        expected_views = {
            "phylogeny": "RaxmlTreeView",
            "gene_presence_absence": "PanHeatmapWidget",
            "pan_metadata": "_MetadataDonutPanel",
            "pan_force": "CgmlstForceGraphView",
            "plasmidome": "PlasmidomePanel",
        }
        for kind in PAN_KINDS:
            with self.subTest(kind=kind):
                panel = self._open(win, kind)
                self.assertIsInstance(panel, PanViewPanel)
                self.assertIs(panel.ctx, win.ctx)
                self.assertEqual(type(panel.view()).__name__, expected_views[kind])
                tab = panel.tab
                self.assertEqual(tab.single_view(), panel.VIEW)
                self.assertFalse(hasattr(tab, "_detail_view_combo"))  # no view selector
                self.assertEqual(len(tab._all_views()), 1)
                # Only the tree view has the tree toolbar and the legend.
                self.assertEqual(tab._tree_toolbar is not None, kind == "phylogeny")
                self.assertEqual(tab._legend is not None, kind == "phylogeny")
                self.assertEqual(tab._sample_filter, set(win.ctx.result.samples))
        # A subset narrows every visible pan panel to the same samples.
        win._set_sample_subset_filter(label="Three", samples=["S1", "S3", "S4"])
        self._idle(win)
        manager = win._panel_manager
        for kind in PAN_KINDS:
            manager.open_panel(kind)
            _spin(self.app, 0.4)
            self.assertEqual(manager.panel(kind).tab._sample_filter, {"S1", "S3", "S4"}, kind)
        heatmap = manager.panel("gene_presence_absence").view()
        # The tree greys the others and the linked heatmap dims their rows
        # (its rows are the tree's leaves).
        leaves = sorted(manager.panel("phylogeny").tab.tree.leaf_order())
        self.assertEqual(leaves, ["S1", "S2", "S3", "S4", "S5", "S6"])
        self.assertEqual(sorted(heatmap.display_samples()), leaves)
        self.assertEqual(sorted(heatmap.dimmed_rows()), ["S2", "S5", "S6"])

    def test_gene_presence_absence_alone_orders_its_rows_like_the_phylogeny(self) -> None:
        win = self._window()
        manager = win._panel_manager
        tree = self._open(win, "phylogeny").tab.tree
        leaves = [str(s) for s in tree.leaf_order()]
        manager.close_panel("phylogeny")
        _spin(self.app, 0.3)
        panel = self._open(win, "gene_presence_absence")
        self.assertFalse(panel.is_linked())
        self.assertIsNone(panel.tab.tree)  # no tree kept for the order
        # Rows in the tree's leaf order, fastGEAR clades grouped.
        self.assertEqual(panel.view().display_samples(), leaves)
        lineages = [CLADES[s][0] for s in leaves]
        self.assertEqual(lineages, sorted(lineages, key=lineages.index))

    def test_gene_presence_absence_rows_follow_a_snippy_tree_with_a_reference_leaf(self) -> None:
        # A snippy tree names the reference sample's leaf "Reference"; the
        # hidden tree behind Gene presence/absence must map it like the tree.
        write_pan_outputs(self.output_dir, snippy_reference="S6")
        win = self._window()
        manager = win._panel_manager
        self.assertEqual(win.run.snippy_reference_sample, "S6")
        tab = self._open(win, "phylogeny").tab
        leaves = [str(s) for s in tab.tree.leaf_order()]
        self.assertIn("S6", leaves)
        self.assertNotIn("Reference", leaves)
        panel = self._open(win, "gene_presence_absence")
        self.assertTrue(self._wait_for(panel.is_linked))
        self.assertEqual(panel.view().display_samples(), leaves)
        manager.close_panel("phylogeny")
        self.assertTrue(self._wait_for(lambda: not panel.is_linked()))
        self.assertEqual(panel.view().display_samples(), leaves)

    def test_hover_goes_out_with_the_panel_uid_and_is_never_echoed(self) -> None:
        win = self._window()
        panel = self._open(win, "phylogeny")
        self.assertTrue(win._panel_manager.is_visible("sample_tree"))
        tab = panel.tab
        with mock.patch.object(tab, "highlight_sample") as one, mock.patch.object(tab, "highlight_samples") as several:
            tab.tree.leafHovered.emit("S3")
            win.ctx.flush_hover()
            _spin(self.app, 0.15)
            self.assertEqual((win.ctx.hover_state.samples, win.ctx.hover_state.origin), (("S3",), "phylogeny"))
            self.assertEqual(win._tree_hover_shown, ("S3",))
            one.assert_not_called()
            several.assert_not_called()
            # A subtree hovers several samples.
            tab.tree.subtreeHovered.emit(["S3", "S4"])
            win.ctx.flush_hover()
            self.assertEqual(win.ctx.hover_state.samples, ("S3", "S4"))
            self.assertEqual(win.ctx.hover_state.kind, "samples")
            _spin(self.app, 0.15)
            several.assert_not_called()

    def test_hover_from_other_panels_reaches_every_visible_pan_panel(self) -> None:
        win = self._window()
        manager = win._panel_manager
        self._open(win, "phylogeny")
        self._open(win, "gene_presence_absence")
        self._open(win, "pan_force", to_side=True)
        for kind in ("phylogeny", "gene_presence_absence", "pan_force"):
            self.assertTrue(manager.is_visible(kind), kind)
        tabs = {kind: manager.panel(kind).tab for kind in ("phylogeny", "gene_presence_absence", "pan_force")}
        patches = [mock.patch.object(tab, "highlight_sample") for tab in tabs.values()]
        spies = dict(zip(tabs, [p.start() for p in patches]))
        try:
            win.ctx.hover(["S4"], origin="sample_tree")
            win.ctx.flush_hover()
            for kind, spy in spies.items():
                self.assertTrue(self._wait_for(lambda s=spy: s.called), kind)
                spy.assert_called_with("S4")
        finally:
            for p in patches:
                p.stop()
        # The real views draw it: the heatmap marks the row, the graph the node.
        win.ctx.hover(["S2"], origin="sample_tree")
        win.ctx.flush_hover()
        heatmap = manager.panel("gene_presence_absence").view()
        self.assertTrue(self._wait_for(lambda: heatmap._highlight_sample == "S2"))
        graph = manager.panel("pan_force").view()
        self.assertTrue(self._wait_for(lambda: graph._highlight_sample == "S2"))

    def test_an_external_hover_emphasizes_the_phylogeny_sample(self) -> None:
        # Feedback 2026-10-08: a hover from another panel gets the same
        # emphasis as one in the tree (thick amber path, dot and chip, the
        # rest dimmed), and clearing it clears everything.
        win = self._window()
        tree = self._open(win, "phylogeny").tab.tree
        self.assertTrue(self._wait_for(lambda: bool(tree.leaf_order())))
        win.ctx.hover(["S4"], origin="sample_tree")
        win.ctx.flush_hover()
        self.assertTrue(self._wait_for(tree.hover_dimmed))
        self.assertTrue(all(item.pen().isCosmetic() for item in tree._edge_items_by_leaf["S4"]))
        overlay = tree.hover_overlay()
        self.assertEqual([leaf for leaf, _r in overlay.dots], ["S4"])
        self.assertEqual([chip[0] for chip in overlay.chips], ["S4"])
        win.ctx.hover(["S1", "S2", "S3", "S4", "S5"], origin="sample_tree")
        win.ctx.flush_hover()
        self.assertTrue(self._wait_for(lambda: len(tree._active_leaf_set) == 5))
        overlay = tree.hover_overlay()
        # A dot for each tip on screen (the panel may scroll some away); the
        # chips stop at three and the rest is counted.
        self.assertGreaterEqual(len(overlay.dots), 3)
        self.assertEqual(len(overlay.chips), 3)
        # "+N more" counts the rest on screen; the tips scrolled away get an
        # edge marker instead, so all five are counted once.
        more = int(overlay.more[0].split()[0]) if overlay.more else 0
        self.assertEqual(more, len(overlay.dots) - 3)
        off = sum(int(text.split()[0]) for _s, text, _r in overlay.edges)
        self.assertEqual(len(overlay.dots) + off, 5)
        win.ctx.hover([], origin="sample_tree")
        win.ctx.flush_hover()
        self.assertTrue(self._wait_for(lambda: not tree.hover_dimmed()))
        self.assertEqual(tree.hover_overlay().dots, [])
        self.assertFalse(any(item.pen().isCosmetic() for items in tree._edge_items_by_leaf.values() for item in items))

    def test_a_heatmap_shows_its_own_hover_and_shares_it(self) -> None:
        win = self._window()
        panel = self._open(win, "gene_presence_absence")
        heatmap = panel.view()
        heatmap.sampleHovered.emit("S5")
        win.ctx.flush_hover()
        self.assertEqual(heatmap._highlight_sample, "S5")
        self.assertEqual(
            (win.ctx.hover_state.samples, win.ctx.hover_state.origin), (("S5",), "gene_presence_absence")
        )
        heatmap.sampleHovered.emit(None)
        win.ctx.flush_hover()
        self.assertIsNone(heatmap._highlight_sample)
        self.assertEqual(win.ctx.hover_state.samples, ())

    def test_plasmidome_hover_shares_every_sample_of_a_cluster(self) -> None:
        win = self._window()
        panel = self._open(win, "plasmidome", to_side=True)
        self.assertTrue(win._panel_manager.is_visible("sample_tree"))
        panel.view().samplesHovered.emit(["S1", "S4"])
        win.ctx.flush_hover()
        self.assertEqual((win.ctx.hover_state.samples, win.ctx.hover_state.origin), (("S1", "S4"), "plasmidome"))
        self.assertEqual(win._tree_hover_shown, ("S1", "S4"))

    def test_clade_cluster_and_plasmid_clicks_narrow_the_window(self) -> None:
        win = self._window()
        origins: List[str] = []
        win.ctx.subsetChanged.connect(lambda _subset, origin: origins.append(origin))
        tree = self._open(win, "phylogeny").tab.tree
        tree.cladeClicked.emit({"name": "2", "samples": ["S3", "S4", "S5"]})
        self._idle(win)
        self.assertEqual(win.ctx.applied.subset.label, "Clade 2")
        self.assertEqual(win.ctx.result.samples, ("S3", "S4", "S5"))
        self.assertEqual(origins[-1], "phylogeny")
        tree.clusterClicked.emit({"name": "3", "samples": ["S5", "S6"]})
        self._idle(win)
        self.assertEqual(win.ctx.applied.subset.label, "Cluster 3")
        plasmidome = self._open(win, "plasmidome").view()
        plasmidome.sampleSubsetRequested.emit({"name": "AA001", "samples": ["S1", "S4"]})
        self._idle(win)
        self.assertEqual((win.ctx.applied.subset.label, win.ctx.applied.subset.kind), ("Plasmid AA001", "plasmid"))
        self.assertEqual(win.ctx.result.samples, ("S1", "S4"))
        self.assertEqual(origins[-1], "plasmidome")

    def test_open_requests_go_to_the_window(self) -> None:
        win = self._window()
        requests: List[Any] = []
        win.ctx.openRequested.connect(lambda request, origin: requests.append((request, origin)))
        tree = self._open(win, "phylogeny").tab.tree
        with mock.patch.object(win, "_open_sample_card_for_sample") as card:
            tree.leafDoubleClicked.emit("S2")
            card.assert_called_once_with("S2")
        self.assertIsInstance(requests[-1][0], OpenSampleCard)
        self.assertEqual(requests[-1][1], "phylogeny")
        with mock.patch.object(win.app_controller.commands, "run") as run:
            tree.cladeWindowRequested.emit({"name": "1", "samples": ["S1", "S2", "S6"]})
            self.assertTrue(run.called)
        self.assertIsInstance(requests[-1][0], OpenSubsetInNewWindow)
        self.assertEqual(requests[-1][0].subset.samples, ("S1", "S2", "S6"))
        plasmidome = self._open(win, "plasmidome").view()
        with mock.patch.object(win, "_on_pan_plasmid_context_requested") as context:
            plasmidome.contextRequested.emit({"name": "AA001", "rows": [{"sample": "S1", "contig": "3"}]})
            self.assertTrue(context.called)
        self.assertIsInstance(requests[-1][0], OpenPlasmidContext)
        self.assertEqual(requests[-1][1], "plasmidome")
        # The window carries it out with the neighbourhood window command
        # (the requesting panel does not apply in the new window).
        with mock.patch.object(win.app_controller.commands, "run") as run:
            plasmidome.contextRequested.emit({"name": "AA001", "rows": [{"sample": "S1", "contig": "3"}]})
        self.assertEqual(run.call_args.args[0], "window.open_neighbourhood")
        self.assertIs(run.call_args.args[1], win)
        self.assertNotIn("beside", run.call_args.kwargs)
        # The window's sample set did not change.
        self.assertIsNone(win.ctx.applied.subset)

    def test_unavailable_panels_are_disabled_with_the_reason(self) -> None:
        empty = tempfile.mkdtemp(prefix="bactolomew_pan_empty_")
        self.addCleanup(shutil.rmtree, empty, True)
        win = _Explorer.build(empty)
        self.windows.append(win)
        _Explorer.show_and_settle(self.app, win, seconds=0.3)
        manager = win._panel_manager
        missing = {
            "phylogeny": "Not available: this run has no phylogenetic tree.",
            "gene_presence_absence": "Not available: this run has no Panaroo output.",
            "pan_force": "Not available: this run has no chewBBACA (cgMLST) output.",
        }
        for kind, reason in missing.items():
            self.assertEqual(manager.is_available(kind), (False, reason))
            self.assertFalse(manager.open_panel(kind))
            self.assertIsNone(manager.open_new_instance(kind))
        # The synthetic metadata and plasmid annotations are there.
        self.assertTrue(manager.is_available("pan_metadata")[0])
        self.assertTrue(manager.is_available("plasmidome")[0])
        # The Panels menu lists them disabled, with the reason as tooltip.
        manager._fill_panels_menu()
        actions = {a.text(): a for a in manager._panels_menu.actions()}
        for kind, reason in missing.items():
            action = actions[manager.title(kind)]
            self.assertFalse(action.isEnabled(), kind)
            self.assertEqual(action.toolTip(), reason)
        from annotation_explorer_polars.panels.palette import PanelPalette

        palette = PanelPalette(manager.specs(), is_open=manager.is_open)
        try:
            entries = {item.data(0, QtCore.Qt.ItemDataRole.UserRole + 1): item for item in palette._items}
            for kind in missing:
                self.assertTrue(entries[kind].isDisabled(), kind)
        finally:
            palette.deleteLater()


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class HoverHubRetiredTests(_PanWindowCase):
    """PanCoreAccessoryTab links nothing itself; the Phylogeny tree and the
    Gene presence/absence heatmap meet on the window's bus."""

    def test_the_tab_has_no_hover_hub_and_no_view_selector(self) -> None:
        for name in (
            "_hover_update_timer",
            "_schedule_hover_highlight",
            "_flush_hover_highlight",
            "_pending_hover_samples",
            "sampleHovered",
            # neoPan's two view selectors and embedded maps are gone.
            "viewSamplesHovered",
            "view_for_key",
            "view_key",
            "all_geo_views",
            "update_geo_data",
            "_create_selector_panel",
            "_ensure_detail_view_ready",
        ):
            self.assertFalse(hasattr(PanCoreAccessoryTab, name), name)
        import importlib.util

        self.assertIsNone(importlib.util.find_spec("annotation_explorer_polars.panels.neopan"))
        win = self._window()
        tab = self._open(win, "phylogeny").tab
        self.assertFalse(tab.findChildren(QtCore.QTimer) and any(
            t.interval() == 24 for t in tab.findChildren(QtCore.QTimer)
        ))

    def test_a_tree_hover_reaches_the_heatmap_through_the_bus(self) -> None:
        win = self._window()
        tree = self._open(win, "phylogeny").tab.tree
        heatmap = self._open(win, "gene_presence_absence").view()
        self.assertTrue(tree.isVisible() and heatmap.isVisible())
        published: List[Any] = []
        win.ctx.bus.subscribe("probe", "hover_changed", published.append)
        with mock.patch.object(heatmap, "highlight_samples") as heat, mock.patch.object(
            tree, "highlight_samples"
        ) as own_tree:
            tree.subtreeHovered.emit(["S1", "S2"])
            # Nothing is linked inside the panels: the heatmap waits for the bus.
            heat.assert_not_called()
            win.ctx.flush_hover()
            self.assertEqual(win.ctx.hover_state.origin, "phylogeny")
            self.assertTrue(self._wait_for(lambda: any(c.args and c.args[0] == ["S1", "S2"] for c in heat.call_args_list)))
            # The tree it came from is skipped (it marks its own hover).
            own_tree.assert_not_called()
        self.assertTrue(published)

    def test_a_heatmap_hover_reaches_the_tree(self) -> None:
        win = self._window()
        tree = self._open(win, "phylogeny").tab.tree
        heatmap = self._open(win, "gene_presence_absence").view()
        with mock.patch.object(tree, "highlight_sample") as tree_one:
            heatmap.sampleHovered.emit("S6")
            win.ctx.flush_hover()
            self.assertEqual(win.ctx.hover_state.origin, "gene_presence_absence")
            self.assertTrue(self._wait_for(lambda: tree_one.called))
            tree_one.assert_called_with("S6")
        # The heatmap shows its own hover at once.
        self.assertEqual(heatmap._highlight_sample, "S6")


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class SharedPanDataTests(unittest.TestCase):
    """One PanData per run, shared by every pan panel of every window."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.app = _app()
        import annotation_explorer_polars.geo_experimental_views as gev

        cls._saved_webengine_cache = gev._WEBENGINE_RUNTIME_CACHE
        gev._WEBENGINE_RUNTIME_CACHE = (False, "disabled in tests")

    @classmethod
    def tearDownClass(cls) -> None:
        import annotation_explorer_polars.geo_experimental_views as gev

        gev._WEBENGINE_RUNTIME_CACHE = cls._saved_webengine_cache

    def setUp(self) -> None:
        self.output_dir = tempfile.mkdtemp(prefix="bactolomew_pan_shared_")
        write_pan_outputs(self.output_dir)
        settings = QtCore.QSettings()
        settings.remove("explore")
        settings.sync()
        self.windows: List[Any] = []

    def tearDown(self) -> None:
        for win in self.windows:
            win.close()
            win.deleteLater()
        _spin(self.app, 0.05)
        shutil.rmtree(self.output_dir, ignore_errors=True)
        settings = QtCore.QSettings()
        settings.remove("explore")
        settings.sync()

    def _settled(self, win: Any) -> Any:
        self.windows.append(win)
        win.resize(1500, 950)
        win.show()
        done = {"v": win.ctx.first_result_done}
        win.firstRefreshFinished.connect(lambda: done.__setitem__("v", True))
        end = QtCore.QDeadlineTimer(20000)
        while not done["v"] and not end.hasExpired():
            _spin(self.app, 0.02)
        _spin(self.app, 1.0)
        return win

    def test_two_windows_share_the_run_pan_data(self) -> None:
        run = _run(self.output_dir)
        controller = AppController(run)
        calls: List[Any] = []
        real_memo = run.pan.memo

        def counting_memo(key: Any, factory: Any) -> Any:
            def counted() -> Any:
                calls.append(key)
                return factory()

            return real_memo(key, counted)

        tabs = []
        with mock.patch.object(run.pan, "memo", side_effect=counting_memo):
            # Classic shows Phylogeny and Gene presence/absence in both windows.
            first = self._settled(controller.open_first_window())
            second = self._settled(controller.open_explore(source=first))
            for win in (first, second):
                manager = win._panel_manager
                for kind in ("gene_presence_absence", "plasmidome", "phylogeny"):
                    self.assertTrue(manager.open_panel(kind))
                    end = QtCore.QDeadlineTimer(10000)
                    while manager.panel(kind) is None and not end.hasExpired():
                        _spin(self.app, 0.05)
                    tabs.append(manager.panel(kind).tab)
        # The leaf order of the Gene presence/absence rows and the plasmid gene rows
        # are each computed once for the run, whatever the number of windows.
        self.assertEqual(sorted(key[0] for key in calls), ["leaf_order", "plasmid_gene_rows"])
        for tab in tabs:
            self.assertIs(tab._pan_data, run.pan)
            self.assertIs(tab._plasmid_cluster_index, run.plasmid_cluster_index)
            self.assertIs(tab._ann_df, run.ann_df)  # shared, not copied
        # The windows stay independent: a subset in one leaves the other alone.
        first._set_sample_subset_filter(label="Two", samples=["S1", "S2"])
        end = QtCore.QDeadlineTimer(10000)
        while first.ctx.is_busy() and not end.hasExpired():
            _spin(self.app, 0.05)
        _spin(self.app, 0.5)
        for win in (first, second):
            # Show Gene presence/absence (it catches up if hidden).
            win._panel_manager.open_panel("gene_presence_absence")
            _spin(self.app, 0.6)
        self.assertEqual(first._panel_manager.panel("gene_presence_absence").tab._sample_filter, {"S1", "S2"})
        self.assertEqual(first._panel_manager.panel("plasmidome").tab._sample_filter, {"S1", "S2"})
        self.assertEqual(len(second.ctx.result.samples), 6)
        self.assertEqual(
            second._panel_manager.panel("gene_presence_absence").tab._sample_filter, set(second.ctx.result.samples)
        )


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class PangenomePresetTests(_PanWindowCase):
    def test_pangenome_puts_the_tree_beside_its_heatmap_and_the_map_beside_plasmidome(self) -> None:
        win = self._window()
        manager = win._panel_manager
        computes = []
        win.ctx.bus.subscribe("probe", "busy_changed", computes.append)
        applied, result = win.ctx.applied, win.ctx.result
        self.assertTrue(manager.apply_preset(layouts.PRESET_PANGENOME))
        _spin(self.app, 0.5)
        self.assertEqual(
            _sorted_groups(manager.host.groups()),
            [["gene_presence_absence"], ["geo"], ["phylogeny"], ["plasmidome"]],
        )
        self.assertEqual(manager.host.floating_count(), 0)
        frames = {kind: manager.frame(kind) for kind in ("phylogeny", "gene_presence_absence", "plasmidome", "geo")}
        pos = {kind: f.mapTo(win, QtCore.QPoint(0, 0)) for kind, f in frames.items()}
        self.assertLess(pos["phylogeny"].x(), pos["gene_presence_absence"].x())  # the tree on the left
        self.assertEqual(pos["phylogeny"].y(), pos["gene_presence_absence"].y())  # one row, same height
        self.assertEqual(frames["phylogeny"].height(), frames["gene_presence_absence"].height())
        self.assertLess(pos["geo"].x(), pos["plasmidome"].x())  # the map left of Plasmidome
        self.assertGreater(pos["geo"].y(), pos["phylogeny"].y())  # below them
        self.assertEqual(pos["geo"].y(), pos["plasmidome"].y())
        self.assertGreater(frames["phylogeny"].height(), frames["geo"].height())
        # Panels shown by the preset catch up; nothing was computed again.
        for kind in ("phylogeny", "gene_presence_absence", "plasmidome"):
            frame = frames[kind]
            self.assertTrue(self._wait_for(lambda f=frame: not f.is_dirty() and f.panel() is not None), kind)
        self.assertIs(win.ctx.applied, applied)
        self.assertIs(win.ctx.result, result)
        self.assertEqual(computes, [])
        # Their rows are linked (tests/test_row_link.py has the details).
        heatmap = manager.panel("gene_presence_absence")
        self.assertTrue(self._wait_for(heatmap.is_linked))

    def _fake(self, available: tuple = ("phylogeny", "gene_presence_absence", "plasmidome")) -> "PanelManager":
        specs = [spec for spec in _fake_specs() if spec.kind not in ("phylogeny", "gene_presence_absence")]
        for kind in ("phylogeny", "gene_presence_absence", "plasmidome"):
            ok = kind in available
            specs.append(
                PanelSpec(
                    kind,
                    kind.replace("_", " ").title(),
                    "Pangenome",
                    placeholder=kind,
                    available=(lambda ok=ok: (ok, "" if ok else "Not available: test reason")),
                    single_instance=False,
                )
            )
        window = QtWidgets.QMainWindow()
        manager = PanelManager(window, specs)
        splitter = QtWidgets.QSplitter()
        splitter.addWidget(QtWidgets.QLabel("rail"))
        splitter.addWidget(manager.host)
        window.setCentralWidget(splitter)
        window.resize(1300, 850)
        window.show()
        manager.build_default()
        _spin(self.app, 0.1)
        self.addCleanup(window.deleteLater)
        self.addCleanup(window.close)
        return manager

    def test_pangenome_leaves_out_what_is_missing(self) -> None:
        cases = {
            ("phylogeny", "gene_presence_absence", "plasmidome"): [
                ["gene_presence_absence"], ["geo"], ["phylogeny"], ["plasmidome"]
            ],
            ("phylogeny", "gene_presence_absence"): [["gene_presence_absence"], ["geo"], ["phylogeny"]],
            ("gene_presence_absence", "plasmidome"): [["gene_presence_absence"], ["geo"], ["plasmidome"]],
            ("plasmidome",): [["geo"], ["plasmidome"]],
        }
        for available, groups in cases.items():
            with self.subTest(available=available):
                manager = self._fake(available)
                self.assertTrue(manager.apply_preset(layouts.PRESET_PANGENOME))
                _spin(self.app, 0.1)
                self.assertEqual(_sorted_groups(manager.host.groups()), groups)
        # Nothing of it: refused, with the panels it needs named.
        manager = self._fake(())
        available, reason = manager.preset_available(layouts.PRESET_PANGENOME)
        self.assertFalse(available)
        self.assertIn("Phylogeny", reason)
        self.assertIn("Gene Presence Absence", reason)

    def test_a_layout_saved_before_the_pan_panels_restores_with_them_closed(self) -> None:
        old = self._fake_without_pan_kinds()
        state, panels = old.host.save_layout(), old.panel_list()
        groups = _sorted_groups(old.host.groups())
        new = self._fake()
        new.apply_preset(layouts.PRESET_PANGENOME)
        self.assertTrue(new.apply_saved_layout(state, panels))
        _spin(self.app, 0.2)
        self.assertEqual(_sorted_groups(new.host.groups()), groups)
        for kind in ("phylogeny", "gene_presence_absence", "plasmidome"):
            self.assertFalse(new.is_open(kind), kind)
            self.assertTrue(new.open_panel(kind), kind)  # and they open fine afterwards

    def _fake_without_pan_kinds(self) -> "PanelManager":
        window = QtWidgets.QMainWindow()
        specs = [spec for spec in _fake_specs() if spec.kind not in ("phylogeny", "gene_presence_absence")]
        manager = PanelManager(window, specs)
        splitter = QtWidgets.QSplitter()
        splitter.addWidget(QtWidgets.QLabel("rail"))
        splitter.addWidget(manager.host)
        window.setCentralWidget(splitter)
        window.resize(1300, 850)
        window.show()
        manager.build_default()
        _spin(self.app, 0.1)
        self.addCleanup(window.deleteLater)
        self.addCleanup(window.close)
        return manager


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class PanPanelLinkingTests(_PanWindowCase):
    """Multi-instance hover, replay on show, hover in, dispose and layouts of
    the pan panels (review round of Phase 9)."""

    def _duplicate(self, win: Any, kind: str, where: str = "right") -> Any:
        manager = win._panel_manager
        uid = manager.duplicate_panel(kind, where)
        self.assertIsNotNone(uid, kind)
        self.assertTrue(self._wait_for(lambda: manager.panel(uid) is not None and manager.is_visible(uid)), uid)
        _spin(self.app, 0.4)
        return manager.panel(uid)

    def test_two_phylogeny_instances_link_their_hover_both_ways(self) -> None:
        win = self._window()
        first = self._open(win, "phylogeny")
        second = self._duplicate(win, "phylogeny")
        self.assertNotEqual(second.uid, first.uid)
        self.assertTrue(win._panel_manager.is_visible(first.uid))
        for source, target in ((second, first), (first, second)):
            with mock.patch.object(target.tab, "highlight_sample") as other, mock.patch.object(
                source.tab, "highlight_sample"
            ) as own:
                source.tab.tree.leafHovered.emit("S4")
                win.ctx.flush_hover()
                self.assertEqual(win.ctx.hover_state.origin, source.uid)
                self.assertTrue(self._wait_for(lambda: other.called), source.uid)
                other.assert_called_with("S4")
                own.assert_not_called()  # never echoed
                source.tab.tree.leafHovered.emit(None)
                win.ctx.flush_hover()
                _spin(self.app, 0.1)

    def test_a_hidden_standalone_pan_panel_catches_up_on_hover_when_shown(self) -> None:
        win = self._window()
        manager = win._panel_manager
        panel = self._open(win, "gene_presence_absence")
        manager.host.move_view("gene_presence_absence", Placement("tab", relative_to="geo"))
        manager.host.set_current_tab("geo")
        _spin(self.app, 0.3)
        self.assertFalse(manager.is_visible("gene_presence_absence"))
        with mock.patch.object(panel.tab, "highlight_sample") as one:
            win.ctx.hover(["S5"], origin="sample_tree")
            win.ctx.flush_hover()
            _spin(self.app, 0.1)
            one.assert_not_called()
            manager.host.set_current_tab("gene_presence_absence")
            _spin(self.app, 0.3)
            one.assert_any_call("S5")

    def test_plasmidome_tints_the_rows_of_hovered_samples_and_pan_metadata_ignores_hover(self) -> None:
        win = self._window()
        plasmidome = self._open(win, "plasmidome", to_side=True)
        metadata = self._open(win, "pan_metadata", to_side=True)
        view = plasmidome.view()
        rows = {row["name"]: i for i, row in view._row_by_table_row.items()}
        self.assertEqual(view.linked_rows(), [])
        win.ctx.hover(["S4"], origin="sample_tree")
        win.ctx.flush_hover()
        self.assertTrue(self._wait_for(lambda: bool(view.linked_rows())))
        # Exactly the rows of the plasmids S4 carries are tinted.
        carrying = sorted(i for i, row in view._row_by_table_row.items() if "S4" in row["samples"])
        self.assertTrue(carrying)
        self.assertLess(len(carrying), len(rows))
        self.assertEqual(sorted(view.linked_rows()), carrying)
        item = view.table.item(view.linked_rows()[0], 1)
        self.assertNotEqual(item.background().style(), QtCore.Qt.BrushStyle.NoBrush)
        win.ctx.hover([], origin="sample_tree")
        win.ctx.flush_hover()
        self.assertTrue(self._wait_for(lambda: not view.linked_rows()))
        self.assertEqual(view.table.item(0, 1).background().style(), QtCore.Qt.BrushStyle.NoBrush)
        # Metadata summary (aggregate donuts) has no per-sample mark: hover in is
        # taken without effect (a documented deviation from plan 4.6).
        self.assertFalse(hasattr(metadata.view(), "highlight_samples"))
        metadata.show_hover(SimpleNamespace(samples=("S1",), origin="sample_tree"))

    def test_plasmidome_hover_is_not_echoed_into_its_own_rows(self) -> None:
        win = self._window()
        panel = self._open(win, "plasmidome", to_side=True)
        view = panel.view()
        view._on_cell_entered(0, 0)
        win.ctx.flush_hover()
        self.assertEqual(win.ctx.hover_state.origin, panel.uid)
        _spin(self.app, 0.2)
        self.assertEqual(view.linked_rows(), [])

    def test_closing_a_duplicate_pan_force_or_phylogeny_keeps_the_first(self) -> None:
        from annotation_explorer_polars.panaroo_force import CgmlstForceGraphView

        win = self._window()
        manager = win._panel_manager
        first = self._open(win, "pan_force")
        dup = self._duplicate(win, "pan_force")
        self.assertIsInstance(dup.view(), CgmlstForceGraphView)
        manager.close_panel(dup.uid)
        _spin(self.app, 0.1)
        self.assertIsNone(manager.panel(dup.uid))
        self.assertIsInstance(first.view(), CgmlstForceGraphView)
        tree = self._open(win, "phylogeny")
        second = self._duplicate(win, "phylogeny")
        manager.close_panel(second.uid)
        _spin(self.app, 0.1)
        self.assertIsNone(manager.panel(second.uid))
        self.assertIsNotNone(tree.tab)

    def test_pan_panels_and_their_options_persist_with_the_layout(self) -> None:
        win = self._window()
        manager = win._panel_manager
        self.assertTrue(manager.apply_preset(layouts.PRESET_PANGENOME))
        _spin(self.app, 0.5)
        tree = self._open(win, "phylogeny")
        dup = self._duplicate(win, "phylogeny", "tab")
        combo = dup.tab._tree_lengths_combo
        combo.setCurrentIndex(combo.findData(LENGTHS_TOPOLOGY))
        _spin(self.app, 0.1)
        self.assertEqual(dup.save_state(), {"layout": "rectangular", "branch_scale": "topology"})
        self.assertEqual(tree.save_state()["branch_scale"], "distance")
        uid = dup.uid
        groups = _sorted_groups(manager.host.groups())
        win.close()
        _spin(self.app, 0.1)
        again = self._window()
        restored_manager = again._panel_manager
        self.assertEqual(_sorted_groups(restored_manager.host.groups()), groups)
        for kind in ("gene_presence_absence", "plasmidome", "phylogeny", uid):
            self.assertTrue(restored_manager.is_open(kind), kind)
        restored_manager.host.set_current_tab(uid)
        self.assertTrue(self._wait_for(lambda: restored_manager.panel(uid) is not None, timeout=5.0))
        _spin(self.app, 0.3)
        restored = restored_manager.panel(uid)
        self.assertIsInstance(restored, PanViewPanel)
        self.assertEqual(restored.tab.tree_display_options()["branch_scale"], "topology")
        self.assertFalse(restored.tab._preserve_branch_lengths)
        restored_manager.host.set_current_tab("phylogeny")
        self.assertTrue(self._wait_for(lambda: restored_manager.panel("phylogeny") is not None, timeout=5.0))
        self.assertEqual(restored_manager.panel("phylogeny").tab.tree_display_options()["branch_scale"], "distance")
        self.assertEqual(restored_manager.host.floating_count(), 0)

    def test_the_window_searches_the_run_folder_for_the_tree_once(self) -> None:
        with mock.patch.object(pan_data_mod, "tree_sources", wraps=pan_data_mod.tree_sources) as spy:
            win = self._window()
            for kind in ("phylogeny", "gene_presence_absence"):
                self._open(win, kind)
            self._duplicate(win, "phylogeny")
            self.assertIs(win.run.pan.tree_index(), win.run.pan.tree_index())
        self.assertEqual(spy.call_count, 1)


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class HeavyPanelHoverBudgetTests(_PanWindowCase):
    """Plan 8.12: hover reaches all visible highlights within 30 ms with four
    heavy panels visible. These are the real Map, Phylogeny and Gene
    presence/absence of Classic's bottom row and a cgMLST MST (compact) beside them
    (WebEngine is off in tests, so the Map is the offline SVG map, whose hover
    used to rebuild the whole map)."""

    def _four_heavy_panels(self, win: Any) -> List[str]:
        manager = win._panel_manager
        geo = manager.focus_panel("geo")
        self.assertEqual(geo, "geo")
        self.assertTrue(manager.open_panel("pan_force", target="geo", to_side=True))
        _spin(self.app, 0.2)
        uids = ["geo", "phylogeny", "gene_presence_absence", "pan_force"]

        def built() -> bool:
            return all(
                manager.is_visible(u) and manager.panel(u) is not None and not manager.frame(u).is_dirty()
                for u in uids
            )

        self.assertTrue(self._wait_for(built, timeout=15.0), {u: manager.is_visible(u) for u in uids})
        self._idle(win)
        return uids

    def test_every_heavy_panel_takes_hover_on_the_heavy_timer(self) -> None:
        win = self._window()
        manager = win._panel_manager
        uids = self._four_heavy_panels(win)
        self.assertEqual(len(uids), 4)
        for uid in uids:
            frame = manager.frame(uid)
            self.assertTrue(manager.registry.spec(frame.kind).heavy, uid)
            hover = [s for s in frame._subs if s.event == "hover_changed"]
            self.assertEqual(len(hover), 1, uid)
            self.assertEqual(hover[0].delivery, Delivery.COALESCED, uid)
            self.assertEqual(hover[0].coalesce_ms, bus_events.HEAVY_HOVER_COALESCE_MS, uid)

    def test_a_hover_reaches_four_real_heavy_panels_within_30_ms(self) -> None:
        win = self._window()
        manager = win._panel_manager
        uids = self._four_heavy_panels(win)
        self.assertEqual(len(uids), 4)
        done: dict = {}
        for uid in uids:
            panel = manager.panel(uid)
            real = panel.show_hover

            def record(hover: Any, _uid: str = uid, _real: Any = real) -> None:
                _real(hover)  # the real highlight; the time is taken after it
                done[_uid] = time.perf_counter()

            panel.show_hover = record
        samples = ["S1", "S4", "S2", "S6", "S3", "S5", "S1"]
        from annotation_explorer_polars.geo_svg_map import SvgGeoView

        geo_view = manager.panel("geo").view._active
        self.assertIsInstance(geo_view, SvgGeoView)

        def median_ms() -> float:
            latencies = []
            for sample in samples:
                win.ctx.hover([], origin="sample_tree")
                _spin(self.app, 0.06)
                done.clear()
                start = time.perf_counter()
                win.ctx.hover([sample], origin="sample_tree", kind="sample")
                loop = QtCore.QEventLoop()
                QtCore.QTimer.singleShot(500, loop.quit)
                poll = QtCore.QTimer()
                poll.setInterval(1)
                poll.timeout.connect(lambda: loop.quit() if len(done) == len(uids) else None)
                poll.start()
                loop.exec()
                poll.stop()
                self.assertEqual(sorted(done), sorted(uids))
                latencies.append((max(done.values()) - start) * 1000)
            return sorted(latencies)[len(latencies) // 2]

        ok, seen = within_budget(median_ms, 30)
        self.assertTrue(ok, f"median ms from ctx.hover to the last heavy panel: {seen}")
        # The highlight really drew: S1 has coordinates in Romania.
        win.ctx.hover(["S1"], origin="sample_tree", kind="sample")
        self.assertTrue(self._wait_for(lambda: bool(geo_view.highlighted_countries()), timeout=1.0))


if __name__ == "__main__":
    unittest.main()
