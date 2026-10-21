"""Gene presence/absence rows beside the Phylogeny (feedback 2026-10-07
item (d), 2026-10-09).

The heatmap always fits every row top to bottom: no vertical scroll or zoom,
whatever the tree does (``pan_view.PanHeatmapWidget``). In one explore window
its rows are the tree's leaves in the tree's order (``panels/row_link.py``);
only the order is linked, and the rows fill the heatmap from its top (no band
left blank to line them up with the tree's leaves). Rows too thin to label
name their sample in the tooltip.

Offscreen on the synthetic run with the pan outputs of ``tests/_pan_fixture.py``.
"""
from __future__ import annotations

import os
import unittest
from pathlib import Path
from typing import Any, List, Tuple

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

from PyQt6 import QtCore, QtGui, QtWidgets

try:
    from tests._pan_fixture import CLADES, SAMPLES
    from tests.test_dock_host import _spin
    from tests.test_pan_panels import _PanWindowCase
except ImportError:
    from _pan_fixture import CLADES, SAMPLES
    from test_dock_host import _spin
    from test_pan_panels import _PanWindowCase

from annotation_explorer_polars.docking import docking_available

if docking_available():
    from annotation_explorer_polars.panaroo_core_accessory import BRANCH_COLOR_CLUSTER
    from annotation_explorer_polars.pan_view import MIN_LABELED_ROW_PX
    from annotation_explorer_polars.tree_view import ROOT_MIDPOINT


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class _RowLinkCase(_PanWindowCase):
    """A window on the synthetic run with pan outputs, and the row helpers."""

    def _pair(self, win: Any) -> Tuple[Any, Any]:
        """The Phylogeny and Gene presence/absence of Classic, linked."""
        manager = win._panel_manager
        tree_panel = self._open(win, "phylogeny")
        heat_panel = self._open(win, "gene_presence_absence")
        self.assertTrue(manager.is_visible("phylogeny") and manager.is_visible("gene_presence_absence"))
        self.assertTrue(self._wait_for(heat_panel.is_linked), "not linked")
        _spin(self.app, 0.2)
        return tree_panel, heat_panel

    def _assert_linked_rows_from_top(self, tree_panel: Any, heat_panel: Any) -> None:
        """The tree's order, and the first row at the top of the plot: the
        rows fill the heatmap, nothing above them is left blank."""
        heatmap = heat_panel.heatmap()
        self.assertEqual(heatmap.display_samples(), tree_panel.tree().leaf_order())
        plot = heatmap._plot
        top = plot.viewportTransform().map(plot.getPlotItem().vb.mapViewToScene(QtCore.QPointF(0.0, 0.0))).y()
        self.assertLessEqual(top, 2.0)

    def _assert_all_rows_fitted(self, heat_panel: Any) -> None:
        """Every row in view: the vertical range is exactly all rows."""
        heatmap = heat_panel.heatmap()
        (_x0, _x1), (y0, y1) = heatmap._plot.getPlotItem().vb.viewRange()
        self.assertAlmostEqual(y0, 0.0, places=4)
        self.assertAlmostEqual(y1, float(len(heatmap.display_samples())), places=4)

    @staticmethod
    def _wheel(heat_panel: Any, angle: int, ctrl: bool = False) -> None:
        plot = heat_panel.heatmap()._plot
        viewport = plot.viewport()
        pos = QtCore.QPointF(viewport.width() / 2.0, viewport.height() / 2.0)
        event = QtGui.QWheelEvent(
            pos,
            QtCore.QPointF(viewport.mapToGlobal(pos.toPoint())),
            QtCore.QPoint(0, 0),
            QtCore.QPoint(0, angle),
            QtCore.Qt.MouseButton.NoButton,
            QtCore.Qt.KeyboardModifier.ControlModifier if ctrl else QtCore.Qt.KeyboardModifier.NoModifier,
            QtCore.Qt.ScrollPhase.NoScrollPhase,
            False,
        )
        QtWidgets.QApplication.sendEvent(viewport, event)

    def _settle(self, seconds: float = 0.2) -> None:
        _spin(self.app, seconds)


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class TreeRowLinkTests(_RowLinkCase):
    def test_phylogeny_draws_the_grouped_tree_with_cluster_branches_and_a_legend(self) -> None:
        win = self._window()
        panel = self._open(win, "phylogeny")
        tab, tree = panel.tab, panel.tree()
        # An old run's unrooted tree file is drawn at its midpoint.
        self.assertEqual(tab.tree_root_mode(), ROOT_MIDPOINT)
        # Genomes outside the sample set are greyed, as neoPan drew them.
        self.assertFalse(tab.prune_tree_to_sample_filter())
        self.assertEqual(tab.branch_color_mode(), BRANCH_COLOR_CLUSTER)
        # Branches colored by cluster, labels neutral, clade labels replaced
        # by the legend, cluster brackets kept.
        self.assertEqual(tree._leaf_display_group_map, {s: CLADES[s][1] for s in SAMPLES})
        self.assertFalse(tree._display_groups_color_labels)
        self.assertFalse(tree._show_group_labels)
        self.assertTrue(tree._show_cluster_annotations)
        # The clades are grouped (contiguous) in the drawn order.
        lineages = [CLADES[s][0] for s in tree.leaf_order()]
        self.assertEqual(lineages, sorted(lineages, key=lineages.index))
        legend = tab._legend
        self.assertTrue(legend.isVisible())
        clades, clusters = tab.tree_legend_entries()
        self.assertEqual([(name, count) for name, count, _color in clades], [("1", 3), ("2", 3)])
        self.assertEqual([name for name, _color in clusters], ["1", "2", "3"])
        # One color source: the legend chips are the tree's group colors.
        for name, _count, color in clades:
            self.assertEqual(color, tree._group_color_map[name].name())
        for name, color in clusters:
            self.assertEqual(color, tree._display_group_color_map[name].name())
        self.assertIn("Clade", legend.text())
        self.assertIn("Cluster", legend.text())
        # The legend counts the genomes of the sample set: it follows a
        # subset (the others stay drawn, greyed).
        win._set_sample_subset_filter(label="Two", samples=["S1", "S3", "S4"])
        self._idle(win)
        self.assertEqual(sorted(tree.leaf_order()), list(SAMPLES))
        clades, clusters = tab.tree_legend_entries()
        self.assertEqual([(name, count) for name, count, _color in clades], [("1", 1), ("2", 2)])
        self.assertEqual([name for name, _color in clusters], ["1", "2"])

    def test_the_heatmap_keeps_every_row_whatever_the_tree_does(self) -> None:
        win = self._window()
        tree_panel, heat_panel = self._pair(win)
        self.assertEqual(win._row_link.leader_uid(), "phylogeny")
        self._assert_linked_rows_from_top(tree_panel, heat_panel)
        self._assert_all_rows_fitted(heat_panel)
        tree = tree_panel.tree()
        order = tree.leaf_order()
        # Zooming and scrolling the tree never moves the heatmap's rows.
        tree.zoom_by(3.0)
        self._settle()
        bar = tree.verticalScrollBar()
        self.assertGreater(bar.maximum(), 0)
        bar.setValue(bar.maximum())
        self._settle()
        self.assertEqual(heat_panel.heatmap().display_samples(), order)
        self._assert_all_rows_fitted(heat_panel)
        # A wheel over the heatmap zooms its genes only and leaves the tree alone.
        vb = heat_panel.heatmap()._plot.getPlotItem().vb
        genes_before = vb.viewRange()[0]
        scale, value = tree.transform().m22(), bar.value()
        self._wheel(heat_panel, 240)
        self._wheel(heat_panel, 240, ctrl=True)
        self._settle()
        self.assertNotEqual(vb.viewRange()[0], genes_before)
        self._assert_all_rows_fitted(heat_panel)
        self.assertEqual((tree.transform().m22(), bar.value()), (scale, value))

    def test_a_filter_dims_rows_and_keeps_them_all_in_view(self) -> None:
        win = self._window()
        tree_panel, heat_panel = self._pair(win)
        win._set_sample_subset_filter(label="Four", samples=["S1", "S2", "S5", "S6"])
        self._idle(win)
        self._settle()
        # The tree greys the others; the heatmap dims their rows.
        self.assertEqual(sorted(tree_panel.tree().leaf_order()), list(SAMPLES))
        self.assertEqual(sorted(heat_panel.heatmap().dimmed_rows()), ["S3", "S4"])
        self._assert_all_rows_fitted(heat_panel)
        self._assert_linked_rows_from_top(tree_panel, heat_panel)
        # Without the tree the heatmap shows its own rows, all fitted.
        win._panel_manager.close_panel("phylogeny")
        self.assertTrue(self._wait_for(lambda: not heat_panel.is_linked()))
        self.assertEqual(sorted(heat_panel.heatmap().display_samples()), ["S1", "S2", "S5", "S6"])
        self._assert_all_rows_fitted(heat_panel)

@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class LinkedRowsDataTests(_PanWindowCase):
    """Linked rows are exactly the tree's leaves, whatever the matrix has."""

    def _drop_from_matrix(self, sample: str) -> None:
        path = Path(self.output_dir) / "panaroo" / "out" / "gene_presence_absence.csv"
        rows: List[List[str]] = [line.split(",") for line in path.read_text(encoding="utf-8").splitlines()]
        index = rows[0].index(sample)
        for row in rows:
            del row[index]
        path.write_text("\n".join(",".join(row) for row in rows) + "\n", encoding="utf-8")

    def _drop_from_tree(self, sample: str) -> None:
        path = Path(self.output_dir) / "trees" / "core_gene_tree.nwk"
        text = path.read_text(encoding="utf-8")
        self.assertIn("((S3:0.10,S4:0.12):0.04,S5:0.30):0.03", text)
        self.assertEqual(sample, "S5")
        path.write_text(text.replace("((S3:0.10,S4:0.12):0.04,S5:0.30):0.03", "(S3:0.10,S4:0.12):0.07"), encoding="utf-8")

    def test_a_leaf_missing_from_the_matrix_is_a_blank_row_and_samples_off_the_tree_are_left_out(self) -> None:
        self._drop_from_matrix("S6")
        self._drop_from_tree("S5")
        win = self._window()
        tree_panel = self._open(win, "phylogeny")
        heat_panel = self._open(win, "gene_presence_absence")
        self.assertTrue(self._wait_for(heat_panel.is_linked))
        _spin(self.app, 0.2)
        heatmap = heat_panel.heatmap()
        order = tree_panel.tree().leaf_order()
        self.assertIn("S6", order)
        self.assertEqual(heatmap.display_samples(), order)
        image = heatmap._image_item.image
        blank = order.index("S6")
        self.assertTrue((image[blank] == 2).all())  # the blank cell value
        self.assertTrue(all((image[i] != 2).all() for i in range(len(order)) if i != blank))
        # S5 is in the matrix and the sample set but not a leaf: no row, and
        # the axis says so.
        self.assertNotIn("S5", order)
        self.assertNotIn("S5", heatmap.display_samples())
        self.assertEqual(heatmap.linked_omitted(), 1)
        self.assertIn("1 not in tree", heatmap._plot.getPlotItem().getAxis("left").labelText)
        # Unlinked, the heatmap lists the matrix's own samples again.
        win._panel_manager.close_panel("phylogeny")
        self.assertTrue(self._wait_for(lambda: not heat_panel.is_linked()))
        self.assertIn("S5", heatmap.display_samples())
        self.assertNotIn("S6", heatmap.display_samples())
        self.assertEqual(heatmap._plot.getPlotItem().getAxis("left").labelText, "Samples")


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class HeatmapPanelTests(_RowLinkCase):
    def test_dragging_the_genes_stays_within_the_matrix(self) -> None:
        win = self._window()
        _tree_panel, heat_panel = self._pair(win)
        heatmap = heat_panel.heatmap()
        vb = heatmap._plot.getPlotItem().vb
        cols = heatmap._image_item.image.shape[1]
        vb.setXRange(-500.0, cols + 500.0, padding=0.0)
        self._settle()
        x0, x1 = vb.viewRange()[0]
        self.assertGreaterEqual(x0, -1e-6)
        self.assertLessEqual(x1, cols + 1e-6)
        vb.translateBy(x=-1000.0)
        self._settle()
        self.assertGreaterEqual(vb.viewRange()[0][0], -1e-6)

    def test_the_palette_opens_the_pair_side_by_side_so_they_link(self) -> None:
        for first, second in (("phylogeny", "gene_presence_absence"), ("gene_presence_absence", "phylogeny")):
            with self.subTest(first=first):
                win = self._window()
                manager = win._panel_manager
                # Genome focus has neither: like a neighbourhood window, the
                # pair was never placed, so the palette's plain open applies.
                self.assertTrue(manager.apply_preset("genome_focus"))
                self._settle(0.3)
                self.assertFalse(manager.host.is_placed(first) or manager.host.is_placed(second))
                self.assertTrue(manager.open_panel(first))
                self._settle(0.3)
                # The second one goes beside the first, not into its tabs.
                self.assertTrue(manager.open_panel(second))
                self.assertTrue(self._wait_for(manager.panel("gene_presence_absence").is_linked), "not linked")
                self.assertNotIn(second, manager.host.group_uids(first))
                self.assertTrue(manager.is_visible(first) and manager.is_visible(second))


# Leaves added to the synthetic tree so it has more rows than a panel's
# height at the tree's own row pitch (the reference run has 153).
_EXTRA = tuple(f"X{i:03d}" for i in range(1, 115))


@unittest.skipUnless(docking_available(), "PyQt6Ads is not installed")
class ManyRowsTests(_RowLinkCase):
    """A tree with more rows than fit at a readable pitch (the reference run
    has 153): the heatmap still shows every row at any window size, and drops
    its row labels for the tooltip when the rows are too thin."""

    def setUp(self) -> None:
        super().setUp()
        out = Path(self.output_dir)
        # A caterpillar of extra leaves hung beside the fixture's tree.
        extra = "S6:0.40"
        for name in _EXTRA:
            extra = f"({extra},{name}:0.05):0.01"
        tree_path = out / "trees" / "core_gene_tree.nwk"
        text = tree_path.read_text(encoding="utf-8")
        self.assertIn(",S6:0.40);", text)
        tree_path.write_text(text.replace(",S6:0.40);", f",{extra});"), encoding="utf-8")
        matrix = out / "panaroo" / "out" / "gene_presence_absence.csv"
        lines = matrix.read_text(encoding="utf-8").splitlines()
        lines[0] += "".join(f",{name}" for name in _EXTRA)
        for i in range(1, len(lines)):
            lines[i] += "".join(f",g{i}_{name}" if (i + n) % 3 == 0 else "," for n, name in enumerate(_EXTRA))
        matrix.write_text("\n".join(lines) + "\n", encoding="utf-8")
        self.total = len(SAMPLES) + len(_EXTRA)

    def test_every_row_stays_in_view_at_any_size_and_thin_rows_lose_their_labels(self) -> None:
        win = self._window()
        tree_panel, heat_panel = self._pair(win)
        heatmap = heat_panel.heatmap()
        self.assertEqual(len(heatmap.display_samples()), self.total)
        for height in (617, 1013):
            with self.subTest(height=height):
                win.resize(1301, height)
                self._settle(0.4)
                self._assert_all_rows_fitted(heat_panel)
                tree_panel.tree().zoom_by(4.0)
                self._settle()
                self._assert_all_rows_fitted(heat_panel)
                tree_panel.tree().reset_view()
                self._settle()
        # 120 rows in a few hundred pixels: too thin to label.
        pitch = heatmap._plot.getPlotItem().vb.height() / self.total
        self.assertLess(pitch, MIN_LABELED_ROW_PX)
        self.assertFalse(heatmap.row_labels_shown())
        # The hovered row names its sample in the tooltip.
        vb = heatmap._plot.getPlotItem().vb
        row = heatmap.display_samples().index("S3")
        heatmap._update_hover_from_scene(vb.mapViewToScene(QtCore.QPointF(1.0, row + 0.5)))
        self.assertEqual(heatmap._plot.viewport().toolTip(), "S3")


if __name__ == "__main__":
    unittest.main()
