from __future__ import annotations

import os
import shutil
import tempfile
import unittest
from unittest import mock
from pathlib import Path
from typing import Any, FrozenSet, List, Set

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

from PyQt6 import QtCore, QtGui, QtTest, QtWidgets

from annotation_explorer_polars.panaroo_common import (
    _FASTGEAR_LINEAGE_PALETTE,
    _METADATA_STRIP_PALETTE,
)
from annotation_explorer_polars.panaroo_core_accessory import LENGTHS_DISTANCE, LENGTHS_TOPOLOGY, PanCoreAccessoryTab
from annotation_explorer_polars import phylo
from annotation_explorer_polars.phylo import (
    GroupOptimizer,
    clone_tree,
    group_root,
    groups_contiguous,
    merge_unary,
)
from annotation_explorer_polars.panaroo_style import _fastgear_lineage_colors
from annotation_explorer_polars.theme import COLORS
from annotation_explorer_polars.row_range_note import row_range_text, visible_rows
from annotation_explorer_polars.tree_view import (
    _ActionTextItem,
    _CLADE_FONT_SIZE,
    _CLADE_RAIL_WIDTH,
    _CLUSTER_FONT_SIZE,
    _LEAF_FONT_SIZE,
    _SCENE_DATA_KIND,
    _SCENE_DATA_NAME,
    _SCENE_DATA_RUN_COUNT,
    _SCENE_DATA_RUN_INDEX,
    _SCENE_DATA_VALUE,
    _TREE_BRANCH_COLOR,
    _TREE_CLUSTER_COLOR,
    _TREE_INACTIVE_COLOR,
    _TREE_LABEL_COLOR,
    _TREE_REFERENCE_COLOR,
    _HOVER_DIM_OPACITY,
    _HOVER_PATH_WIDTH_PX,
    _TreeEdgeItem,
    RaxmlTreeView,
    _NewickParser,
    _TreeNode,
    _midpoint_root,
    _prune_to,
)


TREE_TEXT = (
    "(((A:0.1,B:0.2):0.3,(C:0.1,D:0.1):0.2):0.1,"
    "((E:10,F:0.1):0.1,G:0.1):0.1);"
)
GROUPS = {
    "A": "1",
    "B": "1",
    "C": "1",
    "D": "1",
    "E": "2",
    "F": "2",
    "G": "2",
}
GROUP_COLORS = {"1": "#3f6ea5", "2": "#35826f"}
CLUSTERS = {
    "A": "3",
    "B": "3",
    "C": "1",
    "D": "1",
    "E": "3",
    "F": "3",
    "G": "2",
}


class TreeViewTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def setUp(self) -> None:
        self._views: list[RaxmlTreeView] = []

    def tearDown(self) -> None:
        for view in self._views:
            view.close()
        self.app.processEvents()

    def _make_view(self, *, preserve_branch_lengths: bool = False) -> RaxmlTreeView:
        view = RaxmlTreeView(
            layout_mode="rectangular",
            preserve_branch_lengths=preserve_branch_lengths,
            optimize_groups=False,
            use_group_backgrounds=True,
            show_group_labels=True,
            show_cluster_annotations=True,
        )
        self._views.append(view)
        view.resize(1000, 700)
        view.set_reference_sample("E")
        view.set_leaf_groups(GROUPS, group_colors=GROUP_COLORS)
        view.set_leaf_clusters(CLUSTERS)
        view.set_tree_text(TREE_TEXT)
        view.show()
        self.app.processEvents()
        return view

    @staticmethod
    def _scene_items(view: RaxmlTreeView, kind: str) -> list[QtWidgets.QGraphicsItem]:
        return [item for item in view.scene().items() if item.data(_SCENE_DATA_KIND) == kind]

    def test_clade_bands_are_restrained_label_lane_accents(self) -> None:
        view = self._make_view()
        bands = self._scene_items(view, "clade-band")
        rails = self._scene_items(view, "clade-rail")
        labels = self._scene_items(view, "clade-label")
        self.assertEqual(len(bands), 2)
        self.assertEqual(len(rails), 2)
        self.assertEqual(len(labels), 2)

        for band in bands:
            group_name = str(band.data(_SCENE_DATA_NAME))
            run_index = int(band.data(_SCENE_DATA_RUN_INDEX))
            self.assertLessEqual(band.brush().color().alpha(), 20)
            self.assertEqual(band.pen().style(), QtCore.Qt.PenStyle.NoPen)
            group_leaves = [leaf for leaf, group in GROUPS.items() if group == group_name]
            self.assertGreater(
                band.sceneBoundingRect().left(),
                max(view._leaf_pos_scene[leaf].x() for leaf in group_leaves),
            )
            for leaf_name in group_leaves:
                leaf_point = view._leaf_pos_scene[leaf_name]
                topology_probe = QtCore.QPointF(
                    leaf_point.x() + 0.5,
                    leaf_point.y(),
                )
                self.assertFalse(band.contains(topology_probe))

            rail = next(
                item
                for item in rails
                if item.data(_SCENE_DATA_NAME) == group_name
                and int(item.data(_SCENE_DATA_RUN_INDEX)) == run_index
            )
            self.assertAlmostEqual(rail.rect().width(), _CLADE_RAIL_WIDTH)
            self.assertEqual(rail.brush().color().alpha(), 255)
            self.assertAlmostEqual(
                rail.sceneBoundingRect().height(),
                band.sceneBoundingRect().height(),
            )

        expected_labels = {"Clade 1 · 4 genomes", "Clade 2 · 3 genomes"}
        self.assertEqual({item.text() for item in labels}, expected_labels)
        self.assertEqual(
            {item.font().pointSizeF() for item in labels},
            {_CLADE_FONT_SIZE},
        )
        self.assertTrue(
            all(item.font().weight() == QtGui.QFont.Weight.DemiBold for item in labels)
        )

        distance_view = self._make_view(preserve_branch_lengths=True)
        for band in self._scene_items(distance_view, "clade-band"):
            group_name = str(band.data(_SCENE_DATA_NAME))
            group_leaves = [leaf for leaf, group in GROUPS.items() if group == group_name]
            self.assertGreater(
                band.sceneBoundingRect().left(),
                max(
                    distance_view._leaf_pos_scene[leaf].x()
                    for leaf in group_leaves
                ),
            )
            for leaf_name in group_leaves:
                leaf_point = distance_view._leaf_pos_scene[leaf_name]
                self.assertFalse(
                    band.contains(
                        QtCore.QPointF(leaf_point.x() + 0.5, leaf_point.y())
                    )
                )

    def test_split_cluster_labels_are_compact_and_keep_full_payload(self) -> None:
        view = self._make_view()
        labels = self._scene_items(view, "cluster-label")
        c3_labels = sorted(
            [item for item in labels if item.data(_SCENE_DATA_NAME) == "3"],
            key=lambda item: int(item.data(_SCENE_DATA_RUN_INDEX)),
        )
        self.assertEqual(len(c3_labels), 2)
        self.assertTrue(all(item.text() == "C3" for item in c3_labels))
        self.assertEqual(
            [int(item.data(_SCENE_DATA_RUN_COUNT)) for item in c3_labels],
            [2, 2],
        )
        self.assertIn("block 1 of 2", c3_labels[0].toolTip())
        self.assertIn("block 2 of 2", c3_labels[1].toolTip())
        self.assertTrue(all("4 genomes" in item.toolTip() for item in c3_labels))
        self.assertTrue(
            all(item.font().pointSizeF() == _CLUSTER_FONT_SIZE for item in labels)
        )

        singleton = next(item for item in labels if item.data(_SCENE_DATA_NAME) == "2")
        self.assertNotIn("block", singleton.toolTip())

        brackets = self._scene_items(view, "cluster-bracket")
        self.assertTrue(brackets)
        self.assertTrue(
            all(item.pen().color().name() == _TREE_CLUSTER_COLOR for item in brackets)
        )

        spy = QtTest.QSignalSpy(view.clusterClicked)
        device_rect = c3_labels[0].deviceTransform(view.viewportTransform()).mapRect(
            c3_labels[0].boundingRect()
        )
        QtTest.QTest.mouseClick(
            view.viewport(),
            QtCore.Qt.MouseButton.LeftButton,
            QtCore.Qt.KeyboardModifier.NoModifier,
            device_rect.center().toPoint(),
        )
        self.app.processEvents()
        self.assertEqual(len(spy), 1)
        payload = spy[0][0]
        self.assertEqual(payload["name"], "3")
        self.assertEqual(set(payload["samples"]), {"A", "B", "E", "F"})

    def test_tree_uses_semantic_ink_interaction_and_reference_colors(self) -> None:
        view = self._make_view()
        edge_items = [
            item for item in view.scene().items() if isinstance(item, _TreeEdgeItem)
        ]
        self.assertTrue(edge_items)
        self.assertTrue(
            all(item.pen().color().name() == _TREE_BRANCH_COLOR for item in edge_items)
        )
        self.assertEqual(view._leaf_label_items["A"].brush().color().name(), _TREE_LABEL_COLOR)
        self.assertEqual(
            view._leaf_label_items["E"].brush().color().name(),
            _TREE_REFERENCE_COLOR,
        )
        self.assertEqual(
            view._leaf_label_items["E"].font().weight(),
            QtGui.QFont.Weight.DemiBold,
        )
        self.assertEqual(
            {item.font().pointSizeF() for item in view._leaf_label_items.values()},
            {_LEAF_FONT_SIZE},
        )

        view.highlight_sample("A")
        # Hover is the shared amber (the Map's hovered countries), its
        # label the readable darker shade.
        self.assertEqual(
            view._leaf_label_items["A"].brush().color().name(),
            COLORS["hover_ink"],
        )
        self.assertTrue(
            all(
                item.pen().color().name() == COLORS["hover"]
                for item in view._edge_items_by_leaf["A"]
            )
        )

        view.highlight_sample(None)
        view.set_active_samples(["A"])
        self.assertEqual(
            view._leaf_label_items["B"].brush().color().name(),
            _TREE_INACTIVE_COLOR,
        )

    def test_topology_distance_toggle_preserves_order_highlight_and_view(self) -> None:
        view = self._make_view()
        # Every edge to scale (E's edge would be broken, see the long branch test).
        view.set_show_true_branch_lengths(True)
        topology_positions = {
            leaf: QtCore.QPointF(point)
            for leaf, point in view._leaf_pos_scene.items()
        }
        self.assertAlmostEqual(topology_positions["E"].x(), topology_positions["F"].x())

        view.zoom_by(1.2)
        manual_scale = view.transform().m11()
        view.highlight_sample("A")
        view.set_preserve_branch_lengths(True)
        self.app.processEvents()

        self.assertTrue(view.preserve_branch_lengths_enabled())
        self.assertGreater(view._leaf_pos_scene["E"].x(), view._leaf_pos_scene["F"].x())
        # Cumulative distances are E=10.2, F=0.3, G=0.2, so the renderer's
        # linear x scale must preserve the 10.0:0.1 separation ratio exactly.
        distance_ratio = (
            view._leaf_pos_scene["E"].x() - view._leaf_pos_scene["G"].x()
        ) / (
            view._leaf_pos_scene["F"].x() - view._leaf_pos_scene["G"].x()
        )
        self.assertAlmostEqual(distance_ratio, 100.0, places=6)

        scale_bars = self._scene_items(view, "distance-scale-bar")
        scale_labels = self._scene_items(view, "distance-scale-label")
        self.assertEqual(len(scale_bars), 1)
        self.assertEqual(len(scale_labels), 1)
        self.assertGreater(float(scale_bars[0].data(_SCENE_DATA_VALUE)), 0.0)
        self.assertIn("branch-length units", scale_labels[0].text())
        for leaf, point in topology_positions.items():
            self.assertAlmostEqual(view._leaf_pos_scene[leaf].y(), point.y())
        self.assertAlmostEqual(view.transform().m11(), manual_scale)
        self.assertEqual(view._active_leaf_set, {"A"})

        view.set_preserve_branch_lengths(False)
        self.app.processEvents()
        for leaf, point in topology_positions.items():
            self.assertAlmostEqual(view._leaf_pos_scene[leaf].x(), point.x())
            self.assertAlmostEqual(view._leaf_pos_scene[leaf].y(), point.y())
        self.assertFalse(self._scene_items(view, "distance-scale-bar"))

    def test_invalid_branch_lengths_cannot_masquerade_as_distance(self) -> None:
        view = RaxmlTreeView(
            layout_mode="rectangular",
            preserve_branch_lengths=True,
        )
        self._views.append(view)
        availability_spy = QtTest.QSignalSpy(view.branchDistanceAvailabilityChanged)
        mode_spy = QtTest.QSignalSpy(view.branchScaleModeChanged)

        view.set_tree_text("(A,B);")
        self.app.processEvents()

        self.assertFalse(view.branch_distance_available())
        self.assertFalse(view.preserve_branch_lengths_enabled())
        self.assertIn("missing", view.branch_distance_status_message())
        self.assertFalse(self._scene_items(view, "distance-scale-bar"))
        self.assertTrue(availability_spy)
        self.assertEqual(availability_spy[-1][0], False)
        self.assertTrue(mode_spy)
        self.assertEqual(mode_spy[-1][0], False)

    def test_clade_and_cluster_badges_remain_screen_stable_and_fit(self) -> None:
        view = self._make_view(preserve_branch_lengths=True)
        fixed_flag = (
            QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations
        )
        labels = self._scene_items(view, "clade-label") + self._scene_items(
            view,
            "cluster-label",
        )
        self.assertTrue(labels)
        self.assertTrue(all(item.flags() & fixed_flag for item in labels))
        behind_flag = QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemStacksBehindParent
        for label in labels:
            chips = [
                child
                for child in label.childItems()
                if child.data(_SCENE_DATA_KIND) in {"clade-chip", "cluster-chip"}
            ]
            self.assertEqual(len(chips), 1)
            self.assertTrue(chips[0].flags() & behind_flag)

        def device_size(item: QtWidgets.QGraphicsItem) -> QtCore.QSizeF:
            transform = item.deviceTransform(view.viewportTransform())
            return transform.mapRect(item.boundingRect()).size()

        sizes_before = [device_size(item) for item in labels]
        scale_before = view.transform().m11()
        # A stretch leaves them as they are; the camera zoom scales them.
        view.stretch_by(1.8, 1.8)
        self.app.processEvents()
        sizes_after = [device_size(item) for item in labels]
        self.assertGreater(view.transform().m11(), scale_before)
        for before, after in zip(sizes_before, sizes_after):
            self.assertAlmostEqual(before.width(), after.width(), delta=0.5)
            self.assertAlmostEqual(before.height(), after.height(), delta=0.5)
        view.zoom_by(1.5)
        self.app.processEvents()
        for before, after in zip(sizes_before, [device_size(item) for item in labels]):
            self.assertAlmostEqual(before.width() * 1.5, after.width(), delta=0.5)

        view.fit_to_window()
        self.app.processEvents()
        viewport_rect = QtCore.QRectF(view.viewport().rect())
        for label in labels:
            local_bounds = label.boundingRect().united(label.childrenBoundingRect())
            device_rect = label.deviceTransform(view.viewportTransform()).mapRect(
                local_bounds
            )
            self.assertGreaterEqual(device_rect.left(), viewport_rect.left() - 1.0)
            self.assertLessEqual(device_rect.right(), viewport_rect.right() + 1.0)

        image = QtGui.QImage(
            view.viewport().size(),
            QtGui.QImage.Format.Format_RGB32,
        )
        image.fill(QtGui.QColor("#ffffff"))
        painter = QtGui.QPainter(image)
        view.viewport().render(painter)
        painter.end()
        for label in labels:
            text_rect = label.deviceTransform(view.viewportTransform()).mapRect(
                label.boundingRect()
            ).toAlignedRect()
            dark_pixels = 0
            for y in range(max(0, text_rect.top()), min(image.height(), text_rect.bottom() + 1)):
                for x in range(max(0, text_rect.left()), min(image.width(), text_rect.right() + 1)):
                    color = image.pixelColor(x, y)
                    luminance = (
                        0.2126 * color.red()
                        + 0.7152 * color.green()
                        + 0.0722 * color.blue()
                    )
                    if luminance < 150.0:
                        dark_pixels += 1
            self.assertGreater(dark_pixels, 3, label.text())

    def test_clade_palette_is_calm_separate_and_naturally_stable(self) -> None:
        self.assertEqual(
            _FASTGEAR_LINEAGE_PALETTE[:5],
            ["#3f6ea5", "#35826f", "#7564a7", "#a97832", "#a45f73"],
        )
        self.assertNotEqual(
            _FASTGEAR_LINEAGE_PALETTE[:5],
            _METADATA_STRIP_PALETTE[:5],
        )
        first = _fastgear_lineage_colors(
            {"sample_10": "10", "sample_2": "2", "sample_1": "1"}
        )
        second = _fastgear_lineage_colors(
            {"sample_1": "1", "sample_10": "10", "sample_2": "2"}
        )
        self.assertEqual(list(first), ["1", "2", "10"])
        self.assertEqual(
            {key: color.name() for key, color in first.items()},
            {key: color.name() for key, color in second.items()},
        )
        self.assertEqual(first["1"].name(), _FASTGEAR_LINEAGE_PALETTE[0])
        self.assertEqual(first["2"].name(), _FASTGEAR_LINEAGE_PALETTE[1])
        self.assertEqual(first["10"].name(), _FASTGEAR_LINEAGE_PALETTE[9])

    def test_png_and_pdf_exports_include_the_modern_scene(self) -> None:
        view = self._make_view()
        fixed_flag = (
            QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations
        )
        stable_groups = list(view._device_stable_items)
        export_flag_states: list[list[bool]] = []
        original_export_image = view._export_image

        export_rects: list[QtCore.QRectF] = []

        def recording_export_image(*args: object, **kwargs: object) -> None:
            export_flag_states.append(
                [bool(item.flags() & fixed_flag) for item in stable_groups]
            )
            export_rects.append(QtCore.QRectF(args[1]))
            original_export_image(*args, **kwargs)

        view._export_image = recording_export_image
        with tempfile.TemporaryDirectory() as directory:
            png_path = Path(directory) / "tree.png"
            pdf_path = Path(directory) / "tree.pdf"
            view.export_figure(png_path, image_scale=1.0, image_dpi=96)
            view.export_figure(pdf_path, image_dpi=96)
            self.assertGreater(png_path.stat().st_size, 1_000)
            self.assertGreater(pdf_path.stat().st_size, 1_000)

            image = QtGui.QImage(str(png_path))
            self.assertFalse(image.isNull())
            rail = self._scene_items(view, "clade-rail")[0]
            # The exported rect (the view's scene rect also covers the
            # viewport when the fit leaves the tree smaller than it).
            center = rail.sceneBoundingRect().center() - export_rects[0].topLeft()
            pixel = image.pixelColor(
                max(0, min(image.width() - 1, int(round(center.x())))),
                max(0, min(image.height() - 1, int(round(center.y())))),
            )
            self.assertNotEqual(pixel.name(), "#ffffff")
        self.assertTrue(export_flag_states)
        self.assertTrue(all(not state for state in export_flag_states[0]))
        self.assertTrue(all(item.flags() & fixed_flag for item in stable_groups))

    def test_annotation_rich_toolbar_does_not_offer_lossy_circular_mode(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            tree_path = Path(directory) / "tree.nwk"
            tree_path.write_text(TREE_TEXT, encoding="utf-8")
            annotated = PanCoreAccessoryTab(
                Path(directory),
                None,
                reference_sample=None,
                tree_file=tree_path,
                preserve_branch_lengths=True,
                color_tree_by_fastgear_clade=True,
                annotate_fastgear_clusters=True,
                use_fastgear_clade_backgrounds=True,
                show_info_label=False,
            )
            try:
                self.assertIsNone(annotated._tree_layout_combo)
                self.assertIsNotNone(annotated.tree)
                self.assertEqual(annotated.tree.layout_mode(), "rectangular")
                # One Lengths dropdown: Distance (a tree without branch
                # length sets) or Topology.
                lengths = annotated._tree_lengths_combo
                self.assertIsNotNone(lengths)
                self.assertEqual(lengths.count(), 2)
                self.assertEqual(lengths.itemData(0), LENGTHS_DISTANCE)
                self.assertEqual(lengths.itemData(1), LENGTHS_TOPOLOGY)
                self.assertEqual(lengths.currentData(), LENGTHS_DISTANCE)
                self.assertTrue(lengths.isEnabled())
                self.assertEqual(
                    len(self._scene_items(annotated.tree, "distance-scale-bar")),
                    1,
                )
            finally:
                annotated.close()

            invalid_tree_path = Path(directory) / "tree_without_lengths.nwk"
            invalid_tree_path.write_text("(A,B);", encoding="utf-8")
            invalid = PanCoreAccessoryTab(
                Path(directory),
                None,
                reference_sample=None,
                tree_file=invalid_tree_path,
                preserve_branch_lengths=True,
                color_tree_by_fastgear_clade=True,
                annotate_fastgear_clusters=True,
                use_fastgear_clade_backgrounds=True,
                show_info_label=False,
            )
            try:
                self.assertEqual(invalid._tree_lengths_combo.currentData(), LENGTHS_TOPOLOGY)
                self.assertFalse(invalid._tree_lengths_combo.isEnabled())
                self.assertIn("missing", invalid._tree_lengths_combo.toolTip())
            finally:
                invalid.close()

            plain = PanCoreAccessoryTab(
                Path(directory),
                None,
                reference_sample=None,
                tree_file=None,
                preserve_branch_lengths=False,
                show_info_label=False,
            )
            try:
                self.assertIsNotNone(plain._tree_layout_combo)
                self.assertEqual(plain._tree_layout_combo.count(), 2)
            finally:
                plain.close()


# ---------------------------------------------------------------- phylo
TREE_LEAVES = ("A", "B", "C", "D", "E", "F", "G")


def _parse(text: str) -> _TreeNode:
    return _NewickParser(text).parse()


def _write_tree(test: unittest.TestCase, text: str) -> Path:
    folder = Path(tempfile.mkdtemp(prefix="bactolomew_tree_"))
    test.addCleanup(shutil.rmtree, folder, True)
    path = folder / "tree.nwk"
    path.write_text(text, encoding="utf-8")
    return path


def _leaves(node: _TreeNode) -> List[str]:
    if node.is_leaf:
        return [node.name]
    out: List[str] = []
    for child in node.children:
        out.extend(_leaves(child))
    return out


def _root_split(node: _TreeNode) -> Set[FrozenSet[str]]:
    """The leaf sets under the root's children (what the root position means)."""
    return {frozenset(_leaves(child)) for child in node.children}


def _depths(node: _TreeNode, base: float = 0.0) -> dict:
    if node.is_leaf:
        return {node.name: base}
    out: dict = {}
    for child in node.children:
        out.update(_depths(child, base + float(child.branch_length or 0.0)))
    return out


def _path_length(tree: _TreeNode, left: str, right: str) -> float:
    """Leaf to leaf distance through their most recent common ancestor."""

    def walk(node: _TreeNode, depth: float) -> Any:
        if node.is_leaf:
            found = {node.name: depth} if node.name in (left, right) else {}
            return found, None
        found: dict = {}
        result = None
        for child in node.children:
            child_found, child_result = walk(child, depth + float(child.branch_length or 0.0))
            found.update(child_found)
            result = result if result is not None else child_result
        if result is None and len(found) == 2:
            result = found[left] + found[right] - 2.0 * depth
        return found, result

    return float(walk(tree, 0.0)[1])


class PruneAndRootTests(unittest.TestCase):
    """Pruning (ape ``drop.tip``) and midpoint rooting (``phylo``)."""

    def test_prune_keeps_topology_and_sums_branch_lengths(self) -> None:
        tree = _parse(TREE_TEXT)
        pruned = _prune_to(tree, {"A", "C", "D", "G"})
        self.assertIsNotNone(pruned)
        assert pruned is not None
        self.assertEqual(sorted(_leaves(pruned)), ["A", "C", "D", "G"])
        # (A,B) lost B: A's edge absorbs the merged 0.3 edge above it.
        self.assertEqual(_root_split(pruned), {frozenset({"A", "C", "D"}), frozenset({"G"})})
        depths = _depths(pruned)
        self.assertAlmostEqual(depths["A"], 0.1 + 0.3 + 0.1)
        self.assertAlmostEqual(depths["C"], 0.1 + 0.2 + 0.1)
        # G's edge absorbs the 0.1 edge of the (E,F,G) node.
        self.assertAlmostEqual(depths["G"], 0.1 + 0.1)
        # Distances between kept leaves are unchanged.
        for left, right in (("A", "C"), ("C", "D"), ("A", "G"), ("D", "G")):
            self.assertAlmostEqual(_path_length(pruned, left, right), _path_length(tree, left, right))
        # The file tree itself is untouched.
        self.assertEqual(sorted(_leaves(tree)), list(TREE_LEAVES))

    def test_prune_moves_the_root_to_the_kept_leaves_ancestor(self) -> None:
        tree = _parse(TREE_TEXT)
        pruned = _prune_to(tree, {"C", "D"})
        assert pruned is not None
        self.assertEqual(_root_split(pruned), {frozenset({"C"}), frozenset({"D"})})
        single = _prune_to(tree, {"E"})
        assert single is not None
        self.assertTrue(single.is_leaf)
        self.assertEqual(single.name, "E")
        self.assertIsNone(_prune_to(tree, {"Z"}))
        self.assertIsNone(_prune_to(tree, set()))

    def test_midpoint_root_on_a_small_newick(self) -> None:
        # Longest path A..C = 8; its midpoint lies 2 past the file root on C's edge.
        rooted = _midpoint_root(_parse("((A:1,B:1):1,C:6);"))
        self.assertEqual(_root_split(rooted), {frozenset({"A", "B"}), frozenset({"C"})})
        depths = _depths(rooted)
        self.assertAlmostEqual(depths["A"], 4.0)
        self.assertAlmostEqual(depths["B"], 4.0)
        self.assertAlmostEqual(depths["C"], 4.0)

        tree = _parse(TREE_TEXT)
        rooted = _midpoint_root(tree)
        depths = _depths(rooted)
        # E (10.2 deep) and B (0.6 deep) are the farthest pair (10.8 apart):
        # both sit 5.4 from the midpoint, and no leaf is farther.
        self.assertAlmostEqual(depths["E"], 5.4)
        self.assertAlmostEqual(depths["B"], 5.4)
        self.assertLessEqual(max(depths.values()), 5.4 + 1e-9)
        self.assertEqual(_root_split(rooted), {frozenset({"E"}), frozenset(set(TREE_LEAVES) - {"E"})})
        for left, right in (("A", "G"), ("E", "F"), ("C", "B")):
            self.assertAlmostEqual(_path_length(rooted, left, right), _path_length(tree, left, right))

    def test_scale_bar_value_is_plain_decimal(self) -> None:
        fmt = RaxmlTreeView._format_branch_scale_value
        self.assertEqual(fmt(2e-06), "0.000002")
        self.assertEqual(fmt(0.0005), "0.0005")
        self.assertEqual(fmt(0.002), "0.002")
        self.assertEqual(fmt(5.0), "5")
        self.assertEqual(fmt(50.0), "50")
        self.assertEqual(fmt(3 * 1e-07), "0.0000003")


class RootAndPruneViewTests(unittest.TestCase):
    """RaxmlTreeView's root modes, pruning and scale bar unit, and a tab
    that prunes its canonical tree."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def setUp(self) -> None:
        self._views: List[QtWidgets.QWidget] = []
        self._dir = tempfile.mkdtemp(prefix="bactolomew_tree_view_")

    def tearDown(self) -> None:
        for view in self._views:
            view.close()
            view.deleteLater()
        self.app.processEvents()
        shutil.rmtree(self._dir, ignore_errors=True)

    def _view(self, **kwargs: Any) -> RaxmlTreeView:
        options = dict(layout_mode="rectangular", preserve_branch_lengths=True)
        options.update(kwargs)
        view = RaxmlTreeView(**options)
        self._views.append(view)
        view.resize(900, 700)
        return view

    @staticmethod
    def _items(view: RaxmlTreeView, kind: str) -> list:
        return [item for item in view.scene().items() if item.data(_SCENE_DATA_KIND) == kind]

    @staticmethod
    def _texts(view: RaxmlTreeView) -> List[str]:
        return [item.text() for item in view.scene().items() if isinstance(item, QtWidgets.QGraphicsSimpleTextItem)]

    def test_filter_never_yields_zero_leaves_when_samples_exist(self) -> None:
        view = self._view(prune_inactive=True)
        view.set_tree_text(TREE_TEXT)
        subsets = [
            ["A"],
            ["E", "G"],
            ["A", "C", "F"],
            ["B", "D", "E", "G"],
            list(TREE_LEAVES),
            # Names outside the tree are ignored; the ones inside stay.
            ["A", "not_in_tree", "F"],
        ]
        for subset in subsets:
            view.set_active_samples(subset)
            expected = sorted(name for name in subset if name in TREE_LEAVES)
            self.assertEqual(sorted(view.leaf_order()), expected, subset)
            self.assertTrue(view.leaf_order())
            self.assertFalse(any("Failed" in text or "No tree" in text for text in self._texts(view)))
        # Without a filter every leaf is back.
        view.set_active_samples(None)
        self.assertEqual(sorted(view.leaf_order()), list(TREE_LEAVES))
        # An empty filter says so instead of drawing a blank tree.
        view.set_active_samples([])
        self.assertEqual(view.leaf_order(), [])
        self.assertIn("No genomes match the current filters.", self._texts(view))

    def test_genomes_outside_the_tree_say_so(self) -> None:
        view = self._view(prune_inactive=True)
        view.set_tree_text(TREE_TEXT)
        view.set_active_samples(["Z"])
        self.assertEqual(view.leaf_order(), [])
        self.assertIn("The filtered genome is not in the tree.", self._texts(view))
        view.set_active_samples(["Y", "Z"])
        self.assertIn("None of the 2 filtered genomes are in the tree.", self._texts(view))
        self.assertFalse(any("Failed" in text for text in self._texts(view)))
        # The message keeps its own size (no fit scaling).
        self.assertAlmostEqual(view.transform().m11(), 1.0)
        view.set_active_samples(["A"])
        self.assertEqual(view.leaf_order(), ["A"])

    def test_clearing_the_filter_redraws_the_whole_tree(self) -> None:
        # The render cache must not keep the pruned tree once the filter is gone.
        view = self._view(prune_inactive=True, reroot_for_groups=False, optimize_groups=True)
        view.set_leaf_groups(
            {"A": "1", "B": "1", "C": "2", "D": "2", "E": "1", "F": "2", "G": "2"},
            group_colors={"1": "#3f6ea5", "2": "#35826f"},
        )
        view.set_tree_text(TREE_TEXT)
        view.set_active_samples(["A", "C", "E"])
        self.assertEqual(sorted(view.leaf_order()), ["A", "C", "E"])
        view.set_active_samples(None)
        self.assertEqual(sorted(view.leaf_order()), list(TREE_LEAVES))
        view.set_active_samples(["A", "C", "E"])
        self.assertEqual(sorted(view.leaf_order()), ["A", "C", "E"])
        view.set_prune_inactive(False)
        self.assertEqual(sorted(view.leaf_order()), list(TREE_LEAVES))

    def test_midpoint_root_is_kept_when_pruning(self) -> None:
        # Midpoint of the full tree splits E from the rest; pruning to
        # A, B, C, E keeps that root (the subset's own midpoint would differ
        # only by position, so check the split and the depths too).
        view = self._view(prune_inactive=True)
        view.set_tree_text(TREE_TEXT)
        view.set_root_mode("midpoint")
        view.set_active_samples(["A", "C", "E", "F"])
        rendered = view._resolved_render_root()
        assert rendered is not None
        expected = _prune_to(_midpoint_root(_parse(TREE_TEXT)), {"A", "C", "E", "F"})
        assert expected is not None
        self.assertEqual(_root_split(rendered), _root_split(expected))
        self.assertEqual(_root_split(rendered), {frozenset({"E"}), frozenset({"A", "C", "F"})})
        for leaf, depth in _depths(expected).items():
            self.assertAlmostEqual(_depths(rendered)[leaf], depth)
        # The subset's own midpoint (E to A, 10.7 long) would put E 5.35 deep.
        self.assertAlmostEqual(_depths(rendered)["E"], 5.4)

    def test_greyed_mode_keeps_every_leaf(self) -> None:
        # Pan tree (no prune) keeps greying the filtered-out leaves.
        view = self._view()
        view.set_tree_text(TREE_TEXT)
        view.set_active_samples(["A", "B"])
        self.assertEqual(sorted(view.leaf_order()), list(TREE_LEAVES))

    def test_displayed_root_equals_the_file_root(self) -> None:
        # Groups split by the file root, so a reroot would make them contiguous.
        text = "((A:1,(B:1,C:1):1):1,((D:1,E:1):1,F:1):1);"
        groups = {"A": "1", "B": "2", "C": "2", "D": "2", "E": "2", "F": "1"}
        view = self._view(optimize_groups=True, use_group_backgrounds=True, reroot_for_groups=False)
        view.set_leaf_groups(groups, group_colors={"1": "#3f6ea5", "2": "#35826f"})
        view.set_tree_text(text)
        rendered = view._resolved_render_root()
        assert rendered is not None
        self.assertEqual(_root_split(rendered), _root_split(_parse(text)))
        # Pruned and rotated, the root is still where the file has it.
        view.set_active_samples(["A", "B", "D", "F"])
        view.set_prune_inactive(True)
        rendered = view._resolved_render_root()
        assert rendered is not None
        self.assertEqual(_root_split(rendered), {frozenset({"A", "B"}), frozenset({"D", "F"})})

    def test_midpoint_root_mode(self) -> None:
        view = self._view()
        view.set_tree_text("((A:1,B:1):1,C:6);")
        view.set_root_mode("midpoint")
        rendered = view._resolved_render_root()
        assert rendered is not None
        self.assertAlmostEqual(_depths(rendered)["C"], 4.0)
        view.set_root_mode("file")
        self.assertEqual(view._resolved_render_root(), view._root)
        with self.assertRaises(ValueError):
            view.set_root_mode("outgroup")

    def test_scale_bar_names_its_unit_under_the_tree(self) -> None:
        view = self._view(scale_bar_position="bottom")
        tree_path = Path(self._dir) / "unit.nwk"
        tree_path.write_text(TREE_TEXT, encoding="utf-8")
        view.load_tree_file(tree_path, root_mode="file", branch_length_unit="substitutions/site")
        labels = self._items(view, "distance-scale-label")
        bars = self._items(view, "distance-scale-bar")
        self.assertEqual(len(labels), 1)
        self.assertIn(" substitutions/site", labels[0].text())
        lowest_leaf = max(point.y() for point in view._leaf_pos_scene.values())
        leftmost_leaf = min(point.x() for point in view._leaf_pos_scene.values())
        self.assertGreater(bars[0].sceneBoundingRect().top(), lowest_leaf)
        self.assertLessEqual(bars[0].sceneBoundingRect().left(), leftmost_leaf)

    def test_tab_prunes_the_tree_file_it_is_given(self) -> None:
        tree_path = Path(self._dir) / "core.snippy_raxml.bestTree"
        tree_path.write_text(TREE_TEXT, encoding="utf-8")
        tab = PanCoreAccessoryTab(
            Path(self._dir),
            None,
            reference_sample=None,
            tree_file=tree_path,
            preserve_branch_lengths=True,
            show_info_label=False,
            prune_tree_to_sample_filter=True,
            reroot_tree_for_clades=False,
            scale_bar_position="bottom",
        )
        self._views.append(tab)
        self.assertEqual(tab.current_tree_source().path, tree_path)
        self.assertEqual(tab.tree_root_mode(), "file")
        tree = tab.tree
        assert tree is not None
        tab.set_sample_filter(["B", "C", "G"])
        self.assertEqual(sorted(tree.leaf_order()), ["B", "C", "G"])
        tab.set_sample_filter(None)
        self.assertEqual(sorted(tree.leaf_order()), list(TREE_LEAVES))
        tab.set_tree_root_mode("midpoint")
        self.assertEqual(tree.root_mode(), "midpoint")
        self.assertEqual(sorted(tree.leaf_order()), list(TREE_LEAVES))

    def test_panaroo_tree_draws_generic_units(self) -> None:
        tab = PanCoreAccessoryTab(
            Path(self._dir),
            None,
            reference_sample=None,
            tree_file=None,
            preserve_branch_lengths=True,
            show_info_label=False,
        )
        self._views.append(tab)
        panaroo_tree = Path(self._dir) / "panaroo_tree.nwk"
        panaroo_tree.write_text(TREE_TEXT, encoding="utf-8")
        tab._find_panaroo_tree = lambda: panaroo_tree  # type: ignore[method-assign]
        tab._tree_index_cache = None
        tab._loaded_tree_key = None
        tab._apply_tree_source()
        source = tab.current_tree_source()
        self.assertEqual((source.id, source.path), ("panaroo", panaroo_tree))
        # An old, unrooted file is drawn at its midpoint.
        self.assertEqual(tab.tree_root_mode(), "midpoint")
        tree = tab.tree
        assert tree is not None
        labels = self._items(tree, "distance-scale-label")
        self.assertTrue(labels and labels[0].text().endswith(" branch-length units"), labels)

    def test_branch_length_set_keeps_rows_root_and_support(self) -> None:
        """Another fit's branch lengths (Phylogeny "Branch lengths") change
        only the drawn lengths: the leaf order (every root and order mode,
        pruned too), the support hover and the root edge stay; a tree with
        another topology is refused."""
        main = Path(self._dir) / "mutations.nwk"
        main.write_text("(((A:1,B:1)90:1,C:3)85:1,(D:1,E:1)70:4);", encoding="utf-8")
        other = Path(self._dir) / "with_recombination.nwk"
        other.write_text("(((A:6,B:1)90:9,C:0.5)85:1,(D:1,E:8)70:0.2);", encoding="utf-8")
        view = self._view()
        view.load_tree_file(main, root_mode="file", support_kind="UFBoot", support_threshold=95)
        for prune in (None, ["A", "C", "D", "E"]):
            view.set_prune_inactive(prune is not None)
            view.set_active_samples(prune)
            for root in ("file", "midpoint"):
                view.set_root_mode(root)
                for order in phylo.ORDER_MODES:
                    view.set_order_mode(order)
                    self.assertTrue(view.set_branch_lengths(None))
                    before = view.leaf_order()
                    self.assertTrue(view.set_branch_lengths(other, label="Including recombination"))
                    self.assertEqual(view.leaf_order(), before, (prune, root, order))
                    self.assertEqual(view.support_text(["A", "B"]), "UFBoot 90, not well supported (< 95)")
        view.set_prune_inactive(False)
        view.set_active_samples(None)
        view.set_root_mode("file")
        view.set_order_mode(phylo.ORDER_AS_COMPUTED)
        x = {leaf: point.x() for leaf, point in view._leaf_pos_scene.items()}
        # Edge lengths of the other fit: A-B 6 vs 1, E-D 8 vs 1.
        self.assertAlmostEqual((x["E"] - x["D"]) / (x["A"] - x["B"]), 7.0 / 5.0, places=6)
        self.assertTrue(any("including recombination" in text for text in self._texts(view)))
        # The root stays on the same edge (A, B, C | D, E), splitting the
        # other fit's 1.2 as the tree's own 1 | 4: A at 0.24 + 9 + 6, D at 0.96 + 1.
        self.assertAlmostEqual((x["A"] - x["D"]) / (x["E"] - x["D"]), (15.24 - 1.96) / 7.0, places=6)
        other.write_text("(((A:1,C:1)90:1,B:1)85:1,(D:1,E:1)70:4);", encoding="utf-8")
        self.assertFalse(view.set_branch_lengths(other, label="Including recombination", own_label="Mutations only"))
        self.assertFalse(view.branch_lengths_replaced())
        self.assertIn("different tree", view.branch_lengths_error())
        # The scale bar names the lengths drawn: the tree's own.
        self.assertFalse(any("including recombination" in text for text in self._texts(view)))
        self.assertTrue(any("mutations only" in text for text in self._texts(view)))


GROUP_SPLIT_TEXT = "((A:1,(B:1,C:1):1):1,((D:1,E:1):1,F:1):1);"
GROUP_SPLIT_GROUPS = {"A": "1", "B": "2", "C": "2", "D": "2", "E": "2", "F": "1"}

def _runs(order: List[str], groups: dict) -> int:
    return sum(1 for index, name in enumerate(order) if index == 0 or groups[name] != groups[order[index - 1]])


def _no_unary(node: _TreeNode) -> bool:
    return node.is_leaf or (len(node.children) >= 2 and all(_no_unary(child) for child in node.children))


def _group_rooted(text: str, groups: dict) -> _TreeNode:
    """The tree rooted on the edge that keeps ``groups`` together, with the
    one-child node a bifurcating file root leaves merged away."""
    return merge_unary(group_root(_parse(text), lambda name: groups.get(name, ""), str))


class GroupRootTests(unittest.TestCase):
    """The clade reroot over every edge (``phylo.group_root``) that the tree
    panels use for the "groups" root."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def test_group_root_makes_the_groups_contiguous(self) -> None:
        tree = _parse(GROUP_SPLIT_TEXT)
        # The file root splits both groups; rotation alone cannot join them.
        rotated = clone_tree(tree)
        GroupOptimizer(lambda name: GROUP_SPLIT_GROUPS.get(name, ""), str).optimize(rotated)
        self.assertGreater(_runs(_leaves(rotated), GROUP_SPLIT_GROUPS), 2)
        self.assertFalse(groups_contiguous(rotated, GROUP_SPLIT_GROUPS, str))
        rooted = _group_rooted(GROUP_SPLIT_TEXT, GROUP_SPLIT_GROUPS)
        self.assertEqual(_runs(_leaves(rooted), GROUP_SPLIT_GROUPS), 2)
        self.assertTrue(groups_contiguous(rooted, GROUP_SPLIT_GROUPS, str))
        self.assertTrue(_no_unary(rooted))
        # The file tree is never changed.
        group_root(tree, lambda name: GROUP_SPLIT_GROUPS.get(name, ""), str)
        self.assertEqual(_root_split(tree), _root_split(_parse(GROUP_SPLIT_TEXT)))

    def test_a_filter_keeps_the_full_tree_group_root(self) -> None:
        keep = {"A", "B", "D", "F"}
        expected = _prune_to(_group_rooted(GROUP_SPLIT_TEXT, GROUP_SPLIT_GROUPS), keep)
        assert expected is not None
        # The subset's own group root would differ: the full tree's is kept.
        view = RaxmlTreeView(
            layout_mode="rectangular", preserve_branch_lengths=True, optimize_groups=True,
            reroot_for_groups=False, root_mode="groups", prune_inactive=True,
        )
        self.addCleanup(view.deleteLater)
        view.set_leaf_groups(GROUP_SPLIT_GROUPS)
        view.set_tree_text(GROUP_SPLIT_TEXT)
        view.set_active_samples(sorted(keep))
        drawn = view._resolved_render_root()
        assert drawn is not None
        self.assertEqual(_root_split(drawn), _root_split(expected))

    def test_fewer_than_two_groups_keep_the_file_root(self) -> None:
        tree = _parse(GROUP_SPLIT_TEXT)
        for groups in ({}, {name: "1" for name in GROUP_SPLIT_GROUPS}):
            rooted = group_root(tree, lambda name, groups=groups: groups.get(name, ""), str)
            self.assertEqual(_root_split(rooted), _root_split(tree), groups)

    def test_a_bifurcating_file_root_leaves_no_unary_node(self) -> None:
        text = "((A:1,(B:1,C:1):1):1,(((D:1,E:1):1,F:1):1,G:2):1);"
        groups = {"A": "1", "B": "2", "C": "2", "D": "2", "E": "2", "F": "1", "G": "1"}
        rooted = _group_rooted(text, groups)
        self.assertTrue(_no_unary(rooted))
        self.assertEqual(sorted(_leaves(rooted)), sorted(groups))
        self.assertEqual(_runs(_leaves(rooted), groups), 2)
        # Leaf to leaf distances are unchanged by the reroot and the merge.
        tree = _parse(text)
        for left, right in (("A", "G"), ("B", "F"), ("D", "E"), ("C", "G")):
            self.assertAlmostEqual(_path_length(rooted, left, right), _path_length(tree, left, right))

    def test_root_mode_groups_order_equals_the_clade_reroot_order(self) -> None:
        """One search: root_mode "groups" gives the order of the clade
        reroot (neoPan), with the snippy "Reference" leaf too."""
        text = "((Reference:1,(B:1,C:1):1):1,((D:1,E:1):1,F:1):1);"
        groups = {"REF1": "1", "B": "2", "C": "2", "D": "2", "E": "2", "F": "1"}
        pan = RaxmlTreeView(layout_mode="rectangular", preserve_branch_lengths=True, optimize_groups=True, reroot_for_groups=True)
        by_mode = RaxmlTreeView(
            layout_mode="rectangular", preserve_branch_lengths=True, optimize_groups=True,
            reroot_for_groups=False, root_mode="groups",
        )
        for view in (pan, by_mode):
            self.addCleanup(view.deleteLater)
            view.set_reference_sample("REF1")
            view.set_leaf_groups(groups)
            view.set_tree_text(text)
        self.assertEqual(by_mode.leaf_order(), pan.leaf_order())
        self.assertIn("REF1", [str(name) for name in pan.leaf_order()])
        self.assertEqual(_runs([str(name) for name in pan.leaf_order()], groups), 2)

    def test_phylogeny_view_group_root_does_not_move_with_the_filter(self) -> None:
        """RaxmlTreeView with root_mode "groups" and pruning roots the full
        tree, then prunes it; the rows keep the full order."""
        keep = {"A", "B", "D", "F"}
        view = RaxmlTreeView(
            layout_mode="rectangular", preserve_branch_lengths=True, optimize_groups=True,
            reroot_for_groups=False, root_mode="groups", prune_inactive=True,
        )
        self.addCleanup(view.deleteLater)
        view.set_leaf_groups(GROUP_SPLIT_GROUPS)
        view.set_tree_text(GROUP_SPLIT_TEXT)
        full = list(view.leaf_order())
        view.set_active_samples(sorted(keep))
        drawn = view._resolved_render_root()
        expected = _prune_to(_group_rooted(GROUP_SPLIT_TEXT, GROUP_SPLIT_GROUPS), keep)
        self.assertEqual(_root_split(drawn), _root_split(expected))
        self.assertEqual(list(view.leaf_order()), [name for name in full if name in keep])


class PhyloLayoutTests(unittest.TestCase):
    """The Qt free tree preparation in ``phylo`` (root_tree, prepare_tree,
    the long branch rule) against the Phylogeny panel's drawing."""

    def test_root_tree_groups_and_fallback(self) -> None:
        tree = _parse(GROUP_SPLIT_TEXT)
        rotated = phylo.prepare_tree(tree, leaf_name=str, root_mode=phylo.ROOT_AS_IN_FILE, groups=GROUP_SPLIT_GROUPS)
        self.assertGreater(_runs(_leaves(rotated), GROUP_SPLIT_GROUPS), 2)
        rooted, drawn = phylo.root_tree(tree, leaf_name=str, root_mode=phylo.ROOT_GROUPS, groups=GROUP_SPLIT_GROUPS)
        self.assertEqual(drawn, phylo.ROOT_GROUPS)
        self.assertEqual(_leaves(rooted), _leaves(_group_rooted(GROUP_SPLIT_TEXT, GROUP_SPLIT_GROUPS)))
        self.assertTrue(_no_unary(rooted))
        self.assertEqual(_root_split(tree), _root_split(_parse(GROUP_SPLIT_TEXT)))
        prepared = phylo.prepare_tree(tree, leaf_name=str, root_mode=phylo.ROOT_GROUPS, groups=GROUP_SPLIT_GROUPS)
        self.assertEqual(_leaves(prepared), _leaves(rooted))
        # Fewer than two clades: the file tree itself, drawn as in the file.
        for groups in ({}, {name: "1" for name in GROUP_SPLIT_GROUPS}):
            same, drawn = phylo.root_tree(tree, leaf_name=str, root_mode=phylo.ROOT_GROUPS, groups=groups)
            self.assertEqual(drawn, phylo.ROOT_AS_IN_FILE)
            self.assertIs(same, tree)

    def test_prepare_tree_prunes_after_rooting_the_full_tree(self) -> None:
        tree = _parse(GROUP_SPLIT_TEXT)
        keep = {"A", "B", "D", "F"}
        rooted, _drawn = phylo.root_tree(tree, leaf_name=str, root_mode=phylo.ROOT_GROUPS, groups=GROUP_SPLIT_GROUPS)
        expected = _prune_to(rooted, keep)
        prepared = phylo.prepare_tree(tree, leaf_name=str, root_mode=phylo.ROOT_GROUPS, keep=keep, groups=GROUP_SPLIT_GROUPS)
        self.assertEqual(_root_split(prepared), _root_split(expected))
        # A pruned file rooted tree is rotated again only when clades split.
        full = phylo.prepare_tree(tree, leaf_name=str, root_mode=phylo.ROOT_AS_IN_FILE, groups=GROUP_SPLIT_GROUPS)
        self.assertFalse(groups_contiguous(full, GROUP_SPLIT_GROUPS, str))
        keep = {"A", "B", "C", "F"}
        pruned = phylo.prepare_tree(tree, leaf_name=str, root_mode=phylo.ROOT_AS_IN_FILE, keep=keep, groups=GROUP_SPLIT_GROUPS)
        self.assertTrue(groups_contiguous(pruned, GROUP_SPLIT_GROUPS, str))
        self.assertEqual(sorted(_leaves(pruned)), sorted(keep))
        self.assertIsNone(phylo.prepare_tree(tree, leaf_name=str, keep={"Z"}))
        midpoint = phylo.prepare_tree(_parse(TREE_TEXT), leaf_name=str, root_mode=phylo.ROOT_MIDPOINT)
        self.assertIn("E", {child.name for child in midpoint.children if child.is_leaf})

    def test_ladderize_only_rotates_and_keeps_every_band_together(self) -> None:
        import random

        rng = random.Random(7)
        nodes = [_TreeNode(name=f"L{i:02d}", branch_length=rng.random()) for i in range(60)]
        while len(nodes) > 1:
            a = nodes.pop(rng.randrange(len(nodes)))
            b = nodes.pop(rng.randrange(len(nodes)))
            nodes.append(_TreeNode(children=[a, b], branch_length=rng.random(), support=rng.randrange(100)))
        tree = nodes[0]
        order = _leaves(tree)
        # Lineages contiguous in the computed order; cluster "c" in two runs
        # (not monophyletic), the other clusters contiguous; a few genomes
        # without a cluster (one inside "c"), which must not be rotated into
        # another cluster's strip.
        lineages = {name: "1" if index < 25 else "2" for index, name in enumerate(order)}
        clusters = {
            name: "a" if index < 10 else "c" if index < 20 or 40 <= index < 50 else "b"
            for index, name in enumerate(order)
            if index not in (30, 31, 45)
        }

        def clades(node: _TreeNode) -> dict:
            out: dict = {}

            def walk(current: _TreeNode) -> frozenset:
                leaves = frozenset([current.name]) if current.is_leaf else frozenset().union(*map(walk, current.children))
                out[leaves] = (current.branch_length, current.support)
                return leaves

            walk(node)
            return out

        def run_counts(names: List[str], groups: dict) -> dict:
            counts: dict = {}
            for index, name in enumerate(names):
                group = groups.get(name, "")
                if group and (index == 0 or group != groups.get(names[index - 1], "")):
                    counts[group] = counts.get(group, 0) + 1
            return counts

        for mode in (phylo.ORDER_LADDER_SIZE, phylo.ORDER_LADDER_DEPTH, phylo.ORDER_AS_COMPUTED):
            drawn = phylo.ladderize(tree, leaf_name=str, mode=mode, keep_together=[lineages, clusters])
            names = _leaves(drawn)
            self.assertEqual(clades(drawn), clades(tree), mode)
            self.assertEqual(run_counts(names, lineages), {"1": 1, "2": 1}, mode)
            before = run_counts(order, clusters)
            for cluster, runs in run_counts(names, clusters).items():
                self.assertLessEqual(runs, before[cluster], (mode, cluster))
            if mode != phylo.ORDER_AS_COMPUTED:
                self.assertNotEqual(names, order, mode)
        self.assertEqual(_leaves(tree), order)
        # Without bands to keep: every node has its larger clade on top.
        free = phylo.ladderize(tree, leaf_name=str)

        def size(node: _TreeNode) -> int:
            return 1 if node.is_leaf else sum(size(child) for child in node.children)

        stack = [free]
        while stack:
            node = stack.pop()
            if not node.is_leaf:
                self.assertGreaterEqual(size(node.children[0]), size(node.children[1]))
                stack.extend(node.children)

    def test_layout_matches_the_tree_view_order(self) -> None:
        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        del app
        view = RaxmlTreeView(layout_mode="rectangular", preserve_branch_lengths=True, optimize_groups=True, reroot_for_groups=False)
        self.addCleanup(view.deleteLater)
        view.set_leaf_groups(GROUPS, group_colors=GROUP_COLORS)
        view.set_tree_text(TREE_TEXT)
        root = _parse(TREE_TEXT)
        prepared = phylo.prepare_tree(root, leaf_name=str, groups=GROUPS)
        self.assertEqual(_leaves(prepared), [str(name) for name in view.leaf_order()])
        self.assertEqual(_leaves(root), list(TREE_LEAVES))

    def test_long_branch_is_drawn_short_with_its_true_length(self) -> None:
        # One edge (E) 300 times the rest: drawn short with a break mark and
        # its true length; the rows stay as the unbroken drawing has them.
        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        del app
        text = "((A:1,B:1):1,((C:1,D:1):1,(E:300,F:1):1):1);"
        view = RaxmlTreeView(layout_mode="rectangular", preserve_branch_lengths=True)
        self.addCleanup(view.deleteLater)
        view.load_tree_file(_write_tree(self, text), root_mode="file", branch_length_unit="substitutions/site")
        self.assertEqual(view.broken_branches(), [(("E",), 300.0)])
        broken_order = view.leaf_order()
        def ratio() -> float:
            # At true length E is 300 branch units right of A, F one unit.
            x = {leaf: view._leaf_pos_scene[leaf].x() for leaf in ("A", "E", "F")}
            return (x["E"] - x["A"]) / (x["F"] - x["A"])

        self.assertLess(ratio(), 5.0)
        marks = [item for item in view.scene().items() if item.data(_SCENE_DATA_KIND) == "branch-break-label"]
        self.assertEqual([item.text() for item in marks], ["300 subst/site"])
        self.assertIn("true length 300 substitutions/site", view.branch_break_note(["E"]))
        scale = [item for item in view.scene().items() if item.data(_SCENE_DATA_KIND) == "distance-scale-label"]
        self.assertIn("broken branches not to scale", scale[0].text())
        # The toggle draws every edge at its true length, rows unchanged.
        view.set_show_true_branch_lengths(True)
        self.assertEqual(view.broken_branches(), [])
        self.assertTrue(view.branch_breaks_available())
        self.assertEqual(view.leaf_order(), broken_order)
        self.assertAlmostEqual(ratio(), 300.0, places=6)
        self.assertFalse([item for item in view.scene().items() if item.data(_SCENE_DATA_KIND) == "branch-break"])
        # A midpoint root on a long edge: both root edges are broken.
        two = phylo.long_branch_breaks(_parse("((A:1,B:2):100,(C:1,(D:1,E:3):1):90);"), 10.0)
        self.assertEqual(sorted(two.true_lengths), [90.0, 100.0])
        # A balanced tree breaks nothing.
        self.assertFalse(phylo.long_branch_breaks(_parse("((A:1,B:1):1,(C:1,D:1):1);"), 10.0))


class TreeHoverEmphasisTests(unittest.TestCase):
    """A hovered sample stands out (feedback 2026-10-08): a thick amber
    cosmetic path above the other branches, a tip dot and a name chip in
    viewport pixels, and the rest of the tree dimmed, with no rebuild."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def setUp(self) -> None:
        self._views: List[RaxmlTreeView] = []

    def tearDown(self) -> None:
        for view in self._views:
            view.close()
            view.deleteLater()
        self.app.processEvents()

    def _make_view(self, layout_mode: str = "rectangular", size=(600, 400)) -> RaxmlTreeView:
        view = RaxmlTreeView(layout_mode=layout_mode, use_group_backgrounds=False)
        self._views.append(view)
        view.resize(*size)
        view.set_leaf_groups(GROUPS, group_colors=GROUP_COLORS)
        view.set_leaf_clusters(CLUSTERS)
        view.set_tree_text(TREE_TEXT)
        view.show()
        self.app.processEvents()
        return view

    @staticmethod
    def _edges(view: RaxmlTreeView) -> List[_TreeEdgeItem]:
        return [item for item in view.scene().items() if isinstance(item, _TreeEdgeItem)]

    def test_the_hover_pen_is_cosmetic_amber_and_at_least_3_px(self) -> None:
        view = self._make_view()
        view.highlight_sample("A")
        for item in view._edge_items_by_leaf["A"]:
            pen = item.pen()
            self.assertTrue(pen.isCosmetic())
            self.assertGreaterEqual(pen.widthF(), 3.0)
            self.assertEqual(pen.color().name(), COLORS["hover"])
        self.assertEqual(_HOVER_PATH_WIDTH_PX, 3.0)

    def test_the_hovered_path_draws_above_the_other_edges(self) -> None:
        view = self._make_view()
        view.highlight_sample("E")
        hovered = set(view._edge_items_by_leaf["E"])
        others = [item for item in self._edges(view) if item not in hovered]
        self.assertTrue(others)
        top_other = max(item.zValue() for item in others)
        parents = {item.parentItem() for item in self._edges(view)}
        self.assertEqual(len(parents), 1)  # siblings: z decides the order
        self.assertTrue(all(item.zValue() > top_other for item in hovered))
        # In the scene's stacking order every hovered edge is above the rest.
        order = [item for item in view.scene().items() if isinstance(item, _TreeEdgeItem)]
        last_hovered = max(order.index(item) for item in hovered)
        first_other = min(order.index(item) for item in others)
        self.assertLess(last_hovered, first_other)

    def test_dot_and_chip_name_the_sample_and_do_not_scale_with_zoom(self) -> None:
        view = self._make_view()
        view.highlight_sample("A")
        overlay = view.hover_overlay()
        self.assertEqual([leaf for leaf, _r in overlay.dots], ["A"])
        self.assertEqual(len(overlay.chips), 1)
        name, detail, chip = overlay.chips[0]
        self.assertEqual((name, detail), ("A", "Clade 1 \u00b7 Cluster 3"))
        self.assertIsNone(overlay.more)
        dot = overlay.dots[0][1]
        self.assertAlmostEqual(dot.width(), 7.0)
        viewport = QtCore.QRectF(view.viewport().rect())
        self.assertTrue(viewport.contains(chip))
        view.zoom_by(1.8, dot.center())
        self.app.processEvents()
        zoomed = view.hover_overlay()
        self.assertAlmostEqual(zoomed.dots[0][1].width(), dot.width())
        self.assertAlmostEqual(zoomed.chips[0][2].width(), chip.width())
        self.assertAlmostEqual(zoomed.chips[0][2].height(), chip.height())
        # The overlay is painted: the viewport grab shows the amber dot.
        image = view.viewport().grab().toImage()
        center = zoomed.dots[0][1].center()
        color = image.pixelColor(int(center.x()), int(center.y()))
        self.assertEqual(color.name(), COLORS["hover"])

    def test_chips_stay_inside_the_viewport_near_its_edges(self) -> None:
        view = self._make_view(size=(420, 300))
        view.zoom_by(2.5)
        self.app.processEvents()
        viewport = QtCore.QRectF(view.viewport().rect())
        for leaf in ("A", "G"):
            tip = view._leaf_pos_scene[leaf]
            scale = view.transform().m11()
            # The tip a few pixels inside the right edge (and the top or bottom one).
            view.centerOn(
                tip.x() - (viewport.width() / 2.0 - 6.0) / scale,
                tip.y() + (viewport.height() / 2.0 - 6.0) / scale * (1 if leaf == "A" else -1),
            )
            self.app.processEvents()
            view.highlight_sample(leaf)
            overlay = view.hover_overlay()
            self.assertEqual(len(overlay.chips), 1, leaf)
            chip = overlay.chips[0][2]
            dot = overlay.dots[0][1]
            self.assertTrue(viewport.contains(chip), (leaf, chip, viewport))
            if dot.center().x() + 8.0 + chip.width() > viewport.right():
                self.assertLessEqual(chip.right(), dot.center().x())  # flipped left
            view.highlight_sample(None)

    def test_the_rest_dims_while_hovering_and_is_restored_after(self) -> None:
        view = self._make_view()
        render = mock.patch.object(view, "_render_tree", wraps=view._render_tree)
        spy = render.start()
        self.addCleanup(render.stop)
        content = view._content_item
        count = len(view.scene().items())
        self.assertIsNotNone(content)
        self.assertEqual(content.opacity(), 1.0)
        view.highlight_sample("A")
        self.assertTrue(view.hover_dimmed())
        self.assertAlmostEqual(content.opacity(), _HOVER_DIM_OPACITY)
        self.assertAlmostEqual(view._leaf_label_items["B"].effectiveOpacity(), _HOVER_DIM_OPACITY)
        self.assertAlmostEqual(view._leaf_label_items["A"].effectiveOpacity(), 1.0)
        self.assertEqual(view._leaf_label_items["A"].brush().color().name(), COLORS["hover_ink"])
        hovered = set(view._edge_items_by_leaf["A"])
        for item in self._edges(view):
            expected = 1.0 if item in hovered else _HOVER_DIM_OPACITY
            self.assertAlmostEqual(item.effectiveOpacity(), expected)
        # Moving to another leaf and clearing: no scene rebuild, same items.
        view.highlight_sample("G")
        view.highlight_samples(["C", "D"])
        view.highlight_sample(None)
        spy.assert_not_called()
        self.assertIs(view._content_item, content)
        self.assertEqual(len(view.scene().items()), count)
        # Cleared: nothing dimmed, base pens, no overlay.
        self.assertFalse(view.hover_dimmed())
        self.assertEqual(content.opacity(), 1.0)
        for item in self._edges(view):
            self.assertAlmostEqual(item.effectiveOpacity(), 1.0)
            self.assertFalse(item.pen().isCosmetic())
            self.assertEqual(item.zValue(), 0.0)
        for label in view._leaf_label_items.values():
            self.assertAlmostEqual(label.effectiveOpacity(), 1.0)
        overlay = view.hover_overlay()
        self.assertEqual((overlay.dots, overlay.chips, overlay.more), ([], [], None))

    def test_a_multi_sample_hover_caps_the_chips(self) -> None:
        view = self._make_view(size=(700, 600))
        view.highlight_samples(["A", "B", "C", "D", "E"])
        overlay = view.hover_overlay()
        self.assertEqual(len(overlay.dots), 5)
        self.assertEqual([name for name, _d, _r in overlay.chips], ["A", "B", "C"])
        self.assertIsNotNone(overlay.more)
        self.assertEqual(overlay.more[0], "+2 more")
        rects = [rect for _n, _d, rect in overlay.chips] + [overlay.more[1]]
        viewport = QtCore.QRectF(view.viewport().rect())
        for prev, rect in zip(rects, rects[1:]):
            self.assertGreaterEqual(rect.top(), prev.bottom())  # stacked, no overlap
        self.assertTrue(all(viewport.contains(rect) for rect in rects))
        # Every hovered path is thick and above, the rest dimmed.
        for leaf in ("A", "B", "C", "D", "E"):
            self.assertTrue(all(item.pen().isCosmetic() for item in view._edge_items_by_leaf[leaf]))
        self.assertAlmostEqual(view._leaf_label_items["G"].effectiveOpacity(), _HOVER_DIM_OPACITY)

    def test_the_emphasis_survives_a_redraw_and_is_left_out_of_an_export(self) -> None:
        view = self._make_view()
        view.highlight_sample("A")
        view.set_layout_mode("circular")
        self.app.processEvents()
        self.assertTrue(view.hover_dimmed())
        self.assertAlmostEqual(view._content_item.opacity(), _HOVER_DIM_OPACITY)
        self.assertEqual(len(view.hover_overlay().chips), 1)
        seen = []
        real = view._export_image

        def probe(*args, **kwargs):
            seen.append((view.hover_dimmed(), view._content_item.opacity()))
            return real(*args, **kwargs)

        with tempfile.TemporaryDirectory() as tmp, mock.patch.object(view, "_export_image", probe):
            view.export_figure(Path(tmp) / "tree.png", image_scale=1.0)
        self.assertEqual(seen, [(False, 1.0)])
        self.assertTrue(view.hover_dimmed())  # restored after the export

    def test_the_chip_never_reads_a_display_group_as_a_cluster(self) -> None:
        # Clade branch coloring sets the display groups to clades; with the
        # cluster map empty the chip must not call that clade a cluster.
        view = RaxmlTreeView(layout_mode="rectangular", use_group_backgrounds=False)
        self._views.append(view)
        view.resize(600, 400)
        view.set_leaf_groups(GROUPS, group_colors=GROUP_COLORS)
        view.set_leaf_display_groups(GROUPS, group_colors=GROUP_COLORS)
        view.set_tree_text(TREE_TEXT)
        view.show()
        self.app.processEvents()
        view.highlight_sample("A")
        _name, detail, _rect = view.hover_overlay().chips[0]
        self.assertEqual(detail, "Clade 1")

    def test_a_hovered_clade_or_cluster_label_uses_the_hover_ink(self) -> None:
        view = self._make_view()
        labels = [item for item in view.scene().items() if isinstance(item, _ActionTextItem)]
        self.assertTrue(labels)
        label = labels[0]
        label.hoverEnterEvent(None)
        self.assertEqual(label.brush().color().name(), COLORS["hover_ink"])
        label.hoverLeaveEvent(None)
        self.assertNotEqual(label.brush().color().name(), COLORS["hover_ink"])

    def test_a_short_viewport_keeps_only_the_chips_that_fit(self) -> None:
        leaves = ["A", "B", "C", "D", "E"]
        for height, max_rows in ((150, 4), (90, 3), (70, 2)):
            view = self._make_view(size=(300, height))
            view.highlight_samples(leaves)
            overlay = view.hover_overlay()
            viewport = QtCore.QRectF(view.viewport().rect())
            rects = [rect for _n, _d, rect in overlay.chips] + ([overlay.more[1]] if overlay.more else [])
            self.assertTrue(rects, height)
            self.assertLessEqual(len(rects), max_rows, height)
            for rect in rects:
                self.assertTrue(viewport.contains(rect), (height, rect, viewport))
            for prev, rect in zip(rects, rects[1:]):
                self.assertGreaterEqual(rect.top(), prev.bottom())
            # Every hovered sample is still counted: chips, "+N more" for
            # the tips on screen without a chip, and the off-screen markers.
            shown = len(overlay.chips)
            on_screen = len(overlay.dots)
            if on_screen > shown:
                self.assertIsNotNone(overlay.more)
                expected = f"+{on_screen - shown} more" if shown else f"{on_screen} samples"
                self.assertEqual(overlay.more[0], expected)
            else:
                self.assertIsNone(overlay.more)
            off = sum(int(text.split()[0]) for _s, text, _r in overlay.edges)
            self.assertEqual(on_screen + off, 5, height)
            for _s, _t, rect in overlay.edges:
                self.assertTrue(viewport.contains(rect), (height, rect, viewport))
                for chip in rects:
                    self.assertFalse(chip.intersects(rect), (height, chip, rect))

    def test_hovered_tips_out_of_view_get_an_edge_marker(self) -> None:
        # The fit zoom of a large tree, or a neighbourhood window whose rows
        # belong to tips scrolled away: the dimmed tree must still say where
        # the hovered samples are (final check of round 3).
        view = self._make_view(size=(500, 300))
        view.zoom_by(3.0)
        self.app.processEvents()
        view.centerOn(view._leaf_pos_scene["A"])
        self.app.processEvents()
        view.highlight_sample("G")
        self.assertTrue(view.hover_dimmed())
        overlay = view.hover_overlay()
        self.assertEqual((overlay.dots, overlay.chips, overlay.more), ([], [], None))
        self.assertEqual([(side, text) for side, text, _r in overlay.edges], [("below", "1 hovered sample below")])
        viewport = QtCore.QRectF(view.viewport().rect())
        marker = overlay.edges[0][2]
        self.assertTrue(viewport.contains(marker))
        self.assertGreater(marker.center().y(), viewport.center().y())
        # It is painted (amber border over the dimmed tree).
        image = view.viewport().grab().toImage()
        border = image.pixelColor(int(marker.center().x()), int(marker.top() + 0.5))
        self.assertGreater(border.red(), 200)
        self.assertLess(border.blue(), 120)
        # Samples on both sides and one on screen: a dot and chip for that
        # one, a marker above and below, and no "+N more" double count.
        view.centerOn(view._leaf_pos_scene["D"])
        self.app.processEvents()
        view.highlight_samples(["A", "D", "G"])
        overlay = view.hover_overlay()
        self.assertEqual([leaf for leaf, _r in overlay.dots], ["D"])
        self.assertEqual([name for name, _d, _r in overlay.chips], ["D"])
        self.assertIsNone(overlay.more)
        sides = {side: text for side, text, _r in overlay.edges}
        self.assertEqual(sides, {"above": "1 hovered sample above", "below": "1 hovered sample below"})
        chip = overlay.chips[0][2]
        for _s, _t, rect in overlay.edges:
            self.assertFalse(chip.intersects(rect))
        # Scrolling the tip into view replaces the marker with the dot.
        view.centerOn(view._leaf_pos_scene["G"])
        self.app.processEvents()
        view.highlight_sample("G")
        overlay = view.hover_overlay()
        self.assertEqual([leaf for leaf, _r in overlay.dots], ["G"])
        self.assertEqual(overlay.edges, [])
        view.highlight_sample(None)
        self.assertEqual(view.hover_overlay().edges, [])

    def test_a_viewport_too_short_for_one_chip_keeps_the_dots(self) -> None:
        view = self._make_view(size=(300, 26))
        view.highlight_samples(["A", "B"])
        self.assertLess(view.viewport().height(), 28)
        view.centerOn(view._leaf_pos_scene["A"])
        self.app.processEvents()
        overlay = view.hover_overlay()
        self.assertTrue(overlay.dots)
        self.assertEqual((overlay.chips, overlay.more), ([], None))

    def test_a_chip_moved_off_its_tip_gets_a_leader_line(self) -> None:
        view = self._make_view(size=(500, 300))
        view.zoom_by(0.25)  # rows closer than a chip is tall, as at the fit zoom
        self.app.processEvents()
        view.highlight_samples(["A", "B"])
        overlay = view.hover_overlay()
        tips = {leaf: rect.center() for leaf, rect in overlay.dots}
        self.assertEqual(len(overlay.chips), 2)
        for (name, _d, rect) in overlay.chips:
            off = abs(rect.center().y() - tips[name].y()) > 2.0
            led = [line for line in overlay.leaders if line[0] == tips[name]]
            self.assertEqual(bool(led), off, name)
            if led:
                end = led[0][1]
                self.assertAlmostEqual(end.y(), rect.center().y())
                self.assertIn(round(end.x(), 3), (round(rect.left(), 3), round(rect.right(), 3)))
        self.assertTrue(overlay.leaders)
        view.highlight_sample("A")
        self.assertEqual(view.hover_overlay().leaders, [])


class FitAllRowsTests(unittest.TestCase):
    """The rectangular fit puts every row in the viewport's height
    (feedback 2026-10-08: a linked heatmap shows the tree's rows, and a
    width only fit of a 153 leaf tree showed 24 of them); zoomed in, the
    row range note says which rows are in view."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def setUp(self) -> None:
        self._views: List[RaxmlTreeView] = []

    def tearDown(self) -> None:
        for view in self._views:
            view.close()
            view.deleteLater()
        self.app.processEvents()

    def _make_view(self, leaves: int = 153, size=(520, 300), layout_mode: str = "rectangular") -> RaxmlTreeView:
        view = RaxmlTreeView(layout_mode=layout_mode, use_group_backgrounds=False)
        self._views.append(view)
        view.resize(*size)
        text = "L001:0.1"
        for i in range(2, leaves + 1):
            text = f"({text},L{i:03d}:0.{i % 9 + 1}):0.01"
        view.set_tree_text(text + ";")
        view.set_row_note_enabled(True)
        view.show()
        self.app.processEvents()
        return view

    def _assert_all_rows_in_view(self, view: RaxmlTreeView) -> None:
        order = view.leaf_order()
        self.assertEqual(view.visible_row_range(), (0, len(order) - 1))
        height = view.viewport().height()
        for name in (order[0], order[-1]):
            y = view.leaf_viewport_y(name)
            self.assertTrue(0.0 <= y <= height, (name, y, height))
        self.assertEqual(view.row_range_text(), "")
        self.assertFalse(view.row_note().isVisible())

    def test_row_range_helpers(self) -> None:
        self.assertEqual(visible_rows(0.0, 10.0, 10), (0, 9))
        self.assertEqual(visible_rows(-3.0, 20.0, 10), (0, 9))
        self.assertEqual(visible_rows(4.2, 3.0, 10), (4, 6))
        self.assertIsNone(visible_rows(20.0, 3.0, 10))
        self.assertEqual(row_range_text((0, 9), 10), "")
        self.assertEqual(row_range_text((0, 59), 153), "Rows 1 to 60 of 153")
        self.assertEqual(row_range_text((4, 4), 153), "Row 5 of 153")
        self.assertEqual(row_range_text(None, 153), "No rows in view of 153")
        self.assertEqual(row_range_text(None, 0), "")

    def test_fit_shows_every_row_at_odd_heights(self) -> None:
        view = self._make_view()
        for height in (97, 233, 301, 517, 811):
            with self.subTest(height=height):
                view.resize(520, height)
                self.app.processEvents()
                self._assert_all_rows_in_view(view)

    def test_a_small_tree_still_fits_its_width(self) -> None:
        view = self._make_view(leaves=6, size=(520, 700))
        self._assert_all_rows_in_view(view)
        order = view.leaf_order()
        pitch = (view.leaf_viewport_y(order[-1]) - view.leaf_viewport_y(order[0])) / (len(order) - 1)
        # Width bound: the rows keep a readable pitch, not stretched to 700 px.
        self.assertLess(pitch * len(order), 600.0)

    def test_zoom_then_fit_returns_to_every_row_and_the_note_names_the_rows(self) -> None:
        view = self._make_view()
        self._assert_all_rows_in_view(view)
        view.zoom_by(4.0)
        self.app.processEvents()
        first, last = view.visible_row_range()
        self.assertTrue(first > 0 or last < 152)
        text = f"Rows {first + 1} to {last + 1} of 153"
        self.assertEqual(view.row_range_text(), text)
        note = view.row_note()
        self.assertTrue(note.isVisible())
        self.assertEqual(note.text(), text)
        # Bottom left of the viewport, not part of the scene (no export).
        viewport_rect = QtCore.QRect(view.viewport().mapTo(view, QtCore.QPoint(0, 0)), view.viewport().size())
        self.assertTrue(viewport_rect.contains(note.geometry()), (viewport_rect, note.geometry()))
        self.assertLess(note.geometry().left() - viewport_rect.left(), 20)
        self.assertLess(viewport_rect.bottom() - note.geometry().bottom(), 20)
        self.assertNotIn(note, view.scene().items())
        bar = view.verticalScrollBar()
        bar.setValue(bar.maximum())
        self.app.processEvents()
        self.assertEqual(view.visible_row_range()[1], 152)
        self.assertTrue(view.row_range_text().endswith("to 153 of 153"))
        view.fit_to_window()
        self.app.processEvents()
        self._assert_all_rows_in_view(view)

    def test_a_fit_below_the_zoom_limit_still_zooms_in(self) -> None:
        view = self._make_view(leaves=400, size=(520, 120))
        self._assert_all_rows_in_view(view)
        scale = view.transform().m22()
        self.assertLess(scale, 0.03)
        view.zoom_by(1.15)
        self.assertAlmostEqual(view.transform().m22(), scale * 1.15, places=6)
        for _ in range(12):
            view.zoom_by(1.0 / 1.15)  # no further out than about half the fit
        self.assertLess(view.transform().m22(), scale * 0.7)
        self.assertGreater(view.transform().m22(), scale * 0.4)

    def test_a_tall_tree_keeps_its_width_at_fit(self) -> None:
        # Review 2026-10-08: a uniform scale left a 153 leaf tree in a short
        # panel a sliver about 10 px wide. The rows fit the height, the
        # tree still fits the width.
        view = self._make_view(size=(520, 150))
        self._assert_all_rows_in_view(view)
        transform = view.transform()
        self.assertGreater(transform.m11(), 3.0 * transform.m22())
        content = view.mapFromScene(view._fit_content_rect()).boundingRect()
        self.assertGreater(content.width(), 0.6 * view.viewport().width())
        self.assertLessEqual(content.right(), view.viewport().width() + 1)
        # The branches are drawn visibly on screen at that scale.
        image = view.viewport().grab().toImage()
        ink_columns = 0
        for x in range(0, image.width(), 4):
            if any(QtGui.QColor(image.pixel(x, y)).lightness() < 200 for y in range(0, image.height(), 2)):
                ink_columns += 1
        self.assertGreater(ink_columns * 4, 0.4 * image.width())
        # Alt+wheel spreads the rows at the same width.
        width_scale = transform.m11()
        view.stretch_by(1.0, 2.0)
        self.assertAlmostEqual(view.transform().m11(), width_scale, places=6)
        view.fit_to_window()
        self.app.processEvents()
        self._assert_all_rows_in_view(view)
        self.assertAlmostEqual(view.transform().m11(), width_scale, places=6)

    def test_overlapping_cluster_chips_are_hidden_until_there_is_room(self) -> None:
        # Review 2026-10-08: at the fit of a tall tree in a short panel the
        # screen stable cluster chips of small runs piled on each other.
        view = RaxmlTreeView(layout_mode="rectangular", use_group_backgrounds=False, show_cluster_annotations=True)
        self._views.append(view)
        view.resize(520, 150)
        text = "L001:0.1"
        for i in range(2, 154):
            text = f"({text},L{i:03d}:0.{i % 9 + 1}):0.01"
        clusters = {f"L{i:03d}": ("1" if i <= 4 else "2" if i <= 10 else "3" if i <= 70 else "4") for i in range(1, 154)}
        view.set_leaf_clusters(clusters)
        view.set_tree_text(text + ";")
        view.show()
        self.app.processEvents()
        chips = [item for item in view.scene().items() if item.data(_SCENE_DATA_KIND) == "cluster-label"]
        self.assertEqual(len(chips), 4)

        def screen_rect(item: QtWidgets.QGraphicsItem) -> QtCore.QRectF:
            anchor = view.viewportTransform().map(item.scenePos())
            return item.boundingRect().united(item.childrenBoundingRect()).translated(anchor)

        shown = [item for item in chips if item.isVisible()]
        self.assertLess(len(shown), len(chips))
        for i, first in enumerate(shown):
            for second in shown[i + 1:]:
                self.assertFalse(screen_rect(first).intersects(screen_rect(second)))
        biggest = max(chips, key=lambda item: int(item.data(5)))
        self.assertEqual(biggest.data(_SCENE_DATA_NAME), "4")
        self.assertTrue(biggest.isVisible())
        # Exports show every chip, and the screen keeps its choice.
        with view._scene_scaled_annotations_for_export():
            self.assertTrue(all(item.isVisible() for item in chips))
        self.assertEqual([item for item in chips if item.isVisible()], shown)
        # Zoomed in far enough, every chip has room again.
        for _ in range(30):
            view.zoom_by(1.3)
        self.assertTrue(all(item.isVisible() for item in chips))

    def test_zoom_out_after_the_panel_shrinks_reaches_every_row(self) -> None:
        # Review 2026-10-08: the zoom out limit was the last fit's scale, so
        # after a zoom in and a shorter panel, zooming out stopped short of
        # the rows (Rows 1 to 207 of 400) until Fit.
        view = self._make_view(leaves=400, size=(520, 400))
        self._assert_all_rows_in_view(view)
        for _ in range(3):
            view.zoom_by(1.15)
        view.resize(520, 150)
        self.app.processEvents()
        self.assertNotEqual(view.row_range_text(), "")
        for _ in range(30):
            view.zoom_by(1.0 / 1.15)
        self.app.processEvents()
        self._assert_all_rows_in_view(view)
        view.fit_to_window()
        self.app.processEvents()
        self._assert_all_rows_in_view(view)

    def test_wheel_navigation_camera_zoom_scroll_and_fit(self) -> None:
        # Feedback 2026-10-15: Ctrl+wheel zooms the whole canvas like an
        # image (both axes and the text, under the pointer), Shift+wheel
        # scrolls sideways, Fit returns to the fitted mode.
        view = self._make_view()
        viewport = view.viewport()
        modifiers = QtCore.Qt.KeyboardModifier

        def wheel(pos: QtCore.QPointF, delta: int, mods) -> None:
            event = QtGui.QWheelEvent(
                pos, QtCore.QPointF(viewport.mapToGlobal(pos.toPoint())), QtCore.QPoint(), QtCore.QPoint(0, delta),
                QtCore.Qt.MouseButton.NoButton, mods, QtCore.Qt.ScrollPhase.NoScrollPhase, False,
            )
            QtWidgets.QApplication.sendEvent(viewport, event)

        def label_height() -> float:
            label = view._leaf_label_items[view.leaf_order()[0]]
            return label.deviceTransform(view.viewportTransform()).mapRect(label.boundingRect()).height()

        fit = (view.transform().m11(), view.transform().m22(), label_height())
        self.assertTrue(view.is_fitted())
        anchor = QtCore.QPointF(viewport.width() * 0.6, viewport.height() * 0.4)
        under = view._scene_point_at(anchor)
        for _ in range(3):
            wheel(anchor, 120, modifiers.ControlModifier)
        self.assertFalse(view.is_fitted())
        self.assertAlmostEqual(view.camera_zoom(), 1.2 ** 3, places=6)
        self.assertAlmostEqual(view.transform().m11(), fit[0] * 1.2 ** 3, places=6)
        self.assertAlmostEqual(view.transform().m22(), fit[1] * 1.2 ** 3, places=6)
        self.assertAlmostEqual(label_height(), fit[2] * 1.2 ** 3, delta=0.01)
        moved = view.viewportTransform().map(under) - anchor
        self.assertLess(max(abs(moved.x()), abs(moved.y())), 1.0)
        bar_x, bar_y = view.horizontalScrollBar(), view.verticalScrollBar()
        before = (bar_x.value(), bar_y.value())
        wheel(anchor, -120, modifiers.ShiftModifier)
        self.assertGreater(bar_x.value(), before[0])
        self.assertEqual(bar_y.value(), before[1])
        view.fit_to_window()
        self.assertTrue(view.is_fitted())
        self.assertEqual(view.camera_zoom(), 1.0)
        self.assertAlmostEqual(view.transform().m11(), fit[0], places=6)
        self.assertAlmostEqual(view.transform().m22(), fit[1], places=6)
        self._assert_all_rows_in_view(view)

    def test_circular_layout_has_no_row_note(self) -> None:
        view = self._make_view(layout_mode="circular")
        view.zoom_by(4.0)
        self.app.processEvents()
        self.assertIsNone(view.visible_row_range())
        self.assertEqual(view.row_range_text(), "")
        self.assertFalse(view.row_note().isVisible())


if __name__ == "__main__":
    unittest.main()
