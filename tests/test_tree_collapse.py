"""Collapsing a Phylogeny subtree into a triangle (tree_view, phylo)."""
from __future__ import annotations

import os
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

from PyQt6 import QtWidgets

from annotation_explorer_polars.phylo import ORDER_LADDER_DEPTH, NewickParser, group_starts, leaf_names
from annotation_explorer_polars.tree_view import RaxmlTreeView

# Cluster 3 is not monophyletic: cluster 1 sits between its two blocks.
TREE_TEXT = "(((A:0.1,B:0.2):0.3,(C:0.1,D:0.1):0.2):0.1,((E:1,F:0.1):0.1,G:0.1):0.1);"
CLADES = {"A": "1", "B": "1", "C": "1", "D": "1", "E": "2", "F": "2", "G": "2"}
CLUSTERS = {"A": "3", "B": "3", "C": "1", "D": "1", "E": "3", "F": "3", "G": "2"}


class TreeCollapseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def _view(self) -> RaxmlTreeView:
        view = RaxmlTreeView(layout_mode="rectangular", show_cluster_annotations=True)
        self.addCleanup(view.close)
        view.resize(800, 600)
        view.set_leaf_groups(CLADES)
        view.set_leaf_clusters(CLUSTERS)
        view.set_tree_text(TREE_TEXT)
        return view

    def test_collapsing_keeps_the_leaf_order_and_never_swallows_another_group(self) -> None:
        root = NewickParser(TREE_TEXT).parse()
        starts = group_starts(root, CLUSTERS, str)
        # A group that is not monophyletic starts at each of its pure subtrees.
        self.assertEqual(
            [(group, sorted(leaf_names(node, str))) for group, node in starts],
            [("3", ["A", "B"]), ("1", ["C", "D"]), ("3", ["E", "F"])],
        )
        view = self._view()
        order = view.leaf_order()
        bottom = max(point.y() for point in view._leaf_pos_scene.values())
        handles = view.collapse_handles()
        for _name, members, _collapsed in handles:
            groups = {CLADES[m] for m in members} if len(members) > 2 else {CLUSTERS[m] for m in members}
            self.assertEqual(len(groups), 1, members)
        # The 2 genome blocks of the split cluster 3 get no handle (a right
        # click collapses them); the whole clades and cluster 1 do.
        self.assertEqual(
            sorted(members for _name, members, _collapsed in handles),
            [("A", "B", "C", "D"), ("C", "D"), ("E", "F", "G")],
        )
        for members in (("A", "B"), ("C", "D"), ("E", "F")):
            view.collapse_subtree(members)
        self.assertEqual(
            view.collapsed_rows(),
            [
                # A block of a group that is not monophyletic says how much of it it is.
                ("Cluster 3 · 2 of 4 genomes", ("A", "B")),
                ("Cluster 1 · 2 genomes", ("C", "D")),
                ("Cluster 3 · 2 of 4 genomes", ("E", "F")),
            ],
        )
        self.assertEqual(view.leaf_order(), order)
        # A whole clade over two collapsed clusters: one triangle, the
        # clusters inside stay collapsed. Every leaf is still listed, in its
        # order (a linked heatmap keeps its rows), while the tree is shorter.
        view.collapse_subtree(["A", "B", "C", "D"])
        view.expand_subtree(["E", "F"])
        self.assertEqual([label for label, _m in view.collapsed_rows()], ["Clade 1 · 4 genomes"])
        self.assertEqual(view.leaf_order(), order)
        self.assertLess(max(point.y() for point in view._leaf_pos_scene.values()), bottom)
        # Another order keeps the collapsed subtrees.
        view.set_order_mode(ORDER_LADDER_DEPTH)
        self.assertEqual(len(view.collapsed_rows()), 1)
        self.assertEqual(len(view.collapsed_groups()), 3)
        view.expand_all()
        self.assertEqual(view.collapsed_rows(), [])
        # A clade the root sits in (a grade at the root, as the Groups root
        # can put it) gets a handle on the root and collapses all the same.
        view.set_tree_text("(A:1,(B:1,(C:1,(D:1,E:1):1):1):1);")
        order = view.leaf_order()
        grade = next(members for name, members, _c in view.collapse_handles() if name == "Clade 1")
        self.assertEqual(sorted(grade), ["A", "B", "C", "D"])
        view.collapse_subtree(grade)
        self.assertEqual([label for label, _m in view.collapsed_rows()], ["Clade 1 · 4 genomes"])
        self.assertEqual(view.leaf_order(), order)
        self.assertTrue(view.is_collapsed(grade))

    def test_collapsed_state_round_trips_and_drops_what_the_tree_lacks(self) -> None:
        view = self._view()
        view.collapse_subtree(["E", "F"])
        view.collapse_subtree(["A", "B", "C", "D"])
        saved = [list(group) for group in view.collapsed_groups()]
        other = self._view()
        other.set_collapsed_groups(saved)
        self.assertEqual(other.collapsed_groups(), [("E", "F"), ("A", "B", "C", "D")])
        self.assertEqual([members for _label, members in other.collapsed_rows()], [("A", "B", "C", "D"), ("E", "F")])
        # Another tree without the subtree E, F leaves it out; A to D stays.
        # Back on the first tree, E, F is collapsed again.
        other.set_tree_text("(((A:0.1,B:0.2):0.3,(C:0.1,D:0.1):0.2):0.1,((E:1,G:0.1):0.1,F:0.1):0.1);")
        self.assertEqual(other.collapsed_groups(), [("A", "B", "C", "D")])
        other.set_tree_text(TREE_TEXT)
        self.assertEqual(len(other.collapsed_rows()), 2)
        # A tree pruned to the sample set: its triangle expands and
        # collapses the whole subtree, and keeps it for another filter.
        other.expand_all()
        other.set_prune_inactive(True)
        other.set_active_samples(["A", "B", "C", "E", "F", "G"])
        other.collapse_subtree(["A", "B", "C"])
        self.assertEqual(other.collapsed_groups(), [("A", "B", "C", "D")])
        self.assertEqual(other.collapsed_rows(), [("Clade 1 · 3 genomes", ("A", "B", "C"))])
        other.set_active_samples(["A", "B", "D", "G"])
        self.assertEqual(other.collapsed_rows(), [("Clade 1 · 3 genomes", ("A", "B", "D"))])
        other.toggle_collapsed(["A", "B", "D"])
        self.assertEqual((other.collapsed_groups(), other.collapsed_rows()), ([], []))


if __name__ == "__main__":
    unittest.main()
