from __future__ import annotations

import os
import time
import unittest
from types import SimpleNamespace

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

import polars as pl
from PyQt6 import QtWidgets

from annotation_explorer_polars.model.filter_state import FilterState
from annotation_explorer_polars.tree_controller import TreeControllerMixin


_HEADERS = [
    "Name", "Card", "Origin", "n_total", "ATB", "Country", "db_total_hits",
    "AMRFinderPlus", "ST", "Quality", "Complete", "Contam.",
]

_DEFAULT_RECORDS = {
    "S1": [
        {"Sample": "S1", "Contig": "1", "n_hits": 2, "n_total": 6},
        {"Sample": "S1", "Contig": "2", "n_hits": 0, "n_total": 4},
    ],
    "S2": [{"Sample": "S2", "Contig": "1", "n_hits": 1, "n_total": 12}],
    "S3": [{"Sample": "S3", "Contig": "1", "n_hits": 0, "n_total": 8}],
}


class _SampleTreeHarness(TreeControllerMixin, QtWidgets.QWidget):
    """Minimal host for the sample tree with the app's column layout. The
    tree code reads the run's data from ``self.run`` (RunData); this host
    gives it only the fields the tree uses."""

    def __init__(self) -> None:
        super().__init__()
        self.meta_columns_for_tree = ["Country"]
        self.run = SimpleNamespace(
            annotation_db_columns=["AMRFinderPlus"],
            sample_qc_columns=[],
            sample_total_counts={"S1": 10, "S2": 12, "S3": 8},
            sample_atb_counts={},
            contig_atb_counts={},
            contig_total_count_map={},
            mobsuite_plasmid_members={},
            plasmid_cluster_index=None,
            marker_db_category_map={},
            annotation_hit_counts_pl=pl.DataFrame(
                schema={"Sample": pl.Utf8, "Contig": pl.Utf8, "DB": pl.Utf8, "Hits": pl.Int64}
            ),
        )

        layout = QtWidgets.QVBoxLayout(self)
        self.gene_query = QtWidgets.QLineEdit()
        self.product_query = QtWidgets.QLineEdit()
        self.tree = QtWidgets.QTreeWidget()
        self.tree.setHeaderLabels(_HEADERS)
        self._col_st, self._col_quality, self._col_completeness, self._col_contamination = 8, 9, 10, 11
        self._fixed_tree_columns = (0, 2, 8, 9, 10, 11)
        self.tree.setSortingEnabled(True)
        self.tree.itemExpanded.connect(self._on_item_expanded)
        self.tree.itemCollapsed.connect(self._on_item_collapsed)
        layout.addWidget(self.gene_query)
        layout.addWidget(self.product_query)
        layout.addWidget(self.tree)

    # The window's context holds the applied filters; None reads the boxes.
    applied = None

    def _applied_filters(self) -> FilterState:
        return self.applied if self.applied is not None else FilterState.from_widgets(self)

    def _apply_tree_sample_hover(self, _sample) -> None:
        pass

    def populate(self, counts=None, contig_records=None, samples=("S1", "S2", "S3")) -> None:
        samples = list(samples)
        countries = {"S1": "Romania", "S2": "Italy", "S3": "Spain"}
        meta = pl.DataFrame(
            {"Sample": samples, "Country": [countries.get(s, "Spain") for s in samples]}
        )
        records = contig_records if contig_records is not None else _DEFAULT_RECORDS
        if counts is None:
            searching = bool(self.gene_query.text().strip() or self.product_query.text().strip())
            counts = {"S1": 2, "S2": 1} if searching else dict(self.run.sample_total_counts)
        self._populate_metadata_tree(meta, records, counts, samples)
        self._rerender_visible_table_columns()


class SampleTreeSearchExpandTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def setUp(self) -> None:
        self.host = _SampleTreeHarness()
        self.host.resize(900, 400)
        self.host.show()
        self.app.processEvents()

    def tearDown(self) -> None:
        self.host.close()
        self.app.processEvents()

    def _item(self, sample: str) -> QtWidgets.QTreeWidgetItem:
        return self.host._sample_tree_items[sample]

    def _expanded(self) -> set:
        return {s for s, item in self.host._sample_tree_items.items() if item.isExpanded()}

    def _search(self, gene: str = "", product: str = "") -> None:
        self.host.gene_query.setText(gene)
        self.host.product_query.setText(product)
        self.host.populate()

    def test_no_copies_column_during_search(self) -> None:
        self._search("blaOXA-23")
        tree = self.host.tree
        labels = [tree.headerItem().text(i) for i in range(tree.columnCount())]
        self.assertEqual(labels, _HEADERS)
        self.assertNotIn("Copies", labels)
        self.assertFalse(hasattr(self.host, "_col_copies"))

    def test_names_are_plain_during_search(self) -> None:
        self._search("blaOXA-23")
        for sample in ("S1", "S2", "S3"):
            item = self._item(sample)
            self.assertEqual(item.text(0), sample)
            self.assertNotEqual(item.foreground(0).color().name(), "#2563eb")

    def test_zero_match_sample_greyed_and_collapsed(self) -> None:
        self._search("blaOXA-23")
        s3 = self._item("S3")
        self.assertEqual(s3.foreground(0).color().name(), "#888888")
        self.assertFalse(s3.isExpanded())

    def test_matching_samples_expanded_with_children_visible(self) -> None:
        self._search("blaOXA-23")
        self.assertEqual(self._expanded(), {"S1", "S2"})
        tree = self.host.tree
        s1 = self._item("S1")
        self.assertEqual(s1.childCount(), 2)
        for i in range(s1.childCount()):
            child = s1.child(i)
            self.assertFalse(child.isHidden())
            self.assertTrue(tree.visualItemRect(child).isValid())
        self.assertEqual(s1.child(0).text(0), "Contig 1")

    def test_product_search_also_expands(self) -> None:
        self._search(product="beta-lactamase")
        self.assertEqual(self._expanded(), {"S1", "S2"})

    def test_manual_collapse_sticks_across_refresh_with_same_query(self) -> None:
        self._search("blaOXA-23")
        self._item("S1").setExpanded(False)
        self.app.processEvents()
        self.host.populate()
        self.assertEqual(self._expanded(), {"S2"})
        self.host.populate()
        self.assertEqual(self._expanded(), {"S2"})

    def test_changed_query_reapplies_expansion(self) -> None:
        self._search("blaOXA")
        self._item("S1").setExpanded(False)
        self.host.gene_query.setText("blaOXA-23")
        self.host.populate(counts={"S1": 2, "S3": 1})
        self.assertEqual(self._expanded(), {"S1", "S3"})

    def test_same_query_new_match_opens(self) -> None:
        # Another filter can change which samples match under the same text;
        # newly matching samples open, the manual collapse stays.
        self._search("blaOXA-23")
        self._item("S1").setExpanded(False)
        self.host.populate(counts={"S1": 2, "S3": 1})
        self.assertEqual(self._expanded(), {"S3"})

    def test_manual_collapse_sticks_when_sample_briefly_stops_matching(self) -> None:
        self._search("blaOXA-23")
        self._item("S1").setExpanded(False)
        self.host.populate(counts={"S2": 1})
        self.assertEqual(self._expanded(), {"S2"})
        self.host.populate(counts={"S1": 2, "S2": 1})
        self.assertEqual(self._expanded(), {"S2"})

    def test_reopened_sample_leaves_manual_collapse_set(self) -> None:
        self._search("blaOXA-23")
        self._item("S1").setExpanded(False)
        self.host.populate()
        self._item("S1").setExpanded(True)
        self.host.populate(counts={"S2": 1})
        self.host.populate()
        self.assertEqual(self._expanded(), {"S1", "S2"})

    def test_case_only_edit_keeps_manual_collapse(self) -> None:
        self._search("blaOXA-23")
        self._item("S1").setExpanded(False)
        self._search("BLAoxa-23")
        self.assertEqual(self._expanded(), {"S2"})

    def test_tree_follows_the_applied_search_not_the_box(self) -> None:
        # The tree's counts belong to the applied search; the box can be ahead
        # of them while the user types (it is the context's draft): the old
        # search's state (manual collapse, tooltip text) is kept until the
        # new search is applied.
        self._search("blaOXA-23")
        self._item("S1").setExpanded(False)
        self.host.applied = FilterState(gene="blaOXA-23")
        self.host.gene_query.setText("blaOXA-58")
        self.host.populate(counts={"S1": 2, "S2": 1})
        self.assertEqual(self._expanded(), {"S2"})
        self.assertIn("1 copy of blaOXA-23 in this sample", self._item("S2").toolTip(0))
        self.host.applied = FilterState(gene="blaOXA-58")
        self.host.populate(counts={"S1": 2, "S3": 1})
        self.assertEqual(self._expanded(), {"S1", "S3"})
        self.assertIn("2 copies of blaOXA-58 in this sample", self._item("S1").toolTip(0))

    def test_clearing_search_restores_pre_search_expansion(self) -> None:
        self.host.populate()
        self._item("S3").setExpanded(True)
        self.assertEqual(self._expanded(), {"S3"})
        self._search("blaOXA-23")
        self.assertEqual(self._expanded(), {"S1", "S2"})
        self._search("blaOXA-2")
        self._search("")
        self.assertEqual(self._expanded(), {"S3"})
        # Ordinary refreshes after the search keep the user's state.
        self.host.populate()
        self.assertEqual(self._expanded(), {"S3"})

    def test_tooltip_lines_kept(self) -> None:
        self._search("blaOXA-23")
        self.assertIn("2 copies of blaOXA-23 in this sample", self._item("S1").toolTip(0))
        self.assertIn("1 copy of blaOXA-23 in this sample", self._item("S2").toolTip(0))
        s1 = self._item("S1")
        self.assertIn("2 copies of blaOXA-23 in this contig", s1.child(0).toolTip(0))
        self.assertIn("0 copies of blaOXA-23 in this contig", s1.child(1).toolTip(0))
        self._search("")
        self.assertNotIn("copies of", self._item("S1").toolTip(0))

    def test_plasmid_group_tooltip_sums_members(self) -> None:
        self.host.run.mobsuite_plasmid_members = {("S1", "AA001"): ["2", "3"]}
        plasmid = {"MobSuiteIsPlasmid": "true", "MobSuitePlasmidID": "AA001"}
        records = {
            "S1": [
                {"Sample": "S1", "Contig": "1", "n_hits": 1, "n_total": 6},
                {"Sample": "S1", "Contig": "2", "n_hits": 2, "n_total": 4, **plasmid},
                {"Sample": "S1", "Contig": "3", "n_hits": 1, "n_total": 3, **plasmid},
            ],
        }
        self.host.gene_query.setText("blaOXA")
        self.host.populate(counts={"S1": 4}, contig_records=records, samples=["S1"])
        s1 = self._item("S1")
        self.assertTrue(s1.isExpanded())
        children = {s1.child(i).text(0): s1.child(i) for i in range(s1.childCount())}
        self.assertIn(
            "3 copies of blaOXA in this plasmid group", children["Plasmid AA001"].toolTip(0)
        )

    def test_many_samples_stay_fast(self) -> None:
        samples = [f"X{i:03d}" for i in range(150)]
        records = {
            s: [{"Sample": s, "Contig": str(c), "n_hits": 1 if c == 1 else 0, "n_total": 5}
                for c in range(1, 6)]
            for s in samples
        }
        self.host.run.sample_total_counts = {s: 25 for s in samples}
        counts = {s: 1 for s in samples[::2]}
        expanded_signals = []
        self.host.tree.itemExpanded.connect(lambda item: expanded_signals.append(item))
        self.host.gene_query.setText("blaOXA-23")
        start = time.perf_counter()
        self.host.populate(counts=counts, contig_records=records, samples=samples)
        self.app.processEvents()
        elapsed = time.perf_counter() - start
        self.assertEqual(self._expanded(), set(samples[::2]))
        # Expansion is restored with signals blocked: no per-item signal storm.
        self.assertEqual(expanded_signals, [])
        self.assertLess(elapsed, 3.0)


if __name__ == "__main__":
    unittest.main()
