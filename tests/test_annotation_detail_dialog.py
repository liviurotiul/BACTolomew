from __future__ import annotations

import json
import os
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

import pandas as pd
from PyQt6 import QtWidgets

from annotation_explorer_polars.annotation_detail_dialog import AnnotationDetailDialog


class AnnotationDetailDialogTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def setUp(self) -> None:
        self.nt_sequence = "ATG" + ("ACGT" * 85) + "TAA"
        self.aa_sequence = "MKKTLALAVAVAAC" + ("Q" * 260)
        self.row = {
            "Sample": "S1",
            "Contig": "NODE_1",
            "Gene": "blaOXA-23",
            "GeneID": "family-23",
            "Product": "OXA-23 family carbapenem-hydrolyzing class D beta-lactamase",
            "AnnotationID": "BKT0001",
            "Start": 100,
            "Stop": 399,
            "StartInt": 100,
            "StopInt": 399,
            "Strand": "-",
            "NCBI_NRP_ID": "WP_000941335.1",
            "Refseq_ID": "YP_001713491.1",
            "Nucleotide_Sequence": self.nt_sequence,
            "Amino_Acid_Sequence": self.aa_sequence,
            "Predictors": "Bakta, AMRFinderPlus",
            "Tags": "AM",
            "custom_field": ["preserved", 2],
            "ResistanceMechanisms": "hidden internal value",
        }
        marker_hits = pd.DataFrame(
            [
                {
                    "AnnotationID": "BKT0001",
                    "DB": "AMRFinderPlus",
                    "Method": "protein",
                    "Subject_coverage": 98.4,
                    "Identity": 97.8,
                    "Product": "AMR beta-lactamase",
                    "Bitscore": 300,
                },
                {
                    "AnnotationID": "BKT0001",
                    "DB": "CARD",
                    "Method": "diamond",
                    "Subject_coverage": 100.0,
                    "Identity": 99.2,
                    "Product": "OXA-23 family carbapenemase",
                    "Bitscore": 500,
                },
                {
                    "AnnotationID": "BKT9999",
                    "DB": "Unrelated",
                    "Method": "diamond",
                    "Product": "must not appear",
                    "Bitscore": 900,
                },
                {
                    "AnnotationID": "BKT0001",
                    "DB": "",
                    "Method": "",
                    "Product": "",
                    "Bitscore": 999,
                },
            ]
        )
        ann_df = pd.DataFrame(
            [
                self.row,
                {
                    **self.row,
                    "Sample": "S2",
                    "AnnotationID": "BKT0002",
                    "Amino_Acid_Sequence": "M" + ("A" * 50),
                },
            ]
        )
        self.dialog = AnnotationDetailDialog(
            None,
            self.row,
            list(self.row),
            ann_df=ann_df,
            meta_df=pd.DataFrame([{"Sample": "S1", "Country": "Romania"}]),
            marker_hits_df=marker_hits,
        )
        self.dialog.show()
        self.app.processEvents()

    def tearDown(self) -> None:
        self.dialog.close()
        self.app.processEvents()

    def test_modern_card_hierarchy_is_responsive_and_keeps_identity_visible(self) -> None:
        self.assertEqual(self.dialog.objectName(), "AnnotationDetailDialog")
        self.assertEqual(self.dialog._tabs.tabText(0), "Overview")
        self.assertEqual(self.dialog._tabs.tabText(1), "All fields")
        self.assertEqual(
            self.dialog.findChild(QtWidgets.QLabel, "AnnotationTitle").text(),
            "blaOXA-23",
        )
        self.assertIn(
            "carbapenem-hydrolyzing",
            self.dialog.findChild(QtWidgets.QLabel, "AnnotationProduct").text(),
        )
        self.assertIsNotNone(
            self.dialog.findChild(QtWidgets.QScrollArea, "AnnotationSummaryScroll")
        )
        self.assertFalse(self.dialog.findChildren(QtWidgets.QGroupBox))
        self.assertLessEqual(self.dialog.minimumSizeHint().width(), 760)
        self.assertLessEqual(self.dialog.width(), 940)
        self.assertEqual(self.dialog._context_button.property("role"), "primary")
        self.assertEqual(self.dialog._align_button.text(), "Align sequences…")

    def test_sequences_copy_exact_text_without_inflating_the_window(self) -> None:
        self.assertEqual(len(self.dialog._sequence_edits), 2)
        self.assertEqual(self.dialog._sequence_edits[0].text(), self.nt_sequence)
        self.assertEqual(self.dialog._sequence_edits[1].text(), self.aa_sequence)
        copy_buttons = self.dialog.findChildren(QtWidgets.QToolButton, "AnnotationCopyButton")
        self.assertEqual(len(copy_buttons), 2)
        copy_buttons[1].click()
        self.app.processEvents()
        self.assertEqual(self.app.clipboard().text(), self.aa_sequence)

    def test_marker_linking_sort_and_raw_payload_are_preserved(self) -> None:
        table = self.dialog._marker_table
        self.assertIsNotNone(table)
        self.assertEqual(table.rowCount(), 2)
        self.assertEqual(table.item(0, 0).text(), "CARD")
        self.assertEqual(table.item(1, 0).text(), "AMRFinderPlus")
        self.assertLessEqual(table.height(), 42 + (6 * 32))
        self.assertEqual(table.cellWidget(0, 2).value(), 100)
        self.assertEqual(table.cellWidget(0, 3).value(), 99)

        payload = json.loads(self.dialog._raw_editor.toPlainText())
        self.assertEqual(payload["custom_field"], ["preserved", 2])
        self.assertNotIn("ResistanceMechanisms", payload)
        self.assertEqual([hit["DB"] for hit in payload["MarkerHits"]], ["CARD", "AMRFinderPlus"])

    def test_position_and_orientation_semantics_are_unchanged(self) -> None:
        self.assertEqual(self.dialog._format_position(), "100 ← 399")
        self.assertEqual(self.dialog._format_orientation(), "Reverse (-)")
        badge_text = [
            label.text()
            for label in self.dialog.findChildren(QtWidgets.QLabel)
            if label.property("role") == "annotationChip"
        ]
        self.assertIn("AMR", badge_text)
        self.assertIn("Mobile", badge_text)


if __name__ == "__main__":
    unittest.main()
