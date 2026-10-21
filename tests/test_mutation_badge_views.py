"""A gene with a known resistance mutation (tag X) is drawn in the AMR red
with the mutation badge in both gene views: the neighborhood plot and the
genome track."""
from __future__ import annotations

import os
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

from PyQt6 import QtGui, QtWidgets

from annotation_explorer_polars.context_canvas import ContextCanvas
from annotation_explorer_polars.context_model import CATEGORY_COLORS, gene_type_fill
from annotation_explorer_polars.igv_view import GeneSegment, IgvTrackView, _MutationBadgeItem

PLAIN = {"Gene": "fdhD", "Tags": "", "Start": 100, "Stop": 900, "Sample": "S1", "Contig": "2"}
GYRA = {
    "Gene": "gyrA", "Tags": "X", "Start": 1000, "Stop": 3700, "Sample": "S1", "Contig": "2",
    "ConsensusCategory": "OTHER",
    "ResistanceMutations": "gyrA_S81L", "ResistanceMutationClasses": "QUINOLONE",
}
OXA = {"Gene": "blaOXA-66", "Tags": "AX", "Start": 4000, "Stop": 4800, "Sample": "S1", "Contig": "2",
       "ResistanceMutations": "blaOXA-66_W222L"}


class MutationBadgeViewsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def test_neighborhood_gene_is_red_with_badge(self) -> None:
        canvas = ContextCanvas()
        canvas.resize(1200, 500)
        genes = [
            {"start": row["Start"], "end": row["Stop"], "strand": "+", "label": row["Gene"],
             "row_data": dict(row), "is_center": row is GYRA}
            for row in (PLAIN, GYRA, OXA)
        ]
        canvas.set_rows([{"sample": "S1", "contig": "2", "genes": genes, "is_source": True}])
        layout_genes = {gene.gene.label: gene for gene in canvas.layout_result.rows[0].genes}
        amr = QtGui.QColor(CATEGORY_COLORS["AMR"]).name()
        self.assertEqual(canvas.gene_fill(layout_genes["gyrA"]).name(), amr)
        self.assertIn("Usually not a resistance gene; carries gyrA_S81L (quinolone)", layout_genes["gyrA"].gene.mutation_note)
        self.assertTrue(layout_genes["blaOXA-66"].gene.mutation_note.startswith("Resistance gene that also carries"))
        self.assertEqual(layout_genes["fdhD"].gene.mutation_note, "")
        legend = [kind for group in canvas._legend_groups("1 kb") for kind, _text, _color in group]
        self.assertIn("mutation", legend)

    def test_genome_track_gene_is_red_with_badge(self) -> None:
        view = IgvTrackView()
        view.resize(900, 240)
        segments = [
            GeneSegment(row["Start"], row["Stop"], "+", row["Gene"], QtGui.QColor(gene_type_fill(row)), 1, dict(row))
            for row in (PLAIN, GYRA)
        ]
        view.set_track(segments, 5000, bp_per_px=5.0, hit_bp_positions=[2300])
        view.set_selected_row(GYRA)
        badges = [item for item in view.scene().items() if isinstance(item, _MutationBadgeItem)]
        self.assertEqual(len(badges), 1)
        self.assertIn("gyrA_S81L", badges[0].toolTip())
        fills = {gene.label: item.brush().color().name() for item, gene in view._gene_items}
        self.assertEqual(fills["gyrA"], QtGui.QColor(CATEGORY_COLORS["AMR"]).name())
        self.assertEqual(fills["fdhD"], QtGui.QColor(CATEGORY_COLORS["OTHER"]).name())
        self.assertTrue(view._selection_outline.isVisible())


if __name__ == "__main__":
    unittest.main()
