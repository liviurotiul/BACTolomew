"""P20 concepts: AMRFinderPlus POINT calls land on the right CDS as tag X with the badge
data, acquired-gene rows are ignored, and the organism comes from the MLST species.

Qt-free; tiny fixtures only.
"""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pipeline import AmrFinderWorkflow as AF  # noqa: E402
from annotation_explorer_polars.data import _finalize_annotations_df  # noqa: E402
from annotation_explorer_polars.resistance_mutation import mutation_note  # noqa: E402

HEADER = [
    "Name", "Protein id", "Contig id", "Start", "Stop", "Strand", "Element symbol",
    "Element name", "Scope", "Type", "Subtype", "Class", "Subclass", "Method",
]


def _report_row(protein, contig, start, stop, symbol, subtype, cls, subclass, kind="AMR"):
    return ["S1", protein, contig, str(start), str(stop), "+", symbol, symbol,
            "core", kind, subtype, cls, subclass, "POINTP"]


class PointMutationsOnCdsTest(unittest.TestCase):
    def test_point_rows_tag_the_matching_cds(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            (out / "amrfinder").mkdir()
            rows = [
                HEADER,
                # matched by locus tag
                _report_row("ABC_00437", "contig_2", 43837, 46551, "gyrA_S81L", "POINT", "QUINOLONE", "QUINOLONE"),
                # no protein id: matched by contig and coordinates
                _report_row("NA", "contig_32", 10756, 12975, "parC_S84L", "POINT", "QUINOLONE", "QUINOLONE"),
                # an acquired gene: stays with the DIAMOND screen, never tagged X
                _report_row("ABC_01527", "contig_7", 82960, 83784, "blaOXA-66", "AMR", "BETA-LACTAM", "CARBAPENEM"),
                # a STRESS point mutation (--plus): not a resistance mutation, never tagged X
                _report_row("ABC_00440", "contig_2", 50000, 51000, "xyz_A1V", "POINT", "ACID", "ACID", "STRESS"),
            ]
            (out / "amrfinder" / "S1.tsv").write_text("\n".join("\t".join(r) for r in rows) + "\n")

            ann = pd.DataFrame(
                {
                    "Sample": ["S1", "S1", "S1", "S1"],
                    "Contig": [2, 32, 7, 2],
                    "Start": [43837, 10756, 82960, 50000],
                    "Stop": [46551, 12975, 83784, 51000],
                    "Strand": ["+", "+", "-", "+"],
                    "Locus": ["ABC_00437", "ABC_03473", "ABC_01527", "ABC_00440"],
                    "Gene": ["gyrA", "parC", "blaOXA-66", "xyz"],
                    "MarkerCategory": [None, None, "AMR", None],
                }
            )
            ann = _finalize_annotations_df(AF.attach_resistance_mutations(ann, out))

        self.assertEqual(ann["ResistanceMutations"].tolist(), ["gyrA_S81L", "parC_S84L", "", ""])
        self.assertEqual(ann["ResistanceMutationClasses"].tolist()[:2], ["quinolone", "quinolone"])
        self.assertEqual(ann["Tags"].tolist(), ["X", "X", "A", ""])
        note = mutation_note(ann.iloc[0].to_dict())
        self.assertIn("Usually not a resistance gene", note)
        self.assertIn("gyrA_S81L (quinolone)", note)

    def test_organism_from_mlst_and_skip_without_one(self):
        with tempfile.TemporaryDirectory() as tmp:
            out, db = Path(tmp) / "run", Path(tmp) / "db"
            (out / "mlst").mkdir(parents=True)
            db.mkdir()
            (db / "taxgroup.tsv").write_text(
                "#taxgroup\tgpipe_taxgroup\tnumber_of_nucl_ref_genes\n"
                "Acinetobacter_baumannii\tAcinetobacter\t1\nEscherichia\tEscherichia_coli_Shigella\t6\n"
            )
            mlst = out / "mlst" / "results.json"
            with mock.patch.object(AF, "configured_amrfinder_organism", return_value=""):
                mlst.write_text(json.dumps([{"scheme": "abaumannii_2"}] * 3 + [{"scheme": "ecoli"}]))
                self.assertEqual(AF.resolve_amrfinder_organism(out, db)[0], "Acinetobacter_baumannii")
                mlst.write_text(json.dumps([{"scheme": "bsubtilis"}]))
                organism, note = AF.resolve_amrfinder_organism(out, db)
        self.assertIsNone(organism)
        self.assertIn("no AMRFinderPlus organism", note)


if __name__ == "__main__":
    unittest.main()
