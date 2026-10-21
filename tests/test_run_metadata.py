"""A run opened without a metadata CSV: its own sample table (pipeline/RunMetadata.py).

Every genome of the data gets a row with what the run knows (species from the MLST
scheme, ST, fastGEAR clade and cluster, assembly QC), typed like a CSV with the
labels (ST, Clade, Cluster) kept as text.
"""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from annotation_explorer_polars.data import attach_assembly_quality_metadata, prepare_run_metadata


class RunMetadataTest(unittest.TestCase):
    def test_run_table_from_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp)
            (run / "mlst").mkdir()
            (run / "mlst" / "results.json").write_text(json.dumps([
                {"sample": "a", "scheme": "abaumannii_2", "sequence_type": "492"},
                {"sample": "b", "scheme": "abaumannii_2", "sequence_type": "-"},
            ]))
            (run / "fastgear").mkdir()
            (run / "fastgear" / "sample_clades.tsv").write_text(
                "Sample\tFastGearLineage\tFastGearCluster\nb\t1\t3\n"
            )
            (run / "assembly_qc").mkdir()
            (run / "assembly_qc" / "summary.tsv").write_text(
                "Sample\tQuality\tQuastN50\tCheckM2Completeness\na\t99\t108812\t100.0\n"
            )

            meta, warnings, types = prepare_run_metadata(run, ["a", "b"], ["c"])
            meta = attach_assembly_quality_metadata(meta, run)

        self.assertEqual(warnings, [])
        self.assertEqual(meta["Sample"].tolist(), ["a", "b", "c"])
        self.assertEqual(meta["Species"].tolist()[:2], ["Acinetobacter baumannii"] * 2)
        st = meta["ST"].tolist()
        self.assertEqual(st[0], "492")
        self.assertTrue(meta["ST"].iloc[1:].isna().all())
        self.assertEqual(meta.loc[1, "Clade"], "1")
        self.assertEqual(types["ST"], "string")
        self.assertEqual(types["Clade"], "string")
        self.assertEqual(int(meta.loc[0, "QuastN50"]), 108812)
        self.assertNotIn("MLSTScheme", meta.columns)


if __name__ == "__main__":
    unittest.main()
