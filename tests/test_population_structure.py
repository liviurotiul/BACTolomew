"""Population structure by fastGEAR or FastBAPS (pipeline/PopulationStructure.py).

The critical rules: which method a run uses from what the tools image holds (the default
is FastBAPS, auto prefers FastBAPS, a requested method the image lacks stops the run), how an image's methods are
read (its record, else the fastGEAR files of an older image), where every reader finds a
run's sample_clades.tsv (runs before FastBAPS: fastgear/, as before), and the FastBAPS
mapping (lineage = Level 1, cluster = Level 2) in the shared sample_clades.tsv format.
"""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from pipeline.FastBapsWorkflow import sample_clades_from_clusters
from pipeline.FastGearWorkflow import validate_sample_clades_tsv
from pipeline.PopulationStructure import (
    PopulationMethodUnavailable,
    _methods_from_probe,
    resolve_population_method,
    resume_phylo_overrides,
    sample_clades_path,
)


class PopulationStructureTest(unittest.TestCase):
    def test_method_choice(self) -> None:
        both = ("fastgear", "fastbaps")
        self.assertEqual(resolve_population_method("auto", both).method, "fastbaps")
        self.assertEqual(resolve_population_method("auto", ("fastgear",)).method, "fastgear")
        self.assertEqual(resolve_population_method("", both).method, "fastbaps")
        self.assertEqual(resolve_population_method("fastgear", both).method, "fastgear")
        self.assertIsNone(resolve_population_method("auto", ()).method)
        with self.assertRaises(PopulationMethodUnavailable):
            resolve_population_method("fastbaps", ("fastgear",))
        # The default (nobody set a method) on an image without FastBAPS: fastGEAR.
        self.assertEqual(resolve_population_method("fastbaps", ("fastgear",), explicit=False).method, "fastgear")
        # An image that could not be checked keeps the behavior from before: fastGEAR.
        self.assertEqual(resolve_population_method("auto", None).method, "fastgear")

    def test_image_methods(self) -> None:
        recorded = json.dumps({"mobsuite_db": {}, "population_methods": ["fastbaps"]})
        self.assertEqual(_methods_from_probe(recorded + "\n@@fastgear-files\n"), ("fastbaps",))
        # Images before the record: fastGEAR only when its files are there.
        self.assertEqual(_methods_from_probe('{"mobsuite_db": {}}\n@@fastgear-files\n'), ("fastgear",))
        self.assertEqual(_methods_from_probe("\n"), ())

    def test_run_folder(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp)
            (run / "run_config.json").write_text(json.dumps({"inputs": {"genomes": []}}))
            self.assertEqual(sample_clades_path(run), run / "fastgear" / "sample_clades.tsv")
            # A resume keeps the run's method, whatever the configured default: fastGEAR
            # for a run from before the record.
            self.assertEqual(resume_phylo_overrides(run), {"population_method": "fastgear"})
            record = {"requested": "auto", "method": "fastbaps", "available": ["fastbaps"]}
            (run / "run_config.json").write_text(json.dumps({"inputs": {"population_structure": record}}))
            self.assertEqual(sample_clades_path(run), run / "fastbaps" / "sample_clades.tsv")
            self.assertEqual(resume_phylo_overrides(run), {"population_method": "fastbaps"})
            record["method"] = None
            (run / "run_config.json").write_text(json.dumps({"inputs": {"population_structure": record}}))
            self.assertIsNone(sample_clades_path(run))

    def test_fastbaps_mapping(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            clusters = tmp / "fastbaps_clusters.tsv"
            clusters.write_text("Isolates\tLevel 1\tLevel 2\nb\t2\t3\na\t1\t1\nc\t2\t4\n")
            table = sample_clades_from_clusters(clusters, ["a", "b", "c"])
            dest = tmp / "sample_clades.tsv"
            table.to_csv(dest, sep="\t", index=False)
            rows = validate_sample_clades_tsv(dest, expected_samples=["a", "b", "c"])
        self.assertEqual(rows["Sample"].tolist(), ["a", "b", "c"])
        self.assertEqual(rows["FastGearLineage"].tolist(), ["1", "2", "2"])
        self.assertEqual(rows["FastGearCluster"].tolist(), ["1", "3", "4"])
        self.assertEqual(rows["FastGearMethod"].tolist(), ["fastbaps"] * 3)


if __name__ == "__main__":
    unittest.main()
