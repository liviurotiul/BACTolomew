"""Optional marker and compound databases: a run searches the installed ones, records the
others as not installed, and a run without any still validates; the GUI derives the
"not installed" tags from that record."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from pipeline.DataInput import DataInput
from pipeline.DbCatalog import run_database_record
from pipeline.RunValidation import STEP_DIAMOND_MARKERS, STEP_DIAMOND_REFERENCE, validate_completed_run
from pipeline.Setup import _catalog_targets_dir


def _write(path: Path, text: str = ">a\nACGT\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


class OptionalDatabasesTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self.markers = self.root / "DB" / "seq" / "markers"
        self.catalog = _write(self.root / "DB" / "catalogs" / "marker_db_catalog.json", json.dumps({
            "CARD_AMR": {"category": "AMR", "db_path": "../seq/markers/CARD_AMR.fasta"},
            "VFDB-setB-Prot_VIRULENCE": {"category": "VIRULENCE", "db_path": "../seq/markers/VFDB-setB-Prot_VIRULENCE.fasta"},
        }))
        self.card = _write(self.markers / "CARD_AMR.fasta", ">p\nMKV\n")  # VFDB is not installed

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def test_setup_accepts_a_catalog_whose_databases_are_missing(self) -> None:
        self.card.unlink()
        self.assertTrue(_catalog_targets_dir(self.catalog, self.markers))
        self.assertFalse(_catalog_targets_dir(self.catalog, self.root / "elsewhere"))

    def test_run_records_used_and_not_installed_databases(self) -> None:
        record = run_database_record([self.card], self.catalog)
        (used,) = record["used"]
        self.assertEqual((used["name"], used["category"], used["size"]), ("CARD_AMR", "AMR", self.card.stat().st_size))
        self.assertEqual(len(used["sha256"]), 64)
        self.assertEqual([db["name"] for db in record["not_installed"]], ["VFDB-setB-Prot_VIRULENCE"])

    def test_run_without_marker_or_compound_databases(self) -> None:
        genome = _write(self.root / "g" / "s1.fna")
        bakta_db = self.root / "bakta"
        bakta_db.mkdir()
        cfg = DataInput(
            genomes=[genome], marker_sequences=[], reference_locus_sequences=[], output=self.root / "out",
            snippy_reference=None, bakta_db=str(bakta_db), dmnd_threads=1, bakta_threads=1, window_length=500,
            diamond=True, bakta=True, total_cores=1,
        )
        self.assertEqual((cfg.marker_sequences, cfg.reference_locus_sequences), ([], []))
        run_cfg = {"inputs": {"genomes": [str(genome)], "marker_sequences": [], "reference_locus_sequences": []}}
        steps = {name: {"status": "skipped"} for name in (STEP_DIAMOND_MARKERS, STEP_DIAMOND_REFERENCE)}
        validations, _errors = validate_completed_run(self.root / "out", run_cfg=run_cfg, steps=steps)
        self.assertTrue(validations[STEP_DIAMOND_MARKERS]["outputs_validated"])
        self.assertTrue(validations[STEP_DIAMOND_REFERENCE]["outputs_validated"])
        # With a marker database the step must have run.
        run_cfg["inputs"]["marker_sequences"] = [str(self.card)]
        validations, _errors = validate_completed_run(self.root / "out", run_cfg=run_cfg, steps=steps)
        self.assertFalse(validations[STEP_DIAMOND_MARKERS]["outputs_validated"])

    def test_gui_tags_follow_the_recorded_databases(self) -> None:
        from annotation_explorer_polars.model.run_data import load_run_tag_databases

        out = self.root / "run"
        _write(out / "run_config.json", json.dumps({"inputs": {"databases": {
            "marker": {
                "used": [{"name": "CARD_AMR", "category": "AMR"}, {"name": "CGE-VF", "category": "VIRULENCE"}],
                "not_installed": [{"name": "VFDB-setB-Prot_VIRULENCE", "category": "VIRULENCE"},
                                  {"name": "Phaster-Prophage-Virus-Prot_PHAGE", "category": "PHAGE"}],
            },
            "compound": {"used": [], "not_installed": [{"name": "PAIDB", "category": ""}]},
        }}}))
        tags = load_run_tag_databases(out)
        self.assertEqual(tags["A"].note(), "")
        self.assertEqual(tags["V"].note(), "partly installed")
        self.assertEqual(tags["P"].note(), "not installed")
        self.assertEqual(tags["I"].note(), "not installed")
        self.assertIn("VFDB-setB-Prot_VIRULENCE", tags["V"].tooltip())
        # A run made before the databases were optional records none: nothing is marked.
        _write(out / "run_config.json", json.dumps({"inputs": {}}))
        self.assertEqual(load_run_tag_databases(out), {})


if __name__ == "__main__":
    unittest.main()
