"""Genome sample names must be unique (P05) and a run folder keeps its genome set (P06).

Qt-free: DataInput and the Runner's start-of-run guard, on tiny fixture files in a
temporary folder. No tool or container is started.
"""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pipeline.DataInput import DataInput  # noqa: E402

FASTA = ">c1\nACGTACGT\n"


def _write(path: Path, text: str = FASTA) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


class _TmpCase(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(prefix="bactolomew_genome_set_")
        self.root = Path(self._tmp.name)

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def data_input(self, genomes, output: Path | None = None) -> DataInput:
        marker = _write(self.root / "db" / "markers" / "m_AMR.fasta", ">p\nMKV\n")
        locus = _write(self.root / "db" / "loci" / "l_IS.fasta", ">n\nACGT\n")
        bakta_db = self.root / "db" / "bakta"
        bakta_db.mkdir(parents=True, exist_ok=True)
        return DataInput(
            genomes=list(genomes),
            marker_sequences=[marker],
            reference_locus_sequences=[locus],
            output=output or (self.root / "out"),
            snippy_reference=None,
            bakta_db=str(bakta_db),
            dmnd_threads=1,
            bakta_threads=1,
            window_length=500,
            diamond=True,
            bakta=True,
            total_cores=1,
        )


class GenomeNamesTest(_TmpCase):
    def test_duplicate_sample_names_are_refused(self) -> None:
        a = _write(self.root / "one" / "s1.fna")
        b = _write(self.root / "two" / "s1.fasta")
        c = _write(self.root / "s2.fna")
        with self.assertRaises(ValueError) as cm:
            self.data_input([a, b, c])
        msg = str(cm.exception)
        self.assertIn(str(a), msg)
        self.assertIn(str(b), msg)
        self.assertNotIn(str(c), msg)
        self.assertEqual([p.stem for p in self.data_input([a, c]).genomes], ["s1", "s2"])
        # Names that differ only by case would share a folder on macOS and Windows.
        upper = _write(self.root / "three" / "S2.fna")
        with self.assertRaises(ValueError) as cm:
            self.data_input([a, c, upper])
        self.assertIn(str(upper), str(cm.exception))


class GenomeSetTest(_TmpCase):
    def test_runner_refuses_a_changed_genome_set_before_writing_anything(self) -> None:
        from pipeline import Runner

        genomes = [_write(self.root / "g" / f"{s}.fna") for s in ("a", "b", "c")]
        out = self.root / "run"
        out.mkdir()
        payload = {"inputs": {"genomes": [str(g) for g in genomes[:2]]}}
        (out / "run_config.json").write_text(json.dumps(payload), encoding="utf-8")
        before = (out / "run_config.json").read_bytes()
        with self.assertRaises(Runner.GenomeSetChangedError) as cm:
            Runner._run_workflow_impl(self.data_input(genomes, output=out))
        self.assertIn("Added: c.", str(cm.exception))
        self.assertEqual((out / "run_config.json").read_bytes(), before)
        self.assertEqual(sorted(p.name for p in out.iterdir()), ["run_config.json"])


if __name__ == "__main__":
    unittest.main()
