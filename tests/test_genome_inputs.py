import gzip
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from pipeline.GenomeInputs import prepare_genome_inputs, read_contig_map, uses_normalized_inputs

SPADES = (
    ">NODE_1_length_12_cov_10.5 len=12 cov=10.5 sw=x\r\nACGTACGTACGT\r\n\r\n"
    ">NODE_1_length_12_cov_10.5 duplicate\r\nacgt\r\n"
    ">contig|2 (plasmid)\nGGGG\nCC\n"
)


class GenomeInputsTest(unittest.TestCase):
    def test_new_run_reads_normalized_copies_and_resume_reuses_them(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            a = root / "A.fasta"
            a.write_text(SPADES, encoding="utf-8")
            b = root / "B.fna.gz"
            with gzip.open(b, "wt") as fh:
                fh.write(">c1 desc\nAAAA\n")
            out = root / "run"
            cfg = SimpleNamespace(genomes=[a, b], snippy_reference=a, output=out, genome_inputs={})

            run = prepare_genome_inputs(cfg, {})
            self.assertEqual(run.genomes, [out / "inputs" / "A.fna", out / "inputs" / "B.fna.fna"])
            self.assertEqual([p.stem for p in run.genomes], [p.stem for p in cfg.genomes])
            self.assertEqual(run.snippy_reference, out / "inputs" / "A.fna")
            self.assertEqual(cfg.genomes, [a, b])  # the given config is not changed
            self.assertEqual(
                (out / "inputs" / "A.fna").read_text(),
                ">NODE_1_length_12_cov_10.5\nACGTACGTACGT\n>NODE_1_length_12_cov_10.5_2\nacgt\n"
                ">contig_2\nGGGG\nCC\n",
            )
            self.assertEqual((out / "inputs" / "B.fna.fna").read_text(), ">c1\nAAAA\n")
            self.assertEqual(read_contig_map(out / "inputs" / "A.contig_map.tsv")["contig|2"], "contig_2")
            sample = run.genome_inputs["samples"]["A"]
            self.assertEqual((sample["contigs"], sample["renamed"]), (3, 2))  # the first only loses its description
            self.assertEqual(len(sample["original_sha256"]), 64)

            # A resume with unchanged originals keeps the copies as they are.
            mtime = (out / "inputs" / "A.fna").stat().st_mtime_ns
            recorded = {"inputs": {"genome_inputs": run.genome_inputs}}
            self.assertTrue(uses_normalized_inputs(recorded, folder_has_run=True))
            again = prepare_genome_inputs(cfg, recorded)
            self.assertEqual((out / "inputs" / "A.fna").stat().st_mtime_ns, mtime)
            self.assertEqual(again.genome_inputs, run.genome_inputs)

    def test_run_made_before_keeps_its_original_inputs(self):
        self.assertTrue(uses_normalized_inputs({}, folder_has_run=False))
        self.assertFalse(uses_normalized_inputs({"inputs": {"genomes": ["/x/A.fasta"]}}, folder_has_run=True))


if __name__ == "__main__":
    unittest.main()
