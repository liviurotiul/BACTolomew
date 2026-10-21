"""cgMLST concepts: schemas are adapted with PrepExternalSchema so every allele is a
forward-strand CDS (P01), and an AlleleCall without a self-score for every
representative fails (P03).

The guard test mocks AlleleCall; the schema test runs the real PrepExternalSchema in
the tools image and is skipped without docker or the image.
"""
from __future__ import annotations

import logging
import pickle
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from pipeline.ChewBBACAWorkflow import load_self_scores
from pipeline.ChewbbacaSchemaBuild import (
    BuildOptions,
    adapted_schema_issues,
    build_schema,
    reverse_complement,
    schema_stats,
)

# A 12-codon CDS (ATG ... TAA), a second allele of the same locus, and a non-CDS.
CDS_A = "ATGAAACCCGGGTTTAAACCCGGGAAATTTGCTTAA"
CDS_B = "ATGAAACCCGGGTTTAAACCCGGGAAATTTGCATAA"
NOT_CDS = "ATGAAATAGCCCGGGTTTAAACCCGGGAAATTTTAA"  # internal stop TAG

def _fasta(records: list[tuple[str, str]]) -> str:
    return "".join(f">{h}\n{s}\n" for h, s in records)


def _write_schema(root: Path, loci: dict[str, list[tuple[str, str]]], reps: dict[str, list[tuple[str, str]]],
                  *, table: int = 11) -> Path:
    (root / "short").mkdir(parents=True, exist_ok=True)
    for locus, records in loci.items():
        (root / f"{locus}.fasta").write_text(_fasta(records), encoding="utf-8")
    for locus, records in reps.items():
        (root / "short" / f"{locus}_short.fasta").write_text(_fasta(records), encoding="utf-8")
    (root / ".genes_list").write_bytes(pickle.dumps([f"{locus}.fasta" for locus in loci]))
    (root / ".schema_config").write_bytes(pickle.dumps({"translation_table": [table]}))
    return root


def _write_self_scores(schema: Path, scores: dict) -> None:
    (schema / "short" / "self_scores").write_bytes(pickle.dumps(scores))


class _TmpCase(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)

    def tearDown(self) -> None:
        self._tmp.cleanup()


class RunChewbbacaGuardTests(_TmpCase):
    """run_chewbbaca applies the guard around a (mocked) AlleleCall."""

    def _run(self, stdout: str, scores_after: dict | None):
        from unittest import mock

        import pipeline.ChewBBACAWorkflow as wf

        base = Path(tempfile.mkdtemp(dir=self.tmp))
        schema = _write_schema(base / "schema", {"L1": [("L1_1", CDS_A)]}, {"L1": [("L1_1", CDS_A)]})
        _write_self_scores(schema, {"1": (813, 1451.0)})  # degenerate cache from an old build
        genome = base / "g1.fna"
        genome.write_text(">c\nACGT\n")
        seen_before: list = []

        def fake_run(cmd, **_kw):
            seen_before.append(load_self_scores(schema))
            out = Path(cmd[cmd.index("-o") + 1])
            out.mkdir(parents=True)
            (out / "results_alleles.tsv").write_text("FILE\tL1.fasta\ng1\t1\n")
            if scores_after is not None:
                _write_self_scores(schema, scores_after)
            return subprocess.CompletedProcess(cmd, 0, stdout, "")

        logger = logging.getLogger("test_run_chewbbaca_guard")
        logger.addHandler(logging.NullHandler())
        logger.propagate = False
        config = {"schema_dir": str(schema), "auto_select_from_mlst": False}
        with mock.patch.object(wf, "_resolve_chewbbaca_base_cmd", return_value=["chewBBACA.py"]), \
                mock.patch.object(wf.subprocess, "run", side_effect=fake_run):
            try:
                return wf.run_chewbbaca([genome], base / "run", logger, config=config), seen_before, schema
            except RuntimeError as exc:
                return exc, seen_before, schema

    def test_self_score_mismatch_fails_the_step(self) -> None:
        # A degenerate cache from an old build is dropped before AlleleCall.
        result, seen_before, schema = self._run("ok\n", {"L1_1": (12, 60.0)})
        self.assertEqual(seen_before, [None])
        self.assertIsInstance(result, Path)
        # Fewer self-scores than representatives: the step fails, the cache is dropped.
        result, _seen, schema = self._run("ok\n", {"1": (813, 1451.0)})
        self.assertIsInstance(result, RuntimeError)
        self.assertIn("1 representatives have none", str(result))
        self.assertIsNone(load_self_scores(schema))
        # chewBBACA logging a missing self-score fails the step too.
        result, _seen, schema = self._run(
            "Could not get the self-score for the representative allele L1_1 x\n", {"L1_1": (12, 60.0)}
        )
        self.assertIsInstance(result, RuntimeError)
        self.assertIsNone(load_self_scores(schema))


def _docker_image_available() -> bool:
    if shutil.which("docker") is None:
        return False
    try:
        res = subprocess.run(["docker", "image", "inspect", "bactolomew-tools:0.1.0"],
                             capture_output=True, check=False, timeout=30)
    except Exception:
        return False
    return res.returncode == 0


@unittest.skipUnless(_docker_image_available(), "tools image not available")
class RealPrepExternalSchemaTests(_TmpCase):
    def test_reverse_strand_alleles_are_oriented(self) -> None:
        ext = self.tmp / "ridom"
        ext.mkdir()
        (ext / "L1.fasta").write_text(_fasta([("1", CDS_A), ("2", reverse_complement(CDS_B)), ("3", NOT_CDS)]))
        (ext / "L2.fasta").write_text(_fasta([("1", reverse_complement(CDS_A))]))
        target = self.tmp / "db" / "species"
        result = build_schema(ext, target, options=BuildOptions(cpu=1))
        self.assertEqual(result.written_loci, 2)
        self.assertEqual(result.rejected_alleles, 1)
        stats = schema_stats(target)
        self.assertEqual((stats.alleles, stats.alleles_forward, stats.reps_reverse_only), (3, 3, 0))
        self.assertEqual(adapted_schema_issues(target), [])


if __name__ == "__main__":
    unittest.main()
