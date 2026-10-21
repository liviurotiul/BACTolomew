"""run_state.json lists the planned steps and counts a per-genome step's genomes.

Qt-free: the Runner's real run loop on tiny fixture files, with the tools stubbed out
(Bakta reports its genomes, MLST fails so the run stops there). No container is started.
"""
from __future__ import annotations

import json
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pipeline import Runner  # noqa: E402
from pipeline.DataInput import DataInput  # noqa: E402
from pipeline.RunSteps import PIPELINE_STEPS  # noqa: E402


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


class RunStateProgressTest(unittest.TestCase):
    def _run(self, report_genomes, flush_interval: float) -> tuple[dict, list[dict]]:
        """Run the real loop with Bakta calling report_genomes(progress, genomes, out);
        returns the final run_state and every run_state written while running."""
        with tempfile.TemporaryDirectory(prefix="bactolomew_run_state_") as tmp:
            root = Path(tmp)
            genomes = [_write(root / "genomes" / f"g{i}.fna", ">c1\nACGTACGT\n") for i in range(5)]
            (root / "db" / "bakta").mkdir(parents=True)
            out = root / "out"
            cfg = DataInput(
                genomes=genomes,
                marker_sequences=[_write(root / "db" / "m_AMR.fasta", ">p\nMKV\n")],
                reference_locus_sequences=[_write(root / "db" / "l_IS.fasta", ">n\nACGT\n")],
                output=out,
                snippy_reference=None,
                bakta_db=str(root / "db" / "bakta"),
                dmnd_threads=1,
                bakta_threads=1,
                window_length=500,
                diamond=True,
                bakta=True,
                total_cores=1,
            )

            def fake_bakta(CFG, *, bakta_threads, progress=None):
                report_genomes(progress, genomes, out)
                return True

            state_writes: list[dict] = []
            real_write_json = Runner._write_json

            def recording_write_json(path, data):
                if Path(path).name == Runner.RUN_STATE_FILENAME:
                    state_writes.append(json.loads(json.dumps(data, default=str)))
                real_write_json(path, data)

            with mock.patch.multiple(
                Runner,
                PROGRESS_FLUSH_INTERVAL=flush_interval,
                _write_json=recording_write_json,
                _collect_versions=lambda tools: {},
                _validate_assembly_qc_reuse_outputs=lambda CFG: True,
                _validate_bakta_reuse_outputs=lambda CFG: False,
                run_bakta=fake_bakta,
                validate_bakta_outputs=lambda *a: [],
                ensure_panaroo_out=lambda *a: None,
                ensure_chewbbaca_out=lambda *a: None,
                _validate_mlst_reuse_outputs=lambda CFG: False,
                run_mlst=mock.Mock(side_effect=RuntimeError("mlst stub failure")),
            ):
                with self.assertRaises(RuntimeError):
                    Runner._run_workflow_impl(cfg)

            state = json.loads((out / Runner.RUN_STATE_FILENAME).read_text(encoding="utf-8"))
            return state, state_writes

    @staticmethod
    def _bakta_counts(state_writes: list[dict]) -> list[int]:
        return [
            w["steps"]["bakta"]["progress"]["done"]
            for w in state_writes
            if w.get("current_step") == "bakta" and "progress" in w["steps"]["bakta"]
        ]

    def test_planned_steps_and_genome_progress(self) -> None:
        def report(progress, genomes, out):
            for i, genome in enumerate(genomes):
                progress(i, len(genomes))
                _write(out / "bakta" / genome.stem / f"{genome.stem}.json", "{}")
            progress(len(genomes), len(genomes))

        state, state_writes = self._run(report, flush_interval=3600.0)
        self.assertEqual(
            [s["name"] for s in state["planned_steps"]], [name for name, _ in PIPELINE_STEPS]
        )
        self.assertTrue(all(s["label"] for s in state["planned_steps"]))
        self.assertEqual(state["status"], "failed")
        self.assertEqual(state["steps"]["assembly_qc"]["reason"], "reused_existing_output")
        self.assertEqual(state["steps"]["bakta"]["status"], "success")
        self.assertEqual(state["steps"]["bakta"]["progress"], {"done": 5, "total": 5, "unit": "genomes"})
        self.assertEqual(state["steps"]["mlst"]["status"], "failed")
        # While Bakta ran, the throttle let through only its first and last count.
        self.assertEqual(self._bakta_counts(state_writes), [0, 5])

    def test_held_back_count_is_written_while_a_slow_genome_runs(self) -> None:
        # Three quick (reused) genomes, then a slow one: the held-back count 3/5 must reach
        # run_state.json during the slow genome, not only after it.
        def report(progress, genomes, out):
            for i, genome in enumerate(genomes):
                progress(i, len(genomes))
                if i == 3:
                    time.sleep(0.6)
                _write(out / "bakta" / genome.stem / f"{genome.stem}.json", "{}")
            progress(len(genomes), len(genomes))

        state, state_writes = self._run(report, flush_interval=0.2)
        # 1 and 2 stay held back (superseded by 3); 4 comes more than an interval later.
        self.assertEqual(self._bakta_counts(state_writes), [0, 3, 4, 5])
        self.assertEqual(state["steps"]["bakta"]["progress"]["done"], 5)


if __name__ == "__main__":
    unittest.main()
