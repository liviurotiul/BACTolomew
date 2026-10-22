"""The start window's run cards and the Run progress window show steps, not files.

Synthetic run folders: a running run with planned_steps and a per-genome step
mid-way, a finished one, a stopped one, and an old run_state without planned_steps.
"""
from __future__ import annotations

import json
import os
import socket
import tempfile
import unittest
from datetime import datetime, timedelta, timezone
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

from PyQt6 import QtWidgets

from annotation_explorer_polars import startup_dialog
from annotation_explorer_polars.run_progress import RunProgressDialog
from annotation_explorer_polars.run_steps import format_duration, summarize_run_steps
from annotation_explorer_polars.theme import apply_application_theme
from pipeline.RunSteps import planned_steps

NOW = datetime.now(timezone.utc)


def _iso(minutes_ago: float) -> str:
    return (NOW - timedelta(minutes=minutes_ago)).isoformat()


def _done(start: float, end: float, **extra) -> dict:
    step = {"status": "success", "start_time": _iso(start), "end_time": _iso(end), "outputs_validated": True}
    step.update(extra)
    return step


def _running_state() -> dict:
    return {
        "start_time": _iso(72),
        "status": "running",
        "current_step": "bakta",
        "planned_steps": planned_steps(),
        "steps": {
            "assembly_qc": _done(72, 60, reason="reused_existing_output"),
            "bakta": {"status": "running", "start_time": _iso(14), "progress": {"done": 37, "total": 153, "unit": "genomes"}},
        },
    }


class RunStepsModelTests(unittest.TestCase):
    def test_running_step_line_and_bar(self) -> None:
        steps = summarize_run_steps(_running_state(), "running", now=NOW)
        self.assertEqual(steps.total, 16)
        self.assertEqual(steps.finished, 1)
        self.assertEqual(steps.step_line(), "Step 2 of 16 · Bakta annotation · 14 min · 37/153 genomes")
        self.assertAlmostEqual(steps.fraction, (1 + 37 / 153) / 16)
        self.assertEqual(steps.rows[0].text, "reused · 12 min")
        self.assertEqual(steps.rows[2].text, "waiting")
        self.assertEqual(steps.header_text(), "Step 2 of 16 · 1 h 12 min")

    def test_stopped_run_marks_the_started_step(self) -> None:
        state = _running_state()
        state.update(status="stopped", end_time=_iso(1))
        steps = summarize_run_steps(state, "stopped", now=NOW)
        self.assertEqual(steps.rows[1].text, "stopped · 37/153 genomes")
        self.assertEqual(steps.rows[2].text, "not run")
        self.assertEqual(steps.header_text(), "1 of 16 steps · 1 h 11 min")

    def test_failed_and_skipped_details(self) -> None:
        state = {
            "status": "failed",
            "planned_steps": planned_steps(),
            "steps": {
                "amrfinder_mutations": {"status": "skipped", "reason": "skipped: no database"},
                "gubbins": {"status": "failed", "error": "Gubbins exited with status 1"},
                "phylogeny": {"status": "failed", "error": "snippy failed",
                              "progress": {"done": 77, "total": 152, "unit": "genomes"}},
                "core_tree": _done(12, 5, reason="reused_existing_output"),
                "core_alignment": _done(13, 12.5, reason="reused_existing_output"),
            },
        }
        rows = {row.name: row for row in summarize_run_steps(state, "failed").rows}
        self.assertEqual((rows["amrfinder_mutations"].text, rows["amrfinder_mutations"].detail), ("skipped", "no database"))
        self.assertEqual((rows["gubbins"].text, rows["gubbins"].detail), ("failed", "Gubbins exited with status 1"))
        # A failed per-genome step keeps its count, like a stopped one.
        self.assertEqual(rows["phylogeny"].text, "failed · 77/152 genomes")
        # A reused step that took a while (the core tree's refit) shows how long.
        self.assertEqual(rows["core_tree"].text, "reused · 7 min")
        self.assertEqual(rows["core_alignment"].text, "reused")

    def test_old_run_state_without_planned_steps(self) -> None:
        state = {"steps": {"mlst": _done(5, 4), "bakta": _done(10, 6), "assembly_qc": _done(12, 10)}}
        finished = summarize_run_steps(state, "success")
        self.assertEqual([row.name for row in finished.rows], ["assembly_qc", "bakta", "mlst"])
        self.assertEqual(finished.fraction, 1.0)
        stopped = summarize_run_steps(state, "stopped")
        self.assertEqual(stopped.total, 16)  # the steps it has not reached are listed too
        self.assertEqual(stopped.finished, 3)
        # A known step before the last one reached was not in that run's pipeline.
        older = {"steps": {"mlst": _done(5, 4), "mobsuite": _done(4, 3)}}
        names = [row.name for row in summarize_run_steps(older, "stopped").rows]
        self.assertEqual(names[:2], ["mlst", "mobsuite"])
        self.assertNotIn("amrfinder_mutations", names)
        self.assertIn("panaroo", names)

    def test_format_duration(self) -> None:
        self.assertEqual(format_duration(8), "8 s")
        self.assertEqual(format_duration(14 * 60 + 20), "14 min")
        self.assertEqual(format_duration(72 * 60), "1 h 12 min")
        self.assertEqual(format_duration(120 * 60), "2 h")


class RunProgressViewsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        apply_application_theme(cls.app)

    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(prefix="bactolomew_run_progress_")
        self.root = Path(self._tmp.name)
        self._runs_root = startup_dialog.runs_root
        startup_dialog.runs_root = lambda: self.root

    def tearDown(self) -> None:
        startup_dialog.runs_root = self._runs_root
        self._tmp.cleanup()

    def _run(self, name: str, state: dict, *, lock: bool = False, recorded: str = "") -> Path:
        out = self.root / name / f"{name}__run"
        out.mkdir(parents=True)
        (out / "run_config.json").write_text(
            json.dumps({"run_id": out.name, "sample_set_id": name, "output_dir": recorded or str(out),
                        "inputs": {"genomes": ["/x/a.fna", "/x/b.fna"]}})
        )
        (out / "run_state.json").write_text(json.dumps(state))
        if lock:
            (out / "run.lock.json").write_text(json.dumps({"pid": os.getpid(), "host": socket.gethostname()}))
        return out

    def test_running_run_card_and_progress_window(self) -> None:
        out = self._run("sample_running", _running_state(), lock=True)
        dlg = startup_dialog.StartupDialog()
        try:
            dlg._refresh_timer.stop()
            entry = dlg._run_entries[0]
            self.assertTrue(entry.running)
            self.assertEqual(entry.status_text, "Running: bakta")
            card = dlg.run_list.itemWidget(dlg.run_list.item(0))
            texts = [label.text() for label in card.findChildren(QtWidgets.QLabel)]
            # No status pill: the status lives only in the text line under the bar.
            self.assertNotIn("Running", texts)
            self.assertTrue(any(t.startswith("Step 2 of 16 · Bakta annotation · ") and t.endswith("37/153 genomes") for t in texts))
            buttons = {b.text() for b in card.findChildren(QtWidgets.QPushButton)}
            self.assertEqual(buttons, {"View", "Monitor", "Resume", "Stop", "Force Stop", "More"})
        finally:
            dlg.close()

        progress = RunProgressDialog(out)
        try:
            self.assertEqual(progress.title.text(), "Run in progress")
            self.assertEqual(progress.status.text(), "Running: bakta")
            self.assertTrue(progress.summary.text().startswith("Step 2 of 16 · "))
            self.assertEqual(len(progress._step_rows), 16)
            self.assertTrue(progress._step_rows[1].state.text().endswith("37/153 genomes"))
            self.assertTrue(progress.timer.isActive())

            # The Runner moves on: the next refresh reads run_state.json again.
            state = _running_state()
            state["steps"]["bakta"] = _done(14, 1, progress={"done": 153, "total": 153, "unit": "genomes"})
            state["steps"]["mlst"] = {"status": "running", "start_time": _iso(1)}
            state["current_step"] = "mlst"
            (out / "run_state.json").write_text(json.dumps(state))
            progress._refresh()
            self.assertTrue(progress.summary.text().startswith("Step 3 of 16 · "))
            self.assertTrue(progress._step_rows[1].state.text().startswith("done · 13 min"))
            self.assertEqual(progress.progress.value(), round(2 / 16 * 1000))
        finally:
            progress.close()

    def test_launch_error_stays_after_refresh(self) -> None:
        # The runner failed before run_state.json said so: its error is not replaced.
        out = self._run("sample_launch", _running_state())
        progress = RunProgressDialog(out)
        try:
            progress.mark_failed("Could not start the container\nDocker is not running")
            progress._refresh()
            self.assertEqual(progress.title.text(), "Run failed")
            self.assertEqual(progress.status.text(), "Could not start the container")
            self.assertIn("Docker is not running", progress.details.toPlainText())
            self.assertFalse(progress.timer.isActive())
        finally:
            progress.close()

    def test_long_run_title_is_elided_with_a_tooltip(self) -> None:
        name = "sample_" + "a_very_long_sample_set_identifier_" * 3
        self._run(name, _running_state())
        dlg = startup_dialog.StartupDialog()
        try:
            dlg._refresh_timer.stop()
            dlg.resize(600, 600)
            dlg.show()
            QtWidgets.QApplication.processEvents()
            card = dlg._card(dlg.run_list.item(0))
            title = next(l for l in card.findChildren(QtWidgets.QLabel) if l.property("role") == "runTitle")
            self.assertNotEqual(title.text(), name)
            self.assertIn("…", title.text())
            self.assertEqual(title.toolTip(), name)
            self.assertEqual(card.width(), dlg.run_list.viewport().width())
        finally:
            dlg.close()

    def test_terminal_runs_keep_their_status_texts(self) -> None:
        stopped = _running_state()
        stopped.update(status="stopped", end_time=_iso(1))
        out = self._run("sample_stopped", stopped)
        progress = RunProgressDialog(out)
        try:
            self.assertEqual(progress.title.text(), "Run stopped")
            self.assertEqual(progress.status.text(), "Stopped by user. Click Resume to continue.")
            self.assertFalse(progress.timer.isActive())
        finally:
            progress.close()

        self._run("sample_old", {"status": "stopped", "steps": {"assembly_qc": _done(10, 5)}})
        dlg = startup_dialog.StartupDialog()
        try:
            dlg._refresh_timer.stop()
            for entry in dlg._run_entries:
                self.assertEqual(entry.status_text, "Stopped. Click Resume to continue.")
            old = next(e for e in dlg._run_entries if e.sample_set_id == "sample_old")
            self.assertEqual((old.steps.finished, old.steps.total), (1, 16))
        finally:
            dlg.close()

    def test_a_moved_run_is_listed_where_it_is_now(self) -> None:
        # run_config.json names the folder the run was made in, on another machine.
        out = self._run("sample_moved", {"status": "stopped"}, recorded="/elsewhere/sample_moved__run")
        dlg = startup_dialog.StartupDialog()
        try:
            dlg._refresh_timer.stop()
            entry = dlg._run_entries[0]
            self.assertEqual(entry.path, out.resolve())
            self.assertFalse(entry.shared)
            self.assertTrue(entry.metadata_label.startswith("Metadata: "))
            card = dlg._card(dlg.run_list.item(0))
            more = next(b for b in card.findChildren(QtWidgets.QPushButton) if b.text() == "More")
            texts = [a.text() for a in more.menu().actions()]
            self.assertEqual(texts, ["Run Details...", "Change Metadata...", "Export for Sharing...", "Show in Folder"])
        finally:
            dlg.close()


if __name__ == "__main__":
    unittest.main()
