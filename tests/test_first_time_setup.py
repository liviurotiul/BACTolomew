"""First-time setup steps (pipeline.Setup.SetupJob, no window): what is installed is
found and kept (a light Bakta DB installed earlier too), a missing Bakta DB is installed
(full) through DbDownload, DB/ is made with empty marker and compound catalogs and our
cgMLST scheme catalog (no database bundle), and too little free space for the full Bakta
DB fails its row before anything downloads, then Retry goes on from there. Also: the bar
counts a failed step as nothing, and stopping a streamed command ends what it started too."""
from __future__ import annotations

import json
import sys
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest import mock

from pipeline import Setup as fts
from pipeline import ContainerRuntime as CR
from pipeline import DbDownload
from pipeline.ToolEnv import ContainerContext


def _bakta_db(folder: Path, kind: str = "full") -> Path:
    db = folder / ("db" if kind == "full" else "db-light")
    db.mkdir(parents=True)
    (db / "version.json").write_text(json.dumps({"major": 5, "minor": 1, "type": kind}))
    (db / "bakta.db").write_bytes(b"x")
    return db


def _fake_stream(cmd, log, **_kw):
    """The tools image: bakta_db download writes a full DB, everything else succeeds."""
    if "bakta_db" in cmd:
        assert cmd[cmd.index("--type") + 1] == "full"
        db = _bakta_db(Path(cmd[cmd.index("--output") + 1]))
        rel = db / "amrfinderplus-db" / "2025-07-16.1"
        rel.mkdir(parents=True)
        (rel / "version.txt").write_text("2025-07-16.1\n")
        (db / "amrfinderplus-db" / "latest").symlink_to("2025-07-16.1")
        log(" 45%|####      | 630M/1.39G [00:20<00:25, 31.5MB/s]")
    return 0


class SetupJobTest(unittest.TestCase):
    def setUp(self) -> None:
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)  # a fresh project: no DB/ yet
        self.checkm2 = self.root / "checkm2_db" / fts.CHECKM2_DB_RELATIVE_PATH
        self.checkm2.parent.mkdir(parents=True)
        self.checkm2.write_bytes(b"dmnd")
        self.config = {"runtime": "docker", "container": {"image": "img:1"}, "keep": "me"}
        (self.root / "bactolomew_config.json").write_text(json.dumps(self.config))
        ctx = ContainerContext(runtime="docker", image="img:1", bind_roots=(str(self.root),), home_dir=str(self.root))
        for patcher in (
            mock.patch.object(CR, "runtime_status", return_value=CR.RuntimeStatus("docker", True)),
            mock.patch.object(CR, "docker_image_present", return_value=True),
            mock.patch.object(DbDownload, "_context", return_value=ctx),
            mock.patch.object(DbDownload, "_stream", side_effect=_fake_stream),
        ):
            patcher.start()
            self.addCleanup(patcher.stop)
        self.fastgear = mock.patch.object(fts, "_ensure_fastgear_image").start()
        self.addCleanup(mock.patch.stopall)
        self.steps: list[fts.SetupStep] = []

    def _job(self) -> fts.SetupJob:
        return fts.SetupJob(self.root, on_step=self.steps.append, log=lambda _line: None)

    def _states(self, job: fts.SetupJob) -> dict[str, str]:
        return {s.key: s.state for s in job.steps}

    def test_installed_parts_are_kept_and_bakta_full_and_db_folder_are_made(self) -> None:
        job = self._job()
        self.assertEqual([s.key for s in job.steps], ["image", "bakta", "checkm2", "check"])
        with mock.patch.object(DbDownload, "_free_bytes", return_value=1 << 50), \
                mock.patch.object(DbDownload, "install_checkm2") as checkm2:
            self.assertEqual(job.run(), "done")
        checkm2.assert_not_called()
        self.assertEqual(set(self._states(job).values()), {"done"})
        self.assertIn("45%", [s.text for s in self.steps if s.key == "bakta"])
        record = DbDownload.read_install_record(self.root / "bakta_db")
        self.assertEqual(record["database"], "bakta-full")
        config = json.loads((self.root / "bactolomew_config.json").read_text())
        self.assertEqual(config["bakta_db"], str(self.root / "bakta_db" / "db"))
        self.assertEqual(config["assembly_qc"]["checkm2_database_path"], str(self.checkm2))
        self.assertEqual(config["keep"], "me")
        self.fastgear.assert_not_called()  # container runtime: fastGEAR is in the tools image
        db = self.root / "DB"
        for catalog in ("marker_db_catalog.json", "compound_db_catalog.json"):
            self.assertEqual(json.loads((db / "catalogs" / catalog).read_text()), {})
        schemes = json.loads((db / "cgmlst" / fts.CGMLST_SCHEME_CATALOG_FILENAME).read_text())["schemes"]
        self.assertIn("listeria_monocytogenes", {e["schema_id"] for e in schemes})
        self.assertEqual(config["chewbbaca"]["catalog_path"], str(db / "cgmlst" / fts.CGMLST_SCHEME_CATALOG_FILENAME))
        self.assertEqual(fts.collect_setup_missing_paths(self.root, config), [])

        # Everything installed: found again, nothing downloaded.
        with mock.patch.object(DbDownload, "install_bakta") as bakta:
            self.assertEqual(self._job().run(), "done")
        bakta.assert_not_called()

        # A light Bakta DB installed earlier and configured stays the one runs use.
        light = _bakta_db(self.root / "old_bakta", "light")
        config["bakta_db"] = str(light)
        (self.root / "bactolomew_config.json").write_text(json.dumps(config))
        self.assertEqual(self._job()._find_bakta(), light)

    def test_short_of_space_fails_before_any_download_and_retry_goes_on_from_there(self) -> None:
        job = self._job()
        with mock.patch.object(DbDownload, "_free_bytes", return_value=2 << 30), \
                mock.patch.object(CR, "docker_image_present", return_value=False), \
                mock.patch.object(fts.SetupJob, "_run_image", return_value="img:1") as image:
            self.assertEqual(job.run(), "failed")  # the missing tool image is not pulled first
        image.assert_not_called()
        self.assertEqual(self._states(job), {"image": "waiting", "bakta": "failed", "checkm2": "done",
                                             "check": "waiting"})
        bakta = next(s for s in job.steps if s.key == "bakta")
        self.assertEqual(bakta.note, "Needs about 121 GB free to unpack (2.1 GB free)")
        self.assertFalse((self.root / "bakta_db").exists())

        self.assertLess(fts.setup_fraction(job.steps), 1.0)
        failed_bar = fts.setup_fraction([fts.SetupStep("a", "A", state=fts.STEP_FAILED, fraction=0.9)])
        self.assertEqual(failed_bar, 0.0)

        with mock.patch.object(DbDownload, "_free_bytes", return_value=1 << 50), \
                mock.patch.object(fts.SetupJob, "_run_image", return_value="img:1"):
            self.assertEqual(job.run(), "done")
        self.assertEqual(set(self._states(job).values()), {"done"})


class StopStreamTest(unittest.TestCase):
    def test_stop_ends_a_child_that_ignores_sigterm(self) -> None:
        """`conda run` (or docker) dies on SIGTERM but the tool it started ignores it: the
        stop reaches the whole process group and SIGKILLs what is left after the grace."""
        script = ("import signal, subprocess, sys, time\n"
                  "subprocess.Popen([sys.executable, '-c', 'import signal, time; "
                  "import os; signal.signal(signal.SIGTERM, signal.SIG_IGN); print(\"child\", os.getpid(), flush=True); "
                  "time.sleep(60)'])\n"
                  "time.sleep(60)\n")
        lines: list[str] = []
        result: dict[str, int] = {}

        def run() -> None:
            result["rc"] = CR._stream([sys.executable, "-c", script], lines.append)

        worker = threading.Thread(target=run)
        worker.start()
        deadline = time.monotonic() + 20
        while not any(line.startswith("child ") for line in lines) and time.monotonic() < deadline:
            time.sleep(0.02)
        child = next(int(line.split()[1]) for line in lines if line.startswith("child "))
        started = time.monotonic()
        self.assertTrue(CR.stop_stream(worker.ident, grace=0.5))
        worker.join(15)
        self.assertFalse(worker.is_alive())
        self.assertLess(time.monotonic() - started, 10)
        self.assertNotEqual(result["rc"], 0)
        self.assertFalse(_running(child))


def _running(pid: int) -> bool:
    """`pid` still runs (a zombie waiting for its new parent to reap it does not count)."""
    try:
        state = Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()[0]
    except (OSError, IndexError):
        return False
    return state != "Z"


if __name__ == "__main__":
    unittest.main()
