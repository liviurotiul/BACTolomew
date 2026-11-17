import json
import os
import subprocess
import sys
import tempfile
import unittest
from datetime import datetime, timedelta
from pathlib import Path

from annotation_explorer_polars.run_expectations import load_run_status, run_is_active
from pipeline.RunLock import lock_payload


class StaleRunLockTest(unittest.TestCase):
    """A run.lock.json whose process is gone, or whose pid now belongs to another process,
    is removed and the run reads as interrupted; a live lock or another host's is kept."""

    def _run_dir(self, lock: dict) -> Path:
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        run_dir = Path(tmp.name)
        (run_dir / "run_state.json").write_text(json.dumps({"status": "running"}))
        (run_dir / "run.lock.json").write_text(json.dumps(lock))
        return run_dir

    def test_stale_and_live_locks(self):
        live = lock_payload(run_id="r", start_time="t")
        dead = subprocess.Popen([sys.executable, "-c", "pass"])
        dead.wait()
        other = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
        self.addCleanup(other.kill)
        hour_ago = (datetime.now() - timedelta(hours=1)).astimezone().isoformat()
        old_live = {"pid": other.pid, "host": live["host"], "start_time": datetime.now().astimezone().isoformat()}
        stale = {
            "dead pid": {**live, "pid": dead.pid},
            "old lock with this app's pid": {"pid": os.getpid(), "host": live["host"]},
            "reused pid": {**live, "pid": other.pid, "process_start": live["process_start"] - 3600},
            # 0.3.0 locks (no process_start): the pid's process started after the run began.
            "0.3.0 lock, reused pid": {"pid": other.pid, "host": live["host"], "start_time": hour_ago},
        }
        for label, lock in stale.items():
            with self.subTest(label):
                run_dir = self._run_dir(lock)
                with self.assertLogs("annotation_explorer_polars.run_expectations", "WARNING"):
                    self.assertEqual(load_run_status(run_dir)["status"], "interrupted")
                self.assertFalse((run_dir / "run.lock.json").exists())
        for label, lock in {"live": live, "0.3.0 live": old_live, "other host": {**live, "pid": dead.pid, "host": "elsewhere"}}.items():
            with self.subTest(label):
                run_dir = self._run_dir(lock)
                self.assertTrue(run_is_active(run_dir))
                self.assertTrue((run_dir / "run.lock.json").exists())


if __name__ == "__main__":
    unittest.main()
