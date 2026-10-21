import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from pipeline import DbDownload
from pipeline.ToolEnv import ContainerContext


def _fake_bakta_stream(cmd, log):
    """Stands in for the container: `bakta_db download` writes a full DB with its AMRF DB."""
    if "bakta_db" in cmd:
        out = Path(cmd[cmd.index("--output") + 1])
        db = out / "db"
        rel = db / "amrfinderplus-db" / "2025-07-16.1"
        rel.mkdir(parents=True)
        (db / "version.json").write_text(json.dumps({"major": 5, "minor": 1, "type": "full"}))
        (db / "bakta.db").write_bytes(b"x" * 10)
        (rel / "version.txt").write_text("2025-07-16.1\n")
        (rel / "database_format_version.txt").write_text("4.0.0\n")
        (db / "amrfinderplus-db" / "latest").symlink_to("2025-07-16.1")
    log("progress line")
    return 0


class DbDownloadTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        ctx = ContainerContext(runtime="docker", image="img:1", bind_roots=(str(self.root),),
                               home_dir=str(self.root))
        patcher = mock.patch.object(DbDownload, "_context", return_value=ctx)
        patcher.start()
        self.addCleanup(patcher.stop)
        space = mock.patch.object(DbDownload, "_free_bytes", return_value=1 << 50)
        space.start()
        self.addCleanup(space.stop)
        self.addCleanup(self._tmp.cleanup)

    def test_bakta_installs_via_partial_and_keeps_existing(self):
        lines = []
        target = self.root / "bakta-full"
        with mock.patch.object(DbDownload, "_stream", side_effect=_fake_bakta_stream) as stream:
            rec = DbDownload.install_bakta(target, log=lines.append)
        cmds = [c.args[0] for c in stream.call_args_list]
        self.assertIn("bakta_db", cmds[0])
        self.assertEqual(cmds[0][-4:], ["--output", str(target.with_name(".bakta-full.partial")), "--type", "full"])
        self.assertIn("--database_version", cmds[1])  # the image's amrfinder accepts the bundled AMRF DB
        self.assertEqual(rec["path"], str(target / "db"))
        self.assertEqual(rec["amrfinderplus_db"]["release"], "2025-07-16.1")
        self.assertIn("progress line", lines)
        self.assertFalse(DbDownload.partial_dir(target).exists())
        self.assertEqual(DbDownload.read_install_record(target)["version"]["major"], 5)

        # An existing database is kept (nothing runs) or a new one goes beside it.
        with mock.patch.object(DbDownload, "_stream", side_effect=_fake_bakta_stream) as stream:
            kept = DbDownload.install_bakta(target, log=lines.append)
            self.assertEqual(stream.call_count, 0)
            beside = DbDownload.install_bakta(target, if_exists="beside", log=lines.append)
        self.assertEqual(kept["path"], rec["path"])
        self.assertNotEqual(Path(beside["path"]).parent, target)
        self.assertTrue((target / "db" / "bakta.db").is_file())

    def test_failed_download_leaves_no_target_and_dry_run_writes_nothing(self):
        target = self.root / "checkm2"
        with mock.patch.object(DbDownload, "_stream", return_value=1):
            with self.assertRaises(DbDownload.DbInstallError):
                DbDownload.install_checkm2(target, log=lambda _m: None)
        self.assertFalse(target.exists())
        self.assertFalse(DbDownload.partial_dir(target).exists())

        lines = []
        with mock.patch.object(DbDownload, "_stream") as stream:
            plan = DbDownload.install_checkm2(target, dry_run=True, log=lines.append)
        stream.assert_not_called()
        self.assertEqual(plan["commands"][0][-6:-2], ["checkm2", "database", "--download", "--path"])
        self.assertEqual(list(self.root.iterdir()), [])


if __name__ == "__main__":
    unittest.main()
