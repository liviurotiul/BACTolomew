"""Run provenance (P35): the tool environment a run used is recorded in run_config.json.

Qt-free and Docker-free: `docker image inspect` and the read of the image's
tool_versions.json are faked, so no container starts.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pipeline import ContainerRuntime as CR  # noqa: E402
from pipeline import Runner  # noqa: E402
from pipeline.ToolEnv import ContainerContext, configured_image_ref  # noqa: E402

DIGEST = "sha256:" + "a" * 64
IMAGE_ID = "sha256:" + "b" * 64
REF = "ghcr.io/example/bactolomew-tools:0.2.0"
TOOL_VERSIONS = {
    "fastgear": {"fastGEAR_sha256": "c" * 64},
    "envs": {"bakta": {"lock": "bakta.lock", "tools": {
        "bakta": {"version": "1.10.4", "build": "pyhdfd78af_1", "channel": "bioconda"}}}},
}
INSPECT = {
    "Id": IMAGE_ID,
    "RepoDigests": [f"ghcr.io/example/bactolomew-tools@{DIGEST}"],
    "RepoTags": [REF],
    "Created": "2026-10-13T10:00:00Z",
    "Config": {"Labels": {"org.opencontainers.image.version": "0.2.0",
                          "org.opencontainers.image.revision": "abc123",
                          "maintainer": "ignored"}},
}


IMAGE_CONTENTS = {"mobsuite_db": {"download_date": "2026-10-13", "taxa_sqlite_sha256": "d" * 64}}


def _image_read_stdout(tool_versions, contents, status_line) -> str:
    """What the image read script prints: the two JSON files and mob_init's status line."""
    sec = CR._SECTION
    return (f"{json.dumps(tool_versions) if tool_versions else ''}\n{sec}\n"
            f"{json.dumps(contents) if contents else ''}\n{sec}\n{status_line}\n")


def _fake_run(cmd, **_kwargs):
    if cmd[:3] == ["docker", "image", "inspect"]:
        return SimpleNamespace(returncode=0, stdout=json.dumps(INSPECT) + "\n", stderr="")
    if cmd[:2] == ["docker", "run"] and cmd[-1] == CR._IMAGE_READ_SCRIPT:
        stdout = _image_read_stdout(TOOL_VERSIONS, IMAGE_CONTENTS, "Download date: 2026-10-13. Removing lock file.")
        return SimpleNamespace(returncode=0, stdout=stdout, stderr="")
    raise AssertionError(f"unexpected command: {cmd}")


class RunProvenanceTest(unittest.TestCase):
    def test_environment_block_in_run_config(self) -> None:
        ctx = ContainerContext(runtime="docker", image=f"{REF}@{DIGEST}")
        with mock.patch.object(subprocess, "run", side_effect=_fake_run), \
                mock.patch.dict("os.environ", {}, clear=False) as env:
            env.pop("BACTOLOMEW_NATIVE", None)
            environment = CR.environment_provenance(ctx, configured_image=f"{REF}@{DIGEST}")

        self.assertEqual(environment["runtime"], "docker")
        self.assertEqual(environment["image"], f"{REF}@{DIGEST}")
        self.assertEqual(environment["image_id"], IMAGE_ID)
        self.assertEqual(environment["repo_digests"], INSPECT["RepoDigests"])
        self.assertEqual(environment["image_labels"], {"version": "0.2.0", "revision": "abc123"})
        self.assertEqual(environment["pinned_digest"], DIGEST)
        self.assertTrue(environment["pinned_digest_matches"])
        self.assertEqual(environment["tool_versions"]["envs"], {"bakta": {"bakta": "1.10.4 pyhdfd78af_1"}})
        self.assertEqual(environment["tool_versions"]["fastgear_sha256"], "c" * 64)
        self.assertNotIn("BACTOLOMEW_NATIVE", environment)
        mobsuite = environment["databases"]["mobsuite"]
        self.assertEqual((mobsuite["location"], mobsuite["download_date"]), ("image", "2026-10-13"))
        self.assertEqual(mobsuite["taxa_sqlite_sha256"], "d" * 64)

        cfg = SimpleNamespace(output=Path("/tmp/out"))
        with mock.patch.object(Runner, "_build_run_inputs", return_value={}):
            run_config = Runner._build_run_config(
                cfg, sample_set_id="s", run_id="r", annotation_id_strategy="a",
                thread_plan={}, environment=environment,
            )
        self.assertEqual(json.loads(json.dumps(run_config["environment"])), environment)

    def test_older_image_without_records(self) -> None:
        """0.1.0 has neither tool_versions.json nor image_contents.json: the mob_suite
        date still comes from mob_init's status line."""
        def fake(cmd, **_kwargs):
            if cmd[:3] == ["docker", "image", "inspect"]:
                return SimpleNamespace(returncode=0, stdout=json.dumps(INSPECT) + "\n", stderr="")
            return SimpleNamespace(returncode=0, stdout=_image_read_stdout(
                None, None, "Download date: 2026-07-08. Removing lock file."), stderr="")

        with mock.patch.object(subprocess, "run", side_effect=fake):
            environment = CR.environment_provenance(ContainerContext(runtime="docker", image="bactolomew-tools:0.1.0"))
        self.assertIsNone(environment["tool_versions"])
        self.assertIn("tool_versions_source", environment)
        self.assertEqual(environment["databases"]["mobsuite"]["download_date"], "2026-07-08")

    def test_database_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            bakta = Path(tmp) / "bakta" / "db"
            (bakta / "amrfinderplus-db" / "latest").mkdir(parents=True)
            (bakta / "version.json").write_text(json.dumps(
                {"date": "2024-07-01", "major": 6, "minor": 0, "type": "light", "doi": "10.5281/x"}))
            (bakta / "amrfinderplus-db" / "latest" / "version.txt").write_text("2024-07-22.1\n")
            checkm2 = Path(tmp) / "uniref100.KO.1.dmnd"
            checkm2.write_bytes(b"x" * 10)
            mob = Path(tmp) / "mob"
            mob.mkdir()
            (mob / "status.txt").write_text("Download date: 2026-01-02. Removing lock file.")

            dbs = CR.database_provenance(bakta_db=bakta.parent, checkm2_db=checkm2, mobsuite_db_dir=mob)

        self.assertEqual(dbs["bakta"]["type"], "light")
        self.assertEqual((dbs["bakta"]["major"], dbs["bakta"]["minor"]), (6, 0))
        self.assertEqual(dbs["bakta"]["amrfinderplus_db"], "2024-07-22.1")
        self.assertEqual(dbs["checkm2"]["size"], 10)
        self.assertEqual(dbs["mobsuite"], {"location": "external", "path": str(mob), "download_date": "2026-01-02"})
        self.assertEqual(CR.database_provenance(bakta_db="/nonexistent")["bakta"]["version"], None)

    def test_environment_history_keeps_earlier_environments(self) -> None:
        env_a = {"runtime": "docker", "image": "bactolomew-tools:0.1.0", "image_id": "sha256:1"}
        env_b = {"runtime": "docker", "image": "bactolomew-tools:0.2.0", "image_id": "sha256:2"}
        # Same environment again: nothing to keep.
        self.assertEqual(Runner._environment_history({"created_at": "t0", "environment": env_a}, env_a), [])
        # Resumed under another image: the earlier block is kept, once.
        first = Runner._environment_history({"created_at": "t0", "environment": env_a}, env_b)
        self.assertEqual(first, [{"recorded_at": "t0", "environment": env_a}])
        again = Runner._environment_history(
            {"created_at": "t1", "environment": env_b, "environment_history": first}, env_b)
        self.assertEqual(again, first)
        back = Runner._environment_history(
            {"created_at": "t1", "environment": env_b, "environment_history": first}, env_a)
        self.assertEqual(back, [*first, {"recorded_at": "t1", "environment": env_b}])
        # A folder from before provenance existed.
        old = Runner._environment_history({"created_at": "t0"}, env_b)
        self.assertEqual(len(old), 1)
        self.assertIsNone(old[0]["environment"])

    def test_sif_identity_hash_is_cached(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sif = Path(tmp) / "tools.sif"
            sif.write_bytes(b"sif content")
            with mock.patch.dict(os.environ, {"XDG_CACHE_HOME": str(Path(tmp) / "cache")}):
                first = CR.sif_sha256(str(sif))
                with mock.patch.object(CR, "_sha256_file", side_effect=AssertionError("not cached")):
                    self.assertEqual(CR.sif_sha256(str(sif)), first)
                sif.write_bytes(b"other content!")
                self.assertNotEqual(CR.sif_sha256(str(sif)), first)
        self.assertEqual(first, hashlib.sha256(b"sif content").hexdigest())


class ImageAcquisitionTest(unittest.TestCase):
    """container.build_on_missing and which references are pulled (no Docker started)."""

    PINNED = f"{REF}@{DIGEST}"

    def test_repository_path_and_pull_ref(self) -> None:
        self.assertEqual(CR.pull_ref(self.PINNED), self.PINNED)
        self.assertEqual(CR.pull_ref("localhost:5000/tools:1"), "localhost:5000/tools:1")
        self.assertEqual(CR.pull_ref("bactolomew-tools:0.2.0"), "")
        self.assertEqual(CR.pull_ref("bactolomew-tools:0.2.0", "ghcr.io/x/"), "ghcr.io/x/bactolomew-tools:0.2.0")
        self.assertEqual(configured_image_ref({"container": {"image": "localhost:5000/tools:1", "registry": "r.io"}}),
                         "localhost:5000/tools:1")
        self.assertEqual(configured_image_ref({"container": {"image": "tools:1", "registry": "r.io"}}), "r.io/tools:1")

    def _ensure(self, image, *, pull_ok, build_on_missing, local_name_present=False):
        calls = []

        def present(ref, engine="docker"):
            return local_name_present and ref == CR.split_image_digest(image)[0]

        def pull(img, registry="", log=None, engine="docker"):
            calls.append(("pull", img))
            if not pull_ok:
                raise RuntimeError("docker pull failed")

        def build(img, **_kw):
            calls.append(("build", img))

        with mock.patch.object(CR, "docker_image_present", side_effect=present), \
                mock.patch.object(CR, "pull_image", side_effect=pull), \
                mock.patch.object(CR, "build_image_local", side_effect=build), \
                mock.patch.object(CR, "pinned_digest_mismatch", return_value=""):
            CR.ensure_docker_image(image, log=lambda _m: None, build_on_missing=build_on_missing)
        return calls

    def test_pinned_pull_failure_is_an_error_without_build_on_missing(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "build_on_missing is false"):
            self._ensure(self.PINNED, pull_ok=False, build_on_missing=False)
        self.assertEqual(self._ensure(self.PINNED, pull_ok=False, build_on_missing=True),
                         [("pull", self.PINNED), ("build", REF)])
        self.assertEqual(self._ensure(self.PINNED, pull_ok=True, build_on_missing=False), [("pull", self.PINNED)])

    def test_bare_name_builds_only_with_build_on_missing(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "names no registry"):
            self._ensure("bactolomew-tools:0.2.0", pull_ok=True, build_on_missing=False)
        self.assertEqual(self._ensure("bactolomew-tools:0.2.0", pull_ok=True, build_on_missing=True),
                         [("build", "bactolomew-tools:0.2.0")])

    def test_singularity_pulls_a_full_reference_without_registry(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sif = str(Path(tmp) / "tools.sif")
            streamed = []

            def stream(cmd, log, **_kw):
                streamed.append(cmd)
                return 0

            with mock.patch.object(CR, "_stream", side_effect=stream), \
                    mock.patch.object(CR, "singularity_bin", return_value="apptainer"), \
                    mock.patch.object(CR, "ensure_docker_image", side_effect=AssertionError("needs Docker")):
                CR.ensure_singularity_image(sif, self.PINNED, registry="", log=lambda _m: None)
            self.assertEqual(streamed, [["apptainer", "pull", sif, f"docker://{self.PINNED}"]])

            def failing(cmd, log, **_kw):
                return 1

            with mock.patch.object(CR, "_stream", side_effect=failing), \
                    mock.patch.object(CR, "singularity_bin", return_value="apptainer"), \
                    mock.patch.object(CR, "docker_status",
                                      return_value=CR.RuntimeStatus("docker", False, "not_running", "not running")):
                with self.assertRaisesRegex(RuntimeError, "build_on_missing is false"):
                    CR.ensure_singularity_image(sif, self.PINNED, log=lambda _m: None, build_on_missing=False)
                with self.assertRaisesRegex(RuntimeError, "needs a Docker daemon"):
                    CR.ensure_singularity_image(sif, self.PINNED, log=lambda _m: None, build_on_missing=True)


if __name__ == "__main__":
    unittest.main()
