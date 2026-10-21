"""Container engines: the user flags each engine needs so bind-mount writes belong to the
host user, Apptainer's --home, and telling Docker's failures apart. No container starts."""
from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pipeline import ContainerRuntime as CR  # noqa: E402
from pipeline.ToolEnv import ContainerContext, container_image_cmd  # noqa: E402


class EngineArgvTest(unittest.TestCase):
    def test_user_flags_image_and_home_per_engine(self) -> None:
        base = dict(image="owner/tools:1", bind_roots=("/data",), uid=1000, gid=1000, home_dir="/data/home")
        engines = {
            "rootful docker": (ContainerContext(runtime="docker", **base), ["--user", "1000:1000"], "owner/tools:1"),
            "rootless docker": (ContainerContext(runtime="docker", rootless=True, **base), [], "owner/tools:1"),
            "rootless podman, SELinux": (
                ContainerContext(runtime="podman", engine_exe="podman", podman=True, rootless=True, selinux=True, **base),
                ["--userns=keep-id", "--security-opt", "label=disable"], "docker.io/owner/tools:1"),
            "rootful podman": (ContainerContext(runtime="podman", engine_exe="podman", podman=True, **base),
                               ["--user", "1000:1000"], "docker.io/owner/tools:1"),
        }
        for name, (ctx, flags, image) in engines.items():
            with self.subTest(name):
                cmd = container_image_cmd(ctx, ["run.sh"], workdir="/data", extra_env={"HOME": "/data/c"})
                self.assertEqual(cmd[:3 + len(flags)], [ctx.engine_exe or "docker", "run", "--rm", *flags])
                self.assertEqual(cmd[-2:], [image, "run.sh"])
        sing = ContainerContext(runtime="singularity", sif_path="/s.sif", engine_exe="apptainer", **base)
        cmd = container_image_cmd(sing, ["run.sh"], extra_env={"HOME": "/data/c", "MCR_CACHE_ROOT": "/data/c"})
        self.assertEqual(cmd[:5], ["apptainer", "exec", "--cleanenv", "--home", "/data/c"])
        self.assertNotIn("HOME=/data/c", cmd)
        self.assertIn("MCR_CACHE_ROOT=/data/c", cmd)


class DockerStatusTest(unittest.TestCase):
    def test_failures_are_told_apart(self) -> None:
        cases = {
            "permission": "permission denied while trying to connect to the docker API at unix:///var/run/docker.sock",
            "not_running": "Cannot connect to the Docker daemon at unix:///var/run/docker.sock. Is the docker daemon running?",
            "not_usable": "error during connect: something else",
        }
        for problem, stderr in cases.items():
            with self.subTest(problem), mock.patch.object(CR.shutil, "which", return_value="/usr/bin/docker"), \
                    mock.patch.object(CR, "engine_info", return_value=CR.EngineInfo("docker")), \
                    mock.patch.object(CR.subprocess, "run",
                                      return_value=subprocess.CompletedProcess([], 1, "", stderr)):
                status = CR.docker_status()
                self.assertEqual((status.ok, status.problem), (False, problem))
                self.assertIn(stderr, status.details)
        with mock.patch.object(CR.shutil, "which", return_value=None):
            self.assertEqual(CR.docker_status().problem, "not_installed")


class EngineProbeTest(unittest.TestCase):
    def setUp(self) -> None:
        CR._ENGINE_INFO.clear()
        self.addCleanup(CR._ENGINE_INFO.clear)

    def test_failed_probe_is_not_kept(self) -> None:
        """Docker down at the first probe, rootless once it runs: the second probe sees it."""
        daemon = {"up": False}

        def probe(cmd):
            if cmd[1] == "--version":
                return 0, "Docker version 27.0"
            return (0, '["name=seccomp","name=rootless"]') if daemon["up"] else (1, "Cannot connect")

        with mock.patch.object(CR.shutil, "which", return_value="/usr/bin/docker"), \
                mock.patch.object(CR, "_probe", side_effect=probe), \
                mock.patch.object(CR, "selinux_enforcing", return_value=False):
            self.assertFalse(CR.engine_info("docker").rootless)
            daemon["up"] = True
            self.assertTrue(CR.engine_info("docker").rootless)
            daemon["up"] = False
            self.assertTrue(CR.engine_info("docker").rootless)  # a good probe is kept

    def test_native_context_uses_the_probed_engine(self) -> None:
        ctx = ContainerContext(runtime="native", image="owner/tools:1", uid=1000, gid=1000)
        podman = CR.EngineInfo("podman", podman=True, rootless=True)
        with mock.patch.object(CR, "engine_info", return_value=podman):
            cmd = container_image_cmd(ctx, ["run.sh"], workdir="/data")
        self.assertEqual(cmd[:4], ["podman", "run", "--rm", "--userns=keep-id"])
        self.assertEqual(cmd[-2:], ["docker.io/owner/tools:1", "run.sh"])

    def test_fastgear_base_by_name_must_be_the_pinned_image(self) -> None:
        digest = next(iter(CR.KNOWN_IMAGE_IDS))
        pinned = f"pgcbioinfo/fastgear:latest@{digest}"
        loaded = {"id": CR.KNOWN_IMAGE_IDS[digest], "repo_digests": []}  # a `docker load` copy
        other = {"id": "sha256:" + "0" * 64, "repo_digests": []}
        for info, ok in ((loaded, True), (other, False)):
            with self.subTest(ok=ok), \
                    mock.patch.object(CR, "docker_image_present", side_effect=lambda ref, _e: "@" not in ref), \
                    mock.patch.object(CR, "docker_image_info", return_value=info):
                if ok:
                    self.assertEqual(CR.ensure_fastgear_image(pinned, log=lambda _l: None),
                                     "pgcbioinfo/fastgear:latest")
                else:
                    with mock.patch.object(CR, "_stream", return_value=1):
                        with self.assertRaises(RuntimeError):
                            CR.ensure_fastgear_image(pinned, log=lambda _l: None)


if __name__ == "__main__":
    unittest.main()
