"""The headless path: every pipeline module and the entry point import without PyQt6 and
matplotlib, and the core and RAM budget follows the cgroup (v2, then v1) and SLURM limits."""
from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

from pipeline.Resources import detect_budget

SRC = Path(__file__).resolve().parents[1] / "src"

_IMPORT_ALL_BLOCKED = textwrap.dedent("""
    import importlib, pkgutil, sys

    class Block:
        def find_spec(self, name, path=None, target=None):
            if name.split(".")[0] in ("PyQt6", "PyQt6Ads", "matplotlib", "pyqtgraph"):
                raise ImportError("blocked " + name)

    sys.meta_path.insert(0, Block())
    import pipeline
    for mod in pkgutil.iter_modules(pipeline.__path__):
        importlib.import_module("pipeline." + mod.name)
    import bactolomew_main
    bactolomew_main.build_parser().parse_args(["run", "--out", "x", "--resume"])
    assert bactolomew_main.main(["--version"]) == 0
""")


class HeadlessImportTest(unittest.TestCase):
    def test_pipeline_imports_without_qt_and_matplotlib(self):
        res = subprocess.run([sys.executable, "-c", _IMPORT_ALL_BLOCKED], cwd=str(SRC),
                             capture_output=True, text=True, timeout=300)
        self.assertEqual(res.returncode, 0, res.stderr[-2000:])
        self.assertIn("BACTolomew", res.stdout)


class ResourceBudgetTest(unittest.TestCase):
    def _tree(self, files: dict[str, str]) -> Path:
        root = Path(tempfile.mkdtemp(prefix="cgroup-"))
        self.addCleanup(shutil.rmtree, root, True)
        for rel, text in files.items():
            path = root / rel
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text)
        return root

    def test_cgroup_v2_limits_from_the_process_cgroup_and_its_parents(self):
        root = self._tree({
            "job/cpu.max": "max 100000\n",
            "job/memory.max": str(6 << 30),
            "job/step/cpu.max": "250000 100000\n",  # 2.5 CPUs -> 3 cores
            "job/step/memory.max": "max\n",
        })
        budget = detect_budget({}, cgroup_root=root, proc_cgroup="0::/job/step\n")
        self.assertEqual(budget.core_limits["cgroup CPU quota"], 3)
        self.assertLessEqual(budget.cores, 3)
        self.assertAlmostEqual(budget.memory_limits["cgroup memory limit"], 6.0)
        self.assertLessEqual(budget.memory_gb, 6.0)

    def test_cgroup_v1_and_slurm(self):
        root = self._tree({
            "cpu,cpuacct/slurm/cpu.cfs_quota_us": "400000",
            "cpu,cpuacct/slurm/cpu.cfs_period_us": "100000",
            "memory/slurm/memory.limit_in_bytes": str(9223372036854771712),  # no limit
        })
        proc = "4:memory:/slurm\n3:cpu,cpuacct:/slurm\n"
        env = {"SLURM_CPUS_PER_TASK": "2", "SLURM_MEM_PER_CPU": "3G"}
        budget = detect_budget(env, cgroup_root=root, proc_cgroup=proc)
        self.assertEqual(budget.core_limits["cgroup CPU quota"], 4)
        self.assertEqual(budget.cores, min(2, *budget.core_limits.values()))
        self.assertNotIn("cgroup memory limit", budget.memory_limits)
        self.assertAlmostEqual(budget.memory_limits["SLURM_MEM_PER_CPU"], 6.0)
        self.assertIn("SLURM_CPUS_PER_TASK", budget.describe())


class SetupPlanTest(unittest.TestCase):
    def test_setup_plans_a_sif_and_names_a_runtime_that_cannot_work(self):
        import argparse
        from unittest import mock

        import bactolomew_main as M
        from pipeline import ContainerRuntime as CR

        root = Path("/opt/bactolomew")
        args = argparse.Namespace(runtime="apptainer", sif=None)
        missing = CR.RuntimeStatus("singularity", False, "not_installed", "Apptainer is not installed.")
        with mock.patch.object(CR, "runtime_status", return_value=missing):
            _config, updates, problem = M._setup_config(args, {}, root)
        self.assertEqual(updates, {"runtime": "singularity",
                                   "container": {"sif_path": "/opt/bactolomew/bactolomew-tools-0.2.0.sif"}})
        self.assertIn("not installed", problem)
        # No runtime anywhere: still a plan (for docker), with the reason it cannot run.
        args = argparse.Namespace(runtime=None, sif=None)
        with mock.patch.object(CR, "suggest_runtime", return_value=None):
            config, _updates, problem = M._setup_config(args, {}, root)
        self.assertEqual((config["runtime"], problem), ("docker", M.NO_RUNTIME))


if __name__ == "__main__":
    unittest.main()
