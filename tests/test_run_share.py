"""Shared runs (pipeline.RunShare): a package round trip opens in another runs
folder with its own tables, the compatibility check reports and never blocks,
unsafe package names are refused, and a moved run is found where it is now.

Tiny synthetic runs in a temporary folder (a few JSON files, a tree and three
small parquet tables); never a real run.
"""
from __future__ import annotations

import json
import tempfile
import unittest
import zipfile
from pathlib import Path

import pandas as pd

from annotation_explorer_polars import cache
from annotation_explorer_polars.run_expectations import load_run_status, run_allows_view
from pipeline import RunShare
from pipeline.RunValidation import RESULT_REQUIRED_STEPS


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def make_run(root: Path, *, run_id: str = "sample_abc__run_1", recorded: Path | None = None) -> Path:
    """A completed run folder at ``root`` that recorded ``recorded`` as its output_dir."""
    recorded = recorded or root
    steps = {name: {"status": "success", "outputs_validated": True} for name in RESULT_REQUIRED_STEPS}
    for name in ("assembly_qc", "bakta", "mlst", "mobsuite", "panaroo", "gubbins", "phylogeny",
                 "fastgear_clades", "chewbbaca", "diamond_markers", "diamond_reference_loci", "collect_annotations"):
        steps.setdefault(name, {"status": "success", "outputs_validated": True})
    inputs = {"genomes": ["/genomes/A.fna", "/genomes/B.fna"], "marker_sequences": [], "reference_locus_sequences": []}
    _write_json(root / "run_config.json", {
        "run_id": run_id, "sample_set_id": "sample_abc", "output_dir": str(recorded), "inputs": inputs,
    })
    _write_json(root / "run_report.json", {
        "status": "success", "outputs_validated": True, "run_id": run_id, "steps": steps,
        "output_dir": str(recorded),
    })
    _write_json(root / "run_state.json", {"status": "success", "outputs_validated": True, "run_id": run_id, "steps": steps})
    (root / "trees").mkdir(parents=True, exist_ok=True)
    (root / "trees" / "core.nwk").write_text("(A:0.1,B:0.2);\n")
    _write_json(root / "trees" / "trees.json", {"schema_version": 2, "trees": []})
    (root / "mlst").mkdir(exist_ok=True)
    (root / "mlst" / "results.csv").write_text("sample,scheme,sequence_type\nA,abaumannii_2,2\nB,abaumannii_2,1\n")
    # Working data a package leaves out.
    (root / "bakta" / "A").mkdir(parents=True, exist_ok=True)
    (root / "bakta" / "A" / "A.json").write_text("{}")
    (root / "bakta" / "A" / "A.fna").write_text(">c1\nACGT\n")
    for rel in ("chewbbaca/out/results_alleles.tsv", "chewbbaca/out/results_statistics.tsv",
                "chewbbaca/out/loci_summary_stats.tsv", "chewbbaca/out/logging_info.txt", "mlst/results.json",
                "mobsuite/contigs.tsv", "mobsuite/plasmids.tsv", "mobsuite/summary.json", "assembly_qc/summary.tsv",
                "assembly_qc/summary.json", "assembly_qc/quast/report.tsv", "assembly_qc/checkm2/quality_report.tsv"):
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        (root / rel).write_text("")
    # expected_outputs.json lists absolute paths under the recorded folder.
    _write_json(root / "expected_outputs.json", {
        "output_dir": str(recorded),
        "files": [str(recorded / "run_config.json"), str(recorded / "trees" / "core.nwk")],
    })
    return root


def make_tables(folder: Path, schema: int = cache.TABLE_CACHE_SCHEMA_VERSION) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"Sample": ["A", "B"], "Gene": ["g1", "g2"]}).to_parquet(folder / cache.ANNOTATIONS_CACHE_FILENAME)
    pd.DataFrame({"Sample": ["A"]}).to_parquet(folder / cache.REFERENCE_LOCUS_CACHE_FILENAME)
    pd.DataFrame({"Sample": ["B"]}).to_parquet(folder / cache.MARKER_HITS_CACHE_FILENAME)
    _write_json(folder / cache.TABLE_CACHE_META_FILENAME, {"schema_version": schema, "output_dir": "/home/someone/run"})
    return folder


class RunShareTests(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(prefix="bactolomew_share_")
        self.tmp = Path(self._tmp.name)

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def test_package_round_trip_opens_in_another_runs_folder(self) -> None:
        run = make_run(self.tmp / "made_here" / "run")
        tables = make_tables(self.tmp / "cache_of_maker")
        dest = RunShare.export_run(
            run, tables, self.tmp / "out" / "run.bactolomew.zip",
            include_assemblies=False, scheme_species={"abaumannii_2": "Acinetobacter baumannii"}, app_version="0.2.0",
        )
        names = zipfile.ZipFile(dest).namelist()
        self.assertTrue(all(name.startswith("sample_abc__run_1/") for name in names))
        self.assertNotIn("sample_abc__run_1/bakta/A/A.json", names)  # build inputs stay home
        self.assertNotIn("sample_abc__run_1/bakta/A/A.fna", names)  # assemblies only when asked
        manifest = RunShare.read_manifest(dest)
        self.assertEqual(manifest["format_version"], RunShare.SHARE_FORMAT_VERSION)
        self.assertNotIn(str(self.tmp), json.dumps(manifest))  # no paths of the maker's machine

        runs_root = self.tmp / "viewer" / "runs"
        result = RunShare.import_archive(dest, runs_root)
        self.assertEqual(result.path.parent.parent, runs_root / RunShare.SHARED_RUNS_DIR)
        self.assertEqual((result.damaged, result.missing, result.renamed_from), ((), (), None))
        self.assertTrue(RunShare.is_shared_run(result.path) and RunShare.is_imported_run(result.path))
        self.assertTrue(run_allows_view(load_run_status(result.path)))
        paths = cache.table_paths(result.run_id, result.path)
        self.assertEqual(paths[0].parent, result.path / RunShare.TABLES_DIR)
        self.assertEqual(cache.shared_tables_status(result.path), (True, cache.TABLE_CACHE_SCHEMA_VERSION))
        ann_df, _ref, _hits = cache.read_cached_tables(result.run_id, paths)
        self.assertEqual(list(ann_df["Gene"]), ["g1", "g2"])
        from pipeline.RunMetadata import mlst_scheme_species

        self.assertEqual(mlst_scheme_species(result.path)["abaumannii_2"], "Acinetobacter baumannii")

        # The same run again: its run id is taken in this runs folder.
        again = RunShare.import_archive(dest, runs_root)
        self.assertEqual(again.renamed_from, "sample_abc__run_1")
        self.assertNotEqual(again.run_id, result.run_id)
        self.assertEqual(json.loads((again.path / "run_config.json").read_text())["run_id"], again.run_id)

    def test_package_keeps_no_paths_of_the_exporting_machine(self) -> None:
        import getpass

        run = make_run(self.tmp / "made_here" / "run", recorded=self.tmp / "made_first" / "run")
        home, root = str(Path.home()), str(run)
        config = json.loads((run / "run_config.json").read_text())
        config["inputs"].update(genomes=[f"{home}/genomes/A.fna", f"{home}/genomes/B.fna"],
                                snippy_reference=f"{home}/genomes/A.fna", bakta_db="/Shared/bakta/db")
        config["environment"] = {"databases": {"bakta": {"path": f"{home}/db/bakta"}}}
        _write_json(run / "run_config.json", config)
        (run / "mlst" / "results.csv").write_text(
            f"sample,input_file,scheme,sequence_type\nA,{home}/genomes/A.fna,abaumannii_2,2\n"
        )
        (run / "chewbbaca" / "out" / "logging_info.txt").write_text(
            f"Command: chewBBACA.py -o {root}/chewbbaca/out --schema {self.tmp}/made_first/run/schema\n"
        )
        _write_json(run / "trees" / "trees.json", {
            "schema_version": 2, "trees": [{"id": "core", "path": f"{root}/trees/core.nwk"}],
        })
        dest = RunShare.export_run(run, make_tables(self.tmp / "tables"), self.tmp / "out.zip", app_version="0.2.0")

        user = getpass.getuser()
        with zipfile.ZipFile(dest) as archive:
            for name in archive.namelist():
                text = archive.read(name).decode("utf-8", "replace")
                for private in (home, root, str(self.tmp)) + ((user,) if len(user) >= 3 else ()):
                    self.assertNotIn(private, text, f"{name} holds {private}")
        result = RunShare.import_archive(dest, self.tmp / "viewer" / "runs")
        moved = json.loads((result.path / "run_config.json").read_text())
        self.assertEqual(moved["output_dir"], ".")
        self.assertEqual(moved["inputs"]["snippy_reference"], "A.fna")  # the reference stem still reads
        trees = json.loads((result.path / "trees" / "trees.json").read_text())
        self.assertEqual(trees["trees"][0]["path"], "trees/core.nwk")
        self.assertTrue(run_allows_view(load_run_status(result.path)))

    def test_newer_package_format_imports_with_a_warning(self) -> None:
        run = make_run(self.tmp / "made_here" / "run")
        dest = RunShare.export_run(run, make_tables(self.tmp / "tables"), self.tmp / "out.zip", app_version="")
        newer = self.tmp / "newer.zip"
        with zipfile.ZipFile(dest) as src, zipfile.ZipFile(newer, "w") as out:
            for info in src.infolist():
                data = src.read(info)
                if info.filename.endswith(RunShare.MANIFEST):
                    data = json.dumps({**json.loads(data), "format_version": RunShare.SHARE_FORMAT_VERSION + 1}).encode()
                out.writestr(info, data)
        runs_root = self.tmp / "viewer" / "runs"
        result = RunShare.import_archive(newer, runs_root)
        self.assertTrue(RunShare.run_id_taken(runs_root, RunShare.import_run_id(newer)))
        (result.path / "mlst" / "results.csv").unlink()
        issues = RunShare.check_compatibility(result.path, app_version="0.2.0")
        codes = {issue.code: issue for issue in issues}
        self.assertEqual(codes["format_newer"].level, "warning")
        self.assertIn("1 file missing: mlst/results.csv.", codes["missing"].text)
        self.assertEqual(codes["missing"].short, "1 missing file")
        self.assertIn("version_unknown", codes)  # the run recorded none; the exporter's is not used
        self.assertNotIn("app_older", codes)

    def test_compatibility_warns_and_notes_without_blocking(self) -> None:
        manifest = {
            "format_version": 1, "app_version": "0.3.1", "made_with": "0.3.1", "status": "success",
            "schemas": {"tables": cache.TABLE_CACHE_SCHEMA_VERSION + 1, "trees": 2},
            "files": [{"path": "share/tables/annotations.parquet", "size": 1}], "environment": None,
        }
        issues = RunShare.check_compatibility(
            manifest, app_version="0.2.0", table_schema=cache.TABLE_CACHE_SCHEMA_VERSION, trees_schema=2,
        )
        codes = {issue.code: issue.level for issue in issues}
        self.assertEqual(codes.get("app_newer"), "warning")
        self.assertEqual(codes.get("table_schema"), "warning")
        self.assertEqual(codes.get("no_environment"), "info")
        self.assertEqual(issues[0].level, "warning")  # warnings first
        self.assertTrue(RunShare.warning_summary(issues).startswith("2 warnings: made with BACTolomew 0.3.1"))

    def test_image_note_compares_the_configured_image(self) -> None:
        from unittest import mock

        manifest = {"made_with": "0.2.0", "status": "success", "environment": {"image": "bactolomew-tools:0.1.0"}}

        def note(config: dict) -> list[str]:
            with mock.patch("pipeline.ToolEnv._load_app_config", return_value=config):
                context = RunShare.current_versions()
            context.pop("image_tool_versions")
            issues = RunShare.check_compatibility(manifest, **context)
            return [issue.text for issue in issues if issue.code == "image_differs"]

        docker = {"runtime": "docker", "container": {"image": "bactolomew-tools:0.1.0"}}
        self.assertEqual(note(docker), [])
        other = {"runtime": "docker", "container": {"image": "ghcr.io/x/bactolomew-tools:9.9.9"}}
        self.assertEqual(note(other), ["Tools image bactolomew-tools:0.1.0; this install uses ghcr.io/x/bactolomew-tools:9.9.9."])
        self.assertEqual(note({**other, "mode": "viewer"}), [])  # a viewer has no image

    def test_unsafe_member_names_are_refused(self) -> None:
        bad = self.tmp / "bad.zip"
        with zipfile.ZipFile(bad, "w") as archive:
            archive.writestr("run/share/manifest.json", json.dumps({"format_version": 1, "files": []}))
            archive.writestr("run/../../escape.txt", "x")
        runs_root = self.tmp / "runs"
        with self.assertRaises(RunShare.ShareError):
            RunShare.import_archive(bad, runs_root)
        self.assertFalse((self.tmp / "escape.txt").exists())
        self.assertEqual(list(runs_root.glob(RunShare.STAGING_PREFIX + "*")) if runs_root.exists() else [], [])

    def test_a_moved_run_is_viewable_where_it_is_now(self) -> None:
        moved = make_run(self.tmp / "now_here" / "run", recorded=Path("/somewhere/else/run"))
        status = load_run_status(moved)
        self.assertTrue(run_allows_view(status), status.get("validation_errors"))


if __name__ == "__main__":
    unittest.main()
