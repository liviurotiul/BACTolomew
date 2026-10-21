"""cgMLST schema import (pipeline.CgmlstInstall): sources are detected and registered
without touching the other catalog entries, and an installed schema is replaced only
on request, the old one kept under .replaced/.

The chewBBACA build (PrepExternalSchema and the smoke AlleleCall in the tools image)
is mocked; the real build is covered by test_pipeline_cgmlst_schema.
"""
from __future__ import annotations

import io
import json
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock

import pipeline.CgmlstInstall as ci
from pipeline.ChewbbacaSchemaBuild import (
    BUILD_INFO_FILE,
    CHEWBBACA_SCHEME_CATALOG,
    REPLACED_DIR,
    SCHEMA_BUILD_REPORT_DIR,
    SchemaBuildResult,
    _stream_zip_members,
)

CDS = "ATGAAACCCGGGTTTAAACCCGGGAAATTTGCTTAA"
CATALOG = {
    "version": 1,
    "schemes": [
        {"schema_id": "acinetobacter_baumannii", "schema_dir": "acinetobacter_baumannii",
         "species": ["Acinetobacter baumannii"], "mlst_schemes": ["abaumannii", "abaumannii_2"]},
        {"schema_id": "bacillus_anthracis", "schema_dir": "bacillus_anthracis",
         "species": ["Bacillus anthracis"], "mlst_schemes": ["banthracis"]},
    ],
}


def _fake_build(source, built, **_kw):
    """Stands in for PrepExternalSchema + the smoke AlleleCall: a minimal schema folder."""
    built = Path(built)
    (built / "short").mkdir(parents=True)
    (built / ".genes_list").write_bytes(b"")
    (built / ".schema_config").write_bytes(b"")
    (built / "L1.fasta").write_text(f">L1_1\n{CDS}\n")
    (built / "short" / "L1_short.fasta").write_text(f">L1_1\n{CDS}\n")
    (built / "marker.txt").write_text(f"built from {source}\n")
    (built / SCHEMA_BUILD_REPORT_DIR).mkdir()
    (built / SCHEMA_BUILD_REPORT_DIR / BUILD_INFO_FILE).write_text(json.dumps({"source": str(source)}))
    return SchemaBuildResult(built, built.name, 1, 1, 0, "test", "", alleles=1, representatives=1)


class CgmlstImportTests(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        root = Path(self._tmp.name)
        self.dest = root / "cgmlst"
        self.dest.mkdir()
        (self.dest / CHEWBBACA_SCHEME_CATALOG).write_text(json.dumps(CATALOG, indent=2))
        # Installed A. baumannii schema (it has grown with inferred alleles since).
        self.installed = self.dest / "acinetobacter_baumannii"
        (self.installed / "short").mkdir(parents=True)
        (self.installed / ".genes_list").write_bytes(b"")
        (self.installed / "marker.txt").write_text("installed schema with novel alleles\n")
        # What the user downloaded by hand: a cgMLST.org zip, a folder of loci, an unnamed zip.
        self.downloads = root / "downloads"
        self.downloads.mkdir()
        with zipfile.ZipFile(self.downloads / "Banthracis_cgMLST_alleles.zip", "w") as zf:
            zf.writestr("BAS0001.fasta", f">1\n{CDS}\n")
        (self.downloads / "Abaumannii").mkdir()
        (self.downloads / "Abaumannii" / "ACICU_RS00005.fasta").write_text(f">1\n{CDS}\n")
        with zipfile.ZipFile(self.downloads / "export.zip", "w") as zf:
            zf.writestr("XYZ0001.fasta", f">1\n{CDS}\n")

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def _import(self, paths, **kw):
        with mock.patch.object(ci, "_build_and_validate", side_effect=_fake_build):
            results = ci.import_paths(paths, dest=self.dest, update_config=False, cpu=1, index=[], **kw)
        return {Path(r.source).name: r for r in results}

    def test_import_detects_species_and_keeps_other_entries(self) -> None:
        results = self._import([self.downloads])

        imported = results["Banthracis_cgMLST_alleles.zip"]
        self.assertEqual((imported.status, imported.schema_id), ("imported", "bacillus_anthracis"))
        self.assertTrue(imported.registered)
        self.assertTrue((self.dest / "bacillus_anthracis" / ".genes_list").is_file())
        self.assertEqual(results["Abaumannii"].status, "skipped")
        self.assertIn("already installed", results["Abaumannii"].reason)
        self.assertEqual(results["export.zip"].status, "failed")
        self.assertIn("species not detected", results["export.zip"].reason)

        # The catalog entries are those of before, the imported one with its MLST aliases.
        catalog = json.loads((self.dest / CHEWBBACA_SCHEME_CATALOG).read_text())
        self.assertEqual(catalog["schemes"], CATALOG["schemes"])
        self.assertEqual((self.installed / "marker.txt").read_text(), "installed schema with novel alleles\n")
        self.assertFalse((self.dest / ci.STAGING_DIR).exists())

        # A species given by the user for an unnamed source makes a new catalog entry.
        results = self._import([self.downloads / "export.zip"], species="Xenorhabdus testii", mlst_schemes=["xeno"])
        self.assertEqual((results["export.zip"].status, results["export.zip"].schema_id),
                         ("imported", "xenorhabdus_testii"))
        catalog = json.loads((self.dest / CHEWBBACA_SCHEME_CATALOG).read_text())
        self.assertEqual(catalog["schemes"][:2], CATALOG["schemes"])
        self.assertEqual(catalog["schemes"][2]["species"], ["Xenorhabdus testii"])
        self.assertEqual(catalog["schemes"][2]["mlst_schemes"], ["xeno"])

    def test_replace_keeps_the_old_schema(self) -> None:
        source = self.downloads / "Abaumannii"
        self.assertEqual(self._import([source])["Abaumannii"].status, "skipped")

        result = self._import([source], replace=True)["Abaumannii"]
        self.assertEqual(result.status, "imported")
        backup = Path(result.backup)
        self.assertEqual(backup.parent, self.dest / REPLACED_DIR)
        self.assertTrue(backup.name.startswith("acinetobacter_baumannii."))
        self.assertEqual((backup / "marker.txt").read_text(), "installed schema with novel alleles\n")
        self.assertIn("built from", (self.installed / "marker.txt").read_text())
        catalog = json.loads((self.dest / CHEWBBACA_SCHEME_CATALOG).read_text())
        self.assertEqual(catalog["schemes"], CATALOG["schemes"])
        # The kept copy is not an installed schema.
        self.assertEqual([s.schema_id for s in ci.installed(self.dest)],
                         ["acinetobacter_baumannii"])


class ChewbbacaStepInstallTests(unittest.TestCase):
    """The chewBBACA step's schema: a missing cgMLST.org scheme is installed inside the step
    (its progress reported as the step's status); a species without one fails the step."""

    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        root = Path(self._tmp.name)
        self.dest = root / "cgmlst"
        self.dest.mkdir()
        (self.dest / CHEWBBACA_SCHEME_CATALOG).write_text(json.dumps(CATALOG, indent=2))
        ci.accept_terms(self.dest)
        self.run_dir = root / "run"
        (self.run_dir / "mlst").mkdir(parents=True)
        self.config = {"catalog_path": str(self.dest / CHEWBBACA_SCHEME_CATALOG), "auto_select_from_mlst": True}
        patches = [
            mock.patch("pipeline.ToolEnv._load_app_config", return_value={}),  # not the user's config
            mock.patch("pipeline.ChewBBACAWorkflow.validate_chewbbaca_settings", side_effect=lambda s: s),
            mock.patch.object(ci, "_build_and_validate", side_effect=_fake_build),
        ]
        for patch in patches:
            patch.start()
            self.addCleanup(patch.stop)

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def _mlst(self, scheme: str) -> None:
        rows = [{"sample": f"g{i}", "scheme": scheme, "st": "1"} for i in range(3)]
        (self.run_dir / "mlst" / "results.json").write_text(json.dumps(rows))

    def test_missing_scheme_installs_inside_the_step_and_no_scheme_fails_it(self) -> None:
        import logging

        from pipeline.ChewBBACAWorkflow import CgmlstSchemaMissing, prepare_chewbbaca_schema

        index = [ci.OrgScheme("Banthracis", "Bacillus anthracis cgMLST", "bacillus_anthracis",
                              ["Bacillus anthracis"], "https://www.cgmlst.org/ncs/schema/Banthracis/", "", "")]

        def fake_download(url, out_path, *, on_bytes=None, **_kw):
            Path(out_path).parent.mkdir(parents=True, exist_ok=True)
            Path(out_path).write_bytes(b"x" * 400)
            for received in (100, 200, 400):
                if on_bytes is not None:
                    on_bytes(received, 400)
            return Path(out_path)

        import threading

        from pipeline.RunSteps import StepProgress

        # The Runner's status callback: the step's run_state.json "progress" as it is written.
        steps = {"chewbbaca": {"status": "running"}}
        written: list[object] = []
        store = StepProgress(steps, threading.RLock(),
                             lambda: written.append(dict(steps["chewbbaca"].get("progress") or {})), interval=0)
        runner_status = store.step_status("chewbbaca")
        statuses: list[tuple[str, object]] = []

        def status(text, share):
            statuses.append((text, share))
            runner_status(text, share)

        logger = logging.getLogger("test.chewbbaca.install")
        self._mlst("banthracis")
        with mock.patch.object(ci, "fetch_org_index", return_value=index), \
                mock.patch.object(ci, "_resolve_download_links", return_value=("csv", "zip")), \
                mock.patch.object(ci, "download_file", side_effect=fake_download):
            settings = prepare_chewbbaca_schema(self.run_dir, logger, config=self.config, status=status)
        self.assertEqual(Path(settings["schema_dir"]), (self.dest / "bacillus_anthracis").resolve())
        self.assertTrue((self.dest / "bacillus_anthracis" / ".genes_list").is_file())
        shares = [share for text, share in statuses if text]
        self.assertTrue(all(text == "Installing cgMLST scheme · B. anthracis" for text, _ in statuses[:-1]))
        self.assertEqual(shares, sorted(shares))
        self.assertGreater(len(set(shares)), 4)  # download bytes and build phases, not only 0 and 1
        self.assertEqual(shares[-1], 1.0)
        self.assertEqual(statuses[-1], ("", None))  # the step's status is cleared after the install
        self.assertEqual(written[0], {"label": "Installing cgMLST scheme · B. anthracis", "unit": "%",
                                      "done": 0, "total": 100})
        self.assertEqual(written[-2]["done"], 100)
        self.assertEqual(written[-1], {})  # removed from run_state.json at the end
        self.assertNotIn("progress", steps["chewbbaca"])

        # An MLST scheme no catalog entry (and so no cgMLST.org scheme) maps to: the step fails.
        self._mlst("vparahaemolyticus")
        with self.assertRaises(CgmlstSchemaMissing) as caught:
            prepare_chewbbaca_schema(self.run_dir, logger, config=self.config)
        self.assertTrue(str(caught.exception).startswith(
            "No cgMLST scheme for the MLST scheme vparahaemolyticus: import one in Settings > Add schema"))
        self.assertIn("--mlst vparahaemolyticus", str(caught.exception))


class _Unseekable(io.RawIOBase):
    """A write-only stream: zipfile writes a data descriptor after each member, as
    cgMLST.org's streamed downloads have them."""

    def __init__(self) -> None:
        super().__init__()
        self.data = bytearray()

    def writable(self) -> bool:
        return True

    def write(self, b) -> int:  # noqa: ANN001
        self.data += bytes(b)
        return len(b)


class ZipStreamTests(unittest.TestCase):
    def test_truncated_or_damaged_zip_fails_with_and_without_data_descriptors(self) -> None:
        loci = {f"L{i}.fasta": f">1\n{CDS * (i + 3)}\n".encode() for i in range(4)}
        ordinary, streamed = io.BytesIO(), _Unseekable()
        for sink in (ordinary, streamed):
            with zipfile.ZipFile(sink, "w", zipfile.ZIP_DEFLATED) as zf:
                for name, body in loci.items():
                    zf.writestr(name, body)
        with tempfile.TemporaryDirectory() as tmp:
            for kind, data in (("no descriptors", ordinary.getvalue()), ("descriptors", bytes(streamed.data))):
                members = [i for i in range(len(data)) if data.startswith(b"PK\x03\x04", i)]
                directory = data.find(b"PK\x01\x02")
                damaged = bytearray(data)
                damaged[members[1] + 60] ^= 0xFF
                cases = {
                    "cut in the directory": data[:directory + 20],
                    "cut between members": data[:members[2]],
                    "cut inside a member": data[:members[2] + 45],
                    "damaged member": bytes(damaged),
                }
                for case, cut in cases.items():
                    with self.subTest(kind=kind, case=case):
                        path = Path(tmp) / "alleles.zip"
                        path.write_bytes(cut)
                        if case != "damaged member":  # zipfile cannot open these: the stream reader does
                            self.assertFalse(zipfile.is_zipfile(path))
                        read: dict[str, bytes] = {}
                        if case == "cut in the directory":
                            _stream_zip_members(path, read.__setitem__)
                            self.assertEqual(read, loci)
                        else:
                            with self.assertRaisesRegex(ValueError, "truncated or damaged"):
                                _stream_zip_members(path, read.__setitem__)


if __name__ == "__main__":
    unittest.main()
