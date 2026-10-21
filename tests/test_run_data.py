"""RunData and QueryEngine (window architecture plan 4.3, Phase 2).

The run's data is built once (``RunData.build``, in the loader thread in the
app) and the filter queries run on ``QueryEngine`` without the window. The
window keeps read-only views of the data under the old attribute names, so the
mixins are unchanged. tests/test_refresh_parity.py records the results through
the window; here the same recording is checked against a bare engine.
"""
from __future__ import annotations

import concurrent.futures
import dataclasses
import inspect
import json
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from pathlib import Path
from typing import Any, Dict, List, Optional
from unittest import mock

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
    from tests.test_refresh_parity import GOLDEN_PATH, build_case_matrix, build_synthetic_run, case_filters, record_all
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings
    from test_refresh_parity import GOLDEN_PATH, build_case_matrix, build_synthetic_run, case_filters, record_all

isolate_qt_settings()

import pandas as pd
import polars as pl
from PyQt6 import QtCore, QtWidgets

from annotation_explorer_polars.model import QueryEngine, RunData
from annotation_explorer_polars.constants import MISSING_FILTER_TOKEN
from annotation_explorer_polars.model import run_data as run_data_mod
from annotation_explorer_polars.model.filter_state import AppliedResult, FilterState, SampleSubset


# Attribute the window exposes -> RunData field it reads.
WINDOW_VIEWS: Dict[str, str] = {
    "ann_df": "ann_df",
    "meta_df": "meta_df",
    "ref_df": "ref_df",
    "marker_hits_df": "marker_hits_df",
    "reference_hits_df": "reference_hits_df",
    "_ann_pl": "ann_pl",
    "_meta_pl": "meta_pl",
    "_contig_totals_pl": "contig_totals_pl",
    "annotation_hit_counts_pl": "annotation_hit_counts_pl",
    "_mobsuite_plasmid_members": "mobsuite_plasmid_members",
    "_plasmid_cluster_index": "plasmid_cluster_index",
    "_contig_total_count_map": "contig_total_count_map",
    "sample_total_counts": "sample_total_counts",
    "sample_atb_counts": "sample_atb_counts",
    "contig_atb_counts": "contig_atb_counts",
    "_ann_public_columns": "ann_public_columns",
    "_meta_public_columns": "meta_public_columns",
    "meta_columns": "meta_columns",
    "_sample_qc_columns": "sample_qc_columns",
    "annotation_db_columns": "annotation_db_columns",
    "predictor_options": "predictor_options",
    "_configured_marker_dbs": "configured_marker_dbs",
    "_configured_reference_dbs": "configured_reference_dbs",
    "_snippy_reference_sample": "snippy_reference_sample",
    "_marker_db_category_map": "marker_db_category_map",
    "_neo_tree_has_resistance_traits": "has_resistance_traits",
    "_neo_tree_has_mge_traits": "has_mge_traits",
    "_neo_tree_has_virulence_traits": "has_virulence_traits",
    "output_dir": "output_dir",
    "run_id": "run_id",
}


def build_marker_and_reference_hits() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Marker hits for three trait sources (one row without a DB is dropped) and one reference hit."""
    marker_hits = pd.DataFrame(
        [
            {
                "AnnotationID": "S1:1:0",
                "Sample": "S1",
                "Contig": "1",
                "DB": "AMRFinderPlus_core",
                "Gene": "blaOXA-23",
                "Product": "class D beta-lactamase OXA-23",
                "ResistanceDrugClasses": "BETA-LACTAM; carbapenem",
            },
            {
                "AnnotationID": "S1:3:0",
                "Sample": "S1",
                "Contig": "3",
                "DB": "ISFinder-Prot_MOBILE",
                "Gene": "tnpA",
                "Product": "IS4 family // ISAba1~~~transposase",
                "ResistanceDrugClasses": "",
            },
            {
                "AnnotationID": "S2:2:0",
                "Sample": "S2",
                "Contig": "2",
                "DB": "VFDB-setB",
                "Gene": "ompA",
                "Product": "outer membrane protein A",
                "ResistanceDrugClasses": "",
            },
            {
                "AnnotationID": "S2:1:0",
                "Sample": "S2",
                "Contig": "1",
                "DB": " ",
                "Gene": "lysR",
                "Product": "LysR family transcriptional regulator",
                "ResistanceDrugClasses": "aminoglycoside",
            },
        ]
    )
    ref_hits = pd.DataFrame([{"Sample": "S3", "Contig": "2", "DB": "PAIDB", "Start": 1, "Stop": 900}])
    return marker_hits, ref_hits


def build_parity_run(output_dir: Path) -> RunData:
    """The run the golden recording was made on (no marker or reference hits)."""
    ann_df, meta_df = build_synthetic_run()
    return RunData.build(
        ann_df, meta_df, ref_df=pd.DataFrame(), marker_hits_df=pd.DataFrame(), output_dir=output_dir
    )


def _filters(**changes: Any) -> FilterState:
    """A FilterState from the former snapshot's pieces (no view mode since D3)."""
    return FilterState.from_query(**changes)


class _TempDirTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.output_dir = Path(tempfile.mkdtemp(prefix="bactolomew_run_data_"))
        self.addCleanup(shutil.rmtree, self.output_dir, True)


class RunDataBuildTests(_TempDirTestCase):
    def _build(self) -> RunData:
        ann_df, meta_df = build_synthetic_run()
        marker_hits, ref_hits = build_marker_and_reference_hits()
        return RunData.build(
            ann_df, meta_df, ref_df=ref_hits, marker_hits_df=marker_hits, output_dir=self.output_dir, run_id="run_x"
        )

    def test_build_derives_the_run_tables_and_counts(self) -> None:
        run = self._build()
        self.assertEqual(run.run_id, "run_x")
        self.assertEqual(run.output_dir, self.output_dir)
        # The marker hit without a DB is dropped; every text column exists.
        self.assertEqual(run.marker_hits_df["AnnotationID"].tolist(), ["S1:1:0", "S1:3:0", "S2:2:0"])
        for column in ("SubjectIdentifier", "ResistanceMechanisms", "Method"):
            self.assertIn(column, run.marker_hits_df.columns)
        # Marker search text is attached to the annotation the hit belongs to.
        oxa = run.ann_df.loc[run.ann_df["AnnotationID"] == "S1:1:0", "MarkerSearchText"].iloc[0]
        self.assertIn("amrfinderplus_core", str(oxa).lower())
        self.assertIn("MarkerSearchText", run.ann_public_columns)
        self.assertEqual(run.reference_hits_df.values.tolist(), [["S3", "2", "PAIDB"]])
        self.assertEqual(
            run.annotation_db_columns,
            ("AMRFinderPlus_core", "ISFinder-Prot_MOBILE", "PAIDB", "VFDB-setB"),
        )
        self.assertEqual(
            run.annotation_hit_counts_pl.rows(),
            [("S1", "1", "AMRFinderPlus_core", 1), ("S1", "3", "ISFinder-Prot_MOBILE", 1),
             ("S2", "2", "VFDB-setB", 1), ("S3", "2", "PAIDB", 1)],
        )
        self.assertEqual(dict(run.sample_atb_counts), {"S1": 2})
        self.assertEqual(dict(run.contig_atb_counts), {("S1", "1"): 2})
        self.assertEqual((run.has_resistance_traits, run.has_mge_traits, run.has_virulence_traits), (True, True, True))
        # Counts over the annotations and contigs.
        self.assertEqual(dict(run.sample_total_counts), {"S1": 7, "S2": 6, "S3": 5, "S4": 5, "S5": 4, "S6": 4})
        self.assertEqual(run.contig_total_count_map[("S1", "1")], 3)
        self.assertEqual(run.contig_totals_pl.height, 18)
        self.assertEqual(dict(run.mobsuite_plasmid_members), {("S1", "AA001"): ("3",), ("S2", "AA002"): ("3",), ("S4", "AA001"): ("3",)})
        # Metadata columns: Sample and QC columns are not tree metadata columns.
        self.assertEqual(run.meta_columns, ("Country", "Year", "Source", "Latitude", "Longitude"))
        self.assertEqual(run.sample_qc_columns, ())
        self.assertIn("__search_text", run.meta_pl.columns)
        self.assertIn("__gene_norm", run.ann_pl.columns)
        self.assertIn("CARD_AMR", run.predictor_options)
        # No run_config.json: no configured databases and no snippy reference.
        self.assertEqual((run.configured_marker_dbs, run.configured_reference_dbs), ((), ()))
        self.assertIsNone(run.snippy_reference_sample)

    def test_run_config_gives_configured_databases_and_reference_sample(self) -> None:
        (self.output_dir / "run_config.json").write_text(
            json.dumps(
                {
                    "inputs": {
                        "marker_sequences": ["/db/ZZZ_markers.fasta", "/db/ZZZ_markers.fasta"],
                        "reference_locus_sequences": ["/db/PAIDB_islands.fasta"],
                        "snippy_reference": "/refs/S4.gbk",
                    }
                }
            ),
            encoding="utf-8",
        )
        run = self._build()
        self.assertEqual(run.configured_marker_dbs, ("ZZZ_markers",))
        self.assertEqual(run.configured_reference_dbs, ("PAIDB_islands",))
        self.assertEqual(run.snippy_reference_sample, "S4")
        # Configured databases come first; the reference hit "PAIDB" maps to the one configured "PAIDB_*".
        self.assertEqual(run.reference_hits_df["DB"].tolist(), ["PAIDB_islands"])
        self.assertEqual(
            run.annotation_db_columns,
            ("ZZZ_markers", "PAIDB_islands", "AMRFinderPlus_core", "ISFinder-Prot_MOBILE", "VFDB-setB"),
        )

    def test_run_data_is_read_only(self) -> None:
        run = self._build()
        with self.assertRaises(dataclasses.FrozenInstanceError):
            run.ann_df = pd.DataFrame()  # type: ignore[misc]
        for name in ("sample_total_counts", "sample_atb_counts", "contig_atb_counts", "contig_total_count_map",
                     "mobsuite_plasmid_members", "marker_db_category_map"):
            with self.subTest(field=name), self.assertRaises(TypeError):
                getattr(run, name)["new"] = 1
        for name in ("annotation_db_columns", "predictor_options", "meta_columns", "ann_public_columns"):
            with self.subTest(field=name):
                self.assertIsInstance(getattr(run, name), tuple)

    def test_marker_hits_are_required(self) -> None:
        ann_df, meta_df = build_synthetic_run()
        with self.assertRaises(ValueError):
            RunData.build(ann_df, meta_df, marker_hits_df=None)

    def test_build_reports_progress(self) -> None:
        ann_df, meta_df = build_synthetic_run()
        seen: List[str] = []
        RunData.build(ann_df, meta_df, marker_hits_df=pd.DataFrame(), status=seen.append)
        self.assertEqual(seen, ["Preparing annotations…", "Indexing contigs and plasmids…", "Counting database hits…"])

    def test_build_in_a_worker_thread_matches_the_main_thread(self) -> None:
        main = self._build()
        box: Dict[str, Any] = {}

        def _work() -> None:
            try:
                box["run"] = self._build()
            except BaseException as exc:  # reported below
                box["error"] = exc

        thread = threading.Thread(target=_work)
        thread.start()
        thread.join(60)
        self.assertNotIn("error", box)
        other = box["run"]
        self.assertTrue(main.ann_pl.equals(other.ann_pl))
        self.assertTrue(main.meta_pl.equals(other.meta_pl))
        self.assertTrue(main.ann_df.equals(other.ann_df))
        self.assertTrue(
            main.contig_totals_pl.sort(["Sample", "Contig"]).equals(other.contig_totals_pl.sort(["Sample", "Contig"]))
        )
        for name in ("sample_total_counts", "sample_atb_counts", "contig_atb_counts", "contig_total_count_map"):
            self.assertEqual(dict(getattr(main, name)), dict(getattr(other, name)), name)
        self.assertEqual(main.annotation_db_columns, other.annotation_db_columns)


    def test_plasmid_cluster_index_is_frozen(self) -> None:
        index = self._build().plasmid_cluster_index
        self.assertTrue(index)
        self.assertEqual(sorted(index.presence_by_sample), ["S1", "S2", "S4"])
        for name in ("clusters", "presence_by_sample", "instance_to_cluster"):
            with self.subTest(field=name), self.assertRaises(TypeError):
                getattr(index, name)["new"] = 1
        self.assertIsInstance(index.presence_by_sample["S1"], frozenset)
        with self.assertRaises(RuntimeError):
            index.add_row({"Sample": "S9", "MobSuiteIsPlasmid": "True", "MobSuitePlasmidID": "AA009"})
        self.assertNotIn("S9", index.presence_by_sample)
        # Reads still work on the frozen index.
        key = index.instance_to_cluster[("S1", "AA001")]
        self.assertTrue(index.display_name_for_cluster(key))
        self.assertIn(key, index.cluster_order())

    def test_contig_order_and_plasmid_members_are_deterministic(self) -> None:
        # Phase 3: contigs keep their annotation order (polars group_by used to
        # hand them out in a different order per run, and with them the member
        # order of a multi-contig plasmid, its first table row and its track).
        ann_pl = pl.DataFrame(
            {
                "Sample": ["S1"] * 6 + ["S2"] * 2,
                "Contig": ["7", "7", "2", "12", "2", "9", "1", "1"],
                "MobSuiteIsPlasmid": ["True", "True", "True", "True", "True", "False", "False", "False"],
                "MobSuitePlasmidID": ["AB1", "AB1", "AB1", "AB1", "AB1", "", "", ""],
            }
        )
        for _round in range(5):
            totals = run_data_mod.build_contig_totals(ann_pl)
            self.assertEqual(
                totals.select(["Sample", "Contig"]).rows(),
                [("S1", "7"), ("S1", "2"), ("S1", "12"), ("S1", "9"), ("S2", "1")],
            )
            self.assertEqual(totals.get_column("n_total").to_list(), [2, 2, 1, 1, 2])
            members = run_data_mod.build_mobsuite_plasmid_members(totals)
            self.assertEqual(members, {("S1", "AB1"): ("7", "2", "12")})

    def test_model_package_loads_without_qt(self) -> None:
        # The loader thread builds RunData; importing it must not pull in Qt.
        # Since Phase 3 the query engine and the filter types load no Qt
        # either, nor (Phase 9) the run's PanData.
        code = (
            "import sys\n"
            "from annotation_explorer_polars.model.run_data import RunData\n"
            "from annotation_explorer_polars.model import QueryEngine, FilterState, AppliedResult, PanData\n"
            "QueryEngine, FilterState, AppliedResult, PanData\n"
            "print(sorted(m for m in sys.modules if m.startswith('PyQt6')))\n"
        )
        env = dict(os.environ)
        src = str(Path(__file__).resolve().parents[1] / "src")
        env["PYTHONPATH"] = src + os.pathsep + env.get("PYTHONPATH", "")
        out = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True, timeout=120)
        self.assertEqual(out.returncode, 0, out.stderr)
        self.assertEqual(out.stdout.strip(), "[]")


class DatabaseClassifierTests(unittest.TestCase):
    """The classifiers and trait detection moved from the tree and filter mixins unchanged."""

    def test_resistance_databases(self) -> None:
        categories = {"MyDB": "AMR", "OtherDB": "bmrg", "Plain": "VIRULENCE"}
        cases = {
            "AMRFinderPlus_core": True,
            "card_variants": True,  # CARD_ prefix, any case
            "XCARD_db": False,  # CARD_ only as a prefix
            "BacMet_BMRG": True,
            "MyDB": True,  # by its category only
            "OtherDB": True,
            "Plain": False,
            "VFDB-setB": False,
            "": False,
            "  ": False,
        }
        for name, expected in cases.items():
            with self.subTest(db=name):
                self.assertIs(run_data_mod.db_is_resistance(name, categories), expected)

    def test_mobile_element_databases(self) -> None:
        categories = {"MobileMGE": "MOBILE", "MobileOther": "MOBILE", "MGE_amr": "AMR"}
        cases = {
            "ISFinder-Prot_MOBILE": True,
            "cge-mge_db": True,
            "MobileMGE": True,  # MOBILE category and "MGE" in the name
            "MobileOther": False,  # MOBILE category alone is not enough
            "MGE_amr": False,  # "MGE" in the name without the MOBILE category
            "VFDB-setB": False,
            "": False,
        }
        for name, expected in cases.items():
            with self.subTest(db=name):
                self.assertIs(run_data_mod.db_is_mge_source(name, categories), expected)

    def test_virulence_databases(self) -> None:
        categories = {"MyVF": "virulence", "Other": "AMR"}
        cases = {"VFDB-setB": True, "virulence_finder": True, "MyVF": True, "Other": False, "ISFinder": False, "": False}
        for name, expected in cases.items():
            with self.subTest(db=name):
                self.assertIs(run_data_mod.db_is_virulence_source(name, categories), expected)

    def test_trait_detection(self) -> None:
        empty = pd.DataFrame()
        for detect in (run_data_mod.detect_resistance_traits, run_data_mod.detect_mge_traits, run_data_mod.detect_virulence_traits):
            with self.subTest(detect=detect.__name__):
                self.assertFalse(detect(empty, {}))
                self.assertFalse(detect(None, {}))  # type: ignore[arg-type]
        cog_only = pd.DataFrame([{"DB": "COG", "Gene": "lysR", "Product": "IS4 family transposase", "ResistanceDrugClasses": ""}])
        self.assertFalse(run_data_mod.detect_resistance_traits(cog_only, {}))
        # A transposase named in a non-mobile database is no mobile element trait.
        self.assertFalse(run_data_mod.detect_mge_traits(cog_only, {}))
        self.assertFalse(run_data_mod.detect_virulence_traits(cog_only, {}))
        # Resistance from a drug class column alone, whatever the database.
        drug_class = cog_only.assign(ResistanceDrugClasses="aminoglycoside")
        self.assertTrue(run_data_mod.detect_resistance_traits(drug_class, {}))
        # A resistance database by its category.
        self.assertTrue(run_data_mod.detect_resistance_traits(cog_only, {"COG": "AMR"}))
        # A mobile element database needs a label in its rows.
        isfinder = pd.DataFrame([{"DB": "ISFinder-Prot_MOBILE", "Gene": "", "Product": ""}])
        self.assertFalse(run_data_mod.detect_mge_traits(isfinder, {}))
        self.assertTrue(run_data_mod.detect_mge_traits(isfinder.assign(Product="IS4 family // ISAba1~~~transposase"), {}))
        # A virulence database by its category, with a gene label; a generic label is no trait.
        vf = pd.DataFrame([{"DB": "MyVF", "Gene": "ompA", "Product": ""}])
        self.assertTrue(run_data_mod.detect_virulence_traits(vf, {"MyVF": "VIRULENCE"}))
        self.assertFalse(run_data_mod.detect_virulence_traits(vf, {}))
        generic = pd.DataFrame([{"DB": "VFDB", "Gene": "hypothetical protein", "GeneID": "", "Product": "hypothetical protein"}])
        self.assertFalse(run_data_mod.detect_virulence_traits(generic, {}))
        # Without a DB column every row is read for labels (a mobile element label needs a DB name).
        no_db = pd.DataFrame([{"Gene": "ompA", "Product": "IS26 transposase"}])
        self.assertTrue(run_data_mod.detect_virulence_traits(no_db, {}))
        self.assertFalse(run_data_mod.detect_mge_traits(no_db, {}))

    def test_build_without_trait_sources_has_no_traits(self) -> None:
        ann_df, meta_df = build_synthetic_run()
        cog_only = pd.DataFrame([{"AnnotationID": "S1:1:1", "Sample": "S1", "Contig": "1", "DB": "COGorg24", "Gene": "lysR",
                                  "Product": "LysR family transcriptional regulator"}])
        with mock.patch.object(run_data_mod, "load_marker_db_category_map", return_value={}):
            run = RunData.build(ann_df, meta_df, marker_hits_df=cog_only)
        self.assertEqual((run.has_resistance_traits, run.has_mge_traits, run.has_virulence_traits), (False, False, False))

    def test_mobsuite_plasmid_group_key(self) -> None:
        key = run_data_mod.mobsuite_plasmid_group_key
        self.assertEqual(key({"MobSuiteIsPlasmid": "False", "MobSuitePlasmidID": "AA001"}), "")
        self.assertEqual(key({"MobSuiteIsPlasmid": None, "MobSuitePlasmidID": "AA001"}), "")
        self.assertEqual(key({"MobSuiteIsPlasmid": " Yes ", "MobSuitePlasmidID": "AA001", "Sample": "S1"}), "AA001")
        self.assertEqual(
            key({"MobSuiteIsPlasmid": "plasmid", "MobSuitePlasmidID": "nan", "Sample": "S1", "MobSuitePrimaryClusterID": "AA100"}),
            "S1:AA100",
        )
        self.assertEqual(
            key({"MobSuiteIsPlasmid": "1", "MobSuitePlasmidID": "", "Sample": "S1", "MobSuitePrimaryClusterID": "<NA>",
                 "MobSuiteSecondaryClusterID": "AB200"}),
            "S1:secondary:AB200",
        )
        # Clusters without a sample give no key.
        self.assertEqual(key({"MobSuiteIsPlasmid": "true", "MobSuitePrimaryClusterID": "AA100"}), "")
        self.assertEqual(key({"MobSuiteIsPlasmid": "true", "Sample": "S1"}), "")


class QueryEngineTests(_TempDirTestCase):
    def setUp(self) -> None:
        super().setUp()
        self.run = build_parity_run(self.output_dir)
        self.engine = QueryEngine(self.run)

    def test_engine_reproduces_the_golden_recording_without_a_window(self) -> None:
        golden = json.loads(GOLDEN_PATH.read_text(encoding="utf-8"))
        recorded = record_all(self.engine)
        self.assertEqual(sorted(recorded), sorted(golden))
        for case_id, expected in golden.items():
            with self.subTest(case=case_id):
                self.assertEqual(recorded[case_id], expected)

    def test_result_cache_shares_one_immutable_result(self) -> None:
        calls: List[FilterState] = []
        compute = self.engine._compute

        def counting(filters: FilterState) -> Any:
            calls.append(filters)
            return compute(filters)

        self.engine._compute = counting  # type: ignore[method-assign]
        ann = {"gene_text": "blaOXA", "product_text": None, "molecule": "", "predictors": ["CARD_AMR", "COGorg24_COG"], "tags": [], "hide_hypothetical": False}
        first = self.engine.compute(_filters(ann_filters=dict(ann), meta_filters={"Country": ["Romania", "Italy"]}))
        # Same filters, keys and values in another order: the same cached object.
        reordered = dict(reversed(list(ann.items())), predictors=["COGorg24_COG", "CARD_AMR"])
        second = self.engine.compute(_filters(ann_filters=reordered, meta_filters={"Country": ["Italy", "Romania"]}))
        self.assertEqual(len(calls), 1)
        self.assertIs(first, second)
        self.assertIsInstance(first, AppliedResult)
        # Nothing a reader holds can change the cached result.
        self.assertIsInstance(first.samples, tuple)
        self.assertIsInstance(first.allowed_samples, tuple)
        with self.assertRaises(dataclasses.FrozenInstanceError):
            first.samples = ()  # type: ignore[misc]
        with self.assertRaises(TypeError):
            first.filtered_counts["S1"] = 0  # type: ignore[index]
        with self.assertRaises(TypeError):
            first.contig_records["S1"] = ()  # type: ignore[index]
        self.assertIsInstance(first.contig_records["S1"], tuple)
        # Any other field is another entry; a subset label is not a field of the query.
        self.engine.compute(_filters(ann_filters=dict(ann), hide_empty_contigs=True))
        self.engine.compute(_filters(ann_filters=dict(ann), sample_query="s1"))
        self.assertEqual(len(calls), 3)
        labelled = FilterState(subset=SampleSubset(("S1", "S2"), label="Clade A"))
        relabelled = FilterState(subset=SampleSubset(("S1", "S2"), label="Clade B"))
        self.assertIs(self.engine.compute(labelled), self.engine.compute(relabelled))
        self.assertEqual(len(calls), 4)

    def test_result_cache_is_bounded(self) -> None:
        for token, query in enumerate(["s1", "s2", "s3", "s4", "s5", "s6", "romania", "italy", "spain", "2019"], start=1):
            self.engine.compute(_filters(sample_query=query))
        self.assertEqual(len(self.engine._result_cache), 8)

    def test_compute_is_thread_safe(self) -> None:
        snapshots = [case_filters(case) for case in build_case_matrix()[::9]]

        def summary(result: Any) -> tuple:
            return (
                result.allowed_samples,
                result.samples,
                result.filtered_ann.get_column("AnnotationID").to_list(),
                result.stats_text,
            )

        expected = [summary(QueryEngine(self.run).compute(s)) for s in snapshots]
        shared = QueryEngine(self.run)
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
            for _round in range(3):
                got = list(pool.map(lambda s: summary(shared.compute(s)), snapshots))
                self.assertEqual(got, expected)

    def test_contig_annotations_reads_all_or_filtered_rows(self) -> None:
        all_s1 = self.engine.contig_annotations(allowed_samples=["S1", "S2"], sample="S1", contig="1")
        self.assertEqual(all_s1["AnnotationID"].tolist(), ["S1:1:0", "S1:1:1", "S1:1:2"])
        self.assertEqual(list(all_s1.columns), list(self.run.ann_public_columns))
        merged = self.engine.contig_annotations(allowed_samples=["S1"], sample="S1", contigs=["1", "3"], contig="2")
        self.assertEqual(sorted(merged["Contig"].unique()), ["1", "3"])
        self.assertTrue(self.engine.contig_annotations(allowed_samples=[], sample="S1", contig="1").empty)
        self.assertTrue(self.engine.contig_annotations(allowed_samples=["S2"], sample="S1", contig="1").empty)
        result = self.engine.compute(_filters(ann_filters={"gene_text": "blaOXA"}))
        hits = self.engine.contig_annotations(filtered=result.filtered_ann, sample="S1", contig="3")
        self.assertEqual(hits["AnnotationID"].tolist(), ["S1:3:1"])
        empty = self.engine.contig_annotations(filtered=self.engine.empty_annotations(), sample="S1", contig="3")
        self.assertTrue(empty.empty)
        self.assertEqual(list(empty.columns), list(self.run.ann_public_columns))

    def test_map_payload_for_a_sample_set(self) -> None:
        self.assertEqual(self.engine.map_payload([]), {"points": [], "countries": {}, "maxCountryCount": 1})
        payload = self.engine.map_payload(["S1", "S2", "S5", "S9"])
        self.assertEqual([p["sample"] for p in payload["points"]], ["S1"])
        self.assertEqual(payload["countries"], {"Italy": 1})
        self.assertEqual(payload["total"], 4)

    def test_card_bundles_share_the_run_frames(self) -> None:
        bundle = self.engine.sample_bundle("S1")
        self.assertEqual(bundle["sample"], "S1")
        self.assertIs(bundle["ann_df"], self.run.ann_df)
        self.assertIs(bundle["reference_hits_df"], self.run.reference_hits_df)
        gene = self.engine.annotation_bundle()
        self.assertIs(gene["marker_hits_df"], self.run.marker_hits_df)
        self.assertEqual(set(gene), {"ann_df", "meta_df", "marker_hits_df"})


    def test_result_cache_evicts_the_least_recently_used(self) -> None:
        computes: List[Optional[str]] = []
        compute = self.engine._compute

        def counting(filters: FilterState) -> Any:
            computes.append(filters.sample_query)
            return compute(filters)

        self.engine._compute = counting  # type: ignore[method-assign]
        queries = [f"s{index}" for index in range(1, 9)]
        for token, query in enumerate(queries, start=1):
            self.engine.compute(_filters(sample_query=query))
        # Use the oldest entry again, then add two more: s2 and s3 go, s1 stays.
        self.engine.compute(_filters(sample_query="s1"))
        self.engine.compute(_filters(sample_query="romania"))
        self.engine.compute(_filters(sample_query="italy"))
        self.assertEqual(len(self.engine._result_cache), 8)
        computes.clear()
        self.engine.compute(_filters(sample_query="s1"))
        self.engine.compute(_filters(sample_query="s4"))
        self.assertEqual(computes, [])
        self.engine.compute(_filters(sample_query="s2"))
        self.assertEqual(computes, ["s2"])

    def test_contig_record_rows_are_read_only(self) -> None:
        first = self.engine.compute(_filters())
        row = first.contig_records["S1"][0]
        with self.assertRaises(TypeError):
            row["n_hits"] = 99  # type: ignore[index]
        # A copy is what a reader changes (as the sample tree does).
        copy = dict(row)
        copy["__display_label"] = "Contig 1"
        second = self.engine.compute(_filters())
        self.assertIs(second.contig_records["S1"][0], row)
        self.assertNotIn("__display_label", second.contig_records["S1"][0])

    def test_metadata_subset(self) -> None:
        def samples(query: Optional[str], filters: Dict[str, List[str]]) -> List[str]:
            return self.engine.metadata_subset(query, filters).get_column("Sample").to_list()

        self.assertEqual(samples(None, {}), ["S1", "S2", "S3", "S4", "S5", "S6"])
        self.assertEqual(samples("  ROMANIA ", {}), ["S1", "S4"])
        self.assertEqual(samples(None, {"Country": [MISSING_FILTER_TOKEN]}), ["S5"])
        self.assertEqual(samples(None, {"Country": ["Italy", MISSING_FILTER_TOKEN]}), ["S2", "S5", "S6"])
        self.assertEqual(samples("2020", {"Source": ["blood"]}), ["S3"])
        # Unknown columns and empty value lists do not filter.
        self.assertEqual(samples(None, {"Nope": ["x"], "Country": []}), ["S1", "S2", "S3", "S4", "S5", "S6"])
        subset = self.engine.metadata_subset(None, {"Year": ["2019"]})
        self.assertEqual(subset.columns, list(self.run.meta_public_columns))
        # The value order does not matter for the cache.
        self.assertIs(
            self.engine.metadata_subset(None, {"Country": ["Italy", "Spain"]}),
            self.engine.metadata_subset(None, {"Country": ["Spain", "Italy"]}),
        )

    def test_filter_annotations(self) -> None:
        def ids(filters: Dict[str, Any], samples: List[str], **kwargs: Any) -> List[str]:
            return self.engine.filter_annotations(filters, samples, **kwargs).get_column("AnnotationID").to_list()

        self.assertEqual(ids({"gene_text": "blaOXA"}, ["S1", "S3"]), ["S1:1:0", "S1:3:1", "S3:2:0"])
        self.assertEqual(
            ids({"molecule": "plasmid", "hide_hypothetical": True}, ["S1", "S2"]),
            ["S1:3:0", "S1:3:1", "S2:3:0", "S2:3:1"],
        )
        self.assertEqual(ids({"tags": ["v"]}, ["S1", "S2", "S9"]), ["S1:2:1", "S2:2:0"])
        self.assertEqual(ids({"product_text": "TRANSPOSASE"}, ["S1", "S2"]), ["S1:3:0", "S2:3:0"])
        self.assertEqual(ids({"predictors": ["card_amr"]}, ["S2"]), ["S2:3:1"])
        # The marker search text is searched by the gene text too.
        self.assertEqual(ids({"gene_text": "aac(6')"}, ["S2"]), ["S2:3:1"])
        # No annotation filter: every annotation of the samples (D3: there is no
        # sample-query search of the annotations any more).
        self.assertEqual(ids({}, ["S1"]), ["S1:1:0", "S1:1:1", "S1:1:2", "S1:2:0", "S1:2:1", "S1:3:0", "S1:3:1"])
        # No allowed samples: an empty frame with the public columns.
        for samples in ([], ["", "  "]):
            empty = self.engine.filter_annotations({"gene_text": "blaOXA"}, samples)
            self.assertTrue(empty.is_empty())
            self.assertEqual(empty.columns, list(self.run.ann_public_columns))


def _wait_first_refresh(app: QtWidgets.QApplication, win: Any, timeout: float = 15.0) -> None:
    end = time.monotonic() + timeout
    while time.monotonic() < end:
        app.processEvents(QtCore.QEventLoop.ProcessEventsFlag.AllEvents, 20)
        if win.ctx.first_result_done and not win.ctx.scheduler.is_running() and win.ctx.scheduler.pending_filters() is None:
            return
        time.sleep(0.005)
    raise AssertionError("first refresh did not finish in time")


class WindowOverRunDataTests(_TempDirTestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        import annotation_explorer_polars.geo_experimental_views as gev

        cls._saved_webengine_cache = gev._WEBENGINE_RUNTIME_CACHE
        gev._WEBENGINE_RUNTIME_CACHE = (False, "disabled in tests")

    @classmethod
    def tearDownClass(cls) -> None:
        import annotation_explorer_polars.geo_experimental_views as gev

        gev._WEBENGINE_RUNTIME_CACHE = cls._saved_webengine_cache

    def _window(self, **kwargs: Any) -> Any:
        from annotation_explorer_polars.windows.explore_window import ExploreWindow

        win = ExploreWindow(**kwargs)
        self.addCleanup(self._close, win)
        _wait_first_refresh(self.app, win)
        return win

    def _close(self, win: Any) -> None:
        win.close()
        win.deleteLater()
        self.app.processEvents()

    def test_the_window_reads_run_data_through_run_only(self) -> None:
        # Phase 10: the window's delegating views (ann_df, _ann_pl, ...) are
        # gone; the mixins and every caller read win.run.
        run = build_parity_run(self.output_dir)
        win = self._window(run_data=run)
        self.assertIs(win.run, run)
        for attr in WINDOW_VIEWS:
            with self.subTest(attr=attr):
                self.assertFalse(hasattr(win, attr))
        # The tree's metadata columns are a window list (a tree option), seeded from the run.
        self.assertEqual(win.meta_columns_for_tree, list(run.meta_columns))

    def test_window_from_tables_builds_its_own_run_data(self) -> None:
        ann_df, meta_df = build_synthetic_run()
        win = self._window(ann_df=ann_df, meta_df=meta_df, ref_df=pd.DataFrame(), marker_hits_df=pd.DataFrame(), output_dir=self.output_dir)
        self.assertIsInstance(win.run, RunData)
        self.assertIsInstance(win.engine, QueryEngine)
        self.assertIs(win.engine.run, win.run)
        self.assertEqual(win.run.output_dir, self.output_dir)
        self.assertEqual(win.status.text(), "Samples: 6 | Annotations after filters: 31")

    def test_constructor_takes_tables_or_run_data(self) -> None:
        from annotation_explorer_polars.windows.explore_window import ExploreWindow

        run = build_parity_run(self.output_dir)
        ann_df, meta_df = build_synthetic_run()
        with self.assertRaises(ValueError):
            ExploreWindow(ann_df, meta_df, marker_hits_df=pd.DataFrame(), run_data=run)
        with self.assertRaises(ValueError):
            ExploreWindow()
        with self.assertRaises(ValueError):
            ExploreWindow(ann_df, meta_df)  # marker hits are required

    def test_windows_on_one_run_share_its_data_and_engine(self) -> None:
        run = build_parity_run(self.output_dir)
        engine = QueryEngine(run)
        computes: List[str] = []
        compute = engine._compute

        def recording(filters: FilterState) -> Any:
            computes.append(threading.current_thread().name)
            return compute(filters)

        engine._compute = recording  # type: ignore[method-assign]
        first = self._window(run_data=run, engine=engine)
        second = self._window(run_data=run, engine=engine)
        self.assertIs(first.run, second.run)
        self.assertIs(first.run.ann_df, second.run.ann_df)
        # The second window's first refresh has the same filters: served from the engine's cache.
        self.assertEqual(len(computes), 1)
        self.assertEqual(first.status.text(), second.status.text())
        self.assertEqual(sorted(first._sample_tree_items), sorted(second._sample_tree_items))
        # The compute ran on a pool thread, never on the GUI thread.
        self.assertNotEqual(computes[0], threading.main_thread().name)
        # Filters stay per window: a search in one leaves the other alone.
        first.gene_query.setText("blaOXA")
        end = time.monotonic() + 10
        while time.monotonic() < end and first.status.text() == second.status.text():
            self.app.processEvents(QtCore.QEventLoop.ProcessEventsFlag.AllEvents, 20)
            time.sleep(0.005)
        self.assertEqual(first.status.text(), "Samples: 3 | Annotations after filters: 4")
        self.assertEqual(second.status.text(), "Samples: 6 | Annotations after filters: 31")


    def test_engine_must_belong_to_the_run(self) -> None:
        from annotation_explorer_polars.windows.explore_window import ExploreWindow

        run = build_parity_run(self.output_dir)
        other = QueryEngine(build_parity_run(self.output_dir))
        with self.assertRaises(ValueError):
            ExploreWindow(run_data=run, engine=other)
        ann_df, meta_df = build_synthetic_run()
        with self.assertRaises(ValueError):
            ExploreWindow(ann_df, meta_df, marker_hits_df=pd.DataFrame(), engine=other)

    def test_window_on_ready_run_data_does_no_data_prep(self) -> None:
        run = build_parity_run(self.output_dir)
        with mock.patch.object(RunData, "build", side_effect=AssertionError("RunData.build on the GUI thread")):
            win = self._window(run_data=run)
        self.assertIs(win.run, run)

    def test_refresh_never_runs_window_code_on_the_pool_thread(self) -> None:
        from annotation_explorer_polars.windows.explore_window import ExploreWindow

        run = build_parity_run(self.output_dir)
        real_compute = QueryEngine.compute
        offenders: list = []

        def guarded_compute(engine_self, filters):
            # The pool thread runs the engine with nothing of the window on its stack.
            if threading.current_thread() is threading.main_thread():
                offenders.append("compute on the GUI thread")
            for info in inspect.stack():
                if isinstance(info.frame.f_locals.get("self"), ExploreWindow):
                    offenders.append(f"window code on the pool thread: {info.function}")
            return real_compute(engine_self, filters)

        with mock.patch.object(QueryEngine, "compute", guarded_compute):
            win = self._window(run_data=run)
            win.gene_query.setText("blaOXA")
            end = time.monotonic() + 10
            while time.monotonic() < end and "after filters: 4" not in win.status.text():
                self.app.processEvents(QtCore.QEventLoop.ProcessEventsFlag.AllEvents, 20)
                time.sleep(0.005)
        self.assertEqual(offenders, [])
        self.assertEqual(win.status.text(), "Samples: 3 | Annotations after filters: 4")

    def test_contig_rows_read_all_or_filtered_annotations(self) -> None:
        run = build_parity_run(self.output_dir)
        win = self._window(run_data=run)
        base = win.engine.result_annotations(win.ctx.result, hits=False, sample="S1", contig="3")
        self.assertEqual(base["AnnotationID"].tolist(), ["S1:3:0", "S1:3:1"])
        win.gene_query.setText("blaOXA")
        end = time.monotonic() + 10
        while time.monotonic() < end and "after filters: 4" not in win.status.text():
            self.app.processEvents(QtCore.QEventLoop.ProcessEventsFlag.AllEvents, 20)
            time.sleep(0.005)
        self.assertEqual(win.engine.result_annotations(win.ctx.result, hits=True, sample="S1", contig="3")["AnnotationID"].tolist(), ["S1:3:1"])
        self.assertEqual(win.engine.result_annotations(win.ctx.result, hits=False, sample="S1", contig="3")["AnnotationID"].tolist(), ["S1:3:0", "S1:3:1"])
        merged = win.engine.result_annotations(win.ctx.result, hits=False, sample="S1", contigs=["1", "3"])
        self.assertEqual(sorted(merged["Contig"].unique()), ["1", "3"])
        # The applied samples bound the base rows.
        win._set_sample_subset_filter(label="S2 only", samples=["S2"], kind="test")
        end = time.monotonic() + 10
        while time.monotonic() < end and list(win.ctx.result.allowed_samples) != ["S2"]:
            self.app.processEvents(QtCore.QEventLoop.ProcessEventsFlag.AllEvents, 20)
            time.sleep(0.005)
        self.assertEqual(list(win.ctx.result.allowed_samples), ["S2"])
        self.assertTrue(win.engine.result_annotations(win.ctx.result, hits=False, sample="S1", contig="3").empty)

    def test_gene_cards_open_with_the_run_data(self) -> None:
        from annotation_explorer_polars.annotation_detail_dialog import AnnotationDetailDialog

        (self.output_dir / "run_config.json").write_text(
            json.dumps({"inputs": {"snippy_reference": "/refs/S4.gbk"}}), encoding="utf-8"
        )
        run = build_parity_run(self.output_dir)
        self.assertEqual(run.snippy_reference_sample, "S4")
        win = self._window(run_data=run)
        windows = win.app_controller.windows
        self.addCleanup(lambda: [dlg.close() for dlg in windows.info_windows()])
        row = run.ann_df.loc[run.ann_df["AnnotationID"] == "S1:1:0"].iloc[0].to_dict()

        # From the track and from a table row: the same gene, so one card.
        win._on_igv_gene_double_clicked(dict(row))
        shown = windows.info_windows(win)
        self.assertEqual(len(shown), 1)
        win._set_table_model(
            [{"type": "data", "data": dict(row), "is_hit": False}], [(column, column) for column in ("Sample", "Gene", "Product")]
        )
        win._open_annotation_details_for_row(0)
        self.assertEqual(windows.info_windows(win), shown)
        dlg = shown[0]
        self.assertIsInstance(dlg, AnnotationDetailDialog)
        self.assertEqual(dlg._row_data["AnnotationID"], "S1:1:0")
        self.assertIs(dlg._ann_df, run.ann_df)
        self.assertIs(dlg._meta_df, run.meta_df)
        self.assertIs(dlg._marker_hits_df, run.marker_hits_df)
        # Its actions run on this window (the owner). Compare neighborhoods
        # opens a new window on the same run, whose neighborhood viewer gets
        # the run folder (for figure exports).
        self.assertIs(dlg.owner(), win)
        self.assertTrue(dlg.can_run_commands())
        self.assertTrue(dlg.run_command("window.open_neighbourhood", rows=[{"sample": "S1"}], title="t", source="gene_card"))
        new = [w for w in win.app_controller.windows.explore_windows() if w is not win]
        self.assertEqual(len(new), 1)
        self.addCleanup(new[0].close)
        self.assertIs(new[0].run, run)
        manager = new[0]._panel_manager
        viewer = next(manager.panel(f.uid) for f in manager.frames("neighbourhood"))
        options = viewer._view_options
        self.assertEqual(options["output_dir"], self.output_dir)
        # No tree copy in the viewer: the snippy reference is not passed.
        self.assertNotIn("pan_reference_sample", options)


class LoaderBuildsRunDataTests(_TempDirTestCase):
    """app._RunDataWorker builds the RunData off the GUI thread and hands it to the window."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def test_loader_emits_run_data(self) -> None:
        from annotation_explorer_polars import app as app_mod

        ann_df, meta_df = build_synthetic_run()
        marker_hits, ref_hits = build_marker_and_reference_hits()
        build_threads: List[str] = []
        real_build = RunData.build.__func__

        def recording_build(cls: type, *args: Any, **kwargs: Any) -> RunData:
            build_threads.append(threading.current_thread().name)
            return real_build(cls, *args, **kwargs)

        with mock.patch("annotation_explorer_polars.cache.cached_tables_exist", return_value=True), mock.patch(
            "annotation_explorer_polars.cache.read_cached_tables", return_value=(ann_df, ref_hits, marker_hits)
        ), mock.patch(
            "annotation_explorer_polars.data.prepare_metadata", return_value=(meta_df, ["a warning"], {})
        ), mock.patch(
            "annotation_explorer_polars.data.attach_assembly_quality_metadata", side_effect=lambda df, _dir: df
        ), mock.patch.object(RunData, "build", classmethod(recording_build)):
            worker = app_mod._RunDataWorker("run_x", self.output_dir, self.output_dir / "meta.csv")
            statuses: List[str] = []
            loaded: List[Any] = []
            failed: List[Any] = []
            worker.status.connect(statuses.append)
            worker.loaded.connect(loaded.append)
            worker.failed.connect(lambda *args: failed.append(args))
            worker.start()
            end = time.monotonic() + 30
            while time.monotonic() < end and not (loaded or failed):
                self.app.processEvents(QtCore.QEventLoop.ProcessEventsFlag.AllEvents, 20)
                time.sleep(0.005)
            worker.wait(10000)
        self.assertEqual(failed, [])
        self.assertEqual(len(loaded), 1)
        payload = loaded[0]
        self.assertEqual(sorted(payload), ["run_data", "warnings"])
        self.assertEqual(payload["warnings"], ["a warning"])
        run = payload["run_data"]
        self.assertIsInstance(run, RunData)
        self.assertEqual((run.run_id, run.output_dir), ("run_x", self.output_dir))
        self.assertEqual(run.annotation_db_columns, ("AMRFinderPlus_core", "ISFinder-Prot_MOBILE", "PAIDB", "VFDB-setB"))
        self.assertEqual(build_threads, [build_threads[0]])
        self.assertNotEqual(build_threads[0], threading.main_thread().name)
        self.assertEqual(
            statuses,
            ["Reading cached data…", "Preparing metadata…", "Preparing annotations…",
             "Indexing contigs and plasmids…", "Counting database hits…",
             # Phase 9: the pangenome folder lookups, also on the worker thread.
             "Finding pangenome outputs…"],
        )


if __name__ == "__main__":
    unittest.main()
