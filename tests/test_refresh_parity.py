"""Golden tests for the shared compute (window architecture plan, Phase 0 on).

They record ``QueryEngine.compute`` and the map payload on a small synthetic
run, over a matrix of view mode x sample query x annotation filters x subset x
hide-empty (plus a few metadata filter cases). The real ``ExploreWindow``
is built, so the phases that move the data preparation and the compute are
checked against the same recording.

Phase 3 (decision D3) took the annotation table's view mode out of the query:
the case ids keep their view mode prefix, and the Contig, Context and Global
rows of one case are now identical apart from the token (the prefix is only a
label: the recorder passes no view mode). The former Global rows that used the
sample query as an annotation query (and narrowed the contigs to the hits)
were re-recorded deliberately; their old records are kept in
golden/refresh_parity_d3_old_global.json and
``test_d3_old_global_cases_map_to_the_new_fields`` checks the old-to-new
mapping. In the same phase the sample tree's rows and the "Samples: N" count
became the window's one sample set (plan 4.4, D4), which re-recorded the tree
fields and the stats text of every case with an annotation filter that leaves
some samples without hits. Regenerate the recording with BACTOLOMEW_REGEN_GOLDEN=1 only when a
change is intended, and list the changed cases.
"""
from __future__ import annotations

import json
import math
import os
import shutil
import sys
import tempfile
import time
import unittest
from pathlib import Path
from typing import Any, Dict, List, Optional

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

import pandas as pd
from PyQt6 import QtCore, QtWidgets

from annotation_explorer_polars.constants import EXPECTED_ANNOTATION_COLUMNS, MISSING_FILTER_TOKEN
from annotation_explorer_polars.data import _finalize_annotations_df
from annotation_explorer_polars.model.filter_state import FilterState


GOLDEN_PATH = Path(__file__).resolve().parent / "golden" / "refresh_parity.json"
# The Global rows as recorded before D3 (Phase 3), for the old-to-new mapping.
D3_OLD_GLOBAL_PATH = Path(__file__).resolve().parent / "golden" / "refresh_parity_d3_old_global.json"
REGEN_ENV = "BACTOLOMEW_REGEN_GOLDEN"

# Synthetic run: six samples, three contigs each (S5 has contig "10" so the
# numeric contig sort is recorded too). Contig 3 of S1, S2 and S4 is a plasmid.
_GENES: Dict[str, Dict[str, Any]] = {
    "oxa": {
        "Gene": "blaOXA-23",
        "Product": "class D beta-lactamase OXA-23",
        "MarkerCategory": "AMR",
        "Predictors": "AMRF_core-AMR,CARD_AMR,COGorg24_COG",
        "Method": "AMRFinderPlus",
    },
    "lysr": {
        "Gene": "lysR",
        "Product": "LysR family transcriptional regulator",
        "Predictors": "COGorg24_COG",
        "Method": "COG",
    },
    "hypo": {
        "Gene": "",
        "Product": "hypothetical protein",
        "Predictors": "",
        "Method": "prodigal",
        "Hypothetical": True,
    },
    "ompa": {
        "Gene": "ompA",
        "Product": "outer membrane protein A",
        "MarkerCategory": "VIRULENCE",
        "Predictors": "VFDB-setB-Prot_VIRULENCE,COGorg24_COG",
        "Method": "VFDB",
    },
    "tnpa": {
        "Gene": "tnpA",
        "Product": "IS4 family transposase",
        "MarkerCategory": "MOBILE",
        "Predictors": "ISFinder-Prot_MOBILE",
        "Method": "ISFinder",
    },
    "aac": {
        "Gene": "",
        "Product": "aminoglycoside N-acetyltransferase",
        "MarkerCategory": "AMR",
        "Predictors": "CARD_AMR",
        "Method": "CARD",
        "MarkerSearchText": "aac(6')-ib card_amr",
    },
}

_LAYOUT: Dict[str, Dict[str, List[str]]] = {
    "S1": {"1": ["oxa", "lysr", "hypo"], "2": ["lysr", "ompa"], "3": ["tnpa", "oxa"]},
    "S2": {"1": ["lysr", "hypo"], "2": ["ompa", "hypo"], "3": ["tnpa", "aac"]},
    "S3": {"1": ["lysr"], "2": ["oxa", "hypo"], "3": ["lysr", "hypo"]},
    "S4": {"1": ["lysr", "ompa"], "2": ["hypo"], "3": ["oxa", "tnpa"]},
    "S5": {"1": ["lysr", "hypo"], "2": ["ompa"], "10": ["hypo"]},
    "S6": {"1": ["hypo"], "2": ["lysr"], "3": ["lysr", "hypo"]},
}

_PLASMIDS = {("S1", "3"): "AA001", ("S2", "3"): "AA002", ("S4", "3"): "AA001"}

_META = [
    # Sample, Country, Year, Source, Latitude, Longitude
    ("S1", "Romania", "2019", "blood", 44.43, 26.10),
    ("S2", "Italy", "2020", "urine", None, None),
    ("S3", "Spain", "2020", "blood", None, None),
    ("S4", "Romania", "2021", "wound", 45.75, 21.23),
    ("S5", "", "2021", "urine", None, None),
    ("S6", "Italy", "2019", "sputum", 41.90, 12.50),
]


def build_synthetic_run() -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: List[Dict[str, Any]] = []
    for sample, contigs in _LAYOUT.items():
        for contig, genes in contigs.items():
            for idx, key in enumerate(genes):
                row: Dict[str, Any] = {col: "" for col in EXPECTED_ANNOTATION_COLUMNS}
                row.update(
                    {
                        "Sample": sample,
                        "Contig": contig,
                        "Start": 1 + idx * 1000,
                        "Stop": 900 + idx * 1000,
                        "Strand": "+" if idx % 2 == 0 else "-",
                        "Hypothetical": False,
                        "AnnotationID": f"{sample}:{contig}:{idx}",
                        "Locus": f"{sample}_{contig}_{idx:03d}",
                        "MobSuiteIsPlasmid": "True" if (sample, contig) in _PLASMIDS else "False",
                        "MobSuitePlasmidID": _PLASMIDS.get((sample, contig), ""),
                        "MobSuiteMoleculeType": "plasmid" if (sample, contig) in _PLASMIDS else "chromosome",
                    }
                )
                row.update(_GENES[key])
                rows.append(row)
    # The loader's own finalisation adds AnnUID, StartInt, StopInt and
    # MarkerDisplay and derives Tags from the categories, as for a real run.
    ann_df = _finalize_annotations_df(pd.DataFrame(rows))
    ann_df["Hypothetical"] = ann_df["Hypothetical"].astype(bool)
    meta_df = pd.DataFrame(
        _META,
        columns=["Sample", "Country", "Year", "Source", "Latitude", "Longitude"],
    )
    return ann_df, meta_df


def _ann_filters(**changes: Any) -> Dict[str, Any]:
    base: Dict[str, Any] = {
        "gene_text": None,
        "product_text": None,
        "molecule": "",
        "predictors": [],
        "tags": [],
        "hide_hypothetical": False,
    }
    base.update(changes)
    return base


_VIEW_MODES = ("Contig", "Context", "Global")
# "romania" only matches metadata, "blaoxa" only matches annotations, "s2"
# matches the sample name in both.
_SAMPLE_QUERIES: tuple[Optional[str], ...] = (None, "romania", "blaoxa", "s2")
_ANN_FILTERS: Dict[str, Dict[str, Any]] = {
    "none": _ann_filters(),
    "gene": _ann_filters(gene_text="blaOXA"),
    "gene_marker_text": _ann_filters(gene_text="aac(6')"),
    "product": _ann_filters(product_text="transposase"),
    "plasmid": _ann_filters(molecule="plasmid"),
    "chromosome_tag_v": _ann_filters(molecule="chromosome", tags=["V"]),
    "predictor_hide_hypo": _ann_filters(predictors=["CARD_AMR"], hide_hypothetical=True),
}
_SUBSETS: Dict[str, List[str]] = {"all": [], "s1_s3": ["S1", "S2", "S3"]}
_META_FILTERS: Dict[str, Dict[str, List[str]]] = {
    "country_romania": {"Country": ["Romania"]},
    "country_missing": {"Country": [MISSING_FILTER_TOKEN]},
    "country_italy_or_missing": {"Country": ["Italy", MISSING_FILTER_TOKEN]},
}
# Combinations the main matrix does not cross. The sample query together with a
# metadata filter guards D3 (Phase 3 moves the sample query into the metadata
# filter path).
_EXTRA_CASES: Dict[str, Dict[str, Any]] = {
    "q=romania|meta=country_romania": {"sample_query": "romania", "meta_filters": {"Country": ["Romania"]}},
    "q=s2|meta=country_italy_or_missing": {
        "sample_query": "s2",
        "meta_filters": {"Country": ["Italy", MISSING_FILTER_TOKEN]},
    },
    "q=s2|meta=country_romania": {"sample_query": "s2", "meta_filters": {"Country": ["Romania"]}},
    "subset=s1_s3|meta=country_italy_or_missing": {
        "subset": ["S1", "S2", "S3"],
        "meta_filters": {"Country": ["Italy", MISSING_FILTER_TOKEN]},
    },
    "ann=gene_and_product": {"ann_filters": _ann_filters(gene_text="blaOXA", product_text="beta-lactamase")},
    "ann=hide_hypo_only": {"ann_filters": _ann_filters(hide_hypothetical=True)},
    "meta=year_2020": {"meta_filters": {"Year": ["2020"]}},
    "meta=source_blood|ann=gene": {
        "meta_filters": {"Source": ["blood"]},
        "ann_filters": _ann_filters(gene_text="blaOXA"),
    },
    "subset=with_unknown_sample": {"subset": ["S1", "S9"]},
    "subset=with_unknown_sample|ann=gene": {"subset": ["S1", "S9"], "ann_filters": _ann_filters(gene_text="blaOXA")},
}


def build_case_matrix() -> List[Dict[str, Any]]:
    cases: List[Dict[str, Any]] = []
    for view_mode in _VIEW_MODES:
        for sample_query in _SAMPLE_QUERIES:
            for ann_name, ann in _ANN_FILTERS.items():
                for subset_name, subset in _SUBSETS.items():
                    for hide_empty in (False, True):
                        cases.append(
                            {
                                "id": f"{view_mode}|q={sample_query or ''}|ann={ann_name}|subset={subset_name}|hide_empty={int(hide_empty)}",
                                "view_mode": view_mode,
                                "sample_query": sample_query,
                                "meta_filters": {},
                                "subset": list(subset),
                                "ann_filters": dict(ann),
                                "hide_empty": hide_empty,
                            }
                        )
    for view_mode in _VIEW_MODES:
        for meta_name, meta_filters in _META_FILTERS.items():
            for ann_name in ("none", "gene"):
                cases.append(
                    {
                        "id": f"{view_mode}|meta={meta_name}|ann={ann_name}",
                        "view_mode": view_mode,
                        "sample_query": None,
                        "meta_filters": {k: list(v) for k, v in meta_filters.items()},
                        "subset": [],
                        "ann_filters": dict(_ANN_FILTERS[ann_name]),
                        "hide_empty": False,
                    }
                )
    for view_mode in _VIEW_MODES:
        for extra_name, extra in _EXTRA_CASES.items():
            case = {
                "id": f"{view_mode}|extra|{extra_name}",
                "view_mode": view_mode,
                "sample_query": None,
                "meta_filters": {},
                "subset": [],
                "ann_filters": _ann_filters(),
                "hide_empty": False,
            }
            case.update(
                {
                    key: (dict(value) if isinstance(value, dict) else list(value) if isinstance(value, list) else value)
                    for key, value in extra.items()
                }
            )
            cases.append(case)
    return cases


def _jsonable(value: Any) -> Any:
    if isinstance(value, float) and math.isnan(value):
        return None
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        return value.item()
    return value


def case_filters(case: Dict[str, Any]) -> FilterState:
    """The filters of one case (its view mode is a table display option, D3)."""
    return FilterState.from_query(
        sample_query=case["sample_query"],
        meta_filters={k: list(v) for k, v in case["meta_filters"].items()},
        subset_samples=list(case["subset"]),
        ann_filters=dict(case["ann_filters"]),
        hide_empty_contigs=bool(case["hide_empty"]),
    )


def record_result(engine: Any, filters: FilterState, token: int) -> Dict[str, Any]:
    result = engine.compute(filters)
    # D4: the map shows the window's one sample set.
    map_payload = engine.map_payload(result.samples)
    contig_records = {
        sample: [
            [
                str(row.get("Contig")),
                int(row.get("n_hits") or 0),
                int(row.get("n_total") or 0),
                bool(row.get("MobSuiteIsPlasmid")),
                str(row.get("MobSuitePlasmidID") or ""),
            ]
            for row in rows
        ]
        for sample, rows in result.contig_records.items()
    }
    # Sample lists and map points are in metadata order since Phase 3 (they
    # were recorded sorted before, when their order was not stable; the
    # synthetic metadata is in sample order, so the records did not change).
    return _jsonable(
        {
            "token": token,
            "allowed_samples": list(result.allowed_samples),
            "panel_allowed_samples": list(result.samples),
            # The sample tree lists the window's one sample set (plan 4.4, D4),
            # in every view mode (D3); the window test
            # test_widget_snapshot_matches_golden_case checks the real tree.
            "tree_allowed_samples": list(result.samples),
            "meta_for_tree_samples": result.metadata.get_column("Sample").to_list()
            if "Sample" in result.metadata.columns
            else [],
            "filtered_ann_ids": result.filtered_ann.get_column("AnnotationID").to_list()
            if "AnnotationID" in result.filtered_ann.columns
            else [],
            "filtered_counts": dict(sorted(result.filtered_counts.items())),
            "contig_records": contig_records,
            "stats_text": result.stats_text,
            "ann_filters": filters.ann_filters(),
            "map": map_payload,
        }
    )


def record_payload(engine: Any, case: Dict[str, Any], token: int) -> Dict[str, Any]:
    return record_result(engine, case_filters(case), token)


def record_all(engine: Any) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    for token, case in enumerate(build_case_matrix(), start=1):
        out[case["id"]] = record_payload(engine, case, token)
    return out


def write_golden(recorded: Dict[str, Dict[str, Any]]) -> None:
    GOLDEN_PATH.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"  {json.dumps(case_id)}: {json.dumps(value, sort_keys=True, ensure_ascii=False)}"
        for case_id, value in recorded.items()
    ]
    GOLDEN_PATH.write_text("{\n" + ",\n".join(lines) + "\n}\n", encoding="utf-8")


class _ExplorerFixture:
    """Builds the real ExploreWindow on the synthetic run (QSettings are
    redirected for the whole process by isolate_qt_settings at import)."""

    app: QtWidgets.QApplication
    output_dir: str

    @classmethod
    def build(cls, testcase: type) -> Any:
        testcase.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        import annotation_explorer_polars.geo_experimental_views as gev

        testcase._saved_webengine_cache = gev._WEBENGINE_RUNTIME_CACHE
        gev._WEBENGINE_RUNTIME_CACHE = (False, "disabled in tests")
        from annotation_explorer_polars.windows.explore_window import ExploreWindow

        testcase.output_dir = tempfile.mkdtemp(prefix="bactolomew_parity_run_")
        ann_df, meta_df = build_synthetic_run()
        win = ExploreWindow(
            ann_df,
            meta_df,
            ref_df=pd.DataFrame(),
            marker_hits_df=pd.DataFrame(),
            output_dir=Path(testcase.output_dir),
        )
        return win

    @classmethod
    def teardown(cls, testcase: type, win: Any) -> None:
        import annotation_explorer_polars.geo_experimental_views as gev

        if win is not None:
            win.close()
            win.deleteLater()
            testcase.app.processEvents()
        gev._WEBENGINE_RUNTIME_CACHE = testcase._saved_webengine_cache
        shutil.rmtree(testcase.output_dir, ignore_errors=True)


def _wait_idle(app: QtWidgets.QApplication, win: Any, timeout: float = 15.0) -> None:
    """Process events until no refresh is debounced, queued or running."""
    end = time.monotonic() + timeout
    quiet_since: Optional[float] = None
    while time.monotonic() < end:
        app.processEvents(QtCore.QEventLoop.ProcessEventsFlag.AllEvents, 20)
        busy = (
            win.ctx.scheduler.debounce_timer.isActive()
            or win.ctx.scheduler.is_running()
            or win.ctx.scheduler.pending_filters() is not None
            or not win.ctx.first_result_done
        )
        if busy:
            quiet_since = None
        elif quiet_since is None:
            quiet_since = time.monotonic()
        elif time.monotonic() - quiet_since > 0.05:
            return
        time.sleep(0.005)
    raise AssertionError("refresh did not finish in time")


def _ann_active(case: Dict[str, Any]) -> bool:
    ann = case["ann_filters"]
    return bool(
        ann.get("gene_text")
        or ann.get("product_text")
        or ann.get("molecule")
        or ann.get("predictors")
        or ann.get("tags")
        or ann.get("hide_hypothetical")
    )


def _samples_of(ann_ids: List[str]) -> List[str]:
    return sorted({str(ann_id).split(":", 1)[0] for ann_id in ann_ids})


class RefreshParityGoldenTests(unittest.TestCase):
    win: Any = None
    regenerated = False

    @classmethod
    def setUpClass(cls) -> None:
        cls.win = _ExplorerFixture.build(cls)
        _wait_idle(cls.app, cls.win)
        cls.recorded = record_all(cls.win.engine)
        cls.cases = {case["id"]: case for case in build_case_matrix()}
        if os.environ.get(REGEN_ENV) == "1":
            write_golden(cls.recorded)
            cls.regenerated = True
            print(
                f"\nWARNING: {REGEN_ENV}=1: rewrote {GOLDEN_PATH} ({len(cls.recorded)} cases); "
                "the golden comparison is skipped in this run.",
                file=sys.stderr,
            )

    @classmethod
    def tearDownClass(cls) -> None:
        _ExplorerFixture.teardown(cls, cls.win)

    def test_golden_file_exists(self) -> None:
        self.assertTrue(
            GOLDEN_PATH.is_file(),
            f"{GOLDEN_PATH} is missing; generate it with {REGEN_ENV}=1",
        )

    def test_matrix_matches_golden(self) -> None:
        if self.regenerated:
            self.skipTest(f"golden file regenerated ({REGEN_ENV}=1); nothing to compare")
        golden = json.loads(GOLDEN_PATH.read_text(encoding="utf-8"))
        self.assertEqual(sorted(golden), sorted(self.recorded), "case ids differ from the golden file")
        for case_id, expected in golden.items():
            with self.subTest(case=case_id):
                self.assertEqual(self.recorded[case_id], expected)

    def test_compute_is_deterministic_with_warm_caches(self) -> None:
        # A second pass hits the metadata and annotation LRU caches.
        self.assertEqual(record_all(self.win.engine), self.recorded)

    def test_d3_view_mode_does_not_change_the_result(self) -> None:
        # D3: Contig, Around matches and All matches only change how the
        # annotation table displays the result, never the result. The recorder
        # cannot pass a view mode at all (the case id prefix is only a label),
        # so the guard is the shape of the API: FilterState has no view mode
        # and compute takes nothing else. The window tests in
        # test_explore_context.py check that a mode switch starts no compute.
        import dataclasses
        import inspect

        self.assertNotIn("view_mode", {f.name for f in dataclasses.fields(FilterState)})
        self.assertEqual(list(inspect.signature(self.win.engine.compute).parameters), ["filters"])
        by_case: Dict[str, Dict[str, Dict[str, Any]]] = {}
        for case_id, rec in self.recorded.items():
            mode, rest = case_id.split("|", 1)
            by_case.setdefault(rest, {})[mode] = {k: v for k, v in rec.items() if k != "token"}
        self.assertEqual(len(by_case), len(self.recorded) // 3)
        for rest, modes in by_case.items():
            with self.subTest(case=rest):
                self.assertEqual(sorted(modes), sorted(_VIEW_MODES))
                self.assertEqual(modes["Global"], modes["Contig"])
                self.assertEqual(modes["Context"], modes["Contig"])

    def test_d3_sample_query_filters_samples_only(self) -> None:
        # The sample box is a sample and metadata search in every mode.
        romania = self.recorded["Global|q=romania|ann=none|subset=all|hide_empty=0"]
        self.assertEqual(romania["allowed_samples"], ["S1", "S4"])
        self.assertEqual(romania["panel_allowed_samples"], ["S1", "S4"])
        self.assertEqual(_samples_of(romania["filtered_ann_ids"]), ["S1", "S4"])
        # A gene name in the sample box matches no sample (the hint case).
        oxa = self.recorded["Global|q=blaoxa|ann=none|subset=all|hide_empty=0"]
        self.assertEqual((oxa["allowed_samples"], oxa["panel_allowed_samples"], oxa["filtered_ann_ids"]), ([], [], []))
        self.assertEqual(oxa["stats_text"], "Samples: 0 | Annotations after filters: 0")
        self.assertTrue(self.win.engine.compute(FilterState(sample_query="blaoxa")).sample_query_matches_annotations)
        self.assertFalse(self.win.engine.compute(FilterState(sample_query="romania")).sample_query_matches_annotations)
        self.assertFalse(self.win.engine.compute(FilterState(sample_query="no-such-text")).sample_query_matches_annotations)
        # The gene search is the way to find genes.
        gene = self.recorded["Global|q=|ann=gene|subset=all|hide_empty=0"]
        self.assertEqual(gene["panel_allowed_samples"], ["S1", "S3", "S4"])

    def test_d3_old_global_cases_map_to_the_new_fields(self) -> None:
        # Old-to-new mapping (plan 4.3 parity note): what the former Global
        # sample query found is found by the same text in "Gene/annotation
        # contains" when it named a gene (blaoxa; romania found nothing), and
        # by the sample search itself when it named a sample (s2).
        old = json.loads(D3_OLD_GLOBAL_PATH.read_text(encoding="utf-8"))
        checked = 0
        for case_id, old_rec in old.items():
            case = self.cases[case_id]
            query = case["sample_query"]
            if not query:
                continue
            ann = dict(case["ann_filters"])
            if query == "s2":
                filters = case_filters(case)
            elif ann.get("gene_text"):
                continue  # the case already uses the gene field
            else:
                ann["gene_text"] = query
                filters = case_filters({**case, "sample_query": None, "ann_filters": ann})
            with self.subTest(case=case_id):
                new = self.win.engine.compute(filters)
                self.assertEqual(new.filtered_ann.get_column("AnnotationID").to_list(), old_rec["filtered_ann_ids"])
                checked += 1
        # 87 former Global cases with a sample query, less the 16 that already
        # had a gene text.
        self.assertEqual(checked, 71)


    # Named rules from plan 4.3. The golden equality already covers them; these
    # say which rule moved when a later phase changes a recording on purpose.

    def test_rule_panel_samples_shrink_only_with_annotation_filter(self) -> None:
        for case_id, rec in self.recorded.items():
            case = self.cases[case_id]
            with self.subTest(case=case_id):
                if _ann_active(case):
                    self.assertEqual(rec["panel_allowed_samples"], _samples_of(rec["filtered_ann_ids"]))
                    self.assertLessEqual(set(rec["panel_allowed_samples"]), set(rec["allowed_samples"]))
                else:
                    self.assertEqual(rec["panel_allowed_samples"], rec["allowed_samples"])

    def test_rule_tree_and_status_show_the_sample_set(self) -> None:
        # Plan 4.4 and D4: the tree's rows and the status bar count are the
        # window's one sample set, in every view mode (D3).
        for case_id, rec in self.recorded.items():
            with self.subTest(case=case_id):
                self.assertEqual(rec["tree_allowed_samples"], rec["panel_allowed_samples"])
                self.assertEqual(rec["meta_for_tree_samples"], rec["panel_allowed_samples"])
                self.assertTrue(rec["stats_text"].startswith(f"Samples: {len(rec['panel_allowed_samples'])} |"))

    def test_rule_map_shows_the_sample_set(self) -> None:
        # D4: the map counts exactly the window's sample set.
        for case_id, rec in self.recorded.items():
            with self.subTest(case=case_id):
                if rec["panel_allowed_samples"]:
                    self.assertEqual(rec["map"]["total"], len(rec["panel_allowed_samples"]))
                points = [p["sample"] for p in rec["map"]["points"]]
                self.assertLessEqual(set(points), set(rec["panel_allowed_samples"]))

    def test_rule_hide_empty_drops_contigs_without_hits(self) -> None:
        for case_id, rec in self.recorded.items():
            case = self.cases[case_id]
            if not _ann_active(case):
                continue
            with self.subTest(case=case_id):
                counts = [hits for contigs in rec["contig_records"].values() for _c, hits, *_r in contigs]
                if case["hide_empty"]:
                    self.assertTrue(all(hits > 0 for hits in counts))
        shown = self.recorded["Contig|q=|ann=gene|subset=all|hide_empty=0"]["contig_records"]
        hidden = self.recorded["Contig|q=|ann=gene|subset=all|hide_empty=1"]["contig_records"]
        self.assertIn(["2", 0, 2, False, ""], shown["S1"])
        self.assertEqual(hidden["S1"], [["1", 1, 3, False, ""], ["3", 1, 2, True, "AA001"]])
        self.assertNotIn("S2", hidden)
        # Only "Hide contigs without hits" drops contigs (D3: All matches keeps them).
        self.assertEqual(self.recorded["Global|q=|ann=gene|subset=all|hide_empty=0"]["contig_records"], shown)

    def test_global_status_text_with_search(self) -> None:
        self.assertEqual(
            self.recorded["Global|q=|ann=gene|subset=all|hide_empty=0"]["stats_text"],
            "Samples: 3 | Annotations after filters: 4",
        )


class RefreshParityWindowTests(unittest.TestCase):
    """Window behaviours that 4.3 of the plan lists next to the compute."""

    win: Any = None

    @classmethod
    def setUpClass(cls) -> None:
        cls.win = _ExplorerFixture.build(cls)
        cls.win.resize(1400, 900)
        _wait_idle(cls.app, cls.win)

    @classmethod
    def tearDownClass(cls) -> None:
        _ExplorerFixture.teardown(cls, cls.win)

    def setUp(self) -> None:
        win = self.win
        win.sample_query.clear()
        win.gene_query.clear()
        win.product_query.clear()
        win._panel_manager.panel("sample_tree").hide_empty_contigs_chk.setChecked(False)
        win._clear_sample_subset_filter()
        if win.view_mode != "Contig":
            win.annotations_panel.set_view_mode("Contig")
        win.refresh()
        _wait_idle(self.app, win)

    def _tree_samples(self) -> List[str]:
        return sorted(self.win._sample_tree_items)

    def _table_samples_and_contigs(self) -> set:
        model = self.win.table.model()
        out = set()
        for row in range(model.rowCount()):
            data = model.row_data(row) or {}
            if data.get("Sample"):
                out.add((str(data.get("Sample")), str(data.get("Contig"))))
        return out

    def test_global_without_search_clears_table(self) -> None:
        from annotation_explorer_polars.table_controller import ALL_MATCHES_EMPTY_TEXT

        win = self.win
        stats = win.status.text()
        win.annotations_panel.set_view_mode("Global")
        _wait_idle(self.app, win)
        self.assertEqual(win.view_mode, "Global")
        self.assertEqual(win.table.model().rowCount(), 0)
        # The hint is the Annotations panel's own line (D3), not the shared status bar.
        self.assertEqual(win.annotations_panel.status_text(), ALL_MATCHES_EMPTY_TEXT)
        self.assertEqual(win.status.text(), stats)
        # The tree still lists every sample: nothing narrows it without a search.
        self.assertEqual(self._tree_samples(), ["S1", "S2", "S3", "S4", "S5", "S6"])

    def test_d3_gene_name_in_the_sample_search_filters_samples(self) -> None:
        # Before D3 this searched the annotations in All matches; now the
        # sample box filters samples only, in every mode.
        win = self.win
        win.annotations_panel.set_view_mode("Global")
        win.sample_query.setText("blaoxa")
        _wait_idle(self.app, win)
        self.assertEqual(win.table.model().rowCount(), 0)
        self.assertEqual(self._tree_samples(), [])
        self.assertEqual(win.status.text(), "Samples: 0 | Annotations after filters: 0")
        # The gene field is where genes are searched.
        win.sample_query.clear()
        win.gene_query.setText("blaoxa")
        _wait_idle(self.app, win)
        self.assertEqual(win.table.model().rowCount(), 4)
        self.assertEqual(win.annotations_panel.status_text(), "4 annotations across 3 samples")
        self.assertEqual(win.status.text(), "Samples: 3 | Annotations after filters: 4")
        # The tree shows the window's one sample set in every mode (plan 4.4).
        self.assertEqual(self._tree_samples(), ["S1", "S3", "S4"])
        win.annotations_panel.set_view_mode("Contig")
        self.assertEqual(self._tree_samples(), ["S1", "S3", "S4"])

    def test_widget_snapshot_matches_golden_case(self) -> None:
        win = self.win
        win.annotations_panel.set_view_mode("Context")
        win.gene_query.setText("blaOXA")
        win._panel_manager.panel("sample_tree").hide_empty_contigs_chk.setChecked(True)
        win._set_sample_subset_filter(label="Three", samples=["S1", "S2", "S3"], kind="test")
        _wait_idle(self.app, win)
        golden = json.loads(GOLDEN_PATH.read_text(encoding="utf-8"))
        expected = golden["Context|q=|ann=gene|subset=s1_s3|hide_empty=1"]
        self.assertEqual(sorted(list(win.ctx.result.allowed_samples)), expected["allowed_samples"])
        self.assertEqual(win.ctx.result.filtered_ann.get_column("AnnotationID").to_list(), expected["filtered_ann_ids"])
        self.assertEqual(self._tree_samples(), sorted(expected["tree_allowed_samples"]))
        map_payload = win.engine.map_payload(win.ctx.sample_set())
        map_payload["points"] = sorted(map_payload["points"], key=lambda p: p["sample"])
        self.assertEqual(map_payload, expected["map"])

    def test_open_in_context_view_switches_mode_and_loads_contig(self) -> None:
        win = self.win
        self.assertEqual(win.view_mode, "Contig")
        win._open_tree_contig_payload({"sample": "S3", "contig": "2"}, mode="Context")
        _wait_idle(self.app, win)
        self.assertEqual(win.view_mode, "Context")
        self.assertEqual(win.view_mode, "Context")
        self.assertEqual((win.ctx.selection.sample, win.ctx.selection.contig), ("S3", "2"))
        self.assertEqual(self._table_samples_and_contigs(), {("S3", "2")})

    def test_sample_card_shows_summary_and_raw_only(self) -> None:
        # D10: the card is Summary and Raw; the deleted resistance and markers
        # tab builders must not come back, and the caller's signature must hold.
        win = self.win
        dlg = win._open_sample_card_for_sample("S1")
        self.app.processEvents()
        self.assertEqual(win.app_controller.windows.info_windows(win), [dlg])
        try:
            tabs = dlg.findChild(QtWidgets.QTabWidget)
            self.assertEqual([tabs.tabText(i) for i in range(tabs.count())], ["Summary", "Raw"])
            raw = dlg.findChild(QtWidgets.QPlainTextEdit)
            self.assertIn('"S1"', raw.toPlainText())
        finally:
            dlg.close()
            self.app.processEvents()

    def test_sample_card_database_summary_table(self) -> None:
        from annotation_explorer_polars.sample_card_dialog import SampleCardDialog

        ann_df, meta_df = build_synthetic_run()
        hits = pd.DataFrame(
            {
                "Sample": ["S1", "S1", "S1", "S2"],
                "Contig": ["1", "3", "3", "1"],
                "DB": ["CARD", "CARD", "VFDB", "CARD"],
            }
        )
        dlg = SampleCardDialog(None, sample="S1", ann_df=ann_df, meta_df=meta_df, marker_hits_df=hits)
        try:
            tabs = dlg.findChild(QtWidgets.QTabWidget)
            self.assertEqual([tabs.tabText(i) for i in range(tabs.count())], ["Summary", "Raw"])
            tables = [
                t for t in dlg.findChildren(QtWidgets.QTableWidget)
                if [t.horizontalHeaderItem(i).text() for i in range(t.columnCount())] == ["Type", "DB", "Hits", "Contigs"]
            ]
            self.assertEqual(len(tables), 1)
            table = tables[0]
            rows = [
                [table.item(r, c).text() for c in range(table.columnCount())]
                for r in range(table.rowCount())
            ]
            self.assertEqual(sorted(rows), [["Marker", "CARD", "2", "2"], ["Marker", "VFDB", "1", "1"]])
        finally:
            dlg.deleteLater()
            self.app.processEvents()

    def test_only_one_sample_subset_at_a_time(self) -> None:
        win = self.win
        win._set_sample_subset_filter(label="First", samples=["S1", "S2"], kind="test")
        _wait_idle(self.app, win)
        self.assertEqual(self._tree_samples(), ["S1", "S2"])
        win._set_sample_subset_filter(label="Second", samples=["S4", "S5"], kind="test")
        _wait_idle(self.app, win)
        self.assertEqual(list(win.ctx.draft.subset.samples), ["S4", "S5"])
        self.assertEqual(win.ctx.draft.subset.label, "Second")
        self.assertEqual(self._tree_samples(), ["S4", "S5"])
        chips = []
        for idx in range(win.filters_layout.count()):
            widget = win.filters_layout.itemAt(idx).widget()
            if widget is None:
                continue
            labels = [label.text() for label in widget.findChildren(QtWidgets.QLabel)]
            if "Sample set" in labels:
                chips.append(labels)
        self.assertEqual(len(chips), 1)
        self.assertIn("Second | n=2", chips[0])
        win._clear_sample_subset_filter()
        _wait_idle(self.app, win)
        self.assertEqual(self._tree_samples(), ["S1", "S2", "S3", "S4", "S5", "S6"])


if __name__ == "__main__":
    unittest.main()
