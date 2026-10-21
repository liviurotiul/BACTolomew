"""Gene evidence from Bakta and AMRFinderPlus (step bakta_evidence) gives the annotation
tags without any marker database; a run made before the step loads as before, and a hit
table left by a database removed from the run is not read."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from pipeline.BaktaEvidence import build_evidence, evidence_record


def _cds(locus: str, start: int, product: str, **extra) -> dict:
    return {"type": "cds", "sequence": "contig_1", "start": start, "stop": start + 299, "strand": "+",
            "locus": locus, "gene": extra.pop("gene", None), "product": product, **extra}


BAKTA = {
    "version": {"bakta": "1.10.4", "db": {"version": "6.0", "type": "full"}},
    "features": [
        _cds("S1_0001", 1, "general secretion pathway protein E", gene="gspE",
             expert=[{"type": "expert_proteins", "source": "VFDB", "id": "VFG050757", "gene": "gspE1",
                      "product": "general secretion pathway protein E", "identity": 1.0, "subject_cov": 1.0,
                      "db_xrefs": ["VFDB:VFG050757", "VFDB:VFC0086"]}]),
        _cds("S1_0002", 401, "IS66 family ISAba17 transposase", db_xrefs=["IS:ISAba17"]),
        _cds("S1_0003", 801, "phage tail tape measure protein"),
        _cds("S1_0004", 1201, "phage integrase family protein"),  # a generic integrase: no tag
        _cds("S1_0005", 1601, "class D beta-lactamase OXA-23", gene="blaOXA-23"),
        _cds("S1_0006", 2001, "type VI secretion system baseplate subunit TssK"),  # not phage
    ],
}
AMRFINDER_HEADER = ("Name\tProtein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\t"
                    "Subtype\tClass\tSubclass\tMethod\t% Coverage of reference\t% Identity to reference\t"
                    "Closest reference accession\n")
AMRFINDER_ROWS = (
    "S1\tS1_0005\tcontig_1\t1601\t1900\t+\tblaOXA-23\tOXA-23\tcore\tAMR\tAMR\tBETA-LACTAM\tCARBAPENEM\tEXACTP\t100.00\t100.00\tWP_1\n"
)


class BaktaEvidenceTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.run = Path(self._tmp.name) / "run"
        (self.run / "bakta" / "S1").mkdir(parents=True)
        (self.run / "bakta" / "S1" / "S1.json").write_text(json.dumps(BAKTA), encoding="utf-8")
        (self.run / "amrfinder").mkdir()
        (self.run / "amrfinder" / "S1.tsv").write_text(AMRFINDER_HEADER + AMRFINDER_ROWS, encoding="utf-8")
        (self.run / "amrfinder" / "summary.json").write_text(
            json.dumps({"organism": "Acinetobacter_baumannii", "samples": ["S1"]}), encoding="utf-8")
        # MOB-suite found a relaxase on contig_1, which it calls chromosome: no PLASMID row.
        mob = self.run / "mobsuite" / "per_sample" / "S1"
        mob.mkdir(parents=True)
        (mob / "biomarkers.blast.txt").write_text(
            "qseqid\tsseqid\tsstart\tsend\tpident\tqcovhsp\tbiomarker\n"
            "r1|MOBP\tcontig_1\t1601\t1900\t99.0\t100.0\trelaxase\n", encoding="utf-8")
        (mob / "contig_report.txt").write_text(
            "sample_id\tmolecule_type\tcontig_id\nS1\tchromosome\tcontig_1\n", encoding="utf-8")
        # The old marker search: a CARD hit on the OXA gene, and the table of a database
        # removed from the run since.
        hits = self.run / "annotations" / "bakta_vs_queries"
        hits.mkdir(parents=True)
        (hits / "bakta_annotations_nt-vs-CARD_AMR.tsv").write_text(
            "BKT0000005\tHM|ARO:3001396|OXA-23\t99.0\t273\t0\t0\t1\t273\t1\t273\t1e-100\t500\t273\t273\n", encoding="utf-8")
        (hits / "bakta_annotations_nt-vs-VFDB-setB-Prot_VIRULENCE.tsv").write_text(
            "BKT0000003\tVFG1\t99.0\t99\t0\t0\t1\t99\t1\t99\t1e-50\t200\t99\t99\n", encoding="utf-8")
        self.run_cfg = {"inputs": {
            "genomes": [str(self.run / "S1.fna")],
            "marker_sequences": ["/nowhere/CARD_AMR.fasta"],
            "reference_locus_sequences": [],
            # No catalog lists CARD any more: the run's recorded category is used.
            "marker_categories": {"CARD_AMR": "AMR"},
        }}

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def _load(self):
        from annotation_explorer_polars import data

        original = data.MARKER_DB_CATALOG_PATH
        data.MARKER_DB_CATALOG_PATH = Path(self._tmp.name) / "no_catalog.json"
        try:
            ann, hits = data.build_annotations_from_output(self.run, self.run_cfg)
        finally:
            data.MARKER_DB_CATALOG_PATH = original
        return dict(zip(ann["Locus"], ann["Tags"].fillna(""))), hits

    def test_evidence_rows_and_tags(self) -> None:
        # A run made before the step: the CARD hit tags the OXA gene, the removed VFDB table is not read.
        tags, hits = self._load()
        self.assertEqual(tags, {"S1_0001": "", "S1_0002": "", "S1_0003": "", "S1_0004": "", "S1_0005": "A", "S1_0006": ""})
        self.assertEqual(set(hits["DB"]), {"CARD_AMR"})

        summary = build_evidence(self.run, ["S1"])
        rows = (self.run / "evidence" / "S1.tsv").read_text(encoding="utf-8").splitlines()
        found = {(r.split("\t")[5], r.split("\t")[8]) for r in rows[1:]}
        self.assertEqual(found, {
            ("S1_0001", "Bakta VFDB"), ("S1_0002", "Bakta IS"), ("S1_0002", "Bakta product name"),
            ("S1_0003", "Bakta product name"), ("S1_0005", "AMRFinderPlus"),
        })
        self.assertIn("Effector delivery system", rows[1])
        self.assertEqual(summary["tags"]["P"], ["Bakta product name"])
        self.assertEqual(evidence_record(self.run)["tags"]["A"], ["AMRFinderPlus"])

        tags, hits = self._load()
        self.assertEqual(tags, {"S1_0001": "V", "S1_0002": "M", "S1_0003": "P", "S1_0004": "", "S1_0005": "A", "S1_0006": ""})
        product_rows = hits[hits["DB"] == "Bakta product name"]
        self.assertEqual(set(product_rows["Method"]), {"product name match"})

    def test_gui_tags_follow_the_evidence_sources(self) -> None:
        from annotation_explorer_polars.model.run_data import load_run_tag_databases

        build_evidence(self.run, ["S1"])
        record = {"marker": {"used": [], "not_installed": []}, "compound": {"used": [], "not_installed": []},
                  "evidence": evidence_record(self.run)}
        (self.run / "run_config.json").write_text(json.dumps({"inputs": {"databases": record}}), encoding="utf-8")
        tags = load_run_tag_databases(self.run)
        self.assertEqual(tags["A"].note(), "")
        self.assertEqual(tags["V"].note(), "")
        self.assertEqual(tags["P"].note(), "Bakta product names")
        self.assertEqual(tags["I"].note(), "not installed")  # no compound database: no source

    def test_calls_on_no_cds_and_disrupted_genes(self) -> None:
        from annotation_explorer_polars import data

        (self.run / "amrfinder" / "S1.tsv").write_text(AMRFINDER_HEADER + AMRFINDER_ROWS + (
            "S1\tS1_0006\tcontig_1\t2001\t2300\t+\tabaF\tfosfomycin efflux AbaF\tcore\tAMR\tAMR\tFOSFOMYCIN\t"
            "FOSFOMYCIN\tINTERNAL_STOP\t100.00\t99.00\tWP_2\n"
            # tet(B) cut by the contig end: Bakta called no CDS there; it overlaps S1_0006 by 50 bp
            "S1\tNA\tcontig_1\t2251\t3450\t-\ttet(B)\ttetracycline efflux Tet(B)\tcore\tAMR\tAMR\t"
            "TETRACYCLINE\tTETRACYCLINE\tBLASTX\t100.00\t100.00\tWP_3\n"), encoding="utf-8")
        summary = build_evidence(self.run, ["S1"])
        self.assertEqual(summary["amrfinder"]["unplaced"], {"S1": 1})

        ann, hits = data.build_annotations_from_output(self.run, {"inputs": {**self.run_cfg["inputs"], "marker_sequences": []}})
        tet = ann[ann["Gene"] == "tet(B)"].iloc[0]
        self.assertEqual((tet["AnnotationID"], tet["Start"], tet["Stop"], tet["Strand"], tet["Tags"]),
                         ("AFP0000001", 2251, 3450, "-", "A"))
        self.assertEqual(ann.loc[ann["Locus"] == "S1_0006", "Tags"].iloc[0], "A")
        disrupted = hits[hits["SubjectIdentifier"] == "WP_2"].iloc[0]
        self.assertIn("disrupted: internal stop codon", disrupted["Product"])

    def test_user_database_hits_are_labelled(self) -> None:
        from annotation_explorer_polars import data
        from annotation_explorer_polars.model.run_data import (
            label_user_databases, load_run_tag_databases, load_run_user_databases)

        hits_dir = self.run / "annotations" / "bakta_vs_queries"
        (hits_dir / "bakta_annotations_nt-vs-VFDB_Ab30.tsv").write_text(
            "BKT0000001\tVFG2\t99.0\t99\t0\t0\t1\t99\t1\t99\t1e-50\t200\t99\t99\n", encoding="utf-8")
        build_evidence(self.run, ["S1"])
        db = {"name": "VFDB_Ab30", "category": "VIRULENCE", "path": "/DB/seq/markers/VFDB_Ab30.fasta"}
        inputs = {**self.run_cfg["inputs"], "marker_sequences": [db["path"]],
                  "marker_categories": {"VFDB_Ab30": "VIRULENCE"},
                  "databases": {"marker": {"used": [db], "not_installed": []},
                                "compound": {"used": [], "not_installed": []}, "evidence": evidence_record(self.run)}}
        (self.run / "run_config.json").write_text(json.dumps({"inputs": inputs}), encoding="utf-8")
        ann, hits = data.build_annotations_from_output(self.run, {"inputs": inputs})
        self.assertEqual(set(hits.loc[hits["DB"] == "VFDB_Ab30", "Method"]), {"VFDB_Ab30"})  # the whole name

        user_dbs = load_run_user_databases(self.run)
        ann, hits, _ref = label_user_databases(ann, hits, ann.iloc[0:0], user_dbs)
        self.assertIn("VFDB_Ab30 (U)", set(hits["DB"]))
        self.assertEqual(ann.loc[ann["Locus"] == "S1_0001", "Predictors"].iloc[0], "VFDB_Ab30 (U),Bakta VFDB")
        self.assertIn("VFDB_Ab30 (U)", ann.loc[ann["Locus"] == "S1_0001", "MarkerDisplay"].iloc[0])  # the Marker column
        self.assertIn("VFDB_Ab30 (U)", load_run_tag_databases(self.run)["V"].tooltip())

        # A run made before the step (no evidence record): its databases are built-in ones.
        del inputs["databases"]["evidence"]
        (self.run / "run_config.json").write_text(json.dumps({"inputs": inputs}), encoding="utf-8")
        self.assertEqual(load_run_user_databases(self.run), frozenset())


if __name__ == "__main__":
    unittest.main()
