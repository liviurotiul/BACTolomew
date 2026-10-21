"""P02 / P04 concepts: target models are never reported as acquired AMR, and the metadata
index never guesses.

Qt-free; tiny fixtures only.
"""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pipeline.ResistanceMetadata as RM  # noqa: E402
from pipeline import DiamondWorkflow as DW  # noqa: E402
from pipeline import MarkerDbBuild as MB  # noqa: E402

AMRF_COLUMNS = [
    "allele", "gene_family", "whitelisted_taxa", "product_name", "scope", "type", "subtype",
    "class", "subclass", "refseq_protein_accession", "refseq_nucleotide_accession",
    "genbank_protein_accession", "genbank_nucleotide_accession",
]


def _amrf_row(allele, family, subtype, cls, subclass, prot, gb_prot="", nuc="", gb_nuc=""):
    return {
        "allele": allele, "gene_family": family, "whitelisted_taxa": "", "product_name": family,
        "scope": "core", "type": "AMR", "subtype": subtype, "class": cls, "subclass": subclass,
        "refseq_protein_accession": prot, "refseq_nucleotide_accession": nuc,
        "genbank_protein_accession": gb_prot, "genbank_nucleotide_accession": gb_nuc,
    }


AMRF_ROWS = [
    # adeS: point-mutation references only; the first row used to be attached to every hit.
    _amrf_row("adeS_F12S", "adeS", "POINT", "TETRACYCLINE", "TETRACYCLINE", "WP_031379793.1"),
    _amrf_row("adeS_N125K", "adeS", "POINT", "TETRACYCLINE", "TETRACYCLINE", "WP_031379793.1"),
    _amrf_row("gyrA_G79C", "gyrA", "POINT", "QUINOLONE", "QUINOLONE", "WP_000116450.1"),
    _amrf_row("ompK36_D135DGD", "ompK36", "POINT_DISRUPT", "BETA-LACTAM", "CARBAPENEM", "WP_000000001.1"),
    _amrf_row("blaOXA-23", "blaOXA", "AMR", "BETA-LACTAM", "CARBAPENEM", "WP_001046004.1", nuc="NG_049525.1"),
    # fexA: acquired-gene row plus point rows on the same accession; the acquired row wins.
    _amrf_row("fexA_G33A", "fexA", "POINT", "PHENICOL", "PHENICOL", "WP_015585966.1"),
    _amrf_row("", "fexA", "AMR", "PHENICOL", "CHLORAMPHENICOL/FLORFENICOL", "WP_015585966.1"),
    # Two acquired rows that disagree on the allele: keep the family, leave the allele empty.
    _amrf_row("blaX-1", "blaX", "AMR", "BETA-LACTAM", "BETA-LACTAM", "WP_999999999.1"),
    _amrf_row("blaX-2", "blaX", "AMR", "BETA-LACTAM", "BETA-LACTAM", "WP_999999999.1"),
]

CARD_INDEX_ROWS = [
    # ARO, protein, ARO name, family, drug class
    ("ARO:3001418", "AAV65289.1", "OXA-23", "OXA beta-lactamase", "carbapenem;penicillin beta-lactam"),
    ("ARO:3003817", "AJF82744.1", "Acinetobacter baumannii gyrA conferring resistance to fluoroquinolones",
     "fluoroquinolone resistant gyrA", "fluoroquinolone antibiotic"),
    ("ARO:3000549", "ADM92606.1", "adeS", "RND efflux regulator", "tetracycline antibiotic"),
    # Same protein accession under two homolog AROs with different classes (ambiguous by protein).
    ("ARO:3000001", "AAA00001.1", "NfxB type A", "RND pump", "fluoroquinolone antibiotic"),
    ("ARO:3000002", "AAA00001.1", "NfxB type B", "RND pump", "aminoglycoside antibiotic"),
]

CARD_JSON = {
    "1": {"model_type": "protein homolog model", "ARO_accession": "3001418"},
    "2": {"model_type": "protein variant model", "ARO_accession": "3003817"},
    "3": {"model_type": "protein homolog model", "ARO_accession": "3000549"},
    "4": {"model_type": "protein homolog model", "ARO_accession": "3000001"},
    "5": {"model_type": "protein homolog model", "ARO_accession": "3000002"},
    "6": {"model_type": "protein knockout model", "ARO_accession": "3003809"},
    # Listed as HM in an older CARD FASTA, only an overexpression model in card.json (baeS).
    "7": {"model_type": "protein overexpression model", "ARO_accession": "3000829"},
    "_version": "4.0.1",
}


def _write_amrf_catalog(path: Path) -> Path:
    pd.DataFrame(AMRF_ROWS, columns=AMRF_COLUMNS).to_csv(path, sep="\t", index=False)
    return path


def _write_card(dir_: Path, with_json: bool = True) -> Path:
    rows = [
        {"ARO Accession": a, "Protein Accession": p, "ARO Name": n, "CARD Short Name": n,
         "AMR Gene Family": f, "Drug Class": c, "Resistance Mechanism": "x"}
        for a, p, n, f, c in CARD_INDEX_ROWS
    ]
    path = dir_ / "aro_index.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    if with_json:
        (dir_ / "card.json").write_text(json.dumps(CARD_JSON), encoding="utf-8")
    return path


def _hit(aid, db, product, identity=99.0, bitscore=500.0, category="AMR"):
    return {
        "AnnotationID": aid, "Product": product, "Identity": identity, "Align_len": 300,
        "Mismatches": 0, "Gap_openings": 0, "Q_start": 1, "Q_end": 900, "S_start": 1, "S_end": 300,
        "Evalue": 1e-100, "Bitscore": bitscore, "Q_len": 900, "S_len": 300, "Subject_coverage": 100.0,
        "Method": db.split("_")[0], "MarkerCategory": category, "Category": category, "DB": db,
    }


class _TmpCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        RM._load_amrf_catalog.cache_clear()
        RM._load_card_index.cache_clear()

    def tearDown(self):
        self._tmp.cleanup()

    def meta(self, with_card_json=True):
        card_dir = self.tmp / "card"
        card_dir.mkdir(exist_ok=True)
        return {
            "AMRF_core-AMR": {"source": "AMRF", "files": [str(_write_amrf_catalog(self.tmp / "catalog.tsv"))]},
            "CARD_AMR": {"source": "CARD", "files": [str(_write_card(card_dir, with_card_json))]},
        }


CARD_FASTA = """>HM|gb|AAV65289.1|ARO:3001418|OXA-23 [Acinetobacter baumannii]
MNKYFTCYVVASLFLSGCTVQHNLINETPSQIVQGHNQVIHQYFDEKNTSGVLVIQTDKK
>VM|gb|AJF82744.1|ARO:3003817|Abau_gyrA_FLO [Acinetobacter baumannii]
MSDSLNTQNLHELLEDTPV
>KM|gb|AHB92962.1|ARO:3003809|Abau_OprD_IPM [Acinetobacter baumannii]
MKVMKWSAIALAVSAGSTQ
>OM|gb|X00001.1|ARO:3000999|tetR [Escherichia coli]
MTKLQPNTVIRAALDLLNE
>HM|gb|BAA15934.1|ARO:3000829|baeS [Escherichia coli]
MRIGMRLLLGYFLLVAVAA
"""

AMRF_FASTA = """>WP_001046004.1 MULTISPECIES: carbapenem-hydrolyzing class D beta-lactamase OXA-23 [Bacteria]
MNKYFTCYVVASLFLSGCTV
>WP_031379793.1 MULTISPECIES: two-component sensor histidine kinase AdeS [Acinetobacter]
MKSKLGISKQLFIALTIVN
>WP_015585966.1 MULTISPECIES: chloramphenicol/florfenicol efflux MFS transporter FexA [Staphylococcaceae]
MKKDINLWGLLLGVLSIAL
>WP_000000001.1 porin OmpK36 [Klebsiella]
MKVKVLSLLVPALLVAGAA
"""


class TargetModelsTests(_TmpCase):
    """CARD VM/KM/OM and AMRF POINT/POINT_DISRUPT/AMR-SUSCEPTIBLE are never acquired AMR."""

    def test_target_models_are_dropped_from_db_and_hits(self):
        root = self.tmp / "DB"
        (root / "seq" / "markers").mkdir(parents=True)
        (root / "metadata" / "AMRF").mkdir(parents=True)
        (root / "metadata" / "CARD").mkdir(parents=True)
        (root / "seq" / "markers" / "CARD_AMR.fasta").write_text(CARD_FASTA)
        (root / "seq" / "markers" / "AMRF_core-AMR.fasta").write_text(AMRF_FASTA)
        _write_amrf_catalog(root / "metadata" / "AMRF" / "ReferenceGeneCatalog.txt")
        (root / "metadata" / "CARD" / "card.json").write_text(json.dumps(CARD_JSON))

        MB.build_marker_dbs(root)
        markers = root / "seq" / "markers"
        card = [h.split()[0] for h, _ in MB.iter_fasta_records(markers / "CARD_AMR.fasta")]
        self.assertEqual(card, ["HM|gb|AAV65289.1|ARO:3001418|OXA-23"])
        amrf = [h.split()[0] for h, _ in MB.iter_fasta_records(markers / "AMRF_core-AMR.fasta")]
        self.assertEqual(amrf, ["WP_001046004.1", "WP_015585966.1"])
        # The unfiltered source is archived, so a rebuild starts from it again.
        self.assertTrue((root / "seq" / MB.SOURCES_DIRNAME / "CARD_AMR.unfiltered.fasta.gz").is_file())

        # Hits from a DB built before the fix are dropped by the metadata guard.
        hits = pd.DataFrame([
            _hit("BKT1", "CARD_AMR", "VM|gb|AJF82744.1|ARO:3003817|Abau_gyrA_FLO", identity=99.9),
            _hit("BKT1", "CARD_AMR", "HM|gb|AAV65289.1|ARO:3001418|OXA-23", identity=85.0),
            _hit("BKT2", "CARD_AMR", "KM|gb|AHB92962.1|ARO:3003809|Abau_OprD_IPM"),
            _hit("BKT3", "CARD_AMR", "HM|gb|BAA15934.1|ARO:3000829|baeS"),
            _hit("BKT1", "AMRF_core-AMR", "WP_000116450.1", category="core-AMR"),
            _hit("BKT5", "AMRF_core-AMR", "WP_031379793.1", category="core-AMR"),
            _hit("BKT4", "AMRF_core-AMR", "WP_001046004.1", category="core-AMR"),
        ])
        kept = DW.collapse_diamond_annotations(RM.drop_resistance_target_hits(hits, self.meta()))
        self.assertEqual(
            sorted(zip(kept["AnnotationID"], kept["DB"], kept["Product"])),
            [("BKT1", "CARD_AMR", "HM|gb|AAV65289.1|ARO:3001418|OXA-23"),
             ("BKT4", "AMRF_core-AMR", "WP_001046004.1")],
        )


class MetadataIndexTests(_TmpCase):
    """P04: a hit gets its own reference's metadata, never another catalog row's."""

    def test_index_never_guesses(self):
        idx = RM._load_amrf_catalog(str(_write_amrf_catalog(self.tmp / "c.tsv")))
        self.assertEqual(idx.by_accession["WP_001046004.1"]["ResistanceAlleles"], "blaOXA-23")
        # Acquired row wins over point rows on the same accession.
        self.assertEqual(idx.by_accession["WP_015585966.1"]["ResistanceAntibiotics"], "CHLORAMPHENICOL/FLORFENICOL")
        # Rows that disagree: keep the shared family, leave the allele empty.
        self.assertTrue(pd.isna(idx.by_accession["WP_999999999.1"]["ResistanceAlleles"]))
        self.assertEqual(idx.by_accession["WP_999999999.1"]["ResistanceGeneFamilies"], "blaX")
        card_dir = self.tmp / "card"
        card_dir.mkdir()
        card = RM._load_card_index(str(_write_card(card_dir)))
        self.assertTrue(pd.isna(card.by_protein["AAA00001.1"]["ResistanceAlleles"]))
        self.assertEqual(card.by_aro["ARO:3000002"]["ResistanceAlleles"], "NfxB type B")


if __name__ == "__main__":
    unittest.main()
