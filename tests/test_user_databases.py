import gzip
import tempfile
import unittest
from pathlib import Path

from pipeline import UserDatabases as UD
from pipeline.DbCatalog import load_db_catalog_index


class UserDatabasesTests(unittest.TestCase):
    def test_add_validate_remove(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            root = tmp / "DB"
            protein = tmp / "vir.faa"
            protein.write_text(">p1 toxin\nMKLVVLLAGSTRQW\n>p2\nMSTNPKPQRKTKRNTNRRPQDVKFPGG\n")
            nucleotide = tmp / "nuc.fasta"
            nucleotide.write_text(">n1\nATGAAACGCATTAGCACCACCATTACCGGTAACGGTGCGTGA\n")
            duplicate = tmp / "dup.faa"
            duplicate.write_text(">p1\nMKLV\n>p1\nMSTN\n")

            with self.assertRaisesRegex(UD.DatabaseError, "nucleotide"):
                UD.add_database(nucleotide, "nuc", "VIRULENCE", db_root=root)
            with self.assertRaisesRegex(UD.DatabaseError, "more than once"):
                UD.add_database(duplicate, "dup", "VIRULENCE", db_root=root)
            with self.assertRaises(UD.DatabaseError):
                UD.add_database(protein, "bad name", "VIRULENCE", db_root=root)
            dotted = tmp / "dot.faa"
            dotted.write_text(">d1\nMKV.LA\n")  # DIAMOND makedb refuses '.'
            with self.assertRaisesRegex(UD.DatabaseError, "no amino acids"):
                UD.add_database(dotted, "dot", "AMR", db_root=root)

            db = UD.add_database(protein, "MyVir", "virulence", note="CC BY", db_root=root)
            self.assertEqual((db.sequences, db.category, db.status), (2, "VIRULENCE", "ok"))
            self.assertTrue((root / "seq/markers/MyVir.fasta").is_file())
            with self.assertRaisesRegex(UD.DatabaseError, "exists already"):
                UD.add_database(protein, "myvir", "AMR", db_root=root)

            # The pipeline's catalog reader finds it with its category.
            index = load_db_catalog_index(UD.catalog_path("marker", root))
            self.assertEqual(index.category_map([db.path]), {"MyVir": "VIRULENCE"})

            # A gzipped FASTA with a UTF-8 BOM is unpacked into the DB folder; IUPAC codes
            # count as nucleotides; the evidence sources' names are taken.
            packed = tmp / "vir2.faa.gz"
            packed.write_bytes(gzip.compress(b"\xef\xbb\xbf>q1\nMKLVVLLAGSTRQW\n"))
            db2 = UD.add_database(packed, "Packed", "AMR", db_root=root)
            self.assertEqual((root / "seq/markers/Packed.fasta").read_bytes(), b">q1\nMKLVVLLAGSTRQW\n")
            self.assertEqual(db2.sequences, 1)
            iupac = tmp / "iupac.fasta"
            iupac.write_text(">n2\nATGRYKMSWNNACGTBDHVATGCATGCA\n")
            with self.assertRaisesRegex(UD.DatabaseError, "nucleotide"):
                UD.add_database(iupac, "iupac", "AMR", db_root=root)
            for taken in ("AMRFinderPlus", "MOB-suite", "Bakta_VFDB"):
                with self.assertRaisesRegex(UD.DatabaseError, "evidence sources"):
                    UD.add_database(protein, taken, "AMR", db_root=root)
            with self.assertRaisesRegex(UD.DatabaseError, "-vs-"):
                UD.add_database(protein, "a-vs-b", "AMR", db_root=root)
            UD.remove_database("Packed", root)

            moved = UD.remove_database("MyVir", root)
            self.assertTrue(moved.is_file() and moved.parent.parent == root / "_removed")
            self.assertFalse(db.path.exists())
            self.assertEqual(UD.list_databases(root), [])


if __name__ == "__main__":
    unittest.main()
