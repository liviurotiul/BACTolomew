"""Phylogeny rework, the critical concepts only: the ordered core alignment follows the
ordering genome and maps back to genes, a Gubbins run that did not converge is surfaced
as a warning, the automatic rule picks one Gubbins run or one per lineage, the
merged masked alignment masks each genome only with its own lineage's calls, and the
"Including recombination" branch lengths keep the core tree's topology, root edge and
leaf order."""
from __future__ import annotations

import json
import logging
import tempfile
import unittest
from pathlib import Path

from pipeline.CoreAlignment import (
    build_gene_map,
    file_sha256,
    iter_fasta_records,
    map_alignment_interval,
    parse_embl_header,
    place_from_cells,
    read_gff_index,
    write_ordered_alignment,
)
from pipeline.GubbinsWorkflow import (
    MaskSource,
    merge_masked_alignment,
    not_converged_warning,
    parse_convergence,
    resolve_gubbins_groups,
)
from pipeline.Newick import leaf_names, read_newick, root_split, split_lengths, split_supports
from pipeline.ToolEnv import configured_phylo_settings
from pipeline.TreesIndex import write_trees_index

# Panaroo order: geneA (4 bp), geneB (3 bp), geneC (5 bp), geneD (2 bp, absent from S1).
EMBL = """ID   Genome standard; DNA; PRO; 14 BP.
FT   feature         1..4
FT                   /label=geneA.aln
FT   feature         5..7
FT                   /label=geneB.aln
FT   feature         8..12
FT                   /label=geneC.aln
FT   feature         13..14
FT                   /label=geneD.aln
"""
# S1: geneC first on the longest contig (minus strand), then geneA; geneB on contig 2,
# split in two fragments (";" cell); geneD only refound.
GFF = """##gff-version 3
##sequence-region c2 1 500
##sequence-region c1 1 1000
c1\tx\tCDS\t100\t104\t.\t-\t0\tID=S1_3
c1\tx\tCDS\t300\t303\t.\t+\t0\tID=S1_1
c2\tx\tCDS\t10\t11\t.\t+\t0\tID=S1_2a
c2\tx\tCDS\t40\t41\t.\t+\t0\tID=S1_2b
"""
CELLS = {"geneA": "S1_1", "geneB": "S1_2a;S1_2b", "geneC": "S1_3", "geneD": "S1_refound_9"}
ALIGNMENT = {
    "S1": "ACGT" "GGA" "AACGT" "TT",
    "S2": "ACGA" "GG-" "AACGN" "TC",
}


class OrderedCoreAlignmentTest(unittest.TestCase):
    def test_genes_follow_the_ordering_genome_and_map_back(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "h.embl").write_text(EMBL)
            (tmp / "s1.gff3").write_text(GFF)
            (tmp / "src.aln").write_text("".join(f">{k}\n{v}\n" for k, v in ALIGNMENT.items()))
            blocks = parse_embl_header(tmp / "h.embl")
            placements, counts = place_from_cells([b.cluster for b in blocks], CELLS, read_gff_index(tmp / "s1.gff3"))
            self.assertEqual((counts.placed, counts.fragmented, counts.refound_only), (3, 1, 1))
            rows = build_gene_map(blocks, placements, {"geneA", "geneC"})
            # Longest contig first (c1: geneC at 100, geneA at 300), then c2, then unplaced.
            self.assertEqual([r.cluster for r in rows], ["geneC", "geneA", "geneB", "geneD"])
            self.assertEqual([r.reversed for r in rows], [True, False, False, False])
            self.assertEqual((rows[2].ordering_start, rows[2].ordering_end), (10, 41))
            self.assertFalse(rows[3].placed)

            write_ordered_alignment(tmp / "src.aln", rows, tmp / "out.aln", expected_length=14)
            out = {name: seq.decode() for name, seq in iter_fasta_records(tmp / "out.aln")}
            comp = str.maketrans("ACGT", "TGCA")
            for name, src in ALIGNMENT.items():
                for row in rows:
                    panaroo = src[row.panaroo_aln_start - 1: row.panaroo_aln_end]
                    if row.reversed:
                        panaroo = panaroo.translate(comp)[::-1]
                    self.assertEqual(out[name][row.aln_start - 1: row.aln_end], panaroo, (name, row.cluster))

            # Columns 2..6 of the ordered alignment: geneC (reversed, gene columns 4..1)
            # and the first column of geneA.
            hits = map_alignment_interval(rows, 2, 6)
            self.assertEqual([h["cluster"] for h in hits], ["geneC", "geneA"])
            self.assertEqual((hits[0]["gene_start"], hits[0]["gene_end"]), (1, 4))
            self.assertEqual((hits[0]["ordering_start"], hits[0]["ordering_end"]), (101, 104))
            self.assertEqual((hits[1]["gene_start"], hits[1]["ordering_start"]), (1, 300))


class GubbinsConvergenceTest(unittest.TestCase):
    def test_non_convergence_is_a_warning(self) -> None:
        stuck = parse_convergence("...\nChecking for convergence...\nMaximum number of iterations (20) reached.\n")
        self.assertEqual((stuck["converged"], stuck["iterations_run"]), (False, 20))
        warning = not_converged_warning(stuck, 20)
        self.assertEqual(warning["code"], "gubbins_not_converged")

        done = parse_convergence("Convergence after 7 iterations: Tree observed before.\n")
        self.assertEqual((done["converged"], done["iterations_run"]), (True, 7))
        self.assertIsNone(not_converged_warning(done, 20))
        # A log that says neither is treated as not converged.
        self.assertIsNotNone(not_converged_warning(parse_convergence(""), 20))



class GubbinsModeTest(unittest.TestCase):
    def test_auto_picks_one_run_or_one_per_lineage(self) -> None:
        settings = {**configured_phylo_settings(), "gubbins_mode": "auto", "gubbins_min_lineage_size": 5}
        samples = [f"S{i}" for i in range(12)]
        one = resolve_gubbins_groups(samples, settings, {s: "1" for s in samples})
        self.assertEqual((one.resolved_mode, [g.id for g in one.groups]), ("single", ["all"]))

        lineages = {s: ("1" if i < 4 else "2") for i, s in enumerate(samples)}
        per = resolve_gubbins_groups(samples, settings, lineages)
        self.assertEqual(per.resolved_mode, "per_lineage")
        self.assertEqual([(g.id, len(g.members), g.skip_reason) for g in per.groups],
                         [("1", 4, "below_min_lineage_size"), ("2", 8, None)])

        fallback = resolve_gubbins_groups(samples, settings, None, "fastGEAR exited 1")
        self.assertEqual((fallback.resolved_mode, fallback.mode_reason), ("single", "fastgear_unavailable"))
        self.assertEqual([w["code"] for w in fallback.warnings], ["gubbins_lineages_unavailable"])
        with self.assertRaises(RuntimeError):
            resolve_gubbins_groups(samples, {**settings, "gubbins_mode": "per_lineage"}, None)


class MaskedAlignmentMergeTest(unittest.TestCase):
    def test_each_genome_is_masked_only_by_its_own_lineage(self) -> None:
        rows = {"A1": "ACGTACGTAC", "A2": "ACGTACGTAA", "B1": "TTGTACGNAC", "B2": "TTGTACGTAC", "C1": "ACGAACGTAC"}
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "ordered.aln").write_text("".join(f">{k}\n{v}\n" for k, v in rows.items()))
            # Lineage A: block 3..5 in A1. Lineage B: a block at the same columns in B2 and
            # 6..9 in B1 (column 8 already N). B2 was left out by Gubbins in this test.
            gff_a = '##gff-version 3\nSEQUENCE\tGUBBINS\tCDS\t3\t5\t0\t.\t0\ttaxa=" A1";\n'
            gff_b = (
                '##gff-version 3\nSEQUENCE\tGUBBINS\tCDS\t3\t5\t0\t.\t0\ttaxa=" B2";\n'
                'SEQUENCE\tGUBBINS\tCDS\t6\t9\t0\t.\t0\ttaxa=" B1";\n'
            )
            (tmp / "a.gff").write_text(gff_a)
            (tmp / "b.gff").write_text(gff_b)
            masked_a = {"A1": "ACNNNCGTAC", "A2": rows["A2"]}
            masked_b = {"B1": "TTGTANNNNC", "B2": "TTNNNCGTAC"}
            (tmp / "a.masked").write_text("".join(f">{k}\n{v}\n" for k, v in masked_a.items()))
            (tmp / "b.masked").write_text("".join(f">{k}\n{v}\n" for k, v in masked_b.items()))
            sources = [
                MaskSource("A", frozenset(masked_a), frozenset(masked_a), tmp / "a.masked", tmp / "a.gff"),
                MaskSource("B", frozenset(masked_b), frozenset({"B1"}), tmp / "b.masked", tmp / "b.gff"),
            ]
            # B2's masked row carries a block although Gubbins dropped it: refused.
            with self.assertRaises(RuntimeError):
                merge_masked_alignment(tmp / "ordered.aln", sources, {"C1": ("C", "lineage_skipped")}, tmp / "out.aln")
            masked_b["B2"] = rows["B2"]
            (tmp / "b.masked").write_text("".join(f">{k}\n{v}\n" for k, v in masked_b.items()))
            result = merge_masked_alignment(
                tmp / "ordered.aln", sources, {"C1": ("C", "lineage_skipped")}, tmp / "out.aln"
            )
            out = {name: seq.decode() for name, seq in iter_fasta_records(tmp / "out.aln")}
            self.assertEqual(list(out), list(rows))
            self.assertEqual(out["A1"], "ACNNNCGTAC")
            self.assertEqual(out["A2"], rows["A2"])
            self.assertEqual(out["B1"], "TTGTANNNNC")
            self.assertEqual(out["B2"], rows["B2"])  # dropped by Gubbins: original row
            self.assertEqual(out["C1"], rows["C1"])  # lineage not run: original row
            by_sample = {r["sample"]: r for r in result["rows"]}
            self.assertEqual(by_sample["B1"]["masked_bases"], 3)  # column 8 was already N
            self.assertEqual(by_sample["B2"]["unmasked_reason"], "dropped_by_gubbins")
            self.assertEqual(by_sample["C1"]["unmasked_reason"], "lineage_skipped")

            # A masked row that differs from its own blocks fails the merge.
            masked_a["A2"] = "ACGTNCGTAA"
            (tmp / "a.masked").write_text("".join(f">{k}\n{v}\n" for k, v in masked_a.items()))
            with self.assertRaises(RuntimeError):
                merge_masked_alignment(tmp / "ordered.aln", sources, {"C1": ("C", "lineage_skipped")}, tmp / "out.aln")


class RecombinationBranchLengthsTest(unittest.TestCase):
    def test_refit_keeps_topology_root_edge_and_leaf_order(self) -> None:
        # The core tree is midpoint rooted on the edge above (A,B). The refit (IQ-TREE -te
        # output: other child order, no support) has a very long D, so a fresh midpoint
        # root would sit on D's edge instead.
        core = "((A:0.1,B:0.1)90:0.5,(C:0.1,D:0.1)80:0.05,(E:0.1,F:0.2)70:0.05);\n"
        refit = "((D:3.0,C:0.1):0.05,(B:0.1,A:0.1):0.02,(F:0.2,E:0.1):0.05);\n"
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp)

            def write(rel: str, text: str) -> Path:
                path = run / rel
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text)
                return path

            raw = write("trees/core_recomb_iqtree/core_recomb.treefile", core)
            write("trees/core_recomb_refit/core_recomb_with_recombination.treefile", refit)
            write("core_alignment/ordering.json", json.dumps({"output": {"sha256": "ordered"}}))
            write(
                "gubbins/summary.json",
                json.dumps({"schema_version": 3, "resolved_mode": "per_lineage", "masked_alignment": {"sha256": "masked"}}),
            )
            core_stamp = {
                "id": "core_recomb",
                "raw_file": "trees/core_recomb_iqtree/core_recomb.treefile",
                "input_sha256": "masked",
                "alignment": "gubbins/core_gene_alignment.masked.aln",
                "model": "K3Pu+F+I+G4",
                "branch_length_unit": "substitutions/site",
                "recombination_filtered": True,
            }
            write("trees/core_recomb.json", json.dumps(core_stamp))
            refit_stamp = {
                "status": "built",
                "input_tree_sha256": file_sha256(raw),
                "alignment": "core_alignment/core_gene_alignment.ordered.aln",
                "alignment_sha256": "ordered",
                "model": "K3Pu+F+I+G4",
            }
            write("trees/core_recomb_refit/run.json", json.dumps(refit_stamp))
            settings = configured_phylo_settings({"outgroup": ""})
            payload, _ = write_trees_index(run, list("ABCDEF"), settings, logging.getLogger("test"))

            sets = payload["trees"][0]["branch_length_sets"]
            self.assertEqual([s["id"] for s in sets], ["mutations", "with_recombination"])
            shown = read_newick(run / "trees" / sets[0]["file"])
            other = read_newick(run / "trees" / sets[1]["file"])
            (side, share), (other_side, other_share) = root_split(shown), root_split(other)
            self.assertIn(side, (frozenset("AB"), frozenset("CDEF")))
            self.assertEqual(other_side, side)
            self.assertAlmostEqual(other_share, share)
            self.assertEqual(leaf_names(other), leaf_names(shown))
            self.assertEqual(split_supports(other), split_supports(shown))
            self.assertEqual(split_lengths(other)[frozenset("D")], 3.0)

            # A refit of another tree is not listed and its file goes.
            write("trees/core_recomb_refit/run.json", json.dumps({**refit_stamp, "input_tree_sha256": "old"}))
            payload, _ = write_trees_index(run, list("ABCDEF"), settings, logging.getLogger("test"))
            self.assertNotIn("branch_length_sets", payload["trees"][0])
            self.assertFalse((run / "trees" / sets[1]["file"]).exists())


if __name__ == "__main__":
    unittest.main()
