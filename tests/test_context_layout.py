from __future__ import annotations

import math
import unittest

from PyQt6 import QtGui

from annotation_explorer_polars.context_layout import ContextLayoutEngine, link_payload
from annotation_explorer_polars.context_model import (
    ContextRenderModel,
    annotation_category,
    family_color,
    homology_key_for_gene,
)


def gene(start: int, end: int, name: str, strand: str = "+", **extra):
    value = {
        "start": start,
        "end": end,
        "strand": strand,
        "label": name,
        "row_data": {"Gene": name, "PanarooCluster": f"cluster_{name}"},
    }
    value.update(extra)
    return value


class ContextModelTests(unittest.TestCase):
    def test_homology_priority_and_stable_family_colour(self) -> None:
        item = gene(10, 90, "blaOXA-23")
        item["row_data"]["Gene"] = "different-name"
        self.assertEqual(homology_key_for_gene(item), "cluster:cluster blaoxa 23")
        key = homology_key_for_gene(item)
        self.assertEqual(family_color(key), family_color(key))

    def test_unknown_fields_and_gene_order_are_preserved(self) -> None:
        rows = [{
            "sample": "S1",
            "unknown_row_field": {"kept": True},
            "genes": [gene(900, 990, "late", custom="first"), gene(100, 190, "early", custom="second")],
        }]
        model = ContextRenderModel.from_rows(rows)
        self.assertEqual([value.index for value in model.rows[0].genes], [0, 1])
        self.assertEqual(model.rows[0].source["unknown_row_field"], {"kept": True})
        self.assertEqual(model.rows[0].genes[0].source["custom"], "first")

    def test_center_bp_is_the_exact_alignment_anchor(self) -> None:
        rows = [{"sample": "S1", "center_bp": 175, "genes": [gene(100, 200, "focus", is_center=True)]}]
        layout = ContextLayoutEngine(px_per_bp=1.0).layout(ContextRenderModel.from_rows(rows))
        placed = layout.rows[0].genes[0]
        self.assertEqual((placed.left, placed.right), (-75.0, 25.0))

    def test_explicit_category_wins_over_generic_resistance_source(self) -> None:
        for explicit in ("BMRG", "VIRULENCE", "MOBILE"):
            with self.subTest(explicit=explicit):
                self.assertEqual(
                    annotation_category({
                        "ConsensusCategory": explicit,
                        "ResistanceMetadataSource": "CARD resistance",
                    }),
                    explicit,
                )

    def test_marker_search_bag_uses_stable_informative_token(self) -> None:
        first = {
            "label": "gene",
            "row_data": {"MarkerSearchText": "beta-lactamase family, CARD, BLAST, AMR"},
        }
        second = {
            "label": "gene",
            "row_data": {"MarkerSearchText": "beta-lactamase family, AMRFinder, exact hit"},
        }
        self.assertEqual(homology_key_for_gene(first), "marker:beta lactamase family")
        self.assertEqual(homology_key_for_gene(first), homology_key_for_gene(second))

    def test_useful_wrapper_label_wins_over_generated_unknown_id(self) -> None:
        model = ContextRenderModel.from_rows([{
            "genes": [{
                "start": 1,
                "end": 100,
                "label": "membrane-associated enzyme",
                "row_data": {
                    "Gene": "",
                    "GeneID": "UNK:abcdef123456",
                    "Product": "membrane-associated enzyme",
                },
            }]
        }])
        self.assertEqual(model.rows[0].genes[0].label, "membrane-associated enzyme")

    def test_qcolor_and_gene_level_molecule_metadata_are_normalized(self) -> None:
        model = ContextRenderModel.from_rows([{
            "genes": [{
                "start": 1,
                "end": 100,
                "label": "abc",
                "color": QtGui.QColor("#12ab34"),
                "row_data": {"Gene": "abc", "MobSuiteIsPlasmid": "true"},
            }]
        }])
        self.assertEqual(model.rows[0].genes[0].annotation_fill, "#12ab34")
        self.assertEqual(model.rows[0].molecule, "plasmid")

    def test_non_finite_coordinates_and_row_metrics_are_rejected(self) -> None:
        model = ContextRenderModel.from_rows([{
            "context_similarity": float("inf"),
            "center_bp": float("nan"),
            "genes": [
                gene(0, 100, "valid"),
                gene(float("nan"), 200, "nan"),
                gene(300, float("inf"), "infinite"),
            ],
        }])
        self.assertEqual([value.label for value in model.rows[0].genes], ["valid"])
        self.assertIsNone(model.rows[0].similarity)
        self.assertIsNone(model.rows[0].center_bp)
        layout = ContextLayoutEngine().layout(model)
        self.assertTrue(all(math.isfinite(value) for value in (
            layout.min_x,
            layout.max_x,
            layout.min_y,
            layout.max_y,
        )))


class ContextLayoutTests(unittest.TestCase):
    def test_reverse_match_is_flipped_around_focal_coordinate(self) -> None:
        rows = [
            {"sample": "reference", "is_source": True, "center_bp": 150, "genes": [gene(100, 200, "focus", "+", is_center=True)]},
            {"sample": "reverse", "context_orientation": "reverse", "center_bp": 1150, "genes": [gene(1100, 1200, "focus", "-", is_center=True)]},
        ]
        layout = ContextLayoutEngine(px_per_bp=1.0).layout(ContextRenderModel.from_rows(rows))
        self.assertFalse(layout.rows[0].flipped)
        self.assertTrue(layout.rows[1].flipped)
        self.assertEqual(layout.rows[1].genes[0].display_strand, "+")
        self.assertEqual(layout.rows[1].genes[0].center_x, 0.0)

    def test_gene_widths_remain_proportional_to_base_pairs(self) -> None:
        rows = [{"genes": [gene(0, 100, "short"), gene(200, 600, "long")]}]
        layout = ContextLayoutEngine(px_per_bp=0.1).layout(ContextRenderModel.from_rows(rows))
        widths = [placed.width for placed in layout.rows[0].genes]
        self.assertAlmostEqual(widths[1] / widths[0], 4.0)

    def test_repeated_families_are_matched_one_to_one(self) -> None:
        rows = []
        for sample, offset in (("a", 0), ("b", 1000)):
            rows.append({
                "sample": sample,
                "center_bp": offset + 500,
                "genes": [
                    gene(offset + 100, offset + 200, "dup"),
                    gene(offset + 300, offset + 400, "dup"),
                    gene(offset + 500, offset + 600, "single", is_center=True),
                ],
            })
        layout = ContextLayoutEngine().layout(ContextRenderModel.from_rows(rows))
        links = layout.pairs[0].links
        self.assertEqual(len(links), 3)
        duplicate_links = [link for link in links if "dup" in link.family_key]
        self.assertEqual(len(duplicate_links), 2)
        self.assertEqual(len({link.lower_gene.visual_index for link in duplicate_links}), 2)

    def test_unequal_repeats_choose_the_nearest_collinear_copy(self) -> None:
        upper = {
            "center_bp": 1000,
            "genes": [
                gene(100, 200, "dup"),
                gene(950, 1050, "dup", is_center=True),
            ],
        }
        lower = {
            "center_bp": 1000,
            "genes": [gene(950, 1050, "dup", is_center=True)],
        }
        layout = ContextLayoutEngine(px_per_bp=1.0).layout(
            ContextRenderModel.from_rows([upper, lower])
        )
        self.assertEqual(len(layout.pairs[0].links), 1)
        self.assertEqual(layout.pairs[0].links[0].upper_gene.visual_index, 1)

    def test_weak_product_text_does_not_create_default_links(self) -> None:
        weak = {"start": 0, "end": 100, "strand": "+", "label": "hypothetical protein", "row_data": {"Product": "membrane-associated protein"}}
        layout = ContextLayoutEngine().layout(ContextRenderModel.from_rows([{"genes": [weak]}, {"genes": [dict(weak)]}]))
        self.assertEqual(layout.pairs[0].links, ())

    def test_link_payload_retains_dialog_schema(self) -> None:
        rows = [{"genes": [gene(0, 100, "focus", is_center=True)]}, {"genes": [gene(0, 100, "focus", is_center=True)]}]
        layout = ContextLayoutEngine().layout(ContextRenderModel.from_rows(rows))
        payload = link_payload(layout.pairs[0].links[0])
        self.assertEqual(payload["matched_gene_count"], 1)
        self.assertIn("evidence_summary", payload)
        self.assertEqual(len(payload["matches"]), 1)
        self.assertEqual(payload["matches"][0]["upper_row"], 0)

    def test_short_rows_fold_expand_or_omit(self) -> None:
        rows = [
            {"sample": "ref", "is_source": True, "center_bp": 150,
             "genes": [gene(0, 100, "left"), gene(100, 200, "focus", is_center=True)]},
            {"sample": "short_a", "molecule": "plasmid", "center_bp": 150,
             "genes": [gene(100, 200, "focus", is_center=True)]},
            {"sample": "full", "molecule": "chromosome", "center_bp": 150,
             "genes": [gene(0, 100, "left"), gene(100, 200, "focus", is_center=True)]},
            {"sample": "short_b", "molecule": "chromosome", "collapsed_track_count": 3,
             "center_bp": 150, "genes": [gene(100, 200, "focus", is_center=True)]},
        ]
        model = ContextRenderModel.from_rows(rows)
        self.assertEqual([row.is_short for row in model.rows], [False, True, False, True])
        engine = ContextLayoutEngine()
        inline = engine.layout(model)
        self.assertEqual(len(inline.rows), 4)
        self.assertIsNone(inline.fold)

        folded = engine.layout(model, short_rows="folded")
        self.assertEqual([row.row.sample for row in folded.rows], ["ref", "full"])
        self.assertEqual(folded.fold.contig_count, 4)
        self.assertEqual(folded.fold.molecule_counts(), {"plasmid": 1, "chromosome": 3})
        self.assertEqual(folded.max_y, folded.fold.center_y + folded.fold.height / 2.0 + engine.bottom_margin)

        expanded = engine.layout(model, short_rows="expanded")
        thin = [row for row in expanded.rows if row.short]
        self.assertEqual([row.row.sample for row in thin], ["short_a", "short_b"])
        self.assertTrue(all(row.pitch == engine.short_pitch for row in thin))
        self.assertTrue(all(row.center_y > expanded.fold.center_y for row in thin))
        self.assertEqual(expanded.pairs[-1].links, ())
        self.assertEqual(expanded.row_at(thin[0].center_y), thin[0])
        self.assertIsNone(expanded.row_at(expanded.fold.center_y))

        omitted = engine.layout(model, short_rows="omitted")
        self.assertEqual(len(omitted.rows), 2)
        self.assertIsNone(omitted.fold)

    def test_reference_row_is_never_short(self) -> None:
        model = ContextRenderModel.from_rows([
            {"sample": "ref", "is_source": True, "genes": [gene(0, 100, "focus", is_center=True)]},
        ])
        self.assertFalse(model.rows[0].is_short)
        self.assertEqual(len(ContextLayoutEngine().layout(model, short_rows="folded").rows), 1)

    def test_comparison_reference_wins_over_source(self) -> None:
        model = ContextRenderModel.from_rows([
            {"sample": "source", "is_source": True, "is_comparison_reference": False, "genes": [gene(0, 100, "a")]},
            {"sample": "chosen", "is_comparison_reference": True, "genes": [gene(0, 100, "a")]},
        ])
        self.assertEqual([row.is_reference for row in model.rows], [False, True])

    def test_split_molecules_adds_labelled_dividers(self) -> None:
        def row(sample, molecule, **extra):
            value = {"sample": sample, "molecule": molecule, "center_bp": 150,
                     "genes": [gene(0, 100, "left"), gene(100, 200, "focus", is_center=True)]}
            value.update(extra)
            return value

        rows = [
            row("ref", "chromosome", is_source=True),
            row("p1", "plasmid"),
            row("p2", "plasmid", collapsed_track_count=2),
            row("c1", "chromosome"),
            {"sample": "short", "molecule": "plasmid", "center_bp": 150,
             "genes": [gene(100, 200, "focus", is_center=True)]},
        ]
        engine = ContextLayoutEngine()
        layout = engine.layout(ContextRenderModel.from_rows(rows), short_rows="folded", split_molecules=True)
        self.assertEqual(
            [divider.text for divider in layout.dividers],
            ["Plasmid contigs, 3 (+1 folded)", "Chromosome contigs, 1"],
        )
        self.assertEqual(layout.dividers[0].folded, 1)
        # Each divider sits between the blocks it separates.
        ref, p1, _p2, c1 = layout.rows
        self.assertLess(ref.center_y, layout.dividers[0].center_y)
        self.assertLess(layout.dividers[0].center_y, p1.center_y)
        self.assertLess(layout.dividers[1].center_y, c1.center_y)
        self.assertIsNone(layout.row_at(layout.dividers[0].center_y))

    def test_contig_ends_and_partial_come_from_the_rows(self) -> None:
        model = ContextRenderModel.from_rows([
            {"contig_ends_low": True, "contig_ends_high": False, "context_partial": True,
             "genes": [gene(0, 100, "a")]},
            {"genes": [gene(0, 100, "a")]},
        ])
        first, second = model.rows
        self.assertTrue(first.ends_low and first.partial)
        self.assertFalse(first.ends_high)
        self.assertFalse(second.ends_low or second.ends_high or second.partial)
        layout = ContextLayoutEngine().layout(model)
        self.assertTrue(layout.rows[0].cap_left)
        self.assertFalse(layout.rows[0].cap_right)

    def test_molecule_without_any_information_is_unknown(self) -> None:
        model = ContextRenderModel.from_rows([
            {"genes": [gene(0, 100, "a")]},
            {"molecule": "unknown", "is_plasmid": False, "genes": [gene(0, 100, "a")]},
        ])
        self.assertEqual([row.molecule for row in model.rows], ["unknown", "unknown"])

    def test_merged_rows_count_their_members_by_molecule(self) -> None:
        focal_only = [gene(100, 200, "focus", is_center=True)]
        rows = [
            {"sample": "ref", "is_source": True, "molecule": "chromosome", "center_bp": 150,
             "genes": [gene(0, 100, "left"), gene(100, 200, "focus", is_center=True)]},
            {"sample": "merged", "molecule": "mixed", "center_bp": 150, "collapsed_track_count": 5,
             "collapsed_molecule_counts": {"chromosome": 3, "plasmid": 2}, "genes": focal_only},
        ]
        layout = ContextLayoutEngine().layout(
            ContextRenderModel.from_rows(rows), short_rows="folded", split_molecules=True
        )
        self.assertEqual(layout.fold.molecule_counts(), {"chromosome": 3, "plasmid": 2})
        model = ContextRenderModel.from_rows(rows + [
            {"sample": "p1", "molecule": "plasmid", "center_bp": 150,
             "genes": [gene(0, 100, "left"), gene(100, 200, "focus", is_center=True)]},
        ])
        layout = ContextLayoutEngine().layout(model, short_rows="folded", split_molecules=True)
        self.assertEqual([divider.text for divider in layout.dividers], ["Plasmid contigs, 1 (+2 folded)"])

    def test_rows_without_a_focal_gene_start_at_their_first_gene(self) -> None:
        rows = [
            {"sample": "long", "genes": [gene(0, 100, "repA"), gene(200, 300, "traA"), gene(400, 500, "parA")]},
            {"sample": "short", "genes": [gene(5000, 5100, "repA"), gene(5200, 5300, "traA")]},
        ]
        layout = ContextLayoutEngine(px_per_bp=1.0).layout(ContextRenderModel.from_rows(rows))
        self.assertEqual([row.genes[0].left for row in layout.rows], [0.0, 0.0])


if __name__ == "__main__":
    unittest.main()
