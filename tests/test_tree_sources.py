"""The run's trees (PHYLOGENY_REWORK_2026-10-10.md section 4): every tree
consumer draws the recombination-filtered core tree when the run has one, an
old run falls back to its old tree files, and a Gubbins run that did not
converge is said in the Phylogeny header."""
from __future__ import annotations

import json
import shutil
import tempfile
import unittest
from pathlib import Path

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

from PyQt6 import QtWidgets

from annotation_explorer_polars.model.pan_data import PanData
from annotation_explorer_polars.panaroo_core_accessory import LENGTHS_TOPOLOGY, PanCoreAccessoryTab

# Rooted at the midpoint by the pipeline; UFBoot values as internal labels.
CORE_RECOMB = "((A:0.1,B:0.1)97:0.2,(C:0.1,(D:0.1,E:0.1)72:0.1):0.2);"
SNIPPY = "((A:1,C:1):1,(B:1,(D:1,E:1):1):1);"


def _entry(tree_id: str, label: str, *, filtered: bool, support: str = "none") -> dict:
    return {
        "id": tree_id,
        "label": label,
        "file": f"{tree_id}.nwk",
        "support": {"type": support, "replicates": 1000 if support != "none" else 0, "threshold": 95 if support == "ufboot" else None},
        "rooting": {"method": "midpoint", "outgroup": None, "note": None},
        "recombination_filtered": filtered,
        "branch_length_unit": "substitutions/site",
        "samples_missing": [] if filtered else ["F"],
    }


class TreeSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def setUp(self) -> None:
        self.run = Path(tempfile.mkdtemp(prefix="bactolomew_trees_"))
        self.addCleanup(shutil.rmtree, self.run, True)
        (self.run / "trees").mkdir()
        (self.run / "gubbins").mkdir()
        (self.run / "gubbins" / "core_gene.final_tree.tre").write_text(SNIPPY, encoding="utf-8")
        (self.run / "trees" / "snippy_raxml.tree").write_text(SNIPPY, encoding="utf-8")
        (self.run / "trees" / "core_gene_tree.nwk").write_text(SNIPPY, encoding="utf-8")

    def _write_index(self, default: str, converged: bool = True) -> None:
        (self.run / "trees" / "core_recomb.nwk").write_text(CORE_RECOMB, encoding="utf-8")
        (self.run / "trees" / "snippy_core.nwk").write_text(SNIPPY, encoding="utf-8")
        index = {
            "schema_version": 1,
            "default_tree": default,
            "gubbins": {"ran": True, "converged": converged, "iterations_run": 20, "max_iterations": 20},
            "trees": [
                _entry("snippy_core", "Whole-genome SNPs (Snippy), not recombination-corrected", filtered=False),
                _entry("core_recomb", "Core genome, recombination-filtered", filtered=True, support="ufboot"),
            ],
        }
        (self.run / "trees" / "trees.json").write_text(json.dumps(index), encoding="utf-8")

    def test_the_recombination_filtered_tree_is_the_default(self) -> None:
        # An old run: the Gubbins final tree before the unfiltered RAxML trees.
        legacy = PanData(self.run).tree_index()
        self.assertFalse(legacy.from_index)
        self.assertEqual(legacy.default_id, "gubbins_final")
        self.assertEqual([s.id for s in legacy.sources], ["gubbins_final", "snippy_raxml", "core_gene_raxml"])
        # trees.json is the only source once present; a default naming no
        # listed tree falls back to the recombination-filtered one.
        self._write_index(default="missing")
        index = PanData(self.run).tree_index()
        self.assertTrue(index.from_index)
        self.assertEqual(index.default_id, "core_recomb")
        self.assertEqual([s.id for s in index.sources], ["snippy_core", "core_recomb"])
        self.assertEqual(PanData(self.run).tree_path(), self.run / "trees" / "core_recomb.nwk")

    def test_the_header_says_the_root_the_missing_genomes_and_non_convergence(self) -> None:
        self._write_index(default="core_recomb", converged=False)
        tab = PanCoreAccessoryTab(self.run, None, None, single_view="tree", pan_data=PanData(self.run))
        self.addCleanup(tab.close)
        self.assertEqual(tab.tree_source_id(), "core_recomb")
        self.assertEqual(tab.tree_root_mode(), "file")
        self.assertEqual(tab.tree_root_text(), "Root: midpoint")
        self.assertFalse(tab._gubbins_chip.isHidden())
        self.assertTrue(tab._tree_missing_label.isHidden())
        # Support stays with its edge whatever root is drawn, in the hover only.
        tree = tab.tree
        self.assertEqual(tree.support_text(["A", "B"]), "UFBoot 97, well supported at ≥ 95")
        tab.set_tree_root_mode("midpoint")
        self.assertEqual(tree.support_text(["D", "E"]), "UFBoot 72, not well supported (< 95)")
        # The Snippy tree is not recombination-filtered: no chip, its missing genome said.
        tab.select_tree_source("snippy_core")
        self.assertTrue(tab._gubbins_chip.isHidden())
        self.assertEqual(tab._tree_missing_label.text(), "1 genome not in this tree")
        self.assertEqual(tree.support_text(["D", "E"]), "")
        self.assertEqual(tab.tree_display_options()["source"], "snippy_core")

    def test_a_saved_branch_scale_and_lengths_restore_into_the_settings_column(self) -> None:
        """Panel states saved with the old toolbar (Branch scale, the Branch
        lengths switch, True lengths, Root, Order) restore into the Settings
        column's controls, and are saved again under the same keys."""
        self._write_index(default="core_recomb")
        (self.run / "trees" / "core_recomb_rec.nwk").write_text(CORE_RECOMB.replace("0.1", "0.3"), encoding="utf-8")
        index = json.loads((self.run / "trees" / "trees.json").read_text(encoding="utf-8"))
        index["trees"][1]["branch_length_sets"] = [
            {"id": "mutations", "label": "Mutations only", "file": "core_recomb.nwk", "default": True},
            {"id": "with_recombination", "label": "Including recombination", "file": "core_recomb_rec.nwk"},
        ]
        (self.run / "trees" / "trees.json").write_text(json.dumps(index), encoding="utf-8")
        tab = PanCoreAccessoryTab(
            self.run, None, None, single_view="tree", pan_data=PanData(self.run), preserve_branch_lengths=True
        )
        self.addCleanup(tab.close)
        lengths = tab._tree_lengths_combo
        self.assertEqual(
            [lengths.itemText(i) for i in range(lengths.count())],
            ["Mutations only", "Including recombination", "Topology"],
        )
        self.assertEqual(lengths.currentData(), "mutations")
        old_state = {
            "layout": "rectangular",
            "branch_scale": "distance",
            "branch_lengths": "with_recombination",
            "true_lengths": "on",
            "order": "ladder_depth",
            "root": "midpoint",
            "source": "core_recomb",
        }
        tab.set_tree_display_options(old_state)
        self.assertEqual(lengths.currentData(), "with_recombination")
        self.assertEqual(tab.branch_length_set_id(), "with_recombination")
        self.assertTrue(tab._tree_true_lengths_check.isChecked())
        self.assertEqual(tab._tree_order_combo.currentData(), "ladder_depth")
        self.assertEqual(tab._tree_root_combo.currentData(), "midpoint")
        # Saved under the same keys (the default source is not saved).
        self.assertEqual(tab.tree_display_options(), {k: v for k, v in old_state.items() if k != "source"})
        # Topology keeps the lengths picked for when Distance comes back.
        tab.set_tree_display_options({**old_state, "branch_scale": "topology"})
        self.assertEqual(lengths.currentData(), LENGTHS_TOPOLOGY)
        self.assertEqual(tab.tree_display_options()["branch_lengths"], "with_recombination")
        # Picking a set in the one dropdown draws it to scale.
        lengths.setCurrentIndex(lengths.findData("mutations"))
        options = tab.tree_display_options()
        self.assertEqual(options["branch_scale"], "distance")
        self.assertNotIn("branch_lengths", options)
        # The fold is saved once the user picks it.
        self.assertNotIn("settings", options)
        tab.set_tree_settings_open(False)
        self.assertEqual(tab.tree_display_options()["settings"], "folded")
        # Not shown yet: opened when shown in a panel with room for it.
        tab.set_tree_display_options({"settings": "open"})
        self.assertEqual(tab.tree_display_options()["settings"], "open")


if __name__ == "__main__":
    unittest.main()
