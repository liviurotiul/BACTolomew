from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

import pandas as pd
from PyQt6 import QtCore, QtGui, QtTest, QtWidgets

from annotation_explorer_polars.context_canvas import (
    ALL_LABELS_MIN_SCALE,
    CANVAS_BACKGROUND,
    KEY_LABELS_MIN_SCALE,
    MATCH_TOOLTIP,
    MIN_TEXT_PX,
    REFERENCE_TOOLTIP,
    ROW_HOVER_COLOR,
    SELECTED_COLOR,
    TEXT_COLOR,
    ContextCanvas,
    _GeneNameItem,
    _primary_row_label,
    _molecule_word,
    _paint_rail_row,
    _rail_row_geometry,
    _rail_row_layout,
    _row_sentence,
    _secondary_row_label,
)
from annotation_explorer_polars.context_model import CATEGORY_COLORS
from annotation_explorer_polars import context_canvas
from annotation_explorer_polars import igv_collection_view
from annotation_explorer_polars.igv_collection_view import (
    _MOBILE_NAME_PATTERN,
    CONTEXT_FOOTER_HINT,
    IgvCollectionView,
    LegacyIgvCollectionView,
    NeighbourhoodView,
)
from annotation_explorer_polars.theme import APP_STYLESHEET, COLORS


def _viewer(parent=None, **kwargs) -> NeighbourhoodView:
    """The neighborhood viewer as a standalone window, at the size the old
    ``IgvCollectionDialog`` shim opened with (the shim went in Phase 10)."""
    view = NeighbourhoodView(parent, **kwargs)
    view.resize(1760, 860)
    return view


def _choose(box: QtWidgets.QComboBox, value) -> None:
    """Pick the item holding ``value``, as a user choosing it from the list."""
    index = box.findData(value)
    assert index >= 0, value
    box.setCurrentIndex(index)


def fixture_rows(row_count: int = 4, gene_count: int = 9):
    rows = []
    for row_index in range(row_count):
        genes = []
        base = 100_000 * (row_index + 1)
        for gene_index in range(gene_count):
            start = base + gene_index * 800
            name = ("dnaA", "gyrB", "blaOXA-23", "tnpA", "recA")[gene_index % 5]
            genes.append({
                "start": start,
                "end": start + 520,
                "strand": "-" if (gene_index + row_index) % 4 == 0 else "+",
                "label": name,
                "is_center": gene_index == gene_count // 2,
                "row_data": {
                    "Gene": name,
                    "Product": name,
                    "PanarooCluster": f"cluster_{gene_index}",
                    "ConsensusCategory": "AMR" if name.startswith("bla") else "",
                    "unknown_gene_field": "preserved",
                },
            })
        rows.append({
            "sample": "REFERENCE" if row_index == 0 else f"sample_{row_index}",
            "contig": f"contig_{row_index}",
            "is_source": row_index == 0,
            "center_bp": genes[gene_count // 2]["start"] + 260,
            "context_similarity": 1.0 if row_index == 0 else 0.8,
            "context_orientation": "reverse" if row_index == 2 else "forward",
            "molecule": "plasmid" if row_index % 2 else "chromosome",
            "genes": genes,
            "sample_meta": {"Country": "Romania"},
        })
    return rows


def short_contig_rows():
    """Reference plus one partial neighbor row and three focal-only contigs."""
    rows = fixture_rows(2, 9)
    rows[0].update(contig_genes_before=4, contig_genes_after=40)
    rows[1].update(contig_genes_before=4, contig_genes_after=4)
    for index, molecule in enumerate(("plasmid", "chromosome", "chromosome")):
        focal = dict(rows[0]["genes"][4], start=1000, end=1520)
        rows.append({
            "sample": f"short_{index}",
            "contig": f"{90 + index}",
            "center_bp": 1260,
            "molecule": molecule,
            "context_similarity": 0.06,
            "contig_genes_before": 0,
            "contig_genes_after": 0,
            "genes": [focal],
        })
    return rows


class ContextCanvasTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def setUp(self) -> None:
        self.canvas = ContextCanvas()
        self.canvas.resize(1100, 500)
        self.canvas.set_rows(fixture_rows())
        self.canvas.show()
        self.app.processEvents()

    def tearDown(self) -> None:
        self.canvas.close()
        self.app.processEvents()

    def test_public_compatibility_methods_and_modes(self) -> None:
        self.canvas.set_display_mode("figure")
        self.canvas.set_density_mode("comfortable")
        self.canvas.set_normalize_orientation(False)
        self.canvas.set_synteny_blocks_enabled(False)
        self.canvas.set_static_row_labels(True, track_offset=220)
        self.canvas.fit_all_contexts()
        self.canvas.reset_zoom()
        self.assertEqual(self.canvas._display_mode, "figure")
        self.assertEqual(self.canvas._density, "comfortable")
        self.assertFalse(self.canvas._normalize_orientation)

    def test_publication_style_changes_only_the_figure_and_keeps_choices(self) -> None:
        self.canvas.set_color_mode("family")
        self.canvas.set_label_mode("all")
        before = self.canvas.graphics_view.viewport().grab().toImage()
        before_bytes = before.bits().asstring(before.sizeInBytes())
        export_dir = tempfile.TemporaryDirectory()
        self.addCleanup(export_dir.cleanup)
        explore_export = Path(export_dir.name) / "explore.png"
        publication_export = Path(export_dir.name) / "publication.png"
        self.assertTrue(self.canvas.export_png(explore_export, scale=1.0))
        self.canvas.set_export_publication_style(True)
        self.assertTrue(self.canvas.export_png(publication_export, scale=1.0))
        self.assertNotEqual(explore_export.read_bytes(), publication_export.read_bytes())
        self.app.processEvents()
        after = self.canvas.graphics_view.viewport().grab().toImage()
        self.assertEqual(before_bytes, after.bits().asstring(after.sizeInBytes()))
        self.assertEqual(self.canvas._display_mode, "exploratory")
        self.assertEqual(self.canvas.color_mode, "family")
        self.assertEqual(self.canvas.label_mode, "all")
        # The compatibility display mode restyles but never resets choices.
        self.canvas.set_display_mode("figure")
        self.assertEqual(self.canvas.color_mode, "family")
        self.assertEqual(self.canvas.label_mode, "all")
        self.canvas.set_display_mode("exploratory")
        self.assertEqual(self.canvas.graphics_view.backgroundBrush().color().name(), COLORS["surface"])

    def test_default_color_is_gene_type_with_table_tag_colors(self) -> None:
        self.assertEqual(ContextCanvas().color_mode, "annotation")
        rows = fixture_rows(2, 7)
        rows[0]["genes"][0]["color"] = "#ffe082"
        self.canvas.set_rows(rows)
        self.app.processEvents()
        genes = self.canvas.layout_result.rows[0].genes
        uncategorized = next(gene for gene in genes if not gene.gene.category)
        amr = next(gene for gene in genes if gene.gene.category == "AMR")
        fill = self.canvas._row_items[0]._fill_for_gene
        self.assertEqual(fill(uncategorized).name(), CATEGORY_COLORS["OTHER"])
        self.assertEqual(fill(amr).name(), "#c62828")
        legend = dict(self.canvas._color_legend_entries())
        self.assertEqual(legend["Other"].name(), CATEGORY_COLORS["OTHER"])
        self.assertIn("AMR", legend)

        all_amr = fixture_rows(2, 7)
        for row in all_amr:
            for value in row["genes"]:
                value["row_data"]["ConsensusCategory"] = "AMR"
        self.canvas.set_rows(all_amr)
        self.app.processEvents()
        self.assertNotIn("Other", dict(self.canvas._color_legend_entries()))
        self.canvas.set_color_mode("family")
        self.assertEqual(fill(uncategorized).name(), uncategorized.gene.family_fill)

    def test_legend_card_is_always_visible_and_folds(self) -> None:
        self.canvas.graphics_view.viewport().grab()
        rect = self.canvas._legend_rect
        self.assertFalse(rect.isEmpty())
        self.assertTrue(QtCore.QRectF(self.canvas.graphics_view.viewport().rect()).contains(rect))
        kinds = [kind for group in self.canvas._legend_groups("2 kb") for kind, _text, _color in group]
        for kind in ("swatch", "link", "focal", "scale"):
            self.assertIn(kind, kinds)
        texts = [text for group in self.canvas._legend_groups("2 kb") for _kind, text, _color in group]
        self.assertIn("same direction", texts)
        self.assertIn("inverted", texts)
        point = QtCore.QPointF(rect.center())
        press = QtGui.QMouseEvent(
            QtCore.QEvent.Type.MouseButtonPress, point, point,
            QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.MouseButton.LeftButton,
            QtCore.Qt.KeyboardModifier.NoModifier,
        )
        self.assertTrue(self.canvas._handle_mouse_press(press))
        self.assertTrue(self.canvas.legend_collapsed)
        self.canvas.graphics_view.viewport().grab()
        self.assertLess(self.canvas._legend_rect.width(), rect.width())
        self.assertFalse(self.canvas.has_selection())
        # In a figure the key sits below the plot instead of on top of genes.
        _width, height, _scale = self.canvas._document_size()
        self.assertGreater(height, 42.0 + self.canvas.layout_result.height + 30.0)

    def test_focal_band_replaces_focal_to_focal_ribbons(self) -> None:
        self.assertIsNotNone(self.canvas._focal_band_rect)
        self.assertEqual(self.canvas.focal_label, "recA")
        for item in self.canvas._pair_items:
            for link, _path in item._paths:
                self.assertFalse(link.upper_gene.gene.is_focal and link.lower_gene.gene.is_focal)
        self.assertTrue(any(item._paths for item in self.canvas._pair_items))
        header = self.canvas.plot_caption
        font = header.font()
        font.setPixelSize(11)
        tab = header._focal_tab_rect(font)
        self.assertIsNotNone(tab)
        self.assertIn("recA", header.tooltip_at(tab.center().x()))

    def test_zoom_does_not_rebuild_scene_items(self) -> None:
        before_rows = tuple(self.canvas._row_items)
        before_pairs = tuple(self.canvas._pair_items)
        before_count = len(self.canvas._scene.items())
        self.canvas.zoom_by(1.4)
        self.app.processEvents()
        self.assertEqual(before_rows, tuple(self.canvas._row_items))
        self.assertEqual(before_pairs, tuple(self.canvas._pair_items))
        self.assertEqual(before_count, len(self.canvas._scene.items()))

    def test_translucent_links_do_not_use_unstable_device_cache(self) -> None:
        self.assertTrue(self.canvas._pair_items)
        self.assertEqual(
            self.canvas._row_items[0].cacheMode(),
            QtWidgets.QGraphicsItem.CacheMode.DeviceCoordinateCache,
        )
        self.assertEqual(
            self.canvas._pair_items[0].cacheMode(),
            QtWidgets.QGraphicsItem.CacheMode.ItemCoordinateCache,
        )

    def test_selected_link_cache_is_invalidated_even_when_offscreen(self) -> None:
        self.canvas.set_rows(fixture_rows(30, 9))
        self.app.processEvents()
        self.canvas.set_link_mode("selected")
        offscreen_pair = self.canvas._pair_items[-1]
        calls = []
        offscreen_pair.update = lambda *_args: calls.append(True)
        self.canvas._select_gene(self.canvas.layout_result.rows[0].genes[0])
        self.assertTrue(calls)

    def test_keyboard_open_emits_original_annotation_payload(self) -> None:
        spy = QtTest.QSignalSpy(self.canvas.geneDoubleClicked)
        target = self.canvas.layout_result.rows[0].genes[0]
        self.canvas._select_gene(target)
        self.assertTrue(self.canvas.open_selected_gene())
        self.assertEqual(len(spy), 1)
        emitted = spy[0][0]
        self.assertEqual(emitted["unknown_gene_field"], "preserved")
        self.assertNotIn("__context_family__", emitted)

    def test_selection_survives_data_refresh_when_family_is_visible(self) -> None:
        target = self.canvas.layout_result.rows[0].genes[0]
        self.canvas._select_gene(target)
        selected_family = self.canvas.selected_family
        self.canvas.set_rows(fixture_rows())
        self.app.processEvents()
        self.assertEqual(self.canvas.selected_family, selected_family)
        self.assertIsNotNone(self.canvas._selected_gene)

    def test_selection_stays_on_same_locus_after_row_reordering(self) -> None:
        target = self.canvas.layout_result.rows[1].genes[0]
        self.canvas._select_gene(target)
        self.canvas.set_rows(list(reversed(fixture_rows())))
        self.app.processEvents()
        selected = self.canvas._selected_gene
        self.assertIsNotNone(selected)
        selected_row = self.canvas.layout_result.rows[selected.row_index].row
        self.assertEqual(selected_row.sample, "sample_1")

    def test_selection_is_rebound_after_density_and_orientation_relayout(self) -> None:
        target = self.canvas.layout_result.rows[2].genes[2]
        self.canvas._select_gene(target)
        self.canvas.set_density_mode("comfortable")
        self.canvas.set_normalize_orientation(False)
        self.assertTrue(self.canvas.navigate_selection(1, 0))
        current_genes = {
            id(gene)
            for row in self.canvas.layout_result.rows
            for gene in row.genes
        }
        self.assertIn(id(self.canvas._selected_gene), current_genes)

    def test_complete_png_and_pdf_exports_include_document(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            png = Path(tmp) / "context.png"
            pdf = Path(tmp) / "context.pdf"
            self.assertTrue(self.canvas.export_png(png, scale=1.0))
            self.assertTrue(self.canvas.export_pdf(pdf))
            self.assertGreater(png.stat().st_size, 1000)
            self.assertGreater(pdf.stat().st_size, 1000)

    def test_export_excludes_interactive_selection_chrome(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            before = Path(tmp) / "before.png"
            after = Path(tmp) / "after.png"
            self.assertTrue(self.canvas.export_png(before, scale=1.0))
            self.canvas._select_gene(self.canvas.layout_result.rows[0].genes[0])
            self.assertTrue(self.canvas.export_png(after, scale=1.0))
            self.assertEqual(before.read_bytes(), after.read_bytes())

    def test_empty_export_returns_false(self) -> None:
        self.canvas.clear()
        self.app.processEvents()
        with tempfile.TemporaryDirectory() as tmp:
            self.assertFalse(self.canvas.export_png(Path(tmp) / "empty.png"))

    def test_visible_empty_canvas_paints_without_exception(self) -> None:
        empty = ContextCanvas()
        empty.resize(700, 360)
        empty.set_rows(fixture_rows())
        empty.show()
        self.app.processEvents()
        empty.zoom_by(4.0)
        empty.graphics_view.horizontalScrollBar().setValue(
            empty.graphics_view.horizontalScrollBar().maximum()
        )
        empty.set_rows([])
        self.app.processEvents()
        self.assertTrue(empty.isVisible())
        self.assertAlmostEqual(empty.graphics_view.transform().m11(), 1.0)
        mapped_center = empty.graphics_view.mapFromScene(empty._scene.sceneRect().center())
        self.assertTrue(empty.graphics_view.viewport().rect().contains(mapped_center))
        message = next(
            item
            for item in empty._scene.items()
            if hasattr(item, "text") and item.text().startswith("No neighborhoods")
        )
        message_rect = empty.graphics_view.mapFromScene(message.sceneBoundingRect()).boundingRect()
        self.assertTrue(empty.graphics_view.viewport().rect().intersects(message_rect))
        empty.close()

    def test_hit_tolerance_is_device_stable_across_zoom_levels(self) -> None:
        target = self.canvas.layout_result.rows[0].genes[0]
        point = target.center_x, target.center_y
        self.canvas._set_horizontal_scale(0.05)
        self.assertIsNotNone(
            self.canvas._gene_at_scene(
                self.canvas.graphics_view.mapToScene(
                    self.canvas.graphics_view.mapFromScene(*point)
                )
            )
        )
        self.canvas._set_horizontal_scale(12.0)
        far_outside = type(self.canvas._scene.sceneRect().center())(
            target.right + 2.0,
            target.center_y,
        )
        self.assertIsNone(self.canvas._gene_at_scene(far_outside))

    def _hover_at(self, gene) -> None:
        point = QtCore.QPointF(self.canvas.graphics_view.mapFromScene(gene.center_x, gene.center_y))
        event = QtGui.QMouseEvent(
            QtCore.QEvent.Type.MouseMove, point, point,
            QtCore.Qt.MouseButton.NoButton, QtCore.Qt.MouseButton.NoButton,
            QtCore.Qt.KeyboardModifier.NoModifier,
        )
        self.canvas._handle_mouse_move(event)

    def test_gene_hover_spotlights_family_and_tints_row(self) -> None:
        row = self.canvas.layout_result.rows[1]
        gene = next(g for g in row.genes if not g.gene.is_focal and g.gene.family_key)
        self._hover_at(gene)
        self.assertIs(self.canvas._hovered_gene, gene)
        self.assertEqual(self.canvas.hovered_family, gene.gene.family_key)
        self.assertTrue(self.canvas._veil_item.isVisible())
        self.assertTrue(self.canvas._row_hover_item.isVisible())
        self.assertEqual(self.canvas._row_hover_item.brush().color().name(), ROW_HOVER_COLOR.name())
        self.assertEqual(ROW_HOVER_COLOR.name(), COLORS["row_hover"])
        self.assertGreater(len(self.canvas._dynamic_highlight_items), 1)
        tooltip = self.canvas.graphics_view.toolTip()
        self.assertIn("Found in", tooltip)
        self.assertIn("Double-click for the full gene card", tooltip)

    def test_focal_hover_only_tints_the_row(self) -> None:
        row = self.canvas.layout_result.rows[1]
        focal = next(g for g in row.genes if g.gene.is_focal)
        self._hover_at(focal)
        self.assertIs(self.canvas._hovered_gene, focal)
        self.assertIsNone(self.canvas.hovered_family)
        self.assertFalse(self.canvas._veil_item.isVisible())
        self.assertTrue(self.canvas._row_hover_item.isVisible())
        self.assertEqual(self.canvas._dynamic_highlight_items, [])

    def test_spotlight_follows_the_pointer_through_scroll_and_zoom(self) -> None:
        row = self.canvas.layout_result.rows[0]
        gene = next(g for g in row.genes if not g.gene.is_focal)
        self._hover_at(gene)
        self.canvas.zoom_by(1.0)
        self.assertIs(self.canvas._hovered_gene, gene)
        self.assertTrue(self.canvas._veil_item.isVisible())
        self.canvas.zoom_by(3.0)
        under = self.canvas._pointer_target()[0]
        self.assertIs(self.canvas._hovered_gene, under)
        self.assertEqual(self.canvas._veil_item.isVisible(), under is not None and not under.gene.is_focal)
        # Scrolling one gene along puts the next gene under the still pointer.
        self.canvas.zoom_by(2.0)
        self._hover_at(gene)
        self.assertIs(self.canvas._hovered_gene, gene)
        neighbor = row.genes[row.genes.index(gene) + 1]
        bar = self.canvas.graphics_view.horizontalScrollBar()
        scale_x = self.canvas.graphics_view.transform().m11()
        bar.setValue(bar.value() + round((neighbor.center_x - gene.center_x) * scale_x))
        self.assertIs(self.canvas._hovered_gene, neighbor)
        self.assertIsNot(neighbor, gene)
        self.canvas.graphics_view.leaveEvent(QtCore.QEvent(QtCore.QEvent.Type.Leave))
        self.assertIsNone(self.canvas._hovered_gene)
        self.assertFalse(self.canvas._veil_item.isVisible())

    def test_key_card_avoids_genes_or_folds(self) -> None:
        for size in ((1100, 500), (520, 260)):
            self.canvas.resize(*size)
            self.app.processEvents()
            self.canvas.graphics_view.viewport().grab()
            rect = self.canvas._legend_rect
            self.assertFalse(rect.isEmpty())
            if self.canvas._legend_auto_folded:
                continue
            for gene_rect in self.canvas._visible_gene_rects():
                self.assertFalse(rect.intersects(gene_rect), size)
        # A key the viewer opens on purpose stays open, even over genes.
        self.canvas.set_legend_collapsed(False)
        self.canvas.graphics_view.viewport().grab()
        self.assertFalse(self.canvas._legend_auto_folded)

    def test_key_card_keeps_off_the_fold_line(self) -> None:
        # Long rows fill the right corners, so the key looks to the left,
        # where the no-neighbors sentence starts.
        rows = fixture_rows(3, 31)
        for index in range(3):
            focal = dict(rows[0]["genes"][15], start=1000, end=1520)
            rows.append({
                "sample": f"short_{index}", "contig": f"{90 + index}", "center_bp": 1260,
                "molecule": "chromosome", "context_similarity": 0.06, "genes": [focal],
            })
        self.canvas.set_rows(rows)
        checked = 0
        for size in ((600, 240), (600, 320), (700, 200)):
            self.canvas.resize(*size)
            self.app.processEvents()
            self.canvas.reset_zoom()
            self.app.processEvents()
            self.canvas.graphics_view.viewport().grab()
            link = self.canvas._fold_link_rect
            if link.isEmpty():
                continue
            checked += 1
            if self.canvas._legend_auto_folded:
                continue
            sentence = QtCore.QRectF(12.0, link.top(), link.right() - 12.0, link.height())
            self.assertFalse(self.canvas._legend_rect.intersects(sentence), size)
        self.assertGreater(checked, 0)

    def test_view_state_round_trip(self) -> None:
        self.canvas.set_rows(fixture_rows(40, 31))
        self.canvas.resize(900, 400)
        self.app.processEvents()
        view = self.canvas.graphics_view
        self.canvas.zoom_by(1.6)
        view.verticalScrollBar().setValue(view.verticalScrollBar().maximum() // 2)
        self.app.processEvents()
        state = self.canvas.view_state()
        self.assertGreater(state[2], 0)
        self.canvas.reset_zoom()
        view.verticalScrollBar().setValue(0)
        self.app.processEvents()
        self.assertNotAlmostEqual(view.transform().m11(), state[0])
        self.canvas.restore_view_state(state)
        self.app.processEvents()
        self.assertAlmostEqual(view.transform().m11(), state[0])
        self.assertEqual(view.verticalScrollBar().value(), state[2])

    def test_gene_under_the_key_card_wins_hover_and_click(self) -> None:
        gene = next(g for g in self.canvas.layout_result.rows[1].genes if not g.gene.is_focal)
        center = QtCore.QPointF(self.canvas.graphics_view.mapFromScene(gene.center_x, gene.center_y))
        # Pretend the card is drawn right over the gene.
        self.canvas._legend_rect = QtCore.QRectF(center.x() - 40.0, center.y() - 30.0, 80.0, 60.0)
        self._hover_at(gene)
        self.assertIs(self.canvas._hovered_gene, gene)
        press = QtGui.QMouseEvent(
            QtCore.QEvent.Type.MouseButtonPress, center, center,
            QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.MouseButton.LeftButton,
            QtCore.Qt.KeyboardModifier.NoModifier,
        )
        collapsed = self.canvas.legend_collapsed
        self.assertTrue(self.canvas._handle_mouse_press(press))
        self.assertIs(self.canvas._selected_gene, gene)
        self.assertEqual(self.canvas.legend_collapsed, collapsed)

    def test_every_focal_gene_takes_the_reference_focal_color(self) -> None:
        rows = fixture_rows()
        for index, row in enumerate(rows):
            row["genes"][4]["row_data"]["ConsensusCategory"] = "" if index == 2 else "AMR"
        self.canvas.set_rows(rows)
        self.app.processEvents()
        untyped = next(g for g in self.canvas.layout_result.rows[2].genes if g.gene.is_focal)
        self.assertEqual(untyped.gene.category, "")
        self.assertEqual(self.canvas.gene_fill(untyped).name(), CATEGORY_COLORS["AMR"])
        self.canvas.set_color_mode("family")
        reference_focal = next(g for g in self.canvas.layout_result.rows[0].genes if g.gene.is_focal)
        self.assertEqual(self.canvas.gene_fill(untyped).name(), self.canvas.gene_fill(reference_focal).name())

    def test_spotlit_family_keeps_its_names(self) -> None:
        self.canvas.zoom_by(4.0)
        gene = next(g for g in self.canvas.layout_result.rows[0].genes if g.gene.label == "tnpA")
        self._hover_at(gene)
        names = [item.text for item in self.canvas._dynamic_highlight_items if isinstance(item, _GeneNameItem)]
        self.assertIn("tnpA", names)
        self.canvas.set_label_mode("off")
        self.canvas._rebuild_interaction_highlights()
        names = [item for item in self.canvas._dynamic_highlight_items if isinstance(item, _GeneNameItem)]
        self.assertEqual(names, [])

    def test_presence_text_says_shown_when_rows_are_filtered(self) -> None:
        gene = next(g for g in self.canvas.layout_result.rows[0].genes if not g.gene.is_focal)
        self.assertTrue(self.canvas.presence_text(gene).endswith("genomes with neighbors"))
        self.canvas.set_rows_filtered(True)
        self.assertTrue(self.canvas.presence_text(gene).endswith("genomes with neighbors shown"))

    def test_figure_key_is_drawn_below_the_plot(self) -> None:
        origins = []
        original = self.canvas._paint_legend

        def recording(painter, origin, **kwargs):
            if painter is not None:
                origins.append(QtCore.QPointF(origin))
            return original(painter, origin, **kwargs)

        self.canvas._paint_legend = recording
        export_dir = tempfile.TemporaryDirectory()
        self.addCleanup(export_dir.cleanup)
        self.assertTrue(self.canvas.export_png(Path(export_dir.name) / "key.png", scale=1.0))
        self.assertEqual(len(origins), 1)
        self.assertGreaterEqual(origins[0].y(), 42.0 + self.canvas.layout_result.height)

    def test_pin_is_solid_primary_and_differs_from_hover(self) -> None:
        gene = next(g for g in self.canvas.layout_result.rows[0].genes if not g.gene.is_focal)
        self._hover_at(gene)
        hover_pens = {item.pen().color().name() for item in self.canvas._dynamic_highlight_items if hasattr(item, "pen")}
        self.canvas._select_gene(gene)
        pin_pens = {item.pen().color().name() for item in self.canvas._dynamic_highlight_items if hasattr(item, "pen")}
        self.assertIn(SELECTED_COLOR.name(), pin_pens)
        self.assertNotIn(SELECTED_COLOR.name(), hover_pens)
        self.assertTrue(self.canvas._veil_item.isVisible())
        self.canvas.clear_selection()
        self.canvas.graphics_view.leaveEvent(QtCore.QEvent(QtCore.QEvent.Type.Leave))
        self.assertFalse(self.canvas._veil_item.isVisible())

    def test_key_genes_label_typed_genes_on_every_row(self) -> None:
        row_item = next(item for item in self.canvas._row_items if not item.row.row.is_reference)
        amr = next(g for g in row_item.row.genes if g.gene.category == "AMR" and not g.gene.is_focal)
        plain = next(g for g in row_item.row.genes if not g.gene.category and not g.gene.is_focal)
        self.assertTrue(row_item._should_label(amr, KEY_LABELS_MIN_SCALE))
        self.assertFalse(row_item._should_label(plain, KEY_LABELS_MIN_SCALE))
        # The reference row is named at every zoom, typed or not.
        reference_item = next(item for item in self.canvas._row_items if item.row.row.is_reference)
        for gene in reference_item.row.genes:
            if not gene.gene.is_focal:
                self.assertTrue(reference_item._should_label(gene, 0.1), gene.gene.label)
        self.canvas.set_label_mode("off")
        self.assertFalse(row_item._should_label(amr, 4.0))

    def test_viewport_motion_clears_stale_hover_row(self) -> None:
        row = self.canvas.layout_result.rows[1]
        self.canvas._set_hover(None, row)
        self.assertEqual(self.canvas.hovered_row_index, 1)
        self.canvas.graphics_view.verticalScrollBar().setValue(
            self.canvas.graphics_view.verticalScrollBar().value() + 1
        )
        self.canvas._viewport_changed()
        self.assertIsNone(self.canvas.hovered_row_index)

    def test_generic_gene_has_exact_hover_and_selection_state(self) -> None:
        rows = [{
            "sample": "generic",
            "genes": [{
                "start": 10,
                "end": 11,
                "strand": "+",
                "label": "hypothetical protein",
                "row_data": {"Product": "hypothetical protein"},
            }],
        }]
        self.canvas.set_rows(rows)
        self.app.processEvents()
        target = self.canvas.layout_result.rows[0].genes[0]
        self.canvas._set_hover(target, self.canvas.layout_result.rows[0])
        self.assertIs(self.canvas._hovered_gene, target)
        self.canvas._select_gene(target)
        self.assertIs(self.canvas._selected_gene, target)
        self.assertIsNone(self.canvas.selected_family)

    def test_all_label_mode_never_shows_fewer_labels_than_default(self) -> None:
        for scale_x in (0.3, 0.55, 0.7, 1.0):
            for item in self.canvas._row_items:
                self.canvas.set_label_mode("reference")
                default = {id(g) for g in item.row.genes if item._should_label(g, scale_x)}
                self.canvas.set_label_mode("all")
                everything = {id(g) for g in item.row.genes if item._should_label(g, scale_x)}
                self.assertLessEqual(default, everything, f"scale {scale_x}")

    def test_horizontal_scrollbar_only_when_needed(self) -> None:
        self.assertEqual(
            self.canvas.graphics_view.horizontalScrollBarPolicy(),
            QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded,
        )
        bar = self.canvas.graphics_view.horizontalScrollBar()
        self.canvas.fit_all_contexts()
        self.app.processEvents()
        self.assertFalse(bar.isVisible())
        self.canvas.zoom_by(4.0)
        self.app.processEvents()
        self.assertTrue(bar.isVisible())

    def test_rows_are_top_aligned(self) -> None:
        alignment = self.canvas.graphics_view.alignment()
        self.assertTrue(alignment & QtCore.Qt.AlignmentFlag.AlignTop)

    def test_theme_tokens_and_minimum_text(self) -> None:
        self.assertGreaterEqual(MIN_TEXT_PX, 10)
        self.assertEqual(TEXT_COLOR.name(), COLORS["text"])
        self.assertEqual(SELECTED_COLOR.name(), COLORS["primary"])
        self.assertEqual(CANVAS_BACKGROUND.name(), COLORS["surface"])

    def test_all_labels_adds_names_at_opening_zoom(self) -> None:
        row_item = next(item for item in self.canvas._row_items if not item.row.row.is_reference)
        gene = next(g for g in row_item.row.genes if g.gene.label and not g.gene.is_focal)
        scale = max(ALL_LABELS_MIN_SCALE, 0.3)
        self.assertFalse(row_item._should_label(gene, scale))
        self.canvas.set_label_mode("all")
        self.assertTrue(row_item._should_label(gene, scale))

    def test_selected_gene_stays_visible_when_plot_narrows(self) -> None:
        view = self.canvas.graphics_view
        self.canvas.zoom_by(4.0)
        self.app.processEvents()
        row = self.canvas.layout_result.rows[0]
        gene = row.genes[0]
        view.centerOn(gene.center_x, gene.center_y)
        self.canvas._select_gene(gene)
        self.canvas.resize(500, 500)
        self.app.processEvents()
        visible = view.mapToScene(view.viewport().rect()).boundingRect()
        self.assertTrue(visible.contains(QtCore.QPointF(gene.center_x, gene.center_y)))

    def test_merged_row_second_line_does_not_name_one_member(self) -> None:
        rows = fixture_rows()
        rows[0].update(collapsed_sample_count=3, collapsed_track_count=3, display_label="3 genomes")
        rows[1].update(collapsed_sample_count=2, collapsed_track_count=2, display_label="2 genomes")
        self.canvas.set_rows(rows)
        self.app.processEvents()
        reference, merged = self.canvas.layout_result.rows[:2]
        self.assertEqual(_secondary_row_label(reference), "includes REFERENCE")
        self.assertEqual(_secondary_row_label(merged), "2 contigs")

    def test_gene_tooltip_names_minus_strand(self) -> None:
        row = self.canvas.layout_result.rows[0]
        minus = next(gene for gene in row.genes if gene.gene.strand == "-")
        tooltip = self.canvas._gene_tooltip(minus, row)
        self.assertIn("minus strand", tooltip)
        self.assertIn("contig_0, ", tooltip)
        self.assertIn("Found in", tooltip)
        self.assertNotIn("Family evidence", tooltip)
        self.assertNotIn("\u2013", tooltip)

    def test_rail_tooltips_explain_match_reference_and_badges(self) -> None:
        rail = self.canvas.row_header
        rows = self.canvas.layout_result.rows
        view = self.canvas.graphics_view

        def rail_y(row) -> float:
            return float(view.mapFromScene(QtCore.QPointF(0.0, row.center_y)).y())

        match_x = rail.width() - 20.0
        self.assertEqual(rail.tooltip_at(match_x, rail_y(rows[0])), REFERENCE_TOOLTIP)
        self.assertEqual(rail.tooltip_at(match_x, rail_y(rows[1])), MATCH_TOOLTIP)
        plasmid = next(row for row in rows if row.row.molecule == "plasmid")
        parts = _rail_row_geometry(plasmid, rail_y(plasmid), float(rail.width()), True)
        for key in ("dot", "molecule"):
            self.assertIn("Plasmid contig", rail.tooltip_at(parts[key].center().x(), parts[key].center().y()))
        # The name explains the row in one sentence, then its metadata.
        sentence = rail.tooltip_at(30.0, rail_y(plasmid) - 10.0)
        self.assertTrue(sentence.startswith("Plasmid contig_1 of sample_1. 80% match to the reference."))
        self.assertNotIn("Sample:", sentence)
        self.assertEqual(self.canvas.rail_caption.tooltip_at(self.canvas.rail_caption.width() - 5.0), MATCH_TOOLTIP)

    def test_ribbon_click_selects_family_instead_of_link_panel(self) -> None:
        spy = QtTest.QSignalSpy(self.canvas.syntenyBlockClicked)
        item = next(item for item in self.canvas._pair_items if item._paths)
        link, _path = item._paths[0]
        upper = link.upper_gene
        point = QtCore.QPointF(upper.center_x, (upper.center_y + link.lower_gene.center_y) / 2.0)
        self.assertIsNotNone(self.canvas._link_at(point))
        self.assertIsNone(self.canvas._gene_at_scene(point))
        view_point = QtCore.QPointF(self.canvas.graphics_view.mapFromScene(point))
        event = QtGui.QMouseEvent(
            QtCore.QEvent.Type.MouseButtonPress,
            view_point,
            view_point,
            QtCore.Qt.MouseButton.LeftButton,
            QtCore.Qt.MouseButton.LeftButton,
            QtCore.Qt.KeyboardModifier.NoModifier,
        )
        self.assertTrue(self.canvas._handle_mouse_press(event))
        self.assertEqual(len(spy), 0)
        self.assertEqual(self.canvas.selected_family, link.family_key)
        self.assertTrue(self.canvas._handle_mouse_press(event))
        self.assertFalse(self.canvas.has_selection())

    def test_ribbon_without_family_still_toggles_off(self) -> None:
        from types import SimpleNamespace

        gene = self.canvas.layout_result.rows[0].genes[0]
        link = SimpleNamespace(family_key="", upper_gene=gene)
        self.canvas._gene_at_scene = lambda _point: None
        self.canvas._link_at = lambda _point: link
        event = QtGui.QMouseEvent(
            QtCore.QEvent.Type.MouseButtonPress,
            QtCore.QPointF(5.0, 5.0),
            QtCore.QPointF(5.0, 5.0),
            QtCore.Qt.MouseButton.LeftButton,
            QtCore.Qt.MouseButton.LeftButton,
            QtCore.Qt.KeyboardModifier.NoModifier,
        )
        self.canvas._handle_mouse_press(event)
        self.assertTrue(self.canvas.has_selection())
        self.canvas._handle_mouse_press(event)
        self.assertFalse(self.canvas.has_selection())

    def test_escape_clears_selection_then_passes_through(self) -> None:
        self.canvas._select_gene(self.canvas.layout_result.rows[0].genes[0])
        first = QtGui.QKeyEvent(QtCore.QEvent.Type.KeyPress, QtCore.Qt.Key.Key_Escape, QtCore.Qt.KeyboardModifier.NoModifier)
        self.canvas.graphics_view.keyPressEvent(first)
        self.assertTrue(first.isAccepted())
        self.assertFalse(self.canvas.has_selection())
        second = QtGui.QKeyEvent(QtCore.QEvent.Type.KeyPress, QtCore.Qt.Key.Key_Escape, QtCore.Qt.KeyboardModifier.NoModifier)
        self.canvas.graphics_view.keyPressEvent(second)
        self.assertFalse(second.isAccepted())

    def test_default_alias_points_to_v2(self) -> None:
        self.assertIs(IgvCollectionView, ContextCanvas)

    def test_collapsed_row_uses_group_label_not_representative_sample(self) -> None:
        rows = fixture_rows(1, 5)
        rows[0]["sample"] = "arbitrary_representative"
        rows[0]["display_label"] = "12 samples"
        rows[0]["collapsed_samples"] = [f"sample_{index}" for index in range(12)]
        rows[0]["collapsed_sample_count"] = 12
        self.canvas.set_rows(rows)
        self.app.processEvents()
        row = self.canvas.layout_result.rows[0].row
        self.assertEqual(_primary_row_label(row), "12 samples")

    def test_short_contigs_fold_into_one_summary_line(self) -> None:
        rows = short_contig_rows()
        self.canvas.set_rows(rows)
        self.app.processEvents()
        layout = self.canvas.layout_result
        self.assertEqual([row.row.sample for row in layout.rows], ["REFERENCE", "sample_1"])
        self.assertIsNotNone(layout.fold)
        self.assertFalse(layout.fold.expanded)
        self.assertEqual(
            self.canvas.fold_sentence(),
            "3 contigs carry recA with no neighbors (2 chromosome, 1 plasmid)",
        )
        # No ribbons run to folded contigs, and the fold sits below the rows.
        self.assertGreater(layout.fold.center_y, layout.rows[-1].center_y)
        self.assertEqual(len(layout.pairs), 1)
        # Presence counts still cover every genome, folded ones included.
        focal = next(gene for gene in layout.rows[0].genes if gene.gene.is_focal)
        self.assertIn("of 5 genomes", self.canvas.presence_text(focal))
        # Clicking the line in the plot shows the contigs as thin rows.
        view = self.canvas.graphics_view
        point = QtCore.QPointF(view.mapFromScene(QtCore.QPointF(0.0, layout.fold.center_y)))
        event = QtGui.QMouseEvent(
            QtCore.QEvent.Type.MouseButtonPress, point, point,
            QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.MouseButton.LeftButton,
            QtCore.Qt.KeyboardModifier.NoModifier,
        )
        self.assertTrue(self.canvas._handle_mouse_press(event))
        layout = self.canvas.layout_result
        self.assertTrue(layout.fold.expanded)
        thin = [row for row in layout.rows if row.short]
        self.assertEqual(len(thin), 3)
        self.assertTrue(all(row.pitch < layout.rows[0].pitch for row in thin))
        self.assertTrue(all(row.center_y > layout.fold.center_y for row in thin))
        rail = self.canvas.row_header
        y = float(view.mapFromScene(QtCore.QPointF(0.0, thin[0].center_y)).y())
        self.assertEqual(rail.tooltip_at(rail.width() - 10.0, y), "Contig too short to compare")
        self.assertIn("too short to compare", _row_sentence(thin[0], "blaOXA-23"))
        self.canvas.toggle_short_rows()
        self.assertFalse(self.canvas.layout_result.fold.expanded)

    def test_reference_row_is_never_folded(self) -> None:
        rows = short_contig_rows()
        rows[0]["genes"] = [gene for gene in rows[0]["genes"] if gene["is_center"]]
        self.canvas.set_rows(rows)
        self.app.processEvents()
        self.assertEqual(self.canvas.layout_result.rows[0].row.sample, "REFERENCE")

    def test_partial_rows_and_contig_end_caps(self) -> None:
        dialog = _viewer(
            None, rows=short_contig_rows(), title="Gene neighborhoods: blaOXA-23",
            flank_sizes=(5, 10, 15, 30), initial_flank=5,
        )
        rows = dialog.view.layout_result.rows
        reference, partial = rows[0], rows[1]
        # 4 genes before the focal gene and a flank of 5: the contig ends.
        self.assertTrue(reference.row.partial)
        self.assertTrue(reference.cap_left)
        self.assertFalse(reference.cap_right)
        self.assertTrue(partial.cap_left and partial.cap_right)
        self.assertIn("ends within the shown flank", _row_sentence(partial))
        # At a flank of 30 the reference's long side is cut by the flank, not
        # by the contig, so that side has no cap.
        dialog._set_context_flank(30)
        reference = dialog.view.layout_result.rows[0]
        self.assertTrue(reference.cap_left)
        self.assertFalse(reference.cap_right)
        self.assertIn(("contig_end", "contig ends here", None), [
            entry for group in dialog.view._legend_groups("1 kb") for entry in group
        ])
        dialog.close()

    def test_mirroring_swaps_contig_end_caps(self) -> None:
        rows = short_contig_rows()[:2]
        rows[1].update(context_orientation="reverse", contig_ends_low=True, contig_ends_high=False)
        self.canvas.set_rows(rows)
        row = self.canvas.layout_result.rows[1]
        self.assertTrue(row.flipped)
        self.assertFalse(row.cap_left)
        self.assertTrue(row.cap_right)

    def test_figures_respect_include_short_contigs(self) -> None:
        self.canvas.set_rows(short_contig_rows())
        self.app.processEvents()
        with tempfile.TemporaryDirectory() as tmp:
            self.canvas.set_export_include_short_contigs(True)
            with self.canvas._figure_layout():
                self.assertEqual(sum(row.short for row in self.canvas.layout_result.rows), 3)
            self.canvas.set_export_include_short_contigs(False)
            with self.canvas._figure_layout():
                self.assertEqual(len(self.canvas.layout_result.rows), 2)
                self.assertIn("(not shown)", self.canvas._omitted_short_note())
            self.assertTrue(self.canvas.export_png(Path(tmp) / "figure.png", scale=1.0))
        # The screen keeps its folded line afterwards.
        self.assertEqual(len(self.canvas.layout_result.rows), 2)
        self.assertFalse(self.canvas.layout_result.fold.expanded)

    def test_rail_uses_contig_words_and_quiet_hints(self) -> None:
        rows = fixture_rows()
        rows[1]["contig"] = "48"
        rows[3]["molecule"] = "unknown"
        self.canvas.set_rows(rows)
        layout = self.canvas.layout_result
        self.assertEqual(_secondary_row_label(layout.rows[1]), "contig 48")
        mirrored = layout.rows[2]
        self.assertTrue(mirrored.flipped)
        self.assertNotIn("reversed", _secondary_row_label(mirrored))
        parts, words = _rail_row_layout(mirrored, 100.0, 280.0, True)
        self.assertIn("mirrored", parts)
        self.assertIn(" \u00b7 mirrored", [text for _key, text, _rect, _color in words])
        self.assertNotIn("mirrored", _rail_row_geometry(layout.rows[1], 100.0, 280.0, True))
        # Line two of a contig row reads 'Plasmid \u00b7 contig 48'.
        _parts, words = _rail_row_layout(layout.rows[1], 100.0, 280.0, True)
        self.assertEqual(
            "".join(text for _key, text, _rect, _color in words),
            f"{_molecule_word(layout.rows[1].row.molecule)} \u00b7 contig 48",
        )
        self.assertTrue(_row_sentence(mirrored).endswith("Shown mirrored."))
        self.assertIn("molecule type unknown", _row_sentence(layout.rows[3]))
        self.assertEqual(_molecule_word("unknown"), "Unknown")
        self.assertEqual(_molecule_word("chromosome"), "Chromosome")
        self.assertEqual(_molecule_word("plasmid"), "Plasmid")
        self.assertIn("Reference:", _row_sentence(layout.rows[0]))

    def test_rail_paints_no_chips_or_bars(self) -> None:
        for name in ("_paint_chip", "_paint_match_bar", "_chip_width", "MOLECULE_CHIPS", "METRIC_CHUNK"):
            self.assertFalse(hasattr(context_canvas, name), name)
        rows = fixture_rows()
        rows[1]["context_partial"] = True
        self.canvas.set_rows(rows)
        image = QtGui.QImage(280, 60 * 4, QtGui.QImage.Format.Format_ARGB32)
        image.fill(QtGui.QColor(COLORS["surface"]))
        painter = QtGui.QPainter(image)
        for index, row in enumerate(self.canvas.layout_result.rows[:4]):
            _paint_rail_row(painter, row, 30.0 + 60.0 * index, 280.0)
        painter.end()
        colors = {
            QtGui.QColor(image.pixel(x, y)).name()
            for x in range(image.width())
            for y in range(image.height())
        }
        # No pill fills: the old chip and match-bar tints never appear.
        for tint in ("#e3edff", "#e7f5ec", "#eef2f6", "#86acee", "#edf2f7"):
            self.assertNotIn(tint, colors)
        self.assertIn(COLORS["primary"], colors)
        # The match is plain text in the right column.
        _parts, words = _rail_row_layout(self.canvas.layout_result.rows[1], 30.0, 280.0, True)
        self.assertEqual(words[-1][1], " \u00b7 partial")

    def test_row_right_click_asks_for_a_menu(self) -> None:
        spy = QtTest.QSignalSpy(self.canvas.rowMenuRequested)
        row = self.canvas.layout_result.rows[1]
        self.canvas._request_row_menu(row, QtCore.QPoint(3, 4))
        self.assertEqual(len(spy), 1)
        self.assertEqual(spy[0][0]["sample"], "sample_1")

    def test_plasmid_tooltips_use_contig_coordinates_not_offsets(self) -> None:
        rows = fixture_rows(1, 3)
        rows[0].update(is_plasmid=True, cluster_key="primary:1", joined_contigs=True,
                       title="REFERENCE, plasmid pX1 (IncFII, MOBF)")
        rows[0]["genes"][0]["row_data"].update(Contig="7", StartInt=10, StopInt=530)
        self.canvas.set_rows(rows)
        row = self.canvas.layout_result.rows[0]
        first = next(gene for gene in row.genes if gene.gene.payload.get("Contig") == "7")
        other = next(gene for gene in row.genes if gene is not first)
        self.assertEqual(_secondary_row_label(row), "plasmid pX1 (IncFII, MOBF)")
        self.assertIn("plasmid pX1 (IncFII, MOBF), contig 7, 10 to 530", self.canvas._gene_tooltip(first, row))
        self.assertIn("position in plasmid", self.canvas._gene_tooltip(other, row))


class ContextDialogTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def setUp(self) -> None:
        # The panel's open state is remembered for the session; start each
        # test from the default (no choice made yet).
        NeighbourhoodView.settings_panel_open = None

    def test_gene_neighborhood_constructor_contract(self) -> None:
        dialog = _viewer(
            None,
            rows=fixture_rows(),
            title="Gene neighborhoods: blaOXA-23",
            flank_sizes=(5, 10, 15, 30),
            initial_flank=15,
        )
        dialog.show()
        self.app.processEvents()
        self.assertIsInstance(dialog.view, ContextCanvas)
        self.assertEqual(dialog.flank_box.currentData(), 15)
        self.assertEqual(dialog.flank_box.currentText(), "±15")
        self.assertEqual(dialog.settings_title.text(), "Settings")
        self.assertTrue(dialog.settings_open())
        self.assertTrue(dialog.filter_rail.isVisible())
        dialog.settings_collapse_btn.click()
        self.app.processEvents()
        self.assertFalse(dialog.filter_rail.isVisible())
        dialog._focus_context_search()
        self.app.processEvents()
        self.assertTrue(dialog.filter_rail.isVisible())
        dialog.sample_query.setText("definitely-no-such-sample")
        QtTest.QTest.qWait(190)
        self.app.processEvents()
        self.assertEqual(len(dialog.view.layout_result.rows), 0)
        self.assertEqual(dialog.settings_title.text(), "Settings")
        self.assertTrue(dialog.settings_filter_badge.isVisible())
        self.assertEqual(dialog.settings_filter_badge.text(), "1 filter")
        dialog._clear_context_filters()
        self.assertGreater(len(dialog.view.layout_result.rows), 0)
        dialog.close()

    def test_no_display_menu_and_publication_style_only_affects_export(self) -> None:
        dialog = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x")
        dialog.show()
        self.app.processEvents()
        for name in ("display_btn", "_context_display_actions", "context_flank_box",
                     "related_analysis_btn", "_context_design"):
            self.assertFalse(hasattr(dialog, name), name)
        _choose(dialog.labels_box, "all")
        dialog.publication_style_action.setChecked(True)
        self.assertTrue(dialog.view.export_publication_style)
        self.assertEqual(dialog.view._display_mode, "exploratory")
        self.assertEqual(dialog.view.label_mode, "all")
        self.assertEqual(dialog.view.color_mode, "annotation")
        menu_texts = [action.text() for action in dialog.export_btn.menu().actions()]
        for text in ("Figure as PNG (2x)", "Figure as PDF (vector)", "Publication style",
                     "Compact rows", "Comfortable rows", "Include short contigs"):
            self.assertIn(text, menu_texts)
        self.assertTrue(dialog.include_short_contigs_action.isChecked())
        dialog.include_short_contigs_action.setChecked(False)
        self.assertFalse(dialog._export_include_short_contigs)
        dialog.close()

    def test_zoom_buttons_do_what_they_say_and_more_menu_is_gone(self) -> None:
        dialog = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x")
        dialog.resize(1400, 700)
        dialog.show()
        self.app.processEvents()
        view = dialog.view
        view.fit_all_contexts()
        dialog.fit_reference_btn.click()
        fitted_reference = view.graphics_view.transform().m11()
        view.reset_zoom()
        self.assertAlmostEqual(fitted_reference, view.graphics_view.transform().m11())
        dialog.show_all_btn.click()
        shown_all = view.graphics_view.transform().m11()
        view.fit_all_contexts()
        self.assertAlmostEqual(shown_all, view.graphics_view.transform().m11())
        for name in ("more_btn", "fit_contexts_btn", "context_design_box", "context_density_box",
                     "reset_zoom_btn", "export_png_btn", "export_pdf_btn"):
            self.assertFalse(hasattr(dialog, name), name)
        dialog.close()

    def test_settings_dropdowns_drive_the_canvas(self) -> None:
        dialog = _viewer(
            None, rows=fixture_rows(4, 31), title="Gene neighborhoods: x",
            flank_sizes=(5, 10, 15, 30), initial_flank=15,
        )
        dialog.show()
        self.app.processEvents()
        self.assertEqual(dialog.color_by_box.currentData(), "annotation")
        self.assertEqual(
            [dialog.color_by_box.itemText(i) for i in range(dialog.color_by_box.count())],
            ["Gene type", "Gene family"],
        )
        self.assertEqual(dialog.view.color_mode, "annotation")
        self.assertEqual(dialog.labels_box.currentData(), "reference")
        self.assertEqual(dialog.labels_box.currentText(), "Key genes")
        self.assertEqual(
            [dialog.labels_box.itemText(i) for i in range(dialog.labels_box.count())],
            ["Key genes", "All", "Off"],
        )
        self.assertEqual(
            [dialog.flank_box.itemText(i) for i in range(dialog.flank_box.count())],
            ["±5", "±10", "±15", "±30"],
        )
        _choose(dialog.color_by_box, "family")
        self.assertEqual(dialog.view.color_mode, "family")
        _choose(dialog.labels_box, "off")
        self.assertEqual(dialog.view.label_mode, "off")
        # The keyboard steps through the choices.
        dialog.labels_box.setFocus()
        QtTest.QTest.keyClick(dialog.labels_box, QtCore.Qt.Key.Key_Up)
        self.assertEqual(dialog.view.label_mode, "all")
        _choose(dialog.flank_box, 5)
        self.app.processEvents()
        self.assertEqual(dialog._context_flank, 5)
        self.assertEqual(max(len(row.genes) for row in dialog.view.layout_result.rows), 11)
        dialog._context_density_actions["comfortable"].trigger()
        self.assertEqual(dialog.view._density, "comfortable")
        dialog._context_link_actions["off"].trigger()
        self.assertEqual(dialog.view.link_mode, "off")
        self.assertEqual(dialog.normalize_orientation_check.text(), "Mirror rows to match reference")
        dialog.close()

    def test_plasmid_context_flank_is_a_disabled_full_contig_item(self) -> None:
        rows = fixture_rows(2, 7)
        for row in rows:
            for value in row["genes"]:
                value.pop("is_center", None)
        dialog = _viewer(None, rows=rows, title="Plasmid context: pX1")
        box = dialog.flank_box
        self.assertEqual(box.count(), 1)
        self.assertEqual(box.currentText(), "Full contig")
        self.assertFalse(box.isEnabled())
        dialog.close()

    def test_title_block_chips_and_footer_hint(self) -> None:
        rows = fixture_rows()
        rows[3]["genes"][2]["row_data"]["ConsensusCategory"] = "mobile"
        rows[3]["genes"][2]["row_data"]["Gene"] = "IS26"
        dialog = _viewer(None, rows=rows, title="Gene neighborhoods: blaOXA-23")
        self.assertEqual(dialog.title_label.text(), "blaOXA-23")
        chips = {key: chip.text() for key, chip in dialog.title_chips.items()}
        self.assertEqual(chips["reference"], "Reference: REFERENCE, contig_0")
        self.assertIn("gene you opened", dialog.title_chips["reference"].toolTip())
        self.assertEqual(chips["genomes"], "4 genomes")
        self.assertEqual(chips["neighbors"], "4 with neighboring genes")
        self.assertEqual(chips["molecules"], "4 contigs: 2 plasmid, 2 chromosome")
        self.assertIn("IS26", chips["flanked"])
        self.assertTrue(chips["flanked"].startswith("Flanked by: "))
        self.assertEqual(dialog.footer_hint.text(), CONTEXT_FOOTER_HINT)
        self.assertIn("Arrow keys follow the family", CONTEXT_FOOTER_HINT)
        for text in list(chips.values()) + [CONTEXT_FOOTER_HINT]:
            self.assertNotIn("\u2014", text)
            self.assertNotIn("\u2013", text)
        dialog.close()

    def test_pin_shows_chip_and_labelled_fields(self) -> None:
        dialog = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x")
        dialog.show()
        self.app.processEvents()
        self.assertFalse(dialog.pin_chip.isVisible())
        gene = next(g for g in dialog.view.layout_result.rows[0].genes if g.gene.label == "gyrB")
        dialog.view._select_gene(gene)
        self.app.processEvents()
        self.assertTrue(dialog.pin_chip.isVisible())
        self.assertEqual(dialog.pin_chip_label.text(), "Pinned: gyrB")
        fields = {name: label.text() for name, label in dialog._detail_fields.items()}
        self.assertEqual(fields["Gene"], "gyrB")
        self.assertIn("contig_0", fields["Position"])
        self.assertEqual(fields["Gene type"], "Other")
        self.assertEqual(fields["Found in"], "4 of 4 genomes with neighbors")
        spy = QtTest.QSignalSpy(dialog.geneDoubleClicked)
        dialog.open_gene_card_btn.click()
        self.assertEqual(len(spy), 1)
        dialog.pin_chip_clear_btn.click()
        self.assertFalse(dialog.view.has_selection())
        self.assertFalse(dialog.pin_chip.isVisible())
        self.assertFalse(dialog.synteny_detail_panel.isVisible())
        dialog.close()

    def test_default_order_is_one_ranking_by_match(self) -> None:
        def gene(index, label, strand="+", center=False):
            return {
                "start": index * 1000,
                "end": index * 1000 + 800,
                "strand": strand,
                "label": label,
                "is_center": center,
                "row_data": {"Gene": label},
            }

        labels = [f"g{index}" for index in range(31)]
        labels[15] = "aph"
        reference = {
            "sample": "REF1", "contig": "1", "molecule": "chromosome", "is_source": True,
            "center_bp": 15400,
            "genes": [gene(i, label, center=(i == 15)) for i, label in enumerate(labels)],
        }
        same = dict(reference, sample="CHR2", is_source=False)
        plasmid = {
            "sample": "PL1", "contig": "48", "molecule": "plasmid", "center_bp": 400,
            "genes": [gene(0, "aph", center=True)],
        }
        dialog = _viewer(None, rows=[reference, plasmid, same], title="Order")
        dialog.show()
        self.app.processEvents()
        self.assertFalse(dialog.group_molecule_check.isChecked())
        self.assertEqual(dialog.group_molecule_check.text(), "Split plasmid and chromosome")
        self.assertEqual([row["sample"] for row in dialog._visible_rows], ["REF1", "CHR2", "PL1"])
        dialog.group_molecule_check.setChecked(True)
        self.assertEqual([row["sample"] for row in dialog._visible_rows], ["REF1", "PL1", "CHR2"])
        self.assertIn("plasmid contigs before chromosome", dialog.view.plot_caption.caption)
        dialog.close()

    def test_merge_has_one_name_and_header_is_plain(self) -> None:
        dialog = _viewer(None, rows=fixture_rows(), title="Merge")
        dialog.show()
        self.app.processEvents()
        self.assertEqual(dialog.collapse_btn.text(), "Merge identical neighborhoods")
        self.assertTrue(dialog.collapse_btn.isCheckable())
        # The chips describe the comparison; the status line only shows changes.
        self.assertEqual(dialog.info.text(), "")
        self.assertFalse(dialog.info.isVisible())
        dialog.collapse_btn.setChecked(True)
        self.assertEqual(dialog.collapse_btn.text(), "Merge identical neighborhoods")
        info = dialog.info.text()
        self.assertIn("merged into", info)
        for jargon in ("(s)", "|", "track", "signature", "\u2014", "\u2013"):
            self.assertNotIn(jargon, info)
        dialog.close()

    def test_flanked_by_counts_only_named_mobile_genes(self) -> None:
        for text in ("IS26", "ISAba1", "tnpA", "Tn3", "intI1", "IS6 family transposase", "integrase"):
            self.assertIsNotNone(_MOBILE_NAME_PATTERN.search(text), text)
        for text in ("is2", "isc1", "this is it", "spoIVCA", "site-specific recombinase"):
            self.assertIsNone(_MOBILE_NAME_PATTERN.search(text), text)
        rows = fixture_rows()
        payload = rows[1]["genes"][3]["row_data"]
        rows[1]["genes"][3]["label"] = "spoIVCA"
        payload.update({"Gene": "spoIVCA", "Product": "site-specific recombinase", "ConsensusCategory": "MOBILE"})
        dialog = _viewer(None, rows=rows, title="Gene neighborhoods: x")
        self.assertNotIn("spoIVCA", dialog.title_chips["flanked"].toolTip())
        self.assertIn("tnpA", dialog.title_chips["flanked"].text())
        dialog.close()

    def test_gene_row_on_a_plasmid_contig_names_its_contig(self) -> None:
        rows = fixture_rows()
        rows[0]["is_plasmid"] = True
        rows[0]["title"] = "REFERENCE | contig_0"
        dialog = _viewer(None, rows=rows, title="Gene neighborhoods: x")
        self.assertEqual(dialog.title_chips["reference"].text(), "Reference: REFERENCE, contig_0")
        dialog.close()

    def test_filtered_rows_tell_tooltips_they_count_rows_shown(self) -> None:
        dialog = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x")
        self.assertFalse(dialog.view.rows_filtered)
        dialog.sample_query.setText("sample_1")
        dialog._apply_filters()
        self.assertTrue(dialog.view.rows_filtered)
        gene = next(g for g in dialog.view.layout_result.rows[0].genes if not g.gene.is_focal)
        self.assertIn("shown", dialog.view.presence_text(gene))
        dialog.close()

    def test_classic_renderer_hides_switches_it_ignores(self) -> None:
        original = igv_collection_view.IgvCollectionView
        igv_collection_view.IgvCollectionView = LegacyIgvCollectionView
        try:
            dialog = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x")
        finally:
            igv_collection_view.IgvCollectionView = original
        self.assertIsInstance(dialog.view, LegacyIgvCollectionView)
        self.assertTrue(dialog.color_by_group.isHidden())
        self.assertTrue(dialog.labels_group.isHidden())
        self.assertFalse(dialog.publication_style_action.isVisible())
        self.assertFalse(dialog.include_short_contigs_action.isVisible())
        dialog.close()
        dialog = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x")
        self.assertFalse(dialog.color_by_group.isHidden())
        self.assertTrue(dialog.publication_style_action.isVisible())
        dialog.close()

    def test_no_top_command_bar_and_title_stays_put_when_settings_fold(self) -> None:
        dialog = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x")
        dialog.resize(1500, 700)
        dialog.show()
        self.app.processEvents()
        self.assertFalse(hasattr(dialog, "command_bar"))
        self.assertFalse(hasattr(dialog, "filter_toggle_btn"))
        self.assertEqual(dialog.findChildren(QtWidgets.QFrame, "ContextCommandBar"), [])
        self.assertLess(dialog.minimumSizeHint().width(), 1000)
        # Every button outside the plot lives in the Settings column, the
        # strip or the pinned gene panel.
        homes = (dialog.filter_rail, dialog.settings_strip, dialog.synteny_detail_panel, dialog.view)
        for button in dialog.findChildren(QtWidgets.QAbstractButton):
            self.assertTrue(any(home.isAncestorOf(button) for home in homes), button.text())
        title_geometry = dialog.title_label.parentWidget().geometry()
        dialog.settings_collapse_btn.click()
        self.app.processEvents()
        self.assertFalse(dialog.filter_rail.isVisible())
        self.assertEqual(dialog.title_label.parentWidget().geometry(), title_geometry)
        dialog.settings_expand_btn.click()
        self.app.processEvents()
        self.assertTrue(dialog.filter_rail.isVisible())
        self.assertEqual(dialog.title_label.parentWidget().geometry(), title_geometry)
        self.assertTrue(dialog.footer_hint.isVisible())
        self.assertGreater(dialog.footer_hint.mapTo(dialog, QtCore.QPoint(0, 0)).y(), dialog.filter_rail.mapTo(dialog, QtCore.QPoint(0, 0)).y())
        dialog.close()

    def test_actions_card_holds_zoom_and_export(self) -> None:
        dialog = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x", output_dir=Path(tempfile.gettempdir()))
        dialog.setStyleSheet(APP_STYLESHEET)
        dialog.resize(1600, 900)
        dialog.show()
        self.app.processEvents()
        card = dialog.actions_card
        texts = [b.text() for b in card.findChildren(QtWidgets.QAbstractButton) if b.isVisibleTo(card)]
        # No Phylogeny button: a Phylogeny panel of the window shows the tree.
        self.assertEqual(sorted(texts), sorted(["Fit reference", "Show all", "Export figure"]))
        self.assertFalse(hasattr(dialog, "phylogeny_btn"))
        self.assertFalse(hasattr(dialog, "tree_link_chip"))
        # The last card of the column.
        cards = [dialog.view_card, dialog.order_card, card]
        tops = [c.mapTo(dialog.filter_rail, QtCore.QPoint(0, 0)).y() for c in cards]
        self.assertEqual(tops, sorted(tops))
        fit, show_all = dialog.fit_reference_btn, dialog.show_all_btn
        self.assertEqual(fit.width(), show_all.width())
        self.assertEqual(fit.y(), show_all.y())
        self.assertLess(fit.x(), show_all.x())
        # Export spans the card, from Fit's left edge to Show all's right.
        self.assertEqual(dialog.export_btn.x(), fit.x())
        self.assertEqual(dialog.export_btn.geometry().right(), show_all.geometry().right())
        # The Export figure menu is unchanged.
        menu = dialog.export_btn.menu()
        self.assertIsNotNone(menu)
        texts = [action.text() for action in menu.actions()]
        for text in ("Figure as PNG (2x)", "Figure as PDF (vector)", "Publication style", "Include short contigs"):
            self.assertIn(text, texts)
        dialog.close()

    def test_pinned_chip_sits_in_the_plot_header_clear_of_the_focal_tab(self) -> None:
        for width, inspector in ((1600, False), (1180, True)):
            NeighbourhoodView.settings_panel_open = None
            dialog = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x")
            dialog.setStyleSheet(APP_STYLESHEET)
            dialog.resize(width, 700)
            dialog.show()
            self.app.processEvents()
            header = dialog.view.plot_caption
            self.assertIs(dialog.pin_chip.parentWidget(), header)
            gene = next(g for g in dialog.view.layout_result.rows[0].genes if g.gene.label == "gyrB")
            dialog.view._select_gene(gene)
            self.app.processEvents()
            self.assertTrue(dialog.pin_chip.isVisible())
            self.assertEqual(dialog.synteny_detail_panel.isVisible(), True)
            chip = QtCore.QRectF(dialog.pin_chip.geometry())
            self.assertTrue(QtCore.QRectF(header.rect()).contains(chip))
            tab = header._focal_tab_rect(header._tab_font(header.font()))
            if tab is not None and not header._corner_hides_tab:
                self.assertFalse(chip.intersects(tab))
            dialog.pin_chip_clear_btn.click()
            self.app.processEvents()
            self.assertFalse(dialog.view.has_selection())
            self.assertFalse(dialog.pin_chip.isVisible())
            dialog.close()

    def test_pinned_chip_moves_beside_a_focal_tab_at_the_right_end(self) -> None:
        dialog = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x")
        dialog.resize(1500, 700)
        dialog.show()
        self.app.processEvents()
        header = dialog.view.plot_caption
        dialog.pin_chip_label.setText("Pinned: gyrB")
        dialog.pin_chip.show()
        right_end = QtCore.QRectF(dialog.pin_chip.geometry())
        # A focal tab under the chip's right-end spot pushes the chip left of it.
        fake_tab = QtCore.QRectF(right_end.center().x() - 30.0, 4.0, 60.0, header.height() - 4.0)
        header._focal_tab_rect = lambda _font: fake_tab
        header.refresh()
        moved = QtCore.QRectF(dialog.pin_chip.geometry())
        self.assertFalse(moved.intersects(fake_tab))
        self.assertLess(moved.right(), fake_tab.left())
        self.assertFalse(header._corner_hides_tab)
        dialog.close()

    def test_pinned_chip_hides_the_tab_name_and_stays_inside_a_narrow_header(self) -> None:
        dialog = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x")
        dialog.setStyleSheet(APP_STYLESHEET)
        dialog.resize(1500, 700)
        dialog.show()
        self.app.processEvents()
        header = dialog.view.plot_caption
        # A long name is cut short; the tooltip keeps it in full.
        long_name = "WP_" + "0123456789" * 6
        dialog._show_gene_details({"__context_label__": long_name})
        self.app.processEvents()
        label = dialog.pin_chip_label
        self.assertTrue(label.text().startswith("Pinned: WP_"))
        self.assertTrue(label.text().endswith("\u2026"))
        self.assertLessEqual(label.fontMetrics().horizontalAdvance(label.text()), NeighbourhoodView.PIN_CHIP_MAX_TEXT_WIDTH)
        self.assertIn(long_name, dialog.pin_chip.toolTip())
        # A tab filling the strip leaves no room for both: the tab's name
        # yields and the chip stays whole, at or right of the strip start.
        fake_tab = QtCore.QRectF(context_canvas.RAIL_BADGE_X, 4.0, header.width() - 20.0, header.height() - 4.0)
        header._focal_tab_rect = lambda _font: fake_tab
        header.refresh()
        chip = dialog.pin_chip.geometry()
        self.assertTrue(header._corner_hides_tab)
        self.assertGreaterEqual(chip.left(), context_canvas.RAIL_BADGE_X)
        self.assertTrue(header.rect().contains(chip))
        self.assertEqual(header.tooltip_at(fake_tab.center().x()), header.caption)
        # Too narrow even for the chip: it shrinks to the strip, never off its left end.
        header.resize(120, header.height())
        header.refresh()
        chip = dialog.pin_chip.geometry()
        self.assertGreaterEqual(chip.left(), context_canvas.RAIL_BADGE_X)
        self.assertLessEqual(chip.right(), header.width())
        dialog.close()

    def test_settings_fold_to_a_slim_strip_and_ctrl_f_opens_them(self) -> None:
        dialog = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x")
        dialog.resize(1600, 800)
        dialog.show()
        self.app.processEvents()
        self.assertTrue(dialog.settings_open())
        self.assertFalse(dialog.settings_strip.isVisible())
        self.assertTrue(dialog.filter_rail.isAncestorOf(dialog.settings_collapse_btn))
        # The fold button sits at the right of the column header, level with the title.
        title = dialog.settings_title
        self.assertGreater(dialog.settings_collapse_btn.x(), title.geometry().right())
        self.assertLess(
            abs(dialog.settings_collapse_btn.geometry().center().y() - title.geometry().center().y()), 4
        )
        dialog.settings_collapse_btn.click()
        self.app.processEvents()
        self.assertFalse(dialog.settings_open())
        self.assertTrue(dialog.settings_strip.isVisible())
        self.assertEqual(dialog.settings_strip.width(), NeighbourhoodView.SETTINGS_STRIP_WIDTH)
        self.assertEqual(dialog.settings_expand_btn.toolTip(), "Show settings (Ctrl+F)")
        self.assertIs(NeighbourhoodView.settings_panel_open, False)
        shortcut = next(
            item for item in dialog.findChildren(QtGui.QShortcut)
            if item.key() == QtGui.QKeySequence("Ctrl+F")
        )
        shortcut.activated.emit()
        self.app.processEvents()
        self.assertTrue(dialog.settings_open())
        self.assertFalse(dialog.settings_strip.isVisible())
        self.assertTrue(dialog.sample_query.hasFocus())
        self.assertIs(NeighbourhoodView.settings_panel_open, True)
        # The strip's caption opens the column as well.
        dialog.settings_collapse_btn.click()
        self.app.processEvents()
        self.assertEqual(dialog.settings_strip_caption.text(), "Settings")
        dialog.settings_strip_caption.clicked.emit()
        self.app.processEvents()
        self.assertTrue(dialog.settings_open())
        dialog.close()

    def test_settings_captions_sit_above_full_width_controls(self) -> None:
        dialog = _viewer(
            None, rows=fixture_rows(), title="Gene neighborhoods: x", flank_sizes=(5, 10, 15, 30), initial_flank=15,
        )
        dialog.setStyleSheet(APP_STYLESHEET)
        dialog.resize(1600, 900)
        dialog.show()
        self.app.processEvents()
        rail = dialog.filter_rail
        self.assertFalse(rail.findChildren(QtWidgets.QFormLayout))

        def box(widget: QtWidgets.QWidget) -> QtCore.QRect:
            return QtCore.QRect(widget.mapTo(rail, QtCore.QPoint(0, 0)), widget.size())

        controls = {
            "Flank": dialog.flank_box,
            "Color by": dialog.color_by_box,
            "Labels": dialog.labels_box,
            "Molecule": dialog.molecule_box,
            "Sort by": dialog.sort_column_box,
            "Order": dialog.sort_order_box,
        }
        reference = box(dialog.sample_query)
        for caption_text, control in controls.items():
            group = control.parentWidget()
            captions = [
                label for label in group.findChildren(QtWidgets.QLabel)
                if label.property("role") == "settingsCaption"
            ]
            self.assertEqual([label.text() for label in captions], [caption_text])
            caption, field = box(captions[0]), box(control)
            self.assertLess(caption.bottom(), field.top(), caption_text)
            self.assertEqual(caption.left(), field.left(), caption_text)
            # Every control spans its card, like the search field.
            self.assertEqual(field.left(), reference.left(), caption_text)
            self.assertEqual(field.width(), reference.width(), caption_text)
        # One height for every control.
        height = NeighbourhoodView.SETTINGS_CONTROL_HEIGHT
        for widget in (
            dialog.sample_query, dialog.sort_column_box, dialog.sort_order_box, dialog.add_filter_btn,
            dialog.fit_reference_btn, dialog.export_btn, dialog.flank_box, dialog.color_by_box,
            dialog.labels_box, dialog.molecule_box,
        ):
            self.assertEqual(widget.height(), height, widget)
        # Add filter is a normal button; Clear sits on its line at the right.
        self.assertLess(dialog.add_filter_btn.width(), reference.width() // 2)
        self.assertEqual(box(dialog.add_filter_btn).left(), reference.left())
        clear, add = box(dialog.clear_filters_btn), box(dialog.add_filter_btn)
        self.assertLess(abs(clear.center().y() - add.center().y()), 2)
        self.assertLessEqual(abs(clear.right() - reference.right()), 1)
        # The Clear link's text ends flush with the controls' right edge.
        image = dialog.clear_filters_btn.grab().toImage()
        background = image.pixelColor(0, 0)
        ink = [
            x for x in range(image.width()) for y in range(image.height())
            if abs(image.pixelColor(x, y).red() - background.red()) > 40
        ]
        self.assertLessEqual(clear.right() - (clear.left() + max(ink)), 4)
        # No filters: the chip list takes no room under the button row.
        self.assertFalse(dialog.filters_container.isVisibleTo(dialog.view_card))
        self.assertLessEqual(
            dialog.view_card.geometry().bottom() - dialog.add_filter_btn.mapTo(dialog.view_card.parentWidget(), QtCore.QPoint(0, 0)).y()
            - dialog.add_filter_btn.height(),
            12 + 1,
        )
        dialog._metadata_filters["Country"] = ["Romania"]
        dialog._render_metadata_filters()
        self.app.processEvents()
        self.assertTrue(dialog.filters_container.isVisibleTo(dialog.view_card))
        self.assertGreater(box(dialog.filters_container).top(), add.bottom())
        dialog._clear_context_filters()
        self.app.processEvents()
        self.assertFalse(dialog.filters_container.isVisibleTo(dialog.view_card))
        # Checkboxes start at the captions' left edge.
        for check in (dialog.normalize_orientation_check, dialog.group_molecule_check, dialog.collapse_btn):
            self.assertEqual(box(check).left(), reference.left())
        dialog.close()

    def test_settings_column_fits_its_cards_without_sideways_scrolling(self) -> None:
        dialog = _viewer(
            None, rows=fixture_rows(), title="Gene neighborhoods: x", flank_sizes=(5, 10, 15, 30), initial_flank=15,
        )
        dialog.setStyleSheet(APP_STYLESHEET)
        dialog.resize(1600, 1100)
        dialog.show()
        self.app.processEvents()
        scroll = dialog.settings_scroll
        content = scroll.widget()
        self.assertEqual(scroll.horizontalScrollBar().maximum(), 0)
        self.assertLessEqual(content.minimumSizeHint().width(), scroll.viewport().width())
        margin = NeighbourhoodView.SETTINGS_MARGIN
        rail = dialog.filter_rail
        # A tall window fits every card: no stray vertical scroll bar either.
        self.assertFalse(scroll.verticalScrollBar().isVisible())
        cards = [dialog.view_card, dialog.order_card, dialog.actions_card]
        boxes = [QtCore.QRect(card.mapTo(rail, QtCore.QPoint(0, 0)), card.size()) for card in cards]
        for card, box in zip(cards, boxes):
            # Same margin on both sides: no card clipped at either edge.
            self.assertEqual(box.x(), margin)
            self.assertEqual(rail.width() - box.right() - 1, margin)
            self.assertEqual(card.width(), cards[0].width())
            self.assertEqual(card.layout().contentsMargins(), QtCore.QMargins(12, 12, 12, 12))
        gaps = [boxes[i + 1].y() - boxes[i].bottom() - 1 for i in range(len(boxes) - 1)]
        # Actions follows Order with the usual gap; the column's background
        # simply continues below it.
        self.assertEqual(gaps, [margin] * len(gaps))
        self.assertGreater(rail.height() - boxes[-1].bottom() - 1, margin)
        dialog.close()

    def test_actions_card_stays_in_view_at_short_window_heights(self) -> None:
        for size in ((1600, 900), (1180, 760)):
            NeighbourhoodView.settings_panel_open = True
            dialog = _viewer(
                None, rows=fixture_rows(), title="Gene neighborhoods: x", flank_sizes=(5, 10, 15, 30), initial_flank=15,
            )
            dialog.setStyleSheet(APP_STYLESHEET)
            dialog.resize(*size)
            dialog.show()
            self.app.processEvents()
            rail, scroll = dialog.filter_rail, dialog.settings_scroll
            # Pinned under the scrolling cards, not one of them.
            self.assertFalse(scroll.isAncestorOf(dialog.actions_card))
            self.assertTrue(dialog.settings_footer.isAncestorOf(dialog.actions_card))
            for button in (dialog.fit_reference_btn, dialog.show_all_btn, dialog.export_btn):
                box = QtCore.QRect(button.mapTo(rail, QtCore.QPoint(0, 0)), button.size())
                self.assertTrue(rail.rect().contains(box), (size, button.text()))
            bar = scroll.verticalScrollBar()
            if bar.isVisible():
                # The cards scroll under a rule, and the Actions card keeps
                # their right edge beside the scroll bar.
                self.assertIs(dialog.settings_footer.property("divider"), True)
                order_right = dialog.order_card.mapTo(rail, QtCore.QPoint(0, 0)).x() + dialog.order_card.width()
                actions_right = dialog.actions_card.mapTo(rail, QtCore.QPoint(0, 0)).x() + dialog.actions_card.width()
                self.assertEqual(order_right, actions_right)
            dialog.close()
        NeighbourhoodView.settings_panel_open = None

    def test_settings_dropdowns_menu_button_and_search_clear_icon(self) -> None:
        dialog = _viewer(
            None, rows=fixture_rows(), title="Gene neighborhoods: x", flank_sizes=(5, 10, 15, 30), initial_flank=15,
        )
        dialog.setStyleSheet(APP_STYLESHEET)
        dialog.resize(1600, 900)
        dialog.show()
        dialog.sample_query.setText("x")
        self.app.processEvents()
        # One shared drop-down style: no boxed arrow, a plain chevron image.
        self.assertIn("QComboBox::down-arrow", APP_STYLESHEET)
        self.assertNotIn("QWidget#SettingsContent QComboBox", APP_STYLESHEET)
        muted = QtGui.QColor(COLORS["text_muted"])
        for widget in (
            dialog.flank_box, dialog.color_by_box, dialog.labels_box, dialog.molecule_box,
            dialog.sort_column_box, dialog.sort_order_box, dialog.export_btn,
        ):
            # The plain chevron is drawn in the muted text color at the right.
            image = widget.grab().toImage()
            x0, y0 = widget.width() - 22, widget.height() // 2 - 4
            near = [
                image.pixelColor(x, y) for x in range(x0, x0 + 16) for y in range(y0, y0 + 8)
            ]
            self.assertTrue(
                any(abs(c.red() - muted.red()) + abs(c.green() - muted.green()) + abs(c.blue() - muted.blue()) < 60 for c in near),
                widget,
            )
        # Export reads like the other buttons: no bold, no default triangle.
        self.assertLess(dialog.export_btn.font().weight(), QtGui.QFont.Weight.DemiBold)
        self.assertIn("ContextExportButton::menu-indicator", APP_STYLESHEET)
        # The search field's clear icon is centered and inside the field.
        icons = [b for b in dialog.sample_query.findChildren(QtWidgets.QToolButton) if b.isVisible()]
        self.assertTrue(icons)
        for icon in icons:
            geometry = icon.geometry()
            self.assertGreaterEqual(geometry.top(), 0)
            self.assertLessEqual(geometry.bottom(), dialog.sample_query.height() - 1)
            self.assertLessEqual(abs(geometry.center().y() - dialog.sample_query.rect().center().y()), 1)
        dialog.close()

    def test_toggling_settings_keeps_zoom_and_rows_in_view(self) -> None:
        dialog = _viewer(None, rows=fixture_rows(40, 31), title="Gene neighborhoods: x")
        dialog.resize(1500, 700)
        dialog.show()
        QtTest.QTest.qWait(30)
        self.app.processEvents()
        view = dialog.view
        view.zoom_by(1.7)
        bar = view.graphics_view.verticalScrollBar()
        bar.setValue(bar.maximum() // 2)
        self.app.processEvents()
        self.assertGreater(bar.value(), 0)
        scale = view.graphics_view.transform().m11()
        vertical = view.graphics_view.verticalScrollBar().value()
        width = view.graphics_view.viewport().width()
        for button, checked in ((dialog.settings_collapse_btn, False), (dialog.settings_expand_btn, True)):
            button.click()
            QtTest.QTest.qWait(30)
            self.app.processEvents()
            self.assertEqual(dialog.settings_open(), checked)
            self.assertAlmostEqual(view.graphics_view.transform().m11(), scale)
            self.assertEqual(view.graphics_view.verticalScrollBar().value(), vertical)
            self.assertEqual(view.graphics_view.viewport().width() > width, not checked)
        # The choice lasts for the session: a new window opens the same way.
        dialog.settings_collapse_btn.click()
        dialog.close()
        again = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x")
        again.show()
        self.app.processEvents()
        self.assertFalse(again.settings_open())
        self.assertFalse(again.filter_rail.isVisible())
        self.assertTrue(again.settings_strip.isVisible())
        again.close()

    def test_settings_toggle_keeps_left_edge_with_focal_near_it(self) -> None:
        dialog = _viewer(None, rows=fixture_rows(40, 31), title="Gene neighborhoods: x")
        dialog.resize(1500, 700)
        dialog.show()
        QtTest.QTest.qWait(30)
        view = dialog.view
        graphics = view.graphics_view
        view.zoom_by(3.0)
        self.app.processEvents()
        # Focal column just inside the left edge, scrolled clear of the ends.
        band = view._focal_band_rect
        visible = graphics.mapToScene(graphics.viewport().rect()).boundingRect()
        shift = (band.left() - visible.left() - visible.width() * 0.11) * graphics.transform().m11()
        graphics.horizontalScrollBar().setValue(graphics.horizontalScrollBar().value() + int(shift))
        self.app.processEvents()
        left = view.view_state()[1]
        for checked in (False, True):
            dialog._toggle_filter_rail(checked)
            QtTest.QTest.qWait(30)
            self.app.processEvents()
            self.assertAlmostEqual(view.view_state()[1], left, delta=2.0)
        dialog.close()

    def test_settings_panel_default_follows_window_width(self) -> None:
        wide = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x")
        wide.resize(1600, 800)
        wide.show()
        self.app.processEvents()
        self.assertTrue(wide.filter_rail.isVisible())
        self.assertFalse(wide.settings_strip.isVisible())
        wide.close()
        narrow = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x")
        narrow.resize(1180, 800)
        narrow.show()
        self.app.processEvents()
        self.assertFalse(narrow.filter_rail.isVisible())
        self.assertTrue(narrow.settings_strip.isVisible())
        # The automatic fold is not a choice: the next wide window opens it.
        self.assertIsNone(NeighbourhoodView.settings_panel_open)
        narrow.settings_expand_btn.click()
        self.app.processEvents()
        self.assertTrue(narrow.filter_rail.isVisible())
        narrow.close()
        # Once chosen, the choice wins over the width.
        again = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x")
        again.resize(1180, 800)
        again.show()
        self.app.processEvents()
        self.assertTrue(again.filter_rail.isVisible())
        again.close()

    def test_settings_panel_cannot_be_dragged_shut(self) -> None:
        dialog = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x")
        dialog.resize(1600, 800)
        dialog.show()
        self.app.processEvents()
        dialog.main_split.setSizes([0, 1600])
        self.app.processEvents()
        self.assertGreater(dialog.main_split.sizes()[0], 0)
        dialog._focus_context_search()
        self.app.processEvents()
        self.assertTrue(dialog.settings_open())
        self.assertGreater(dialog.main_split.sizes()[0], 0)
        self.assertTrue(dialog.sample_query.isVisible())
        # One panel object per name, and merging shows its on state.
        self.assertEqual(len(dialog.findChildren(QtWidgets.QWidget, "FilterRail")), 1)
        self.assertIsInstance(dialog.collapse_btn, QtWidgets.QCheckBox)
        dialog.close()

    def test_no_tree_of_its_own(self) -> None:
        # The viewer never builds a tree copy, even with a run folder: a
        # Phylogeny panel of the window shows the tree instead.
        dialog = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x", output_dir=Path(tempfile.gettempdir()))
        dialog.show()
        self.app.processEvents()
        from annotation_explorer_polars.panaroo_views import PanCoreAccessoryTab

        self.assertEqual(dialog.findChildren(PanCoreAccessoryTab), [])
        for name in ("phylogeny_btn", "tree_link_chip", "pan_host", "panCopyCreated", "content_split"):
            self.assertFalse(hasattr(dialog, name), name)
        dialog.close()

    def test_pinning_keeps_the_focal_band_in_place(self) -> None:
        dialog = _viewer(None, rows=fixture_rows(4, 31), title="Gene neighborhoods: x")
        dialog.resize(1500, 700)
        dialog.show()
        self.app.processEvents()
        QtTest.QTest.qWait(20)
        view = dialog.view
        # Zoomed in, so the plot scrolls instead of centering the scene.
        view.zoom_by(2.0)
        self.app.processEvents()
        band = view._focal_band_rect
        before = view.graphics_view.mapFromScene(band.center()).x()
        focal = next(g for g in view.layout_result.rows[0].genes if g.gene.is_focal)
        view._select_gene(focal)
        self.app.processEvents()
        self.assertTrue(dialog.synteny_detail_panel.isVisible())
        self.assertEqual(view.graphics_view.mapFromScene(band.center()).x(), before)
        dialog.close()

    def test_plasmid_context_header_has_no_focal_gene_or_contig_key(self) -> None:
        rows = []
        for index in range(3):
            rows.append({
                "sample": f"S{min(index, 1)}",
                "contig": f"S{index}:primary:12",
                "title": f"S{min(index, 1)} | pX1",
                "display_label": f"S{min(index, 1)} | pX1",
                "molecule": "plasmid",
                "is_plasmid": True,
                "cluster_key": "primary:12",
                "genes": [
                    {"start": g * 900, "end": g * 900 + 700, "strand": "+", "label": f"g{g}",
                     "row_data": {"Gene": f"g{g}"}}
                    for g in range(4)
                ],
            })
        dialog = _viewer(None, rows=rows, title="Plasmid context")
        dialog.show()
        self.app.processEvents()
        chips = {key: chip.text() for key, chip in dialog.title_chips.items()}
        self.assertEqual(chips["reference"], "Reference: S0, plasmid pX1")
        self.assertEqual(chips["genomes"], "2 genomes")
        self.assertEqual(chips["molecules"], "3 plasmids")
        self.assertNotIn("neighbors", chips)
        self.assertNotIn("flanked", chips)
        for text in chips.values():
            self.assertNotIn("primary:12", text)
        # The rail names the plasmid as the chip does, not its instance key.
        dialog.show()
        self.app.processEvents()
        reference = dialog.view.layout_result.rows[0]
        self.assertEqual(_secondary_row_label(reference), "plasmid pX1")
        dialog.collapse_btn.setChecked(True)
        tooltip = dialog._collapse_rows(dialog._visible_rows)[0]["meta_tooltip"]
        self.assertIn("pX1", tooltip)
        self.assertNotIn("primary:12", tooltip)
        dialog.close()

    def test_export_file_named_after_window_title(self) -> None:
        dialog = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: blaOXA-23")
        path = dialog._default_export_path("png")
        self.assertTrue(path.endswith("bactolomew_gene_neighborhoods_blaoxa_23.png"))
        dialog.close()

    def test_inspector_shows_minus_strand(self) -> None:
        dialog = _viewer(None, rows=fixture_rows(), title="Strand")
        dialog.show()
        self.app.processEvents()
        row = dialog.view.layout_result.rows[0]
        minus = next(gene for gene in row.genes if gene.gene.strand == "-")
        dialog.view._select_gene(minus)
        self.assertEqual(dialog._detail_fields["Strand"].text(), "minus")
        dialog.close()

    def test_escape_clears_selection_before_closing_from_anywhere(self) -> None:
        # In a container that closes on Esc (a dialog), as the old window shim did.
        host = QtWidgets.QDialog()
        host_layout = QtWidgets.QVBoxLayout(host)
        dialog = NeighbourhoodView(host, rows=fixture_rows(), title="Esc")
        host_layout.addWidget(dialog)
        host.resize(1760, 860)
        host.show()
        self.app.processEvents()
        dialog.view._select_gene(dialog.view.layout_result.rows[0].genes[0])
        self.assertTrue(dialog.synteny_detail_panel.isVisible())
        self.assertTrue(dialog.pin_chip.isVisible())
        # Focus outside the plot: the viewer must still clear first.
        dialog.fit_reference_btn.setFocus()
        QtTest.QTest.keyClick(dialog.fit_reference_btn, QtCore.Qt.Key.Key_Escape)
        self.assertTrue(host.isVisible())
        self.assertFalse(dialog.view.has_selection())
        self.assertFalse(dialog.synteny_detail_panel.isVisible())
        self.assertFalse(dialog.pin_chip.isVisible())
        # From inside the plot as well.
        dialog.view._select_gene(dialog.view.layout_result.rows[0].genes[0])
        QtTest.QTest.keyClick(dialog.view.graphics_view, QtCore.Qt.Key.Key_Escape)
        self.assertTrue(host.isVisible())
        self.assertFalse(dialog.synteny_detail_panel.isVisible())
        # Nothing left to clear: the Esc passes to the container.
        QtTest.QTest.keyClick(dialog.view.graphics_view, QtCore.Qt.Key.Key_Escape)
        self.assertFalse(host.isVisible())
        host.deleteLater()

    def test_plasmid_constructor_without_focal_metadata(self) -> None:
        rows = fixture_rows(2, 7)
        for row in rows:
            row.pop("center_bp", None)
            row.pop("is_source", None)
            for value in row["genes"]:
                value.pop("is_center", None)
        dialog = _viewer(None, rows=rows, title="Plasmid context")
        dialog.show()
        self.app.processEvents()
        self.assertEqual(len(dialog.view.layout_result.rows), 2)
        self.assertIsNotNone(dialog.view.layout_result.rows[0].row.focal_gene)
        dialog.close()

    def test_missing_pandas_values_do_not_break_tooltip_or_inspector(self) -> None:
        rows = fixture_rows(2, 7)
        rows[0]["molecule"] = pd.NA
        rows[0]["molecule_kind"] = "chromosome"
        rows[0]["genes"][0]["label"] = "membrane enzyme"
        rows[0]["genes"][0]["row_data"]["Gene"] = pd.NA
        rows[0]["genes"][0]["row_data"]["Product"] = pd.NA
        dialog = _viewer(None, rows=rows, title="Plasmid context")
        dialog.show()
        self.app.processEvents()
        target = dialog.view.layout_result.rows[0].genes[0]
        tooltip = dialog.view._gene_tooltip(target, dialog.view.layout_result.rows[0])
        self.assertIn("Found in", tooltip)
        dialog.view._select_gene(target)
        self.app.processEvents()
        self.assertTrue(dialog.synteny_detail_panel.isVisible())
        dialog.close()

    def test_settings_panel_has_view_and_order_cards(self) -> None:
        dialog = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x")
        dialog.show()
        self.app.processEvents()
        self.assertTrue(dialog.filter_rail.isVisible())
        rail = dialog.filter_rail
        for widget in (
            dialog.flank_box,
            dialog.color_by_box,
            dialog.labels_box,
            dialog.normalize_orientation_check,
            dialog.group_molecule_check,
            dialog.collapse_btn,
            dialog.molecule_box,
        ):
            self.assertTrue(dialog.view_card.isAncestorOf(widget))
            self.assertTrue(rail.isAncestorOf(widget))
        cards = [dialog.view_card, dialog.order_card, dialog.actions_card]
        tops = [card.mapTo(rail, QtCore.QPoint(0, 0)).y() for card in cards]
        self.assertEqual(tops, sorted(tops))
        self.assertIs(dialog.sample_query.parent(), dialog.view_card)
        self.assertTrue(dialog.order_card.isAncestorOf(dialog.sort_column_box))
        self.assertEqual(dialog.sort_column_box.itemText(0), "Match to reference")
        molecule = dialog.molecule_box
        self.assertEqual([molecule.itemText(i) for i in range(molecule.count())], ["All", "Plasmid", "Chromosome"])
        self.assertEqual(dialog.clear_filters_btn.text(), "Clear")
        _choose(molecule, "plasmid")
        self.assertTrue(all(r["molecule"] == "plasmid" for r in dialog._visible_rows))
        self.assertEqual(dialog.settings_filter_badge.text(), "1 filter")
        self.assertEqual(dialog.settings_strip_badge.text(), "1")
        self.assertIn("1 filter on", dialog.settings_expand_btn.toolTip())
        dialog._clear_context_filters()
        self.assertTrue(dialog.settings_filter_badge.isHidden())
        self.assertEqual(dialog.settings_expand_btn.toolTip(), "Show settings (Ctrl+F)")
        self.assertEqual(molecule.currentData(), "")
        self.assertEqual(len(dialog._visible_rows), 4)
        dialog.close()

    def test_plasmids_first_draws_labelled_dividers(self) -> None:
        dialog = _viewer(None, rows=short_contig_rows() + fixture_rows()[2:], title="Gene neighborhoods: x")
        self.assertEqual(dialog.view.layout_result.dividers, ())
        dialog.group_molecule_check.setChecked(True)
        texts = [divider.text for divider in dialog.view.layout_result.dividers]
        self.assertEqual(texts, ["Plasmid contigs, 2 (+1 folded)", "Chromosome contigs, 1 (+2 folded)"])
        for text in texts:
            self.assertNotIn("\u2014", text)
        dialog.close()

    def test_merge_keeps_molecule_counts_blocks_and_shared_caps(self) -> None:
        rows = short_contig_rows() + fixture_rows()[2:]
        dialog = _viewer(None, rows=rows, title="Gene neighborhoods: x")
        dialog.collapse_btn.setChecked(True)
        # Every focal-only contig merges into one row, of mixed molecules.
        self.assertIn("(2 chromosome, 1 plasmid)", dialog.view.fold_sentence())
        merged = [r for r in dialog._collapse_rows(dialog._visible_rows) if r.get("collapsed_track_count") == 3][0]
        self.assertEqual(merged["molecule"], "mixed")
        self.assertEqual(merged["collapsed_molecule_counts"], {"plasmid": 1, "chromosome": 2})
        dialog.group_molecule_check.setChecked(True)
        texts = [divider.text for divider in dialog.view.layout_result.dividers]
        self.assertEqual(texts, ["Plasmid contigs, 2 (+1 folded)", "Chromosome contigs, 1 (+2 folded)"])
        a = {"sample": "a", "contig": "1", "contig_ends_low": True, "contig_ends_high": True,
             "genes": [{"start": 1, "end": 9, "label": "x", "is_center": True}]}
        b = dict(a, sample="b", contig_ends_high=False)
        collapsed = dialog._collapse_rows([a, b])
        self.assertEqual((collapsed[0]["contig_ends_low"], collapsed[0]["contig_ends_high"]), (True, False))
        dialog.close()

    def test_split_after_merge_keeps_each_molecule_in_one_block(self) -> None:
        dialog = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x")
        rows = [
            {"sample": "ref", "molecule": "chromosome", "is_comparison_reference": True},
            {"sample": "p1", "molecule": "plasmid"},
            {"sample": "m", "molecule": "mixed"},
            {"sample": "p2", "molecule": "plasmid"},
            {"sample": "c1", "molecule": "chromosome"},
        ]
        self.assertEqual(
            [row["sample"] for row in dialog._split_by_molecule(rows)],
            ["ref", "p1", "p2", "m", "c1"],
        )
        dialog.close()

    def test_export_keeps_a_pinned_short_row_gene_and_stays_quiet(self) -> None:
        canvas = ContextCanvas()
        canvas.resize(1100, 500)
        canvas.set_rows(short_contig_rows())
        canvas.set_short_rows_expanded(True)
        short = [row for row in canvas.layout_result.rows if row.short][0]
        canvas._select_gene(short.genes[0])
        emitted = []
        canvas.geneSelected.connect(lambda payload: emitted.append(payload))
        canvas.set_export_include_short_contigs(False)
        with tempfile.TemporaryDirectory() as folder:
            self.assertTrue(canvas.export_png(Path(folder) / "figure.png", scale=1.0))
        self.assertIsNotNone(canvas._selected_gene)
        self.assertEqual(canvas._selected_gene.gene.label, short.genes[0].gene.label)
        self.assertEqual(emitted, [])
        canvas.close()

    def test_second_locus_on_one_contig_names_its_position(self) -> None:
        rows = fixture_rows(3)
        rows[2].update(sample=rows[1]["sample"], contig=rows[1]["contig"])
        canvas = ContextCanvas()
        canvas.set_rows(rows)
        texts = [_secondary_row_label(row) for row in canvas.layout_result.rows]
        self.assertEqual(texts[0], "contig_0")
        self.assertRegex(texts[1], r"^contig_1 at [0-9,.]+ kb$")
        self.assertNotEqual(texts[1], texts[2])
        canvas.close()

    def test_merged_row_is_named_by_its_count_and_partial_reference_keeps_its_hint(self) -> None:
        rows = fixture_rows(2)
        rows[1].update(collapsed_samples=["a", "b"], collapsed_sample_count=2, display_label="2 genomes")
        rows[0].update(context_partial=True)
        canvas = ContextCanvas()
        canvas.resize(1200, 500)
        canvas.show()
        canvas.set_rows(rows)
        self.app.processEvents()
        reference, merged = canvas.layout_result.rows
        self.assertEqual(_primary_row_label(merged.row), "2 genomes")
        self.assertNotIn("count", _rail_row_geometry(merged, 30.0, 280.0, True))
        ref_parts = _rail_row_geometry(reference, 30.0, 280.0, True)
        self.assertIn("partial", ref_parts)
        rail = canvas.row_header
        y = float(canvas.graphics_view.mapFromScene(QtCore.QPointF(0.0, reference.center_y)).y())
        ref_parts = _rail_row_geometry(reference, y, float(rail.width()), True)
        partial = ref_parts["partial"]
        self.assertIn("every match covers only", rail.tooltip_at(partial.center().x(), partial.center().y()))
        self.assertIn("every match covers only the side it shows", _row_sentence(reference))
        canvas.close()

    def test_gene_card_rows_count_genes_on_each_side_of_the_focal_gene(self) -> None:
        from annotation_explorer_polars.annotation_detail_dialog import AnnotationDetailDialog

        genes = ["dnaA", "aph", "gyrB", "recA", "tnpA"]
        ann = pd.DataFrame([
            {"Sample": "S1", "Contig": "7", "Gene": name, "Product": name, "AnnotationID": f"A{index}",
             "Start": 100 + index * 1000, "Stop": 900 + index * 1000, "StartInt": 100 + index * 1000,
             "StopInt": 900 + index * 1000, "Strand": "+", "MobSuiteIsPlasmid": "True"}
            for index, name in enumerate(genes)
        ])
        focal = ann.iloc[1].to_dict()
        dialog = AnnotationDetailDialog(None, focal, list(focal), ann_df=ann)
        rows = dialog._build_context_collection_rows(window_genes=30)
        dialog.close()
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual((row["contig_genes_before"], row["contig_genes_after"]), (1, 3))
        self.assertEqual(row["title"], "S1, contig 7")
        self.assertEqual(row["molecule"], "plasmid")

    def test_use_as_reference_reranks_and_updates_the_chip(self) -> None:
        dialog = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x")
        target = dict(dialog.view.layout_result.rows[2].row.source)
        dialog._use_row_as_reference(target)
        first = dialog.view.layout_result.rows[0].row
        self.assertEqual(first.sample, target["sample"])
        self.assertTrue(first.is_reference)
        self.assertEqual(sum(row.row.is_reference for row in dialog.view.layout_result.rows), 1)
        self.assertIn(target["sample"], dialog.title_chips["reference"].text())
        self.assertIn("Use as reference", dialog.title_chips["reference"].toolTip())
        dialog.close()

    def test_row_menu_offers_open_in_table_only_with_a_callback(self) -> None:
        def menu_texts(dialog):
            texts = []
            original = QtWidgets.QMenu.exec
            QtWidgets.QMenu.exec = lambda menu, *args: texts.extend(a.text() for a in menu.actions())
            try:
                dialog._show_row_menu(dict(dialog.view.layout_result.rows[1].row.source), QtCore.QPoint())
            finally:
                QtWidgets.QMenu.exec = original
            return texts

        dialog = _viewer(None, rows=fixture_rows(), title="Gene neighborhoods: x")
        self.assertEqual(menu_texts(dialog), ["Use as reference", "Copy accession"])
        dialog.close()
        opened = []
        dialog = _viewer(
            None, rows=fixture_rows(), title="Gene neighborhoods: x", open_samples_callback=opened.append
        )
        self.assertEqual(menu_texts(dialog), ["Use as reference", "Copy accession", "Open sample in table"])
        dialog._open_samples(["sample_1"])
        self.assertEqual(opened, [["sample_1"]])
        dialog.close()

    def test_plasmid_context_reference_has_most_neighboring_genes(self) -> None:
        rows = []
        for index, count in enumerate((2, 5, 3)):
            rows.append({
                "sample": f"S{index}", "contig": f"S{index}:primary:1",
                "title": f"S{index}, plasmid pX1", "molecule": "plasmid", "is_plasmid": True,
                "cluster_key": "primary:1", "joined_contigs": True,
                "genes": [
                    {"start": g * 900, "end": g * 900 + 700, "strand": "+", "label": f"g{g}",
                     "row_data": {"Gene": f"g{g}"}}
                    for g in range(count)
                ],
            })
        dialog = _viewer(None, rows=rows, title="Plasmid context: pX1")
        self.assertEqual(dialog.title_chips["reference"].text(), "Reference: S1, plasmid pX1")
        self.assertEqual(dialog.view.layout_result.rows[0].row.sample, "S1")
        dialog.close()

    def test_plain_plasmid_names_have_no_bar_separators(self) -> None:
        from annotation_explorer_polars.plasmid_model import plain_plasmid_name

        self.assertEqual(plain_plasmid_name("p12 | IncFII | MOBF"), "p12 (IncFII, MOBF)")
        self.assertEqual(plain_plasmid_name("p12"), "p12")
        self.assertEqual(plain_plasmid_name(None), "")


if __name__ == "__main__":
    unittest.main()
