from __future__ import annotations

import os
import unittest
from types import SimpleNamespace

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

from PyQt6 import QtCore, QtGui, QtTest, QtWidgets

from annotation_explorer_polars.explorer_ui import _FilterMultiSelectComboBox
from annotation_explorer_polars.filters_controller import FiltersControllerMixin
from annotation_explorer_polars.model.filter_state import FilterState
from annotation_explorer_polars.panels.annotations import (
    AnnotationsPanel,
    build_view_mode_switch,
    build_window_options,
)
from annotation_explorer_polars.panels.sample_tree import SampleTreePanel
from annotation_explorer_polars.settings_column import settings_card, settings_control, stacked_field
from annotation_explorer_polars import theme
from annotation_explorer_polars.theme import apply_application_theme


class _FilterControlHarness(FiltersControllerMixin, QtWidgets.QWidget):
    def __init__(self) -> None:
        super().__init__()
        # The rail reads the run's predictors from self.run (RunData).
        self.run = SimpleNamespace(predictor_options=["AMRFinderPlus", "Bakta", "eggNOG"], tag_databases={})
        self._predictor_model = None
        self._tag_model = None
        self.refresh_count = 0

        layout = QtWidgets.QVBoxLayout(self)
        self.predictor_box = _FilterMultiSelectComboBox(self)
        self.predictor_box.view().pressed.connect(self._toggle_predictor_item)
        self._init_predictor_model()
        layout.addWidget(self.predictor_box)

        self.tag_box = _FilterMultiSelectComboBox(self)
        self.tag_box.view().pressed.connect(self._toggle_tag_item)
        self._init_tag_model()
        layout.addWidget(self.tag_box)

    def _rail_filters_changed(self, *, debounce: bool = False) -> None:
        # The window hands the new draft to its context; count the edits.
        self.refresh_count += 1


class AnnotationFilterControlTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        apply_application_theme(cls.app)

    def setUp(self) -> None:
        self.widget = _FilterControlHarness()
        self.widget.resize(280, 120)
        self.widget.show()
        self.app.processEvents()

    def tearDown(self) -> None:
        self.widget.predictor_box.hidePopup()
        self.widget.tag_box.hidePopup()
        self.widget.close()
        self.app.processEvents()

    def _open_from_field_body(self, box: _FilterMultiSelectComboBox) -> None:
        editor = box.lineEdit()
        self.assertIsNotNone(editor)
        QtTest.QTest.mouseClick(
            editor,
            QtCore.Qt.MouseButton.LeftButton,
            pos=editor.rect().center(),
        )
        self.app.processEvents()
        self.assertTrue(box.view().window().isVisible())

    def _click_popup_row(self, box: _FilterMultiSelectComboBox, row: int) -> None:
        index = box.model().index(row, 0)
        rect = box.view().visualRect(index)
        self.assertTrue(rect.isValid())
        QtTest.QTest.mouseClick(
            box.view().viewport(),
            QtCore.Qt.MouseButton.LeftButton,
            pos=rect.center(),
        )
        self.app.processEvents()

    def test_predictor_and_tag_fields_open_from_the_visible_body(self) -> None:
        self.assertEqual(self.widget.predictor_box.model().rowCount(), 4)
        self.assertEqual(self.widget.tag_box.model().rowCount(), 8)  # "Any (clear)" plus the 7 tags

        self._open_from_field_body(self.widget.predictor_box)
        self.widget.predictor_box.hidePopup()
        self.app.processEvents()
        self._open_from_field_body(self.widget.tag_box)

    def test_popup_rows_update_multiselect_state_and_summary(self) -> None:
        self._open_from_field_body(self.widget.predictor_box)
        self._click_popup_row(self.widget.predictor_box, 1)
        self.assertEqual(self.widget._selected_predictors(), ["AMRFinderPlus"])
        self.assertEqual(self.widget.predictor_box.lineEdit().text(), "Predictor (AMRFinderPlus)")

        self.widget.predictor_box.hidePopup()
        self.app.processEvents()
        self._open_from_field_body(self.widget.tag_box)
        self._click_popup_row(self.widget.tag_box, 1)
        self.assertEqual(self.widget._selected_tags(), ["A"])
        self.assertEqual(self.widget.tag_box.lineEdit().text(), "Tags (A)")
        self.assertGreaterEqual(self.widget.refresh_count, 2)

    def test_keyboard_opens_the_multiselect_popup(self) -> None:
        editor = self.widget.predictor_box.lineEdit()
        self.assertIsNotNone(editor)
        editor.setFocus()
        QtTest.QTest.keyClick(editor, QtCore.Qt.Key.Key_Space)
        self.app.processEvents()
        self.assertTrue(self.widget.predictor_box.view().window().isVisible())


    def test_predictor_list_is_wide_enough_for_long_names(self) -> None:
        box = self.widget.predictor_box
        model = box.model()
        long_name = "CGE-PlasmidFinder_Enterobacteriales_Prot_PLASMID"
        item = QtGui.QStandardItem(long_name)
        item.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
        item.setData(QtCore.Qt.CheckState.Unchecked, QtCore.Qt.ItemDataRole.CheckStateRole)
        model.appendRow(item)
        box.showPopup()
        self.app.processEvents()
        text_width = box.view().fontMetrics().horizontalAdvance(long_name)
        self.assertGreater(box.view().viewport().width(), text_width)


class TableViewDropdownTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        apply_application_theme(cls.app)

    def test_table_shows_dropdown_drives_and_follows_the_view_mode(self) -> None:
        # The "Table shows" box belongs to its Annotations panel (Phase 10:
        # the window's hidden view_mode_box combo is gone).
        strip, box = build_view_mode_switch()
        table = QtWidgets.QTableView()
        panel = AnnotationsPanel(strip, table, table=table, mode_box=box)
        changes = []
        panel.viewModeChanged.connect(changes.append)
        self.assertIs(panel.mode_box(), box)
        self.assertEqual(
            [box.itemText(i) for i in range(box.count())], ["Contig", "Around matches", "All matches"]
        )
        self.assertEqual(box.currentData(), "Contig")
        box.setCurrentIndex(box.findData("Global"))
        self.assertEqual(panel.view_mode(), "Global")
        self.assertEqual(changes, ["Global"])
        self.assertEqual(box.toolTip(), "Table: matching annotations across all samples")
        panel.set_view_mode("Context")
        self.assertEqual(box.currentText(), "Around matches")
        self.assertEqual(changes, ["Global", "Context"])
        # A restored mode shows in the box without a redraw request.
        panel.restore_state({"view_mode": "Contig"})
        self.assertEqual(box.currentData(), "Contig")
        self.assertEqual(changes, ["Global", "Context"])
        self.assertEqual(panel.save_state(), {"view_mode": "Contig"})
        panel.deleteLater()

    def test_window_filter_checkbox_sits_in_the_strip_and_follows_the_draft(self) -> None:
        # Feedback 2026-10-07: "Hide hypothetical" moved from the rail to the
        # Annotations panel ("Hide contigs without hits" moved on to the
        # Sample tree). It is a window filter: the panel shows ctx.draft and
        # reports user toggles.
        class _Ctx(QtCore.QObject):
            draftChanged = QtCore.pyqtSignal(object)

            def __init__(self) -> None:
                super().__init__()
                self.draft = FilterState(hide_empty_contigs=True)

        strip, box = build_view_mode_switch()
        options = build_window_options(strip)
        self.assertEqual(list(options), ["hide_hypothetical"])
        hypo = strip.findChild(QtWidgets.QCheckBox, "HideHypotheticalCheck")
        self.assertIs(hypo, options["hide_hypothetical"])
        self.assertIsNone(strip.findChild(QtWidgets.QCheckBox, "HideEmptyContigsCheck"))
        self.assertEqual(hypo.text(), "Hide hypothetical")
        # The tooltip says this is a window filter, not an option of the table.
        self.assertIn("any panel of this window", hypo.toolTip())
        self.assertIn("still list them", hypo.toolTip())
        self.assertFalse(hasattr(AnnotationsPanel, "window_options"))
        # One row: caption, box, the checkbox; it wraps when narrow, so the
        # checkbox never forces a wide panel.
        layout = strip.layout()
        order = [layout.itemAt(i).widget() for i in range(layout.count())]
        self.assertEqual(order[1:], [box, hypo])
        one_line = strip.sizeHint().height()
        self.assertLess(strip.minimumSizeHint().width(), 300)
        self.assertGreater(layout.heightForWidth(180), one_line)
        self.assertEqual(layout.heightForWidth(2000), one_line)

        table = QtWidgets.QTableView()
        panel = AnnotationsPanel(strip, table, table=table, mode_box=box, window_options=options)
        self.assertIs(panel.hide_hypo_chk, hypo)
        self.assertFalse(hasattr(panel, "hide_empty_contigs_chk"))
        toggles = []
        panel.windowOptionToggled.connect(lambda field, value: toggles.append((field, value)))
        ctx = _Ctx()
        panel.bind(ctx, "annotations")
        # Bound: shows the draft at once, without reporting a toggle.
        self.assertFalse(hypo.isChecked())
        self.assertEqual(toggles, [])
        # The draft changes elsewhere (another instance, a cloned window).
        ctx.draft = FilterState(hide_hypothetical=True)
        ctx.draftChanged.emit(ctx.draft)
        self.assertTrue(hypo.isChecked())
        self.assertEqual(toggles, [])
        # The user toggles: the panel reports it, the window writes the draft.
        hypo.click()
        self.assertEqual(toggles, [("hide_hypothetical", False)])
        # Window filters are not panel options: save_state is the view mode only.
        self.assertEqual(panel.save_state(), {"view_mode": "Contig"})
        panel.unbind()
        ctx.draftChanged.emit(FilterState(hide_hypothetical=True))
        self.assertFalse(hypo.isChecked())
        panel.deleteLater()

    def test_sample_tree_panel_has_hide_contigs_checkbox_that_follows_the_draft(self) -> None:
        # Feedback 2026-10-07: "Hide contigs without hits" sits on the Sample
        # tree panel (the only panel it changes), in a strip above the tree.
        class _Ctx(QtCore.QObject):
            draftChanged = QtCore.pyqtSignal(object)

            def __init__(self) -> None:
                super().__init__()
                self.draft = FilterState(hide_empty_contigs=True)

        class _Window:
            def _build_sample_tree_panel(self) -> QtWidgets.QWidget:
                holder = QtWidgets.QWidget()
                self.tree = QtWidgets.QTreeWidget(holder)
                QtWidgets.QVBoxLayout(holder).addWidget(self.tree)
                return holder

        panel = SampleTreePanel(_Window())
        check = panel.findChild(QtWidgets.QCheckBox, "HideEmptyContigsCheck")
        self.assertIs(check, panel.hide_empty_contigs_chk)
        self.assertEqual(check.text(), "Hide contigs without hits")
        self.assertIn("Window filter", check.toolTip())
        self.assertIn("this tree", check.toolTip())
        self.assertIsNone(panel.findChild(QtWidgets.QCheckBox, "HideHypotheticalCheck"))
        # The strip sits above the tree.
        layout = panel.layout()
        self.assertIs(layout.itemAt(0).widget(), panel.options_strip)
        self.assertEqual(panel.options_strip.objectName(), "SampleTreeOptions")
        self.assertIs(check.parentWidget(), panel.options_strip)
        toggles = []
        panel.windowOptionToggled.connect(lambda field, value: toggles.append((field, value)))
        ctx = _Ctx()
        panel.bind(ctx, "sample_tree")
        self.assertTrue(check.isChecked())
        self.assertEqual(toggles, [])
        ctx.draft = FilterState(hide_empty_contigs=False)
        ctx.draftChanged.emit(ctx.draft)
        self.assertFalse(check.isChecked())
        self.assertEqual(toggles, [])
        check.click()
        self.assertEqual(toggles, [("hide_empty_contigs", True)])
        self.assertEqual(panel.save_state(), {})
        panel.unbind()
        ctx.draftChanged.emit(FilterState(hide_empty_contigs=False))
        self.assertTrue(check.isChecked())
        panel.deleteLater()

    def test_stacked_field_puts_a_muted_caption_above_its_control(self) -> None:
        # The rail's fields (Sort by, Order, the annotation filters).
        card, layout = settings_card(None, "Sorting")
        combo = settings_control(QtWidgets.QComboBox())
        field = stacked_field(card, layout, "Sort by", combo)
        card.resize(240, 120)
        card.show()
        self.app.processEvents()
        label = next(lab for lab in field.findChildren(QtWidgets.QLabel) if lab.text() == "Sort by")
        self.assertEqual(label.property("role"), "settingsCaption")
        self.assertIs(label.buddy(), combo)
        self.assertLess(label.geometry().bottom(), combo.geometry().top())
        self.assertEqual(combo.height(), 30)
        card.close()


if __name__ == "__main__":
    unittest.main()


class SharedComboStyleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        apply_application_theme(cls.app)

    def test_style_sheet_images_exist(self) -> None:
        for path in theme._IMAGES.values():
            self.assertTrue(os.path.isfile(path), path)
            self.assertIn(f'url("{path}")', theme.APP_STYLESHEET)

    def test_wheel_scrolls_the_page_instead_of_an_unfocused_combo(self) -> None:
        area = QtWidgets.QScrollArea()
        content = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(content)
        other = QtWidgets.QLineEdit()
        layout.addWidget(other)
        combos = []
        for _ in range(30):
            combo = QtWidgets.QComboBox()
            combo.addItems(["a", "b", "c"])
            layout.addWidget(combo)
            combos.append(combo)
        area.setWidget(content)
        area.resize(200, 200)
        area.show()
        other.setFocus()
        self.app.processEvents()
        combo = combos[2]
        self.assertEqual(combo.focusPolicy(), QtCore.Qt.FocusPolicy.StrongFocus)

        def wheel() -> None:
            pos = QtCore.QPointF(10, 10)
            event = QtGui.QWheelEvent(
                pos,
                combo.mapToGlobal(pos),
                QtCore.QPoint(),
                QtCore.QPoint(0, -120),
                QtCore.Qt.MouseButton.NoButton,
                QtCore.Qt.KeyboardModifier.NoModifier,
                QtCore.Qt.ScrollPhase.NoScrollPhase,
                False,
            )
            QtWidgets.QApplication.sendEvent(combo, event)
            self.app.processEvents()

        wheel()
        self.assertEqual(combo.currentIndex(), 0)
        self.assertGreater(area.verticalScrollBar().value(), 0)
        combo.setFocus()
        self.app.processEvents()
        wheel()
        self.assertEqual(combo.currentIndex(), 1)
        area.close()
