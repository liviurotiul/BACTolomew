"""The cgMLST MST's Settings column (the shared ``settings_column``
look): captions above 30 px drop-downs in a card, folded to its strip while
the panel is too narrow for the graph and opened again when it grows, unless
the user folded it (or opened it in the narrow panel)."""
from __future__ import annotations

import os
import tempfile
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

from PyQt6 import QtWidgets

from annotation_explorer_polars.cgmlst_force_plot import CgmlstForcePlotPanel
from annotation_explorer_polars.theme import apply_application_theme


class CgmlstForceSettingsColumnTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        apply_application_theme(cls.app)

    def setUp(self) -> None:
        self._dir = tempfile.TemporaryDirectory()
        self.panel = CgmlstForcePlotPanel(self._dir.name)
        self.panel.set_trait_values({"Clade": {"a": "1", "b": "2"}, "Country": {"a": "X"}})

    def tearDown(self) -> None:
        self.panel.close()
        self.panel.deleteLater()
        self.app.processEvents()
        self._dir.cleanup()

    def _show_at(self, width: int) -> None:
        self.panel.resize(width, 600)
        self.panel.show()
        for _ in range(5):
            self.app.processEvents()

    def test_controls_are_captioned_drop_downs_in_a_card(self) -> None:
        panel = self.panel
        captions = {
            label.text()
            for label in panel._settings_column.findChildren(QtWidgets.QLabel)
            if label.property("role") == "settingsCaption"
        }
        self.assertEqual(captions, {"Color nodes by", "Shape nodes by"})
        self.assertEqual(panel._color_combo.height(), 30)
        self.assertEqual(
            [panel._shape_combo.itemText(i) for i in range(panel._shape_combo.count())], ["None", "Clade", "Country"]
        )

    def test_a_narrow_panel_folds_the_column_until_it_grows(self) -> None:
        panel = self.panel
        self._show_at(1100)
        self.assertTrue(panel.settings_open())
        self.assertFalse(panel._settings_strip.isVisible())
        self._show_at(500)
        self.assertFalse(panel.settings_open())
        self.assertTrue(panel._settings_strip.isVisible())
        self._show_at(1100)
        self.assertTrue(panel.settings_open())
        # The user's fold stands when the panel grows.
        panel._settings_column.collapse_button.click()
        self._show_at(1300)
        self.assertFalse(panel.settings_open())
        panel._settings_strip.expand_button.click()
        self.assertTrue(panel.settings_open())

    def test_the_users_open_stands_in_a_narrow_panel(self) -> None:
        panel = self.panel
        self._show_at(500)
        self.assertFalse(panel.settings_open())
        panel._settings_strip.expand_button.click()
        self._show_at(490)
        self.assertTrue(panel.settings_open())
        # Once the panel has room for both, the size rule applies again.
        self._show_at(1100)
        self._show_at(500)
        self.assertFalse(panel.settings_open())


if __name__ == "__main__":
    unittest.main()
