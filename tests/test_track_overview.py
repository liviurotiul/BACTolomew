"""The genome track's overview is its scroll bar: dragging the window pans,
a click centers there, dragging an edge zooms with the other edge kept; a
resize keeps the view where it was and the overview in view."""
from __future__ import annotations

import os
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtTest import QTest

from annotation_explorer_polars.igv_view import GeneSegment, IgvTrackView

LEFT = QtCore.Qt.MouseButton.LeftButton
NO_MOD = QtCore.Qt.KeyboardModifier.NoModifier


class TrackOverviewTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def _view(self, height: int) -> IgvTrackView:
        view = IgvTrackView()
        view.resize(700, height)
        view.show()
        genes = [
            GeneSegment(200 + i * 900, 900 + i * 900, "+", f"g{i}", QtGui.QColor("#f08c00"), 1 + i % 2,
                        {"AnnUID": str(i)})
            for i in range(200)
        ]
        view.set_track(genes, 180_000, bp_per_px=10.0, hit_bp_positions=[50_000])
        view.center_on_bp(50_000)
        self.app.processEvents()
        self.addCleanup(view.close)
        return view

    def _drag(self, overview: QtWidgets.QWidget, x_from: float, x_to: float) -> None:
        y = overview.height() // 2
        QTest.mousePress(overview, LEFT, NO_MOD, QtCore.QPoint(round(x_from), y))
        for step in range(1, 6):
            QTest.mouseMove(overview, QtCore.QPoint(round(x_from + (x_to - x_from) * step / 5), y))
        QTest.mouseRelease(overview, LEFT, NO_MOD, QtCore.QPoint(round(x_to), y))

    def test_drag_click_and_edge_zoom(self) -> None:
        view = self._view(300)
        overview = view._overview
        self.assertFalse(view.horizontalScrollBar().isVisible())
        start, end = view.visible_bp_range()
        window = overview.window_rect()
        self._drag(overview, window.center().x(), window.center().x() + 50)
        moved = view.visible_bp_range()
        shift = overview.bp_at(50) - overview.bp_at(0)
        self.assertAlmostEqual(moved[0] - start, shift, delta=view._bp_per_px * 2)
        self.assertAlmostEqual(moved[1] - moved[0], end - start, delta=1)

        axis = overview.axis()
        QTest.mouseClick(overview, LEFT, NO_MOD, QtCore.QPoint(round(axis.left() + axis.width() * 0.75), 10))
        left, right = view.visible_bp_range()
        self.assertAlmostEqual((left + right) / 2, 135_000, delta=overview.bp_at(1) - overview.bp_at(0))

        # Drag the right edge out: zoom out with the left edge kept.
        window = overview.window_rect()
        self._drag(overview, window.right() + 2, window.right() + 60)
        zoomed = view.visible_bp_range()
        self.assertAlmostEqual(zoomed[0], left, delta=view._bp_per_px * 2)
        self.assertGreater(zoomed[1] - zoomed[0], (right - left) * 1.5)

    def test_resize_keeps_position_and_overview(self) -> None:
        # Shrinking across a lane pitch step re-lays the lanes out; the view
        # stays on the same stretch of contig and the overview stays whole.
        view = self._view(300)
        before = view.visible_bp_range()
        for height in (133, 120, 60):
            view.resize(700, height)
            QTest.qWait(250)
            middle = sum(view.visible_bp_range()) / 2
            self.assertAlmostEqual(middle, sum(before) / 2, delta=view._bp_per_px * 2, msg=height)
            self.assertEqual(view._overview.geometry().bottom(), view.contentsRect().bottom())


if __name__ == "__main__":
    unittest.main()
