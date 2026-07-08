"""Live, interactive SVG map widget for the BACTolomew "Map Lab".

Reuses backend_svg's SVG-string generation and renders it into a
QGraphicsView via a QGraphicsSvgItem, with hand-drag panning and
wheel-zoom. No OpenGL.
"""
from __future__ import annotations

from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtSvg import QSvgRenderer
from PyQt6.QtSvgWidgets import QGraphicsSvgItem

import backend_svg


class _MapView(QtWidgets.QGraphicsView):
    """QGraphicsView with hand-drag pan and cursor-anchored wheel zoom."""

    _MIN_SCALE = 0.2
    _MAX_SCALE = 40.0

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setDragMode(QtWidgets.QGraphicsView.DragMode.ScrollHandDrag)
        self.setTransformationAnchor(
            QtWidgets.QGraphicsView.ViewportAnchor.AnchorUnderMouse
        )
        self.setResizeAnchor(
            QtWidgets.QGraphicsView.ViewportAnchor.AnchorViewCenter
        )
        self.setRenderHints(
            QtGui.QPainter.RenderHint.Antialiasing
            | QtGui.QPainter.RenderHint.SmoothPixmapTransform
        )
        self.setHorizontalScrollBarPolicy(
            QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded
        )
        self.setVerticalScrollBarPolicy(
            QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded
        )
        self.setBackgroundBrush(QtGui.QColor(geo_ocean()))
        self._zoom = 1.0

    def wheelEvent(self, event):
        delta = event.angleDelta().y()
        if delta == 0:
            return
        factor = 1.25 if delta > 0 else 1.0 / 1.25
        new_zoom = self._zoom * factor
        if new_zoom < self._MIN_SCALE:
            factor = self._MIN_SCALE / self._zoom
            new_zoom = self._MIN_SCALE
        elif new_zoom > self._MAX_SCALE:
            factor = self._MAX_SCALE / self._zoom
            new_zoom = self._MAX_SCALE
        self._zoom = new_zoom
        self.scale(factor, factor)
        event.accept()


def geo_ocean() -> str:
    try:
        import geo_data
        return geo_data.OCEAN
    except Exception:
        return "#f4f8fc"


def build_widget(parent=None) -> QtWidgets.QWidget:
    """Return a live, interactive QWidget showing the ST492 fixture map."""
    # Reuse backend_svg's SVG generation (same data as the static backend).
    svg = backend_svg._build_svg()
    svg_bytes = QtCore.QByteArray(svg.encode("utf-8"))

    renderer = QSvgRenderer(svg_bytes)
    item = QGraphicsSvgItem()
    item.setSharedRenderer(renderer)
    # keep the renderer alive for the item's lifetime
    item._renderer = renderer

    scene = QtWidgets.QGraphicsScene(parent)
    scene.setBackgroundBrush(QtGui.QColor(geo_ocean()))
    scene.addItem(item)
    scene.setSceneRect(item.boundingRect())

    view = _MapView(parent)
    view.setScene(scene)
    # keep references alive on the view
    view._scene = scene
    view._item = item

    def _fit():
        view.fitInView(item.boundingRect(), QtCore.Qt.AspectRatioMode.KeepAspectRatio)
        view._zoom = 1.0

    # fit on first show
    QtCore.QTimer.singleShot(0, _fit)
    view._fit = _fit
    return view


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    w = build_widget()
    w.resize(1000, 560)
    w.show()
    sys.exit(app.exec())
