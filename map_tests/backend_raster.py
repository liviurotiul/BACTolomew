"""Pure Qt raster map backend: QImage + QPainter, no widget, no GL."""
from __future__ import annotations

from PyQt6 import QtGui, QtCore, QtWidgets

import fixture
import geo_data

SCALE = 2
W = fixture.WIDTH * SCALE
H = fixture.HEIGHT * SCALE


def _qcolor(hexstr: str) -> QtGui.QColor:
    return QtGui.QColor(hexstr)


def render(out_png: str) -> None:
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    img = QtGui.QImage(W, H, QtGui.QImage.Format.Format_ARGB32)
    img.fill(_qcolor(geo_data.OCEAN))

    p = QtGui.QPainter(img)
    p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
    p.setRenderHint(QtGui.QPainter.RenderHint.TextAntialiasing, True)

    counts = geo_data.counts_by_geojson_name(fixture.COUNTRIES)

    border_pen = QtGui.QPen(_qcolor(geo_data.BORDER))
    border_pen.setWidthF(0.6 * SCALE)
    border_pen.setJoinStyle(QtCore.Qt.PenJoinStyle.RoundJoin)

    # Countries as choropleth
    for name, rings in geo_data.load_world():
        path = QtGui.QPainterPath()
        path.setFillRule(QtCore.Qt.FillRule.OddEvenFill)
        for ring in rings:
            if len(ring) < 2:
                continue
            x0, y0 = geo_data.project(ring[0][0], ring[0][1], W, H)
            path.moveTo(x0, y0)
            for lon, lat in ring[1:]:
                x, y = geo_data.project(lon, lat, W, H)
                path.lineTo(x, y)
            path.closeSubpath()
        c = geo_data.match_count(name, counts)
        fill = _qcolor(geo_data.count_color(c, fixture.MAX_COUNT))
        p.setPen(border_pen)
        p.setBrush(fill)
        p.drawPath(path)

    # Sample points
    r = 3.2 * SCALE
    pt_pen = QtGui.QPen(_qcolor(geo_data.POINT_EDGE))
    pt_pen.setWidthF(0.9 * SCALE)
    p.setPen(pt_pen)
    p.setBrush(_qcolor(geo_data.POINT_FILL))
    for pt in fixture.POINTS:
        x, y = geo_data.project(pt["lon"], pt["lat"], W, H)
        p.drawEllipse(QtCore.QPointF(x, y), r, r)

    # Title
    ink = _qcolor(geo_data.INK)
    title_font = QtGui.QFont("Sans Serif", 15 * SCALE)
    title_font.setBold(True)
    p.setFont(title_font)
    p.setPen(ink)
    margin = 14 * SCALE
    fm = QtGui.QFontMetrics(title_font)
    p.drawText(margin, margin + fm.ascent(), fixture.TITLE)

    # Legend: horizontal colour scale
    _draw_legend(p, ink)

    p.end()
    img.save(out_png)


def _draw_legend(p: QtGui.QPainter, ink: QtGui.QColor) -> None:
    maxc = fixture.MAX_COUNT
    # representative counts low..high
    stops = [1, max(2, maxc // 8), max(4, maxc // 4), max(8, maxc // 2), maxc]
    seen = []
    for s in stops:
        if s not in seen and s >= 1:
            seen.append(s)
    stops = seen

    sw = 34 * SCALE          # swatch width
    sh = 16 * SCALE          # swatch height
    gap = 2 * SCALE
    n = len(stops)
    total_w = n * sw + (n - 1) * gap
    x0 = W - total_w - 18 * SCALE
    y0 = H - 46 * SCALE

    label_font = QtGui.QFont("Sans Serif", 8 * SCALE)
    label_font.setBold(True)
    p.setFont(label_font)
    fm = QtGui.QFontMetrics(label_font)

    # header
    p.setPen(ink)
    p.drawText(x0, y0 - 6 * SCALE, "samples per country")

    # panel background for readability
    pad = 6 * SCALE
    panel = QtCore.QRectF(
        x0 - pad,
        y0 - 6 * SCALE - fm.height(),
        total_w + 2 * pad,
        sh + fm.height() + 6 * SCALE + 2 * pad,
    )
    bg = QtGui.QColor(255, 255, 255, 205)
    p.setPen(QtCore.Qt.PenStyle.NoPen)
    p.setBrush(bg)
    p.drawRoundedRect(panel, 5 * SCALE, 5 * SCALE)

    # re-draw header on top of panel
    p.setPen(ink)
    p.setFont(label_font)
    p.drawText(x0, y0 - 6 * SCALE, "samples per country")

    swatch_pen = QtGui.QPen(_qcolor(geo_data.BORDER))
    swatch_pen.setWidthF(0.6 * SCALE)
    for i, s in enumerate(stops):
        x = x0 + i * (sw + gap)
        rect = QtCore.QRectF(x, y0, sw, sh)
        p.setPen(swatch_pen)
        p.setBrush(_qcolor(geo_data.count_color(s, maxc)))
        p.drawRect(rect)
        p.setPen(ink)
        label = str(s)
        tw = fm.horizontalAdvance(label)
        p.drawText(int(x + (sw - tw) / 2), int(y0 + sh + fm.ascent() + 2 * SCALE), label)


if __name__ == "__main__":
    import sys
    render(sys.argv[1])
