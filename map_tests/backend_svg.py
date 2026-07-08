"""SVG map backend for the BACTolomew ST492 sample-map comparison.

Builds an SVG string (viewBox "0 0 1100 620") with an ocean rect, a choropleth
path per country ring, sample-point circles, a title and a colour-scale legend.
The SVG is written to <out>.svg and then rasterised to PNG with PyQt6 QtSvg.
"""
from __future__ import annotations

import fixture
import geo_data

W = fixture.WIDTH   # 1100
H = fixture.HEIGHT  # 620


def _fnum(v: float) -> str:
    """Round to 1 decimal and drop a trailing '.0' to keep the SVG small."""
    s = f"{v:.1f}"
    if s.endswith(".0"):
        s = s[:-2]
    return s


def _ring_path(ring) -> str:
    """One 'M x y L x y ... Z' subpath for a lon/lat ring (skip degenerate)."""
    if len(ring) < 3:
        return ""
    pts = []
    for lon, lat in ring:
        x, y = geo_data.project(lon, lat, W, H)
        pts.append((x, y))
    d = "M" + f"{_fnum(pts[0][0])} {_fnum(pts[0][1])}"
    for x, y in pts[1:]:
        d += "L" + f"{_fnum(x)} {_fnum(y)}"
    return d + "Z"


def _build_svg() -> str:
    counts_norm = geo_data.counts_by_geojson_name(fixture.COUNTRIES)
    world = geo_data.load_world()

    parts = []
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{W}" height="{H}" viewBox="0 0 {W} {H}">'
    )

    # 1. ocean background
    parts.append(
        f'<rect x="0" y="0" width="{W}" height="{H}" fill="{geo_data.OCEAN}"/>'
    )

    # 2+3. choropleth countries with thin borders
    parts.append(
        f'<g stroke="{geo_data.BORDER}" stroke-width="0.4" '
        f'stroke-linejoin="round">'
    )
    for name, rings in world:
        count = geo_data.match_count(name, counts_norm)
        fill = geo_data.count_color(count, fixture.MAX_COUNT)
        d = "".join(_ring_path(r) for r in rings)
        if not d:
            continue
        parts.append(f'<path d="{d}" fill="{fill}" fill-rule="evenodd"/>')
    parts.append("</g>")

    # 4. sample points
    parts.append(
        f'<g fill="{geo_data.POINT_FILL}" stroke="{geo_data.POINT_EDGE}" '
        f'stroke-width="0.7">'
    )
    for p in fixture.POINTS:
        x, y = geo_data.project(p["lon"], p["lat"], W, H)
        parts.append(
            f'<circle cx="{_fnum(x)}" cy="{_fnum(y)}" r="2.6"/>'
        )
    parts.append("</g>")

    # 5. title
    parts.append(
        f'<text x="{W/2:.0f}" y="34" text-anchor="middle" '
        f'font-family="DejaVu Sans, Arial, sans-serif" font-size="22" '
        f'font-weight="700" fill="{geo_data.INK}">{_escape(fixture.TITLE)}</text>'
    )

    # 6. legend: colour scale for the counts
    parts.append(_legend())

    parts.append("</svg>")
    return "".join(parts)


def _escape(s: str) -> str:
    return (
        s.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )


def _legend() -> str:
    """Discrete colour swatches + labels, lower-left panel."""
    maxc = fixture.MAX_COUNT
    # sample counts across the log ramp
    stops = [0, 1, 3, 8, 15, maxc]
    seen = set()
    ticks = []
    for c in stops:
        if c not in seen and c <= maxc:
            seen.add(c)
            ticks.append(c)
    if maxc not in seen:
        ticks.append(maxc)

    sw = 26            # swatch width
    sh = 14            # swatch height
    pad = 12
    x0 = 20
    panel_w = pad * 2 + sw * len(ticks)
    panel_h = 62
    y0 = H - panel_h - 20

    g = [f'<g font-family="DejaVu Sans, Arial, sans-serif">']
    # panel background
    g.append(
        f'<rect x="{x0}" y="{y0}" width="{panel_w}" height="{panel_h}" '
        f'rx="6" fill="#ffffff" fill-opacity="0.88" '
        f'stroke="{geo_data.BORDER}" stroke-width="1"/>'
    )
    g.append(
        f'<text x="{x0 + pad}" y="{y0 + 18}" font-size="12" '
        f'font-weight="700" fill="{geo_data.INK}">Samples per country</text>'
    )
    sy = y0 + 26
    for i, c in enumerate(ticks):
        sx = x0 + pad + i * sw
        fill = geo_data.count_color(c, maxc)
        g.append(
            f'<rect x="{sx}" y="{sy}" width="{sw}" height="{sh}" '
            f'fill="{fill}" stroke="{geo_data.BORDER}" stroke-width="0.5"/>'
        )
        g.append(
            f'<text x="{sx + sw/2:.0f}" y="{sy + sh + 12}" '
            f'text-anchor="middle" font-size="10" fill="{geo_data.INK}">'
            f'{c}</text>'
        )
    g.append("</g>")
    return "".join(g)


def render(out_png: str) -> None:
    svg = _build_svg()

    # write raw SVG alongside the PNG
    if out_png.lower().endswith(".png"):
        out_svg = out_png[:-4] + ".svg"
    else:
        out_svg = out_png + ".svg"
    with open(out_svg, "w", encoding="utf-8") as fh:
        fh.write(svg)

    # rasterise SVG -> PNG via PyQt6 QtSvg
    from PyQt6 import QtWidgets, QtGui, QtCore
    from PyQt6.QtSvg import QSvgRenderer

    app = QtWidgets.QApplication.instance()
    if app is None:
        app = QtWidgets.QApplication([])

    scale = 2  # higher-DPI raster
    img = QtGui.QImage(
        W * scale, H * scale, QtGui.QImage.Format.Format_ARGB32_Premultiplied
    )
    img.fill(QtGui.QColor("white"))

    renderer = QSvgRenderer(QtCore.QByteArray(svg.encode("utf-8")))
    painter = QtGui.QPainter(img)
    painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
    renderer.render(painter, QtCore.QRectF(0, 0, W * scale, H * scale))
    painter.end()

    img.save(out_png, "PNG")


if __name__ == "__main__":
    import sys
    render(sys.argv[1])
