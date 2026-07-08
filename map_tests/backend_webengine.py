"""QtWebEngine map backend.

Renders the ST492 world choropleth as a self-contained (no-network) HTML page
whose body is a full-page inline <svg> built the SAME way as the svg backend
(ocean rect + country <path> choropleth + point <circle>s + title + legend).
The page is loaded into an off-screen QWebEngineView and grabbed to a PNG.

This backend is EXPECTED to possibly crash on broken-GL machines; that failure
is exactly the signal the harness wants. GPU is disabled and software OpenGL is
forced up front so it at least renders under xvfb.
"""
from __future__ import annotations

# MUST be set before importing PyQt6 / WebEngine.
import os
os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = (
    "--disable-gpu --disable-software-rasterizer --no-sandbox --in-process-gpu"
)
os.environ.setdefault("QT_OPENGL", "software")

import html as _html

from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWebEngineWidgets import QWebEngineView

import fixture
import geo_data


def _build_svg() -> str:
    """Build the full-page SVG markup (same construction as the svg backend)."""
    W, H = fixture.WIDTH, fixture.HEIGHT
    counts_norm = geo_data.counts_by_geojson_name(fixture.COUNTRIES)
    world = geo_data.load_world()

    parts: list[str] = []
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
        f'width="{W}" height="{H}">'
    )

    # 1. ocean-tinted background
    parts.append(
        f'<rect x="0" y="0" width="{W}" height="{H}" fill="{geo_data.OCEAN}"/>'
    )

    # 2 + 3. choropleth countries with thin borders
    for name, rings in world:
        count = geo_data.match_count(name, counts_norm)
        fill = geo_data.count_color(count, fixture.MAX_COUNT)
        subpaths: list[str] = []
        for ring in rings:
            if len(ring) < 3:
                continue
            coords = []
            for lon, lat in ring:
                x, y = geo_data.project(lon, lat, W, H)
                coords.append(f"{x:.2f},{y:.2f}")
            if coords:
                subpaths.append("M" + "L".join(coords) + "Z")
        if not subpaths:
            continue
        d = "".join(subpaths)
        parts.append(
            f'<path d="{d}" fill="{fill}" stroke="{geo_data.BORDER}" '
            f'stroke-width="0.4" fill-rule="evenodd"/>'
        )

    # 4. sample points
    for p in fixture.POINTS:
        x, y = geo_data.project(p["lon"], p["lat"], W, H)
        parts.append(
            f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.1" '
            f'fill="{geo_data.POINT_FILL}" stroke="{geo_data.POINT_EDGE}" '
            f'stroke-width="0.9"/>'
        )

    # 5. title
    title = _html.escape(fixture.TITLE)
    parts.append(
        f'<rect x="0" y="0" width="{W}" height="34" fill="#ffffff" '
        f'fill-opacity="0.72"/>'
    )
    parts.append(
        f'<text x="18" y="23" font-family="Segoe UI, Helvetica, Arial, sans-serif" '
        f'font-size="17" font-weight="600" fill="{geo_data.INK}">{title}</text>'
    )

    # 6. legend / colour scale
    parts.append(_legend_svg(W, H))

    parts.append("</svg>")
    return "".join(parts)


def _legend_svg(W: int, H: int) -> str:
    """A compact horizontal count colour-scale in the lower-left."""
    maxc = fixture.MAX_COUNT
    lx, ly = 18, H - 58
    lw, lh = 220, 12
    steps = 60
    parts: list[str] = []

    # panel background
    parts.append(
        f'<rect x="{lx - 10}" y="{ly - 26}" width="{lw + 74}" height="60" rx="7" '
        f'fill="#ffffff" fill-opacity="0.82" stroke="{geo_data.BORDER}" '
        f'stroke-width="1"/>'
    )
    parts.append(
        f'<text x="{lx}" y="{ly - 10}" font-family="Segoe UI, Helvetica, Arial, '
        f'sans-serif" font-size="12" font-weight="600" fill="{geo_data.INK}">'
        f'Samples per country</text>'
    )

    # gradient bar via discrete rects following the same log ramp
    seg_w = lw / steps
    for i in range(steps):
        t = i / (steps - 1)
        # invert log1p ramp so colour matches count_color across the scale
        import math
        count = math.expm1(t * math.log1p(maxc))
        col = geo_data.count_color(max(1, round(count)), maxc)
        parts.append(
            f'<rect x="{lx + i * seg_w:.2f}" y="{ly}" width="{seg_w + 0.6:.2f}" '
            f'height="{lh}" fill="{col}"/>'
        )
    parts.append(
        f'<rect x="{lx}" y="{ly}" width="{lw}" height="{lh}" fill="none" '
        f'stroke="{geo_data.BORDER}" stroke-width="0.8"/>'
    )

    # tick labels: 1 and max
    parts.append(
        f'<text x="{lx}" y="{ly + lh + 13}" font-family="Segoe UI, Helvetica, '
        f'Arial, sans-serif" font-size="11" fill="{geo_data.INK}">1</text>'
    )
    parts.append(
        f'<text x="{lx + lw}" y="{ly + lh + 13}" text-anchor="end" '
        f'font-family="Segoe UI, Helvetica, Arial, sans-serif" font-size="11" '
        f'fill="{geo_data.INK}">{maxc}</text>'
    )

    # no-data swatch
    sx = lx + lw + 12
    parts.append(
        f'<rect x="{sx}" y="{ly}" width="{lh}" height="{lh}" '
        f'fill="{geo_data.NO_DATA_FILL}" stroke="{geo_data.BORDER}" '
        f'stroke-width="0.8"/>'
    )
    parts.append(
        f'<text x="{sx + lh + 5}" y="{ly + lh - 1}" font-family="Segoe UI, '
        f'Helvetica, Arial, sans-serif" font-size="11" fill="{geo_data.INK}">'
        f'no data</text>'
    )
    return "".join(parts)


def _build_html() -> str:
    svg = _build_svg()
    return (
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<style>html,body{margin:0;padding:0;background:"
        f"{geo_data.OCEAN};" "}"
        "svg{display:block}</style></head><body>"
        f"{svg}</body></html>"
    )


def render(out_png: str) -> None:
    W, H = fixture.WIDTH, fixture.HEIGHT

    app = QtWidgets.QApplication.instance()
    if app is None:
        # non-empty argv is REQUIRED by WebEngine
        app = QtWidgets.QApplication(["x"])

    view = QWebEngineView()
    view.setAttribute(QtCore.Qt.WidgetAttribute.WA_DontShowOnScreen, True)
    view.resize(W, H)
    view.show()

    state = {"loaded": False, "ok": False}
    loop = QtCore.QEventLoop()

    def _on_load(ok: bool) -> None:
        state["loaded"] = True
        state["ok"] = bool(ok)
        if loop.isRunning():
            loop.quit()

    view.loadFinished.connect(_on_load)

    # safety timeout (~8s)
    timer = QtCore.QTimer()
    timer.setSingleShot(True)
    timer.timeout.connect(lambda: loop.quit() if loop.isRunning() else None)
    timer.start(8000)

    view.setHtml(_build_html())

    if not state["loaded"]:
        loop.exec()
    timer.stop()

    # let the compositor paint the freshly-loaded content
    for _ in range(6):
        app.processEvents()
        t = QtCore.QElapsedTimer()
        t.start()
        while t.elapsed() < 120:
            app.processEvents()

    image = QtGui.QImage(W, H, QtGui.QImage.Format.Format_ARGB32)
    image.fill(QtGui.QColor(geo_data.OCEAN))
    painter = QtGui.QPainter(image)
    view.render(painter)
    painter.end()

    # grab() as a fallback / primary; prefer the grab if it is non-empty
    pix = view.grab()
    if not pix.isNull():
        image = pix.toImage()

    if not image.save(out_png, "PNG"):
        raise RuntimeError(f"QImage.save failed for {out_png}")


if __name__ == "__main__":
    import sys
    render(sys.argv[1])
