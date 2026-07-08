"""Live interactive matplotlib (Agg) map widget for the BACTolomew Map Lab.

Reuses backend_matplotlib's figure-drawing logic but embeds a live Figure via
FigureCanvasQTAgg plus a NavigationToolbar2QT so the user gets pan / zoom /
home / save. Same ST492 fixture data as backend_matplotlib.py.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")  # Agg backend, no GL

import math

from matplotlib.figure import Figure
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg,
    NavigationToolbar2QT,
)

from PyQt6 import QtWidgets

import fixture
import geo_data


def _draw_map(fig, points, countries, max_count):
    """Draw the ST492 choropleth + sample points onto a matplotlib Figure.

    This mirrors backend_matplotlib.render() but targets a supplied Figure so
    the same drawing can be embedded live in a Qt canvas.
    """
    counts_norm = geo_data.counts_by_geojson_name(countries)
    maxc = max_count
    world = geo_data.load_world()

    ax = fig.add_subplot(111)

    # Ocean-tinted background.
    ax.set_facecolor(geo_data.OCEAN)
    fig.patch.set_facecolor("#ffffff")

    # Draw every country as a filled polygon coloured by its sample count.
    patches = []
    facecolors = []
    for name, rings in world:
        count = geo_data.match_count(name, counts_norm)
        fill = geo_data.count_color(count, maxc)
        for ring in rings:
            if len(ring) < 3:
                continue
            patches.append(MplPolygon(ring, closed=True))
            facecolors.append(fill)

    coll = PatchCollection(patches, match_original=False)
    coll.set_facecolor(facecolors)
    coll.set_edgecolor(geo_data.BORDER)
    coll.set_linewidth(0.35)
    ax.add_collection(coll)

    # Sample points (lon, lat).
    lons = [p["lon"] for p in points]
    lats = [p["lat"] for p in points]
    ax.scatter(
        lons, lats,
        s=18,
        c=geo_data.POINT_FILL,
        edgecolors=geo_data.POINT_EDGE,
        linewidths=0.6,
        zorder=5,
        alpha=0.95,
    )

    # Equirectangular framing.
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_aspect("equal")

    # Clean look: no spines, no ticks.
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])

    # Title.
    ax.set_title(fixture.TITLE, fontsize=14, fontweight="bold",
                 color=geo_data.INK, pad=12)

    # Colour scale legend using the same log-ramp as count_color.
    ramp = LinearSegmentedColormap.from_list(
        "counts",
        [tuple(c / 255.0 for c in geo_data._RAMP_LO),
         tuple(c / 255.0 for c in geo_data._RAMP_HI)],
    )

    class _LogCountNorm(Normalize):
        def __call__(self, value, clip=None):
            import numpy as np
            v = np.asarray(value, dtype=float)
            denom = math.log1p(maxc) if maxc > 0 else 1.0
            return np.clip(np.log1p(np.clip(v, 0, None)) / denom, 0.0, 1.0)

    sm = ScalarMappable(cmap=ramp, norm=_LogCountNorm(vmin=0, vmax=maxc))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.02,
                        aspect=30, shrink=0.6)
    ticks = [0, 1, 5, 10, maxc]
    ticks = sorted({t for t in ticks if 0 <= t <= maxc})
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([str(t) for t in ticks])
    cbar.ax.tick_params(labelsize=8, color=geo_data.BORDER, labelcolor=geo_data.INK)
    cbar.outline.set_edgecolor(geo_data.BORDER)
    cbar.outline.set_linewidth(0.5)
    cbar.set_label("samples per country", fontsize=8, color=geo_data.INK)

    return ax


def build_widget(parent=None):
    """Return a live, interactive Qt widget showing the ST492 fixture map.

    Pan / wheel-zoom / home / save are provided by the matplotlib
    NavigationToolbar2QT. The map itself is drawn with the same logic as
    backend_matplotlib.render().
    """
    points = [
        {"sample": p["sample"], "lat": p["lat"], "lon": p["lon"],
         "country": p["country"]}
        for p in fixture.POINTS
    ]
    countries = fixture.COUNTRIES
    max_count = fixture.MAX_COUNT

    container = QtWidgets.QWidget(parent)
    layout = QtWidgets.QVBoxLayout(container)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(0)

    fig = Figure(figsize=(fixture.WIDTH / 150.0, fixture.HEIGHT / 150.0))
    canvas = FigureCanvasQTAgg(fig)

    _draw_map(fig, points, countries, max_count)

    toolbar = NavigationToolbar2QT(canvas, container)

    # Wheel-zoom on top of the toolbar's pan/zoom-rect tools.
    def _on_scroll(event):
        ax = event.inaxes
        if ax is None or event.xdata is None or event.ydata is None:
            return
        scale = 0.8 if event.button == "up" else 1.25
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        xd, yd = event.xdata, event.ydata
        ax.set_xlim(xd - (xd - x0) * scale, xd + (x1 - xd) * scale)
        ax.set_ylim(yd - (yd - y0) * scale, yd + (y1 - yd) * scale)
        canvas.draw_idle()

    canvas.mpl_connect("scroll_event", _on_scroll)

    layout.addWidget(toolbar)
    layout.addWidget(canvas, 1)

    return container


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    w = build_widget()
    w.resize(1000, 560)
    w.show()
    sys.exit(app.exec())
