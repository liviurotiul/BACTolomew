"""Live, interactive Cartopy map widget for the BACTolomew Map Lab.

Reuses backend_cartopy's rendering logic (offline geo_data polygons on a
PlateCarree axes) but draws onto a live FigureCanvasQTAgg embedded in a
QWidget, with a NavigationToolbar2QT for pan + wheel-zoom. Fully offline:
no Natural Earth downloads, no GL, no network.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("QtAgg")

from matplotlib.figure import Figure
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
import cartopy.crs as ccrs
from shapely.geometry import Polygon

from PyQt6 import QtWidgets

import fixture
import geo_data


def _draw_map(fig, points, countries, max_count):
    """Render the ST492 choropleth + sample points onto a cartopy axes.

    Mirrors backend_cartopy.render(), but paints onto a live Figure/axes
    instead of saving a PNG.
    """
    counts_norm = geo_data.counts_by_geojson_name(countries)
    maxc = max_count

    proj = ccrs.PlateCarree()
    ax = fig.add_subplot(1, 1, 1, projection=proj)
    ax.set_global()

    # Ocean-tinted background.
    ax.set_facecolor(geo_data.OCEAN)
    try:
        ax.background_patch.set_facecolor(geo_data.OCEAN)
    except AttributeError:
        pass

    # Optional graticule (pure matplotlib, no network).
    try:
        ax.gridlines(draw_labels=False, linewidth=0.3,
                     color="#d7e0ea", alpha=0.7)
    except Exception:
        pass

    # Choropleth: build shapely polygons from the offline rings and add them.
    for name, rings in geo_data.load_world():
        count = geo_data.match_count(name, counts_norm)
        fill = geo_data.count_color(count, maxc)
        geoms = []
        for ring in rings:
            if len(ring) >= 3:
                geoms.append(Polygon(ring))
        if not geoms:
            continue
        try:
            ax.add_geometries(
                geoms, crs=proj,
                facecolor=fill, edgecolor=geo_data.BORDER,
                linewidth=0.3, zorder=2,
            )
        except Exception:
            pass

    # Sample points.
    lons = [p["lon"] for p in points]
    lats = [p["lat"] for p in points]
    ax.scatter(
        lons, lats, transform=proj,
        s=26, facecolor=geo_data.POINT_FILL,
        edgecolor=geo_data.POINT_EDGE, linewidth=0.6,
        alpha=0.95, zorder=5,
    )

    ax.set_title(fixture.TITLE, fontsize=15, color=geo_data.INK,
                 fontweight="bold", pad=10)

    # Legend / colour scale for the counts (log ramp, matching count_color).
    stops = [0, 1, 2, 5, 10, maxc]
    seen = set()
    handles = []
    labels = []
    for c in stops:
        if c in seen:
            continue
        seen.add(c)
        col = geo_data.count_color(c, maxc)
        handles.append(Patch(facecolor=col, edgecolor=geo_data.BORDER,
                             linewidth=0.3))
        labels.append("no data" if c == 0 else f"{c}")
    handles.append(Line2D([0], [0], marker="o", linestyle="none",
                          markerfacecolor=geo_data.POINT_FILL,
                          markeredgecolor=geo_data.POINT_EDGE,
                          markersize=7))
    labels.append("sample")

    leg = ax.legend(
        handles, labels, title="samples / country",
        loc="lower left", fontsize=8, title_fontsize=9,
        frameon=True, framealpha=0.9, ncol=1,
        borderpad=0.6, handlelength=1.4, labelspacing=0.35,
    )
    leg.get_frame().set_edgecolor(geo_data.BORDER)
    leg.get_frame().set_facecolor("#ffffff")
    for txt in leg.get_texts():
        txt.set_color(geo_data.INK)
    leg.get_title().set_color(geo_data.INK)
    leg.set_zorder(10)

    return ax


class _WheelZoomFilter:
    """Attaches scroll-to-zoom to a cartopy axes (toolbar gives pan)."""

    def __init__(self, ax, canvas):
        self.ax = ax
        self.canvas = canvas
        self.cid = canvas.mpl_connect("scroll_event", self._on_scroll)

    def _on_scroll(self, event):
        ax = self.ax
        if event.inaxes is not ax:
            return
        base = 1.2
        scale = 1.0 / base if event.button == "up" else base
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        xd = event.xdata if event.xdata is not None else (x0 + x1) / 2.0
        yd = event.ydata if event.ydata is not None else (y0 + y1) / 2.0
        nx0 = xd - (xd - x0) * scale
        nx1 = xd + (x1 - xd) * scale
        ny0 = yd - (yd - y0) * scale
        ny1 = yd + (y1 - yd) * scale
        ax.set_xlim(nx0, nx1)
        ax.set_ylim(ny0, ny1)
        self.canvas.draw_idle()


def build_widget(parent=None):
    """Build a live, interactive ST492 cartopy map widget.

    Returns a QWidget embedding a matplotlib FigureCanvas (cartopy PlateCarree
    axes) plus a NavigationToolbar. Pan via the toolbar's pan tool; zoom via
    the toolbar's zoom-rect tool or the mouse wheel.
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

    fig = Figure(figsize=(fixture.WIDTH / 100.0, fixture.HEIGHT / 100.0),
                 facecolor=geo_data.OCEAN)
    canvas = FigureCanvasQTAgg(fig)
    ax = _draw_map(fig, points, countries, max_count)

    toolbar = NavigationToolbar2QT(canvas, container)

    layout.addWidget(toolbar)
    layout.addWidget(canvas, 1)

    # Keep references alive on the widget.
    container._fig = fig
    container._canvas = canvas
    container._ax = ax
    container._toolbar = toolbar
    container._wheel = _WheelZoomFilter(ax, canvas)

    canvas.draw_idle()
    return container


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    w = build_widget()
    w.resize(1000, 560)
    w.show()
    sys.exit(app.exec())
