"""Cartopy (matplotlib Agg) backend for the ST492 sample map.

Draws offline polygons from geo_data.load_world() as a choropleth, plus sample
points. Does NOT rely on Natural Earth downloads -- any cartopy feature-download
call is wrapped in try/except so it still renders fully offline.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
from shapely.geometry import Polygon

import fixture
import geo_data


def render(out_png: str) -> None:
    counts_norm = geo_data.counts_by_geojson_name(fixture.COUNTRIES)
    maxc = fixture.MAX_COUNT

    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(fixture.WIDTH / 100.0, fixture.HEIGHT / 100.0))
    ax = plt.axes(projection=proj)
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
    lons = [p["lon"] for p in fixture.POINTS]
    lats = [p["lat"] for p in fixture.POINTS]
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

    fig.savefig(out_png, dpi=150, bbox_inches="tight",
                facecolor=geo_data.OCEAN)
    plt.close(fig)


if __name__ == "__main__":
    import sys; render(sys.argv[1])
