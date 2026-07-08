"""Matplotlib (Agg) world-map choropleth backend for the ST492 comparison."""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")  # must be set before pyplot import

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import math

import fixture
import geo_data


def render(out_png: str) -> None:
    world = geo_data.load_world()
    counts_norm = geo_data.counts_by_geojson_name(fixture.COUNTRIES)
    maxc = fixture.MAX_COUNT

    # Figure sized to the target pixel geometry.
    fig, ax = plt.subplots(figsize=(fixture.WIDTH / 150.0, fixture.HEIGHT / 150.0))

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
    lons = [p["lon"] for p in fixture.POINTS]
    lats = [p["lat"] for p in fixture.POINTS]
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

    fig.savefig(out_png, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)


if __name__ == "__main__":
    import sys; render(sys.argv[1])
