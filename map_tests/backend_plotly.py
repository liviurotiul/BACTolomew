"""Plotly (kaleido) sample-map backend for the BACTolomew ST492 comparison.

Renders a choropleth of per-country sample counts plus one marker per sample,
exported to PNG via kaleido (headless Chromium through choreographer).
"""
from __future__ import annotations

import plotly.graph_objects as go

import fixture
import geo_data

# Map fixture display names -> plotly "country names" locationmode labels.
_PLOTLY_NAME = {
    "USA": "United States",
    "UK (Scotland)": "United Kingdom",
    "Kosovo": "Kosovo",
    "Serbia": "Serbia",
    "Romania": "Romania",
    "China": "China",
    "Switzerland": "Switzerland",
    "Montenegro": "Montenegro",
    "Croatia": "Croatia",
    "Hungary": "Hungary",
    "Germany": "Germany",
    "Saudi Arabia": "Saudi Arabia",
    "Slovenia": "Slovenia",
    "France": "France",
    "Lithuania": "Lithuania",
    "Lebanon": "Lebanon",
    "Thailand": "Thailand",
    "Italy": "Italy",
    "Canada": "Canada",
    "Japan": "Japan",
    "Sudan": "Sudan",
    "India": "India",
    "Singapore": "Singapore",
    "Brazil": "Brazil",
    "Fiji": "Fiji",
}


def render(out_png: str) -> None:
    # --- aggregate counts onto plotly country names -------------------------
    agg: dict[str, int] = {}
    for name, count in fixture.COUNTRIES.items():
        plotly_name = _PLOTLY_NAME.get(name, name)
        agg[plotly_name] = agg.get(plotly_name, 0) + int(count)

    locations = list(agg.keys())
    counts = [agg[k] for k in locations]

    # --- choropleth ---------------------------------------------------------
    choropleth = go.Choropleth(
        locations=locations,
        locationmode="country names",
        z=counts,
        zmin=0,
        zmax=fixture.MAX_COUNT,
        colorscale="Blues",
        marker_line_color=geo_data.BORDER,
        marker_line_width=0.4,
        colorbar=dict(
            title=dict(text="Samples", side="right"),
            thickness=14,
            len=0.7,
            x=1.005,
            tickfont=dict(size=11, color=geo_data.INK),
        ),
        hovertemplate="%{location}: %{z}<extra></extra>",
    )

    # --- sample points ------------------------------------------------------
    lons = [p["lon"] for p in fixture.POINTS]
    lats = [p["lat"] for p in fixture.POINTS]
    names = [p["sample"] for p in fixture.POINTS]
    scatter = go.Scattergeo(
        lon=lons,
        lat=lats,
        text=names,
        mode="markers",
        marker=dict(
            color=geo_data.POINT_FILL,
            size=6,
            line=dict(color=geo_data.POINT_EDGE, width=0.8),
            opacity=0.95,
        ),
        hovertemplate="%{text}<extra></extra>",
        showlegend=False,
    )

    fig = go.Figure(data=[choropleth, scatter])
    fig.update_layout(
        title=dict(
            text=fixture.TITLE,
            x=0.5,
            xanchor="center",
            font=dict(size=20, color=geo_data.INK),
        ),
        geo=dict(
            scope="world",
            showcountries=True,
            countrycolor=geo_data.BORDER,
            countrywidth=0.4,
            showland=True,
            landcolor=geo_data.NO_DATA_FILL,
            showocean=True,
            oceancolor=geo_data.OCEAN,
            showcoastlines=True,
            coastlinecolor=geo_data.BORDER,
            coastlinewidth=0.4,
            showframe=False,
            projection_type="equirectangular",
            lataxis=dict(range=[-60, 85]),
            lonaxis=dict(range=[-180, 180]),
            bgcolor=geo_data.OCEAN,
        ),
        width=fixture.WIDTH,
        height=fixture.HEIGHT,
        margin=dict(l=8, r=8, t=52, b=8),
        paper_bgcolor="#ffffff",
    )

    fig.write_image(out_png, scale=2)


if __name__ == "__main__":
    import sys
    render(sys.argv[1])
