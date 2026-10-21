"""LIVE interactive Plotly map widget for the BACTolomew "Map Lab".

Builds the same ST492 fixture figure as backend_plotly.py and renders it fully
interactive (pan / wheel-zoom / hover) inside a PyQt6 QWebEngineView by feeding
it the self-contained Plotly HTML.

NOTE: WebEngine may crash on broken-GL hosts (expected; the Lab isolates it).
Under xvfb it renders fine.
"""
from __future__ import annotations

import os
import tempfile

# Use the real GL stack (GPU once the driver is fixed; software under xvfb).
os.environ.setdefault("QTWEBENGINE_CHROMIUM_FLAGS", "--no-sandbox")

import plotly.graph_objects as go
from PyQt6 import QtWidgets

# Keep temp HTML files alive for the lifetime of the process.
_TMP_FILES: list = []

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


def build_figure() -> go.Figure:
    """Build the ST492 choropleth + sample-point figure (same data as backend_plotly).

    The static backend uses Plotly's geo subplot (Choropleth/Scattergeo), whose
    base geometry is fetched as topojson from cdn.plot.ly -- unreachable offline,
    so nothing draws. This LIVE version renders the exact same data on a plain
    *cartesian* plot instead: it draws the country polygons from the local
    Natural Earth geojson (via geo_data) as filled traces in plate-carree
    (linear lon/lat) space. Fully offline, and still pan/wheel-zoom interactive.
    """
    points = [
        {"sample": p["sample"], "lat": p["lat"], "lon": p["lon"], "country": p["country"]}
        for p in fixture.POINTS
    ]
    countries = fixture.COUNTRIES
    max_count = fixture.MAX_COUNT

    counts_norm = geo_data.counts_by_geojson_name(countries)

    # --- one filled polygon trace per country -------------------------------
    traces: list = []
    for name, rings in geo_data.load_world():
        count = geo_data.match_count(name, counts_norm)
        xs: list = []
        ys: list = []
        for ring in rings:
            for lon, lat in ring:
                xs.append(lon)
                ys.append(lat)
            xs.append(None)  # break between rings/subpaths
            ys.append(None)
        traces.append(
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines",
                fill="toself",
                fillcolor=geo_data.count_color(count, max_count),
                line=dict(color=geo_data.BORDER, width=0.4),
                hoverinfo="text",
                text=f"{name}: {count}",
                name=name,
                showlegend=False,
            )
        )

    # --- colorbar proxy (cartesian fills carry no scale of their own) --------
    colorbar_proxy = go.Scatter(
        x=[None],
        y=[None],
        mode="markers",
        marker=dict(
            colorscale=[[0.0, "#deebf7"], [1.0, "#084594"]],
            cmin=0,
            cmax=max_count,
            color=[0],
            showscale=True,
            colorbar=dict(
                title=dict(text="Samples", side="right"),
                thickness=14,
                len=0.7,
                x=1.005,
                tickfont=dict(size=11, color=geo_data.INK),
            ),
        ),
        hoverinfo="skip",
        showlegend=False,
    )
    traces.append(colorbar_proxy)

    # --- sample points ------------------------------------------------------
    lons = [p["lon"] for p in points]
    lats = [p["lat"] for p in points]
    names = [p["sample"] for p in points]
    traces.append(
        go.Scatter(
            x=lons,
            y=lats,
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
    )

    fig = go.Figure(data=traces)
    fig.update_layout(
        title=dict(
            text=fixture.TITLE,
            x=0.5,
            xanchor="center",
            font=dict(size=20, color=geo_data.INK),
        ),
        xaxis=dict(
            range=[-180, 180],
            visible=False,
            constrain="domain",
        ),
        yaxis=dict(
            range=[-60, 85],
            visible=False,
            scaleanchor="x",  # equal degree aspect => plate carree
            scaleratio=1,
            constrain="domain",
        ),
        plot_bgcolor=geo_data.OCEAN,  # ocean
        paper_bgcolor="#ffffff",
        margin=dict(l=8, r=8, t=52, b=8),
        dragmode="pan",
        hovermode="closest",
        autosize=True,
    )
    return fig


def build_widget(parent=None):
    """Return a live, interactive QWidget showing the ST492 fixture map."""
    from PyQt6 import QtCore

    # QtWebEngine requires the shared-OpenGL-contexts attribute. Ideally set
    # before the QApplication is constructed; set it defensively here too so we
    # work even when the app already exists (the Map Lab host sets it early).
    QtWidgets.QApplication.setAttribute(
        QtCore.Qt.ApplicationAttribute.AA_ShareOpenGLContexts, True
    )

    # Import lazily so import-time failures surface here (WebEngine is heavy).
    from PyQt6 import QtWebEngineWidgets

    fig = build_figure()
    # Responsive: fill the WebEngineView; keep pan+zoom interactivity.
    html = fig.to_html(
        include_plotlyjs="inline",
        full_html=True,
        config={
            "scrollZoom": True,
            "displayModeBar": True,
            "displaylogo": False,
            "responsive": True,
        },
    )

    view = QtWebEngineWidgets.QWebEngineView(parent)

    # The inline Plotly bundle is ~4.9 MB, far over QWebEngineView.setHtml()'s
    # ~2 MB data-URL limit, so setHtml() would silently render nothing. Write to
    # a temp file and load it as a real URL instead (still fully live/interactive).
    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".html", delete=False, encoding="utf-8"
    )
    tmp.write(html)
    tmp.flush()
    tmp.close()
    _TMP_FILES.append(tmp.name)

    view.load(QtCore.QUrl.fromLocalFile(tmp.name))

    # Block until the page (and Plotly) has finished loading so a subsequent
    # show()/grab() captures a drawn frame rather than a blank widget.
    loop = QtCore.QEventLoop()
    done = {"ok": False}

    def _finished(ok):
        done["ok"] = ok
        loop.quit()

    view.loadFinished.connect(_finished)
    timeout = QtCore.QTimer()
    timeout.setSingleShot(True)
    timeout.timeout.connect(loop.quit)
    timeout.start(15000)
    loop.exec()

    # WebEngine composites asynchronously in a separate process; a caller that
    # only pumps the event loop a handful of times right after show() would grab
    # a blank frame. Pre-size + show here and spin the loop for a short spell so
    # Plotly has actually drawn at least one frame by the time we return.
    view.resize(1000, 560)
    view.show()
    _pump(1500)

    return view


def _pump(ms: int) -> None:
    """Process events for roughly `ms` milliseconds of wall-clock time."""
    import time

    end = time.monotonic() + ms / 1000.0
    while time.monotonic() < end:
        QtWidgets.QApplication.processEvents()
        time.sleep(0.01)


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    w = build_widget()
    w.resize(1000, 560)
    w.show()
    sys.exit(app.exec())
