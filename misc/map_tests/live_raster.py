"""Live interactive Map Lab widget (raster backend).

Returns the REAL shipping RasterGeoMapView from the annotation_explorer_polars
package, fed with the shared ST492 fixture data. Pan/zoom/highlight are built in.
"""
from __future__ import annotations

import os
import sys

from PyQt6 import QtWidgets  # noqa: F401  (ensures a Qt binding is imported)

import fixture

_SRC = "/home/bals/Documents/dev/BACTolomew/src"
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from annotation_explorer_polars.geo_raster_map import RasterGeoMapView


def build_widget(parent=None):
    view = RasterGeoMapView(parent)

    points = [
        {
            "sample": p["sample"],
            "lat": p["lat"],
            "lon": p["lon"],
            "country": p["country"],
        }
        for p in fixture.POINTS
    ]
    view.update_data(
        points=points,
        countries=fixture.COUNTRIES,
        max_count=fixture.MAX_COUNT,
    )
    return view


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    w = build_widget()
    w.resize(1000, 560)
    w.show()
    sys.exit(app.exec())
