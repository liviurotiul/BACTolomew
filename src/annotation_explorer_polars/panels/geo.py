"""The Map panel (plan 4.6, D4): the window's one sample set on a map.

Multi-instance and heavy: built on first show, redrawn only while visible.
Two maps of the same window share one payload through the engine's cache.
It shows the countries of hovered samples in the hover color, and hovering a
country with samples hovers that country's samples in the window (origin: this
panel), so the other panels emphasize them; leaving the country or the map
clears it. The samples are those the map counts for the country (the window's
sample set, the map's country normalization).
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

from PyQt6 import QtCore, QtWidgets

from .base import GEO, ExplorePanel
from .pan import hover_kind


class GeoPanel(ExplorePanel):
    KIND = GEO
    TITLE = "Map"
    HEAVY = True

    def __init__(self, window: Any, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setObjectName("GeoPanel")
        self._window = window
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.view: Any = None
        # Normalized country -> its samples in the current sample set.
        self._country_samples: Dict[str, List[str]] = {}
        self._hover_country = ""
        self._sent_hover: Tuple[str, ...] = ()
        try:
            from ..geo_experimental_views import GeoMapLibreView

            self.view = GeoMapLibreView(parent=self)
        except Exception as exc:
            message = QtWidgets.QLabel(f"Failed to initialize Geo map: {exc}", self)
            message.setWordWrap(True)
            message.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
            layout.addWidget(message)
            return
        layout.addWidget(self.view)
        self.view.countryHovered.connect(self._on_country_hovered)

    def apply_result(self, filters: Any, result: Any) -> None:
        del filters
        if self.view is None:
            return
        engine = self._window.engine
        payload = engine.map_payload(result.samples)
        by_country = getattr(engine, "map_country_samples", None)
        self._country_samples = _samples_by_country(by_country(result.samples) if by_country else {})
        self.view.update_data(
            points=payload.get("points", []),
            countries=payload.get("countries", {}),
            max_count=payload.get("maxCountryCount", 1),
            total=payload.get("total"),
        )
        if self._hover_country:  # a hovered country follows the new sample set
            self._send_country_hover()

    def show_hover(self, hover: Any) -> None:
        if self.view is not None:
            self.view.highlight_samples(list(getattr(hover, "samples", ()) or ()))

    # Country hover out (the map is the origin; the frame never echoes it back)
    def country_samples(self, country: str) -> List[str]:
        """The samples the map counts for this normalized country name."""
        return list(self._country_samples.get(str(country or ""), ()))

    def _on_country_hovered(self, country: str) -> None:
        self._hover_country = str(country or "")
        self._send_country_hover()

    def _send_country_hover(self) -> None:
        samples = tuple(self.country_samples(self._hover_country)) if self._hover_country else ()
        # A country without samples sends nothing; it only ends a hover this map started.
        if samples == self._sent_hover:
            return
        self._sent_hover = samples
        if self.ctx is not None:
            self.ctx.hover(samples, origin=self.uid or None, kind=hover_kind(samples))

    def _clear_country_hover(self) -> None:
        # The map backend forgets its country too: it reports only changes, so
        # a stale key would swallow the next hover of the same country.
        if self.view is not None:
            self.view.reset_hover()
        self._on_country_hovered("")

    def leaveEvent(self, event: Any) -> None:  # type: ignore[override]
        super().leaveEvent(event)
        self._clear_country_hover()

    def on_hidden(self) -> None:
        self._clear_country_hover()

    def unbind(self) -> None:
        self._clear_country_hover()
        super().unbind()


def _samples_by_country(country_samples: Dict[str, List[str]]) -> Dict[str, List[str]]:
    """Group the engine's per-country samples by the name the map counts them under."""
    from ..geo_raster_map import _norm

    out: Dict[str, List[str]] = {}
    for name, samples in (country_samples or {}).items():
        key = _norm(name)
        if key:
            out.setdefault(key, []).extend(str(s) for s in samples)
    return out
