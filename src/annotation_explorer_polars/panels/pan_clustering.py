"""The Pangenome clustering panel: how a part of the pangenome (a
prevalence range picked in the band under the plot) makes the window's
isolates cluster, as a classical MDS of their Jaccard distances in 2D or 3D
with k-means clusters (``pan_clustering_view.PanClusteringView``).

Heavy and multi-instance, like the pan panels. Hover goes both ways through
the window's bus; a k-means cluster can narrow the window's sample set or
open in a new explore window; a double-click opens the sample card.
"""
from __future__ import annotations

from typing import Any, Dict, Optional

from PyQt6 import QtWidgets

from ..prevalence_band import percent_text
from ..session.requests import OpenSampleCard
from .base import PAN_MDS, ExplorePanel
from .pan import hover_kind, hovered_samples, request_pan_subset, request_pan_subset_window


class PanClusteringPanel(ExplorePanel):
    KIND = PAN_MDS
    TITLE = "Pangenome clustering"
    GROUP = "Pangenome"
    HEAVY = True

    def __init__(self, window: Any, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setObjectName("PanClusteringPanel")
        self._window = window
        self.view: Any = None
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

    def set_view(self, view: Any) -> None:
        """The PanClusteringView (built by the window's factory as this panel's child)."""
        self.view = view
        self.layout().addWidget(view)
        view.samplesHovered.connect(self._on_samples_hovered)
        view.sampleOpenRequested.connect(self._on_sample_open_requested)
        view.subsetRequested.connect(self._on_subset_requested)
        view.band.rangeCommitted.connect(self.titleDetailChanged)

    def title_detail(self) -> str:
        """The prevalence range ("15 to 95%")."""
        if self.view is None:
            return ""
        return self.state_detail({"low": self.view.band.low(), "high": self.view.band.high()})

    @staticmethod
    def state_detail(state: Dict[str, Any]) -> str:
        """``title_detail`` from a saved state (PanelSpec.state_detail); ""
        without a range (the view's default is not known before the build)."""
        try:
            low, high = percent_text(float(state["low"])), percent_text(float(state["high"]))
        except (KeyError, TypeError, ValueError):
            return ""
        return f"{low[:-1]} to {high}"

    def apply_result(self, filters: Any, result: Any) -> None:
        del filters
        self.view.set_window_samples(list(result.samples))

    def show_hover(self, hover: Any) -> None:
        self.view.highlight_samples(tuple(getattr(hover, "samples", ()) or ()))

    def save_state(self) -> Dict[str, Any]:
        return self.view.save_state() if self.view is not None else {}

    def restore_state(self, state: Dict[str, Any]) -> None:
        if self.view is not None and state:
            self.view.restore_state(dict(state))
            self.titleDetailChanged.emit()

    def dispose(self) -> None:
        if self.view is not None:
            self.view.dispose()

    def _on_samples_hovered(self, samples: object) -> None:
        values = hovered_samples(samples)
        if self.ctx is not None:
            self.ctx.hover(values, origin=self.uid or None, kind=hover_kind(values))

    def _on_sample_open_requested(self, sample: str) -> None:
        if sample and self.ctx is not None:
            self.ctx.request(OpenSampleCard(sample), origin=self.uid or "")

    def _on_subset_requested(self, payload: object, new_window: bool) -> None:
        if new_window:
            request_pan_subset_window(self, payload, "cluster")
        else:
            request_pan_subset(self, payload, "cluster")
