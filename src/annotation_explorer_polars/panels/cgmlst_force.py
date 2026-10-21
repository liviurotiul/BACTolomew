"""The cgMLST MST panel (plan 4.6): the window's sample set as a minimum
spanning tree of the cgMLST distances, drawn as a force graph.

Multi-instance and heavy. Hover works both ways. Its display options are
the "Color nodes by" and "Shape nodes by" traits.
"""
from __future__ import annotations

from typing import Any, Dict, Optional

from PyQt6 import QtWidgets

from .base import CGMLST_FORCE, ExplorePanel


class CgmlstForcePanel(ExplorePanel):
    KIND = CGMLST_FORCE
    TITLE = "cgMLST MST"
    HEAVY = True

    def __init__(self, window: Any, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setObjectName("CgmlstForcePanel")
        self._window = window
        self.plot: Any = None
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

    def set_plot(self, plot: Any) -> None:
        """The CgmlstForcePlotPanel (built by the window as this panel's child)."""
        self.plot = plot
        self.layout().addWidget(plot)
        plot.sampleHovered.connect(self._on_sample_hovered)
        plot.styleChanged.connect(self.titleDetailChanged)

    def title_detail(self) -> str:
        """The trait the nodes are colored by ("by Country")."""
        return self.state_detail({"color_by": self.plot.color_trait() if self.plot is not None else ""})

    @staticmethod
    def state_detail(state: Dict[str, Any]) -> str:
        """``title_detail`` from a saved state (PanelSpec.state_detail)."""
        trait = str(state.get("color_by") or "")
        return f"by {trait}" if trait else ""

    def save_state(self) -> Dict[str, Any]:
        if self.plot is None:
            return {}
        return {"color_by": self.plot.color_trait(), "shape_by": self.plot.shape_trait()}

    def restore_state(self, state: Dict[str, Any]) -> None:
        if self.plot is None or not state:
            return
        color, shape = state.get("color_by"), state.get("shape_by")
        self.plot.set_style_traits(
            str(color) if color is not None else None,
            str(shape) if shape is not None else None,
        )

    def on_shown(self) -> None:
        # A graph whose data changed while hidden is built when shown.
        if self.plot is not None:
            self.plot.ensure_fresh()

    def apply_result(self, filters: Any, result: Any) -> None:
        del filters
        self._window._refresh_cgmlst_force_plot(self.plot, result.samples)

    def show_hover(self, hover: Any) -> None:
        samples = list(getattr(hover, "samples", ()) or ())
        if len(samples) > 1:
            self.plot.highlight_samples(samples)
        else:
            self.plot.highlight_sample(samples[0] if samples else None)

    def _on_sample_hovered(self, sample: object) -> None:
        value = str(sample).strip() if sample is not None else ""
        if self.ctx is not None:
            self.ctx.hover([value] if value else [], origin=self.uid or None)
