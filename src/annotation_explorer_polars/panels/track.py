"""The Genome track panel (plan 4.6): paired with the first Annotations panel.

Single instance. It cannot exist without that Annotations panel (feedback
2026-10-07, ``PanelSpec.requires``): opening the track opens Annotations too,
closing Annotations closes the track, and a restored layout or preset that
has the track alone closes it.

The track draws the contig of the paired Annotations panel's selection
(Contig, Around matches) or the context of its selected row (All matches);
the pairing and the table and track scroll sync are the window's
(``IgvControllerMixin``). The panel itself only hosts the view: results,
selection and hover reach the track through the paired Annotations panel.
The optional dimming of other samples' hits (plan 4.6) is not done: the
track only ever shows one sample's contig.
"""
from __future__ import annotations

from typing import Any, Optional

from PyQt6 import QtWidgets

from .base import TRACK, ExplorePanel


class TrackPanel(ExplorePanel):
    KIND = TRACK
    TITLE = "Genome track"
    SINGLE_INSTANCE = True

    def __init__(self, window: Any, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setObjectName("TrackPanel")
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(window._build_track_widget())
