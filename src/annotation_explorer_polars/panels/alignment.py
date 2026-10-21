"""The alignment panel (plan 4.6, 8.7): a gene card's alignment result docked
in the card's owner explore window.

A document panel: each result is its own instance, opened with the aligned
sequences by the ``panel.open_alignment`` command and never saved with the
layout (3.11). Its rows follow the window's one sample set (the consensus row
stays that of the whole alignment); a hovered row hovers its sample in the
window and rows of samples hovered elsewhere are tinted.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional

from PyQt6 import QtWidgets

from .base import ALIGNMENT, ExplorePanel


@dataclass
class AlignmentDocument:
    """One alignment result (kept in memory, never saved)."""

    alignment_text: str
    status_text: str = ""
    title: str = "Sequence alignment"
    key: Any = None


class AlignmentPanel(ExplorePanel):
    KIND = ALIGNMENT
    TITLE = "Sequence alignment"
    HEAVY = True
    NEEDS_PARAMS = True

    def __init__(self, document: AlignmentDocument, parent: Optional[QtWidgets.QWidget] = None) -> None:
        from ..annotation_detail_dialog import AlignmentView

        super().__init__(parent)
        self.setObjectName("AlignmentPanel")
        self.document = document
        self._hovering = False
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(0)
        self.view = AlignmentView(
            self,
            alignment_text=document.alignment_text,
            status_text=document.status_text,
            sample_hover_callback=self._on_sample_hovered,
            panel_mode=True,
        )
        layout.addWidget(self.view)

    def title(self) -> str:
        return self.document.title or self.TITLE

    def apply_result(self, filters: Any, result: Any) -> None:
        del filters
        self.view.set_sample_scope(result.samples)

    def show_hover(self, hover: Any) -> None:
        self.view.set_linked_samples(tuple(getattr(hover, "samples", ()) or ()))

    def _on_sample_hovered(self, sample: Optional[str]) -> None:
        if self.ctx is None:
            return
        value = str(sample or "").strip()
        self._hovering = bool(value)
        self.ctx.hover([value] if value else [], origin=self.uid or None)

    def unbind(self) -> None:
        # A row hovered when the panel closes is not left hovered.
        ctx = self.ctx
        if ctx is not None and self._hovering:
            ctx.clear_hover(origin=self.uid)
        super().unbind()

    def save_state(self) -> Dict[str, Any]:
        return {}
