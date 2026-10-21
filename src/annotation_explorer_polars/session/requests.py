"""Typed open requests of an explore window (plan 4.4, ``ctx.request``).

Panels never call each other or the window directly to open something: they
hand one of these to their context, and the window that owns the context
carries it out. Requests always target the context's own window.

Phase 4 has the requests that today's panels make, plus ``OpenPanel``;
Phase 6 adds ``ClosePanel`` (a neighbourhood panel's Ctrl+W). Neighbourhood
panels are opened by the ``window.open_neighbourhood`` command (a gene card)
or ``OpenPlasmidContext`` (a pan view), so no ``OpenNeighbourhood`` request
is needed. Phase 7 adds ``OpenSubsetInNewWindow``: a clade, cluster or
plasmid opened in a new explore window of its own. It and
``OpenPlasmidContext`` (which shows the comparison in a new neighbourhood
window) are the requests that do not act on the context's own window; they
start another, independent one. ``ShowAnnotationsTable`` makes a window
without an Annotations panel show one (a neighbourhood panel's "Open
sample(s) in table").
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Tuple

from ..model.filter_state import SampleSubset


@dataclass(frozen=True)
class OpenSampleCard:
    sample: str


@dataclass(frozen=True)
class OpenGeneCard:
    """A gene card for one annotation row (the row's fields as a dict)."""

    row: Tuple[Tuple[str, Any], ...]

    @classmethod
    def from_row(cls, row: Dict[str, Any]) -> "OpenGeneCard":
        return cls(tuple((str(k), v) for k, v in dict(row).items()))

    def row_dict(self) -> Dict[str, Any]:
        return dict(self.row)


@dataclass(frozen=True)
class OpenPanel:
    """Open a panel of ``kind`` (a new instance for multi-instance kinds when
    one is open already). ``params`` are panel params such as a view mode;
    ``beside`` is the uid of the panel it opens next to (placement rule 3)."""

    kind: str
    params: Tuple[Tuple[str, Any], ...] = ()
    beside: Optional[str] = None

    def params_dict(self) -> Dict[str, Any]:
        return dict(self.params)


@dataclass(frozen=True)
class FocusPanel:
    """Show and focus a panel of ``kind`` (the one used last), applying ``params``."""

    kind: str
    params: Tuple[Tuple[str, Any], ...] = ()

    def params_dict(self) -> Dict[str, Any]:
        return dict(self.params)


@dataclass(frozen=True)
class ClosePanel:
    """Close the panel ``uid`` of this window (a panel's own Ctrl+W, plan 3.7)."""

    uid: str


@dataclass(frozen=True)
class OpenPlasmidContext:
    """The plasmid context viewer for a plasmid cluster (a pan view's request payload)."""

    payload: Tuple[Tuple[str, Any], ...] = field(default=())

    @classmethod
    def from_payload(cls, payload: Dict[str, Any]) -> "OpenPlasmidContext":
        return cls(tuple((str(k), v) for k, v in dict(payload).items()))

    def payload_dict(self) -> Dict[str, Any]:
        return dict(self.payload)


@dataclass(frozen=True)
class ShowAnnotationsTable:
    """Show the window's sample set in a table: when no Annotations panel is
    open, open the Sample tree and Annotations beside the panel ``beside``.
    A window that shows Annotations already does not change."""

    beside: Optional[str] = None


@dataclass(frozen=True)
class OpenSubsetInNewWindow:
    """Open a new explore window whose rail carries ``subset`` as its locked
    window subset (plan 3.10). The new window is independent from then on."""

    subset: SampleSubset


def params_tuple(params: Optional[Dict[str, Any]]) -> Tuple[Tuple[str, Any], ...]:
    return tuple(sorted((str(k), v) for k, v in dict(params or {}).items()))
