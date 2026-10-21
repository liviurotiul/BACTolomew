"""The filter rail as window chrome bound to the context (plan 4.4, D11).

The rail is not a panel: it is never in the registry, never docked and never
closed while its window is open. It writes the context's draft (the rail
widgets and their handlers are still the window's ``FiltersControllerMixin``;
moving them here is a follow-up of the plan's 7.1) and re-renders from the context: the gene search summary and
the one-time D3 hint from each applied result, and the subset chip whenever
``set_subset`` replaces it (a clade, cluster or plasmid from the pan views, or
"Open sample(s) in table").
"""
from __future__ import annotations

from typing import Any, Optional

from PyQt6 import QtCore

from ..session import bus as bus_events

OWNER = "rail"


class FilterRail(QtCore.QObject):
    def __init__(self, window: Any, ctx: Any, parent: Optional[QtCore.QObject] = None) -> None:
        super().__init__(parent if parent is not None else window)
        self._window = window
        self.ctx = ctx
        ctx.bus.subscribe(OWNER, bus_events.RESULT_APPLIED, self._on_result_applied)
        ctx.subsetChanged.connect(self._on_subset_changed)

    def _on_result_applied(self, payload: Any) -> None:
        filters, result = payload
        self._window._update_table_search_hint(result.filtered_ann)
        self._window._update_sample_search_hint(filters, result)

    def _on_subset_changed(self, _subset: Any, _origin: str) -> None:
        self._window._render_metadata_filters()
