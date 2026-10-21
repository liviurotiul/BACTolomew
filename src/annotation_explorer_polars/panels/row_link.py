"""Row order of a window's Phylogeny and Gene presence/absence panels.

When an explore window shows a Phylogeny panel and a Gene presence/absence
panel, the heatmap's rows are the tree's leaves, in the tree's order (genomes
outside the sample set are dimmed rows, as the tree greys them). Only the
order is linked: the heatmap always fits every row top to bottom
(pan_view.PanHeatmapWidget) and never scrolls or zooms with the tree.
``TreeRowLink`` is owned by the window (one per window, created with its
``PanelManager``); it only sees the panels of its own manager, so windows
never link to each other. Panels never call each other: the link reads the
leader's leaf order and hands it to the followers.

- Leader: the first visible, built Phylogeny panel of the window. Further
  Phylogeny panels are not linked.
- Followers: every visible, built Gene presence/absence panel of the window.
  With no leader (or a circular tree) they show their own rows (the order
  of the run's default tree). A linked heatmap follows whatever tree source
  and root its leader draws (it takes the leader's drawn leaf order).
- Heights: each panel fits its own rows. The heatmap's rows fill its panel
  from the top and the tree fits its own viewport (below its toolbar and
  legend), so the rows are not level with the leaves; nothing leaves a band
  blank to line them up. Zooming or scrolling the tree leaves the heatmap as
  it is.
- A sync runs when the tree's rows (leaf order, layout) change or a panel is
  shown, hidden, opened or closed, never on the tree's scroll or zoom alone.
  Updates are coalesced (a zero-delay timer): one sync per event-loop turn.

Hidden panels (a background tab, Maximize, closed) are left alone and synced
again when shown (``linkStateChanged`` of the panels).
"""
from __future__ import annotations

from typing import Any, List, Optional

from PyQt6 import QtCore, sip

from .base import GENE_PRESENCE_ABSENCE, PHYLOGENY


def _alive(obj: Any) -> bool:
    return obj is not None and not sip.isdeleted(obj)


class TreeRowLink(QtCore.QObject):
    """Links the rows of one window's Phylogeny and Gene presence/absence
    panels (see the module docstring)."""

    # After each sync: the leader's uid ("" for none) and the linked followers.
    synced = QtCore.pyqtSignal(str, list)

    def __init__(self, manager: Any, parent: Optional[QtCore.QObject] = None) -> None:
        super().__init__(parent)
        self._manager = manager
        self._timer = QtCore.QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(0)
        self._timer.timeout.connect(self.sync)
        self._leader_uid = ""
        self._followers: List[str] = []
        manager.layoutApplied.connect(self.schedule)
        manager.instancesChanged.connect(self.schedule)
        manager.panelBuilt.connect(self._on_panel_built)
        host = getattr(manager, "host", None)
        if host is not None:
            host.viewOpened.connect(self.schedule)
            host.viewClosed.connect(self.schedule)
            host.activeViewChanged.connect(self.schedule)
            host.maximizedChanged.connect(self.schedule)
        for uid in manager.uids():
            self._on_panel_built(uid)

    # ------------------------------------------------------------ wiring
    def _on_panel_built(self, uid: str) -> None:
        panel = self._manager.panel(uid)
        if panel is None or panel.property("treeRowLinked"):
            return
        kind = self._manager.frame(uid).kind
        if kind == PHYLOGENY:
            panel.rowsViewChanged.connect(self.schedule)
            panel.linkStateChanged.connect(self.schedule)
        elif kind == GENE_PRESENCE_ABSENCE:
            panel.linkStateChanged.connect(self.schedule)
        else:
            return
        panel.setProperty("treeRowLinked", True)
        self.schedule()

    def schedule(self, *_args: object) -> None:
        """Sync on the next event-loop turn (coalesced)."""
        if not self._timer.isActive():
            self._timer.start()

    # ------------------------------------------------------------ roles
    def _visible_panels(self, kind: str) -> List[Any]:
        out = []
        for uid in self._manager.uids(kind):
            panel = self._manager.panel(uid)
            if _alive(panel) and self._manager.is_visible(uid) and getattr(panel, "tab", None) is not None:
                out.append(panel)
        return out

    def leader(self) -> Optional[Any]:
        """The Phylogeny panel the heatmaps follow (None when none is shown)."""
        for panel in self._visible_panels(PHYLOGENY):
            return panel
        return None

    def followers(self) -> List[Any]:
        return self._visible_panels(GENE_PRESENCE_ABSENCE)

    def leader_uid(self) -> str:
        return self._leader_uid

    def linked_uids(self) -> List[str]:
        """The followers linked at the last sync."""
        return list(self._followers)

    # ------------------------------------------------------------ sync
    def sync(self) -> None:
        self._timer.stop()
        leader = self.leader()
        followers = self.followers()
        order = leader.row_order() if leader is not None else None
        for panel in followers:
            panel.follow_rows(order)
        linked = [panel.uid for panel in followers] if order is not None else []
        self._leader_uid = leader.uid if leader is not None and order is not None else ""
        self._followers = linked
        self.synced.emit(self._leader_uid, list(linked))
