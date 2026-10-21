"""Panels of one explore window: instances, placement, menus and the palette.

``PanelManager`` owns the docking area (``DockHost``) of one explore window and
the panel instances in it. It owns the tab and group menus (plan 3.5), the
panel palette and placement rules (3.6), Reopen Closed Panel, Maximize, Reset
Layout, the built-in presets and saved layouts (3.8), the empty state (3.1)
and the panel list used by the layout store (3.11). No keyboard shortcuts are added (D8): everything is reachable through
menus, buttons and drag.

Multi-instance kinds (map, the pan panels and
Annotations; plan 3.6, 8.5 and 8.10): the first
instance has the kind as uid and is hidden on close, so it reopens in place,
like a single-instance panel. Further instances (Duplicate Panel, Split
Right/Down from the tab menu, opening the kind again) get a suffixed uid; a
close disposes them and deletes their dock, and Reopen Closed Panel rebuilds
them from their saved display options beside where they were.

Document panels (neighbourhood and alignment comparisons, plan 8.7) are
opened with their data through ``open_document``: beside their source (plan
3.6 rule 3), one per comparison (opening the same one again focuses it),
disposed on close and reopened from memory by Reopen Closed Panel. They are
not part of saved layouts (3.11) and survive Reset Layout. The palette lists
only the open ones.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

from PyQt6 import QtCore, QtWidgets

from ..docking import layouts
from ..docking.dock_host import DockHost, Placement
from .base import GENE_PRESENCE_ABSENCE, PHYLOGENY, PanelFrame, PanelRegistry, PanelSet, PanelSpec
from .palette import PanelPalette

CLOSED_STACK_LIMIT = 10

# Source panel : a document panel split to its right (open_document).
PLOT_PANEL_SPLIT = (40, 60)

# Panels the empty state offers as one-click buttons (the most used ones).
_EMPTY_STATE_KINDS = ("sample_tree", "annotations", "geo")

_SIDE_LABELS = (("left", "Left"), ("right", "Right"), ("top", "Above"), ("bottom", "Below"))


@dataclass
class ClosedPanel:
    """One entry of the Reopen Closed Panel stack (plan 3.6 rule 5)."""

    uid: str
    kind: str
    state: Dict[str, Any] = field(default_factory=dict)
    title: str = ""
    # An extra instance is deleted on close: reopen it beside this panel,
    # as a tab of its group ("tab") or split on that side of it.
    beside: Optional[str] = None
    where: str = "tab"
    # A document panel's data (in memory only), to reopen it with.
    document: Any = None


class PanelManager(PanelSet):
    """Owns the docking area of one explore window and the panels in it."""

    activePanelChanged = QtCore.pyqtSignal(str)
    # An instance was added (a duplicate) or removed for good (a closed duplicate).
    instancesChanged = QtCore.pyqtSignal()
    # The whole docking area was rearranged (a preset, a saved layout,
    # Classic). The host emits no viewOpened / viewClosed while it restores.
    layoutApplied = QtCore.pyqtSignal()

    def __init__(
        self,
        window: QtWidgets.QWidget,
        specs: "Sequence[PanelSpec] | PanelRegistry",
        *,
        ctx: Any = None,
    ) -> None:
        super().__init__(window, specs, ctx=ctx)
        self._window = window
        self.host = DockHost()
        self._custom_titles: Dict[str, str] = {}
        self._default_titles: Dict[str, str] = {}
        self._closed_stack: List[ClosedPanel] = []
        self._last_active_by_kind: Dict[str, str] = {}
        for spec in self.registry.specs():
            if not spec.document:
                self._create_frame(spec, spec.uid)
        self.host.set_empty_page(self._build_empty_page())
        self.host.viewClosed.connect(self._on_view_closed)
        self.host.viewOpened.connect(self._on_view_opened)
        self.host.viewCloseRequested.connect(self._on_view_close_requested)
        self.host.activeViewChanged.connect(self._on_active_view_changed)
        self.host.tabMenuRequested.connect(self._show_tab_menu)
        self.host.groupMenuRequested.connect(self._show_group_menu)
        self.host.addToGroupRequested.connect(lambda uid: self.show_palette(target=uid))
        # The host is silent while it restores a layout, and a new instance's
        # first placement is no open toggle.
        self.layoutApplied.connect(self.refresh_default_titles)
        self.instancesChanged.connect(self.refresh_default_titles)
        self._view_menu: Optional[QtWidgets.QMenu] = None
        self._panels_menu: Optional[QtWidgets.QMenu] = None
        self._reopen_action = None
        self._maximize_action = None

    def _on_frame_created(self, frame: PanelFrame) -> None:
        extra = frame.uid != frame.kind
        if frame.spec.document:
            title = self.document_title(frame.uid)
        else:
            title = frame.spec.title
            frame.titleDetailChanged.connect(self._on_title_detail_changed)
            frame.panelBuilt.connect(self._on_title_detail_changed)
        self._default_titles[frame.uid] = title
        self.host.register_view(frame.uid, title, frame, custom_close=extra)

    # Default titles: the kind's title while one instance of it is open;
    # with several open, each says what it shows ("Phylogeny · Core
    # genome", ExplorePanel.title_detail), numbered from 2 when two would
    # read the same. Titles the user set (rename_panel) are kept.
    def _on_title_detail_changed(self, uid: str) -> None:
        frame = self.frame(uid)
        if frame is not None:
            self.refresh_default_titles(frame.kind)

    def refresh_default_titles(self, kind: Optional[str] = None) -> None:
        """Retitle the panels of ``kind`` (else of every kind) by the rule above."""
        specs = [self.registry.spec(kind)] if kind is not None else self.registry.specs()
        for spec in specs:
            if spec is None or spec.document:
                continue
            uids = self.uids(spec.kind)
            open_uids = [uid for uid in uids if self.is_open(uid)]
            several = len(open_uids) > 1
            # Only open panels claim a title: a closed one (the kind's first
            # instance stays registered when closed) gets the plain title.
            used = {self._custom_titles[uid] for uid in open_uids if uid in self._custom_titles}
            for uid in open_uids + [uid for uid in uids if uid not in open_uids]:
                frame = self.frame(uid)
                base = spec.title
                title = base
                if uid in open_uids:
                    if several:
                        detail = frame.title_detail() if frame is not None else ""
                        if detail:
                            base = f"{spec.title} · {detail}"
                    title, number = base, 2
                    while title in used:
                        title = f"{base} {number}"
                        number += 1
                if uid in self._custom_titles:
                    self._default_titles[uid] = title
                    continue
                used.add(title)
                if self._default_titles.get(uid) != title:
                    self._default_titles[uid] = title
                    self.host.set_view_title(uid, title)

    # ------------------------------------------------------------- lookups
    def specs(self) -> List[PanelSpec]:
        return self.registry.specs()

    def spec(self, kind: str) -> Optional[PanelSpec]:
        return self.registry.spec(kind)

    def is_open(self, uid: str) -> bool:
        return self.host.has_view(uid) and self.host.is_open(uid)

    def is_available(self, kind: str) -> Tuple[bool, str]:
        """For a kind, or for the uid of one of its instances."""
        frame = self.frame(kind)
        return self.registry.available(frame.kind if frame is not None else kind)

    def is_extra_instance(self, uid: str) -> bool:
        frame = self.frame(uid)
        return frame is not None and frame.uid != frame.kind

    def closed_stack(self) -> List[str]:
        return [entry.uid for entry in self._closed_stack]

    def title(self, uid: str) -> str:
        if uid in self._custom_titles:
            return self._custom_titles[uid]
        if uid in self._default_titles:
            return self._default_titles[uid]
        spec = self.registry.spec(uid)
        return spec.title if spec is not None else uid

    def last_active(self, kind: str) -> Optional[str]:
        """The instance of ``kind`` that was active last (else the first open one)."""
        uid = self._last_active_by_kind.get(kind)
        if uid is not None and self.frame(uid) is not None and self.is_open(uid):
            return uid
        for other in self.uids(kind):
            if self.is_open(other):
                return other
        return kind if self.frame(kind) is not None else None

    # -------------------------------------------------------------- opening
    def open_panel(
        self,
        kind: str,
        *,
        to_side: bool = False,
        target: Optional[str] = None,
        placement: Optional[Placement] = None,
    ) -> bool:
        """Open or focus a panel, following the placement rules of plan 3.6.

        ``kind`` may also be the uid of an existing instance. A single-instance
        kind focuses its panel when it is open; a multi-instance kind whose
        first instance is open gets a new instance."""
        uid = kind
        if self.frame(uid) is None:
            return False
        kind = self.frame(uid).kind
        available, _reason = self.is_available(kind)
        if not available:
            return False
        if self.is_open(uid):
            if uid == kind and self.registry.is_multi_instance(kind) and (to_side or target or placement):
                # Another view of the same kind, where it was asked for.
                return self.new_instance(kind, placement=self._placement_for(to_side, target, placement)) is not None
            # Single-instance kinds (and an open instance) focus the existing panel.
            self.host.focus_view(uid)
            return True
        if self.frame(uid).spec.requires:
            return self._open_with_requirement(uid, to_side=to_side, target=target, placement=placement)
        if placement is None:
            if not to_side and target is None and self.host.is_placed(uid):
                placement = None  # reopen where it was last
            elif not to_side and target is None:
                placement = self._linked_pair_placement(kind) or self._placement_for(False, None, None)
            else:
                placement = self._placement_for(to_side, target, None)
        return self.host.open_view(uid, placement)

    # Kinds whose rows link (panels/row_link.py: the heatmap's rows in the
    # tree's order, level with the leaves at the tree's fit) and the side each
    # opens on next to the other one: the tree on the left, the heatmap on the
    # right.
    _LINKED_PAIR_SIDES = {PHYLOGENY: (GENE_PRESENCE_ABSENCE, "left"), GENE_PRESENCE_ABSENCE: (PHYLOGENY, "right")}

    def _linked_pair_placement(self, kind: str) -> Optional[Placement]:
        """A Phylogeny or Gene presence/absence opened by default while the
        other one is open: split beside it (a tab of the same group would
        never be shown with it, so the two would not link)."""
        pair = self._LINKED_PAIR_SIDES.get(kind)
        if pair is None:
            return None
        partner_kind, side = pair
        for partner in self.uids(partner_kind):
            if self.is_open(partner):
                return Placement(side, relative_to=partner)
        return None

    def _open_with_requirement(
        self,
        uid: str,
        *,
        to_side: bool,
        target: Optional[str],
        placement: Optional[Placement],
    ) -> bool:
        """Open a panel that cannot exist without another one (the Genome
        track without Annotations): the required panel opens first, where
        the request asked (else in place, else by the usual rules), and this
        one goes below it. An explicit request for this panel alone is
        honored when the required panel is open already; a plain reopen of a
        panel placed before goes back where it was."""
        required = self.frame(uid).spec.requires
        explicit = to_side or target is not None or placement is not None
        if self.frame(required) is None:
            return False
        if self.is_open(required):
            if explicit:
                return self.host.open_view(uid, self._placement_for(to_side, target, placement))
            if self.host.is_placed(uid):
                return self.host.open_view(uid, None)
        else:
            if not explicit and self.host.is_placed(uid) and self.host.is_placed(required):
                # Both were closed together (closing the required panel closes
                # this one): reopen both where they were.
                if not self.open_panel(required):
                    return False
                if not self.host.open_view(uid, None):
                    return False
                self._fix_squashed_pair(uid, required)
                return True
            if not self.open_panel(required, to_side=to_side, target=target, placement=placement):
                return False
        if self.host.open_view(uid, Placement("bottom", relative_to=required)):
            self.host.set_split_sizes(uid, list(layouts.CLASSIC_TABLE_TRACK))
            return True
        # The bottom split was refused (both in one group, nothing else in
        # it): a tab beside the required panel instead.
        return self.host.open_view(uid, Placement("tab", relative_to=required))

    def _fix_squashed_pair(self, uid: str, required: str) -> None:
        """After both panels came back where they were: if the stored split
        left one of them a sliver (they were closed while one group was
        maximized), give the pair the Classic proportions again."""
        if self.host.group_uids(uid) == self.host.group_uids(required):
            return
        sizes = self.host.split_sizes(uid)
        if len(sizes) == 2 and sizes == self.host.split_sizes(required) and sum(sizes) > 0:
            if min(sizes) < 0.2 * sum(sizes):
                self.host.set_split_sizes(uid, list(layouts.CLASSIC_TABLE_TRACK))

    def requirement_title(self, uid: str) -> str:
        """The title of the panel ``uid`` cannot exist without, else ""."""
        frame = self.frame(uid)
        required = frame.spec.requires if frame is not None else ""
        return self.title(required) if required and self.frame(required) is not None else ""

    def enforce_requirements(self) -> List[str]:
        """Close the open panels whose required panel is not open (a saved
        layout, a session or a preset that has the Genome track without
        Annotations). Nothing is opened: a restore never adds panels. The
        closed panels are not offered by Reopen Closed Panel. Returns the
        closed uids."""
        closed: List[str] = []
        for frame in list(self.frames()):
            required = frame.spec.requires
            if required and self.is_open(frame.uid) and not self.is_open(required):
                self.close_panel(frame.uid)
                closed.append(frame.uid)
        if closed:
            self._closed_stack = [e for e in self._closed_stack if e.uid not in closed]
        return closed

    def _placement_for(
        self,
        to_side: bool,
        target: Optional[str],
        placement: Optional[Placement],
    ) -> Placement:
        """Placement rules 1 and 2 of plan 3.6 (a tab of the active or target
        group, or split right of the active group). A default open while a
        comparison is the active tab splits beside it rather than over it."""
        if placement is not None:
            return placement
        active = self.host.active_view()
        if to_side:
            return Placement("right", relative_to=active, outer=active is None)
        if target is not None:
            return Placement("tab", relative_to=target)
        frame = self.frame(active) if active is not None else None
        if frame is not None and frame.spec.document:
            # The active tab is a comparison (a neighbourhood window holds
            # nothing else): a tab would hide the rows that feed the new
            # panel's hover, so split to its left instead (the Phylogeny
            # beside the rows, as the tree sits left of its heatmap).
            return Placement("left", relative_to=active)
        return Placement("tab", relative_to=active)

    def open_new_instance(self, kind: str, *, to_side: bool = False, target: Optional[str] = None) -> Optional[str]:
        """Palette and menus: the first instance if it is closed, else a new one."""
        if self.frame(kind) is None or not self.is_available(kind)[0]:
            return None
        if not self.is_open(kind) or not self.registry.is_multi_instance(kind):
            return kind if self.open_panel(kind, to_side=to_side, target=target) else None
        placement = None if to_side or target is not None else self._linked_pair_placement(kind)
        return self.new_instance(kind, placement=placement or self._placement_for(to_side, target, None))

    def new_instance(
        self,
        kind: str,
        *,
        placement: Optional[Placement] = None,
        state: Optional[Dict[str, Any]] = None,
        uid: Optional[str] = None,
        title: str = "",
    ) -> Optional[str]:
        """Add one more instance of a multi-instance kind (placed, shown, focused)."""
        spec = self.registry.spec(kind)
        if spec is None or spec.single_instance or spec.document or not spec.available()[0]:
            return None
        frame = self._create_instance(spec, uid=uid, state=state, title=title)
        if placement is None:
            placement = self._new_instance_placement(kind)
        if not self.host.open_view(frame.uid, placement):
            self.remove_instance(frame.uid)
            return None
        self.instancesChanged.emit()
        return frame.uid

    def _new_instance_placement(self, kind: str) -> Placement:
        """A new instance goes in a tab of the active panel, except over its
        linked partner (a heatmap over the Phylogeny would hide the tree and
        unlink the rows): then in a tab of an open panel of its own kind,
        else beside the partner."""
        active = self.host.active_view()
        pair = self._LINKED_PAIR_SIDES.get(kind)
        frame = self.frame(active) if active else None
        if pair is not None and frame is not None and frame.kind == pair[0]:
            for same in self.uids(kind):
                if self.is_open(same):
                    return Placement("tab", relative_to=same)
            return Placement(pair[1], relative_to=active)
        return Placement("tab", relative_to=active)

    def _create_instance(
        self,
        spec: PanelSpec,
        *,
        uid: Optional[str] = None,
        state: Optional[Dict[str, Any]] = None,
        title: str = "",
    ) -> PanelFrame:
        if uid is None or self.frame(uid) is not None or PanelRegistry.kind_of(uid) != spec.kind or uid == spec.kind:
            uid = PanelRegistry.new_uid(spec.kind, list(self._frames) + [e.uid for e in self._closed_stack])
        frame = self._create_frame(spec, uid, state)
        if title:
            self.rename_panel(uid, title)
        return frame

    # ------------------------------------------------------------ documents
    def open_document(
        self,
        kind: str,
        document: Any,
        *,
        beside: Optional[str] = None,
        state: Optional[Dict[str, Any]] = None,
    ) -> Optional[str]:
        """Open a document panel of ``kind`` beside ``beside`` (else beside
        the active group, plan 3.6 rule 3) and focus it. A comparison that is
        open already (same ``document.key``) is focused instead; when that
        panel's load failed or found nothing, it loads again first. Further
        panels of a kind that is open already join that panel's group as
        tabs, so comparisons collect in one place instead of each splitting
        the window again."""
        existing = self.find_document(kind, getattr(document, "key", None))
        if existing is not None and self.is_open(existing):
            panel = self.panel(existing)
            retry = getattr(panel, "retry_with", None)
            if callable(retry):
                retry(document)
            self.host.focus_view(existing)
            return existing
        frame = self._create_document_frame(kind, document, state=state)
        if frame is None:
            return None
        sibling = self._open_document_of_kind(kind, exclude=frame.uid)
        if sibling is not None:
            placement = Placement("tab", relative_to=sibling)
        else:
            placement = self.beside_placement(beside)
        if not self.host.open_view(frame.uid, placement):
            self.remove_instance(frame.uid)
            return None
        if placement.where == "right":
            # These were large windows: give the new group most of the split.
            self.host.set_split_sizes(frame.uid, list(PLOT_PANEL_SPLIT))
        self.instancesChanged.emit()
        return frame.uid

    def _open_document_of_kind(self, kind: str, *, exclude: str = "") -> Optional[str]:
        """An open document panel of ``kind`` (the active one first)."""
        active = self.host.active_view()
        uids = [uid for uid in self.uids(kind) if uid != exclude and self.is_open(uid)]
        if active in uids:
            return active
        return uids[-1] if uids else None

    def beside_placement(self, source: Optional[str] = None) -> Placement:
        """Placement rule 3 (plan 3.6): beside ``source`` (else the active
        panel): a tab of the group to its right, or a new group split to its
        right when there is none."""
        if source is None or not self.is_open(source):
            source = self.host.active_view()
        if source is None:
            return Placement("tab")
        right = self.host.neighbor_group(source, "right")
        if right is not None:
            return Placement("tab", relative_to=right)
        return Placement("right", relative_to=source)

    def _reopen_document(self, entry: ClosedPanel) -> Optional[str]:
        sibling = self._open_document_of_kind(entry.kind)
        if entry.beside and self.is_open(entry.beside):
            placement = Placement(entry.where, relative_to=entry.beside)
        elif sibling is not None:
            placement = Placement("tab", relative_to=sibling)
        else:
            placement = self.beside_placement(None)
        existing = self.find_document(entry.kind, getattr(entry.document, "key", None))
        if existing is not None:
            self.host.focus_view(existing)
            return existing
        frame = self._create_document_frame(entry.kind, entry.document, uid=entry.uid, state=entry.state)
        if frame is None:
            return None
        if entry.title:
            self.rename_panel(frame.uid, entry.title)
        if not self.host.open_view(frame.uid, placement):
            self.remove_instance(frame.uid)
            return None
        self.instancesChanged.emit()
        return frame.uid

    def duplicate_panel(self, uid: str, where: str = "tab") -> Optional[str]:
        """Duplicate Panel (``where="tab"``: a tab of the same group) and the
        tab menu's Split Right / Split Down for multi-instance kinds: a new
        instance with the same display options, beside ``uid``."""
        frame = self.frame(uid)
        if frame is None or not self.registry.is_multi_instance(frame.kind) or not self.is_open(uid):
            return None
        return self.new_instance(frame.kind, placement=Placement(where, relative_to=uid), state=frame.save_state())

    def remove_instance(self, uid: str) -> None:
        """Dispose an extra instance and delete its dock."""
        if self.frame(uid) is None or not self.is_extra_instance(uid):
            return
        kind = self.frame(uid).kind
        self._custom_titles.pop(uid, None)
        self._default_titles.pop(uid, None)
        self._drop_frame(uid)
        self.host.remove_view(uid)
        self.refresh_default_titles(kind)
        self.instancesChanged.emit()

    def focus_panel(self, kind: str) -> Optional[str]:
        """Show and focus the instance of ``kind`` used last (open it if closed)."""
        uid = self.last_active(kind)
        if uid is None:
            return None
        if not self.open_panel(uid):
            return None
        self.host.focus_view(uid)
        return uid

    def close_panel(self, uid: str) -> None:
        if self.is_extra_instance(uid):
            if self.is_open(uid):
                self._close_extra_instance(uid)
            return
        self.host.close_view(uid)

    def _on_view_close_requested(self, uid: str) -> None:
        if self.is_extra_instance(uid):
            self._close_extra_instance(uid)

    def _close_extra_instance(self, uid: str) -> None:
        frame = self.frame(uid)
        if frame is None:
            return
        beside = next((other for other in self.host.group_uids(uid) if other != uid and self.is_open(other)), None)
        where = "tab"
        if beside is None:
            # Alone in its group: remember a neighbouring group and the side
            # the panel was on, so it reopens in its own split there.
            for side, reopen_side in (("left", "right"), ("right", "left"), ("top", "bottom"), ("bottom", "top")):
                beside = self.host.neighbor_group(uid, side)
                if beside is not None:
                    where = reopen_side
                    break
        entry = ClosedPanel(
            uid=uid,
            kind=frame.kind,
            state=frame.save_state(),
            title=self._custom_titles.get(uid, ""),
            beside=beside,
            where=where,
            document=self.document(uid),
        )
        self.remove_instance(uid)
        self._push_closed(entry)

    def reopen_closed(self) -> Optional[str]:
        while self._closed_stack:
            entry = self._closed_stack.pop()
            if not self.is_available(entry.kind)[0]:
                continue
            if self.frame(entry.uid) is not None:
                if self.is_open(entry.uid):
                    continue
                self.open_panel(entry.uid)
                return entry.uid
            if entry.document is not None:
                uid = self._reopen_document(entry)
                if uid is not None:
                    return uid
                continue
            if self.registry.is_multi_instance(entry.kind):
                if entry.beside and self.is_open(entry.beside):
                    placement = Placement(entry.where, relative_to=entry.beside)
                else:
                    placement = Placement("tab", relative_to=self.host.active_view())
                uid = self.new_instance(
                    entry.kind,
                    placement=placement,
                    state=entry.state,
                    uid=entry.uid,
                    title=entry.title,
                )
                if uid is not None:
                    return uid
        return None

    def _push_closed(self, entry: ClosedPanel) -> None:
        self._closed_stack = [e for e in self._closed_stack if e.uid != entry.uid]
        self._closed_stack.append(entry)
        del self._closed_stack[:-CLOSED_STACK_LIMIT]

    # -------------------------------------------------------------- layouts
    def _available_kinds(self) -> List[str]:
        return [spec.kind for spec in self.registry.specs() if self.registry.available(spec.kind)[0]]

    def build_default(self) -> None:
        layouts.build_classic(self.host, self._available_kinds())
        self.enforce_requirements()
        self.layoutApplied.emit()

    def preset_available(self, key: str) -> Tuple[bool, str]:
        """A preset is offered when at least one of its main panels is."""
        preset = layouts.preset(key)
        if preset is None:
            return False, "Unknown layout"
        available = set(self._available_kinds())
        if any(kind in available for kind in preset.needs_any):
            return True, ""
        names = ", ".join(self.registry.spec(k).title for k in preset.needs_any if self.registry.spec(k) is not None)
        return False, f"Not available for this run: needs {names or 'panels this run does not have'}"

    def reset_layout(self) -> None:
        """Rebuild Classic with the default panel titles. Touches only the
        docking area: never the rail, the filters, the sample set or a panel's
        display options (such as the table's view mode). Extra instances are
        closed (Classic has one of each). Open document panels stay open: a
        column on the right holds them."""
        self.apply_preset(layouts.PRESET_CLASSIC)

    def apply_preset(self, key: str) -> bool:
        """Arrange the docking area as a built-in layout (plan 3.8), like
        Reset Layout does for Classic: default titles, one instance of each
        kind, panels the preset leaves out closed, open document panels kept
        in a column on the right. False when the preset is not available."""
        preset = layouts.preset(key)
        if preset is None or not self.preset_available(key)[0]:
            return False
        documents = self._open_documents()
        for uid in [f.uid for f in self.frames() if f.uid != f.kind and f.uid not in documents]:
            self.remove_instance(uid)
        for uid in list(self._custom_titles):
            self.rename_panel(uid, "")
        self.host.clear_layout()
        preset.build(self.host, self._available_kinds())
        self.enforce_requirements()
        self._place_documents_right(documents)
        self._closed_stack.clear()
        self.layoutApplied.emit()
        return True

    def apply_saved_layout(self, dock_state: bytes, panels: Sequence[Dict[str, object]]) -> bool:
        """Arrange the docking area as a saved layout (a user layout or the
        default for new windows): its panels, their display options and
        titles. Open document panels stay open in a column on the right.
        Panels of kinds this run lacks are closed. If the layout cannot be
        applied, the current one is put back and False is returned."""
        if not dock_state:
            return False
        # Nothing open yet (a window starting with the default layout): on
        # failure the area is left empty and the caller builds Classic.
        had_layout = bool(self.host.open_uids())
        before_state = self.host.save_layout() if had_layout else b""
        before_panels = self.panel_list()
        documents = self._open_documents()
        ok = self._restore_saved(dock_state, panels, documents)
        if not ok:
            if not (had_layout and self._restore_saved(before_state, before_panels, documents)):
                self.host.clear_layout()
                self._drop_unplaced_instances()
                if had_layout:
                    self.build_default()
                    self._place_documents_right(documents)
            self.layoutApplied.emit()
            return False
        self._closed_stack.clear()
        self.layoutApplied.emit()
        return True

    def _restore_saved(
        self,
        dock_state: bytes,
        panels: Sequence[Dict[str, object]],
        documents: Sequence[str],
    ) -> bool:
        for uid in list(self._custom_titles):
            if uid not in documents:
                self.rename_panel(uid, "")
        ok = self.restore_layout(dock_state, panels)
        if ok:
            self.close_unavailable()
            ok = self.has_open_panels()
        if ok:
            self._place_documents_right([uid for uid in documents if not self.is_open(uid)])
        return ok

    def close_unavailable(self) -> None:
        """Close the open panels of kinds this run does not have (a layout
        saved with another run may name them)."""
        for frame in list(self.frames()):
            if self.is_open(frame.uid) and not self.registry.available(frame.kind)[0]:
                self.close_panel(frame.uid)

    def has_open_panels(self) -> bool:
        """True when a panel other than a document panel is open: a layout
        made of document panels alone (or of nothing) cannot be applied."""
        return any(self.is_open(f.uid) for f in self.frames() if not f.spec.document)

    def _open_documents(self) -> List[str]:
        return [f.uid for f in self.frames() if f.spec.document and self.is_open(f.uid)]

    def _place_documents_right(self, documents: Sequence[str]) -> None:
        first: Optional[str] = None
        for uid in documents:
            if self.frame(uid) is None:
                continue
            if first is None:
                self.host.open_view(uid, Placement("right", outer=True))
                first = uid
            else:
                self.host.open_view(uid, Placement("tab", relative_to=first))
        if first is not None:
            # ADS gives a group added beside a splitter's existing groups
            # only its minimum width; a fixed share keeps it usable.
            self.host.set_split_share(first, layouts.DOCUMENT_COLUMN_SHARE)

    def panel_list(self) -> List[Dict[str, object]]:
        out: List[Dict[str, object]] = []
        for frame in self.frames():
            if frame.spec.document:
                continue  # computed at run time, never restored (plan 3.11)
            params: Dict[str, object] = dict(frame.save_state())
            if frame.uid in self._custom_titles:
                params["title"] = self._custom_titles[frame.uid]
            out.append({"uid": frame.uid, "kind": frame.kind, "params": params})
        return out

    def apply_panel_list(self, panels: Sequence[Dict[str, object]], *, create_instances: bool = True) -> None:
        """Panel params from a saved layout; unknown kinds are dropped. With
        ``create_instances`` the extra instances it lists are created (not yet
        placed: the dock state places them)."""
        for entry in panels:
            kind = str(entry.get("kind") or "")
            spec = self.registry.spec(kind)
            if spec is None or spec.document:
                continue
            params = entry.get("params")
            params = dict(params) if isinstance(params, dict) else {}
            uid = str(entry.get("uid") or kind)
            if uid != kind:
                if not create_instances or spec.single_instance or PanelRegistry.kind_of(uid) != kind:
                    continue
                if self.frame(uid) is None:
                    self._create_frame(spec, uid)
            frame = self.frame(uid)
            if frame is None:
                continue
            title = str(params.pop("title", "") or "").strip()
            if title:
                self.rename_panel(uid, title)
            frame.restore_state(params)

    def restore_layout(self, dock_state: bytes, panels: Sequence[Dict[str, object]]) -> bool:
        """Saved layout first; False means the caller should build Classic."""
        self.apply_panel_list(panels)
        if not dock_state:
            self._drop_unplaced_instances()
            return False
        # A layout with nothing open is valid too (the empty state shows).
        ok = self.host.restore_layout(dock_state)
        if ok and self.enforce_requirements() and not self.has_open_panels():
            # The layout held nothing but a panel that cannot exist alone (the
            # Genome track without Annotations): treat it as unusable so the
            # caller builds Classic or keeps the previous layout, rather than
            # showing an empty window the user never saved.
            ok = False
        if ok:
            self._place_listed_panels(panels)
        if not ok:
            self.host.clear_layout()
        self._drop_unplaced_instances()
        self.layoutApplied.emit()
        return ok

    def _place_listed_panels(self, panels: Sequence[Dict[str, object]]) -> None:
        """Open the panels a layout lists with a ``place`` ({"where",
        "relative_to"}) that its dock state left closed or out: a migrated
        neoPan brings its Gene presence/absence beside its Phylogeny
        (shell/layout_store.migrate_layout)."""
        for entry in panels:
            place = entry.get("place")
            if not isinstance(place, dict):
                continue
            uid = str(entry.get("uid") or entry.get("kind") or "")
            anchor = str(place.get("relative_to") or "")
            where = str(place.get("where") or "right")
            if self.frame(uid) is None or self.is_open(uid) or not self.is_open(anchor):
                continue
            if not self.is_available(uid)[0]:
                continue
            self.host.open_view(uid, Placement(where, relative_to=anchor))

    def _drop_unplaced_instances(self) -> None:
        # Document panels are kept: a layout never holds them, and the caller
        # puts open ones back (apply_saved_layout).
        for uid in [f.uid for f in self.frames() if f.uid != f.kind and not f.spec.document]:
            if not self.is_open(uid):
                self.remove_instance(uid)

    def rename_panel(self, uid: str, title: str) -> None:
        title = title.strip()
        frame = self.frame(uid)
        if frame is None:
            return
        default = self._default_titles.get(uid) or frame.spec.title
        if not title or title == default:
            self._custom_titles.pop(uid, None)
            self.host.set_view_title(uid, default)
        else:
            self._custom_titles[uid] = title
            self.host.set_view_title(uid, title)
        if not frame.spec.document:
            # Another instance may no longer need (or now need) a number.
            self.refresh_default_titles(frame.kind)

    # -------------------------------------------------------------- palette
    def _palette_hint(self, target: Optional[str]) -> str:
        """The palette's footer: where a default Open puts the panel here
        (it follows _placement_for)."""
        side = "Open to the Side splits it to the right."
        if target is not None:
            return f"Open adds the panel as a tab of this group; {side}"
        active = self.host.active_view()
        frame = self.frame(active) if active is not None else None
        if frame is not None and frame.spec.document:
            return f"Open splits the panel to the left of the comparison; {side}"
        return f"Open adds the panel as a tab of the active group; {side}"

    def show_palette(self, target: Optional[str] = None) -> None:
        """Open the panel palette. ``target``: open into that panel's group."""
        palette = PanelPalette(
            [spec for spec in self.specs() if not spec.document],
            is_open=self.is_open,
            available=self.is_available,
            title="New panel here" if target else "New panel",
            parent=self._window,
            hint=self._palette_hint(target),
            documents=[
                (frame.uid, self.title(frame.uid), frame.spec.group)
                for frame in self.frames()
                if frame.spec.document and self.is_open(frame.uid)
            ],
        )

        def _open(kind: str, to_side: bool) -> None:
            if self.is_document(kind):
                # An open comparison: show it (in the target group when the
                # palette came from a group's "+").
                if target is not None and target != kind and self.is_open(target):
                    self.host.move_view(kind, Placement("right" if to_side else "tab", relative_to=target))
                self.host.focus_view(kind)
            elif to_side and target is not None and self.is_open(target):
                # From a group's "+": to the side of that group.
                if self.is_open(kind) and self.registry.is_multi_instance(kind):
                    self.new_instance(kind, placement=Placement("right", relative_to=target))
                else:
                    self.open_panel(kind, placement=Placement("right", relative_to=target))
            else:
                self.open_new_instance(kind, to_side=to_side, target=target)

        palette.opened.connect(_open)
        palette.exec()
        palette.deleteLater()

    # ---------------------------------------------------------------- menus
    def install_menus(self, menu_bar: QtWidgets.QMenuBar) -> QtWidgets.QMenu:
        view = menu_bar.addMenu("View")
        view.setObjectName("ViewMenu")
        # The panel picker is Window > New Panel... (Ctrl+Shift+P), owned by
        # the window.
        self._panels_menu = view.addMenu("Panels")
        self._panels_menu.aboutToShow.connect(self._fill_panels_menu)
        self._reopen_action = view.addAction("Reopen Closed Panel")
        self._reopen_action.triggered.connect(lambda _c=False: self.reopen_closed())
        self._maximize_action = view.addAction("Maximize Active Group")
        self._maximize_action.setCheckable(True)
        self._maximize_action.triggered.connect(self._toggle_maximize_active)
        view.addSeparator()
        reset_action = view.addAction("Reset Layout")
        reset_action.setToolTip("Restore the default arrangement of panels; filters are kept")
        reset_action.triggered.connect(lambda _c=False: self.reset_layout())
        view.aboutToShow.connect(self._update_view_menu)
        self._view_menu = view
        return view

    def _update_view_menu(self) -> None:
        if self._reopen_action is not None:
            self._reopen_action.setEnabled(bool(self._reopenable()))
        if self._maximize_action is not None:
            self._maximize_action.setChecked(self.host.is_maximized())
            self._maximize_action.setEnabled(
                self.host.is_maximized() or len(self.host.groups()) > 1
            )

    def _reopenable(self) -> List[str]:
        return [
            e.uid
            for e in self._closed_stack
            if not self.is_open(e.uid) and self.is_available(e.kind)[0]
        ]

    def _fill_panels_menu(self) -> None:
        menu = self._panels_menu
        if menu is None:
            return
        menu.clear()
        for frame in self.frames():
            uid = frame.uid
            action = menu.addAction(self.title(uid))
            action.setCheckable(True)
            action.setChecked(self.is_open(uid))
            available, reason = self.is_available(uid)
            action.setEnabled(available)
            required = self.requirement_title(uid)
            if not available:
                action.setToolTip(reason)
            elif required:
                action.setToolTip(f"Opens with {required}; closing {required} closes it")
            action.toggled.connect(
                lambda checked, _uid=uid: self.open_panel(_uid) if checked else self.close_panel(_uid)
            )

    def _toggle_maximize_active(self, checked: bool) -> None:
        if not checked:
            self.host.end_maximize()
            return
        active = self.host.active_view()
        if active is not None:
            self.host.maximize_group(active)

    def build_tab_menu(self, uid: str) -> QtWidgets.QMenu:
        """Tab context menu (plan 3.5). No Float or Detach entry exists.

        For multi-instance kinds Split Right / Split Down open a copy of the
        panel in a new group on that side (as an editor split does) and
        Duplicate Panel adds a copy as a tab of the same group; for the other
        kinds the split entries move the panel."""
        menu = QtWidgets.QMenu(self._window)
        menu.setObjectName("PanelTabMenu")
        group = self.host.group_uids(uid)
        alone = len(group) <= 1
        frame = self.frame(uid)
        multi = frame is not None and self.registry.is_multi_instance(frame.kind)
        split_right = menu.addAction("Split Right")
        split_down = menu.addAction("Split Down")
        if multi:
            split_right.setToolTip("Open a copy of this panel in a new group on the right")
            split_down.setToolTip("Open a copy of this panel in a new group below")
            split_right.triggered.connect(lambda _c=False: self.duplicate_panel(uid, "right"))
            split_down.triggered.connect(lambda _c=False: self.duplicate_panel(uid, "bottom"))
        else:
            split_right.setEnabled(not alone)
            split_down.setEnabled(not alone)
            split_right.triggered.connect(lambda _c=False: self.host.move_view(uid, Placement("right", relative_to=uid)))
            split_down.triggered.connect(lambda _c=False: self.host.move_view(uid, Placement("bottom", relative_to=uid)))
        move_menu = menu.addMenu("Move to Group")
        any_move = False
        for side, label in _SIDE_LABELS:
            neighbor = self.host.neighbor_group(uid, side)
            action = move_menu.addAction(label)
            if neighbor is not None:
                action.triggered.connect(
                    lambda _c=False, n=neighbor: self.host.move_view(uid, Placement("tab", relative_to=n))
                )
                any_move = True
            elif not alone:
                action.triggered.connect(
                    lambda _c=False, s=side: self.host.move_view(uid, Placement(s, relative_to=uid))
                )
                any_move = True
            else:
                action.setEnabled(False)
        move_menu.addSeparator()
        new_group = move_menu.addAction("New Group")
        new_group.setEnabled(not alone)
        new_group.triggered.connect(lambda _c=False: self.host.move_view(uid, Placement("right", relative_to=uid)))
        move_menu.setEnabled(any_move or not alone)
        if multi:
            duplicate = menu.addAction("Duplicate Panel")
            duplicate.setToolTip("Open a copy of this panel as a tab of the same group")
            duplicate.triggered.connect(lambda _c=False: self.duplicate_panel(uid, "tab"))
        rename = menu.addAction("Rename...")
        rename.triggered.connect(lambda _c=False: self._rename_interactive(uid))
        menu.addSeparator()
        close = menu.addAction("Close")
        close.triggered.connect(lambda _c=False: self.close_panel(uid))
        close_others = menu.addAction("Close Others")
        # The panel kept open also keeps the panel it requires (closing
        # Annotations would close the Genome track the user chose to keep).
        keep = {uid, frame.spec.requires if frame is not None else ""}
        close_others.setEnabled(any(other not in keep for other in group))
        close_others.triggered.connect(
            lambda _c=False: [self.close_panel(other) for other in self.host.group_uids(uid) if other not in keep]
        )
        close_all = menu.addAction("Close All in Group")
        close_all.triggered.connect(lambda _c=False: [self.close_panel(k) for k in self.host.group_uids(uid)])
        return menu

    def build_group_menu(self, uid: str) -> QtWidgets.QMenu:
        """Group menu (the group strip's button or a right-click on its empty
        strip; plan 3.5). ``uid`` is the group's current panel."""
        menu = QtWidgets.QMenu(self._window)
        menu.setObjectName("PanelGroupMenu")
        group = self.host.group_uids(uid)
        alone = len(group) <= 1
        open_here = menu.addAction("Open Panel Here...")
        open_here.triggered.connect(lambda _c=False: self.show_palette(target=uid))
        reopen = menu.addAction("Reopen Closed Panel")
        reopen.setEnabled(bool(self._reopenable()))
        reopen.triggered.connect(lambda _c=False: self.reopen_closed())
        menu.addSeparator()
        split_right = menu.addAction("Split Group Right")
        split_right.setToolTip("Move the visible panel into a new group on the right")
        split_right.setEnabled(not alone)
        split_right.triggered.connect(lambda _c=False: self.host.move_view(uid, Placement("right", relative_to=uid)))
        split_down = menu.addAction("Split Group Down")
        split_down.setToolTip("Move the visible panel into a new group below")
        split_down.setEnabled(not alone)
        split_down.triggered.connect(lambda _c=False: self.host.move_view(uid, Placement("bottom", relative_to=uid)))
        maximize = menu.addAction("Restore Group" if self.host.is_maximized() else "Maximize Group")
        maximize.setEnabled(self.host.is_maximized() or len(self.host.groups()) > 1)
        maximize.triggered.connect(lambda _c=False: self.host.toggle_maximize(uid))
        menu.addSeparator()
        close_group = menu.addAction("Close Group")
        close_group.triggered.connect(lambda _c=False: [self.close_panel(k) for k in self.host.group_uids(uid)])
        return menu

    def _show_tab_menu(self, uid: str, pos: QtCore.QPoint) -> None:
        menu = self.build_tab_menu(uid)
        menu.exec(pos)
        menu.deleteLater()

    def _show_group_menu(self, uid: str, pos: QtCore.QPoint) -> None:
        menu = self.build_group_menu(uid)
        menu.exec(pos)
        menu.deleteLater()

    def _rename_interactive(self, uid: str) -> None:
        current = self.title(uid)
        text, ok = QtWidgets.QInputDialog.getText(
            self._window,
            "Rename panel",
            "Panel title (leave empty for the default):",
            text=current,
        )
        if ok:
            self.rename_panel(uid, text)

    # ----------------------------------------------------------- empty state
    def _build_empty_page(self) -> QtWidgets.QWidget:
        page = QtWidgets.QWidget()
        page.setObjectName("DockEmptyState")
        outer = QtWidgets.QVBoxLayout(page)
        outer.addStretch(1)
        box = QtWidgets.QWidget(page)
        box.setFixedWidth(560)
        layout = QtWidgets.QVBoxLayout(box)
        layout.setSpacing(10)
        title = QtWidgets.QLabel("No panels are open", box)
        title.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        title.setStyleSheet("font-size: 16px; font-weight: 600;")
        layout.addWidget(title)
        hint = QtWidgets.QLabel("Open a panel, or restore the default layout. The filters on the left still apply.", box)
        hint.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        hint.setWordWrap(True)
        hint.setProperty("role", "muted")
        layout.addWidget(hint)
        row = QtWidgets.QHBoxLayout()
        row.addStretch(1)
        open_button = QtWidgets.QPushButton("New panel...", box)
        open_button.setObjectName("EmptyStateOpenPanel")
        open_button.clicked.connect(lambda _c=False: self.show_palette())
        row.addWidget(open_button)
        for kind in _EMPTY_STATE_KINDS:
            spec = self.registry.spec(kind)
            if spec is None:
                continue
            button = QtWidgets.QPushButton(spec.title, box)
            button.setObjectName(f"EmptyStateOpen_{kind}")
            button.clicked.connect(lambda _c=False, _k=kind: self.open_panel(_k))
            row.addWidget(button)
        row.addStretch(1)
        layout.addLayout(row)
        restore = QtWidgets.QPushButton("Restore default layout", box)
        restore.setObjectName("EmptyStateRestoreLayout")
        restore.clicked.connect(lambda _c=False: self.reset_layout())
        restore_row = QtWidgets.QHBoxLayout()
        restore_row.addStretch(1)
        restore_row.addWidget(restore)
        restore_row.addStretch(1)
        layout.addLayout(restore_row)
        outer.addWidget(box, alignment=QtCore.Qt.AlignmentFlag.AlignHCenter)
        outer.addStretch(2)
        return page

    # ------------------------------------------------------------ internals
    def _on_view_closed(self, uid: str) -> None:
        frame = self.frame(uid)
        if frame is not None and not frame.spec.document:
            self.refresh_default_titles(frame.kind)
        if frame is None or self.is_extra_instance(uid):
            return  # an extra instance pushed its own entry when it was closed
        self._push_closed(ClosedPanel(uid=uid, kind=frame.kind, state=frame.save_state()))
        # A panel that cannot exist without this one closes with it (closing
        # Annotations closes the Genome track). Pushed after it, so Reopen
        # Closed Panel brings both back. The host is silent while it
        # restores; enforce_requirements covers restores and presets.
        for dependent in self.registry.dependents(uid):
            if self.frame(dependent) is not None and self.is_open(dependent):
                self.close_panel(dependent)

    def _on_view_opened(self, uid: str) -> None:
        self._closed_stack = [e for e in self._closed_stack if e.uid != uid]
        frame = self.frame(uid)
        if frame is not None and not frame.spec.document:
            self.refresh_default_titles(frame.kind)

    def _on_active_view_changed(self, uid: str) -> None:
        frame = self.frame(uid)
        if frame is not None:
            self._last_active_by_kind[frame.kind] = uid
        self.activePanelChanged.emit(uid)
