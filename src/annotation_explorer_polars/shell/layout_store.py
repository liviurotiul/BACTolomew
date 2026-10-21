"""Layout persistence of explore windows (plan 3.11), in QSettings.

Keys:
- ``explore/v2/run/<sha1(output_dir)>/windows/<n>/...``: the saved layout of
  explore window number n of the run: format version, geometry, rail width
  and collapsed state, the docking area state blob and the panel list. Each
  window saves into the slot of its own number.
- ``explore/v2/run/<sha1(output_dir)>/session``: the slots of the windows that
  were open when the run was last left (JSON list, in window order); reopening
  the run opens one window per slot. Windows with a locked window subset are
  not listed (they start from a clade or the like at run time).
- ``explore/v2/lastWindow/...``: geometry, rail state and panel display
  options (no dock state, no panel titles) of the window saved last, for runs
  that have no saved layout yet. The pre-docking build kept these in one
  global group, so every run reopened with them; this keeps that behavior.
- ``explore/v2/default/...``: the arrangement chosen with "Set as Default for
  New Windows" (format version, dock state, panel list; plan 3.8). New
  Explore Window, windows opened from a clade or the like, and runs with no
  saved layout start with it instead of Classic. Only that menu entry writes
  it; "Use Classic for New Windows" removes it. Window size and rail state
  for runs with no saved layout still come from ``lastWindow``, which every
  save updates.
- ``explore/v2/layouts/<name>/...``: user layouts saved with "Save Layout
  As..." (the name, format version, dock state and panel list; no filters,
  no rail state, no geometry). The group name is the layout name,
  percent-encoded except for letters, digits, spaces and "_.-~".
- ``explore/v2/hints/<name>``: one-time hints already shown (``HINT_*`` keys).

Format versions: version 1 named the panels neoPan (``neopan``), Pan tree
(``pan_tree``), neoTree (``neotree``) and Gene presence (``gene_presence``).
Version 2 has Phylogeny and Gene presence/absence instead; a version 1 layout
is migrated in memory when it is loaded (``migrate_layout``) and written as
version 2 at its next save, so loading never writes the settings.

Only the layout is stored, never filters or searches (D5). The pre-docking
``annotationExplorer`` group is no longer read (its one-time migration went
with the legacy layout in Phase 10) and never written.
"""
from __future__ import annotations

import hashlib
import json
import secrets
import xml.etree.ElementTree as ElementTree
from urllib.parse import quote, unquote
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from PyQt6 import QtCore

FORMAT_VERSION = 2
# Panel kinds of format version 1 and the kind that replaces each.
PHYLOGENY_KIND = "phylogeny"
GENE_PRESENCE_ABSENCE_KIND = "gene_presence_absence"
RETIRED_KINDS = {
    "neopan": PHYLOGENY_KIND,
    "pan_tree": PHYLOGENY_KIND,
    "neotree": PHYLOGENY_KIND,
    "gene_presence": GENE_PRESENCE_ABSENCE_KIND,
}
# Their params that still mean something (the Pan tree's tree display
# options are the Phylogeny panel's); a custom title is always kept.
_KEPT_PARAMS = {"pan_tree": ("layout", "branch_scale")}
ROOT = "explore/v2"
LAST_WINDOW_GROUP = f"{ROOT}/lastWindow"
DEFAULT_GROUP = f"{ROOT}/default"
USER_LAYOUTS_GROUP = f"{ROOT}/layouts"
LAYOUT_NAME_MAX = 60
# Panel params that pick the run's data rather than arrange panels (the
# Phylogeny panel's tree source, root and collapsed subtrees, which name the
# run's genomes): kept in a run's window slots, left out of the default, the
# user layouts and the last window, which arrange panels only (plan 3.11,
# PHYLOGENY_REWORK section 4).
RUN_DATA_PARAMS = ("source", "root", "collapsed")

HINTS_GROUP = f"{ROOT}/hints"
# Shown once, the first time the sample search matches no sample while the
# same text would match genes or products (decision D3, plan 4.3).
HINT_SAMPLE_SEARCH_D3 = "sample_search_d3"

# A new explore window's start: the default layout (Classic) and a clean
# rail, instead of a saved slot (plan 3.10, New Explore Window).
START_DEFAULT_LAYOUT = "default"
# A window that starts with an empty docking area (no Classic): the
# neighbourhood window, whose comparison panel is opened into it right after
# it is built.
START_EMPTY_LAYOUT = "empty"

RAIL_MIN_WIDTH = 248
RAIL_MAX_WIDTH = 420
RAIL_DEFAULT_WIDTH = 288


def clamp_rail_width(value: object) -> int:
    try:
        width = int(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        width = RAIL_DEFAULT_WIDTH
    return max(RAIL_MIN_WIDTH, min(width, RAIL_MAX_WIDTH))


def run_key(output_dir: Optional[Path]) -> str:
    """Stable per-run key: sha1 of the run's resolved output folder."""
    text = ""
    if output_dir is not None:
        try:
            text = str(Path(output_dir).expanduser().resolve())
        except OSError:
            text = str(output_dir)
    return hashlib.sha1(text.encode("utf-8")).hexdigest()


def hint_was_shown(name: str) -> bool:
    return bool(QtCore.QSettings().value(f"{HINTS_GROUP}/{name}", False, type=bool))


def mark_hint_shown(name: str) -> None:
    settings = QtCore.QSettings()
    settings.setValue(f"{HINTS_GROUP}/{name}", True)
    settings.sync()


def _byte_array(value: object) -> QtCore.QByteArray:
    if isinstance(value, QtCore.QByteArray):
        return value
    if isinstance(value, (bytes, bytearray)):
        return QtCore.QByteArray(bytes(value))
    return QtCore.QByteArray()


def _kind_of(uid: str) -> str:
    return str(uid).split("-", 1)[0]


def _renamed(uid: str) -> Optional[str]:
    """The version 2 uid of a version 1 uid of a retired kind (else None):
    ``pan_tree`` becomes ``phylogeny``, ``pan_tree-91c0`` ``phylogeny-91c0``."""
    kind = _kind_of(uid)
    new_kind = RETIRED_KINDS.get(kind)
    if new_kind is None:
        return None
    return new_kind + str(uid)[len(kind):]


def _decode_dock_state(state: bytes) -> Tuple[Optional[bytes], bool]:
    """The dock state's XML and whether it was compressed (qCompress)."""
    data = bytes(state or b"")
    if data.lstrip().startswith(b"<"):
        return data, False
    xml = bytes(QtCore.qUncompress(QtCore.QByteArray(data)))
    return (xml, True) if xml.lstrip().startswith(b"<") else (None, True)


def _migrate_dock_state(state: bytes) -> Tuple[bytes, Dict[str, str], List[str], List[str]]:
    """Rename the docks of retired kinds in a version 1 dock state.

    Returns the new state, the new uid of every renamed dock (old uid to new
    uid), the new uids that are open and the old uids that were open. When several old docks become the
    same uid (neoPan and the Pan tree both become ``phylogeny``), the first
    open one keeps it; the others get a fresh suffixed uid and are closed, so
    the restore drops them and one Phylogeny stays open."""
    xml, compressed = _decode_dock_state(state)
    if xml is None:
        return state, {}, [], []
    try:
        root = ElementTree.fromstring(xml)
    except ElementTree.ParseError:
        return state, {}, [], []
    widgets = list(root.iter("Widget"))
    names = [str(w.get("Name") or "") for w in widgets]
    taken = {name for name in names if _renamed(name) is None}
    renames: Dict[str, str] = {}
    # Which old dock keeps each new uid: an already new one never moves; else
    # the first open one, else the first one.
    keeper: Dict[str, int] = {}
    for index, widget in enumerate(widgets):
        new = _renamed(names[index])
        if new is None or new in taken:
            continue
        is_open = str(widget.get("Closed") or "0") != "1"
        if new not in keeper or (is_open and str(widgets[keeper[new]].get("Closed") or "0") == "1"):
            keeper[new] = index
    used = set(taken) | set(keeper)
    opened: List[str] = []
    was_open = [
        names[i] for i, w in enumerate(widgets) if _renamed(names[i]) is not None and str(w.get("Closed") or "0") != "1"
    ]
    for index, widget in enumerate(widgets):
        old = names[index]
        new = _renamed(old)
        if new is None:
            continue
        if keeper.get(new) != index:
            kind = _kind_of(new)
            while True:
                new = f"{kind}-{secrets.token_hex(2)}"
                if new not in used:
                    break
            widget.set("Closed", "1")
        used.add(new)
        widget.set("Name", new)
        renames[old] = new
        if str(widget.get("Closed") or "0") != "1":
            opened.append(new)
    if not renames:
        return state, {}, [], []
    for area in root.iter("Area"):
        current = str(area.get("Current") or "")
        if current not in renames:
            continue
        new = renames[current]
        if new in opened:
            area.set("Current", new)
            continue
        # The area's current dock was closed by the migration: show another.
        others = [str(w.get("Name") or "") for w in area.iter("Widget") if str(w.get("Closed") or "0") != "1"]
        area.set("Current", others[0] if others else new)
    out = ElementTree.tostring(root, encoding="utf-8", xml_declaration=True)
    if compressed:
        out = bytes(QtCore.qCompress(QtCore.QByteArray(out), 9))
    return out, renames, opened, was_open


def migrate_layout(layout: "WindowLayout") -> "WindowLayout":
    """A version 1 layout as version 2: neoPan, the Pan tree and neoTree become
    Phylogeny, Gene presence becomes Gene presence/absence, in the dock state
    and in the panel list. A neoPan that was open showed its tree and its
    heatmap, so Gene presence/absence opens to the right of its Phylogeny
    when the migrated layout has none open (the panel entry's ``place``).
    The Pan tree's display options carry over; neoPan's and neoTree's are
    dropped (titles are kept). Other versions are returned unchanged."""
    if layout.version != 1:
        return layout
    dock_state, renames, opened, was_open = (
        _migrate_dock_state(layout.dock_state) if layout.dock_state else (b"", {}, [], [])
    )
    panels: List[Dict[str, Any]] = []
    seen: set = set()
    entries = [dict(entry) for entry in layout.panels if isinstance(entry, dict)]
    for entry in entries:
        uid = str(entry.get("uid") or entry.get("kind") or "")
        old_kind = str(entry.get("kind") or _kind_of(uid))
        if old_kind not in RETIRED_KINDS:
            if uid not in seen:
                panels.append(entry)
                seen.add(uid)
            continue
        new_uid = renames.get(uid) or _renamed(uid) or uid
        if new_uid in seen:
            continue  # a second old panel for one new uid (no dock state)
        params = entry.get("params")
        params = dict(params) if isinstance(params, dict) else {}
        kept = {key: params[key] for key in _KEPT_PARAMS.get(old_kind, ()) if key in params}
        if str(params.get("title") or "").strip():
            kept["title"] = params["title"]
        panels.append({"uid": new_uid, "kind": RETIRED_KINDS[old_kind], "params": kept})
        seen.add(new_uid)
    # Closed duplicates the dock state names: listed, so their docks exist
    # while the state is restored (and are dropped right after).
    for old, new in renames.items():
        if new not in seen:
            panels.append({"uid": new, "kind": _kind_of(new), "params": {}})
            seen.add(new)
    neopan_open = [old for old in was_open if _kind_of(old) == "neopan"]
    anchor = renames.get(neopan_open[0], "") if neopan_open else ""
    if anchor not in opened:
        # neoPan lost its uid to another open Phylogeny: the heatmap goes there.
        anchor = next((uid for uid in opened if _kind_of(uid) == PHYLOGENY_KIND), "")
    if anchor and not any(_kind_of(uid) == GENE_PRESENCE_ABSENCE_KIND for uid in opened):
        place = {"where": "right", "relative_to": anchor}
        entry = next((p for p in panels if p.get("uid") == GENE_PRESENCE_ABSENCE_KIND), None)
        if entry is None:
            panels.append({"uid": GENE_PRESENCE_ABSENCE_KIND, "kind": GENE_PRESENCE_ABSENCE_KIND, "params": {}, "place": place})
        else:
            entry["place"] = place
    return WindowLayout(
        geometry=layout.geometry,
        rail_width=layout.rail_width,
        rail_collapsed=layout.rail_collapsed,
        dock_state=dock_state,
        panels=panels,
        version=FORMAT_VERSION,
    )


@dataclass
class WindowLayout:
    """What one explore window saves (layout only, D5)."""

    geometry: bytes = b""
    rail_width: int = RAIL_DEFAULT_WIDTH
    rail_collapsed: bool = False
    dock_state: bytes = b""
    panels: List[Dict[str, Any]] = field(default_factory=list)
    version: int = FORMAT_VERSION


class LayoutStore:
    def __init__(self, settings: Optional[QtCore.QSettings] = None) -> None:
        self._settings = settings

    def _qsettings(self) -> QtCore.QSettings:
        return self._settings if self._settings is not None else QtCore.QSettings()

    @staticmethod
    def window_group(key: str, number: int = 1) -> str:
        return f"{ROOT}/run/{key}/windows/{int(number)}"

    def load_window(self, key: str, number: int = 1) -> Optional[WindowLayout]:
        return self._load_group(self.window_group(key, number))

    def _load_group(self, group: str) -> Optional[WindowLayout]:
        settings = self._qsettings()
        settings.beginGroup(group)
        try:
            if not settings.contains("version"):
                return None
            try:
                version = int(settings.value("version", 0))
            except (TypeError, ValueError):
                version = 0
            try:
                panels = json.loads(str(settings.value("panels", "[]") or "[]"))
            except (TypeError, ValueError):
                panels = []
            if not isinstance(panels, list):
                panels = []
            # An older format is migrated in memory only (written as the
            # current one at the next save).
            return migrate_layout(
                WindowLayout(
                    geometry=bytes(_byte_array(settings.value("geometry"))),
                    rail_width=clamp_rail_width(settings.value("railWidth", RAIL_DEFAULT_WIDTH)),
                    rail_collapsed=settings.value("railCollapsed", False, type=bool),
                    dock_state=bytes(_byte_array(settings.value("dockState"))),
                    panels=[p for p in panels if isinstance(p, dict)],
                    version=version,
                )
            )
        finally:
            settings.endGroup()

    def save_window(self, key: str, layout: WindowLayout, number: int = 1) -> None:
        self._save_group(self.window_group(key, number), layout)

    @staticmethod
    def session_key(key: str) -> str:
        return f"{ROOT}/run/{key}/session"

    def load_session(self, key: str) -> List[int]:
        """The window slots to reopen for this run (empty when none is recorded)."""
        try:
            raw = json.loads(str(self._qsettings().value(self.session_key(key), "[]") or "[]"))
        except (TypeError, ValueError):
            return []
        slots: List[int] = []
        for value in raw if isinstance(raw, list) else []:
            try:
                slot = int(value)
            except (TypeError, ValueError):
                continue
            if slot > 0 and slot not in slots:
                slots.append(slot)
        return slots

    def save_session(self, key: str, slots: List[int]) -> None:
        settings = self._qsettings()
        settings.setValue(self.session_key(key), json.dumps([int(slot) for slot in slots]))
        settings.sync()

    def load_last_window(self) -> Optional[WindowLayout]:
        """Fallback for a run with no saved layout: the last saved window's
        geometry, rail state and panel display options (never a dock state)."""
        return self._load_group(LAST_WINDOW_GROUP)

    def save_last_window(self, layout: WindowLayout) -> None:
        panels = []
        for entry in layout.panels:
            params = dict(entry.get("params") or {})
            params.pop("title", None)  # a renamed tab belongs to its run
            for key in RUN_DATA_PARAMS:
                params.pop(key, None)
            panels.append({**entry, "params": params})
        self._save_group(
            LAST_WINDOW_GROUP,
            WindowLayout(
                geometry=layout.geometry,
                rail_width=layout.rail_width,
                rail_collapsed=layout.rail_collapsed,
                dock_state=b"",
                panels=panels,
                version=layout.version,
            ),
        )

    # ------------------------------------------- default and user layouts
    def load_default(self) -> Optional[WindowLayout]:
        """The arrangement set with "Set as Default for New Windows", if any."""
        return self._load_group(DEFAULT_GROUP)

    def has_default(self) -> bool:
        return self._qsettings().contains(f"{DEFAULT_GROUP}/version")

    def save_default(self, layout: WindowLayout) -> None:
        self._save_arrangement(DEFAULT_GROUP, layout)

    def clear_default(self) -> None:
        settings = self._qsettings()
        settings.remove(DEFAULT_GROUP)
        settings.sync()

    @staticmethod
    def clean_layout_name(name: object) -> str:
        """A user layout name: one line, trimmed, at most LAYOUT_NAME_MAX characters."""
        text = " ".join(str(name or "").split())
        return text[:LAYOUT_NAME_MAX].strip()

    @staticmethod
    def _user_layout_group(name: str) -> str:
        return f"{USER_LAYOUTS_GROUP}/{quote(name, safe=' ')}"

    def layout_names(self) -> List[str]:
        """The saved user layouts, sorted by name (case-insensitive)."""
        settings = self._qsettings()
        settings.beginGroup(USER_LAYOUTS_GROUP)
        try:
            names = []
            for group in settings.childGroups():
                stored = settings.value(f"{group}/name")
                name = self.clean_layout_name(stored if stored else unquote(group))
                if name and settings.contains(f"{group}/version"):
                    names.append(name)
        finally:
            settings.endGroup()
        return sorted(set(names), key=lambda n: (n.casefold(), n))

    def existing_layout_name(self, name: str) -> Optional[str]:
        """The saved user layout ``name`` refers to: the same name, or one
        that differs only in case ("Map" and "map" are one layout)."""
        name = self.clean_layout_name(name)
        if not name:
            return None
        if self._qsettings().contains(f"{self._user_layout_group(name)}/version"):
            return name
        folded = name.casefold()
        return next((n for n in self.layout_names() if n.casefold() == folded), None)

    def has_user_layout(self, name: str) -> bool:
        return self.existing_layout_name(name) is not None

    def load_user_layout(self, name: str) -> Optional[WindowLayout]:
        existing = self.existing_layout_name(name)
        return self._load_group(self._user_layout_group(existing)) if existing else None

    def save_user_layout(self, name: str, layout: WindowLayout) -> str:
        """Save (or overwrite) a user layout; returns the name it was saved
        under. A layout whose name differs only in case is replaced, and the
        new spelling is kept."""
        name = self.clean_layout_name(name)
        if not name:
            raise ValueError("A layout needs a name")
        existing = self.existing_layout_name(name)
        if existing is not None and existing != name:
            self._qsettings().remove(self._user_layout_group(existing))
        self._save_arrangement(self._user_layout_group(name), layout, name=name)
        return name

    def delete_user_layout(self, name: str) -> bool:
        existing = self.existing_layout_name(name)
        if existing is None:
            return False
        settings = self._qsettings()
        settings.remove(self._user_layout_group(existing))
        settings.sync()
        return True

    def _save_arrangement(self, group: str, layout: WindowLayout, *, name: str = "") -> None:
        """The docking area only: dock state and panel list (no geometry, no
        rail state, no filters)."""
        settings = self._qsettings()
        settings.remove(group)
        settings.beginGroup(group)
        try:
            if name:
                settings.setValue("name", name)
            settings.setValue("version", int(layout.version))
            settings.setValue("dockState", QtCore.QByteArray(layout.dock_state))
            settings.setValue("panels", json.dumps(_arrangement_panels(layout.panels), sort_keys=True))
        finally:
            settings.endGroup()
        settings.sync()

    def _save_group(self, group: str, layout: WindowLayout) -> None:
        settings = self._qsettings()
        settings.beginGroup(group)
        try:
            settings.setValue("version", int(layout.version))
            settings.setValue("geometry", QtCore.QByteArray(layout.geometry))
            settings.setValue("railWidth", clamp_rail_width(layout.rail_width))
            settings.setValue("railCollapsed", bool(layout.rail_collapsed))
            settings.setValue("dockState", QtCore.QByteArray(layout.dock_state))
            settings.setValue("panels", json.dumps(layout.panels, sort_keys=True))
        finally:
            settings.endGroup()
        settings.sync()


def _arrangement_panels(panels: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """A layout's panel list without the params that pick run data."""
    out = []
    for entry in panels:
        params = entry.get("params")
        if isinstance(params, dict) and any(key in params for key in RUN_DATA_PARAMS):
            entry = {**entry, "params": {k: v for k, v in params.items() if k not in RUN_DATA_PARAMS}}
        out.append(entry)
    return out
