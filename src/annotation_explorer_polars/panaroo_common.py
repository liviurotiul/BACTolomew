from __future__ import annotations

import re
from pathlib import Path
from typing import List, Tuple


_FASTGEAR_LINEAGE_PALETTE = [
    "#1d4ed8",
    "#dc2626",
    "#059669",
    "#7c3aed",
    "#ea580c",
    "#0f766e",
    "#be185d",
    "#4f46e5",
    "#65a30d",
    "#b45309",
    "#0891b2",
    "#4338ca",
    "#b91c1c",
    "#15803d",
    "#a21caf",
    "#0369a1",
]

_METADATA_STRIP_PALETTE = [
    "#1d4ed8",
    "#dc2626",
    "#059669",
    "#7c3aed",
    "#ea580c",
    "#0f766e",
    "#be185d",
    "#4f46e5",
    "#65a30d",
    "#b45309",
    "#0891b2",
    "#4338ca",
    "#b91c1c",
    "#15803d",
    "#a21caf",
    "#0369a1",
    "#f59e0b",
    "#2563eb",
    "#16a34a",
    "#e11d48",
]

_NEO_TREE_RESISTANCE_TOKEN = "__neo_tree_resistance__"
_NEO_TREE_MGE_TOKEN = "__neo_tree_mge__"
_NEO_TREE_VIRULENCE_TOKEN = "__neo_tree_virulence__"
_NEO_TREE_PANAROO_TOKEN = "__neo_tree_panaroo__"
_NEO_TREE_PLASMID_TOKEN = "__neo_tree_plasmids__"


def _natural_sort_key(value: object) -> Tuple[Tuple[int, object], ...]:
    text = str(value or "").strip()
    if not text:
        return ((1, ""),)
    parts = re.split(r"(\d+)", text)
    key: List[Tuple[int, object]] = []
    for part in parts:
        if not part:
            continue
        if part.isdigit():
            key.append((0, int(part)))
        else:
            key.append((1, part.casefold()))
    return tuple(key)


def _file_signature(path: Path) -> Tuple[str, int, int]:
    resolved = path.resolve()
    stat = resolved.stat()
    return str(resolved), int(stat.st_mtime_ns), int(stat.st_size)
