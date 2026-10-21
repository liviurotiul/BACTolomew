from __future__ import annotations

import re
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple


_FASTGEAR_LINEAGE_PALETTE = [
    "#3f6ea5",
    "#35826f",
    "#7564a7",
    "#a97832",
    "#a45f73",
    "#4f7f8a",
    "#8d6b56",
    "#5f7a55",
    "#8a667d",
    "#4f6f9c",
    "#b26b5e",
    "#4e8778",
    "#7b6ba7",
    "#9a7b3f",
    "#557c87",
    "#6b7280",
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

# fastGEAR cluster colors: the app's cluster palette (Cluster 1 blue, 2 red,
# 3 green, 4 purple, ...), the same in every tree panel.
CLUSTER_HUES: Tuple[str, ...] = tuple(_METADATA_STRIP_PALETTE)


def cluster_colors(clusters: Sequence[str], sort_key: Optional[Callable[[str], object]] = None) -> Dict[str, str]:
    """fastGEAR cluster colors: ``CLUSTER_HUES`` in cluster order."""
    ordered = sorted({str(value) for value in clusters if str(value).strip()}, key=sort_key or str)
    return {value: CLUSTER_HUES[index % len(CLUSTER_HUES)] for index, value in enumerate(ordered)}


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
