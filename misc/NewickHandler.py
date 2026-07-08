# phylo_utils.py
from __future__ import annotations
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Callable, Dict, List, Tuple, Optional

from pipeline.Utils import norm_label

# Minimal tree node
@dataclass(slots=True)
class Node:
    name: str = ""
    length: float = 1.0
    children: Optional[List["Node"]] = None
    def __post_init__(self):
        if self.children is None:
            self.children = []

def load_newick_text(newick_path: str | Path | None, sample_cols: List[str]) -> str:
    if newick_path:
        return Path(newick_path).read_text()
    return star_newick(sample_cols)

def star_newick(leaves: Iterable[str]) -> str:
    clean = [norm_label(s) for s in leaves]
    return "(" + ",".join(clean) + ");"

def parse_newick(s: str) -> Node:
    s = re.sub(r"\s+", "", s)
    i = 0
    def name():
        nonlocal i
        j = i
        while i < len(s) and s[i] not in ",():;":
            i += 1
        return s[j:i]
    def length():
        nonlocal i
        if i < len(s) and s[i] == ":":
            i += 1
            j = i
            while i < len(s) and s[i] not in ",():;":
                i += 1
            return float(s[j:i])
        return None
    def subtree():
        nonlocal i
        if s[i] == "(":
            i += 1
            ch = [subtree()]
            while s[i] == ",":
                i += 1
                ch.append(subtree())
            if s[i] != ")":
                raise ValueError("Newick parse error")
            i += 1
            return Node(name(), length(), ch)
        else:
            return Node(name(), length(), [])
    r = subtree()
    if i < len(s) and s[i] == ";":
        i += 1
    return r

def layout_tree(root: Node) -> tuple[Dict[str, float], List[Tuple[float,float,float,float]], float, Dict[str, float]]:
    leaf_y: Dict[str, float] = {}
    segs: List[Tuple[float,float,float,float]] = []
    labels: List[Tuple[str,float]] = []
    yk = [0]

    def walk(node: Node, x0: float = 0.0) -> Tuple[float, float]:
        if not node.children:
            y = yk[0]; yk[0] += 1
            nm = norm_label(node.name or f"L{y}")
            x = x0 + (node.length or 0.0)
            leaf_y[nm] = y
            labels.append((nm, y))
            return x, y
        ys = []; mx = x0
        for ch in node.children:
            x1 = x0 + (ch.length or 0.0)
            xc, yc = walk(ch, x1)
            ys.append((x1, yc))
            segs.append((x0, yc, x1, yc))
            mx = max(mx, xc, x1)
        ymin = min(y for _, y in ys); ymax = max(y for _, y in ys)
        segs.append((x0, ymin, x0, ymax))
        return mx, (ymin + ymax) / 2

    xmax, _ = walk(root, 0.0)
    return leaf_y, segs, xmax, dict(labels)

def align_samples_to_tree(sample_cols: List[str], leaf_y: Dict[str, float], *, norm_label: Callable[[str], str]) -> List[str]:
    norm_map = {norm_label(c): c for c in sample_cols}
    ordered = [norm_map[n] for n in sorted(leaf_y, key=lambda k: leaf_y[k]) if n in norm_map]
    extras = [c for c in sample_cols if c not in ordered]
    return ordered + extras
