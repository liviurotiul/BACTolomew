"""Newick trees without Qt: the node type and parser, a writer, pruning (ape
``drop.tip``), midpoint and outgroup rooting, splits and the Robinson-Foulds
distance.

Branch support: a numeric label on an internal node (IQ-TREE UFBoot or SH-aLRT,
RAxML support) is read into ``TreeNode.support`` and never kept as a name. The
support belongs to the edge above the node, so rerooting moves each value with
its edge: the (split, support) set of a tree is the same before and after
rooting. Leaf labels are always names (a sample may be called "373").

The pipeline writes the rooted trees with it (``pipeline.TreesIndex``) and the
GUI's ``phylo.py`` imports the tree model back from here.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Callable, Dict, FrozenSet, List, Optional, Sequence, Set, Tuple


@dataclass
class TreeNode:
    name: str = ""
    children: List["TreeNode"] = field(default_factory=list)
    branch_length: Optional[float] = None
    # Support of the edge above this node (internal nodes only), as written in the file.
    support: Optional[float] = None

    @property
    def is_leaf(self) -> bool:
        return len(self.children) == 0


_NUMBER = re.compile(r"^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$")


class NewickParser:
    def __init__(self, text: str) -> None:
        self.text = text.strip()
        self.n = len(self.text)
        self.i = 0

    def parse(self) -> TreeNode:
        if not self.text:
            raise ValueError("Empty Newick string.")
        node = self._parse_subtree()
        self._skip_ws()
        if self._peek() == ";":
            self.i += 1
        return node

    def _peek(self) -> str:
        if self.i >= self.n:
            return ""
        return self.text[self.i]

    def _skip_ws(self) -> None:
        while self.i < self.n:
            ch = self.text[self.i]
            if ch.isspace():
                self.i += 1
            elif ch == "[":
                # Newick comment: skipped.
                end = self.text.find("]", self.i)
                self.i = self.n if end < 0 else end + 1
            else:
                break

    def _parse_subtree(self) -> TreeNode:
        self._skip_ws()
        if self._peek() == "(":
            self.i += 1
            children: List[TreeNode] = []
            while True:
                children.append(self._parse_subtree())
                self._skip_ws()
                ch = self._peek()
                if ch == ",":
                    self.i += 1
                    continue
                if ch == ")":
                    self.i += 1
                    break
                raise ValueError(f"Invalid Newick near index {self.i}.")
            label, quoted = self._parse_label()
            branch = self._parse_branch_length()
            node = TreeNode(children=children, branch_length=branch)
            if not quoted and _NUMBER.match(label):
                node.support = float(label)
            else:
                node.name = label
            return node

        label, _quoted = self._parse_label()
        branch = self._parse_branch_length()
        return TreeNode(name=label, branch_length=branch)

    def _parse_label(self) -> Tuple[str, bool]:
        self._skip_ws()
        if self.i >= self.n:
            return "", False
        if self._peek() == "'":
            self.i += 1
            parts: List[str] = []
            while self.i < self.n:
                ch = self.text[self.i]
                if ch == "'":
                    self.i += 1
                    if self.i < self.n and self.text[self.i] == "'":
                        parts.append("'")
                        self.i += 1
                        continue
                    break
                parts.append(ch)
                self.i += 1
            return "".join(parts).strip(), True

        start = self.i
        while self.i < self.n and self.text[self.i] not in ",():;[":
            self.i += 1
        return self.text[start:self.i].strip(), False

    def _parse_branch_length(self) -> Optional[float]:
        self._skip_ws()
        if self._peek() != ":":
            return None
        self.i += 1
        self._skip_ws()
        start = self.i
        while self.i < self.n and self.text[self.i] not in ",() ;[":
            self.i += 1
        raw = self.text[start:self.i].strip()
        if not raw:
            return None
        try:
            return float(raw)
        except Exception:
            return None


def parse_newick(text: str) -> TreeNode:
    return NewickParser(text).parse()


def read_newick(path) -> TreeNode:
    from pathlib import Path

    return parse_newick(Path(path).read_text(encoding="utf-8", errors="ignore"))


# ------------------------------------------------------------------ writing
_UNSAFE_LABEL = re.compile(r"[\s(),:;\[\]']")


def _label(name: str) -> str:
    if not name:
        return ""
    if _UNSAFE_LABEL.search(name):
        return "'" + name.replace("'", "''") + "'"
    return name


def _number(value: float) -> str:
    text = format(float(value), ".10g")
    return text


def _support_text(value: float) -> str:
    return str(int(value)) if float(value).is_integer() else _number(value)


def write_newick(root: TreeNode, *, support: bool = True) -> str:
    """Newick text of the tree (ending with ";\\n"). Internal nodes write their
    support as the label when there is one (else their name)."""

    def walk(node: TreeNode) -> str:
        if node.is_leaf:
            text = _label(node.name)
        else:
            inner = ",".join(walk(child) for child in node.children)
            if support and node.support is not None:
                label = _support_text(node.support)
            else:
                label = _label(node.name)
            text = f"({inner}){label}"
        if node.branch_length is not None and math.isfinite(node.branch_length):
            text += ":" + _number(node.branch_length)
        return text

    return walk(root) + ";\n"


# ------------------------------------------------------------------ helpers
def sum_branch_lengths(upper: Optional[float], lower: Optional[float]) -> Optional[float]:
    """The length of two edges merged into one (a missing length counts as 0
    unless both are missing)."""
    if upper is None and lower is None:
        return None
    return float(upper or 0.0) + float(lower or 0.0)


def _merged_support(upper: Optional[float], lower: Optional[float]) -> Optional[float]:
    # Two edges merged into one stand for the same split; the lower value is the one written.
    return lower if lower is not None else upper


def leaf_names(root: TreeNode) -> List[str]:
    out: List[str] = []
    stack = [root]
    while stack:
        node = stack.pop()
        if node.is_leaf:
            out.append(node.name)
        else:
            stack.extend(reversed(node.children))
    return out


def rename_leaves(root: TreeNode, mapping: Callable[[str], str]) -> TreeNode:
    stack = [root]
    while stack:
        node = stack.pop()
        if node.is_leaf:
            node.name = mapping(node.name)
        else:
            stack.extend(node.children)
    return root


def prune_to(
    root: TreeNode,
    keep: Set[str],
    leaf_key: Optional[Callable[[str], str]] = None,
) -> Optional[TreeNode]:
    """A copy of the tree with only the leaves whose key is in ``keep`` (ape
    ``drop.tip`` semantics): dropped leaves go, internal nodes left with one
    child are merged into it with their branch lengths summed, and the root
    moves down to the most recent common ancestor of the kept leaves. The
    topology and branch lengths between kept leaves are unchanged. None when
    no leaf is kept."""
    key = leaf_key or (lambda name: name)

    def walk(node: TreeNode) -> Optional[TreeNode]:
        if node.is_leaf:
            if key(node.name or "") not in keep:
                return None
            return TreeNode(name=node.name, branch_length=node.branch_length)
        children = [child for child in (walk(child) for child in node.children) if child is not None]
        if not children:
            return None
        if len(children) == 1:
            only = children[0]
            only.branch_length = sum_branch_lengths(node.branch_length, only.branch_length)
            if not only.is_leaf:
                only.support = _merged_support(node.support, only.support)
            return only
        return TreeNode(
            name=node.name,
            children=children,
            branch_length=node.branch_length,
            support=node.support,
        )

    pruned = walk(root)
    if pruned is not None:
        # The edge above the new root is not part of the kept tree.
        pruned.branch_length = root.branch_length
        pruned.support = None
    return pruned


# ------------------------------------------------------------------ rooting
class _Graph:
    """The tree as an undirected graph (nodes in file preorder)."""

    def __init__(self, root: TreeNode) -> None:
        self.nodes: List[TreeNode] = []
        self.parent: Dict[int, Optional[int]] = {}
        index: Dict[int, int] = {}
        stack: List[Tuple[TreeNode, Optional[int]]] = [(root, None)]
        while stack:
            node, parent = stack.pop()
            index[id(node)] = len(self.nodes)
            self.parent[len(self.nodes)] = parent
            self.nodes.append(node)
            me = index[id(node)]
            stack.extend((child, me) for child in reversed(node.children))
        self.adjacency: Dict[int, List[int]] = {i: [] for i in range(len(self.nodes))}
        for i, parent in self.parent.items():
            if parent is not None:
                self.adjacency[parent].append(i)
                self.adjacency[i].append(parent)
        self.leaves = [i for i, node in enumerate(self.nodes) if node.is_leaf]

    def edge_child(self, a: int, b: int) -> int:
        """The node below the edge a-b in the file's rooting."""
        return b if self.parent.get(b) == a else a

    def edge_length(self, a: int, b: int, *, use_counts: bool) -> float:
        if use_counts:
            return 1.0
        value = self.nodes[self.edge_child(a, b)].branch_length
        if value is None or not math.isfinite(value) or value < 0.0:
            return 0.0
        return float(value)

    def edge_support(self, a: int, b: int) -> Optional[float]:
        child = self.nodes[self.edge_child(a, b)]
        return None if child.is_leaf else child.support

    def rooted_on_edge(
        self,
        upper: int,
        lower: int,
        len_upper: Optional[float],
        len_lower: Optional[float],
        *,
        clean_lengths: bool = False,
    ) -> TreeNode:
        """A copy rooted on the edge upper-lower, the root splitting it into
        len_upper (to upper) and len_lower (to lower). Support moves with its edge;
        a node left with one child (the old root) is merged into it. With
        clean_lengths, missing or invalid lengths are written as 0."""
        root_support = self.edge_support(upper, lower)

        def build(node_id: int, from_id: int, branch: Optional[float], support: Optional[float]) -> TreeNode:
            src = self.nodes[node_id]
            built = TreeNode(name=src.name, branch_length=branch)
            if not src.is_leaf:
                built.support = support
            for neighbor in self.adjacency[node_id]:
                if neighbor == from_id:
                    continue
                if clean_lengths:
                    length: Optional[float] = self.edge_length(node_id, neighbor, use_counts=False)
                else:
                    length = self.nodes[self.edge_child(node_id, neighbor)].branch_length
                built.children.append(build(neighbor, node_id, length, self.edge_support(node_id, neighbor)))
            if len(built.children) == 1:
                only = built.children[0]
                only.branch_length = sum_branch_lengths(built.branch_length, only.branch_length)
                if not only.is_leaf:
                    only.support = _merged_support(built.support, only.support)
                return only
            if not built.children:
                built.support = None
            return built

        new_root = TreeNode(name="", branch_length=None)
        new_root.children = [
            build(upper, lower, len_upper, root_support),
            build(lower, upper, len_lower, root_support),
        ]
        return new_root


def midpoint_root(root: TreeNode) -> TreeNode:
    """A copy of the tree rooted at the midpoint of its longest leaf to leaf
    path. Missing or invalid lengths count as 0; a tree without any positive
    length is rooted on edge counts instead (and keeps its file lengths). A
    node left with one child by the reroot (the old root) is merged into it
    with the lengths summed. Support values move with their edges."""
    graph = _Graph(root)
    if len(graph.leaves) < 2:
        return prune_to(root, {graph.nodes[i].name for i in graph.leaves}) or root
    use_counts = not any(
        graph.edge_length(i, p, use_counts=False) > 0.0
        for i, p in graph.parent.items()
        if p is not None
    )

    def farthest(start: int) -> Tuple[int, Dict[int, float], Dict[int, int]]:
        dist = {start: 0.0}
        parent: Dict[int, int] = {}
        todo = [start]
        while todo:
            current = todo.pop()
            for neighbor in graph.adjacency[current]:
                if neighbor not in dist:
                    dist[neighbor] = dist[current] + graph.edge_length(current, neighbor, use_counts=use_counts)
                    parent[neighbor] = current
                    todo.append(neighbor)
        # Ties go to the first leaf in file order, so the result is stable.
        return max(graph.leaves, key=lambda leaf: (dist[leaf], -leaf)), dist, parent

    end_a, _dist, _parent = farthest(graph.leaves[0])
    end_b, dist, parent = farthest(end_a)
    half = dist[end_b] / 2.0
    # Walk from end_b towards end_a to the edge that holds the midpoint.
    lower = end_b
    upper = parent.get(lower)
    while upper is not None and dist[upper] > half:
        lower = upper
        upper = parent.get(lower)
    if upper is None:
        return prune_to(root, {graph.nodes[i].name for i in graph.leaves}) or root
    if use_counts:
        return graph.rooted_on_edge(upper, lower, None, None)
    return graph.rooted_on_edge(upper, lower, half - dist[upper], dist[lower] - half, clean_lengths=True)


def root_on_outgroup(root: TreeNode, outgroup: Sequence[str]) -> Tuple[Optional[TreeNode], Optional[str]]:
    """(tree rooted on the edge that separates the outgroup from the rest, None)
    or (None, reason) when it cannot be: "missing: a, b" (outgroup taxa not in
    the tree) or "not_monophyletic". The root halves that edge, as Gubbins does."""
    wanted = [name for name in dict.fromkeys(str(x).strip() for x in outgroup) if name]
    graph = _Graph(root)
    leaf_set = {graph.nodes[i].name for i in graph.leaves}
    missing = [name for name in wanted if name not in leaf_set]
    if not wanted or missing:
        return None, "missing: " + ", ".join(missing or ["(empty)"])
    target = set(wanted)
    if target == leaf_set:
        return None, "not_monophyletic"
    below = _leaves_below(graph)
    for node_id, parent in graph.parent.items():
        if parent is None:
            continue
        side = below[node_id]
        if side == target or (leaf_set - side) == target:
            length = graph.nodes[node_id].branch_length
            half = None if length is None else length / 2.0
            return graph.rooted_on_edge(parent, node_id, half, half), None
    return None, "not_monophyletic"


def root_split(rooted: TreeNode) -> Optional[Tuple[FrozenSet[str], float]]:
    """(leaves of the root's first child, the share of the root edge on that side) of a
    tree rooted on an edge (two root children); None for any other root."""
    if len(rooted.children) != 2:
        return None
    first, second = (max(0.0, float(child.branch_length or 0.0)) for child in rooted.children)
    total = first + second
    return frozenset(leaf_names(rooted.children[0])), (first / total if total > 0.0 else 0.5)


def root_on_split(root: TreeNode, side: FrozenSet[str], share: float) -> Optional[TreeNode]:
    """A copy of an unrooted tree (IQ-TREE's trifurcating root) rooted on the edge that
    separates ``side`` from the other leaves, ``side`` as the first root child with
    ``share`` of the edge's length: the root of another tree with the same topology (see
    ``root_split``) put on the same edge at the same relative position. None when no edge
    has that split, or when it is one of the two edges of a bifurcating file root."""
    graph = _Graph(root)
    leaf_set = frozenset(graph.nodes[i].name for i in graph.leaves)
    below = _leaves_below(graph)
    for node_id, parent in graph.parent.items():
        if parent is None or (below[node_id] != side and leaf_set - below[node_id] != side):
            continue
        if graph.parent[parent] is None and len(graph.nodes[parent].children) == 2:
            return None
        length = graph.nodes[node_id].branch_length
        lower_share = share if below[node_id] == side else 1.0 - share
        lower = None if length is None else length * lower_share
        upper = None if length is None else length - lower
        rooted = graph.rooted_on_edge(parent, node_id, upper, lower)
        if below[node_id] == side:
            rooted.children.reverse()
        return rooted
    return None


def order_like(tree: TreeNode, reference: TreeNode) -> TreeNode:
    """Sort every node's children in place by their first leaf in ``reference``'s leaf
    order, so a tree with the same rooted topology lists its leaves in the same order."""
    position = {name: i for i, name in enumerate(leaf_names(reference))}

    def first(node: TreeNode) -> int:
        if node.is_leaf:
            return position.get(node.name, len(position))
        keyed = sorted(((first(child), i, child) for i, child in enumerate(node.children)), key=lambda item: item[:2])
        node.children = [child for _, _, child in keyed]
        return keyed[0][0]

    first(tree)
    return tree


def _leaves_below(graph: _Graph) -> Dict[int, FrozenSet[str]]:
    below: Dict[int, FrozenSet[str]] = {}
    for i in reversed(range(len(graph.nodes))):  # children come after parents in preorder
        node = graph.nodes[i]
        if node.is_leaf:
            below[i] = frozenset([node.name])
        else:
            acc: Set[str] = set()
            for j in graph.adjacency[i]:
                if graph.parent.get(j) == i:
                    acc |= below[j]
            below[i] = frozenset(acc)
    return below


# ------------------------------------------------------------------ splits
def _split_key(side: FrozenSet[str], leaf_set: FrozenSet[str], anchor: str) -> FrozenSet[str]:
    return side if anchor not in side else leaf_set - side


def split_supports(root: TreeNode) -> Dict[FrozenSet[str], Optional[float]]:
    """Non-trivial bipartitions of the unrooted tree, each as the side without the
    alphabetically first leaf, with the support of its edge."""
    graph = _Graph(root)
    leaf_set = frozenset(graph.nodes[i].name for i in graph.leaves)
    if not leaf_set:
        return {}
    anchor = min(leaf_set)
    below = _leaves_below(graph)
    out: Dict[FrozenSet[str], Optional[float]] = {}
    for node_id, parent in graph.parent.items():
        if parent is None:
            continue
        side = below[node_id]
        if len(side) < 2 or len(side) > len(leaf_set) - 2:
            continue
        key = _split_key(side, leaf_set, anchor)
        support = graph.nodes[node_id].support
        if key not in out or out[key] is None:
            out[key] = support
    return out


def copy_split_supports(target: TreeNode, source: TreeNode) -> None:
    """Give each internal edge of ``target`` (in place) the support of the edge with the
    same split in ``source``, None where ``source`` has no such split."""
    supports = split_supports(source)
    graph = _Graph(target)
    leaf_set = frozenset(graph.nodes[i].name for i in graph.leaves)
    anchor = min(leaf_set)
    below = _leaves_below(graph)
    for node_id, parent in graph.parent.items():
        node = graph.nodes[node_id]
        if parent is not None and not node.is_leaf:
            node.support = supports.get(_split_key(below[node_id], leaf_set, anchor))


def split_lengths(root: TreeNode) -> Dict[FrozenSet[str], float]:
    """Every edge of the unrooted tree (leaf edges included) as its split (the side
    without the alphabetically first leaf) with its length; the two edges of a
    bifurcating root are one edge."""
    graph = _Graph(root)
    leaf_set = frozenset(graph.nodes[i].name for i in graph.leaves)
    if not leaf_set:
        return {}
    anchor = min(leaf_set)
    below = _leaves_below(graph)
    out: Dict[FrozenSet[str], float] = {}
    for node_id, parent in graph.parent.items():
        if parent is None or below[node_id] == leaf_set:
            continue
        key = _split_key(below[node_id], leaf_set, anchor)
        out[key] = out.get(key, 0.0) + graph.edge_length(parent, node_id, use_counts=False)
    return out


def robinson_foulds(a: TreeNode, b: TreeNode) -> Tuple[int, float, int]:
    """(RF distance, normalized RF, shared leaves) on the leaves both trees have."""
    shared = set(leaf_names(a)) & set(leaf_names(b))
    n = len(shared)
    if n < 4:
        return 0, 0.0, n
    pa = prune_to(a, shared)
    pb = prune_to(b, shared)
    sa = set(split_supports(pa)) if pa else set()
    sb = set(split_supports(pb)) if pb else set()
    rf = len(sa ^ sb)
    return rf, rf / float(2 * (n - 3)), n


def support_values(root: TreeNode) -> List[float]:
    return [v for v in split_supports(root).values() if v is not None]


__all__ = [
    "NewickParser",
    "TreeNode",
    "copy_split_supports",
    "leaf_names",
    "midpoint_root",
    "order_like",
    "parse_newick",
    "prune_to",
    "read_newick",
    "rename_leaves",
    "robinson_foulds",
    "root_on_outgroup",
    "root_on_split",
    "root_split",
    "split_lengths",
    "split_supports",
    "sum_branch_lengths",
    "support_values",
    "write_newick",
]
