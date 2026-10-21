"""Phylogenetic trees without Qt: the clade grouping shared by every
tree panel (``tree_view``): ``GroupOptimizer`` rotates children so the
fastGEAR clades run contiguously and ``group_root`` reroots the tree on the
edge that keeps them together (every edge tried as the root).

``root_tree`` / ``prepare_tree`` give the drawn tree in one call (the root
asked for, chosen on the full tree, with the file root as the fallback when
there are fewer than two clades, then pruned to the window's genomes);
``tree_view`` places and draws it. ``long_branch_breaks`` is the long branch
rule of the Phylogeny panel's Distance drawing (``tree_view``).
``ladderize`` gives the drawn leaf order (the Order dropdown): larger clades
first at every node, the clade and cluster runs kept together;
``follow_order`` draws a pruned tree in its full tree's order.
``group_starts`` finds where each clade or cluster starts (the nodes the
Phylogeny panel offers to collapse) and ``resolve_collapsed`` finds the
collapsed subtrees, kept as leaf sets, in a drawing of the tree."""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Set, Tuple

# The Newick node type and parser, pruning (ape ``drop.tip``) and midpoint
# rooting live in the Qt free pipeline module, which also roots and writes the
# run's trees; a numeric internal label is TreeNode.support, never a name.
from pipeline.Newick import (
    NewickParser,
    TreeNode,
    midpoint_root,
    prune_to,
    sum_branch_lengths,
)


# ---------------------------------------------------------------- rooting
ROOT_AS_IN_FILE = "file"
ROOT_MIDPOINT = "midpoint"
# Rooted on the edge that keeps the fastGEAR clades (and so their clusters)
# together; the file root when there are fewer than two clades.
ROOT_GROUPS = "groups"
ROOT_MODES = (ROOT_GROUPS, ROOT_AS_IN_FILE, ROOT_MIDPOINT)


def edge_length(node: TreeNode) -> float:
    """A node's edge length; missing, negative or non-finite counts as 0."""
    value = node.branch_length
    if value is None or not math.isfinite(value) or value < 0.0:
        return 0.0
    return float(value)


def analyze_branch_lengths(root: TreeNode) -> Tuple[bool, str]:
    """Whether the tree's branch lengths can be drawn (Distance), and the
    reason shown when they cannot."""
    edge_count = 0
    missing_count = 0
    invalid_count = 0
    positive_count = 0
    stack = [root]
    while stack:
        node = stack.pop()
        for child in node.children:
            edge_count += 1
            branch_length = child.branch_length
            if branch_length is None:
                missing_count += 1
            elif not math.isfinite(branch_length) or branch_length < 0.0:
                invalid_count += 1
            elif branch_length > 0.0:
                positive_count += 1
            stack.append(child)

    if edge_count == 0:
        return False, "Distance unavailable: the tree has no branches."
    issues: List[str] = []
    if missing_count:
        issues.append(f"{missing_count} missing")
    if invalid_count:
        issues.append(f"{invalid_count} negative or non-finite")
    if issues:
        details = " and ".join(issues)
        return (
            False,
            f"Distance unavailable: {details} branch lengths out of {edge_count}; "
            "Topology spacing is shown instead.",
        )
    if positive_count == 0:
        return (
            False,
            "Distance unavailable: all branch lengths are zero; Topology spacing is shown instead.",
        )
    return (
        True,
        "Distance uses cumulative Newick branch-length units.",
    )


# ---------------------------------------------------------------- grouping
@dataclass
class GroupSummary:
    counts: Dict[str, int]
    total: int
    dominant: Optional[str]
    purity: float
    leaf_order: List[str]
    group_order: List[str]
    structural_score: float


def score_group_sequence(groups: Sequence[str]) -> float:
    """How contiguous a sequence of groups is (long runs score high)."""
    if not groups:
        return 0.0
    score = 0.0
    transitions = 0
    missing = 0
    prev: Optional[str] = None
    run_len = 0
    for raw_group in groups:
        group = str(raw_group or "")
        if not group:
            missing += 1
        if prev is None:
            prev = group
            run_len = 1
            continue
        if group == prev:
            run_len += 1
            continue
        if prev:
            score += float(run_len * run_len)
        transitions += 1
        prev = group
        run_len = 1
    if prev:
        score += float(run_len * run_len)
    score -= float(transitions) * 1.35
    score -= float(missing) * 0.75
    return score


class GroupOptimizer:
    """Orders every node's children so leaf groups (fastGEAR clades) run
    contiguously. ``group_of`` maps a leaf name to its group ("" for none);
    ``leaf_name`` maps a raw Newick leaf label to the name ``group_of`` and
    the leaf orders use."""

    def __init__(self, group_of: Callable[[str], str], leaf_name: Callable[[str], str]) -> None:
        self._group_of = group_of
        self._leaf_name = leaf_name

    def optimize(self, node: TreeNode) -> GroupSummary:
        """Rotate ``node``'s subtree in place; returns its summary."""
        if node.is_leaf:
            name = self._leaf_name(node.name or "")
            group = self._group_of(name) or ""
            return GroupSummary(
                counts={group: 1} if group else {},
                total=1,
                dominant=group or None,
                purity=1.0 if group else 0.0,
                leaf_order=[name],
                group_order=[group],
                structural_score=0.0,
            )
        pairs = [(child, self.optimize(child)) for child in node.children]
        pairs = self._optimize_child_pairs(pairs)
        node.children = [child for child, _summary in pairs]
        return self._merge(pairs)

    @staticmethod
    def candidate_score(summary: GroupSummary) -> Tuple[float, int, float, int]:
        distinct = len({group for group in summary.group_order if group})
        missing = sum(1 for group in summary.group_order if not group)
        total = score_group_sequence(summary.group_order) + summary.structural_score
        return total, distinct, summary.purity, -missing

    @staticmethod
    def _child_sort_key(summary: GroupSummary) -> Tuple[int, str, float, int, str]:
        dominant = summary.dominant or ""
        first_leaf = summary.leaf_order[0] if summary.leaf_order else ""
        return (1 if not dominant else 0, dominant, -summary.purity, -summary.total, first_leaf)

    @staticmethod
    def _structural_bonus(summaries: Sequence[GroupSummary]) -> float:
        if len(summaries) <= 1:
            return 0.0
        dominants = [summary.dominant for summary in summaries if summary.dominant]
        distinct = len(set(dominants))
        weighted_purity = sum(summary.purity * float(max(1, summary.total)) for summary in summaries)
        mixed_penalty = sum((1.0 - summary.purity) * float(max(1, summary.total)) for summary in summaries)
        bonus = 0.0
        if distinct > 1:
            bonus += weighted_purity * (1.0 + 0.35 * float(distinct - 1))
        else:
            bonus -= mixed_penalty * 0.6
        if (
            len(summaries) == 2
            and summaries[0].dominant
            and summaries[1].dominant
            and summaries[0].dominant != summaries[1].dominant
        ):
            bonus += 0.5 * float(summaries[0].total + summaries[1].total)
        return bonus

    def _merge(self, pairs: Sequence[Tuple[TreeNode, GroupSummary]]) -> GroupSummary:
        counts: Dict[str, int] = {}
        total = 0
        leaf_order: List[str] = []
        group_order: List[str] = []
        structural = 0.0
        summaries = [summary for _child, summary in pairs]
        for summary in summaries:
            total += summary.total
            leaf_order.extend(summary.leaf_order)
            group_order.extend(summary.group_order)
            structural += summary.structural_score
            for group, count in summary.counts.items():
                counts[group] = counts.get(group, 0) + int(count)
        dominant: Optional[str] = None
        purity = 0.0
        if counts:
            dominant = max(sorted(counts), key=lambda key: (counts[key], key))
            purity = float(counts[dominant]) / float(max(1, total))
        structural += self._structural_bonus(summaries)
        return GroupSummary(counts, total, dominant, purity, leaf_order, group_order, structural)

    def _arrangement_score(self, pairs: Sequence[Tuple[TreeNode, GroupSummary]]) -> Tuple[float, int, float, str]:
        merged = self._merge(pairs)
        total = score_group_sequence(merged.group_order) + merged.structural_score
        distinct = len({summary.dominant for _child, summary in pairs if summary.dominant})
        first_leaf = merged.leaf_order[0] if merged.leaf_order else ""
        return total, distinct, merged.purity, first_leaf

    def _optimize_child_pairs(
        self, pairs: Sequence[Tuple[TreeNode, GroupSummary]]
    ) -> List[Tuple[TreeNode, GroupSummary]]:
        if len(pairs) <= 1:
            return list(pairs)
        if len(pairs) <= 6:
            best_pairs: Optional[List[Tuple[TreeNode, GroupSummary]]] = None
            best_score: Optional[Tuple[float, int, float, str]] = None
            for perm in itertools.permutations(pairs):
                candidate = list(perm)
                score = self._arrangement_score(candidate)
                if best_score is None or score > best_score:
                    best_score = score
                    best_pairs = candidate
            return best_pairs if best_pairs is not None else list(pairs)
        ordered = sorted(pairs, key=lambda pair: self._child_sort_key(pair[1]))
        best: List[Tuple[TreeNode, GroupSummary]] = [ordered[0]]
        for pair in ordered[1:]:
            chosen: Optional[List[Tuple[TreeNode, GroupSummary]]] = None
            chosen_score: Optional[Tuple[float, int, float, str]] = None
            for index in range(len(best) + 1):
                candidate = best[:index] + [pair] + best[index:]
                score = self._arrangement_score(candidate)
                if chosen_score is None or score > chosen_score:
                    chosen_score = score
                    chosen = candidate
            best = chosen if chosen is not None else best
        return best


def clone_tree(node: TreeNode) -> TreeNode:
    """A copy of ``node``'s tree (names, branch lengths and support)."""
    top = TreeNode(name=node.name, branch_length=node.branch_length, support=node.support)
    todo = [(node, top)]
    while todo:
        source, built = todo.pop()
        for child in source.children:
            copy = TreeNode(name=child.name, branch_length=child.branch_length, support=child.support)
            built.children.append(copy)
            todo.append((child, copy))
    return top


def group_root(
    root: TreeNode,
    group_of: Callable[[str], str],
    leaf_name: Callable[[str], str],
    *,
    node_name: Optional[Callable[[TreeNode], str]] = None,
) -> TreeNode:
    """The tree rerooted so the groups run contiguously: the file root and
    then every edge (in file preorder, split in halves) is tried as the root,
    each candidate rotated by ``GroupOptimizer``, and the first strictly best
    one is returned, already rotated. ``node_name`` names the rebuilt nodes
    (the raw label by default). Never mutates ``root``. A node left with one
    child by the reroot (a bifurcating file root) is kept as it is: see
    ``merge_unary``."""
    name_of = node_name or (lambda node: str(node.name or ""))
    optimizer = GroupOptimizer(group_of, leaf_name)
    best = clone_tree(root)
    best_score = GroupOptimizer.candidate_score(optimizer.optimize(best))
    # The unrooted graph: nodes in preorder; a node's neighbors are its
    # parent, then its children in file order; edges in preorder of the child.
    names: List[str] = []
    adjacency: List[List[Tuple[int, float]]] = []
    edges: List[Tuple[int, int, float]] = []
    stack: List[Tuple[TreeNode, int]] = [(root, -1)]
    while stack:
        node, parent = stack.pop()
        index = len(names)
        names.append(name_of(node))
        adjacency.append([])
        if parent >= 0:
            length = edge_length(node)
            adjacency[parent].append((index, length))
            adjacency[index].append((parent, length))
            edges.append((parent, index, length))
        stack.extend((child, index) for child in reversed(node.children))

    def build(start: int, came_from: int, branch: float) -> TreeNode:
        top = TreeNode(name=names[start], branch_length=branch)
        todo: List[Tuple[TreeNode, int, int]] = [(top, start, came_from)]
        while todo:
            built, current, previous = todo.pop()
            for neighbor, length in adjacency[current]:
                if neighbor == previous:
                    continue
                child = TreeNode(name=names[neighbor], branch_length=float(length))
                built.children.append(child)
                todo.append((child, neighbor, current))
        return top

    for left, right, length in edges:
        left_branch = max(0.0, length * 0.5)
        candidate = TreeNode(name="", branch_length=None)
        candidate.children = [build(left, right, left_branch), build(right, left, max(0.0, length - left_branch))]
        score = GroupOptimizer.candidate_score(optimizer.optimize(candidate))
        if score > best_score:
            best, best_score = candidate, score
    return best


def merge_unary(root: TreeNode) -> TreeNode:
    """Merge every internal node with one child into that child, the branch
    lengths summed (in place; returns the new top node)."""
    while len(root.children) == 1:
        only = root.children[0]
        only.branch_length = root.branch_length
        root = only
    stack = [root]
    while stack:
        node = stack.pop()
        merged: List[TreeNode] = []
        for child in node.children:
            while len(child.children) == 1:
                only = child.children[0]
                only.branch_length = sum_branch_lengths(child.branch_length, only.branch_length)
                child = only
            merged.append(child)
        node.children = merged
        stack.extend(merged)
    return root


def leaf_names(root: TreeNode, leaf_name: Callable[[str], str]) -> List[str]:
    """The leaves of ``root`` in drawing order (top to bottom)."""
    names: List[str] = []
    stack = [root]
    while stack:
        node = stack.pop()
        if node.is_leaf:
            names.append(leaf_name(node.name or ""))
        else:
            stack.extend(reversed(node.children))
    return names


def groups_contiguous(root: TreeNode, groups: Mapping[str, str], leaf_name: Callable[[str], str]) -> bool:
    """Every group (leaf name to clade) runs as one block of rows; leaves
    without a group are skipped."""
    seen: Set[str] = set()
    previous = ""
    for name in leaf_names(root, leaf_name):
        group = groups.get(name, "")
        if not group or group == previous:
            continue
        if group in seen:
            return False
        seen.add(group)
        previous = group
    return True


def _has_groups(groups: Optional[Mapping[str, str]]) -> bool:
    return len({value for value in (groups or {}).values() if value}) >= 2


def _rotated(tree: TreeNode, groups: Mapping[str, str], leaf_name: Callable[[str], str]) -> TreeNode:
    """A copy of ``tree`` with its children rotated so the groups run
    contiguously."""
    copy = clone_tree(tree)
    group_map = dict(groups)
    GroupOptimizer(lambda name: group_map.get(name, ""), leaf_name).optimize(copy)
    return copy


def root_tree(
    root: TreeNode,
    *,
    leaf_name: Callable[[str], str],
    root_mode: str = ROOT_AS_IN_FILE,
    groups: Optional[Mapping[str, str]] = None,
    rotate: bool = False,
) -> Tuple[TreeNode, str]:
    """The full tree rooted as ``root_mode`` asks, and the root drawn:
    ROOT_GROUPS falls back to the file root with fewer than two ``groups``
    (``group_root`` itself always returns the best candidate it found, the
    file root included). A group rooted tree comes back already rotated;
    ``rotate`` rotates the other roots too (a copy; with two or more
    groups), so the result is what ``prepare_tree(rotated=True)`` takes."""
    if root_mode == ROOT_GROUPS and _has_groups(groups):
        group_map = dict(groups or {})
        rooted = group_root(root, lambda name: group_map.get(name, ""), leaf_name)
        return merge_unary(rooted), ROOT_GROUPS
    if root_mode == ROOT_MIDPOINT:
        tree, drawn = midpoint_root(root), ROOT_MIDPOINT
    else:
        tree, drawn = root, ROOT_AS_IN_FILE
    if rotate and _has_groups(groups):
        tree = _rotated(tree, groups or {}, leaf_name)
    return tree, drawn


def prepare_tree(
    root: TreeNode,
    *,
    leaf_name: Callable[[str], str],
    root_mode: str = ROOT_AS_IN_FILE,
    keep: Optional[Set[str]] = None,
    groups: Optional[Mapping[str, str]] = None,
    rotated: bool = False,
) -> Optional[TreeNode]:
    """The tree to draw: rooted on the full tree as ``root_mode`` asks
    and rotated so the ``groups`` (leaf name to clade) run contiguously when
    there are two or more (``root_tree``), then pruned to ``keep`` (None
    keeps every leaf). Root and rotation are chosen on the full tree, so a
    filter never moves the root and the pruned rows keep the full order (a
    subsequence of it); the pruned tree is rotated again only when its groups
    are no longer contiguous. ``rotated`` says ``root`` is already a
    ``root_tree(rotate=True)`` result (a cached one). Never mutates
    ``root``; None when no leaf is kept."""
    tree = root
    if not rotated:
        tree, _drawn = root_tree(root, leaf_name=leaf_name, root_mode=root_mode, groups=groups, rotate=True)
    if keep is not None:
        pruned = prune_to(tree, set(keep), leaf_name)
        if pruned is None:
            return None
        if _has_groups(groups) and not groups_contiguous(pruned, groups or {}, leaf_name):
            return _rotated(pruned, groups or {}, leaf_name)
        return pruned
    return clone_tree(tree)


# ---------------------------------------------------------------- order
# The leaf order drawn: the tree ladderized (larger clades first, at the top)
# by clade size or by depth, or as computed (the root and clade rotation
# only, no ladderizing).
ORDER_LADDER_SIZE = "ladder_size"
ORDER_LADDER_DEPTH = "ladder_depth"
ORDER_AS_COMPUTED = "computed"
ORDER_MODES = (ORDER_LADDER_SIZE, ORDER_LADDER_DEPTH, ORDER_AS_COMPUTED)

# A child order and the orientation of each child (0: its own ends as
# computed, 1: mirrored), by child index.
_Config = Tuple[Tuple[int, ...], Tuple[int, ...]]
# Per band map: the first and last band of a subtree's rows (None: no band).
_Ends = Tuple[Tuple[Optional[int], Optional[int]], ...]


def _postorder(root: TreeNode) -> List[TreeNode]:
    """Every node after its children (iterative: deep trees are fine)."""
    out: List[TreeNode] = []
    stack = [root]
    while stack:
        node = stack.pop()
        out.append(node)
        stack.extend(node.children)
    out.reverse()
    return out


def _band_runs(order: Sequence[str], groups: Mapping[str, str]) -> Dict[str, int]:
    """Each leaf's band: the run of its group it sits in along ``order``, so
    a group drawn in two runs today is two bands, each kept together. Leaves
    without a group form runs too, so none can be rotated into a band and
    split it (and a subtree with no group is one band, free to ladderize)."""
    bands: Dict[str, int] = {}
    previous: Optional[str] = None
    band = -1
    for name in order:
        group = groups.get(name, "")
        if band < 0 or group != previous:
            band += 1
            previous = group
        bands[name] = band
    return bands


def _inversions(values: Sequence[int]) -> int:
    """The pairs out of order in ``values`` (merge sort, O(k log k))."""
    items = list(values)
    count = 0
    width = 1
    while width < len(items):
        merged: List[int] = []
        for start in range(0, len(items), 2 * width):
            left = items[start:start + width]
            right = items[start + width:start + 2 * width]
            a = b = 0
            while a < len(left) and b < len(right):
                if left[a] <= right[b]:
                    merged.append(left[a])
                    a += 1
                else:
                    merged.append(right[b])
                    count += len(left) - a
                    b += 1
            merged.extend(left[a:])
            merged.extend(right[b:])
        items = merged
        width *= 2
    return count


def _join_bands(
    parts: Sequence[Tuple[frozenset, Optional[int], Optional[int]]]
) -> Optional[Tuple[Optional[int], Optional[int]]]:
    """The first and last band of children drawn one after the other (each
    child's bands already contiguous), or None when a band would be split."""
    closed: Set[int] = set()
    current: Optional[int] = None
    first: Optional[int] = None
    for bands, head, tail in parts:
        if not bands:
            continue
        if first is None:
            first = head
        if head == current:
            if len(bands) > 1:
                if not bands.isdisjoint(closed):
                    return None
                closed.add(current)
                closed.update(bands - {tail})
        else:
            if current is not None and (current in bands or not bands.isdisjoint(closed)):
                return None
            if current is not None:
                closed.add(current)
            closed.update(bands - {tail})
        current = tail
    return first, current


def ladderize(
    root: TreeNode,
    *,
    leaf_name: Callable[[str], str],
    mode: str = ORDER_LADDER_SIZE,
    keep_together: Sequence[Mapping[str, str]] = (),
) -> TreeNode:
    """A copy of ``root`` ladderized: at every node the children are ordered
    larger first (top), by clade size (number of leaves; then depth) or by
    depth (the longest path to a tip, the child's own edge included; then
    size), ties going to the clade with the smallest leaf name. Rotations
    only: topology, branch lengths, support and the root never change.

    ``keep_together`` (leaf name to group, e.g. fastGEAR lineages and
    clusters) has priority: every run of a group in ``root``'s order (its
    bands; one per group when the groups are contiguous) stays one contiguous
    run, and ladderizing decides the remaining choices. Dynamic programming
    over the tree: each node is drawn with its bands' ends as in ``root`` or
    mirrored, and the order with the fewest out of ladder child pairs is
    chosen (``root``'s own order is always possible, so the bands never
    break). ``ORDER_AS_COMPUTED`` returns a plain copy."""
    if mode not in (ORDER_LADDER_SIZE, ORDER_LADDER_DEPTH):
        return clone_tree(root)
    nodes = _postorder(root)
    order = leaf_names(root, leaf_name)
    band_maps = [_band_runs(order, groups) for groups in keep_together if groups]
    count = len(band_maps)

    size: Dict[int, int] = {}
    reach: Dict[int, float] = {}  # height plus the node's own edge
    smallest: Dict[int, str] = {}
    bands: Dict[int, Tuple[frozenset, ...]] = {}
    ends: Dict[int, _Ends] = {}
    # Per node and orientation: (cost, config) of its best drawing.
    best: Dict[int, List[Tuple[int, _Config]]] = {}

    def mirrored(value: _Ends) -> _Ends:
        return tuple((tail, head) for head, tail in value)

    for node in nodes:
        key = id(node)
        if node.is_leaf:
            name = leaf_name(node.name or "")
            size[key] = 1
            reach[key] = edge_length(node)
            smallest[key] = name
            per_map = [band_map.get(name) for band_map in band_maps]
            bands[key] = tuple(frozenset() if band is None else frozenset((band,)) for band in per_map)
            ends[key] = tuple((band, band) for band in per_map)
            best[key] = [(0, ((), ())), (0, ((), ()))]
            continue
        children = node.children
        ids = [id(child) for child in children]
        size[key] = sum(size[i] for i in ids)
        reach[key] = edge_length(node) + max(reach[i] for i in ids)
        smallest[key] = min(smallest[i] for i in ids)
        bands[key] = tuple(frozenset().union(*(bands[i][m] for i in ids)) for m in range(count))
        # The bands' ends as computed (children in order, none mirrored).
        joined = [
            _join_bands([(bands[i][m], ends[i][m][0], ends[i][m][1]) for i in ids]) for m in range(count)
        ]
        ends[key] = tuple(value if value is not None else (None, None) for value in joined)

        def ladder_key(index: int) -> Tuple[float, ...]:
            i = ids[index]
            depth = round(reach[i], 12)
            primary, secondary = (size[i], depth) if mode == ORDER_LADDER_SIZE else (depth, size[i])
            return (-primary, -secondary)

        preferred = sorted(range(len(ids)), key=lambda index: (ladder_key(index), smallest[ids[index]]))
        rank = {index: position for position, index in enumerate(preferred)}

        # A child whose bands are the same either way round takes its
        # cheaper orientation; only the others are tried both ways.
        choices: List[Tuple[int, ...]] = []
        for i in ids:
            costs = best[i]
            if ends[i] == mirrored(ends[i]):
                choices.append((0 if costs[0][0] <= costs[1][0] else 1,))
            else:
                choices.append((0, 1))
        free = all(len(bands[key][m]) <= 1 for m in range(count))
        if free:
            orientation = tuple(choice[0] for choice in choices)
            cost = sum(best[i][o][0] for i, o in zip(ids, orientation))
            config = (tuple(preferred), orientation)
            best[key] = [(cost, config), (cost, config)]
            continue
        k = len(ids)
        identity = tuple(range(k))
        if k <= 3:
            perms = list(itertools.permutations(identity))
            orientations = list(itertools.product(*choices))
        else:
            # Large polytomy: as computed and mirrored, each with every run
            # of neighbors in one band (one band per map, so movable among
            # themselves) put in ladder order, and the ladder order itself.
            signature = [
                tuple(next(iter(bands[i][m])) for m in range(count))
                if all(len(bands[i][m]) == 1 for m in range(count))
                else None
                for i in ids
            ]
            perms = [identity, tuple(reversed(identity)), tuple(preferred)]
            for base in (identity, tuple(reversed(identity))):
                sorted_runs: List[int] = []
                for _, run in itertools.groupby(
                    base, key=lambda index: signature[index] if signature[index] is not None else ("solo", index)
                ):
                    sorted_runs.extend(sorted(run, key=rank.__getitem__))
                perms.append(tuple(sorted_runs))
            orientations = [
                tuple(choice[0] for choice in choices),
                tuple(choice[-1] for choice in choices),
                tuple(0 for _ in choices),
                tuple(1 for _ in choices),
            ]
        target = ends[key]
        target_mirror = mirrored(target)
        found: List[Optional[Tuple[int, _Config]]] = [None, None]
        # As computed and mirrored first: always valid, and ties keep them.
        candidates: List[_Config] = [(identity, tuple(0 for _ in ids)), (tuple(reversed(identity)), tuple(1 for _ in ids))]
        candidates.extend((perm, orient) for perm in perms for orient in orientations)
        for perm, orient in candidates:
            realized: List[Tuple[Optional[int], Optional[int]]] = []
            for m in range(count):
                parts = []
                for index in perm:
                    i = ids[index]
                    head, tail = ends[i][m]
                    if orient[index]:
                        head, tail = tail, head
                    parts.append((bands[i][m], head, tail))
                value = _join_bands(parts)
                if value is None:
                    break
                realized.append(value)
            else:
                drawn = tuple(realized)
                cost = _inversions([rank[index] for index in perm]) + sum(best[ids[index]][orient[index]][0] for index in range(k))
                for side, wanted in ((0, target), (1, target_mirror)):
                    if drawn == wanted and (found[side] is None or cost < found[side][0]):
                        found[side] = (cost, (perm, orient))
        best[key] = [found[0], found[1]]  # type: ignore[list-item]

    def build() -> TreeNode:
        top_key = id(root)
        side = 0 if best[top_key][0][0] <= best[top_key][1][0] else 1
        top = TreeNode(name=root.name, branch_length=root.branch_length, support=root.support)
        todo: List[Tuple[TreeNode, TreeNode, int]] = [(root, top, side)]
        while todo:
            source, built, orientation = todo.pop()
            if source.is_leaf:
                continue
            perm, orient = best[id(source)][orientation][1]
            for index in perm:
                child = source.children[index]
                copy = TreeNode(name=child.name, branch_length=child.branch_length, support=child.support)
                built.children.append(copy)
                todo.append((child, copy, orient[index]))
        return top

    return build()


def _run_counts(order: Sequence[str], groups: Mapping[str, str]) -> Dict[str, int]:
    """How many runs each group has along ``order`` (a leaf without a group
    ends a run, as it splits the group's strip)."""
    counts: Dict[str, int] = {}
    previous = ""
    for name in order:
        group = groups.get(name, "")
        if group and group != previous:
            counts[group] = counts.get(group, 0) + 1
        previous = group
    return counts


def follow_order(
    root: TreeNode,
    order: Sequence[str],
    *,
    leaf_name: Callable[[str], str],
    keep_together: Sequence[Mapping[str, str]] = (),
) -> Optional[TreeNode]:
    """A copy of ``root`` with every node's children sorted by their first
    leaf in ``order`` (rotations only), e.g. a pruned tree drawn in its full
    tree's ladder order so a filter only removes rows. None when a leaf is
    not in ``order`` or a group of ``keep_together`` would get more runs
    than ``root`` draws it in."""
    rank = {name: index for index, name in enumerate(order)}
    first: Dict[int, int] = {}
    for node in _postorder(root):
        if node.is_leaf:
            position = rank.get(leaf_name(node.name or ""))
            if position is None:
                return None
            first[id(node)] = position
        else:
            first[id(node)] = min(first[id(child)] for child in node.children)
    top = TreeNode(name=root.name, branch_length=root.branch_length, support=root.support)
    todo = [(root, top)]
    while todo:
        source, built = todo.pop()
        for child in sorted(source.children, key=lambda child: first[id(child)]):
            copy = TreeNode(name=child.name, branch_length=child.branch_length, support=child.support)
            built.children.append(copy)
            todo.append((child, copy))
    before = leaf_names(root, leaf_name)
    after = leaf_names(top, leaf_name)
    for groups in keep_together:
        was = _run_counts(before, groups)
        if any(count > was.get(group, 0) for group, count in _run_counts(after, groups).items()):
            return None
    return top


# ---------------------------------------------------------------- layout
@dataclass(frozen=True)
class BranchBreaks:
    """The edges a Distance drawing shortens (``long_branch_breaks``): the
    nodes whose parent edge is broken, their true lengths (same order), the
    reference depth H (the 95th percentile root to tip depth with the broken
    edges counted as 0) and the length every broken edge is drawn with."""

    nodes: Tuple[TreeNode, ...] = ()
    true_lengths: Tuple[float, ...] = ()
    reference_depth: float = 0.0
    drawn_length: float = 0.0

    def __bool__(self) -> bool:
        return bool(self.nodes)


# How long a broken edge is drawn, as a share of the reference depth H.
BROKEN_BRANCH_DRAWN_SHARE = 0.5
# At most this many of the longest edges are broken.
MAX_BROKEN_BRANCHES = 4


def _percentile(values: Sequence[float], share: float) -> float:
    """The ``share`` percentile of ``values`` (linear interpolation)."""
    ordered = sorted(values)
    if not ordered:
        return 0.0
    position = share * float(len(ordered) - 1)
    low = int(math.floor(position))
    high = min(low + 1, len(ordered) - 1)
    return ordered[low] + (ordered[high] - ordered[low]) * (position - low)


def long_branch_breaks(root: TreeNode, factor: float) -> BranchBreaks:
    """The long branch rule of a Distance drawing (dataset agnostic): with the
    edges sorted by length (L1 >= L2 >= ...), B_m the m longest and H_m the
    95th percentile root to tip depth with the edges of B_m counted as 0,
    break B_m for the smallest m (1 to ``MAX_BROKEN_BRANCHES``) with
    L_m > factor * H_m and L_(m+1) <= factor * H_m. A midpoint root on a long
    edge splits it into two long root edges, found as m = 2. Empty when no m
    qualifies (or the tree has no depth left without them)."""
    # Preorder: every node after its parent; edge i is the edge above nodes[i].
    nodes: List[TreeNode] = []
    parent: List[int] = []
    stack: List[Tuple[TreeNode, int]] = [(root, -1)]
    while stack:
        node, up = stack.pop()
        parent.append(up)
        nodes.append(node)
        stack.extend((child, len(nodes) - 1) for child in node.children)
    leaves = [index for index, node in enumerate(nodes) if node.is_leaf]
    if len(leaves) < 2:
        return BranchBreaks()
    lengths = [0.0] + [edge_length(node) for node in nodes[1:]]
    ordered = sorted(range(1, len(nodes)), key=lambda index: (-lengths[index], index))
    for m in range(1, min(MAX_BROKEN_BRANCHES, len(ordered)) + 1):
        broken = set(ordered[:m])
        depth = [0.0] * len(nodes)
        for index in range(1, len(nodes)):
            depth[index] = depth[parent[index]] + (0.0 if index in broken else lengths[index])
        reference = _percentile([depth[leaf] for leaf in leaves], 0.95)
        if reference <= 0.0:
            continue
        limit = float(factor) * reference
        next_length = lengths[ordered[m]] if m < len(ordered) else 0.0
        if lengths[ordered[m - 1]] > limit and next_length <= limit:
            return BranchBreaks(
                nodes=tuple(nodes[index] for index in ordered[:m]),
                true_lengths=tuple(lengths[index] for index in ordered[:m]),
                reference_depth=reference,
                drawn_length=BROKEN_BRANCH_DRAWN_SHARE * reference,
            )
    return BranchBreaks()


# ---------------------------------------------------------------- collapsing
#
# A collapsed subtree is kept as the set of its leaf names (not a node), so it
# survives a new drawing of the same tree (another order, a filter) and is
# found again on the way back from another tree source or root.


_SIGNATURE_MASK = (1 << 64) - 1


def _leaf_signature(names: object) -> Tuple[int, int]:
    """A leaf set's signature: its size and the sum of its names' hashes
    (equal sets have equal signatures; a match is checked on the names)."""
    size = 0
    total = 0
    for name in names:  # type: ignore[attr-defined]
        size += 1
        total = (total + hash(name)) & _SIGNATURE_MASK
    return size, total


def _node_signatures(root: TreeNode, leaf_name: Callable[[str], str]) -> Dict[int, Tuple[int, int]]:
    """Every node's leaf set signature (by node id), in linear time."""
    signatures: Dict[int, Tuple[int, int]] = {}
    for node in _postorder(root):
        if node.is_leaf:
            signatures[id(node)] = (1, hash(leaf_name(node.name or "")) & _SIGNATURE_MASK)
            continue
        size = 0
        total = 0
        for child in node.children:
            child_size, child_total = signatures[id(child)]
            size += child_size
            total = (total + child_total) & _SIGNATURE_MASK
        signatures[id(node)] = (size, total)
    return signatures


class _SubtreeFinder:
    """Finds the node of ``root`` (never the root) whose leaves are exactly
    a given set, without building every node's leaf set."""

    def __init__(self, root: TreeNode, leaf_name: Callable[[str], str]) -> None:
        self.root = root
        self.leaf_name = leaf_name
        self.leaves = frozenset(leaf_names(root, leaf_name))
        signatures = _node_signatures(root, leaf_name)
        # The innermost node of a set (a node with one child holds the same).
        self._by_signature: Dict[Tuple[int, int], List[TreeNode]] = {}
        for node in _postorder(root):
            if node is not root:
                self._by_signature.setdefault(signatures[id(node)], []).append(node)

    def node_of(self, names: frozenset) -> Optional[TreeNode]:
        for node in self._by_signature.get(_leaf_signature(names), ()):
            if frozenset(leaf_names(node, self.leaf_name)) == names:
                return node
        return None


def group_starts(
    root: TreeNode,
    groups: Mapping[str, str],
    leaf_name: Callable[[str], str],
    *,
    min_size: int = 2,
) -> List[Tuple[str, TreeNode]]:
    """Where each group (leaf name to group, such as a fastGEAR lineage or
    cluster) starts in the tree: every maximal subtree whose leaves all have
    that one group, with at least ``min_size`` leaves, in drawing order
    (top to bottom). A monophyletic group starts at one node; a group that
    is not (a cluster with another nested among its genomes) starts at
    several, and no start ever holds a genome of another group. The root
    is never a start (collapsing it would hide the whole tree)."""
    if not groups:
        return []
    group_of: Dict[int, Optional[str]] = {}
    size: Dict[int, int] = {}
    for node in _postorder(root):
        if node.is_leaf:
            group_of[id(node)] = groups.get(leaf_name(node.name or "")) or None
            size[id(node)] = 1
            continue
        child_groups = {group_of[id(child)] for child in node.children}
        only = next(iter(child_groups)) if len(child_groups) == 1 else None
        group_of[id(node)] = only
        size[id(node)] = sum(size[id(child)] for child in node.children)
    starts: List[Tuple[str, TreeNode]] = []
    stack: List[Tuple[TreeNode, Optional[str]]] = [(root, None)]
    while stack:
        node, parent_group = stack.pop()
        group = group_of[id(node)]
        if node is not root and group is not None and group != parent_group:
            if size[id(node)] >= min_size:
                starts.append((group, node))
        stack.extend((child, group) for child in reversed(node.children))
    return starts


def root_grades(
    root: TreeNode,
    groups: Mapping[str, str],
    leaf_name: Callable[[str], str],
) -> List[Tuple[str, frozenset]]:
    """The groups that are a grade at the root: not one subtree, while
    every other leaf is (a root placed inside a clade leaves that clade's
    genomes hanging from the path down to the rest of the tree). Each such
    group with at least two leaves, and its leaves. ``collapse_root_grade``
    draws one as a subtree."""
    if not groups:
        return []
    finder = _SubtreeFinder(root, leaf_name)
    members: Dict[str, Set[str]] = {}
    for name in finder.leaves:
        group = groups.get(name)
        if group:
            members.setdefault(group, set()).add(name)
    out: List[Tuple[str, frozenset]] = []
    for group, names in members.items():
        grade = frozenset(names)
        rest = finder.leaves - grade
        if len(grade) < 2 or not rest or finder.node_of(grade) is not None:
            continue
        if finder.node_of(frozenset(rest)) is not None:
            out.append((group, grade))
    return out


def collapse_root_grade(root: TreeNode, grade: frozenset, leaf_name: Callable[[str], str]) -> Optional[TreeNode]:
    """``root`` with the grade ``grade`` (see ``root_grades``) gathered
    under one new node at the root, so it can be drawn collapsed: the path
    down to the rest of the tree keeps its nodes and lengths (the rest is
    drawn where it was), and each piece of the grade hangs from the new node
    with its distance from the root. The new node comes first when the
    grade's leaves were drawn above the rest. None when ``grade`` is not a
    grade of ``root``. Never mutates ``root``."""
    finder = _SubtreeFinder(root, leaf_name)
    grade = frozenset(grade) & finder.leaves
    rest_node = finder.node_of(frozenset(finder.leaves - grade)) if len(grade) >= 2 else None
    if rest_node is None or finder.node_of(grade) is not None:
        return None
    # The path from the root down to the rest.
    parent_of: Dict[int, TreeNode] = {}
    stack = [root]
    while stack:
        node = stack.pop()
        for child in node.children:
            parent_of[id(child)] = node
            stack.append(child)
    path = [rest_node]
    while path[-1] is not root:
        path.append(parent_of[id(path[-1])])
    path.reverse()
    pieces: List[TreeNode] = []
    built_root = TreeNode(name=root.name, branch_length=root.branch_length, support=root.support)
    built = built_root
    distance = 0.0
    for node, below in zip(path, path[1:]):
        for child in node.children:
            if child is not below:
                piece = clone_tree(child)
                piece.branch_length = distance + edge_length(child)
                pieces.append(piece)
        distance += edge_length(below)
        if below is rest_node:
            next_node = rest_node
        else:
            next_node = TreeNode(name=below.name, branch_length=below.branch_length, support=below.support)
        built.children.append(next_node)
        built = next_node
    holder = TreeNode(name="", branch_length=0.0)
    order = leaf_names(root, leaf_name)
    position = {name: index for index, name in enumerate(order)}
    pieces.sort(key=lambda piece: position.get(leaf_names(piece, leaf_name)[0], 0))
    holder.children = pieces
    first_grade = min(position[name] for name in grade)
    first_rest = min(position[name] for name in finder.leaves - grade)
    if first_grade < first_rest:
        built_root.children.insert(0, holder)
    else:
        built_root.children.append(holder)
    return built_root


def full_subtree(root: TreeNode, names: frozenset, drawn: frozenset, leaf_name: Callable[[str], str]) -> frozenset:
    """The leaves of the smallest subtree of ``root`` (the unpruned tree)
    holding ``names`` when its leaves among ``drawn`` (the pruned tree's)
    are exactly ``names``; else ``names``. A subtree collapsed in a pruned
    drawing is kept so, and found again whatever the filter."""
    names = frozenset(names)
    if not names:
        return names
    count: Dict[int, int] = {}
    size: Dict[int, int] = {}
    best: Optional[TreeNode] = None
    for node in _postorder(root):
        if node.is_leaf:
            count[id(node)] = 1 if leaf_name(node.name or "") in names else 0
            size[id(node)] = 1
        else:
            count[id(node)] = sum(count[id(child)] for child in node.children)
            size[id(node)] = sum(size[id(child)] for child in node.children)
        if best is None and count[id(node)] == len(names):
            best = node
    if best is None:
        return names
    full = frozenset(leaf_names(best, leaf_name))
    return full if (full & drawn) == names else names


def resolve_collapsed(
    root: TreeNode,
    collapsed: Sequence[frozenset],
    leaf_name: Callable[[str], str],
    *,
    pruned: bool = False,
) -> Tuple[Dict[int, frozenset], List[frozenset], Dict[frozenset, frozenset]]:
    """The nodes of ``root`` to draw collapsed, by node id (the outermost
    of nested collapsed subtrees; the inner ones stay collapsed inside), each
    with the collapsed set it draws; the sets that do not match this
    drawing; and each matched set's drawn leaves. A set matches when a node
    holds exactly its drawn leaves and, unless ``pruned``, it has all of
    them. A ``pruned`` drawing (the tree cut to the sample set) collapses
    the node of a set's drawn leaves from two on (pruning keeps a clade a
    clade) and never calls a set with fewer drawn leaves unmatched."""
    finder = _SubtreeFinder(root, leaf_name)
    wanted: Dict[int, frozenset] = {}
    unmatched: List[frozenset] = []
    drawn_of: Dict[frozenset, frozenset] = {}
    for members in collapsed:
        members = frozenset(members)
        drawn = members & finder.leaves
        node = finder.node_of(drawn) if len(drawn) >= 2 else None
        if node is not None and (pruned or drawn == members):
            wanted[id(node)] = members
            drawn_of[members] = drawn
        elif not pruned or len(drawn) >= 2:
            unmatched.append(members)
    # Only the outermost: a collapsed node inside another is not drawn.
    out: Dict[int, frozenset] = {}
    stack: List[TreeNode] = [root]
    while stack:
        node = stack.pop()
        if id(node) in wanted:
            out[id(node)] = wanted[id(node)]
            continue
        stack.extend(node.children)
    return out, unmatched, drawn_of
