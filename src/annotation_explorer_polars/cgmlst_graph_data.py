from __future__ import annotations

import heapq
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .panaroo_common import _file_signature  # noqa: F401 - panaroo_force imports it from here


_CHEWBBACA_ALLELES_FILENAME = "results_alleles.tsv"


def _find_chewbbaca_out(output_dir: Path) -> Optional[Path]:
    direct = output_dir / "chewbbaca" / "out"
    if (direct / _CHEWBBACA_ALLELES_FILENAME).is_file():
        return direct
    root = output_dir / "chewbbaca"
    if not root.exists():
        return None
    candidates = [p.parent for p in root.rglob(_CHEWBBACA_ALLELES_FILENAME) if p.is_file()]
    if not candidates:
        return None
    candidates.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return candidates[0]


@lru_cache(maxsize=4)
def _load_chewbbaca_profiles_cached(
    path_str: str,
    _mtime_ns: int,
    _size: int,
) -> Tuple[Tuple[str, ...], np.ndarray, np.ndarray]:
    path = Path(path_str)
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    if df.empty or df.shape[1] < 2:
        raise ValueError(f"Invalid chewBBACA alleles table: {path}")

    sample_series = df.iloc[:, 0].fillna("").astype(str).map(str.strip)
    keep_mask = sample_series.ne("")
    if not keep_mask.any():
        raise ValueError(f"chewBBACA alleles table contains no sample rows: {path}")

    sample_labels = tuple(sample_series.loc[keep_mask].tolist())
    allele_df = df.loc[keep_mask, df.columns[1:]]
    # chewBBACA writes a novel (inferred) allele as "INF-<id>": it is a called allele,
    # not missing data. LNF, NIPH, PLOT3/5, ASM, ALM, ... stay missing.
    allele_df = allele_df.apply(lambda col: col.str.replace(r"^INF-", "", regex=True))
    numeric_df = allele_df.apply(pd.to_numeric, errors="coerce")
    valid = ~numeric_df.isna().to_numpy(dtype=bool)
    numeric_values = numeric_df.fillna(-1).to_numpy(dtype=np.int32, copy=False)
    numeric_values.setflags(write=False)
    valid.setflags(write=False)
    return sample_labels, numeric_values, valid


@lru_cache(maxsize=4)
def _load_chewbbaca_distance_cache(
    path_str: str,
    _mtime_ns: int,
    _size: int,
) -> Tuple[Tuple[str, ...], np.ndarray, np.ndarray, int]:
    sample_labels, values, valid = _load_chewbbaca_profiles_cached(path_str, _mtime_ns, _size)
    n_samples = int(values.shape[0])
    n_loci = int(values.shape[1]) if values.ndim == 2 else 0
    dist = np.zeros((n_samples, n_samples), dtype=np.int32)
    score = np.zeros((n_samples, n_samples), dtype=np.int64)
    if n_samples <= 1 or n_loci <= 0:
        dist.setflags(write=False)
        score.setflags(write=False)
        return sample_labels, dist, score, n_loci

    locus_span = np.int64(n_loci + 1)
    for i in range(n_samples - 1):
        both = valid[i][np.newaxis, :] & valid[i + 1 :]
        comparable = both.sum(axis=1, dtype=np.int32)
        diff = np.count_nonzero((values[i][np.newaxis, :] != values[i + 1 :]) & both, axis=1).astype(np.int32, copy=False)
        missing = n_loci - comparable
        no_overlap = comparable == 0
        edge_dist = np.where(no_overlap, n_loci + 1, diff).astype(np.int32, copy=False)
        edge_score = edge_dist.astype(np.int64) * locus_span + missing.astype(np.int64)
        dist[i, i + 1 :] = edge_dist
        dist[i + 1 :, i] = edge_dist
        score[i, i + 1 :] = edge_score
        score[i + 1 :, i] = edge_score

    dist.setflags(write=False)
    score.setflags(write=False)
    return sample_labels, dist, score, n_loci


def _build_mst_edges(
    score: np.ndarray,
    display_dist: np.ndarray,
    sample_labels: Sequence[str],
    *,
    reference_sample: Optional[str] = None,
) -> List[Tuple[int, int, int]]:
    """Prim's minimum spanning tree over ``score``: (parent, child, distance)
    edges in the order the tree grows, ties broken by sample label; a parent
    always comes before its children.

    Identical profiles (score 0) are all tied, and Prim hangs them off one
    member as a star. Where a node has ``_IDENTICAL_CHAIN_MIN`` or more such
    leaves, they are linked one after another instead: still a minimum
    spanning tree (every edge between them is 0), and drawn as a compact
    block (``_mst_tree_drawing`` folds long chains) rather than a wide ring."""
    n = int(score.shape[0])
    if n <= 1:
        return []
    labels = [str(label) for label in sample_labels]
    ref = str(reference_sample or "").strip()
    start = 0
    if ref and ref in labels:
        start = labels.index(ref)
    else:
        finite = np.where(score > 0, score, np.iinfo(np.int64).max)
        mean_scores = finite.mean(axis=1)
        start = int(np.argmin(mean_scores))

    # Label order as ranks (equal labels share a rank), so ties compare as
    # the labels do.
    _unique, label_rank = np.unique(np.asarray(labels, dtype=str), return_inverse=True)
    used = np.zeros(n, dtype=bool)
    parent = np.full(n, start, dtype=np.int64)
    used[start] = True
    best_score = np.array(score[start], dtype=np.int64, copy=True)
    best_score[start] = np.iinfo(np.int64).max

    edges: List[Tuple[int, int, int]] = []
    while len(edges) < n - 1:
        candidates = np.flatnonzero(~used)
        if candidates.size == 0:
            break
        # The smallest (score, parent label, label), then the lowest index.
        pick = candidates[best_score[candidates] == best_score[candidates].min()]
        pick = pick[label_rank[parent[pick]] == label_rank[parent[pick]].min()]
        pick = pick[label_rank[pick] == label_rank[pick].min()]
        next_idx = int(pick[0])
        parent_idx = int(parent[next_idx])
        edges.append((parent_idx, next_idx, int(display_dist[parent_idx, next_idx])))
        used[next_idx] = True
        remaining = ~used
        row = np.asarray(score[next_idx], dtype=np.int64)
        better = remaining & (
            (row < best_score)
            | ((row == best_score) & (label_rank[next_idx] < label_rank[parent]))
        )
        best_score[better] = row[better]
        parent[better] = next_idx
    return _chain_identical_leaves(edges, score, display_dist)


_IDENTICAL_CHAIN_MIN = 4


def _chain_identical_leaves(
    edges: List[Tuple[int, int, int]],
    score: np.ndarray,
    display_dist: np.ndarray,
) -> List[Tuple[int, int, int]]:
    """Relink the identical leaves of each node as a chain, in edge order."""
    degree: Dict[int, int] = {}
    for left, right, _distance in edges:
        degree[left] = degree.get(left, 0) + 1
        degree[right] = degree.get(right, 0) + 1
    leaves: Dict[int, List[int]] = {}  # parent -> positions of its identical leaf edges
    for position, (left, right, _distance) in enumerate(edges):
        if degree[right] == 1 and int(score[left, right]) == 0:
            leaves.setdefault(left, []).append(position)
    out = list(edges)
    for positions in leaves.values():
        if len(positions) < _IDENTICAL_CHAIN_MIN:
            continue
        previous = edges[positions[0]][1]
        for position in positions[1:]:
            child = edges[position][1]
            out[position] = (previous, child, int(display_dist[previous, child]))
            previous = child
    return out


# Tree layout of the minimum spanning tree, in units of the drawn node
# diameter (the view scales a node to one unit when zoomed out).
_NODE_CLEARANCE = 0.7  # the smallest node to edge distance; nodes stay twice this apart
_EDGE_MIN_LENGTH = 2.0 * _NODE_CLEARANCE  # identical profiles and the shortest edges
# The edge scale (``_mst_edge_scale``), dataset agnostic: with P the 95th
# percentile of the nonzero edge distances, an edge is drawn to one linear
# scale on which P alleles are _EDGE_P95_LENGTH long (at least
# _EDGE_MIN_LENGTH), unless it is longer than _EDGE_BREAK_FACTOR * P alleles:
# then it is drawn about _EDGE_BROKEN_LENGTH long (longer than every other
# edge's length to scale), with a break mark and its true distance. Broken
# edges grow slowly (_EDGE_BROKEN_GROWTH per doubling) so their order shows.
_EDGE_P95_LENGTH = 16.0
_EDGE_BREAK_FACTOR = 4.0
_EDGE_BROKEN_LENGTH = 2.0 * _EDGE_P95_LENGTH
_EDGE_BROKEN_GROWTH = 0.15
_CHILD_FAN = 1.5 * math.pi  # around a node, the angle its children may use
_FOLD_MIN_NODES = 6  # an unbranched end of the tree this long is folded
_FOLD_ASPECT = 1.5  # width to height of a folded chain
_ROOT_GAP_CHOICES = 12  # up to this many children, the root tries its free room in each gap


def _fold_chain(lengths: Sequence[float]) -> List[Tuple[float, float]]:
    """A chain folded back and forth (boustrophedon) into a block.

    The chain starts at the origin, heading along +x; ``lengths`` are its
    edge lengths. Each edge runs along the row while it stays within the
    block's width, else it turns up to start the next row the other way.
    Rows are an edge length (at least two clearances) apart, so no two edges
    cross and every node keeps clear of the edges it is not on; every node
    stays at x >= 0. Returns, per edge, its (heading, length) in the chain's
    frame: 0, pi (along a row) or pi / 2 (to the next row).
    """
    total = float(sum(lengths))
    if not lengths:
        return []
    width = max(2.0 * _EDGE_MIN_LENGTH, math.sqrt(_FOLD_ASPECT * total * total / len(lengths)))
    x = 0.0
    direction = 1.0
    steps: List[Tuple[float, float]] = []
    for length in lengths:
        ahead = x + direction * length
        if -1e-9 <= ahead <= width + 1e-9:
            x = ahead
            steps.append((0.0 if direction > 0 else math.pi, float(length)))
        else:
            direction = -direction
            steps.append((0.5 * math.pi, float(length)))
    return steps


@dataclass(frozen=True)
class MstEdgeScale:
    """How long the layout draws an edge of a given allele distance.

    Linear (``units_per_allele`` layout units per allele, at least
    ``min_length``) up to ``break_above`` alleles; a longer edge is broken:
    drawn ``broken_length`` long, a little more the longer it is (by
    ``_EDGE_BROKEN_GROWTH`` per doubling of its distance), with a break mark
    and its true distance.
    """

    units_per_allele: float
    break_above: float = math.inf
    broken_length: float = 0.0
    min_length: float = _EDGE_MIN_LENGTH

    def is_broken(self, distance: float) -> bool:
        return float(distance) > self.break_above

    def length(self, distance: float) -> float:
        if self.is_broken(distance):
            return self.broken_length * (1.0 + _EDGE_BROKEN_GROWTH * math.log2(float(distance) / self.break_above))
        return max(self.min_length, self.units_per_allele * max(0.0, float(distance)))

    def broken_edges(self, distances: Sequence[float]) -> Tuple[int, ...]:
        """The positions of the broken edges among ``distances``."""
        return tuple(idx for idx, distance in enumerate(distances) if self.is_broken(distance))


def _mst_edge_scale(distances: Sequence[float]) -> MstEdgeScale:
    """The edge scale of a tree with these edge distances (alleles); see
    ``_EDGE_BREAK_FACTOR``."""
    values = np.asarray([float(d) for d in distances], dtype=float)
    nonzero = values[np.isfinite(values) & (values > 0.0)]
    if nonzero.size == 0:
        return MstEdgeScale(units_per_allele=1.0)
    typical = float(np.percentile(nonzero, 95.0))
    unit = _EDGE_P95_LENGTH / typical
    break_above = _EDGE_BREAK_FACTOR * typical
    if float(nonzero.max()) <= break_above:
        return MstEdgeScale(units_per_allele=unit)
    longest_kept = float(nonzero[nonzero <= break_above].max())
    broken_length = max(_EDGE_BROKEN_LENGTH, 1.1 * unit * longest_kept)
    return MstEdgeScale(units_per_allele=unit, break_above=break_above, broken_length=broken_length)


def _convex_hull(points: np.ndarray) -> np.ndarray:
    pts = np.unique(np.round(np.asarray(points, dtype=float), 9), axis=0)
    if len(pts) <= 2:
        return pts

    def turn(o: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
        return float((a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]))

    lower: List[np.ndarray] = []
    upper: List[np.ndarray] = []
    for point in pts:
        while len(lower) >= 2 and turn(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    for point in pts[::-1]:
        while len(upper) >= 2 and turn(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return np.asarray(lower[:-1] + upper[:-1], dtype=float)


def _polygon_origin_distance(polygon: np.ndarray) -> float:
    """Distance from the origin to a convex polygon (0 when inside it)."""
    if len(polygon) == 1:
        return float(np.hypot(*polygon[0]))
    start = polygon
    side = np.roll(polygon, -1, axis=0) - start
    t = np.clip(-(start * side).sum(axis=1) / np.maximum((side * side).sum(axis=1), 1e-12), 0.0, 1.0)
    closest = start + t[:, np.newaxis] * side
    distance = float(np.min(np.hypot(closest[:, 0], closest[:, 1])))
    if len(polygon) >= 3:
        turns = side[:, 0] * -start[:, 1] - side[:, 1] * -start[:, 0]
        if np.all(turns > 0) or np.all(turns < 0):
            return 0.0
    return distance


def _subtree_wedge(hull: np.ndarray, length: float) -> Optional[Tuple[float, float]]:
    """The angles, seen from a parent, of a child's subtree drawn ``length``
    away along the parent's +x axis, widened so that everything outside the
    wedge stays a clearance away. None while the subtree's hull would come
    within two clearances of the parent (with the parent outside the convex
    hull, the angles never wrap around, even for a subtree that bends back
    beside the parent)."""
    x = hull[:, 0] + length
    y = hull[:, 1]
    if _polygon_origin_distance(np.column_stack([x, y])) < 2.0 * _NODE_CLEARANCE:
        return None
    return _hull_wedge(hull, length)


def _hull_wedge(hull: np.ndarray, length: float) -> Tuple[float, float]:
    """``_subtree_wedge`` for a length already known to be valid."""
    x = hull[:, 0] + length
    y = hull[:, 1]
    radius = np.hypot(x, y)
    angle = np.arctan2(y, x)
    pad = np.arcsin(np.minimum(1.0, _NODE_CLEARANCE / radius))
    return float(np.min(angle - pad)), float(np.max(angle + pad))


def _mst_tree_drawing(
    n: int,
    edges: Sequence[Tuple[int, int, int]],
    sample_labels: Sequence[str],
    scale: Optional[MstEdgeScale] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Node positions for a spanning tree (or forest) that keep the drawing clear.

    The seed and the fallback of ``cgmlst_mst_layout.mst_drawing``. A
    balloon layout, built bottom up from the tree's centroid: each child's
    subtree is drawn in its own local frame, then placed around its parent in
    a wedge of its own. The wedges of siblings never overlap and every wedge
    keeps a clearance band, so by construction no two edges cross, nodes stay
    at least ``2 * _NODE_CLEARANCE`` apart and no node comes within
    ``_NODE_CLEARANCE`` of an edge it is not on. Edges get their length from
    the allele distance (``scale``, by default ``_mst_edge_scale`` of the
    edges: linear, very long edges broken) and only grow where the children
    of a node need room (``_edge_scale_parts`` tells how much). An unbranched end of ``_FOLD_MIN_NODES`` or
    more nodes is folded back and forth into a block (``_fold_chain``) rather
    than drawn as one long line. Deterministic: ties go by sample label.
    The tree's long axis ends up horizontal (the view may turn it to fit).
    Returns the positions and, per node, its parent in the layout (-1 for
    the root of each component).
    """
    count = max(0, int(n))
    if count <= 1:
        return np.zeros((count, 2), dtype=float), np.full(count, -1, dtype=np.int64)
    if scale is None:
        scale = _mst_edge_scale([distance for _left, _right, distance in edges])
    labels = [str(label) for label in sample_labels]
    rank = np.empty(count, dtype=np.int64)
    rank[sorted(range(count), key=lambda idx: (labels[idx] if idx < len(labels) else "", idx))] = np.arange(count)
    adjacency: List[List[Tuple[int, float]]] = [[] for _ in range(count)]
    for left, right, distance in edges:
        if 0 <= int(left) < count and 0 <= int(right) < count and int(left) != int(right):
            adjacency[int(left)].append((int(right), float(distance)))
            adjacency[int(right)].append((int(left), float(distance)))
    for neighbors in adjacency:
        neighbors.sort(key=lambda item: rank[item[0]])

    def walk(root: int, parent: np.ndarray) -> Tuple[List[int], Dict[int, float]]:
        parent[root] = root
        order = [root]
        lengths: Dict[int, float] = {}
        for node in order:
            for child, distance in adjacency[node]:
                if parent[child] == -1:
                    parent[child] = node
                    lengths[child] = scale.length(distance)
                    order.append(child)
        return order, lengths

    def component_layout(members: List[int], first_parent: np.ndarray) -> Tuple[Dict[int, np.ndarray], np.ndarray]:
        size = {node: 1 for node in members}
        for node in reversed(members[1:]):
            size[int(first_parent[node])] += size[node]
        total = len(members)

        def heaviest_part(node: int) -> int:
            parts = [total - size[node]]
            parts.extend(size[child] for child, _ in adjacency[node] if int(first_parent[child]) == node and child != node)
            return max(parts)

        root = min(members, key=lambda node: (heaviest_part(node), rank[node]))
        # The walk only follows edges, so it never leaves this component.
        parent = np.full(count, -1, dtype=np.int64)
        order, lengths = walk(root, parent)
        children: Dict[int, List[int]] = {node: [] for node in order}
        for node in order[1:]:
            children[int(parent[node])].append(node)
        # Unbranched ends: the number of nodes from a node down to its leaf
        # (0 where the tree branches below). Long ones are folded into a block.
        chain: Dict[int, int] = {}
        for node in reversed(order):
            kids = children[node]
            chain[node] = 1 if not kids else (chain[kids[0]] + 1 if len(kids) == 1 and chain[kids[0]] else 0)
        folds: Dict[int, List[int]] = {}  # the top of a folded chain -> the chain
        for node in order:
            if chain[node] >= _FOLD_MIN_NODES and (node == root or not chain[int(parent[node])]):
                path = [node]
                while children[path[-1]]:
                    path.append(children[path[-1]][0])
                folds[node] = path
        folded = {node for path in folds.values() for node in path[1:]}
        hulls: Dict[int, np.ndarray] = {}
        local: Dict[int, Tuple[float, float]] = {}  # child -> (angle around its parent, length)
        for node in reversed(order):
            if node in folded:
                continue
            kids = children[node]
            if not kids:
                hulls[node] = np.zeros((1, 2), dtype=float)
                continue
            path = folds.get(node)
            if path is not None:
                points = [(0.0, 0.0)]
                heading = 0.0
                x = y = 0.0
                for child, (direction, length) in zip(path[1:], _fold_chain([lengths[c] for c in path[1:]])):
                    local[child] = (direction - heading, length)
                    heading = direction
                    x += length * math.cos(direction)
                    y += length * math.sin(direction)
                    points.append((x, y))
                hulls[node] = _convex_hull(np.asarray(points, dtype=float))
                continue
            room = 2.0 * math.pi if node == root else _CHILD_FAN
            reach: List[float] = []
            wedges: List[Tuple[float, float]] = []
            for child in kids:
                length = lengths[child]
                wedge = _subtree_wedge(hulls[child], length)
                grow = 0.25
                while wedge is None:
                    length += grow
                    grow *= 1.5
                    wedge = _subtree_wedge(hulls[child], length)
                reach.append(length)
                wedges.append(wedge)
            used = sum(hi - lo for lo, hi in wedges)
            if used > room:
                # Too many children for the room: lengthen, a step at a
                # time, the edge whose wedge narrows the most per unit of
                # length. Leaves all narrow alike, so they are one entry:
                # the shortest leaf edges grow together to a common length.
                leaf = np.array([len(hulls[child]) == 1 and not hulls[child].any() for child in kids], dtype=bool)
                leaf_base = np.asarray(reach, dtype=float)[leaf]
                floor = float(leaf_base.min()) if leaf_base.size else 0.0

                def leaf_width(length: float) -> float:
                    return float(np.sum(2.0 * np.arcsin(_NODE_CLEARANCE / np.maximum(leaf_base, length))))

                def lengthened(idx: int) -> Tuple[float, float, Tuple[float, float]]:
                    step = max(0.2, 0.08 * reach[idx])
                    wedge = _hull_wedge(hulls[kids[idx]], reach[idx] + step)  # longer only narrows a valid wedge
                    gain = ((wedges[idx][1] - wedges[idx][0]) - (wedge[1] - wedge[0])) / step
                    return gain, step, wedge

                def leaves_lengthened() -> Tuple[float, float]:
                    step = max(0.2, 0.08 * floor)
                    single = 2.0 * (math.asin(_NODE_CLEARANCE / floor) - math.asin(_NODE_CLEARANCE / (floor + step)))
                    return single / step, step

                queue: List[Tuple[float, int, int, float, Optional[Tuple[float, float]]]] = []
                for idx in np.flatnonzero(~leaf).tolist():
                    gain, step, wedge = lengthened(idx)
                    queue.append((-gain, int(rank[kids[idx]]), idx, step, wedge))
                if leaf_base.size:
                    gain, step = leaves_lengthened()
                    queue.append((-gain, -1, -1, step, None))
                heapq.heapify(queue)
                while used > room:
                    _gain, order_rank, idx, step, wedge = heapq.heappop(queue)
                    if idx < 0:
                        before = leaf_width(floor)
                        floor += step
                        used += leaf_width(floor) - before
                        gain, step = leaves_lengthened()
                        heapq.heappush(queue, (-gain, -1, -1, step, None))
                        continue
                    assert wedge is not None
                    used += (wedge[1] - wedge[0]) - (wedges[idx][1] - wedges[idx][0])
                    reach[idx] += step
                    wedges[idx] = wedge
                    gain, step, wedge = lengthened(idx)
                    heapq.heappush(queue, (-gain, order_rank, idx, step, wedge))
                for idx in np.flatnonzero(leaf).tolist():
                    if reach[idx] < floor:
                        reach[idx] = floor
                        pad = math.asin(_NODE_CLEARANCE / floor)
                        wedges[idx] = (-pad, pad)
            # The widest subtree straight ahead, the others alternately on
            # either side; around the root, evenly spread.
            by_width = sorted(range(len(kids)), key=lambda idx: (-(wedges[idx][1] - wedges[idx][0]), rank[kids[idx]]))
            sequence: List[int] = []
            for position, idx in enumerate(by_width):
                if position % 2 == 0:
                    sequence.append(idx)
                else:
                    sequence.insert(0, idx)
            # Around a node, the children's wedges side by side, centered on
            # the way ahead. Around the root, the room left over goes between
            # them: evenly, or all of it in one place when that gives a more
            # compact picture (two big halves at an angle rather than in a
            # long line).
            gaps = [[0.0] * len(kids)]
            if node == root:
                free = max(0.0, room - used)
                gaps = [[free / len(kids)] * len(kids)]
                if 1 < len(kids) <= _ROOT_GAP_CHOICES and free > 0.05:
                    gaps.extend([free if slot == place else 0.0 for slot in range(len(kids))] for place in range(len(kids)))
            best: Optional[Tuple[float, List[Tuple[int, float]], List[np.ndarray]]] = None
            for spacing in gaps:
                cursor = 0.0 if node == root else -0.5 * used
                angles: List[Tuple[int, float]] = []
                points = [np.zeros((1, 2), dtype=float)]
                for idx, gap in zip(sequence, spacing):
                    lo, hi = wedges[idx]
                    angle = cursor - lo
                    cursor += (hi - lo) + gap
                    angles.append((idx, angle))
                    cos_a, sin_a = math.cos(angle), math.sin(angle)
                    shifted = hulls[kids[idx]] + np.array([reach[idx], 0.0])
                    points.append(shifted @ np.array([[cos_a, sin_a], [-sin_a, cos_a]]))
                if len(gaps) == 1:
                    best = (0.0, angles, points)
                    break
                side = _square_side(np.vstack(points))
                if best is None or side < 0.97 * best[0]:
                    best = (side, angles, points)
            assert best is not None
            _side, angles, points = best
            for idx, angle in angles:
                child = kids[idx]
                local[child] = (angle, reach[idx])
                del hulls[child]
            hulls[node] = _convex_hull(np.vstack(points))
        heading = {root: 0.0}
        placed = {root: np.zeros(2, dtype=float)}
        for node in order[1:]:
            angle, length = local[node]
            up = int(parent[node])
            heading[node] = heading[up] + angle
            placed[node] = placed[up] + length * np.array([math.cos(heading[node]), math.sin(heading[node])])
        return placed, parent

    coords = np.zeros((count, 2), dtype=float)
    parents = np.full(count, -1, dtype=np.int64)
    first_parent = np.full(count, -1, dtype=np.int64)
    offset_x = 0.0
    for start in sorted(range(count), key=lambda node: rank[node]):
        if first_parent[start] != -1:
            continue
        members, _lengths = walk(start, first_parent)
        placed, parent = component_layout(members, first_parent)
        parents[members] = parent[members]
        part = np.array([placed[node] for node in members], dtype=float)
        part = _long_axis_horizontal(part)
        part[:, 0] += offset_x - float(np.min(part[:, 0]))
        part[:, 1] -= 0.5 * (float(np.min(part[:, 1])) + float(np.max(part[:, 1])))
        offset_x = float(np.max(part[:, 0])) + 4.0 * _NODE_CLEARANCE
        coords[members] = part
    coords -= 0.5 * (coords.min(axis=0) + coords.max(axis=0))
    roots = parents == np.arange(count)
    parents[roots] = -1
    return coords, parents


def _edge_scale_parts(
    coords: np.ndarray,
    edges: Sequence[Tuple[int, int, int]],
    parents: np.ndarray,
    scale: MstEdgeScale,
) -> Tuple[np.ndarray, np.ndarray]:
    """Where each edge is drawn to scale: the layout lengthens some edges to
    give a node's children room, and the view draws only the first
    ``scale.length`` of such an edge solid, from the parent's end, and the
    rest dotted. Returns, per edge, the node its solid part starts from and
    the solid share of the drawn edge (1.0: drawn to scale)."""
    count = len(edges)
    start = np.zeros(count, dtype=np.int64)
    share = np.ones(count, dtype=float)
    for idx, (left, right, distance) in enumerate(edges):
        left, right = int(left), int(right)
        start[idx] = right if 0 <= right < len(parents) and int(parents[left]) == right else left
        drawn = float(np.hypot(*(coords[right] - coords[left]))) if max(left, right) < len(coords) else 0.0
        wanted = scale.length(distance)
        if drawn > 1.02 * wanted + 1e-9:
            share[idx] = wanted / drawn
    return start, share


def _square_side(points: np.ndarray) -> float:
    """The side of the smallest square, at any of 36 turns, around the points
    (the picture's size when fitted into a square-ish plot)."""
    hull = _convex_hull(points)
    angles = np.radians(np.arange(0.0, 90.0, 2.5))
    x = np.outer(np.cos(angles), hull[:, 0]) + np.outer(np.sin(angles), hull[:, 1])
    y = np.outer(-np.sin(angles), hull[:, 0]) + np.outer(np.cos(angles), hull[:, 1])
    return float(np.min(np.maximum(np.ptp(x, axis=1), np.ptp(y, axis=1))))


def _long_axis_horizontal(coords: np.ndarray) -> np.ndarray:
    """Rotate about the centroid so the main axis of the points is horizontal."""
    centered = coords - coords.mean(axis=0)
    if len(centered) < 2:
        return centered
    _values, vectors = np.linalg.eigh(np.cov(centered.T))
    major = vectors[:, 1]
    if major[0] < 0 or (abs(major[0]) < 1e-12 and major[1] < 0):
        major = -major
    rotation = np.array([[major[0], major[1]], [-major[1], major[0]]])
    return centered @ rotation.T
