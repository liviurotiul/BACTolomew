"""The drawing of the cgMLST minimum spanning tree.

``mst_drawing`` starts from the balloon layout (``_mst_tree_drawing``: clear
by construction, but mechanical: even fans, straight runs, shared
directions) and relaxes it into a natural picture that keeps it clear:

* up to ``_SPRINGS_MAX_NODES`` nodes, constrained springs (``_spring_layout``):
  repulsion between nodes and from foreign edges turns each edge about its
  parent, and lengthened edges shrink back toward their length to scale;
* larger trees, up to ``_STRESS_MAX_NODES``, stress majorization
  (``_stress_layout``): faster, as natural in its angles, a few more dotted
  edges;
* above that, for long fans that would cost seconds and gain nothing (the
  cost guards in ``mst_drawing``), or when the result fails the final
  check, the balloon itself.

Both relaxations work on the edges' headings and lengths, never draw an
edge shorter than its length to scale (``MstEdgeScale.length``), never move
a node more than a fraction of the clearance per step and take back the
part of a step that breaks a rule, so every accepted state keeps the
balloon's guarantees: nodes ``2 * _NODE_CLEARANCE`` apart, no node within
``_NODE_CLEARANCE`` of an edge it is not on, no crossings. ``_clarity_problems``
checks the result in full anyway. Deterministic: no random numbers, a fixed
order everywhere. Units are node diameters, as in ``cgmlst_graph_data``.
"""
from __future__ import annotations

import math
from typing import List, Optional, Sequence, Tuple

import numpy as np
from scipy.spatial import cKDTree

from .cgmlst_graph_data import (
    _NODE_CLEARANCE,
    MstEdgeScale,
    _long_axis_horizontal,
    _mst_edge_scale,
    _mst_tree_drawing,
)

_NODE_MIN = 2.0 * _NODE_CLEARANCE
_SPRINGS_MAX_NODES = 1000
_STRESS_MAX_NODES = 2500  # the stress targets are an n x n matrix
# Cost guards (about 3.5 s at most): the balloon's edges' total length, in
# node diameters (a 1000 node cohort draws about 7000), and for stress the
# most edges at one node.
_SPRINGS_MAX_LENGTH = 60_000
_STRESS_MAX_LENGTH = 20_000
_STRESS_MAX_FAN = 150

Edges = Sequence[Tuple[int, int, int]]


def mst_drawing(
    n: int,
    edges: Edges,
    sample_labels: Sequence[str],
    scale: Optional[MstEdgeScale] = None,
) -> Tuple[np.ndarray, np.ndarray, str]:
    """Node positions for the minimum spanning tree, the layout's parents
    (as ``_mst_tree_drawing``) and a line saying which layout was used and,
    for a fallback, why. The long axis ends up horizontal."""
    if scale is None:
        scale = _mst_edge_scale([distance for _left, _right, distance in edges])
    seed, parents = _mst_tree_drawing(n, edges, sample_labels, scale)
    count = int(seed.shape[0])
    if count < 3:
        return seed, parents, "balloon (fewer than 3 nodes)"
    if count > _STRESS_MAX_NODES:
        return seed, parents, f"balloon ({count} nodes, above {_STRESS_MAX_NODES})"
    method = "springs" if count <= _SPRINGS_MAX_NODES else "stress"
    # Both relaxations sample every edge at each step, and the stress check
    # compares every pair of edges at a node: long fans (a hub with many
    # isolates, all far from each other) cost seconds and gain nothing, as
    # their edges cannot get any shorter.
    nr = np.flatnonzero(parents >= 0)
    drawn = float(np.linalg.norm(seed[nr] - seed[parents[nr]], axis=1).sum())
    max_length = _SPRINGS_MAX_LENGTH if method == "springs" else _STRESS_MAX_LENGTH
    if drawn > max_length:
        return seed, parents, f"balloon (edges {drawn:.0f} node diameters long in all, above {max_length} for {method})"
    fan = int(np.bincount(parents[nr], minlength=count).max()) + 1
    if method == "stress" and fan > _STRESS_MAX_FAN:
        return seed, parents, f"balloon (a node with {fan} edges, above {_STRESS_MAX_FAN} for stress)"
    try:
        coords = (_spring_layout if method == "springs" else _stress_layout)(seed, parents, edges, scale)
    except Exception as exc:  # pragma: no cover - a safety net, logged
        return seed, parents, f"balloon ({method} layout failed: {type(exc).__name__}: {exc})"
    problems = _clarity_problems(coords, parents)
    if problems:
        return seed, parents, f"balloon ({method} layout failed the clarity check: {problems})"
    coords = _long_axis_horizontal(coords)
    coords -= 0.5 * (coords.min(axis=0) + coords.max(axis=0))
    return coords, parents, method


def _clarity_problems(coords: np.ndarray, parents: np.ndarray) -> str:
    """'' when the drawing keeps every rule (nodes 2 clearances apart, no
    node within a clearance of an edge it is not on, no crossing), else what
    breaks."""
    if not np.all(np.isfinite(coords)):
        return "positions not finite"
    tree = _Tree(len(coords), parents)
    overlaps, on_edges, crossings = _Geometry(tree).problems(
        coords, _NODE_MIN - 1e-6, _NODE_CLEARANCE - 1e-6, 0.5
    )
    found = [
        f"{count} {what}"
        for count, what in ((overlaps, "node overlaps"), (on_edges, "nodes on edges"), (crossings, "crossings"))
        if count
    ]
    return ", ".join(found)


# ---------------------------------------------------------------- the tree

class _Tree:
    """A forest given by its parents (-1 at a root), in DFS preorder: a
    subtree is the range tin[v]:tout[v] of ``order``."""

    def __init__(self, n: int, parents: np.ndarray):
        par = np.asarray(parents, dtype=np.int64)
        children: List[List[int]] = [[] for _ in range(n)]
        for v in range(n):
            if par[v] >= 0:
                children[int(par[v])].append(v)
        roots = [v for v in range(n) if par[v] < 0]
        order: List[int] = []
        tin = np.zeros(n, np.int64)
        depth = np.zeros(n, np.int64)
        for root in roots:
            stack = [root]
            while stack:
                v = stack.pop()
                tin[v] = len(order)
                order.append(v)
                for c in reversed(children[v]):
                    depth[c] = depth[v] + 1
                    stack.append(c)
        size = np.ones(n, np.int64)
        for v in order[::-1]:
            if par[v] >= 0:
                size[par[v]] += size[v]
        self.n = n
        self.par = par
        self.order = np.asarray(order, np.int64)
        self.tin = tin
        self.tout = tin + size
        self.size = size
        self.depth = depth
        self.roots = np.asarray(roots, np.int64)
        self.nonroot = np.flatnonzero(par >= 0)

    def positions(self, theta: np.ndarray, length: np.ndarray, root_pos: np.ndarray) -> np.ndarray:
        """Node positions from each edge's heading and length (prefix sums
        over the preorder)."""
        n = self.n
        steps = np.zeros((n + 1, 2))
        v = self.nonroot
        vec = np.column_stack([length[v] * np.cos(theta[v]), length[v] * np.sin(theta[v])])
        np.add.at(steps, self.tin[v], vec)
        np.add.at(steps, self.tout[v], -vec)
        r = self.roots
        np.add.at(steps, self.tin[r], root_pos[r])
        np.add.at(steps, self.tout[r], -root_pos[r])
        return np.cumsum(steps[:n], axis=0)[self.tin]

    def subtree_sum(self, values: np.ndarray) -> np.ndarray:
        """Sum of ``values`` (per node, any trailing shape) over each subtree."""
        pre = values[self.order]
        cs = np.concatenate([np.zeros((1,) + pre.shape[1:]), np.cumsum(pre, axis=0)])
        return cs[self.tout] - cs[self.tin]

    def subtree_reach(self, P: np.ndarray) -> np.ndarray:
        """Per node, about the farthest its subtree reaches from its parent
        (the bounding box corner farthest from it)."""
        n = self.n
        pre = P[self.order]
        lo_t, hi_t = [pre], [pre]
        k = 1
        while (1 << k) <= n:
            h = 1 << (k - 1)
            lo, hi = lo_t[-1], hi_t[-1]
            lo_t.append(np.minimum(lo[: len(lo) - h], lo[h:]))
            hi_t.append(np.maximum(hi[: len(hi) - h], hi[h:]))
            k += 1
        level = np.floor(np.log2(self.size)).astype(np.int64)
        lo_box = np.empty((n, 2))
        hi_box = np.empty((n, 2))
        for k in np.unique(level).tolist():
            sel = np.flatnonzero(level == k)
            a = self.tin[sel]
            b = self.tout[sel] - (1 << k)
            lo_box[sel] = np.minimum(lo_t[k][a], lo_t[k][b])
            hi_box[sel] = np.maximum(hi_t[k][a], hi_t[k][b])
        pp = P[np.maximum(self.par, 0)]
        return np.hypot(
            np.maximum(np.abs(lo_box[:, 0] - pp[:, 0]), np.abs(hi_box[:, 0] - pp[:, 0])),
            np.maximum(np.abs(lo_box[:, 1] - pp[:, 1]), np.abs(hi_box[:, 1] - pp[:, 1])),
        )

    def path_nodes(self, a: int, b: int) -> List[int]:
        """Nodes on the paths from lca(a, b) down to a and b (lca excluded)."""
        out = []
        par, depth = self.par, self.depth
        while depth[a] > depth[b]:
            out.append(a)
            a = int(par[a])
        while depth[b] > depth[a]:
            out.append(b)
            b = int(par[b])
        while a != b:
            out.append(a)
            out.append(b)
            a, b = int(par[a]), int(par[b])
        return out


def _segment_samples(P, el, er, step):
    """Points along each edge at most ``step`` apart, and their edge."""
    L = np.linalg.norm(P[er] - P[el], axis=1)
    cnt = np.maximum(2, np.ceil(L / step).astype(np.int64) + 1)
    eid = np.repeat(np.arange(len(el)), cnt)
    start = np.concatenate([[0], np.cumsum(cnt)[:-1]])
    t = (np.arange(eid.size) - np.repeat(start, cnt)) / np.repeat(cnt - 1, cnt)
    return P[el[eid]] + t[:, None] * (P[er[eid]] - P[el[eid]]), eid


def _point_segment(P, el, er, nodes, eids):
    """Distance of each node to its edge, the offset from the closest point
    and that point's place along the edge (0 to 1)."""
    a = P[el[eids]]
    ab = P[er[eids]] - a
    t = np.clip(((P[nodes] - a) * ab).sum(1) / np.maximum((ab * ab).sum(1), 1e-12), 0.0, 1.0)
    diff = P[nodes] - (a + t[:, None] * ab)
    return np.linalg.norm(diff, axis=1), diff, t


def _node_edge_pairs(tree, stree, eid, n_edges, r, el, er):
    """Unique (node, edge) pairs with a sample of the edge within r of the
    node, the node not on the edge; sorted."""
    m = tree.sparse_distance_matrix(stree, r, output_type="ndarray")
    if m.size == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    key = np.unique(m["i"].astype(np.int64) * n_edges + eid[m["j"]])
    nodes, eids = key // n_edges, key % n_edges
    keep = (el[eids] != nodes) & (er[eids] != nodes)
    return nodes[keep], eids[keep]


def _orient(a, b, c):
    return (b[..., 0] - a[..., 0]) * (c[..., 1] - a[..., 1]) - (b[..., 1] - a[..., 1]) * (c[..., 0] - a[..., 0])


def _cross(p1, p2, q1, q2):
    return (_orient(p1, p2, q1) * _orient(p1, p2, q2) < 0) & (_orient(q1, q2, p1) * _orient(q1, q2, p2) < 0)


def _seg_dist(p, a, b):
    """Distance of points p (..., 2) to segments a-b (..., 2), broadcast."""
    ab = b - a
    t = np.clip(((p - a) * ab).sum(-1) / np.maximum((ab * ab).sum(-1), 1e-12), 0.0, 1.0)
    return np.linalg.norm(p - (a + t[..., None] * ab), axis=-1)


class _Geometry:
    """The tree's edges, edge i joining par(child[i]) and child[i]."""

    def __init__(self, tree: _Tree):
        self.child = tree.nonroot
        self.el = tree.par[self.child]
        self.er = self.child

    def _search(self, P, node_min, clear, step):
        kd = cKDTree(P)
        close = kd.query_pairs(node_min, output_type="ndarray")
        pts, eid = _segment_samples(P, self.el, self.er, step)
        skd = cKDTree(pts)
        # A node within ``clear`` of an edge is within this of a sample.
        r = math.sqrt(clear * clear + 0.25 * step * step) + 1e-9
        nodes, eids = _node_edge_pairs(kd, skd, eid, len(self.el), r, self.el, self.er)
        if nodes.size:
            d, _diff, _t = _point_segment(P, self.el, self.er, nodes, eids)
            nodes, eids = nodes[d < clear], eids[d < clear]
        return close, nodes, eids, skd, eid

    def _crossings(self, P, skd, eid, step):
        """Pairs of crossing edges: two crossing edges have a sample each
        within step / 2 of the crossing."""
        el, er = self.el, self.er
        sp = skd.query_pairs(step + 1e-9, output_type="ndarray")
        if not sp.size:
            return np.zeros(0, np.int64), np.zeros(0, np.int64)
        x, y = eid[sp[:, 0]], eid[sp[:, 1]]
        other = x != y
        key = np.unique(np.minimum(x[other], y[other]) * len(el) + np.maximum(x[other], y[other]))
        i, j = key // len(el), key % len(el)
        share = (el[i] == el[j]) | (el[i] == er[j]) | (er[i] == el[j]) | (er[i] == er[j])
        i, j = i[~share], j[~share]
        cross = _cross(P[el[i]], P[er[i]], P[el[j]], P[er[j]])
        return i[cross], j[cross]

    def violations(self, P, node_min, clear, step, crossings=True):
        """(a, b) node pairs whose relative placement breaks a rule."""
        if len(self.el) == 0:
            return [tuple(pair) for pair in cKDTree(P).query_pairs(node_min, output_type="ndarray").tolist()]
        close, nodes, eids, skd, eid = self._search(P, node_min, clear, step)
        out = [tuple(pair) for pair in close.tolist()]
        for u, e in zip(nodes.tolist(), eids.tolist()):
            out.append((u, int(self.child[e])))
            out.append((u, int(self.el[e])))
        if crossings:
            i, j = self._crossings(P, skd, eid, step)
            for a, b in zip(i.tolist(), j.tolist()):
                out.append((int(self.child[a]), int(self.child[b])))
                out.append((int(self.el[a]), int(self.el[b])))
        return out

    def problems(self, P, node_min, clear, step) -> Tuple[int, int, int]:
        """Counts of node overlaps, (node, foreign edge) pairs and crossings."""
        if len(self.el) == 0:
            return len(cKDTree(P).query_pairs(node_min)), 0, 0
        close, nodes, _eids, skd, eid = self._search(P, node_min, clear, step)
        i, _j = self._crossings(P, skd, eid, step)
        return len(close), int(nodes.size), int(i.size)


def _edge_state(seed, parents, edges, scale, tree):
    """Per node, the length to scale of the edge to its parent, and the
    edge's heading and length in ``seed``."""
    n = tree.n
    distance = {}
    for a, b, d in edges:
        distance[(int(a), int(b))] = distance[(int(b), int(a))] = d
    target = np.zeros(n)
    nr = tree.nonroot
    target[nr] = [scale.length(distance[(int(parents[v]), int(v))]) for v in nr.tolist()]
    vec = seed[nr] - seed[tree.par[nr]]
    theta = np.zeros(n)
    length = np.zeros(n)
    theta[nr] = np.arctan2(vec[:, 1], vec[:, 0])
    length[nr] = np.hypot(vec[:, 0], vec[:, 1])
    return target, theta, length


# ---------------------------------------------------------------- springs

# Forces (layout units): Fruchterman-Reingold repulsion between nodes closer
# than _CUTOFF (tapered to 0 there) plus a stiff barrier below _BARRIER, a
# softened all-pairs repulsion (_K_FAR), a push of nodes off foreign edges
# closer than _EDGE_REACH and a weak pull toward the root.
_K = 2.2
_CUTOFF = 7.0
_BARRIER = 2.2
_EDGE_REACH = 1.6
_GRAVITY = 0.0005
_K_FAR = 1.0
_SOFT = 2.0
_EXACT_MAX = 400  # above this many nodes, the far field uses a 32 x 32 grid
_FAR_EVERY = 4  # steps between far field updates
# Steps: a turn is _ETA * torque / inertia, at most _MAX_TURN rad and a
# _MAX_MOVE shift of the subtree at full temperature, cooling to _T_END; a
# lengthened edge loses _SHRINK of its excess per step. No node moves more
# than _MAX_STEP (below half a clearance) per step.
_ETA = 0.3
_MAX_TURN = 0.12
_MAX_MOVE = 0.6
_T_END = 0.02
_SHRINK = 0.06
_MAX_STEP = 0.3
_SAMPLE = 1.4  # edge sampling of the per step check
_GUARD_ROUNDS = 8
_HOLD, _HOLD_MAX = 2, 24  # a held back edge stays frozen _HOLD * 2^(fails - 1) steps
# Shortening passes: a stretched edge's subtree turned rigidly about the
# parent (in _SHORTEN_DEG steps, the smallest turn first) at its length to
# scale, or failing that at _SHORTEN_LEVELS - 1 lengths between.
_SHORTEN_DEG = 10.0
_SHORTEN_LEVELS = 3
_SHORTEN_MAX_SUBTREE = 60
_SHORTEN_FULL = 8  # a parent whose children miss this many times in a row is full for the pass
_SHORTEN_BATCH = 400_000  # candidates x nodes x edges tested at once


def _spring_settings(n: int) -> Tuple[int, int]:
    """(steps, final shortening passes) for a tree of n nodes: the full
    schedule up to 400 nodes; above, half the steps and one pass give the
    same share of edges to scale (within 0.5 %) in about 60 % of the time."""
    if n <= 400:
        return 120, 2
    return 60, 1


def _far_repulsion(P):
    """_K_FAR^2 * sum_j (p_i - p_j) / (|p_i - p_j|^2 + soft^2): all pairs up
    to _EXACT_MAX nodes, else against grid cell centroids (masses)."""
    n = len(P)
    z = P[:, 0] + 1j * P[:, 1]
    if n <= _EXACT_MAX:
        src, mass, s2 = z, None, _SOFT * _SOFT
    else:
        lo = P.min(0)
        size = max(float(np.ptp(P, axis=0).max()), 1e-6) / 32 * (1 + 1e-9)
        cell = np.floor((P - lo) / size).astype(np.int64)
        _uniq, inv = np.unique(cell[:, 0] * 33 + cell[:, 1], return_inverse=True)
        mass = np.bincount(inv).astype(float)
        src = (np.bincount(inv, weights=P[:, 0]) + 1j * np.bincount(inv, weights=P[:, 1])) / mass
        s2 = _SOFT * _SOFT + size * size
    out = np.empty(n, dtype=complex)
    for s0 in range(0, n, 256):
        diff = z[s0 : s0 + 256, None] - src[None, :]
        w = 1.0 / (diff.real * diff.real + diff.imag * diff.imag + s2)
        if mass is not None:
            w *= mass
        out[s0 : s0 + 256] = (diff * w).sum(1)
    return _K_FAR * _K_FAR * np.column_stack([out.real, out.imag])


def _forces(P, geom: _Geometry, center, far):
    F = far.copy()
    kd = cKDTree(P)
    pairs = kd.query_pairs(_CUTOFF, output_type="ndarray")
    if pairs.size:
        a, b = pairs[:, 0], pairs[:, 1]
        diff = P[a] - P[b]
        d = np.maximum(np.linalg.norm(diff, axis=1), 1e-6)
        mag = _K * _K / d * (1.0 - d / _CUTOFF) ** 2 + np.where(d < _BARRIER, 8.0 * (_BARRIER - d), 0.0)
        f = diff / d[:, None] * mag[:, None]
        np.add.at(F, a, f)
        np.add.at(F, b, -f)
    el, er = geom.el, geom.er
    if len(el):
        pts, eid = _segment_samples(P, el, er, 1.0)
        nodes, eids = _node_edge_pairs(kd, cKDTree(pts), eid, len(el), _EDGE_REACH + 0.5, el, er)
        if nodes.size:
            d, diff, t = _point_segment(P, el, er, nodes, eids)
            near = d < _EDGE_REACH
            nodes, eids, d, diff, t = nodes[near], eids[near], np.maximum(d[near], 1e-6), diff[near], t[near]
            f = diff / d[:, None] * (6.0 * (_EDGE_REACH - d))[:, None]
            np.add.at(F, nodes, f)
            np.add.at(F, el[eids], -f * (1.0 - t)[:, None])
            np.add.at(F, er[eids], -f * t[:, None])
    F -= _GRAVITY * (P - center)
    return F


def _shorten_pass(P, tree: _Tree, geom: _Geometry, target, node_min, clear):
    """Bring stretched edges back to their length to scale where there is
    room: the child's subtree turned rigidly about the parent, each
    candidate checked against the rest of the drawing (nodes, nodes on
    edges both ways, crossings). Small subtrees first. Returns the new
    positions and how many edges it shortened."""
    P = P.copy()
    el, er = geom.el, geom.er
    n = len(P)
    cur = np.linalg.norm(P[tree.nonroot] - P[tree.par[tree.nonroot]], axis=1)
    stretched = tree.nonroot[cur > 1.02 * target[tree.nonroot] + 1e-9]
    stretched = stretched[tree.size[stretched] <= _SHORTEN_MAX_SUBTREE]
    stretched = stretched[np.lexsort((stretched, tree.size[stretched]))]
    steps = np.radians(np.arange(0.0, 180.01, _SHORTEN_DEG))
    turns = np.concatenate([[0.0], np.ravel(np.column_stack([steps[1:], -steps[1:]]))])
    c, s = np.cos(turns), np.sin(turns)
    rot = np.stack([np.column_stack([c, -s]), np.column_stack([s, c])], axis=1)
    done = 0
    in_sub = np.zeros(n, bool)
    pos_in_s = np.full(n, -1, np.int64)
    misses = np.zeros(n, np.int64)  # per parent, children in a row that found no room
    for v in stretched.tolist():
        p = int(tree.par[v])
        S = tree.order[tree.tin[v] : tree.tout[v]]
        L_now = float(np.hypot(*(P[v] - P[p])))
        tgt = float(target[v])
        if L_now <= 1.02 * tgt or misses[p] >= _SHORTEN_FULL:
            continue
        in_sub[:] = False
        in_sub[S] = True
        heading = math.atan2(*(P[v] - P[p])[::-1])
        rel = P[S] - P[v]
        reach = float(np.max(np.linalg.norm(rel, axis=1))) + L_now + 2.0
        # The rest of the drawing near the parent: its nodes (R) and edges (fe).
        R = np.flatnonzero(~in_sub & (np.linalg.norm(P - P[p], axis=1) <= reach + 1.5))
        nodes_near = cKDTree(P[R]) if R.size else None
        moving_e = in_sub[er]  # (p, v) and the subtree's own edges
        fe = np.flatnonzero(~moving_e)
        fe = fe[_seg_dist(P[p], P[el[fe]], P[er[fe]]) <= reach + 1.0]
        fa, fb = P[el[fe]], P[er[fe]]
        f_lo, f_hi = np.minimum(fa, fb), np.maximum(fa, fb)
        fe_at_p = (el[fe] == p) | (er[fe] == p)
        me = np.flatnonzero(moving_e)
        pos_in_s[:] = -1
        pos_in_s[S] = np.arange(len(S))
        me_a = pos_in_s[el[me]]  # -1: the parent
        me_b = pos_in_s[er[me]]
        share = (me_a == -1)[:, None] & fe_at_p[None, :]  # edges meeting at the parent
        Rq = P[R[R != p]]  # p is on (p, v) only; the rest are foreign to every moving edge

        def ends(cands):
            """The moving edges' ends for each candidate placement."""
            return np.where(me_a[None, :, None] >= 0, cands[:, np.maximum(me_a, 0)], P[p]), cands[:, me_b]

        def first_fit(cands: np.ndarray) -> int:
            """The first of the candidate placements (C, len(S), 2) whose
            subtree keeps clear of the rest (its nodes already keep clear of
            the rest's nodes), or -1. All C are tested at once: edges outside
            the box around all of them (plus 1 > clear) can neither come
            within a clearance nor cross."""
            lo = np.minimum(cands.min(axis=(0, 1)), P[p]) - 1.0
            hi = np.maximum(cands.max(axis=(0, 1)), P[p]) + 1.0
            alive = np.arange(len(cands))  # the candidates still fitting; cheapest test first
            near = np.flatnonzero(np.all((f_lo <= hi) & (f_hi >= lo), axis=1))
            if near.size:
                na, nb = fa[near][None, None], fb[near][None, None]
                alive = alive[_seg_dist(cands[:, :, None, :], na, nb).min(axis=(1, 2)) >= clear]
            close = Rq[np.all((Rq >= lo) & (Rq <= hi), axis=1)]
            if alive.size and close.size:
                A, B = ends(cands[alive])
                alive = alive[_seg_dist(close[None, None], A[:, :, None, :], B[:, :, None, :]).min(axis=(1, 2)) >= clear]
            if alive.size and near.size:
                A, B = ends(cands[alive])
                X = _cross(A[:, :, None, :], B[:, :, None, :], na, nb)
                alive = alive[~(X & ~share[None, :, near]).any(axis=(1, 2))]
            return int(alive[0]) if alive.size else -1

        levels = [tgt] + [tgt * (L_now / tgt) ** (k / _SHORTEN_LEVELS) for k in range(1, _SHORTEN_LEVELS)]
        # Every placement in the order they are tried: each length, each turn.
        phi = heading + turns
        u = np.column_stack([np.cos(phi), np.sin(phi)])
        Snew = (
            P[p] + np.asarray(levels)[:, None, None, None] * u[None, :, None, :] + np.einsum("tij,kj->tki", rot, rel)[None]
        ).reshape(-1, len(S), 2)
        ok = np.ones(len(Snew), bool)
        if nodes_near is not None:
            d, _idx = nodes_near.query(Snew.reshape(-1, 2), distance_upper_bound=node_min)
            ok = (d.reshape(len(Snew), -1) >= node_min).all(axis=1)
        # The first candidate alone (it often fits), then the rest in
        # batches of bounded size.
        cand = np.flatnonzero(ok)
        batch = max(1, _SHORTEN_BATCH // (len(S) * max(1, len(fe), len(Rq))))
        start = 0
        while start < cand.size:
            group = cand[start : start + (1 if start == 0 else batch)]
            start += len(group)
            k = first_fit(Snew[group])
            if k >= 0:
                P[S] = Snew[group[k]]
                done += 1
                misses[p] = 0
                break
        else:
            misses[p] += 1
    return P, done


def _spring_layout(seed: np.ndarray, parents: np.ndarray, edges: Edges, scale: MstEdgeScale) -> np.ndarray:
    """Constrained springs in angle space, from the balloon ``seed``.

    The state is, per non-root node, the heading and length of the edge to
    its parent. The forces on a subtree turn its edge by their torque about
    the parent over the subtree's moment of inertia (prefix sums over the
    preorder, so a step is O(n)); lengths only shrink toward their length to
    scale. After each step the rules are checked (no crossing test needed:
    no node moves more than _MAX_STEP, so an edge can only start crossing
    another by a node coming within a clearance of it); each violation
    takes back the step on the two tree paths from the pair's lowest common
    ancestor and freezes them a while; if a few rounds do not clear the
    step, it is undone. Shortening passes before and after."""
    n = len(seed)
    steps_total, shorten_passes = _spring_settings(n)
    tree = _Tree(n, parents)
    geom = _Geometry(tree)
    target, theta, length = _edge_state(seed, parents, edges, scale, tree)
    root_pos = np.array(seed, dtype=float)
    center = seed[tree.roots].mean(axis=0)
    node_min = _NODE_MIN - 1e-6
    clear = _NODE_CLEARANCE - 1e-6
    nr = tree.nonroot

    def from_positions(P):
        vec = P[nr] - P[tree.par[nr]]
        theta[nr] = np.arctan2(vec[:, 1], vec[:, 0])
        length[nr] = np.hypot(vec[:, 0], vec[:, 1])
        return tree.positions(theta, length, root_pos)

    P, _done = _shorten_pass(tree.positions(theta, length, root_pos), tree, geom, target, node_min, clear)
    P = from_positions(P)
    decay = _T_END ** (1.0 / max(1, steps_total - 1))
    temp = 1.0
    hold = np.zeros(n, np.int64)
    fails = np.zeros(n, np.int64)
    far = np.zeros((n, 2))
    for it in range(steps_total):
        if it % _FAR_EVERY == 0:
            far = _far_repulsion(P)
        F = _forces(P, geom, center, far)
        sumF = tree.subtree_sum(F)
        sumC = tree.subtree_sum(P[:, 0] * F[:, 1] - P[:, 1] * F[:, 0])
        sumP = tree.subtree_sum(P)
        sumP2 = tree.subtree_sum((P * P).sum(1))
        pp = P[np.maximum(tree.par, 0)]
        torque = sumC - (pp[:, 0] * sumF[:, 1] - pp[:, 1] * sumF[:, 0])
        inertia = sumP2 - 2.0 * (pp * sumP).sum(1) + tree.size * (pp * pp).sum(1)
        cap = np.minimum(_MAX_TURN, _MAX_MOVE / np.maximum(tree.subtree_reach(P), 1.0)) * temp
        new_theta = theta.copy()
        new_theta[nr] += np.clip(_ETA * torque[nr] / np.maximum(inertia[nr], 1e-6), -cap[nr], cap[nr])
        new_length = length.copy()
        new_length[nr] = target[nr] + (length[nr] - target[nr]) * (1.0 - _SHRINK)
        held = hold > 0
        new_length[held] = length[held]
        new_theta[held] = theta[held]
        hold[held] -= 1
        snap = new_length - target < 1e-4
        new_length[snap] = target[snap]
        Pn = tree.positions(new_theta, new_length, root_pos)
        for _shrink in range(3):
            disp = float(np.max(np.linalg.norm(Pn - P, axis=1)))
            if disp <= _MAX_STEP:
                break
            f = 0.95 * _MAX_STEP / disp
            new_theta = theta + f * (new_theta - theta)
            new_length = length + f * (new_length - length)
            Pn = tree.positions(new_theta, new_length, root_pos)
        else:
            temp *= decay
            continue
        for rnd in range(_GUARD_ROUNDS):
            bad = geom.violations(Pn, node_min, clear, _SAMPLE, crossings=False)
            if not bad:
                theta, length, P = new_theta, new_length, Pn
                break
            revert = set()
            for a, b in bad:
                path = tree.path_nodes(int(a), int(b))
                revert.update(path)
                if rnd >= 2 and path:
                    # The lowest common ancestor's whole subtree goes back.
                    lca = int(tree.par[path[-1]])
                    revert.update(tree.order[tree.tin[lca] + 1 : tree.tout[lca]].tolist())
            idx = np.fromiter(sorted(revert), dtype=np.int64)
            fails[idx] += 1
            hold[idx] = np.minimum(_HOLD_MAX, _HOLD * 2 ** np.minimum(fails[idx] - 1, 6))
            new_theta[idx] = theta[idx]
            new_length[idx] = length[idx]
            Pn = tree.positions(new_theta, new_length, root_pos)
        temp *= decay
    for _pass in range(shorten_passes):
        P, done = _shorten_pass(P, tree, geom, target, node_min, clear)
        if not done:
            break
    return P


# ---------------------------------------------------------------- stress

_STRESS_STEP_MAX = 0.34  # a node moving less cannot carry an edge across another (< clearance / 2)
_SOFT_NN = 2.4  # soft spacing targets, above the hard ones
_SOFT_NE = 1.15
_SOFT_K = 0.6
_LEN_K = 0.5  # pull of a lengthened edge back to its length to scale
_LOCAL_PAIRS = 40  # stress partners per node (nearest in the tree), plus _PIVOTS far ones
_PIVOTS = 24
_MIN_ADJ_DEG = 10.0  # smallest angle between two edges at a node (less at a busy node)


def _tree_distances(tree: _Tree, target: np.ndarray) -> np.ndarray:
    """Path lengths between all nodes at the edges' lengths to scale (the
    stress targets), from one root."""
    n = tree.n
    order = tree.order
    D = np.zeros((n, n))
    root = int(order[0])
    for v in order[1:].tolist():
        D[root, v] = D[root, tree.par[v]] + target[v]
    for v in order[1:].tolist():
        D[v] = D[tree.par[v]] + target[v]
        D[v, order[tree.tin[v] : tree.tout[v]]] -= 2.0 * target[v]
    np.fill_diagonal(D, 0.0)
    return D


def _smacof(n, pi, pj, Dt, W, X, root, iters, tol=3e-4):
    """SMACOF (Guttman transform) over the given pairs, ``root`` held in place."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import splu

    diag = np.bincount(pi, weights=W, minlength=n) + np.bincount(pj, weights=W, minlength=n)
    rows = np.concatenate([pi, pj, np.arange(n)])
    cols = np.concatenate([pj, pi, np.arange(n)])
    Lw = coo_matrix((np.concatenate([-W, -W, diag]), (rows, cols)), shape=(n, n)).tocsc()
    free = np.flatnonzero(np.arange(n) != root)
    lu = splu(Lw[free][:, free].tocsc())
    col = Lw[free][:, [root]].toarray()
    X = X.copy()
    prev = None
    for _it in range(iters):
        diff = X[pi] - X[pj]
        d = np.maximum(np.hypot(diff[:, 0], diff[:, 1]), 1e-9)
        stress = float(np.sum(W * (d - Dt) ** 2))
        if prev is not None and prev - stress < tol * prev:
            break
        prev = stress
        c = (W * Dt / d)[:, None] * diff
        Y = np.zeros((n, 2))
        for k in range(2):
            Y[:, k] = np.bincount(pi, weights=c[:, k], minlength=n) - np.bincount(pj, weights=c[:, k], minlength=n)
        X[free] = lu.solve(np.ascontiguousarray(Y[free] - col @ X[[root]]))
    return X


def _angle_pairs(tree: _Tree):
    """Pairs of edges at a node: the node, the two other ends and the
    smallest angle allowed between them."""
    adj: List[List[int]] = [[] for _ in range(tree.n)]
    for c in tree.nonroot.tolist():
        adj[int(tree.par[c])].append(c)
        adj[c].append(int(tree.par[c]))
    V, O1, O2, M = [], [], [], []
    for v, nb in enumerate(adj):
        if len(nb) < 2:
            continue
        lim = min(math.radians(_MIN_ADJ_DEG), math.pi / len(nb))
        for x in range(len(nb)):
            for y in range(x + 1, len(nb)):
                V.append(v)
                O1.append(nb[x])
                O2.append(nb[y])
                M.append(lim)
    return np.asarray(V, np.int64), np.asarray(O1, np.int64), np.asarray(O2, np.int64), np.asarray(M)


def _pair_angles(P, V, O1, O2):
    u = P[O1] - P[V]
    w = P[O2] - P[V]
    return np.abs(np.arctan2(u[:, 0] * w[:, 1] - u[:, 1] * w[:, 0], np.einsum("ij,ij->i", u, w)))


def _wrap(a):
    return (a + np.pi) % (2.0 * np.pi) - np.pi


def _stress_layout(seed: np.ndarray, parents: np.ndarray, edges: Edges, scale: MstEdgeScale) -> np.ndarray:
    """Stress majorization of the tree, kept clear at every step.

    1. Unconstrained stress (SMACOF; targets: tree path lengths at the
       lengths to scale, weights 1 / D^2, over each node's nearest in the
       tree and a few pivots) from the balloon.
    2. Move each edge's heading toward it and its length to scale, as far
       as stays clear (a held back edge halves its share, a free one
       doubles it).
    3. Polish: stress descent plus a soft spacing term, in heading / length
       space, preconditioned by each subtree's moment.
    4. Snap lengthened edges to their length to scale where it stays clear.
    Each step keeps the rules: the edges below the lowest common ancestor
    of the nodes that break one are held back."""
    n = len(seed)
    tree = _Tree(n, parents)
    if len(tree.roots) != 1:
        raise ValueError("the stress layout needs a connected tree")
    geom = _Geometry(tree)
    target, theta, length = _edge_state(seed, parents, edges, scale, tree)
    length = np.maximum(length, target)
    root_pos = np.array(seed, dtype=float)
    root = int(tree.roots[0])
    nr = tree.nonroot
    P = tree.positions(theta, length, root_pos)
    node_min = _NODE_MIN - 1e-9
    clear = _NODE_CLEARANCE - 1e-9
    V, O1, O2, min_angle = _angle_pairs(tree)

    def step(theta, length, P, d_theta, d_length, rounds=8, priority=False):
        """The step with the edges that break a rule held back. With
        ``priority``, in the first round the first of those edges (in
        preorder) still moves. Returns the new state and the held edges."""
        blocked = (d_theta == 0.0) & (d_length == 0.0)
        for rnd in range(rounds):
            th = theta + np.where(blocked, 0.0, d_theta)
            ln = length + np.where(blocked, 0.0, d_length)
            Q = tree.positions(th, ln, root_pos)
            moved = float(np.max(np.abs(Q - P)))
            bad = geom.violations(Q, node_min, clear, 2.0, crossings=moved * 1.4143 > _STRESS_STEP_MAX)
            if V.size:
                before = _pair_angles(P, V, O1, O2)
                after = _pair_angles(Q, V, O1, O2)
                for k in np.flatnonzero((after < min_angle) & (after < before - 1e-9)).tolist():
                    bad.append((int(O1[k]), int(O2[k])))
            if not bad:
                return th, ln, Q, blocked
            count = int(blocked.sum())
            for a, b in bad:
                path = sorted((c for c in tree.path_nodes(int(a), int(b)) if not blocked[c]), key=lambda c: tree.tin[c])
                if priority and rnd == 0 and len(path) > 1:
                    path = path[1:]
                blocked[path] = True
            if int(blocked.sum()) == count and not (priority and rnd == 0):
                break
        return theta, length, P, np.ones(n, dtype=bool)

    # Stress pairs: each node's nearest in the tree, plus pivots spread over it.
    D = _tree_distances(tree, target)
    k = min(n - 1, _LOCAL_PAIRS)
    near = np.concatenate([np.argpartition(D[r0 : r0 + 256], k, axis=1)[:, : k + 1] for r0 in range(0, n, 256)])
    rows = np.repeat(np.arange(n), near.shape[1])
    cols = near.ravel()
    pivots = [root]
    far = D[root].copy()
    for _ in range(_PIVOTS - 1):
        pivots.append(int(np.argmax(far)))
        far = np.minimum(far, D[pivots[-1]])
    rows = np.concatenate([rows, np.repeat(np.asarray(pivots), n)])
    cols = np.concatenate([cols, np.tile(np.arange(n), len(pivots))])
    keep = rows != cols
    key = np.unique(np.minimum(rows[keep], cols[keep]) * n + np.maximum(rows[keep], cols[keep]))
    pi, pj = key // n, key % n
    Dt = D[pi, pj]
    W = 1.0 / (Dt * Dt)
    del D

    # 1. and 2.
    X = _smacof(n, pi, pj, Dt, W, P, root, 60)
    vec = X[nr] - X[tree.par[nr]]
    theta_t = theta.copy()
    theta_t[nr] = np.arctan2(vec[:, 1], vec[:, 0])
    share = np.ones(n)

    def advance(theta, length, P, d_theta, d_length, todo, **kw):
        theta2, length2, P2, blocked = step(theta, length, P, d_theta, d_length, **kw)
        went = todo & ~blocked
        share[todo & blocked] *= 0.5
        share[went] = np.minimum(1.0, share[went] * 2.0)
        return (theta2, length2, P2) if went.any() else (theta, length, P)

    for _it in range(200):
        gap_theta = _wrap(theta_t - theta)
        gap_length = target - length
        gap_theta[root] = gap_length[root] = 0.0
        reach = np.maximum(tree.subtree_reach(P), 1e-6)
        todo = ((np.abs(gap_theta) * reach > 1e-3) | (np.abs(gap_length) > 1e-6)) & (share >= 1.0 / 64.0)
        if not todo.any():
            break
        theta, length, P = advance(theta, length, P, np.where(todo, gap_theta * share, 0.0), np.where(todo, gap_length * share, 0.0), todo)

    # 3. Polish.
    weight = 2.0 * (np.bincount(pi, weights=W, minlength=n) + np.bincount(pj, weights=W, minlength=n)) + 2.0 * _SOFT_K
    edge_keys = np.unique(np.minimum(geom.el, geom.er) * n + np.maximum(geom.el, geom.er))

    def gradient(P):
        diff = P[pi] - P[pj]
        d = np.maximum(np.hypot(diff[:, 0], diff[:, 1]), 1e-9)
        f = (2.0 * W * (d - Dt) / d)[:, None] * diff
        g = np.zeros((n, 2))
        for k in range(2):
            g[:, k] = np.bincount(pi, weights=f[:, k], minlength=n) - np.bincount(pj, weights=f[:, k], minlength=n)
        close = cKDTree(P).query_pairs(r=_SOFT_NN, output_type="ndarray")
        if len(close):
            close = close[~np.isin(np.minimum(close[:, 0], close[:, 1]) * n + np.maximum(close[:, 0], close[:, 1]), edge_keys)]
        if len(close):
            dv = P[close[:, 0]] - P[close[:, 1]]
            dd = np.maximum(np.hypot(dv[:, 0], dv[:, 1]), 1e-9)
            ff = (-2.0 * _SOFT_K * (_SOFT_NN - dd) / dd)[:, None] * dv
            np.add.at(g, close[:, 0], ff)
            np.add.at(g, close[:, 1], -ff)
        pts, eid = _segment_samples(P, geom.el, geom.er, 2.0)
        nodes, eids = _node_edge_pairs(cKDTree(P), cKDTree(pts), eid, len(geom.el), _SOFT_NE + 1.0, geom.el, geom.er)
        if nodes.size:
            dd, diff, tt = _point_segment(P, geom.el, geom.er, nodes, eids)
            m = dd < _SOFT_NE
            nodes, eids, dd, diff, tt = nodes[m], eids[m], dd[m], diff[m], tt[m]
            ff = (2.0 * _SOFT_K * (_SOFT_NE - dd) / np.maximum(dd, 1e-9))[:, None] * diff
            np.add.at(g, nodes, -ff)
            np.add.at(g, geom.el[eids], ff * (1 - tt)[:, None])
            np.add.at(g, geom.er[eids], ff * tt[:, None])
        return g

    share = np.ones(n)
    for _it in range(15):
        g = gradient(P)
        S = tree.subtree_sum(g)
        Xs = tree.subtree_sum(P[:, 0] * g[:, 1] - P[:, 1] * g[:, 0])
        pp = P[np.maximum(tree.par, 0)]
        tau = Xs - (pp[:, 0] * S[:, 1] - pp[:, 1] * S[:, 0])
        g_length = np.einsum("ij,ij->i", S, np.column_stack([np.cos(theta), np.sin(theta)])) + 2.0 * _LEN_K * (length - target)
        A = tree.subtree_sum(weight)
        Ar2 = tree.subtree_sum(weight[:, None] * np.column_stack([(P * P).sum(1), P[:, 0], P[:, 1]]))
        H = Ar2[:, 0] - 2.0 * (pp[:, 0] * Ar2[:, 1] + pp[:, 1] * Ar2[:, 2]) + (pp * pp).sum(1) * A
        reach = np.maximum(tree.subtree_reach(P), 1e-6)
        d_theta = np.clip(-tau / np.maximum(H, 1e-6), -1.0 / reach, 1.0 / reach) * share
        d_length = np.minimum(0.0, np.maximum(target, length - g_length / (A + 2.0 * _LEN_K)) - length) * share
        d_theta[root] = d_length[root] = 0.0
        todo = (np.abs(d_theta) * reach > 2e-3) | (np.abs(d_length) > 2e-3)
        if not todo.any():
            break
        theta, length, P = advance(theta, length, P, np.where(todo, d_theta, 0.0), np.where(todo, d_length, 0.0), todo)

    # 4. Snap.
    for _round in range(8):
        d_length = np.maximum(target - length, -1.0)
        d_length[root] = 0.0
        if float(np.max(np.abs(d_length))) < 1e-9:
            break
        theta2, length2, P2, blocked = step(theta, length, P, np.zeros(n), d_length, priority=True)
        if blocked.all() or np.array_equal(length2, length):
            break
        theta, length, P = theta2, length2, P2
    return P
