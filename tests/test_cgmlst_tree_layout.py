"""The cgMLST MST's layout (cgmlst_mst_layout.mst_drawing: the balloon of
cgmlst_graph_data._mst_tree_drawing relaxed by constrained springs, or by
stress on large trees): whatever the tree, no two edges cross, nodes stay
apart and off the edges they are not on, and the same input gives the same
picture. Identical profiles form a chain, and long unbranched chains stay a
compact block. Edge lengths are linear in the allele distance, outliers drawn
shortened (broken); an edge lengthened to make room is to scale in its solid
part.
"""
from __future__ import annotations

import unittest

import numpy as np

from annotation_explorer_polars.cgmlst_graph_data import (
    _NODE_CLEARANCE,
    _build_mst_edges,
    _edge_scale_parts,
    _mst_edge_scale,
    _mst_tree_drawing,
)
from annotation_explorer_polars.cgmlst_mst_layout import (
    _clarity_problems,
    _spring_layout,
    _stress_layout,
    mst_drawing,
)


def _layout(n, edges, labels, scale=None):
    """The app's layout; it must not have fallen back to the balloon."""
    pos, parents, method = mst_drawing(n, edges, labels, scale)
    if n >= 3 and method not in ("springs", "stress"):
        raise AssertionError(method)
    return pos, parents


def _random_trees():
    rng = np.random.default_rng(492)
    for shape in ("random", "star", "chain"):
        for n in (2, 3, 40, 160):
            edges = []
            for child in range(1, n):
                if shape == "star":
                    parent = 0 if rng.random() < 0.7 else int(rng.integers(0, child))
                elif shape == "chain":
                    parent = child - 1 if rng.random() < 0.8 else int(rng.integers(0, child))
                else:
                    parent = int(rng.integers(0, child))
                # Identical profiles, small distances and a few very long edges.
                roll = rng.random()
                distance = 0 if roll < 0.3 else int(rng.integers(1, 40)) if roll < 0.95 else int(rng.integers(500, 2400))
                edges.append((parent, child, distance))
            yield f"{shape}{n}", n, edges


def _crossing(p1, p2, q1, q2) -> bool:
    def turn(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    return turn(p1, p2, q1) * turn(p1, p2, q2) < 0 and turn(q1, q2, p1) * turn(q1, q2, p2) < 0


def _segment_distance(point, start, end) -> float:
    span = end - start
    t = float(np.clip(np.dot(point - start, span) / max(float(np.dot(span, span)), 1e-12), 0.0, 1.0))
    return float(np.linalg.norm(point - (start + t * span)))


def _clear_layout_checks(test: unittest.TestCase, n, edges, pos) -> None:
    tolerance = 1e-6
    gaps = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=2)
    np.fill_diagonal(gaps, np.inf)
    test.assertGreaterEqual(float(gaps.min()), 2.0 * _NODE_CLEARANCE - tolerance)
    for idx, (left, right, _distance) in enumerate(edges):
        for node in range(n):
            if node not in (left, right):
                test.assertGreaterEqual(_segment_distance(pos[node], pos[left], pos[right]), _NODE_CLEARANCE - tolerance)
        for left2, right2, _distance2 in edges[idx + 1:]:
            if {left, right} & {left2, right2}:
                continue
            test.assertFalse(_crossing(pos[left], pos[right], pos[left2], pos[right2]))


class MstTreeLayoutTest(unittest.TestCase):
    def test_identical_profiles_form_a_block_not_a_ring(self) -> None:
        # 40 identical genomes and one 25 alleles away: Prim makes a star of 39
        # leaves around one of them; they become a chain, drawn as a block.
        n = 41
        score = np.zeros((n, n), dtype=np.int64)
        score[40, :40] = score[:40, 40] = 25 * 2391
        dist = (score // 2391).astype(np.int32)
        labels = [f"S{idx:03d}" for idx in range(n)]
        edges = _build_mst_edges(score, dist, labels)
        self.assertEqual(len(edges), n - 1)
        self.assertEqual(sum(distance for _l, _r, distance in edges), 25)
        degree = np.bincount([node for left, right, _d in edges for node in (left, right)], minlength=n)
        self.assertLessEqual(int(degree.max()), 3)
        for left, right, distance in edges:
            self.assertEqual(int(dist[left, right]), distance)
        pos, _parents = _layout(n, edges, labels)
        _clear_layout_checks(self, n, edges, pos)
        extent = np.ptp(pos[:40], axis=0)
        # The star drawn as is: a ring about 19 units across. A line: 55 long.
        # (The balloon folds the chain into a block; the springs loosen it
        # into a curled strand, still well short of a line.)
        self.assertLess(float(extent[0] * extent[1]), 19.0 * 19.0)
        self.assertLess(float(extent.max()), 0.6 * 55.0)

    def test_a_long_chain_stays_compact(self) -> None:
        n = 300
        edges = [(idx - 1, idx, idx % 7) for idx in range(1, n)]
        labels = [f"S{idx:03d}" for idx in range(n)]
        pos, _parents = _layout(n, edges, labels)
        _clear_layout_checks(self, n, edges, pos)
        extent = np.ptp(pos, axis=0)
        self.assertLess(float(extent.max() / extent.min()), 4.0)

    def test_layout_has_no_crossings_and_keeps_nodes_apart(self) -> None:
        tolerance = 1e-6
        for name, n, edges in _random_trees():
            labels = [f"S{idx:03d}" for idx in range(n)]
            pos, _parents = _layout(n, edges, labels)
            with self.subTest(tree=name):
                self.assertEqual(pos.shape, (n, 2))
                self.assertTrue(np.array_equal(pos, _layout(n, edges, labels)[0]))
                gaps = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=2)
                np.fill_diagonal(gaps, np.inf)
                self.assertGreaterEqual(float(gaps.min()), 2.0 * _NODE_CLEARANCE - tolerance)
                for idx, (left, right, _distance) in enumerate(edges):
                    for node in range(n):
                        if node not in (left, right):
                            self.assertGreaterEqual(
                                _segment_distance(pos[node], pos[left], pos[right]),
                                _NODE_CLEARANCE - tolerance,
                            )
                    for left2, right2, _distance2 in edges[idx + 1:]:
                        if {left, right} & {left2, right2}:
                            continue
                        self.assertFalse(_crossing(pos[left], pos[right], pos[left2], pos[right2]))

    def test_edge_lengths_are_linear_and_outliers_broken(self) -> None:
        # A tree of mostly 0 to 60 allele edges with a busy hub (its short
        # edges must be lengthened to fit), and three far outliers.
        rng = np.random.default_rng(1454)
        n = 160
        edges = []
        for child in range(1, n):
            parent = 0 if child % 5 == 0 else int(rng.integers(max(0, child - 6), child))
            distance = int(rng.integers(0, 61)) if child not in (41, 91, 141) else int(rng.integers(1400, 2000))
            edges.append((parent, child, distance))
        labels = [f"S{idx:03d}" for idx in range(n)]
        distances = np.array([distance for _l, _r, distance in edges], dtype=float)
        scale = _mst_edge_scale(distances)
        broken = np.zeros(len(edges), dtype=bool)
        broken[list(scale.broken_edges(distances))] = True
        self.assertEqual(set(np.flatnonzero(distances >= 1400)), set(np.flatnonzero(broken)))
        unit = scale.units_per_allele
        self.assertAlmostEqual(scale.length(40), 2.0 * scale.length(20))
        self.assertEqual(scale.length(0), scale.min_length)
        pos, parents = _layout(n, edges, labels, scale)
        _clear_layout_checks(self, n, edges, pos)
        drawn = np.array([np.linalg.norm(pos[right] - pos[left]) for left, right, _d in edges])
        wanted = np.array([scale.length(d) for d in distances])
        self.assertTrue(np.all(drawn >= wanted - 1e-6))
        # A lengthened edge is solid for its length to scale, from its
        # parent's end: every solid part is to scale (within the 2 % the
        # view leaves undotted), on one linear scale above the minimum.
        start, share = _edge_scale_parts(pos, edges, parents, scale)
        self.assertTrue(np.any(share < 1.0))
        for (left, right, _d), node in zip(edges, start):
            other = right if node == left else left
            self.assertIn(node, (left, right))
            self.assertEqual(int(parents[other]), int(node))
        solid = drawn * share
        np.testing.assert_allclose(solid, wanted, rtol=0.021)
        linear = ~broken & (distances * unit >= scale.min_length)
        np.testing.assert_allclose(solid[linear] / distances[linear], unit, rtol=0.021)
        # A broken edge is drawn shorter than to scale, longer than every
        # other edge's length to scale, and longer the farther it is.
        self.assertTrue(np.all(drawn[broken] < unit * distances[broken]))
        self.assertGreater(float(solid[broken].min()), float(solid[~broken].max()))
        order = np.argsort(distances[broken])
        self.assertTrue(np.all(np.diff(wanted[broken][order]) > 0.0))

    def test_springs_and_stress_keep_the_balloon_clear_and_draw_more_to_scale(self) -> None:
        # Both relaxations start from the balloon; each result passes the
        # full clarity check, is the same on every run and draws at least
        # as many edges at their length to scale as the balloon did.
        for name, n, edges in _random_trees():
            if n < 40:
                continue
            labels = [f"S{idx:03d}" for idx in range(n)]
            scale = _mst_edge_scale([distance for _l, _r, distance in edges])
            seed, parents = _mst_tree_drawing(n, edges, labels, scale)
            wanted = np.array([scale.length(distance) for _l, _r, distance in edges])
            left = np.array([edge[0] for edge in edges])
            right = np.array([edge[1] for edge in edges])

            def to_scale(pos):
                drawn = np.linalg.norm(pos[right] - pos[left], axis=1)
                self.assertTrue(np.all(drawn >= wanted - 1e-6))
                return int(np.count_nonzero(drawn <= 1.01 * wanted))

            for layout in (_spring_layout, _stress_layout):
                with self.subTest(tree=name, layout=layout.__name__):
                    pos = layout(seed, parents, edges, scale)
                    self.assertEqual(_clarity_problems(pos, parents), "")
                    self.assertTrue(np.array_equal(pos, layout(seed, parents, edges, scale)))
                    self.assertGreaterEqual(to_scale(pos), to_scale(seed))


if __name__ == "__main__":
    unittest.main()
