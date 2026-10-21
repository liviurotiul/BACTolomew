"""The Pangenome clustering compute (model/pan_clustering.py) against plain numpy."""
from __future__ import annotations

import unittest

import numpy as np

from annotation_explorer_polars.model import pan_clustering as pc


class PanClusteringComputeTest(unittest.TestCase):
    def test_range_jaccard_and_mds_match_numpy(self) -> None:
        rng = np.random.default_rng(7)
        n, g = 40, 300
        # Genes x samples, as Panaroo writes it; gene prevalences spread over 0..1.
        presence = (rng.random((g, n)) < rng.random(g)[:, None]).astype(np.uint8)
        # Genes at exactly 30% and 50% stay out (strict bounds), and two
        # isolates carry none of the genes in range.
        presence[0] = 0
        presence[0, :12] = 1  # 12 / 40 = 30%
        presence[1] = 0
        presence[1, :20] = 1  # 20 / 40 = 50%
        names = [f"S{i:02d}" for i in range(n)]
        genes = [f"g{i:03d}" for i in range(g)]
        prevalence = presence.mean(axis=1)
        in_range = (prevalence > 0.3) & (prevalence < 0.5)
        presence[np.ix_(in_range, [38, 39])] = 0
        prevalence = presence.mean(axis=1)
        in_range = (prevalence > 0.3) & (prevalence < 0.5)

        slice_ = pc.presence_slice(presence, names, genes, [""] * g)
        mask = pc.genes_in_range(slice_.prevalence, 0.3, 0.5)
        self.assertEqual(sorted(np.asarray(slice_.genes)[mask]), sorted(np.asarray(genes)[in_range]))
        self.assertNotIn("g000", np.asarray(slice_.genes)[mask])
        self.assertNotIn("g001", np.asarray(slice_.genes)[mask])
        # A bound off by rounding (0.35000000000000003) still leaves out a gene at 35%.
        self.assertFalse(pc.genes_in_range(np.array([7 / 20]), 0.05, 0.35000000000000003)[0])
        self.assertFalse(pc.genes_in_range(np.array([7 / 20]), 0.34999999999999997, 0.5)[0])

        result = pc.ClusteringCache().compute(slice_, "all", 0.3, 0.5, k=4)
        self.assertEqual(result.dropped, ("S38", "S39"))
        self.assertEqual(result.n_genes, int(in_range.sum()))

        # Jaccard against a pairwise loop over the kept isolates.
        x = presence[in_range][:, :38].T.astype(bool)
        expected = np.zeros((38, 38))
        for i in range(38):
            for j in range(38):
                union = np.logical_or(x[i], x[j]).sum()
                expected[i, j] = 1.0 - np.logical_and(x[i], x[j]).sum() / union if union else 0.0
        got = pc.jaccard_distance(slice_.matrix[:38][:, mask])
        np.testing.assert_allclose(got, expected, atol=1e-12)

        # Classical MDS against a full eigendecomposition: same distances
        # between the points on the three axes, same explained shares.
        d2 = expected ** 2
        j = np.eye(38) - 1.0 / 38
        b = -0.5 * j @ d2 @ j
        values, vectors = np.linalg.eigh(b)
        top = np.argsort(values)[::-1][:3]
        ref = vectors[:, top] * np.sqrt(values[top])
        np.testing.assert_allclose(np.abs(result.coords), np.abs(ref), atol=1e-8)
        np.testing.assert_allclose(result.explained, values[top] / np.trace(b), atol=1e-10)

        # A Euclidean configuration is recovered exactly (up to rotation).
        points = rng.normal(size=(25, 3)) * [3.0, 2.0, 1.0]
        euclid = np.sqrt(((points[:, None] - points[None]) ** 2).sum(axis=-1))
        coords, _explained = pc.classical_mds(euclid, 3)
        recovered = np.sqrt(((coords[:, None] - coords[None]) ** 2).sum(axis=-1))
        np.testing.assert_allclose(recovered, euclid, atol=1e-9)

    def test_kmeans_is_deterministic_and_numbered_by_size(self) -> None:
        rng = np.random.default_rng(3)
        centers = np.array([[0, 0, 0], [10, 0, 0], [0, 10, 0], [0, 0, 10], [10, 10, 0], [10, 0, 10]], dtype=float)
        sizes = [30, 25, 20, 15, 10, 5]
        points = np.vstack([c + 0.3 * rng.normal(size=(s, 3)) for c, s in zip(centers, sizes)])
        truth = np.repeat(np.arange(6), sizes)
        order = rng.permutation(len(points))
        points, truth = points[order], truth[order]

        first = pc.kmeans(points, k=6, seed=1)
        np.testing.assert_array_equal(first, pc.kmeans(points, k=6, seed=1))
        self.assertEqual(sorted(set(first.tolist())), [1, 2, 3, 4, 5, 6])
        # The blobs come back whole, numbered from the largest.
        for blob, size in enumerate(sizes):
            labels = set(first[truth == blob].tolist())
            self.assertEqual(labels, {blob + 1})
            self.assertEqual(int((first == blob + 1).sum()), size)
        # k is capped at the number of distinct points.
        self.assertEqual(set(pc.kmeans(np.zeros((5, 3)), k=6).tolist()), {1})


if __name__ == "__main__":
    unittest.main()
