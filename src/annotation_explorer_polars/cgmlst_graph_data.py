from __future__ import annotations

import math
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .panaroo_common import _file_signature, _natural_sort_key


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


def _circle_layout(n_samples: int) -> np.ndarray:
    if n_samples <= 0:
        return np.zeros((0, 2), dtype=float)
    if n_samples == 1:
        return np.zeros((1, 2), dtype=float)
    angles = np.linspace(0.0, 2.0 * np.pi, num=n_samples, endpoint=False, dtype=float)
    coords = np.column_stack([np.cos(angles), np.sin(angles)])
    return coords.astype(float, copy=False)


def _build_mst_edges(
    score: np.ndarray,
    display_dist: np.ndarray,
    sample_labels: Sequence[str],
    *,
    reference_sample: Optional[str] = None,
) -> List[Tuple[int, int, int]]:
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

    used = np.zeros(n, dtype=bool)
    best_score = np.full(n, np.iinfo(np.int64).max, dtype=np.int64)
    parent = np.full(n, start, dtype=np.int32)
    used[start] = True
    best_score[:] = score[start]
    best_score[start] = np.iinfo(np.int64).max

    edges: List[Tuple[int, int, int]] = []
    while len(edges) < n - 1:
        candidates = np.flatnonzero(~used)
        if candidates.size == 0:
            break
        next_idx = min(
            (int(idx) for idx in candidates),
            key=lambda idx: (
                int(best_score[idx]),
                labels[int(parent[idx])],
                labels[idx],
            ),
        )
        parent_idx = int(parent[next_idx])
        edges.append((parent_idx, next_idx, int(display_dist[parent_idx, next_idx])))
        used[next_idx] = True
        remaining = np.flatnonzero(~used)
        for idx in remaining:
            candidate_score = int(score[next_idx, idx])
            current_parent = int(parent[idx])
            current_score = int(best_score[idx])
            if candidate_score < current_score or (
                candidate_score == current_score
                and labels[next_idx] < labels[current_parent]
            ):
                best_score[idx] = candidate_score
                parent[idx] = next_idx
    return edges


def _build_mst_adjacency(
    n: int,
    edges: Sequence[Tuple[int, int, int]],
) -> List[List[Tuple[int, float]]]:
    adjacency: List[List[Tuple[int, float]]] = [[] for _ in range(max(0, int(n)))]
    for left, right, dist in edges:
        if left < 0 or right < 0 or left >= n or right >= n:
            continue
        weight = max(0.0, float(dist))
        adjacency[left].append((right, weight))
        adjacency[right].append((left, weight))
    return adjacency


def _rooted_radial_mst_layout(
    sample_labels: Sequence[str],
    edges: Sequence[Tuple[int, int, int]],
    *,
    reference_idx: Optional[int] = None,
) -> np.ndarray:
    n = len(sample_labels)
    if n <= 1:
        return np.zeros((n, 2), dtype=float)

    adjacency = _build_mst_adjacency(n, edges)
    indegree = np.zeros(n, dtype=np.int32)
    for _left, right, _dist in edges:
        if 0 <= int(right) < n:
            indegree[int(right)] += 1
    if reference_idx is not None and 0 <= int(reference_idx) < n:
        root = int(reference_idx)
    else:
        roots = np.flatnonzero(indegree == 0)
        root = int(roots[0]) if roots.size else 0

    parent = np.full(n, -1, dtype=np.int32)
    level = np.zeros(n, dtype=np.int32)
    depth = np.zeros(n, dtype=float)
    order: List[int] = [root]
    parent[root] = root
    for node in order:
        neighbors = sorted(
            adjacency[node],
            key=lambda item: _natural_sort_key(sample_labels[int(item[0])]),
        )
        for child, weight in neighbors:
            if child == int(parent[node]) or parent[child] != -1:
                continue
            parent[child] = node
            level[child] = level[node] + 1
            depth[child] = depth[node] + 0.9 + math.sqrt(max(0.0, float(weight)))
            order.append(child)
    if len(order) != n:
        return _circle_layout(n)

    leaf_weight = np.ones(n, dtype=float)
    for node in reversed(order):
        children = [child for child, _weight in adjacency[node] if parent[child] == node]
        if children:
            leaf_weight[node] = float(sum(leaf_weight[child] for child in children))

    angles = np.zeros(n, dtype=float)

    def assign_angles(node: int, start: float, end: float) -> None:
        children = [child for child, _weight in adjacency[node] if parent[child] == node]
        if not children:
            angles[node] = 0.5 * (start + end)
            return
        children.sort(key=lambda child: (-leaf_weight[child], _natural_sort_key(sample_labels[child])))
        span = max(1e-6, float(end - start))
        gap = min(0.08 * span, 0.18)
        usable_span = max(1e-6, span - (gap * max(0, len(children) - 1)))
        total_weight = max(1.0, float(sum(leaf_weight[child] for child in children)))
        cursor = float(start)
        first_start = cursor
        last_end = cursor
        for idx, child in enumerate(children):
            fraction = float(leaf_weight[child]) / total_weight
            child_span = usable_span * fraction
            child_start = cursor
            child_end = child_start + child_span
            if idx == len(children) - 1:
                child_end = end
            assign_angles(child, child_start, child_end)
            if idx == 0:
                first_start = child_start
            last_end = child_end
            cursor = child_end + gap
        angles[node] = 0.5 * (first_start + last_end)

    assign_angles(root, -math.pi, math.pi)
    angles -= angles[root]

    max_level = max(1, int(np.max(level)) if level.size else 0)
    max_depth = max(1.0, float(np.max(depth)) if depth.size else 0.0)
    coords = np.zeros((n, 2), dtype=float)
    for node in range(n):
        level_component = float(level[node]) / float(max_level)
        depth_component = float(depth[node]) / float(max_depth)
        radius = (0.32 * level_component) + (0.68 * depth_component)
        if node == root:
            radius = 0.0
        coords[node] = radius * np.array(
            [math.cos(float(angles[node])), math.sin(float(angles[node]))],
            dtype=float,
        )
    return coords


def _spread_duplicate_profile_nodes(
    coords: np.ndarray,
    sample_labels: Sequence[str],
    allele_values: np.ndarray,
    allele_valid: np.ndarray,
) -> np.ndarray:
    if int(coords.shape[0]) <= 1:
        return coords
    grouped: Dict[Tuple[Tuple[int, ...], Tuple[bool, ...]], List[int]] = {}
    for idx, _label in enumerate(sample_labels):
        key = (
            tuple(int(value) for value in allele_values[idx].tolist()),
            tuple(bool(value) for value in allele_valid[idx].tolist()),
        )
        grouped.setdefault(key, []).append(idx)

    updated = np.array(coords, dtype=float, copy=True)
    for members in grouped.values():
        if len(members) <= 1:
            continue
        members_sorted = sorted(members, key=lambda idx: _natural_sort_key(sample_labels[idx]))
        center = np.mean(updated[members_sorted], axis=0)
        radius = 0.14 + (0.03 * max(0, len(members_sorted) - 2))
        angles = np.linspace(0.0, 2.0 * np.pi, num=len(members_sorted), endpoint=False, dtype=float)
        for idx, angle in zip(members_sorted, angles):
            updated[idx] = center + (radius * np.array([math.cos(angle), math.sin(angle)], dtype=float))
    return updated
