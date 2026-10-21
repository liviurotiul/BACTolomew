"""How a part of the pangenome makes the isolates cluster (no Qt).

The Pangenome clustering panel's compute, after the user's R script (vegan
and plotly): the Panaroo presence/absence matrix (an empty cell is absent),
the genes whose prevalence is strictly between ``low`` and ``high`` (fractions
of the analyzed isolates), the isolates with no gene in that range dropped
(and reported), Jaccard distances between the binary profiles, classical
(Torgerson) MDS on three axes and k-means on the MDS coordinates with a fixed
seed. Per cluster, the genes present in the cluster and absent from every
other analyzed isolate ("unique genes").

Everything here is pure numpy (scipy for the top eigenpairs) and thread
safe: the panel runs ``ClusteringCache.compute`` on a pool thread.
"""
from __future__ import annotations

import threading
import time
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Dict, Hashable, List, Optional, Sequence, Tuple

import numpy as np
import polars as pl
from scipy import linalg as sla

# The script's range and k.
DEFAULT_LOW = 0.30
DEFAULT_HIGH = 0.50
DEFAULT_K = 6
DEFAULT_SEED = 1
# k-means++ starts; 10 left the script's 30% to 50% clustering seed dependent
# (a local optimum about 20% worse for some seeds), 50 reaches the best one.
DEFAULT_N_INIT = 50
MDS_AXES = 3
MIN_ISOLATES = 3
# Prevalences are count / n; two of them differ by far more than this.
_PREVALENCE_EPS = 1e-9


@dataclass(frozen=True, eq=False)
class PresenceSlice:
    """The presence matrix of the analyzed isolates, genes in prevalence order
    (highest first, then by name). Genes absent from every analyzed isolate
    are left out. ``matrix`` is a read-only bool isolates x genes array."""

    samples: Tuple[str, ...]
    genes: Tuple[str, ...]
    annotations: Tuple[str, ...]
    matrix: np.ndarray
    prevalence: np.ndarray

    @property
    def n_isolates(self) -> int:
        return len(self.samples)


def presence_slice(
    presence: np.ndarray,
    sample_names: Sequence[str],
    genes: Sequence[str],
    annotations: Sequence[str],
    samples: Optional[Sequence[str]] = None,
) -> PresenceSlice:
    """The slice of ``presence`` (genes x samples, as Panaroo writes it) for
    ``samples`` (None: every sample of the matrix; names the matrix lacks are
    skipped)."""
    column = {str(name): idx for idx, name in enumerate(sample_names)}
    if samples is None:
        chosen = [str(name) for name in sample_names]
    else:
        seen = set()
        chosen = []
        for name in samples:
            name = str(name)
            if name in column and name not in seen:
                seen.add(name)
                chosen.append(name)
    cols = np.fromiter((column[name] for name in chosen), dtype=np.intp, count=len(chosen))
    full = np.asarray(presence)
    if full.dtype == np.uint8 and full.size and int(full.max()) <= 1:
        full = full.view(bool)  # Panaroo's 0/1 matrix, without a bool copy
    sub = full[:, cols]  # genes x isolates
    counts = np.count_nonzero(sub, axis=1)
    n = max(1, len(chosen))
    prevalence = counts / float(n)
    present = np.flatnonzero(counts > 0)
    names = np.asarray([str(genes[i]) for i in present], dtype=object)
    # Highest prevalence first, then by gene name.
    order = present[np.lexsort((names.astype(str), -prevalence[present]))] if present.size else present
    picked = sub[order]
    del sub
    matrix = np.ascontiguousarray(picked.T, dtype=bool)
    del picked
    matrix.setflags(write=False)
    prev = np.ascontiguousarray(prevalence[order])
    prev.setflags(write=False)
    return PresenceSlice(
        samples=tuple(chosen),
        genes=tuple(str(genes[i]) for i in order),
        annotations=tuple(str(annotations[i]) if i < len(annotations) else "" for i in order),
        matrix=matrix,
        prevalence=prev,
    )


def genes_in_range(prevalence: np.ndarray, low: float, high: float) -> np.ndarray:
    """Genes with ``low < prevalence < high`` (strict, as the script, with a
    tolerance so a bound such as 0.35000000000000003 still leaves out a gene
    at 35%); a ``high`` of 1 (100%) keeps the genes of every isolate too."""
    prevalence = np.asarray(prevalence, dtype=float)
    upper = (
        prevalence < float(high) - _PREVALENCE_EPS if float(high) < 1.0 else np.ones(prevalence.shape, dtype=bool)
    )
    return (prevalence > float(low) + _PREVALENCE_EPS) & upper


def jaccard_distance(x: np.ndarray) -> np.ndarray:
    """Jaccard distances between the rows of a binary matrix (0 between two
    empty rows), with one matrix product and two n x n buffers."""
    a = np.asarray(x, dtype=np.float32)
    counts = a.sum(axis=1, dtype=np.float64)
    dist = (a @ a.T).astype(np.float64)  # the intersections
    del a
    union = np.add.outer(counts, counts)
    union -= dist
    empty = union <= 0
    union[empty] = 1.0
    dist /= union
    del union
    np.subtract(1.0, dist, out=dist)
    dist[empty] = 0.0
    np.fill_diagonal(dist, 0.0)
    return dist


def classical_mds(dist: np.ndarray, k: int = MDS_AXES) -> Tuple[np.ndarray, np.ndarray]:
    """Classical (Torgerson) MDS of a distance matrix, as R's ``cmdscale``:
    the coordinates on the top ``k`` axes (an axis with a negative eigenvalue
    gets 0 coordinates) and each axis' share of the total inertia
    (eigenvalue over the trace of the double-centered matrix). Each axis is
    signed so its largest absolute coordinate is positive.

    The trace counts the negative eigenvalues a non Euclidean distance such
    as Jaccard gives, as R's ``eig / sum(eig)`` does, so the shares are those
    of the script."""
    b = np.square(np.asarray(dist, dtype=np.float64))
    n = b.shape[0]
    k = int(max(1, min(k, n)))
    rows = b.mean(axis=1)
    # Double centering in place: b = -0.5 * (d2 - row - col + mean).
    b -= rows[:, None]
    b -= rows[None, :]
    b += rows.mean()
    b *= -0.5
    values, vectors = sla.eigh(b, subset_by_index=[n - k, n - 1])
    values, vectors = values[::-1], vectors[:, ::-1]
    coords = vectors * np.sqrt(np.clip(values, 0.0, None))[None, :]
    for axis in range(k):
        col = coords[:, axis]
        if col.size and col[int(np.argmax(np.abs(col)))] < 0:
            coords[:, axis] = -col
    trace = float(np.trace(b))
    explained = np.clip(values, 0.0, None) / trace if trace > 0 else np.zeros(k)
    return coords, explained


def kmeans(
    coords: np.ndarray,
    k: int = DEFAULT_K,
    seed: int = DEFAULT_SEED,
    n_init: int = DEFAULT_N_INIT,
    max_iter: int = 100,
) -> np.ndarray:
    """Lloyd k-means with k-means++ starts from ``seed``; the best of
    ``n_init`` starts. Labels 1..k numbered by cluster size (largest first,
    ties by first member), so the same input always gives the same labels.
    ``k`` is capped at the number of distinct points."""
    x = np.asarray(coords, dtype=np.float64)
    n = x.shape[0]
    if n == 0:
        return np.zeros(0, dtype=np.int32)
    k = int(max(1, min(int(k), np.unique(x, axis=0).shape[0])))
    rng = np.random.default_rng(int(seed))
    best_labels: Optional[np.ndarray] = None
    best_inertia = np.inf
    sq = (x * x).sum(axis=1)
    for _start in range(max(1, int(n_init))):
        centers = np.empty((k, x.shape[1]))
        centers[0] = x[rng.integers(n)]
        closest = ((x - centers[0]) ** 2).sum(axis=1)
        for c in range(1, k):
            total = closest.sum()
            idx = int(rng.choice(n, p=closest / total)) if total > 0 else int(rng.integers(n))
            centers[c] = x[idx]
            closest = np.minimum(closest, ((x - centers[c]) ** 2).sum(axis=1))
        labels = np.full(n, -1)
        for _it in range(max_iter):
            d = sq[:, None] - 2.0 * x @ centers.T + (centers * centers).sum(axis=1)[None, :]
            new = np.argmin(d, axis=1)
            if np.array_equal(new, labels):
                break
            labels = new
            # Emptied clusters restart at the points farthest from their
            # centers, each at a different one.
            sizes = np.bincount(labels, minlength=k)
            sums = np.stack([np.bincount(labels, weights=x[:, j], minlength=k) for j in range(x.shape[1])], axis=1)
            filled = sizes > 0
            centers[filled] = sums[filled] / sizes[filled, None]
            if not filled.all():
                far = np.argsort(-d[np.arange(n), labels], kind="stable")
                for reseeded, c in enumerate(np.flatnonzero(~filled)):
                    centers[c] = x[int(far[min(reseeded, n - 1)])]
        inertia = float(((x - centers[labels]) ** 2).sum())
        if inertia < best_inertia - 1e-12:
            best_inertia, best_labels = inertia, labels.copy()
    assert best_labels is not None
    used = np.unique(best_labels)
    sizes = [(-(best_labels == c).sum(), int(np.flatnonzero(best_labels == c)[0]), int(c)) for c in used]
    renumber = {c: rank + 1 for rank, (_size, _first, c) in enumerate(sorted(sizes))}
    return np.asarray([renumber[int(c)] for c in best_labels], dtype=np.int32)


@dataclass(frozen=True)
class UniqueGene:
    cluster: int
    gene: str
    annotation: str
    n_present: int
    cluster_size: int


def unique_genes(
    matrix: np.ndarray, labels: np.ndarray, genes: Sequence[str], annotations: Sequence[str]
) -> List[UniqueGene]:
    """Per cluster, the genes (columns of ``matrix``, isolates x genes, rows
    in ``labels`` order) present in at least one member and in no isolate of
    another cluster."""
    x = np.asarray(matrix, dtype=bool)
    labels = np.asarray(labels)
    out: List[UniqueGene] = []
    if x.size == 0:
        return out
    per_cluster = {int(c): x[labels == c].sum(axis=0) for c in np.unique(labels)}
    for cluster in sorted(per_cluster):
        inside = per_cluster[cluster]
        outside = x[labels != cluster].any(axis=0)
        size = int((labels == cluster).sum())
        for idx in np.flatnonzero((inside > 0) & ~outside):
            out.append(UniqueGene(cluster, str(genes[idx]), str(annotations[idx]), int(inside[idx]), size))
    return out


@dataclass(frozen=True, eq=False)
class MdsStage:
    """The range's genes and the MDS of the isolates that have any of them."""

    low: float
    high: float
    gene_idx: np.ndarray  # columns of the slice in range
    kept: np.ndarray  # rows of the slice kept (any gene in range)
    dropped: Tuple[str, ...]
    coords: np.ndarray  # len(kept) x 3
    explained: np.ndarray
    message: str
    seconds: Dict[str, float] = field(default_factory=dict)


@dataclass(frozen=True, eq=False)
class ClusteringResult:
    """One clustering of a slice: what the panel draws and exports."""

    slice: PresenceSlice
    low: float
    high: float
    k: int
    seed: int
    samples: Tuple[str, ...]  # the kept isolates, rows of coords and labels
    dropped: Tuple[str, ...]
    gene_idx: np.ndarray
    rows: np.ndarray  # the slice rows of ``samples``
    coords: np.ndarray
    explained: np.ndarray
    labels: np.ndarray
    message: str = ""
    seconds: Dict[str, float] = field(default_factory=dict)

    @property
    def n_genes(self) -> int:
        return int(self.gene_idx.size)

    @property
    def ok(self) -> bool:
        return not self.message and self.coords.shape[0] > 0

    def range_matrix(self) -> np.ndarray:
        """Kept isolates x genes in range."""
        return self.slice.matrix[np.ix_(self.rows, self.gene_idx)]

    def unique_genes(self) -> List[UniqueGene]:
        if not self.ok:
            return []
        genes = [self.slice.genes[i] for i in self.gene_idx]
        annotations = [self.slice.annotations[i] for i in self.gene_idx]
        return unique_genes(self.range_matrix(), self.labels, genes, annotations)

    def row_order(self) -> List[int]:
        """Slice rows for the band: by cluster, then MDS1; dropped isolates last."""
        if self.ok:
            rows = [int(self.rows[i]) for i in np.lexsort((self.coords[:, 0], self.labels))]
        else:
            rows = [int(r) for r in self.rows]
        taken = set(rows)
        return rows + [i for i in range(len(self.slice.samples)) if i not in taken]


def mds_stage(slice_: PresenceSlice, low: float, high: float) -> MdsStage:
    started = time.perf_counter()
    in_range = np.flatnonzero(genes_in_range(slice_.prevalence, low, high))
    empty = np.zeros((0, MDS_AXES))
    if slice_.n_isolates < MIN_ISOLATES:
        # Every gene of one or two isolates is at 100% or 50%: the range is not the reason.
        count = f"{slice_.n_isolates} isolate" + ("" if slice_.n_isolates == 1 else "s")
        return MdsStage(
            low, high, in_range, np.zeros(0, dtype=np.intp), (), empty, np.zeros(MDS_AXES),
            f"The clustering needs at least {MIN_ISOLATES} isolates with a Panaroo profile; this set has {count}.",
        )
    sub = slice_.matrix[:, in_range]
    has_any = sub.any(axis=1) if in_range.size else np.zeros(slice_.n_isolates, dtype=bool)
    kept = np.flatnonzero(has_any)
    dropped = tuple(slice_.samples[i] for i in np.flatnonzero(~has_any))
    if in_range.size == 0:
        return MdsStage(low, high, in_range, kept[:0], dropped, empty, np.zeros(MDS_AXES), "No genes in this range.")
    if kept.size < 3:
        return MdsStage(
            low, high, in_range, kept, dropped, empty, np.zeros(MDS_AXES),
            f"Fewer than {MIN_ISOLATES} isolates have genes in this range.",
        )
    dist = jaccard_distance(sub[kept])
    t_jaccard = time.perf_counter()
    coords, explained = classical_mds(dist, MDS_AXES)
    del dist
    if coords.shape[1] < MDS_AXES:
        coords = np.hstack([coords, np.zeros((coords.shape[0], MDS_AXES - coords.shape[1]))])
        explained = np.concatenate([explained, np.zeros(MDS_AXES - explained.size)])
    t_mds = time.perf_counter()
    return MdsStage(
        low, high, in_range, kept, dropped, coords, explained, "",
        {"jaccard": t_jaccard - started, "mds": t_mds - t_jaccard},
    )


class ClusteringCache:
    """The MDS of the latest (slice, range) pairs and their k-means labels,
    least recently used out first (thread safe). Distance matrices are never
    kept."""

    def __init__(self, size: int = 8) -> None:
        self._size = int(size)
        self._lock = threading.Lock()
        self._mds: "OrderedDict[Hashable, MdsStage]" = OrderedDict()
        self._labels: "OrderedDict[Hashable, np.ndarray]" = OrderedDict()

    def _put(self, store: "OrderedDict[Hashable, object]", key: Hashable, value: object) -> None:
        store[key] = value
        store.move_to_end(key)
        while len(store) > self._size * (4 if store is self._labels else 1):
            store.popitem(last=False)

    def compute(
        self,
        slice_: PresenceSlice,
        slice_key: Hashable,
        low: float,
        high: float,
        k: int = DEFAULT_K,
        seed: int = DEFAULT_SEED,
    ) -> ClusteringResult:
        key = (slice_key, round(float(low), 6), round(float(high), 6))
        with self._lock:
            stage = self._mds.get(key)
            if stage is not None:
                self._mds.move_to_end(key)
        if stage is None:
            stage = mds_stage(slice_, low, high)
            with self._lock:
                self._put(self._mds, key, stage)
        seconds = dict(stage.seconds)
        label_key = (key, int(k), int(seed))
        with self._lock:
            labels = self._labels.get(label_key)
        if labels is None:
            started = time.perf_counter()
            labels = kmeans(stage.coords, k, seed) if stage.coords.shape[0] else np.zeros(0, dtype=np.int32)
            seconds["kmeans"] = time.perf_counter() - started
            with self._lock:
                self._put(self._labels, label_key, labels)
        return ClusteringResult(
            slice=slice_,
            low=float(low),
            high=float(high),
            k=int(labels.max()) if labels.size else int(k),
            seed=int(seed),
            samples=tuple(slice_.samples[i] for i in stage.kept) if stage.coords.shape[0] else (),
            dropped=stage.dropped,
            gene_idx=stage.gene_idx,
            rows=stage.kept if stage.coords.shape[0] else np.zeros(0, dtype=np.intp),
            coords=stage.coords,
            explained=stage.explained,
            labels=labels,
            message=stage.message,
            seconds=seconds,
        )


def band_image(slice_: PresenceSlice, row_order: Sequence[int], max_cols: int = 1600, max_rows: int = 96) -> np.ndarray:
    """The slice as a small image for the band: rows (isolates in
    ``row_order``) and columns (genes in prevalence order) averaged into at
    most ``max_rows`` x ``max_cols`` cells; values are the share present.
    The genes are binned first, so the full matrix is never copied."""
    x = slice_.matrix
    if x.size == 0:
        return np.zeros((1, 1), dtype=np.float32)

    def starts(n: int, limit: int) -> np.ndarray:
        bins = (np.arange(n) * min(n, limit)) // n
        return np.flatnonzero(np.r_[True, np.diff(bins) > 0])

    rows, cols = x.shape
    r_starts, c_starts = starts(rows, max_rows), starts(cols, max_cols)
    binned = np.add.reduceat(x, c_starts, axis=1, dtype=np.int32)  # isolates x gene bins
    if len(row_order):
        binned = binned[np.asarray(list(row_order), dtype=np.intp)]
    out = np.add.reduceat(binned, r_starts, axis=0, dtype=np.int64).astype(np.float32)
    r_sizes = np.diff(np.r_[r_starts, rows]).astype(np.float32)
    c_sizes = np.diff(np.r_[c_starts, cols]).astype(np.float32)
    out /= r_sizes[:, None] * c_sizes[None, :]
    return out


def export_tables(
    result: ClusteringResult, color_label: str = "", color_by_sample: Optional[Dict[str, str]] = None
) -> Tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """The cluster assignments (dropped isolates with an empty cluster), the
    samples of each cluster and the unique genes of each cluster."""
    colors = color_by_sample or {}
    names = list(result.samples) + list(result.dropped)
    n_kept = len(result.samples)
    coords = result.coords if result.ok else np.zeros((0, MDS_AXES))

    def axis(i: int) -> List[Optional[float]]:
        return [float(coords[r, i]) for r in range(n_kept)] + [None] * len(result.dropped)

    assignments: Dict[str, list] = {
        "Genome_Name": names,
        "Cluster_Number": [int(c) for c in result.labels] + [None] * len(result.dropped)
        if result.ok
        else [None] * len(names),
    }
    if color_label:
        assignments[color_label] = [colors.get(name, "") for name in names]
    for i in range(MDS_AXES):
        assignments[f"MDS{i + 1}"] = axis(i)
    clusters_df = pl.DataFrame(assignments, schema_overrides={"Cluster_Number": pl.Int32})

    members: Dict[int, List[str]] = {}
    if result.ok:
        for name, label in zip(result.samples, result.labels):
            members.setdefault(int(label), []).append(name)
    samples_df = pl.DataFrame(
        {
            "Cluster_Number": sorted(members),
            "N_isolates": [len(members[c]) for c in sorted(members)],
            "Genome_Names": [";".join(members[c]) for c in sorted(members)],
        },
        schema={"Cluster_Number": pl.Int32, "N_isolates": pl.Int32, "Genome_Names": pl.Utf8},
    )
    unique = result.unique_genes()
    unique_df = pl.DataFrame(
        {
            "Cluster_Number": [u.cluster for u in unique],
            "Gene": [u.gene for u in unique],
            "Annotation": [u.annotation for u in unique],
            "N_present_in_cluster": [u.n_present for u in unique],
            "Cluster_size": [u.cluster_size for u in unique],
            "In_every_member": [u.n_present == u.cluster_size for u in unique],
        },
        schema={
            "Cluster_Number": pl.Int32,
            "Gene": pl.Utf8,
            "Annotation": pl.Utf8,
            "N_present_in_cluster": pl.Int32,
            "Cluster_size": pl.Int32,
            "In_every_member": pl.Boolean,
        },
    )
    return clusters_df, samples_df, unique_df
