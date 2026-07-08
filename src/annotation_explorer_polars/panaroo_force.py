from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import html
import math
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np
from PyQt6 import QtCore, QtGui, QtWidgets

from .cgmlst_graph_data import (
    _CHEWBBACA_ALLELES_FILENAME,
    _build_mst_adjacency,
    _build_mst_edges,
    _circle_layout,
    _file_signature,
    _find_chewbbaca_out,
    _load_chewbbaca_distance_cache,
    _load_chewbbaca_profiles_cached,
    _rooted_radial_mst_layout,
    _spread_duplicate_profile_nodes,
)
from .panaroo_widgets import _PgUnavailableWidget
from .panaroo_style import _apply_analysis_panel_style, _log_view_event
from .refresh_worker import RefreshWorker

try:
    import pyqtgraph as pg
except Exception:
    pg = None


@dataclass(frozen=True)
class _CgmlstForcePayload:
    token: int
    request_key: Tuple[object, ...]
    sample_labels: Tuple[str, ...]
    sample_called_loci: Tuple[int, ...]
    initial_coords: np.ndarray
    edges: Tuple[Tuple[int, int, int], ...]
    rest_lengths: Tuple[float, ...]
    root_index: int
    n_loci: int
    mean_pairwise_distance: float
    duplicate_profile_groups: int


@dataclass(frozen=True)
class CgmlstForceConfig:
    charge_strength: float = 0.006
    collision_radius: float = 0.23
    collision_strength: float = 2.35
    spring_strength: float = 2.0
    rest_length_min: float = 0.22
    rest_length_span: float = 3.5
    anchor_strength: float = 0.0
    gravity_strength: float = 0.0
    reference_anchor_strength: float = 0.095
    node_edge_radius: float = 0.5
    node_edge_strength: float = 0.5
    dt_base: float = 0.004
    dt_alpha_scale: float = 0.020
    damping: float = 0.96
    max_speed_base: float = 0.072
    max_speed_alpha_scale: float = 0.046
    cooling: float = 0.999
    stop_alpha: float = 0.001
    stop_speed: float = 0.00005
    max_ticks: int = 460
    node_size: float = 21.0
    focus_node_size: float = 13.0
    hover_node_size: float = 26.0
    reference_node_size: float = 21.0

    def to_dict(self) -> Dict[str, float | int]:
        return dict(asdict(self))


DEFAULT_CGMLST_FORCE_CONFIG = CgmlstForceConfig()


def _normalize_coords(coords: np.ndarray, *, anchor_idx: Optional[int] = None) -> np.ndarray:
    arr = np.array(coords, dtype=float, copy=True)
    if arr.size == 0:
        return arr
    if anchor_idx is not None and 0 <= int(anchor_idx) < int(arr.shape[0]):
        arr -= arr[int(anchor_idx)]
    else:
        arr -= np.mean(arr, axis=0, keepdims=True)
    scale = max(float(np.ptp(arr[:, 0])), float(np.ptp(arr[:, 1])), 1e-6)
    arr /= scale
    if anchor_idx is not None and 0 <= int(anchor_idx) < int(arr.shape[0]):
        arr -= arr[int(anchor_idx)]
    arr.setflags(write=False)
    return arr


def _infer_mst_root_index(
    sample_labels: Sequence[str],
    edges: Sequence[Tuple[int, int, int]],
    reference_sample: Optional[str],
) -> int:
    ref = str(reference_sample or "").strip()
    labels = [str(label).strip() for label in sample_labels]
    if ref and ref in labels:
        return labels.index(ref)
    indegree = np.zeros(len(labels), dtype=np.int32)
    for _left, right, _dist in edges:
        if 0 <= int(right) < len(labels):
            indegree[int(right)] += 1
    roots = np.flatnonzero(indegree == 0)
    if roots.size:
        return int(roots[0])
    return 0


def _rooted_layered_mst_layout(
    sample_labels: Sequence[str],
    edges: Sequence[Tuple[int, int, int]],
    *,
    reference_idx: Optional[int] = None,
) -> np.ndarray:
    n = len(sample_labels)
    if n <= 1:
        return np.zeros((n, 2), dtype=float)

    adjacency = _build_mst_adjacency(n, edges)
    root = int(reference_idx) if reference_idx is not None and 0 <= int(reference_idx) < n else 0
    parent = np.full(n, -1, dtype=np.int32)
    depth = np.zeros(n, dtype=float)
    order: List[int] = [root]
    parent[root] = root

    for node in order:
        neighbors = sorted(
            adjacency[node],
            key=lambda item: (float(item[1]), str(sample_labels[int(item[0])]).casefold()),
        )
        for child, weight in neighbors:
            if child == int(parent[node]) or parent[child] != -1:
                continue
            parent[child] = node
            depth[child] = depth[node] + 1.0 + (0.55 * math.sqrt(max(0.0, float(weight))))
            order.append(child)

    if len(order) != n:
        return _rooted_radial_mst_layout(sample_labels, edges, reference_idx=reference_idx)

    x = np.zeros(n, dtype=float)
    next_leaf = 0.0

    def assign_x(node: int) -> None:
        nonlocal next_leaf
        children = [child for child, _weight in adjacency[node] if parent[child] == node]
        if not children:
            x[node] = next_leaf
            next_leaf += 1.0
            return
        children.sort(key=lambda child: (depth[child], str(sample_labels[child]).casefold()))
        for child in children:
            assign_x(child)
        x[node] = float(np.mean([x[child] for child in children]))

    assign_x(root)
    max_x = max(1.0, float(np.ptp(x)))
    max_depth = max(1.0, float(np.max(depth)))
    coords = np.zeros((n, 2), dtype=float)
    coords[:, 0] = (x - x[root]) / max_x
    coords[:, 1] = -depth / max_depth
    return coords


def _force_simulation_step(
    pos: np.ndarray,
    vel: np.ndarray,
    alpha: float,
    ticks: int,
    cfg: CgmlstForceConfig,
    anchor_positions: np.ndarray,
    edge_left_idx: np.ndarray,
    edge_right_idx: np.ndarray,
    rest_lengths: np.ndarray,
    node_edge_node_idx: np.ndarray,
    node_edge_left_idx: np.ndarray,
    node_edge_right_idx: np.ndarray,
    root_index: int,
) -> Tuple[np.ndarray, np.ndarray, float, int, bool]:
    n = pos.shape[0]
    if n <= 1:
        return pos, vel, 0.0, ticks, False

    forces = np.zeros_like(pos)
    delta = pos[:, np.newaxis, :] - pos[np.newaxis, :, :]
    dist = np.linalg.norm(delta, axis=2)
    np.fill_diagonal(dist, np.inf)
    safe_dist = np.maximum(dist, 1e-4)
    unit = np.divide(
        delta,
        safe_dist[:, :, np.newaxis],
        out=np.zeros_like(delta),
        where=np.isfinite(safe_dist)[:, :, np.newaxis],
    )

    charge_strength = cfg.charge_strength * math.sqrt(max(2, n))
    repulsion = np.where(np.isfinite(dist), charge_strength / np.square(safe_dist), 0.0)
    forces += np.sum(unit * repulsion[:, :, np.newaxis], axis=1)

    overlap = np.maximum(0.0, cfg.collision_radius - dist)
    forces += np.sum(unit * (overlap * cfg.collision_strength)[:, :, np.newaxis], axis=1)

    if edge_left_idx.size and rest_lengths.size:
        edge_vec = pos[edge_right_idx] - pos[edge_left_idx]
        edge_len = np.linalg.norm(edge_vec, axis=1)
        safe_edge_len = np.maximum(edge_len, 1e-6)
        direction = edge_vec / safe_edge_len[:, np.newaxis]
        zero_mask = edge_len <= 1e-6
        if np.any(zero_mask):
            direction[zero_mask] = np.array([1.0, 0.0], dtype=float)
        force_mag = cfg.spring_strength * (safe_edge_len - rest_lengths)
        force_vec = direction * force_mag[:, np.newaxis]
        np.add.at(forces, edge_left_idx, force_vec)
        np.add.at(forces, edge_right_idx, -force_vec)

    if cfg.node_edge_radius > 0.0 and cfg.node_edge_strength > 0.0 and node_edge_node_idx.size:
        threshold = float(cfg.node_edge_radius)
        edge_k = float(cfg.node_edge_strength)
        node_pos = pos[node_edge_node_idx]
        start = pos[node_edge_left_idx]
        end = pos[node_edge_right_idx]
        seg = end - start
        seg_len_sq = np.sum(seg * seg, axis=1)
        seg_safe = np.maximum(seg_len_sq, 1e-12)
        t = np.clip(
            np.sum((node_pos - start) * seg, axis=1) / seg_safe,
            0.0,
            1.0,
        )
        closest = start + (t[:, np.newaxis] * seg)
        away = node_pos - closest
        gap = np.linalg.norm(away, axis=1)
        active = gap < threshold
        if np.any(active):
            away = away[active]
            gap_active = gap[active]
            seg_active = seg[active]
            t_active = t[active]
            node_idx = node_edge_node_idx[active]
            left_idx = node_edge_left_idx[active]
            right_idx = node_edge_right_idx[active]
            direction = np.zeros_like(away)
            nonzero = gap_active > 1e-6
            if np.any(nonzero):
                direction[nonzero] = away[nonzero] / gap_active[nonzero, np.newaxis]
            zero_mask = ~nonzero
            if np.any(zero_mask):
                zero_seg = seg_active[zero_mask]
                perp = np.column_stack([-zero_seg[:, 1], zero_seg[:, 0]])
                perp_norm = np.linalg.norm(perp, axis=1)
                valid_perp = perp_norm > 1e-9
                if np.any(valid_perp):
                    perp[valid_perp] /= perp_norm[valid_perp, np.newaxis]
                if np.any(~valid_perp):
                    perp[~valid_perp] = np.array([1.0, 0.0], dtype=float)
                direction[zero_mask] = perp
            strength = edge_k * (1.0 - (gap_active / threshold))
            push = direction * (1.35 * strength)[:, np.newaxis]
            offset = direction * (0.70 * strength)[:, np.newaxis]
            np.add.at(forces, node_idx, push)
            np.add.at(forces, left_idx, -offset * (1.0 - t_active)[:, np.newaxis])
            np.add.at(forces, right_idx, -offset * t_active[:, np.newaxis])

    if anchor_positions.shape == pos.shape and pos.size and cfg.anchor_strength > 0.0:
        forces += (anchor_positions - pos) * cfg.anchor_strength
    if cfg.gravity_strength > 0.0:
        forces -= pos * cfg.gravity_strength
    if 0 <= root_index < n and cfg.reference_anchor_strength > 0.0:
        forces[root_index] += (anchor_positions[root_index] - pos[root_index]) * cfg.reference_anchor_strength

    dt = (cfg.dt_alpha_scale * alpha) + cfg.dt_base
    vel = (vel + (forces * dt)) * cfg.damping
    speed = np.linalg.norm(vel, axis=1)
    max_speed = cfg.max_speed_base + (cfg.max_speed_alpha_scale * alpha)
    speed_scale = np.minimum(1.0, max_speed / np.maximum(speed, 1e-9))
    vel *= speed_scale[:, np.newaxis]
    pos = pos + vel

    ticks += 1
    alpha *= cfg.cooling
    mean_speed = float(np.mean(np.linalg.norm(vel, axis=1))) if vel.size else 0.0
    keep_running = not (
        alpha <= cfg.stop_alpha
        or mean_speed <= cfg.stop_speed
        or ticks >= int(cfg.max_ticks)
    )
    if not keep_running:
        alpha = 0.0
    return pos, vel, alpha, ticks, keep_running


class _ForceSimulationWorker(QtCore.QObject):
    frameReady = QtCore.pyqtSignal(int, object, int, float)
    finished = QtCore.pyqtSignal(int, object, int)

    def __init__(
        self,
        token: int,
        positions: np.ndarray,
        velocities: np.ndarray,
        alpha: float,
        cfg: CgmlstForceConfig,
        anchor_positions: np.ndarray,
        edge_left_idx: np.ndarray,
        edge_right_idx: np.ndarray,
        rest_lengths: np.ndarray,
        node_edge_node_idx: np.ndarray,
        node_edge_left_idx: np.ndarray,
        node_edge_right_idx: np.ndarray,
        root_index: int,
        *,
        substeps: int,
        frame_interval_ms: int = 33,
    ) -> None:
        super().__init__()
        self._token = int(token)
        self._positions = np.array(positions, dtype=float, copy=True)
        self._velocities = np.array(velocities, dtype=float, copy=True)
        self._alpha = float(alpha)
        self._cfg = cfg
        self._anchor_positions = np.array(anchor_positions, dtype=float, copy=True)
        self._edge_left_idx = np.array(edge_left_idx, dtype=np.int32, copy=True)
        self._edge_right_idx = np.array(edge_right_idx, dtype=np.int32, copy=True)
        self._rest_lengths = np.array(rest_lengths, dtype=float, copy=True)
        self._node_edge_node_idx = np.array(node_edge_node_idx, dtype=np.int32, copy=True)
        self._node_edge_left_idx = np.array(node_edge_left_idx, dtype=np.int32, copy=True)
        self._node_edge_right_idx = np.array(node_edge_right_idx, dtype=np.int32, copy=True)
        self._root_index = int(root_index)
        self._substeps = max(1, int(substeps))
        self._frame_interval_ms = max(1, int(frame_interval_ms))
        self._ticks = 0
        self._stop_requested = False

    def stop(self) -> None:
        self._stop_requested = True

    @QtCore.pyqtSlot()
    def run(self) -> None:
        keep_running = self._positions.shape[0] > 1 and self._alpha > 0.0
        while keep_running and not self._stop_requested:
            for _ in range(self._substeps):
                self._positions, self._velocities, self._alpha, self._ticks, keep_running = _force_simulation_step(
                    self._positions,
                    self._velocities,
                    self._alpha,
                    self._ticks,
                    self._cfg,
                    self._anchor_positions,
                    self._edge_left_idx,
                    self._edge_right_idx,
                    self._rest_lengths,
                    self._node_edge_node_idx,
                    self._node_edge_left_idx,
                    self._node_edge_right_idx,
                    self._root_index,
                )
                if not keep_running or self._stop_requested:
                    break
            self.frameReady.emit(self._token, self._positions.copy(), int(self._ticks), float(self._alpha))
            if keep_running and not self._stop_requested:
                QtCore.QThread.msleep(self._frame_interval_ms)
        self.finished.emit(self._token, self._positions.copy(), int(self._ticks))


def _build_cgmlst_force_payload(
    token: int,
    request_key: Tuple[object, ...],
    output_dir: Path,
    sample_filter: Optional[Tuple[str, ...]],
    reference_sample: Optional[str],
) -> _CgmlstForcePayload:
    out = _find_chewbbaca_out(output_dir)
    if out is None:
        raise FileNotFoundError("chewBBACA output not found.")
    alleles_path = out / _CHEWBBACA_ALLELES_FILENAME
    sig = _file_signature(alleles_path)
    _all_profile_labels, all_values, all_valid = _load_chewbbaca_profiles_cached(*sig)
    all_labels, full_dist, full_score, n_loci = _load_chewbbaca_distance_cache(*sig)

    _selected = set(sample_filter or ())
    keep_idx = list(range(len(all_labels)))
    if not keep_idx:
        raise ValueError("No selected samples available for cgMLST force graph view.")

    sample_labels = tuple(str(all_labels[idx]).strip() for idx in keep_idx)
    dist = full_dist[np.ix_(keep_idx, keep_idx)].astype(float, copy=False)
    score = full_score[np.ix_(keep_idx, keep_idx)]
    values = all_values[keep_idx]
    valid = all_valid[keep_idx]
    sample_called_loci = tuple(int(np.count_nonzero(valid[idx])) for idx in range(valid.shape[0]))
    edges = tuple(
        _build_mst_edges(
            score,
            dist.astype(np.int32, copy=False),
            sample_labels,
            reference_sample=reference_sample,
        )
    )
    if len(sample_labels) <= 1:
        initial_coords = np.zeros((len(sample_labels), 2), dtype=float)
        initial_coords.setflags(write=False)
        rest_lengths: Tuple[float, ...] = tuple()
    else:
        root_index = _infer_mst_root_index(sample_labels, edges, reference_sample)
        initial_coords = _rooted_layered_mst_layout(
            sample_labels,
            edges,
            reference_idx=root_index,
        )
        if not np.isfinite(initial_coords).all() or np.max(np.ptp(initial_coords, axis=0)) <= 1e-8:
            initial_coords = _rooted_radial_mst_layout(
                sample_labels,
                edges,
                reference_idx=root_index,
            )
        if not np.isfinite(initial_coords).all() or np.max(np.ptp(initial_coords, axis=0)) <= 1e-8:
            initial_coords = _circle_layout(len(sample_labels))
        initial_coords = _spread_duplicate_profile_nodes(initial_coords, sample_labels, values, valid)
        initial_coords = _normalize_coords(initial_coords, anchor_idx=root_index)
        rest_lengths = tuple()
    if len(sample_labels) <= 1:
        root_index = 0

    upper_tri = dist[np.triu_indices(dist.shape[0], k=1)]
    finite_upper = upper_tri[np.isfinite(upper_tri)]
    mean_pairwise_distance = float(np.mean(finite_upper)) if finite_upper.size else 0.0

    profile_counts: Dict[Tuple[Tuple[int, ...], Tuple[bool, ...]], int] = {}
    for idx, _label in enumerate(sample_labels):
        key = (
            tuple(int(value) for value in values[idx].tolist()),
            tuple(bool(value) for value in valid[idx].tolist()),
        )
        profile_counts[key] = profile_counts.get(key, 0) + 1
    duplicate_profile_groups = sum(1 for count in profile_counts.values() if count > 1)

    return _CgmlstForcePayload(
        token=int(token),
        request_key=request_key,
        sample_labels=sample_labels,
        sample_called_loci=sample_called_loci,
        initial_coords=initial_coords,
        edges=edges,
        rest_lengths=rest_lengths,
        root_index=int(root_index),
        n_loci=int(n_loci),
        mean_pairwise_distance=mean_pairwise_distance,
        duplicate_profile_groups=int(duplicate_profile_groups),
    )


class CgmlstForceGraphView(QtWidgets.QWidget):
    sampleHovered = QtCore.pyqtSignal(object)
    _MAX_NODE_EDGE_CONSTRAINTS = 180_000

    def __init__(
        self,
        output_dir: Path,
        *,
        reference_sample: Optional[str] = None,
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        _apply_analysis_panel_style(self)
        self._output_dir = Path(output_dir)
        _log_view_event(f"created cgMLST force graph view for {self._output_dir}")
        self._reference_sample = str(reference_sample or "").strip()
        self._sample_filter: Optional[set[str]] = None
        self._highlight_sample: Optional[str] = None
        self._highlight_samples: Set[str] = set()
        self._hover_index: Optional[int] = None
        self._lineage_by_sample: Dict[str, str] = {}
        self._lineage_colors: Dict[str, QtGui.QColor] = {}
        self._cluster_by_sample: Dict[str, str] = {}
        self._custom_color_by_sample: Dict[str, QtGui.QColor] = {}
        self._symbol_by_sample: Dict[str, str] = {}
        self._sample_hover_metadata: Dict[str, Dict[str, str]] = {}
        self._sample_labels: List[str] = []
        self._coords_by_index: Dict[int, Tuple[float, float]] = {}
        self._point_by_index: Dict[int, object] = {}
        self._base_color_by_index: Dict[int, QtGui.QColor] = {}
        self._sample_by_index: Dict[int, str] = {}
        self._called_loci_by_index: Dict[int, int] = {}
        self._sample_called_loci: Tuple[int, ...] = tuple()
        self._edges: List[Tuple[int, int, int]] = []
        self._rest_lengths = np.zeros(0, dtype=float)
        self._edge_left_idx = np.zeros(0, dtype=np.int32)
        self._edge_right_idx = np.zeros(0, dtype=np.int32)
        self._node_edge_node_idx = np.zeros(0, dtype=np.int32)
        self._node_edge_left_idx = np.zeros(0, dtype=np.int32)
        self._node_edge_right_idx = np.zeros(0, dtype=np.int32)
        self._config = DEFAULT_CGMLST_FORCE_CONFIG
        self._refresh_token = 0
        self._active_refresh_token: Optional[int] = None
        self._worker_refs: Dict[int, RefreshWorker] = {}
        self._last_request_key: Optional[Tuple[object, ...]] = None
        self._last_payload: Optional[_CgmlstForcePayload] = None
        self._needs_refresh = True
        self._anchor_positions = np.zeros((0, 2), dtype=float)
        self._positions = np.zeros((0, 2), dtype=float)
        self._velocities = np.zeros((0, 2), dtype=float)
        self._root_index = 0
        self._simulation_alpha = 0.0
        self._simulation_ticks = 0
        self._simulation_token = 0
        self._simulation_thread: Optional[QtCore.QThread] = None
        self._simulation_worker: Optional[_ForceSimulationWorker] = None
        self._positions_revision = 0
        self._point_style_key_by_index: Dict[int, Tuple[object, ...]] = {}
        self._brush_cache: Dict[Tuple[int, int, int, int], QtGui.QBrush] = {}
        self._pen_cache: Dict[Tuple[int, int, int, int, float], QtGui.QPen] = {}
        self._last_edge_highlight_key: Optional[Tuple[Tuple[str, ...], int]] = None
        self._scatter_item: Optional["pg.ScatterPlotItem"] = None
        self._edge_item: Optional["pg.PlotDataItem"] = None
        self._edge_highlight_item: Optional["pg.PlotDataItem"] = None
        self._hover_text: Optional["pg.TextItem"] = None

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self._info = QtWidgets.QLabel()
        self._info.setWordWrap(True)
        self._info.setStyleSheet("color: #607183; font-weight: 600; padding: 6px 8px;")
        self._info.setVisible(False)
        layout.addWidget(self._info)

        if pg is None:
            layout.addWidget(_PgUnavailableWidget("PyQtGraph is required for the cgMLST force-graph view."))
            return

        self._plot = pg.PlotWidget()
        self._plot.setBackground("#ffffff")
        self._plot.getPlotItem().setMenuEnabled(False)
        self._plot.hideAxis("left")
        self._plot.hideAxis("bottom")
        self._plot.setAspectLocked(True, ratio=1.0)
        layout.addWidget(self._plot, stretch=1)
        self._hover_proxy = pg.SignalProxy(
            self._plot.scene().sigMouseMoved,
            rateLimit=60,
            slot=self._on_mouse_moved,
        )
        self.destroyed.connect(lambda _obj=None: self._stop_simulation(wait=True))

    def force_config(self) -> CgmlstForceConfig:
        return replace(self._config)

    def set_force_config(self, config: CgmlstForceConfig, *, restart: bool = True) -> None:
        if config == self._config:
            return
        self._config = config
        self._rebuild_rest_lengths()
        if restart and self._last_payload is not None:
            self._apply_payload(self._last_payload)
            self._update_info_text(self._last_payload, running=self._simulation_worker is not None)
        else:
            self._apply_point_styles()
            self._apply_edge_highlight()

    def restart_simulation(self) -> None:
        if self._last_payload is None:
            return
        self._apply_payload(self._last_payload)
        self._update_info_text(self._last_payload, running=self._simulation_worker is not None)

    def set_reference_sample(self, sample: Optional[str], *, refresh: bool = True) -> None:
        value = str(sample or "").strip()
        if value == self._reference_sample:
            return
        self._reference_sample = value
        self._needs_refresh = True
        if refresh and self.isVisible():
            self.refresh()
        else:
            self._apply_point_styles()
            self._apply_edge_highlight()

    def set_sample_filter(self, samples: Optional[Sequence[str]], *, refresh: bool = True) -> None:
        next_filter = {str(s).strip() for s in samples if str(s).strip()} if samples is not None else None
        if next_filter == self._sample_filter:
            return
        self._sample_filter = next_filter
        self._apply_point_styles()
        self._apply_edge_highlight()

    def set_fastgear_assignments(
        self,
        lineage_by_sample: Mapping[str, str],
        lineage_colors: Mapping[str, QtGui.QColor],
        cluster_by_sample: Mapping[str, str],
    ) -> None:
        self._lineage_by_sample = {str(k).strip(): str(v).strip() for k, v in lineage_by_sample.items() if str(k).strip()}
        self._lineage_colors = {
            str(k).strip(): QtGui.QColor(v)
            for k, v in lineage_colors.items()
            if str(k).strip()
        }
        self._cluster_by_sample = {str(k).strip(): str(v).strip() for k, v in cluster_by_sample.items() if str(k).strip()}
        if self._last_payload is not None:
            self._render_scene()
        else:
            self._apply_point_styles()
            self._apply_edge_highlight()

    def set_node_style_maps(
        self,
        *,
        color_by_sample: Optional[Mapping[str, object]] = None,
        symbol_by_sample: Optional[Mapping[str, object]] = None,
    ) -> None:
        normalized_colors: Dict[str, QtGui.QColor] = {}
        for raw_sample, raw_color in (color_by_sample or {}).items():
            sample = str(raw_sample).strip()
            color = raw_color if isinstance(raw_color, QtGui.QColor) else QtGui.QColor(raw_color)
            if not sample or not color.isValid():
                continue
            normalized_colors[sample] = QtGui.QColor(color)

        normalized_symbols: Dict[str, str] = {}
        for raw_sample, raw_symbol in (symbol_by_sample or {}).items():
            sample = str(raw_sample).strip()
            symbol = str(raw_symbol).strip()
            if not sample or not symbol:
                continue
            normalized_symbols[sample] = symbol

        self._custom_color_by_sample = normalized_colors
        self._symbol_by_sample = normalized_symbols
        if self._last_payload is not None:
            self._render_scene()
        else:
            self._apply_point_styles()
            self._apply_edge_highlight()

    def set_sample_hover_metadata(self, metadata: Optional[Mapping[str, Mapping[str, object]]]) -> None:
        normalized: Dict[str, Dict[str, str]] = {}
        for raw_sample, raw_fields in (metadata or {}).items():
            sample = str(raw_sample).strip()
            if not sample or not isinstance(raw_fields, Mapping):
                continue
            fields: Dict[str, str] = {}
            for raw_key, raw_value in raw_fields.items():
                key = str(raw_key).strip()
                value = str(raw_value).strip() if raw_value is not None else ""
                if not key or not value:
                    continue
                fields[key] = value
            if fields:
                normalized[sample] = fields
        if normalized == self._sample_hover_metadata:
            return
        self._sample_hover_metadata = normalized
        if self._hover_index is not None and self._hover_text is not None:
            self._update_hover_text(self._hover_index)

    def ensure_fresh(self) -> None:
        if self._needs_refresh:
            self.refresh()

    def highlight_sample(self, sample: Optional[str]) -> None:
        value = str(sample).strip() if sample is not None else ""
        if self._highlight_sample == (value or None) and not self._highlight_samples:
            return
        self._highlight_sample = value or None
        self._highlight_samples = set()
        self._apply_point_styles()
        self._apply_edge_highlight()

    def highlight_samples(self, samples: Optional[Sequence[str]]) -> None:
        normalized = {str(s).strip() for s in (samples or []) if str(s).strip()}
        if self._highlight_sample is None and normalized == self._highlight_samples:
            return
        self._highlight_samples = normalized
        self._highlight_sample = None
        self._apply_point_styles()
        self._apply_edge_highlight()

    def _current_request_key(self) -> Tuple[object, ...]:
        out = _find_chewbbaca_out(self._output_dir)
        if out is None:
            return ("missing",)
        alleles_path = out / _CHEWBBACA_ALLELES_FILENAME
        if not alleles_path.is_file():
            return ("missing",)
        sig = _file_signature(alleles_path)
        return ("ok", *sig, str(self._reference_sample or "").strip())

    def _apply_payload(self, payload: _CgmlstForcePayload) -> None:
        if pg is None:
            return
        self._stop_simulation(wait=True)
        self._plot.clear()
        self._scatter_item = None
        self._edge_item = None
        self._edge_highlight_item = None
        self._hover_text = None
        self._hover_index = None
        self._sample_labels = list(payload.sample_labels)
        self._coords_by_index = {}
        self._point_by_index = {}
        self._base_color_by_index = {}
        self._sample_by_index = {}
        self._called_loci_by_index = {}
        self._sample_called_loci = tuple(int(value) for value in payload.sample_called_loci)
        self._edges = list(payload.edges)
        self._root_index = int(payload.root_index) if 0 <= int(payload.root_index) < len(self._sample_labels) else 0
        self._anchor_positions = np.array(payload.initial_coords, dtype=float, copy=True)
        self._positions = np.array(payload.initial_coords, dtype=float, copy=True)
        self._velocities = np.zeros_like(self._positions)
        self._rebuild_rest_lengths()
        self._rebuild_topology_cache()
        self._simulation_alpha = 0.38 if self._positions.shape[0] > 1 else 0.0
        self._simulation_ticks = 0
        self._positions_revision = 0
        self._point_style_key_by_index = {}
        self._last_edge_highlight_key = None

        self._edge_item = pg.PlotDataItem([], [], pen=pg.mkPen(QtGui.QColor(148, 163, 184, 170), width=1.3), connect="finite")
        self._plot.addItem(self._edge_item)
        self._edge_highlight_item = pg.PlotDataItem([], [], pen=pg.mkPen(QtGui.QColor(15, 23, 42, 235), width=2.4), connect="finite")
        self._plot.addItem(self._edge_highlight_item)
        self._scatter_item = pg.ScatterPlotItem()
        self._plot.addItem(self._scatter_item)
        self._hover_text = pg.TextItem(anchor=(0.0, 1.4))
        self._hover_text.setZValue(30)
        self._hover_text.setVisible(False)
        self._plot.addItem(self._hover_text)

        self._render_scene()
        self._set_initial_view_range()
        if self._simulation_alpha > 0.0:
            self._start_simulation()
            self._update_info_text(payload, running=True)
        else:
            self._update_info_text(payload, running=False)

    def _update_info_text(self, payload: Optional[_CgmlstForcePayload], *, running: bool) -> None:
        if payload is None:
            return
        lineage_count = len({value for value in self._lineage_by_sample.values() if value})
        status = "running" if running else "settled"
        _log_view_event(
            f"cgMLST force graph | samples: {len(payload.sample_labels)} | loci: {payload.n_loci} | "
            f"MST edges: {len(payload.edges)} | mean pairwise allele distance: {payload.mean_pairwise_distance:.2f} | "
            f"fastGEAR clades: {lineage_count} | duplicate cgMLST profiles: {payload.duplicate_profile_groups} | "
            f"simulation: {status}"
        )
        self._info.clear()
        self._info.setVisible(False)

    def _set_initial_view_range(self) -> None:
        if pg is None or self._positions.size == 0:
            return
        vb = self._plot.getViewBox()
        if vb is None:
            return
        coords = np.asarray(self._positions, dtype=float)
        min_x = float(np.min(coords[:, 0]))
        max_x = float(np.max(coords[:, 0]))
        min_y = float(np.min(coords[:, 1]))
        max_y = float(np.max(coords[:, 1]))
        center_x = 0.5 * (min_x + max_x)
        center_y = 0.5 * (min_y + max_y)
        half_extent = max(
            1.6,
            1.9 * max(
                0.5 * (max_x - min_x),
                0.5 * (max_y - min_y),
            ) + 0.35,
        )
        vb.disableAutoRange()
        vb.setRange(
            xRange=(center_x - half_extent, center_x + half_extent),
            yRange=(center_y - half_extent, center_y + half_extent),
            padding=0.0,
        )

    def _render_scene(self) -> None:
        if pg is None or self._edge_item is None or self._edge_highlight_item is None or self._scatter_item is None:
            return
        self._coords_by_index = {}
        self._base_color_by_index = {}
        self._sample_by_index = {}
        self._positions_revision += 1
        self._last_edge_highlight_key = None
        edge_x: List[float] = []
        edge_y: List[float] = []
        for left, right, _dist in self._edges:
            if left >= self._positions.shape[0] or right >= self._positions.shape[0]:
                continue
            edge_x.extend([float(self._positions[left, 0]), float(self._positions[right, 0]), np.nan])
            edge_y.extend([float(self._positions[left, 1]), float(self._positions[right, 1]), np.nan])
        self._edge_item.setData(edge_x, edge_y)

        spots = []
        style_keys: Dict[int, Tuple[object, ...]] = {}
        active_filter = self._sample_filter
        emphasis = self._emphasis_samples()
        for idx, sample in enumerate(self._sample_labels):
            x = float(self._positions[idx, 0])
            y = float(self._positions[idx, 1])
            base = self._base_color_for_sample(sample)
            brush, pen, size, style_key = self._point_style_components(
                idx,
                sample,
                base,
                emphasis,
                active_filter,
            )
            spots.append(
                {
                    "pos": (x, y),
                    "data": {"index": int(idx), "sample": sample},
                    "brush": brush,
                    "pen": pen,
                    "size": size,
                    "symbol": self._symbol_by_sample.get(sample, "o"),
                }
            )
            style_keys[idx] = style_key
            self._coords_by_index[idx] = (x, y)
            self._base_color_by_index[idx] = QtGui.QColor(base)
            self._sample_by_index[idx] = sample
            self._called_loci_by_index[idx] = self._sample_called_loci[idx] if idx < len(self._sample_called_loci) else 0
        self._scatter_item.setData(spots)
        self._point_by_index = {}
        self._point_style_key_by_index = style_keys
        for point in self._scatter_item.points():
            data = point.data()
            if not isinstance(data, Mapping):
                continue
            idx = data.get("index")
            if idx is None:
                continue
            try:
                self._point_by_index[int(idx)] = point
            except Exception:
                continue

        if self._hover_index is not None and self._hover_text is not None:
            coords = self._coords_by_index.get(self._hover_index)
            if coords is not None:
                self._hover_text.setPos(float(coords[0]), float(coords[1]))

        self._apply_edge_highlight()

    def _update_hover_text(self, idx: int) -> None:
        if self._hover_text is None:
            return
        sample = self._sample_by_index.get(idx, "")
        if not sample:
            self._hover_text.setVisible(False)
            return
        coords = self._coords_by_index.get(idx, (0.0, 0.0))
        lineage = self._lineage_by_sample.get(sample, "")
        cluster = self._cluster_by_sample.get(sample, "")
        lines = [html.escape(sample)]
        if lineage:
            lines.append(f"<span style='color:#475569'>Clade: {html.escape(lineage)}</span>")
        if cluster:
            lines.append(f"<span style='color:#475569'>Cluster: {html.escape(cluster)}</span>")
        for label, value in self._sample_hover_metadata.get(sample, {}).items():
            lines.append(
                f"<span style='color:#475569'>{html.escape(label)}: {html.escape(value)}</span>"
            )
        called_loci = self._called_loci_by_index.get(idx)
        total_loci = int(self._last_payload.n_loci) if self._last_payload is not None else 0
        if called_loci is not None and total_loci > 0:
            pct = (100.0 * float(called_loci)) / float(total_loci)
            lines.append(
                f"<span style='color:#0f172a'>cgMLST called loci: {pct:.1f}% ({called_loci}/{total_loci})</span>"
            )
        self._hover_text.setHtml("<br>".join(lines))
        self._hover_text.setPos(float(coords[0]), float(coords[1]))
        self._hover_text.setVisible(True)

    def _start_simulation(self) -> None:
        if self._positions.shape[0] <= 1 or self._simulation_alpha <= 0.0:
            return
        self._stop_simulation(wait=True)
        self._simulation_token += 1
        token = int(self._simulation_token)
        substeps = 4 if self._positions.shape[0] <= 256 else 2
        thread = QtCore.QThread(self)
        worker = _ForceSimulationWorker(
            token,
            self._positions,
            self._velocities,
            self._simulation_alpha,
            self._config,
            self._anchor_positions,
            self._edge_left_idx,
            self._edge_right_idx,
            self._rest_lengths,
            self._node_edge_node_idx,
            self._node_edge_left_idx,
            self._node_edge_right_idx,
            self._root_index,
            substeps=substeps,
            frame_interval_ms=33,
        )
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.frameReady.connect(self._on_simulation_frame)
        worker.finished.connect(self._on_simulation_finished)
        worker.finished.connect(thread.quit)
        worker.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        self._simulation_thread = thread
        self._simulation_worker = worker
        thread.start()

    def _stop_simulation(self, *, wait: bool) -> None:
        worker = self._simulation_worker
        thread = self._simulation_thread
        self._simulation_worker = None
        self._simulation_thread = None
        if worker is not None:
            self._simulation_token += 1
            worker.stop()
        if wait and thread is not None and thread.isRunning():
            thread.quit()
            thread.wait(1000)

    def _on_simulation_frame(self, token: int, positions_obj: object, ticks: int, alpha: float) -> None:
        if token != self._simulation_token or self._simulation_worker is None:
            return
        self._positions = np.asarray(positions_obj, dtype=float)
        self._simulation_ticks = int(ticks)
        self._simulation_alpha = float(alpha)
        sample_count = len(self._sample_labels)
        if sample_count > 1200:
            render_stride = 4
        elif sample_count > 600:
            render_stride = 3
        elif sample_count > 300:
            render_stride = 2
        else:
            render_stride = 1
        if render_stride > 1 and (self._simulation_ticks % render_stride) != 0:
            return
        self._render_scene()

    def _on_simulation_finished(self, token: int, positions_obj: object, ticks: int) -> None:
        if token != self._simulation_token:
            return
        self._positions = np.asarray(positions_obj, dtype=float)
        self._simulation_ticks = int(ticks)
        self._simulation_alpha = 0.0
        self._simulation_worker = None
        self._simulation_thread = None
        self._render_scene()
        self._update_info_text(self._last_payload, running=False)

    def _rebuild_rest_lengths(self) -> None:
        if not self._edges:
            self._rest_lengths = np.zeros(0, dtype=float)
            return
        max_edge = max((int(distance) for _left, _right, distance in self._edges), default=1)
        max_edge = max(max_edge, 1)
        self._rest_lengths = np.asarray(
            [
                self._config.rest_length_min
                + (
                    self._config.rest_length_span
                    * (max(0.0, float(distance)) / float(max_edge))
                )
                for _left, _right, distance in self._edges
            ],
            dtype=float,
        )

    def _rebuild_topology_cache(self) -> None:
        n = int(len(self._sample_labels))
        if not self._edges or n <= 1:
            self._edge_left_idx = np.zeros(0, dtype=np.int32)
            self._edge_right_idx = np.zeros(0, dtype=np.int32)
            self._node_edge_node_idx = np.zeros(0, dtype=np.int32)
            self._node_edge_left_idx = np.zeros(0, dtype=np.int32)
            self._node_edge_right_idx = np.zeros(0, dtype=np.int32)
            return
        edge_left = np.asarray([int(left) for left, _right, _dist in self._edges], dtype=np.int32)
        edge_right = np.asarray([int(right) for _left, right, _dist in self._edges], dtype=np.int32)
        self._edge_left_idx = edge_left
        self._edge_right_idx = edge_right

        def node_at_rank_excluding(rank: int, left: int, right: int) -> int:
            first, second = sorted((int(left), int(right)))
            node = int(rank)
            if node >= first:
                node += 1
            if node >= second:
                node += 1
            return node

        available_per_edge = max(0, n - 2)
        edge_count = max(1, len(self._edges))
        total_constraints = available_per_edge * edge_count
        if total_constraints <= self._MAX_NODE_EDGE_CONSTRAINTS:
            ranks_by_edge = None
        else:
            ranks_per_edge = max(1, self._MAX_NODE_EDGE_CONSTRAINTS // edge_count)
            ranks_per_edge = min(available_per_edge, ranks_per_edge)
            ranks_by_edge = np.linspace(
                0,
                available_per_edge - 1,
                num=ranks_per_edge,
                dtype=np.int32,
            )

        node_edge_node: List[int] = []
        node_edge_left: List[int] = []
        node_edge_right: List[int] = []
        for left, right, _dist in self._edges:
            left_int = int(left)
            right_int = int(right)
            if ranks_by_edge is None:
                nodes = (
                    node for node in range(n)
                    if node != left_int and node != right_int
                )
            else:
                nodes = (
                    node_at_rank_excluding(int(rank), left_int, right_int)
                    for rank in ranks_by_edge
                )
            for node in nodes:
                node_edge_node.append(node)
                node_edge_left.append(left_int)
                node_edge_right.append(right_int)
        self._node_edge_node_idx = np.asarray(node_edge_node, dtype=np.int32)
        self._node_edge_left_idx = np.asarray(node_edge_left, dtype=np.int32)
        self._node_edge_right_idx = np.asarray(node_edge_right, dtype=np.int32)

    def _on_refresh_finished(self, payload_obj: object) -> None:
        if not isinstance(payload_obj, _CgmlstForcePayload):
            return
        self._worker_refs.pop(int(payload_obj.token), None)
        if self._active_refresh_token != int(payload_obj.token):
            return
        self._active_refresh_token = None
        self._last_request_key = payload_obj.request_key
        self._last_payload = payload_obj
        self._needs_refresh = False
        self._apply_payload(payload_obj)

    def _on_refresh_failed(self, token: int, error_text: str) -> None:
        self._worker_refs.pop(int(token), None)
        if self._active_refresh_token != int(token):
            return
        self._active_refresh_token = None
        self._needs_refresh = True
        self._stop_simulation(wait=False)
        last_line = str(error_text or "").strip().splitlines()[-1] if str(error_text or "").strip() else "Unknown error"
        self._info.setText(f"Failed to build cgMLST force graph: {last_line}")
        self._info.setVisible(True)

    def _set_hover_index(self, idx: Optional[int]) -> None:
        if idx == self._hover_index:
            return
        self._hover_index = idx
        if idx is None:
            if self._hover_text is not None:
                self._hover_text.setVisible(False)
            self._apply_point_styles()
            self._apply_edge_highlight()
            self.sampleHovered.emit(None)
            return
        sample = self._sample_by_index.get(idx, "")
        if self._hover_text is not None and sample:
            self._update_hover_text(idx)
        self._apply_point_styles()
        self._apply_edge_highlight()
        self.sampleHovered.emit(sample or None)

    def _on_mouse_moved(self, evt: object) -> None:
        if pg is None or self._scatter_item is None:
            return
        try:
            pos = evt[0] if isinstance(evt, (tuple, list)) else evt
        except Exception:
            pos = evt
        if pos is None or not self._plot.sceneBoundingRect().contains(pos):
            self._set_hover_index(None)
            return
        vb = self._plot.getViewBox()
        if vb is None:
            self._set_hover_index(None)
            return
        mouse_point = vb.mapSceneToView(pos)
        points = self._scatter_item.pointsAt(mouse_point)
        if points is None:
            self._set_hover_index(None)
            return
        try:
            points_count = len(points)
        except Exception:
            points_count = int(np.asarray(points).size)
        if points_count == 0:
            self._set_hover_index(None)
            return
        best_idx: Optional[int] = None
        best_d2 = float("inf")
        mx = float(mouse_point.x())
        my = float(mouse_point.y())
        for point in points:
            data = point.data()
            if not isinstance(data, Mapping):
                continue
            idx = data.get("index")
            if idx is None:
                continue
            try:
                idx_int = int(idx)
            except Exception:
                continue
            xy = self._coords_by_index.get(idx_int)
            if xy is None:
                continue
            dx = float(xy[0]) - mx
            dy = float(xy[1]) - my
            d2 = (dx * dx) + (dy * dy)
            if d2 < best_d2:
                best_d2 = d2
                best_idx = idx_int
        self._set_hover_index(best_idx)

    def _cached_brush(self, color: QtGui.QColor) -> QtGui.QBrush:
        key = (int(color.red()), int(color.green()), int(color.blue()), int(color.alpha()))
        brush = self._brush_cache.get(key)
        if brush is None:
            brush = pg.mkBrush(*key)
            self._brush_cache[key] = brush
        return brush

    def _cached_pen(self, color: QtGui.QColor, width: float) -> QtGui.QPen:
        rounded_width = round(float(width), 3)
        key = (
            int(color.red()),
            int(color.green()),
            int(color.blue()),
            int(color.alpha()),
            rounded_width,
        )
        pen = self._pen_cache.get(key)
        if pen is None:
            pen = pg.mkPen(QtGui.QColor(key[0], key[1], key[2], key[3]), width=rounded_width)
            self._pen_cache[key] = pen
        return pen

    def _base_color_for_sample(self, sample: str) -> QtGui.QColor:
        custom = self._custom_color_by_sample.get(sample)
        if custom is not None and custom.isValid():
            return QtGui.QColor(custom)
        lineage = self._lineage_by_sample.get(sample, "")
        return QtGui.QColor(self._lineage_colors.get(lineage, QtGui.QColor("#2563eb")))

    def _emphasis_samples(self) -> Set[str]:
        emphasis = set(self._highlight_samples)
        if self._highlight_sample:
            emphasis.add(self._highlight_sample)
        hover_sample = self._sample_by_index.get(self._hover_index or -1, "")
        if hover_sample:
            emphasis.add(hover_sample)
        return emphasis

    def _point_style_components(
        self,
        idx: int,
        sample: str,
        base: QtGui.QColor,
        emphasis: Set[str],
        active_filter: Optional[set[str]],
    ) -> Tuple[QtGui.QBrush, QtGui.QPen, float, Tuple[object, ...]]:
        in_filter = active_filter is None or sample in active_filter
        dim_others = bool(emphasis)
        active = in_filter and (not dim_others or sample in emphasis)
        if active:
            alpha = 215
        elif in_filter:
            alpha = 72
        else:
            alpha = 24
        brush_color = QtGui.QColor(base)
        brush_color.setAlpha(alpha)
        is_hover = self._hover_index == idx
        is_focus = in_filter and sample in emphasis
        if is_hover:
            size = self._config.hover_node_size
        else:
            size = self._config.node_size
        pen_color = QtGui.QColor("#0f172a" if (is_hover or is_focus) else "#334155")
        if not in_filter:
            pen_color.setAlpha(60)
        pen_width = 2.4 if is_hover else (2.0 if is_focus else 1.0)
        style_key = (
            int(brush_color.red()),
            int(brush_color.green()),
            int(brush_color.blue()),
            int(brush_color.alpha()),
            round(float(size), 3),
            int(pen_color.red()),
            int(pen_color.green()),
            int(pen_color.blue()),
            int(pen_color.alpha()),
            round(float(pen_width), 3),
        )
        return (
            self._cached_brush(brush_color),
            self._cached_pen(pen_color, pen_width),
            float(size),
            style_key,
        )

    def _apply_point_styles(self) -> None:
        if pg is None:
            return
        active_filter = self._sample_filter
        emphasis = self._emphasis_samples()
        for idx, point in self._point_by_index.items():
            base = self._base_color_by_index.get(idx, QtGui.QColor("#475569"))
            sample = self._sample_by_index.get(idx, "")
            brush, pen, size, style_key = self._point_style_components(
                idx,
                sample,
                base,
                emphasis,
                active_filter,
            )
            if self._point_style_key_by_index.get(idx) == style_key:
                continue
            self._point_style_key_by_index[idx] = style_key
            point.setBrush(brush)
            point.setPen(pen)
            point.setSize(size)

    def _apply_edge_highlight(self) -> None:
        if pg is None or self._edge_highlight_item is None:
            return
        emphasis = set(self._highlight_samples)
        if self._highlight_sample:
            emphasis.add(self._highlight_sample)
        hover_sample = self._sample_by_index.get(self._hover_index or -1, "")
        if hover_sample:
            emphasis.add(hover_sample)
        edge_key = (tuple(sorted(emphasis)), int(self._positions_revision))
        if self._last_edge_highlight_key == edge_key:
            return
        self._last_edge_highlight_key = edge_key
        if not emphasis:
            self._edge_highlight_item.setData([], [])
            return
        emphasized_idx = {idx for idx, sample in self._sample_by_index.items() if sample in emphasis}
        xs: List[float] = []
        ys: List[float] = []
        for left, right, _dist in self._edges:
            if left not in emphasized_idx and right not in emphasized_idx:
                continue
            if left >= self._positions.shape[0] or right >= self._positions.shape[0]:
                continue
            xs.extend([float(self._positions[left, 0]), float(self._positions[right, 0]), np.nan])
            ys.extend([float(self._positions[left, 1]), float(self._positions[right, 1]), np.nan])
        self._edge_highlight_item.setData(xs, ys)

    def refresh(self) -> None:
        if pg is None:
            return
        request_key = self._current_request_key()
        if request_key == ("missing",):
            self._info.setText("chewBBACA output not found.")
            self._info.setVisible(True)
            return

        if self._last_request_key == request_key and self._last_payload is not None:
            self._needs_refresh = False
            self._apply_payload(self._last_payload)
            return

        if self._active_refresh_token is not None and self._last_request_key == request_key:
            return

        self._refresh_token += 1
        token = int(self._refresh_token)
        self._active_refresh_token = token
        self._last_request_key = request_key
        self._needs_refresh = False
        self._stop_simulation(wait=False)
        if len(request_key) >= 2 and request_key[0] == "ok":
            print(f"[file] using cgMLST allele calls: {request_key[1]}")
        self._info.setText("Building cgMLST force graph...")
        self._info.setVisible(True)
        worker = RefreshWorker(
            token,
            lambda token=token, request_key=request_key, reference_sample=self._reference_sample:
                _build_cgmlst_force_payload(
                    token,
                    request_key,
                    self._output_dir,
                    None,
                    reference_sample,
                ),
        )
        worker.signals.finished.connect(self._on_refresh_finished)
        worker.signals.failed.connect(self._on_refresh_failed)
        self._worker_refs[token] = worker
        QtCore.QThreadPool.globalInstance().start(worker)
