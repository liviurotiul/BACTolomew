from __future__ import annotations

from dataclasses import dataclass
import html
import math
from pathlib import Path
import time
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np
from PyQt6 import QtCore, QtGui, QtWidgets

from .cgmlst_graph_data import (
    _CHEWBBACA_ALLELES_FILENAME,
    MstEdgeScale,
    _build_mst_edges,
    _file_signature,
    _find_chewbbaca_out,
    _load_chewbbaca_distance_cache,
    _load_chewbbaca_profiles_cached,
    _edge_scale_parts,
    _mst_edge_scale,
)
from .cgmlst_mst_layout import mst_drawing
from .panaroo_widgets import _PgUnavailableWidget
from .panaroo_style import _apply_analysis_panel_style, _log_view_event
from .refresh_worker import RefreshWorker
from .theme import COLORS, OUT_OF_SET_OPACITY, out_of_set_color

# Distance labels and break marks of genomes outside the window's samples:
# faded with them, but a little stronger than the nodes so they stay legible.
_FADED_TEXT_OPACITY = min(1.0, 1.8 * OUT_OF_SET_OPACITY)

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
    # Node positions in node diameters (see cgmlst_mst_layout.mst_drawing).
    coords: np.ndarray
    edges: Tuple[Tuple[int, int, int], ...]
    # How edge lengths follow the allele distance (linear, very long edges
    # broken: drawn short with a break mark).
    scale: MstEdgeScale
    # Per edge, the node its part drawn to scale starts from and that part's
    # share of the drawn edge (below 1 where the layout lengthened the edge
    # to make room; the rest is drawn dotted). See _edge_scale_parts.
    scale_start: np.ndarray
    scale_share: np.ndarray
    n_loci: int
    mean_pairwise_distance: float
    duplicate_profile_groups: int


def _build_cgmlst_force_payload(
    token: int,
    request_key: Tuple[object, ...],
    output_dir: Path,
    reference_sample: Optional[str],
) -> _CgmlstForcePayload:
    out = _find_chewbbaca_out(output_dir)
    if out is None:
        raise FileNotFoundError("chewBBACA output not found.")
    alleles_path = out / _CHEWBBACA_ALLELES_FILENAME
    sig = _file_signature(alleles_path)
    _all_profile_labels, values, valid = _load_chewbbaca_profiles_cached(*sig)
    all_labels, full_dist, score, n_loci = _load_chewbbaca_distance_cache(*sig)
    if not all_labels:
        raise ValueError("No selected samples available for the cgMLST MST view.")

    sample_labels = tuple(str(label).strip() for label in all_labels)
    dist = full_dist.astype(float, copy=False)
    sample_called_loci = tuple(int(np.count_nonzero(valid[idx])) for idx in range(valid.shape[0]))
    edges = tuple(
        _build_mst_edges(
            score,
            dist.astype(np.int32, copy=False),
            sample_labels,
            reference_sample=reference_sample,
        )
    )
    scale = _mst_edge_scale([distance for _left, _right, distance in edges])
    started = time.perf_counter()
    coords, parents, layout = mst_drawing(len(sample_labels), edges, sample_labels, scale)
    _log_view_event(
        f"cgMLST MST layout: {layout}, {len(sample_labels)} nodes, {time.perf_counter() - started:.2f} s"
    )
    coords.setflags(write=False)
    scale_start, scale_share = _edge_scale_parts(coords, edges, parents, scale)

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
        coords=coords,
        edges=edges,
        scale=scale,
        scale_start=scale_start,
        scale_share=scale_share,
        n_loci=int(n_loci),
        mean_pairwise_distance=mean_pairwise_distance,
        duplicate_profile_groups=int(duplicate_profile_groups),
    )


# A label's screen rectangle: left, top, right, bottom.
_Box = Tuple[float, float, float, float]


def _clip_segments(
    x0: np.ndarray, y0: np.ndarray, dx: np.ndarray, dy: np.ndarray, box: _Box
) -> Tuple[np.ndarray, np.ndarray]:
    """The part of each segment (x0 + t dx, y0 + t dy for t in 0..1) inside
    ``box``, as its t range (Liang-Barsky, all at once); t0 > t1 where the
    segment misses the box."""
    t0 = np.zeros(np.size(x0))
    t1 = np.ones(np.size(x0))
    for p, q in ((-dx, x0 - box[0]), (dx, box[2] - x0), (-dy, y0 - box[1]), (dy, box[3] - y0)):
        p = np.asarray(p, dtype=float)
        flat = np.abs(p) < 1e-12
        with np.errstate(divide="ignore", invalid="ignore"):
            t = q / np.where(flat, 1.0, p)
        t0 = np.where(~flat & (p < 0.0), np.maximum(t0, t), t0)
        t1 = np.where(~flat & (p > 0.0), np.minimum(t1, t), t1)
        t1 = np.where(flat & (q < 0.0), -1.0, t1)
    return t0, t1


def _segments_cross_box(
    x0: np.ndarray, y0: np.ndarray, dx: np.ndarray, dy: np.ndarray, box: _Box, *, skip: int = -1
) -> bool:
    """Whether a segment other than number ``skip`` passes through ``box``."""
    t0, t1 = _clip_segments(x0, y0, dx, dy, box)
    hit = t0 <= t1
    if 0 <= skip < hit.size:
        hit[skip] = False
    return bool(hit.any())


class _EdgeDistanceLabels(QtWidgets.QGraphicsObject):
    """Allele distance labels on the force graph's edges, drawn in screen pixels.

    Zoomed in (the median edge long enough on screen), an edge shows its
    distance at its midpoint once the label fits between its two nodes; labels
    fade in as the view zooms in and the gap opens, and are hidden zoomed out.
    The edges of the focused nodes (the hovered node) show theirs at any zoom,
    emphasized: at the midpoint when it is free, else along the edge just
    outside the focused node. No label covers another one or leaves its edge:
    a focused label with no free place on its edge (short edges around a busy
    node, zoomed out) is left out. Plain labels never cover a node either;
    focused ones are placed first and avoid nodes when they can. A broken
    edge (drawn shortened, see ``MstEdgeScale``) shows a break mark, two
    slanted strokes across a gap in the line as in the Phylogeny, with its
    true distance beside it at any zoom; both are placed before every other
    label, clear of the nodes wherever the edge has room. One item paints
    them all.
    """

    _FONT_PX = 11
    _PAD_X = 5.0
    _PAD_Y = 2.0
    _FADE_PX = 18.0  # gap growth over which a midpoint label fades in
    _CELL_PX = 48.0  # collision grid cell
    _BREAK_HALF_PX = 5.0  # the break mark's half height (across the edge)
    _BREAK_GAP_PX = 3.0  # half the gap in the line

    def __init__(self, view_box: "pg.ViewBox") -> None:
        super().__init__()
        self._vb = view_box
        self._rect = QtCore.QRectF()
        self._positions = np.zeros((0, 2), dtype=float)
        self._left = np.zeros(0, dtype=np.int32)
        self._right = np.zeros(0, dtype=np.int32)
        self._texts: List[str] = []
        self._broken = np.zeros(0, dtype=bool)
        self._solid: Tuple[np.ndarray, np.ndarray] = (np.zeros(0), np.ones(0))
        self._focus: Set[int] = set()
        self._enlarged: Set[int] = set()
        self._faded = np.zeros(0, dtype=bool)
        self._node_radius_px = 10.5
        self._focus_radius_px = 13.0
        self._obstacle: Optional[QtWidgets.QGraphicsItem] = None
        self._font = QtGui.QFont(QtWidgets.QApplication.font())
        self._font.setPixelSize(self._FONT_PX)
        self._font_bold = QtGui.QFont(self._font)
        self._font_bold.setWeight(QtGui.QFont.Weight.DemiBold)
        self._sizes: Dict[Tuple[str, bool], Tuple[float, float]] = {}
        self._fill = QtGui.QColor(COLORS["surface"])
        self._border = QtGui.QColor(COLORS["border_strong"])
        self._ink = QtGui.QColor(COLORS["text_muted"])
        self._focus_border = QtGui.QColor(COLORS["hover"])
        self._focus_ink = QtGui.QColor(COLORS["hover_ink"])
        self._break_ink = QtGui.QColor(COLORS["text"])
        self.setZValue(25)
        self.setAcceptedMouseButtons(QtCore.Qt.MouseButton.NoButton)
        view_box.sigRangeChanged.connect(self._on_view_changed)
        view_box.sigResized.connect(self._on_view_changed)
        self._on_view_changed()

    def set_graph(
        self,
        positions: np.ndarray,
        left: np.ndarray,
        right: np.ndarray,
        texts: Sequence[str],
        broken: Sequence[int] = (),
        solid: Optional[Tuple[np.ndarray, np.ndarray]] = None,
    ) -> None:
        """``broken``: the indices of the edges drawn shortened; ``solid``:
        per edge, the part drawn solid as (t0, t1) from ``left`` (0) to
        ``right`` (1), where break marks go (default: all of it)."""
        self._positions = positions
        self._left = left
        self._right = right
        self._texts = list(texts)
        self._broken = np.zeros(int(left.size), dtype=bool)
        self._broken[[int(idx) for idx in broken if 0 <= int(idx) < left.size]] = True
        if solid is None or np.size(solid[0]) != left.size:
            solid = (np.zeros(left.size), np.ones(left.size))
        self._solid = solid
        self.update()

    def detach(self) -> None:
        """Stop following the view (the item is about to be removed)."""
        for signal in (self._vb.sigRangeChanged, self._vb.sigResized):
            try:
                signal.disconnect(self._on_view_changed)
            except (TypeError, RuntimeError):
                pass

    def set_obstacle(self, item: Optional[QtWidgets.QGraphicsItem]) -> None:
        """An item (the hover card) that labels keep clear of while it shows."""
        self._obstacle = item

    def set_positions(self, positions: np.ndarray) -> None:
        self._positions = positions
        self.update()

    def set_focus(
        self,
        focus: Set[int],
        faded: np.ndarray,
        *,
        enlarged: Set[int],
        node_radius_px: float,
        focus_radius_px: float,
    ) -> None:
        """``focus``: node indices whose edges always show their labels;
        ``faded``: per node, True when the node is dimmed (outside the
        window's samples, or not focused while something is);
        ``enlarged``: nodes drawn at the hover size (``focus_radius_px``)."""
        self._focus = set(focus)
        self._faded = faded
        self._enlarged = set(enlarged)
        self._node_radius_px = float(node_radius_px)
        self._focus_radius_px = float(focus_radius_px)
        self.update()

    def _on_view_changed(self, *_args: object) -> None:
        rect = self._vb.viewRect()
        if rect != self._rect:
            self.prepareGeometryChange()
            self._rect = QtCore.QRectF(rect)

    def boundingRect(self) -> QtCore.QRectF:  # type: ignore[override]
        return QtCore.QRectF(self._rect)

    def _label_size(self, text: str, bold: bool) -> Tuple[float, float]:
        key = (text, bold)
        size = self._sizes.get(key)
        if size is None:
            metrics = QtGui.QFontMetricsF(self._font_bold if bold else self._font)
            size = (
                metrics.horizontalAdvance(text) + 2.0 * self._PAD_X,
                metrics.height() + 2.0 * self._PAD_Y,
            )
            self._sizes[key] = size
        return size

    def _place_labels(
        self, tr: QtGui.QTransform
    ) -> Tuple[List[Tuple[float, bool, _Box, str]], List[Tuple[float, float, float, float, Optional[_Box]]]]:
        """The labels to draw through the painter transform ``tr``, in screen
        pixels: (alpha, emphasized, box, text), and the break marks: (alpha,
        x, y, angle of the edge in degrees, its label's box or None when
        there was no room for it)."""
        n_edges = int(self._left.size)
        if n_edges == 0 or self._positions.shape[0] == 0 or len(self._texts) != n_edges:
            return [], []
        px = self._positions[:, 0]
        py = self._positions[:, 1]
        sx = tr.m11() * px + tr.m21() * py + tr.dx()
        sy = tr.m12() * px + tr.m22() * py + tr.dy()
        x0, y0 = sx[self._left], sy[self._left]
        x1, y1 = sx[self._right], sy[self._right]
        dx = x1 - x0
        dy = y1 - y0
        length = np.hypot(dx, dy)
        visible = tr.mapRect(self._rect)
        vx0, vy0 = visible.left(), visible.top()
        vx1, vy1 = visible.right(), visible.bottom()
        r_node = self._node_radius_px
        r_focus = self._focus_radius_px
        focus = self._focus
        faded = self._faded if self._faded.size == self._positions.shape[0] else None
        median_px = float(np.median(length))
        zoom_gate = min(1.0, max(0.0, (median_px - (2.0 * r_node + 14.0)) / 20.0))
        broken = self._broken if self._broken.size == n_edges else np.zeros(n_edges, dtype=bool)
        if not focus and zoom_gate <= 0.0 and not broken.any():
            return [], []
        cell = self._CELL_PX
        enlarged = self._enlarged

        def radius(node: int) -> float:
            return r_focus if node in enlarged else r_node

        def cells_of(box: _Box) -> List[Tuple[int, int]]:
            return [
                (gx, gy)
                for gx in range(int(box[0] // cell), int(box[2] // cell) + 1)
                for gy in range(int(box[1] // cell), int(box[3] // cell) + 1)
            ]

        def add(grid: Dict[Tuple[int, int], List[_Box]], box: _Box) -> None:
            for key in cells_of(box):
                grid.setdefault(key, []).append(box)

        def hits(grid: Dict[Tuple[int, int], List[_Box]], box: _Box) -> bool:
            for key in cells_of(box):
                for other in grid.get(key, ()):
                    if box[0] < other[2] and box[2] > other[0] and box[1] < other[3] and box[3] > other[1]:
                        return True
            return False

        def label_box(cx: float, cy: float, w: float, h: float) -> _Box:
            return (cx - w / 2.0, cy - h / 2.0, cx + w / 2.0, cy + h / 2.0)

        def in_view(box: _Box) -> bool:
            return box[2] >= vx0 and box[0] <= vx1 and box[3] >= vy0 and box[1] <= vy1

        def inside_view(box: _Box) -> bool:
            return box[0] >= vx0 and box[2] <= vx1 and box[1] >= vy0 and box[3] <= vy1

        # Nodes in view block labels (a slightly smaller box than the disc,
        # so a label may touch a node's rim).
        node_grid: Dict[Tuple[int, int], List[_Box]] = {}
        label_grid: Dict[Tuple[int, int], List[_Box]] = {}
        on_screen = np.nonzero(
            (sx >= vx0 - r_focus) & (sx <= vx1 + r_focus) & (sy >= vy0 - r_focus) & (sy <= vy1 + r_focus)
        )[0]
        for node in on_screen.tolist():
            r = 0.8 * radius(node)
            add(node_grid, (sx[node] - r, sy[node] - r, sx[node] + r, sy[node] + r))

        obstacle = self._obstacle
        if obstacle is not None and obstacle.isVisible():
            card = tr.mapRect(self.mapRectFromScene(obstacle.sceneBoundingRect()))
            add(label_grid, (card.left(), card.top(), card.right(), card.bottom()))

        def along_extent(i: int, w: float, h: float) -> float:
            if length[i] <= 1e-9:
                return w
            return abs(w * dx[i] / length[i]) + abs(h * dy[i] / length[i])

        if focus:
            focus_arr = np.fromiter(focus, dtype=np.int64, count=len(focus))
            focus_edge = np.isin(self._left, focus_arr) | np.isin(self._right, focus_arr)
        else:
            focus_edge = np.zeros(n_edges, dtype=bool)

        # (alpha, emphasized, box, text): broken edges first, then focused.
        placed: List[Tuple[float, bool, _Box, str]] = []
        # Break marks: (alpha, x, y, angle in degrees, the label's box or None).
        marks: List[Tuple[float, float, float, float, Optional[_Box]]] = []
        half = self._BREAK_HALF_PX
        for i in np.nonzero(broken)[0].tolist():
            a, b = int(self._left[i]), int(self._right[i])
            if length[i] <= 1e-9:
                continue
            # The solid part of the edge inside the view.
            clipped = _clip_segments(x0[i : i + 1], y0[i : i + 1], dx[i : i + 1], dy[i : i + 1], (vx0, vy0, vx1, vy1))
            t0 = max(float(clipped[0][0]), float(self._solid[0][i]))
            t1 = min(float(clipped[1][0]), float(self._solid[1][i]))
            if t0 > t1:
                continue
            ux, uy = dx[i] / length[i], dy[i] / length[i]
            nx, ny = -uy, ux
            # The mark keeps off the end nodes; tried along the visible part,
            # from its middle outwards.
            lo = max(t0, (radius(a) + half + 3.0) / length[i])
            hi = min(t1, 1.0 - (radius(b) + half + 3.0) / length[i])
            if lo > hi:
                lo = hi = 0.5 * (t0 + t1)
            middle = 0.5 * (lo + hi)
            steps = [middle]
            for k in range(1, 7):
                steps.extend((middle - k * (hi - lo) / 12.0, middle + k * (hi - lo) / 12.0))
            emphasized = bool(focus_edge[i])
            text = self._texts[i]
            w, h = self._label_size(text, emphasized)
            off = half + 3.0 + 0.5 * (abs(w * nx) + abs(h * ny))
            # The label above the edge first (as in the Phylogeny), else below.
            sides = (1.0, -1.0) if ny < 0.0 or (ny == 0.0 and nx < 0.0) else (-1.0, 1.0)
            options: List[Tuple[float, float, _Box, _Box]] = []
            for t in steps:
                bx, by = x0[i] + t * dx[i], y0[i] + t * dy[i]
                mark = (bx - half, by - half, bx + half, by + half)
                for side in sides:
                    options.append((bx, by, mark, label_box(bx + side * nx * off, by + side * ny * off, w, h)))
            # Inside the view, clear of labels and nodes, and off the other
            # edges (so it reads as this edge's); else crossing an edge;
            # else only clear of labels; else (a crowded view, zoomed out)
            # the mark alone, in the middle.
            clear = [
                o for o in options
                if inside_view(o[3]) and not hits(label_grid, o[3]) and not hits(label_grid, o[2])
                and not hits(node_grid, o[3]) and not hits(node_grid, o[2])
            ]
            choice = (
                next((o for o in clear if not _segments_cross_box(x0, y0, dx, dy, o[3], skip=i)), None)
                or (clear[0] if clear else None)
                or next((o for o in options if inside_view(o[3]) and not hits(label_grid, o[3])), None)
            )
            bx, by, mark, box = choice if choice is not None else options[0]
            add(label_grid, mark)
            # Dimmed with a dimmed end, or, focused, with a dimmed focused end.
            ends = [node for node in (a, b) if node in focus] if emphasized else [a, b]
            alpha = _FADED_TEXT_OPACITY if faded is not None and any(faded[node] for node in ends) else 1.0
            labeled = choice is not None
            marks.append((alpha, float(bx), float(by), math.degrees(math.atan2(dy[i], dx[i])), box if labeled else None))
            if labeled:
                add(label_grid, (box[0] - 1.0, box[1] - 1.0, box[2] + 1.0, box[3] + 1.0))
                placed.append((alpha, emphasized, box, text))

        for i in np.nonzero(focus_edge & ~broken)[0].tolist():
            a, b = int(self._left[i]), int(self._right[i])
            text = self._texts[i]
            w, h = self._label_size(text, True)
            if a in focus:
                ox, oy, ux, uy, near = x0[i], y0[i], dx[i], dy[i], a
            else:
                ox, oy, ux, uy, near = x1[i], y1[i], -dx[i], -dy[i], b
            norm = math.hypot(ux, uy)
            ux, uy = (ux / norm, uy / norm) if norm > 1e-9 else (0.0, -1.0)
            candidates: List[_Box] = []
            if length[i] - radius(a) - radius(b) >= along_extent(i, w, h) + 4.0:
                candidates.append(label_box(x0[i] + 0.5 * dx[i], y0[i] + 0.5 * dy[i], w, h))
            # Otherwise (or when the midpoint is taken): along the edge,
            # just outside the focused node, then further out up to the
            # other node.
            reach = radius(near) + 3.0 + 0.5 * (abs(w * ux) + abs(h * uy))
            step = 0.6 * h + 2.0
            far_end = float(length[i])
            for k in range(10):
                along = reach + k * step
                if along > far_end:
                    break
                candidates.append(label_box(ox + ux * along, oy + uy * along, w, h))
            # Inside the view and free of nodes; else over a node (dimmed
            # or the far one); partly outside the view only when nothing
            # else fits. Never over another label.
            free = [c for c in candidates if not hits(label_grid, c)]
            choice = (
                next((c for c in free if inside_view(c) and not hits(node_grid, c)), None)
                or next((c for c in free if inside_view(c)), None)
                or next((c for c in free if in_view(c)), None)
            )
            if choice is not None:
                add(label_grid, choice)
                # A focused node outside the window's samples is greyed out;
                # so are its labels.
                alpha = _FADED_TEXT_OPACITY if faded is not None and faded[near] else 1.0
                placed.append((alpha, True, choice, text))

        # Every other edge: only zoomed in (by the on-screen length of the
        # median edge), then once its label fits between its nodes. Edges
        # with no room or whose midpoint is off screen are dropped up front.
        if zoom_gate > 0.0:
            mx = x0 + 0.5 * dx
            my = y0 + 0.5 * dy
            margin = 40.0
            plain = np.nonzero(
                ~focus_edge
                & ~broken
                & (length > 2.0 * r_node + 4.0)
                & (mx >= vx0 - margin) & (mx <= vx1 + margin)
                & (my >= vy0 - margin) & (my <= vy1 + margin)
            )[0]
            for i in plain.tolist():
                a, b = int(self._left[i]), int(self._right[i])
                text = self._texts[i]
                w, h = self._label_size(text, False)
                gap = float(length[i]) - 2.0 * r_node - along_extent(i, w, h) - 4.0
                if gap <= 0.0:
                    continue
                alpha = zoom_gate * min(1.0, gap / self._FADE_PX)
                if faded is not None and (faded[a] or faded[b]):
                    alpha *= _FADED_TEXT_OPACITY
                if alpha < 0.05:
                    continue
                box = label_box(float(mx[i]), float(my[i]), w, h)
                if in_view(box) and not hits(label_grid, box) and not hits(node_grid, box):
                    add(label_grid, (box[0] - 1.0, box[1] - 1.0, box[2] + 1.0, box[3] + 1.0))
                    placed.append((alpha, False, box, text))
        return placed, marks

    def paint(self, painter: QtGui.QPainter, _option: object, _widget: object = None) -> None:  # type: ignore[override]
        placed, marks = self._place_labels(painter.transform())
        if not placed and not marks:
            return
        painter.save()
        painter.resetTransform()
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        painter.setRenderHint(QtGui.QPainter.RenderHint.TextAntialiasing, True)
        gap = self._BREAK_GAP_PX
        half = self._BREAK_HALF_PX
        slant = 0.6 * half
        for alpha, px, py, angle, _box in marks:
            painter.save()
            painter.translate(px, py)
            painter.rotate(angle)
            # The gap: the line under it is covered with the plot's color.
            painter.fillRect(QtCore.QRectF(-gap, -2.5, 2.0 * gap, 5.0), self._fill)
            ink = QtGui.QColor(self._break_ink)
            ink.setAlphaF(alpha)
            pen = QtGui.QPen(ink, 1.3)
            pen.setCapStyle(QtCore.Qt.PenCapStyle.RoundCap)
            painter.setPen(pen)
            for x in (-gap, gap):
                painter.drawLine(QtCore.QPointF(x - slant, half), QtCore.QPointF(x + slant, -half))
            painter.restore()
        current_font: Optional[QtGui.QFont] = None
        # Plain labels under the emphasized ones.
        for alpha, emphasized, box, text in reversed(placed):
            fill = QtGui.QColor(self._fill)
            border = QtGui.QColor(self._focus_border if emphasized else self._border)
            ink = QtGui.QColor(self._focus_ink if emphasized else self._ink)
            fill.setAlphaF(0.92 * alpha)
            border.setAlphaF(alpha)
            ink.setAlphaF(alpha)
            rect = QtCore.QRectF(box[0], box[1], box[2] - box[0], box[3] - box[1])
            painter.setPen(QtGui.QPen(border, 1.4 if emphasized else 1.0))
            painter.setBrush(fill)
            painter.drawRoundedRect(rect, rect.height() / 2.0, rect.height() / 2.0)
            font = self._font_bold if emphasized else self._font
            if font is not current_font:
                painter.setFont(font)
                current_font = font
            painter.setPen(ink)
            painter.drawText(rect, int(QtCore.Qt.AlignmentFlag.AlignCenter), text)
        painter.restore()


def _paint_node_symbol(painter: QtGui.QPainter, symbol: str, rect: QtCore.QRectF) -> None:
    """A scatter symbol ("o", "s", "t", "d", "star", "+", "x", "p", "h",
    "t1", "t2", "t3") in ``rect`` with the painter's pen and brush."""
    cx, cy = rect.center().x(), rect.center().y()
    radius = 0.5 * min(rect.width(), rect.height())

    def polygon(points: Sequence[Tuple[float, float]]) -> None:
        painter.drawPolygon(QtGui.QPolygonF([QtCore.QPointF(x, y) for x, y in points]))

    def regular(sides: int, inner: Optional[float] = None) -> None:
        corners = sides * (2 if inner is not None else 1)
        points = []
        for idx in range(corners):
            angle = math.radians(-90.0 + 360.0 * idx / corners)
            r = radius if inner is None or idx % 2 == 0 else radius * inner
            points.append((cx + r * math.cos(angle), cy + r * math.sin(angle)))
        polygon(points)

    left, right, top, bottom = rect.left(), rect.right(), rect.top(), rect.bottom()
    if symbol == "s":
        painter.drawRect(rect)
    elif symbol == "t":
        # As pyqtgraph draws them: "t" points down, "t1" up.
        polygon([(cx, bottom), (left, top), (right, top)])
    elif symbol == "t1":
        polygon([(cx, top), (left, bottom), (right, bottom)])
    elif symbol == "t2":
        polygon([(left, top), (left, bottom), (right, cy)])
    elif symbol == "t3":
        polygon([(left, cy), (right, top), (right, bottom)])
    elif symbol == "d":
        polygon([(cx, top), (left, cy), (cx, bottom), (right, cy)])
    elif symbol == "p":
        regular(5)
    elif symbol == "h":
        regular(6)
    elif symbol == "star":
        regular(5, inner=0.45)
    elif symbol == "+":
        painter.drawLine(QtCore.QPointF(cx, top), QtCore.QPointF(cx, bottom))
        painter.drawLine(QtCore.QPointF(left, cy), QtCore.QPointF(right, cy))
    elif symbol == "x":
        painter.drawLine(QtCore.QPointF(left, top), QtCore.QPointF(right, bottom))
        painter.drawLine(QtCore.QPointF(left, bottom), QtCore.QPointF(right, top))
    else:
        painter.drawEllipse(rect)


@dataclass(frozen=True)
class GraphKeyEntry:
    """One key line: a node drawn with ``symbol`` in ``color`` (``color``
    None: the neutral shape key), its ``label`` and how many samples have it."""

    label: str
    count: int
    color: Optional[QtGui.QColor] = None
    symbol: str = "o"


def _nodes_in_rect(view_box: "pg.ViewBox", positions: np.ndarray, radius_px: float, rect: QtCore.QRectF) -> int:
    """How many nodes (discs of ``radius_px``) touch ``rect``, a rectangle in
    the ViewBox's own coordinates."""
    if positions.size == 0:
        return 0
    # From the view range rather than the child transform, which the
    # ViewBox only updates when it next paints.
    (x0, x1), (y0, y1) = view_box.viewRange()
    area = view_box.boundingRect()
    if x1 <= x0 or y1 <= y0:
        return 0
    sx = area.left() + (positions[:, 0] - x0) * (area.width() / (x1 - x0))
    sy = area.top() + (y1 - positions[:, 1]) * (area.height() / (y1 - y0))
    r = float(radius_px)
    inside = (
        (sx >= rect.left() - r) & (sx <= rect.right() + r)
        & (sy >= rect.top() - r) & (sy <= rect.bottom() + r)
    )
    return int(np.count_nonzero(inside))


def _view_corners(view_box: "pg.ViewBox", width: float, height: float, margin: float) -> List[QtCore.QRectF]:
    """A ``width`` by ``height`` rectangle in each corner of the ViewBox,
    ``margin`` in: bottom right, top right, bottom left, top left."""
    area = view_box.boundingRect()
    left, right = area.left() + margin, area.right() - margin - width
    top, bottom = area.top() + margin, area.bottom() - margin - height
    return [
        QtCore.QRectF(right, bottom, width, height),
        QtCore.QRectF(right, top, width, height),
        QtCore.QRectF(left, bottom, width, height),
        QtCore.QRectF(left, top, width, height),
    ]


class _GraphKey(QtWidgets.QGraphicsObject):
    """The graph's key, a card over the plot like the neighborhood view's.

    A child of the ViewBox, so it stays put while the graph pans and zooms.
    Each section is a caption line, then its entries side by side, wrapping
    like a figure legend. ``place`` puts it in the plot corner that covers the
    fewest nodes; the view calls it when it fits the graph, and the key
    follows resizes and new content, never a pan or zoom. When every corner
    would cover nodes it folds to a small "Key" pill (as does a click on the
    card), and a click on the pill opens it again, then for good.
    """

    # The key's captions or values changed (the view refits to keep room
    # for it); new counts alone do not move the graph.
    changed = QtCore.pyqtSignal()
    # The key picked its corner (the scale bar then picks another one).
    placed = QtCore.pyqtSignal()

    _FONT_PX = 11
    _MARGIN = 10.0
    _PAD_X = 11.0
    _PAD_Y = 8.0
    _SWATCH = 10.0
    _ENTRY_GAP = 10.0
    _MAX_WIDTH = 300.0

    def __init__(self, view_box: "pg.ViewBox") -> None:
        super().__init__(view_box)
        self._vb = view_box
        self._sections: List[Tuple[str, List[GraphKeyEntry]]] = []
        self._positions = np.zeros((0, 2), dtype=float)
        self._node_radius_px = 10.5
        self.collapsed = False
        self._forced_open = False
        self._auto_folded = False
        self._corner = 0
        self._rect = QtCore.QRectF()
        # Laid out items: (kind, x, y, entry or text), in card coordinates.
        self._rows: List[Tuple[str, float, float, object]] = []
        self._font = QtGui.QFont(QtWidgets.QApplication.font())
        self._font.setPixelSize(self._FONT_PX)
        self._bold = QtGui.QFont(self._font)
        self._bold.setBold(True)
        self._metrics = QtGui.QFontMetricsF(self._font)
        self._line = max(16.0, self._metrics.height() + 4.0)
        self.setZValue(40)
        self.setAcceptHoverEvents(True)
        self.setAcceptedMouseButtons(QtCore.Qt.MouseButton.LeftButton)
        self.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.setVisible(False)
        view_box.sigResized.connect(self.place)

    def set_sections(self, sections: Sequence[Tuple[str, Sequence[GraphKeyEntry]]]) -> None:
        before = [(title, [entry.label for entry in entries]) for title, entries in self._sections]
        self._sections = [(str(title), list(entries)) for title, entries in sections if entries]
        self.setVisible(bool(self._sections))
        self.place()
        if [(title, [entry.label for entry in entries]) for title, entries in self._sections] != before:
            self.changed.emit()

    def set_nodes(self, positions: np.ndarray, radius_px: float) -> None:
        """The nodes the key keeps clear of (used at the next ``place``)."""
        self._positions = positions
        self._node_radius_px = float(radius_px)

    def set_collapsed(self, collapsed: bool) -> None:
        self.collapsed = bool(collapsed)
        self._forced_open = not self.collapsed
        self.place()

    def reserved_size(self, view_width: float, view_height: float) -> Optional[QtCore.QSizeF]:
        """The room the open card takes with its margins (None: no card, or
        folded by the viewer)."""
        if not self._sections or self.collapsed:
            return None
        size = self._layout_card(view_width - 2.0 * self._MARGIN, view_height - 2.0 * self._MARGIN)
        return QtCore.QSizeF(size.width() + 2.0 * self._MARGIN, size.height() + 2.0 * self._MARGIN)

    def set_faded(self, hover_card: Optional[QtCore.QRectF]) -> None:
        """Fade while the hover card (a scene rectangle) overlaps the card."""
        faded = hover_card is not None and self.mapRectToScene(self._rect).intersects(hover_card)
        self.setOpacity(0.15 if faded else 1.0)

    def covers(self, scene_pos: QtCore.QPointF) -> bool:
        return self.isVisible() and self.mapRectToScene(self._rect).contains(scene_pos)

    def folded(self) -> bool:
        return self.collapsed or self._auto_folded

    def _entry_width(self, entry: GraphKeyEntry) -> float:
        return (
            self._SWATCH + 5.0 + self._metrics.horizontalAdvance(entry.label)
            + 4.0 + self._metrics.horizontalAdvance(f"({entry.count})")
        )

    def _layout_card(self, max_width: float, max_height: float) -> QtCore.QSizeF:
        """Lay out the open card within ``max_width`` by ``max_height``; the
        entries that do not fit are counted on a last "+N more" line."""
        line = self._line
        inner = max(80.0, min(self._MAX_WIDTH, max_width) - 2.0 * self._PAD_X)
        # Lines of (kind, x, value), entries flowing within the inner width.
        lines: List[List[Tuple[str, float, object]]] = []
        for title, entries in self._sections:
            lines.append([("caption", 0.0, title)])
            current: List[Tuple[str, float, object]] = []
            x = 0.0
            for entry in entries:
                width = self._entry_width(entry)
                if current and x + width > inner:
                    lines.append(current)
                    current, x = [], 0.0
                current.append(("entry", x, entry))
                x += width + self._ENTRY_GAP
            if current:
                lines.append(current)
        room = max(2, int((max_height - 2.0 * self._PAD_Y) // line) - 1)  # under the title
        hidden = 0
        if len(lines) > room:
            hidden = sum(1 for kept in lines[room - 1:] for kind, _x, _value in kept if kind == "entry")
            lines = lines[: room - 1]
            while lines and lines[-1][0][0] == "caption":
                lines.pop()  # a caption with none of its entries
            lines.append([("more", 0.0, f"+{hidden} more")])
        self._rows = []
        width = QtGui.QFontMetricsF(self._bold).horizontalAdvance("Key")
        y = self._PAD_Y + line
        for row in lines:
            for kind, x, value in row:
                self._rows.append((kind, self._PAD_X + x, y, value))
                if kind == "entry":
                    width = max(width, x + self._entry_width(value))
                else:
                    width = max(width, x + self._metrics.horizontalAdvance(str(value)))
            y += line
        return QtCore.QSizeF(width + 2.0 * self._PAD_X, y + self._PAD_Y)

    def _corners(self, width: float, height: float) -> List[QtCore.QRectF]:
        return _view_corners(self._vb, width, height, self._MARGIN)

    def _emptiest(self, corners: List[QtCore.QRectF]) -> Tuple[int, int]:
        counts = [_nodes_in_rect(self._vb, self._positions, self._node_radius_px, corner) for corner in corners]
        best = min(range(len(corners)), key=lambda idx: (counts[idx], idx))
        return best, counts[best]

    def place(self, *_args: object) -> None:
        """Pick the corner and whether the card fits there open."""
        area = self._vb.boundingRect()
        rect: Optional[QtCore.QRectF] = None
        self._auto_folded = False
        if self._sections and not self.collapsed:
            size = self._layout_card(area.width() - 2.0 * self._MARGIN, area.height() - 2.0 * self._MARGIN)
            corners = self._corners(size.width(), size.height())
            best, count = self._emptiest(corners)
            if count == 0 or self._forced_open:
                rect = corners[best]
                self._corner = best
            else:
                # Every corner would cover nodes: fold until asked.
                self._auto_folded = True
        if rect is None:
            width = QtGui.QFontMetricsF(self._bold).horizontalAdvance("Key") + 22.0
            corners = self._corners(width, 22.0)
            self._corner, _count = self._emptiest(corners)
            rect = corners[self._corner]
        if rect != self._rect:
            self.prepareGeometryChange()
            self._rect = rect
        self.setToolTip("Click to show the key" if self.folded() else "Click to fold the key")
        self.update()
        self.placed.emit()

    def boundingRect(self) -> QtCore.QRectF:  # type: ignore[override]
        return self._rect.adjusted(0.0, 0.0, 0.0, 2.0)

    def paint(self, painter: QtGui.QPainter, _option: object, _widget: object = None) -> None:  # type: ignore[override]
        rect = self._rect
        if rect.isEmpty():
            return
        painter.save()
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        painter.setRenderHint(QtGui.QPainter.RenderHint.TextAntialiasing, True)
        shadow = QtGui.QColor(COLORS["text_strong"])
        shadow.setAlpha(28)
        painter.setPen(QtCore.Qt.PenStyle.NoPen)
        painter.setBrush(shadow)
        painter.drawRoundedRect(rect.translated(0.0, 1.5), 6.0, 6.0)
        plate = QtGui.QColor(COLORS["surface"])
        plate.setAlpha(235)
        painter.setPen(QtGui.QPen(QtGui.QColor(COLORS["border_soft"]), 1.0))
        painter.setBrush(plate)
        painter.drawRoundedRect(rect, 6.0, 6.0)
        muted = QtGui.QColor(COLORS["text_muted"])
        painter.setFont(self._bold)
        painter.setPen(muted)
        if self.folded():
            painter.drawText(rect, int(QtCore.Qt.AlignmentFlag.AlignCenter), "Key")
            painter.restore()
            return
        line = self._line
        ox, oy = rect.left(), rect.top()
        painter.drawText(QtCore.QPointF(ox + self._PAD_X, oy + self._PAD_Y + line - 5.0), "Key")
        painter.setFont(self._font)
        text = QtGui.QColor(COLORS["text"])
        outline = QtGui.QColor(COLORS["text_strong"])
        neutral = QtGui.QColor(COLORS["border"])
        for kind, x, y, value in self._rows:
            baseline = oy + y + line - 5.0
            if kind in ("caption", "more"):
                painter.setPen(muted)
                painter.drawText(QtCore.QPointF(ox + x, baseline), str(value))
                continue
            entry: GraphKeyEntry = value  # type: ignore[assignment]
            middle = oy + y + 0.5 * line
            swatch = QtCore.QRectF(ox + x, middle - 0.5 * self._SWATCH, self._SWATCH, self._SWATCH)
            if entry.color is not None:
                painter.setPen(QtGui.QPen(QtGui.QColor(entry.color).darker(125), 0.8))
                painter.setBrush(QtGui.QColor(entry.color))
                painter.drawRoundedRect(swatch, 2.0, 2.0)
            else:
                painter.setPen(QtGui.QPen(outline, 1.1))
                painter.setBrush(neutral)
                _paint_node_symbol(painter, entry.symbol, swatch)
            label_x = ox + x + self._SWATCH + 5.0
            painter.setPen(text)
            painter.drawText(QtCore.QPointF(label_x, baseline), entry.label)
            painter.setPen(muted)
            painter.drawText(
                QtCore.QPointF(label_x + self._metrics.horizontalAdvance(entry.label) + 4.0, baseline),
                f"({entry.count})",
            )
        painter.restore()

    def mousePressEvent(self, event: QtWidgets.QGraphicsSceneMouseEvent) -> None:  # type: ignore[override]
        # Taking the press keeps the view from panning under the card.
        event.accept()

    def mouseReleaseEvent(self, event: QtWidgets.QGraphicsSceneMouseEvent) -> None:  # type: ignore[override]
        event.accept()
        if self._rect.contains(event.pos()):
            # A key folded to make room for nodes opens on click, like one
            # the viewer folded.
            self.set_collapsed(not self.folded())


class _ScaleBar(QtWidgets.QGraphicsObject):
    """The graph's scale bar, a small card in a corner of the plot: a bar a
    round number of alleles long ("10 alleles") on the drawing's linear
    scale, its length following the zoom.

    A child of the ViewBox, like the key. ``place`` puts it in the corner
    with the fewest nodes, never the key's (bottom left first, a little
    above the plot's auto range button); it runs whenever the key is placed
    (a fit, a resize, folding or opening the key). Pans and zooms only
    change the bar's length, in the same corner. Hidden while the graph has no scale (every
    profile identical) or one allele is longer than the widest bar.
    """

    _FONT_PX = 11
    _MARGIN = 10.0
    _PAD_X = 8.0
    _PAD_Y = 5.0
    _TARGET_PX = 70.0  # the bar's preferred length
    _MAX_PX = 150.0
    _TICK = 4.0
    _BUTTON_ROOM = 6.0  # extra room under the bar for the auto range button
    # Bottom left, bottom right, top left, top right (as _view_corners numbers them).
    _CORNER_ORDER = (2, 0, 3, 1)

    def __init__(self, view_box: "pg.ViewBox", key: Optional[_GraphKey]) -> None:
        super().__init__(view_box)
        self._vb = view_box
        self._key = key
        self._units_per_allele = 0.0
        self._positions = np.zeros((0, 2), dtype=float)
        self._node_radius_px = 10.5
        self._corner = 2
        self._alleles = 0
        self._bar_px = 0.0
        self._rect = QtCore.QRectF()
        self._font = QtGui.QFont(QtWidgets.QApplication.font())
        self._font.setPixelSize(self._FONT_PX)
        self._metrics = QtGui.QFontMetricsF(self._font)
        self.setZValue(40)
        self.setAcceptedMouseButtons(QtCore.Qt.MouseButton.NoButton)
        self.setVisible(False)
        if key is not None:
            key.placed.connect(self.place)
        else:
            view_box.sigResized.connect(self.place)
        view_box.sigRangeChanged.connect(self._on_range_changed)

    def set_scale(self, units_per_allele: float) -> None:
        """Layout units per allele (0: no scale to show)."""
        self._units_per_allele = max(0.0, float(units_per_allele))
        self.place()

    def set_nodes(self, positions: np.ndarray, radius_px: float) -> None:
        self._positions = positions
        self._node_radius_px = float(radius_px)

    def set_faded(self, hover_card: Optional[QtCore.QRectF]) -> None:
        """Fade while the hover card (a scene rectangle) overlaps the bar."""
        faded = hover_card is not None and self.mapRectToScene(self._rect).intersects(hover_card)
        self.setOpacity(0.15 if faded else 1.0)

    def _measure(self) -> Optional[Tuple[int, float]]:
        """The bar's alleles and length in pixels at the current zoom."""
        if self._units_per_allele <= 0.0:
            return None
        (x0, x1), _y_range = self._vb.viewRange()
        width = self._vb.boundingRect().width()
        if x1 <= x0 or width <= 0.0:
            return None
        px_per_allele = self._units_per_allele * width / (x1 - x0)
        if px_per_allele <= 0.0 or not math.isfinite(px_per_allele):
            return None
        # The round number (1, 2 or 5 times a power of ten) whose bar is
        # closest to the preferred length.
        best: Optional[Tuple[float, int, float]] = None
        for magnitude in (10 ** k for k in range(0, 7)):
            for step in (1, 2, 5):
                bar = step * magnitude * px_per_allele
                miss = abs(math.log(bar / self._TARGET_PX))
                if bar <= self._MAX_PX and (best is None or miss < best[0]):
                    best = (miss, step * magnitude, bar)
        return None if best is None else (best[1], best[2])

    def _text(self) -> str:
        return "1 allele" if self._alleles == 1 else f"{self._alleles} alleles"

    def _size(self) -> QtCore.QSizeF:
        width = max(self._bar_px, self._metrics.horizontalAdvance(self._text()))
        return QtCore.QSizeF(width + 2.0 * self._PAD_X, self._metrics.height() + self._TICK + 3.0 + 2.0 * self._PAD_Y)

    def place(self, *_args: object) -> None:
        """Measure the bar and pick its corner."""
        self._layout(choose_corner=True)

    def _on_range_changed(self, *_args: object) -> None:
        self._layout(choose_corner=False)

    def _layout(self, *, choose_corner: bool) -> None:
        measured = self._measure()
        if measured is None:
            self.setVisible(False)
            return
        self._alleles, self._bar_px = measured
        size = self._size()
        corners = _view_corners(self._vb, size.width(), size.height(), self._MARGIN)
        corners[2].translate(0.0, -self._BUTTON_ROOM)
        if choose_corner:
            key = self._key
            taken = key._corner if key is not None and key.isVisible() else None
            free = [idx for idx in self._CORNER_ORDER if idx != taken]
            counts = {
                idx: _nodes_in_rect(self._vb, self._positions, self._node_radius_px, corners[idx]) for idx in free
            }
            self._corner = min(free, key=lambda idx: (counts[idx], free.index(idx)))
        rect = corners[self._corner]
        if rect != self._rect:
            self.prepareGeometryChange()
            self._rect = rect
        self.setVisible(True)
        self.update()

    def boundingRect(self) -> QtCore.QRectF:  # type: ignore[override]
        return self._rect.adjusted(0.0, 0.0, 0.0, 2.0)

    def paint(self, painter: QtGui.QPainter, _option: object, _widget: object = None) -> None:  # type: ignore[override]
        rect = self._rect
        if rect.isEmpty():
            return
        painter.save()
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        painter.setRenderHint(QtGui.QPainter.RenderHint.TextAntialiasing, True)
        plate = QtGui.QColor(COLORS["surface"])
        plate.setAlpha(225)
        painter.setPen(QtGui.QPen(QtGui.QColor(COLORS["border_soft"]), 1.0))
        painter.setBrush(plate)
        painter.drawRoundedRect(rect, 5.0, 5.0)
        ink = QtGui.QColor(COLORS["text_muted"])
        # The bar with end ticks, centered over its caption.
        x0 = rect.center().x() - 0.5 * self._bar_px
        x1 = x0 + self._bar_px
        y = rect.top() + self._PAD_Y + self._TICK
        pen = QtGui.QPen(ink, 1.5)
        pen.setCapStyle(QtCore.Qt.PenCapStyle.FlatCap)
        painter.setPen(pen)
        painter.drawLine(QtCore.QPointF(x0, y), QtCore.QPointF(x1, y))
        painter.drawLine(QtCore.QPointF(x0, y - self._TICK), QtCore.QPointF(x0, y))
        painter.drawLine(QtCore.QPointF(x1, y - self._TICK), QtCore.QPointF(x1, y))
        painter.setFont(self._font)
        text_rect = QtCore.QRectF(rect.left(), y + 3.0, rect.width(), self._metrics.height())
        painter.drawText(text_rect, int(QtCore.Qt.AlignmentFlag.AlignCenter), self._text())
        painter.restore()


class CgmlstForceGraphView(QtWidgets.QWidget):
    """The cgMLST minimum spanning tree of the run's genomes.

    The layout (``mst_drawing``) is computed once per allele table, in
    node diameters; a node is drawn ``_NODE_PX`` wide, or at most one layout
    unit (less its outline and the edge line) when the view is zoomed out
    further, so drawn nodes never overlap and never touch an edge at any zoom.
    Edge lengths are linear in the allele distance, very long edges drawn
    shortened with a break mark and their true distance (``MstEdgeScale``);
    a scale bar shows the scale. Where the layout lengthened an edge to make
    room for a node's children, only its length to scale is solid (from the
    parent's end) and the rest is dotted, so solid lines read on the bar. Genomes outside the window's samples, and
    their edges, are drawn in the shared out of set grey (``out_of_set_color``).
    While the view fits the graph, the tree's long axis follows the plot's.
    """

    sampleHovered = QtCore.pyqtSignal(object)
    _NODE_PX = 21.0
    _HOVER_NODE_PX = 26.0
    _EDGE_COLOR = QtGui.QColor(148, 163, 184, 170)

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
        _log_view_event(f"created cgMLST MST view for {self._output_dir}")
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
        self._edge_left_idx = np.zeros(0, dtype=np.int32)
        self._edge_right_idx = np.zeros(0, dtype=np.int32)
        self._refresh_token = 0
        self._active_refresh_token: Optional[int] = None
        self._worker_refs: Dict[int, RefreshWorker] = {}
        self._last_request_key: Optional[Tuple[object, ...]] = None
        self._last_payload: Optional[_CgmlstForcePayload] = None
        self._needs_refresh = True
        self._positions = np.zeros((0, 2), dtype=float)
        # The tree turned a quarter turn (the long axis vertical).
        self._upright = False
        # Node size as a fraction of _NODE_PX: below 1 while the view shows
        # a layout unit smaller than _NODE_PX pixels.
        self._node_scale = 1.0
        self._hover_node_px = self._HOVER_NODE_PX
        self._line_factor = 1.0
        self._positions_revision = 0
        self._point_style_key_by_index: Dict[int, Tuple[object, ...]] = {}
        self._brush_cache: Dict[Tuple[int, int, int, int], QtGui.QBrush] = {}
        self._pen_cache: Dict[Tuple[int, int, int, int, float, bool], QtGui.QPen] = {}
        self._last_edge_highlight_key: Optional[Tuple[Tuple[str, ...], int]] = None
        self._scatter_item: Optional["pg.ScatterPlotItem"] = None
        self._edge_item: Optional["pg.PlotDataItem"] = None
        self._edge_faded_item: Optional["pg.PlotDataItem"] = None
        self._edge_highlight_item: Optional["pg.PlotDataItem"] = None
        # The dotted parts of the lengthened edges, per layer: "faded",
        # "plain" and "highlight", each over the matching solid item.
        self._edge_dotted_items: Dict[str, "pg.PlotDataItem"] = {}
        self._scale_start = np.zeros(0, dtype=np.int64)
        self._scale_share = np.zeros(0, dtype=float)
        self._hover_text: Optional["pg.TextItem"] = None
        self._distance_labels: Optional[_EdgeDistanceLabels] = None
        self._label_focus_key: Optional[Tuple[object, ...]] = None
        self._key: Optional[_GraphKey] = None
        self._scale_bar: Optional[_ScaleBar] = None

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
        view_box = self._plot.getViewBox()
        view_box.sigRangeChanged.connect(self._on_view_scale_changed)
        view_box.sigResized.connect(self._on_view_scale_changed)
        self._key = _GraphKey(view_box)
        self._key.changed.connect(self._on_key_changed)
        # After the key: on a resize the bar picks its corner once the key has.
        self._scale_bar = _ScaleBar(view_box, self._key)
        # Fit the graph to the view when it is built and when the view is
        # resized (docked panels change size often), until the user pans or
        # zooms.
        self._auto_fit = True
        self._fit_timer = QtCore.QTimer(self)
        self._fit_timer.setSingleShot(True)
        self._fit_timer.setInterval(30)
        self._fit_timer.timeout.connect(self._fit_view_to_nodes)
        self._plot.getViewBox().sigRangeChangedManually.connect(self._on_range_changed_manually)
        self._hover_proxy = pg.SignalProxy(
            self._plot.scene().sigMouseMoved,
            rateLimit=60,
            slot=self._on_mouse_moved,
        )
        # The scene reports moves only inside the plot: leaving it clears
        # the hover (else the last hovered node stays highlighted).
        self._plot.viewport().installEventFilter(self)

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
        self._apply_edge_fade()
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

    def set_key_sections(self, sections: Sequence[Tuple[str, Sequence[GraphKeyEntry]]]) -> None:
        """The key card's sections: (caption, entries); none hides the key."""
        if self._key is not None:
            self._key.set_sections(sections)

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
        if self._distance_labels is not None:
            self._distance_labels.detach()
        self._plot.clear()
        self._scatter_item = None
        self._edge_item = None
        self._edge_faded_item = None
        self._edge_highlight_item = None
        self._edge_dotted_items = {}
        self._hover_text = None
        self._distance_labels = None
        self._label_focus_key = None
        self._hover_index = None
        self._sample_labels = list(payload.sample_labels)
        self._coords_by_index = {}
        self._point_by_index = {}
        self._base_color_by_index = {}
        self._sample_by_index = {}
        self._called_loci_by_index = {}
        self._sample_called_loci = tuple(int(value) for value in payload.sample_called_loci)
        self._edges = list(payload.edges)
        self._upright = False
        self._positions = np.array(payload.coords, dtype=float, copy=True)
        self._edge_left_idx = np.asarray([int(left) for left, _right, _dist in self._edges], dtype=np.int32)
        self._edge_right_idx = np.asarray([int(right) for _left, right, _dist in self._edges], dtype=np.int32)
        self._scale_start = np.asarray(payload.scale_start, dtype=np.int64)
        self._scale_share = np.asarray(payload.scale_share, dtype=float)
        self._positions_revision = 0
        self._point_style_key_by_index = {}
        self._last_edge_highlight_key = None

        self._edge_faded_item = pg.PlotDataItem([], [], connect="finite")
        self._edge_item = pg.PlotDataItem([], [], connect="finite")
        self._edge_highlight_item = pg.PlotDataItem([], [], connect="finite")
        for layer, item in (("faded", self._edge_faded_item), ("plain", self._edge_item), ("highlight", self._edge_highlight_item)):
            self._plot.addItem(item)
            self._edge_dotted_items[layer] = pg.PlotDataItem([], [], connect="finite")
            self._plot.addItem(self._edge_dotted_items[layer])
        self._apply_edge_pens()
        self._scatter_item = pg.ScatterPlotItem()
        self._plot.addItem(self._scatter_item)
        self._distance_labels = _EdgeDistanceLabels(self._plot.getViewBox())
        self._plot.addItem(self._distance_labels, ignoreBounds=True)
        n_loci = int(payload.n_loci)
        self._distance_labels.set_graph(
            self._positions,
            self._edge_left_idx,
            self._edge_right_idx,
            [str(int(distance)) if int(distance) <= n_loci else "n/a" for _left, _right, distance in self._edges],
            broken=payload.scale.broken_edges([distance for _left, _right, distance in self._edges]),
            solid=self._solid_ranges(),
        )
        if self._scale_bar is not None:
            self._scale_bar.set_scale(payload.scale.units_per_allele if any(d > 0 for _l, _r, d in self._edges) else 0.0)
        # A card behind the hover text keeps it readable over edges and labels.
        card = QtGui.QColor(COLORS["surface"])
        card.setAlpha(235)
        self._hover_text = pg.TextItem(
            color=COLORS["text_strong"],
            anchor=(0.0, 1.4),
            border=pg.mkPen(COLORS["border_strong"]),
            fill=pg.mkBrush(card),
        )
        self._hover_text.setZValue(30)
        self._hover_text.setVisible(False)
        if self._key is not None:
            self._key.set_faded(None)
        if self._scale_bar is not None:
            self._scale_bar.set_faded(None)
        self._plot.addItem(self._hover_text)
        self._distance_labels.set_obstacle(self._hover_text)

        self._auto_fit = True
        self._fit_view_to_nodes(render=False)
        self._node_scale, self._hover_node_px, self._line_factor = self._current_node_sizes()
        self._render_scene()
        self._update_info_text(payload)

    def _update_info_text(self, payload: _CgmlstForcePayload) -> None:
        lineage_count = len({value for value in self._lineage_by_sample.values() if value})
        _log_view_event(
            f"cgMLST MST | samples: {len(payload.sample_labels)} | loci: {payload.n_loci} | "
            f"MST edges: {len(payload.edges)} | mean pairwise allele distance: {payload.mean_pairwise_distance:.2f} | "
            f"fastGEAR clades: {lineage_count} | duplicate cgMLST profiles: {payload.duplicate_profile_groups}"
        )
        self._info.clear()
        self._info.setVisible(False)

    def _render_scene(self) -> None:
        if pg is None or self._edge_item is None or self._edge_highlight_item is None or self._scatter_item is None:
            return
        self._coords_by_index = {}
        self._base_color_by_index = {}
        self._sample_by_index = {}
        self._positions_revision += 1
        self._last_edge_highlight_key = None
        self._apply_edge_fade()
        self._apply_edge_pens()
        if self._distance_labels is not None:
            self._distance_labels.set_positions(self._positions)
        self._place_cards()

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
        self._apply_label_focus(emphasis, active_filter)
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
                self._place_hover_card(coords)

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
        muted = COLORS["text_muted"]
        lines = [html.escape(sample)]
        if lineage:
            lines.append(f"<span style='color:{muted}'>Clade: {html.escape(lineage)}</span>")
        if cluster:
            lines.append(f"<span style='color:{muted}'>Cluster: {html.escape(cluster)}</span>")
        for label, value in self._sample_hover_metadata.get(sample, {}).items():
            lines.append(
                f"<span style='color:{muted}'>{html.escape(label)}: {html.escape(value)}</span>"
            )
        called_loci = self._called_loci_by_index.get(idx)
        total_loci = int(self._last_payload.n_loci) if self._last_payload is not None else 0
        if called_loci is not None and total_loci > 0:
            pct = (100.0 * float(called_loci)) / float(total_loci)
            lines.append(
                f"cgMLST called loci: {pct:.1f}% ({called_loci}/{total_loci})"
            )
        self._hover_text.setHtml("<br>".join(lines))
        self._place_hover_card(coords)
        self._hover_text.setVisible(True)

    def _place_hover_card(self, coords: Tuple[float, float]) -> None:
        """Put the hover card above and right of its node, flipped left or
        below when that side would leave the plot."""
        card = self._hover_text
        if card is None:
            return
        x, y = float(coords[0]), float(coords[1])
        vb = self._plot.getViewBox()
        node = vb.mapViewToScene(QtCore.QPointF(x, y))
        area = vb.sceneBoundingRect()
        size = card.boundingRect()
        w, h = float(size.width()), float(size.height())
        right_room = area.right() - node.x()
        left_room = node.x() - area.left()
        ax = 0.0 if right_room >= w or right_room >= left_room else 1.0
        above_room = node.y() - area.top()
        below_room = area.bottom() - node.y()
        ay = 1.4 if above_room >= 1.4 * h or above_room >= below_room else -0.4
        if tuple(card.anchor) != (ax, ay):
            card.setAnchor((ax, ay))
        card.setPos(x, y)
        # The key and the scale bar are drawn above the graph: they fade
        # under the hover card.
        hover_card = QtCore.QRectF(node.x() - ax * w, node.y() - ay * h, w, h)
        if self._key is not None:
            self._key.set_faded(hover_card)
        if self._scale_bar is not None:
            self._scale_bar.set_faded(hover_card)

    def _on_range_changed_manually(self, *_args: object) -> None:
        self._auto_fit = False

    def _on_key_changed(self) -> None:
        # New key content: the fit keeps room for the open key.
        if pg is not None and self._auto_fit:
            self._fit_timer.start()

    @staticmethod
    def _edge_width(scale: float) -> float:
        # Thinner lines with small nodes, so they stay lighter than the nodes.
        return 1.3 * min(1.0, max(0.6, 0.5 + scale))

    @staticmethod
    def _outline_width(scale: float, *, hover: bool = False, focus: bool = False) -> float:
        return (2.4 if hover else (2.0 if focus else 1.0)) * (0.5 + 0.5 * scale)

    def _current_node_sizes(self) -> Tuple[float, float, float]:
        """At the current zoom: the node size as a fraction of _NODE_PX, the
        hovered node's size in pixels, and the factor on the line widths.

        The layout keeps node centers 1.4 units apart and 0.7 units from the
        edges they are not on. A node is at most one unit wide, and its
        outline (an emphasized one) and the edge line get thinner while the
        unit is a few pixels, so a drawn node never touches another node or
        an edge."""
        if pg is None or self._positions.size == 0:
            return 1.0, self._HOVER_NODE_PX, 1.0
        try:
            unit_px = 1.0 / float(self._plot.getViewBox().viewPixelSize()[0])
        except (ZeroDivisionError, TypeError, ValueError):
            return 1.0, self._HOVER_NODE_PX, 1.0
        if not math.isfinite(unit_px) or unit_px <= 0.0:
            return 1.0, self._HOVER_NODE_PX, 1.0
        node = min(self._NODE_PX, unit_px)
        scale = node / self._NODE_PX
        edge = self._edge_width(scale)
        # The emphasized edges are drawn twice as wide.
        lines = min(1.0, (1.4 * unit_px - node) / (self._outline_width(scale, focus=True) + 2.0 * edge))
        hover = min(
            self._HOVER_NODE_PX * scale,
            1.4 * unit_px - lines * (self._outline_width(scale, hover=True) + 2.0 * edge),
        )
        return scale, max(node, hover), lines

    def _on_view_scale_changed(self, *_args: object) -> None:
        scale, hover, lines = self._current_node_sizes()
        if abs(scale - self._node_scale) <= 0.02 * self._node_scale:
            return
        self._node_scale = scale
        self._hover_node_px = hover
        self._line_factor = lines
        self._apply_node_scale()

    def _apply_node_scale(self) -> None:
        """Resize the nodes, their outlines and the edge lines after a zoom,
        without rebuilding the scene."""
        if self._scatter_item is None or self._edge_item is None or not self._point_by_index:
            return
        self._apply_edge_pens()
        radius = 0.5 * self._NODE_PX * self._node_scale
        if self._key is not None:
            self._key.set_nodes(self._positions, radius)
        if self._scale_bar is not None:
            self._scale_bar.set_nodes(self._positions, radius)
        active_filter = self._sample_filter
        emphasis = self._emphasis_samples()
        pens: List[QtGui.QPen] = []
        sizes = np.empty(len(self._sample_labels), dtype=float)
        for idx, sample in enumerate(self._sample_labels):
            base = self._base_color_by_index.get(idx, QtGui.QColor("#475569"))
            _brush, pen, size, style_key = self._point_style_components(idx, sample, base, emphasis, active_filter)
            pens.append(pen)
            sizes[idx] = size
            self._point_style_key_by_index[idx] = style_key
        self._scatter_item.setPen(pens, update=False)
        self._scatter_item.setSize(sizes)
        self._apply_label_focus(emphasis, active_filter)

    def _apply_edge_pens(self) -> None:
        """The edge pens at the current node size: plain, faded (the edges of
        genomes outside the window's samples) and highlighted, each solid
        and dotted (the lengthened part of an edge)."""
        width = self._line_factor * self._edge_width(self._node_scale)
        layers = (
            ("faded", self._edge_faded_item, out_of_set_color(self._EDGE_COLOR), width),
            ("plain", self._edge_item, QtGui.QColor(self._EDGE_COLOR), width),
            ("highlight", self._edge_highlight_item, QtGui.QColor(COLORS["hover"]), 2.0 * width),
        )
        for layer, item, color, layer_width in layers:
            if item is not None:
                item.setPen(self._cached_pen(color, layer_width))
            dotted = self._edge_dotted_items.get(layer)
            if dotted is not None:
                dotted.setPen(self._cached_pen(color, layer_width, dotted=True))

    def _solid_ranges(self) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        """Per edge, its part drawn solid as (t0, t1) along left to right."""
        if self._scale_share.size != self._edge_left_idx.size:
            return None
        from_left = self._scale_start == self._edge_left_idx
        share = self._scale_share
        return np.where(from_left, 0.0, 1.0 - share), np.where(from_left, share, 1.0)

    def _edge_segments(self, mask: np.ndarray) -> Tuple[Tuple[np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray]]:
        """The masked edges as line data (x, y with a gap after each
        segment): their parts drawn to scale, and the lengthened rest."""
        left, right = self._edge_left_idx[mask], self._edge_right_idx[mask]
        if self._scale_share.size == self._edge_left_idx.size:
            start, share = self._scale_start[mask], self._scale_share[mask]
        else:
            start, share = left, np.ones(left.size, dtype=float)
        end = np.where(start == left, right, left)
        a, b = self._positions[start], self._positions[end]
        knee = a + share[:, np.newaxis] * (b - a)
        longer = share < 1.0

        def lines(p: np.ndarray, q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
            gap = np.full(p.shape[0], np.nan)
            return np.column_stack([p[:, 0], q[:, 0], gap]).ravel(), np.column_stack([p[:, 1], q[:, 1], gap]).ravel()

        return lines(a, knee), lines(knee[longer], b[longer])

    def _place_cards(self) -> None:
        """Place the key; the scale bar follows (never in the key's corner)."""
        radius = 0.5 * self._NODE_PX * self._node_scale
        if self._scale_bar is not None:
            self._scale_bar.set_nodes(self._positions, radius)
        if self._key is not None:
            self._key.set_nodes(self._positions, radius)
            self._key.place()
        elif self._scale_bar is not None:
            self._scale_bar.place()

    def _apply_edge_fade(self) -> None:
        """Split the edges: those of a genome outside the window's samples
        are drawn in the faint out of set grey, under the others. An edge the
        layout lengthened is solid for its length to scale, from its parent
        node, and dotted beyond."""
        if self._edge_item is None:
            return
        n = int(self._positions.shape[0])
        left, right = self._edge_left_idx, self._edge_right_idx
        valid = (left < n) & (right < n)
        faded = np.zeros(left.size, dtype=bool)
        if self._sample_filter is not None and n:
            outside = ~np.fromiter((s in self._sample_filter for s in self._sample_labels[:n]), dtype=bool, count=n)
            faded[valid] = outside[left[valid]] | outside[right[valid]]
        for layer, item, mask in (("plain", self._edge_item, valid & ~faded), ("faded", self._edge_faded_item, valid & faded)):
            solid, dotted = self._edge_segments(mask)
            if item is not None:
                item.setData(*solid)
            if layer in self._edge_dotted_items:
                self._edge_dotted_items[layer].setData(*dotted)

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:  # type: ignore[override]
        super().resizeEvent(event)
        if pg is not None and getattr(self, "_auto_fit", False):
            self._fit_timer.start()

    def _fit_ranges(self, coords: np.ndarray, area: QtCore.QRectF) -> Tuple[float, List[float], List[float]]:
        """The view ranges that show every node with a small margin, and the
        pixels per layout unit they give (aspect locked)."""
        min_x, max_x = float(np.min(coords[:, 0])), float(np.max(coords[:, 0]))
        min_y, max_y = float(np.min(coords[:, 1])), float(np.max(coords[:, 1]))
        # A node (a layout unit wide) plus a little air around the tree.
        margin = 1.0 + 0.03 * max(max_x - min_x, max_y - min_y)
        x_range = [min_x - margin, max_x + margin]
        y_range = [min_y - margin, max_y + margin]
        if area.width() <= 0 or area.height() <= 0:
            return 0.0, x_range, y_range
        span_x, span_y = x_range[1] - x_range[0], y_range[1] - y_range[0]
        full = min(area.width() / span_x, area.height() / span_y)
        # Room for the open key beside or below the tree, unless that costs
        # the graph more than a third of its size.
        key = self._key.reserved_size(area.width(), area.height()) if self._key is not None else None
        if key is None:
            return full, x_range, y_range
        beside = min(max(0.0, area.width() - key.width()) / span_x, area.height() / span_y)
        below = min(area.width() / span_x, max(0.0, area.height() - key.height()) / span_y)
        if max(beside, below) < (2.0 / 3.0) * full:
            return full, x_range, y_range
        if beside >= below:
            x_range[1] += key.width() / beside
            return beside, x_range, y_range
        y_range[0] -= key.height() / below
        return below, x_range, y_range

    def _fit_view_to_nodes(self, *, render: bool = True) -> None:
        """Show every node with a small margin (aspect stays locked), the
        tree's long axis along the plot's when that shows it larger."""
        if pg is None or self._positions.size == 0:
            return
        vb = self._plot.getViewBox()
        if vb is None:
            return
        area = vb.boundingRect()
        coords = np.asarray(self._positions, dtype=float)
        # A quarter turn either way: (x, y) -> (-y, x) and back.
        turned = np.column_stack([coords[:, 1], -coords[:, 0]]) if self._upright else np.column_stack([-coords[:, 1], coords[:, 0]])
        scale, x_range, y_range = self._fit_ranges(coords, area)
        turned_scale, turned_x, turned_y = self._fit_ranges(turned, area)
        if turned_scale > 1.1 * scale:
            self._upright = not self._upright
            self._positions = turned
            x_range, y_range = turned_x, turned_y
        vb.disableAutoRange()
        vb.setRange(xRange=tuple(x_range), yRange=tuple(y_range), padding=0.0)
        if self._positions is turned and render:
            self._render_scene()
        else:
            self._place_cards()

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
        last_line = str(error_text or "").strip().splitlines()[-1] if str(error_text or "").strip() else "Unknown error"
        self._info.setText(f"Failed to build the cgMLST MST: {last_line}")
        self._info.setVisible(True)

    def eventFilter(self, obj: QtCore.QObject, event: QtCore.QEvent) -> bool:  # type: ignore[override]
        if event.type() == QtCore.QEvent.Type.Leave and pg is not None and obj is self._plot.viewport():
            self._clear_hover()
        return super().eventFilter(obj, event)

    def hideEvent(self, event: QtGui.QHideEvent) -> None:  # type: ignore[override]
        super().hideEvent(event)
        if pg is not None:
            self._clear_hover()

    def _clear_hover(self) -> None:
        # Drop a mouse move still waiting in the rate limiter, so it cannot
        # bring the hover back after the pointer left.
        self._hover_proxy.timer.stop()
        self._hover_proxy.args = None
        self._set_hover_index(None)

    def _set_hover_index(self, idx: Optional[int]) -> None:
        if idx == self._hover_index:
            return
        self._hover_index = idx
        if idx is None:
            if self._hover_text is not None:
                self._hover_text.setVisible(False)
            if self._key is not None:
                self._key.set_faded(None)
            if self._scale_bar is not None:
                self._scale_bar.set_faded(None)
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
        if pos is None or not self._plot.sceneBoundingRect().contains(pos) or (self._key is not None and self._key.covers(pos)):
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

    def _cached_pen(self, color: QtGui.QColor, width: float, *, dotted: bool = False) -> QtGui.QPen:
        rounded_width = round(float(width), 3)
        key = (
            int(color.red()),
            int(color.green()),
            int(color.blue()),
            int(color.alpha()),
            rounded_width,
            dotted,
        )
        pen = self._pen_cache.get(key)
        if pen is None:
            pen = pg.mkPen(QtGui.QColor(key[0], key[1], key[2], key[3]), width=rounded_width)
            if dotted:
                # Short dashes, about 2 px on and 3 px off at any line width.
                unit = max(0.5, rounded_width)
                pen.setDashPattern([2.0 / unit, 3.0 / unit])
                pen.setCapStyle(QtCore.Qt.PenCapStyle.FlatCap)
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
        hover_sample = self._sample_by_index.get(self._hover_index if self._hover_index is not None else -1, "")
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
        is_hover = self._hover_index == idx
        is_focus = in_filter and sample in emphasis
        size = self._hover_node_px if is_hover else self._NODE_PX * self._node_scale
        if in_filter:
            brush_color = QtGui.QColor(base)
            brush_color.setAlpha(215 if active else 72)
            pen_color = QtGui.QColor("#0f172a" if (is_hover or is_focus) else "#334155")
        else:
            # Outside the window's samples: the faint grey of the Phylogeny
            # and the heatmap.
            brush_color = out_of_set_color()
            pen_color = out_of_set_color("#0f172a" if is_hover else None)
        pen_width = self._line_factor * self._outline_width(self._node_scale, hover=is_hover, focus=is_focus)
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
        self._apply_label_focus(emphasis, active_filter)

    def _apply_label_focus(self, emphasis: Set[str], active_filter: Optional[set[str]]) -> None:
        """The hovered node (here, or a single sample hovered in another
        panel) shows the distances of all its edges."""
        labels = self._distance_labels
        if labels is None:
            return
        focus_samples: Set[str] = set()
        if self._highlight_sample:
            focus_samples.add(self._highlight_sample)
        hover_sample = self._sample_by_index.get(self._hover_index if self._hover_index is not None else -1, "")
        if hover_sample:
            focus_samples.add(hover_sample)
        key = (
            frozenset(focus_samples),
            frozenset(emphasis),
            active_filter,
            self._hover_index,
            self._node_scale,
        )
        if key == self._label_focus_key:
            return
        self._label_focus_key = key
        focus = {idx for idx, sample in self._sample_by_index.items() if sample in focus_samples}
        n = len(self._sample_labels)
        faded = np.zeros(n, dtype=bool)
        if active_filter is not None:
            faded |= ~np.fromiter((s in active_filter for s in self._sample_labels), dtype=bool, count=n)
        if emphasis:
            faded |= ~np.fromiter((s in emphasis for s in self._sample_labels), dtype=bool, count=n)
        labels.set_focus(
            focus,
            faded,
            enlarged={self._hover_index} if self._hover_index is not None else set(),
            node_radius_px=0.5 * self._NODE_PX * self._node_scale,
            focus_radius_px=0.5 * self._hover_node_px,
        )

    def _apply_edge_highlight(self) -> None:
        if pg is None or self._edge_highlight_item is None:
            return
        emphasis = set(self._highlight_samples)
        if self._highlight_sample:
            emphasis.add(self._highlight_sample)
        hover_sample = self._sample_by_index.get(self._hover_index if self._hover_index is not None else -1, "")
        if hover_sample:
            emphasis.add(hover_sample)
        edge_key = (tuple(sorted(emphasis)), int(self._positions_revision))
        if self._last_edge_highlight_key == edge_key:
            return
        self._last_edge_highlight_key = edge_key
        dotted_item = self._edge_dotted_items.get("highlight")
        if not emphasis:
            self._edge_highlight_item.setData([], [])
            if dotted_item is not None:
                dotted_item.setData([], [])
            return
        n = int(self._positions.shape[0])
        emphasized = np.zeros(n, dtype=bool)
        emphasized[[idx for idx, sample in self._sample_by_index.items() if sample in emphasis and idx < n]] = True
        left, right = self._edge_left_idx, self._edge_right_idx
        valid = (left < n) & (right < n)
        mask = np.zeros(left.size, dtype=bool)
        mask[valid] = emphasized[left[valid]] | emphasized[right[valid]]
        solid, dotted = self._edge_segments(mask)
        self._edge_highlight_item.setData(*solid)
        if dotted_item is not None:
            dotted_item.setData(*dotted)

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
        if len(request_key) >= 2 and request_key[0] == "ok":
            print(f"[file] using cgMLST allele calls: {request_key[1]}")
        self._info.setText("Building the cgMLST MST...")
        self._info.setVisible(True)
        worker = RefreshWorker(
            token,
            lambda token=token, request_key=request_key, reference_sample=self._reference_sample:
                _build_cgmlst_force_payload(
                    token,
                    request_key,
                    self._output_dir,
                    reference_sample,
                ),
        )
        worker.signals.finished.connect(self._on_refresh_finished)
        worker.signals.failed.connect(self._on_refresh_failed)
        self._worker_refs[token] = worker
        QtCore.QThreadPool.globalInstance().start(worker)
