"""A 3D point cloud drawn with plain pyqtgraph (no OpenGL).

``Projection3D`` rotates the points with a 3 x 3 matrix and projects them
with a mild perspective onto the plane; the view draws them back to front in
an ordinary ``ScatterPlotItem``. ``OrbitViewBox`` is a ViewBox whose mouse
turns the cloud while its orbit mode is on (drag rotates around the center,
the wheel zooms, a double-click resets) and pans and zooms as usual
otherwise, so one plot shows the 2D and the 3D view.

Two depth cues besides the drawing order and the depth levels: a reference
box around the cloud (``ReferenceBox``: a wireframe, a faint grid on the
three walls away from the viewer and the axis titles on the near outline,
as matplotlib's 3D panes) and sphere-like shading of the markers
(``sphere_brush``, rendered once per style into the scatter's atlas).
"""
from __future__ import annotations

from typing import List, NamedTuple, Optional, Sequence, Tuple

import numpy as np
from PyQt6 import QtCore, QtGui

import pyqtgraph as pg

# Radians per pixel of drag.
ROTATE_PER_PX = 0.01
ZOOM_STEP = 1.15
ZOOM_RANGE = (0.3, 8.0)
# The camera sits this far from the center (radius 1): f = D / (D - z).
CAMERA_DISTANCE = 3.0
# The fixed view range of the 3D view, both axes.
VIEW_EXTENT = 1.2
# Depth cue levels, far to near: size factor and share of the background.
DEPTH_SIZES = (0.55, 0.66, 0.78, 0.89, 1.0, 1.11, 1.23, 1.35)
DEPTH_FADES = (0.6, 0.5, 0.4, 0.3, 0.2, 0.12, 0.05, 0.0)
# The reference box: padding around the cloud (share of its largest span)
# and about this many grid steps along that span.
BOX_PAD = 0.05
BOX_STEPS = 6
# Marker shading (shares toward white and black): the highlight at the top
# left and the rim.
SHADE_HIGHLIGHT = 0.8
SHADE_RIM = 0.5


def _rot_x(angle: float) -> np.ndarray:
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def _rot_y(angle: float) -> np.ndarray:
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


# A slightly turned start, seen a little from above (the reference box's
# floor is behind the points), so the three axes are all seen.
DEFAULT_ROTATION = _rot_x(0.45) @ _rot_y(0.6)


class Projection3D:
    """Points centered and scaled into the unit ball (one scale for every
    axis, so the MDS axes keep their proportions), a rotation and a zoom."""

    def __init__(self) -> None:
        self.points = np.zeros((0, 3))
        # The cloud's center and scale: points = (coords - center) / scale.
        self.center = np.zeros(3)
        self.scale = 1.0
        self.rotation = DEFAULT_ROTATION.copy()
        self.zoom = 1.0
        self._updates = 0

    def set_points(self, points: np.ndarray) -> None:
        p = np.asarray(points, dtype=float).reshape(-1, 3)
        self.center = np.zeros(3)
        self.scale = 1.0
        if p.shape[0]:
            self.center = p.mean(axis=0)
            p = p - self.center
            radius = float(np.sqrt((p * p).sum(axis=1)).max())
            if radius > 0:
                p = p / radius
                self.scale = radius
        self.points = p

    def reset(self) -> None:
        self.rotation = DEFAULT_ROTATION.copy()
        self.zoom = 1.0

    def rotate(self, dx: float, dy: float) -> None:
        """Turn by a drag of ``dx``, ``dy`` pixels (screen axes, through the center)."""
        self.rotation = _rot_y(dx * ROTATE_PER_PX) @ _rot_x(dy * ROTATE_PER_PX) @ self.rotation
        self._updates += 1
        if self._updates % 100 == 0:
            # Keep it a rotation despite rounding.
            q, r = np.linalg.qr(self.rotation)
            self.rotation = q * np.sign(np.diag(r))[None, :]

    def zoom_by(self, steps: float) -> None:
        self.zoom = float(np.clip(self.zoom * ZOOM_STEP ** float(steps), *ZOOM_RANGE))

    def project(self, points: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
        """Screen xy (view units) and depth z (larger is nearer) of ``points``
        (default: the cloud)."""
        p = self.points if points is None else np.asarray(points, dtype=float).reshape(-1, 3)
        r = p @ self.rotation.T
        # The view's y axis points up, as the rotated y.
        f = CAMERA_DISTANCE / (CAMERA_DISTANCE - r[:, 2])
        return r[:, :2] * (f * self.zoom)[:, None], r[:, 2]

    def state(self) -> dict:
        return {"rotation": [float(v) for v in self.rotation.ravel()], "zoom": float(self.zoom)}

    def set_state(self, state: object) -> None:
        if not isinstance(state, dict):
            return
        values = state.get("rotation")
        if isinstance(values, (list, tuple)) and len(values) == 9:
            try:
                matrix = np.asarray([float(v) for v in values]).reshape(3, 3)
            except (TypeError, ValueError):
                matrix = None
            if matrix is not None and np.allclose(matrix @ matrix.T, np.eye(3), atol=1e-3):
                self.rotation = matrix
        try:
            self.zoom = float(np.clip(float(state.get("zoom", self.zoom)), *ZOOM_RANGE))
        except (TypeError, ValueError):
            pass


def depth_levels(z: np.ndarray) -> np.ndarray:
    """Depth (-1 far to 1 near) in len(DEPTH_SIZES) levels, 0 the farthest."""
    n = len(DEPTH_SIZES)
    return np.clip(((np.asarray(z) + 1.0) * 0.5 * n).astype(int), 0, n - 1)


def nice_step(span: float) -> float:
    """The 1, 2 or 5 times a power of ten at or just above ``span``."""
    if not np.isfinite(span) or span <= 0:
        return 1.0
    power = 10.0 ** np.floor(np.log10(span))
    for factor in (1.0, 2.0, 5.0, 10.0):
        if factor * power >= span * (1 - 1e-9):
            return float(factor * power)
    return float(10.0 * power)


def _hull_exit(points: np.ndarray, origin: np.ndarray, direction: np.ndarray) -> float:
    """How far from ``origin`` along the unit ``direction`` the ray leaves
    the convex hull of ``points`` (0 when it starts outside or on it)."""
    pts = sorted(map(tuple, np.asarray(points, dtype=float)))

    def chain(seq: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
        out: List[Tuple[float, float]] = []
        for q in seq:
            while len(out) >= 2 and (
                (out[-1][0] - out[-2][0]) * (q[1] - out[-2][1]) - (out[-1][1] - out[-2][1]) * (q[0] - out[-2][0])
            ) <= 0:
                out.pop()
            out.append(q)
        return out[:-1]

    hull = np.asarray(chain(pts) + chain(pts[::-1]))  # counterclockwise
    if hull.shape[0] < 3:
        return 0.0
    edges = np.roll(hull, -1, axis=0) - hull
    normals = np.column_stack([edges[:, 1], -edges[:, 0]])  # outward
    slack = ((hull - origin) * normals).sum(axis=1)  # >= 0 inside
    if (slack < 0).any():
        return 0.0
    speed = normals @ direction
    leaving = speed > 1e-12
    return float((slack[leaving] / speed[leaving]).min()) if leaving.any() else 0.0


class BoxFrame(NamedTuple):
    """One frame of the reference box, projected (view units)."""

    panes: List[np.ndarray]  # the three far walls, 4 x 2 each
    grid: np.ndarray  # pairs of points (``connect="pairs"``): the grid and the three edges at the near corner
    edges: np.ndarray  # pairs: the nine other edges and the tick marks
    labels: np.ndarray  # 3 x 2, the axis titles' positions
    anchors: List[Tuple[float, float]]  # their TextItem anchors
    outward: np.ndarray  # 3 x 2, the unit screen direction each title sits in
    corners: np.ndarray  # 8 x 2


class ReferenceBox:
    """A box around the normalized cloud with grid lines at round values of
    the original coordinates (one step for every axis)."""

    def __init__(self) -> None:
        self.lo = np.zeros(3)
        self.hi = np.zeros(3)
        self.ticks: List[np.ndarray] = [np.zeros(0)] * 3

    def fit(self, projection: Projection3D) -> None:
        p = projection.points
        if p.shape[0] == 0:
            self.lo, self.hi = np.zeros(3), np.zeros(3)
            self.ticks = [np.zeros(0)] * 3
            return
        lo, hi = p.min(axis=0), p.max(axis=0)
        pad = BOX_PAD * float(np.max(hi - lo))
        self.lo, self.hi = lo - pad, hi + pad
        scale, center = projection.scale, projection.center
        step = nice_step(float(np.max(hi - lo)) * scale / BOX_STEPS)
        self.ticks = []
        for axis in range(3):
            first = np.ceil((self.lo[axis] * scale + center[axis]) / step)
            last = np.floor((self.hi[axis] * scale + center[axis]) / step)
            values = np.arange(first, last + 1) * step
            self.ticks.append((values - center[axis]) / scale)

    def radius(self) -> float:
        """The farthest corner from the center of rotation."""
        corners = np.stack(np.meshgrid(*zip(self.lo, self.hi), indexing="ij"), axis=-1).reshape(-1, 3)
        return float(np.sqrt((corners * corners).sum(axis=1)).max())

    def frame(self, projection: Projection3D, tick: float, gap: float) -> BoxFrame:
        """The box as turned now; tick marks ``tick`` and titles ``gap`` view
        units out from the box's outline."""
        lo, hi = self.lo, self.hi
        mid = 0.5 * (lo + hi)
        # The camera in the cloud's coordinates: the far wall of each axis is
        # the one on the other side of the center.
        camera = projection.rotation.T @ np.array([0.0, 0.0, CAMERA_DISTANCE])
        far_is_hi = camera < mid

        def side(axis: int, far: bool) -> float:
            return float(hi[axis] if far == bool(far_is_hi[axis]) else lo[axis])

        segments: List[Tuple[np.ndarray, np.ndarray]] = []  # 3D pairs
        grid: List[Tuple[np.ndarray, np.ndarray]] = []
        panes3d: List[np.ndarray] = []
        for axis in range(3):
            j, k = [a for a in range(3) if a != axis]
            for vj in (lo[j], hi[j]):
                for vk in (lo[k], hi[k]):
                    a = np.empty(3)
                    a[j], a[k], a[axis] = vj, vk, lo[axis]
                    b = a.copy()
                    b[axis] = hi[axis]
                    # The edges at the corner nearest the viewer cross the
                    # cloud in front of points: drawn as faint as the grid.
                    near = vj == side(j, False) and vk == side(k, False)
                    (grid if near else segments).append((a, b))
            wall = side(axis, True)
            pane = np.empty((4, 3))
            pane[:, axis] = wall
            pane[:, j] = (lo[j], hi[j], hi[j], lo[j])
            pane[:, k] = (lo[k], lo[k], hi[k], hi[k])
            panes3d.append(pane)
            for along, across in ((j, k), (k, j)):
                for t in self.ticks[along]:
                    a = np.empty(3)
                    a[axis], a[along], a[across] = wall, t, lo[across]
                    b = a.copy()
                    b[across] = hi[across]
                    grid.append((a, b))

        def project_pairs(pairs: List[Tuple[np.ndarray, np.ndarray]]) -> np.ndarray:
            if not pairs:
                return np.zeros((0, 2))
            xy, _z = projection.project(np.asarray(pairs).reshape(-1, 3))
            return xy

        center_xy = projection.project(mid)[0][0]
        corners = np.stack(np.meshgrid(*zip(lo, hi), indexing="ij"), axis=-1).reshape(-1, 3)
        corners_xy = projection.project(corners)[0]
        labels = np.zeros((3, 2))
        outward = np.zeros((3, 2))
        anchors: List[Tuple[float, float]] = []
        ticks: List[np.ndarray] = []
        for axis in range(3):
            j, k = [a for a in range(3) if a != axis]
            # The edges on the box's outline: on one far wall and one near
            # wall. The title goes on the nearer of the two.
            best = None
            for fj, fk in ((True, False), (False, True)):
                a = np.empty(3)
                a[j], a[k], a[axis] = side(j, fj), side(k, fk), lo[axis]
                b = a.copy()
                b[axis] = hi[axis]
                xy, z = projection.project(np.vstack([a, b, 0.5 * (a + b)]))
                if best is None or z[2] > best[2]:
                    best = (xy, a, z[2])
            assert best is not None
            xy, a, _depth = best
            along = xy[1] - xy[0]
            length = float(np.hypot(*along))
            direction = along / length if length > 1e-9 else np.zeros(2)
            out = xy[2] - center_xy
            out = out - direction * float(out @ direction)
            norm = float(np.hypot(*out))
            if norm < 1e-9:
                out = xy[2] - center_xy
                norm = float(np.hypot(*out))
            out = out / norm if norm > 1e-9 else np.array([0.0, -1.0])
            outward[axis] = out
            # Past the box's outline: an edge seen end on (looking down its
            # axis) is a short line inside the outline, not on it.
            reach = _hull_exit(corners_xy, xy[2], out)
            labels[axis] = xy[2] + out * (reach + gap)
            anchors.append((float(0.5 - 0.5 * out[0]), float(0.5 + 0.5 * out[1])))
            if self.ticks[axis].size:
                pts = np.repeat(a[None, :], self.ticks[axis].size, axis=0)
                pts[:, axis] = self.ticks[axis]
                txy, _z = projection.project(pts)
                ticks.append(np.stack([txy, txy + out * tick], axis=1).reshape(-1, 2))
        return BoxFrame(
            panes=[projection.project(p)[0] for p in panes3d],
            grid=project_pairs(grid),
            edges=np.vstack([project_pairs(segments)] + ticks),
            labels=labels,
            anchors=anchors,
            outward=outward,
            corners=corners_xy,
        )


def _mix(color: QtGui.QColor, target: QtGui.QColor, share: float) -> QtGui.QColor:
    """``color`` moved ``share`` of the way to ``target`` (alpha kept)."""
    return QtGui.QColor(
        int(round(color.red() + (target.red() - color.red()) * share)),
        int(round(color.green() + (target.green() - color.green()) * share)),
        int(round(color.blue() + (target.blue() - color.blue()) * share)),
        color.alpha(),
    )


# Marker shapes whose fill leaves most of their bounding box empty (thin
# arms, triangles): their shading is drawn nearer the middle.
def sphere_brush(color: QtGui.QColor) -> QtGui.QBrush:
    """A shaded sphere of ``color``: a bright highlight toward the top left,
    the color itself around the middle and a darker rim, so a round marker
    reads as a ball lit from the top left (a radial gradient over the
    marker's bounding box, so one brush fits every size)."""
    white, black = QtGui.QColor(255, 255, 255), QtGui.QColor(0, 0, 0)
    gradient = QtGui.QRadialGradient(QtCore.QPointF(0.5, 0.5), 0.55, QtCore.QPointF(0.34, 0.3))
    gradient.setCoordinateMode(QtGui.QGradient.CoordinateMode.ObjectBoundingMode)
    gradient.setColorAt(0.0, _mix(color, white, SHADE_HIGHLIGHT))
    gradient.setColorAt(0.18, _mix(color, white, 0.35))
    gradient.setColorAt(0.5, color)
    gradient.setColorAt(0.82, _mix(color, black, 0.2))
    gradient.setColorAt(1.0, _mix(color, black, SHADE_RIM))
    return QtGui.QBrush(gradient)


class OrbitViewBox(pg.ViewBox):
    """A ViewBox that turns a 3D cloud while ``orbit`` is on."""

    orbitDragged = QtCore.pyqtSignal(float, float)
    orbitZoomed = QtCore.pyqtSignal(float)

    def __init__(self, *args: object, **kwargs: object) -> None:
        super().__init__(*args, **kwargs)
        self.orbit = False

    def set_orbit(self, on: bool) -> None:
        self.orbit = bool(on)
        self.setMouseEnabled(not self.orbit, not self.orbit)

    def mouseDragEvent(self, ev: object, axis: Optional[int] = None) -> None:  # type: ignore[override]
        if not self.orbit:
            super().mouseDragEvent(ev, axis)
            return
        if ev.button() != QtCore.Qt.MouseButton.LeftButton:  # type: ignore[attr-defined]
            ev.ignore()  # type: ignore[attr-defined]
            return
        ev.accept()  # type: ignore[attr-defined]
        delta = ev.screenPos() - ev.lastScreenPos()  # type: ignore[attr-defined]
        self.orbitDragged.emit(float(delta.x()), float(delta.y()))

    def wheelEvent(self, ev: object, axis: Optional[int] = None) -> None:  # type: ignore[override]
        if not self.orbit:
            super().wheelEvent(ev, axis)
            return
        ev.accept()  # type: ignore[attr-defined]
        delta = ev.delta() if hasattr(ev, "delta") else ev.angleDelta().y()  # type: ignore[attr-defined]
        self.orbitZoomed.emit(float(delta) / 120.0)


def fade_toward(rgba: Sequence[int], background: Sequence[int], share: float) -> Tuple[int, int, int, int]:
    """``rgba`` moved ``share`` of the way to ``background`` (alpha kept)."""
    r, g, b = (int(round(c + (bg - c) * share)) for c, bg in zip(rgba[:3], background[:3]))
    return r, g, b, int(rgba[3]) if len(rgba) > 3 else 255


class DepthScatterItem(pg.ScatterPlotItem):
    """A scatter whose points are redrawn in a new order and at new places
    every frame without restyling them.

    ``set_level_styles`` styles every point once per depth level (one
    ``setData``, which renders the symbols into the item's atlas);
    ``show_frame`` then picks each point's styled record at its level, in
    drawing order, and only moves it. Restyling thousands of points per frame
    is what makes ``setData`` slow; picking records is a numpy copy.
    """

    def __init__(self, *args: object, **kwargs: object) -> None:
        super().__init__(*args, **kwargs)
        self._template: Optional[np.ndarray] = None
        self._styles: Optional[Tuple[object, ...]] = None
        self._n = 0
        self._dpr = 1.0

    def has_level_styles(self) -> bool:
        return self._template is not None

    def drop_level_styles(self) -> None:
        self._template = None
        self._styles = None

    def set_level_styles(
        self, brushes: np.ndarray, sizes: np.ndarray, symbols: Sequence[object], pens: Sequence[object]
    ) -> None:
        """``brushes`` and ``sizes`` are levels x points; ``symbols`` and
        ``pens`` one per point."""
        levels, n = brushes.shape
        self._styles = (brushes, sizes, symbols, pens)
        self.setData(
            x=np.zeros(levels * n),
            y=np.zeros(levels * n),
            brush=list(brushes.ravel()),
            size=np.asarray(sizes, dtype=float).ravel(),
            symbol=list(symbols) * levels,
            pen=list(pens) * levels,
        )
        self._template = self.data.copy()
        self._n = n
        self._dpr = self.fragmentAtlas.devicePixelRatio()

    def show_frame(self, xy: np.ndarray, points: np.ndarray, levels: np.ndarray) -> None:
        """Draw ``points`` (in this order) at ``xy``, each at its depth level."""
        if self._template is None:
            return
        if self.fragmentAtlas.devicePixelRatio() != self._dpr and self._styles is not None:
            # A paint on a screen of another pixel ratio re-rendered the atlas.
            self.set_level_styles(*self._styles)  # type: ignore[arg-type]
        data = self._template[np.asarray(levels, dtype=np.intp) * self._n + np.asarray(points, dtype=np.intp)]
        data["x"] = xy[:, 0]
        data["y"] = xy[:, 1]
        self.data = data
        self.prepareGeometryChange()
        self.informViewBoundsChanged()
        self.bounds = [None, None]
        self.invalidate()

    def clear(self) -> None:  # type: ignore[override]
        self.drop_level_styles()
        super().clear()
