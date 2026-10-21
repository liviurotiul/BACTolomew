"""The band under the Pangenome clustering plot: the gene presence/absence
of the analyzed isolates as a small image, genes by prevalence (core on the
left, rare genes on the right), with a two handle range over it that picks
the part of the pangenome the clustering uses.

A handle at prevalence p sits after the genes the range leaves out on its
side, so the range covers exactly the genes in it. The genes outside are
veiled; the Panaroo categories (core, soft core, shell, cloud) are marked
along the top. Exact values go in the two spin boxes beside the band.
"""
from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

import numpy as np
from PyQt6 import QtCore, QtGui, QtWidgets

from .model.pan_clustering import genes_in_range
from .settings_column import ElidedHint
from .theme import COLORS

BAND_HEIGHT = 76
STEP = 0.005  # prevalence step of the handles (0.5%)
_STEPS = 200  # 1 / STEP, so a snapped value is an exact division
MIN_GAP = STEP
_HANDLE_GRAB_PX = 8
_GRIP_W = 9.0
# Shift while dragging: one step of 0.5% per this many pixels.
_FINE_PX = 4.0
# Panaroo's categories: core >= 99%, soft core 95% to 99%, shell 15% to 95%, cloud below.
_CATEGORIES = (("core", 0.99), ("soft core", 0.95), ("shell", 0.15), ("cloud", 0.0))


def snap(value: float) -> float:
    return float(np.clip(round(float(value) * _STEPS) / _STEPS, 0.0, 1.0))


def percent_text(value: float) -> str:
    text = f"{100.0 * float(value):.1f}"
    return (text[:-2] if text.endswith(".0") else text) + "%"


def _rgb(name: str) -> np.ndarray:
    color = QtGui.QColor(COLORS[name])
    return np.array([color.red(), color.green(), color.blue()], dtype=np.float32)


def image_from_shares(shares: np.ndarray) -> QtGui.QImage:
    """The band image: each cell's share of present genes, from the surface
    color (absent) to the text color (present)."""
    shares = np.clip(np.asarray(shares, dtype=np.float32), 0.0, 1.0)
    low, high = _rgb("surface"), _rgb("text")
    rgb = (low[None, None, :] + (high - low)[None, None, :] * shares[:, :, None]).astype(np.uint8)
    rgb = np.ascontiguousarray(rgb)
    h, w = rgb.shape[:2]
    image = QtGui.QImage(rgb.data, w, h, 3 * w, QtGui.QImage.Format.Format_RGB888)
    return image.copy()


class _BandCanvas(QtWidgets.QWidget):
    """The image, the veils and the two handles."""

    # (low, high) while dragging, then once more on release.
    rangeEdited = QtCore.pyqtSignal(float, float)
    rangeCommitted = QtCore.pyqtSignal(float, float)
    handleDoubleClicked = QtCore.pyqtSignal(str)

    _TOP = 13  # the category strip
    _PAD = 6

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setMouseTracking(True)
        self.setFocusPolicy(QtCore.Qt.FocusPolicy.StrongFocus)
        self.setMinimumWidth(120)
        self.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Expanding)
        self.prevalence = np.zeros(0)
        self.low, self.high = 0.30, 0.50
        self._image: Optional[QtGui.QImage] = None
        self._scaled: Optional[QtGui.QPixmap] = None
        self._drag: Optional[str] = None  # "low", "high" or "range"
        self._drag_origin: Tuple[float, float, float] = (0.0, 0.0, 0.0)
        self._fine: Optional[Tuple[float, float, float]] = None  # Shift drag: (x, low, high)
        self._active = "low"
        self._font = QtGui.QFont(self.font())
        self._font.setPixelSize(10)
        self.setAccessibleName("Prevalence range")

    # ---------------------------------------------------------- data
    def set_prevalence(self, prevalence: np.ndarray) -> None:
        self.prevalence = np.asarray(prevalence, dtype=float)
        self.update()

    def set_image(self, image: Optional[QtGui.QImage]) -> None:
        self._image = image
        self._scaled = None
        self.update()

    def set_range(self, low: float, high: float) -> None:
        self.low, self.high = float(low), float(high)
        self.update()

    # ---------------------------------------------------------- geometry
    def _area(self) -> QtCore.QRectF:
        return QtCore.QRectF(self.rect()).adjusted(self._PAD, self._TOP, -self._PAD, -2)

    def _genes(self) -> int:
        return int(self.prevalence.size)

    def _x_of_index(self, index: float) -> float:
        area = self._area()
        n = max(1, self._genes())
        return area.left() + area.width() * float(index) / n

    def _x_of(self, which: str) -> float:
        p = self.prevalence
        if which == "high":
            index = int(np.count_nonzero(p >= self.high)) if self.high < 1.0 else 0
        else:
            index = int(np.count_nonzero(p > self.low))
        return self._x_of_index(index)

    def _value_at(self, x: float) -> float:
        """The prevalence of the gene under ``x`` (0 past the last gene)."""
        n = self._genes()
        if n == 0:
            return 0.0
        area = self._area()
        index = int(np.floor((x - area.left()) / max(1.0, area.width()) * n))
        if index < 0:
            return 1.0
        if index >= n:
            return 0.0
        return float(self.prevalence[index])

    # ---------------------------------------------------------- painting
    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        self._scaled = None
        super().resizeEvent(event)

    def paintEvent(self, event: QtGui.QPaintEvent) -> None:
        del event
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        area = self._area()
        painter.fillRect(area, QtGui.QColor(COLORS["surface"]))
        if self._image is not None and not self._image.isNull():
            size = area.size().toSize()
            if self._scaled is None or self._scaled.size() != size:
                self._scaled = QtGui.QPixmap.fromImage(
                    self._image.scaled(
                        size,
                        QtCore.Qt.AspectRatioMode.IgnoreAspectRatio,
                        QtCore.Qt.TransformationMode.SmoothTransformation,
                    )
                )
            painter.drawPixmap(area.topLeft(), self._scaled)
        self._paint_categories(painter, area)
        if self._genes():
            x_high, x_low = self._x_of("high"), self._x_of("low")
            veil = QtGui.QColor(COLORS["window"])
            veil.setAlpha(150)
            painter.fillRect(QtCore.QRectF(area.left(), area.top(), x_high - area.left(), area.height()), veil)
            painter.fillRect(QtCore.QRectF(x_low, area.top(), area.right() - x_low, area.height()), veil)
            painter.setPen(QtGui.QPen(QtGui.QColor(COLORS["primary"]), 1.5))
            painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
            painter.drawRect(QtCore.QRectF(x_high, area.top() + 0.75, max(1.0, x_low - x_high), area.height() - 1.5))
            for which, x in (("high", x_high), ("low", x_low)):
                self._paint_handle(painter, x, area, which)
        painter.setPen(QtGui.QPen(QtGui.QColor(COLORS["border"]), 1.0))
        painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
        painter.drawRect(area)
        painter.end()

    def _paint_categories(self, painter: QtGui.QPainter, area: QtCore.QRectF) -> None:
        if not self._genes():
            return
        painter.setFont(self._font)
        metrics = QtGui.QFontMetricsF(self._font)
        muted = QtGui.QColor(COLORS["text_muted"])
        edges = [self._x_of_index(np.count_nonzero(self.prevalence >= bound)) for _name, bound in _CATEGORIES]
        start = area.left()
        for (name, _bound), end in zip(_CATEGORIES, edges):
            width = end - start
            if width > 2:
                if start > area.left() + 0.5:
                    painter.setPen(QtGui.QPen(QtGui.QColor(COLORS["border_strong"]), 1.0, QtCore.Qt.PenStyle.DotLine))
                    painter.drawLine(QtCore.QPointF(start, 1.0), QtCore.QPointF(start, area.top()))
                if metrics.horizontalAdvance(name) + 6 < width:
                    painter.setPen(muted)
                    painter.drawText(
                        QtCore.QRectF(start, 0.0, width, area.top()),
                        int(QtCore.Qt.AlignmentFlag.AlignCenter),
                        name,
                    )
            start = end

    def _grip_center(self, x: float, which: str) -> float:
        """Grips sit outside the range (the high one left of its line, the
        low one right of it), so a narrow range keeps both in reach; at the
        band's ends they stay inside it."""
        center = x - _GRIP_W / 2.0 if which == "high" else x + _GRIP_W / 2.0
        area = self._area()
        return float(np.clip(center, area.left() + _GRIP_W / 2.0, area.right() - _GRIP_W / 2.0))

    def _paint_handle(self, painter: QtGui.QPainter, x: float, area: QtCore.QRectF, which: str) -> None:
        active = (self.hasFocus() and self._active == which) or self._drag == which
        color = QtGui.QColor(COLORS["primary_hover"] if active else COLORS["primary"])
        painter.setPen(QtCore.Qt.PenStyle.NoPen)
        painter.setBrush(color)
        painter.drawRect(QtCore.QRectF(x - 1.0, area.top(), 2.0, area.height()))
        center = self._grip_center(x, which)
        grip = QtCore.QRectF(center - _GRIP_W / 2.0, area.center().y() - 10.0, _GRIP_W, 20.0)
        painter.setPen(QtGui.QPen(QtGui.QColor(COLORS["surface"]), 1.2))
        painter.drawRoundedRect(grip, 3.0, 3.0)
        painter.setPen(QtGui.QPen(QtGui.QColor(COLORS["surface"]), 1.0))
        for dx in (-1.5, 1.5):
            painter.drawLine(
                QtCore.QPointF(center + dx, grip.top() + 6.0), QtCore.QPointF(center + dx, grip.bottom() - 6.0)
            )

    # ---------------------------------------------------------- mouse
    def _handle_at(self, x: float) -> Optional[str]:
        d_high = abs(x - self._grip_center(self._x_of("high"), "high"))
        d_low = abs(x - self._grip_center(self._x_of("low"), "low"))
        if min(d_high, d_low) > _HANDLE_GRAB_PX:
            return None
        return "high" if d_high < d_low else "low"

    def _set_values(self, low: float, high: float) -> None:
        low, high = snap(low), snap(high)
        if high - low < MIN_GAP - 1e-9:
            return
        if (low, high) != (self.low, self.high):
            self.low, self.high = low, high
            self.update()
            self.rangeEdited.emit(low, high)

    def _move_handle(self, which: str, value: float) -> None:
        if which == "low":
            self._set_values(min(snap(value), self.high - MIN_GAP), self.high)
        else:
            self._set_values(self.low, max(snap(value), self.low + MIN_GAP))

    def mousePressEvent(self, event: QtGui.QMouseEvent) -> None:
        if event.button() != QtCore.Qt.MouseButton.LeftButton or not self._genes():
            super().mousePressEvent(event)
            return
        x = event.position().x()
        handle = self._handle_at(x)
        if handle is None:
            x_high, x_low = self._x_of("high"), self._x_of("low")
            if x_high < x < x_low:
                self._drag = "range"
                self._drag_origin = (x, x_high, x_low)
            else:
                handle = "high" if x < x_high else "low"
                self._move_handle(handle, self._value_at(x))
        if handle is not None:
            self._drag = handle
            self._active = handle
        self._fine = None
        self.update()
        event.accept()

    def _fine_drag(self, x: float) -> None:
        """Shift held: the dragged handle or range moves 0.5% per few pixels
        from where Shift was first seen, whatever the genes under it."""
        if self._fine is None:
            self._fine = (x, self.low, self.high)
        start, low, high = self._fine
        delta = -round((x - start) / _FINE_PX) * STEP  # right is toward rare genes
        if self._drag == "low":
            self._move_handle("low", low + delta)
        elif self._drag == "high":
            self._move_handle("high", high + delta)
        else:
            delta = float(np.clip(delta, -low, 1.0 - high))
            self._set_values(low + delta, high + delta)

    def _drag_tip(self, event: QtGui.QMouseEvent) -> None:
        genes = int(np.count_nonzero(genes_in_range(self.prevalence, self.low, self.high)))
        QtWidgets.QToolTip.showText(
            event.globalPosition().toPoint(),
            f"{percent_text(self.low)} to {percent_text(self.high)} · {genes:,} genes"
            "\nHold Shift for 0.5% steps",
            self,
        )

    def mouseMoveEvent(self, event: QtGui.QMouseEvent) -> None:
        x = event.position().x()
        if self._drag is not None and event.modifiers() & QtCore.Qt.KeyboardModifier.ShiftModifier:
            self._fine_drag(x)
            self._drag_tip(event)
        elif self._drag in ("low", "high"):
            self._fine = None
            self._move_handle(self._drag, self._value_at(x))
            self._drag_tip(event)
        elif self._drag == "range":
            start, x_high, x_low = self._drag_origin
            dx = x - start
            high, low = self._value_at(x_high + dx), self._value_at(x_low + dx - 0.01)
            self._fine = None
            if high - low >= MIN_GAP:
                self._set_values(low, high)
            self._drag_tip(event)
        else:
            handle = self._handle_at(x)
            inside = self._x_of("high") < x < self._x_of("low")
            self.setCursor(
                QtCore.Qt.CursorShape.SizeHorCursor
                if handle
                else (QtCore.Qt.CursorShape.OpenHandCursor if inside else QtCore.Qt.CursorShape.PointingHandCursor)
            )
            if self._genes():
                n = self._genes()
                area = self._area()
                index = int(np.clip((x - area.left()) / max(1.0, area.width()) * n, 0, n - 1))
                QtWidgets.QToolTip.showText(
                    event.globalPosition().toPoint(),
                    f"Prevalence {percent_text(self.prevalence[index])} (gene {index + 1:,} of {n:,})",
                    self,
                )
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event: QtGui.QMouseEvent) -> None:
        if self._drag is not None and event.button() == QtCore.Qt.MouseButton.LeftButton:
            self._drag = None
            self._fine = None
            self.update()
            self.rangeCommitted.emit(self.low, self.high)
        super().mouseReleaseEvent(event)

    def mouseDoubleClickEvent(self, event: QtGui.QMouseEvent) -> None:
        handle = self._handle_at(event.position().x())
        if handle is not None:
            self._drag = None
            self.handleDoubleClicked.emit(handle)
            event.accept()
            return
        super().mouseDoubleClickEvent(event)

    def keyPressEvent(self, event: QtGui.QKeyEvent) -> None:
        key = event.key()
        # Left is toward the core (higher prevalence), as the band is drawn.
        steps = {QtCore.Qt.Key.Key_Left: 1, QtCore.Qt.Key.Key_Right: -1}.get(key)
        if steps is None:
            super().keyPressEvent(event)
            return
        if event.modifiers() & QtCore.Qt.KeyboardModifier.ShiftModifier:
            steps *= 10
        value = self.low if self._active == "low" else self.high
        self._move_handle(self._active, value + steps * STEP)
        self.rangeCommitted.emit(self.low, self.high)
        event.accept()

    def focusInEvent(self, event: QtGui.QFocusEvent) -> None:
        self.update()
        super().focusInEvent(event)

    def focusOutEvent(self, event: QtGui.QFocusEvent) -> None:
        self.update()
        super().focusOutEvent(event)


class PrevalenceRangeBand(QtWidgets.QFrame):
    """The band with its readout ("1,234 genes · 30% to 50%"), the two spin
    boxes and the dropped isolates hint."""

    # (low, high): edited while dragging, committed on release or entry.
    rangeEdited = QtCore.pyqtSignal(float, float)
    rangeCommitted = QtCore.pyqtSignal(float, float)

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setObjectName("PrevalenceRangeBand")
        self.setFixedHeight(BAND_HEIGHT)
        self.setStyleSheet(
            f"QFrame#PrevalenceRangeBand {{ background: {COLORS['surface_subtle']};"
            f" border-top: 1px solid {COLORS['border']}; }}"
        )
        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(8, 4, 10, 4)
        layout.setSpacing(10)
        self.canvas = _BandCanvas(self)
        self.canvas.setToolTip(
            "Genes by prevalence, core on the left. Drag a handle or the range"
            " (hold Shift for 0.5% steps); arrow keys move the focused handle."
        )
        layout.addWidget(self.canvas, stretch=1)

        side = QtWidgets.QWidget(self)
        side.setMaximumWidth(206)
        side.setSizePolicy(QtWidgets.QSizePolicy.Policy.Preferred, QtWidgets.QSizePolicy.Policy.Preferred)
        side_layout = QtWidgets.QVBoxLayout(side)
        side_layout.setContentsMargins(0, 0, 0, 0)
        side_layout.setSpacing(3)
        self.readout = ElidedHint("", side)
        self.readout.setStyleSheet(f"color: {COLORS['text_strong']}; font-weight: 600;")
        side_layout.addWidget(self.readout)
        spins = QtWidgets.QHBoxLayout()
        spins.setContentsMargins(0, 0, 0, 0)
        spins.setSpacing(6)
        self.low_spin = self._spin(side, "Lowest prevalence (genes in more isolates than this)")
        self.high_spin = self._spin(side, "Highest prevalence (genes in fewer isolates than this; 100% keeps the core)")
        to_label = QtWidgets.QLabel("to", side)
        to_label.setProperty("role", "muted")
        spins.addWidget(self.low_spin)
        spins.addWidget(to_label)
        spins.addWidget(self.high_spin)
        spins.addStretch(1)
        side_layout.addLayout(spins)
        self.hint = ElidedHint("", side)
        self.hint.setStyleSheet(f"color: {COLORS['text_muted']}; font-size: 11px;")
        side_layout.addWidget(self.hint)
        layout.addWidget(side)

        self.canvas.rangeEdited.connect(self._on_canvas_edited)
        self.canvas.rangeCommitted.connect(self._on_canvas_committed)
        self.canvas.handleDoubleClicked.connect(self._focus_spin)
        self.low_spin.editingFinished.connect(self._on_spin_entered)
        self.high_spin.editingFinished.connect(self._on_spin_entered)
        self._sync_spins()
        self._update_readout()

    def _spin(self, parent: QtWidgets.QWidget, tooltip: str) -> QtWidgets.QDoubleSpinBox:
        spin = QtWidgets.QDoubleSpinBox(parent)
        spin.setRange(0.0, 100.0)
        spin.setDecimals(1)
        spin.setSingleStep(100.0 * STEP)
        spin.setSuffix("%")
        spin.setKeyboardTracking(False)
        spin.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        spin.setStyleSheet("QDoubleSpinBox { padding: 0px 4px; border-radius: 6px; }")
        # The widest text plus the arrows, the padding and the frame.
        text = QtGui.QFontMetrics(spin.font()).horizontalAdvance("100.0%")
        spin.setFixedSize(text + 34, 24)
        spin.setToolTip(tooltip)
        return spin

    # ---------------------------------------------------------- API
    def low(self) -> float:
        return self.canvas.low

    def high(self) -> float:
        return self.canvas.high

    def set_prevalence(self, prevalence: np.ndarray) -> None:
        self.canvas.set_prevalence(prevalence)
        self._update_readout()

    def set_image(self, image: Optional[QtGui.QImage]) -> None:
        self.canvas.set_image(image)

    def set_range(self, low: float, high: float) -> None:
        low, high = snap(low), snap(high)
        if high - low < MIN_GAP:
            high = min(1.0, low + MIN_GAP)
            low = high - MIN_GAP
        self.canvas.set_range(low, high)
        self._sync_spins()
        self._update_readout()

    def genes_in_range(self) -> int:
        return int(np.count_nonzero(genes_in_range(self.canvas.prevalence, self.low(), self.high())))

    def set_dropped(self, dropped: Sequence[str]) -> None:
        names: List[str] = [str(s) for s in dropped]
        if names:
            noun = "isolate has" if len(names) == 1 else "isolates have"
            self.hint.setText(f"{len(names):,} {noun} no genes in range")
            shown = names[:40]
            more = f"\n... and {len(names) - len(shown):,} more" if len(names) > len(shown) else ""
            self.hint.setToolTip("Left out of the clustering:\n" + "\n".join(shown) + more)
        else:
            self.hint.setText("")
            self.hint.setToolTip("")

    # ---------------------------------------------------------- internals
    def _sync_spins(self) -> None:
        for spin, value in ((self.low_spin, self.low()), (self.high_spin, self.high())):
            with QtCore.QSignalBlocker(spin):
                spin.setValue(round(100.0 * value, 1))

    def _update_readout(self) -> None:
        text = f"{self.genes_in_range():,} genes · {percent_text(self.low())} to {percent_text(self.high())}"
        self.readout.setText(text)
        self.readout.setToolTip(text)

    def _on_canvas_edited(self, low: float, high: float) -> None:
        self._sync_spins()
        self._update_readout()
        self.rangeEdited.emit(low, high)

    def _on_canvas_committed(self, low: float, high: float) -> None:
        self._sync_spins()
        self._update_readout()
        self.rangeCommitted.emit(low, high)

    def _on_spin_entered(self) -> None:
        low, high = snap(self.low_spin.value() / 100.0), snap(self.high_spin.value() / 100.0)
        if high - low < MIN_GAP:
            # The other bound gives way to the one just entered.
            if self.sender() is self.low_spin:
                high = min(1.0, low + MIN_GAP)
                low = high - MIN_GAP
            else:
                low = max(0.0, high - MIN_GAP)
                high = low + MIN_GAP
        if (low, high) == (self.low(), self.high()):
            self._sync_spins()
            return
        self.set_range(low, high)
        self.rangeCommitted.emit(self.low(), self.high())

    def _focus_spin(self, which: str) -> None:
        spin = self.low_spin if which == "low" else self.high_spin
        spin.setFocus(QtCore.Qt.FocusReason.OtherFocusReason)
        spin.selectAll()
