"""Borderless, straight-cornered loading splash shown while a run view loads.

A bigger splash with an interesting blue gradient background, gently animated
waves, and a thin blue "loading line" sweeping left to right along the bottom
(the sweeping accent seen in Adobe / Geneious splash screens).

The line and waves derive their phase from elapsed wall-clock time (not a
per-tick counter), so they keep their correct speed even when the frame rate
drops under load. The heavy run load runs on a background thread (see app.py) so
the event loop stays live and the animation actually plays; only the final
on-GUI-thread window construction briefly pauses it just before the splash
closes.
"""
from __future__ import annotations

import math
from typing import Optional

from PyQt6 import QtCore, QtGui, QtWidgets


class LoadingLine(QtWidgets.QWidget):
    """An indeterminate marquee: a soft blue segment sweeping left to right."""

    def __init__(
        self,
        parent: Optional[QtWidgets.QWidget] = None,
        *,
        track: str = "#0b1533",
        track_alpha: int = 90,
        accent: str = "#cfe3ff",
        height: int = 4,
        period_ms: int = 1150,
    ) -> None:
        super().__init__(parent)
        self._track = QtGui.QColor(track)
        self._track.setAlpha(track_alpha)
        self._accent = QtGui.QColor(accent)
        self.setFixedHeight(height)
        self.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed,
        )
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_TransparentForMouseEvents, True)
        self._period_ms = max(1, period_ms)
        self._phase = 0.0  # manual fallback (used for static captures)
        self._elapsed = QtCore.QElapsedTimer()
        self._timer = QtCore.QTimer(self)
        self._timer.setInterval(16)
        self._timer.timeout.connect(self.update)

    def start(self) -> None:
        if not self._elapsed.isValid():
            self._elapsed.start()
        if not self._timer.isActive():
            self._timer.start()

    def stop(self) -> None:
        self._timer.stop()

    def _current_phase(self) -> float:
        if self._elapsed.isValid():
            return (self._elapsed.elapsed() % self._period_ms) / self._period_ms
        return self._phase

    def paintEvent(self, _event: QtGui.QPaintEvent) -> None:
        w = self.width()
        h = self.height()
        if w <= 0 or h <= 0:
            return
        painter = QtGui.QPainter(self)
        painter.fillRect(0, 0, w, h, self._track)

        seg = max(80, int(w * 0.30))
        # Segment travels fully off both ends, so the sweep reads continuous.
        x = -seg + self._current_phase() * (w + seg)
        gradient = QtGui.QLinearGradient(x, 0.0, x + seg, 0.0)
        edge = QtGui.QColor(self._accent)
        edge.setAlpha(0)
        core = QtGui.QColor(self._accent)
        gradient.setColorAt(0.0, edge)
        gradient.setColorAt(0.5, core)
        gradient.setColorAt(1.0, edge)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        painter.fillRect(QtCore.QRectF(x, 0.0, float(seg), float(h)), QtGui.QBrush(gradient))
        painter.end()


class LoadingSplash(QtWidgets.QWidget):
    """Frameless, opaque, straight-cornered splash: blue gradient + waves + line."""

    _WAVE_SPEED = 0.9  # radians per second of drift

    def __init__(
        self,
        *,
        title: str = "BACTolomew",
        subtitle: str = "ANNOTATION EXPLORER",
        status: str = "Loading…",
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(
            parent,
            QtCore.Qt.WindowType.FramelessWindowHint
            | QtCore.Qt.WindowType.Dialog
            | QtCore.Qt.WindowType.WindowStaysOnTopHint,
        )
        self.setObjectName("LoadingSplash")
        # No stylesheet background on the splash itself: paintEvent draws the
        # gradient + waves. Labels are styled individually (transparent bg) so
        # the painted background shows through.
        self._closing = False
        self._wave_phase = 0.0  # manual fallback (used for static captures)
        self._wave_elapsed = QtCore.QElapsedTimer()

        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        body = QtWidgets.QWidget()
        body.setAttribute(QtCore.Qt.WidgetAttribute.WA_TranslucentBackground, True)
        body_layout = QtWidgets.QVBoxLayout(body)
        body_layout.setContentsMargins(46, 42, 46, 34)
        body_layout.setSpacing(8)

        self._title_label = QtWidgets.QLabel(title)
        self._title_label.setStyleSheet(
            "color: #ffffff; font-size: 36px; font-weight: 800; background: transparent;"
        )
        body_layout.addWidget(self._title_label)

        self._subtitle_label = QtWidgets.QLabel(subtitle)
        self._subtitle_label.setStyleSheet(
            "color: #bcd0f7; font-size: 11px; font-weight: 700; background: transparent;"
        )
        sub_font = self._subtitle_label.font()
        sub_font.setLetterSpacing(QtGui.QFont.SpacingType.AbsoluteSpacing, 3.0)
        self._subtitle_label.setFont(sub_font)
        body_layout.addWidget(self._subtitle_label)

        body_layout.addStretch(1)

        self._status_label = QtWidgets.QLabel(status)
        self._status_label.setStyleSheet(
            "color: #d3ddf7; font-size: 13px; background: transparent;"
        )
        body_layout.addWidget(self._status_label)

        outer.addWidget(body, 1)

        self._line = LoadingLine(self)
        outer.addWidget(self._line)

        self.setFixedSize(700, 410)

        # Repaint the waves ~30fps; the actual drift is time-based (see paintEvent).
        self._wave_timer = QtCore.QTimer(self)
        self._wave_timer.setInterval(33)
        self._wave_timer.timeout.connect(self.update)

    # ---- API ----
    def start(self) -> None:
        self._center_on_screen()
        self.show()
        self.raise_()
        self._line.start()
        if not self._wave_elapsed.isValid():
            self._wave_elapsed.start()
        self._wave_timer.start()
        QtWidgets.QApplication.processEvents()

    def set_status(self, text: str) -> None:
        """Update the status line. The event loop repaints it (and the animation)."""
        self._status_label.setText(text)

    def finish(self) -> None:
        if self._closing:
            return
        self._closing = True
        self._line.stop()
        self._wave_timer.stop()
        self.close()
        self.deleteLater()

    # ---- internals ----
    def _center_on_screen(self) -> None:
        screen = self.screen() or QtWidgets.QApplication.primaryScreen()
        if screen is None:
            return
        geo = screen.availableGeometry()
        self.move(
            geo.center().x() - self.width() // 2,
            geo.center().y() - self.height() // 2,
        )

    def _current_wave_phase(self) -> float:
        if self._wave_elapsed.isValid():
            return (self._wave_elapsed.elapsed() / 1000.0) * self._WAVE_SPEED
        return self._wave_phase

    def paintEvent(self, _event: QtGui.QPaintEvent) -> None:
        w = self.width()
        h = self.height()
        if w <= 0 or h <= 0:
            return
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)

        # Diagonal blue gradient background.
        grad = QtGui.QLinearGradient(0.0, 0.0, float(w), float(h))
        grad.setColorAt(0.0, QtGui.QColor("#101a4a"))
        grad.setColorAt(0.55, QtGui.QColor("#1e3aa8"))
        grad.setColorAt(1.0, QtGui.QColor("#2f6fed"))
        painter.fillRect(self.rect(), QtGui.QBrush(grad))

        # Layered, gently drifting waves in the lower portion.
        wave_phase = self._current_wave_phase()
        layers = (
            # (baseline_frac, amplitude, wavelength, phase_mult, phase_off, alpha)
            (0.60, 20.0, 430.0, 1.00, 0.0, 16),
            (0.72, 28.0, 340.0, 1.35, 1.2, 22),
            (0.85, 22.0, 260.0, 0.80, 2.4, 30),
        )
        for base_frac, amp, wavelength, pmult, poff, alpha in layers:
            painter.fillPath(
                self._wave_path(w, h, base_frac * h, amp, wavelength, wave_phase * pmult + poff),
                QtGui.QBrush(QtGui.QColor(255, 255, 255, alpha)),
            )

        # Subtle 1px inner border for definition on light desktops (straight corners).
        pen = QtGui.QPen(QtGui.QColor(0, 0, 0, 60))
        pen.setWidth(1)
        painter.setPen(pen)
        painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
        painter.drawRect(0, 0, w - 1, h - 1)
        painter.end()

    @staticmethod
    def _wave_path(
        w: int, h: int, baseline: float, amplitude: float, wavelength: float, phase: float
    ) -> QtGui.QPainterPath:
        path = QtGui.QPainterPath()
        path.moveTo(0.0, float(h))
        path.lineTo(0.0, baseline)
        step = 6.0
        x = 0.0
        two_pi = 2.0 * math.pi
        while x <= w:
            y = baseline + amplitude * math.sin(two_pi * (x / wavelength) + phase)
            path.lineTo(x, y)
            x += step
        path.lineTo(float(w), baseline + amplitude * math.sin(two_pi * (w / wavelength) + phase))
        path.lineTo(float(w), float(h))
        path.closeSubpath()
        return path
