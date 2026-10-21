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
            (0.60, 20.0, 430.0, 1.00, 0.0, 20),
            # The middle wave is drawn over the ship, so it needs enough opacity to
            # actually read as water lapping over the hull when the ship sinks.
            (0.72, 28.0, 340.0, 1.35, 1.2, 70),
            (0.85, 22.0, 260.0, 0.80, 2.4, 95),
        )

        def draw_wave(layer: tuple) -> None:
            base_frac, amp, wavelength, pmult, poff, alpha = layer
            painter.fillPath(
                self._wave_path(w, h, base_frac * h, amp, wavelength, wave_phase * pmult + poff),
                QtGui.QBrush(QtGui.QColor(255, 255, 255, alpha)),
            )

        # Back wave, then a little navigator's ship riding it, then the nearer
        # waves - so the hull nestles into the water instead of floating above it.
        draw_wave(layers[0])
        self._draw_ship(painter, w, h, wave_phase, layers[0])
        for layer in layers[1:]:
            draw_wave(layer)

        # Subtle 1px inner border for definition on light desktops (straight corners).
        pen = QtGui.QPen(QtGui.QColor(0, 0, 0, 60))
        pen.setWidth(1)
        painter.setPen(pen)
        painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
        painter.drawRect(0, 0, w - 1, h - 1)
        painter.end()

    def _draw_ship(
        self, painter: QtGui.QPainter, w: int, h: int, wave_phase: float, base_layer: tuple
    ) -> None:
        """A small three-masted galleon riding the back wave.

        Drawn in the same flat, translucent-white idiom as the waves. It bobs on
        the wave surface and rocks with the wave's local slope, and the sails
        gently breathe with the wind.
        """
        _base_frac, _amp, wavelength, pmult, poff, _alpha = base_layer
        # Asymmetric sway: a quick lurch forward (~0.5 s) then a slow drift back
        # (~3 s), repeating. The ship faces right-to-left, so "forward" is left.
        t = wave_phase / self._WAVE_SPEED  # elapsed seconds
        t_fwd, t_back = 0.5, 3.0
        u = t % (t_fwd + t_back)
        if u < t_fwd:
            s = u / t_fwd
            pos = s * s * (3.0 - 2.0 * s)          # 0 -> 1, eased (fast forward)
        else:
            s = (u - t_fwd) / t_back
            pos = 1.0 - s * s * (3.0 - 2.0 * s)    # 1 -> 0, eased (slow back)
        sway = 24.0  # px to either side of centre
        drift = (0.5 - pos) * 2.0 * sway           # pos 1 -> left (forward), pos 0 -> right (back)
        ship_x = w * 0.60 + drift
        # Bob on the wave it rides, at the same average height (between the back
        # 0.60 wave and the middle 0.72 wave). The bob is deep enough that at the
        # bottom of a swell the hull sinks behind the middle wave (drawn after it).
        angle = 2.0 * math.pi * (ship_x / wavelength) + wave_phase * pmult + poff
        surface_y = 0.63 * h + 26.0 * math.sin(angle)
        rock_deg = 5.0 * math.cos(angle)

        sail_col = QtGui.QColor(255, 255, 255, 236)
        hull_col = QtGui.QColor(232, 240, 255, 232)
        mast_col = QtGui.QColor(255, 255, 255, 205)
        rig_col = QtGui.QColor(255, 255, 255, 90)
        flag_col = QtGui.QColor(255, 255, 255, 236)
        # Gunports and the stern windows read as openings: the background blue.
        port_col = QtGui.QColor(30, 58, 168, 150)
        lantern_col = QtGui.QColor(255, 236, 170, 235)

        painter.save()
        painter.translate(ship_x, surface_y)
        painter.rotate(rock_deg)
        painter.scale(0.9, 0.9)
        painter.scale(-1.0, 1.0)  # mirror so the ship sails right-to-left

        # Mast positions (x) and heads (y); the bow points to +x before the mirror.
        fore_x, main_x, mizzen_x = 50.0, 2.0, -56.0
        fore_top, main_top, mizzen_top = -128.0, -150.0, -112.0
        bowsprit_tip = QtCore.QPointF(124, -36)

        # Standing rigging (behind everything).
        rig_pen = QtGui.QPen(rig_col)
        rig_pen.setWidthF(1.1)
        painter.setPen(rig_pen)
        painter.drawLine(bowsprit_tip, QtCore.QPointF(fore_x, fore_top))      # bowsprit -> fore top
        painter.drawLine(QtCore.QPointF(fore_x, fore_top), QtCore.QPointF(main_x, main_top))
        painter.drawLine(QtCore.QPointF(main_x, main_top), QtCore.QPointF(mizzen_x, mizzen_top))
        painter.drawLine(QtCore.QPointF(mizzen_x, mizzen_top), QtCore.QPointF(-92, -36))  # mizzen top -> stern
        painter.drawLine(QtCore.QPointF(main_x, main_top), QtCore.QPointF(-30, -10))     # main shrouds
        painter.drawLine(QtCore.QPointF(fore_x, fore_top), QtCore.QPointF(70, -18))      # fore shrouds

        # Masts + bowsprit.
        mast_pen = QtGui.QPen(mast_col)
        mast_pen.setWidthF(3.2)
        mast_pen.setCapStyle(QtCore.Qt.PenCapStyle.RoundCap)
        painter.setPen(mast_pen)
        painter.drawLine(QtCore.QPointF(fore_x, -14), QtCore.QPointF(fore_x, fore_top))
        painter.drawLine(QtCore.QPointF(main_x, -2), QtCore.QPointF(main_x, main_top))
        mast_pen.setWidthF(2.6)
        painter.setPen(mast_pen)
        painter.drawLine(QtCore.QPointF(mizzen_x, -22), QtCore.QPointF(mizzen_x, mizzen_top))
        painter.drawLine(QtCore.QPointF(78, -16), bowsprit_tip)                 # bowsprit

        # Sails, billowing to leeward; the belly breathes with the wind.
        bulge = 10.0 + 2.5 * math.sin(wave_phase * 1.6)
        painter.setPen(QtCore.Qt.PenStyle.NoPen)
        painter.setBrush(sail_col)
        painter.drawPath(self._square_sail(36, 66, -120, -84, bulge * 0.8))    # fore topsail
        painter.drawPath(self._square_sail(30, 72, -78, -30, bulge))           # fore course
        painter.drawPath(self._square_sail(-16, 22, -142, -102, bulge * 0.85)) # main topsail
        painter.drawPath(self._square_sail(-26, 32, -96, -26, bulge * 1.15))   # main course
        painter.drawPath(self._lateen_sail(mizzen_x, bulge))                   # mizzen lateen
        painter.drawPath(self._square_sail(100, 118, -32, -14, bulge * 0.45))  # spritsail

        # Hull (drawn over the mast feet): stepped stern castle, waist, forecastle.
        painter.setBrush(hull_col)
        painter.drawPath(self._hull_path())
        painter.drawPath(self._beakhead_path())

        # Gunports along the side and windows in the stern castle.
        painter.setBrush(port_col)
        for x in range(-40, 56, 15):
            painter.drawRoundedRect(QtCore.QRectF(x, 3.0, 6.0, 4.5), 1.0, 1.0)
        for x in (-82, -74, -66):
            painter.drawRoundedRect(QtCore.QRectF(x, -27.0, 5.0, 5.0), 1.0, 1.0)

        # Stern lantern.
        painter.setBrush(lantern_col)
        painter.drawEllipse(QtCore.QPointF(-96, -44), 3.2, 3.2)

        # Fluttering pennants at every masthead (the main one longest).
        flutter = 4.0 * math.sin(wave_phase * 3.0)
        painter.setBrush(flag_col)
        for x, top, length in ((main_x, main_top, 32.0), (fore_x, fore_top, 22.0), (mizzen_x, mizzen_top, 18.0)):
            flag = QtGui.QPainterPath()
            flag.moveTo(x, top)
            flag.lineTo(x + length + flutter, top + 4 + flutter * 0.4)
            flag.lineTo(x, top + 9)
            flag.closeSubpath()
            painter.drawPath(flag)

        painter.restore()

    @staticmethod
    def _hull_path() -> QtGui.QPainterPath:
        """A galleon hull: high stepped stern castle (poop and quarterdeck), low
        waist, raised forecastle and a sharp stem."""
        p = QtGui.QPainterPath()
        p.moveTo(-98, -40)                       # poop deck, overhanging stern
        p.lineTo(-88, -6)                        # raked transom
        p.lineTo(-80, 12)                        # stern down into the water
        p.cubicTo(-40, 32, 40, 32, 74, 14)       # deep belly
        p.lineTo(86, -10)                        # stem up to the bow
        p.lineTo(82, -22)                        # forecastle front
        p.lineTo(58, -24)                        # forecastle deck
        p.lineTo(54, -10)                        # step down to the waist
        p.cubicTo(24, -6, -16, -6, -40, -10)     # waist sheer
        p.lineTo(-44, -26)                       # step up to the quarterdeck
        p.lineTo(-70, -30)                       # quarterdeck
        p.lineTo(-74, -40)                       # step up to the poop
        p.closeSubpath()
        return p

    @staticmethod
    def _beakhead_path() -> QtGui.QPainterPath:
        p = QtGui.QPainterPath()
        p.moveTo(84, -12)
        p.lineTo(104, -8)
        p.lineTo(80, 0)
        p.closeSubpath()
        return p

    @staticmethod
    def _lateen_sail(mast_x: float, bulge: float) -> QtGui.QPainterPath:
        """The mizzen's triangular lateen sail on its slanted yard."""
        top = QtCore.QPointF(mast_x + 22, -110)
        low = QtCore.QPointF(mast_x - 30, -48)
        clew = QtCore.QPointF(mast_x + 6, -28)
        p = QtGui.QPainterPath()
        p.moveTo(top)
        p.lineTo(low)                                                       # yard
        p.quadTo(QtCore.QPointF(mast_x - 4 + bulge * 0.4, -30 + bulge * 0.3), clew)   # foot
        p.quadTo(QtCore.QPointF(mast_x + 22 + bulge * 0.6, -66), top)       # leech
        p.closeSubpath()
        return p

    @staticmethod
    def _square_sail(
        xl: float, xr: float, ytop: float, ybot: float, bulge: float
    ) -> QtGui.QPainterPath:
        """A billowing square sail bulging to leeward (down/right)."""
        span = ybot - ytop
        p = QtGui.QPainterPath()
        p.moveTo(xl, ytop)
        p.lineTo(xr, ytop)                                                  # yard (top)
        p.cubicTo(xr + bulge * 0.5, ytop + span * 0.4,
                  xr + bulge * 0.5, ytop + span * 0.7,
                  xr - bulge * 0.1, ybot)                                   # leech (right)
        p.cubicTo((xl + xr) * 0.5 + bulge, ybot + bulge * 0.55,
                  (xl + xr) * 0.5 - bulge * 0.2, ybot + bulge * 0.35,
                  xl + bulge * 0.1, ybot)                                   # foot (bottom belly)
        p.cubicTo(xl + bulge * 0.35, ytop + span * 0.6,
                  xl + bulge * 0.25, ytop + span * 0.3,
                  xl, ytop)                                                 # luff (left)
        p.closeSubpath()
        return p

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
