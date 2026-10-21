#!/usr/bin/env python
"""BACTolomew Map Lab - a tiny app to compare the candidate map renderers LIVE and
to crash-test them the same way the real app builds them.

    python map_tests/map_lab.py

Each renderer opens in its OWN subprocess (lab_viewer.py), so a GPU/GL segfault
only kills that window - the Lab keeps running and flags it 💥.
"""
from __future__ import annotations
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# The live widgets need the APP's env (PyQt6/matplotlib/etc.) and must run exactly
# as the app builds them - that's where a GL crash surfaces. Re-exec under it if
# the current interpreter can't even import PyQt6 (e.g. launched from base).
BACTOLOMEW_PY = os.environ.get("BACTOLOMEW_PY") or "/home/bals/miniforge3/envs/bactolomew-qt/bin/python"
if not os.path.exists(BACTOLOMEW_PY):
    BACTOLOMEW_PY = sys.executable
try:
    import PyQt6  # noqa: F401
except Exception:
    if os.path.exists(BACTOLOMEW_PY) and os.path.realpath(sys.executable) != os.path.realpath(BACTOLOMEW_PY):
        os.execv(BACTOLOMEW_PY, [BACTOLOMEW_PY, os.path.abspath(__file__), *sys.argv[1:]])
    raise

from PyQt6 import QtCore, QtGui, QtWidgets  # noqa: E402
import live_backends  # noqa: E402

VIEWER = os.path.join(HERE, "lab_viewer.py")
OUT = os.path.join(HERE, "out")

TAGS = {
    "safe": ("#0f9d6b", "No GL · safe"),
    "gl": ("#d97706", "Needs GPU GL"),
}


class Card(QtWidgets.QFrame):
    def __init__(self, name: str) -> None:
        super().__init__()
        self.name = name
        self.info = live_backends.meta(name) or {"label": name, "tag": "safe", "desc": ""}
        self.proc: QtCore.QProcess | None = None
        self._auto = False
        self.setObjectName("card")
        self.setStyleSheet(
            "QFrame#card{background:#ffffff;border:1px solid #d9e2ec;border-radius:14px;}"
        )
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(14, 14, 14, 14)
        lay.setSpacing(9)

        thumb = QtWidgets.QLabel()
        thumb.setFixedHeight(150)
        thumb.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        thumb.setStyleSheet("background:#f4f8fc;border:1px solid #e7ebf1;border-radius:8px;")
        png = os.path.join(OUT, f"{name}.png")
        if os.path.exists(png):
            pix = QtGui.QPixmap(png).scaledToHeight(148, QtCore.Qt.TransformationMode.SmoothTransformation)
            thumb.setPixmap(pix)
        else:
            thumb.setText("(run test_maps.py\nfor a thumbnail)")
            thumb.setStyleSheet(thumb.styleSheet() + "color:#94a3b8;")
        lay.addWidget(thumb)

        header = QtWidgets.QHBoxLayout()
        title = QtWidgets.QLabel(self.info["label"])
        title.setStyleSheet("font-size:14px;font-weight:700;color:#172433;")
        header.addWidget(title)
        header.addStretch(1)
        color, text = TAGS.get(self.info["tag"], TAGS["safe"])
        chip = QtWidgets.QLabel(text)
        chip.setStyleSheet(
            f"color:#fff;background:{color};border-radius:9px;padding:2px 9px;font-size:10px;font-weight:700;"
        )
        header.addWidget(chip)
        lay.addLayout(header)

        desc = QtWidgets.QLabel(self.info["desc"])
        desc.setWordWrap(True)
        desc.setStyleSheet("color:#566577;font-size:12px;")
        lay.addWidget(desc)

        row = QtWidgets.QHBoxLayout()
        self.open_btn = QtWidgets.QPushButton("Open live  ▶")
        self.open_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.open_btn.setStyleSheet(
            "QPushButton{background:#2f6fed;color:#fff;border:none;border-radius:8px;padding:7px 14px;font-weight:600;}"
            "QPushButton:hover{background:#295fd0;}"
        )
        self.open_btn.clicked.connect(self.open_live)
        row.addWidget(self.open_btn)
        self.status = QtWidgets.QLabel("not run")
        self.status.setStyleSheet("color:#64748b;font-size:12px;")
        row.addWidget(self.status)
        row.addStretch(1)
        lay.addLayout(row)

    def _set(self, text: str, color: str) -> None:
        self.status.setText(text)
        self.status.setStyleSheet(f"color:{color};font-size:12px;font-weight:600;")

    def open_live(self, *, auto: bool = False, grab: str | None = None) -> None:
        if self.proc is not None and self.proc.state() != QtCore.QProcess.ProcessState.NotRunning:
            return
        self._auto = auto
        self.proc = QtCore.QProcess(self)
        args = [VIEWER, self.name]
        if grab:
            args += ["--grab", grab]
        self.proc.setProgram(BACTOLOMEW_PY)
        self.proc.setArguments(args)
        self.proc.setWorkingDirectory(HERE)
        self.proc.finished.connect(self._finished)
        self.proc.errorOccurred.connect(self._errored)
        self._set("opening …", "#64748b")
        self.proc.start()
        if not auto:
            QtCore.QTimer.singleShot(3500, self._still_open)

    def _still_open(self) -> None:
        if self.proc is not None and self.proc.state() == QtCore.QProcess.ProcessState.Running:
            self._set("✅ open - interact, then close", "#0f9d6b")

    def _finished(self, code: int, status: QtCore.QProcess.ExitStatus) -> None:
        crashed = status == QtCore.QProcess.ExitStatus.CrashExit
        if crashed:
            self._set("💥 crashed (GL)", "#dc2626")
        elif self._auto:
            self._set("✅ renders - no crash", "#0f9d6b")
        else:
            self._set("window closed", "#64748b")
        done_cb = getattr(self, "_on_auto_done", None)
        if self._auto and callable(done_cb):
            done_cb(self.name, crashed)

    def _errored(self, _err: QtCore.QProcess.ProcessError) -> None:
        self._set("⚠️ failed to launch", "#d97706")


class MapLab(QtWidgets.QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("BACTolomew Map Lab")
        self.resize(880, 720)
        central = QtWidgets.QWidget()
        central.setStyleSheet("background:#eef3f8;")
        self.setCentralWidget(central)
        outer = QtWidgets.QVBoxLayout(central)
        outer.setContentsMargins(22, 20, 22, 20)
        outer.setSpacing(14)

        title = QtWidgets.QLabel("Map renderer lab")
        title.setStyleSheet("font-size:22px;font-weight:800;color:#172433;")
        outer.addWidget(title)
        sub = QtWidgets.QLabel(
            "Open each renderer live (own window, real pan/zoom). Each runs in its own "
            "process, so a GPU/GL crash only closes that window and is flagged below."
        )
        sub.setWordWrap(True)
        sub.setStyleSheet("color:#566577;font-size:13px;")
        outer.addWidget(sub)

        controls = QtWidgets.QHBoxLayout()
        self.safe_btn = QtWidgets.QPushButton("Open all safe maps")
        self.crash_btn = QtWidgets.QPushButton("Auto crash-test all  ⚡")
        for b in (self.safe_btn, self.crash_btn):
            b.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
            b.setStyleSheet(
                "QPushButton{background:#ffffff;border:1px solid #cdd6e0;border-radius:8px;padding:8px 14px;font-weight:600;color:#172433;}"
                "QPushButton:hover{background:#f4f7fb;}"
            )
        self.safe_btn.clicked.connect(self._open_all_safe)
        self.crash_btn.clicked.connect(self._crash_test_all)
        controls.addWidget(self.safe_btn)
        controls.addWidget(self.crash_btn)
        controls.addStretch(1)
        outer.addLayout(controls)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        holder = QtWidgets.QWidget()
        grid = QtWidgets.QGridLayout(holder)
        grid.setSpacing(16)
        grid.setContentsMargins(0, 4, 0, 0)
        self.cards: dict[str, Card] = {}
        for i, name in enumerate(live_backends.NAMES):
            card = Card(name)
            self.cards[name] = card
            grid.addWidget(card, i // 2, i % 2)
        scroll.setWidget(holder)
        outer.addWidget(scroll, stretch=1)

        self._auto_queue: list[str] = []

    def _open_all_safe(self) -> None:
        for name in live_backends.NAMES:
            if (live_backends.meta(name) or {}).get("tag") == "safe":
                self.cards[name].open_live()

    def _crash_test_all(self) -> None:
        # sequential: run one, then the next when it finishes (via _on_auto_done)
        self._auto_queue = list(live_backends.NAMES)
        self.crash_btn.setEnabled(False)
        for card in self.cards.values():
            card._on_auto_done = self._auto_next
        self._auto_next(None, False)

    def _auto_next(self, _name, _crashed) -> None:
        if not self._auto_queue:
            self.crash_btn.setEnabled(True)
            return
        name = self._auto_queue.pop(0)
        grab = os.path.join(OUT, f"live_{name}.png")
        self.cards[name].open_live(auto=True, grab=grab)


def main() -> None:
    app = QtWidgets.QApplication(sys.argv)
    lab = MapLab()
    lab.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
