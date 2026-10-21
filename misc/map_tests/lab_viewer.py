"""Show ONE live map renderer in its own window - this is the isolated process
the Map Lab spawns, so a GL segfault here only kills this window.

Usage:
    python lab_viewer.py <name> [--grab OUT.png]
"""
from __future__ import annotations
import faulthandler
import os
import sys

faulthandler.enable()
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
sys.path.insert(0, HERE)

from PyQt6 import QtCore, QtWidgets  # noqa: E402
import live_backends  # noqa: E402


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: lab_viewer.py <name> [--grab OUT.png]", file=sys.stderr)
        return 2
    name = sys.argv[1]
    grab = None
    if "--grab" in sys.argv:
        idx = sys.argv.index("--grab")
        if idx + 1 < len(sys.argv):
            grab = sys.argv[idx + 1]

    # The QtWebEngine renderers need AA_ShareOpenGLContexts set BEFORE the
    # QApplication is created (they import QtWebEngine at module load). Use the
    # real GL stack (GPU once the driver is fixed; software under xvfb) - no
    # --disable-gpu, which only painted black.
    if name in ("webengine", "plotly", "maplibre"):
        # --ignore-gpu-blocklist lets WebGL (MapLibre) run even on software GL,
        # so it also works on a GPU-less box; a real GPU is used when present.
        os.environ.setdefault("QTWEBENGINE_CHROMIUM_FLAGS", "--no-sandbox --ignore-gpu-blocklist")
        QtWidgets.QApplication.setAttribute(
            QtCore.Qt.ApplicationAttribute.AA_ShareOpenGLContexts, True
        )

    # non-empty argv is required by WebEngine backends
    app = QtWidgets.QApplication(["map-lab-viewer", name])
    info = live_backends.meta(name) or {"label": name}
    win = QtWidgets.QMainWindow()
    win.setWindowTitle(f"Map Lab - {info['label']}")
    win.resize(1120, 660)

    try:
        widget = live_backends.build(name, win)
        win.setCentralWidget(widget)
    except Exception as exc:  # noqa: BLE001
        import traceback

        traceback.print_exc()
        label = QtWidgets.QLabel(f"'{name}' failed to build:\n\n{exc}")
        label.setWordWrap(True)
        label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        label.setStyleSheet("padding:24px; color:#b91c1c; font-size:13px;")
        win.setCentralWidget(label)

    win.show()

    if grab:
        # auto-screenshot mode (thumbnails / automated crash test): render, save, quit
        def _shot() -> None:
            app.processEvents()
            win.grab().save(grab)
            app.quit()

        QtCore.QTimer.singleShot(2600, _shot)

    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
