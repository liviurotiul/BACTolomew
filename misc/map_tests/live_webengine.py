"""LIVE QtWebEngine map widget for the BACTolomew Map Lab.

Reuses backend_webengine's self-contained inline-SVG choropleth (no network)
and wraps it in a QWebEngineView, adding client-side JavaScript so the map is
INTERACTIVE: drag to pan, mouse-wheel to zoom (both implemented by mutating the
SVG viewBox, so no external map tiles / network are ever needed).

Unlike the static backend this is a live, on-screen widget: build_widget()
returns a visible QWidget. WA_DontShowOnScreen is deliberately NOT set.

May crash on broken-GL hosts (expected; isolated by the Lab). GPU is disabled
and software OpenGL forced up front so it renders under xvfb.
"""
from __future__ import annotations

# Use the real GL stack (GPU once the driver is fixed; software under xvfb).
import os
os.environ.setdefault("QTWEBENGINE_CHROMIUM_FLAGS", "--no-sandbox")

from PyQt6 import QtCore, QtWidgets

# WebEngine requires AA_ShareOpenGLContexts to be set (Qt only checks that the
# flag IS set, in 6.10 even after the QApplication exists). The harness test
# creates the QApplication *before* importing this module, so we must set the
# attribute here BEFORE any QtWebEngine import (backend_webengine imports it).
QtWidgets.QApplication.setAttribute(
    QtCore.Qt.ApplicationAttribute.AA_ShareOpenGLContexts, True
)

from PyQt6.QtWebEngineWidgets import QWebEngineView  # noqa: E402

import fixture  # noqa: E402
import geo_data  # noqa: E402,F401  (used transitively by backend)
import backend_webengine  # noqa: E402


# JS that turns the static full-page <svg> into a pan + wheel-zoom surface by
# mutating its viewBox. Self-contained, no network.
_INTERACT_JS = """
<script>
(function () {
  function init() {
    var svg = document.querySelector('svg');
    if (!svg) { return; }
    var W = %d, H = %d;
    var vb = { x: 0, y: 0, w: W, h: H };
    function apply() {
      svg.setAttribute('viewBox',
        vb.x + ' ' + vb.y + ' ' + vb.w + ' ' + vb.h);
    }
    apply();
    svg.style.cursor = 'grab';

    // wheel zoom, centred on the cursor
    svg.addEventListener('wheel', function (e) {
      e.preventDefault();
      var rect = svg.getBoundingClientRect();
      var mx = (e.clientX - rect.left) / rect.width;
      var my = (e.clientY - rect.top) / rect.height;
      var factor = e.deltaY < 0 ? 0.88 : 1.0 / 0.88;
      var nw = vb.w * factor, nh = vb.h * factor;
      // clamp zoom range
      var minW = W * 0.05, maxW = W * 3.0;
      if (nw < minW) { factor = minW / vb.w; nw = minW; nh = vb.h * factor; }
      if (nw > maxW) { factor = maxW / vb.w; nw = maxW; nh = vb.h * factor; }
      vb.x += (vb.w - nw) * mx;
      vb.y += (vb.h - nh) * my;
      vb.w = nw; vb.h = nh;
      apply();
    }, { passive: false });

    // drag to pan
    var dragging = false, lastX = 0, lastY = 0;
    svg.addEventListener('mousedown', function (e) {
      dragging = true; lastX = e.clientX; lastY = e.clientY;
      svg.style.cursor = 'grabbing';
    });
    window.addEventListener('mouseup', function () {
      dragging = false; svg.style.cursor = 'grab';
    });
    window.addEventListener('mousemove', function (e) {
      if (!dragging) { return; }
      var rect = svg.getBoundingClientRect();
      var dx = (e.clientX - lastX) / rect.width * vb.w;
      var dy = (e.clientY - lastY) / rect.height * vb.h;
      vb.x -= dx; vb.y -= dy;
      lastX = e.clientX; lastY = e.clientY;
      apply();
    });
  }
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else { init(); }
})();
</script>
"""


def _build_interactive_html() -> str:
    """The backend's self-contained HTML + interaction JS, full-viewport."""
    svg = backend_webengine._build_svg()
    js = _INTERACT_JS % (fixture.WIDTH, fixture.HEIGHT)
    return (
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<style>html,body{margin:0;padding:0;height:100%;background:"
        f"{geo_data.OCEAN};" "}"
        "svg{display:block;width:100vw;height:100vh}</style></head><body>"
        f"{svg}{js}</body></html>"
    )


def _pump(app, ms: int) -> None:
    """Process events for a real wall-clock duration so WebEngine can paint."""
    t = QtCore.QElapsedTimer()
    t.start()
    while t.elapsed() < ms:
        app.processEvents(QtCore.QEventLoop.ProcessEventsFlag.AllEvents, 20)


def build_widget(parent=None):
    """Return a live, on-screen, interactive QWebEngineView showing the map.

    WebEngine loads HTML asynchronously, so we block here until loadFinished
    (with a safety timeout) and pump real-time events, guaranteeing the map is
    painted by the time this returns -- otherwise a caller that only spins
    processEvents() instantaneously would grab a blank page.
    """
    app = QtWidgets.QApplication.instance()
    if app is None:  # WebEngine requires a QApplication with non-empty argv
        app = QtWidgets.QApplication(["x"])

    view = QWebEngineView(parent)
    view.setMinimumSize(480, 300)
    # Lay the page out at a full default size up front, so a caller that grabs
    # immediately after an instantaneous resize still gets a full-frame render.
    view.resize(1000, 560)

    loop = QtCore.QEventLoop()
    done = {"v": False}

    def _on_load(_ok: bool) -> None:
        done["v"] = True
        if loop.isRunning():
            loop.quit()

    view.loadFinished.connect(_on_load)

    timer = QtCore.QTimer()
    timer.setSingleShot(True)
    timer.timeout.connect(lambda: loop.quit() if loop.isRunning() else None)
    timer.start(8000)

    view.setHtml(_build_interactive_html())

    # Show so the compositor is active, then wait for the load to finish.
    view.show()
    if not done["v"]:
        loop.exec()
    timer.stop()

    # Let the freshly-loaded content actually paint to the surface.
    _pump(app, 600)
    return view


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv or ["x"])
    w = build_widget()
    w.resize(1000, 560)
    w.show()
    sys.exit(app.exec())
