#!/usr/bin/env python
"""Empirically test whether ANY software-GL configuration makes QtWebEngine
render on this machine (instead of crashing or painting black).

    python map_tests/webengine_gl_test.py

Each config runs in its own subprocess (crash-isolated). We load a bright test
page in a real QWebEngineView, screenshot it, and classify the result:
  CRASH - the process segfaulted (GL context creation killed it)
  BLACK - it ran but painted (near) black  (Chromium software-compositing bug)
  RENDERS - it painted the bright page  ->  WebEngine is usable with this config!
"""
from __future__ import annotations
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BACTOLOMEW_PY = os.environ.get("BACTOLOMEW_PY") or "/home/bals/miniforge3/envs/bactolomew-qt/bin/python"
if not os.path.exists(BACTOLOMEW_PY):
    BACTOLOMEW_PY = sys.executable

# name -> {env: {..extra env..}, attrs: [Qt ApplicationAttribute names]}
CONFIGS = {
    "default (as the app does)": {},
    "chromium --disable-gpu": {"env": {"QTWEBENGINE_CHROMIUM_FLAGS": "--disable-gpu --no-sandbox"}},
    "Mesa llvmpipe (LIBGL_ALWAYS_SOFTWARE=1)": {
        "env": {"LIBGL_ALWAYS_SOFTWARE": "1", "GALLIUM_DRIVER": "llvmpipe"}},
    "Mesa llvmpipe + Qt software GL": {
        "env": {"LIBGL_ALWAYS_SOFTWARE": "1", "GALLIUM_DRIVER": "llvmpipe", "QT_OPENGL": "software"},
        "attrs": ["AA_UseSoftwareOpenGL"]},
    "Chromium SwiftShader (--use-gl=swiftshader)": {
        "env": {"QTWEBENGINE_CHROMIUM_FLAGS": "--use-gl=swiftshader --no-sandbox --in-process-gpu"}},
    "Chromium ANGLE+SwiftShader": {
        "env": {"QTWEBENGINE_CHROMIUM_FLAGS":
                "--use-gl=angle --use-angle=swiftshader --enable-unsafe-swiftshader --no-sandbox --in-process-gpu"}},
    "Mesa llvmpipe + SwiftShader": {
        "env": {"LIBGL_ALWAYS_SOFTWARE": "1", "GALLIUM_DRIVER": "llvmpipe",
                "QTWEBENGINE_CHROMIUM_FLAGS": "--use-gl=angle --use-angle=swiftshader --enable-unsafe-swiftshader --no-sandbox"}},
}

_HTML = (
    "<html><body style='margin:0;background:#2563eb;color:#fff;font:28px sans-serif;'>"
    "<div style='padding:40px'>WEBENGINE RENDER TEST<br><span style='font-size:16px'>"
    "if you can read this, WebEngine renders</span></div></body></html>"
)


def _worker(config_name: str, out_png: str) -> int:
    import faulthandler
    faulthandler.enable()
    from PyQt6 import QtCore, QtWidgets
    cfg = CONFIGS.get(config_name, {})
    # WebEngine always needs shared GL contexts; plus any config-specific attrs.
    QtWidgets.QApplication.setAttribute(QtCore.Qt.ApplicationAttribute.AA_ShareOpenGLContexts, True)
    for attr in cfg.get("attrs", []):
        QtWidgets.QApplication.setAttribute(getattr(QtCore.Qt.ApplicationAttribute, attr), True)
    app = QtWidgets.QApplication(["gl-test", config_name])
    from PyQt6.QtWebEngineWidgets import QWebEngineView
    view = QWebEngineView()
    view.resize(420, 300)
    view.show()
    view.setHtml(_HTML)

    def _shot() -> None:
        app.processEvents()
        view.grab().save(out_png)
        app.quit()

    QtCore.QTimer.singleShot(3500, _shot)
    return app.exec()


def _classify_png(path: str) -> str:
    try:
        from PIL import Image
        import statistics
        im = Image.open(path).convert("L")
        px = list(im.getdata())
        mean = statistics.fmean(px)
        return "RENDERS" if mean > 40 else "BLACK"
    except Exception:
        return "BLACK"


def main() -> None:
    if len(sys.argv) >= 4 and sys.argv[1] == "--worker":
        sys.exit(_worker(sys.argv[2], sys.argv[3]))

    import subprocess
    outdir = os.path.join(HERE, "out")
    os.makedirs(outdir, exist_ok=True)
    print("\nQtWebEngine software-GL probe  (each config isolated in a subprocess)")
    print("=" * 70)
    results = []
    for name, cfg in CONFIGS.items():
        out_png = os.path.join(outdir, "glprobe.png")
        if os.path.exists(out_png):
            os.remove(out_png)
        env = dict(os.environ)
        env.update(cfg.get("env", {}))
        print(f"  testing: {name} …", flush=True)
        try:
            proc = subprocess.run(
                [BACTOLOMEW_PY, os.path.abspath(__file__), "--worker", name, out_png],
                env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=60,
            )
            if proc.returncode < 0:
                verdict = "CRASH"
            elif os.path.exists(out_png) and os.path.getsize(out_png) > 500:
                verdict = _classify_png(out_png)
            else:
                verdict = "NO-OUTPUT"
        except subprocess.TimeoutExpired:
            verdict = "TIMEOUT"
        results.append((name, verdict))

    print("\n" + "=" * 70 + "\nRESULTS\n" + "=" * 70)
    icon = {"RENDERS": "✅", "CRASH": "💥", "BLACK": "⬛", "TIMEOUT": "⏱", "NO-OUTPUT": "⚠️"}
    for name, verdict in results:
        print(f"  {icon.get(verdict,'?')} {verdict:9} {name}")
    winners = [n for n, v in results if v == "RENDERS"]
    print("\nWebEngine RENDERS with:", ", ".join(winners) if winners else "(none - no software-GL config worked)")
    if winners:
        print("→ We can use the real interactive web map (MapLibre/plotly) with that config.\n")
    else:
        print("→ WebEngine is not viable on this machine; the native tile map is the way.\n")


if __name__ == "__main__":
    main()
