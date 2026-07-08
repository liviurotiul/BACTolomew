#!/usr/bin/env python
"""Run every candidate map backend in an ISOLATED subprocess and report which
ones render vs. crash on THIS machine — without opening the whole app.

Usage:
    python map_tests/test_maps.py

Each backend writes map_tests/out/<name>.png. A backend that segfaults (e.g. on
broken GPU OpenGL) only kills its own subprocess; this harness keeps going.
"""
from __future__ import annotations
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)

# Run the backends with the APP's env (bactolomew-qt) so they have PyQt6 /
# matplotlib / cartopy etc. AND run in the same environment the app uses -- that
# is where a GL crash would actually happen. Falls back to whatever launched us.
# Override with BACTOLOMEW_PY=/path/to/python if your env lives elsewhere.
_CANDIDATE_PY = os.environ.get("BACTOLOMEW_PY") or "/home/bals/miniforge3/envs/bactolomew-qt/bin/python"
BACKEND_PY = _CANDIDATE_PY if os.path.exists(_CANDIDATE_PY) else sys.executable

# (name, script, human description)
BACKENDS = [
    ("raster", "backend_raster.py", "Qt QPainter map — no GL, no deps"),
    ("svg", "backend_svg.py", "QtSvg vector map — no GL, no deps"),
    ("matplotlib", "backend_matplotlib.py", "matplotlib Agg — no GL"),
    ("cartopy", "backend_cartopy.py", "matplotlib + Cartopy projection"),
    ("plotly", "backend_plotly.py", "plotly + kaleido static PNG (headless Chromium)"),
    ("webengine", "backend_webengine.py", "QtWebEngine MapLibre, GPU-disabled"),
]


def run_one(name: str, script: str, timeout: int = 150):
    out_png = os.path.join(OUT, f"{name}.png")
    if os.path.exists(out_png):
        os.remove(out_png)
    t0 = time.time()
    try:
        proc = subprocess.run(
            [BACKEND_PY, os.path.join(HERE, script), out_png],
            cwd=HERE,
            env=dict(os.environ),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return ("TIMEOUT", time.time() - t0, None, f"exceeded {timeout}s")
    except Exception as exc:  # noqa: BLE001
        return ("ERROR", time.time() - t0, None, str(exc))
    dt = time.time() - t0
    made = os.path.exists(out_png) and os.path.getsize(out_png) > 1000
    tail = (proc.stderr or proc.stdout or "").strip().replace("\n", " ")[-260:]
    if proc.returncode == 0 and made:
        return ("OK", dt, out_png, "")
    if proc.returncode < 0:
        return ("CRASH", dt, None, f"killed by signal {-proc.returncode}. {tail}")
    if not made:
        return ("NO-OUTPUT", dt, None, tail or f"exit {proc.returncode}, no PNG")
    return (f"EXIT {proc.returncode}", dt, out_png, tail)


def main() -> None:
    print(f"\nBACTolomew map-backend test  ·  outputs -> {OUT}")
    print(f"backend python: {BACKEND_PY}")
    print("=" * 74)
    results = []
    for name, script, desc in BACKENDS:
        if not os.path.exists(os.path.join(HERE, script)):
            print(f"  {name:12} SKIP  (missing {script})")
            continue
        print(f"  {name:12} running …", flush=True)
        status, dt, png, note = run_one(name, script)
        results.append((name, status, dt, png, note, desc))

    print("\n" + "=" * 74)
    print("RESULTS  (✅ works here · 💥 crashed · ⚠️ error)")
    print("=" * 74)
    for name, status, dt, png, note, desc in results:
        mark = "✅" if status == "OK" else ("💥" if status == "CRASH" else "⚠️ ")
        print(f"{mark} {name:12} {status:11} {dt:5.1f}s   {desc}")
        if note:
            print(f"      └─ {note}")
    ok = [r[0] for r in results if r[1] == "OK"]
    print("\nRendered successfully on this machine:", ", ".join(ok) if ok else "(none)")
    print(f"Compare the look by opening the PNGs in:  {OUT}\n")


if __name__ == "__main__":
    main()
