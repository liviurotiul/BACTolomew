"""The pipeline's steps in run order, with their short labels (no heavy imports).

The Runner writes them into run_state.json as "planned_steps" at the start of each run;
the GUI reads that list, and falls back to PIPELINE_STEPS for run_state files written
before it existed.
"""

from __future__ import annotations

import threading
import time
from typing import Callable, Optional

# (name in run_state.json "steps", short label) in the order the Runner runs them.
PIPELINE_STEPS: tuple[tuple[str, str], ...] = (
    ("assembly_qc", "Assembly QC"),
    ("bakta", "Bakta annotation"),
    ("mlst", "MLST typing"),
    ("amrfinder_mutations", "AMRFinderPlus"),
    ("mobsuite", "MOB-suite plasmids"),
    ("bakta_evidence", "Bakta evidence"),
    ("panaroo", "Panaroo pangenome"),
    ("core_alignment", "Core alignment"),
    ("fastgear_clades", "fastGEAR lineages"),
    ("gubbins", "Gubbins recombination"),
    ("core_tree", "Core tree"),
    ("phylogeny", "Snippy SNP tree"),
    ("chewbbaca", "chewBBACA cgMLST"),
    ("diamond_markers", "Marker databases"),
    ("diamond_reference_loci", "Reference loci search"),
    ("collect_annotations", "Collect annotations"),
)

STEP_LABELS: dict[str, str] = dict(PIPELINE_STEPS)

# Per-genome workflows report progress(done, total) after each genome finishes
# (whether it succeeded or failed); the Runner stores it as the step's "progress".
ProgressCallback = Optional[Callable[[int, int], None]]
# Work inside a step reports status(text, fraction): what it is doing ("Installing cgMLST
# scheme · A. baumannii") and how far (0 to 1, None when unknown); status("", None) ends it.
# The Runner stores it as the step's "progress" with the text as its "label".
StatusCallback = Optional[Callable[[str, Optional[float]], None]]


def step_label(name: str) -> str:
    return STEP_LABELS.get(name) or name.replace("_", " ").capitalize()


def planned_steps() -> list[dict[str, str]]:
    return [{"name": name, "label": label} for name, label in PIPELINE_STEPS]


def report_progress(progress: ProgressCallback, done: int, total: int) -> None:
    """Call a workflow's progress callback; a failing callback never stops the tool."""
    if progress is None:
        return
    try:
        progress(int(done), int(total))
    except Exception:
        pass


def status_progress(text: str, fraction: Optional[float]) -> Optional[dict]:
    """The step "progress" entry of a status(text, fraction) report (None ends it):
    {"label": text, "done": percent, "total": 100, "unit": "%"}, without done/total when
    the fraction is unknown."""
    text = str(text or "").strip()
    if not text:
        return None
    entry: dict = {"label": text, "unit": "%"}
    if fraction is not None:
        entry.update(done=int(max(0.0, min(1.0, float(fraction))) * 100), total=100)
    return entry


def report_status(status: StatusCallback, text: str, fraction: Optional[float] = None) -> None:
    """Call a status callback; a failing callback never stops the work."""
    if status is None:
        return
    try:
        status(text, fraction)
    except Exception:
        pass


class StepProgress:
    """The "progress" of a run's running steps, kept in ``steps`` (the run_state.json
    "steps", guarded by ``lock``) and written with ``write()`` (run_state.json with the run
    still running): at once when a progress first appears, changes its label or is
    removed, else at most every ``interval`` s; a value held back by that limit is written
    by a one-shot timer once the interval has passed, so a slow genome after several quick
    ones never leaves an old count on display. The owner calls ``written()`` after each
    run_state.json write (the throttle counts from it)."""

    def __init__(self, steps: dict, lock, write: Callable[[], None], *, interval: float = 1.0) -> None:
        self._steps = steps
        self._lock = lock
        self._write = write
        self._interval = float(interval)
        self._last_write = 0.0
        self._trailing: Optional[threading.Timer] = None

    def written(self) -> None:
        with self._lock:
            self._last_write = time.monotonic()
            if self._trailing is not None:
                self._trailing.cancel()
                self._trailing = None

    def _running(self, name: str) -> bool:
        step = self._steps.get(name)
        return isinstance(step, dict) and step.get("status") == "running"

    def _write_held_back(self, name: str) -> None:
        with self._lock:
            self._trailing = None
            if self._running(name):
                self._write()

    def store(self, name: str, entry: Optional[dict], *, flush_now: bool = False) -> None:
        """Store (with None remove) step ``name``'s progress while it runs."""
        with self._lock:
            if not self._running(name):
                return
            step = self._steps[name]
            if entry is None:
                if step.pop("progress", None) is not None:
                    self._write()
                return
            previous = step.get("progress")
            step["progress"] = entry
            waited = time.monotonic() - self._last_write
            if previous is None or flush_now or waited >= self._interval:
                self._write()
            elif self._trailing is None:
                self._trailing = threading.Timer(self._interval - waited, self._write_held_back, (name,))
                self._trailing.daemon = True
                self._trailing.start()

    def genome_progress(self, name: str) -> ProgressCallback:
        """progress(done, total) of a per-genome step (written when it starts, when done,
        and throttled in between)."""

        def _progress(done: int, total: int) -> None:
            self.store(name, {"done": int(done), "total": int(total), "unit": "genomes"},
                       flush_now=done >= total)

        return _progress

    def step_status(self, name: str) -> StatusCallback:
        """status(text, fraction) of work inside a step (the cgMLST install of the chewBBACA
        step): stored with the text as its label (a new text is written at once);
        status("", None) removes it."""

        def _status(text: str, fraction: Optional[float]) -> None:
            entry = status_progress(text, fraction)
            with self._lock:
                step = self._steps.get(name)
                previous = step.get("progress") if isinstance(step, dict) else None
                new_label = entry is not None and (previous or {}).get("label") != entry["label"]
                self.store(name, entry, flush_now=new_label)

        return _status


__all__ = [
    "PIPELINE_STEPS",
    "STEP_LABELS",
    "ProgressCallback",
    "StatusCallback",
    "StepProgress",
    "planned_steps",
    "report_progress",
    "report_status",
    "status_progress",
    "step_label",
]
