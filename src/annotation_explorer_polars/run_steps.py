"""Run progress by pipeline steps, read from run_state.json (no Qt).

The Runner lists the steps a run will consider in ``planned_steps`` (name and
label, in run order) and keeps one entry per started step in ``steps`` (status,
start and end time, the reason of a reused or skipped step, and a
``progress`` {done, total, unit} for the steps that go genome by genome, or
{label, done, total, unit "%"} for work inside a step, such as the cgMLST scheme
install of the chewBBACA step). The
start window's run cards and the Run progress window both show this summary;
neither counts output files.

Runs written before ``planned_steps`` existed list the steps present in
``steps``, in the pipeline's usual order.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Optional

from pipeline.RunSteps import PIPELINE_STEPS, step_label

# For run_state files without planned_steps: the pipeline's step order.
_KNOWN_ORDER = {name: index for index, (name, _label) in enumerate(PIPELINE_STEPS)}

# Row states. Finished ones count toward the bar.
WAITING = "waiting"
RUNNING = "running"
DONE = "done"
REUSED = "reused"
SKIPPED = "skipped"
FAILED = "failed"
PARTIAL = "partial"
HALTED = "halted"  # started, but the run stopped or ended before it finished
FINISHED_STATES = frozenset({DONE, REUSED, SKIPPED})


@dataclass(frozen=True)
class StepRow:
    name: str
    label: str
    state: str
    text: str  # e.g. "running · 14 min · 37/153 genomes"
    detail: str  # tooltip: the error of a failed step, the reason of a skipped one
    seconds: Optional[float]
    done: Optional[int] = None
    total: Optional[int] = None
    unit: str = ""
    note: str = ""  # work inside the step: "Installing cgMLST scheme · A. baumannii"

    @property
    def count_text(self) -> str:
        """'37/153 genomes', or '38%' for work inside the step."""
        if self.done is None or not self.total:
            return ""
        if self.unit == "%":
            return f"{self.done}%"
        return f"{self.done}/{self.total} {self.unit}".rstrip()

    @property
    def progress_text(self) -> str:
        """The count, after the text of work inside the step ('Installing ... · 38%')."""
        return " · ".join(b for b in (self.note, self.count_text) if b)

    @property
    def fraction(self) -> float:
        # Work inside a step (an install before the tool runs) does not fill the step's share.
        if self.note or self.done is None or not self.total:
            return 0.0
        return max(0.0, min(1.0, self.done / self.total))


@dataclass(frozen=True)
class RunSteps:
    rows: tuple[StepRow, ...]
    status: str
    elapsed_seconds: Optional[float]
    current_index: Optional[int]  # the running step, else the first unfinished one

    @property
    def total(self) -> int:
        return len(self.rows)

    @property
    def finished(self) -> int:
        return sum(1 for row in self.rows if row.state in FINISHED_STATES)

    @property
    def current(self) -> Optional[StepRow]:
        if self.current_index is None:
            return None
        return self.rows[self.current_index]

    @property
    def fraction(self) -> float:
        """Finished steps over planned steps, plus the running step's share."""
        if not self.rows:
            return 0.0
        value = float(self.finished)
        current = self.current
        if current is not None and current.state == RUNNING:
            value += current.fraction
        return max(0.0, min(1.0, value / len(self.rows)))

    def step_line(self) -> str:
        """The run card's line while running: 'Step 9 of 15 · Gubbins · 14 min'."""
        current = self.current
        if current is None:
            return ""
        bits = [f"Step {self.current_index + 1} of {self.total}", current.label]
        if current.state == RUNNING and current.seconds is not None:
            bits.append(format_duration(current.seconds))
        if current.progress_text:
            bits.append(current.progress_text)
        return " · ".join(bits)

    def header_text(self) -> str:
        """The Run progress window's header: 'Step 9 of 15 · 1 h 12 min'."""
        if not self.rows:
            return format_duration(self.elapsed_seconds) if self.elapsed_seconds is not None else ""
        if self.status == "running" and self.current is not None:
            bits = [f"Step {self.current_index + 1} of {self.total}"]
        else:
            bits = [f"{self.finished} of {self.total} steps"]
        if self.elapsed_seconds is not None:
            bits.append(format_duration(self.elapsed_seconds))
        return " · ".join(bits)


def format_duration(seconds: Optional[float]) -> str:
    """'8 s', '14 min', '1 h 12 min'."""
    if seconds is None:
        return ""
    total = max(0, int(round(seconds)))
    if total < 60:
        return f"{total} s"
    minutes = total // 60
    if minutes < 60:
        return f"{minutes} min"
    hours, minutes = divmod(minutes, 60)
    return f"{hours} h {minutes} min" if minutes else f"{hours} h"


def _parse_time(value: object) -> Optional[datetime]:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        stamp = datetime.fromisoformat(text)
    except ValueError:
        return None
    if stamp.tzinfo is None:
        stamp = stamp.astimezone()
    return stamp


def _seconds_between(start: object, end: object) -> Optional[float]:
    first = _parse_time(start)
    second = _parse_time(end)
    if first is None or second is None:
        return None
    return max(0.0, (second - first).total_seconds())


def _as_int(value: object) -> Optional[int]:
    try:
        return int(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def _planned(state: dict, status: str) -> list[tuple[str, str, dict]]:
    """(name, label, planned entry) in run order.

    Without planned_steps (runs from before it), the steps present in ``steps`` in
    the pipeline's order; a run that did not succeed also lists the known steps
    after the last one it reached, so its bar does not read as complete (a known
    step before that one is one its pipeline version did not have).
    """
    steps = state.get("steps") if isinstance(state.get("steps"), dict) else {}
    raw = state.get("planned_steps")
    planned: list[tuple[str, str, dict]] = []
    seen: set[str] = set()
    if isinstance(raw, list) and raw:
        for item in raw:
            entry = item if isinstance(item, dict) else {"name": item}
            name = str(entry.get("name") or "").strip()
            if not name or name in seen:
                continue
            seen.add(name)
            label = str(entry.get("label") or "").strip() or step_label(name)
            planned.append((name, label, entry))
        # A step the Runner started without planning it still shows, after the plan.
        extra = [name for name in steps if name not in seen]
    else:
        names = set(steps)
        if status != "success":
            reached = max((_KNOWN_ORDER[name] for name in names if name in _KNOWN_ORDER), default=-1)
            names.update(name for name, index in _KNOWN_ORDER.items() if index > reached)
        extra = sorted(names, key=lambda name: (_KNOWN_ORDER.get(name, len(_KNOWN_ORDER)), name))
    for name in extra:
        planned.append((str(name), step_label(str(name)), {}))
    return planned


def _row(name: str, label: str, planned: dict, step: Optional[dict], *, run_status: str, now: datetime) -> StepRow:
    progress = step.get("progress") if isinstance(step, dict) else None
    done = total = None
    unit = note = ""
    if isinstance(progress, dict):
        done = _as_int(progress.get("done"))
        total = _as_int(progress.get("total"))
        unit = str(progress.get("unit") or "").strip()
        note = str(progress.get("label") or "").strip()
    count = StepRow(name, label, "", "", "", None, done, total, unit, note).count_text

    if not isinstance(step, dict):
        skip = planned.get("skip") or planned.get("skipped")
        if skip:
            reason = skip if isinstance(skip, str) else str(planned.get("reason") or "")
            return StepRow(name, label, SKIPPED, "skipped", _skip_reason(reason), None)
        text = "waiting" if run_status == "running" else "not run"
        return StepRow(name, label, WAITING, text, "", None)

    status = str(step.get("status") or "").strip().lower()
    reason = str(step.get("reason") or "").strip()
    start = step.get("start_time")
    if status == "running":
        if run_status == "running":
            seconds = _seconds_between(start, now.isoformat())
            bits = ["running", format_duration(seconds) if seconds is not None else "", count]
            return StepRow(name, label, RUNNING, " · ".join(b for b in bits if b), "", seconds, done, total, unit,
                           note)
        halted = "stopped" if run_status == "stopped" else "interrupted"
        text = " · ".join(b for b in (halted, "" if note else count) if b)
        return StepRow(name, label, HALTED, text, "", None, done, total, unit, note)

    seconds = _seconds_between(start, step.get("end_time"))
    if status == "success":
        if reason == "reused_existing_output":
            # Checking reused outputs is quick; one that took a while (the core tree's
            # refit) shows how long.
            bits = ["reused", format_duration(seconds) if seconds is not None and seconds >= 60 else ""]
            return StepRow(name, label, REUSED, " · ".join(b for b in bits if b), "", seconds, done, total, unit)
        bits = ["done", format_duration(seconds) if seconds is not None else ""]
        if total and not note:
            bits.append(f"{total} {unit}".rstrip())
        return StepRow(name, label, DONE, " · ".join(b for b in bits if b), "", seconds, done, total, unit)
    if status == "skipped":
        return StepRow(name, label, SKIPPED, "skipped", _skip_reason(reason), seconds)
    errors = step.get("validation_errors")
    error_lines = [str(e).strip() for e in errors if str(e).strip()] if isinstance(errors, list) else []
    error = str(step.get("error") or "").strip()
    detail = "\n".join([error] + error_lines if error else error_lines) or reason
    if note:  # a failed install: its text and share say nothing about the step's outputs
        count = ""
    if status == "partial":
        text = " · ".join(b for b in ("incomplete", count) if b)
        return StepRow(name, label, PARTIAL, text, detail, seconds, done, total, unit)
    text = " · ".join(b for b in ("failed", count) if b)
    return StepRow(name, label, FAILED, text, detail, seconds, done, total, unit)


def _skip_reason(reason: str) -> str:
    text = str(reason or "").strip()
    if text.lower().startswith("skipped:"):
        text = text[len("skipped:"):].strip()
    return text


def summarize_run_steps(
    state: Optional[dict],
    status: str,
    *,
    now: Optional[datetime] = None,
) -> RunSteps:
    """The step rows of a run from its run_state.json payload.

    ``status`` is the run's status as ``run_expectations.load_run_status`` reports it
    ("running" only while the run's lock is held).
    """
    state = state if isinstance(state, dict) else {}
    status = str(status or "").strip().lower()
    now = now or datetime.now(timezone.utc)
    steps = state.get("steps") if isinstance(state.get("steps"), dict) else {}
    rows = tuple(
        _row(name, label, entry, steps.get(name), run_status=status, now=now)
        for name, label, entry in _planned(state, status)
    )

    current_index: Optional[int] = None
    current_name = str(state.get("current_step") or "").strip()
    for index, row in enumerate(rows):
        if row.state == RUNNING or (current_name and row.name == current_name and row.state not in FINISHED_STATES):
            current_index = index
            break
    if current_index is None and status == "running":
        current_index = next((i for i, row in enumerate(rows) if row.state not in FINISHED_STATES), None)

    if status == "running":
        elapsed = _seconds_between(state.get("start_time"), now.isoformat())
    else:
        elapsed = _seconds_between(state.get("start_time"), state.get("end_time") or state.get("updated_at"))
    return RunSteps(rows=rows, status=status, elapsed_seconds=elapsed, current_index=current_index)


__all__ = [
    "FINISHED_STATES",
    "RunSteps",
    "StepRow",
    "format_duration",
    "summarize_run_steps",
]
