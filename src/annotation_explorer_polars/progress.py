from __future__ import annotations

import time
from typing import Optional


def log_step(label: str, prefix: str = "[load]") -> float:
    print(f"{prefix} {label}...")
    return time.perf_counter()


def log_done(
    label: str,
    start: Optional[float] = None,
    *,
    prefix: str = "[load]",
    extra: Optional[str] = None,
) -> None:
    elapsed = None
    if start is not None:
        elapsed = time.perf_counter() - start
    if elapsed is not None:
        base = f"{prefix} {label} done in {elapsed:.2f}s"
    else:
        base = f"{prefix} {label} done"
    if extra:
        base = f"{base} ({extra})"
    print(base)


class ProgressLogger:
    def __init__(self, total: Optional[int] = None, prefix: str = "[progress]") -> None:
        self.total = total
        self.prefix = prefix
        self.current = 0

    def step(self, label: str) -> float:
        self.current += 1
        if self.total:
            print(f"{self.prefix} {self.current}/{self.total} {label}...")
        else:
            print(f"{self.prefix} step {self.current}: {label}...")
        return time.perf_counter()

    def done(self, label: str, start: Optional[float] = None, extra: Optional[str] = None) -> None:
        elapsed = None
        if start is not None:
            elapsed = time.perf_counter() - start
        if self.total:
            base = f"{self.prefix} {self.current}/{self.total} {label}"
        else:
            base = f"{self.prefix} step {self.current}: {label}"
        if elapsed is not None:
            base = f"{base} done in {elapsed:.2f}s"
        else:
            base = f"{base} done"
        if extra:
            base = f"{base} ({extra})"
        print(base)
