"""Wall-clock budgets for the timing tests (window plan 8.12).

A busy machine can make one measurement slow although the code is fine, so a
timing test measures a median, and ``within_budget`` measures once more before
it fails. ``BACTOLOMEW_TEST_TIMING_FACTOR`` (default 1) scales every budget,
for a slow or heavily loaded test machine.
"""
from __future__ import annotations

import os
from typing import Callable, List, Tuple

_FACTOR_ENV = "BACTOLOMEW_TEST_TIMING_FACTOR"


def budget_ms(base_ms: float) -> float:
    try:
        factor = float(os.environ.get(_FACTOR_ENV, "1") or 1)
    except ValueError:
        factor = 1.0
    return base_ms * max(1.0, factor)


def within_budget(measure_ms: Callable[[], float], base_ms: float, tries: int = 2) -> Tuple[bool, List[float]]:
    """Run ``measure_ms`` up to ``tries`` times; True once a value fits the budget."""
    limit = budget_ms(base_ms)
    seen: List[float] = []
    for _ in range(max(1, tries)):
        value = measure_ms()
        seen.append(value)
        if value < limit:
            return True, seen
    return False, seen
