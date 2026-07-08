from __future__ import annotations

from dataclasses import dataclass
import traceback
from typing import Any, Callable, Dict, List, Optional

import polars as pl
from PyQt6 import QtCore


@dataclass(frozen=True)
class RefreshSnapshot:
    token: int
    sample_query: Optional[str]
    meta_filters: Dict[str, List[str]]
    sample_subset_samples: List[str]
    ann_filters: Dict[str, Any]
    hide_empty_contigs: bool
    view_mode: str


@dataclass
class RefreshPayload:
    token: int
    ann_filters: Dict[str, Any]
    allowed_samples: List[str]
    panel_allowed_samples: List[str]
    filtered_ann_pl: pl.DataFrame
    meta_for_tree_pl: pl.DataFrame
    contig_records: Dict[str, List[Dict[str, Any]]]
    filtered_counts: Dict[str, int]
    tree_allowed_samples: List[str]
    stats_text: str


class RefreshWorkerSignals(QtCore.QObject):
    finished = QtCore.pyqtSignal(object)
    failed = QtCore.pyqtSignal(int, str)


class RefreshWorker(QtCore.QRunnable):
    def __init__(self, token: int, compute: Callable[[], RefreshPayload]) -> None:
        super().__init__()
        self._token = int(token)
        self._compute = compute
        self.signals = RefreshWorkerSignals()
        self.setAutoDelete(True)

    def run(self) -> None:
        try:
            payload = self._compute()
        except Exception:
            self.signals.failed.emit(self._token, traceback.format_exc())
            return
        self.signals.finished.emit(payload)
