from __future__ import annotations

import traceback
from typing import Any, Callable

from PyQt6 import QtCore


class RefreshWorkerSignals(QtCore.QObject):
    finished = QtCore.pyqtSignal(object)
    failed = QtCore.pyqtSignal(int, str)


class RefreshWorker(QtCore.QRunnable):
    """Runs one compute on a pool thread and reports it with its token.

    The filter query and its result types (``FilterState``, ``AppliedResult``)
    live in model/filter_state.py, so the model package loads no Qt module.
    """

    def __init__(self, token: int, compute: Callable[[], Any]) -> None:
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
