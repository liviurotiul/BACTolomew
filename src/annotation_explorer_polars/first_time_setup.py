from __future__ import annotations

import threading
from dataclasses import replace
from pathlib import Path
from typing import Optional

from PyQt6 import QtCore, QtGui, QtWidgets

from pipeline.Setup import (
    STEP_DONE,
    STEP_FAILED,
    SetupJob,
    SetupStep,
    setup_fraction,
)

from .run_progress import _set_state


# --------------------------------------------------------------------------- #
# The window
# --------------------------------------------------------------------------- #
class _ElidedLabel(QtWidgets.QLabel):
    """One line that ends in "..." when it does not fit; the full text is its tooltip."""

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self._full = ""
        self.setSizePolicy(QtWidgets.QSizePolicy.Policy.Ignored, QtWidgets.QSizePolicy.Policy.Preferred)

    def set_full_text(self, text: str) -> None:
        self._full = text
        self._elide()

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:  # noqa: N802 - Qt override
        super().resizeEvent(event)
        self._elide()

    def _elide(self) -> None:
        shown = self.fontMetrics().elidedText(self._full, QtCore.Qt.TextElideMode.ElideRight, max(0, self.width()))
        if shown != self.text():
            self.setText(shown)
        self.setToolTip(self._full if shown != self._full else "")


class _SetupRow(QtWidgets.QWidget):
    """One step: state marker, name (with the failure, or a hint, below it), size and state."""

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setObjectName("RunStepRow")
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_StyledBackground, True)
        row = QtWidgets.QHBoxLayout(self)
        row.setContentsMargins(12, 7, 12, 7)
        row.setSpacing(10)
        self.marker = QtWidgets.QLabel()
        self.marker.setProperty("role", "stepMarker")
        row.addWidget(self.marker, 0, QtCore.Qt.AlignmentFlag.AlignVCenter)
        names = QtWidgets.QVBoxLayout()
        names.setSpacing(2)
        self.name = QtWidgets.QLabel()
        self.name.setProperty("role", "stepName")
        names.addWidget(self.name)
        self.note = _ElidedLabel()
        self.note.setProperty("role", "stepNote")
        self.note.setVisible(False)
        names.addWidget(self.note)
        self.hint = _ElidedLabel()
        self.hint.setProperty("role", "stepDetail")
        self.hint.setVisible(False)
        names.addWidget(self.hint)
        row.addLayout(names, 1)
        self.size = QtWidgets.QLabel()
        self.size.setProperty("role", "stepState")
        self.size.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
        row.addWidget(self.size, 0, QtCore.Qt.AlignmentFlag.AlignVCenter)
        self.state = QtWidgets.QLabel()
        self.state.setProperty("role", "stepState")
        self.state.setMinimumWidth(78)
        self.state.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
        row.addWidget(self.state, 0, QtCore.Qt.AlignmentFlag.AlignVCenter)

    def set_step(self, step: SetupStep) -> None:
        self.name.setText(step.label)
        self.size.setText(step.size)
        self.state.setText(step.text)
        self.note.set_full_text(step.note)
        self.note.setVisible(bool(step.note))
        self.hint.set_full_text(step.hint)
        self.hint.setVisible(bool(step.hint) and not step.note and step.state != STEP_DONE)
        for widget in (self.marker, self.name, self.state):
            _set_state(widget, step.state)
        self.setToolTip(step.detail)


class _SetupBridge(QtCore.QObject):
    step = QtCore.pyqtSignal(object)
    line = QtCore.pyqtSignal(str)


class _SetupWorker(QtCore.QThread):
    outcome = QtCore.pyqtSignal(str)

    def __init__(self, job: SetupJob, parent: Optional[QtCore.QObject] = None) -> None:
        super().__init__(parent)
        self.job = job

    def run(self) -> None:  # worker thread
        self.job.thread_id = threading.get_ident()
        try:
            result = self.job.run()
        except Exception as exc:  # noqa: BLE001 - a bug in a step: shown as a failure
            self.job._log_fn(f"ERROR: {exc}")
            result = "failed"
        finally:
            self.job.thread_id = None
        self.outcome.emit(result)


class FirstTimeSetupDialog(QtWidgets.QDialog):
    """First-time setup in one window: it starts by itself, shows one row per step and
    ends with Start (or Retry after a failure). The work runs on a worker thread."""

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None, *, dev_mode: bool = False,
                 root: Optional[Path] = None) -> None:
        super().__init__(parent)
        self.setObjectName("FirstTimeSetupDialog")
        self.setWindowTitle("First-Time Setup")
        self.resize(540, 0)
        self._bridge = _SetupBridge(self)
        self._bridge.step.connect(self._on_step)
        self._bridge.line.connect(self._append_log)
        self.job = SetupJob(root, dev_mode=dev_mode, on_step=self._bridge.step.emit, log=self._bridge.line.emit)
        self._worker: Optional[_SetupWorker] = None
        self._stopping = False
        self._outcome = ""  # how the last run ended: "done", "failed", "cancelled"
        self._rows: dict[str, _SetupRow] = {}
        self._view = {step.key: replace(step) for step in self.job.steps}  # what the rows show
        self._build_ui()
        for step in list(self._view.values()):
            self._on_step(step)
        QtCore.QTimer.singleShot(0, self._start)

    # ---- UI ----
    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(20, 18, 20, 14)
        layout.setSpacing(10)
        layout.setSizeConstraint(QtWidgets.QLayout.SizeConstraint.SetMinimumSize)

        header = QtWidgets.QHBoxLayout()
        title = QtWidgets.QLabel("First-Time Setup")
        title.setProperty("role", "dialogTitle")
        header.addWidget(title, 1)
        self.summary = QtWidgets.QLabel("")
        self.summary.setProperty("role", "runCount")
        self.summary.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
        header.addWidget(self.summary, 0)
        layout.addLayout(header)
        self.status = QtWidgets.QLabel("Getting BACTolomew ready")
        self.status.setProperty("role", "muted")
        layout.addWidget(self.status)

        self.progress = QtWidgets.QProgressBar()
        self.progress.setProperty("role", "stepBar")
        self.progress.setRange(0, 1000)
        self.progress.setTextVisible(False)
        layout.addSpacing(2)
        layout.addWidget(self.progress)
        layout.addSpacing(4)

        card = QtWidgets.QFrame()
        card.setObjectName("RunStepsCard")
        rows = QtWidgets.QVBoxLayout(card)
        rows.setContentsMargins(1, 4, 1, 4)
        rows.setSpacing(0)
        for step in self._view.values():
            row = _SetupRow()
            rows.addWidget(row)
            self._rows[step.key] = row
        layout.addWidget(card)

        self.details_toggle = QtWidgets.QToolButton()
        self.details_toggle.setObjectName("SetupDetailsToggle")
        self.details_toggle.setText("Details")
        self.details_toggle.setCheckable(True)
        self.details_toggle.setArrowType(QtCore.Qt.ArrowType.RightArrow)
        self.details_toggle.setToolButtonStyle(QtCore.Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.details_toggle.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.details_toggle.toggled.connect(self._toggle_details)
        layout.addWidget(self.details_toggle, 0, QtCore.Qt.AlignmentFlag.AlignLeft)
        self.details = QtWidgets.QPlainTextEdit()
        self.details.setObjectName("SetupLog")
        self.details.setReadOnly(True)
        self.details.setMaximumBlockCount(3000)
        self.details.setLineWrapMode(QtWidgets.QPlainTextEdit.LineWrapMode.NoWrap)
        mono = QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.SystemFont.FixedFont)
        mono.setPointSize(max(8, mono.pointSize() - 1))
        self.details.setFont(mono)
        self.details.setMinimumHeight(170)
        self.details.setVisible(False)
        layout.addWidget(self.details, 1)
        layout.addStretch(0)  # a taller window grows the log, or this space above the buttons

        footer = QtWidgets.QFrame()
        footer.setObjectName("ProgressFooter")
        buttons = QtWidgets.QHBoxLayout(footer)
        buttons.setContentsMargins(0, 10, 0, 0)
        buttons.setSpacing(8)
        buttons.addStretch(1)
        self.cancel_btn = QtWidgets.QPushButton("Cancel")
        self.cancel_btn.setProperty("role", "flat")
        self.cancel_btn.clicked.connect(self.reject)
        buttons.addWidget(self.cancel_btn)
        self.retry_btn = QtWidgets.QPushButton("Retry")
        self.retry_btn.setProperty("role", "primary")
        self.retry_btn.clicked.connect(self._start)
        self.retry_btn.setVisible(False)
        buttons.addWidget(self.retry_btn)
        self.start_btn = QtWidgets.QPushButton("Start")
        self.start_btn.setProperty("role", "primary")
        self.start_btn.setEnabled(False)
        self.start_btn.clicked.connect(self.accept)
        buttons.addWidget(self.start_btn)
        layout.addWidget(footer)

    def _toggle_details(self, shown: bool) -> None:
        self.details_toggle.setArrowType(QtCore.Qt.ArrowType.DownArrow if shown else QtCore.Qt.ArrowType.RightArrow)
        self.details.setVisible(shown)
        self._fit()

    def _fit(self) -> None:
        """Without the log, the window is as tall as its rows (a failure line comes and goes)."""
        if not self.details.isVisible():
            QtCore.QTimer.singleShot(0, self._fit_now)

    def _fit_now(self) -> None:
        for _ in range(4):  # row, card, dialog: each passes its new size up one level
            QtWidgets.QApplication.sendPostedEvents(None, QtCore.QEvent.Type.LayoutRequest.value)
        self.resize(self.width(), self._needed_height())

    def _needed_height(self) -> int:
        """The height the content needs at the current width (wrapped labels included)."""
        need = self.minimumSizeHint().height()
        layout = self.layout()
        if layout is not None and layout.hasHeightForWidth():
            need = max(need, layout.totalHeightForWidth(self.width()))
        return need

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:  # noqa: N802 - Qt override
        super().resizeEvent(event)
        if event.size().width() != event.oldSize().width():
            QtCore.QTimer.singleShot(0, self._grow_to_content)

    def _grow_to_content(self) -> None:
        """A narrower window can wrap a line: grow (never shrink) to keep every row whole."""
        need = self._needed_height()
        if self.height() < need:
            self.resize(self.width(), need)

    def _append_log(self, line: str) -> None:
        self.details.appendPlainText(line)

    def _on_step(self, step: SetupStep) -> None:
        old = self._view[step.key]
        noted = (bool(old.note), old.state == STEP_DONE) != (bool(step.note), step.state == STEP_DONE)
        self._view[step.key] = step
        self._rows[step.key].set_step(step)
        if noted:
            self._fit()
        steps = list(self._view.values())
        done = sum(1 for s in steps if s.state == STEP_DONE)
        self.summary.setText(f"{done} of {len(steps)} steps")
        self.progress.setValue(int(round(setup_fraction(steps) * 1000)))
        failed = any(s.state == STEP_FAILED for s in steps)
        _set_state(self.progress, "failed" if failed else ("success" if done == len(steps) else "running"))

    # ---- run ----
    def _start(self) -> None:
        if self._worker is not None and self._worker.isRunning():
            return
        self._stopping = False
        self._outcome = ""
        self.job.rearm()
        self.status.setText("Getting BACTolomew ready")
        self.retry_btn.setVisible(False)
        self.start_btn.setVisible(True)
        self.start_btn.setEnabled(False)
        self.cancel_btn.setVisible(True)
        self.cancel_btn.setEnabled(True)
        self._worker = _SetupWorker(self.job, self)
        self._worker.outcome.connect(self._finished)
        self._worker.start()

    def is_running(self) -> bool:
        return self._worker is not None and self._worker.isRunning()

    def _finished(self, outcome: str) -> None:
        if self._worker is not None:
            self._worker.wait()
        self._outcome = outcome
        if outcome == "cancelled" or self._stopping:
            QtWidgets.QDialog.reject(self)
            return
        if outcome == "done":
            self.status.setText("BACTolomew is ready")
            self.cancel_btn.setVisible(False)
            self.start_btn.setEnabled(True)
            self.start_btn.setFocus()
            return
        self.status.setText("Setup could not finish")
        self.start_btn.setVisible(False)
        self.retry_btn.setVisible(True)
        self.retry_btn.setFocus()

    # ---- closing ----
    def reject(self) -> None:
        """Cancel, Escape and the window's close button: during the setup this asks once,
        then stops the setup (the download in progress is discarded) and closes."""
        if not self.is_running():
            super().reject()
            return
        if self._stopping:
            return
        box = QtWidgets.QMessageBox(self)
        box.setObjectName("SetupStopQuestion")
        box.setIcon(QtWidgets.QMessageBox.Icon.NoIcon)
        box.setWindowTitle("First-Time Setup")
        box.setText("Stop the setup?")
        box.setInformativeText("The download in progress is discarded. Setup continues the next time "
                               "BACTolomew starts.")
        stop = box.addButton("Stop setup", QtWidgets.QMessageBox.ButtonRole.DestructiveRole)
        keep = box.addButton("Keep going", QtWidgets.QMessageBox.ButtonRole.RejectRole)
        for button, role in ((stop, "flat"), (keep, "primary")):
            button.setProperty("role", role)
            button.style().unpolish(button)
            button.style().polish(button)
        box.setDefaultButton(keep)
        box.exec()
        if box.clickedButton() is not stop:
            return
        if not self.is_running():  # it ended while the question was open
            if self._outcome != "done":  # finished setups stay open for Start
                super().reject()
            return
        self._stopping = True
        self.status.setText("Stopping...")
        self.cancel_btn.setEnabled(False)
        self.job.cancel()
