from __future__ import annotations

from pathlib import Path
from typing import Optional

from PyQt6 import QtCore, QtWidgets

from .run_expectations import (
    describe_run_status,
    load_run_status,
    pending_outputs_text,
)
from .run_steps import RUNNING, StepRow, summarize_run_steps

TERMINAL_STATUSES = frozenset({"success", "failed", "stopped", "interrupted", "partial"})


def _compact_status_text(status_info: dict[str, object], *, pending_text: Optional[str], full_status: str) -> str:
    status = str(status_info.get("status") or "").strip().lower()
    current_step = str(status_info.get("current_step") or "").strip()
    if status == "running":
        return f"Running: {current_step}" if current_step else "Running"
    if status == "success":
        return "Completed."
    if status == "stopped":
        return "Stopped by user. Click Resume to continue."
    if status == "interrupted":
        return "Interrupted. Click Resume to continue."
    if status == "failed":
        return "Failed. Click Resume to retry."
    if status == "partial":
        return pending_text or "Partial run. Click Resume to continue."
    return full_status or pending_text or "Status unknown. Click Resume to continue."


def _compact_failure_text(error: str) -> str:
    text = str(error or "").strip()
    if not text:
        return "Failed. Click Resume to retry."
    first_line = text.splitlines()[0].strip()
    if len(first_line) > 140:
        first_line = first_line[:137].rstrip() + "..."
    return first_line or "Failed. Click Resume to retry."


def status_bar_state(status: str) -> str:
    """The step bar's color state for a run status (theme: QProgressBar[role="stepBar"])."""
    status = str(status or "").strip().lower()
    if status == "success":
        return "success"
    if status == "failed":
        return "failed"
    if status in {"stopped", "interrupted", "partial"}:
        return "stopped"
    return "running"


def pending_text_for(status: str, output_dir: Path) -> Optional[str]:
    """The partial run's explanation from its expected outputs; only partial and
    unknown runs show it, so no other run checks files."""
    if str(status or "").strip().lower() in {"partial", "unknown"}:
        return pending_outputs_text(output_dir)
    return None


def _set_state(widget: QtWidgets.QWidget, state: str) -> None:
    if widget.property("state") == state:
        return
    widget.setProperty("state", state)
    widget.style().unpolish(widget)
    widget.style().polish(widget)


class _StepRowWidget(QtWidgets.QWidget):
    """One planned step: a state marker, its label, and its state text."""

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setObjectName("RunStepRow")
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_StyledBackground, True)
        row = QtWidgets.QHBoxLayout(self)
        row.setContentsMargins(12, 7, 12, 7)
        row.setSpacing(10)
        self.marker = QtWidgets.QLabel()
        self.marker.setProperty("role", "stepMarker")
        self.name = QtWidgets.QLabel()
        self.name.setProperty("role", "stepName")
        # Work inside the step, under its name: "Installing cgMLST scheme · A. baumannii".
        self.note = QtWidgets.QLabel()
        self.note.setProperty("role", "stepDetail")
        self.note.setVisible(False)
        names = QtWidgets.QVBoxLayout()
        names.setContentsMargins(0, 0, 0, 0)
        names.setSpacing(1)
        names.addWidget(self.name)
        names.addWidget(self.note)
        self.state = QtWidgets.QLabel()
        self.state.setProperty("role", "stepState")
        self.state.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
        row.addWidget(self.marker, 0, QtCore.Qt.AlignmentFlag.AlignVCenter)
        row.addLayout(names, 1)
        row.addWidget(self.state, 0)

    def set_row(self, step: StepRow) -> None:
        # The state first: it changes the font, and the texts are then measured with it.
        for widget in (self.marker, self.name, self.state):
            _set_state(widget, step.state)
        self.name.setText(step.label)
        self.note.setText(step.note)
        self.note.setVisible(bool(step.note) and step.state == RUNNING)
        self.state.setText(step.text)
        self.setToolTip(step.detail)


class RunProgressDialog(QtWidgets.QDialog):
    closed = QtCore.pyqtSignal()

    def __init__(self, output_dir: Path, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.output_dir = Path(output_dir)
        self._step_rows: list[_StepRowWidget] = []
        self._step_keys: tuple[tuple[str, str], ...] = ()
        # The runner's error from mark_failed; run_state.json may never have recorded it.
        self._launch_error: Optional[str] = None
        self.setObjectName("RunProgressDialog")
        self.setWindowTitle("Run progress")
        self.resize(560, 680)
        self._build_ui()

        self.timer = QtCore.QTimer(self)
        self.timer.setInterval(2000)
        self.timer.timeout.connect(self._refresh)
        self.timer.start()
        self._refresh()

    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(20, 18, 20, 14)
        layout.setSpacing(10)

        header = QtWidgets.QHBoxLayout()
        header.setSpacing(12)
        self.title = QtWidgets.QLabel("Run in progress")
        self.title.setProperty("role", "dialogTitle")
        header.addWidget(self.title, 1)
        self.summary = QtWidgets.QLabel("")
        self.summary.setProperty("role", "runCount")
        self.summary.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
        header.addWidget(self.summary, 0)
        layout.addLayout(header)

        self.status = QtWidgets.QLabel("")
        self.status.setWordWrap(True)
        self.status.setProperty("role", "muted")
        self.status.setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignTop)
        self.status.setMaximumHeight(self.fontMetrics().lineSpacing() * 3)
        layout.addWidget(self.status)

        self.details = QtWidgets.QPlainTextEdit()
        self.details.setObjectName("RunErrorDetails")
        self.details.setReadOnly(True)
        self.details.setVisible(False)
        self.details.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed)
        self.details.setLineWrapMode(QtWidgets.QPlainTextEdit.LineWrapMode.WidgetWidth)
        layout.addWidget(self.details)

        self.progress = QtWidgets.QProgressBar()
        self.progress.setProperty("role", "stepBar")
        self.progress.setRange(0, 1000)
        self.progress.setValue(0)
        self.progress.setTextVisible(False)
        layout.addSpacing(2)
        layout.addWidget(self.progress)
        layout.addSpacing(4)

        card = QtWidgets.QFrame()
        card.setObjectName("RunStepsCard")
        card_layout = QtWidgets.QVBoxLayout(card)
        card_layout.setContentsMargins(1, 4, 1, 4)
        card_layout.setSpacing(0)
        scroll = QtWidgets.QScrollArea()
        scroll.setObjectName("RunStepsScroll")
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        content = QtWidgets.QWidget()
        content.setObjectName("RunStepsContent")
        self._rows_layout = QtWidgets.QVBoxLayout(content)
        self._rows_layout.setContentsMargins(0, 0, 0, 0)
        self._rows_layout.setSpacing(0)
        self.empty = QtWidgets.QLabel("No steps recorded yet.")
        self.empty.setProperty("role", "muted")
        self.empty.setContentsMargins(12, 8, 12, 8)
        self._rows_layout.addWidget(self.empty)
        self._rows_layout.addStretch(1)
        scroll.setWidget(content)
        card_layout.addWidget(scroll)
        layout.addWidget(card, 1)

        footer = QtWidgets.QFrame()
        footer.setObjectName("ProgressFooter")
        btn_row = QtWidgets.QHBoxLayout(footer)
        btn_row.setContentsMargins(0, 10, 0, 0)
        btn_row.addStretch(1)
        self.close_btn = QtWidgets.QPushButton("Close")
        self.close_btn.setEnabled(True)
        self.close_btn.clicked.connect(self.close)
        btn_row.addWidget(self.close_btn)
        layout.addWidget(footer)

    def set_output_dir(self, output_dir: Path) -> None:
        self.output_dir = Path(output_dir)
        self._launch_error = None
        self.title.setText("Run in progress")
        self.status.setText("")
        self.status.setToolTip("")
        self.details.clear()
        self.details.setVisible(False)
        self.progress.setValue(0)
        self._set_steps(())
        if not self.timer.isActive():
            self.timer.start()
        self._refresh()

    def _show_details(self, text: str) -> None:
        """The error box, as tall as its text up to six lines (then it scrolls)."""
        self.details.setPlainText(text)
        self.details.setVisible(bool(text))
        if not text:
            return
        metrics = self.details.fontMetrics()
        width = max(120, self.details.viewport().width() or self.width() - 64)
        lines = sum(
            max(1, -(-metrics.horizontalAdvance(line) // width)) for line in text.splitlines() or [""]
        )
        self.details.setFixedHeight(metrics.lineSpacing() * min(6, max(1, lines)) + 30)

    def mark_failed(self, error: str) -> None:
        self._launch_error = str(error or "")
        self.timer.stop()
        self._refresh()

    def _show_launch_error(self) -> None:
        error = self._launch_error or ""
        self.title.setText("Run failed")
        compact = _compact_failure_text(error)
        self.status.setText(compact)
        self.status.setToolTip(error)
        detail = error.strip()
        self._show_details(detail if detail != compact else "")

    def mark_finished(self) -> None:
        self._refresh()

    def _refresh(self) -> None:
        status_info = load_run_status(self.output_dir)
        status = str(status_info.get("status") or "")
        state = status_info.get("state")
        steps = summarize_run_steps(state if isinstance(state, dict) else None, status)

        if self._launch_error is not None:
            # The runner's own error stays; run_state.json may still say running.
            self._show_launch_error()
        else:
            pending = pending_text_for(status, self.output_dir)
            full_status = describe_run_status(status_info, pending) or ""
            status_text = _compact_status_text(status_info, pending_text=pending, full_status=full_status)
            self.title.setText(self._title_for_status(status))
            self.status.setText(status_text)
            self.status.setToolTip(full_status)
            self._set_details(full_status, status_text, status)

        self.summary.setText(steps.header_text())
        fraction = steps.fraction
        if not steps.rows and status == "success":
            fraction = 1.0
        self.progress.setValue(int(round(fraction * 1000)))
        _set_state(self.progress, status_bar_state("failed" if self._launch_error is not None else status))
        self._set_steps(steps.rows)

        if status in TERMINAL_STATUSES and self.timer.isActive():
            self.timer.stop()

    def _set_steps(self, rows: tuple[StepRow, ...]) -> None:
        keys = tuple((row.name, row.label) for row in rows)
        if keys != self._step_keys:
            for widget in self._step_rows:
                self._rows_layout.removeWidget(widget)
                widget.deleteLater()
            self._step_rows = []
            for index in range(len(rows)):
                widget = _StepRowWidget()
                self._rows_layout.insertWidget(index, widget)
                self._step_rows.append(widget)
            self._step_keys = keys
        for widget, row in zip(self._step_rows, rows):
            widget.set_row(row)
        self.empty.setVisible(not rows)

    def _title_for_status(self, status: str) -> str:
        status = str(status or "").strip().lower()
        if status == "success":
            return "Run completed"
        if status == "failed":
            return "Run failed"
        if status == "stopped":
            return "Run stopped"
        if status == "interrupted":
            return "Run interrupted"
        if status == "partial":
            return "Run partial"
        return "Run in progress"

    def _set_details(self, full_status: str, compact_status: str, status: str) -> None:
        status = str(status or "").strip().lower()
        detail = str(full_status or "").strip()
        should_show = status in {"failed", "partial"} and bool(detail) and detail != compact_status
        self._show_details(detail if should_show else "")

    def closeEvent(self, event: QtCore.QEvent) -> None:
        if self.timer.isActive():
            self.timer.stop()
        self.closed.emit()
        super().closeEvent(event)
