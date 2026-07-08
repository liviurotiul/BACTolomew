from __future__ import annotations

from pathlib import Path
from typing import Optional

from PyQt6 import QtCore, QtWidgets

from .run_expectations import (
    classify_expected_path,
    describe_run_status,
    expected_output_exists,
    load_run_status,
    resolve_expected_files,
    summarize_pending_outputs,
)


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


class RunProgressDialog(QtWidgets.QDialog):
    closed = QtCore.pyqtSignal()

    def __init__(self, output_dir: Path, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.output_dir = Path(output_dir)
        self.expected_files: list[Path] = []
        self.group_order = [
            "Run files",
            "Bakta",
            "MLST",
            "MOB-suite",
            "Panaroo",
            "Gubbins",
            "Phylogeny",
            "fastGEAR",
            "chewBBACA",
            "DIAMOND marker DBs",
            "DIAMOND marker hits",
            "DIAMOND loci DBs",
            "DIAMOND loci hits",
        ]
        self.group_labels: dict[str, QtWidgets.QLabel] = {}
        self.setWindowTitle("Run progress")
        self.resize(560, 520)
        self._build_ui()

        self.timer = QtCore.QTimer(self)
        self.timer.setInterval(2000)
        self.timer.timeout.connect(self._refresh)
        self.timer.start()
        self._refresh()

    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(10)

        self.title = QtWidgets.QLabel("Run in progress")
        title_font = self.title.font()
        title_font.setPointSize(max(10, title_font.pointSize() + 2))
        title_font.setBold(True)
        self.title.setFont(title_font)
        layout.addWidget(self.title)

        self.status = QtWidgets.QLabel("")
        self.status.setWordWrap(True)
        self.status.setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignTop)
        self.status.setMinimumHeight(self.fontMetrics().lineSpacing() * 2)
        self.status.setMaximumHeight(self.fontMetrics().lineSpacing() * 3)
        layout.addWidget(self.status)

        self.details = QtWidgets.QPlainTextEdit()
        self.details.setReadOnly(True)
        self.details.setVisible(False)
        self.details.setMaximumHeight(120)
        self.details.setLineWrapMode(QtWidgets.QPlainTextEdit.LineWrapMode.WidgetWidth)
        self.details.setStyleSheet(
            "QPlainTextEdit {"
            " background: #fff7f7;"
            " border: 1px solid #e4b3b3;"
            " border-radius: 6px;"
            " color: #842020;"
            " font-size: 11px;"
            "}"
        )
        layout.addWidget(self.details)

        self.progress = QtWidgets.QProgressBar()
        self.progress.setMinimum(0)
        self.progress.setMaximum(1)
        self.progress.setValue(0)
        layout.addWidget(self.progress)

        self.counts = QtWidgets.QLabel("")
        self.counts.setWordWrap(True)
        self.counts.setVisible(False)
        layout.addWidget(self.counts)

        grid = QtWidgets.QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(12)
        grid.setVerticalSpacing(8)
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 0)
        header_left = QtWidgets.QLabel("Step")
        header_left.setStyleSheet("font-weight: 600; color: #2b3c4d;")
        header_right = QtWidgets.QLabel("Completed")
        header_right.setStyleSheet("font-weight: 600; color: #2b3c4d;")
        grid.addWidget(header_left, 0, 0)
        grid.addWidget(header_right, 0, 1)
        for row, name in enumerate(self.group_order, start=1):
            label = QtWidgets.QLabel(name)
            label.setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignVCenter)
            value = QtWidgets.QLabel("0/0")
            value.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
            grid.addWidget(label, row, 0)
            grid.addWidget(value, row, 1)
            self.group_labels[name] = value
        layout.addLayout(grid)

        btn_row = QtWidgets.QHBoxLayout()
        btn_row.addStretch(1)
        self.close_btn = QtWidgets.QPushButton("Close")
        self.close_btn.setEnabled(True)
        self.close_btn.clicked.connect(self.close)
        btn_row.addWidget(self.close_btn)
        layout.addLayout(btn_row)

    def set_output_dir(self, output_dir: Path) -> None:
        self.output_dir = Path(output_dir)
        self.expected_files = []
        self.title.setText("Run in progress")
        self.status.setText("")
        self.status.setToolTip("")
        self.details.clear()
        self.details.setVisible(False)
        self.progress.setValue(0)
        self.counts.setText("")
        self.counts.setVisible(False)
        for label in self.group_labels.values():
            label.setText("0/0")
        if not self.timer.isActive():
            self.timer.start()
        self._refresh()

    def mark_failed(self, error: str) -> None:
        self.title.setText("Run failed")
        compact = _compact_failure_text(error)
        self.status.setText(compact)
        self.status.setToolTip(error)
        detail = str(error or "").strip()
        self.details.setPlainText(detail)
        self.details.setVisible(bool(detail and detail != compact))
        self.timer.stop()
        self._refresh()

    def mark_finished(self) -> None:
        self._refresh()

    def _mark_completed(self) -> None:
        self.title.setText("Run completed")
        self.status.setText("All expected files have been created.")
        self.status.setToolTip("")
        self.details.clear()
        self.details.setVisible(False)
        self.timer.stop()

    def _load_expected_files(self) -> bool:
        resolved = resolve_expected_files(self.output_dir)
        if not resolved:
            return False
        self.expected_files = resolved
        self.progress.setMaximum(len(self.expected_files))
        return True

    def _classify(self, path: Path) -> Optional[str]:
        return classify_expected_path(path)

    def _expected_exists(self, path: Path) -> bool:
        return expected_output_exists(self.output_dir, path)

    def _refresh(self) -> None:
        status_info = load_run_status(self.output_dir)
        if not self.expected_files:
            if not self._load_expected_files():
                full_status = describe_run_status(status_info, None) or ""
                status_text = _compact_status_text(status_info, pending_text=None, full_status=full_status)
                self.status.setText(status_text)
                self.status.setToolTip(full_status)
                self._set_details(full_status, status_text, str(status_info.get("status") or ""))
                self.title.setText(self._title_for_status(str(status_info.get("status") or "")))
                self.progress.setValue(0)
                self.counts.setText("")
                for label in self.group_labels.values():
                    label.setText("0/0")
                if str(status_info.get("status") or "") in {"success", "failed", "stopped", "interrupted", "partial"}:
                    self.timer.stop()
                return

        existing = 0
        total = len(self.expected_files)
        group_totals: dict[str, int] = {k: 0 for k in self.group_order}
        group_found: dict[str, int] = {k: 0 for k in self.group_order}

        for p in self.expected_files:
            exists = self._expected_exists(p)
            if exists:
                existing += 1
            group = self._classify(p)
            if group and group in group_totals:
                group_totals[group] += 1
                if exists:
                    group_found[group] += 1

        self.progress.setValue(existing)
        self.counts.setText("")
        self.counts.setVisible(False)

        for name in self.group_order:
            total_n = group_totals.get(name, 0)
            found_n = group_found.get(name, 0)
            if total_n == 0:
                self.group_labels[name].setText("—")
            else:
                self.group_labels[name].setText(f"{found_n}/{total_n}")

        pending = summarize_pending_outputs(self.output_dir, self.expected_files) if total > 0 else None
        full_status = describe_run_status(status_info, pending) or ""
        status_text = _compact_status_text(status_info, pending_text=pending, full_status=full_status)
        title = self._title_for_status(str(status_info.get("status") or ""))
        self.title.setText(title)
        self.status.setText(status_text)
        self.status.setToolTip(full_status)
        self._set_details(full_status, status_text, str(status_info.get("status") or ""))

        terminal_statuses = {"success", "failed", "stopped", "interrupted", "partial"}
        if str(status_info.get("status") or "") in terminal_statuses and self.timer.isActive():
            self.timer.stop()

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
        self.details.setPlainText(detail if should_show else "")
        self.details.setVisible(should_show)

    def closeEvent(self, event: QtCore.QEvent) -> None:
        if self.timer.isActive():
            self.timer.stop()
        self.closed.emit()
        super().closeEvent(event)
