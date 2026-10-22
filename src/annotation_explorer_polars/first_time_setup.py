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


# --------------------------------------------------------------------------- #
# Full install or viewer only, and the way from a viewer install to the pipeline
# --------------------------------------------------------------------------- #
CHOICE_FULL = "full"
CHOICE_VIEWER = "viewer"

PIPELINE_NEEDED_TIP = "Needs the pipeline setup"


class _ChoiceCard(QtWidgets.QFrame):
    """One option in the run card style; a click selects it, a double click chooses it."""

    clicked = QtCore.pyqtSignal()
    activated = QtCore.pyqtSignal()

    def __init__(self, title: str, meta: str, text: str, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setObjectName("RunCard")
        self.setProperty("selected", False)
        self.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(16, 12, 16, 12)
        layout.setSpacing(4)
        title_label = QtWidgets.QLabel(title)
        title_label.setProperty("role", "runTitle")
        layout.addWidget(title_label)
        meta_label = QtWidgets.QLabel(meta)
        meta_label.setProperty("role", "runMeta")
        layout.addWidget(meta_label)
        body = QtWidgets.QLabel(text)
        body.setProperty("role", "runStep")
        body.setWordWrap(True)
        layout.addSpacing(2)
        layout.addWidget(body)

    def set_selected(self, selected: bool) -> None:
        self.setProperty("selected", bool(selected))
        self.style().unpolish(self)
        self.style().polish(self)

    def mousePressEvent(self, event: QtGui.QMouseEvent) -> None:  # noqa: N802 - Qt override
        if event.button() == QtCore.Qt.MouseButton.LeftButton:
            self.clicked.emit()
        super().mousePressEvent(event)

    def mouseDoubleClickEvent(self, event: QtGui.QMouseEvent) -> None:  # noqa: N802 - Qt override
        if event.button() == QtCore.Qt.MouseButton.LeftButton:
            self.activated.emit()
        super().mouseDoubleClickEvent(event)


class SetupChoiceDialog(QtWidgets.QDialog):
    """The first start without a config: a full install (tools image and databases) or a
    viewer that only opens runs made elsewhere. `choice` is CHOICE_FULL or CHOICE_VIEWER."""

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setObjectName("SetupChoiceDialog")
        self.setWindowTitle("First-Time Setup")
        self.resize(540, 0)
        self.choice = CHOICE_FULL
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(20, 18, 20, 14)
        layout.setSpacing(10)
        layout.setSizeConstraint(QtWidgets.QLayout.SizeConstraint.SetMinimumSize)

        title = QtWidgets.QLabel("First-Time Setup")
        title.setProperty("role", "dialogTitle")
        layout.addWidget(title)
        hint = QtWidgets.QLabel("Choose what this install does. A viewer can set up the pipeline later.")
        hint.setProperty("role", "muted")
        hint.setWordWrap(True)
        layout.addWidget(hint)
        layout.addSpacing(4)

        self.cards = {
            CHOICE_FULL: _ChoiceCard(
                "Full (run analyses)",
                "Tools image and databases, over 40 GB",
                "Run predictions on your genomes and open any run. Setup downloads the tools image "
                "and the Bakta and CheckM2 databases (about 125 GB free while they unpack).",
            ),
            CHOICE_VIEWER: _ChoiceCard(
                "Viewer only (open results from others)",
                "No downloads",
                "Open runs shared with you (Import Run... in the run list) and explore every view. "
                "Settings can set up the pipeline later.",
            ),
        }
        for key, card in self.cards.items():
            card.setAccessibleName(key)
            card.clicked.connect(lambda key=key: self._select(key))
            card.activated.connect(lambda key=key: (self._select(key), self.accept()))
            layout.addWidget(card)
        self._select(CHOICE_FULL)

        layout.addStretch(0)
        footer = QtWidgets.QFrame()
        footer.setObjectName("ProgressFooter")
        buttons = QtWidgets.QHBoxLayout(footer)
        buttons.setContentsMargins(0, 10, 0, 0)
        buttons.setSpacing(8)
        buttons.addStretch(1)
        cancel = QtWidgets.QPushButton("Cancel")
        cancel.setProperty("role", "flat")
        cancel.clicked.connect(self.reject)
        buttons.addWidget(cancel)
        self.continue_btn = QtWidgets.QPushButton("Continue")
        self.continue_btn.setObjectName("SetupChoiceContinue")
        self.continue_btn.setProperty("role", "primary")
        self.continue_btn.setDefault(True)
        self.continue_btn.clicked.connect(self.accept)
        buttons.addWidget(self.continue_btn)
        layout.addWidget(footer)

    def _select(self, key: str) -> None:
        self.choice = key
        for name, card in self.cards.items():
            card.set_selected(name == key)


def _root_config(root: Path) -> dict:
    from pipeline.Setup import _read_config

    return _read_config(root)


def _pick_runtime(root: Path, *, dev_mode: bool) -> None:
    """A config without a runtime gets the first usable one here, as `bactolomew setup`
    does, so the setup that follows pulls the tools image (a dev setup stays native)."""
    import json
    import os

    from pipeline import ContainerRuntime as CR

    config = _root_config(root)
    if dev_mode or config.get("runtime") or os.environ.get("BACTOLOMEW_NATIVE"):
        return
    runtime = CR.suggest_runtime()
    if not runtime:
        return  # the setup's image row says what is missing
    config["runtime"] = runtime
    path = root / "bactolomew_config.json"
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(config, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(path)


def _run_full_setup(parent: Optional[QtWidgets.QWidget], root: Path, *, dev_mode: bool) -> bool:
    _pick_runtime(root, dev_mode=dev_mode)
    dlg = FirstTimeSetupDialog(parent, dev_mode=dev_mode, root=root)
    return dlg.exec() == QtWidgets.QDialog.DialogCode.Accepted


def ensure_app_ready(*, dev_mode: bool = False, parent: Optional[QtWidgets.QWidget] = None,
                     root: Optional[Path] = None) -> bool:
    """Before the run list: without a config, ask full or viewer first. A viewer install
    writes its config and is ready (nothing to download); a full one runs First-Time Setup
    when anything is missing. False when the user cancels."""
    from pipeline.Setup import collect_setup_missing_paths, write_viewer_config
    from pipeline.ToolEnv import APP_MODE_VIEWER, configured_app_mode, project_root

    root = Path(root) if root is not None else project_root()
    if not (root / "bactolomew_config.json").exists():
        choice = SetupChoiceDialog(parent)
        if choice.exec() != QtWidgets.QDialog.DialogCode.Accepted:
            return False
        if choice.choice == CHOICE_VIEWER:
            write_viewer_config(root)
            return True
        return _run_full_setup(parent, root, dev_mode=dev_mode)
    config = _root_config(root)
    if configured_app_mode(config) == APP_MODE_VIEWER:
        write_viewer_config(root)  # makes runs/ and cache/ when they are missing
        return True
    if not collect_setup_missing_paths(root, config, dev_mode=dev_mode):
        return True
    return _run_full_setup(parent, root, dev_mode=dev_mode)


def offer_pipeline_setup(parent: Optional[QtWidgets.QWidget] = None, *, dev_mode: bool = False,
                         root: Optional[Path] = None, ask: bool = True) -> bool:
    """From a viewer install to the pipeline: ask (unless `ask` is False), then run
    First-Time Setup, which makes the install full. True when the pipeline is set up."""
    from pipeline.ToolEnv import APP_MODE_FULL, configured_app_mode, project_root

    root = Path(root) if root is not None else project_root()
    if configured_app_mode(_root_config(root)) == APP_MODE_FULL:
        return True
    if ask:
        box = QtWidgets.QMessageBox(parent)
        box.setObjectName("PipelineSetupQuestion")
        box.setIcon(QtWidgets.QMessageBox.Icon.NoIcon)
        box.setWindowTitle("Set Up Pipeline")
        box.setText("Running predictions needs the pipeline setup (tools image and databases, over 40 GB). "
                    "Set it up now?")
        go = box.addButton("Set Up", QtWidgets.QMessageBox.ButtonRole.AcceptRole)
        cancel = box.addButton("Cancel", QtWidgets.QMessageBox.ButtonRole.RejectRole)
        for button, role in ((cancel, "flat"), (go, "primary")):
            button.setProperty("role", role)
            button.style().unpolish(button)
            button.style().polish(button)
        box.setDefaultButton(go)
        box.exec()
        if box.clickedButton() is not go:
            return False
    if not _run_full_setup(parent, root, dev_mode=dev_mode):
        return False
    return configured_app_mode(_root_config(root)) == APP_MODE_FULL
