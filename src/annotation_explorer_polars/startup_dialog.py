from __future__ import annotations

import json
import hashlib
import os
import weakref
import socket
import signal
import time
from datetime import datetime
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Callable, List, Optional

from PyQt6 import QtCore, QtWidgets

from .run_progress import RunProgressDialog
from .constants import RUN_CONFIG_FILENAME, RUN_LOCK_FILENAME, RUN_STOP_FILENAME
from .run_expectations import (
    classify_expected_path,
    describe_run_status,
    expected_output_exists,
    load_run_status,
    resolve_expected_files,
    run_allows_resume,
    run_allows_view,
    run_is_active,
    summarize_pending_outputs,
)

if TYPE_CHECKING:
    from pipeline.DataInput import DataInput


@dataclass
class RunEntry:
    run_id: str
    sample_set_id: str
    sample_count: Optional[int]
    path: Path
    run_config: Optional[dict]
    expected_total: int
    expected_found: int
    running: bool
    status: str
    status_text: Optional[str]
    can_view: bool
    can_resume: bool
    pending_text: Optional[str] = None
    status_detail: Optional[str] = None


class PipelineRunner(QtCore.QThread):
    finished_ok = QtCore.pyqtSignal()
    failed = QtCore.pyqtSignal(str)

    def __init__(self, cfg: "DataInput", run_id: str, output_dir: Path, parent: Optional[QtCore.QObject] = None) -> None:
        super().__init__(parent)
        self.cfg = cfg
        self.run_id = run_id
        self.output_dir = Path(output_dir)

    def run(self) -> None:
        from pipeline.Runner import run_workflow

        try:
            run_workflow(self.cfg)
        except BaseException as exc:
            self.failed.emit(str(exc) or exc.__class__.__name__)
            return
        self.finished_ok.emit()


_ACTIVE_RUNNERS: list[PipelineRunner] = []


def _project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def runs_root() -> Path:
    return _project_root() / "runs"


def _load_app_config() -> dict:
    path = _project_root() / "bactolomew_config.json"
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text())
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


def _ensure_app_config(
    parent: Optional[QtWidgets.QWidget] = None,
    *,
    dev_mode: bool = False,
) -> Optional[dict]:
    from .first_time_setup import FirstTimeSetupDialog, apply_first_time_setup, collect_setup_missing_paths

    root = _project_root()
    config_path = root / "bactolomew_config.json"
    config = _load_app_config()

    missing: list[Path] = []
    if not config_path.exists():
        missing.append(config_path)
    missing.extend(collect_setup_missing_paths(root, config, dev_mode=dev_mode))

    if not missing:
        return config

    dlg = FirstTimeSetupDialog(missing=missing, parent=parent, dev_mode=dev_mode)
    if dlg.exec() != QtWidgets.QDialog.DialogCode.Accepted:
        return None
    result = dlg.result()
    if result is None:
        return None
    try:
        apply_first_time_setup(result)
    except Exception as exc:
        QtWidgets.QMessageBox.critical(parent, "Setup failed", str(exc))
        return None
    return _load_app_config()


def _compute_run_id(genomes: list[Path]) -> str:
    samples = sorted({p.stem for p in genomes if p is not None})
    joined = "|".join(samples)
    digest = hashlib.sha1(joined.encode("utf-8")).hexdigest()[:12]
    return f"sample_{digest}"


def _new_attempt_id(sample_set_id: str) -> str:
    stamp = datetime.now().astimezone().strftime("%Y%m%d_%H%M%S_%f")
    return f"{sample_set_id}__run_{stamp}"


def _pid_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    return True


def _load_catalog_entries(path: Path) -> list[tuple[Path, str]]:
    from pipeline.DbCatalog import load_db_catalog

    try:
        catalog = load_db_catalog(path)
    except Exception as exc:
        QtWidgets.QMessageBox.critical(None, "Catalog error", f"Failed to load catalog: {path}\n{exc}")
        return []
    entries: list[tuple[Path, str]] = []
    for _name, info in (catalog or {}).items():
        if not isinstance(info, dict):
            continue
        p_raw = info.get("db_path") or info.get("path") or info.get("db-path")
        if not p_raw:
            continue
        p = Path(p_raw).expanduser().resolve()
        cat = str(info.get("category") or "").strip()
        entries.append((p, cat))
    return entries


def _find_active_run(exclude: Optional[Path] = None) -> Optional[Path]:
    runs_dir = runs_root()
    if not runs_dir.exists():
        return None
    for lock_path in runs_dir.rglob(RUN_LOCK_FILENAME):
        run_dir = lock_path.parent
        if exclude and run_dir.resolve() == exclude.resolve():
            continue
        if run_is_active(run_dir):
            return run_dir
    return None


def _compact_run_status(
    status_info: dict[str, object],
    *,
    pending_text: Optional[str],
    full_status: Optional[str],
) -> str:
    status = str(status_info.get("status") or "").strip().lower()
    current_step = str(status_info.get("current_step") or "").strip()
    if status == "running":
        return f"Running: {current_step}" if current_step else "Running"
    if status == "success":
        return "Completed"
    if status == "stopped":
        return "Stopped. Click Resume to continue."
    if status == "interrupted":
        return "Interrupted. Click Resume to continue."
    if status == "failed":
        return "Failed. Click Resume to retry."
    if status == "partial":
        if pending_text:
            return pending_text
        return "Partial run. Click Resume to continue."
    return full_status or pending_text or "Status unknown. Click Resume to continue."


# --- shared control styling (matches the app's light / #2f6fed design language) ---
_PRIMARY_BTN = (
    "QPushButton{background:#2f6fed;color:#ffffff;border:none;border-radius:8px;"
    "padding:6px 18px;font-weight:600;}"
    "QPushButton:hover{background:#295fd0;}"
    "QPushButton:disabled{background:#c6d3e6;color:#eef3f8;}"
)
_SECONDARY_BTN = (
    "QPushButton{background:#ffffff;color:#1f2c3b;border:1px solid #cdd6e0;border-radius:8px;"
    "padding:6px 14px;font-weight:600;}"
    "QPushButton:hover{background:#f2f6fb;border-color:#b9c9de;}"
    "QPushButton:disabled{color:#aab6c4;border-color:#e6ecf3;}"
)
_GHOST_BTN = (
    "QPushButton{background:transparent;color:#64748b;border:none;border-radius:8px;"
    "padding:6px 14px;font-weight:600;}"
    "QPushButton:hover{background:#e9eef4;color:#1f2c3b;}"
)


def _pill_style(fg: str, r: int, g: int, b: int) -> str:
    return (
        f"color:{fg};background:rgba({r},{g},{b},0.13);border:1px solid rgba({r},{g},{b},0.32);"
        "border-radius:9px;padding:2px 10px;font-weight:700;font-size:10px;"
    )


class RunCard(QtWidgets.QFrame):
    def __init__(
        self,
        entry: RunEntry,
        on_view: Optional[Callable[[RunEntry], None]] = None,
        on_monitor: Optional[Callable[[RunEntry], None]] = None,
        on_resume: Optional[Callable[[RunEntry], None]] = None,
        on_stop: Optional[Callable[[RunEntry], None]] = None,
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.setObjectName("RunCard")
        self.setStyleSheet(
            """
            QFrame#RunCard {
                background: #ffffff;
                border: 1px solid #dbe3ec;
                border-radius: 14px;
            }
            QFrame#RunCard:hover {
                border: 1px solid #b9c9de;
                background: #fcfdff;
            }
            """
        )
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(16, 14, 16, 14)
        layout.setSpacing(9)
        title = QtWidgets.QLabel(entry.sample_set_id or entry.run_id)
        title.setStyleSheet("font-weight: 700; font-size: 14px; color: #172433;")
        attempt_label = entry.run_id
        if entry.sample_count is None:
            sample_text = "Samples: unknown"
        else:
            sample_text = f"{entry.sample_count} samples"
        subtitle_bits = [attempt_label, sample_text]
        subtitle = QtWidgets.QLabel(" • ".join(bit for bit in subtitle_bits if bit))
        subtitle.setStyleSheet("color: #64748b; font-size: 11px;")
        layout.addWidget(title)
        layout.addWidget(subtitle)

        progress = QtWidgets.QProgressBar()
        progress.setMinimum(0)
        progress.setMaximum(100)
        pct = 0
        if entry.expected_total > 0:
            if entry.expected_found >= entry.expected_total:
                pct = 100
            else:
                pct = min(99, int((entry.expected_found * 100) / entry.expected_total))
        elif entry.can_view:
            pct = 100
        progress.setValue(max(0, min(100, pct)))
        progress.setTextVisible(False)
        progress.setFixedHeight(8)
        progress.setStyleSheet(
            "QProgressBar {"
            "  background: #e6ecf4;"
            "  border: none;"
            "  border-radius: 4px;"
            "}"
            "QProgressBar::chunk {"
            "  background: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #5b8bf2, stop:1 #2f6fed);"
            "  border-radius: 4px;"
            "}"
        )
        layout.addSpacing(2)
        layout.addWidget(progress)

        if entry.expected_total > 0:
            frac_text = f"{entry.expected_found} / {entry.expected_total}"
        else:
            frac_text = f"{pct}%"

        status_text = ""
        pill_style = ""
        if entry.running:
            status_text = entry.status_text or "Running"
            pill_style = _pill_style("#2f6fed", 47, 111, 237)
        elif entry.status == "success":
            status_text = entry.status_text or "Completed"
            pill_style = _pill_style("#0f7a37", 15, 157, 107)
        elif entry.status in {"failed", "interrupted"}:
            status_text = entry.status_text or "Failed"
            pill_style = _pill_style("#c0322b", 220, 38, 38)
        elif entry.status == "stopped":
            status_text = entry.status_text or "Stopped"
            pill_style = _pill_style("#a86612", 217, 119, 6)
        elif entry.status_text:
            status_text = entry.status_text
            pill_style = _pill_style("#a86612", 217, 119, 6)

        status_row = QtWidgets.QHBoxLayout()
        status_row.setSpacing(6)
        if status_text:
            status = QtWidgets.QLabel(status_text)
            status.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.NoTextInteraction)
            status.setSizePolicy(QtWidgets.QSizePolicy.Policy.Maximum, QtWidgets.QSizePolicy.Policy.Fixed)
            status.setStyleSheet(pill_style)
            detail = entry.status_detail or entry.status_text
            if detail and detail != status_text:
                status.setToolTip(detail)
            status_row.addWidget(status)
        status_row.addStretch(1)
        frac = QtWidgets.QLabel(frac_text)
        frac.setStyleSheet("color:#8a99ab; font-size:11px; font-weight:600;")
        status_row.addWidget(frac)
        layout.addLayout(status_row)

        btn_row = QtWidgets.QHBoxLayout()
        btn_row.addStretch(1)
        view_btn = QtWidgets.QPushButton("View")
        view_btn.setEnabled(entry.can_view)
        view_btn.setStyleSheet(_PRIMARY_BTN)
        view_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        if not entry.can_view:
            view_btn.setToolTip(
                entry.status_detail
                or entry.status_text
                or entry.pending_text
                or "This run is not ready to open."
            )
        view_btn.clicked.connect(lambda: on_view(entry) if on_view else None)
        btn_row.addWidget(view_btn)

        if not entry.can_view:
            monitor_btn = QtWidgets.QPushButton("Monitor")
            monitor_btn.clicked.connect(lambda: on_monitor(entry) if on_monitor else None)
            resume_btn = QtWidgets.QPushButton("Resume")
            resume_btn.setEnabled(entry.can_resume)
            resume_btn.clicked.connect(lambda: on_resume(entry) if on_resume else None)
            stop_btn = QtWidgets.QPushButton("Stop")
            stop_btn.setEnabled(entry.running)
            stop_btn.clicked.connect(lambda: on_stop(entry) if on_stop else None)
            force_btn = QtWidgets.QPushButton("Force Stop")
            force_btn.setEnabled(entry.running)
            force_btn.clicked.connect(lambda: on_stop(entry, True) if on_stop else None)
            for _b in (monitor_btn, resume_btn, stop_btn, force_btn):
                _b.setStyleSheet(_SECONDARY_BTN)
                _b.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
            btn_row.addWidget(monitor_btn)
            btn_row.addWidget(resume_btn)
            btn_row.addWidget(stop_btn)
            btn_row.addWidget(force_btn)
        layout.addLayout(btn_row)


class NewRunDialog(QtWidgets.QDialog):
    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.genome_paths: List[Path] = []

        self.setWindowTitle("Start New Run")
        self.resize(720, 520)
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)

        intro = QtWidgets.QLabel(
            "Select genomes to run the full required pipeline."
        )
        intro.setWordWrap(True)
        layout.addWidget(intro)

        genomes_group = QtWidgets.QGroupBox("Genomes")
        genomes_layout = QtWidgets.QVBoxLayout(genomes_group)
        self.genome_list = QtWidgets.QListWidget()
        genomes_layout.addWidget(self.genome_list)

        btn_row = QtWidgets.QHBoxLayout()
        add_btn = QtWidgets.QPushButton("Add Genomes…")
        add_btn.clicked.connect(self._add_genomes)
        remove_btn = QtWidgets.QPushButton("Remove Selected")
        remove_btn.clicked.connect(self._remove_selected_genomes)
        clear_btn = QtWidgets.QPushButton("Clear")
        clear_btn.clicked.connect(self._clear_genomes)
        btn_row.addWidget(add_btn)
        btn_row.addWidget(remove_btn)
        btn_row.addWidget(clear_btn)
        btn_row.addStretch(1)
        genomes_layout.addLayout(btn_row)
        layout.addWidget(genomes_group)

        btns = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel
        )
        btns.accepted.connect(self._validate_and_accept)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)

    def _add_genomes(self) -> None:
        paths, _ = QtWidgets.QFileDialog.getOpenFileNames(
            self,
            "Select genome FASTA files",
            "",
            "FASTA files (*.fa *.fasta *.fna *.fa.gz *.fasta.gz *.fna.gz);;All files (*)",
        )
        if not paths:
            return
        for p in paths:
            path = Path(p)
            if path not in self.genome_paths:
                self.genome_paths.append(path)
                self.genome_list.addItem(str(path))

    def _remove_selected_genomes(self) -> None:
        selected = self.genome_list.selectedItems()
        if not selected:
            return
        for item in selected:
            row = self.genome_list.row(item)
            self.genome_list.takeItem(row)
        self.genome_paths = [Path(self.genome_list.item(i).text()) for i in range(self.genome_list.count())]

    def _clear_genomes(self) -> None:
        self.genome_list.clear()
        self.genome_paths = []

    def _validate_and_accept(self) -> None:
        if not self.genome_paths:
            QtWidgets.QMessageBox.critical(self, "Validation error", "Please select at least one genome.")
            return
        self.accept()


class StartupDialog(QtWidgets.QDialog):
    def __init__(self, parent: Optional[QtWidgets.QWidget] = None, *, dev_mode: bool = False) -> None:
        super().__init__(parent)
        self._dev_mode = bool(dev_mode)
        self.output_dir: Optional[Path] = None
        self.run_config: Optional[dict] = None
        self.run_id: Optional[str] = None
        self._run_entries: List[RunEntry] = []
        self._pipeline_runner: Optional[PipelineRunner] = None
        self._progress_dialog: Optional[RunProgressDialog] = None
        self._active_run_dir: Optional[Path] = None
        self._active_run_id: Optional[str] = None

        self.setWindowTitle("Select Run")
        self.resize(640, 520)
        self.setObjectName("StartupDialog")
        self.setStyleSheet("QDialog#StartupDialog { background: #eef3f8; }")
        self._build_ui()
        self._load_runs()

        self._refresh_timer = QtCore.QTimer(self)
        self._refresh_timer.setInterval(5000)
        self._refresh_timer.timeout.connect(self._load_runs)
        self._refresh_timer.start()

    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(22, 20, 22, 18)
        layout.setSpacing(12)

        title = QtWidgets.QLabel("Select a run")
        title.setStyleSheet("font-size: 18px; font-weight: 800; color: #172433;")
        layout.addWidget(title)
        intro = QtWidgets.QLabel("Open a completed run in the annotation explorer, or start a new one.")
        intro.setWordWrap(True)
        intro.setStyleSheet("color: #64748b; font-size: 12px;")
        layout.addWidget(intro)

        list_row = QtWidgets.QHBoxLayout()
        list_row.addStretch(1)
        self.run_list = QtWidgets.QListWidget()
        self.run_list.setObjectName("RunList")
        self.run_list.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.run_list.setSpacing(8)
        self.run_list.setFixedSize(520, 360)
        self.run_list.setStyleSheet(
            """
            QListWidget#RunList {
                background: transparent;
                border: none;
                padding: 4px;
            }
            QListWidget#RunList::item:selected {
                background: transparent;
            }
            QListWidget#RunList QScrollBar:vertical {
                background: transparent;
                width: 9px;
                margin: 2px 0;
            }
            QListWidget#RunList QScrollBar::handle:vertical {
                background: rgba(31,45,61,0.24);
                border-radius: 4px;
                min-height: 28px;
            }
            QListWidget#RunList QScrollBar::handle:vertical:hover {
                background: rgba(31,45,61,0.40);
            }
            QListWidget#RunList QScrollBar::add-line:vertical,
            QListWidget#RunList QScrollBar::sub-line:vertical {
                height: 0;
            }
            QListWidget#RunList QScrollBar::add-page:vertical,
            QListWidget#RunList QScrollBar::sub-page:vertical {
                background: transparent;
            }
            """
        )
        self.run_list.currentItemChanged.connect(self._on_run_selected)
        list_row.addWidget(self.run_list, stretch=0)
        list_row.addStretch(1)
        layout.addLayout(list_row)

        btn_row = QtWidgets.QHBoxLayout()
        new_run_btn = QtWidgets.QPushButton("Start New Run")
        new_run_btn.setStyleSheet(_SECONDARY_BTN)
        new_run_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        new_run_btn.clicked.connect(self._open_new_run_dialog)
        btn_row.addWidget(new_run_btn)
        btn_row.addStretch(1)
        cancel_btn = QtWidgets.QPushButton("Cancel")
        cancel_btn.setStyleSheet(_GHOST_BTN)
        cancel_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        cancel_btn.clicked.connect(self.reject)
        btn_row.addWidget(cancel_btn)
        layout.addLayout(btn_row)

    def _load_runs(self) -> None:
        prev_run_id = self.run_id
        prev_path = str(self.output_dir.resolve()) if self.output_dir else ""
        self.run_list.clear()
        self._run_entries = []
        runs_dir = runs_root()
        if not runs_dir.exists():
            msg = QtWidgets.QListWidgetItem(f"No runs found (missing {runs_dir}).")
            msg.setFlags(QtCore.Qt.ItemFlag.NoItemFlags)
            self.run_list.addItem(msg)
            return

        run_cfg_paths = sorted(runs_dir.rglob(RUN_CONFIG_FILENAME), reverse=True)
        seen_dirs: set[str] = set()
        for run_cfg_path in run_cfg_paths:
            run_cfg = None
            try:
                payload = json.loads(run_cfg_path.read_text(encoding="utf-8"))
            except Exception:
                payload = None
            if isinstance(payload, dict):
                run_cfg = payload
                output_dir = Path(payload.get("output_dir") or run_cfg_path.parent).expanduser().resolve()
            else:
                output_dir = run_cfg_path.parent.resolve()
            if str(output_dir) in seen_dirs:
                continue
            seen_dirs.add(str(output_dir))

            expected_files = resolve_expected_files(output_dir)
            expected_total = len(expected_files)
            expected_found = sum(1 for fp in expected_files if expected_output_exists(output_dir, fp))

            sample_count = None
            if isinstance(run_cfg, dict):
                genomes = run_cfg.get("inputs", {}).get("genomes") if isinstance(run_cfg.get("inputs"), dict) else None
                if isinstance(genomes, list):
                    sample_count = len(genomes)

            status_info = load_run_status(output_dir)
            running = bool(status_info.get("active"))
            pending_text = summarize_pending_outputs(output_dir, expected_files) if expected_files else None
            can_view = run_allows_view(status_info)
            can_resume = run_allows_resume(status_info)
            status_detail = describe_run_status(status_info, pending_text)
            status_text = _compact_run_status(
                status_info,
                pending_text=pending_text,
                full_status=status_detail,
            )
            run_id = str(
                (run_cfg or {}).get("run_id")
                or status_info.get("run_id")
                or output_dir.name
            ).strip() or output_dir.name
            sample_set_id = str(
                (run_cfg or {}).get("sample_set_id")
                or status_info.get("sample_set_id")
                or output_dir.parent.name
            ).strip() or output_dir.parent.name
            entry = RunEntry(
                run_id=run_id,
                sample_set_id=sample_set_id,
                sample_count=sample_count,
                path=output_dir,
                run_config=run_cfg if isinstance(run_cfg, dict) else None,
                expected_total=expected_total,
                expected_found=expected_found,
                running=running,
                status=str(status_info.get("status") or "unknown"),
                status_text=status_text,
                can_view=can_view,
                can_resume=can_resume,
                pending_text=pending_text,
                status_detail=status_detail,
            )
            self._run_entries.append(entry)

        if not self._run_entries:
            msg = QtWidgets.QListWidgetItem("No runs found in the runs/ folder.")
            msg.setFlags(QtCore.Qt.ItemFlag.NoItemFlags)
            self.run_list.addItem(msg)
            return

        for entry in self._run_entries:
            item = QtWidgets.QListWidgetItem()
            item.setData(QtCore.Qt.ItemDataRole.UserRole, entry)
            card = RunCard(
                entry,
                self._open_run_from_entry,
                self._open_progress_for_entry,
                self._resume_run,
                self._stop_run,
                self.run_list,
            )
            card.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed)
            card_width = max(240, self.run_list.viewport().width() - 8)
            card.setFixedWidth(card_width)
            item.setSizeHint(card.sizeHint())
            self.run_list.addItem(item)
            self.run_list.setItemWidget(item, card)

        if self.run_list.count() > 0:
            if prev_path:
                for i in range(self.run_list.count()):
                    item = self.run_list.item(i)
                    entry = item.data(QtCore.Qt.ItemDataRole.UserRole)
                    if isinstance(entry, RunEntry) and str(entry.path.resolve()) == prev_path:
                        self.run_list.setCurrentRow(i)
                        break
                else:
                    self.run_list.setCurrentRow(0)
            elif prev_run_id:
                for i in range(self.run_list.count()):
                    item = self.run_list.item(i)
                    entry = item.data(QtCore.Qt.ItemDataRole.UserRole)
                    if isinstance(entry, RunEntry) and entry.run_id == prev_run_id:
                        self.run_list.setCurrentRow(i)
                        break
                else:
                    self.run_list.setCurrentRow(0)
            else:
                self.run_list.setCurrentRow(0)

    def _on_run_selected(
        self,
        current: Optional[QtWidgets.QListWidgetItem],
        _previous: Optional[QtWidgets.QListWidgetItem],
    ) -> None:
        entry = current.data(QtCore.Qt.ItemDataRole.UserRole) if current else None
        if not isinstance(entry, RunEntry):
            self.output_dir = None
            self.run_config = None
            self.run_id = None
            return
        self.output_dir = entry.path
        self.run_config = entry.run_config
        self.run_id = entry.run_id

    def _selected_run_entry(self) -> Optional[RunEntry]:
        if not self.output_dir:
            return None
        out = self.output_dir.resolve()
        for entry in self._run_entries:
            if entry.path.resolve() == out:
                return entry
        return None

    def _open_run_from_entry(self, entry: RunEntry) -> None:
        if not entry.can_view:
            QtWidgets.QMessageBox.information(
                self,
                "Run not ready",
                entry.status_text or entry.pending_text or "Only successful runs can be opened in the visualizer.",
            )
            return
        self.output_dir = entry.path
        self.run_config = entry.run_config
        self.run_id = entry.run_id
        self._validate_and_accept()

    def _validate_and_accept(self) -> None:
        errors = self._validate_inputs()
        if errors:
            QtWidgets.QMessageBox.critical(self, "Validation error", "\n".join(errors))
            return
        self.accept()

    def _open_progress_dialog(self, output_dir: Path) -> RunProgressDialog:
        if self._progress_dialog and self._progress_dialog.isVisible():
            self._progress_dialog.set_output_dir(output_dir)
            self._progress_dialog.raise_()
            self._progress_dialog.activateWindow()
            return self._progress_dialog
        dlg = RunProgressDialog(output_dir, parent=self)
        self._progress_dialog = dlg
        dlg.closed.connect(self._clear_progress_dialog)
        dlg.show()
        return dlg

    def _clear_progress_dialog(self) -> None:
        self._progress_dialog = None

    def _open_progress_for_entry(self, entry: RunEntry) -> None:
        self._open_progress_dialog(entry.path)

    def _open_progress_for_selected(self) -> None:
        if self.output_dir:
            self._open_progress_dialog(self.output_dir)
            return
        if self._active_run_dir:
            self._open_progress_dialog(self._active_run_dir)
            return
        QtWidgets.QMessageBox.information(self, "Monitor run", "Select a run first.")

    def _stop_run(self, entry: RunEntry, force: bool = False) -> None:
        if not entry.running:
            QtWidgets.QMessageBox.information(self, "Stop run", "This run is not currently active.")
            return
        if force:
            lock_path = entry.path / RUN_LOCK_FILENAME
            try:
                data = json.loads(lock_path.read_text(encoding="utf-8")) if lock_path.exists() else {}
            except Exception as exc:
                QtWidgets.QMessageBox.critical(self, "Force stop", f"Failed to read lock file: {exc}")
                return
            host = str(data.get("host") or "").strip()
            pid = int(data.get("pid") or 0)
            if host and host != socket.gethostname():
                QtWidgets.QMessageBox.critical(self, "Force stop", "Run is on a different host.")
                return
            if pid <= 0:
                QtWidgets.QMessageBox.critical(self, "Force stop", "Missing PID in lock file.")
                return
            try:
                os.kill(pid, signal.SIGTERM)
                time.sleep(0.5)
                if _pid_alive(pid):
                    os.kill(pid, signal.SIGKILL)
            except Exception as exc:
                QtWidgets.QMessageBox.critical(self, "Force stop", f"Failed to stop process: {exc}")
                return
            QtWidgets.QMessageBox.information(self, "Force stop", "Force stop signal sent.")
            return

        stop_path = entry.path / RUN_STOP_FILENAME
        try:
            stop_path.write_text("stop\n", encoding="utf-8")
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "Stop run", f"Failed to request stop: {exc}")
            return
        QtWidgets.QMessageBox.information(self, "Stop run", "Stop requested. The run will stop after the current step.")

    def _data_input_from_run_config(self, run_cfg: dict, output_dir: Path) -> Optional["DataInput"]:
        from pipeline.ChewBBACAWorkflow import validate_chewbbaca_settings
        from pipeline.DataInput import DataInput
        from pipeline.ToolEnv import (
            configured_assembly_qc_settings,
            configured_chewbbaca_settings,
            configured_fastgear_settings,
            configured_mobsuite_settings,
        )

        inputs = run_cfg.get("inputs") if isinstance(run_cfg, dict) else None
        if not isinstance(inputs, dict):
            QtWidgets.QMessageBox.critical(self, "Run config error", "run_config.json is missing inputs.")
            return None
        expected_files = resolve_expected_files(output_dir)
        missing_groups = {
            group
            for path in expected_files
            if not expected_output_exists(output_dir, path)
            for group in [classify_expected_path(path)]
            if group is not None
        }
        try:
            cfg = DataInput(
                genomes=[Path(p) for p in inputs.get("genomes") or []],
                marker_sequences=[Path(p) for p in inputs.get("marker_sequences") or []],
                reference_locus_sequences=[Path(p) for p in inputs.get("reference_locus_sequences") or []],
                output=output_dir,
                snippy_reference=inputs.get("snippy_reference"),
                bakta_db=inputs.get("bakta_db"),
                dmnd_threads=inputs.get("dmnd_threads", inputs.get("cores", 1)),
                bakta_threads=inputs.get("bakta_threads", inputs.get("cores", 1)),
                window_length=int(inputs.get("window_length", 500)),
                diamond=inputs.get("diamond", True),
                bakta=inputs.get("bakta", True),
                total_cores=int(inputs.get("cores", 1)),
                marker_categories=inputs.get("marker_categories") or {},
                reference_locus_categories=inputs.get("reference_locus_categories") or {},
                fastgear_config=inputs.get("fastgear") or configured_fastgear_settings(),
                mobsuite_config=configured_mobsuite_settings(inputs.get("mobsuite")),
                chewbbaca_config=configured_chewbbaca_settings(inputs.get("chewbbaca")),
                assembly_qc_config=configured_assembly_qc_settings(inputs.get("assembly_qc")),
                run_id=str(run_cfg.get("run_id") or output_dir.name),
            )
            validate_chewbbaca_settings(cfg.chewbbaca_config)
            if "fastGEAR" in missing_groups and not cfg.fastgear_config.get("image"):
                raise ValueError(
                    "fastGEAR Docker image not configured. Set bactolomew_config.json fastgear.image before resuming."
                )
            return cfg
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "Run config error", str(exc))
            return None

    def _resume_run(self, entry: RunEntry) -> None:
        active = _find_active_run(exclude=entry.path)
        if active is not None:
            QtWidgets.QMessageBox.warning(
                self,
                "Run in progress",
                f"A pipeline run is already in progress in:\n{active}\nWait for it to finish before resuming another.",
            )
            return
        if self._pipeline_runner and self._pipeline_runner.isRunning():
            QtWidgets.QMessageBox.warning(
                self,
                "Run in progress",
                "A pipeline run is already in progress. Wait for it to finish before resuming another.",
            )
            return

        run_cfg = entry.run_config
        if not isinstance(run_cfg, dict):
            run_cfg_path = entry.path / RUN_CONFIG_FILENAME
            if not run_cfg_path.exists():
                QtWidgets.QMessageBox.critical(
                    self,
                    "Missing config",
                    f"{RUN_CONFIG_FILENAME} not found for this run. Cannot resume.",
                )
                return
            try:
                run_cfg = json.loads(run_cfg_path.read_text(encoding="utf-8"))
            except Exception as exc:
                QtWidgets.QMessageBox.critical(self, "Run config error", str(exc))
                return

        cfg = self._data_input_from_run_config(run_cfg, entry.path)
        if cfg is None:
            return

        self._active_run_dir = entry.path
        self._active_run_id = entry.run_id

        progress = self._open_progress_dialog(entry.path)
        runner = PipelineRunner(cfg, run_id=entry.run_id, output_dir=entry.path)
        self._pipeline_runner = runner
        _ACTIVE_RUNNERS.append(runner)
        runner_ref = weakref.ref(runner)

        def _cleanup_runner() -> None:
            r = runner_ref()
            if r and r in _ACTIVE_RUNNERS:
                _ACTIVE_RUNNERS.remove(r)

        runner.finished_ok.connect(_cleanup_runner)
        runner.failed.connect(lambda _msg: _cleanup_runner())
        runner.failed.connect(progress.mark_failed)
        runner.finished_ok.connect(progress.mark_finished)
        runner.failed.connect(self._on_pipeline_failed)
        runner.finished_ok.connect(self._on_pipeline_finished)
        runner.start()

    def _on_pipeline_finished(self) -> None:
        runner = self.sender()
        if isinstance(runner, PipelineRunner):
            run_id = runner.run_id
            output_dir = runner.output_dir
        else:
            run_id = self._active_run_id or ""
            output_dir = self._active_run_dir or Path()
        self._pipeline_runner = None
        self._active_run_id = run_id
        self._active_run_dir = output_dir if output_dir else self._active_run_dir
        self._load_runs()
        for i in range(self.run_list.count()):
            item = self.run_list.item(i)
            entry = item.data(QtCore.Qt.ItemDataRole.UserRole)
            if isinstance(entry, RunEntry) and entry.run_id == run_id:
                self.run_list.setCurrentRow(i)
                break

    def _on_pipeline_failed(self, error: str) -> None:
        runner = self.sender()
        if isinstance(runner, PipelineRunner):
            output_dir = runner.output_dir
            self._active_run_id = runner.run_id
        else:
            output_dir = self._active_run_dir
        self._pipeline_runner = None
        if output_dir:
            self._active_run_dir = output_dir
        self._load_runs()
        if error.strip() == "Run stopped by user":
            QtWidgets.QMessageBox.information(self, "Run stopped", error)
        else:
            QtWidgets.QMessageBox.critical(self, "Run failed", error)

    def _validate_inputs(self) -> List[str]:
        errors: List[str] = []
        if not self.output_dir:
            errors.append("Please select a run.")

        selected = self._selected_run_entry()
        if selected is not None:
            if not selected.can_view:
                errors.append(
                    selected.status_text
                    or selected.pending_text
                    or "Run is not ready yet. Only successful runs can be opened."
                )

        if errors:
            return errors

        output_dir = self.output_dir
        if not output_dir.exists() or not output_dir.is_dir():
            errors.append("Output folder does not exist or is not a directory.")
            return errors

        from .data import validate_output_dir

        # Fast checks only, so "View" accepts instantly and the loading splash
        # appears immediately. The deep per-file output scan (slow on large runs)
        # runs in the background loader when building from outputs.
        run_cfg, out_errors = validate_output_dir(output_dir, deep=False)
        errors.extend(out_errors)
        if not errors:
            self.run_config = run_cfg
        return errors

    def _open_new_run_dialog(self) -> None:
        config = _ensure_app_config(self, dev_mode=self._dev_mode)
        if config is None:
            return
        dlg = NewRunDialog(self)
        if dlg.exec() != QtWidgets.QDialog.DialogCode.Accepted:
            return
        self._start_pipeline_run(dlg, config)

    def _start_pipeline_run(self, dlg: NewRunDialog, config: dict) -> None:
        active = _find_active_run()
        if active is not None:
            QtWidgets.QMessageBox.warning(
                self,
                "Run in progress",
                f"A pipeline run is already in progress in:\n{active}\nWait for it to finish before starting another.",
            )
            return
        if self._pipeline_runner and self._pipeline_runner.isRunning():
            QtWidgets.QMessageBox.warning(
                self,
                "Run in progress",
                "A pipeline run is already in progress. Wait for it to finish before starting another.",
            )
            return
        genomes = list(dlg.genome_paths or [])
        if not genomes:
            QtWidgets.QMessageBox.warning(self, "Missing", "Add ≥1 genome file.")
            return
        for p in genomes:
            if not p.exists() or not p.is_file():
                QtWidgets.QMessageBox.critical(self, "Error", f"Genome file not found: {p}")
                return

        sample_set_id = _compute_run_id(genomes)
        run_id = _new_attempt_id(sample_set_id)
        output_dir = runs_root() / sample_set_id / run_id

        from pipeline.ChewBBACAWorkflow import validate_chewbbaca_settings
        from pipeline.DataInput import DataInput
        from pipeline.DbCatalog import COMPOUND_DB_CATALOG_PATH, MARKER_DB_CATALOG_PATH
        from pipeline.ToolEnv import (
            configured_assembly_qc_settings,
            configured_bakta_db_path,
            configured_chewbbaca_settings,
            configured_fastgear_settings,
            configured_mobsuite_settings,
        )

        marker_entries = _load_catalog_entries(MARKER_DB_CATALOG_PATH)
        locus_entries = _load_catalog_entries(COMPOUND_DB_CATALOG_PATH)
        markers = [p for p, _ in marker_entries]
        loci = [p for p, _ in locus_entries]

        if not markers:
            QtWidgets.QMessageBox.warning(self, "Missing", "No marker sequences found in the catalog.")
            return
        if not loci:
            QtWidgets.QMessageBox.warning(self, "Missing", "No reference locus sequences found in the catalog.")
            return

        missing_markers = [p for p in markers if not p.exists() or not p.is_file()]
        if missing_markers:
            QtWidgets.QMessageBox.critical(
                self,
                "Missing files",
                "Marker sequences not found:\n" + "\n".join(str(p) for p in missing_markers),
            )
            return
        missing_loci = [p for p in loci if not p.exists() or not p.is_file()]
        if missing_loci:
            QtWidgets.QMessageBox.critical(
                self,
                "Missing files",
                "Reference locus sequences not found:\n" + "\n".join(str(p) for p in missing_loci),
            )
            return

        marker_cats = {p.stem: cat for p, cat in marker_entries}
        locus_cats = {p.stem: cat for p, cat in locus_entries}

        bdb = (config.get("bakta_db") or "").strip()
        if self._dev_mode and not bdb:
            bdb = configured_bakta_db_path()
        if not bdb:
            QtWidgets.QMessageBox.critical(
                self,
                "Missing Bakta DB",
                "Bakta DB path not found in bactolomew_config.json.",
            )
            return
        if not Path(bdb).expanduser().is_dir():
            QtWidgets.QMessageBox.critical(self, "Error", f"Bakta DB folder does not exist: {bdb}")
            return

        fastgear_config = configured_fastgear_settings(config.get("fastgear"))
        if not fastgear_config.get("image"):
            QtWidgets.QMessageBox.critical(
                self,
                "Missing fastGEAR config",
                "fastGEAR Docker image not found in bactolomew_config.json (fastgear.image).",
            )
            return
        mobsuite_config = configured_mobsuite_settings(config.get("mobsuite"))
        chewbbaca_config = configured_chewbbaca_settings(config.get("chewbbaca"))
        assembly_qc_config = configured_assembly_qc_settings(config.get("assembly_qc"))
        try:
            validate_chewbbaca_settings(chewbbaca_config)
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "Missing chewBBACA config", str(exc))
            return

        cpu = os.cpu_count() or 1

        try:
            cfg = DataInput(
                genomes=genomes,
                marker_sequences=markers,
                reference_locus_sequences=loci,
                output=output_dir.expanduser().resolve(),
                snippy_reference=None,
                bakta_db=bdb,
                dmnd_threads=cpu,
                bakta_threads=cpu,
                window_length=500,
                diamond=True,
                bakta=True,
                total_cores=cpu,
                marker_categories=marker_cats,
                reference_locus_categories=locus_cats,
                fastgear_config=fastgear_config,
                mobsuite_config=mobsuite_config,
                chewbbaca_config=chewbbaca_config,
                assembly_qc_config=assembly_qc_config,
                run_id=run_id,
            )
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "Invalid input", str(exc))
            return

        self._active_run_dir = output_dir
        self._active_run_id = run_id

        progress = self._open_progress_dialog(output_dir)
        runner = PipelineRunner(cfg, run_id=run_id, output_dir=output_dir)
        self._pipeline_runner = runner
        _ACTIVE_RUNNERS.append(runner)
        runner_ref = weakref.ref(runner)

        def _cleanup_runner() -> None:
            r = runner_ref()
            if r and r in _ACTIVE_RUNNERS:
                _ACTIVE_RUNNERS.remove(r)

        runner.finished_ok.connect(_cleanup_runner)
        runner.failed.connect(lambda _msg: _cleanup_runner())
        runner.failed.connect(progress.mark_failed)
        runner.finished_ok.connect(progress.mark_finished)
        runner.failed.connect(self._on_pipeline_failed)
        runner.finished_ok.connect(self._on_pipeline_finished)
        runner.start()

        return
