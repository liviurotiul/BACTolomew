"""The run list's sharing dialogs: Import Run..., Export for Sharing... and Run
Details. The work (pipeline.RunShare) runs on a worker thread; the dialogs only
show it. Same look as the run list: captions above, white cards, theme roles."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any, Callable, Mapping, Optional, Sequence

from PyQt6 import QtCore, QtWidgets

from pipeline import RunShare
from pipeline.RunShare import Issue, format_size

# A built annotation table takes about 3 MB per genome (448 MB for 153 genomes).
_TABLE_BYTES_PER_GENOME = 3_000_000


def compatibility_context() -> dict[str, Any]:
    """This install's versions for RunShare.check_compatibility (the table format
    is the explorer's)."""
    from .cache import TABLE_CACHE_SCHEMA_VERSION

    context = RunShare.current_versions()
    context["table_schema"] = TABLE_CACHE_SCHEMA_VERSION
    return context


def run_issues(source: Any, context: Optional[Mapping[str, Any]] = None) -> list[Issue]:
    """The compatibility findings of a run folder or of a package's manifest."""
    context = dict(context or compatibility_context())
    return RunShare.check_compatibility(source, **context)


def _caption(text: str) -> QtWidgets.QLabel:
    label = QtWidgets.QLabel(text)
    label.setProperty("role", "settingsCaption")
    return label


def _card(rows: Sequence[tuple[str, str]]) -> QtWidgets.QFrame:
    """A white card of label / value rows."""
    card = QtWidgets.QFrame()
    card.setObjectName("RunDetailsSection")
    grid = QtWidgets.QGridLayout(card)
    grid.setContentsMargins(14, 10, 14, 10)
    grid.setHorizontalSpacing(14)
    grid.setVerticalSpacing(5)
    for row, (label, value) in enumerate(rows):
        name = QtWidgets.QLabel(label)
        name.setProperty("role", "runMeta")
        name.setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignTop)
        text = QtWidgets.QLabel(value)
        text.setProperty("role", "runStep")
        text.setWordWrap(True)
        text.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        grid.addWidget(name, row, 0)
        grid.addWidget(text, row, 1)
    grid.setColumnStretch(1, 1)
    return card


def _issue_rows(issues: Sequence[Issue]) -> list[tuple[str, str]]:
    return [("Warning" if issue.level == "warning" else "Note", issue.text) for issue in issues]


def _issues_card(issues: Sequence[Issue]) -> QtWidgets.QFrame:
    card = _card(_issue_rows(issues))
    grid = card.layout()
    for row, issue in enumerate(issues):
        if issue.level == "warning":
            item = grid.itemAtPosition(row, 0)
            if item is not None and item.widget() is not None:
                item.widget().setProperty("role", "runWarning")
    return card


def _dialog_header(layout: QtWidgets.QVBoxLayout, title: str, hint: str) -> QtWidgets.QLabel:
    heading = QtWidgets.QLabel(title)
    heading.setProperty("role", "windowTitle")
    layout.addWidget(heading)
    sub = QtWidgets.QLabel(hint)
    sub.setWordWrap(True)
    sub.setProperty("role", "windowHint")
    layout.addWidget(sub)
    layout.addSpacing(8)
    return sub


def _footer(layout: QtWidgets.QVBoxLayout) -> QtWidgets.QHBoxLayout:
    footer = QtWidgets.QFrame()
    footer.setObjectName("StartFooter")
    row = QtWidgets.QHBoxLayout(footer)
    row.setContentsMargins(0, 12, 0, 0)
    row.setSpacing(8)
    layout.addWidget(footer)
    return row


def _button(text: str, role: Optional[str] = None) -> QtWidgets.QPushButton:
    button = QtWidgets.QPushButton(text)
    button.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
    if role:
        button.setProperty("role", role)
    return button


def _repolish(widget: QtWidgets.QWidget) -> None:
    widget.style().unpolish(widget)
    widget.style().polish(widget)


class _Worker(QtCore.QThread):
    """Runs ``job(progress, cancelled)`` off the GUI thread."""

    progressed = QtCore.pyqtSignal(int, int, str)
    done = QtCore.pyqtSignal(object)
    failed = QtCore.pyqtSignal(str)

    def __init__(self, job: Callable[..., Any], parent: Optional[QtCore.QObject] = None) -> None:
        super().__init__(parent)
        self._job = job
        self._cancel = False

    def cancel(self) -> None:
        self._cancel = True

    def run(self) -> None:  # the worker thread
        try:
            result = self._job(self.progressed.emit, lambda: self._cancel)
        except RunShare.ShareCancelled:
            self.failed.emit("")
            return
        except Exception as exc:
            self.failed.emit(str(exc) or exc.__class__.__name__)
            return
        self.done.emit(result)


class _WorkDialog(QtWidgets.QDialog):
    """The shared frame of Import and Export: a progress line, Cancel stops the
    work (or closes), the action button starts it."""

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setObjectName("RunShareDialog")
        self._worker: Optional[_Worker] = None
        self.layout_ = QtWidgets.QVBoxLayout(self)
        self.layout_.setContentsMargins(22, 20, 22, 14)
        self.layout_.setSpacing(4)

    def _add_progress(self) -> None:
        self.progress = QtWidgets.QProgressBar()
        self.progress.setProperty("role", "stepBar")
        self.progress.setRange(0, 1000)
        self.progress.setTextVisible(False)
        self.progress.hide()
        self.progress_text = QtWidgets.QLabel("")
        self.progress_text.setProperty("role", "runMeta")
        self.progress_text.hide()
        self.layout_.addSpacing(6)
        self.layout_.addWidget(self.progress)
        self.layout_.addWidget(self.progress_text)

    def _start(self, job: Callable[..., Any], *, busy_text: str) -> None:
        self.progress.setValue(0)
        self.progress.setProperty("state", "running")
        _repolish(self.progress)
        self.progress.show()
        self.progress_text.setText(busy_text)
        self.progress_text.show()
        self.action_button.setEnabled(False)
        self._set_inputs_enabled(False)
        worker = _Worker(job, self)
        self._worker = worker
        worker.progressed.connect(self._on_progress)
        worker.done.connect(self._on_done)
        worker.failed.connect(self._on_failed)
        worker.start()

    def _set_inputs_enabled(self, enabled: bool) -> None:
        pass

    def _on_progress(self, done: int, total: int, text: str) -> None:
        if total > 0:
            self.progress.setValue(int(1000 * min(done, total) / total))
            self.progress_text.setText(f"{format_size(done)} of {format_size(total)}")

    def _on_done(self, result: object) -> None:
        raise NotImplementedError

    def _on_failed(self, message: str) -> None:
        self._worker = None
        self.action_button.setEnabled(True)
        self._set_inputs_enabled(True)
        self.progress.hide()
        self.progress_text.hide()
        if message:
            QtWidgets.QMessageBox.critical(self, self.windowTitle(), message)

    def reject(self) -> None:  # type: ignore[override]
        if self._worker is not None and self._worker.isRunning():
            self._worker.cancel()
            self.progress_text.setText("Stopping…")
            self._worker.wait(30000)
            return
        super().reject()

    def closeEvent(self, event) -> None:  # type: ignore[override]
        if self._worker is not None and self._worker.isRunning():
            self._worker.cancel()
            self._worker.wait(30000)
        super().closeEvent(event)


# --------------------------------------------------------------------- import


class ImportRunDialog(_WorkDialog):
    """Import Run...: a share package (.zip) or a run folder; shows what it is and
    what may not match, and imports it anyway (nothing found blocks it)."""

    def __init__(self, runs_root: Path, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Import Run")
        self.resize(620, 520)
        self._runs_root = Path(runs_root)
        self._source: Optional[Path] = None
        self.result: Optional[RunShare.ImportResult] = None
        self.issues: list[Issue] = []
        layout = self.layout_
        _dialog_header(
            layout,
            "Import a shared run",
            "Open a run someone shared with you (.zip), or copy a run folder into the run list. "
            "Viewing it needs no pipeline setup.",
        )
        layout.addWidget(_caption("Package or run folder"))
        row = QtWidgets.QHBoxLayout()
        row.setSpacing(6)
        self.path_edit = QtWidgets.QLineEdit()
        self.path_edit.setReadOnly(True)
        self.path_edit.setPlaceholderText("No package chosen")
        row.addWidget(self.path_edit, 1)
        self.file_button = _button("Choose File...")
        self.file_button.clicked.connect(self._choose_file)
        row.addWidget(self.file_button)
        self.folder_button = _button("Choose Folder...")
        self.folder_button.clicked.connect(self._choose_folder)
        row.addWidget(self.folder_button)
        layout.addLayout(row)
        layout.addSpacing(10)

        self.body = QtWidgets.QVBoxLayout()
        self.body.setSpacing(4)
        layout.addLayout(self.body)
        layout.addStretch(1)
        self._add_progress()

        buttons = _footer(layout)
        buttons.addStretch(1)
        cancel = _button("Cancel", "flat")
        cancel.clicked.connect(self.reject)
        buttons.addWidget(cancel)
        self.action_button = _button("Import", "primary")
        self.action_button.setEnabled(False)
        self.action_button.clicked.connect(self._import)
        buttons.addWidget(self.action_button)

    def _set_inputs_enabled(self, enabled: bool) -> None:
        self.file_button.setEnabled(enabled)
        self.folder_button.setEnabled(enabled)

    def _choose_file(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Choose a shared run", str(Path.home()), "Shared runs (*.zip);;All files (*)"
        )
        if path:
            self.set_source(Path(path))

    def _choose_folder(self) -> None:
        path = QtWidgets.QFileDialog.getExistingDirectory(self, "Choose a run folder", str(Path.home()))
        if path:
            self.set_source(Path(path))

    def _clear_body(self) -> None:
        while self.body.count():
            item = self.body.takeAt(0)
            if item.widget() is not None:
                item.widget().deleteLater()

    def set_source(self, path: Path) -> None:
        """Read the package (or folder) and show its summary and findings."""
        self._clear_body()
        self._source = None
        self.issues = []
        self.path_edit.setText(str(path))
        summary, issues, error = self._describe(path)
        if error:
            note = QtWidgets.QLabel(error)
            note.setProperty("role", "runWarning")
            note.setWordWrap(True)
            self.body.addWidget(note)
            self.action_button.setEnabled(False)
            self.action_button.setText("Import")
            return
        self._source = path
        self.issues = issues
        self.body.addWidget(_caption("Run"))
        self.body.addWidget(_card(summary))
        if issues:
            self.body.addSpacing(8)
            self.body.addWidget(_caption("Before you import"))
            self.body.addWidget(_issues_card(issues))
        warnings = [issue for issue in issues if issue.level == "warning"]
        self.action_button.setText("Import Anyway" if warnings else "Import")
        self.action_button.setEnabled(True)

    def _describe(self, path: Path) -> tuple[list[tuple[str, str]], list[Issue], str]:
        if path.is_dir():
            config_path = path / "run_config.json"
            if not config_path.is_file():
                return [], [], "This folder is not a BACTolomew run: run_config.json is missing."
            info = RunShare.run_info(path)
            rows = [row for row in info.get("Run", []) if row[0] in {"Run", "Sample set", "Genomes", "Status"}]
            rows += [row for row in info.get("Made with", []) if row[0] in {"BACTolomew", "Tools image"}]
            size = sum(p.stat().st_size for p in path.rglob("*") if p.is_file() and not p.is_symlink())
            rows.append(("Size", format_size(size)))
            return rows, run_issues(path), ""
        manifest = RunShare.read_manifest(path)
        if not manifest:
            return [], [], "This file is not a shared BACTolomew run."
        files = [f for f in manifest.get("files") or [] if isinstance(f, Mapping)]
        size = sum(int(f.get("size") or 0) for f in files)
        rows = [
            ("Run", str(manifest.get("run_id") or "")),
            ("Sample set", str(manifest.get("sample_set_id") or "")),
            ("Genomes", str(manifest.get("sample_count") or "")),
            ("Made with", f"BACTolomew {manifest.get('made_with')}" if manifest.get("made_with") else "not recorded"),
            ("Exported", f"{manifest.get('exported_at') or ''} (BACTolomew {manifest.get('app_version') or '?'})"),
            ("Assemblies", "included" if manifest.get("includes_assemblies") else "not included"),
            ("Metadata", "CSV included" if manifest.get("metadata") == "csv" else "run data only"),
            ("Size", f"{format_size(size)} unpacked"),
        ]
        free = shutil.disk_usage(self._runs_root if self._runs_root.exists() else Path.home()).free
        if size > free:
            return rows, [], f"Not enough disk space: {format_size(size)} needed, {format_size(free)} free."
        return [row for row in rows if row[1]], run_issues(manifest), ""

    def _import(self) -> None:
        source = self._source
        if source is None:
            return
        runs_root = self._runs_root
        if RunShare.run_id_taken(runs_root, RunShare.import_run_id(source)) and not confirm_import_copy(self):
            return

        def job(progress: Callable[..., None], cancelled: Callable[[], bool]) -> RunShare.ImportResult:
            if source.is_dir():
                return RunShare.import_folder(source, runs_root, progress=progress, cancelled=cancelled)
            return RunShare.import_archive(source, runs_root, progress=progress, cancelled=cancelled)

        self._start(job, busy_text="Importing…")

    def _on_done(self, result: object) -> None:
        self._worker = None
        self.result = result if isinstance(result, RunShare.ImportResult) else None
        self.accept()


# --------------------------------------------------------------------- export


def _tables_dir(run_id: str, run_dir: Path) -> tuple[Optional[Path], bool]:
    """(folder with the run's tables, they are up to date) for an export."""
    from .cache import cache_root, cached_tables_exist, run_cache_slug

    if RunShare.is_shared_run(run_dir):
        return run_dir / RunShare.TABLES_DIR, True
    if not (cache_root() / run_cache_slug(run_id)).is_dir():
        return None, False
    if cached_tables_exist(run_id, run_dir):
        return cache_root() / run_cache_slug(run_id), True
    return None, False


def _metadata_csv(run_id: str) -> Optional[Path]:
    from .cache import load_run_cache

    cached = load_run_cache(run_id, create=False)
    path = Path(cached["metadata_csv"]).expanduser() if cached.get("metadata_csv") else None
    return path if path is not None and path.is_file() else None


class ExportRunDialog(_WorkDialog):
    """Export for Sharing...: one package with what the explorer needs to show the
    run (its tables built first when they are not yet)."""

    def __init__(
        self,
        run_id: str,
        run_dir: Path,
        *,
        title: str = "",
        sample_count: Optional[int] = None,
        app_version: str = "",
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Export for Sharing")
        self.resize(600, 460)
        self._run_id = run_id
        self._run_dir = Path(run_dir)
        self._app_version = app_version
        self.result_path: Optional[Path] = None
        self._tables_dir, self._tables_ready = _tables_dir(run_id, self._run_dir)
        sizes = RunShare.estimate_share_size(self._run_dir, self._tables_dir)
        if not self._tables_ready:
            sizes["tables"] = int(sample_count or 0) * _TABLE_BYTES_PER_GENOME
        self._sizes = sizes
        self._csv = _metadata_csv(run_id)
        layout = self.layout_
        _dialog_header(
            layout,
            "Export for sharing",
            f"{title or run_id}: one file with what the explorer needs to show this run. "
            "Others open it with Import Run..., without the pipeline.",
        )
        layout.addWidget(_caption("Save as"))
        row = QtWidgets.QHBoxLayout()
        row.setSpacing(6)
        self.path_edit = QtWidgets.QLineEdit(str(Path.home() / f"{run_id}{RunShare.ARCHIVE_SUFFIX}"))
        row.addWidget(self.path_edit, 1)
        self.choose_button = _button("Choose...")
        self.choose_button.clicked.connect(self._choose_dest)
        row.addWidget(self.choose_button)
        layout.addLayout(row)
        layout.addSpacing(10)

        layout.addWidget(_caption("Include"))
        self.assemblies_box = QtWidgets.QCheckBox(
            f"Include assemblies (FASTA), {format_size(sizes['assemblies'])}"
        )
        self.assemblies_box.setEnabled(sizes["assemblies"] > 0)
        self.assemblies_box.setChecked(sizes["assemblies"] > 0)
        self.assemblies_box.toggled.connect(self._update_total)
        layout.addWidget(self.assemblies_box)
        self.metadata_box = QtWidgets.QCheckBox(
            f"Include metadata CSV ({self._csv.name})" if self._csv else "Include metadata CSV"
        )
        self.metadata_box.setEnabled(self._csv is not None)
        self.metadata_box.setChecked(self._csv is not None)
        self.metadata_box.toggled.connect(self._update_total)
        layout.addWidget(self.metadata_box)
        note = QtWidgets.QLabel(
            "The CSV may hold private sample data. Without it, the run opens with its own data."
            if self._csv
            else "Without a metadata CSV, the run opens with its own data."
        )
        note.setProperty("role", "windowHint")
        note.setWordWrap(True)
        layout.addWidget(note)
        layout.addSpacing(10)

        layout.addWidget(_caption("Size"))
        self.total_label = QtWidgets.QLabel("")
        self.total_label.setProperty("role", "runStep")
        self.total_label.setWordWrap(True)
        layout.addWidget(self.total_label)
        self._update_total()
        layout.addStretch(1)
        self._add_progress()

        buttons = _footer(layout)
        buttons.addStretch(1)
        cancel = _button("Cancel", "flat")
        cancel.clicked.connect(self.reject)
        buttons.addWidget(cancel)
        self.action_button = _button("Export", "primary")
        self.action_button.clicked.connect(self._export)
        buttons.addWidget(self.action_button)

    def _set_inputs_enabled(self, enabled: bool) -> None:
        for widget in (self.path_edit, self.choose_button, self.assemblies_box, self.metadata_box):
            widget.setEnabled(enabled)
        if enabled:
            self.assemblies_box.setEnabled(self._sizes["assemblies"] > 0)
            self.metadata_box.setEnabled(self._csv is not None)

    def _choose_dest(self) -> None:
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export for sharing", self.path_edit.text(), "Shared runs (*.zip)"
        )
        if path:
            if not path.endswith(".zip"):
                path += RunShare.ARCHIVE_SUFFIX
            self.path_edit.setText(path)

    def total_bytes(self) -> int:
        total = self._sizes["run_files"] + self._sizes["tables"]
        if self.assemblies_box.isChecked():
            total += self._sizes["assemblies"]
        if self.metadata_box.isChecked() and self._csv is not None:
            total += self._csv.stat().st_size
        return total

    def _update_total(self) -> None:
        text = f"About {format_size(self.total_bytes())}"
        if not self._tables_ready:
            text += ". The data tables are built first, from the run's Bakta results (a few minutes)."
        self.total_label.setText(text)

    def _export(self) -> None:
        dest = Path(self.path_edit.text().strip()).expanduser()
        if not dest.name:
            return
        if dest.exists():
            answer = QtWidgets.QMessageBox.question(self, "Export for Sharing", f"{dest.name} exists. Replace it?")
            if answer != QtWidgets.QMessageBox.StandardButton.Yes:
                return
        run_id, run_dir = self._run_id, self._run_dir
        include_assemblies = self.assemblies_box.isChecked()
        csv = self._csv if self.metadata_box.isChecked() else None
        tables_dir, ready = self._tables_dir, self._tables_ready
        app_version = self._app_version

        def job(progress: Callable[..., None], cancelled: Callable[[], bool]) -> Path:
            from pipeline.RunMetadata import mlst_scheme_species

            nonlocal tables_dir
            if not ready:
                tables_dir = _build_tables(run_id, run_dir, progress)
            if cancelled():
                raise RunShare.ShareCancelled()
            return RunShare.export_run(
                run_dir,
                tables_dir,
                dest,
                include_assemblies=include_assemblies,
                metadata_csv=csv,
                scheme_species=mlst_scheme_species(run_dir),
                app_version=app_version,
                progress=progress,
                cancelled=cancelled,
            )

        self._start(job, busy_text="Building data tables…" if not ready else "Exporting…")

    def _on_done(self, result: object) -> None:
        self._worker = None
        self.result_path = Path(result) if isinstance(result, (str, Path)) else None
        self.progress.setValue(1000)
        self.progress.setProperty("state", "success")
        _repolish(self.progress)
        size = format_size(self.result_path.stat().st_size) if self.result_path and self.result_path.exists() else ""
        self.progress_text.setText(f"Saved {self.result_path.name if self.result_path else ''} ({size})")
        self.action_button.setText("Done")
        self.action_button.setEnabled(True)
        self.action_button.clicked.disconnect()
        self.action_button.clicked.connect(self.accept)


def _build_tables(run_id: str, run_dir: Path, progress: Callable[..., None]) -> Path:
    """Build and cache a full run's tables (as opening it does), for an export."""
    from .cache import cache_run_dir, write_cached_tables
    from .data import load_output_data, validate_output_dir

    progress(0, 0, "Building data tables…")
    _cfg, errors = validate_output_dir(run_dir, deep=True)
    if errors:
        raise RunShare.ShareError("\n".join(errors))
    ann_df, ref_df, marker_hits_df = load_output_data(run_dir)
    write_cached_tables(run_id, run_dir, ann_df, ref_df, marker_hits_df)
    return cache_run_dir(run_id)


# --------------------------------------------------------------- run details


class RunDetailsDialog(QtWidgets.QDialog):
    """Run Details...: the run, what made it (BACTolomew, the tools image, the
    databases), the analysis, the findings of the compatibility check and, for
    an imported run, where it came from."""

    SECTIONS = ("Run", "Made with", "Analysis", "Warnings", "Shared")

    def __init__(
        self,
        run_dir: Path,
        issues: Sequence[Issue] = (),
        *,
        title: str = "",
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self.setObjectName("RunShareDialog")
        self.setWindowTitle("Run Details")
        self.resize(640, 640)
        info = RunShare.run_info(Path(run_dir))
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(22, 20, 22, 14)
        layout.setSpacing(4)
        _dialog_header(layout, "Run details", title or Path(run_dir).name)

        scroll = QtWidgets.QScrollArea()
        scroll.setObjectName("RunStepsScroll")
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        content = QtWidgets.QWidget()
        content.setObjectName("RunStepsContent")
        body = QtWidgets.QVBoxLayout(content)
        body.setContentsMargins(0, 0, 4, 0)
        body.setSpacing(4)
        self.section_titles: list[str] = []
        for section in self.SECTIONS:
            if section == "Warnings":
                if not issues:
                    continue
                card = _issues_card(list(issues))
            elif info.get(section):
                card = _card(info[section])
            else:
                continue
            self.section_titles.append(section)
            body.addWidget(_caption(section))
            body.addWidget(card)
            body.addSpacing(8)
        body.addStretch(1)
        scroll.setWidget(content)
        layout.addWidget(scroll, 1)

        buttons = _footer(layout)
        buttons.addStretch(1)
        close = _button("Close", "primary")
        close.clicked.connect(self.accept)
        buttons.addWidget(close)


def confirm_import_copy(parent: QtWidgets.QWidget) -> bool:
    """Before importing a run whose run id is already in the run list: "Import Copy"
    imports it as a separate run (its run id gets a _shared suffix)."""
    box = QtWidgets.QMessageBox(parent)
    box.setIcon(QtWidgets.QMessageBox.Icon.Warning)
    box.setWindowTitle("Import Run")
    box.setText("This run is already in your run list. Import it as a separate copy?")
    copy_button = box.addButton("Import Copy", QtWidgets.QMessageBox.ButtonRole.AcceptRole)
    box.addButton(QtWidgets.QMessageBox.StandardButton.Cancel)
    box.setDefaultButton(copy_button)
    box.exec()
    return box.clickedButton() is copy_button


def confirm_open_with_warnings(parent: QtWidgets.QWidget, title: str, issues: Sequence[Issue]) -> bool:
    """Before View on a run with warnings: they are listed, "Open anyway" opens."""
    warnings = [issue for issue in issues if issue.level == "warning"]
    if not warnings:
        return True
    box = QtWidgets.QMessageBox(parent)
    box.setIcon(QtWidgets.QMessageBox.Icon.Warning)
    box.setWindowTitle("Open run")
    box.setText(f"{title} may not show correctly in this version.")
    box.setInformativeText("\n".join(f"- {issue.text}" for issue in warnings))
    open_button = box.addButton("Open Anyway", QtWidgets.QMessageBox.ButtonRole.AcceptRole)
    box.addButton(QtWidgets.QMessageBox.StandardButton.Cancel)
    box.setDefaultButton(open_button)
    box.exec()
    return box.clickedButton() is open_button


__all__ = [
    "ExportRunDialog",
    "ImportRunDialog",
    "RunDetailsDialog",
    "compatibility_context",
    "confirm_import_copy",
    "confirm_open_with_warnings",
    "run_issues",
]
