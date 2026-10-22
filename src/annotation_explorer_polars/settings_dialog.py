"""Settings dialog reachable from the runs window's gear button.

Phase 1 exposes the global compute limits (max CPU threads, max container RAM)
persisted to bactolomew_config.json's `resources` block via ToolEnv.save_app_config,
the default population structure method (`phylo.population_method`), plus the installed
cgMLST schemas, with "Add schema…" (import local files, one or
many) and "Download from cgMLST.org…" (pipeline.CgmlstInstall). The Databases page
lists, adds and removes the marker and compound databases (pipeline.UserDatabases).
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

from PyQt6 import QtCore, QtWidgets

from pipeline import UserDatabases
from pipeline.Resources import detect_budget

from .settings_column import settings_combo



def _caption(text: str) -> QtWidgets.QLabel:
    """A field caption above its control (theme role settingsCaption)."""
    label = QtWidgets.QLabel(text)
    label.setProperty("role", "settingsCaption")
    return label


def _hint(text: str) -> QtWidgets.QLabel:
    label = QtWidgets.QLabel(text)
    label.setWordWrap(True)
    label.setProperty("role", "windowHint")
    return label


_DB_COLUMNS = ("Name", "Category", "Sequences", "Size", "Added", "Status")


class SettingsDialog(QtWidgets.QDialog):
    """Edit global settings. Reads/writes bactolomew_config.json (merge-based)."""

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None, *, db_root: Optional[Path] = None) -> None:
        super().__init__(parent)
        self._db_root = db_root  # None: the project's DB folder
        self.setWindowTitle("Settings")
        self.setObjectName("SettingsDialog")
        self.resize(560, 840)
        self._cpu_total = detect_budget().cores
        self._build_ui()
        self._apply_app_mode()
        self._load_current()

    # ---- UI ----
    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(22, 20, 22, 18)
        layout.setSpacing(12)
        # The window never gets smaller than its content (no overlapping controls).
        layout.setSizeConstraint(QtWidgets.QLayout.SizeConstraint.SetMinimumSize)

        title = QtWidgets.QLabel("Settings")
        title.setProperty("role", "windowTitle")
        layout.addWidget(title)

        self.pages = QtWidgets.QTabWidget()
        general = QtWidgets.QWidget()
        general_layout = QtWidgets.QVBoxLayout(general)
        general_layout.setContentsMargins(12, 12, 12, 12)
        general_layout.setSpacing(12)
        general_layout.addWidget(self._build_resources_group())
        general_layout.addWidget(self._build_population_group())
        general_layout.addWidget(self._build_schemas_group(), stretch=1)
        self.pages.addTab(general, "General")
        self.pages.addTab(self._build_databases_page(), "Databases")
        layout.addWidget(self.pages, stretch=1)

        footer = QtWidgets.QFrame()
        footer.setObjectName("StartFooter")
        btn_row = QtWidgets.QHBoxLayout(footer)
        btn_row.setContentsMargins(0, 12, 0, 0)
        btn_row.setSpacing(8)
        self.pipeline_btn = QtWidgets.QPushButton("Set Up Pipeline…")
        self.pipeline_btn.setObjectName("SetUpPipelineButton")
        self.pipeline_btn.setToolTip("Download the tools image and databases (over 40 GB) to run predictions")
        self.pipeline_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.pipeline_btn.clicked.connect(self._set_up_pipeline)
        btn_row.addWidget(self.pipeline_btn)
        btn_row.addStretch(1)
        cancel_btn = QtWidgets.QPushButton("Cancel")
        cancel_btn.setProperty("role", "flat")
        cancel_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        cancel_btn.clicked.connect(self.reject)
        btn_row.addWidget(cancel_btn)
        save_btn = QtWidgets.QPushButton("Save")
        save_btn.setProperty("role", "primary")
        save_btn.setDefault(True)
        save_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        save_btn.clicked.connect(self._save)
        btn_row.addWidget(save_btn)
        layout.addWidget(footer)

    def _build_resources_group(self) -> QtWidgets.QGroupBox:
        group = QtWidgets.QGroupBox("Resources")
        # Each caption above its control, the wrapped hint below them.
        box = QtWidgets.QVBoxLayout(group)
        box.setContentsMargins(14, 22, 14, 14)
        box.setSpacing(4)

        # Max CPU threads
        self.threads_spin = QtWidgets.QSpinBox()
        self.threads_spin.setRange(8, max(8, self._cpu_total))
        self.threads_spin.setSuffix(" threads")
        box.addWidget(_caption("Max CPU threads"))
        threads_row = QtWidgets.QHBoxLayout()
        threads_row.addWidget(self.threads_spin)
        threads_row.addStretch(1)
        box.addLayout(threads_row)
        box.addSpacing(8)

        # Max container RAM (opt-in; 0 in config == no limit)
        ram_row = QtWidgets.QHBoxLayout()
        self.ram_check = QtWidgets.QCheckBox("Limit to")
        self.ram_spin = QtWidgets.QSpinBox()
        self.ram_spin.setRange(16, 4096)
        self.ram_spin.setSuffix(" GB")
        self.ram_check.toggled.connect(self.ram_spin.setEnabled)
        ram_row.addWidget(self.ram_check)
        ram_row.addWidget(self.ram_spin)
        ram_row.addStretch(1)
        box.addWidget(_caption("Max RAM"))
        box.addLayout(ram_row)
        box.addSpacing(8)

        box.addWidget(
            _hint(
                f"Physical cores available: {self._cpu_total}. The thread cap applies to every new "
                "run (a per-run --cores still overrides). The RAM limit is applied to Docker tool "
                "containers (--memory) and enables CheckM2 low-memory mode when tight."
            )
        )
        return group

    def _build_population_group(self) -> QtWidgets.QGroupBox:
        group = QtWidgets.QGroupBox("Population structure")
        box = QtWidgets.QVBoxLayout(group)
        box.setContentsMargins(14, 22, 14, 14)
        box.setSpacing(4)
        self._population_method = "fastbaps"
        options = [
            ("fastbaps", "FastBAPS", "FastBAPS clusters: lineage = Level 1, cluster = Level 2."),
            ("fastgear", "fastGEAR", "fastGEAR lineages and clusters (x86_64 tools images only)."),
            ("auto", "Automatic", "FastBAPS when the tools image has it, else fastGEAR."),
        ]
        self.population_combo = settings_combo(group, options, self._population_method, self._set_population_method)
        self.population_combo.setObjectName("DefaultPopulationMethodCombo")
        box.addWidget(_caption("Default method"))
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.population_combo)
        row.addStretch(1)
        box.addLayout(row)
        box.addSpacing(8)
        box.addWidget(
            _hint(
                "A new run offers the choice when the tools image has both; a resumed run keeps its method."
            )
        )
        return group

    def _set_population_method(self, method: str) -> None:
        self._population_method = str(method)

    def _build_schemas_group(self) -> QtWidgets.QGroupBox:
        group = QtWidgets.QGroupBox("cgMLST schemas")
        v = QtWidgets.QVBoxLayout(group)
        v.setContentsMargins(14, 22, 14, 14)
        v.setSpacing(8)

        self.schema_list = QtWidgets.QListWidget()
        self.schema_list.setMinimumHeight(132)  # about six schemas without scrolling
        v.addWidget(self.schema_list, stretch=1)

        row = QtWidgets.QHBoxLayout()
        add_btn = QtWidgets.QPushButton("Add schema…")
        self.schema_add_btn = add_btn
        add_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        add_btn.clicked.connect(self._add_schema)
        row.addWidget(add_btn)
        download_btn = QtWidgets.QPushButton("Download from cgMLST.org…")
        self.schema_download_btn = download_btn
        download_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        download_btn.clicked.connect(self._download_schemas)
        row.addWidget(download_btn)
        row.addStretch(1)
        v.addLayout(row)

        v.addWidget(
            _hint(
                "Installed species schemas under DB/cgmlst. A run auto-selects the matching schema "
                "from its MLST result."
            )
        )
        return group

    def _build_databases_page(self) -> QtWidgets.QWidget:
        page = QtWidgets.QWidget()
        box = QtWidgets.QVBoxLayout(page)
        box.setContentsMargins(14, 14, 14, 14)
        box.setSpacing(8)
        box.addWidget(_caption("Marker and compound databases"))

        self.db_table = QtWidgets.QTableWidget(0, len(_DB_COLUMNS))
        self.db_table.setObjectName("UserDatabaseTable")
        self.db_table.setHorizontalHeaderLabels(list(_DB_COLUMNS))
        self.db_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.db_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.db_table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.db_table.setAlternatingRowColors(True)
        self.db_table.setWordWrap(False)
        self.db_table.verticalHeader().setVisible(False)
        self.db_table.setMinimumHeight(180)
        header = self.db_table.horizontalHeader()
        header.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.db_table.setTextElideMode(QtCore.Qt.TextElideMode.ElideRight)
        self.db_table.itemSelectionChanged.connect(self._sync_db_buttons)
        box.addWidget(self.db_table, stretch=1)

        self.db_empty = _hint("")
        box.addWidget(self.db_empty)
        self.db_empty_space = QtWidgets.QWidget()
        box.addWidget(self.db_empty_space, stretch=1)

        row = QtWidgets.QHBoxLayout()
        add_btn = QtWidgets.QPushButton("Add database…")
        self.db_add_btn = add_btn
        add_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        add_btn.clicked.connect(self._add_database)
        row.addWidget(add_btn)
        self.db_remove_btn = QtWidgets.QPushButton("Remove")
        self.db_remove_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.db_remove_btn.clicked.connect(self._remove_database)
        row.addWidget(self.db_remove_btn)
        row.addStretch(1)
        box.addLayout(row)

        box.addWidget(
            _hint(
                "Protein FASTAs searched with DIAMOND after Bakta, from the next run. You are "
                "responsible for each database's license."
            )
        )
        return page

    # ---- data ----
    def _load_current(self) -> None:
        from pipeline.ToolEnv import configured_resources

        res = configured_resources()
        self.threads_spin.setValue(max(8, min(res["max_threads"], max(8, self._cpu_total))))
        ram = res["max_ram_gb"]
        self.ram_check.setChecked(ram > 0)
        self.ram_spin.setEnabled(ram > 0)
        self.ram_spin.setValue(ram if ram >= 16 else 16)
        try:
            from pipeline.ToolEnv import configured_phylo_settings

            method = str(configured_phylo_settings().get("population_method") or "fastbaps")
        except ValueError:
            method = "fastbaps"
        index = self.population_combo.findData(method)
        self.population_combo.setCurrentIndex(index if index >= 0 else 0)
        # Saved only when changed: an unset method stays the default, which an image without
        # FastBAPS replaces by fastGEAR instead of refusing the run.
        self._population_loaded = self._population_method

        self._reload_schema_list()
        self._reload_databases()

    def _reload_databases(self) -> None:
        try:
            databases = UserDatabases.list_databases(self._db_root)
            error = ""
        except UserDatabases.DatabaseError as exc:
            databases, error = [], str(exc)
        table = self.db_table
        table.setRowCount(len(databases))
        for r, db in enumerate(databases):
            values = (
                db.name,
                UserDatabases.CATEGORY_LABELS.get(db.category, db.category),
                "" if db.sequences is None else f"{db.sequences:,}",
                UserDatabases.human_size(db.size) if db.size else "",
                db.added[:10],
                db.status,
            )
            for c, value in enumerate(values):
                item = QtWidgets.QTableWidgetItem(value)
                if c in (2, 3):
                    item.setTextAlignment(
                        QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter
                    )
                item.setToolTip("\n".join(t for t in (str(db.path), db.note) if t))
                table.setItem(r, c, item)
        table.setVisible(bool(databases))
        self.db_empty.setText(
            error or "No databases added. Runs annotate with Bakta and the pipeline's own tools."
        )
        self.db_empty.setVisible(not databases)
        self.db_empty_space.setVisible(not databases)
        self._sync_db_buttons()

    def _sync_db_buttons(self) -> None:
        self.db_remove_btn.setEnabled(bool(self.db_table.selectedItems()))

    def _reload_schema_list(self) -> None:
        self.schema_list.clear()
        for name in self._installed_schema_names():
            self.schema_list.addItem(name)
        if self.schema_list.count() == 0:
            item = QtWidgets.QListWidgetItem("No schemas installed yet.")
            item.setFlags(QtCore.Qt.ItemFlag.NoItemFlags)
            self.schema_list.addItem(item)

    def _installed_schema_names(self) -> list[str]:
        from pipeline.ToolEnv import configured_database_dirs

        cgmlst_dir = configured_database_dirs().get("cgmlst_dir") or ""
        if not cgmlst_dir:
            return []
        root = Path(cgmlst_dir)
        if not root.is_dir():
            return []
        names: list[str] = []
        for child in sorted(root.iterdir(), key=lambda p: p.name.lower()):
            if child.is_dir() and (child / ".genes_list").is_file():
                names.append(child.name)
        return names

    def _apply_app_mode(self) -> None:
        """A viewer install has no tools image: what builds schemas or databases is off,
        and Set Up Pipeline… offers the full setup."""
        from pipeline.ToolEnv import viewer_mode

        viewer = viewer_mode()
        self.pipeline_btn.setVisible(viewer)
        for button in (self.schema_add_btn, self.schema_download_btn, self.db_add_btn):
            button.setEnabled(not viewer)
            button.setToolTip("Needs the pipeline setup" if viewer else "")

    def _set_up_pipeline(self) -> None:
        from .first_time_setup import offer_pipeline_setup

        offer_pipeline_setup(self)
        self._apply_app_mode()

    # ---- actions ----
    def _add_schema(self) -> None:
        from .schema_import_dialog import AddSchemaDialog

        dlg = AddSchemaDialog(self)
        if dlg.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            self._reload_schema_list()

    def _download_schemas(self) -> None:
        from pipeline.CgmlstInstall import default_dest

        from .cgmlst_download_dialog import CgmlstDownloadDialog

        dlg = CgmlstDownloadDialog(default_dest(), self)
        dlg.exec()
        if dlg.installed_any:
            self._reload_schema_list()

    def _add_database(self) -> None:
        dlg = AddDatabaseDialog(self, db_root=self._db_root)
        if dlg.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            self._reload_databases()

    def _remove_database(self) -> None:
        row = self.db_table.currentRow()
        item = self.db_table.item(row, 0) if row >= 0 else None
        if item is None:
            return
        name = item.text()
        answer = QtWidgets.QMessageBox.question(
            self,
            "Remove database",
            f"Remove {name}? Runs stop searching it; its file is moved to DB/_removed.",
        )
        if answer != QtWidgets.QMessageBox.StandardButton.Yes:
            return
        try:
            UserDatabases.remove_database(name, self._db_root)
        except Exception as exc:  # noqa: BLE001 - surfaced to the user
            QtWidgets.QMessageBox.critical(self, "Remove failed", str(exc))
        self._reload_databases()

    def _save(self) -> None:
        from pipeline.ToolEnv import save_app_config

        resources = {
            "max_threads": int(self.threads_spin.value()),
            "max_ram_gb": int(self.ram_spin.value()) if self.ram_check.isChecked() else 0,
        }
        updates: dict = {"resources": resources}
        if self._population_method != getattr(self, "_population_loaded", None):
            updates["phylo"] = {"population_method": self._population_method}
        try:
            save_app_config(updates)
        except Exception as exc:  # pragma: no cover - surfaced to the user
            QtWidgets.QMessageBox.critical(self, "Save failed", f"Could not save settings:\n{exc}")
            return
        self.accept()


class AddDatabaseDialog(QtWidgets.QDialog):
    """Add a protein FASTA as a marker or compound database (pipeline.UserDatabases)."""

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None, *, db_root: Optional[Path] = None) -> None:
        super().__init__(parent)
        self._db_root = db_root
        self.setWindowTitle("Add database")
        self.setObjectName("AddDatabaseDialog")
        self.resize(520, 0)
        self.added: Optional[UserDatabases.UserDatabase] = None
        self._worker: Optional[_AddDatabaseWorker] = None
        self._category = UserDatabases.MARKER_CATEGORIES[0]

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(22, 20, 22, 18)
        layout.setSpacing(10)
        layout.setSizeConstraint(QtWidgets.QLayout.SizeConstraint.SetMinimumSize)
        title = QtWidgets.QLabel("Add database")
        title.setProperty("role", "windowTitle")
        layout.addWidget(title)

        file_row = QtWidgets.QHBoxLayout()
        self.file_edit = QtWidgets.QLineEdit()
        self.file_edit.setPlaceholderText("Protein FASTA (amino acids)")
        browse = QtWidgets.QPushButton("Browse…")
        browse.clicked.connect(self._browse)
        file_row.addWidget(self.file_edit, 1)
        file_row.addWidget(browse)
        layout.addWidget(_caption("FASTA file"))
        layout.addLayout(file_row)

        self.name_edit = QtWidgets.QLineEdit()
        self.name_edit.setPlaceholderText("letters, digits, - and _")
        layout.addWidget(_caption("Name"))
        layout.addWidget(self.name_edit)

        options = [(c, UserDatabases.CATEGORY_LABELS[c], "") for c in UserDatabases.CATEGORIES]
        self.category_combo = settings_combo(self, options, self._category, self._set_category)
        layout.addWidget(_caption("Category"))
        layout.addWidget(self.category_combo)

        self.note_edit = QtWidgets.QPlainTextEdit()
        self.note_edit.setPlaceholderText("Source, version or license (optional)")
        self.note_edit.setFixedHeight(64)
        layout.addWidget(_caption("Note"))
        layout.addWidget(self.note_edit)

        self.error_label = _hint("")
        self.error_label.setProperty("role", "warningChip")
        self.error_label.hide()
        layout.addWidget(self.error_label)

        footer = QtWidgets.QFrame()
        footer.setObjectName("StartFooter")
        btns = QtWidgets.QHBoxLayout(footer)
        btns.setContentsMargins(0, 12, 0, 0)
        btns.setSpacing(8)
        btns.addStretch(1)
        self.cancel_btn = QtWidgets.QPushButton("Cancel")
        self.cancel_btn.setProperty("role", "flat")
        self.cancel_btn.clicked.connect(self.reject)
        btns.addWidget(self.cancel_btn)
        self.add_btn = QtWidgets.QPushButton("Add")
        self.add_btn.setProperty("role", "primary")
        self.add_btn.setDefault(True)
        self.add_btn.clicked.connect(self._add)
        btns.addWidget(self.add_btn)
        layout.addWidget(footer)

    def _set_category(self, value: str) -> None:
        self._category = value

    def _browse(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Protein FASTA", "", "FASTA (*.fasta *.fa *.faa *.fas *.gz);;All files (*)"
        )
        if path:
            self.set_file(path)

    def set_file(self, path: str) -> None:
        self.file_edit.setText(path)
        if not self.name_edit.text().strip():
            self.name_edit.setText(re.sub(r"[^A-Za-z0-9_-]+", "_", Path(path).stem).strip("_-")[:64])

    def _add(self) -> None:
        # Checking, copying and hashing a large FASTA takes minutes: run it off the GUI thread.
        self.error_label.hide()
        self._set_busy(True)
        self._worker = _AddDatabaseWorker(
            Path(self.file_edit.text().strip()),
            self.name_edit.text(),
            self._category,
            self.note_edit.toPlainText(),
            self._db_root,
        )
        self._worker.done.connect(self._add_done)
        self._worker.failed.connect(self._add_failed)
        self._worker.start()

    def _set_busy(self, busy: bool) -> None:
        for widget in (self.file_edit, self.name_edit, self.category_combo, self.note_edit, self.add_btn, self.cancel_btn):
            widget.setEnabled(not busy)
        self.add_btn.setText("Adding…" if busy else "Add")

    def _add_done(self, db: UserDatabases.UserDatabase) -> None:
        self._worker.wait()
        self._worker = None
        self.added = db
        self.accept()

    def _add_failed(self, message: str) -> None:
        self._worker.wait()
        self._worker = None
        self._set_busy(False)
        self.error_label.setText(message)
        self.error_label.show()

    def reject(self) -> None:
        if self._worker is None:  # the copy cannot be stopped halfway
            super().reject()


class _AddDatabaseWorker(QtCore.QThread):
    """Runs UserDatabases.add_database off the GUI thread."""

    done = QtCore.pyqtSignal(object)
    failed = QtCore.pyqtSignal(str)

    def __init__(self, fasta: Path, name: str, category: str, note: str, db_root: Optional[Path]) -> None:
        super().__init__()
        self._args = (fasta, name, category)
        self._note = note
        self._db_root = db_root

    def run(self) -> None:
        try:
            db = UserDatabases.add_database(*self._args, note=self._note, db_root=self._db_root)
        except (UserDatabases.DatabaseError, OSError) as exc:
            self.failed.emit(str(exc))
        except Exception as exc:  # noqa: BLE001 - surfaced to the user
            self.failed.emit(f"Could not add the database: {exc}")
        else:
            self.done.emit(db)
