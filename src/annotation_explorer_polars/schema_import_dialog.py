"""Dialog to import a new cgMLST schema into DB/cgmlst and register it.

Accepts a folder of per-locus allele FASTAs, a .zip of them, a single combined
FASTA, a single-locus FASTA, or an already-built chewBBACA schema. Conversion
runs on a background thread; on success the schema is registered in
scheme_species_catalog.json (species + MLST aliases) so it participates in the
existing MLST-based auto-selection.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

from PyQt6 import QtCore, QtWidgets

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
)
_HINT_STYLE = "color:#64748b;font-size:11px;"

_FASTA_EXTS = (".fasta", ".fa", ".fna", ".ffn", ".faa")


class _ImportWorker(QtCore.QThread):
    progress = QtCore.pyqtSignal(str)
    done = QtCore.pyqtSignal(object)
    failed = QtCore.pyqtSignal(str)

    def __init__(self, source: Path, target_dir: Path, single_fasta_mode: str, overwrite: bool) -> None:
        super().__init__()
        self._source = source
        self._target_dir = target_dir
        self._mode = single_fasta_mode
        self._overwrite = overwrite
        self._cancel = False

    def request_cancel(self) -> None:
        self._cancel = True

    def run(self) -> None:  # pragma: no cover - exercised via the GUI
        try:
            from pipeline.ChewbbacaSchemaBuild import build_schema

            def cb(msg: str) -> None:
                if self._cancel:
                    raise RuntimeError("cancelled")
                self.progress.emit(msg)

            result = build_schema(
                self._source,
                self._target_dir,
                single_fasta_mode=self._mode,
                overwrite=self._overwrite,
                progress=cb,
            )
            self.done.emit(result)
        except Exception as exc:  # noqa: BLE001 - surfaced to the user
            self.failed.emit(str(exc))


class AddSchemaDialog(QtWidgets.QDialog):
    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Add cgMLST schema")
        self.setObjectName("AddSchemaDialog")
        self.setStyleSheet("QDialog#AddSchemaDialog { background: #eef3f8; }")
        self.resize(560, 340)
        self.imported_schema_id: Optional[str] = None
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(22, 20, 22, 18)
        layout.setSpacing(10)

        title = QtWidgets.QLabel("Add cgMLST schema")
        title.setStyleSheet("font-size: 16px; font-weight: 800; color: #172433;")
        layout.addWidget(title)

        form = QtWidgets.QFormLayout()
        form.setSpacing(10)

        # Source (file or folder)
        src_row = QtWidgets.QHBoxLayout()
        self.source_edit = QtWidgets.QLineEdit()
        self.source_edit.setPlaceholderText("FASTA file, folder of allele FASTAs, .zip, or a chewBBACA schema")
        browse_file = QtWidgets.QPushButton("File…")
        browse_file.setStyleSheet(_SECONDARY_BTN)
        browse_file.clicked.connect(self._browse_file)
        browse_dir = QtWidgets.QPushButton("Folder…")
        browse_dir.setStyleSheet(_SECONDARY_BTN)
        browse_dir.clicked.connect(self._browse_dir)
        src_row.addWidget(self.source_edit, 1)
        src_row.addWidget(browse_file)
        src_row.addWidget(browse_dir)
        form.addRow("Source", src_row)

        # Species
        self.species_edit = QtWidgets.QLineEdit()
        self.species_edit.setPlaceholderText("e.g. Acinetobacter baumannii")
        self.species_edit.textChanged.connect(self._sync_target_hint)
        form.addRow("Species", self.species_edit)

        # MLST aliases
        self.mlst_edit = QtWidgets.QLineEdit()
        self.mlst_edit.setPlaceholderText("comma-separated, e.g. abaumannii, abaumannii_2 (optional)")
        form.addRow("MLST scheme(s)", self.mlst_edit)

        # Single-FASTA interpretation
        self.mode_combo = QtWidgets.QComboBox()
        self.mode_combo.addItems([
            "Auto-detect",
            "Set of loci (split by header)",
            "Alleles of a single locus",
        ])
        form.addRow("If a single FASTA", self.mode_combo)

        layout.addLayout(form)

        self.target_hint = QtWidgets.QLabel("")
        self.target_hint.setWordWrap(True)
        self.target_hint.setStyleSheet(_HINT_STYLE)
        layout.addWidget(self.target_hint)

        hint = QtWidgets.QLabel(
            "A folder may hold per-locus FASTAs (a Ridom combined all-loci file is ignored). "
            "The schema is written under DB/cgmlst and auto-selected when a run's MLST matches."
        )
        hint.setWordWrap(True)
        hint.setStyleSheet(_HINT_STYLE)
        layout.addWidget(hint)

        layout.addStretch(1)
        btns = QtWidgets.QHBoxLayout()
        btns.addStretch(1)
        cancel = QtWidgets.QPushButton("Cancel")
        cancel.setStyleSheet(_SECONDARY_BTN)
        cancel.clicked.connect(self.reject)
        btns.addWidget(cancel)
        self.import_btn = QtWidgets.QPushButton("Import")
        self.import_btn.setStyleSheet(_PRIMARY_BTN)
        self.import_btn.clicked.connect(self._do_import)
        btns.addWidget(self.import_btn)
        layout.addLayout(btns)

    # ---- helpers ----
    def _browse_file(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Select schema source file", "",
            "Schema sources (*.fasta *.fa *.fna *.ffn *.faa *.zip);;All files (*)",
        )
        if path:
            self.source_edit.setText(path)

    def _browse_dir(self) -> None:
        path = QtWidgets.QFileDialog.getExistingDirectory(self, "Select schema source folder")
        if path:
            self.source_edit.setText(path)

    def _cgmlst_dir(self) -> Optional[Path]:
        from pipeline.ToolEnv import configured_database_dirs

        raw = configured_database_dirs().get("cgmlst_dir") or ""
        return Path(raw) if raw else None

    def _sync_target_hint(self) -> None:
        from pipeline.ChewbbacaSchemaBuild import slugify_species

        species = self.species_edit.text().strip()
        cg = self._cgmlst_dir()
        if not species or cg is None:
            self.target_hint.setText("")
            return
        schema_id = slugify_species(species)
        self.target_hint.setText(f"Will be installed as: {cg / schema_id}")

    # ---- import ----
    def _resolve_mode(self, source: Path) -> str:
        from pipeline.ChewbbacaSchemaBuild import is_combined_fasta

        if not (source.is_file() and source.suffix.lower() in _FASTA_EXTS):
            return "loci_set"  # ignored for folders/zips/prebuilt
        choice = self.mode_combo.currentIndex()
        if choice == 1:
            return "loci_set"
        if choice == 2:
            return "single_locus"
        # auto
        try:
            return "loci_set" if is_combined_fasta(source) else "single_locus"
        except Exception:
            return "loci_set"

    def _do_import(self) -> None:
        from pipeline.ChewbbacaSchemaBuild import slugify_species

        source_text = self.source_edit.text().strip()
        species = self.species_edit.text().strip()
        if not source_text:
            QtWidgets.QMessageBox.warning(self, "Missing source", "Choose a schema source file or folder.")
            return
        source = Path(source_text).expanduser()
        if not source.exists():
            QtWidgets.QMessageBox.warning(self, "Missing source", f"Source does not exist:\n{source}")
            return
        if not species:
            QtWidgets.QMessageBox.warning(self, "Missing species", "Enter a species name.")
            return
        cg = self._cgmlst_dir()
        if cg is None:
            QtWidgets.QMessageBox.critical(
                self, "No cgMLST directory",
                "databases.cgmlst_dir is not configured. Run first-time setup first.",
            )
            return

        schema_id = slugify_species(species)
        target_dir = cg / schema_id
        overwrite = False
        if target_dir.exists():
            resp = QtWidgets.QMessageBox.question(
                self, "Replace schema?",
                f"A schema '{schema_id}' already exists.\nReplace it?",
                QtWidgets.QMessageBox.StandardButton.Yes | QtWidgets.QMessageBox.StandardButton.No,
            )
            if resp != QtWidgets.QMessageBox.StandardButton.Yes:
                return
            overwrite = True

        mode = self._resolve_mode(source)
        mlst = [a.strip() for a in self.mlst_edit.text().replace(";", ",").split(",") if a.strip()]

        # Run the conversion on a worker thread with an indeterminate progress dialog.
        worker = _ImportWorker(source, target_dir, mode, overwrite)
        prog = QtWidgets.QProgressDialog("Preparing…", "Cancel", 0, 0, self)
        prog.setWindowTitle("Importing schema")
        prog.setWindowModality(QtCore.Qt.WindowModality.WindowModal)
        prog.setMinimumDuration(0)
        prog.setAutoClose(False)
        prog.setAutoReset(False)

        state: dict = {}
        loop = QtCore.QEventLoop()
        worker.progress.connect(prog.setLabelText)
        worker.done.connect(lambda res: (state.update(result=res), loop.quit()))
        worker.failed.connect(lambda msg: (state.update(error=msg), loop.quit()))
        prog.canceled.connect(worker.request_cancel)

        self.import_btn.setEnabled(False)
        worker.start()
        prog.show()
        loop.exec()  # returns when the worker emits done/failed
        prog.close()
        worker.wait()
        self.import_btn.setEnabled(True)

        if "error" in state:
            if state["error"] == "cancelled":
                return
            QtWidgets.QMessageBox.critical(self, "Import failed", state["error"])
            return
        result = state.get("result")
        if result is None:
            return

        # Register in the catalog (+ config available_schemas).
        try:
            from pipeline.ChewbbacaSchemaBuild import register_schema
            from pipeline.ToolEnv import configured_chewbbaca_settings

            catalog_path = configured_chewbbaca_settings().get("catalog_path") or None
            register_schema(
                cg, schema_id, [species], mlst,
                catalog_path=Path(catalog_path) if catalog_path else None,
            )
        except ValueError as exc:
            QtWidgets.QMessageBox.warning(
                self, "Schema imported, not auto-selectable",
                f"The schema was built ({result.written_loci} loci) but registration reported:\n{exc}\n\n"
                "It is available for manual selection.",
            )
            self.imported_schema_id = schema_id
            self.accept()
            return
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.warning(
                self, "Schema imported, registration issue",
                f"The schema was built but could not be registered:\n{exc}",
            )

        self.imported_schema_id = schema_id
        warn = f"\n\nNote: {result.warning}" if result.warning else ""
        miss = f"\n{result.missing_loci} expected loci were missing." if result.missing_loci else ""
        QtWidgets.QMessageBox.information(
            self, "Schema imported",
            f"Installed '{schema_id}' with {result.written_loci} loci "
            f"(source: {result.source_format}).{miss}{warn}",
        )
        self.accept()
