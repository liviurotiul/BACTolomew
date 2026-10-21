"""Dialog to import cgMLST schemas into DB/cgmlst and register them.

Accepts a folder of per-locus allele FASTAs, a .zip of them, a single combined
FASTA, a single-locus FASTA, or an already-built chewBBACA schema; several files,
or a folder holding many of these (every zip downloaded from cgMLST.org), are
imported one by one with the species taken from each name. The work is
pipeline.CgmlstInstall.import_paths on a background thread; each schema is
registered in scheme_species_catalog.json (species + MLST aliases) so it
participates in the existing MLST-based auto-selection.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, Optional, Sequence

from PyQt6 import QtCore, QtWidgets



def _hint(text: str = "") -> QtWidgets.QLabel:
    label = QtWidgets.QLabel(text)
    label.setWordWrap(True)
    label.setProperty("role", "windowHint")
    return label


def _field(caption: str, control) -> QtWidgets.QVBoxLayout:
    """A caption above its control (a widget or a layout)."""
    box = QtWidgets.QVBoxLayout()
    box.setSpacing(4)
    label = QtWidgets.QLabel(caption)
    label.setProperty("role", "settingsCaption")
    box.addWidget(label)
    if isinstance(control, QtWidgets.QLayout):
        box.addLayout(control)
    else:
        box.addWidget(control)
    return box


_SOURCE_SEPARATOR = "; "


class _InstallWorker(QtCore.QThread):
    """Runs one pipeline.CgmlstInstall job (install or import_paths) off the GUI thread."""

    progress = QtCore.pyqtSignal(str)
    fraction = QtCore.pyqtSignal(str, float)  # (schema id, share of that item done)
    done = QtCore.pyqtSignal(object)
    failed = QtCore.pyqtSignal(str)

    def __init__(self, job: Callable[..., list]) -> None:
        super().__init__()
        self._job = job
        self._cancel = False

    def request_cancel(self) -> None:
        self._cancel = True

    def run(self) -> None:  # pragma: no cover - exercised via the GUI
        try:
            results = self._job(
                progress=self.progress.emit,
                on_result=lambda res, _done, _total: self.fraction.emit(res.schema_id or res.source, 1.0),
                should_cancel=lambda: self._cancel,
                on_fraction=self.fraction.emit,
            )
            self.done.emit(results)
        except Exception as exc:  # noqa: BLE001 - surfaced to the user
            self.failed.emit(str(exc))


def run_install_job(
    parent: QtWidgets.QWidget, title: str, job: Callable[..., list], total: int
) -> tuple[Optional[list], str]:
    """Run ``job`` on a worker thread behind a progress window; the GUI stays responsive.

    The bar follows the items' downloads and build phases (``on_fraction``); with
    ``total`` 0 (not known) it only shows activity. Returns (results, "") or
    (None, error); a cancelled job still returns its results (the unfinished items say
    "cancelled").
    """
    worker = _InstallWorker(job)
    scale = 1000
    prog = QtWidgets.QProgressDialog("Preparing…", "Cancel", 0, scale if total > 0 else 0, parent)
    prog.setWindowTitle(title)
    prog.setWindowModality(QtCore.Qt.WindowModality.WindowModal)
    prog.setMinimumDuration(0)
    prog.setAutoClose(False)
    prog.setAutoReset(False)
    prog.setMinimumWidth(460)

    state: dict = {}
    loop = QtCore.QEventLoop()
    worker.progress.connect(prog.setLabelText)
    shares: dict[str, float] = {}

    def on_fraction(item: str, share: float) -> None:
        shares[item] = max(share, shares.get(item, 0.0))
        if total > 0:
            prog.setValue(int(scale * min(1.0, sum(shares.values()) / total)))

    worker.fraction.connect(on_fraction)
    worker.done.connect(lambda res: (state.update(result=res), loop.quit()))
    worker.failed.connect(lambda msg: (state.update(error=msg), loop.quit()))
    prog.canceled.connect(worker.request_cancel)
    prog.canceled.connect(lambda: prog.setLabelText("Cancelling after the current step…"))

    worker.start()
    prog.show()
    loop.exec()  # returns when the worker emits done/failed
    prog.close()
    worker.wait()
    return state.get("result"), state.get("error", "")


def show_install_summary(parent: QtWidgets.QWidget, title: str, results: Sequence) -> None:
    """One message listing what was installed, skipped and failed (several sources)."""
    lines: list[str] = []
    done = [r for r in results if r.status in ("installed", "imported")]
    skipped = [r for r in results if r.status == "skipped"]
    failed = [r for r in results if r.status == "failed"]
    if done:
        lines.append("Installed: " + ", ".join(r.schema_id for r in done))
        for r in done:
            notes = [n for n in (r.warning, r.registration_error) if n]
            if r.backup:
                notes.append(f"previous version kept in {r.backup}")
            if notes:
                lines.append(f"  {r.schema_id}: " + "; ".join(notes))
    if skipped:
        lines.append("Skipped (already installed): " + ", ".join(r.schema_id or Path(r.source).name for r in skipped))
    if failed:
        lines.append("Failed:")
        lines += [f"  {r.schema_id or Path(r.source).name}: {r.reason}" for r in failed]
    box = QtWidgets.QMessageBox.warning if failed else QtWidgets.QMessageBox.information
    box(parent, title, "\n".join(lines) or "Nothing to do.")


class AddSchemaDialog(QtWidgets.QDialog):
    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Add cgMLST schema")
        self.setObjectName("AddSchemaDialog")
        self.resize(580, 500)
        self.imported_schema_id: Optional[str] = None
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(22, 20, 22, 18)
        layout.setSpacing(10)
        # The window never gets smaller than its content (no overlapping controls).
        layout.setSizeConstraint(QtWidgets.QLayout.SizeConstraint.SetMinimumSize)

        title = QtWidgets.QLabel("Add cgMLST schema")
        title.setProperty("role", "windowTitle")
        layout.addWidget(title)
        layout.addSpacing(4)

        # Source (file or folder)
        src_row = QtWidgets.QHBoxLayout()
        self.source_edit = QtWidgets.QLineEdit()
        self.source_edit.setPlaceholderText("FASTA file, folder of allele FASTAs, .zip, or a chewBBACA schema")
        browse_file = QtWidgets.QPushButton("File…")
        browse_file.clicked.connect(self._browse_file)
        browse_dir = QtWidgets.QPushButton("Folder…")
        browse_dir.clicked.connect(self._browse_dir)
        src_row.addWidget(self.source_edit, 1)
        src_row.addWidget(browse_file)
        src_row.addWidget(browse_dir)
        layout.addLayout(_field("Source", src_row))

        # Species
        self.species_edit = QtWidgets.QLineEdit()
        self.species_edit.setPlaceholderText("e.g. Acinetobacter baumannii")
        self.species_edit.textChanged.connect(self._sync_target_hint)
        layout.addLayout(_field("Species", self.species_edit))

        # MLST aliases
        self.mlst_edit = QtWidgets.QLineEdit()
        self.mlst_edit.setPlaceholderText("comma-separated, e.g. abaumannii, abaumannii_2 (optional)")
        layout.addLayout(_field("MLST scheme(s)", self.mlst_edit))

        # Single-FASTA interpretation
        self.mode_combo = QtWidgets.QComboBox()
        self.mode_combo.addItems([
            "Auto-detect",
            "Set of loci (split by header)",
            "Alleles of a single locus",
        ])
        layout.addLayout(_field("If a single FASTA", self.mode_combo))

        self.target_hint = _hint()
        layout.addWidget(self.target_hint)

        layout.addWidget(_hint(
            "Several files, or a folder holding many schemas (for example every zip downloaded "
            "from cgMLST.org), are imported one by one; each species is then taken from its file "
            "or folder name. An installed schema is replaced only when you confirm it, and the "
            "previous version is kept in the .replaced folder."
        ))

        hint = _hint(
            "A folder may hold per-locus FASTAs (a Ridom combined all-loci file is ignored). "
            "Loci are adapted with chewBBACA PrepExternalSchema: alleles are put in coding "
            "orientation and alleles that are not a complete CDS are dropped (listed in the "
            "schema's schema_build_report folder). "
            "The schema is written under DB/cgmlst and auto-selected when a run's MLST matches."
        )
        layout.addWidget(hint)

        layout.addStretch(1)
        footer = QtWidgets.QFrame()
        footer.setObjectName("StartFooter")
        btns = QtWidgets.QHBoxLayout(footer)
        btns.setContentsMargins(0, 12, 0, 0)
        btns.setSpacing(8)
        btns.addStretch(1)
        cancel = QtWidgets.QPushButton("Cancel")
        cancel.setProperty("role", "flat")
        cancel.clicked.connect(self.reject)
        btns.addWidget(cancel)
        self.import_btn = QtWidgets.QPushButton("Import")
        self.import_btn.setProperty("role", "primary")
        self.import_btn.clicked.connect(self._do_import)
        btns.addWidget(self.import_btn)
        layout.addWidget(footer)

    # ---- helpers ----
    def showEvent(self, event) -> None:  # noqa: N802 - Qt override
        super().showEvent(event)
        QtCore.QTimer.singleShot(0, self._fit_height)

    def _fit_height(self) -> None:
        """Grow the window to the height its wrapped hints need at its width (the minimum
        size constraint does not follow a label's height for width)."""
        layout = self.layout()
        if layout is None:
            return
        need = (layout.totalHeightForWidth(self.width()) if layout.hasHeightForWidth()
                else layout.totalMinimumSize().height())
        if need > self.height():
            self.resize(self.width(), need)

    def _browse_file(self) -> None:
        paths, _ = QtWidgets.QFileDialog.getOpenFileNames(
            self, "Select schema source file", "",
            "Schema sources (*.fasta *.fa *.fna *.ffn *.faa *.zip);;All files (*)",
        )
        if paths:
            self.source_edit.setText(_SOURCE_SEPARATOR.join(paths))

    def _browse_dir(self) -> None:
        path = QtWidgets.QFileDialog.getExistingDirectory(self, "Select schema source folder")
        if path:
            self.source_edit.setText(path)

    def _cgmlst_dir(self) -> Optional[Path]:
        from pipeline.ToolEnv import configured_database_dirs

        raw = configured_database_dirs().get("cgmlst_dir") or ""
        return Path(raw) if raw else None

    def _sync_target_hint(self) -> None:
        from pipeline.CgmlstInstall import schema_id_for_species

        species = self.species_edit.text().strip()
        cg = self._cgmlst_dir()
        if not species or cg is None:
            self.target_hint.setText("")
            return
        try:  # the same resolver as the import (a catalog schema listing the species)
            schema_id = schema_id_for_species(species, cg)
        except Exception:  # noqa: BLE001 - an unreadable catalog fails the import itself
            self.target_hint.setText("")
            return
        self.target_hint.setText(f"Will be installed as: {cg / schema_id}")
        self._fit_height()

    # ---- import ----
    def _sources(self) -> list[Path]:
        return [Path(t.strip()).expanduser() for t in self.source_edit.text().split(";") if t.strip()]

    def _single_fasta_mode(self) -> str:
        return ("auto", "loci_set", "single_locus")[max(0, self.mode_combo.currentIndex())]

    def _do_import(self) -> None:
        from pipeline.CgmlstInstall import import_paths

        sources = self._sources()
        species = self.species_edit.text().strip()
        if not sources:
            QtWidgets.QMessageBox.warning(self, "Missing source", "Choose a schema source file or folder.")
            return
        for source in sources:
            if not source.exists():
                QtWidgets.QMessageBox.warning(self, "Missing source", f"Source does not exist:\n{source}")
                return
        cg = self._cgmlst_dir()
        if cg is None:
            QtWidgets.QMessageBox.critical(
                self, "No cgMLST directory",
                "databases.cgmlst_dir is not configured. Run first-time setup first.",
            )
            return
        mlst = [a.strip() for a in self.mlst_edit.text().replace(";", ",").split(",") if a.strip()]
        common = dict(species=species or None, dest=cg, mlst_schemes=mlst,
                      single_fasta_mode=self._single_fasta_mode())

        # Plan first (no build): what each source is detected as, and what is installed.
        # On the worker too: detection may read every source and the cgMLST.org lists.
        self.import_btn.setEnabled(False)
        plan, error = run_install_job(
            self, "Importing schema", lambda **cb: import_paths(sources, dry_run=True, **common, **cb), 0,
        )
        self.import_btn.setEnabled(True)
        if error == "cancelled":
            return
        if error or plan is None:
            QtWidgets.QMessageBox.critical(self, "Import failed", error or "no result")
            return
        single = len(plan) == 1
        if single and not species and plan[0].status == "failed" and "species" in plan[0].reason:
            QtWidgets.QMessageBox.warning(self, "Missing species", "Enter a species name.")
            return
        installed = [r for r in plan if r.status == "skipped" and r.reason.startswith("already installed")]
        replace = False
        if installed:
            if single:
                resp = QtWidgets.QMessageBox.question(
                    self, "Replace schema?",
                    f"A schema '{installed[0].schema_id}' already exists.\nReplace it?",
                    QtWidgets.QMessageBox.StandardButton.Yes | QtWidgets.QMessageBox.StandardButton.No,
                )
                if resp != QtWidgets.QMessageBox.StandardButton.Yes:
                    return
                replace = True
            else:
                resp = QtWidgets.QMessageBox.question(
                    self, "Replace schemas?",
                    f"{len(installed)} of these schemas are already installed: "
                    f"{', '.join(r.schema_id for r in installed)}.\n"
                    "Replace them? The installed versions are kept in the .replaced folder; "
                    "with No they are skipped.",
                    QtWidgets.QMessageBox.StandardButton.Yes | QtWidgets.QMessageBox.StandardButton.No,
                )
                replace = resp == QtWidgets.QMessageBox.StandardButton.Yes
        if not any(r.status == "planned" for r in plan) and not replace:
            show_install_summary(self, "Add cgMLST schema", plan)
            return

        self.import_btn.setEnabled(False)
        results, error = run_install_job(
            self, "Importing schema",
            lambda **cb: import_paths(sources, replace=replace, **common, **cb),
            sum(r.status == "planned" for r in plan) + (len(installed) if replace else 0),
        )
        self.import_btn.setEnabled(True)
        if error or results is None:
            QtWidgets.QMessageBox.critical(self, "Import failed", error or "no result")
            return
        done = [r for r in results if r.status == "imported"]
        if done:
            self.imported_schema_id = done[0].schema_id
        if not single:
            show_install_summary(self, "Add cgMLST schema", results)
            if done:
                self.accept()
            return

        result = results[0]
        if result.status != "imported":
            if result.reason != "cancelled":
                QtWidgets.QMessageBox.critical(self, "Import failed", result.reason)
            return
        if result.alias_conflict:
            QtWidgets.QMessageBox.warning(
                self, "Schema imported, not auto-selectable",
                f"The schema was built ({result.loci} loci) but registration reported:\n"
                f"{result.registration_error}\n\nIt is available for manual selection.",
            )
            self.accept()
            return
        if result.registration_error:
            QtWidgets.QMessageBox.warning(
                self, "Schema imported, registration issue",
                f"The schema was built but could not be registered:\n{result.registration_error}",
            )
        warn = f"\n\nNote: {result.warning}" if result.warning else ""
        if result.backup:
            warn += f"\n\nThe previous version is kept in {result.backup}."
        miss = (
            f"\n{result.missing_loci} expected loci were missing or had no valid CDS allele."
            if result.missing_loci else ""
        )
        QtWidgets.QMessageBox.information(
            self, "Schema imported",
            f"Installed '{result.schema_id}' with {result.loci} loci "
            f"(source: {result.source_format}).{miss}{warn}",
        )
        self.accept()
