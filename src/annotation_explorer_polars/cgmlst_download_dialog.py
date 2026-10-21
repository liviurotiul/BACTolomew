"""Dialog to download cgMLST schemes from cgMLST.org into DB/cgmlst.

Lists the cgMLST.org schemes (with their download size when the server gives one),
shows the cgMLST.org terms once before the first download and records their
acceptance, then runs pipeline.CgmlstInstall.install on a background thread:
each scheme is downloaded, adapted with chewBBACA PrepExternalSchema, checked with
a test allele call and registered for the MLST-based auto-selection.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

from PyQt6 import QtCore, QtWidgets

from .schema_import_dialog import _hint, run_install_job, show_install_summary


def _size_text(n: Optional[int]) -> str:
    if not n:
        return ""
    value = float(n)
    for unit in ("B", "KB", "MB", "GB"):
        if value < 1024 or unit == "GB":
            return f"{value:.0f} {unit}" if unit in ("B", "KB") else f"{value:.1f} {unit}"
        value /= 1024
    return ""


def confirm_cgmlst_terms(parent: Optional[QtWidgets.QWidget], dest: Path, *, before_run: bool = False) -> bool:
    """The cgMLST.org terms, shown once: True when they are accepted (now or before).

    ``before_run``: asked before a run that may install its species' scheme; "Not now"
    starts the run anyway (it stops at the chewBBACA step if it needs a download).
    """
    from pipeline.CgmlstInstall import TERMS_TEXT, accept_terms, default_dest, terms_accepted

    if terms_accepted(dest) or (before_run and terms_accepted(default_dest())):
        return True
    box = QtWidgets.QMessageBox(parent)
    box.setIcon(QtWidgets.QMessageBox.Icon.Information)
    box.setWindowTitle("cgMLST.org terms of use")
    box.setText(TERMS_TEXT)
    if before_run:
        box.setInformativeText(
            "This run installs the cgMLST scheme of its species from cgMLST.org if it is not "
            "installed yet. Without the terms accepted, such a run stops at the chewBBACA step."
        )
    accept = box.addButton("Accept", QtWidgets.QMessageBox.ButtonRole.AcceptRole)
    box.addButton("Not now" if before_run else "Cancel", QtWidgets.QMessageBox.ButtonRole.RejectRole)
    box.exec()
    if box.clickedButton() is not accept:
        return False
    try:
        accept_terms(dest)
    except OSError as exc:
        QtWidgets.QMessageBox.critical(parent, "cgMLST.org terms of use", f"Could not record the acceptance:\n{exc}")
        return False
    return True


# List workers still fetching when their dialog closed (kept alive until they finish).
_ORPHANED: set = set()


class _ListWorker(QtCore.QThread):
    """Fetches the cgMLST.org index and the download sizes."""

    done = QtCore.pyqtSignal(object, object)  # (index, rows)
    failed = QtCore.pyqtSignal(str)

    def __init__(self, dest: Path) -> None:
        super().__init__()
        self._dest = dest

    def run(self) -> None:  # pragma: no cover - exercised via the GUI
        try:
            from pipeline.CgmlstInstall import available, fetch_org_index

            index = fetch_org_index()
            rows = available(self._dest, index=index, sizes=True)
            self.done.emit(index, rows)
        except Exception as exc:  # noqa: BLE001 - surfaced in the dialog
            self.failed.emit(str(exc))


class CgmlstDownloadDialog(QtWidgets.QDialog):
    def __init__(self, dest: Path, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Download from cgMLST.org")
        self.setObjectName("CgmlstDownloadDialog")
        self.resize(620, 560)
        self._dest = Path(dest)
        self._index: list = []
        self.installed_any = False
        self._build_ui()
        self._worker = _ListWorker(self._dest)
        self._worker.done.connect(self._fill)
        self._worker.failed.connect(self._fetch_failed)
        self._worker.start()

    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(22, 20, 22, 18)
        layout.setSpacing(10)
        layout.setSizeConstraint(QtWidgets.QLayout.SizeConstraint.SetMinimumSize)

        title = QtWidgets.QLabel("Download from cgMLST.org")
        title.setProperty("role", "windowTitle")
        layout.addWidget(title)
        layout.addSpacing(4)

        caption = QtWidgets.QLabel("Species")
        caption.setProperty("role", "settingsCaption")
        layout.addWidget(caption)
        self.scheme_list = QtWidgets.QListWidget()
        self.scheme_list.setMinimumHeight(260)
        self.scheme_list.itemChanged.connect(self._sync_button)
        layout.addWidget(self.scheme_list, stretch=1)
        self.status_hint = _hint("Loading the cgMLST.org list…")
        layout.addWidget(self.status_hint)

        layout.addWidget(_hint(
            "Each scheme is downloaded, adapted with chewBBACA PrepExternalSchema, checked with a "
            f"test allele call and added to the cgMLST schemas in {self._dest}, where a run "
            "auto-selects it from its MLST result. Installed schemes are not downloaded again. "
            "cgMLST.org schemes are for non-commercial use only."
        ))

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
        self.install_btn = QtWidgets.QPushButton("Download and install")
        self.install_btn.setProperty("role", "primary")
        self.install_btn.setEnabled(False)
        self.install_btn.clicked.connect(self._install)
        btns.addWidget(self.install_btn)
        layout.addWidget(footer)

    # ---- list ----
    def _fill(self, index: list, rows: list) -> None:
        self._index = list(index)
        self.scheme_list.blockSignals(True)
        self.scheme_list.clear()
        for row in rows:
            if not row.slug:
                continue  # only in the catalog: not downloadable
            parts = [row.title]
            size = _size_text(row.size_bytes)
            if size:
                parts.append(size)
            if row.installed:
                parts.append("installed")
            elif not row.bacterial:
                parts.append("not bacterial, not supported by chewBBACA")
            item = QtWidgets.QListWidgetItem("  ·  ".join(parts))
            item.setData(QtCore.Qt.ItemDataRole.UserRole, row.slug)
            if row.installed or not row.bacterial:
                item.setFlags(QtCore.Qt.ItemFlag.NoItemFlags)
            else:
                item.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
                item.setCheckState(QtCore.Qt.CheckState.Unchecked)
            self.scheme_list.addItem(item)
        self.scheme_list.blockSignals(False)
        self.status_hint.setText(
            f"{sum(1 for r in rows if r.slug)} schemes on cgMLST.org, "
            f"{sum(1 for r in rows if r.slug and r.installed)} installed."
        )

    def _fetch_failed(self, error: str) -> None:
        self.status_hint.setText(f"cgMLST.org could not be reached: {error}")

    def _checked_slugs(self) -> list[str]:
        out = []
        for i in range(self.scheme_list.count()):
            item = self.scheme_list.item(i)
            if item.flags() & QtCore.Qt.ItemFlag.ItemIsUserCheckable and \
                    item.checkState() == QtCore.Qt.CheckState.Checked:
                out.append(str(item.data(QtCore.Qt.ItemDataRole.UserRole)))
        return out

    def _sync_button(self, *_args) -> None:
        self.install_btn.setEnabled(bool(self._checked_slugs()))

    # ---- install ----
    def _install(self) -> None:
        from pipeline.CgmlstInstall import install

        slugs = self._checked_slugs()
        if not slugs or not confirm_cgmlst_terms(self, self._dest):
            return
        index = list(self._index)
        self.install_btn.setEnabled(False)
        results, error = run_install_job(
            self, "Installing cgMLST schemas",
            lambda **cb: install(slugs, self._dest, index=index, accept=True, **cb),
            len(slugs),
        )
        self.install_btn.setEnabled(True)
        if error or results is None:
            QtWidgets.QMessageBox.critical(self, "Download failed", error or "no result")
            return
        show_install_summary(self, "Download from cgMLST.org", results)
        self.installed_any = any(r.status == "installed" for r in results)
        if self.installed_any:
            self.accept()

    def done(self, result: int) -> None:  # noqa: D401 - QDialog override
        worker = self._worker
        if worker.isRunning():
            worker.done.disconnect()
            worker.failed.disconnect()
            _ORPHANED.add(worker)
            worker.finished.connect(lambda: _ORPHANED.discard(worker))
        super().done(result)
