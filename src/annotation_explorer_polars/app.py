from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Optional
from urllib import request

from PyQt6 import QtCore, QtWidgets

from .first_time_setup import FirstTimeSetupDialog, apply_first_time_setup, collect_setup_missing_paths
from .loading_splash import LoadingSplash
from .progress import ProgressLogger
from .startup_dialog import StartupDialog


WORLD_GEOJSON_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/"
    "ne_110m_admin_0_countries.geojson"
)
WORLD_GEOJSON_50M_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/"
    "ne_50m_admin_0_countries.geojson"
)
WORLD_GEOJSON_10M_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/"
    "ne_10m_admin_0_countries.geojson"
)
WORLD_GEOJSON_MAX_SIZE = 5 * 1024 * 1024


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_app_config() -> dict:
    path = project_root() / "bactolomew_config.json"
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


def _missing_startup_paths(*, dev_mode: bool = False) -> list[Path]:
    root = project_root()
    config_path = root / "bactolomew_config.json"
    config = _load_app_config()
    missing: list[Path] = [config_path] if not config_path.exists() else []
    missing.extend(collect_setup_missing_paths(root, config, dev_mode=dev_mode))
    return missing


def _run_first_time_setup(missing: list[Path], *, dev_mode: bool = False) -> bool:
    dlg = FirstTimeSetupDialog(missing=missing, dev_mode=dev_mode)
    if dlg.exec() != QtWidgets.QDialog.DialogCode.Accepted:
        return False
    result = dlg.result()
    if result is None:
        return False
    try:
        apply_first_time_setup(result)
    except Exception as exc:
        QtWidgets.QMessageBox.critical(
            None,
            "Setup failed",
            str(exc),
        )
        return False
    return True


class _ToolsImageWorker(QtCore.QThread):
    line = QtCore.pyqtSignal(str)
    done = QtCore.pyqtSignal(bool, str)

    def __init__(self, runtime, image, registry, sif, repo_root, parent=None) -> None:
        super().__init__(parent)
        self._args = (runtime, image, registry, sif, repo_root)

    def run(self) -> None:  # worker thread
        from pipeline import ContainerRuntime as CR

        runtime, image, registry, sif, repo_root = self._args
        try:
            CR.ensure_tools_image(
                runtime, image=image, registry=registry, sif_path=sif,
                repo_root=Path(repo_root), log=lambda m: self.line.emit(m),
            )
            self.done.emit(True, "")
        except Exception as exc:  # surfaced in the dialog
            self.done.emit(False, str(exc))


def _ensure_tools_image_ready() -> bool:
    """In container mode, make sure the tools image is present (build/pull with progress)."""
    from pipeline.ToolEnv import runtime_is_container

    config = _load_app_config()
    if not runtime_is_container(config):
        return True

    from pipeline import ContainerRuntime as CR

    container = config.get("container") if isinstance(config.get("container"), dict) else {}
    image = str(container.get("image") or "").strip() or "bactolomew-tools:0.1.0"
    runtime = str(config.get("runtime") or "docker").strip().lower()
    registry = str(container.get("registry") or "").strip()
    sif = str(container.get("sif_path") or "").strip()

    if runtime == "singularity":
        if CR.sif_present(sif):
            return True
    else:
        if not CR.docker_available():
            QtWidgets.QMessageBox.critical(
                None, "Docker not available",
                "Docker is the selected tool runtime but the Docker daemon is not reachable.\n"
                'Start Docker (or edit "runtime" in bactolomew_config.json) and relaunch.',
            )
            return False
        if CR.docker_image_present(image):
            return True

    # Acquire the image with a streaming progress dialog.
    dlg = QtWidgets.QDialog()
    dlg.setWindowTitle("Preparing BACTolomew tools")
    dlg.setModal(True)
    dlg.setWindowFlags(
        QtCore.Qt.WindowType.Dialog
        | QtCore.Qt.WindowType.CustomizeWindowHint
        | QtCore.Qt.WindowType.WindowTitleHint
    )
    dlg.resize(720, 440)
    lay = QtWidgets.QVBoxLayout(dlg)
    hint = "Pulling from registry" if registry else "Building the tool image locally (first time can take many minutes)"
    status = QtWidgets.QLabel(f"Setting up the tool image: {image}\n{hint}…")
    status.setWordWrap(True)
    lay.addWidget(status)
    log = QtWidgets.QPlainTextEdit()
    log.setReadOnly(True)
    log.setStyleSheet("font-family: monospace; font-size: 11px;")
    lay.addWidget(log, 1)

    outcome = {"ok": False}
    worker = _ToolsImageWorker(runtime, image, registry, sif, str(project_root()))

    def _append(text: str) -> None:
        log.appendPlainText(text)
        bar = log.verticalScrollBar()
        bar.setValue(bar.maximum())

    def _finish(ok: bool, err: str) -> None:
        outcome["ok"] = ok
        if ok:
            dlg.accept()
        else:
            status.setText(f"Setup failed: {err}")
            btns = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Close)
            btns.rejected.connect(dlg.reject)
            btns.accepted.connect(dlg.reject)
            lay.addWidget(btns)

    worker.line.connect(_append)
    worker.done.connect(_finish)
    worker.start()
    dlg.exec()
    worker.wait()
    if not outcome["ok"]:
        QtWidgets.QMessageBox.critical(
            None, "Tool image setup failed",
            "Could not build or pull the BACTolomew tools image. See the log for details.",
        )
    return outcome["ok"]


def _ensure_startup_ready(*, dev_mode: bool = False) -> bool:
    if not _ensure_tools_image_ready():
        return False
    missing = _missing_startup_paths(dev_mode=dev_mode)
    if missing:
        return _run_first_time_setup(missing, dev_mode=dev_mode)
    return True


def _download_world_geojson(target: Path, url: str) -> Optional[str]:
    try:
        with request.urlopen(url, timeout=30) as response:
            if response.status != 200:
                return f"HTTP {response.status}"
            data = response.read()
        tmp_path = target.with_suffix(".tmp")
        tmp_path.write_bytes(data)
        tmp_path.replace(target)
    except Exception as exc:
        return str(exc)
    return None


def _ensure_world_geojson() -> None:
    assets_dir = Path(__file__).resolve().parent / "assets"
    target = assets_dir / "world_countries.geojson"
    if target.exists():
        try:
            if target.stat().st_size <= WORLD_GEOJSON_MAX_SIZE:
                pass
            else:
                err = _download_world_geojson(target, WORLD_GEOJSON_URL)
                if err:
                    print(f"Warning: Failed to download world_countries.geojson: {err}")
        except OSError:
            pass
    else:
        assets_dir.mkdir(parents=True, exist_ok=True)
        err = _download_world_geojson(target, WORLD_GEOJSON_URL)
        if err:
            print(f"Warning: Failed to download world_countries.geojson: {err}")

    high_target = assets_dir / "world_countries_50m.geojson"
    if not high_target.exists():
        err = _download_world_geojson(high_target, WORLD_GEOJSON_50M_URL)
        if err:
            print(f"Warning: Failed to download world_countries_50m.geojson: {err}")

    ultra_target = assets_dir / "world_countries_10m.geojson"
    if not ultra_target.exists():
        err = _download_world_geojson(ultra_target, WORLD_GEOJSON_10M_URL)
        if err:
            print(f"Warning: Failed to download world_countries_10m.geojson: {err}")


class _WorldGeojsonWorkerSignals(QtCore.QObject):
    finished = QtCore.pyqtSignal(list)


class _WorldGeojsonWorker(QtCore.QRunnable):
    def __init__(self) -> None:
        super().__init__()
        self.signals = _WorldGeojsonWorkerSignals()
        self.setAutoDelete(True)

    def run(self) -> None:
        warnings: list[str] = []
        try:
            _ensure_world_geojson()
        except Exception as exc:
            warnings.append(str(exc))
        self.signals.finished.emit(warnings)


def _start_world_geojson_worker() -> None:
    worker = _WorldGeojsonWorker()

    def _log_warnings(warnings: list) -> None:
        for warning in warnings:
            print(f"Warning: Failed to prepare world GeoJSON assets: {warning}")

    worker.signals.finished.connect(_log_warnings)
    QtCore.QThreadPool.globalInstance().start(worker)


def _validate_metadata_csv(path: Path) -> list[str]:
    from .polars_io import read_csv_df

    if not path.exists() or not path.is_file():
        return ["Metadata CSV does not exist."]
    if path.stat().st_size == 0:
        return ["Metadata CSV is empty (0 bytes)."]
    try:
        df = read_csv_df(path, nrows=0)
    except Exception as exc:
        return [f"Failed to read metadata CSV: {exc}"]
    if "Sample" not in df.columns:
        return ["Metadata CSV is missing required column: Sample"]
    return []


def _resolve_metadata_csv(
    run_id: str,
    output_dir: Path,
    *,
    force_prompt: bool = False,
) -> Optional[Path]:
    from .cache import load_run_cache, save_run_cache

    if not force_prompt:
        cached = load_run_cache(run_id)
        cached_path = cached.get("metadata_csv")
        if cached_path:
            path = Path(cached_path).expanduser()
            if not _validate_metadata_csv(path):
                return path

    while True:
        path_str, _ = QtWidgets.QFileDialog.getOpenFileName(
            None,
            "Select metadata CSV",
            str(output_dir),
            "CSV files (*.csv);;All files (*)",
        )
        if not path_str:
            return None
        path = Path(path_str).expanduser()
        errors = _validate_metadata_csv(path)
        if errors:
            QtWidgets.QMessageBox.critical(None, "Metadata error", "\n".join(errors))
            continue
        save_run_cache(run_id, path)
        return path


class _RunDataWorker(QtCore.QThread):
    """Loads a run's tables + metadata off the GUI thread.

    Reading the (potentially hundreds-of-MB) cached parquet takes several
    seconds; doing it here keeps the GUI event loop — and thus the loading
    splash animation — alive. Parquet/polars/pandas release the GIL during the
    heavy work, so the main thread keeps painting. The resulting DataFrames are
    handed back to the GUI thread via the `loaded` signal; the window itself is
    still built on the GUI thread (Qt requirement).
    """

    status = QtCore.pyqtSignal(str)
    loaded = QtCore.pyqtSignal(object)
    failed = QtCore.pyqtSignal(str, str)

    def __init__(self, run_id: str, output_dir: Path, metadata_csv: object, parent=None) -> None:
        super().__init__(parent)
        self._run_id = run_id
        self._output_dir = output_dir
        self._metadata_csv = metadata_csv

    def run(self) -> None:  # executes on the worker thread
        try:
            import pandas as pd

            from .cache import cached_tables_exist, read_cached_tables, write_cached_tables
            from .data import (
                attach_assembly_quality_metadata,
                load_output_data,
                prepare_metadata,
                validate_output_dir,
            )

            if cached_tables_exist(self._run_id, self._output_dir):
                self.status.emit("Reading cached data…")
                ann_df, ref_df, marker_hits_df = read_cached_tables(self._run_id)
            else:
                # Only when building from outputs do we need the deep per-file
                # validation (skipped in the dialog to keep "View" instant).
                self.status.emit("Validating run outputs…")
                _cfg, val_errors = validate_output_dir(self._output_dir, deep=True)
                if val_errors:
                    self.failed.emit("Run validation failed", "\n".join(val_errors))
                    return
                self.status.emit("Building data cache…")
                ann_df, ref_df, marker_hits_df = load_output_data(self._output_dir)
                write_cached_tables(self._run_id, self._output_dir, ann_df, ref_df, marker_hits_df)

            self.status.emit("Preparing metadata…")
            ref_samples = (
                ref_df["Sample"].astype(str).str.strip().replace({"": pd.NA}).dropna().unique().tolist()
                if not ref_df.empty and "Sample" in ref_df.columns
                else []
            )
            meta_df, warnings, _types = prepare_metadata(
                self._metadata_csv,
                ann_df["Sample"].unique().tolist(),
                ref_samples,
            )
            meta_df = attach_assembly_quality_metadata(meta_df, self._output_dir)
        except Exception as exc:  # surfaced on the GUI thread
            import traceback

            traceback.print_exc()
            self.failed.emit("Failed to load run", str(exc) or exc.__class__.__name__)
            return

        self.loaded.emit(
            {
                "ann_df": ann_df,
                "ref_df": ref_df,
                "marker_hits_df": marker_hits_df,
                "meta_df": meta_df,
                "warnings": warnings,
            }
        )


def main(*, dev_mode: bool = False) -> None:
    # Print a native stack trace if the process hard-crashes (segfault), which
    # a Python traceback can't capture.
    import faulthandler
    faulthandler.enable()
    # The Geo map defaults to the MapLibre GL web map (QtWebEngine), which needs
    # shared OpenGL contexts configured before the QApplication is created.
    # Harmless if the panel later falls back to the offline SVG map.
    QtCore.QCoreApplication.setAttribute(
        QtCore.Qt.ApplicationAttribute.AA_ShareOpenGLContexts,
        True,
    )
    app = QtWidgets.QApplication(sys.argv)
    if not _ensure_startup_ready(dev_mode=dev_mode):
        return
    (project_root() / "cache").mkdir(parents=True, exist_ok=True)
    (project_root() / "runs").mkdir(parents=True, exist_ok=True)
    app.setStyleSheet(
        """
        QWidget {
            color: #1f2c3b;
            selection-background-color: #e3edff;
            selection-color: #172434;
        }
        QMainWindow {
            background: #eef3f8;
        }
        QDialog {
            background: #f4f7fb;
        }
        QGroupBox {
            font-weight: 600;
            border: 1px solid #d7e0ea;
            border-radius: 12px;
            margin-top: 16px;
            background: #fbfcfe;
            color: #243447;
        }
        QGroupBox::title {
            subcontrol-origin: margin;
            subcontrol-position: top left;
            padding: 1px 8px;
            left: 12px;
            top: 0px;
            background: #eef3f8;
            color: #5a6b7c;
        }
        QLabel {
            color: #223246;
        }
        QLineEdit, QComboBox, QPlainTextEdit, QTextEdit {
            border: 1px solid #ccd6e2;
            border-radius: 8px;
            padding: 6px 8px;
            background: #ffffff;
            color: #172434;
        }
        QLineEdit:focus, QComboBox:focus, QPlainTextEdit:focus, QTextEdit:focus {
            border: 1px solid #4d86f7;
            background: #ffffff;
        }
        QPushButton {
            border: 1px solid #ccd6e2;
            border-radius: 8px;
            padding: 6px 12px;
            background: #ffffff;
            color: #203041;
            font-weight: 500;
        }
        QPushButton:hover {
            background: #f5f9ff;
            border-color: #9ab8f8;
        }
        QPushButton:pressed {
            background: #e4eeff;
            border-color: #7ea3ef;
        }
        QPushButton:disabled {
            background: #f5f7fa;
            border-color: #dde4ec;
            color: #93a1b2;
        }
        QTableView, QTreeView, QListWidget {
            background: #ffffff;
            alternate-background-color: #f8fbfd;
            gridline-color: #e6edf3;
            selection-background-color: #e3edff;
            selection-color: #172434;
            border: 1px solid #d7e0ea;
            border-radius: 12px;
            outline: none;
        }
        QHeaderView::section {
            background: #f4f7fb;
            color: #5a6b7c;
            padding: 8px 10px;
            border: none;
            border-right: 1px solid #e5ebf2;
            border-bottom: 1px solid #d7e0ea;
            font-weight: 600;
        }
        QTableCornerButton::section {
            background: #f4f7fb;
            border: none;
            border-right: 1px solid #e5ebf2;
            border-bottom: 1px solid #d7e0ea;
        }
        QTabWidget::pane {
            border: 1px solid #d7e0ea;
            border-radius: 12px;
            background: #ffffff;
            top: -1px;
        }
        QTabBar::tab {
            background: #f2f5f9;
            border: 1px solid #d7e0ea;
            color: #5b6c7d;
            padding: 8px 14px;
            border-top-left-radius: 10px;
            border-top-right-radius: 10px;
            margin-right: 4px;
        }
        QTabBar::tab:selected {
            background: #ffffff;
            color: #1f2c3b;
            border-bottom-color: #ffffff;
        }
        QTabBar::tab:hover:!selected {
            background: #f7faff;
            color: #32475c;
        }
        QProgressBar {
            border: 1px solid #d3dbe6;
            border-radius: 8px;
            background: #f3f6fa;
            color: #223246;
            text-align: center;
        }
        QProgressBar::chunk {
            border-radius: 8px;
            background: #7ea6f8;
        }
        QMenu {
            background: #ffffff;
            border: 1px solid #d7e0ea;
            padding: 6px;
        }
        QMenu::item {
            padding: 6px 18px;
            border-radius: 6px;
        }
        QMenu::item:selected {
            background: #ebf2ff;
        }
        QScrollBar:vertical {
            background: transparent;
            width: 12px;
            margin: 4px 0 4px 0;
        }
        QScrollBar::handle:vertical {
            background: #ccd6e2;
            border-radius: 6px;
            min-height: 28px;
        }
        QScrollBar::handle:vertical:hover {
            background: #b7c5d5;
        }
        QScrollBar:horizontal {
            background: transparent;
            height: 12px;
            margin: 0 4px 0 4px;
        }
        QScrollBar::handle:horizontal {
            background: #ccd6e2;
            border-radius: 6px;
            min-width: 28px;
        }
        QScrollBar::handle:horizontal:hover {
            background: #b7c5d5;
        }
        QScrollBar::add-line, QScrollBar::sub-line {
            width: 0px;
            height: 0px;
        }
        QStatusBar {
            background: #eef3f8;
            color: #607183;
            border-top: 1px solid #d7e0ea;
        }
        QToolTip {
            background: #1f2c3b;
            color: #ffffff;
            border: none;
            padding: 6px 8px;
        }
        """
    )
    # Make text slightly larger globally
    f = app.font()
    if f.pointSize() > 0:
        f.setPointSize(f.pointSize() + 1)
    app.setFont(f)
    dlg = StartupDialog(dev_mode=dev_mode)
    if dlg.exec() != QtWidgets.QDialog.DialogCode.Accepted:
        return
    if not dlg.output_dir or not dlg.run_id:
        return
    splash = LoadingSplash(status="Loading run…")
    splash.start()

    metadata_csv = _resolve_metadata_csv(
        dlg.run_id,
        dlg.output_dir,
        force_prompt=False,
    )
    if metadata_csv is None:
        splash.finish()
        return

    # Load the (multi-second) tables + metadata off the GUI thread so the splash
    # animation keeps running; build the window on the GUI thread once ready.
    state: dict = {}
    worker = _RunDataWorker(dlg.run_id, dlg.output_dir, metadata_csv)
    state["worker"] = worker  # keep a strong reference for the thread's lifetime

    def _on_status(text: str) -> None:
        splash.set_status(text)

    def _on_failed(title: str, message: str) -> None:
        splash.finish()
        QtWidgets.QMessageBox.critical(None, title, message)
        app.quit()

    def _on_loaded(payload: object) -> None:
        data = payload if isinstance(payload, dict) else {}
        warnings = data.get("warnings") or []
        if warnings:
            # Parented to the splash so this modal box stacks above it.
            QtWidgets.QMessageBox.warning(splash, "Metadata warnings", "\n".join(warnings))
        splash.set_status("Building views…")
        from .explorer_window import AnnotationExplorer

        try:
            win = AnnotationExplorer(
                data["ann_df"],
                data["meta_df"],
                ref_df=data["ref_df"],
                marker_hits_df=data["marker_hits_df"],
                output_dir=dlg.output_dir,
            )
        except BaseException as exc:  # never orphan the always-on-top splash
            splash.finish()
            QtWidgets.QMessageBox.critical(
                None, "Failed to build view", str(exc) or exc.__class__.__name__
            )
            app.quit()
            return
        state["win"] = win
        # Dismiss the splash the moment the first refresh settles (success or failure);
        # keep a generous hard fallback only in case that signal never arrives.
        win.firstRefreshFinished.connect(splash.finish)
        QtCore.QTimer.singleShot(60000, splash.finish)
        win.resize(1400, 900)
        win.show()
        QtCore.QTimer.singleShot(0, _start_world_geojson_worker)

    worker.status.connect(_on_status)
    worker.failed.connect(_on_failed)
    worker.loaded.connect(_on_loaded)
    worker.start()

    app.exec()


if __name__ == "__main__":
    main()
