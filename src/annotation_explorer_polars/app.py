from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Callable, Optional
from urllib import request

from PyQt6 import QtCore, QtWidgets

from .first_time_setup import ensure_app_ready
from .loading_splash import LoadingSplash
from .startup_dialog import StartupDialog
from .theme import apply_application_theme


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


def _ensure_startup_ready(*, dev_mode: bool = False) -> bool:
    """The install is ready to start: without a config, the choice between viewing
    shared runs and the pipeline setup; First-Time Setup when a full install lacks
    anything. False when the user cancels."""
    return ensure_app_ready(dev_mode=dev_mode)


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


def _choose_metadata(run_id: str, output_dir: Path, parent: Optional[QtWidgets.QWidget] = None) -> Optional[object]:
    """Ask for the run's metadata: a CSV (validated) or the run's own data
    (``RUN_METADATA_SOURCE``). The choice is remembered for the run; None
    when the user cancels."""
    from .cache import RUN_METADATA_SOURCE, save_run_cache

    while True:
        box = QtWidgets.QMessageBox(parent)
        box.setIcon(QtWidgets.QMessageBox.Icon.Question)
        box.setWindowTitle("Metadata")
        box.setText("Choose a metadata CSV for this run (it needs a Sample column).")
        box.setInformativeText(
            "Without one, the run's own data is used: species, MLST, fastGEAR clades and assembly QC."
        )
        csv_button = box.addButton("Select metadata CSV...", QtWidgets.QMessageBox.ButtonRole.AcceptRole)
        run_button = box.addButton("Use run data only", QtWidgets.QMessageBox.ButtonRole.AcceptRole)
        box.addButton(QtWidgets.QMessageBox.StandardButton.Cancel)
        box.setDefaultButton(csv_button)
        box.exec()
        clicked = box.clickedButton()
        if clicked is run_button:
            save_run_cache(run_id, None, run_data_only=True)
            return RUN_METADATA_SOURCE
        if clicked is not csv_button:
            return None
        while True:
            path_str, _ = QtWidgets.QFileDialog.getOpenFileName(
                parent,
                "Select metadata CSV",
                str(output_dir),
                "CSV files (*.csv);;All files (*)",
            )
            if not path_str:
                break  # back to the choice
            path = Path(path_str).expanduser()
            errors = _validate_metadata_csv(path)
            if errors:
                QtWidgets.QMessageBox.critical(parent, "Metadata error", "\n".join(errors))
                continue
            save_run_cache(run_id, path)
            return path


def _resolve_metadata_csv(run_id: str, output_dir: Path) -> Optional[object]:
    """The run's remembered metadata (a valid CSV, or its own data only),
    else the user's choice; None when the user cancels."""
    from .cache import RUN_METADATA_SOURCE, load_run_cache

    cached = load_run_cache(run_id)
    if cached.get("metadata_source") == RUN_METADATA_SOURCE:
        return RUN_METADATA_SOURCE
    cached_path = cached.get("metadata_csv")
    if cached_path:
        path = Path(cached_path).expanduser()
        if not _validate_metadata_csv(path):
            return path
    shared = _shared_metadata(run_id, output_dir)
    if shared is not None:
        return shared
    return _choose_metadata(run_id, output_dir)


def _shared_metadata(run_id: str, output_dir: Path) -> Optional[object]:
    """A shared run first opens with the metadata it was shared with: its CSV
    (``share/metadata.csv``), or the run's own data; remembered like a choice."""
    from pipeline.RunShare import METADATA_CSV, is_shared_run, read_manifest

    from .cache import RUN_METADATA_SOURCE, save_run_cache

    if not is_shared_run(output_dir):
        return None
    csv_path = output_dir / METADATA_CSV
    if csv_path.is_file() and not _validate_metadata_csv(csv_path):
        save_run_cache(run_id, csv_path)
        return csv_path
    if read_manifest(output_dir).get("metadata") == "run":
        save_run_cache(run_id, None, run_data_only=True)
        return RUN_METADATA_SOURCE
    return None


SHARED_TABLES_UNREADABLE = "This shared run has no data tables this version can read."


class _RunDataWorker(QtCore.QThread):
    """Loads a run's tables + metadata and builds its RunData off the GUI thread.

    Reading the (potentially hundreds-of-MB) cached parquet and deriving the
    run's indexes and counts (``RunData.build``) take several seconds; doing
    it here keeps the GUI event loop, and thus the loading splash animation,
    alive. The RunData is handed back to the GUI thread via the `loaded`
    signal; the window itself is still built on the GUI thread (Qt
    requirement).
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

            from pipeline.RunShare import is_shared_run

            from .cache import (
                RUN_METADATA_SOURCE,
                cached_tables_exist,
                read_cached_tables,
                shared_table_paths,
                shared_tables_status,
                write_cached_tables,
            )
            from .data import (
                attach_assembly_quality_metadata,
                load_output_data,
                prepare_metadata,
                prepare_run_metadata,
                validate_output_dir,
            )
            from .model import RunData

            if is_shared_run(self._output_dir):
                # A shared run has nothing to build its tables from: read the ones
                # it came with (also in another format; the run card warned).
                if not shared_tables_status(self._output_dir)[0]:
                    self.failed.emit("Failed to load run", SHARED_TABLES_UNREADABLE)
                    return
                self.status.emit("Reading shared data…")
                try:
                    ann_df, ref_df, marker_hits_df = read_cached_tables(
                        self._run_id, shared_table_paths(self._output_dir)
                    )
                except Exception:
                    self.failed.emit("Failed to load run", SHARED_TABLES_UNREADABLE)
                    return
            elif cached_tables_exist(self._run_id, self._output_dir):
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
            if self._metadata_csv == RUN_METADATA_SOURCE:
                meta_df, warnings, _types = prepare_run_metadata(
                    self._output_dir,
                    ann_df["Sample"].unique().tolist(),
                    ref_samples,
                )
            else:
                meta_df, warnings, _types = prepare_metadata(
                    self._metadata_csv,
                    ann_df["Sample"].unique().tolist(),
                    ref_samples,
                )
            meta_df = attach_assembly_quality_metadata(meta_df, self._output_dir)
            run_data = RunData.build(
                ann_df,
                meta_df,
                ref_df=ref_df,
                marker_hits_df=marker_hits_df,
                output_dir=self._output_dir,
                run_id=self._run_id,
                status=self.status.emit,
            )
            # Pangenome folder lookups here, not on the UI thread when the
            # window checks which pan panels are available (Phase 9).
            self.status.emit("Finding pangenome outputs…")
            run_data.pan.resolve_lookups()
        except Exception as exc:  # surfaced on the GUI thread
            import traceback

            traceback.print_exc()
            self.failed.emit("Failed to load run", str(exc) or exc.__class__.__name__)
            return

        self.loaded.emit({"run_data": run_data, "warnings": warnings})


def _docking_ready() -> bool:
    """Explore windows need PyQt6Ads (there is no fallback layout): without
    it, say so with the install command instead of failing later."""
    from .docking import (
        DOCKING_UNAVAILABLE_HEADLINE,
        docking_available,
        docking_install_hint,
        docking_unavailable_message,
    )

    if docking_available():
        return True
    print(docking_unavailable_message(), file=sys.stderr)
    box = QtWidgets.QMessageBox(
        QtWidgets.QMessageBox.Icon.Critical,
        "PyQt6Ads is missing",
        DOCKING_UNAVAILABLE_HEADLINE,
    )
    box.setInformativeText(docking_install_hint())
    box.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
    box.exec()
    return False


def _open_loaded_run(
    payload: object,
    *,
    app: QtWidgets.QApplication,
    splash: Any,
    state: dict,
    open_run_list: Optional[Callable[[Any, QtWidgets.QWidget], None]] = None,
) -> None:
    """The loader finished: build the run's AppController and its first
    explore window. ``state["app"]`` keeps the controller alive; the app quits
    when its last explore window closes."""
    data = payload if isinstance(payload, dict) else {}
    warnings = data.get("warnings") or []
    if warnings:
        # Parented to the splash so this modal box stacks above it.
        QtWidgets.QMessageBox.warning(splash, "Metadata warnings", "\n".join(warnings))
    splash.set_status("Building views…")
    from .shell.app_controller import AppController

    try:
        # The application object owns the run's data, the windows and the
        # commands; the app quits when its last explore window closes.
        controller = AppController(data["run_data"], manage_quit=True, open_run_list=open_run_list)
        win = controller.open_first_window()
    except BaseException as exc:  # never orphan the always-on-top splash
        splash.finish()
        QtWidgets.QMessageBox.critical(
            None, "Failed to build view", str(exc) or exc.__class__.__name__
        )
        app.exit(1)  # nothing to show: a nonzero exit status
        return
    reloaded = "app" in state
    state["app"] = controller
    # Dismiss the splash the moment the first refresh settles (success or failure);
    # keep a generous hard fallback only in case that signal never arrives.
    win.firstRefreshFinished.connect(splash.finish)
    QtCore.QTimer.singleShot(60000, splash.finish)
    if not reloaded:
        QtCore.QTimer.singleShot(0, _start_world_geojson_worker)


def _load_run(app: QtWidgets.QApplication, run_id: str, output_dir: Path, metadata_csv: object, *, state: dict) -> None:
    """Load the run behind the splash and open its windows. File > Run List...
    closes the run's windows (their session is recorded) and shows the run list
    again; Cancel there quits."""
    splash = LoadingSplash(status="Loading run…")
    splash.start()

    # Load the (multi-second) tables + metadata off the GUI thread so the splash
    # animation keeps running; build the window on the GUI thread once ready.
    worker = _RunDataWorker(run_id, output_dir, metadata_csv)
    state["worker"] = worker  # keep a strong reference for the thread's lifetime

    def _on_status(text: str) -> None:
        splash.set_status(text)

    def _on_failed(title: str, message: str) -> None:
        splash.finish()
        QtWidgets.QMessageBox.critical(None, title, message)
        # Back to the run list; Cancel there quits with a nonzero exit status.
        state["load_failed"] = True
        QtCore.QTimer.singleShot(0, lambda: _show_run_list(app, state))

    def _open_run_list(controller: Any, _window: QtWidgets.QWidget) -> None:
        if controller.close_for_reload():
            # After the menu's slot returns: the window it belongs to is closed.
            QtCore.QTimer.singleShot(0, lambda: _show_run_list(app, state))

    def _on_loaded(payload: object) -> None:
        _open_loaded_run(payload, app=app, splash=splash, state=state, open_run_list=_open_run_list)

    worker.status.connect(_on_status)
    worker.failed.connect(_on_failed)
    worker.loaded.connect(_on_loaded)
    worker.start()


def _choose_run(*, dev_mode: bool = False) -> Optional[tuple[str, Path, object]]:
    """The run list (everything about runs lives there), then the run's
    metadata; None when the user cancels. Cancel in the metadata choice goes
    back to the run list."""
    while True:
        dlg = StartupDialog(dev_mode=dev_mode)
        if dlg.exec() != QtWidgets.QDialog.DialogCode.Accepted:
            return None
        if not dlg.output_dir or not dlg.run_id:
            return None
        metadata_csv = _resolve_metadata_csv(dlg.run_id, dlg.output_dir)
        if metadata_csv is not None:
            return dlg.run_id, dlg.output_dir, metadata_csv


def _show_run_list(app: QtWidgets.QApplication, state: dict) -> None:
    """File > Run List...: the run's windows are closed; choose a run (Cancel quits)."""
    choice = _choose_run(dev_mode=bool(state.get("dev_mode")))
    failed = state.pop("load_failed", False)
    if choice is None:
        app.exit(1 if failed else 0)  # 1: the run could not load
        return
    _load_run(app, *choice, state=state)


def main(*, dev_mode: bool = False) -> int:
    """Run the app; the exit status for the entry points (1 when it cannot
    start because PyQt6Ads is missing, 0 when the user cancels startup)."""
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
    # Identity, palette, font, and QSS must be in place before readiness checks:
    # those checks can open the first-time-setup and tool-image dialogs.
    apply_application_theme(app)
    if not _docking_ready():
        return 1
    if not _ensure_startup_ready(dev_mode=dev_mode):
        return 0
    (project_root() / "cache").mkdir(parents=True, exist_ok=True)
    (project_root() / "runs").mkdir(parents=True, exist_ok=True)
    state: dict = {"dev_mode": dev_mode}
    choice = _choose_run(dev_mode=dev_mode)
    if choice is None:
        return 0
    _load_run(app, *choice, state=state)
    return int(app.exec())


if __name__ == "__main__":
    sys.exit(main())
