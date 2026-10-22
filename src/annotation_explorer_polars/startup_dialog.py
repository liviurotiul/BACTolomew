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
from typing import TYPE_CHECKING, Callable, List, Optional, Sequence

from PyQt6 import QtCore, QtWidgets

from .run_progress import RunProgressDialog, pending_text_for, status_bar_state
from .constants import RUN_CONFIG_FILENAME, RUN_LOCK_FILENAME, RUN_STOP_FILENAME
from .run_expectations import (
    classify_expected_path,
    describe_run_status,
    expected_output_exists,
    load_run_state,
    load_run_status,
    resolve_expected_files,
    run_allows_resume,
    run_allows_view,
    run_is_active,
)
from .run_steps import RunSteps, summarize_run_steps

if TYPE_CHECKING:
    from pipeline.DataInput import DataInput


@dataclass
class RunEntry:
    run_id: str
    sample_set_id: str
    sample_count: Optional[int]
    path: Path
    run_config: Optional[dict]
    steps: RunSteps
    running: bool
    status: str
    status_text: Optional[str]
    can_view: bool
    can_resume: bool
    pending_text: Optional[str] = None
    status_detail: Optional[str] = None
    # Imported with Import Run... (a share package or a copied run folder).
    shared: bool = False
    # Compatibility findings (pipeline.RunShare.Issue), warnings first.
    issues: tuple = ()
    metadata_label: str = ""
    version: str = ""

    @property
    def warnings(self) -> list:
        return [issue for issue in self.issues if issue.level == "warning"]


def _viewer_mode() -> bool:
    """A viewer-only install (no tools image or databases)."""
    from pipeline.ToolEnv import viewer_mode

    return viewer_mode()


PIPELINE_SETUP_NEEDED = "Needs the pipeline setup"
IMPORTED_VIEW_ONLY = "Imported runs can be viewed only"


def _metadata_label(run_id: str, output_dir: Path) -> str:
    """The card's metadata line, from the run's remembered choice."""
    from pipeline.RunShare import METADATA_CSV, is_shared_run, read_manifest

    from .cache import RUN_METADATA_SOURCE, load_run_cache

    cached = load_run_cache(run_id, create=False)
    if cached.get("metadata_source") == RUN_METADATA_SOURCE:
        return "Metadata: run data only"
    if cached.get("metadata_csv"):
        return f"Metadata: {Path(cached['metadata_csv']).name}"
    if is_shared_run(output_dir):
        # Opened first with what it was shared with (app._shared_metadata).
        if (output_dir / METADATA_CSV).is_file():
            return "Metadata: the shared CSV"
        if read_manifest(output_dir).get("metadata") == "run":
            return "Metadata: run data only"
    return "Metadata: not chosen"


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
    from pipeline.Setup import collect_setup_missing_paths

    from .first_time_setup import FirstTimeSetupDialog

    root = _project_root()
    config_path = root / "bactolomew_config.json"
    config = _load_app_config()

    missing: list[Path] = []
    if not config_path.exists():
        missing.append(config_path)
    missing.extend(collect_setup_missing_paths(root, config, dev_mode=dev_mode))

    if not missing:
        return config

    dlg = FirstTimeSetupDialog(parent=parent, dev_mode=dev_mode)
    if dlg.exec() != QtWidgets.QDialog.DialogCode.Accepted:
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


def _load_catalog_entries(path: Path) -> list[tuple[Path, str]] | None:
    """The catalog's (path, category) entries, or None when it cannot be read."""
    from pipeline.DbCatalog import load_db_catalog

    try:
        catalog = load_db_catalog(path)
    except Exception as exc:
        QtWidgets.QMessageBox.critical(None, "Catalog error", f"Failed to load catalog: {path}\n{exc}")
        return None
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


# "Not now" before a new run: not asked again before new runs in this session (a run that
# then needs a download stops at its chewBBACA step and its Resume asks again).
_TERMS_DECLINED_FOR_NEW_RUNS = False


def ask_cgmlst_terms_before_run(
    parent: QtWidgets.QWidget, output_dir: Optional[Path], chewbbaca_config: Optional[dict]
) -> None:
    """Before a run (``output_dir`` None for a new one) whose chewBBACA step may install its
    species' scheme from cgMLST.org: the cgMLST.org terms, unless they are accepted ("Not
    now" starts the run anyway; it stops at the chewBBACA step if it needs a download).
    A new run's species is not known yet, so after "Not now" later new runs of this session
    are not asked again; a Resume that needs a download is."""
    global _TERMS_DECLINED_FOR_NEW_RUNS
    from pipeline.ChewBBACAWorkflow import cgmlst_terms_accepted, run_may_install_schema
    from pipeline.CgmlstInstall import default_dest

    from .cgmlst_download_dialog import confirm_cgmlst_terms

    if cgmlst_terms_accepted(config=chewbbaca_config):
        return
    if output_dir is None and _TERMS_DECLINED_FOR_NEW_RUNS:
        return
    if output_dir is not None:
        steps = (load_run_state(output_dir) or {}).get("steps")
        step = steps.get("chewbbaca") if isinstance(steps, dict) else None
        if isinstance(step, dict) and step.get("status") == "success":
            return
    if run_may_install_schema(output_dir, chewbbaca_config):
        accepted = confirm_cgmlst_terms(parent, default_dest(), before_run=True)
        if output_dir is None and not accepted:
            _TERMS_DECLINED_FOR_NEW_RUNS = True


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


def _cgmlst_stop_text(status_info: dict[str, object]) -> str:
    """The card line of a run stopped at its chewBBACA step for want of a cgMLST scheme
    ("" for any other failure)."""
    from pipeline.ChewBBACAWorkflow import CGMLST_NO_SCHEMA_FOR_SPECIES, CGMLST_TERMS_NOT_ACCEPTED

    state = status_info.get("state")
    steps = state.get("steps") if isinstance(state, dict) else None
    step = steps.get("chewbbaca") if isinstance(steps, dict) else None
    code = step.get("code") if isinstance(step, dict) and step.get("status") == "failed" else ""
    if code == CGMLST_NO_SCHEMA_FOR_SPECIES:
        return "Failed: no cgMLST scheme. Import one in Settings > Add schema, then Resume."
    if code == CGMLST_TERMS_NOT_ACCEPTED:
        return "Failed: cgMLST.org terms not accepted. Click Resume to accept them."
    return ""


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
        return _cgmlst_stop_text(status_info) or "Failed. Click Resume to retry."
    if status == "partial":
        if pending_text:
            return pending_text
        return "Partial run. Click Resume to continue."
    return full_status or pending_text or "Status unknown. Click Resume to continue."


def _repolish(widget: QtWidgets.QWidget) -> None:
    widget.style().unpolish(widget)
    widget.style().polish(widget)


# Vertical gap between two run cards in the start window's list.
_CARD_GAP = 8


class _ElidedLabel(QtWidgets.QLabel):
    """A one-line label that elides its text in the middle when it does not fit;
    the full text is then its tooltip."""

    def __init__(self, text: str, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(text, parent)
        self._full_text = text

    def sizeHint(self) -> QtCore.QSize:  # type: ignore[override]
        hint = super().sizeHint()
        return QtCore.QSize(self.fontMetrics().horizontalAdvance(self._full_text) + 2, hint.height())

    def minimumSizeHint(self) -> QtCore.QSize:  # type: ignore[override]
        return QtCore.QSize(40, super().minimumSizeHint().height())

    def resizeEvent(self, event) -> None:  # type: ignore[override]
        super().resizeEvent(event)
        shown = self.fontMetrics().elidedText(self._full_text, QtCore.Qt.TextElideMode.ElideMiddle, self.width())
        if shown != self.text():
            self.setText(shown)
        self.setToolTip(self._full_text if shown != self._full_text else "")


class RunCard(QtWidgets.QFrame):
    """One run in the start window: name, samples, step progress, status, actions.
    Styled by the theme (QFrame#RunCard, stepBar)."""

    def __init__(
        self,
        entry: RunEntry,
        on_view: Optional[Callable[[RunEntry], None]] = None,
        on_monitor: Optional[Callable[[RunEntry], None]] = None,
        on_resume: Optional[Callable[[RunEntry], None]] = None,
        on_stop: Optional[Callable[[RunEntry], None]] = None,
        parent: Optional[QtWidgets.QWidget] = None,
        *,
        on_more: Optional[Callable[[RunEntry, str], None]] = None,
    ) -> None:
        super().__init__(parent)
        self.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.setObjectName("RunCard")
        self.setProperty("selected", False)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(16, 12, 16, 12)
        layout.setSpacing(6)

        title_row = QtWidgets.QHBoxLayout()
        title_row.setSpacing(8)
        title = _ElidedLabel(entry.sample_set_id or entry.run_id)
        title.setProperty("role", "runTitle")
        title_row.addWidget(title, 1)
        if entry.shared:
            chip = QtWidgets.QLabel("Shared")
            chip.setProperty("role", "runChip")
            chip.setToolTip("Imported with Import Run...")
            title_row.addWidget(chip, 0, QtCore.Qt.AlignmentFlag.AlignVCenter)
        layout.addLayout(title_row)

        attempt_label = entry.run_id
        if entry.sample_count is None:
            sample_text = "Samples: unknown"
        else:
            sample_text = f"{entry.sample_count} samples"
        subtitle_bits = [attempt_label, sample_text]
        subtitle = QtWidgets.QLabel(" • ".join(bit for bit in subtitle_bits if bit))
        subtitle.setProperty("role", "runMeta")
        layout.addWidget(subtitle)
        if entry.metadata_label:
            metadata = QtWidgets.QLabel(entry.metadata_label)
            metadata.setProperty("role", "runMeta")
            layout.addWidget(metadata)
        from pipeline.RunShare import warning_summary

        warning_text = warning_summary(entry.issues)
        if warning_text:
            warning = QtWidgets.QLabel(warning_text)
            warning.setProperty("role", "runWarning")
            warning.setWordWrap(True)
            warning.setToolTip("\n".join(issue.text for issue in entry.warnings))
            layout.addWidget(warning)

        steps = entry.steps
        progress = QtWidgets.QProgressBar()
        progress.setProperty("role", "stepBar")
        progress.setRange(0, 1000)
        fraction = steps.fraction
        if not steps.rows and entry.can_view:
            fraction = 1.0
        progress.setValue(int(round(fraction * 1000)))
        progress.setTextVisible(False)
        progress.setProperty("state", "running" if entry.running else status_bar_state(entry.status))
        layout.addSpacing(4)
        layout.addWidget(progress)
        layout.addSpacing(2)

        # While running: "Step 9 of 15 · Gubbins · 14 min"; otherwise the status sentence.
        if entry.running:
            line_text = steps.step_line() or entry.status_text or "Running"
        else:
            line_text = entry.status_text or ""
        info_row = QtWidgets.QHBoxLayout()
        info_row.setSpacing(10)
        info = QtWidgets.QLabel(line_text)
        info.setProperty("role", "runStep")
        info.setWordWrap(True)
        info.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.NoTextInteraction)
        detail = entry.status_detail or entry.status_text
        if detail and detail != line_text:
            info.setToolTip(detail)
        info_row.addWidget(info, 1)
        if steps.rows and not entry.running:
            count = QtWidgets.QLabel(f"{steps.finished} / {steps.total} steps")
            count.setProperty("role", "runCount")
            info_row.addWidget(count, 0, QtCore.Qt.AlignmentFlag.AlignTop)
        layout.addLayout(info_row)

        btn_row = QtWidgets.QHBoxLayout()
        btn_row.setSpacing(6)
        btn_row.addStretch(1)
        view_btn = QtWidgets.QPushButton("View")
        view_btn.setEnabled(entry.can_view)
        view_btn.setProperty("role", "primary")
        view_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        if not entry.can_view:
            view_btn.setToolTip(
                entry.status_detail
                or entry.status_text
                or entry.pending_text
                or "This run is not ready to open."
            )
        view_btn.clicked.connect(lambda: on_view(entry) if on_view else None)

        if not entry.can_view:
            monitor_btn = QtWidgets.QPushButton("Monitor")
            monitor_btn.clicked.connect(lambda: on_monitor(entry) if on_monitor else None)
            resume_btn = QtWidgets.QPushButton("Resume")
            resume_btn.setEnabled(entry.can_resume)
            if entry.shared:
                resume_btn.setToolTip(IMPORTED_VIEW_ONLY)
            elif not entry.can_resume and not entry.running and _viewer_mode():
                resume_btn.setToolTip(PIPELINE_SETUP_NEEDED)
            resume_btn.clicked.connect(lambda: on_resume(entry) if on_resume else None)
            stop_btn = QtWidgets.QPushButton("Stop")
            stop_btn.setEnabled(entry.running)
            stop_btn.clicked.connect(lambda: on_stop(entry) if on_stop else None)
            force_btn = QtWidgets.QPushButton("Force Stop")
            force_btn.setEnabled(entry.running)
            force_btn.setProperty("role", "danger")
            force_btn.clicked.connect(lambda: on_stop(entry, True) if on_stop else None)
            for _b in (monitor_btn, resume_btn, stop_btn, force_btn):
                _b.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
                btn_row.addWidget(_b)
        btn_row.addWidget(self._more_button(entry, on_more))
        btn_row.addWidget(view_btn)
        layout.addSpacing(4)
        layout.addLayout(btn_row)

    def _more_button(
        self, entry: RunEntry, on_more: Optional[Callable[[RunEntry, str], None]]
    ) -> QtWidgets.QPushButton:
        """More: everything else about the run (details, metadata, sharing, its folder)."""
        button = QtWidgets.QPushButton("More")
        button.setObjectName("RunMoreButton")
        button.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        menu = QtWidgets.QMenu(button)
        items = [
            ("details", "Run Details...", True, "Versions, tools, databases and what may not match"),
            ("metadata", "Change Metadata...", not entry.running, "Choose this run's metadata CSV or use the run's own data"),
            ("export", "Export for Sharing...", entry.can_view and entry.status == "success" and not entry.running,
             "One file with what the explorer needs to show this run"),
            ("folder", "Show in Folder", True, "Open the run's folder"),
        ]
        if entry.shared:
            items.append(("remove", "Remove...", not entry.running, "Delete this imported run from the run list"))
        for key, text, enabled, tip in items:
            action = menu.addAction(text)
            action.setEnabled(bool(enabled))
            action.setToolTip(tip)
            action.setObjectName(f"Run{key.capitalize()}Action")
            action.triggered.connect(lambda _c=False, k=key: on_more(entry, k) if on_more else None)
        menu.setToolTipsVisible(True)
        button.setMenu(menu)
        return button

    def set_selected(self, selected: bool) -> None:
        if bool(self.property("selected")) == bool(selected):
            return
        self.setProperty("selected", bool(selected))
        _repolish(self)


def _image_population_methods() -> tuple:
    """The population structure methods of the configured tools image (empty when it
    cannot be checked); cached per image id by pipeline.PopulationStructure, so only the
    first New Run window of a session waits for the image check (in a worker thread)."""
    try:
        from pipeline.PopulationStructure import image_population_methods

        return image_population_methods() or ()
    except Exception:
        return ()


class _PopulationProbe(QtCore.QObject):
    """Runs _image_population_methods in a worker thread; `done` arrives on the GUI thread."""

    done = QtCore.pyqtSignal(object)

    def start(self) -> None:
        import threading

        threading.Thread(target=self._run, daemon=True, name="new-run-image-check").start()

    def _run(self) -> None:
        methods = _image_population_methods()
        try:
            self.done.emit(methods)
        except RuntimeError:  # the window was closed first
            pass


class NewRunDialog(QtWidgets.QDialog):
    def __init__(
        self, parent: Optional[QtWidgets.QWidget] = None, *, population_methods: Optional[Sequence[str]] = None
    ) -> None:
        super().__init__(parent)
        self.genome_paths: List[Path] = []
        # The method picked in the dialog (None: the config's phylo.population_method).
        self.population_method: Optional[str] = None
        self.population_combo: Optional[QtWidgets.QComboBox] = None

        self.setObjectName("NewRunDialog")
        self.setWindowTitle("Start New Run")
        self.resize(720, 520)
        self._build_ui()
        # The tools image's population structure methods (a choice only when it has both):
        # given, or checked in a worker thread so the window opens at once.
        if population_methods is not None:
            self.set_population_methods(population_methods)
        else:
            self._probe = _PopulationProbe(self)
            self._probe.done.connect(self.set_population_methods)
            self._probe.start()

    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(20, 18, 20, 14)
        layout.setSpacing(10)

        intro = QtWidgets.QLabel(
            "Select genomes to run the full required pipeline."
        )
        intro.setWordWrap(True)
        intro.setProperty("role", "muted")
        layout.addWidget(intro)

        genomes_group = QtWidgets.QGroupBox("Assemblies")
        genomes_layout = QtWidgets.QVBoxLayout(genomes_group)
        self.genome_list = QtWidgets.QListWidget()
        genomes_layout.addWidget(self.genome_list)

        btn_row = QtWidgets.QHBoxLayout()
        add_btn = QtWidgets.QPushButton("Add Assemblies…")
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
        self._population_box = QtWidgets.QWidget()
        self._population_box.setObjectName("PopulationMethodBox")
        self._population_layout = QtWidgets.QVBoxLayout(self._population_box)
        self._population_layout.setContentsMargins(0, 0, 0, 0)
        self._population_layout.setSpacing(4)
        layout.addWidget(self._population_box)

        btns = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel
        )
        btns.button(QtWidgets.QDialogButtonBox.StandardButton.Ok).setProperty("role", "primary")
        btns.accepted.connect(self._validate_and_accept)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)

    def set_population_methods(self, available: Sequence[str]) -> None:
        """A "Population structure" dropdown when the tools image has fastGEAR and FastBAPS
        (default: phylo.population_method, auto meaning FastBAPS); with one method a muted
        line names it and the run uses it, whatever the configured default; nothing when the
        image could not be checked."""
        from pipeline.PopulationStructure import AUTO_PREFERENCE, DEFAULT_METHOD, METHOD_LABELS
        from pipeline.ToolEnv import configured_phylo_settings

        from .settings_column import settings_combo

        layout = self._population_layout
        methods = [m for m in AUTO_PREFERENCE if m in tuple(available or ())]
        if len(methods) == 1:
            self.population_method = methods[0]
            note = QtWidgets.QLabel(f"Population structure: {METHOD_LABELS[methods[0]]} (the only method in the tools image).")
            note.setProperty("role", "muted")
            note.setObjectName("PopulationMethodNote")
            layout.addWidget(note)
            return
        if len(methods) < 2:
            return
        try:
            configured = str(configured_phylo_settings().get("population_method") or DEFAULT_METHOD)
        except ValueError:
            configured = DEFAULT_METHOD
        self.population_method = configured if configured in methods else methods[0]
        tooltips = {
            "fastgear": "fastGEAR lineages and clusters.",
            "fastbaps": "FastBAPS clusters: lineage = Level 1, cluster = Level 2.",
        }
        caption = QtWidgets.QLabel("Population structure")
        caption.setProperty("role", "settingsCaption")
        layout.addWidget(caption)
        combo = settings_combo(
            self,
            [(m, METHOD_LABELS[m], tooltips[m]) for m in methods],
            self.population_method,
            self._set_population_method,
        )
        combo.setObjectName("PopulationMethodCombo")
        row = QtWidgets.QHBoxLayout()
        row.addWidget(combo)
        row.addStretch(1)
        layout.addLayout(row)
        self.population_combo = combo

    def _set_population_method(self, method: str) -> None:
        self.population_method = str(method)

    def _add_genomes(self) -> None:
        paths, _ = QtWidgets.QFileDialog.getOpenFileNames(
            self,
            "Select assembly FASTA files",
            "",
            "FASTA files (*.fa *.fasta *.fna *.fa.gz *.fasta.gz *.fna.gz);;All files (*)",
        )
        if not paths:
            return
        self._add_genome_paths([Path(p) for p in paths])

    def _add_genome_paths(self, paths: List[Path]) -> List[Path]:
        """Add genome files, refusing any whose sample name is already taken.

        Names come from the resolved file (symlinks followed), as the pipeline names
        samples, so the picker and DataInput agree. Returns the rejected paths; a
        warning lists each one with the file that already holds its name.
        """
        from pipeline.DataInput import genome_name_keys

        known: set[str] = set()
        taken: dict[str, Path] = {}
        for existing in self.genome_paths:
            ident, keys = genome_name_keys(existing)
            known.add(ident)
            for key in keys:
                taken.setdefault(key, existing)

        rejected: List[tuple[Path, Optional[Path]]] = []
        for path in paths:
            ident, keys = genome_name_keys(path)
            if path in self.genome_paths or ident in known:
                continue
            clash = next((taken[k] for k in sorted(keys) if k in taken), None)
            if clash is not None:
                rejected.append((path, clash))
                continue
            known.add(ident)
            for key in keys:
                taken.setdefault(key, path)
            self.genome_paths.append(path)
            self.genome_list.addItem(str(path))
        if rejected:
            lines = [
                f"{path}\n    same sample name as {other}" if other else str(path)
                for path, other in rejected
            ]
            QtWidgets.QMessageBox.warning(
                self,
                "Duplicate genome names",
                "Each genome needs a unique sample name (the file name without .gz and "
                "the FASTA suffix; for a link, the name of the file it points to). "
                "These files were not added:\n\n" + "\n".join(lines),
            )
        return [path for path, _ in rejected]

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
        from pipeline.DataInput import find_genome_name_collisions, format_genome_name_collisions

        collisions = find_genome_name_collisions(self.genome_paths)
        if collisions:
            QtWidgets.QMessageBox.critical(self, "Validation error", format_genome_name_collisions(collisions))
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
        # A partial run's explanation checks its expected files: once per run state.
        self._pending_cache: dict[str, tuple[tuple[str, str], Optional[str]]] = {}
        from .run_share_dialogs import compatibility_context

        self._compat = compatibility_context()

        self.setWindowTitle("Select Run")
        self.resize(720, 600)
        self.setMinimumWidth(600)  # a card's five buttons in one row
        self.setObjectName("StartupDialog")
        self._build_ui()
        self._load_runs()

        self._refresh_timer = QtCore.QTimer(self)
        self._refresh_timer.setInterval(5000)
        self._refresh_timer.timeout.connect(self._load_runs)
        self._refresh_timer.start()

    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(22, 20, 22, 14)
        layout.setSpacing(4)

        title = QtWidgets.QLabel("Select a run")
        title.setProperty("role", "windowTitle")
        layout.addWidget(title)
        intro = QtWidgets.QLabel("Open a completed run in the annotation explorer, or start a new one.")
        intro.setWordWrap(True)
        intro.setProperty("role", "windowHint")
        layout.addWidget(intro)
        layout.addSpacing(10)

        self.run_list = QtWidgets.QListWidget()
        self.run_list.setObjectName("RunList")
        self.run_list.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.run_list.setSpacing(0)  # cards line up with the header and footer; _CARD_GAP below each
        self.run_list.setVerticalScrollMode(QtWidgets.QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.run_list.verticalScrollBar().setSingleStep(16)
        self.run_list.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.run_list.viewport().installEventFilter(self)
        self.run_list.currentItemChanged.connect(self._on_run_selected)
        layout.addWidget(self.run_list, 1)

        footer = QtWidgets.QFrame()
        footer.setObjectName("StartFooter")
        btn_row = QtWidgets.QHBoxLayout(footer)
        btn_row.setContentsMargins(0, 12, 0, 0)
        btn_row.setSpacing(8)
        settings_btn = QtWidgets.QPushButton("⚙")
        settings_btn.setToolTip("Settings")
        settings_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        settings_btn.clicked.connect(self._open_settings_dialog)
        btn_row.addWidget(settings_btn)
        new_run_btn = QtWidgets.QPushButton("Start New Run")
        new_run_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        new_run_btn.clicked.connect(self._open_new_run_dialog)
        btn_row.addWidget(new_run_btn)
        import_btn = QtWidgets.QPushButton("Import Run...")
        import_btn.setObjectName("ImportRunButton")
        import_btn.setToolTip("Add a run someone shared with you (.zip) or a run folder")
        import_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        import_btn.clicked.connect(self._import_run)
        btn_row.addWidget(import_btn)
        btn_row.addStretch(1)
        cancel_btn = QtWidgets.QPushButton("Cancel")
        cancel_btn.setProperty("role", "flat")
        cancel_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        cancel_btn.clicked.connect(self.reject)
        btn_row.addWidget(cancel_btn)
        layout.addWidget(footer)

    def eventFilter(self, watched: QtCore.QObject, event: QtCore.QEvent) -> bool:  # type: ignore[override]
        if watched is self.run_list.viewport() and event.type() == QtCore.QEvent.Type.Resize:
            self._fit_cards()
        return super().eventFilter(watched, event)

    def _fit_cards(self) -> None:
        """Cards span the list's width; their height follows their wrapped text."""
        width = max(240, self.run_list.viewport().width())
        for i in range(self.run_list.count()):
            item = self.run_list.item(i)
            card = self._card(item)
            if card is None:
                continue
            card.setFixedWidth(width)
            height = card.heightForWidth(width) if card.hasHeightForWidth() else -1
            if height <= 0:
                height = card.sizeHint().height()
            item.setSizeHint(QtCore.QSize(width, height + _CARD_GAP))

    def _card(self, item: Optional[QtWidgets.QListWidgetItem]) -> Optional["RunCard"]:
        holder = self.run_list.itemWidget(item) if item is not None else None
        return holder.findChild(RunCard) if holder is not None else None

    def _pending_text(self, output_dir: Path, status_info: dict[str, object]) -> Optional[str]:
        status = str(status_info.get("status") or "unknown")
        key = (status, str(status_info.get("updated_at") or ""))
        cached = self._pending_cache.get(str(output_dir))
        if cached is not None and cached[0] == key:
            return cached[1]
        text = pending_text_for(status, output_dir)
        self._pending_cache[str(output_dir)] = (key, text)
        return text

    def _load_runs(self) -> None:
        if QtWidgets.QApplication.activePopupWidget() is not None:
            return  # a card's More menu is open: refresh on the next tick
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

        from pipeline.RunShare import STAGING_PREFIX, is_imported_run, is_shared_run

        run_cfg_paths = sorted(runs_dir.rglob(RUN_CONFIG_FILENAME), reverse=True)
        seen_dirs: set[str] = set()
        for run_cfg_path in run_cfg_paths:
            if any(part.startswith(STAGING_PREFIX) for part in run_cfg_path.relative_to(runs_dir).parts):
                continue  # an import being unpacked
            run_cfg = None
            try:
                payload = json.loads(run_cfg_path.read_text(encoding="utf-8"))
            except Exception:
                payload = None
            if isinstance(payload, dict):
                run_cfg = payload
            # The run is where its run_config.json is: the output_dir it recorded
            # names the folder it was made in (another one once moved or shared).
            output_dir = run_cfg_path.parent.resolve()
            if str(output_dir) in seen_dirs:
                continue
            seen_dirs.add(str(output_dir))

            sample_count = None
            if isinstance(run_cfg, dict):
                genomes = run_cfg.get("inputs", {}).get("genomes") if isinstance(run_cfg.get("inputs"), dict) else None
                if isinstance(genomes, list):
                    sample_count = len(genomes)

            status_info = load_run_status(output_dir)
            status = str(status_info.get("status") or "unknown")
            running = bool(status_info.get("active"))
            state = status_info.get("state")
            steps = summarize_run_steps(state if isinstance(state, dict) else None, status)
            pending_text = self._pending_text(output_dir, status_info)
            shared = is_imported_run(output_dir)
            package = is_shared_run(output_dir)
            # A shared run opens whatever the checks found (they are shown first).
            can_view = run_allows_view(status_info) or (package and not running)
            can_resume = run_allows_resume(status_info) and not shared and not _viewer_mode()
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
                steps=steps,
                running=running,
                status=status,
                status_text=status_text,
                can_view=can_view,
                can_resume=can_resume,
                pending_text=pending_text,
                status_detail=status_detail,
                shared=shared,
                issues=tuple(self._run_issues(output_dir, package)),
                metadata_label=_metadata_label(run_id, output_dir),
                version=str((run_cfg or {}).get("bactolomew_version") or ""),
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
                on_more=self._run_action,
            )
            card.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed)
            # The card sits in a transparent holder whose bottom margin is the gap.
            holder = QtWidgets.QWidget()
            holder_layout = QtWidgets.QVBoxLayout(holder)
            holder_layout.setContentsMargins(0, 0, 0, _CARD_GAP)
            holder_layout.setSpacing(0)
            holder_layout.addWidget(card, 0, QtCore.Qt.AlignmentFlag.AlignTop)
            self.run_list.addItem(item)
            self.run_list.setItemWidget(item, holder)
        self._fit_cards()

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
        for item, selected in ((_previous, False), (current, True)):
            card = self._card(item)
            if card is not None:
                card.set_selected(selected)
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
        from .run_share_dialogs import confirm_open_with_warnings

        self._refresh_timer.stop()
        try:
            confirmed = confirm_open_with_warnings(self, entry.sample_set_id or entry.run_id, entry.issues)
        finally:
            self._refresh_timer.start()
        if not confirmed:
            return
        self.output_dir = entry.path
        self.run_config = entry.run_config
        self.run_id = entry.run_id
        self._validate_and_accept()

    def _run_issues(self, output_dir: Path, package: bool) -> list:
        """The run's compatibility findings; a local run's status is on its card already."""
        from .run_share_dialogs import run_issues

        try:
            issues = run_issues(output_dir, self._compat)
        except Exception:
            return []
        return [issue for issue in issues if package or issue.code != "status"]

    # ------------------------------------------------------- More menu

    def _run_action(self, entry: RunEntry, action: str) -> None:
        # After the menu has closed: the refresh rebuilds the cards (and this menu).
        QtCore.QTimer.singleShot(0, lambda: self._run_paused(self._run_action_now, entry, action))

    def _run_paused(self, handler: Callable, *args) -> None:
        """Run ``handler`` with the 5 s refresh stopped (it rebuilds every card)."""
        self._refresh_timer.stop()
        try:
            handler(*args)
        finally:
            if self.isVisible():
                self._refresh_timer.start()

    def _run_action_now(self, entry: RunEntry, action: str) -> None:
        handlers = {
            "details": self._show_run_details,
            "metadata": self._change_metadata,
            "export": self._export_run,
            "folder": self._show_in_folder,
            "remove": self._remove_run,
        }
        handler = handlers.get(action)
        if handler is not None:
            handler(entry)

    def _show_run_details(self, entry: RunEntry) -> None:
        from .run_share_dialogs import RunDetailsDialog

        RunDetailsDialog(entry.path, entry.issues, title=entry.sample_set_id or entry.run_id, parent=self).exec()

    def _change_metadata(self, entry: RunEntry) -> None:
        """The run's metadata CSV or its own data, remembered for the next View."""
        from .app import _choose_metadata

        _choose_metadata(entry.run_id, entry.path, parent=self)
        self._load_runs()

    def _export_run(self, entry: RunEntry) -> None:
        from pipeline import __version__

        from .run_share_dialogs import ExportRunDialog

        ExportRunDialog(
            entry.run_id,
            entry.path,
            title=entry.sample_set_id or entry.run_id,
            sample_count=entry.sample_count,
            app_version=__version__,
            parent=self,
        ).exec()

    def _show_in_folder(self, entry: RunEntry) -> None:
        from PyQt6 import QtGui

        QtGui.QDesktopServices.openUrl(QtCore.QUrl.fromLocalFile(str(entry.path)))

    def _remove_run(self, entry: RunEntry) -> None:
        """Delete an imported run and its cache (only imported runs)."""
        import shutil

        from pipeline.RunShare import ShareError, remove_imported_run

        from .cache import cache_root, run_cache_slug

        answer = QtWidgets.QMessageBox.question(
            self,
            "Remove run",
            f"Delete the imported run {entry.sample_set_id or entry.run_id}?\n{entry.path}",
        )
        if answer != QtWidgets.QMessageBox.StandardButton.Yes:
            return
        try:
            remove_imported_run(entry.path)
        except (OSError, ShareError) as exc:
            QtWidgets.QMessageBox.critical(self, "Remove run", str(exc))
            return
        shutil.rmtree(cache_root() / run_cache_slug(entry.run_id), ignore_errors=True)
        if self.output_dir is not None and self.output_dir == entry.path:
            self.output_dir = None
        self._load_runs()

    def _import_run(self) -> None:
        self._run_paused(self._import_run_now)

    def _import_run_now(self) -> None:
        from .run_share_dialogs import ImportRunDialog

        dlg = ImportRunDialog(runs_root(), parent=self)
        if dlg.exec() != QtWidgets.QDialog.DialogCode.Accepted or dlg.result is None:
            return
        self.output_dir = dlg.result.path.resolve()
        self.run_id = dlg.result.run_id
        self._load_runs()

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
        QtWidgets.QMessageBox.information(
            self, "Stop run",
            "Stop requested. The run will stop after the current step, or during a cgMLST scheme install.",
        )

    def _data_input_from_run_config(self, run_cfg: dict, output_dir: Path) -> Optional["DataInput"]:
        from pipeline.ChewBBACAWorkflow import validate_chewbbaca_settings
        from pipeline.DataInput import DataInput
        from pipeline.PopulationStructure import resume_phylo_overrides
        from pipeline.ToolEnv import (
            configured_assembly_qc_settings,
            configured_bakta_db_path,
            configured_chewbbaca_settings,
            configured_database_dirs,
            configured_fastgear_settings,
            configured_mobsuite_settings,
            relocate_db_path,
        )

        inputs = run_cfg.get("inputs") if isinstance(run_cfg, dict) else None
        if not isinstance(inputs, dict):
            QtWidgets.QMessageBox.critical(self, "Run config error", "run_config.json is missing inputs.")
            return None

        # Re-resolve reference-DB inputs against the current config so a run stays
        # loadable after the project or DB directory was moved/renamed (the stored
        # absolute path is the fallback). Genomes are user inputs outside the DB and
        # are left exactly as recorded.
        # A recorded database removed since is left out (the run records it as not installed).
        db_dirs = configured_database_dirs()
        marker_sequences = [
            path for path in (
                Path(relocate_db_path(p, db_dirs.get("marker_sequences_dir")))
                for p in (inputs.get("marker_sequences") or [])
            ) if path.is_file()
        ]
        reference_locus_sequences = [
            path for path in (
                Path(relocate_db_path(p, db_dirs.get("compound_sequences_dir")))
                for p in (inputs.get("reference_locus_sequences") or [])
            ) if path.is_file()
        ]
        kept = {p.stem for p in (*marker_sequences, *reference_locus_sequences)}
        bakta_db = inputs.get("bakta_db")
        if bakta_db and not Path(str(bakta_db)).expanduser().exists():
            bakta_db = configured_bakta_db_path() or bakta_db

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
                marker_sequences=marker_sequences,
                reference_locus_sequences=reference_locus_sequences,
                output=output_dir,
                snippy_reference=inputs.get("snippy_reference"),
                bakta_db=bakta_db,
                dmnd_threads=inputs.get("dmnd_threads", inputs.get("cores", 1)),
                bakta_threads=inputs.get("bakta_threads", inputs.get("cores", 1)),
                window_length=int(inputs.get("window_length", 500)),
                diamond=inputs.get("diamond", True),
                bakta=inputs.get("bakta", True),
                total_cores=int(inputs.get("cores", 1)),
                marker_categories={k: v for k, v in (inputs.get("marker_categories") or {}).items() if k in kept},
                reference_locus_categories={
                    k: v for k, v in (inputs.get("reference_locus_categories") or {}).items() if k in kept
                },
                fastgear_config=inputs.get("fastgear") or configured_fastgear_settings(),
                mobsuite_config=configured_mobsuite_settings(inputs.get("mobsuite")),
                chewbbaca_config=configured_chewbbaca_settings(inputs.get("chewbbaca")),
                assembly_qc_config=configured_assembly_qc_settings(inputs.get("assembly_qc")),
                # A run keeps the population structure method it used.
                phylo_config=resume_phylo_overrides(output_dir),
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

        ask_cgmlst_terms_before_run(self, entry.path, cfg.chewbbaca_config)

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

        from pipeline.RunShare import is_shared_run

        if is_shared_run(output_dir):
            # A shared run opens with what it has (its warnings were shown).
            if not (output_dir / RUN_CONFIG_FILENAME).is_file():
                errors.append(f"Required file not found: {RUN_CONFIG_FILENAME}")
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

    def _open_settings_dialog(self) -> None:
        from .settings_dialog import SettingsDialog

        # Settings are read from bactolomew_config.json on the next accessor call /
        # next run, so there is nothing to refresh here after the dialog closes.
        SettingsDialog(self).exec()

    def _open_new_run_dialog(self) -> None:
        if _viewer_mode():
            # A viewer install sets up the pipeline first (tools image and databases).
            from .first_time_setup import offer_pipeline_setup

            if not offer_pipeline_setup(self, dev_mode=self._dev_mode):
                return
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
            DEV_BAKTA_DB_PATH,
            configured_assembly_qc_settings,
            configured_chewbbaca_settings,
            configured_fastgear_settings,
            configured_mobsuite_settings,
            effective_core_budget,
        )

        # Only the installed databases are searched; the run logs and records the
        # others as not installed (none skips that search).
        # An unreadable catalog stops the start (as the CLI does).
        marker_catalog = _load_catalog_entries(MARKER_DB_CATALOG_PATH)
        if marker_catalog is None:
            return
        locus_catalog = _load_catalog_entries(COMPOUND_DB_CATALOG_PATH)
        if locus_catalog is None:
            return
        marker_entries = [(p, cat) for p, cat in marker_catalog if p.is_file()]
        locus_entries = [(p, cat) for p, cat in locus_catalog if p.is_file()]
        markers = [p for p, _ in marker_entries]
        loci = [p for p, _ in locus_entries]
        marker_cats = {p.stem: cat for p, cat in marker_entries}
        locus_cats = {p.stem: cat for p, cat in locus_entries}

        bdb = (config.get("bakta_db") or "").strip()
        if self._dev_mode and not bdb:
            bdb = DEV_BAKTA_DB_PATH
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

        cpu = effective_core_budget()  # min(physical cores, configured max_threads)

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
                phylo_config={"population_method": dlg.population_method} if dlg.population_method else {},
                run_id=run_id,
            )
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "Invalid input", str(exc))
            return

        ask_cgmlst_terms_before_run(self, None, cfg.chewbbaca_config)

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
