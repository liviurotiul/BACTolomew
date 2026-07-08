from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path
from typing import Dict, Iterable, Optional

from PyQt6 import QtCore, QtGui, QtWidgets


ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from annotation_explorer_polars.panaroo_force import (  # noqa: E402
    DEFAULT_CGMLST_FORCE_CONFIG,
    CgmlstForceConfig,
    CgmlstForceGraphView,
)
from annotation_explorer_polars.panaroo_data import (  # noqa: E402
    _find_fastgear_sample_clades,
    _load_fastgear_clusters,
    _load_fastgear_lineages,
)
from annotation_explorer_polars.panaroo_style import _fastgear_lineage_colors  # noqa: E402


PREFERRED_RUN_DIR = ROOT / "runs" / "sample_78eaa9126850" / "run_78eaa9126850"


def _discover_default_run_dir() -> Path:
    if PREFERRED_RUN_DIR.is_dir():
        return PREFERRED_RUN_DIR
    candidates = sorted(
        (
            path.parent.parent.parent
            for path in ROOT.glob("runs/*/*/chewbbaca/out/results_alleles.tsv")
            if path.is_file()
        ),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    if candidates:
        return candidates[0]
    return PREFERRED_RUN_DIR


def _load_reference_sample(run_dir: Path, explicit: Optional[str]) -> Optional[str]:
    if explicit:
        return str(explicit).strip() or None
    run_config_path = run_dir / "run_config.json"
    if not run_config_path.is_file():
        return None
    try:
        config = json.loads(run_config_path.read_text(encoding="utf-8"))
    except Exception:
        return None
    ref = str(config.get("inputs", {}).get("snippy_reference", "") or "").strip()
    if not ref:
        return None
    return Path(ref).stem or ref


def _apply_fastgear_colors(view: CgmlstForceGraphView, run_dir: Path) -> None:
    sample_clades = _find_fastgear_sample_clades(run_dir)
    if sample_clades is None:
        return
    lineage_by_sample = _load_fastgear_lineages(sample_clades)
    cluster_by_sample = _load_fastgear_clusters(sample_clades)
    lineage_colors = _fastgear_lineage_colors(lineage_by_sample)
    view.set_fastgear_assignments(lineage_by_sample, lineage_colors, cluster_by_sample)


class _NumericTuningRow(QtWidgets.QWidget):
    valueChanged = QtCore.pyqtSignal()

    def __init__(
        self,
        name: str,
        *,
        minimum: float,
        maximum: float,
        step: float,
        initial: float | int,
        integer: bool = False,
        help_text: str = "",
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self.name = str(name)
        self._minimum = float(minimum)
        self._maximum = float(maximum)
        self._step = float(step)
        self._integer = bool(integer)
        self._slider_min = 0
        self._slider_max = int(round((self._maximum - self._minimum) / self._step))

        layout = QtWidgets.QGridLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setHorizontalSpacing(10)
        layout.setVerticalSpacing(3)

        label = QtWidgets.QLabel(self.name.replace("_", " "))
        font = label.font()
        font.setBold(True)
        label.setFont(font)
        layout.addWidget(label, 0, 0)

        self._slider = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self._slider.setRange(self._slider_min, self._slider_max)
        layout.addWidget(self._slider, 0, 1)

        if self._integer:
            spin: QtWidgets.QAbstractSpinBox = QtWidgets.QSpinBox()
            cast_spin = spin
            cast_spin.setRange(int(round(self._minimum)), int(round(self._maximum)))
            cast_spin.setSingleStep(int(max(1, round(self._step))))
            cast_spin.setValue(int(round(float(initial))))
        else:
            spin = QtWidgets.QDoubleSpinBox()
            cast_spin = spin
            step_text = f"{self._step:.6f}".rstrip("0").rstrip(".")
            decimals = len(step_text.split(".", 1)[1]) if "." in step_text else 0
            cast_spin.setDecimals(max(3, decimals))
            cast_spin.setRange(self._minimum, self._maximum)
            cast_spin.setSingleStep(self._step)
            cast_spin.setValue(float(initial))
        self._spin = cast_spin
        self._spin.setKeyboardTracking(False)
        self._spin.setMinimumWidth(110)
        layout.addWidget(self._spin, 0, 2)

        if help_text:
            help_label = QtWidgets.QLabel(help_text)
            help_label.setWordWrap(True)
            help_label.setStyleSheet("color: #475569; font-size: 11px;")
            layout.addWidget(help_label, 1, 0, 1, 3)

        self._slider.valueChanged.connect(self._on_slider_changed)
        self._spin.valueChanged.connect(self._on_spin_changed)
        self.set_value(initial)

    def _slider_to_value(self, slider_value: int) -> float | int:
        value = self._minimum + (float(slider_value) * self._step)
        if self._integer:
            return int(round(value))
        return round(value, 6)

    def _value_to_slider(self, value: float | int) -> int:
        return int(round((float(value) - self._minimum) / self._step))

    def _on_slider_changed(self, slider_value: int) -> None:
        value = self._slider_to_value(slider_value)
        blocker = QtCore.QSignalBlocker(self._spin)
        try:
            if self._integer:
                self._spin.setValue(int(value))
            else:
                self._spin.setValue(float(value))
        finally:
            del blocker
        self.valueChanged.emit()

    def _on_spin_changed(self, value: object) -> None:
        slider_value = self._value_to_slider(float(value))
        blocker = QtCore.QSignalBlocker(self._slider)
        try:
            self._slider.setValue(slider_value)
        finally:
            del blocker
        self.valueChanged.emit()

    def set_value(self, value: float | int) -> None:
        value_f = min(self._maximum, max(self._minimum, float(value)))
        slider_value = self._value_to_slider(value_f)
        blocker_slider = QtCore.QSignalBlocker(self._slider)
        blocker_spin = QtCore.QSignalBlocker(self._spin)
        try:
            self._slider.setValue(slider_value)
            if self._integer:
                self._spin.setValue(int(round(value_f)))
            else:
                self._spin.setValue(float(value_f))
        finally:
            del blocker_slider
            del blocker_spin

    def value(self) -> float | int:
        value = self._spin.value()
        return int(value) if self._integer else float(value)


class ForceGraphTunerWindow(QtWidgets.QMainWindow):
    def __init__(self, run_dir: Path, reference_sample: Optional[str]) -> None:
        super().__init__()
        self._run_dir = Path(run_dir)
        self._reference_sample = str(reference_sample or "").strip() or None
        self._apply_timer = QtCore.QTimer(self)
        self._apply_timer.setInterval(140)
        self._apply_timer.setSingleShot(True)
        self._apply_timer.timeout.connect(self._apply_config)
        self._controls: Dict[str, _NumericTuningRow] = {}

        self.setWindowTitle("BACTolomew cgMLST Force Graph Tuner")
        self.resize(1680, 1020)

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        root = QtWidgets.QHBoxLayout(central)
        root.setContentsMargins(10, 10, 10, 10)
        root.setSpacing(12)

        self._view = CgmlstForceGraphView(
            self._run_dir,
            reference_sample=self._reference_sample,
        )
        _apply_fastgear_colors(self._view, self._run_dir)
        root.addWidget(self._view, stretch=3)

        right_panel = QtWidgets.QWidget()
        right_layout = QtWidgets.QVBoxLayout(right_panel)
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.setSpacing(10)
        root.addWidget(right_panel, stretch=2)

        summary = QtWidgets.QLabel(
            f"Run: {self._run_dir}\n"
            f"Reference: {self._reference_sample or 'none'}\n"
            "Move sliders to tune the live cgMLST force layout. The graph restarts after short input pauses."
        )
        summary.setWordWrap(True)
        summary.setStyleSheet("font-size: 12px;")
        right_layout.addWidget(summary)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        right_layout.addWidget(scroll, stretch=1)

        scroll_body = QtWidgets.QWidget()
        self._form_layout = QtWidgets.QVBoxLayout(scroll_body)
        self._form_layout.setContentsMargins(0, 0, 0, 0)
        self._form_layout.setSpacing(12)
        scroll.setWidget(scroll_body)

        self._build_control_groups()

        buttons = QtWidgets.QHBoxLayout()
        right_layout.addLayout(buttons)
        reset_btn = QtWidgets.QPushButton("Reset Defaults")
        reset_btn.clicked.connect(self._reset_defaults)
        buttons.addWidget(reset_btn)
        restart_btn = QtWidgets.QPushButton("Restart Simulation")
        restart_btn.clicked.connect(self._view.restart_simulation)
        buttons.addWidget(restart_btn)
        copy_btn = QtWidgets.QPushButton("Copy Params")
        copy_btn.clicked.connect(self._copy_params)
        buttons.addWidget(copy_btn)
        buttons.addStretch(1)

        self._params_box = QtWidgets.QPlainTextEdit()
        self._params_box.setReadOnly(True)
        self._params_box.setMinimumHeight(240)
        right_layout.addWidget(self._params_box)

        self._update_params_box()
        QtCore.QTimer.singleShot(0, self._view.refresh)

    def _build_control_groups(self) -> None:
        self._add_group(
            "Main forces",
            [
                ("charge_strength", 0.0005, 0.0500, 0.0005, False, "Node-node repulsion. Higher values push samples apart."),
                ("collision_radius", 0.0200, 0.4000, 0.0050, False, "Minimum personal space for nodes before hard collision push kicks in."),
                ("collision_strength", 0.0000, 5.0000, 0.0500, False, "Strength of the close-range collision separation."),
                ("spring_strength", 0.0000, 2.0000, 0.0200, False, "How hard MST edges pull toward their target lengths."),
                ("rest_length_min", 0.0000, 1.5000, 0.0100, False, "Base edge length for zero or minimal allelic distance."),
                ("rest_length_span", 0.0000, 3.5000, 0.0200, False, "Extra edge length added as allelic distance increases."),
                ("anchor_strength", 0.0000, 0.3000, 0.0050, False, "Pull back toward the initial layout. Lower values allow freer motion."),
                ("gravity_strength", 0.0000, 0.1000, 0.0020, False, "Soft pull toward the graph center."),
                ("reference_anchor_strength", 0.0000, 0.6000, 0.0050, False, "How strongly the root/reference sample is held near its initial position."),
                ("node_edge_radius", 0.0000, 0.5000, 0.0050, False, "How close a node can approach a non-incident edge before being pushed away."),
                ("node_edge_strength", 0.0000, 0.5000, 0.0050, False, "Strength of node-edge separation."),
            ],
        )
        self._add_group(
            "Motion",
            [
                ("dt_base", 0.0010, 0.0600, 0.0010, False, "Base integration step."),
                ("dt_alpha_scale", 0.0000, 0.0800, 0.0010, False, "Extra integration step while the simulation is still hot."),
                ("damping", 0.5000, 0.9900, 0.0050, False, "Velocity retention per tick. Higher is livelier."),
                ("max_speed_base", 0.0010, 0.1000, 0.0010, False, "Hard speed cap once the layout settles."),
                ("max_speed_alpha_scale", 0.0000, 0.1200, 0.0010, False, "Extra speed allowance while the layout is still hot."),
                ("cooling", 0.9000, 0.9990, 0.0010, False, "How quickly the simulation cools per tick."),
                ("stop_alpha", 0.0010, 0.2000, 0.0010, False, "Stop when alpha falls below this."),
                ("stop_speed", 0.00005, 0.01000, 0.00005, False, "Stop when mean node speed falls below this."),
                ("max_ticks", 50, 2000, 10, True, "Hard stop limit for simulation ticks."),
            ],
        )
        self._add_group(
            "Node sizes",
            [
                ("node_size", 4, 48, 1, True, "Default node size."),
                ("focus_node_size", 4, 56, 1, True, "Size for highlighted nodes."),
                ("hover_node_size", 4, 64, 1, True, "Size for the hovered node."),
                ("reference_node_size", 4, 56, 1, True, "Minimum size for the reference node."),
            ],
        )

    def _add_group(
        self,
        title: str,
        specs: Iterable[tuple[str, float, float, float, bool, str]],
    ) -> None:
        box = QtWidgets.QGroupBox(title)
        layout = QtWidgets.QVBoxLayout(box)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(10)
        defaults = DEFAULT_CGMLST_FORCE_CONFIG.to_dict()
        for name, minimum, maximum, step, integer, help_text in specs:
            row = _NumericTuningRow(
                name,
                minimum=minimum,
                maximum=maximum,
                step=step,
                initial=defaults[name],
                integer=integer,
                help_text=help_text,
            )
            row.valueChanged.connect(self._schedule_apply)
            self._controls[name] = row
            layout.addWidget(row)
        self._form_layout.addWidget(box)

    def _schedule_apply(self) -> None:
        self._update_params_box()
        self._apply_timer.start()

    def _current_config(self) -> CgmlstForceConfig:
        cfg = DEFAULT_CGMLST_FORCE_CONFIG
        for name, control in self._controls.items():
            cfg = replace(cfg, **{name: control.value()})
        return cfg

    def _apply_config(self) -> None:
        self._view.set_force_config(self._current_config(), restart=True)
        self._update_params_box()

    def _reset_defaults(self) -> None:
        defaults = DEFAULT_CGMLST_FORCE_CONFIG.to_dict()
        for name, control in self._controls.items():
            control.set_value(defaults[name])
        self._apply_config()

    def _copy_params(self) -> None:
        payload = json.dumps(self._current_config().to_dict(), indent=2, sort_keys=True)
        QtWidgets.QApplication.clipboard().setText(payload)

    def _update_params_box(self) -> None:
        config = self._current_config().to_dict()
        self._params_box.setPlainText(json.dumps(config, indent=2, sort_keys=True))


def main() -> None:
    parser = argparse.ArgumentParser(description="Tune BACTolomew cgMLST force-graph parameters.")
    parser.add_argument(
        "--run-dir",
        default=str(_discover_default_run_dir()),
        help="Run directory containing chewBBACA outputs. Defaults to the current ST492 run if present.",
    )
    parser.add_argument(
        "--reference-sample",
        default="",
        help="Override the reference sample label used to pin the graph.",
    )
    args = parser.parse_args()

    run_dir = Path(args.run_dir).expanduser().resolve()
    if not run_dir.is_dir():
        raise SystemExit(f"Run directory not found: {run_dir}")

    app = QtWidgets.QApplication(sys.argv)
    app.setApplicationName("BACTolomew cgMLST Force Tuner")
    reference_sample = _load_reference_sample(run_dir, args.reference_sample)
    window = ForceGraphTunerWindow(run_dir, reference_sample)
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
