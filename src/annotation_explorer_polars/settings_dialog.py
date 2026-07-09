"""Settings dialog reachable from the runs window's gear button.

Phase 1 exposes the global compute limits (max CPU threads, max container RAM)
persisted to bactolomew_config.json's `resources` block via ToolEnv.save_app_config,
plus a read-only view of the installed cgMLST schemas. The "Add schema…" import
flow is wired in a later phase.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

from PyQt6 import QtCore, QtGui, QtWidgets

# Local copies of the app's control styles (duplicated to avoid importing
# startup_dialog, which imports this module).
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

_GROUP_STYLE = (
    "QGroupBox{font-weight:700;color:#1f2c3b;border:1px solid #d7dfe8;border-radius:10px;"
    "margin-top:10px;padding:12px 12px 12px 12px;background:#ffffff;}"
    "QGroupBox::title{subcontrol-origin:margin;left:12px;padding:0 4px;}"
)
_HINT_STYLE = "color:#64748b;font-size:11px;"


class SettingsDialog(QtWidgets.QDialog):
    """Edit global settings. Reads/writes bactolomew_config.json (merge-based)."""

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Settings")
        self.setObjectName("SettingsDialog")
        self.setStyleSheet("QDialog#SettingsDialog { background: #eef3f8; }")
        self.resize(560, 500)
        self._cpu_total = os.cpu_count() or 1
        self._build_ui()
        self._load_current()

    # ---- UI ----
    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(22, 20, 22, 18)
        layout.setSpacing(12)

        title = QtWidgets.QLabel("Settings")
        title.setStyleSheet("font-size: 18px; font-weight: 800; color: #172433;")
        layout.addWidget(title)

        layout.addWidget(self._build_resources_group())
        layout.addWidget(self._build_schemas_group(), stretch=1)

        btn_row = QtWidgets.QHBoxLayout()
        btn_row.addStretch(1)
        cancel_btn = QtWidgets.QPushButton("Cancel")
        cancel_btn.setStyleSheet(_GHOST_BTN)
        cancel_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        cancel_btn.clicked.connect(self.reject)
        btn_row.addWidget(cancel_btn)
        save_btn = QtWidgets.QPushButton("Save")
        save_btn.setStyleSheet(_PRIMARY_BTN)
        save_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        save_btn.clicked.connect(self._save)
        btn_row.addWidget(save_btn)
        layout.addLayout(btn_row)

    def _build_resources_group(self) -> QtWidgets.QGroupBox:
        group = QtWidgets.QGroupBox("Resources")
        group.setStyleSheet(_GROUP_STYLE)
        form = QtWidgets.QFormLayout(group)
        form.setContentsMargins(12, 16, 12, 12)
        form.setSpacing(10)

        # Max CPU threads
        self.threads_spin = QtWidgets.QSpinBox()
        self.threads_spin.setRange(8, max(8, self._cpu_total))
        self.threads_spin.setSuffix(" threads")
        form.addRow("Max CPU threads", self.threads_spin)

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
        form.addRow("Max RAM", ram_row)

        hint = QtWidgets.QLabel(
            f"Physical cores available: {self._cpu_total}. The thread cap applies to every new "
            "run (a per-run --cores still overrides). The RAM limit is applied to Docker tool "
            "containers (--memory) and enables CheckM2 low-memory mode when tight."
        )
        hint.setWordWrap(True)
        hint.setStyleSheet(_HINT_STYLE)
        form.addRow(hint)
        return group

    def _build_schemas_group(self) -> QtWidgets.QGroupBox:
        group = QtWidgets.QGroupBox("cgMLST schemas")
        group.setStyleSheet(_GROUP_STYLE)
        v = QtWidgets.QVBoxLayout(group)
        v.setContentsMargins(12, 16, 12, 12)
        v.setSpacing(8)

        self.schema_list = QtWidgets.QListWidget()
        self.schema_list.setStyleSheet(
            "QListWidget{background:#f7f9fc;border:1px solid #e2e8f0;border-radius:8px;}"
        )
        v.addWidget(self.schema_list, stretch=1)

        row = QtWidgets.QHBoxLayout()
        add_btn = QtWidgets.QPushButton("Add schema…")
        add_btn.setStyleSheet(_SECONDARY_BTN)
        add_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        add_btn.clicked.connect(self._add_schema)
        row.addWidget(add_btn)
        row.addStretch(1)
        v.addLayout(row)

        hint = QtWidgets.QLabel(
            "Installed species schemas under DB/cgmlst. A run auto-selects the matching schema "
            "from its MLST result."
        )
        hint.setWordWrap(True)
        hint.setStyleSheet(_HINT_STYLE)
        v.addWidget(hint)
        return group

    # ---- data ----
    def _load_current(self) -> None:
        from pipeline.ToolEnv import configured_resources

        res = configured_resources()
        self.threads_spin.setValue(max(8, min(res["max_threads"], max(8, self._cpu_total))))
        ram = res["max_ram_gb"]
        self.ram_check.setChecked(ram > 0)
        self.ram_spin.setEnabled(ram > 0)
        self.ram_spin.setValue(ram if ram >= 16 else 16)

        self._reload_schema_list()

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

    # ---- actions ----
    def _add_schema(self) -> None:
        QtWidgets.QMessageBox.information(
            self,
            "Add schema",
            "Importing a new cgMLST schema (single FASTA, a folder of allele FASTAs, or a "
            "prebuilt chewBBACA schema) is being added in the next update.",
        )

    def _save(self) -> None:
        from pipeline.ToolEnv import save_app_config

        resources = {
            "max_threads": int(self.threads_spin.value()),
            "max_ram_gb": int(self.ram_spin.value()) if self.ram_check.isChecked() else 0,
        }
        try:
            save_app_config({"resources": resources})
        except Exception as exc:  # pragma: no cover - surfaced to the user
            QtWidgets.QMessageBox.critical(self, "Save failed", f"Could not save settings:\n{exc}")
            return
        self.accept()
