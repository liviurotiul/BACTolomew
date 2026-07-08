from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set

import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets

from .panaroo_common import _natural_sort_key
from .plasmid_model import (
    PlasmidClusterIndex,
    build_plasmid_cluster_context_rows,
    plasmid_cluster_key_from_row,
)


_PLASMID_TRUE_VALUES = {"true", "1", "yes", "y", "p", "plasmid"}


def _clean_text(value: Any) -> str:
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    text = str(value or "").strip()
    if text.casefold() in {"", "-", ".", "na", "n/a", "nan", "none", "null", "<na>"}:
        return ""
    return text


def _truthy(value: Any) -> bool:
    return str(value or "").strip().casefold() in _PLASMID_TRUE_VALUES


def _gene_label(row: Dict[str, Any]) -> str:
    for col in ("Gene", "Product", "ConsensusDonorProduct", "GeneID"):
        value = _clean_text(row.get(col))
        if value and not value.upper().startswith("UNK:"):
            return value.casefold()
    return _clean_text(row.get("AnnotationID")) or _clean_text(row.get("Id"))


class _PlasmidTable(QtWidgets.QTableWidget):
    leftTable = QtCore.pyqtSignal()

    def leaveEvent(self, event: QtCore.QEvent) -> None:
        self.leftTable.emit()
        super().leaveEvent(event)


class PlasmidomePanel(QtWidgets.QWidget):
    samplesHovered = QtCore.pyqtSignal(object)
    sampleSubsetRequested = QtCore.pyqtSignal(object)
    contextRequested = QtCore.pyqtSignal(object)

    _COLUMNS = ["Plasmid", "Samples", "Mobility", "Unique genes", "ARG", "Mobile"]

    def __init__(
        self,
        ann_df: Optional[pd.DataFrame],
        output_dir: Optional[Path],
        plasmid_index: PlasmidClusterIndex,
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self._ann_df = ann_df.copy() if ann_df is not None else pd.DataFrame()
        self._output_dir = Path(output_dir) if output_dir is not None else None
        self._plasmid_index = plasmid_index
        self._sample_filter: Optional[Set[str]] = None
        self._rows: List[Dict[str, Any]] = []
        self._row_by_table_row: Dict[int, Dict[str, Any]] = {}
        self._build_ui()
        self.refresh()

    def _build_ui(self) -> None:
        root = QtWidgets.QVBoxLayout(self)
        root.setContentsMargins(10, 10, 10, 10)
        root.setSpacing(8)

        header = QtWidgets.QHBoxLayout()
        title = QtWidgets.QLabel("Plasmidome", self)
        title.setStyleSheet("font-weight: 800; color: #0f172a; font-size: 15px;")
        header.addWidget(title)
        header.addStretch(1)
        self.summary_label = QtWidgets.QLabel("", self)
        self.summary_label.setStyleSheet("color: #64748b; font-weight: 600;")
        header.addWidget(self.summary_label)
        root.addLayout(header)

        self.table = _PlasmidTable(self)
        self.table.setColumnCount(len(self._COLUMNS))
        self.table.setHorizontalHeaderLabels(self._COLUMNS)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setMouseTracking(True)
        self.table.setContextMenuPolicy(QtCore.Qt.ContextMenuPolicy.CustomContextMenu)
        self.table.setAlternatingRowColors(True)
        self.table.setShowGrid(False)
        self.table.verticalHeader().setVisible(False)
        self.table.verticalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Fixed)
        self.table.verticalHeader().setDefaultSectionSize(34)
        self.table.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.Stretch)
        for col in range(1, len(self._COLUMNS)):
            self.table.horizontalHeader().setSectionResizeMode(col, QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.table.cellEntered.connect(self._on_cell_entered)
        self.table.customContextMenuRequested.connect(self._show_row_context_menu)
        self.table.leftTable.connect(lambda: self.samplesHovered.emit(None))
        root.addWidget(self.table, stretch=1)

    def set_sample_filter(self, samples: Optional[Sequence[str]]) -> None:
        next_filter = {str(sample).strip() for sample in samples or [] if str(sample).strip()} or None
        if next_filter == self._sample_filter:
            return
        self._sample_filter = next_filter
        self.refresh()

    def _cluster_samples(self, cluster_key: str) -> List[str]:
        summary = self._plasmid_index.clusters.get(cluster_key)
        samples = set(summary.samples) if summary is not None else set()
        if self._sample_filter is not None:
            samples &= self._sample_filter
        return sorted(samples, key=_natural_sort_key)

    def _cluster_gene_stats(self, cluster_key: str, samples: Set[str]) -> Dict[str, int]:
        if self._ann_df.empty or "MobSuiteIsPlasmid" not in self._ann_df.columns:
            return {"unique": 0, "arg": 0, "mobile": 0}
        df = self._ann_df.loc[self._ann_df["MobSuiteIsPlasmid"].map(_truthy)].copy()
        if df.empty:
            return {"unique": 0, "arg": 0, "mobile": 0}
        if samples:
            df = df.loc[df["Sample"].fillna("").astype(str).str.strip().isin(samples)].copy()
        if df.empty:
            return {"unique": 0, "arg": 0, "mobile": 0}
        records = df.to_dict(orient="records")
        keys = [plasmid_cluster_key_from_row(record) for record in records]
        df["__cluster_key"] = keys
        df = df.loc[df["__cluster_key"].astype(str).str.strip() == cluster_key]
        if df.empty:
            return {"unique": 0, "arg": 0, "mobile": 0}
        unique_genes = {
            label
            for label in (_gene_label(row) for row in df.to_dict(orient="records"))
            if label
        }
        tags = df.get("Tags", pd.Series("", index=df.index)).fillna("").astype(str)
        return {
            "unique": len(unique_genes),
            "arg": int(tags.str.contains("A", regex=False).sum()),
            "mobile": int(tags.str.contains("M", regex=False).sum()),
        }

    def _mobility_text(self, cluster_key: str) -> str:
        summary = self._plasmid_index.clusters.get(cluster_key)
        if summary is None:
            return ""
        mobility = sorted(summary.mobility, key=_natural_sort_key)
        if mobility:
            return ", ".join(mobility[:3])
        if summary.mpf_types:
            return "conjugative?"
        if summary.relaxase_types:
            return "mobilizable?"
        return ""

    def _build_rows(self) -> List[Dict[str, Any]]:
        rows: List[Dict[str, Any]] = []
        for cluster_key in self._plasmid_index.cluster_order(self._sample_filter):
            samples = self._cluster_samples(cluster_key)
            if not samples:
                continue
            stats = self._cluster_gene_stats(cluster_key, set(samples))
            summary = self._plasmid_index.clusters.get(cluster_key)
            rows.append(
                {
                    "cluster_key": cluster_key,
                    "name": self._plasmid_index.display_name_for_cluster(cluster_key),
                    "samples": samples,
                    "sample_count": len(samples),
                    "mobility": self._mobility_text(cluster_key),
                    "unique_genes": stats["unique"],
                    "arg_genes": stats["arg"],
                    "mobile_genes": stats["mobile"],
                    "tooltip": self._plasmid_index.tooltip_for_cluster(cluster_key),
                    "prevalence": summary.prevalence if summary is not None else len(samples),
                }
            )
        rows.sort(key=lambda row: (-int(row.get("sample_count") or 0), _natural_sort_key(str(row.get("name") or ""))))
        return rows

    def refresh(self) -> None:
        self._rows = self._build_rows()
        self._render_table()
        self.summary_label.setText(f"{len(self._rows)} plasmid cluster(s)")

    def _render_table(self) -> None:
        self._row_by_table_row = {}
        self.table.setRowCount(len(self._rows))
        for row_idx, row in enumerate(self._rows):
            self._row_by_table_row[row_idx] = row
            values = [
                str(row["name"]),
                str(row["sample_count"]),
                str(row["mobility"] or ""),
                str(row["unique_genes"]),
                str(row["arg_genes"]),
                str(row["mobile_genes"]),
            ]
            tooltip = str(row.get("tooltip") or row["name"])
            for col_idx, value in enumerate(values):
                item = QtWidgets.QTableWidgetItem(value)
                item.setToolTip(tooltip)
                if col_idx == 0:
                    font = item.font()
                    font.setBold(True)
                    item.setFont(font)
                    item.setForeground(QtGui.QBrush(QtGui.QColor("#14532d")))
                else:
                    item.setTextAlignment(int(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter))
                self.table.setItem(row_idx, col_idx, item)
            self.table.setRowHeight(row_idx, 34)

    def _on_cell_entered(self, row: int, _column: int) -> None:
        payload = self._row_by_table_row.get(int(row))
        if not payload:
            self.samplesHovered.emit(None)
            return
        self.samplesHovered.emit(list(payload.get("samples") or []))

    def _show_row_context_menu(self, position: QtCore.QPoint) -> None:
        index = self.table.indexAt(position)
        if not index.isValid():
            return
        row = int(index.row())
        payload = self._row_by_table_row.get(int(row))
        if not payload:
            return
        self.table.selectRow(row)
        menu = QtWidgets.QMenu(self.table)
        filter_action = menu.addAction("Filter by this plasmid")
        context_action = menu.addAction("Open context")
        chosen = menu.exec(self.table.viewport().mapToGlobal(position))
        if chosen == filter_action:
            self._emit_filter_requested(payload)
        elif chosen == context_action:
            self._emit_context_requested(str(payload.get("cluster_key") or ""))

    def _emit_filter_requested(self, payload: Dict[str, Any]) -> None:
        self.sampleSubsetRequested.emit(
            {
                "name": str(payload.get("name") or ""),
                "samples": list(payload.get("samples") or []),
            }
        )

    def _emit_context_requested(self, cluster_key: str) -> None:
        rows = build_plasmid_cluster_context_rows(
            self._ann_df,
            self._plasmid_index,
            cluster_key,
            sample_filter=self._sample_filter,
        )
        name = self._plasmid_index.display_name_for_cluster(cluster_key)
        self.contextRequested.emit(
            {
                "name": name,
                "cluster_key": cluster_key,
                "samples": self._cluster_samples(cluster_key),
                "rows": rows,
            }
        )
