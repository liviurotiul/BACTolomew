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
from .theme import COLORS


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


def plasmid_gene_rows(ann_df: Optional[pd.DataFrame]) -> pd.DataFrame:
    """The annotations on MOB-suite plasmid contigs, one row each, with the
    columns the plasmid cluster table counts from: ``__sample``,
    ``__cluster_key``, ``__gene_label`` and ``__tags``. Built once per run
    (it reads every annotation) and then only filtered, per cluster."""
    columns = ["__sample", "__cluster_key", "__gene_label", "__tags"]
    if ann_df is None or ann_df.empty or "MobSuiteIsPlasmid" not in ann_df.columns:
        return pd.DataFrame(columns=columns)
    df = ann_df.loc[ann_df["MobSuiteIsPlasmid"].map(_truthy)]
    if df.empty:
        return pd.DataFrame(columns=columns)
    records = df.to_dict(orient="records")
    sample = df["Sample"] if "Sample" in df.columns else pd.Series("", index=df.index)
    tags = df["Tags"] if "Tags" in df.columns else pd.Series("", index=df.index)
    out = pd.DataFrame(
        {
            "__sample": sample.fillna("").astype(str).str.strip().to_numpy(),
            "__cluster_key": [str(plasmid_cluster_key_from_row(record)).strip() for record in records],
            "__gene_label": [_gene_label(record) for record in records],
            "__tags": tags.fillna("").astype(str).to_numpy(),
        }
    )
    return out


class _PlasmidTable(QtWidgets.QTableWidget):
    leftTable = QtCore.pyqtSignal()

    def leaveEvent(self, event: QtCore.QEvent) -> None:
        self.leftTable.emit()
        super().leaveEvent(event)


class PlasmidomePanel(QtWidgets.QWidget):
    samplesHovered = QtCore.pyqtSignal(object)
    sampleSubsetRequested = QtCore.pyqtSignal(object)
    # The row menu's "Open plasmid in new explore window" (plan 3.10).
    sampleSubsetWindowRequested = QtCore.pyqtSignal(object)
    contextRequested = QtCore.pyqtSignal(object)

    _COLUMNS = ["Plasmid", "Samples", "Mobility", "Unique genes", "ARG", "Mobile"]

    def __init__(
        self,
        ann_df: Optional[pd.DataFrame],
        output_dir: Optional[Path],
        plasmid_index: PlasmidClusterIndex,
        parent: Optional[QtWidgets.QWidget] = None,
        *,
        gene_rows: Optional[pd.DataFrame] = None,
    ) -> None:
        super().__init__(parent)
        # Read only (shared with the window); ``gene_rows`` is
        # ``plasmid_gene_rows(ann_df)`` when the caller already has it.
        self._ann_df = ann_df if ann_df is not None else pd.DataFrame()
        self._gene_rows = gene_rows
        self._output_dir = Path(output_dir) if output_dir is not None else None
        self._plasmid_index = plasmid_index
        self._sample_filter: Optional[Set[str]] = None
        # Samples hovered elsewhere in the window: the rows of their plasmids are tinted.
        self._linked_samples: Set[str] = set()
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
        # None shows every cluster; an empty sample set shows none.
        next_filter = (
            {str(sample).strip() for sample in samples if str(sample).strip()} if samples is not None else None
        )
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

    def _plasmid_gene_rows(self) -> pd.DataFrame:
        if self._gene_rows is None:
            self._gene_rows = plasmid_gene_rows(self._ann_df)
        return self._gene_rows

    def _cluster_gene_stats(self, cluster_key: str, samples: Set[str]) -> Dict[str, int]:
        df = self._plasmid_gene_rows()
        if df.empty:
            return {"unique": 0, "arg": 0, "mobile": 0}
        mask = df["__cluster_key"] == cluster_key
        if samples:
            mask &= df["__sample"].isin(samples)
        df = df.loc[mask]
        if df.empty:
            return {"unique": 0, "arg": 0, "mobile": 0}
        unique_genes = {label for label in df["__gene_label"] if label}
        tags = df["__tags"]
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
        if self._linked_samples:
            self._apply_linked_tint()

    def set_linked_samples(self, samples: Optional[Sequence[str]]) -> None:
        """Tint the rows of the plasmids that the samples hovered elsewhere
        carry (none clears)."""
        linked = {str(sample).strip() for sample in samples or [] if str(sample).strip()}
        if linked == self._linked_samples:
            return
        self._linked_samples = linked
        self._apply_linked_tint()

    def linked_rows(self) -> List[int]:
        """The table rows tinted for the hovered samples."""
        return [
            row_idx
            for row_idx, row in self._row_by_table_row.items()
            if self._linked_samples.intersection(row.get("samples") or ())
        ]

    def _apply_linked_tint(self) -> None:
        tinted = set(self.linked_rows())
        tint = QtGui.QBrush(QtGui.QColor(COLORS["primary_soft"]))
        for row_idx in range(self.table.rowCount()):
            brush = tint if row_idx in tinted else QtGui.QBrush()
            for col_idx in range(self.table.columnCount()):
                item = self.table.item(row_idx, col_idx)
                if item is not None:
                    item.setBackground(brush)

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
        menu = self.row_menu(payload)
        menu.exec(self.table.viewport().mapToGlobal(position))
        menu.deleteLater()

    def row_menu(self, payload: Dict[str, Any]) -> QtWidgets.QMenu:
        """The menu of one plasmid row."""
        menu = QtWidgets.QMenu(self.table)
        filter_action = menu.addAction("Filter by this plasmid")
        window_action = menu.addAction("Open plasmid in new explore window")
        window_action.setObjectName("PlasmidWindowAction")
        window_action.setToolTip(
            "A new explore window that shows only the samples carrying this plasmid; this window stays as it is"
        )
        context_action = menu.addAction("Open context")
        filter_action.triggered.connect(lambda _c=False: self._emit_filter_requested(payload))
        window_action.triggered.connect(lambda _c=False: self.sampleSubsetWindowRequested.emit(self._subset_payload(payload)))
        context_action.triggered.connect(
            lambda _c=False: self._emit_context_requested(str(payload.get("cluster_key") or ""))
        )
        return menu

    @staticmethod
    def _subset_payload(payload: Dict[str, Any]) -> Dict[str, Any]:
        return {
            "name": str(payload.get("name") or ""),
            "samples": list(payload.get("samples") or []),
        }

    def _emit_filter_requested(self, payload: Dict[str, Any]) -> None:
        self.sampleSubsetRequested.emit(self._subset_payload(payload))

    def _emit_context_requested(self, cluster_key: str) -> None:
        # Every instance of the plasmid, whatever this view's sample filter:
        # the comparison opens in a neighbourhood window that narrows its rows
        # to that window's own sample set, which the user may widen later.
        rows = build_plasmid_cluster_context_rows(
            self._ann_df,
            self._plasmid_index,
            cluster_key,
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
