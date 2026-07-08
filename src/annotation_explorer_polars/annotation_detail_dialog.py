from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
import traceback
from collections import Counter
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets

from .constants import PLASMID_TRUE_VALUES
from .igv_collection_view import IgvCollectionDialog

_HIDDEN_INTERNAL_MARKER_FIELDS = {
    "SubjectIdentifier",
    "SubjectProteinAccession",
    "SubjectDNAAccession",
    "SubjectAROAccession",
    "SubjectModelName",
    "ResistanceMetadataSource",
    "ResistanceAlleles",
    "ResistanceGeneFamilies",
    "ResistanceDrugClasses",
    "ResistanceAntibiotics",
    "ResistanceMechanisms",
}


def _normalize_detail_value(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, np.ndarray):
        return [_normalize_detail_value(v) for v in value.tolist()]
    if isinstance(value, dict):
        return {str(k): _normalize_detail_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_normalize_detail_value(v) for v in value]
    if isinstance(value, (pd.Series, pd.Index)):
        return [_normalize_detail_value(v) for v in value.tolist()]
    if hasattr(value, "as_py"):
        try:
            return _normalize_detail_value(value.as_py())
        except Exception:
            pass
    if hasattr(value, "to_pylist"):
        try:
            return _normalize_detail_value(value.to_pylist())
        except Exception:
            pass
    try:
        if pd.isna(value):
            return None
    except Exception:
        pass
    return value


def _safe_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and pd.isna(value):
        return ""
    s = str(value).strip()
    if not s or s.lower() == "nan":
        return ""
    return s


def _parse_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, float) and pd.isna(value):
        return None
    try:
        return float(value)
    except Exception:
        return None


def _clean_sequence_text(value: Any) -> str:
    seq = _safe_value(value)
    if not seq:
        return ""
    return "".join(ch for ch in seq if not ch.isspace())


def _safe_header_part(value: Any) -> str:
    raw = _safe_value(value)
    if not raw:
        return ""
    return "".join(ch if (ch.isalnum() or ch in "-_.|") else "_" for ch in raw)


def _parse_pasted_sequence(text: str) -> str:
    raw = (text or "").strip()
    if not raw:
        return ""
    lines = [ln.strip() for ln in raw.splitlines() if ln.strip()]
    if not lines:
        return ""
    if lines[0].startswith(">"):
        seq_lines = [ln for ln in lines[1:] if not ln.startswith(">")]
        return _clean_sequence_text("".join(seq_lines))
    return _clean_sequence_text("".join(lines))


def _show_floating_dialog(
    owner_store: List[QtWidgets.QDialog],
    dlg: QtWidgets.QDialog,
) -> None:
    dlg.setWindowModality(QtCore.Qt.WindowModality.NonModal)
    dlg.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose, True)
    owner_store.append(dlg)

    def _cleanup(*_args: object, _dlg: QtWidgets.QDialog = dlg) -> None:
        if _dlg in owner_store:
            owner_store.remove(_dlg)

    dlg.destroyed.connect(_cleanup)
    dlg.show()
    dlg.raise_()
    dlg.activateWindow()


class _ContextCollectionSignals(QtCore.QObject):
    finished = QtCore.pyqtSignal(int, object)
    failed = QtCore.pyqtSignal(int, str)


class _ContextCollectionWorker(QtCore.QRunnable):
    def __init__(self, token: int, compute: Callable[[], List[Dict[str, Any]]]) -> None:
        super().__init__()
        self._token = int(token)
        self._compute = compute
        self.signals = _ContextCollectionSignals()
        self.setAutoDelete(True)

    def run(self) -> None:
        try:
            rows = self._compute()
        except Exception:
            self.signals.failed.emit(self._token, traceback.format_exc())
            return
        self.signals.finished.emit(self._token, rows)


class ZoomableAlignmentTable(QtWidgets.QTableWidget):
    zoomRequested = QtCore.pyqtSignal(int)

    def __init__(self, rows: int, cols: int, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(rows, cols, parent)
        # On many Qt styles, wheel events are delivered to the viewport first.
        self.viewport().installEventFilter(self)

    def _handle_wheel_zoom(self, event: QtGui.QWheelEvent) -> bool:
        delta = event.angleDelta().y()
        if delta == 0:
            delta = event.pixelDelta().y()
        if delta == 0:
            return False
        step = 2 if delta > 0 else -2
        self.zoomRequested.emit(step)
        event.accept()
        return True

    def eventFilter(self, obj: QtCore.QObject, event: QtCore.QEvent) -> bool:
        if obj is self.viewport() and event.type() == QtCore.QEvent.Type.Wheel:
            wheel_event = event  # type: ignore[assignment]
            if self._handle_wheel_zoom(wheel_event):
                return True
        return super().eventFilter(obj, event)

    def wheelEvent(self, event: QtGui.QWheelEvent) -> None:  # type: ignore[override]
        if self._handle_wheel_zoom(event):
            return
        super().wheelEvent(event)


class AlignmentResultDialog(QtWidgets.QDialog):
    def __init__(
        self,
        parent: Optional[QtWidgets.QWidget],
        *,
        alignment_text: str,
        status_text: str,
        sample_hover_callback: Optional[Callable[[Optional[str]], None]] = None,
    ) -> None:
        super().__init__(parent)
        self._alignment_text = alignment_text
        self._sample_hover_callback = sample_hover_callback
        self._syncing_scroll = False
        self._col_size = 22
        self._row_height = 22
        self._info_table: Optional[QtWidgets.QTableWidget] = None
        self._seq_table: Optional[ZoomableAlignmentTable] = None
        self._row_count = 0
        self._seq_col_count = 0
        self._row_samples: List[Optional[str]] = []
        self._hovered_sample: Optional[str] = None
        self._zoom_event_sources: set[QtCore.QObject] = set()
        self._zoom_pending_step = 0
        self._zoom_timer = QtCore.QTimer(self)
        self._zoom_timer.setSingleShot(True)
        self._zoom_timer.setInterval(40)
        self._zoom_timer.timeout.connect(self._flush_zoom_steps)
        self._row_stretch_timer = QtCore.QTimer(self)
        self._row_stretch_timer.setSingleShot(True)
        self._row_stretch_timer.setInterval(60)
        self._row_stretch_timer.timeout.connect(self._stretch_rows_to_fill_viewport)
        self._last_applied_col_size = self._col_size
        self.setWindowTitle("Alignment Result")
        self.resize(1200, 560)
        self._build_ui(status_text)

    def _build_ui(self, status_text: str) -> None:
        _ = status_text
        layout = QtWidgets.QVBoxLayout(self)
        records = self._parse_fasta_alignment(self._alignment_text)

        if records:
            left_table, right_table = self._build_alignment_tables(records)
            self._info_table = left_table
            self._seq_table = right_table
            self._register_zoom_sources(left_table, right_table)
            splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
            splitter.setChildrenCollapsible(False)
            splitter.setHandleWidth(1)
            splitter.addWidget(left_table)
            splitter.addWidget(right_table)
            splitter.setStretchFactor(0, 0)
            splitter.setStretchFactor(1, 1)
            splitter.setSizes([380, 940])
            layout.addWidget(splitter, stretch=1)
            self._fit_tables_height(left_table, right_table)
            self._apply_zoom()
        else:
            raw = QtWidgets.QPlainTextEdit()
            raw.setReadOnly(True)
            raw.setLineWrapMode(QtWidgets.QPlainTextEdit.LineWrapMode.NoWrap)
            raw_font = QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.SystemFont.FixedFont)
            raw.setFont(raw_font)
            raw.setPlainText(self._alignment_text)
            layout.addWidget(raw, stretch=1)

        btn_row = QtWidgets.QHBoxLayout()
        btn_row.addStretch(1)
        save_btn = QtWidgets.QPushButton("Save FASTA")
        save_btn.clicked.connect(self._save_alignment)
        close_btn = QtWidgets.QPushButton("Close")
        close_btn.clicked.connect(self.accept)
        btn_row.addWidget(save_btn)
        btn_row.addWidget(close_btn)
        layout.addLayout(btn_row)
        self.adjustSize()

    def _parse_fasta_alignment(self, text: str) -> List[Tuple[str, str]]:
        records: List[Tuple[str, str]] = []
        header: Optional[str] = None
        seq_parts: List[str] = []
        for raw_line in (text or "").splitlines():
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = _clean_sequence_text("".join(seq_parts))
                    if seq:
                        records.append((header, seq))
                header = line[1:].strip() or f"seq_{len(records) + 1}"
                seq_parts = []
            else:
                if header is None:
                    # Not FASTA; abort colored rendering.
                    return []
                seq_parts.append(line)
        if header is not None:
            seq = _clean_sequence_text("".join(seq_parts))
            if seq:
                records.append((header, seq))
        return records

    def _build_alignment_tables(
        self, records: List[Tuple[str, str]]
    ) -> Tuple[ZoomableAlignmentTable, ZoomableAlignmentTable]:
        max_len = max(len(seq) for _, seq in records)
        padded = [(hdr, seq.ljust(max_len, "-")) for hdr, seq in records]
        consensus_marks = self._consensus_marks([seq for _, seq in padded])
        row_labels = [hdr for hdr, _ in padded] + ["consensus"]
        self._row_samples = [self._sample_from_alignment_header(lbl) for lbl in row_labels]
        row_count = len(row_labels)
        self._row_count = row_count
        self._seq_col_count = max_len

        left = ZoomableAlignmentTable(row_count, 1)
        left.setHorizontalHeaderLabels(["Sequence"])
        left.verticalHeader().setVisible(False)
        left.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        left.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.NoSelection)
        left.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        left.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        left.setVerticalScrollMode(QtWidgets.QAbstractItemView.ScrollMode.ScrollPerPixel)
        left.setStyleSheet(
            "QTableWidget {"
            " background: #f8fbff;"
            " border: 1px solid #d8e2ed;"
            " border-radius: 12px;"
            " gridline-color: #e7edf3;"
            "}"
            "QHeaderView::section {"
            " background: #f1f5fb;"
            " border: 1px solid #d8e2ed;"
            " padding: 4px;"
            " font-weight: 600;"
            " color: #58697b;"
            "}"
        )
        info_font = left.font()
        info_metrics = QtGui.QFontMetrics(info_font)
        max_label_width = max(info_metrics.horizontalAdvance(lbl) for lbl in row_labels)
        desired_info_width = max(260, min(520, max_label_width + 28))
        left.setColumnWidth(0, desired_info_width - 18)
        left.setMinimumWidth(desired_info_width)
        left.setMaximumWidth(desired_info_width + 18)

        right = ZoomableAlignmentTable(row_count, max_len)
        right.setHorizontalHeaderLabels([str(i + 1) for i in range(max_len)])
        right_hh = right.horizontalHeader()
        right_hh.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Fixed)
        right_hh.setMinimumSectionSize(4)
        right_hh.setDefaultSectionSize(self._col_size)
        right_hh.setStretchLastSection(False)
        right.verticalHeader().setVisible(False)
        right.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        right.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.NoSelection)
        right.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        right.setVerticalScrollMode(QtWidgets.QAbstractItemView.ScrollMode.ScrollPerPixel)
        right.setHorizontalScrollMode(QtWidgets.QAbstractItemView.ScrollMode.ScrollPerPixel)
        right.setStyleSheet(
            "QTableWidget {"
            " background: #ffffff;"
            " border: 1px solid #d8e2ed;"
            " border-radius: 12px;"
            " gridline-color: #e7edf3;"
            "}"
            "QHeaderView::section {"
            " background: #f4f7fb;"
            " border: 1px solid #e4ebf2;"
            " padding: 4px;"
            " font-weight: 600;"
            " color: #5a6b7c;"
            "}"
        )
        left.horizontalHeader().setFixedHeight(24)
        right.horizontalHeader().setFixedHeight(24)
        palette = self._residue_palette()
        left.setUpdatesEnabled(False)
        right.setUpdatesEnabled(False)
        try:
            for row_idx, label in enumerate(row_labels):
                elided = info_metrics.elidedText(
                    label,
                    QtCore.Qt.TextElideMode.ElideRight,
                    max(80, desired_info_width - 26),
                )
                item = QtWidgets.QTableWidgetItem(elided)
                item.setToolTip(label)
                if label == "consensus":
                    item.setForeground(QtGui.QBrush(QtGui.QColor("#475569")))
                    f = item.font()
                    f.setBold(True)
                    item.setFont(f)
                left.setItem(row_idx, 0, item)
                left.setRowHeight(row_idx, self._row_height)
                right.setRowHeight(row_idx, self._row_height)

            for row_idx, (_label, seq) in enumerate(padded):
                for col_idx, res in enumerate(seq):
                    bg = QtGui.QColor(palette.get(res.upper(), "#e5e7eb"))
                    item = QtWidgets.QTableWidgetItem(res)
                    item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                    item.setBackground(QtGui.QBrush(bg))
                    item.setForeground(QtGui.QBrush(QtGui.QColor("#0f172a")))
                    if res != "-" and consensus_marks[col_idx] != "*" and res != consensus_marks[col_idx]:
                        f = item.font()
                        f.setBold(True)
                        item.setFont(f)
                    right.setItem(row_idx, col_idx, item)

            consensus_row = len(padded)
            for col_idx, mark in enumerate(consensus_marks):
                item = QtWidgets.QTableWidgetItem(mark)
                item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                item.setForeground(QtGui.QBrush(QtGui.QColor("#64748b")))
                right.setItem(consensus_row, col_idx, item)
        finally:
            left.setUpdatesEnabled(True)
            right.setUpdatesEnabled(True)

        left.verticalScrollBar().valueChanged.connect(
            lambda value: self._sync_vertical_scroll(left, right, value)
        )
        right.verticalScrollBar().valueChanged.connect(
            lambda value: self._sync_vertical_scroll(right, left, value)
        )
        left.setMouseTracking(True)
        right.setMouseTracking(True)
        left.cellEntered.connect(lambda row, _col: self._on_hover_row(row))
        right.cellEntered.connect(lambda row, _col: self._on_hover_row(row))
        left.zoomRequested.connect(self._apply_zoom_step)
        right.zoomRequested.connect(self._apply_zoom_step)
        return left, right

    def _register_zoom_sources(
        self,
        left: QtWidgets.QTableWidget,
        right: QtWidgets.QTableWidget,
    ) -> None:
        # Some platforms route wheel events to headers/viewport/scrollbars rather than table.
        # Register all relevant surfaces to make wheel-zoom reliable.
        objs: List[QtCore.QObject] = [
            left,
            right,
            left.viewport(),
            right.viewport(),
            left.horizontalHeader(),
            right.horizontalHeader(),
            left.verticalHeader(),
            right.verticalHeader(),
            left.horizontalHeader().viewport(),
            right.horizontalHeader().viewport(),
            left.verticalHeader().viewport(),
            right.verticalHeader().viewport(),
            left.horizontalScrollBar(),
            right.horizontalScrollBar(),
            left.verticalScrollBar(),
            right.verticalScrollBar(),
        ]
        for obj in objs:
            if obj is None or obj in self._zoom_event_sources:
                continue
            obj.installEventFilter(self)
            self._zoom_event_sources.add(obj)

    def _sync_vertical_scroll(
        self,
        source: QtWidgets.QTableWidget,
        target: QtWidgets.QTableWidget,
        value: int,
    ) -> None:
        if self._syncing_scroll:
            return
        self._syncing_scroll = True
        try:
            if target.verticalScrollBar().value() != value:
                target.verticalScrollBar().setValue(value)
        finally:
            self._syncing_scroll = False

    def _apply_zoom_step(self, step: int) -> None:
        step_i = int(step)
        if step_i == 0:
            return
        self._zoom_pending_step += step_i
        self._zoom_timer.start()

    def _flush_zoom_steps(self) -> None:
        if self._zoom_pending_step == 0:
            return
        pending = self._zoom_pending_step
        self._zoom_pending_step = 0
        min_size = self._min_col_size()
        new_size = max(min_size, min(64, self._col_size + pending))
        if new_size == self._col_size:
            return
        self._col_size = new_size
        self._apply_zoom()

    def _min_col_size(self) -> int:
        if self._seq_table is None or self._seq_col_count <= 0:
            return 4
        viewport_w = max(1, self._seq_table.viewport().width())
        # Do not shrink so much that the full alignment becomes narrower than the visible container.
        return max(4, int((viewport_w + self._seq_col_count - 1) / self._seq_col_count))

    def _wheel_step(self, event: QtGui.QWheelEvent) -> int:
        angle = event.angleDelta()
        primary = angle.y() if abs(angle.y()) >= abs(angle.x()) else angle.x()
        if primary == 0:
            pixel = event.pixelDelta()
            primary = pixel.y() if abs(pixel.y()) >= abs(pixel.x()) else pixel.x()
        if primary == 0:
            return 0
        # Standard wheel notch is ~120; scale smoothly for trackpads.
        magnitude = max(1, int(abs(primary) / 60))
        return magnitude * (2 if primary > 0 else -2)

    def eventFilter(self, obj: QtCore.QObject, event: QtCore.QEvent) -> bool:
        if obj in {self._info_table.viewport() if self._info_table else None, self._seq_table.viewport() if self._seq_table else None}:
            if event.type() == QtCore.QEvent.Type.Leave:
                self._set_hover_sample(None)
        if obj in self._zoom_event_sources and event.type() == QtCore.QEvent.Type.Wheel:
            wheel_event = event  # type: ignore[assignment]
            step = self._wheel_step(wheel_event)
            if step != 0:
                self._apply_zoom_step(step)
                wheel_event.accept()
                return True
        return super().eventFilter(obj, event)

    def _apply_zoom(self) -> None:
        if self._info_table is None or self._seq_table is None:
            return
        self._col_size = max(self._min_col_size(), self._col_size)
        if self._last_applied_col_size == self._col_size:
            return
        hh = self._seq_table.horizontalHeader()
        hh.setDefaultSectionSize(self._col_size)
        hh.setMinimumSectionSize(4)
        # Resize all sections in C++ once instead of Python-looping all columns.
        hh.resizeSections(QtWidgets.QHeaderView.ResizeMode.Fixed)
        self._seq_table.viewport().update()
        self._last_applied_col_size = self._col_size

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        super().resizeEvent(event)
        # Re-apply bounds after resize because the minimum horizontal scale depends on viewport width.
        self._apply_zoom()
        self._schedule_row_stretch()

    def _schedule_row_stretch(self) -> None:
        self._row_stretch_timer.start()

    def _fit_tables_height(
        self,
        left: QtWidgets.QTableWidget,
        right: QtWidgets.QTableWidget,
    ) -> None:
        row_total = sum(left.rowHeight(r) for r in range(left.rowCount()))
        header_h = left.horizontalHeader().height()
        frame_h = left.frameWidth() * 2
        content_h = row_total + header_h + frame_h + 4
        screen_h = QtGui.QGuiApplication.primaryScreen().availableGeometry().height()
        cap_h = int(screen_h * 0.72)
        base_h = min(content_h, cap_h)

        # Keep a sensible minimum so empty space is reduced initially,
        # but allow both panes to expand when the dialog is resized taller.
        left.setMinimumHeight(base_h)
        right.setMinimumHeight(base_h)
        left.setMaximumHeight(16777215)
        right.setMaximumHeight(16777215)

        if content_h <= cap_h:
            left.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded)
            right.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        else:
            left.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded)
            right.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._schedule_row_stretch()

    def _stretch_rows_to_fill_viewport(self) -> None:
        if self._info_table is None or self._seq_table is None or self._row_count <= 0:
            return
        left_vh = max(0, self._info_table.viewport().height())
        right_vh = max(0, self._seq_table.viewport().height())
        viewport_h = min(left_vh, right_vh)
        if viewport_h <= 0:
            return

        # Only stretch upward when there is extra room; keep current baseline otherwise.
        fill_row_h = max(self._row_height, int(viewport_h / self._row_count))
        fill_row_h = min(64, fill_row_h)
        for row in range(self._row_count):
            self._info_table.setRowHeight(row, fill_row_h)
            self._seq_table.setRowHeight(row, fill_row_h)

    def _consensus_marks(self, seqs: List[str]) -> List[str]:
        if not seqs:
            return []
        max_len = max(len(s) for s in seqs)
        marks: List[str] = []
        for i in range(max_len):
            col = [s[i] if i < len(s) else "-" for s in seqs]
            non_gap = [c for c in col if c != "-"]
            if not non_gap:
                marks.append("-")
                continue
            if len(set(non_gap)) == 1 and len(non_gap) == len(col):
                marks.append("*")
            else:
                marks.append(Counter(non_gap).most_common(1)[0][0])
        return marks

    def _residue_palette(self) -> Dict[str, str]:
        return {
            "A": "#b7e4c7",
            "C": "#ffd6a5",
            "D": "#a9def9",
            "E": "#cdb4db",
            "F": "#ffc8dd",
            "G": "#d9ed92",
            "H": "#bde0fe",
            "I": "#fde4cf",
            "K": "#f1c0e8",
            "L": "#caffbf",
            "M": "#ffcad4",
            "N": "#e4c1f9",
            "P": "#fef9c3",
            "Q": "#bee1e6",
            "R": "#fecdd3",
            "S": "#c7d2fe",
            "T": "#bae6fd",
            "V": "#dcfce7",
            "W": "#fef08a",
            "Y": "#fde68a",
            "-": "#f3f4f6",
        }

    def _sample_from_alignment_header(self, header: str) -> Optional[str]:
        text = (header or "").strip()
        if not text or text.lower() == "consensus":
            return None
        if text.lower().startswith("wild_type"):
            return None
        parts = [p.strip() for p in text.split("|")]
        if len(parts) >= 2 and parts[1]:
            return parts[1]
        return None

    def _on_hover_row(self, row: int) -> None:
        if row < 0 or row >= len(self._row_samples):
            self._set_hover_sample(None)
            return
        self._set_hover_sample(self._row_samples[row])

    def _set_hover_sample(self, sample: Optional[str]) -> None:
        if self._hovered_sample == sample:
            return
        self._hovered_sample = sample
        if self._sample_hover_callback is not None:
            self._sample_hover_callback(sample)

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:  # type: ignore[override]
        self._set_hover_sample(None)
        super().closeEvent(event)

    def _save_alignment(self) -> None:
        dest, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save alignment",
            "alignment.fasta",
            "FASTA files (*.fasta *.fa *.faa);;All files (*)",
        )
        if not dest:
            return
        try:
            Path(dest).write_text(self._alignment_text, encoding="utf-8")
            QtWidgets.QMessageBox.information(self, "Alignment", f"Saved alignment to {dest}")
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, "Alignment", f"Could not save alignment: {exc}")


class GeneCardAlignmentDialog(QtWidgets.QDialog):
    def __init__(
        self,
        parent: Optional[QtWidgets.QWidget],
        *,
        row_data: Dict[str, Any],
        ann_df: Optional[pd.DataFrame],
        sample_hover_callback: Optional[Callable[[Optional[str]], None]] = None,
    ) -> None:
        super().__init__(parent)
        self._row_data = row_data
        self._ann_df = ann_df
        self._sample_hover_callback = sample_hover_callback
        self._floating_dialogs: List[QtWidgets.QDialog] = []
        self.setWindowTitle("Align Sequences")
        self.resize(680, 420)
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)

        desc = QtWidgets.QLabel(
            "Choose which amino acid sequences to align for this annotation."
        )
        desc.setWordWrap(True)
        layout.addWidget(desc)

        self.available_chk = QtWidgets.QCheckBox("Align with available sequences")
        self.available_chk.setChecked(True)
        layout.addWidget(self.available_chk)

        self.wt_chk = QtWidgets.QCheckBox("Add wild type")
        self.wt_chk.setChecked(False)
        self.wt_chk.toggled.connect(self._toggle_wild_type)
        layout.addWidget(self.wt_chk)

        self.wt_label = QtWidgets.QLabel("Paste wild type amino acid sequence")
        self.wt_edit = QtWidgets.QPlainTextEdit()
        self.wt_edit.setPlaceholderText("Paste raw amino acid sequence or FASTA entry")
        self.wt_edit.setLineWrapMode(QtWidgets.QPlainTextEdit.LineWrapMode.NoWrap)
        self.wt_edit.setFixedHeight(120)
        self.wt_label.setVisible(False)
        self.wt_edit.setVisible(False)
        layout.addWidget(self.wt_label)
        layout.addWidget(self.wt_edit)

        self.status_label = QtWidgets.QLabel("")
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)

        btn_row = QtWidgets.QHBoxLayout()
        btn_row.addStretch(1)
        run_btn = QtWidgets.QPushButton("Run Alignment")
        run_btn.clicked.connect(self._run_alignment)
        close_btn = QtWidgets.QPushButton("Close")
        close_btn.clicked.connect(self.reject)
        btn_row.addWidget(run_btn)
        btn_row.addWidget(close_btn)
        layout.addLayout(btn_row)

    def _toggle_wild_type(self, checked: bool) -> None:
        self.wt_label.setVisible(bool(checked))
        self.wt_edit.setVisible(bool(checked))

    def _collect_available_sequences(self) -> List[Tuple[str, str]]:
        out: List[Tuple[str, str]] = []
        seen: set[Tuple[str, str]] = set()

        def _append_row(rec: Dict[str, Any]) -> None:
            seq = _clean_sequence_text(rec.get("Amino_Acid_Sequence"))
            if not seq:
                return
            ann_id = _safe_header_part(rec.get("AnnotationID")) or "ANN"
            sample = _safe_header_part(rec.get("Sample")) or "sample"
            contig = _safe_header_part(rec.get("Contig")) or "contig"
            gene = _safe_header_part(rec.get("Gene")) or "gene"
            header = f"{ann_id}|{sample}|{contig}|{gene}"
            key = (header, seq)
            if key in seen:
                return
            seen.add(key)
            out.append((header, seq))

        # Always include the selected annotation sequence first when available.
        _append_row(self._row_data)

        if self._ann_df is None or self._ann_df.empty:
            return out
        if "Amino_Acid_Sequence" not in self._ann_df.columns:
            return out

        gene_id = _safe_value(self._row_data.get("GeneID"))
        gene = _safe_value(self._row_data.get("Gene"))
        ann = self._ann_df
        subset = ann.iloc[0:0]
        if gene_id and "GeneID" in ann.columns:
            mask = ann["GeneID"].fillna("").astype(str).str.strip() == gene_id
            subset = ann[mask]
        elif gene and "Gene" in ann.columns:
            mask = ann["Gene"].fillna("").astype(str).str.strip().str.lower() == gene.lower()
            subset = ann[mask]

        for _, row in subset.iterrows():
            _append_row(row.to_dict())
        return out

    def _find_aligner(self) -> Optional[str]:
        from pipeline.ToolEnv import runtime_is_container

        if runtime_is_container():
            return "mafft"  # provided by the tools image
        for cmd in ["mafft", "muscle", "clustalo"]:
            if shutil.which(cmd):
                return cmd
        return None

    def _perform_alignment(self, records: List[Tuple[str, str]]) -> Tuple[str, str]:
        tmpdir = Path(tempfile.mkdtemp(prefix="gene_card_align_"))
        in_path = tmpdir / "input.fasta"
        out_path = tmpdir / "alignment.fasta"
        with in_path.open("w", encoding="utf-8") as fh:
            for header, seq in records:
                fh.write(f">{header}\n{seq}\n")

        aligner = self._find_aligner()
        if not aligner:
            text = in_path.read_text(encoding="utf-8")
            return text, "No aligner (mafft/muscle/clustalo) found. Showing unaligned FASTA."

        try:
            if aligner == "mafft":
                from pipeline.ToolEnv import gui_tool_cmd

                cmd = [*gui_tool_cmd("mafft", "mafft", mount_paths=[str(in_path.parent)]), "--auto", str(in_path)]
                res = subprocess.run(cmd, capture_output=True, text=True)
                if res.returncode == 0 and res.stdout.strip():
                    return res.stdout, "Alignment completed with mafft."
                return in_path.read_text(encoding="utf-8"), f"mafft failed: {res.stderr.strip() or 'unknown error'}"

            if aligner == "muscle":
                res = subprocess.run(
                    [aligner, "-in", str(in_path), "-out", str(out_path)],
                    capture_output=True,
                    text=True,
                )
                if res.returncode == 0 and out_path.exists():
                    return out_path.read_text(encoding="utf-8"), "Alignment completed with muscle."
                return in_path.read_text(encoding="utf-8"), f"muscle failed: {res.stderr.strip() or 'unknown error'}"

            if aligner == "clustalo":
                res = subprocess.run(
                    [aligner, "-i", str(in_path), "-o", str(out_path), "--force"],
                    capture_output=True,
                    text=True,
                )
                if res.returncode == 0 and out_path.exists():
                    return out_path.read_text(encoding="utf-8"), "Alignment completed with clustalo."
                return in_path.read_text(encoding="utf-8"), f"clustalo failed: {res.stderr.strip() or 'unknown error'}"

        except Exception as exc:
            return in_path.read_text(encoding="utf-8"), f"Alignment failed: {exc}"

        return in_path.read_text(encoding="utf-8"), "Alignment failed unexpectedly."

    def _run_alignment(self) -> None:
        records: List[Tuple[str, str]] = []
        if self.available_chk.isChecked():
            records.extend(self._collect_available_sequences())

        if self.wt_chk.isChecked():
            wt_seq = _parse_pasted_sequence(self.wt_edit.toPlainText())
            if not wt_seq:
                QtWidgets.QMessageBox.warning(self, "Align", "Wild type is enabled, but no valid sequence was pasted.")
                return
            records.append(("wild_type", wt_seq))

        if not records:
            QtWidgets.QMessageBox.warning(self, "Align", "No sequences selected for alignment.")
            return

        if len(records) == 1:
            QtWidgets.QMessageBox.information(
                self,
                "Align",
                "Only one sequence is available. Add wild type or enable available sequences for MSA.",
            )
            return

        alignment_text, status = self._perform_alignment(records)
        self.status_label.setText(f"{status} ({len(records)} sequences)")
        dlg = AlignmentResultDialog(
            None,
            alignment_text=alignment_text,
            status_text=status,
            sample_hover_callback=self._sample_hover_callback,
        )
        _show_floating_dialog(self._floating_dialogs, dlg)


class AnnotationDetailDialog(QtWidgets.QDialog):
    def __init__(
        self,
        parent: Optional[QtWidgets.QWidget],
        row_data: Dict[str, Any],
        column_order: Sequence[str],
        ann_df: Optional[pd.DataFrame] = None,
        meta_df: Optional[pd.DataFrame] = None,
        marker_hits_df: Optional[pd.DataFrame] = None,
        output_dir: Optional[Path] = None,
        pan_reference_sample: Optional[str] = None,
        pan_fallback_tree: Optional[Path] = None,
        sample_hover_callback: Optional[Callable[[Optional[str]], None]] = None,
    ) -> None:
        super().__init__(parent)
        self._row_data = row_data
        self._column_order = [c for c in column_order if c in row_data]
        self._ann_df = ann_df
        self._meta_df = meta_df if meta_df is not None else pd.DataFrame()
        self._marker_hits_df = marker_hits_df if marker_hits_df is not None else pd.DataFrame()
        self._output_dir = Path(output_dir) if output_dir is not None else None
        self._pan_reference_sample = str(pan_reference_sample).strip() if pan_reference_sample else None
        self._pan_fallback_tree = Path(pan_fallback_tree) if pan_fallback_tree is not None else None
        self._sample_hover_callback = sample_hover_callback
        self._descriptor_labels: List[QtWidgets.QLabel] = []
        self._floating_dialogs: List[QtWidgets.QDialog] = []
        self._context_request_token = 0
        self._context_worker: Optional[_ContextCollectionWorker] = None
        self._context_progress: Optional[QtWidgets.QProgressDialog] = None
        self._ann_sample_norm_cache: Optional[pd.Series] = None
        self._ann_contig_norm_cache: Optional[pd.Series] = None
        self._ann_gene_norm_cache: Optional[pd.Series] = None
        self._ann_annid_norm_cache: Optional[pd.Series] = None
        self._ann_start_num_cache: Optional[pd.Series] = None
        self._ann_stop_num_cache: Optional[pd.Series] = None
        self._context_contig_cache: Dict[Tuple[str, str], pd.DataFrame] = {}
        self._sample_metadata_cache: Dict[str, Dict[str, str]] = {}
        self.setWindowTitle("Annotation details")
        self.resize(820, 620)
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)

        tabs = QtWidgets.QTabWidget()
        tabs.addTab(self._build_summary_tab(), "Summary")
        tabs.addTab(self._build_raw_tab(), "Raw")
        layout.addWidget(tabs)

        btns = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Close)
        context_btn = btns.addButton("Context", QtWidgets.QDialogButtonBox.ButtonRole.ActionRole)
        align_btn = btns.addButton("Align", QtWidgets.QDialogButtonBox.ButtonRole.ActionRole)
        context_btn.clicked.connect(self._open_context_dialog)
        align_btn.clicked.connect(self._open_alignment_dialog)
        btns.rejected.connect(self.reject)
        btns.accepted.connect(self.accept)
        layout.addWidget(btns)

    def _build_summary_tab(self) -> QtWidgets.QWidget:
        widget = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(widget)

        header = QtWidgets.QLabel(self._format_header())
        header.setWordWrap(True)
        header.setStyleSheet("font-weight: 600;")
        layout.addWidget(header)

        pos_group = QtWidgets.QGroupBox("Position")
        pos_layout = QtWidgets.QFormLayout(pos_group)
        pos_layout.addRow(self._make_descriptor_label("Position"), QtWidgets.QLabel(self._format_position()))
        pos_layout.addRow(self._make_descriptor_label("Orientation"), QtWidgets.QLabel(self._format_orientation()))
        layout.addWidget(pos_group)

        ids_group = QtWidgets.QGroupBox("Identifiers")
        ids_layout = QtWidgets.QFormLayout(ids_group)
        for label, key in [("NCBI", "NCBI_NRP_ID"), ("RefSeq", "Refseq_ID")]:
            value = _safe_value(self._row_data.get(key))
            if value:
                ids_layout.addRow(self._make_descriptor_label(label), QtWidgets.QLabel(value))
        if ids_layout.rowCount() == 0:
            ids_layout.addRow(self._make_descriptor_label("(none)"), QtWidgets.QLabel(""))
        layout.addWidget(ids_group)

        seq_group = QtWidgets.QGroupBox("Sequences")
        seq_layout = QtWidgets.QFormLayout(seq_group)
        nt_seq = _safe_value(self._row_data.get("Nucleotide_Sequence"))
        aa_seq = _safe_value(self._row_data.get("Amino_Acid_Sequence"))
        if nt_seq:
            seq_layout.addRow(self._make_descriptor_label("Nucleotide"), self._build_sequence_row(nt_seq))
        if aa_seq:
            seq_layout.addRow(self._make_descriptor_label("Amino acid"), self._build_sequence_row(aa_seq))
        if seq_layout.rowCount() == 0:
            seq_layout.addRow(self._make_descriptor_label("(none)"), QtWidgets.QLabel(""))
        layout.addWidget(seq_group)

        marker_hits = self._marker_hits_for_row()
        if not marker_hits.empty:
            markers_group = QtWidgets.QGroupBox("Markers")
            markers_layout = QtWidgets.QVBoxLayout(markers_group)
            table = QtWidgets.QTableWidget()
            table.setColumnCount(5)
            table.setHorizontalHeaderLabels(["DB", "Method", "Coverage", "Identity", "Product"])
            table.setRowCount(len(marker_hits))
            table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
            table.horizontalHeader().setStretchLastSection(True)
            table.verticalHeader().setVisible(False)
            for row_idx, marker in enumerate(marker_hits.to_dict("records")):
                db = str(marker.get("DB") or "")
                method = str(marker.get("Method") or "")
                product = str(marker.get("Product") or "")
                coverage = marker.get("Subject_coverage") or ""
                identity = marker.get("Identity") or ""
                table.setItem(row_idx, 0, QtWidgets.QTableWidgetItem(db))
                table.setItem(row_idx, 1, QtWidgets.QTableWidgetItem(method))
                table.setItem(row_idx, 4, QtWidgets.QTableWidgetItem(product))
                self._set_metric_cell(table, row_idx, 2, coverage, "#7bb6ff")
                self._set_metric_cell(table, row_idx, 3, identity, "#ff8f8f")
            table.resizeColumnsToContents()
            markers_layout.addWidget(table)
            layout.addWidget(markers_group)

        layout.addStretch(1)
        self._update_descriptor_widths()
        return widget

    def _build_raw_tab(self) -> QtWidgets.QWidget:
        widget = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(widget)
        raw = QtWidgets.QPlainTextEdit()
        raw.setReadOnly(True)
        normalized = {
            k: _normalize_detail_value(v)
            for k, v in self._row_data.items()
            if k not in _HIDDEN_INTERNAL_MARKER_FIELDS
        }
        marker_hits = self._marker_hits_for_row()
        if not marker_hits.empty:
            visible_cols = [c for c in marker_hits.columns if c not in _HIDDEN_INTERNAL_MARKER_FIELDS]
            normalized["MarkerHits"] = _normalize_detail_value(marker_hits[visible_cols].to_dict("records"))
        raw.setPlainText(json.dumps(normalized, indent=2, ensure_ascii=False))
        layout.addWidget(raw)
        return widget

    def _format_header(self) -> str:
        sample = _safe_value(self._row_data.get("Sample")) or "Unknown sample"
        contig_raw = _safe_value(self._row_data.get("Contig")) or "Unknown contig"
        contig = contig_raw if contig_raw.lower().startswith("contig") else f"Contig {contig_raw}"
        gene = _safe_value(self._row_data.get("Gene")) or "Unknown gene"
        return f"{sample} > {contig} > {gene}"

    def _format_position(self) -> str:
        start = _safe_value(self._row_data.get("Start"))
        stop = _safe_value(self._row_data.get("Stop"))
        if start and stop:
            arrow = "→"
            strand = _safe_value(self._row_data.get("Strand"))
            if strand == "-":
                arrow = "←"
            return f"{start} {arrow} {stop}"
        return start or stop or ""

    def _format_orientation(self) -> str:
        strand = _safe_value(self._row_data.get("Strand")) or "?"
        if strand == "+":
            return "Forward (+)"
        if strand == "-":
            return "Reverse (-)"
        return strand

    def _build_sequence_row(self, seq: str) -> QtWidgets.QWidget:
        widget = QtWidgets.QWidget()
        layout = QtWidgets.QHBoxLayout(widget)
        layout.setContentsMargins(0, 0, 0, 0)
        preview = self._sequence_preview(seq)
        label = QtWidgets.QLabel(preview)
        label.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        label.setToolTip(seq)
        copy_btn = QtWidgets.QPushButton("Copy")
        copy_btn.clicked.connect(lambda _checked=False, s=seq: self._copy_text(s))
        layout.addWidget(label, stretch=1)
        layout.addWidget(copy_btn)
        return widget

    def _sequence_preview(self, seq: str, max_len: int = 80) -> str:
        if len(seq) <= max_len:
            return f"{seq} (len={len(seq)})"
        return f"{seq[:max_len]}… (len={len(seq)})"

    def _copy_text(self, text: str) -> None:
        QtGui.QGuiApplication.clipboard().setText(text)

    def _marker_hits_for_row(self) -> pd.DataFrame:
        df = self._marker_hits_df
        if df is None or df.empty:
            return pd.DataFrame()
        if "AnnotationID" not in df.columns:
            return pd.DataFrame()

        ann_id = _safe_value(self._row_data.get("AnnotationID"))
        if not ann_id:
            # Strict linking: marker hits are keyed by Bakta AnnotationID (BKT...).
            return pd.DataFrame()

        out = df[df["AnnotationID"].astype(str).str.strip() == ann_id].copy()
        if out.empty:
            return pd.DataFrame()

        # Guard against stale/invalid cached rows that contain no marker payload.
        payload_cols = [c for c in ("DB", "Method", "Product") if c in out.columns]
        if payload_cols:
            keep = pd.Series(False, index=out.index)
            for col in payload_cols:
                vals = out[col].astype(str).str.strip()
                keep |= vals.ne("") & vals.ne("nan") & vals.ne("<NA>")
            out = out[keep]

        if out.empty:
            return pd.DataFrame()

        sort_cols: List[str] = []
        ascending: List[bool] = []
        if "Bitscore" in out.columns:
            sort_cols.append("Bitscore")
            ascending.append(False)
        if "Identity" in out.columns:
            sort_cols.append("Identity")
            ascending.append(False)
        if sort_cols:
            out = out.sort_values(sort_cols, ascending=ascending, na_position="last")
        return out.reset_index(drop=True)

    def _set_metric_cell(
        self, table: QtWidgets.QTableWidget, row: int, col: int, value: Any, color: str
    ) -> None:
        metric = _parse_float(value)
        if metric is None:
            table.setItem(row, col, QtWidgets.QTableWidgetItem(""))
            return
        metric = max(0.0, min(100.0, metric))
        bar = QtWidgets.QProgressBar()
        bar.setRange(0, 100)
        bar.setValue(int(round(metric)))
        bar.setFormat(f"{metric:.1f}")
        bar.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        bar.setTextVisible(True)
        bar.setStyleSheet(
            "QProgressBar {"
            " border: 1px solid #d3dbe6;"
            " border-radius: 7px;"
            " background: #f3f6fa;"
            " color: #223246;"
            " text-align: center;"
            "}"
            f"QProgressBar::chunk {{ background-color: {color}; border-radius: 7px; }}"
        )
        table.setCellWidget(row, col, bar)
        table.setItem(row, col, QtWidgets.QTableWidgetItem(f"{metric:.1f}"))

    def _make_descriptor_label(self, text: str) -> QtWidgets.QLabel:
        label = QtWidgets.QLabel(text)
        label.setStyleSheet("font-weight: 600; color: #2b3c4d;")
        label.setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignVCenter)
        self._descriptor_labels.append(label)
        return label

    def _update_descriptor_widths(self) -> None:
        if not self._descriptor_labels:
            return
        width = max(120, int(self.width() * 0.15))
        for label in self._descriptor_labels:
            label.setFixedWidth(width)

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        super().resizeEvent(event)
        self._update_descriptor_widths()

    def _open_alignment_dialog(self) -> None:
        dlg = GeneCardAlignmentDialog(
            None,
            row_data=self._row_data,
            ann_df=self._ann_df,
            sample_hover_callback=self._sample_hover_callback,
        )
        _show_floating_dialog(self._floating_dialogs, dlg)

    def _open_context_dialog(self) -> None:
        if self._context_worker is not None:
            if self._context_progress is not None:
                self._context_progress.raise_()
                self._context_progress.activateWindow()
            return
        self._context_request_token += 1
        token = self._context_request_token

        progress = QtWidgets.QProgressDialog("Loading gene context…", "", 0, 0, self)
        progress.setWindowTitle("Context")
        progress.setCancelButton(None)
        progress.setMinimumDuration(0)
        progress.setAutoClose(False)
        progress.setAutoReset(False)
        progress.setWindowModality(QtCore.Qt.WindowModality.WindowModal)
        progress.show()
        progress.raise_()
        progress.activateWindow()
        self._context_progress = progress

        worker = _ContextCollectionWorker(
            token,
            lambda: self._build_context_collection_rows(window_genes=15),
        )
        worker.signals.finished.connect(self._on_context_collection_ready)
        worker.signals.failed.connect(self._on_context_collection_failed)
        self._context_worker = worker
        QtCore.QThreadPool.globalInstance().start(worker)

    def _finish_context_collection_request(self, token: int) -> bool:
        if token != self._context_request_token:
            return False
        self._context_worker = None
        if self._context_progress is not None:
            self._context_progress.close()
            self._context_progress.deleteLater()
            self._context_progress = None
        return True

    def _on_context_collection_ready(self, token: int, rows_obj: object) -> None:
        if not self._finish_context_collection_request(token):
            return
        rows = rows_obj if isinstance(rows_obj, list) else []
        if not rows:
            QtWidgets.QMessageBox.information(
                self,
                "Context",
                "No contexts found for this gene.",
            )
            return
        gene = _safe_value(self._row_data.get("Gene")) or "gene"
        dlg = IgvCollectionDialog(
            None,
            rows=rows,
            title=f"IGV_collection - {gene}",
            output_dir=self._output_dir,
            pan_reference_sample=self._pan_reference_sample,
            pan_fallback_tree=self._pan_fallback_tree,
        )
        dlg.geneDoubleClicked.connect(self._open_collection_gene_card)
        _show_floating_dialog(self._floating_dialogs, dlg)

    def _on_context_collection_failed(self, token: int, error_text: str) -> None:
        if not self._finish_context_collection_request(token):
            return
        QtWidgets.QMessageBox.warning(
            self,
            "Context",
            "Failed to build the gene context view.\n\n"
            f"{error_text.strip().splitlines()[-1] if error_text.strip() else 'Unknown error.'}",
        )

    def _ensure_context_lookup_caches(self) -> bool:
        df = self._ann_df
        if df is None or df.empty:
            return False
        if "Sample" not in df.columns or "Contig" not in df.columns:
            return False

        if self._ann_sample_norm_cache is None:
            self._ann_sample_norm_cache = df["Sample"].astype(str).str.strip()
        if self._ann_contig_norm_cache is None:
            self._ann_contig_norm_cache = df["Contig"].astype(str).str.strip()
        if self._ann_gene_norm_cache is None:
            if "Gene" in df.columns:
                self._ann_gene_norm_cache = df["Gene"].astype(str).str.strip().str.lower()
            else:
                self._ann_gene_norm_cache = pd.Series("", index=df.index, dtype=str)
        if self._ann_annid_norm_cache is None:
            if "AnnotationID" in df.columns:
                self._ann_annid_norm_cache = df["AnnotationID"].astype(str).str.strip()
            else:
                self._ann_annid_norm_cache = pd.Series("", index=df.index, dtype=str)
        if self._ann_start_num_cache is None:
            start_source = df["StartInt"] if "StartInt" in df.columns else df.get("Start")
            self._ann_start_num_cache = pd.to_numeric(start_source, errors="coerce")
        if self._ann_stop_num_cache is None:
            stop_source = df["StopInt"] if "StopInt" in df.columns else df.get("Stop")
            self._ann_stop_num_cache = pd.to_numeric(stop_source, errors="coerce")
        return True

    def _context_contig_rows(self, sample: str, contig: str) -> pd.DataFrame:
        key = (str(sample).strip(), str(contig).strip())
        cached = self._context_contig_cache.get(key)
        if cached is not None:
            return cached
        if not self._ensure_context_lookup_caches():
            empty = pd.DataFrame()
            self._context_contig_cache[key] = empty
            return empty

        assert self._ann_df is not None
        assert self._ann_sample_norm_cache is not None
        assert self._ann_contig_norm_cache is not None
        assert self._ann_start_num_cache is not None
        assert self._ann_stop_num_cache is not None
        assert self._ann_gene_norm_cache is not None
        assert self._ann_annid_norm_cache is not None

        mask = (self._ann_sample_norm_cache == key[0]) & (self._ann_contig_norm_cache == key[1])
        if not bool(mask.any()):
            empty = pd.DataFrame()
            self._context_contig_cache[key] = empty
            return empty

        contig_df = self._ann_df.loc[mask].copy()
        contig_df["__start"] = self._ann_start_num_cache.loc[contig_df.index]
        contig_df["__stop"] = self._ann_stop_num_cache.loc[contig_df.index]
        contig_df["__gene_norm"] = self._ann_gene_norm_cache.loc[contig_df.index]
        contig_df["__annid_norm"] = self._ann_annid_norm_cache.loc[contig_df.index]
        contig_df = contig_df.dropna(subset=["__start", "__stop"])
        if contig_df.empty:
            self._context_contig_cache[key] = contig_df
            return contig_df
        contig_df["__start"] = contig_df["__start"].astype(int)
        contig_df["__stop"] = contig_df["__stop"].astype(int)
        sort_cols = ["__start", "__stop"]
        if "__gene_norm" in contig_df.columns:
            sort_cols.append("__gene_norm")
        if "Product" in contig_df.columns:
            sort_cols.append("Product")
        contig_df = contig_df.sort_values(sort_cols, na_position="last").reset_index(drop=True)
        self._context_contig_cache[key] = contig_df
        return contig_df

    def _build_context_collection_rows(self, window_genes: int = 15) -> List[Dict[str, Any]]:
        df = self._ann_df
        if df is None or df.empty:
            return []
        if "Sample" not in df.columns or "Contig" not in df.columns:
            return []
        if not self._ensure_context_lookup_caches():
            return []

        assert self._ann_gene_norm_cache is not None
        assert self._ann_annid_norm_cache is not None

        gene_key = _safe_value(self._row_data.get("Gene")).lower()
        ann_id_key = _safe_value(self._row_data.get("AnnotationID"))
        sample_key = _safe_value(self._row_data.get("Sample"))
        contig_key = _safe_value(self._row_data.get("Contig"))
        start_key = self._to_int(self._row_data.get("StartInt"), self._row_data.get("Start"))
        stop_key = self._to_int(self._row_data.get("StopInt"), self._row_data.get("Stop"))

        candidate_index = pd.Index([])
        if gene_key:
            candidate_index = self._ann_gene_norm_cache[self._ann_gene_norm_cache == gene_key].index
        if candidate_index.empty and ann_id_key:
            candidate_index = self._ann_annid_norm_cache[self._ann_annid_norm_cache == ann_id_key].index
        if candidate_index.empty:
            return []
        candidates = df.loc[candidate_index]
        if candidates.empty:
            return []

        rows: List[Tuple[int, Dict[str, Any]]] = []
        seen: set[Tuple[str, str, int, int, str]] = set()
        for _, hit_row in candidates.iterrows():
            sample = _safe_value(hit_row.get("Sample"))
            contig = _safe_value(hit_row.get("Contig"))
            if not sample or not contig:
                continue

            contig_df = self._context_contig_rows(sample, contig)
            if contig_df.empty:
                continue

            hit_start = self._to_int(hit_row.get("StartInt"), hit_row.get("Start"))
            hit_stop = self._to_int(hit_row.get("StopInt"), hit_row.get("Stop"))
            if hit_start is None or hit_stop is None:
                continue
            hit_ann = _safe_value(hit_row.get("AnnotationID"))
            hit_gene = _safe_value(hit_row.get("Gene")).lower()

            mask = (contig_df["__start"] == int(hit_start)) & (contig_df["__stop"] == int(hit_stop))
            if hit_ann:
                mask &= contig_df["__annid_norm"] == hit_ann
            if not bool(mask.any()) and hit_gene:
                mask = (
                    (contig_df["__gene_norm"] == hit_gene)
                    & (contig_df["__start"] == int(hit_start))
                    & (contig_df["__stop"] == int(hit_stop))
                )
            if bool(mask.any()):
                hit_idx = int(mask[mask].index[0])
            else:
                mid = int(round((int(hit_start) + int(hit_stop)) / 2.0))
                mid_series = ((contig_df["__start"] + contig_df["__stop"]) / 2.0 - float(mid)).abs()
                hit_idx = int(mid_series.idxmin())

            center_row = contig_df.iloc[hit_idx]
            center_start = int(center_row["__start"])
            center_stop = int(center_row["__stop"])
            center_ann = _safe_value(center_row.get("AnnotationID"))
            key = (sample, contig, center_start, center_stop, center_ann)
            if key in seen:
                continue
            seen.add(key)

            lo = max(0, hit_idx - int(window_genes))
            hi = min(len(contig_df) - 1, hit_idx + int(window_genes))
            center_bp = int(round((center_start + center_stop) / 2.0))

            genes: List[Dict[str, Any]] = []
            for pos in range(lo, hi + 1):
                row = contig_df.iloc[pos]
                start = int(row["__start"])
                stop = int(row["__stop"])
                label = _safe_value(row.get("Gene")) or _safe_value(row.get("Product")) or "gene"
                genes.append(
                    {
                        "start": start,
                        "end": stop,
                        "strand": _safe_value(row.get("Strand")) or "+",
                        "label": label,
                        "color": self._igv_like_color(_safe_value(row.get("Tags"))),
                        "row_data": {
                            k: _normalize_detail_value(v)
                            for k, v in row.drop(
                                labels=["__start", "__stop", "__gene_norm", "__annid_norm"],
                                errors="ignore",
                            ).to_dict().items()
                        },
                        "is_center": pos == hit_idx,
                    }
                )

            is_source = (
                sample == sample_key
                and contig == contig_key
                and start_key is not None
                and stop_key is not None
                and center_start == int(start_key)
                and center_stop == int(stop_key)
            )
            sample_meta = self._sample_metadata_map(sample)
            is_plasmid = str(center_row.get("MobSuiteIsPlasmid") or "").strip().lower() in PLASMID_TRUE_VALUES
            row_payload = {
                "title": f"{sample} | {contig}",
                "sample": sample,
                "contig": contig,
                "molecule": "plasmid" if is_plasmid else "chromosome",
                "is_plasmid": is_plasmid,
                "is_source": is_source,
                "sample_meta": sample_meta,
                "meta_tooltip": self._sample_metadata_tooltip(sample, contig, sample_meta),
                "center_bp": center_bp,
                "genes": genes,
            }
            rows.append((0 if is_source else 1, row_payload))

        rows.sort(key=lambda pair: (pair[0], str(pair[1].get("title") or "")))
        return [payload for _, payload in rows]

    @staticmethod
    def _to_int(primary: Any, fallback: Any = None) -> Optional[int]:
        for value in (primary, fallback):
            if value is None:
                continue
            try:
                if isinstance(value, float) and pd.isna(value):
                    continue
            except Exception:
                pass
            try:
                return int(value)
            except Exception:
                try:
                    return int(float(value))
                except Exception:
                    continue
        return None

    @staticmethod
    def _igv_like_color(tags: str) -> str:
        value = str(tags or "")
        if "A" in value:
            return "#ef9a9a"
        if "B" in value:
            return "#a5d6a7"
        if "V" in value:
            return "#ce93d8"
        if "M" in value:
            return "#ffcc80"
        if "P" in value:
            return "#9ec3ff"
        if "I" in value:
            return "#ffe082"
        return "#90caf9"

    def _sample_metadata_map(self, sample: str) -> Dict[str, str]:
        cached = self._sample_metadata_cache.get(sample)
        if cached is not None:
            return dict(cached)
        out: Dict[str, str] = {}
        df = self._meta_df
        if df is None or df.empty or "Sample" not in df.columns:
            return out
        rows = df[df["Sample"].astype(str).str.strip() == str(sample).strip()]
        if rows.empty:
            return out
        row = rows.iloc[0]
        for col in rows.columns:
            if col == "Sample":
                continue
            value = _safe_value(row.get(col))
            if value:
                out[str(col)] = value
        self._sample_metadata_cache[sample] = dict(out)
        return out

    def _sample_metadata_tooltip(self, sample: str, contig: str, sample_meta: Optional[Dict[str, str]] = None) -> str:
        header = [f"Sample: {sample}", f"Contig: {contig}"]
        meta_map = sample_meta if sample_meta is not None else self._sample_metadata_map(sample)
        if not meta_map:
            return "\n".join(header)
        details: List[str] = []
        for col in sorted(meta_map.keys()):
            details.append(f"{col}: {meta_map[col]}")
        if not details:
            return "\n".join(header)
        return "\n".join(header + [""] + details)

    def _open_collection_gene_card(self, payload: Dict[str, Any]) -> None:
        if not payload:
            return
        row_data = {str(k): _normalize_detail_value(v) for k, v in dict(payload).items()}
        dlg = AnnotationDetailDialog(
            None,
            row_data,
            list(row_data.keys()),
            ann_df=self._ann_df,
            meta_df=self._meta_df,
            marker_hits_df=self._marker_hits_df,
            output_dir=self._output_dir,
            pan_reference_sample=self._pan_reference_sample,
            pan_fallback_tree=self._pan_fallback_tree,
            sample_hover_callback=self._sample_hover_callback,
        )
        _show_floating_dialog(self._floating_dialogs, dlg)
