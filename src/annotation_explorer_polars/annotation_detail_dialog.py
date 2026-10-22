from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
import threading
from collections import Counter
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets

from pipeline.BaktaEvidence import KEYWORD_METHOD

from .constants import PLASMID_TRUE_VALUES
from .shell.commands import PANEL_OPEN_ALIGNMENT, WINDOW_OPEN_NEIGHBOURHOOD
from .resistance_mutation import (
    MUTATION_BADGE_SVG,
    MUTATION_TAG,
    has_resistance_mutation,
    mutation_note,
    resistance_mutations,
)
from .theme import COLORS
from .windows.info_window import InfoWindow

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


def normalize_detail_row(payload: Dict[str, Any]) -> Dict[str, Any]:
    """A gene row from the neighborhood viewer, as plain values for a gene card."""
    return {str(k): _normalize_detail_value(v) for k, v in dict(payload or {}).items()}


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


class AlignmentView(QtWidgets.QWidget):
    """A gene card's alignment result: the sequence names beside the aligned
    residues, with a consensus row. The body of the ``AlignmentPanel`` of an
    explore window (``panel_mode=True``).

    ``sample_hover_callback`` gets the sample of the hovered row (None when
    the pointer leaves). In a panel the rows follow the window's sample set
    (``set_sample_scope``; the consensus is that of the whole alignment) and
    rows of samples hovered elsewhere are tinted (``set_linked_samples``).
    """

    # The Close button (the dialog closes; a panel has its tab's close button).
    closeRequested = QtCore.pyqtSignal()

    def __init__(
        self,
        parent: Optional[QtWidgets.QWidget],
        *,
        alignment_text: str,
        status_text: str,
        sample_hover_callback: Optional[Callable[[Optional[str]], None]] = None,
        panel_mode: bool = False,
    ) -> None:
        super().__init__(parent)
        self._alignment_text = alignment_text
        self._sample_hover_callback = sample_hover_callback
        self._panel_mode = bool(panel_mode)
        self._scope_samples: Optional[frozenset] = None
        self._linked_samples: frozenset = frozenset()
        # Sequence names carry each sample as _safe_header_part wrote it;
        # the window's own names, by that form, for hovering out.
        self._window_names: Dict[str, str] = {}
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
        btn_row.addWidget(save_btn)
        if not self._panel_mode:
            close_btn = QtWidgets.QPushButton("Close")
            close_btn.clicked.connect(self.closeRequested.emit)
            btn_row.addWidget(close_btn)
        layout.addLayout(btn_row)

    def alignment_text(self) -> str:
        return self._alignment_text

    def row_samples(self) -> List[Optional[str]]:
        """The sample of each row (None for the consensus and wild type rows)."""
        return list(self._row_samples)

    def set_sample_scope(self, samples: Optional[Sequence[str]]) -> None:
        """Show only the rows of ``samples`` (the window's sample set; None
        shows every row). Rows without a sample (consensus, wild type) stay."""
        scope = self._header_forms(samples) if samples is not None else None
        if scope == self._scope_samples:
            return
        self._scope_samples = scope
        for row, sample in enumerate(self._row_samples):
            hidden = scope is not None and sample is not None and sample not in scope
            for table in (self._info_table, self._seq_table):
                if table is not None:
                    table.setRowHidden(row, hidden)
        self._schedule_row_stretch()

    def is_row_shown(self, row: int) -> bool:
        table = self._info_table
        return table is not None and 0 <= row < table.rowCount() and not table.isRowHidden(row)

    def set_linked_samples(self, samples: Optional[Sequence[str]]) -> None:
        """Tint the name cells of rows whose sample is hovered elsewhere."""
        linked = self._header_forms(samples)
        if linked == self._linked_samples or self._info_table is None:
            return
        self._linked_samples = linked
        tint = QtGui.QBrush(QtGui.QColor(COLORS["primary_soft"]))
        for row, sample in enumerate(self._row_samples):
            item = self._info_table.item(row, 0)
            if item is not None:
                item.setBackground(tint if sample is not None and sample in linked else QtGui.QBrush())

    def _header_forms(self, samples: Optional[Sequence[str]]) -> frozenset:
        """``samples`` as the sequence names spell them (special characters
        become '_'), so 'ST492 isolate#1' finds its row 'ST492_isolate_1'."""
        forms = set()
        for sample in samples or ():
            raw = str(sample).strip()
            form = _safe_header_part(raw)
            if form:
                forms.add(form)
                self._window_names.setdefault(form, raw)
        return frozenset(forms)

    def linked_rows(self) -> List[int]:
        return [row for row, sample in enumerate(self._row_samples) if sample is not None and sample in self._linked_samples]

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
        shown = sum(1 for row in range(self._row_count) if not self._info_table.isRowHidden(row)) or 1
        fill_row_h = max(self._row_height, int(viewport_h / shown))
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
            # The window's spelling of the sample, when it is known.
            self._sample_hover_callback(self._window_names.get(sample, sample) if sample else sample)

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


ALIGN_NEEDS_TOOLS_TIP = "Needs the pipeline setup, or mafft installed"


def find_aligner() -> Optional[str]:
    """The aligner the gene card uses: mafft from the tools image on a container setup,
    else one on PATH. A viewer install has no tools image, so only one on PATH."""
    from pipeline.ToolEnv import runtime_is_container, viewer_mode

    if not viewer_mode() and runtime_is_container():
        return "mafft"  # provided by the tools image
    for cmd in ["mafft", "muscle", "clustalo"]:
        if shutil.which(cmd):
            return cmd
    return None


class GeneCardAlignmentDialog(QtWidgets.QDialog):
    """The gene card's alignment parameter form (attached to the card). The
    result goes to ``on_result(alignment_text, status_text)``; the card shows
    it in its owner explore window."""

    def __init__(
        self,
        parent: Optional[QtWidgets.QWidget],
        *,
        row_data: Dict[str, Any],
        ann_df: Optional[pd.DataFrame],
        on_result: Optional[Callable[[str, str], None]] = None,
    ) -> None:
        super().__init__(parent)
        self._row_data = row_data
        self._ann_df = ann_df
        self._on_result = on_result
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
        return find_aligner()

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
                from pipeline.ToolEnv import gui_tool_cmd, viewer_mode

                base = ["mafft"] if viewer_mode() else gui_tool_cmd("mafft", "mafft", mount_paths=[str(in_path.parent)])
                cmd = [*base, "--auto", str(in_path)]
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
        if self._on_result is not None:
            self._on_result(alignment_text, status)


class AnnotationDetailDialog(InfoWindow):
    """The gene card (an info window, plan section 5). It keeps its actions,
    which run as commands on its owner explore window: Compare neighborhoods
    opens the neighborhood viewer in a new explore window with the owner's
    filters, Align shows the parameter form and then the result in the owner
    (with the result's sample hover in that window)."""

    COMMANDS = (WINDOW_OPEN_NEIGHBOURHOOD, PANEL_OPEN_ALIGNMENT)

    def __init__(
        self,
        parent: Optional[QtWidgets.QWidget],
        row_data: Dict[str, Any],
        column_order: Sequence[str],
        ann_df: Optional[pd.DataFrame] = None,
        meta_df: Optional[pd.DataFrame] = None,
        marker_hits_df: Optional[pd.DataFrame] = None,
    ) -> None:
        super().__init__(parent)
        self._row_data = row_data
        self._column_order = [c for c in column_order if c in row_data]
        self._ann_df = ann_df
        self._meta_df = meta_df if meta_df is not None else pd.DataFrame()
        self._marker_hits_df = marker_hits_df if marker_hits_df is not None else pd.DataFrame()
        self._alignment_form: Optional[GeneCardAlignmentDialog] = None
        self._sequence_edits: List[QtWidgets.QLineEdit] = []
        self._marker_table: Optional[QtWidgets.QTableWidget] = None
        self._raw_editor: Optional[QtWidgets.QPlainTextEdit] = None
        # One collection at a time over this card's lookup caches (they are
        # filled lazily, on a worker thread).
        self._context_rows_lock = threading.Lock()
        self._ann_sample_norm_cache: Optional[pd.Series] = None
        self._ann_contig_norm_cache: Optional[pd.Series] = None
        self._ann_gene_norm_cache: Optional[pd.Series] = None
        self._ann_annid_norm_cache: Optional[pd.Series] = None
        self._ann_start_num_cache: Optional[pd.Series] = None
        self._ann_stop_num_cache: Optional[pd.Series] = None
        self._context_contig_cache: Dict[Tuple[str, str], pd.DataFrame] = {}
        self._sample_metadata_cache: Dict[str, Dict[str, str]] = {}
        gene = _safe_value(self._row_data.get("Gene")) or "Unnamed annotation"
        self.setObjectName("AnnotationDetailDialog")
        self.set_card_title(f"Gene annotation: {gene}")
        self.setMinimumSize(720, 540)
        self.resize(920, 700)
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 14)
        layout.setSpacing(12)

        layout.addWidget(self._build_hero())

        tabs = QtWidgets.QTabWidget()
        tabs.setObjectName("AnnotationDetailTabs")
        tabs.setDocumentMode(True)
        tabs.addTab(self._build_summary_tab(), "Overview")
        tabs.addTab(self._build_raw_tab(), "All fields")
        self._tabs = tabs
        layout.addWidget(tabs, stretch=1)

        btns = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Close)
        btns.setObjectName("AnnotationActionButtons")
        context_btn = btns.addButton(
            "Compare neighborhoods…",
            QtWidgets.QDialogButtonBox.ButtonRole.ActionRole,
        )
        context_btn.setProperty("role", "primary")
        context_btn.setToolTip(
            "Compare this gene with matching neighborhoods across samples, "
            "in a new explore window with this window's filters"
        )
        align_btn = btns.addButton(
            "Align sequences…",
            QtWidgets.QDialogButtonBox.ButtonRole.ActionRole,
        )
        align_btn.setToolTip("Align this protein with matching available sequences")
        self._context_button = context_btn
        self._align_button = align_btn
        context_btn.clicked.connect(self._open_context_dialog)
        align_btn.clicked.connect(self._open_alignment_dialog)
        self._update_action_buttons()
        btns.rejected.connect(self.reject)
        btns.accepted.connect(self.accept)

        actions = QtWidgets.QFrame()
        actions.setObjectName("AnnotationActions")
        actions_layout = QtWidgets.QHBoxLayout(actions)
        actions_layout.setContentsMargins(0, 11, 0, 0)
        actions_layout.setSpacing(10)
        action_hint = QtWidgets.QLabel("Explore this annotation")
        action_hint.setProperty("role", "annotationActionHint")
        actions_layout.addWidget(action_hint)
        actions_layout.addStretch(1)
        actions_layout.addWidget(btns)
        layout.addWidget(actions)

    def _build_hero(self) -> QtWidgets.QFrame:
        hero = QtWidgets.QFrame()
        hero.setObjectName("AnnotationHero")
        layout = QtWidgets.QVBoxLayout(hero)
        layout.setContentsMargins(18, 15, 18, 15)
        layout.setSpacing(5)

        eyebrow = QtWidgets.QLabel("GENE ANNOTATION")
        eyebrow.setProperty("role", "annotationEyebrow")
        layout.addWidget(eyebrow)

        gene = _safe_value(self._row_data.get("Gene")) or "Unnamed annotation"
        title = QtWidgets.QLabel(gene)
        title.setObjectName("AnnotationTitle")
        title.setProperty("role", "annotationTitle")
        title.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(title)

        product = _safe_value(self._row_data.get("Product")) or "No product description available"
        product_label = QtWidgets.QLabel(product)
        product_label.setObjectName("AnnotationProduct")
        product_label.setProperty("role", "annotationProduct")
        product_label.setWordWrap(True)
        product_label.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(product_label)

        sample = _safe_value(self._row_data.get("Sample")) or "Unknown sample"
        contig_raw = _safe_value(self._row_data.get("Contig")) or "Unknown contig"
        contig = contig_raw if contig_raw.lower().startswith("contig") else f"Contig {contig_raw}"
        source = QtWidgets.QLabel(f"{sample}  •  {contig}")
        source.setProperty("role", "annotationSource")
        source.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(source)

        badges = self._annotation_badges()
        if badges:
            badge_wrap = QtWidgets.QWidget()
            badge_layout = QtWidgets.QHBoxLayout(badge_wrap)
            badge_layout.setContentsMargins(0, 5, 0, 0)
            badge_layout.setSpacing(6)
            for text, variant, tooltip in badges:
                badge_layout.addWidget(self._make_annotation_badge(text, variant, tooltip))
            badge_layout.addStretch(1)
            layout.addWidget(badge_wrap)
        return hero

    def _annotation_badges(self) -> List[Tuple[str, str, str]]:
        badges: List[Tuple[str, str, str]] = []
        position = self._format_position()
        if position:
            badges.append((position, "fact", "Genomic coordinates"))
        orientation = self._format_orientation()
        if orientation:
            badges.append((orientation, "fact", "Gene orientation"))
        annotation_id = _safe_value(self._row_data.get("AnnotationID"))
        if annotation_id:
            badges.append((annotation_id, "identifier", "Annotation ID"))

        tag_labels = {
            "A": ("AMR", "amr"),
            "B": ("BMRG", "bmrg"),
            "V": ("Virulence", "virulence"),
            "M": ("Mobile", "mobile"),
            "P": ("Phage", "phage"),
            "I": ("Overlap", "overlap"),
            MUTATION_TAG: ("Resistance mutation", "amr"),
        }
        tag_text = _safe_value(self._row_data.get("Tags")).upper()
        # A tag from a keyword in Bakta's product (Mobile, Phage) says so.
        by_name = _safe_value(self._row_data.get("ConsensusSource")) == KEYWORD_METHOD
        for token, (label, variant) in tag_labels.items():
            if token in tag_text:
                tooltip = mutation_note(self._row_data) if token == MUTATION_TAG else ""
                if not tooltip and by_name and token in "MP":
                    tooltip = f"Annotation tag: {label}, from the Bakta product name (not a database search)"
                badges.append((label, variant, tooltip or f"Annotation tag: {label}"))

        predictor_text = _safe_value(self._row_data.get("Predictors"))
        predictors: List[str] = []
        for raw in predictor_text.replace(";", ",").split(","):
            predictor = raw.strip()
            if predictor and predictor not in predictors:
                predictors.append(predictor)
        for predictor in predictors:
            badges.append((predictor, "predictor", f"Predictor: {predictor}"))

        max_badges = 7
        if len(badges) > max_badges:
            hidden = badges[max_badges - 1 :]
            hidden_names = ", ".join(text for text, _variant, _tooltip in hidden)
            badges = badges[: max_badges - 1] + [
                (f"+{len(hidden)} more", "more", hidden_names)
            ]
        return badges

    @staticmethod
    def _make_annotation_badge(text: str, variant: str, tooltip: str) -> QtWidgets.QLabel:
        label = QtWidgets.QLabel(text)
        label.setProperty("role", "annotationChip")
        label.setProperty("variant", variant)
        label.setToolTip(tooltip or text)
        return label

    @staticmethod
    def _build_section_card(
        title: str,
        summary: str = "",
    ) -> Tuple[QtWidgets.QFrame, QtWidgets.QVBoxLayout]:
        card = QtWidgets.QFrame()
        card.setObjectName("AnnotationCard")
        card_layout = QtWidgets.QVBoxLayout(card)
        card_layout.setContentsMargins(15, 13, 15, 14)
        card_layout.setSpacing(10)

        heading_layout = QtWidgets.QHBoxLayout()
        heading_layout.setContentsMargins(0, 0, 0, 0)
        heading_layout.setSpacing(8)
        heading = QtWidgets.QLabel(title)
        heading.setProperty("role", "annotationSectionTitle")
        heading_layout.addWidget(heading)
        heading_layout.addStretch(1)
        if summary:
            meta = QtWidgets.QLabel(summary)
            meta.setProperty("role", "annotationSectionMeta")
            heading_layout.addWidget(meta)
        card_layout.addLayout(heading_layout)
        return card, card_layout

    @staticmethod
    def _build_detail_field(label_text: str, value_text: str) -> QtWidgets.QWidget:
        field = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(field)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(3)
        label = QtWidgets.QLabel(label_text)
        label.setProperty("role", "annotationFieldLabel")
        value = QtWidgets.QLabel(value_text or "-")
        value.setProperty("role", "annotationFieldValue")
        value.setWordWrap(True)
        value.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(label)
        layout.addWidget(value)
        return field

    def _build_detail_grid(self, fields: Sequence[Tuple[str, str]]) -> QtWidgets.QWidget:
        grid_widget = QtWidgets.QWidget()
        grid = QtWidgets.QGridLayout(grid_widget)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(28)
        grid.setVerticalSpacing(12)
        for index, (label, value) in enumerate(fields):
            grid.addWidget(self._build_detail_field(label, value), index // 2, index % 2)
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 1)
        return grid_widget

    def _build_summary_tab(self) -> QtWidgets.QWidget:
        page = QtWidgets.QWidget()
        page_layout = QtWidgets.QVBoxLayout(page)
        page_layout.setContentsMargins(0, 0, 0, 0)

        scroll = QtWidgets.QScrollArea()
        scroll.setObjectName("AnnotationSummaryScroll")
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        content = QtWidgets.QWidget()
        content.setObjectName("AnnotationSummaryContent")
        layout = QtWidgets.QVBoxLayout(content)
        layout.setContentsMargins(2, 10, 2, 4)
        layout.setSpacing(10)

        location_card, location_layout = self._build_section_card("Location")
        location_fields = [
            ("Sample", _safe_value(self._row_data.get("Sample")) or "Unknown sample"),
            ("Contig", _safe_value(self._row_data.get("Contig")) or "Unknown contig"),
            ("Coordinates", self._format_position()),
            ("Strand", self._format_orientation()),
        ]
        location_layout.addWidget(self._build_detail_grid(location_fields))
        layout.addWidget(location_card)

        if has_resistance_mutation(self._row_data):
            layout.addWidget(self._build_mutation_card())

        identifier_fields: List[Tuple[str, str]] = []
        for label, key in [
            ("Annotation ID", "AnnotationID"),
            ("NCBI protein", "NCBI_NRP_ID"),
            ("RefSeq protein", "Refseq_ID"),
            ("Locus", "Locus"),
            ("Gene ID", "GeneID"),
        ]:
            value = _safe_value(self._row_data.get(key))
            if value:
                identifier_fields.append((label, value))
        if identifier_fields:
            ids_card, ids_layout = self._build_section_card("Identifiers")
            ids_layout.addWidget(self._build_detail_grid(identifier_fields))
            layout.addWidget(ids_card)

        nt_seq = _safe_value(self._row_data.get("Nucleotide_Sequence"))
        aa_seq = _safe_value(self._row_data.get("Amino_Acid_Sequence"))
        if nt_seq or aa_seq:
            seq_count = int(bool(nt_seq)) + int(bool(aa_seq))
            seq_card, seq_layout = self._build_section_card(
                "Sequences",
                f"{seq_count} available",
            )
            if nt_seq:
                seq_layout.addWidget(self._build_sequence_row(nt_seq, "Nucleotide", "nt"))
            if aa_seq:
                seq_layout.addWidget(self._build_sequence_row(aa_seq, "Amino acid", "aa"))
            layout.addWidget(seq_card)

        marker_hits = self._marker_hits_for_row()
        if not marker_hits.empty:
            marker_count = len(marker_hits)
            marker_card, marker_layout = self._build_section_card(
                "Marker evidence",
                f"{marker_count} {'hit' if marker_count == 1 else 'hits'}",
            )
            table = self._build_marker_table(marker_hits)
            self._marker_table = table
            marker_layout.addWidget(table)
            layout.addWidget(marker_card)

        layout.addStretch(1)
        scroll.setWidget(content)
        page_layout.addWidget(scroll)
        return page

    def _build_mutation_card(self) -> QtWidgets.QFrame:
        """AMRFinderPlus point mutations of this CDS: the badge, the note, symbols and drug classes."""
        mutations = resistance_mutations(self._row_data)
        card, card_layout = self._build_section_card(
            "Resistance mutations",
            "AMRFinderPlus",
        )
        note_row = QtWidgets.QHBoxLayout()
        note_row.setContentsMargins(0, 0, 0, 0)
        note_row.setSpacing(8)
        icon = QtWidgets.QLabel()
        icon.setPixmap(QtGui.QIcon(str(MUTATION_BADGE_SVG)).pixmap(QtCore.QSize(16, 16), self.devicePixelRatioF()))
        icon.setAlignment(QtCore.Qt.AlignmentFlag.AlignTop)
        note_row.addWidget(icon)
        note = QtWidgets.QLabel(mutation_note(self._row_data))
        note.setWordWrap(True)
        note_row.addWidget(note, stretch=1)
        card_layout.addLayout(note_row)
        fields = [("Mutations", ", ".join(mutations) or "Tagged X")]
        classes = _safe_value(self._row_data.get("ResistanceMutationClasses")).replace(";", ", ")
        if classes:
            fields.append(("Drug classes", classes))
        card_layout.addWidget(self._build_detail_grid(fields))
        return card

    def _build_marker_table(self, marker_hits: pd.DataFrame) -> QtWidgets.QTableWidget:
        table = QtWidgets.QTableWidget()
        table.setObjectName("AnnotationMarkerTable")
        table.setColumnCount(5)
        table.setHorizontalHeaderLabels(["DB", "Method", "Coverage", "Identity", "Product"])
        table.setRowCount(len(marker_hits))
        table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        table.setAlternatingRowColors(True)
        table.setHorizontalScrollMode(QtWidgets.QAbstractItemView.ScrollMode.ScrollPerPixel)
        table.setVerticalScrollMode(QtWidgets.QAbstractItemView.ScrollMode.ScrollPerPixel)
        table.verticalHeader().setVisible(False)
        table.verticalHeader().setDefaultSectionSize(32)

        header = table.horizontalHeader()
        header.setStretchLastSection(False)
        header.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(2, QtWidgets.QHeaderView.ResizeMode.Fixed)
        header.setSectionResizeMode(3, QtWidgets.QHeaderView.ResizeMode.Fixed)
        header.setSectionResizeMode(4, QtWidgets.QHeaderView.ResizeMode.Stretch)
        table.setColumnWidth(2, 90)
        table.setColumnWidth(3, 90)

        for row_idx, marker in enumerate(marker_hits.to_dict("records")):
            db = str(marker.get("DB") or "")
            method = str(marker.get("Method") or "")
            product = str(marker.get("Product") or "")
            coverage = marker.get("Subject_coverage") or ""
            identity = marker.get("Identity") or ""
            table.setItem(row_idx, 0, QtWidgets.QTableWidgetItem(db))
            table.setItem(row_idx, 1, QtWidgets.QTableWidgetItem(method))
            product_item = QtWidgets.QTableWidgetItem(product)
            product_item.setToolTip(product)
            table.setItem(row_idx, 4, product_item)
            self._set_metric_cell(table, row_idx, 2, coverage, "#86acee")
            self._set_metric_cell(table, row_idx, 3, identity, "#86acee")

        visible_rows = min(max(len(marker_hits), 1), 6)
        table.setFixedHeight(42 + (visible_rows * 32))
        return table

    def _build_raw_tab(self) -> QtWidgets.QWidget:
        widget = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(widget)
        layout.setContentsMargins(2, 10, 2, 4)
        layout.setSpacing(8)

        title = QtWidgets.QLabel("Complete annotation record")
        title.setProperty("role", "annotationSectionTitle")
        layout.addWidget(title)
        description = QtWidgets.QLabel(
            "Normalized source fields and marker evidence for this annotation."
        )
        description.setProperty("role", "muted")
        description.setWordWrap(True)
        layout.addWidget(description)

        raw = QtWidgets.QPlainTextEdit()
        raw.setObjectName("AnnotationRawJson")
        raw.setReadOnly(True)
        raw.setLineWrapMode(QtWidgets.QPlainTextEdit.LineWrapMode.NoWrap)
        raw.setFont(QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.SystemFont.FixedFont))
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
        raw.setAccessibleName("Complete annotation record as JSON")
        self._raw_editor = raw
        layout.addWidget(raw, stretch=1)
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

    def _build_sequence_row(self, seq: str, label_text: str, unit: str) -> QtWidgets.QWidget:
        widget = QtWidgets.QFrame()
        widget.setObjectName("AnnotationSequenceRow")
        layout = QtWidgets.QVBoxLayout(widget)
        layout.setContentsMargins(10, 8, 10, 9)
        layout.setSpacing(6)

        heading_layout = QtWidgets.QHBoxLayout()
        heading_layout.setContentsMargins(0, 0, 0, 0)
        heading = QtWidgets.QLabel(label_text)
        heading.setProperty("role", "annotationSequenceTitle")
        length = QtWidgets.QLabel(f"{len(seq):,} {unit}")
        length.setProperty("role", "annotationSequenceLength")
        heading_layout.addWidget(heading)
        heading_layout.addStretch(1)
        heading_layout.addWidget(length)
        layout.addLayout(heading_layout)

        sequence_layout = QtWidgets.QHBoxLayout()
        sequence_layout.setContentsMargins(0, 0, 0, 0)
        sequence_layout.setSpacing(8)
        sequence_edit = QtWidgets.QLineEdit(seq)
        sequence_edit.setObjectName("AnnotationSequenceText")
        sequence_edit.setReadOnly(True)
        sequence_edit.setMinimumWidth(0)
        sequence_edit.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Ignored,
            QtWidgets.QSizePolicy.Policy.Fixed,
        )
        sequence_edit.setFont(
            QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.SystemFont.FixedFont)
        )
        sequence_edit.setCursorPosition(0)
        sequence_edit.setAccessibleName(f"{label_text} sequence")
        self._sequence_edits.append(sequence_edit)

        copy_btn = QtWidgets.QToolButton()
        copy_btn.setObjectName("AnnotationCopyButton")
        copy_btn.setText("Copy")
        copy_btn.setToolTip(f"Copy the complete {label_text.lower()} sequence")
        copy_btn.setAccessibleName(f"Copy {label_text.lower()} sequence")
        copy_btn.clicked.connect(lambda _checked=False, s=seq: self._copy_text(s))
        sequence_layout.addWidget(sequence_edit, stretch=1)
        sequence_layout.addWidget(copy_btn)
        layout.addLayout(sequence_layout)
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
        bar.setFixedHeight(20)
        bar.setStyleSheet(
            "QProgressBar {"
            " border: none;"
            " border-radius: 6px;"
            " background: #edf2f7;"
            " color: #172434;"
            " text-align: center;"
            " font-weight: 600;"
            "}"
            f"QProgressBar::chunk {{ background-color: {color}; border-radius: 6px; }}"
        )
        table.setCellWidget(row, col, bar)
        table.setItem(row, col, QtWidgets.QTableWidgetItem(f"{metric:.1f}"))

    # Owner window and actions

    def menu_label(self) -> str:
        """Names the sample and contig too: cards of one gene in several
        samples have the same title."""
        sample = _safe_value(self._row_data.get("Sample"))
        contig = _safe_value(self._row_data.get("Contig"))
        where = ", ".join(part for part in (sample, f"contig {contig}" if contig else "") if part)
        label = f"{self.card_title()} ({where})" if where else self.card_title()
        badge = self.owner_badge()
        return f"{label} [{badge}]" if badge else label

    def _owner_changed(self) -> None:
        self._update_action_buttons()

    def _update_action_buttons(self) -> None:
        """The actions need an owner explore window to show their results in."""
        enabled = self.can_run_commands()
        tip_off = "Open this card from an explore window to use its actions"
        from pipeline.ToolEnv import viewer_mode

        no_aligner = viewer_mode() and find_aligner() is None
        for button, tip_on in (
            (
                getattr(self, "_context_button", None),
                "Compare this gene with matching neighborhoods across samples, "
                "in a new explore window with this window's filters",
            ),
            (getattr(self, "_align_button", None), "Align this protein with matching available sequences"),
        ):
            if button is None:
                continue
            if no_aligner and button is getattr(self, "_align_button", None):
                button.setEnabled(False)
                button.setToolTip(ALIGN_NEEDS_TOOLS_TIP)
                continue
            button.setEnabled(enabled)
            button.setToolTip(tip_on if enabled else tip_off)

    def _open_alignment_dialog(self) -> None:
        """The parameter form, attached to this card (one at a time)."""
        if not self.can_run_commands():
            return  # no owner window to show the result in
        form = self._alignment_form
        if form is None:
            form = GeneCardAlignmentDialog(
                self,
                row_data=self._row_data,
                ann_df=self._ann_df,
                on_result=self._show_alignment_result,
            )
            form.setWindowModality(QtCore.Qt.WindowModality.NonModal)
            form.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose, True)
            form.destroyed.connect(self._forget_alignment_form)
            self._alignment_form = form
        form.show()
        form.raise_()
        form.activateWindow()

    def _forget_alignment_form(self, *_args: object) -> None:
        self._alignment_form = None

    def _show_alignment_result(self, alignment_text: str, status_text: str) -> None:
        gene = _safe_value(self._row_data.get("Gene")) or "gene"
        self.run_command(
            PANEL_OPEN_ALIGNMENT,
            alignment_text=alignment_text,
            status_text=status_text,
            title=f"Alignment: {gene}",
        )

    def _open_context_dialog(self) -> None:
        """Compare neighborhoods: a neighbourhood panel in a new explore
        window cloned from the owner's filters, busy (inline) until its rows
        are collected on a worker thread. Asking again for the same gene shows
        that panel (in whichever window holds it)."""
        if not self.can_run_commands():
            return  # no owner window to show the neighborhoods in
        gene = _safe_value(self._row_data.get("Gene")) or "gene"
        self.run_command(
            WINDOW_OPEN_NEIGHBOURHOOD,
            # Load the largest supported neighborhood once; the panel can then
            # switch between ±5/±10/±15/±30 genes without another scan.
            loader=self._collect_context_rows,
            title=f"Gene neighborhoods: {gene}",
            source="gene_card",
            key=("gene", self.info_key()) if self.info_key() is not None else None,
            reference_sample=_safe_value(self._row_data.get("Sample")),
            flank_sizes=(5, 10, 15, 30),
            initial_flank=15,
        )

    def _collect_context_rows(self) -> List[Dict[str, Any]]:
        """The Compare neighborhoods loader (runs on a worker thread)."""
        with self._context_rows_lock:
            return self._build_context_collection_rows(window_genes=30)

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
            plasmid_call = _safe_value(center_row.get("MobSuiteIsPlasmid")).lower()
            if plasmid_call in {"<na>", "none", "null", "-"}:
                plasmid_call = ""
            is_plasmid = plasmid_call in PLASMID_TRUE_VALUES
            # Without a MOB-suite call the molecule type is unknown, not
            # chromosome.
            molecule = "plasmid" if is_plasmid else "chromosome" if plasmid_call else "unknown"
            contig_label = contig if contig.lower().startswith("contig") else f"contig {contig}"
            row_payload = {
                "title": f"{sample}, {contig_label}",
                "sample": sample,
                "contig": contig,
                "molecule": molecule,
                "is_plasmid": is_plasmid,
                "is_source": is_source,
                "sample_meta": sample_meta,
                "meta_tooltip": self._sample_metadata_tooltip(sample, contig, sample_meta),
                "center_bp": center_bp,
                "genes": genes,
                # Genes the contig holds on each side of the focal gene, so the
                # comparison can tell 'the contig ends here' from 'nothing
                # similar here' at any flank.
                "contig_genes_before": int(hit_idx),
                "contig_genes_after": int(len(contig_df) - 1 - hit_idx),
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
        if "A" in value or MUTATION_TAG in value:
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
