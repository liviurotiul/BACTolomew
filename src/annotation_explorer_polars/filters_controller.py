from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

import polars as pl
from PyQt6 import QtCore, QtGui, QtWidgets

from .constants import MISSING_FILTER_TOKEN
from .metadata_filter_dialog import MetadataValueDialog
from .model.filter_state import AppliedResult, FilterState, SampleSubset
from .settings_column import ElidedHint
from .theme import COLORS
from .shell.layout_store import HINT_SAMPLE_SEARCH_D3, hint_was_shown, mark_hint_shown

MISSING_FILTER_LABEL = "(missing)"
# A filter row's values keep this much room in the narrowest rail; a long
# column name elides first.
FILTER_SUMMARY_MIN_WIDTH = 56
SAMPLE_SEARCH_HINT_TEXT = (
    "The sample search filters samples only. To find genes, use Gene/annotation contains or Product contains."
)

TAG_OPTIONS: List[Tuple[str, str]] = [
    ("A", "AMR"),
    ("B", "BMRG"),
    ("V", "Virulence"),
    ("M", "Mobile"),
    ("P", "Phage"),
    ("I", "Compound overlap"),
    ("X", "Resistance mutation"),
]


class FiltersControllerMixin:
    def _normalize_sample_subset_values(self, samples: Sequence[str]) -> List[str]:
        seen: set[str] = set()
        ordered: List[str] = []
        for raw in samples:
            sample = str(raw).strip()
            if not sample or sample in seen:
                continue
            seen.add(sample)
            ordered.append(sample)
        return ordered

    def _set_sample_subset_filter(
        self,
        *,
        label: str,
        samples: Sequence[str],
        kind: Optional[str] = None,
        origin: str = "rail",
    ) -> None:
        """Replace the one sample-subset filter (the rail's chip) through the
        context (plan 4.4 ``set_subset``); every panel follows."""
        normalized = self._normalize_sample_subset_values(samples)
        if not normalized:
            return
        subset = SampleSubset(tuple(normalized), label=str(label).strip(), kind=str(kind or "").strip())
        self.ctx.set_subset(subset, origin=origin)

    def _clear_sample_subset_filter(self) -> None:
        if self.ctx.draft.subset is None:
            return
        self.ctx.set_subset(None, origin="rail")

    def _show_add_filter_menu(self) -> None:
        menu = QtWidgets.QMenu(self)
        for col in sorted(self.run.meta_columns):
            label = f"{col} (filtered)" if col in self.metadata_filters else col
            action = menu.addAction(label)
            action.triggered.connect(lambda _checked=False, c=col: self._edit_metadata_filter(c))
        anchor = self.add_filter_btn.mapToGlobal(QtCore.QPoint(0, self.add_filter_btn.height()))
        menu.exec(anchor)


    def _init_predictor_model(self) -> None:
        model = QtGui.QStandardItemModel(self.predictor_box)
        clear_item = QtGui.QStandardItem("Any (clear)")
        clear_item.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled)
        clear_item.setData("__clear__", QtCore.Qt.ItemDataRole.UserRole)
        model.appendRow(clear_item)
        for pred in self.run.predictor_options:
            item = QtGui.QStandardItem(pred)
            item.setFlags(
                QtCore.Qt.ItemFlag.ItemIsEnabled
                | QtCore.Qt.ItemFlag.ItemIsUserCheckable
            )
            item.setData(QtCore.Qt.CheckState.Unchecked, QtCore.Qt.ItemDataRole.CheckStateRole)
            item.setData(pred, QtCore.Qt.ItemDataRole.UserRole)
            model.appendRow(item)
        model.itemChanged.connect(self._on_predictor_model_changed)
        self._predictor_model = model
        self.predictor_box.setModel(model)
        self._update_predictor_box_text()


    def _selected_predictors(self) -> List[str]:
        preds: List[str] = []
        if not self._predictor_model:
            return preds
        for row in range(self._predictor_model.rowCount()):
            item = self._predictor_model.item(row)
            if not item:
                continue
            token = item.data(QtCore.Qt.ItemDataRole.UserRole)
            if not token or token == "__clear__":
                continue
            if item.checkState() == QtCore.Qt.CheckState.Checked:
                preds.append(str(token))
        return preds


    def _update_predictor_box_text(self) -> None:
        selected = self._selected_predictors()
        if selected:
            txt = f"Predictor ({','.join(sorted(selected))})"
        else:
            txt = "Predictor (any)"
        if self.predictor_box.isEditable() and self.predictor_box.lineEdit() is not None:
            line_edit = self.predictor_box.lineEdit()
            line_edit.setText(txt)
            line_edit.setCursorPosition(0)
            line_edit.deselect()
        else:
            self.predictor_box.setCurrentText(txt)


    def _toggle_predictor_item(self, index: QtCore.QModelIndex) -> None:
        if not self._predictor_model:
            return
        item = self._predictor_model.itemFromIndex(index)
        if not item:
            return
        token = item.data(QtCore.Qt.ItemDataRole.UserRole)
        if token == "__clear__":
            for row in range(self._predictor_model.rowCount()):
                it = self._predictor_model.item(row)
                if not it:
                    continue
                if it.data(QtCore.Qt.ItemDataRole.UserRole) == "__clear__":
                    continue
                it.setCheckState(QtCore.Qt.CheckState.Unchecked)
            self._update_predictor_box_text()
            self._apply_rail_filters()
            return
        state = item.checkState()
        item.setCheckState(
            QtCore.Qt.CheckState.Unchecked
            if state == QtCore.Qt.CheckState.Checked
            else QtCore.Qt.CheckState.Checked
        )


    def _on_predictor_model_changed(self, _item: QtGui.QStandardItem) -> None:
        self._update_predictor_box_text()
        self._apply_rail_filters()


    def _init_tag_model(self) -> None:
        model = QtGui.QStandardItemModel(self.tag_box)
        clear_item = QtGui.QStandardItem("Any (clear)")
        clear_item.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled)
        clear_item.setData("__clear__", QtCore.Qt.ItemDataRole.UserRole)
        model.appendRow(clear_item)
        for tag, label in TAG_OPTIONS:
            # A tag whose databases were not installed for this run says so, so that
            # no hits is not read as nothing found; one found only by Bakta product
            # names says that. Its tooltip names its sources (a run with gene evidence).
            databases = self.run.tag_databases.get(tag)
            note = databases.note() if databases is not None else ""
            item = QtGui.QStandardItem(f"{tag}: {label}" + (f" ({note})" if note else ""))
            item.setFlags(
                QtCore.Qt.ItemFlag.ItemIsEnabled
                | QtCore.Qt.ItemFlag.ItemIsUserCheckable
            )
            item.setData(QtCore.Qt.CheckState.Unchecked, QtCore.Qt.ItemDataRole.CheckStateRole)
            item.setData(tag, QtCore.Qt.ItemDataRole.UserRole)
            if databases is not None and databases.tooltip():
                item.setToolTip(databases.tooltip())
            if note:
                muted = COLORS["text_disabled" if note == "not installed" else "text_muted"]
                item.setForeground(QtGui.QBrush(QtGui.QColor(muted)))
            model.appendRow(item)
        model.itemChanged.connect(self._on_tag_model_changed)
        self._tag_model = model
        self.tag_box.setModel(model)
        self._update_tag_box_text()


    def _selected_tags(self) -> List[str]:
        tags: List[str] = []
        if not self._tag_model:
            return tags
        for row in range(self._tag_model.rowCount()):
            item = self._tag_model.item(row)
            if not item:
                continue
            tag = item.data(QtCore.Qt.ItemDataRole.UserRole)
            if not tag or tag == "__clear__":
                continue
            if item.checkState() == QtCore.Qt.CheckState.Checked:
                tags.append(str(tag))
        return tags


    def _update_tag_box_text(self) -> None:
        selected = self._selected_tags()
        if selected:
            txt = f"Tags ({','.join(sorted(selected))})"
        else:
            txt = "Tags (any)"
        if self.tag_box.isEditable() and self.tag_box.lineEdit() is not None:
            line_edit = self.tag_box.lineEdit()
            line_edit.setText(txt)
            line_edit.setCursorPosition(0)
            line_edit.deselect()
        else:
            self.tag_box.setCurrentText(txt)


    def _toggle_tag_item(self, index: QtCore.QModelIndex) -> None:
        if not self._tag_model:
            return
        item = self._tag_model.itemFromIndex(index)
        if not item:
            return
        tag = item.data(QtCore.Qt.ItemDataRole.UserRole)
        if tag == "__clear__":
            for row in range(self._tag_model.rowCount()):
                it = self._tag_model.item(row)
                if not it:
                    continue
                if it.data(QtCore.Qt.ItemDataRole.UserRole) == "__clear__":
                    continue
                it.setCheckState(QtCore.Qt.CheckState.Unchecked)
            self._update_tag_box_text()
            self._apply_rail_filters()
            return
        state = item.checkState()
        item.setCheckState(
            QtCore.Qt.CheckState.Unchecked
            if state == QtCore.Qt.CheckState.Checked
            else QtCore.Qt.CheckState.Checked
        )


    def _on_tag_model_changed(self, _item: QtGui.QStandardItem) -> None:
        self._update_tag_box_text()
        self._apply_rail_filters()


    def _edit_metadata_filter(self, column: str) -> None:
        values = self._metadata_value_options(column)
        current = set(self.metadata_filters.get(column, []))
        dlg = MetadataValueDialog(self, column, values, current)
        if dlg.exec() != QtWidgets.QDialog.DialogCode.Accepted:
            return
        selected = dlg.selected_values()
        if not selected:
            if column in self.metadata_filters:
                del self.metadata_filters[column]
        else:
            self.metadata_filters[column] = selected
        self._render_metadata_filters()
        self._apply_rail_filters()


    def _metadata_value_options(self, column: str) -> List[Tuple[str, str]]:
        if column not in self.run.meta_pl.columns:
            return []
        values_lf = self.run.meta_pl.lazy().select(pl.col(column).cast(pl.Utf8).str.strip_chars().alias("v"))
        has_missing = (
            values_lf
            .filter(pl.col("v").is_null() | (pl.col("v") == ""))
            .select(pl.len().alias("n"))
            .collect()
            .item(0, "n")
            > 0
        )
        values = (
            values_lf
            .filter(pl.col("v").is_not_null() & (pl.col("v") != ""))
            .unique()
            .sort("v")
            .collect()
            .get_column("v")
            .to_list()
        )
        values = sorted({str(v) for v in values if str(v).strip()})
        options: List[Tuple[str, str]] = [(v, v) for v in values]
        if has_missing:
            options.insert(0, (MISSING_FILTER_LABEL, MISSING_FILTER_TOKEN))
        return options


    def _render_metadata_filters(self) -> None:
        while self.filters_layout.count():
            item = self.filters_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
        window_subset = getattr(self, "_window_subset", None)
        if window_subset is not None and window_subset.samples:
            # The locked subset of a window opened from a clade, cluster or
            # plasmid (plan 3.10): shown first, never removable.
            # The label goes on its own wrapping line under "Window ...
            # Locked", so a long subset name never hides the count.
            row = QtWidgets.QWidget()
            row.setObjectName("WindowSubsetChip")
            row_layout = QtWidgets.QGridLayout(row)
            row_layout.setContentsMargins(0, 0, 0, 0)
            row_layout.setVerticalSpacing(2)
            label = QtWidgets.QLabel("Window")
            count = len(window_subset.samples)
            summary_text = window_subset.label or f"{count} samples"
            summary = QtWidgets.QLabel(f"{summary_text} | n={count}")
            summary.setObjectName("WindowSubsetSummary")
            summary.setWordWrap(True)
            summary.setSizePolicy(QtWidgets.QSizePolicy.Policy.Ignored, QtWidgets.QSizePolicy.Policy.Preferred)
            summary.setToolTip(", ".join(window_subset.samples))
            locked = QtWidgets.QLabel("Locked")
            locked.setProperty("role", "muted")
            locked.setToolTip(
                "This window shows only these samples. Open a new explore window "
                "(Window > New Explore Window) to see all samples."
            )
            row_layout.addWidget(label, 0, 0)
            row_layout.addWidget(locked, 0, 1, QtCore.Qt.AlignmentFlag.AlignRight)
            row_layout.addWidget(summary, 1, 0, 1, 2)
            row_layout.setColumnStretch(0, 1)
            self.filters_layout.addWidget(row)
        subset = self.ctx.draft.subset
        subset_samples = list(subset.samples) if subset is not None else []
        subset_label = str(subset.label or "").strip() if subset is not None else ""
        if subset_samples:
            row = QtWidgets.QWidget()
            row_layout = QtWidgets.QHBoxLayout(row)
            row_layout.setContentsMargins(0, 0, 0, 0)
            label = QtWidgets.QLabel("Sample set")
            summary_text = subset_label or f"{len(subset_samples)} samples"
            summary_line = f"{summary_text} | n={len(subset_samples)}"
            summary = ElidedHint(summary_line)
            # The tooltip keeps what a narrow rail elides, then the samples.
            summary.setToolTip(f"{summary_line}\n{', '.join(subset_samples)}")
            remove_btn = QtWidgets.QPushButton("Remove")
            remove_btn.clicked.connect(self._clear_sample_subset_filter)
            row_layout.addWidget(label)
            row_layout.addWidget(summary, stretch=1)
            row_layout.addWidget(remove_btn)
            self.filters_layout.addWidget(row)
        for column, values in self.metadata_filters.items():
            row = QtWidgets.QWidget()
            row_layout = QtWidgets.QHBoxLayout(row)
            row_layout.setContentsMargins(0, 0, 0, 0)
            # Both elide in a narrow rail, so a long column name leaves room
            # for its values (the tooltips show them whole).
            label = ElidedHint(column)
            label.setToolTip(column)
            summary = ElidedHint(self._format_filter_values(values))
            summary.setMinimumWidth(FILTER_SUMMARY_MIN_WIDTH)
            summary.setToolTip(self._format_filter_values(values, full=True))
            remove_btn = QtWidgets.QPushButton("Remove")
            remove_btn.clicked.connect(lambda _checked=False, c=column: self._remove_metadata_filter(c))
            row_layout.addWidget(label)
            row_layout.addWidget(summary, stretch=1)
            row_layout.addWidget(remove_btn)
            self.filters_layout.addWidget(row)
        # No empty gap under "Add filter" in the Sample card.
        self.filters_container.setVisible(self.filters_layout.count() > 0)
        self.filters_layout.addStretch(1)


    def _remove_metadata_filter(self, column: str) -> None:
        if column in self.metadata_filters:
            del self.metadata_filters[column]
            self._render_metadata_filters()
            self._apply_rail_filters()


    def _format_filter_values(self, values: Sequence[str], full: bool = False) -> str:
        display: List[str] = []
        for v in values:
            if v == MISSING_FILTER_TOKEN:
                display.append(MISSING_FILTER_LABEL)
            else:
                display.append(str(v))
        if full or len(display) <= 5:
            return ", ".join(display)
        return ", ".join(display[:5]) + f", … (+{len(display) - 5} more)"


    # The rail writes the window's draft filters (plan 4.4); text fields
    # debounce, the other controls apply immediately.

    def _schedule_refresh(self, *_args: object) -> None:
        """A rail text field changed: apply after the debounce delay."""
        self._rail_filters_changed(debounce=True)

    def _apply_rail_filters(self, *_args: object) -> None:
        """A rail chip, checkbox or list changed: apply now."""
        self._rail_filters_changed(debounce=False)

    def _rail_filters_changed(self, *, debounce: bool = False) -> None:
        self.ctx.update_draft(FilterState.from_widgets(self), debounce=debounce)

    def _applied_filters(self) -> FilterState:
        """The filters of the result on screen."""
        return self.ctx.applied


    def _on_annotations_view_mode_changed(self, mode: str) -> None:
        """The paired Annotations panel's mode changed: redraw it and its track."""
        if mode == "Contig":
            self._igv_dirty = True
        self._apply_view_mode()
        self._update_active_panel_status()


    def _has_active_annotation_filter(self) -> bool:
        """An annotation filter is applied (gene, product, molecule, predictors, tags, hide hypothetical)."""
        return self._applied_filters().has_annotation_filter()


    def _should_show_global(self) -> bool:
        """All matches lists rows only with an applied annotation filter (D3:
        the sample search no longer counts as an annotation search)."""
        return self._has_active_annotation_filter()


    def _update_sample_search_hint(self, filters: FilterState, result: AppliedResult) -> None:
        """One-time hint (D3, plan 4.3): the first time the sample search matches
        no sample while the same text is in gene or product fields, say where
        genes are searched. Shown once per installation, never again."""
        label = getattr(self, "sample_search_hint", None)
        if label is None:
            return
        shown_for = getattr(self, "_sample_search_hint_query", None)
        if shown_for is not None and filters.sample_query != shown_for:
            label.setVisible(False)
            self._sample_search_hint_query = None
        if not result.sample_query_matches_annotations or hint_was_shown(HINT_SAMPLE_SEARCH_D3):
            return
        mark_hint_shown(HINT_SAMPLE_SEARCH_D3)
        self._sample_search_hint_query = filters.sample_query
        label.setText(SAMPLE_SEARCH_HINT_TEXT)
        label.setVisible(True)
