from __future__ import annotations

from typing import Optional

from PyQt6 import QtCore

from .annotation_table import AnnotationTableModel


class HoverSyncMixin:
    def _apply_pan_visual_highlight(self) -> None:
        active_sample = self._hover_active_sample
        active_source = self._hover_active_source

        self._apply_tree_sample_hover(active_sample)
        if self.pan_neo_tab is not None and self.pan_neo_tab is not active_source:
            self.pan_neo_tab.highlight_sample(active_sample)
        if self._phandango_pan_heatmap is not None and self._phandango_pan_heatmap.isVisible() and self._phandango_pan_heatmap is not active_source:
            self._phandango_pan_heatmap.highlight_sample(active_sample)
        if self._cgmlst_force_plot_panel is not None and self._cgmlst_force_plot_panel.isVisible() and self._cgmlst_force_plot_panel is not active_source:
            self._cgmlst_force_plot_panel.highlight_sample(active_sample)
        if self.pan_similarity_tab is not None and self.pan_similarity_tab.isVisible():
            self.pan_similarity_tab.highlight_sample(active_sample)
        if self.pan_structvar_tab is not None and self.pan_structvar_tab.isVisible():
            self.pan_structvar_tab.highlight_sample(active_sample)
        for geo_view in self._all_geo_views():
            if geo_view.isVisible():
                geo_view.highlight_sample(active_sample)

    def eventFilter(self, obj: QtCore.QObject, event: QtCore.QEvent) -> bool:
        if obj is self.tree.viewport():
            if event.type() == QtCore.QEvent.Type.MouseMove:
                pos = event.position().toPoint() if hasattr(event, "position") else event.pos()
                item = self.tree.itemAt(pos)
                sample = self._sample_from_tree_item(item)
                if sample != self._hover_active_sample:
                    if sample:
                        self._set_pan_highlight(sample)
                    else:
                        self._clear_pan_highlight()
            elif event.type() == QtCore.QEvent.Type.Leave:
                self._clear_pan_highlight()
        elif obj is self.table.viewport():
            if event.type() == QtCore.QEvent.Type.MouseMove:
                pos = event.position().toPoint() if hasattr(event, "position") else event.pos()
                index = self.table.indexAt(pos)
                sample = self._sample_from_table_index(index)
                if sample != self._hover_active_sample:
                    if sample:
                        self._set_pan_highlight(sample)
                    else:
                        self._clear_pan_highlight()
            elif event.type() == QtCore.QEvent.Type.Leave:
                self._clear_pan_highlight()
        return super().eventFilter(obj, event)

    def _sample_from_table_index(self, index: QtCore.QModelIndex) -> Optional[str]:
        if not index.isValid():
            return None
        model = self.table.model()
        if not isinstance(model, AnnotationTableModel):
            return None
        row_data = model.row_data(index.row())
        if not row_data:
            return None
        raw = row_data.get("Sample", self.current_sample)
        sample = str(raw).strip() if raw is not None else ""
        return sample or None

    def _on_pan_sample_hovered(self, sample: object) -> None:
        value = str(sample).strip() if sample is not None else ""
        source = self.sender()
        if value:
            self._set_pan_highlight(value, source=source)
        else:
            self._clear_pan_highlight(source=source)

    def _set_pan_highlight(self, sample: Optional[str], source: Optional[QtCore.QObject] = None) -> None:
        if self._hover_active_sample == sample and self._hover_active_source is source:
            return
        self._hover_active_sample = sample
        self._hover_active_source = source
        self._apply_pan_visual_highlight()

    def _clear_pan_highlight(self, source: Optional[QtCore.QObject] = None) -> None:
        self._hover_active_sample = None
        self._hover_active_source = source
        self._apply_pan_visual_highlight()
