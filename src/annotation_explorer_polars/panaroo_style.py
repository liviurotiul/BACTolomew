from __future__ import annotations

from typing import Dict, Mapping

from PyQt6 import QtCore, QtGui, QtWidgets

from .panaroo_common import (
    _FASTGEAR_LINEAGE_PALETTE,
    _METADATA_STRIP_PALETTE,
    _natural_sort_key,
)


_ANALYSIS_PANEL_STYLESHEET = """
#NeoAnalysisPanel {
  background: #ffffff;
  border: 1px solid #d8e2ed;
  border-radius: 8px;
}
#NeoAnalysisPanel QStackedWidget {
  background: transparent;
  border: 0;
}
QGroupBox#NeoAnalysisGroup {
  background: #ffffff;
  border: 1px solid #d8e2ed;
  border-radius: 8px;
  margin-top: 14px;
  color: #334155;
  font-weight: 700;
}
QGroupBox#NeoAnalysisGroup::title {
  subcontrol-origin: margin;
  left: 10px;
  padding: 0 4px;
}
"""


def _log_view_event(message: str) -> None:
    print(f"[view] {message}")


def _apply_analysis_panel_style(widget: QtWidgets.QWidget) -> None:
    widget.setObjectName("NeoAnalysisPanel")
    widget.setAttribute(QtCore.Qt.WidgetAttribute.WA_StyledBackground, True)
    widget.setStyleSheet(_ANALYSIS_PANEL_STYLESHEET)
    widget.setSizePolicy(
        QtWidgets.QSizePolicy.Policy.Expanding,
        QtWidgets.QSizePolicy.Policy.Expanding,
    )


def _apply_analysis_group_style(group: QtWidgets.QGroupBox) -> None:
    group.setObjectName("NeoAnalysisGroup")
    group.setAttribute(QtCore.Qt.WidgetAttribute.WA_StyledBackground, True)
    group.setStyleSheet(_ANALYSIS_PANEL_STYLESHEET)


def _fastgear_lineage_colors(lineage_by_sample: Mapping[str, str]) -> Dict[str, QtGui.QColor]:
    unique = sorted(
        {str(lineage).strip() for lineage in lineage_by_sample.values() if str(lineage).strip()}
    )
    colors: Dict[str, QtGui.QColor] = {}
    for idx, lineage in enumerate(unique):
        colors[lineage] = QtGui.QColor(_FASTGEAR_LINEAGE_PALETTE[idx % len(_FASTGEAR_LINEAGE_PALETTE)])
    return colors


def _fastgear_cluster_colors(cluster_by_sample: Mapping[str, str]) -> Dict[str, QtGui.QColor]:
    unique = sorted(
        {str(cluster).strip() for cluster in cluster_by_sample.values() if str(cluster).strip()},
        key=_natural_sort_key,
    )
    colors: Dict[str, QtGui.QColor] = {}
    for idx, cluster in enumerate(unique):
        colors[cluster] = QtGui.QColor(_METADATA_STRIP_PALETTE[idx % len(_METADATA_STRIP_PALETTE)])
    return colors
