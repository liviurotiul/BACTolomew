from __future__ import annotations

from typing import Dict, Mapping

from PyQt6 import QtCore, QtGui, QtWidgets

from .panaroo_common import (
    _FASTGEAR_LINEAGE_PALETTE,
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


def _fastgear_lineage_colors(lineage_by_sample: Mapping[str, str]) -> Dict[str, QtGui.QColor]:
    unique = sorted(
        {str(lineage).strip() for lineage in lineage_by_sample.values() if str(lineage).strip()},
        key=_natural_sort_key,
    )
    colors: Dict[str, QtGui.QColor] = {}
    for idx, lineage in enumerate(unique):
        numeric_label = lineage
        if numeric_label.lower().startswith("clade"):
            numeric_label = numeric_label[len("clade"):].strip(" _-:")
        palette_idx = int(numeric_label) - 1 if numeric_label.isdigit() and int(numeric_label) > 0 else idx
        colors[lineage] = QtGui.QColor(
            _FASTGEAR_LINEAGE_PALETTE[palette_idx % len(_FASTGEAR_LINEAGE_PALETTE)]
        )
    return colors
