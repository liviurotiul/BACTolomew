from __future__ import annotations

from typing import List, Sequence, Tuple

from PyQt6 import QtCore, QtWidgets


class MetadataValueDialog(QtWidgets.QDialog):
    def __init__(
        self,
        parent: QtWidgets.QWidget,
        column: str,
        values: Sequence[Tuple[str, str]],
        selected: Sequence[str],
    ) -> None:
        super().__init__(parent)
        self._selected = set(selected)
        self.setWindowTitle(f"Filter: {column}")
        self.resize(420, 360)
        self._build_ui(values)

    def _build_ui(self, values: Sequence[Tuple[str, str]]) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        label = QtWidgets.QLabel("Select values to include:")
        layout.addWidget(label)

        self.list_widget = QtWidgets.QListWidget()
        for display, token in values:
            item = QtWidgets.QListWidgetItem(display)
            item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            state = QtCore.Qt.CheckState.Checked if token in self._selected else QtCore.Qt.CheckState.Unchecked
            item.setCheckState(state)
            item.setData(QtCore.Qt.ItemDataRole.UserRole, token)
            self.list_widget.addItem(item)
        layout.addWidget(self.list_widget)

        btns = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel
        )
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)

    def selected_values(self) -> List[str]:
        selected: List[str] = []
        for i in range(self.list_widget.count()):
            item = self.list_widget.item(i)
            if item.checkState() == QtCore.Qt.CheckState.Checked:
                token = item.data(QtCore.Qt.ItemDataRole.UserRole)
                selected.append(str(token))
        return selected

