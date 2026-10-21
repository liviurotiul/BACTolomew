from __future__ import annotations

import atexit
from pathlib import Path
from string import Template
from types import MappingProxyType

from PyQt6 import QtCore, QtGui, QtWidgets, sip


APP_ORGANIZATION = "BACTolomew"
APP_DOMAIN = "bactolomew.app"
APP_NAME = "BACTolomew"
APP_DISPLAY_NAME = "BACTolomew Annotation Explorer"


# Semantic UI colors. Plot, taxonomy, resistance, and other scientific colors
# deliberately remain with their owning views; these values describe app chrome.
COLORS = MappingProxyType(
    {
        "window": "#eef3f8",
        "dialog": "#f4f7fb",
        "surface": "#ffffff",
        "surface_subtle": "#f8fbfd",
        "surface_hover": "#f5f9ff",
        # Table and tree row hover; also the hovered row in the neighborhood plot.
        "row_hover": "#eef4fc",
        "surface_pressed": "#e4eeff",
        "border": "#d7e0ea",
        "border_strong": "#ccd6e2",
        "border_soft": "#e5ebf2",
        "text": "#1f2c3b",
        "text_strong": "#172434",
        "text_muted": "#607183",
        "text_disabled": "#93a1b2",
        "primary": "#2f6fed",
        "primary_hover": "#295fd0",
        "primary_soft": "#e3edff",
        "primary_soft_hover": "#ebf2ff",
        "focus": "#4d86f7",
        "focus_soft": "#c8dafe",
        "danger": "#c0322b",
        "success": "#0f7a37",
        "warning": "#a86612",
        # Soft fills and borders of status pills and boxes (run cards, the
        # Run progress window's error box); the same shades as the chips.
        "danger_soft": "#fff0ed",
        "danger_border": "#f2c9c2",
        "success_soft": "#edf8f2",
        "success_border": "#c8e7d4",
        "warning_soft": "#fff7e8",
        "warning_border": "#efd9ad",
        # Primary text on a primary_soft fill, and the border of a hovered card.
        "primary_ink": "#1e4f9f",
        "focus_hover": "#9ab8f8",
        # A disabled colored (primary or danger) button.
        "disabled_fill": "#f5f7fa",
        "disabled_border": "#dde4ec",
        # Tree chrome: neutral branches, and the grey of
        # genomes outside the window's filter in a greyed tree.
        "tree_branch": "#64748b",
        "tree_inactive": "#94a3b8",
        # Hover emphasis shared by every panel (the Map's hovered countries,
        # the Phylogeny's hovered path and tip): an amber no cluster or
        # window color uses. ``hover_ink`` is the readable text shade.
        "hover": "#f59e0b",
        "hover_ink": "#b45309",
    }
)

# One color per explore window (by window number), shown as a stripe along the
# top of the window and of its info cards while several windows are open
# (plan 3.10). Distinct hues of similar weight; the first is the app's primary.
WINDOW_COLORS = (
    "#2f6fed",
    "#e8590c",
    "#0f9d58",
    "#9c36b5",
    "#d6336c",
    "#0b7285",
    "#b8860b",
    "#5c6b7a",
)


def window_color(number: int) -> str:
    """The color of explore window ``number`` (1, 2, ...)."""
    index = max(0, int(number or 1) - 1)
    return WINDOW_COLORS[index % len(WINDOW_COLORS)]


# Genomes outside the window's sample set, in the views that keep them
# greyed instead of hiding them (the Phylogeny tree's branches, labels and
# collapsed clades, the linked Gene presence/absence rows): drawn in the
# ``tree_inactive`` grey at this opacity, so the sample set stands out.
OUT_OF_SET_OPACITY = 0.25


def out_of_set_color(color: object = None) -> QtGui.QColor:
    """``color`` (default: the ``tree_inactive`` grey) as an out of set
    genome is drawn: at OUT_OF_SET_OPACITY of its own alpha."""
    out = QtGui.QColor(COLORS["tree_inactive"] if color is None else color)
    out.setAlphaF(out.alphaF() * OUT_OF_SET_OPACITY)
    return out


def out_of_set_rgb(rgb: tuple, background: tuple = (255, 255, 255)) -> tuple:
    """An opaque ``rgb`` at OUT_OF_SET_OPACITY over ``background`` (for the
    views drawn as images, such as the heatmap's cells)."""
    share = OUT_OF_SET_OPACITY
    return tuple(int(round(b + (c - b) * share)) for c, b in zip(rgb, background))


# Drop-down chevrons and the tick of checkable list rows. The SVGs carry the
# text_muted, text_disabled and surface colors above; keep them in step.
_ASSETS = Path(__file__).resolve().parent / "assets"
_IMAGES = MappingProxyType(
    {
        "chevron_down": (_ASSETS / "chevron_down.svg").as_posix(),
        "chevron_down_disabled": (_ASSETS / "chevron_down_disabled.svg").as_posix(),
        "check_tick": (_ASSETS / "check_tick.svg").as_posix(),
    }
)

# Toolbar icons (text_muted strokes), by name.
_ICONS = MappingProxyType(
    {
        "zoom_in": _ASSETS / "zoom_in.svg",
        "zoom_out": _ASSETS / "zoom_out.svg",
        "zoom_fit": _ASSETS / "zoom_fit.svg",
        "rows_readable": _ASSETS / "rows_readable.svg",
        "search": _ASSETS / "search.svg",
    }
)


def icon(name: str) -> QtGui.QIcon:
    """A toolbar icon of the theme ("zoom_in", "zoom_out", "zoom_fit",
    "rows_readable", "search")."""
    return QtGui.QIcon(str(_ICONS[name]))


APP_STYLESHEET = Template(
    """
    QWidget {
        color: $text;
        selection-background-color: $primary_soft;
        selection-color: $text_strong;
    }
    QMainWindow {
        background: $window;
    }
    QDialog {
        background: $dialog;
    }
    QWidget#ContextTitleBlock {
        border-bottom: 1px solid $border;
    }
    QLabel[role="contextTitle"] {
        color: $text_strong;
        font-size: 17px;
        font-weight: 700;
    }
    QLabel[role="contextHint"] {
        color: $text_muted;
        font-size: 11px;
    }
    QFrame#ContextFooter {
        background: transparent;
        border: none;
        border-top: 1px solid $border;
    }
    QFrame#ContextPinChip {
        background: $primary;
        border: 1px solid $primary;
        border-radius: 9px;
    }
    QFrame#ContextPinChip QLabel {
        color: $surface;
        font-size: 11px;
        font-weight: 600;
    }
    QFrame#ContextPinChip QToolButton {
        min-height: 12px;
        min-width: 14px;
        padding: 0px 2px;
        border: none;
        color: $surface;
        font-weight: 700;
    }
    QFrame#ContextPinChip QToolButton:hover {
        background: $primary_hover;
    }
    QFrame#contextCanvasToolbar, QFrame#TreeToolbar {
        background: $surface_subtle;
        border: 1px solid $border;
        border-radius: 8px;
    }
    QFrame#ContextInspector {
        background: $surface;
        border: 1px solid $border;
        border-radius: 10px;
    }
    QFrame#ContextInspector QTextBrowser {
        background: $surface;
        border: none;
        padding: 4px;
    }
    QFrame#ContextFilterChip {
        background: $primary_soft_hover;
        border: 1px solid #c8dafe;
        border-radius: 8px;
    }
    QDialog#AnnotationDetailDialog {
        background: $window;
    }
    QDialog#AnnotationDetailDialog QFrame#AnnotationHero,
    QDialog#AnnotationDetailDialog QFrame#AnnotationCard {
        background: $surface;
        border: 1px solid $border;
        border-radius: 12px;
    }
    QLabel[role="annotationEyebrow"] {
        color: $primary;
        font-size: 10px;
        font-weight: 700;
        letter-spacing: 1px;
    }
    QLabel[role="annotationTitle"] {
        color: $text_strong;
        font-size: 21px;
        font-weight: 700;
    }
    QLabel[role="annotationProduct"] {
        color: #33465a;
        font-size: 13px;
    }
    QLabel[role="annotationSource"] {
        color: $text_muted;
        font-size: 11px;
    }
    QLabel[role="annotationChip"] {
        background: #edf3fa;
        border: 1px solid #d6e1ec;
        border-radius: 9px;
        color: #40556b;
        font-size: 10px;
        font-weight: 600;
        padding: 3px 8px;
    }
    QLabel[role="annotationChip"][variant="amr"] {
        background: #fff0ed;
        border-color: #f2c9c2;
        color: #963b34;
    }
    QLabel[role="annotationChip"][variant="bmrg"] {
        background: #edf8f2;
        border-color: #c8e7d4;
        color: #2f7350;
    }
    QLabel[role="annotationChip"][variant="virulence"] {
        background: #f6effb;
        border-color: #dfcdec;
        color: #76508c;
    }
    QLabel[role="annotationChip"][variant="primary"] {
        background: $primary_soft;
        border-color: $focus_soft;
        color: #1e4f9f;
    }
    QLabel[role="annotationChip"][variant="mobile"],
    QLabel[role="annotationChip"][variant="overlap"] {
        background: #fff7e8;
        border-color: #efd9ad;
        color: #826027;
    }
    QLabel[role="annotationChip"][variant="phage"] {
        background: #edf5ff;
        border-color: #cbdff5;
        color: #3c6590;
    }
    QLabel[role="annotationSectionTitle"] {
        color: $text_strong;
        font-size: 13px;
        font-weight: 700;
    }
    QLabel[role="annotationSectionMeta"],
    QLabel[role="annotationActionHint"] {
        color: $text_muted;
        font-size: 11px;
    }
    QLabel[role="annotationFieldLabel"] {
        color: $text_muted;
        font-size: 10px;
        font-weight: 600;
    }
    QLabel[role="annotationFieldValue"] {
        color: #293b4e;
        font-size: 12px;
    }
    QFrame#AnnotationSequenceRow {
        background: $surface_subtle;
        border: 1px solid $border_soft;
        border-radius: 8px;
    }
    QLabel[role="annotationSequenceTitle"] {
        color: #34475a;
        font-weight: 600;
    }
    QLabel[role="annotationSequenceLength"] {
        color: $text_muted;
        font-size: 10px;
    }
    QLineEdit#AnnotationSequenceText {
        min-height: 18px;
        border: 1px solid #dce5ee;
        border-radius: 6px;
        padding: 5px 7px;
        background: $surface;
        color: #314154;
    }
    QToolButton#AnnotationCopyButton {
        min-height: 18px;
        border: 1px solid $border_strong;
        background: $surface;
        padding: 5px 9px;
        font-weight: 600;
    }
    QToolButton#AnnotationCopyButton:hover {
        background: $surface_hover;
        border-color: #9ab8f8;
    }
    QTabWidget#AnnotationDetailTabs::pane {
        border: none;
        background: transparent;
    }
    QScrollArea#AnnotationSummaryScroll,
    QWidget#AnnotationSummaryContent {
        border: none;
        background: transparent;
    }
    QTableWidget#AnnotationMarkerTable {
        border-radius: 8px;
    }
    QPlainTextEdit#AnnotationRawJson {
        border-radius: 9px;
        background: $surface;
    }
    QFrame#AnnotationActions {
        background: transparent;
        border: none;
        border-top: 1px solid $border;
    }
    QWidget#ContextCanvas, QWidget#ContextRowHeader {
        background: $surface;
    }
    QFrame#contextCanvasToolbar QLabel {
        color: $text_muted;
        font-weight: 600;
    }
    QFrame#TreeZoomPill {
        background: $surface;
        border: 1px solid $border_strong;
        border-radius: 8px;
    }
    QFrame#TreeZoomPill QToolButton {
        min-height: 0px;
        min-width: 0px;
        padding: 4px 7px;
        border: none;
        border-radius: 7px;
        background: transparent;
    }
    QFrame#TreeZoomPill QToolButton:hover {
        background: $primary_soft_hover;
    }
    QFrame#contextCanvasToolbar QToolButton,
    QFrame#TreeToolbar QToolButton {
        background: $surface;
        border-color: $border_strong;
    }
    QFrame#contextCanvasToolbar QToolButton:hover,
    QFrame#TreeToolbar QToolButton:hover {
        background: $primary_soft_hover;
        border-color: #9ab8f8;
    }
    QWidget#FilterRail {
        background: transparent;
    }
    /* The Phylogeny's Settings column laid over the tree in a narrow panel. */
    QWidget#FilterRail[overlay="true"] {
        background: $window;
        border-right: 1px solid $border_strong;
    }
    QScrollArea#SettingsScroll, QScrollArea#SettingsScroll > QWidget > QWidget,
    QWidget#SettingsContent {
        background: transparent;
    }
    QWidget#SettingsFooter {
        background: transparent;
        border-top: 1px solid transparent;
    }
    QWidget#SettingsFooter[divider="true"] {
        border-top: 1px solid $border;
    }
    QScrollArea#SettingsScroll QScrollBar:vertical {
        width: 6px;
        margin: 2px 0px;
    }
    QScrollArea#SettingsScroll QScrollBar::handle:vertical {
        background: $border;
        border-radius: 3px;
    }
    QScrollArea#SettingsScroll QScrollBar::handle:vertical:hover {
        background: $border_strong;
    }
    QGroupBox {
        font-weight: 600;
        border: 1px solid $border;
        border-radius: 12px;
        margin-top: 16px;
        background: #fbfcfe;
        color: #243447;
    }
    QGroupBox::title {
        subcontrol-origin: margin;
        subcontrol-position: top left;
        padding: 1px 8px;
        left: 12px;
        top: 0px;
        background: $window;
        color: #5a6b7c;
    }
    QLabel {
        color: #223246;
    }
    QLabel[role="muted"] {
        color: $text_muted;
    }
    QLabel[role="warningChip"] {
        background: #fff7e8;
        border: 1px solid #efd9ad;
        border-radius: 9px;
        color: $warning;
        font-size: 11px;
        font-weight: 600;
        padding: 2px 8px;
    }
    QLabel[role="dialogTitle"] {
        color: $text_strong;
        font-size: 15px;
        font-weight: 700;
    }
    QLabel[role="sectionTitle"], QLabel[role="chipTitle"] {
        color: $text_strong;
        font-weight: 600;
    }
    QLineEdit, QComboBox, QPlainTextEdit, QTextEdit,
    QSpinBox, QDoubleSpinBox, QDateEdit, QDateTimeEdit, QTimeEdit {
        border: 1px solid $border_strong;
        border-radius: 8px;
        padding: 6px 8px;
        background: $surface;
        color: $text_strong;
        selection-background-color: $primary_soft;
    }
    /* One drop-down look for the whole app: a 30 px white field with a
       plain muted chevron, and a list that drops below it. */
    QComboBox {
        min-height: 28px;
        padding: 0px 8px;
        combobox-popup: 0;
    }
    /* The "Table shows" switch above the annotation table stays compact. */
    QComboBox#TableViewBox {
        min-height: 22px;
    }
    /* The editor of an editable combo is part of the field, not a second box.
       Qt insets line edit text by 2 px, so the field pads 2 px less to keep
       its text in line with the other drop-downs. */
    QComboBox:editable {
        padding-left: 6px;
    }
    QComboBox QLineEdit {
        border: none;
        border-radius: 0px;
        padding: 0px;
        margin: 0px;
        background: transparent;
    }
    QComboBox::drop-down {
        subcontrol-origin: padding;
        subcontrol-position: center right;
        width: 22px;
        border: none;
        background: transparent;
    }
    QComboBox::down-arrow {
        image: url("$chevron_down");
        width: 10px;
        height: 6px;
    }
    QComboBox::down-arrow:disabled {
        image: url("$chevron_down_disabled");
    }
    QComboBox QAbstractItemView {
        background: $surface;
        border: 1px solid $border;
        border-radius: 8px;
        padding: 4px;
        outline: none;
        color: $text_strong;
    }
    QComboBox QAbstractItemView::item {
        min-height: 18px;
        padding: 6px 10px;
        border: none;
        border-radius: 6px;
    }
    /* No item color: the list's text_strong is the default, so a row's own
       foreground (a muted, not installed tag) shows. */
    QComboBox QAbstractItemView::item:hover,
    QComboBox QAbstractItemView::item:selected,
    QComboBox QAbstractItemView::item:focus {
        background: $primary_soft;
    }
    QComboBox QAbstractItemView::item:disabled {
        color: $text_disabled;
    }
    QComboBox QAbstractItemView::indicator {
        width: 14px;
        height: 14px;
        border: 1px solid $border_strong;
        border-radius: 4px;
        background: $surface;
    }
    QComboBox QAbstractItemView::indicator:checked {
        border-color: $primary;
        background: $primary;
        image: url("$check_tick");
    }
    QLineEdit:focus, QComboBox:focus, QPlainTextEdit:focus, QTextEdit:focus,
    QSpinBox:focus, QDoubleSpinBox:focus, QDateEdit:focus,
    QDateTimeEdit:focus, QTimeEdit:focus {
        border: 1px solid $focus;
        background: $surface;
    }
    QLineEdit:disabled, QComboBox:disabled, QPlainTextEdit:disabled,
    QTextEdit:disabled, QSpinBox:disabled, QDoubleSpinBox:disabled {
        background: #f5f7fa;
        border-color: #dde4ec;
        color: $text_disabled;
    }
    QSpinBox::up-button, QSpinBox::down-button,
    QDoubleSpinBox::up-button, QDoubleSpinBox::down-button {
        border: none;
        background: transparent;
    }
    QPushButton {
        min-height: 20px;
        border: 1px solid $border_strong;
        border-radius: 8px;
        padding: 6px 12px;
        background: $surface;
        color: #203041;
        font-weight: 500;
    }
    QPushButton:hover {
        background: $surface_hover;
        border-color: #9ab8f8;
    }
    QPushButton:pressed {
        background: $surface_pressed;
        border-color: #7ea3ef;
    }
    QPushButton:focus {
        border: 1px solid $focus;
    }
    QPushButton:disabled {
        background: #f5f7fa;
        border-color: #dde4ec;
        color: $text_disabled;
    }
    QPushButton[role="primary"] {
        background: $primary;
        border-color: $primary;
        color: $surface;
        font-weight: 600;
    }
    QPushButton[role="primary"]:hover {
        background: $primary_hover;
        border-color: $primary_hover;
    }
    QPushButton[role="danger"] {
        color: $danger;
        border-color: #e4b3b3;
    }
    QToolButton {
        min-height: 20px;
        border: 1px solid transparent;
        border-radius: 7px;
        padding: 5px 8px;
        background: transparent;
        color: #314154;
    }
    QToolButton:hover {
        background: $surface_hover;
        border-color: $border;
    }
    QToolButton:pressed, QToolButton:checked {
        background: $primary_soft;
        border-color: #9ab8f8;
        color: $text_strong;
    }
    QToolButton:focus {
        border-color: $focus;
    }
    QToolButton:disabled {
        color: $text_disabled;
    }
    QToolButton#SettingsChevron {
        min-height: 0px;
        min-width: 0px;
        padding: 0px;
        border: 1px solid transparent;
        border-radius: 6px;
        color: $text_muted;
        font-size: 17px;
        font-weight: 600;
    }
    QToolButton#SettingsChevron:hover {
        background: $surface_hover;
        border-color: $border;
        color: $text_strong;
    }
    QFrame#ContextSettingsStrip {
        background: transparent;
        border-right: 1px solid $border;
    }
    QLabel#SettingsStripBadge {
        padding: 0px 4px;
        border-radius: 8px;
        background: $primary;
        color: $surface;
        font-size: 10px;
        font-weight: 700;
    }
    QLabel[role="settingsCaption"] {
        color: $text_muted;
        font-size: 11px;
        font-weight: 600;
    }
    QWidget#SettingsContent QPushButton, QWidget#SettingsFooter QPushButton {
        min-height: 28px;
        max-height: 28px;
        padding: 0px 12px;
    }
    QWidget#SettingsContent QPushButton:checked, QWidget#SettingsFooter QPushButton:checked {
        background: $primary_soft;
        border-color: #9ab8f8;
        color: #1e4f9f;
    }
    QWidget#SettingsFooter QPushButton#ContextExportButton {
        padding: 0px 26px 0px 12px;
    }
    QWidget#SettingsFooter QPushButton#ContextExportButton::menu-indicator {
        image: url("$chevron_down");
        subcontrol-origin: border;
        subcontrol-position: center right;
        right: 9px;
        width: 10px;
        height: 6px;
    }
    QWidget#SettingsFooter QPushButton#ContextExportButton::menu-indicator:disabled {
        image: url("$chevron_down_disabled");
    }
    QWidget#SettingsContent QLineEdit {
        padding: 0px 8px;
    }
    QWidget#SettingsContent QLineEdit QToolButton {
        min-height: 0px;
        padding: 0px;
        border: none;
    }
    QWidget#SettingsContent QCheckBox {
        padding: 2px 0px;
    }
    QLabel[role="panelTitle"] {
        color: $text_strong;
        font-size: 14px;
        font-weight: 700;
        padding: 2px 2px 0px 2px;
    }
    QFrame#ContextRowsCard {
        background: $surface;
        border: 1px solid $border;
        border-radius: 10px;
    }
    QToolButton#ContextLinkButton {
        border: none;
        background: transparent;
        color: $primary;
        font-weight: 600;
        padding: 2px 4px;
    }
    QWidget#SettingsContent QToolButton#ContextLinkButton {
        padding: 2px 0px;
    }
    QToolButton#ContextLinkButton:hover {
        color: $primary_hover;
        text-decoration: underline;
    }
    QToolButton#ContextLinkButton:disabled {
        color: $text_disabled;
    }
    QToolButton[role="filterChip"] {
        min-width: 18px;
        max-width: 18px;
        padding: 2px;
        color: $text_muted;
    }
    QCheckBox, QRadioButton {
        spacing: 7px;
        padding: 3px 1px;
        color: #314154;
    }
    QCheckBox:hover, QRadioButton:hover,
    QCheckBox:focus, QRadioButton:focus {
        color: #1e4f9f;
    }
    QCheckBox:disabled, QRadioButton:disabled {
        color: $text_disabled;
    }
    QTableView, QTreeView, QListWidget {
        background: $surface;
        alternate-background-color: $surface_subtle;
        gridline-color: #e6edf3;
        selection-background-color: $primary_soft;
        selection-color: $text_strong;
        border: 1px solid $border;
        border-radius: 12px;
        outline: none;
    }
    QTableView:focus, QTreeView:focus, QListWidget:focus {
        border-color: $focus;
    }
    QHeaderView::section {
        background: $dialog;
        color: #5a6b7c;
        padding: 8px 10px;
        border: none;
        border-right: 1px solid $border_soft;
        border-bottom: 1px solid $border;
        font-weight: 600;
    }
    QHeaderView::section:hover {
        background: #edf3fa;
    }
    QTableCornerButton::section {
        background: $dialog;
        border: none;
        border-right: 1px solid $border_soft;
        border-bottom: 1px solid $border;
    }
    QTabWidget::pane {
        border: 1px solid $border;
        border-radius: 12px;
        background: $surface;
        top: -1px;
    }
    QTabBar::tab {
        background: #f2f5f9;
        border: 1px solid $border;
        color: #5b6c7d;
        padding: 8px 14px;
        border-top-left-radius: 10px;
        border-top-right-radius: 10px;
        margin-right: 4px;
    }
    QTabBar::tab:selected {
        background: $surface;
        color: $text;
        border-bottom-color: $surface;
    }
    QTabBar::tab:hover:!selected {
        background: #f7faff;
        color: #32475c;
    }
    QProgressBar {
        border: 1px solid #d3dbe6;
        border-radius: 8px;
        background: #f3f6fa;
        color: #223246;
        text-align: center;
    }
    QProgressBar::chunk {
        border-radius: 8px;
        background: #7ea6f8;
    }
    QMenuBar {
        background: $surface;
        border-bottom: 1px solid $border;
        padding: 2px 5px;
    }
    QMenuBar::item {
        border-radius: 6px;
        padding: 5px 9px;
        background: transparent;
    }
    QMenuBar::item:selected, QMenuBar::item:pressed {
        background: $primary_soft_hover;
    }
    QMenu {
        background: $surface;
        border: 1px solid $border;
        padding: 6px;
    }
    QMenu::separator {
        height: 1px;
        margin: 5px 8px;
        background: $border_soft;
    }
    QMenu::item {
        padding: 6px 18px;
        border-radius: 6px;
    }
    QMenu::item:selected {
        background: $primary_soft_hover;
    }
    QToolBar {
        background: $surface;
        border: none;
        border-bottom: 1px solid $border;
        spacing: 4px;
        padding: 4px 6px;
    }
    QToolBar::separator {
        background: $border;
        width: 1px;
        margin: 4px 6px;
    }
    QSplitter::handle {
        background: transparent;
    }
    QSplitter::handle:hover, QSplitter::handle:pressed {
        background: $focus_soft;
    }
    QSplitter::handle:horizontal {
        width: 6px;
        margin: 0px 1px;
    }
    QSplitter::handle:vertical {
        height: 6px;
        margin: 1px 0px;
    }
    QScrollArea {
        background: transparent;
        border: none;
    }
    QScrollBar:vertical {
        background: transparent;
        width: 12px;
        margin: 4px 0px;
    }
    QScrollBar::handle:vertical {
        background: $border_strong;
        border-radius: 6px;
        min-height: 28px;
    }
    QScrollBar::handle:vertical:hover {
        background: #b7c5d5;
    }
    QScrollBar:horizontal {
        background: transparent;
        height: 12px;
        margin: 0px 4px;
    }
    QScrollBar::handle:horizontal {
        background: $border_strong;
        border-radius: 6px;
        min-width: 28px;
    }
    QScrollBar::handle:horizontal:hover {
        background: #b7c5d5;
    }
    QScrollBar::add-line, QScrollBar::sub-line {
        width: 0px;
        height: 0px;
    }
    /* Docking area (Qt Advanced Docking System, its own stylesheet off).
       Panels sit in groups: a tab strip over the visible panel. */
    ads--CDockContainerWidget {
        background: $window;
    }
    ads--CDockContainerWidget ads--CDockSplitter::handle {
        background: transparent;
    }
    ads--CDockContainerWidget ads--CDockSplitter::handle:hover,
    ads--CDockContainerWidget ads--CDockSplitter::handle:pressed {
        background: $focus_soft;
    }
    ads--CDockAreaWidget {
        background: $surface;
        border: 1px solid $border;
        border-radius: 6px;
    }
    ads--CDockAreaWidget #tabsMenuButton::menu-indicator {
        image: none;
    }
    ads--CDockAreaTitleBar {
        background: $dialog;
        border-bottom: 1px solid $border;
        border-top-left-radius: 6px;
        border-top-right-radius: 6px;
        padding: 0px;
    }
    ads--CDockAreaWidget[focused="true"] ads--CDockAreaTitleBar {
        border-bottom: 1px solid $focus;
    }
    ads--CDockWidgetTab {
        background: transparent;
        border: none;
        border-right: 1px solid $border_soft;
        border-top: 2px solid transparent;
        padding: 0px 2px 0px 6px;
        min-height: 28px;
    }
    ads--CDockWidgetTab:hover {
        background: $surface_hover;
    }
    ads--CDockWidgetTab[activeTab="true"] {
        background: $surface;
        border-top: 2px solid $border_strong;
    }
    ads--CDockWidgetTab[focused="true"] {
        background: $surface;
        border-top: 2px solid $primary;
    }
    ads--CDockWidgetTab QLabel {
        color: $text_muted;
        background: transparent;
    }
    ads--CDockWidgetTab[activeTab="true"] QLabel,
    ads--CDockWidgetTab[focused="true"] QLabel {
        color: $text_strong;
    }
    ads--CDockWidget {
        background: $surface;
        border: none;
    }
    #tabCloseButton {
        margin-top: 2px;
        background: none;
        border: none;
        border-radius: 4px;
        padding: 0px;
        qproperty-icon: url(:/ads/images/close-button.svg),
            url(:/ads/images/close-button-disabled.svg) disabled;
        qproperty-iconSize: 14px;
    }
    #tabCloseButton:hover {
        background: $primary_soft_hover;
    }
    #tabsMenuButton {
        qproperty-icon: url(:/ads/images/tabs-menu-button.svg);
        qproperty-iconSize: 16px;
        border: none;
        background: none;
    }
    QToolButton#DockGroupAddButton, QToolButton#DockGroupMenuButton {
        border: none;
        border-radius: 4px;
        background: transparent;
        color: $text_muted;
        padding: 2px 7px;
        font-weight: 600;
    }
    QToolButton#DockGroupAddButton:hover, QToolButton#DockGroupMenuButton:hover {
        background: $primary_soft_hover;
        color: $text_strong;
    }
    QWidget#DockEmptyState {
        background: $window;
    }
    QStatusBar {
        min-height: 24px;
        background: $window;
        color: $text_muted;
        border-top: 1px solid $border;
    }
    QStatusBar::item {
        border: none;
    }
    QStatusBar QLabel {
        color: $text_muted;
        padding: 2px 8px;
    }
    /* Start window (run selection) and the Run progress window. */
    QDialog#StartupDialog, QDialog#RunProgressDialog, QDialog#NewRunDialog,
    QDialog#SettingsDialog, QDialog#AddSchemaDialog, QDialog#AddDatabaseDialog, QDialog#FirstTimeSetupDialog,
    QMessageBox#SetupStopQuestion {
        background: $window;
    }
    QLabel[role="windowTitle"] {
        color: $text_strong;
        font-size: 17px;
        font-weight: 700;
    }
    QLabel[role="windowHint"] {
        color: $text_muted;
        font-size: 11px;
    }
    QListWidget#RunList {
        background: transparent;
        border: none;
        padding: 0px;
    }
    QListWidget#RunList::item,
    QListWidget#RunList::item:hover,
    QListWidget#RunList::item:selected {
        background: transparent;
        border: none;
    }
    QFrame#RunCard, QFrame#RunStepsCard {
        background: $surface;
        border: 1px solid $border;
        border-radius: 10px;
    }
    QFrame#RunCard:hover {
        border-color: $focus_hover;
    }
    QFrame#RunCard[selected="true"] {
        border: 1px solid $focus;
    }
    QFrame#RunCard QLabel[role="runTitle"] {
        color: $text_strong;
        font-size: 13px;
        font-weight: 700;
    }
    QLabel[role="runMeta"] {
        color: $text_muted;
        font-size: 11px;
    }
    QLabel[role="runStep"] {
        color: $text;
        font-size: 11px;
    }
    QLabel[role="runCount"] {
        color: $text_muted;
        font-size: 11px;
        font-weight: 600;
    }
    QProgressBar[role="stepBar"] {
        min-height: 6px;
        max-height: 6px;
        border: none;
        border-radius: 3px;
        background: $border_soft;
    }
    QProgressBar[role="stepBar"]::chunk {
        border-radius: 3px;
        background: $primary;
    }
    QProgressBar[role="stepBar"][state="success"]::chunk {
        background: $success;
    }
    QProgressBar[role="stepBar"][state="failed"]::chunk {
        background: $danger;
    }
    QProgressBar[role="stepBar"][state="stopped"]::chunk {
        background: $tree_inactive;
    }
    QFrame#StartFooter, QFrame#ProgressFooter {
        background: transparent;
        border: none;
        border-top: 1px solid $border;
    }
    QFrame#RunCard QPushButton[role="primary"]:disabled,
    QFrame#RunCard QPushButton[role="danger"]:disabled {
        background: $disabled_fill;
        border-color: $disabled_border;
        color: $text_disabled;
    }
    QPushButton[role="flat"] {
        background: transparent;
        border-color: transparent;
        color: $text_muted;
    }
    QPushButton[role="flat"]:hover {
        background: $surface_hover;
        border-color: $border;
        color: $text_strong;
    }
    QPlainTextEdit#RunErrorDetails {
        background: $danger_soft;
        border: 1px solid $danger_border;
        color: $danger;
        font-size: 11px;
    }
    QScrollArea#RunStepsScroll, QWidget#RunStepsContent {
        background: transparent;
    }
    QWidget#RunStepRow {
        background: transparent;
        border-bottom: 1px solid $border_soft;
    }
    QLabel[role="stepMarker"] {
        min-width: 8px;
        max-width: 8px;
        min-height: 8px;
        max-height: 8px;
        border-radius: 5px;
        border: 1px solid $border_strong;
        background: $surface;
    }
    QLabel[role="stepMarker"][state="running"] {
        background: $primary;
        border-color: $primary;
    }
    QLabel[role="stepMarker"][state="done"] {
        background: $success;
        border-color: $success;
    }
    QLabel[role="stepMarker"][state="reused"] {
        background: $success_border;
        border-color: $success;
    }
    QLabel[role="stepMarker"][state="skipped"] {
        background: $border_strong;
        border-color: $border_strong;
    }
    QLabel[role="stepMarker"][state="failed"] {
        background: $danger;
        border-color: $danger;
    }
    QLabel[role="stepMarker"][state="partial"],
    QLabel[role="stepMarker"][state="halted"] {
        background: $warning;
        border-color: $warning;
    }
    QLabel[role="stepName"] {
        color: $text;
    }
    QLabel[role="stepName"][state="running"] {
        color: $text_strong;
        font-weight: 600;
    }
    QLabel[role="stepName"][state="waiting"] {
        color: $text_muted;
    }
    QLabel[role="stepState"] {
        color: $text_muted;
        font-size: 11px;
    }
    QLabel[role="stepState"][state="running"] {
        color: $primary;
        font-weight: 600;
    }
    QLabel[role="stepState"][state="failed"] {
        color: $danger;
        font-weight: 600;
    }
    QLabel[role="stepState"][state="partial"],
    QLabel[role="stepState"][state="halted"] {
        color: $warning;
    }
    QLabel[role="stepDetail"] {
        color: $text_muted;
        font-size: 11px;
    }
    QLabel[role="stepNote"] {
        color: $danger;
        font-size: 11px;
    }
    QToolButton#SetupDetailsToggle {
        border: none;
        background: transparent;
        color: $text_muted;
        padding: 2px 0px;
    }
    QToolButton#SetupDetailsToggle:hover {
        color: $text_strong;
    }
    QPlainTextEdit#SetupLog {
        background: $surface;
        color: $text;
    }
    QToolTip {
        background: $text;
        color: $surface;
        border: none;
        padding: 6px 8px;
    }
    """
).substitute(COLORS, **_IMAGES)


def _application_palette() -> QtGui.QPalette:
    palette = QtGui.QPalette()
    palette.setColor(QtGui.QPalette.ColorRole.Window, QtGui.QColor(COLORS["window"]))
    palette.setColor(QtGui.QPalette.ColorRole.WindowText, QtGui.QColor(COLORS["text"]))
    palette.setColor(QtGui.QPalette.ColorRole.Base, QtGui.QColor(COLORS["surface"]))
    palette.setColor(QtGui.QPalette.ColorRole.AlternateBase, QtGui.QColor(COLORS["surface_subtle"]))
    palette.setColor(QtGui.QPalette.ColorRole.ToolTipBase, QtGui.QColor(COLORS["text"]))
    palette.setColor(QtGui.QPalette.ColorRole.ToolTipText, QtGui.QColor(COLORS["surface"]))
    palette.setColor(QtGui.QPalette.ColorRole.Text, QtGui.QColor(COLORS["text_strong"]))
    palette.setColor(QtGui.QPalette.ColorRole.Button, QtGui.QColor(COLORS["surface"]))
    palette.setColor(QtGui.QPalette.ColorRole.ButtonText, QtGui.QColor(COLORS["text"]))
    palette.setColor(QtGui.QPalette.ColorRole.Highlight, QtGui.QColor(COLORS["primary_soft"]))
    palette.setColor(QtGui.QPalette.ColorRole.HighlightedText, QtGui.QColor(COLORS["text_strong"]))
    palette.setColor(QtGui.QPalette.ColorRole.Link, QtGui.QColor(COLORS["primary"]))
    palette.setColor(QtGui.QPalette.ColorRole.PlaceholderText, QtGui.QColor(COLORS["text_muted"]))
    palette.setColor(
        QtGui.QPalette.ColorGroup.Disabled,
        QtGui.QPalette.ColorRole.Text,
        QtGui.QColor(COLORS["text_disabled"]),
    )
    palette.setColor(
        QtGui.QPalette.ColorGroup.Disabled,
        QtGui.QPalette.ColorRole.ButtonText,
        QtGui.QColor(COLORS["text_disabled"]),
    )
    return palette


def _x11_compositor_running() -> bool:
    """True when an X11 compositing manager owns the _NET_WM_CM_S<screen>
    selection; without one, see-through windows are drawn black."""

    import ctypes
    import ctypes.util

    try:
        name = ctypes.util.find_library("X11")
        if not name:
            return False
        xlib = ctypes.cdll.LoadLibrary(name)
        xlib.XOpenDisplay.restype = ctypes.c_void_p
        xlib.XOpenDisplay.argtypes = [ctypes.c_char_p]
        xlib.XDefaultScreen.argtypes = [ctypes.c_void_p]
        xlib.XInternAtom.restype = ctypes.c_ulong
        xlib.XInternAtom.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int]
        xlib.XGetSelectionOwner.restype = ctypes.c_ulong
        xlib.XGetSelectionOwner.argtypes = [ctypes.c_void_p, ctypes.c_ulong]
        xlib.XCloseDisplay.argtypes = [ctypes.c_void_p]
        display = xlib.XOpenDisplay(None)
        if not display:
            return False
        try:
            screen = xlib.XDefaultScreen(display)
            atom = xlib.XInternAtom(display, f"_NET_WM_CM_S{screen}".encode(), 0)
            return bool(atom) and xlib.XGetSelectionOwner(display, atom) != 0
        finally:
            xlib.XCloseDisplay(display)
    except (OSError, AttributeError):
        return False


_SEE_THROUGH_POPUPS: bool | None = None


def _see_through_popups_supported() -> bool:
    """Whether combo lists can sit on a see-through window so their rounded
    corners show. Wayland, Windows and macOS always composite; plain X11 only
    does with a compositing manager running."""

    global _SEE_THROUGH_POPUPS
    if _SEE_THROUGH_POPUPS is None:
        if QtGui.QGuiApplication.platformName() == "xcb":
            _SEE_THROUGH_POPUPS = _x11_compositor_running()
        else:
            _SEE_THROUGH_POPUPS = True
    return _SEE_THROUGH_POPUPS


class _ComboWheelGuard(QtCore.QObject):
    """Lets the mouse wheel scroll the page past a combo box instead of
    changing its value, unless the combo has keyboard focus."""

    def eventFilter(self, watched: QtCore.QObject, event: QtCore.QEvent) -> bool:  # type: ignore[override]
        if (
            event.type() != QtCore.QEvent.Type.Wheel
            or not isinstance(watched, QtWidgets.QComboBox)
            or watched.hasFocus()
        ):
            return False
        area = watched.parentWidget()
        while area is not None and not isinstance(area, QtWidgets.QAbstractScrollArea):
            area = area.parentWidget()
        if area is None:
            event.ignore()
            return True
        viewport = area.viewport()
        position = viewport.mapFromGlobal(event.globalPosition())
        forwarded = QtGui.QWheelEvent(
            position,
            event.globalPosition(),
            event.pixelDelta(),
            event.angleDelta(),
            event.buttons(),
            event.modifiers(),
            event.phase(),
            event.inverted(),
        )
        QtCore.QCoreApplication.sendEvent(viewport, forwarded)
        event.accept()
        return True


_COMBO_WHEEL_GUARD = _ComboWheelGuard()


class _AppStyle(QtWidgets.QProxyStyle):
    """Fusion, with two app-wide combo box rules: the wheel only changes a
    focused combo, and the drop-down list is drawn on a see-through window
    (where the platform composites) so the style sheet's rounded corners show."""

    def polish(self, target):  # type: ignore[override]
        if isinstance(target, QtWidgets.QComboBox):
            if target.focusPolicy() == QtCore.Qt.FocusPolicy.WheelFocus:
                target.setFocusPolicy(QtCore.Qt.FocusPolicy.StrongFocus)
            target.installEventFilter(_COMBO_WHEEL_GUARD)
            if _see_through_popups_supported():
                # The list's window exists once view() is asked for, and must
                # be see-through before it is first shown.
                target.view().window().setAttribute(
                    QtCore.Qt.WidgetAttribute.WA_TranslucentBackground, True
                )
        return super().polish(target)


_EXIT_CLEANUP_REGISTERED = False


def _delete_python_owned_windows() -> None:
    """Delete the parentless widgets Python owns before PyQt's exit cleanup does.

    PyQt deletes them at exit while it walks its table of wrapped objects.
    Deleting one runs Python code (this style's polish, event filters, slots)
    and each object that code wraps is added to that table; when the table
    grows it is reallocated under the walk and the process crashes on exit
    (pyqtgraph's plot menus and closed dock widgets do this). Deleting them
    here, before the walk starts, is the same cleanup without the crash.
    """
    if QtWidgets.QApplication.instance() is None:
        return
    for widget in QtWidgets.QApplication.topLevelWidgets():
        if not sip.isdeleted(widget) and sip.ispyowned(widget) and widget.parent() is None:
            sip.delete(widget)


def _register_exit_cleanup() -> None:
    # atexit runs last registered first, and PyQt registered its own cleanup
    # when QtCore was imported, so this one runs before it.
    global _EXIT_CLEANUP_REGISTERED
    if not _EXIT_CLEANUP_REGISTERED:
        atexit.register(_delete_python_owned_windows)
        _EXIT_CLEANUP_REGISTERED = True


def apply_application_theme(app: QtWidgets.QApplication) -> None:
    """Install app identity and the shared visual system before any dialog opens."""

    QtCore.QCoreApplication.setOrganizationName(APP_ORGANIZATION)
    QtCore.QCoreApplication.setOrganizationDomain(APP_DOMAIN)
    QtCore.QCoreApplication.setApplicationName(APP_NAME)
    app.setApplicationDisplayName(APP_DISPLAY_NAME)

    # Installed once: replacing the style repolishes every widget alive.
    if "Fusion" in QtWidgets.QStyleFactory.keys() and not isinstance(app.style(), _AppStyle):
        app.setStyle(_AppStyle("Fusion"))
        _register_exit_cleanup()
    app.setPalette(_application_palette())

    # One point larger than the platform font, once: applying the theme again
    # (as the tests do per test class) must not keep growing the text.
    if not app.property("bactolomewFontRaised"):
        font = app.font()
        if font.pointSize() > 0:
            font.setPointSize(font.pointSize() + 1)
        app.setFont(font)
        app.setProperty("bactolomewFontRaised", True)
    app.setStyleSheet(APP_STYLESHEET)
