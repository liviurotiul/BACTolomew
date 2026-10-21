"""Qt Advanced Docking System configuration for every explore window (plan 3.12).

This module holds data only and never imports PyQt6Ads itself: ``dock_host``
(the one module that imports the library) passes the module in through
``apply_global_config`` and ``apply_dock_features``. The flag names below were
confirmed against PyQt6Ads 5.1.1 in the Phase 0 docking spike, under Xvfb and
then on the user's real display (36 of 36 checks). Phase 1's checklist C on the
real display is still open, so keep every change to them in this file.

PyQt6Ads specifics:
- the config flags are class level and must be set before the first
  ``CDockManager`` exists;
- ``eConfigFlag`` values do not combine with ``|`` in PyQt6Ads, so start from a
  known base and set each flag with its own ``setConfigFlag`` call;
- the default config at import (DefaultNonOpaqueConfig) turns on
  ``DoubleClickUndocksWidget`` and ``DockAreaHasUndockButton``; both must be off.
"""
from __future__ import annotations

from typing import Any, Tuple

# Starting point for the explicit flags below.
BASE_CONFIG = "DefaultBaseConfig"

CONFIG_FLAGS_ON: Tuple[str, ...] = (
    "OpaqueSplitterResize",  # dividers resize live
    "FocusHighlighting",  # the active panel's tab is marked
    "AllTabsHaveCloseButton",
    "DockAreaHasTabsMenuButton",  # ADS's own tab list, for overflowing strips
    "DockAreaDynamicTabsMenuButtonVisibility",  # ... shown only when tabs overflow
    "MiddleMouseButtonClosesTab",
    "EqualSplitOnInsertion",
    "DragPreviewIsDynamic",  # non-opaque drag: a preview moves, never the widget
    "DragPreviewShowsContentPixmap",
    "DisableStylesheet",  # theme.py styles the ads-- selectors
)

CONFIG_FLAGS_OFF: Tuple[str, ...] = (
    "DockAreaHasCloseButton",
    "DockAreaHasUndockButton",
    "DoubleClickUndocksWidget",  # on by default; double-click maximizes instead
    "ActiveTabHasCloseButton",  # redundant with AllTabsHaveCloseButton
    "FloatingContainerHasWidgetTitle",  # irrelevant without floating
)

# Per panel (every CDockWidget). Movable, Closable and Focusable stay on;
# FocusHighlighting needs Focusable.
DOCK_FEATURES_OFF: Tuple[str, ...] = (
    "DockWidgetFloatable",  # no panel can ever leave its explore window
    "DockWidgetPinnable",  # no auto-hide
)
DOCK_FEATURES_ON: Tuple[str, ...] = (
    "DockWidgetMovable",
    "DockWidgetClosable",
    "DockWidgetFocusable",
)

# Extra instances of multi-instance panels (plan 3.12): the close asks the
# owner first (it disposes the panel) and the dock is deleted, not hidden.
# Single-instance panels (and the first instance of each kind) are hidden on
# close and reopen in place.
CUSTOM_CLOSE_FEATURES: Tuple[str, ...] = (
    "CustomCloseHandling",
    "DockWidgetDeleteOnClose",
)

_applied = False


def apply_global_config(ads: Any) -> None:
    """Set the manager-wide flags once, before the first CDockManager."""
    global _applied
    if _applied:
        return
    manager = ads.CDockManager
    flags = manager.eConfigFlag
    manager.setConfigFlags(getattr(flags, BASE_CONFIG))
    for name in CONFIG_FLAGS_ON:
        manager.setConfigFlag(getattr(flags, name), True)
    for name in CONFIG_FLAGS_OFF:
        manager.setConfigFlag(getattr(flags, name), False)
    # Auto-hide stays off: autoHideConfigFlags() is empty by default and
    # AutoHideFeatureEnabled is never set.
    _applied = True


def apply_dock_features(ads: Any, dock: Any, *, custom_close: bool = False) -> None:
    """Give one CDockWidget the fixed per-panel features. ``custom_close``: the
    owner handles the close (disposes the panel, then deletes the dock)."""
    features = ads.CDockWidget.DockWidgetFeature
    for name in DOCK_FEATURES_OFF:
        dock.setFeature(getattr(features, name), False)
    for name in DOCK_FEATURES_ON:
        dock.setFeature(getattr(features, name), True)
    for name in CUSTOM_CLOSE_FEATURES:
        dock.setFeature(getattr(features, name), bool(custom_close))
