"""Docking area of the explore window (plan 3 and 6.3).

Only ``dock_host`` imports PyQt6Ads. The explore window needs it: there is no
fallback layout any more (the one-release legacy layout went in Phase 10), so
the app checks ``docking_available()`` at startup and says how to install it.
"""
from __future__ import annotations

import shlex
import sys

# What requirements.txt asks for (D1).
ADS_REQUIREMENT = "PyQt6Ads>=5.1,<6"


class DockingUnavailableError(RuntimeError):
    """PyQt6Ads (the docking library of every explore window) is missing."""


def docking_available() -> bool:
    from .dock_host import ads_available

    return ads_available()


def docking_import_error() -> str:
    """Why PyQt6Ads could not be imported ("" when it was)."""
    from .dock_host import ads_import_error

    return ads_import_error()


# The first line of the startup error; the message box shows it as its main
# text and docking_install_hint() below it.
DOCKING_UNAVAILABLE_HEADLINE = (
    "BACTolomew needs the docking library PyQt6Ads (Qt Advanced Docking System) "
    "for its explore windows, and it could not be loaded."
)


def docking_install_hint() -> str:
    """How to install PyQt6Ads into this Python, plus the import error."""
    reason = docking_import_error()
    # Quoted so the line still works when pasted for a path with spaces.
    lines = [
        "Install it into the Python environment that runs BACTolomew:",
        f"    {shlex.quote(sys.executable)} -m pip install {shlex.quote(ADS_REQUIREMENT)}",
        "(or: pip install -r requirements.txt)",
    ]
    if reason:
        lines += ["", f"Import error: {reason}"]
    return "\n".join(lines)


def docking_unavailable_message() -> str:
    """The startup error text: what is missing and how to install it."""
    return DOCKING_UNAVAILABLE_HEADLINE + "\n\n" + docking_install_hint()


def require_docking() -> None:
    """Raise ``DockingUnavailableError`` (with the install hint) without PyQt6Ads."""
    if not docking_available():
        raise DockingUnavailableError(docking_unavailable_message())
