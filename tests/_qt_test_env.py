"""Keeps the test suite away from the user's real QSettings.

Every test module that creates a QApplication or a window calls
``isolate_qt_settings()`` at import time, before any QSettings is created.
All modules of one test process share one scratch folder (its path is kept in
an environment variable, so a second import of this helper under another
module name reuses it); the folder is removed when the process exits.
"""
from __future__ import annotations

import atexit
import os
import shutil
import tempfile

from PyQt6 import QtCore

_ENV_KEY = "BACTOLOMEW_TEST_SETTINGS_DIR"


def isolate_qt_settings() -> str:
    """Point every QSettings format and scope at a per-process scratch folder."""
    path = os.environ.get(_ENV_KEY)
    if not path or not os.path.isdir(path):
        path = tempfile.mkdtemp(prefix="bactolomew_test_settings_")
        os.environ[_ENV_KEY] = path
        atexit.register(shutil.rmtree, path, True)
    fmt = QtCore.QSettings.Format
    scope = QtCore.QSettings.Scope
    QtCore.QSettings.setDefaultFormat(fmt.IniFormat)
    for each_format in (fmt.IniFormat, fmt.NativeFormat):
        for each_scope in (scope.UserScope, scope.SystemScope):
            QtCore.QSettings.setPath(each_format, each_scope, path)
    return path
