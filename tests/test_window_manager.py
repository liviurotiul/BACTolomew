"""AppController, WindowManager, commands and info windows (window plan Phase 5, 8.6).

- Info cards are separate top-level windows, one per entity: opening the same
  card again raises the open one (D6).
- A card is bound to the explore window that opened it; when that window
  closes, the card rebinds to the most recently active explore window.
- The app quits after the last explore window closes; cards never keep it
  alive.
- The sample card is only a card: no commands, no action buttons (D10).
- The gene card keeps its actions, which run as commands on its owner window.
"""
from __future__ import annotations

import os
import shutil
import tempfile
import time
import unittest
from pathlib import Path
from typing import Any, Callable, List, Optional
from unittest import mock

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from tests._qt_test_env import isolate_qt_settings
except ImportError:  # "unittest discover -s tests" puts tests/ itself on sys.path
    from _qt_test_env import isolate_qt_settings

isolate_qt_settings()

import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets, sip

try:
    from tests.test_dock_host import _app, _Explorer, _spin
    from tests.test_refresh_parity import build_synthetic_run
except ImportError:
    from test_dock_host import _app, _Explorer, _spin
    from test_refresh_parity import build_synthetic_run

from annotation_explorer_polars import annotation_detail_dialog
from annotation_explorer_polars.windows import explore_window
from annotation_explorer_polars.annotation_detail_dialog import (
    AnnotationDetailDialog,
    GeneCardAlignmentDialog,
)
from annotation_explorer_polars.panels.alignment import AlignmentPanel
from annotation_explorer_polars.panels.neighbourhood import NeighbourhoodPanel
from annotation_explorer_polars.model import RunData
from annotation_explorer_polars.sample_card_dialog import SampleCardDialog
from annotation_explorer_polars.session.requests import OpenPlasmidContext, OpenSampleCard
from annotation_explorer_polars.shell import commands as cmd
from annotation_explorer_polars.shell.app_controller import AppController
from annotation_explorer_polars.shell.commands import CommandRegistry
from annotation_explorer_polars.shell.window_manager import WindowManager
from annotation_explorer_polars.windows.info_window import InfoKey, InfoWindow, gene_key, sample_key


# ---------------------------------------------------------------- fakes


class _FakeExplore(QtWidgets.QMainWindow):
    """Stands in for an explore window: emits windowClosed after closing."""

    windowClosed = QtCore.pyqtSignal()

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        super().closeEvent(event)
        if event.isAccepted():
            self.windowClosed.emit()


class _FakeCard(InfoWindow):
    COMMANDS = ("demo.run",)

    def __init__(self, title: str = "Card") -> None:
        super().__init__(None)
        self.set_card_title(title)
        self.owner_changes = 0

    def _owner_changed(self) -> None:
        self.owner_changes += 1


def _alive(obj: Any) -> bool:
    return obj is not None and not sip.isdeleted(obj)


class WindowManagerTests(unittest.TestCase):
    """The WindowManager on its own, with fake explore windows and cards."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.app = _app()

    def setUp(self) -> None:
        self.quits: List[int] = []
        self.commands = CommandRegistry()
        self.ran: List[Any] = []
        self.commands.register("demo.run", lambda window, **kw: self.ran.append((window, kw)))
        self.wm = WindowManager(self.commands, quit_app=lambda: self.quits.append(1))
        self.windows: List[_FakeExplore] = []

    def tearDown(self) -> None:
        for card in self.wm.info_windows():
            card.close()
        for window in self.windows:
            if _alive(window):
                window.deleteLater()
        _spin(self.app, 0.05)

    def _explore(self, title: str = "Explore") -> _FakeExplore:
        window = _FakeExplore()
        window.setWindowTitle(title)
        window.resize(400, 300)
        window.show()
        self.windows.append(window)
        self.wm.add_explore(window)
        return window

    def _open(self, key: InfoKey, owner: Optional[QtWidgets.QWidget] = None, title: str = "Card") -> Any:
        built: List[_FakeCard] = []

        def factory() -> _FakeCard:
            built.append(_FakeCard(title))
            return built[0]

        card = self.wm.open_info(key, factory, owner=owner)
        return card, built

    def test_windows_are_numbered_and_the_newest_is_active(self) -> None:
        one = self._explore("one")
        two = self._explore("two")
        self.assertEqual((self.wm.number(one), self.wm.number(two)), (1, 2))
        self.assertEqual(self.wm.explore_windows(), [one, two])
        self.assertIs(self.wm.active_explore(), two)
        two.close()
        three = self._explore("three")
        self.assertEqual(self.wm.number(three), 2)  # the lowest free number

    def test_dedupe_raises_the_existing_card(self) -> None:
        owner = self._explore()
        first, built_first = self._open(sample_key("S1"), owner)
        second, built_second = self._open(sample_key("S1"), owner)
        self.assertIs(first, second)
        self.assertEqual((len(built_first), len(built_second)), (1, 0))
        self.assertEqual(self.wm.info_windows(), [first])
        self.assertTrue(first.isVisible())
        other, _ = self._open(sample_key("S2"), owner)
        self.assertIsNot(other, first)
        self.assertEqual(len(self.wm.info_windows()), 2)

    def test_a_closed_card_opens_anew(self) -> None:
        owner = self._explore()
        first, _ = self._open(sample_key("S1"), owner)
        first.close()
        _spin(self.app, 0.05)
        self.assertEqual(self.wm.info_windows(), [])
        self.assertFalse(_alive(first))  # deleted on close
        second, built = self._open(sample_key("S1"), owner)
        self.assertEqual(len(built), 1)
        self.assertIs(self.wm.info_window(sample_key("S1")), second)

    def test_close_info(self) -> None:
        owner = self._explore()
        card, _ = self._open(sample_key("S1"), owner)
        self.assertTrue(self.wm.close_info(sample_key("S1")))
        _spin(self.app, 0.05)
        self.assertFalse(self.wm.close_info(sample_key("S1")))
        self.assertFalse(_alive(card))

    def test_cards_are_top_level_windows_never_children_of_an_explore_window(self) -> None:
        owner = self._explore()
        card, _ = self._open(sample_key("S1"), owner)
        self.assertTrue(card.isWindow())
        self.assertIsNone(card.parent())
        self.assertEqual(card.windowType(), QtCore.Qt.WindowType.Window)
        self.assertEqual(card.windowModality(), QtCore.Qt.WindowModality.NonModal)

    def test_reopening_from_another_window_rebinds_the_card(self) -> None:
        one = self._explore("one")
        two = self._explore("two")
        card, _ = self._open(sample_key("S1"), one)
        self.assertIs(card.owner(), one)
        again, built = self._open(sample_key("S1"), two)
        self.assertIs(again, card)
        self.assertEqual(built, [])
        self.assertIs(card.owner(), two)

    def test_orphaned_cards_rebind_to_the_most_recently_active_window(self) -> None:
        one = self._explore("one")
        two = self._explore("two")
        three = self._explore("three")
        card, _ = self._open(sample_key("S1"), two)
        other, _ = self._open(sample_key("S2"), one)
        self.assertEqual(card.windowTitle(), "Card [2]")  # several windows: owner badge
        # Focus moves to window one, then three: three is the most recent.
        self.app.focusChanged.emit(None, one)
        self.app.focusChanged.emit(None, three)
        self.assertIs(self.wm.active_explore(), three)
        two.close()
        _spin(self.app, 0.05)
        self.assertIs(card.owner(), three)
        self.assertEqual(card.windowTitle(), "Card [3]")
        self.assertIs(other.owner(), one)  # cards of other windows keep their owner
        self.assertTrue(card.isVisible())
        three.close()
        _spin(self.app, 0.05)
        self.assertIs(card.owner(), one)
        # One window left: no badge.
        self.assertEqual(card.windowTitle(), "Card")
        self.assertEqual(self.quits, [])

    def test_focus_inside_a_window_marks_it_active(self) -> None:
        one = self._explore("one")
        two = self._explore("two")
        child = QtWidgets.QLineEdit(one)
        self.app.focusChanged.emit(None, child)
        self.assertIs(self.wm.active_explore(), one)
        # Focus in a card does not change the active explore window.
        card, _ = self._open(sample_key("S1"), two)
        self.app.focusChanged.emit(None, card)
        self.assertIs(self.wm.active_explore(), one)

    def test_the_app_quits_after_the_last_explore_window_closes(self) -> None:
        one = self._explore("one")
        two = self._explore("two")
        card, _ = self._open(sample_key("S1"), one)
        one.close()
        _spin(self.app, 0.05)
        self.assertEqual(self.quits, [])
        self.assertTrue(_alive(card) and card.isVisible())
        two.close()
        _spin(self.app, 0.05)
        self.assertEqual(self.quits, [1])
        self.assertFalse(_alive(card))
        self.assertEqual(self.wm.info_windows(), [])
        self.assertEqual(self.wm.explore_windows(), [])

    def test_cards_never_keep_the_app_alive_and_closing_them_never_quits(self) -> None:
        owner = self._explore()
        card, _ = self._open(sample_key("S1"), owner)
        card.close()
        _spin(self.app, 0.05)
        self.assertEqual(self.quits, [])
        self._open(sample_key("S2"), owner)
        owner.close()
        _spin(self.app, 0.05)
        self.assertEqual(self.quits, [1])
        self.assertEqual(self.wm.info_windows(), [])

    def test_a_destroyed_explore_window_counts_as_closed(self) -> None:
        one = self._explore("one")
        two = self._explore("two")
        card, _ = self._open(sample_key("S1"), one)
        sip.delete(one)
        self.assertEqual(self.wm.explore_windows(), [two])
        self.assertIs(card.owner(), two)

    def test_there_are_no_transient_windows_any_more(self) -> None:
        # Phase 10: document panels are always panels of the docking area.
        self.assertFalse(hasattr(self.wm, "show_transient"))
        self.assertFalse(hasattr(self.wm, "transient_windows"))

    def test_card_commands_run_on_the_owner_window(self) -> None:
        one = self._explore("one")
        two = self._explore("two")
        card, _ = self._open(sample_key("S1"), one)
        self.assertTrue(card.run_command("demo.run", value=3))
        one.close()
        _spin(self.app, 0.05)
        card.run_command("demo.run", value=4)
        self.assertEqual(self.ran, [(one, {"value": 3}), (two, {"value": 4})])
        with self.assertRaises(ValueError):
            card.run_command("info.open_sample", sample="S1")  # not one of its commands

    def test_window_menu_lists_windows_and_cards_without_shortcuts(self) -> None:
        one = self._explore("one")
        two = self._explore("two")
        self._open(sample_key("S1"), one, title="Sample card: S1")
        menu = QtWidgets.QMenu()
        self.wm.fill_window_menu(menu, two)
        texts = [a.text() for a in menu.actions()]
        self.assertEqual(texts, ["Explore windows", "one [1]", "two [2]", "", "Info cards", "Sample card: S1 [1]"])
        checked = [a.text() for a in menu.actions() if a.isChecked()]
        self.assertEqual(checked, ["two [2]"])
        # The section labels are drawn (disabled bold entries, not addSection,
        # whose text the app's style does not show).
        headers = [a for a in menu.actions() if a.text() in ("Explore windows", "Info cards")]
        self.assertTrue(all(not a.isEnabled() and a.font().bold() and not a.isSeparator() for a in headers))
        self.assertTrue(all(a.shortcut().isEmpty() for a in menu.actions()))
        for card in self.wm.info_windows():
            card.close()
        _spin(self.app, 0.05)
        self.wm.fill_window_menu(menu, two)
        empty = [a for a in menu.actions() if a.text() == "No open cards"]
        self.assertEqual(len(empty), 1)
        self.assertFalse(empty[0].isEnabled())

    def test_open_explore_uses_the_factory(self) -> None:
        made: List[_FakeExplore] = []

        def factory() -> _FakeExplore:
            window = _FakeExplore()
            made.append(window)
            self.windows.append(window)
            return window

        wm = WindowManager(self.commands, window_factory=factory)
        window = wm.open_explore()
        self.assertEqual(made, [window])
        self.assertEqual(wm.explore_windows(), [window])
        with self.assertRaises(RuntimeError):
            WindowManager(self.commands).open_explore()


class InfoKeyAndCommandTests(unittest.TestCase):
    def test_gene_keys_are_the_same_for_the_same_gene_from_any_opener(self) -> None:
        table_row = {"Sample": "S1", "Contig": "3", "Locus": "S1_3_000", "Start": 1, "Stop": 900, "Strand": "+"}
        self.assertEqual(gene_key(table_row), InfoKey("gene", ("S1", "3", "S1_3_000")))
        no_locus = {"Sample": "S1", "Contig": "3", "Start": "1", "Stop": 900.0, "Strand": "+", "Locus": float("nan")}
        self.assertEqual(gene_key(no_locus), InfoKey("gene", ("S1", "3", "1-900+")))
        self.assertEqual(gene_key({**no_locus, "Start": 1.0}), gene_key(no_locus))
        self.assertNotEqual(gene_key({**no_locus, "Strand": "-"}), gene_key(no_locus))
        self.assertEqual(sample_key(" S1 "), InfoKey("sample", ("S1",)))

    def test_rows_without_locus_or_coordinates_never_share_a_gene_key(self) -> None:
        by_id = {"Sample": "S1", "AnnotationID": "S1:1:7", "Gene": "blaOXA"}
        self.assertEqual(gene_key(by_id), InfoKey("gene", ("S1", "", "id:S1:1:7")))
        self.assertNotEqual(gene_key({**by_id, "AnnotationID": "S1:1:8"}), gene_key(by_id))
        loose = {"Sample": "S1", "Gene": "blaOXA"}
        self.assertNotEqual(gene_key(loose), gene_key({"Sample": "S1", "Gene": "blaKPC"}))
        self.assertEqual(gene_key(dict(reversed(list(loose.items())))), gene_key(loose))

    def test_registry(self) -> None:
        registry = CommandRegistry()
        registry.register("a.b", lambda window, **kw: (window, kw), "A B")
        self.assertEqual(registry.run("a.b", "W", x=1), ("W", {"x": 1}))
        self.assertEqual(registry.names(), ["a.b"])
        self.assertEqual(registry.title("a.b"), "A B")
        with self.assertRaises(ValueError):
            registry.register("a.b", lambda window: None)
        with self.assertRaises(KeyError):
            registry.run("missing", "W")
        with self.assertRaises(ValueError):
            registry.run("a.b", None)


# ---------------------------------------------------------------- real windows


def _wait(app: QtWidgets.QApplication, predicate: Callable[[], bool], timeout: float = 10.0) -> bool:
    end = time.monotonic() + timeout
    while time.monotonic() < end:
        if predicate():
            return True
        _spin(app, 0.02)
    return bool(predicate())


def _document_panels(win: Any, kind: str) -> List[Any]:
    """The open document panels (neighbourhood or alignment) of ``win``."""
    manager = win._panel_manager
    return [manager.panel(f.uid) for f in manager.frames(kind) if manager.is_open(f.uid) and manager.panel(f.uid)]


class _RealWindowCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = _app()
        import annotation_explorer_polars.geo_experimental_views as gev

        cls._saved_webengine_cache = gev._WEBENGINE_RUNTIME_CACHE
        gev._WEBENGINE_RUNTIME_CACHE = (False, "disabled in tests")

    @classmethod
    def tearDownClass(cls) -> None:
        import annotation_explorer_polars.geo_experimental_views as gev

        gev._WEBENGINE_RUNTIME_CACHE = cls._saved_webengine_cache

    def setUp(self) -> None:
        self.output_dir = tempfile.mkdtemp(prefix="bactolomew_wm_run_")
        self.windows: List[Any] = []
        settings = QtCore.QSettings()
        settings.remove("explore")
        settings.sync()

    def tearDown(self) -> None:
        for win in self.windows:
            if _alive(win):
                win.close()
                win.deleteLater()
        _spin(self.app, 0.05)
        shutil.rmtree(self.output_dir, ignore_errors=True)
        settings = QtCore.QSettings()
        settings.remove("explore")
        settings.sync()

    def _window(self) -> Any:
        win = _Explorer.build(self.output_dir)
        self.windows.append(win)
        _Explorer.show_and_settle(self.app, win, seconds=0.2)
        return win

    def _run_data(self) -> RunData:
        ann_df, meta_df = build_synthetic_run()
        return RunData.build(
            ann_df, meta_df, ref_df=pd.DataFrame(), marker_hits_df=pd.DataFrame(), output_dir=Path(self.output_dir)
        )

    def _settle(self, win: Any) -> None:
        done = {"v": win.ctx.first_result_done}
        win.firstRefreshFinished.connect(lambda: done.__setitem__("v", True))
        _wait(self.app, lambda: done["v"], 20)
        _spin(self.app, 0.2)

    @staticmethod
    def _row(win: Any, annotation_id: str = "S1:1:0") -> dict:
        return win.run.ann_df.loc[win.run.ann_df["AnnotationID"] == annotation_id].iloc[0].to_dict()


class AppControllerTests(_RealWindowCase):
    def test_the_controller_owns_the_run_and_its_windows(self) -> None:
        run = self._run_data()
        controller = AppController(run)
        first = controller.open_first_window()
        self.windows.append(first)
        self._settle(first)
        second = controller.open_explore()
        self.windows.append(second)
        self._settle(second)
        for win in (first, second):
            self.assertIs(win.app_controller, controller)
            self.assertIs(win.run, run)
            self.assertIs(win.engine, controller.engine)
            self.assertIs(win._layout_store, controller.layouts)
        self.assertIsNot(first.ctx, second.ctx)  # independent windows (D7)
        self.assertEqual(controller.windows.explore_windows(), [first, second])
        self.assertEqual([controller.windows.number(w) for w in (first, second)], [1, 2])
        for name in (
            "info.open_sample", "info.open_annotation", "panel.open", "panel.reopen_closed",
            "panel.open_neighbourhood", "window.open_neighbourhood", "panel.open_alignment", "filter.set_subset",
            "filter.clear_subset", "nav.show_contig", "layout.reset", "app.run_list",
        ):
            self.assertTrue(controller.commands.has(name), name)

    def test_the_app_quits_after_its_last_explore_window_and_closes_the_cards(self) -> None:
        app = self.app
        before = app.quitOnLastWindowClosed()
        self.addCleanup(app.setQuitOnLastWindowClosed, before)
        with mock.patch.object(QtCore.QCoreApplication, "quit") as quit_app:
            controller = AppController(self._run_data(), manage_quit=True)
        self.assertFalse(app.quitOnLastWindowClosed())
        first = controller.open_first_window()
        second = controller.open_explore()
        self.windows += [first, second]
        self._settle(first)
        self._settle(second)
        card = first._open_sample_card_for_sample("S1")
        gene = first._open_gene_card(self._row(first))
        self.assertEqual(card.windowTitle(), "Sample card: S1 [1]")
        first.close()
        _spin(app, 0.1)
        quit_app.assert_not_called()
        # Orphaned cards rebind to the window that is left.
        self.assertIs(card.owner(), second)
        self.assertIs(gene.owner(), second)
        self.assertEqual(card.windowTitle(), "Sample card: S1")
        self.assertEqual(controller.windows.info_windows(second), [card, gene])
        second.close()
        _spin(app, 0.1)
        quit_app.assert_called_once_with()
        self.assertFalse(_alive(card))
        self.assertFalse(_alive(gene))

    def test_a_window_built_on_its_own_has_a_private_controller_without_quit_policy(self) -> None:
        win = self._window()
        controller = win.app_controller
        self.assertIsInstance(controller, AppController)
        self.assertIs(controller.run, win.run)
        self.assertIs(controller.engine, win.engine)
        self.assertEqual(controller.windows.explore_windows(), [win])
        card = win._open_sample_card_for_sample("S2")
        with mock.patch.object(QtCore.QCoreApplication, "quit") as quit_app:
            win.close()
            _spin(self.app, 0.1)
        quit_app.assert_not_called()
        self.assertFalse(_alive(card))

    def test_the_controller_refuses_another_runs_engine(self) -> None:
        run = self._run_data()
        other = self._run_data()
        controller = AppController(run)
        with self.assertRaises(ValueError):
            explore_window.ExploreWindow(run_data=other, app_controller=controller)
        with self.assertRaises(ValueError):
            AppController(run, AppController(other).engine)

    def test_no_floating_dialog_helper_is_left(self) -> None:
        self.assertFalse(hasattr(explore_window.ExploreWindow, "_show_floating_dialog"))
        self.assertFalse(hasattr(annotation_detail_dialog, "_show_floating_dialog"))
        win = self._window()
        self.assertFalse(hasattr(win, "_floating_dialogs"))


class InfoCardTests(_RealWindowCase):
    def test_opening_the_same_sample_card_twice_raises_it(self) -> None:
        win = self._window()
        first = win._open_sample_card_for_sample("S1")
        second = win._open_sample_card_for_sample("S1")
        # neoPan asks through the context: the same card.
        win.ctx.request(OpenSampleCard("S1"), origin="phylogeny")
        self.assertIs(first, second)
        self.assertEqual(win.app_controller.windows.info_windows(), [first])
        self.assertIsInstance(first, SampleCardDialog)
        self.assertIs(first.owner(), win)
        self.assertEqual(first.windowTitle(), "Sample card: S1")
        other = win._open_sample_card_for_sample("S2")
        self.assertIsNot(other, first)

    def test_the_sample_card_exposes_no_commands(self) -> None:
        win = self._window()
        card = win._open_sample_card_for_sample("S1")
        self.assertEqual(card.commands(), ())
        self.assertEqual(SampleCardDialog.COMMANDS, ())
        for name in win.app_controller.commands.names():
            with self.assertRaises(ValueError):
                card.run_command(name)
        tabs = card.findChild(QtWidgets.QTabWidget)
        self.assertEqual([tabs.tabText(i) for i in range(tabs.count())], ["Summary", "Raw"])
        # No buttons but Close: nothing that filters, reveals or opens.
        buttons = [b.text() for b in card.findChildren(QtWidgets.QAbstractButton) if not isinstance(b, QtWidgets.QTabBar)]
        self.assertEqual([t for t in buttons if t], ["Close"])
        # No hover output either: the card has no hover signal or callback.
        self.assertFalse(any(name for name in dir(card) if "hover" in name.lower()))

    def test_the_same_gene_from_the_track_and_the_table_is_one_card(self) -> None:
        win = self._window()
        row = self._row(win)
        win._on_igv_gene_double_clicked(dict(row))
        win._open_gene_card({k: (str(v) if k in ("Start", "Stop") else v) for k, v in row.items()})
        cards = win.app_controller.windows.info_windows(win)
        self.assertEqual(len(cards), 1)
        self.assertIsInstance(cards[0], AnnotationDetailDialog)
        # The Window menu tells cards of one gene in several samples apart.
        self.assertEqual(cards[0].menu_label(), f"{cards[0].windowTitle()} (S1, contig 1)")
        self.assertEqual(cards[0].objectName(), "AnnotationDetailDialog")
        self.assertEqual(cards[0].info_key(), InfoKey("gene", ("S1", "1", row["Locus"])))
        self.assertEqual(cards[0].commands(), (cmd.WINDOW_OPEN_NEIGHBOURHOOD, cmd.PANEL_OPEN_ALIGNMENT))

    def _other_windows(self, win: Any) -> List[Any]:
        others = [w for w in win.app_controller.windows.explore_windows() if w is not win]
        self.windows += [w for w in others if w not in self.windows]
        return others

    def test_compare_neighborhoods_opens_in_a_new_window_cloned_from_the_owner(self) -> None:
        win = self._window()
        win.gene_query.setText("blaOXA")
        win.refresh()
        _spin(self.app, 0.3)
        card = win._open_gene_card(self._row(win))
        self.assertTrue(card._context_button.isEnabled())
        self.assertIn("in a new explore window", card._context_button.toolTip())
        ran: List[str] = []
        original_run = win.app_controller.commands.run

        def recording_run(name: str, window: Any, **kwargs: Any) -> Any:
            ran.append(name)
            self.assertIs(window, win)
            return original_run(name, window, **kwargs)

        with mock.patch.object(win.app_controller.commands, "run", side_effect=recording_run):
            card._context_button.click()
        self.assertEqual(ran, [cmd.WINDOW_OPEN_NEIGHBOURHOOD])
        # A new explore window holds the comparison; the owner has none.
        others = self._other_windows(win)
        self.assertEqual(len(others), 1)
        new = others[0]
        self.assertEqual(_document_panels(win, "neighbourhood"), [])
        viewers = _document_panels(new, "neighbourhood")
        self.assertEqual(len(viewers), 1)
        viewer = viewers[0]
        self.assertIsInstance(viewer, NeighbourhoodPanel)
        self.assertTrue(new.isAncestorOf(viewer))
        self.assertTrue(new.isVisible())
        # Only the comparison in its docking area; the owner's filters applied.
        self.assertEqual(new._panel_manager.host.open_uids(), [viewer.uid])
        self.assertTrue(_wait(self.app, lambda: new.ctx.first_result_done))
        self.assertEqual(new.ctx.applied, win.ctx.draft)
        self.assertEqual(new.gene_query.text(), "blaOXA")
        # The filled panel (plan 8.7: an inline busy state, no progress dialog).
        self.assertFalse(card.findChildren(QtWidgets.QProgressDialog))
        self.assertTrue(_wait(self.app, lambda: viewer.view is not None))
        self.assertTrue(viewer.document.title.startswith("Gene neighborhoods: "))
        self.assertTrue(viewer.title().endswith(" neighborhoods"))  # the gene first, for narrow tabs
        self.assertEqual(new._panel_manager.title(viewer.uid), viewer.title())
        self.assertTrue(new.windowTitle().endswith(f" - {viewer.title()}"), new.windowTitle())
        # Asked again for the same gene, the card shows the open panel: no third window.
        card._context_button.click()
        self.assertEqual(len(self._other_windows(win)), 1)
        self.assertEqual(_document_panels(new, "neighbourhood"), [viewer])
        # The card stays with its owner; a gene double-clicked in the viewer
        # opens its card owned by the neighbourhood window.
        self.assertIs(card.owner(), win)
        viewer.view.geneDoubleClicked.emit(self._row(win, "S1:1:1"))
        gene_cards = new.app_controller.windows.info_windows(new)
        self.assertEqual(len(gene_cards), 1)
        # The viewer is part of its window: it goes with it.
        new.close()
        _spin(self.app, 0.1)
        self.assertFalse(_alive(viewer) and viewer.isVisible())
        self.assertTrue(win.isVisible())

    def test_align_shows_the_form_on_the_card_and_the_result_in_the_owner_window(self) -> None:
        win = self._window()
        card = win._open_gene_card(self._row(win))
        card._align_button.click()
        form = card._alignment_form
        self.assertIsInstance(form, GeneCardAlignmentDialog)
        self.assertIs(form.parent(), card)  # the parameter form is attached to the card
        card._align_button.click()
        self.assertIs(card._alignment_form, form)  # one form per card
        records = [("a|S1|1|g", "MKA"), ("b|S2|1|g", "MKA")]
        with mock.patch.object(GeneCardAlignmentDialog, "_collect_available_sequences", return_value=records), \
                mock.patch.object(GeneCardAlignmentDialog, "_perform_alignment", return_value=(">a|S1|1|g\nMKA\n>b|S2|1|g\nMKA\n", "Aligned.")):
            form._run_alignment()
        # The result is a panel of the owner window (plan 8.7).
        results = _document_panels(win, "alignment")
        self.assertEqual(len(results), 1)
        result = results[0]
        self.assertIsInstance(result, AlignmentPanel)
        self.assertTrue(result.title().startswith("Alignment: "))
        result.view._set_hover_sample("S2")
        win.ctx.flush_hover()
        self.assertEqual(win.ctx.hover_state.samples, ("S2",))
        self.assertEqual(win.ctx.hover_state.origin, result.uid)
        win._panel_manager.close_panel(result.uid)
        win.ctx.flush_hover()
        self.assertEqual(win.ctx.hover_state.samples, ())  # cleared on close
        card.close()
        _spin(self.app, 0.05)
        self.assertFalse(_alive(form))  # the form goes with its card

    def test_a_card_without_an_owner_has_its_actions_off(self) -> None:
        win = self._window()
        card = AnnotationDetailDialog(None, self._row(win), list(self._row(win)), ann_df=win.run.ann_df)
        self.addCleanup(card.deleteLater)
        self.assertFalse(card.can_run_commands())
        self.assertFalse(card._context_button.isEnabled())
        self.assertFalse(card._align_button.isEnabled())
        self.assertFalse(card.run_command(cmd.WINDOW_OPEN_NEIGHBOURHOOD, rows=[], title="x"))
        # Called directly, the actions start nothing whose result would be lost.
        with mock.patch.object(card, "_build_context_collection_rows") as collect:
            card._open_context_dialog()
        collect.assert_not_called()
        self.assertEqual(_document_panels(win, "neighbourhood"), [])
        card._open_alignment_dialog()
        self.assertIsNone(card._alignment_form)

    def test_plasmid_context_opens_through_the_command_in_a_new_window(self) -> None:
        win = self._window()
        card = AnnotationDetailDialog(None, self._row(win, "S1:1:1"), list(self._row(win)), ann_df=win.run.ann_df)
        rows = card._build_context_collection_rows(window_genes=30)
        card.deleteLater()
        self.assertTrue(rows)
        ran: List[Any] = []
        original_run = win.app_controller.commands.run

        def recording_run(name: str, window: Any, **kwargs: Any) -> Any:
            ran.append((name, window, dict(kwargs)))
            return original_run(name, window, **kwargs)

        with mock.patch.object(win.app_controller.commands, "run", side_effect=recording_run):
            win.ctx.request(OpenPlasmidContext.from_payload({"rows": rows, "name": "pX"}), origin="sample_tree")
        self.assertEqual(ran[0][0], cmd.WINDOW_OPEN_NEIGHBOURHOOD)
        self.assertIs(ran[0][1], win)
        self.assertEqual(ran[0][2]["source"], "plasmid")
        self.assertNotIn("beside", ran[0][2])  # a panel of this window does not apply there
        self.assertEqual(_document_panels(win, "neighbourhood"), [])
        others = self._other_windows(win)
        self.assertEqual(len(others), 1)
        new = others[0]
        viewers = _document_panels(new, "neighbourhood")
        self.assertEqual(len(viewers), 1)
        self.assertEqual(viewers[0].title(), "pX plasmid context")
        self.assertEqual(viewers[0].document.title, "Plasmid context: pX")
        self.assertTrue(new.isAncestorOf(viewers[0]))
        # The same plasmid again (from either window) shows the open panel.
        win.ctx.request(OpenPlasmidContext.from_payload({"rows": rows, "name": "pX"}), origin="sample_tree")
        new.ctx.request(OpenPlasmidContext.from_payload({"rows": rows, "name": "pX"}), origin=viewers[0].uid)
        self.assertEqual(len(self._other_windows(win)), 1)
        self.assertEqual(_document_panels(new, "neighbourhood"), viewers)
        # Another plasmid from the neighbourhood window opens another window,
        # cloned from that one.
        new.ctx.request(OpenPlasmidContext.from_payload({"rows": rows[:1], "name": "pY"}), origin=viewers[0].uid)
        self.assertEqual(len(self._other_windows(win)), 2)


class ShellMenuTests(_RealWindowCase):
    def test_file_run_list_is_the_only_run_action_and_returns_to_the_run_list(self) -> None:
        # Everything about runs (metadata, sharing, details, settings) lives in the
        # run list: the explore window keeps only the way back to it.
        opened: List[Any] = []
        controller = AppController(self._run_data(), open_run_list=lambda c, w: opened.append((c, w)))
        win = controller.open_first_window()
        self.windows.append(win)
        self._settle(win)
        titles = [a.text() for a in win.menuBar().actions()]
        self.assertEqual(titles[0], "File")
        self.assertNotIn("Edit", titles)
        self.assertEqual(titles[-1], "Window")
        texts = [a.text() for a in win._file_menu.actions() if a.text()]
        self.assertEqual(texts, ["Run List...", "Close Window", "Quit"])
        run_list = next(a for a in win._file_menu.actions() if a.objectName() == "RunListAction")
        self.assertTrue(run_list.shortcut().isEmpty())
        run_list.trigger()
        self.assertEqual(opened, [(controller, win)])

    def test_window_menu_lists_this_window_and_its_cards(self) -> None:
        win = self._window()
        win._open_sample_card_for_sample("S3")
        menu = win._window_menu
        # Filled once the window has registered: it lists itself before the
        # menu is ever shown.
        # Phase 8 added the Layouts submenu (plan 3.8, Window > Layouts); the
        # 2026-10-08 round New Panel... (Ctrl+Shift+P) and Show Filters, which
        # replace the top bar's Panels button and Filters toggle.
        head = ["New Explore Window", "Duplicate Window", "", "New Panel...", "Layouts", "Show Filters", ""]
        self.assertEqual([a.text() for a in menu.actions()][:9], head + ["Explore windows", win.windowTitle()])
        menu.aboutToShow.emit()
        texts = [a.text() for a in menu.actions()]
        self.assertEqual(texts, head + ["Explore windows", win.windowTitle(), "", "Info cards", "Sample card: S3"])
        shortcuts = {a.text(): a.shortcut().toString() for a in menu.actions() if not a.shortcut().isEmpty()}
        self.assertEqual(shortcuts, {"New Panel...": "Ctrl+Shift+P"})
        for action in win._file_menu.actions():
            self.assertTrue(action.shortcut().isEmpty())


class TwoWindowTests(_RealWindowCase):
    """Two explore windows of one AppController: shared engine, independent
    windows (D7), commands that run against one window only."""

    def _two(self, *, manage_quit: bool = False) -> Any:
        controller = AppController(self._run_data(), manage_quit=manage_quit)
        first = controller.open_first_window()
        second = controller.open_explore()
        self.windows += [first, second]
        self._settle(first)
        self._settle(second)
        return controller, first, second

    def test_commands_run_against_one_window_only(self) -> None:
        controller, first, second = self._two()
        run = controller.commands.run
        all_samples = tuple(first.ctx.sample_set())
        self.assertGreater(len(all_samples), 2)
        # filter.set_subset / filter.clear_subset
        run(cmd.FILTER_SET_SUBSET, second, label="Two samples", samples=["S1", "S2"])
        self.assertTrue(_wait(self.app, lambda: second.ctx.applied.subset is not None and not second.ctx.is_busy()))
        _spin(self.app, 0.1)
        self.assertEqual(second.ctx.applied.subset.samples, ("S1", "S2"))
        self.assertEqual(set(second.ctx.sample_set()), {"S1", "S2"})
        self.assertIsNone(first.ctx.draft.subset)
        self.assertEqual(tuple(first.ctx.sample_set()), all_samples)
        run(cmd.FILTER_CLEAR_SUBSET, second)
        self.assertTrue(_wait(self.app, lambda: second.ctx.applied.subset is None and not second.ctx.is_busy()))
        _spin(self.app, 0.1)
        self.assertEqual(tuple(second.ctx.sample_set()), all_samples)
        # nav.show_contig: "Open annotations" (Contig) by default, or Context.
        run(cmd.NAV_SHOW_CONTIG, second, sample="S1", contig="1")
        _spin(self.app, 0.2)
        self.assertEqual((second.ctx.selection.sample, second.ctx.selection.contig), ("S1", "1"))
        self.assertEqual(second._panel_manager.panel("annotations").view_mode(), "Contig")
        self.assertNotEqual(first.ctx.selection.sample, "S1")
        run(cmd.NAV_SHOW_CONTIG, second, sample="S1", contig="1", mode="Context")
        _spin(self.app, 0.2)
        self.assertEqual(second._panel_manager.panel("annotations").view_mode(), "Context")
        self.assertEqual(first._panel_manager.panel("annotations").view_mode(), "Contig")
        # panel.open, panel.reopen_closed, layout.reset
        manager = second._panel_manager
        uid = run(cmd.PANEL_OPEN, second, kind="annotations", beside="annotations")
        _spin(self.app, 0.2)
        self.assertEqual(manager.uids("annotations"), ["annotations", uid])
        self.assertEqual(first._panel_manager.uids("annotations"), ["annotations"])
        manager.host._docks[uid].requestCloseDockWidget()
        _spin(self.app, 0.1)
        self.assertEqual(manager.uids("annotations"), ["annotations"])
        self.assertEqual(run(cmd.PANEL_REOPEN_CLOSED, second), uid)
        _spin(self.app, 0.2)
        self.assertEqual(manager.uids("annotations"), ["annotations", uid])
        run(cmd.LAYOUT_RESET, second)
        _spin(self.app, 0.2)
        self.assertEqual(manager.uids("annotations"), ["annotations"])

    def test_hover_stays_in_its_window(self) -> None:
        _controller, first, second = self._two()
        second.ctx.hover(["S2"], origin="sample_tree")
        second.ctx.flush_hover()
        first.ctx.flush_hover()
        self.assertEqual(second.ctx.hover_state.samples, ("S2",))
        self.assertEqual(first.ctx.hover_state.samples, ())
        first.ctx.hover(["S3"], origin="sample_tree")
        first.ctx.flush_hover()
        self.assertEqual(second.ctx.hover_state.samples, ("S2",))
        self.assertEqual(first.ctx.hover_state.samples, ("S3",))

    def test_a_card_reopened_from_the_other_window_works_for_that_window(self) -> None:
        controller, first, second = self._two()
        card = first._open_sample_card_for_sample("S1")
        self.assertEqual(card.windowTitle(), "Sample card: S1 [1]")
        self.assertIs(second._open_sample_card_for_sample("S1"), card)
        self.assertIs(card.owner(), second)
        self.assertEqual(card.windowTitle(), "Sample card: S1 [2]")
        row = self._row(first)
        gene = first._open_gene_card(row)
        self.assertIs(second._open_gene_card(dict(row)), gene)
        self.assertIs(gene.owner(), second)
        self.assertEqual(gene.menu_label(), f"{gene.card_title()} (S1, contig 1) [2]")
        # Its actions now run on window 2: the comparison opens in a new
        # window cloned from window 2, not from window 1.
        second.gene_query.setText("blaOXA")
        gene._context_button.click()
        third = [w for w in controller.windows.explore_windows() if w not in (first, second)]
        self.windows += third
        self.assertEqual(len(third), 1)
        self.assertEqual(len(_document_panels(third[0], "neighbourhood")), 1)
        self.assertEqual(third[0].gene_query.text(), "blaOXA")
        self.assertEqual(_document_panels(first, "neighbourhood"), [])
        self.assertEqual(_document_panels(second, "neighbourhood"), [])
        gene._align_button.click()
        records = [("a|S1|1|g", "MKA"), ("b|S2|1|g", "MKA")]
        with mock.patch.object(GeneCardAlignmentDialog, "_collect_available_sequences", return_value=records), \
                mock.patch.object(GeneCardAlignmentDialog, "_perform_alignment", return_value=(">a|S1|1|g\nMKA\n>b|S2|1|g\nMKA\n", "Aligned.")):
            gene._alignment_form._run_alignment()
        self.assertEqual(_document_panels(first, "alignment"), [])
        result = _document_panels(second, "alignment")[0]
        result.view._set_hover_sample("S2")
        second.ctx.flush_hover()
        first.ctx.flush_hover()
        self.assertEqual(second.ctx.hover_state.samples, ("S2",))
        self.assertEqual(second.ctx.hover_state.origin, result.uid)
        self.assertEqual(first.ctx.hover_state.samples, ())

    def test_a_closed_explore_window_stops_autosaving_and_is_deleted(self) -> None:
        controller, first, second = self._two()
        timer = second._layout_autosave_timer
        self.assertTrue(timer.isActive())
        second.close()
        self.assertFalse(timer.isActive())
        QtCore.QCoreApplication.sendPostedEvents(None, QtCore.QEvent.Type.DeferredDelete.value)
        _spin(self.app, 0.05)
        self.assertFalse(_alive(second))
        self.assertEqual(controller.windows.explore_windows(), [first])

    def test_layouts_are_saved_before_the_app_quits(self) -> None:
        app = self.app
        before = app.quitOnLastWindowClosed()
        self.addCleanup(app.setQuitOnLastWindowClosed, before)
        with mock.patch.object(QtCore.QCoreApplication, "quit") as quit_app:
            controller = AppController(self._run_data(), manage_quit=True)
        win = controller.open_first_window()
        self.windows.append(win)
        self._settle(win)
        key = win._layout_run_key
        extra = controller.commands.run(cmd.PANEL_OPEN, win, kind="annotations", beside="annotations")
        _spin(app, 0.2)
        stored: List[Any] = []
        quit_app.side_effect = lambda *_a: stored.append(controller.layouts.load_window(key))
        win.close()
        _spin(app, 0.1)
        quit_app.assert_called_once_with()
        self.assertEqual(len(stored), 1)
        self.assertIsNotNone(stored[0])
        self.assertIn(extra, [entry["uid"] for entry in stored[0].panels])


class _Splash:
    def __init__(self) -> None:
        self.statuses: List[str] = []
        self.finished = 0

    def set_status(self, text: str) -> None:
        self.statuses.append(text)

    def finish(self) -> None:
        self.finished += 1


class AppLoaderTests(_RealWindowCase):
    """app.main's loaded handler: the AppController replaces state["win"]."""

    def test_the_loaded_run_gets_an_app_controller_that_quits_after_its_last_window(self) -> None:
        from annotation_explorer_polars import app as app_module

        app = self.app
        before = app.quitOnLastWindowClosed()
        self.addCleanup(app.setQuitOnLastWindowClosed, before)
        run = self._run_data()
        splash = _Splash()
        state: dict = {}
        with mock.patch.object(QtCore.QCoreApplication, "quit") as quit_app, \
                mock.patch.object(app_module, "_start_world_geojson_worker") as geo:
            app_module._open_loaded_run({"run_data": run, "warnings": []}, app=app, splash=splash, state=state)
            controller = state["app"]
            self.assertNotIn("win", state)
            self.assertIsInstance(controller, AppController)
            self.assertIs(controller.run, run)
            self.assertFalse(app.quitOnLastWindowClosed())
            windows = controller.windows.explore_windows()
            self.assertEqual(len(windows), 1)
            win = windows[0]
            self.windows.append(win)
            self.assertIs(win.app_controller, controller)
            self.assertTrue(win.isVisible())
            self.assertEqual(splash.statuses, ["Building views…"])
            # The splash goes when the first refresh settles.
            self.assertTrue(_wait(app, lambda: splash.finished > 0, 20))
            _spin(app, 0.05)
            geo.assert_called_once_with()
            card = win._open_sample_card_for_sample("S1")
            win.close()
            _spin(app, 0.1)
            quit_app.assert_called_once_with()
            self.assertFalse(_alive(card))

    def test_a_failed_build_closes_the_splash_and_quits(self) -> None:
        from annotation_explorer_polars import app as app_module

        splash = _Splash()
        state: dict = {}
        with mock.patch.object(QtCore.QCoreApplication, "exit") as exit_app, \
                mock.patch.object(QtWidgets.QMessageBox, "critical") as critical:
            app_module._open_loaded_run({}, app=self.app, splash=splash, state=state)
        self.assertEqual(splash.finished, 1)
        critical.assert_called_once()
        exit_app.assert_called_once_with(1)
        self.assertNotIn("app", state)


class ShortcutTests(_RealWindowCase):
    """D8: no new keyboard shortcut in the window, the cards or their forms,
    except the window's New Panel... (Ctrl+Shift+P, asked for by the user)."""

    @staticmethod
    def _shortcuts(widget: QtWidgets.QWidget) -> List[str]:
        keys = [s.key().toString() for s in widget.findChildren(QtGui.QShortcut)]
        actions = list(widget.findChildren(QtGui.QAction)) + list(widget.actions())
        keys += [a.shortcut().toString() for a in actions if not a.shortcut().isEmpty()]
        return sorted(set(keys))

    def test_no_window_card_or_form_has_a_new_shortcut(self) -> None:
        win = self._window()
        sample = win._open_sample_card_for_sample("S1")
        gene = win._open_gene_card(self._row(win))
        gene._align_button.click()
        self.assertEqual(self._shortcuts(win), ["Ctrl+Shift+P"])
        self.assertEqual(win._new_panel_action.shortcut().toString(), "Ctrl+Shift+P")
        for widget in (sample, gene, gene._alignment_form):
            self.assertEqual(self._shortcuts(widget), [], type(widget).__name__)
        # The neighbourhood viewer keeps its two shortcuts from before, now
        # acting only inside its panel (plan 3.7).
        gene._context_button.click()
        new = [w for w in win.app_controller.windows.explore_windows() if w is not win]
        self.windows += new
        self.assertEqual(self._shortcuts(new[0]), ["Ctrl+Shift+P"])
        viewer = _document_panels(new[0], "neighbourhood")[0]
        self.assertTrue(_wait(self.app, lambda: viewer.view is not None))
        self.assertEqual(self._shortcuts(viewer), ["Ctrl+F", "Ctrl+W"])
        for shortcut in viewer.findChildren(QtGui.QShortcut):
            self.assertEqual(shortcut.context(), QtCore.Qt.ShortcutContext.WidgetWithChildrenShortcut)


if __name__ == "__main__":
    unittest.main()
