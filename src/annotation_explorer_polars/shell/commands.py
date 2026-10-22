"""Named commands (plan 4.2, ``CommandRegistry``).

Every cross-cutting action is a named command. A command always runs against
one explore window: the caller's, or for a card its owner's. Menus, buttons
and cards run commands by name, so a keyboard binding can be added later
without restructuring (D8: none is added, except Ctrl+Shift+P on Window >
New Panel..., which runs ``panel.open`` without a kind: the panel picker).

Phase 5 registers the commands whose actions existed then; Phase 7 adds the
window commands (``window.new``, ``window.duplicate``, ``window.open_subset``,
``window.close``, ``app.quit``; later ``window.open_neighbourhood``); Phase 8 the layout commands (``layout.save``,
``layout.apply``, ``layout.set_default``, ``layout.clear_default``,
``layout.delete``).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Mapping, Optional, Sequence

if TYPE_CHECKING:  # pragma: no cover
    from .app_controller import AppController

INFO_OPEN_SAMPLE = "info.open_sample"
INFO_OPEN_ANNOTATION = "info.open_annotation"
PANEL_OPEN = "panel.open"
PANEL_REOPEN_CLOSED = "panel.reopen_closed"
PANEL_OPEN_NEIGHBOURHOOD = "panel.open_neighbourhood"
PANEL_OPEN_ALIGNMENT = "panel.open_alignment"
FILTER_SET_SUBSET = "filter.set_subset"
FILTER_CLEAR_SUBSET = "filter.clear_subset"
NAV_SHOW_CONTIG = "nav.show_contig"
LAYOUT_RESET = "layout.reset"
LAYOUT_SAVE = "layout.save"
LAYOUT_APPLY = "layout.apply"
LAYOUT_SET_DEFAULT = "layout.set_default"
LAYOUT_CLEAR_DEFAULT = "layout.clear_default"
LAYOUT_DELETE = "layout.delete"
APP_RUN_LIST = "app.run_list"
WINDOW_NEW = "window.new"
WINDOW_DUPLICATE = "window.duplicate"
WINDOW_OPEN_SUBSET = "window.open_subset"
WINDOW_OPEN_NEIGHBOURHOOD = "window.open_neighbourhood"
WINDOW_CLOSE = "window.close"
APP_QUIT = "app.quit"

Handler = Callable[..., Any]


@dataclass(frozen=True)
class Command:
    name: str
    title: str
    handler: Handler


class CommandRegistry:
    """Name to handler. ``run(name, window, **kwargs)`` calls
    ``handler(window, **kwargs)`` and returns its result."""

    def __init__(self) -> None:
        self._commands: Dict[str, Command] = {}

    def register(self, name: str, handler: Handler, title: str = "") -> None:
        if name in self._commands:
            raise ValueError(f"Command {name!r} is registered already")
        self._commands[name] = Command(name, title or name, handler)

    def has(self, name: str) -> bool:
        return name in self._commands

    def names(self) -> List[str]:
        return list(self._commands)

    def title(self, name: str) -> str:
        return self._commands[name].title

    def run(self, name: str, window: Any, **kwargs: Any) -> Any:
        command = self._commands.get(name)
        if command is None:
            raise KeyError(f"Unknown command {name!r}")
        if window is None:
            raise ValueError(f"Command {name!r} needs an explore window")
        return command.handler(window, **kwargs)


def register_explore_commands(registry: CommandRegistry, controller: "AppController") -> None:
    """The commands of an explore window. The handlers carry out the action in
    the given window (most delegate to the window's existing code); the info
    card commands go through the controller's WindowManager."""
    windows = controller.windows

    def open_sample(window: Any, *, sample: str) -> Any:
        from ..sample_card_dialog import SampleCardDialog
        from ..windows.info_window import sample_key

        sample = str(sample or "").strip()
        if not sample:
            return None
        return windows.open_info(
            sample_key(sample),
            lambda: SampleCardDialog(None, **window.engine.sample_bundle(sample)),
            owner=window,
        )

    def open_annotation(window: Any, *, row: Mapping[str, Any]) -> Any:
        from ..annotation_detail_dialog import AnnotationDetailDialog
        from ..windows.info_window import gene_key

        row_data = dict(row or {})
        if not row_data:
            return None
        return windows.open_info(
            gene_key(row_data),
            lambda: AnnotationDetailDialog(None, row_data, list(row_data.keys()), **window.engine.annotation_bundle()),
            owner=window,
        )

    def open_panel(
        window: Any,
        *,
        kind: Optional[str] = None,
        params: Optional[Mapping[str, Any]] = None,
        beside: Optional[str] = None,
    ) -> Any:
        if not kind:
            window._show_panel_palette()
            return None
        from ..session.requests import OpenPanel, params_tuple

        return window._open_panel_request(OpenPanel(kind, params_tuple(dict(params or {})), beside))

    def reopen_closed(window: Any) -> Any:
        manager = getattr(window, "_panel_manager", None)
        return manager.reopen_closed() if manager is not None else None

    def open_neighbourhood(window: Any, **kwargs: Any) -> Any:
        # In this window (the neighbourhood window places its panel with it).
        return window._open_neighbourhood_view(**kwargs)

    def open_alignment(window: Any, **kwargs: Any) -> Any:
        return window._open_alignment_view(**kwargs)

    def set_subset(
        window: Any,
        *,
        label: str,
        samples: Sequence[str],
        kind: str = "sample",
        origin: str = "rail",
    ) -> None:
        window._set_sample_subset_filter(label=label, samples=list(samples), kind=kind, origin=origin)

    def clear_subset(window: Any) -> None:
        window._clear_sample_subset_filter()

    def show_contig(window: Any, *, sample: str, contig: str, mode: str = "Contig", **payload: Any) -> None:
        # The tree's "Open annotations" (mode "Contig") and "Open in Context view" (mode "Context").
        window._open_tree_contig_payload({"sample": sample, "contig": contig, **payload}, mode=mode)

    def reset_layout(window: Any) -> None:
        manager = getattr(window, "_panel_manager", None)
        if manager is not None:
            manager.reset_layout()

    # Layouts (plan 3.8): the docking area only, never the rail or filters.
    # ``layout`` is a user layout's name (``name`` is taken by ``run``).
    def save_layout(window: Any, *, layout: Optional[str] = None) -> Optional[str]:
        # Without a name the window asks for one (Save Layout As...).
        saver = getattr(window, "_save_layout_as", None)
        return saver(layout) if callable(saver) else None

    def apply_layout(window: Any, *, preset: Optional[str] = None, layout: Optional[str] = None) -> bool:
        applier = getattr(window, "_apply_layout", None)
        return bool(applier(preset=preset, name=layout)) if callable(applier) else False

    def set_default_layout(window: Any) -> bool:
        setter = getattr(window, "_set_default_layout", None)
        return bool(setter()) if callable(setter) else False

    def clear_default_layout(window: Any) -> None:
        del window  # the default is shared by every window
        controller.layouts.clear_default()

    def delete_layout(window: Any, *, layout: str) -> bool:
        del window  # user layouts are shared by every window
        return controller.layouts.delete_user_layout(layout)

    def run_list(window: Any) -> None:
        # Close this run's windows and choose a run in the run list.
        controller.open_run_list(window)

    def new_window(window: Any) -> Any:
        # A clean rail and the default layout (plan 3.10).
        return controller.open_explore(source=window)

    def duplicate_window(window: Any) -> Any:
        # A copy of this window's filters and layout, once; independent after.
        return controller.duplicate_explore(window)

    def open_subset_window(window: Any, *, subset: Any) -> Any:
        # A new window whose rail carries the subset as its locked window subset.
        return controller.open_explore(source=window, window_subset=subset)

    def open_neighbourhood_window(window: Any, **kwargs: Any) -> Any:
        # Compare neighborhoods (a gene card, a plasmid context): a new explore
        # window with this window's filters, holding the comparison panel.
        return controller.open_neighbourhood_window(window, **kwargs)

    def close_window(window: Any) -> bool:
        return bool(window.close())

    def quit_app(window: Any) -> None:
        del window  # every explore window of the run closes
        controller.quit_all()

    registry.register(INFO_OPEN_SAMPLE, open_sample, "Open sample card")
    registry.register(INFO_OPEN_ANNOTATION, open_annotation, "Open gene card")
    registry.register(PANEL_OPEN, open_panel, "New panel")
    registry.register(PANEL_REOPEN_CLOSED, reopen_closed, "Reopen closed panel")
    registry.register(PANEL_OPEN_NEIGHBOURHOOD, open_neighbourhood, "Compare neighborhoods in this window")
    registry.register(PANEL_OPEN_ALIGNMENT, open_alignment, "Show alignment")
    registry.register(FILTER_SET_SUBSET, set_subset, "Filter to samples")
    registry.register(FILTER_CLEAR_SUBSET, clear_subset, "Clear sample subset")
    registry.register(NAV_SHOW_CONTIG, show_contig, "Show contig in Annotations")
    registry.register(LAYOUT_RESET, reset_layout, "Reset layout")
    registry.register(LAYOUT_SAVE, save_layout, "Save layout as")
    registry.register(LAYOUT_APPLY, apply_layout, "Apply layout")
    registry.register(LAYOUT_SET_DEFAULT, set_default_layout, "Set as default for new windows")
    registry.register(LAYOUT_CLEAR_DEFAULT, clear_default_layout, "Use Classic for new windows")
    registry.register(LAYOUT_DELETE, delete_layout, "Delete saved layout")
    registry.register(APP_RUN_LIST, run_list, "Run list")
    registry.register(WINDOW_NEW, new_window, "New explore window")
    registry.register(WINDOW_DUPLICATE, duplicate_window, "Duplicate window")
    registry.register(WINDOW_OPEN_SUBSET, open_subset_window, "Open in new explore window")
    registry.register(WINDOW_OPEN_NEIGHBOURHOOD, open_neighbourhood_window, "Compare neighborhoods")
    registry.register(WINDOW_CLOSE, close_window, "Close window")
    registry.register(APP_QUIT, quit_app, "Quit")
