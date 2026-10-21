"""Built-in docking layouts (plan 3.8): Classic and the presets.

Every preset arranges the docking area only (the rail is outside the layout)
and never touches the filters or the sample set.

Classic (the default) reproduces the pre-docking window, with the tree and
the gene presence/absence heatmap side by side below:

    +---------------------+-----------------------------+
    | Sample tree         | Annotations (table)         |
    |                     +-----------------------------+
    |                     | Genome track                |
    +----------------+----+------------+----------------+
    | Map            | Phylogeny       | Gene presence/ |
    |                |                 | absence        |
    +----------------+-----------------+----------------+

The cgMLST MST is not open by default.

Genome focus:              Pangenome:                       Geo + table:

    +------+-------------+     +----------+----------------+     +---------+-------------+
    |Sample| Annotations |     | Phylogeny| Gene presence/ |     | Map     | Annotations |
    |tree  +-------------+     |          | absence        |     |         |             |
    |      | Genome track|     +----------+-------+--------+     |         |             |
    +------+-------------+     | Map              | Plasmid|     +---------+-------------+
                               |                  | -ome   |
                               +------------------+--------+

The Genome track is placed only with Annotations: it cannot exist without
it (``PanelSpec.requires``, enforced by ``PanelManager``). Phylogeny and Gene
presence/absence sit side by side at the same height, so at the tree's fit
the window's row link (panels/row_link.py) puts each tree leaf level with its
heatmap row (the heatmap always fits every row).

A preset leaves out the panels a run does not have; it is offered only when
at least one of its main panels is available.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

from ..panels.base import ANNOTATIONS, GENE_PRESENCE_ABSENCE, GEO, PHYLOGENY, PLASMIDOME, SAMPLE_TREE, TRACK
from .dock_host import DockHost, Placement

CLASSIC_PANELS = (SAMPLE_TREE, ANNOTATIONS, TRACK, GEO, PHYLOGENY, GENE_PRESENCE_ABSENCE)

# Proportions of the Classic splitters. The bottom row gets 30 percent of the
# height, more than the old fixed layout's fifth, so the pan panels are
# usable without resizing first (Phase 0 follow-up). The genome track keeps
# about the height it had before (about 210 px of content at 1800x1050): with
# a tab strip above both the table and the track it needs a larger share
# than the old 75:25 table split.
CLASSIC_ROWS = (70, 30)  # top row : bottom row
CLASSIC_BOTTOM_COLUMNS = (34, 30, 36)  # map : phylogeny : gene presence/absence
CLASSIC_TOP_COLUMNS = (42, 58)  # sample tree : annotations column
CLASSIC_TABLE_TRACK = (62, 38)  # annotations : genome track

GENOME_FOCUS_COLUMNS = (26, 74)  # sample tree : annotations column
PANGENOME_ROWS = (62, 38)  # phylogeny and heatmap row : map and plasmidome row
PANGENOME_COLUMNS = (35, 65)  # phylogeny : gene presence/absence
PANGENOME_BOTTOM_COLUMNS = (50, 50)  # map : plasmidome
GEO_TABLE_COLUMNS = (50, 50)  # map : annotations
# Open document panels (gene neighborhoods and the like) kept by a preset or
# a saved layout sit in a column on the right with this percent of the width
# (about what Classic gives them), whatever the layout's own splits.
DOCUMENT_COLUMN_SHARE = 33

PRESET_CLASSIC = "classic"
PRESET_GENOME_FOCUS = "genome_focus"
PRESET_PANGENOME = "pangenome"
PRESET_GEO_TABLE = "geo_table"


@dataclass(frozen=True)
class Preset:
    """One built-in layout of the Window > Layouts menu."""

    key: str
    title: str
    description: str
    # It is offered when at least one of these is available for the run.
    needs_any: Tuple[str, ...]
    build: Callable[[DockHost, Sequence[str]], None]


def _present(host: DockHost, panels: Iterable[str], available: Iterable[str]) -> List[str]:
    allowed = set(available)
    return [uid for uid in panels if uid in allowed and host.has_view(uid)]


def build_classic(host: DockHost, available: Iterable[str] = CLASSIC_PANELS) -> None:
    """Place the Classic panels into an empty docking area."""
    present = _present(host, CLASSIC_PANELS, available)

    def has(uid: str) -> bool:
        return uid in present

    if has(SAMPLE_TREE):
        host.place_view(SAMPLE_TREE, Placement("left", outer=True))
    if has(ANNOTATIONS):
        if has(SAMPLE_TREE):
            host.place_view(ANNOTATIONS, Placement("right", relative_to=SAMPLE_TREE))
        else:
            host.place_view(ANNOTATIONS, Placement("right", outer=True))
    if has(TRACK) and has(ANNOTATIONS):
        # The track never exists without Annotations (PanelSpec.requires).
        host.place_view(TRACK, Placement("bottom", relative_to=ANNOTATIONS))
    bottom = [uid for uid in (GEO, PHYLOGENY, GENE_PRESENCE_ABSENCE) if has(uid)]
    if bottom:
        top_row = any(has(uid) for uid in (SAMPLE_TREE, ANNOTATIONS))
        host.place_view(bottom[0], Placement("bottom", outer=True) if top_row else Placement("left", outer=True))
        for left, uid in zip(bottom, bottom[1:]):
            host.place_view(uid, Placement("right", relative_to=left))
        if len(bottom) == len(CLASSIC_BOTTOM_COLUMNS):
            host.set_split_sizes(bottom[0], list(CLASSIC_BOTTOM_COLUMNS))
        elif len(bottom) == 2:
            host.set_split_sizes(bottom[0], [50, 50])
        if top_row:
            # The rows: the bottom groups sit in a row of their own when there
            # are several of them (one splitter further out).
            host.set_split_sizes(bottom[0], list(CLASSIC_ROWS), outer=1 if len(bottom) > 1 else 0)
    if has(SAMPLE_TREE) and has(ANNOTATIONS):
        host.set_split_sizes(SAMPLE_TREE, list(CLASSIC_TOP_COLUMNS))
    if has(TRACK) and has(ANNOTATIONS):
        host.set_split_sizes(TRACK, list(CLASSIC_TABLE_TRACK))


def build_genome_focus(host: DockHost, available: Iterable[str]) -> None:
    """A narrow sample tree beside the annotation table over the genome track."""
    present = _present(host, (SAMPLE_TREE, ANNOTATIONS, TRACK), available)
    if SAMPLE_TREE in present:
        host.place_view(SAMPLE_TREE, Placement("left", outer=True))
    if ANNOTATIONS in present:
        # The track goes only with Annotations (PanelSpec.requires).
        host.place_view(ANNOTATIONS, Placement("right", outer=True))
        if TRACK in present:
            host.place_view(TRACK, Placement("bottom", relative_to=ANNOTATIONS))
            host.set_split_sizes(TRACK, list(CLASSIC_TABLE_TRACK))
        if SAMPLE_TREE in present:
            host.set_split_sizes(SAMPLE_TREE, list(GENOME_FOCUS_COLUMNS))


def build_pangenome(host: DockHost, available: Iterable[str]) -> None:
    """Phylogeny beside Gene presence/absence (linked rows, same height), the
    map beside Plasmidome below them."""
    present = _present(host, (PHYLOGENY, GENE_PRESENCE_ABSENCE, GEO, PLASMIDOME), available)
    top = [uid for uid in (PHYLOGENY, GENE_PRESENCE_ABSENCE) if uid in present]
    bottom = [uid for uid in (GEO, PLASMIDOME) if uid in present]
    if top:
        host.place_view(top[0], Placement("left", outer=True))
        if len(top) == 2:
            host.place_view(top[1], Placement("right", relative_to=top[0]))
            host.set_split_sizes(top[0], list(PANGENOME_COLUMNS))
    if bottom:
        host.place_view(bottom[0], Placement("bottom", outer=True) if top else Placement("left", outer=True))
        if len(bottom) == 2:
            host.place_view(bottom[1], Placement("right", relative_to=bottom[0]))
            host.set_split_sizes(bottom[0], list(PANGENOME_BOTTOM_COLUMNS))
    if top and bottom:
        # The two rows (one splitter further out when the bottom row holds
        # two groups side by side).
        host.set_split_sizes(bottom[0], list(PANGENOME_ROWS), outer=1 if len(bottom) > 1 else 0)


def build_geo_table(host: DockHost, available: Iterable[str]) -> None:
    """The map and the annotation table side by side."""
    present = _present(host, (GEO, ANNOTATIONS), available)
    if GEO in present:
        host.place_view(GEO, Placement("left", outer=True))
    if ANNOTATIONS in present:
        host.place_view(ANNOTATIONS, Placement("right", outer=True))
        if GEO in present:
            # ADS keeps the sizes of the splitter it reuses; split evenly.
            host.set_split_sizes(GEO, list(GEO_TABLE_COLUMNS))


PRESETS: Tuple[Preset, ...] = (
    Preset(
        PRESET_CLASSIC,
        "Classic",
        "The default: sample tree, annotations and genome track, with the map, the phylogeny and gene presence/absence below",
        CLASSIC_PANELS,
        build_classic,
    ),
    Preset(
        PRESET_GENOME_FOCUS,
        "Genome focus",
        "A narrow sample tree beside the annotation table and the genome track",
        (SAMPLE_TREE, ANNOTATIONS),
        build_genome_focus,
    ),
    Preset(
        PRESET_PANGENOME,
        "Pangenome",
        "Phylogeny beside gene presence/absence, with the map and Plasmidome below",
        (PHYLOGENY, GENE_PRESENCE_ABSENCE, PLASMIDOME),
        build_pangenome,
    ),
    Preset(
        PRESET_GEO_TABLE,
        "Geo + table",
        "The map and the annotation table side by side",
        (GEO, ANNOTATIONS),
        build_geo_table,
    ),
)

_PRESETS_BY_KEY: Dict[str, Preset] = {preset.key: preset for preset in PRESETS}


def preset(key: str) -> Optional[Preset]:
    return _PRESETS_BY_KEY.get(key)
