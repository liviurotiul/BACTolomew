"""PanData: the pangenome outputs of one run (window plan 8.10, Phase 9).

One per run, on ``RunData`` as ``run.pan``, shared by every pan panel of every
explore window of the run (Phylogeny, Gene presence/absence, Metadata summary,
cgMLST MST (compact), Plasmidome; sharing read-only data is not linking). It says which pan
outputs the run has (the palette disables the panels whose data is missing)
and hands out what the pan views would otherwise each look up or build again:
the run's trees (``trees/trees.json``, or the old tree files), the Panaroo, fastGEAR and chewBBACA folders, the fastGEAR
clades and clusters, the per-sample hover lines and the run's frozen plasmid
cluster index.

No Qt. Folder lookups are done once, on first use, and cached (a run's output
does not change while it is open); the tables themselves are read through the
lru caches of ``panaroo_data``. ``memo`` caches values the pan views compute
from these (for example the tree's leaf order), so a second panel or window
reuses them. Thread safe.
"""
from __future__ import annotations

import json
import math
import threading
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any, Callable, Dict, Hashable, List, Mapping, Optional, Sequence, Tuple

import pandas as pd

from ..cgmlst_graph_data import _find_chewbbaca_out
from ..panaroo_data import (
    _find_fastgear_sample_clades,
    _find_panaroo_out,
    _load_fastgear_clusters,
    _load_fastgear_lineages,
)
from ..plasmid_model import PlasmidClusterIndex

# The pan views a standalone pan panel shows (``PanCoreAccessoryTab(single_view=...)``).
VIEW_TREE = "tree"
VIEW_PRESENCE = "presence"
VIEW_METADATA = "metadata"
VIEW_FORCE = "force"
VIEW_PLASMIDOME = "plasmidome"
SINGLE_VIEWS = (VIEW_TREE, VIEW_PRESENCE, VIEW_METADATA, VIEW_FORCE, VIEW_PLASMIDOME)

COUNTRY_TOOLTIP_CANDIDATES = ("Country", "_country")
ST_TOOLTIP_CANDIDATES = (
    "ST",
    "Sequence Type",
    "SequenceType",
    "sequence_type",
    "MLST ST",
    "MLST_ST",
)
COLLECTION_DATE_TOOLTIP_CANDIDATES = (
    "Collection date",
    "Collection Date",
    "Collection_date",
    "collection_date",
    "CollectionDate",
    "Date",
    "date",
)

_MISSING_TEXT = {"nan", "none", "null", "na", "n/a", "<na>", "missing", "unknown"}


def normalize_metadata_value(value: object) -> str:
    """A metadata cell as text; empty for missing values ("NA", "unknown", ...)."""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    text = str(value).strip()
    if not text or text.casefold() in _MISSING_TEXT:
        return ""
    return text


def resolve_column(columns: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    """The first candidate column present (exact name first, then ignoring case)."""
    available = {str(col): str(col) for col in columns}
    lowered = {str(col).casefold(): str(col) for col in columns}
    for candidate in candidates:
        if candidate in available:
            return available[candidate]
        match = lowered.get(str(candidate).casefold())
        if match:
            return match
    return None


def load_mlst_sequence_types(output_dir: Optional[Path]) -> Dict[str, str]:
    """Sample to MLST sequence type from ``mlst/results.csv`` (empty if absent)."""
    if output_dir is None:
        return {}
    path = Path(output_dir) / "mlst" / "results.csv"
    if not path.is_file():
        return {}
    try:
        df = pd.read_csv(path, dtype=str, low_memory=False)
    except Exception:
        return {}
    if not {"sample", "sequence_type"}.issubset({str(column).strip() for column in df.columns}):
        return {}
    sample_col = resolve_column(list(df.columns), ("sample",))
    st_col = resolve_column(list(df.columns), ("sequence_type",))
    if not sample_col or not st_col:
        return {}
    sequence_types: Dict[str, str] = {}
    for sample_raw, st_raw in df[[sample_col, st_col]].itertuples(index=False):
        sample = normalize_metadata_value(sample_raw)
        st_value = normalize_metadata_value(st_raw)
        if sample and st_value and sample not in sequence_types:
            sequence_types[sample] = st_value
    return sequence_types


def build_sample_hover_metadata(
    metadata_df: Optional[pd.DataFrame],
    output_dir: Optional[Path],
) -> Dict[str, Dict[str, str]]:
    """The lines a pan view shows when a sample is hovered: country, ST and
    collection date from the metadata, the ST from MLST when the metadata has none."""
    metadata: Dict[str, Dict[str, str]] = {}
    if metadata_df is not None and not metadata_df.empty and "Sample" in metadata_df.columns:
        columns = list(metadata_df.columns)
        country_col = resolve_column(columns, COUNTRY_TOOLTIP_CANDIDATES)
        st_col = resolve_column(columns, ST_TOOLTIP_CANDIDATES)
        date_col = resolve_column(columns, COLLECTION_DATE_TOOLTIP_CANDIDATES)
        for _idx, row in metadata_df.iterrows():
            sample = normalize_metadata_value(row.get("Sample"))
            if not sample:
                continue
            fields = metadata.setdefault(sample, {})
            for column, name in ((country_col, "Country"), (st_col, "ST"), (date_col, "Collection date")):
                if column and name not in fields:
                    value = normalize_metadata_value(row.get(column))
                    if value:
                        fields[name] = value
    for sample, sequence_type in load_mlst_sequence_types(output_dir).items():
        if sequence_type:
            metadata.setdefault(sample, {}).setdefault("ST", sequence_type)
    return {sample: fields for sample, fields in metadata.items() if fields}


# ---------------------------------------------------------------- trees
TREES_INDEX = Path("trees") / "trees.json"
TREES_INDEX_SCHEMA = 1

ROOTING_MIDPOINT = "midpoint"
ROOTING_OUTGROUP = "outgroup"

SUPPORT_NONE = "none"
SUPPORT_NAMES = {"ufboot": "UFBoot", "sh_alrt": "SH-aLRT"}

UNIT_SUBSTITUTIONS = "substitutions/site"
UNIT_GUBBINS_SNPS = "SNPs (reconstructed point mutations)"


@dataclass(frozen=True)
class BranchLengthSet:
    """One set of branch lengths of a tree (trees.json
    ``branch_length_sets``): the same topology, support and root, with the
    lengths of another fit (for the core tree: mutations only, the masked
    alignment, or including recombination, the unmasked one)."""

    id: str
    label: str
    path: Path
    default: bool = False
    branch_length_unit: Optional[str] = None
    alignment: Optional[str] = None
    model: Optional[str] = None


@dataclass(frozen=True)
class TreeSource:
    """One tree the pan views can draw: an entry of ``trees/trees.json``, or
    for a run without it a tree file at a known location (``rooted`` False:
    the file root is arbitrary, the views root it at the midpoint)."""

    id: str
    label: str
    path: Path
    branch_length_unit: Optional[str] = None
    # The root the file has: ``rooting`` ("midpoint" or "outgroup") and the
    # outgroup, from trees.json; ``rooted`` is False for a legacy file.
    rooted: bool = False
    rooting: Optional[str] = None
    outgroup: Optional[str] = None
    rooting_note: Optional[str] = None
    support_type: str = SUPPORT_NONE
    support_threshold: Optional[float] = None
    recombination_filtered: bool = False
    samples_missing: Tuple[str, ...] = ()
    # Genomes in the tree whose recombination was not masked (a skipped
    # lineage, dropped by Gubbins; trees.json schema 2).
    samples_unmasked: Tuple[str, ...] = ()
    method: Optional[str] = None
    model: Optional[str] = None
    notes: Tuple[str, ...] = ()
    # Its branch length sets, the default first; empty unless trees.json
    # lists two or more whose files exist (only then the Phylogeny toolbar
    # offers the choice).
    branch_length_sets: Tuple[BranchLengthSet, ...] = ()

    def branch_length_set(self, set_id: Optional[str]) -> Optional[BranchLengthSet]:
        """The set ``set_id``, else the default one (None without sets)."""
        for entry in self.branch_length_sets:
            if entry.id == set_id:
                return entry
        return self.branch_length_sets[0] if self.branch_length_sets else None

    def support_name(self) -> Optional[str]:
        """"UFBoot" or "SH-aLRT" when the file's internal labels are support
        values; None without support."""
        return SUPPORT_NAMES.get(self.support_type)

    @property
    def short_label(self) -> str:
        """A short name ("Core genome", "Snippy SNPs"), for a panel title."""
        for prefix, name in _TREE_SHORT_LABELS:
            if self.id.startswith(prefix):
                return name
        return self.label.split(",", 1)[0].strip() or self.id


# Short tree names by id prefix (pipeline/CoreTreeWorkflow.TREE_LABELS and
# the old file names of _LEGACY_TREES).
_TREE_SHORT_LABELS = (
    ("core", "Core genome"),
    ("gubbins", "Gubbins tree"),
    ("snippy", "Snippy SNPs"),
    ("panaroo", "Panaroo tree"),
)


@dataclass(frozen=True)
class TreeIndex:
    """The trees of a run (``tree_sources``): ``sources`` in display order,
    ``default_id`` the one every tree consumer draws unless a Phylogeny panel
    picks another, ``gubbins`` the Gubbins summary of trees.json (None
    without one) and ``from_index`` whether trees.json listed them."""

    sources: Tuple[TreeSource, ...] = ()
    default_id: Optional[str] = None
    gubbins: Optional[Mapping[str, Any]] = None
    from_index: bool = False

    def source(self, source_id: Optional[str]) -> Optional[TreeSource]:
        for source in self.sources:
            if source.id == source_id:
                return source
        return None

    def default(self) -> Optional[TreeSource]:
        return self.source(self.default_id)

    def gubbins_not_converged(self) -> bool:
        """Gubbins ran and did not report convergence (per lineage: a
        lineage that ran did not converge)."""
        gubbins = self.gubbins or {}
        if self.gubbins_per_lineage():
            return bool(self.gubbins_lineages_not_converged())
        return bool(gubbins.get("ran")) and gubbins.get("converged") is not True

    def gubbins_per_lineage(self) -> bool:
        """Gubbins ran once per fastGEAR lineage (trees.json schema 2)."""
        gubbins = self.gubbins or {}
        return gubbins.get("resolved_mode") == "per_lineage" and isinstance(gubbins.get("lineages"), (list, tuple))

    def gubbins_lineages(self) -> List[Mapping[str, Any]]:
        """The per-lineage Gubbins entries (empty for one run on all genomes)."""
        if not self.gubbins_per_lineage():
            return []
        return [entry for entry in (self.gubbins or {}).get("lineages", ()) if isinstance(entry, Mapping)]

    def gubbins_lineages_not_converged(self) -> List[Mapping[str, Any]]:
        """The lineages whose Gubbins run ran and did not converge."""
        return [entry for entry in self.gubbins_lineages() if entry.get("ran") and entry.get("converged") is not True]


def _text(value: object) -> Optional[str]:
    text = str(value).strip() if value is not None else ""
    return text or None


def _names_text(value: object) -> Optional[str]:
    """trees.json ``rooting.outgroup`` (a list of names, or null) as display text."""
    if isinstance(value, (list, tuple)):
        return ", ".join(str(item).strip() for item in value if str(item).strip()) or None
    return _text(value)


def _number(value: object) -> Optional[float]:
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _branch_length_sets(entry: Mapping[str, Any], trees_dir: Path, unit: Optional[str]) -> Tuple[BranchLengthSet, ...]:
    """A trees.json entry's ``branch_length_sets`` whose files exist, the
    default first; empty unless two or more remain (a missing file or an old
    run: no choice)."""
    listed = entry.get("branch_length_sets")
    if not isinstance(listed, list):
        return ()
    sets: List[BranchLengthSet] = []
    for item in listed:
        if not isinstance(item, Mapping):
            continue
        set_id, file_name = _text(item.get("id")), _text(item.get("file"))
        if not set_id or not file_name or any(known.id == set_id for known in sets):
            continue
        path = trees_dir / file_name
        if not path.is_file():
            continue
        sets.append(
            BranchLengthSet(
                id=set_id,
                label=_text(item.get("label")) or set_id,
                path=path,
                default=item.get("default") is True,
                branch_length_unit=_text(item.get("branch_length_unit")) or unit,
                alignment=_text(item.get("alignment")),
                model=_text(item.get("model")),
            )
        )
    if len(sets) < 2:
        return ()
    sets.sort(key=lambda item: not item.default)
    return tuple(sets)


def _index_source(entry: Mapping[str, Any], trees_dir: Path, top_rooting: Mapping[str, Any]) -> Optional[TreeSource]:
    """A trees.json entry as a TreeSource; None when it names no file that exists."""
    source_id = _text(entry.get("id"))
    file_name = _text(entry.get("file"))
    if not source_id or not file_name:
        return None
    path = trees_dir / file_name
    if not path.is_file():
        return None
    rooting = entry.get("rooting") if isinstance(entry.get("rooting"), Mapping) else top_rooting
    support = entry.get("support") if isinstance(entry.get("support"), Mapping) else {}
    support_type = _text(support.get("type")) or SUPPORT_NONE
    if support_type not in SUPPORT_NAMES:
        support_type = SUPPORT_NONE
    missing = entry.get("samples_missing")
    unmasked = entry.get("samples_unmasked")
    notes = entry.get("notes")
    unit = _text(entry.get("branch_length_unit"))
    return TreeSource(
        id=source_id,
        label=_text(entry.get("label")) or source_id,
        path=path,
        branch_length_unit=unit,
        rooted=True,
        rooting=_text(rooting.get("method")) or ROOTING_MIDPOINT,
        outgroup=_names_text(rooting.get("outgroup")),
        rooting_note=_text(rooting.get("note")),
        support_type=support_type,
        support_threshold=_number(support.get("threshold")) if support_type != SUPPORT_NONE else None,
        recombination_filtered=entry.get("recombination_filtered") is True,
        samples_missing=tuple(str(s) for s in missing) if isinstance(missing, list) else (),
        samples_unmasked=tuple(str(s) for s in unmasked) if isinstance(unmasked, list) else (),
        method=_text(entry.get("method")),
        model=_text(entry.get("model")),
        notes=tuple(str(n) for n in notes if str(n).strip()) if isinstance(notes, list) else (),
        branch_length_sets=_branch_length_sets(entry, trees_dir, unit),
    )


def _read_trees_index(output_dir: Path) -> Optional[TreeIndex]:
    """The trees of ``trees/trees.json`` (contract: PHYLOGENY_REWORK section
    4); None when the run has no readable index or it lists no tree file that
    exists (the caller then falls back to the old file names)."""
    path = output_dir / TREES_INDEX
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if not isinstance(data, Mapping) or not isinstance(data.get("trees"), list):
        return None
    top_rooting = data.get("rooting") if isinstance(data.get("rooting"), Mapping) else {}
    sources: List[TreeSource] = []
    for entry in data["trees"]:
        if not isinstance(entry, Mapping):
            continue
        source = _index_source(entry, path.parent, top_rooting)
        if source is not None and all(source.id != known.id for known in sources):
            sources.append(source)
    if not sources:
        return None
    # Default: ``default_tree`` when listed, else the first
    # recombination-filtered tree, else the first.
    default_id = _text(data.get("default_tree"))
    if default_id is None or all(source.id != default_id for source in sources):
        filtered = [source for source in sources if source.recombination_filtered]
        default_id = (filtered or sources)[0].id
    gubbins = data.get("gubbins")
    return TreeIndex(
        sources=tuple(sources),
        default_id=default_id,
        gubbins=MappingProxyType(dict(gubbins)) if isinstance(gubbins, Mapping) else None,
        from_index=True,
    )


# A run without trees.json: its tree files at their old names, in this order
# (the first recombination-filtered one is the default). The files are
# unrooted (the views root them at the midpoint).
_LEGACY_TREES = (
    (
        "gubbins_final",
        "Gubbins final tree, recombination-filtered",
        (Path("gubbins") / "core_gene.final_tree.tre",),
        UNIT_GUBBINS_SNPS,
        True,
    ),
    (
        "snippy_raxml",
        "Snippy SNPs, not recombination-corrected",
        (
            Path("trees") / "snippy_raxml.tree",
            Path("trees") / "snippy_raxml.newick",
            Path("snippy") / "core" / "snippy_raxml.tree",
        ),
        UNIT_SUBSTITUTIONS,
        False,
    ),
    (
        "core_gene_raxml",
        "Core genome, not recombination-corrected",
        (Path("trees") / "core_gene_tree.nwk",),
        UNIT_SUBSTITUTIONS,
        False,
    ),
)
PANAROO_TREE_ID = "panaroo"
PANAROO_TREE_LABEL = "Panaroo tree"


def _legacy_tree_index(output_dir: Path, panaroo_tree: Optional[Path]) -> TreeIndex:
    sources: List[TreeSource] = []
    for source_id, label, candidates, unit, filtered in _LEGACY_TREES:
        path = next((output_dir / c for c in candidates if (output_dir / c).is_file()), None)
        if path is not None:
            sources.append(
                TreeSource(id=source_id, label=label, path=path, branch_length_unit=unit, recombination_filtered=filtered)
            )
    if panaroo_tree is not None:
        sources.append(TreeSource(id=PANAROO_TREE_ID, label=PANAROO_TREE_LABEL, path=panaroo_tree))
    filtered_sources = [source for source in sources if source.recombination_filtered]
    default = (filtered_sources or sources)[0].id if sources else None
    return TreeIndex(sources=tuple(sources), default_id=default)


def tree_sources(output_dir: Optional[Path], panaroo_tree: Callable[[], Optional[Path]]) -> TreeIndex:
    """The run's trees: ``trees/trees.json`` when it lists any, else the tree
    files at the old names (Gubbins final tree, Snippy RAxML, core gene RAxML)
    and Panaroo's tree (``panaroo_tree()``, searched only then)."""
    if output_dir is None or not Path(output_dir).exists():
        return TreeIndex()
    output_dir = Path(output_dir)
    index = _read_trees_index(output_dir)
    if index is not None:
        return index
    return _legacy_tree_index(output_dir, panaroo_tree())


def find_panaroo_tree(output_dir: Optional[Path]) -> Optional[Path]:
    """The newest tree file in the Panaroo output, if any."""
    if output_dir is None:
        return None
    out = _find_panaroo_out(Path(output_dir))
    if out is None:
        return None
    patterns = ["*.nwk", "*.newick", "*tree*.nwk", "*tree*.newick", "*bestTree*"]
    matches: List[Path] = []
    for pattern in patterns:
        matches.extend([p for p in out.rglob(pattern) if p.is_file()])
    if not matches:
        return None
    matches.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return matches[0]


_UNSET = object()


class PanData:
    """The pangenome outputs of one run (immutable from the outside)."""

    def __init__(
        self,
        output_dir: Optional[Path],
        meta_df: Optional[pd.DataFrame] = None,
        plasmid_cluster_index: Optional[PlasmidClusterIndex] = None,
        reference_sample: Optional[str] = None,
    ) -> None:
        self._output_dir = Path(output_dir) if output_dir is not None else None
        self._meta_df = meta_df if meta_df is not None else pd.DataFrame()
        self._plasmid_cluster_index = (
            plasmid_cluster_index if plasmid_cluster_index is not None else PlasmidClusterIndex().freeze()
        )
        self._reference_sample = reference_sample
        self._lock = threading.RLock()
        self._cache: Dict[Hashable, Any] = {}

    # ------------------------------------------------------------- basics
    @property
    def output_dir(self) -> Optional[Path]:
        return self._output_dir

    @property
    def meta_df(self) -> pd.DataFrame:
        """The run's metadata table (shared; read only)."""
        return self._meta_df

    @property
    def reference_sample(self) -> Optional[str]:
        return self._reference_sample

    @property
    def plasmid_cluster_index(self) -> PlasmidClusterIndex:
        """The run's frozen MOB-suite plasmid cluster index (``RunData``'s)."""
        return self._plasmid_cluster_index

    def memo(self, key: Hashable, factory: Callable[[], Any]) -> Any:
        """``factory()`` once per ``key`` for the whole run (thread safe)."""
        with self._lock:
            value = self._cache.get(("memo", key), _UNSET)
            if value is _UNSET:
                value = factory()
                self._cache[("memo", key)] = value
            return value

    def _cached(self, name: str, factory: Callable[[], Any]) -> Any:
        with self._lock:
            value = self._cache.get(name, _UNSET)
            if value is _UNSET:
                value = factory() if self._output_dir is not None else None
                self._cache[name] = value
            return value

    # ------------------------------------------------------------- files
    def panaroo_tree_path(self) -> Optional[Path]:
        return self._cached("panaroo_tree", lambda: find_panaroo_tree(self._output_dir))

    def tree_index(self) -> TreeIndex:
        """The run's trees and the default one (``tree_sources``)."""
        value = self._cached("tree_index", lambda: tree_sources(self._output_dir, self.panaroo_tree_path))
        return value if value is not None else TreeIndex()

    def tree_path(self) -> Optional[Path]:
        """The default tree's file (every tree consumer draws it unless a
        Phylogeny panel picks another source)."""
        default = self.tree_index().default()
        return default.path if default is not None else None

    def panaroo_out(self) -> Optional[Path]:
        return self._cached("panaroo_out", lambda: _find_panaroo_out(self._output_dir))

    def fastgear_clades_path(self) -> Optional[Path]:
        return self._cached("fastgear", lambda: _find_fastgear_sample_clades(self._output_dir))

    def chewbbaca_out(self) -> Optional[Path]:
        return self._cached("chewbbaca", lambda: _find_chewbbaca_out(self._output_dir))

    def resolve_lookups(self) -> None:
        """Do the folder lookups now (they may walk the whole run folder).
        The loader calls this on its worker thread, so the palette's
        availability checks at window start do not search on the UI thread."""
        self.tree_path()
        self.panaroo_out()
        self.fastgear_clades_path()
        self.chewbbaca_out()

    # ------------------------------------------------------------- tables
    def fastgear_assignments(self) -> Tuple[Mapping[str, str], Mapping[str, str]]:
        """fastGEAR clade (lineage) and cluster of each sample (empty without fastGEAR)."""

        def load() -> Tuple[Mapping[str, str], Mapping[str, str]]:
            path = self.fastgear_clades_path()
            if path is None:
                return MappingProxyType({}), MappingProxyType({})
            try:
                lineages = _load_fastgear_lineages(path)
            except Exception:
                lineages = {}
            try:
                clusters = _load_fastgear_clusters(path)
            except Exception:
                clusters = {}
            return MappingProxyType(dict(lineages)), MappingProxyType(dict(clusters))

        with self._lock:
            value = self._cache.get("fastgear_assignments", _UNSET)
            if value is _UNSET:
                value = load()
                self._cache["fastgear_assignments"] = value
            return value

    def sample_hover_metadata(self) -> Dict[str, Dict[str, str]]:
        """Sample to its hover lines (a fresh copy for the caller)."""
        with self._lock:
            value = self._cache.get("hover_metadata", _UNSET)
            if value is _UNSET:
                value = build_sample_hover_metadata(self._meta_df, self._output_dir)
                self._cache["hover_metadata"] = value
        return {sample: dict(fields) for sample, fields in value.items()}

    # ------------------------------------------------------------- availability
    def has_metadata(self) -> bool:
        df = self._meta_df
        return bool(not df.empty and "Sample" in df.columns and any(str(c) != "Sample" for c in df.columns))

    def availability(self, view: str) -> Tuple[bool, str]:
        """Whether the run has what a standalone pan panel of ``view`` shows,
        and if not, why (the palette's tooltip)."""
        if self._output_dir is None:
            return False, "Not available: this run has no output folder."
        if view == VIEW_TREE:
            if self.tree_path() is None:
                return False, "Not available: this run has no phylogenetic tree."
        elif view == VIEW_PRESENCE:
            if self.panaroo_out() is None:
                return False, "Not available: this run has no Panaroo output."
        elif view == VIEW_METADATA:
            if not self.has_metadata():
                return False, "Not available: this run has no sample metadata."
        elif view == VIEW_FORCE:
            if self.chewbbaca_out() is None:
                return False, "Not available: this run has no chewBBACA (cgMLST) output."
        elif view == VIEW_PLASMIDOME:
            if not self._plasmid_cluster_index:
                return False, "Not available: this run has no MOB-suite plasmid clusters."
        else:
            return False, "Unknown pan view"
        return True, ""
