"""Gubbins per lineage and the recombination-masked core alignment (step ``gubbins``).

Input: ``core_alignment/core_gene_alignment.ordered.aln`` (``pipeline.CoreAlignment``)
and the lineages of the population structure step (``fastgear/sample_clades.tsv`` from
fastGEAR, or ``fastbaps/sample_clades.tsv`` from FastBAPS, ``pipeline.PopulationStructure``).
The same rule applies to every dataset (``phylo.gubbins_mode``):

- ``auto`` (default): one Gubbins run per lineage when there are two or more,
  else one run on all genomes; ``single``: always one run; ``per_lineage``: per
  lineage, an error when the lineages are unavailable.
- A lineage with fewer than ``phylo.gubbins_min_lineage_size`` genomes, or fewer than
  three variable columns (Gubbins' ``--min-snps``), is not run; its genomes stay
  unmasked and the summary says why.
- Each run reads the lineage's rows of the ordered alignment with every column, so
  Gubbins coordinates are ordered-alignment coordinates. Runs go in parallel within the
  step's thread budget, each into a temp folder that replaces ``lineage_<id>/`` only on
  success; a lineage whose genomes, input and parameters did not change is reused.
- Every lineage is masked with Gubbins' own ``mask_gubbins_aln.py`` (``N``), then the
  masked rows are merged into ``gubbins/core_gene_alignment.masked.aln`` in the ordered
  alignment's row order, each genome masked only by its own lineage's calls (checked
  row by row). The core tree is built on that file.
- Recombination between lineages (imports on the long inter-lineage branch) is not
  masked; this is recorded in the summary's notes.

Single mode keeps Gubbins' outputs at the top of ``gubbins/`` (``core_gene.*``) as before.
``gubbins/summary.json`` (schema 3) records the mode, every group with its genomes,
convergence, recombination statistics and why a group did not run, and the masked
alignment. Warnings never stop the run; a Gubbins failure does.
"""
from __future__ import annotations

import csv
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import tempfile
import threading
import time
from concurrent.futures import FIRST_EXCEPTION, ThreadPoolExecutor, wait
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator, Mapping, Sequence

import numpy as np

from pipeline.CoreAlignment import (
    BASE_CODE,
    ORDERED_ALIGNMENT_REL,
    fasta_record,
    file_sha256,
    iter_fasta_records,
    map_alignment_interval,
    predicted_gubbins_dropouts,
    read_gene_map,
    read_ordering_summary,
)
from pipeline.Newick import leaf_names, read_newick
from pipeline.PopulationStructure import sample_clades_path
from pipeline.ToolEnv import phylo_outgroup_names, tool_cmd

ALGORITHM_VERSION = 3
SUMMARY_SCHEMA_VERSION = 3
GUBBINS_PREFIX = "core_gene"
GUBBINS_SUMMARY = "summary.json"
GUBBINS_FINAL_TREE = f"{GUBBINS_PREFIX}.final_tree.tre"
GUBBINS_NODE_TREE = f"{GUBBINS_PREFIX}.node_labelled.final_tree.tre"
GUBBINS_FILTERED_ALIGNMENT = f"{GUBBINS_PREFIX}.filtered_polymorphic_sites.fasta"
GUBBINS_RECOMBINATION_GFF = f"{GUBBINS_PREFIX}.recombination_predictions.gff"
GUBBINS_PER_BRANCH = f"{GUBBINS_PREFIX}.per_branch_statistics.csv"
BLOCK_GENES_TSV = "recombination_blocks_genes.tsv"
COMBINED_GFF = "recombination_predictions.gff"
MASKED_ALIGNMENT = "core_gene_alignment.masked.aln"
MASKED_ALIGNMENT_REL = Path("gubbins") / MASKED_ALIGNMENT
MASKED_BASES_TSV = "masked_bases.tsv"
GROUP_STAMP = "run.json"
SINGLE_GROUP = "all"
MASK_CHAR = "N"
# Gubbins' --min-snps default: with fewer variable columns no block can be called.
MIN_VARIABLE_SITES = 3
# Files of one Gubbins run, in its folder (lineage_<id>/, or gubbins/ in single mode).
GROUP_OUT_FILES: tuple[str, ...] = (
    GUBBINS_FINAL_TREE,
    GUBBINS_FILTERED_ALIGNMENT,
    GUBBINS_RECOMBINATION_GFF,
    GUBBINS_PER_BRANCH,
)
# Files of a schema 2 run (one cohort run, before per-lineage runs).
GUBBINS_OUT_FILES: tuple[Path, ...] = (
    Path(GUBBINS_FINAL_TREE),
    Path(GUBBINS_FILTERED_ALIGNMENT),
    Path(GUBBINS_RECOMBINATION_GFF),
    Path(GUBBINS_SUMMARY),
)

# Why a group did not run.
SKIP_FEW_GENOMES = "below_min_genomes"
SKIP_LINEAGE_SIZE = "below_min_lineage_size"
SKIP_VARIABLE_SITES = "too_few_variable_sites"
# Why a genome of the masked alignment is not masked.
UNMASKED_LINEAGE_SKIPPED = "lineage_skipped"
UNMASKED_DROPPED = "dropped_by_gubbins"
UNMASKED_NO_LINEAGE = "no_lineage"

INTER_LINEAGE_NOTE = (
    "Recombination between lineages is not masked: each lineage is analyzed on its own, so "
    "imports from one lineage into another on the long branch between them (or before they "
    "split) are invisible. They lengthen the branches between lineages, not the branching "
    "within a lineage."
)

_CONVERGED = re.compile(r"Convergence after (\d+) iterations?: (.+?)\s*$", re.MULTILINE)
_MAX_REACHED = re.compile(r"Maximum number of iterations \((\d+)\) reached", re.MULTILINE)


@dataclass(frozen=True)
class GubbinsResult:
    output_dir: Path
    ran: bool
    summary: dict[str, Any]
    warnings: list[dict[str, str]] = field(default_factory=list)


@dataclass
class GubbinsGroup:
    """One Gubbins run: a lineage (fastGEAR or FastBAPS), or ``all`` genomes in single mode."""

    id: str
    members: list[str]
    skip_reason: str | None = None

    @property
    def key(self) -> tuple[str, ...]:
        return tuple(sorted(self.members))


@dataclass
class GroupPlan:
    requested_mode: str
    resolved_mode: str  # single | per_lineage
    mode_reason: str  # configured | two_or_more_lineages | one_lineage | fastgear_unavailable
    groups: list[GubbinsGroup]
    unassigned: list[str] = field(default_factory=list)
    warnings: list[dict[str, str]] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    source: str = "fastGEAR"  # the method that found the lineages (fastGEAR or FastBAPS)

    @property
    def single(self) -> bool:
        return self.resolved_mode == "single"


def _warning(code: str, message: str) -> dict[str, str]:
    return {"step": "gubbins", "code": code, "message": message}


def _group_label(group_id: str) -> str:
    return "the cohort" if group_id == SINGLE_GROUP else f"lineage {group_id}"


def _lineage_sort_key(value: str) -> tuple[int, Any]:
    return (0, int(value)) if value.isdigit() else (1, value)


# ------------------------------------------------------------ mode and groups
def read_fastgear_lineages(output_path: Path) -> dict[str, str]:
    """Sample -> lineage from the run's sample_clades.tsv (fastgear/ or fastbaps/, by the
    run's population method; raises when unreadable or when the run has no method)."""
    path = sample_clades_path(output_path)
    if path is None:
        raise RuntimeError("the run has no population structure method")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "Sample" not in reader.fieldnames or "FastGearLineage" not in reader.fieldnames:
            raise RuntimeError(f"{path} has no Sample and FastGearLineage columns.")
        return {
            str(row["Sample"]).strip(): str(row["FastGearLineage"] or "").strip()
            for row in reader
            if str(row.get("Sample") or "").strip()
        }


def resolve_gubbins_groups(
    samples: Sequence[str],
    settings: Mapping[str, Any],
    lineages: Mapping[str, str] | None,
    lineage_error: str | None = None,
    *,
    lineage_source: str = "fastGEAR",
) -> GroupPlan:
    """The Gubbins groups of a cohort by phylo.gubbins_mode (see the module docstring).

    lineages: sample -> lineage, or None when the population structure method
    (lineage_source: fastGEAR or FastBAPS) is unavailable or failed (lineage_error says
    why). The variable-column rule is applied later, on the staged alignments.
    """
    samples = list(samples)
    requested = str(settings["gubbins_mode"])
    min_genomes = int(settings["gubbins_min_genomes"])
    min_size = int(settings["gubbins_min_lineage_size"])
    warnings: list[dict[str, str]] = []
    notes: list[str] = []

    def single(reason: str) -> GroupPlan:
        skip = SKIP_FEW_GENOMES if len(samples) < min_genomes else None
        return GroupPlan(requested, "single", reason, [GubbinsGroup(SINGLE_GROUP, samples, skip)], [], warnings, notes,
                         lineage_source)

    if requested == "single":
        return single("configured")
    if lineages is None:
        detail = lineage_error or "fastgear/sample_clades.tsv is missing"
        if requested == "per_lineage":
            raise RuntimeError(
                f"phylo.gubbins_mode is per_lineage, but the {lineage_source} lineages are unavailable ({detail})."
            )
        warnings.append(
            _warning(
                "gubbins_lineages_unavailable",
                f"The {lineage_source} lineages are unavailable ({detail}), so Gubbins runs once on all genomes. "
                "On a cohort of distinct lineages its recombination calls are less reliable; a resume "
                f"that completes {lineage_source} reruns Gubbins per lineage.",
            )
        )
        return single("fastgear_unavailable")
    by_lineage: dict[str, list[str]] = {}
    unassigned: list[str] = []
    for sample in samples:
        lineage = str(lineages.get(sample) or "").strip()
        if lineage:
            by_lineage.setdefault(lineage, []).append(sample)
        else:
            unassigned.append(sample)
    if len(by_lineage) < 2:
        if requested == "per_lineage":
            notes.append(
                f"phylo.gubbins_mode is per_lineage, but {lineage_source} found one lineage, so Gubbins ran once "
                "on all genomes."
            )
        return single("one_lineage")
    groups = [
        GubbinsGroup(lineage, members, SKIP_LINEAGE_SIZE if len(members) < min_size else None)
        for lineage, members in sorted(by_lineage.items(), key=lambda item: _lineage_sort_key(item[0]))
    ]
    reason = "two_or_more_lineages" if requested == "auto" else "configured"
    return GroupPlan(requested, "per_lineage", reason, groups, unassigned, warnings, notes, lineage_source)


def group_folder_rel(group_id: str, single: bool) -> Path:
    """A group's output folder relative to the run (gubbins/ itself in single mode)."""
    return Path("gubbins") if single else Path("gubbins") / f"lineage_{group_id}"


def expected_gubbins_files(summary: Mapping[str, Any] | None) -> list[Path]:
    """Run-relative output files of a Gubbins step from its summary (any schema)."""
    if summary is None or summary.get("schema_version") != SUMMARY_SCHEMA_VERSION:
        if summary is not None and summary.get("ran") is False:
            return [Path("gubbins") / GUBBINS_SUMMARY]
        return [Path("gubbins") / rel for rel in GUBBINS_OUT_FILES]
    single = summary.get("resolved_mode") == "single"
    ran = [str(g.get("id")) for g in summary.get("groups") or [] if isinstance(g, Mapping) and g.get("ran")]
    return planned_gubbins_files(ran, single)


def planned_gubbins_files(ran_group_ids: Sequence[str], single: bool) -> list[Path]:
    files = [Path("gubbins") / GUBBINS_SUMMARY]
    if not ran_group_ids:
        return files
    files += [MASKED_ALIGNMENT_REL, Path("gubbins") / MASKED_BASES_TSV, Path("gubbins") / BLOCK_GENES_TSV]
    if not single:
        files.append(Path("gubbins") / COMBINED_GFF)
    for group_id in ran_group_ids:
        folder = group_folder_rel(group_id, single)
        files += [folder / name for name in (*GROUP_OUT_FILES, GROUP_STAMP)]
    return files


# ------------------------------------------------------------ Gubbins output parsing
def parse_convergence(stdout: str) -> dict[str, Any]:
    """Convergence from Gubbins' stdout: converged True / False / None (not stated)."""
    found = _CONVERGED.search(stdout or "")
    if found:
        return {
            "converged": True,
            "iterations_run": int(found.group(1)),
            "convergence_reason": found.group(2).rstrip("."),
        }
    found = _MAX_REACHED.search(stdout or "")
    if found:
        return {
            "converged": False,
            "iterations_run": int(found.group(1)),
            "convergence_reason": "maximum number of iterations reached",
        }
    return {"converged": None, "iterations_run": None, "convergence_reason": None}


def not_converged_warning(
    convergence: Mapping[str, Any], max_iterations: int, group_id: str = SINGLE_GROUP
) -> dict[str, str] | None:
    if convergence.get("converged") is True:
        return None
    if convergence.get("converged") is False:
        detail = f"it stopped at the maximum of {max_iterations} iterations"
    else:
        detail = "its log does not state convergence"
    where = "" if group_id == SINGLE_GROUP else f" on lineage {group_id}"
    return _warning(
        "gubbins_not_converged",
        f"Gubbins did not converge{where} ({detail}): the recombination calls and the "
        "recombination-filtered tree come from the last iteration and may still change "
        "with more iterations (phylo.gubbins_iterations).",
    )


def gubbins_version() -> str | None:
    try:
        result = subprocess.run(
            [*tool_cmd("run_gubbins.py"), "--version"], capture_output=True, text=True, check=False
        )
    except FileNotFoundError:
        return None
    for line in reversed(((result.stdout or "") + "\n" + (result.stderr or "")).splitlines()):
        if re.fullmatch(r"\d+(\.\d+)+", line.strip()):
            return line.strip()
    return None


def gubbins_params(settings: Mapping[str, Any], *, threads: int, outgroup: str | None) -> dict[str, Any]:
    return {
        "seed": int(settings["seed"]),
        "iterations": int(settings["gubbins_iterations"]),
        "converge_method": str(settings["gubbins_converge_method"]),
        "filter_percentage": float(settings["gubbins_filter_percentage"]),
        "tree_builder": str(settings["gubbins_tree_builder"]),
        "outgroup": outgroup,
        "threads": int(max(1, threads)),
    }


def _without_threads(params: Any) -> Any:
    """Parameters compared for reuse: the thread count is recorded but does not make a
    finished run stale (a resume on another machine or core budget reuses it)."""
    if not isinstance(params, Mapping):
        return params
    return {key: value for key, value in params.items() if key != "threads"}


def _predicted_dropouts(ordering: Mapping[str, Any], settings: Mapping[str, Any]) -> list[str]:
    """Taxa Gubbins will leave out with the current phylo.gubbins_filter_percentage (the
    gap fraction of a row is the same in a lineage's input, which keeps every column)."""
    missing = ordering.get("missing_fraction")
    if isinstance(missing, Mapping):
        return predicted_gubbins_dropouts(
            {str(k): float(v) for k, v in missing.items()}, float(settings["gubbins_filter_percentage"])
        )
    return list(ordering.get("predicted_gubbins_dropouts") or [])


def _resolve_outgroup(settings: Mapping[str, Any], dropouts: Sequence[str]) -> tuple[str | None, list[dict[str, str]], list[str]]:
    """(Gubbins --outgroup value, warnings, notes); single mode only."""
    names = phylo_outgroup_names(settings)
    if not names:
        return None, [], []
    dropped = [name for name in names if name in set(dropouts)]
    if dropped:
        message = (
            f"Outgroup sample(s) {', '.join(dropped)} have too many N or gap characters and will be "
            "left out by Gubbins, so Gubbins runs without an outgroup."
        )
        return None, [_warning("outgroup_unusable", message)], []
    notes = []
    if len(names) > 1:
        notes.append(
            f"Gubbins takes one outgroup name; it was given {names[0]}, the first of the configured "
            "outgroup. The trees are rooted on the whole outgroup afterwards."
        )
    return names[0], [], notes


def _read_per_branch(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def root_edge_nodes(node_labelled_tree: Path) -> list[str]:
    """Names of the root's two children in Gubbins' node-labelled tree ([] otherwise).

    Gubbins roots its tree (midpoint or outgroup), which splits one unrooted edge into
    the root's two branches; the lineage-structure share counts them as one edge, so
    it does not depend on where the root falls on that edge.
    """
    try:
        root = read_newick(node_labelled_tree)
    except Exception:
        return []
    if len(root.children) != 2 or not all(child.name for child in root.children):
        return []
    return [child.name for child in root.children]


def recombination_statistics(per_branch: Path, taxa: set[str], root_edge: Sequence[str] = ()) -> dict[str, Any]:
    rows = _read_per_branch(per_branch)
    inside = sum(_float(r.get("Number of SNPs Inside Recombinations")) for r in rows)
    outside = sum(_float(r.get("Number of SNPs Outside Recombinations")) for r in rows)
    # One row per unrooted edge: the root's two branches are merged into one.
    edges: list[tuple[str, float]] = []
    root_snps = 0.0
    for r in rows:
        snps = _float(r.get("Number of SNPs Outside Recombinations"))
        if root_edge and r.get("Node") in root_edge:
            root_snps += snps
        else:
            edges.append((str(r.get("Node") or ""), snps))
    if root_edge:
        edges.append((" + ".join(root_edge) + " (root edge)", root_snps))
    longest = max(edges, key=lambda edge: edge[1], default=None)
    share = None
    if longest is not None and outside > 0:
        share = longest[1] / outside
    tips = []
    for r in rows:
        if r.get("Node") in taxa:
            length = _float(r.get("Genome Length"))
            frac = _float(r.get("Bases in Recombinations Excluding Gaps")) / length if length else 0.0
            tips.append((frac, r["Node"]))
    tips.sort(reverse=True)
    return {
        "snps_inside_recombination": int(inside),
        "snps_outside_recombination": int(outside),
        "r_m_snp_ratio": (inside / outside) if outside else None,
        "longest_branch_share": {
            "node": longest[0] if longest else None,
            "share": share,
            "snps_outside_recombination": int(longest[1]) if longest else None,
        },
        "tip_recombination_fraction_top5": [{"sample": name, "fraction": round(frac, 5)} for frac, name in tips[:5]],
    }


def gff_blocks(path: Path) -> list[tuple[int, int, list[str]]]:
    """(start, end, taxa) of every recombination block of a Gubbins GFF (1-based, inclusive)."""
    blocks = []
    with Path(path).open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            taxa = re.search(r'taxa="([^"]*)"', cols[8])
            names = [x for x in (taxa.group(1).split() if taxa else []) if x]
            blocks.append((int(cols[3]), int(cols[4]), names))
    return blocks


def write_block_gene_table(gffs: Sequence[tuple[str, Path]], gene_map: Path, dest: Path) -> dict[str, dict[str, int]]:
    """One row per recombination block (lineage, coordinates, genes it overlaps) from
    each group's GFF, ordered alignment coordinates mapped back through gene_map.tsv.

    Genes the ordering genome lacks are appended in Panaroo order, so their neighbors
    in the alignment are arbitrary: a block touching one is flagged
    (touches_appended_gene) and left out of blocks_spanning_several_genes. The
    genome_span column (same contig only) shows how far apart in the ordering genome
    the first and last gene are; intergenic and accessory DNA is absent from the
    alignment, so a short block can bridge genes that are far apart.
    """
    rows = read_gene_map(gene_map)
    stats: dict[str, dict[str, int]] = {}
    tmp = dest.with_name(dest.name + ".tmp")
    with tmp.open("w", encoding="utf-8") as out:
        out.write(
            "lineage\tblock_start\tblock_end\tlength\ttaxa_count\tgenes_overlapped\tgenes\t"
            "first_gene_contig\tfirst_gene_genome_start\tlast_gene_contig\tlast_gene_genome_end\t"
            "genome_span\ttouches_appended_gene\tany_gene_outside_panaroo_filtered\n"
        )
        for group_id, gff in gffs:
            blocks = gff_blocks(gff)
            multi_gene = 0
            touching_appended = 0
            for start, end, taxa in blocks:
                genes = map_alignment_interval(rows, start, end)
                appended = any(not g.get("placed") for g in genes)
                touching_appended += int(appended)
                multi_gene += int(len(genes) > 1 and not appended)
                first, last = (genes[0], genes[-1]) if genes else ({}, {})
                span = ""
                if (
                    not appended
                    and first.get("ordering_contig")
                    and first.get("ordering_contig") == last.get("ordering_contig")
                    and first.get("ordering_start") is not None
                    and last.get("ordering_end") is not None
                ):
                    span = abs(int(last["ordering_end"]) - int(first["ordering_start"])) + 1
                outside_filtered = any(g.get("in_panaroo_filtered") is False for g in genes)
                values = (
                    group_id,
                    start,
                    end,
                    end - start + 1,
                    len(taxa),
                    len(genes),
                    ";".join(g["cluster"] for g in genes),
                    first.get("ordering_contig") or "",
                    first.get("ordering_start") or "",
                    last.get("ordering_contig") or "",
                    last.get("ordering_end") or "",
                    span,
                    "true" if appended else "false",
                    "true" if outside_filtered else "false",
                )
                out.write("\t".join(str(v) for v in values) + "\n")
            stats[group_id] = {
                "blocks": len(blocks),
                "blocks_spanning_several_genes": multi_gene,
                "blocks_touching_appended_genes": touching_appended,
            }
    os.replace(tmp, dest)
    return stats


def write_combined_gff(gffs: Sequence[tuple[str, Path]], dest: Path) -> int:
    """Every group's recombination blocks in one GFF, column 9 extended by lineage="<id>"."""
    count = 0
    tmp = dest.with_name(dest.name + ".tmp")
    with tmp.open("w", encoding="utf-8") as out:
        out.write("##gff-version 3\n")
        for group_id, gff in gffs:
            with Path(gff).open("r", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    if line.startswith("#") or not line.strip():
                        continue
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 9:
                        continue
                    attrs = cols[8].rstrip()
                    cols[8] = f'{attrs}{"" if attrs.endswith(";") else ";"}lineage="{group_id}";'
                    out.write("\t".join(cols) + "\n")
                    count += 1
    os.replace(tmp, dest)
    return count


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    os.replace(tmp, path)


def _format_subprocess_error(stdout: str, stderr: str, returncode: int) -> str:
    text = "\n".join(part for part in ((stdout or "").strip(), (stderr or "").strip()) if part)
    return "\n".join(text.splitlines()[-30:]) or f"exit code {returncode}"


def diversity_warnings(diversity: Mapping[str, Any], settings: Mapping[str, Any], group_id: str = SINGLE_GROUP) -> list[dict[str, str]]:
    mean = float(diversity.get("mean_pairwise_p") or 0.0)
    maximum = float(diversity.get("max_pairwise_p") or 0.0)
    mean_limit = float(settings["gubbins_diversity_warn"])
    max_limit = float(settings["gubbins_max_diversity_warn"])
    if mean <= mean_limit and maximum <= max_limit:
        return []
    who = "The cohort's core genome is" if group_id == SINGLE_GROUP else f"The core genome of lineage {group_id} is"
    advice = (
        " Consider running it per lineage (phylo.gubbins_mode)." if group_id == SINGLE_GROUP else ""
    )
    return [
        _warning(
            "gubbins_high_diversity",
            f"{who} diverse for Gubbins (mean pairwise p-distance {mean:.4f}, max {maximum:.4f}; "
            f"warning above {mean_limit} / {max_limit}). Gubbins is designed for closely related "
            f"isolates: on distant genomes recombination calls and the filtered tree are less reliable.{advice}",
        )
    ]


# ------------------------------------------------------------ staging
@dataclass
class StagedGroup:
    """A group's Gubbins input: its rows of the ordered alignment, every column."""

    path: Path
    sha256: str
    variable_sites: int
    diversity: dict[str, Any]


class _GroupCounter:
    """Per-column state of one group while its rows stream by: the first called base,
    whether a second state was seen, and the called bases (capped at 2)."""

    def __init__(self, length: int) -> None:
        self.first = np.full(length, 4, dtype=np.uint8)
        self.variable = np.zeros(length, dtype=bool)
        self.called = np.zeros(length, dtype=np.uint8)

    def add(self, codes: np.ndarray) -> None:
        called = codes < 4
        np.copyto(self.first, codes, where=called & (self.first == 4))
        self.variable |= called & (codes != self.first)
        self.called += called
        np.minimum(self.called, 2, out=self.called)


def _pairwise_diversity(path: Path, variable_idx: np.ndarray, columns_compared: int) -> dict[str, Any]:
    """Mean and max pairwise p-distance of a staged group (same definition as the
    cohort's in core_alignment/ordering.json: over columns with >= 2 called bases)."""
    names: list[str] = []
    rows: list[np.ndarray] = []
    for name, raw in iter_fasta_records(path):
        rows.append(BASE_CODE[np.frombuffer(raw, dtype=np.uint8)][variable_idx])
        names.append(name)
    result: dict[str, Any] = {
        "columns_compared": int(columns_compared),
        "variable_sites": int(variable_idx.size),
        "mean_pairwise_p": 0.0,
        "max_pairwise_p": 0.0,
        "max_pair": None,
        "median_pairwise_snps": 0.0,
    }
    if not variable_idx.size or len(rows) < 2 or not columns_compared:
        return result
    matrix = np.vstack(rows)
    counts = np.stack([(matrix == k).sum(axis=0) for k in range(4)]).astype(np.float64)
    called = counts.sum(axis=0)
    usable = called >= 2
    same = (counts[:, usable] * (counts[:, usable] - 1.0)).sum(axis=0) / (called[usable] * (called[usable] - 1.0))
    result["mean_pairwise_p"] = float((1.0 - same).sum() / float(columns_compared))
    n = matrix.shape[0]
    both = np.zeros((n, n), dtype=np.float64)
    same_pairs = np.zeros((n, n), dtype=np.float64)
    for lo in range(0, matrix.shape[1], 20000):
        block = matrix[:, lo: lo + 20000]
        called_block = (block < 4).astype(np.float32)
        both += called_block @ called_block.T
        for k in range(4):
            onehot = (block == k).astype(np.float32)
            same_pairs += onehot @ onehot.T
    iu = np.triu_indices(n, k=1)
    pair_diffs = (both - same_pairs)[iu]
    best = int(np.argmax(pair_diffs))
    result["max_pairwise_p"] = float(pair_diffs[best]) / float(columns_compared)
    result["max_pair"] = [names[iu[0][best]], names[iu[1][best]]]
    result["median_pairwise_snps"] = float(np.median(pair_diffs))
    return result


def stage_group_alignments(
    alignment: Path, groups: Sequence[GubbinsGroup], dest_dir: Path
) -> dict[str, StagedGroup]:
    """Write each group's rows of the ordered alignment (all columns, the alignment's
    row order and line format) to dest_dir/<id>.aln in one pass, with its sha256 and
    variable columns; then the group's pairwise diversity."""
    group_of = {name: group.id for group in groups for name in group.members}
    handles: dict[str, Any] = {}
    digests: dict[str, Any] = {}
    counters: dict[str, _GroupCounter] = {}
    paths = {group.id: Path(dest_dir) / f"{group.id}.aln" for group in groups}
    try:
        for group in groups:
            handles[group.id] = paths[group.id].open("wb")
            digests[group.id] = hashlib.sha256()
        for name, raw in iter_fasta_records(alignment):
            group_id = group_of.get(name)
            if group_id is None:
                continue
            seq = np.frombuffer(raw, dtype=np.uint8)
            counter = counters.get(group_id)
            if counter is None:
                counter = counters[group_id] = _GroupCounter(seq.size)
            counter.add(BASE_CODE[seq])
            chunk = fasta_record(name, seq)
            digests[group_id].update(chunk)
            handles[group_id].write(chunk)
    finally:
        for handle in handles.values():
            handle.close()
    staged: dict[str, StagedGroup] = {}
    for group in groups:
        counter = counters.get(group.id)
        if counter is None:
            raise RuntimeError(f"No row of {_group_label(group.id)} is in {alignment}.")
        variable_idx = np.flatnonzero(counter.variable)
        compared = int(np.count_nonzero(counter.called >= 2))
        staged[group.id] = StagedGroup(
            path=paths[group.id],
            sha256=digests[group.id].hexdigest(),
            variable_sites=int(variable_idx.size),
            diversity=_pairwise_diversity(paths[group.id], variable_idx, compared),
        )
    return staged


def split_threads(total: int, weights: Mapping[str, int]) -> dict[str, int]:
    """Threads per group: 1 each, the rest of the budget in proportion to the group's
    variable columns, the remainder to the largest group. The sum is the budget when
    there are at most that many groups; with more groups each gets 1 and
    ``_run_parallel`` queues them so the running ones stay within the budget."""
    if not weights:
        return {}
    total = max(1, int(total))
    extra = max(0, total - len(weights))
    weight_sum = sum(max(1, int(w)) for w in weights.values())
    threads = {key: 1 + (extra * max(1, int(w))) // weight_sum for key, w in weights.items()}
    largest = max(weights, key=lambda key: (int(weights[key]), key))
    threads[largest] += max(0, total - sum(threads.values()))
    return threads


# ------------------------------------------------------------ masking
@dataclass(frozen=True)
class MaskSource:
    """One group's masked alignment for the merge: rows of all its members in the
    ordered alignment's order, masked by its own GFF over the taxa Gubbins kept."""

    group_id: str
    members: frozenset[str]
    kept: frozenset[str]
    masked_alignment: Path
    gff: Path


def _blocks_by_taxon(source: MaskSource) -> dict[str, list[tuple[int, int]]]:
    by_taxon: dict[str, list[tuple[int, int]]] = {}
    for start, end, taxa in gff_blocks(source.gff):
        outside = sorted(set(taxa) - source.members)
        if outside:
            raise RuntimeError(
                f"A recombination block of {_group_label(source.group_id)} ({start}..{end}) names taxa "
                f"outside it: {', '.join(outside[:5])}."
            )
        for taxon in taxa:
            by_taxon.setdefault(taxon, []).append((start, end))
    return by_taxon


def merge_masked_alignment(
    alignment: Path,
    sources: Sequence[MaskSource],
    unmasked_reasons: Mapping[str, tuple[str | None, str]],
    dest: Path,
) -> dict[str, Any]:
    """Stream the ordered alignment once and write dest: a genome of a group that ran
    and that Gubbins kept takes its row from that group's masked alignment, every other
    genome keeps its original row (unmasked_reasons: name -> (group, reason)).

    Each masked row is checked against the original: it differs only where it now holds
    N, and exactly at that genome's own blocks in its group's GFF (minus columns already
    N). So a genome is masked only by its own lineage's calls; any mismatch raises.

    Per genome: masked_bases counts masked A/C/G/T (masked_positions also counts the
    gaps Gubbins' masking turns into N) and masked_fraction is masked_bases over the
    genome's called bases. A block on an internal branch masks every genome below it,
    so a high fraction means sitting below such a branch, not recombining more.
    """
    mask_byte = ord(MASK_CHAR)
    source_of = {name: source for source in sources for name in source.members}
    blocks = {source.group_id: _blocks_by_taxon(source) for source in sources}
    iterators: dict[str, Iterator[tuple[str, bytes]]] = {
        source.group_id: iter_fasta_records(source.masked_alignment) for source in sources
    }
    digest = hashlib.sha256()
    rows: list[dict[str, Any]] = []
    length = None
    tmp = dest.with_name(dest.name + ".tmp")
    try:
        with tmp.open("wb") as out:
            for name, raw in iter_fasta_records(alignment):
                original = np.frombuffer(raw, dtype=np.uint8)
                if length is None:
                    length = original.size
                elif original.size != length:
                    raise RuntimeError(f"Row {name} of {alignment} has {original.size} columns, not {length}.")
                source = source_of.get(name)
                row = original
                group_id: str | None = None
                reason: str | None = None
                masked_count = 0
                masked_positions = 0
                if source is not None:
                    group_id = source.group_id
                    masked_name, masked_raw = next(iterators[group_id], (None, b""))
                    if masked_name != name:
                        raise RuntimeError(
                            f"The masked alignment of {_group_label(group_id)} has {masked_name!r} where "
                            f"{name!r} was expected."
                        )
                    masked = np.frombuffer(masked_raw.upper(), dtype=np.uint8)
                    if masked.size != original.size:
                        raise RuntimeError(f"The masked row of {name} has {masked.size} columns, not {original.size}.")
                    expected = np.zeros(original.size, dtype=bool)
                    if name in source.kept:
                        for start, end in blocks[group_id].get(name, ()):
                            expected[start - 1: end] = True
                    expected &= original != mask_byte
                    changed = masked != original
                    called = BASE_CODE[original] < 4
                    if not np.array_equal(changed, expected) or np.any(masked[changed] != mask_byte):
                        raise RuntimeError(
                            f"The masked row of {name} does not match {_group_label(group_id)}'s recombination "
                            f"blocks ({int(np.count_nonzero(changed))} columns changed, "
                            f"{int(np.count_nonzero(expected))} expected)."
                        )
                    if name in source.kept:
                        row = masked
                        masked_positions = int(np.count_nonzero(changed))
                        # Gubbins' masking also turns gaps (absent gene, indel) into N.
                        masked_count = int(np.count_nonzero(called[changed]))
                    else:
                        reason = UNMASKED_DROPPED
                else:
                    if name not in unmasked_reasons:
                        raise RuntimeError(f"Genome {name} is neither in a Gubbins group nor recorded as unmasked.")
                    group_id, reason = unmasked_reasons[name]
                called_bases = int(np.count_nonzero(BASE_CODE[original] < 4))
                chunk = fasta_record(name, row)
                digest.update(chunk)
                out.write(chunk)
                rows.append(
                    {
                        "sample": name,
                        "group": group_id,
                        "masked_bases": masked_count,
                        "masked_positions": masked_positions,
                        "masked_fraction": round(masked_count / float(called_bases), 6) if called_bases else 0.0,
                        "unmasked_reason": reason,
                    }
                )
        for group_id, iterator in iterators.items():
            extra = next(iterator, None)
            if extra is not None:
                raise RuntimeError(f"The masked alignment of {_group_label(group_id)} has an extra row {extra[0]!r}.")
        os.replace(tmp, dest)
    finally:
        tmp.unlink(missing_ok=True)
    return {"sha256": digest.hexdigest(), "length": int(length or 0), "rows": rows}


def _write_masked_bases(rows: Sequence[Mapping[str, Any]], dest: Path) -> None:
    tmp = dest.with_name(dest.name + ".tmp")
    with tmp.open("w", encoding="utf-8") as out:
        out.write("sample\tgroup\tmasked_bases\tmasked_positions\tmasked_fraction\tunmasked_reason\n")
        for row in rows:
            out.write(
                f"{row['sample']}\t{row['group'] or ''}\t{row['masked_bases']}\t{row['masked_positions']}\t"
                f"{row['masked_fraction']}\t"
                f"{row['unmasked_reason'] or ''}\n"
            )
    os.replace(tmp, dest)


# ------------------------------------------------------------ running
# mask_gubbins_aln.py loads its whole group alignment in memory: at most this many at once.
MAX_PARALLEL_MASK_JOBS = 4
_CANCELLED = "cancelled after another Gubbins run failed"


class _Processes:
    """Tool processes of the parallel group runs; a failure stops the others.

    Terminating the ``docker run`` client does not stop its container (the signal
    reaches ``conda run`` as PID 1, which ignores it), so a docker command gets a
    ``--cidfile`` and cancel kills the container itself.
    """

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._procs: list[tuple[subprocess.Popen, Path | None]] = []
        self._cid_dir = Path(tempfile.mkdtemp(prefix="bactolomew_gubbins_cid."))
        self.cancelled = False

    def run(self, cmd: Sequence[str], cwd: Path) -> tuple[int, str, str]:
        cmd = list(cmd)
        with self._lock:
            if self.cancelled:
                raise RuntimeError(_CANCELLED)
            cidfile = None
            if cmd[:2] == ["docker", "run"]:
                cidfile = self._cid_dir / f"{len(self._procs)}.cid"
                cmd[2:2] = ["--cidfile", str(cidfile)]
            proc = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            self._procs.append((proc, cidfile))
        stdout, stderr = proc.communicate()
        return proc.returncode, stdout or "", stderr or ""

    def cancel(self) -> None:
        with self._lock:
            self.cancelled = True
            for proc, cidfile in self._procs:
                if proc.poll() is not None:
                    continue
                cid = cidfile.read_text().strip() if cidfile is not None and cidfile.is_file() else ""
                if cid:
                    subprocess.run(["docker", "kill", cid], capture_output=True, check=False)
                proc.terminate()

    def close(self) -> None:
        shutil.rmtree(self._cid_dir, ignore_errors=True)


class _Budget:
    """Threads in use by the running jobs: a job starts, in submission order, once its
    threads fit in the budget (a job larger than the budget runs alone)."""

    def __init__(self, total: int) -> None:
        self.total = max(1, int(total))
        self.used = 0
        self.next_ticket = 0
        self.cond = threading.Condition()

    def run(self, ticket: int, cost: int, fn: Any) -> Any:
        cost = min(max(1, int(cost)), self.total)
        with self.cond:
            self.cond.wait_for(lambda: self.next_ticket == ticket and self.used + cost <= self.total)
            self.used += cost
            self.next_ticket += 1
            self.cond.notify_all()
        try:
            return fn()
        finally:
            with self.cond:
                self.used -= cost
                self.cond.notify_all()


def _run_parallel(
    jobs: Mapping[str, Any], processes: _Processes, budget: int, costs: Mapping[str, int] | None = None
) -> dict[str, Any]:
    """Run callables (id -> fn) in parallel, in the given order, so the summed cost
    (threads; 1 each by default) of the running ones stays within budget. On the first
    failure the others are stopped and the failure is raised."""
    if not jobs:
        return {}
    gate = _Budget(budget)
    results: dict[str, Any] = {}
    with ThreadPoolExecutor(max_workers=len(jobs)) as pool:
        futures = {
            pool.submit(gate.run, ticket, (costs or {}).get(key, 1), fn): key
            for ticket, (key, fn) in enumerate(jobs.items())
        }
        done, _pending = wait(futures, return_when=FIRST_EXCEPTION)
        if any(f.exception() is not None for f in done):
            processes.cancel()
        wait(futures)
        errors = [f.exception() for f in futures if f.exception() is not None]
        real = [e for e in errors if _CANCELLED not in str(e)]
        if errors:
            raise (real or errors)[0]
        for future, key in futures.items():
            results[key] = future.result()
    return results


def _gubbins_command(run_dir: Path, input_path: Path, params: Mapping[str, Any], threads: int) -> list[str]:
    cmd = [
        *tool_cmd("run_gubbins.py", workdir=str(run_dir)),
        "--prefix", GUBBINS_PREFIX,
        "--threads", str(threads),
        "--seed", str(params["seed"]),
        "--iterations", str(params["iterations"]),
        "--converge-method", str(params["converge_method"]),
        "--filter-percentage", str(params["filter_percentage"]),
        "--tree-builder", str(params["tree_builder"]),
    ]
    if params.get("outgroup"):
        cmd += ["--outgroup", str(params["outgroup"])]
    cmd.append(str(input_path))
    return cmd


def _run_one_group(
    group: GubbinsGroup,
    input_path: Path,
    run_dir: Path,
    params: Mapping[str, Any],
    processes: _Processes,
    logger: logging.Logger,
) -> dict[str, Any]:
    cmd = _gubbins_command(run_dir, input_path, params, int(params["threads"]))
    logger.info(
        "Running Gubbins on %s (%d genomes, %d threads).", _group_label(group.id), len(group.members), params["threads"]
    )
    logger.info("CMD: " + " ".join(cmd))
    started = time.monotonic()
    returncode, stdout, stderr = processes.run(cmd, run_dir)
    wall = round(time.monotonic() - started, 1)
    (run_dir / "command.stdout.log").write_text(stdout, encoding="utf-8")
    (run_dir / "command.stderr.log").write_text(stderr, encoding="utf-8")
    if returncode != 0:
        raise RuntimeError(f"Gubbins failed on {_group_label(group.id)}: {_format_subprocess_error(stdout, stderr, returncode)}")
    missing = [
        name for name in GROUP_OUT_FILES
        if not (run_dir / name).is_file() or (run_dir / name).stat().st_size <= 0
    ]
    if missing:
        raise RuntimeError(f"Gubbins on {_group_label(group.id)} did not produce: " + ", ".join(missing))
    for leftover in run_dir.glob("tmp*"):
        if leftover.is_dir():
            shutil.rmtree(leftover, ignore_errors=True)
    logger.info("Gubbins on %s finished in %.0f s.", _group_label(group.id), wall)
    return {"command": [str(x) for x in cmd], "stdout": stdout, "wall_time_s": wall}


def group_warnings(entry: Mapping[str, Any], settings: Mapping[str, Any]) -> list[dict[str, str]]:
    """Warnings of a group that ran, from its summary entry and the current settings
    (rebuilt when a run is reused, so they name the group's current lineage id)."""
    group_id = str(entry.get("id"))
    warnings: list[dict[str, str]] = diversity_warnings(entry.get("diversity") or {}, settings, group_id)
    item = not_converged_warning(entry, int(entry.get("max_iterations") or settings["gubbins_iterations"]), group_id)
    if item:
        warnings.append(item)
    left_out = list(entry.get("samples_left_out") or [])
    if left_out:
        filter_percentage = float((entry.get("params") or {}).get("filter_percentage") or settings["gubbins_filter_percentage"])
        warnings.append(
            _warning(
                "gubbins_samples_unmasked",
                f"Gubbins left {len(left_out)} genome(s) of {_group_label(group_id)} out (more than "
                f"{filter_percentage:g}% N or gap in the core alignment): {', '.join(left_out)}. "
                "They stay in the core tree without recombination masking.",
            )
        )
    longest = entry.get("longest_branch_share") or {}
    share = longest.get("share")
    if share is not None and share > float(settings["gubbins_longest_branch_warn"]):
        if group_id == SINGLE_GROUP:
            text = (
                "the cohort holds distinct lineages. Gubbins is designed for closely related isolates; its "
                "recombination calls are most reliable within one lineage (phylo.gubbins_mode per_lineage, "
                "which needs fastGEAR or FastBAPS lineages)."
            )
        else:
            text = (
                f"lineage {group_id} itself holds distinct sublineages, so its recombination calls are less "
                "reliable."
            )
        warnings.append(
            _warning(
                "gubbins_lineage_structure",
                f"{share:.0%} of the clonal-frame SNPs of {_group_label(group_id)} sit on one branch "
                f"({longest.get('node')}): {text}",
            )
        )
    return warnings


def _group_entry(
    group: GubbinsGroup,
    folder: Path,
    staged: StagedGroup,
    params: Mapping[str, Any],
    settings: Mapping[str, Any],
    *,
    ordered_sha: str,
    version: str | None,
    run: Mapping[str, Any],
) -> dict[str, Any]:
    """Summary entry (and reuse stamp) of a group that ran, parsed from its folder."""
    kept = set(leaf_names(read_newick(folder / GUBBINS_FINAL_TREE)))
    entry = {
        "id": group.id,
        "members": list(group.members),
        "genomes": len(group.members),
        "ran": True,
        "reused": False,
        "skip_reason": None,
        "folder": None,
        "variable_sites": staged.variable_sites,
        "diversity": staged.diversity,
        "threads": int(params["threads"]),
        "command": list(run.get("command") or []),
        "input_sha256": staged.sha256,
        "ordered_alignment_sha256": ordered_sha,
        "params": dict(params),
        "tool_version": version,
        "algorithm_version": ALGORITHM_VERSION,
        **parse_convergence(str(run.get("stdout") or "")),
        "max_iterations": int(params["iterations"]),
        "recombination_blocks": len(gff_blocks(folder / GUBBINS_RECOMBINATION_GFF)),
        **recombination_statistics(folder / GUBBINS_PER_BRANCH, kept, root_edge_nodes(folder / GUBBINS_NODE_TREE)),
        "longest_branch_warn_above": float(settings["gubbins_longest_branch_warn"]),
        "samples_left_out": [name for name in group.members if name not in kept],
        "wall_time_s": run.get("wall_time_s"),
    }
    entry["warnings"] = group_warnings(entry, settings)
    return entry


def _skipped_entry(group: GubbinsGroup, staged: StagedGroup | None, settings: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "id": group.id,
        "members": list(group.members),
        "genomes": len(group.members),
        "ran": False,
        "reused": False,
        "skip_reason": group.skip_reason,
        "folder": None,
        "variable_sites": staged.variable_sites if staged else None,
        "diversity": staged.diversity if staged else None,
        "threads": None,
        "command": None,
        "input_sha256": staged.sha256 if staged else None,
        "converged": None,
        "iterations_run": None,
        "max_iterations": int(settings["gubbins_iterations"]),
        "convergence_reason": None,
        "recombination_blocks": None,
        "snps_inside_recombination": None,
        "snps_outside_recombination": None,
        "r_m_snp_ratio": None,
        "longest_branch_share": None,
        "tip_recombination_fraction_top5": [],
        "samples_left_out": [],
        "wall_time_s": None,
        "warnings": [],
    }


def _skip_warning(group: GubbinsGroup, settings: Mapping[str, Any], staged: StagedGroup | None) -> dict[str, str]:
    names = ", ".join(group.members)
    if group.skip_reason == SKIP_FEW_GENOMES:
        return _warning(
            "gubbins_skipped_few_genomes",
            f"Gubbins was skipped: {len(group.members)} genome(s), fewer than phylo.gubbins_min_genomes "
            f"({int(settings['gubbins_min_genomes'])}). The core tree is built on the ordered core "
            "alignment without recombination filtering.",
        )
    if group.skip_reason == SKIP_LINEAGE_SIZE:
        why = (
            f"it has {len(group.members)} genome(s), fewer than phylo.gubbins_min_lineage_size "
            f"({int(settings['gubbins_min_lineage_size'])})"
        )
    else:
        why = (
            f"it has {staged.variable_sites if staged else 0} variable column(s), fewer than the "
            f"{MIN_VARIABLE_SITES} Gubbins needs to call a block"
        )
    who = f"lineage {group.id}" if group.id != SINGLE_GROUP else "the cohort"
    return _warning(
        "gubbins_lineage_skipped",
        f"Gubbins was not run on {who}: {why}. Its genomes stay in the core tree without "
        f"recombination masking: {names}.",
    )


def _previous_group_stamps(gubbins_dir: Path) -> dict[tuple[str, ...], tuple[dict[str, Any], Path]]:
    """Stamps of the groups of the previous Gubbins step (schema 3): member key ->
    (entry, folder holding its files). Each run writes its stamp (run.json, in
    lineage_<id>/ or at the top of gubbins/ in single mode) as soon as it finishes, so a
    step that fails later (masking, merge) still reuses it on resume."""
    stamps: dict[tuple[str, ...], tuple[dict[str, Any], Path]] = {}
    for folder in [gubbins_dir, *sorted(gubbins_dir.glob("lineage_*"))]:
        try:
            entry = json.loads((folder / GROUP_STAMP).read_text(encoding="utf-8"))
        except Exception:
            continue
        if isinstance(entry, dict) and entry.get("ran"):
            stamps[tuple(sorted(entry.get("members") or []))] = (entry, folder)
    return stamps


def _stamp_reusable(
    stamp: Mapping[str, Any],
    folder: Path,
    group: GubbinsGroup,
    params: Mapping[str, Any],
    *,
    ordered_sha: str,
    version: str | None,
    input_sha: str | None,
) -> bool:
    return (
        stamp.get("algorithm_version") == ALGORITHM_VERSION
        and stamp.get("ordered_alignment_sha256") == ordered_sha
        and tuple(sorted(stamp.get("members") or [])) == group.key
        and _without_threads(stamp.get("params")) == _without_threads(params)
        and stamp.get("tool_version") == version
        and (input_sha is None or stamp.get("input_sha256") == input_sha)
        and all((folder / name).is_file() and (folder / name).stat().st_size > 0 for name in GROUP_OUT_FILES)
    )


def _single_layout_files(gubbins_dir: Path) -> list[Path]:
    return [
        path for path in gubbins_dir.iterdir()
        if path.is_file()
        and (path.name.startswith(f"{GUBBINS_PREFIX}.") or path.name.startswith("command.") or path.name == GROUP_STAMP)
    ]


def _move_into(run_dir: Path, gubbins_dir: Path) -> None:
    """Single mode: Gubbins' files go to the top of gubbins/."""
    for path in run_dir.iterdir():
        target = gubbins_dir / path.name
        if target.is_dir():
            shutil.rmtree(target)
        os.replace(path, target)
    run_dir.rmdir()


def _mask_group(input_path: Path, gff: Path, out: Path, workdir: Path, processes: _Processes) -> Path:
    cmd = [
        *tool_cmd("mask_gubbins_aln.py", workdir=str(workdir)),
        "--aln", str(input_path),
        "--gff", str(gff),
        "--out", str(out),
        "--missing-char", MASK_CHAR,
    ]
    returncode, stdout, stderr = processes.run(cmd, workdir)
    if returncode != 0 or not out.is_file():
        raise RuntimeError(f"mask_gubbins_aln.py failed: {_format_subprocess_error(stdout, stderr, returncode)}")
    return out


def run_gubbins(
    output_path: Path,
    settings: Mapping[str, Any],
    logger: logging.Logger,
    *,
    threads: int,
    lineages: Mapping[str, str] | None,
    lineage_error: str | None = None,
    lineage_source: str = "fastGEAR",
) -> GubbinsResult:
    """The gubbins step: groups by phylo.gubbins_mode, Gubbins per group (reusing the
    unchanged ones), Gubbins' masking per group merged into the masked core alignment,
    the combined block tables and gubbins/summary.json (written last)."""
    output_path = Path(output_path).expanduser().resolve()
    gubbins_dir = output_path / "gubbins"
    alignment = output_path / ORDERED_ALIGNMENT_REL
    ordering = read_ordering_summary(output_path)
    if ordering is None or not alignment.is_file():
        raise RuntimeError("The ordered core alignment (step core_alignment) is required before Gubbins.")
    ordered_sha = str((ordering.get("output") or {}).get("sha256") or "") or file_sha256(alignment)
    taxa = [name for name, _seq in iter_fasta_records(alignment)]
    plan = resolve_gubbins_groups(taxa, settings, lineages, lineage_error, lineage_source=lineage_source)
    single = plan.single
    gubbins_dir.mkdir(parents=True, exist_ok=True)
    previous = _previous_group_stamps(gubbins_dir)
    # Written last: a failure below leaves no summary, so the next run starts over.
    (gubbins_dir / GUBBINS_SUMMARY).unlink(missing_ok=True)
    warnings: list[dict[str, str]] = list(plan.warnings)
    notes: list[str] = list(plan.notes)
    dropouts = _predicted_dropouts(ordering, settings)
    outgroup = None
    if single:
        outgroup, outgroup_warnings, outgroup_notes = _resolve_outgroup(settings, dropouts)
        warnings += outgroup_warnings
        notes += outgroup_notes
    elif phylo_outgroup_names(settings):
        notes.append(
            "Gubbins runs per lineage without an outgroup (a cohort outgroup sits in another lineage); "
            "the cohort tree is rooted on phylo.outgroup afterwards."
        )
    base_params = gubbins_params(settings, threads=threads, outgroup=outgroup)
    candidates = [group for group in plan.groups if group.skip_reason is None]
    version = gubbins_version() if candidates else None
    logger.info(
        "Gubbins mode %s (requested %s, %s): %s.",
        plan.resolved_mode,
        plan.requested_mode,
        plan.mode_reason,
        ", ".join(f"{_group_label(g.id)} {len(g.members)} genomes" for g in plan.groups),
    )

    staging = Path(tempfile.mkdtemp(prefix=".stage.", dir=gubbins_dir))
    processes = _Processes()
    try:
        if single:
            diversity = dict(ordering.get("diversity") or {})
            staged = {
                g.id: StagedGroup(alignment, ordered_sha, int(diversity.get("variable_sites") or 0), diversity)
                for g in candidates
            }
        else:
            staged = stage_group_alignments(alignment, candidates, staging)
        for group in candidates:
            if staged[group.id].variable_sites < MIN_VARIABLE_SITES:
                group.skip_reason = SKIP_VARIABLE_SITES
        to_run = [group for group in plan.groups if group.skip_reason is None]
        for group in plan.groups:
            if group.skip_reason:
                item = _skip_warning(group, settings, staged.get(group.id))
                warnings.append(item)
                logger.warning(item["message"])
        if plan.unassigned:
            item = _warning(
                "gubbins_lineage_skipped",
                f"{len(plan.unassigned)} genome(s) have no {plan.source} lineage and stay in the core tree without "
                f"recombination masking: {', '.join(plan.unassigned)}.",
            )
            warnings.append(item)
            logger.warning(item["message"])
        predicted = [name for name in dropouts if any(name in g.members for g in to_run)]
        if predicted:
            logger.warning(
                "Gubbins' missing-data filter flags %d genome(s) with more than %g%% N or gap: %s. A genome "
                "missing from its group's final tree stays in the core tree unmasked.",
                len(predicted),
                float(settings["gubbins_filter_percentage"]),
                ", ".join(predicted),
            )

        # Reuse: a group with the same genomes, input and parameters keeps its run.
        group_params = {g.id: {**base_params} for g in to_run}
        reused: dict[str, tuple[dict[str, Any], Path]] = {}
        for group in to_run:
            stamp, folder = previous.get(group.key, (None, None))
            if stamp is not None and _stamp_reusable(
                stamp,
                folder,
                group,
                group_params[group.id],
                ordered_sha=ordered_sha,
                version=version,
                input_sha=None if single else staged[group.id].sha256,
            ):
                reused[group.id] = (stamp, folder)
        fresh = [group for group in to_run if group.id not in reused]
        _prepare_layout(gubbins_dir, single, reused, keep_single_files=bool(single and reused))
        split = split_threads(threads, {g.id: staged[g.id].variable_sites for g in fresh})
        for group in fresh:
            group_params[group.id]["threads"] = split[group.id]
        if reused:
            logger.info("Reusing Gubbins on %s (genomes, input and parameters unchanged).",
                        ", ".join(_group_label(g) for g in reused))

        # Run the changed groups in parallel (largest first, within the thread budget),
        # each into its own temp folder that becomes its final folder, with its stamp, as
        # soon as it finishes (a later failure keeps it for the resume).

        def run_group(group: GubbinsGroup, run_dir: Path) -> dict[str, Any]:
            run = _run_one_group(group, staged[group.id].path, run_dir, group_params[group.id], processes, logger)
            entry = _group_entry(
                group, run_dir, staged[group.id], group_params[group.id], settings,
                ordered_sha=ordered_sha, version=version, run=run,
            )
            if single:
                entry["folder"] = "."
                _move_into(run_dir, gubbins_dir)
                _write_json(gubbins_dir / GROUP_STAMP, entry)
            else:
                final = gubbins_dir / f"lineage_{group.id}"
                entry["folder"] = final.name
                _write_json(run_dir / GROUP_STAMP, entry)
                if final.exists():
                    shutil.rmtree(final)
                os.replace(run_dir, final)
            return entry

        jobs = {}
        for group in sorted(fresh, key=lambda g: -group_params[g.id]["threads"]):
            run_dir = Path(tempfile.mkdtemp(prefix=f".run.{group.id}.", dir=gubbins_dir))
            os.chmod(run_dir, 0o755)
            jobs[group.id] = lambda g=group, d=run_dir: run_group(g, d)
        entries: dict[str, dict[str, Any]] = _run_parallel(
            jobs, processes, threads, {g.id: group_params[g.id]["threads"] for g in fresh}
        )
        for group in to_run:
            if group.id in reused:
                entry = dict(reused[group.id][0])
                entry.update({"id": group.id, "reused": True, "members": list(group.members)})
                entry["warnings"] = group_warnings(entry, settings)
                entry["folder"] = "." if single else f"lineage_{group.id}"
                _write_json(gubbins_dir / ("" if single else entry["folder"]) / GROUP_STAMP, entry)
                entries[group.id] = entry
        for group in plan.groups:
            if group.skip_reason:
                entries[group.id] = _skipped_entry(group, staged.get(group.id), settings)
                if not single:
                    folder = gubbins_dir / f"lineage_{group.id}"
                    folder.mkdir(exist_ok=True)
                    _write_json(folder / GROUP_STAMP, entries[group.id])

        for group in to_run:
            warnings.extend(entries[group.id].get("warnings") or [])
            for item in entries[group.id].get("warnings") or []:
                if group.id not in reused:
                    logger.warning(item["message"])

        def folder_of(group_id: str) -> Path:
            return gubbins_dir if single else gubbins_dir / f"lineage_{group_id}"

        gffs = [(g.id, folder_of(g.id) / GUBBINS_RECOMBINATION_GFF) for g in to_run]
        masked_record = None
        block_stats: dict[str, dict[str, int]] = {}
        if to_run:
            # Gubbins' own masking per group, then one merged, checked cohort alignment.
            mask_jobs = {
                g.id: (
                    lambda g=g: _mask_group(
                        staged[g.id].path,
                        folder_of(g.id) / GUBBINS_RECOMBINATION_GFF,
                        staging / f"{g.id}.masked.aln",
                        staging,
                        processes,
                    )
                )
                for g in to_run
            }
            logger.info("Masking recombination with mask_gubbins_aln.py (%d group(s)).", len(mask_jobs))
            masked_files = _run_parallel(mask_jobs, processes, min(threads, MAX_PARALLEL_MASK_JOBS))
            sources = [
                MaskSource(
                    group_id=g.id,
                    members=frozenset(g.members),
                    kept=frozenset(leaf_names(read_newick(folder_of(g.id) / GUBBINS_FINAL_TREE))),
                    masked_alignment=masked_files[g.id],
                    gff=folder_of(g.id) / GUBBINS_RECOMBINATION_GFF,
                )
                for g in to_run
            ]
            reasons: dict[str, tuple[str | None, str]] = {
                name: (g.id, UNMASKED_LINEAGE_SKIPPED) for g in plan.groups if g.skip_reason for name in g.members
            }
            reasons.update({name: (None, UNMASKED_NO_LINEAGE) for name in plan.unassigned})
            merged = merge_masked_alignment(alignment, sources, reasons, gubbins_dir / MASKED_ALIGNMENT)
            _write_masked_bases(merged["rows"], gubbins_dir / MASKED_BASES_TSV)
            unmasked = [
                {"sample": r["sample"], "group": r["group"], "reason": r["unmasked_reason"]}
                for r in merged["rows"] if r["unmasked_reason"]
            ]
            masked_record = {
                "path": MASKED_ALIGNMENT_REL.as_posix(),
                "sha256": merged["sha256"],
                "length": merged["length"],
                "missing_char": MASK_CHAR,
                "inputs": {
                    "ordered_alignment_sha256": ordered_sha,
                    "gff_sha256": {group_id: file_sha256(gff) for group_id, gff in gffs},
                },
                "samples_unmasked": unmasked,
                "masked_bases": (Path("gubbins") / MASKED_BASES_TSV).as_posix(),
            }
            logger.info(
                "Masked core alignment: %d genomes, %d unmasked; %s.",
                len(merged["rows"]),
                len(unmasked),
                MASKED_ALIGNMENT_REL,
            )
            if not single:
                write_combined_gff(gffs, gubbins_dir / COMBINED_GFF)
            block_stats = write_block_gene_table(
                gffs, output_path / "core_alignment" / "gene_map.tsv", gubbins_dir / BLOCK_GENES_TSV
            )
            for group_id, stats in block_stats.items():
                entries[group_id].update(
                    {
                        "blocks_spanning_several_genes": stats["blocks_spanning_several_genes"],
                        "blocks_touching_appended_genes": stats["blocks_touching_appended_genes"],
                    }
                )
        else:
            for stale in (MASKED_ALIGNMENT, MASKED_BASES_TSV, BLOCK_GENES_TSV, COMBINED_GFF):
                (gubbins_dir / stale).unlink(missing_ok=True)
    finally:
        processes.close()
        shutil.rmtree(staging, ignore_errors=True)
        for leftover in gubbins_dir.glob(".run.*"):
            shutil.rmtree(leftover, ignore_errors=True)

    if not single and to_run:
        notes.append(INTER_LINEAGE_NOTE)
    ordered_entries = [entries[g.id] for g in plan.groups]
    ran = [e for e in ordered_entries if e.get("ran")]
    not_converged = [e["id"] for e in ran if e.get("converged") is not True]
    inside = sum(int(e.get("snps_inside_recombination") or 0) for e in ran)
    outside = sum(int(e.get("snps_outside_recombination") or 0) for e in ran)
    single_entry = ordered_entries[0] if single else {}
    summary = {
        "schema_version": SUMMARY_SCHEMA_VERSION,
        "algorithm_version": ALGORITHM_VERSION,
        "tool": "gubbins",
        "tool_version": version,
        "mode": plan.requested_mode,
        "resolved_mode": plan.resolved_mode,
        "mode_reason": plan.mode_reason,
        "lineage_source": plan.source,
        "ran": bool(ran),
        "skip_reason": None if ran else "no group ran (" + ", ".join(
            f"{_group_label(e['id'])}: {e['skip_reason']}" for e in ordered_entries
        ) + ")",
        "input_alignment": ORDERED_ALIGNMENT_REL.as_posix(),
        "input_sha256": ordered_sha,
        "ordering_genome": ordering.get("ordering_genome"),
        "params": base_params,
        "min_genomes": int(settings["gubbins_min_genomes"]),
        "min_lineage_size": int(settings["gubbins_min_lineage_size"]),
        "sequence_count": len(taxa),
        "converged": (not not_converged) if ran else None,
        "not_converged": not_converged,
        "iterations_run": max((int(e["iterations_run"]) for e in ran if e.get("iterations_run")), default=None),
        "max_iterations": int(settings["gubbins_iterations"]),
        "convergence_reason": single_entry.get("convergence_reason") if single else None,
        "diversity": dict(ordering.get("diversity") or {}),
        "recombination_blocks": sum(int(e.get("recombination_blocks") or 0) for e in ran),
        "snps_inside_recombination": inside,
        "snps_outside_recombination": outside,
        "r_m_snp_ratio": (inside / outside) if outside else None,
        # Genomes missing from the Gubbins tree; per lineage there is no cohort Gubbins tree,
        # and dropped genomes stay in the core tree unmasked (masked_alignment.samples_unmasked).
        "samples_left_out": list(single_entry.get("samples_left_out") or []) if single else [],
        "unassigned_samples": list(plan.unassigned),
        "groups": ordered_entries,
        "masked_alignment": masked_record,
        "warnings": warnings,
        "notes": notes,
        "outputs": {
            "masked_alignment": MASKED_ALIGNMENT_REL.as_posix() if masked_record else None,
            "masked_bases": (Path("gubbins") / MASKED_BASES_TSV).as_posix() if masked_record else None,
            "blocks_genes": (Path("gubbins") / BLOCK_GENES_TSV).as_posix() if masked_record else None,
            "recombination_gff": (
                (Path("gubbins") / (GUBBINS_RECOMBINATION_GFF if single else COMBINED_GFF)).as_posix()
                if masked_record else None
            ),
            "final_tree": (Path("gubbins") / GUBBINS_FINAL_TREE).as_posix() if single and ran else None,
        },
    }
    _write_json(gubbins_dir / GUBBINS_SUMMARY, summary)
    for entry in ran:
        logger.info(
            "Gubbins on %s %s after %s iteration(s): %s recombination blocks, r/m %.2f%s.",
            _group_label(entry["id"]),
            "converged" if entry.get("converged") else "did NOT converge",
            entry.get("iterations_run"),
            entry.get("recombination_blocks"),
            entry.get("r_m_snp_ratio") or 0.0,
            " (reused)" if entry.get("reused") else "",
        )
    return GubbinsResult(output_dir=gubbins_dir, ran=bool(ran), summary=summary, warnings=warnings)


def _prepare_layout(
    gubbins_dir: Path,
    single: bool,
    reused: Mapping[str, tuple[dict[str, Any], Path]],
    *,
    keep_single_files: bool,
) -> None:
    """Remove the previous layout's outputs that this run does not keep: the other
    mode's files, lineage folders of changed groups, legacy staging folders. Reused
    lineage folders are renamed to their current lineage id."""
    for name in ("input",):
        legacy = gubbins_dir / name
        if legacy.is_dir():
            shutil.rmtree(legacy)
    for leftover in gubbins_dir.glob("tmp*"):
        if leftover.is_dir():
            shutil.rmtree(leftover, ignore_errors=True)
    for leftover in gubbins_dir.glob(".run.*"):
        shutil.rmtree(leftover, ignore_errors=True)
    if single:
        for folder in gubbins_dir.glob("lineage_*"):
            shutil.rmtree(folder)
        (gubbins_dir / COMBINED_GFF).unlink(missing_ok=True)
        if not keep_single_files:
            for path in _single_layout_files(gubbins_dir):
                path.unlink()
        return
    for path in _single_layout_files(gubbins_dir):
        path.unlink()
    parked: dict[str, Path] = {}
    for group_id, (_stamp, folder) in reused.items():
        park = gubbins_dir / f".reuse.{group_id}"
        if park.exists():
            shutil.rmtree(park)
        os.replace(folder, park)
        parked[group_id] = park
    for folder in gubbins_dir.glob("lineage_*"):
        shutil.rmtree(folder)
    for group_id, park in parked.items():
        os.replace(park, gubbins_dir / f"lineage_{group_id}")


def read_gubbins_summary(output_path: Path) -> dict[str, Any] | None:
    try:
        data = json.loads((Path(output_path) / "gubbins" / GUBBINS_SUMMARY).read_text(encoding="utf-8"))
    except Exception:
        return None
    return data if isinstance(data, dict) else None


def masked_alignment_sha256(summary: Mapping[str, Any] | None) -> str | None:
    """The recorded sha256 of the masked core alignment (None when there is none)."""
    record = (summary or {}).get("masked_alignment")
    if isinstance(record, Mapping) and record.get("sha256"):
        return str(record["sha256"])
    return None


def gubbins_reuse_issue(
    output_path: Path,
    settings: Mapping[str, Any],
    *,
    threads: int,
    lineages: Mapping[str, str] | None,
    lineage_error: str | None = None,
    lineage_source: str = "fastGEAR",
) -> str | None:
    """Why the existing Gubbins output cannot be reused as a whole (None when it can).
    When it cannot, run_gubbins still reuses every group whose genomes, input and
    parameters did not change."""
    output_path = Path(output_path)
    summary = read_gubbins_summary(output_path)
    if summary is None:
        return "no Gubbins summary"
    if summary.get("schema_version") != SUMMARY_SCHEMA_VERSION:
        return "the Gubbins output predates per-lineage runs and the masked core alignment"
    if summary.get("algorithm_version") != ALGORITHM_VERSION:
        return "the Gubbins step changed"
    ordering = read_ordering_summary(output_path) or {}
    current_sha = (ordering.get("output") or {}).get("sha256")
    if not current_sha or summary.get("input_sha256") != current_sha:
        return "input alignment changed"
    alignment = output_path / ORDERED_ALIGNMENT_REL
    taxa = [name for name, _seq in iter_fasta_records(alignment)] if alignment.is_file() else []
    try:
        plan = resolve_gubbins_groups(taxa, settings, lineages, lineage_error, lineage_source=lineage_source)
    except RuntimeError as exc:
        return str(exc)
    if summary.get("resolved_mode") != plan.resolved_mode:
        return f"Gubbins mode changed ({summary.get('resolved_mode')} to {plan.resolved_mode})"
    recorded = {
        tuple(sorted(e.get("members") or [])): e for e in summary.get("groups") or [] if isinstance(e, Mapping)
    }
    if len(recorded) != len(plan.groups) or sorted(summary.get("unassigned_samples") or []) != sorted(plan.unassigned):
        return "the lineages changed"
    size_rules = (SKIP_FEW_GENOMES, SKIP_LINEAGE_SIZE)
    outgroup = None
    if plan.single:
        outgroup, _warnings, _notes = _resolve_outgroup(settings, _predicted_dropouts(ordering, settings))
    params = gubbins_params(settings, threads=threads, outgroup=outgroup)
    version_needed = False
    for group in plan.groups:
        entry = recorded.get(group.key)
        if entry is None:
            return "the lineages changed"
        if entry.get("id") != group.id:
            return "the lineage ids changed"
        recorded_size_skip = entry.get("skip_reason") if entry.get("skip_reason") in size_rules else None
        if recorded_size_skip != group.skip_reason:
            return "a lineage crossed phylo.gubbins_min_lineage_size or phylo.gubbins_min_genomes"
        if entry.get("ran"):
            if _without_threads(entry.get("params")) != _without_threads(params):
                return "parameters changed"
            version_needed = True
    if version_needed and summary.get("tool_version") != gubbins_version():
        return "Gubbins version changed"
    record = summary.get("masked_alignment")
    if summary.get("ran"):
        if not isinstance(record, Mapping):
            return "the masked core alignment record is missing"
        masked = output_path / str(record.get("path") or MASKED_ALIGNMENT_REL)
        if not masked.is_file() or file_sha256(masked) != record.get("sha256"):
            return "the masked core alignment is missing or was modified"
        single = summary.get("resolved_mode") == "single"
        for group_id, sha in ((record.get("inputs") or {}).get("gff_sha256") or {}).items():
            gff = output_path / group_folder_rel(str(group_id), single) / GUBBINS_RECOMBINATION_GFF
            if not gff.is_file() or file_sha256(gff) != sha:
                return f"the recombination calls of {_group_label(str(group_id))} changed"
    return None


__all__ = [
    "ALGORITHM_VERSION",
    "BLOCK_GENES_TSV",
    "COMBINED_GFF",
    "GROUP_OUT_FILES",
    "GUBBINS_FILTERED_ALIGNMENT",
    "GUBBINS_FINAL_TREE",
    "GUBBINS_OUT_FILES",
    "GUBBINS_PREFIX",
    "GUBBINS_RECOMBINATION_GFF",
    "GUBBINS_SUMMARY",
    "GroupPlan",
    "GubbinsGroup",
    "GubbinsResult",
    "INTER_LINEAGE_NOTE",
    "MASKED_ALIGNMENT",
    "MASKED_ALIGNMENT_REL",
    "MASKED_BASES_TSV",
    "MaskSource",
    "SINGLE_GROUP",
    "SUMMARY_SCHEMA_VERSION",
    "expected_gubbins_files",
    "group_folder_rel",
    "group_warnings",
    "gubbins_reuse_issue",
    "masked_alignment_sha256",
    "merge_masked_alignment",
    "not_converged_warning",
    "parse_convergence",
    "planned_gubbins_files",
    "read_fastgear_lineages",
    "read_gubbins_summary",
    "resolve_gubbins_groups",
    "run_gubbins",
    "split_threads",
    "stage_group_alignments",
]
