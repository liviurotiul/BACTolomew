"""Rooted tree files and ``trees/trees.json`` (the GUI's only source of tree choice).

Producers leave raw, unrooted trees and a stamp ``trees/<id>.json`` each:

- ``core_recomb`` / ``core``: ``pipeline.CoreTreeWorkflow`` (IQ-TREE on the
  recombination-masked core alignment, or the ordered one when Gubbins did not run);
- ``gubbins_final``: written here from ``gubbins/summary.json`` (the Gubbins final tree),
  only when Gubbins ran once on the whole cohort (per-lineage Gubbins trees hold one
  lineage each and are never listed as cohort trees);
- ``snippy_core``: ``pipeline.PhylogenyWorkflow`` (RAxML-NG on Snippy's core.full.aln),
  or ``adopt_legacy_snippy_tree`` for a run made before stamps.

``write_trees_index`` roots every raw tree by one rule (``phylo.outgroup`` when set and
usable, else midpoint; support moves with its edge), renames Snippy's ``Reference`` tip to
the reference sample, checks the tips against the run's genomes, writes ``trees/<id>.nwk``
and rebuilds ``trees.json`` from the stamps whose input still matches its upstream file.
It is cheap and runs after ``core_tree`` and ``phylogeny`` on every run. With one Gubbins
run on the whole cohort, changing ``phylo.outgroup`` is not only a re-root: Gubbins takes
the outgroup (``--outgroup``), so Gubbins and then the core tree rerun; per lineage it is
only a re-root. ``rooting.outgroup`` in ``trees.json`` is the list of outgroup names (or
null). The contract is section 4 of PHYLOGENY_REWORK_2026-10-10.md (schema 2: R2.7).

The core_recomb entry lists ``branch_length_sets`` when its refit on the unmasked
alignment (``CoreTreeWorkflow.run_recombination_refit``) matches the current raw tree and
ordered alignment: "Mutations only" (``core_recomb.nwk``) and "Including recombination"
(``core_recomb_with_recombination.nwk``, rooted on the same edge at the same relative
position as core_recomb, with its support values and children in the same order, so
both files list the leaves in the same order). Without a current refit the key is absent.
"""
from __future__ import annotations

import json
import logging
import os
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any, Mapping, Sequence

from pipeline.CoreAlignment import file_sha256, read_ordering_summary
from pipeline.CoreTreeWorkflow import (
    BRANCH_LENGTH_SET_LABELS,
    REFIT_DIR,
    REFIT_TREE_ID,
    TREE_LABELS,
    current_refit_stamp,
    tree_stamp_path,
)
from pipeline.GubbinsWorkflow import (
    GUBBINS_FINAL_TREE,
    SUMMARY_SCHEMA_VERSION as GUBBINS_SCHEMA_VERSION,
    masked_alignment_sha256,
    read_gubbins_summary,
)
from pipeline.Newick import (
    copy_split_supports,
    leaf_names,
    midpoint_root,
    order_like,
    read_newick,
    rename_leaves,
    root_on_outgroup,
    root_on_split,
    root_split,
    write_newick,
)
from pipeline.ToolEnv import phylo_outgroup_names

SCHEMA_VERSION = 2
TREES_JSON = "trees.json"
TREE_ORDER = ("core_recomb", "core", "gubbins_final", "snippy_core")
SNIPPY_RAW_REL = Path("trees") / "raw" / "snippy_core.raxml.nwk"
LEGACY_SNIPPY_TREE_REL = Path("trees") / "snippy_raxml.tree"
SNIPPY_ALIGNMENT_REL = Path("snippy") / "core" / "core.full.aln"
SNIPPY_ALGORITHM_VERSION = 1
GUBBINS_FINAL_ALGORITHM_VERSION = 1
RUN_WARNING_STEPS = ("core_alignment", "gubbins", "core_tree", "phylogeny")

_ENTRY_KEYS = (
    "id",
    "label",
    "file",
    "method",
    "tool_version",
    "alignment",
    "alignment_kind",
    "sites",
    "invariant_correction",
    "model",
    "log_likelihood",
    "support",
    "rooting",
    "recombination_filtered",
    "ordering_genome",
    "branch_length_unit",
    "samples",
    "samples_missing",
    "samples_unmasked",
    "reference_sample",
    "seed",
    "threads",
    "fast_mode",
    "notes",
)


def _warning(code: str, message: str) -> dict[str, str]:
    return {"step": "phylogeny", "code": code, "message": message}


def _read_json(path: Path) -> dict[str, Any] | None:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except Exception:
        return None
    return data if isinstance(data, dict) else None


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    os.replace(tmp, path)


def snippy_reference_stem(output_path: Path) -> str | None:
    try:
        text = (Path(output_path) / "snippy" / "reference.txt").read_text(encoding="utf-8").strip()
    except OSError:
        return None
    return Path(text).stem if text else None


def write_snippy_stamp(
    output_path: Path,
    *,
    tool_version: str | None,
    threads: int | None,
    small_tree_note: str | None = None,
) -> dict[str, Any]:
    """Stamp of the Snippy tree (trees/snippy_core.json) from the files on disk.
    small_tree_note: the tree was written without RAxML-NG (fewer than 4 sequences)."""
    output_path = Path(output_path)
    alignment = output_path / SNIPPY_ALIGNMENT_REL
    stamp = {
        "id": "snippy_core",
        "label": TREE_LABELS["snippy_core"],
        "raw_file": str(SNIPPY_RAW_REL),
        "input_path": str(SNIPPY_ALIGNMENT_REL),
        "input_sha256": file_sha256(alignment) if alignment.is_file() else None,
        "algorithm_version": SNIPPY_ALGORITHM_VERSION,
        "params": {"model": "GTR+G", "search": "raxml-ng --search", "seed": 1},
        "tool_version": tool_version,
        "threads": threads,
        "seed": 1,
        "method": f"RAxML-NG{(' ' + tool_version) if tool_version else ''}, GTR+G, best tree (no support)",
        "alignment": str(SNIPPY_ALIGNMENT_REL),
        "alignment_kind": "snippy_core_full",
        "sites": None,
        "invariant_correction": None,
        "model": "GTR+G",
        "log_likelihood": None,
        "support": {"type": "none", "replicates": 0, "threshold": None},
        "recombination_filtered": False,
        "ordering_genome": None,
        "branch_length_unit": "substitutions/site",
        "samples_missing": [],
        "reference_sample": snippy_reference_stem(output_path),
        "fast_mode": False,
        "notes": [],
    }
    if small_tree_note:
        stamp.update(
            {
                "method": "only possible tree (fewer than 4 sequences, no RAxML-NG)",
                "model": None,
                "params": {"search": "none"},
                "notes": [small_tree_note],
            }
        )
    _write_json(tree_stamp_path(output_path, "snippy_core"), stamp)
    return stamp


def adopt_legacy_snippy_tree(output_path: Path, logger: logging.Logger) -> bool:
    """Copy an old run's trees/snippy_raxml.tree to trees/raw/ and stamp it (no Snippy
    call). The legacy file stays where it was for the scripts that still read it."""
    output_path = Path(output_path)
    legacy = output_path / LEGACY_SNIPPY_TREE_REL
    raw = output_path / SNIPPY_RAW_REL
    if legacy.is_file() and not raw.is_file():
        raw.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(legacy, raw)
        logger.info("Copied the Snippy tree %s to %s.", legacy, raw)
    if raw.is_file() and not tree_stamp_path(output_path, "snippy_core").is_file():
        write_snippy_stamp(output_path, tool_version=None, threads=None)
        return True
    return False


def write_gubbins_final_stamp(output_path: Path) -> dict[str, Any] | None:
    output_path = Path(output_path)
    summary = read_gubbins_summary(output_path)
    stamp_path = tree_stamp_path(output_path, "gubbins_final")
    single = bool(summary) and (
        summary.get("schema_version") != GUBBINS_SCHEMA_VERSION or summary.get("resolved_mode") == "single"
    )
    if (
        not summary
        or not single
        or not summary.get("ran", True)
        or not (output_path / "gubbins" / GUBBINS_FINAL_TREE).is_file()
    ):
        stamp_path.unlink(missing_ok=True)
        (output_path / "trees" / "gubbins_final.nwk").unlink(missing_ok=True)
        return None
    params = summary.get("params") or {}
    version = summary.get("tool_version")
    stamp = {
        "id": "gubbins_final",
        "label": TREE_LABELS["gubbins_final"],
        "raw_file": str(Path("gubbins") / GUBBINS_FINAL_TREE),
        "input_path": str(Path("core_alignment") / "core_gene_alignment.ordered.aln"),
        "input_sha256": summary.get("input_sha256"),
        "algorithm_version": GUBBINS_FINAL_ALGORITHM_VERSION,
        "params": params,
        "tool_version": version,
        "threads": params.get("threads"),
        "seed": params.get("seed"),
        "method": (
            f"Gubbins{(' ' + version) if version else ''} final tree "
            f"({params.get('tree_builder') or 'raxml'} inside Gubbins), no support"
        ),
        "alignment": str(Path("gubbins") / "core_gene.filtered_polymorphic_sites.fasta"),
        "alignment_kind": "core_snps",
        "sites": None,
        "invariant_correction": None,
        "model": None,
        "log_likelihood": None,
        "support": {"type": "none", "replicates": 0, "threshold": None},
        "recombination_filtered": True,
        "ordering_genome": summary.get("ordering_genome"),
        "branch_length_unit": "SNPs (reconstructed point mutations)",
        "samples_missing": list(summary.get("samples_left_out") or []),
        "samples_unmasked": [],
        "reference_sample": None,
        "fast_mode": False,
        "notes": [],
    }
    _write_json(stamp_path, stamp)
    return stamp


def _current_input_sha(output_path: Path, tree_id: str, stamp: Mapping[str, Any]) -> str | None:
    if tree_id == "core_recomb":
        # The masked alignment's sha256 as recorded and checked by the Gubbins step.
        return masked_alignment_sha256(read_gubbins_summary(output_path))
    if tree_id in ("core", "gubbins_final"):
        return ((read_ordering_summary(output_path) or {}).get("output") or {}).get("sha256")
    if tree_id == "snippy_core":
        path = output_path / SNIPPY_ALIGNMENT_REL
        return file_sha256(path) if path.is_file() else None
    return None


_LINEAGE_KEYS = (
    "id",
    "genomes",
    "ran",
    "skip_reason",
    "converged",
    "iterations_run",
    "max_iterations",
    "convergence_reason",
    "recombination_blocks",
    "r_m_snp_ratio",
)


def gubbins_index_block(gubbins: Mapping[str, Any] | None) -> dict[str, Any] | None:
    """trees.json ``gubbins``: the schema 1 keys with values aggregated over the Gubbins
    groups (converged only when every group that ran converged), plus the mode and one
    entry per group (``lineages``; the one group ``all`` in single mode)."""
    if gubbins is None:
        return None
    params = gubbins.get("params") or {}
    groups = [g for g in gubbins.get("groups") or [] if isinstance(g, Mapping)]
    return {
        "ran": bool(gubbins.get("ran", True)),
        "converged": gubbins.get("converged"),
        "iterations_run": gubbins.get("iterations_run"),
        "max_iterations": gubbins.get("max_iterations") or params.get("iterations"),
        "converge_method": params.get("converge_method"),
        "samples_left_out": list(gubbins.get("samples_left_out") or []),
        "mode": gubbins.get("mode"),
        "resolved_mode": gubbins.get("resolved_mode") or ("single" if gubbins else None),
        "lineages": [{key: group.get(key) for key in _LINEAGE_KEYS} for group in groups],
    }


def root_tree(tree, outgroup: Sequence[str]) -> tuple[Any, dict[str, Any], dict[str, str] | None]:
    """(rooted tree, rooting record, warning or None) by the run's rule."""
    if outgroup:
        rooted, reason = root_on_outgroup(tree, outgroup)
        if rooted is not None:
            return rooted, {"method": "outgroup", "outgroup": list(outgroup), "note": None}, None
        if reason == "not_monophyletic":
            code = "outgroup_not_monophyletic"
            note = f"The outgroup {', '.join(outgroup)} is not monophyletic in this tree; rooted at midpoint."
        else:
            code = "outgroup_unusable"
            note = f"Outgroup taxa are not in this tree ({reason.split(': ', 1)[-1]}); rooted at midpoint."
        return midpoint_root(tree), {"method": "midpoint", "outgroup": list(outgroup), "note": note}, _warning(code, note)
    return midpoint_root(tree), {"method": "midpoint", "outgroup": None, "note": None}, None


def _write_tree(path: Path, tree) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(write_newick(tree), encoding="utf-8")
    os.replace(tmp, path)


def _branch_length_sets(
    output_path: Path, stamp: Mapping[str, Any], raw: Path, rooted, logger: logging.Logger
) -> list[dict[str, Any]] | None:
    """The core_recomb entry's branch_length_sets, writing the rooted refit tree; None
    (and no refit file) when there is no current refit."""
    out_file = output_path / "trees" / f"{REFIT_TREE_ID}.nwk"
    refit_raw = output_path / REFIT_DIR / f"{REFIT_TREE_ID}.treefile"
    split = root_split(rooted)
    refit = current_refit_stamp(output_path, raw)
    refit_rooted = root_on_split(read_newick(refit_raw), *split) if refit is not None and split else None
    if refit is None or refit_rooted is None:
        if refit is not None:
            logger.warning("The refit core tree has no edge for the core tree's root; it is not listed.")
        out_file.unlink(missing_ok=True)
        return None
    # Each split takes the core tree's support, whatever labels IQ-TREE's -te output kept.
    copy_split_supports(refit_rooted, rooted)
    _write_tree(out_file, order_like(refit_rooted, rooted))
    return [
        {
            "id": "mutations",
            "label": BRANCH_LENGTH_SET_LABELS["mutations"],
            "file": f"{stamp.get('id')}.nwk",
            "default": True,
            "alignment": stamp.get("alignment"),
            "model": stamp.get("model"),
            "branch_length_unit": stamp.get("branch_length_unit"),
            "tree_length": (refit.get("tree_length") or {}).get("mutations"),
        },
        {
            "id": "with_recombination",
            "label": BRANCH_LENGTH_SET_LABELS["with_recombination"],
            "file": out_file.name,
            "default": False,
            "alignment": refit.get("alignment"),
            "model": refit.get("model"),
            "branch_length_unit": refit.get("branch_length_unit") or "substitutions/site",
            "method": (
                f"IQ-TREE{(' ' + refit['tool_version']) if refit.get('tool_version') else ''}, topology of the "
                "core tree fixed, branch lengths re-estimated on the unmasked core alignment"
            ),
            "log_likelihood": refit.get("log_likelihood"),
            "tree_length": (refit.get("tree_length") or {}).get("with_recombination"),
        },
    ]


def write_trees_index(
    output_path: Path,
    genomes: Sequence[Any],
    settings: Mapping[str, Any],
    logger: logging.Logger,
    *,
    run_warnings: Sequence[Mapping[str, str]] = (),
) -> tuple[dict[str, Any], list[dict[str, str]]]:
    """Root the raw trees, write trees/<id>.nwk and trees/trees.json. Returns
    (trees.json payload, warnings raised here)."""
    output_path = Path(output_path)
    trees_dir = output_path / "trees"
    trees_dir.mkdir(parents=True, exist_ok=True)
    expected = {Path(str(g)).stem for g in genomes}
    outgroup = phylo_outgroup_names(settings)
    write_gubbins_final_stamp(output_path)
    warnings: list[dict[str, str]] = []
    entries: list[dict[str, Any]] = []
    for tree_id in TREE_ORDER:
        stamp = _read_json(tree_stamp_path(output_path, tree_id))
        out_file = trees_dir / f"{tree_id}.nwk"
        if stamp is None:
            out_file.unlink(missing_ok=True)
            continue
        raw = output_path / str(stamp.get("raw_file") or "")
        current = _current_input_sha(output_path, tree_id, stamp)
        if not raw.is_file() or not current or stamp.get("input_sha256") != current:
            message = f"Tree \"{TREE_LABELS[tree_id]}\" is stale (its input changed) and is not listed."
            warnings.append(_warning("tree_stale", message))
            logger.warning(message)
            out_file.unlink(missing_ok=True)
            continue
        tree = read_newick(raw)
        reference = stamp.get("reference_sample")
        if tree_id == "snippy_core" and reference:
            leaves = set(leaf_names(tree))
            if "Reference" in leaves and reference not in leaves:
                rename_leaves(tree, lambda name: reference if name == "Reference" else name)
        tips = set(leaf_names(tree))
        missing = sorted(expected - tips)
        allowed_missing = set(stamp.get("samples_missing") or [])
        unexpected = sorted(tips - expected)
        if unexpected or set(missing) - allowed_missing:
            message = (
                f"The tips of tree \"{TREE_LABELS[tree_id]}\" do not match the run's genomes "
                f"(missing {len(set(missing) - allowed_missing)}, unexpected {len(unexpected)}: "
                f"{', '.join((sorted(set(missing) - allowed_missing) + unexpected)[:5])}); it is not listed."
            )
            warnings.append(_warning("tree_tip_mismatch", message))
            logger.warning(message)
            out_file.unlink(missing_ok=True)
            continue
        rooted, rooting, root_warning = root_tree(tree, outgroup)
        if root_warning:
            root_warning = {**root_warning, "message": f"{TREE_LABELS[tree_id]}: {root_warning['message']}"}
            warnings.append(root_warning)
            logger.warning(root_warning["message"])
        _write_tree(out_file, rooted)
        entry = {key: stamp.get(key) for key in _ENTRY_KEYS}
        entry.update(
            {
                "id": tree_id,
                "label": TREE_LABELS[tree_id],
                "file": out_file.name,
                "rooting": rooting,
                "samples": len(tips),
                "samples_missing": missing,
                "notes": list(stamp.get("notes") or []),
            }
        )
        if tree_id == "core_recomb":
            sets = _branch_length_sets(output_path, stamp, raw, rooted, logger)
            if sets:
                entry["branch_length_sets"] = sets
        entries.append(entry)
    if not any("branch_length_sets" in entry for entry in entries):
        (trees_dir / f"{REFIT_TREE_ID}.nwk").unlink(missing_ok=True)

    gubbins_block = gubbins_index_block(read_gubbins_summary(output_path))
    default = next((e["id"] for e in entries if e["recombination_filtered"]), entries[0]["id"] if entries else None)
    default_entry = next((e for e in entries if e["id"] == default), None)
    ordering = read_ordering_summary(output_path) or {}
    combined: list[dict[str, str]] = []
    for item in [*run_warnings, *warnings]:
        if item.get("step") in RUN_WARNING_STEPS and dict(item) not in combined:
            combined.append(dict(item))
    payload = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "default_tree": default,
        "rooting": default_entry["rooting"] if default_entry else {"method": "midpoint", "outgroup": None, "note": None},
        "ordering_genome": ordering.get("ordering_genome"),
        "gubbins": gubbins_block,
        "warnings": combined,
        "trees": entries,
    }
    _write_json(trees_dir / TREES_JSON, payload)
    logger.info(
        "Wrote %s: %s (default %s, rooted %s).",
        trees_dir / TREES_JSON,
        ", ".join(e["id"] for e in entries) or "no trees",
        default,
        payload["rooting"]["method"],
    )
    return payload, warnings


def invalidate_trees(output_path: Path, tree_ids: Sequence[str]) -> None:
    """Delete the listed trees (rooted file and stamp) before an upstream rerun, so a
    downstream failure never leaves an old tree beside new upstream results."""
    trees_dir = Path(output_path) / "trees"
    for tree_id in tree_ids:
        for suffix in (".nwk", ".json"):
            (trees_dir / f"{tree_id}{suffix}").unlink(missing_ok=True)
        if tree_id == "core_recomb":
            (trees_dir / f"{REFIT_TREE_ID}.nwk").unlink(missing_ok=True)


__all__ = [
    "LEGACY_SNIPPY_TREE_REL",
    "SNIPPY_RAW_REL",
    "TREES_JSON",
    "adopt_legacy_snippy_tree",
    "gubbins_index_block",
    "invalidate_trees",
    "snippy_reference_stem",
    "write_snippy_stamp",
    "write_trees_index",
]
