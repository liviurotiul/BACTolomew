"""Core genome tree with IQ-TREE (step ``core_tree``).

The default tree of a run (``core_recomb``) is IQ-TREE on the recombination-masked core
alignment (``gubbins/core_gene_alignment.masked.aln``: the ordered core alignment with
each genome's recombinant columns set to N by its own Gubbins group, see
``pipeline.GubbinsWorkflow``). The invariant sites are in the alignment, so no constant
site correction is needed and branch lengths are substitutions per core site.
ModelFinder picks the model (``-m MFP``) and 1000 ultrafast bootstraps (``-B 1000
--bnni``) give branch support, with a fixed seed (multithreaded IQ-TREE is still not
bit-reproducible: on the same input the log-likelihood and support can differ in the
last digits and poorly supported splits (UFBoot well under 95) can come out
differently, which can move leaves in the drawn order; the thread count is recorded, not part of the
reuse check). A reused tree takes its Gubbins-derived notes and unmasked genomes fresh
from the Gubbins summary (``refresh_gubbins_notes``).

When no Gubbins group ran (too few genomes, or every lineage too small) the tree is
built on the ordered core alignment itself (``core``, not recombination-filtered).
Fewer than four genomes: no bootstrap; two or fewer: a trivial tree. Above
``phylo.iqtree_fast_above_genomes`` IQ-TREE runs ``--fast`` with SH-aLRT support. Every
choice goes into ``trees/core_recomb_iqtree/run.json`` (the reuse stamp) and the tree's
``trees/<id>.json`` stamp, from which ``pipeline.TreesIndex`` roots it and lists it in
``trees/trees.json``.

Branch lengths including recombination (``run_recombination_refit``, after the core
tree): the core_recomb topology refit with IQ-TREE ``-te`` on the unmasked ordered
alignment, same model, no tree search, no bootstrap. It gives the same tree a second set
of branch lengths (``trees/core_recomb_refit/``, listed by ``pipeline.TreesIndex`` as the
core_recomb entry's ``branch_length_sets``). It has its own reuse stamp and never makes
the core tree rerun.
"""
from __future__ import annotations

import json
import logging
import os
import re
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from pipeline.CoreAlignment import (
    BASE_CODE,
    ORDERED_ALIGNMENT_REL,
    file_sha256,
    iter_fasta_records,
    ordered_alignment_sha256,
    read_ordering_summary,
)
from pipeline.GubbinsWorkflow import (
    GUBBINS_FINAL_TREE,
    INTER_LINEAGE_NOTE,
    MASKED_ALIGNMENT_REL,
    MASKED_BASES_TSV,
    SUMMARY_SCHEMA_VERSION,
    group_folder_rel,
    masked_alignment_sha256,
    read_gubbins_summary,
)
from pipeline.Newick import (
    leaf_names,
    read_newick,
    robinson_foulds,
    split_lengths,
    split_supports,
    support_values,
)
from pipeline.ToolEnv import tool_cmd

ALGORITHM_VERSION = 2
IQTREE_DIR = Path("trees") / "core_recomb_iqtree"
RUN_STAMP = "run.json"
KEPT_SUFFIXES = (".iqtree", ".log", ".treefile", ".contree", ".model.gz")
MIN_GENOMES_FOR_SUPPORT = 4

TREE_LABELS = {
    "core_recomb": "Core genome, recombination-filtered",
    "core": "Core genome (Gubbins skipped)",
    "gubbins_final": "Gubbins final tree, recombination-filtered (no support)",
    "snippy_core": "Whole-genome SNPs (Snippy), not recombination-corrected",
}
SUPPORT_THRESHOLDS = {"ufboot": 95, "sh_alrt": 80}


@dataclass(frozen=True)
class CoreTreeResult:
    tree_id: str
    raw_tree: Path
    stamp: dict[str, Any]
    warnings: list[dict[str, str]] = field(default_factory=list)


def _warning(code: str, message: str) -> dict[str, str]:
    return {"step": "core_tree", "code": code, "message": message}


def iqtree_version() -> str | None:
    try:
        result = subprocess.run([*tool_cmd("iqtree2"), "--version"], capture_output=True, text=True, check=False)
    except FileNotFoundError:
        return None
    found = re.search(r"version\s+(\d+(?:\.\d+)+)", (result.stdout or "") + (result.stderr or ""))
    return found.group(1) if found else None


@dataclass(frozen=True)
class CoreTreePlan:
    tree_id: str
    input_path: Path
    input_sha256: str
    taxa: int
    fast_mode: bool
    bootstraps: int
    support_type: str
    model: str


def plan_core_tree(output_path: Path, settings: Mapping[str, Any]) -> CoreTreePlan:
    """What the core tree step will build, from the Gubbins summary and settings."""
    output_path = Path(output_path)
    gubbins = read_gubbins_summary(output_path)
    if gubbins is None or gubbins.get("schema_version") != SUMMARY_SCHEMA_VERSION:
        raise RuntimeError("The Gubbins step (schema 3 summary, or its skip record) is required before the core tree.")
    masked_sha = masked_alignment_sha256(gubbins)
    if gubbins.get("ran") and masked_sha:
        tree_id = "core_recomb"
        input_path = output_path / MASKED_ALIGNMENT_REL
        input_sha = masked_sha
    else:
        tree_id = "core"
        input_path = output_path / ORDERED_ALIGNMENT_REL
        input_sha = str(((read_ordering_summary(output_path) or {}).get("output") or {}).get("sha256") or "")
    if not input_path.is_file() or not input_sha:
        raise RuntimeError(f"Core tree input is missing: {input_path}")
    taxa = int(gubbins.get("sequence_count") or 0) or sum(1 for _ in iter_fasta_records(input_path))
    fast = taxa > int(settings["iqtree_fast_above_genomes"])
    bootstraps = int(settings["iqtree_bootstraps"]) if taxa >= MIN_GENOMES_FOR_SUPPORT else 0
    support_type = "none" if bootstraps == 0 else ("sh_alrt" if fast else "ufboot")
    return CoreTreePlan(
        tree_id=tree_id,
        input_path=input_path,
        input_sha256=input_sha,
        taxa=taxa,
        fast_mode=fast,
        bootstraps=bootstraps,
        support_type=support_type,
        model=str(settings["iqtree_model"]),
    )


def _stamp_params(plan: CoreTreePlan, settings: Mapping[str, Any], *, threads: int) -> dict[str, Any]:
    return {
        "tree_id": plan.tree_id,
        "model": plan.model,
        "bootstraps": plan.bootstraps,
        "support_type": plan.support_type,
        "bnni": bool(settings["iqtree_bnni"]) and plan.support_type == "ufboot",
        "fast_mode": plan.fast_mode,
        "seed": int(settings["seed"]),
        "threads": int(max(1, threads)),
    }


# ------------------------------------------------------------ IQ-TREE parsing
def parse_iqtree_outputs(prefix: Path) -> dict[str, Any]:
    report = Path(str(prefix) + ".iqtree")
    log = Path(str(prefix) + ".log")
    text = report.read_text(encoding="utf-8", errors="ignore") if report.is_file() else ""
    log_text = log.read_text(encoding="utf-8", errors="ignore") if log.is_file() else ""
    best = re.search(r"Best-fit model:\s*(\S+)\s+chosen according to (\w+)", text + "\n" + log_text)
    model_line = re.search(r"^Model of substitution:\s*(\S+)", text, flags=re.MULTILINE)
    loglik = re.search(r"Log-likelihood of the tree:\s*(-?[\d.]+)", text)
    composition_failed = re.findall(r"^\s*\d+\s+(\S+)\s+[\d.]+%\s+failed\s", log_text, flags=re.MULTILINE)
    identical = re.findall(r"NOTE:\s*(\S+) is identical to (\S+) but kept", log_text)
    ignored = re.search(r"NOTE:\s*(\d+) identical sequences .*? will be ignored", log_text)
    return {
        "model": best.group(1) if best else (model_line.group(1) if model_line else None),
        "model_criterion": best.group(2) if best else None,
        "log_likelihood": float(loglik.group(1)) if loglik else None,
        "ufboot_converged": not bool(re.search(r"bootstrap analysis did not converge", log_text, flags=re.IGNORECASE)),
        "composition_failed": composition_failed,
        "identical_kept": len(identical),
        "identical_ignored": int(ignored.group(1)) if ignored else 0,
    }


def _alignment_shape(path: Path) -> tuple[list[str], int, int]:
    """(taxa, columns, variable columns) of an alignment: a column is variable when it
    holds two different bases among A, C, G and T."""
    names: list[str] = []
    first: np.ndarray | None = None
    variable: np.ndarray | None = None
    for name, raw in iter_fasta_records(path):
        names.append(name)
        codes = BASE_CODE[np.frombuffer(raw.upper(), dtype=np.uint8)]
        called = codes < 4
        if first is None:
            first = np.where(called, codes, 4).astype(np.uint8)
            variable = np.zeros(codes.size, dtype=bool)
            continue
        np.copyto(first, codes, where=called & (first == 4))
        variable |= called & (codes != first)
    if first is None or variable is None:
        return names, 0, 0
    return names, int(first.size), int(np.count_nonzero(variable))


def _unmasked_names(gubbins: Mapping[str, Any]) -> list[str]:
    record = gubbins.get("masked_alignment") or {}
    return [str(item.get("sample")) for item in record.get("samples_unmasked") or [] if isinstance(item, Mapping)]


_UNMASKED_TEXT = {
    "lineage_skipped": "lineage not run by Gubbins",
    "dropped_by_gubbins": "left out by Gubbins (too many N or gaps)",
    "no_lineage": "no fastGEAR lineage",
}


def _masking_notes(gubbins: Mapping[str, Any]) -> list[str]:
    notes: list[str] = []
    if gubbins.get("resolved_mode") == "per_lineage":
        notes.append(INTER_LINEAGE_NOTE)
    record = gubbins.get("masked_alignment") or {}
    by_reason: dict[str, list[str]] = {}
    for item in record.get("samples_unmasked") or []:
        if isinstance(item, Mapping):
            by_reason.setdefault(str(item.get("reason")), []).append(str(item.get("sample")))
    for reason, names in by_reason.items():
        notes.append(
            f"{len(names)} genome(s) not recombination-masked ({_UNMASKED_TEXT.get(reason, reason)}): "
            + ", ".join(names)
        )
    return notes


def _rf_to_gubbins_trees(output_path: Path, gubbins: Mapping[str, Any], tree) -> list[dict[str, Any]]:
    """RF of the core tree to each Gubbins tree on their shared genomes (one cohort tree
    in single mode, one per lineage otherwise)."""
    single = gubbins.get("resolved_mode") == "single"
    out = []
    for group in gubbins.get("groups") or []:
        if not isinstance(group, Mapping) or not group.get("ran"):
            continue
        path = output_path / group_folder_rel(str(group.get("id")), single) / GUBBINS_FINAL_TREE
        if not path.is_file():
            continue
        distance, normalized, shared = robinson_foulds(tree, read_newick(path))
        if shared >= 4:
            out.append(
                {"group": str(group.get("id")), "rf": distance, "normalized_rf": round(normalized, 4), "shared_leaves": shared}
            )
    return out


def _gubbins_notes(output_path: Path, gubbins: Mapping[str, Any], tree) -> tuple[list[dict[str, Any]], list[str]]:
    """(RF to the Gubbins trees, notes) of a core_recomb tree: the parts of its stamp that
    come from the Gubbins step, so a reused tree can take them fresh (lineage ids,
    unmasked genomes) without rerunning IQ-TREE."""
    rf = _rf_to_gubbins_trees(output_path, gubbins, tree)
    notes = []
    for item in rf:
        where = "the Gubbins final tree" if item["group"] == "all" else f"the lineage {item['group']} Gubbins tree"
        notes.append(
            f"Robinson-Foulds distance to {where}: {item['rf']} (normalized {item['normalized_rf']:.3f}, "
            f"{item['shared_leaves']} shared genomes); recombination was reconstructed on the Gubbins "
            "tree, not this one."
        )
    return rf, notes + _masking_notes(gubbins)


def refresh_gubbins_notes(output_path: Path) -> bool:
    """A reused core_recomb tree: rewrite its Gubbins-derived notes, RF distances and
    unmasked genomes from the current gubbins/summary.json (the masked alignment can be
    byte-identical while the lineage ids or the unmasked set changed). True when the
    stamps changed."""
    output_path = Path(output_path)
    run_stamp = read_core_tree_stamp(output_path)
    if not run_stamp or (run_stamp.get("params") or {}).get("tree_id") != "core_recomb" or "iqtree_notes" not in run_stamp:
        return False
    tree_path = tree_stamp_path(output_path, "core_recomb")
    try:
        tree_stamp = json.loads(tree_path.read_text(encoding="utf-8"))
        tree = read_newick(output_path / IQTREE_DIR / "core_recomb.treefile")
    except Exception:
        return False
    gubbins = read_gubbins_summary(output_path) or {}
    rf, gubbins_notes = _gubbins_notes(output_path, gubbins, tree)
    notes = list(run_stamp["iqtree_notes"]) + gubbins_notes
    unmasked = _unmasked_names(gubbins)
    if (
        run_stamp.get("notes") == notes
        and run_stamp.get("rf_to_gubbins") == rf
        and tree_stamp.get("notes") == notes
        and tree_stamp.get("samples_unmasked") == unmasked
    ):
        return False
    _write_json(output_path / IQTREE_DIR / RUN_STAMP, {**run_stamp, "notes": notes, "rf_to_gubbins": rf})
    _write_json(tree_path, {**tree_stamp, "notes": notes, "samples_unmasked": unmasked})
    return True


def _write_trivial_tree(names: Sequence[str], dest: Path) -> None:
    if not names:
        raise RuntimeError("Cannot write a tree without taxa.")
    text = f"{names[0]};\n" if len(names) == 1 else "(" + ",".join(f"{n}:0" for n in names) + ");\n"
    dest.write_text(text, encoding="utf-8")


def _format_error(result: subprocess.CompletedProcess[str]) -> str:
    text = "\n".join(part for part in ((result.stdout or "").strip(), (result.stderr or "").strip()) if part)
    return "\n".join(text.splitlines()[-30:]) or f"exit code {result.returncode}"


def _method_text(version: str | None, params: Mapping[str, Any], fixed_model: bool) -> str:
    tool = f"IQ-TREE {version}" if version else "IQ-TREE"
    model = "ModelFinder (BIC)" if not fixed_model else f"model {params['model']}"
    if params["support_type"] == "ufboot":
        return f"{tool}, {model}, ultrafast bootstrap"
    if params["support_type"] == "sh_alrt":
        return f"{tool} fast search, {model}, SH-aLRT"
    return f"{tool}, {model}, no support"


def tree_stamp_path(output_path: Path, tree_id: str) -> Path:
    return Path(output_path) / "trees" / f"{tree_id}.json"


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    os.replace(tmp, path)


def run_core_tree(
    output_path: Path,
    settings: Mapping[str, Any],
    logger: logging.Logger,
    *,
    threads: int,
) -> CoreTreeResult:
    output_path = Path(output_path).expanduser().resolve()
    plan = plan_core_tree(output_path, settings)
    params = _stamp_params(plan, settings, threads=threads)
    ordering = read_ordering_summary(output_path) or {}
    gubbins = read_gubbins_summary(output_path) or {}
    version = iqtree_version()
    warnings: list[dict[str, str]] = []
    notes: list[str] = []
    final_dir = output_path / IQTREE_DIR
    final_dir.parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=".core_tree.", dir=final_dir.parent))
    prefix = work / plan.tree_id
    taxa_names, sites, variable_sites = _alignment_shape(plan.input_path)
    parsed: dict[str, Any] = {}
    try:
        raw_tree = prefix.with_suffix(".treefile")
        if plan.taxa <= 2:
            _write_trivial_tree(taxa_names, raw_tree)
            notes.append(f"Trivial tree: {plan.taxa} genome(s).")
        elif not variable_sites:
            _write_trivial_tree(taxa_names, raw_tree)
            notes.append("no_variable_sites: star tree, no IQ-TREE call.")
        else:
            cmd = [
                *tool_cmd("iqtree2", workdir=str(work)),
                "-s", str(plan.input_path),
                "-m", plan.model,
                "--seed", str(params["seed"]),
                "-T", str(params["threads"]),
                "--prefix", str(prefix),
                "--quiet",
            ]
            if params["support_type"] == "ufboot":
                cmd += ["-B", str(plan.bootstraps)]
                if params["bnni"]:
                    cmd.append("--bnni")
            elif params["support_type"] == "sh_alrt":
                cmd += ["--fast", "--alrt", str(plan.bootstraps)]
            logger.info(
                "Running IQ-TREE %s on %s (%d genomes, %s sites, model %s, support %s, %d threads).",
                version or "",
                plan.input_path.name,
                plan.taxa,
                sites,
                plan.model,
                params["support_type"],
                params["threads"],
            )
            logger.info("CMD: " + " ".join(cmd))
            result = subprocess.run(cmd, capture_output=True, text=True, check=False)
            if result.returncode != 0 or not raw_tree.is_file():
                log_tail = ""
                log_path = Path(str(prefix) + ".log")
                if log_path.is_file():
                    log_tail = "\n".join(log_path.read_text(errors="ignore").splitlines()[-30:])
                raise RuntimeError(f"IQ-TREE failed: {log_tail or _format_error(result)}")
            parsed = parse_iqtree_outputs(prefix)
        if plan.support_type == "none" or plan.taxa < MIN_GENOMES_FOR_SUPPORT:
            warnings.append(
                _warning(
                    "core_tree_no_support_few_genomes",
                    f"The core tree has no branch support: {plan.taxa} genome(s), bootstrap needs at least "
                    f"{MIN_GENOMES_FOR_SUPPORT}.",
                )
            )
        if parsed and params["support_type"] == "ufboot" and not parsed.get("ufboot_converged", True):
            warnings.append(
                _warning(
                    "core_tree_ufboot_not_converged",
                    "IQ-TREE reports that the ultrafast bootstrap did not converge; support values are less "
                    "reliable.",
                )
            )
        if parsed.get("composition_failed"):
            notes.append(
                f"{len(parsed['composition_failed'])} sequence(s) failed IQ-TREE's composition chi2 test: "
                + ", ".join(parsed["composition_failed"][:20])
            )
        identical = int(parsed.get("identical_kept") or 0) + int(parsed.get("identical_ignored") or 0)
        if identical:
            notes.append(
                f"{identical} sequence(s) are identical to another one; support among identical sequences "
                "is meaningless."
            )
        # A trivial or star tree written without IQ-TREE carries no support.
        built = plan.taxa > 2 and bool(variable_sites)
        support_type = plan.support_type if built else "none"
        tree = read_newick(raw_tree)
        supports = support_values(tree)
        iqtree_notes = list(notes)
        rf = None
        if plan.tree_id == "core_recomb":
            rf, gubbins_notes = _gubbins_notes(output_path, gubbins, tree)
            notes.extend(gubbins_notes)
        for item in warnings:
            logger.warning(item["message"])

        for path in list(work.iterdir()):
            if path.is_file() and not any(path.name.endswith(s) for s in KEPT_SUFFIXES):
                path.unlink()
        stamp_run = {
            "algorithm_version": ALGORITHM_VERSION,
            "tool_version": version,
            "input_path": str(plan.input_path.relative_to(output_path)),
            "input_sha256": plan.input_sha256,
            "params": params,
            "sites": sites,
            "variable_sites": variable_sites,
            **{k: parsed.get(k) for k in ("model", "model_criterion", "log_likelihood", "ufboot_converged")},
            "composition_failed": parsed.get("composition_failed") or [],
            "identical_collapsed": identical,
            "support_median": float(np.median(supports)) if supports else None,
            "support_share_at_threshold": (
                round(sum(1 for v in supports if v >= SUPPORT_THRESHOLDS[support_type]) / len(supports), 4)
                if supports and support_type in SUPPORT_THRESHOLDS else None
            ),
            "rf_to_gubbins": rf,
            "warnings": warnings,
            "iqtree_notes": iqtree_notes,
            "notes": notes,
        }
        _write_json(work / RUN_STAMP, stamp_run)
        if final_dir.exists():
            shutil.rmtree(final_dir)
        os.replace(work, final_dir)
    finally:
        if work.exists():
            shutil.rmtree(work, ignore_errors=True)

    raw_final = final_dir / f"{plan.tree_id}.treefile"
    fixed_model = not str(settings["iqtree_model"]).upper().startswith("MF")
    stamp = {
        "id": plan.tree_id,
        "label": TREE_LABELS[plan.tree_id],
        "raw_file": str(raw_final.relative_to(output_path)),
        "input_path": str(plan.input_path.relative_to(output_path)),
        "input_sha256": stamp_run["input_sha256"],
        "algorithm_version": ALGORITHM_VERSION,
        "params": params,
        "tool_version": version,
        "threads": params["threads"],
        "seed": params["seed"],
        "method": (
            _method_text(version, params, fixed_model)
            if built
            else ("trivial tree" if plan.taxa <= 2 else "star tree (no variable sites)")
        ),
        "alignment": str(plan.input_path.relative_to(output_path)),
        "alignment_kind": "core_alignment_masked" if plan.tree_id == "core_recomb" else "core_alignment",
        "sites": sites,
        "invariant_correction": None,
        "model": parsed.get("model"),
        "log_likelihood": parsed.get("log_likelihood"),
        "support": {
            "type": support_type,
            "replicates": plan.bootstraps if support_type != "none" else 0,
            "threshold": SUPPORT_THRESHOLDS.get(support_type),
        },
        "recombination_filtered": plan.tree_id == "core_recomb",
        "ordering_genome": ordering.get("ordering_genome"),
        "branch_length_unit": "substitutions/site",
        # Every genome is in the masked alignment; the ones Gubbins did not mask are listed.
        "samples_missing": [],
        "samples_unmasked": _unmasked_names(gubbins) if plan.tree_id == "core_recomb" else [],
        "reference_sample": None,
        "fast_mode": plan.fast_mode,
        "notes": notes,
    }
    for other in ("core_recomb", "core"):
        if other != plan.tree_id:
            for suffix in (".json", ".nwk"):
                (output_path / "trees" / f"{other}{suffix}").unlink(missing_ok=True)
    _write_json(tree_stamp_path(output_path, plan.tree_id), stamp)
    logger.info(
        "Core tree %s: model %s, log-likelihood %s, %d genomes, support %s (median %s).",
        plan.tree_id,
        stamp["model"],
        stamp["log_likelihood"],
        len(leaf_names(read_newick(raw_final))),
        support_type,
        stamp_run["support_median"],
    )
    return CoreTreeResult(tree_id=plan.tree_id, raw_tree=raw_final, stamp=stamp, warnings=warnings)


def read_core_tree_stamp(output_path: Path) -> dict[str, Any] | None:
    try:
        data = json.loads((Path(output_path) / IQTREE_DIR / RUN_STAMP).read_text(encoding="utf-8"))
    except Exception:
        return None
    return data if isinstance(data, dict) else None


def core_tree_reuse_issue(output_path: Path, settings: Mapping[str, Any], *, threads: int) -> str | None:
    """Why the existing core tree cannot be reused (None when it can)."""
    output_path = Path(output_path)
    stamp = read_core_tree_stamp(output_path)
    if stamp is None:
        return "no input stamp"
    if stamp.get("algorithm_version") != ALGORITHM_VERSION:
        return "the core tree step changed"
    try:
        plan = plan_core_tree(output_path, settings)
    except RuntimeError as exc:
        return str(exc)
    # The thread count is recorded but does not make a finished tree stale.
    current = _stamp_params(plan, settings, threads=threads)
    recorded = {k: v for k, v in (stamp.get("params") or {}).items() if k != "threads"}
    if recorded != {k: v for k, v in current.items() if k != "threads"}:
        return "parameters changed"
    # The masked alignment's sha256 is the one recorded (and checked) by the Gubbins step.
    if stamp.get("input_sha256") != plan.input_sha256:
        return "input alignment changed"
    if stamp.get("tool_version") != iqtree_version():
        return "IQ-TREE version changed"
    raw = output_path / IQTREE_DIR / f"{plan.tree_id}.treefile"
    if not raw.is_file() or not tree_stamp_path(output_path, plan.tree_id).is_file():
        return "tree files are missing"
    return None


# ------------------------------------------------ branch lengths including recombination
REFIT_DIR = Path("trees") / "core_recomb_refit"
REFIT_TREE_ID = "core_recomb_with_recombination"
REFIT_ALGORITHM_VERSION = 1
REFIT_KEPT_SUFFIXES = (".iqtree", ".log", ".treefile")
BRANCH_LENGTH_SET_LABELS = {"mutations": "Mutations only", "with_recombination": "Including recombination"}
_REFIT_SKIP_TEXT = {
    "disabled": "phylo.recombination_branch_lengths is off",
    "gubbins_skipped": "Gubbins did not run, so the core tree already includes recombination",
    "nothing_masked": "Gubbins masked no base, so both sets of branch lengths would be the same",
    "no_iqtree_tree": "the core tree was not built by IQ-TREE (too few genomes or no variable sites)",
}


def read_refit_stamp(output_path: Path) -> dict[str, Any] | None:
    try:
        data = json.loads((Path(output_path) / REFIT_DIR / RUN_STAMP).read_text(encoding="utf-8"))
    except Exception:
        return None
    return data if isinstance(data, dict) else None


def current_refit_stamp(output_path: Path, raw_tree: Path) -> dict[str, Any] | None:
    """The built refit stamp when it is the refit of ``raw_tree`` (the raw core_recomb
    tree) on the current unmasked alignment; None otherwise (none, skipped, failed, or
    a refit kept from before the core tree or the alignment changed)."""
    output_path = Path(output_path)
    refit = read_refit_stamp(output_path) or {}
    current = (
        refit.get("status") == "built"
        and (output_path / REFIT_DIR / f"{REFIT_TREE_ID}.treefile").is_file()
        and Path(raw_tree).is_file()
        and refit.get("input_tree_sha256") == file_sha256(raw_tree)
        and refit.get("alignment_sha256") == ordered_alignment_sha256(output_path)
    )
    return refit if current else None


def _masked_base_total(output_path: Path) -> int | None:
    """A/C/G/T bases set to N by the Gubbins masking, over every genome (None when the
    table is missing)."""
    try:
        lines = (output_path / "gubbins" / MASKED_BASES_TSV).read_text(encoding="utf-8").splitlines()
    except OSError:
        return None
    header = lines[0].split("\t") if lines else []
    if "masked_bases" not in header:
        return None
    column = header.index("masked_bases")
    return sum(int(row.split("\t")[column] or 0) for row in lines[1:] if row.strip())


def _refit_inputs(output_path: Path, settings: Mapping[str, Any]) -> tuple[dict[str, Any] | None, str | None]:
    """(what the refit reads, None), or (None, why it is skipped)."""
    if not settings.get("recombination_branch_lengths", True):
        return None, "disabled"
    core = read_core_tree_stamp(output_path) or {}
    if (core.get("params") or {}).get("tree_id") != "core_recomb":
        return None, "gubbins_skipped"
    raw = output_path / IQTREE_DIR / "core_recomb.treefile"
    if not core.get("model") or not raw.is_file():
        return None, "no_iqtree_tree"
    ordered_sha = ordered_alignment_sha256(output_path)
    alignment = output_path / ORDERED_ALIGNMENT_REL
    if not ordered_sha or not alignment.is_file():
        raise RuntimeError(f"The unmasked core alignment is missing: {alignment}")
    if _masked_base_total(output_path) == 0 or masked_alignment_sha256(read_gubbins_summary(output_path)) == ordered_sha:
        return None, "nothing_masked"
    return {
        "input_tree": str(raw.relative_to(output_path)),
        "input_tree_sha256": file_sha256(raw),
        "alignment": str(ORDERED_ALIGNMENT_REL),
        # Recorded and checked against the file by the core_alignment step.
        "alignment_sha256": ordered_sha,
        "model": str(core["model"]),
        "seed": int(settings["seed"]),
    }, None


def _refit_reuse_key(inputs: Mapping[str, Any], version: str | None) -> dict[str, Any]:
    return {
        **{k: inputs[k] for k in ("input_tree_sha256", "alignment_sha256", "model", "seed")},
        "tool_version": version,
        "algorithm_version": REFIT_ALGORITHM_VERSION,
    }


def _branch_length_change(mutations, with_recombination) -> dict[str, Any]:
    """Tree lengths and the edges that grew most, matched by split."""
    before = split_lengths(mutations)
    after = split_lengths(with_recombination)
    ratios = sorted(after[k] / before[k] for k in before if before[k] > 1e-5)
    grown = sorted(before, key=lambda k: after[k] - before[k], reverse=True)[:5]
    leaves = frozenset(leaf_names(mutations))
    smaller = {k: (k if 2 * len(k) <= len(leaves) else leaves - k) for k in grown}
    total_before, total_after = sum(before.values()), sum(after.values())
    return {
        "tree_length": {"mutations": round(total_before, 6), "with_recombination": round(total_after, 6)},
        "tree_length_ratio": round(total_after / total_before, 4) if total_before > 0 else None,
        "median_edge_ratio": round(ratios[len(ratios) // 2], 4) if ratios else None,
        "largest_increases": [
            {
                "genomes": len(smaller[k]),
                "leaf": next(iter(smaller[k])) if len(smaller[k]) == 1 else None,
                "mutations": round(before[k], 6),
                "with_recombination": round(after[k], 6),
            }
            for k in grown
        ],
    }


def _replace_refit_dir(final_dir: Path, record: Mapping[str, Any]) -> None:
    """Leave only the record (a skip or a failure) in the refit folder."""
    if final_dir.exists():
        shutil.rmtree(final_dir)
    _write_json(final_dir / RUN_STAMP, record)


def _refit_failed(final_dir: Path, stamp: dict[str, Any], exc: Exception, logger: logging.Logger) -> dict[str, Any]:
    """Record a failed refit. A refit built before is kept (trees.json lists it only
    while it still matches the current tree and alignment); a later run retries."""
    message = f"Branch lengths including recombination could not be computed: {exc}"
    logger.warning(message)
    stamp.update({"status": "failed", "error": str(exc), "warnings": [_warning("core_tree_refit_failed", message)]})
    try:
        previous = json.loads((final_dir / RUN_STAMP).read_text(encoding="utf-8"))
    except Exception:
        previous = None
    kept = isinstance(previous, dict) and previous.get("status") == "built"
    if not (kept and (final_dir / f"{REFIT_TREE_ID}.treefile").is_file()):
        try:
            _replace_refit_dir(final_dir, stamp)
        except OSError as record_exc:
            logger.warning("The failed refit could not be recorded: %s", record_exc)
    return stamp


def run_recombination_refit(
    output_path: Path,
    settings: Mapping[str, Any],
    logger: logging.Logger,
    *,
    threads: int,
) -> dict[str, Any]:
    """Refit the core_recomb branch lengths on the unmasked ordered alignment (topology
    fixed, ``-te``) into ``trees/core_recomb_refit/``, or reuse the last refit, or
    record why there is none. Returns the refit stamp (``status``: built, reused,
    skipped or failed). A failure is a warning, never a failed core tree."""
    output_path = Path(output_path).expanduser().resolve()
    final_dir = output_path / REFIT_DIR
    stamp: dict[str, Any] = {"algorithm_version": REFIT_ALGORITHM_VERSION}
    try:
        return _recombination_refit(output_path, final_dir, settings, logger, threads=threads, stamp=stamp)
    except Exception as exc:
        return _refit_failed(final_dir, stamp, exc, logger)


def _recombination_refit(
    output_path: Path,
    final_dir: Path,
    settings: Mapping[str, Any],
    logger: logging.Logger,
    *,
    threads: int,
    stamp: dict[str, Any],
) -> dict[str, Any]:
    """``run_recombination_refit``'s body; ``stamp`` collects what a failure records."""
    inputs, skip = _refit_inputs(output_path, settings)
    if skip or inputs is None:
        record = {
            "algorithm_version": REFIT_ALGORITHM_VERSION,
            "status": "skipped",
            "skip_reason": skip,
            "message": _REFIT_SKIP_TEXT[skip],
            "warnings": [],
        }
        _replace_refit_dir(final_dir, record)
        logger.info("Branch lengths including recombination: not computed (%s).", record["message"])
        return record
    version = iqtree_version()
    key = _refit_reuse_key(inputs, version)
    previous = read_refit_stamp(output_path) or {}
    previous_key = dict(previous.get("key") or {})
    if version is None:
        # IQ-TREE's version could not be read (container briefly unavailable): match
        # the rest of the key rather than rerun and lose a good refit.
        previous_key["tool_version"] = None
    if (
        previous.get("status") == "built"
        and previous_key == key
        and (final_dir / f"{REFIT_TREE_ID}.treefile").is_file()
    ):
        logger.info("Branch lengths including recombination: reused (same tree, alignment, model and IQ-TREE).")
        return {**previous, "status": "reused"}

    threads = int(max(1, threads))
    final_dir.parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=".core_recomb_refit.", dir=final_dir.parent))
    prefix = work / REFIT_TREE_ID
    cmd = [
        *tool_cmd("iqtree2", workdir=str(work)),
        "-s", str(output_path / inputs["alignment"]),
        "-te", str(output_path / inputs["input_tree"]),
        "-m", inputs["model"],
        # Identical sequences stay in: the topology is fixed, IQ-TREE must not regroup them.
        "-keep-ident",
        "--seed", str(inputs["seed"]),
        "-T", str(threads),
        "--prefix", str(prefix),
        "--quiet",
    ]
    stamp.update({
        "algorithm_version": REFIT_ALGORITHM_VERSION,
        "key": key,
        "tool_version": version,
        **inputs,
        "threads": threads,
        "command": cmd,
        "warnings": [],
    })
    logger.info(
        "Refitting the core_recomb branch lengths on the unmasked alignment (IQ-TREE %s, model %s, fixed "
        "topology, %d threads).",
        version or "",
        inputs["model"],
        threads,
    )
    logger.info("CMD: " + " ".join(cmd))
    started = time.monotonic()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        raw = prefix.with_suffix(".treefile")
        if result.returncode != 0 or not raw.is_file():
            log_path = Path(str(prefix) + ".log")
            log_tail = (
                "\n".join(log_path.read_text(errors="ignore").splitlines()[-30:]) if log_path.is_file() else ""
            )
            raise RuntimeError(f"IQ-TREE failed: {log_tail or _format_error(result)}")
        source = read_newick(output_path / inputs["input_tree"])
        refit = read_newick(raw)
        if set(leaf_names(refit)) != set(leaf_names(source)) or set(split_supports(refit)) != set(
            split_supports(source)
        ):
            raise RuntimeError("IQ-TREE returned a different topology than the core tree it was given.")
        parsed = parse_iqtree_outputs(prefix)
        for path in list(work.iterdir()):
            if path.is_file() and not any(path.name.endswith(suffix) for suffix in REFIT_KEPT_SUFFIXES):
                path.unlink()
        stamp.update(
            {
                "status": "built",
                "wall_time_s": round(time.monotonic() - started, 1),
                "log_likelihood": parsed.get("log_likelihood"),
                "branch_length_unit": "substitutions/site",
                **_branch_length_change(source, refit),
            }
        )
        _write_json(work / RUN_STAMP, stamp)
        if final_dir.exists():
            shutil.rmtree(final_dir)
        os.replace(work, final_dir)
    except Exception:
        stamp["wall_time_s"] = round(time.monotonic() - started, 1)
        raise
    finally:
        if work.exists():
            shutil.rmtree(work, ignore_errors=True)
    logger.info(
        "Branch lengths including recombination: tree length %.4g -> %.4g substitutions/site (x%s) in %.0f s.",
        stamp["tree_length"]["mutations"],
        stamp["tree_length"]["with_recombination"],
        stamp["tree_length_ratio"],
        stamp["wall_time_s"],
    )
    return stamp


__all__ = [
    "BRANCH_LENGTH_SET_LABELS",
    "CoreTreeResult",
    "IQTREE_DIR",
    "REFIT_DIR",
    "REFIT_TREE_ID",
    "TREE_LABELS",
    "core_tree_reuse_issue",
    "parse_iqtree_outputs",
    "plan_core_tree",
    "current_refit_stamp",
    "read_core_tree_stamp",
    "read_refit_stamp",
    "refresh_gubbins_notes",
    "run_core_tree",
    "run_recombination_refit",
    "tree_stamp_path",
]
