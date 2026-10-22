"""FastBAPS lineages and clusters on the ordered core alignment (step ``fastgear_clades``
with phylo.population_method fastbaps; see ``pipeline.PopulationStructure``).

FastBAPS (Tonkin-Hill et al. 2019, R package fastbaps in the image's ``fastbaps`` env)
reads ``core_alignment/core_gene_alignment.ordered.aln`` (the file fastGEAR, Gubbins and
the core tree read), fits the Dirichlet prior (``optimise_prior``, type
``optimise.symmetric``) and clusters hierarchically at two levels (``multi_res_baps``,
levels 2). The run's ``sample_clades.tsv`` maps them as:

- lineage (FastGearLineage, the groups Gubbins runs per) = Level 1,
- cluster (FastGearCluster) = Level 2, nested in Level 1 (cluster ids are unique across
  lineages, as fastGEAR's are),
- strain index = the genome's row in the alignment (1 based), method = ``fastbaps``.

On the ST492 reference cohort Level 1 matches fastGEAR's clusters (one of them split in
two) and isolates fastGEAR's lineage 1 exactly, so a FastBAPS run has more, smaller
lineages than a fastGEAR run of the same genomes.

Outputs in ``fastbaps/``: ``run_fastbaps.R`` (the script that ran), ``fastbaps_clusters.tsv``
(FastBAPS' own table: Isolates, Level 1, Level 2), ``sample_clades.tsv``, ``fastbaps.log``
and ``input.json`` (the input stamp: alignment sha256, parameters, script, image; the
step reruns when any of them changes).
"""
from __future__ import annotations

import hashlib
import json
import logging
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import pandas as pd

from pipeline.CoreAlignment import ORDERED_ALIGNMENT_REL, read_ordering_summary
from pipeline.FastGearWorkflow import _image_identity, validate_sample_clades_tsv
from pipeline.PopulationStructure import SAMPLE_CLADES_NAME, method_dir

ALGORITHM_VERSION = 1
INPUT_STAMP = "input.json"
CLUSTERS_TSV = "fastbaps_clusters.tsv"
SCRIPT_NAME = "run_fastbaps.R"
LOG_NAME = "fastbaps.log"
PARAMS = {"prior": "optimise.symmetric", "levels": 2}

# The R script, written into fastbaps/ for each run (its sha256 is in the stamp).
R_SCRIPT = """\
# FastBAPS on the ordered core alignment (written by BACTolomew, pipeline/FastBapsWorkflow.py).
# Usage: Rscript run_fastbaps.R <alignment.fasta> <clusters.tsv> <seed>
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages({library(fastbaps)})
set.seed(as.integer(args[3]))
sparse <- import_fasta_sparse_nt(args[1])
sparse <- optimise_prior(sparse, type = "optimise.symmetric")
clusters <- multi_res_baps(sparse, levels = 2)
write.table(clusters, args[2], sep = "\\t", quote = FALSE, row.names = FALSE)
cat("fastbaps", as.character(packageVersion("fastbaps")), "\\n")
"""


@dataclass
class FastBapsResult:
    sample_clades_tsv: Path
    clusters_tsv: Path


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _fastbaps_stamp(output_path: Path, seed: int, ctx) -> dict[str, Any]:
    ordering = read_ordering_summary(output_path) or {}
    return {
        "algorithm_version": ALGORITHM_VERSION,
        "method": "fastbaps",
        "input_path": str(Path(output_path) / ORDERED_ALIGNMENT_REL),
        "input_sha256": (ordering.get("output") or {}).get("sha256"),
        "params": {**PARAMS, "seed": int(seed)},
        "script_sha256": _sha256_text(R_SCRIPT),
        "image": _image_identity(ctx),
    }


def fastbaps_reuse_issue(output_path: Path, seed: int) -> str | None:
    """Why the existing FastBAPS output cannot be reused (None when it can)."""
    from pipeline.ToolEnv import current_container_context, make_container_context

    output_path = Path(output_path)
    try:
        recorded = json.loads((method_dir(output_path, "fastbaps") / INPUT_STAMP).read_text(encoding="utf-8"))
    except Exception:
        return "no input stamp"
    ctx = current_container_context() or make_container_context(output_dir=output_path)
    current = _fastbaps_stamp(output_path, seed, ctx)
    if recorded.get("input_sha256") != current["input_sha256"] or not current["input_sha256"]:
        return "input alignment changed"
    for key in ("algorithm_version", "params", "script_sha256", "image"):
        if recorded.get(key) != current[key]:
            return f"{key.replace('_', ' ')} changed"
    return None


def alignment_names(path: Path) -> list[str]:
    """The record names of a FASTA file, in file order."""
    names: list[str] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                names.append(line[1:].strip().split()[0] if line[1:].strip() else "")
    return names


def sample_clades_from_clusters(clusters_tsv: Path, alignment_order: Sequence[str]) -> pd.DataFrame:
    """FastBAPS' table (Isolates, Level 1, Level 2) as sample_clades.tsv rows: lineage =
    Level 1, cluster = Level 2, strain index = the row in the alignment (1 based)."""
    df = pd.read_csv(clusters_tsv, sep="\t", dtype=str)
    columns = {str(c).strip().lower().replace(".", " "): c for c in df.columns}
    try:
        sample_col, level1, level2 = columns["isolates"], columns["level 1"], columns["level 2"]
    except KeyError:
        raise RuntimeError(
            f"{clusters_tsv} has no Isolates, Level 1 and Level 2 columns (found {', '.join(map(str, df.columns))})"
        ) from None
    index = {name: i + 1 for i, name in enumerate(alignment_order)}
    out = pd.DataFrame({
        "Sample": df[sample_col].astype(str).str.strip(),
        "FastGearLineage": df[level1].astype(str).str.strip(),
        "FastGearCluster": df[level2].astype(str).str.strip(),
    })
    out["FastGearStrainIndex"] = [index.get(name, "") for name in out["Sample"]]
    out["FastGearMethod"] = "fastbaps"
    order = pd.to_numeric(out["FastGearStrainIndex"], errors="coerce")
    return out.assign(_order=order).sort_values(["_order", "Sample"], kind="stable").drop(columns="_order")


def run_fastbaps_clades(
    output_path: Path,
    *,
    expected_samples: Sequence[str],
    logger: logging.Logger,
    seed: int = 1,
) -> FastBapsResult:
    from pipeline.ToolEnv import current_container_context, make_container_context, tool_cmd

    output_path = Path(output_path)
    alignment = output_path / ORDERED_ALIGNMENT_REL
    if not alignment.is_file():
        raise RuntimeError("The ordered core alignment (step core_alignment) is required before FastBAPS.")
    out_dir = method_dir(output_path, "fastbaps")
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)
    script = out_dir / SCRIPT_NAME
    script.write_text(R_SCRIPT, encoding="utf-8")
    clusters_tsv = out_dir / CLUSTERS_TSV
    partial = out_dir / f".{CLUSTERS_TSV}.partial"

    cmd = [*tool_cmd("fastbaps", "Rscript", workdir=str(out_dir)), str(script), str(alignment), str(partial), str(int(seed))]
    logger.info("Running FastBAPS.")
    logger.info("CMD: " + " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    (out_dir / LOG_NAME).write_text((result.stdout or "") + (result.stderr or ""), encoding="utf-8")
    if result.returncode != 0 or not partial.is_file():
        details = "\n".join(line for line in ((result.stderr or "") + (result.stdout or "")).strip().splitlines()[-5:])
        partial.unlink(missing_ok=True)
        raise RuntimeError(
            f"FastBAPS failed with exit code {result.returncode}" + (f": {details}" if details else ".")
        )
    os.replace(partial, clusters_tsv)

    table = sample_clades_from_clusters(clusters_tsv, alignment_names(alignment))
    sample_clades_tsv = out_dir / SAMPLE_CLADES_NAME
    table.to_csv(sample_clades_tsv, sep="\t", index=False)
    # The same checks as a fastGEAR result (columns, one row per genome of the run).
    validate_sample_clades_tsv(sample_clades_tsv, expected_samples=expected_samples)
    ctx = current_container_context() or make_container_context(output_dir=output_path)
    (out_dir / INPUT_STAMP).write_text(json.dumps(_fastbaps_stamp(output_path, seed, ctx), indent=2), encoding="utf-8")
    lineages = table["FastGearLineage"].nunique()
    clusters = table["FastGearCluster"].nunique()
    logger.info(f"FastBAPS: {lineages} lineage(s) (Level 1), {clusters} cluster(s) (Level 2): {sample_clades_tsv}")
    return FastBapsResult(sample_clades_tsv=sample_clades_tsv, clusters_tsv=clusters_tsv)


__all__ = [
    "CLUSTERS_TSV",
    "FastBapsResult",
    "fastbaps_reuse_issue",
    "run_fastbaps_clades",
    "sample_clades_from_clusters",
]
