from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from pipeline.Utils import CATEGORIES, cat_priority

try:
    CATEGORY_PRIORITY_MAP = dict(cat_priority)
except Exception:
    CATEGORY_PRIORITY_MAP = {}

if not CATEGORY_PRIORITY_MAP:
    CATEGORY_PRIORITY_MAP = {cat: i for i, cat in enumerate(CATEGORIES)}
if "OTHER" not in CATEGORY_PRIORITY_MAP:
    CATEGORY_PRIORITY_MAP["OTHER"] = max(CATEGORY_PRIORITY_MAP.values(), default=0) + 1


def infer_annotation_category_and_source(markers: list[dict]) -> tuple[str, Any, Any]:
    """
    Given a list of marker dicts, return:
      (ConsensusCategory, ConsensusSourceMethod, ConsensusDonorProduct)
    """
    if not markers:
        return "OTHER", pd.NA, pd.NA

    best_cat = None
    best_rank = float("inf")
    best_source = pd.NA
    best_donor = pd.NA

    for m in markers:
        cat = m.get("Category") or m.get("MarkerCategory")
        if not cat:
            continue

        rank = CATEGORY_PRIORITY_MAP.get(cat, CATEGORY_PRIORITY_MAP.get("OTHER", float("inf")))
        if rank < best_rank:
            best_rank = rank
            best_cat = cat
            best_source = m.get("Method")
            best_donor = m.get("Product")

    if best_cat is None:
        return "OTHER", pd.NA, pd.NA

    return best_cat, best_source, best_donor


def assign_annotation_ids(df: pd.DataFrame, prefix: str = "ANN") -> pd.DataFrame:
    df = df.reset_index(drop=True).copy()
    df["AnnotationID"] = [f"{prefix}{i:07d}" for i in range(1, len(df) + 1)]
    return df


def write_annotation_fasta(
    df: pd.DataFrame,
    output_path: Path,
    logger: logging.Logger,
    id_col: str = "AnnotationID",
    seq_col: str = "Nucleotide_Sequence",
    filename: str = "annotations_nt.fasta",
) -> Path | None:
    if df.empty or id_col not in df.columns or seq_col not in df.columns:
        logger.warning("No annotations with sequences to export to FASTA.")
        return None

    df = df.dropna(subset=[id_col, seq_col])
    if df.empty:
        logger.warning("All annotations lack nucleotide sequences; FASTA not written.")
        return None

    fa_path = output_path / filename
    with fa_path.open("w", encoding="utf-8") as fh:
        for _, row in df.iterrows():
            aid = str(row[id_col])
            seq = str(row[seq_col]).replace(" ", "").replace("\n", "")
            if not seq:
                continue
            fh.write(f">{aid}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")

    logger.info(f"Wrote nucleotide FASTA with annotations: {fa_path}")
    return fa_path


def add_compound_clusters(comp_df: pd.DataFrame, clusters_df: pd.DataFrame) -> pd.DataFrame:
    """
    For each compound annotation, find all overlapping gene/MarkerSequence clusters
    (same Sample+Contig, any interval overlap) and write them as a comma-separated string.
    """
    comp_df = comp_df.copy()

    if comp_df.empty or clusters_df.empty:
        comp_df["ClusterID"] = pd.NA
        return comp_df

    cdf = clusters_df[["Sample", "Contig", "Start", "Stop", "ClusterID"]].dropna(subset=["Start", "Stop"])
    cdf["Start"] = pd.to_numeric(cdf["Start"], errors="coerce")
    cdf["Stop"] = pd.to_numeric(cdf["Stop"], errors="coerce")
    cdf = cdf.dropna(subset=["Start", "Stop"])
    cdf["Start"] = cdf["Start"].astype(int)
    cdf["Stop"] = cdf["Stop"].astype(int)

    comp_df["Start"] = pd.to_numeric(comp_df["Start"], errors="coerce")
    comp_df["Stop"] = pd.to_numeric(comp_df["Stop"], errors="coerce")
    comp_df.loc[comp_df["Start"] > comp_df["Stop"], ["Start", "Stop"]] = comp_df.loc[
        comp_df["Start"] > comp_df["Stop"], ["Stop", "Start"]
    ].values
    comp_df["ClusterID"] = ""

    for (samp, cont), gcomp in comp_df.groupby(["Sample", "Contig"], dropna=False):
        cg = cdf[(cdf["Sample"] == samp) & (cdf["Contig"] == cont)]
        if cg.empty:
            continue

        starts = cg["Start"].to_numpy()
        stops = cg["Stop"].to_numpy()
        cids = cg["ClusterID"].astype(str).to_numpy()

        for idx, row in gcomp.iterrows():
            s = row["Start"]
            e = row["Stop"]
            if pd.isna(s) or pd.isna(e):
                continue
            s = int(s)
            e = int(e)
            mask = (starts <= e) & (stops >= s)
            if mask.any():
                comp_df.at[idx, "ClusterID"] = ",".join(sorted(set(cids[mask].tolist())))

    comp_df.loc[comp_df["ClusterID"] == "", "ClusterID"] = pd.NA
    return comp_df


def compute_category(row) -> str:
    product = row["Product"]

    if not pd.isna(product):
        if "virulence" in product.lower():
            return "VIRULENCE"
        if "transposase" in product.lower() or "integrase" in product.lower():
            return "MOBILE"
        if "replicon" in product.lower() or "plasmid" in product.lower():
            return "PLASMID"
        if "phage" in product.lower():
            return "PHAGE"
        if "resistance" in product.lower():
            return "AMR"
        if "efflux" in product.lower():
            return "EFFLUX"

    return None


def compute_gene(annotation: str) -> str:
    try:
        return annotation.split("|")[1].strip()
    except Exception:
        return annotation


__all__ = [
    "CATEGORY_PRIORITY_MAP",
    "infer_annotation_category_and_source",
    "assign_annotation_ids",
    "write_annotation_fasta",
    "add_compound_clusters",
    "compute_category",
    "compute_gene",
]
