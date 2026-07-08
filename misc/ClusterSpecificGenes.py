import numpy as np
import pandas as pd
from typing import Tuple, Dict, Optional
from pathlib import Path
import re

def compute_cluster_specific_from_clusters(
    gpa_df: pd.DataFrame,
    clusters_by_method: Dict[str, Dict[str, str] | None],
    *,
    gene_col: Optional[str] = None,
    name_cols: tuple[str, ...] = ("Non-unique Gene name", "Annotation"),
    min_cluster_size: int = 2,
    drop_never_present: bool = True,         # drop genes with zero copies across considered samples (per method)
    drop_globally_core: bool = True,         # drop genes present in ALL genomes (≥1 copy in every sample)
    add_fisher_from_presence: bool = True,   # add p/q from binary presence (per-genome)
    alpha: float = 0.05,
    drop_insignificant: bool = True,         # if True and p/q available, drop non-significant rows
    use_fdr: bool = True,                    # BH-FDR if p-values are computed
    drop_unspecific: bool = True,            # drop genes with low specificity by copies
    unspecific_threshold: float = 0.5,       # threshold for copy_specificity (keep >= this)
    drop_hypothetical: bool = True,          # drop annotation == 'hypothetical protein'
    out_dir: Optional[str | Path] = None,    # if provided, save per-cluster CSVs here (e.g., 'out/')
    save_empty: bool = False,                # save even empty tables
) -> Tuple[Dict[str, Dict[str, pd.DataFrame]], pd.DataFrame]:
    """
    Count *occurrences* (copies) of each gene inside vs outside each cluster,
    for every clustering method in `clusters_by_method`.

    Returns
    -------
    by_method : dict[method][cluster] -> DataFrame with columns:
        ['method','cluster','gene','gene_id','gene_name','annotation',
         'in_copies','out_copies','total_copies',
         'is_unique_by_copy','copy_share_in','copy_specificity',
         (optional) 'p_value','q_value']
    long_df   : concatenated DataFrame of all clusters/methods.
    """

    # --- detect gene + sample columns ---
    cols = gpa_df.columns.astype(str).tolist()
    if gene_col is None:
        gene_col = "Gene" if "Gene" in cols else cols[0]

    meta_like = {
        gene_col, "Gene ID", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences",
        "Avg sequences per isolate", "Genome Fragment", "Order within Fragment", "Accessory Fragment",
        "Accessory Order with Fragment", "QC", "Min group size nuc", "Max group size nuc", "Avg group size nuc"
    }
    sample_cols = [c for c in cols if c not in meta_like]
    if len(sample_cols) < 2:
        raise ValueError("Fewer than 2 sample columns detected in gpa_df.")

    # --- parse per-sample cells to (a) copy counts and (b) binary presence ---
    def _count_entries(cell) -> int:
        if cell is None:
            return 0
        s = str(cell).strip()
        if s == "" or s == "0" or s.lower() == "nan" or s.lower() == "none":
            return 0
        toks = [t for t in re.split(r'[;,|\s]+', s) if t]
        return len(toks)

    def _present(cell) -> int:
        return 1 if _count_entries(cell) > 0 else 0

    tbl = gpa_df[[gene_col] + sample_cols].copy()
    C = tbl[sample_cols].map(_count_entries).astype(np.int32)   # copies per sample
    M = tbl[sample_cols].map(_present).astype(np.int8)          # presence per sample
    genes = tbl[gene_col].astype(str).tolist()

    # Optional descriptors
    gene_id_ser = gpa_df.get("Gene ID")
    name_ser    = next((gpa_df[c] for c in name_cols if c in gpa_df.columns), None)
    annot_ser   = gpa_df.get("Annotation")

    # --- drop globally core genes (present in every genome) ---
    if drop_globally_core:
        core_mask = (M.sum(axis=1).to_numpy(dtype=np.int64) == M.shape[1])
        if core_mask.any():
            keep = ~core_mask
            C = C.loc[keep];  M = M.loc[keep]
            genes = [g for g, k in zip(genes, keep) if k]
            if gene_id_ser is not None: gene_id_ser = gene_id_ser.loc[keep]
            if name_ser    is not None: name_ser    = name_ser.loc[keep]
            if annot_ser   is not None: annot_ser   = annot_ser.loc[keep]

    # Optional Fisher exact (presence-based)
    fisher_exact = None
    if add_fisher_from_presence:
        try:
            from scipy.stats import fisher_exact as _fisher
            fisher_exact = _fisher
        except Exception:
            fisher_exact = None

    def _bh_fdr(p: pd.Series) -> pd.Series:
        p = pd.Series(p, dtype=float)
        v = p.dropna()
        m = len(v)
        if m == 0:
            return p
        order = v.sort_values().index
        ranks = pd.Series(np.arange(1, m+1, dtype=float), index=order)
        q = (v.loc[order] * m / ranks).iloc[::-1].cummin().iloc[::-1].clip(upper=1.0)
        out = pd.Series(index=p.index, dtype=float); out.loc[order] = q
        return out

    # Prepare output directory if requested
    out_path: Optional[Path] = None
    if out_dir is not None:
        out_path = Path(out_dir)
        out_path.mkdir(parents=True, exist_ok=True)

    by_method: Dict[str, Dict[str, pd.DataFrame]] = {}
    all_rows: list[pd.DataFrame] = []

    for method, label_map in (clusters_by_method or {}).items():
        if not isinstance(label_map, dict) or not label_map:
            continue

        labeled_samples = [s for s in sample_cols if s in label_map]
        if len(labeled_samples) < 2:
            continue

        # clusters -> sample lists (drop tiny clusters)
        cl_to_samples: Dict[str, list[str]] = {}
        for s in labeled_samples:
            cl_to_samples.setdefault(str(label_map[s]), []).append(s)
        cl_to_samples = {cl: ss for cl, ss in cl_to_samples.items() if len(ss) >= min_cluster_size}
        if not cl_to_samples:
            continue

        # Optionally drop genes with zero copies across considered samples (per method)
        if drop_never_present:
            keep = (C[labeled_samples].sum(axis=1).to_numpy(dtype=np.int64) > 0)
            C_use = C.loc[keep, labeled_samples]
            M_use = M.loc[keep, labeled_samples]
            genes_use = [g for g, k in zip(genes, keep) if k]
            gene_id_use = gene_id_ser.loc[keep].to_numpy(object) if gene_id_ser is not None else np.array([None]*keep.sum(), dtype=object)
            gene_name_use = (name_ser.loc[keep].astype(str).replace("nan","").to_numpy(object)
                             if name_ser is not None else np.array([None]*keep.sum(), dtype=object))
            annot_use = (annot_ser.loc[keep].astype(str).replace("nan","").to_numpy(object)
                         if annot_ser is not None else np.array([None]*keep.sum(), dtype=object))
        else:
            C_use, M_use = C[labeled_samples], M[labeled_samples]
            genes_use = genes
            gene_id_use = gene_id_ser.to_numpy(object) if gene_id_ser is not None else np.array([None]*len(genes_use), dtype=object)
            gene_name_use = (name_ser.astype(str).replace("nan","").to_numpy(object)
                             if name_ser is not None else np.array([None]*len(genes_use), dtype=object))
            annot_use = (annot_ser.astype(str).replace("nan","").to_numpy(object)
                         if annot_ser is not None else np.array([None]*len(genes_use), dtype=object))

        by_cluster: Dict[str, pd.DataFrame] = {}

        # Precompute total copies across all labeled samples (per gene)
        total_copies_all = C_use.sum(axis=1).to_numpy(dtype=np.int64)

        for cl, in_samples in cl_to_samples.items():
            out_samples = [s for s in labeled_samples if s not in in_samples]

            # copies per gene
            in_copies  = C_use[in_samples].sum(axis=1).to_numpy(dtype=np.int64)
            out_copies = C_use[out_samples].sum(axis=1).to_numpy(dtype=np.int64) if out_samples else np.zeros_like(in_copies)
            total_copies = in_copies + out_copies  # equals total_copies_all

            is_unique_by_copy = (in_copies > 0) & (out_copies == 0)

            with np.errstate(divide="ignore", invalid="ignore"):
                copy_share_in = np.where(total_copies > 0, in_copies / total_copies, np.nan)
                copy_specificity = np.where(total_copies > 0, 2.0*copy_share_in - 1.0, np.nan)  # [-1,1]

            df = pd.DataFrame({
                "method": method,
                "cluster": cl,
                "gene": genes_use,
                "gene_id": gene_id_use,
                "gene_name": gene_name_use,
                "annotation": annot_use,
                "in_copies": in_copies,
                "out_copies": out_copies,
                "total_copies": total_copies,
                "is_unique_by_copy": is_unique_by_copy,
                "copy_share_in": copy_share_in,
                "copy_specificity": copy_specificity,
            })

            # Presence-based Fisher (optional)
            have_p = False
            if add_fisher_from_presence and fisher_exact is not None and len(out_samples) > 0:
                a = M_use[in_samples].sum(axis=1).to_numpy(dtype=np.int64)                 # genomes IN with gene
                b = len(in_samples) - a                                                     # genomes IN without
                c = M_use[out_samples].sum(axis=1).to_numpy(dtype=np.int64)                # genomes OUT with gene
                d = len(out_samples) - c                                                    # genomes OUT without
                pvals = np.full(df.shape[0], np.nan, dtype=float)
                for i, (ai, bi, ci, di) in enumerate(zip(a, b, c, d)):
                    try:
                        pvals[i] = fisher_exact([[int(ai), int(bi)], [int(ci), int(di)]], alternative="greater")[1]
                    except Exception:
                        pvals[i] = np.nan
                df["p_value"] = pvals
                df["q_value"] = _bh_fdr(df["p_value"])
                have_p = True

            # --- filter out unwanted rows ---
            if drop_insignificant and have_p:
                use_col = "q_value" if (use_fdr and "q_value" in df.columns) else "p_value"
                mask = df[use_col] <= float(alpha)
                if mask.notna().any():
                    df = df[mask.fillna(False)]

            if drop_unspecific and "copy_specificity" in df.columns:
                mask = df["copy_specificity"] >= float(unspecific_threshold)
                if mask.notna().any():
                    df = df[mask.fillna(False)]

            if drop_hypothetical and "annotation" in df.columns:
                ann = df["annotation"].astype(str).str.strip().str.casefold()
                df = df[~ann.eq("hypothetical protein")]

            df = df.reset_index(drop=True)

            # Sort: unique-by-copy first, then most concentrated inside
            if not df.empty:
                df = df.sort_values(
                    by=["is_unique_by_copy", "copy_specificity", "copy_share_in", "total_copies"],
                    ascending=[False, False, False, False]
                ).reset_index(drop=True)

            by_cluster[cl] = df
            all_rows.append(df)

            # --- Save per-cluster CSV if requested ---
            if out_path is not None and (save_empty or not df.empty):
                method_tag = re.sub(r"[^A-Za-z0-9]+", "", str(method).upper())
                cluster_tag = re.sub(r"[^A-Za-z0-9]+", "", str(cl))
                fname = f"{method_tag}_{cluster_tag}.csv"
                df.to_csv(out_path / fname, index=False)

        by_method[method] = by_cluster

    long_df = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame(
        columns=["method","cluster","gene","gene_id","gene_name","annotation",
                 "in_copies","out_copies","total_copies",
                 "is_unique_by_copy","copy_share_in","copy_specificity","p_value","q_value"]
    )
    return by_method, long_df
