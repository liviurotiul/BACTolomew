import numpy as np
import pandas as pd

from skbio import DistanceMatrix
from skbio.stats.distance import permanova as sk_permanova, anosim as sk_anosim
try:
    from skbio.stats.distance import permdisp as sk_permdisp
    _HAS_PERMDISP = True
except Exception:
    _HAS_PERMDISP = False

from sklearn.metrics import silhouette_score, adjusted_rand_score, normalized_mutual_info_score, v_measure_score

def _coerce_dm(D_or_dm, samples=None) -> tuple[DistanceMatrix, np.ndarray, list[str]]:
    """
    Return (DistanceMatrix, ndarray D, ids). If D_or_dm is ndarray, require samples.
    """
    if isinstance(D_or_dm, DistanceMatrix):
        dm = D_or_dm
        return dm, np.asarray(dm.data, dtype=float), list(dm.ids)
    # ndarray path
    D = np.asarray(D_or_dm, dtype=float)
    if D.ndim != 2 or D.shape[0] != D.shape[1]:
        raise ValueError("D must be a square distance matrix.")
    if samples is None:
        raise ValueError("When passing a NumPy distance matrix, you must also provide `samples` to define row/col order.")
    if len(samples) != D.shape[0]:
        raise ValueError("Length of `samples` must match D.shape[0].")
    return DistanceMatrix(D, ids=list(map(str, samples))), D, list(map(str, samples))

def _sanitize_labels(col: pd.Series, other_label: str) -> pd.Series:
    return (
        col.astype("string")
           .str.strip()
           .replace({"": pd.NA, "nan": pd.NA, "NaN": pd.NA, "NONE": pd.NA, "None": pd.NA})
           .fillna(other_label)
           .astype(str)
    )

def _labels_from_metadata(ids: list[str], metadata_df: pd.DataFrame, column: str, other_label="Other") -> np.ndarray:
    lower = {c.lower(): c for c in metadata_df.columns}
    if column.lower() == "sample" or column.lower() not in lower:
        raise ValueError(f"Metadata column '{column}' not found (or is 'Sample').")
    col = lower[column.lower()]
    md = metadata_df[["Sample", col]].copy()
    md[col] = _sanitize_labels(md[col], other_label)
    smap = dict(zip(md["Sample"].astype(str), md[col]))
    return np.array([smap.get(s, other_label) for s in ids], dtype=object)

def _pseudo_r2_from_distance_labels(D: np.ndarray, y: np.ndarray) -> float:
    """
    Pseudo R² used in PERMANOVA context: (SSt - SSw) / SSt with squared distances.
    """
    Dsq = D * D
    n = D.shape[0]
    if n < 2:
        return np.nan
    SSt = Dsq.sum() / n
    if not np.isfinite(SSt) or SSt <= 0:
        return np.nan
    SSw = 0.0
    cats = pd.Categorical(y).categories
    for lvl in cats:
        idx = np.where(y == lvl)[0]
        if len(idx):
            SSw += Dsq[np.ix_(idx, idx)].sum() / len(idx)
    return float(max(SSt - SSw, 0.0) / SSt)

def _extract_stat_p_from_mapping(obj):
    """
    Accepts a pandas Series/row or dict-like, and searches common key aliases.
    """
    # Normalize keys to a compact form for matching
    if hasattr(obj, "to_dict"):
        d = obj.to_dict()
    else:
        d = dict(obj)

    norm = {str(k).strip().lower().replace(" ", "").replace("-", "").replace("_", ""): v
            for k, v in d.items()}

    # Possible aliases for statistic and p-value across skbio versions
    stat_keys = ("teststatistic", "statistic", "f", "r", "pseudof", "pseudor")
    pval_keys = ("pvalue", "p", "pval", "pvalueperm", "pvalue(permutation)", "pvaluepermutation")

    stat = next((norm[k] for k in stat_keys if k in norm), None)
    pval = next((norm[k] for k in pval_keys if k in norm), None)

    if stat is None or pval is None:
        # Some DataFrames store columns with spaces/hyphens; try looser lookup
        # Look for any key that contains 'p' and 'value'
        if pval is None:
            pval = next((v for k, v in norm.items() if "p" in k and "value" in k), None)
        if stat is None:
            # prefer keys containing 'stat'
            stat = next((v for k, v in norm.items() if "stat" in k), None)
            # or fall back to 'f'/'r' if present
            if stat is None:
                stat = norm.get("f", norm.get("r", None))

    if stat is None or pval is None:
        raise KeyError(f"Could not locate statistic/p-value keys in {list(d.keys())}")

    return float(stat), float(pval)

def _group_counts(labels):
    vc = pd.Series(labels).value_counts()
    n_groups = int((vc > 0).sum())
    n_singletons = int((vc == 1).sum())
    n_ge2 = int((vc >= 2).sum())
    return vc, n_groups, n_singletons, n_ge2

def _srget(sr: pd.Series, *keys, default=np.nan):
    """Robustly read stat/p from a scikit-bio results Series using several aliases."""
    for k in keys:
        if k in sr:
            return sr[k]
    return default

def _sanitize_labels(col: pd.Series, other_label: str) -> pd.Series:
    return (
        col.astype("string")
           .str.strip()
           .replace({"": pd.NA, "nan": pd.NA, "NaN": pd.NA, "NONE": pd.NA, "None": pd.NA})
           .fillna(other_label)
           .astype(str)
    )

def _labels_from_metadata(ids: list[str], metadata_df: pd.DataFrame, column: str, other_label="Other") -> np.ndarray:
    lower = {c.lower(): c for c in metadata_df.columns}
    if column.lower() == "sample" or column.lower() not in lower:
        raise ValueError(f"Metadata column '{column}' not found (or is 'Sample').")
    col = lower[column.lower()]
    md = metadata_df[["Sample", col]].copy()
    md[col] = _sanitize_labels(md[col], other_label)
    smap = dict(zip(md["Sample"].astype(str), md[col]))
    return np.array([smap.get(s, other_label) for s in ids], dtype=object)

def _pseudo_r2_from_distance_labels(D: np.ndarray, y: np.ndarray) -> float:
    """
    Pseudo R² used in PERMANOVA context: (SSt - SSw) / SSt with squared distances.
    """
    Dsq = D * D
    n = D.shape[0]
    if n < 2:
        return np.nan
    SSt = Dsq.sum() / n
    if not np.isfinite(SSt) or SSt <= 0:
        return np.nan
    SSw = 0.0
    cats = pd.Categorical(y).categories
    for lvl in cats:
        idx = np.where(y == lvl)[0]
        if len(idx):
            SSw += Dsq[np.ix_(idx, idx)].sum() / len(idx)
    return float(max(SSt - SSw, 0.0) / SSt)

def _safe_dm(D_or_dm, samples=None):
    """Coerce to DistanceMatrix + ndarray; fix obvious issues and assert finite."""
    dm, D, ids = _coerce_dm(D_or_dm, samples)
    # force symmetry and zero diagonal (tolerant to tiny FP error)
    D = np.asarray(D, float)
    if not np.allclose(D, D.T, atol=1e-12, equal_nan=False):
        D = 0.5 * (D + D.T)
    np.fill_diagonal(D, 0.0)
    if not np.isfinite(D).all():
        raise ValueError("Distance matrix contains NaN/inf.")
    return DistanceMatrix(D, ids=ids), D, ids

def _extract_stat_p(res):
    """Return (statistic, pvalue) from scikit-bio results across versions/shapes."""
    # 1) Newer scikit-bio: attributes on the result object
    stat = getattr(res, "statistic", None)
    pval = getattr(res, "p_value", None)
    if stat is not None and pval is not None:
        return float(stat), float(pval)

    # 2) Result exposes .to_series()
    if hasattr(res, "to_series"):
        sr = res.to_series()
        # Fall through to the generic Series/DataFrame handler below
        try:
            return _extract_stat_p_from_mapping(sr)
        except Exception:
            pass

    # 3) It might already be a pandas object or a namedtuple/dict-like
    try:
        import pandas as pd
        if isinstance(res, pd.Series):
            return _extract_stat_p_from_mapping(res)
        if isinstance(res, pd.DataFrame):
            # typical scikit-bio older output: 1-row DataFrame
            if len(res.index) >= 1:
                return _extract_stat_p_from_mapping(res.iloc[0])
        if hasattr(res, "_asdict"):  # namedtuple
            return _extract_stat_p_from_mapping(res._asdict())
        if isinstance(res, dict):
            return _extract_stat_p_from_mapping(res)
    except Exception:
        pass

    # 4) Give up
    return float("nan"), float("nan")

def _extract_stat_p_from_mapping(obj):
    """
    Accepts a pandas Series/row or dict-like, and searches common key aliases.
    """
    # Normalize keys to a compact form for matching
    if hasattr(obj, "to_dict"):
        d = obj.to_dict()
    else:
        d = dict(obj)

    norm = {str(k).strip().lower().replace(" ", "").replace("-", "").replace("_", ""): v
            for k, v in d.items()}

    # Possible aliases for statistic and p-value across skbio versions
    stat_keys = ("teststatistic", "statistic", "f", "r", "pseudof", "pseudor")
    pval_keys = ("pvalue", "p", "pval", "pvalueperm", "pvalue(permutation)", "pvaluepermutation")

    stat = next((norm[k] for k in stat_keys if k in norm), None)
    pval = next((norm[k] for k in pval_keys if k in norm), None)

    if stat is None or pval is None:
        # Some DataFrames store columns with spaces/hyphens; try looser lookup
        # Look for any key that contains 'p' and 'value'
        if pval is None:
            pval = next((v for k, v in norm.items() if "p" in k and "value" in k), None)
        if stat is None:
            # prefer keys containing 'stat'
            stat = next((v for k, v in norm.items() if "stat" in k), None)
            # or fall back to 'f'/'r' if present
            if stat is None:
                stat = norm.get("f", norm.get("r", None))

    if stat is None or pval is None:
        raise KeyError(f"Could not locate statistic/p-value keys in {list(d.keys())}")

    return float(stat), float(pval)

def analyze_metadata_traits(
    D_or_dm,
    metadata_df: pd.DataFrame,
    *,
    samples: list[str] | None = None,
    permutations: int = 999,
    other_label: str = "Other",
    min_per_group: int = 3,
    drop_other: bool = True,
    compute_silhouette: bool = True,
    max_silhouette_n: int = 200,
    n_jobs: int = 1,   # parallel over traits
) -> pd.DataFrame:
    # ---- coerce DM once
    if isinstance(D_or_dm, DistanceMatrix):
        dm_full = D_or_dm
        D_full = np.asarray(dm_full.data, float)
        ids_full = list(map(str, dm_full.ids))
    else:
        D_full = np.asarray(D_or_dm, float)
        if D_full.ndim != 2 or D_full.shape[0] != D_full.shape[1]:
            raise ValueError("D must be square.")
        if samples is None or len(samples) != D_full.shape[0]:
            raise ValueError("Provide samples matching D.")
        ids_full = list(map(str, samples))
        dm_full = DistanceMatrix(D_full, ids=ids_full)

    # symmetrize/clean once
    if not np.allclose(D_full, D_full.T, atol=1e-12):
        D_full = 0.5 * (D_full + D_full.T)
        dm_full = DistanceMatrix(D_full, ids=ids_full)
    np.fill_diagonal(D_full, 0.0)
    if not np.isfinite(D_full).all():
        raise ValueError("Distance matrix contains NaN/inf.")

    # ---- label matrix once, aligned to ids_full
    md = metadata_df.copy()
    colmap = {c.lower(): c for c in md.columns}
    if "sample" not in colmap:
        raise ValueError("metadata_df must have 'Sample' column")
    md = md.rename(columns={colmap["sample"]: "Sample"})
    md["Sample"] = md["Sample"].astype(str)

    cols = [c for c in md.columns if c != "Sample"]
    if not cols:
        return pd.DataFrame(columns=[
            "trait","n_samples","n_groups","n_singletons","group_sizes",
            "permanova_F","permanova_p","R2","permdisp_F","permdisp_p",
            "anosim_R","anosim_p","silhouette","why_perm","why_disp"
        ])

    md_idx = (
        md.set_index("Sample")[cols]
          .apply(lambda s: s.astype("string").str.strip()
                 .replace({"": pd.NA, "nan": pd.NA, "NaN": pd.NA, "NONE": pd.NA, "None": pd.NA}))
    )
    md_idx = md_idx.reindex(ids_full)
    md_idx = md_idx.fillna(other_label).astype(str)
    Y = md_idx.values  # shape (n_samples, n_traits)
    traits = list(md_idx.columns)

    def _pseudo_r2(D, y):
        D2 = D * D
        n = D.shape[0]
        if n < 2: return np.nan
        SSt = D2.sum() / n
        if not np.isfinite(SSt) or SSt <= 0: return np.nan
        SSw = 0.0
        for lvl in pd.unique(y):
            idx = np.where(y == lvl)[0]
            if idx.size:
                SSw += D2[np.ix_(idx, idx)].sum() / idx.size
        return float(max(SSt - SSw, 0.0) / SSt)

    def _one(trait_idx: int):
        trait = traits[trait_idx]
        y_full = Y[:, trait_idx]

        keep = np.ones(len(y_full), dtype=bool)
        if drop_other:
            keep &= (y_full != other_label)

        # min group size filter
        if keep.any():
            vc = pd.Series(y_full[keep]).value_counts()
            big = set(vc.index[vc.values >= min_per_group])
            keep &= np.isin(y_full, list(big))

        idx = np.where(keep)[0]
        y = y_full[idx]
        if idx.size == 0:
            return {
                "trait": trait, "n_samples": 0, "n_groups": 0, "n_singletons": 0, "group_sizes": "",
                "permanova_F": np.nan, "permanova_p": np.nan, "R2": np.nan,
                "permdisp_F": np.nan, "permdisp_p": np.nan, "anosim_R": np.nan, "anosim_p": np.nan,
                "silhouette": np.nan, "why_perm": "no data after filtering", "why_disp": "no data after filtering"
            }

        D = D_full[np.ix_(idx, idx)]
        ids = [ids_full[i] for i in idx]

        vc = pd.Series(y).value_counts()
        n_groups = int((vc > 0).sum())
        n_singletons = int((vc == 1).sum())
        n_ge2 = int((vc >= 2).sum())
        group_sizes_str = "; ".join([f"{k}:{int(v)}" for k, v in vc.items()])

        permanova_F = permanova_p = permdisp_F = permdisp_p = anosim_R = anosim_p = sil = np.nan
        why_perm = why_disp = None
        R2 = np.nan

        if len(y) >= 3 and n_groups >= 2:
            try:
                res = sk_permanova(DistanceMatrix(D, ids=ids), y, permutations=permutations)
                if hasattr(res, "statistic"):
                    permanova_F = float(res.statistic); permanova_p = float(res.p_value)
                else:
                    sr = res.to_series() if hasattr(res, "to_series") else None
                    if sr is not None:
                        permanova_F = float(sr.get("test statistic", sr.get("pseudo-F", np.nan)))
                        permanova_p = float(sr.get("p-value", np.nan))
            except Exception as e:
                why_perm = f"permanova error: {e}"

            try:
                R2 = _pseudo_r2(D, y)
            except Exception as e:
                why_perm = f"{why_perm}; R2 error: {e}" if why_perm else f"R2 error: {e}"

            try:
                res = sk_anosim(DistanceMatrix(D, ids=ids), y, permutations=permutations)
                if hasattr(res, "statistic"):
                    anosim_R = float(res.statistic); anosim_p = float(res.p_value)
                else:
                    sr = res.to_series() if hasattr(res, "to_series") else None
                    if sr is not None:
                        anosim_R = float(sr.get("test statistic", sr.get("R", np.nan)))
                        anosim_p = float(sr.get("p-value", np.nan))
            except Exception as e:
                pass

            if compute_silhouette and n_ge2 >= 2 and len(y) <= max_silhouette_n:
                try:
                    sil = float(silhouette_score(D, y, metric="precomputed"))
                except Exception:
                    pass
        else:
            why_perm = (f"after filtering: kept {len(y)} samples across {n_groups} groups; "
                        f"need ≥3 samples and ≥2 groups (min_per_group={min_per_group})")

        if _HAS_PERMDISP:
            if n_ge2 >= 2:
                try:
                    res = sk_permdisp(DistanceMatrix(D, ids=ids), y, permutations=permutations)
                    if hasattr(res, "statistic"):
                        permdisp_F = float(res.statistic); permdisp_p = float(res.p_value)
                    else:
                        sr = res.to_series() if hasattr(res, "to_series") else None
                        if sr is not None:
                            permdisp_F = float(sr.get("test statistic", np.nan))
                            permdisp_p = float(sr.get("p-value", np.nan))
                except Exception as e:
                    why_disp = f"permdisp error: {e}"
            else:
                why_disp = (f"after filtering: need ≥2 groups with ≥2 samples each; "
                            f"groups(n≥2)={n_ge2}")

        return {
            "trait": trait,
            "n_samples": int(len(y)),
            "n_groups": int(n_groups),
            "n_singletons": int(n_singletons),
            "group_sizes": group_sizes_str,
            "permanova_F": permanova_F,
            "permanova_p": permanova_p,
            "R2": R2,
            "permdisp_F": permdisp_F,
            "permdisp_p": permdisp_p,
            "anosim_R": anosim_R,
            "anosim_p": anosim_p,
            "silhouette": sil,
            "why_perm": why_perm,
            "why_disp": why_disp,
        }

    if n_jobs == 1:
        rows = [_one(i) for i in range(len(traits))]
    else:
        from joblib import Parallel, delayed
        rows = Parallel(n_jobs=n_jobs, prefer="processes")(delayed(_one)(i) for i in range(len(traits)))

    out = pd.DataFrame(rows).sort_values(
        by=["R2", "permanova_p"], ascending=[False, True], na_position="last"
    ).reset_index(drop=True)
    return out


def run_permanova(D_or_dm, metadata_df, column: str, *, samples=None, permutations=999):
    dm, D, ids = _coerce_dm(D_or_dm, samples)
    y = _labels_from_metadata(ids, metadata_df, column)
    res = sk_permanova(dm, y, permutations=permutations)  # <-- positional
    out = res.to_series().to_dict() if hasattr(res, "to_series") else {
        "test statistic": getattr(res, "statistic", np.nan),
        "p-value": getattr(res, "p_value", np.nan),
        "permutations": permutations,
        "method name": "PERMANOVA",
    }
    out["column"] = column
    try:
        out["R2"] = _pseudo_r2_from_distance_labels(D, y)
    except Exception as e:
        print("R² error:", e)
    return out

def run_permdisp(D_or_dm, metadata_df, column: str, *, samples=None, permutations=999):
    if not _HAS_PERMDISP:
        raise RuntimeError("scikit-bio permdisp not available; upgrade scikit-bio >=0.6.")
    dm, _, ids = _coerce_dm(D_or_dm, samples)
    y = _labels_from_metadata(ids, metadata_df, column)
    res = sk_permdisp(dm, y, permutations=permutations)  # <-- positional
    return res.to_series().to_dict() if hasattr(res, "to_series") else {
        "test statistic": getattr(res, "statistic", np.nan),
        "p-value": getattr(res, "p_value", np.nan),
        "permutations": permutations,
        "method name": "PERMDISP",
    }

def run_anosim(D_or_dm, metadata_df, column: str, *, samples=None, permutations=999):
    dm, _, ids = _coerce_dm(D_or_dm, samples)
    y = _labels_from_metadata(ids, metadata_df, column)
    res = sk_anosim(dm, y, permutations=permutations)  # <-- positional
    return res.to_series().to_dict() if hasattr(res, "to_series") else {
        "test statistic": getattr(res, "statistic", np.nan),
        "p-value": getattr(res, "p_value", np.nan),
        "permutations": permutations,
        "method name": "ANOSIM",
    }

# Optional: label agreement helpers (unchanged)
def silhouette_from_labels(D: np.ndarray, labels: np.ndarray) -> float:
    return float(silhouette_score(D, labels, metric="precomputed"))

def label_agreements(true_labels: np.ndarray, pred_labels: np.ndarray) -> dict:
    return {
        "ARI": adjusted_rand_score(true_labels, pred_labels),
        "NMI": normalized_mutual_info_score(true_labels, pred_labels),
        "V_measure": v_measure_score(true_labels, pred_labels),
    }
