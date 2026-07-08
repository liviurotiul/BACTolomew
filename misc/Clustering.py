import pandas as pd
import numpy as np

def pcoa_from_gene_absence_jaccard(
    df: "pd.DataFrame",
    *,
    gene_col: str | None = None,
    invert: bool = False,        # Jaccard usually on presence; set True to use absence
    n_components: int | None = 2,
    standardize: bool = False,   # unused; kept for API symmetry
) -> dict:
    warnings: list[str] = []
    ok = True

    if not isinstance(df, pd.DataFrame) or df.empty:
        return {"ok": False, "warnings": ["Input is empty or not a DataFrame."]}

    cols = df.columns.astype(str).tolist()
    if gene_col is None:
        gene_col = "Gene" if "Gene" in cols else cols[0]
    if gene_col not in cols:
        return {"ok": False, "warnings": [f"Gene column '{gene_col}' not found."]}

    meta_like = {
        gene_col, "Gene ID", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences",
        "Avg sequences per isolate", "Genome Fragment", "Order within Fragment", "Accessory Fragment",
        "Accessory Order with Fragment", "QC", "Min group size nuc", "Max group size nuc", "Avg group size nuc"
    }
    sample_cols = [c for c in cols if c not in meta_like]
    if len(sample_cols) < 2:
        return {"ok": False, "warnings": ["Fewer than 2 sample columns detected."]}

    g = df[[gene_col] + sample_cols].copy()
    def _present(x):
        if x is None: return 0
        s = str(x).strip().lower()
        return 0 if (s == "" or s == "0" or s == "nan" or s == "none") else 1

    P = g[sample_cols].map(_present).astype(np.int8)
    genes_all = g[gene_col].astype(str).tolist()
    A = (1 - P) if invert else P

    is_all_missing = g[sample_cols].isna().all(axis=1) | (
        g[sample_cols].astype(str).apply(lambda s: s.str.strip() == "").all(axis=1)
    )
    removed_all_missing_idx = np.where(is_all_missing.values)[0].tolist()
    removed_all_missing = [genes_all[i] for i in removed_all_missing_idx]

    keep_mask_initial = ~is_all_missing.values
    A0 = A.loc[keep_mask_initial].to_numpy(dtype=np.int8)
    genes0 = [genes_all[i] for i, k in enumerate(keep_mask_initial) if k]
    if A0.shape[0] == 0:
        return {"ok": False, "warnings": ["All genes were missing; nothing to analyze."]}

    var_per_gene = A0.var(axis=1)
    keep_mask_var = var_per_gene > 0.0
    removed_zero_var = [g for g, k in zip(genes0, keep_mask_var) if not k]

    A1 = A0[keep_mask_var, :]
    genes_kept = [g for g, k in zip(genes0, keep_mask_var) if k]
    if A1.shape[0] < 1:
        return {
            "ok": False,
            "warnings": ["No variable genes remain after filtering."],
            "removed_genes_all_missing": removed_all_missing,
            "removed_genes_zero_variance": removed_zero_var,
        }

    X = A1.T.astype(np.float64)     # samples x genes (binary)
    samples = sample_cols
    n = X.shape[0]
    if n < 2:
        return {"ok": False, "warnings": ["Need at least 2 samples."]}

    # Jaccard distance matrix (vectorized)
    ones_per_sample = X.sum(axis=1)                      # (n,)
    intersection = X @ X.T                               # (n,n)
    union = ones_per_sample[:, None] + ones_per_sample[None, :] - intersection
    with np.errstate(divide="ignore", invalid="ignore"):
        J = np.where(union > 0, intersection / union, 1.0)  # similarity
    D = 1.0 - J                                            # distance
    np.fill_diagonal(D, 0.0)

    # PCoA (classical MDS)
    D2 = D * D
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * (H @ D2 @ H)
    B = (B + B.T) * 0.5

    eigvals, eigvecs = np.linalg.eigh(B)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    pos = eigvals > 1e-12
    if not np.any(pos):
        warnings.append("All eigenvalues non-positive; coordinates are zero.")
        ok = False
        k = 0 if n_components is None else max(0, int(n_components))
        scores = np.zeros((n, k))
        explained = np.zeros((k,))
        ratio = np.zeros((k,))
        return {
            "ok": ok,
            "warnings": warnings,
            "scores": scores,
            "explained_variance": explained,
            "explained_variance_ratio": ratio,
            "samples": samples,
            "genes": genes_kept,
            "removed_genes_zero_variance": removed_zero_var,
            "removed_genes_all_missing": removed_all_missing,
            "distance_matrix": D,
            "params": {
                "gene_col": gene_col,
                "sample_cols": sample_cols,
                "invert": invert,
                "n_components": k,
                "metric": "jaccard",
                "method": "pcoa",
            },
        }

    eigvals_pos = eigvals[pos]
    eigvecs_pos = eigvecs[:, pos]
    r = eigvals_pos.size
    k = r if n_components is None else max(1, min(int(n_components), r))

    coords = eigvecs_pos[:, :k] * np.sqrt(eigvals_pos[:k])
    explained = eigvals_pos[:k]
    ratio = explained / eigvals_pos.sum()

    return {
        "ok": ok,
        "warnings": warnings,
        "scores": coords,                             # n_samples x k
        "explained_variance": explained,
        "explained_variance_ratio": ratio,
        "samples": samples,
        "genes": genes_kept,
        "removed_genes_zero_variance": removed_zero_var,
        "removed_genes_all_missing": removed_all_missing,
        "distance_matrix": D,
        "params": {
            "gene_col": gene_col,
            "sample_cols": sample_cols,
            "invert": invert,
            "n_components": k,
            "metric": "jaccard",
            "method": "pcoa",
        },
    }

def pca_from_gene_absence(
    df: "pd.DataFrame",
    *,
    gene_col: str | None = None,
    invert: bool = True,
    n_components: int | None = 2,
    standardize: bool = False,
) -> dict:
    """
    Perform PCA on a gene (presence/absence) style dataframe with robust prechecks.

    Parameters
    ----------
    df : pd.DataFrame
        Table like Roary/Panaroo 'gene_presence_absence.csv'-style, where one column
        holds the gene/group name and the per-sample columns contain either blanks/NaN/0
        for absence and non-empty strings/numbers for presence.
        If `invert=True` (default), 1 means ABSENCE and 0 means presence in the matrix used for PCA.
    gene_col : str | None
        Name of the gene column. If None, the function will auto-detect it (prefers 'Gene').
    invert : bool
        If True, invert presence to model ABSENCE (1=absent, 0=present). If False, 1=present.
    n_components : int | None
        Number of principal components to return. If None, uses min(n_samples, n_features).
    standardize : bool
        If True, z-score each gene feature after centering (useful when features have unequal variance).

    Returns
    -------
    dict
        {
          'X': np.ndarray [n_samples, n_features]   # centered (and optionally standardized) matrix
          'scores': np.ndarray [n_samples, k],      # sample coordinates in PC space
          'components': np.ndarray [k, n_features], # PC loadings (like sklearn.components_)
          'feature_loadings': np.ndarray [n_features, k], # loadings scaled by sqrt(expl.var)
          'explained_variance': np.ndarray [k],
          'explained_variance_ratio': np.ndarray [k],
          'singular_values': np.ndarray [k],
          'feature_means': np.ndarray [n_features],  # means used for centering
          'feature_stds': np.ndarray [n_features],   # stds used (1 if standardize=False)
          'samples': list[str],
          'genes': list[str],                        # genes kept after filtering
          'removed_genes_zero_variance': list[str],
          'removed_genes_all_missing': list[str],
          'params': dict(...),
          'warnings': list[str],
          'ok': bool
        }
    """

    warnings: list[str] = []
    ok = True

    if not isinstance(df, pd.DataFrame) or df.empty:
        return {
            "ok": False,
            "warnings": ["Input is empty or not a DataFrame."],
        }

    # --- Detect columns: gene column + sample columns ---
    cols = df.columns.astype(str).tolist()
    if gene_col is None:
        gene_col = "Gene" if "Gene" in cols else cols[0]

    if gene_col not in cols:
        return {
            "ok": False,
            "warnings": [f"Gene column '{gene_col}' not found in DataFrame."],
        }

    meta_like = {
        gene_col, "Gene ID", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences",
        "Avg sequences per isolate", "Genome Fragment", "Order within Fragment", "Accessory Fragment",
        "Accessory Order with Fragment", "QC", "Min group size nuc", "Max group size nuc", "Avg group size nuc"
    }
    sample_cols = [c for c in cols if c not in meta_like]
    if len(sample_cols) < 2:
        return {
            "ok": False,
            "warnings": ["Fewer than 2 sample columns detected; PCA requires ≥2 samples."],
        }

    # --- Build binary presence matrix from typical gene_presence_absence tables ---
    # Map each cell: non-empty/nonzero/non-NaN => 1 (present), else 0 (absent)
    g = df[[gene_col] + sample_cols].copy()
    def _present(x):
        if x is None:
            return 0
        s = str(x).strip().lower()
        if s == "" or s == "0" or s == "nan" or s == "none":
            return 0
        return 1

    P = g[sample_cols].map(_present).astype(np.int8)  # presence matrix (genes x samples)
    genes_all = g[gene_col].astype(str).tolist()

    # Invert if we want absence=1 (common for "gene absence" analysis)
    A = (1 - P) if invert else P  # A: genes x samples (binary)

    # --- Remove genes that are all missing (all NaN before mapping) OR all-constant (zero variance) ---
    # "All missing" after mapping would be indistinguishable from "all absent" (=1) if invert=True,
    # so we check missing separately from the raw table:
    is_all_missing = g[sample_cols].isna().all(axis=1) | (g[sample_cols].astype(str).apply(lambda s: s.str.strip()=="").all(axis=1))
    removed_all_missing_idx = np.where(is_all_missing.values)[0].tolist()
    removed_all_missing = [genes_all[i] for i in removed_all_missing_idx]

    # Filter them out first
    keep_mask_initial = ~is_all_missing.values
    A0 = A.loc[keep_mask_initial].to_numpy(dtype=float)  # genes_kept0 x samples
    genes0 = [genes_all[i] for i, k in enumerate(keep_mask_initial) if k]

    if A0.shape[0] == 0:
        return {
            "ok": False,
            "warnings": ["All genes were missing; nothing to analyze."],
        }

    # Zero variance across samples (constant 0 or 1 after (inversion) mapping)
    var_per_gene = A0.var(axis=1)
    keep_mask_var = var_per_gene > 0.0
    removed_zero_var = [g for g, k in zip(genes0, keep_mask_var) if not k]

    A1 = A0[keep_mask_var, :]
    genes_kept = [g for g, k in zip(genes0, keep_mask_var) if k]

    if A1.shape[0] < 2:
        return {
            "ok": False,
            "warnings": ["Fewer than 2 variable genes remain after filtering; PCA is not meaningful."],
            "removed_genes_all_missing": removed_all_missing,
            "removed_genes_zero_variance": removed_zero_var,
        }

    # --- Arrange matrix as samples x features (genes) for PCA ---
    # Current A1 is genes x samples; transpose to samples x genes
    X = A1.T  # [n_samples, n_features]
    samples = sample_cols
    n_samples, n_features = X.shape

    if n_samples < 2:
        return {
            "ok": False,
            "warnings": ["PCA requires at least 2 samples."],
            "removed_genes_all_missing": removed_all_missing,
            "removed_genes_zero_variance": removed_zero_var,
        }

    # --- Center (and optionally standardize) features ---
    feature_means = X.mean(axis=0)
    Xc = X - feature_means

    if standardize:
        feature_stds = Xc.std(axis=0, ddof=1)
        # avoid divide-by-zero (shouldn't happen due to variance filtering, but be safe)
        safe_stds = np.where(feature_stds == 0, 1.0, feature_stds)
        Xcs = Xc / safe_stds
    else:
        feature_stds = np.ones(n_features, dtype=float)
        Xcs = Xc

    # --- PCA via SVD ---
    # Xcs = U * S * Vt  (U: n_samples x r, S: r, Vt: r x n_features)
    U, S, Vt = np.linalg.svd(Xcs, full_matrices=False)
    r = min(n_samples, n_features)
    if n_components is None:
        k = r
    else:
        if n_components < 1:
            warnings.append("Requested n_components < 1; adjusted to 1.")
            n_components = 1
        k = min(n_components, r)

    # Explained variance per component: S^2 / (n_samples - 1)
    denom = max(n_samples - 1, 1)
    explained_variance_all = (S ** 2) / denom
    total_var = explained_variance_all.sum() if explained_variance_all.size else 0.0
    if total_var <= 0:
        warnings.append("Total variance is zero after centering/standardization.")
        ok = False

    explained_variance = explained_variance_all[:k]
    explained_variance_ratio = (explained_variance / total_var) if total_var > 0 else np.zeros_like(explained_variance)

    # Components (like sklearn): k x n_features
    components = Vt[:k, :]

    # Scores (sample coordinates): n_samples x k
    scores = Xcs @ components.T

    # Feature loadings scaled by sqrt(explained variance) (one common convention)
    # shape: n_features x k
    feature_loadings = components.T * np.sqrt(explained_variance[np.newaxis, :])

    # Package results
    result = {
        "ok": ok,
        "warnings": warnings,
        "X": Xcs,  # centered (and possibly standardized) samples x features
        "scores": scores,
        "components": components,
        "feature_loadings": feature_loadings,
        "explained_variance": explained_variance,
        "explained_variance_ratio": explained_variance_ratio,
        "singular_values": S[:k],
        "feature_means": feature_means,
        "feature_stds": feature_stds,
        "samples": samples,
        "genes": genes_kept,
        "removed_genes_zero_variance": removed_zero_var,
        "removed_genes_all_missing": removed_all_missing,
        "params": {
            "gene_col": gene_col,
            "sample_cols": sample_cols,
            "invert": invert,
            "n_components": k,
            "standardize": standardize,
        },
    }
    return result

def kmeans(X: np.ndarray, k: int, *, n_init: int = 5, max_iter: int = 100, rng: int | None = None):
    rng_obj = np.random.default_rng(rng)
    best_inertia = np.inf
    best_labels = None
    best_centers = None

    for _ in range(n_init):
        # k unique random initial centers
        init_idx = rng_obj.choice(X.shape[0], size=k, replace=False)
        centers = X[init_idx].copy()

        for _ in range(max_iter):
            # assign
            d2 = ((X[:, None, :] - centers[None, :, :])**2).sum(axis=2)  # (n,k)
            labels = d2.argmin(axis=1)

            # recompute centers
            new_centers = np.vstack([
                X[labels == j].mean(axis=0) if np.any(labels == j) else centers[j]
                for j in range(k)
            ])

            if np.allclose(new_centers, centers):
                centers = new_centers
                break
            centers = new_centers

        # inertia
        d2 = ((X - centers[labels])**2).sum(axis=1).sum()
        if d2 < best_inertia:
            best_inertia = d2
            best_labels = labels
            best_centers = centers

    return best_labels, best_centers, float(best_inertia)

def calinski_harabasz(X: np.ndarray, labels: np.ndarray, centers: np.ndarray) -> float:
    """Compute Calinski–Harabasz index."""
    n, p = X.shape
    k = centers.shape[0]
    if k <= 1 or k >= n:
        return 0.0
    overall_mean = X.mean(axis=0)
    # Between-cluster dispersion
    Bk = 0.0
    # Within-cluster dispersion
    Wk = 0.0
    for j in range(k):
        Xj = X[labels == j]
        nj = Xj.shape[0]
        if nj == 0:
            continue
        diff_c = centers[j] - overall_mean
        Bk += nj * float(diff_c @ diff_c)
        # within
        diffs = Xj - centers[j]
        Wk += float((diffs * diffs).sum())
    # CH index
    return (Bk / (k - 1)) / (Wk / (n - k)) if (Wk > 0 and n > k) else 0.0

def best_kmeans(X: np.ndarray):
    """Choose k in [2..min(6, n-1)] maximizing CH index; return labels, k, centers."""
    n = X.shape[0]
    k_max = max(2, min(6, n - 1))
    best = {"score": -np.inf, "labels": None, "k": 2, "centers": None}
    for k in range(2, k_max + 1):
        labels, centers, _ = kmeans(X, k, n_init=8, max_iter=100, rng=42)
        score = calinski_harabasz(X, labels, centers)
        if score > best["score"]:
            best.update({"score": score, "labels": labels, "k": k, "centers": centers})
    return best["labels"], best["k"], best["centers"]