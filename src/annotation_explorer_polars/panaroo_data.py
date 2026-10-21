from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from .panaroo_common import _file_signature
from .polars_io import read_csv_df

_META_COLUMNS = {
    "Gene",
    "Non-unique Gene name",
    "Annotation",
    "No. isolates",
    "No. sequences",
    "Avg sequences per isolate",
    "Min group size nuc",
    "Max group size nuc",
    "Avg group size nuc",
    "Genome Fragment",
    "Order within Fragment",
    "Accessory Fragment",
    "Accessory Order with Fragment",
    "QC",
}

def _find_panaroo_out(output_dir: Path) -> Optional[Path]:
    direct = output_dir / "panaroo" / "out"
    if (direct / "gene_presence_absence.csv").is_file():
        return direct
    root = output_dir / "panaroo"
    if not root.exists():
        return None
    candidates = [p.parent for p in root.rglob("gene_presence_absence.csv") if p.is_file()]
    if not candidates:
        return None
    candidates.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return candidates[0]


@lru_cache(maxsize=16)
def _load_presence_matrix_cached(
    path_str: str,
    _mtime_ns: int,
    _size: int,
) -> Tuple[np.ndarray, Tuple[str, ...], Tuple[str, ...], Tuple[str, ...]]:
    path = Path(path_str)
    df = read_csv_df(path, low_memory=False)
    sample_cols = [c for c in df.columns if c not in _META_COLUMNS]
    if not sample_cols:
        return np.zeros((0, 0), dtype=np.uint8), tuple(), tuple(), tuple()
    raw = df[sample_cols].fillna("").astype(str)
    presence = (raw.ne("")) & (raw.ne("0"))
    genes = df["Gene"].astype(str).tolist() if "Gene" in df.columns else [f"gene_{i}" for i in range(len(df))]
    products = df["Annotation"].fillna("").astype(str).tolist() if "Annotation" in df.columns else [""] * len(df)
    presence_arr = presence.to_numpy(dtype=np.uint8)
    presence_arr.setflags(write=False)
    return presence_arr, tuple(str(s) for s in sample_cols), tuple(genes), tuple(products)


def _load_presence_matrix(panaroo_out: Path) -> Tuple[np.ndarray, List[str], List[str], List[str]]:
    path = panaroo_out / "gene_presence_absence.csv"
    presence, sample_cols, genes, products = _load_presence_matrix_cached(*_file_signature(path))
    return presence, list(sample_cols), list(genes), list(products)


def _find_fastgear_sample_clades(output_dir: Path) -> Optional[Path]:
    direct = output_dir / "fastgear" / "sample_clades.tsv"
    if direct.is_file():
        return direct
    fastgear_root = output_dir / "fastgear"
    if not fastgear_root.exists():
        return None
    matches = [p for p in fastgear_root.rglob("sample_clades.tsv") if p.is_file()]
    if not matches:
        return None
    matches.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return matches[0]


@lru_cache(maxsize=16)
def _load_fastgear_assignments_cached(
    path_str: str,
    _mtime_ns: int,
    _size: int,
) -> Tuple[Tuple[str, str, str], ...]:
    path = Path(path_str)
    df = read_csv_df(path, sep="\t", dtype=str)
    required = {"Sample", "FastGearLineage", "FastGearCluster"}
    if df.empty or not required.issubset(df.columns):
        return tuple()
    subset = (
        df.loc[:, ["Sample", "FastGearLineage", "FastGearCluster"]]
        .fillna("")
        .astype(str)
        .assign(
            Sample=lambda d: d["Sample"].str.strip(),
            FastGearLineage=lambda d: d["FastGearLineage"].str.strip(),
            FastGearCluster=lambda d: d["FastGearCluster"].str.strip(),
        )
    )
    subset = subset[(subset["Sample"] != "") & (subset["FastGearLineage"] != "")]
    if subset.empty:
        return tuple()
    subset = subset.drop_duplicates(subset=["Sample"], keep="first")
    return tuple(
        (row["Sample"], row["FastGearLineage"], row["FastGearCluster"])
        for _, row in subset.iterrows()
    )


def _load_fastgear_lineages(sample_clades_path: Path) -> Dict[str, str]:
    pairs = _load_fastgear_assignments_cached(*_file_signature(sample_clades_path))
    return {sample: lineage for sample, lineage, _cluster in pairs}


def _load_fastgear_clusters(sample_clades_path: Path) -> Dict[str, str]:
    pairs = _load_fastgear_assignments_cached(*_file_signature(sample_clades_path))
    return {sample: cluster for sample, _lineage, cluster in pairs if str(cluster).strip()}


def _normalize_resistance_token(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    if "," in text:
        text = text.split(",", 1)[0].strip()
    return re.sub(r"\s+", " ", text)


def _extract_resistance_gene_labels(row: Mapping[str, object]) -> List[str]:
    preferred_gene = _normalize_resistance_token(row.get("Gene"))
    if preferred_gene:
        return [preferred_gene]

    out: List[str] = []
    seen: set[str] = set()

    def _add(value: object) -> None:
        token = _normalize_resistance_token(value)
        if not token:
            return
        key = token.casefold()
        if key in seen:
            return
        seen.add(key)
        out.append(token)

    for field in ("ResistanceAlleles", "ResistanceGeneFamilies"):
        raw = str(row.get(field) or "").strip()
        if not raw:
            continue
        for part in raw.split(";"):
            _add(part)
        if out:
            return out

    fallback = _normalize_resistance_token(row.get("Product"))
    return [fallback] if fallback else []


def _normalize_mge_token(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    text = re.sub(r"\s+", " ", text)
    return text.strip(" |;,:")


def _extract_mge_labels(row: Mapping[str, object]) -> List[str]:
    db_name = str(row.get("DB") or "").strip()
    product = str(row.get("Product") or "").strip()
    if not db_name or not product:
        return []

    candidates: List[str] = []
    db_name_upper = db_name.upper()
    if "ISFINDER" in db_name_upper:
        token = product.split("//")[-1].strip()
        token = token.split("~~~", 1)[0].strip()
        token = token.split("|", 1)[0].strip()
        if token:
            candidates.append(token)
    if "CGE-MGE" in db_name_upper:
        token = product.split("|", 1)[0].strip()
        if token:
            candidates.append(token)

    if not candidates:
        match = re.search(r"\b(?:IS[A-Za-z0-9._-]+|Tn[A-Za-z0-9._-]+)\b", product)
        if match:
            candidates.append(match.group(0))

    out: List[str] = []
    seen: set[str] = set()
    for candidate in candidates:
        token = _normalize_mge_token(candidate)
        if not token:
            continue
        if not (token.startswith("IS") or token.startswith("Tn")):
            continue
        key = token.casefold()
        if key in seen:
            continue
        seen.add(key)
        out.append(token)
    return out


def _extract_virulence_labels(row: Mapping[str, object]) -> List[str]:
    generic_labels = {
        "hypothetical protein",
        "putative protein",
        "putative",
        "protein of unknown function",
        "uncharacterized protein",
    }

    def _normalize(value: object) -> str:
        text = str(value or "").strip()
        if not text:
            return ""
        text = re.sub(r"\s+", " ", text)
        return text.strip(" |;,:")

    def _valid(label: str) -> bool:
        if not label:
            return False
        return label.casefold() not in generic_labels

    gene = _normalize(row.get("Gene"))
    if _valid(gene):
        return [gene]

    gene_id = _normalize(row.get("GeneID"))
    if _valid(gene_id):
        return [gene_id]

    product = _normalize(row.get("Product"))
    if _valid(product):
        return [product]

    return []


def _select_samples_from_presence(
    presence: np.ndarray,
    sample_labels: Sequence[str],
    sample_filter: Optional[set[str]],
) -> Tuple[np.ndarray, List[str]]:
    if sample_filter is None:
        return presence, list(sample_labels)
    keep_idx = [i for i, name in enumerate(sample_labels) if str(name).strip() in sample_filter]
    if not keep_idx:
        return np.zeros((presence.shape[0], 0), dtype=presence.dtype), []
    sel = presence[:, keep_idx]
    labels = [str(sample_labels[i]) for i in keep_idx]
    return sel, labels
