from __future__ import annotations

import gzip
import re
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from .panaroo_common import _file_signature
from .polars_io import read_csv_df

try:
    from scipy.spatial.distance import pdist, squareform
except Exception:
    pdist = None
    squareform = None

try:
    from skbio import DistanceMatrix
    from skbio.tree import nj
except Exception:
    DistanceMatrix = None
    nj = None

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

_PANAROO_ALIGNMENT_CANDIDATES = (
    "core_gene_alignment.aln",
    "core_gene_alignment.fa",
    "core_gene_alignment.fasta",
    "core_gene_alignment.aln.gz",
    "core_gene_alignment.fa.gz",
    "core_gene_alignment.fasta.gz",
)


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


def _find_panaroo_alignment(output_dir: Path) -> Optional[Path]:
    out = _find_panaroo_out(output_dir)
    if out is None:
        return None
    for candidate in _PANAROO_ALIGNMENT_CANDIDATES:
        path = out / candidate
        if path.is_file():
            return path
    matches: List[Path] = []
    for candidate in _PANAROO_ALIGNMENT_CANDIDATES:
        matches.extend([p for p in out.rglob(candidate) if p.is_file()])
    if not matches:
        return None
    matches.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return matches[0]


def _open_text_auto(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return path.open("r", encoding="utf-8", errors="ignore")


@lru_cache(maxsize=8)
def _load_core_alignment_cached(
    path_str: str,
    _mtime_ns: int,
    _size: int,
) -> Tuple[Tuple[str, ...], np.ndarray]:
    path = Path(path_str)
    labels: List[str] = []
    sequences: List[np.ndarray] = []
    current_name: Optional[str] = None
    current_parts: List[str] = []

    def _flush_current() -> None:
        nonlocal current_name, current_parts
        if current_name is None:
            return
        sequence = "".join(current_parts).strip()
        if sequence:
            labels.append(current_name)
            sequences.append(np.frombuffer(sequence.encode("ascii", errors="ignore"), dtype=np.uint8))
        current_name = None
        current_parts = []

    with _open_text_auto(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                _flush_current()
                current_name = line[1:].strip()
                continue
            current_parts.append(line)
    _flush_current()

    if not labels:
        matrix = np.zeros((0, 0), dtype=np.uint8)
        matrix.setflags(write=False)
        return tuple(), matrix

    lengths = {int(seq.shape[0]) for seq in sequences}
    if len(lengths) != 1:
        raise ValueError(f"Alignment has inconsistent sequence lengths: {sorted(lengths)}")

    matrix = np.vstack(sequences).astype(np.uint8, copy=False)
    matrix.setflags(write=False)
    return tuple(labels), matrix


def _quote_newick_label(label: str) -> str:
    return "'" + str(label).replace("'", "''") + "'"


def _replace_newick_leaf_tokens(newick: str, token_to_label: Mapping[str, str]) -> str:
    pattern = re.compile(r"(?<=\(|,)([^:(),;]+)(?=[:),;])")

    def _repl(match: re.Match[str]) -> str:
        token = match.group(1)
        label = token_to_label.get(token)
        if label is None:
            return token
        return _quote_newick_label(label)

    return pattern.sub(_repl, newick)


@lru_cache(maxsize=64)
def _build_subset_tree_newick_cached(
    path_str: str,
    _mtime_ns: int,
    _size: int,
    selected_samples: Tuple[str, ...],
) -> str:
    if pdist is None or squareform is None or DistanceMatrix is None or nj is None:
        raise RuntimeError("Subset tree recomputation requires scipy and scikit-bio.")

    labels, matrix = _load_core_alignment_cached(path_str, _mtime_ns, _size)
    if matrix.shape[0] == 0:
        return ""

    selected_set = {str(sample).strip() for sample in selected_samples if str(sample).strip()}
    if not selected_set:
        return ""

    keep_idx = [idx for idx, label in enumerate(labels) if label in selected_set]
    if not keep_idx:
        return ""

    subset_labels = [labels[idx] for idx in keep_idx]
    if len(subset_labels) == 1:
        return f"{_quote_newick_label(subset_labels[0])};"

    subset_matrix = matrix[keep_idx]
    if len(subset_labels) == 2:
        distance = float(np.mean(subset_matrix[0] != subset_matrix[1]))
        branch = max(0.0, distance * 0.5)
        return (
            f"({_quote_newick_label(subset_labels[0])}:{branch:.10f},"
            f"{_quote_newick_label(subset_labels[1])}:{branch:.10f});"
        )

    token_ids = [f"leaf_{idx}" for idx in range(len(subset_labels))]
    distances = squareform(pdist(subset_matrix, metric="hamming"))
    dm = DistanceMatrix(distances, token_ids)
    tree = nj(dm)
    newick = str(tree).strip()
    return _replace_newick_leaf_tokens(newick, dict(zip(token_ids, subset_labels)))


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
