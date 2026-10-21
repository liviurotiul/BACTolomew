from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

import pandas as pd

from .constants import PLASMID_TRUE_VALUES
from .panaroo_common import _natural_sort_key
from .polars_io import read_csv_df
from pipeline.MobSuiteWorkflow import mobsuite_row_has_direct_plasmid_evidence


_EMPTY_VALUES = {"", "-", ".", "na", "n/a", "nan", "none", "null", "<na>"}
_GENERIC_GENE_LABELS = {
    "cds",
    "gene",
    "protein",
    "hypothetical protein",
    "putative protein",
    "uncharacterized protein",
    "unknown protein",
    "conserved hypothetical protein",
}


def _clean_text(value: Any) -> str:
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    text = str(value or "").strip()
    return "" if text.casefold() in _EMPTY_VALUES else text


def _normalized_text(value: Any) -> str:
    text = _clean_text(value).casefold()
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def _is_generated_unknown_id(value: Any) -> bool:
    text = _clean_text(value)
    if not text:
        return False
    if text.upper().startswith("UNK:"):
        return True
    return bool(re.fullmatch(r"(?i)UNK[\s:_-]+[a-f0-9]{6,}", text))


def _is_informative_gene_label(value: Any) -> bool:
    text = _clean_text(value)
    if not text or _is_generated_unknown_id(text):
        return False
    norm = _normalized_text(text)
    if len(norm) < 2 or norm in _GENERIC_GENE_LABELS:
        return False
    if norm.startswith("hypothetical") or norm.startswith("uncharacterized"):
        return False
    if norm.startswith("duf"):
        return False
    return True


def _truthy_text(value: Any) -> bool:
    return str(value or "").strip().casefold() in PLASMID_TRUE_VALUES


def _short_hash(*parts: Any, length: int = 8) -> str:
    digest = hashlib.sha1("|".join(_clean_text(part) for part in parts).encode("utf-8")).hexdigest()
    return digest[: max(4, int(length))]


def _split_values(value: Any) -> List[str]:
    text = _clean_text(value)
    if not text:
        return []
    return [
        token.strip()
        for token in re.split(r"[,;|]\s*", text)
        if _clean_text(token)
    ]


def _first_value(values: Iterable[Any]) -> str:
    for value in values:
        text = _clean_text(value)
        if text:
            return text
    return ""


def _compact_descriptor(value: Any) -> str:
    text = _clean_text(value)
    if not text:
        return ""
    token = _split_values(text)[0] if _split_values(text) else text
    token = token.strip()
    cluster_match = re.fullmatch(r"(?i)(rep|mob|mpf)_?cluster_?([A-Za-z0-9_.-]+)", token)
    if cluster_match:
        return f"{cluster_match.group(1).lower()}{cluster_match.group(2)}"
    token = re.sub(r"(?i)^rep_cluster_", "rep", token)
    return token[:24]


def plasmid_cluster_key_from_row(row: Mapping[str, Any]) -> str:
    primary = _clean_text(row.get("MobSuitePrimaryClusterID") or row.get("PrimaryClusterID"))
    if primary:
        return f"primary:{primary}"
    secondary = _clean_text(row.get("MobSuiteSecondaryClusterID") or row.get("SecondaryClusterID"))
    if secondary:
        return f"secondary:{secondary}"
    neighbor = _clean_text(row.get("MobSuiteMashNearestNeighbor") or row.get("MashNearestNeighbor"))
    md5 = _clean_text(row.get("MD5") or row.get("MobSuiteMD5"))
    sample = _clean_text(row.get("Sample"))
    plasmid_id = _clean_text(row.get("MobSuitePlasmidID") or row.get("PlasmidID"))
    if neighbor:
        return f"neighbor:{neighbor}"
    return f"unclustered:{_short_hash(sample, plasmid_id, md5, row.get('Contig'))}"


def plasmid_instance_key_from_row(row: Mapping[str, Any]) -> str:
    plasmid_id = _clean_text(row.get("MobSuitePlasmidID") or row.get("PlasmidID"))
    sample = _clean_text(row.get("Sample"))
    if plasmid_id:
        return plasmid_id
    cluster_key = plasmid_cluster_key_from_row(row)
    if sample and cluster_key:
        return f"{sample}:{cluster_key}"
    return cluster_key


def _display_base_for_cluster_key(cluster_key: str) -> str:
    prefix, _, value = str(cluster_key or "").partition(":")
    clean_value = _clean_text(value)
    if prefix == "primary" and clean_value:
        return f"p{clean_value}"
    if prefix == "secondary" and clean_value:
        return f"pSec-{clean_value}"
    if prefix == "neighbor" and clean_value:
        return f"pNN-{clean_value}"
    if clean_value:
        return f"pUnclustered-{clean_value[:8]}"
    return "pUnclustered"


def plain_plasmid_name(display_name: Any) -> str:
    """'p123 (IncFII, MOBF)' for the cluster name 'p123 | IncFII | MOBF'.

    Row labels and window titles use this form, so they read as words
    rather than as '|'-separated fields.
    """
    parts = [part.strip() for part in _clean_text(display_name).split("|") if part.strip()]
    if not parts:
        return ""
    return f"{parts[0]} ({', '.join(parts[1:])})" if len(parts) > 1 else parts[0]


@dataclass
class PlasmidClusterSummary:
    key: str
    display_name: str
    primary_cluster_id: str = ""
    secondary_cluster_ids: Set[str] = field(default_factory=set)
    samples: Set[str] = field(default_factory=set)
    plasmid_ids: Set[str] = field(default_factory=set)
    replicon_types: Set[str] = field(default_factory=set)
    relaxase_types: Set[str] = field(default_factory=set)
    mpf_types: Set[str] = field(default_factory=set)
    orit_types: Set[str] = field(default_factory=set)
    mobility: Set[str] = field(default_factory=set)
    mash_neighbors: Set[str] = field(default_factory=set)
    sizes: List[int] = field(default_factory=list)
    num_contigs: List[int] = field(default_factory=list)

    @property
    def prevalence(self) -> int:
        return len(self.samples)

    @property
    def min_size(self) -> int:
        return min(self.sizes) if self.sizes else 0

    @property
    def max_size(self) -> int:
        return max(self.sizes) if self.sizes else 0

    def tooltip(self) -> str:
        lines = [
            f"Plasmid cluster: {self.display_name}",
            f"Samples: {self.prevalence}",
        ]
        if self.primary_cluster_id:
            lines.append(f"MOB primary cluster: {self.primary_cluster_id}")
        if self.secondary_cluster_ids:
            lines.append("Secondary cluster(s): " + ", ".join(sorted(self.secondary_cluster_ids, key=_natural_sort_key)[:8]))
        if self.replicon_types:
            lines.append("Replicon: " + ", ".join(sorted(self.replicon_types, key=_natural_sort_key)[:8]))
        if self.relaxase_types:
            lines.append("Relaxase: " + ", ".join(sorted(self.relaxase_types, key=_natural_sort_key)[:8]))
        if self.mpf_types:
            lines.append("MPF: " + ", ".join(sorted(self.mpf_types, key=_natural_sort_key)[:8]))
        if self.mobility:
            lines.append("Mobility: " + ", ".join(sorted(self.mobility, key=_natural_sort_key)[:8]))
        if self.min_size or self.max_size:
            if self.min_size == self.max_size:
                lines.append(f"Size: {self.max_size:,} bp")
            else:
                lines.append(f"Size range: {self.min_size:,}-{self.max_size:,} bp")
        if self.mash_neighbors:
            lines.append("Mash neighbor: " + ", ".join(sorted(self.mash_neighbors, key=_natural_sort_key)[:5]))
        return "\n".join(lines)


class PlasmidClusterIndex:
    def __init__(self) -> None:
        self.clusters: Dict[str, PlasmidClusterSummary] = {}
        self.presence_by_sample: Dict[str, Set[str]] = {}
        self.instance_to_cluster: Dict[Tuple[str, str], str] = {}
        self._frozen = False

    def __bool__(self) -> bool:
        return bool(self.clusters)

    def freeze(self) -> "PlasmidClusterIndex":
        """Make the index read only (a run's shared index, window plan 4.3).

        The maps become ``MappingProxyType`` views, the per-sample cluster sets
        frozensets, and ``add_row`` raises. The cluster summaries stay plain
        objects: readers treat them as read only.
        """
        if not self._frozen:
            self.clusters = MappingProxyType(dict(self.clusters))  # type: ignore[assignment]
            self.presence_by_sample = MappingProxyType(  # type: ignore[assignment]
                {sample: frozenset(keys) for sample, keys in self.presence_by_sample.items()}
            )
            self.instance_to_cluster = MappingProxyType(dict(self.instance_to_cluster))  # type: ignore[assignment]
            self._frozen = True
        return self

    def add_row(self, row: Mapping[str, Any]) -> None:
        if self._frozen:
            raise RuntimeError("PlasmidClusterIndex is frozen (read only).")
        sample = _clean_text(row.get("Sample"))
        cluster_key = plasmid_cluster_key_from_row(row)
        if not sample or not cluster_key:
            return
        primary = _clean_text(row.get("MobSuitePrimaryClusterID") or row.get("PrimaryClusterID"))
        display_name = self._make_display_name(cluster_key, row)
        summary = self.clusters.get(cluster_key)
        if summary is None:
            summary = PlasmidClusterSummary(
                key=cluster_key,
                display_name=display_name,
                primary_cluster_id=primary,
            )
            self.clusters[cluster_key] = summary
        elif len(display_name) > len(summary.display_name):
            summary.display_name = display_name
        if primary and not summary.primary_cluster_id:
            summary.primary_cluster_id = primary

        summary.samples.add(sample)
        self.presence_by_sample.setdefault(sample, set()).add(cluster_key)

        plasmid_id = _clean_text(row.get("MobSuitePlasmidID") or row.get("PlasmidID"))
        if plasmid_id:
            summary.plasmid_ids.add(plasmid_id)
            self.instance_to_cluster[(sample, plasmid_id)] = cluster_key
        instance_key = plasmid_instance_key_from_row(row)
        if instance_key:
            self.instance_to_cluster[(sample, instance_key)] = cluster_key

        secondary = _clean_text(row.get("MobSuiteSecondaryClusterID") or row.get("SecondaryClusterID"))
        if secondary:
            summary.secondary_cluster_ids.add(secondary)
        for attr, column_names in (
            ("replicon_types", ("MobSuiteRepliconTypes", "RepliconTypes")),
            ("relaxase_types", ("MobSuiteRelaxaseTypes", "RelaxaseTypes")),
            ("mpf_types", ("MobSuiteMPFType", "MPFType")),
            ("orit_types", ("MobSuiteOriTTypes", "OriTTypes")),
            ("mobility", ("MobSuitePredictedMobility", "PredictedMobility")),
            ("mash_neighbors", ("MobSuiteMashNearestNeighbor", "MashNearestNeighbor")),
        ):
            value = _first_value(row.get(column) for column in column_names)
            getattr(summary, attr).update(_split_values(value))
        try:
            size = int(float(_first_value(row.get(column) for column in ("MobSuiteLength", "Size"))))
            if size > 0:
                summary.sizes.append(size)
        except Exception:
            pass
        try:
            count = int(float(_clean_text(row.get("NumContigs"))))
            if count > 0:
                summary.num_contigs.append(count)
        except Exception:
            pass

    def _make_display_name(self, cluster_key: str, row: Mapping[str, Any]) -> str:
        base = _display_base_for_cluster_key(cluster_key)
        descriptors = []
        replicon = _compact_descriptor(_first_value(row.get(column) for column in ("MobSuiteRepliconTypes", "RepliconTypes")))
        relaxase = _compact_descriptor(_first_value(row.get(column) for column in ("MobSuiteRelaxaseTypes", "RelaxaseTypes")))
        if replicon:
            descriptors.append(replicon)
        if relaxase and relaxase not in descriptors:
            descriptors.append(relaxase)
        return " | ".join([base, *descriptors[:2]])

    def display_name_for_cluster(self, cluster_key: str) -> str:
        summary = self.clusters.get(cluster_key)
        return summary.display_name if summary is not None else _display_base_for_cluster_key(cluster_key)

    def tooltip_for_cluster(self, cluster_key: str) -> str:
        summary = self.clusters.get(cluster_key)
        return summary.tooltip() if summary is not None else self.display_name_for_cluster(cluster_key)

    def display_name_for_row(self, row: Mapping[str, Any]) -> str:
        return self.display_name_for_cluster(plasmid_cluster_key_from_row(row))

    def cluster_order(self, sample_filter: Optional[Set[str]] = None) -> List[str]:
        if sample_filter:
            allowed = {str(sample).strip() for sample in sample_filter if str(sample).strip()}
            keys = {
                key
                for sample, values in self.presence_by_sample.items()
                if sample in allowed
                for key in values
            }
        else:
            keys = set(self.clusters)
        return sorted(
            keys,
            key=lambda key: (
                -len(self.clusters.get(key, PlasmidClusterSummary(key, "")).samples),
                _natural_sort_key(self.display_name_for_cluster(key)),
            ),
        )


def _read_tsv_if_present(path: Path) -> pd.DataFrame:
    if not path.is_file():
        return pd.DataFrame()
    try:
        return read_csv_df(path, sep="\t", dtype=str).fillna("")
    except Exception:
        try:
            return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False).fillna("")
        except Exception:
            return pd.DataFrame()


def build_plasmid_cluster_index(
    *,
    output_dir: Optional[Path] = None,
    ann_df: Optional[pd.DataFrame] = None,
) -> PlasmidClusterIndex:
    index = PlasmidClusterIndex()
    plasmids_df = pd.DataFrame()
    if output_dir is not None:
        plasmids_df = _read_tsv_if_present(Path(output_dir) / "mobsuite" / "plasmids.tsv")
    if not plasmids_df.empty:
        for row in plasmids_df.to_dict(orient="records"):
            if not mobsuite_row_has_direct_plasmid_evidence(row):
                continue
            index.add_row(row)

    if ann_df is not None and not ann_df.empty and "Sample" in ann_df.columns:
        cols = [
            col
            for col in (
                "Sample",
                "Contig",
                "MobSuiteIsPlasmid",
                "MobSuitePlasmidID",
                "MobSuitePrimaryClusterID",
                "MobSuiteSecondaryClusterID",
                "MobSuiteRepliconTypes",
                "MobSuiteRelaxaseTypes",
                "MobSuiteMPFType",
                "MobSuiteOriTTypes",
                "MobSuitePredictedMobility",
                "MobSuiteMashNearestNeighbor",
                "MobSuiteMashNeighborDistance",
                "MobSuiteLength",
            )
            if col in ann_df.columns
        ]
        if cols and "MobSuiteIsPlasmid" in cols:
            subset = ann_df.loc[:, cols].copy()
            plasmid_mask = subset["MobSuiteIsPlasmid"].map(_truthy_text)
            subset = subset.loc[plasmid_mask].drop_duplicates()
            for row in subset.to_dict(orient="records"):
                index.add_row(row)
    return index


def _gene_label(row: Mapping[str, Any]) -> str:
    for col in ("Gene", "Product", "ConsensusDonorProduct", "GeneID"):
        value = _clean_text(row.get(col))
        if _is_informative_gene_label(value):
            return value
    return "gene"


def _gene_key(row: Mapping[str, Any]) -> str:
    for col in ("PanarooCluster", "PanarooGeneCluster", "PanarooGene", "PanarooID", "GeneCluster", "Cluster"):
        value = _clean_text(row.get(col))
        if value:
            return f"cluster:{value.casefold()}"
    for col in ("Gene", "Product", "ConsensusDonorProduct", "GeneID"):
        value = _clean_text(row.get(col))
        if _is_informative_gene_label(value):
            return f"{col.casefold()}:{_normalized_text(value)}"
    aa_seq = _clean_text(row.get("Amino_Acid_Sequence")).replace("*", "")
    if len(aa_seq) >= 20:
        digest = hashlib.sha1(aa_seq.encode("utf-8")).hexdigest()[:16]
        return f"aa:{digest}"
    gene_id = _clean_text(row.get("GeneID"))
    if gene_id:
        return f"locus:{gene_id.casefold()}"
    return f"locus:{_short_hash(row.get('Sample'), row.get('Contig'), row.get('Start'), row.get('Stop'))}"


def _gene_color(row: Mapping[str, Any]) -> str:
    tags = _clean_text(row.get("Tags"))
    if "A" in tags or "X" in tags:  # X: resistance mutation, drawn as AMR
        return "#ef9a9a"
    if "B" in tags:
        return "#a5d6a7"
    if "V" in tags:
        return "#ce93d8"
    if "M" in tags:
        return "#ffcc80"
    if "P" in tags:
        return "#9ec3ff"
    if "I" in tags:
        return "#ffe082"
    return "#90caf9"


def _safe_int(value: Any) -> Optional[int]:
    try:
        if pd.isna(value):
            return None
    except Exception:
        pass
    try:
        return int(float(str(value).strip()))
    except Exception:
        return None


def _contig_sort_key(value: Any) -> Tuple[Tuple[int, object], ...]:
    text = _clean_text(value)
    return _natural_sort_key(text)


def _plasmid_instance_genes(df: pd.DataFrame) -> List[Dict[str, Any]]:
    if df.empty:
        return []
    rows: List[Dict[str, Any]] = []
    offset = 1
    gap = 1000
    for _contig, contig_df in sorted(df.groupby("Contig", dropna=False), key=lambda item: _contig_sort_key(item[0])):
        ordered = contig_df.copy()
        start_values = ordered["StartInt"] if "StartInt" in ordered.columns else ordered.get("Start")
        stop_values = ordered["StopInt"] if "StopInt" in ordered.columns else ordered.get("Stop")
        ordered = ordered.assign(
            __start=pd.to_numeric(start_values, errors="coerce"),
            __stop=pd.to_numeric(stop_values, errors="coerce"),
        )
        ordered = ordered.dropna(subset=["__start", "__stop"]).sort_values(["__start", "__stop"], kind="stable")
        contig_max = 0
        for _idx, raw_row in ordered.iterrows():
            row = raw_row.to_dict()
            start = _safe_int(row.get("__start"))
            stop = _safe_int(row.get("__stop"))
            if start is None or stop is None:
                continue
            left = min(start, stop)
            right = max(start, stop)
            contig_max = max(contig_max, right)
            rows.append(
                {
                    "start": int(offset + left),
                    "end": int(offset + right),
                    "strand": _clean_text(row.get("Strand")) or "+",
                    "label": _gene_label(row),
                    "color": _gene_color(row),
                    "row_data": row,
                }
            )
        offset += max(1, contig_max) + gap
    return rows


def _instance_similarity(a: Sequence[str], b: Sequence[str]) -> float:
    set_a = set(a)
    set_b = set(b)
    if not set_a and not set_b:
        return 0.0
    inter = len(set_a & set_b)
    union = len(set_a | set_b) or 1
    jaccard = inter / union
    order_bonus = 0.0
    if a and b and inter:
        b_pos = {value: idx for idx, value in enumerate(b)}
        projected = [b_pos[value] for value in a if value in b_pos]
        if len(projected) >= 2:
            ordered_pairs = sum(1 for left, right in zip(projected, projected[1:]) if left <= right)
            order_bonus = ordered_pairs / max(1, len(projected) - 1)
    return (0.75 * jaccard) + (0.25 * order_bonus)


def build_plasmid_representative_context_rows(
    ann_df: pd.DataFrame,
    index: PlasmidClusterIndex,
    *,
    sample_filter: Optional[Set[str]] = None,
) -> List[Dict[str, Any]]:
    if ann_df is None or ann_df.empty or not index:
        return []
    required = {"Sample", "Contig", "MobSuiteIsPlasmid"}
    if not required.issubset(set(ann_df.columns)):
        return []
    sample_filter = {str(sample).strip() for sample in sample_filter or set() if str(sample).strip()} or None
    df = ann_df.loc[ann_df["MobSuiteIsPlasmid"].map(_truthy_text)].copy()
    if sample_filter is not None:
        df = df.loc[df["Sample"].fillna("").astype(str).str.strip().isin(sample_filter)].copy()
    if df.empty:
        return []
    records = df.to_dict(orient="records")
    cluster_values = [plasmid_cluster_key_from_row(row) for row in records]
    instance_values = [plasmid_instance_key_from_row(row) for row in records]
    df["__plasmid_cluster_key"] = cluster_values
    df["__plasmid_instance_key"] = instance_values
    df = df[df["__plasmid_cluster_key"].astype(str).str.strip().ne("")]

    instances_by_cluster: Dict[str, List[Dict[str, Any]]] = {}
    for (cluster_key, sample, instance_key), group in df.groupby(
        ["__plasmid_cluster_key", "Sample", "__plasmid_instance_key"],
        dropna=False,
    ):
        cluster_key = _clean_text(cluster_key)
        sample = _clean_text(sample)
        instance_key = _clean_text(instance_key)
        if not cluster_key or not sample:
            continue
        genes = _plasmid_instance_genes(group)
        if not genes:
            continue
        gene_keys = [_gene_key(gene.get("row_data", {})) for gene in genes]
        instances_by_cluster.setdefault(cluster_key, []).append(
            {
                "sample": sample,
                "instance_key": instance_key,
                "genes": genes,
                "gene_keys": gene_keys,
                "gene_count": len(genes),
            }
        )

    out_rows: List[Dict[str, Any]] = []
    for cluster_key in index.cluster_order(sample_filter):
        instances = instances_by_cluster.get(cluster_key) or []
        if not instances:
            continue
        if len(instances) == 1:
            representative = instances[0]
        else:
            scored: List[Tuple[float, int, Dict[str, Any]]] = []
            for idx, candidate in enumerate(instances):
                others = [other for other in instances if other is not candidate]
                avg = sum(
                    _instance_similarity(candidate["gene_keys"], other["gene_keys"])
                    for other in others
                ) / max(1, len(others))
                scored.append((avg, int(candidate.get("gene_count") or 0), candidate))
            scored.sort(key=lambda item: (item[0], item[1]), reverse=True)
            representative = scored[0][2]
        summary = index.clusters.get(cluster_key)
        tooltip = index.tooltip_for_cluster(cluster_key)
        tooltip += f"\nRepresentative: {representative['sample']}"
        plain_name = plain_plasmid_name(index.display_name_for_cluster(cluster_key))
        out_rows.append(
            {
                "sample": representative["sample"],
                "contig": representative.get("instance_key") or cluster_key,
                "display_label": f"{representative['sample']}, plasmid {plain_name}",
                "title": f"{representative['sample']}, plasmid {plain_name}",
                "molecule": "plasmid",
                "is_plasmid": True,
                # Gene positions are offsets into the plasmid's contigs joined
                # end to end, not contig coordinates.
                "joined_contigs": True,
                "meta_tooltip": tooltip,
                "cluster_key": cluster_key,
                "prevalence": summary.prevalence if summary is not None else 1,
                "genes": representative["genes"],
            }
        )
    return out_rows


def build_plasmid_cluster_context_rows(
    ann_df: pd.DataFrame,
    index: PlasmidClusterIndex,
    cluster_key: str,
    *,
    sample_filter: Optional[Set[str]] = None,
) -> List[Dict[str, Any]]:
    if ann_df is None or ann_df.empty or not index:
        return []
    cluster_key = _clean_text(cluster_key)
    if not cluster_key:
        return []
    required = {"Sample", "Contig", "MobSuiteIsPlasmid"}
    if not required.issubset(set(ann_df.columns)):
        return []
    sample_filter = {str(sample).strip() for sample in sample_filter or set() if str(sample).strip()} or None
    df = ann_df.loc[ann_df["MobSuiteIsPlasmid"].map(_truthy_text)].copy()
    if sample_filter is not None:
        df = df.loc[df["Sample"].fillna("").astype(str).str.strip().isin(sample_filter)].copy()
    if df.empty:
        return []
    records = df.to_dict(orient="records")
    df["__plasmid_cluster_key"] = [plasmid_cluster_key_from_row(row) for row in records]
    df["__plasmid_instance_key"] = [plasmid_instance_key_from_row(row) for row in records]
    df = df[df["__plasmid_cluster_key"].astype(str).str.strip() == cluster_key]
    if df.empty:
        return []

    rows: List[Dict[str, Any]] = []
    display_name = plain_plasmid_name(index.display_name_for_cluster(cluster_key))
    for (sample, instance_key), group in df.groupby(["Sample", "__plasmid_instance_key"], dropna=False):
        sample = _clean_text(sample)
        instance_key = _clean_text(instance_key)
        if not sample:
            continue
        genes = _plasmid_instance_genes(group)
        if not genes:
            continue
        rows.append(
            {
                "sample": sample,
                "contig": instance_key or cluster_key,
                "display_label": f"{sample}, plasmid {display_name}",
                "title": f"{sample}, plasmid {display_name}",
                "molecule": "plasmid",
                "is_plasmid": True,
                # Gene positions are offsets into the plasmid's contigs joined
                # end to end, not contig coordinates.
                "joined_contigs": True,
                "meta_tooltip": index.tooltip_for_cluster(cluster_key),
                "cluster_key": cluster_key,
                "genes": genes,
            }
        )
    rows.sort(key=lambda row: _natural_sort_key(str(row.get("sample") or "")))
    return rows
