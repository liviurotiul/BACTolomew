from __future__ import annotations

import hashlib
import math
import re
from dataclasses import dataclass
from typing import Any, Mapping, Optional, Sequence, Tuple

from .resistance_mutation import has_resistance_mutation, mutation_note


EMPTY_TEXT_VALUES = {"", "-", ".", "na", "n/a", "nan", "none", "null", "<na>"}
GENERIC_LABELS = {
    "cds",
    "gene",
    "protein",
    "hypothetical protein",
    "putative protein",
    "unknown protein",
    "uncharacterized protein",
    "conserved hypothetical protein",
}
GENERIC_MARKER_TOKENS = {
    "amr",
    "amrf",
    "amrfinder",
    "bakta",
    "blast",
    "card",
    "cge mge",
    "isfinder",
    "marker",
    "mobile",
    "protein",
    "vfdb",
    "virulence",
}

# A deterministic, high-contrast categorical palette.  The fill encodes a
# homologous family; interaction state is always carried by the outline.
FAMILY_PALETTE: Tuple[str, ...] = (
    "#4e79a7",
    "#f28e2b",
    "#59a14f",
    "#e15759",
    "#76b7b2",
    "#b07aa1",
    "#edc948",
    "#ff9da7",
    "#9c755f",
    "#bab0ab",
    "#2f6b9a",
    "#d47418",
    "#3f8f55",
    "#c84d59",
    "#568e8b",
    "#9561a1",
    "#b89b20",
    "#d97b91",
)

# Gene type colors: the same AMR, BMRG, virulence and mobile colors as the
# annotation table's tag chips (annotation_table._tag_color).  Everything else
# is a quiet grey so the typed genes stand out around the focal gene.
CATEGORY_COLORS = {
    "AMR": "#c62828",
    "VIRULENCE": "#6a1b9a",
    "MOBILE": "#ff8f00",
    "BMRG": "#1b5e20",
    "OTHER": "#c5ced8",
}
CATEGORY_LABELS = {
    "AMR": "AMR",
    "VIRULENCE": "Virulence",
    "MOBILE": "Mobile",
    "BMRG": "BMRG",
    "OTHER": "Other",
}


def clean_text(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() in EMPTY_TEXT_VALUES else text


def truthy_value(value: Any) -> bool:
    """Return a conservative bool without asking pandas missing values to coerce."""
    text = clean_text(value).lower()
    return text in {"1", "1.0", "true", "yes", "y", "on", "plasmid"}


def normalized_color_text(value: Any) -> str:
    if value is None:
        return ""
    name_method = getattr(value, "name", None)
    if callable(name_method):
        try:
            named = str(name_method()).strip()
            if named:
                return named
        except Exception:
            pass
    text = clean_text(value)
    # Avoid treating a Python object representation as a CSS colour.
    return "" if text.startswith("<") and text.endswith(">") else text


def normalized_text(value: Any) -> str:
    text = clean_text(value).lower()
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def informative_text(value: Any) -> bool:
    text = normalized_text(value)
    if len(text) < 2 or text in GENERIC_LABELS:
        return False
    return not (
        text.startswith("hypothetical")
        or text.startswith("uncharacterized")
        or text.startswith("duf")
    )


def canonical_marker_symbol(value: Any) -> str:
    text = clean_text(value)
    if not text:
        return ""
    for pattern in (
        r"\bintI\d+[A-Za-z0-9_-]*\b",
        r"\bIS(?:[A-Za-z]+[0-9][A-Za-z0-9_-]*|[0-9]+[A-Za-z0-9_-]*)\b",
        r"\bTn(?:[A-Za-z]+[0-9][A-Za-z0-9_-]*|[0-9]+[A-Za-z0-9_-]*)\b",
    ):
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match and match.group(0).lower() != "isfinder":
            return match.group(0)
    return ""


def first_informative_marker_token(value: Any) -> str:
    text = clean_text(value)
    if not text:
        return ""
    for token in re.split(r"[,;|]\s*", text):
        token = token.strip()
        normalized = normalized_text(token)
        if (
            normalized
            and normalized not in GENERIC_MARKER_TOKENS
            and informative_text(token)
        ):
            return token
    return ""


def _payload_for_gene(gene: Mapping[str, Any]) -> Mapping[str, Any]:
    payload = gene.get("row_data")
    return payload if isinstance(payload, dict) else {}


def homology_key_for_gene(gene: Mapping[str, Any]) -> str:
    payload = _payload_for_gene(gene)
    for column in (
        "PanarooCluster",
        "PanarooGeneCluster",
        "PanarooGene",
        "PanarooID",
        "GeneCluster",
        "Cluster",
    ):
        value = clean_text(payload.get(column))
        if informative_text(value):
            return f"cluster:{normalized_text(value)}"

    for column in ("Gene", "GeneID"):
        value = clean_text(payload.get(column))
        if value.upper().startswith("UNK:"):
            continue
        marker = canonical_marker_symbol(value)
        if marker:
            return f"gene:{normalized_text(marker)}"
        if informative_text(value):
            return f"gene:{normalized_text(value)}"

    for column in ("MarkerSearchText", "ConsensusDonorProduct", "Product"):
        marker = canonical_marker_symbol(payload.get(column))
        if marker:
            return f"gene:{normalized_text(marker)}"

    marker_token = first_informative_marker_token(payload.get("MarkerSearchText"))
    if marker_token:
        return f"marker:{normalized_text(marker_token)}"

    for column in ("ConsensusDonorProduct", "Product"):
        value = clean_text(payload.get(column))
        if informative_text(value):
            return f"product:{normalized_text(value)}"

    amino_acids = clean_text(payload.get("Amino_Acid_Sequence")).replace("*", "")
    if len(amino_acids) >= 20:
        digest = hashlib.sha1(amino_acids.encode("utf-8")).hexdigest()[:16]
        return f"aa:{digest}"

    label = clean_text(gene.get("label"))
    marker = canonical_marker_symbol(label)
    if marker:
        return f"gene:{normalized_text(marker)}"
    if informative_text(label):
        return f"label:{normalized_text(label)}"
    return ""


def homology_evidence(key: str) -> Tuple[str, float]:
    prefix = str(key or "").split(":", 1)[0].lower()
    return {
        "cluster": ("Panaroo/cluster", 1.30),
        "aa": ("protein sequence", 1.18),
        "gene": ("gene", 1.00),
        "marker": ("database annotation", 0.92),
        "product": ("product", 0.72),
        "label": ("label", 0.58),
    }.get(prefix, ("unknown", 0.45))


def annotation_category(payload: Mapping[str, Any]) -> str:
    aliases = {
        "amr": "AMR",
        "arg": "AMR",
        "args": "AMR",
        "core amr": "AMR",
        "plus amr": "AMR",
        "bmrg": "BMRG",
        "core bmrg": "BMRG",
        "plus bmrg": "BMRG",
        "virulence": "VIRULENCE",
        "plus virulence": "VIRULENCE",
        "mobile": "MOBILE",
    }
    # A known resistance mutation (tag X) makes the gene read as AMR, whatever
    # the gene is otherwise (gyrA, parC, ...); the views add the badge.
    if has_resistance_mutation(payload):
        return "AMR"
    # Explicit biological categories are authoritative.  Source names such as
    # CARD describe evidence and must not turn a BMRG/mobile/virulence call into
    # AMR merely because resistance metadata is also present.
    for column in ("ConsensusCategory", "MarkerCategory", "Category"):
        for token in re.split(r"[,;|]\s*", clean_text(payload.get(column))):
            category = aliases.get(normalized_text(token))
            if category:
                return category
    # A phage or plasmid gene has no gene type color; its source and product
    # names must not make it read as AMR ("large terminase" holds "arg").
    if normalized_text(clean_text(payload.get("ConsensusCategory"))) in {"phage", "plasmid"}:
        return ""

    legacy_tags = clean_text(payload.get("Tags")).upper()
    for tag, category in (("A", "AMR"), ("B", "BMRG"), ("V", "VIRULENCE"), ("M", "MOBILE")):
        if tag in legacy_tags:
            return category

    values = [
        clean_text(payload.get(column))
        for column in (
            "Predictors",
            "MarkerDisplay",
            "ConsensusSource",
            "ResistanceMetadataSource",
            "MarkerSearchText",
        )
        if clean_text(payload.get(column))
    ]
    blob = normalized_text(" ".join(values))
    if any(token in blob for token in ("amr", "arg", "card", "amrfinder", "resistance")):
        return "AMR"
    if any(token in blob for token in ("virulence", "vfdb")):
        return "VIRULENCE"
    if any(token in blob for token in ("mobile", "isfinder", "mge", "transpos")):
        return "MOBILE"
    if any(token in blob for token in ("bmrg", "biocide", "metal")):
        return "BMRG"
    return ""


def gene_type_fill(payload: Mapping[str, Any]) -> str:
    """The gene type color of an annotation row (the neighborhood's "Gene
    type" coloring, also used by the genome track)."""
    return CATEGORY_COLORS.get(annotation_category(payload), CATEGORY_COLORS["OTHER"])


def display_label_for_gene(gene: Mapping[str, Any]) -> str:
    payload = _payload_for_gene(gene)
    candidates = (
        payload.get("Gene"),
        gene.get("label"),
        payload.get("ConsensusDonorProduct"),
        payload.get("Product"),
        payload.get("GeneID"),
    )
    for candidate in candidates:
        value = clean_text(candidate)
        if value.upper().startswith("UNK:"):
            continue
        marker = canonical_marker_symbol(candidate)
        if marker:
            return marker
        if informative_text(value):
            return value if len(value) <= 34 else value[:31].rstrip() + "…"
    return ""


def family_color(key: str) -> str:
    if not key:
        return "#c5ced8"
    digest = hashlib.sha1(key.encode("utf-8")).digest()
    return FAMILY_PALETTE[int.from_bytes(digest[:2], "big") % len(FAMILY_PALETTE)]


def molecule_kind(row: Mapping[str, Any]) -> str:
    """'plasmid', 'chromosome', 'mixed' or 'unknown' (no plasmid prediction)."""
    explicit = clean_text(row.get("molecule")) or clean_text(row.get("molecule_kind"))
    explicit = explicit.lower()
    if explicit in {"plasmid", "chromosome", "mixed", "unknown"}:
        return explicit
    if "is_plasmid" in row or "MobSuiteIsPlasmid" in row:
        value = row.get("is_plasmid", row.get("MobSuiteIsPlasmid", False))
        return "plasmid" if truthy_value(value) else "chromosome"
    states = []
    genes = row.get("genes")
    if isinstance(genes, list):
        for gene in genes:
            if not isinstance(gene, dict):
                continue
            payload = gene.get("row_data")
            if isinstance(payload, dict) and "MobSuiteIsPlasmid" in payload:
                states.append(truthy_value(payload.get("MobSuiteIsPlasmid")))
    if states:
        if all(states):
            return "plasmid"
        if any(states):
            return "mixed"
    return "unknown"


def optional_count(value: Any) -> Optional[int]:
    """A non-negative gene count, or None when the row does not say."""
    if isinstance(value, bool) or not clean_text(value):
        return None
    try:
        count = int(float(value))
    except (TypeError, ValueError, OverflowError):
        return None
    return count if count >= 0 else None


def normalized_strand(value: Any) -> str:
    # A literal '-' is a biological value here, not the missing-value token it
    # may represent in annotation text fields.
    text = "" if value is None else str(value).strip().lower()
    return "-" if text in {"-", "-1", "reverse", "rev"} else "+"


@dataclass(frozen=True)
class ContextGene:
    index: int
    start: float
    end: float
    strand: str
    label: str
    family_key: str
    evidence: str
    evidence_weight: float
    family_fill: str
    annotation_fill: str
    category: str
    is_focal: bool
    label_priority: int
    payload: Mapping[str, Any]
    source: Mapping[str, Any]
    # The resistance mutation badge's sentence ('' when the gene has none).
    mutation_note: str = ""

    @property
    def low(self) -> float:
        return min(self.start, self.end)

    @property
    def high(self) -> float:
        return max(self.start, self.end)

    @property
    def center(self) -> float:
        return (self.low + self.high) / 2.0


@dataclass(frozen=True)
class ContextRow:
    index: int
    sample: str
    contig: str
    title: str
    molecule: str
    similarity: Optional[float]
    comparison_orientation: str
    is_reference: bool
    sample_count: int
    samples: Tuple[str, ...]
    center_bp: Optional[float]
    genes: Tuple[ContextGene, ...]
    source: Mapping[str, Any]
    # The contig really ends after the first (low) or last (high) gene shown,
    # in contig coordinates.  Only set when the rows say how many genes the
    # contig holds on each side of the focal gene.
    ends_low: bool = False
    ends_high: bool = False
    # The contig ends inside the chosen flank on at least one side.
    partial: bool = False

    @property
    def has_neighbors(self) -> bool:
        return any(not gene.is_focal for gene in self.genes)

    @property
    def is_short(self) -> bool:
        """Only the focal gene: the contig is too short to compare.

        The reference row is never short, so it is never folded away.
        """
        return (
            not self.is_reference
            and any(gene.is_focal for gene in self.genes)
            and not self.has_neighbors
        )

    @property
    def focal_gene(self) -> Optional[ContextGene]:
        for gene in self.genes:
            if gene.is_focal:
                return gene
        if self.center_bp is not None and self.genes:
            return min(self.genes, key=lambda gene: abs(gene.center - float(self.center_bp)))
        return self.genes[len(self.genes) // 2] if self.genes else None

    @property
    def anchor_bp(self) -> float:
        if self.center_bp is not None:
            return float(self.center_bp)
        focal = self.focal_gene
        if focal is not None:
            return focal.center
        if self.genes:
            return (min(g.low for g in self.genes) + max(g.high for g in self.genes)) / 2.0
        return 0.0


@dataclass(frozen=True)
class ContextRenderModel:
    rows: Tuple[ContextRow, ...]

    @property
    def reference_row(self) -> Optional[ContextRow]:
        for row in self.rows:
            if row.is_reference:
                return row
        return self.rows[0] if self.rows else None

    @classmethod
    def from_rows(cls, rows: Sequence[Mapping[str, Any]]) -> "ContextRenderModel":
        normalized_rows = []
        for row_index, raw_row in enumerate(rows):
            if not isinstance(raw_row, dict):
                continue
            genes = []
            raw_genes = raw_row.get("genes")
            if not isinstance(raw_genes, (list, tuple)):
                raw_genes = []
            for gene_index, raw_gene in enumerate(raw_genes):
                if not isinstance(raw_gene, dict):
                    continue
                try:
                    start = float(raw_gene.get("start"))
                    end = float(raw_gene.get("end"))
                except (TypeError, ValueError):
                    continue
                if not (math.isfinite(start) and math.isfinite(end)):
                    continue
                if start == end:
                    end = start + 1.0
                payload = _payload_for_gene(raw_gene)
                key = homology_key_for_gene(raw_gene)
                evidence, weight = homology_evidence(key)
                category = annotation_category(payload)
                label = display_label_for_gene(raw_gene)
                focal = truthy_value(raw_gene.get("is_center"))
                priority = 100 if focal else 80 if category else 60 if label else 0
                annotation_fill = CATEGORY_COLORS.get(
                    category,
                    normalized_color_text(raw_gene.get("color")) or "#90caf9",
                )
                genes.append(
                    ContextGene(
                        index=gene_index,
                        start=start,
                        end=end,
                        strand=normalized_strand(raw_gene.get("strand")),
                        label=label,
                        family_key=key,
                        evidence=evidence,
                        evidence_weight=float(weight),
                        family_fill=family_color(key),
                        annotation_fill=annotation_fill,
                        category=category,
                        is_focal=focal,
                        label_priority=priority,
                        payload=dict(payload),
                        source=dict(raw_gene),
                        mutation_note=mutation_note(payload),
                    )
                )
            try:
                similarity: Optional[float] = float(raw_row.get("context_similarity"))
            except (TypeError, ValueError):
                similarity = None
            if similarity is not None and not math.isfinite(similarity):
                similarity = None
            try:
                center_bp: Optional[float] = float(raw_row.get("center_bp"))
            except (TypeError, ValueError):
                center_bp = None
            if center_bp is not None and not math.isfinite(center_bp):
                center_bp = None
            sample = clean_text(raw_row.get("sample"))
            contig = clean_text(raw_row.get("contig"))
            title = clean_text(raw_row.get("display_label")) or clean_text(raw_row.get("title"))
            if not title:
                title = sample or contig or f"Neighborhood {row_index + 1}"
            raw_samples = raw_row.get("samples")
            if not isinstance(raw_samples, (list, tuple, set)):
                raw_samples = raw_row.get("collapsed_samples")
            if isinstance(raw_samples, (list, tuple, set)):
                samples = tuple(dict.fromkeys(clean_text(value) for value in raw_samples if clean_text(value)))
            else:
                samples = (sample,) if sample else ()
            count_value = raw_row.get("sample_count")
            if not clean_text(count_value):
                count_value = raw_row.get("collapsed_sample_count")
            if not clean_text(count_value):
                count_value = len(samples) or 1
            try:
                sample_count = max(1, int(count_value))
            except (TypeError, ValueError, OverflowError):
                sample_count = max(1, len(samples))
            # 'is_comparison_reference' is set by the dialog for every row it
            # scores, and wins over 'is_source' so another row can be chosen
            # as the reference.
            if "is_comparison_reference" in raw_row:
                is_reference = truthy_value(raw_row.get("is_comparison_reference"))
            else:
                is_reference = truthy_value(raw_row.get("is_source"))
            normalized_rows.append(
                ContextRow(
                    index=row_index,
                    sample=sample,
                    contig=contig,
                    title=title,
                    molecule=molecule_kind(raw_row),
                    similarity=similarity,
                    comparison_orientation=clean_text(raw_row.get("context_orientation")).lower(),
                    is_reference=is_reference,
                    sample_count=sample_count,
                    samples=samples,
                    center_bp=center_bp,
                    genes=tuple(genes),
                    source=dict(raw_row),
                    ends_low=truthy_value(raw_row.get("contig_ends_low")),
                    ends_high=truthy_value(raw_row.get("contig_ends_high")),
                    partial=truthy_value(raw_row.get("context_partial")),
                )
            )
        return cls(rows=tuple(normalized_rows))
