from __future__ import annotations

import hashlib
import html
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from PyQt6 import QtCore, QtGui, QtWidgets

from .constants import MISSING_FILTER_TOKEN, PLASMID_TRUE_VALUES
from .metadata_filter_dialog import MetadataValueDialog
from .panaroo_views import PanCoreAccessoryTab


def _truthy_text(value: Any) -> bool:
    text = str(value or "").strip().lower()
    return text in PLASMID_TRUE_VALUES


def _row_molecule_kind(row: Dict[str, Any]) -> str:
    explicit = str(row.get("molecule") or row.get("molecule_kind") or "").strip().lower()
    if explicit in {"plasmid", "chromosome", "mixed"}:
        return explicit
    if "is_plasmid" in row:
        return "plasmid" if _truthy_text(row.get("is_plasmid")) else "chromosome"
    if "MobSuiteIsPlasmid" in row:
        return "plasmid" if _truthy_text(row.get("MobSuiteIsPlasmid")) else "chromosome"

    states: List[bool] = []
    genes = row.get("genes")
    if isinstance(genes, list):
        for gene in genes:
            if not isinstance(gene, dict):
                continue
            payload = gene.get("row_data")
            if not isinstance(payload, dict) or "MobSuiteIsPlasmid" not in payload:
                continue
            states.append(_truthy_text(payload.get("MobSuiteIsPlasmid")))
    if states:
        if all(states):
            return "plasmid"
        if any(states):
            return "mixed"
    return "chromosome"


def _molecule_badge(row: Dict[str, Any]) -> Tuple[str, str]:
    kind = _row_molecule_kind(row)
    if kind == "plasmid":
        return "P", "#16a34a"
    if kind == "mixed":
        return "P/C", "#64748b"
    return "C", "#2563eb"


_EMPTY_TEXT_VALUES = {"", "-", ".", "na", "n/a", "nan", "none", "null", "<na>"}
_GENERIC_HOMOLOGY_LABELS = {
    "cds",
    "gene",
    "protein",
    "hypothetical protein",
    "putative protein",
    "uncharacterized protein",
    "unknown protein",
    "conserved hypothetical protein",
}
_GENERIC_MARKER_LABELS = {
    "amrf",
    "amrfinder",
    "bakta",
    "blast",
    "card",
    "cog",
    "diamond",
    "isfinder",
    "marker",
    "marker gene",
    "mobile",
    "protein",
    "vfdb",
}
_FIGURE_MARKER_CATEGORIES = {
    "AMR",
    "ARG",
    "ARGS",
    "BMRG",
    "CORE-AMR",
    "CORE-BMRG",
    "MOBILE",
    "PLUS-AMR",
    "PLUS-BMRG",
    "PLUS-VIRULENCE",
    "VIRULENCE",
}
_FIGURE_CATEGORY_ALIASES = {
    "AMR": "AMR",
    "ARG": "AMR",
    "ARGS": "AMR",
    "CORE-AMR": "AMR",
    "PLUS-AMR": "AMR",
    "BMRG": "BMRG",
    "CORE-BMRG": "BMRG",
    "PLUS-BMRG": "BMRG",
    "VIRULENCE": "VIRULENCE",
    "PLUS-VIRULENCE": "VIRULENCE",
    "MOBILE": "MOBILE",
}
_FIGURE_CATEGORY_SORT = {"AMR": 0, "BMRG": 1, "VIRULENCE": 2, "MOBILE": 3}
_FIGURE_CATEGORY_HUES = {
    "AMR": (356, 14),
    "BMRG": (138, 20),
    "VIRULENCE": (280, 18),
    "MOBILE": (40, 18),
}
_FIGURE_FAMILY_SHAPES = {
    "transposase": "triangle",
    "integrase": "diamond",
    "recombinase": "down_triangle",
    "replication": "pentagon",
    "mobilization": "hexagon",
    "conjugation": "plus",
    "beta_lactamase": "pentagon",
    "aminoglycoside": "hexagon",
    "folate": "diamond",
    "macrolide": "hexagon",
    "tetracycline": "plus",
    "quinolone": "down_triangle",
    "colistin": "pentagon",
    "efflux": "plus",
    "metal_resistance": "circle",
    "toxin": "hexagon",
    "adhesion": "circle",
    "secretion": "plus",
    "capsule": "hexagon",
    "siderophore": "diamond",
    "marker_other": "square",
}
_FIGURE_FAMILY_LABELS = {
    "transposase": "transposase / IS / Tn",
    "integrase": "integrase / integron",
    "recombinase": "recombinase / resolvase",
    "replication": "replication",
    "mobilization": "mobilization / relaxase",
    "conjugation": "conjugation / transfer",
    "beta_lactamase": "beta-lactamase",
    "aminoglycoside": "aminoglycoside resistance",
    "folate": "folate-pathway resistance",
    "macrolide": "macrolide resistance",
    "tetracycline": "tetracycline resistance",
    "quinolone": "quinolone resistance",
    "colistin": "colistin/lipid A resistance",
    "efflux": "efflux/biocide resistance",
    "metal_resistance": "metal resistance",
    "toxin": "toxin",
    "adhesion": "adhesion",
    "secretion": "secretion system",
    "capsule": "capsule/surface",
    "siderophore": "siderophore/iron uptake",
    "marker_other": "other marker",
}


def _clean_context_text(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in _EMPTY_TEXT_VALUES:
        return ""
    return text


def _normalized_homology_text(value: Any) -> str:
    text = _clean_context_text(value).lower()
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def _is_informative_homology_text(value: Any) -> bool:
    norm = _normalized_homology_text(value)
    if len(norm) < 2 or norm in _GENERIC_HOMOLOGY_LABELS:
        return False
    if norm.startswith("hypothetical") or norm.startswith("uncharacterized"):
        return False
    if norm.startswith("duf"):
        return False
    return True


def _first_informative_marker_token(value: Any) -> str:
    text = _clean_context_text(value)
    if not text:
        return ""
    for token in re.split(r"[,;|]\s*", text):
        token = token.strip()
        if _is_informative_homology_text(token):
            return token
    return text if _is_informative_homology_text(text) else ""


def _canonical_marker_symbol(value: Any) -> str:
    text = _clean_context_text(value)
    if not text:
        return ""
    int_match = re.search(r"\bintI\d+[A-Za-z0-9_-]*\b", text, flags=re.IGNORECASE)
    if int_match:
        return int_match.group(0)
    for match in re.finditer(
        r"\bIS(?:[A-Za-z]+[0-9][A-Za-z0-9_-]*|[0-9]+[A-Za-z0-9_-]*)\b",
        text,
        flags=re.IGNORECASE,
    ):
        symbol = match.group(0)
        if symbol.lower() in {"isfinder"}:
            continue
        return symbol
    for match in re.finditer(
        r"\bTn(?:[A-Za-z]+[0-9][A-Za-z0-9_-]*|[0-9]+[A-Za-z0-9_-]*)\b",
        text,
        flags=re.IGNORECASE,
    ):
        return match.group(0)
    return ""


def _iter_context_value_tokens(value: Any) -> Sequence[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        out: List[str] = []
        for item in value:
            out.extend(_iter_context_value_tokens(item))
        return out
    text = _clean_context_text(value)
    if not text:
        return []
    return [token.strip() for token in re.split(r"[,;|]\s*", text) if token.strip()]


def _is_figure_marker_category(value: Any) -> bool:
    for token in _iter_context_value_tokens(value):
        normalized = re.sub(r"\s+", "-", token.strip().upper())
        if normalized in _FIGURE_MARKER_CATEGORIES:
            return True
    return False


def _canonical_figure_category(value: Any) -> str:
    for token in _iter_context_value_tokens(value):
        normalized = re.sub(r"\s+", "-", token.strip().upper())
        category = _FIGURE_CATEGORY_ALIASES.get(normalized)
        if category:
            return category
    return ""


def _best_figure_category(values: Sequence[Any]) -> str:
    categories = [
        category
        for value in values
        if (category := _canonical_figure_category(value))
    ]
    if not categories:
        return ""
    return min(categories, key=lambda category: _FIGURE_CATEGORY_SORT.get(category, 99))


def _is_useful_marker_shape_label(value: Any) -> bool:
    text = _clean_context_text(value)
    if _is_generated_marker_label(text):
        return False
    if not _is_informative_homology_text(text):
        return False
    norm = _normalized_homology_text(text)
    if norm in _GENERIC_MARKER_LABELS:
        return False
    if _is_figure_marker_category(text):
        return False
    return True


def _is_generated_marker_label(value: Any) -> bool:
    text = _clean_context_text(value)
    if not text:
        return False
    if text.upper().startswith("UNK:"):
        return True
    return bool(re.fullmatch(r"(?i)UNK[\s:_-]+[a-f0-9]{6,}", text))


def _first_useful_marker_shape_token(value: Any) -> str:
    for token in _iter_context_value_tokens(value):
        if _is_useful_marker_shape_label(token):
            return token
    text = _clean_context_text(value)
    return text if _is_useful_marker_shape_label(text) else ""


def _stable_hash_int(value: str) -> int:
    digest = hashlib.sha1(str(value or "").encode("utf-8")).hexdigest()
    return int(digest[:12], 16)


def _stable_figure_category_color(category: str, label: str) -> QtGui.QColor:
    base_hue, spread = _FIGURE_CATEGORY_HUES.get(category, (215, 10))
    seed = _stable_hash_int(f"{category}:{_normalized_homology_text(label)}")
    hue = int((base_hue + (seed % (spread * 2 + 1)) - spread) % 360)
    saturation = 166 + ((seed >> 9) % 58)
    value = 176 + ((seed >> 18) % 54)
    if category == "MOBILE":
        saturation = 178 + ((seed >> 9) % 50)
        value = 194 + ((seed >> 18) % 42)
    return QtGui.QColor.fromHsv(hue, saturation, value)


def _compact_marker_text(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", _clean_context_text(value).lower())


def _is_cge_mge_source(value: Any) -> bool:
    text = _clean_context_text(value).lower()
    return "cge-mge" in text or "cge mge" in _normalized_homology_text(text)


def _cge_mge_subject_parts(value: Any) -> Tuple[str, str]:
    text = _clean_context_text(value)
    if not text:
        return "", ""
    for chunk in re.split(r"[,;]\s*", text):
        parts = [part.strip() for part in chunk.split("|")]
        if len(parts) < 4:
            continue
        if not parts[0] or not parts[1].isdigit() or not parts[2]:
            continue
        element = parts[0]
        product = parts[3]
        return element, product
    return "", ""


def _cge_mge_product_label(value: Any) -> str:
    _element, product = _cge_mge_subject_parts(value)
    if not product:
        return ""
    marker_symbol = _canonical_marker_symbol(product)
    if marker_symbol:
        return marker_symbol
    return product if _is_useful_marker_shape_label(product) else ""


def _iter_cge_mge_subject_chunks(value: Any) -> Sequence[str]:
    text = _clean_context_text(value)
    if not text:
        return []
    chunks = [chunk.strip() for chunk in re.split(r"[,;]\s*", text) if chunk.strip()]
    subjects: List[str] = []
    for idx, chunk in enumerate(chunks):
        element, _product = _cge_mge_subject_parts(chunk)
        if not element:
            continue
        lookahead = " ".join(chunks[idx + 1 : idx + 4])
        if _is_cge_mge_source(lookahead):
            subjects.append(chunk)
    return subjects


def _cge_mge_product_label_from_hit(hit: Dict[str, Any]) -> str:
    if not isinstance(hit, dict):
        return ""
    source_blob = " ".join(
        _clean_context_text(hit.get(col))
        for col in ("DB", "Method")
        if _clean_context_text(hit.get(col))
    )
    if not _is_cge_mge_source(source_blob):
        return ""
    for col in ("SubjectIdentifier", "Product", "SubjectModelName"):
        label = _cge_mge_product_label(hit.get(col))
        if label:
            return label
    return ""


def _iter_cge_mge_subject_values(payload: Dict[str, Any]) -> Sequence[str]:
    values: List[str] = []
    for hit in _marker_hit_records(payload):
        if not isinstance(hit, dict):
            continue
        source_blob = " ".join(
            _clean_context_text(hit.get(col))
            for col in ("DB", "Method")
            if _clean_context_text(hit.get(col))
        )
        if not _is_cge_mge_source(source_blob):
            continue
        for col in ("SubjectIdentifier", "Product", "SubjectModelName"):
            value = _clean_context_text(hit.get(col))
            if value:
                values.append(value)

    source_blob = " ".join(
        _clean_context_text(payload.get(col))
        for col in ("ConsensusSource", "MarkerDisplay", "MarkerSearchText")
        if _clean_context_text(payload.get(col))
    )
    if _is_cge_mge_source(source_blob):
        values.extend(_iter_cge_mge_subject_chunks(payload.get("MarkerSearchText")))
        donor = _clean_context_text(payload.get("ConsensusDonorProduct"))
        if donor and _cge_mge_subject_parts(donor)[0]:
            values.append(donor)
    return values


def _cge_mge_product_label_for_gene(payload: Optional[Dict[str, Any]]) -> str:
    if not isinstance(payload, dict):
        return ""
    for value in _iter_cge_mge_subject_values(payload):
        label = _cge_mge_product_label(value)
        if label:
            return label
    return ""


def _is_cge_ice_marker_gene(payload: Optional[Dict[str, Any]]) -> bool:
    if not isinstance(payload, dict):
        return False
    for value in _iter_cge_mge_subject_values(payload):
        element, _product = _cge_mge_subject_parts(value)
        if _compact_marker_text(element).startswith("ice"):
            return True
    return False


def _figure_marker_text_blob(label: str, payload: Optional[Dict[str, Any]]) -> str:
    values: List[Any] = [label]
    if isinstance(payload, dict):
        values.extend(payload.get(col) for col in _MARKER_LABEL_COLUMNS)
        values.extend(
            payload.get(col)
            for col in (
                "Predictors",
                "MarkerDisplay",
                "ConsensusSource",
                "ResistanceDrugClasses",
                "ResistanceAntibiotics",
                "ResistanceMechanisms",
            )
        )
        for hit in _marker_hit_records(payload):
            values.extend(
                hit.get(col)
                for col in (
                    "Product",
                    "SubjectIdentifier",
                    "SubjectModelName",
                    "DB",
                    "Method",
                )
            )
    return " ".join(_clean_context_text(value) for value in values if _clean_context_text(value))


def _figure_shape_family_for_marker(
    display_label: str,
    payload: Optional[Dict[str, Any]],
    category: str,
) -> str:
    blob = _normalized_homology_text(_figure_marker_text_blob(display_label, payload))
    label_blob = _normalized_homology_text(display_label)
    label_words = label_blob.split()
    compact_label = _compact_marker_text(display_label)
    marker_symbol = _canonical_marker_symbol(display_label)
    marker_symbol_compact = _compact_marker_text(marker_symbol)

    def label_starts_with(*prefixes: str) -> bool:
        return compact_label.startswith(prefixes) or any(
            word.startswith(prefixes) for word in label_words
        )

    def label_matches(pattern: str) -> bool:
        return bool(re.search(pattern, label_blob))

    # Exact label-derived families win over broad annotation text. This keeps a
    # merT-style resistance marker from becoming a transposase because one of
    # its annotations also mentions an IS/Tn/mobile element.
    if compact_label.startswith("inti") or marker_symbol_compact.startswith("inti"):
        return "integrase"
    if (
        marker_symbol_compact.startswith("is")
        or re.match(r"^(?:tnp|tn\d)", compact_label)
    ):
        return "transposase"
    if label_matches(r"\b(?:recombinase|resolvase|invertase|xer[cd])\b"):
        return "recombinase"
    if label_starts_with("rep") or label_matches(r"\b(?:replication|replicase)\b"):
        return "replication"
    if label_starts_with("mob") or label_matches(
        r"\b(?:relaxase|orit|mobilization|mobilisation)\b"
    ):
        return "mobilization"
    if re.match(r"^(?:tra|trb)[a-z0-9]*", compact_label) or label_matches(
        r"\b(?:conjug|type iv secretion|pilus assembly|transfer protein)\b"
    ):
        return "conjugation"

    if label_starts_with("bla") or label_matches(
        r"\b(?:beta lactamase|carbapenemase|oxacillinase)\b"
    ):
        return "beta_lactamase"
    if (
        label_starts_with("aad", "aac", "aph", "ant", "arm", "rmt", "str")
        or "aminoglycoside" in label_blob
    ):
        return "aminoglycoside"
    if label_starts_with("dfr", "sul") or label_matches(
        r"\b(?:trimethoprim|sulfonamide|folate)\b"
    ):
        return "folate"
    if label_starts_with("erm", "mph", "msr", "mef", "lnu", "lsa", "cfr") or label_matches(
        r"\b(?:macrolide|lincosamide|streptogramin|ribosomal methyltransferase)\b"
    ):
        return "macrolide"
    if label_starts_with("tet") or "tetracycline" in label_blob:
        return "tetracycline"
    if (
        label_starts_with("qnr", "gyr", "par")
        or "quinolone" in label_blob
        or "fluoroquinolone" in label_blob
    ):
        return "quinolone"
    if label_starts_with("mcr", "pmr", "ept") or label_matches(
        r"\b(?:colistin|polymyxin|lipid a)\b"
    ):
        return "colistin"
    if label_starts_with(
        "qac", "sug", "emr", "mdf", "acr", "mex", "ade", "smr"
    ) or label_matches(r"\b(?:efflux|biocide|disinfectant|multidrug transporter)\b"):
        return "efflux"
    if label_starts_with(
        "ars", "mer", "czc", "cop", "cue", "cus", "sil", "znt", "cad", "chr", "nre", "nik"
    ) or label_matches(
        r"\b(?:metal|copper|arsenic|mercury|mercuric|zinc|cadmium|chromate|nickel)\b"
    ):
        return "metal_resistance"

    if label_matches(r"\b(?:toxin|hemolysin|haemolysin|cytotoxin|enterotoxin)\b"):
        return "toxin"
    if label_starts_with("fim", "afa", "pap", "sfa", "pil") or label_matches(
        r"\b(?:adhesin|fimbrial|fimbria|pilus|attachment)\b"
    ):
        return "adhesion"
    if label_matches(r"\b(?:secretion system|type iii|type vi|type ii|secreted effector)\b"):
        return "secretion"
    if label_starts_with("kps", "wzi", "wzc", "cps") or label_matches(
        r"\b(?:capsule|capsular|surface polysaccharide|lps)\b"
    ):
        return "capsule"
    if label_matches(
        r"\b(?:siderophore|yersiniabactin|aerobactin|enterobactin|iron uptake)\b"
    ):
        return "siderophore"

    if category == "MOBILE":
        if re.search(r"\b(?:integrase|integron)\b", blob):
            return "integrase"
        if re.search(r"\b(?:recombinase|resolvase|invertase|xer[cd])\b", blob):
            return "recombinase"
        if re.search(r"\b(?:replication|replicase)\b", blob):
            return "replication"
        if re.search(r"\b(?:relaxase|orit|mobilization|mobilisation)\b", blob):
            return "mobilization"
        if re.search(r"\b(?:conjug|type iv secretion|pilus assembly|transfer protein)\b", blob):
            return "conjugation"
        if re.search(r"\b(?:transposase|insertion sequence|transposon)\b", blob):
            return "transposase"

    if re.search(r"\b(?:beta lactamase|carbapenemase|oxacillinase)\b", blob):
        return "beta_lactamase"
    if "aminoglycoside" in blob:
        return "aminoglycoside"
    if re.search(r"\b(?:trimethoprim|sulfonamide|folate)\b", blob):
        return "folate"
    if re.search(r"\b(?:macrolide|lincosamide|streptogramin|ribosomal methyltransferase)\b", blob):
        return "macrolide"
    if "tetracycline" in blob:
        return "tetracycline"
    if "quinolone" in blob or "fluoroquinolone" in blob:
        return "quinolone"
    if re.search(r"\b(?:colistin|polymyxin|lipid a)\b", blob):
        return "colistin"
    if re.search(r"\b(?:efflux|biocide|disinfectant|multidrug transporter)\b", blob):
        return "efflux"
    if re.search(r"\b(?:metal|copper|arsenic|mercury|mercuric|zinc|cadmium|chromate|nickel)\b", blob):
        return "metal_resistance"
    if re.search(r"\b(?:toxin|hemolysin|haemolysin|cytotoxin|enterotoxin)\b", blob):
        return "toxin"
    if re.search(r"\b(?:adhesin|fimbrial|fimbria|pilus|attachment)\b", blob):
        return "adhesion"
    if re.search(r"\b(?:secretion system|type iii|type vi|type ii|secreted effector)\b", blob):
        return "secretion"
    if re.search(r"\b(?:capsule|capsular|surface polysaccharide|lps)\b", blob):
        return "capsule"
    if re.search(r"\b(?:siderophore|yersiniabactin|aerobactin|enterobactin|iron uptake)\b", blob):
        return "siderophore"

    if category == "BMRG":
        return "efflux" if "biocide" in blob or "efflux" in blob else "metal_resistance"
    if category == "MOBILE":
        return "transposase"
    return "marker_other"


def _homology_key_for_gene(gene: Dict[str, Any]) -> str:
    payload = gene.get("row_data") if isinstance(gene.get("row_data"), dict) else {}

    for col in (
        "PanarooCluster",
        "PanarooGeneCluster",
        "PanarooGene",
        "PanarooID",
        "GeneCluster",
        "Cluster",
    ):
        value = _clean_context_text(payload.get(col))
        if _is_informative_homology_text(value):
            return f"cluster:{_normalized_homology_text(value)}"

    for col in ("Gene", "GeneID"):
        value = _clean_context_text(payload.get(col))
        if value.startswith("UNK:"):
            continue
        marker_symbol = _canonical_marker_symbol(value)
        if marker_symbol:
            return f"gene:{_normalized_homology_text(marker_symbol)}"
        if _is_informative_homology_text(value):
            return f"gene:{_normalized_homology_text(value)}"

    for col in ("MarkerSearchText", "ConsensusDonorProduct", "Product"):
        marker_symbol = _canonical_marker_symbol(payload.get(col))
        if marker_symbol:
            return f"gene:{_normalized_homology_text(marker_symbol)}"

    marker = _first_informative_marker_token(payload.get("MarkerSearchText"))
    if marker:
        return f"marker:{_normalized_homology_text(marker)}"

    for col in ("ConsensusDonorProduct", "Product"):
        value = _clean_context_text(payload.get(col))
        if _is_informative_homology_text(value):
            return f"product:{_normalized_homology_text(value)}"

    aa_seq = _clean_context_text(payload.get("Amino_Acid_Sequence")).replace("*", "")
    if len(aa_seq) >= 20:
        digest = hashlib.sha1(aa_seq.encode("utf-8")).hexdigest()[:16]
        return f"aa:{digest}"

    label = _clean_context_text(gene.get("label"))
    marker_symbol = _canonical_marker_symbol(label)
    if marker_symbol:
        return f"gene:{_normalized_homology_text(marker_symbol)}"
    if _is_informative_homology_text(label):
        return f"label:{_normalized_homology_text(label)}"
    return ""


def _homology_evidence_for_key(key: str) -> Tuple[str, float]:
    prefix = str(key or "").split(":", 1)[0].strip().lower()
    if prefix == "cluster":
        return "Panaroo/cluster", 1.30
    if prefix == "aa":
        return "protein sequence", 1.18
    if prefix == "gene":
        return "gene", 1.00
    if prefix == "marker":
        return "database annotation", 0.92
    if prefix == "product":
        return "product", 0.72
    if prefix == "label":
        return "label", 0.58
    return "unknown", 0.45


def _homology_display_for_key(key: str) -> str:
    text = str(key or "").split(":", 1)[-1].strip()
    return text or "unknown"


def _gene_homology_profile(gene: Dict[str, Any]) -> Dict[str, Any]:
    key = _homology_key_for_gene(gene)
    evidence, weight = _homology_evidence_for_key(key)
    return {
        "key": key,
        "label": _homology_display_for_key(key),
        "evidence": evidence,
        "weight": float(weight),
    }


_MARKER_LABEL_COLUMNS = (
    "ResistanceAlleles",
    "Gene",
    "GeneID",
    "MarkerSearchText",
    "ConsensusDonorProduct",
    "Product",
    "ResistanceGeneFamilies",
)
_MARKER_CATEGORY_COLUMNS = ("ConsensusCategory", "MarkerCategory", "Category")
_MARKER_HIT_COLLECTION_COLUMNS = ("MarkerHits", "Markers")
_RESISTANCE_EVIDENCE_COLUMNS = (
    "ResistanceMetadataSource",
    "ResistanceAlleles",
    "ResistanceGeneFamilies",
    "ResistanceDrugClasses",
    "ResistanceAntibiotics",
    "ResistanceMechanisms",
)


def _marker_hit_records(payload: Dict[str, Any]) -> List[Dict[str, Any]]:
    records: List[Dict[str, Any]] = []
    for col in _MARKER_HIT_COLLECTION_COLUMNS:
        value = payload.get(col)
        if not isinstance(value, list):
            continue
        records.extend(item for item in value if isinstance(item, dict))
    return records


def _figure_marker_category_for_gene(payload: Dict[str, Any]) -> str:
    category_values: List[Any] = []
    for col in _MARKER_CATEGORY_COLUMNS:
        value = payload.get(col)
        if _clean_context_text(value):
            category_values.append(value)

    for hit in _marker_hit_records(payload):
        for col in _MARKER_CATEGORY_COLUMNS:
            value = hit.get(col)
            if _clean_context_text(value):
                category_values.append(value)

    category = _best_figure_category(category_values)
    if category:
        return category
    if category_values:
        return ""

    source_values = [
        payload.get(col)
        for col in ("Predictors", "MarkerDisplay", "ConsensusSource", "MarkerSearchText")
    ]
    for hit in _marker_hit_records(payload):
        source_values.extend(hit.get(col) for col in ("DB", "Method", "Product"))
    source_blob = _normalized_homology_text(",".join(_clean_context_text(value) for value in source_values))

    if "bmrg" in source_blob:
        return "BMRG"
    if "vfdb" in source_blob or "virulence" in source_blob:
        return "VIRULENCE"
    if "isfinder" in source_blob or "mge" in source_blob or "mobile" in source_blob:
        return "MOBILE"
    if "card" in source_blob or "amrf" in source_blob or re.search(r"\bamr\b", source_blob):
        return "AMR"

    # Older caches may only carry resistance metadata summaries, which are still
    # the ARG classes we want in figure mode.
    if any(_clean_context_text(payload.get(col)) for col in _RESISTANCE_EVIDENCE_COLUMNS):
        return "AMR"
    return ""


def _has_figure_marker_category(payload: Dict[str, Any]) -> bool:
    return bool(_figure_marker_category_for_gene(payload))


def _marker_annotation_label_for_gene(label: str, payload: Optional[Dict[str, Any]]) -> str:
    if not isinstance(payload, dict):
        return ""
    if not _has_figure_marker_category(payload):
        return ""
    marker_hit_list = _marker_hit_records(payload)
    category = _figure_marker_category_for_gene(payload)

    for col in _MARKER_LABEL_COLUMNS:
        if category == "MOBILE" and col == "MarkerSearchText":
            label_from_cge = _cge_mge_product_label_for_gene(payload)
            if label_from_cge:
                return label_from_cge
            continue
        value = payload.get(col)
        marker_symbol = _canonical_marker_symbol(value)
        if marker_symbol:
            return marker_symbol

    for hit in marker_hit_list:
        if not isinstance(hit, dict):
            continue
        if category == "MOBILE":
            label_from_cge = _cge_mge_product_label_from_hit(hit)
            if label_from_cge:
                return label_from_cge
        for col in ("Product", "SubjectIdentifier", "SubjectModelName", "Method", "DB"):
            marker_symbol = _canonical_marker_symbol(hit.get(col))
            if marker_symbol:
                return marker_symbol

    if category == "MOBILE":
        label_from_cge = _cge_mge_product_label_for_gene(payload)
        if label_from_cge:
            return label_from_cge

    for col in _MARKER_LABEL_COLUMNS:
        if category == "MOBILE" and col == "MarkerSearchText":
            label_from_cge = _cge_mge_product_label_for_gene(payload)
            if label_from_cge:
                return label_from_cge
            continue
        token = _first_useful_marker_shape_token(payload.get(col))
        if token:
            return token

    for hit in marker_hit_list:
        if not isinstance(hit, dict):
            continue
        for col in ("Product", "SubjectIdentifier", "SubjectModelName"):
            token = _first_useful_marker_shape_token(hit.get(col))
            if token:
                return token

    marker_symbol = _canonical_marker_symbol(label)
    if marker_symbol:
        return marker_symbol

    label_token = _first_useful_marker_shape_token(label)
    if label_token:
        return label_token

    return ""


class IgvCollectionView(QtWidgets.QGraphicsView):
    contextSamplesHovered = QtCore.pyqtSignal(object)
    geneDoubleClicked = QtCore.pyqtSignal(dict)
    syntenyBlockClicked = QtCore.pyqtSignal(dict)
    viewportChanged = QtCore.pyqtSignal()

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self._scene = QtWidgets.QGraphicsScene(self)
        self.setScene(self._scene)
        self.setBackgroundBrush(QtGui.QBrush(QtGui.QColor("#f8fafc")))
        self.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        self.setDragMode(QtWidgets.QGraphicsView.DragMode.ScrollHandDrag)
        self.setTransformationAnchor(QtWidgets.QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setResizeAnchor(QtWidgets.QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setViewportUpdateMode(QtWidgets.QGraphicsView.ViewportUpdateMode.BoundingRectViewportUpdate)
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.setMinimumHeight(320)
        self.setMouseTracking(True)
        self.viewport().setMouseTracking(True)
        self._rows: List[Dict[str, Any]] = []
        self._min_rel_bp = -1000.0
        self._max_rel_bp = 1000.0
        self._row_height = 22.0
        self._row_gap = 150.0
        self._top_margin = 86.0
        self._bottom_margin = 28.0
        self._label_width = 40.0
        self._arrow_bp = 90.0
        self._gene_label_gap_px = 8.0
        self._row_label_left_padding_px = 12.0
        self._row_label_track_gap_px = 34.0
        self._row_label_min_lane_px = 300.0
        self._row_label_max_lane_px = 560.0
        self._min_zoom_x = 0.03
        self._figure_min_zoom_x = 0.045
        self._max_zoom_x = 20.0
        self._min_zoom_y = 0.16
        self._max_zoom_y = 1.0
        self._max_synteny_gap_genes = 3
        self._show_synteny_blocks = True
        self._left_align_contexts = True
        self._show_molecule_badges = True
        self._static_row_labels = True
        # Pixel lane before the first gene track. It is converted to scene units
        # on render so zooming cannot collapse labels into the track.
        self._track_x_offset = self._row_label_min_lane_px
        self._computed_track_x_offset = self._track_x_offset
        self._display_mode = "exploratory"
        self._row_label_items: List[Tuple[QtWidgets.QGraphicsSimpleTextItem, float, float]] = []
        self._gene_label_items: List[Tuple[QtWidgets.QGraphicsSimpleTextItem, str, float, float, float]] = []
        self._gene_visuals: List[Dict[str, Any]] = []
        self._synteny_band_items: List[QtWidgets.QGraphicsPathItem] = []
        self._hover_highlight_items: List[QtWidgets.QGraphicsItem] = []
        self._figure_styles: Dict[str, Dict[str, Any]] = {}
        self._figure_legend_entries: Dict[str, Dict[str, Any]] = {}
        self._figure_legend_items: List[Tuple[QtWidgets.QGraphicsItem, float, float]] = []
        self._figure_legend_size = QtCore.QSizeF(0.0, 0.0)
        self._selected_gene_key: Optional[str] = None
        self._hover_homology_key: Optional[str] = None
        self._selection_line_items: List[QtWidgets.QGraphicsLineItem] = []
        self._hovered_samples_key: Optional[Tuple[str, ...]] = None

        self.horizontalScrollBar().valueChanged.connect(self._on_scrollbar_changed)
        self.verticalScrollBar().valueChanged.connect(self._on_scrollbar_changed)

    def clear(self) -> None:
        self._emit_context_hover(None)
        self._rows = []
        self._row_label_items = []
        self._gene_label_items = []
        self._gene_visuals = []
        self._synteny_band_items = []
        self._hover_highlight_items = []
        self._figure_styles = {}
        self._figure_legend_entries = {}
        self._figure_legend_items = []
        self._selected_gene_key = None
        self._hover_homology_key = None
        self._selection_line_items = []
        self._scene.clear()
        self.resetTransform()
        self.viewportChanged.emit()

    def set_synteny_blocks_enabled(self, enabled: bool) -> None:
        self._show_synteny_blocks = bool(enabled)
        self._render_scene()

    def set_molecule_badges_visible(self, visible: bool) -> None:
        next_value = bool(visible)
        if next_value == self._show_molecule_badges:
            return
        self._show_molecule_badges = next_value
        self._render_scene()

    def set_static_row_labels(self, enabled: bool, *, track_offset: float = 0.0) -> None:
        next_enabled = bool(enabled)
        next_offset = max(0.0, float(track_offset or 0.0))
        if next_enabled and next_offset <= 0.0:
            next_offset = self._row_label_min_lane_px
        if next_enabled == self._static_row_labels and abs(next_offset - self._track_x_offset) < 1e-6:
            return
        self._static_row_labels = next_enabled
        self._track_x_offset = next_offset if next_enabled else 0.0
        self._render_scene()

    def _scale_x(self) -> float:
        return max(abs(float(self.transform().m11())), 1e-6)

    def _scale_y(self) -> float:
        return max(abs(float(self.transform().m22())), 1e-6)

    def _static_track_offset_scene(self) -> float:
        if not self._static_row_labels:
            return 0.0
        return float(self._computed_track_x_offset) / self._scale_x()

    def _row_label_text_for_row(self, row: Dict[str, Any], row_idx: int) -> Tuple[str, str]:
        sample = str(row.get("sample") or "").strip()
        contig = str(row.get("contig") or "").strip()
        display_label = str(row.get("display_label") or "").strip()
        title = display_label or str(row.get("title") or f"Context {row_idx + 1}")
        if not display_label and sample and contig:
            title = f"{sample} | {contig}"
        badge_text, _badge_color = _molecule_badge(row)
        return str(title), str(badge_text or "")

    def _compute_track_lane_px(self) -> float:
        if not self._static_row_labels:
            return 0.0
        title_font = QtGui.QFont()
        title_font.setPointSize(8)
        title_font.setBold(True)
        badge_font = QtGui.QFont(title_font)
        title_metrics = QtGui.QFontMetrics(title_font)
        badge_metrics = QtGui.QFontMetrics(badge_font)
        widest_label_px = 0.0
        for row_idx, row in enumerate(self._rows):
            if not isinstance(row, dict):
                continue
            title, badge_text = self._row_label_text_for_row(row, row_idx)
            width = float(title_metrics.horizontalAdvance(title))
            if badge_text and self._show_molecule_badges:
                width += float(badge_metrics.horizontalAdvance(badge_text) + 8)
            widest_label_px = max(widest_label_px, width)
        lane_px = self._row_label_left_padding_px + widest_label_px + self._row_label_track_gap_px
        lane_px = max(float(self._track_x_offset), lane_px, self._row_label_min_lane_px)
        return min(lane_px, self._row_label_max_lane_px)

    def set_display_mode(self, mode: str) -> None:
        normalized = str(mode or "").strip().lower()
        if normalized not in {"exploratory", "figure"}:
            normalized = "exploratory"
        if normalized == self._display_mode:
            return
        self._display_mode = normalized
        self._hover_homology_key = None
        self._render_scene()

    def set_rows(self, rows: Sequence[Dict[str, Any]]) -> None:
        self._emit_context_hover(None)
        self._rows = [dict(r) for r in rows if isinstance(r, dict)]
        self.resetTransform()
        self._render_scene()
        QtCore.QTimer.singleShot(0, self._center_on_origin)

    def _center_on_origin(self) -> None:
        if self._scene is None or self._scene.itemsBoundingRect().isNull():
            return
        if self._left_align_contexts:
            self.horizontalScrollBar().setValue(self.horizontalScrollBar().minimum())
            self._reposition_row_labels()
            return
        self.centerOn(self._scene.sceneRect().center())

    def _compute_row_anchors(self) -> Tuple[List[int], float, float]:
        row_anchors: List[int] = []
        min_rel = None
        max_rel = None
        for row in self._rows:
            genes = row.get("genes") or []
            if not isinstance(genes, list):
                row_anchors.append(0)
                continue
            spans: List[Tuple[int, int]] = []
            for gene in genes:
                if not isinstance(gene, dict):
                    continue
                try:
                    start = int(gene.get("start") or 0)
                    end = int(gene.get("end") or 0)
                except Exception:
                    continue
                spans.append((min(start, end), max(start, end)))
            if not spans:
                row_anchors.append(0)
                continue
            row_anchor = min(start for start, _end in spans)
            row_anchors.append(row_anchor)
            for start, end in spans:
                rel_start = float(start - row_anchor)
                rel_end = float(end - row_anchor)
                min_rel = rel_start if min_rel is None else min(min_rel, rel_start)
                max_rel = rel_end if max_rel is None else max(max_rel, rel_end)

        if min_rel is None or max_rel is None:
            return row_anchors, -1000.0, 1000.0
        if max_rel <= min_rel:
            max_rel = min_rel + 1.0
        return row_anchors, float(min_rel), float(max_rel)

    def _render_scene(self) -> None:
        self._scene.clear()
        self._row_label_items = []
        self._gene_label_items = []
        self._gene_visuals = []
        self._synteny_band_items = []
        self._hover_highlight_items = []
        self._figure_styles = {}
        self._figure_legend_entries = {}
        self._figure_legend_items = []
        self._figure_legend_size = QtCore.QSizeF(0.0, 0.0)
        self._selection_line_items = []
        if not self._rows:
            return
        selected_before = self._selected_gene_key

        row_anchors, min_rel, max_rel = self._compute_row_anchors()
        self._min_rel_bp = min_rel
        self._max_rel_bp = max_rel
        self._computed_track_x_offset = self._compute_track_lane_px()
        track_offset = self._static_track_offset_scene()
        x0 = track_offset + self._min_rel_bp - 18.0
        x1 = track_offset + self._max_rel_bp + 36.0

        for row_idx, row in enumerate(self._rows):
            y_top = self._top_margin + row_idx * (self._row_height + self._row_gap)
            y_center = y_top + self._row_height / 2.0
            title, _row_badge_text = self._row_label_text_for_row(row, row_idx)
            tooltip = str(row.get("meta_tooltip") or "").strip() or title

            baseline = QtWidgets.QGraphicsLineItem(track_offset, y_center, x1, y_center)
            baseline.setPen(QtGui.QPen(QtGui.QColor("#d8dee7"), 0.9))
            baseline.setToolTip(tooltip)
            baseline.setData(1, int(row_idx))
            baseline.setData(2, "track")
            baseline.setZValue(2)
            self._scene.addItem(baseline)
            badge_text, badge_color = _molecule_badge(row)
            title_x_offset = 0.0
            if badge_text and self._show_molecule_badges:
                badge_item = QtWidgets.QGraphicsSimpleTextItem(badge_text)
                badge_font = QtGui.QFont()
                badge_font.setPointSize(8)
                badge_font.setBold(True)
                badge_item.setFont(badge_font)
                badge_item.setBrush(QtGui.QBrush(QtGui.QColor(badge_color)))
                badge_item.setToolTip(tooltip)
                badge_item.setData(1, int(row_idx))
                badge_item.setData(2, "label")
                badge_item.setData(4, "badge")
                badge_item.setFlag(
                    QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations, True
                )
                badge_item.setZValue(41)
                self._scene.addItem(badge_item)
                self._row_label_items.append((badge_item, y_top, 0.0))
                title_x_offset = float(QtGui.QFontMetrics(badge_font).horizontalAdvance(badge_text) + 8)

            title_item = QtWidgets.QGraphicsSimpleTextItem(title)
            title_font = QtGui.QFont()
            title_font.setPointSize(8)
            title_font.setBold(True)
            title_item.setFont(title_font)
            title_item.setBrush(QtGui.QBrush(QtGui.QColor("#4b5563")))
            title_item.setToolTip(tooltip)
            title_item.setData(1, int(row_idx))
            title_item.setData(2, "label")
            title_item.setData(4, "title")
            title_item.setData(5, title)
            title_item.setFlag(
                QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations, True
            )
            title_item.setZValue(40)
            self._scene.addItem(title_item)
            self._row_label_items.append((title_item, y_top, title_x_offset))
            if self._static_row_labels:
                scale_x = self._scale_x()
                left_scene = self.mapToScene(QtCore.QPoint(0, 0)).x()
                title_item.setPos(
                    left_scene + ((self._row_label_left_padding_px + title_x_offset) / scale_x),
                    float(y_top - 12.0),
                )

            genes = row.get("genes") or []
            if not isinstance(genes, list):
                continue
            row_anchor = row_anchors[row_idx] if row_idx < len(row_anchors) else 0
            for gene_idx, gene in enumerate(genes):
                if not isinstance(gene, dict):
                    continue
                try:
                    start = int(gene.get("start") or 0)
                    end = int(gene.get("end") or 0)
                except Exception:
                    continue
                strand = str(gene.get("strand") or "+").strip() or "+"
                label = str(gene.get("label") or "gene")
                color_value = gene.get("color")
                color = color_value if isinstance(color_value, QtGui.QColor) else QtGui.QColor(str(color_value or "#90caf9"))
                rel_start = track_offset + float(min(start, end) - row_anchor)
                rel_end = track_offset + float(max(start, end) - row_anchor)
                payload = gene.get("row_data") if isinstance(gene.get("row_data"), dict) else None
                is_center = bool(gene.get("is_center"))
                self._add_gene_item(
                    rel_start,
                    rel_end,
                    y_top,
                    row_idx=row_idx,
                    gene_idx=gene_idx,
                    strand=strand,
                    label=label,
                    color=color,
                    payload=payload,
                    is_center=is_center,
                )

        self._render_synteny_layer()

        total_height = self._top_margin + len(self._rows) * (self._row_height + self._row_gap) + self._bottom_margin
        self._scene.setSceneRect(
            0.0 if self._static_row_labels else x0 - self._label_width,
            0.0,
            (x1 + 24.0) if self._static_row_labels else (x1 - x0) + self._label_width + 24.0,
            total_height,
        )
        self._render_figure_legend()
        self._update_zoom_dependent_items()
        self._apply_gene_selection(selected_before)
        self.viewportChanged.emit()

    def _figure_style_for_gene(
        self,
        homology_key: str,
        label: str,
        payload: Optional[Dict[str, Any]],
    ) -> Dict[str, Any]:
        category = _figure_marker_category_for_gene(payload) if isinstance(payload, dict) else ""
        marker_label = _marker_annotation_label_for_gene(label, payload)
        if category == "MOBILE" and not marker_label:
            marker_label = "mobile element"
        if not category or not marker_label:
            return {
                "shape": "square",
                "color": QtGui.QColor("#cbd5e1"),
                "label": "",
                "legend": False,
            }
        display_label = marker_label
        family = _figure_shape_family_for_marker(display_label, payload, category)
        shape = "triangle" if category == "MOBILE" else _FIGURE_FAMILY_SHAPES.get(family, "square")
        ice_associated = bool(category == "MOBILE" and _is_cge_ice_marker_gene(payload))
        ice_key = "ice" if ice_associated else "plain"
        style_key = f"{category}:{family}:{_normalized_homology_text(display_label)}:{ice_key}"
        if style_key in self._figure_styles:
            self._figure_styles[style_key]["count"] = int(self._figure_styles[style_key].get("count") or 0) + 1
            return self._figure_styles[style_key]

        style = {
            "shape": shape,
            "color": _stable_figure_category_color(category, display_label),
            "label": str(display_label),
            "category": category,
            "family": family,
            "family_label": _FIGURE_FAMILY_LABELS.get(family, family.replace("_", " ")),
            "ice": ice_associated,
            "legend": True,
            "count": 1,
        }
        self._figure_styles[style_key] = style
        self._figure_legend_entries[style_key] = style
        return style

    @staticmethod
    def _figure_symbol_path(shape: str, rect: QtCore.QRectF) -> QtGui.QPainterPath:
        path = QtGui.QPainterPath()
        cx = rect.center().x()
        cy = rect.center().y()
        left = rect.left()
        right = rect.right()
        top = rect.top()
        bottom = rect.bottom()
        if shape == "circle":
            path.addEllipse(rect)
        elif shape == "square":
            path.addRoundedRect(rect, 2.0, 2.0)
        elif shape == "diamond":
            path.addPolygon(
                QtGui.QPolygonF(
                    [
                        QtCore.QPointF(cx, top),
                        QtCore.QPointF(right, cy),
                        QtCore.QPointF(cx, bottom),
                        QtCore.QPointF(left, cy),
                    ]
                )
            )
        elif shape == "triangle":
            path.addPolygon(
                QtGui.QPolygonF(
                    [
                        QtCore.QPointF(cx, top),
                        QtCore.QPointF(right, bottom),
                        QtCore.QPointF(left, bottom),
                    ]
                )
            )
        elif shape == "down_triangle":
            path.addPolygon(
                QtGui.QPolygonF(
                    [
                        QtCore.QPointF(left, top),
                        QtCore.QPointF(right, top),
                        QtCore.QPointF(cx, bottom),
                    ]
                )
            )
        elif shape == "hexagon":
            path.addPolygon(
                QtGui.QPolygonF(
                    [
                        QtCore.QPointF(left + rect.width() * 0.25, top),
                        QtCore.QPointF(right - rect.width() * 0.25, top),
                        QtCore.QPointF(right, cy),
                        QtCore.QPointF(right - rect.width() * 0.25, bottom),
                        QtCore.QPointF(left + rect.width() * 0.25, bottom),
                        QtCore.QPointF(left, cy),
                    ]
                )
            )
        elif shape == "pentagon":
            path.addPolygon(
                QtGui.QPolygonF(
                    [
                        QtCore.QPointF(cx, top),
                        QtCore.QPointF(right, top + rect.height() * 0.38),
                        QtCore.QPointF(right - rect.width() * 0.18, bottom),
                        QtCore.QPointF(left + rect.width() * 0.18, bottom),
                        QtCore.QPointF(left, top + rect.height() * 0.38),
                    ]
                )
            )
        else:
            bar_w = max(2.0, rect.width() * 0.32)
            bar_h = max(2.0, rect.height() * 0.32)
            path.addRoundedRect(QtCore.QRectF(cx - bar_w / 2.0, top, bar_w, rect.height()), 1.5, 1.5)
            path.addRoundedRect(QtCore.QRectF(left, cy - bar_h / 2.0, rect.width(), bar_h), 1.5, 1.5)
        return path

    def _add_figure_gene_item(
        self,
        *,
        x0: float,
        x1: float,
        y_top: float,
        row_idx: int,
        gene_idx: int,
        strand: str,
        label: str,
        payload: Optional[Dict[str, Any]],
        is_center: bool,
        gene_key: str,
        homology_key: str,
        homology_profile: Dict[str, Any],
    ) -> None:
        total_w = max(1.0, x1 - x0)
        style = self._figure_style_for_gene(homology_key, label, payload)
        cx = (x0 + x1) / 2.0
        cy = y_top + self._row_height / 2.0
        symbol_size = 28.0 if bool(style.get("legend")) else 20.0
        symbol_rect = QtCore.QRectF(-symbol_size / 2.0, -symbol_size / 2.0, symbol_size, symbol_size)
        item = QtWidgets.QGraphicsPathItem(self._figure_symbol_path(str(style.get("shape") or "circle"), symbol_rect))
        item.setBrush(QtGui.QBrush(style.get("color") if isinstance(style.get("color"), QtGui.QColor) else QtGui.QColor("#cbd5e1")))
        pen_color = QtGui.QColor("#f97316") if is_center else QtGui.QColor("#0f172a")
        item.setPen(QtGui.QPen(pen_color, 1.5 if is_center else 0.9))
        item.setFlag(QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations, True)
        item.setPos(cx, cy)
        item.setToolTip(label)
        if payload is not None:
            item.setData(0, payload)
        item.setData(1, int(row_idx))
        item.setData(2, "track")
        item.setData(3, homology_key)
        item.setZValue(24)
        self._scene.addItem(item)
        overlay_items: List[QtWidgets.QGraphicsItem] = []
        if bool(style.get("ice")):
            ice_font = QtGui.QFont()
            ice_font.setPointSize(6)
            ice_font.setBold(True)
            ice_badge = QtWidgets.QGraphicsSimpleTextItem("ICE")
            ice_badge.setFont(ice_font)
            ice_badge.setBrush(QtGui.QBrush(QtGui.QColor("#713f12")))
            ice_badge.setFlag(QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations, True)
            ice_badge.setToolTip(f"{label}\nCGE MGE hit belongs to an ICE reference element")
            if payload is not None:
                ice_badge.setData(0, payload)
            ice_badge.setData(1, int(row_idx))
            ice_badge.setData(2, "track")
            ice_badge.setData(3, homology_key)
            ice_badge.setZValue(35)
            metrics = QtGui.QFontMetrics(ice_font)
            ice_badge.setPos(
                cx - (float(metrics.horizontalAdvance("ICE")) / 2.0),
                cy - (float(metrics.height()) / 2.0) + 1.0,
            )
            self._scene.addItem(ice_badge)
            overlay_items.append(ice_badge)
        ribbon_half_width = max(12.0, min(45.0, total_w * 0.18))
        self._gene_visuals.append(
            {
                "gene_key": gene_key,
                "homology_key": homology_key,
                "homology_label": str(homology_profile.get("label") or ""),
                "homology_evidence": str(homology_profile.get("evidence") or ""),
                "homology_weight": float(homology_profile.get("weight") or 0.0),
                "row_idx": int(row_idx),
                "gene_idx": int(gene_idx),
                "label": str(label),
                "strand": str(strand or "+"),
                "x0": float(cx - ribbon_half_width),
                "x1": float(cx + ribbon_half_width),
                "x_center": float(cx),
                "y_top": float(y_top),
                "shape_items": [item],
                "base_pens": [QtGui.QPen(item.pen())],
                "overlay_items": overlay_items,
                "text_item": None,
                "text_base_color": QtGui.QColor("#0f172a"),
                "fixed_symbol": True,
                "fixed_symbol_size": float(symbol_size),
            }
        )

    def _add_gene_item(
        self,
        x_start: float,
        x_end: float,
        y_top: float,
        *,
        row_idx: int,
        gene_idx: int,
        strand: str,
        label: str,
        color: QtGui.QColor,
        payload: Optional[Dict[str, Any]],
        is_center: bool,
    ) -> None:
        x0 = min(x_start, x_end)
        x1 = max(x_start, x_end)
        if x1 <= x0:
            return
        total_w = x1 - x0
        arrow_w = min(
            total_w,
            max(8.0, min(self._arrow_bp, total_w * 0.35)),
        )
        outline_color = QtGui.QColor("#111827") if is_center else QtGui.QColor("#374151")
        outline_w = 1.2 if is_center else 0.8
        shape_items: List[QtWidgets.QGraphicsItem] = []
        base_pens: List[QtGui.QPen] = []
        gene_key = self._gene_key(payload, label)
        homology_profile = _gene_homology_profile({"label": label, "row_data": payload or {}})
        homology_key = str(homology_profile.get("key") or "")
        if self._display_mode == "figure":
            self._add_figure_gene_item(
                x0=float(x0),
                x1=float(x1),
                y_top=float(y_top),
                row_idx=int(row_idx),
                gene_idx=int(gene_idx),
                strand=strand,
                label=label,
                payload=payload,
                is_center=is_center,
                gene_key=gene_key,
                homology_key=homology_key,
                homology_profile=homology_profile,
            )
            return

        body_rect: Optional[QtCore.QRectF] = None
        is_reverse = strand == "-"
        if is_reverse:
            head_base_x = x0 + arrow_w
            body_w = max(0.0, x1 - head_base_x)
            if body_w > 0:
                body_rect = QtCore.QRectF(head_base_x, y_top, body_w, self._row_height)
            arrow = QtGui.QPolygonF(
                [
                    QtCore.QPointF(head_base_x, y_top),
                    QtCore.QPointF(x0, y_top + self._row_height / 2.0),
                    QtCore.QPointF(head_base_x, y_top + self._row_height),
                ]
            )
        else:
            head_base_x = x1 - arrow_w
            body_w = max(0.0, head_base_x - x0)
            if body_w > 0:
                body_rect = QtCore.QRectF(x0, y_top, body_w, self._row_height)
            arrow = QtGui.QPolygonF(
                [
                    QtCore.QPointF(head_base_x, y_top),
                    QtCore.QPointF(x1, y_top + self._row_height / 2.0),
                    QtCore.QPointF(head_base_x, y_top + self._row_height),
                ]
            )

        parent_item: Optional[QtWidgets.QGraphicsRectItem] = None
        if body_rect is not None:
            parent_item = QtWidgets.QGraphicsRectItem(body_rect)
            parent_item.setBrush(QtGui.QBrush(color))
            parent_item.setPen(QtGui.QPen(outline_color, outline_w))
            parent_item.setToolTip(label)
            if payload is not None:
                parent_item.setData(0, payload)
            parent_item.setData(1, int(row_idx))
            parent_item.setData(2, "track")
            parent_item.setData(3, homology_key)
            parent_item.setZValue(24)
            self._scene.addItem(parent_item)
            shape_items.append(parent_item)
            base_pens.append(QtGui.QPen(parent_item.pen()))

        arrow_item = QtWidgets.QGraphicsPolygonItem(arrow)
        arrow_item.setBrush(QtGui.QBrush(color))
        arrow_item.setPen(QtGui.QPen(outline_color, outline_w))
        arrow_item.setToolTip(label)
        if payload is not None:
            arrow_item.setData(0, payload)
        arrow_item.setData(1, int(row_idx))
        arrow_item.setData(2, "track")
        arrow_item.setData(3, homology_key)
        arrow_item.setZValue(24)
        self._scene.addItem(arrow_item)
        shape_items.append(arrow_item)
        base_pens.append(QtGui.QPen(arrow_item.pen()))

        text_item: Optional[QtWidgets.QGraphicsSimpleTextItem] = None
        if label and total_w >= 8.0 and _is_informative_homology_text(label):
            font = QtGui.QFont()
            font.setPointSize(9)
            text_item = QtWidgets.QGraphicsSimpleTextItem(label)
            text_item.setFont(font)
            text_item.setBrush(QtGui.QBrush(QtGui.QColor("#1f2937")))
            text_item.setFlag(
                QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations, True
            )
            text_item.setRotation(-45.0)
            text_item.setData(1, int(row_idx))
            text_item.setData(2, "track")
            text_item.setData(3, homology_key)
            text_item.setZValue(30)
            self._position_gene_label_item(text_item, (x0 + x1) / 2.0, y_top)
            self._scene.addItem(text_item)
            self._gene_label_items.append(
                (
                    text_item,
                    str(label),
                    float(max(1.0, total_w)),
                    float((x0 + x1) / 2.0),
                    float(y_top),
                )
            )

        self._gene_visuals.append(
            {
                "gene_key": gene_key,
                "homology_key": homology_key,
                "homology_label": str(homology_profile.get("label") or ""),
                "homology_evidence": str(homology_profile.get("evidence") or ""),
                "homology_weight": float(homology_profile.get("weight") or 0.0),
                "row_idx": int(row_idx),
                "gene_idx": int(gene_idx),
                "label": str(label),
                "strand": str(strand or "+"),
                "x0": float(x0),
                "x1": float(x1),
                "x_center": float((x0 + x1) / 2.0),
                "y_top": float(y_top),
                "shape_items": shape_items,
                "base_pens": base_pens,
                "text_item": text_item,
                "text_base_color": QtGui.QColor("#0f172a"),
            }
        )

    @staticmethod
    def _normalized_strand(value: Any) -> str:
        strand = str(value or "+").strip()
        return "-" if strand == "-" else "+"

    def _row_visuals_by_index(self) -> Dict[int, List[Dict[str, Any]]]:
        by_row: Dict[int, List[Dict[str, Any]]] = {}
        for visual in self._gene_visuals:
            try:
                row_idx = int(visual.get("row_idx") or 0)
            except Exception:
                continue
            by_row.setdefault(row_idx, []).append(visual)
        for visuals in by_row.values():
            visuals.sort(key=lambda v: int(v.get("gene_idx") or 0))
        return by_row

    def _match_adjacent_row_visuals(
        self,
        upper: Sequence[Dict[str, Any]],
        lower: Sequence[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        def relative_index(visual: Dict[str, Any], count: int) -> float:
            if count <= 1:
                return 0.5
            try:
                return float(int(visual.get("gene_idx") or 0)) / float(max(1, count - 1))
            except Exception:
                return 0.5

        upper_by_key: Dict[str, List[Dict[str, Any]]] = {}
        lower_by_key: Dict[str, List[Dict[str, Any]]] = {}
        for visual in upper:
            key = str(visual.get("homology_key") or "")
            if key:
                upper_by_key.setdefault(key, []).append(visual)
        for visual in lower:
            key = str(visual.get("homology_key") or "")
            if key:
                lower_by_key.setdefault(key, []).append(visual)

        upper_count = len(upper)
        lower_count = len(lower)
        matches: List[Dict[str, Any]] = []
        used_lower: set[int] = set()
        for key, upper_candidates in upper_by_key.items():
            lower_candidates = lower_by_key.get(key, [])
            if not lower_candidates:
                continue
            pair_costs: List[Tuple[float, Dict[str, Any], Dict[str, Any]]] = []
            for up in upper_candidates:
                up_rel = relative_index(up, upper_count)
                up_idx = int(up.get("gene_idx") or 0)
                for lw in lower_candidates:
                    if id(lw) in used_lower:
                        continue
                    low_rel = relative_index(lw, lower_count)
                    low_idx = int(lw.get("gene_idx") or 0)
                    strand_penalty = (
                        0.12
                        if self._normalized_strand(up.get("strand"))
                        != self._normalized_strand(lw.get("strand"))
                        else 0.0
                    )
                    order_penalty = abs(up_rel - low_rel)
                    index_penalty = abs(up_idx - low_idx) / float(max(1, max(upper_count, lower_count) - 1))
                    center_penalty = abs(
                        float(up.get("x_center") or 0.0) - float(lw.get("x_center") or 0.0)
                    ) / max(1.0, float(self._max_rel_bp - self._min_rel_bp))
                    pair_costs.append(
                        (
                            (order_penalty * 2.0) + (index_penalty * 0.6) + (center_penalty * 0.35) + strand_penalty,
                            up,
                            lw,
                        )
                    )

            used_upper: set[int] = set()
            for _cost, up, lw in sorted(pair_costs, key=lambda item: item[0]):
                if id(up) in used_upper or id(lw) in used_lower:
                    continue
                used_upper.add(id(up))
                used_lower.add(id(lw))
                matches.append(
                    {
                        "upper": up,
                        "lower": lw,
                        "key": key,
                        "evidence": str(up.get("homology_evidence") or lw.get("homology_evidence") or ""),
                        "weight": max(float(up.get("homology_weight") or 0.0), float(lw.get("homology_weight") or 0.0)),
                    }
                )

        matches.sort(
            key=lambda match: (
                int(match["upper"].get("gene_idx") or 0),
                int(match["lower"].get("gene_idx") or 0),
            )
        )
        return matches

    def _collapse_synteny_matches(self, matches: Sequence[Dict[str, Any]]) -> List[Tuple[List[Dict[str, Any]], int]]:
        blocks: List[Tuple[List[Dict[str, Any]], int]] = []
        current: List[Dict[str, Any]] = []
        direction = 0

        def flush() -> None:
            nonlocal current, direction
            if current:
                blocks.append((current, direction or 1))
            current = []
            direction = 0

        for match in matches:
            if not current:
                current = [match]
                direction = 0
                continue

            previous = current[-1]
            prev_up = int(previous["upper"].get("gene_idx") or 0)
            prev_low = int(previous["lower"].get("gene_idx") or 0)
            cur_up = int(match["upper"].get("gene_idx") or 0)
            cur_low = int(match["lower"].get("gene_idx") or 0)
            up_delta = cur_up - prev_up
            low_delta = cur_low - prev_low
            next_direction = 1 if low_delta > 0 else -1 if low_delta < 0 else 0
            can_extend = (
                1 <= up_delta <= self._max_synteny_gap_genes + 1
                and 1 <= abs(low_delta) <= self._max_synteny_gap_genes + 1
                and next_direction != 0
                and (direction == 0 or direction == next_direction)
            )
            if can_extend:
                current.append(match)
                direction = direction or next_direction
            else:
                flush()
                current = [match]
                direction = 0
        flush()
        return blocks

    @staticmethod
    def _visual_span(visuals: Sequence[Dict[str, Any]]) -> Tuple[float, float]:
        xs: List[float] = []
        for visual in visuals:
            xs.extend([float(visual.get("x0") or 0.0), float(visual.get("x1") or 0.0)])
        if not xs:
            return 0.0, 1.0
        x0 = min(xs)
        x1 = max(xs)
        if x1 - x0 < 8.0:
            center = (x0 + x1) / 2.0
            x0 = center - 4.0
            x1 = center + 4.0
        return x0, x1

    @staticmethod
    def _synteny_labels_for_block(block: Sequence[Dict[str, Any]]) -> List[str]:
        labels: List[str] = []
        seen: set[str] = set()
        for match in block:
            label = str(match["upper"].get("label") or "").strip()
            if not label:
                label = str(match["key"]).split(":", 1)[-1]
            if label and label not in seen:
                seen.add(label)
                labels.append(label)
        return labels

    def _add_synteny_ribbon(self, block: Sequence[Dict[str, Any]], direction: int) -> None:
        if not block:
            return
        upper_visuals = [match["upper"] for match in block]
        lower_visuals = [match["lower"] for match in block]
        upper_left, upper_right = self._visual_span(upper_visuals)
        lower_left, lower_right = self._visual_span(lower_visuals)
        block_len = len(block)
        pad = 1.2 if block_len >= 2 else 0.6
        upper_left -= pad
        upper_right += pad
        lower_left -= pad
        lower_right += pad

        try:
            upper_row = int(upper_visuals[0].get("row_idx") or 0)
            lower_row = int(lower_visuals[0].get("row_idx") or upper_row + 1)
        except Exception:
            return
        y_upper = self._top_margin + upper_row * (self._row_height + self._row_gap) + self._row_height - 0.4
        y_lower = self._top_margin + lower_row * (self._row_height + self._row_gap) + 0.4
        if y_lower <= y_upper:
            return

        strand_mismatch = any(
            self._normalized_strand(match["upper"].get("strand"))
            != self._normalized_strand(match["lower"].get("strand"))
            for match in block
        )
        inverted = direction < 0 or strand_mismatch
        if inverted:
            brush_color = QtGui.QColor(239, 68, 68, 72 if block_len >= 2 else 44)
        else:
            brush_color = QtGui.QColor(100, 116, 139, 54 if block_len >= 2 else 34)

        path = QtGui.QPainterPath()
        mid_y = (y_upper + y_lower) / 2.0
        if inverted:
            path.moveTo(upper_left, y_upper)
            path.cubicTo(upper_left, mid_y, lower_right, mid_y, lower_right, y_lower)
            path.lineTo(lower_left, y_lower)
            path.cubicTo(lower_left, mid_y, upper_right, mid_y, upper_right, y_upper)
        else:
            path.moveTo(upper_left, y_upper)
            path.cubicTo(upper_left, mid_y, lower_left, mid_y, lower_left, y_lower)
            path.lineTo(lower_right, y_lower)
            path.cubicTo(lower_right, mid_y, upper_right, mid_y, upper_right, y_upper)
        path.closeSubpath()

        item = QtWidgets.QGraphicsPathItem(path)
        item.setBrush(QtGui.QBrush(brush_color))
        item.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
        item.setZValue(9 if block_len >= 2 else 8)
        item.setData(2, "synteny")
        labels = self._synteny_labels_for_block(block)
        label_preview = ", ".join(labels[:6])
        if len(labels) > 6:
            label_preview += f", +{len(labels) - 6} more"
        if block_len >= 2:
            title = "Inverted synteny block" if inverted else "Synteny block"
        else:
            title = "Inverted homology link" if inverted else "Homology link"
        evidence_counts: Dict[str, int] = {}
        matched_genes: List[Dict[str, Any]] = []
        for match in block:
            evidence = str(match.get("evidence") or match["upper"].get("homology_evidence") or "unknown")
            evidence_counts[evidence] = evidence_counts.get(evidence, 0) + 1
            matched_genes.append(
                {
                    "key": str(match.get("key") or ""),
                    "label": str(match["upper"].get("label") or _homology_display_for_key(str(match.get("key") or ""))),
                    "upper_row": int(match["upper"].get("row_idx") or 0),
                    "upper_gene": int(match["upper"].get("gene_idx") or 0),
                    "upper_strand": self._normalized_strand(match["upper"].get("strand")),
                    "lower_row": int(match["lower"].get("row_idx") or 0),
                    "lower_gene": int(match["lower"].get("gene_idx") or 0),
                    "lower_strand": self._normalized_strand(match["lower"].get("strand")),
                    "evidence": evidence,
                }
            )
        evidence_text = ", ".join(
            f"{name} x{count}" if count > 1 else name
            for name, count in sorted(evidence_counts.items(), key=lambda item: (-item[1], item[0]))
        )
        block_info = {
            "title": title,
            "matched_gene_count": int(block_len),
            "orientation": "inverted" if inverted else "same orientation",
            "direction": int(direction),
            "evidence_counts": evidence_counts,
            "evidence_summary": evidence_text,
            "labels": labels,
            "label_preview": label_preview,
            "matches": matched_genes,
        }
        tooltip = f"{title}: {block_len} matched genes"
        if evidence_text:
            tooltip += f"\nEvidence: {evidence_text}"
        if label_preview:
            tooltip += f"\n{label_preview}"
        item.setToolTip(tooltip)
        item.setData(4, block_info)
        self._scene.addItem(item)
        self._synteny_band_items.append(item)

    def _render_synteny_layer(self) -> None:
        for item in self._synteny_band_items:
            try:
                self._scene.removeItem(item)
            except Exception:
                pass
        self._synteny_band_items = []
        if not self._show_synteny_blocks or len(self._rows) < 2:
            return

        by_row = self._row_visuals_by_index()
        for row_idx in range(len(self._rows) - 1):
            upper = by_row.get(row_idx, [])
            lower = by_row.get(row_idx + 1, [])
            if not upper or not lower:
                continue
            matches = self._match_adjacent_row_visuals(upper, lower)
            for block, direction in self._collapse_synteny_matches(matches):
                self._add_synteny_ribbon(block, direction)

    @staticmethod
    def _gene_key(payload: Optional[Dict[str, Any]], fallback_label: str = "") -> str:
        gene_name = ""
        if isinstance(payload, dict):
            gene_name = str(payload.get("Gene") or "").strip()
            if not gene_name:
                gene_name = str(payload.get("Product") or "").strip()
            if not gene_name:
                gene_name = str(payload.get("label") or "").strip()
        if not gene_name:
            gene_name = str(fallback_label).strip()
        return gene_name.lower()

    @staticmethod
    def _ordered_hover_samples(raw: Optional[Sequence[str]]) -> List[str]:
        if not raw:
            return []
        out: List[str] = []
        seen: set[str] = set()
        for value in raw:
            sample = str(value).strip()
            if not sample or sample in seen:
                continue
            seen.add(sample)
            out.append(sample)
        return out

    def _hover_target_from_item(
        self,
        item: Optional[QtWidgets.QGraphicsItem],
    ) -> Tuple[Optional[int], Optional[str]]:
        cur = item
        while cur is not None:
            try:
                row_idx_raw = cur.data(1)
            except Exception:
                row_idx_raw = None
            try:
                hover_kind = cur.data(2)
            except Exception:
                hover_kind = None
            if row_idx_raw is not None:
                try:
                    return int(row_idx_raw), str(hover_kind or "").strip() or None
                except Exception:
                    pass
            cur = cur.parentItem()
        return None, None

    def _row_index_for_scene_y(self, scene_y: float) -> Optional[int]:
        rel_y = float(scene_y) - self._top_margin
        if rel_y < 0.0:
            return None
        row_span = self._row_height + self._row_gap
        if row_span <= 0.0:
            return None
        row_idx = int(rel_y // row_span)
        if row_idx < 0 or row_idx >= len(self._rows):
            return None
        offset = rel_y - float(row_idx) * row_span
        if 0.0 <= offset <= self._row_height:
            return row_idx
        return None

    def _hover_samples_for_row(self, row: Dict[str, Any], hover_kind: Optional[str]) -> List[str]:
        collapsed_samples = self._ordered_hover_samples(row.get("collapsed_samples"))
        if len(collapsed_samples) > 1:
            return collapsed_samples
        sample = str(row.get("sample") or "").strip()
        if sample:
            return [sample]
        return collapsed_samples

    def _emit_context_hover(self, samples: Optional[Sequence[str]]) -> None:
        normalized = self._ordered_hover_samples(samples)
        key = tuple(normalized) if normalized else None
        if key == self._hovered_samples_key:
            return
        self._hovered_samples_key = key
        self.contextSamplesHovered.emit(normalized or None)

    def _on_scrollbar_changed(self, _value: int) -> None:
        self._reposition_row_labels()
        self._reposition_figure_legend()
        self.viewportChanged.emit()

    def _reposition_row_labels(self) -> None:
        if not self._row_label_items:
            return
        scale_x = self._scale_x()
        scale_y = self._scale_y()
        label_x_px = self._row_label_left_padding_px if self._static_row_labels else 8.0
        track_lane_px = self._computed_track_x_offset if self._static_row_labels else 0.0
        label_gap_px = self._row_label_track_gap_px if self._static_row_labels else 8.0
        left_scene = self.mapToScene(QtCore.QPoint(0, 0)).x()
        x = float(left_scene + (label_x_px / scale_x))
        if self._static_row_labels:
            for item, _y_top, x_offset in self._row_label_items:
                is_title = item.data(4) == "title" or float(x_offset) > 0.0
                item.setVisible(scale_y >= (0.32 if is_title else 0.12))
                if is_title:
                    full_text = str(item.data(5) or item.text())
                    available_px = max(
                        36.0,
                        float(track_lane_px) - label_x_px - float(x_offset) - label_gap_px,
                    )
                    metrics = QtGui.QFontMetrics(item.font())
                    item.setText(
                        metrics.elidedText(
                            full_text,
                            QtCore.Qt.TextElideMode.ElideRight,
                            int(available_px),
                        )
                    )
                item.setPos(x + (float(x_offset) / scale_x), float(_y_top - 12.0))
            return
        for item, y_top, x_offset in self._row_label_items:
            is_title = item.data(4) == "title" or float(x_offset) > 0.0
            item.setVisible(scale_y >= (0.32 if is_title else 0.12))
            item.setPos(x + (float(x_offset) / scale_x), float(y_top - 12.0))

    def _position_gene_label_item(
        self,
        item: QtWidgets.QGraphicsSimpleTextItem,
        x_center: float,
        y_top: float,
        *,
        scale_x: Optional[float] = None,
        scale_y: Optional[float] = None,
    ) -> None:
        if scale_x is None:
            scale_x = abs(float(self.transform().m11()))
        if scale_y is None:
            scale_y = abs(float(self.transform().m22()))
        scale_x = max(float(scale_x or 1.0), 1e-6)
        scale_y = max(float(scale_y or 1.0), 1e-6)
        rect = item.boundingRect()
        item.setTransformOriginPoint(rect.width() / 2.0, rect.height())
        item.setPos(
            float(x_center) - (float(rect.width()) / 2.0 / scale_x),
            float(y_top) - ((float(rect.height()) + self._gene_label_gap_px) / scale_y),
        )

    def _update_gene_label_elision(self) -> None:
        if not self._gene_label_items:
            return
        scale_x = max(abs(float(self.transform().m11())), self._current_min_zoom_x())
        scale_y = abs(float(self.transform().m22()))
        diagonal = 0.70710678118
        for item, full_text, scene_width, x_center, y_top in self._gene_label_items:
            metrics = QtGui.QFontMetrics(item.font())
            try:
                row_idx = int(item.data(1) or 0)
            except Exception:
                row_idx = 0
            available_scene = self._top_margin if row_idx <= 0 else self._row_gap
            available_above_px = max(0.0, (available_scene * scale_y) - 12.0)
            horizontal_anchor_px = float(scene_width) * scale_x
            max_from_vertical_px = (available_above_px - (metrics.height() * diagonal)) / diagonal
            max_px = int(
                max(
                    0.0,
                    min(max_from_vertical_px, max(28.0, horizontal_anchor_px * 3.5)),
                )
            )
            if available_above_px < 28.0 or horizontal_anchor_px < 3.0 or max_px < 20:
                item.setText("")
                continue
            item.setText(metrics.elidedText(full_text, QtCore.Qt.TextElideMode.ElideRight, max_px))
            self._position_gene_label_item(
                item,
                float(x_center),
                float(y_top),
                scale_x=scale_x,
                scale_y=scale_y,
            )

    def _update_synteny_overview_opacity(self) -> None:
        if not self._synteny_band_items:
            return
        scale_y = max(0.0, min(self._max_zoom_y, abs(float(self.transform().m22()))))
        opacity = max(0.48, min(1.0, 0.42 + (scale_y * 0.58)))
        for item in self._synteny_band_items:
            item.setOpacity(opacity)

    def _add_figure_legend_item(self, item: QtWidgets.QGraphicsItem, x_offset: float, y_offset: float) -> None:
        item.setFlag(QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations, True)
        self._scene.addItem(item)
        self._figure_legend_items.append((item, float(x_offset), float(y_offset)))

    def _render_figure_legend(self) -> None:
        self._figure_legend_items = []
        self._figure_legend_size = QtCore.QSizeF(0.0, 0.0)
        if self._display_mode != "figure" or not self._figure_legend_entries:
            return
        entries = sorted(
            self._figure_legend_entries.values(),
            key=lambda style: (
                _FIGURE_CATEGORY_SORT.get(str(style.get("category") or ""), 99),
                str(style.get("family_label") or style.get("family") or "").casefold(),
                str(style.get("label") or "").casefold(),
            ),
        )

        title_font = QtGui.QFont()
        title_font.setPointSize(8)
        title_font.setBold(True)
        label_font = QtGui.QFont()
        label_font.setPointSize(8)
        title_metrics = QtGui.QFontMetrics(title_font)
        label_metrics = QtGui.QFontMetrics(label_font)
        max_rows_per_column = 30
        column_count = max(1, (len(entries) + max_rows_per_column - 1) // max_rows_per_column)
        column_widths: List[float] = []
        for column_idx in range(column_count):
            col_entries = entries[
                column_idx * max_rows_per_column : (column_idx + 1) * max_rows_per_column
            ]
            widest_label = max(
                [label_metrics.horizontalAdvance(str(style.get("label") or "")) for style in col_entries] or [0]
            )
            column_widths.append(max(120.0, float(widest_label + 42.0)))
        title_width = float(title_metrics.horizontalAdvance(f"Figure legend ({len(entries)})"))
        column_gap = 12.0
        width = max(title_width + 20.0, sum(column_widths) + column_gap * float(max(0, column_count - 1)) + 20.0)
        row_h = 24.0
        visible_rows = min(max_rows_per_column, len(entries))
        height = 32.0 + row_h * float(visible_rows)
        self._figure_legend_size = QtCore.QSizeF(width, height)

        bg = QtWidgets.QGraphicsRectItem(0.0, 0.0, width, height)
        bg.setBrush(QtGui.QBrush(QtGui.QColor(255, 255, 255, 226)))
        bg.setPen(QtGui.QPen(QtGui.QColor("#cbd5e1"), 1.0))
        bg.setZValue(86)
        self._add_figure_legend_item(bg, 0.0, 0.0)

        title = QtWidgets.QGraphicsSimpleTextItem(f"Figure legend ({len(entries)})")
        title.setFont(title_font)
        title.setBrush(QtGui.QBrush(QtGui.QColor("#334155")))
        title.setZValue(87)
        self._add_figure_legend_item(title, 10.0, 7.0)

        x = 10.0
        for column_idx in range(column_count):
            y = 28.0
            col_entries = entries[
                column_idx * max_rows_per_column : (column_idx + 1) * max_rows_per_column
            ]
            for style in col_entries:
                color = style.get("color") if isinstance(style.get("color"), QtGui.QColor) else QtGui.QColor("#94a3b8")
                tooltip = (
                    f"{style.get('label') or 'gene'} | "
                    f"{style.get('category') or 'marker'} | "
                    f"{style.get('family_label') or style.get('family') or 'family'} "
                    f"({int(style.get('count') or 1)} occurrence(s))"
                )
                if bool(style.get("ice")):
                    tooltip += " | CGE ICE-associated"
                symbol = QtWidgets.QGraphicsPathItem(
                    self._figure_symbol_path(str(style.get("shape") or "circle"), QtCore.QRectF(0.0, 0.0, 16.0, 16.0))
                )
                symbol.setBrush(QtGui.QBrush(color))
                symbol.setPen(QtGui.QPen(QtGui.QColor("#0f172a"), 0.8))
                symbol.setZValue(87)
                symbol.setToolTip(tooltip)
                self._add_figure_legend_item(symbol, x + 1.0, y + 1.0)

                if bool(style.get("ice")):
                    ice_font = QtGui.QFont()
                    ice_font.setPointSize(5)
                    ice_font.setBold(True)
                    ice_text = QtWidgets.QGraphicsSimpleTextItem("ICE")
                    ice_text.setFont(ice_font)
                    ice_text.setBrush(QtGui.QBrush(QtGui.QColor("#713f12")))
                    ice_text.setZValue(88)
                    ice_text.setToolTip(tooltip)
                    ice_metrics = QtGui.QFontMetrics(ice_font)
                    self._add_figure_legend_item(
                        ice_text,
                        x + 9.0 - (float(ice_metrics.horizontalAdvance("ICE")) / 2.0),
                        y + 9.0 - (float(ice_metrics.height()) / 2.0),
                    )

                text = QtWidgets.QGraphicsSimpleTextItem(str(style.get("label") or "gene"))
                text.setFont(label_font)
                text.setBrush(QtGui.QBrush(QtGui.QColor("#334155")))
                text.setZValue(87)
                text.setToolTip(tooltip)
                self._add_figure_legend_item(text, x + 28.0, y + 1.0)
                y += row_h
            x += column_widths[column_idx] + column_gap

    def _reposition_figure_legend(self) -> None:
        if not self._figure_legend_items or self._figure_legend_size.isEmpty():
            return
        scale_x = abs(float(self.transform().m11()))
        scale_y = abs(float(self.transform().m22()))
        if scale_x <= 1e-6:
            scale_x = 1.0
        if scale_y <= 1e-6:
            scale_y = 1.0
        viewport = self.viewport().rect()
        anchor = self.mapToScene(QtCore.QPoint(max(0, viewport.width() - 14), 14))
        base_x = float(anchor.x()) - (float(self._figure_legend_size.width()) / scale_x)
        base_y = float(anchor.y())
        for item, x_offset, y_offset in self._figure_legend_items:
            item.setPos(base_x + (x_offset / scale_x), base_y + (y_offset / scale_y))

    def _update_zoom_dependent_items(self) -> None:
        self._reposition_row_labels()
        self._update_gene_label_elision()
        self._update_synteny_overview_opacity()
        self._reposition_figure_legend()

    def _current_min_zoom_x(self) -> float:
        return self._figure_min_zoom_x if self._display_mode == "figure" else self._min_zoom_x

    def _apply_zoom_step(self, *, zoom_out: bool) -> bool:
        zoom_factor = 1.12
        current_x = max(abs(float(self.transform().m11())), 1e-6)
        current_y = max(abs(float(self.transform().m22())), 1e-6)
        min_zoom_x = self._current_min_zoom_x()

        if zoom_out:
            if current_x > min_zoom_x + 1e-4:
                target_x = max(min_zoom_x, current_x / zoom_factor)
                self.scale(target_x / current_x, 1.0)
                if self._static_row_labels:
                    self._render_scene()
                    self._center_on_origin()
                else:
                    self._update_zoom_dependent_items()
                self.viewportChanged.emit()
                return True
            if current_y > self._min_zoom_y + 1e-4:
                target_y = max(self._min_zoom_y, current_y / zoom_factor)
                self.scale(1.0, target_y / current_y)
                self._update_zoom_dependent_items()
                self.viewportChanged.emit()
                return True
            return False

        if current_y < self._max_zoom_y - 1e-4:
            target_y = min(self._max_zoom_y, current_y * zoom_factor)
            self.scale(1.0, target_y / current_y)
            self._update_zoom_dependent_items()
            self.viewportChanged.emit()
            return True
        if current_x < self._max_zoom_x - 1e-4:
            target_x = min(self._max_zoom_x, current_x * zoom_factor)
            self.scale(target_x / current_x, 1.0)
            if self._static_row_labels:
                self._render_scene()
                self._center_on_origin()
            else:
                self._update_zoom_dependent_items()
            self.viewportChanged.emit()
            return True
        return False

    def _row_center_y(self, row_idx: int) -> float:
        y_top = self._top_margin + int(row_idx) * (self._row_height + self._row_gap)
        return y_top + self._row_height / 2.0

    def _clear_hover_highlights(self) -> None:
        for item in self._hover_highlight_items:
            try:
                self._scene.removeItem(item)
            except Exception:
                pass
        self._hover_highlight_items = []

    def _add_hover_highlight(self, visual: Dict[str, Any]) -> None:
        if bool(visual.get("fixed_symbol")):
            cx = float(visual.get("x_center") or 0.0)
            cy = self._row_center_y(int(visual.get("row_idx") or 0))
            size = float(visual.get("fixed_symbol_size") or 34.0) + 12.0
            rect = QtCore.QRectF(-size / 2.0, -size / 2.0, size, size)
            halo = QtWidgets.QGraphicsRectItem(rect)
            halo.setFlag(QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations, True)
            halo.setPos(cx, cy)
            halo.setBrush(QtGui.QBrush(QtGui.QColor(6, 182, 212, 54)))
            halo.setPen(QtGui.QPen(QtGui.QColor("#06b6d4"), 3.2))
            halo.setZValue(43)
            halo.setData(2, "hover-highlight")
            halo.setData(3, str(visual.get("homology_key") or ""))
            halo.setData(1, int(visual.get("row_idx") or 0))
            self._scene.addItem(halo)
            self._hover_highlight_items.append(halo)
            return
        x0 = float(visual.get("x0") or 0.0)
        x1 = float(visual.get("x1") or x0)
        y_top = float(visual.get("y_top") or 0.0)
        if x1 <= x0:
            return
        rect = QtCore.QRectF(x0 - 7.0, y_top - 7.0, (x1 - x0) + 14.0, self._row_height + 14.0)
        halo = QtWidgets.QGraphicsRectItem(rect)
        halo.setBrush(QtGui.QBrush(QtGui.QColor(6, 182, 212, 54)))
        halo.setPen(QtGui.QPen(QtGui.QColor("#06b6d4"), 3.2))
        halo.setZValue(43)
        halo.setData(2, "hover-highlight")
        halo.setData(3, str(visual.get("homology_key") or ""))
        halo.setData(1, int(visual.get("row_idx") or 0))
        self._scene.addItem(halo)
        self._hover_highlight_items.append(halo)

    def _refresh_visual_emphasis(self) -> None:
        self._clear_hover_highlights()
        gene_key = self._selected_gene_key
        homology_key = self._hover_homology_key
        has_selection = bool(gene_key)
        has_hover = bool(homology_key)
        for visual in self._gene_visuals:
            is_selected = has_selection and visual.get("gene_key") == gene_key
            is_hovered = has_hover and visual.get("homology_key") == homology_key
            if is_hovered:
                self._add_hover_highlight(visual)
            shape_items = visual.get("shape_items") or []
            base_pens = visual.get("base_pens") or []
            for idx, item in enumerate(shape_items):
                if idx < len(base_pens):
                    pen = QtGui.QPen(base_pens[idx])
                else:
                    try:
                        pen = QtGui.QPen(item.pen())  # type: ignore[attr-defined]
                    except Exception:
                        pen = QtGui.QPen()
                if is_hovered:
                    pen.setColor(QtGui.QColor("#00e5ff"))
                    pen.setWidthF(max(3.4, pen.widthF() + 2.4))
                    item.setOpacity(1.0)
                    item.setZValue(44)
                elif is_selected:
                    pen.setColor(QtGui.QColor("#f97316"))
                    pen.setWidthF(max(1.2, pen.widthF() + 0.8))
                    item.setOpacity(1.0)
                    item.setZValue(42)
                elif has_selection:
                    item.setOpacity(0.18)
                    item.setZValue(24)
                else:
                    item.setOpacity(1.0)
                    item.setZValue(24)
                try:
                    item.setPen(pen)  # type: ignore[attr-defined]
                except Exception:
                    pass
            text_item = visual.get("text_item")
            if text_item is not None:
                if is_hovered:
                    hover_font = QtGui.QFont(text_item.font())
                    hover_font.setBold(True)
                    text_item.setFont(hover_font)
                    text_item.setBrush(QtGui.QBrush(QtGui.QColor("#075985")))
                    text_item.setOpacity(1.0)
                    text_item.setZValue(50)
                elif is_selected:
                    base_font = QtGui.QFont(text_item.font())
                    base_font.setBold(False)
                    text_item.setFont(base_font)
                    text_item.setBrush(QtGui.QBrush(QtGui.QColor("#c2410c")))
                    text_item.setOpacity(1.0)
                    text_item.setZValue(48)
                elif has_selection:
                    base_font = QtGui.QFont(text_item.font())
                    base_font.setBold(False)
                    text_item.setFont(base_font)
                    text_item.setBrush(QtGui.QBrush(QtGui.QColor("#94a3b8")))
                    text_item.setOpacity(0.55)
                    text_item.setZValue(30)
                else:
                    base_font = QtGui.QFont(text_item.font())
                    base_font.setBold(False)
                    text_item.setFont(base_font)
                    text_item.setBrush(QtGui.QBrush(visual.get("text_base_color") or QtGui.QColor("#0f172a")))
                    text_item.setOpacity(1.0)
                    text_item.setZValue(30)
            for overlay_item in visual.get("overlay_items") or []:
                if is_hovered:
                    overlay_item.setOpacity(1.0)
                    overlay_item.setZValue(50)
                elif is_selected:
                    overlay_item.setOpacity(1.0)
                    overlay_item.setZValue(48)
                elif has_selection:
                    overlay_item.setOpacity(0.35)
                    overlay_item.setZValue(35)
                else:
                    overlay_item.setOpacity(1.0)
                    overlay_item.setZValue(35)

    def _apply_gene_selection(self, gene_key: Optional[str]) -> None:
        for line in self._selection_line_items:
            try:
                self._scene.removeItem(line)
            except Exception:
                pass
        self._selection_line_items = []
        self._selected_gene_key = gene_key

        has_selection = bool(gene_key)
        if not has_selection:
            self._refresh_visual_emphasis()
            return

        matches = [v for v in self._gene_visuals if v.get("gene_key") == gene_key]
        if len(matches) < 2:
            self._refresh_visual_emphasis()
            return
        by_row: Dict[int, List[Dict[str, Any]]] = {}
        for visual in matches:
            by_row.setdefault(int(visual.get("row_idx") or 0), []).append(visual)
        rows = sorted(by_row.keys())
        if len(rows) < 2:
            self._refresh_visual_emphasis()
            return
        for i in range(len(rows) - 1):
            upper = by_row[rows[i]]
            lower = by_row[rows[i + 1]]
            if not upper or not lower:
                continue
            for up in upper:
                x_up = float(up.get("x_center") or 0.0)
                best = min(lower, key=lambda lw: abs(float(lw.get("x_center") or 0.0) - x_up))
                x_down = float(best.get("x_center") or 0.0)
                y_up = self._row_center_y(rows[i])
                y_down = self._row_center_y(rows[i + 1])
                line = QtWidgets.QGraphicsLineItem(x_up, y_up, x_down, y_down)
                line.setPen(QtGui.QPen(QtGui.QColor("#fb7185"), 1.5, QtCore.Qt.PenStyle.SolidLine))
                line.setZValue(26)
                self._scene.addItem(line)
                self._selection_line_items.append(line)
        self._refresh_visual_emphasis()

    def _apply_homology_hover(self, homology_key: Optional[str]) -> None:
        # Gene-hover emphasis was too visually disruptive. Keep hover handling
        # for track/sample metadata, but never highlight/fade genes on hover.
        if self._hover_homology_key is None and not self._hover_highlight_items:
            return
        self._hover_homology_key = None
        self._refresh_visual_emphasis()

    def _homology_key_from_item(self, item: Optional[QtWidgets.QGraphicsItem]) -> str:
        cur = item
        while cur is not None:
            try:
                key = str(cur.data(3) or "").strip()
            except Exception:
                key = ""
            if key:
                return key
            cur = cur.parentItem()
        return ""

    def _synteny_payload_from_item(self, item: Optional[QtWidgets.QGraphicsItem]) -> Optional[Dict[str, Any]]:
        cur = item
        while cur is not None:
            try:
                kind = str(cur.data(2) or "").strip()
            except Exception:
                kind = ""
            if kind == "synteny":
                try:
                    payload = cur.data(4)
                except Exception:
                    payload = None
                if isinstance(payload, dict):
                    return dict(payload)
            cur = cur.parentItem()
        return None

    def _payload_from_item(self, item: Optional[QtWidgets.QGraphicsItem]) -> Optional[Dict[str, Any]]:
        cur = item
        while cur is not None:
            try:
                payload = cur.data(0)
            except Exception:
                payload = None
            if isinstance(payload, dict):
                return payload
            cur = cur.parentItem()
        return None

    def mousePressEvent(self, event: QtGui.QMouseEvent) -> None:
        pos = event.position().toPoint() if hasattr(event, "position") else event.pos()
        item = self._scene.itemAt(self.mapToScene(pos), QtGui.QTransform())
        synteny_payload = self._synteny_payload_from_item(item)
        if synteny_payload is not None:
            self.syntenyBlockClicked.emit(synteny_payload)
            event.accept()
            return
        self._apply_gene_selection(None)
        super().mousePressEvent(event)

    def mouseDoubleClickEvent(self, event: QtGui.QMouseEvent) -> None:
        pos = event.position().toPoint() if hasattr(event, "position") else event.pos()
        item = self._scene.itemAt(self.mapToScene(pos), QtGui.QTransform())
        payload = self._payload_from_item(item)
        if payload is not None:
            self.geneDoubleClicked.emit(payload)
            event.accept()
            return
        super().mouseDoubleClickEvent(event)

    def mouseMoveEvent(self, event: QtGui.QMouseEvent) -> None:
        pos = event.position().toPoint() if hasattr(event, "position") else event.pos()
        scene_pos = self.mapToScene(pos)
        item = self._scene.itemAt(scene_pos, QtGui.QTransform())
        self._apply_homology_hover(None)
        row_idx, hover_kind = self._hover_target_from_item(item)
        if row_idx is None:
            row_idx = self._row_index_for_scene_y(scene_pos.y())
            if row_idx is not None:
                hover_kind = "track"
        if row_idx is None or row_idx < 0 or row_idx >= len(self._rows):
            self._emit_context_hover(None)
        else:
            self._emit_context_hover(self._hover_samples_for_row(self._rows[row_idx], hover_kind))
        super().mouseMoveEvent(event)

    def reset_zoom(self) -> None:
        self.resetTransform()
        if self._static_row_labels:
            self._render_scene()
            self._center_on_origin()
        else:
            self._center_on_origin()
            self._update_zoom_dependent_items()
        self.viewportChanged.emit()

    def fit_all_contexts(self) -> None:
        rect = self._scene.sceneRect()
        if rect.isNull():
            return
        self.fitInView(rect.adjusted(-12.0, -12.0, 12.0, 12.0), QtCore.Qt.AspectRatioMode.KeepAspectRatio)
        if self._static_row_labels:
            self._render_scene()
        else:
            self._update_zoom_dependent_items()
        self.viewportChanged.emit()

    def _export_scene_rect(self) -> QtCore.QRectF:
        rect = self._scene.itemsBoundingRect()
        if rect.isNull():
            rect = self._scene.sceneRect()
        return rect.adjusted(-24.0, -24.0, 24.0, 24.0)

    def export_png(self, path: Path, *, scale: float = 2.0) -> bool:
        rect = self._export_scene_rect()
        if rect.isNull():
            return False
        width = max(1, int(round(rect.width() * float(scale))))
        height = max(1, int(round(rect.height() * float(scale))))
        image = QtGui.QImage(width, height, QtGui.QImage.Format.Format_ARGB32_Premultiplied)
        image.fill(QtGui.QColor("#f8fafc"))
        painter = QtGui.QPainter(image)
        try:
            painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
            self._scene.render(painter, QtCore.QRectF(0.0, 0.0, float(width), float(height)), rect)
        finally:
            painter.end()
        return bool(image.save(str(path)))

    def export_pdf(self, path: Path) -> bool:
        rect = self._export_scene_rect()
        if rect.isNull():
            return False
        width = max(1.0, float(rect.width()))
        height = max(1.0, float(rect.height()))
        writer = QtGui.QPdfWriter(str(path))
        writer.setResolution(96)
        page_size = QtGui.QPageSize(
            QtCore.QSizeF(width, height),
            QtGui.QPageSize.Unit.Point,
            "BACTolomew synteny context",
        )
        writer.setPageSize(page_size)
        writer.setPageMargins(
            QtCore.QMarginsF(0.0, 0.0, 0.0, 0.0),
            QtGui.QPageLayout.Unit.Point,
        )
        painter = QtGui.QPainter(writer)
        try:
            painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
            self._scene.render(painter, QtCore.QRectF(0.0, 0.0, width, height), rect)
        finally:
            painter.end()
        return True

    def wheelEvent(self, event: QtGui.QWheelEvent) -> None:
        delta = event.angleDelta().y()
        if delta == 0:
            return super().wheelEvent(event)
        self._apply_zoom_step(zoom_out=delta < 0)
        event.accept()

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        super().resizeEvent(event)
        self._update_zoom_dependent_items()
        self.viewportChanged.emit()

    def leaveEvent(self, event: QtCore.QEvent) -> None:
        self._emit_context_hover(None)
        self._apply_homology_hover(None)
        super().leaveEvent(event)


class IgvCollectionDialog(QtWidgets.QDialog):
    geneDoubleClicked = QtCore.pyqtSignal(dict)

    def __init__(
        self,
        parent: Optional[QtWidgets.QWidget],
        *,
        rows: Sequence[Dict[str, Any]],
        title: str = "IGV_collection",
        output_dir: Optional[Path] = None,
        pan_reference_sample: Optional[str] = None,
        pan_fallback_tree: Optional[Path] = None,
    ) -> None:
        super().__init__(parent)
        self._all_rows: List[Dict[str, Any]] = [dict(r) for r in rows if isinstance(r, dict)]
        self._metadata_filters: Dict[str, List[str]] = {}
        self._metadata_columns = self._collect_metadata_columns()
        self._missing_label = "(missing)"
        self._collapsed = False
        self._output_dir = Path(output_dir) if output_dir is not None else None
        self._pan_reference_sample = str(pan_reference_sample).strip() if pan_reference_sample else None
        self._pan_fallback_tree = Path(pan_fallback_tree) if pan_fallback_tree is not None else None
        self._pan_core_copies: List[PanCoreAccessoryTab] = []
        self._pan_core_host: Optional[QtWidgets.QWidget] = None
        self._pan_core_host_layout: Optional[QtWidgets.QVBoxLayout] = None
        self._pan_core_placeholder: Optional[QtWidgets.QLabel] = None
        self._visible_rows: List[Dict[str, Any]] = []
        self._highlighted_context_samples: Optional[Tuple[str, ...]] = None

        self.setWindowTitle(title)
        self.resize(1760, 860)
        self._build_ui()
        self._render_metadata_filters()
        self._apply_filters()
        if self._pan_core_host is not None:
            QtCore.QTimer.singleShot(0, self._materialize_pan_core_copy)

    def _build_ui(self) -> None:
        root = QtWidgets.QVBoxLayout(self)
        main_split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        main_split.setChildrenCollapsible(False)
        root.addWidget(main_split, stretch=1)

        left = QtWidgets.QWidget(self)
        left.setMinimumWidth(260)
        left.setMaximumWidth(360)
        left_layout = QtWidgets.QVBoxLayout(left)
        left_layout.setContentsMargins(8, 8, 8, 8)
        left_layout.setSpacing(8)
        left_layout.setAlignment(QtCore.Qt.AlignmentFlag.AlignTop)

        title = QtWidgets.QLabel("Sample Filters")
        title.setStyleSheet("font-weight: 600; color: #2b3c4d;")
        left_layout.addWidget(title)

        self.sample_query = QtWidgets.QLineEdit(left)
        self.sample_query.setPlaceholderText("Search samples and metadata...")
        self.sample_query.textChanged.connect(self._apply_filters)
        left_layout.addWidget(self.sample_query)

        molecule_row = QtWidgets.QWidget(left)
        molecule_layout = QtWidgets.QHBoxLayout(molecule_row)
        molecule_layout.setContentsMargins(0, 0, 0, 0)
        molecule_layout.setSpacing(6)
        molecule_label = QtWidgets.QLabel("Molecule", molecule_row)
        self.molecule_filter_box = QtWidgets.QComboBox(molecule_row)
        self.molecule_filter_box.addItem("All", "")
        self.molecule_filter_box.addItem("Plasmid only", "plasmid")
        self.molecule_filter_box.addItem("Chromosome only", "chromosome")
        self.molecule_filter_box.currentIndexChanged.connect(self._apply_filters)
        molecule_layout.addWidget(molecule_label)
        molecule_layout.addWidget(self.molecule_filter_box, stretch=1)
        left_layout.addWidget(molecule_row)

        self.add_filter_btn = QtWidgets.QPushButton("Add filter", left)
        self.add_filter_btn.setEnabled(bool(self._metadata_columns))
        self.add_filter_btn.clicked.connect(self._show_add_filter_menu)
        left_layout.addWidget(self.add_filter_btn)

        self.filters_container = QtWidgets.QWidget(left)
        self.filters_container.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Preferred,
            QtWidgets.QSizePolicy.Policy.Maximum,
        )
        self.filters_layout = QtWidgets.QVBoxLayout(self.filters_container)
        self.filters_layout.setContentsMargins(0, 0, 0, 0)
        self.filters_layout.setSpacing(6)
        self.filters_layout.setAlignment(QtCore.Qt.AlignmentFlag.AlignTop)
        left_layout.addWidget(self.filters_container)

        sort_group = QtWidgets.QGroupBox("Sorting", left)
        sort_layout = QtWidgets.QFormLayout(sort_group)
        self.sort_column_box = QtWidgets.QComboBox(sort_group)
        self.sort_column_box.addItem("Sample", "__sample__")
        self.sort_column_box.addItem("Contig", "__contig__")
        self.sort_column_box.addItem("Context similarity", "__context_similarity__")
        for col in self._metadata_columns:
            self.sort_column_box.addItem(col, col)
        self.sort_order_box = QtWidgets.QComboBox(sort_group)
        self.sort_order_box.addItem("Ascending", "asc")
        self.sort_order_box.addItem("Descending", "desc")
        self.group_molecule_check = QtWidgets.QCheckBox("Group plasmid/chromosome", sort_group)
        self.group_molecule_check.setChecked(True)
        sort_col_idx = self.sort_column_box.findData("__context_similarity__")
        if sort_col_idx >= 0:
            self.sort_column_box.setCurrentIndex(sort_col_idx)
        sort_order_idx = self.sort_order_box.findData("desc")
        if sort_order_idx >= 0:
            self.sort_order_box.setCurrentIndex(sort_order_idx)
        self.sort_column_box.currentIndexChanged.connect(self._on_sort_column_changed)
        self.sort_order_box.currentIndexChanged.connect(self._apply_filters)
        self.group_molecule_check.toggled.connect(self._apply_filters)
        sort_layout.addRow("Sort by", self.sort_column_box)
        sort_layout.addRow("Order", self.sort_order_box)
        sort_layout.addRow("", self.group_molecule_check)
        left_layout.addWidget(sort_group)

        view_group = QtWidgets.QGroupBox("Context View", left)
        view_controls = QtWidgets.QVBoxLayout(view_group)
        view_controls.setContentsMargins(8, 8, 8, 8)
        view_controls.setSpacing(6)

        design_row = QtWidgets.QWidget(view_group)
        design_layout = QtWidgets.QHBoxLayout(design_row)
        design_layout.setContentsMargins(0, 0, 0, 0)
        design_layout.setSpacing(6)
        mode_label = QtWidgets.QLabel("Design", design_row)
        self.context_design_box = QtWidgets.QComboBox(design_row)
        self.context_design_box.addItem("Exploratory", "exploratory")
        self.context_design_box.addItem("Figure", "figure")
        design_layout.addWidget(mode_label)
        design_layout.addWidget(self.context_design_box, stretch=1)
        view_controls.addWidget(design_row)

        nav_row = QtWidgets.QWidget(view_group)
        nav_layout = QtWidgets.QHBoxLayout(nav_row)
        nav_layout.setContentsMargins(0, 0, 0, 0)
        nav_layout.setSpacing(6)
        self.reset_zoom_btn = QtWidgets.QPushButton("Reset zoom", nav_row)
        self.fit_contexts_btn = QtWidgets.QPushButton("Fit all", nav_row)
        nav_layout.addWidget(self.reset_zoom_btn)
        nav_layout.addWidget(self.fit_contexts_btn)
        view_controls.addWidget(nav_row)

        export_row = QtWidgets.QWidget(view_group)
        export_layout = QtWidgets.QHBoxLayout(export_row)
        export_layout.setContentsMargins(0, 0, 0, 0)
        export_layout.setSpacing(6)
        self.export_png_btn = QtWidgets.QPushButton("Export PNG", export_row)
        self.export_pdf_btn = QtWidgets.QPushButton("Export PDF", export_row)
        export_layout.addWidget(self.export_png_btn)
        export_layout.addWidget(self.export_pdf_btn)
        view_controls.addWidget(export_row)
        left_layout.addWidget(view_group)

        self.collapse_btn = QtWidgets.QPushButton("Collapse", left)
        self.collapse_btn.setCheckable(True)
        self.collapse_btn.toggled.connect(self._toggle_collapsed)
        left_layout.addWidget(self.collapse_btn)
        left_layout.addStretch(1)

        main_split.addWidget(left)

        right = QtWidgets.QWidget(self)
        right_layout = QtWidgets.QVBoxLayout(right)
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.setSpacing(6)
        self.info = QtWidgets.QLabel("")
        self.info.setStyleSheet("color:#6a7b8d; padding-left: 4px;")
        right_layout.addWidget(self.info)

        content_split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal, right)
        content_split.setChildrenCollapsible(False)
        right_layout.addWidget(content_split, stretch=1)

        view_panel = QtWidgets.QWidget(content_split)
        view_panel_layout = QtWidgets.QVBoxLayout(view_panel)
        view_panel_layout.setContentsMargins(0, 0, 0, 0)
        view_panel_layout.setSpacing(6)

        view_split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal, view_panel)
        view_split.setChildrenCollapsible(False)
        self.view = IgvCollectionView(view_split)
        self.view.geneDoubleClicked.connect(self.geneDoubleClicked.emit)
        self.view.contextSamplesHovered.connect(self._on_context_samples_hovered)
        self.view.syntenyBlockClicked.connect(self._show_synteny_block_details)
        self.context_design_box.currentIndexChanged.connect(self._on_context_design_changed)
        self.reset_zoom_btn.clicked.connect(self.view.reset_zoom)
        self.fit_contexts_btn.clicked.connect(self.view.fit_all_contexts)
        self.export_png_btn.clicked.connect(self._export_context_png)
        self.export_pdf_btn.clicked.connect(self._export_context_pdf)
        view_split.addWidget(self.view)

        self.synteny_detail_panel = QtWidgets.QFrame(view_split)
        self.synteny_detail_panel.setFrameShape(QtWidgets.QFrame.Shape.StyledPanel)
        self.synteny_detail_panel.setMinimumWidth(280)
        self.synteny_detail_panel.setMaximumWidth(420)
        detail_layout = QtWidgets.QVBoxLayout(self.synteny_detail_panel)
        detail_layout.setContentsMargins(10, 8, 10, 10)
        detail_layout.setSpacing(6)
        detail_header = QtWidgets.QHBoxLayout()
        detail_title = QtWidgets.QLabel("Synteny Block", self.synteny_detail_panel)
        detail_title.setStyleSheet("font-weight: 700; color: #334155;")
        detail_close = QtWidgets.QToolButton(self.synteny_detail_panel)
        detail_close.setText("x")
        detail_close.clicked.connect(self.synteny_detail_panel.hide)
        detail_header.addWidget(detail_title)
        detail_header.addStretch(1)
        detail_header.addWidget(detail_close)
        detail_layout.addLayout(detail_header)
        self.synteny_detail_text = QtWidgets.QTextBrowser(self.synteny_detail_panel)
        self.synteny_detail_text.setOpenExternalLinks(False)
        self.synteny_detail_text.setStyleSheet(
            "QTextBrowser { background: #ffffff; border: 1px solid #d8dee7; border-radius: 6px; padding: 6px; }"
        )
        detail_layout.addWidget(self.synteny_detail_text, stretch=1)
        view_split.addWidget(self.synteny_detail_panel)
        self.synteny_detail_panel.hide()
        view_split.setStretchFactor(0, 1)
        view_split.setStretchFactor(1, 0)
        view_split.setSizes([980, 300])
        view_panel_layout.addWidget(view_split, stretch=1)

        content_split.addWidget(view_panel)

        close_row = QtWidgets.QHBoxLayout()
        close_row.addStretch(1)
        close_btn = QtWidgets.QPushButton("Close")
        close_btn.clicked.connect(self.accept)
        close_row.addWidget(close_btn)
        right_layout.addLayout(close_row)
        main_split.addWidget(right)

        if self._output_dir is not None:
            pan_host = QtWidgets.QWidget(content_split)
            pan_host_layout = QtWidgets.QVBoxLayout(pan_host)
            pan_host_layout.setContentsMargins(0, 0, 0, 0)
            pan_host_layout.setSpacing(0)
            placeholder = QtWidgets.QLabel("Loading neoPan…", pan_host)
            placeholder.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
            placeholder.setStyleSheet("color: #6a7b8d; font-weight: 600;")
            pan_host_layout.addWidget(placeholder, stretch=1)
            self._pan_core_host = pan_host
            self._pan_core_host_layout = pan_host_layout
            self._pan_core_placeholder = placeholder
            content_split.addWidget(pan_host)
            main_split.setStretchFactor(0, 0)
            main_split.setStretchFactor(1, 1)
            content_split.setStretchFactor(0, 1)
            content_split.setStretchFactor(1, 1)
            main_split.setSizes([300, 1460])
            content_split.setSizes([720, 720])
        else:
            main_split.setStretchFactor(0, 0)
            main_split.setStretchFactor(1, 1)
            content_split.setStretchFactor(0, 1)
            main_split.setSizes([300, 1460])

    def _new_pan_core_copy(self, parent: Optional[QtWidgets.QWidget] = None) -> Optional[PanCoreAccessoryTab]:
        if self._output_dir is None:
            return None
        return PanCoreAccessoryTab(
            self._output_dir,
            None,
            reference_sample=self._pan_reference_sample,
            fallback_tree=self._pan_fallback_tree,
            preserve_branch_lengths=True,
            color_tree_by_fastgear_clade=True,
            optimize_tree_by_fastgear_clade=True,
            annotate_fastgear_clusters=True,
            use_fastgear_clade_backgrounds=True,
            show_heatmap=False,
            show_info_label=False,
            parent=parent,
        )

    def _materialize_pan_core_copy(self) -> None:
        if self._pan_core_copies or self._pan_core_host is None or self._pan_core_host_layout is None:
            return
        pan_panel = self._new_pan_core_copy(self._pan_core_host)
        if pan_panel is None:
            if self._pan_core_placeholder is not None:
                self._pan_core_placeholder.setText("neoPan unavailable.")
            return
        if self._pan_core_placeholder is not None:
            self._pan_core_placeholder.deleteLater()
            self._pan_core_placeholder = None
        self._pan_core_copies = [pan_panel]
        self._pan_core_host_layout.addWidget(pan_panel)
        self._set_pan_core_sample_filter([str(row.get("sample") or "") for row in self._visible_rows])
        if self._highlighted_context_samples is not None:
            self._set_pan_core_highlight(list(self._highlighted_context_samples))

    def _collect_metadata_columns(self) -> List[str]:
        cols: set[str] = set()
        for row in self._all_rows:
            meta = row.get("sample_meta")
            if not isinstance(meta, dict):
                continue
            for col in meta.keys():
                name = str(col).strip()
                if name:
                    cols.add(name)
        return sorted(cols)

    def _show_add_filter_menu(self) -> None:
        menu = QtWidgets.QMenu(self)
        for col in self._metadata_columns:
            label = f"{col} (filtered)" if col in self._metadata_filters else col
            action = menu.addAction(label)
            action.triggered.connect(lambda _checked=False, c=col: self._edit_metadata_filter(c))
        anchor = self.add_filter_btn.mapToGlobal(QtCore.QPoint(0, self.add_filter_btn.height()))
        menu.exec(anchor)

    def _on_sort_column_changed(self) -> None:
        if self.sort_column_box.currentData() == "__context_similarity__":
            desc_idx = self.sort_order_box.findData("desc")
            if desc_idx >= 0 and self.sort_order_box.currentData() != "desc":
                self.sort_order_box.blockSignals(True)
                self.sort_order_box.setCurrentIndex(desc_idx)
                self.sort_order_box.blockSignals(False)
        self._apply_filters()

    def _metadata_value_options(self, column: str) -> List[Tuple[str, str]]:
        values: set[str] = set()
        has_missing = False
        for row in self._all_rows:
            meta = row.get("sample_meta")
            if not isinstance(meta, dict):
                has_missing = True
                continue
            value = str(meta.get(column, "")).strip()
            if value:
                values.add(value)
            else:
                has_missing = True
        opts: List[Tuple[str, str]] = [(v, v) for v in sorted(values)]
        if has_missing:
            opts.insert(0, (self._missing_label, MISSING_FILTER_TOKEN))
        return opts

    def _edit_metadata_filter(self, column: str) -> None:
        values = self._metadata_value_options(column)
        current = list(self._metadata_filters.get(column, []))
        dlg = MetadataValueDialog(self, column, values, current)
        if dlg.exec() != QtWidgets.QDialog.DialogCode.Accepted:
            return
        selected = dlg.selected_values()
        if selected:
            self._metadata_filters[column] = selected
        else:
            self._metadata_filters.pop(column, None)
        self._render_metadata_filters()
        self._apply_filters()

    def _render_metadata_filters(self) -> None:
        while self.filters_layout.count():
            item = self.filters_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
        for column, values in self._metadata_filters.items():
            row = QtWidgets.QWidget(self.filters_container)
            row.setSizePolicy(
                QtWidgets.QSizePolicy.Policy.Preferred,
                QtWidgets.QSizePolicy.Policy.Fixed,
            )
            row_layout = QtWidgets.QHBoxLayout(row)
            row_layout.setContentsMargins(0, 0, 0, 0)
            row_layout.setSpacing(6)
            label = QtWidgets.QLabel(column)
            summary = QtWidgets.QLabel(self._format_filter_values(values))
            summary.setToolTip(self._format_filter_values(values, full=True))
            remove_btn = QtWidgets.QPushButton("Remove")
            remove_btn.clicked.connect(lambda _checked=False, c=column: self._remove_metadata_filter(c))
            row_layout.addWidget(label)
            row_layout.addWidget(summary, stretch=1)
            row_layout.addWidget(remove_btn)
            self.filters_layout.addWidget(row)
        self.filters_layout.addStretch(1)

    def _remove_metadata_filter(self, column: str) -> None:
        if column in self._metadata_filters:
            self._metadata_filters.pop(column, None)
            self._render_metadata_filters()
            self._apply_filters()

    @staticmethod
    def _normalized_samples(samples: Optional[Sequence[str]]) -> List[str]:
        out: List[str] = []
        seen: set[str] = set()
        for raw in samples or []:
            sample = str(raw).strip()
            if not sample or sample in seen:
                continue
            seen.add(sample)
            out.append(sample)
        return out

    def _set_pan_core_sample_filter(self, samples: Optional[Sequence[str]]) -> None:
        normalized = self._normalized_samples(samples)
        if not self._pan_core_copies:
            return
        for panel in self._pan_core_copies:
            panel.set_sample_filter(normalized or None)

    def _set_pan_core_highlight(self, samples: Optional[Sequence[str]]) -> None:
        normalized = self._normalized_samples(samples)
        key = tuple(normalized) if normalized else None
        if key == self._highlighted_context_samples:
            return
        self._highlighted_context_samples = key
        for panel in self._pan_core_copies:
            panel.highlight_sample(None)
            panel.highlight_samples(None)
            if not normalized:
                continue
            if len(normalized) == 1:
                panel.highlight_sample(normalized[0])
            else:
                panel.highlight_samples(normalized)

    def _on_context_samples_hovered(self, samples: object) -> None:
        if not isinstance(samples, (list, tuple)):
            self._set_pan_core_highlight(None)
            return
        self._set_pan_core_highlight([str(s) for s in samples])

    def _on_context_design_changed(self) -> None:
        box = getattr(self, "context_design_box", None)
        mode = "exploratory"
        if isinstance(box, QtWidgets.QComboBox):
            mode = str(box.currentData() or "exploratory")
        self.view.set_display_mode(mode)

    def _show_synteny_block_details(self, payload: Dict[str, Any]) -> None:
        title = html.escape(str(payload.get("title") or "Synteny block"))
        orientation = html.escape(str(payload.get("orientation") or ""))
        evidence = html.escape(str(payload.get("evidence_summary") or "unknown"))
        count = int(payload.get("matched_gene_count") or 0)
        labels = [str(v) for v in payload.get("labels") or [] if str(v).strip()]
        matches = payload.get("matches") if isinstance(payload.get("matches"), list) else []
        rows = [
            "<style>"
            "table{border-collapse:collapse;width:100%;font-size:11px;}"
            "td,th{border-bottom:1px solid #e2e8f0;padding:3px 4px;text-align:left;}"
            "th{color:#475569;background:#f8fafc;}"
            ".muted{color:#64748b;}"
            "</style>",
            f"<h3 style='margin:0 0 6px 0;color:#0f172a'>{title}</h3>",
            f"<p class='muted' style='margin:0 0 8px 0'>Genes: <b>{count}</b><br>Orientation: <b>{orientation}</b><br>Evidence: <b>{evidence}</b></p>",
        ]
        if labels:
            rows.append(
                "<p style='margin:0 0 8px 0'><b>Labels:</b> "
                + html.escape(", ".join(labels[:12]))
                + (f" <span class='muted'>(+{len(labels) - 12} more)</span>" if len(labels) > 12 else "")
                + "</p>"
            )
        if matches:
            rows.append("<table><tr><th>Gene</th><th>Rows</th><th>Strand</th><th>Evidence</th></tr>")
            for match in matches[:40]:
                label = html.escape(str(match.get("label") or _homology_display_for_key(str(match.get("key") or ""))))
                upper_row = int(match.get("upper_row") or 0) + 1
                lower_row = int(match.get("lower_row") or 0) + 1
                strand = html.escape(f"{match.get('upper_strand') or '?'} -> {match.get('lower_strand') or '?'}")
                ev = html.escape(str(match.get("evidence") or ""))
                rows.append(f"<tr><td>{label}</td><td>{upper_row} -> {lower_row}</td><td>{strand}</td><td>{ev}</td></tr>")
            rows.append("</table>")
            if len(matches) > 40:
                rows.append(f"<p class='muted'>Showing first 40 of {len(matches)} matched genes.</p>")
        self.synteny_detail_text.setHtml("".join(rows))
        self.synteny_detail_panel.show()

    def _default_export_path(self, suffix: str) -> str:
        base = self._output_dir if self._output_dir is not None else Path.home()
        return str(base / f"bactolomew_synteny_context.{suffix}")

    def _export_context_png(self) -> None:
        path_text, _selected = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Export synteny context as PNG",
            self._default_export_path("png"),
            "PNG image (*.png)",
        )
        if not path_text:
            return
        path = Path(path_text)
        if path.suffix.lower() != ".png":
            path = path.with_suffix(".png")
        if not self.view.export_png(path):
            QtWidgets.QMessageBox.warning(self, "Export failed", "Could not write the PNG export.")

    def _export_context_pdf(self) -> None:
        path_text, _selected = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Export synteny context as PDF",
            self._default_export_path("pdf"),
            "PDF document (*.pdf)",
        )
        if not path_text:
            return
        path = Path(path_text)
        if path.suffix.lower() != ".pdf":
            path = path.with_suffix(".pdf")
        if not self.view.export_pdf(path):
            QtWidgets.QMessageBox.warning(self, "Export failed", "Could not write the PDF export.")

    def _format_filter_values(self, values: Sequence[str], *, full: bool = False) -> str:
        display: List[str] = []
        for value in values:
            if value == MISSING_FILTER_TOKEN:
                display.append(self._missing_label)
            else:
                display.append(str(value))
        if full or len(display) <= 5:
            return ", ".join(display)
        return ", ".join(display[:5]) + f", … (+{len(display) - 5} more)"

    def _selected_molecule_filter(self) -> str:
        box = getattr(self, "molecule_filter_box", None)
        if not isinstance(box, QtWidgets.QComboBox):
            return ""
        data = box.currentData()
        molecule = str(data if data is not None else "").strip().lower()
        return molecule if molecule in {"plasmid", "chromosome"} else ""

    def _group_by_molecule_enabled(self) -> bool:
        checkbox = getattr(self, "group_molecule_check", None)
        return bool(isinstance(checkbox, QtWidgets.QCheckBox) and checkbox.isChecked())

    def _row_matches_filters(self, row: Dict[str, Any], query: str) -> bool:
        sample = str(row.get("sample") or "").strip()
        contig = str(row.get("contig") or "").strip()
        meta = row.get("sample_meta")
        meta_map = meta if isinstance(meta, dict) else {}

        molecule = self._selected_molecule_filter()
        if molecule:
            row_kind = _row_molecule_kind(row)
            if molecule == "plasmid" and row_kind != "plasmid":
                return False
            if molecule == "chromosome" and row_kind != "chromosome":
                return False

        if query:
            haystack: List[str] = [sample, contig]
            for k, v in meta_map.items():
                haystack.append(str(k))
                haystack.append(str(v))
            if not any(query in str(v).lower() for v in haystack if str(v)):
                return False

        for column, accepted in self._metadata_filters.items():
            value = str(meta_map.get(column, "")).strip()
            include_missing = MISSING_FILTER_TOKEN in accepted
            allowed = [v for v in accepted if v != MISSING_FILTER_TOKEN]
            if value:
                if allowed and value in allowed:
                    continue
                if include_missing and not value:
                    continue
                return False
            if include_missing:
                continue
            return False
        return True

    def _apply_filters(self) -> None:
        query = self.sample_query.text().strip().lower()
        visible = [row for row in self._all_rows if self._row_matches_filters(row, query)]
        visible = self._sort_rows(visible)
        self._visible_rows = list(visible)
        self._set_pan_core_sample_filter(
            [str(row.get("sample") or "") for row in visible]
        )
        display_rows = self._collapse_rows(visible) if self._collapsed else visible
        self.view.set_rows(display_rows)
        if self._collapsed:
            unique_samples = {
                str(row.get("sample") or "").strip()
                for row in visible
                if str(row.get("sample") or "").strip()
            }
            self.info.setText(
                f"{len(display_rows)} signature context(s) from {len(visible)}/{len(self._all_rows)} "
                f"track(s) across {len(unique_samples)} sample(s) | +/-15 genes around the selected gene"
            )
        else:
            self.info.setText(
                f"{len(visible)}/{len(self._all_rows)} context track(s) | +/-15 genes around the selected gene"
            )

    def _toggle_collapsed(self, checked: bool) -> None:
        self._collapsed = bool(checked)
        self.collapse_btn.setText("Uncollapse" if self._collapsed else "Collapse")
        self._apply_filters()

    def _reference_row(self) -> Optional[Dict[str, Any]]:
        for row in self._all_rows:
            if bool(row.get("is_source")):
                return row
        return self._all_rows[0] if self._all_rows else None

    @staticmethod
    def _context_signature(row: Dict[str, Any]) -> Tuple[Tuple[str, str, bool], ...]:
        genes = row.get("genes")
        if not isinstance(genes, list):
            return tuple()
        signature: List[Tuple[str, str, bool]] = []
        for gene in genes:
            if not isinstance(gene, dict):
                continue
            profile = _gene_homology_profile(gene)
            label = str(profile.get("key") or gene.get("label") or "").strip().lower()
            strand = str(gene.get("strand") or "+").strip() or "+"
            signature.append((label, strand, bool(gene.get("is_center"))))
        return tuple(signature)

    @staticmethod
    def _context_gene_profiles(row: Dict[str, Any]) -> List[Dict[str, Any]]:
        genes = row.get("genes")
        if not isinstance(genes, list):
            return []
        center_idx = 0
        for idx, gene in enumerate(genes):
            if isinstance(gene, dict) and bool(gene.get("is_center")):
                center_idx = idx
                break
        else:
            center_idx = max(0, len(genes) // 2)

        profiles: List[Dict[str, Any]] = []
        for idx, gene in enumerate(genes):
            if not isinstance(gene, dict):
                continue
            profile = _gene_homology_profile(gene)
            key = str(profile.get("key") or "").strip()
            if not key:
                fallback = _normalized_homology_text(gene.get("label"))
                if fallback:
                    key = f"raw:{fallback}"
            if not key:
                continue
            distance = abs(int(idx) - int(center_idx))
            local_weight = 1.0 + (1.35 / (1.0 + float(distance)))
            profiles.append(
                {
                    "idx": int(idx),
                    "center_distance": int(distance),
                    "key": key,
                    "strand": IgvCollectionDialog._normalized_context_strand(gene.get("strand")),
                    "weight": float(profile.get("weight") or 0.5) * local_weight,
                }
            )
        return profiles

    @staticmethod
    def _normalized_context_strand(value: Any) -> str:
        return "-" if str(value or "+").strip() == "-" else "+"

    @staticmethod
    def _weighted_collinear_score(
        reference: Sequence[Dict[str, Any]],
        candidate: Sequence[Dict[str, Any]],
    ) -> Tuple[float, int]:
        if not reference or not candidate:
            return 0.0, 0
        m = len(reference)
        n = len(candidate)
        dp: List[List[Tuple[float, int]]] = [[(0.0, 0) for _ in range(n + 1)] for _ in range(m + 1)]
        for i in range(m):
            ref = reference[i]
            for j in range(n):
                best = dp[i][j + 1] if dp[i][j + 1][0] >= dp[i + 1][j][0] else dp[i + 1][j]
                if str(ref.get("key") or "") == str(candidate[j].get("key") or ""):
                    strand_factor = (
                        1.0
                        if str(ref.get("strand") or "+") == str(candidate[j].get("strand") or "+")
                        else 0.82
                    )
                    center_delta = abs(int(ref.get("center_distance") or 0) - int(candidate[j].get("center_distance") or 0))
                    proximity_factor = max(0.70, 1.0 - (0.045 * float(center_delta)))
                    match_score = min(float(ref.get("weight") or 0.0), float(candidate[j].get("weight") or 0.0))
                    match_score *= strand_factor * proximity_factor
                    linked = (dp[i][j][0] + match_score, dp[i][j][1] + 1)
                    if linked[0] > best[0]:
                        best = linked
                dp[i + 1][j + 1] = best
        return dp[m][n]

    def _context_similarity_score(
        self,
        row: Dict[str, Any],
        reference_profiles: Optional[Sequence[Dict[str, Any]]] = None,
    ) -> float:
        ref_profiles = list(reference_profiles or [])
        if not ref_profiles:
            ref_row = self._reference_row()
            ref_profiles = self._context_gene_profiles(ref_row) if ref_row is not None else []
        cand_profiles = self._context_gene_profiles(row)
        if not ref_profiles and not cand_profiles:
            return 1.0
        if not ref_profiles or not cand_profiles:
            return 0.0

        forward_score, forward_count = self._weighted_collinear_score(ref_profiles, cand_profiles)
        reverse_candidate = [
            {
                **profile,
                "strand": "-" if str(profile.get("strand") or "+") == "+" else "+",
            }
            for profile in reversed(cand_profiles)
        ]
        reverse_score, reverse_count = self._weighted_collinear_score(ref_profiles, reverse_candidate)
        reverse_score *= 0.90
        raw_score = max(forward_score, reverse_score)
        match_count = max(forward_count, reverse_count)
        if match_count >= 3:
            raw_score *= 1.0 + min(0.28, 0.035 * float(match_count - 2))
        denom = max(
            sum(float(profile.get("weight") or 0.0) for profile in ref_profiles),
            sum(float(profile.get("weight") or 0.0) for profile in cand_profiles),
            1.0,
        )
        return max(0.0, min(1.0, raw_score / denom))

    @staticmethod
    def _sample_label(sample_count: int) -> str:
        return f"{sample_count} sample" if sample_count == 1 else f"{sample_count} samples"

    @staticmethod
    def _ordered_unique_strings(values: Sequence[str]) -> List[str]:
        out: List[str] = []
        seen: set[str] = set()
        for value in values:
            text = str(value).strip()
            if not text or text in seen:
                continue
            seen.add(text)
            out.append(text)
        return out

    def _collapsed_tooltip(
        self,
        *,
        sample_names: Sequence[str],
        track_titles: Sequence[str],
        sample_count: int,
        track_count: int,
    ) -> str:
        lines = [
            f"{self._sample_label(sample_count)} share this context signature",
            f"{track_count} context track{'s' if track_count != 1 else ''}",
        ]
        if sample_names:
            lines.extend(["", "Samples:"])
            lines.extend(sample_names[:20])
            if len(sample_names) > 20:
                lines.append(f"... (+{len(sample_names) - 20} more)")
        if track_count > sample_count and track_titles:
            lines.extend(["", "Tracks:"])
            lines.extend(track_titles[:20])
            if len(track_titles) > 20:
                lines.append(f"... (+{len(track_titles) - 20} more)")
        return "\n".join(lines)

    def _collapse_rows(self, rows: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
        if not rows:
            return []

        grouped: Dict[Tuple[Tuple[str, str, bool], ...], List[Dict[str, Any]]] = {}
        ordered_signatures: List[Tuple[Tuple[str, str, bool], ...]] = []
        for row in rows:
            signature = self._context_signature(row)
            if signature not in grouped:
                grouped[signature] = []
                ordered_signatures.append(signature)
            grouped[signature].append(row)

        collapsed_rows: List[Dict[str, Any]] = []
        for signature in ordered_signatures:
            members = grouped.get(signature, [])
            if not members:
                continue
            representative_row = min(
                members,
                key=lambda row: (
                    0 if bool(row.get("is_source")) else 1,
                    str(row.get("sample") or "").strip().lower(),
                    str(row.get("contig") or "").strip().lower(),
                ),
            )
            representative = dict(representative_row)
            sample_names = self._ordered_unique_strings(
                [str(member.get("sample") or "") for member in members]
            )
            track_titles = self._ordered_unique_strings(
                [
                    str(member.get("title") or f"{member.get('sample') or ''} | {member.get('contig') or ''}")
                    for member in members
                ]
            )
            sample_count = len(sample_names)
            track_count = len(members)
            molecule_kinds = {
                _row_molecule_kind(member)
                for member in members
                if isinstance(member, dict)
            }
            representative["display_label"] = self._sample_label(sample_count)
            if len(molecule_kinds) == 1:
                representative["molecule_kind"] = next(iter(molecule_kinds))
            elif molecule_kinds:
                representative["molecule_kind"] = "mixed"
            representative["meta_tooltip"] = self._collapsed_tooltip(
                sample_names=sample_names,
                track_titles=track_titles,
                sample_count=sample_count,
                track_count=track_count,
            )
            representative["collapsed_samples"] = list(sample_names)
            representative["collapsed_sample_count"] = sample_count
            representative["collapsed_track_count"] = track_count
            collapsed_rows.append(representative)

        return collapsed_rows

    def _sort_rows(self, rows: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
        if not rows:
            return []
        sort_col = "__sample__"
        ascending = True
        data = self.sort_column_box.currentData() if hasattr(self, "sort_column_box") else None
        if data is not None:
            sort_col = str(data)
        order_data = self.sort_order_box.currentData() if hasattr(self, "sort_order_box") else None
        if order_data is not None:
            ascending = str(order_data).strip().lower().startswith("asc")

        def _sample_key(row: Dict[str, Any]) -> str:
            return str(row.get("sample") or "").strip()

        def _contig_key(row: Dict[str, Any]) -> str:
            return str(row.get("contig") or "").strip()

        def _sort_value(row: Dict[str, Any]) -> str:
            if sort_col == "__sample__":
                return _sample_key(row)
            if sort_col == "__contig__":
                return _contig_key(row)
            meta = row.get("sample_meta")
            if isinstance(meta, dict):
                return str(meta.get(sort_col, "")).strip()
            return ""

        if sort_col == "__context_similarity__":
            ref_row = self._reference_row()
            ref_profiles = self._context_gene_profiles(ref_row) if ref_row is not None else []
            group_by_molecule = self._group_by_molecule_enabled()

            def _molecule_rank(row: Dict[str, Any]) -> Tuple[int, str]:
                if not group_by_molecule:
                    return (0, "")
                kind = _row_molecule_kind(row)
                if kind == "plasmid":
                    return (0, kind)
                if kind == "mixed":
                    return (1, kind)
                if kind == "chromosome":
                    return (2, kind)
                return (3, kind)

            scored_rows = [
                (self._context_similarity_score(row, ref_profiles), row)
                for row in rows
            ]
            if ascending:
                return [
                    row
                    for _score, row in sorted(
                        scored_rows,
                        key=lambda scored: (
                            _molecule_rank(scored[1]),
                            scored[0],
                            _sample_key(scored[1]).lower(),
                            _contig_key(scored[1]).lower(),
                        ),
                    )
                ]
            return [
                row
                for _score, row in sorted(
                    scored_rows,
                    key=lambda scored: (
                        _molecule_rank(scored[1]),
                        0 if bool(scored[1].get("is_source")) else 1,
                        -scored[0],
                        _sample_key(scored[1]).lower(),
                        _contig_key(scored[1]).lower(),
                    ),
                )
            ]

        present: List[Tuple[str, Dict[str, Any]]] = []
        missing: List[Dict[str, Any]] = []
        for row in rows:
            value = _sort_value(row)
            if value:
                present.append((value, row))
            else:
                missing.append(row)

        present_sorted = sorted(
            present,
            key=lambda vr: (
                str(vr[0]).lower(),
                _sample_key(vr[1]).lower(),
                _contig_key(vr[1]).lower(),
            ),
            reverse=not ascending,
        )
        missing_sorted = sorted(
            missing,
            key=lambda r: (_sample_key(r).lower(), _contig_key(r).lower()),
            reverse=not ascending,
        )
        return [row for _value, row in present_sorted] + missing_sorted
