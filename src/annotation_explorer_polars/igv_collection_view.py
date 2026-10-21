from __future__ import annotations

import hashlib
import html
import os
import re
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from PyQt6 import QtCore, QtGui, QtWidgets

from .constants import MISSING_FILTER_TOKEN, PLASMID_TRUE_VALUES
from .context_model import CATEGORY_LABELS, display_label_for_gene, optional_count
from .model.run_data import load_run_tag_databases, not_installed_key_note
from .metadata_filter_dialog import MetadataValueDialog
from .settings_column import (
    SETTINGS_CONTROL_HEIGHT,
    SETTINGS_MARGIN,
    SETTINGS_STRIP_WIDTH,
    ElidedHint,
    SettingsColumn,
    SettingsStrip,
    WrapLayout,
    auto_fold,
    make_chip,
    select_combo_data,
    settings_card,
    settings_combo,
    settings_control,
    stacked_field,
)


def _scalar_text(value: Any) -> str:
    return "" if value is None else str(value).strip()


def _truthy_text(value: Any) -> bool:
    text = _scalar_text(value).lower()
    return text in PLASMID_TRUE_VALUES


def _row_molecule_kind(row: Dict[str, Any]) -> str:
    explicit = _scalar_text(row.get("molecule")).lower()
    if explicit not in {"plasmid", "chromosome", "mixed", "unknown"}:
        explicit = _scalar_text(row.get("molecule_kind")).lower()
    if explicit in {"plasmid", "chromosome", "mixed", "unknown"}:
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
    return "unknown"


def _molecule_badge(row: Dict[str, Any]) -> Tuple[str, str]:
    kind = _row_molecule_kind(row)
    if kind == "plasmid":
        return "P", "#16a34a"
    if kind == "mixed":
        return "P/C", "#64748b"
    return "C", "#2563eb"


_EMPTY_TEXT_VALUES = {"", "-", ".", "na", "n/a", "nan", "none", "null", "<na>"}


def _contig_text(value: Any) -> str:
    """'contig 48' for a bare ID; leave IDs that already say 'contig' alone."""
    text = str(value if value is not None else "").strip()
    if not text or text.lower().startswith("contig"):
        return text
    return f"contig {text}"


def _strand_text(value: Any) -> str:
    """Name a strand in words; '-' is a real strand here, not a missing value."""
    text = str(value if value is not None else "").strip()
    if text == "-":
        return "minus"
    if text == "+":
        return "plus"
    return "unknown"
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
        # Compact is the default because comparative context views can contain
        # dozens of rows.  Comfortable mode remains available from the canvas
        # toolbar when longer rotated gene labels are more important than
        # seeing many neighborhoods at once.
        self._density_mode = "compact"
        self._row_gap = 62.0
        self._top_margin = 70.0
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
        self._left_align_contexts = False
        self._show_molecule_badges = True
        self._static_row_labels = True
        self._normalize_orientation = True
        # Pixel lane before the first gene track. It is converted to scene units
        # on render so zooming cannot collapse labels into the track.
        self._track_x_offset = self._row_label_min_lane_px
        self._computed_track_x_offset = self._track_x_offset
        self._focal_x_scene = 0.0
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
        self._focal_x_scene = 0.0
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
        comparison_bits: List[str] = []
        if bool(row.get("is_source") or row.get("is_comparison_reference")):
            comparison_bits.append("REFERENCE")
        else:
            try:
                comparison_bits.append(f"{float(row.get('context_similarity')):.0%} similar")
            except Exception:
                pass
            orientation = str(row.get("context_orientation") or "").strip().lower()
            if orientation == "reverse":
                comparison_bits.append("reverse match")
            elif orientation == "forward":
                comparison_bits.append("forward match")
        _anchor, flipped = self._row_anchor_and_flip(row)
        if flipped:
            comparison_bits.append("flipped")
        if comparison_bits:
            title = " · ".join(comparison_bits + [str(title)])
        badge_text, _badge_color = _molecule_badge(row)
        return str(title), str(badge_text or "")

    @staticmethod
    def _comparison_tooltip(row: Dict[str, Any]) -> str:
        if bool(row.get("is_source") or row.get("is_comparison_reference")):
            return "Reference neighborhood used for similarity ranking"
        lines: List[str] = []
        try:
            lines.append(f"Similarity to reference: {float(row.get('context_similarity')):.1%}")
        except Exception:
            pass
        orientation = str(row.get("context_orientation") or "").strip().lower()
        if orientation:
            lines.append(f"Best collinear orientation: {orientation}")
        try:
            lines.append(f"Matched genes: {int(row.get('context_match_count'))}")
        except Exception:
            pass
        return "\n".join(lines)

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

    def set_density_mode(self, mode: str) -> None:
        normalized = str(mode or "").strip().lower()
        if normalized not in {"compact", "comfortable"}:
            normalized = "compact"
        if normalized == self._density_mode:
            return
        self._density_mode = normalized
        if normalized == "comfortable":
            self._row_gap = 108.0
            self._top_margin = 82.0
        else:
            self._row_gap = 62.0
            self._top_margin = 70.0
        self._render_scene()
        self._center_on_origin()

    def set_normalize_orientation(self, enabled: bool) -> None:
        normalized = bool(enabled)
        if normalized == self._normalize_orientation:
            return
        self._normalize_orientation = normalized
        self._render_scene()
        self._center_on_origin()

    def set_rows(self, rows: Sequence[Dict[str, Any]]) -> None:
        self._emit_context_hover(None)
        self._rows = [dict(r) for r in rows if isinstance(r, dict)]
        self.resetTransform()
        self._render_scene()
        QtCore.QTimer.singleShot(0, self._center_on_origin)

    def _center_on_origin(self) -> None:
        if self._scene is None or self._scene.itemsBoundingRect().isNull():
            return
        if self._focal_x_scene > 0.0:
            scale_x = self._scale_x()
            viewport_width = max(1.0, float(self.viewport().width()))
            label_lane = (
                min(float(self._computed_track_x_offset), viewport_width * 0.45)
                if self._static_row_labels
                else 0.0
            )
            target_view_x = label_lane + max(1.0, viewport_width - label_lane) / 2.0
            hbar = self.horizontalScrollBar()
            target = int(round((self._focal_x_scene * scale_x) - target_view_x))
            hbar.setValue(max(hbar.minimum(), min(hbar.maximum(), target)))
            self.verticalScrollBar().setValue(self.verticalScrollBar().minimum())
            self._reposition_row_labels()
            self._reposition_figure_legend()
            return
        if self._left_align_contexts:
            self.horizontalScrollBar().setValue(self.horizontalScrollBar().minimum())
            self._reposition_row_labels()
            return
        self.centerOn(self._scene.sceneRect().center())

    @staticmethod
    def _gene_bounds(gene: Dict[str, Any]) -> Optional[Tuple[float, float]]:
        try:
            start = float(gene.get("start") or 0)
            end = float(gene.get("end") or 0)
        except Exception:
            return None
        return min(start, end), max(start, end)

    def _focal_gene_for_row(self, row: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        genes = row.get("genes")
        if not isinstance(genes, list):
            return None
        valid_genes = [gene for gene in genes if isinstance(gene, dict) and self._gene_bounds(gene)]
        if not valid_genes:
            return None
        for gene in valid_genes:
            if bool(gene.get("is_center")):
                return gene
        try:
            center_bp = float(row.get("center_bp"))
        except Exception:
            center_bp = None  # type: ignore[assignment]
        if center_bp is not None:
            containing = [
                gene
                for gene in valid_genes
                if (bounds := self._gene_bounds(gene)) is not None
                and bounds[0] <= center_bp <= bounds[1]
            ]
            if containing:
                return containing[0]
            return min(
                valid_genes,
                key=lambda gene: abs(
                    sum(self._gene_bounds(gene) or (0.0, 0.0)) / 2.0 - center_bp
                ),
            )
        return valid_genes[len(valid_genes) // 2]

    def _row_anchor_and_flip(self, row: Dict[str, Any]) -> Tuple[float, bool]:
        focal_gene = self._focal_gene_for_row(row)
        focal_bounds = self._gene_bounds(focal_gene) if focal_gene is not None else None
        try:
            anchor = float(row.get("center_bp"))
        except Exception:
            anchor = None  # type: ignore[assignment]
        if anchor is None and focal_bounds is not None:
            anchor = (focal_bounds[0] + focal_bounds[1]) / 2.0
        if anchor is None:
            genes = row.get("genes")
            spans: List[Tuple[float, float]] = []
            if isinstance(genes, list):
                spans = [
                    bounds
                    for gene in genes
                    if isinstance(gene, dict)
                    if (bounds := self._gene_bounds(gene)) is not None
                ]
            anchor = min((start for start, _end in spans), default=0.0)
        focal_strand = self._normalized_strand(focal_gene.get("strand")) if focal_gene else "+"
        return float(anchor), bool(self._normalize_orientation and focal_strand == "-")

    @staticmethod
    def _relative_bounds(start: float, end: float, anchor: float, flip: bool) -> Tuple[float, float]:
        rel_start = min(start, end) - anchor
        rel_end = max(start, end) - anchor
        if flip:
            return -rel_end, -rel_start
        return rel_start, rel_end

    def _compute_row_anchors(self) -> Tuple[List[Tuple[float, bool]], float, float]:
        row_layouts: List[Tuple[float, bool]] = []
        min_rel = None
        max_rel = None
        for row in self._rows:
            genes = row.get("genes") or []
            if not isinstance(genes, list):
                row_layouts.append((0.0, False))
                continue
            anchor, flip = self._row_anchor_and_flip(row)
            row_layouts.append((anchor, flip))
            for gene in genes:
                if not isinstance(gene, dict):
                    continue
                bounds = self._gene_bounds(gene)
                if bounds is None:
                    continue
                rel_start, rel_end = self._relative_bounds(bounds[0], bounds[1], anchor, flip)
                min_rel = rel_start if min_rel is None else min(min_rel, rel_start)
                max_rel = rel_end if max_rel is None else max(max_rel, rel_end)

        if min_rel is None or max_rel is None:
            return row_layouts, -1000.0, 1000.0
        if max_rel <= min_rel:
            max_rel = min_rel + 1.0
        return row_layouts, float(min_rel), float(max_rel)

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
        self._focal_x_scene = 0.0
        if not self._rows:
            return
        selected_before = self._selected_gene_key

        row_layouts, min_rel, max_rel = self._compute_row_anchors()
        self._min_rel_bp = min_rel
        self._max_rel_bp = max_rel
        self._computed_track_x_offset = self._compute_track_lane_px()
        track_offset = self._static_track_offset_scene()
        focal_x = track_offset + max(0.0, -self._min_rel_bp) + 18.0
        self._focal_x_scene = float(focal_x)
        x0 = track_offset
        x1 = focal_x + self._max_rel_bp + 36.0

        for row_idx, row in enumerate(self._rows):
            y_top = self._top_margin + row_idx * (self._row_height + self._row_gap)
            y_center = y_top + self._row_height / 2.0
            title, _row_badge_text = self._row_label_text_for_row(row, row_idx)
            tooltip = str(row.get("meta_tooltip") or "").strip() or title
            comparison_tooltip = self._comparison_tooltip(row)
            if comparison_tooltip:
                tooltip = f"{tooltip}\n\n{comparison_tooltip}"

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
            row_anchor, flip_row = row_layouts[row_idx] if row_idx < len(row_layouts) else (0.0, False)
            for gene_idx, gene in enumerate(genes):
                if not isinstance(gene, dict):
                    continue
                try:
                    start = int(gene.get("start") or 0)
                    end = int(gene.get("end") or 0)
                except Exception:
                    continue
                strand = self._normalized_strand(gene.get("strand"))
                label = str(gene.get("label") or "gene")
                color_value = gene.get("color")
                color = color_value if isinstance(color_value, QtGui.QColor) else QtGui.QColor(str(color_value or "#90caf9"))
                raw_rel_start, raw_rel_end = self._relative_bounds(
                    float(start),
                    float(end),
                    float(row_anchor),
                    bool(flip_row),
                )
                rel_start = focal_x + raw_rel_start
                rel_end = focal_x + raw_rel_end
                if flip_row:
                    strand = "-" if strand == "+" else "+"
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

        # Keep the focal alignment explicit in interactive use and exports.
        last_row_y = self._top_margin + (len(self._rows) - 1) * (self._row_height + self._row_gap)
        focal_guide = QtWidgets.QGraphicsLineItem(
            focal_x,
            max(18.0, self._top_margin - 38.0),
            focal_x,
            last_row_y + self._row_height + 8.0,
        )
        focal_pen = QtGui.QPen(QtGui.QColor("#f59e0b"), 1.0, QtCore.Qt.PenStyle.DashLine)
        focal_guide.setPen(focal_pen)
        focal_guide.setZValue(3)
        focal_guide.setToolTip("Focal gene alignment")
        self._scene.addItem(focal_guide)

        focal_label = QtWidgets.QGraphicsSimpleTextItem("focal gene")
        focal_font = QtGui.QFont()
        focal_font.setPointSize(8)
        focal_font.setBold(True)
        focal_label.setFont(focal_font)
        focal_label.setBrush(QtGui.QBrush(QtGui.QColor("#b45309")))
        focal_label.setPos(focal_x + 6.0, max(2.0, self._top_margin - 52.0))
        focal_label.setZValue(39)
        self._scene.addItem(focal_label)

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
            if payload is not None:
                text_item.setData(0, payload)
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
        if event.button() == QtCore.Qt.MouseButton.LeftButton:
            payload = self._payload_from_item(item)
            if payload is not None:
                fallback = str(payload.get("Gene") or payload.get("Product") or "")
                gene_key = self._gene_key(payload, fallback)
                if gene_key:
                    # A second click on the pinned homologue family clears it;
                    # double-click still proceeds to the full annotation card.
                    self._apply_gene_selection(
                        None if gene_key == self._selected_gene_key else gene_key
                    )
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


# Keep the established import and dialog contract while the new renderer rolls
# out.  Setting BACTOLOMEW_CONTEXT_RENDERER=classic is an intentionally
# temporary diagnostic escape hatch; the high-quality canvas is the default.
LegacyIgvCollectionView = IgvCollectionView
if os.environ.get("BACTOLOMEW_CONTEXT_RENDERER", "v2").strip().lower() == "classic":
    IgvCollectionView = LegacyIgvCollectionView
else:
    from .context_canvas import ContextCanvas as IgvCollectionView


# The annotation tags of the key's gene types (its colors).
KEY_GENE_TYPE_TAGS = tuple(
    (tag, CATEGORY_LABELS[category])
    for tag, category in (("A", "AMR"), ("V", "VIRULENCE"), ("M", "MOBILE"), ("B", "BMRG"))
)


# How far from the focal gene a mobile gene counts as flanking it.
MOBILE_FLANK_GENES = 5
# Transposases, IS elements, transposons and integrases.  IS and Tn names
# are matched case-sensitively so words such as 'is2' never count.
_MOBILE_NAME_PATTERN = re.compile(
    r"\b((?i:tnp\w*|transpos\w*|integrase)|IS\d+\w*|IS[A-Z][A-Za-z]*\d+\w*|Tn\d+\w*|intI\d*)\b"
)
DETAIL_FIELD_NAMES = ("Gene", "Product", "Position", "Strand", "Gene type", "Found in")
CONTEXT_FOOTER_HINT = (
    "Hover a gene to trace it. Click to pin. Arrow keys follow the family. "
    "Double-click or Enter for the full gene card. Ctrl+scroll to zoom."
)


class NeighbourhoodView(QtWidgets.QWidget):
    """The gene neighborhood (and plasmid context) comparison: the title line,
    the Settings column, the plot and its pinned gene details.

    It is the body of the ``NeighbourhoodPanel`` of an explore window
    (``panel_mode=True``): there the rows follow the window's one sample set
    (``set_sample_scope``) instead of the viewer's own metadata filters (D4),
    and its two shortcuts act only inside the panel (plan 3.7). Without panel
    mode it is a standalone viewer with its own metadata filters and
    window-wide shortcuts (the tests build it that way since the
    ``IgvCollectionDialog`` window shim went in Phase 10).
    """

    geneDoubleClicked = QtCore.pyqtSignal(dict)
    # Ctrl+W: a panel asks its window to close it (a standalone viewer's
    # owner decides).
    closeRequested = QtCore.pyqtSignal()
    # None until the viewer opens or closes the Settings panel; that choice
    # is remembered for the rest of the session (never written to QSettings)
    # on the ``settings_memory`` class the opener passes (the panel's class;
    # a standalone viewer uses its own). Until then the panel starts open in windows wide
    # enough to keep a usable plot beside it.
    settings_panel_open: Optional[bool] = None
    SETTINGS_PANEL_AUTO_MIN_WIDTH = 1400
    # Panel mode: the plot's minimum height, and the room the plot keeps
    # beside an open Settings column; a panel smaller than that folds the
    # column while it is small and opens it again when it grows.
    PANEL_PLOT_MIN_HEIGHT = 140
    OUT_OF_SCOPE_TEXT = "None of these neighborhoods are in this window's samples.\nChange the filters in the rail to show them."
    PANEL_PLOT_MIN_WIDTH = 360
    # One height for every control in the Settings column, one gap between
    # its cards and one margin around them.
    SETTINGS_CONTROL_HEIGHT = SETTINGS_CONTROL_HEIGHT
    SETTINGS_MARGIN = SETTINGS_MARGIN
    SETTINGS_STRIP_WIDTH = SETTINGS_STRIP_WIDTH
    PIN_CHIP_MAX_TEXT_WIDTH = 220

    def __init__(
        self,
        parent: Optional[QtWidgets.QWidget],
        *,
        rows: Sequence[Dict[str, Any]],
        title: str = "IGV_collection",
        output_dir: Optional[Path] = None,
        flank_sizes: Optional[Sequence[int]] = None,
        initial_flank: Optional[int] = None,
        open_samples_callback: Optional[Callable[[List[str]], None]] = None,
        panel_mode: bool = False,
        settings_memory: Optional[type] = None,
    ) -> None:
        super().__init__(parent)
        self._title = str(title or "")
        self._panel_mode = bool(panel_mode)
        self._settings_memory = settings_memory if settings_memory is not None else type(self)
        self._all_rows: List[Dict[str, Any]] = [dict(r) for r in rows if isinstance(r, dict)]
        # Panel mode: the samples of the window's sample set (None: every row).
        self._scope_samples: Optional[frozenset] = None
        self._scope_notes: List[str] = []
        self._linked_samples: Tuple[str, ...] = ()
        # Shows samples in the main window's table; the row menu offers
        # 'Open sample in table' only when the opener passes one.
        self._open_samples_callback = open_samples_callback
        # Set by 'Use as reference'; None keeps the default reference.
        self._reference_identity: Optional[Tuple[str, str, str]] = None
        self._molecule_filter = ""
        self._metadata_filters: Dict[str, List[str]] = {}
        self._metadata_columns = self._collect_metadata_columns()
        self._missing_label = "(missing)"
        self._collapsed = False
        self._context_color_mode = "annotation"
        self._context_label_mode = "reference"
        # Figure option: draw contigs with only the focal gene as thin rows.
        self._export_include_short_contigs = True
        self._has_focal_gene = any(
            isinstance(gene, dict) and bool(gene.get("is_center"))
            for row in self._all_rows
            for gene in (row.get("genes") or [])
        )
        self._output_dir = Path(output_dir) if output_dir is not None else None
        self._visible_rows: List[Dict[str, Any]] = []
        self._filter_timer = QtCore.QTimer(self)
        self._filter_timer.setSingleShot(True)
        self._filter_timer.setInterval(170)
        self._filter_timer.timeout.connect(self._apply_filters)
        self._flank_sizes = sorted(
            {
                int(value)
                for value in (flank_sizes or [])
                if isinstance(value, (int, float)) and int(value) > 0
            }
        )
        requested_flank = int(initial_flank) if initial_flank is not None else None
        if self._flank_sizes:
            self._context_flank = (
                requested_flank
                if requested_flank in self._flank_sizes
                else self._flank_sizes[min(len(self._flank_sizes) - 1, 2)]
            )
        else:
            self._context_flank = None

        # Panel mode: the Settings column folded because the panel is small
        # (it opens again when the panel grows; the user's choice stands).
        self._settings_auto_folded = False
        self._settings_fit_timer = QtCore.QTimer(self)
        self._settings_fit_timer.setSingleShot(True)
        self._settings_fit_timer.setInterval(0)
        self._settings_fit_timer.timeout.connect(self._fit_settings_to_panel)

        self._build_ui()
        self._render_metadata_filters()
        self._apply_filters()

    def title_text(self) -> str:
        """The comparison's title, e.g. 'Gene neighborhoods: blaOXA-23'."""
        return self._title

    def is_panel_mode(self) -> bool:
        return self._panel_mode

    def _build_ui(self) -> None:
        root = QtWidgets.QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)
        # The title line spans the window, so opening the Settings column
        # below it never shoves it sideways.
        root.addWidget(self._build_title_block())

        body = QtWidgets.QWidget(self)
        body_layout = QtWidgets.QHBoxLayout(body)
        body_layout.setContentsMargins(0, 0, 0, 0)
        body_layout.setSpacing(0)
        # The folded column: a slim strip with the expand button.
        body_layout.addWidget(self._build_settings_strip(body))
        main_split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal, body)
        # Never dragged shut: the header and strip buttons are the one way to
        # fold the column, so the strip shows exactly when the column is shut.
        main_split.setChildrenCollapsible(False)
        self.main_split = main_split
        body_layout.addWidget(main_split, stretch=1)
        root.addWidget(body, stretch=1)

        left = self._build_rows_panel()
        main_split.addWidget(left)

        right = QtWidgets.QWidget(self)
        right_layout = QtWidgets.QVBoxLayout(right)
        right_layout.setContentsMargins(self.SETTINGS_MARGIN, 8, self.SETTINGS_MARGIN, 4)
        right_layout.setSpacing(0)
        view_panel = QtWidgets.QWidget(right)
        view_panel_layout = QtWidgets.QVBoxLayout(view_panel)
        view_panel_layout.setContentsMargins(0, 0, 0, 0)
        view_panel_layout.setSpacing(6)

        view_split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal, view_panel)
        view_split.setChildrenCollapsible(True)
        self.view_split = view_split
        self.view = IgvCollectionView(view_split)
        if self._panel_mode:
            # A docked panel can sit in a short group: the plot scrolls
            # instead of pushing the footer out of the panel.
            self.view.setMinimumHeight(self.PANEL_PLOT_MIN_HEIGHT)
        self.view.geneDoubleClicked.connect(self.geneDoubleClicked.emit)
        if hasattr(self.view, "set_not_installed_note"):
            self.view.set_not_installed_note(
                *not_installed_key_note(load_run_tag_databases(self._output_dir), KEY_GENE_TYPE_TAGS)
            )
        self.view.syntenyBlockClicked.connect(self._show_synteny_block_details)
        if hasattr(self.view, "geneSelected"):
            self.view.geneSelected.connect(self._show_gene_details)
        if hasattr(self.view, "rowMenuRequested"):
            self.view.rowMenuRequested.connect(self._show_row_menu)
        self.normalize_orientation_check.toggled.connect(self.view.set_normalize_orientation)
        self.fit_reference_btn.clicked.connect(self.view.reset_zoom)
        self.show_all_btn.clicked.connect(self.view.fit_all_contexts)
        # The Pinned chip sits in the plot's header strip; the classic
        # renderer has none, so there it heads the Actions card.
        attach_chip = getattr(self.view, "set_header_corner_widget", None)
        if callable(attach_chip):
            attach_chip(self.pin_chip)
        else:
            self._actions_layout.insertWidget(1, self.pin_chip, alignment=QtCore.Qt.AlignmentFlag.AlignLeft)
        self._set_context_view_option("set_export_title", self.title_text())
        self._set_context_view_option("set_color_mode", self._context_color_mode)
        self._set_context_view_option("set_label_mode", self._context_label_mode)
        self._set_context_view_option("set_export_include_short_contigs", self._export_include_short_contigs)
        self._hide_controls_the_view_lacks()
        view_split.addWidget(self.view)
        view_split.addWidget(self._build_detail_panel(view_split))
        self.synteny_detail_panel.hide()
        view_split.setStretchFactor(0, 1)
        view_split.setStretchFactor(1, 0)
        view_split.setSizes([980, 300])
        view_panel_layout.addWidget(view_split, stretch=1)

        right_layout.addWidget(view_panel, stretch=1)
        main_split.addWidget(right)

        main_split.setStretchFactor(0, 0)
        main_split.setStretchFactor(1, 1)
        main_split.setSizes([self._filter_rail_last_width, 1480])

        # The Settings column holds every control; its header button folds
        # it to the strip for a wider plot.
        is_open = self._settings_memory.settings_panel_open is not False
        left.setVisible(is_open)
        self.settings_strip.setVisible(not is_open)
        footer = QtWidgets.QFrame(self)
        footer.setObjectName("ContextFooter")
        footer_layout = QtWidgets.QHBoxLayout(footer)
        footer_layout.setContentsMargins(16, 5, 16, 6)
        if self._panel_mode:
            self.footer_hint = ElidedHint(CONTEXT_FOOTER_HINT, footer)
            self.footer_hint.setToolTip(CONTEXT_FOOTER_HINT)
        else:
            self.footer_hint = QtWidgets.QLabel(CONTEXT_FOOTER_HINT, footer)
        self.footer_hint.setProperty("role", "contextHint")
        footer_layout.addWidget(self.footer_hint, stretch=1)
        root.addWidget(footer)
        # In a panel both act only while the focus is inside it (plan 3.7):
        # Ctrl+F in another panel never jumps here, and Ctrl+W closes only
        # this panel. A standalone viewer keeps window-wide shortcuts.
        context = (
            QtCore.Qt.ShortcutContext.WidgetWithChildrenShortcut
            if self._panel_mode
            else QtCore.Qt.ShortcutContext.WindowShortcut
        )
        for keys, slot in (("Ctrl+W", self.closeRequested.emit), ("Ctrl+F", self._focus_context_search)):
            shortcut = QtGui.QShortcut(QtGui.QKeySequence(keys), self)
            shortcut.setContext(context)
            shortcut.activated.connect(slot)

    def _build_settings_strip(self, parent: QtWidgets.QWidget) -> QtWidgets.QWidget:
        """The folded Settings column: an expand button and a filter count
        (active filters stay visible while the column is folded)."""
        strip = SettingsStrip(parent, "Show settings (Ctrl+F)")
        strip.expandRequested.connect(lambda: self._toggle_filter_rail(True))
        self.settings_strip = strip
        self.settings_expand_btn = strip.expand_button
        self.settings_strip_badge = strip.badge
        self.settings_strip_caption = strip.caption
        return strip

    def _build_rows_panel(self) -> QtWidgets.QWidget:
        """The Settings column: a header with the fold button, then the
        'View' (looks and which rows), 'Order' and 'Actions' cards."""
        left = SettingsColumn(self, object_name="FilterRail", max_width=420)
        left.collapseRequested.connect(lambda: self._toggle_filter_rail(False))
        self.filter_rail = left
        self.settings_title = left.title
        self.settings_collapse_btn = left.collapse_button
        self.settings_scroll = left.scroll
        self.settings_footer = left.footer
        self.settings_filter_badge = make_chip("", "Rows are filtered: Clear in 'View' resets them", "primary")
        left.insert_header_widget(self.settings_filter_badge)
        self.settings_filter_badge.hide()
        content = left.content
        left.add_card(self._build_view_card(content))
        left.add_card(self._build_order_card(content))
        # Actions stay pinned under the scrolling cards, so zoom and Export
        # are in view at any window height.
        left.footer_layout.addWidget(self._build_actions_card(left.footer))
        left.finish(min_width=248)
        self._filter_rail_last_width = max(300, left.minimumWidth())
        return left

    def _build_view_card(self, parent: QtWidgets.QWidget) -> QtWidgets.QFrame:
        """Flank, color, labels, mirroring, the plasmid/chromosome split,
        merging, and which rows are shown."""
        card, layout = settings_card(parent, "View")
        self.view_card = card

        if self._flank_sizes:
            flank_options = [
                (int(flank), f"±{flank}", f"Show up to {flank} genes on each side of the focal gene")
                for flank in self._flank_sizes
            ]
        else:
            flank_options = [("full", "Full contig", "Each row shows every loaded gene of the plasmid")]
        self.flank_box = settings_combo(card, flank_options, self._context_flank, self._set_context_flank)
        if not self._flank_sizes:
            self.flank_box.setEnabled(False)
        self.flank_group = stacked_field(card, layout, "Flank", self.flank_box)

        self.color_by_box = settings_combo(
            card,
            [
                ("annotation", "Gene type", "Color genes by type: AMR, virulence, mobile, BMRG or other"),
                ("family", "Gene family", "Same color means the same gene family (colors can repeat)"),
            ],
            self._context_color_mode,
            self._set_context_color_mode,
        )
        self.color_by_group = stacked_field(card, layout, "Color by", self.color_by_box)

        self.labels_box = settings_combo(
            card,
            [
                (
                    "reference",
                    "Key genes",
                    "Name the reference row, plus AMR, mobile and virulence genes on every row once zoomed in",
                ),
                ("all", "All", "Name every gene that fits"),
                ("off", "Off", "Name only the gene under the pointer"),
            ],
            self._context_label_mode,
            self._set_context_label_mode,
        )
        self.labels_group = stacked_field(card, layout, "Labels", self.labels_box)

        self.normalize_orientation_check = QtWidgets.QCheckBox("Mirror rows to match reference", card)
        self.normalize_orientation_check.setChecked(True)
        self.normalize_orientation_check.setToolTip(
            "Flip other rows so their matching genes run the same way as the reference; "
            "the reference is never flipped"
        )
        layout.addWidget(self.normalize_orientation_check)
        self.group_molecule_check = QtWidgets.QCheckBox("Split plasmid and chromosome", card)
        # Off by default so the list is one ranking: reference first, then
        # match to the reference, highest first.
        self.group_molecule_check.setChecked(False)
        self.group_molecule_check.setToolTip(
            "Split the rows into plasmid and chromosome blocks, plasmids first, "
            "each ordered on its own under a labelled divider"
        )
        self.group_molecule_check.toggled.connect(self._apply_filters)
        layout.addWidget(self.group_molecule_check)
        # A checkbox like the split above, so the on state shows at a glance.
        self.collapse_btn = QtWidgets.QCheckBox("Merge identical neighborhoods", card)
        self.collapse_btn.setToolTip(
            "Show neighborhoods with exactly the same genes, order and strands as one row"
        )
        self.collapse_btn.toggled.connect(self._toggle_collapsed)
        layout.addWidget(self.collapse_btn)

        # Which rows: search, molecule switch and metadata filters (in a
        # panel: the window's sample set instead of the metadata filters, D4).
        self.sample_query = QtWidgets.QLineEdit(card)
        self.sample_query.setPlaceholderText("Search samples and metadata...")
        self.sample_query.setClearButtonEnabled(True)
        self.sample_query.textChanged.connect(self._schedule_filter_apply)
        layout.addWidget(settings_control(self.sample_query))
        if self._panel_mode:
            # A panel has no search of its own: the window's search (the
            # filter rail) sets which samples it shows. Kept, empty and
            # hidden, so the filter code stays one.
            self.sample_query.setVisible(False)
        self.molecule_box = settings_combo(
            card,
            [
                ("", "All", "Show every contig"),
                ("plasmid", "Plasmid", "Show plasmid contigs only"),
                ("chromosome", "Chromosome", "Show chromosome contigs only"),
            ],
            self._molecule_filter,
            self._set_molecule_filter,
        )
        self.molecule_group = stacked_field(card, layout, "Molecule", self.molecule_box)

        self.add_filter_btn = QtWidgets.QPushButton("Add filter", card)
        self.add_filter_btn.setEnabled(bool(self._metadata_columns) and not self._panel_mode)
        # A panel has no metadata filters of its own: the rail's apply to
        # every panel of the window (D4).
        self.add_filter_btn.setVisible(not self._panel_mode)
        self.add_filter_btn.setToolTip(
            "Filter rows by sample metadata" if self._metadata_columns else "No sample metadata to filter by"
        )
        self.add_filter_btn.clicked.connect(self._show_add_filter_menu)
        settings_control(self.add_filter_btn)
        self.clear_filters_btn = QtWidgets.QToolButton(card)
        self.clear_filters_btn.setObjectName("ContextLinkButton")
        self.clear_filters_btn.setText("Clear")
        self.clear_filters_btn.setToolTip(
            "Clear the search and the molecule switch" if self._panel_mode else "Clear the search and every filter"
        )
        self.clear_filters_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.clear_filters_btn.clicked.connect(self._clear_context_filters)
        # Just as wide as its text, so the link ends flush with the controls.
        self.clear_filters_btn.ensurePolished()
        self.clear_filters_btn.setFixedWidth(self.clear_filters_btn.fontMetrics().horizontalAdvance("Clear") + 4)
        filter_row = QtWidgets.QHBoxLayout()
        filter_row.setContentsMargins(0, 0, 0, 0)
        filter_row.setSpacing(8)
        filter_row.addWidget(self.add_filter_btn)
        filter_row.addStretch(1)
        filter_row.addWidget(self.clear_filters_btn, alignment=QtCore.Qt.AlignmentFlag.AlignVCenter)
        layout.addLayout(filter_row)

        # Filter chips, one per line under the button.
        self.filters_container = QtWidgets.QWidget(card)
        self.filters_container.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Preferred,
            QtWidgets.QSizePolicy.Policy.Maximum,
        )
        self.filters_layout = QtWidgets.QVBoxLayout(self.filters_container)
        self.filters_layout.setContentsMargins(0, 0, 0, 0)
        self.filters_layout.setSpacing(6)
        self.filters_layout.setAlignment(QtCore.Qt.AlignmentFlag.AlignTop)
        # Hidden while empty, so no blank band sits under the button row.
        self.filters_container.hide()
        layout.addWidget(self.filters_container)
        self.scope_label: Optional[QtWidgets.QLabel] = None
        self._scope_chip_layout: Optional[WrapLayout] = None
        if self._panel_mode:
            # Read only: why the list may be shorter than the comparison.
            self.scope_label = QtWidgets.QLabel("", card)
            self.scope_label.setObjectName("ContextScopeLabel")
            self.scope_label.setProperty("role", "muted")
            self.scope_label.setWordWrap(True)
            layout.addWidget(self.scope_label)
            self.scope_chips = QtWidgets.QWidget(card)
            self._scope_chip_layout = WrapLayout(self.scope_chips, spacing=6, line_spacing=4)
            self._scope_chip_layout.setContentsMargins(0, 0, 0, 0)
            self.scope_chips.hide()
            layout.addWidget(self.scope_chips)
            self._render_scope()
        return card

    def _build_order_card(self, parent: QtWidgets.QWidget) -> QtWidgets.QFrame:
        """Sort key and direction."""
        card, layout = settings_card(parent, "Order")
        self.order_card = card
        self.sort_column_box = QtWidgets.QComboBox(card)
        self.sort_column_box.addItem("Match to reference", "__context_similarity__")
        self.sort_column_box.addItem("Sample", "__sample__")
        self.sort_column_box.addItem("Contig", "__contig__")
        for col in self._metadata_columns:
            self.sort_column_box.addItem(col, col)
        self.sort_order_box = QtWidgets.QComboBox(card)
        self.sort_order_box.addItem("Ascending", "asc")
        self.sort_order_box.addItem("Descending", "desc")
        sort_order_idx = self.sort_order_box.findData("desc")
        if sort_order_idx >= 0:
            self.sort_order_box.setCurrentIndex(sort_order_idx)
        self.sort_column_box.currentIndexChanged.connect(self._on_sort_column_changed)
        self.sort_order_box.currentIndexChanged.connect(self._apply_filters)
        for box in (self.sort_column_box, self.sort_order_box):
            # Long metadata names never widen the column.
            box.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
            box.setMinimumContentsLength(12)
            settings_control(box)
        self.sort_column_group = stacked_field(card, layout, "Sort by", self.sort_column_box)
        self.sort_order_group = stacked_field(card, layout, "Order", self.sort_order_box)
        return card

    def _build_actions_card(self, parent: QtWidgets.QWidget) -> QtWidgets.QFrame:
        """Zoom and figure export."""
        card, layout = settings_card(parent, "Actions")
        self.actions_card = card
        self._actions_layout = layout

        # Shown in the plot's header strip once the view exists.
        self.pin_chip = QtWidgets.QFrame(card)
        self.pin_chip.setObjectName("ContextPinChip")
        pin_layout = QtWidgets.QHBoxLayout(self.pin_chip)
        pin_layout.setContentsMargins(9, 1, 3, 1)
        pin_layout.setSpacing(2)
        self.pin_chip_label = QtWidgets.QLabel("", self.pin_chip)
        self.pin_chip_clear_btn = QtWidgets.QToolButton(self.pin_chip)
        self.pin_chip_clear_btn.setText("×")
        self.pin_chip_clear_btn.setToolTip("Unpin (Esc)")
        self.pin_chip_clear_btn.clicked.connect(self._clear_context_selection)
        pin_layout.addWidget(self.pin_chip_label)
        pin_layout.addWidget(self.pin_chip_clear_btn)
        self.pin_chip.hide()

        self.fit_reference_btn = QtWidgets.QPushButton("Fit reference", card)
        self.fit_reference_btn.setToolTip(
            "Zoom to the reference neighborhood and center the focal gene (F)"
        )
        self.show_all_btn = QtWidgets.QPushButton("Show all", card)
        self.show_all_btn.setToolTip("Zoom out until every row fits in full")
        zoom_row = QtWidgets.QHBoxLayout()
        zoom_row.setContentsMargins(0, 0, 0, 0)
        zoom_row.setSpacing(8)
        for button in (self.fit_reference_btn, self.show_all_btn):
            settings_control(button)
            # Equal widths: the layout splits the row by stretch alone.
            button.setSizePolicy(QtWidgets.QSizePolicy.Policy.Ignored, QtWidgets.QSizePolicy.Policy.Fixed)
            zoom_row.addWidget(button, stretch=1)
        layout.addLayout(zoom_row)

        self.export_btn = QtWidgets.QPushButton("Export figure", card)
        self.export_btn.setObjectName("ContextExportButton")
        self.export_btn.setToolTip("Save the plot as a PNG or PDF figure, or set the figure options")
        self.export_btn.setMenu(self._build_export_menu(self.export_btn))
        layout.addWidget(settings_control(self.export_btn))
        return card

    def _set_molecule_filter(self, molecule: Any) -> None:
        value = str(molecule or "").strip().lower()
        self._molecule_filter = value if value in {"plasmid", "chromosome"} else ""
        select_combo_data(self.molecule_box, self._molecule_filter)
        self._apply_filters()

    def _build_title_block(self) -> QtWidgets.QWidget:
        """Slim title line in the gene card style: eyebrow, name, answer chips."""
        block = QtWidgets.QWidget(self)
        block.setObjectName("ContextTitleBlock")
        # Styled, so its bottom rule sets it off from the Settings column.
        block.setAttribute(QtCore.Qt.WidgetAttribute.WA_StyledBackground, True)
        layout = QtWidgets.QVBoxLayout(block)
        layout.setContentsMargins(16, 11, 16, 9)
        layout.setSpacing(2)
        heading, _sep, subject = self.title_text().partition(":")
        if not subject.strip():
            heading, subject = "Gene neighborhoods", heading
        eyebrow = QtWidgets.QLabel(heading.strip().upper(), block)
        eyebrow.setProperty("role", "annotationEyebrow")
        layout.addWidget(eyebrow)
        self.title_label = QtWidgets.QLabel(subject.strip(), block)
        self.title_label.setProperty("role", "contextTitle")
        self.title_label.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(self.title_label)
        chips = QtWidgets.QWidget(block)
        chip_row = WrapLayout(chips, spacing=6, line_spacing=4)
        chip_row.setContentsMargins(0, 4, 0, 0)
        self._title_chip_layout = chip_row
        self.title_chips: Dict[str, QtWidgets.QLabel] = {}
        self.info = QtWidgets.QLabel("", chips)
        self.info.setProperty("role", "muted")
        layout.addWidget(chips)
        self._render_title_chips()
        return block

    def _render_title_chips(self) -> None:
        layout = self._title_chip_layout
        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget is not None and widget is not self.info:
                widget.deleteLater()
        self.title_chips = {}
        for key, text, tooltip, variant in self._title_chip_specs():
            chip = make_chip(text, tooltip, variant)
            self.title_chips[key] = chip
            layout.addWidget(chip)
        layout.addWidget(self.info)

    def _title_chip_specs(self) -> List[Tuple[str, str, str, str]]:
        """(key, text, tooltip, variant) for each chip under the title.

        In a panel the counts cover the rows in the window's samples; the
        reference stays the comparison's own (every row is scored against it).
        """
        all_rows = self._all_rows
        scoped = self._scope_samples is not None
        rows = [row for row in all_rows if self._row_in_scope(row)] if scoped else all_rows
        narrowed = scoped and len(rows) != len(all_rows)
        specs: List[Tuple[str, str, str, str]] = []
        reference = self._reference_row(all_rows)
        if reference is not None:
            sample = _clean_context_text(reference.get("sample"))
            where = self._row_place_text(reference)
            text = "Reference: " + ", ".join(part for part in (sample, where) if part)
            if self._reference_identity is not None:
                tooltip = (
                    "Chosen with 'Use as reference' (right-click a row). "
                    "Every other row is scored against it."
                )
            elif bool(reference.get("is_source")):
                tooltip = "The neighborhood of the gene you opened. Every other row is scored against it."
            else:
                tooltip = (
                    "The row with the most neighboring genes. Every other row is scored against it; "
                    "right-click a row to use it instead."
                )
            if scoped and not self._row_in_scope(reference):
                tooltip += " Its sample is not in this window's samples, so its row is hidden."
            specs.append(("reference", text, tooltip, "primary"))
        genomes = {_clean_context_text(row.get("sample")) for row in rows} - {""}
        if narrowed:
            every_genome = {_clean_context_text(row.get("sample")) for row in all_rows} - {""}
            specs.append((
                "genomes",
                f"{len(genomes)} of {self._count_text(len(every_genome), 'genome')}",
                "Genomes of this comparison in this window's samples (set in the filter rail)",
                "",
            ))
        else:
            specs.append((
                "genomes",
                self._count_text(len(genomes), "genome"),
                "Genomes with at least one row in this comparison",
                "",
            ))
        kinds = [_row_molecule_kind(row) for row in rows]
        if self._has_focal_gene:
            with_neighbors = {
                _clean_context_text(row.get("sample"))
                for row in rows
                if any(
                    isinstance(gene, dict) and not bool(gene.get("is_center"))
                    for gene in (row.get("genes") or [])
                )
            } - {""}
            specs.append((
                "neighbors",
                f"{len(with_neighbors)} with neighboring genes",
                "Genomes where the focal gene has other genes beside it on the same contig. "
                "In the others it sits on a contig too short to compare; those contigs are "
                "folded into one line at the bottom of the plot.",
                "",
            ))
            parts = [
                f"{count} {label}"
                for label, count in (
                    ("plasmid", kinds.count("plasmid")),
                    ("chromosome", kinds.count("chromosome")),
                    ("mixed", kinds.count("mixed")),
                    ("unknown", kinds.count("unknown")),
                )
                if count
            ]
            if parts:
                # Rows are contigs, and a genome can have several, so this
                # counts contigs to stay apart from the genome chip.
                specs.append((
                    "molecules",
                    f"{self._count_text(len(kinds), 'contig')}: " + ", ".join(parts),
                    "Contigs carrying the focal gene, by molecule type (one row each)",
                    "",
                ))
            flanking = self._flanking_mobile_genes()
            if flanking:
                names = ", ".join(name for name, _count in flanking[:3])
                detail = ", ".join(
                    f"{name} ({self._count_text(count, 'contig')})" for name, count in flanking
                )
                specs.append((
                    "flanked",
                    f"Flanked by: {names}",
                    f"Mobile genes within {MOBILE_FLANK_GENES} genes of the focal gene: {detail}",
                    "mobile",
                ))
        else:
            specs.append((
                "molecules",
                self._count_text(kinds.count("plasmid") or len(rows), "plasmid"),
                "Plasmids in this comparison",
                "",
            ))
        return specs

    def _flanking_mobile_genes(self) -> List[Tuple[str, int]]:
        """Mobile genes near the focal gene, most frequent first, as (name, contigs).

        Only transposases, IS elements, transposons and integrases count,
        found by name or product; other 'mobile' tags such as recombinases
        are left out.  The window is fixed so the chip does not change with
        the Flank switch.
        """
        counts: Dict[str, int] = {}
        for row in self._all_rows:
            genes = row.get("genes")
            center = self._focal_gene_index(row) if isinstance(genes, list) else None
            if center is None or not any(isinstance(g, dict) and g.get("is_center") for g in genes):
                continue
            names: set[str] = set()
            low = max(0, center - MOBILE_FLANK_GENES)
            for gene in genes[low:center + MOBILE_FLANK_GENES + 1]:
                if not isinstance(gene, dict) or gene.get("is_center"):
                    continue
                payload = gene.get("row_data") if isinstance(gene.get("row_data"), dict) else {}
                label = display_label_for_gene(gene)
                product = _clean_context_text(payload.get("Product"))
                if not (_MOBILE_NAME_PATTERN.search(label) or _MOBILE_NAME_PATTERN.search(product)):
                    continue
                if not label or label.lower() == product.lower():
                    label = "transposase" if "transpos" in product.lower() else (label or "mobile gene")
                names.add(label)
            for name in names:
                counts[name] = counts.get(name, 0) + 1
        return sorted(counts.items(), key=lambda item: (-item[1], item[0].lower()))

    def _menu_section_label(self, menu: QtWidgets.QMenu, text: str) -> None:
        label = QtWidgets.QLabel(text, menu)
        label.setProperty("role", "muted")
        label.setContentsMargins(12, 6, 12, 2)
        action = QtWidgets.QWidgetAction(menu)
        action.setDefaultWidget(label)
        action.setEnabled(False)
        menu.addAction(action)

    def _build_export_menu(self, parent: QtWidgets.QWidget) -> QtWidgets.QMenu:
        menu = QtWidgets.QMenu(parent)
        menu.setToolTipsVisible(True)
        menu.addAction("Figure as PNG (2x)", self._export_context_png)
        menu.addAction("Figure as PDF (vector)", self._export_context_pdf)
        menu.addSeparator()
        self._menu_section_label(menu, "Figure options")
        self.publication_style_action = menu.addAction("Publication style")
        self.publication_style_action.setCheckable(True)
        self.publication_style_action.setChecked(False)
        self.publication_style_action.setToolTip(
            "Figures only: white rows, dark gene outlines and quieter ribbons"
        )
        self.publication_style_action.toggled.connect(
            lambda checked: self._set_context_view_option("set_export_publication_style", checked)
        )
        self._context_density_actions: Dict[str, QtGui.QAction] = {}
        self._context_link_actions: Dict[str, QtGui.QAction] = {}
        for attr, method, choices in (
            (
                "_context_density_actions",
                "set_density_mode",
                (("Compact rows", "compact", True), ("Comfortable rows", "comfortable", False)),
            ),
            (
                "_context_link_actions",
                "set_link_mode",
                (
                    ("Links between neighboring rows", "blocks", True),
                    ("All links (dense)", "all", False),
                    ("No links", "off", False),
                ),
            ),
        ):
            menu.addSeparator()
            group = QtGui.QActionGroup(menu)
            group.setExclusive(True)
            for text, value, checked in choices:
                action = menu.addAction(text)
                action.setCheckable(True)
                action.setChecked(checked)
                group.addAction(action)
                getattr(self, attr)[value] = action
                action.triggered.connect(
                    lambda _checked=False, m=method, v=value: self._set_context_view_option(m, v)
                )
        menu.addSeparator()
        self.include_short_contigs_action = menu.addAction("Include short contigs")
        self.include_short_contigs_action.setCheckable(True)
        self.include_short_contigs_action.setChecked(self._export_include_short_contigs)
        self.include_short_contigs_action.setToolTip(
            "Include rows whose contig holds only the focal gene"
        )
        self.include_short_contigs_action.toggled.connect(self._set_export_include_short_contigs)
        return menu

    def _hide_controls_the_view_lacks(self) -> None:
        """The classic renderer has its own coloring and labels: hide switches
        and figure options it would silently ignore."""
        supported = callable(getattr(self.view, "set_color_mode", None))
        for control in (
            self.color_by_group,
            self.labels_group,
            self.publication_style_action,
            self.include_short_contigs_action,
        ):
            control.setVisible(supported)

    def _set_export_include_short_contigs(self, checked: bool) -> None:
        self._export_include_short_contigs = bool(checked)
        self._set_context_view_option("set_export_include_short_contigs", self._export_include_short_contigs)

    def _show_row_menu(self, row: object, global_pos: QtCore.QPoint) -> None:
        """Right-click on a row: re-rank against it, copy or open its samples."""
        if not isinstance(row, dict):
            return
        samples = [
            str(sample).strip()
            for sample in (row.get("collapsed_samples") or [row.get("sample")])
            if str(sample or "").strip()
        ]
        merged = len(samples) > 1
        menu = QtWidgets.QMenu(self)
        menu.setToolTipsVisible(True)
        use = menu.addAction("Use as reference")
        is_reference = self._is_reference_row(row)
        short = bool(row.get("genes")) and self._neighbor_gene_count(row) == 0 and any(
            isinstance(gene, dict) and gene.get("is_center") for gene in row.get("genes") or []
        )
        use.setEnabled(not is_reference and not short)
        if is_reference:
            use.setToolTip("This row is already the reference")
        elif short:
            use.setToolTip("Contig too short to compare: it holds only the focal gene")
        elif row.get("context_partial"):
            use.setToolTip(
                "Score and order every other row against this one. Its contig ends within"
                " the shown flank, so matches would only cover the side it shows"
            )
        else:
            use.setToolTip("Score and order every other row against this one")
        use.triggered.connect(lambda _checked=False, r=dict(row): self._use_row_as_reference(r))
        copy = menu.addAction("Copy accessions" if merged else "Copy accession")
        copy.setEnabled(bool(samples))
        copy.triggered.connect(lambda _checked=False, values=list(samples): self._copy_accessions(values))
        if self._open_samples_callback is not None and samples:
            if self._panel_mode:
                opener = menu.addAction(
                    "Show only these samples in the window" if merged else "Show only this sample in the window"
                )
            else:
                opener = menu.addAction("Open samples in table" if merged else "Open sample in table")
            if self._panel_mode:
                opener.setToolTip(
                    "Narrow this window's samples to " + ("these samples" if merged else "this sample")
                    + " (a chip in the filter rail; every panel follows)"
                )
            else:
                opener.setToolTip("Filter the main table to " + ("these samples" if merged else "this sample"))
            opener.triggered.connect(lambda _checked=False, values=list(samples): self._open_samples(values))
        menu.exec(global_pos)

    def _use_row_as_reference(self, row: Dict[str, Any]) -> None:
        identity = self._context_row_identity(row)
        if identity is None:
            return
        self._reference_identity = identity
        self._render_title_chips()
        self._apply_filters()

    @staticmethod
    def _copy_accessions(samples: Sequence[str]) -> None:
        clipboard = QtWidgets.QApplication.clipboard()
        if clipboard is not None and samples:
            clipboard.setText("\n".join(samples))

    def _open_samples(self, samples: Sequence[str]) -> None:
        callback = self._open_samples_callback
        if callback is not None and samples:
            callback(list(samples))

    def _build_detail_panel(self, parent: QtWidgets.QWidget) -> QtWidgets.QFrame:
        """Pinned gene details as labelled fields, in the gene card's field style."""
        panel = QtWidgets.QFrame(parent)
        panel.setObjectName("ContextInspector")
        panel.setFrameShape(QtWidgets.QFrame.Shape.StyledPanel)
        panel.setMinimumWidth(280)
        panel.setMaximumWidth(420)
        self.synteny_detail_panel = panel
        detail_layout = QtWidgets.QVBoxLayout(panel)
        detail_layout.setContentsMargins(14, 10, 10, 12)
        detail_layout.setSpacing(8)
        detail_header = QtWidgets.QHBoxLayout()
        self.context_inspector_title = QtWidgets.QLabel("Pinned gene", panel)
        self.context_inspector_title.setProperty("role", "sectionTitle")
        detail_close = QtWidgets.QToolButton(panel)
        detail_close.setText("×")
        detail_close.setToolTip("Unpin and close (Esc)")
        detail_close.clicked.connect(self._clear_context_selection)
        detail_header.addWidget(self.context_inspector_title)
        detail_header.addStretch(1)
        detail_header.addWidget(detail_close)
        detail_layout.addLayout(detail_header)

        self.gene_fields = QtWidgets.QWidget(panel)
        fields_layout = QtWidgets.QVBoxLayout(self.gene_fields)
        fields_layout.setContentsMargins(0, 0, 0, 0)
        fields_layout.setSpacing(7)
        self._detail_fields: Dict[str, QtWidgets.QLabel] = {}
        for name in DETAIL_FIELD_NAMES:
            label = QtWidgets.QLabel(name.upper(), self.gene_fields)
            label.setProperty("role", "annotationFieldLabel")
            value = QtWidgets.QLabel("", self.gene_fields)
            value.setProperty("role", "annotationFieldValue")
            value.setWordWrap(True)
            value.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
            box = QtWidgets.QVBoxLayout()
            box.setSpacing(1)
            box.addWidget(label)
            box.addWidget(value)
            fields_layout.addLayout(box)
            self._detail_fields[name] = value
        self.open_gene_card_btn = QtWidgets.QPushButton("Open gene card", self.gene_fields)
        self.open_gene_card_btn.setProperty("role", "primary")
        self.open_gene_card_btn.setToolTip("Open the full gene card (double-click the gene or press Enter)")
        self.open_gene_card_btn.clicked.connect(self._open_pinned_gene_card)
        fields_layout.addSpacing(4)
        fields_layout.addWidget(self.open_gene_card_btn, alignment=QtCore.Qt.AlignmentFlag.AlignLeft)
        detail_layout.addWidget(self.gene_fields)

        # Kept for the classic renderer's homology link details.
        self.synteny_detail_text = QtWidgets.QTextBrowser(panel)
        self.synteny_detail_text.setOpenExternalLinks(False)
        self.synteny_detail_text.hide()
        detail_layout.addWidget(self.synteny_detail_text, stretch=1)
        detail_layout.addStretch(1)
        panel.hide()
        return panel

    def _open_pinned_gene_card(self) -> None:
        opener = getattr(self.view, "open_selected_gene", None)
        if callable(opener):
            opener()

    def _set_context_flank(self, flank: Any) -> None:
        try:
            value = int(flank)
        except (TypeError, ValueError):
            return
        if value <= 0 or value == self._context_flank:
            return
        self._context_flank = value
        self._apply_filters()

    def _set_context_color_mode(self, mode: str) -> None:
        self._context_color_mode = mode
        self._set_context_view_option("set_color_mode", mode)

    def _set_context_label_mode(self, mode: str) -> None:
        self._context_label_mode = mode
        self._set_context_view_option("set_label_mode", mode)

    def _clear_context_selection(self) -> bool:
        """Clear the selected gene and close the details panel.

        Returns True when there was something to clear, so Esc can close the
        window only on a second press.
        """
        cleared = False
        has_selection = getattr(self.view, "has_selection", None)
        if callable(has_selection) and has_selection():
            self.view.clear_selection()
            cleared = True
        if self.synteny_detail_panel.isVisible():
            self.synteny_detail_panel.hide()
            cleared = True
        return cleared

    def keyPressEvent(self, event: QtGui.QKeyEvent) -> None:
        if event.key() == QtCore.Qt.Key.Key_Escape and self._clear_context_selection():
            event.accept()
            return
        super().keyPressEvent(event)

    def _schedule_filter_apply(self) -> None:
        self._filter_timer.start()

    def _set_context_view_option(self, method_name: str, value: str) -> None:
        method = getattr(self.view, method_name, None)
        if callable(method):
            method(value)

    def settings_open(self) -> bool:
        """True while the Settings column is open, False while it is the strip."""
        return not self.filter_rail.isHidden()

    def _toggle_filter_rail(self, visible: bool) -> None:
        self._settings_auto_folded = False
        rail = getattr(self, "filter_rail", None)
        splitter = getattr(self, "main_split", None)
        if not isinstance(rail, QtWidgets.QWidget) or not isinstance(splitter, QtWidgets.QSplitter):
            return
        self._settings_memory.settings_panel_open = bool(visible)
        if bool(visible) != rail.isHidden():
            if not visible or splitter.sizes()[0] > 0:
                return
        # Opening or closing the column only widens or narrows the plot: the
        # zoom and the rows in view stay as they were.
        capture = getattr(self.view, "view_state", None)
        state = capture() if callable(capture) else None
        if visible:
            self.settings_strip.hide()
            rail.show()
            width = max(rail.minimumWidth(), min(420, int(self._filter_rail_last_width)))
            splitter.setSizes([width, max(600, self.width() - width)])
        else:
            sizes = splitter.sizes()
            if sizes and sizes[0] > 0:
                self._filter_rail_last_width = sizes[0]
            rail.hide()
            self.settings_strip.show()
        restore = getattr(self.view, "restore_view_state", None)
        if state is not None and callable(restore):
            # Again once the layout has settled, as the new width lands later.
            # The timer belongs to the view, so closing first cancels it.
            restore(state)
            timer = QtCore.QTimer(self.view)
            timer.setSingleShot(True)
            timer.timeout.connect(lambda: restore(state))
            timer.timeout.connect(timer.deleteLater)
            timer.start(0)

    def apply_settings_default(self, width: int) -> None:
        """Once, before the viewer is first shown at ``width`` pixels.

        No choice made yet this session: a narrow viewer starts with the
        Settings column folded, so the plot keeps its room. Called when the
        viewer is first shown (``showEvent``), so the plot fits its width once.
        """
        if getattr(self, "_settings_default_applied", False):
            return
        self._settings_default_applied = True
        if self._settings_memory.settings_panel_open is None and width < self.SETTINGS_PANEL_AUTO_MIN_WIDTH:
            self.filter_rail.hide()
            self.settings_strip.show()

    def set_settings_open(self, visible: bool) -> None:
        """Open or fold the Settings column without remembering it as the
        session's choice (a panel restoring its saved option)."""
        self._settings_default_applied = True
        memory = self._settings_memory
        remembered = memory.settings_panel_open
        self._toggle_filter_rail(bool(visible))
        memory.settings_panel_open = remembered

    def showEvent(self, event: QtGui.QShowEvent) -> None:
        self.apply_settings_default(self.width())
        super().showEvent(event)
        if self._panel_mode:
            self._settings_fit_timer.start()

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        super().resizeEvent(event)
        if self._panel_mode:
            # Once the layout has settled at the new size.
            self._settings_fit_timer.start()

    def settings_wanted(self) -> bool:
        """The Settings column's chosen state: open, or folded only because
        the panel is small for now (what a panel saves)."""
        return self.settings_open() or self._settings_auto_folded

    def _panel_too_small_for_settings(self) -> bool:
        """A docked panel too small for the column beside a usable plot: the
        cards would clip and the Actions buttons overlap."""
        rail = self.filter_rail
        if self.width() < rail.minimumWidth() + self.PANEL_PLOT_MIN_WIDTH:
            return True
        # Room for the column's header, the pinned Actions card and a card.
        needed = self.settings_footer.sizeHint().height() + 160
        return self.main_split.height() < needed

    def _fit_settings_to_panel(self) -> None:
        if not self._panel_mode or not self.isVisible():
            return
        self._settings_auto_folded = auto_fold(
            self._panel_too_small_for_settings(),
            self.settings_open(),
            self._settings_auto_folded,
            self.set_settings_open,
        )

    def _focus_context_search(self) -> None:
        if not self.settings_open():
            self._toggle_filter_rail(True)
        if self._panel_mode:
            # No search box in a panel: Ctrl+F only opens the settings.
            self.molecule_box.setFocus(QtCore.Qt.FocusReason.ShortcutFocusReason)
            return
        self.sample_query.setFocus(QtCore.Qt.FocusReason.ShortcutFocusReason)
        self.sample_query.selectAll()

    def _show_gene_details(self, payload: Dict[str, Any]) -> None:
        if not payload:
            self.synteny_detail_panel.hide()
            self.pin_chip.hide()
            return
        label = (
            _clean_context_text(payload.get("__context_label__"))
            or _clean_context_text(payload.get("Gene"))
            or "Unnamed gene"
        )
        product = _clean_context_text(payload.get("Product"))
        start = payload.get("__context_start__")
        end = payload.get("__context_end__")
        position = ""
        try:
            low, high = sorted((int(float(start)), int(float(end))))
            position = f"{low:,} to {high:,}"
        except (TypeError, ValueError):
            pass
        place = _clean_context_text(payload.get("__context_place__"))
        sample = _clean_context_text(payload.get("__context_sample__"))
        # The canvas words the position itself, so plasmid offsets never
        # read as contig coordinates.
        where = _clean_context_text(payload.get("__context_position__")) or ", ".join(
            part for part in (place, position) if part
        )
        if sample:
            where = f"{sample}, {where}" if where else sample
        found = _clean_context_text(payload.get("__context_found_in__"))
        if found.lower().startswith("found in "):
            found = found[len("found in "):]
        values = {
            "Gene": label,
            "Product": product or "No product description",
            "Position": where or "Unknown",
            "Strand": _strand_text(payload.get("__context_strand__")),
            "Gene type": _clean_context_text(payload.get("__context_type__"))
            or _clean_context_text(payload.get("__context_category__"))
            or "Other",
            "Found in": found or "Unknown",
        }
        for name, value in values.items():
            self._detail_fields[name].setText(value)
        self.context_inspector_title.setText("Pinned gene")
        self.synteny_detail_text.hide()
        self.gene_fields.show()
        self.synteny_detail_panel.show()
        # A long name (a bare accession, say) is cut short; the tooltip has it.
        self.pin_chip_label.setText(self.pin_chip_label.fontMetrics().elidedText(
            f"Pinned: {label}", QtCore.Qt.TextElideMode.ElideRight, self.PIN_CHIP_MAX_TEXT_WIDTH
        ))
        self.pin_chip.setToolTip(f"{label} and its gene family are highlighted in every row")
        self.pin_chip.show()

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
            row = QtWidgets.QFrame(self.filters_container)
            row.setObjectName("ContextFilterChip")
            row.setSizePolicy(
                QtWidgets.QSizePolicy.Policy.Preferred,
                QtWidgets.QSizePolicy.Policy.Fixed,
            )
            row_layout = QtWidgets.QHBoxLayout(row)
            row_layout.setContentsMargins(0, 0, 0, 0)
            row_layout.setSpacing(6)
            label = QtWidgets.QLabel(column)
            label.setProperty("role", "chipTitle")
            summary = QtWidgets.QLabel(self._format_filter_values(values))
            summary.setProperty("role", "muted")
            summary.setToolTip(self._format_filter_values(values, full=True))
            remove_btn = QtWidgets.QToolButton(row)
            remove_btn.setText("×")
            remove_btn.setProperty("role", "filterChip")
            remove_btn.setToolTip(f"Remove {column} filter")
            remove_btn.clicked.connect(lambda _checked=False, c=column: self._remove_metadata_filter(c))
            row_layout.addWidget(label)
            row_layout.addWidget(summary, stretch=1)
            row_layout.addWidget(remove_btn)
            self.filters_layout.addWidget(row)
        self.filters_layout.addStretch(1)
        self.filters_container.setVisible(bool(self._metadata_filters))
        self._update_filter_badge()

    def _clear_context_filters(self) -> None:
        self._metadata_filters.clear()
        with QtCore.QSignalBlocker(self.sample_query):
            self.sample_query.clear()
        self._molecule_filter = ""
        select_combo_data(self.molecule_box, "")
        self._render_metadata_filters()
        self._apply_filters()

    def _update_filter_badge(self) -> None:
        count = len(self._metadata_filters)
        if self.sample_query.text().strip():
            count += 1
        if self._molecule_filter:
            count += 1
        text = self._count_text(count, "filter") if count else ""
        self.settings_filter_badge.setText(text)
        self.settings_filter_badge.setVisible(bool(count))
        self.settings_strip_badge.setText(str(count) if count else "")
        self.settings_strip_badge.setVisible(bool(count))
        self.settings_strip_badge.setToolTip(f"{text} on" if count else "")
        self.settings_expand_btn.setToolTip(
            f"Show settings (Ctrl+F). {text[:1].upper()}{text[1:]} on" if count else "Show settings (Ctrl+F)"
        )
        self.clear_filters_btn.setEnabled(bool(count))

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

    def _show_synteny_block_details(self, payload: Dict[str, Any]) -> None:
        # Only the classic renderer emits ribbon clicks.
        title = html.escape(str(payload.get("title") or "Synteny block"))
        self.context_inspector_title.setText("Homology link")
        self.gene_fields.hide()
        self.synteny_detail_text.show()
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
        # Name the file after the window title, e.g. the gene name.
        slug = re.sub(r"[^A-Za-z0-9]+", "_", self.title_text()).strip("_").lower()
        return str(base / f"bactolomew_{slug or 'gene_neighborhoods'}.{suffix}")

    def _export_context_png(self) -> None:
        path_text, _selected = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Export gene neighborhoods as PNG",
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
            "Export gene neighborhoods as PDF",
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
        molecule = str(getattr(self, "_molecule_filter", "") or "").strip().lower()
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

    @staticmethod
    def _focal_gene_index(row: Dict[str, Any]) -> Optional[int]:
        genes = row.get("genes")
        if not isinstance(genes, list) or not genes:
            return None
        for idx, gene in enumerate(genes):
            if isinstance(gene, dict) and bool(gene.get("is_center")):
                return idx
        try:
            center_bp = float(row.get("center_bp"))
        except Exception:
            center_bp = None  # type: ignore[assignment]
        if center_bp is not None:
            candidates: List[Tuple[float, int]] = []
            for idx, gene in enumerate(genes):
                if not isinstance(gene, dict):
                    continue
                try:
                    start = float(gene.get("start") or 0)
                    end = float(gene.get("end") or 0)
                except Exception:
                    continue
                low, high = min(start, end), max(start, end)
                if low <= center_bp <= high:
                    return idx
                candidates.append((abs(((low + high) / 2.0) - center_bp), idx))
            if candidates:
                return min(candidates)[1]
        return len(genes) // 2

    def _row_with_selected_flank(self, row: Dict[str, Any]) -> Dict[str, Any]:
        out = dict(row)
        flank = self._context_flank
        genes = row.get("genes")
        if not isinstance(genes, list) or not genes:
            return out
        center_idx = self._focal_gene_index(row)
        if center_idx is None:
            return out
        if flank is None:
            start_idx, end_idx = 0, len(genes)
        else:
            start_idx = max(0, int(center_idx) - int(flank))
            end_idx = min(len(genes), int(center_idx) + int(flank) + 1)
            out["genes"] = [dict(gene) if isinstance(gene, dict) else gene for gene in genes[start_idx:end_idx]]
            out["display_flank"] = int(flank)
        self._mark_contig_ends(
            out,
            shown_before=int(center_idx) - start_idx,
            shown_after=end_idx - 1 - int(center_idx),
            flank=flank,
        )
        return out

    @staticmethod
    def _mark_contig_ends(
        row: Dict[str, Any],
        *,
        shown_before: int,
        shown_after: int,
        flank: Optional[int],
    ) -> None:
        """Mark where the contig really ends, from the gene card's counts.

        'contig_genes_before' and 'contig_genes_after' say how many genes the
        contig holds on each side of the focal gene.  A row whose shown genes
        reach that count stops because the assembly ends there; a contig with
        fewer genes than the flank on a side makes the row 'partial'.  Rows
        without the counts (plasmid context) get neither.
        """
        before = optional_count(row.get("contig_genes_before"))
        after = optional_count(row.get("contig_genes_after"))
        if before is None or after is None or not any(
            isinstance(gene, dict) and gene.get("is_center") for gene in (row.get("genes") or [])
        ):
            return
        row["contig_ends_low"] = shown_before >= before
        row["contig_ends_high"] = shown_after >= after
        has_neighbors = (before + after) > 0
        row["context_partial"] = bool(
            flank is not None and has_neighbors and (before < int(flank) or after < int(flank))
        )

    def _annotate_comparison_metadata(
        self,
        rows: Sequence[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        prepared = [dict(row) for row in rows]
        # Keep similarity scores anchored to the original focal neighborhood.
        # A sample search may hide that row, but it must not silently turn the
        # first remaining result into a new scientific reference.
        reference_source = self._reference_row(self._all_rows)
        reference_row = (
            self._row_with_selected_flank(reference_source)
            if reference_source is not None
            else None
        )
        reference_key = self._context_row_identity(reference_row)
        reference_profiles = (
            self._context_gene_profiles(reference_row) if reference_row is not None else []
        )
        for row in prepared:
            is_reference = bool(
                reference_key is not None
                and self._context_row_identity(row) == reference_key
            )
            score, orientation, match_count = self._context_similarity_details(
                row,
                reference_profiles,
            )
            row["is_comparison_reference"] = bool(is_reference)
            row["context_similarity"] = float(score)
            row["context_orientation"] = "reference" if is_reference else orientation
            row["context_match_count"] = int(match_count)
        return prepared

    @staticmethod
    def _context_row_identity(row: Optional[Dict[str, Any]]) -> Optional[Tuple[str, str, str]]:
        if not isinstance(row, dict):
            return None
        sample = str(row.get("sample") or "").strip()
        contig = str(row.get("contig") or "").strip()
        center = str(row.get("center_bp") or "").strip()
        if not sample and not contig and not center:
            return None
        return sample, contig, center

    def _apply_filters(self) -> None:
        self._update_filter_badge()
        query = self.sample_query.text().strip().lower()
        visible = [
            self._row_with_selected_flank(row)
            for row in self._all_rows
            if self._row_in_scope(row) and self._row_matches_filters(row, query)
        ]
        visible = self._annotate_comparison_metadata(visible)
        visible = self._sort_rows(visible)
        self._visible_rows = list(visible)
        display_rows = self._collapse_rows(visible) if self._collapsed else visible
        if self._collapsed and self._group_by_molecule_enabled():
            # Merging can turn a group 'mixed'; keep each molecule in one block.
            display_rows = self._split_by_molecule(display_rows)
        # Presence counts in tooltips cover the rows shown; say so when row
        # filters hide some, so they never read as the title chip's totals.
        self._set_context_view_option("set_rows_filtered", len(visible) != len(self._all_rows))
        self._set_context_view_option("set_molecule_split", self._group_by_molecule_enabled())
        if self._panel_mode:
            # The window's samples leave no row: say where to change them.
            out_of_scope = bool(self._all_rows) and self._scope_samples is not None and not self.scope_row_count()
            self._set_context_view_option("set_empty_text", self.OUT_OF_SCOPE_TEXT if out_of_scope else "")
        self.view.set_rows(display_rows)
        self._set_context_view_option("set_sort_caption", self._sort_caption())
        status = self._info_text(visible, display_rows)
        self.info.setText(status)
        self.info.setVisible(bool(status))

    # The window's sample set (panel mode, D4)

    def _row_in_scope(self, row: Dict[str, Any]) -> bool:
        scope = self._scope_samples
        if scope is None:
            return True
        samples = row.get("collapsed_samples") or [row.get("sample")]
        return any(str(sample or "").strip() in scope for sample in samples)

    def set_sample_scope(self, samples: Optional[Sequence[str]], notes: Sequence[str] = ()) -> None:
        """Show only the rows of ``samples`` (the window's one sample set;
        None shows every row). ``notes`` name the window's active filters,
        shown read only in the Show card so a short list explains itself."""
        scope = frozenset(str(s).strip() for s in samples or () if str(s).strip()) if samples is not None else None
        notes = [str(note) for note in notes if str(note).strip()]
        if scope == self._scope_samples and notes == self._scope_notes:
            return
        self._scope_samples = scope
        self._scope_notes = notes
        self._render_scope()
        self._render_title_chips()
        self._apply_filters()

    def scope_row_count(self) -> int:
        return sum(1 for row in self._all_rows if self._row_in_scope(row))

    def _render_scope(self) -> None:
        label = self.scope_label
        layout = self._scope_chip_layout
        if label is None or layout is None:
            return
        total = len(self._all_rows)
        shown = self.scope_row_count()
        genomes = {_clean_context_text(row.get("sample")) for row in self._all_rows} - {""}
        in_scope = {
            _clean_context_text(row.get("sample")) for row in self._all_rows if self._row_in_scope(row)
        } - {""}
        if self._scope_samples is None or shown == total:
            text = "Every row is in this window's samples."
        else:
            text = (
                f"{self._count_text(shown, 'row')} of {total} ({len(in_scope)} of "
                f"{self._count_text(len(genomes), 'genome')}) are in this window's samples."
            )
        if self._scope_notes:
            text += " Filters in the rail:"
        label.setText(text)
        label.setToolTip(
            "The rows follow the samples of this window, set in the filter rail. "
            "Change or clear the filters there to show more rows."
        )
        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
        for note in self._scope_notes:
            chip = make_chip(note, "A filter of this window (change it in the filter rail)")
            chip.setParent(self.scope_chips)
            layout.addWidget(chip)
        self.scope_chips.setVisible(bool(self._scope_notes))

    def set_linked_samples(self, samples: Optional[Sequence[str]]) -> None:
        """Samples hovered elsewhere in the window: their rows are tinted."""
        key = tuple(self._normalized_samples(samples))
        if key == self._linked_samples:
            return
        self._linked_samples = key
        self._set_context_view_option("set_linked_samples", list(key))

    @staticmethod
    def _count_text(count: int, singular: str, plural: Optional[str] = None) -> str:
        return f"{count} {singular if count == 1 else (plural or singular + 's')}"

    def _info_text(
        self,
        visible: Sequence[Dict[str, Any]],
        display_rows: Sequence[Dict[str, Any]],
    ) -> str:
        """What the plot shows, when filters or merging change it; else ''.

        The title chips already describe the whole comparison.
        """
        total = len(self._all_rows)
        text = ""
        if len(visible) != total:
            text = f"Showing {len(visible)} of {self._count_text(total, 'neighborhood')}"
        if self._collapsed:
            merged = f"merged into {self._count_text(len(display_rows), 'row')}"
            text = f"{text}, {merged}" if text else f"{self._count_text(total, 'neighborhood')} {merged}"
        return text

    @staticmethod
    def _row_place_text(row: Dict[str, Any]) -> str:
        """Where a row sits: 'contig 48', or the plasmid name for plasmid rows.

        Plasmid-context rows (one per plasmid, with a 'cluster_key') carry an
        instance key in 'contig', which is not a contig name, so their display
        title is used instead.  A gene row on a plasmid contig keeps its contig.
        """
        if bool(row.get("is_plasmid")) and _clean_context_text(row.get("cluster_key")):
            title = _clean_context_text(row.get("title"))
            sample = _clean_context_text(row.get("sample"))
            # 'sample, plasmid pX1' (or the older 'sample | pX1').
            for separator in ("|", ","):
                head, found, tail = title.partition(separator)
                if found and (not sample or head.strip() == sample):
                    title = tail.strip()
                    break
            return f"plasmid {title}" if title and not title.lower().startswith("plasmid") else title
        return _contig_text(_clean_context_text(row.get("contig")))

    def _sort_caption(self) -> str:
        data = self.sort_column_box.currentData() if hasattr(self, "sort_column_box") else None
        order = self.sort_order_box.currentData() if hasattr(self, "sort_order_box") else None
        ascending = str(order or "asc").strip().lower().startswith("asc")
        sort_col = str(data or "__sample__")
        if sort_col == "__context_similarity__":
            caption = "Reference first, then by match to the reference, " + (
                "lowest first" if ascending else "highest first"
            )
        else:
            name = {"__sample__": "sample", "__contig__": "contig"}.get(sort_col, sort_col)
            caption = f"Sorted by {name}, " + ("A to Z" if ascending else "Z to A")
        if self._group_by_molecule_enabled():
            caption += "; plasmid contigs before chromosome contigs"
        if self._collapsed:
            caption += ". Identical neighborhoods merged"
        return caption

    def _toggle_collapsed(self, checked: bool) -> None:
        # The checkable button shows the state; its name stays the same.
        self._collapsed = bool(checked)
        self._apply_filters()

    def _reference_row(
        self,
        rows: Optional[Sequence[Dict[str, Any]]] = None,
    ) -> Optional[Dict[str, Any]]:
        """The row every other row is scored against.

        A row picked with 'Use as reference' wins, then the gene the viewer
        opened ('is_source').  Without either (plasmid context) the row with
        the most neighboring genes is the fairest yardstick.
        """
        candidates = list(rows) if rows is not None else self._all_rows
        if self._reference_identity is not None:
            for row in candidates:
                if self._context_row_identity(row) == self._reference_identity:
                    return row
        for row in candidates:
            if bool(row.get("is_source")):
                return row
        if not candidates:
            return None
        return max(
            enumerate(candidates),
            key=lambda item: (self._neighbor_gene_count(item[1]), -item[0]),
        )[1]

    @staticmethod
    def _neighbor_gene_count(row: Dict[str, Any]) -> int:
        genes = row.get("genes")
        if not isinstance(genes, list):
            return 0
        return sum(1 for gene in genes if isinstance(gene, dict) and not gene.get("is_center"))

    @staticmethod
    def _is_reference_row(row: Dict[str, Any]) -> bool:
        """Scored rows carry 'is_comparison_reference', which wins over 'is_source'."""
        if "is_comparison_reference" in row:
            return bool(row.get("is_comparison_reference"))
        return bool(row.get("is_source"))

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
                    "strand": NeighbourhoodView._normalized_context_strand(gene.get("strand")),
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

    def _context_similarity_details(
        self,
        row: Dict[str, Any],
        reference_profiles: Optional[Sequence[Dict[str, Any]]] = None,
    ) -> Tuple[float, str, int]:
        ref_profiles = list(reference_profiles or [])
        if not ref_profiles:
            ref_row = self._reference_row()
            ref_profiles = self._context_gene_profiles(ref_row) if ref_row is not None else []
        cand_profiles = self._context_gene_profiles(row)
        if not ref_profiles and not cand_profiles:
            return 1.0, "forward", 0
        if not ref_profiles or not cand_profiles:
            return 0.0, "forward", 0

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
        if reverse_score > forward_score:
            raw_score = reverse_score
            match_count = reverse_count
            orientation = "reverse"
        else:
            raw_score = forward_score
            match_count = forward_count
            orientation = "forward"
        if match_count >= 3:
            raw_score *= 1.0 + min(0.28, 0.035 * float(match_count - 2))
        denom = max(
            sum(float(profile.get("weight") or 0.0) for profile in ref_profiles),
            sum(float(profile.get("weight") or 0.0) for profile in cand_profiles),
            1.0,
        )
        return max(0.0, min(1.0, raw_score / denom)), orientation, int(match_count)

    def _context_similarity_score(
        self,
        row: Dict[str, Any],
        reference_profiles: Optional[Sequence[Dict[str, Any]]] = None,
    ) -> float:
        score, _orientation, _match_count = self._context_similarity_details(
            row,
            reference_profiles,
        )
        return score

    @staticmethod
    def _sample_label(sample_count: int) -> str:
        return f"{sample_count} genome" if sample_count == 1 else f"{sample_count} genomes"

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
            f"{self._sample_label(sample_count)} with identical neighborhoods",
            f"{track_count} contig{'s' if track_count != 1 else ''} merged into this row",
        ]
        if sample_names:
            lines.extend(["", "Genomes:"])
            lines.extend(sample_names[:20])
            if len(sample_names) > 20:
                lines.append(f"... (+{len(sample_names) - 20} more)")
        if track_count > sample_count and track_titles:
            lines.extend(["", "Contigs:"])
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
                    0 if self._is_reference_row(row) else 1,
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
                    (
                        f"{member.get('sample') or ''}, {self._row_place_text(member)}"
                        if self._row_place_text(member)
                        else str(member.get("sample") or "")
                    )
                    for member in members
                ]
            )
            sample_count = len(sample_names)
            track_count = len(members)
            molecule_counts: Dict[str, int] = {}
            for member in members:
                kind = _row_molecule_kind(member)
                molecule_counts[kind] = molecule_counts.get(kind, 0) + 1
            representative["display_label"] = self._sample_label(sample_count)
            # 'molecule' wins over 'molecule_kind', so both are set.  A row
            # merging several types is 'mixed'; the per-type counts let the
            # fold line, dividers and chip tooltip say what it holds.
            if len(molecule_counts) == 1:
                representative["molecule"] = representative["molecule_kind"] = next(iter(molecule_counts))
            elif molecule_counts:
                representative["molecule"] = representative["molecule_kind"] = "mixed"
            representative["collapsed_molecule_counts"] = dict(molecule_counts)
            # Members share their shown genes, so they are all partial or
            # none is; a contig-end cap holds only when every member ends there.
            for end in ("contig_ends_low", "contig_ends_high"):
                representative[end] = all(bool(member.get(end)) for member in members)
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

    @staticmethod
    def _molecule_block_rank(row: Dict[str, Any]) -> int:
        return {"plasmid": 0, "mixed": 1, "chromosome": 2}.get(_row_molecule_kind(row), 3)

    def _split_by_molecule(self, rows: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Stable reorder into molecule blocks, the reference row first."""
        return sorted(
            rows,
            key=lambda row: (0 if self._is_reference_row(row) else 1, self._molecule_block_rank(row)),
        )

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
            ref_row = self._reference_row(rows)
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

            scored_rows: List[Tuple[float, Dict[str, Any]]] = []
            for row in rows:
                try:
                    score = float(row.get("context_similarity"))
                except Exception:
                    score = self._context_similarity_score(row, ref_profiles)
                scored_rows.append((score, row))
            if ascending:
                return [
                    row
                    for _score, row in sorted(
                        scored_rows,
                        key=lambda scored: (
                            0 if self._is_reference_row(scored[1]) else 1,
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
                        0 if self._is_reference_row(scored[1]) else 1,
                        _molecule_rank(scored[1]),
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
        ordered = [row for _value, row in present_sorted] + missing_sorted
        if self._group_by_molecule_enabled():
            # 'Split plasmid and chromosome' splits any order into blocks.
            ordered.sort(key=self._molecule_block_rank)
        return ordered
