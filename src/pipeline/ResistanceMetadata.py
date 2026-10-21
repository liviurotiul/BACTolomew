from __future__ import annotations

import json
import logging
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd


SUBJECT_METADATA_COLUMNS = [
    "SubjectIdentifier",
    "SubjectProteinAccession",
    "SubjectDNAAccession",
    "SubjectAROAccession",
    "SubjectModelName",
]

RESISTANCE_MARKER_COLUMNS = [
    "ResistanceMetadataSource",
    "ResistanceAlleles",
    "ResistanceGeneFamilies",
    "ResistanceDrugClasses",
    "ResistanceAntibiotics",
    "ResistanceMechanisms",
]

ACCESSION_RE = re.compile(r"^[A-Za-z]{1,6}_?\d+(?:\.\d+)?$")
ARO_RE = re.compile(r"\bARO:\d+\b")

# CARD protein models, as the header prefix of CARD_AMR.fasta records. Only protein homolog
# models (HM) mean "this acquired gene is present"; variant (VM), knockout (KM) and
# overexpression (OM) models describe a resistance-relevant target whose wild type is
# susceptible, so a homology hit to them is not a resistance call (RGI needs a mutation check).
CARD_HOMOLOG_MODEL_PREFIX = "HM"
CARD_TARGET_MODEL_PREFIXES = frozenset({"VM", "KM", "OM"})
CARD_HOMOLOG_MODEL_TYPE = "protein homolog model"

# AMRFinderPlus catalog subtypes that describe point mutations or susceptible alleles of a
# target gene, not an acquired gene. Hits to these references are never reported as resistance.
AMRF_TARGET_SUBTYPES = frozenset({"POINT", "POINT_DISRUPT", "AMR-SUSCEPTIBLE"})

AMRF_ACCESSION_COLUMNS = (
    "refseq_protein_accession",
    "genbank_protein_accession",
    "refseq_nucleotide_accession",
    "genbank_nucleotide_accession",
)


def _safe_text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    return text


def _normalize_resistance_value(value: Any) -> Any:
    text = _safe_text(value)
    return text or pd.NA


def _join_unique(values: list[Any]) -> str:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        text = _safe_text(value)
        if not text or text in seen:
            continue
        seen.add(text)
        out.append(text)
    return ",".join(out)


def _parse_subject_metadata(product: Any) -> dict[str, Any]:
    header = _safe_text(product)
    subject_id = header.split()[0] if header else ""
    aro = ""
    protein = ""
    dna = ""
    model = ""

    if header:
        m = ARO_RE.search(header)
        if m:
            aro = m.group(0)

    if subject_id:
        parts = [p.strip() for p in subject_id.split("|") if p.strip()]
        if parts:
            candidates = [p for p in parts if ACCESSION_RE.match(p)]
            if candidates:
                protein = candidates[0]
                if len(candidates) > 1:
                    dna = candidates[1]
            if len(parts) >= 5:
                model = parts[4]
        if not protein and ACCESSION_RE.match(subject_id):
            protein = subject_id

    return {
        "SubjectIdentifier": subject_id or pd.NA,
        "SubjectProteinAccession": protein or pd.NA,
        "SubjectDNAAccession": dna or pd.NA,
        "SubjectAROAccession": aro or pd.NA,
        "SubjectModelName": model or pd.NA,
    }


def card_model_prefix(product: Any) -> str:
    """Return the CARD model prefix (HM, VM, KM, OM) of a CARD_AMR header, or ''."""
    subject_id = _safe_text(product).split()
    if not subject_id:
        return ""
    head = subject_id[0].split("|", 1)
    if len(head) < 2:
        return ""
    prefix = head[0].strip().upper()
    return prefix if len(prefix) == 2 and prefix.endswith("M") else ""


def parse_marker_subject_metadata(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in SUBJECT_METADATA_COLUMNS:
        if col not in out.columns:
            out[col] = pd.NA
    if out.empty:
        return out

    parsed = out.get("Product", pd.Series(dtype="object")).apply(_parse_subject_metadata)
    parsed_df = pd.DataFrame(parsed.tolist(), index=out.index)
    for col in SUBJECT_METADATA_COLUMNS:
        if col in parsed_df.columns:
            out[col] = parsed_df[col]
    return out


def _payload_key(payload: dict[str, Any]) -> tuple[str, ...]:
    return tuple(_safe_text(payload.get(col)) for col in RESISTANCE_MARKER_COLUMNS)


def _resolve_payloads(payloads: list[dict[str, Any]]) -> tuple[dict[str, Any], bool]:
    """Collapse the catalog rows of one key into one payload without guessing.

    One distinct row: that row. Several distinct rows: keep only the fields every row agrees
    on (for example the gene family) and leave the allele and every disputed field empty.
    Returns (payload, ambiguous).
    """
    distinct: dict[tuple[str, ...], dict[str, Any]] = {}
    for payload in payloads:
        distinct.setdefault(_payload_key(payload), payload)
    if len(distinct) == 1:
        return next(iter(distinct.values())), False
    merged: dict[str, Any] = {}
    for col in RESISTANCE_MARKER_COLUMNS:
        values = {_safe_text(p.get(col)) for p in distinct.values()}
        merged[col] = _normalize_resistance_value(values.pop()) if len(values) == 1 else pd.NA
    merged["ResistanceAlleles"] = pd.NA
    return merged, True


def _resolve_index(rows_by_key: dict[str, list[dict[str, Any]]]) -> tuple[dict[str, dict[str, Any]], frozenset[str]]:
    index: dict[str, dict[str, Any]] = {}
    ambiguous: set[str] = set()
    for key, payloads in rows_by_key.items():
        payload, is_ambiguous = _resolve_payloads(payloads)
        index[key] = payload
        if is_ambiguous:
            ambiguous.add(key)
    return index, frozenset(ambiguous)


@dataclass(frozen=True)
class AmrfCatalogIndex:
    by_accession: dict[str, dict[str, Any]]
    # Accessions whose catalog rows are all point-mutation / susceptible-allele references.
    target_accessions: frozenset[str]
    # Accessions with several distinct reportable rows (allele left empty).
    ambiguous_accessions: frozenset[str]


@lru_cache(maxsize=8)
def _load_amrf_catalog(path_str: str) -> AmrfCatalogIndex:
    path = Path(path_str)
    if not path.exists() or not path.is_file():
        return AmrfCatalogIndex({}, frozenset(), frozenset())

    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    reportable: dict[str, list[dict[str, Any]]] = {}
    target_only: set[str] = set()
    for row in df.to_dict(orient="records"):
        subtype = _safe_text(row.get("subtype")).upper()
        accessions = {_safe_text(row.get(col)) for col in AMRF_ACCESSION_COLUMNS} - {""}
        if subtype in AMRF_TARGET_SUBTYPES:
            target_only.update(accessions)
            continue
        payload = {
            "ResistanceMetadataSource": "AMRF",
            "ResistanceAlleles": _normalize_resistance_value(row.get("allele")),
            "ResistanceGeneFamilies": _normalize_resistance_value(row.get("gene_family")),
            "ResistanceDrugClasses": _normalize_resistance_value(row.get("class")),
            "ResistanceAntibiotics": _normalize_resistance_value(row.get("subclass")),
            "ResistanceMechanisms": pd.NA,
        }
        for accession in accessions:
            reportable.setdefault(accession, []).append(payload)
    by_accession, ambiguous = _resolve_index(reportable)
    return AmrfCatalogIndex(
        by_accession=by_accession,
        target_accessions=frozenset(target_only - set(by_accession)),
        ambiguous_accessions=ambiguous,
    )


def load_card_target_aros(card_json: Path) -> frozenset[str]:
    """AROs whose CARD models are all non-homolog (variant, knockout, overexpression, rRNA, meta)."""
    if not card_json.is_file():
        return frozenset()
    try:
        data = json.loads(card_json.read_text(encoding="utf-8"))
    except Exception:
        return frozenset()
    types_by_aro: dict[str, set[str]] = {}
    for model in data.values() if isinstance(data, dict) else []:
        if not isinstance(model, dict) or "model_type" not in model:
            continue
        aro = _safe_text(model.get("ARO_accession"))
        if aro:
            types_by_aro.setdefault(f"ARO:{aro}", set()).add(_safe_text(model.get("model_type")).lower())
    return frozenset(aro for aro, types in types_by_aro.items() if CARD_HOMOLOG_MODEL_TYPE not in types)


@dataclass(frozen=True)
class CardIndex:
    by_aro: dict[str, dict[str, Any]]
    by_protein: dict[str, dict[str, Any]]
    # AROs that only have variant/knockout/overexpression models (from card.json next to the index).
    target_aros: frozenset[str]
    ambiguous_keys: frozenset[str]


@lru_cache(maxsize=8)
def _load_card_index(path_str: str) -> CardIndex:
    path = Path(path_str)
    if not path.exists() or not path.is_file():
        return CardIndex({}, {}, frozenset(), frozenset())

    target_aros = load_card_target_aros(path.parent / "card.json")
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    rows_by_aro: dict[str, list[dict[str, Any]]] = {}
    rows_by_protein: dict[str, list[dict[str, Any]]] = {}
    for row in df.to_dict(orient="records"):
        aro = _safe_text(row.get("ARO Accession"))
        if aro in target_aros:
            continue
        payload = {
            "ResistanceMetadataSource": "CARD",
            "ResistanceAlleles": _normalize_resistance_value(row.get("ARO Name") or row.get("CARD Short Name")),
            "ResistanceGeneFamilies": _normalize_resistance_value(row.get("AMR Gene Family")),
            "ResistanceDrugClasses": _normalize_resistance_value(row.get("Drug Class")),
            "ResistanceAntibiotics": pd.NA,
            "ResistanceMechanisms": _normalize_resistance_value(row.get("Resistance Mechanism")),
        }
        protein = _safe_text(row.get("Protein Accession"))
        if aro:
            rows_by_aro.setdefault(aro, []).append(payload)
        if protein:
            rows_by_protein.setdefault(protein, []).append(payload)
    by_aro, ambiguous_aro = _resolve_index(rows_by_aro)
    by_protein, ambiguous_protein = _resolve_index(rows_by_protein)
    return CardIndex(
        by_aro=by_aro,
        by_protein=by_protein,
        target_aros=target_aros,
        ambiguous_keys=ambiguous_aro | ambiguous_protein,
    )


def _normalize_metadata_by_db(metadata_by_db: dict[str, dict[str, Any]] | None) -> dict[str, dict[str, Any]]:
    normalized: dict[str, dict[str, Any]] = {}
    for db_name, raw_meta in (metadata_by_db or {}).items():
        if not isinstance(raw_meta, dict):
            continue
        source = _safe_text(raw_meta.get("source")).upper()
        files_raw = raw_meta.get("files") or []
        if isinstance(files_raw, (str, Path)):
            files = [str(Path(files_raw).expanduser())]
        elif isinstance(files_raw, (list, tuple, set)):
            files = [str(Path(v).expanduser()) for v in files_raw if _safe_text(v)]
        else:
            files = []
        if source and files:
            normalized[str(db_name)] = {"source": source, "files": files}
    return normalized


def _load_source_index(meta: dict[str, Any], logger: logging.Logger | None) -> Any:
    path_str = meta["files"][0]
    if not Path(path_str).exists() and logger is not None:
        logger.warning("%s resistance metadata file not found: %s", meta["source"], path_str)
    if meta["source"] == "AMRF":
        return _load_amrf_catalog(path_str)
    if meta["source"] == "CARD":
        return _load_card_index(path_str)
    return None


def _is_target_reference(source: str, index: Any, product: Any) -> bool:
    if source == "CARD":
        if card_model_prefix(product) in CARD_TARGET_MODEL_PREFIXES:
            return True
        aro = _safe_text(_parse_subject_metadata(product)["SubjectAROAccession"])
        return bool(index is not None and aro and aro in index.target_aros)
    if source == "AMRF" and index is not None:
        accession = _safe_text(_parse_subject_metadata(product)["SubjectProteinAccession"])
        return bool(accession and accession in index.target_accessions)
    return False


def drop_resistance_target_hits(
    df: pd.DataFrame,
    metadata_by_db: dict[str, dict[str, Any]] | None,
    *,
    logger: logging.Logger | None = None,
) -> pd.DataFrame:
    """Remove marker hits to resistance *target* references.

    CARD variant/knockout/overexpression models and AMRFinderPlus POINT, POINT_DISRUPT and
    AMR-SUSCEPTIBLE references match the wild-type housekeeping gene (gyrA, parC, oprD, ...),
    so a homology hit says nothing about resistance. The marker databases no longer contain
    them; this guard keeps hits from an older database or older DIAMOND TSVs out as well.
    Run it before collapsing hits per (AnnotationID, DB) so a valid lower-ranked hit survives.
    """
    if df.empty or "Product" not in df.columns or "DB" not in df.columns:
        return df
    normalized = _normalize_metadata_by_db(metadata_by_db)
    if not normalized:
        return df

    drop = pd.Series(False, index=df.index)
    for db_name, meta in normalized.items():
        if meta["source"] not in ("AMRF", "CARD"):
            continue
        in_db = df["DB"].astype("string") == db_name
        if not in_db.any():
            continue
        index = _load_source_index(meta, logger)
        products = df.loc[in_db.fillna(False), "Product"]
        verdict = {p: _is_target_reference(meta["source"], index, p) for p in products.dropna().unique()}
        drop.loc[products.index] = products.map(verdict).fillna(False).astype(bool)
    n_drop = int(drop.sum())
    if not n_drop:
        return df
    if logger is not None:
        per_db = df.loc[drop, "DB"].value_counts().to_dict()
        logger.warning(
            "Dropped %d marker hit(s) to resistance target models (CARD variant/knockout/"
            "overexpression, AMRFinderPlus point-mutation/susceptible references), which are not "
            "resistance calls: %s. Rebuild the marker databases (python -m pipeline.MarkerDbBuild) "
            "and rerun the DIAMOND marker step if this run used an older database.",
            n_drop,
            per_db,
        )
    return df.loc[~drop].copy()


def enrich_marker_hits_with_resistance(
    df: pd.DataFrame,
    metadata_by_db: dict[str, dict[str, Any]] | None,
    *,
    logger: logging.Logger | None = None,
) -> pd.DataFrame:
    if metadata_by_db:
        df = drop_resistance_target_hits(df, metadata_by_db, logger=logger)
    out = parse_marker_subject_metadata(df)
    for col in RESISTANCE_MARKER_COLUMNS:
        if col not in out.columns:
            out[col] = pd.NA
    if out.empty or not metadata_by_db:
        return out

    normalized_meta = _normalize_metadata_by_db(metadata_by_db)
    index_cache: dict[str, Any] = {}

    payload_by_index: dict[Any, dict[str, Any]] = {}
    for idx, row in out.iterrows():
        db_name = _safe_text(row.get("DB"))
        meta = normalized_meta.get(db_name)
        if not meta:
            continue
        if db_name not in index_cache:
            index_cache[db_name] = _load_source_index(meta, logger)
        index = index_cache[db_name]
        if index is None:
            continue

        payload: dict[str, Any] | None = None
        if meta["source"] == "AMRF":
            protein = _safe_text(row.get("SubjectProteinAccession"))
            dna = _safe_text(row.get("SubjectDNAAccession"))
            payload = index.by_accession.get(protein) or index.by_accession.get(dna)
        elif meta["source"] == "CARD":
            # The ARO in the header names the hit's own model, so it decides first; the protein
            # accession is a fallback whose ambiguous entries carry no allele.
            aro = _safe_text(row.get("SubjectAROAccession"))
            protein = _safe_text(row.get("SubjectProteinAccession"))
            payload = index.by_aro.get(aro) or index.by_protein.get(protein)

        if payload:
            payload_by_index[idx] = payload

    if not payload_by_index:
        return out

    payload_df = pd.DataFrame.from_dict(payload_by_index, orient="index")
    for col in RESISTANCE_MARKER_COLUMNS:
        if col in payload_df.columns:
            out.loc[payload_df.index, col] = payload_df[col]
    return out


def build_annotation_resistance_maps(marker_hits_df: pd.DataFrame) -> dict[str, dict[str, str]]:
    maps = {col: {} for col in RESISTANCE_MARKER_COLUMNS}
    if marker_hits_df.empty or "AnnotationID" not in marker_hits_df.columns:
        return maps

    grouped = marker_hits_df.groupby("AnnotationID", sort=False, dropna=False)
    for col in RESISTANCE_MARKER_COLUMNS:
        if col not in marker_hits_df.columns:
            continue
        maps[col] = grouped[col].agg(lambda s: _join_unique(s.tolist())).to_dict()
    return maps


__all__ = [
    "SUBJECT_METADATA_COLUMNS",
    "RESISTANCE_MARKER_COLUMNS",
    "CARD_HOMOLOG_MODEL_PREFIX",
    "CARD_TARGET_MODEL_PREFIXES",
    "AMRF_TARGET_SUBTYPES",
    "card_model_prefix",
    "parse_marker_subject_metadata",
    "drop_resistance_target_hits",
    "enrich_marker_hits_with_resistance",
    "build_annotation_resistance_maps",
]
