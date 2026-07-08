from __future__ import annotations

import logging
import re
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


@lru_cache(maxsize=8)
def _load_amrf_accession_index(path_str: str) -> dict[str, dict[str, Any]]:
    path = Path(path_str)
    if not path.exists() or not path.is_file():
        return {}

    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    index: dict[str, dict[str, Any]] = {}
    accession_cols = [
        "refseq_protein_accession",
        "genbank_protein_accession",
        "refseq_nucleotide_accession",
        "genbank_nucleotide_accession",
    ]
    for _, row in df.iterrows():
        payload = {
            "ResistanceMetadataSource": "AMRF",
            "ResistanceAlleles": _normalize_resistance_value(row.get("allele")),
            "ResistanceGeneFamilies": _normalize_resistance_value(row.get("gene_family")),
            "ResistanceDrugClasses": _normalize_resistance_value(row.get("class")),
            "ResistanceAntibiotics": _normalize_resistance_value(row.get("subclass")),
            "ResistanceMechanisms": pd.NA,
        }
        for col in accession_cols:
            accession = _safe_text(row.get(col))
            if accession and accession not in index:
                index[accession] = payload
    return index


@lru_cache(maxsize=8)
def _load_card_index(path_str: str) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    path = Path(path_str)
    if not path.exists() or not path.is_file():
        return {}, {}

    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    by_aro: dict[str, dict[str, Any]] = {}
    by_protein: dict[str, dict[str, Any]] = {}
    for _, row in df.iterrows():
        payload = {
            "ResistanceMetadataSource": "CARD",
            "ResistanceAlleles": _normalize_resistance_value(row.get("ARO Name") or row.get("CARD Short Name")),
            "ResistanceGeneFamilies": _normalize_resistance_value(row.get("AMR Gene Family")),
            "ResistanceDrugClasses": _normalize_resistance_value(row.get("Drug Class")),
            "ResistanceAntibiotics": pd.NA,
            "ResistanceMechanisms": _normalize_resistance_value(row.get("Resistance Mechanism")),
        }
        aro = _safe_text(row.get("ARO Accession"))
        protein = _safe_text(row.get("Protein Accession"))
        if aro and aro not in by_aro:
            by_aro[aro] = payload
        if protein and protein not in by_protein:
            by_protein[protein] = payload
    return by_aro, by_protein


def enrich_marker_hits_with_resistance(
    df: pd.DataFrame,
    metadata_by_db: dict[str, dict[str, Any]] | None,
    *,
    logger: logging.Logger | None = None,
) -> pd.DataFrame:
    out = parse_marker_subject_metadata(df)
    for col in RESISTANCE_MARKER_COLUMNS:
        if col not in out.columns:
            out[col] = pd.NA
    if out.empty or not metadata_by_db:
        return out

    normalized_meta: dict[str, dict[str, Any]] = {}
    for db_name, raw_meta in (metadata_by_db or {}).items():
        if not isinstance(raw_meta, dict):
            continue
        source = _safe_text(raw_meta.get("source"))
        files_raw = raw_meta.get("files") or []
        if isinstance(files_raw, (str, Path)):
            files = [str(Path(files_raw).expanduser())]
        elif isinstance(files_raw, (list, tuple, set)):
            files = [str(Path(v).expanduser()) for v in files_raw if _safe_text(v)]
        else:
            files = []
        if source and files:
            normalized_meta[str(db_name)] = {
                "source": source,
                "files": files,
            }

    amrf_cache: dict[str, dict[str, dict[str, Any]]] = {}
    card_cache: dict[str, tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]] = {}

    payload_by_index: dict[int, dict[str, Any]] = {}
    for idx, row in out.iterrows():
        db_name = _safe_text(row.get("DB"))
        meta = normalized_meta.get(db_name)
        if not meta:
            continue

        source = _safe_text(meta.get("source")).upper()
        files = [str(Path(p).expanduser()) for p in meta.get("files", []) if _safe_text(p)]
        if not files:
            continue

        payload: dict[str, Any] | None = None
        if source == "AMRF":
            path_str = files[0]
            if path_str not in amrf_cache:
                if not Path(path_str).exists() and logger is not None:
                    logger.warning("AMRF resistance metadata file not found: %s", path_str)
                amrf_cache[path_str] = _load_amrf_accession_index(path_str)
            index_map = amrf_cache[path_str]
            protein = _safe_text(row.get("SubjectProteinAccession"))
            dna = _safe_text(row.get("SubjectDNAAccession"))
            payload = index_map.get(protein) or index_map.get(dna)
        elif source == "CARD":
            path_str = files[0]
            if path_str not in card_cache:
                if not Path(path_str).exists() and logger is not None:
                    logger.warning("CARD resistance metadata file not found: %s", path_str)
                card_cache[path_str] = _load_card_index(path_str)
            by_aro, by_protein = card_cache[path_str]
            aro = _safe_text(row.get("SubjectAROAccession"))
            protein = _safe_text(row.get("SubjectProteinAccession"))
            payload = by_aro.get(aro) or by_protein.get(protein)

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
    "parse_marker_subject_metadata",
    "enrich_marker_hits_with_resistance",
    "build_annotation_resistance_maps",
]
