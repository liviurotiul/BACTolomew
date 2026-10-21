from __future__ import annotations
from collections import defaultdict
import re
import time
import json
import logging
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .constants import (
    EXPECTED_ANNOTATION_COLUMNS,
    MARKER_SEARCH_TEXT_COLUMN,
    MISSING_FILTER_TOKEN,
    PLASMID_TRUE_VALUES,
    RUN_CONFIG_FILENAME,
    RUN_REPORT_FILENAME,
)
from .run_expectations import load_run_report, load_run_status, run_allows_view
from .progress import log_done, log_step
from .polars_io import read_csv_df
from pipeline.RunValidation import RESULT_REQUIRED_STEPS, skippable_steps, validate_completed_run
from pipeline.AmrFinderWorkflow import attach_resistance_mutations
from pipeline.AnnotationUtils import (
    assign_annotation_ids,
    infer_annotation_category_and_source,
)
from pipeline.BaktaEvidence import (
    SOURCE_AMRFINDER,
    contig_key,
    evidence_marker_hits,
    read_evidence,
    unplaced_evidence,
    unplaced_locus,
)
from pipeline.BaktaWorkflow import parse_multiple_bakta_json
from pipeline.AssemblyQualityWorkflow import ASSEMBLY_QC_COLUMNS, load_assembly_quality_summary
from pipeline.DbCatalog import (
    MARKER_DB_CATALOG_PATH,
    COMPOUND_DB_CATALOG_PATH,
    load_db_catalog_index,
    normalize_resistance_metadata,
)
from pipeline.DiamondWorkflow import (
    collapse_diamond_annotations,
    load_diamond_tsvs,
    read_annotations,
)
from pipeline.MobSuiteWorkflow import load_mobsuite_contig_predictions
from pipeline.RunMetadata import RUN_METADATA_LABEL_COLUMNS, build_run_metadata
from pipeline.ResistanceMetadata import (
    RESISTANCE_MARKER_COLUMNS,
    SUBJECT_METADATA_COLUMNS,
    build_annotation_resistance_maps,
    drop_resistance_target_hits,
    enrich_marker_hits_with_resistance,
)
from pipeline.FastaUtils import subseq
from pipeline.UserDatabases import run_user_databases
from pipeline.Utils import CATEGORIES, assign_gene_ids

# Metadata preparation guidelines:
# - Keep the required sample identifier column named exactly `Sample`.
# - Keep one metadata row per sample; do not duplicate `Sample` values.
# - Use empty cells for missing values. Do not mix placeholders such as
#   `unknown`, `NA`, `null`, `-`, or `?`.
# - Keep categorical identifiers as strings even if they look numeric.
#   In this project that includes values such as `subtype`, `k_locus`,
#   and `complex_type`.
# - Keep `collection_date` as text unless the whole column is a true
#   machine-readable date series with consistent precision. Year-only
#   values should stay as strings.
# - Prefer stable, lowercase snake_case column names for new metadata
#   files so they can be matched predictably in the GUI.
_ST492_METADATA_PATH = Path("/Data/ST492/metadata/metadata_ALL_new_merged_with_clusters_with_complex.csv")
_ST492_METADATA_TYPE_OVERRIDES: Dict[str, str] = {
    "Sample": "string",
    "subtype": "string",
    "isolation_source": "string",
    "collection_date": "string",
    "country": "string",
    "city": "string",
    "host": "string",
    "disease": "string",
    "k_locus": "string",
    "complex_type": "string",
}


def _to_int(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _make_ann_uid(row: pd.Series) -> str:
    parts = [
        str(row.get("Sample") or "").strip(),
        str(row.get("Contig") or "").strip(),
        str(row.get("StartInt") if pd.notna(row.get("StartInt")) else row.get("Start") or "").strip(),
        str(row.get("StopInt") if pd.notna(row.get("StopInt")) else row.get("Stop") or "").strip(),
        str(row.get("Strand") or "").strip(),
    ]
    for fb in [row.get("AnnotationID"), row.get("GeneID"), row.get("Gene"), row.get("Product"), row.get("Id")]:
        if pd.notna(fb) and str(fb).strip():
            parts.append(str(fb).strip())
            break
    else:
        parts.append("ann")
    return "|".join(parts)


def _build_marker_map(df: pd.DataFrame, marker_cols: List[str]) -> Dict[str, List[Dict[str, Any]]]:
    if df.empty or "AnnotationID" not in df.columns:
        return {}
    ann = df["AnnotationID"]
    mask = ann.notna()
    if not mask.any():
        return {}
    sub = df.loc[mask, ["AnnotationID"] + marker_cols]
    records = sub[marker_cols].to_dict("records")
    ann_ids = sub["AnnotationID"].tolist()
    out: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for ann_id, rec in zip(ann_ids, records):
        out[ann_id].append(rec)
    return dict(out)


def _join_unique(values: Sequence[Any]) -> str:
    seen: set[str] = set()
    out: List[str] = []
    for value in values:
        if value is None or (isinstance(value, float) and pd.isna(value)):
            continue
        text = str(value).strip()
        if not text or text.lower() == "nan" or text in seen:
            continue
        seen.add(text)
        out.append(text)
    return ",".join(out)


MARKER_SEARCH_TEXT_COLUMNS: Tuple[str, ...] = (
    "Product",
    "DB",
    "Method",
    "MarkerCategory",
    "Category",
    "SubjectIdentifier",
    "SubjectProteinAccession",
    "SubjectDNAAccession",
    "SubjectAROAccession",
    "SubjectModelName",
    *RESISTANCE_MARKER_COLUMNS,
)


def build_marker_search_text_map(df: pd.DataFrame) -> Dict[str, str]:
    if df.empty or "AnnotationID" not in df.columns:
        return {}
    search_cols = [col for col in MARKER_SEARCH_TEXT_COLUMNS if col in df.columns]
    if not search_cols:
        return {}

    # This text is only ever used for lowercase substring search, so token order
    # and duplication don't matter -- only the set of tokens per AnnotationID.
    # A native polars group-by is ~4x faster than the old pandas row-wise
    # apply + python group agg (which dominated window load on large runs).
    import polars as pl

    pdf = pl.from_pandas(df[["AnnotationID", *search_cols]])
    pdf = pdf.with_columns(
        [pl.col(col).cast(pl.Utf8, strict=False) for col in ["AnnotationID", *search_cols]]
    )
    pdf = pdf.with_columns(pl.col("AnnotationID").str.strip_chars()).filter(
        pl.col("AnnotationID") != ""
    )
    if pdf.height == 0:
        return {}

    long = pdf.unpivot(index="AnnotationID", on=search_cols, value_name="val")
    long = long.with_columns(pl.col("val").str.strip_chars()).filter(
        pl.col("val").is_not_null()
        & (pl.col("val") != "")
        & (pl.col("val").str.to_lowercase() != "nan")
    )
    grouped = long.group_by("AnnotationID").agg(
        pl.col("val").unique(maintain_order=True).str.join(",")
    )
    return {key: value for key, value in grouped.iter_rows() if value}


def _finalize_marker_hits_df(df: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "AnnotationID",
        "Sample",
        "Contig",
        "Start",
        "Stop",
        "Strand",
        "Gene",
        "GeneID",
        "Method",
        "DB",
        "MarkerCategory",
        "Category",
        "Product",
        "Identity",
        "Align_len",
        "Mismatches",
        "Gap_openings",
        "Q_start",
        "Q_end",
        "S_start",
        "S_end",
        "Evalue",
        "Bitscore",
        "Q_len",
        "S_len",
        "Subject_coverage",
        *SUBJECT_METADATA_COLUMNS,
        *RESISTANCE_MARKER_COLUMNS,
    ]
    out = df.copy()
    for col in cols:
        if col not in out.columns:
            out[col] = pd.NA
    return out[cols]


# The annotation tags of marker categories ("core-AMR", "plus-VIRULENCE", ...).
_CATEGORY_TAGS = (("AMR", "A"), ("BMRG", "B"), ("VIRULENCE", "V"), ("MOBILE", "M"), ("PHAGE", "P"))


def category_tags(categories: str) -> str:
    """The tag letters (A, B, V, M, P) of a text naming marker categories."""
    text = str(categories or "").upper()
    return "".join(tag for name, tag in _CATEGORY_TAGS if name in text)


def _finalize_annotations_df(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "Sample" in df.columns:
        df["Sample"] = df["Sample"].astype(str).str.strip()
    else:
        df["Sample"] = pd.NA
    if "Contig" in df.columns:
        df["Contig"] = df["Contig"].astype(str).str.strip()
    else:
        df["Contig"] = pd.NA

    df["StartInt"] = _to_int(df.get("Start"))
    df["StopInt"] = _to_int(df.get("Stop"))
    df["AnnUID"] = df.apply(_make_ann_uid, axis=1)

    if "MarkerDisplay" in df.columns:
        df["MarkerDisplay"] = df["MarkerDisplay"].replace({"": pd.NA})
    else:
        df["MarkerDisplay"] = pd.NA
    # fall back to existing MarkerCategory if parsed is empty
    if "MarkerCategory" in df.columns:
        df["MarkerDisplay"] = df["MarkerDisplay"].fillna(df["MarkerCategory"])

    if "Predictors" in df.columns:
        df["Predictors"] = df["Predictors"].replace({"": pd.NA})
    else:
        df["Predictors"] = pd.NA

    for col in ("ConsensusCategory", "MarkerCategory"):
        if col not in df.columns:
            df[col] = pd.NA
    cat_series = (
        df["ConsensusCategory"].fillna("").astype(str)
        + " "
        + df["MarkerCategory"].fillna("").astype(str)
    ).str.upper()

    df["Tags"] = cat_series.apply(category_tags)
    # X: a known resistance mutation (AMRFinderPlus POINT) in this CDS.
    if "ResistanceMutations" in df.columns:
        for col in ("ResistanceMutations", "ResistanceMutationClasses"):
            df[col] = df[col].fillna("").astype(str) if col in df.columns else ""
        mutated = df["ResistanceMutations"].str.strip() != ""
        df.loc[mutated, "Tags"] = df.loc[mutated, "Tags"] + "X"

    needed_cols = [
        *EXPECTED_ANNOTATION_COLUMNS,
        "AnnUID",
        "MarkerDisplay",
        "StartInt",
        "StopInt",
    ]
    for col in needed_cols:
        if col not in df.columns:
            df[col] = pd.NA
    return df[needed_cols]


def load_annotations(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Annotations CSV not found: {path}")
    df = read_csv_df(path, dtype=str)
    return _finalize_annotations_df(df)


def load_metadata(path: Path, known_samples: Sequence[str]) -> pd.DataFrame:
    df, _warnings, _types = prepare_metadata(path, known_samples, [])
    return df


def load_reference_locus_samples(path: Path) -> List[str]:
    if not path.exists():
        raise FileNotFoundError(f"Reference locus CSV not found: {path}")
    t0 = time.perf_counter()
    print(f"[file] reading {path}")
    df = read_csv_df(path, usecols=["Sample"], dtype=str)
    print(f"[file] read {path} in {time.perf_counter() - t0:.2f}s")
    samples = df["Sample"].astype(str).str.strip()
    samples = samples[samples != ""]
    return samples.dropna().unique().tolist()


def prepare_metadata(
    path: Path,
    annotation_samples: Sequence[str],
    reference_samples: Sequence[str],
) -> Tuple[pd.DataFrame, List[str], Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Metadata CSV not found: {path}")
    t0 = time.perf_counter()
    print(f"[file] reading {path}")
    df = read_csv_df(path, dtype=str)
    print(f"[file] read {path} in {time.perf_counter() - t0:.2f}s")
    if "Sample" not in df.columns:
        raise ValueError("Metadata CSV is missing required column: Sample")

    warnings: List[str] = []
    df = df.copy()
    df["Sample"] = df["Sample"].astype(str).str.strip()
    df.loc[df["Sample"] == "", "Sample"] = pd.NA

    empty_sample_rows = df["Sample"].isna().sum()
    if empty_sample_rows:
        warnings.append(f"Metadata has {int(empty_sample_rows)} row(s) with empty Sample; they were dropped.")
        df = df[df["Sample"].notna()].copy()

    ann_samples = _normalize_samples(annotation_samples)
    ref_samples = _normalize_samples(reference_samples)
    target_samples = sorted(set(ann_samples) | set(ref_samples))
    meta_samples = set(df["Sample"].dropna().unique().tolist())

    missing_in_meta = [s for s in target_samples if s not in meta_samples]
    extra_in_meta = [s for s in meta_samples if s not in set(target_samples)]

    if missing_in_meta:
        warnings.append(
            "Metadata missing samples present in data: "
            f"{_format_sample_list(missing_in_meta)}. Added empty rows for them."
        )
        df = pd.concat([df, pd.DataFrame({"Sample": missing_in_meta})], ignore_index=True)

    if target_samples:
        df = df[df["Sample"].isin(target_samples)].copy()

    if extra_in_meta:
        warnings.append(
            "Metadata has samples not found in annotations/reference locus data: "
            f"{_format_sample_list(extra_in_meta)}. These were ignored."
        )

    df = _coerce_empty_to_na(df)
    df, inferred_types = _infer_metadata_types(df, overrides=_metadata_type_overrides_for_path(path))
    df.attrs["inferred_types"] = inferred_types
    return df, warnings, inferred_types


def prepare_run_metadata(
    output_dir: Path,
    annotation_samples: Sequence[str],
    reference_samples: Sequence[str],
) -> Tuple[pd.DataFrame, List[str], Dict[str, str]]:
    """``prepare_metadata`` for a run opened without a metadata CSV: the run's
    own sample table (species, MLST, fastGEAR, assembly QC) for every genome
    in the data, typed like a CSV."""
    samples = sorted(set(_normalize_samples(annotation_samples)) | set(_normalize_samples(reference_samples)))
    df = build_run_metadata(Path(output_dir), samples)
    df = df.replace("", pd.NA)
    overrides = {col: "string" for col in RUN_METADATA_LABEL_COLUMNS if col in df.columns}
    df, inferred_types = _infer_metadata_types(df, overrides=overrides)
    df.attrs["inferred_types"] = inferred_types
    return df, [], inferred_types


def attach_assembly_quality_metadata(meta_df: pd.DataFrame, output_dir: Optional[Path]) -> pd.DataFrame:
    if output_dir is None or meta_df.empty or "Sample" not in meta_df.columns:
        return meta_df
    summary = load_assembly_quality_summary(Path(output_dir))
    if summary.empty or "Sample" not in summary.columns or "Quality" not in summary.columns:
        return meta_df

    qc_columns = [col for col in ASSEMBLY_QC_COLUMNS if col in summary.columns]
    qc = summary[qc_columns].copy()
    qc["Sample"] = qc["Sample"].astype(str).str.strip()
    qc = qc[qc["Sample"] != ""].drop_duplicates(subset=["Sample"], keep="first")
    integer_cols = {
        "Quality",
        "QuastTotalLength",
        "QuastContigs",
        "QuastLargestContig",
        "QuastN50",
        "QuastL50",
        "QuastNs",
    }
    float_cols = {
        "QuastGC",
        "QuastNsPer100kbp",
        "CheckM2Completeness",
        "CheckM2Contamination",
    }
    for col in integer_cols:
        if col in qc.columns:
            qc[col] = pd.to_numeric(qc[col], errors="coerce").astype("Int64")
    for col in float_cols:
        if col in qc.columns:
            qc[col] = pd.to_numeric(qc[col], errors="coerce")

    out = meta_df.copy()
    attrs = dict(getattr(out, "attrs", {}) or {})
    inferred_types = dict(attrs.get("inferred_types") or {})
    drop_cols = [col for col in qc.columns if col != "Sample" and col in out.columns]
    if drop_cols:
        out = out.drop(columns=drop_cols)
    out = out.merge(qc, on="Sample", how="left")
    qc_non_sample = [col for col in qc.columns if col != "Sample"]
    ordered = ["Sample", *qc_non_sample, *[col for col in out.columns if col not in {"Sample", *qc_non_sample}]]
    out = out[ordered]
    for col in integer_cols & set(qc.columns):
        inferred_types[col] = "int"
    for col in float_cols & set(qc.columns):
        inferred_types[col] = "float"
    for col in set(qc.columns) - integer_cols - float_cols - {"Sample"}:
        inferred_types[col] = "string"
    out.attrs["assembly_qc_columns"] = qc_non_sample
    out.attrs["inferred_types"] = inferred_types
    return out


def _normalize_samples(samples: Sequence[str]) -> List[str]:
    out: List[str] = []
    for s in samples:
        if s is None:
            continue
        val = str(s).strip()
        if not val or val.lower() == "nan":
            continue
        if val:
            out.append(val)
    return out


def _format_sample_list(samples: Sequence[str], limit: int = 10) -> str:
    samples = list(samples)
    if len(samples) <= limit:
        return ", ".join(samples)
    return ", ".join(samples[:limit]) + f", … (+{len(samples) - limit} more)"


def _coerce_empty_to_na(df: pd.DataFrame) -> pd.DataFrame:
    for col in df.columns:
        if col == "Sample":
            continue
        df[col] = df[col].replace(r"^\\s*$", pd.NA, regex=True)
    return df


def _metadata_type_overrides_for_path(path: Path) -> Dict[str, str]:
    try:
        normalized = path.expanduser().resolve()
    except Exception:
        normalized = path.expanduser().absolute()
    try:
        st492_path = _ST492_METADATA_PATH.expanduser().resolve()
    except Exception:
        st492_path = _ST492_METADATA_PATH.expanduser().absolute()
    if normalized == st492_path:
        return dict(_ST492_METADATA_TYPE_OVERRIDES)
    return {}


def _coerce_series_to_declared_type(series: pd.Series, declared_type: str) -> pd.Series:
    s = series.copy()
    s = s.replace(r"^\\s*$", pd.NA, regex=True)
    if declared_type == "string":
        return s.astype("string")
    if declared_type == "int":
        return pd.to_numeric(s, errors="coerce").astype("Int64")
    if declared_type == "float":
        return pd.to_numeric(s, errors="coerce").astype("float")
    if declared_type == "date":
        return pd.to_datetime(s, errors="coerce")
    raise ValueError(f"Unsupported metadata declared type: {declared_type}")


def _infer_metadata_types(
    df: pd.DataFrame,
    *,
    overrides: Optional[Dict[str, str]] = None,
) -> Tuple[pd.DataFrame, Dict[str, str]]:
    inferred: Dict[str, str] = {}
    overrides = overrides or {}
    for col in df.columns:
        declared_type = overrides.get(col)
        if declared_type is not None:
            df[col] = _coerce_series_to_declared_type(df[col], declared_type)
            inferred[col] = declared_type
            continue
        if col == "Sample":
            df[col] = df[col].astype(str).str.strip()
            inferred[col] = "string"
            continue
        df[col], inferred[col] = _infer_series_type(df[col])
    return df, inferred


def _infer_series_type(series: pd.Series) -> Tuple[pd.Series, str]:
    s = series.copy()
    s = s.replace(r"^\\s*$", pd.NA, regex=True)
    non_null = s.dropna()
    if non_null.empty:
        return s.astype("string"), "string"

    lowered = non_null.astype(str).str.strip().str.lower()
    bool_values = {"true", "false", "yes", "no", "0", "1"}
    if lowered.isin(bool_values).all():
        mapped = lowered.map(
            {
                "true": True,
                "yes": True,
                "1": True,
                "false": False,
                "no": False,
                "0": False,
            }
        )
        return mapped.reindex(s.index), "bool"

    numeric = pd.to_numeric(non_null, errors="coerce")
    if numeric.notna().all():
        numeric_full = pd.to_numeric(s, errors="coerce")
        if (numeric.dropna() % 1 == 0).all():
            return numeric_full.astype("Int64"), "int"
        return numeric_full.astype("float"), "float"

    if _looks_like_date(non_null):
        parsed = pd.to_datetime(non_null, errors="coerce")
        if parsed.notna().mean() >= 0.9:
            return pd.to_datetime(s, errors="coerce"), "date"

    return s.astype("string"), "string"


def _looks_like_date(values: pd.Series) -> bool:
    s = values.astype(str).str.strip()
    has_sep = s.str.contains(r"[-/]", regex=True)
    is_yyyymmdd = s.str.match(r"^\\d{8}$")
    return (has_sep | is_yyyymmdd).mean() >= 0.8


def _contains(series: pd.Series, text: Optional[str]) -> pd.Series:
    if not text:
        return pd.Series(True, index=series.index)
    text = text.strip()
    if not text:
        return pd.Series(True, index=series.index)
    return series.astype(str).str.contains(text, case=False, na=False)


def _contains_any(df: pd.DataFrame, text: Optional[str]) -> pd.Series:
    if not text:
        return pd.Series(True, index=df.index)
    text = text.strip()
    if not text:
        return pd.Series(True, index=df.index)
    mask = pd.Series(False, index=df.index)
    for col in df.columns:
        mask = mask | df[col].astype(str).str.contains(text, case=False, na=False)
    return mask


def apply_metadata_filters(
    df: pd.DataFrame,
    sample_query: Optional[str],
    column_filters: Dict[str, Sequence[str]],
) -> pd.DataFrame:
    out = df.copy()
    if sample_query:
        out = out[_contains_any(out, sample_query)]
    for col, raw_values in column_filters.items():
        if col not in out.columns or not raw_values:
            continue
        values = {str(v) for v in raw_values if v is not None}
        include_missing = MISSING_FILTER_TOKEN in values
        values.discard(MISSING_FILTER_TOKEN)
        if include_missing and values:
            mask = out[col].isna() | out[col].astype(str).isin(values)
        elif include_missing:
            mask = out[col].isna()
        elif values:
            mask = out[col].astype(str).isin(values)
        else:
            continue
        out = out[mask]
    return out


def apply_annotation_filters(df: pd.DataFrame, filters: Dict[str, Any]) -> pd.DataFrame:
    out = df
    if filters.get("predictors"):
        preds = [str(p).strip() for p in filters["predictors"] if p]
        if preds and "Predictors" in out.columns:
            pattern = r"(^|,\\s*)(" + "|".join(map(re.escape, preds)) + r")(,|$)"
            out = out[out["Predictors"].astype(str).str.contains(pattern, regex=True, na=False)]
    if filters.get("tags"):
        tags = [str(t).upper() for t in filters["tags"] if t]
        if tags:
            pattern = "[" + "".join(tags) + "]"
            out = out[out["Tags"].astype(str).str.upper().str.contains(pattern, regex=True, na=False)]
    if filters.get("gene_text"):
        query = filters["gene_text"]
        if "Gene" in out.columns:
            mask = _contains(out["Gene"], query)
        else:
            mask = pd.Series(False, index=out.index)
        for col in (MARKER_SEARCH_TEXT_COLUMN, "MarkerDisplay", "Predictors", "Method"):
            if col in out.columns:
                mask = mask | _contains(out[col], query)
        out = out[mask]
    if filters.get("product_text"):
        out = out[_contains(out["Product"], filters["product_text"])]
    if filters.get("hide_hypothetical"):
        out = out[~out["Hypothetical"].astype(str).str.lower().isin(["true", "1", "yes"])]
    molecule = str(filters.get("molecule") or "").strip().lower()
    if molecule in {"plasmid", "chromosome"}:
        if "MobSuiteIsPlasmid" in out.columns:
            plasmid_mask = (
                out["MobSuiteIsPlasmid"]
                .astype(str)
                .str.strip()
                .str.lower()
                .isin(PLASMID_TRUE_VALUES)
            )
        else:
            plasmid_mask = pd.Series(False, index=out.index)
        out = out[plasmid_mask] if molecule == "plasmid" else out[~plasmid_mask]
    return out


def _read_json_file(path: Path) -> Dict[str, Any]:
    start = time.perf_counter()
    print(f"[file] reading {path}")
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise ValueError(f"Failed to read {path.name}: {exc}")
    if not isinstance(data, dict):
        raise ValueError(f"{path.name} must be a JSON object.")
    print(f"[file] read {path} in {time.perf_counter() - start:.2f}s")
    return data


def load_run_config(output_dir: Path) -> Dict[str, Any]:
    path = Path(output_dir) / RUN_CONFIG_FILENAME
    if not path.exists() or not path.is_file():
        raise FileNotFoundError(f"Required file not found: {path}")
    return _read_json_file(path)


def validate_output_dir(
    output_dir: Path, *, deep: bool = True
) -> Tuple[Optional[Dict[str, Any]], List[str]]:
    # deep=False skips validate_completed_run, which scans every expected output
    # file across all samples/steps (seconds on large runs). The fast JSON checks
    # (finalized status, run_report, required steps) always run.
    errors: List[str] = []
    output_dir = Path(output_dir)
    if not output_dir.exists() or not output_dir.is_dir():
        errors.append("Output folder does not exist or is not a directory.")
        return None, errors

    run_cfg = None
    run_cfg_path = output_dir / RUN_CONFIG_FILENAME
    if not run_cfg_path.exists():
        errors.append(f"Required file not found: {RUN_CONFIG_FILENAME}")
        return None, errors
    try:
        run_cfg = _read_json_file(run_cfg_path)
    except Exception as exc:
        errors.append(str(exc))
        return None, errors

    status_info = load_run_status(output_dir)
    if not run_allows_view(status_info):
        status = str(status_info.get("status") or "").strip() or "unknown"
        detail = str(status_info.get("error") or "").strip()
        msg = f"Run is not finalized successfully (status: {status})."
        if detail:
            msg += f" {detail}"
        errors.append(msg)
        return run_cfg, errors

    run_report = load_run_report(output_dir)
    if not isinstance(run_report, dict):
        errors.append(f"Required file not found or invalid: {RUN_REPORT_FILENAME}")
        return run_cfg, errors
    if str(run_report.get("status") or "").strip().lower() != "success":
        errors.append(f"{RUN_REPORT_FILENAME} is not finalized successfully.")
        return run_cfg, errors
    if not bool(run_report.get("outputs_validated")):
        errors.append(f"{RUN_REPORT_FILENAME} says outputs were not validated successfully.")
        return run_cfg, errors

    inputs = run_cfg.get("inputs") if isinstance(run_cfg, dict) else None
    if not isinstance(inputs, dict):
        errors.append("run_config.json is missing the 'inputs' object.")
        return run_cfg, errors
    marker_inputs = inputs.get("marker_sequences")
    locus_inputs = inputs.get("reference_locus_sequences")
    # An empty list is a run without any installed marker (or compound) database.
    if not isinstance(marker_inputs, list):
        errors.append("run_config.json is missing inputs.marker_sequences.")
    if not isinstance(locus_inputs, list):
        errors.append("run_config.json is missing inputs.reference_locus_sequences.")
    report_steps = run_report.get("steps")
    if not isinstance(report_steps, dict):
        errors.append(f"{RUN_REPORT_FILENAME} is missing the 'steps' object.")
        return run_cfg, errors
    skippable = skippable_steps(run_cfg)
    for step_name in RESULT_REQUIRED_STEPS:
        step_info = report_steps.get(step_name)
        step_status = str(step_info.get("status") or "").strip().lower() if isinstance(step_info, dict) else ""
        if step_status == "skipped" and step_name in skippable:
            continue
        if step_status != "success":
            errors.append(f"Required step is not successful: {step_name} ({step_status or 'missing'}).")

    if deep:
        validations, validation_errors = validate_completed_run(
            output_dir,
            run_cfg=run_cfg,
            steps=report_steps,
        )
        for step_name, info in validations.items():
            if not bool(info.get("outputs_validated")):
                step_errors = info.get("validation_errors") or []
                if isinstance(step_errors, list):
                    errors.extend(str(item) for item in step_errors if str(item).strip())
                else:
                    errors.append(f"Run validation failed for step: {step_name}")
        errors.extend(str(item) for item in validation_errors if str(item).strip())

    deduped_errors: List[str] = []
    seen_errors: set[str] = set()
    for error in errors:
        text = str(error).strip()
        if not text or text in seen_errors:
            continue
        seen_errors.add(text)
        deduped_errors.append(text)

    return run_cfg, deduped_errors


def _get_inputs_list(run_cfg: Dict[str, Any], key: str) -> List[str]:
    inputs = run_cfg.get("inputs")
    if not isinstance(inputs, dict):
        raise ValueError("run_config.json is missing the 'inputs' object.")
    vals = inputs.get(key)
    if not isinstance(vals, list):
        raise ValueError(f"run_config.json is missing inputs.{key}.")
    return [str(v) for v in vals]


def _db_name_from_tsv(path: Path) -> str:
    stem = path.stem
    if "-vs-" not in stem:
        raise ValueError(f"Unexpected DIAMOND TSV name: {path.name}")
    return stem.split("-vs-", 1)[1]


def _run_hit_tables(
    tsv_files: List[Path], run_dbs: set, *, logger: logging.Logger, label: str
) -> List[Path]:
    """The hit tables of the databases the run searched; a table left by a database
    removed since (a folder resumed with fewer databases) is not read."""
    kept = [p for p in tsv_files if _db_name_from_tsv(p) in run_dbs]
    stale = sorted({_db_name_from_tsv(p) for p in tsv_files} - run_dbs)
    if stale:
        logger.warning("Ignoring %s hit tables of databases not in this run: %s", label, ", ".join(stale))
    return kept


def _ensure_annotation_columns(df: pd.DataFrame) -> pd.DataFrame:
    for col in EXPECTED_ANNOTATION_COLUMNS:
        if col not in df.columns:
            df[col] = pd.NA
    return df


def _normalize_contig_for_mobsuite(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    if text.lower().startswith("contig_"):
        text = text.split("_", 1)[1]
    try:
        fval = float(text)
    except Exception:
        return text
    if fval.is_integer():
        return str(int(fval))
    return text


def _attach_mobsuite_predictions(df: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    predictions = load_mobsuite_contig_predictions(output_dir)
    for col in (
        "MobSuitePrediction",
        "MobSuiteIsPlasmid",
        "MobSuiteSourceContig",
        "MobSuitePlasmidID",
        "MobSuiteMoleculeType",
        "MobSuitePrimaryClusterID",
        "MobSuiteSecondaryClusterID",
        "MobSuiteLength",
        "MobSuiteRepliconTypes",
        "MobSuiteRelaxaseTypes",
        "MobSuiteMPFType",
        "MobSuiteOriTTypes",
        "MobSuitePredictedMobility",
        "MobSuiteMashNearestNeighbor",
        "MobSuiteMashNeighborDistance",
    ):
        if col not in df.columns:
            df[col] = pd.NA
    if not predictions or df.empty:
        if "MobSuiteIsPlasmid" in df.columns:
            df["MobSuiteIsPlasmid"] = df["MobSuiteIsPlasmid"].fillna("false")
        return df

    def _lookup(row: pd.Series) -> dict[str, str]:
        sample = str(row.get("Sample") or "").strip()
        contig = _normalize_contig_for_mobsuite(row.get("Contig"))
        return predictions.get((sample, contig), {})

    hits = df.apply(_lookup, axis=1)
    df["MobSuitePrediction"] = hits.apply(lambda item: item.get("Prediction") or pd.NA)
    df["MobSuiteIsPlasmid"] = hits.apply(lambda item: item.get("IsPlasmid") or "false")
    df["MobSuiteSourceContig"] = hits.apply(lambda item: item.get("SourceContig") or pd.NA)
    df["MobSuitePlasmidID"] = hits.apply(lambda item: item.get("PlasmidID") or pd.NA)
    df["MobSuiteMoleculeType"] = hits.apply(lambda item: item.get("MoleculeType") or pd.NA)
    df["MobSuitePrimaryClusterID"] = hits.apply(lambda item: item.get("PrimaryClusterID") or pd.NA)
    df["MobSuiteSecondaryClusterID"] = hits.apply(lambda item: item.get("SecondaryClusterID") or pd.NA)
    df["MobSuiteLength"] = hits.apply(lambda item: item.get("Size") or pd.NA)
    df["MobSuiteRepliconTypes"] = hits.apply(lambda item: item.get("RepliconTypes") or pd.NA)
    df["MobSuiteRelaxaseTypes"] = hits.apply(lambda item: item.get("RelaxaseTypes") or pd.NA)
    df["MobSuiteMPFType"] = hits.apply(lambda item: item.get("MPFType") or pd.NA)
    df["MobSuiteOriTTypes"] = hits.apply(lambda item: item.get("OriTTypes") or pd.NA)
    df["MobSuitePredictedMobility"] = hits.apply(lambda item: item.get("PredictedMobility") or pd.NA)
    df["MobSuiteMashNearestNeighbor"] = hits.apply(lambda item: item.get("MashNearestNeighbor") or pd.NA)
    df["MobSuiteMashNeighborDistance"] = hits.apply(lambda item: item.get("MashNeighborDistance") or pd.NA)
    return df


def _build_db_entry_map(
    db_paths: List[str],
    catalog_path: Path,
    *,
    logger: logging.Logger,
    label: str,
    recorded_categories: Optional[Dict[str, Any]] = None,
) -> Dict[str, Dict[str, Any]]:
    """The catalog entry of each database the run searched; a database the catalog no
    longer lists (or no catalog) keeps the category the run recorded for it."""
    recorded = recorded_categories if isinstance(recorded_categories, dict) else {}
    index = None
    if db_paths:
        try:
            index = load_db_catalog_index(catalog_path, logger=logger)
        except FileNotFoundError:
            index = None
    out: Dict[str, Dict[str, Any]] = {}
    missing: List[str] = []
    for p_str in db_paths:
        p = Path(p_str).expanduser()
        entry = index.entry_for_path(p) if index is not None else None
        if entry is None and p.stem in recorded:
            entry = {"category": recorded[p.stem], "db_path": str(p)}
        if entry is None:
            missing.append(str(p))
            continue
        stem = p.stem
        if stem in out:
            raise ValueError(f"Duplicate {label} DB stem in run_config: {stem}")
        out[stem] = dict(entry)
    if missing:
        raise FileNotFoundError(
            f"Missing {label} catalog entries for: {', '.join(missing)}"
        )
    return out


def _build_db_category_map(
    db_paths: List[str],
    catalog_path: Path,
    *,
    logger: logging.Logger,
    label: str,
    recorded_categories: Optional[Dict[str, Any]] = None,
) -> Dict[str, str]:
    entry_map = _build_db_entry_map(
        db_paths,
        catalog_path,
        logger=logger,
        label=label,
        recorded_categories=recorded_categories,
    )
    return {stem: str(entry.get("category") or "").strip() for stem, entry in entry_map.items()}


def _build_db_resistance_metadata_map(
    run_cfg: Dict[str, Any],
    db_entry_map: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    inputs = run_cfg.get("inputs") if isinstance(run_cfg, dict) else None
    raw = inputs.get("marker_resistance_metadata") if isinstance(inputs, dict) else None

    out: Dict[str, Dict[str, Any]] = {}
    if isinstance(raw, dict):
        for stem, payload in raw.items():
            meta = payload.get("resistance_metadata") if isinstance(payload, dict) and "resistance_metadata" in payload else payload
            normalized = normalize_resistance_metadata(meta)
            if normalized is not None:
                out[str(stem)] = normalized

    for stem, entry in db_entry_map.items():
        if stem in out:
            continue
        normalized = normalize_resistance_metadata(entry.get("resistance_metadata"))
        if normalized is not None:
            out[stem] = normalized
    return out


# AnnotationIDs of the genes that are an AMRFinderPlus call on no Bakta CDS (after every
# Bakta row's BKT id, so those stay the ones the marker search used).
UNPLACED_ID_PREFIX = "AFP"
_UNPLACED_KEY_COLUMN = "__evidence_key"


def _with_unplaced_genes(bakta_df: pd.DataFrame, evidence: Optional[pd.DataFrame], output_dir: Path) -> pd.DataFrame:
    """The Bakta rows plus one row per AMRFinderPlus gene call that lies on no Bakta CDS
    (its contig coordinates, Method "AMRFinderPlus"), so it is found, tagged and drawn
    like any gene; the run's evidence keeps those calls with an empty locus tag."""
    calls = unplaced_evidence(evidence)
    if calls.empty or bakta_df.empty:
        return bakta_df
    rows: List[Dict[str, Any]] = []
    for call in calls.to_dict("records"):
        sample = str(call["sample"])
        contig = contig_key(call["contig"])
        try:
            start, stop = int(call["start"]), int(call["stop"])
        except (TypeError, ValueError):
            continue
        fasta = output_dir / "bakta" / sample / f"{sample}.fna"
        rows.append({
            "Sample": sample,
            "Contig": int(contig) if contig.isdigit() else pd.NA,
            "Start": start,
            "Stop": stop,
            "Strand": call.get("strand") or pd.NA,
            "Gene": call.get("gene") or pd.NA,
            "Product": call.get("product") or pd.NA,
            "Hypothetical": False,
            "Method": SOURCE_AMRFINDER,
            "Nucleotide_Sequence": subseq(fasta, str(call["contig"]), start, stop, call.get("strand") or "+")
            if fasta.is_file() else pd.NA,
            _UNPLACED_KEY_COLUMN: unplaced_locus(call),
        })
    if not rows:
        return bakta_df
    extra = pd.DataFrame(rows)
    extra["AnnotationID"] = [f"{UNPLACED_ID_PREFIX}{i:07d}" for i in range(1, len(extra) + 1)]
    return pd.concat([bakta_df, extra], ignore_index=True)


def _evidence_marker_hits(evidence: Optional[pd.DataFrame], bakta_df: pd.DataFrame) -> pd.DataFrame:
    """The run's gene evidence (evidence/*.tsv) as marker hits of the Bakta rows (and of
    the rows `_with_unplaced_genes` added); empty for a run made before the bakta_evidence
    step."""
    if evidence is None or bakta_df.empty:
        return pd.DataFrame()
    t_ev = log_step("Reading the gene evidence")
    locus = bakta_df["Locus"].astype(str)
    if _UNPLACED_KEY_COLUMN in bakta_df.columns:
        locus = bakta_df[_UNPLACED_KEY_COLUMN].fillna(locus).astype(str)
    keys = zip(bakta_df["Sample"].astype(str), locus)
    hits = evidence_marker_hits(evidence, dict(zip(keys, bakta_df["AnnotationID"])))
    log_done("Reading the gene evidence", t_ev, extra=f"{len(hits)} rows")
    return hits


def build_annotations_from_output(output_dir: Path, run_cfg: Dict[str, Any]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    logger = logging.getLogger("annotation_explorer")
    output_dir = Path(output_dir)

    annotation_id_strategy = str(run_cfg.get("annotation_id_strategy") or "").strip().lower()
    if annotation_id_strategy and annotation_id_strategy != "sorted_bakta_json_v1":
        raise ValueError(f"Unsupported annotation_id_strategy: {annotation_id_strategy}")
    bakta_json_files = sorted((output_dir / "bakta").rglob("*.json"), key=lambda p: str(p))
    if not bakta_json_files:
        raise FileNotFoundError("No Bakta JSON files found under output/bakta.")

    marker_paths = _get_inputs_list(run_cfg, "marker_sequences")
    marker_entry_map = _build_db_entry_map(
        marker_paths,
        MARKER_DB_CATALOG_PATH,
        logger=logger,
        label="marker",
        recorded_categories=run_cfg["inputs"].get("marker_categories"),
    )
    marker_cats = {stem: str(entry.get("category") or "").strip() for stem, entry in marker_entry_map.items()}
    marker_resistance_metadata = _build_db_resistance_metadata_map(run_cfg, marker_entry_map)

    tsv_dir = output_dir / "annotations" / "bakta_vs_queries"
    all_tsv_files = sorted(tsv_dir.rglob("*.tsv"))
    tsv_files = _run_hit_tables(all_tsv_files, set(marker_cats), logger=logger, label="marker")
    if not tsv_files and marker_paths:
        raise FileNotFoundError("No DIAMOND TSVs found under output/annotations/bakta_vs_queries.")

    t_bakta = log_step(f"Parsing Bakta JSONs ({len(bakta_json_files)})")
    bakta_df = parse_multiple_bakta_json(bakta_json_files, logger)
    log_done("Parsing Bakta JSONs", t_bakta, extra=f"{len(bakta_df)} rows")
    bakta_df = _ensure_annotation_columns(bakta_df)

    evidence = read_evidence(output_dir)
    if not bakta_df.empty:
        bakta_df["Method"] = "bakta"
        bakta_df = assign_annotation_ids(bakta_df, prefix="BKT")
        bakta_df = _with_unplaced_genes(bakta_df, evidence, output_dir)
        bakta_df = _ensure_annotation_columns(bakta_df)

    # A user database's hits keep its whole name as their Method.
    user_dbs = run_user_databases(run_cfg.get("inputs"))
    t_diamond = log_step(f"Reading marker DIAMOND TSVs ({len(tsv_files)})")
    if len(tsv_files) == len(all_tsv_files):
        diamond_annotations = read_annotations(tsv_dir, user_dbs)
    elif tsv_files:
        diamond_annotations = pd.concat([read_annotations(p, user_dbs) for p in tsv_files], ignore_index=True)
    else:  # only tables of databases no longer in the run are left
        diamond_annotations = pd.DataFrame()
    log_done("Reading marker DIAMOND TSVs", t_diamond, extra=f"{len(diamond_annotations)} hits")
    if not diamond_annotations.empty:
        diamond_annotations["MarkerCategory"] = diamond_annotations["DB"].map(marker_cats)
        missing_map = diamond_annotations["MarkerCategory"].isna() & diamond_annotations["DB"].notna()
        if missing_map.any():
            bad = sorted(diamond_annotations.loc[missing_map, "DB"].unique().tolist())
            raise ValueError(
                "Missing marker categories for DBs: " + ", ".join([str(b) for b in bad])
            )

        diamond_annotations["MarkerCategory"] = diamond_annotations["MarkerCategory"].where(
            diamond_annotations["MarkerCategory"].isin(CATEGORIES)
        )
        diamond_annotations["Category"] = diamond_annotations["MarkerCategory"]
        diamond_annotations = drop_resistance_target_hits(
            diamond_annotations,
            marker_resistance_metadata,
            logger=logger,
        )
        t_collapse = log_step("Collapsing marker DIAMOND hits")
        diamond_annotations = collapse_diamond_annotations(diamond_annotations)
        log_done("Collapsing marker DIAMOND hits", t_collapse, extra=f"{len(diamond_annotations)} rows")
        t_res = log_step("Attaching resistance metadata to marker hits")
        diamond_annotations = enrich_marker_hits_with_resistance(
            diamond_annotations,
            marker_resistance_metadata,
            logger=logger,
        )
        log_done("Attaching resistance metadata to marker hits", t_res, extra=f"{len(diamond_annotations)} rows")

    # A run with the bakta_evidence step: its rows join the marker hits (DB = their source).
    evidence_hits = _evidence_marker_hits(evidence, bakta_df)
    bakta_df = bakta_df.drop(columns=[_UNPLACED_KEY_COLUMN], errors="ignore")
    if not evidence_hits.empty:
        diamond_annotations = (
            pd.concat([diamond_annotations, evidence_hits], ignore_index=True)
            if not diamond_annotations.empty
            else evidence_hits
        )
    if not diamond_annotations.empty:
        marker_cols = [c for c in diamond_annotations.columns if c != "AnnotationID"]
        t_map = log_step("Building marker map (linear scan)")
        marker_map = _build_marker_map(diamond_annotations, marker_cols)
        log_done("Building marker map (linear scan)", t_map, extra=f"{len(marker_map)} annotations")
        marker_hits_df = diamond_annotations.copy()
    else:
        marker_map = {}
        marker_hits_df = pd.DataFrame()

    if not marker_hits_df.empty and not bakta_df.empty:
        join_cols = [
            c
            for c in ["AnnotationID", "Sample", "Contig", "Start", "Stop", "Strand", "Gene", "GeneID"]
            if c in bakta_df.columns
        ]
        if "AnnotationID" in join_cols:
            ann_lookup = bakta_df[join_cols].drop_duplicates(subset=["AnnotationID"])
            marker_hits_df = marker_hits_df.merge(ann_lookup, on="AnnotationID", how="left")
    marker_hits_df = _finalize_marker_hits_df(marker_hits_df)

    method_map: Dict[str, str] = {}
    predictor_map: Dict[str, str] = {}
    resistance_maps = {col: {} for col in RESISTANCE_MARKER_COLUMNS}
    if not marker_hits_df.empty and "AnnotationID" in marker_hits_df.columns:
        grouped = marker_hits_df.groupby("AnnotationID", sort=False, dropna=False)
        method_map = grouped["Method"].agg(_join_unique).to_dict()
        predictor_map = grouped["DB"].agg(_join_unique).to_dict()
        resistance_maps = build_annotation_resistance_maps(marker_hits_df)

    if not bakta_df.empty:
        if "AnnotationID" not in bakta_df.columns:
            raise ValueError("Bakta annotations missing AnnotationID; cannot attach markers.")
        t_attach = log_step("Attaching marker summaries to Bakta rows")
        bakta_df["MarkerDisplay"] = bakta_df["AnnotationID"].map(method_map).replace({"": pd.NA})
        bakta_df["Predictors"] = bakta_df["AnnotationID"].map(predictor_map).replace({"": pd.NA})
        bakta_df[MARKER_SEARCH_TEXT_COLUMN] = (
            bakta_df["AnnotationID"].map(build_marker_search_text_map(marker_hits_df)).replace({"": pd.NA})
        )
        for col, col_map in resistance_maps.items():
            bakta_df[col] = bakta_df["AnnotationID"].map(col_map).replace({"": pd.NA})
        log_done("Attaching marker summaries to Bakta rows", t_attach, extra=f"{len(bakta_df)} rows")

        markers_for_consensus = bakta_df["AnnotationID"].map(marker_map)
        markers_for_consensus = markers_for_consensus.apply(
            lambda v: v if isinstance(v, list) else []
        )
        t_consensus = log_step("Inferring consensus category/source")
        consensus = markers_for_consensus.apply(infer_annotation_category_and_source)
        bakta_df["ConsensusCategory"] = consensus.str[0]
        bakta_df["ConsensusSource"] = consensus.str[1]
        bakta_df["ConsensusDonorProduct"] = consensus.str[2]
        log_done("Inferring consensus category/source", t_consensus)

        t_gene_ids = log_step("Assigning gene IDs")
        bakta_df = assign_gene_ids(bakta_df)
        log_done("Assigning gene IDs", t_gene_ids, extra=f"{len(bakta_df)} rows")

        t_mobsuite = log_step("Attaching MOB-suite plasmid contig predictions")
        bakta_df = _attach_mobsuite_predictions(bakta_df, output_dir)
        plasmid_rows = int((bakta_df["MobSuiteIsPlasmid"].astype(str).str.lower() == "true").sum())
        log_done("Attaching MOB-suite plasmid contig predictions", t_mobsuite, extra=f"{plasmid_rows} annotation rows")

        t_mut = log_step("Attaching AMRFinderPlus point mutations")
        bakta_df = attach_resistance_mutations(bakta_df, output_dir)
        mutated = int((bakta_df["ResistanceMutations"] != "").sum())
        log_done("Attaching AMRFinderPlus point mutations", t_mut, extra=f"{mutated} CDS")

    t_final = log_step("Finalizing annotation table")
    bakta_df = _ensure_annotation_columns(bakta_df)
    bakta_df = _finalize_annotations_df(bakta_df)
    log_done("Finalizing annotation table", t_final, extra=f"{len(bakta_df)} rows")
    return bakta_df, marker_hits_df


def build_reference_locus_from_output(output_dir: Path, run_cfg: Dict[str, Any]) -> pd.DataFrame:
    logger = logging.getLogger("annotation_explorer")
    output_dir = Path(output_dir)

    tsv_dir = output_dir / "annotations" / "reference_loci"
    tsv_files = sorted(tsv_dir.rglob("*.tsv"))
    locus_paths = _get_inputs_list(run_cfg, "reference_locus_sequences")
    if not tsv_files and locus_paths:
        raise FileNotFoundError("No DIAMOND TSVs found under output/annotations/reference_loci.")

    locus_cats = _build_db_category_map(
        locus_paths,
        COMPOUND_DB_CATALOG_PATH,
        logger=logger,
        label="reference locus",
        recorded_categories=run_cfg["inputs"].get("reference_locus_categories"),
    )
    tsv_files = _run_hit_tables(tsv_files, set(locus_cats), logger=logger, label="reference locus")

    t_locus = log_step(f"Reading reference locus DIAMOND TSVs ({len(tsv_files)})")
    ref_df = load_diamond_tsvs(tsv_files, locus_cats, logger, run_user_databases(run_cfg.get("inputs")))
    log_done("Reading reference locus DIAMOND TSVs", t_locus, extra=f"{len(ref_df)} hits")
    return ref_df


def load_output_data(
    output_dir: Path,
    progress_callback: Optional[Callable[[int, int], None]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    run_cfg = load_run_config(output_dir)
    ann_df, marker_hits_df = build_annotations_from_output(output_dir, run_cfg)
    ref_df = build_reference_locus_from_output(output_dir, run_cfg)
    t_overlap = log_step("Computing compound overlap tags")
    ann_df = add_compound_overlap_tags(ann_df, ref_df, progress_callback=progress_callback)
    log_done("Computing compound overlap tags", t_overlap)
    return ann_df, ref_df, marker_hits_df


def add_compound_overlap_tags(
    ann_df: pd.DataFrame,
    comp_df: pd.DataFrame,
    progress_callback: Optional[Callable[[int, int], None]] = None,
) -> pd.DataFrame:
    if ann_df.empty or comp_df.empty:
        return ann_df

    out = ann_df.copy()

    for col in ("Sample", "Contig"):
        if col in out.columns:
            out[col] = out[col].astype(str).str.strip()
        else:
            out[col] = pd.NA

    comp = comp_df.copy()
    for col in ("Sample", "Contig"):
        if col in comp.columns:
            comp[col] = comp[col].astype(str).str.strip()
        else:
            comp[col] = pd.NA

    out["StartInt"] = _to_int(out.get("Start"))
    out["StopInt"] = _to_int(out.get("Stop"))
    comp["StartInt"] = _to_int(comp.get("Start"))
    comp["StopInt"] = _to_int(comp.get("Stop"))

    out["Tags"] = out.get("Tags", "").fillna("").astype(str)

    comp = comp.dropna(subset=["Sample", "Contig", "StartInt", "StopInt"])
    if comp.empty:
        return out

    groups = list(comp.groupby(["Sample", "Contig"], dropna=False))
    total_groups = len(groups)
    for group_idx, ((sample, contig), comp_g) in enumerate(groups, start=1):
        idx = out.index[(out["Sample"] == sample) & (out["Contig"] == contig)]
        if idx.empty:
            if progress_callback is not None and (group_idx % 25 == 0 or group_idx == total_groups):
                try:
                    progress_callback(group_idx, total_groups)
                except Exception:
                    pass
            continue

        comp_starts = comp_g["StartInt"].to_numpy(dtype=float)
        comp_stops = comp_g["StopInt"].to_numpy(dtype=float)
        comp_valid = ~np.isnan(comp_starts) & ~np.isnan(comp_stops)
        comp_starts = comp_starts[comp_valid]
        comp_stops = comp_stops[comp_valid]
        if comp_starts.size == 0:
            if progress_callback is not None and (group_idx % 25 == 0 or group_idx == total_groups):
                try:
                    progress_callback(group_idx, total_groups)
                except Exception:
                    pass
            continue

        comp_min = np.minimum(comp_starts, comp_stops)
        comp_max = np.maximum(comp_starts, comp_stops)
        order_comp = np.argsort(comp_min, kind="mergesort")
        merged: List[Tuple[float, float]] = []
        for s, e in zip(comp_min[order_comp], comp_max[order_comp]):
            if not merged:
                merged.append((float(s), float(e)))
                continue
            prev_s, prev_e = merged[-1]
            if s <= prev_e:
                merged[-1] = (prev_s, max(prev_e, float(e)))
            else:
                merged.append((float(s), float(e)))

        ann_starts = out.loc[idx, "StartInt"].to_numpy(dtype=float)
        ann_stops = out.loc[idx, "StopInt"].to_numpy(dtype=float)
        ann_valid = ~np.isnan(ann_starts) & ~np.isnan(ann_stops)
        if not ann_valid.any():
            if progress_callback is not None and (group_idx % 25 == 0 or group_idx == total_groups):
                try:
                    progress_callback(group_idx, total_groups)
                except Exception:
                    pass
            continue

        overlap = np.zeros(len(ann_starts), dtype=bool)
        valid_pos = np.where(ann_valid)[0]
        ann_min = np.minimum(ann_starts[valid_pos], ann_stops[valid_pos])
        ann_max = np.maximum(ann_starts[valid_pos], ann_stops[valid_pos])
        order_ann = np.argsort(ann_min, kind="mergesort")
        j = 0
        merged_count = len(merged)
        for ordered_idx in order_ann:
            pos = int(valid_pos[ordered_idx])
            a_s = float(ann_min[ordered_idx])
            a_e = float(ann_max[ordered_idx])
            while j < merged_count and merged[j][1] < a_s:
                j += 1
            if j < merged_count and merged[j][0] <= a_e:
                overlap[pos] = True

        hit_idx = idx[overlap & ann_valid]
        if hit_idx.empty:
            if progress_callback is not None and (group_idx % 25 == 0 or group_idx == total_groups):
                try:
                    progress_callback(group_idx, total_groups)
                except Exception:
                    pass
            continue

        def _add_i_tag(s: str) -> str:
            s = s or ""
            return s if "I" in s else s + "I"

        out.loc[hit_idx, "Tags"] = out.loc[hit_idx, "Tags"].apply(_add_i_tag)

        if progress_callback is not None and (group_idx % 25 == 0 or group_idx == total_groups):
            try:
                progress_callback(group_idx, total_groups)
            except Exception:
                pass

    return out
