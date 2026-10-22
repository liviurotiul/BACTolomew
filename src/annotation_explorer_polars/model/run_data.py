"""RunData: the immutable data of one loaded run (window plan 4.3, Phase 2).

Everything here is derived from the run's tables once, when the run is
opened, and is then only read. ``RunData.build`` runs in the run loader's
worker thread (``app._RunDataWorker``), so neither this module nor what it
imports loads Qt (tests/test_run_data.py checks this).

Shared by every explore window of the run (sharing read-only data is not
linking). The pandas and polars frames are shared objects: treat them as read
only. The maps are ``MappingProxyType`` views, the lists are tuples, and the
plasmid cluster index is frozen (``PlasmidClusterIndex.freeze``). ``pan``
(``PanData``, Phase 9) describes the pangenome outputs for the pan panels and
loads them lazily, on first use.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import pandas as pd
import polars as pl

from pipeline.DbCatalog import MARKER_DB_CATALOG_PATH, load_db_catalog
from pipeline.UserDatabases import USER_LABEL_SUFFIX, run_user_databases, user_label

from ..constants import MARKER_SEARCH_TEXT_COLUMN, PLASMID_TRUE_VALUES
from ..data import build_marker_search_text_map, category_tags
from ..panaroo_data import _extract_mge_labels, _extract_virulence_labels
from ..plasmid_model import PlasmidClusterIndex, build_plasmid_cluster_index
from .pan_data import PanData


SAMPLE_QC_META_COLUMNS = frozenset(
    {
        "Quality",
        "QuastTotalLength",
        "QuastContigs",
        "QuastLargestContig",
        "QuastN50",
        "QuastL50",
        "QuastGC",
        "QuastNs",
        "QuastNsPer100kbp",
        "CheckM2Completeness",
        "CheckM2Contamination",
        "CheckM2CompletenessModel",
        "QualityNotes",
    }
)

_MARKER_HITS_COLUMNS = (
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
    "Product",
    "Identity",
    "Bitscore",
    "Subject_coverage",
    "SubjectIdentifier",
    "SubjectProteinAccession",
    "SubjectDNAAccession",
    "SubjectAROAccession",
    "SubjectModelName",
    "ResistanceMetadataSource",
    "ResistanceAlleles",
    "ResistanceGeneFamilies",
    "ResistanceDrugClasses",
    "ResistanceAntibiotics",
    "ResistanceMechanisms",
)

_MARKER_HITS_TEXT_COLUMNS = (
    "AnnotationID",
    "Sample",
    "Contig",
    "DB",
    "Method",
    "Product",
    "SubjectIdentifier",
    "SubjectProteinAccession",
    "SubjectDNAAccession",
    "SubjectAROAccession",
    "SubjectModelName",
    "ResistanceMetadataSource",
    "ResistanceAlleles",
    "ResistanceGeneFamilies",
    "ResistanceDrugClasses",
    "ResistanceAntibiotics",
    "ResistanceMechanisms",
)

_MOBSUITE_TEXT_COLUMNS = (
    "MobSuitePlasmidID",
    "MobSuitePrimaryClusterID",
    "MobSuiteSecondaryClusterID",
    "MobSuiteMoleculeType",
    "MobSuiteRepliconTypes",
    "MobSuiteRelaxaseTypes",
    "MobSuiteMPFType",
    "MobSuiteOriTTypes",
    "MobSuitePredictedMobility",
    "MobSuiteMashNearestNeighbor",
    "MobSuiteMashNeighborDistance",
)


# Text helpers

def clean_text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and pd.isna(value):
        return ""
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    return text


def clean_series(s: pd.Series) -> pd.Series:
    """Vectorised equivalent of applying clean_text over a whole column."""
    cleaned = s.where(s.notna(), "").astype(str).str.strip()
    return cleaned.mask(cleaned.str.lower() == "nan", "")


def norm_text_expr(column: str) -> pl.Expr:
    return pl.col(column).cast(pl.Utf8, strict=False).fill_null("").str.strip_chars()


def _nonempty_text(value: Any) -> str:
    text = str(value or "").strip()
    return "" if text.lower() in {"", "nan", "none", "null", "<na>"} else text


def mobsuite_plasmid_group_key(row: Mapping[str, Any]) -> str:
    """Key that groups the contigs of one MOB-suite plasmid within a sample."""
    if str(row.get("MobSuiteIsPlasmid") or "").strip().lower() not in PLASMID_TRUE_VALUES:
        return ""
    plasmid_id = _nonempty_text(row.get("MobSuitePlasmidID"))
    if plasmid_id:
        return plasmid_id
    sample = _nonempty_text(row.get("Sample"))
    primary = _nonempty_text(row.get("MobSuitePrimaryClusterID"))
    if sample and primary:
        return f"{sample}:{primary}"
    secondary = _nonempty_text(row.get("MobSuiteSecondaryClusterID"))
    if sample and secondary:
        return f"{sample}:secondary:{secondary}"
    return ""


# Database classification (used by the sample tree columns)

def db_is_resistance(db_name: str, category_map: Mapping[str, str]) -> bool:
    name = clean_text(db_name)
    if not name:
        return False
    category_upper = clean_text(category_map.get(name, "")).upper()
    if "AMR" in category_upper or "BMRG" in category_upper:
        return True
    name_upper = name.upper()
    return "AMR" in name_upper or "BMRG" in name_upper or name_upper.startswith("CARD_")


def db_is_mge_source(db_name: str, category_map: Mapping[str, str]) -> bool:
    name = clean_text(db_name)
    if not name:
        return False
    name_upper = name.upper()
    if "ISFINDER" in name_upper or "CGE-MGE" in name_upper or name_upper == "BAKTA IS":
        return True
    category_upper = clean_text(category_map.get(name, "")).upper()
    return category_upper == "MOBILE" and ("ISFINDER" in name_upper or "MGE" in name_upper)


def db_is_virulence_source(db_name: str, category_map: Mapping[str, str]) -> bool:
    name = clean_text(db_name)
    if not name:
        return False
    name_upper = name.upper()
    if "VIRULENCE" in name_upper or "VFDB" in name_upper:
        return True
    return clean_text(category_map.get(name, "")).upper() == "VIRULENCE"


# Build steps (each was a step of the explore window's __init__ before Phase 2)

def normalize_marker_hits_table(marker_hits_df: pd.DataFrame) -> pd.DataFrame:
    if marker_hits_df.empty:
        return pd.DataFrame(columns=list(_MARKER_HITS_COLUMNS))
    required = {"Sample", "Contig", "DB"}
    if not required.issubset(set(marker_hits_df.columns)):
        return pd.DataFrame(columns=list(marker_hits_df.columns))

    out = marker_hits_df.copy()
    for col in _MARKER_HITS_TEXT_COLUMNS:
        if col not in out.columns:
            out[col] = pd.NA
    for col in ["AnnotationID", "Sample", "Contig", "DB"]:
        out[col] = clean_series(out[col])
    out = out[(out["Sample"] != "") & (out["Contig"] != "") & (out["DB"] != "")]
    return out.reset_index(drop=True)


def attach_marker_search_text(ann_df: pd.DataFrame, marker_hits_df: pd.DataFrame) -> pd.DataFrame:
    out = ann_df.copy()
    if MARKER_SEARCH_TEXT_COLUMN not in out.columns:
        out[MARKER_SEARCH_TEXT_COLUMN] = pd.NA
    if out.empty or marker_hits_df.empty or "AnnotationID" not in out.columns:
        return out
    marker_text_map = build_marker_search_text_map(marker_hits_df)
    if not marker_text_map:
        return out
    mapped = out["AnnotationID"].fillna("").astype(str).str.strip().map(marker_text_map)
    mapped = mapped.replace({"": pd.NA})
    out[MARKER_SEARCH_TEXT_COLUMN] = mapped.fillna(out[MARKER_SEARCH_TEXT_COLUMN])
    return out


def prepare_metadata_polars(df: pl.DataFrame) -> pl.DataFrame:
    search_parts = [norm_text_expr(str(col)).str.to_lowercase() for col in df.columns]
    return (
        df.lazy()
        .with_columns(
            [
                norm_text_expr("Sample").alias("__sample_norm"),
                (
                    pl.concat_str(search_parts, separator=" ")
                    if search_parts
                    else pl.lit("", dtype=pl.Utf8)
                ).alias("__search_text"),
            ]
        )
        .collect()
    )


def prepare_annotation_polars(df: pl.DataFrame) -> pl.DataFrame:
    def _expr_or_empty(column: str) -> pl.Expr:
        if column in df.columns:
            return norm_text_expr(column)
        return pl.lit("", dtype=pl.Utf8)

    return (
        df.lazy()
        .with_columns(
            [
                _expr_or_empty("Sample").alias("__sample_norm"),
                _expr_or_empty("Contig").alias("__contig_norm"),
                _expr_or_empty("Gene").str.to_lowercase().alias("__gene_norm"),
                _expr_or_empty("Product").str.to_lowercase().alias("__product_norm"),
                _expr_or_empty("Predictors").str.to_lowercase().alias("__predictors_norm"),
                _expr_or_empty("Method").str.to_lowercase().alias("__method_norm"),
                _expr_or_empty("MarkerDisplay").str.to_lowercase().alias("__marker_display_norm"),
                _expr_or_empty(MARKER_SEARCH_TEXT_COLUMN).str.to_lowercase().alias("__marker_search_norm"),
                _expr_or_empty("Hypothetical").str.to_lowercase().alias("__hypothetical_norm"),
                _expr_or_empty("MobSuiteIsPlasmid").str.to_lowercase().alias("__mobsuite_plasmid_norm"),
                _expr_or_empty("Tags").str.to_uppercase().alias("__tags_upper"),
            ]
        )
        .with_columns(
            [
                pl.concat_str(
                    [
                        pl.col("__sample_norm").str.to_lowercase(),
                        pl.col("__contig_norm").str.to_lowercase(),
                        pl.col("__gene_norm"),
                        pl.col("__product_norm"),
                        pl.col("__predictors_norm"),
                        pl.col("__method_norm"),
                        pl.col("__marker_display_norm"),
                        pl.col("__marker_search_norm"),
                        pl.col("__tags_upper").str.to_lowercase(),
                    ],
                    separator=" ",
                ).alias("__search_text"),
                pl.col("__hypothetical_norm").is_in(["true", "1", "yes"]).alias("__is_hypothetical"),
                pl.col("__mobsuite_plasmid_norm").is_in(list(PLASMID_TRUE_VALUES)).alias("__is_mobsuite_plasmid"),
            ]
        )
        .collect()
    )


def build_contig_totals(ann_pl: pl.DataFrame) -> pl.DataFrame:
    mobsuite_is_plasmid_expr = (
        pl.col("MobSuiteIsPlasmid")
        .cast(pl.Utf8, strict=False)
        .str.strip_chars()
        .str.to_lowercase()
        .is_in(["true", "1", "yes", "y", "p", "plasmid"])
        .fill_null(False)
        .any()
        .alias("MobSuiteIsPlasmid")
        if "MobSuiteIsPlasmid" in ann_pl.columns
        else pl.lit(False).alias("MobSuiteIsPlasmid")
    )

    def _first_clean_text_expr(column: str) -> pl.Expr:
        if column not in ann_pl.columns:
            return pl.lit("").alias(column)
        text = pl.col(column).cast(pl.Utf8, strict=False).str.strip_chars()
        return (
            pl.when(text.str.to_lowercase().is_in(["", "nan", "none", "null", "<na>"]))
            .then(None)
            .otherwise(text)
            .drop_nulls()
            .first()
            .fill_null("")
            .alias(column)
        )

    # maintain_order: contigs (and so the member order of a multi-contig
    # plasmid) in annotation order, the same in every run.
    return (
        ann_pl.lazy()
        .group_by(["Sample", "Contig"], maintain_order=True)
        .agg(
            [
                pl.len().alias("n_total"),
                mobsuite_is_plasmid_expr,
                *[_first_clean_text_expr(column) for column in _MOBSUITE_TEXT_COLUMNS],
            ]
        )
        .collect()
    )


def build_mobsuite_plasmid_members(contigs_pl: pl.DataFrame) -> Dict[Tuple[str, str], Tuple[str, ...]]:
    members: Dict[Tuple[str, str], List[str]] = {}
    if contigs_pl.is_empty():
        return {}
    for row in contigs_pl.iter_rows(named=True):
        sample = _nonempty_text(row.get("Sample"))
        contig = _nonempty_text(row.get("Contig"))
        group_key = mobsuite_plasmid_group_key(row)
        if not sample or not contig or not group_key:
            continue
        bucket = members.setdefault((sample, group_key), [])
        if contig not in bucket:
            bucket.append(contig)
    return {key: tuple(values) for key, values in members.items()}


def _read_run_inputs(output_dir: Optional[Path]) -> Dict[str, Any]:
    if output_dir is None:
        return {}
    run_cfg_path = output_dir / "run_config.json"
    if not run_cfg_path.exists():
        return {}
    try:
        data = json.loads(run_cfg_path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    if not isinstance(data, dict):
        return {}
    inputs = data.get("inputs")
    return inputs if isinstance(inputs, dict) else {}


def load_run_user_databases(output_dir: Optional[Path]) -> frozenset[str]:
    """The run's user databases (pipeline.UserDatabases.run_user_databases): their hits
    are shown as "<name> (U)"; none for a run made before the bakta_evidence step."""
    return run_user_databases(_read_run_inputs(output_dir))


def _join_labels(value: Any, labels: Mapping[str, str]) -> Any:
    """A comma separated list of database names with each user database labelled."""
    text = clean_text(value)
    if not text:
        return value
    return ",".join(labels.get(part.strip(), part.strip()) for part in text.split(","))


def label_user_databases(
    ann_df: pd.DataFrame,
    marker_hits_df: pd.DataFrame,
    ref_df: pd.DataFrame,
    user_dbs: frozenset[str],
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """The tables with every user database named "<name> (U)" where its hits are shown
    (DB, Method, the table's Marker column, Predictors, the reference loci's Method). The
    run's files and the table cache keep the plain names."""
    if not user_dbs:
        return ann_df, marker_hits_df, ref_df
    labels = {name: user_label(name) for name in user_dbs}
    hits = marker_hits_df
    if not hits.empty and "DB" in hits.columns:
        hits = hits.copy()
        user_rows = hits["DB"].isin(user_dbs)
        if user_rows.any():
            hits.loc[user_rows, "DB"] = hits.loc[user_rows, "DB"].map(labels)
            if "Method" in hits.columns:
                hits.loc[user_rows, "Method"] = hits.loc[user_rows, "DB"]
            if "MarkerDisplay" in ann_df.columns and "Method" in hits.columns:
                ids = set(hits.loc[user_rows, "AnnotationID"].astype(str))
                methods = (
                    hits[hits["AnnotationID"].astype(str).isin(ids)]
                    .groupby("AnnotationID", sort=False)["Method"]
                    .agg(lambda values: ",".join(dict.fromkeys(clean_text(v) for v in values if clean_text(v))))
                )
                ann_df = ann_df.copy()
                touched = ann_df["AnnotationID"].astype(str).isin(ids)
                ann_df.loc[touched, "MarkerDisplay"] = ann_df.loc[touched, "AnnotationID"].astype(str).map(methods)
    if "Predictors" in ann_df.columns:
        ann_df = ann_df.copy()
        has = ann_df["Predictors"].notna()
        ann_df.loc[has, "Predictors"] = ann_df.loc[has, "Predictors"].map(lambda v: _join_labels(v, labels))
    if not ref_df.empty:
        ref_df = ref_df.copy()
        for col in ("DB", "Method"):
            if col in ref_df.columns:
                ref_df[col] = ref_df[col].map(lambda v: labels.get(clean_text(v), v))
    return ann_df, hits, ref_df


def load_run_database_names(output_dir: Optional[Path]) -> Tuple[Tuple[str, ...], Tuple[str, ...]]:
    inputs = _read_run_inputs(output_dir)

    def _stems(key: str) -> Tuple[str, ...]:
        vals = inputs.get(key)
        if not isinstance(vals, list):
            return ()
        out: List[str] = []
        seen: set[str] = set()
        for raw in vals:
            stem = Path(str(raw)).stem.strip()
            if not stem or stem in seen:
                continue
            seen.add(stem)
            out.append(stem)
        return tuple(out)

    return _stems("marker_sequences"), _stems("reference_locus_sequences")


# The tag of the compound (reference locus) databases: their hits tag overlapping genes.
COMPOUND_TAG = "I"


# The evidence source that is a keyword rule on Bakta's product, not a database search.
PRODUCT_NAME_SOURCE = "Bakta product name"
# Every annotation tag a run with gene evidence accounts for.
EVIDENCE_RUN_TAGS = ("A", "B", "V", "M", "P", COMPOUND_TAG, "X")


@dataclass(frozen=True)
class TagDatabases:
    """The databases behind one annotation tag in a run: the marker databases it searched,
    the catalog's ones (or, with gene evidence, the sources) that were not installed, and
    the gene evidence sources it had (Bakta VFDB, AMRFinderPlus, ...; ``evidence`` is True
    for a run made with the bakta_evidence step)."""

    used: Tuple[str, ...] = ()
    not_installed: Tuple[str, ...] = ()
    sources: Tuple[str, ...] = ()
    evidence: bool = False

    def _keyword_only(self) -> bool:
        return not self.used and self.sources == (PRODUCT_NAME_SOURCE,)

    def note(self) -> str:
        """"not installed" when the tag had no source at all, "Bakta product names" when
        its only source is the product name match, "partly installed" when some of its
        databases or sources were missing, else ""."""
        if self.evidence:
            if not self.used and not self.sources:
                return "not installed"
            if self._keyword_only():
                return "Bakta product names"
            return "partly installed" if self.not_installed else ""
        if not self.not_installed:
            return ""
        return "partly installed" if self.used else "not installed"

    def tooltip(self) -> str:
        if self.evidence:
            parts: List[str] = []
            if self._keyword_only():
                parts.append("From Bakta product names (keywords), not from a database search.")
            elif self.sources or self.used:
                parts.append("Sources: " + ", ".join((*self.sources, *self.used)) + ".")
            if self.not_installed:
                parts.append("Not installed for this run: " + ", ".join(self.not_installed) + ".")
            if not self.used and not self.sources:
                parts.append("No source for it in this run.")
            return " ".join(parts)
        if not self.not_installed:
            return ""
        text = "Not installed for this run: " + ", ".join(self.not_installed) + "."
        if self.used:
            return text + " Searched: " + ", ".join(self.used) + "."
        return text + " No database was searched for it."


def not_installed_key_note(
    tag_databases: Mapping[str, TagDatabases], tags: Sequence[Tuple[str, str]]
) -> Tuple[str, str]:
    """(short text, tooltip) for a view's key: the text names the not installed databases
    of these (tag, label)s ("" when every one was installed), the tooltip each tag's
    sources and missing databases (a run with gene evidence always names its sources;
    an older run only the missing ones, "" when none was)."""
    names: List[str] = []
    lines: List[str] = []
    for tag, label in tags:
        databases = tag_databases.get(tag)
        if databases is None:
            continue
        missing = list(databases.not_installed)
        if not missing and databases.note() == "not installed":
            missing = [label]
        for name in missing:
            # A user database keeps its whole "<name> (U)".
            short = name if name.endswith(USER_LABEL_SUFFIX) else (re.split(r"[-_]", name, maxsplit=1)[0] or name)
            if short not in names:
                names.append(short)
        tooltip = databases.tooltip()
        if tooltip and (missing or databases.evidence):
            lines.append(f"{label}: {tooltip}")
    text = "not installed: " + ", ".join(names) if names else ""
    return text, "\n".join(lines)


def load_run_tag_databases(output_dir: Optional[Path]) -> Dict[str, TagDatabases]:
    """Tag letter -> its databases, from the run's recorded `inputs.databases` (with the
    gene evidence sources of the bakta_evidence step); empty for a run made before the
    databases were optional (those had every one installed)."""
    inputs = _read_run_inputs(output_dir)
    record = inputs.get("databases")
    if not isinstance(record, dict):
        return {}
    user_dbs = run_user_databases(inputs)
    used: Dict[str, List[str]] = {}
    missing: Dict[str, List[str]] = {}
    for kind in ("marker", "compound"):
        entries = record.get(kind) if isinstance(record.get(kind), dict) else {}
        for key, bucket in (("used", used), ("not_installed", missing)):
            for db in entries.get(key) or []:
                if not isinstance(db, dict):
                    continue
                name = clean_text(db.get("name"))
                if name in user_dbs:
                    name = user_label(name)
                tags = COMPOUND_TAG if kind == "compound" else category_tags(clean_text(db.get("category")))
                for tag in tags:
                    if name and name not in bucket.setdefault(tag, []):
                        bucket[tag].append(name)
    evidence = record.get("evidence")
    if not isinstance(evidence, dict):
        return {
            tag: TagDatabases(tuple(used.get(tag, ())), tuple(missing.get(tag, ())))
            for tag in sorted(set(used) | set(missing))
        }

    def _names(key: str, tag: str) -> Tuple[str, ...]:
        values = (evidence.get(key) or {}).get(tag) if isinstance(evidence.get(key), dict) else None
        return tuple(clean_text(v) for v in values or [] if clean_text(v))

    return {
        tag: TagDatabases(
            tuple(used.get(tag, ())),
            (*_names("missing", tag), *missing.get(tag, ())),
            _names("tags", tag),
            evidence=True,
        )
        for tag in EVIDENCE_RUN_TAGS
    }


def load_snippy_reference_sample(output_dir: Optional[Path]) -> Optional[str]:
    ref_path = _read_run_inputs(output_dir).get("snippy_reference")
    if not ref_path:
        return None
    stem = Path(str(ref_path)).stem.strip()
    return stem or None


def load_marker_db_category_map() -> Dict[str, str]:
    try:
        raw = load_db_catalog(MARKER_DB_CATALOG_PATH)
    except Exception:
        return {}
    if not isinstance(raw, dict):
        return {}

    out: Dict[str, str] = {}
    for name, info in raw.items():
        if not isinstance(info, dict):
            continue
        category = clean_text(info.get("category"))
        if not category:
            continue
        key = clean_text(name)
        if key:
            out[key] = category
        db_path = clean_text(info.get("db_path") or info.get("path"))
        if db_path:
            stem = Path(db_path).stem.strip()
            if stem:
                out.setdefault(stem, category)
    return out


def map_reference_db_name(value: str, configured_reference_dbs: Sequence[str]) -> str:
    db = clean_text(value)
    if not db:
        return ""
    if not configured_reference_dbs:
        return db
    if db in configured_reference_dbs:
        return db
    candidates = [name for name in configured_reference_dbs if name.startswith(f"{db}_")]
    if len(candidates) == 1:
        return candidates[0]
    return db


def build_reference_hits_table(ref_df: pd.DataFrame, configured_reference_dbs: Sequence[str]) -> pd.DataFrame:
    cols = ["Sample", "Contig", "DB"]
    if ref_df.empty or "Sample" not in ref_df.columns or "Contig" not in ref_df.columns:
        return pd.DataFrame(columns=cols)

    db_col = "DB" if "DB" in ref_df.columns else "Method" if "Method" in ref_df.columns else None
    if db_col is None:
        return pd.DataFrame(columns=cols)

    rows: List[Tuple[str, str, str]] = []
    for sample, contig, db_raw in ref_df[["Sample", "Contig", db_col]].itertuples(index=False, name=None):
        sample_name = clean_text(sample)
        contig_name = clean_text(contig)
        if not sample_name or not contig_name:
            continue
        mapped_db = map_reference_db_name(clean_text(db_raw), configured_reference_dbs)
        if not mapped_db:
            continue
        rows.append((sample_name, contig_name, mapped_db))
    if not rows:
        return pd.DataFrame(columns=cols)
    return pd.DataFrame(rows, columns=cols)


def collect_annotation_db_columns(
    configured_marker_dbs: Sequence[str],
    configured_reference_dbs: Sequence[str],
    marker_hits_df: pd.DataFrame,
    reference_hits_df: pd.DataFrame,
) -> Tuple[str, ...]:
    ordered: List[str] = []
    seen: set[str] = set()

    for db in [*configured_marker_dbs, *configured_reference_dbs]:
        name = clean_text(db)
        if not name or name in seen:
            continue
        seen.add(name)
        ordered.append(name)

    observed: List[str] = []
    if not marker_hits_df.empty and "DB" in marker_hits_df.columns:
        observed.extend(marker_hits_df["DB"].astype(str).tolist())
    if not reference_hits_df.empty and "DB" in reference_hits_df.columns:
        observed.extend(reference_hits_df["DB"].astype(str).tolist())

    for db in sorted({clean_text(v) for v in observed if clean_text(v)}):
        if db in seen:
            continue
        seen.add(db)
        ordered.append(db)

    return tuple(ordered)


def format_resistance_drug_class_label(value: Any) -> str:
    label = clean_text(value).replace("_", " ")
    label = re.sub(r"\s+", " ", label).strip()
    if not label:
        return ""
    if label.upper() == label:
        return label.title()
    return label[0].upper() + label[1:]


def split_resistance_drug_classes(value: Any) -> List[str]:
    text = clean_text(value)
    if not text:
        return []
    out: List[str] = []
    seen: set[str] = set()
    for part in re.split(r"[;,]", text):
        label = format_resistance_drug_class_label(part)
        if not label:
            continue
        key = label.casefold()
        if key in seen:
            continue
        seen.add(key)
        out.append(label)
    return out


def build_resistance_atb_count_maps(
    marker_hits_df: pd.DataFrame,
) -> Tuple[Dict[str, int], Dict[Tuple[str, str], int]]:
    if marker_hits_df.empty or "ResistanceDrugClasses" not in marker_hits_df.columns:
        return {}, {}

    sample_classes: Dict[str, set[str]] = {}
    contig_classes: Dict[Tuple[str, str], set[str]] = {}
    cols = ["Sample", "Contig", "ResistanceDrugClasses"]
    for sample_raw, contig_raw, classes_raw in marker_hits_df[cols].itertuples(index=False, name=None):
        sample_name = clean_text(sample_raw)
        contig_name = clean_text(contig_raw)
        if not sample_name:
            continue
        classes = split_resistance_drug_classes(classes_raw)
        if not classes:
            continue
        sample_classes.setdefault(sample_name, set()).update(classes)
        if contig_name:
            contig_classes.setdefault((sample_name, contig_name), set()).update(classes)

    sample_counts = {sample: len(values) for sample, values in sample_classes.items()}
    contig_counts = {key: len(values) for key, values in contig_classes.items()}
    return sample_counts, contig_counts


def build_annotation_hit_counts(marker_hits_df: pd.DataFrame, reference_hits_df: pd.DataFrame) -> pl.DataFrame:
    empty = pl.DataFrame(schema={"Sample": pl.Utf8, "Contig": pl.Utf8, "DB": pl.Utf8, "Hits": pl.Int64})
    frames: List[pl.DataFrame] = []
    if not marker_hits_df.empty:
        frames.append(pl.from_pandas(marker_hits_df[["Sample", "Contig", "DB"]]))
    if not reference_hits_df.empty:
        frames.append(pl.from_pandas(reference_hits_df[["Sample", "Contig", "DB"]]))
    if not frames:
        return empty

    combined = pl.concat(frames, how="vertical_relaxed")
    counts = (
        combined.lazy()
        .with_columns(
            [
                pl.col("Sample").cast(pl.Utf8).str.strip_chars().alias("Sample"),
                pl.col("Contig").cast(pl.Utf8).str.strip_chars().alias("Contig"),
                pl.col("DB").cast(pl.Utf8).str.strip_chars().alias("DB"),
            ]
        )
        .filter(
            (pl.col("Sample").is_not_null())
            & (pl.col("Sample") != "")
            & (pl.col("Contig").is_not_null())
            & (pl.col("Contig") != "")
            & (pl.col("DB").is_not_null())
            & (pl.col("DB") != "")
        )
        .group_by(["Sample", "Contig", "DB"])
        .agg(pl.len().alias("Hits"))
        .sort(["Sample", "Contig", "DB"], nulls_last=True)
        .collect()
    )
    if counts.is_empty():
        return empty
    return counts.select(["Sample", "Contig", "DB", "Hits"])


def collect_predictor_options(ann_pl: pl.DataFrame) -> Tuple[str, ...]:
    if "Predictors" not in ann_pl.columns:
        return ()
    values = (
        ann_pl.lazy()
        .select(pl.col("Predictors").cast(pl.Utf8).fill_null("").str.split(",").alias("__pred"))
        # Explicit: polars 2 changed the default to False (1.x warns).
        .explode("__pred", empty_as_null=True)
        .with_columns(pl.col("__pred").cast(pl.Utf8).str.strip_chars().alias("__pred"))
        .filter(pl.col("__pred").is_not_null() & (pl.col("__pred") != ""))
        .select(pl.col("__pred").unique().sort())
        .collect()
        .get_column("__pred")
        .to_list()
    )
    return tuple(str(v) for v in values)


def _distinct_dbs(marker_hits_df: pd.DataFrame) -> List[str]:
    return sorted({str(db).strip() for db in marker_hits_df["DB"].dropna().astype(str).tolist() if str(db).strip()})


def detect_resistance_traits(marker_hits_df: pd.DataFrame, category_map: Mapping[str, str]) -> bool:
    if marker_hits_df is None or marker_hits_df.empty:
        return False
    if "DB" in marker_hits_df.columns:
        if any(db_is_resistance(db, category_map) for db in _distinct_dbs(marker_hits_df)):
            return True
    for col in ("ResistanceDrugClasses", "ResistanceAlleles", "ResistanceGeneFamilies"):
        if col in marker_hits_df.columns and marker_hits_df[col].fillna("").astype(str).str.strip().ne("").any():
            return True
    return False


def _detect_labelled_traits(
    marker_hits_df: pd.DataFrame,
    is_source: Callable[[str], bool],
    extract_labels: Callable[[Dict[str, Any]], Any],
) -> bool:
    if marker_hits_df is None or marker_hits_df.empty:
        return False
    if "DB" in marker_hits_df.columns:
        source_dbs = [db for db in _distinct_dbs(marker_hits_df) if is_source(db)]
        if not source_dbs:
            return False
        db_values = marker_hits_df["DB"].fillna("").astype(str).str.strip()
        subset = marker_hits_df[db_values.isin(source_dbs)]
    else:
        subset = marker_hits_df
    if subset.empty:
        return False
    return any(extract_labels(row) for row in subset.to_dict(orient="records"))


def detect_mge_traits(marker_hits_df: pd.DataFrame, category_map: Mapping[str, str]) -> bool:
    return _detect_labelled_traits(
        marker_hits_df, lambda db: db_is_mge_source(db, category_map), _extract_mge_labels
    )


def detect_virulence_traits(marker_hits_df: pd.DataFrame, category_map: Mapping[str, str]) -> bool:
    return _detect_labelled_traits(
        marker_hits_df, lambda db: db_is_virulence_source(db, category_map), _extract_virulence_labels
    )


@dataclass(frozen=True, eq=False)
class RunData:
    # pandas
    ann_df: pd.DataFrame
    meta_df: pd.DataFrame
    ref_df: pd.DataFrame
    marker_hits_df: pd.DataFrame
    reference_hits_df: pd.DataFrame
    # polars (ann_pl and meta_pl carry the "__*_norm" and "__search_text" columns)
    ann_pl: pl.DataFrame
    meta_pl: pl.DataFrame
    contig_totals_pl: pl.DataFrame
    annotation_hit_counts_pl: pl.DataFrame
    # maps
    mobsuite_plasmid_members: Mapping[Tuple[str, str], Tuple[str, ...]]
    plasmid_cluster_index: PlasmidClusterIndex
    contig_total_count_map: Mapping[Tuple[str, str], int]
    sample_total_counts: Mapping[str, int]
    sample_atb_counts: Mapping[str, int]
    contig_atb_counts: Mapping[Tuple[str, str], int]
    # column lists
    ann_public_columns: Tuple[str, ...]
    meta_public_columns: Tuple[str, ...]
    meta_columns: Tuple[str, ...]
    sample_qc_columns: Tuple[str, ...]
    annotation_db_columns: Tuple[str, ...]
    predictor_options: Tuple[str, ...]
    # run config
    configured_marker_dbs: Tuple[str, ...]
    configured_reference_dbs: Tuple[str, ...]
    snippy_reference_sample: Optional[str]
    marker_db_category_map: Mapping[str, str]
    # Annotation tag -> the databases behind it (shown as not installed when missing).
    tag_databases: Mapping[str, TagDatabases]
    # Trait sources present in the marker hits
    has_resistance_traits: bool
    has_mge_traits: bool
    has_virulence_traits: bool
    output_dir: Optional[Path]
    run_id: Optional[str]
    # The pangenome outputs, shared by every pan panel (plan 8.10, Phase 9).
    pan: PanData

    @classmethod
    def build(
        cls,
        ann_df: pd.DataFrame,
        meta_df: pd.DataFrame,
        *,
        ref_df: Optional[pd.DataFrame] = None,
        marker_hits_df: Optional[pd.DataFrame] = None,
        output_dir: Optional[Path] = None,
        run_id: Optional[str] = None,
        status: Optional[Callable[[str], None]] = None,
    ) -> "RunData":
        """Derive everything the explore windows read from the run's tables.

        Runs off the GUI thread in the app; ``status`` reports progress for
        the loading splash.
        """
        if marker_hits_df is None:
            raise ValueError("marker_hits_df is required. Rebuild cache using the new schema.")

        def _status(text: str) -> None:
            if status is not None:
                status(text)

        output_path = Path(output_dir) if output_dir is not None else None
        _status("Preparing annotations…")
        user_dbs = load_run_user_databases(output_path)
        ann_df, marker_hits, ref = label_user_databases(
            ann_df,
            normalize_marker_hits_table(marker_hits_df),
            ref_df.copy() if ref_df is not None else pd.DataFrame(),
            user_dbs,
        )
        ann_df = attach_marker_search_text(ann_df, marker_hits)
        ann_pl = prepare_annotation_polars(pl.from_pandas(ann_df))
        meta_pl = prepare_metadata_polars(pl.from_pandas(meta_df))

        _status("Indexing contigs and plasmids…")
        contig_totals_pl = build_contig_totals(ann_pl)
        mobsuite_members = build_mobsuite_plasmid_members(contig_totals_pl)
        plasmid_index = build_plasmid_cluster_index(output_dir=output_path, ann_df=ann_df).freeze()
        contig_total_count_map: Dict[Tuple[str, str], int] = {
            (str(row.get("Sample") or "").strip(), str(row.get("Contig") or "").strip()): int(row.get("n_total") or 0)
            for row in contig_totals_pl.iter_rows(named=True)
            if str(row.get("Sample") or "").strip() and str(row.get("Contig") or "").strip()
        }
        sample_total_counts: Dict[str, int] = {
            str(row.get("Sample") or "").strip(): int(row.get("n_total") or 0)
            for row in (
                ann_pl.lazy()
                .group_by("Sample")
                .agg(pl.len().alias("n_total"))
                .collect()
                .iter_rows(named=True)
            )
            if str(row.get("Sample") or "").strip()
        }

        _status("Counting database hits…")
        configured_marker_dbs, configured_reference_dbs = (
            tuple(user_label(name) if name in user_dbs else name for name in names)
            for names in load_run_database_names(output_path)
        )
        # The gene evidence sources' categories, then the marker catalog's.
        category_map = {"Bakta VFDB": "VIRULENCE", "Bakta IS": "MOBILE", **load_marker_db_category_map()}
        recorded = _read_run_inputs(output_path).get("marker_categories") or {}
        # A database this install's catalog lacks (a run made elsewhere, no DB/)
        # keeps the category the run recorded.
        for name, category in recorded.items() if isinstance(recorded, dict) else ():
            if clean_text(name) and clean_text(category):
                category_map.setdefault(clean_text(name), clean_text(category))
        for name in user_dbs:
            category = category_map.get(name) or (recorded.get(name) if isinstance(recorded, dict) else None)
            if category:
                category_map[user_label(name)] = clean_text(category)
        reference_hits_df = build_reference_hits_table(ref, configured_reference_dbs)
        sample_atb_counts, contig_atb_counts = build_resistance_atb_count_maps(marker_hits)
        annotation_db_columns = collect_annotation_db_columns(
            configured_marker_dbs, configured_reference_dbs, marker_hits, reference_hits_df
        )
        annotation_hit_counts_pl = build_annotation_hit_counts(marker_hits, reference_hits_df)
        has_resistance = detect_resistance_traits(marker_hits, category_map)
        has_mge = detect_mge_traits(marker_hits, category_map)
        has_virulence = detect_virulence_traits(marker_hits, category_map)
        snippy_reference_sample = load_snippy_reference_sample(output_path)

        return cls(
            ann_df=ann_df,
            meta_df=meta_df,
            ref_df=ref,
            marker_hits_df=marker_hits,
            reference_hits_df=reference_hits_df,
            ann_pl=ann_pl,
            meta_pl=meta_pl,
            contig_totals_pl=contig_totals_pl,
            annotation_hit_counts_pl=annotation_hit_counts_pl,
            mobsuite_plasmid_members=MappingProxyType(mobsuite_members),
            plasmid_cluster_index=plasmid_index,
            contig_total_count_map=MappingProxyType(contig_total_count_map),
            sample_total_counts=MappingProxyType(sample_total_counts),
            sample_atb_counts=MappingProxyType(sample_atb_counts),
            contig_atb_counts=MappingProxyType(contig_atb_counts),
            ann_public_columns=tuple(ann_df.columns),
            meta_public_columns=tuple(meta_df.columns),
            meta_columns=tuple(
                c for c in meta_df.columns if c != "Sample" and c not in SAMPLE_QC_META_COLUMNS
            ),
            sample_qc_columns=tuple(c for c in meta_df.columns if c in SAMPLE_QC_META_COLUMNS),
            annotation_db_columns=annotation_db_columns,
            predictor_options=collect_predictor_options(ann_pl),
            configured_marker_dbs=configured_marker_dbs,
            configured_reference_dbs=configured_reference_dbs,
            snippy_reference_sample=snippy_reference_sample,
            marker_db_category_map=MappingProxyType(category_map),
            tag_databases=MappingProxyType(load_run_tag_databases(output_path)),
            has_resistance_traits=has_resistance,
            has_mge_traits=has_mge,
            has_virulence_traits=has_virulence,
            output_dir=output_path,
            run_id=run_id,
            # Folder lookups and tables load lazily, on a pan panel's first use.
            pan=PanData(output_path, meta_df, plasmid_index, snippy_reference_sample),
        )
