from __future__ import annotations

from collections import OrderedDict
import bisect
import json
import re
import threading
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import polars as pl
from PyQt6 import QtCore, QtGui, QtWidgets

from pipeline.DbCatalog import MARKER_DB_CATALOG_PATH, load_db_catalog

from .annotation_detail_dialog import AnnotationDetailDialog
from .annotation_table import AnnotationTableModel, DetailsDelegate, TagsDelegate
from .constants import ANNOTATION_COLUMNS_TO_SHOW, MARKER_SEARCH_TEXT_COLUMN
from .data import build_marker_search_text_map
from .geo_experimental_views import GeoMapLibreView
from .geo_view import GeoMapView
from .panaroo_views import (
    EntropyTrackTab,
    PanCoreAccessoryTab,
    StructPresenceHeatmapTab,
)
from .pan_view import PanHeatmapWidget
from .cgmlst_force_plot import CgmlstForcePlotPanel
from .igv_view import GeneSegment, IgvTrackView
from .metadata_filter_dialog import MetadataValueDialog
from .explorer_ui import ExplorerUIMixin
from .filters_controller import FiltersControllerMixin
from .tree_controller import TreeControllerMixin
from .table_controller import TableControllerMixin
from .igv_controller import IgvControllerMixin
from .panels_controller import PanelsControllerMixin
from .hover_sync import HoverSyncMixin
from .plasmid_model import PlasmidClusterIndex, build_plasmid_cluster_index
from .polars_queries import PolarsQueriesMixin
from .refresh_worker import RefreshPayload, RefreshSnapshot, RefreshWorker


_SAMPLE_QC_META_COLUMNS = {
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


class AnnotationExplorer(
    ExplorerUIMixin,
    FiltersControllerMixin,
    TreeControllerMixin,
    TableControllerMixin,
    IgvControllerMixin,
    PanelsControllerMixin,
    HoverSyncMixin,
    PolarsQueriesMixin,
    QtWidgets.QMainWindow,
):
    # Emitted once, when the first refresh has populated the sample tree — used
    # to dismiss the loading splash exactly when the view is actually ready.
    firstRefreshFinished = QtCore.pyqtSignal()

    def __init__(
        self,
        ann_df: pd.DataFrame,
        meta_df: pd.DataFrame,
        ref_df: Optional[pd.DataFrame] = None,
        marker_hits_df: Optional[pd.DataFrame] = None,
        output_dir: Optional[Path] = None,
    ) -> None:
        super().__init__()
        if marker_hits_df is None:
            raise ValueError("marker_hits_df is required. Rebuild cache using the new schema.")
        self.marker_hits_df = self._normalize_marker_hits_table(marker_hits_df)
        self._yield_to_ui()
        ann_df = self._attach_marker_search_text(ann_df, self.marker_hits_df)
        self._yield_to_ui()
        self.ann_df = ann_df
        self.meta_df = meta_df
        self.ref_df = ref_df.copy() if ref_df is not None else pd.DataFrame()
        self.output_dir = Path(output_dir) if output_dir is not None else None
        self._sample_qc_columns = [c for c in meta_df.columns if c in _SAMPLE_QC_META_COLUMNS]
        self.meta_columns = [c for c in meta_df.columns if c != "Sample" and c not in _SAMPLE_QC_META_COLUMNS]
        # show all metadata columns in the tree
        self.meta_columns_for_tree = list(self.meta_columns)
        self.metadata_filters: Dict[str, List[str]] = {}
        self._sample_subset_filter_samples: List[str] = []
        self._sample_subset_filter_label: Optional[str] = None
        self._sample_subset_filter_kind: Optional[str] = None
        self._ann_public_columns = list(ann_df.columns)
        self._meta_public_columns = list(meta_df.columns)
        self._polars_cache_lock = threading.Lock()
        self._metadata_subset_cache: "OrderedDict[object, object]" = OrderedDict()
        self._annotation_filter_cache: "OrderedDict[object, object]" = OrderedDict()
        self._table_column_selected: List[str] = []
        self._table_column_group_checkboxes: Dict[str, QtWidgets.QCheckBox] = {}
        self._table_column_checkboxes: Dict[str, QtWidgets.QCheckBox] = {}
        self._table_column_group_members: Dict[str, List[str]] = {}
        self._table_column_ui_updating = False
        self._neo_tree_has_resistance_traits = False
        self._neo_tree_has_mge_traits = False
        self._neo_tree_has_virulence_traits = False
        self._neo_tree_trait_available: List[str] = []
        self._neo_tree_trait_selected: List[str] = []
        self._neo_tree_trait_checkboxes: Dict[str, QtWidgets.QCheckBox] = {}
        self._neo_tree_trait_ui_updating = False
        self._neo_tree_trait_select_all_checkbox: Optional[QtWidgets.QCheckBox] = None
        self._neo_tree_trait_color_signals_connected = False
        self._neo_tree_amr_color = QtGui.QColor("#f97316")
        self._neo_tree_mge_color = QtGui.QColor("#0f766e")
        self._neo_tree_virulence_color = QtGui.QColor("#be185d")
        self._ann_pl = self._prepare_annotation_polars(pl.from_pandas(ann_df))
        self._meta_pl = self._prepare_metadata_polars(pl.from_pandas(meta_df))
        mobsuite_is_plasmid_expr = (
            pl.col("MobSuiteIsPlasmid")
            .cast(pl.Utf8, strict=False)
            .str.strip_chars()
            .str.to_lowercase()
            .is_in(["true", "1", "yes", "y", "p", "plasmid"])
            .fill_null(False)
            .any()
            .alias("MobSuiteIsPlasmid")
            if "MobSuiteIsPlasmid" in self._ann_pl.columns
            else pl.lit(False).alias("MobSuiteIsPlasmid")
        )
        def _first_clean_text_expr(column: str) -> pl.Expr:
            if column not in self._ann_pl.columns:
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

        self._contig_totals_pl = (
            self._ann_pl.lazy()
            .group_by(["Sample", "Contig"])
            .agg(
                [
                    pl.len().alias("n_total"),
                    mobsuite_is_plasmid_expr,
                    _first_clean_text_expr("MobSuitePlasmidID"),
                    _first_clean_text_expr("MobSuitePrimaryClusterID"),
                    _first_clean_text_expr("MobSuiteSecondaryClusterID"),
                    _first_clean_text_expr("MobSuiteMoleculeType"),
                    _first_clean_text_expr("MobSuiteRepliconTypes"),
                    _first_clean_text_expr("MobSuiteRelaxaseTypes"),
                    _first_clean_text_expr("MobSuiteMPFType"),
                    _first_clean_text_expr("MobSuiteOriTTypes"),
                    _first_clean_text_expr("MobSuitePredictedMobility"),
                    _first_clean_text_expr("MobSuiteMashNearestNeighbor"),
                    _first_clean_text_expr("MobSuiteMashNeighborDistance"),
                ]
            )
            .collect()
        )
        self._yield_to_ui()
        self._mobsuite_plasmid_members = self._build_mobsuite_plasmid_members(self._contig_totals_pl)
        self._plasmid_cluster_index: PlasmidClusterIndex = build_plasmid_cluster_index(
            output_dir=self.output_dir,
            ann_df=ann_df,
        )
        self._yield_to_ui()
        self._contig_total_count_map: Dict[Tuple[str, str], int] = {
            (str(row.get("Sample") or "").strip(), str(row.get("Contig") or "").strip()): int(row.get("n_total") or 0)
            for row in self._contig_totals_pl.iter_rows(named=True)
            if str(row.get("Sample") or "").strip() and str(row.get("Contig") or "").strip()
        }
        self.sample_total_counts: Dict[str, int] = {
            str(row.get("Sample") or "").strip(): int(row.get("n_total") or 0)
            for row in (
                self._ann_pl.lazy()
                .group_by("Sample")
                .agg(pl.len().alias("n_total"))
                .collect()
                .iter_rows(named=True)
            )
            if str(row.get("Sample") or "").strip()
        }
        self._yield_to_ui()
        self._base_allowed_samples: List[str] = []
        self._base_ann_pl = self._empty_ann_public_pl()
        self._filtered_ann_pl = self._empty_ann_public_pl()
        self.base_ann = ann_df.iloc[0:0].copy()
        self.filtered_ann = ann_df.iloc[0:0].copy()

        self._configured_marker_dbs, self._configured_reference_dbs = self._load_run_database_names()
        self._marker_db_category_map = self._load_marker_db_category_map()
        self._snippy_reference_sample = self._load_snippy_reference_sample()
        self._yield_to_ui()
        self._neo_tree_has_resistance_traits = self._detect_neo_tree_resistance_traits()
        self._neo_tree_has_mge_traits = self._detect_neo_tree_mge_traits()
        self._neo_tree_has_virulence_traits = self._detect_neo_tree_virulence_traits()
        self._yield_to_ui()
        self._neo_tree_trait_available = self._collect_neo_tree_trait_available_columns()
        self._neo_tree_trait_selected = self._neo_tree_default_selected_traits()
        self.reference_hits_df = self._build_reference_hits_table(self.ref_df)
        self._yield_to_ui()
        self.sample_atb_counts, self.contig_atb_counts = self._build_resistance_atb_count_maps()
        self.annotation_db_columns = self._collect_annotation_db_columns()
        self._yield_to_ui()
        self.annotation_hit_counts_pl = self._build_annotation_hit_counts_polars()
        self._yield_to_ui()

        self.predictor_options = self._collect_predictor_options()
        self._yield_to_ui()
        self._tag_model: Optional[QtGui.QStandardItemModel] = None
        self._predictor_model: Optional[QtGui.QStandardItemModel] = None
        self.view_mode = "Contig"
        self.context_window = 30
        self.current_sample: Optional[str] = None
        self.current_contig: Optional[str] = None
        self.current_contigs: Optional[List[str]] = None
        self.current_contig_label: Optional[str] = None
        self._table_model: Optional[AnnotationTableModel] = None
        self._tags_delegate: Optional[TagsDelegate] = None
        self._details_delegate: Optional[DetailsDelegate] = None
        self._delegate_reset_column_count = 0
        self.igv_view: Optional[IgvTrackView] = None
        self._igv_key: Optional[Tuple[Any, ...]] = None
        self._igv_dirty = False
        self._igv_bp_segments: List[Tuple[int, int, float, float]] = []
        self._igv_raw_contig_len = 0
        self._igv_display_contig_len = 0
        self._syncing_table_to_igv = False
        self._syncing_igv_to_table = False
        self._pending_igv_center_bp: Optional[int] = None
        self._table_row_bp_index: List[Tuple[int, int]] = []
        self._table_bp_values: List[int] = []
        self.map_view: Optional[GeoMapView] = None
        self.map_web_view: Optional[GeoMapLibreView] = None
        self.pan_neo_tab: Optional[PanCoreAccessoryTab] = None
        self.pan_neo_tree_tab: Optional[PanCoreAccessoryTab] = None
        self.pan_similarity_tab: Optional[QtWidgets.QWidget] = None
        self.pan_structvar_tab: Optional[StructPresenceHeatmapTab] = None
        self.pan_entropy_tab: Optional[EntropyTrackTab] = None
        self._plot_menu: Optional[QtWidgets.QMenu] = None
        self._phandango_plot_action: Optional[QtGui.QAction] = None
        self._cgmlst_force_plot_action: Optional[QtGui.QAction] = None
        self._phandango_plot_dialog: Optional[QtWidgets.QDialog] = None
        self._phandango_plot_splitter: Optional[QtWidgets.QSplitter] = None
        self._phandango_presence_toggle: Optional[QtWidgets.QCheckBox] = None
        self._phandango_pan_heatmap: Optional[PanHeatmapWidget] = None
        self._cgmlst_force_plot_dialog: Optional[QtWidgets.QDialog] = None
        self._cgmlst_force_plot_panel: Optional[CgmlstForcePlotPanel] = None
        self._top_splitter: Optional[QtWidgets.QSplitter] = None
        self._filters_panel: Optional[QtWidgets.QWidget] = None
        self._sample_tree_items: Dict[str, QtWidgets.QTreeWidgetItem] = {}
        self._tree_hover_sample: Optional[str] = None
        self._tree_bulk_toggle_in_progress = False
        self._hover_active_sample: Optional[str] = None
        self._hover_active_source: Optional[QtCore.QObject] = None
        self._table_scroll_timer = QtCore.QTimer(self)
        self._table_scroll_timer.setSingleShot(True)
        self._table_scroll_timer.setInterval(40)
        self._table_scroll_timer.timeout.connect(self._sync_igv_from_table_scroll)
        self._igv_scroll_timer = QtCore.QTimer(self)
        self._igv_scroll_timer.setSingleShot(True)
        self._igv_scroll_timer.setInterval(40)
        self._igv_scroll_timer.timeout.connect(self._sync_table_from_igv_center)
        self._sync_animation_ms = 200
        self._table_scroll_animation: Optional[QtCore.QPropertyAnimation] = None
        self._table_sync_token = 0
        self._igv_sync_token = 0
        self._floating_dialogs: List[QtWidgets.QDialog] = []
        self._refresh_generation = 0
        self._first_refresh_done = False
        self._refresh_running = False
        self._refresh_pending_snapshot: Optional[RefreshSnapshot] = None
        self._refresh_worker: Optional[RefreshWorker] = None
        self._refresh_pool = QtCore.QThreadPool.globalInstance()
        self.ann_filters: Dict[str, Any] = {}
        self._deferred_panel_allowed_samples: List[str] = []
        self._deferred_panel_ann_filters: Dict[str, Any] = {}
        self._deferred_geo_dirty = False
        self._deferred_pan_tab_dirty: set[str] = set()
        self._deferred_panel_timer = QtCore.QTimer(self)
        self._deferred_panel_timer.setSingleShot(True)
        self._deferred_panel_timer.setInterval(250)
        self._deferred_panel_timer.timeout.connect(self._flush_deferred_panel_refresh)

        self._build_ui()
        self._render_metadata_filters()
        self.refresh()

    # UI setup

    def _show_floating_dialog(self, dlg: QtWidgets.QDialog) -> None:
        dlg.setWindowModality(QtCore.Qt.WindowModality.NonModal)
        dlg.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self._floating_dialogs.append(dlg)

        def _cleanup(*_args: object, _dlg: QtWidgets.QDialog = dlg) -> None:
            if _dlg in self._floating_dialogs:
                self._floating_dialogs.remove(_dlg)

        dlg.destroyed.connect(_cleanup)
        dlg.show()
        dlg.raise_()
        dlg.activateWindow()
    def _load_run_database_names(self) -> Tuple[List[str], List[str]]:
        if self.output_dir is None:
            return [], []
        run_cfg_path = self.output_dir / "run_config.json"
        if not run_cfg_path.exists():
            return [], []
        try:
            data = json.loads(run_cfg_path.read_text(encoding="utf-8"))
        except Exception:
            return [], []
        if not isinstance(data, dict):
            return [], []
        inputs = data.get("inputs")
        if not isinstance(inputs, dict):
            return [], []

        def _stems(key: str) -> List[str]:
            vals = inputs.get(key)
            if not isinstance(vals, list):
                return []
            out: List[str] = []
            seen: set[str] = set()
            for raw in vals:
                stem = Path(str(raw)).stem.strip()
                if not stem or stem in seen:
                    continue
                seen.add(stem)
                out.append(stem)
            return out

        return _stems("marker_sequences"), _stems("reference_locus_sequences")

    def _load_snippy_reference_sample(self) -> Optional[str]:
        if self.output_dir is None:
            return None
        run_cfg_path = self.output_dir / "run_config.json"
        if not run_cfg_path.exists():
            return None
        try:
            data = json.loads(run_cfg_path.read_text(encoding="utf-8"))
        except Exception:
            return None
        if not isinstance(data, dict):
            return None
        inputs = data.get("inputs")
        if not isinstance(inputs, dict):
            return None
        ref_path = inputs.get("snippy_reference")
        if not ref_path:
            return None
        stem = Path(str(ref_path)).stem.strip()
        return stem or None

    @staticmethod
    def _clean_text(value: Any) -> str:
        if value is None:
            return ""
        if isinstance(value, float) and pd.isna(value):
            return ""
        text = str(value).strip()
        if not text or text.lower() == "nan":
            return ""
        return text

    @staticmethod
    def _clean_series(s: pd.Series) -> pd.Series:
        """Vectorised equivalent of applying _clean_text over a whole column."""
        cleaned = s.where(s.notna(), "").astype(str).str.strip()
        return cleaned.mask(cleaned.str.lower() == "nan", "")

    @staticmethod
    def _yield_to_ui() -> None:
        # Let the loading splash repaint during the heavy synchronous data
        # preparation in __init__ (the standard Qt splash + processEvents idiom).
        # Safe here: this all runs before any of this window's own widgets or
        # signals exist, so there is nothing to re-enter.
        app = QtWidgets.QApplication.instance()
        if app is not None:
            app.processEvents()

    @staticmethod
    def _attach_marker_search_text(
        ann_df: pd.DataFrame,
        marker_hits_df: pd.DataFrame,
    ) -> pd.DataFrame:
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

    def _build_sample_metadata_lookup(self, meta_df: pd.DataFrame) -> Dict[str, Dict[str, str]]:
        out: Dict[str, Dict[str, str]] = {}
        if meta_df.empty or "Sample" not in meta_df.columns:
            return out
        for _, row in meta_df.iterrows():
            sample_name = self._clean_text(row.get("Sample"))
            if not sample_name or sample_name in out:
                continue
            out[sample_name] = {
                str(col): self._clean_text(row.get(col))
                for col in meta_df.columns
                if str(col) != "Sample"
            }
        return out

    def _load_marker_db_category_map(self) -> Dict[str, str]:
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
            category = self._clean_text(info.get("category"))
            if not category:
                continue
            key = self._clean_text(name)
            if key:
                out[key] = category
            db_path = self._clean_text(info.get("db_path") or info.get("path"))
            if db_path:
                stem = Path(db_path).stem.strip()
                if stem:
                    out.setdefault(stem, category)
        return out

    def _normalize_marker_hits_table(self, marker_hits_df: pd.DataFrame) -> pd.DataFrame:
        if marker_hits_df.empty:
            return pd.DataFrame(
                columns=[
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
                ]
            )
        required = {"Sample", "Contig", "DB"}
        if not required.issubset(set(marker_hits_df.columns)):
            return pd.DataFrame(columns=list(marker_hits_df.columns))

        out = marker_hits_df.copy()
        for col in [
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
        ]:
            if col not in out.columns:
                out[col] = pd.NA
        for col in ["AnnotationID", "Sample", "Contig", "DB"]:
            out[col] = self._clean_series(out[col])
        out = out[(out["Sample"] != "") & (out["Contig"] != "") & (out["DB"] != "")]
        return out.reset_index(drop=True)

    def _map_reference_db_name(self, value: str) -> str:
        db = self._clean_text(value)
        if not db:
            return ""
        if not self._configured_reference_dbs:
            return db
        if db in self._configured_reference_dbs:
            return db
        candidates = [name for name in self._configured_reference_dbs if name.startswith(f"{db}_")]
        if len(candidates) == 1:
            return candidates[0]
        return db

    def _build_reference_hits_table(self, ref_df: pd.DataFrame) -> pd.DataFrame:
        cols = ["Sample", "Contig", "DB"]
        if ref_df.empty or "Sample" not in ref_df.columns or "Contig" not in ref_df.columns:
            return pd.DataFrame(columns=cols)

        db_col = "DB" if "DB" in ref_df.columns else "Method" if "Method" in ref_df.columns else None
        if db_col is None:
            return pd.DataFrame(columns=cols)

        rows: List[Tuple[str, str, str]] = []
        for sample, contig, db_raw in ref_df[["Sample", "Contig", db_col]].itertuples(index=False, name=None):
            sample_name = self._clean_text(sample)
            contig_name = self._clean_text(contig)
            if not sample_name or not contig_name:
                continue
            mapped_db = self._map_reference_db_name(self._clean_text(db_raw))
            if not mapped_db:
                continue
            rows.append((sample_name, contig_name, mapped_db))
        if not rows:
            return pd.DataFrame(columns=cols)
        return pd.DataFrame(rows, columns=cols)

    def _collect_annotation_db_columns(self) -> List[str]:
        ordered: List[str] = []
        seen: set[str] = set()

        for db in [*self._configured_marker_dbs, *self._configured_reference_dbs]:
            name = self._clean_text(db)
            if not name or name in seen:
                continue
            seen.add(name)
            ordered.append(name)

        observed: List[str] = []
        if not self.marker_hits_df.empty and "DB" in self.marker_hits_df.columns:
            observed.extend(self.marker_hits_df["DB"].astype(str).tolist())
        if not self.reference_hits_df.empty and "DB" in self.reference_hits_df.columns:
            observed.extend(self.reference_hits_df["DB"].astype(str).tolist())

        for db in sorted({self._clean_text(v) for v in observed if self._clean_text(v)}):
            if db in seen:
                continue
            seen.add(db)
            ordered.append(db)

        return ordered

    def _build_resistance_atb_count_maps(self) -> Tuple[Dict[str, int], Dict[Tuple[str, str], int]]:
        sample_counts: Dict[str, int] = {}
        contig_counts: Dict[Tuple[str, str], int] = {}
        if self.marker_hits_df.empty or "ResistanceDrugClasses" not in self.marker_hits_df.columns:
            return sample_counts, contig_counts

        sample_classes: Dict[str, set[str]] = {}
        contig_classes: Dict[Tuple[str, str], set[str]] = {}
        cols = ["Sample", "Contig", "ResistanceDrugClasses"]
        for _i, (sample_raw, contig_raw, classes_raw) in enumerate(
            self.marker_hits_df[cols].itertuples(index=False, name=None)
        ):
            if (_i & 0x3FFF) == 0:  # keep the splash animating through this 500k+ row loop
                self._yield_to_ui()
            sample_name = self._clean_text(sample_raw)
            contig_name = self._clean_text(contig_raw)
            if not sample_name:
                continue
            classes = self._split_resistance_drug_classes(classes_raw)
            if not classes:
                continue
            sample_bucket = sample_classes.setdefault(sample_name, set())
            sample_bucket.update(classes)
            if contig_name:
                contig_bucket = contig_classes.setdefault((sample_name, contig_name), set())
                contig_bucket.update(classes)

        sample_counts = {sample: len(values) for sample, values in sample_classes.items()}
        contig_counts = {key: len(values) for key, values in contig_classes.items()}
        return sample_counts, contig_counts

    def _split_resistance_drug_classes(self, value: Any) -> List[str]:
        text = self._clean_text(value)
        if not text:
            return []
        out: List[str] = []
        seen: set[str] = set()
        for part in re.split(r"[;,]", text):
            label = self._format_resistance_drug_class_label(part)
            if not label:
                continue
            key = label.casefold()
            if key in seen:
                continue
            seen.add(key)
            out.append(label)
        return out

    def _format_resistance_drug_class_label(self, value: Any) -> str:
        label = self._clean_text(value).replace("_", " ")
        label = re.sub(r"\s+", " ", label).strip()
        if not label:
            return ""
        if label.upper() == label:
            return label.title()
        return label[0].upper() + label[1:]


    # Actions
    def refresh(self) -> None:
        self._refresh_generation += 1
        snapshot = RefreshSnapshot(
            token=self._refresh_generation,
            sample_query=self.sample_query.text().strip() or None,
            meta_filters={k: list(v) for k, v in self.metadata_filters.items()},
            sample_subset_samples=list(self._sample_subset_filter_samples),
            ann_filters=self._current_annotation_filters(),
            hide_empty_contigs=self.hide_empty_contigs_chk.isChecked(),
            view_mode=str(self.view_mode),
        )
        self._refresh_pending_snapshot = snapshot
        if not self._refresh_running:
            self._start_next_refresh()


    def _start_next_refresh(self) -> None:
        snapshot = self._refresh_pending_snapshot
        if snapshot is None:
            return
        self._refresh_pending_snapshot = None
        self._refresh_running = True
        worker = RefreshWorker(snapshot.token, lambda s=snapshot: self._compute_refresh_payload(s))
        worker.signals.finished.connect(self._on_refresh_payload_ready)
        worker.signals.failed.connect(self._on_refresh_payload_failed)
        self._refresh_worker = worker
        self._refresh_pool.start(worker)


    def _on_refresh_payload_ready(self, payload_obj: object) -> None:
        self._refresh_running = False
        self._refresh_worker = None
        if isinstance(payload_obj, RefreshPayload) and payload_obj.token == self._refresh_generation:
            self._apply_refresh_payload(payload_obj)
        self._start_next_refresh()


    def _on_refresh_payload_failed(self, token: int, error_text: str) -> None:
        self._refresh_running = False
        self._refresh_worker = None
        print(error_text)
        if token == self._refresh_generation:
            self.status.setText("Filter refresh failed. See terminal for traceback.")
        if not self._first_refresh_done:
            self._first_refresh_done = True
            self.firstRefreshFinished.emit()
        self._start_next_refresh()


    def _apply_refresh_payload(self, payload: RefreshPayload) -> None:
        self._base_allowed_samples = list(payload.allowed_samples)
        self._filtered_ann_pl = payload.filtered_ann_pl
        self._update_table_search_hint(payload.filtered_ann_pl)
        self._populate_metadata_tree(
            payload.meta_for_tree_pl,
            payload.contig_records,
            payload.filtered_counts,
            payload.tree_allowed_samples,
        )
        self._rerender_visible_table_columns()
        self.status.setText(payload.stats_text)
        self.ann_filters = dict(payload.ann_filters)
        self._apply_view_mode()
        # Keep the core results responsive, then let the heavier bottom panels catch up.
        self._schedule_deferred_panel_refresh(payload.panel_allowed_samples, self.ann_filters)
        self._apply_pan_visual_highlight()
        if not self._first_refresh_done:
            self._first_refresh_done = True
            self.firstRefreshFinished.emit()
