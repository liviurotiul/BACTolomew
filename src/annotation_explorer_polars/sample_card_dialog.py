from __future__ import annotations

import html
import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import pandas as pd
from PyQt6 import QtCore, QtWidgets

from pipeline.UserDatabases import USER_LABEL_SUFFIX

from .annotation_detail_dialog import _normalize_detail_value, _parse_float, _safe_value
from .resistance_mutation import MUTATION_TAG, resistance_mutations
from .windows.info_window import InfoWindow

_HIDDEN_INTERNAL_MARKER_FIELDS = {
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
}


class _AtbClassCard(QtWidgets.QFrame):
    def __init__(
        self,
        label: str,
        support_count: int,
        tooltip: str,
        warning: bool = False,
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self.setObjectName("AtbClassCard")
        self.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed,
        )
        self.setMinimumWidth(180)
        self.setMinimumHeight(62)
        self.setToolTip(tooltip)
        if warning:
            self.setStyleSheet(
                "#AtbClassCard {"
                "  border: 1px solid #efc1bc;"
                "  border-radius: 12px;"
                "  background: #fff2f0;"
                "}"
                "#AtbClassCard:hover {"
                "  border: 1px solid #e29b93;"
                "  background: #ffe7e3;"
                "}"
            )
        else:
            self.setStyleSheet(
                "#AtbClassCard {"
                "  border: 1px solid #d9e4ef;"
                "  border-radius: 12px;"
                "  background: #f8fbff;"
                "}"
                "#AtbClassCard:hover {"
                "  border: 1px solid #89aceb;"
                "  background: #eef4ff;"
                "}"
            )

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(10, 8, 10, 8)
        layout.setSpacing(2)

        title = QtWidgets.QLabel()
        title.setWordWrap(True)
        title.setStyleSheet("font-weight: 600; color: #203041;")
        title.setToolTip(tooltip)
        title.setAttribute(QtCore.Qt.WidgetAttribute.WA_TransparentForMouseEvents, True)
        title.setTextFormat(QtCore.Qt.TextFormat.RichText)
        title_html = html.escape(label)
        if warning:
            title_html = (
                f"{title_html} "
                "<span style='color:#d92d20; font-weight:700;'>!</span>"
            )
        title.setText(title_html)
        layout.addWidget(title)

        support_label = QtWidgets.QLabel(
            f"{int(support_count)} supporting gene{'s' if int(support_count) != 1 else ''}"
        )
        support_label.setStyleSheet(
            "color: #8b5249; font-size: 11px;" if warning else "color: #64778a; font-size: 11px;"
        )
        support_label.setToolTip(tooltip)
        support_label.setAttribute(QtCore.Qt.WidgetAttribute.WA_TransparentForMouseEvents, True)
        layout.addWidget(support_label)


class SampleCardDialog(InfoWindow):
    """The sample card (an info window, plan section 5, D10): a read-only card
    with its Summary and Raw tabs and nothing else. It has no commands, no
    buttons that filter, reveal or open anything, and no hover output."""

    COMMANDS = ()

    def __init__(
        self,
        parent: Optional[QtWidgets.QWidget],
        *,
        sample: str,
        ann_df: Optional[pd.DataFrame] = None,
        meta_df: Optional[pd.DataFrame] = None,
        marker_hits_df: Optional[pd.DataFrame] = None,
        reference_hits_df: Optional[pd.DataFrame] = None,
        output_dir: Optional[Path] = None,
    ) -> None:
        super().__init__(parent)
        self._sample = str(sample or "").strip()
        self._ann_df = ann_df if ann_df is not None else pd.DataFrame()
        self._meta_df = meta_df if meta_df is not None else pd.DataFrame()
        self._marker_hits_df = marker_hits_df if marker_hits_df is not None else pd.DataFrame()
        self._reference_hits_df = reference_hits_df if reference_hits_df is not None else pd.DataFrame()
        self._output_dir = Path(output_dir) if output_dir is not None else None
        self._descriptor_labels: List[QtWidgets.QLabel] = []

        self._sample_ann = self._filter_sample_frame(self._ann_df)
        self._sample_marker_hits = self._filter_sample_frame(self._marker_hits_df)
        self._sample_reference_hits = self._filter_sample_frame(self._reference_hits_df)
        self._sample_meta = self._sample_metadata_map()
        self._resistance_ann = self._build_resistance_annotations()
        self._mutation_ann = self._build_mutation_annotations()
        self._atb_cards = self._build_atb_cards()

        self.set_card_title(f"Sample card: {self._sample}" if self._sample else "Sample card")
        self.resize(980, 720)
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)

        tabs = QtWidgets.QTabWidget()
        tabs.addTab(self._build_summary_tab(), "Summary")
        tabs.addTab(self._build_raw_tab(), "Raw")
        layout.addWidget(tabs)

        btns = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Close)
        btns.rejected.connect(self.reject)
        btns.accepted.connect(self.accept)
        layout.addWidget(btns)

    def _build_summary_tab(self) -> QtWidgets.QWidget:
        widget = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(widget)

        header = QtWidgets.QLabel(self._format_header())
        header.setWordWrap(True)
        header.setStyleSheet("font-size: 18px; font-weight: 600; color: #203041;")
        layout.addWidget(header)

        subtitle = QtWidgets.QLabel(self._format_subtitle())
        subtitle.setWordWrap(True)
        subtitle.setStyleSheet("color: #6a7b8d;")
        layout.addWidget(subtitle)

        stats_group = QtWidgets.QGroupBox("Overview")
        stats_layout = QtWidgets.QFormLayout(stats_group)
        for label, value in self._summary_rows():
            stats_layout.addRow(self._make_descriptor_label(label), QtWidgets.QLabel(value))
        layout.addWidget(stats_group)

        meta_group = QtWidgets.QGroupBox("Metadata")
        meta_layout = QtWidgets.QFormLayout(meta_group)
        if self._sample_meta:
            for key, value in self._sample_meta.items():
                meta_layout.addRow(self._make_descriptor_label(key), self._build_wrapped_value(value))
        else:
            meta_layout.addRow(self._make_descriptor_label("(none)"), QtWidgets.QLabel("No metadata available"))
        layout.addWidget(meta_group)

        layout.addWidget(self._build_atb_group())

        db_summary_group = QtWidgets.QGroupBox("Database Hit Summary")
        db_summary_layout = QtWidgets.QVBoxLayout(db_summary_group)
        summary_table = self._build_db_summary_table()
        if summary_table is None:
            db_summary_layout.addWidget(QtWidgets.QLabel("No marker or reference hits found for this sample."))
        else:
            db_summary_layout.addWidget(summary_table)
        layout.addWidget(db_summary_group)

        layout.addStretch(1)
        self._update_descriptor_widths()
        return self._wrap_scroll(widget)

    def _build_raw_tab(self) -> QtWidgets.QWidget:
        widget = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(widget)
        raw = QtWidgets.QPlainTextEdit()
        raw.setReadOnly(True)
        raw.setPlainText(
            json.dumps(self._build_raw_payload(), indent=2, ensure_ascii=False)
        )
        layout.addWidget(raw)
        return widget

    def _build_wrapped_value(self, value: str) -> QtWidgets.QLabel:
        label = QtWidgets.QLabel(value)
        label.setWordWrap(True)
        label.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        return label

    def _format_header(self) -> str:
        sample = self._sample or "Unknown sample"
        if self._output_dir is None:
            return sample
        return f"{sample} | {self._output_dir.name}"

    def _format_subtitle(self) -> str:
        parts = [
            f"{len(self._sample_ann)} annotations",
            f"{self._count_unique_values(self._sample_ann, 'Contig')} contigs",
            f"{len(self._resistance_ann)} resistance-tagged annotations",
        ]
        if not self._mutation_ann.empty:
            parts.append(f"{len(self._mutation_symbols())} resistance mutations")
        if self._atb_cards:
            parts.append(f"{len(self._atb_cards)} ATB classes")
        if not self._sample_marker_hits.empty:
            parts.append(f"{len(self._sample_marker_hits)} marker hits")
        if not self._sample_reference_hits.empty:
            parts.append(f"{len(self._sample_reference_hits)} reference-locus hits")
        return " | ".join(parts)

    def _summary_rows(self) -> List[Tuple[str, str]]:
        return [
            ("Sample", self._sample or "Unknown"),
            ("Annotations", str(len(self._sample_ann))),
            ("Contigs", str(self._count_unique_values(self._sample_ann, "Contig"))),
            ("Unique genes", str(self._count_unique_values(self._sample_ann, "Gene"))),
            ("Resistance genes", str(len(self._resistance_ann))),
            ("Resistance mutations", self._mutation_summary()),
            ("ATB classes", str(len(self._atb_cards))),
            ("Marker hits", str(len(self._sample_marker_hits))),
            ("Marker DBs", str(self._count_unique_values(self._sample_marker_hits, "DB"))),
            ("Reference loci hits", str(len(self._sample_reference_hits))),
            ("Reference DBs", str(self._count_unique_values(self._sample_reference_hits, "DB"))),
            ("Metadata fields", str(len(self._sample_meta))),
        ]

    def _build_atb_group(self) -> QtWidgets.QGroupBox:
        group = QtWidgets.QGroupBox("ATB")
        layout = QtWidgets.QVBoxLayout(group)

        hint = QtWidgets.QLabel("Hover a resistance class to see the supporting genes and marker evidence.")
        hint.setWordWrap(True)
        hint.setStyleSheet("color: #6a7b8d;")
        layout.addWidget(hint)

        if not self._atb_cards:
            empty = QtWidgets.QLabel("No antibiotic resistance classes identified for this sample.")
            empty.setWordWrap(True)
            layout.addWidget(empty)
            return group

        cards_wrap = QtWidgets.QWidget()
        cards_grid = QtWidgets.QGridLayout(cards_wrap)
        cards_grid.setContentsMargins(0, 0, 0, 0)
        cards_grid.setHorizontalSpacing(8)
        cards_grid.setVerticalSpacing(8)
        max_cols = 4
        for idx, card_data in enumerate(self._atb_cards):
            card = _AtbClassCard(
                card_data["label"],
                int(card_data["support_count"]),
                str(card_data["tooltip"]),
                bool(card_data.get("warning")),
            )
            row = idx // max_cols
            col = idx % max_cols
            cards_grid.addWidget(card, row, col)
        for col in range(max_cols):
            cards_grid.setColumnStretch(col, 1)
        cards_grid.setRowStretch((len(self._atb_cards) // max_cols) + 1, 1)
        layout.addWidget(cards_wrap)
        return group

    def _build_resistance_annotations(self) -> pd.DataFrame:
        df = self._sample_ann
        if df.empty or "Tags" not in df.columns:
            return pd.DataFrame(columns=list(df.columns))
        tags = df["Tags"].astype(str).str.upper()
        out = df[tags.str.contains("A", na=False)].copy()
        if out.empty:
            return out
        sort_cols = [c for c in ["Contig", "StartInt", "StopInt", "Gene", "Product"] if c in out.columns]
        if sort_cols:
            out = out.sort_values(sort_cols, na_position="last")
        return out.reset_index(drop=True)

    def _build_mutation_annotations(self) -> pd.DataFrame:
        """CDS of this sample with an AMRFinderPlus point mutation (tag X)."""
        df = self._sample_ann
        if df.empty:
            return pd.DataFrame(columns=list(df.columns))
        mask = pd.Series(False, index=df.index)
        if "ResistanceMutations" in df.columns:
            mask |= df["ResistanceMutations"].fillna("").astype(str).str.strip() != ""
        if "Tags" in df.columns:
            mask |= df["Tags"].astype(str).str.upper().str.contains(MUTATION_TAG, regex=False, na=False)
        return df[mask].reset_index(drop=True)

    def _mutation_symbols(self) -> List[str]:
        symbols: List[str] = []
        for row in self._mutation_ann.to_dict("records"):
            for symbol in resistance_mutations(row) or [_safe_value(row.get("Gene")) or "unnamed"]:
                if symbol not in symbols:
                    symbols.append(symbol)
        return symbols

    def _mutation_summary(self) -> str:
        symbols = self._mutation_symbols()
        return f"{len(symbols)} ({', '.join(symbols)})" if symbols else "0"

    def _atb_sources(self) -> List[Tuple[List[str], Dict[str, Any]]]:
        """(drug classes, support row) from the marker hits and from the resistance mutations."""
        out: List[Tuple[List[str], Dict[str, Any]]] = []
        source_df = self._sample_marker_hits
        if not source_df.empty and "ResistanceDrugClasses" in source_df.columns:
            for row in source_df.to_dict("records"):
                classes = self._split_atb_classes(row.get("ResistanceDrugClasses"))
                support_row = self._build_atb_support_row(row) if classes else None
                if support_row is not None:
                    out.append((classes, support_row))
        for row in self._mutation_ann.to_dict("records"):
            classes = self._split_atb_classes(row.get("ResistanceMutationClasses"))
            for symbol in resistance_mutations(row) if classes else []:
                out.append((classes, {
                    "GeneSymbol": symbol,
                    "MarkerAnnotation": "AMRFinderPlus mutation",
                    "Subject_coverage": "",
                    "Identity": "",
                    "_CoverageSort": -1.0,
                    "_IdentitySort": -1.0,
                }))
        return out

    def _build_atb_cards(self) -> List[Dict[str, Any]]:
        grouped: Dict[str, Dict[str, Any]] = {}
        for classes, support_row in self._atb_sources():
            for class_name in classes:
                display = self._format_atb_label(class_name)
                if not display:
                    continue
                key = display.casefold()
                bucket = grouped.setdefault(
                    key,
                    {
                        "label": display,
                        "genes": [],
                        "gene_keys": set(),
                        "support_rows": {},
                    },
                )
                gene = str(support_row["GeneSymbol"])
                gene_key = gene.casefold()
                if gene_key not in bucket["gene_keys"]:
                    bucket["gene_keys"].add(gene_key)
                    bucket["genes"].append(gene)
                row_key = (
                    str(support_row["GeneSymbol"]).casefold(),
                    str(support_row["MarkerAnnotation"]).casefold(),
                )
                current = bucket["support_rows"].get(row_key)
                if current is None or self._atb_support_row_sort_key(support_row) > self._atb_support_row_sort_key(current):
                    bucket["support_rows"][row_key] = support_row

        cards: List[Dict[str, Any]] = []
        for bucket in grouped.values():
            genes = sorted(bucket["genes"], key=str.casefold)
            if not genes:
                continue
            support_rows = sorted(
                bucket["support_rows"].values(),
                key=lambda row: (
                    str(row.get("GeneSymbol") or "").casefold(),
                    str(row.get("MarkerAnnotation") or "").casefold(),
                    -float(row.get("_CoverageSort") or -1.0),
                    -float(row.get("_IdentitySort") or -1.0),
                ),
            )
            cards.append(
                {
                    "label": bucket["label"],
                    "support_count": len(genes),
                    "supporting_genes": genes,
                    "support_rows": support_rows,
                    "warning": self._atb_card_has_warning(len(genes), support_rows),
                    "tooltip": self._format_atb_tooltip(str(bucket["label"]), support_rows),
                }
            )
        cards.sort(key=lambda card: (-int(card["support_count"]), str(card["label"]).casefold()))
        return cards

    def _split_atb_classes(self, value: Any) -> List[str]:
        text = _safe_value(value)
        if not text:
            return []
        out: List[str] = []
        seen = set()
        for part in re.split(r"[;,]", text):
            label = self._format_atb_label(part)
            if not label:
                continue
            key = label.casefold()
            if key in seen:
                continue
            seen.add(key)
            out.append(label)
        return out

    def _canonical_atb_label(self, value: Any) -> str:
        label = _safe_value(value)
        if not label:
            return ""
        label = re.sub(r"\s+", " ", label.replace("_", " ")).strip()
        if not label:
            return ""
        # CARD often uses labels like "X antibiotic" while AMRF uses just "X".
        # Collapse those into one canonical class name before grouping/counting.
        label = re.sub(r"(?i)\s+antibiotics?$", "", label).strip()
        label = re.sub(r"(?i)\s+antimicrobial(s)?$", "", label).strip()
        label = re.sub(r"(?i)\s+agents?$", "", label).strip()
        label = re.sub(r"\s+", " ", label).strip(" -_/;,")
        return label

    def _format_atb_label(self, value: Any) -> str:
        label = self._canonical_atb_label(value)
        if not label:
            return ""
        if label.upper() == label:
            return label.title()
        return label[0].upper() + label[1:]

    def _supporting_genes_for_row(self, row: Dict[str, Any]) -> List[str]:
        genes = self._extract_gene_tokens(row.get("ResistanceAlleles"))
        if not genes:
            genes = self._extract_gene_tokens(row.get("ResistanceGeneFamilies"))

        gene_name = _safe_value(row.get("Gene"))
        if gene_name:
            genes.insert(0, gene_name)

        out: List[str] = []
        seen = set()
        for gene in genes:
            label = re.sub(r"\s+", " ", str(gene).strip())
            if not label:
                continue
            key = label.casefold()
            if key in seen:
                continue
            seen.add(key)
            out.append(label)
        return out

    def _preferred_gene_symbol(self, row: Dict[str, Any]) -> str:
        gene_name = _safe_value(row.get("Gene"))
        if gene_name:
            return gene_name
        genes = self._supporting_genes_for_row(row)
        if genes:
            return str(genes[0])
        product = _safe_value(row.get("Product"))
        if product:
            return product
        return ""

    def _format_marker_db_label(self, value: Any) -> str:
        db = _safe_value(value)
        if not db:
            return ""
        mapping = {
            "CARD_AMR": "CARD",
            "AMRF_core-AMR": "AMRFinder core",
            "AMRF_plus-AMR": "AMRFinder plus",
        }
        if db in mapping:
            return mapping[db]
        if db.endswith(USER_LABEL_SUFFIX):  # a user database keeps its name
            return db
        return db.replace("_", " ")

    def _build_atb_support_row(self, row: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        gene_symbol = self._preferred_gene_symbol(row)
        if not gene_symbol:
            return None
        marker_annotation = self._format_marker_db_label(row.get("DB"))
        coverage = _parse_float(row.get("Subject_coverage"))
        identity = _parse_float(row.get("Identity"))
        return {
            "GeneSymbol": gene_symbol,
            "MarkerAnnotation": marker_annotation or _safe_value(row.get("DB")),
            "Subject_coverage": "" if coverage is None else f"{coverage:.1f}",
            "Identity": "" if identity is None else f"{identity:.1f}",
            "_CoverageSort": -1.0 if coverage is None else float(coverage),
            "_IdentitySort": -1.0 if identity is None else float(identity),
        }

    def _atb_support_row_sort_key(self, row: Dict[str, Any]) -> Tuple[float, float]:
        return (
            float(row.get("_CoverageSort") or -1.0),
            float(row.get("_IdentitySort") or -1.0),
        )

    def _support_row_is_under_full_match(self, row: Dict[str, Any]) -> bool:
        coverage = _parse_float(row.get("_CoverageSort"))
        identity = _parse_float(row.get("_IdentitySort"))
        cov_under = coverage is not None and 0.0 <= coverage < 100.0
        id_under = identity is not None and 0.0 <= identity < 100.0
        return cov_under or id_under

    def _atb_card_has_warning(
        self,
        support_count: int,
        support_rows: Sequence[Dict[str, Any]],
    ) -> bool:
        if int(support_count) != 1:
            return False
        return any(self._support_row_is_under_full_match(row) for row in support_rows)

    def _extract_gene_tokens(self, value: Any) -> List[str]:
        text = _safe_value(value)
        if not text:
            return []
        out: List[str] = []
        for part in text.split(";"):
            token = str(part).strip()
            if not token:
                continue
            if "," in token:
                token = token.split(",", 1)[0].strip()
            token = re.sub(r"\s+", " ", token)
            if token:
                out.append(token)
        return out

    def _format_atb_tooltip(
        self,
        label: str,
        support_rows: Sequence[Dict[str, Any]],
    ) -> str:
        shown = list(support_rows[:24])
        hidden_count = max(0, len(support_rows) - len(shown))
        lines = [
            f"<b>{html.escape(label)}</b><br>",
            f"{len(support_rows)} supporting marker annotation{'s' if len(support_rows) != 1 else ''}<br><br>",
            "<table cellspacing='0' cellpadding='3'>",
            "<tr>"
            "<td><b>Gene</b></td>"
            "<td><b>DB</b></td>"
            "<td><b>ID</b></td>"
            "<td><b>COV</b></td>"
            "</tr>",
        ]
        for row in shown:
            row_warn = self._support_row_is_under_full_match(row)
            gene_symbol = self._tooltip_metric_span(
                str(row.get("GeneSymbol") or ""),
                highlight=row_warn,
            )
            marker_annotation = html.escape(str(row.get("MarkerAnnotation") or ""))
            identity = self._tooltip_metric_span(
                self._format_percent_text(row.get("Identity")),
                highlight=self._metric_text_under_100(row.get("_IdentitySort")),
            )
            coverage = self._tooltip_metric_span(
                self._format_percent_text(row.get("Subject_coverage")),
                highlight=self._metric_text_under_100(row.get("_CoverageSort")),
            )
            lines.append(
                "<tr>"
                f"<td>{gene_symbol}</td>"
                f"<td>{marker_annotation}</td>"
                f"<td>{identity}</td>"
                f"<td>{coverage}</td>"
                "</tr>"
            )
        lines.append("</table>")
        if hidden_count:
            lines.append(f"<br>+ {hidden_count} more")
        return "".join(lines)

    def _metric_text_under_100(self, value: Any) -> bool:
        metric = _parse_float(value)
        return metric is not None and 0.0 <= metric < 100.0

    def _format_percent_text(self, value: Any) -> str:
        metric = _parse_float(value)
        if metric is None:
            return html.escape(_safe_value(value))
        return f"{metric:.1f}%"

    def _tooltip_metric_span(self, text: str, *, highlight: bool) -> str:
        safe = html.escape(str(text or ""))
        if not highlight:
            return safe
        return f"<span style='color:#f0c94d; font-weight:700;'>{safe}</span>"

    def _build_db_summary_table(self) -> Optional[QtWidgets.QTableWidget]:
        records: List[Dict[str, Any]] = []
        if not self._sample_marker_hits.empty and "DB" in self._sample_marker_hits.columns:
            marker_summary = (
                self._sample_marker_hits.groupby("DB", dropna=False)
                .agg(Hits=("DB", "size"), Contigs=("Contig", "nunique"))
                .reset_index()
            )
            for row in marker_summary.to_dict("records"):
                records.append(
                    {
                        "Type": "Marker",
                        "DB": _safe_value(row.get("DB")),
                        "Hits": int(row.get("Hits") or 0),
                        "Contigs": int(row.get("Contigs") or 0),
                    }
                )
        if not self._sample_reference_hits.empty and "DB" in self._sample_reference_hits.columns:
            reference_summary = (
                self._sample_reference_hits.groupby("DB", dropna=False)
                .agg(Hits=("DB", "size"), Contigs=("Contig", "nunique"))
                .reset_index()
            )
            for row in reference_summary.to_dict("records"):
                records.append(
                    {
                        "Type": "Reference",
                        "DB": _safe_value(row.get("DB")),
                        "Hits": int(row.get("Hits") or 0),
                        "Contigs": int(row.get("Contigs") or 0),
                    }
                )
        if not records:
            return None
        records.sort(key=lambda row: (str(row.get("Type") or ""), str(row.get("DB") or "")))
        columns = [
            ("Type", "Type"),
            ("DB", "DB"),
            ("Hits", "Hits"),
            ("Contigs", "Contigs"),
        ]
        return self._build_table_from_records(records, columns)

    def _build_table_from_records(
        self,
        records: Sequence[Dict[str, Any]],
        columns: Sequence[Tuple[str, str]],
    ) -> QtWidgets.QTableWidget:
        table = QtWidgets.QTableWidget()
        table.setColumnCount(len(columns))
        table.setHorizontalHeaderLabels([label for _, label in columns])
        table.setRowCount(len(records))
        table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        table.setAlternatingRowColors(True)
        table.verticalHeader().setVisible(False)
        table.horizontalHeader().setStretchLastSection(True)

        for row_idx, record in enumerate(records):
            for col_idx, (key, _label) in enumerate(columns):
                value = record.get(key)
                item = QtWidgets.QTableWidgetItem(_safe_value(value))
                table.setItem(row_idx, col_idx, item)

        header = table.horizontalHeader()
        for idx in range(max(0, len(columns) - 1)):
            header.setSectionResizeMode(idx, QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        if columns:
            header.setSectionResizeMode(len(columns) - 1, QtWidgets.QHeaderView.ResizeMode.Stretch)
        table.resizeRowsToContents()
        table.setSortingEnabled(True)
        return table

    def _build_raw_payload(self) -> Dict[str, Any]:
        ann_cols = [c for c in self._sample_ann.columns if c not in _HIDDEN_INTERNAL_MARKER_FIELDS]
        resistance_cols = [c for c in self._resistance_ann.columns if c not in _HIDDEN_INTERNAL_MARKER_FIELDS]
        marker_cols = [c for c in self._sample_marker_hits.columns if c not in _HIDDEN_INTERNAL_MARKER_FIELDS]
        return {
            "sample": self._sample,
            "output_dir": str(self._output_dir) if self._output_dir is not None else None,
            "summary": {
                "annotations": len(self._sample_ann),
                "contigs": self._count_unique_values(self._sample_ann, "Contig"),
                "unique_genes": self._count_unique_values(self._sample_ann, "Gene"),
                "resistance_annotations": len(self._resistance_ann),
                "resistance_mutations": self._mutation_symbols(),
                "atb_classes": len(self._atb_cards),
                "marker_hits": len(self._sample_marker_hits),
                "reference_hits": len(self._sample_reference_hits),
            },
            "atb": _normalize_detail_value(
                [
                    {
                        "class": card["label"],
                        "support_count": card["support_count"],
                        "supporting_genes": card["supporting_genes"],
                        "support_rows": [
                            {
                                "GeneSymbol": row.get("GeneSymbol"),
                                "MarkerAnnotation": row.get("MarkerAnnotation"),
                                "Subject_coverage": row.get("Subject_coverage"),
                                "Identity": row.get("Identity"),
                            }
                            for row in card.get("support_rows", [])
                        ],
                    }
                    for card in self._atb_cards
                ]
            ),
            "metadata": _normalize_detail_value(self._sample_meta),
            "annotations": _normalize_detail_value(self._sample_ann[ann_cols].to_dict("records")),
            "resistance_annotations": _normalize_detail_value(self._resistance_ann[resistance_cols].to_dict("records")),
            "marker_hits": _normalize_detail_value(self._sample_marker_hits[marker_cols].to_dict("records")),
            "reference_hits": _normalize_detail_value(self._sample_reference_hits.to_dict("records")),
        }

    def _wrap_scroll(self, widget: QtWidgets.QWidget) -> QtWidgets.QScrollArea:
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        scroll.setWidget(widget)
        return scroll

    def _filter_sample_frame(self, df: Optional[pd.DataFrame]) -> pd.DataFrame:
        if df is None or df.empty or "Sample" not in df.columns:
            return pd.DataFrame(columns=list(df.columns) if isinstance(df, pd.DataFrame) else [])
        mask = df["Sample"].astype(str).str.strip() == self._sample
        out = df[mask].copy()
        if out.empty:
            return out
        sort_cols = [c for c in ["Contig", "StartInt", "StopInt", "Gene", "Product", "DB"] if c in out.columns]
        if sort_cols:
            out = out.sort_values(sort_cols, na_position="last")
        return out.reset_index(drop=True)

    def _sample_metadata_map(self) -> Dict[str, str]:
        out: Dict[str, str] = {}
        df = self._meta_df
        if df is None or df.empty or "Sample" not in df.columns:
            return out
        rows = df[df["Sample"].astype(str).str.strip() == self._sample]
        if rows.empty:
            return out
        row = rows.iloc[0]
        for col in rows.columns:
            if col == "Sample":
                continue
            value = _safe_value(row.get(col))
            if value:
                out[str(col)] = value
        return out

    def _count_unique_values(self, df: pd.DataFrame, column: str) -> int:
        if df.empty or column not in df.columns:
            return 0
        values = {
            _safe_value(value)
            for value in df[column].tolist()
            if _safe_value(value)
        }
        return len(values)

    def _make_descriptor_label(self, text: str) -> QtWidgets.QLabel:
        label = QtWidgets.QLabel(text)
        label.setStyleSheet("font-weight: 600; color: #2b3c4d;")
        label.setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignVCenter)
        self._descriptor_labels.append(label)
        return label

    def _update_descriptor_widths(self) -> None:
        if not self._descriptor_labels:
            return
        width = max(120, int(self.width() * 0.16))
        for label in self._descriptor_labels:
            label.setFixedWidth(width)

    def resizeEvent(self, event) -> None:  # type: ignore[override]
        super().resizeEvent(event)
        self._update_descriptor_widths()
