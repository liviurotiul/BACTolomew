from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from .context_model import ContextGene, ContextRenderModel, ContextRow


@dataclass(frozen=True)
class GeneLayout:
    row_index: int
    visual_index: int
    gene: ContextGene
    left: float
    right: float
    center_y: float
    display_strand: str

    @property
    def center_x(self) -> float:
        return (self.left + self.right) / 2.0

    @property
    def width(self) -> float:
        return max(0.0, self.right - self.left)


@dataclass(frozen=True)
class RowLayout:
    row: ContextRow
    row_index: int
    center_y: float
    pitch: float
    gene_height: float
    flipped: bool
    genes: Tuple[GeneLayout, ...]
    min_x: float
    max_x: float
    # A thin row for a contig with only the focal gene (shown on request).
    short: bool = False
    # The contig really ends at this end of the row (after mirroring).
    cap_left: bool = False
    cap_right: bool = False

    def gene_at(
        self,
        x: float,
        y: float,
        *,
        x_tolerance: float = 3.0,
        y_tolerance: float = 3.0,
    ) -> Optional[GeneLayout]:
        if abs(y - self.center_y) > (self.gene_height / 2.0 + y_tolerance):
            return None
        # Rows generally contain only tens of genes.  Keeping the tuple sorted
        # makes this inexpensive and avoids thousands of QGraphicsItems.
        for gene in self.genes:
            if gene.left - x_tolerance <= x <= gene.right + x_tolerance:
                return gene
            if gene.left > x + x_tolerance:
                break
        return None


@dataclass(frozen=True)
class LinkLayout:
    upper_row: int
    lower_row: int
    upper_gene: GeneLayout
    lower_gene: GeneLayout
    family_key: str
    evidence: str
    evidence_weight: float
    inverted: bool

    @property
    def label(self) -> str:
        return self.upper_gene.gene.label or self.family_key.split(":", 1)[-1]


@dataclass(frozen=True)
class PairLayout:
    upper_row: int
    lower_row: int
    links: Tuple[LinkLayout, ...]


@dataclass(frozen=True)
class DividerLayout:
    """Labelled band between the plasmid and chromosome blocks."""

    center_y: float
    height: float
    molecule: str
    text: str
    # Contigs of this type folded into the summary line at the bottom.
    folded: int = 0


@dataclass(frozen=True)
class FoldLayout:
    """The summary line standing for contigs with only the focal gene."""

    center_y: float
    height: float
    rows: Tuple[ContextRow, ...]
    expanded: bool

    @property
    def contig_count(self) -> int:
        return sum(contig_count(row) for row in self.rows)

    def molecule_counts(self) -> Dict[str, int]:
        counts: Dict[str, int] = {}
        for row in self.rows:
            for kind, count in row_molecule_counts(row).items():
                counts[kind] = counts.get(kind, 0) + count
        return counts


def contig_count(row: ContextRow) -> int:
    """Contigs a row stands for: several when identical rows are merged."""
    try:
        merged = int(row.source.get("collapsed_track_count") or 0)
    except (TypeError, ValueError):
        merged = 0
    return max(1, merged)


def row_molecule_counts(row: ContextRow) -> Dict[str, int]:
    """Contigs per molecule type behind a row.

    A merged row of plasmid and chromosome contigs is drawn as 'mixed' but
    carries its members' types, so counts never lump them under 'mixed'.
    """
    raw = row.source.get("collapsed_molecule_counts")
    counts: Dict[str, int] = {}
    if isinstance(raw, dict):
        for kind, value in raw.items():
            try:
                number = int(value)
            except (TypeError, ValueError):
                continue
            if number > 0:
                counts[str(kind)] = counts.get(str(kind), 0) + number
    return counts or {row.molecule: contig_count(row)}


MOLECULE_BLOCK_NAMES = {
    "plasmid": "Plasmid contigs",
    "chromosome": "Chromosome contigs",
    "mixed": "Mixed rows",
    "unknown": "Contigs of unknown type",
}


@dataclass(frozen=True)
class ContextLayout:
    model: ContextRenderModel
    rows: Tuple[RowLayout, ...]
    pairs: Tuple[PairLayout, ...]
    min_x: float
    max_x: float
    min_y: float
    max_y: float
    px_per_bp: float
    focal_x: float = 0.0
    dividers: Tuple[DividerLayout, ...] = ()
    fold: Optional[FoldLayout] = None

    @property
    def width(self) -> float:
        return max(1.0, self.max_x - self.min_x)

    @property
    def height(self) -> float:
        return max(1.0, self.max_y - self.min_y)

    def row_at(self, y: float) -> Optional[RowLayout]:
        for row in self.rows:
            if abs(y - row.center_y) <= row.pitch / 2.0:
                return row
        return None

    def gene_at(
        self,
        x: float,
        y: float,
        *,
        x_tolerance: float = 3.0,
        y_tolerance: float = 3.0,
    ) -> Optional[GeneLayout]:
        row = self.row_at(y)
        return (
            row.gene_at(
                x,
                y,
                x_tolerance=x_tolerance,
                y_tolerance=y_tolerance,
            )
            if row is not None
            else None
        )

    @property
    def reference_row(self) -> Optional[RowLayout]:
        for row in self.rows:
            if row.row.is_reference:
                return row
        return self.rows[0] if self.rows else None


class ContextLayoutEngine:
    """Pure layout for the interactive canvas and document exports."""

    def __init__(
        self,
        *,
        px_per_bp: float = 0.065,
        compact_pitch: float = 50.0,
        comfortable_pitch: float = 66.0,
        compact_gene_height: float = 16.0,
        comfortable_gene_height: float = 18.0,
        top_margin: float = 34.0,
        bottom_margin: float = 28.0,
        horizontal_margin: float = 34.0,
        minimum_link_evidence: float = 0.90,
        short_pitch: float = 24.0,
        short_gene_height: float = 10.0,
        divider_height: float = 26.0,
        fold_height: float = 30.0,
    ) -> None:
        self.px_per_bp = float(px_per_bp)
        self.compact_pitch = float(compact_pitch)
        self.comfortable_pitch = float(comfortable_pitch)
        self.compact_gene_height = float(compact_gene_height)
        self.comfortable_gene_height = float(comfortable_gene_height)
        self.top_margin = float(top_margin)
        self.bottom_margin = float(bottom_margin)
        self.horizontal_margin = float(horizontal_margin)
        self.minimum_link_evidence = float(minimum_link_evidence)
        self.short_pitch = float(short_pitch)
        self.short_gene_height = float(short_gene_height)
        self.divider_height = float(divider_height)
        self.fold_height = float(fold_height)

    @staticmethod
    def _should_flip(row: ContextRow, normalize_orientation: bool) -> bool:
        if not normalize_orientation or row.is_reference:
            return False
        if row.comparison_orientation == "reverse":
            return True
        if row.comparison_orientation == "forward":
            return False
        focal = row.focal_gene
        return bool(focal is not None and focal.strand == "-")

    @staticmethod
    def _display_strand(strand: str, flipped: bool) -> str:
        if not flipped:
            return strand
        return "+" if strand == "-" else "-"

    def _layout_row(
        self,
        row: ContextRow,
        row_index: int,
        *,
        pitch: float,
        gene_height: float,
        normalize_orientation: bool,
        center_y: Optional[float] = None,
        short: bool = False,
    ) -> RowLayout:
        flipped = self._should_flip(row, normalize_orientation)
        anchor = row.anchor_bp
        if row.center_bp is None and row.genes and not any(gene.is_focal for gene in row.genes):
            # No focal gene (plasmid context): start every row at its first
            # gene, instead of at a middle gene that differs between rows.
            anchor = max(gene.high for gene in row.genes) if flipped else min(gene.low for gene in row.genes)
        if center_y is None:
            center_y = self.top_margin + row_index * pitch
        genes = []
        for gene in row.genes:
            if flipped:
                rel_left = anchor - gene.high
                rel_right = anchor - gene.low
            else:
                rel_left = gene.low - anchor
                rel_right = gene.high - anchor
            genes.append(
                GeneLayout(
                    row_index=row_index,
                    visual_index=gene.index,
                    gene=gene,
                    left=rel_left * self.px_per_bp,
                    right=rel_right * self.px_per_bp,
                    center_y=center_y,
                    display_strand=self._display_strand(gene.strand, flipped),
                )
            )
        genes.sort(key=lambda gene: (gene.left, gene.right, gene.visual_index))
        min_x = min((gene.left for gene in genes), default=-1.0)
        max_x = max((gene.right for gene in genes), default=1.0)
        return RowLayout(
            row=row,
            row_index=row_index,
            center_y=center_y,
            pitch=pitch,
            gene_height=gene_height,
            flipped=flipped,
            genes=tuple(genes),
            min_x=min_x,
            max_x=max_x,
            short=short,
            # Mirroring swaps which contig end sits on the left.
            cap_left=row.ends_high if flipped else row.ends_low,
            cap_right=row.ends_low if flipped else row.ends_high,
        )

    @staticmethod
    def _divider_text(molecule: str, rows: Sequence[ContextRow], folded: int) -> str:
        name = MOLECULE_BLOCK_NAMES.get(molecule, "Contigs")
        text = f"{name}, {sum(contig_count(row) for row in rows)}"
        if folded:
            text += f" (+{folded} folded)"
        return text

    def _molecule_blocks(
        self,
        rows: Sequence[ContextRow],
        folded: Sequence[ContextRow],
    ) -> Dict[int, Tuple[str, int]]:
        """Divider text keyed by the index of the row that opens each block.

        Rows arrive already ordered by molecule.  A reference row at the top
        stays above every block.
        """
        dividers: Dict[int, Tuple[str, int]] = {}
        start = 1 if rows and rows[0].is_reference else 0
        index = start
        while index < len(rows):
            molecule = rows[index].molecule
            end = index
            while end < len(rows) and rows[end].molecule == molecule:
                end += 1
            hidden = sum(row_molecule_counts(row).get(molecule, 0) for row in folded)
            dividers[index] = (self._divider_text(molecule, rows[index:end], hidden), hidden)
            index = end
        return dividers

    @staticmethod
    def _pair_occurrences(
        upper_genes: Sequence[GeneLayout],
        lower_genes: Sequence[GeneLayout],
    ) -> Tuple[Tuple[GeneLayout, GeneLayout], ...]:
        """Choose a monotonic, minimum-distance one-to-one repeat matching."""
        if not upper_genes or not lower_genes:
            return ()
        upper_sorted = sorted(upper_genes, key=lambda gene: gene.center_x)
        lower_sorted = sorted(lower_genes, key=lambda gene: gene.center_x)
        upper_is_longer = len(upper_sorted) >= len(lower_sorted)
        longer = upper_sorted if upper_is_longer else lower_sorted
        shorter = lower_sorted if upper_is_longer else upper_sorted
        n, m = len(longer), len(shorter)
        infinity = float("inf")
        costs = [[infinity] * (m + 1) for _ in range(n + 1)]
        matched = [[False] * (m + 1) for _ in range(n + 1)]
        for i in range(n + 1):
            costs[i][0] = 0.0
        for i in range(1, n + 1):
            for j in range(1, min(i, m) + 1):
                skip_cost = costs[i - 1][j]
                match_cost = costs[i - 1][j - 1] + abs(
                    longer[i - 1].center_x - shorter[j - 1].center_x
                )
                if match_cost <= skip_cost:
                    costs[i][j] = match_cost
                    matched[i][j] = True
                else:
                    costs[i][j] = skip_cost
        pairs = []
        i, j = n, m
        while j > 0 and i > 0:
            if matched[i][j]:
                long_gene = longer[i - 1]
                short_gene = shorter[j - 1]
                pairs.append(
                    (long_gene, short_gene)
                    if upper_is_longer
                    else (short_gene, long_gene)
                )
                i -= 1
                j -= 1
            else:
                i -= 1
        pairs.reverse()
        return tuple(pairs)

    def _match_pair(self, upper: RowLayout, lower: RowLayout) -> PairLayout:
        upper_by_key: Dict[str, list[GeneLayout]] = {}
        lower_by_key: Dict[str, list[GeneLayout]] = {}
        for layout, target in ((upper, upper_by_key), (lower, lower_by_key)):
            for gene in layout.genes:
                key = gene.gene.family_key
                if not key or gene.gene.evidence_weight < self.minimum_link_evidence:
                    continue
                target.setdefault(key, []).append(gene)

        links = []
        for key in sorted(set(upper_by_key).intersection(lower_by_key)):
            # Pair repeated families in visual order.  This is deterministic,
            # one-to-one, and preserves local collinearity after normalization.
            for upper_gene, lower_gene in self._pair_occurrences(
                upper_by_key[key], lower_by_key[key]
            ):
                links.append(
                    LinkLayout(
                        upper_row=upper.row_index,
                        lower_row=lower.row_index,
                        upper_gene=upper_gene,
                        lower_gene=lower_gene,
                        family_key=key,
                        evidence=upper_gene.gene.evidence,
                        evidence_weight=min(
                            upper_gene.gene.evidence_weight,
                            lower_gene.gene.evidence_weight,
                        ),
                        inverted=upper_gene.display_strand != lower_gene.display_strand,
                    )
                )
        links.sort(key=lambda link: (link.upper_gene.center_x, link.lower_gene.center_x, link.family_key))
        return PairLayout(upper.row_index, lower.row_index, tuple(links))

    def layout(
        self,
        model: ContextRenderModel,
        *,
        density: str = "compact",
        normalize_orientation: bool = True,
        short_rows: str = "inline",
        split_molecules: bool = False,
    ) -> ContextLayout:
        """Place the rows top to bottom.

        short_rows decides what happens to contigs that hold only the focal
        gene: "inline" draws them like any row, "folded" replaces them with
        one summary line at the bottom, "expanded" draws that line with thin
        rows under it, and "omitted" leaves them out.  split_molecules puts a
        labelled divider above each plasmid or chromosome block.
        """
        comfortable = str(density).lower() == "comfortable"
        pitch = self.comfortable_pitch if comfortable else self.compact_pitch
        gene_height = self.comfortable_gene_height if comfortable else self.compact_gene_height
        mode = str(short_rows or "inline").strip().lower()
        if mode not in {"inline", "folded", "expanded", "omitted"}:
            mode = "inline"
        if mode == "inline":
            main_rows: List[ContextRow] = list(model.rows)
            short_list: List[ContextRow] = []
        else:
            main_rows = [row for row in model.rows if not row.is_short]
            short_list = [row for row in model.rows if row.is_short]
        folded_below = short_list if mode in {"folded", "expanded"} else []
        dividers_at = self._molecule_blocks(main_rows, folded_below) if split_molecules else {}

        placed: List[RowLayout] = []
        dividers: List[DividerLayout] = []
        # Top edge of the next band; the first row is centered on top_margin.
        cursor = self.top_margin - pitch / 2.0
        for row in main_rows:
            block = dividers_at.get(len(placed))
            if block is not None:
                dividers.append(DividerLayout(
                    center_y=cursor + self.divider_height / 2.0,
                    height=self.divider_height,
                    molecule=row.molecule,
                    text=block[0],
                    folded=block[1],
                ))
                cursor += self.divider_height
            placed.append(self._layout_row(
                row,
                len(placed),
                pitch=pitch,
                gene_height=gene_height,
                normalize_orientation=normalize_orientation,
                center_y=cursor + pitch / 2.0,
            ))
            cursor += pitch
        fold: Optional[FoldLayout] = None
        if folded_below:
            fold = FoldLayout(
                center_y=cursor + self.fold_height / 2.0,
                height=self.fold_height,
                rows=tuple(folded_below),
                expanded=mode == "expanded",
            )
            cursor += self.fold_height
            if mode == "expanded":
                for row in folded_below:
                    placed.append(self._layout_row(
                        row,
                        len(placed),
                        pitch=self.short_pitch,
                        gene_height=self.short_gene_height,
                        normalize_orientation=normalize_orientation,
                        center_y=cursor + self.short_pitch / 2.0,
                        short=True,
                    ))
                    cursor += self.short_pitch
        rows = tuple(placed)
        # No ribbons run into or between the thin rows of short contigs.
        pairs = tuple(
            self._match_pair(rows[index], rows[index + 1])
            if not (rows[index].short or rows[index + 1].short)
            else PairLayout(rows[index].row_index, rows[index + 1].row_index, ())
            for index in range(max(0, len(rows) - 1))
        )
        raw_min_x = min((row.min_x for row in rows), default=-500.0)
        raw_max_x = max((row.max_x for row in rows), default=500.0)
        min_x = raw_min_x - self.horizontal_margin
        max_x = raw_max_x + self.horizontal_margin
        min_y = 0.0
        max_y = cursor + self.bottom_margin if (rows or fold is not None) else 240.0
        return ContextLayout(
            model=model,
            rows=rows,
            pairs=pairs,
            min_x=min_x,
            max_x=max_x,
            min_y=min_y,
            max_y=max_y,
            px_per_bp=self.px_per_bp,
            dividers=tuple(dividers),
            fold=fold,
        )


def link_payload(link: LinkLayout) -> dict:
    evidence_counts = {link.evidence: 1}
    orientation = "inverted" if link.inverted else "same orientation"
    match = {
        "key": link.family_key,
        "label": link.label,
        "upper_row": link.upper_row,
        "upper_gene": link.upper_gene.visual_index,
        "upper_strand": link.upper_gene.display_strand,
        "lower_row": link.lower_row,
        "lower_gene": link.lower_gene.visual_index,
        "lower_strand": link.lower_gene.display_strand,
        "evidence": link.evidence,
    }
    return {
        "title": "Inverted homology link" if link.inverted else "Homology link",
        "matched_gene_count": 1,
        "orientation": orientation,
        "direction": -1 if link.inverted else 1,
        "evidence_counts": evidence_counts,
        "evidence_summary": link.evidence,
        "labels": [link.label] if link.label else [],
        "label_preview": link.label,
        "matches": [match],
    }
