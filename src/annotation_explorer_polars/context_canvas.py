from __future__ import annotations

import contextlib
import math
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

from PyQt6 import QtCore, QtGui, QtWidgets

from .context_layout import (
    ContextLayout,
    ContextLayoutEngine,
    DividerLayout,
    FoldLayout,
    GeneLayout,
    LinkLayout,
    PairLayout,
    RowLayout,
    row_molecule_counts,
)
from .context_model import (
    CATEGORY_COLORS,
    CATEGORY_LABELS,
    FAMILY_PALETTE,
    ContextRenderModel,
    clean_text,
    truthy_value,
)
from .resistance_mutation import MUTATION_LEGEND_TEXT, paint_mutation_badge
from .theme import COLORS


# Chrome colors come from the app theme; gene and ribbon fills stay with the
# scientific palettes in context_model.
CANVAS_BACKGROUND = QtGui.QColor(COLORS["surface"])
ROW_ALT_BACKGROUND = QtGui.QColor(COLORS["surface_subtle"])
BASELINE_COLOR = QtGui.QColor(COLORS["border"])
SEPARATOR_COLOR = QtGui.QColor(COLORS["border_soft"])
HEADER_BACKGROUND = QtGui.QColor(COLORS["dialog"])
TEXT_COLOR = QtGui.QColor(COLORS["text"])
STRONG_TEXT_COLOR = QtGui.QColor(COLORS["text_strong"])
MUTED_TEXT_COLOR = QtGui.QColor(COLORS["text_muted"])
FOCAL_COLOR = QtGui.QColor(COLORS["text_strong"])
SELECTED_COLOR = QtGui.QColor(COLORS["primary"])
HOVER_COLOR = QtGui.QColor(COLORS["focus"])
HOVER_HALO_COLOR = QtGui.QColor(COLORS["focus_soft"])
ROW_HOVER_COLOR = QtGui.QColor(COLORS["row_hover"])
FOCAL_BAND_COLOR = QtGui.QColor(COLORS["primary_soft"])
REFERENCE_COLOR = QtGui.QColor(COLORS["primary"])
PUBLICATION_BACKGROUND = QtGui.QColor(COLORS["surface"])
# Ribbon fills: grey joins genes running the same way, red joins inverted ones.
LINK_RGB = (100, 116, 139)
INVERTED_LINK_RGB = (213, 94, 94)
# White veil laid over the plot during a spotlight, so everything outside the
# traced family fades slightly.
SPOTLIGHT_VEIL_ALPHA = 105
# Lowest zoom at which the "All" labels mode adds names; at "Fit reference"
# in a 1000 px window the scale is about 0.3.
ALL_LABELS_MIN_SCALE = 0.25
# "Key genes" labels these gene types on every row once zoomed in this far.
# With the Settings panel open, "Fit reference" in a 1600 px window lands
# near 0.43, so the opening view still names them (labels never overlap).
KEY_LABEL_CATEGORIES = frozenset({"AMR", "MOBILE", "VIRULENCE"})
KEY_LABELS_MIN_SCALE = 0.4
# No theme token is dark enough for a printed baseline.
PUBLICATION_BASELINE = QtGui.QColor("#b8c2cc")
# Smallest text the canvas and row rail draw, in device-independent pixels.
MIN_TEXT_PX = 10

RAIL_HEADER_HEIGHT = 24
RAIL_BADGE_X = 10.0
RAIL_MATCH_WIDTH = 56.0
RAIL_RIGHT_MARGIN = 8.0
RAIL_WIDTH = 280
MATCH_TOOLTIP = "Match to the reference neighborhood"
REFERENCE_TOOLTIP = "Reference neighborhood: every other row is scored against this one"
SHORT_TOOLTIP = "Contig too short to compare"
PARTIAL_TOOLTIP = (
    "Partial: the contig ends within the shown flank, so part of the neighborhood is missing"
)
PARTIAL_REFERENCE_TOOLTIP = (
    "Partial: its contig ends within the shown flank, so every match covers only the side it shows"
)
MIRRORED_TOOLTIP = "Shown mirrored, so its genes run the same way as the reference"
# Short-contig summary line, in the style of the annotation table's
# "N genes not shown" gap rows: light band, muted italic text.
GAP_BACKGROUND = QtGui.QColor(COLORS["dialog"])
GAP_TEXT_COLOR = QtGui.QColor(COLORS["text_muted"])
LINK_TEXT_COLOR = QtGui.QColor(COLORS["primary"])
CONTIG_END_COLOR = QtGui.QColor(COLORS["text_muted"])
# Rail second line: a small colored dot and a word name the molecule type.
MOLECULE_DOTS = {
    "plasmid": ("Plasmid", COLORS["success"]),
    "chromosome": ("Chromosome", COLORS["primary"]),
    "mixed": ("Mixed", COLORS["text_disabled"]),
    "unknown": ("Unknown", COLORS["text_disabled"]),
}
RAIL_DOT_SIZE = 7.0
RAIL_DETAIL_PX = 11
# 'partial' is a hint, not an alarm: muted text with a little warmth.
PARTIAL_TEXT_COLOR = QtGui.QColor("#846c4b")


def gene_arrow_path(rect: QtCore.QRectF, strand: str, scale_x: float) -> QtGui.QPainterPath:
    """Return a proportional five-point gene arrow with a screen-stable head
    (the neighborhood plot's and the genome track's gene glyph)."""
    scale_x = max(abs(float(scale_x)), 1e-6)
    minimum_width = 2.5 / scale_x
    if rect.width() < minimum_width:
        center = rect.center().x()
        rect = QtCore.QRectF(center - minimum_width / 2.0, rect.top(), minimum_width, rect.height())
    head = min(rect.width() * 0.42, 8.0 / scale_x)
    head = max(min(rect.width() * 0.42, head), min(rect.width() * 0.28, 2.0 / scale_x))
    path = QtGui.QPainterPath()
    if strand == "-":
        shoulder = rect.left() + head
        path.moveTo(rect.right(), rect.top())
        path.lineTo(shoulder, rect.top())
        path.lineTo(rect.left(), rect.center().y())
        path.lineTo(shoulder, rect.bottom())
        path.lineTo(rect.right(), rect.bottom())
    else:
        shoulder = rect.right() - head
        path.moveTo(rect.left(), rect.top())
        path.lineTo(shoulder, rect.top())
        path.lineTo(rect.right(), rect.center().y())
        path.lineTo(shoulder, rect.bottom())
        path.lineTo(rect.left(), rect.bottom())
    path.closeSubpath()
    return path


def _nice_scale_bp(target_bp: float) -> float:
    target_bp = max(1.0, float(target_bp))
    magnitude = 10.0 ** math.floor(math.log10(target_bp))
    normalized = target_bp / magnitude
    nice = 1.0 if normalized < 1.5 else 2.0 if normalized < 3.5 else 5.0 if normalized < 7.5 else 10.0
    return nice * magnitude


def _primary_row_label(row: Any) -> str:
    """Use a collapsed group's label instead of its arbitrary representative."""
    if row.sample_count > 1:
        return row.title or f"{row.sample_count} genomes"
    return row.sample or row.title


def _secondary_row_label(row_layout: Any) -> str:
    """Second rail line: 'contig 48', or a merged group's summary.

    A merged row stands for several contigs, so one member's contig would be
    misleading; it names the reference genome instead when the group holds
    it, otherwise how many contigs were merged.  Mirroring and 'partial'
    follow it on the rail as words of their own.
    """
    row = row_layout.row
    if row.sample_count > 1:
        if row.is_reference and row.sample:
            return f"includes {row.sample}"
        try:
            contigs = int(row.source.get("collapsed_track_count") or 0)
        except (TypeError, ValueError):
            contigs = 0
        return f"{contigs} contigs" if contigs > 1 else ""
    return _row_place_text(row)


def _molecule_word(molecule: str) -> str:
    return MOLECULE_DOTS.get(molecule, MOLECULE_DOTS["unknown"])[0]


def _molecule_dot_color(molecule: str) -> QtGui.QColor:
    return QtGui.QColor(MOLECULE_DOTS.get(molecule, MOLECULE_DOTS["unknown"])[1])


def _molecule_tooltip(molecule: str, counts: Optional[Dict[str, int]] = None) -> str:
    if molecule == "plasmid":
        return "Plasmid contig (MOB-suite prediction)"
    if molecule == "chromosome":
        return "Chromosome contig (MOB-suite prediction)"
    if molecule == "mixed":
        names = {"plasmid": "plasmid", "chromosome": "chromosome", "unknown": "of unknown type"}
        parts = [
            f"{counts[kind]} {names[kind]}"
            for kind in ("plasmid", "chromosome", "unknown")
            if counts and counts.get(kind)
        ]
        if parts:
            return "Merged row of contigs of more than one type: " + ", ".join(parts)
        return "Merged row with both plasmid and chromosome contigs"
    return "Molecule type unknown: no plasmid prediction for this contig"


def _strand_words(strand: Any) -> str:
    text = str(strand or "").strip()
    if text == "-":
        return "minus strand"
    if text == "+":
        return "plus strand"
    return "unknown strand"


def _is_plasmid_context_row(row: Any) -> bool:
    source = row.source if isinstance(row.source, dict) else {}
    plasmid = str(source.get("is_plasmid")).strip().lower() in {"1", "true", "yes"}
    return plasmid and bool(clean_text(source.get("cluster_key")))


def _plasmid_name_from_title(title: str, sample: str) -> str:
    """The plasmid part of 'sample, plasmid pX1' (or the older 'sample | pX1')."""
    title = clean_text(title)
    for separator in ("|", ","):
        head, found, tail = title.partition(separator)
        if found and (not sample or head.strip() == sample):
            title = tail.strip()
            break
    return title


def _row_place_text(row: Any, *, with_locus: bool = True) -> str:
    """'contig 48' for a contig row, the plasmid name for a plasmid-context row.

    Plasmid-context rows (one per plasmid, with a 'cluster_key') carry an
    instance key in 'contig', which is not a contig name, so their display
    title is used instead.  A gene row on a plasmid contig keeps its contig.
    """
    if _is_plasmid_context_row(row):
        title = _plasmid_name_from_title(row.title, row.sample)
        if title and not title.lower().startswith("plasmid"):
            return f"plasmid {title}"
        return title
    contig = row.contig
    if contig and not contig.lower().startswith("contig"):
        contig = f"contig {contig}"
    if with_locus and contig and truthy_value(row.source.get("repeated_contig")):
        # Two copies of the focal gene on one contig: say where each sits.
        contig = f"{contig} at {_locus_text(row)}"
    return contig


def _locus_text(row: Any) -> str:
    """Where the focal gene starts on its contig, as '1,234 kb'."""
    focal = row.focal_gene
    position = float(focal.low) if focal is not None else float(row.anchor_bp)
    kb = max(0.0, position) / 1000.0
    return f"{kb:,.0f} kb" if kb >= 10 else f"{kb:.1f} kb"


def _mark_repeated_contigs(rows: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Flag rows that share a sample and contig with another row.

    Such rows are second copies of the focal gene on the same contig; their
    rail text adds the locus position so they never read as duplicates.
    Merged rows and plasmid-context rows are left alone.
    """
    def key(row: Dict[str, Any]) -> Optional[Tuple[str, str]]:
        if not isinstance(row, dict) or row.get("collapsed_samples"):
            return None
        if clean_text(row.get("cluster_key")) and truthy_value(row.get("is_plasmid")):
            return None
        sample, contig = clean_text(row.get("sample")), clean_text(row.get("contig"))
        return (sample, contig) if sample and contig else None

    counts: Dict[Tuple[str, str], int] = {}
    for row in rows:
        found = key(row)
        if found is not None:
            counts[found] = counts.get(found, 0) + 1
    marked: List[Dict[str, Any]] = []
    for row in rows:
        found = key(row)
        if found is not None and counts[found] > 1:
            row = {**row, "repeated_contig": True}
        marked.append(row)
    return marked


def _row_sentence(row_layout: Any, focal_label: str = "") -> str:
    """One plain sentence about a row, for the tooltip on its name.

    For example 'Plasmid contig 48 of GCA_006366115.1. 67% match to the
    reference. Shown mirrored.'
    """
    row = row_layout.row
    parts: List[str] = []
    if row.sample_count > 1:
        parts.append(f"{row.sample_count} genomes with identical neighborhoods.")
    else:
        place = _row_place_text(row)
        kind = {"plasmid": "Plasmid", "chromosome": "Chromosome"}.get(row.molecule, "")
        if place.lower().startswith("plasmid"):
            kind = ""
        where = f"{kind} {place}".strip() if kind else (place[:1].upper() + place[1:])
        if where and row.sample:
            where = f"{where} of {row.sample}"
        elif row.sample:
            where = row.sample
        if row.molecule == "unknown":
            where += " (molecule type unknown)"
        if where:
            parts.append(where + ".")
    if row.is_reference:
        parts.append("Reference: every other row is scored against this one.")
    elif row.is_short:
        gene = f" {focal_label}" if focal_label else " the focal gene"
        parts.append(f"Contig too short to compare: it carries only{gene}.")
    elif row.similarity is not None:
        parts.append(f"{row.similarity:.0%} match to the reference.")
    if row.partial and not row.is_short:
        parts.append(
            "The contig ends within the shown flank, so part of the neighborhood is missing"
            + (": every match covers only the side it shows." if row.is_reference else ".")
        )
    if row_layout.flipped:
        parts.append("Shown mirrored.")
    sentence = " ".join(parts)
    # The row's own metadata follows; its 'Sample:' and 'Contig:' lines
    # would only repeat the sentence.
    meta = [
        line
        for line in clean_text(row.source.get("meta_tooltip")).splitlines()
        if not line.startswith(("Sample: ", "Contig: "))
    ]
    while meta and not meta[0].strip():
        meta.pop(0)
    return "\n\n".join(part for part in (sentence, "\n".join(meta)) if part)


def _gene_position_text(gene: Any, row: Any, *, with_strand: bool = True) -> str:
    """Where a gene sits, in words that never pretend to be contig coordinates.

    Plasmid-context rows join a plasmid's contigs end to end, so their
    drawing positions are offsets into the plasmid; the gene's own contig
    coordinates are used when the annotation has them.
    """
    strand = f", {_strand_words(gene.gene.strand)}" if with_strand else ""
    # The gene's own coordinates follow, so no locus hint is needed.
    place = _row_place_text(row, with_locus=False)
    if _is_plasmid_context_row(row) or bool(row.source.get("joined_contigs")):
        payload = gene.gene.payload
        contig = clean_text(payload.get("Contig"))
        try:
            start = float(clean_text(payload.get("StartInt")) or clean_text(payload.get("Start")))
            stop = float(clean_text(payload.get("StopInt")) or clean_text(payload.get("Stop")))
        except (TypeError, ValueError):
            start = stop = float("nan")
        prefix = f"{place}, " if place else ""
        if contig and math.isfinite(start) and math.isfinite(stop):
            low, high = sorted((int(start), int(stop)))
            contig = contig if contig.lower().startswith("contig") else f"contig {contig}"
            return f"{prefix}{contig}, {low:,} to {high:,}{strand}"
        low, high = int(gene.gene.low), int(gene.gene.high)
        return f"{prefix}position in plasmid {low:,} to {high:,}{strand}"
    position = f"{int(gene.gene.low):,} to {int(gene.gene.high):,}{strand}"
    return f"{place}, {position}" if place else position


def _fold_sentence(fold: Any, focal_label: str) -> str:
    """'46 contigs carry aph(3')-Ia with no neighbors (40 chromosome, 6 plasmid)'."""
    count = fold.contig_count
    gene = focal_label or "the focal gene"
    verb = "carries" if count == 1 else "carry"
    noun = "contig" if count == 1 else "contigs"
    counts = fold.molecule_counts()
    parts = [
        f"{counts[kind]} {kind}"
        for kind in ("chromosome", "plasmid", "mixed", "unknown")
        if counts.get(kind)
    ]
    detail = f" ({', '.join(parts)})" if parts else ""
    return f"{count} {noun} {verb} {gene} with no neighbors{detail}"


def _scale_label(scale_bp: float) -> str:
    return f"{scale_bp / 1000.0:g} kb" if scale_bp >= 1000 else f"{scale_bp:g} bp"


class _GenomeView(QtWidgets.QGraphicsView):
    def __init__(self, owner: "ContextCanvas") -> None:
        super().__init__(owner)
        self.owner = owner
        self.setRenderHints(
            QtGui.QPainter.RenderHint.Antialiasing
            | QtGui.QPainter.RenderHint.TextAntialiasing
            | QtGui.QPainter.RenderHint.SmoothPixmapTransform
        )
        self.setBackgroundBrush(QtGui.QBrush(CANVAS_BACKGROUND))
        self.setDragMode(QtWidgets.QGraphicsView.DragMode.ScrollHandDrag)
        self.setTransformationAnchor(QtWidgets.QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        # Resizing keeps the left edge in place: opening the details panel
        # narrows the plot without sliding the focal band sideways.
        self.setResizeAnchor(QtWidgets.QGraphicsView.ViewportAnchor.NoAnchor)
        self.setViewportUpdateMode(QtWidgets.QGraphicsView.ViewportUpdateMode.BoundingRectViewportUpdate)
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        # Rows start under the header strip even when they do not fill the
        # viewport (for example with merged neighborhoods).
        self.setAlignment(QtCore.Qt.AlignmentFlag.AlignHCenter | QtCore.Qt.AlignmentFlag.AlignTop)
        self.setMouseTracking(True)
        self.viewport().setMouseTracking(True)
        self.setFocusPolicy(QtCore.Qt.FocusPolicy.StrongFocus)

    def mouseMoveEvent(self, event: QtGui.QMouseEvent) -> None:
        self.owner._handle_mouse_move(event)
        super().mouseMoveEvent(event)

    def mousePressEvent(self, event: QtGui.QMouseEvent) -> None:
        if self.owner._handle_mouse_press(event):
            event.accept()
            return
        super().mousePressEvent(event)

    def contextMenuEvent(self, event: QtGui.QContextMenuEvent) -> None:
        layout = self.owner.layout_result
        if layout is not None:
            row = layout.row_at(self.mapToScene(event.pos()).y())
            if row is not None:
                self.owner._request_row_menu(row, event.globalPos())
                event.accept()
                return
        super().contextMenuEvent(event)

    def mouseDoubleClickEvent(self, event: QtGui.QMouseEvent) -> None:
        if self.owner._handle_mouse_double_click(event):
            event.accept()
            return
        super().mouseDoubleClickEvent(event)

    def wheelEvent(self, event: QtGui.QWheelEvent) -> None:
        # Plain wheel remains useful for long comparisons.  Ctrl+wheel zooms
        # horizontally around the pointer without changing row height.
        if event.modifiers() & QtCore.Qt.KeyboardModifier.ControlModifier:
            self.owner.zoom_by(1.16 if event.angleDelta().y() > 0 else 1.0 / 1.16)
            event.accept()
            return
        super().wheelEvent(event)

    def keyPressEvent(self, event: QtGui.QKeyEvent) -> None:
        key = event.key()
        if key == QtCore.Qt.Key.Key_Escape:
            # The first Esc clears the selection; with nothing selected the
            # event goes up to the dialog, which closes on the next Esc.
            if self.owner.has_selection():
                self.owner.clear_selection()
                event.accept()
                return
            event.ignore()
            return
        if key == QtCore.Qt.Key.Key_F:
            self.owner.reset_zoom()
            event.accept()
            return
        if key in (QtCore.Qt.Key.Key_Plus, QtCore.Qt.Key.Key_Equal):
            self.owner.zoom_by(1.18)
            event.accept()
            return
        if key == QtCore.Qt.Key.Key_Minus:
            self.owner.zoom_by(1.0 / 1.18)
            event.accept()
            return
        directions = {
            QtCore.Qt.Key.Key_Left: (-1, 0),
            QtCore.Qt.Key.Key_Right: (1, 0),
            QtCore.Qt.Key.Key_Up: (0, -1),
            QtCore.Qt.Key.Key_Down: (0, 1),
        }
        if key in directions and self.owner.navigate_selection(*directions[key]):
            event.accept()
            return
        if key in (QtCore.Qt.Key.Key_Return, QtCore.Qt.Key.Key_Enter):
            if self.owner.open_selected_gene():
                event.accept()
                return
        super().keyPressEvent(event)

    def leaveEvent(self, event: QtCore.QEvent) -> None:
        self.owner._pointer = None
        self.owner._set_hover(None, None)
        super().leaveEvent(event)

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        super().resizeEvent(event)
        # A wider plot keeps its left edge: whatever was in view still is.
        if event.oldSize().width() <= 0 or event.size().width() <= event.oldSize().width():
            self.owner._keep_focal_column_visible()
        # Opening the details panel narrows the plot; keep the clicked gene
        # on screen instead of letting it slide out at the edge.
        self.owner.ensure_selection_visible()
        self.owner._viewport_changed()

    def drawForeground(self, painter: QtGui.QPainter, rect: QtCore.QRectF) -> None:
        super().drawForeground(painter, rect)
        self.owner._draw_view_overlay(painter)


class _ColumnHeader(QtWidgets.QWidget):
    """Thin caption strip that sits above the row rail or the plot."""

    def __init__(self, owner: "ContextCanvas", *, rail: bool) -> None:
        super().__init__(owner)
        self.owner = owner
        self._rail = bool(rail)
        self.caption = ""
        # A small widget (the dialog's 'Pinned' chip) kept at the right end
        # of the strip, or beside the focal tab when the tab sits there.
        self._corner: Optional[QtWidgets.QWidget] = None
        # Set when the strip is too short for both: the tab's name then
        # yields, as the focal band and the title still name the gene.
        self._corner_hides_tab = False
        self.setFixedHeight(RAIL_HEADER_HEIGHT)
        self.setMouseTracking(True)

    def set_caption(self, text: str) -> None:
        self.caption = str(text or "")
        self.update()

    def set_corner_widget(self, widget: QtWidgets.QWidget) -> None:
        self._corner = widget
        widget.setParent(self)
        widget.installEventFilter(self)
        self.refresh()

    def refresh(self) -> None:
        """Repaint after the plot scrolled, zoomed or was rebuilt."""
        self.place_corner()
        self.update()

    def eventFilter(self, watched: QtCore.QObject, event: QtCore.QEvent) -> bool:
        if watched is self._corner and event.type() in (
            QtCore.QEvent.Type.Show,
            QtCore.QEvent.Type.Hide,
            QtCore.QEvent.Type.LayoutRequest,
        ):
            # New text or shown: place it again and keep the caption clear of it.
            self.refresh()
        return super().eventFilter(watched, event)

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        super().resizeEvent(event)
        self.place_corner()

    @staticmethod
    def _tab_font(base: QtGui.QFont) -> QtGui.QFont:
        font = QtGui.QFont(base)
        font.setPixelSize(11)
        font.setWeight(QtGui.QFont.Weight.DemiBold)
        return font

    def place_corner(self) -> None:
        """Right end of the strip, else just left or right of the focal tab:
        never on top of the tab."""
        widget = self._corner
        self._corner_hides_tab = False
        if widget is None or widget.isHidden():
            return
        hint = widget.sizeHint()
        right = float(self.width() - RAIL_RIGHT_MARGIN)
        # Never wider than the strip, so its left part is never cut off.
        width = max(1, min(hint.width(), int(right - RAIL_BADGE_X)))
        height = max(1, min(hint.height(), self.height() - 4))
        tab = self._focal_tab_rect(self._tab_font(self.font()))
        candidates = [right - width]
        if tab is not None:
            candidates += [tab.left() - 10.0 - width, tab.right() + 10.0]

        def clear_of_tab(x: float) -> bool:
            return tab is None or x + width + 8.0 <= tab.left() or x >= tab.right() + 8.0

        x = next(
            (c for c in candidates if c >= RAIL_BADGE_X and c + width <= right + 0.5 and clear_of_tab(c)),
            None,
        )
        if x is None:
            x = max(RAIL_BADGE_X, candidates[0])
            self._corner_hides_tab = True
        widget.setGeometry(int(round(x)), (self.height() - height) // 2, width, height)

    def _focal_tab_rect(self, font: QtGui.QFont) -> Optional[QtCore.QRectF]:
        """Where the focal gene's name sits, centered over the focal band."""
        band = self.owner._focal_band_rect
        label = self.owner.focal_label
        if self._rail or band is None or not label:
            return None
        viewport = self.owner.graphics_view.viewport()
        view_point = self.owner.graphics_view.mapFromScene(QtCore.QPointF(band.center().x(), 0.0))
        x = float(self.mapFromGlobal(viewport.mapToGlobal(view_point)).x())
        width = float(QtGui.QFontMetrics(font).horizontalAdvance(label)) + 16.0
        rect = QtCore.QRectF(x - width / 2.0, 4.0, width, float(self.height()) - 4.0)
        if rect.right() < 0.0 or rect.left() > float(self.width()):
            return None
        return rect

    def tooltip_at(self, x: float) -> str:
        if not self._rail:
            tab = None if self._corner_hides_tab else self._focal_tab_rect(self._tab_font(self.font()))
            if tab is not None and tab.left() <= x <= tab.right():
                return f"Focal gene: every row is lined up on {self.owner.focal_label}"
            return self.caption
        if x >= self.width() - RAIL_MATCH_WIDTH - RAIL_RIGHT_MARGIN:
            return MATCH_TOOLTIP
        return "Genome accession, molecule type and contig"

    def event(self, event: QtCore.QEvent) -> bool:
        if event.type() == QtCore.QEvent.Type.ToolTip:
            text = self.tooltip_at(float(event.pos().x()))
            if text:
                QtWidgets.QToolTip.showText(event.globalPos(), text, self)
            else:
                QtWidgets.QToolTip.hideText()
            return True
        return super().event(event)

    def paintEvent(self, event: QtGui.QPaintEvent) -> None:
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.TextAntialiasing, True)
        painter.fillRect(self.rect(), HEADER_BACKGROUND)
        font = painter.font()
        font.setPixelSize(11)
        font.setWeight(QtGui.QFont.Weight.DemiBold)
        painter.setFont(font)
        painter.setPen(MUTED_TEXT_COLOR)
        text_rect = QtCore.QRectF(self.rect()).adjusted(RAIL_BADGE_X, 0.0, -RAIL_RIGHT_MARGIN, 0.0)
        align_v = QtCore.Qt.AlignmentFlag.AlignVCenter
        if self._rail:
            painter.drawText(text_rect, QtCore.Qt.AlignmentFlag.AlignLeft | align_v, "Genome")
            painter.drawText(text_rect, QtCore.Qt.AlignmentFlag.AlignRight | align_v, "Match")
        else:
            # The focal band's name sits here, on top of the band, so it
            # stays visible however far the rows are scrolled.
            tab = None if self._corner_hides_tab else self._focal_tab_rect(font)
            corner = self._corner
            if corner is not None and not corner.isHidden() and corner.x() > text_rect.left():
                text_rect.setRight(min(text_rect.right(), float(corner.x()) - 10.0))
            if tab is not None:
                text_rect.setRight(min(text_rect.right(), tab.left() - 10.0))
                painter.setPen(QtCore.Qt.PenStyle.NoPen)
                painter.setBrush(FOCAL_BAND_COLOR)
                painter.drawRoundedRect(tab.adjusted(0.0, 0.0, 0.0, 6.0), 5.0, 5.0)
                painter.setPen(REFERENCE_COLOR)
                painter.drawText(tab, QtCore.Qt.AlignmentFlag.AlignCenter, self.owner.focal_label)
            if self.caption and text_rect.width() > 40:
                font.setWeight(QtGui.QFont.Weight.Normal)
                painter.setFont(font)
                painter.setPen(MUTED_TEXT_COLOR)
                caption = QtGui.QFontMetrics(font).elidedText(
                    self.caption,
                    QtCore.Qt.TextElideMode.ElideRight,
                    max(20, int(text_rect.width())),
                )
                painter.drawText(text_rect, QtCore.Qt.AlignmentFlag.AlignLeft | align_v, caption)
        painter.setPen(QtGui.QPen(BASELINE_COLOR, 1.0))
        painter.drawLine(0, self.height() - 1, self.width(), self.height() - 1)
        if self._rail:
            painter.drawLine(self.width() - 1, 0, self.width() - 1, self.height())


def _rail_font(pixel_size: int, *, bold: bool = False, italic: bool = False) -> QtGui.QFont:
    font = QtGui.QFont()
    font.setPixelSize(pixel_size)
    if bold:
        font.setWeight(QtGui.QFont.Weight.DemiBold)
    font.setItalic(italic)
    return font


def _rail_match_text(row_layout: RowLayout) -> Tuple[str, QtGui.QFont, QtGui.QColor]:
    """The right-hand text of a rail row: 'Reference', '91%' or '-'."""
    row = row_layout.row
    if row.is_reference and not row_layout.short:
        return "Reference", _rail_font(12, bold=True), REFERENCE_COLOR
    if row_layout.short or row.similarity is None:
        return "-", _rail_font(12), MUTED_TEXT_COLOR
    return f"{row.similarity:.0%}", _rail_font(12), MUTED_TEXT_COLOR


def _rail_detail_segments(row_layout: RowLayout, show_molecule: bool) -> List[Tuple[str, str, QtGui.QColor]]:
    """The muted words after the dot: (key, text, color), in reading order.

    For example 'Plasmid', ' · contig 43', ' · mirrored', ' · partial'.
    A thin short-contig row keeps only its contig.
    """
    row = row_layout.row
    segments: List[Tuple[str, str, QtGui.QColor]] = []
    if show_molecule and not row_layout.short:
        segments.append(("molecule", _molecule_word(row.molecule), MUTED_TEXT_COLOR))
    place = _secondary_row_label(row_layout)
    if place:
        segments.append(("place", place, MUTED_TEXT_COLOR))
    if not row_layout.short:
        if row_layout.flipped:
            segments.append(("mirrored", "mirrored", MUTED_TEXT_COLOR))
        if row.partial:
            segments.append(("partial", "partial", PARTIAL_TEXT_COLOR))
    return [
        (key, text if index == 0 else f" \u00b7 {text}", color)
        for index, (key, text, color) in enumerate(segments)
    ]


def _rail_row_layout(
    row_layout: RowLayout, center: float, width: float, show_molecule: bool
) -> Tuple[Dict[str, QtCore.QRectF], List[Tuple[str, str, QtCore.QRectF, QtGui.QColor]]]:
    """Where the parts of one rail row sit; shared by painting and tooltips.

    Returns the named rects ('name', 'match', 'dot', and one per detail
    word) and the detail words as drawn, elided from the right to fit.
    """
    row = row_layout.row
    # Plasmid-context rows are all plasmids, and their place text says so.
    show_molecule = show_molecule and not _is_plasmid_context_row(row)
    right = width - RAIL_RIGHT_MARGIN
    parts: Dict[str, QtCore.QRectF] = {}
    match_text, match_font, _color = _rail_match_text(row_layout)
    match_width = float(QtGui.QFontMetricsF(match_font).horizontalAdvance(match_text))
    # Never narrower than a short percentage, so the tooltip is easy to find.
    match_width = max(match_width, 28.0)
    if row_layout.short:
        # Thin row: accession, dot and contig on one line.
        line_top, line_height = center - 8.0, 16.0
        name_font = _rail_font(11, bold=True)
        detail_top = line_top
    else:
        line_top, line_height = center - 17.0, 16.0
        name_font = _rail_font(12, bold=True)
        detail_top = center + 1.0
    parts["match"] = QtCore.QRectF(right - match_width, line_top, match_width, line_height)
    name_right = parts["match"].left() - 8.0
    name_width = float(QtGui.QFontMetricsF(name_font).horizontalAdvance(_primary_row_label(row)))
    parts["name"] = QtCore.QRectF(RAIL_BADGE_X, line_top, max(20.0, min(name_width, name_right - RAIL_BADGE_X)), line_height)

    x = RAIL_BADGE_X
    detail_right = right
    if row_layout.short:
        x = parts["name"].right() + 8.0
        detail_right = name_right
    if show_molecule and x + RAIL_DOT_SIZE < detail_right:
        parts["dot"] = QtCore.QRectF(
            x, detail_top + (line_height - RAIL_DOT_SIZE) / 2.0, RAIL_DOT_SIZE, RAIL_DOT_SIZE
        )
        x += RAIL_DOT_SIZE + 5.0
    metrics = QtGui.QFontMetricsF(_rail_font(RAIL_DETAIL_PX))
    words: List[Tuple[str, str, QtCore.QRectF, QtGui.QColor]] = []
    for key, text, color in _rail_detail_segments(row_layout, show_molecule):
        room = detail_right - x
        if room < 12.0:
            break
        text_width = float(metrics.horizontalAdvance(text))
        elided = text_width > room
        if elided:
            text = metrics.elidedText(text, QtCore.Qt.TextElideMode.ElideRight, room)
            text_width = float(metrics.horizontalAdvance(text))
        rect = QtCore.QRectF(x, detail_top, text_width, line_height)
        parts[key] = rect
        words.append((key, text, rect, color))
        x = rect.right()
        if elided:
            break
    return parts, words


def _rail_row_geometry(row_layout: RowLayout, center: float, width: float, show_molecule: bool) -> Dict[str, QtCore.QRectF]:
    """The named rects of one rail row (see _rail_row_layout)."""
    return _rail_row_layout(row_layout, center, width, show_molecule)[0]


def _paint_rail_row(
    painter: QtGui.QPainter,
    row_layout: RowLayout,
    center: float,
    width: float,
    *,
    show_molecule: bool = True,
) -> None:
    """Draw one row's two quiet lines; the rail and figures share it.

    Line one: the accession (or 'N genomes') and the match on the right.
    Line two: a molecule dot, then muted words such as 'Plasmid · contig 43'.
    """
    row = row_layout.row
    parts, words = _rail_row_layout(row_layout, center, width, show_molecule)
    painter.save()
    painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
    painter.setRenderHint(QtGui.QPainter.RenderHint.TextAntialiasing, True)
    align_left = QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignVCenter
    name_font = _rail_font(11 if row_layout.short else 12, bold=True)
    name_rect = parts["name"]
    name = QtGui.QFontMetricsF(name_font).elidedText(
        _primary_row_label(row), QtCore.Qt.TextElideMode.ElideRight, name_rect.width() + 1.0
    )
    painter.setFont(name_font)
    painter.setPen(TEXT_COLOR)
    painter.drawText(name_rect.adjusted(0.0, 0.0, 4.0, 0.0), align_left, name)

    match_text, match_font, match_color = _rail_match_text(row_layout)
    painter.setFont(match_font)
    painter.setPen(match_color)
    painter.drawText(
        parts["match"],
        QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter,
        match_text,
    )

    if "dot" in parts:
        painter.setPen(QtCore.Qt.PenStyle.NoPen)
        painter.setBrush(_molecule_dot_color(row.molecule))
        painter.drawEllipse(parts["dot"])
    painter.setFont(_rail_font(RAIL_DETAIL_PX))
    for _key, text, rect, color in words:
        painter.setPen(color)
        painter.drawText(rect.adjusted(0.0, 0.0, 4.0, 0.0), align_left, text)
    painter.restore()


def _paint_band_label(
    painter: QtGui.QPainter,
    rect: QtCore.QRectF,
    text: str,
    *,
    link: str = "",
    italic: bool = True,
) -> QtCore.QRectF:
    """Muted (italic) text with an optional link word; returns the link rect."""
    painter.save()
    font = _rail_font(11, italic=italic)
    if not italic:
        font.setWeight(QtGui.QFont.Weight.DemiBold)
    metrics = QtGui.QFontMetrics(font)
    link_font = _rail_font(11, bold=True)
    link_width = float(QtGui.QFontMetrics(link_font).horizontalAdvance(link)) if link else 0.0
    text_width = max(20, int(rect.width() - (link_width + 10.0 if link else 0.0)))
    shown = metrics.elidedText(text, QtCore.Qt.TextElideMode.ElideRight, text_width)
    painter.setFont(font)
    painter.setPen(GAP_TEXT_COLOR)
    align = QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignVCenter
    painter.drawText(rect, align, shown)
    link_rect = QtCore.QRectF()
    if link:
        left = rect.left() + float(metrics.horizontalAdvance(shown)) + 10.0
        link_rect = QtCore.QRectF(left, rect.top(), link_width, rect.height())
        painter.setFont(link_font)
        painter.setPen(LINK_TEXT_COLOR)
        painter.drawText(link_rect, align, link)
    painter.restore()
    return link_rect


class _RowHeaderRail(QtWidgets.QWidget):
    def __init__(self, owner: "ContextCanvas") -> None:
        super().__init__(owner)
        self.owner = owner
        self.setObjectName("ContextRowHeader")
        self.setMouseTracking(True)
        self.setAutoFillBackground(True)
        palette = self.palette()
        palette.setColor(QtGui.QPalette.ColorRole.Window, CANVAS_BACKGROUND)
        self.setPalette(palette)

    def sizeHint(self) -> QtCore.QSize:
        return QtCore.QSize(RAIL_WIDTH, 480)

    def _scene_y(self, y: float) -> float:
        return self.owner.graphics_view.mapToScene(QtCore.QPoint(0, int(y))).y()

    def _row_at_widget_y(self, y: float) -> Optional[RowLayout]:
        layout = self.owner.layout_result
        if layout is None:
            return None
        return layout.row_at(self._scene_y(y))

    def _widget_center(self, scene_y: float) -> float:
        return float(self.owner.graphics_view.mapFromScene(QtCore.QPointF(0.0, scene_y)).y())

    def tooltip_at(self, x: float, y: float) -> str:
        """Explain the match, the molecule dot, a hint word or the row under the pointer."""
        owner = self.owner
        if owner._fold_at_scene_y(self._scene_y(y)):
            return owner.fold_tooltip()
        divider = owner._divider_at_scene_y(self._scene_y(y))
        if divider is not None:
            return owner.divider_tooltip(divider)
        row = self._row_at_widget_y(y)
        if row is None:
            return ""
        parts = _rail_row_geometry(row, self._widget_center(row.center_y), float(self.width()), owner.show_molecule_badges)
        point = QtCore.QPointF(x, y)

        def hit(key: str) -> bool:
            return key in parts and parts[key].adjusted(-3.0, -3.0, 3.0, 3.0).contains(point)

        if hit("match"):
            if row.row.is_reference and not row.short:
                return f"{REFERENCE_TOOLTIP}. {PARTIAL_REFERENCE_TOOLTIP}" if row.row.partial else REFERENCE_TOOLTIP
            if row.short:
                return SHORT_TOOLTIP
            if row.row.partial:
                return f"{MATCH_TOOLTIP}. {PARTIAL_TOOLTIP}"
            if row.row.similarity is not None:
                return MATCH_TOOLTIP
        if hit("dot") or hit("molecule"):
            return _molecule_tooltip(row.row.molecule, row_molecule_counts(row.row))
        if hit("mirrored"):
            return MIRRORED_TOOLTIP
        if hit("partial"):
            return PARTIAL_REFERENCE_TOOLTIP if row.row.is_reference else PARTIAL_TOOLTIP
        return _row_sentence(row, owner.focal_label)

    def event(self, event: QtCore.QEvent) -> bool:
        if event.type() == QtCore.QEvent.Type.ToolTip:
            pos = event.pos()
            text = self.tooltip_at(float(pos.x()), float(pos.y()))
            if text:
                QtWidgets.QToolTip.showText(event.globalPos(), text, self)
            else:
                QtWidgets.QToolTip.hideText()
            return True
        return super().event(event)

    def mouseMoveEvent(self, event: QtGui.QMouseEvent) -> None:
        self.owner._pointer = ("rail", QtCore.QPointF(event.position()))
        row = self._row_at_widget_y(event.position().y())
        self.owner._set_hover(None, row)
        on_fold = self.owner._fold_at_scene_y(self._scene_y(event.position().y()))
        self.setCursor(
            QtCore.Qt.CursorShape.PointingHandCursor if on_fold else QtCore.Qt.CursorShape.ArrowCursor
        )
        super().mouseMoveEvent(event)

    def mousePressEvent(self, event: QtGui.QMouseEvent) -> None:
        if (
            event.button() == QtCore.Qt.MouseButton.LeftButton
            and self.owner._fold_at_scene_y(self._scene_y(event.position().y()))
        ):
            self.owner.toggle_short_rows()
            event.accept()
            return
        super().mousePressEvent(event)

    def contextMenuEvent(self, event: QtGui.QContextMenuEvent) -> None:
        row = self._row_at_widget_y(float(event.pos().y()))
        if row is not None:
            self.owner._request_row_menu(row, event.globalPos())
            event.accept()
            return
        super().contextMenuEvent(event)

    def leaveEvent(self, event: QtCore.QEvent) -> None:
        self.owner._pointer = None
        self.owner._set_hover(None, None)
        super().leaveEvent(event)

    def paintEvent(self, event: QtGui.QPaintEvent) -> None:
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.TextAntialiasing, True)
        painter.fillRect(self.rect(), CANVAS_BACKGROUND)
        layout = self.owner.layout_result
        if layout is None or (not layout.rows and layout.fold is None):
            painter.setPen(MUTED_TEXT_COLOR)
            painter.drawText(
                self.rect().adjusted(14, 16, -10, -10),
                QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignTop,
                "No rows",
            )
            return

        width = float(self.width())
        for row in layout.rows:
            center = self._widget_center(row.center_y)
            row_rect = QtCore.QRectF(0.0, center - row.pitch / 2.0, width, row.pitch)
            if row_rect.bottom() < 0 or row_rect.top() > self.height():
                continue
            if row.row_index % 2 and self.owner._display_mode == "exploratory":
                painter.fillRect(row_rect, ROW_ALT_BACKGROUND)
            if self.owner.hovered_row_index == row.row_index or self.owner.row_is_linked(row):
                # Same tint as the plot row, so a name can be followed to its genes.
                painter.fillRect(row_rect, ROW_HOVER_COLOR)
            _paint_rail_row(painter, row, center, width, show_molecule=self.owner.show_molecule_badges)
            painter.setPen(QtGui.QPen(SEPARATOR_COLOR, 1.0))
            painter.drawLine(QtCore.QPointF(0.0, row_rect.bottom()), QtCore.QPointF(width, row_rect.bottom()))
        for divider in layout.dividers:
            center = self._widget_center(divider.center_y)
            rect = QtCore.QRectF(0.0, center - divider.height / 2.0, width, divider.height)
            if rect.bottom() < 0 or rect.top() > self.height():
                continue
            self.owner._paint_divider_band(painter, rect, divider.text)
        fold = layout.fold
        if fold is not None:
            center = self._widget_center(fold.center_y)
            rect = QtCore.QRectF(0.0, center - fold.height / 2.0, width, fold.height)
            if rect.bottom() >= 0 and rect.top() <= self.height():
                painter.fillRect(rect, GAP_BACKGROUND)
                painter.setPen(QtGui.QPen(SEPARATOR_COLOR, 1.0))
                painter.drawLine(rect.topLeft(), rect.topRight())
                painter.drawLine(rect.bottomLeft(), rect.bottomRight())
                _paint_band_label(
                    painter,
                    rect.adjusted(RAIL_BADGE_X, 0.0, -RAIL_RIGHT_MARGIN, 0.0),
                    "Contigs with no neighbors",
                )

        painter.setPen(QtGui.QPen(BASELINE_COLOR, 1.0))
        painter.drawLine(int(width) - 1, 0, int(width) - 1, self.height())


def _draw_contig_end(
    painter: QtGui.QPainter,
    x: float,
    center_y: float,
    gene_height: float,
    scale_x: float,
    direction: float,
    device_scale: float = 1.0,
) -> None:
    """A bracket-shaped end cap: the bar marks the end, the serifs point
    back into the contig."""
    half = gene_height / 2.0 + 3.0
    serif = 4.0 / max(scale_x, 1e-6)
    pen = QtGui.QPen(CONTIG_END_COLOR, 2.0 * max(device_scale, 1.0))
    pen.setCosmetic(True)
    pen.setCapStyle(QtCore.Qt.PenCapStyle.FlatCap)
    painter.save()
    painter.setPen(pen)
    painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
    path = QtGui.QPainterPath(QtCore.QPointF(x + direction * serif, center_y - half))
    path.lineTo(x, center_y - half)
    path.lineTo(x, center_y + half)
    path.lineTo(x + direction * serif, center_y + half)
    painter.drawPath(path)
    painter.restore()


class _ContextRowItem(QtWidgets.QGraphicsItem):
    def __init__(self, owner: "ContextCanvas", row: RowLayout, scene_min_x: float, scene_max_x: float) -> None:
        super().__init__()
        self.owner = owner
        self.row = row
        self._bounds = QtCore.QRectF(
            scene_min_x,
            row.center_y - row.pitch / 2.0,
            max(1.0, scene_max_x - scene_min_x),
            row.pitch,
        )
        self.setZValue(3.0)
        self.setFlag(
            QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemUsesExtendedStyleOption,
            True,
        )
        # Rows are opaque enough for a device cache and account for most of
        # the cost of shared-family highlighting in dense comparisons.
        self.setCacheMode(QtWidgets.QGraphicsItem.CacheMode.DeviceCoordinateCache)

    def boundingRect(self) -> QtCore.QRectF:
        return self._bounds

    def _fill_for_gene(self, gene: GeneLayout) -> QtGui.QColor:
        return self.owner.gene_fill(gene)

    def _should_label(self, gene: GeneLayout, scale_x: float) -> bool:
        # The hovered and pinned genes get their names from the spotlight
        # layer, so row labels never pop in and out under the pointer.
        if gene.gene.is_focal:
            # On screen the focal name sits on top of the focal band; a
            # figure names it on the reference row instead.
            return bool(self.owner.exporting and self.row.row.is_reference)
        if self.owner.label_mode == "off":
            return False
        is_reference = self.row.row.is_reference
        # The reference row is named at every zoom (typed genes first, then
        # as many others as fit); other rows name their AMR, mobile and
        # virulence genes once zoomed in.
        key_label = bool(
            (is_reference and gene.gene.label_priority >= 60)
            or (gene.gene.category in KEY_LABEL_CATEGORIES and scale_x >= KEY_LABELS_MIN_SCALE)
        )
        if self.owner.label_mode == "all":
            # "All" adds every named gene from the opening zoom onward (labels
            # are collision checked), and never drops a key gene label.
            return key_label or (scale_x >= ALL_LABELS_MIN_SCALE and bool(gene.gene.label))
        return key_label

    def _draw_fixed_text(
        self,
        painter: QtGui.QPainter,
        x: float,
        y: float,
        text: str,
        font: QtGui.QFont,
        color: QtGui.QColor,
        *,
        centered: bool = True,
    ) -> QtCore.QRectF:
        transform = painter.worldTransform()
        device_scale = max(float(self.owner._paint_device_scale), 1e-6)
        scale_x = max(abs(float(transform.m11())) / device_scale, 1e-6)
        scale_y = max(abs(float(transform.m22())) / device_scale, 1e-6)
        metrics = QtGui.QFontMetrics(font)
        width = float(metrics.horizontalAdvance(text))
        height = float(metrics.height())
        painter.save()
        painter.translate(x, y)
        painter.scale(1.0 / scale_x, 1.0 / scale_y)
        painter.setFont(font)
        painter.setPen(color)
        left = -width / 2.0 if centered else 0.0
        rect = QtCore.QRectF(left, -height / 2.0, width, height)
        painter.drawText(rect, QtCore.Qt.AlignmentFlag.AlignCenter, text)
        painter.restore()
        device_center = transform.map(QtCore.QPointF(x, y))
        return QtCore.QRectF(
            device_center.x() + left * device_scale,
            device_center.y() - height * device_scale / 2.0,
            width * device_scale,
            height * device_scale,
        )

    def paint(
        self,
        painter: QtGui.QPainter,
        option: QtWidgets.QStyleOptionGraphicsItem,
        widget: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        del widget
        transform = painter.worldTransform()
        device_scale = max(float(self.owner._paint_device_scale), 1e-6)
        scale_x = max(abs(float(transform.m11())) / device_scale, 1e-6)
        exposed = option.exposedRect.adjusted(-12.0 / scale_x, -3.0, 12.0 / scale_x, 3.0)
        baseline = (
            PUBLICATION_BASELINE
            if self.owner._display_mode in {"figure", "publication"}
            else BASELINE_COLOR
        )
        painter.setPen(QtGui.QPen(baseline, 1.0))
        painter.drawLine(
            QtCore.QPointF(self.row.min_x, self.row.center_y),
            QtCore.QPointF(self.row.max_x, self.row.center_y),
        )
        # An end cap where the contig really ends, so 'the assembly stops
        # here' reads differently from 'nothing similar here'.
        gap = min(5.0 / scale_x, 20.0)
        for edge, present, direction in (
            (self.row.min_x - gap, self.row.cap_left, 1.0),
            (self.row.max_x + gap, self.row.cap_right, -1.0),
        ):
            if present:
                _draw_contig_end(painter, edge, self.row.center_y, self.row.gene_height, scale_x, direction, device_scale)

        for gene in self.row.genes:
            top = self.row.center_y - self.row.gene_height / 2.0
            rect = QtCore.QRectF(gene.left, top, max(gene.width, 0.01), self.row.gene_height)
            if not rect.intersects(exposed):
                continue
            path = gene_arrow_path(rect, gene.display_strand, scale_x)
            fill = self._fill_for_gene(gene)
            fill.setAlpha(245)
            painter.setBrush(QtGui.QBrush(fill))
            # Hover and pin outlines live in the spotlight layer on top.
            if gene.gene.is_focal:
                pen = QtGui.QPen(FOCAL_COLOR, 2.0 * device_scale)
            elif self.owner._display_mode in {"figure", "publication"}:
                outline = QtGui.QColor(COLORS["text"])
                outline.setAlpha(205)
                pen = QtGui.QPen(outline, 0.8 * device_scale)
            else:
                outline = QtGui.QColor(COLORS["surface"])
                outline.setAlpha(215)
                pen = QtGui.QPen(outline, 0.8 * device_scale)
            pen.setCosmetic(True)
            painter.setPen(pen)
            painter.drawPath(path)

        # Labels are priority sorted and collision checked in device pixels.
        candidates = [gene for gene in self.row.genes if gene.gene.label and self._should_label(gene, scale_x)]
        candidates.sort(key=lambda gene: (-gene.gene.label_priority, gene.center_x))
        occupied: List[QtCore.QRectF] = []
        for gene in candidates:
            gene_bounds = QtCore.QRectF(
                gene.left,
                self.row.center_y - self.row.gene_height / 2.0 - 16.0,
                max(gene.width, 0.01),
                self.row.gene_height + 18.0,
            )
            if not gene_bounds.intersects(exposed):
                continue
            font = QtGui.QFont()
            font.setPixelSize(MIN_TEXT_PX)
            font.setBold(gene.gene.is_focal)
            metrics = QtGui.QFontMetrics(font)
            text_width = float(metrics.horizontalAdvance(gene.gene.label))
            gene_screen_width = gene.width * scale_x
            if gene_screen_width >= text_width + 12.0:
                fill = self._fill_for_gene(gene)
                text_color = QtGui.QColor(COLORS["surface"]) if fill.lightness() < 145 else TEXT_COLOR
                label_rect = self._draw_fixed_text(
                    painter,
                    gene.center_x,
                    self.row.center_y,
                    gene.gene.label,
                    font,
                    text_color,
                )
            else:
                label_y = self.row.center_y - self.row.gene_height / 2.0 - 8.0
                device_center = transform.map(QtCore.QPointF(gene.center_x, label_y))
                candidate_rect = QtCore.QRectF(
                    device_center.x() - (text_width / 2.0 + 3.0) * device_scale,
                    device_center.y() - metrics.height() * device_scale / 2.0,
                    (text_width + 6.0) * device_scale,
                    float(metrics.height()) * device_scale,
                )
                if any(candidate_rect.intersects(existing) for existing in occupied):
                    continue
                label_rect = self._draw_fixed_text(
                    painter,
                    gene.center_x,
                    label_y,
                    gene.gene.label,
                    font,
                    TEXT_COLOR,
                )
            occupied.append(label_rect.adjusted(-3.0, -1.0, 3.0, 1.0))


# Resistance mutation badge: screen pixels, on the gene's lower right corner
# (labels sit inside or above a gene, so the corner below stays free).
MUTATION_BADGE_PX = 12.0


def badge_inset_px(width_px: float, plus_strand: bool) -> float:
    """How far left of a gene's right end its badge sits, in screen pixels:
    under the arrow's shoulder on the plus strand, by the tail's corner on
    the minus strand, never past the middle of a short gene."""
    head = min(width_px * 0.42, 8.0)
    return min(head if plus_strand else MUTATION_BADGE_PX / 2.0, width_px / 2.0)


class _MutationBadgeLayer(QtWidgets.QGraphicsItem):
    """The resistance mutation badges of every row, above the spotlight.

    One item for the plot: it draws a fixed-size badge at each gene with a
    ``mutation_note`` (the gene's tooltip carries the sentence), on screen and
    in exports alike. Only the exposed badges are painted; on screen they
    come from one cached pixmap, exports render the SVG. While a family is
    spotlit the badges of the other genes fade like the veiled genes.
    """

    def __init__(self, owner: "ContextCanvas", layout: ContextLayout) -> None:
        super().__init__()
        self.owner = owner
        self._genes = [
            (gene, row.gene_height)
            for row in layout.rows
            for gene in row.genes
            if gene.gene.mutation_note
        ]
        self._bounds = QtCore.QRectF(layout.min_x, layout.min_y, layout.width, layout.height)
        # ids of the spotlit genes while the veil shows, else None.
        self._lit: Optional[set] = None
        self._pixmap: Optional[QtGui.QPixmap] = None
        self.setZValue(6.2)
        self.setAcceptedMouseButtons(QtCore.Qt.MouseButton.NoButton)

    def boundingRect(self) -> QtCore.QRectF:
        return self._bounds

    def set_spotlit(self, genes: Optional[Iterable[GeneLayout]]) -> None:
        lit = None if genes is None else {id(gene) for gene in genes}
        if lit != self._lit:
            self._lit = lit
            self.update()

    def _badge_pixmap(self, ratio: float) -> QtGui.QPixmap:
        side = max(1, int(round(MUTATION_BADGE_PX * ratio)))
        pixmap = self._pixmap
        if pixmap is None or pixmap.width() != side:
            image = QtGui.QImage(side, side, QtGui.QImage.Format.Format_ARGB32_Premultiplied)
            image.fill(QtCore.Qt.GlobalColor.transparent)
            image_painter = QtGui.QPainter(image)
            image_painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
            paint_mutation_badge(image_painter, QtCore.QRectF(0.0, 0.0, side, side))
            image_painter.end()
            pixmap = QtGui.QPixmap.fromImage(image)
            self._pixmap = pixmap
        return pixmap

    def paint(
        self,
        painter: QtGui.QPainter,
        option: QtWidgets.QStyleOptionGraphicsItem,
        widget: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        del widget
        transform = painter.worldTransform()
        device_scale = max(float(self.owner._paint_device_scale), 1e-6)
        scale_x = max(abs(float(transform.m11())) / device_scale, 1e-6)
        scale_y = max(abs(float(transform.m22())) / device_scale, 1e-6)
        size = MUTATION_BADGE_PX
        half_x = size / 2.0 / scale_x
        half_y = size / 2.0 / scale_y
        exposed = option.exposedRect.adjusted(-half_x, -half_y, half_x, half_y)
        bounds = self._bounds
        pixmap = None
        if not self.owner.exporting:
            pixmap = self._badge_pixmap(device_scale * painter.device().devicePixelRatioF())
        veiled_opacity = 1.0 - SPOTLIGHT_VEIL_ALPHA / 255.0
        lit = self._lit
        target = QtCore.QRectF(-size / 2.0, -size / 2.0, size, size)
        for gene, gene_height in self._genes:
            inset = badge_inset_px(gene.width * scale_x, gene.display_strand != "-")
            # Kept inside the plot, so nothing is drawn past the bounds.
            x = min(gene.left + gene.width - inset / scale_x, bounds.right() - half_x)
            y = min(gene.center_y + gene_height / 2.0, bounds.bottom() - half_y)
            if not exposed.contains(QtCore.QPointF(x, y)):
                continue
            painter.save()
            if lit is not None and id(gene) not in lit:
                painter.setOpacity(veiled_opacity)
            painter.translate(x, y)
            painter.scale(1.0 / scale_x, 1.0 / scale_y)
            if pixmap is not None:
                painter.drawPixmap(target, pixmap, QtCore.QRectF(pixmap.rect()))
            else:
                paint_mutation_badge(painter, target)
            painter.restore()


class _SyntenyPairItem(QtWidgets.QGraphicsItem):
    def __init__(self, owner: "ContextCanvas", pair: PairLayout) -> None:
        super().__init__()
        self.owner = owner
        self.pair = pair
        self._paths: List[Tuple[LinkLayout, QtGui.QPainterPath]] = []
        bounds = QtCore.QRectF()
        for link in pair.links:
            upper = link.upper_gene
            lower = link.lower_gene
            if upper.gene.is_focal and lower.gene.is_focal:
                # Focal genes already line up under the focal band; stacking
                # their ribbons only drew an accidental grey column.
                continue
            y_upper = upper.center_y + owner.gene_height_for_row(upper.row_index) / 2.0 + 0.5
            y_lower = lower.center_y - owner.gene_height_for_row(lower.row_index) / 2.0 - 0.5
            mid_y = (y_upper + y_lower) / 2.0
            path = QtGui.QPainterPath(QtCore.QPointF(upper.left, y_upper))
            path.cubicTo(upper.left, mid_y, lower.left, mid_y, lower.left, y_lower)
            path.lineTo(lower.right, y_lower)
            path.cubicTo(lower.right, mid_y, upper.right, mid_y, upper.right, y_upper)
            path.closeSubpath()
            self._paths.append((link, path))
            bounds = bounds.united(path.boundingRect()) if not bounds.isNull() else path.boundingRect()
        self._bounds = bounds.adjusted(-2.0, -2.0, 2.0, 2.0) if not bounds.isNull() else QtCore.QRectF()
        self.setZValue(1.0)
        self.setFlag(
            QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemUsesExtendedStyleOption,
            True,
        )
        # A device cache on these translucent ribbons can produce opaque black
        # tiles on the Qt offscreen/raster backend.  Item-coordinate caching is
        # stable and still avoids replaying every ribbon during interaction.
        self.setCacheMode(QtWidgets.QGraphicsItem.CacheMode.ItemCoordinateCache)

    def boundingRect(self) -> QtCore.QRectF:
        return self._bounds

    def _link_is_visible(self, link: LinkLayout) -> bool:
        selected = bool(
            self.owner.selected_family
            and link.family_key == self.owner.selected_family
            and (not self.owner.exporting or self.owner.link_mode == "selected")
        )
        hovered = bool(
            not self.owner.exporting
            and self.owner.hovered_family
            and link.family_key == self.owner.hovered_family
        )
        if self.owner.link_mode == "off" or not self.owner.show_links:
            return False
        if self.owner.link_mode == "selected":
            return selected or hovered
        if self.owner.link_mode == "blocks":
            layout = self.owner.layout_result
            scene_width = layout.width if layout is not None else 800.0
            max_local_shift = max(90.0, scene_width * 0.18)
            if abs(link.upper_gene.center_x - link.lower_gene.center_x) > max_local_shift:
                return False
        return True

    def paint(
        self,
        painter: QtGui.QPainter,
        option: QtWidgets.QStyleOptionGraphicsItem,
        widget: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        del widget
        for link, path in self._paths:
            if not path.boundingRect().intersects(option.exposedRect) or not self._link_is_visible(link):
                continue
            publication = self.owner._display_mode in {"figure", "publication"}
            if link.inverted:
                color = QtGui.QColor(*INVERTED_LINK_RGB, 26 if publication else 32)
                outline = QtGui.QColor(183, 70, 75, 50)
            else:
                color = QtGui.QColor(*LINK_RGB, 26 if publication else 34)
                outline = QtGui.QColor(*LINK_RGB, 0)
            if self.owner.link_mode == "selected":
                color.setAlpha(105)
                outline.setAlpha(145)
            painter.setBrush(QtGui.QBrush(color))
            pen = QtGui.QPen(outline, 0.8 * max(float(self.owner._paint_device_scale), 1.0))
            pen.setCosmetic(True)
            painter.setPen(pen)
            painter.drawPath(path)

    def link_at(self, scene_point: QtCore.QPointF) -> Optional[LinkLayout]:
        for link, path in reversed(self._paths):
            if self._link_is_visible(link) and path.contains(scene_point):
                return link
        return None


class _DynamicHighlightPathItem(QtWidgets.QGraphicsPathItem):
    """Small interaction-only path that is omitted from document exports."""

    def __init__(self, owner: "ContextCanvas", path: QtGui.QPainterPath) -> None:
        super().__init__(path)
        self.owner = owner
        self.setZValue(6.0)
        self.setAcceptedMouseButtons(QtCore.Qt.MouseButton.NoButton)

    def paint(
        self,
        painter: QtGui.QPainter,
        option: QtWidgets.QStyleOptionGraphicsItem,
        widget: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        if not self.owner.exporting:
            super().paint(painter, option, widget)


class _InteractionRectItem(QtWidgets.QGraphicsRectItem):
    """Row hover tint or spotlight veil; never part of a figure."""

    def __init__(self, owner: "ContextCanvas", color: QtGui.QColor, z_value: float) -> None:
        super().__init__()
        self.owner = owner
        self.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
        self.setBrush(QtGui.QBrush(color))
        self.setZValue(z_value)
        self.setAcceptedMouseButtons(QtCore.Qt.MouseButton.NoButton)
        self.hide()

    def paint(
        self,
        painter: QtGui.QPainter,
        option: QtWidgets.QStyleOptionGraphicsItem,
        widget: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        if not self.owner.exporting:
            super().paint(painter, option, widget)


class _SpotlightLabelItem(QtWidgets.QGraphicsItem):
    """Screen-sized name pill above the hovered or pinned gene."""

    def __init__(self, owner: "ContextCanvas", text: str, anchor: QtCore.QPointF, *, pinned: bool) -> None:
        super().__init__()
        self.owner = owner
        self.text = text
        self.pinned = bool(pinned)
        self.font = QtGui.QFont()
        self.font.setPixelSize(11)
        self.font.setBold(True)
        metrics = QtGui.QFontMetricsF(self.font)
        width = metrics.horizontalAdvance(text) + 12.0
        height = metrics.height() + 4.0
        self._rect = QtCore.QRectF(-width / 2.0, -height - 3.0, width, height)
        self.setFlag(QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations, True)
        self.setAcceptedMouseButtons(QtCore.Qt.MouseButton.NoButton)
        self.setPos(anchor)
        self.setZValue(7.0)

    def boundingRect(self) -> QtCore.QRectF:
        return self._rect.adjusted(-1.0, -1.0, 1.0, 1.0)

    def paint(
        self,
        painter: QtGui.QPainter,
        option: QtWidgets.QStyleOptionGraphicsItem,
        widget: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        del option, widget
        if self.owner.exporting:
            return
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        if self.pinned:
            painter.setPen(QtGui.QPen(SELECTED_COLOR, 1.0))
            painter.setBrush(SELECTED_COLOR)
            text_color = QtGui.QColor(COLORS["surface"])
        else:
            painter.setPen(QtGui.QPen(QtGui.QColor(COLORS["border_strong"]), 1.0))
            painter.setBrush(QtGui.QColor(COLORS["surface"]))
            text_color = STRONG_TEXT_COLOR
        painter.drawRoundedRect(self._rect, 5.0, 5.0)
        painter.setFont(self.font)
        painter.setPen(text_color)
        painter.drawText(self._rect, QtCore.Qt.AlignmentFlag.AlignCenter, self.text)


class _GeneNameItem(QtWidgets.QGraphicsItem):
    """Screen-sized gene name redrawn above the spotlight veil."""

    def __init__(self, owner: "ContextCanvas", text: str, anchor: QtCore.QPointF, color: QtGui.QColor) -> None:
        super().__init__()
        self.owner = owner
        self.text = text
        self.color = QtGui.QColor(color)
        self.font = QtGui.QFont()
        self.font.setPixelSize(MIN_TEXT_PX)
        metrics = QtGui.QFontMetricsF(self.font)
        width = metrics.horizontalAdvance(text)
        height = metrics.height()
        self._rect = QtCore.QRectF(-width / 2.0, -height / 2.0, width, height)
        self.setFlag(QtWidgets.QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations, True)
        self.setAcceptedMouseButtons(QtCore.Qt.MouseButton.NoButton)
        self.setPos(anchor)
        self.setZValue(6.5)

    def boundingRect(self) -> QtCore.QRectF:
        return self._rect.adjusted(-1.0, -1.0, 1.0, 1.0)

    def paint(
        self,
        painter: QtGui.QPainter,
        option: QtWidgets.QStyleOptionGraphicsItem,
        widget: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        del option, widget
        if self.owner.exporting:
            return
        painter.setFont(self.font)
        painter.setPen(self.color)
        painter.drawText(self._rect, QtCore.Qt.AlignmentFlag.AlignCenter, self.text)


class ContextCanvas(QtWidgets.QWidget):
    """A plot-first comparative neighborhood canvas with a frozen row rail."""

    contextSamplesHovered = QtCore.pyqtSignal(object)
    geneDoubleClicked = QtCore.pyqtSignal(dict)
    geneSelected = QtCore.pyqtSignal(dict)
    # Kept for the dialog contract shared with the classic renderer; ribbon
    # clicks now select the family instead of emitting this.
    syntenyBlockClicked = QtCore.pyqtSignal(dict)
    viewportChanged = QtCore.pyqtSignal()
    # Right-click on a row: (the row's source dict, global position).  The
    # dialog builds the menu, since it owns the reference and the samples.
    rowMenuRequested = QtCore.pyqtSignal(object, QtCore.QPoint)

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(parent)
        self.setObjectName("ContextCanvas")
        self.setMinimumHeight(320)
        self._model = ContextRenderModel(rows=())
        self.layout_result: Optional[ContextLayout] = None
        self._engine = ContextLayoutEngine()
        self._density = "compact"
        self._display_mode = "exploratory"
        self._normalize_orientation = True
        self.show_links = True
        self.link_mode = "blocks"
        # "reference" is the "Key genes" labels mode.
        self.label_mode = "reference"
        # Gene type ("annotation") is the default; "family" colors by homolog.
        self.color_mode = "annotation"
        self.show_molecule_badges = True
        # Shown in the plot when no row is left ('' keeps the default).
        self.empty_text = ""
        # Contigs holding only the focal gene fold into one summary line;
        # "Show" draws them as thin rows under it.
        self.short_rows_expanded = False
        # Export > Include short contigs: figures draw them as thin rows.
        self.export_include_short_contigs = True
        # Labelled dividers between plasmid and chromosome blocks.
        self.split_molecules = False
        self._short_rows_mode_override: Optional[str] = None
        self._fold_link_rect = QtCore.QRectF()
        self.exporting = False
        # Publication style applies to figures only and never to the screen.
        self.export_publication_style = False
        self.legend_collapsed = False
        # Set once the viewer opens the key on purpose: it then stays open
        # even where it has to cover genes.
        self._legend_forced_open = False
        self._legend_auto_folded = False
        self._legend_rect = QtCore.QRectF()
        # Gene type databases not installed for the run: a muted key line and tooltip.
        self._not_installed_text = ""
        self._not_installed_tooltip = ""
        self._paint_device_scale = 1.0
        self.hovered_family: Optional[str] = None
        self._hovered_gene: Optional[GeneLayout] = None
        self.selected_family: Optional[str] = None
        self.hovered_row_index: Optional[int] = None
        self._selected_gene: Optional[GeneLayout] = None
        # Last pointer position over the plot or the rail, so a spotlight can
        # follow whatever is under it after a scroll or zoom.
        self._pointer: Optional[Tuple[str, QtCore.QPointF]] = None
        self._row_items: List[_ContextRowItem] = []
        self._row_background_items: List[QtWidgets.QGraphicsRectItem] = []
        self._pair_items: List[_SyntenyPairItem] = []
        self._family_pair_items: Dict[str, List[_SyntenyPairItem]] = {}
        self._family_link_paths: Dict[str, List[Tuple[LinkLayout, QtGui.QPainterPath]]] = {}
        self._family_gene_layouts: Dict[str, List[GeneLayout]] = {}
        self._family_samples: Dict[str, set] = {}
        self._all_genomes: set = set()
        self._genomes_with_neighbors: set = set()
        # True while row filters hide some genomes: the counts below then
        # cover only the rows shown, and the text says so.
        self.rows_filtered = False
        self._has_focal_genes = False
        self._dynamic_highlight_items: List[QtWidgets.QGraphicsItem] = []
        self._row_hover_item: Optional[_InteractionRectItem] = None
        # Samples hovered in another panel of the window: their rows get the
        # row hover tint (set_linked_samples).
        self.linked_samples: frozenset = frozenset()
        self._linked_row_items: List[_InteractionRectItem] = []
        self._veil_item: Optional[_InteractionRectItem] = None
        self._badge_layer: Optional[_MutationBadgeLayer] = None
        self._focal_band_rect: Optional[QtCore.QRectF] = None
        # Every focal gene is drawn in the reference focal gene's color, so a
        # partial or untyped hit still reads as the gene being compared.
        self._focal_source: Optional[GeneLayout] = None
        self.focal_label = ""
        self._last_hover_samples: Optional[Tuple[str, ...]] = None
        self.export_title = "Gene neighborhood comparison"

        root = QtWidgets.QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)
        self.splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal, self)
        self.splitter.setChildrenCollapsible(False)
        # Both columns carry a header strip of the same height so rail rows
        # stay level with the plot rows.
        self.rail_column = QtWidgets.QWidget(self)
        self.rail_column.setMinimumWidth(200)
        self.rail_column.setMaximumWidth(360)
        rail_layout = QtWidgets.QVBoxLayout(self.rail_column)
        rail_layout.setContentsMargins(0, 0, 0, 0)
        rail_layout.setSpacing(0)
        self.rail_caption = _ColumnHeader(self, rail=True)
        self.row_header = _RowHeaderRail(self)
        rail_layout.addWidget(self.rail_caption)
        rail_layout.addWidget(self.row_header, stretch=1)
        plot_column = QtWidgets.QWidget(self)
        plot_layout = QtWidgets.QVBoxLayout(plot_column)
        plot_layout.setContentsMargins(0, 0, 0, 0)
        plot_layout.setSpacing(0)
        self.plot_caption = _ColumnHeader(self, rail=False)
        self.graphics_view = _GenomeView(self)
        self._scene = QtWidgets.QGraphicsScene(self.graphics_view)
        self.graphics_view.setScene(self._scene)
        plot_layout.addWidget(self.plot_caption)
        plot_layout.addWidget(self.graphics_view, stretch=1)
        self.splitter.addWidget(self.rail_column)
        self.splitter.addWidget(plot_column)
        self.splitter.setStretchFactor(0, 0)
        self.splitter.setStretchFactor(1, 1)
        self.splitter.setSizes([RAIL_WIDTH, 1100])
        root.addWidget(self.splitter)

        self.graphics_view.horizontalScrollBar().valueChanged.connect(self._viewport_changed)
        self.graphics_view.verticalScrollBar().valueChanged.connect(self._viewport_changed)
        self.splitter.splitterMoved.connect(lambda _pos, _index: self._viewport_changed())

    def clear(self) -> None:
        self._emit_hover_samples(None)
        self._model = ContextRenderModel(rows=())
        self.layout_result = None
        self.hovered_family = None
        self._hovered_gene = None
        self.hovered_row_index = None
        self.selected_family = None
        self._selected_gene = None
        self._reset_scene_bookkeeping()
        self._scene.clear()
        self.graphics_view.resetTransform()
        self.row_header.update()
        self.viewportChanged.emit()

    def set_rows(self, rows: Sequence[Dict[str, Any]]) -> None:
        self._emit_hover_samples(None)
        self._model = ContextRenderModel.from_rows(_mark_repeated_contigs(rows))
        self.hovered_family = None
        self._hovered_gene = None
        self.hovered_row_index = None
        self._rebuild_layout()
        QtCore.QTimer.singleShot(0, self.reset_zoom)

    def _reset_scene_bookkeeping(self) -> None:
        self._row_items = []
        self._row_background_items = []
        self._pair_items = []
        self._family_pair_items = {}
        self._family_link_paths = {}
        self._family_gene_layouts = {}
        self._family_samples = {}
        self._all_genomes = set()
        self._genomes_with_neighbors = set()
        self._has_focal_genes = False
        self._badge_layer = None
        self._dynamic_highlight_items = []
        self._row_hover_item = None
        self._linked_row_items = []
        self._veil_item = None
        self._focal_band_rect = None
        self._focal_source: Optional[GeneLayout] = None
        self.focal_label = ""

    def _rebuild_layout(self, *, restore_identity: Optional[Tuple[str, str, str, float, float]] = None) -> None:
        previous_gene = self._selected_gene
        previous_identity = restore_identity or self._gene_identity(previous_gene, self.layout_result)
        self.layout_result = self._engine.layout(
            self._model,
            density=self._density,
            normalize_orientation=self._normalize_orientation,
            short_rows=self._short_rows_mode(),
            split_molecules=self.split_molecules,
        )
        self._scene.clear()
        self._reset_scene_bookkeeping()
        layout = self.layout_result
        self._scene.setSceneRect(layout.min_x, layout.min_y, layout.width, layout.height)
        self._collect_presence(layout)
        if not layout.rows:
            self.selected_family = None
            self._selected_gene = None
            if previous_gene is not None:
                self.geneSelected.emit({})
            message = self._scene.addSimpleText(self.empty_text or "No neighborhoods match the current filters")
            font = message.font()
            font.setPointSize(11)
            message.setFont(font)
            message.setBrush(MUTED_TEXT_COLOR)
            bounds = message.boundingRect()
            center = self._scene.sceneRect().center()
            message.setPos(
                center.x() - bounds.width() / 2.0,
                center.y() - bounds.height() / 2.0,
            )
            self.row_header.update()
            self.plot_caption.refresh()
            return

        for row in layout.rows:
            rect = QtCore.QRectF(
                layout.min_x,
                row.center_y - row.pitch / 2.0,
                layout.width,
                row.pitch,
            )
            background = self._scene.addRect(
                rect,
                QtGui.QPen(QtCore.Qt.PenStyle.NoPen),
                QtGui.QBrush(self._row_background_color(row.row_index)),
            )
            background.setZValue(0.0)
            self._row_background_items.append(background)
        # Divider and short-contig bands span the plot; their text sits in
        # the rail and, for the summary line, in the view overlay.
        bands = [(divider.center_y, divider.height, HEADER_BACKGROUND) for divider in layout.dividers]
        if layout.fold is not None:
            bands.append((layout.fold.center_y, layout.fold.height, GAP_BACKGROUND))
        for center_y, height, color in bands:
            band = self._scene.addRect(
                QtCore.QRectF(layout.min_x, center_y - height / 2.0, layout.width, height),
                QtGui.QPen(SEPARATOR_COLOR, 0.0),
                QtGui.QBrush(color),
            )
            band.setZValue(0.6)
        self._row_hover_item = _InteractionRectItem(self, ROW_HOVER_COLOR, 0.3)
        self._scene.addItem(self._row_hover_item)

        # One light band behind every row marks the focal gene column.
        focal_genes = [gene for row in layout.rows for gene in row.genes if gene.gene.is_focal]
        if focal_genes:
            left = min(gene.left for gene in focal_genes)
            right = max(gene.right for gene in focal_genes)
            pad = max(4.0, (right - left) * 0.12)
            # The band stops at the last row or the summary line, not in the
            # empty margin under them.
            bottom = max(row.center_y + row.pitch / 2.0 for row in layout.rows)
            if layout.fold is not None:
                bottom = max(bottom, layout.fold.center_y + layout.fold.height / 2.0)
            self._focal_band_rect = QtCore.QRectF(
                left - pad, layout.min_y, (right - left) + 2.0 * pad, bottom - layout.min_y
            )
            band = self._scene.addRect(
                self._focal_band_rect,
                QtGui.QPen(QtCore.Qt.PenStyle.NoPen),
                QtGui.QBrush(FOCAL_BAND_COLOR),
            )
            band.setZValue(0.5)
            reference = layout.reference_row
            reference_focal = next(
                (gene for gene in (reference.genes if reference is not None else ()) if gene.gene.is_focal),
                focal_genes[0],
            )
            self.focal_label = reference_focal.gene.label
            self._focal_source = (
                reference_focal
                if reference_focal.gene.category
                else next((gene for gene in focal_genes if gene.gene.category), reference_focal)
            )
            self._has_focal_genes = True

        for pair in layout.pairs:
            item = _SyntenyPairItem(self, pair)
            self._scene.addItem(item)
            self._pair_items.append(item)
            for link, path in item._paths:
                if link.family_key:
                    self._family_link_paths.setdefault(link.family_key, []).append((link, path))
            for family in {link.family_key for link, _path in item._paths if link.family_key}:
                self._family_pair_items.setdefault(family, []).append(item)
        for row in layout.rows:
            item = _ContextRowItem(self, row, layout.min_x, layout.max_x)
            self._scene.addItem(item)
            self._row_items.append(item)
            for gene in row.genes:
                if gene.gene.family_key:
                    self._family_gene_layouts.setdefault(gene.gene.family_key, []).append(gene)
        self._veil_item = _InteractionRectItem(
            self,
            QtGui.QColor(*CANVAS_BACKGROUND.getRgb()[:3], SPOTLIGHT_VEIL_ALPHA),
            4.0,
        )
        self._veil_item.setRect(QtCore.QRectF(layout.min_x, layout.min_y, layout.width, layout.height))
        self._scene.addItem(self._veil_item)
        badges = _MutationBadgeLayer(self, layout)
        self._badge_layer = badges if badges._genes else None
        if self._badge_layer is not None:
            self._scene.addItem(badges)
        restored = None
        if previous_identity is not None:
            candidates = [
                gene
                for row in layout.rows
                for gene in row.genes
                if self._gene_identity(gene, layout) == previous_identity
            ]
            if candidates:
                restored = candidates[0]
        self.selected_family = (restored.gene.family_key or None) if restored is not None else None
        self._selected_gene = restored
        if previous_gene is not None and restored is None:
            self.geneSelected.emit({})
        elif restored is not None:
            self._emit_selected_gene(restored)
        self._rebuild_interaction_highlights()
        self._update_linked_row_items()
        self.row_header.update()
        self.plot_caption.refresh()
        self.viewportChanged.emit()

    def row_is_linked(self, row: RowLayout) -> bool:
        """True when a sample of ``row`` is hovered in another panel."""
        linked = self.linked_samples
        return bool(linked) and any(str(sample).strip() in linked for sample in row.row.samples)

    def set_linked_samples(self, samples: Optional[Sequence[str]]) -> None:
        """Tint the rows of samples hovered elsewhere in the window (never in
        a figure; the canvas's own hover is unchanged)."""
        linked = frozenset(str(s).strip() for s in samples or () if str(s).strip())
        if linked == self.linked_samples:
            return
        self.linked_samples = linked
        self._update_linked_row_items()
        self.row_header.update()

    def _update_linked_row_items(self) -> None:
        for item in self._linked_row_items:
            if item.scene() is self._scene:
                self._scene.removeItem(item)
        self._linked_row_items = []
        layout = self.layout_result
        if layout is None or not self.linked_samples:
            return
        for row in layout.rows:
            if not self.row_is_linked(row):
                continue
            item = _InteractionRectItem(self, ROW_HOVER_COLOR, 0.3)
            item.setRect(QtCore.QRectF(layout.min_x, row.center_y - row.pitch / 2.0, layout.width, row.pitch))
            item.show()
            self._scene.addItem(item)
            self._linked_row_items.append(item)

    def _row_background_color(self, row_index: int, *, publication: Optional[bool] = None) -> QtGui.QColor:
        if publication is None:
            publication = self._display_mode in {"figure", "publication"}
        if publication:
            return PUBLICATION_BACKGROUND
        return ROW_ALT_BACKGROUND if row_index % 2 else CANVAS_BACKGROUND

    def _short_rows_mode(self) -> str:
        """How the layout treats contigs that hold only the focal gene."""
        if self._short_rows_mode_override is not None:
            mode = self._short_rows_mode_override
        else:
            mode = "expanded" if self.short_rows_expanded else "folded"
        rows = self._model.rows
        if rows and all(row.is_short for row in rows) and mode != "inline":
            # Nothing else to show (for example after a search): draw them.
            return "expanded"
        return mode

    def _collect_presence(self, layout: ContextLayout) -> None:
        """Genome and family counts over every row, folded ones included."""
        del layout
        for row in self._model.rows:
            samples = row.samples or ((row.sample,) if row.sample else ())
            self._all_genomes.update(samples)
            if row.has_neighbors:
                self._genomes_with_neighbors.update(samples)
            for gene in row.genes:
                if gene.family_key:
                    self._family_samples.setdefault(gene.family_key, set()).update(samples)

    def set_short_rows_expanded(self, expanded: bool) -> None:
        expanded = bool(expanded)
        if expanded == self.short_rows_expanded:
            return
        self.short_rows_expanded = expanded
        scroll = self.graphics_view.verticalScrollBar().value()
        transform = self.graphics_view.transform()
        horizontal = self.graphics_view.horizontalScrollBar().value()
        self._rebuild_layout()
        # Keep the view where it was: the summary line stays under the pointer.
        self.graphics_view.setTransform(transform)
        self.graphics_view.horizontalScrollBar().setValue(horizontal)
        self.graphics_view.verticalScrollBar().setValue(scroll)
        self._viewport_changed()

    def toggle_short_rows(self) -> None:
        self.set_short_rows_expanded(not self.short_rows_expanded)

    def set_export_include_short_contigs(self, include: bool) -> None:
        """Figures only: draw short contigs as thin rows, or leave them out."""
        self.export_include_short_contigs = bool(include)

    def set_molecule_split(self, enabled: bool) -> None:
        """Label the plasmid and chromosome blocks with dividers.

        The dialog orders the rows; this only draws the dividers, and takes
        effect with the next set_rows.
        """
        self.split_molecules = bool(enabled)

    def _fold_at_scene_y(self, y: float) -> bool:
        fold = self.layout_result.fold if self.layout_result is not None else None
        return bool(fold is not None and abs(y - fold.center_y) <= fold.height / 2.0)

    def _divider_at_scene_y(self, y: float) -> Optional[DividerLayout]:
        layout = self.layout_result
        if layout is None:
            return None
        for divider in layout.dividers:
            if abs(y - divider.center_y) <= divider.height / 2.0:
                return divider
        return None

    @staticmethod
    def divider_tooltip(divider: DividerLayout) -> str:
        tooltip = f"{divider.text}: each block is ordered on its own"
        if divider.folded:
            noun = "contig" if divider.folded == 1 else "contigs"
            tooltip += (
                f". {divider.folded} more {noun} of this type carry only the focal gene"
                " and are folded into the line at the bottom"
            )
        return tooltip

    def fold_sentence(self) -> str:
        fold = self.layout_result.fold if self.layout_result is not None else None
        return _fold_sentence(fold, self.focal_label) if fold is not None else ""

    def fold_tooltip(self) -> str:
        fold = self.layout_result.fold if self.layout_result is not None else None
        if fold is None:
            return ""
        noun = "this contig" if fold.contig_count == 1 else f"these {fold.contig_count} contigs"
        if fold.expanded:
            return f"Click to fold {noun} into one line again"
        return f"Contigs too short to compare. Click to show {noun} as thin rows"

    def _paint_divider_band(self, painter: QtGui.QPainter, rect: QtCore.QRectF, text: str) -> None:
        painter.save()
        painter.fillRect(rect, HEADER_BACKGROUND)
        painter.setPen(QtGui.QPen(SEPARATOR_COLOR, 1.0))
        painter.drawLine(rect.topLeft(), rect.topRight())
        painter.drawLine(rect.bottomLeft(), rect.bottomRight())
        painter.restore()
        _paint_band_label(
            painter,
            rect.adjusted(RAIL_BADGE_X, 0.0, -RAIL_RIGHT_MARGIN, 0.0),
            text,
            italic=False,
        )

    def _request_row_menu(self, row: RowLayout, global_pos: QtCore.QPoint) -> None:
        self.rowMenuRequested.emit(dict(row.row.source), QtCore.QPoint(global_pos))

    @staticmethod
    def _gene_identity(
        gene: Optional[GeneLayout],
        layout: Optional[ContextLayout],
    ) -> Optional[Tuple[str, str, str, float, float]]:
        if gene is None or layout is None or not (0 <= gene.row_index < len(layout.rows)):
            return None
        row = layout.rows[gene.row_index].row
        payload = gene.gene.payload
        annotation_id = next(
            (
                value
                for column in ("AnnotationID", "LocusTag", "GeneID")
                if (value := clean_text(payload.get(column)))
            ),
            "",
        )
        row_identity = "\x1f".join((row.sample, row.contig, row.title))
        return (
            row_identity,
            annotation_id,
            gene.gene.label,
            float(gene.gene.start),
            float(gene.gene.end),
        )

    def gene_height_for_row(self, row_index: int) -> float:
        layout = self.layout_result
        if layout is None or not (0 <= row_index < len(layout.rows)):
            return 16.0
        return layout.rows[row_index].gene_height

    def set_synteny_blocks_enabled(self, enabled: bool) -> None:
        self.show_links = bool(enabled)
        for item in self._pair_items:
            item.update()
        self._rebuild_interaction_highlights()
        self.graphics_view.viewport().update()

    def set_link_mode(self, mode: str) -> None:
        normalized = str(mode or "").strip().lower()
        self.link_mode = normalized if normalized in {"blocks", "selected", "off", "all"} else "blocks"
        self.show_links = self.link_mode != "off"
        for item in self._pair_items:
            item.update()
        self._rebuild_interaction_highlights()
        self.graphics_view.viewport().update()

    def set_label_mode(self, mode: str) -> None:
        normalized = str(mode or "").strip().lower()
        normalized = "reference" if normalized == "key" else normalized
        self.label_mode = normalized if normalized in {"reference", "all", "off"} else "reference"
        for item in self._row_items:
            item.update()

    def set_color_mode(self, mode: str) -> None:
        normalized = str(mode or "").strip().lower()
        normalized = "annotation" if normalized == "type" else normalized
        self.color_mode = normalized if normalized in {"family", "annotation"} else "annotation"
        for item in self._row_items:
            item.update()
        # Spotlighted genes are redrawn with their fill, and the key changes.
        self._rebuild_interaction_highlights()
        self.graphics_view.viewport().update()

    def set_export_publication_style(self, enabled: bool) -> None:
        """Figures only: white rows, dark outlines and quieter ribbons."""
        self.export_publication_style = bool(enabled)

    def set_not_installed_note(self, text: str, tooltip: str = "") -> None:
        """Name the gene type databases the run did not have in the key, so a type
        with no genes is not read as none found."""
        self._not_installed_text = str(text or "")
        self._not_installed_tooltip = str(tooltip or "")
        self.graphics_view.viewport().update()

    def set_legend_collapsed(self, collapsed: bool) -> None:
        self.legend_collapsed = bool(collapsed)
        self._legend_forced_open = not self.legend_collapsed
        self.graphics_view.viewport().update()

    def _toggle_legend(self) -> None:
        # A key folded to make room for genes opens on click, like one the
        # viewer folded.
        self.set_legend_collapsed(not (self.legend_collapsed or self._legend_auto_folded))

    def set_molecule_badges_visible(self, visible: bool) -> None:
        self.show_molecule_badges = bool(visible)
        self.row_header.update()

    def set_static_row_labels(self, enabled: bool, *, track_offset: float = 0.0) -> None:
        self.rail_column.setVisible(bool(enabled))
        if enabled and track_offset > 0:
            self.splitter.setSizes([int(max(200.0, min(360.0, track_offset))), max(1, self.width() - int(track_offset))])

    def set_header_corner_widget(self, widget: QtWidgets.QWidget) -> None:
        """Keep a small widget at the right end of the plot's header strip."""
        self.plot_caption.set_corner_widget(widget)

    def set_sort_caption(self, text: str) -> None:
        """Describe the current row order above the plot."""
        self.plot_caption.set_caption(text)

    def set_export_title(self, text: str) -> None:
        self.export_title = clean_text(text) or "Gene neighborhood comparison"

    def set_display_mode(self, mode: str) -> None:
        normalized = str(mode or "").strip().lower()
        normalized = normalized if normalized in {"exploratory", "figure", "publication"} else "exploratory"
        normalized = "figure" if normalized == "publication" else normalized
        self._display_mode = normalized
        publication = normalized == "figure"
        # Styling only: color, label and link choices are never reset here.
        self.graphics_view.setBackgroundBrush(
            QtGui.QBrush(PUBLICATION_BACKGROUND if publication else CANVAS_BACKGROUND)
        )
        self._rebuild_layout()
        self.row_header.update()
        self.graphics_view.viewport().update()

    def set_density_mode(self, mode: str) -> None:
        normalized = str(mode or "").strip().lower()
        normalized = normalized if normalized in {"compact", "comfortable"} else "compact"
        if normalized == self._density:
            return
        self._density = normalized
        self._rebuild_layout()
        QtCore.QTimer.singleShot(0, self.reset_zoom)

    def set_normalize_orientation(self, enabled: bool) -> None:
        normalized = bool(enabled)
        if normalized == self._normalize_orientation:
            return
        self._normalize_orientation = normalized
        self._rebuild_layout()
        QtCore.QTimer.singleShot(0, self.reset_zoom)

    def gene_fill(self, gene: GeneLayout) -> QtGui.QColor:
        """Fill for the current Color by choice: gene type or gene family."""
        if gene.gene.is_focal and self._focal_source is not None:
            gene = self._focal_source
        if self.color_mode == "annotation":
            # Untyped genes are "Other" so the gene type key stays true.
            value = gene.gene.annotation_fill if gene.gene.category else CATEGORY_COLORS["OTHER"]
        else:
            value = gene.gene.family_fill
        color = QtGui.QColor(value)
        return color if color.isValid() else QtGui.QColor(CATEGORY_COLORS["OTHER"])

    def _spotlight_groups(self) -> List[Tuple[str, Optional[str], GeneLayout]]:
        """(state, family, gene) for the pinned and the hovered spotlight.

        Hovering a focal gene only tints its row: the focal genes already line
        up under the focal band.
        """
        groups: List[Tuple[str, Optional[str], GeneLayout]] = []
        pinned = self._selected_gene
        if pinned is not None:
            groups.append(("pinned", self.selected_family, pinned))
        hovered = self._hovered_gene
        if (
            hovered is not None
            and not hovered.gene.is_focal
            and hovered is not pinned
            and not (self.hovered_family and self.hovered_family == self.selected_family)
        ):
            groups.append(("hovered", self.hovered_family, hovered))
        return groups

    def _update_row_hover_item(self) -> None:
        item = self._row_hover_item
        layout = self.layout_result
        if item is None or layout is None:
            return
        index = self.hovered_row_index
        if index is None or not (0 <= index < len(layout.rows)):
            item.hide()
            return
        row = layout.rows[index]
        item.setRect(QtCore.QRectF(layout.min_x, row.center_y - row.pitch / 2.0, layout.width, row.pitch))
        item.show()

    def _rebuild_interaction_highlights(self) -> None:
        for item in self._dynamic_highlight_items:
            if item.scene() is self._scene:
                self._scene.removeItem(item)
        self._dynamic_highlight_items = []
        self._update_row_hover_item()
        groups = self._spotlight_groups() if not self.exporting and self.layout_result is not None else []
        if self._veil_item is not None:
            self._veil_item.setVisible(bool(groups))
        if self._badge_layer is not None:
            self._badge_layer.set_spotlit(
                [
                    gene
                    for _state, family, anchor in groups
                    for gene in (self._family_gene_layouts.get(family, ()) if family else (anchor,))
                ]
                if groups
                else None
            )
        if not groups:
            return

        visible = self.graphics_view.mapToScene(
            self.graphics_view.viewport().rect()
        ).boundingRect()
        scale_x = max(abs(float(self.graphics_view.transform().m11())), 1e-6)

        def add(item: QtWidgets.QGraphicsItem) -> None:
            self._scene.addItem(item)
            self._dynamic_highlight_items.append(item)

        name_font = QtGui.QFont()
        name_font.setPixelSize(MIN_TEXT_PX)
        name_metrics = QtGui.QFontMetricsF(name_font)
        # Names placed above small genes, per row, in screen pixels.
        names_above: Dict[int, List[Tuple[float, float]]] = {}

        def add_name(gene: GeneLayout, top: float) -> None:
            """Keep a spotlit gene's name readable above the veil."""
            text = gene.gene.label
            if not text or self.label_mode == "off":
                return
            text_width = name_metrics.horizontalAdvance(text)
            if gene.width * scale_x >= text_width + 12.0:
                fill = self.gene_fill(gene)
                color = QtGui.QColor(COLORS["surface"]) if fill.lightness() < 145 else TEXT_COLOR
                add(_GeneNameItem(self, text, QtCore.QPointF(gene.center_x, gene.center_y), color))
                return
            center = gene.center_x * scale_x
            span = (center - text_width / 2.0 - 3.0, center + text_width / 2.0 + 3.0)
            taken = names_above.setdefault(gene.row_index, [])
            if any(span[0] < right and left < span[1] for left, right in taken):
                return
            taken.append(span)
            add(_GeneNameItem(self, text, QtCore.QPointF(gene.center_x, top - 8.0), TEXT_COLOR))

        for state, family, anchor in groups:
            pinned = state == "pinned"
            # The ribbons joining the family darken above the veil.
            if family and self.link_mode != "off":
                for link, path in self._family_link_paths.get(family, ()):
                    if not path.boundingRect().intersects(visible):
                        continue
                    rgb = INVERTED_LINK_RGB if link.inverted else LINK_RGB
                    ribbon = _DynamicHighlightPathItem(self, path)
                    ribbon.setZValue(4.5)
                    ribbon.setBrush(QtGui.QBrush(QtGui.QColor(*rgb, 120)))
                    ribbon.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
                    add(ribbon)
            genes = self._family_gene_layouts.get(family, ()) if family else (anchor,)
            for gene in genes:
                height = self.gene_height_for_row(gene.row_index)
                rect = QtCore.QRectF(
                    gene.left,
                    gene.center_y - height / 2.0,
                    max(gene.width, 0.01),
                    height,
                )
                if not rect.intersects(visible):
                    continue
                path = gene_arrow_path(rect, gene.display_strand, scale_x)
                is_anchor = gene is anchor
                if not pinned:
                    # Hover is a soft halo; a pin is a solid primary outline.
                    halo = _DynamicHighlightPathItem(self, path)
                    halo.setZValue(5.5)
                    halo.setBrush(QtGui.QBrush(QtCore.Qt.BrushStyle.NoBrush))
                    halo_pen = QtGui.QPen(HOVER_HALO_COLOR, 7.0 if is_anchor else 5.0)
                    halo_pen.setCosmetic(True)
                    halo_pen.setJoinStyle(QtCore.Qt.PenJoinStyle.RoundJoin)
                    halo.setPen(halo_pen)
                    add(halo)
                arrow = _DynamicHighlightPathItem(self, path)
                arrow.setZValue(6.0)
                arrow.setBrush(QtGui.QBrush(self.gene_fill(gene)))
                if pinned:
                    pen = QtGui.QPen(SELECTED_COLOR, 3.0 if is_anchor else 2.2)
                elif is_anchor:
                    pen = QtGui.QPen(HOVER_COLOR, 1.8)
                elif gene.gene.is_focal:
                    pen = QtGui.QPen(FOCAL_COLOR, 2.0)
                else:
                    pen = QtGui.QPen(CANVAS_BACKGROUND, 0.8)
                pen.setCosmetic(True)
                pen.setJoinStyle(QtCore.Qt.PenJoinStyle.MiterJoin)
                arrow.setPen(pen)
                add(arrow)
                if not is_anchor:
                    add_name(gene, rect.top())
            if anchor.gene.label and visible.contains(QtCore.QPointF(anchor.center_x, anchor.center_y)):
                top = anchor.center_y - self.gene_height_for_row(anchor.row_index) / 2.0
                add(_SpotlightLabelItem(
                    self,
                    anchor.gene.label,
                    QtCore.QPointF(anchor.center_x, top),
                    pinned=pinned,
                ))

    def _refresh_interaction_items(self, *families: Optional[str]) -> None:
        """Rebuild the spotlight layer; ribbons only repaint in 'selected' mode."""
        pair_targets: Dict[int, _SyntenyPairItem] = {}
        self._rebuild_interaction_highlights()
        if self.link_mode == "selected":
            for family in families:
                if not family:
                    continue
                for item in self._family_pair_items.get(family, ()):
                    pair_targets[id(item)] = item
        for item in pair_targets.values():
            item.update()

    def _pointer_target(self) -> Tuple[Optional[GeneLayout], Optional[RowLayout]]:
        """Gene and row under the last known pointer position."""
        layout = self.layout_result
        if self._pointer is None or layout is None:
            return None, None
        source, position = self._pointer
        if source == "rail":
            return None, self.row_header._row_at_widget_y(position.y())
        point = self.graphics_view.mapToScene(position.toPoint())
        gene = self._gene_at_scene(point)
        if gene is None and self._legend_rect.contains(position):
            return None, None
        return gene, layout.row_at(point.y())

    def _viewport_changed(self, *_args: Any) -> None:
        # Scrolling or zooming moves the plot under a still pointer: keep the
        # spotlight on whatever is under it now instead of wiping it.
        gene, row = self._pointer_target()
        self._set_hover(gene, row, refresh=False)
        self.row_header.update()
        self.plot_caption.refresh()
        self._rebuild_interaction_highlights()
        self.graphics_view.viewport().update()
        self.viewportChanged.emit()

    def zoom_by(self, factor: float) -> None:
        current = abs(float(self.graphics_view.transform().m11()))
        target = max(0.05, min(12.0, current * float(factor)))
        if current <= 0:
            return
        self.graphics_view.scale(target / current, 1.0)
        self._viewport_changed()

    def _set_horizontal_scale(self, scale_x: float) -> None:
        scale_x = max(0.05, min(12.0, float(scale_x)))
        transform = QtGui.QTransform()
        transform.scale(scale_x, 1.0)
        self.graphics_view.setTransform(transform)

    def reset_zoom(self) -> None:
        layout = self.layout_result
        if layout is None:
            return
        if not layout.rows:
            self.graphics_view.resetTransform()
            self.graphics_view.centerOn(self._scene.sceneRect().center())
            self._viewport_changed()
            return
        self._set_horizontal_scale(self._reference_fit_scale(layout))
        self.graphics_view.centerOn(0.0, layout.rows[0].center_y)
        self.graphics_view.verticalScrollBar().setValue(self.graphics_view.verticalScrollBar().minimum())
        self._viewport_changed()

    def view_state(self) -> Tuple[float, float, int, int]:
        """Zoom, the scene x at the plot's left edge, the vertical scroll and the plot width."""
        view = self.graphics_view
        left = float(view.mapToScene(QtCore.QPoint(0, 0)).x())
        return (
            float(view.transform().m11()),
            left,
            int(view.verticalScrollBar().value()),
            int(view.viewport().width()),
        )

    def restore_view_state(self, state: Tuple[float, float, int, int]) -> None:
        """Put back a view_state() after the plot changed width.

        The zoom and the rows in view stay; the left edge stays too, unless
        the narrower plot would push the focal column out of sight.
        """
        scale_x, left, vertical, width = state
        view = self.graphics_view
        if abs(float(view.transform().m11()) - scale_x) > 1e-9:
            self._set_horizontal_scale(scale_x)
        current_left = float(view.mapToScene(QtCore.QPoint(0, 0)).x())
        bar = view.horizontalScrollBar()
        bar.setValue(bar.value() + int(round((left - current_left) * scale_x)))
        view.verticalScrollBar().setValue(vertical)
        if view.viewport().width() < width:
            self._keep_focal_column_visible()
        self._viewport_changed()

    def _reference_fit_scale(self, layout: ContextLayout) -> float:
        """Horizontal scale that fits the reference row in the view."""
        reference = layout.reference_row
        span = (reference.max_x - reference.min_x) if reference is not None else layout.width
        span = max(240.0, span + 80.0)
        available = max(240.0, float(self.graphics_view.viewport().width()) - 36.0)
        return max(0.18, min(2.5, available / span))

    def _keep_focal_column_visible(self) -> None:
        """After a resize, bring the focal column back if it left the view.

        The view keeps its left edge on resize, so the focal band stays put
        when the details panel opens; only a window narrowed until the band
        nears an edge recenters on it.
        """
        band = getattr(self, "_focal_band_rect", None)
        if band is None:
            return
        view = self.graphics_view
        visible = view.mapToScene(view.viewport().rect()).boundingRect()
        # A band squeezed against an edge counts as gone, too.
        margin = visible.width() * 0.1
        if visible.left() + margin <= band.left() and band.right() <= visible.right() - margin:
            return
        view.centerOn(band.center().x(), visible.center().y())

    def fit_all_contexts(self) -> None:
        layout = self.layout_result
        if layout is None or not layout.rows:
            return
        available = max(240.0, float(self.graphics_view.viewport().width()) - 36.0)
        self._set_horizontal_scale(max(0.05, min(4.0, available / max(1.0, layout.width))))
        self.graphics_view.centerOn((layout.min_x + layout.max_x) / 2.0, layout.rows[0].center_y)
        self.graphics_view.verticalScrollBar().setValue(self.graphics_view.verticalScrollBar().minimum())
        self._viewport_changed()

    def presence_text(self, gene: GeneLayout) -> str:
        """'Found in 12 of 14 genomes with neighbors' for the gene's family."""
        family = gene.gene.family_key
        if not family:
            return "No matching genes in other rows"
        carriers = self._family_samples.get(family, set())
        shown = " shown" if self.rows_filtered else ""
        if gene.gene.is_focal or not self._has_focal_genes:
            return f"Found in {len(carriers)} of {len(self._all_genomes)} genomes{shown}"
        # Only genomes with genes beside the focal gene can carry a neighbor.
        found = len(carriers & self._genomes_with_neighbors)
        return f"Found in {found} of {len(self._genomes_with_neighbors)} genomes with neighbors{shown}"

    def set_empty_text(self, text: str) -> None:
        """What the plot says when no row is left ('' for the default)."""
        self.empty_text = str(text or "")

    def set_rows_filtered(self, filtered: bool) -> None:
        """Row filters hide genomes: presence counts then name the rows shown."""
        self.rows_filtered = bool(filtered)

    def _gene_tooltip(self, gene: GeneLayout, row: RowLayout) -> str:
        data = gene.gene.payload
        product = clean_text(data.get("Product"))
        lines = [gene.gene.label or "Unnamed gene"]
        if product and product != gene.gene.label:
            lines.append(product)
        if gene.gene.mutation_note:
            lines.append(gene.gene.mutation_note)
        lines.append(_gene_position_text(gene, row.row))
        lines.append(self.presence_text(gene))
        lines.append("Double-click for the full gene card")
        return "\n".join(lines)

    def _set_hover(
        self,
        gene: Optional[GeneLayout],
        row: Optional[RowLayout],
        *,
        refresh: bool = True,
    ) -> None:
        family = (gene.gene.family_key or None) if gene is not None and not gene.gene.is_focal else None
        row_index = row.row_index if row is not None else (gene.row_index if gene is not None else None)
        if gene is self._hovered_gene and family == self.hovered_family and row_index == self.hovered_row_index:
            return
        previous_family = self.hovered_family
        self._hovered_gene = gene
        self.hovered_family = family
        self.hovered_row_index = row_index
        if refresh:
            self._refresh_interaction_items(previous_family, family)
        self.row_header.update()
        hovered_row = row
        if hovered_row is None and gene is not None and self.layout_result is not None:
            hovered_row = self.layout_result.rows[gene.row_index]
        self._emit_hover_samples(hovered_row.row.samples if hovered_row is not None else None)

    def _emit_hover_samples(self, samples: Optional[Sequence[str]]) -> None:
        key = tuple(str(sample).strip() for sample in (samples or []) if str(sample).strip())
        normalized: Optional[Tuple[str, ...]] = key or None
        if normalized == self._last_hover_samples:
            return
        self._last_hover_samples = normalized
        self.contextSamplesHovered.emit(list(normalized) if normalized else None)

    def _handle_mouse_move(self, event: QtGui.QMouseEvent) -> None:
        layout = self.layout_result
        if layout is None:
            return
        position = QtCore.QPointF(event.position())
        self._pointer = ("plot", position)
        point = self.graphics_view.mapToScene(position.toPoint())
        gene = self._gene_at_scene(point)
        # A gene under the pointer wins over the key card drawn above it.
        if gene is None and self._legend_rect.contains(position):
            self._set_hover(None, None)
            self.graphics_view.viewport().setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
            folded = self.legend_collapsed or self._legend_auto_folded
            tooltip = "Click to show the key" if folded else "Click to fold the key"
            if self._not_installed_tooltip:
                tooltip = f"{self._not_installed_tooltip}\n{tooltip}"
            self.graphics_view.setToolTip(tooltip)
            return
        row = layout.row_at(point.y())
        self._set_hover(gene, row)
        on_fold = gene is None and self._fold_at_scene_y(point.y())
        self.graphics_view.viewport().setCursor(
            QtCore.Qt.CursorShape.PointingHandCursor
            if gene is not None or on_fold
            else QtCore.Qt.CursorShape.ArrowCursor
        )
        if gene is not None and row is not None:
            tooltip = self._gene_tooltip(gene, row)
        elif on_fold:
            tooltip = self.fold_tooltip()
        else:
            divider = self._divider_at_scene_y(point.y()) if gene is None else None
            tooltip = self.divider_tooltip(divider) if divider is not None else ""
        self.graphics_view.setToolTip(tooltip)

    def _link_at(self, point: QtCore.QPointF) -> Optional[LinkLayout]:
        for item in reversed(self._pair_items):
            link = item.link_at(point)
            if link is not None:
                return link
        return None

    def _gene_at_scene(self, point: QtCore.QPointF) -> Optional[GeneLayout]:
        layout = self.layout_result
        if layout is None:
            return None
        transform = self.graphics_view.transform()
        scale_x = max(abs(float(transform.m11())), 1e-6)
        scale_y = max(abs(float(transform.m22())), 1e-6)
        return layout.gene_at(
            point.x(),
            point.y(),
            x_tolerance=4.0 / scale_x,
            y_tolerance=3.0 / scale_y,
        )

    def _emit_selected_gene(self, gene: GeneLayout) -> None:
        payload = dict(gene.gene.payload)
        payload["__context_label__"] = gene.gene.label
        payload["__context_family__"] = gene.gene.family_key
        payload["__context_category__"] = gene.gene.category
        payload["__context_start__"] = gene.gene.start
        payload["__context_end__"] = gene.gene.end
        payload["__context_strand__"] = gene.gene.strand
        gene_type = CATEGORY_LABELS.get(gene.gene.category or "OTHER", "Other")
        if gene.gene.mutation_note:
            mutations = clean_text(gene.gene.payload.get("ResistanceMutations")).replace(";", ", ")
            gene_type += f" (resistance mutation{': ' + mutations if mutations else ''})"
        payload["__context_type__"] = gene_type
        payload["__context_found_in__"] = self.presence_text(gene)
        if self.layout_result is not None and 0 <= gene.row_index < len(self.layout_result.rows):
            row = self.layout_result.rows[gene.row_index].row
            payload["__context_place__"] = _row_place_text(row)
            payload["__context_sample__"] = row.sample
            payload["__context_position__"] = _gene_position_text(gene, row, with_strand=False)
        self.geneSelected.emit(payload)

    def _select_gene(self, gene: GeneLayout) -> None:
        previous_family = self.selected_family
        self._selected_gene = gene
        self.selected_family = gene.gene.family_key or None
        self._refresh_interaction_items(previous_family, self.selected_family)
        self._emit_selected_gene(gene)

    def ensure_selection_visible(self) -> None:
        gene = self._selected_gene
        if gene is None or self.layout_result is None:
            return
        rect = QtCore.QRectF(gene.left, gene.center_y - 8.0, max(gene.width, 1.0), 16.0)
        self.graphics_view.ensureVisible(rect, 40, 20)

    def has_selection(self) -> bool:
        return self.selected_family is not None or self._selected_gene is not None

    def clear_selection(self) -> None:
        if self.selected_family is None and self._selected_gene is None:
            return
        previous_family = self.selected_family
        self.selected_family = None
        self._selected_gene = None
        self._refresh_interaction_items(previous_family)
        self.geneSelected.emit({})

    def _handle_mouse_press(self, event: QtGui.QMouseEvent) -> bool:
        if event.button() != QtCore.Qt.MouseButton.LeftButton or self.layout_result is None:
            return False
        point = self.graphics_view.mapToScene(event.position().toPoint())
        gene = self._gene_at_scene(point)
        if gene is None and self._legend_rect.contains(QtCore.QPointF(event.position())):
            self._toggle_legend()
            return True
        if gene is None and self._fold_at_scene_y(point.y()):
            self.toggle_short_rows()
            return True
        if gene is not None:
            if self._selected_gene == gene:
                self.clear_selection()
            else:
                self._select_gene(gene)
            return True
        link = self._link_at(point)
        if link is not None:
            # A ribbon stands for its gene family: clicking it selects that
            # family just like clicking either of the genes it joins.
            already = (
                link.family_key == self.selected_family
                if link.family_key
                else self._selected_gene == link.upper_gene
            )
            if already:
                self.clear_selection()
            else:
                self._select_gene(link.upper_gene)
            return True
        self.clear_selection()
        return False

    def _handle_mouse_double_click(self, event: QtGui.QMouseEvent) -> bool:
        if event.button() != QtCore.Qt.MouseButton.LeftButton or self.layout_result is None:
            return False
        point = self.graphics_view.mapToScene(event.position().toPoint())
        gene = self._gene_at_scene(point)
        if gene is None:
            return self._legend_rect.contains(QtCore.QPointF(event.position()))
        self.geneDoubleClicked.emit(dict(gene.gene.payload))
        return True

    def navigate_selection(self, dx: int, dy: int) -> bool:
        layout = self.layout_result
        if layout is None or not layout.rows:
            return False
        current = self._selected_gene
        if current is None:
            row = layout.rows[0]
            current = next((gene for gene in row.genes if gene.gene.is_focal), row.genes[0] if row.genes else None)
            if current is None:
                return False
        target = current
        if dx:
            row = layout.rows[current.row_index]
            index = row.genes.index(current)
            target = row.genes[max(0, min(len(row.genes) - 1, index + dx))]
        elif dy:
            row_index = max(0, min(len(layout.rows) - 1, current.row_index + dy))
            row = layout.rows[row_index]
            same_family = [gene for gene in row.genes if gene.gene.family_key == current.gene.family_key]
            target = min(same_family or list(row.genes), key=lambda gene: abs(gene.center_x - current.center_x)) if row.genes else current
        self._select_gene(target)
        self.graphics_view.centerOn(target.center_x, target.center_y)
        return True

    def open_selected_gene(self) -> bool:
        if self._selected_gene is None:
            return False
        self.geneDoubleClicked.emit(dict(self._selected_gene.gene.payload))
        return True

    def _color_legend_entries(self) -> List[Tuple[str, QtGui.QColor]]:
        """Gene type key entries for the types present in the rows."""
        if self.color_mode != "annotation":
            return []
        present = {
            gene.category
            for row in self._model.rows
            for gene in row.genes
            if gene.category
        }
        entries = [
            (CATEGORY_LABELS[category], QtGui.QColor(CATEGORY_COLORS[category]))
            for category in ("AMR", "VIRULENCE", "MOBILE", "BMRG")
            if category in present
        ]
        has_other = any(
            not gene.category
            for row in self._model.rows
            for gene in row.genes
        )
        if has_other:
            entries.append((CATEGORY_LABELS["OTHER"], QtGui.QColor(CATEGORY_COLORS["OTHER"])))
        return entries

    def _legend_groups(self, scale_text: str) -> List[List[Tuple[str, str, Optional[QtGui.QColor]]]]:
        """Key content as (swatch kind, text, color) tokens, one group per line."""
        groups: List[List[Tuple[str, str, Optional[QtGui.QColor]]]] = []
        if self.color_mode == "annotation":
            groups.append([("swatch", label, color) for label, color in self._color_legend_entries()])
            if self._not_installed_text and not self.exporting:  # screen only, not figures
                groups.append([("note", self._not_installed_text, None)])
        else:
            groups.append([("families", "Same color: same gene family", None)])
        if self._badge_layer is not None:
            groups.append([("mutation", MUTATION_LEGEND_TEXT, None)])
        if self.link_mode != "off" and self.show_links:
            groups.append([
                ("link", "same direction", QtGui.QColor(*LINK_RGB, 150)),
                ("link", "inverted", QtGui.QColor(*INVERTED_LINK_RGB, 175)),
            ])
        if self._has_focal_genes:
            focal = next(
                (gene for row in (self.layout_result.rows if self.layout_result else ()) for gene in row.genes if gene.gene.is_focal),
                None,
            )
            fill = self.gene_fill(focal) if focal is not None else QtGui.QColor(CATEGORY_COLORS["OTHER"])
            text = f"focal gene ({self.focal_label})" if self.focal_label else "focal gene"
            groups.append([("focal", text, fill)])
        if self.layout_result is not None and any(
            row.cap_left or row.cap_right for row in self.layout_result.rows
        ):
            groups.append([("contig_end", "contig ends here", None)])
        groups.append([("scale", scale_text, None)])
        return [group for group in groups if group]

    @staticmethod
    def _legend_font() -> QtGui.QFont:
        font = QtGui.QFont()
        font.setPixelSize(11)
        return font

    def _paint_legend(
        self,
        painter: Optional[QtGui.QPainter],
        origin: QtCore.QPointF,
        *,
        max_width: float,
        scale_width: float,
        scale_text: str,
        card: bool,
    ) -> QtCore.QSizeF:
        """Lay out (and, with a painter, draw) the key; returns its size.

        As a card each group gets its own line under a 'Key' title; as a
        figure strip the groups run on one line and wrap only when needed.
        """
        font = self._legend_font()
        metrics = QtGui.QFontMetricsF(font)
        line_height = max(16.0, metrics.height() + 4.0)
        swatch_widths = {"swatch": 10.0, "families": 25.0, "link": 16.0, "focal": 20.0, "contig_end": 16.0, "mutation": 12.0,
                         "note": 0.0}
        x0, y0 = origin.x(), origin.y()
        x, y = x0, y0
        used_width = 0.0
        if painter is not None:
            painter.save()
            painter.setFont(font)
        if card:
            if painter is not None:
                title_font = QtGui.QFont(font)
                title_font.setBold(True)
                painter.setFont(title_font)
                painter.setPen(MUTED_TEXT_COLOR)
                painter.drawText(QtCore.QPointF(x0, y0 + line_height - 5.0), "Key")
                painter.setFont(font)
            y += line_height
        for group_index, group in enumerate(self._legend_groups(scale_text)):
            if card and group_index:
                x = x0
                y += line_height
            elif group_index:
                x += 14.0
            for kind, text, color in group:
                swatch = scale_width if kind == "scale" else swatch_widths.get(kind, 10.0)
                width = swatch + 5.0 + metrics.horizontalAdvance(text)
                if x > x0 and x + width > x0 + max_width:
                    x = x0
                    y += line_height
                if painter is not None:
                    self._paint_legend_swatch(painter, kind, color, QtCore.QRectF(x, y, swatch, line_height))
                    painter.setPen(MUTED_TEXT_COLOR if kind in ("scale", "note") else TEXT_COLOR)
                    painter.drawText(QtCore.QPointF(x + swatch + 5.0, y + line_height - 5.0), text)
                x += width + 10.0
                used_width = max(used_width, x - x0 - 10.0)
        if painter is not None:
            painter.restore()
        return QtCore.QSizeF(used_width, y + line_height - y0)

    def _paint_legend_swatch(
        self,
        painter: QtGui.QPainter,
        kind: str,
        color: Optional[QtGui.QColor],
        rect: QtCore.QRectF,
    ) -> None:
        middle = rect.center().y()
        if kind == "swatch" and color is not None:
            painter.setPen(QtGui.QPen(QtGui.QColor(color).darker(125), 0.8))
            painter.setBrush(color)
            painter.drawRoundedRect(QtCore.QRectF(rect.left(), middle - 5.0, 10.0, 10.0), 2.0, 2.0)
        elif kind == "mutation":
            paint_mutation_badge(painter, QtCore.QRectF(rect.left(), middle - 6.0, 12.0, 12.0))
        elif kind == "families":
            painter.setPen(QtCore.Qt.PenStyle.NoPen)
            for index, value in enumerate(FAMILY_PALETTE[:3]):
                painter.setBrush(QtGui.QColor(value))
                painter.drawRoundedRect(QtCore.QRectF(rect.left() + index * 9.0, middle - 4.0, 7.0, 8.0), 1.5, 1.5)
        elif kind == "link" and color is not None:
            ribbon = QtGui.QPolygonF([
                QtCore.QPointF(rect.left(), middle - 5.0),
                QtCore.QPointF(rect.left() + 10.0, middle - 5.0),
                QtCore.QPointF(rect.left() + 16.0, middle + 5.0),
                QtCore.QPointF(rect.left() + 6.0, middle + 5.0),
            ])
            painter.setPen(QtCore.Qt.PenStyle.NoPen)
            painter.setBrush(color)
            painter.drawPolygon(ribbon)
        elif kind == "focal" and color is not None:
            painter.setPen(QtCore.Qt.PenStyle.NoPen)
            painter.setBrush(FOCAL_BAND_COLOR)
            painter.drawRect(QtCore.QRectF(rect.left() + 3.0, middle - 7.0, 14.0, 14.0))
            painter.setPen(QtGui.QPen(FOCAL_COLOR, 1.4))
            painter.setBrush(color)
            painter.drawPath(gene_arrow_path(QtCore.QRectF(rect.left(), middle - 4.0, 20.0, 8.0), "+", 1.0))
        elif kind == "contig_end":
            painter.setPen(QtGui.QPen(BASELINE_COLOR, 1.0))
            painter.drawLine(QtCore.QPointF(rect.left(), middle), QtCore.QPointF(rect.left() + 11.0, middle))
            _draw_contig_end(painter, rect.left() + 13.0, middle, 8.0, 1.0, -1.0)
        elif kind == "scale":
            painter.setPen(QtGui.QPen(MUTED_TEXT_COLOR, 1.2))
            painter.drawLine(QtCore.QPointF(rect.left(), middle), QtCore.QPointF(rect.right(), middle))
            painter.drawLine(QtCore.QPointF(rect.left(), middle - 4.0), QtCore.QPointF(rect.left(), middle + 4.0))
            painter.drawLine(QtCore.QPointF(rect.right(), middle - 4.0), QtCore.QPointF(rect.right(), middle + 4.0))

    def _draw_view_overlay(self, painter: QtGui.QPainter) -> None:
        layout = self.layout_result
        if layout is None or not layout.rows:
            self._legend_rect = QtCore.QRectF()
            return
        scale_x = max(abs(float(self.graphics_view.transform().m11())), 1e-6)
        scale_bp = _nice_scale_bp(100.0 / (layout.px_per_bp * scale_x))
        scale_width = scale_bp * layout.px_per_bp * scale_x
        viewport = self.graphics_view.viewport()
        painter.save()
        painter.resetTransform()
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        self._paint_fold_line(painter)
        right = float(viewport.width()) - 12.0
        bottom = float(viewport.height()) - 12.0
        # Map legend card style: white, rounded, muted text, soft shadow.
        rect: Optional[QtCore.QRectF] = None
        self._legend_auto_folded = False
        if not self.legend_collapsed:
            padding_x, padding_y = 11.0, 8.0
            size = self._paint_legend(
                None,
                QtCore.QPointF(0.0, 0.0),
                max_width=300.0,
                scale_width=scale_width,
                scale_text=_scale_label(scale_bp),
                card=True,
            )
            rect = self._legend_card_rect(
                size.width() + 2.0 * padding_x,
                size.height() + 2.0 * padding_y,
            )
            if rect is None:
                # Every corner would cover genes: fold to the pill instead.
                self._legend_auto_folded = True
        if rect is None:
            font = self._legend_font()
            font.setBold(True)
            width = QtGui.QFontMetricsF(font).horizontalAdvance("Key") + 22.0
            rect = QtCore.QRectF(right - width, bottom - 22.0, width, 22.0)
        folded = self.legend_collapsed or self._legend_auto_folded
        shadow = QtGui.QColor(COLORS["text_strong"])
        shadow.setAlpha(28)
        painter.setPen(QtCore.Qt.PenStyle.NoPen)
        painter.setBrush(shadow)
        painter.drawRoundedRect(rect.translated(0.0, 1.5), 6.0, 6.0)
        plate = QtGui.QColor(COLORS["surface"])
        plate.setAlpha(235)
        painter.setPen(QtGui.QPen(SEPARATOR_COLOR, 1.0))
        painter.setBrush(plate)
        painter.drawRoundedRect(rect, 6.0, 6.0)
        if folded:
            font = self._legend_font()
            font.setBold(True)
            painter.setFont(font)
            painter.setPen(MUTED_TEXT_COLOR)
            painter.drawText(rect, QtCore.Qt.AlignmentFlag.AlignCenter, "Key")
        else:
            self._paint_legend(
                painter,
                QtCore.QPointF(rect.left() + 11.0, rect.top() + 8.0),
                max_width=300.0,
                scale_width=scale_width,
                scale_text=_scale_label(scale_bp),
                card=True,
            )
        painter.restore()
        self._legend_rect = rect

    def _paint_fold_line(self, painter: QtGui.QPainter) -> None:
        """The short-contig sentence, pinned to the plot's left edge."""
        self._fold_link_rect = QtCore.QRectF()
        fold = self.layout_result.fold if self.layout_result is not None else None
        if fold is None:
            return
        view = self.graphics_view
        center = float(view.mapFromScene(QtCore.QPointF(0.0, fold.center_y)).y())
        rect = QtCore.QRectF(12.0, center - fold.height / 2.0, float(view.viewport().width()) - 24.0, fold.height)
        if rect.bottom() < 0 or rect.top() > view.viewport().height():
            return
        self._fold_link_rect = _paint_band_label(
            painter,
            rect,
            self.fold_sentence() + ".",
            link="Hide" if fold.expanded else "Show",
        )

    def _visible_gene_rects(self) -> List[QtCore.QRectF]:
        """Viewport rectangles of the genes on screen, labels included."""
        layout = self.layout_result
        if layout is None:
            return []
        view = self.graphics_view
        visible = view.mapToScene(view.viewport().rect()).boundingRect()
        rects: List[QtCore.QRectF] = []
        for row in layout.rows:
            top = row.center_y - row.gene_height / 2.0
            if row.center_y + row.gene_height / 2.0 < visible.top() or top - 14.0 > visible.bottom():
                continue
            for gene in row.genes:
                if gene.right < visible.left() or gene.left > visible.right():
                    continue
                scene_rect = QtCore.QRectF(gene.left, top, max(gene.width, 0.01), row.gene_height)
                # A little headroom covers the name drawn above small genes.
                rects.append(view.mapFromScene(scene_rect).boundingRect().toRectF().adjusted(0.0, -12.0, 0.0, 0.0))
        return rects

    def _legend_card_rect(self, width: float, height: float) -> Optional[QtCore.QRectF]:
        """First plot corner where the open key covers no gene or fold line, or None.

        Bottom right first, then the other corners.  A key the viewer opened
        on purpose stays bottom right even over genes.
        """
        viewport = self.graphics_view.viewport()
        margin = 12.0
        left = margin
        right = float(viewport.width()) - margin - width
        top = margin
        bottom = float(viewport.height()) - margin - height
        corners = [
            QtCore.QRectF(right, bottom, width, height),
            QtCore.QRectF(right, top, width, height),
            QtCore.QRectF(left, bottom, width, height),
            QtCore.QRectF(left, top, width, height),
        ]
        if self._legend_forced_open:
            return corners[0]
        if left > right or top > bottom:
            return None
        genes = self._visible_gene_rects()
        if not self._fold_link_rect.isEmpty():
            # The no-neighbors sentence counts as content too; it was painted
            # just before the key, from the left edge to its link.
            link = self._fold_link_rect
            genes.append(QtCore.QRectF(left, link.top(), link.right() - left + 4.0, link.height()))
        for corner in corners:
            if not any(corner.intersects(gene) for gene in genes):
                return corner
        return None

    def _export_scale(self, x_scale: float) -> Tuple[float, float]:
        layout = self.layout_result
        px_per_bp = layout.px_per_bp if layout is not None else 0.065
        scale_bp = _nice_scale_bp(130.0 / max(px_per_bp * x_scale, 1e-6))
        return scale_bp, scale_bp * px_per_bp * x_scale

    def _document_size(self) -> Tuple[float, float, float]:
        layout = self.layout_result
        if layout is None:
            return 1200.0, 500.0, 1.0
        header_width = float(RAIL_WIDTH)
        plot_width = max(900.0, min(3000.0, layout.width))
        x_scale = plot_width / max(1.0, layout.width)
        width = header_width + plot_width + 48.0
        scale_bp, scale_width = self._export_scale(x_scale)
        # The key sits under the plot, never on top of genes.
        legend = self._paint_legend(
            None,
            QtCore.QPointF(0.0, 0.0),
            max_width=width - 36.0,
            scale_width=scale_width,
            scale_text=_scale_label(scale_bp),
            card=False,
        )
        footnote = 18.0 if self._omitted_short_note() else 0.0
        return width, 42.0 + layout.height + 14.0 + legend.height() + footnote + 16.0, x_scale

    def _omitted_short_note(self) -> str:
        """Figure footnote naming the short contigs a figure leaves out."""
        layout = self.layout_result
        if layout is None or self._short_rows_mode() != "omitted":
            return ""
        shown = {id(row.row) for row in layout.rows}
        omitted = tuple(row for row in self._model.rows if row.is_short and id(row) not in shown)
        if not omitted:
            return ""
        fold = FoldLayout(center_y=0.0, height=0.0, rows=omitted, expanded=False)
        return _fold_sentence(fold, self.focal_label) + " (not shown)."

    @contextlib.contextmanager
    def _figure_layout(self) -> Iterator[None]:
        """Lay the rows out as Export > Include short contigs asks, then
        restore the screen layout and view."""
        wanted = "expanded" if self.export_include_short_contigs else "omitted"
        has_short = any(row.is_short for row in self._model.rows)
        if not has_short or self._short_rows_mode() == wanted:
            yield
            return
        view = self.graphics_view
        transform = view.transform()
        scroll = (view.horizontalScrollBar().value(), view.verticalScrollBar().value())
        # The figure layout may leave out the row of the pinned gene: keep
        # the pin and keep the dialog's details panel out of it.
        pinned = self._gene_identity(self._selected_gene, self.layout_result)
        was_blocked = self.blockSignals(True)
        self._short_rows_mode_override = wanted
        try:
            self._rebuild_layout()
            yield
        finally:
            self._short_rows_mode_override = None
            self._rebuild_layout(restore_identity=pinned)
            view.setTransform(transform)
            view.horizontalScrollBar().setValue(scroll[0])
            view.verticalScrollBar().setValue(scroll[1])
            self.blockSignals(was_blocked)
            self._viewport_changed()

    def _paint_document(self, painter: QtGui.QPainter) -> None:
        layout = self.layout_result
        if layout is None:
            return
        # Publication style is a figure option: it applies only while the
        # figure is painted and leaves the screen as it was.
        screen_mode = self._display_mode
        if self.export_publication_style:
            self._display_mode = "figure"
        try:
            self._paint_document_content(painter, layout)
        finally:
            self._display_mode = screen_mode
            for index, item in enumerate(self._row_background_items):
                item.setBrush(QtGui.QBrush(self._row_background_color(index)))

    def _paint_document_content(self, painter: QtGui.QPainter, layout: ContextLayout) -> None:
        width, height, x_scale = self._document_size()
        header_width = float(RAIL_WIDTH)
        top = 42.0
        publication = self._display_mode in {"figure", "publication"}
        document_background = PUBLICATION_BACKGROUND if publication else CANVAS_BACKGROUND
        painter.fillRect(QtCore.QRectF(0.0, 0.0, width, height), document_background)
        for index, item in enumerate(self._row_background_items):
            item.setBrush(QtGui.QBrush(self._row_background_color(index, publication=publication)))
        title_font = painter.font()
        title_font.setPixelSize(15)
        title_font.setBold(True)
        painter.setFont(title_font)
        painter.setPen(STRONG_TEXT_COLOR)
        painter.drawText(QtCore.QPointF(18.0, 26.0), self.export_title)

        self.exporting = True
        cached_items = self._row_items + self._pair_items
        cache_modes = [item.cacheMode() for item in cached_items]
        for item in cached_items:
            # Screen caches make dense interaction fast; exports deliberately
            # bypass them so PDF paths remain vector and high-DPI PNG is crisp.
            item.setCacheMode(QtWidgets.QGraphicsItem.CacheMode.NoCache)
            item.update()
        source = QtCore.QRectF(layout.min_x, layout.min_y, layout.width, layout.height)
        target = QtCore.QRectF(header_width + 18.0, top, layout.width * x_scale, layout.height)
        try:
            self._scene.render(
                painter,
                target,
                source,
                QtCore.Qt.AspectRatioMode.IgnoreAspectRatio,
            )
        finally:
            self.exporting = False
            for item, cache_mode in zip(cached_items, cache_modes):
                item.setCacheMode(cache_mode)
                item.update()

        painter.fillRect(QtCore.QRectF(0.0, top, header_width, layout.height), CANVAS_BACKGROUND)
        for row in layout.rows:
            center = top + row.center_y
            if row.row_index % 2 and self._display_mode == "exploratory":
                painter.fillRect(QtCore.QRectF(0.0, center - row.pitch / 2.0, header_width, row.pitch), ROW_ALT_BACKGROUND)
            _paint_rail_row(painter, row, center, header_width, show_molecule=self.show_molecule_badges)
        for divider in layout.dividers:
            rect = QtCore.QRectF(0.0, top + divider.center_y - divider.height / 2.0, header_width, divider.height)
            self._paint_divider_band(painter, rect, divider.text)
        if layout.fold is not None:
            fold = layout.fold
            band = QtCore.QRectF(0.0, top + fold.center_y - fold.height / 2.0, header_width, fold.height)
            painter.fillRect(band, GAP_BACKGROUND)
            # The sentence runs over the plot band, without the Show link.
            _paint_band_label(
                painter,
                QtCore.QRectF(header_width + 26.0, band.top(), width - header_width - 44.0, band.height()),
                self.fold_sentence() + ".",
            )
        painter.setPen(QtGui.QPen(BASELINE_COLOR, 1.0))
        painter.drawLine(QtCore.QPointF(header_width, top), QtCore.QPointF(header_width, top + layout.height))

        scale_bp, scale_width = self._export_scale(x_scale)
        legend = self._paint_legend(
            painter,
            QtCore.QPointF(18.0, top + layout.height + 14.0),
            max_width=width - 36.0,
            scale_width=scale_width,
            scale_text=_scale_label(scale_bp),
            card=False,
        )
        note = self._omitted_short_note()
        if note:
            _paint_band_label(
                painter,
                QtCore.QRectF(18.0, top + layout.height + 14.0 + legend.height() + 2.0, width - 36.0, 16.0),
                note,
            )

    def export_png(self, path: Path, *, scale: float = 2.0) -> bool:
        if self.layout_result is None or not self.layout_result.rows:
            return False
        with self._figure_layout():
            return self._export_png(path, scale=scale)

    def _export_png(self, path: Path, *, scale: float) -> bool:
        width, height, _x_scale = self._document_size()
        requested_scale = max(0.5, float(scale))
        max_pixels = 80_000_000.0
        safe_scale = min(requested_scale, math.sqrt(max_pixels / max(1.0, width * height)))
        pixel_width = max(1, int(math.ceil(width * safe_scale)))
        pixel_height = max(1, int(math.ceil(height * safe_scale)))
        image = QtGui.QImage(pixel_width, pixel_height, QtGui.QImage.Format.Format_ARGB32_Premultiplied)
        if image.isNull():
            return False
        image.fill(CANVAS_BACKGROUND)
        painter = QtGui.QPainter(image)
        if not painter.isActive():
            return False
        try:
            painter.setRenderHints(QtGui.QPainter.RenderHint.Antialiasing | QtGui.QPainter.RenderHint.TextAntialiasing)
            self._paint_device_scale = safe_scale
            painter.scale(safe_scale, safe_scale)
            self._paint_document(painter)
        finally:
            self._paint_device_scale = 1.0
            painter.end()
        return bool(image.save(str(path)))

    def export_pdf(self, path: Path) -> bool:
        if self.layout_result is None or not self.layout_result.rows:
            return False
        with self._figure_layout():
            return self._export_pdf(path)

    def _export_pdf(self, path: Path) -> bool:
        width, height, _x_scale = self._document_size()
        writer = QtGui.QPdfWriter(str(path))
        # At 72 DPI one logical document unit maps to one PDF point.  Using 96
        # here would reserve a page 4/3 larger than the painted content.
        writer.setResolution(72)
        writer.setPageSize(
            QtGui.QPageSize(
                QtCore.QSizeF(width, height),
                QtGui.QPageSize.Unit.Point,
                "BACTolomew gene neighborhoods",
            )
        )
        writer.setPageMargins(QtCore.QMarginsF(), QtGui.QPageLayout.Unit.Point)
        painter = QtGui.QPainter(writer)
        if not painter.isActive():
            return False
        try:
            painter.setRenderHints(QtGui.QPainter.RenderHint.Antialiasing | QtGui.QPainter.RenderHint.TextAntialiasing)
            self._paint_device_scale = 1.0
            self._paint_document(painter)
        finally:
            painter.end()
        return bool(Path(path).exists() and Path(path).stat().st_size > 0)
