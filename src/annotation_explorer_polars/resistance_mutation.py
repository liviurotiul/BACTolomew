"""Resistance mutations on genes that are not resistance genes themselves.

AMRFinderPlus point mutations (the ``amrfinder_mutations`` step) give a CDS
the tag ``X`` and the columns ``ResistanceMutations`` (';'-joined element
symbols, e.g. ``gyrA_S81L``) and ``ResistanceMutationClasses`` (their drug
classes). Such a gene is drawn in the AMR red, plus a small badge
(``assets/mutation_badge.svg``) with the sentence of :func:`mutation_note`.
A gene that is a resistance gene and also carries a mutation keeps the AMR
color and gets the badge too. Without the columns only the tag counts.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, List, Mapping, Optional


MUTATION_TAG = "X"
MUTATION_BADGE_SVG = Path(__file__).resolve().parent / "assets" / "mutation_badge.svg"
MUTATION_LEGEND_TEXT = "resistance mutation"

# The badge's QSvgRenderer, made on first paint (this module stays importable
# without Qt, like context_model, which uses it).
_renderer: Any = None


def clean_text(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"", "-", ".", "na", "n/a", "nan", "none", "null", "<na>"} else text


def _split(value: Any) -> List[str]:
    parts = (clean_text(part) for part in clean_text(value).split(";"))
    return list(dict.fromkeys(part for part in parts if part))


def resistance_mutations(row: Optional[Mapping[str, Any]]) -> List[str]:
    """The AMRFinderPlus element symbols of the row (``['gyrA_S81L']``)."""
    return _split(row.get("ResistanceMutations")) if row else []


def has_resistance_mutation(row: Optional[Mapping[str, Any]]) -> bool:
    """The row carries a known resistance mutation (tag X or the column)."""
    if not row:
        return False
    return MUTATION_TAG in clean_text(row.get("Tags")).upper() or bool(resistance_mutations(row))


def mutation_note(row: Optional[Mapping[str, Any]]) -> str:
    """The badge's tooltip sentence, or '' when the row has no mutation."""
    if not has_resistance_mutation(row):
        return ""
    assert row is not None
    mutations = resistance_mutations(row)
    classes = [value.lower() for value in _split(row.get("ResistanceMutationClasses"))]
    if mutations:
        what = ", ".join(mutations)
        if classes:
            what += f" ({', '.join(classes)})"
        what += ", known resistance mutations" if len(mutations) > 1 else ", a known resistance mutation"
    else:
        what = "a known resistance mutation"
    source = "AMRFinderPlus"
    # AMRFinderPlus names the gene it assessed in the symbol (pmrB_I164F); when Bakta named
    # the CDS otherwise (baeS), say which gene AMRFinderPlus read it as.
    assessed = list(dict.fromkeys(m.split("_", 1)[0] for m in mutations if "_" in m))
    gene = clean_text(row.get("Gene")).lower()
    if assessed and gene not in {a.lower() for a in assessed}:
        source += f", which reads this CDS as {' / '.join(assessed)}"
    if "A" in clean_text(row.get("Tags")).upper():
        return f"Resistance gene that also carries {what} ({source})"
    return f"Usually not a resistance gene; carries {what} ({source})"


def paint_mutation_badge(painter: Any, rect: Any) -> None:
    """Draw the shared badge into ``rect`` (a QRectF) with a QPainter."""
    global _renderer
    if _renderer is None:
        from PyQt6 import QtSvg

        _renderer = QtSvg.QSvgRenderer(str(MUTATION_BADGE_SVG))
    _renderer.render(painter, rect)
