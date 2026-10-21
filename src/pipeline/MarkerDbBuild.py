"""Build the AMR marker FASTAs in DB/seq/markers from their unfiltered sources.

One-time DB build step, not part of a run. The marker databases are searched with DIAMOND for
homology, which can only call *acquired* genes. Target-gene references are therefore removed:

* CARD_AMR: keep protein homolog models (HM), as RGI does for presence calls; drop protein
  variant (VM), knockout (KM) and overexpression (OM) models, which match the wild-type
  gyrA, parC, oprD, adeRS, ... of every genome.
* AMRF_*: drop references whose ReferenceGeneCatalog rows are all POINT, POINT_DISRUPT or
  AMR-SUSCEPTIBLE (point-mutation and susceptible-allele references). An accession that also
  has an acquired-gene row (AMR, BIOCIDE, METAL, ...) is kept; accessions missing from the
  catalog are kept and listed in the report.

The unfiltered FASTA of each database is kept as ``<sources>/<stem>.unfiltered.fasta.gz``; the
first build archives the current marker FASTA there, so rebuilding is idempotent. Each build
writes ``<sources>/reports/<stem>.removed.tsv`` and ``<sources>/marker_build_report.json``.
The per-run ``dmnd/queries/<stem>.dmnd`` files are rebuilt by the pipeline when the FASTA
changes (see DiamondWorkflow.dmnd_is_current).

Usage:
    PYTHONPATH=src python -m pipeline.MarkerDbBuild [--db-root DB] [--dry-run]
"""
from __future__ import annotations

import argparse
import datetime as _dt
import gzip
import hashlib
import json
import os
import shutil
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterator, Optional

import pandas as pd

from pipeline.ResistanceMetadata import (
    AMRF_ACCESSION_COLUMNS,
    AMRF_TARGET_SUBTYPES,
    CARD_HOMOLOG_MODEL_PREFIX,
    load_card_target_aros,
    card_model_prefix,
)

CARD_DB_STEMS = ("CARD_AMR",)
AMRF_DB_STEMS = (
    "AMRF_core-AMR",
    "AMRF_plus-AMR",
    "AMRF_core-BMRG",
    "AMRF_plus-BMRG",
    "AMRF_plus-VIRULENCE",
)
CARD_MODEL_LABELS = {
    "HM": "CARD protein homolog model (HM)",
    "VM": "CARD protein variant model (VM)",
    "KM": "CARD protein knockout model (KM)",
    "OM": "CARD protein overexpression model (OM)",
}
SOURCES_DIRNAME = "marker_sources"

Classifier = Callable[[str], Optional[str]]


@dataclass
class FilterResult:
    db: str
    source: Path
    output: Path
    entries_before: int = 0
    entries_after: int = 0
    dropped: list[tuple[str, str]] = field(default_factory=list)  # (header, reason)
    kept_not_in_catalog: list[str] = field(default_factory=list)
    changed: bool = False

    def dropped_by_reason(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for _, reason in self.dropped:
            counts[reason] = counts.get(reason, 0) + 1
        return dict(sorted(counts.items()))


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def iter_fasta_records(path: Path) -> Iterator[tuple[str, list[str]]]:
    """Yield (header without '>', sequence lines with newline) for each record."""
    header: str | None = None
    lines: list[str] = []
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, lines
                header = line[1:].rstrip("\n")
                lines = []
            elif header is not None:
                lines.append(line if line.endswith("\n") else line + "\n")
    if header is not None:
        yield header, lines


def _accession(header: str) -> str:
    return header.split()[0] if header.split() else ""


def card_classifier(target_aros: frozenset[str] = frozenset()) -> Classifier:
    """Keep CARD protein homolog models only.

    Every record is also checked by ARO against card.json: an HM record whose ARO has only
    non-homolog models in card.json (for example baeS/baeR, overexpression models in CARD 4.0.1)
    is dropped too, so the DB agrees with the runtime guard in ResistanceMetadata.
    """

    def classify(header: str) -> Optional[str]:
        prefix = card_model_prefix(header)
        if prefix and prefix != CARD_HOMOLOG_MODEL_PREFIX:
            return CARD_MODEL_LABELS.get(prefix, f"CARD non-homolog model ({prefix})")
        for part in _accession(header).split("|"):
            if part.startswith("ARO:") and part in target_aros:
                return "CARD non-homolog model in card.json (by ARO)"
        return None

    return classify


def load_amrf_subtypes(catalog: Path) -> dict[str, frozenset[str]]:
    df = pd.read_csv(catalog, sep="\t", dtype=str).fillna("")
    subtypes: dict[str, set[str]] = {}
    for row in df.to_dict(orient="records"):
        subtype = str(row.get("subtype", "")).strip().upper()
        for col in AMRF_ACCESSION_COLUMNS:
            acc = str(row.get(col, "")).strip()
            if acc:
                subtypes.setdefault(acc, set()).add(subtype)
    return {acc: frozenset(v) for acc, v in subtypes.items()}


def amrf_classifier(subtypes_by_accession: dict[str, frozenset[str]]) -> Classifier:
    """Drop references whose catalog rows are all target subtypes. Unknown accessions are kept."""

    def classify(header: str) -> Optional[str]:
        subtypes = subtypes_by_accession.get(_accession(header))
        if subtypes and subtypes <= AMRF_TARGET_SUBTYPES:
            return "AMRF " + "+".join(sorted(subtypes))
        return None

    return classify


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def filter_fasta(
    db: str,
    source: Path,
    output: Path,
    classify: Classifier,
    *,
    known_accessions: Optional[set[str]] = None,
    dry_run: bool = False,
) -> FilterResult:
    """Write the records of ``source`` that ``classify`` keeps to ``output``.

    ``output`` is replaced atomically and only when its content changes, so an unchanged
    database keeps its mtime and its per-run .dmnd stays current.
    """
    result = FilterResult(db=db, source=source, output=output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{output.name}.", dir=str(output.parent))
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as out:
            for header, lines in iter_fasta_records(source):
                result.entries_before += 1
                reason = classify(header)
                if reason:
                    result.dropped.append((header, reason))
                    continue
                if known_accessions is not None and _accession(header) not in known_accessions:
                    result.kept_not_in_catalog.append(header)
                result.entries_after += 1
                out.write(f">{header}\n")
                out.writelines(lines)
        result.changed = not output.exists() or _sha256(tmp) != _sha256(output)
        if result.changed and not dry_run:
            if output.exists():
                shutil.copymode(output, tmp)
            os.replace(tmp, output)
        else:
            tmp.unlink()
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    return result


def archive_source(fasta: Path, archived: Path) -> bool:
    """Gzip ``fasta`` to ``archived`` (deterministic header) unless the archive exists."""
    if archived.exists():
        return False
    archived.parent.mkdir(parents=True, exist_ok=True)
    tmp = archived.with_name(archived.name + ".tmp")
    with fasta.open("rb") as src, open(tmp, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
            shutil.copyfileobj(src, gz)
    os.replace(tmp, archived)
    return True


def write_removed_report(result: FilterResult, reports_dir: Path) -> Path:
    reports_dir.mkdir(parents=True, exist_ok=True)
    path = reports_dir / f"{result.db}.removed.tsv"
    rows = [{"header": h, "accession": _accession(h), "reason": r} for h, r in result.dropped]
    pd.DataFrame(rows, columns=["header", "accession", "reason"]).to_csv(path, sep="\t", index=False)
    return path


def build_marker_dbs(db_root: Path, *, dry_run: bool = False) -> dict:
    markers_dir = db_root / "seq" / "markers"
    sources_dir = db_root / "seq" / SOURCES_DIRNAME
    reports_dir = sources_dir / "reports"
    amrf_catalog = db_root / "metadata" / "AMRF" / "ReferenceGeneCatalog.txt"
    card_json = db_root / "metadata" / "CARD" / "card.json"

    subtypes = load_amrf_subtypes(amrf_catalog) if amrf_catalog.is_file() else {}
    card_classify = card_classifier(load_card_target_aros(card_json))
    amrf_classify = amrf_classifier(subtypes)

    summary: dict = {
        "built_at": _dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "amrf_catalog": str(amrf_catalog),
        "amrf_catalog_version": _read_first_line(db_root / "metadata" / "AMRF" / "version.txt"),
        "card_json": str(card_json),
        "rules": {
            "CARD": "keep protein homolog models (HM); drop variant (VM), knockout (KM), overexpression (OM), and any record whose ARO card.json lists only non-homolog models",
            "AMRF": "drop accessions whose catalog subtypes are all in " + ",".join(sorted(AMRF_TARGET_SUBTYPES)),
        },
        "databases": {},
    }
    for stem in (*CARD_DB_STEMS, *AMRF_DB_STEMS):
        fasta = markers_dir / f"{stem}.fasta"
        source = sources_dir / f"{stem}.unfiltered.fasta.gz"
        if not source.exists():
            if not fasta.exists():
                continue
            if dry_run:
                source = fasta
            else:
                archive_source(fasta, source)
        is_card = stem in CARD_DB_STEMS
        result = filter_fasta(
            stem,
            source,
            fasta,
            card_classify if is_card else amrf_classify,
            known_accessions=None if is_card else set(subtypes),
            dry_run=dry_run,
        )
        report = None if dry_run else write_removed_report(result, reports_dir)
        summary["databases"][stem] = {
            "source": str(source),
            "source_sha256": _sha256(source),
            "output": str(fasta),
            "output_sha256": None if dry_run else _sha256(fasta),
            "entries_before": result.entries_before,
            "entries_after": result.entries_after,
            "dropped": len(result.dropped),
            "dropped_by_reason": result.dropped_by_reason(),
            "kept_not_in_catalog": result.kept_not_in_catalog,
            "filtered": bool(result.dropped),
            "changed_by_this_build": result.changed,
            "removed_report": str(report) if report else None,
        }
    if not dry_run:
        sources_dir.mkdir(parents=True, exist_ok=True)
        (sources_dir / "marker_build_report.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )
    return summary


def _read_first_line(path: Path) -> str | None:
    try:
        return path.read_text(encoding="utf-8").splitlines()[0].strip()
    except Exception:
        return None


def main(argv: list[str] | None = None) -> int:
    default_root = Path(__file__).resolve().parents[2] / "DB"
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--db-root", type=Path, default=default_root)
    parser.add_argument("--dry-run", action="store_true", help="report counts without writing anything")
    args = parser.parse_args(argv)
    summary = build_marker_dbs(args.db_root.resolve(), dry_run=args.dry_run)
    for stem, info in summary["databases"].items():
        print(
            f"{stem}: {info['entries_before']} -> {info['entries_after']} "
            f"(dropped {info['dropped']}: {info['dropped_by_reason']}; written={info['changed_by_this_build']})"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
