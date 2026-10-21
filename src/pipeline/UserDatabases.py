"""The marker and compound databases a user adds (Settings > Databases, `bactolomew db`).

BACTolomew ships no marker or compound databases. A user adds a protein FASTA (plain or
gzipped): it is checked (protein, unique sequence ids), copied (unpacked) into
DB/seq/markers (DB/seq/compound for a reference locus) and listed in the catalog with its
category, so the DIAMOND search after Bakta uses it from the next run. The explorer shows
a user database's hits as "<name> (U)" (`user_label`). Removing one drops its catalog entry and moves its
file into DB/_removed/<date>/ (nothing is deleted). Catalog writes are atomic and locked.

    python -m pipeline.UserDatabases --list
    python -m pipeline.UserDatabases --add FASTA --name NAME --category VIRULENCE [--note TEXT]
    python -m pipeline.UserDatabases --remove NAME
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import sys
import tempfile
from zlib import error as zlib_error
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, BinaryIO, Mapping, Optional, Sequence

from pipeline import DbCatalog
from pipeline.BaktaEvidence import EVIDENCE_SOURCES
from pipeline.ChewbbacaSchemaBuild import _atomic_write_json, _catalog_lock


MARKER_CATEGORIES = ("AMR", "BMRG", "VIRULENCE", "MOBILE", "PHAGE", "PLASMID")
COMPOUND_CATEGORY = "COMPOUND"  # a reference locus (several genes), searched as a compound element
CATEGORIES = (*MARKER_CATEGORIES, COMPOUND_CATEGORY)
CATEGORY_LABELS = {
    "AMR": "AMR",
    "BMRG": "BMRG (stress, metal, biocide)",
    "VIRULENCE": "Virulence",
    "MOBILE": "Mobile element",
    "PHAGE": "Phage",
    "PLASMID": "Plasmid",
    COMPOUND_CATEGORY: "Reference locus (compound)",
}

# Letters, digits, '-' and '_'; DIAMOND hit files are named <genome>-vs-<database>.
NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]{0,63}$")
_PROTEIN_LETTERS = frozenset("ABCDEFGHIJKLMNOPQRSTUVWXYZ*-")  # DIAMOND makedb rejects "."
# A, C, G, T, U, N and the other IUPAC nucleotide codes (several are amino acids too).
_NUCLEOTIDE_LETTERS = frozenset("ACGTUNRYSWKMBDHV")
NUCLEOTIDE_SHARE = 0.9  # at least this share of nucleotide codes: a nucleotide FASTA
_GZIP_MAGIC = b"\x1f\x8b"
_UTF8_BOM = b"\xef\xbb\xbf"

# What the explorer appends to a user database's name wherever its hits are shown.
USER_LABEL_SUFFIX = " (U)"


def user_label(name: str) -> str:
    """'VFDB_Ab30' -> 'VFDB_Ab30 (U)', the name the explorer shows for a user database."""
    return f"{name}{USER_LABEL_SUFFIX}"


def run_user_databases(inputs: Optional[Mapping[str, Any]]) -> frozenset[str]:
    """The marker and compound databases (file stems) of a run that are user databases:
    every one of a run made with the bakta_evidence step (its run_config records
    `inputs.databases.evidence`; BACTolomew ships none since then), none of an older run
    (those searched the built-in catalog)."""
    if not isinstance(inputs, Mapping):
        return frozenset()
    record = inputs.get("databases")
    if not isinstance(record, Mapping) or not isinstance(record.get("evidence"), Mapping):
        return frozenset()
    names: set[str] = set()
    for key in ("marker_sequences", "reference_locus_sequences"):
        for raw in inputs.get(key) or []:
            stem = Path(str(raw)).stem.strip()
            if stem:
                names.add(stem)
    for kind in ("marker", "compound"):
        entries = record.get(kind) if isinstance(record.get(kind), Mapping) else {}
        for db in [*(entries.get("used") or []), *(entries.get("not_installed") or [])]:
            name = str(db.get("name") or "").strip() if isinstance(db, Mapping) else ""
            if name:
                names.add(name)
    return frozenset(names)


def _reserved_name(name: str) -> bool:
    """A name of the pipeline's own evidence sources (AMRFinderPlus, MOB-suite, Bakta ...)."""
    key = re.sub(r"[-_ ]", "", name).lower()
    reserved = {re.sub(r"[-_ ]", "", source).lower() for source in EVIDENCE_SOURCES}
    return key.startswith("bakta") or key in reserved


class DatabaseError(ValueError):
    """A database cannot be added or removed; the message says why."""


@dataclass(frozen=True)
class UserDatabase:
    name: str
    kind: str  # "marker" or "compound"
    category: str
    path: Path
    size: int
    sequences: Optional[int]
    added: str
    note: str
    source_path: str
    sha256: str

    @property
    def ok(self) -> bool:
        return self.path.is_file()

    @property
    def status(self) -> str:
        return "ok" if self.ok else "file missing"


def default_db_root() -> Path:
    return DbCatalog.DB_ROOT


def catalog_path(kind: str, db_root: Optional[Path] = None) -> Path:
    root = Path(db_root) if db_root else default_db_root()
    name = DbCatalog.COMPOUND_DB_CATALOG_PATH.name if kind == "compound" else DbCatalog.MARKER_DB_CATALOG_PATH.name
    return root / "catalogs" / name


def _seq_dir(kind: str, db_root: Path) -> Path:
    return db_root / "seq" / ("compound" if kind == "compound" else "markers")


def human_size(size: int) -> str:
    value = float(size)
    for unit in ("B", "KB", "MB", "GB"):
        if value < 1024 or unit == "GB":
            return f"{value:.0f} {unit}" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024
    return f"{size} B"


# ---------------------------------------------------------------- FASTA checks
def _is_gzip(path: Path) -> bool:
    with Path(path).open("rb") as fh:
        return fh.read(2) == _GZIP_MAGIC


def _open_fasta(path: Path) -> BinaryIO:
    """The FASTA's bytes, unpacked when it is gzipped."""
    return gzip.open(path, "rb") if _is_gzip(path) else Path(path).open("rb")  # type: ignore[return-value]


def _unpacked_size(path: Path) -> int:
    """About the bytes the FASTA takes once unpacked (gzip records it modulo 4 GiB)."""
    size = Path(path).stat().st_size
    if not _is_gzip(path):
        return size
    with Path(path).open("rb") as fh:
        fh.seek(-4, os.SEEK_END)
        recorded = int.from_bytes(fh.read(4), "little")
    return max(recorded, size * 4)


def validate_protein_fasta(path: Path) -> int:
    """The number of sequences of a protein FASTA (plain or gzipped, a UTF-8 BOM is
    ignored); DatabaseError when it is not one (nucleotide, empty sequences, characters
    that are no residues, repeated ids)."""
    path = Path(path).expanduser()
    if not path.is_file():
        raise DatabaseError(f"File not found: {path}")
    try:
        return _validate_protein_lines(path)
    except (gzip.BadGzipFile, EOFError, zlib_error) as exc:
        raise DatabaseError(f"The gzipped file cannot be unpacked ({exc}).") from exc


def _validate_protein_lines(path: Path) -> int:
    ids: set[str] = set()
    current: Optional[str] = None
    current_len = 0
    letters = 0
    nucleotides = 0
    line_no = 0

    def close(line_no: int) -> None:
        if current is not None and current_len == 0:
            raise DatabaseError(f"Sequence '{current}' (line {line_no}) is empty.")

    with _open_fasta(path) as raw_fh, io.TextIOWrapper(raw_fh, encoding="utf-8-sig", errors="replace") as fh:
        for line_no, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                close(line_no)
                seq_id = line[1:].split(None, 1)[0] if line[1:].strip() else ""
                if not seq_id:
                    raise DatabaseError(f"Line {line_no} has a header without a sequence id.")
                if seq_id in ids:
                    raise DatabaseError(f"Sequence id '{seq_id}' appears more than once (line {line_no}).")
                ids.add(seq_id)
                current, current_len = seq_id, 0
                continue
            if current is None:
                raise DatabaseError("Not a FASTA file: the first line must be a '>' header.")
            residues = line.upper().replace(" ", "")
            bad = set(residues) - _PROTEIN_LETTERS
            if bad:
                raise DatabaseError(
                    f"Line {line_no} holds characters that are no amino acids: {''.join(sorted(bad))[:10]}"
                )
            current_len += len(residues)
            letters += len(residues)
            nucleotides += sum(residues.count(c) for c in _NUCLEOTIDE_LETTERS)
        close(line_no)
    if not ids:
        raise DatabaseError("Not a FASTA file: no '>' header found.")
    if letters and nucleotides / letters >= NUCLEOTIDE_SHARE:
        raise DatabaseError(
            "This looks like a nucleotide FASTA. The search after Bakta compares proteins: "
            "add the protein sequences (amino acids) instead."
        )
    return len(ids)


_COUNT_CACHE: dict[tuple[str, int, int], int] = {}


def _count_sequences(path: Path) -> Optional[int]:
    try:
        st = path.stat()
    except OSError:
        return None
    key = (str(path), st.st_size, st.st_mtime_ns)
    if key not in _COUNT_CACHE:
        count = 0
        with path.open("rb") as fh:
            for line in fh:
                if line.startswith(b">"):
                    count += 1
        _COUNT_CACHE[key] = count
    return _COUNT_CACHE[key]


def _stored_count(info: dict[str, Any]) -> Optional[int]:
    """The sequence count the catalog entry recorded when it was added."""
    count = info.get("sequences")
    return count if isinstance(count, int) and not isinstance(count, bool) and count > 0 else None


# ---------------------------------------------------------------- catalogs
def _read_raw(catalog: Path) -> dict[str, Any]:
    if not catalog.is_file():
        return {}
    try:
        data = json.loads(catalog.read_text(encoding="utf-8"))
    except Exception as exc:
        raise DatabaseError(f"Cannot read the catalog {catalog}: {exc}") from exc
    if not isinstance(data, dict):
        raise DatabaseError(f"The catalog {catalog} is not a JSON object.")
    return data


def _entry_path(catalog: Path, info: dict[str, Any]) -> Path:
    return Path(DbCatalog._resolve_catalog_path(catalog, info.get("db_path") or info.get("path") or ""))


def list_databases(db_root: Optional[Path] = None) -> list[UserDatabase]:
    """Every entry of the marker and compound catalogs (an unreadable catalog raises)."""
    root = Path(db_root) if db_root else default_db_root()
    out: list[UserDatabase] = []
    for kind in ("marker", "compound"):
        catalog = catalog_path(kind, root)
        for name, info in _read_raw(catalog).items():
            if not isinstance(info, dict) or not (info.get("db_path") or info.get("path")):
                continue
            path = _entry_path(catalog, info)
            try:
                size = path.stat().st_size
            except OSError:
                size = 0
            category = str(info.get("category") or "") or (COMPOUND_CATEGORY if kind == "compound" else "")
            out.append(UserDatabase(
                name=str(name),
                kind=kind,
                category=category,
                path=path,
                size=size,
                sequences=(_stored_count(info) or _count_sequences(path)) if size else None,
                added=str(info.get("added") or ""),
                note=str(info.get("note") or ""),
                source_path=str(info.get("source_path") or ""),
                sha256=str(info.get("sha256") or ""),
            ))
    return out


def check_name(name: str, db_root: Optional[Path] = None) -> str:
    """``name`` stripped, when it is a safe name no database (or file) uses yet."""
    name = str(name or "").strip()
    if not NAME_RE.match(name) or "-vs-" in name:
        raise DatabaseError(
            "Use a name of letters, digits, '-' and '_' (up to 64, starting with a letter or digit) "
            "that does not contain '-vs-' (the hit files are named <genome>-vs-<database>)."
        )
    if _reserved_name(name):
        raise DatabaseError(
            f"'{name}' is taken by the pipeline's own evidence sources (AMRFinderPlus, MOB-suite, "
            "Bakta ...); choose another name."
        )
    root = Path(db_root) if db_root else default_db_root()
    taken: set[str] = set()
    for db in list_databases(root):
        taken |= {db.name.lower(), db.path.stem.lower()}
    for kind in ("marker", "compound"):
        folder = _seq_dir(kind, root)
        if folder.is_dir():
            taken |= {p.stem.lower() for p in folder.iterdir()}
    if name.lower() in taken:
        raise DatabaseError(f"A database named '{name}' exists already.")
    return name


def add_database(
    fasta: Path,
    name: str,
    category: str,
    *,
    note: str = "",
    db_root: Optional[Path] = None,
) -> UserDatabase:
    """Check ``fasta``, copy it into the DB folder and list it in its catalog."""
    root = Path(db_root) if db_root else default_db_root()
    category = str(category or "").strip().upper()
    if category not in CATEGORIES:
        raise DatabaseError(f"Unknown category '{category}'; use one of {', '.join(CATEGORIES)}.")
    kind = "compound" if category == COMPOUND_CATEGORY else "marker"
    source = Path(fasta).expanduser().resolve()
    sequences = validate_protein_fasta(source)
    catalog = catalog_path(kind, root)
    with _catalog_lock(catalog):
        name = check_name(name, root)
        data = _read_raw(catalog)
        folder = _seq_dir(kind, root)
        folder.mkdir(parents=True, exist_ok=True)
        target = folder / f"{name}.fasta"
        needed = _unpacked_size(source)
        free = shutil.disk_usage(folder).free
        if free < needed:
            raise DatabaseError(
                f"Not enough free space in {folder}: the copy needs {human_size(needed)}, "
                f"{human_size(free)} free."
            )
        digest = hashlib.sha256()
        fd, tmp = tempfile.mkstemp(prefix=f".{target.name}.", suffix=".tmp", dir=str(folder))
        try:
            with os.fdopen(fd, "wb") as out, _open_fasta(source) as src:
                first = True
                for chunk in iter(lambda: src.read(1 << 20), b""):
                    if first and chunk.startswith(_UTF8_BOM):
                        chunk = chunk[len(_UTF8_BOM):]
                    first = False
                    digest.update(chunk)
                    out.write(chunk)
            os.replace(tmp, target)
        except BaseException:
            Path(tmp).unlink(missing_ok=True)
            raise
        entry: dict[str, Any] = {
            "category": "" if kind == "compound" else category,
            "db_path": os.path.relpath(target, catalog.parent),
            "user_added": True,
            "added": datetime.now().isoformat(timespec="seconds"),
            "source_path": str(source),
            "sha256": digest.hexdigest(),
            "sequences": sequences,
        }
        if note.strip():
            entry["note"] = note.strip()
        data[name] = entry
        try:
            _atomic_write_json(catalog, data)
        except BaseException:
            target.unlink(missing_ok=True)
            raise
    return next(db for db in list_databases(root) if db.name == name and db.kind == kind)


def remove_database(name: str, db_root: Optional[Path] = None) -> Optional[Path]:
    """Drop ``name`` from its catalog and move its file into DB/_removed/<date>/; the
    file's new path (None when it was missing or lies outside the DB folder, then it is
    left where it is)."""
    root = Path(db_root) if db_root else default_db_root()
    for kind in ("marker", "compound"):
        catalog = catalog_path(kind, root)
        with _catalog_lock(catalog):
            data = _read_raw(catalog)
            if name not in data:
                continue
            info = data.pop(name)
            path = _entry_path(catalog, info) if isinstance(info, dict) else None
            moved: Optional[Path] = None
            if path is not None and path.is_file() and path.resolve().is_relative_to((root / "seq").resolve()):
                dest_dir = root / "_removed" / datetime.now().strftime("%Y-%m-%d")
                dest_dir.mkdir(parents=True, exist_ok=True)
                moved = dest_dir / path.name
                n = 1
                while moved.exists():
                    moved = dest_dir / f"{path.stem}-{n}{path.suffix}"
                    n += 1
                shutil.move(str(path), str(moved))
            try:
                _atomic_write_json(catalog, data)
            except BaseException:
                if moved is not None:
                    shutil.move(str(moved), str(path))
                raise
            return moved
    raise DatabaseError(f"No database named '{name}'.")


# ---------------------------------------------------------------- CLI
def add_cli_arguments(parser: argparse.ArgumentParser) -> None:
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--list", action="store_true", help="list the marker and compound databases")
    action.add_argument("--add", metavar="FASTA", help="add a protein FASTA (plain or .gz) as a database")
    action.add_argument("--remove", metavar="NAME", help="remove a database (its file moves to DB/_removed)")
    parser.add_argument("--name", help="name of the added database (letters, digits, - _)")
    parser.add_argument("--category", type=str.upper, choices=CATEGORIES,
                        help="category of the added database (COMPOUND: a reference locus)")
    parser.add_argument("--note", default="", help="optional note, for example the database's license")
    parser.add_argument("--db-root", metavar="DIR", help="DB folder (default: the project's DB)")


def run_cli(args: argparse.Namespace) -> int:
    root = Path(args.db_root).expanduser() if args.db_root else None
    try:
        if args.list:
            databases = list_databases(root)
            if not databases:
                print("No databases added. Runs search with Bakta and the pipeline's own tools only.")
            for db in databases:
                count = "" if db.sequences is None else f"{db.sequences} seq"
                print(f"{db.name:<40} {db.category:<14} {count:>12} {human_size(db.size):>10} "
                      f"{db.added[:10]:<10} {db.status}")
            return 0
        if args.add:
            if not args.name or not args.category:
                print("--add needs --name and --category.", file=sys.stderr)
                return 2
            db = add_database(Path(args.add), args.name, args.category, note=args.note, db_root=root)
            print(f"Added {db.name} ({db.category}, {db.sequences} sequences): {db.path}")
            print("It is searched after Bakta from the next run.")
            return 0
        moved = remove_database(args.remove, root)
        print(f"Removed {args.remove}" + (f"; its file is in {moved}" if moved else ""))
        return 0
    except DatabaseError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m pipeline.UserDatabases", description=__doc__.split("\n\n")[0])
    add_cli_arguments(parser)
    return run_cli(parser.parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())
