"""Build a chewBBACA-format cgMLST schema from user-provided sources and register it.

This reuses the conversion approach proven in
``scripts/build_chewbbaca_schemes_from_ncs.py`` (the cgMLST.org / Ridom NCS
importer): each locus becomes a normalized per-allele FASTA at the schema root,
with a ``short/`` representative dir and the ``.schema_config`` / ``.genes_list``
pickles chewBBACA expects — so no chewBBACA subcommand is required to adapt a
schema.

Supported sources (see :func:`build_schema`):
  * a folder of per-locus allele FASTAs (Ridom style; a stray combined
    "all-loci" FASTA is ignored via the ``targets.csv`` manifest or header
    detection),
  * a ``.zip`` of per-locus FASTAs,
  * a single combined multi-locus FASTA (locus encoded in each header),
  * a single FASTA holding all alleles of ONE locus,
  * an already-built chewBBACA schema directory (validated and copied).
"""
from __future__ import annotations

import csv
import json
import os
import pickle
import re
import shutil
import struct
import tempfile
import zipfile
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Optional

CHEWBBACA_SHORT_DIR = "short"
CHEWBBACA_GENE_LIST = ".genes_list"
CHEWBBACA_SCHEMA_CONFIG = ".schema_config"
CHEWBBACA_SCHEME_CATALOG = "scheme_species_catalog.json"
FASTA_EXTS = (".fasta", ".fa", ".fna", ".ffn", ".faa")

# chewBBACA 3.5.4 schema parameters (matches the existing NCS-built schemas).
SCHEMA_CONFIG = {
    "bsr": [0.6],
    "prodigal_training_file": [None],
    "translation_table": [11],
    "minimum_locus_length": [0],
    "chewBBACA_version": ["3.5.4"],
    "size_threshold": [0.2],
    "word_size": [5],
    "window_size": [5],
    "cluster_sim": [0.2],
    "representative_filter": [0.9],
    "intraCluster_filter": [0.9],
}

ZIP_LOCAL_FILE_SIG = b"PK\x03\x04"
ZIP_DATA_DESCRIPTOR_SIG = b"PK\x07\x08"
HTML_PREFIXES = (b"<html", b"<!doctype html")

Progress = Optional[Callable[[str], None]]


def _emit(progress: Progress, message: str) -> None:
    if progress is not None:
        try:
            progress(message)
        except Exception:
            pass


# ---------------------------------------------------------------------------
# FASTA primitives (adapted from scripts/build_chewbbaca_schemes_from_ncs.py)
# ---------------------------------------------------------------------------
def wrapped_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width))


def parse_fasta_records(data: bytes, *, source_name: str) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    parts: list[str] = []
    text = data.decode("utf-8", errors="replace")
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None and parts:
                records.append((header, "".join(parts)))
            header = line[1:].strip()
            parts = []
            continue
        if header is None:
            continue
        parts.append(line)
    if header is not None and parts:
        records.append((header, "".join(parts)))
    if not records:
        raise ValueError(f"{source_name} does not contain readable FASTA records")
    return records


def first_fasta_record(data: bytes, *, source_name: str) -> tuple[str, str]:
    header: str | None = None
    parts: list[str] = []
    text = data.decode("utf-8", errors="replace")
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                break
            header = line[1:].strip()
            continue
        if header is None:
            continue
        parts.append(line)
    if not header or not parts:
        raise ValueError(f"{source_name} does not contain a readable FASTA record")
    return header, "".join(parts)


def normalized_allele_header(locus_name: str, raw_header: str, seen: set[str]) -> str:
    locus_stem = Path(locus_name).stem
    token = raw_header.split()[0].strip()
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", token).strip("_") or "1"
    header = f"{locus_stem}_{token}"
    if header in seen:
        suffix = 2
        while f"{header}_{suffix}" in seen:
            suffix += 1
        header = f"{header}_{suffix}"
    seen.add(header)
    return header


def normalize_locus_fasta(locus_name: str, locus_data: bytes) -> bytes:
    records = parse_fasta_records(locus_data, source_name=locus_name)
    seen: set[str] = set()
    out: list[str] = []
    for raw_header, sequence in records:
        header = normalized_allele_header(locus_name, raw_header, seen)
        out.append(f">{header}")
        out.append(wrapped_sequence(sequence))
    return ("\n".join(out) + "\n").encode("utf-8")


def write_short_fasta(short_path: Path, locus_data: bytes, *, source_name: str) -> None:
    header, sequence = first_fasta_record(locus_data, source_name=source_name)
    short_path.write_text(f">{header}\n{wrapped_sequence(sequence)}\n", encoding="utf-8")


class SchemeStager:
    """Stages normalized loci into a temp dir, plus short/ representatives."""

    def __init__(self, temp_dir: Path, expected_names: Iterable[str] | None) -> None:
        self.temp_dir = temp_dir
        self.short_dir = temp_dir / CHEWBBACA_SHORT_DIR
        self.short_dir.mkdir(parents=True, exist_ok=True)
        self.expected_names = set(expected_names or [])
        self.written_names: list[str] = []
        self.written_set: set[str] = set()
        self.unexpected_names: list[str] = []
        self.duplicate_names: list[str] = []

    def write_locus(self, locus_name: str, locus_data: bytes) -> None:
        if self.expected_names and locus_name not in self.expected_names:
            self.unexpected_names.append(locus_name)
            return
        if locus_name in self.written_set:
            self.duplicate_names.append(locus_name)
            return
        normalized = normalize_locus_fasta(locus_name, locus_data)
        (self.temp_dir / locus_name).write_bytes(normalized)
        short_target = self.short_dir / f"{Path(locus_name).stem}_short.fasta"
        write_short_fasta(short_target, normalized, source_name=locus_name)
        self.written_set.add(locus_name)
        self.written_names.append(locus_name)

    def finalize_schema_metadata(self, expected_order: list[str]) -> None:
        ordered = [name for name in expected_order if name in self.written_set]
        if not ordered:
            ordered = list(self.written_names)
        with (self.temp_dir / CHEWBBACA_SCHEMA_CONFIG).open("wb") as fh:
            pickle.dump(SCHEMA_CONFIG, fh, protocol=pickle.HIGHEST_PROTOCOL)
        with (self.temp_dir / CHEWBBACA_GENE_LIST).open("wb") as fh:
            pickle.dump(ordered, fh, protocol=pickle.HIGHEST_PROTOCOL)


# ---------------------------------------------------------------------------
# Source readers
# ---------------------------------------------------------------------------
def read_expected_loci(targets_csv: Path) -> list[str]:
    expected: list[str] = []
    seen: set[str] = set()
    with targets_csv.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            locus_name = next(
                (field.strip() for field in reversed(row) if field.strip().lower().endswith(".fasta")),
                "",
            )
            if not locus_name or locus_name in seen:
                continue
            seen.add(locus_name)
            expected.append(locus_name)
    return expected


def convert_zip_archive(alleles_path: Path, stager: SchemeStager, progress: Progress = None) -> str:
    with zipfile.ZipFile(alleles_path) as archive:
        infos = [i for i in archive.infolist() if not i.is_dir()]
        for idx, info in enumerate(infos):
            stager.write_locus(Path(info.filename).name, archive.read(info))
            if progress and idx % 200 == 0:
                _emit(progress, f"Staging loci… {idx + 1}/{len(infos)}")
    return "zip_archive"


def _stream_zip_members(alleles_path: Path, on_member: Callable[[str, bytes], None]) -> str:
    data = alleles_path.read_bytes()
    pos = 0
    while pos + 30 <= len(data) and data[pos:pos + 4] == ZIP_LOCAL_FILE_SIG:
        header = data[pos + 4:pos + 30]
        (_v, _f, compression, _mt, _md, _crc, _cs, _us, name_len, extra_len) = struct.unpack(
            "<HHHHHIIIHH", header
        )
        name_start = pos + 30
        name = data[name_start:name_start + name_len].decode("utf-8", errors="replace")
        payload_start = name_start + name_len + extra_len
        if compression != 8:
            raise ValueError(f"Unsupported ZIP stream compression method {compression} for {name}")
        decompressor = zlib.decompressobj(-15)
        chunks: list[bytes] = []
        read_pos = payload_start
        while True:
            chunk = data[read_pos:read_pos + 1024 * 1024]
            if not chunk:
                raise EOFError(f"truncated while decoding {name}")
            out = decompressor.decompress(chunk)
            if out:
                chunks.append(out)
            if decompressor.unused_data:
                read_pos += len(chunk) - len(decompressor.unused_data)
                break
            read_pos += len(chunk)
        on_member(Path(name).name, b"".join(chunks))
        pos = read_pos + (16 if data[read_pos:read_pos + 4] == ZIP_DATA_DESCRIPTOR_SIG else 12)
    return ""


def convert_zip_stream(alleles_path: Path, stager: SchemeStager) -> str:
    _stream_zip_members(alleles_path, stager.write_locus)
    return "zip_stream"


def locus_name_from_header(header: str, expected_names: set[str]) -> str:
    token = header.split()[0]
    candidates: list[str] = []
    if token.lower().endswith(".fasta"):
        candidates.append(token)
    candidates.append(f"{token}.fasta")
    if "_" in token:
        prefix, suffix = token.rsplit("_", 1)
        if suffix.isdigit():
            candidates.append(f"{prefix}.fasta")
    for candidate in candidates:
        if candidate in expected_names:
            return candidate
    if "_" in token:
        prefix, suffix = token.rsplit("_", 1)
        if suffix.isdigit():
            return f"{prefix}.fasta"
    return f"{token}.fasta"


def convert_combined_fasta(
    alleles_path: Path, stager: SchemeStager, expected_names: list[str]
) -> str:
    groups: dict[str, list[str]] = {}
    current_header: str | None = None
    current_seq: list[str] = []

    def flush() -> None:
        nonlocal current_header, current_seq
        if current_header is None:
            return
        locus_name = locus_name_from_header(current_header, set(expected_names))
        groups.setdefault(locus_name, []).append(
            f">{current_header}\n{wrapped_sequence(''.join(current_seq))}\n"
        )
        current_header = None
        current_seq = []

    with alleles_path.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                current_header = line[1:].strip()
                continue
            if current_header is not None:
                current_seq.append(line)
    flush()
    for locus_name, records in groups.items():
        stager.write_locus(locus_name, "".join(records).encode("utf-8"))
    return "combined_fasta"


def _looks_like_html(path: Path) -> bool:
    try:
        with path.open("rb") as fh:
            prefix = fh.read(256).lstrip().lower()
        return any(prefix.startswith(m) for m in HTML_PREFIXES)
    except Exception:
        return False


def detect_file_format(path: Path) -> str:
    """Classify a FILE source: zip_archive / zip_stream / combined_fasta / html_error / unknown."""
    if _looks_like_html(path):
        return "html_error"
    if zipfile.is_zipfile(path):
        return "zip_archive"
    try:
        with path.open("rb") as fh:
            prefix = fh.read(8)
    except Exception:
        return "unknown"
    if prefix.startswith(ZIP_LOCAL_FILE_SIG):
        return "zip_stream"
    if prefix.lstrip().startswith(b">"):
        return "combined_fasta"
    return "unknown"


# ---------------------------------------------------------------------------
# Folder handling + the Ridom "per-locus + combined file" quirk
# ---------------------------------------------------------------------------
def _list_fastas(folder: Path) -> list[Path]:
    return sorted(
        (p for p in folder.iterdir() if p.is_file() and p.suffix.lower() in FASTA_EXTS),
        key=lambda p: p.name.lower(),
    )


def _find_targets_csv(folder: Path) -> Optional[Path]:
    for name in ("targets.csv", "targets.tsv"):
        cand = folder / name
        if cand.is_file():
            return cand
    return None


def is_combined_fasta(path: Path, *, sample: int = 60) -> bool:
    """True when a FASTA holds >=2 distinct, *named* loci (a Ridom combined file).

    Per-locus allele files have bare numeric headers (>1, >2, …) whose derived
    locus stems are numeric -> not combined. A combined file has locus-named
    headers (>ACICU_RS00005_1, >ACICU_RS00010_1, …) -> combined.
    """
    stems: set[str] = set()
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                token = line[1:].split()[0].strip() if len(line) > 1 else ""
                if not token:
                    continue
                stem = token
                if "_" in token:
                    prefix, suffix = token.rsplit("_", 1)
                    if suffix.isdigit():
                        stem = prefix
                stems.add(stem)
                if len(stems) > sample:
                    break
    except Exception:
        return False
    lettered = {s for s in stems if re.search(r"[A-Za-z]", s)}
    return len(lettered) >= 2


def is_prebuilt_schema(path: Path) -> bool:
    return (
        path.is_dir()
        and (path / CHEWBBACA_GENE_LIST).is_file()
        and (path / CHEWBBACA_SHORT_DIR).is_dir()
    )


def slugify_species(name: str, *, fallback: str = "schema") -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", str(name).strip().lower()).strip("_")
    return slug or fallback


# ---------------------------------------------------------------------------
# High-level build
# ---------------------------------------------------------------------------
@dataclass
class SchemaBuildResult:
    target_dir: Path
    schema_id: str
    written_loci: int
    expected_loci: int
    missing_loci: int
    source_format: str
    warning: str


def _replace_output_dir(target_dir: Path, temp_dir: Path, *, overwrite: bool) -> None:
    if target_dir.exists():
        if not overwrite:
            raise FileExistsError(f"schema already exists: {target_dir}")
        shutil.rmtree(target_dir)
    target_dir.parent.mkdir(parents=True, exist_ok=True)
    temp_dir.replace(target_dir)


def _new_temp(target_dir: Path) -> Path:
    target_dir.parent.mkdir(parents=True, exist_ok=True)
    return Path(tempfile.mkdtemp(prefix=f".{target_dir.name}_", dir=str(target_dir.parent)))


def _finalize(
    stager: SchemeStager,
    expected_order: list[str],
    target_dir: Path,
    temp_root: Path,
    *,
    overwrite: bool,
    source_format: str,
    warning: str,
) -> SchemaBuildResult:
    if not stager.written_names:
        raise ValueError("no loci were produced from the source")
    stager.finalize_schema_metadata(expected_order or list(stager.written_names))
    _replace_output_dir(target_dir, temp_root, overwrite=overwrite)
    missing = [n for n in expected_order if n not in stager.written_set]
    dup = sorted(set(stager.duplicate_names))
    parts = [warning] if warning else []
    if dup:
        parts.append(f"ignored {len(dup)} duplicate locus file(s)")
    return SchemaBuildResult(
        target_dir=target_dir,
        schema_id=target_dir.name,
        written_loci=len(stager.written_names),
        expected_loci=len(expected_order) if expected_order else len(stager.written_names),
        missing_loci=len(missing),
        source_format=source_format,
        warning="; ".join(parts),
    )


def build_schema(
    source: Path,
    target_dir: Path,
    *,
    single_fasta_mode: str = "loci_set",  # "loci_set" | "single_locus"
    overwrite: bool = False,
    progress: Progress = None,
) -> SchemaBuildResult:
    """Build a chewBBACA schema at ``target_dir`` from ``source``.

    ``source`` may be a directory (prebuilt schema, or a folder of per-locus /
    combined FASTAs) or a file (.zip, or a single FASTA). ``single_fasta_mode``
    disambiguates a lone FASTA: a multi-locus set vs. all alleles of one locus.
    """
    source = Path(source).expanduser()
    target_dir = Path(target_dir).expanduser()
    if not source.exists():
        raise FileNotFoundError(f"source not found: {source}")

    # 1) already-built chewBBACA schema directory -> validate + copy
    if source.is_dir() and is_prebuilt_schema(source):
        _emit(progress, "Copying prebuilt chewBBACA schema…")
        if target_dir.exists():
            if not overwrite:
                raise FileExistsError(f"schema already exists: {target_dir}")
            shutil.rmtree(target_dir)
        target_dir.parent.mkdir(parents=True, exist_ok=True)
        tmp = _new_temp(target_dir)
        shutil.rmtree(tmp)  # copytree needs a non-existent dest
        try:
            shutil.copytree(source, tmp)
            tmp.replace(target_dir)
        except Exception:
            shutil.rmtree(tmp, ignore_errors=True)
            raise
        try:
            genes = pickle.loads((target_dir / CHEWBBACA_GENE_LIST).read_bytes())
            n = len(genes) if isinstance(genes, list) else 0
        except Exception:
            n = sum(1 for p in target_dir.glob("*.fasta"))
        return SchemaBuildResult(target_dir, target_dir.name, n, n, 0, "prebuilt_schema", "")

    # 2) directory of loci FASTAs (Ridom)
    if source.is_dir():
        return _build_from_folder(source, target_dir, overwrite=overwrite, progress=progress)

    # 3) single file
    fmt = detect_file_format(source)
    if fmt == "html_error":
        raise ValueError("source is an HTML error page, not a FASTA/zip")
    temp_root = _new_temp(target_dir)
    try:
        if fmt in ("zip_archive", "zip_stream"):
            _emit(progress, "Reading zip archive…")
            stager = SchemeStager(temp_root, None)
            sfmt = convert_zip_archive(source, stager, progress) if fmt == "zip_archive" \
                else convert_zip_stream(source, stager)
            return _finalize(stager, list(stager.written_names), target_dir, temp_root,
                             overwrite=overwrite, source_format=sfmt, warning="")
        # a lone FASTA
        if single_fasta_mode == "single_locus":
            _emit(progress, "Building single-locus schema…")
            stager = SchemeStager(temp_root, None)
            locus_name = f"{slugify_species(source.stem, fallback='locus')}.fasta"
            stager.write_locus(locus_name, source.read_bytes())
            return _finalize(stager, [locus_name], target_dir, temp_root,
                             overwrite=overwrite, source_format="single_locus_fasta", warning="")
        _emit(progress, "Splitting combined FASTA into loci…")
        stager = SchemeStager(temp_root, None)
        convert_combined_fasta(source, stager, [])
        return _finalize(stager, list(stager.written_names), target_dir, temp_root,
                         overwrite=overwrite, source_format="combined_fasta", warning="")
    except Exception:
        shutil.rmtree(temp_root, ignore_errors=True)
        raise


def _build_from_folder(
    folder: Path, target_dir: Path, *, overwrite: bool, progress: Progress
) -> SchemaBuildResult:
    fastas = _list_fastas(folder)
    if not fastas:
        raise ValueError(f"no FASTA files found in {folder}")

    targets_csv = _find_targets_csv(folder)
    temp_root = _new_temp(target_dir)
    try:
        if targets_csv is not None:
            # Manifest-driven: only files named in targets.csv become loci; anything
            # else (e.g. a stray combined all-loci FASTA) is skipped by the stager.
            expected = read_expected_loci(targets_csv)
            _emit(progress, f"Staging {len(expected)} loci (targets.csv)…")
            stager = SchemeStager(temp_root, expected)
            for idx, f in enumerate(fastas):
                stager.write_locus(f.name, f.read_bytes())
                if idx % 200 == 0:
                    _emit(progress, f"Staging loci… {idx + 1}/{len(fastas)}")
            return _finalize(stager, expected, target_dir, temp_root,
                             overwrite=overwrite, source_format="folder+targets", warning="")

        # No manifest: split per-locus files from any combined all-loci file(s).
        combined = [f for f in fastas if is_combined_fasta(f)]
        per_locus = [f for f in fastas if f not in combined]
        stager = SchemeStager(temp_root, None)
        if per_locus:
            _emit(progress, f"Staging {len(per_locus)} loci…")
            for idx, f in enumerate(per_locus):
                stager.write_locus(f.name, f.read_bytes())
                if idx % 200 == 0:
                    _emit(progress, f"Staging loci… {idx + 1}/{len(per_locus)}")
            warning = f"ignored {len(combined)} combined-style FASTA file(s)" if combined else ""
            return _finalize(stager, [f.name for f in per_locus], target_dir, temp_root,
                             overwrite=overwrite, source_format="folder", warning=warning)
        if combined:
            _emit(progress, "Splitting combined FASTA into loci…")
            for f in combined:
                convert_combined_fasta(f, stager, [])
            return _finalize(stager, list(stager.written_names), target_dir, temp_root,
                             overwrite=overwrite, source_format="folder_combined", warning="")
        raise ValueError(f"no usable loci FASTAs found in {folder}")
    except Exception:
        shutil.rmtree(temp_root, ignore_errors=True)
        raise


# ---------------------------------------------------------------------------
# Registration
# ---------------------------------------------------------------------------
def _atomic_write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def register_schema(
    cgmlst_dir: Path,
    schema_id: str,
    species: list[str],
    mlst_schemes: list[str],
    *,
    catalog_path: Optional[Path] = None,
    update_config: bool = True,
) -> Path:
    """Add/replace a catalog entry and (optionally) the config available_schemas map.

    Raises ValueError on an MLST-scheme alias already claimed by another schema.
    """
    cgmlst_dir = Path(cgmlst_dir)
    catalog_path = Path(catalog_path) if catalog_path else cgmlst_dir / CHEWBBACA_SCHEME_CATALOG

    data: dict = {"version": 1, "schemes": []}
    if catalog_path.is_file():
        try:
            loaded = json.loads(catalog_path.read_text(encoding="utf-8"))
            if isinstance(loaded, dict):
                data = loaded
        except Exception:
            pass
    schemes = data.get("schemes")
    if not isinstance(schemes, list):
        schemes = []

    mlst_clean = [a.strip() for a in mlst_schemes if a and a.strip()]
    alias_norm = {a.lower() for a in mlst_clean}
    for entry in schemes:
        if not isinstance(entry, dict) or str(entry.get("schema_id")) == schema_id:
            continue
        for alias in entry.get("mlst_schemes") or []:
            if str(alias).strip().lower() in alias_norm:
                raise ValueError(
                    f"MLST scheme '{alias}' already maps to schema '{entry.get('schema_id')}'"
                )

    new_entry = {
        "schema_id": schema_id,
        "schema_dir": schema_id,
        "species": [s.strip() for s in species if s and s.strip()],
        "mlst_schemes": mlst_clean,
    }
    schemes = [e for e in schemes if not (isinstance(e, dict) and str(e.get("schema_id")) == schema_id)]
    schemes.append(new_entry)
    schemes.sort(key=lambda e: str(e.get("schema_id")))
    data["schemes"] = schemes
    data.setdefault("version", 1)
    _atomic_write_json(catalog_path, data)

    if update_config:
        try:
            from pipeline.ToolEnv import save_app_config

            save_app_config({
                "chewbbaca": {
                    "available_schemas": {schema_id: str((cgmlst_dir / schema_id).resolve())}
                }
            })
        except Exception:
            pass
    return catalog_path
