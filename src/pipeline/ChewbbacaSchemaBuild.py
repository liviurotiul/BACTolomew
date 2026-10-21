"""Build a chewBBACA-format cgMLST schema from user-provided sources and register it.

Source loci are staged as one allele FASTA per locus (headers normalized to
``<locus>_<allele>``) and then adapted with chewBBACA ``PrepExternalSchema``
(run in the tools image like every other tool). PrepExternalSchema translates
every allele, writes it in coding orientation (cgMLST.org / Ridom alleles are
stored on the reference strand, so about half of them are reverse
complemented), drops alleles that are not a complete CDS, selects the
representatives in ``short/`` and writes ``.schema_config`` / ``.genes_list``.
The rejected alleles and loci are kept in ``schema_build_report/`` inside the
schema, next to a ``build_info.json`` with the before/after counts.

Supported sources (see :func:`build_schema`):
  * a folder of per-locus allele FASTAs (Ridom style; a stray combined
    "all-loci" FASTA is ignored via the ``targets.csv`` manifest or header
    detection),
  * a ``.zip`` of per-locus FASTAs (optionally filtered by a ``targets.csv``),
  * a single combined multi-locus FASTA (locus encoded in each header),
  * a single FASTA holding all alleles of ONE locus,
  * an already-built chewBBACA schema directory (copied when every
    representative is a forward CDS, otherwise re-adapted from its locus
    FASTAs).
"""
from __future__ import annotations

import csv
import fcntl
import json
import os
import re
import shutil
import struct
import subprocess
import tempfile
import zipfile
import zlib
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable, Iterator, Optional

CHEWBBACA_SHORT_DIR = "short"
CHEWBBACA_GENE_LIST = ".genes_list"
CHEWBBACA_SCHEMA_CONFIG = ".schema_config"
CHEWBBACA_SCHEME_CATALOG = "scheme_species_catalog.json"
FASTA_EXTS = (".fasta", ".fa", ".fna", ".ffn", ".faa")

# PrepExternalSchema settings (chewBBACA defaults; stored in .schema_config).
PREP_BSR = "0.6"
PREP_TRANSLATION_TABLE = 11
# Report written inside every built schema (not a .fasta, so chewBBACA ignores it).
SCHEMA_BUILD_REPORT_DIR = "schema_build_report"
BUILD_INFO_FILE = "build_info.json"
INVALID_ALLELES_FILE = "invalid_alleles.txt"
INVALID_LOCI_FILE = "invalid_loci.txt"
SUMMARY_STATS_FILE = "summary_stats.tsv"
# Runtime caches chewBBACA AlleleCall writes into a schema; never copied on import.
CHEWBBACA_RUNTIME_CACHES = ("pre_computed", "loci_modes", "short/self_scores")
# Held (shared) in a schema by a run's chewBBACA step; a replace needs it exclusively.
SCHEMA_USE_LOCK = ".bactolomew_in_use.lock"

ZIP_LOCAL_FILE_SIG = b"PK\x03\x04"
ZIP_DATA_DESCRIPTOR_SIG = b"PK\x07\x08"
ZIP_CENTRAL_DIR_SIG = b"PK\x01\x02"
ZIP_END_OF_DIR_SIG = b"PK\x05\x06"
HTML_PREFIXES = (b"<html", b"<!doctype html")

Progress = Optional[Callable[[str], None]]


ShouldCancel = Optional[Callable[[], bool]]


class Cancelled(RuntimeError):
    """A build (or install) was cancelled by the user."""


def _emit(progress: Progress, message: str) -> None:
    if progress is not None:
        try:
            progress(message)
        except Exception:
            pass


def run_watched(
    cmd: list[str],
    *,
    watch: Optional[Callable[[], None]] = None,
    should_cancel: ShouldCancel = None,
    interval: float = 1.0,
) -> subprocess.CompletedProcess:
    """subprocess.run(cmd, capture_output=True, text=True) that calls ``watch`` every
    ``interval`` s while the command runs, and kills it and raises Cancelled once
    ``should_cancel`` says so."""
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    while True:
        try:
            stdout, stderr = proc.communicate(timeout=interval)
            return subprocess.CompletedProcess(cmd, proc.returncode, stdout, stderr)
        except subprocess.TimeoutExpired:
            pass
        if should_cancel is not None and should_cancel():
            proc.kill()
            proc.communicate()
            raise Cancelled("cancelled")
        if watch is not None:
            try:
                watch()
            except Exception:  # noqa: BLE001 - only progress
                pass


# ---------------------------------------------------------------------------
# FASTA primitives
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


def normalize_locus_fasta(locus_name: str, locus_data: bytes) -> tuple[bytes, int]:
    """Rename headers to ``<locus>_<allele>`` and rewrap; returns (FASTA bytes, allele count)."""
    records = parse_fasta_records(locus_data, source_name=locus_name)
    seen: set[str] = set()
    out: list[str] = []
    for raw_header, sequence in records:
        header = normalized_allele_header(locus_name, raw_header, seen)
        out.append(f">{header}")
        out.append(wrapped_sequence(sequence))
    return ("\n".join(out) + "\n").encode("utf-8"), len(records)


def read_fasta_file(path: Path) -> list[tuple[str, str]]:
    """All (header, sequence) records of a FASTA file; empty list when it has none."""
    try:
        return parse_fasta_records(Path(path).read_bytes(), source_name=str(path))
    except ValueError:
        return []


# ---------------------------------------------------------------------------
# CDS orientation (translation table 11)
# ---------------------------------------------------------------------------
# Stop and start codons per NCBI genetic code (11: bacteria; 4: Mycoplasma, TGA = Trp).
_STOP_CODONS = {
    11: frozenset({"TAA", "TAG", "TGA"}),
    4: frozenset({"TAA", "TAG"}),
}
_START_CODONS = {
    11: frozenset({"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"}),
    4: frozenset({"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA", "TTA"}),
}
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(_COMPLEMENT)[::-1]


def is_forward_cds(sequence: str, table: int = 11) -> bool:
    """True when ``sequence`` is a complete CDS as written: start..stop, no internal stop."""
    stops = _STOP_CODONS.get(int(table), _STOP_CODONS[11])
    starts = _START_CODONS.get(int(table), _START_CODONS[11])
    seq = sequence.upper()
    if len(seq) < 6 or len(seq) % 3 or set(seq) - set("ACGT"):
        return False
    if seq[:3] not in starts or seq[-3:] not in stops:
        return False
    return all(seq[i:i + 3] not in stops for i in range(3, len(seq) - 3, 3))


def cds_orientation(sequence: str, table: int = 11) -> str:
    """'forward' (valid as written), 'reverse' (valid only reverse complemented) or 'invalid'."""
    if is_forward_cds(sequence, table):
        return "forward"
    if is_forward_cds(reverse_complement(sequence), table):
        return "reverse"
    return "invalid"


def read_build_info(schema_dir: Path) -> dict:
    """The schema's schema_build_report/build_info.json ({} when absent or unreadable)."""
    try:
        data = json.loads((Path(schema_dir) / SCHEMA_BUILD_REPORT_DIR / BUILD_INFO_FILE).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def schema_translation_table(schema_dir: Path, default: int = 11) -> int:
    """The genetic code recorded in the schema's build_info.json (``default`` when unknown).

    The schema's own .schema_config is a pickle; it is never loaded on the host
    (an imported schema could carry any pickle). :func:`read_schema_metadata`
    reads it in the tools image when a schema has no build record yet.
    """
    try:
        value = read_build_info(schema_dir).get("translation_table")
        return int(value) if value else default
    except (TypeError, ValueError):
        return default


# Runs in the chewbbaca env of the tools image: the schema's pickles are loaded there,
# never on the host. Prints {"genes_list": <loci listed or -1>, "translation_table": <n or null>}.
_SCHEMA_METADATA_SCRIPT = r"""
import json, os, pickle, sys
d = sys.argv[1]
out = {"genes_list": -1, "translation_table": None}
try:
    with open(os.path.join(d, ".genes_list"), "rb") as fh:
        genes = pickle.load(fh)
    out["genes_list"] = len(genes) if isinstance(genes, list) else -1
except Exception:
    pass
try:
    with open(os.path.join(d, ".schema_config"), "rb") as fh:
        config = pickle.load(fh)
    value = config.get("translation_table")
    value = value[0] if isinstance(value, list) else value
    out["translation_table"] = int(value) if value else None
except Exception:
    pass
print(json.dumps(out))
"""


@dataclass
class SchemaMetadata:
    genes_list: int = -1  # loci listed in .genes_list (-1: missing or unreadable)
    translation_table: Optional[int] = None


def read_schema_metadata(schema_dir: Path, *, python_cmd: list[str] | None = None) -> SchemaMetadata:
    """Read a schema's .genes_list and .schema_config pickles in the tools image."""
    from pipeline.ToolEnv import gui_tool_cmd

    schema_dir = Path(schema_dir).resolve()
    cmd = list(python_cmd) if python_cmd else gui_tool_cmd("chewbbaca", "python", mount_paths=[schema_dir])
    res = subprocess.run(cmd + ["-c", _SCHEMA_METADATA_SCRIPT, str(schema_dir)],
                         capture_output=True, text=True, check=False)
    lines = [ln for ln in (res.stdout or "").splitlines() if ln.strip().startswith("{")]
    if res.returncode != 0 or not lines:
        tail = "\n".join((res.stderr or res.stdout or "").strip().splitlines()[-5:])
        raise RuntimeError(f"could not read the schema settings of {schema_dir} in the tools image: {tail}")
    data = json.loads(lines[-1])
    table = data.get("translation_table")
    return SchemaMetadata(int(data.get("genes_list", -1)), int(table) if table else None)


# ---------------------------------------------------------------------------
# Schema statistics and validation
# ---------------------------------------------------------------------------
@dataclass
class SchemaStats:
    loci: int = 0
    genes_list: int = -1
    alleles: int = 0
    representatives: int = 0
    self_scores: int = -1
    reps_missing_self_score: int = -1
    alleles_forward: int = 0
    alleles_reverse_only: int = 0
    alleles_invalid: int = 0
    reps_forward: int = 0
    reps_reverse_only: int = 0
    reps_invalid: int = 0
    loci_without_short: int = 0


def _locus_fastas(schema_dir: Path) -> list[Path]:
    return sorted(
        (p for p in Path(schema_dir).iterdir() if p.is_file() and p.name.lower().endswith(".fasta")),
        key=lambda p: p.name,
    )


def schema_stats(
    schema_dir: Path,
    *,
    orientation: bool = True,
    metadata: SchemaMetadata | None = None,
) -> SchemaStats:
    """Count loci, alleles, representatives, self-scores and CDS orientation of a schema.

    ``genes_list`` stays -1 unless ``metadata`` (read in the tools image) is given;
    the genetic code comes from ``metadata`` or else the schema's build_info.json.
    """
    from pipeline.ChewBBACAWorkflow import load_self_scores, schema_representative_ids

    schema_dir = Path(schema_dir)
    stats = SchemaStats()
    table = (metadata.translation_table if metadata else None) or schema_translation_table(schema_dir)
    loci = _locus_fastas(schema_dir)
    stats.loci = len(loci)
    if metadata is not None:
        stats.genes_list = metadata.genes_list
    short_dir = schema_dir / CHEWBBACA_SHORT_DIR
    for locus in loci:
        records = read_fasta_file(locus)
        stats.alleles += len(records)
        if orientation:
            for _h, seq in records:
                kind = cds_orientation(seq, table)
                if kind == "forward":
                    stats.alleles_forward += 1
                elif kind == "reverse":
                    stats.alleles_reverse_only += 1
                else:
                    stats.alleles_invalid += 1
        short = short_dir / f"{locus.stem}_short.fasta"
        reps = read_fasta_file(short) if short.is_file() else []
        if not reps:
            stats.loci_without_short += 1
        stats.representatives += len(reps)
        if orientation:
            for _h, seq in reps:
                kind = cds_orientation(seq, table)
                if kind == "forward":
                    stats.reps_forward += 1
                elif kind == "reverse":
                    stats.reps_reverse_only += 1
                else:
                    stats.reps_invalid += 1
    scores = load_self_scores(schema_dir)
    if scores is not None:
        stats.self_scores = len(scores)
        rep_ids = schema_representative_ids(schema_dir)
        stats.reps_missing_self_score = len(rep_ids - set(scores))
    return stats


def adapted_schema_issues(schema_dir: Path, *, metadata: SchemaMetadata | None = None) -> list[str]:
    """Problems that make an adapted schema unusable (empty list when it is valid).

    The .genes_list pickle is checked in the tools image (``metadata`` when the
    caller already read it).
    """
    schema_dir = Path(schema_dir)
    issues: list[str] = []
    for name in (CHEWBBACA_GENE_LIST, CHEWBBACA_SCHEMA_CONFIG):
        if not (schema_dir / name).is_file():
            issues.append(f"missing {name}")
    if issues:
        return issues
    stats = schema_stats(schema_dir, metadata=metadata or read_schema_metadata(schema_dir))
    if stats.loci == 0:
        issues.append("no locus FASTA files")
    if stats.genes_list != stats.loci:
        issues.append(f".genes_list lists {stats.genes_list} loci but {stats.loci} locus FASTAs exist")
    if stats.loci_without_short:
        issues.append(f"{stats.loci_without_short} loci have no representative in short/")
    if stats.reps_reverse_only:
        issues.append(f"{stats.reps_reverse_only} representatives are valid CDS only when reverse complemented")
    if stats.reps_invalid:
        issues.append(f"{stats.reps_invalid} representatives are not a valid CDS")
    if stats.alleles_reverse_only:
        issues.append(f"{stats.alleles_reverse_only} alleles are valid CDS only when reverse complemented")
    return issues


class SchemeStager:
    """Stages normalized per-locus allele FASTAs (the PrepExternalSchema input)."""

    def __init__(self, temp_dir: Path, expected_names: Iterable[str] | None) -> None:
        self.temp_dir = temp_dir
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        self.expected_names = set(expected_names or [])
        self.written_names: list[str] = []
        self.written_set: set[str] = set()
        self.unexpected_names: list[str] = []
        self.duplicate_names: list[str] = []
        self.allele_count = 0

    def write_locus(self, locus_name: str, locus_data: bytes) -> None:
        if self.expected_names and locus_name not in self.expected_names:
            self.unexpected_names.append(locus_name)
            return
        if locus_name in self.written_set:
            self.duplicate_names.append(locus_name)
            return
        normalized, count = normalize_locus_fasta(locus_name, locus_data)
        stem = Path(locus_name).stem
        (self.temp_dir / f"{stem}.fasta").write_bytes(normalized)
        self.allele_count += count
        self.written_set.add(locus_name)
        self.written_names.append(locus_name)


# ---------------------------------------------------------------------------
# chewBBACA PrepExternalSchema
# ---------------------------------------------------------------------------
def _chewbbaca_base_cmd(mount_paths: Iterable[Path]) -> list[str]:
    from pipeline.ToolEnv import gui_tool_cmd

    return gui_tool_cmd("chewbbaca", "chewBBACA.py", mount_paths=[str(p) for p in mount_paths])


def _default_cpu() -> int:
    from pipeline.ToolEnv import effective_core_budget

    return max(1, min(32, effective_core_budget() - 1))


def run_prep_external_schema(
    external_dir: Path,
    out_dir: Path,
    *,
    cpu: int | None = None,
    base_cmd: list[str] | None = None,
    translation_table: int = PREP_TRANSLATION_TABLE,
    progress: Progress = None,
    should_cancel: ShouldCancel = None,
) -> Path:
    """Adapt ``external_dir`` (one FASTA per locus) into a chewBBACA schema at ``out_dir``.

    ``out_dir`` must not exist. PrepExternalSchema writes its invalid-allele,
    invalid-locus and summary reports beside ``out_dir``; they are moved into
    ``out_dir/schema_build_report``. Returns that report directory. While it runs,
    ``progress`` gets "Adapting loci… i/n" from the loci written so far; ``should_cancel``
    stops it (Cancelled).
    """
    external_dir = Path(external_dir).resolve()
    out_dir = Path(out_dir).resolve()
    if out_dir.exists():
        raise FileExistsError(f"PrepExternalSchema output already exists: {out_dir}")
    cmd = list(base_cmd) if base_cmd else _chewbbaca_base_cmd([external_dir, out_dir.parent])
    cmd += [
        "PrepExternalSchema",
        "-g", str(external_dir),
        "-o", str(out_dir),
        "--cpu", str(max(1, int(cpu or _default_cpu()))),
        "--bsr", PREP_BSR,
        "--t", str(int(translation_table)),
    ]
    _emit(progress, "Adapting loci with chewBBACA PrepExternalSchema…")
    staged = sum(1 for f in external_dir.iterdir() if f.name.lower().endswith(".fasta"))
    adapted = [0]

    def watch() -> None:
        count = sum(1 for f in out_dir.iterdir() if f.name.lower().endswith(".fasta")) if out_dir.is_dir() else 0
        if staged and count > adapted[0]:
            adapted[0] = count
            _emit(progress, f"Adapting loci… {min(count, staged)}/{staged}")

    res = run_watched(cmd, watch=watch, should_cancel=should_cancel, interval=2.0)
    log_text = f"$ {' '.join(cmd)}\n\n--- stdout ---\n{res.stdout or ''}\n--- stderr ---\n{res.stderr or ''}"
    if res.returncode != 0 or not (out_dir / CHEWBBACA_GENE_LIST).is_file():
        tail = "\n".join((res.stderr or res.stdout or "").strip().splitlines()[-15:])
        raise RuntimeError(
            f"chewBBACA PrepExternalSchema failed (exit code {res.returncode}).\n{tail}"
        )
    report_dir = out_dir / SCHEMA_BUILD_REPORT_DIR
    report_dir.mkdir(parents=True, exist_ok=True)
    (report_dir / "prep_external_schema.log").write_text(log_text, encoding="utf-8")
    prefix = out_dir.name
    for suffix in (INVALID_ALLELES_FILE, INVALID_LOCI_FILE, SUMMARY_STATS_FILE):
        produced = out_dir.parent / f"{prefix}_{suffix}"
        if produced.is_file():
            shutil.move(str(produced), str(report_dir / suffix))
        else:
            (report_dir / suffix).write_text("", encoding="utf-8")
    return report_dir


def _chewbbaca_version(log_text: str) -> str:
    match = re.search(r"chewBBACA version:\s*(\S+)", log_text)
    return match.group(1) if match else ""


def _nonempty_lines(path: Path) -> list[str]:
    try:
        return [ln for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    except OSError:
        return []


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
    """Read a zip's members in order from their local headers (for an archive whose
    central directory is damaged), with or without data descriptors. Raises ValueError
    when the archive is truncated or damaged: a member cut off or failing its CRC, or the
    data ending (or going on with something else) before the central directory."""
    data = alleles_path.read_bytes()
    view = memoryview(data)
    pos = 0
    count = 0

    def broken(why: str) -> ValueError:
        return ValueError(f"the zip archive is truncated or damaged: {why} ({count} members read)")

    while data[pos:pos + 4] == ZIP_LOCAL_FILE_SIG:
        if pos + 30 > len(data):
            raise broken("it ends inside a member header")
        (_v, flags, compression, _mt, _md, crc, comp_size, size, name_len, extra_len) = struct.unpack(
            "<HHHHHIIIHH", data[pos + 4:pos + 30]
        )
        name_start = pos + 30
        payload_start = name_start + name_len + extra_len
        if payload_start > len(data):
            raise broken("it ends inside a member header")
        name = data[name_start:name_start + name_len].decode("utf-8", errors="replace")
        has_descriptor = bool(flags & 0x08)  # sizes and CRC follow the data
        if compression == 8:
            decompressor = zlib.decompressobj(-15)
            chunks: list[bytes] = []
            read_pos = payload_start
            while not decompressor.eof:
                chunk = view[read_pos:read_pos + 1024 * 1024]
                if not chunk:
                    raise broken(f"it ends inside {name}")
                try:
                    chunks.append(decompressor.decompress(chunk))
                except zlib.error as exc:
                    raise broken(f"{name} is damaged ({exc})") from exc
                read_pos += len(chunk)
            end = read_pos - len(decompressor.unused_data)
            payload = b"".join(chunks)
        elif compression == 0 and not has_descriptor:
            end = payload_start + comp_size
            if end > len(data):
                raise broken(f"it ends inside {name}")
            payload = data[payload_start:end]
        else:
            raise ValueError(f"Unsupported ZIP stream compression method {compression} for {name}")
        if has_descriptor:
            if data[end:end + 4] == ZIP_DATA_DESCRIPTOR_SIG:
                end += 4
            if end + 12 > len(data):
                raise broken(f"it ends inside {name}")
            crc, _comp_size, size = struct.unpack("<III", data[end:end + 12])
            end += 12
        if zlib.crc32(payload) != crc or len(payload) & 0xFFFFFFFF != size:
            raise broken(f"{name} fails its checksum")
        on_member(Path(name).name, payload)
        count += 1
        pos = end
    if data[pos:pos + 4] not in (ZIP_CENTRAL_DIR_SIG, ZIP_END_OF_DIR_SIG):
        raise broken("it ends before the zip directory" if pos >= len(data) else f"unexpected data at byte {pos}")
    if not count:
        raise ValueError("the zip archive holds no members")
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
    source_loci: int = 0
    source_alleles: int = 0
    alleles: int = 0
    representatives: int = 0
    rejected_alleles: int = 0
    rejected_loci: list[str] = field(default_factory=list)


@dataclass
class BuildOptions:
    cpu: int | None = None
    base_cmd: list[str] | None = None  # chewBBACA command prefix (default: tools image)
    python_cmd: list[str] | None = None  # python of the chewbbaca env (default: tools image)
    translation_table: int = PREP_TRANSLATION_TABLE  # 4 for Mycoplasma
    should_cancel: ShouldCancel = None  # stops PrepExternalSchema and the test allele call


REPLACED_DIR = ".replaced"


class SchemaInUse(RuntimeError):
    """A run's chewBBACA step is using the schema that would be replaced."""


def _open_lock(path: Path) -> Optional[int]:
    """A descriptor of the lock file ``path`` (None when it cannot be created)."""
    try:
        return os.open(str(path), os.O_RDWR | os.O_CREAT, 0o664)
    except OSError:
        return None


@contextmanager
def schema_in_use(schema_dir: Path) -> Iterator[None]:
    """Mark ``schema_dir`` as used by a chewBBACA run (a shared lock) while the block runs.

    A replace (:func:`install_schema_dir`) refuses a schema marked so: AlleleCall in
    inferred mode writes novel alleles into it. A schema folder that is not writable
    is used without the mark.
    """
    fd = _open_lock(Path(schema_dir) / SCHEMA_USE_LOCK)
    try:
        if fd is not None:
            fcntl.flock(fd, fcntl.LOCK_SH)
        yield
    finally:
        if fd is not None:
            os.close(fd)


def install_schema_dir(built_dir: Path, target_dir: Path, *, replace: bool = False) -> Optional[Path]:
    """Move ``built_dir`` to ``target_dir``; returns where a replaced schema was kept.

    A schema is never deleted: chewBBACA adds the novel alleles of every inferred-mode
    run to it, so an existing ``target_dir`` is replaced only with ``replace`` and then
    moved to ``<parent>/.replaced/<name>.<UTC timestamp>`` first. Raises SchemaInUse
    while a run's chewBBACA step uses it (:func:`schema_in_use`).
    """
    built_dir, target_dir = Path(built_dir), Path(target_dir)
    target_dir.parent.mkdir(parents=True, exist_ok=True)
    if not target_dir.exists():
        built_dir.replace(target_dir)
        return None
    if not replace:
        raise FileExistsError(f"schema already exists: {target_dir}")
    fd = _open_lock(target_dir / SCHEMA_USE_LOCK)
    try:
        if fd is not None:
            try:
                fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                raise SchemaInUse(f"{target_dir} is in use by a running chewBBACA step; "
                                  "replace it when the run is done") from None
        backup_root = target_dir.parent / REPLACED_DIR
        backup_root.mkdir(exist_ok=True)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        backup = backup_root / f"{target_dir.name}.{stamp}"
        n = 2
        while backup.exists():
            backup = backup_root / f"{target_dir.name}.{stamp}_{n}"
            n += 1
        target_dir.replace(backup)
        try:
            built_dir.replace(target_dir)
        except BaseException:
            backup.replace(target_dir)  # put the previous schema back
            raise
    finally:
        if fd is not None:
            os.close(fd)
    return backup


def _new_temp(target_dir: Path) -> Path:
    target_dir.parent.mkdir(parents=True, exist_ok=True)
    return Path(tempfile.mkdtemp(prefix=f".{target_dir.name}_", dir=str(target_dir.parent)))


def _finalize(
    stager: SchemeStager,
    expected_order: list[str],
    target_dir: Path,
    temp_root: Path,
    *,
    source: Path,
    source_format: str,
    warning: str,
    options: BuildOptions,
    progress: Progress,
) -> SchemaBuildResult:
    """Adapt the staged loci with PrepExternalSchema, validate, then move it to ``target_dir``."""
    if not stager.written_names:
        raise ValueError("no loci were produced from the source")
    adapted = temp_root / target_dir.name
    report_dir = run_prep_external_schema(
        stager.temp_dir, adapted, cpu=options.cpu, base_cmd=options.base_cmd,
        translation_table=options.translation_table, progress=progress, should_cancel=options.should_cancel,
    )
    _emit(progress, "Validating the adapted schema…")
    metadata = read_schema_metadata(adapted, python_cmd=options.python_cmd)
    issues = adapted_schema_issues(adapted, metadata=metadata)
    if issues:
        raise RuntimeError("adapted schema failed validation: " + "; ".join(issues))

    stats = schema_stats(adapted, orientation=False, metadata=metadata)
    rejected_alleles = _nonempty_lines(report_dir / INVALID_ALLELES_FILE)
    rejected_loci = _nonempty_lines(report_dir / INVALID_LOCI_FILE)
    log_text = (report_dir / "prep_external_schema.log").read_text(encoding="utf-8")
    expected = expected_order or list(stager.written_names)
    missing = [n for n in expected if n not in stager.written_set]
    dup = sorted(set(stager.duplicate_names))
    parts = [warning] if warning else []
    if dup:
        parts.append(f"ignored {len(dup)} duplicate locus file(s)")
    if rejected_loci:
        parts.append(f"{len(rejected_loci)} locus/loci had no valid CDS allele and were dropped")
    if rejected_alleles:
        parts.append(f"{len(rejected_alleles)} allele(s) were not a valid CDS and were dropped")
    build_info = {
        "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "builder": "chewBBACA PrepExternalSchema",
        "chewbbaca_version": _chewbbaca_version(log_text),
        "bsr": float(PREP_BSR),
        "translation_table": int(options.translation_table),
        "source": str(source),
        "source_format": source_format,
        "expected_loci": len(expected),
        "source_loci": len(stager.written_names),
        "source_alleles": stager.allele_count,
        "missing_source_loci": missing,
        "loci": stats.loci,
        "alleles": stats.alleles,
        "representatives": stats.representatives,
        "rejected_alleles": len(rejected_alleles),
        "rejected_loci": rejected_loci,
    }
    (report_dir / BUILD_INFO_FILE).write_text(json.dumps(build_info, indent=2) + "\n", encoding="utf-8")
    install_schema_dir(adapted, target_dir)
    shutil.rmtree(temp_root, ignore_errors=True)
    return SchemaBuildResult(
        target_dir=target_dir,
        schema_id=target_dir.name,
        written_loci=stats.loci,
        expected_loci=len(expected),
        missing_loci=len(missing) + len(rejected_loci),
        source_format=source_format,
        warning="; ".join(parts),
        source_loci=len(stager.written_names),
        source_alleles=stager.allele_count,
        alleles=stats.alleles,
        representatives=stats.representatives,
        rejected_alleles=len(rejected_alleles),
        rejected_loci=rejected_loci,
    )


def _copy_prebuilt(source: Path, target_dir: Path, metadata: SchemaMetadata) -> SchemaBuildResult:
    """Copy a chewBBACA-made schema (without its runtime caches) and record its build."""
    tmp = _new_temp(target_dir)
    shutil.rmtree(tmp)  # copytree needs a non-existent dest
    try:
        shutil.copytree(source, tmp, ignore=shutil.ignore_patterns(*[Path(c).name for c in CHEWBBACA_RUNTIME_CACHES],
                                                                 SCHEMA_USE_LOCK))
        stats = schema_stats(tmp, orientation=False, metadata=metadata)
        previous = read_build_info(source)
        build_info = {
            **previous,
            "imported_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "builder": previous.get("builder") or "chewBBACA (prebuilt schema)",
            "translation_table": int(metadata.translation_table or PREP_TRANSLATION_TABLE),
            "source": str(source),
            "source_format": "prebuilt_schema",
            "loci": stats.loci,
            "alleles": stats.alleles,
            "representatives": stats.representatives,
        }
        build_info.setdefault("built_at", build_info["imported_at"])
        report_dir = tmp / SCHEMA_BUILD_REPORT_DIR
        report_dir.mkdir(exist_ok=True)
        (report_dir / BUILD_INFO_FILE).write_text(json.dumps(build_info, indent=2) + "\n", encoding="utf-8")
        install_schema_dir(tmp, target_dir)
    except Exception:
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    return SchemaBuildResult(
        target_dir, target_dir.name, stats.loci, stats.loci, 0, "prebuilt_schema", "",
        source_loci=stats.loci, source_alleles=stats.alleles, alleles=stats.alleles,
        representatives=stats.representatives,
    )


def build_schema(
    source: Path,
    target_dir: Path,
    *,
    single_fasta_mode: str = "loci_set",  # "loci_set" | "single_locus"
    progress: Progress = None,
    targets_csv: Path | None = None,
    options: BuildOptions | None = None,
) -> SchemaBuildResult:
    """Build a chewBBACA schema at ``target_dir`` from ``source``.

    ``source`` may be a directory (prebuilt schema, or a folder of per-locus /
    combined FASTAs) or a file (.zip, or a single FASTA). ``single_fasta_mode``
    disambiguates a lone FASTA: a multi-locus set vs. all alleles of one locus.
    ``targets_csv`` (a cgMLST.org locus manifest) restricts a zip source to the
    scheme's loci. Every non-prebuilt source is adapted with PrepExternalSchema.
    ``target_dir`` must not exist yet (see :func:`install_schema_dir` to replace one).
    """
    source = Path(source).expanduser()
    target_dir = Path(target_dir).expanduser()
    options = options or BuildOptions()
    if not source.exists():
        raise FileNotFoundError(f"source not found: {source}")
    if target_dir.exists():
        raise FileExistsError(f"schema already exists: {target_dir}")

    # 1) already-built chewBBACA schema directory: copy it when chewBBACA made it (its
    #    representatives are forward CDS), otherwise re-adapt its locus FASTAs (an old
    #    hand-built schema). Its pickles are read in the tools image.
    if source.is_dir() and is_prebuilt_schema(source):
        metadata = (
            read_schema_metadata(source, python_cmd=options.python_cmd)
            if (source / CHEWBBACA_SCHEMA_CONFIG).is_file() else None
        )
        stats = schema_stats(source, metadata=metadata) if metadata else None
        if stats and stats.genes_list == stats.loci and not (
            stats.reps_reverse_only or stats.reps_invalid or stats.loci_without_short
        ):
            _emit(progress, "Copying prebuilt chewBBACA schema…")
            return _copy_prebuilt(source, target_dir, metadata)
        _emit(progress, "Prebuilt schema is not in coding orientation; re-adapting its loci…")
        return _build_from_folder(source, target_dir, progress=progress,
                                  options=options, source_format="prebuilt_schema_readapted")

    # 2) directory of loci FASTAs (Ridom)
    if source.is_dir():
        return _build_from_folder(source, target_dir, progress=progress, options=options)

    # 3) single file
    fmt = detect_file_format(source)
    if fmt == "html_error":
        raise ValueError("source is an HTML error page, not a FASTA/zip")
    temp_root = _new_temp(target_dir)
    common = dict(source=source, options=options, progress=progress)
    try:
        external = temp_root / "external"
        if fmt in ("zip_archive", "zip_stream"):
            _emit(progress, "Reading zip archive…")
            expected = read_expected_loci(Path(targets_csv)) if targets_csv else []
            stager = SchemeStager(external, expected or None)
            sfmt = convert_zip_archive(source, stager, progress) if fmt == "zip_archive" \
                else convert_zip_stream(source, stager)
            return _finalize(stager, expected or list(stager.written_names), target_dir, temp_root,
                             source_format=sfmt + ("+targets" if expected else ""), warning="",
                             **common)
        # a lone FASTA
        if single_fasta_mode == "single_locus":
            _emit(progress, "Building single-locus schema…")
            stager = SchemeStager(external, None)
            locus_name = f"{slugify_species(source.stem, fallback='locus')}.fasta"
            stager.write_locus(locus_name, source.read_bytes())
            return _finalize(stager, [locus_name], target_dir, temp_root,
                             source_format="single_locus_fasta", warning="", **common)
        _emit(progress, "Splitting combined FASTA into loci…")
        stager = SchemeStager(external, None)
        convert_combined_fasta(source, stager, [])
        return _finalize(stager, list(stager.written_names), target_dir, temp_root,
                         source_format="combined_fasta", warning="", **common)
    except Exception:
        shutil.rmtree(temp_root, ignore_errors=True)
        raise


def _build_from_folder(
    folder: Path,
    target_dir: Path,
    *,
    progress: Progress,
    options: BuildOptions,
    source_format: str = "",
) -> SchemaBuildResult:
    fastas = _list_fastas(folder)
    if not fastas:
        raise ValueError(f"no FASTA files found in {folder}")

    targets_csv = _find_targets_csv(folder)
    temp_root = _new_temp(target_dir)
    external = temp_root / "external"
    common = dict(source=folder, options=options, progress=progress)
    try:
        if targets_csv is not None:
            # Manifest-driven: only files named in targets.csv become loci; anything
            # else (e.g. a stray combined all-loci FASTA) is skipped by the stager.
            expected = read_expected_loci(targets_csv)
            _emit(progress, f"Staging {len(expected)} loci (targets.csv)…")
            stager = SchemeStager(external, expected)
            for idx, f in enumerate(fastas):
                stager.write_locus(f.name, f.read_bytes())
                if idx % 200 == 0:
                    _emit(progress, f"Staging loci… {idx + 1}/{len(fastas)}")
            return _finalize(stager, expected, target_dir, temp_root,
                             source_format=source_format or "folder+targets", warning="", **common)

        # No manifest: split per-locus files from any combined all-loci file(s).
        combined = [f for f in fastas if is_combined_fasta(f)]
        per_locus = [f for f in fastas if f not in combined]
        stager = SchemeStager(external, None)
        if per_locus:
            _emit(progress, f"Staging {len(per_locus)} loci…")
            for idx, f in enumerate(per_locus):
                stager.write_locus(f.name, f.read_bytes())
                if idx % 200 == 0:
                    _emit(progress, f"Staging loci… {idx + 1}/{len(per_locus)}")
            warning = f"ignored {len(combined)} combined-style FASTA file(s)" if combined else ""
            return _finalize(stager, [f.name for f in per_locus], target_dir, temp_root,
                             source_format=source_format or "folder", warning=warning, **common)
        if combined:
            _emit(progress, "Splitting combined FASTA into loci…")
            for f in combined:
                convert_combined_fasta(f, stager, [])
            return _finalize(stager, list(stager.written_names), target_dir, temp_root,
                             source_format=source_format or "folder_combined", warning="", **common)
        raise ValueError(f"no usable loci FASTAs found in {folder}")
    except Exception:
        shutil.rmtree(temp_root, ignore_errors=True)
        raise


# ---------------------------------------------------------------------------
# Registration
# ---------------------------------------------------------------------------
def _atomic_write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(json.dumps(payload, indent=2, sort_keys=False) + "\n")
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


@contextmanager
def _catalog_lock(catalog_path: Path) -> Iterator[None]:
    """Serialize catalog and config updates across processes (a GUI import and a run's
    auto-install) with a lock file next to the catalog."""
    catalog_path.parent.mkdir(parents=True, exist_ok=True)
    fd = _open_lock(catalog_path.with_name(f".{catalog_path.name}.lock"))
    try:
        if fd is not None:
            fcntl.flock(fd, fcntl.LOCK_EX)
        yield
    finally:
        if fd is not None:
            os.close(fd)


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

    The other entries are kept as they are (a catalog in the bare-list format is kept
    and written in the ``{"version", "schemes"}`` format). Raises ValueError on an
    MLST-scheme alias already claimed by another schema, RuntimeError when the catalog
    cannot be read or the config cannot be updated (the catalog entry is then written).
    """
    cgmlst_dir = Path(cgmlst_dir)
    catalog_path = Path(catalog_path) if catalog_path else cgmlst_dir / CHEWBBACA_SCHEME_CATALOG

    with _catalog_lock(catalog_path):
        data: dict = {"version": 1, "schemes": []}
        if catalog_path.is_file():
            try:
                loaded = json.loads(catalog_path.read_text(encoding="utf-8"))
            except Exception as exc:
                raise RuntimeError(f"cannot read the cgMLST catalog {catalog_path}: {exc}") from exc
            if isinstance(loaded, dict):
                data = loaded
            elif isinstance(loaded, list):
                data = {"version": 1, "schemes": loaded}
            else:
                raise RuntimeError(f"cgMLST catalog {catalog_path} is neither a list nor an object")
        schemes = data.get("schemes")
        if not isinstance(schemes, list):
            raise RuntimeError(f"cgMLST catalog {catalog_path} has no list of schemes")

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
        schemes.sort(key=lambda e: str(e.get("schema_id")) if isinstance(e, dict) else "")
        data["schemes"] = schemes
        data.setdefault("version", 1)
        _atomic_write_json(catalog_path, data)

        if update_config:
            from pipeline.ToolEnv import save_app_config

            try:
                save_app_config({
                    "chewbbaca": {
                        "available_schemas": {schema_id: str((cgmlst_dir / schema_id).resolve())}
                    }
                })
            except Exception as exc:
                raise RuntimeError(f"added to the catalog, but the config could not be updated: {exc}") from exc
    return catalog_path
