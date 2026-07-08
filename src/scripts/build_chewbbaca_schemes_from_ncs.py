#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import pickle
import re
import shutil
import struct
import sys
import tempfile
import zipfile
import zlib
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Iterable


DEFAULT_INPUT_DIR = Path("/Data/cgmlst/download_tmp")
DEFAULT_OUTPUT_DIR = Path("/Data/cgmlst/schemes")
CHEWBBACA_SHORT_DIR = "short"
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


@dataclass
class SchemeInput:
    slug: str
    title: str
    input_dir: Path
    targets_csv: Path
    alleles_path: Path
    output_name: str


@dataclass
class ConversionResult:
    slug: str
    title: str
    output_name: str
    output_dir: str
    source_format: str
    status: str
    expected_loci: int
    written_loci: int
    missing_loci: int
    unexpected_loci: int
    error: str = ""
    warning: str = ""
    missing_preview: list[str] = field(default_factory=list)
    unexpected_preview: list[str] = field(default_factory=list)


def announce(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def slugify_title(title: str, fallback_slug: str) -> str:
    label = re.sub(r"\bcgmlst\b", "", title, flags=re.IGNORECASE)
    label = label.replace("&", " and ")
    label = label.replace("/", " ")
    label = label.replace("(", " ")
    label = label.replace(")", " ")
    label = re.sub(r"[^A-Za-z0-9]+", "_", label).strip("_").lower()
    return label or fallback_slug.lower()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert cgMLST.org NCS downloads into chewBBACA schema directories.",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=f"Directory containing per-scheme NCS downloads (default: {DEFAULT_INPUT_DIR})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory where chewBBACA schemas will be created (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--match",
        action="append",
        default=[],
        help="Case-insensitive substring filter for scheme slug/title; may be passed multiple times.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing output scheme directory if it already exists.",
    )
    parser.add_argument(
        "--allow-partial",
        action="store_true",
        help="Write partial schemes when the source archive is truncated or incomplete.",
    )
    parser.add_argument(
        "--report-only",
        action="store_true",
        help="Inspect inputs and write a conversion report without creating schema directories.",
    )
    return parser.parse_args()


def load_scheme_titles(input_root: Path) -> dict[str, str]:
    manifest = input_root / "schemes.json"
    if not manifest.is_file():
        return {}
    data = json.loads(manifest.read_text(encoding="utf-8"))
    out: dict[str, str] = {}
    if isinstance(data, list):
        for row in data:
            if not isinstance(row, dict):
                continue
            slug = str(row.get("slug") or "").strip()
            title = str(row.get("title") or "").strip()
            if slug:
                out[slug] = title or slug
    return out


def find_alleles_source(scheme_dir: Path) -> Path | None:
    for name in (
        "alleles.zip",
        "alleles.fasta",
        "alleles.truncated.zip",
        "alleles.html",
    ):
        path = scheme_dir / name
        if path.is_file():
            return path
    return None


def discover_scheme_inputs(input_root: Path, match_terms: list[str]) -> list[SchemeInput]:
    titles = load_scheme_titles(input_root)
    schemes: list[SchemeInput] = []
    for child in sorted(input_root.iterdir(), key=lambda path: path.name.lower()):
        if not child.is_dir():
            continue
        targets_csv = child / "targets.csv"
        alleles_path = find_alleles_source(child)
        if not targets_csv.is_file() or alleles_path is None:
            continue
        slug = child.name
        title = titles.get(slug, slug)
        haystack = " ".join([slug.lower(), title.lower()])
        if match_terms and not all(term in haystack for term in match_terms):
            continue
        schemes.append(
            SchemeInput(
                slug=slug,
                title=title,
                input_dir=child,
                targets_csv=targets_csv,
                alleles_path=alleles_path,
                output_name=slugify_title(title, slug),
            )
        )
    return schemes


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


def looks_like_html_error(path: Path) -> bool:
    with path.open("rb") as fh:
        prefix = fh.read(256).lstrip().lower()
    return any(prefix.startswith(marker) for marker in HTML_PREFIXES)


def wrapped_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width))


def first_fasta_record(data: bytes, *, source_name: str) -> tuple[str, str]:
    header: str | None = None
    sequence_parts: list[str] = []
    text = data.decode("utf-8", errors="replace")
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                break
            header = line[1:].strip()
            continue
        if header is None:
            continue
        sequence_parts.append(line)
    if not header or not sequence_parts:
        raise ValueError(f"{source_name} does not contain a readable FASTA record")
    return header, "".join(sequence_parts)


def parse_fasta_records(data: bytes, *, source_name: str) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence_parts: list[str] = []
    text = data.decode("utf-8", errors="replace")
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None and sequence_parts:
                records.append((header, "".join(sequence_parts)))
            header = line[1:].strip()
            sequence_parts = []
            continue
        if header is None:
            continue
        sequence_parts.append(line)
    if header is not None and sequence_parts:
        records.append((header, "".join(sequence_parts)))
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


def normalize_locus_fasta(locus_name: str, locus_data: bytes) -> bytes:
    records = parse_fasta_records(locus_data, source_name=locus_name)
    seen_headers: set[str] = set()
    out_lines: list[str] = []
    for raw_header, sequence in records:
        header = normalized_allele_header(locus_name, raw_header, seen_headers)
        out_lines.append(f">{header}")
        out_lines.append(wrapped_sequence(sequence))
    return ("\n".join(out_lines) + "\n").encode("utf-8")


def write_short_fasta(short_path: Path, locus_data: bytes, *, source_name: str) -> None:
    header, sequence = first_fasta_record(locus_data, source_name=source_name)
    short_path.write_text(
        f">{header}\n{wrapped_sequence(sequence)}\n",
        encoding="utf-8",
    )


class SchemeStager:
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
        target = self.temp_dir / locus_name
        normalized_data = normalize_locus_fasta(locus_name, locus_data)
        target.write_bytes(normalized_data)
        short_target = self.short_dir / f"{Path(locus_name).stem}_short.fasta"
        write_short_fasta(short_target, normalized_data, source_name=locus_name)
        self.written_set.add(locus_name)
        self.written_names.append(locus_name)

    def finalize_schema_metadata(self, expected_order: list[str]) -> None:
        ordered = [name for name in expected_order if name in self.written_set]
        if not ordered:
            ordered = list(self.written_names)
        with (self.temp_dir / ".schema_config").open("wb") as fh:
            pickle.dump(SCHEMA_CONFIG, fh, protocol=pickle.HIGHEST_PROTOCOL)
        with (self.temp_dir / ".genes_list").open("wb") as fh:
            pickle.dump(ordered, fh, protocol=pickle.HIGHEST_PROTOCOL)


def convert_zip_archive(alleles_path: Path, stager: SchemeStager) -> tuple[str, str]:
    with zipfile.ZipFile(alleles_path) as archive:
        for info in archive.infolist():
            if info.is_dir():
                continue
            stager.write_locus(Path(info.filename).name, archive.read(info))
    return "zip_archive", ""


def _stream_zip_members(
    alleles_path: Path,
    on_member: Callable[[str, bytes], None],
) -> str:
    data = alleles_path.read_bytes()
    pos = 0
    while pos + 30 <= len(data) and data[pos:pos + 4] == ZIP_LOCAL_FILE_SIG:
        header = data[pos + 4:pos + 30]
        (
            _version,
            _flags,
            compression,
            _mtime,
            _mdate,
            _crc32,
            _compressed_size,
            _uncompressed_size,
            name_len,
            extra_len,
        ) = struct.unpack("<HHHHHIIIHH", header)
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
                consumed = len(chunk) - len(decompressor.unused_data)
                read_pos += consumed
                break
            read_pos += len(chunk)

        on_member(Path(name).name, b"".join(chunks))
        if data[read_pos:read_pos + 4] == ZIP_DATA_DESCRIPTOR_SIG:
            pos = read_pos + 16
        else:
            pos = read_pos + 12

    trailing = len(data) - pos
    if trailing not in (0,):
        return f"ZIP stream ended with {trailing} trailing byte(s)"
    return ""


def convert_zip_stream(alleles_path: Path, stager: SchemeStager) -> tuple[str, str]:
    warning = _stream_zip_members(alleles_path, stager.write_locus)
    return "zip_stream", warning


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
    alleles_path: Path,
    stager: SchemeStager,
    expected_names: list[str],
) -> tuple[str, str]:
    groups: dict[str, list[str]] = {}
    current_header: str | None = None
    current_seq: list[str] = []

    def flush_current() -> None:
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
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush_current()
                current_header = line[1:].strip()
                continue
            if current_header is not None:
                current_seq.append(line)
    flush_current()

    for locus_name, records in groups.items():
        stager.write_locus(locus_name, "".join(records).encode("utf-8"))
    return "combined_fasta", ""


def detect_source_format(alleles_path: Path) -> str:
    if looks_like_html_error(alleles_path):
        return "html_error"
    if zipfile.is_zipfile(alleles_path):
        return "zip_archive"
    with alleles_path.open("rb") as fh:
        prefix = fh.read(8)
    if prefix.startswith(ZIP_LOCAL_FILE_SIG):
        return "zip_stream"
    if prefix.lstrip().startswith(b">"):
        return "combined_fasta"
    return "unknown"


def convert_source(
    alleles_path: Path,
    stager: SchemeStager,
    expected_names: list[str],
) -> tuple[str, str]:
    source_format = detect_source_format(alleles_path)
    if source_format == "html_error":
        preview = alleles_path.read_text(encoding="utf-8", errors="replace").strip()
        raise ValueError(f"HTML error page instead of FASTA download: {preview[:160]}")
    if source_format == "zip_archive":
        return convert_zip_archive(alleles_path, stager)
    if source_format == "zip_stream":
        return convert_zip_stream(alleles_path, stager)
    if source_format == "combined_fasta":
        return convert_combined_fasta(alleles_path, stager, expected_names)

    raise ValueError("unsupported alleles source format")


def replace_output_dir(target_dir: Path, temp_dir: Path, *, overwrite: bool) -> None:
    if target_dir.exists():
        if not overwrite:
            raise FileExistsError(f"output directory already exists: {target_dir}")
        shutil.rmtree(target_dir)
    target_dir.parent.mkdir(parents=True, exist_ok=True)
    temp_dir.replace(target_dir)


def skipped_result(
    scheme: SchemeInput,
    *,
    target_dir: Path,
    source_format: str,
    expected_loci: int,
    written_loci: int,
    unexpected_loci: int,
    error: str,
    warning: str = "",
    missing_loci: int | None = None,
    missing_preview: list[str] | None = None,
    unexpected_preview: list[str] | None = None,
) -> ConversionResult:
    if missing_loci is None:
        missing_loci = max(0, expected_loci - written_loci)
    return ConversionResult(
        slug=scheme.slug,
        title=scheme.title,
        output_name=scheme.output_name,
        output_dir=str(target_dir),
        source_format=source_format,
        status="skipped",
        expected_loci=expected_loci,
        written_loci=written_loci,
        missing_loci=missing_loci,
        unexpected_loci=unexpected_loci,
        error=error,
        warning=warning,
        missing_preview=missing_preview or [],
        unexpected_preview=unexpected_preview or [],
    )


def convert_scheme(
    scheme: SchemeInput,
    output_root: Path,
    *,
    overwrite: bool,
    allow_partial: bool,
    report_only: bool,
) -> ConversionResult:
    expected_names = read_expected_loci(scheme.targets_csv)
    target_dir = output_root / scheme.output_name
    if target_dir.exists() and not overwrite and not report_only:
        return skipped_result(
            scheme,
            target_dir=target_dir,
            source_format="unknown",
            expected_loci=len(expected_names),
            written_loci=0,
            unexpected_loci=0,
            error=f"output directory already exists: {target_dir}",
        )

    temp_root = Path(tempfile.mkdtemp(prefix=f"{scheme.output_name}_", dir=str(output_root)))
    warning = ""
    source_format = detect_source_format(scheme.alleles_path)
    try:
        stager = SchemeStager(temp_root, expected_names)
        try:
            source_format, warning = convert_source(scheme.alleles_path, stager, expected_names)
        except Exception as exc:
            missing_preview = expected_names[:10]
            return skipped_result(
                scheme,
                target_dir=target_dir,
                source_format=source_format,
                expected_loci=len(expected_names),
                written_loci=len(stager.written_names),
                unexpected_loci=len(stager.unexpected_names),
                error=str(exc),
                warning=warning,
                missing_preview=missing_preview,
                unexpected_preview=stager.unexpected_names[:10],
            )

        missing = [name for name in expected_names if name not in stager.written_set]
        unexpected = sorted(set(stager.unexpected_names))
        duplicate = sorted(set(stager.duplicate_names))
        warning_parts = [part for part in [warning] if part]
        if duplicate:
            warning_parts.append(f"ignored duplicate loci: {', '.join(duplicate[:5])}")
        warning_text = "; ".join(warning_parts)

        if missing and not allow_partial:
            return skipped_result(
                scheme,
                target_dir=target_dir,
                source_format=source_format,
                expected_loci=len(expected_names),
                written_loci=len(stager.written_names),
                missing_loci=len(missing),
                unexpected_loci=len(unexpected),
                error="source download is incomplete; skipping this species",
                warning=warning_text,
                missing_preview=missing[:10],
                unexpected_preview=unexpected[:10],
            )

        if not report_only:
            stager.finalize_schema_metadata(expected_names)
            replace_output_dir(target_dir, temp_root, overwrite=overwrite)
        status = "partial" if missing else "converted"
        return ConversionResult(
            slug=scheme.slug,
            title=scheme.title,
            output_name=scheme.output_name,
            output_dir=str(target_dir),
            source_format=source_format,
            status=status,
            expected_loci=len(expected_names),
            written_loci=len(stager.written_names),
            missing_loci=len(missing),
            unexpected_loci=len(unexpected),
            warning=warning_text,
            missing_preview=missing[:10],
            unexpected_preview=unexpected[:10],
        )
    finally:
        if temp_root.exists():
            shutil.rmtree(temp_root, ignore_errors=True)


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")


def write_tsv(path: Path, rows: list[dict[str, object]], field_order: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(field_order)
        for row in rows:
            writer.writerow([row.get(field, "") for field in field_order])


def main() -> int:
    args = parse_args()
    input_root = args.input_dir.expanduser()
    output_root = args.output_dir.expanduser()
    match_terms = [term.strip().lower() for term in args.match if term.strip()]

    announce(f"Input directory: {input_root}")
    announce(f"Output directory: {output_root}")
    if match_terms:
        announce(f"Filtering schemes with: {', '.join(match_terms)}")

    if not input_root.exists() or not input_root.is_dir():
        announce(f"Input directory not found: {input_root}")
        return 1

    output_root.mkdir(parents=True, exist_ok=True)
    schemes = discover_scheme_inputs(input_root, match_terms)
    if not schemes:
        announce("No matching downloaded schemes were found.")
        return 1

    announce(f"Discovered {len(schemes)} downloaded scheme directories.")
    if args.report_only:
        announce("Report-only mode enabled. No schema directories will be written.")

    results: list[ConversionResult] = []
    for idx, scheme in enumerate(schemes, start=1):
        announce(f"[{idx}/{len(schemes)}] Converting {scheme.slug} -> {scheme.output_name}")
        result = convert_scheme(
            scheme,
            output_root,
            overwrite=args.overwrite,
            allow_partial=args.allow_partial,
            report_only=args.report_only,
        )
        results.append(result)
        details = (
            f"status={result.status}, format={result.source_format}, "
            f"written={result.written_loci}/{result.expected_loci}"
        )
        if result.warning:
            details += f", warning={result.warning}"
        if result.error:
            details += f", error={result.error}"
        announce(f"  {details}")

    report_rows = [asdict(result) for result in results]
    announce("Writing conversion reports.")
    write_json(output_root / "conversion_report.json", report_rows)
    write_tsv(
        output_root / "conversion_report.tsv",
        report_rows,
        [
            "slug",
            "title",
            "output_name",
            "output_dir",
            "source_format",
            "status",
            "expected_loci",
            "written_loci",
            "missing_loci",
            "unexpected_loci",
            "error",
            "warning",
        ],
    )

    partial = sum(1 for result in results if result.status == "partial")
    converted = sum(1 for result in results if result.status == "converted")
    skipped = sum(1 for result in results if result.status == "skipped")

    announce(
        f"Summary: converted={converted}, partial={partial}, skipped={skipped}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
