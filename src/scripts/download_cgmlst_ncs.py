#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from dataclasses import asdict, dataclass
from html import unescape
from html.parser import HTMLParser
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Iterable, Optional


INDEX_URL = "https://www.cgmlst.org/ncs/"
DEFAULT_OUTPUT_DIR = Path("/Data/cgmlst/download_tmp")
DEFAULT_TIMEOUT = 120
DEFAULT_RETRIES = 5
DEFAULT_WORKERS = 4
USER_AGENT = "BACTolomew-cgMLST-NCS-downloader/1.0"
HTML_PREFIXES = (b"<html", b"<!doctype html")


def announce(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def _normalize_space(text: str) -> str:
    return re.sub(r"\s+", " ", unescape(text or "")).strip()


def _scheme_slug_from_url(url: str) -> str:
    parsed = urllib.parse.urlparse(url)
    slug = parsed.path.rstrip("/").split("/")[-1]
    return slug or "scheme"


def _request(url: str, *, timeout: int) -> urllib.request.Request:
    return urllib.request.Request(
        url,
        headers={
            "User-Agent": USER_AGENT,
            "Accept": "*/*",
        },
    )


def _fetch_text(url: str, *, timeout: int, retries: int) -> str:
    last_exc: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(_request(url, timeout=timeout), timeout=timeout) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except Exception as exc:  # urllib raises a mix of exception types
            last_exc = exc
            if attempt == retries:
                break
            time.sleep(min(2 ** (attempt - 1), 5))
    raise RuntimeError(f"Failed to fetch {url}: {last_exc}")


class LinkParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[tuple[str, str]] = []
        self._href: Optional[str] = None
        self._text_parts: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, Optional[str]]]) -> None:
        if tag.lower() != "a":
            return
        attr_map = dict(attrs)
        href = attr_map.get("href")
        if href is None:
            return
        self._href = href
        self._text_parts = []

    def handle_data(self, data: str) -> None:
        if self._href is None:
            return
        self._text_parts.append(data)

    def handle_endtag(self, tag: str) -> None:
        if tag.lower() != "a" or self._href is None:
            return
        text = _normalize_space("".join(self._text_parts))
        self.links.append((self._href, text))
        self._href = None
        self._text_parts = []


def _parse_links(html_text: str) -> list[tuple[str, str]]:
    parser = LinkParser()
    parser.feed(html_text)
    return parser.links


@dataclass
class Scheme:
    title: str
    scheme_url: str
    slug: str
    locus_url: str
    csv_url: str
    fasta_url: str


@dataclass
class DownloadTask:
    scheme_slug: str
    scheme_title: str
    kind: str
    source_url: str
    output_path: Path


@dataclass
class DownloadResult:
    scheme_slug: str
    scheme_title: str
    kind: str
    source_url: str
    output_path: str
    status: str
    attempts: int = 0
    error: str = ""


def _fasta_candidate_paths(directory: Path) -> list[Path]:
    return [
        directory / "alleles.zip",
        directory / "alleles.fasta",
        directory / "alleles.truncated.zip",
        directory / "alleles.html",
    ]


def _classify_fasta_path(path: Path) -> str:
    if _looks_like_html(path):
        return "html"

    with path.open("rb") as fh:
        prefix = fh.read(8)

    if zipfile.is_zipfile(path):
        return "zip"
    if prefix.startswith(b"PK\x03\x04"):
        return "truncated_zip_like"
    if prefix.lstrip().startswith(b">"):
        return "fasta"
    return "unknown"


def _fasta_filename_for_kind(kind: str) -> str:
    if kind == "zip":
        return "alleles.zip"
    if kind == "fasta":
        return "alleles.fasta"
    if kind == "html":
        return "alleles.html"
    if kind == "truncated_zip_like":
        return "alleles.truncated.zip"
    return "alleles.bin"


def _looks_like_html(path: Path) -> bool:
    with path.open("rb") as fh:
        prefix = fh.read(256).lstrip().lower()
    return any(prefix.startswith(marker) for marker in HTML_PREFIXES)


def _validate_csv_file(path: Path) -> str | None:
    if not path.exists() or not path.is_file():
        return "CSV file is missing"
    if path.stat().st_size == 0:
        return "CSV file is empty"
    if _looks_like_html(path):
        return "received an HTML error page instead of CSV data"

    try:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader, None)
            if not header:
                return "CSV file has no header row"

            header_fields = {field.strip() for field in header if field.strip()}
            if "Gene" not in header_fields or "Alleles" not in header_fields:
                return "CSV header does not look like an NCS targets table"

            data_rows = 0
            fasta_refs = 0
            for row in reader:
                if any(field.strip() for field in row):
                    data_rows += 1
                if any(field.strip().lower().endswith(".fasta") for field in row):
                    fasta_refs += 1
    except Exception as exc:
        return f"CSV validation failed: {exc}"

    if data_rows == 0:
        return "CSV file has no data rows"
    if fasta_refs == 0:
        return "CSV file contains no locus FASTA references"
    return None


def _expected_fasta_names_from_csv(path: Path) -> set[str] | None:
    validation_error = _validate_csv_file(path)
    if validation_error is not None:
        return None

    expected: set[str] = set()
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            locus_name = next(
                (field.strip() for field in reversed(row) if field.strip().lower().endswith(".fasta")),
                "",
            )
            if locus_name:
                expected.add(locus_name)
    return expected or None


def _validate_plain_fasta_file(path: Path) -> str | None:
    header_seen = False
    sequence_seen = False
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            for raw_line in fh:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    header_seen = True
                    continue
                if header_seen:
                    sequence_seen = True
                    break
    except Exception as exc:
        return f"plain FASTA validation failed: {exc}"

    if not header_seen:
        return "plain FASTA file contains no header lines"
    if not sequence_seen:
        return "plain FASTA file contains no sequence data"
    return None


def _validate_fasta_file(path: Path) -> str | None:
    if not path.exists() or not path.is_file():
        return "FASTA file is missing"
    if path.stat().st_size == 0:
        return "FASTA file is empty"
    if _looks_like_html(path):
        return "received an HTML error page instead of FASTA data"

    classification = _classify_fasta_path(path)

    if classification == "zip":
        try:
            with zipfile.ZipFile(path) as archive:
                infos = [info for info in archive.infolist() if not info.is_dir()]
                if not infos:
                    return "ZIP archive contains no files"
                if not all(info.filename.lower().endswith(".fasta") for info in infos):
                    return "ZIP archive contains non-FASTA members"

                bad_member = archive.testzip()
                if bad_member is not None:
                    return f"ZIP archive member failed CRC validation: {bad_member}"

                sample = archive.read(infos[0]).decode("utf-8", errors="replace")
                if ">" not in sample:
                    return "ZIP archive members do not look like FASTA files"

                expected_names = _expected_fasta_names_from_csv(path.parent / "targets.csv")
                if expected_names is not None:
                    archive_names = {Path(info.filename).name for info in infos}
                    missing = expected_names - archive_names
                    extras = archive_names - expected_names
                    if missing:
                        return f"ZIP archive is missing {len(missing)} loci expected from targets.csv"
                    if extras:
                        return f"ZIP archive contains {len(extras)} unexpected loci not listed in targets.csv"
        except Exception as exc:
            return f"ZIP FASTA validation failed: {exc}"
        return None

    if classification == "truncated_zip_like":
        return "file looks like a truncated ZIP archive"

    if classification == "html":
        return "received an HTML error page instead of FASTA data"

    return _validate_plain_fasta_file(path)


def _validate_download(task: DownloadTask, path: Path) -> str | None:
    if task.kind == "csv":
        return _validate_csv_file(path)
    if task.kind == "fasta":
        return _validate_fasta_file(path)
    return f"unknown download kind: {task.kind}"


def discover_schemes(index_url: str, *, timeout: int, retries: int, match_terms: list[str]) -> list[Scheme]:
    announce(f"Fetching scheme index: {index_url}")
    index_html = _fetch_text(index_url, timeout=timeout, retries=retries)
    links = _parse_links(index_html)
    candidate_links = [
        (href, text)
        for href, text in links
        if text and "cgMLST" in text and "/ncs/schema/" in urllib.parse.urljoin(index_url, href)
    ]
    announce(f"Found {len(candidate_links)} scheme candidates on the index page.")
    schemes: list[Scheme] = []
    seen: set[str] = set()

    for href, text in candidate_links:
        full_url = urllib.parse.urljoin(index_url, href)
        if full_url in seen:
            continue

        slug = _scheme_slug_from_url(full_url)
        haystack = " ".join([text.lower(), slug.lower(), full_url.lower()])
        if match_terms and not all(term in haystack for term in match_terms):
            continue

        locus_url = urllib.parse.urljoin(full_url.rstrip("/") + "/", "locus/")
        announce(f"Inspecting scheme page: {slug}")
        csv_url, fasta_url = discover_scheme_downloads(
            locus_url,
            timeout=timeout,
            retries=retries,
        )
        schemes.append(
            Scheme(
                title=text,
                scheme_url=full_url,
                slug=slug,
                locus_url=locus_url,
                csv_url=csv_url,
                fasta_url=fasta_url,
            )
        )
        seen.add(full_url)

    announce(f"Resolved download links for {len(schemes)} schemes.")
    return schemes


def discover_scheme_downloads(locus_url: str, *, timeout: int, retries: int) -> tuple[str, str]:
    locus_html = _fetch_text(locus_url, timeout=timeout, retries=retries)
    csv_url = ""
    fasta_url = ""

    for href, text in _parse_links(locus_html):
        if text == "Download table as CSV":
            csv_url = urllib.parse.urljoin(locus_url, href)
        elif text == "Download alleles as FASTA":
            fasta_url = urllib.parse.urljoin(locus_url, href)

    if not csv_url:
        separator = "&" if "?" in locus_url else "?"
        csv_url = f"{locus_url}{separator}content-type=csv"
    if not fasta_url:
        fasta_url = urllib.parse.urljoin(locus_url, "../alleles/")

    return csv_url, fasta_url


def build_tasks(schemes: Iterable[Scheme], output_dir: Path) -> list[DownloadTask]:
    tasks: list[DownloadTask] = []
    for scheme in schemes:
        scheme_dir = output_dir / scheme.slug
        tasks.append(
            DownloadTask(
                scheme_slug=scheme.slug,
                scheme_title=scheme.title,
                kind="csv",
                source_url=scheme.csv_url,
                output_path=scheme_dir / "targets.csv",
            )
        )
        tasks.append(
            DownloadTask(
                scheme_slug=scheme.slug,
                scheme_title=scheme.title,
                kind="fasta",
                source_url=scheme.fasta_url,
                output_path=scheme_dir / "alleles.zip",
            )
        )
    return tasks


def download_one(task: DownloadTask, *, timeout: int, retries: int, overwrite: bool) -> DownloadResult:
    task.output_path.parent.mkdir(parents=True, exist_ok=True)
    if not overwrite:
        existing_paths = (
            [path for path in _fasta_candidate_paths(task.output_path.parent) if path.exists()]
            if task.kind == "fasta"
            else ([task.output_path] if task.output_path.exists() else [])
        )
        for existing_path in existing_paths:
            existing_error = _validate_download(task, existing_path)
            if existing_error is None:
                return DownloadResult(
                    scheme_slug=task.scheme_slug,
                    scheme_title=task.scheme_title,
                    kind=task.kind,
                    source_url=task.source_url,
                    output_path=str(existing_path),
                    status="skipped",
                    attempts=0,
                )

    last_exc: Exception | None = None
    for attempt in range(1, retries + 1):
        tmp_path: Path | None = None
        try:
            with urllib.request.urlopen(_request(task.source_url, timeout=timeout), timeout=timeout) as resp:
                with NamedTemporaryFile(
                    mode="wb",
                    prefix=task.output_path.name + ".",
                    suffix=".part",
                    dir=str(task.output_path.parent),
                    delete=False,
                ) as tmp:
                    while True:
                        chunk = resp.read(1024 * 1024)
                        if not chunk:
                            break
                        tmp.write(chunk)
                    tmp_path = Path(tmp.name)
            validation_error = _validate_download(task, tmp_path)
            if validation_error is not None:
                raise RuntimeError(validation_error)
            final_path = task.output_path
            if task.kind == "fasta":
                final_kind = _classify_fasta_path(tmp_path)
                final_path = task.output_path.parent / _fasta_filename_for_kind(final_kind)
            tmp_path.replace(final_path)
            if task.kind == "fasta":
                for stale_path in _fasta_candidate_paths(task.output_path.parent):
                    if stale_path == final_path:
                        continue
                    try:
                        stale_path.unlink(missing_ok=True)
                    except Exception:
                        pass
            return DownloadResult(
                scheme_slug=task.scheme_slug,
                scheme_title=task.scheme_title,
                kind=task.kind,
                source_url=task.source_url,
                output_path=str(final_path),
                status="downloaded",
                attempts=attempt,
            )
        except Exception as exc:
            last_exc = exc
            if tmp_path is not None:
                try:
                    tmp_path.unlink(missing_ok=True)
                except Exception:
                    pass
            if attempt < retries:
                time.sleep(min(2 ** (attempt - 1), 5))

    return DownloadResult(
        scheme_slug=task.scheme_slug,
        scheme_title=task.scheme_title,
        kind=task.kind,
        source_url=task.source_url,
        output_path=str(task.output_path),
        status="failed",
        attempts=retries,
        error=str(last_exc or "unknown error"),
    )


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")


def write_tsv(path: Path, rows: list[dict[str, object]], field_order: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(field_order) + "\n")
        for row in rows:
            values = [str(row.get(field, "")) for field in field_order]
            fh.write("\t".join(values) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download cgMLST.org NCS per-scheme target CSVs and allele FASTAs.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Root output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--index-url",
        default=INDEX_URL,
        help=f"Scheme index URL (default: {INDEX_URL})",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=DEFAULT_WORKERS,
        help=f"Parallel download workers (default: {DEFAULT_WORKERS})",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=DEFAULT_TIMEOUT,
        help=f"Per-request timeout in seconds (default: {DEFAULT_TIMEOUT})",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=DEFAULT_RETRIES,
        help=f"Maximum attempts per fetch/download including validation (default: {DEFAULT_RETRIES})",
    )
    parser.add_argument(
        "--match",
        action="append",
        default=[],
        help="Case-insensitive substring filter; may be passed multiple times.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-download files even if they already exist.",
    )
    parser.add_argument(
        "--list-only",
        action="store_true",
        help="Discover schemes and write manifests without downloading files.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    match_terms = [_normalize_space(term).lower() for term in args.match if _normalize_space(term)]
    output_dir: Path = args.output_dir.expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    announce(f"Output directory: {output_dir}")
    if match_terms:
        announce(f"Filtering schemes with: {', '.join(match_terms)}")

    schemes = discover_schemes(
        args.index_url,
        timeout=args.timeout,
        retries=args.retries,
        match_terms=match_terms,
    )
    if not schemes:
        announce("No matching cgMLST schemes were discovered.")
        return 1

    scheme_payload = [asdict(scheme) for scheme in schemes]
    announce("Writing scheme manifest files.")
    write_json(output_dir / "schemes.json", scheme_payload)
    write_tsv(
        output_dir / "schemes.tsv",
        scheme_payload,
        ["slug", "title", "scheme_url", "locus_url", "csv_url", "fasta_url"],
    )

    announce(f"Discovered {len(schemes)} schemes.")
    if args.list_only:
        announce("List-only mode requested. Skipping downloads.")
        return 0

    tasks = build_tasks(schemes, output_dir)
    results: list[DownloadResult] = []
    failures = 0

    announce(f"Preparing to download {len(tasks)} files across {len(schemes)} schemes.")
    announce(f"Starting downloads with {max(1, args.workers)} worker(s).")

    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.workers)) as executor:
        future_map = {
            executor.submit(
                download_one,
                task,
                timeout=args.timeout,
                retries=args.retries,
                overwrite=args.overwrite,
            ): task
            for task in tasks
        }
        total = len(future_map)
        completed = 0
        for future in concurrent.futures.as_completed(future_map):
            completed += 1
            result = future.result()
            results.append(result)
            if result.status == "failed":
                failures += 1
                print(
                    f"[{completed}/{total}] FAILED {result.scheme_slug} {result.kind} after {result.attempts} attempt(s): {result.error}",
                    file=sys.stderr,
                )
            else:
                attempt_note = ""
                if result.status == "downloaded":
                    attempt_note = f" (attempt {result.attempts}/{args.retries})"
                print(
                    f"[{completed}/{total}] {result.status.upper()} {result.scheme_slug} {result.kind} -> {result.output_path}{attempt_note}",
                    file=sys.stderr,
                )

    result_rows = [asdict(result) for result in results]
    announce("Writing download result manifests.")
    write_json(output_dir / "download_results.json", result_rows)
    write_tsv(
        output_dir / "download_results.tsv",
        result_rows,
        ["scheme_slug", "scheme_title", "kind", "source_url", "output_path", "status", "attempts", "error"],
    )

    if failures:
        announce(f"Finished with {failures} failed download(s).")
        return 1

    announce("Finished without download failures.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
