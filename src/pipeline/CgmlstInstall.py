"""Install cgMLST schemas for chewBBACA: download them from cgMLST.org or import local files.

Every schema is built with :mod:`pipeline.ChewbbacaSchemaBuild` (chewBBACA
PrepExternalSchema in the tools image) into ``<dest>/.staging``, validated (every
representative a forward CDS, then a smoke AlleleCall that must leave a self-score
for every representative), moved into ``<dest>/<schema_id>`` and registered in the
existing ``scheme_species_catalog.json`` (and, for the configured cgMLST folder, in
the config's ``chewbbaca.available_schemas``). Other catalog entries are kept as they
are.

A schema is never overwritten silently: chewBBACA adds every novel allele of an
inferred-mode run to it, so an installed schema is skipped unless ``replace`` is
given, and a replaced one is kept under ``<dest>/.replaced/<schema_id>.<timestamp>``.

cgMLST.org (Ridom) allows its schemes to be downloaded for non-commercial use only and
does not allow redistributing them: :data:`TERMS_TEXT` is shown once before the first
download and the acceptance is recorded in ``<dest>/.cgmlst_org_terms.json``.

Command line (``--dest`` defaults to the configured cgMLST folder)::

    python -m pipeline.CgmlstInstall --list [--sizes]
    python -m pipeline.CgmlstInstall --status
    python -m pipeline.CgmlstInstall --install "Acinetobacter baumannii" Lmonocytogenes
    python -m pipeline.CgmlstInstall --install-all
    python -m pipeline.CgmlstInstall --import ~/Downloads/cgmlst/*.zip [--species S]
    ... [--replace] [--dest DIR] [--dry-run] [--jobs N] [--cpu N] [--accept-terms]
"""
from __future__ import annotations

import argparse
import csv
import getpass
import json
import re
import shutil
import sys
import tempfile
import threading
import time
import urllib.parse
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from html import unescape
from html.parser import HTMLParser
from pathlib import Path
from typing import Callable, Iterable, Optional, Sequence

from pipeline.ChewbbacaSchemaBuild import (
    _atomic_write_json,
    BUILD_INFO_FILE,
    CHEWBBACA_GENE_LIST,
    CHEWBBACA_SCHEME_CATALOG,
    CHEWBBACA_SHORT_DIR,
    FASTA_EXTS,
    SCHEMA_BUILD_REPORT_DIR,
    BuildOptions,
    Cancelled,
    SchemaBuildResult,
    build_schema,
    detect_file_format,
    install_schema_dir,
    is_combined_fasta,
    is_prebuilt_schema,
    read_build_info,
    read_fasta_file,
    register_schema,
    run_watched,
    slugify_species,
)

INDEX_URL = "https://www.cgmlst.org/ncs/"
USER_AGENT = "BACTolomew-cgMLST-installer/1.0"
DEFAULT_TIMEOUT = 120
DEFAULT_RETRIES = 5
DEFAULT_JOBS = 2
STAGING_DIR = ".staging"
TERMS_FILE = ".cgmlst_org_terms.json"
HTML_PREFIXES = (b"<html", b"<!doctype html")
# cgMLST.org schemes whose organisms use another genetic code than 11.
TRANSLATION_TABLES = {"mycoplasma_gallisepticum": 4}
# Genera of cgMLST.org schemes that are not bacteria (Candidozyma auris, a yeast):
# chewBBACA predicts genes with Prodigal, which is for prokaryotes, so they are not installed.
NON_BACTERIAL_GENERA = frozenset({"candidozyma", "candida", "aspergillus", "cryptococcus", "saccharomyces"})
NON_BACTERIAL_REASON = "not a bacterial scheme: chewBBACA predicts genes with Prodigal, which is for prokaryotes"
# Share of an import's loci that must be named like an installed schema's loci to
# take that schema's species when the file name does not tell it.
LOCUS_OVERLAP_MIN = 0.8

TERMS_TEXT = (
    "cgMLST.org (Ridom GmbH) provides the schemes of its Nomenclature Server for "
    "non-commercial use only. The downloaded allele data must not be redistributed. "
    "The full terms of use are on https://www.cgmlst.org. Accept these terms to "
    "download schemes from cgMLST.org."
)

Progress = Optional[Callable[[str], None]]
OnResult = Optional[Callable[["ItemResult", int, int], None]]
ShouldCancel = Optional[Callable[[], bool]]
# on_fraction(schema_id, share): how far one item is (0 to 1, never decreasing), from the
# bytes of its download and the phases of its build.
OnFraction = Optional[Callable[[str, float], None]]

# Share of a cgMLST.org install reached once its files are downloaded; the build does the rest.
DOWNLOAD_SHARE = 0.4
# The build's messages (ChewbbacaSchemaBuild and this module) and the share of the build
# reached when each one is reported; "Staging loci… i/n" moves between 0.05 and 0.3.
_BUILD_PHASES = (
    ("Reading zip archive", 0.02),
    ("Splitting combined FASTA", 0.04),
    ("Copying prebuilt", 0.05),
    ("Staging", 0.05),
    ("Adapting loci with chewBBACA PrepExternalSchema", 0.3),
    ("Validating the adapted schema", 0.6),
    ("Checking the schema with a test allele call", 0.7),
    ("Registering the schema", 0.95),
)
_STAGING_COUNT_RE = re.compile(r"^Staging loci\D*(\d+)/(\d+)")
# "Adapting loci… i/n" (PrepExternalSchema's loci written so far) moves between 0.3 and 0.6.
_ADAPTING_COUNT_RE = re.compile(r"^Adapting loci\D*(\d+)/(\d+)")

_REGISTER_LOCK = threading.Lock()


class TermsNotAccepted(RuntimeError):
    pass


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------
@dataclass
class OrgScheme:
    """One scheme of the cgMLST.org Nomenclature Server index."""

    slug: str
    title: str
    schema_id: str
    species: list[str]
    scheme_url: str
    csv_url: str
    alleles_url: str
    size_bytes: Optional[int] = None

    @property
    def bacterial(self) -> bool:
        return _norm(self.schema_id).split("_")[0] not in NON_BACTERIAL_GENERA


@dataclass
class AvailableScheme:
    schema_id: str
    title: str
    species: list[str]
    mlst_schemes: list[str]
    slug: str  # cgMLST.org slug ("" when the scheme is only in the catalog)
    in_catalog: bool
    installed: bool
    size_bytes: Optional[int] = None
    bacterial: bool = True


@dataclass
class InstalledScheme:
    schema_id: str
    path: str
    loci: int
    alleles: Optional[int]  # at build time (build_info.json), or counted with count_alleles
    built_at: str
    source: str
    builder: str
    in_catalog: bool
    species: list[str] = field(default_factory=list)
    mlst_schemes: list[str] = field(default_factory=list)


@dataclass
class ItemResult:
    source: str
    schema_id: str = ""
    status: str = "failed"  # installed | imported | skipped | failed | planned
    reason: str = ""
    path: str = ""
    backup: str = ""
    loci: int = 0
    alleles: int = 0
    missing_loci: int = 0  # expected loci missing from the source or without a valid CDS allele
    source_format: str = ""
    warning: str = ""
    species: list[str] = field(default_factory=list)
    registered: bool = False
    registration_error: str = ""
    alias_conflict: bool = False  # registration_error: an MLST scheme already maps to another schema


# ---------------------------------------------------------------------------
# Folders, catalog and config
# ---------------------------------------------------------------------------
def default_dest() -> Path:
    """The configured cgMLST folder (``databases.cgmlst_dir``), else ``DB/cgmlst``."""
    from pipeline.ToolEnv import configured_database_dirs, project_root

    raw = configured_database_dirs().get("cgmlst_dir") or ""
    return Path(raw).expanduser() if raw else project_root() / "DB" / "cgmlst"


def _is_configured_dest(dest: Path) -> bool:
    from pipeline.ToolEnv import configured_database_dirs

    raw = configured_database_dirs().get("cgmlst_dir") or ""
    try:
        return bool(raw) and Path(raw).expanduser().resolve() == Path(dest).resolve()
    except OSError:
        return False


def catalog_path_for(dest: Path) -> Path:
    """The catalog runs read for ``dest``: the configured one for the configured folder."""
    if _is_configured_dest(dest):
        from pipeline.ToolEnv import configured_chewbbaca_settings

        configured = str(configured_chewbbaca_settings().get("catalog_path") or "").strip()
        if configured:
            return Path(configured).expanduser()
    return Path(dest) / CHEWBBACA_SCHEME_CATALOG


def load_catalog(catalog_path: Path) -> list[dict]:
    """The catalog's scheme entries ([] when there is no catalog yet)."""
    from pipeline.ChewBBACAWorkflow import _load_scheme_catalog

    return _load_scheme_catalog(Path(catalog_path)) if Path(catalog_path).is_file() else []


def terms_accepted(dest: Path | None = None) -> bool:
    try:
        data = json.loads((Path(dest or default_dest()) / TERMS_FILE).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return isinstance(data, dict) and bool(data.get("accepted_at"))


def accept_terms(dest: Path | None = None) -> Path:
    """Record that the cgMLST.org terms were accepted for downloads into ``dest``."""
    path = Path(dest or default_dest()) / TERMS_FILE
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        user = getpass.getuser()
    except Exception:  # noqa: BLE001
        user = ""
    payload = {
        "accepted_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "user": user,
        "terms": TERMS_TEXT,
    }
    _atomic_write_json(path, payload)
    return path


def _installed_dir(dest: Path, schema_id: str) -> Path:
    return Path(dest) / schema_id


# ---------------------------------------------------------------------------
# cgMLST.org download
# ---------------------------------------------------------------------------
def _normalize_space(text: str) -> str:
    return re.sub(r"\s+", " ", unescape(text or "")).strip()


def _open_url(url: str, *, timeout: int, method: str = "GET"):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "*/*"}, method=method)
    return urllib.request.urlopen(request, timeout=timeout)  # noqa: S310 - https URLs of cgmlst.org


def fetch_text(url: str, *, timeout: int = DEFAULT_TIMEOUT, retries: int = DEFAULT_RETRIES) -> str:
    last_exc: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with _open_url(url, timeout=timeout) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except Exception as exc:  # urllib raises a mix of exception types
            last_exc = exc
            if attempt < retries:
                time.sleep(min(2 ** (attempt - 1), 5))
    raise RuntimeError(f"Failed to fetch {url}: {last_exc}")


class _LinkParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[tuple[str, str]] = []
        self._href: Optional[str] = None
        self._text: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, Optional[str]]]) -> None:
        href = dict(attrs).get("href") if tag.lower() == "a" else None
        if href is not None:
            self._href, self._text = href, []

    def handle_data(self, data: str) -> None:
        if self._href is not None:
            self._text.append(data)

    def handle_endtag(self, tag: str) -> None:
        if tag.lower() == "a" and self._href is not None:
            self.links.append((self._href, _normalize_space("".join(self._text))))
            self._href = None


def _parse_links(html_text: str) -> list[tuple[str, str]]:
    parser = _LinkParser()
    parser.feed(html_text)
    return parser.links


def slugify_title(title: str, fallback_slug: str = "") -> str:
    """The schema folder name of a cgMLST.org title ("Burkholderia mallei (RKI) cgMLST" ->
    "burkholderia_mallei_rki"); the names of the installed schemas come from it."""
    label = re.sub(r"\bcgmlst\b", "", title, flags=re.IGNORECASE)
    for old in ("&", "/", "(", ")"):
        label = label.replace(old, " and " if old == "&" else " ")
    label = re.sub(r"[^A-Za-z0-9]+", "_", label).strip("_").lower()
    return label or fallback_slug.lower()


def species_from_title(title: str) -> list[str]:
    """Species of a cgMLST.org title ("Campylobacter jejuni/coli cgMLST" -> both species)."""
    label = re.sub(r"\bcgmlst\b", "", title, flags=re.IGNORECASE)
    label = _normalize_space(re.sub(r"\([^)]*\)", " ", label))
    words = label.split(" ", 1)
    if len(words) == 2 and "/" in words[1]:
        return [f"{words[0]} {epithet.strip()}" for epithet in words[1].split("/") if epithet.strip()]
    return [label] if label else []


def fetch_org_index(index_url: str = INDEX_URL, *, timeout: int = 60, retries: int = 3) -> list[OrgScheme]:
    """The schemes listed on the cgMLST.org Nomenclature Server index (one request)."""
    links = _parse_links(fetch_text(index_url, timeout=timeout, retries=retries))
    out: list[OrgScheme] = []
    seen: set[str] = set()
    for href, text in links:
        url = urllib.parse.urljoin(index_url, href)
        if not text or "cgMLST" not in text or "/ncs/schema/" not in url or url in seen:
            continue
        seen.add(url)
        slug = urllib.parse.urlparse(url).path.rstrip("/").split("/")[-1] or "scheme"
        scheme_url = url.rstrip("/") + "/"
        out.append(OrgScheme(
            slug=slug,
            title=text,
            schema_id=slugify_title(text, slug),
            species=species_from_title(text),
            scheme_url=scheme_url,
            csv_url=urllib.parse.urljoin(scheme_url, "locus/") + "?content-type=csv",
            alleles_url=urllib.parse.urljoin(scheme_url, "alleles/"),
        ))
    if not out:
        raise RuntimeError(f"no cgMLST schemes found on {index_url}")
    return out


def probe_sizes(schemes: Sequence[OrgScheme], *, timeout: int = 20, workers: int = 8) -> None:
    """Fill ``size_bytes`` from the allele download's Content-Length when the server sends one."""

    def _probe(scheme: OrgScheme) -> None:
        try:
            with _open_url(scheme.alleles_url, timeout=timeout, method="HEAD") as resp:
                length = resp.headers.get("Content-Length")
                scheme.size_bytes = int(length) if length and int(length) > 0 else None
        except Exception:  # noqa: BLE001 - the size is optional
            scheme.size_bytes = None

    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        list(pool.map(_probe, schemes))


def _resolve_download_links(scheme: OrgScheme, *, timeout: int, retries: int) -> tuple[str, str]:
    """The table and allele links of the scheme's locus page (the derived URLs as fallback)."""
    locus_url = urllib.parse.urljoin(scheme.scheme_url, "locus/")
    csv_url, alleles_url = scheme.csv_url, scheme.alleles_url
    try:
        for href, text in _parse_links(fetch_text(locus_url, timeout=timeout, retries=retries)):
            if text == "Download table as CSV":
                csv_url = urllib.parse.urljoin(locus_url, href)
            elif text == "Download alleles as FASTA":
                alleles_url = urllib.parse.urljoin(locus_url, href)
    except RuntimeError:
        pass
    return csv_url, alleles_url


def _looks_like_html(path: Path) -> bool:
    with Path(path).open("rb") as fh:
        prefix = fh.read(256).lstrip().lower()
    return any(prefix.startswith(marker) for marker in HTML_PREFIXES)


def targets_csv_issue(path: Path) -> str:
    """Why ``path`` is not a cgMLST.org locus table ("" when it is one)."""
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return "the locus table is missing or empty"
    if _looks_like_html(path):
        return "received an HTML page instead of the locus table"
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = {f.strip() for f in next(reader, None) or [] if f.strip()}
        if "Gene" not in header or "Alleles" not in header:
            return "the locus table header is not a cgMLST.org targets table"
        if not any(any(f.strip().lower().endswith(".fasta") for f in row) for row in reader):
            return "the locus table names no locus FASTA"
    return ""


def alleles_zip_issue(path: Path, targets_csv: Path | None = None) -> str:
    """Why ``path`` is not a complete cgMLST.org allele archive ("" when it is one)."""
    from pipeline.ChewbbacaSchemaBuild import read_expected_loci

    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return "the allele archive is missing or empty"
    if _looks_like_html(path):
        return "received an HTML page instead of the allele archive"
    if not zipfile.is_zipfile(path):
        return "the allele download is not a zip archive (truncated?)"
    try:
        with zipfile.ZipFile(path) as archive:
            names = {Path(i.filename).name for i in archive.infolist() if not i.is_dir()}
            if not names:
                return "the allele archive is empty"
            bad = archive.testzip()
            if bad is not None:
                return f"the allele archive member {bad} is corrupt"
    except Exception as exc:  # noqa: BLE001
        return f"the allele archive cannot be read: {exc}"
    if targets_csv is not None and Path(targets_csv).is_file():
        missing = set(read_expected_loci(Path(targets_csv))) - names
        if missing:
            return f"the allele archive lacks {len(missing)} loci of the locus table"
    return ""


def _content_length(resp) -> Optional[int]:
    try:
        value = int(resp.getheader("Content-Length") or 0)
    except (AttributeError, TypeError, ValueError):
        return None
    return value or None


def download_file(
    url: str,
    out_path: Path,
    *,
    check: Callable[[Path], str],
    timeout: int = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
    should_cancel: ShouldCancel = None,
    on_bytes: Optional[Callable[[int, Optional[int]], None]] = None,
    expected_size: Optional[int] = None,
) -> Path:
    """Download ``url`` to ``out_path`` (through a .part file), retried until ``check`` passes.

    ``on_bytes(received, total)`` follows the download; ``total`` is the server's
    Content-Length, else ``expected_size``, else None.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    last = ""
    for attempt in range(1, retries + 1):
        part = out_path.with_name(out_path.name + ".part")
        try:
            with _open_url(url, timeout=timeout) as resp, part.open("wb") as fh:
                total = _content_length(resp) or expected_size
                received = 0
                while True:
                    if should_cancel is not None and should_cancel():
                        raise Cancelled("cancelled")
                    chunk = resp.read(1024 * 1024)
                    if not chunk:
                        break
                    fh.write(chunk)
                    received += len(chunk)
                    if on_bytes is not None:
                        on_bytes(received, total)
            issue = check(part)
            if not issue:
                part.replace(out_path)
                return out_path
            last = issue
        except Cancelled:
            part.unlink(missing_ok=True)
            raise
        except Exception as exc:  # noqa: BLE001 - retried
            last = str(exc)
        part.unlink(missing_ok=True)
        if attempt < retries:
            time.sleep(min(2 ** (attempt - 1), 5))
    raise RuntimeError(f"download of {url} failed after {retries} attempt(s): {last}")


# ---------------------------------------------------------------------------
# Validation: smoke AlleleCall (P03 self-scores)
# ---------------------------------------------------------------------------
# Stops in all six frames (TTAA is its own reverse complement): no CDS spans two loci.
_SPACER = "TTAA" * 30


def write_synthetic_genome(schema_dir: Path, out_path: Path) -> Path:
    """A genome made of the first representative of every locus, one spacer apart.

    AlleleCall on it computes short/self_scores for every representative and calls
    each locus, without needing a real genome of the species.
    """
    parts: list[str] = []
    for short in sorted((Path(schema_dir) / CHEWBBACA_SHORT_DIR).glob("*_short.fasta")):
        records = read_fasta_file(short)
        if records:
            parts.append(records[0][1])
    if not parts:
        raise RuntimeError(f"{schema_dir} has no representative alleles")
    sequence = _SPACER.join(parts)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    width = 80
    out_path.write_text(
        ">synthetic_representatives\n"
        + "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width)) + "\n",
        encoding="utf-8",
    )
    return out_path


def smoke_allele_call(
    schema: Path,
    out_dir: Path,
    *,
    genomes: Sequence[Path] = (),
    cpu: int = 1,
    base_cmd: list[str] | None = None,
    should_cancel: ShouldCancel = None,
) -> dict:
    """AlleleCall without inferred alleles (the schema gains none); checks the self-scores.
    ``should_cancel`` stops it (Cancelled).

    Without ``genomes`` a synthetic genome of the schema's representatives is used.
    Raises RuntimeError when the call fails or a representative has no self-score.
    """
    from pipeline.ChewBBACAWorkflow import load_self_scores, self_score_issue
    from pipeline.ToolEnv import gui_tool_cmd

    schema, out_dir = Path(schema).resolve(), Path(out_dir).resolve()
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)
    synthetic = not genomes
    inputs = [write_synthetic_genome(schema, out_dir / "synthetic.fasta")] if synthetic \
        else [Path(g).resolve() for g in genomes]
    listing = out_dir / "genomes.txt"
    listing.write_text("".join(f"{g}\n" for g in inputs), encoding="utf-8")
    cmd = list(base_cmd) if base_cmd else gui_tool_cmd(
        "chewbbaca", "chewBBACA.py", mount_paths=[schema, out_dir, *inputs])
    cmd += ["AlleleCall", "-i", str(listing), "-g", str(schema), "-o", str(out_dir / "out"),
            "--cpu", str(max(1, int(cpu))), "--no-inferred"]
    if synthetic:
        cmd += ["--pm", "meta"]
    t0 = time.time()
    res = run_watched(cmd, should_cancel=should_cancel)
    (out_dir / "command.stdout.log").write_text(res.stdout or "", encoding="utf-8")
    (out_dir / "command.stderr.log").write_text(res.stderr or "", encoding="utf-8")
    missing_lines = (res.stdout or "").count("Could not get the self-score")
    scores = load_self_scores(schema)
    summary = {
        "returncode": res.returncode,
        "seconds": round(time.time() - t0, 1),
        "genomes": "synthetic" if synthetic else [str(g) for g in inputs],
        "missing_self_score_lines": missing_lines,
        "self_scores": len(scores) if scores is not None else -1,
    }
    if res.returncode != 0:
        tail = "\n".join((res.stderr or res.stdout or "").strip().splitlines()[-10:])
        raise RuntimeError(f"smoke AlleleCall failed (exit {res.returncode}): {tail}")
    issue = self_score_issue(schema) if scores is not None else "AlleleCall wrote no short/self_scores"
    if issue or missing_lines:
        raise RuntimeError(f"self-scores: {issue or ''} ({missing_lines} missing-score lines)")
    return summary


def _build_and_validate(
    source: Path,
    built: Path,
    *,
    targets_csv: Path | None,
    single_fasta_mode: str,
    options: BuildOptions,
    smoke_genomes: Sequence[Path],
    smoke_dir: Path,
    progress: Progress,
) -> SchemaBuildResult:
    """Build ``source`` into ``built`` (orientation checked by build_schema), then the smoke call."""
    result = build_schema(source, built, single_fasta_mode=single_fasta_mode, targets_csv=targets_csv,
                          options=options, progress=progress)
    if progress:
        progress("Checking the schema with a test allele call…")
    smoke_allele_call(built, smoke_dir, genomes=smoke_genomes, cpu=options.cpu or 1, base_cmd=options.base_cmd,
                      should_cancel=options.should_cancel)
    return result


# ---------------------------------------------------------------------------
# Listing
# ---------------------------------------------------------------------------
def available(
    dest: Path | None = None,
    *,
    index: Sequence[OrgScheme] | None = None,
    fetch: bool = True,
    sizes: bool = False,
    progress: Progress = None,
) -> list[AvailableScheme]:
    """The schemes known from cgMLST.org and the catalog, with their install state.

    ``index`` is a cgMLST.org index already fetched; with ``fetch`` it is fetched here,
    and when cgMLST.org cannot be reached only the catalog entries are listed.
    """
    dest = Path(dest or default_dest())
    entries = load_catalog(catalog_path_for(dest))
    if index is None and fetch:
        try:
            index = fetch_org_index()
        except Exception as exc:  # noqa: BLE001 - listed from the catalog only
            if progress:
                progress(f"cgMLST.org could not be reached ({exc}); listing the catalog only.")
            index = []
    index = list(index or [])
    if sizes and index:
        probe_sizes(index)
    by_id = {e["schema_id"]: e for e in entries}
    rows: list[AvailableScheme] = []
    for scheme in index:
        entry = by_id.pop(scheme.schema_id, None)
        rows.append(AvailableScheme(
            schema_id=scheme.schema_id,
            title=scheme.title,
            species=list(entry["species"]) if entry and entry["species"] else list(scheme.species),
            mlst_schemes=list(entry["mlst_schemes"]) if entry else [],
            slug=scheme.slug,
            in_catalog=entry is not None,
            installed=_installed_dir(dest, scheme.schema_id).is_dir(),
            size_bytes=scheme.size_bytes,
            bacterial=scheme.bacterial,
        ))
    for entry in by_id.values():
        rows.append(AvailableScheme(
            schema_id=entry["schema_id"],
            title=(entry["species"] or [entry["schema_id"]])[0],
            species=list(entry["species"]),
            mlst_schemes=list(entry["mlst_schemes"]),
            slug="",
            in_catalog=True,
            installed=_installed_dir(dest, entry["schema_dir"]).is_dir(),
        ))
    rows.sort(key=lambda r: r.title.lower())
    return rows


def _count_alleles(schema_dir: Path) -> int:
    total = 0
    for fasta in Path(schema_dir).glob("*.fasta"):
        with fasta.open("rb") as fh:
            total += sum(1 for line in fh if line.startswith(b">"))
    return total


def installed(dest: Path | None = None, *, count_alleles: bool = False) -> list[InstalledScheme]:
    """The schemas installed in ``dest`` (folders with a .genes_list)."""
    dest = Path(dest or default_dest())
    entries = {e["schema_dir"]: e for e in load_catalog(catalog_path_for(dest))}
    out: list[InstalledScheme] = []
    if not dest.is_dir():
        return out
    for child in sorted(dest.iterdir(), key=lambda p: p.name.lower()):
        if child.name.startswith(".") or not (child / CHEWBBACA_GENE_LIST).is_file():
            continue
        info = read_build_info(child)
        entry = entries.get(child.name)
        alleles = _count_alleles(child) if count_alleles else info.get("alleles")
        out.append(InstalledScheme(
            schema_id=child.name,
            path=str(child),
            loci=sum(1 for p in child.glob("*.fasta")),
            alleles=int(alleles) if alleles is not None else None,
            built_at=str(info.get("built_at") or ""),
            source=str(info.get("source") or ""),
            builder=str(info.get("builder") or ""),
            in_catalog=entry is not None,
            species=list(entry["species"]) if entry else [],
            mlst_schemes=list(entry["mlst_schemes"]) if entry else [],
        ))
    return out


# ---------------------------------------------------------------------------
# Species detection
# ---------------------------------------------------------------------------
def _norm(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(text).lower()).strip("_")


# Names of a source that say nothing about its species.
_NOISE_TOKENS = frozenset({"cgmlst", "alleles", "allele", "schema", "scheme", "ncs", "download", "fasta",
                           "zip", "chewbbaca", "ridom", "loci", "truncated"})


@dataclass
class _Target:
    schema_id: str
    species: list[str]
    mlst_schemes: list[str]
    # Names that identify the schema, matched as whole tokens of a source name, with a
    # rank: 3 for its id, cgMLST.org slug and title and their "Abaumannii" style short
    # forms, 2 for the short form of a complex's first species ("cjejuni").
    keys: dict[str, int] = field(default_factory=dict)
    # Names of its species: matched only when they are the whole source name (less its
    # noise tokens), so "Bmallei_fli" is not taken for a "Burkholderia mallei" schema.
    species_keys: set[str] = field(default_factory=set)

    def add(self, key: str, rank: int) -> None:
        if key:
            self.keys[key] = max(rank, self.keys.get(key, 0))


def _species_keys(species: str) -> set[str]:
    """"Acinetobacter baumannii" -> {"acinetobacter_baumannii", "abaumannii"}."""
    words = _norm(species).split("_")
    keys = {"_".join(words)} if words and words[0] else set()
    if len(words) >= 2 and words[0] and words[1] not in ("spp", "sp"):
        keys.add(words[0][0] + words[1])
    return keys


def _add_schema_id_keys(target: _Target, schema_id: str) -> None:
    """"burkholderia_mallei_rki" -> itself, "bmallei_rki"; "brucella_spp" -> "brucella";
    "campylobacter_jejuni_coli" -> "cjejuni_coli", "cjejuni_complex" and (rank 2) "cjejuni"."""
    words = _norm(schema_id).split("_")
    target.add("_".join(words), 3)
    if len(words) < 2:
        return
    if words[1] in ("spp", "sp"):
        target.add(words[0], 3)
        return
    short = words[0][0] + words[1]
    target.add("_".join([short, *words[2:]]), 3)
    if len(words) > 2 and not any(len(w) <= 3 for w in words[2:]):  # a complex, not a variant (rki)
        target.add(f"{short}_complex", 3)
        target.add(short, 2)


def _targets(entries: Iterable[dict], index: Iterable[OrgScheme]) -> list[_Target]:
    """Every schema a source can be detected as, with the names that identify it."""
    targets: dict[str, _Target] = {}
    for e in entries:
        t = targets[e["schema_id"]] = _Target(e["schema_id"], list(e["species"]), list(e["mlst_schemes"]))
        _add_schema_id_keys(t, e["schema_id"])
        for sp in e["species"]:
            t.species_keys |= _species_keys(sp)
    for s in index:
        t = targets.get(s.schema_id)
        if t is None:
            t = targets[s.schema_id] = _Target(s.schema_id, list(s.species), [])
            _add_schema_id_keys(t, s.schema_id)
            for sp in s.species:
                t.species_keys |= _species_keys(sp)
        t.add(_norm(s.slug), 3)
        t.add(_norm(s.title), 3)
    for t in targets.values():
        # A genus alone ("brucella" of brucella_spp) names the schema only when it is the
        # whole name: "Brucella_abortus" is another species, not the genus-wide schema.
        genera = {_norm(sp).split("_")[0] for sp in t.species} | {_norm(t.schema_id).split("_")[0]}
        for key in [k for k in t.keys if k in genera]:
            del t.keys[key]
            t.species_keys.add(key)
    return list(targets.values())


def _match_name(name: str, targets: Sequence[_Target]) -> tuple[Optional[_Target], str]:
    """The target ``name`` identifies (None, reason when several do; None, "" when none)."""
    norm = _norm(name)
    padded = f"_{norm}_"
    bare = "_".join(t for t in norm.split("_") if t and t not in _NOISE_TOKENS and not t.isdigit())
    best: dict[str, tuple[int, int]] = {}
    for t in targets:
        scores = [(rank, len(k)) for k, rank in t.keys.items() if f"_{k}_" in padded]
        scores += [(1, len(k)) for k in t.species_keys if k == bare]
        if scores:
            best[t.schema_id] = max(scores)
    if not best:
        return None, ""
    top = max(best.values())
    winners = sorted(sid for sid, score in best.items() if score == top)
    if len(winners) > 1:
        return None, f"the name '{name}' matches several schemas ({', '.join(winners)})"
    return next(t for t in targets if t.schema_id == winners[0]), ""


def _source_locus_names(item: "ImportItem") -> set[str]:
    try:
        if item.kind in ("zip", "ncs_folder"):
            with zipfile.ZipFile(item.source) as archive:
                return {Path(i.filename).stem for i in archive.infolist() if not i.is_dir()}
        if item.source.is_dir():
            return {p.stem for p in item.source.iterdir() if p.is_file() and p.suffix.lower() in FASTA_EXTS}
    except Exception:  # noqa: BLE001 - detection is best effort
        pass
    return set()


def _match_content(item: "ImportItem", dest: Path, targets: Sequence[_Target]) -> Optional[_Target]:
    """The installed schema whose locus names cover the source's (unique), else None."""
    loci = _source_locus_names(item)
    if not loci:
        return None
    hits = []
    for t in targets:
        schema = _installed_dir(dest, t.schema_id)
        if not schema.is_dir():
            continue
        names = {p.stem for p in schema.glob("*.fasta")}
        if names and len(loci & names) / len(loci) >= LOCUS_OVERLAP_MIN:
            hits.append(t)
    return hits[0] if len(hits) == 1 else None


# Locus names of each cgMLST.org scheme (from its locus table), fetched once per process.
_ORG_LOCI: dict[str, frozenset[str]] = {}
_ORG_LOCI_LOCK = threading.Lock()


def _org_locus_names(scheme: OrgScheme, *, timeout: int = 60) -> frozenset[str]:
    """The locus names of a cgMLST.org scheme (empty when its table cannot be fetched)."""
    with _ORG_LOCI_LOCK:
        cached = _ORG_LOCI.get(scheme.csv_url)
    if cached is not None:
        return cached
    try:
        text = fetch_text(scheme.csv_url, timeout=timeout, retries=2)
        rows = csv.reader(text.splitlines(), delimiter="\t")
        next(rows, None)
        names = frozenset(
            Path(f.strip()).stem for row in rows for f in row if f.strip().lower().endswith(".fasta")
        )
    except Exception:  # noqa: BLE001 - detection is best effort
        return frozenset()
    with _ORG_LOCI_LOCK:
        _ORG_LOCI[scheme.csv_url] = names
    return names


def _match_org_content(item: "ImportItem", index: Sequence[OrgScheme], targets: Sequence[_Target],
                       progress: Progress = None) -> Optional[_Target]:
    """The cgMLST.org scheme whose locus table names the source's loci (unique), else None.

    For downloads whose file name says nothing: a browser saves every cgMLST.org allele
    download as "null" (the server sends ``filename=null``).
    """
    loci = _source_locus_names(item)
    if not loci or not index:
        return None
    if progress:
        progress(f"Identifying {item.source.name} from the cgMLST.org locus lists…")
    with ThreadPoolExecutor(max_workers=8) as pool:
        tables = list(pool.map(_org_locus_names, index))
    hits = [scheme for scheme, names in zip(index, tables)
            if names and len(loci & names) / len(loci) >= LOCUS_OVERLAP_MIN]
    if len(hits) != 1:
        return None
    scheme = hits[0]
    found = next((t for t in targets if t.schema_id == scheme.schema_id), None)
    return found or _Target(scheme.schema_id, list(scheme.species), [])


# ---------------------------------------------------------------------------
# Import sources
# ---------------------------------------------------------------------------
@dataclass
class ImportItem:
    source: Path  # what build_schema reads
    kind: str  # zip | ncs_folder | prebuilt | loci_folder | fasta
    names: list[str]  # names the species is detected from, most specific first
    targets_csv: Optional[Path] = None
    label: str = ""
    title: str = ""  # cgMLST.org title from the download folder's schemes.json


def _manifest_titles(folder: Path) -> dict[str, str]:
    """slug -> title from a cgMLST.org download folder's schemes.json, when there is one."""
    try:
        data = json.loads((folder / "schemes.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    if not isinstance(data, list):
        return {}
    return {str(r.get("slug")): str(r.get("title") or "") for r in data if isinstance(r, dict) and r.get("slug")}


def _is_archive(path: Path) -> bool:
    return detect_file_format(path) in ("zip_archive", "zip_stream")


def _item_for_file(path: Path, *, allow_fasta: bool = True) -> Optional[ImportItem]:
    fmt = detect_file_format(path)
    titles = _manifest_titles(path.parent.parent)
    if fmt in ("zip_archive", "zip_stream"):
        targets = path.parent / "targets.csv"
        generic = path.stem.lower() in ("alleles", "alleles.truncated")
        title = titles.get(path.parent.name, "") if generic else ""
        names = [title, path.parent.name] if generic else [path.stem, path.parent.name]
        return ImportItem(path, "zip", [n for n in names if n],
                          targets if generic and targets.is_file() else None, label=str(path), title=title)
    if allow_fasta and fmt == "combined_fasta" and path.suffix.lower() in FASTA_EXTS:
        return ImportItem(path, "fasta", [path.stem, path.parent.name], label=str(path))
    return None


def discover_import_items(paths: Iterable[Path], *, max_depth: int = 3) -> tuple[list[ImportItem], list[ItemResult]]:
    """The schema sources among ``paths``: cgMLST.org zips, folders of per-locus FASTAs,
    chewBBACA schemas, single FASTAs, or folders holding many of these."""
    items: list[ImportItem] = []
    errors: list[ItemResult] = []
    seen: set[Path] = set()

    def add(item: ImportItem) -> None:
        key = item.source.resolve()
        if key not in seen:
            seen.add(key)
            items.append(item)

    def visit(path: Path, depth: int, top: bool) -> bool:
        if path.name.startswith(".") and not top:
            return False
        if path.is_file():
            # In a folder of downloads only archives are sources; a FASTA there (a genome,
            # a note) is not a schema. A FASTA given by itself is one.
            item = _item_for_file(path, allow_fasta=top)
            if item is not None:
                add(item)
                return True
            if top:
                errors.append(ItemResult(str(path), reason="not a zip archive, FASTA file or schema folder"))
            return False
        if not path.is_dir():
            if top:
                errors.append(ItemResult(str(path), reason="not found"))
            return False
        title = _manifest_titles(path.parent).get(path.name, "")
        names = [n for n in (title, path.name) if n]
        if is_prebuilt_schema(path):
            add(ImportItem(path, "prebuilt", names, label=str(path), title=title))
            return True
        if (path / "alleles.zip").is_file():
            targets = path / "targets.csv"
            add(ImportItem(path / "alleles.zip", "ncs_folder", names,
                           targets if targets.is_file() else None, label=str(path), title=title))
            return True
        children = [p for p in path.iterdir() if not p.name.startswith(".")]
        has_fastas = any(p.is_file() and p.suffix.lower() in FASTA_EXTS for p in children)
        # A folder of per-locus FASTAs, unless it also holds archives or schema folders:
        # then it is a folder of downloads and its stray FASTAs are not loci.
        holds_sources = has_fastas and any(
            (p.is_file() and p.suffix.lower() not in FASTA_EXTS and _is_archive(p))
            or (p.is_dir() and (is_prebuilt_schema(p) or (p / "alleles.zip").is_file()))
            for p in children
        )
        if has_fastas and not holds_sources:
            add(ImportItem(path, "loci_folder", names, label=str(path), title=title))
            return True
        if (path / "targets.csv").is_file():
            errors.append(ItemResult(str(path), reason="has a cgMLST.org locus table but no complete allele archive (alleles.zip)"))
            return True
        found = False
        if depth < max_depth:
            for child in sorted(path.iterdir(), key=lambda p: p.name.lower()):
                found = visit(child, depth + 1, False) or found
        if top and not found:
            errors.append(ItemResult(str(path), reason="no schema sources found in this folder"))
        return found

    for raw in paths:
        visit(Path(raw).expanduser(), 0, True)
    return items, errors


# ---------------------------------------------------------------------------
# Install / import
# ---------------------------------------------------------------------------
def _register(dest: Path, catalog: Path, target: _Target, extra_mlst: Sequence[str],
              update_config: bool, result: ItemResult) -> None:
    """Add or update this schema's catalog entry; the other entries are kept as they are."""
    with _REGISTER_LOCK:
        current = {e["schema_id"]: e for e in load_catalog(catalog)}.get(target.schema_id)
        species = list(current["species"]) if current else []
        species += [s for s in target.species if s not in species]
        mlst = list(current["mlst_schemes"]) if current else []
        mlst += [m for m in list(target.mlst_schemes) + list(extra_mlst) if m and m not in mlst]
        try:
            register_schema(dest, target.schema_id, species, mlst, catalog_path=catalog,
                            update_config=update_config)
            result.registered = True
        except ValueError as exc:
            result.registration_error, result.alias_conflict = str(exc), True
        except Exception as exc:  # noqa: BLE001 - the schema is installed; reported
            result.registration_error = str(exc) or type(exc).__name__
        if not mlst:
            result.warning = "; ".join(filter(None, [
                result.warning, "no MLST scheme maps to it, so runs do not auto-select it",
            ]))


def _run_items(
    dest: Path,
    jobs: int,
    work: list[tuple[ItemResult, Callable[[Path], object]]],
    on_result: OnResult,
    should_cancel: ShouldCancel,
) -> list[ItemResult]:
    """Run each item's body (given its own staging folder) with at most ``jobs`` at a time."""
    done = 0
    lock = threading.Lock()

    def run(entry: tuple[ItemResult, Callable[[Path], object]]) -> None:
        nonlocal done
        res, body = entry
        if should_cancel is not None and should_cancel():
            res.status, res.reason = "failed", "cancelled"
        else:
            staging = _new_staging(dest, res.schema_id)
            try:
                body(staging)
            except Cancelled:
                res.status, res.reason = "failed", "cancelled"
            except Exception as exc:  # noqa: BLE001 - reported per item
                res.status, res.reason = "failed", str(exc) or type(exc).__name__
            finally:
                _remove_staging(staging)
        with lock:
            done += 1
            if on_result is not None:
                on_result(res, done, len(work))

    with ThreadPoolExecutor(max_workers=max(1, int(jobs))) as pool:
        list(pool.map(run, work))
    return [res for res, _body in work]


def _report_settled(results: Sequence[ItemResult], on_fraction: OnFraction) -> None:
    """Items settled without a build (skipped, planned or failed before it) are done."""
    if on_fraction is not None:
        for res in results:
            on_fraction(res.schema_id or res.source, 1.0)


def _cpu_per_job(cpu: int | None, jobs: int) -> int:
    if cpu:
        return max(1, int(cpu))
    from pipeline.ToolEnv import effective_core_budget

    return max(1, effective_core_budget() // max(1, jobs))


def _build_share(message: str) -> Optional[float]:
    """The share of a build reached when it reports ``message`` (None: not a phase)."""
    for pattern, start, end in ((_STAGING_COUNT_RE, 0.05, 0.3), (_ADAPTING_COUNT_RE, 0.3, 0.6)):
        count = pattern.match(message)
        if count and int(count.group(2)):
            return start + (end - start) * min(1.0, int(count.group(1)) / int(count.group(2)))
    return next((share for prefix, share in _BUILD_PHASES if message.startswith(prefix)), None)


class _ItemProgress:
    """Reports one item's step and how far it is (``on_fraction``); the build fills the
    shares after ``build_from`` (what the download took). Raises Cancelled when the user
    cancelled (unless ``check`` is False: once the schema is installed it is always
    registered)."""

    def __init__(self, progress: Progress, name: str, should_cancel: ShouldCancel,
                 on_fraction: OnFraction = None, *, build_from: float = 0.0) -> None:
        self._progress = progress
        self._name = name
        self._should_cancel = should_cancel
        self._on_fraction = on_fraction
        self._build_from = build_from
        self._reached = -1.0

    def __call__(self, message: str, *, check: bool = True) -> None:
        if check and self._should_cancel is not None and self._should_cancel():
            raise Cancelled("cancelled")
        if self._progress is not None:
            self._progress(f"[{self._name}] {message}")
        share = _build_share(message)
        if share is not None:
            self.fraction(self._build_from + (1.0 - self._build_from) * share)

    def fraction(self, value: float) -> None:
        value = max(0.0, min(1.0, float(value)))
        if value <= self._reached:
            return
        self._reached = value
        if self._on_fraction is not None:
            self._on_fraction(self._name, value)


def _build_one(
    *,
    source: Path,
    target: _Target,
    dest: Path,
    catalog: Path,
    replace: bool,
    result: ItemResult,
    targets_csv: Optional[Path],
    single_fasta_mode: str,
    options: BuildOptions,
    smoke_genomes: Sequence[Path],
    extra_mlst: Sequence[str],
    update_config: bool,
    emit: _ItemProgress,
    staging: Path,
    done_status: str,
    origin: str = "",
) -> ItemResult:
    """Build into ``staging``, validate, install (keeping a replaced schema) and register.

    ``origin`` replaces the build record's source (the staging copy of a download).
    """
    built = staging / target.schema_id
    build = _build_and_validate(source, built, targets_csv=targets_csv, single_fasta_mode=single_fasta_mode,
                                options=options, smoke_genomes=smoke_genomes, smoke_dir=staging / "smoke",
                                progress=emit)
    if origin:
        info = read_build_info(built)
        info["source"] = origin
        (built / SCHEMA_BUILD_REPORT_DIR / BUILD_INFO_FILE).write_text(json.dumps(info, indent=2) + "\n",
                                                                       encoding="utf-8")
    final = _installed_dir(dest, target.schema_id)
    with _REGISTER_LOCK:
        if final.exists() and not replace:
            # Another install (a second run needing the same scheme) finished it first.
            result.status, result.reason, result.path = "skipped", "installed meanwhile by another install", str(final)
            emit.fraction(1.0)
            return result
        backup = install_schema_dir(built, final, replace=replace)
    result.backup = str(backup or "")
    result.path = str(final)
    result.loci = build.written_loci
    result.alleles = build.alleles
    result.missing_loci = build.missing_loci
    result.source_format = build.source_format
    result.warning = build.warning
    result.status = done_status
    emit("Registering the schema…", check=False)
    _register(dest, catalog, target, extra_mlst, update_config, result)
    emit.fraction(1.0)
    return result


_STAGING_LOCK = threading.Lock()


def _new_staging(dest: Path, schema_id: str) -> Path:
    root = Path(dest) / STAGING_DIR
    with _STAGING_LOCK:
        root.mkdir(parents=True, exist_ok=True)
        return Path(tempfile.mkdtemp(prefix=f"{schema_id}.", dir=str(root)))


def _remove_staging(staging: Path) -> None:
    """Remove an item's staging folder, and ``<dest>/.staging`` once it is empty."""
    shutil.rmtree(staging, ignore_errors=True)
    with _STAGING_LOCK:
        try:
            staging.parent.rmdir()
        except OSError:  # another item is still building there
            pass


def _resolve_org(name: str, index: Sequence[OrgScheme], entries: Sequence[dict]) -> tuple[Optional[OrgScheme], str]:
    targets = _targets(entries, index)
    t, why = _match_name(name, targets)
    if t is None:
        return None, why or f"'{name}' is not a cgMLST.org scheme"
    scheme = next((s for s in index if s.schema_id == t.schema_id), None)
    if scheme is None:
        return None, f"'{t.schema_id}' is in the catalog but not on cgMLST.org"
    return scheme, ""


def install(
    species: Sequence[str],
    dest: Path | None = None,
    *,
    replace: bool = False,
    dry_run: bool = False,
    install_all: bool = False,
    jobs: int = DEFAULT_JOBS,
    cpu: int | None = None,
    smoke_genomes: Sequence[Path] = (),
    accept: bool = False,
    update_config: bool | None = None,
    index: Sequence[OrgScheme] | None = None,
    timeout: int = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
    progress: Progress = None,
    on_result: OnResult = None,
    should_cancel: ShouldCancel = None,
    on_fraction: OnFraction = None,
) -> list[ItemResult]:
    """Download ``species`` (names, schema ids or cgMLST.org slugs) from cgMLST.org,
    build, validate and register them; ``install_all`` takes every scheme not installed.

    Raises TermsNotAccepted (before any download) unless the cgMLST.org terms were
    accepted for ``dest`` or ``accept`` records the acceptance now.
    """
    dest = Path(dest or default_dest())
    catalog = catalog_path_for(dest)
    update_config = _is_configured_dest(dest) if update_config is None else update_config
    entries = load_catalog(catalog)
    index = list(index) if index is not None else fetch_org_index()
    jobs = max(1, int(jobs))
    options_cpu = _cpu_per_job(cpu, jobs)

    planned: list[tuple[OrgScheme, ItemResult]] = []
    results: list[ItemResult] = []
    names = ([s.slug for s in index if s.bacterial and not _installed_dir(dest, s.schema_id).exists()]
             if install_all else list(species))
    seen: set[str] = set()
    for name in names:
        scheme, why = _resolve_org(name, index, entries)
        if scheme is None:
            results.append(ItemResult(name, reason=why))
            continue
        if scheme.schema_id in seen:
            continue
        seen.add(scheme.schema_id)
        res = ItemResult(name, schema_id=scheme.schema_id, species=list(scheme.species),
                         path=str(_installed_dir(dest, scheme.schema_id)))
        if not scheme.bacterial:
            res.reason = NON_BACTERIAL_REASON
            results.append(res)
        elif _installed_dir(dest, scheme.schema_id).exists() and not replace:
            res.status, res.reason = "skipped", "already installed (use replace to rebuild it)"
            results.append(res)
        elif dry_run:
            res.status = "planned"
            res.reason = f"download {scheme.alleles_url}" + (" and replace the installed schema" if
                                                             _installed_dir(dest, scheme.schema_id).exists() else "")
            results.append(res)
        else:
            planned.append((scheme, res))
    if planned and not terms_accepted(dest):
        if not accept:
            raise TermsNotAccepted(TERMS_TEXT)
        accept_terms(dest)
    _report_settled(results, on_fraction)

    def body(scheme: OrgScheme, res: ItemResult, staging: Path) -> None:
        emit = _ItemProgress(progress, scheme.schema_id, should_cancel, on_fraction, build_from=DOWNLOAD_SHARE)
        entry = {e["schema_id"]: e for e in entries}.get(scheme.schema_id)
        target = _Target(scheme.schema_id,
                         list(entry["species"]) if entry and entry["species"] else list(scheme.species),
                         list(entry["mlst_schemes"]) if entry else [])
        emit("Downloading from cgMLST.org…")
        emit.fraction(0.0)
        csv_url, alleles_url = _resolve_download_links(scheme, timeout=timeout, retries=retries)
        download = staging / "download"
        targets_csv = download_file(csv_url, download / "targets.csv", check=targets_csv_issue,
                                    timeout=timeout, retries=retries, should_cancel=should_cancel)
        emit.fraction(0.02)

        def on_bytes(received: int, total: Optional[int]) -> None:
            if total:
                emit.fraction(0.02 + (DOWNLOAD_SHARE - 0.02) * min(1.0, received / total))

        alleles = download_file(alleles_url, download / "alleles.zip",
                                check=lambda p: alleles_zip_issue(p, targets_csv),
                                timeout=timeout, retries=retries, should_cancel=should_cancel,
                                on_bytes=on_bytes, expected_size=scheme.size_bytes)
        options = BuildOptions(cpu=options_cpu, translation_table=TRANSLATION_TABLES.get(scheme.schema_id, 11),
                               should_cancel=should_cancel)
        _build_one(source=alleles, target=target, dest=dest, catalog=catalog, replace=replace, result=res,
                   targets_csv=targets_csv, single_fasta_mode="loci_set", options=options,
                   smoke_genomes=smoke_genomes, extra_mlst=(), update_config=update_config, emit=emit,
                   staging=staging, done_status="installed", origin=scheme.scheme_url)

    work = [(res, lambda staging, s=scheme, r=res: body(s, r, staging)) for scheme, res in planned]
    return results + _run_items(dest, jobs, work, on_result, should_cancel)


def _target_for_species(species: str, targets: Sequence[_Target]) -> _Target:
    """The schema a species typed by the user names: the schema of that id, else the one its
    name identifies, else the only one listing that species, else a new schema."""
    species = species.strip()
    schema_id = slugify_species(species)
    found = next((t for t in targets if t.schema_id == schema_id), None)
    if found is None:
        found, _why = _match_name(species, targets)
    if found is None:
        listing = [t for t in targets if species.lower() in (s.lower() for s in t.species)]
        found = listing[0] if len(listing) == 1 else None
    return found or _Target(schema_id, [species], [])


def schema_id_for_species(species: str, dest: Path | None = None) -> str:
    """The schema id an import with this species (typed by the user) is installed as."""
    dest = Path(dest or default_dest())
    return _target_for_species(species, _targets(load_catalog(catalog_path_for(dest)), [])).schema_id


def import_paths(
    paths: Sequence[Path],
    species: str | None = None,
    dest: Path | None = None,
    *,
    mlst_schemes: Sequence[str] = (),
    replace: bool = False,
    dry_run: bool = False,
    single_fasta_mode: str = "auto",  # auto | loci_set | single_locus
    jobs: int = 1,
    cpu: int | None = None,
    smoke_genomes: Sequence[Path] = (),
    update_config: bool | None = None,
    index: Sequence[OrgScheme] | None = None,
    progress: Progress = None,
    on_result: OnResult = None,
    should_cancel: ShouldCancel = None,
    on_fraction: OnFraction = None,
) -> list[ItemResult]:
    """Import schema sources the user downloaded: files and/or folders, many at once.

    Each source (a cgMLST.org zip or download folder, a folder of per-locus FASTAs,
    a chewBBACA schema, a single FASTA, or a folder holding many of these) gets its
    species from ``species``, else from its file or folder name, else from its locus
    names when they match an installed schema or a cgMLST.org scheme's locus table
    (``index``; fetched only when needed, never when an empty list is given). Reported
    per source: imported, skipped (already installed) or failed with the reason;
    ``dry_run`` only plans ("planned"). ``mlst_schemes`` needs a single source (an MLST
    scheme maps to one schema): ValueError otherwise.
    """
    dest = Path(dest or default_dest())
    catalog = catalog_path_for(dest)
    update_config = _is_configured_dest(dest) if update_config is None else update_config
    entries = load_catalog(catalog)
    org_index = list(index) if index is not None else None
    targets = _targets(entries, org_index or [])
    items, results = discover_import_items(paths)
    jobs = max(1, int(jobs))
    options_cpu = _cpu_per_job(cpu, jobs)

    forced = _target_for_species(species, targets) if species and species.strip() else None

    planned: list[tuple[ImportItem, _Target, ItemResult]] = []
    claimed: dict[str, str] = {}
    for item in items:
        if should_cancel is not None and should_cancel():
            raise Cancelled("cancelled")
        res = ItemResult(item.label or str(item.source))
        target, why = forced, ""
        if target is None:
            for name in item.names:
                target, why = _match_name(name, targets)
                if target is not None or why:
                    break
        if target is None and not why:
            target = _match_content(item, dest, targets)
        if target is None and not why and not item.title:
            if org_index is None:
                try:
                    org_index = fetch_org_index()
                except Exception:  # noqa: BLE001 - offline: the name decides
                    org_index = []
            target = _match_org_content(item, org_index, targets, progress)
        if target is None and not why and item.title:
            # A cgMLST.org scheme without a catalog entry yet: named after its title.
            target = _Target(slugify_title(item.title), species_from_title(item.title), [])
        if target is None:
            res.reason = (why + "; " if why else "") + "species not detected from the name; give the species"
            results.append(res)
            continue
        res.schema_id, res.species = target.schema_id, list(target.species)
        final = _installed_dir(dest, target.schema_id)
        res.path = str(final)
        if target.schema_id in claimed:
            res.reason = f"another source of this import ({claimed[target.schema_id]}) is installed as {target.schema_id}"
            results.append(res)
            continue
        claimed[target.schema_id] = res.source
        try:
            inside = item.source.resolve().is_relative_to(final.resolve())
        except OSError:
            inside = False
        if inside:
            res.status, res.reason = "skipped", "this is the installed schema"
        elif final.exists() and not replace:
            res.status, res.reason = "skipped", "already installed (use replace to rebuild it)"
        elif dry_run:
            res.status = "planned"
            res.reason = f"import as {target.schema_id}" + (" and replace the installed schema" if final.exists() else "")
        else:
            planned.append((item, target, res))
            continue
        results.append(res)

    building = len(planned) + sum(r.status == "planned" for r in results)
    if mlst_schemes and building > 1:
        raise ValueError(f"MLST schemes can be given for one schema at a time; this import has {building}. "
                         "Import the sources one by one, or leave the MLST schemes empty.")

    _report_settled(results, on_fraction)

    def mode_for(item: ImportItem) -> str:
        if item.kind != "fasta" or single_fasta_mode in ("loci_set", "single_locus"):
            return "loci_set" if item.kind != "fasta" else single_fasta_mode
        try:
            return "loci_set" if is_combined_fasta(item.source) else "single_locus"
        except Exception:  # noqa: BLE001
            return "loci_set"

    def body(item: ImportItem, target: _Target, res: ItemResult, staging: Path) -> None:
        emit = _ItemProgress(progress, target.schema_id, should_cancel, on_fraction)
        emit(f"Importing {item.source.name}…")
        emit.fraction(0.0)
        options = BuildOptions(cpu=options_cpu, translation_table=TRANSLATION_TABLES.get(target.schema_id, 11),
                               should_cancel=should_cancel)
        _build_one(source=item.source, target=target, dest=dest, catalog=catalog, replace=replace, result=res,
                   targets_csv=item.targets_csv, single_fasta_mode=mode_for(item), options=options,
                   smoke_genomes=smoke_genomes, extra_mlst=mlst_schemes, update_config=update_config,
                   emit=emit, staging=staging, done_status="imported")

    work = [(res, lambda staging, i=item, t=target, r=res: body(i, t, r, staging)) for item, target, res in planned]
    return results + _run_items(dest, jobs, work, on_result, should_cancel)


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------
def _human_size(n: Optional[int]) -> str:
    if not n:
        return ""
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return ""


def _print_results(results: Sequence[ItemResult]) -> None:
    for r in results:
        detail = r.reason or (f"{r.loci} loci, {r.alleles} alleles -> {r.path}" if r.path else "")
        extra = []
        if r.backup:
            extra.append(f"previous schema kept in {r.backup}")
        if r.warning:
            extra.append(r.warning)
        if r.registration_error:
            extra.append(f"not registered: {r.registration_error}")
        tail = f" ({'; '.join(extra)})" if extra else ""
        print(f"{r.status:<9} {r.schema_id or '-':<45} {r.source}\n          {detail}{tail}")


def _ask_terms(dest: Path, accept_flag: bool) -> bool:
    if terms_accepted(dest):
        return True
    print(TERMS_TEXT, file=sys.stderr)
    if accept_flag:
        accept_terms(dest)
        return True
    if sys.stdin.isatty():
        answer = input("Accept these terms? [y/N] ").strip().lower()
        if answer in ("y", "yes"):
            accept_terms(dest)
            return True
    print("Not accepted: pass --accept-terms to accept them.", file=sys.stderr)
    return False


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m pipeline.CgmlstInstall",
        description="Install chewBBACA cgMLST schemas from cgMLST.org or from local files.",
    )
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--list", action="store_true", help="Schemes on cgMLST.org and in the catalog.")
    action.add_argument("--status", action="store_true", help="Installed schemas.")
    action.add_argument("--install", nargs="+", metavar="SPECIES",
                        help="Download, build and register these schemes (species, schema id or slug).")
    action.add_argument("--install-all", action="store_true", help="Every cgMLST.org scheme not installed.")
    action.add_argument("--import", dest="import_paths", nargs="+", type=Path, metavar="PATH",
                        help="Zips, folders or schemas downloaded by hand (folders may hold many).")
    parser.add_argument("--species", help="Species of the imported source(s) when the names do not tell it.")
    parser.add_argument("--mlst", default="", help="Comma-separated MLST scheme(s) mapping to an imported schema.")
    parser.add_argument("--replace", action="store_true",
                        help="Rebuild installed schemas; the old one is kept under .replaced/.")
    parser.add_argument("--dest", type=Path, help="cgMLST folder (default: the configured one).")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be done.")
    parser.add_argument("--jobs", type=int, default=DEFAULT_JOBS, help="Schemas built in parallel.")
    parser.add_argument("--cpu", type=int, help="CPU threads per build (default: budget / jobs).")
    parser.add_argument("--smoke-genome", type=Path, action="append", default=[],
                        help="Genome for the test allele call (default: one made of the representatives).")
    parser.add_argument("--sizes", action="store_true", help="With --list: ask cgMLST.org for download sizes.")
    parser.add_argument("--count-alleles", action="store_true", help="With --status: count the current alleles.")
    parser.add_argument("--no-config", action="store_true",
                        help="Do not add the schemas to bactolomew_config.json (catalog only).")
    parser.add_argument("--accept-terms", action="store_true",
                        help="Accept the cgMLST.org terms (shown once) without asking; alone, only records it.")
    parser.add_argument("--report", type=Path, help="Write the per-item results as JSON.")
    args = parser.parse_args(argv)

    dest = (args.dest or default_dest()).expanduser()
    update_config = False if args.no_config else None
    progress = lambda msg: print(msg, file=sys.stderr, flush=True)  # noqa: E731

    if not (args.list or args.status or args.install or args.install_all or args.import_paths):
        if args.accept_terms:
            return 0 if _ask_terms(dest, True) else 1
        parser.error("one of --list, --status, --install, --install-all, --import or --accept-terms is required")
    if args.list:
        for row in available(dest, sizes=args.sizes, progress=progress):
            state = "installed" if row.installed else ("" if row.slug else "catalog only")
            print(f"{row.schema_id:<50} {row.title:<60} {_human_size(row.size_bytes):>9} {state}")
        return 0
    if args.status:
        for row in installed(dest, count_alleles=args.count_alleles):
            print(f"{row.schema_id:<50} {row.loci:>6} loci {row.alleles if row.alleles is not None else '?':>9} "
                  f"alleles  built {row.built_at or '?'}  {'' if row.in_catalog else '(not in catalog) '}"
                  f"{row.source}")
        return 0

    try:
        if args.import_paths:
            mlst = [a.strip() for a in args.mlst.replace(";", ",").split(",") if a.strip()]
            results = import_paths(args.import_paths, args.species, dest, mlst_schemes=mlst, replace=args.replace,
                                   dry_run=args.dry_run, jobs=args.jobs, cpu=args.cpu,
                                   smoke_genomes=args.smoke_genome, update_config=update_config,
                                   progress=progress)
        else:
            index = fetch_org_index()

            def run(accept: bool) -> list[ItemResult]:
                return install(args.install or [], dest, replace=args.replace, dry_run=args.dry_run,
                               install_all=args.install_all, jobs=args.jobs, cpu=args.cpu,
                               smoke_genomes=args.smoke_genome, accept=accept, update_config=update_config,
                               index=index, progress=progress)

            try:
                results = run(args.accept_terms)
            except TermsNotAccepted:
                # Asked only when something is to be downloaded.
                if not _ask_terms(dest, False):
                    return 1
                results = run(True)
    except (TermsNotAccepted, ValueError) as exc:
        print(exc, file=sys.stderr)
        return 1
    _print_results(results)
    if args.report:
        args.report.write_text(json.dumps([asdict(r) for r in results], indent=1) + "\n", encoding="utf-8")
    return 1 if any(r.status == "failed" for r in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
