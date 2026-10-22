#!/usr/bin/env python
"""BACTolomew install demo: the app's real install windows on a simulated network and tools.

Run from the repository root:

    PYTHONPATH=src .venv/bin/python src/scripts/demo_install.py [--speed 4] [--keep]

A launcher window walks through the install flows in order. Every step opens the app's
own windows (First-Time Setup, Settings with "Download from cgMLST.org..." and
"Add schema...", the start window and the Run progress window), and their real code runs:
threads, progress windows, validation, the .partial folders, install records, the
cgMLST catalog and the config. Only the lowest layer is replaced:

* the network: ``urllib.request.urlopen`` answers for www.cgmlst.org and NCBI's
  AMRFinderPlus folder from generated data (any other address is refused);
* the tools: ``DbDownload._stream`` (bakta_db, checkm2, amrfinder_update, amrfinder and
  the ``--version`` checks in the tools image), ``ContainerRuntime._stream`` (docker pull
  and docker tag of the tools image) and the ``subprocess.run`` of ChewbbacaSchemaBuild and
  CgmlstInstall (chewBBACA PrepExternalSchema, AlleleCall, the schema metadata reader)
  print realistic output and write small placeholder files (large database files are
  sparse, so their sizes are realistic and they take no disk space);
* the Docker checks: the daemon counts as running; the tools image is missing until
  First-Time Setup pulls it from a simulated registry (fastGEAR's image counts as present).

Any other attempt to start docker, singularity, conda or a bioinformatics tool raises.
Everything is written to a fresh temporary folder (printed at start, removed on exit
unless --keep): its own bactolomew_config.json, DB/cgmlst (with a copy of the cgMLST
catalog and a few installed schemes; First-Time Setup makes the rest of DB/), runs/, a
folder of fake browser downloads, and the QSettings. The user's config, DB, runs, caches and ~/.config are never touched.

Durations: each simulated operation takes about the square root of its real duration
in seconds (a 15 minute download takes 30 s), divided by --speed.
"""
from __future__ import annotations

import argparse
import io
import json
import logging
import os
import pickle
import random
import re
import shlex
import shutil
import signal
import socket
import subprocess
import sys
import tempfile
import threading
import time
import traceback
import urllib.error
import urllib.parse
import urllib.request
import zipfile
import zlib
from dataclasses import dataclass
from datetime import datetime, timedelta
from email.message import Message
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable, Optional

_SRC = Path(__file__).resolve().parents[1]
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from PyQt6 import QtCore, QtGui, QtWidgets  # noqa: E402

from annotation_explorer_polars.theme import COLORS  # noqa: E402

REPO_ROOT = _SRC.parent
DEMO_TAG = "Demo: nothing is downloaded"
CHEWBBACA_VERSION = "3.5.4"
AMRF_RELEASE = "2025-07-16.1"
AMRF_SOFTWARE = "4.0.23"
GB = 1 << 30
MB = 1 << 20

# --------------------------------------------------------------------------- #
# Timing
# --------------------------------------------------------------------------- #
SPEED = 1.0
_SHUTDOWN = threading.Event()


def demo_seconds(real_seconds: float) -> float:
    """How long the demo spends on an operation that takes ``real_seconds`` for real."""
    real = max(0.0, float(real_seconds))
    scaled = real if real <= 1 else real ** 0.5
    return scaled / max(0.05, SPEED)


def pause(seconds: float) -> None:
    """Sleep ``seconds`` of demo time (returns early when the demo shuts down)."""
    _SHUTDOWN.wait(max(0.0, seconds))


def human_size(n: Optional[float]) -> str:
    if n is None:
        return "?"
    value = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024 or unit == "TB":
            return f"{value:.0f} {unit}" if unit in ("B", "KB") else f"{value:.1f} {unit}"
        value /= 1024
    return "?"


def _now_iso(delta: float = 0.0) -> str:
    return (datetime.now().astimezone() + timedelta(seconds=delta)).isoformat(timespec="seconds")


# --------------------------------------------------------------------------- #
# The cgMLST.org schemes the simulated server lists
# --------------------------------------------------------------------------- #
# slug | title | loci | allele download size in bytes ("" = the server sends no size) | first locus.
# Titles, slugs, loci counts and sizes are those of the cgMLST.org Nomenclature Server
# downloads of 2026-05; the allele data the demo serves is generated.
_ORG_TABLE = """\
Abaumannii|Acinetobacter baumannii cgMLST|2390|75100848|ACICU_RS00005
Banthracis|Bacillus anthracis cgMLST|3803|15052264|GBAA_RS00015
Bpertussis|Bordetella pertussis cgMLST|2983|2694425|BP0001
Bmelitensis|Brucella melitensis cgMLST|2704|4525104|BME_RS00005
Brucella|Brucella spp. cgMLST|1764|2854552|BME_RS00005
Bmallei_fli|Burkholderia mallei (FLI) cgMLST|2838|6111038|BMA0002
Bmallei_rki|Burkholderia mallei (RKI) cgMLST|3328|7480182|BMA_RS00010
Bpseudomallei|Burkholderia pseudomallei cgMLST|4221|19615389|BPSL0001
Cjejuni_complex|Campylobacter jejuni/coli cgMLST|637|42895329|Cj0001
Cauris|Candidozyma auris cgMLST|1248|9830112|B9J08_000001
Cfreundii|Citrobacter freundii cgMLST|3250|33429250|G4551_RS00005
Cfreundii_complex|Citrobacter freundii/portucalensis/braakii/europaeus cgMLST|2307|18902319|G4551_RS00005
Cdifficile|Clostridioides difficile cgMLST|2147|21084001|CD630_00010
Cperfringens|Clostridium perfringens cgMLST|1431|7510382|CPF_RS00015
Cdiphtheriae|Corynebacterium diphtheriae cgMLST|1553|5861223|DIP_RS11520
Cpseudotuberculosis|Corynebacterium pseudotuberculosis cgMLST|1577|1191266|CpMEX29_RS00005
Csakazakii_complex|Cronobacter sakazakii/malonaticus cgMLST|2831|19731268|ESA_RS00005
Ehormaechei|Enterobacter hormaechei cgMLST|2178|22376515|LZ581_RS01210
Efaecalis|Enterococcus faecalis cgMLST|1972|15220160|OG1RF_RS00015
Efaecium|Enterococcus faecium cgMLST|1423|17644397|EFAU004_00001
Ecoli|Escherichia coli cgMLST|2513||EC_RS00005
Ftularensis|Francisella tularensis cgMLST|1147|1125687|FTL_0001
Koxytoca_complex|Klebsiella oxytoca/grimontii/michiganensis/pasteurii cgMLST|2536|24811678|EGY21_RS00010
Kpneumoniae_complex|Klebsiella pneumoniae/variicola/quasipneumoniae cgMLST|2358||KP1_RS00005
Lpneumophila|Legionella pneumophila cgMLST|1521|8467643|lpg0004
Lmonocytogenes|Listeria monocytogenes cgMLST|1701|42447079|lmo0001
Mmorganii|Morganella morganii cgMLST|2462|5570111|FX120_RS00010
Mtuberculosis_complex|Mycobacterium tuberculosis/bovis/africanum/canettii cgMLST|2891|35401267|Rv0014
Mabscessus|Mycobacteroides abscessus cgMLST|2904|9934902|MAB_0002
Mgallisepticum|Mycoplasma gallisepticum cgMLST|425|578818|MGA_RS00005
Plarvae|Paenibacillus larvae cgMLST|2419|1625270|ERICI_RS00005
Pmirabilis|Proteus mirabilis cgMLST|2574|4736707|PMI_RS00005
Pstuartii|Providencia stuartii cgMLST|3079|2067828|AL507_00010
Paeruginosa|Pseudomonas aeruginosa cgMLST|3867|45581240|PA0001
Senterica|Salmonella enterica cgMLST|3002|157237945|SEN_RS00005
Smarcescens|Serratia marcescens cgMLST|2692|57508605|SMDB11_RS00005
Sargenteus|Staphylococcus argenteus cgMLST|2168|1536834|IAPLMMLD_00001
Saureus|Staphylococcus aureus cgMLST|1861||SACOL0001
Scapitis|Staphylococcus capitis cgMLST|1492|2499275|SCAP_RS00005
Spyogenes|Streptococcus pyogenes cgMLST|1095|4960429|ETT58_RS00005
Yenterocolitica|Yersinia enterocolitica cgMLST|2582||YE_RS00005
"""

# Schemes installed in the demo's DB/cgmlst from the start ("already installed").
PREINSTALLED = ("Abaumannii", "Lmonocytogenes", "Paeruginosa")
# Species a simulated run is made of (the first one whose schema is not installed).
RUN_SPECIES_ORDER = ("enterococcus_faecium", "streptococcus_pyogenes", "clostridium_perfringens",
                     "proteus_mirabilis", "providencia_stuartii", "serratia_marcescens",
                     "enterococcus_faecalis")

# A species of step 4 that has no scheme on cgMLST.org and none in the catalog: (species,
# MLST scheme, first locus). Its schema files (a folder of per-locus FASTAs) are imported
# with Add schema.
NO_SCHEME_SPECIES = ("Vibrio parahaemolyticus", "vparahaemolyticus", "VP_RS00005")
NO_SCHEME_LOCI = 60

_SENSE_CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
                 if a + b + c not in ("TAA", "TAG", "TGA")]


@dataclass(frozen=True)
class OrgSchemeData:
    slug: str
    title: str
    loci: int
    size: Optional[int]
    first_locus: str

    @property
    def schema_id(self) -> str:
        from pipeline.CgmlstInstall import slugify_title

        return slugify_title(self.title, self.slug)

    @property
    def scheme_url(self) -> str:
        return f"https://www.cgmlst.org/ncs/schema/{self.slug}/"

    def locus_names(self) -> list[str]:
        return _locus_names(self.first_locus, self.loci)

    def alleles(self, locus: str, *, max_alleles: int = 4) -> list[str]:
        return _alleles(f"{self.slug}/{locus}", max_alleles)


def _org_schemes() -> list[OrgSchemeData]:
    out = []
    for line in _ORG_TABLE.strip().splitlines():
        slug, title, loci, size, first = line.split("|")
        out.append(OrgSchemeData(slug, title, int(loci), int(size) if size else None, first))
    return out


ORG = {s.slug: s for s in _org_schemes()}


def org_by_schema_id(schema_id: str) -> Optional[OrgSchemeData]:
    return next((s for s in ORG.values() if s.schema_id == schema_id), None)


def _locus_names(first: str, count: int) -> list[str]:
    found = re.match(r"^(.*?)(\d+)$", first)
    prefix, digits = (found.group(1), found.group(2)) if found else (first + "_", "1")
    start, width = int(digits), len(digits)
    step = 5 if start % 5 == 0 else 1
    return [f"{prefix}{start + i * step:0{width}d}" for i in range(count)]


def _alleles(key: str, max_alleles: int) -> list[str]:
    """Deterministic alleles of one locus: complete forward CDSs, now and then one with an
    internal stop codon (PrepExternalSchema drops those)."""
    rng = random.Random(zlib.crc32(key.encode()))
    codons = [rng.choice(_SENSE_CODONS) for _ in range(rng.randint(50, 140))]
    base = ["ATG", *codons, "TAA"]
    out = ["".join(base)]
    for _ in range(rng.randint(0, max(0, max_alleles - 1))):
        variant = list(base)
        for _m in range(rng.randint(1, 3)):
            variant[rng.randrange(1, len(variant) - 1)] = rng.choice(_SENSE_CODONS)
        if rng.random() < 0.01:
            variant[len(variant) // 2] = "TGA"
        out.append("".join(variant))
    return out


class _StreamSink(io.RawIOBase):
    """A write-only stream that cannot seek, so zipfile writes data descriptors after each
    member, as cgMLST.org's streamed downloads have them."""

    def __init__(self) -> None:
        super().__init__()
        self.data = bytearray()

    def writable(self) -> bool:
        return True

    def write(self, b) -> int:  # noqa: ANN001
        self.data += bytes(b)
        return len(b)


@lru_cache(maxsize=64)
def org_alleles_zip(slug: str) -> bytes:
    """The allele archive cgMLST.org serves for a scheme: one FASTA per locus, alleles >1, >2..."""
    scheme = ORG[slug]
    sink = _StreamSink()
    with zipfile.ZipFile(sink, "w", zipfile.ZIP_DEFLATED) as archive:
        for locus in scheme.locus_names():
            text = "".join(f">{i}\n{seq}\n" for i, seq in enumerate(scheme.alleles(locus), 1))
            archive.writestr(f"{locus}.fasta", text)
    return bytes(sink.data)


@lru_cache(maxsize=64)
def org_targets_csv(slug: str) -> str:
    scheme = ORG[slug]
    rows = ["Locus\tGene\tStart\tLength\tProduct\tAlleles\t"]
    pos = 76
    for locus in scheme.locus_names():
        length = len(scheme.alleles(locus, max_alleles=1)[0])
        rows.append(f"\t{locus}\t\t{pos}\t{length}\thypothetical protein\t{locus}.fasta\t")
        pos += length + 120
    return "\n".join(rows) + "\n"


# --------------------------------------------------------------------------- #
# Console: everything printed (by the demo, the app's code and the fake tools)
# --------------------------------------------------------------------------- #
class _ConsoleTee(io.TextIOBase):
    """sys.stdout / sys.stderr replacement: writes through, and hands complete lines to the
    demo console (from any thread)."""

    def __init__(self, stream, sink: Callable[[str], None]) -> None:
        super().__init__()
        self._stream = stream
        self._sink = sink
        self._buf = ""
        self._lock = threading.Lock()

    def write(self, text: str) -> int:
        try:
            self._stream.write(text)
            self._stream.flush()
        except Exception:  # noqa: BLE001 - a closed terminal must not stop the demo
            pass
        with self._lock:
            self._buf += text
            *lines, self._buf = self._buf.split("\n")
        for line in lines:
            self._sink(line)
        return len(text)

    def flush(self) -> None:
        try:
            self._stream.flush()
        except Exception:  # noqa: BLE001
            pass

    def isatty(self) -> bool:
        return False


def say(text: str) -> None:
    """A demo line in the console (and the terminal)."""
    print(text, flush=True)


# --------------------------------------------------------------------------- #
# Simulated network
# --------------------------------------------------------------------------- #
class FakeResponse(io.BytesIO):
    """An HTTP response whose body arrives at the simulated rate (read() returns at most
    64 KB per call, as a socket would; callers loop until it returns b"")."""

    def __init__(self, body: bytes, url: str, *, seconds: float = 0.0,
                 headers: Optional[dict[str, str]] = None, status: int = 200) -> None:
        super().__init__(body)
        self.url = url
        self.status = status
        self.code = status
        self.headers = Message()
        for key, value in (headers or {}).items():
            self.headers[key] = value
        self._total = max(1, len(body))
        self._seconds = seconds
        if not body and seconds:
            pause(seconds)

    def read(self, size: Optional[int] = -1) -> bytes:  # type: ignore[override]
        if _SHUTDOWN.is_set():
            raise urllib.error.URLError("demo closed")
        if size is None or size < 0:
            data = super().read()
        else:
            data = super().read(min(size, 64 * 1024))
        if data and self._seconds:
            pause(self._seconds * len(data) / self._total)
        return data

    def getcode(self) -> int:
        return self.status

    def getheader(self, name: str, default: Optional[str] = None) -> Optional[str]:
        return self.headers.get(name, default)

    def info(self) -> Message:
        return self.headers


_CGMLST_RATE = 6 * MB  # bytes per second from cgMLST.org
_NCBI_RATE = 12 * MB


def _cgmlst_index_html() -> str:
    links = "\n".join(f'<li><a href="/ncs/schema/{s.slug}/">{s.title}</a></li>' for s in ORG.values())
    return (f"<!doctype html><html><head><title>cgMLST.org Nomenclature Server</title></head><body>"
            f"<h1>Nomenclature Server</h1><ul>\n{links}\n</ul></body></html>")


def _cgmlst_locus_html(scheme: OrgSchemeData) -> str:
    return (f"<!doctype html><html><body><h1>{scheme.title}</h1>"
            f'<a href="?content-type=csv">Download table as CSV</a> '
            f'<a href="../alleles/">Download alleles as FASTA</a></body></html>')


def _http_error(url: str, code: int, reason: str) -> urllib.error.HTTPError:
    return urllib.error.HTTPError(url, code, reason, Message(), None)


def _cgmlst_respond(url: str, parsed: urllib.parse.ParseResult, method: str) -> FakeResponse:
    path = parsed.path
    if path.rstrip("/") == "/ncs":
        return FakeResponse(_cgmlst_index_html().encode(), url, seconds=demo_seconds(1.5))
    found = re.match(r"^/ncs/schema/([^/]+)/(locus/|alleles/)?$", path)
    scheme = ORG.get(found.group(1)) if found else None
    if scheme is None:
        raise _http_error(url, 404, "Not Found")
    part = found.group(2) or ""
    if part == "locus/" and "content-type=csv" in parsed.query:
        return FakeResponse(org_targets_csv(scheme.slug).encode(), url, seconds=demo_seconds(1.0),
                            headers={"Content-Type": "text/csv"})
    if part == "locus/":
        return FakeResponse(_cgmlst_locus_html(scheme).encode(), url, seconds=demo_seconds(0.6))
    if part == "alleles/":
        headers = {"Content-Type": "application/zip",
                   # cgMLST.org names every allele download "null" (a browser saves null.zip).
                   "Content-Disposition": "attachment; filename=null"}
        if scheme.size:
            headers["Content-Length"] = str(scheme.size)
        if method == "HEAD":
            return FakeResponse(b"", url, seconds=demo_seconds(0.3), headers=headers)
        size = scheme.size or 40 * MB
        say(f"[simulated network] GET {url}  ({human_size(size)}, about {size / _CGMLST_RATE:.0f} s "
            f"at {human_size(_CGMLST_RATE)}/s)")
        body = org_alleles_zip(scheme.slug)
        if scheme.size:  # the small stand-in's own length, so a progress bar follows its bytes
            headers["Content-Length"] = str(len(body))
        return FakeResponse(body, url, seconds=demo_seconds(size / _CGMLST_RATE), headers=headers)
    return FakeResponse(f"<html><body><h1>{scheme.title}</h1></body></html>".encode(), url,
                        seconds=demo_seconds(0.6))


_AMRF_FILES = {
    "version.txt": AMRF_RELEASE + "\n",
    "database_format_version.txt": "4.0\n",
    "AMRProt.fa": ">0|WP_000027057.1|1|1|blaTEM-156|blaTEM|hydrolase|2|BETA-LACTAM|BETA-LACTAM|class A\nMSIQHFRVALIPFFAAFCLPVFA\n",
    "AMR_CDS.fa": ">WP_000027057.1:1-861|blaTEM-156\nATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTTAA\n",
    "ReferenceGeneCatalog.txt": "allele\tgene_family\twhitelisted_taxa\tproduct_name\nblaTEM-156\tblaTEM\t\tclass A beta-lactamase TEM-156\n",
    "taxgroup.tsv": "#taxgroup\tgpipe_taxgroup\tnumber_of_nucl_ref_genes\nAcinetobacter_baumannii\tAcinetobacter\t0\n",
    "changes.txt": f"{AMRF_RELEASE}: demo release notes\n",
}


def _ncbi_respond(url: str, parsed: urllib.parse.ParseResult, method: str) -> FakeResponse:
    found = re.match(r"^/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/([^/]+)/([^/]+)/(.*)$",
                     parsed.path)
    if not found:
        raise _http_error(url, 404, "Not Found")
    name = found.group(3)
    if not name:
        links = "\n".join(f'<a href="{n}">{n}</a>' for n in sorted(_AMRF_FILES))
        return FakeResponse(f"<html><body><pre>{links}</pre></body></html>".encode(), url,
                            seconds=demo_seconds(1.0))
    body = _AMRF_FILES.get(name)
    if body is None:
        raise _http_error(url, 404, "Not Found")
    say(f"[simulated network] GET {url}")
    return FakeResponse(body.encode(), url, seconds=demo_seconds(2.0))


def fake_urlopen(url: Any, data: Any = None, timeout: Any = None, *args: Any, **kwargs: Any) -> FakeResponse:
    """urllib.request.urlopen of the demo: cgMLST.org and NCBI's AMRFinderPlus folder are
    simulated, every other address is refused."""
    if isinstance(url, urllib.request.Request):
        full, method = url.get_full_url(), url.get_method()
    else:
        full, method = str(url), "GET"
    parsed = urllib.parse.urlparse(full)
    host = parsed.netloc.lower()
    if host.endswith("cgmlst.org"):
        return _cgmlst_respond(full, parsed, method)
    if host == "ftp.ncbi.nlm.nih.gov":
        return _ncbi_respond(full, parsed, method)
    raise urllib.error.URLError(f"{DEMO_TAG}: network access to {host or full} is blocked")


# --------------------------------------------------------------------------- #
# Simulated tools: the commands DbDownload runs in the tools image
# --------------------------------------------------------------------------- #
def _split_tool(cmd: list[str]) -> tuple[str, list[str]]:
    """(executable, its arguments) of a ToolEnv command (docker/singularity + conda run)."""
    if "conda" in cmd:
        tail = cmd[cmd.index("conda"):]
        for i, token in enumerate(tail):
            if token in ("-n", "-p") and i + 2 < len(tail):
                return tail[i + 2], tail[i + 3:]
    return (cmd[0] if cmd else ""), cmd[1:]


def _arg(args: list[str], flag: str, default: str = "") -> str:
    return args[args.index(flag) + 1] if flag in args and args.index(flag) + 1 < len(args) else default


def _sparse(path: Path, size: int, head: bytes = b"") -> None:
    """A file of ``size`` bytes that takes no disk space (a hole after ``head``)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(head)
        fh.truncate(max(size, len(head)))


def _tqdm_line(done: float, total: float, elapsed: float, rate: float) -> str:
    pct = int(100 * done / total) if total else 100
    filled = int(pct / 10)
    bar = "█" * filled + " " * (10 - filled)
    eta = (total - done) / rate if rate else 0

    def clock(sec: float) -> str:
        sec = int(sec)
        return f"{sec // 60:02d}:{sec % 60:02d}" if sec < 3600 else f"{sec // 3600}:{sec % 3600 // 60:02d}:{sec % 60:02d}"

    return f"{pct:3d}%|{bar}| {human_size(done)}/{human_size(total)} [{clock(elapsed)}<{clock(eta)}, {human_size(rate)}/s]"


def _simulate_transfer(emit: Callable[[str], None], total: int, rate: float, *, label: str = "") -> None:
    """Progress lines of a download of ``total`` bytes at ``rate`` bytes/s, in demo time."""
    real = total / rate
    wall = demo_seconds(real)
    ticks = max(4, min(60, int(wall / 0.5)))
    for i in range(1, ticks + 1):
        pause(wall / ticks)
        if _SHUTDOWN.is_set():
            return
        done = total * i / ticks
        emit((label + " " if label else "") + _tqdm_line(done, total, real * i / ticks, rate))


def _amrf_update_lines(emit: Callable[[str], None], db_root: Path, *, indent: str = "") -> Path:
    """amrfinder_update --database db_root: the newest release of the image's format."""
    emit(f"{indent}Running: amrfinder_update --database {db_root} --force_update")
    emit(f"{indent}Software directory: '/opt/conda/envs/amrfinder/bin/'")
    emit(f"{indent}Software version: {AMRF_SOFTWARE}")
    pause(demo_seconds(3))
    emit(f"{indent}Looking up the published databases at "
         "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/4.0/")
    release_dir = db_root / AMRF_RELEASE
    emit(f"{indent}Downloading AMRFinder database version {AMRF_RELEASE} into '{release_dir}/'")
    for name, body in _AMRF_FILES.items():
        pause(demo_seconds(2.5))
        emit(f"{indent}  {name}")
        release_dir.mkdir(parents=True, exist_ok=True)
        (release_dir / name).write_text(body, encoding="utf-8")
    _sparse(release_dir / "AMR.LIB", 186 * MB, b"HMMER3/f\n")
    _sparse(release_dir / "AMR_DNA-Escherichia", 31 * MB, b">")
    pause(demo_seconds(20))
    emit(f"{indent}Indexing")
    pause(demo_seconds(12))
    latest = db_root / "latest"
    if latest.is_symlink() or latest.exists():
        latest.unlink()
    latest.symlink_to(AMRF_RELEASE)
    emit(f"{indent}Database directory: '{release_dir}'")
    return release_dir


# The full v5.1 database of the image's bakta 1.10.4 (sizes as in DbDownload).
_BAKTA = {"dir": "db", "tarball": 39_401_553_410, "unpacked": 76_302_905_767,
          "extract": 780, "checksum": 210, "record": "10522951"}


def _bakta_db_download(args: list[str], emit: Callable[[str], None]) -> int:
    out = Path(_arg(args, "--output"))
    kind = _arg(args, "--type", "full")
    spec = _BAKTA
    if kind != "full" or not out.is_dir():
        emit(f"ERROR: invalid output path or DB type ({out}, {kind})")
        return 1
    emit("Bakta software version: 1.10.4")
    emit("Required database schema version: 5")
    emit("")
    emit("Fetch DB versions...")
    pause(demo_seconds(3))
    emit("\ttested DB versions: 3")
    emit("\tcompatible DB versions: 1 (v5.1)")
    name = f"{spec['dir']}.tar.gz"
    url = f"https://zenodo.org/record/{spec['record']}/files/{name}"
    emit(f"Download database: v5.1, type={kind}, 2024-01-19, DOI: 10.5281/zenodo.{spec['record']}, URL: {url}")
    _simulate_transfer(emit, spec["tarball"], 32 * MB)
    emit(f"\ttarball: {out / name}")
    emit("Check database checksum...")
    pause(demo_seconds(spec["checksum"]))
    emit("\tdatabase checksum successfully validated")
    emit(f"Extract DB tarball: file={out / name}, output={out}")
    pause(demo_seconds(spec["extract"]))
    db = out / spec["dir"]
    unpacked = spec["unpacked"]
    shares = {"bakta.db": 0.448, "psc.dmnd": 0.535, "pfam.h3f": 0.002, "pfam.h3m": 0.0032,
              "pfam.h3p": 0.0038, "sorf.dmnd": 0.00003, "antifam.h3m": 0.00008,
              "expert-protein-sequences.dmnd": 0.0002, "ncRNA-genes.i1m": 0.0007}
    scale = 1 / sum(shares.values())
    for fname, share in shares.items():
        head = b"SQLite format 3\x00" if fname.endswith(".db") else b"\x00"
        _sparse(db / fname, int(unpacked * share * scale), head)
    (db / "version.json").write_text(json.dumps({
        "date": "2024-01-19", "major": 5, "minor": 1, "type": kind, "doi": f"10.5281/zenodo.{spec['record']}",
    }, indent=2), encoding="utf-8")
    emit("\tremoving tarball file...")
    emit("Update AMRFinderPlus database...")
    _amrf_update_lines(emit, db / "amrfinderplus-db", indent="\t")
    emit("\tAMRFinderPlus database successfully updated")
    emit(f"Bakta database successfully downloaded, extracted and updated: {db}")
    return 0


def _checkm2_download(args: list[str], emit: Callable[[str], None]) -> int:
    path = Path(_arg(args, "--path"))
    stamp = lambda: datetime.now().strftime("[%m/%d/%Y %I:%M:%S %p]")  # noqa: E731
    emit(f"{stamp()} INFO: Command: Download database. Checking internal path information.")
    pause(demo_seconds(2))
    emit(f"{stamp()} INFO: Downloading CheckM2 database from Zenodo: "
         "https://zenodo.org/api/records/14897628/files/checkm2_database.tar.gz/content")
    if FAIL_NEXT_CHECKM2.is_set():
        FAIL_NEXT_CHECKM2.clear()
        _simulate_transfer(emit, 690_000_000, 24 * MB)
        emit(f"{stamp()} ERROR: Connection reset by peer while downloading checkm2_database.tar.gz")
        return 1
    _simulate_transfer(emit, 1_708_000_000, 24 * MB)
    emit(f"{stamp()} INFO: Extracting files from archive...")
    pause(demo_seconds(45))
    _sparse(path / "CheckM2_database" / "uniref100.KO.1.dmnd", 2_946_000_000, b"\x6d\x18\x05\xd4")
    emit(f"{stamp()} INFO: Verifying version and checksums...")
    pause(demo_seconds(25))
    emit(f"{stamp()} INFO: Verification success.")
    emit(f"{stamp()} INFO: Database check successful! (not written to the image: --no_write_json_db)")
    return 0


def _amrfinder_version(args: list[str], emit: Callable[[str], None]) -> int:
    db = Path(_arg(args, "--database"))
    emit("Software directory: '/opt/conda/envs/amrfinder/bin/'")
    emit(f"Software version: {AMRF_SOFTWARE}")
    pause(demo_seconds(4))
    version = db / "version.txt"
    if not version.is_file():
        emit(f"ERROR: Database directory does not exist or is not an AMRFinderPlus database: {db}")
        return 1
    fmt = (db / "database_format_version.txt").read_text(encoding="utf-8").strip() \
        if (db / "database_format_version.txt").is_file() else "?"
    emit(f"Database directory: '{db.resolve()}'")
    emit(f"Database version: {version.read_text(encoding='utf-8').strip()} (format {fmt})")
    return 0


# The tools image: missing until First-Time Setup pulls it (docker pull, then docker tag).
IMAGE_PRESENT = threading.Event()
# Set by the launcher: the next CheckM2 download breaks off half way (First-Time Setup's Retry).
FAIL_NEXT_CHECKM2 = threading.Event()
DEMO_REGISTRY = "registry.example.org/bactolomew"
_IMAGE_LAYERS = (("3c645031de29", 29), ("0a1d3b7c9f21", 412), ("5d2b8e44a9c0", 2_310), ("e71f0c8a6b13", 1_840),
                 ("9b04d6e2f7a5", 1_120), ("c38a5f1e0d92", 960), ("47e9b2d1c6f8", 655), ("b5c0a7e3d419", 88),
                 ("1f6d8c2b0e7a", 3))


def _docker_pull(args: list[str], emit: Callable[[str], None]) -> int:
    """docker pull without a terminal: one status line per layer and event."""
    ref = args[-1]
    name, _, tag = ref.rpartition(":")
    emit(f"{tag}: Pulling from {name.split('/', 1)[-1]}")
    for layer, _mb in _IMAGE_LAYERS:
        emit(f"{layer}: Pulling fs layer")
    total_mb = sum(mb for _l, mb in _IMAGE_LAYERS)
    wall = demo_seconds(total_mb / 45)  # 45 MB/s from the registry
    for layer, mb in _IMAGE_LAYERS:
        pause(wall * mb / total_mb)
        if _SHUTDOWN.is_set():
            return 1
        emit(f"{layer}: Verifying Checksum")
        emit(f"{layer}: Download complete")
        emit(f"{layer}: Pull complete")
    emit("Digest: sha256:7d4e0c1b9a8f6e5d4c3b2a19087f6e5d4c3b2a19087f6e5d4c3b2a19087f6e5d")
    emit(f"Status: Downloaded newer image for {ref}")
    emit(ref)
    return 0


def fake_tool_stream(cmd: list[str], log: Optional[Callable[[str], None]] = None, *,
                     cwd: Optional[Path] = None, env: Optional[dict] = None) -> int:
    """pipeline.ContainerRuntime._stream as DbDownload and the image setup call it: the
    command line, then the tool's output line by line; returns the exit status. The tool
    (or docker) is simulated."""
    emit = log or print
    emit("$ " + " ".join(str(c) for c in cmd))
    exe, args = _split_tool([str(c) for c in cmd])
    name = Path(exe).name
    try:
        if name == "docker" and args[:1] == ["pull"]:
            return _docker_pull(args, emit)
        if name == "docker" and args[:1] == ["tag"]:
            IMAGE_PRESENT.set()
            return 0
        if name in ("bakta", "checkm2") and args == ["--version"]:
            pause(demo_seconds(4))
            emit("bakta 1.11.4" if name == "bakta" else "1.1.0")
            return 0
        if name == "bakta_db" and args[:1] == ["download"]:
            return _bakta_db_download(args, emit)
        if name == "checkm2" and args[:1] == ["database"] and "--download" in args:
            return _checkm2_download(args, emit)
        if name == "amrfinder_update":
            _amrf_update_lines(emit, Path(_arg(args, "--database")))
            return 0
        if name == "amrfinder" and "--database_version" in args:
            return _amrfinder_version(args, emit)
        if name == "amrfinder_index":
            emit(f"Indexing {args[0] if args else ''}")
            pause(demo_seconds(30))
            emit("Indexing done")
            return 0
    except Exception as exc:  # noqa: BLE001 - reported like a failing tool
        emit(f"ERROR: {exc}")
        return 1
    emit(f"{DEMO_TAG}: no simulation for {name}; nothing was run")
    return 127


# --------------------------------------------------------------------------- #
# Simulated tools: chewBBACA (the subprocess.run of ChewbbacaSchemaBuild / CgmlstInstall)
# --------------------------------------------------------------------------- #
def write_chewbbaca_schema(out: Path, loci: dict[str, list[str]], *, table: int = 11) -> None:
    """The files of a chewBBACA schema: locus FASTAs, short/ representatives, .genes_list
    and .schema_config (pickles, as chewBBACA writes them)."""
    (out / "short").mkdir(parents=True, exist_ok=True)
    for name, seqs in loci.items():
        (out / f"{name}.fasta").write_text("".join(f">{name}_{i}\n{s}\n" for i, s in enumerate(seqs, 1)),
                                           encoding="utf-8")
        (out / "short" / f"{name}_short.fasta").write_text(f">{name}_1\n{seqs[0]}\n", encoding="utf-8")
    (out / ".genes_list").write_bytes(pickle.dumps([f"{n}.fasta" for n in loci]))
    (out / ".schema_config").write_bytes(pickle.dumps({
        "bsr": [0.6], "translation_table": [table], "chewBBACA_version": [CHEWBBACA_VERSION],
        "prodigal_training_file": [None], "minimum_locus_length": [0], "size_threshold": [None],
        "word_size": [None], "window_size": [None], "cluster_sim": [None],
        "representative_filter": [None], "intraCluster_filter": [None],
    }))


def _chewbbaca_banner() -> list[str]:
    return [f"chewBBACA version: {CHEWBBACA_VERSION}",
            "Authors: Rafael Mamede, Pedro Cerqueira, Mickael Silva, Joao Carrico, Mario Ramirez",
            "Github: https://github.com/B-UMMI/chewBBACA",
            "Documentation: https://chewbbaca.readthedocs.io/en/latest/index.html",
            "Contacts: imm-bioinfo@medicina.ulisboa.pt", ""]


_KILLED = -9


def _prep_external_schema(args: list[str], killed: Optional[threading.Event] = None) -> subprocess.CompletedProcess:
    """Adapts the loci in ten batches (their FASTAs appear one batch at a time, as the
    real PrepExternalSchema writes them); ``killed`` (the Popen stand-in's kill) ends it."""
    from pipeline.ChewbbacaSchemaBuild import is_forward_cds, read_fasta_file

    src = Path(_arg(args, "-g"))
    out = Path(_arg(args, "-o"))
    table = int(_arg(args, "--t", "11") or 11)
    fastas = sorted(src.glob("*.fasta"))
    invalid_alleles: list[str] = []
    invalid_loci: list[str] = []
    stats: list[str] = ["Gene\tTotal_alleles\tValid_alleles\tProportion_valid"]
    loci: dict[str, list[str]] = {}
    batches = 10
    for batch in range(batches):
        pause(demo_seconds(20 + 0.045 * len(fastas)) / batches)
        if killed is not None and killed.is_set():
            return subprocess.CompletedProcess(args, _KILLED, "", "killed")
        adapted: dict[str, list[str]] = {}
        for fasta in fastas[batch::batches]:
            records = read_fasta_file(fasta)
            valid = [seq for _h, seq in records if is_forward_cds(seq, table)]
            invalid_alleles += [f"{h}\tinternal stop codon" for h, seq in records if not is_forward_cds(seq, table)]
            stats.append(f"{fasta.stem}\t{len(records)}\t{len(valid)}\t{len(valid) / max(1, len(records)):.3f}")
            if valid:
                adapted[fasta.stem] = valid
            else:
                invalid_loci.append(fasta.stem)
        loci.update(adapted)
        write_chewbbaca_schema(out, adapted, table=table)
    write_chewbbaca_schema(out, dict(sorted(loci.items())), table=table)
    for suffix, lines in (("invalid_alleles.txt", invalid_alleles), ("invalid_loci.txt", invalid_loci),
                          ("summary_stats.tsv", stats)):
        (out.parent / f"{out.name}_{suffix}").write_text("".join(f"{ln}\n" for ln in lines), encoding="utf-8")
    stdout = "\n".join(_chewbbaca_banner() + [
        f"Number of loci: {len(fastas)}",
        "Adapting loci...",
        f"Number of invalid loci: {len(invalid_loci)}",
        f"Number of invalid alleles: {len(invalid_alleles)}",
        f"Successfully adapted {len(loci)}/{len(fastas)} loci present in the input schema.",
    ]) + "\n"
    return subprocess.CompletedProcess(args, 0, stdout, "")


def _allele_call(args: list[str], killed: Optional[threading.Event] = None) -> subprocess.CompletedProcess:
    from pipeline.ChewBBACAWorkflow import schema_representative_ids

    schema = Path(_arg(args, "-g"))
    out = Path(_arg(args, "-o"))
    loci = sorted(schema.glob("*.fasta"))
    seconds = demo_seconds(30 + 0.02 * len(loci))
    if killed is None:
        pause(seconds)
    elif killed.wait(seconds):
        return subprocess.CompletedProcess(args, _KILLED, "", "killed")
    reps = schema_representative_ids(schema)
    scores = {rep: (400, 812.0) for rep in reps}
    (schema / "short" / "self_scores").write_bytes(pickle.dumps(scores))
    out.mkdir(parents=True, exist_ok=True)
    (out / "results_statistics.tsv").write_text("FILE\tEXC\tINF\n synthetic\t%d\t0\n" % len(loci),
                                                encoding="utf-8")
    stdout = "\n".join(_chewbbaca_banner() + [
        "Number of inputs: 1",
        f"Number of loci: {len(loci)}",
        "Predicting gene coordinates with Prodigal (meta mode)...",
        "Translating coding sequences...",
        "Computing self-scores for the representative alleles...",
        f"Determining BLASTp self-score for {len(reps)} representatives...done.",
        f"Classified {len(loci)} loci (exact matches: {len(loci)}).",
        f"Wrote allele calling results to {out}",
    ]) + "\n"
    return subprocess.CompletedProcess(args, 0, stdout, "")


def _schema_metadata(args: list[str]) -> subprocess.CompletedProcess:
    schema = Path(args[-1])
    loci = sum(1 for p in schema.glob("*.fasta"))
    table = 11
    try:
        info = json.loads((schema / "schema_build_report" / "build_info.json").read_text(encoding="utf-8"))
        table = int(info.get("translation_table") or 11)
    except (OSError, ValueError, TypeError):
        pass
    pause(demo_seconds(3))
    payload = {"genes_list": loci if (schema / ".genes_list").is_file() else -1, "translation_table": table}
    return subprocess.CompletedProcess(args, 0, json.dumps(payload) + "\n", "")


class _FakePopen:
    """Popen of a simulated chewBBACA call (run on a thread): communicate(timeout) and
    kill(), what ChewbbacaSchemaBuild.run_watched uses."""

    def __init__(self, cmd, *args, **kwargs) -> None:  # noqa: ANN001
        self.args = cmd
        self.returncode: Optional[int] = None
        self._killed = threading.Event()
        self._result: Optional[subprocess.CompletedProcess] = None
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def _run(self) -> None:
        try:
            self._result = FakeSubprocess.run(self.args, killed=self._killed)
        except Exception as exc:  # noqa: BLE001 - the call's failure
            self._result = subprocess.CompletedProcess(self.args, 1, "", str(exc))

    def communicate(self, timeout: Optional[float] = None) -> tuple[str, str]:
        self._thread.join(timeout)
        if self._thread.is_alive():
            raise subprocess.TimeoutExpired(self.args, timeout)
        result = self._result or subprocess.CompletedProcess(self.args, 1, "", "")
        self.returncode = _KILLED if self._killed.is_set() else result.returncode
        return result.stdout, result.stderr

    def kill(self) -> None:
        self._killed.set()

    terminate = kill


class FakeSubprocess:
    """Stands in for the ``subprocess`` module of ChewbbacaSchemaBuild: run() and Popen
    simulate chewBBACA in the tools image; anything else comes from the real module."""

    CompletedProcess = subprocess.CompletedProcess
    Popen = _FakePopen

    @staticmethod
    def run(cmd, *args, killed: Optional[threading.Event] = None, **kwargs) -> subprocess.CompletedProcess:  # noqa: ANN001
        argv = [str(c) for c in cmd]
        shown = argv[argv.index("conda"):] if "conda" in argv else argv
        shown = [a if "\n" not in a else "<schema settings reader>" for a in shown]
        say("[simulated tools image] $ " + shlex.join(shown))
        if "PrepExternalSchema" in argv:
            return _prep_external_schema(argv, killed)
        if "AlleleCall" in argv:
            return _allele_call(argv, killed)
        if "-c" in argv:
            return _schema_metadata(argv)
        return subprocess.CompletedProcess(argv, 127, "", f"{DEMO_TAG}: blocked {argv[:3]}")

    def __getattr__(self, name: str) -> Any:
        return getattr(subprocess, name)


# --------------------------------------------------------------------------- #
# Guards: nothing real is ever started or downloaded
# --------------------------------------------------------------------------- #
_BLOCKED_EXES = frozenset({
    "docker", "podman", "singularity", "apptainer", "conda", "mamba", "micromamba",
    "chewBBACA.py", "chewbbaca", "bakta", "bakta_db", "checkm2", "amrfinder", "amrfinder_update",
    "amrfinder_index", "mob_init", "mlst", "curl", "wget", "git",
})
_RealPopen = subprocess.Popen


_RealSystem = os.system


def _blocked_reason(args: Any) -> Optional[str]:
    """Why a child process must not start in the demo, or None. Every word of the command
    is searched, so a tool behind ``sh -c``, ``bash -lc``, ``env`` or ``timeout`` is found
    too; a child python is refused because it would run without the demo's simulation."""
    if isinstance(args, (str, bytes, os.PathLike)):
        words = [os.fsdecode(args)]
    else:
        words = [os.fsdecode(a) if isinstance(a, (bytes, os.PathLike)) else str(a) for a in args]
    for word in words:
        for token in re.split(r"[\s;&|()`$<>'\"]+", word):
            exe = Path(token).name
            if exe in _BLOCKED_EXES:
                return f"blocked a real tool call ({exe})"
            if re.fullmatch(r"python[\d.]*", exe) or (token and token == sys.executable):
                return f"blocked a child python ({token}): it would run without the demo's simulation"
    return None


class _GuardedPopen(_RealPopen):  # type: ignore[misc, valid-type]
    def __init__(self, args, *a, **kw):  # noqa: ANN001
        reason = _blocked_reason(args)
        if reason:
            raise RuntimeError(f"{DEMO_TAG}: {reason}")
        super().__init__(args, *a, **kw)


def _guarded_system(command: str) -> int:
    reason = _blocked_reason(command)
    if reason:
        raise RuntimeError(f"{DEMO_TAG}: {reason}")
    return _RealSystem(command)


# --------------------------------------------------------------------------- #
# The demo folder
# --------------------------------------------------------------------------- #
class DemoEnv:
    """A fresh project root: config, DB/, runs/, cache/, downloads and QSettings."""

    def __init__(self, root: Path) -> None:
        self.root = Path(root).resolve()
        self.db = self.root / "DB"
        self.cgmlst = self.db / "cgmlst"
        self.catalog = self.cgmlst / "scheme_species_catalog.json"
        self.downloads = self.root / "Downloads"
        self.settings_dir = self.root / "qsettings"
        self.tmp = self.root / "tmp"

        def root_fn() -> Path:  # the project_root every app module gets in the demo
            return self.root

        root_fn._demo_root = True  # type: ignore[attr-defined]
        self.root_fn = root_fn

    # ---- build ----
    def build(self) -> "DemoEnv":
        for path in (self.root / "runs", self.root / "cache", self.tmp, self.settings_dir,
                     self.db / "metadata", self.cgmlst):
            path.mkdir(parents=True, exist_ok=True)
        self._catalog()
        for slug in PREINSTALLED:
            self._installed_scheme(ORG[slug])
        self._config()
        self._downloads()
        return self

    def _catalog(self) -> None:
        from pipeline.Setup import CGMLST_SCHEME_CATALOG_SOURCE

        shutil.copyfile(CGMLST_SCHEME_CATALOG_SOURCE, self.catalog)  # what First-Time Setup copies

    def _installed_scheme(self, scheme: OrgSchemeData) -> None:
        out = self.cgmlst / scheme.schema_id
        loci = {name: scheme.alleles(name, max_alleles=1) for name in scheme.locus_names()}
        write_chewbbaca_schema(out, loci)
        report = out / "schema_build_report"
        report.mkdir(exist_ok=True)
        (report / "build_info.json").write_text(json.dumps({
            "built_at": "2026-10-09T07:52:31+00:00",
            "builder": "chewBBACA PrepExternalSchema",
            "chewbbaca_version": CHEWBBACA_VERSION,
            "bsr": 0.6,
            "translation_table": 11,
            "source": scheme.scheme_url,
            "source_format": "zip_archive+targets",
            "expected_loci": scheme.loci,
            "loci": scheme.loci,
            "alleles": scheme.loci,
            "representatives": scheme.loci,
        }, indent=2) + "\n", encoding="utf-8")

    def _config(self) -> None:
        from pipeline.ToolEnv import configured_fastgear_settings

        schemas = {p.name: str(p) for p in sorted(self.cgmlst.iterdir())
                   if p.is_dir() and (p / ".genes_list").is_file()}
        first = next(iter(schemas.values()), "")
        config = {
            "runtime": "docker",
            "container": {"image": "liviurotiul/bactolomew:0.3.0", "registry": DEMO_REGISTRY, "sif_path": "",
                          "build_on_missing": False},
            "databases": {
                "root_dir": str(self.db),
                "marker_sequences_dir": str(self.db / "seq" / "markers"),
                "compound_sequences_dir": str(self.db / "seq" / "compound"),
                "marker_catalog_path": str(self.db / "catalogs" / "marker_db_catalog.json"),
                "compound_catalog_path": str(self.db / "catalogs" / "compound_db_catalog.json"),
                "metadata_dir": str(self.db / "metadata"),
                "cgmlst_dir": str(self.cgmlst),
                "cgmlst_scheme_catalog_path": str(self.catalog),
            },
            "chewbbaca": {
                "auto_select_from_mlst": True,
                "available_schemas": schemas,
                "catalog_path": str(self.catalog),
                "schema_dir": first,
                "genes_list": str(Path(first) / ".genes_list") if first else "",
                "loci_list": "",
                "training_file": "",
            },
            "fastgear": configured_fastgear_settings({}),
            "resources": {"max_threads": 16, "max_ram_gb": 0},
            "dev_mode": False,
        }
        (self.root / "bactolomew_config.json").write_text(json.dumps(config, indent=2, sort_keys=True),
                                                          encoding="utf-8")

    def _downloads(self) -> None:
        """What a user's Downloads folder holds after fetching schemes by hand in a browser."""
        dl = self.downloads
        dl.mkdir(parents=True, exist_ok=True)
        (dl / "null.zip").write_bytes(org_alleles_zip("Smarcescens"))
        (dl / "null (1).zip").write_bytes(org_alleles_zip("Efaecalis"))
        (dl / "Scapitis_cgMLST_alleles.zip").write_bytes(org_alleles_zip("Scapitis"))
        (dl / "Abaumannii_cgMLST_alleles.zip").write_bytes(org_alleles_zip("Abaumannii"))
        (dl / "Cdifficile_cgMLST_alleles.zip").write_bytes(_truncated(org_alleles_zip("Cdifficile"), 0.55))
        lmo = dl / "Lmonocytogenes"
        lmo.mkdir(exist_ok=True)
        (lmo / "alleles.zip").write_bytes(org_alleles_zip("Lmonocytogenes"))
        (lmo / "targets.csv").write_text(org_targets_csv("Lmonocytogenes"), encoding="utf-8")
        (dl / "README.txt").write_text("Saved from https://www.cgmlst.org/ncs/ with a web browser.\n",
                                       encoding="utf-8")

    def describe_downloads(self) -> list[str]:
        notes = {
            "null.zip": "browser name of a cgMLST.org download (Serratia marcescens, found by its loci)",
            "null (1).zip": "the second one (Enterococcus faecalis, found by its loci)",
            "Scapitis_cgMLST_alleles.zip": "Staphylococcus capitis, found by its name",
            "Abaumannii_cgMLST_alleles.zip": "Acinetobacter baumannii: already installed",
            "Cdifficile_cgMLST_alleles.zip": "Clostridioides difficile, truncated download (broken)",
            "Lmonocytogenes": "folder with targets.csv + alleles.zip: already installed",
            "README.txt": "not a schema (ignored)",
        }
        real = {"null.zip": "Smarcescens", "null (1).zip": "Efaecalis", "Scapitis_cgMLST_alleles.zip": "Scapitis",
                "Abaumannii_cgMLST_alleles.zip": "Abaumannii", "Lmonocytogenes": "Lmonocytogenes"}
        lines = ["  (sizes: the demo's small stand-in, then the real cgMLST.org download it stands for)"]
        for path in sorted(self.downloads.iterdir(), key=lambda p: p.name.lower()):
            size = sum(p.stat().st_size for p in path.rglob("*")) if path.is_dir() else path.stat().st_size
            org = ORG.get(real.get(path.name, ""))
            real_size = f"(real {human_size(org.size)})" if org is not None and org.size else ""
            lines.append(f"  {path.name + ('/' if path.is_dir() else ''):30s} {human_size(size):>9} {real_size}")
            if notes.get(path.name):
                lines.append(f"      {notes[path.name]}")
        return lines

    # ---- state ----
    def config(self) -> dict:
        try:
            return json.loads((self.root / "bactolomew_config.json").read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return {}

    def installed_ids(self) -> set[str]:
        return {p.name for p in self.cgmlst.iterdir() if p.is_dir() and (p / ".genes_list").is_file()}


def _truncated(data: bytes, fraction: float) -> bytes:
    """``data`` cut off inside a member's compressed bytes, past ``fraction`` of it: a
    download that stopped half way."""
    pos = data.find(b"PK\x03\x04", int(len(data) * fraction))
    if pos < 0:
        return data[: int(len(data) * fraction)]
    name_len, extra_len = int.from_bytes(data[pos + 26:pos + 28], "little"), int.from_bytes(data[pos + 28:pos + 30], "little")
    return data[: pos + 30 + name_len + extra_len + 40]


def isolate_qsettings(path: Path) -> None:
    """Every QSettings format and scope in the demo folder (tests/_qt_test_env style)."""
    from PyQt6 import QtCore

    fmt, scope = QtCore.QSettings.Format, QtCore.QSettings.Scope
    QtCore.QSettings.setDefaultFormat(fmt.IniFormat)
    for each_format in (fmt.IniFormat, fmt.NativeFormat):
        for each_scope in (scope.UserScope, scope.SystemScope):
            QtCore.QSettings.setPath(each_format, each_scope, str(path))


# --------------------------------------------------------------------------- #
# Wiring the simulated layer into the app's modules
# --------------------------------------------------------------------------- #
def install_backend(env: DemoEnv) -> None:
    """Point the app at the demo folder and replace the network, tool and Docker layer."""
    import annotation_explorer_polars.cache as cache_mod
    import pipeline.Setup as fts
    import annotation_explorer_polars.startup_dialog as sd
    import pipeline.CgmlstInstall as ci
    import pipeline.ChewBBACAWorkflow as cw
    import pipeline.ChewbbacaSchemaBuild as csb
    import pipeline.ContainerRuntime as cr
    import pipeline.DbDownload as dbd
    import pipeline.ToolEnv as te

    env.tmp.mkdir(parents=True, exist_ok=True)
    tempfile.tempdir = str(env.tmp)
    os.environ.pop("BACTOLOMEW_NATIVE", None)
    # Every module of the app that resolves the project folder (and so the config, DB/,
    # runs/ and cache/) by its own project_root, or keeps a DB/ path as a constant: the
    # ones loaded now, and through the import hook every one imported later.
    for mod_name, module in list(sys.modules.items()):
        _redirect_module(env, mod_name, module)
    sys.meta_path.insert(0, _RedirectFinder(env))

    urllib.request.urlopen = fake_urlopen  # cgMLST.org and NCBI, simulated
    subprocess.Popen = _GuardedPopen  # type: ignore[misc]
    os.system = _guarded_system
    dbd._stream = fake_tool_stream
    cr._stream = fake_tool_stream  # docker pull / docker tag of the tools image
    csb.subprocess = FakeSubprocess()
    cr.docker_available = lambda: True
    cr.docker_status = lambda exe="docker": cr.RuntimeStatus("docker", True)
    cr.runtime_status = lambda runtime: cr.RuntimeStatus(runtime, True)
    cr.engine_info = lambda runtime: cr.EngineInfo("docker")  # plain rootful Docker
    cr.docker_image_present = lambda image, engine="docker": IMAGE_PRESENT.is_set()
    fts._docker_client_available = lambda docker_bin: True
    fts._docker_image_available = lambda docker_bin, image: True

    for module, attr in ((te, "project_root"), (cw, "project_root"), (dbd, "project_root"),
                         (cr, "project_root"), (cache_mod, "project_root"), (fts, "project_root"),
                         (sd, "_project_root")):
        if Path(getattr(module, attr)()) != env.root:
            raise SystemExit(f"demo: {module.__name__}.{attr} still points at the repository")
    if Path(ci.default_dest()) != env.cgmlst:
        raise SystemExit(f"demo: the cgMLST folder is {ci.default_dest()}, not the demo's")
    import pipeline.DbCatalog as dbc

    if Path(dbc.MARKER_DB_CATALOG_PATH).is_relative_to(REPO_ROOT):
        raise SystemExit("demo: pipeline.DbCatalog still points at the repository's DB/")


_APP_PACKAGES = ("pipeline", "annotation_explorer_polars")
_REPO_DB = REPO_ROOT / "DB"


def _redirect_module(env: DemoEnv, mod_name: str, module: Any) -> None:
    """Point one app module at the demo folder: its project_root functions and the DB/
    paths it holds as module constants (pipeline.DbCatalog's and the copies other
    modules import from it)."""
    if module is None or mod_name.split(".")[0] not in _APP_PACKAGES:
        return
    for attr in ("project_root", "_project_root"):
        fn = getattr(module, attr, None)
        if callable(fn) and not getattr(fn, "_demo_root", False):
            try:
                if Path(fn()) == REPO_ROOT:
                    setattr(module, attr, env.root_fn)
            except TypeError:
                pass
    for attr, value in list(vars(module).items()):
        if attr.isupper() and isinstance(value, Path):
            if value == REPO_ROOT and attr == "PROJECT_ROOT" and mod_name == "pipeline.DbCatalog":
                setattr(module, attr, env.root)
            elif value == _REPO_DB or value.is_relative_to(_REPO_DB):
                setattr(module, attr, env.db / value.relative_to(_REPO_DB))


class _RedirectFinder:
    """An import hook: every app module imported after the demo started is redirected
    (_redirect_module) as soon as it has run, so the demo never depends on which modules
    happened to be loaded first."""

    def __init__(self, env: DemoEnv) -> None:
        self.env = env

    def find_spec(self, name: str, path: Any = None, target: Any = None) -> Any:
        if name.split(".")[0] not in _APP_PACKAGES:
            return None
        from importlib.machinery import PathFinder

        spec = PathFinder.find_spec(name, path)
        loader = getattr(spec, "loader", None)
        if loader is None or not hasattr(loader, "exec_module"):
            return spec
        original, env = loader.exec_module, self.env

        def exec_module(module: Any) -> None:
            original(module)
            _redirect_module(env, name, module)

        loader.exec_module = exec_module
        return spec

    def invalidate_caches(self) -> None:
        pass


# --------------------------------------------------------------------------- #
# A simulated pipeline run (writes run_state.json like pipeline.Runner)
# --------------------------------------------------------------------------- #
class _Halt(Exception):
    pass


# Steps done before the demo run reaches the phylogeny: (name, real minutes, per genome).
_DONE_BEFORE = (
    ("assembly_qc", 4.2, True), ("bakta", 41.0, True), ("mlst", 0.7, True),
    ("amrfinder_mutations", 3.1, True), ("mobsuite", 9.4, True), ("bakta_evidence", 0.3, False),
    ("panaroo", 14.2, False),
    ("core_alignment", 2.3, False), ("fastgear_clades", 11.0, False), ("gubbins", 6.1, False),
    ("core_tree", 4.4, False),
)
# Real minutes of the steps after the core tree that the demo passes through quickly.
_LATER = {"phylogeny": 23.0, "diamond_markers": 6.5, "diamond_reference_loci": 4.2, "collect_annotations": 1.4}


def _load_prctl() -> Optional[Callable[..., int]]:
    try:
        import ctypes

        return ctypes.CDLL("libc.so.6", use_errno=True).prctl
    except Exception:  # noqa: BLE001 - not Linux
        return None


_PRCTL = _load_prctl()  # loaded here: nothing may be imported between fork and exec


_ALLELE_CALL_MINUTES = 7.5  # chewBBACA AlleleCall of 12 genomes on 8 threads


def _real_install_minutes(schema_id: str) -> float:
    """Roughly how long a cgMLST.org install takes for real: the download at about 4 MB/s,
    PrepExternalSchema at about 4 loci a second, and the test allele call."""
    scheme = org_by_schema_id(schema_id)
    if scheme is None:
        return 10.0
    size = scheme.size or scheme.loci * 12_000
    return size / (4 * MB) / 60 + scheme.loci / 4 / 60 + 1.5


def _die_with_parent() -> None:
    """In the stand-in run process before it starts: it gets SIGKILL when the thread that
    started it ends (the run's thread; also when the demo is killed or crashes), so no
    ``sleep`` outlives the demo (Linux PR_SET_PDEATHSIG; elsewhere shutdown() kills it on
    a normal exit)."""
    if _PRCTL is not None:
        _PRCTL(1, int(signal.SIGKILL))  # PR_SET_PDEATHSIG = 1


class SimulatedRun(threading.Thread):
    """A run of 12 genomes of one species, from the Snippy tree onwards. The chewBBACA
    step runs ChewBBACAWorkflow.prepare_chewbbaca_schema, the real code run_chewbbaca runs
    before AlleleCall (the schema auto-selection and the cgMLST.org install of a missing
    scheme), with the Runner's status callback, and fails like the Runner's step when it
    raises."""

    def __init__(self, env: DemoEnv, run_dir: Path, *, label: str, resume: bool = False) -> None:
        super().__init__(daemon=True, name=f"demo-run-{label}")
        self.env = env
        self.run_dir = Path(run_dir)
        self.label = label
        self.resume = resume
        self.child: Optional[subprocess.Popen] = None
        self.logger = logging.getLogger(f"bactolomew.demo.run.{label}")
        self.logger.setLevel(logging.INFO)
        if not self.logger.handlers:
            handler = logging.StreamHandler(sys.stdout)
            handler.setFormatter(logging.Formatter(f"[run {label}] %(asctime)s %(levelname)s %(message)s",
                                                   "%H:%M:%S"))
            self.logger.addHandler(handler)
            self.logger.propagate = False
        state = json.loads((self.run_dir / "run_state.json").read_text(encoding="utf-8"))
        self.state = state
        self.steps: dict = state["steps"]
        self.warnings: list = state.get("warnings") or []
        from pipeline.RunSteps import StepProgress

        # The Runner's step progress (its run_state.json throttle and the install's status).
        self.state_lock = threading.RLock()
        self.progress = StepProgress(self.steps, self.state_lock, lambda: self._flush("running", "chewbbaca"),
                                     interval=0.5)

    # ---- files ----
    def _write(self, name: str, payload: dict) -> None:
        tmp = self.run_dir / f".{name}.tmp"
        tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        os.replace(tmp, self.run_dir / name)

    def _flush(self, status: str, current: Optional[str], *, error: str = "", report: bool = False) -> None:
        from pipeline.RunSteps import planned_steps

        with self.state_lock:
            self.progress.written()
            self._flush_locked(status, current, error=error, report=report, planned=planned_steps())

    def _flush_locked(self, status: str, current: Optional[str], *, error: str, report: bool,
                      planned: list) -> None:
        payload = dict(self.state)
        payload.update({
            "updated_at": _now_iso(), "status": status, "current_step": current,
            "planned_steps": planned, "steps": self.steps, "warnings": self.warnings,
            "outputs_validated": status == "success", "report_written": report,
            "stop_requested": (self.run_dir / "run.stop").exists(),
        })
        if status in {"success", "failed", "stopped", "interrupted", "partial"}:
            payload["end_time"] = _now_iso()
        if error:
            payload["error"] = error
        self.state = payload
        self._write("run_state.json", payload)

    def _lock(self) -> None:
        # The run's process: a sleeping child, so the start window's Force Stop ends this
        # demo run (as it would end the Runner), never the demo itself.
        self.child = _RealPopen(["sleep", "86400"], stdin=subprocess.DEVNULL,
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                preexec_fn=_die_with_parent)
        self._write("run.lock.json", {"run_id": self.state.get("run_id"), "pid": self.child.pid,
                                      "host": socket.gethostname(), "start_time": self.state.get("start_time")})

    def _unlock(self) -> None:
        (self.run_dir / "run.lock.json").unlink(missing_ok=True)
        if self.child is not None and self.child.poll() is None:
            self.child.kill()

    def stop_process(self) -> None:
        """At demo exit: the run ends as if its process died."""
        if self.child is not None and self.child.poll() is None:
            self.child.kill()

    def _check(self, *, stop: bool = False) -> None:
        if _SHUTDOWN.is_set() or (self.child is not None and self.child.poll() is not None):
            raise _Halt("killed")  # run_state.json stays "running": the GUI shows it interrupted
        if stop and (self.run_dir / "run.stop").exists():  # the Runner stops between steps
            self.logger.info("Stop requested: stopping before the next step.")
            (self.run_dir / "run.stop").unlink(missing_ok=True)
            self._flush("stopped", None, error="Run stopped by user")
            self._unlock()
            raise _Halt("stopped")

    def _wait(self, seconds: float) -> None:
        end = time.monotonic() + seconds
        while time.monotonic() < end:
            self._check()
            pause(min(0.2, end - time.monotonic()))

    def _start(self, name: str) -> None:
        self.steps[name] = {"status": "running", "start_time": _now_iso(), "outputs_validated": False,
                            "validation_errors": []}
        self._flush("running", name)

    def _end(self, name: str, status: str, *, real_minutes: float = 0.0, **info: Any) -> None:
        step = self.steps.get(name, {})
        step.update({"status": status, "end_time": _now_iso(), "outputs_validated": status == "success", **info})
        if real_minutes:  # a step the demo only shows in passing keeps a realistic duration
            step["start_time"] = _now_iso(-real_minutes * 60)
        self.steps[name] = step
        self._flush("running", None)

    # ---- run ----
    def run(self) -> None:
        try:
            self._lock()
            self._flush("running", None)
            self.logger.info("%s run %s (%s genomes)", "Resuming" if self.resume else "Running",
                             self.state.get("run_id"), len(self.state["inputs"]["genomes"]))
            from pipeline.RunSteps import PIPELINE_STEPS

            for name, _label in PIPELINE_STEPS:
                if str((self.steps.get(name) or {}).get("status")) == "success":
                    continue
                self._check(stop=True)
                if name == "chewbbaca":
                    self._chewbbaca()
                else:
                    self._start(name)
                    self.logger.info("Step %s", name)
                    self._wait(3.0 / max(0.05, SPEED))
                    self._end(name, "success", real_minutes=_LATER.get(name, 3.0))
            report = {"status": "success", "outputs_validated": True, "run_id": self.state.get("run_id"),
                      "sample_set_id": self.state.get("sample_set_id"), "end_time": _now_iso(),
                      "steps": self.steps, "warnings": self.warnings, "validation_errors": []}
            self._write("run_report.json", report)
            self._flush("success", None, report=True)
            self.logger.info("Run completed%s.", " with warnings: " + "; ".join(w["code"] for w in self.warnings)
                             if self.warnings else "")
            self._unlock()
        except _Halt as halt:
            self.logger.info("Run %s.", "stopped" if str(halt) == "stopped" else "ended (process gone)")
        except Exception as exc:  # noqa: BLE001 - shown as a failed run
            from pipeline.ChewBBACAWorkflow import CgmlstSchemaMissing

            self.logger.error("Run failed: %s", exc)
            if not isinstance(exc, CgmlstSchemaMissing):  # an expected stop says what to do
                traceback.print_exc()
            self._flush("failed", None, error=str(exc))
            self._unlock()

    def _chewbbaca(self) -> None:
        from pipeline.ChewBBACAWorkflow import CgmlstSchemaMissing, prepare_chewbbaca_schema
        from pipeline.ChewbbacaSchemaBuild import Cancelled

        out = self.run_dir
        self._start("chewbbaca")
        logger = self.logger
        before = self.env.installed_ids()
        try:
            # ChewBBACAWorkflow.run_chewbbaca up to AlleleCall (the run's config adds nothing);
            # a Stop during the install stops the run, as in the Runner.
            settings = prepare_chewbbaca_schema(out, logger, threads=8,
                                                status=self.progress.step_status("chewbbaca"),
                                                should_cancel=lambda: (out / "run.stop").exists())
        except Cancelled:
            self._check(stop=True)
            raise
        except Exception as exc:
            # pipeline.Runner: the step fails and so does the run; Resume runs it again.
            logger.error("chewBBACA workflow failed: %s", exc)
            code = {"code": exc.code} if isinstance(exc, CgmlstSchemaMissing) else {}
            self._end("chewbbaca", "failed", error=str(exc), **code)
            say("[demo] The run stops at the chewBBACA step. Run progress shows the reason; "
                "Resume runs the step again.")
            raise
        installed = sorted(self.env.installed_ids() - before)
        install_minutes = sum(_real_install_minutes(sid) for sid in installed)
        schema = Path(settings["schema_dir"])
        loci = sum(1 for _ in schema.glob("*.fasta"))
        genomes = self.state["inputs"]["genomes"]
        logger.info("[simulated tools image] chewBBACA.py AlleleCall -i genomes.txt -g %s -o chewbbaca/out "
                    "--cpu 8", schema)
        for i in range(1, len(genomes) + 1):
            self._wait(demo_seconds(25) / 3)
            logger.info("chewBBACA AlleleCall: %d/%d genomes", i, len(genomes))
        chew = out / "chewbbaca" / "out"
        chew.mkdir(parents=True, exist_ok=True)
        names = [p.stem for p in sorted(schema.glob("*.fasta"))]
        rows = ["FILE\t" + "\t".join(names)] + [
            Path(g).stem + "\t" + "\t".join(str(1 + (i * 7 + j) % 4) for j in range(len(names)))
            for i, g in enumerate(genomes)]
        (chew / "results_alleles.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
        for fname in ("results_statistics.tsv", "loci_summary_stats.tsv", "logging_info.txt"):
            (chew / fname).write_text("demo\n", encoding="utf-8")
        logger.info("chewBBACA called %d loci in %d genomes with %s.", loci, len(genomes), schema.name)
        # The row shows a realistic duration: AlleleCall, plus the install when it ran here.
        self._end("chewbbaca", "success", real_minutes=_ALLELE_CALL_MINUTES + install_minutes)
        if install_minutes:
            say(f"[demo] For real the install inside this run adds about {install_minutes:.0f} min to the "
                f"chewBBACA step.")


def create_run(env: DemoEnv, species: str, mlst_scheme: str, *, tag: str) -> Path:
    """A run folder in the state a real run has when its Snippy tree is being built: 12
    genomes of ``species`` typed by MLST as ``mlst_scheme``."""
    from annotation_explorer_polars.startup_dialog import _compute_run_id, _new_attempt_id

    short = "".join(w[0] for w in species.split()[:1]) + species.split()[-1][:3]
    genome_dir = env.root / "genomes" / f"{short}_{tag}"
    genome_dir.mkdir(parents=True, exist_ok=True)
    genomes = []
    for i in range(1, 13):
        path = genome_dir / f"{short}_{tag}_isolate_{i:02d}.fasta"
        path.write_text(f">contig_1\n{'ATGAAACCCGGGTTTTAA' * 4}\n", encoding="utf-8")
        genomes.append(path)
    sample_set = _compute_run_id(genomes)
    attempt = _new_attempt_id(sample_set)
    run_dir = env.root / "runs" / sample_set / attempt
    run_dir.mkdir(parents=True)
    alias = mlst_scheme
    (run_dir / "mlst").mkdir()
    (run_dir / "mlst" / "results.json").write_text(json.dumps([
        {"sample": g.stem, "scheme": alias, "st": str(17 + (i % 3) * 63)} for i, g in enumerate(genomes)
    ], indent=1), encoding="utf-8")
    (run_dir / "run_config.json").write_text(json.dumps({
        "run_id": attempt, "sample_set_id": sample_set, "output_dir": str(run_dir),
        "inputs": {"genomes": [str(g) for g in genomes]},
    }, indent=2), encoding="utf-8")
    total = sum(m for _n, m, _g in _DONE_BEFORE) * 60
    t = -total - 20
    steps = {}
    for name, minutes, per_genome in _DONE_BEFORE:
        step = {"status": "success", "start_time": _now_iso(t), "end_time": _now_iso(t + minutes * 60),
                "outputs_validated": True, "validation_errors": []}
        if per_genome:
            step["progress"] = {"done": 12, "total": 12, "unit": "genomes"}
        steps[name] = step
        t += minutes * 60
    state = {"start_time": _now_iso(-total - 20), "status": "running", "sample_set_id": sample_set,
             "run_id": attempt, "inputs": {"genomes": [str(g) for g in genomes]}, "steps": steps,
             "warnings": [], "output_dir": str(run_dir)}
    (run_dir / "run_state.json").write_text(json.dumps(state, indent=2), encoding="utf-8")
    return run_dir


def write_no_scheme_source(env: DemoEnv) -> Path:
    """The schema files of NO_SCHEME_SPECIES a user got elsewhere: a folder of per-locus
    allele FASTAs (outside Downloads, so step 3's listing does not change)."""
    species, _mlst, first = NO_SCHEME_SPECIES
    folder = env.root / "schemes_from_a_colleague" / species.replace(" ", "_")
    folder.mkdir(parents=True, exist_ok=True)
    for locus in _locus_names(first, NO_SCHEME_LOCI):
        seqs = _alleles(f"noscheme/{locus}", 3)
        (folder / f"{locus}.fasta").write_text("".join(f">{i}\n{q}\n" for i, q in enumerate(seqs, 1)),
                                              encoding="utf-8")
    return folder


def pick_run_species(env: DemoEnv, exclude: set[str]) -> Optional[dict]:
    """A catalog species with an MLST scheme, a cgMLST.org scheme and no installed schema."""
    from pipeline.CgmlstInstall import load_catalog

    entries = {e["schema_id"]: e for e in load_catalog(env.catalog)}
    installed = env.installed_ids()
    order = list(RUN_SPECIES_ORDER) + sorted(entries)
    for sid in order:
        entry = entries.get(sid)
        if entry and entry["mlst_schemes"] and sid not in installed and sid not in exclude \
                and org_by_schema_id(sid) is not None:
            return entry
    return None


# --------------------------------------------------------------------------- #
# Launcher window
# --------------------------------------------------------------------------- #
def run_launcher(env: DemoEnv) -> int:
    from annotation_explorer_polars.theme import apply_application_theme

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
    apply_application_theme(app)
    launcher = Launcher(env)
    launcher.show()
    signal.signal(signal.SIGINT, lambda *_a: app.quit())
    heartbeat = QtCore.QTimer()
    heartbeat.timeout.connect(lambda: None)  # lets Python see Ctrl+C
    heartbeat.start(300)
    code = app.exec()
    launcher.shutdown()
    return code


class _Bridge(QtCore.QObject):
    line = QtCore.pyqtSignal(str)
    finished = QtCore.pyqtSignal(object, object)  # (callback, result-or-exception)


class _DemoWindows(QtCore.QObject):
    """Marks every window the app opens as a demo window, and fills in the downloads folder
    when "Add schema..." opens."""

    def __init__(self, env: DemoEnv, launcher: QtWidgets.QWidget) -> None:
        super().__init__()
        self.env = env
        self.launcher = launcher

    def eventFilter(self, obj: QtCore.QObject, event: QtCore.QEvent) -> bool:  # noqa: N802
        if event.type() == QtCore.QEvent.Type.Show and isinstance(obj, QtWidgets.QDialog) \
                and obj is not self.launcher:
            title = obj.windowTitle()
            if DEMO_TAG not in title:
                obj.setWindowTitle(f"{title} - {DEMO_TAG}" if title else DEMO_TAG)
            self._add_banner(obj)
            if obj.objectName() == "AddSchemaDialog":
                edit = getattr(obj, "source_edit", None)
                if edit is not None and not edit.text():
                    edit.setText(str(self.env.downloads))
        return False

    @staticmethod
    def _add_banner(dialog: QtWidgets.QDialog) -> None:
        """A one-line demo notice at the top of an app window (not on message boxes, which
        keep their title suffix only)."""
        if dialog.property("demoBanner") or isinstance(dialog, (QtWidgets.QMessageBox, QtWidgets.QFileDialog)):
            return
        layout = dialog.layout()
        if not isinstance(layout, QtWidgets.QBoxLayout) or layout.direction() != QtWidgets.QBoxLayout.Direction.TopToBottom:
            return
        dialog.setProperty("demoBanner", True)
        banner = QtWidgets.QLabel(f"{DEMO_TAG} (simulated network and tools)")
        banner.setObjectName("DemoBanner")
        banner.setWordWrap(True)
        banner.setStyleSheet(
            f"QLabel#DemoBanner {{ background: {COLORS['warning_soft']}; border: 1px solid "
            f"{COLORS['warning_border']}; border-radius: 6px; color: {COLORS['text']}; padding: 3px 10px; }}")
        layout.insertWidget(0, banner)


_STEPS = (
    ("First-run setup: tools and databases",
     "Opens First-Time Setup. It runs by itself: the tools image is pulled, the Bakta (full) and CheckM2 "
     "databases download into .partial folders, are verified and renamed into place, and a quick check "
     "writes the config. Open it again to see everything found as done. \"With a failed download\" breaks "
     "the CheckM2 download off once (Retry). AMRFinderPlus is command line only."),
    ("cgMLST schemes from cgMLST.org",
     "Opens Settings: click \"Download from cgMLST.org...\", tick a few species and press Download and "
     "install. Terms are asked once; each scheme is downloaded, adapted, tested and registered."),
    ("Import downloaded schemes",
     "Opens Settings: click \"Add schema...\" (source preset to a browser Downloads folder) and Import. "
     "Each item is imported, skipped as already installed, or fails with its reason."),
    ("A run on a species with no installed scheme",
     "A simulated 12-genome run in the real start window and Run progress. \"Missing scheme\": the "
     "scheme installs from cgMLST.org inside the chewBBACA step, with its progress on the step's row and "
     "the run card. Without the cgMLST.org terms accepted the app asks for them first (\"Terms not "
     "accepted\" forgets them); \"Not now\" stops the run at chewBBACA until they are accepted and the run "
     "resumed. \"No scheme\": a species cgMLST.org has no scheme for stops at chewBBACA; import its "
     "schema files, then Resume in the start window."),
    ("The same from the command line",
     "The real command-line tools with --dry-run (the exact tool commands), --list and --status, "
     "against the demo folder."),
)


class Launcher(QtWidgets.QDialog):
    def __init__(self, env: DemoEnv) -> None:
        super().__init__()
        self.env = env
        self.setObjectName("DemoLauncher")
        self.setWindowTitle(f"BACTolomew install demo - {DEMO_TAG}")
        self.resize(1240, 900)
        self.bridge = _Bridge()
        self.bridge.line.connect(self._append)
        self.bridge.finished.connect(lambda cb, res: cb(res))
        self._busy = False
        self._runs: dict[str, SimulatedRun] = {}
        self._run_count = 0
        self._startup = None
        self._step_buttons: list[QtWidgets.QPushButton] = []
        self._build_ui()
        self._windows = _DemoWindows(env, self)
        QtWidgets.QApplication.instance().installEventFilter(self._windows)
        sys.stdout = _ConsoleTee(sys.__stdout__, self.bridge.line.emit)
        sys.stderr = _ConsoleTee(sys.__stderr__, self.bridge.line.emit)
        self._patch_start_window()
        # main() printed the folder on the terminal already: these go to the console only.
        self._append(f"{DEMO_TAG}. Demo folder: {env.root}")
        say(f"Speed {SPEED:g}: a step that takes T seconds for real takes about sqrt(T)/{SPEED:g} s here.")

    # ---- UI ----
    def _build_ui(self) -> None:
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(22, 18, 22, 14)
        outer.setSpacing(8)
        title = QtWidgets.QLabel("BACTolomew install demo")
        title.setProperty("role", "windowTitle")
        outer.addWidget(title)
        banner = QtWidgets.QLabel(
            f"<b>{DEMO_TAG}.</b> The app's real windows run on a simulated cgMLST.org, NCBI and tools "
            f"image; nothing is downloaded and no container is started. Everything is written to "
            f"{self.env.root} and removed when the demo closes. Do the steps in order.")
        banner.setWordWrap(True)
        banner.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        banner.setStyleSheet(
            f"QLabel {{ background: {COLORS['warning_soft']}; border: 1px solid {COLORS['warning_border']};"
            f" border-radius: 8px; color: {COLORS['text']}; padding: 8px 12px; }}")
        outer.addWidget(banner)

        body = QtWidgets.QHBoxLayout()
        body.setSpacing(14)
        outer.addLayout(body, 1)

        steps_host = QtWidgets.QWidget()
        steps_layout = QtWidgets.QVBoxLayout(steps_host)
        steps_layout.setContentsMargins(0, 0, 4, 0)
        steps_layout.setSpacing(8)
        actions = (
            [("Open First-Time Setup", self.step1_first_run),
             ("With a failed download", lambda: self.step1_first_run(fail_checkm2=True)),
             ("AMRFinderPlus DB (command line)", self.step1_amrf)],
            [("Open Settings", self.step2_settings)],
            [("Open Settings", self.step3_settings), ("List the downloads", self.step3_list)],
            [("Missing scheme", lambda: self.step4_run()),
             ("Terms not accepted", lambda: self.step4_run(forget_terms=True)),
             ("No scheme", lambda: self.step4_no_scheme()),
             ("Import its schema", self.step4_import),
             ("Start window", self.step4_start_window)],
            [("Show the commands", self.step5_cli)],
        )
        for number, ((heading, text), buttons) in enumerate(zip(_STEPS, actions), 1):
            steps_layout.addWidget(self._step_card(number, heading, text, buttons))
        steps_layout.addStretch(1)
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setWidget(steps_host)
        scroll.setMinimumWidth(540)
        body.addWidget(scroll, 5)

        right = QtWidgets.QVBoxLayout()
        right.setSpacing(4)
        caption = QtWidgets.QLabel("Console (what the install prints; the run log of step 4)")
        caption.setProperty("role", "settingsCaption")
        right.addWidget(caption)
        self.console = QtWidgets.QPlainTextEdit()
        self.console.setObjectName("DemoConsole")
        self.console.setReadOnly(True)
        self.console.setMaximumBlockCount(6000)
        self.console.setLineWrapMode(QtWidgets.QPlainTextEdit.LineWrapMode.WidgetWidth)
        mono = QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.SystemFont.FixedFont)
        mono.setPointSize(max(8, mono.pointSize() - 1))
        self.console.setFont(mono)
        right.addWidget(self.console, 1)
        body.addLayout(right, 6)

        footer = QtWidgets.QFrame()
        footer.setObjectName("StartFooter")
        row = QtWidgets.QHBoxLayout(footer)
        row.setContentsMargins(0, 10, 0, 0)
        self.status = QtWidgets.QLabel("")
        self.status.setProperty("role", "windowHint")
        row.addWidget(self.status, 1)
        folder_btn = QtWidgets.QPushButton("Open demo folder")
        folder_btn.clicked.connect(lambda: QtGui.QDesktopServices.openUrl(QtCore.QUrl.fromLocalFile(str(self.env.root))))
        row.addWidget(folder_btn)
        clear_btn = QtWidgets.QPushButton("Clear console")
        clear_btn.setProperty("role", "flat")
        clear_btn.clicked.connect(self.console.clear)
        row.addWidget(clear_btn)
        quit_btn = QtWidgets.QPushButton("Quit")
        quit_btn.setProperty("role", "primary")
        quit_btn.clicked.connect(self.close)
        row.addWidget(quit_btn)
        outer.addWidget(footer)
        self.setStyleSheet(
            f"QDialog#DemoLauncher {{ background: {COLORS['window']}; }}"
            f"QDialog#DemoLauncher QPushButton[role=\"primary\"]:disabled {{ background: {COLORS['disabled_fill']};"
            f" border-color: {COLORS['disabled_border']}; color: {COLORS['text_disabled']}; }}"
            f"QFrame#DemoStep {{ background: {COLORS['surface']}; border: 1px solid {COLORS['border']};"
            f" border-radius: 10px; }}"
            f"QLabel#DemoStepNumber {{ background: {COLORS['primary_soft']}; color: {COLORS['primary_ink']};"
            f" border-radius: 12px; font-weight: 700; min-width: 24px; max-width: 24px; min-height: 24px;"
            f" max-height: 24px; qproperty-alignment: AlignCenter; }}")

    def _step_card(self, number: int, heading: str, text: str, buttons) -> QtWidgets.QFrame:
        card = QtWidgets.QFrame()
        card.setObjectName("DemoStep")
        lay = QtWidgets.QVBoxLayout(card)
        lay.setContentsMargins(14, 10, 14, 10)
        lay.setSpacing(5)
        head = QtWidgets.QHBoxLayout()
        head.setSpacing(10)
        num = QtWidgets.QLabel(str(number))
        num.setObjectName("DemoStepNumber")
        head.addWidget(num, 0, QtCore.Qt.AlignmentFlag.AlignTop)
        label = QtWidgets.QLabel(heading)
        label.setProperty("role", "sectionTitle")
        head.addWidget(label, 1)
        lay.addLayout(head)
        hint = QtWidgets.QLabel(text)
        hint.setWordWrap(True)
        hint.setProperty("role", "windowHint")
        lay.addWidget(hint)
        # At most two buttons a row, so a step's buttons never widen the steps column.
        per_row = 2
        for start in range(0, len(buttons), per_row):
            row = QtWidgets.QHBoxLayout()
            row.setSpacing(6)
            for i, (caption, slot) in enumerate(buttons[start:start + per_row], start):
                btn = QtWidgets.QPushButton(caption)
                if i == 0:
                    btn.setProperty("role", "primary")
                btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
                btn.clicked.connect(slot)
                btn.setObjectName(f"DemoStep{number}Button{i + 1}")
                self._step_buttons.append(btn)
                row.addWidget(btn)
            row.addStretch(1)
            lay.addLayout(row)
        return card

    def _append(self, line: str) -> None:
        self.console.appendPlainText(line)
        bar = self.console.verticalScrollBar()
        bar.setValue(bar.maximum())

    def _set_busy(self, text: str) -> None:
        self._busy = bool(text)
        self.status.setText(text)
        for btn in self._step_buttons:
            if not btn.objectName().startswith("DemoStep4"):
                btn.setEnabled(not self._busy)

    def _background(self, what: str, job: Callable[[], Any], done: Callable[[Any], None]) -> None:
        """Run ``job`` off the GUI thread; ``done(result or exception)`` back on it."""
        self._set_busy(f"{what}...")

        def work() -> None:
            try:
                result = job()
            except BaseException as exc:  # noqa: BLE001 - handed to done()
                traceback.print_exc()
                result = exc
            self.bridge.finished.emit(done, result)

        threading.Thread(target=work, daemon=True, name=f"demo-{what}").start()

    # ---- step 1 ----
    def step1_first_run(self, *, fail_checkm2: bool = False) -> None:
        from annotation_explorer_polars.first_time_setup import FirstTimeSetupDialog

        if fail_checkm2:
            FAIL_NEXT_CHECKM2.set()
        started = time.time()
        say("")
        say("== First-Time Setup (output also under Details in the window) ==")
        # startup_dialog._ensure_app_config opens it the same way.
        dlg = FirstTimeSetupDialog(parent=self)
        accepted = dlg.exec() == QtWidgets.QDialog.DialogCode.Accepted
        FAIL_NEXT_CHECKM2.clear()
        say(f"First-Time Setup {'finished (Start)' if accepted else 'closed'}.")
        self._step1_records(started)

    def _step1_records(self, started: float) -> None:
        """What setup installed or kept, from the install records and the config (console)."""
        from pipeline.DbDownload import read_install_record

        config = self.env.config()
        lines = []
        for title, folder, used in (
            ("Bakta DB", self.env.root / "bakta_db", config.get("bakta_db")),
            ("CheckM2 DIAMOND DB", self.env.root / "checkm2_db",
             (config.get("assembly_qc") or {}).get("checkm2_database_path")),
        ):
            rec = read_install_record(folder)
            if not rec:
                lines.append(f"{title}: not installed{' (configured as ' + str(used) + ')' if used else ''}")
                continue
            fresh = rec.get("installed_at", "") >= datetime.fromtimestamp(started).astimezone().isoformat(timespec="seconds")
            version = rec.get("version") or {}
            if isinstance(version, dict) and version.get("major") is not None:
                vtext = f"v{version.get('major')}.{version.get('minor')} {version.get('type', '')} ({version.get('date')})"
            else:
                vtext = ""
            lines.append(f"{title}: {'installed now' if fresh else 'kept (installed before)'}"
                         f"{', ' + vtext if vtext else ''}, {human_size(rec.get('size_bytes'))}")
            amrf = rec.get("amrfinderplus_db") or {}
            if amrf:
                lines.append(f"    bundled AMRFinderPlus DB {amrf.get('release')} (format {amrf.get('format')})")
            if rec.get("sha256"):
                lines.append(f"    sha256 {rec['sha256'][:16]}...")
            lines.append(f"    {rec.get('path')}")
            lines.append(f"    record: {folder / 'bactolomew_install.json'}")
        lines.append(f"bactolomew_config.json: bakta_db = {config.get('bakta_db')}")
        lines.append("assembly_qc.checkm2_database_path = "
                     f"{(config.get('assembly_qc') or {}).get('checkm2_database_path')}")
        say("\n".join(["== Installed records =="] + lines))

    def step1_amrf(self) -> None:
        import pipeline.DbDownload as dbd

        argv = ["--install", "amrfinderplus"]
        say("")
        say(f"$ PYTHONPATH=src python -m pipeline.DbDownload {shlex.join(argv)}")
        self._background("AMRFinderPlus DB install", lambda: dbd.main(argv),
                         lambda res: (self._set_busy(""), say(f"(exit status {res})")))

    # ---- steps 2 and 3 ----
    def _settings(self) -> None:
        from annotation_explorer_polars.settings_dialog import SettingsDialog

        SettingsDialog(self).exec()

    def step2_settings(self) -> None:
        say("")
        say("== Settings: click \"Download from cgMLST.org...\" ==")
        self._settings()

    def step3_list(self) -> None:
        say("")
        say(f"== {self.env.downloads} ==")
        for line in self.env.describe_downloads():
            say(line)

    def step3_settings(self) -> None:
        self.step3_list()
        say("== Settings: click \"Add schema...\" (the source is filled in), then Import ==")
        self._settings()

    # ---- step 4 ----
    def _patch_start_window(self) -> None:
        """The start window's actions that would start the real pipeline act on the demo."""
        from annotation_explorer_polars import startup_dialog as sd

        launcher = self

        def new_run(dialog) -> None:  # noqa: ANN001
            QtWidgets.QMessageBox.information(dialog, "Start New Run",
                                              f"{DEMO_TAG}: a real run is not started here. Use step 4 of the demo.")

        def view(dialog, entry) -> None:  # noqa: ANN001
            if not entry.can_view:
                QtWidgets.QMessageBox.information(dialog, "Run not ready", entry.status_text or "Not ready.")
                return
            QtWidgets.QMessageBox.information(dialog, "View", f"{DEMO_TAG}: the explorer is not opened for a "
                                              "demo run (its outputs are placeholders).")

        def resume(dialog, entry) -> None:  # noqa: ANN001
            launcher.resume_run(entry.path, dialog)

        sd.StartupDialog._open_new_run_dialog = new_run
        sd.StartupDialog._open_run_from_entry = view
        sd.StartupDialog._resume_run = resume

    def _start_run(self, run_dir: Path, heading: str) -> None:
        """Start the simulated run of ``run_dir`` the way the start window starts one: the
        cgMLST.org terms first when its chewBBACA step may install a scheme."""
        from annotation_explorer_polars.startup_dialog import ask_cgmlst_terms_before_run

        startup = self.step4_start_window()
        ask_cgmlst_terms_before_run(startup, run_dir, None)
        self._run_count += 1
        label = str(self._run_count)
        say("")
        say(f"== Run {label}: {heading} ==")
        sim = SimulatedRun(self.env, run_dir, label=label)
        self._runs[str(run_dir)] = sim
        sim.start()
        startup._load_runs()
        startup._open_progress_dialog(run_dir)

    def step4_run(self, *, forget_terms: bool = False) -> Optional[Path]:
        from pipeline.CgmlstInstall import TERMS_FILE

        if forget_terms:
            (self.env.cgmlst / TERMS_FILE).unlink(missing_ok=True)
            say(f"Removed the cgMLST.org terms acceptance ({self.env.cgmlst / TERMS_FILE}).")
        taken = {run.state.get("demo_schema_id", "") for run in self._runs.values() if run.is_alive()}
        entry = pick_run_species(self.env, taken)
        if entry is None:
            QtWidgets.QMessageBox.information(self, "No species left",
                                              "Every catalog species with an MLST scheme is installed already.")
            return None
        species = (entry.get("species") or [entry["schema_id"]])[0]
        run_dir = create_run(self.env, species, entry["mlst_schemes"][0], tag="missing")
        state = json.loads((run_dir / "run_state.json").read_text(encoding="utf-8"))
        state["demo_schema_id"] = entry["schema_id"]
        (run_dir / "run_state.json").write_text(json.dumps(state, indent=2), encoding="utf-8")
        self._start_run(run_dir, f"12 {species} genomes (MLST scheme {entry['mlst_schemes'][0]}), no "
                                 f"'{entry['schema_id']}' schema installed")
        return run_dir

    def step4_no_scheme(self) -> Optional[Path]:
        species, mlst, _first = NO_SCHEME_SPECIES
        run_dir = create_run(self.env, species, mlst, tag=f"noscheme{self._run_count + 1}")
        self._start_run(run_dir, f"12 {species} genomes (MLST scheme {mlst}): no scheme on cgMLST.org, none "
                                 "installed")
        return run_dir

    def step4_import(self) -> None:
        """Settings > Add schema with the no-scheme species' files, its name and MLST scheme
        filled in (the user would type them)."""
        from annotation_explorer_polars.schema_import_dialog import AddSchemaDialog

        species, mlst, _first = NO_SCHEME_SPECIES
        source = write_no_scheme_source(self.env)
        say("")
        say(f"== Add schema: {source} as {species} (MLST scheme {mlst}); then Resume the stopped run ==")
        dlg = AddSchemaDialog(self)
        dlg.source_edit.setText(str(source))
        dlg.species_edit.setText(species)
        dlg.mlst_edit.setText(mlst)
        dlg.exec()

    def resume_run(self, run_dir: Path, parent: QtWidgets.QWidget) -> None:
        from annotation_explorer_polars.startup_dialog import ask_cgmlst_terms_before_run

        sim = self._runs.get(str(Path(run_dir)))
        if sim is None:
            QtWidgets.QMessageBox.information(parent, "Resume", f"{DEMO_TAG}: only the demo's runs can resume.")
            return
        if sim.is_alive():
            QtWidgets.QMessageBox.information(parent, "Resume", "This run is still running.")
            return
        ask_cgmlst_terms_before_run(parent, Path(run_dir), None)
        say(f"== Resuming run {sim.label} ==")
        again = SimulatedRun(self.env, run_dir, label=sim.label, resume=True)
        self._runs[str(Path(run_dir))] = again
        again.start()
        if self._startup is not None:
            self._startup._load_runs()
            self._startup._open_progress_dialog(Path(run_dir))

    def step4_start_window(self):
        from annotation_explorer_polars.startup_dialog import StartupDialog

        if self._startup is None or not self._startup.isVisible():
            self._startup = StartupDialog(self)
            self._startup.setModal(False)
            self._startup.show()
        self._startup.raise_()
        self._startup.activateWindow()
        return self._startup

    # ---- step 5 ----
    def step5_cli(self) -> None:
        import pipeline.CgmlstInstall as ci
        import pipeline.DbDownload as dbd

        commands = [
            (dbd, ["--list"]),
            (dbd, ["--install", "bakta-full", "checkm2", "amrfinderplus", "--dry-run"]),
            (dbd, ["--install", "amrfinderplus", "--amrf-release", AMRF_RELEASE, "--dry-run"]),
            (ci, ["--list", "--sizes"]),
            (ci, ["--install", "Enterococcus faecium", "Smarcescens", "Abaumannii", "--dry-run"]),
            (ci, ["--import", str(self.env.downloads), "--dry-run"]),
            (ci, ["--status"]),
        ]

        def job() -> None:
            for module, argv in commands:
                say("")
                say(f"$ PYTHONPATH=src python -m {module.__name__} {shlex.join(argv)}")
                try:
                    code = module.main(argv)
                except SystemExit as exc:
                    code = exc.code
                say(f"(exit status {code})")

        say("")
        say(f"== Command line, against the demo folder (bactolomew_config.json and DB/ in {self.env.root}) ==")
        self._background("Command line", job, lambda _res: self._set_busy(""))

    # ---- exit ----
    def reject(self) -> None:  # Escape closes the demo like Quit
        self.close()

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:  # noqa: N802
        QtWidgets.QApplication.instance().quit()
        super().closeEvent(event)

    def shutdown(self) -> None:
        _SHUTDOWN.set()
        for sim in self._runs.values():
            sim.stop_process()
        sys.stdout, sys.stderr = sys.__stdout__, sys.__stderr__


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #
def main(argv: Optional[list[str]] = None) -> int:
    global SPEED
    parser = argparse.ArgumentParser(description="BACTolomew install demo (nothing is downloaded).")
    parser.add_argument("--speed", type=float, default=1.0,
                        help="divide every simulated duration by this (default 1; 4 is brisk)")
    parser.add_argument("--keep", action="store_true", help="keep the demo folder on exit")
    parser.add_argument("--root", type=Path, default=None,
                        help="use this empty folder instead of a new temporary one")
    args = parser.parse_args(argv)
    SPEED = max(0.05, float(args.speed))

    if args.root is not None:
        root = args.root.expanduser().resolve()
        if root.exists() and any(root.iterdir()):
            parser.error(f"{root} is not empty")
        root.mkdir(parents=True, exist_ok=True)
    else:
        root = Path(tempfile.mkdtemp(prefix="bactolomew-install-demo-"))
    print(f"{DEMO_TAG}. Demo folder: {root}", flush=True)
    env = DemoEnv(root)
    try:
        install_backend(env)  # before anything reads a config or a project folder
        env.build()
        isolate_qsettings(env.settings_dir)
        return run_launcher(env)
    finally:
        if args.keep:
            print(f"Demo folder kept: {root}", flush=True)
        else:
            shutil.rmtree(root, ignore_errors=True)
            print(f"Demo folder removed: {root}", flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
