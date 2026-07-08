#!/usr/bin/env python3
"""
Download Acinetobacter baumannii GCA assemblies not already present locally.

Reads the GenBank assembly summary for A. baumannii, builds download URLs for
the *_genomic.fna.gz FASTA files, skips any already found in /Data/ST492/fasta,
and downloads the rest in parallel.

Outputs are saved to /Data/ST492/fasta.
"""

import concurrent.futures
import gzip
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import List, Tuple

FASTA_DIR = Path("/Data/ST492/fastas")
ASM_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Acinetobacter_baumannii/assembly_summary.txt"
MAX_WORKERS = 16


def fetch_assembly_summary(url: str) -> List[Tuple[str, str]]:
    """Fetch and parse assembly_summary, returning list of (accession, ftp_path)."""
    try:
        with urllib.request.urlopen(url) as resp:
            text = resp.read().decode("utf-8", errors="replace")
    except urllib.error.URLError as exc:
        raise SystemExit(f"Failed to fetch assembly summary: {exc}")

    accessions = []
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        # assembly_accession is column 1, ftp_path is column 20 (0-based 19)
        if len(parts) < 20:
            continue
        acc = parts[0]
        ftp = parts[19]
        if acc.startswith("GCA_") and ftp and ftp != "na":
            accessions.append((acc, ftp))
    return accessions


def strip_all_suffixes(path: Path) -> str:
    """Strip all suffixes (including .gz) to get the base name with version."""
    name = path.name
    for suf in path.suffixes:
        if name.endswith(suf):
            name = name[: -len(suf)]
    return name


def accession_root(acc: str) -> str:
    """Return accession without version, e.g., GCA_12345.1 -> GCA_12345."""
    return acc.split(".")[0]


def build_download_tasks(entries: List[Tuple[str, str]], existing: set) -> List[Tuple[str, str, str]]:
    tasks = []
    for acc, ftp in entries:
        if accession_root(acc) in existing:
            continue
        base = ftp.rstrip("/").split("/")[-1]
        fname = f"{base}_genomic.fna.gz"
        url = f"{ftp}/{fname}".replace("ftp://", "https://")
        tasks.append((acc, url, base))
    return tasks


def download_one(acc: str, url: str, base: str, dest_dir: Path) -> Tuple[str, bool, str]:
    dest = dest_dir / f"{base}_genomic.fna"
    try:
        with urllib.request.urlopen(url) as resp:
            compressed = resp.read()
        data = gzip.decompress(compressed)
        with open(dest, "wb") as out_fh:
            out_fh.write(data)
        return acc, True, ""
    except urllib.error.URLError as exc:
        return acc, False, str(exc)
    except OSError as exc:
        return acc, False, str(exc)


def main() -> None:
    entries = fetch_assembly_summary(ASM_SUMMARY_URL)
    if not entries:
        raise SystemExit("No entries found in assembly summary.")

    FASTA_DIR.mkdir(parents=True, exist_ok=True)

    existing_roots = set()
    for p in FASTA_DIR.glob("GCA_*"):
        if p.suffix in {".gz", ".fna", ".fa", ".fasta"} or "".join(p.suffixes[-2:]) in {".fa.gz", ".fna.gz", ".fasta.gz"}:
            base = strip_all_suffixes(p)
            existing_roots.add(accession_root(base))
    tasks = build_download_tasks(entries, existing_roots)

    if not tasks:
        print("No new GCA assemblies to download.")
        return

    total = len(tasks)
    print(f"Found {total} new assemblies to download to {FASTA_DIR}")
    import pdb; pdb.set_trace()

    successes = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        fut_to_task = {ex.submit(download_one, acc, url, base, FASTA_DIR): (acc, url) for acc, url, base in tasks}
        for idx, fut in enumerate(concurrent.futures.as_completed(fut_to_task), start=1):
            acc, url = fut_to_task[fut]
            acc, ok, err = fut.result()
            if ok:
                successes += 1
            else:
                print(f"Failed {acc}: {err}", file=sys.stderr)
            if idx % 10 == 0 or idx == total:
                print(f"[{idx}/{total}] downloads attempted", flush=True)

    print(f"Finished: {successes}/{total} downloaded successfully; outputs in {FASTA_DIR}")


if __name__ == "__main__":
    main()
