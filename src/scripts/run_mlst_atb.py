#!/usr/bin/env python3
"""
Run mlst on FASTA files under /Data/ST492/fastas and collect results.
Outputs to /Data/ST492/metadata/subtype.csv (columns: sample, subtype) and
copies subtype 492 FASTAs to /Data/ST492/final_fasta.
Uses up to 100 workers (or CPU count if smaller) to run mlst in parallel.
"""

import csv
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Tuple

FASTA_DIR = Path("/Data/ST492/fastas")
OUT_CSV = Path("/Data/ST492/metadata/subtype.csv")
OUT_FASTA_DIR = Path("/Data/ST492/final_fasta")
MAX_WORKERS = 100


def find_fasta_files(directory: Path) -> List[Path]:
    exts = {".fa", ".fasta", ".fna"}
    files = []
    for path in directory.rglob("*"):
        if path.is_file() and path.suffix.lower() in exts:
            files.append(path)
    return sorted(files)


def run_mlst(path: Path) -> Tuple[str, str]:
    """Run mlst on a single file. Returns (sample, subtype)."""
    sample = path.stem
    cmd = ["mlst", str(path)]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as exc:
        return sample, f"ERROR:{exc.returncode}"

    line = res.stdout.strip().splitlines()
    if not line:
        return sample, "ERROR:empty_output"

    fields = line[0].split()
    # mlst output: FILE SCHEME ST [loci...]
    subtype = fields[2] if len(fields) >= 3 else "ERROR:parse"
    return sample, subtype


def main() -> None:
    if not FASTA_DIR.exists():
        raise SystemExit(f"Input directory does not exist: {FASTA_DIR}")

    files = find_fasta_files(FASTA_DIR)
    if not files:
        print(f"No FASTA files found under {FASTA_DIR}")
        return

    workers = min(MAX_WORKERS, os.cpu_count() or 1, len(files))
    total = len(files)
    print(f"Found {total} FASTA files under {FASTA_DIR}")
    print(f"Using {workers} parallel workers (max {MAX_WORKERS})")
    results: List[Tuple[str, str]] = []
    sample_to_path: Dict[str, Path] = {p.stem: p for p in files}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        future_map = {executor.submit(run_mlst, f): f for f in files}
        for idx, fut in enumerate(as_completed(future_map), start=1):
            results.append(fut.result())
            if idx % 10 == 0 or idx == total:
                print(f"[{idx}/{total}] completed", flush=True)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample", "subtype"])
        for sample, subtype in sorted(results):
            writer.writerow([sample, subtype])

    print(f"Wrote {len(results)} records to {OUT_CSV}")

    # Copy subtype 492 FASTAs into OUT_FASTA_DIR
    OUT_FASTA_DIR.mkdir(parents=True, exist_ok=True)
    copied = 0
    for sample, subtype in results:
        if subtype == "492":
            src = sample_to_path.get(sample)
            if src and src.exists():
                dest = OUT_FASTA_DIR / src.name
                shutil.copy2(src, dest)
                copied += 1
    print(f"Copied {copied} subtype 492 FASTA files to {OUT_FASTA_DIR}")


if __name__ == "__main__":
    main()
