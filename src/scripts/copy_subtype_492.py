#!/usr/bin/env python3
"""
Copy subtype 492 FASTA files listed in mlst_results.csv into /Data/ST492/final_fasta.
Assumes the CSV has columns: sample, subtype (as produced by run_mlst_atb.py).
"""

import csv
import shutil
from pathlib import Path

ATB_DIR = Path("/Data/ST492/ATB")
MLST_CSV = ATB_DIR / "mlst_results.csv"
OUT_DIR = Path("/Data/ST492/final_fasta")


def load_mlst(path: Path):
    if not path.exists():
        raise SystemExit(f"CSV not found: {path}")
    with path.open() as fh:
        reader = csv.DictReader(fh)
        if not {"sample", "subtype"} <= set(reader.fieldnames or []):
            raise SystemExit("CSV must contain 'sample' and 'subtype' columns.")
        return [(row["sample"], row["subtype"]) for row in reader]


def main() -> None:
    rows = load_mlst(MLST_CSV)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    copied = 0
    missing = 0
    for sample, subtype in rows:
        if subtype != "492":
            continue
        src = next((p for p in ATB_DIR.glob(f"{sample}.*") if p.suffix.lower() in {".fa", ".fasta", ".fna"}), None)
        if src and src.exists():
            shutil.copy2(src, OUT_DIR / src.name)
            copied += 1
        else:
            missing += 1

    print(f"Copied {copied} subtype 492 FASTA files to {OUT_DIR}")
    if missing:
        print(f"Skipped {missing} subtype 492 entries because FASTA files were not found.")


if __name__ == "__main__":
    main()
