#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd


DEFAULT_CSV = "/home/bals/Documents/test-Contexar3/annotations.csv"


def load_blaoxa72(csv_path: Path) -> pd.DataFrame:
    if not csv_path.is_file():
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    df = pd.read_csv(csv_path)
    required_cols = {"Gene", "ConsensusCategory", "Markers"}
    missing = sorted(required_cols.difference(df.columns))
    if missing:
        raise ValueError(f"Missing required columns in {csv_path}: {', '.join(missing)}")

    gene_mask = (
        df["Gene"]
        .astype(str)
        .str.strip()
        .str.casefold()
        .eq("blaoxa-72")
    )

    cols = ["Gene", "ConsensusCategory", "Markers"]
    if "Sample" in df.columns:
        cols = ["Sample"] + cols

    return df.loc[gene_mask, cols].reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "List all entries where Gene == blaOXA-72 and show their "
            "ConsensusCategory and Markers."
        )
    )
    parser.add_argument(
        "csv",
        nargs="?",
        default=DEFAULT_CSV,
        help=f"Path to annotations CSV (default: {DEFAULT_CSV})",
    )
    args = parser.parse_args()

    csv_path = Path(args.csv).expanduser()

    try:
        hits = load_blaoxa72(csv_path)
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(str(exc))

    if hits.empty:
        print("No blaOXA-72 entries found.")
        return

    with pd.option_context("display.max_colwidth", None):
        print(hits.to_string(index=False))


if __name__ == "__main__":
    main()
