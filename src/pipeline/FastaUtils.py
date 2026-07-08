from __future__ import annotations

import gzip
import logging
import os
import re
import time
from pathlib import Path

import pandas as pd

from pipeline.Utils import CATEGORIES

FASTA_EXTS = {".fasta", ".fa", ".fna"}

CONTIG_MAP: dict[str, dict[str, int]] = {}
FASTA_CACHE: dict[str, dict[str, str]] = {}

_RC = str.maketrans(
    "ACGTUMRWSYKVHDBNacgtumrwsykvhdbn",
    "TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn",
)


def prepare_query_fasta(
    genome: Path,
    tmp_root: Path,
    logger,
    write_map: bool = False,
) -> Path:
    """
    Create a sanitized FASTA copy of `genome` under `tmp_root` with the SAME filename,
    but with deflines replaced by a compact per-sample mapping: contig_<int>.
    """
    try:
        tmp_root.mkdir(parents=True, exist_ok=True)
        dest = tmp_root / genome.name
        sample = genome.stem

        try:
            if dest.exists() and dest.stat().st_mtime >= genome.stat().st_mtime:
                return dest
        except Exception:
            pass

        is_gz = genome.suffix.lower() == ".gz"
        in_open = (
            lambda p: gzip.open(p, "rt", encoding="utf-8", errors="ignore")
        ) if is_gz else (lambda p: open(p, "r", encoding="utf-8", errors="ignore"))
        out_open = (
            lambda p: gzip.open(p, "wt", encoding="utf-8")
        ) if is_gz else (lambda p: open(p, "w", encoding="utf-8"))

        map_fp = (tmp_root / (genome.name + ".headers.map.tsv")) if write_map else None
        map_fh = None
        if map_fp:
            map_fh = open(map_fp, "w", encoding="utf-8")
            map_fh.write("sample\told\tnew\tid\n")

        smap = CONTIG_MAP.setdefault(sample, {})

        with in_open(genome) as fin, out_open(dest) as fout:
            for line in fin:
                if not line.startswith(">"):
                    fout.write(line)
                    continue

                raw = line[1:].rstrip("\r\n")
                norm = re.sub(r"\s+", "_", raw)

                cid = smap.get(norm)
                if cid is None:
                    cid = len(smap) + 1
                    smap[norm] = cid

                new_header = f"contig_{cid}"
                fout.write(">" + new_header + "\n")

                if map_fh:
                    map_fh.write(f"{sample}\t{norm}\t{new_header}\t{cid}\n")

        if map_fh:
            map_fh.close()

        try:
            os.utime(dest, (time.time(), time.time()))
        except Exception:
            pass

        return dest

    except Exception as e:
        logger.error(f"Failed to prepare sanitized FASTA for {genome}: {e}")
        return genome


def is_fasta(p: Path) -> bool:
    return p.suffix.lower() in FASTA_EXTS


def collect_fastas(root: Path) -> list[Path]:
    return [p for ext in FASTA_EXTS for p in root.rglob(f"*{ext}")]


def extract_category(p: Path, default=None) -> str:
    m = p.name.split("_")[1].split(".")[0]
    return m


def normalize_cat(cat: str) -> str:
    resolved = cat if cat in CATEGORIES else None
    logging.getLogger("diamond_prep").debug("normalize_cat %s -> %s", cat, resolved)
    return resolved


def _load_fasta_dict(path: Path) -> dict[str, str]:
    key = str(path)
    if key in FASTA_CACHE:
        return FASTA_CACHE[key]
    d: dict[str, str] = {}
    if not path.exists():
        FASTA_CACHE[key] = d
        return d
    name = None
    buf = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    d[name] = "".join(buf)
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name is not None:
            d[name] = "".join(buf)
    FASTA_CACHE[key] = d
    return d


def subseq(fasta_path: Path, contig: str, start: int, stop: int, strand: str) -> str | None:
    if contig is None or pd.isna(start) or pd.isna(stop):
        return None
    seqs = _load_fasta_dict(fasta_path)
    s = seqs.get(str(contig))
    if not s:
        return None
    a = int(start) - 1
    b = int(stop)
    if a < 0 or b > len(s) or a >= b:
        return None
    frag = s[a:b]
    if str(strand).strip() == "-":
        frag = frag.translate(_RC)[::-1]
    return frag


__all__ = [
    "FASTA_EXTS",
    "CONTIG_MAP",
    "FASTA_CACHE",
    "prepare_query_fasta",
    "is_fasta",
    "collect_fastas",
    "extract_category",
    "normalize_cat",
    "subseq",
]
