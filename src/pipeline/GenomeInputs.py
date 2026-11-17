"""Normalized genome inputs of a run.

Assemblers write FASTA headers the tools downstream do not all read the same way
(Shovill/SPAdes: ">NODE_1_length_5000_cov_10 len=5000 cov=10 sw=..."). A new run copies
each genome into its folder as ``inputs/<sample>.fna`` (gzipped inputs are decompressed)
with every header cut to its first whitespace-separated token, characters outside
[A-Za-z0-9_.-] replaced by "_" and made unique within the genome; the sequence lines are
copied unchanged (line endings become "\\n", blank lines are dropped). Each genome's
original-to-new contig names are in ``inputs/<sample>.contig_map.tsv``; run_config.json
``inputs.genome_inputs`` records the originals with their sha256, and ``inputs.genomes``
keeps the paths as given (a resume rebuilds the run from them). Every tool of the run
reads the copies.

A run folder made before this (it records a genome set but no ``genome_inputs``) keeps
reading its original files, so a resume finds its steps' inputs unchanged.
"""
from __future__ import annotations

import copy
import gzip
import re
from pathlib import Path
from typing import Any, Mapping

from pipeline.CoreAlignment import file_sha256

INPUTS_DIRNAME = "inputs"
GENOME_INPUTS_KEY = "genome_inputs"
GENOME_INPUTS_VERSION = 1
CONTIG_MAP_HEADER = "original_header\toriginal_id\tcontig\n"
_UNSAFE_CONTIG_CHARS = re.compile(r"[^A-Za-z0-9_.-]")
_GZIP_MAGIC = b"\x1f\x8b"


def _open_text(path: Path):
    with Path(path).open("rb") as fh:
        gz = fh.read(2) == _GZIP_MAGIC
    if gz:
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return Path(path).open("r", encoding="utf-8", errors="replace", newline="")


def safe_contig_name(header: str) -> str:
    """A header's first whitespace-separated token with unsafe characters replaced."""
    token = header.split(maxsplit=1)[0] if header.strip() else ""
    return _UNSAFE_CONTIG_CHARS.sub("_", token) or "contig"


def normalize_genome_fasta(source: Path, dest: Path, contig_map: Path) -> dict[str, int]:
    """Write the normalized copy of one genome FASTA and its contig name map.

    Returns {"contigs": n, "renamed": contigs whose ID changed (a dropped description does not count)}. A file without any
    header is an error (it is not a FASTA)."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    used: set[str] = set()
    rows: list[str] = []
    renamed = 0
    tmp = dest.with_name(dest.name + ".tmp")
    with _open_text(source) as fin, tmp.open("w", encoding="utf-8", newline="\n") as fout:
        for raw in fin:
            line = raw.rstrip("\r\n")
            if line.startswith(">"):
                header = line[1:].strip()
                base = safe_contig_name(header)
                name, n = base, 1
                while name in used:
                    n += 1
                    name = f"{base}_{n}"
                used.add(name)
                original_id = header.split(maxsplit=1)[0] if header else ""
                renamed += int(name != original_id)
                rows.append(f"{header.replace(chr(9), ' ')}\t{original_id}\t{name}\n")
                fout.write(f">{name}\n")
            elif line.strip():
                if not used:
                    raise ValueError(f"{source}: sequence before the first FASTA header")
                fout.write(line + "\n")
    if not used:
        tmp.unlink(missing_ok=True)
        raise ValueError(f"{source}: no FASTA header found")
    map_tmp = contig_map.with_name(contig_map.name + ".tmp")
    map_tmp.write_text(CONTIG_MAP_HEADER + "".join(rows), encoding="utf-8")
    tmp.replace(dest)
    map_tmp.replace(contig_map)
    return {"contigs": len(used), "renamed": renamed}


def read_contig_map(path: Path) -> dict[str, str]:
    """original_id -> normalized contig name of one genome's contig_map.tsv."""
    out: dict[str, str] = {}
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) == 3:
            out.setdefault(parts[1], parts[2])
    return out


def _recorded(run_cfg: Mapping[str, Any] | None) -> dict | None:
    inputs = run_cfg.get("inputs") if isinstance(run_cfg, Mapping) else None
    record = inputs.get(GENOME_INPUTS_KEY) if isinstance(inputs, Mapping) else None
    return record if isinstance(record, dict) else None


def uses_normalized_inputs(run_cfg: Mapping[str, Any] | None, *, folder_has_run: bool) -> bool:
    """True for a new run folder and for one whose run_config records normalized inputs;
    False for a folder made before them (its tools keep reading the original files)."""
    return _recorded(run_cfg) is not None or not folder_has_run


def prepare_genome_inputs(cfg: Any, run_cfg: Mapping[str, Any] | None, logger=None) -> Any:
    """A copy of the run's DataInput reading the normalized genomes, with the record
    (``genome_inputs``) to write to run_config. A copy whose original still has the
    recorded sha256 is reused as it is (not rewritten, so its mtime stays)."""
    output = Path(cfg.output)
    folder = output / INPUTS_DIRNAME
    previous = (_recorded(run_cfg) or {}).get("samples") or {}
    samples: dict[str, dict] = {}
    by_original: dict[Path, Path] = {}
    for genome in cfg.genomes:
        original = Path(genome)
        sample = original.stem
        dest = folder / f"{sample}.fna"
        contig_map = folder / f"{sample}.contig_map.tsv"
        sha = file_sha256(original)
        before = previous.get(sample) if isinstance(previous.get(sample), Mapping) else {}
        if before.get("original_sha256") == sha and dest.is_file() and contig_map.is_file():
            counts = {"contigs": before.get("contigs"), "renamed": before.get("renamed")}
        else:
            counts = normalize_genome_fasta(original, dest, contig_map)
            if logger is not None and counts["renamed"]:
                logger.info(
                    f"Genome {sample}: {counts['renamed']} of {counts['contigs']} contig names "
                    f"normalized (map in {contig_map})"
                )
        samples[sample] = {
            "original": str(original),
            "original_sha256": sha,
            "fasta": str(dest.relative_to(output)),
            "contig_map": str(contig_map.relative_to(output)),
            **counts,
        }
        by_original[original] = dest
    out = copy.copy(cfg)
    out.genomes = [by_original[Path(g)] for g in cfg.genomes]
    if cfg.snippy_reference is not None:
        out.snippy_reference = by_original.get(Path(cfg.snippy_reference), cfg.snippy_reference)
    out.genome_inputs = {
        "version": GENOME_INPUTS_VERSION,
        "dir": INPUTS_DIRNAME,
        "given_genomes": [str(g) for g in cfg.genomes],
        "given_snippy_reference": str(cfg.snippy_reference) if cfg.snippy_reference else None,
        "samples": samples,
    }
    return out
