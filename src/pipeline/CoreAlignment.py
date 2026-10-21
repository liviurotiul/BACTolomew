"""Genome-ordered core gene alignment (step ``core_alignment``).

Panaroo concatenates the core gene alignments in its own cluster order, so
neighbouring columns of ``core_gene_alignment.aln`` come from unrelated genes and
Gubbins' spatial scan finds recombination blocks spanning them (review P39). This
step reorders the gene blocks along one real genome, the "ordering genome":

- automatic: the most contiguous cohort assembly (fewest contigs, then highest N50,
  then name) that places at least ``phylo.ordering_min_placed`` of the core genes;
- ``phylo.order_reference`` = a cohort sample name: that sample;
- ``phylo.order_reference`` = a genome FASTA path: the core genes' representatives
  are mapped onto it with minimap2.

Genes follow the ordering genome's contigs (longest first) and their start; a gene
on its minus strand is reverse-complemented so a tract crossing neighbouring genes
stays one column range. Core genes the ordering genome lacks are appended in
Panaroo order and counted. The same pass sanitizes the alignment (uppercase,
anything outside ACGTN- becomes N) and measures the cohort's diversity.

Outputs (``core_alignment/``): the ordered alignment (the single input of Gubbins,
the core tree and fastGEAR, review P63), ``gene_map.tsv`` (where each gene sits in
the ordered alignment and in the ordering genome; ``map_alignment_interval`` maps a
column range back to genes) and ``ordering.json`` (every choice, the inputs'
sha256 for reuse, the diversity).
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import logging
import os
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence

import numpy as np

from pipeline.AssemblyQualityWorkflow import _fasta_metrics
from pipeline.ToolEnv import tool_cmd

ALGORITHM_VERSION = 1
SCHEMA_VERSION = 1
CORE_ALIGNMENT_DIR = "core_alignment"
ORDERED_ALIGNMENT_NAME = "core_gene_alignment.ordered.aln"
GENE_MAP_NAME = "gene_map.tsv"
ORDERING_JSON_NAME = "ordering.json"
ORDERED_ALIGNMENT_REL = Path(CORE_ALIGNMENT_DIR) / ORDERED_ALIGNMENT_NAME
CORE_ALIGNMENT_OUT_FILES: tuple[Path, ...] = (
    Path(ORDERED_ALIGNMENT_NAME),
    Path(GENE_MAP_NAME),
    Path(ORDERING_JSON_NAME),
)

# Ids of one cluster on one contig and strand closer than this are fragments of one gene.
FRAGMENT_GAP_BP = 5000
MAX_ORDERING_CANDIDATES = 20
EXTERNAL_MIN_QUERY_COVERAGE = 0.8
WRAP = 60

_PANAROO_ALIGNMENTS: dict[str, tuple[str, str]] = {
    "full": ("core_gene_alignment", "core_alignment_header.embl"),
    "filtered": ("core_gene_alignment_filtered", "core_alignment_filtered_header.embl"),
}
_ALIGNMENT_SUFFIXES = (".aln", ".fa", ".fasta", ".aln.gz", ".fa.gz", ".fasta.gz")

GENE_MAP_COLUMNS = (
    "order",
    "cluster",
    "aln_start",
    "aln_end",
    "length",
    "placed",
    "reversed",
    "ordering_contig",
    "ordering_start",
    "ordering_end",
    "ordering_strand",
    "ordering_locus_tag",
    "copies_in_ordering_genome",
    "fragmented",
    "in_panaroo_filtered",
    "panaroo_aln_start",
    "panaroo_aln_end",
)


# --------------------------------------------------------------------- inputs
@dataclass(frozen=True)
class PanarooAlignment:
    mode: str
    alignment: Path
    header: Path
    filtered_header: Path | None


def resolve_panaroo_core_alignment(output_dir: Path, mode: str = "full") -> PanarooAlignment:
    """The Panaroo core alignment every core consumer uses (``phylo.core_alignment``)."""
    panaroo_out = Path(output_dir) / "panaroo" / "out"
    if mode not in _PANAROO_ALIGNMENTS:
        raise ValueError(f"phylo.core_alignment must be one of {', '.join(_PANAROO_ALIGNMENTS)}: {mode!r}")
    stem, header_name = _PANAROO_ALIGNMENTS[mode]
    alignment = next(
        (
            panaroo_out / f"{stem}{suffix}"
            for suffix in _ALIGNMENT_SUFFIXES
            if (panaroo_out / f"{stem}{suffix}").is_file() and (panaroo_out / f"{stem}{suffix}").stat().st_size > 0
        ),
        None,
    )
    if alignment is None:
        raise RuntimeError(
            f"Panaroo {mode} core alignment not found under {panaroo_out} "
            f"(looked for {stem} with {', '.join(_ALIGNMENT_SUFFIXES)})."
        )
    header = panaroo_out / header_name
    if not header.is_file():
        raise RuntimeError(
            f"Panaroo alignment header {header} is missing: the gene blocks of {alignment.name} "
            "cannot be located, so the alignment cannot be ordered."
        )
    filtered_header = panaroo_out / _PANAROO_ALIGNMENTS["filtered"][1]
    return PanarooAlignment(
        mode=mode,
        alignment=alignment,
        header=header,
        filtered_header=filtered_header if filtered_header.is_file() else None,
    )


@dataclass(frozen=True)
class GeneBlock:
    cluster: str
    start: int  # 1-based inclusive, in the Panaroo alignment
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


_EMBL_FEATURE = re.compile(r"^FT\s+feature\s+(\d+)\.\.(\d+)\s*$")
_EMBL_LABEL = re.compile(r'^FT\s+/label="?([^"\s]+)"?\s*$')


def parse_embl_header(path: Path) -> list[GeneBlock]:
    """Gene blocks (cluster, start, end) of a Panaroo core alignment header, in file order."""
    blocks: list[GeneBlock] = []
    pending: tuple[int, int] | None = None
    with Path(path).open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            found = _EMBL_FEATURE.match(line)
            if found:
                pending = (int(found.group(1)), int(found.group(2)))
                continue
            found = _EMBL_LABEL.match(line)
            if found and pending is not None:
                cluster = found.group(1)
                if cluster.endswith(".aln"):
                    cluster = cluster[: -len(".aln")]
                blocks.append(GeneBlock(cluster=cluster, start=pending[0], end=pending[1]))
                pending = None
    if not blocks:
        raise RuntimeError(f"No gene blocks found in Panaroo alignment header {path}.")
    previous_end = 0
    for block in blocks:
        if block.start != previous_end + 1 or block.end < block.start:
            raise RuntimeError(
                f"Panaroo alignment header {path} is not a contiguous block list "
                f"(block {block.cluster} at {block.start}..{block.end} after {previous_end})."
            )
        previous_end = block.end
    return blocks


def read_presence_cells(gpa_csv: Path, samples: Sequence[str]) -> dict[str, dict[str, str]]:
    """{sample: {cluster: cell}} from Panaroo's gene_presence_absence.csv."""
    out: dict[str, dict[str, str]] = {sample: {} for sample in samples}
    with Path(gpa_csv).open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        columns = {name: index for index, name in enumerate(header)}
        wanted = [(sample, columns[sample]) for sample in samples if sample in columns]
        for row in reader:
            if not row:
                continue
            cluster = row[0]
            for sample, index in wanted:
                if index < len(row) and row[index].strip():
                    out[sample][cluster] = row[index].strip()
    return out


def presence_samples(gpa_csv: Path) -> list[str]:
    with Path(gpa_csv).open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        header = next(csv.reader(handle))
    return [name for name in header[3:] if name]


@dataclass
class GffIndex:
    contig_lengths: dict[str, int]
    features: dict[str, tuple[str, int, int, str]]  # ID -> (contig, start, end, strand)

    def contig_rank(self) -> dict[str, int]:
        ordered = sorted(self.contig_lengths, key=lambda c: (-self.contig_lengths[c], c))
        return {contig: rank for rank, contig in enumerate(ordered)}


def read_gff_index(path: Path) -> GffIndex:
    lengths: dict[str, int] = {}
    features: dict[str, tuple[str, int, int, str]] = {}
    with Path(path).open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4:
                    lengths[parts[1]] = int(parts[3])
                continue
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            if cols[2] == "region":
                lengths.setdefault(cols[0], int(cols[4]))
                continue
            attrs = cols[8]
            found = re.search(r"(?:^|;)ID=([^;]+)", attrs)
            if not found:
                continue
            feature_id = found.group(1)
            if feature_id not in features:
                features[feature_id] = (cols[0], int(cols[3]), int(cols[4]), cols[6])
    return GffIndex(contig_lengths=lengths, features=features)


def bakta_gff_path(output_dir: Path, sample: str) -> Path | None:
    path = Path(output_dir) / "bakta" / sample / f"{sample}.gff3"
    return path if path.is_file() else None


# ------------------------------------------------------------------ placement
@dataclass(frozen=True)
class Placement:
    contig: str
    contig_rank: int
    start: int
    end: int
    strand: str
    locus_tag: str
    copies: int
    fragmented: bool


@dataclass
class PlacementCounts:
    placed: int = 0
    multi_copy: int = 0
    fragmented: int = 0
    refound_only: int = 0
    no_gff_entry: int = 0
    absent: int = 0


def place_from_cells(
    clusters: Sequence[str],
    cells: Mapping[str, str],
    gff: GffIndex,
) -> tuple[dict[str, Placement], PlacementCounts]:
    """Where each cluster sits in one cohort genome, from its presence cell and GFF3.

    Ids on one contig and strand within FRAGMENT_GAP_BP of each other are fragments
    of one gene (one placement spanning them). Separate groups make the cluster
    multi-copy: the copy with the lowest (contig rank, start) is used. Refound ids and
    ids without a GFF row do not place a cluster.
    """
    rank = gff.contig_rank()
    placements: dict[str, Placement] = {}
    counts = PlacementCounts()
    for cluster in clusters:
        cell = cells.get(cluster, "")
        ids = [x.strip() for x in cell.split(";") if x.strip()]
        if not ids:
            counts.absent += 1
            continue
        located: list[tuple[int, int, int, str, str, str]] = []
        refound = 0
        for feature_id in ids:
            if "refound" in feature_id:
                refound += 1
                continue
            hit = gff.features.get(feature_id)
            if hit is None:
                continue
            contig, start, end, strand = hit
            located.append((rank.get(contig, len(rank)), start, end, contig, strand, feature_id))
        if not located:
            if refound:
                counts.refound_only += 1
            else:
                counts.no_gff_entry += 1
            continue
        located.sort()
        groups: list[list[tuple[int, int, int, str, str, str]]] = []
        for item in located:
            if groups:
                last = groups[-1][-1]
                if item[3] == last[3] and item[4] == last[4] and item[1] - last[2] <= FRAGMENT_GAP_BP:
                    groups[-1].append(item)
                    continue
            groups.append([item])
        first = groups[0]
        placements[cluster] = Placement(
            contig=first[0][3],
            contig_rank=first[0][0],
            start=min(x[1] for x in first),
            end=max(x[2] for x in first),
            strand=first[0][4],
            locus_tag=";".join(x[5] for x in first),
            copies=len(groups),
            fragmented=len(first) > 1,
        )
        counts.placed += 1
        counts.multi_copy += int(len(groups) > 1)
        counts.fragmented += int(len(first) > 1)
    return placements, counts


def _cluster_representatives(panaroo_out: Path, clusters: Sequence[str], dest: Path) -> int:
    """Write one representative sequence per core cluster (Panaroo's pan genome
    reference, else the first sequence of the cluster's gene alignment)."""
    wanted = set(clusters)
    written: set[str] = set()
    with dest.open("w", encoding="utf-8") as out:
        reference = panaroo_out / "pan_genome_reference.fa"
        if reference.is_file():
            for name, seq in _iter_fasta_text(reference):
                if name in wanted and name not in written and seq:
                    out.write(f">{name}\n{seq}\n")
                    written.add(name)
        aligned_dir = panaroo_out / "aligned_gene_sequences"
        for cluster in clusters:
            if cluster in written:
                continue
            for suffix in (".aln.fas", ".aln", ".fas"):
                path = aligned_dir / f"{cluster}{suffix}"
                if path.is_file():
                    for _name, seq in _iter_fasta_text(path):
                        seq = seq.replace("-", "")
                        if seq:
                            out.write(f">{cluster}\n{seq}\n")
                            written.add(cluster)
                            break
                    break
    return len(written)


def _iter_fasta_text(path: Path) -> Iterator[tuple[str, str]]:
    name: str | None = None
    chunks: list[str] = []
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].strip().split()[0] if line[1:].strip() else ""
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def place_on_external_genome(
    panaroo_out: Path,
    clusters: Sequence[str],
    genome: Path,
    work_dir: Path,
    logger: logging.Logger,
) -> tuple[dict[str, Placement], PlacementCounts, dict[str, int]]:
    """Place core clusters on an external genome by mapping their representatives
    with minimap2 (best hit with query coverage >= 0.8; strand from the hit)."""
    work_dir.mkdir(parents=True, exist_ok=True)
    queries = work_dir / "core_representatives.fa"
    n_queries = _cluster_representatives(panaroo_out, clusters, queries)
    paf = work_dir / "core_representatives.paf"
    cmd = [
        *tool_cmd("minimap2", workdir=str(work_dir)),
        "-c",
        "-x",
        "asm20",
        "--secondary=no",
        "-o",
        str(paf),
        str(genome),
        str(queries),
    ]
    logger.info("Mapping %d core gene representatives onto %s with minimap2.", n_queries, genome)
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"minimap2 failed on {genome}: {(result.stderr or result.stdout).strip()[-2000:]}")
    best: dict[str, tuple[int, str, int, int, str]] = {}
    contig_lengths: dict[str, int] = {}
    with paf.open("r", encoding="utf-8") as handle:
        for line in handle:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch = (
                cols[0], int(cols[1]), int(cols[2]), int(cols[3]), cols[4], cols[5], int(cols[6]),
                int(cols[7]), int(cols[8]), int(cols[9]),
            )
            contig_lengths[tname] = tlen
            if qlen <= 0 or (qend - qstart) / qlen < EXTERNAL_MIN_QUERY_COVERAGE:
                continue
            if qname not in best or nmatch > best[qname][0]:
                best[qname] = (nmatch, tname, tstart + 1, tend, strand)
    for name, seq in _iter_fasta_text(genome):
        contig_lengths.setdefault(name, len(seq))
    rank = {c: i for i, c in enumerate(sorted(contig_lengths, key=lambda c: (-contig_lengths[c], c)))}
    placements = {
        cluster: Placement(
            contig=hit[1],
            contig_rank=rank.get(hit[1], len(rank)),
            start=hit[2],
            end=hit[3],
            strand=hit[4],
            locus_tag="",
            copies=1,
            fragmented=False,
        )
        for cluster, hit in best.items()
    }
    counts = PlacementCounts(placed=len(placements), absent=len(clusters) - len(placements))
    for path in (queries, paf):
        path.unlink(missing_ok=True)
    return placements, counts, contig_lengths


# --------------------------------------------------------------------- choice
@dataclass
class OrderingChoice:
    method: str  # auto_most_contiguous | config_sample | external_genome
    genome: str  # sample name or external path
    genome_path: Path | None
    gff_path: Path | None
    contigs: int | None
    n50: int | None
    placements: dict[str, Placement]
    counts: PlacementCounts
    candidates: list[dict[str, Any]] = field(default_factory=list)
    warnings: list[dict[str, str]] = field(default_factory=list)


def _contiguity(genomes: Sequence[Path], samples: set[str]) -> list[tuple[int, int, str, Path]]:
    ranked: list[tuple[int, int, str, Path]] = []
    for genome in genomes:
        stem = Path(genome).stem
        if stem not in samples:
            continue
        metrics = _fasta_metrics(Path(genome))
        contigs = metrics.get("contigs")
        n50 = metrics.get("n50")
        if contigs is None:
            continue
        ranked.append((int(contigs), -int(n50 or 0), stem, Path(genome)))
    ranked.sort()
    return ranked


def choose_ordering_genome(
    output_dir: Path,
    genomes: Sequence[Path],
    blocks: Sequence[GeneBlock],
    settings: Mapping[str, Any],
    logger: logging.Logger,
) -> OrderingChoice:
    output_dir = Path(output_dir)
    panaroo_out = output_dir / "panaroo" / "out"
    gpa = panaroo_out / "gene_presence_absence.csv"
    clusters = [block.cluster for block in blocks]
    total = len(clusters)
    override = str(settings.get("order_reference") or "").strip()
    stems = {Path(g).stem: Path(g) for g in genomes}
    table_samples = set(presence_samples(gpa))

    def from_sample(sample: str) -> tuple[dict[str, Placement], PlacementCounts, Path]:
        gff_path = bakta_gff_path(output_dir, sample)
        if gff_path is None:
            raise RuntimeError(f"Bakta GFF3 of the ordering genome {sample} not found under {output_dir / 'bakta'}.")
        cells = read_presence_cells(gpa, [sample])[sample]
        placements, counts = place_from_cells(clusters, cells, read_gff_index(gff_path))
        return placements, counts, gff_path

    if override and override in stems:
        placements, counts, gff_path = from_sample(override)
        metrics = _fasta_metrics(stems[override])
        choice = OrderingChoice(
            method="config_sample",
            genome=override,
            genome_path=stems[override],
            gff_path=gff_path,
            contigs=metrics.get("contigs"),
            n50=metrics.get("n50"),
            placements=placements,
            counts=counts,
        )
    elif override:
        external = Path(override).expanduser()
        if not external.is_file():
            raise RuntimeError(
                f"phylo.order_reference {override!r} is neither a cohort sample nor a genome FASTA file."
            )
        work = output_dir / CORE_ALIGNMENT_DIR / "external_mapping"
        placements, counts, contig_lengths = place_on_external_genome(panaroo_out, clusters, external, work, logger)
        try:
            work.rmdir()
        except OSError:
            pass
        lengths = sorted(contig_lengths.values(), reverse=True)
        half, acc, n50 = sum(lengths) / 2.0, 0, None
        for length in lengths:
            acc += length
            if acc >= half:
                n50 = length
                break
        choice = OrderingChoice(
            method="external_genome",
            genome=str(external.resolve()),
            genome_path=external.resolve(),
            gff_path=None,
            contigs=len(lengths),
            n50=n50,
            placements=placements,
            counts=counts,
        )
    else:
        min_placed = float(settings.get("ordering_min_placed", 0.95))
        ranked = [
            item for item in _contiguity(genomes, table_samples)
            if bakta_gff_path(output_dir, item[2]) is not None
        ][:MAX_ORDERING_CANDIDATES]
        if not ranked:
            raise RuntimeError("No cohort genome with a Bakta GFF3 can order the core alignment.")
        cells_by_sample = read_presence_cells(gpa, [item[2] for item in ranked])
        evaluated: list[tuple[tuple[int, int, str, Path], dict[str, Placement], PlacementCounts, Path]] = []
        winner = None
        for item in ranked:
            gff_path = bakta_gff_path(output_dir, item[2])
            assert gff_path is not None
            placements, counts = place_from_cells(clusters, cells_by_sample[item[2]], read_gff_index(gff_path))
            evaluated.append((item, placements, counts, gff_path))
            if counts.placed >= min_placed * total:
                winner = evaluated[-1]
                break
        warnings: list[dict[str, str]] = []
        if winner is None:
            winner = max(evaluated, key=lambda e: (e[2].placed, -e[0][0], -e[0][1]))
            message = (
                f"No candidate among the {len(evaluated)} most contiguous genomes places "
                f"{min_placed:.0%} of the {total} core genes; ordering by {winner[0][2]}, which "
                f"places {winner[2].placed}. The rest are appended in Panaroo order."
            )
            logger.warning(message)
            warnings.append({"step": "core_alignment", "code": "ordering_many_unplaced", "message": message})
        item, placements, counts, gff_path = winner
        choice = OrderingChoice(
            method="auto_most_contiguous",
            genome=item[2],
            genome_path=item[3],
            gff_path=gff_path,
            contigs=item[0],
            n50=-item[1],
            placements=placements,
            counts=counts,
            candidates=[
                {"sample": e[0][2], "contigs": e[0][0], "n50": -e[0][1], "placed": e[2].placed, "total": total}
                for e in evaluated
            ],
            warnings=warnings,
        )
    return choice


# ------------------------------------------------------------------ gene map
@dataclass(frozen=True)
class GeneMapRow:
    order: int
    cluster: str
    aln_start: int
    aln_end: int
    placed: bool
    reversed: bool
    ordering_contig: str
    ordering_start: int | None
    ordering_end: int | None
    ordering_strand: str
    ordering_locus_tag: str
    copies_in_ordering_genome: int
    fragmented: bool
    in_panaroo_filtered: bool | None
    panaroo_aln_start: int
    panaroo_aln_end: int

    @property
    def length(self) -> int:
        return self.aln_end - self.aln_start + 1


def build_gene_map(
    blocks: Sequence[GeneBlock],
    placements: Mapping[str, Placement],
    filtered_clusters: set[str] | None,
) -> list[GeneMapRow]:
    """Ordered gene blocks: placed genes by (contig rank, start), then the unplaced in Panaroo order."""
    indexed = list(enumerate(blocks))
    placed = sorted(
        (item for item in indexed if item[1].cluster in placements),
        key=lambda item: (
            placements[item[1].cluster].contig_rank,
            placements[item[1].cluster].start,
            placements[item[1].cluster].end,
            item[0],
        ),
    )
    unplaced = [item for item in indexed if item[1].cluster not in placements]
    rows: list[GeneMapRow] = []
    position = 1
    for order, (_index, block) in enumerate([*placed, *unplaced], start=1):
        placement = placements.get(block.cluster)
        rows.append(
            GeneMapRow(
                order=order,
                cluster=block.cluster,
                aln_start=position,
                aln_end=position + block.length - 1,
                placed=placement is not None,
                reversed=bool(placement and placement.strand == "-"),
                ordering_contig=placement.contig if placement else "",
                ordering_start=placement.start if placement else None,
                ordering_end=placement.end if placement else None,
                ordering_strand=placement.strand if placement else "",
                ordering_locus_tag=placement.locus_tag if placement else "",
                copies_in_ordering_genome=placement.copies if placement else 0,
                fragmented=bool(placement and placement.fragmented),
                in_panaroo_filtered=None if filtered_clusters is None else block.cluster in filtered_clusters,
                panaroo_aln_start=block.start,
                panaroo_aln_end=block.end,
            )
        )
        position += block.length
    return rows


def write_gene_map(rows: Sequence[GeneMapRow], path: Path) -> None:
    def cell(value: Any) -> str:
        if value is None:
            return ""
        if isinstance(value, bool):
            return "true" if value else "false"
        return str(value)

    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(GENE_MAP_COLUMNS) + "\n")
        for row in rows:
            values = [
                row.order, row.cluster, row.aln_start, row.aln_end, row.length, row.placed, row.reversed,
                row.ordering_contig, row.ordering_start, row.ordering_end, row.ordering_strand,
                row.ordering_locus_tag, row.copies_in_ordering_genome, row.fragmented,
                row.in_panaroo_filtered, row.panaroo_aln_start, row.panaroo_aln_end,
            ]
            handle.write("\t".join(cell(v) for v in values) + "\n")


def read_gene_map(path: Path) -> list[GeneMapRow]:
    def opt_int(text: str) -> int | None:
        return int(text) if text else None

    def opt_bool(text: str) -> bool | None:
        return None if text == "" else text == "true"

    rows: list[GeneMapRow] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for rec in reader:
            rows.append(
                GeneMapRow(
                    order=int(rec["order"]),
                    cluster=rec["cluster"],
                    aln_start=int(rec["aln_start"]),
                    aln_end=int(rec["aln_end"]),
                    placed=rec["placed"] == "true",
                    reversed=rec["reversed"] == "true",
                    ordering_contig=rec["ordering_contig"],
                    ordering_start=opt_int(rec["ordering_start"]),
                    ordering_end=opt_int(rec["ordering_end"]),
                    ordering_strand=rec["ordering_strand"],
                    ordering_locus_tag=rec["ordering_locus_tag"],
                    copies_in_ordering_genome=int(rec["copies_in_ordering_genome"] or 0),
                    fragmented=rec["fragmented"] == "true",
                    in_panaroo_filtered=opt_bool(rec["in_panaroo_filtered"]),
                    panaroo_aln_start=int(rec["panaroo_aln_start"]),
                    panaroo_aln_end=int(rec["panaroo_aln_end"]),
                )
            )
    return rows


def map_alignment_interval(rows: Sequence[GeneMapRow], start: int, end: int) -> list[dict[str, Any]]:
    """Genes a column range (1-based inclusive, ordered alignment) overlaps.

    For each gene: the overlapped columns, the same range in gene coordinates (columns
    of the Panaroo gene block, counted in its coding direction, so inverted for a
    reversed gene) and the approximate ordering-genome coordinates (start of the gene
    plus the column offset; alignment gaps are not subtracted).
    """
    if end < start:
        start, end = end, start
    out: list[dict[str, Any]] = []
    for row in rows:
        if row.aln_end < start or row.aln_start > end:
            continue
        lo, hi = max(start, row.aln_start), min(end, row.aln_end)
        if row.reversed:
            gene_lo, gene_hi = row.aln_end - hi + 1, row.aln_end - lo + 1
        else:
            gene_lo, gene_hi = lo - row.aln_start + 1, hi - row.aln_start + 1
        genome_lo = genome_hi = None
        if row.placed and row.ordering_start is not None and row.ordering_end is not None:
            genome_lo = min(row.ordering_end, row.ordering_start + (lo - row.aln_start))
            genome_hi = min(row.ordering_end, row.ordering_start + (hi - row.aln_start))
        out.append(
            {
                "cluster": row.cluster,
                "order": row.order,
                "aln_start": lo,
                "aln_end": hi,
                "gene_start": gene_lo,
                "gene_end": gene_hi,
                "reversed": row.reversed,
                "placed": row.placed,
                "ordering_contig": row.ordering_contig,
                "ordering_start": genome_lo,
                "ordering_end": genome_hi,
                "in_panaroo_filtered": row.in_panaroo_filtered,
            }
        )
    return out


# ------------------------------------------------------------------ alignment
def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return Path(path).open("r", encoding="utf-8", errors="ignore")


def _open_binary(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rb")
    return Path(path).open("rb")


def iter_fasta_records(path: Path) -> Iterator[tuple[str, bytes]]:
    """(name, sequence bytes with line breaks removed) per record, streamed."""
    name: str | None = None
    chunks: list[bytes] = []
    with _open_binary(path) as handle:
        for line in handle:
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(chunks)
                name = line[1:].strip().decode("utf-8", errors="replace")
                if not name:
                    raise RuntimeError(f"Empty FASTA header in {path}.")
                chunks = []
            else:
                stripped = line.strip()
                if stripped:
                    if name is None:
                        raise RuntimeError(f"FASTA sequence data before the first header in {path}.")
                    chunks.append(stripped)
    if name is not None:
        yield name, b"".join(chunks)


_SANITIZE = np.full(256, ord("N"), dtype=np.uint8)
for _ch in b"ACGTN-":
    _SANITIZE[_ch] = _ch
    _SANITIZE[ord(chr(_ch).lower())] = _ch
_VALID = np.zeros(256, dtype=bool)
for _ch in b"ACGTNacgtn-":
    _VALID[_ch] = True
_COMPLEMENT = np.arange(256, dtype=np.uint8)
for _a, _b in (("A", "T"), ("C", "G"), ("G", "C"), ("T", "A")):
    _COMPLEMENT[ord(_a)] = ord(_b)
# Base codes: A C G T = 0..3, anything else (N, -) = 4.
BASE_CODE = np.full(256, 4, dtype=np.uint8)
for _i, _ch in enumerate(b"ACGT"):
    BASE_CODE[_ch] = _i


def _wrap(seq: np.ndarray) -> bytes:
    n = seq.size
    full = n // WRAP
    parts: list[bytes] = []
    if full:
        body = np.empty((full, WRAP + 1), dtype=np.uint8)
        body[:, :WRAP] = seq[: full * WRAP].reshape(full, WRAP)
        body[:, WRAP] = ord("\n")
        parts.append(body.tobytes())
    if n % WRAP:
        parts.append(seq[full * WRAP:].tobytes() + b"\n")
    return b"".join(parts)


def fasta_record(name: str, seq: np.ndarray) -> bytes:
    """One FASTA record in the ordered alignment's format (60 columns per line)."""
    return f">{name}\n".encode("utf-8") + _wrap(seq)


@dataclass
class OrderedAlignmentStats:
    sequence_count: int
    alignment_length: int
    sanitized_characters: int
    sha256: str
    missing_fraction: dict[str, float]
    diversity: dict[str, Any]


def write_ordered_alignment(
    source: Path,
    rows: Sequence[GeneMapRow],
    dest: Path,
    *,
    expected_length: int,
) -> OrderedAlignmentStats:
    """Stream the source alignment into dest in gene-map order (reversed genes
    reverse-complemented), sanitized, and count the bases per column."""
    perm = np.concatenate(
        [
            np.arange(row.panaroo_aln_end - 1, row.panaroo_aln_start - 2, -1, dtype=np.int64)
            if row.reversed
            else np.arange(row.panaroo_aln_start - 1, row.panaroo_aln_end, dtype=np.int64)
            for row in rows
        ]
    )
    if perm.size != expected_length:
        raise RuntimeError(f"Gene blocks cover {perm.size} columns, the alignment has {expected_length}.")
    rev_mask = np.zeros(perm.size, dtype=bool)
    for row in rows:
        if row.reversed:
            rev_mask[row.aln_start - 1: row.aln_end] = True
    rev_idx = np.flatnonzero(rev_mask)

    digest = hashlib.sha256()
    counts = np.zeros((4, perm.size), dtype=np.int32)
    missing: dict[str, float] = {}
    sanitized = 0
    seen: set[str] = set()
    tmp = dest.with_name(dest.name + ".tmp")
    with tmp.open("wb") as out:
        for name, raw in iter_fasta_records(source):
            if name in seen:
                raise RuntimeError(f"Duplicate sequence name {name!r} in {source}.")
            seen.add(name)
            arr = np.frombuffer(raw, dtype=np.uint8)
            if arr.size != expected_length:
                raise RuntimeError(
                    f"Sequence {name} in {source} has {arr.size} columns, the header says {expected_length}."
                )
            sanitized += int(np.count_nonzero(~_VALID[arr]))
            seq = _SANITIZE[arr][perm]
            if rev_idx.size:
                seq[rev_idx] = _COMPLEMENT[seq[rev_idx]]
            codes = BASE_CODE[seq]
            for k in range(4):
                counts[k] += codes == k
            missing[name] = float(np.count_nonzero(codes == 4)) / float(seq.size)
            chunk = fasta_record(name, seq)
            digest.update(chunk)
            out.write(chunk)
    os.replace(tmp, dest)
    if len(seen) < 1:
        raise RuntimeError(f"No sequences in {source}.")
    diversity = _diversity(dest, counts)
    return OrderedAlignmentStats(
        sequence_count=len(seen),
        alignment_length=int(perm.size),
        sanitized_characters=sanitized,
        sha256=digest.hexdigest(),
        missing_fraction=missing,
        diversity=diversity,
    )


def _diversity(path: Path, counts: np.ndarray) -> dict[str, Any]:
    """Mean and max pairwise p-distance over the columns with >= 2 called bases.

    Mean: average over those columns of 1 - sum c_a (c_a - 1) / (m (m - 1)). Max: the
    largest pairwise count of differences at variable columns (both called), divided by
    the number of columns with >= 2 called bases.
    """
    called = counts.sum(axis=0).astype(np.float64)
    usable = called >= 2
    n_usable = int(np.count_nonzero(usable))
    if n_usable == 0:
        return {"columns_compared": 0, "variable_sites": 0, "mean_pairwise_p": 0.0, "max_pairwise_p": 0.0}
    c = counts[:, usable].astype(np.float64)
    m = called[usable]
    same = (c * (c - 1.0)).sum(axis=0) / (m * (m - 1.0))
    mean_p = float(np.mean(1.0 - same))
    variable_idx = np.flatnonzero((counts > 0).sum(axis=0) >= 2)
    names: list[str] = []
    rows: list[np.ndarray] = []
    for name, raw in iter_fasta_records(path):
        codes = BASE_CODE[np.frombuffer(raw, dtype=np.uint8)]
        rows.append(codes[variable_idx])
        names.append(name)
    max_p, max_pair, median_snps = 0.0, None, 0.0
    if variable_idx.size and len(rows) >= 2:
        matrix = np.vstack(rows)
        n = matrix.shape[0]
        both = np.zeros((n, n), dtype=np.float64)
        same_pairs = np.zeros((n, n), dtype=np.float64)
        step = 20000
        for lo in range(0, matrix.shape[1], step):
            block = matrix[:, lo: lo + step]
            called_block = (block < 4).astype(np.float32)
            both += called_block @ called_block.T
            for k in range(4):
                onehot = (block == k).astype(np.float32)
                same_pairs += onehot @ onehot.T
        diffs = both - same_pairs
        iu = np.triu_indices(n, k=1)
        pair_diffs = diffs[iu]
        best = int(np.argmax(pair_diffs))
        max_p = float(pair_diffs[best]) / float(n_usable)
        max_pair = [names[iu[0][best]], names[iu[1][best]]]
        median_snps = float(np.median(pair_diffs))
    return {
        "columns_compared": n_usable,
        "variable_sites": int(variable_idx.size),
        "mean_pairwise_p": mean_p,
        "max_pairwise_p": max_p,
        "max_pair": max_pair,
        "median_pairwise_snps": median_snps,
    }


# ------------------------------------------------------------------ run step
def file_sha256(path: Path, *, chunk: int = 1 << 22) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def _genome_set_fingerprint(genomes: Iterable[Any]) -> str:
    joined = "\n".join(sorted({Path(str(g)).stem for g in genomes}))
    return hashlib.sha256(joined.encode("utf-8")).hexdigest()


def _settings_snapshot(settings: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "core_alignment": settings.get("core_alignment"),
        "order_reference": str(settings.get("order_reference") or ""),
        "ordering_min_placed": float(settings.get("ordering_min_placed", 0.95)),
    }


def _run_relative(output_dir: Path, path: Path) -> str:
    """A path inside the run folder as written relative to it (so a moved run folder
    keeps its stamps); a path outside it stays absolute."""
    try:
        return Path(path).resolve().relative_to(Path(output_dir).resolve()).as_posix()
    except ValueError:
        return str(path)


def _recorded_path(output_dir: Path, text: Any) -> Path:
    path = Path(str(text or ""))
    return path if path.is_absolute() else Path(output_dir) / path


def _input_fingerprints(
    output_dir: Path, paths: PanarooAlignment, gpa: Path, choice_files: Mapping[str, Path | None]
) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, path in (
        ("source_alignment", paths.alignment),
        ("header", paths.header),
        ("filtered_header", paths.filtered_header),
        ("gene_presence_absence", gpa),
        *choice_files.items(),
    ):
        if path is None:
            out[key] = None
            continue
        out[key] = {"path": _run_relative(output_dir, path), "sha256": file_sha256(path)}
    return out


@dataclass(frozen=True)
class CoreAlignmentResult:
    output_dir: Path
    alignment: Path
    gene_map: Path
    ordering_json: Path
    summary: dict[str, Any]
    warnings: list[dict[str, str]]


def predicted_gubbins_dropouts(missing_fraction: Mapping[str, float], filter_percentage: float) -> list[str]:
    """Taxa Gubbins will filter out (more than filter_percentage % of N and gap)."""
    return sorted(name for name, frac in missing_fraction.items() if frac * 100.0 > float(filter_percentage))


def run_core_alignment(
    output_dir: Path,
    genomes: Sequence[Path],
    settings: Mapping[str, Any],
    logger: logging.Logger,
) -> CoreAlignmentResult:
    output_dir = Path(output_dir).expanduser().resolve()
    out_dir = output_dir / CORE_ALIGNMENT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    panaroo = resolve_panaroo_core_alignment(output_dir, str(settings.get("core_alignment") or "full"))
    gpa = output_dir / "panaroo" / "out" / "gene_presence_absence.csv"
    blocks = parse_embl_header(panaroo.header)
    filtered_clusters = (
        {block.cluster for block in parse_embl_header(panaroo.filtered_header)}
        if panaroo.filtered_header is not None
        else None
    )
    logger.info(
        "Ordering the Panaroo %s core alignment (%d genes, %d columns).",
        panaroo.mode,
        len(blocks),
        blocks[-1].end,
    )
    choice = choose_ordering_genome(output_dir, genomes, blocks, settings, logger)
    rows = build_gene_map(blocks, choice.placements, filtered_clusters)
    alignment = out_dir / ORDERED_ALIGNMENT_NAME
    stats = write_ordered_alignment(panaroo.alignment, rows, alignment, expected_length=blocks[-1].end)
    gene_map = out_dir / GENE_MAP_NAME
    write_gene_map(rows, gene_map)

    filter_percentage = float(settings.get("gubbins_filter_percentage", 25.0))
    dropouts = predicted_gubbins_dropouts(stats.missing_fraction, filter_percentage)
    if dropouts:
        logger.warning(
            "%d sample(s) have more than %.0f%% N or gap in the core alignment and will be left out "
            "by Gubbins: %s",
            len(dropouts),
            filter_percentage,
            ", ".join(dropouts),
        )
    placed = sum(1 for row in rows if row.placed)
    summary: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "algorithm_version": ALGORITHM_VERSION,
        "method": choice.method,
        "ordering_genome": choice.genome,
        "ordering_genome_path": str(choice.genome_path) if choice.genome_path else None,
        "ordering_genome_contigs": choice.contigs,
        "ordering_genome_n50": choice.n50,
        "candidates": choice.candidates,
        "settings": _settings_snapshot(settings),
        "genome_set_fingerprint": _genome_set_fingerprint(genomes),
        "source_alignment_mode": panaroo.mode,
        "genes_total": len(rows),
        "genes_placed": placed,
        "genes_appended": len(rows) - placed,
        "genes_reversed": sum(1 for row in rows if row.reversed),
        "genes_multi_copy": choice.counts.multi_copy,
        "genes_fragmented": choice.counts.fragmented,
        "genes_refound_only": choice.counts.refound_only,
        "genes_without_gff_entry": choice.counts.no_gff_entry,
        "genes_absent_from_ordering_genome": choice.counts.absent,
        "genes_in_panaroo_filtered": None if filtered_clusters is None else sum(
            1 for row in rows if row.in_panaroo_filtered
        ),
        "sequence_count": stats.sequence_count,
        "alignment_length": stats.alignment_length,
        "sanitized_characters": stats.sanitized_characters,
        "diversity": stats.diversity,
        "gubbins_filter_percentage": filter_percentage,
        "predicted_gubbins_dropouts": dropouts,
        "missing_fraction_max": max(stats.missing_fraction.values()) if stats.missing_fraction else 0.0,
        # Every sample's N and gap fraction: Gubbins recomputes its dropouts from it with
        # the current phylo.gubbins_filter_percentage (not part of this step's settings).
        "missing_fraction": {name: round(frac, 6) for name, frac in sorted(stats.missing_fraction.items())},
        "inputs": _input_fingerprints(
            output_dir,
            panaroo,
            gpa,
            {
                "ordering_genome_gff": choice.gff_path,
                "external_genome": choice.genome_path if choice.method == "external_genome" else None,
            },
        ),
        "output": {"path": ORDERED_ALIGNMENT_REL.as_posix(), "sha256": stats.sha256},
        "gene_map": _run_relative(output_dir, gene_map),
        "warnings": choice.warnings,
        "notes": [
            "Genes on the ordering genome's minus strand are reverse-complemented (A/T and C/G "
            "swapped, N and gaps kept).",
            "Intergenic regions are absent: Gubbins window sizes and block lengths are in "
            "concatenated core-gene columns, not genome base pairs.",
            "Core genes the ordering genome lacks are appended at the end in Panaroo order.",
        ],
    }
    ordering_json = out_dir / ORDERING_JSON_NAME
    tmp = ordering_json.with_name(ordering_json.name + ".tmp")
    tmp.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    os.replace(tmp, ordering_json)
    logger.info(
        "Core alignment ordered by %s (%s): %d of %d genes placed, %d appended, %d reversed; "
        "mean pairwise p-distance %.4f, max %.4f, %d variable sites.",
        choice.genome,
        choice.method,
        placed,
        len(rows),
        len(rows) - placed,
        summary["genes_reversed"],
        stats.diversity["mean_pairwise_p"],
        stats.diversity["max_pairwise_p"],
        stats.diversity["variable_sites"],
    )
    return CoreAlignmentResult(
        output_dir=out_dir,
        alignment=alignment,
        gene_map=gene_map,
        ordering_json=ordering_json,
        summary=summary,
        warnings=list(choice.warnings),
    )


def read_ordering_summary(output_dir: Path) -> dict[str, Any] | None:
    try:
        data = json.loads((Path(output_dir) / CORE_ALIGNMENT_DIR / ORDERING_JSON_NAME).read_text(encoding="utf-8"))
    except Exception:
        return None
    return data if isinstance(data, dict) else None


def core_alignment_reuse_issue(
    output_dir: Path,
    genomes: Sequence[Path],
    settings: Mapping[str, Any],
) -> str | None:
    """Why the existing ordered alignment cannot be reused (None when it can)."""
    output_dir = Path(output_dir)
    summary = read_ordering_summary(output_dir)
    if summary is None:
        return "no ordering record (core_alignment/ordering.json)"
    if summary.get("algorithm_version") != ALGORITHM_VERSION:
        return "the ordering algorithm changed"
    if summary.get("settings") != _settings_snapshot(settings):
        return "phylo ordering settings changed"
    if summary.get("genome_set_fingerprint") != _genome_set_fingerprint(genomes):
        return "the genome set changed"
    try:
        panaroo = resolve_panaroo_core_alignment(output_dir, str(settings.get("core_alignment") or "full"))
    except Exception as exc:
        return str(exc)
    inputs = summary.get("inputs") or {}
    for key, path in (
        ("source_alignment", panaroo.alignment),
        ("header", panaroo.header),
        ("filtered_header", panaroo.filtered_header),
        ("gene_presence_absence", output_dir / "panaroo" / "out" / "gene_presence_absence.csv"),
    ):
        recorded = inputs.get(key)
        if path is None and recorded is None:
            continue
        if not isinstance(recorded, Mapping) or path is None or not Path(path).is_file():
            return f"input {key} changed"
        if recorded.get("sha256") != file_sha256(path):
            return f"input {key} changed"
    for key in ("ordering_genome_gff", "external_genome"):
        recorded = inputs.get(key)
        if isinstance(recorded, Mapping):
            path = _recorded_path(output_dir, recorded.get("path"))
            if not path.is_file() or file_sha256(path) != recorded.get("sha256"):
                return f"input {key} changed"
    alignment = output_dir / ORDERED_ALIGNMENT_REL
    output = summary.get("output") or {}
    if not alignment.is_file() or output.get("sha256") != file_sha256(alignment):
        return "the ordered alignment is missing or was modified"
    if not (output_dir / CORE_ALIGNMENT_DIR / GENE_MAP_NAME).is_file():
        return "gene_map.tsv is missing"
    return None


def ordered_alignment_sha256(output_dir: Path) -> str | None:
    summary = read_ordering_summary(output_dir)
    output = (summary or {}).get("output") or {}
    value = output.get("sha256")
    return str(value) if value else None


__all__ = [
    "ALGORITHM_VERSION",
    "BASE_CODE",
    "CORE_ALIGNMENT_DIR",
    "CORE_ALIGNMENT_OUT_FILES",
    "CoreAlignmentResult",
    "GENE_MAP_NAME",
    "GeneBlock",
    "GeneMapRow",
    "ORDERED_ALIGNMENT_NAME",
    "ORDERED_ALIGNMENT_REL",
    "ORDERING_JSON_NAME",
    "build_gene_map",
    "core_alignment_reuse_issue",
    "fasta_record",
    "file_sha256",
    "iter_fasta_records",
    "map_alignment_interval",
    "ordered_alignment_sha256",
    "parse_embl_header",
    "place_from_cells",
    "predicted_gubbins_dropouts",
    "read_gene_map",
    "read_gff_index",
    "read_ordering_summary",
    "resolve_panaroo_core_alignment",
    "run_core_alignment",
    "write_ordered_alignment",
]
