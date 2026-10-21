#!/usr/bin/env python3
"""Assemble the raw per-tool prediction outputs into /Data/ST492/supplementary as
a set of numbered, descriptively named folders suitable for deposition alongside
the manuscript.

The pipeline run directory holds ~25 GB, most of which is intermediates that a
reader neither needs nor can use: DIAMOND database indexes, per-sample Snippy
working directories, the Bakta .json/.embl/.gbff serialisations of data already
present in the GFF3, the un-filtered whole-genome alignments and the tool logs.
This script copies only the result files, compresses the per-genome collections,
and writes a manifest with SHA-256 checksums.

Re-running is safe: each destination folder is rebuilt from scratch.
"""
import hashlib
import json
import shutil
import subprocess
import tarfile
from pathlib import Path

RUN = Path("/home/bals/Documents/dev/BACTolomew/runs/sample_78eaa9126850/run_78eaa9126850")
OUT = Path("/Data/ST492/supplementary")

# Bakta writes 13 files per genome. Keep the annotation itself (GFF3), the
# feature table (TSV), the proteins (FAA), the per-genome summary (TXT) and the
# record of which database each annotation came from (inference.tsv). Drop
# .embl/.gbff/.json (the same content re-serialised), .ffn/.fna (nucleotide
# sequence, already public), .log and .hypotheticals.*.
BAKTA_KEEP = {
    ".gff3": ("bakta_annotations_gff3", "GFF3 genome annotation"),
    ".tsv": ("bakta_feature_tables", "tab-separated feature table"),
    ".faa": ("bakta_protein_sequences", "predicted protein sequences"),
    ".txt": ("bakta_genome_summaries", "per-genome annotation summary"),
    ".inference.tsv": ("bakta_annotation_inference", "evidence and source database per feature"),
}


def fresh(name: str) -> Path:
    d = OUT / name
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    return d


def copy(src: Path, dst_dir: Path, rename: str | None = None) -> Path | None:
    if not src.exists():
        print(f"    MISSING {src}")
        return None
    dst = dst_dir / (rename or src.name)
    shutil.copy2(src, dst)
    return dst


def gzip_copy(src: Path, dst_dir: Path, rename: str | None = None) -> Path | None:
    """Copy then gzip, for large single files that are still worth shipping."""
    dst = copy(src, dst_dir, rename)
    if dst is None:
        return None
    subprocess.run(["gzip", "-f", "-6", str(dst)], check=True)
    return dst.with_suffix(dst.suffix + ".gz")


def tar_gz(files: list[Path], dst: Path, arcdir: str) -> Path:
    """Bundle many small per-genome files into one archive."""
    with tarfile.open(dst, "w:gz", compresslevel=6) as tf:
        for f in sorted(files):
            tf.add(f, arcname=f"{arcdir}/{f.name}")
    return dst


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


built: list[tuple[str, str]] = []   # (relative path, description)


def note(p: Path, desc: str) -> None:
    built.append((str(p.relative_to(OUT)), desc))


# ---------------------------------------------------------------- 01  Bakta
print("01  genome annotation (Bakta)")
d = fresh("01_genome_annotation_bakta")
samples = sorted(x for x in (RUN / "bakta").iterdir() if x.is_dir())
print(f"    {len(samples)} genome directories")
for suffix, (stem, desc) in BAKTA_KEEP.items():
    files = [s / f"{s.name}{suffix}" for s in samples]
    files = [f for f in files if f.exists()]
    arc = tar_gz(files, d / f"{stem}.tar.gz", stem)
    note(arc, f"{desc}, one file per genome ({len(files)} genomes)")
    print(f"    {arc.name:<40} {len(files):>4} files  {arc.stat().st_size/1e6:>8.1f} MB")

# ------------------------------------------------- 02  DIAMOND marker screening
print("02  marker-gene screening (DIAMOND)")
d = fresh("02_marker_screening_diamond")
q = sorted((RUN / "annotations" / "bakta_vs_queries").glob("*.tsv"))
arc = tar_gz(q, d / "diamond_bakta_cds_vs_marker_databases.tar.gz",
             "diamond_bakta_cds_vs_marker_databases")
note(arc, f"DIAMOND blastx hits, Bakta CDS against each marker database ({len(q)} databases, all genomes pooled per file)")
print(f"    {arc.name:<52} {len(q):>4} files  {arc.stat().st_size/1e6:>8.1f} MB")

r = sorted((RUN / "annotations" / "reference_loci").glob("*.tsv"))
arc = tar_gz(r, d / "diamond_vs_paidb_pathogenicity_islands.tar.gz",
             "diamond_vs_paidb_pathogenicity_islands")
note(arc, f"DIAMOND hits against the PAIDB pathogenicity-island reference loci, one file per genome ({len(r)} genomes)")
print(f"    {arc.name:<52} {len(r):>4} files  {arc.stat().st_size/1e6:>8.1f} MB")

# -------------------------------------------------------------------- 03  MLST
print("03  MLST (Pasteur scheme)")
d = fresh("03_mlst_pasteur")
note(copy(RUN / "mlst" / "results.csv", d, "mlst_results.csv"), "Pasteur-scheme MLST allele profile and sequence type per genome")
note(copy(RUN / "mlst" / "results.json", d, "mlst_results.json"), "the same MLST calls with per-locus allele detail")

# ------------------------------------------------------------ 04  assembly QC
print("04  assembly quality")
d = fresh("04_assembly_quality")
note(copy(RUN / "assembly_qc" / "summary.tsv", d, "assembly_qc_summary.tsv"), "merged assembly QC summary per genome")
note(copy(RUN / "assembly_qc" / "quast" / "report.tsv", d, "quast_report.tsv"), "QUAST assembly statistics")
note(copy(RUN / "assembly_qc" / "quast" / "transposed_report.tsv", d, "quast_report_transposed.tsv"), "QUAST statistics, one row per genome")
note(copy(RUN / "assembly_qc" / "checkm2" / "quality_report.tsv", d, "checkm2_quality_report.tsv"), "CheckM2 completeness and contamination estimates")

# ---------------------------------------------------------------- 05  pangenome
print("05  pangenome (Panaroo)")
d = fresh("05_pangenome_panaroo")
p = RUN / "panaroo" / "out"
note(copy(p / "gene_presence_absence.csv", d), "Panaroo gene presence/absence matrix with cluster annotation")
note(copy(p / "gene_presence_absence.Rtab", d), "binary presence/absence matrix")
note(copy(p / "gene_presence_absence_roary.csv", d), "the same matrix in Roary column format")
note(copy(p / "struct_presence_absence.Rtab", d), "structural-variant presence/absence")
note(copy(p / "summary_statistics.txt", d), "core, soft-core, shell and cloud gene counts")
note(copy(p / "alignment_entropy.csv", d), "per-cluster alignment entropy")
note(copy(p / "pan_genome_reference.fa", d), "one representative nucleotide sequence per gene cluster")

# ------------------------------------------------------------------ 06  cgMLST
print("06  cgMLST (chewBBACA)")
d = fresh("06_cgmlst_chewbbaca")
c = RUN / "chewbbaca" / "out"
note(copy(c / "results_alleles.tsv", d), "chewBBACA allele call matrix, genomes by loci")
note(copy(c / "results_statistics.tsv", d), "per-genome allele call classification counts")
note(copy(c / "loci_summary_stats.tsv", d), "per-locus summary statistics")
note(gzip_copy(c / "results_contigsInfo.tsv", d), "contig and coordinates of every allele call")
note(copy(c / "paralogous_counts.tsv", d), "loci flagged as paralogous")

# ------------------------------------------------------------- 07  recombination
print("07  recombination removal (Gubbins)")
d = fresh("07_recombination_gubbins")
g = RUN / "gubbins"
note(copy(g / "core_gene.final_tree.tre", d, "gubbins_final_tree.tre"), "recombination-corrected maximum-likelihood tree")
note(copy(g / "core_gene.node_labelled.final_tree.tre", d, "gubbins_final_tree_node_labelled.tre"), "the same tree with internal nodes labelled")
note(copy(g / "core_gene.per_branch_statistics.csv", d, "gubbins_per_branch_statistics.csv"), "SNPs, recombination blocks and r/m per branch")
note(gzip_copy(g / "core_gene.recombination_predictions.gff", d, "gubbins_recombination_predictions.gff"), "coordinates of every predicted recombination block")
note(gzip_copy(g / "core_gene.summary_of_snp_distribution.vcf", d, "gubbins_snp_distribution.vcf"), "SNP distribution across the filtered alignment")
note(gzip_copy(g / "core_gene.filtered_polymorphic_sites.fasta", d, "gubbins_filtered_polymorphic_sites.fasta"), "recombination-filtered polymorphic-site alignment used for the final tree")
note(copy(g / "summary.json", d, "gubbins_run_summary.json"), "Gubbins run parameters and per-iteration summary")

# ------------------------------------------------------- 08  lineages (fastGEAR)
print("08  lineage and recombination assignment (fastGEAR)")
d = fresh("08_lineages_fastgear")
f = RUN / "fastgear"
note(copy(f / "sample_clades.tsv", d, "fastgear_sample_lineage_assignment.tsv"), "lineage (cluster) assigned to each genome")
note(copy(f / "output" / "lineage_information.txt", d, "fastgear_lineage_information.txt"), "lineage definitions and membership")
note(copy(f / "output" / "recombinations_recent.txt", d, "fastgear_recent_recombinations.txt"), "recent recombination events, with donor lineage and coordinates")
note(copy(f / "output" / "recombinations_ancestral.txt", d, "fastgear_ancestral_recombinations.txt"), "ancestral recombination events")

# ----------------------------------------------------- 09  reference mapping
print("09  reference mapping and variant calling (Snippy)")
d = fresh("09_reference_mapping_snippy")
s = RUN / "snippy"
note(copy(s / "core" / "core.aln", d, "snippy_core_snp_alignment.aln"), "core SNP alignment across all genomes")
note(copy(s / "core" / "core.txt", d, "snippy_core_summary.txt"), "per-genome mapping and variant summary")
note(gzip_copy(s / "core" / "core.tab", d, "snippy_core_variants.tab"), "every core variant site with the allele in each genome")
note(gzip_copy(s / "core" / "core.vcf", d, "snippy_core_variants.vcf"), "the same variants in VCF format")
note(copy(s / "reference.txt", d, "snippy_reference.txt"), "reference genome used for mapping")

# ------------------------------------------------- 10  plasmid reconstruction
print("10  plasmid reconstruction (MOB-suite)")
d = fresh("10_plasmid_reconstruction_mobsuite")
m = RUN / "mobsuite"
note(gzip_copy(m / "contigs.tsv", d, "mobsuite_contig_assignments.tsv"), "every contig assigned to chromosome or a plasmid cluster")
note(copy(m / "plasmids.tsv", d, "mobsuite_plasmid_clusters.tsv"), "reconstructed plasmids with replicon and relaxase typing")
per = sorted((m / "per_sample").glob("*/mobtyper_results.txt"))
if per:
    with tarfile.open(d / "mobtyper_results_per_genome.tar.gz", "w:gz", compresslevel=6) as tf:
        for pth in per:
            tf.add(pth, arcname=f"mobtyper_results_per_genome/{pth.parent.name}.mobtyper_results.txt")
    note(d / "mobtyper_results_per_genome.tar.gz", f"MOB-typer replicon, relaxase and mobility call per reconstructed plasmid ({len(per)} genomes)")
rep = sorted((m / "per_sample").glob("*/contig_report.txt"))
if rep:
    with tarfile.open(d / "contig_reports_per_genome.tar.gz", "w:gz", compresslevel=6) as tf:
        for pth in rep:
            tf.add(pth, arcname=f"contig_reports_per_genome/{pth.parent.name}.contig_report.txt")
    note(d / "contig_reports_per_genome.tar.gz", f"MOB-recon per-contig molecule assignment ({len(rep)} genomes)")

# --------------------------------------------------------------- 11  trees
print("11  phylogenetic trees")
d = fresh("11_phylogenetic_trees")
note(copy(RUN / "trees" / "core_gene_tree.nwk", d, "core_gene_tree.nwk"), "maximum-likelihood tree from the Panaroo core-gene alignment")
note(copy(RUN / "trees" / "snippy_raxml.tree", d, "snippy_core_snp_tree.nwk"), "maximum-likelihood tree from the Snippy core SNP alignment")
note(copy(RUN / "gubbins" / "core_gene.final_tree.tre", d, "gubbins_recombination_filtered_tree.nwk"), "recombination-filtered tree (the tree used for the figures)")

# ------------------------------------------------------------------ manifest
print("\nwriting manifest and checksums")
rows = []
for rel, desc in built:
    if rel is None:
        continue
    p = OUT / rel
    if not p.exists():
        continue
    rows.append({"file": rel, "bytes": p.stat().st_size,
                 "sha256": sha256(p), "description": desc})

(OUT / "MANIFEST.json").write_text(json.dumps(rows, indent=2) + "\n")
(OUT / "CHECKSUMS.sha256").write_text(
    "".join(f"{r['sha256']}  {r['file']}\n" for r in rows))

total = sum(r["bytes"] for r in rows)
print(f"  {len(rows)} files, {total/1e6:.0f} MB total")
print(f"  wrote {OUT}/MANIFEST.json and {OUT}/CHECKSUMS.sha256")
