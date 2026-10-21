#!/usr/bin/env python3
"""Build the ST492 pan-genome presence panel from the Panaroo output.

Reads Panaroo's gene_presence_absence.Rtab (orthologous-cluster x genome binary
matrix, computed over the full ST2/ST492 set), restricts to the 93 ST492 genomes
(the tree tips), keeps clusters present in >=1 ST492 genome, and classifies each
cluster by within-ST492 frequency using Panaroo's own bands:
  Core (>=99%)  Soft-core (95-<99%)  Shell (15-<95%)  Cloud (<15%)
Genes are ordered by descending frequency (rank). Writes:
  * st492_pangenome_genes.csv   (gene, freq, class, rank)
  * st492_pangenome_matrix.csv  (gene + one 0/1 column per ST492 genome)
for plot_st492_tree_pangenome_phandango.R (full-pangenome phandango, no core panel).
"""
import pandas as pd

RTAB = "/home/bals/Documents/dev/BACTolomew/runs/sample_78eaa9126850/run_78eaa9126850/panaroo/out/gene_presence_absence.Rtab"
METADATA = "/Data/ST492/metadata/metadata_ALL_new_merged_with_clusters_with_complex.csv"
GENES_CSV = "/Data/ST492/figures2026/st492_pangenome_genes.csv"
MATRIX_CSV = "/Data/ST492/figures2026/st492_pangenome_matrix.csv"

# --- ST492 samples ---
md = pd.read_csv(METADATA)
st492 = sorted(md.loc[md["subtype"].astype(str).str.strip().isin(["492", "492.0"]), "Sample"])

# --- Panaroo matrix, restricted to ST492 genomes ---
rt = pd.read_csv(RTAB, sep="\t", index_col=0)
cols = [c for c in rt.columns if c in set(st492)]
sub = rt[cols].astype(int)
n = len(cols)
freq = sub.sum(axis=1)
sub = sub.loc[freq >= 1]           # present in >=1 ST492 genome
freq = freq.loc[sub.index]


def pg_class(f: int) -> str:
    r = f / n
    if r >= 0.99:
        return "Core"
    if r >= 0.95:
        return "Soft-core"
    if r >= 0.15:
        return "Shell"
    return "Cloud"


genes = pd.DataFrame({"gene": sub.index, "freq": freq.values})
genes["class"] = genes["freq"].apply(pg_class)
genes = genes.sort_values(["freq", "gene"], ascending=[False, True]).reset_index(drop=True)
genes["rank"] = genes.index + 1
genes.to_csv(GENES_CSV, index=False)

# subset presence matrix in the same gene order
out = sub.loc[genes["gene"]].reset_index()
out = out.rename(columns={out.columns[0]: "gene"})
out.to_csv(MATRIX_CSV, index=False)

# --- report ---
counts = genes["class"].value_counts()
order = ["Core", "Soft-core", "Shell", "Cloud"]
print(f"ST492 genomes: {n}  |  pan-genome clusters present in >=1: {len(genes)}")
for c in order:
    print(f"  {c:10s} {int(counts.get(c, 0)):5d}")
print(f"wrote {GENES_CSV}")
print(f"wrote {MATRIX_CSV}")
