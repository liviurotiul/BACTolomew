#!/usr/bin/env python3
"""Extract the phage (vOTU) results from Phages.xlsx into tidy CSVs for plotting
and for supplementary tables.

Writes to /Data/ST492/figures2026:
  st492_phage_votu_catalog.csv     one row per vOTU (57) with all attributes
  st492_phage_genome_counts.csv    one row per genome (152) with its vOTU load
  st492_phage_genome_votu_long.csv genome x vOTU membership, long format
  st492_phage_species_summary.csv  one row per predicted phage species
"""
import pandas as pd

XLSX = "/Data/ST492/marius-phage/Phages.xlsx"
OUT = "/Data/ST492/figures2026"

FUNC_COLS = ["connector", "head and packaging", "integration and excision", "lysis",
             "moron, auxiliary metabolic gene and host takeover", "other", "tail",
             "transcription regulation", "unknown function",
             "DNA, RNA and nucleotide metabolism", "tRNAs", "tmRNAs", "CRISPRs"]

# ---------------------------------------------------------------- vOTU catalog
m = pd.read_excel(XLSX, sheet_name="vOTU_Master_Classified")
m = m[m["vOTU_ID"].notna() & m["length"].notna()].copy()      # drop max/min summary rows

m["st492_n"] = pd.to_numeric(m["ST492"], errors="coerce").fillna(0).astype(int)
m["st2_n"] = pd.to_numeric(m["ST2"], errors="coerce").fillna(0).astype(int)
m["st_group"] = m.apply(
    lambda r: "ST492 only" if r.st492_n > 0 and r.st2_n == 0
    else ("ST2 only" if r.st2_n > 0 and r.st492_n == 0 else "Shared"), axis=1)
m["species"] = m["Description"].fillna("Unassigned")
m["length_kb"] = (m["length"] / 1000).round(1)

cat = pd.DataFrame({
    "votu_id": m["vOTU_ID"],
    "length_bp": m["length"].astype(int),
    "length_kb": m["length_kb"],
    "checkv_quality": m["checkv_quality"],
    "completeness_pct": pd.to_numeric(m["completeness"], errors="coerce").round(2),
    "provirus": m["provirus"],
    "species": m["species"],
    "genus": m["Genus"].fillna("Unclassified"),
    "novelty": m["Novelty_Classification"].fillna("Unknown"),
    "mash_distance": pd.to_numeric(m["mash_distance"], errors="coerce"),
    "st_group": m["st_group"],
    "n_st492_genomes": m["st492_n"],
    "n_st2_genomes": m["st2_n"],
    "n_genomes_total": pd.to_numeric(m["genome_count"], errors="coerce").fillna(0).astype(int),
    "cds": pd.to_numeric(m["CDS"], errors="coerce"),
    "hypothetical_pct": pd.to_numeric(m["HP % from CDS"], errors="coerce").round(1),
    "n_amr_genes": pd.to_numeric(m["CARD_AMR_Genes"], errors="coerce").fillna(0).astype(int),
    "n_virulence_factors": pd.to_numeric(m["VFDB_Virulence_Factors"], errors="coerce").fillna(0).astype(int),
})
for c in FUNC_COLS:
    if c in m.columns:
        key = "n_" + c.split(",")[0].replace(" ", "_").replace("/", "_").lower()
        cat[key] = pd.to_numeric(m[c], errors="coerce").fillna(0).astype(int)
cat = cat.sort_values("votu_id")
cat.to_csv(f"{OUT}/st492_phage_votu_catalog.csv", index=False)

# ------------------------------------------------------- per-genome vOTU load
g = pd.read_excel(XLSX, sheet_name="Genome_to_vOTUs_mapping")
g = g[g["genome_id"].notna() & g["ST"].isin([2, 492, 2.0, 492.0])].copy()
vcols = ["vOTU_list"] + [c for c in g.columns if c.startswith("Unnamed")]
g["votus"] = g[vcols].apply(lambda r: [x for x in r if isinstance(x, str) and x.startswith("vOTU")], axis=1)
counts = pd.DataFrame({
    "genome_id": g["genome_id"].astype(str),
    "st": "ST" + g["ST"].astype(int).astype(str),
    "n_votus": g["votus"].apply(len),
    "votu_list": g["votus"].apply(lambda v: ";".join(v)),
}).sort_values(["st", "genome_id"])

# attach capsule type and assembly contiguity (ST492 genomes) from the NCBI sheet
nx = pd.read_excel(XLSX, sheet_name="NCBI_METADATA_ST492")
nx["gid"] = nx["accession"].astype(str).str.replace(r"\.\d+$", "", regex=True)
nx = nx[["gid", "K locus", "Number of scaffolds", "Country", "isolation_source"]].rename(
    columns={"K locus": "k_locus", "Number of scaffolds": "n_scaffolds",
             "Country": "country", "isolation_source": "isolation_source"})
counts["gid"] = counts["genome_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
counts = counts.merge(nx, on="gid", how="left").drop(columns="gid")
counts.to_csv(f"{OUT}/st492_phage_genome_counts.csv", index=False)

# membership, long format
long = counts.assign(votu_id=counts["votu_list"].str.split(";")).explode("votu_id")
long = long[long["votu_id"].astype(str).str.startswith("vOTU")][["genome_id", "st", "votu_id"]]
long.to_csv(f"{OUT}/st492_phage_genome_votu_long.csv", index=False)

# ----------------------------------------------------------- species summary
sp = (cat.groupby("species")
        .agg(n_votus=("votu_id", "size"),
             n_genomes=("n_genomes_total", "sum"),
             n_st492=("n_st492_genomes", "sum"),
             n_st2=("n_st2_genomes", "sum"),
             median_length_kb=("length_kb", "median"))
        .reset_index().sort_values(["n_votus", "n_genomes"], ascending=False))
sp.to_csv(f"{OUT}/st492_phage_species_summary.csv", index=False)

# --------------------------------------------------- vOTU vs capsule association
import numpy as np
from scipy import stats as _st
from statsmodels.stats.multitest import multipletests

c492 = counts[(counts["st"] == "ST492") & counts["k_locus"].notna()][["genome_id", "k_locus"]]
l492 = long[long["st"] == "ST492"].merge(c492, on="genome_id")
pa = (l492.assign(v=1).pivot_table(index="genome_id", columns="votu_id", values="v", fill_value=0)
      .reindex(c492["genome_id"]).fillna(0).astype(int))
kl = c492.set_index("genome_id")["k_locus"]

rows = []
for v in pa.columns:
    for k in sorted(kl.unique()):
        a = int(((pa[v] == 1) & (kl == k)).sum()); b = int(((pa[v] == 1) & (kl != k)).sum())
        cc = int(((pa[v] == 0) & (kl == k)).sum()); d = int(((pa[v] == 0) & (kl != k)).sum())
        odds, pv = _st.fisher_exact([[a, b], [cc, d]])
        rows.append(dict(votu_id=v, k_locus=k, n_present=a, kl_n=a + cc,
                         pct=100 * a / (a + cc) if (a + cc) else 0.0,
                         odds_ratio=odds, p=pv))
assoc = pd.DataFrame(rows)
assoc["q"] = multipletests(assoc["p"], method="fdr_bh")[1]
assoc["log2_or"] = np.log2(assoc["odds_ratio"].replace(0, 1 / 128).clip(upper=128))
assoc.to_csv(f"{OUT}/st492_phage_votu_kl_association.csv", index=False)

pa_out = pa.reset_index().merge(c492, on="genome_id")
pa_out.to_csv(f"{OUT}/st492_phage_presence_matrix.csv", index=False)
print(f"vOTU-KL associations : {len(assoc)} tested, {(assoc.q < 0.05).sum()} significant (BH q<0.05)")

print(f"vOTUs               : {len(cat)}")
print(f"  ST492 only        : {(cat.st_group=='ST492 only').sum()}")
print(f"  ST2 only          : {(cat.st_group=='ST2 only').sum()}")
print(f"  Shared            : {(cat.st_group=='Shared').sum()}")
print(f"  carrying an ARG   : {(cat.n_amr_genes>0).sum()}")
print(f"  carrying a VF     : {(cat.n_virulence_factors>0).sum()}")
print(f"  length range      : {cat.length_bp.min():,} to {cat.length_bp.max():,} bp")
print(f"genomes             : {len(counts)}  ({(counts.st=='ST492').sum()} ST492, {(counts.st=='ST2').sum()} ST2)")
print(f"  vOTUs per genome  : ST492 median {counts[counts.st=='ST492'].n_votus.median():.0f} "
      f"(range {counts[counts.st=='ST492'].n_votus.min()}-{counts[counts.st=='ST492'].n_votus.max()}), "
      f"ST2 median {counts[counts.st=='ST2'].n_votus.median():.0f} "
      f"(range {counts[counts.st=='ST2'].n_votus.min()}-{counts[counts.st=='ST2'].n_votus.max()})")
print(f"species             : {len(sp)}")
print(f"\nwrote 4 CSVs to {OUT}")
