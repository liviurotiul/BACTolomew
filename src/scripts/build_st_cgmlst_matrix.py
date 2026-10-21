#!/usr/bin/env python3
"""Build a pairwise cgMLST allelic-similarity matrix for ST2 + ST492, block-ordered
by ST -> region -> country for a heatmap with nested group rectangles.

Reads chewBBACA allele calls, keeps the cgMLST core (loci called in >=95% of the
153 ST2+ST492 genomes), and computes pairwise allelic SIMILARITY = fraction of
co-called loci carrying an identical allele. Genomes are ordered ST492 (by region
then country) -> United Kingdom -> ST2 (by country). Writes:
  * st_cgmlst_matrix_order.csv  (sample, idx, st, region, country)
  * st_cgmlst_matrix_long.csv   (i, j, sim)   -- all 153x153 cells
for plot_st_cgmlst_matrix.R.
"""
import re
import numpy as np
import pandas as pd

ALLELES = "/home/bals/Documents/dev/BACTolomew/runs/sample_78eaa9126850/run_78eaa9126850/chewbbaca/out/results_alleles.tsv"
METADATA = "/Data/ST492/metadata/metadata_ALL_new_merged_with_clusters_with_complex.csv"
ORDER_CSV = "/Data/ST492/figures2026/st_cgmlst_matrix_order.csv"
LONG_CSV = "/Data/ST492/figures2026/st_cgmlst_matrix_long.csv"
CORE_CALLRATE = 0.95

# ST492 region grouping (granular, matching the reference figure)
REGION_492 = {
    "Serbia": "Eastern Europe", "Montenegro": "Eastern Europe", "Kosovo": "Eastern Europe",
    "Croatia": "Eastern Europe", "Slovenia": "Eastern Europe", "Romania": "Eastern Europe",
    "Hungary": "Eastern Europe",
    "Switzerland": "Central Europe", "Germany": "Central Europe", "France": "Central Europe",
    "USA": "North America",
    "Saudi Arabia": "West Asia", "Lebanon": "West Asia",
    "Scotland": "United Kingdom",
    "NA": "Unknown", "Unknown": "Unknown",
}
REGION_ORDER = ["Eastern Europe", "Central Europe", "North America", "West Asia",
                "Unknown", "United Kingdom"]

md = pd.read_csv(METADATA)
md["st"] = md["subtype"].astype(str).str.strip().str.replace(r"\.0$", "", regex=True)
keep = md[md["st"].isin(["2", "492"])].copy()
st_of = dict(zip(keep["Sample"], keep["st"]))
country_of = {s: (c if isinstance(c, str) else "NA")
              for s, c in zip(keep["Sample"], keep["country"])}

# ---- allele calls -> numeric (INF-x -> x; tags -> missing) ----
al = pd.read_csv(ALLELES, sep="\t", dtype=str).rename(columns=lambda c: "Sample" if c.lower() == "file" else c)
sub = al[al["Sample"].isin(st_of)].set_index("Sample")


def clean(v):
    if isinstance(v, str):
        if v.startswith("INF-"):
            v = v[4:]
        if re.fullmatch(r"\d+", v):
            return int(v)
    return np.nan


M = sub.map(clean)
core = M.columns[M.notna().mean(axis=0) >= CORE_CALLRATE]
Mc = M[core]
samples = list(Mc.index)
Mv = Mc.to_numpy()
n, L = Mv.shape

# ---- pairwise allelic similarity = identical / co-called ----
sim = np.eye(n)
for i in range(n):
    a = Mv[i]
    ai = ~np.isnan(a)
    for j in range(i + 1, n):
        b = Mv[j]
        both = ai & ~np.isnan(b)
        cc = int(both.sum())
        s = float(np.sum((a == b) & both)) / cc if cc else np.nan
        sim[i, j] = sim[j, i] = s

# ---- group + order ----
def region_of(s):
    if st_of[s] == "2":
        return "ST2"
    return REGION_492.get(country_of[s], "Unknown")


rows = pd.DataFrame({
    "sample": samples,
    "st": [("ST" + st_of[s]) for s in samples],
    "region": [region_of(s) for s in samples],
    "country": [country_of[s] for s in samples],
})
# order: ST492 regions (fixed order) -> ST2; within a region by country size desc; within country by sample
csize = rows.groupby(["region", "country"]).size().rename("csz").reset_index()
rows = rows.merge(csize, on=["region", "country"])
rows["st_rank"] = np.where(rows["st"] == "ST492", 0, 1)
rows["reg_rank"] = rows["region"].apply(lambda r: REGION_ORDER.index(r) if r in REGION_ORDER else 99)
rows = rows.sort_values(["st_rank", "reg_rank", "region", "csz", "country", "sample"],
                        ascending=[True, True, True, False, True, True]).reset_index(drop=True)
rows["idx"] = np.arange(1, len(rows) + 1)

# reorder similarity to match
pos = {s: k for k, s in enumerate(samples)}
perm = [pos[s] for s in rows["sample"]]
simo = sim[np.ix_(perm, perm)]

rows[["sample", "idx", "st", "region", "country"]].to_csv(ORDER_CSV, index=False)
ii, jj = np.meshgrid(np.arange(1, n + 1), np.arange(1, n + 1), indexing="ij")
pd.DataFrame({"i": ii.ravel(), "j": jj.ravel(), "sim": simo.ravel()}).to_csv(LONG_CSV, index=False)

print(f"genomes={n} (ST492={sum(st_of[s]=='492' for s in samples)}, ST2={sum(st_of[s]=='2' for s in samples)})  core loci={L}")
print(f"similarity range: {np.nanmin(simo[simo<1]):.3f} .. {np.nanmax(simo[simo<1]):.3f} (off-diagonal)")
print("region block sizes:")
print(rows.groupby(["st", "region"], sort=False).size().to_string())
print(f"wrote {ORDER_CSV}\nwrote {LONG_CSV}")
