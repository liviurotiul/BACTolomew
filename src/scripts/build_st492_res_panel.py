#!/usr/bin/env python3
"""Build the ST492 resistance presence panel from the pipeline's AMR/BMRG hits.

Uses the cached per-sample marker hits (cache/<run>/marker_hits.parquet), keeps
resistance hits (AMRFinderPlus + CARD antimicrobial-resistance DBs and the
BacMet2 / AMRFinderPlus biocide-&-metal-resistance-gene DBs) at >=80% identity
and >=80% subject coverage, restricts to the 93 ST492 genomes, and keeps the
VARIABLE resistome (present in 1..92 of 93) for the phandango matrix. The
always-present intrinsic core (gyrA, parC, adeIJK, rpoB, ...) is dropped -- there
is no core-repertoire panel in this figure.

Each gene is assigned a resistance class (curated mapping below, with a
drug-class fallback for anything new). Writes a tidy CSV
(Sample,gene,subcategory,present) for plot_st492_tree_resistance_phandango.R.
"""
import pandas as pd

MARKER_HITS = "/home/bals/Documents/dev/BACTolomew/cache/run_78eaa9126850/marker_hits.parquet"
METADATA = "/Data/ST492/metadata/metadata_ALL_new_merged_with_clusters_with_complex.csv"
OUT_CSV = "/Data/ST492/figures2026/st492_res_panel_presence.csv"

MIN_IDENTITY = 80.0   # DIAMOND blastx of Bakta CDS vs the resistance DBs
MIN_COVERAGE = 80.0   # subject (reference protein) coverage
# resistance DBs: "*-AMR" = antimicrobial resistance, "*-BMRG" = biocide/metal resistance genes
RES_DB_RE = r"AMR|BMRG"

# --- curated gene -> resistance class -----------------------------------------
CAT_GENES = {
    "Beta-lactam": [
        "blaOXA-66", "blaOXA-72", "blaOXA-23", "blaOXA-82", "blaOXA-148", "blaOXA",
        "blaADC-30", "blaADC-361", "blaADC-56", "blaADC-33", "blaADC",
        "blaNDM-1", "blaTEM-1", "blaTEM-116", "blaSHV-5", "ftsI",
    ],
    "Aminoglycoside": [
        "ant(3'')-IIa", "aph(6)-Id", "aph(3'')-Ib", "aph(3')-VI", "aph(3')-VIa",
        "aph(3')-Ia", "aac(3)-Ia", "aac(6')-Ib", "aacA16", "aadA1", "aadA2", "armA", "rpsL",
    ],
    "Fluoroquinolone": ["gyrA", "parC"],
    "Tetracycline": ["tet(B)", "tetR(B)", "tetR", "rpsJ"],
    "Folate pathway": ["sul1", "sul2", "dfrA1", "dfrA12"],
    "Macrolide": ["msr(E)", "mph(E)"],
    "Colistin": ["pmrA", "pmrB", "pmrA_M12I", "lpxA", "lpxC", "lpsB"],
    "Efflux & other": [
        "adeR", "adeS", "adeL", "adeC", "adeF", "adeG", "adeH", "adeI", "adeJ", "adeK",
        "adeA", "adeB", "adeN", "adeT1", "adeT2", "abeM", "abeS", "amvA", "abaF", "abaQ",
        "baeR", "baeS", "oprD", "tolC", "csrA", "ompR", "cxpE", "opgE", "rpoB", "tuf",
    ],
    "Biocide/metal": ["merA", "merR", "nreB", "qacEdelta1"],
}
GENE_CAT = {g: cat for cat, genes in CAT_GENES.items() for g in genes}

# fallback: normalise a messy ResistanceDrugClasses string to a class
def fallback_cat(drug_classes: str, kind: str) -> str:
    if kind == "BMRG":
        return "Biocide/metal"
    s = (drug_classes or "").lower()
    if not s or s == "nan":
        return "Efflux & other"
    n = s.count(";") + 1
    if n >= 4:
        return "Efflux & other"          # broad multidrug efflux pumps
    if "beta-lactam" in s or "carbapenem" in s or "cephalosporin" in s or "penicillin" in s or "monobactam" in s:
        return "Beta-lactam"
    if "aminoglycoside" in s:
        return "Aminoglycoside"
    if "tetracycline" in s or "glycylcycline" in s:
        return "Tetracycline"
    if "sulfonamide" in s or "sulfone" in s or "diaminopyrimidine" in s or "trimethoprim" in s:
        return "Folate pathway"
    if "macrolide" in s or "lincosamide" in s or "streptogramin" in s:
        return "Macrolide"
    if "colistin" in s or "peptide" in s or "polymyxin" in s:
        return "Colistin"
    if "fluoroquinolone" in s or "quinolone" in s:
        return "Fluoroquinolone"
    return "Efflux & other"

# --- ST492 samples ---
md = pd.read_csv(METADATA)
st492 = sorted(md.loc[md["subtype"].astype(str).str.strip().isin(["492", "492.0"]), "Sample"])
n_samples = len(st492)

# --- resistance hits for ST492 ---
mh = pd.read_parquet(MARKER_HITS)
res = mh[
    mh["DB"].astype(str).str.contains(RES_DB_RE, na=False)
    & (pd.to_numeric(mh["Identity"], errors="coerce") >= MIN_IDENTITY)
    & (pd.to_numeric(mh["Subject_coverage"], errors="coerce") >= MIN_COVERAGE)
    & mh["Sample"].isin(st492)
].copy()
res["kind"] = res["DB"].astype(str).apply(lambda d: "BMRG" if "BMRG" in d else "AMR")

# per-gene: sample set + drug classes + kind
present_samples = res.groupby("Gene")["Sample"].agg(lambda s: set(s))
gene_kind = res.groupby("Gene")["kind"].agg(lambda s: "BMRG" if (s == "BMRG").all() else "AMR")
gene_drug = res.groupby("Gene")["ResistanceDrugClasses"].agg(
    lambda s: s.dropna().astype(str).mode().iat[0] if len(s.dropna()) else ""
)

def category_of(gene: str) -> str:
    if gene in GENE_CAT:
        return GENE_CAT[gene]
    return fallback_cat(gene_drug.get(gene, ""), gene_kind.get(gene, "AMR"))

# variable resistome: present in at least one but not all genomes (1..n-1)
var_genes = [g for g, ss in present_samples.items() if 1 <= len(ss) <= n_samples - 1]

rows = []
uncur = []
for g in var_genes:
    ss = present_samples[g]
    cat = category_of(g)
    if g not in GENE_CAT:
        uncur.append((g, cat))
    for s in st492:
        rows.append((s, g, cat, s in ss))
out = pd.DataFrame(rows, columns=["Sample", "gene", "subcategory", "present"])
out.to_csv(OUT_CSV, index=False)

# --- report ---
summ = (
    out[out["present"]].groupby(["subcategory", "gene"]).size().rename("n").reset_index()
    .sort_values(["subcategory", "n"], ascending=[True, False])
)
core_n = sum(1 for g, ss in present_samples.items() if len(ss) == n_samples)
print(f"samples={n_samples}  variable resistance genes kept={len(var_genes)}  "
      f"(of {present_samples.size} genes; {core_n} core/always-present dropped)")
print("by resistance class:")
print(summ.groupby("subcategory")["gene"].nunique().to_string())
if uncur:
    print("\nfallback-categorised (not in curated map):")
    for g, c in uncur:
        print(f"   {g:14s} -> {c}")
print("\nper-gene presence:")
for _, r in summ.iterrows():
    print(f"   {r['subcategory'][:16]:18s} {r['gene']:14s} {r['n']:3d}/{n_samples}")
print(f"\nwrote {OUT_CSV}")
