#!/usr/bin/env python3
"""Build the ST492 virulence-factor presence panel from the pipeline's VFDB hits.

Uses the cached per-sample VFDB-setB marker hits (cache/<run>/marker_hits.parquet),
keeps hits >=80% identity and >=80% subject coverage, restricts to the 93 ST492
genomes, and splits the virulence genes into:
  * VARIABLE  - absent from at least one genome (present in 1..92 of 93) -> phandango matrix
  * CORE      - present in ALL 93 genomes -> inventory strip below the plot
Each gene's virulence subcategory is VFDB's own functional class (Adherence,
Biofilm, Effector delivery system, Nutritional/Metabolic factor, Immune
modulation, Exotoxin, Regulation), parsed from the VFDB FASTA headers. Writes a
tidy CSV (Sample,gene,subcategory,present) for the phandango figure.
"""
import re
import pandas as pd

MARKER_HITS = "/home/bals/Documents/dev/BACTolomew/cache/run_78eaa9126850/marker_hits.parquet"
METADATA = "/Data/ST492/metadata/metadata_ALL_new_merged_with_clusters_with_complex.csv"
VFDB_FASTA = "/home/bals/Documents/dev/BACTolomew/DB/seq/markers/VFDB-setB-Prot_VIRULENCE.fasta"
OUT_CSV = "/Data/ST492/figures2026/st492_vf_panel_presence.csv"

MIN_IDENTITY = 80.0   # DIAMOND blastx of Bakta CDS vs VFDB-setB
MIN_COVERAGE = 80.0   # subject (VFDB protein) coverage
CORE_MIN = 91         # present in >90 genomes -> core repertoire; the rest (1..90) -> phandango

# --- ST492 samples ---
md = pd.read_csv(METADATA)
st492 = sorted(md.loc[md["subtype"].astype(str).str.strip().isin(["492", "492.0"]), "Sample"])
n_samples = len(st492)

# --- VFG -> VFDB functional category (from FASTA headers) ---
# header: >VFG037170(gb|...) (plc1) phospholipase C [Phospholipase C (VF0470) - Exotoxin (VFC0235)] [organism]
vfg_re = re.compile(r"^>(VFG\d+)")
cat_re = re.compile(r"-\s*([A-Za-z/ ,]+?)\s*\(VFC\d+\)")
cat_of = {}
with open(VFDB_FASTA) as fh:
    for line in fh:
        if line.startswith(">"):
            m1 = vfg_re.match(line)
            if not m1:
                continue
            m2 = cat_re.search(line)
            cat_of[m1.group(1)] = m2.group(1).strip() if m2 else "Other"

# --- VFDB virulence hits for ST492 ---
mh = pd.read_parquet(MARKER_HITS)
vir = mh[
    mh["DB"].astype(str).str.contains("VFDB", na=False)
    & (pd.to_numeric(mh["Identity"], errors="coerce") >= MIN_IDENTITY)
    & (pd.to_numeric(mh["Subject_coverage"], errors="coerce") >= MIN_COVERAGE)
    & mh["Sample"].isin(st492)
].copy()
vir["VFG"] = vir["SubjectIdentifier"].astype(str).str.extract(r"(VFG\d+)")
vir["VFcat"] = vir["VFG"].map(cat_of).fillna("Other")

# per-gene: sample set + dominant category
present_samples = vir.groupby("Gene")["Sample"].agg(lambda s: set(s))
gene_cat = vir.groupby("Gene")["VFcat"].agg(lambda s: s.mode().iat[0] if len(s.mode()) else "Other")

# variable genes: present in at least one genome but not "core" (1..CORE_MIN-1)
var_genes = [g for g, ss in present_samples.items() if 1 <= len(ss) < CORE_MIN]

rows = []
for g in var_genes:
    ss = present_samples[g]
    cat = gene_cat[g]
    for s in st492:
        rows.append((s, g, cat, s in ss))
out = pd.DataFrame(rows, columns=["Sample", "gene", "subcategory", "present"])
out.to_csv(OUT_CSV, index=False)

# --- core inventory (present in >90 genomes; not shown per-sample) ---
CORE_CSV = "/Data/ST492/figures2026/st492_vf_core_inventory.csv"
core = [(g, gene_cat[g], len(ss)) for g, ss in present_samples.items() if len(ss) >= CORE_MIN]
core_df = pd.DataFrame(core, columns=["gene", "subcategory", "present_count"])
core_df = core_df.sort_values(["subcategory", "gene"])
core_df.to_csv(CORE_CSV, index=False)
print(f"core inventory: {len(core_df)} genes (present in >={CORE_MIN} of {n_samples}) -> {CORE_CSV}")

# --- report ---
summ = (
    out[out["present"]].groupby(["subcategory", "gene"]).size().rename("n").reset_index()
    .sort_values(["subcategory", "n"], ascending=[True, False])
)
print(f"samples={n_samples}  variable VF genes kept={len(var_genes)}  (of {present_samples.size} VFDB genes)")
print("by VFDB category:")
print(summ.groupby("subcategory")["gene"].nunique().to_string())
print("\nper-gene presence:")
for _, r in summ.iterrows():
    print(f"   {r['subcategory'][:26]:28s} {r['gene']:10s} {r['n']:3d}/{n_samples}")
print(f"\nwrote {OUT_CSV}")
