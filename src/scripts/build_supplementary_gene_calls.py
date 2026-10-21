#!/usr/bin/env python3
"""Export per-sample virulence and resistance gene database calls (with identity &
coverage) for the genes shown in the ST492 phandango figures, as two supplementary
CSVs. Calls use the same thresholds as the figures (>=80% identity, >=80% subject
coverage) over the 93 ST492 genomes.
"""
import pandas as pd

MARKER_HITS = "/home/bals/Documents/dev/BACTolomew/cache/run_78eaa9126850/marker_hits.parquet"
METADATA = "/Data/ST492/metadata/metadata_ALL_new_merged_with_clusters_with_complex.csv"
FIG = "/Data/ST492/figures2026"
OUT = "/Data/ST492/supplementary"
MIN_ID, MIN_COV = 80.0, 80.0

# ST492 samples (the phandango genome set)
md = pd.read_csv(METADATA)
st492 = set(md.loc[md["subtype"].astype(str).str.strip().isin(["492", "492.0"]), "Sample"])

# gene -> category, taken from the figure panels (so we export exactly the genes shown)
vf_cat = pd.concat([
    pd.read_csv(f"{FIG}/st492_vf_panel_presence.csv")[["gene", "subcategory"]],
    pd.read_csv(f"{FIG}/st492_vf_core_inventory.csv")[["gene", "subcategory"]],
]).drop_duplicates("gene").set_index("gene")["subcategory"].to_dict()
res_cat = (pd.read_csv(f"{FIG}/st492_res_panel_presence.csv")[["gene", "subcategory"]]
           .drop_duplicates("gene").set_index("gene")["subcategory"].to_dict())

mh = pd.read_parquet(MARKER_HITS)
mh = mh[mh["Sample"].isin(st492)
        & (pd.to_numeric(mh["Identity"], errors="coerce") >= MIN_ID)
        & (pd.to_numeric(mh["Subject_coverage"], errors="coerce") >= MIN_COV)]

REN = {"Sample": "sample_id", "Gene": "gene", "DB": "database", "Method": "method",
       "Identity": "identity", "Subject_coverage": "coverage",
       "SubjectIdentifier": "subject_identifier", "SubjectProteinAccession": "subject_protein_accession",
       "Product": "product", "Contig": "contig", "Start": "start", "Stop": "stop",
       "Strand": "strand", "AnnotationID": "annotation_id",
       "ResistanceDrugClasses": "resistance_drug_classes",
       "ResistanceGeneFamilies": "resistance_gene_families",
       "ResistanceMechanisms": "resistance_mechanisms"}

# --- virulence (VFDB) ---
vf = mh[mh["DB"].astype(str).str.contains("VFDB", na=False) & mh["Gene"].isin(vf_cat)].rename(columns=REN)
vf["vf_category"] = vf["gene"].map(vf_cat)
vf_cols = ["sample_id", "gene", "vf_category", "database", "method", "identity", "coverage",
           "subject_identifier", "product", "contig", "start", "stop", "strand", "annotation_id"]
vf = vf[vf_cols].round({"identity": 1, "coverage": 1}).sort_values(["vf_category", "gene", "sample_id"])
vf.to_csv(f"{OUT}/st492_virulence_gene_calls.csv", index=False)

# --- resistance (AMR + BMRG) ---
res = mh[mh["DB"].astype(str).str.contains("AMR|BMRG", na=False) & mh["Gene"].isin(res_cat)].rename(columns=REN)
res["resistance_class"] = res["gene"].map(res_cat)
res_cols = ["sample_id", "gene", "resistance_class", "database", "method", "identity", "coverage",
            "resistance_drug_classes", "resistance_gene_families", "resistance_mechanisms",
            "subject_protein_accession", "subject_identifier", "product",
            "contig", "start", "stop", "strand", "annotation_id"]
res = res[res_cols].round({"identity": 1, "coverage": 1}).sort_values(["resistance_class", "gene", "sample_id"])
res.to_csv(f"{OUT}/st492_resistance_gene_calls.csv", index=False)

print(f"virulence calls: {len(vf)} rows, {vf['gene'].nunique()} genes, {vf['sample_id'].nunique()} genomes")
print(f"  databases: {sorted(vf['database'].unique())}")
print(f"resistance calls: {len(res)} rows, {res['gene'].nunique()} genes, {res['sample_id'].nunique()} genomes")
print(f"  databases: {sorted(res['database'].unique())}")
print(f"wrote {OUT}/st492_virulence_gene_calls.csv")
print(f"wrote {OUT}/st492_resistance_gene_calls.csv")
