#!/usr/bin/env python3
"""Build the supplementary Excel workbook linking every database hit back to its
Bakta annotation.

Sheet 1 (Bakta_annotations) is the master feature table: one row per Bakta
feature across all 153 genomes, keyed by AnnotationID (the "BKT" identifier).
Every following sheet holds the hits from one marker database, carrying the same
AnnotationID so a hit can be joined straight back to the Bakta feature it sits
on (VLOOKUP / XLOOKUP on column A, or a pivot).

COG is excluded: it was removed from the pipeline and its databases, but the
cached run predates that removal, so its 477,677 hits are still in the parquet.
"""
import pandas as pd

CACHE = "/home/bals/Documents/dev/BACTolomew/cache/run_78eaa9126850"
METADATA = "/Data/ST492/supplementary/ST2_ST492_isolate_metadata.csv"
OUT = "/Data/ST492/supplementary/ST2_ST492_bakta_annotations_and_database_hits.xlsx"

# DB code in the parquet -> worksheet name (Excel caps sheet names at 31 chars)
DB_SHEETS = {
    "AMRF_core-AMR": "AMRFinderPlus_core_AMR",
    "AMRF_plus-AMR": "AMRFinderPlus_plus_AMR",
    "AMRF_core-BMRG": "AMRFinderPlus_core_biocide",
    "AMRF_plus-BMRG": "AMRFinderPlus_plus_biocide",
    "CARD_AMR": "CARD",
    "BacMet2-EXP-database_BMRG": "BacMet2",
    "VFDB-setB-Prot_VIRULENCE": "VFDB_setB",
    "Phaster-Prophage-Virus-Prot_PHAGE": "PHASTER_prophage",
    "CGE-MGE-Prot_MOBILE": "MobileElementFinder",
    "ISFinder-Prot_MOBILE": "ISfinder",
    "CGE-PlasmidFinder-Prot_PLASMID": "PlasmidFinder",
}

# Refseq_ID is empty for every one of the 574,416 features, so it is not carried.
ANN_COLS = ["AnnotationID", "Sample", "Contig", "Locus", "Start", "Stop", "Strand",
            "Gene", "Product", "Hypothetical", "NCBI_NRP_ID"]

HIT_COLS = ["AnnotationID", "Sample", "Contig", "Start", "Stop", "Strand",
            "Gene", "GeneID", "Product", "Identity", "Subject_coverage",
            "Align_len", "Mismatches", "Evalue", "Bitscore", "Q_len", "S_len",
            "SubjectIdentifier", "SubjectProteinAccession", "SubjectDNAAccession",
            "SubjectAROAccession", "SubjectModelName"]

# Only carried on the sheets where they actually hold something.
RES_COLS = ["ResistanceAlleles", "ResistanceGeneFamilies", "ResistanceDrugClasses",
            "ResistanceAntibiotics", "ResistanceMechanisms"]


def nonempty(df, cols):
    """Keep only the columns that are present and not entirely blank."""
    out = []
    for c in cols:
        if c in df.columns and df[c].notna().any() and (df[c].astype(str).str.strip() != "").any():
            out.append(c)
    return out


print("reading Bakta annotations")
ann = pd.read_parquet(f"{CACHE}/annotations.parquet", columns=ANN_COLS)
ann = ann.sort_values(["Sample", "Contig", "Start"], kind="stable")
print(f"  {len(ann):,} features across {ann.Sample.nunique()} genomes")

print("reading marker hits")
mh = pd.read_parquet(f"{CACHE}/marker_hits.parquet")
before = len(mh)
mh = mh[mh.DB != "COGorg24_COG"]
print(f"  {before:,} hits, {before - len(mh):,} dropped as COG, {len(mh):,} kept")

import xlsxwriter

# constant_memory keeps peak RAM flat over the 574k-row sheet, but it flushes each
# row as soon as the next one starts, so every cell must be written in order:
# header first, then the data rows. Writing back to row 0 afterwards (as
# DataFrame.to_excel followed by a header restyle would) silently drops the data.
book = xlsxwriter.Workbook(OUT, {"constant_memory": True})
hdr = book.add_format({"bold": True, "bg_color": "#132033", "font_color": "white",
                       "border": 1, "align": "left", "valign": "vcenter"})


def write_sheet(df, name, widths):
    ws = book.add_worksheet(name)
    ws.freeze_panes(1, 1)
    for i, col in enumerate(df.columns):
        ws.set_column(i, i, widths.get(col, 14))
        ws.write(0, i, col, hdr)
    clean = df.astype(object).where(pd.notna(df), None)
    for r, row in enumerate(clean.itertuples(index=False, name=None), start=1):
        ws.write_row(r, 0, row)
    print(f"  {name:<28} {len(df):>8,} rows x {len(df.columns)} cols")


print("writing sheets")
write_sheet(ann, "Bakta_annotations",
            {"AnnotationID": 13, "Sample": 18, "Locus": 15, "Product": 52,
             "NCBI_NRP_ID": 16, "Refseq_ID": 14, "Gene": 12})

for db, sheet in DB_SHEETS.items():
    sub = mh[mh.DB == db].copy()
    if sub.empty:
        print(f"  {sheet:<28} no hits, sheet skipped")
        continue
    cols = nonempty(sub, HIT_COLS) + nonempty(sub, RES_COLS)
    sub = sub[cols].sort_values(["Sample", "Contig", "Start"], kind="stable")
    write_sheet(sub, sheet,
                {"AnnotationID": 13, "Sample": 18, "Product": 46, "Gene": 16,
                 "SubjectIdentifier": 34, "SubjectModelName": 22,
                 "ResistanceAlleles": 34, "ResistanceGeneFamilies": 40,
                 "ResistanceDrugClasses": 40, "ResistanceMechanisms": 30,
                 "ResistanceAntibiotics": 30})

# Isolate metadata last, so the Bakta table stays sheet 1. Its sample_id column
# matches the Sample column on every other sheet (153/153), so it joins cleanly.
meta = pd.read_csv(METADATA)
write_sheet(meta, "Isolate_metadata",
            {"sample_id": 20, "isolation_source": 18, "isolation_source_detail": 22,
             "geographic_region": 18, "cgmlst_complex_type": 19, "ncbi_biosample": 16,
             "ncbi_bioproject": 16, "ncbi_assembly_name": 20, "ncbi_organism_name": 24,
             "ncbi_infraspecific_strain": 24, "ncbi_collected_by": 26,
             "ncbi_geo_loc_name": 20, "ncbi_isolation_source": 22,
             "ncbi_host_disease": 20, "ncbi_host_health_state": 22})

book.close()
import os
print(f"\nwrote {OUT}  ({os.path.getsize(OUT)/1e6:.1f} MB)")
