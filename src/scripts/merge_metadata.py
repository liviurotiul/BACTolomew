import re
import pandas as pd

CSV1_PATH = "/Data/ST492/metadata/metadata.csv"
CSV2_PATH = "/Data/ST492/metadata/metadata_492_new.csv"
OUT_PATH  = "/Data/ST492/metadata/metadata_492_new_merged.csv"

def base_and_ver(gca: str):
    m = re.match(r"^(GCA_\d+)(?:\.(\d+))?$", str(gca).strip())
    if not m:
        return None, None
    base, ver = m.group(1), m.group(2)
    return base, int(ver) if ver is not None else None

def year_only(x: str) -> str:
    if x is None:
        return ""
    s = str(x).strip()
    if not s or s.lower() in {"nan", "none"}:
        return ""
    # prefer a leading YYYY
    m = re.match(r"^\s*(\d{4})", s)
    if m:
        return m.group(1)
    # fallback: any 4-digit year
    m = re.search(r"(\d{4})", s)
    return m.group(1) if m else ""

# csv1
c1 = pd.read_csv(CSV1_PATH, dtype=str)

# rename columns in csv1
rename_map_c1 = {}
if "geolocation" in c1.columns:
    rename_map_c1["geolocation"] = "Country"
if "isolation_date" in c1.columns:
    rename_map_c1["isolation_date"] = "collection_date"
c1 = c1.rename(columns=rename_map_c1)

# normalize year in csv1 collection_date
if "collection_date" in c1.columns:
    c1["collection_date"] = c1["collection_date"].map(year_only)

c1["Sample"] = c1["Sample"].astype(str)

# build version map from csv1 Samples
ver_map = {}
for s in c1["Sample"].dropna().astype(str):
    base, ver = base_and_ver(s)
    if base and ver is not None and (base not in ver_map or ver > ver_map[base]):
        ver_map[base] = ver

# csv2 (tsv)
c2 = pd.read_csv(CSV2_PATH, sep="\t", dtype=str).fillna("")
if "accession" not in c2.columns:
    raise SystemExit("csv2 must have an 'accession' column")

def patch_acc(acc):
    base, ver = base_and_ver(acc)
    if not base or ver is not None:
        return acc
    return f"{base}.{ver_map[base]}" if base in ver_map else acc

c2["accession"] = c2["accession"].map(patch_acc)
c2 = c2.rename(columns={"accession": "Sample"})

# add csv1 rows missing from csv2 by Sample
c2_samples = set(c2["Sample"].astype(str))
c1_only = c1[~c1["Sample"].astype(str).isin(c2_samples)]

# merge and clean empty columns
final = pd.concat([c2, c1_only], ignore_index=True, sort=False)
final = final.fillna("")
final = final.dropna(axis=1, how="all")
final = final.loc[:, (final != "").any(axis=0)]

final.to_csv(OUT_PATH, index=False)
print(f"Saved cleaned versioned CSV -> {OUT_PATH}  (rows: {len(final)})")
