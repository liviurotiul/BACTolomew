import re
import pandas as pd

NWK_PATH   = "/home/bals/Documents/test-Contexar3/trees/BestTree_Nov2025.nwk"
META_PATH  = "/Data/ST492/metadata/metadata_ALL_new_merged.csv"  # metadata with Sample column
OUT_PATH   = "/home/bals/Documents/test-Contexar/trees/vBovin_Best_RaxML_Phylo_ST2_ST492_fixed.nwk"

def base_and_ver(s: str):
    m = re.fullmatch(r"(GCA_\d+)(?:\.(\d+))?", s.strip())
    if not m:
        return None, None
    base, ver = m.group(1), m.group(2)
    return base, (int(ver) if ver is not None else None)

def strip_suffix(name: str) -> str:
    return re.sub(r"(?:_ST(?:492|2))$", "", name)

# Build version map (base -> highest version) from metadata
meta = pd.read_csv(META_PATH, dtype=str).fillna("")
sample_col = "Sample" if "Sample" in meta.columns else ("accession" if "accession" in meta.columns else None)
if sample_col is None:
    raise SystemExit("Metadata must contain 'Sample' or 'accession' column.")

ver_map = {}
for s in meta[sample_col].astype(str):
    base, ver = base_and_ver(s)
    if base and ver is not None and (base not in ver_map or ver > ver_map[base]):
        ver_map[base] = ver

def transform_label(lbl: str) -> str:
    lbl2 = strip_suffix(lbl)
    base, ver = base_and_ver(lbl2)
    if base and ver is None and base in ver_map:
        return f"{base}.{ver_map[base]}"
    return lbl2

# Replace labels in Newick: tokens after '(' or ',' up to next ':', ',', ')' or ';'
pattern = re.compile(r'(?<=\(|,)([^:\(\),;]+)')

with open(NWK_PATH, "r", encoding="utf-8") as f:
    newick = f.read()

def repl(m: re.Match) -> str:
    token = m.group(1).strip()
    return transform_label(token)

fixed = pattern.sub(repl, newick)

with open(OUT_PATH, "w", encoding="utf-8") as f:
    f.write(fixed)

print(f"Saved fixed Newick -> {OUT_PATH}")

