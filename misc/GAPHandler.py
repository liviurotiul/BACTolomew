import pandas as pd
from pathlib import Path

def GAP_detect_columns(df: pd.DataFrame):
    cols = df.columns.astype(str).tolist()
    gene_col = "Gene" if "Gene" in cols else cols[0]
    meta_like = {
        gene_col, "Gene ID", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences",
        "Avg sequences per isolate", "Genome Fragment", "Order within Fragment", "Accessory Fragment",
        "Accessory Order with Fragment", "QC", "Min group size nuc", "Max group size nuc", "Avg group size nuc"
    }
    sample_cols = [c for c in cols if c not in meta_like]
    return gene_col, sample_cols

def GAP_load_gene_table(csv_path: Path) -> tuple[pd.DataFrame, str, list[str]]:
    df = pd.read_csv(csv_path)
    gene_col, sample_cols = GAP_detect_columns(df)
    return df, gene_col, sample_cols

def presence_matrix(df: pd.DataFrame, gene_col:str, sample_cols: list[str]) -> pd.DataFrame:
    g = df[[gene_col] + sample_cols].copy()
    m = g[sample_cols].astype(str).map(
        lambda x: 0 if (x is None or str(x).strip()=="" or str(x).strip()=="0" or str(x).strip().lower()=="nan") else 1
    )
    m.index = g[gene_col].astype(str)
    m.index.name = "Gene"
    return m