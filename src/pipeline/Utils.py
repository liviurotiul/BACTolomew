import pandas as pd
import hashlib

cat_specs_no_other = [
    ("AMR", "#c62828"),
    ("BMRG", "#1b5e20"),
    ("EFFLUX", "green"),
    ("VIRULENCE", "#6a1b9a"),
    ("MOBILE", "#ff8f00"),
    ("PHAGE", "#1565c0"),
    ("PLASMID", "#f57c00"),
    ("COMPOUND", "red"),
    ("MULTIPLE", "#ff6f00"),
    ("core-AMR", "#c62828"),
    ("plus-AMR", "#c62828"),
    ("core-BMRG", "#1b5e20"),
    ("plus-BMRG", "#1b5e20"),
    ("plus-VIRULENCE", "#6a1b9a")
]

cat_specs = [
    ("AMR", "#c62828"),
    ("BMRG", "#1b5e20"),
    ("EFFLUX", "green"),
    ("VIRULENCE", "#6a1b9a"),
    ("MOBILE", "#ff8f00"),
    ("PHAGE", "#1565c0"),
    ("PLASMID", "#f57c00"),
    ("COMPOUND", "red"),
    ("MULTIPLE", "#ff6f00"),
    ("core-AMR", "#c62828"),
    ("plus-AMR", "#c62828"),
    ("core-BMRG", "#1b5e20"),
    ("plus-BMRG", "#1b5e20"),
    ("plus-VIRULENCE", "#6a1b9a"),
    ("OTHER", "#9e9e9e")
]

cat_priority = {
    "AMR": 1,
    "core-AMR": 2,
    "plus-AMR": 2,
    "BMRG": 3,
    "plus-BMRG": 3,
    "core-BMRG": 4,
    "VIRULENCE": 5,
    "plus-VIRULENCE": 5,
    "MOBILE": 6,
    "PLASMID": 7,
    "PHAGE": 8,
    "EFFLUX": 9,
    "COMPOUND": 10,
    "MULTIPLE": 11,
    "OTHER": 12,
}

CATEGORIES = [cat[0] for cat in cat_specs]

bakta_allowed_versions = ["bakta 1.10.4"]


def assign_gene_ids(
    genes_df: pd.DataFrame,
    *,
    column: str = "GeneID",
) -> pd.DataFrame:
    """
    Add a stable identifier column that:
      - Uses Gene/ConsensusDonorProduct/Product when available.
      - Replaces DUF/hypothetical entries with unique per-locus tokens.
    """
    df = genes_df.copy()
    for c in ("Sample", "Contig", "Gene", "Product", "ConsensusDonorProduct", "Strand"):
        if c not in df.columns:
            df[c] = pd.NA

    for c in ("Sample", "Contig"):
        df[c] = df[c].astype(str)
    for c in ("Gene", "Product", "ConsensusDonorProduct"):
        df[c] = df[c].astype("string")

    # Keep Gene as provided by Bakta; do not override from marker annotations.

    df["Start"] = pd.to_numeric(df.get("Start"), errors="coerce")
    df["Stop"] = pd.to_numeric(df.get("Stop"), errors="coerce")

    base_key = df["Gene"].astype("string").str.strip()

    unk_mask = df["Product"].astype("string").str.contains(
        r"\b(?:hypothetical|DUF)\b", case=False, na=False
    )
    sig = (
        df["Sample"].astype("string").fillna("") + "|"
        + df["Contig"].astype("string").fillna("") + "|"
        + df["Gene"].astype("string").fillna("") + "|"
        + df["Product"].astype("string").fillna("") + "|"
        + df["Start"].fillna(-1).round().astype("Int64").astype(str) + "|"
        + df["Stop"].fillna(-1).round().astype("Int64").astype(str)
    )

    tokens = sig.apply(
        lambda s: "UNK:" + hashlib.blake2b(str(s).encode("utf-8"), digest_size=6).hexdigest()
    )

    gene_id = base_key.copy()
    gene_id = gene_id.mask(unk_mask, tokens)
    gene_id = gene_id.where(gene_id.notna() & (gene_id != ""), tokens)
    gene_id = gene_id.fillna(tokens)
    df[column] = gene_id
    return df
