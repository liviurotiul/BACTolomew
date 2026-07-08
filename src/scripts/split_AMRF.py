# split_amrf_fasta.py
import sys, re, os
import pandas as pd

ACC_PAT = re.compile(r"[A-Z]{1,3}[_]?\d+(?:\.\d+)?|[A-Z]{3,}\d+(?:\.\d+)?|[A-Z]{1,2}_\d+\.\d+")

def read_fasta(path):
    with open(path, "r") as fh:
        hdr, seq = None, []
        for line in fh:
            if line.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(seq)
                hdr = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if hdr is not None:
            yield hdr, "".join(seq)

def split_multi_acc(x):
    if pd.isna(x) or not str(x).strip():
        return []
    s = str(x).replace(";", " ").replace(",", " ").strip()
    return [t for t in s.split() if t != "-"]

def norm_acc(a):
    a = str(a).strip()
    if not a: return []
    base = a.split()[0]
    nover = base.split(".")[0] if "." in base else base
    return list(dict.fromkeys([base, nover]))

def sanitize_fname(scope, typ):
    s = (scope or "unknown").strip().lower()
    t = (typ or "unknown").strip().upper()
    def clean(x): return re.sub(r"[^\w.+-]+","_", x)
    return f"AMRF_{clean(s)}_{clean(t)}.fasta"

def build_indexes(df):
    cols = {
        "gb_prot": "GenBank protein",
        "gb_nuc":  "GenBank nucleotide",
        "rf_prot": "RefSeq protein",
        "rf_nuc":  "RefSeq nucleotide",
    }
    out = {k:{} for k in cols}
    for _, row in df.iterrows():
        scope = row.get("Scope", "unknown")
        typ   = row.get("Type", "unknown")
        for key, col in cols.items():
            if col not in df.columns: continue
            accs = []
            for raw in split_multi_acc(row.get(col, "")):
                accs.extend(norm_acc(raw))
            for a in accs:
                out[key].setdefault(a, (scope, typ))
    return out  # dicts per source column

def best_match(header, idxs):
    toks = []
    # 1) first token
    first = header.split()[0]
    toks.extend(norm_acc(first))
    # 2) any accession-like tokens in the full header
    for m in ACC_PAT.findall(header):
        toks.extend(norm_acc(m))
    # de-dup, preserve order
    toks = list(dict.fromkeys(toks))
    # priority: GenBank protein → GenBank nucleotide → RefSeq protein → RefSeq nucleotide
    for key in ("gb_prot","gb_nuc","rf_prot","rf_nuc"):
        d = idxs[key]
        for t in toks:
            if t in d:
                return d[t]
    return ("unknown","unknown")

def main(fasta_path, tsv_path, outdir="."):
    os.makedirs(outdir, exist_ok=True)
    meta = pd.read_csv(tsv_path, sep="\t", dtype=str).fillna("")
    idxs = build_indexes(meta)

    buckets = {}  # filename -> list[(header, seq)]
    unmatched = []

    for hdr, seq in read_fasta(fasta_path):
        scope, typ = best_match(hdr, idxs)
        fname = sanitize_fname(scope, typ)
        buckets.setdefault(fname, []).append((hdr, seq))
        if scope == "unknown" and typ == "unknown":
            unmatched.append(hdr)

    for fname, recs in buckets.items():
        with open(os.path.join(outdir, fname), "w") as out:
            for h, s in recs:
                out.write(f">{h}\n")
                for i in range(0, len(s), 60):
                    out.write(s[i:i+60] + "\n")

    rep = []
    for fname, recs in sorted(buckets.items()):
        # parse back scope/type from filename just for the report
        m = re.match(r"AMRF_(.+)_(.+)\.fasta$", fname)
        scope = m.group(1) if m else "unknown"
        typ = m.group(2) if m else "unknown"
        rep.append({"file": fname, "scope": scope, "type": typ, "count": len(recs)})
    pd.DataFrame(rep).sort_values(["scope","type"]).to_csv(os.path.join(outdir, "split_report.tsv"), sep="\t", index=False)

    if unmatched:
        with open(os.path.join(outdir, "unmatched_headers.txt"), "w") as fh:
            fh.write("\n".join(unmatched) + "\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python split_amrf_fasta.py <sequences.fasta> <metadata.tsv> [outdir]")
        sys.exit(1)
    fasta = sys.argv[1]
    tsv   = sys.argv[2]
    outd  = sys.argv[3] if len(sys.argv) > 3 else "."
    main(fasta, tsv, outd)
