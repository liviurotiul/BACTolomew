# --- minimal header-cleanup ----------------------------------------------
path = "/Data/DatabasesFasta/VFDB_setB_Prot.fasta"
out  = path.rsplit(".", 1)[0] + "_simplified.fasta"

with open(path) as fin, open(out, "w") as fout:
    for ln in fin:
        if not ln.startswith(">"):
            fout.write(ln)               # plain sequence line
            continue

        # examples:
        # >VFG037170(gb|WP_001081754) (plc1) phospholipase C [...]
        # >VFG037179(gb|YP_001085084) (plc2) phospholipase C [...]
        try:
            fam      = ln[1:].split("(", 1)[0]                 # VFG037170
            acc      = ln.split("|", 1)[0].split(")")[0]       # WP_001081754
            product  = ln.split(")", 2)[2].split("[", 1)[0].strip()
        except:
            print("⚠ could not parse header:", ln.strip())

        fout.write(f">{fam}|{acc}|{product}\n")

print("✔ headers simplified →", out)
