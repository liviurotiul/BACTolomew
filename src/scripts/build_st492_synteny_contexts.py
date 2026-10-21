#!/usr/bin/env python3
"""For each key ARG/porin, extract its REAL genomic neighbourhood (from Bakta GFF3s)
across ALL ST492 genomes that carry it, group the neighbourhoods into distinct
context variants, and emit the common variants aligned on the focal gene - one block
per gene-of-interest, for per-locus gggenes figures showing within-ST492 variation.
Writes /Data/ST492/figures2026/st492_synteny_variants.csv
"""
import os, re, difflib
from collections import Counter, defaultdict
from urllib.parse import unquote
import pandas as pd

BK = "/home/bals/Documents/dev/BACTolomew/runs/sample_78eaa9126850/run_78eaa9126850/bakta"
MH = "/home/bals/Documents/dev/BACTolomew/cache/run_78eaa9126850/marker_hits.parquet"
MD = "/Data/ST492/metadata/metadata_ALL_new_merged_with_clusters_with_complex.csv"
OUT = "/Data/ST492/figures2026/st492_synteny_variants.csv"
WINDOW = 14         # CDS each side (display) - wide enough to show a full integron /
                    # resistance island even when the focal gene sits at a contig edge and
                    # the whole cassette lies to one side of it
SIG_WINDOW = 9      # CDS each side used to define a context "variant"
MAX_VARIANTS = 5    # top variants shown per gene
MIN_N = 1           # min genomes to show a variant

md = pd.read_csv(MD); md["st"] = md["subtype"].astype(str).str.strip().str.replace(r"\.0$", "", regex=True)
st_of = dict(zip(md["Sample"], md["st"])); ctry = dict(zip(md["Sample"], md["country"].fillna("NA")))
st492 = set(md.loc[md["st"] == "492", "Sample"])
mh = pd.read_parquet(MH)
# COG was removed from this pipeline but the cached run predates that, so its rows are
# still in marker_hits and would be counted as extra copies of a gene.
mh = mh[mh["DB"].astype(str) != "COGorg24_COG"]

def has_gff(s): return os.path.exists(f"{BK}/{s}/{s}.gff3")

MOB = "/home/bals/Documents/dev/BACTolomew/runs/sample_78eaa9126850/run_78eaa9126850/mobsuite/per_sample"
_mob = {}
def mob_call(sample, contig):
    """mob_recon molecule type for a contig -> (molecule_type, primary_cluster_id, rep_type)."""
    if sample not in _mob:
        d, rep_cl = {}, set()
        fp = f"{MOB}/{sample}/contig_report.txt"
        if os.path.exists(fp):
            m = pd.read_csv(fp, sep="\t", dtype=str)
            for _, r in m.iterrows():
                rt = str(r["rep_type(s)"])
                d[str(r["contig_id"])] = (r["molecule_type"], str(r["primary_cluster_id"]),
                                          int(float(r["size"])), rt)
                # which plasmid clusters have replicon evidence somewhere in this genome
                if r["molecule_type"] == "plasmid" and rt not in ("-", "nan", "None", ""):
                    rep_cl.add(str(r["primary_cluster_id"]))
        _mob[sample] = d; _mob_rep[sample] = rep_cl
    return _mob[sample].get(str(contig), (None, None, 0, None))

_mob_rep = {}
def cluster_has_replicon(sample, cluster):
    """Does this plasmid cluster carry a replicon on ANY contig of this genome?

    Checking the focal contig alone is too strict: a plasmid split over several contigs
    carries its rep gene on only one of them. Checking the cluster is enough to
    corroborate the bin, and it is what separates a real plasmid from a pseudo-cluster
    such as AE271, which has no replicon, relaxase, MPF or oriT on any of its contigs
    in any genome and turns out to be chromosome broken into pieces.
    """
    if sample not in _mob_rep: mob_call(sample, None)
    return cluster in _mob_rep.get(sample, set())

_gff = {}
def read_gff(sample):
    if sample in _gff: return _gff[sample]
    rows = []
    with open(f"{BK}/{sample}/{sample}.gff3") as fh:
        for line in fh:
            if line.startswith("#") or "\tCDS\t" not in line: continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9: continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            rows.append((f[0], int(f[3]), int(f[4]), f[6], a.get("product", a.get("Name", ""))))
    d = pd.DataFrame(rows, columns=["contig", "start", "end", "strand", "product"]).sort_values(["contig", "start"])
    _gff[sample] = d; return d

_fna = {}
def read_fna(sample):
    """Contig sequences, for detecting assembly fragments that duplicate a longer contig."""
    if sample in _fna: return _fna[sample]
    d, name = {}, None
    fp = f"{BK}/{sample}/{sample}.fna"
    if os.path.exists(fp):
        with open(fp) as fh:
            for line in fh:
                if line.startswith(">"): name = line[1:].split()[0]; d[name] = []
                elif name: d[name].append(line.strip())
        d = {k: "".join(v).upper() for k, v in d.items()}
    _fna[sample] = d; return d

_RC = str.maketrans("ACGTN", "TGCAN")
def _rc(x): return x.translate(_RC)[::-1]

def drop_redundant(sample, copies):
    """Remove copies sitting on a short contig whose sequence is already contained,
    forwards or reverse-complemented, in a longer contig that carries the same gene.
    Short assembly fragments of one physical locus would otherwise be counted as
    separate copies and can form an entire spurious context variant."""
    if len(copies) < 2: return copies
    seqs = read_fna(sample)
    if not seqs: return copies
    withlen = sorted(copies, key=lambda c: -len(seqs.get(c[0], "")))
    keep = []
    for c in withlen:
        sh = seqs.get(c[0], "")
        if not sh: keep.append(c); continue
        # Test the middle 80 percent, not the whole contig: these fragments typically
        # overhang the parent contig end by a few tens of bases, so exact containment of
        # the full sequence fails even when the locus is plainly the same one.
        # Probe-based containment. A single exact substring test fails on two real cases:
        # a fragment nearly as long as its parent (rejected by a length gate) and one that
        # differs by a single base (rejected by exact matching). Sampling 200 bp probes
        # every 400 bp across the core and requiring 90 percent of them to be found
        # tolerates both. Verified to catch 9 of 9 known-redundant intI1 contigs.
        core = sh[len(sh) // 10: len(sh) - len(sh) // 10]
        probes = [core[i:i + 200] for i in range(0, max(1, len(core) - 200), 400)]
        probes = [q for q in probes if len(q) == 200]
        red = False
        for k in keep:
            if k[0] == c[0] or not probes: continue
            lg = seqs.get(k[0], "")
            if len(sh) > len(lg): continue
            found = sum(1 for q in probes if q in lg or _rc(q) in lg)
            if found >= 0.9 * len(probes):
                red = True; break
        if red: continue
        keep.append(c)
    return sorted(keep)

def label_cat(product, length=None):
    p = unquote(product).replace("%2C", ",").strip(); pl = p.lower()
    # ISCR1 is emitted by Bakta only as "IS91 family transposase". It is 1,542 bp /
    # 513 aa and matches the CGE-MGE ISCR1 reference at 100 percent over 513/513
    # residues, with no ISfinder hit. The generic label hides the rolling-circle
    # mechanism that explains these cassette structures.
    if length and 1500 <= length <= 1600 and "is91" in pl and "transposase" in pl:
        return "IS91/ISCR1", "MGE"
    if "pilin" in pl and "fimbrial" in pl: return "fimA", "other"
    if "usher" in pl or "fimd" in pl: return "fimD", "other"
    if "chaperone papd" in pl or "papd" in pl: return "fimC", "other"
    if "hypothetical" in pl: return "HP", "HP"
    m = re.findall(r"OXA-(\d+)", p)
    if m: return f"blaOXA-{m[-1]}", "ARG"
    if "ndm" in pl: return "blaNDM-1", "ARG"
    if "armа" in pl or "arma" in pl: return "armA", "ARG"
    if "sul1" in pl: return "sul1", "ARG"
    if "sul2" in pl: return "sul2", "ARG"
    if "qace" in pl: return "qacEΔ1", "ARG"
    if re.search(r"aac\(3\)", pl): return "aac(3)-Ia", "ARG"
    if re.search(r"aac\(6", pl): return "aac(6')-Ib", "ARG"
    if "aada2" in pl: return "aadA2", "ARG"
    if "aada1" in pl: return "aadA1", "ARG"
    if "aada5" in pl: return "aadA5", "ARG"
    if re.search(r"aph\(3", pl): return "aph(3')", "ARG"
    if "dfra12" in pl: return "dfrA12", "ARG"
    if "dfra" in pl: return "dfrA", "ARG"
    if re.search(r"tetracycline.*(efflux|resistance)|tet\(", pl): return "tet", "ARG"
    if "mph(e)" in pl or "macrolide 2'-phosphotransferase" in pl: return "mph(E)", "ARG"
    if "msr(e)" in pl or "abc-f" in pl: return "msr(E)", "ARG"
    # class 1 integron integrase - MUST precede any generic integrase/transposase rule
    if "inti1" in pl or "integron integrase" in pl or "class 1 integron" in pl:
        return "intI1", "MGE"
    # named insertion sequences - surface the SPECIFIC element from the Bakta product
    # e.g. "IS4 family ISEc29 transposase" -> ISEc29 ; "IS91 family transposase" -> IS91
    if re.search(r"transposase|insertion sequence", pl) or re.search(r"\bIS[A-Za-z]*\d", p):
        m = (re.search(r"(?:family|element)\s+(IS[A-Za-z]*\d+[A-Za-z]*)", p)
             or re.search(r"\b(IS[A-Za-z]*\d+[A-Za-z]*)\b", p))
        if m: return m.group(1), "MGE"
        mt = re.search(r"\b(Tn\d+)\b", p)
        return (mt.group(1) if mt else "tnp"), "MGE"
    # transposon resolvases / site-specific recombinases (Tn21, Tn1331, Tn3 ...)
    if "resolvase" in pl or "recombinase" in pl or "tnpr" in pl or "transposon" in pl:
        mt = re.search(r"\b(Tn\d+)\b", p)
        return (mt.group(1) if mt else "tnpR"), "MGE"
    # antitoxin MUST be tested first: Bakta's "Toxin-antitoxin system, antitoxin component"
    # starts with "toxin" at position 0, where the (?<!anti) lookbehind succeeds because
    # nothing precedes it, so the toxin rule would swallow the antitoxin.
    if re.search(r"brna|antitoxin|hicb|hica", pl): return "brnA", "TA"
    if re.search(r"brnt|toxin", pl): return "brnT", "TA"
    if "CarO" in p or "carbapenem-associated resistance" in pl: return "carO", "porin"  # case-sensitive: avoid "carOtene"
    if "oprd" in pl:
        # Bakta applies the name oprD to four proteins of 438, 426, 416 and 395 codons.
        # Only "OprD family outer membrane porin" (1,317 bp) is oprD itself; the three
        # "Outer membrane porin, OprD family" products are divergent family members
        # sharing 21 to 26 percent amino-acid identity with it and 44 percent with each
        # other, so they must not be drawn or counted as copies of oprD.
        if pl.startswith("outer membrane porin, oprd"): return "OprD-family", "porin"
        return "oprD", "porin"
    if re.search(r"\bompa\b", pl): return "ompA", "porin"
    for key, lab in [("replication initiation protein repm", "repM"), ("repm", "repM"),
                     ("dna replication", "DNA rep"), ("septicolysin", "septicolysin"),
                     ("znud2", "znuD2"), ("tonb", "tonB"), ("mnat", "mnaT"), ("fxsa", "fxsA"),
                     ("queg", "queG"), ("biob", "bioB"), ("ruvc", "ruvC"), ("hipb", "hipB"),
                     ("permease", "permease"), ("rps27", "rPS27A"), ("nnr", "nnr2"),
                     ("cop", "copG"), ("lysr", "lysR"), ("mg chelatase", "Mg chelatase"),
                     ("mg-chelatase", "Mg chelatase"), ("rhat", "rhaT"), ("dna-binding", "DNA BP"),
                     ("tetrahydropyridine-2,6-dicarboxylate", "dapD"), ("succinyldiaminopimelate", "dapD")]:
        if key in pl: return lab, "other"
    if "duf" in pl:
        m = re.search(r"DUF\d+", p); return (m.group(0) if m else "DUF") + " DCP", "other"
    if "domain-containing protein" in pl:
        return (p.split(" domain")[0][:10] + " DCP"), "other"
    if "methionine adenosyltransferase" in pl: return "metK", "other"
    if "eal domain" in pl or "ggdef" in pl or "diguanylate" in pl: return "eAL", "other"
    m = re.search(r"\b([A-Z][a-z]{1,4}[A-Z][A-Za-z0-9]?)\b\s*$", p)
    if m:
        s = m.group(1); return s[0].lower() + s[1:], "other"
    base = re.sub(r"\bfamily\b", "", p)          # drop only the word "family", keep the function
    base = re.sub(r"\s+protein$", "", base)
    base = re.sub(r"\s*[,(].*", "", base)         # drop from first comma/paren
    base = re.sub(r"\s{2,}", " ", base).strip()
    if len(base) > 26:                            # keep a couple of words, don't stop at the first
        cut = base[:26].rsplit(" ", 1)[0]
        base = cut if len(cut) >= 12 else base[:25] + "…"
    return base, "other"

def locate(goi, sample):
    """Every distinct copy of the focal gene in this genome, as (contig, gff_start).

    marker_hits carries one row per CDS per database, so the same copy appears under
    both AMRF_core-AMR and CARD_AMR; deduplicate on AnnotationID before iterating,
    otherwise a two-copy genome looks like four hits and iloc[0] silently keeps one.
    """
    g = read_gff(sample)
    hit = mh[(mh["Sample"] == sample) & (mh["Gene"].astype(str) == goi)]
    out = []
    if len(hit):
        for _, h in hit.drop_duplicates(subset="AnnotationID").iterrows():
            contig = f"contig_{h['Contig']}"; s0 = min(int(h["Start"]), int(h["Stop"]))
            gc = g[(g["contig"] == contig) & (abs(g["start"] - s0) < 60)]
            if not len(gc): continue
            # A marker hit can land on a paralogue that shares the gene name. Only accept the
            # CDS as a focal copy if its own product resolves to the gene of interest.
            r0 = gc.iloc[0]
            lab, _c = label_cat(r0["product"], int(r0["end"]) - int(r0["start"]) + 1)
            if lab == goi: out.append((contig, int(r0["start"])))
    else:  # carO by product
        gc = g[g["product"].str.contains(r"CarO|[Cc]arbapenem-associated resistance", case=True, na=False, regex=True)]
        for _, r in gc.iterrows(): out.append((r["contig"], int(r["start"])))
    # a copy can be reached from more than one hit row; keep each physical locus once,
    # then drop short contigs that merely duplicate a longer one in the same assembly
    return drop_redundant(sample, sorted(set(out)))

def context(goi, sample, contig, cstart):
    g = read_gff(sample)
    gcon = g[g["contig"] == contig].reset_index(drop=True)
    hits = gcon.index[gcon["start"] == cstart]
    if not len(hits): return None
    idx = hits[0]
    gc = gcon.iloc[[idx]]
    lo, hi = max(0, idx - WINDOW), min(len(gcon), idx + WINDOW + 1)
    # A window that runs off the end of the contig shows a truncation, not a biological
    # boundary. Tracks built from such windows must say so, otherwise a short track reads
    # as "this locus has lost its neighbours" when the assembly simply stopped.
    truncated = (lo == 0) or (hi == len(gcon))
    whole_contig = (lo == 0) and (hi == len(gcon))
    nb = gcon.iloc[lo:hi].copy().reset_index(drop=True)
    cstart, cstrand = gc.iloc[0]["start"], gc.iloc[0]["strand"]
    labs, cats = [], []
    for _, r in nb.iterrows():
        if r["start"] == cstart: labs.append(goi); cats.append("center")
        else:
            l, c = label_cat(r["product"], int(r["end"]) - int(r["start"]) + 1)
            labs.append(l); cats.append(c)
    nb["label"], nb["category"] = labs, cats
    # blaNDM-1 only. Bakta calls a 315 bp opposite-strand ORF overlapping the 3 prime end
    # of blaNDM-1 a "hypothetical protein". Translated on the blaNDM-1 strand it is an exact
    # fusion of the first 58 codons of Ble-MBL (WP_004201167.1, 121 aa) and the C-terminal
    # 40 codons of qacEdelta1 (ARO:3005010, 115 aa): the truncated Tn125 cargo welded to the
    # integron 3 prime conserved segment. Neither half reaches the 80 percent subject-coverage
    # threshold, so no marker database ever calls it. Present in 18 of 21 blaNDM-1 loci.
    if goi == "blaNDM-1":
        c0, c1, cs = gc.iloc[0]["start"], gc.iloc[0]["end"], gc.iloc[0]["strand"]
        m = ((nb["label"] == "HP") & (nb["end"] - nb["start"] + 1).between(300, 330)
             & (nb["start"] < c1 + 60) & (nb["end"] > c0 - 60) & (nb["strand"] != cs))
        nb.loc[m, "label"] = "\u0394bleMBL::\u0394qacE\u03941"
        nb.loc[m, "category"] = "ARG"
    cidx = nb.index[nb["start"] == cstart][0]
    cmid = (nb.loc[cidx, "start"] + nb.loc[cidx, "end"]) / 2
    nb["x0"] = nb["start"] - cmid; nb["x1"] = nb["end"] - cmid
    if cstrand == "-":  # flip so focal points forward
        nb["x0"], nb["x1"] = -nb["x1"], -nb["x0"]
        nb["strand"] = nb["strand"].map({"+": "-", "-": "+"})
    nb["orient"] = (nb["strand"] == "+").astype(int)
    nb = nb.sort_values("x0").reset_index(drop=True)
    # signature = ordered labels within +/- SIG_WINDOW genes of focal.
    # Locate the focal gene by category, not by label: a second copy of the same gene
    # inside the window carries the same label, and [0] would pick the wrong one
    # (head-to-tail plasmid dimers put both copies inside a +/- 14 CDS window).
    c2 = nb.index[nb["category"] == "center"][0]
    s_lo, s_hi = max(0, c2 - SIG_WINDOW), min(len(nb), c2 + SIG_WINDOW + 1)
    sig = tuple(nb["label"].iloc[s_lo:s_hi])
    # When the window covers an entire contig, the contig is very likely a small circular
    # plasmid and the point at which the assembler linearised it is arbitrary. Two assemblies
    # of the same molecule then differ only by rotation, and would otherwise be drawn as
    # separate context variants. Canonicalise such signatures to their smallest rotation so
    # one plasmid is one track.
    if whole_contig and len(sig) > 2:
        sig = min(tuple(sig[i:] + sig[:i]) for i in range(len(sig)))
    return nb, sig, contig, truncated

GOIS = ["blaOXA-72", "blaOXA-148", "blaOXA-23", "blaOXA-66", "blaNDM-1", "intI1", "carO", "oprD"]
frames = []
for goi in GOIS:
    if goi == "carO":
        samples = [s for s in sorted(st492) if has_gff(s)]
    else:
        samples = [s for s in mh.loc[mh["Gene"].astype(str) == goi, "Sample"].unique()
                   if s in st492 and has_gff(s)]
    # One entry per physical gene COPY, not per genome: a genome carrying two copies on
    # two contigs contributes two contexts and can fall into two different variants.
    groups = defaultdict(list); rep = {}; contig_of = {}; trunc_of = {}   # sig -> (rep_sample, nb); copy key -> contig
    n_copies = 0
    for s in samples:
        for contig, cstart in locate(goi, s):
            try:
                res = context(goi, s, contig, cstart)
            except Exception:
                res = None
            if res is None: continue
            nb, sig, contig, trunc = res
            key = (s, contig, cstart); n_copies += 1
            groups[sig].append(key); contig_of[key] = contig; trunc_of[key] = trunc
            if sig not in rep or len(nb) > len(rep[sig][1]):
                rep[sig] = (s, nb)                 # show the member with the most flanking CDS

    # mob_recon molecule type. A "chromosome" call really means "not binned to a plasmid";
    # on a short contig that is not a confident call, so require >= CHROM_MIN bp before
    # trusting it and otherwise report the location as unresolved.
    # A mob_recon "chromosome" call means only "not binned to a plasmid". Below 100 kb
    # that is its default bin for an unplaceable contig: such contigs are routinely
    # 99.7 to 100 percent identical to plasmid-binned contigs in sibling assemblies.
    CHROM_MIN = 100000
    def classify(m, sample=None):
        mt, cl, size, _rep = m
        # A plasmid call needs corroboration, exactly as a chromosome call does. Without it
        # a 2.5 kb contig carrying no replicon prints as a confident plasmid cluster while
        # the same contig called chromosome would be demoted to unresolved.
        if mt == "plasmid":
            if sample is None or cluster_has_replicon(sample, cl): return "plasmid", cl
            return "unresolved", None
        if mt == "chromosome" and size >= CHROM_MIN: return "chromosome", None
        return "unresolved", None
    def loc_counts(kset):
        locs = [classify(mob_call(k[0], contig_of[k]), k[0]) for k in kset]
        n_pl = sum(1 for l, _ in locs if l == "plasmid")
        n_ch = sum(1 for l, _ in locs if l == "chromosome")
        return locs, n_pl, n_ch, len(locs) - n_pl - n_ch

    # Rank by size, then force-include any group holding a copy whose location differs from
    # the dominant one for this locus. Ranking on size alone drops every rare context, and
    # the minority location is exactly the informative one: for a plasmid-borne gene the
    # chromosomal insertions, for an intrinsically chromosomal gene the plasmid escapes.
    # Without this the figure silently asserts the gene is always in one place.
    # Only a CONFIDENT minority forces inclusion: plasmid where the locus is normally
    # chromosomal, or chromosomal where it is normally plasmid-borne. "unresolved" means
    # only that the contig was too short to call, so it must not drag a track into the
    # figure; otherwise every fragmented assembly becomes its own variant.
    all_locs = [classify(mob_call(k[0], contig_of[k]), k[0])[0] for ks in groups.values() for k in ks]
    confident = [l for l in all_locs if l != "unresolved"]
    dominant = Counter(confident).most_common(1)[0][0] if confident else "plasmid"
    minority = "chromosome" if dominant == "plasmid" else "plasmid"
    by_size = sorted(groups.items(), key=lambda kv: -len(kv[1]))
    ranked = by_size[:MAX_VARIANTS]
    chosen = {sig for sig, _ in ranked}
    forced = 0
    for sig, kset in by_size[MAX_VARIANTS:]:
        if sig in chosen: continue
        if any(classify(mob_call(k[0], contig_of[k]), k[0])[0] == minority for k in kset):
            ranked.append((sig, kset)); chosen.add(sig); forced += 1

    # A track whose gene order is a contiguous run inside a longer track, and whose window
    # is cut by a contig end, is the same locus caught on a shorter contig rather than a
    # different arrangement. Mark it, so a reader does not read assembly fragmentation as
    # biological variation. Runs shorter than MIN_FRAG labels are not assigned a parent:
    # they are contained in almost everything and so identify nothing.
    MIN_FRAG = 3
    lab_of = {sig: list(rep[sig][1]["label"]) for sig, _ in ranked}
    def sublist(a, b):
        return len(a) < len(b) and any(b[i:i + len(a)] == a for i in range(len(b) - len(a) + 1))
    def mol_of(kset):
        locs, n_pl, n_ch, n_un = loc_counts(kset)
        if n_pl >= max(n_ch, 1) and n_pl:
            cl = Counter(c for l, c in locs if l == "plasmid" and c and c != "-")
            return ("plasmid", cl.most_common(1)[0][0] if cl else None)
        return ("chromosome" if n_ch > n_pl else "unresolved", None)
    mol_v = {i: mol_of(k) for i, (_, k) in enumerate(ranked, 1)}
    def compatible(a, b):
        """Same molecule, treating "unresolved" as compatible with anything.
        Unresolved records that the binning could not decide, not that the molecule
        differs, so refusing to link it to a resolved track hides real contig splits."""
        if a[0] == "unresolved" or b[0] == "unresolved": return True
        return a == b
    parent = {}; near_only = set()  # near_only retained for the label gate below
    for i, (sig_i, k_i) in enumerate(ranked, 1):
        li = lab_of[sig_i]
        if len(li) < MIN_FRAG: continue
        # Same molecule and cluster only. A free plasmid is trivially contained in its own
        # chromosomally integrated copy, and that is a biological difference, not a split.
        hits = [j for j, (sig_j, k_j) in enumerate(ranked, 1)
                if j != i and compatible(mol_v[i], mol_v[j]) and sublist(li, lab_of[sig_j])]
        # A near-identical rule used to link tracks here as well, but it joined tracks that
        # differ in gene CONTENT (one CDS called as two over the same span) rather than being
        # one locus seen at different contig lengths. Those are separate arrangements and are
        # left as separate tracks. True rotations of a circular plasmid are already merged
        # upstream, by canonicalising whole-contig signatures.
        if hits: parent[i] = hits

    for rank, (sig, kset) in enumerate(ranked, 1):
        if len(kset) < MIN_N: continue
        genomes = sorted({k[0] for k in kset})
        cc = Counter(ctry.get(g, "NA") for g in genomes)
        ctxt = ", ".join(f"{c} {n}" for c, n in cc.most_common(3))
        locs, n_pl, n_ch, n_un = loc_counts(kset)
        resolved = n_pl + n_ch
        if resolved == 0 or n_un > resolved:
            molecule, cluster, loc = "unresolved", "", "location unresolved"
        elif n_pl >= n_ch:
            cl = Counter(c for l, c in locs if l == "plasmid" and c and c != "-")
            cluster = cl.most_common(1)[0][0] if cl else None
            molecule = "plasmid"
            # show the split when the cluster is not unanimous, rather than printing the
            # majority cluster as though every member carried it
            if cluster and len(cl) > 1:
                loc = f"plasmid {cluster} ({cl[cluster]}/{n_pl})"
            else:
                loc = f"plasmid {cluster}" if cluster else "plasmid"
        else:
            molecule, cluster, loc = "chromosome", "", "chromosome"
        if n_ch and molecule == "plasmid": loc += f", {n_ch} chromosomal"
        n_lab = f"n = {len(kset)} copies" + (f" in {len(genomes)} genomes" if len(genomes) != len(kset) else "")
        n_tr = sum(1 for k in kset if trunc_of.get(k))
        # a sub-run link only means "split across contigs" when the window is truncated;
        # a near-identical link is an annotation difference and needs no truncation
        par = parent.get(rank) if (n_tr or rank in near_only) else None
        # truncation is still computed above for the fragment-bracket logic (par),
        # but is no longer printed in the track label
        track = f"variant {rank}  ·  {loc}  ·  {n_lab}  ({ctxt})"
        rep_s, nb = rep[sig]; nb = nb.copy()
        nb["goi"] = goi; nb["variant"] = rank; nb["n"] = len(kset)
        nb["track"] = track; nb["rep_sample"] = rep_s
        nb["molecule"] = molecule; nb["plasmid_cluster"] = cluster or ""
        # variant number this track is a contig-split fragment of, for the figure to bracket
        nb["fragment_of"] = par[0] if par else 0
        frames.append(nb[["goi", "variant", "n", "track", "rep_sample", "molecule",
                          "plasmid_cluster", "fragment_of", "x0", "x1", "orient", "label", "category"]])
    print(f"{goi:11s} {len(samples):3d} genomes / {n_copies:3d} copies -> {len(groups):2d} contexts; "
          f"showing {len(ranked):2d}  (dominant location: {dominant}; {forced} forced in as minority-location)")

allg = pd.concat(frames, ignore_index=True)
allg.to_csv(OUT, index=False)
print(f"\nwrote {OUT}  ({allg['goi'].nunique()} loci, {allg.groupby('goi')['variant'].nunique().sum()} variant tracks)")
