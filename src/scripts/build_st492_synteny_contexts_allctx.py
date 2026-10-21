#!/usr/bin/env python3
"""Build the synteny context tables with EVERY context kept, not just the five most common.

The shipped builder, src/scripts/build_st492_synteny_contexts.py, caps each locus at
MAX_VARIANTS = 5 so the per-locus figures stay legible. The all-contexts figures need the
full set, so this wrapper reads that builder's source, patches it in memory and executes it.
The builder file itself is never modified and stays byte-identical, which matters because
the capped figures in the manuscript are generated from it. Every correction below therefore
lives here, as an anchored in-memory patch that must apply exactly once or the run stops.

Rewritten in memory:
  MAX_VARIANTS = 5    ->  no cap
  OUT                 ->  st492_synteny_variants_allctx8.csv (the capped CSV survives)
  corrections found by the 2026-10 audit of the blaOXA-72 all-contexts figure:
    * two doubled (tandem) AC163 plasmid contigs get a track of their own that draws the
      whole contig, and both blaOXA-72 copies on each are counted (named contigs only, no
      heuristic; user decision of 2026-10-06)
    * whole-contig circular copies whose sequences are identical after trimming the
      assembler's terminal overlap share one signature, so a different opening point or a
      CDS cut by the contig end no longer creates a separate context
    * one sub-100 kb MOB-suite chromosome call is accepted where core chromosomal genes
      flank the gene (named contig only)
    * blaOXA-72 contigs that MOB-suite co-binned into AB082 with a rep_cluster_734 plasmid,
      but which carry the AC163 replicon rep_cluster_1172 themselves, are counted as AC163
    * country breakdowns list every country when a track has at most seven, otherwise the
      top three plus "+K other", ordered by count then name and never cutting inside a tie
    * track labels say "n = 1 copy", not "n = 1 copies"
    * display labels (never the grouping signature) are corrected: TonB-dependent receptor
      fragments are not "tonB", no undefined " DCP" suffix, EAL/GGDEF named properly,
      NCBIfam accessions unlabelled, and for blaOXA-72 the Bakta gene name is used where the
      product itself names that gene (but a domain abbreviation such as RRM stays a domain)
    * a per-copy provenance table is written, which the builder never did
    * gene-level rules for every locus (user decisions of 2026-10-06, see _st492_display_labels):
      the bleMBL::qacEdelta1 read-through ORF is computed from the assembly and named only when
      both halves check out; a grey CDS with an AMR call is an ARG; IS transposases print by
      family unless an element is named; aph(3') carries its subclass; Bakta gene names replace
      placeholders; a CDS cut by a contig end gets a prime
    * named blaNDM-1 corrections read from the assembly: the contig-end msr(E)' fragment and
      sul1' and ISCR1' broken by ambiguity codes; a gene split by an N scaffold gap (the
      qacEdelta1' of the Serbian AC237 assemblies) is found generically and drawn with a gap mark
    * every arrow carries a note (why its name has a prime, a delta or "fam.") and its best
      CGE-MGE / ISFinder match (mobile), every track its scaffold gaps and contig ends, and the
      copies dropped as duplicate contigs are listed
    * fragment_of never names a parent on another molecule or for an unresolved track (the
      plot already refused such links), and is set only when exactly one parent qualifies
    * a whole circle whose contig ends repeat a stretch inexactly is recognised as a circle
  and by the 2026-10-06 intI1 fix round (user decisions of 2026-10-06):
    * two distinct intI1 copies on one chromosome (GCA_964331735.1 contig_1) are a named
      decision: both counted, each its own track marked "same genome as variant N", the
      build stops if they move; any other undecided several-copies contig is reported
    * contigs byte-identical to a plasmid cluster's contig elsewhere count as that cluster
      (named, sequence-checked: two intI1 AC715 contigs are AC237)
    * a cluster is printed only on a strict majority of a track's copies, never as a fraction;
      otherwise the track is U
    * IS elements named by the partner CDS of a two-gene element (IS1326), resolvases printed
      tnpR, a Tn21 tnpM called by CGE-MGE drawn as MGE, intI1 printed only for a class 1
      integrase, a CDS cut by a contig end recognised also by its sequence lying inside a
      longer gene drawn at the locus
  and by the 2026-10-07 stacking round (user decisions of 2026-10-07):
    * circle facts are per contig, computed once, and the same under every locus (columns
      end_overlap_bp, end_repeat_bp, seamless_circle, closed_circle, circle_kind); eight named
      9,830 bp blaOXA-72 contigs are closed circles whose ends join without overlap
      (_ST492_SEAMLESS); truncated / circle_opening / whole_contig stay window facts
    * the bleMBL::qacEdelta1 ORF of a genome with no qacEdelta1 of its own is checked against the
      CARD qacEdelta1 reference; a Bakta-named IS fragment gets a prime (IS15'); the frameshifted
      terminal-repeat IS15 copy of GCA_031156335.1 contig_80 is IS15' (named)
    * "DNA-binding ..., X family" regulators print "X fam." on every locus; the ambiguity-broken
      ISCR1 keeps its name (IS91/ISCR1'); a second copy of the gene of interest is category
      goi_copy; a CDS with only a resistance point-mutation call stays grey with a note
  and by its repair round (2026-10-07):
    * an arrow inside a contig's inexact terminal repeat is noted "terminal repeat" whatever
      the window reaches (one contig, one note under every locus)
    * a remnant renamed from its sequence says so in its note (share of 21-mers, reference
      genome, Bakta's product), and an ISFinder-named or primed IS records its CDS length
      against the element's transposase, so the caption can state the prime rule
  and by the 2026-10-06 repair round after the merge:
    * an arrow read from the assembly at a contig end is drawn only when it follows the
      window's outermost CDS (intI1 dropped msr(E)' beyond an undrawn ISEc29)
    * a window that stops at the opening point of a closed circle is not "truncated" (new
      per-copy column circle_opening; contig_end "R, circle opening")
    * a domain-only product prints Bakta's gene symbol (ycbK, hdeD, rosB), and blaOXA-72 no
      longer prints "Lrp fam." / "AcrR fam." where Bakta gives lrp / acrR
    * the family-only IS note gives ISFinder's coverage of the CDS, not of the element
    * the same-genome mark reads "same genome as variant N" (CSV number), the focal copy's
      CGE-MGE match is recorded too, and per copy every track holding its gene order as a run
      (run_tracks)

Outputs (in /Data/ST492/figures2026, or --outdir):
  st492_synteny_variants_allctx8.csv   one row per gene arrow, per variant track, with the
                                       representative CDS coordinates (contig, start, end)
  st492_synteny_contexts_allctx8.csv   one row per gene copy, with its variant assignment,
                                       contig facts, MOB-suite call and country
  st492_synteny_dropped_allctx8.csv    one row per hit dropped as a duplicate contig, with the
                                       contig holding it, the 31-mer containment and the
                                       builder's own probe share

Usage:  python3 src/scripts/build_st492_synteny_contexts_allctx.py [--outdir DIR]
"""
import argparse
import pathlib
import re
import runpy
import sys
import tempfile

SRC = pathlib.Path(__file__).with_name("build_st492_synteny_contexts.py")
SUFFIX = "_allctx8"

# Injected into the builder just before its main loop, so it can wrap the functions defined
# above that point. Nothing in it changes how a context is grouped except where a comment
# says so and cites the evidence.
PRELUDE = r'''
# ==== allctx corrections, injected by build_st492_synteny_contexts_allctx.py ===========
import re as _re
import sys as _sys
from urllib.parse import unquote as _unquote

_ST492_COPYINFO = {}    # (goi, sample, contig, cstart) -> per-copy facts for the dump
_ST492_COPYROWS = []

# Keep Bakta's gene= attribute next to the product. Used for display labels only, never for
# the grouping signature, so it cannot merge or split a context.
_read_gff_base = read_gff
def read_gff(sample):
    d = _read_gff_base(sample)
    if "gene" not in d.columns:
        gene = {}
        with open(f"{BK}/{sample}/{sample}.gff3") as fh:
            for line in fh:
                if line.startswith("#") or "\tCDS\t" not in line: continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9: continue
                a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
                gene[(f[0], int(f[3]), int(f[4]), f[6])] = a.get("gene", "")
        d["gene"] = [gene.get(k, "") for k in zip(d["contig"], d["start"], d["end"], d["strand"])]
    return d

# Two doubled (tandem) AC163 plasmid contigs, both from Swiss genomes. Each carries blaOXA-72
# twice, exactly 9,830 bp apart (the modal AC163 unit length): GCA_964331335.1 contig_5
# (19,150 bp, copies at 5,885 and 15,715) and GCA_964331905.1 contig_5 (19,266 bp, at 6,010 and
# 15,840). blastn of each contig against itself puts its second half on its first:
#   GCA_964331905.1: 9831-19170 = 1-9340 at 100.000 percent; 31-mer containment of the second
#                    period in the first 0.997
#   GCA_964331335.1: 9829-19065 = 1-9234 at 99.47 percent (its differences lie mostly in the
#                    first 3 kb of the contig); 31-mer containment 0.973
# Both are shorter than 2 x 9,830 bp. 18 other genomes in the same gene-order context carry
# the plasmid as a single 9.8 to 9.9 kb contig (16 MOB-suite AC163 calls and 2 AB082 co-bins,
# see _st492_classify_note), and 1 more such contig differs from it by one substitution.
# Whether these two contigs are a real tandem dimer or the monomer assembled about 1.96 times
# over cannot be told from the assemblies: a perfect tandem repeat gives the same assembly
# graph either way, and only read depth or long reads could settle it. By the user's decision
# (2026-10-06) both copies on each contig are counted, and the two contigs share one track of
# their own that draws the whole doubled contig, so the second copy is visible. The figure
# says the contig holds two copies 9,830 bp apart and that an assembly artefact cannot be
# ruled out; it never calls it a proven dimer. Left to the base builder, each copy's +/- 9 CDS
# window would cut the doubled contig at a different place, giving four copies three or four
# different signatures, so the track is pinned by a fixed signature instead. Deliberately two
# named contigs and no length or spacing heuristic: a general rule would also catch genuine
# same-contig duplications elsewhere.
#
# Not to be confused with the second blaOXA-72 hit that marker_hits reports in each of the
# three Zurich genomes GCA_031875465.1, GCA_031875745.1 and GCA_031876255.1. Those hits sit on
# short contigs (contig_140 2,895 bp, contig_88 2,894 bp, contig_84 2,895 bp) that, apart from
# a 40 to 41 bp overhang, are the reverse complement of the terminal about 2.85 kb of the
# AC163 contig in the same genome (contig_60 6108-8961, contig_51 6110-8962, contig_47
# 6111-8965). They are the same physical locus caught twice, and the base builder's
# drop_redundant already removes them, so each of these genomes contributes one copy. That
# stays so: drop_redundant only compares DIFFERENT contigs, so it never touches the two
# doubled contigs above.
_ST492_TANDEM = {("blaOXA-72", "GCA_964331335.1", "contig_5"),
                 ("blaOXA-72", "GCA_964331905.1", "contig_5")}
_ST492_TANDEM_PERIOD = 9830
_ST492_TANDEM_SPAN = {}     # (goi, sample, contig) -> (start of first copy, start of second)
# One fixed signature for all four copies, so they form one track and nothing else joins it
_ST492_TANDEM_SIG = ("blaOXA-72 x2 on one contig, 9,830 bp apart",)

# Copies drop_redundant removes as the same locus caught on a shorter contig, recorded so the
# caption can list them (guide 2.2) and the verifier can recount them. blaNDM-1: 21 marker-hit
# loci, 4 dropped, 17 counted. GCA_031875465.1 contig_149, GCA_031876255.1 contig_90 and
# GCA_052150845.1 contig_99 (1,049 bp each) are, apart from a 41 bp overhang, the reverse
# complement of the first 1,008 bp of the blaNDM-1 contig of the same genome (contig_76,
# contig_56, contig_62); GCA_031875745.1 contig_92 (1,974 bp) lies wholly inside contig_43.
# Counted once, like the three Zurich blaOXA-72 fragments below. intI1: 9 dropped, 4 near-full
# duplicates (at least 90 percent of the length of the contig holding them; the user decided on
# 2026-10-06 that they stay counted once) and 5 shorter fragments.
_ST492_DROPPED = []
_ST492_LAST_DROP = []
_drop_redundant_base = drop_redundant
def drop_redundant(sample, copies):
    keep = _drop_redundant_base(sample, copies)
    _ST492_LAST_DROP[:] = [c for c in copies if c not in keep]
    return keep

def _st492_kmer_contained(a, b, k=31):
    """Share of a's 31-mers (every 7th) found in b on either strand, as the verifier counts."""
    ka = {a[i:i + k] for i in range(0, len(a) - k + 1, 7)}
    if not ka: return 0.0
    rb = _rc(b)
    return sum(1 for x in ka if x in b or x in rb) / len(ka)

def _st492_probe_frac(sh, lg):
    """The builder's own containment measure, copied: the share of 200 bp probes, taken every
    400 bp across the middle 80 percent of the short contig, found in the long one on either
    strand. The builder drops the short contig at 0.9; the caption prints the share."""
    core = sh[len(sh) // 10: len(sh) - len(sh) // 10]
    probes = [q for q in (core[i:i + 200] for i in range(0, max(1, len(core) - 200), 400))
              if len(q) == 200]
    if not probes: return 0.0
    return sum(1 for q in probes if q in lg or _rc(q) in lg) / len(probes)

# Two DISTINCT copies on one contig, each drawn as its own context (user decision of
# 2026-10-06). GCA_964331735.1 contig_1 is a 3.96 Mb chromosome (MOB-suite chromosome) with
# two complete class 1 integrons 679,550 bp apart on opposite strands: intI1 at 2,160,506 (-)
# and 2,840,056 (+). The two intI1 CDS differ at 3 positions and their flanks share only about
# 1.35 kb before diverging (aadA1/aac(3)-Ia, Tn21 and Tn1331 on one side; armA, ISCR1 and the
# dfrA12 array with dusA/acnA on the other), so this is not a tandem and the whole-contig rule
# of _ST492_TANDEM cannot apply to a chromosome. Both copies are counted (63 copies in 62
# genomes) and each keeps the standard +/- WINDOW CDS track; the figure marks the two tracks
# as the same genome. Named contig and exact copy positions, so the build stops if they move.
_ST492_SEPARATE = {("intI1", "GCA_964331735.1", "contig_1"): ((2160506, "-"), (2840056, "+"))}

_locate_base = locate
def locate(goi, sample):
    _ST492_LAST_DROP[:] = []
    out = _locate_base(goi, sample)
    seqs = read_fna(sample)
    for contig, cstart in _ST492_LAST_DROP:
        best = max(((_st492_kmer_contained(seqs[contig], seqs[c]), c) for c, _s in out if c != contig),
                   default=(0.0, ""))
        sh, lg = seqs[contig], seqs.get(best[1], "")
        _ST492_DROPPED.append(dict(goi=goi, sample=sample, contig=contig, cstart=cstart,
                                   contig_len=len(sh), in_contig=best[1], in_contig_len=len(lg),
                                   containment=round(best[0], 3),
                                   probe_share=round(_st492_probe_frac(sh, lg), 3) if lg else 0.0,
                                   kind=("near-full duplicate" if lg and len(sh) >= 0.9 * len(lg)
                                         else "fragment")))
    g = read_gff(sample)
    for (dg, ds, dc), want in _ST492_SEPARATE.items():
        if (dg, ds) != (goi, sample): continue
        got = tuple(sorted((c[1], g[(g["contig"] == dc) & (g["start"] == c[1])]["strand"].iloc[0])
                           for c in out if c[0] == dc))
        if got != tuple(sorted(want)):
            _sys.exit(f"two-copy decision for {goi} {sample} {dc} no longer matches the data: "
                      f"copies at {got}, expected {sorted(want)}")
    # Any OTHER contig with several copies of the gene has no recorded decision. Reported, not
    # fatal, so a locus nobody has audited yet still builds; the verifier lists them too.
    for contig, k in Counter(c[0] for c in out).items():
        if (k > 1 and (goi, sample, contig) not in _ST492_TANDEM
                and (goi, sample, contig) not in _ST492_SEPARATE):
            print(f"  [allctx] UNDECIDED: {goi} {sample} {contig} carries {k} copies")
    tand = [c for c in out if (goi, sample, c[0]) in _ST492_TANDEM]
    if not tand: return out
    starts = sorted(c[1] for c in tand)
    contig = tand[0][0]
    if len(starts) != 2 or starts[1] - starts[0] != _ST492_TANDEM_PERIOD:
        _sys.exit(f"tandem exception for {sample} {contig} no longer matches the data: "
                  f"copies at {starts}, expected two {_ST492_TANDEM_PERIOD} bp apart")
    _ST492_TANDEM_SPAN[(goi, sample, contig)] = (starts[0], starts[1])
    return out                                      # both copies kept and counted

def _st492_whole(sample, contig, cstart):
    gcon = read_gff(sample); gcon = gcon[gcon["contig"] == contig].reset_index(drop=True)
    idx = gcon.index[gcon["start"] == cstart]
    if not len(idx): return False
    i = idx[0]
    return (i - WINDOW <= 0) and (i + WINDOW + 1 >= len(gcon))

_context_base = context
def context(goi, sample, contig, cstart):
    span = _ST492_TANDEM_SPAN.get((goi, sample, contig))
    if span is None:
        res = _context_base(goi, sample, contig, cstart)
        if res is not None:
            _ST492_COPYINFO[(goi, sample, contig, cstart)] = dict(
                whole_contig=_st492_whole(sample, contig, cstart), tandem=False,
                ncds_window=len(res[0]), sig=res[1], seq_merged_from="")
        return res
    # Both copies of a doubled contig return the SAME drawing: the whole contig, aligned on its
    # first copy, with the second copy drawn as a further blaOXA-72 arrow 9,830 bp away. The
    # first copy is the anchor because, both copies being on the minus strand, that puts the
    # contig at -12.0 to +6.3 kb, inside the x-window the other tracks already span; anchored on
    # the second copy it would run to +16.1 kb and widen the figure. The base window of
    # +/- 14 CDS covers all 21 CDS of either contig from the first copy, which is checked.
    s1, _s2 = span
    if not _st492_whole(sample, contig, s1):
        _sys.exit(f"tandem exception for {sample} {contig}: the +/- {WINDOW} CDS window from "
                  f"the first copy no longer covers the whole contig")
    res = _context_base(goi, sample, contig, s1)
    if res is None: return None
    nb, _sig, _c, trunc = res
    _ST492_COPYINFO[(goi, sample, contig, cstart)] = dict(
        whole_contig=True, tandem=True, ncds_window=len(nb), sig=_ST492_TANDEM_SIG,
        seq_merged_from="")
    return nb, _ST492_TANDEM_SIG, contig, trunc

# Whole-contig circular copies, compared by SEQUENCE before grouping. Rotation
# canonicalisation of the label ring only merges identical label rings, so the same plasmid
# opened inside a CDS (Bakta then calls the remnant something else, e.g. "DNA-binding protein"
# for a cut "DNA replication protein") became a separate context. Here a contig counts as a
# closed circle only when its two ends overlap exactly (the assembler's terminal repeat, 77 bp
# for the 9.8 to 9.9 kb AC163 contigs); the overlap is trimmed and the circle is anchored at
# the focal gene's start on its own strand. Copies with byte-identical anchored circles are
# the same molecule and are given one signature: the one held by the most copies at this
# locus. No identity threshold: a single substitution keeps a copy separate.
def _st492_overlap(q, kmax=1000, kmin=20):
    for k in range(min(kmax, len(q) // 2), kmin - 1, -1):
        if q[:k] == q[-k:]: return k
    return 0

# A whole plasmid can also come as a contig whose two ends repeat the same stretch INEXACTLY,
# so the exact test above misses it. blaNDM-1 csv v6 (GCA_031156335.1 contig_80, 22,849 bp) is
# the 22,080 bp AC237 circle of csv v7 and v8 plus 767 bp of it again: positions 1-766 and
# 22083-22849 align with a single 1 bp indel (blastn 99.87 percent), and the rest aligns to
# csv v7 end to end. Used only to describe such a contig (whole circle, not cut by the
# assembly), never to merge or regroup it.
def _st492_edits(a, b, band=8):
    """Edits aligning all of a to the best-matching prefix of b, within a diagonal band."""
    prev = {j: j for j in range(min(band, len(b)) + 1)}
    for i in range(1, len(a) + 1):
        cur = {}
        for j in range(max(0, i - band), min(len(b), i + band) + 1):
            best = prev.get(j, 10**9) + 1
            if j: best = min(best, cur.get(j - 1, 10**9) + 1,
                             prev.get(j - 1, 10**9) + (a[i - 1] != b[j - 1]))
            cur[j] = best
        prev = cur
    return min(prev.values())

def _st492_near_overlap(q, kmax=2000, seed=30):
    """(length, edits) of an inexact terminal repeat of at least 100 bp with at most 1 percent
    edits; (0, 0) when the contig has none or its overlap is exact."""
    L = len(q)
    if L < 4 * seed or _st492_overlap(q): return (0, 0)
    lo = max(L - kmax, L // 2)
    for m in _re.finditer(_re.escape(q[:seed]), q[lo:]):
        tail = q[lo + m.start():]
        if len(tail) < 100: continue
        e = _st492_edits(tail, q[:len(tail) + 8])
        if e <= max(1, 0.01 * len(tail)): return (len(tail), e)
    return (0, 0)

# Closed circles whose ends join WITHOUT any overlap: the assembler already trimmed its terminal
# repeat, so neither test above can see them. Named, not detected (user preference 2026-10-07):
# the eight 9,830 bp blaOXA-72 AC163 contigs below are each, on one strand or the other, a
# rotation of the 9,907 bp circle GCA_029679805.1 contig_57 with its 77 bp end overlap trimmed
# (0 differences), checked on every build. Used to describe a contig (closed circle, not cut by
# the assembly), never to merge or regroup it, so the tracks and their numbers do not move.
_ST492_SEAMLESS_REF = ("GCA_029679805.1", "contig_57")
_ST492_SEAMLESS = {("GCA_016474335.1", "contig_35"), ("GCA_016474395.1", "contig_35"),
                   ("GCA_022017315.1", "contig_40"), ("GCA_022017345.1", "contig_39"),
                   ("GCA_028577815.3", "contig_47"), ("GCA_048215925.1", "contig_42"),
                   ("GCA_048245695.1", "contig_39"), ("GCA_048245935.1", "contig_44")}

# Circle facts are properties of the CONTIG, so they are computed once per contig and are the
# same under every locus that draws it (blaNDM-1 and intI1 used to disagree about the same three
# AC237 circles because the old fields were only filled when the window covered the contig).
#   overlap: exact end overlap (bp); repeat, edits: inexact terminal repeat; seamless: named
#   above; closed: any of the three; kind: the words the caption uses for it
_ST492_CIRC = {}
def _st492_contig_circle(sample, contig):
    if (sample, contig) in _ST492_CIRC: return _ST492_CIRC[(sample, contig)]
    q = read_fna(sample).get(contig, "")
    ov = _st492_overlap(q) if q else 0
    rep, ed = _st492_near_overlap(q) if q and not ov else (0, 0)
    seam = (sample, contig) in _ST492_SEAMLESS
    if seam:
        r = read_fna(_ST492_SEAMLESS_REF[0])[_ST492_SEAMLESS_REF[1]]
        C = r[:-_st492_overlap(r)]
        if not (len(q) == len(C) and (q in C + C or _rc(q) in C + C)):
            _sys.exit(f"named seamless circle {sample} {contig} is no longer a rotation of "
                      f"{' '.join(_ST492_SEAMLESS_REF)}")
    kind = (f"its two ends overlap by {ov} bp" if ov else
            f"its two ends repeat {rep:,} bp, {ed} edit{'' if ed == 1 else 's'} apart" if rep else
            f"its ends join without overlap: a rotation of the closed circle {' '.join(_ST492_SEAMLESS_REF)}"
            if seam else "")
    _ST492_CIRC[(sample, contig)] = dict(overlap=ov, repeat=rep, edits=ed, seamless=seam,
                                         closed=bool(ov or rep or seam), kind=kind)
    return _ST492_CIRC[(sample, contig)]

def _st492_circle(goi, sample, contig, cstart):
    info = _ST492_COPYINFO.get((goi, sample, contig, cstart), {})
    q = read_fna(sample).get(contig)
    if not q: return None
    g = read_gff(sample)
    row = g[(g["contig"] == contig) & (g["start"] == cstart)]
    if not len(row): return None
    a, b, strand = int(row.iloc[0]["start"]) - 1, int(row.iloc[0]["end"]) - 1, row.iloc[0]["strand"]
    # A doubled contig keeps its own track: one period of it can be byte-identical to a
    # monomer circle, and merging on that would put its copies back into the monomer track.
    if (goi, sample, contig) in _ST492_TANDEM: return None
    if not info.get("whole_contig"): return None
    k = _st492_overlap(q)
    if not k: return None
    C = q[:-k]
    if not (0 <= a <= b < len(C)): return None
    if strand == "+": return C[a:] + C[:a]
    return _rc(C[b + 1:] + C[:b + 1])

def _st492_regroup(goi, groups, rep):
    cls = defaultdict(list)
    for sig, ks in groups.items():
        for k in ks:
            c = _st492_circle(goi, *k)
            if c: cls[c].append((sig, k))
    moves = []
    for c, members in cls.items():
        sigs = sorted({s for s, _ in members})
        if len(sigs) < 2: continue
        tgt = max(sigs, key=lambda s: len(groups[s]))
        moves += [(k, s, tgt) for s, k in members if s != tgt]
    for k, s, t in moves:
        groups[s].remove(k); groups[t].append(k)
        _ST492_COPYINFO[(goi,) + k]["seq_merged_from"] = " ".join(s)
    for s in [s for s, ks in groups.items() if not ks]:
        del groups[s]; rep.pop(s, None)
    if moves:
        print(f"  [allctx] {goi}: {len(moves)} whole-contig copies joined a sequence-identical "
              f"circular context ({len({m[1] for m in moves})} signatures emptied)")
    return groups, rep

# Location calls. Two narrow, documented exceptions to the builder's classify():
#  * GCA_021570155.1 contig_44 (46,943 bp) is called chromosome by MOB-suite but demoted to
#    unresolved by the 100 kb floor. It carries an ISAba31-bracketed AC163 unit flanked on
#    both sides by core chromosomal genes (mutS, lrp, lpxL, fdxA; acrR, araJ), so the call is
#    accepted for this contig. Named, not a rule: other short chromosome calls stay unresolved.
#  * MOB-suite put five blaOXA-72 contigs into AB082 together with a co-resident
#    rep_cluster_734 plasmid, although each carries rep_cluster_1172 (the AC163 replicon)
#    itself and is fully contained in an AC163 contig. For blaOXA-72 only, an AB082 call on a
#    contig whose own replicon is rep_cluster_1172, in a genome whose AB082 bin also holds a
#    rep_cluster_734 contig, is counted as AC163.
#  * a contig whose sequence is byte-identical to a contig that MOB-suite put in a plasmid
#    cluster in another genome is that cluster, whatever bin ID its own genome got. Named
#    contigs only, each checked against its reference contig's sequence on every build:
#      intI1: GCA_022017315.1 contig_62 and GCA_022017345.1 contig_59 (15,062 bp, rep_cluster_1254)
#      were binned AC715, but are identical to GCA_016474395.1 contig_59 (AC237; MOB-suite md5
#      64263cfaa5f436abd4621515369b7c94 for all three) and their bins hold the same companion
#      contigs as AC237 bins elsewhere: a cluster-ID split, not a different plasmid.
_ST492_CHROM_CONFIRMED = {("GCA_021570155.1", "contig_44")}
_ST492_SAME_AS_CLUSTER = {
    ("GCA_022017315.1", "contig_62"): ("AC237", "GCA_016474395.1", "contig_59"),
    ("GCA_022017345.1", "contig_59"): ("AC237", "GCA_016474395.1", "contig_59"),
}

def _st492_same_as(sample, contig):
    ent = _ST492_SAME_AS_CLUSTER.get((sample, contig))
    if not ent: return None
    cl, rs, rc_ = ent
    if (read_fna(sample).get(contig) != read_fna(rs).get(rc_)
            or mob_call(rs, rc_)[:2] != ("plasmid", cl)):
        _sys.exit(f"named cluster exception {sample} {contig} -> {cl}: it is no longer identical "
                  f"to {rs} {rc_}, or that contig is no longer in {cl}")
    return cl, rs, rc_

def _st492_cobin734(sample):
    mob_call(sample, None)
    return any(cl == "AB082" and "rep_cluster_734" in str(rt)
               for _mt, cl, _sz, rt in _mob.get(sample, {}).values())

def _st492_classify_note(goi, sample, contig):
    m = mob_call(sample, contig)
    if (sample, contig) in _ST492_CHROM_CONFIRMED and m[0] == "chromosome":
        return ("chromosome", None), "MOB-suite chromosome call accepted under 100 kb"
    same = _st492_same_as(sample, contig)
    if same:
        return (("plasmid", same[0]),
                f"{m[1]} bin; identical to {same[1]} {same[2]} ({same[0]}); counted as {same[0]}")
    if (goi == "blaOXA-72" and m[0] == "plasmid" and m[1] == "AB082"
            and str(m[3]).strip() == "rep_cluster_1172" and _st492_cobin734(sample)):
        return ("plasmid", "AC163"), "AB082 co-bin with rep_cluster_734; counted as AC163"
    return None, ""

def _st492_classify(goi, sample, contig, classify):
    res, _note = _st492_classify_note(goi, sample, contig)
    return res if res else classify(mob_call(sample, contig), sample)

# Country breakdown: every country when there are at most _ST492_ALL_COUNTRIES_MAX (7),
# otherwise the top three and "+K other", K being the remaining genomes. The builder printed
# the top three only, so a label could silently account for 35 of 42 genomes. A missing
# country prints as "unknown". Countries are ordered by count, ties alphabetically, and the
# cut never splits a tie: every country with the same count as the third one is shown.
# Counter.most_common() broke ties by insertion order, so a 20-copy label printed "USA 3" and
# hid Slovenia 3 and Switzerland 3 in "+9 other", which made three of the seven Swiss
# blaOXA-72 genomes invisible. The limit was raised from four to seven by the user's decision
# of 2026-10-06, so a seven-country track (blaOXA-72, 18 copies) names its single Swiss genome
# instead of folding it into "+4 other".
_ST492_ALL_COUNTRIES_MAX = 7
def _st492_country_text(cc):
    items = sorted((("unknown" if c == "NA" else c, n) for c, n in cc.items()),
                   key=lambda t: (-t[1], t[0]))
    if len(items) <= _ST492_ALL_COUNTRIES_MAX:
        shown, rest = items, 0
    else:
        k = 3
        while k < len(items) and items[k][1] == items[2][1]:
            k += 1
        shown, rest = items[:k], sum(n for _, n in items[k:])
    txt = ", ".join(f"{c} {n}" for c, n in shown)
    return txt + (f", +{rest} other" if rest else "")

# "n = 1 copies" was written into the track column for every single-copy track and only
# patched over in the drawing. Pluralised here, so the table says what the figure says.
def _st492_n_lab(n_copies, n_genomes):
    lab = f"n = {n_copies} cop{'y' if n_copies == 1 else 'ies'}"
    if n_genomes != n_copies:
        lab += f" in {n_genomes} genome{'' if n_genomes == 1 else 's'}"
    return lab

# Track location text (replaces the builder's block). A plasmid cluster is printed only when it
# holds a strict majority of ALL the track's copies, and then without a fraction; otherwise the
# track is lettered U (user decision of 2026-10-06, every locus). The builder printed the most
# common cluster whenever plasmid calls were not outnumbered by unresolved ones, so a 2-of-4 tie
# read "P AB082", and it printed fractions such as "AC237 (21/23)" whose denominator was not the
# track's n. The chromosome branch is the builder's, unchanged: resolved calls must not be
# outnumbered, and chromosome must outnumber plasmid.
def _st492_location(locs, n):
    n_pl = sum(1 for l, _ in locs if l == "plasmid")
    n_ch = sum(1 for l, _ in locs if l == "chromosome")
    n_un = len(locs) - n_pl - n_ch
    if n_pl + n_ch == 0 or n_un > n_pl + n_ch:
        return "unresolved", "", "location unresolved"
    if n_pl < n_ch:
        return "chromosome", "", "chromosome"
    cl = Counter(c for l, c in locs if l == "plasmid" and c and c != "-")
    top = sorted(cl.items(), key=lambda t: (-t[1], t[0]))[0] if cl else None
    if not top or 2 * top[1] <= n:
        return "unresolved", "", "location unresolved"
    loc = f"plasmid {top[0]}"
    if n_ch: loc += f", {n_ch} chromosomal"
    return "plasmid", top[0], loc

def _st492_dump(goi, rank, kset, contig_of, trunc_of, parents, classify):
    # every qualifying parent, also when the link was refused for having more than one
    parents = _ST492_CANDS.get((goi, rank), parents)
    for k in kset:
        s, contig, cstart = k
        m = mob_call(s, contig)
        res, note = _st492_classify_note(goi, s, contig)
        loc, cl = res if res else classify(m, s)
        info = _ST492_COPYINFO.get((goi,) + k, {})
        q = read_fna(s).get(contig, "")
        # contig facts (the same under every locus); a doubled contig is not called a circle
        cf = _st492_contig_circle(s, contig)
        closed = cf["closed"] and not info.get("tandem")
        whole = bool(info.get("whole_contig"))
        _ST492_COPYROWS.append(dict(
            goi=goi, variant=rank, sample=s, contig=contig, cstart=cstart,
            contig_len=len(q),
            # WINDOW facts (they depend on the locus): the builder flags every window that
            # reaches a contig end; a closed circle's ends are only where the assembler opened
            # it, so a window there is not truncated. circle_opening marks a window that stops
            # at a circle's opening point without drawing the circle whole (intI1 csv v15 to
            # v17); a circle drawn whole has whole_contig and closed_circle both True.
            truncated=bool(trunc_of.get(k)) and not closed,
            circle_opening=bool(trunc_of.get(k)) and closed and not whole,
            whole_contig=whole,
            # CONTIG facts: exact overlap between the two contig ends, inexact terminal repeat
            # (_st492_near_overlap), named seamless circle, and whether any of them holds
            end_overlap_bp=cf["overlap"],
            end_repeat_bp=cf["repeat"],
            seamless_circle=cf["seamless"],
            closed_circle=closed,
            circle_kind=cf["kind"] if closed else "",
            # one of the two counted copies on a doubled contig (see _ST492_TANDEM)
            tandem_copy=bool(info.get("tandem")),
            seq_merged=bool(info.get("seq_merged_from")),
            mob_molecule_type=m[0], mob_primary_cluster_id=m[1], mob_rep_type=m[3],
            location=loc, location_cluster=cl or "", location_note=note,
            candidate_parents=" ".join(map(str, parents)), sig_window=SIG_WINDOW,
            # every track whose gene order contains this one's as a contiguous run
            run_tracks=" ".join(map(str, _ST492_RUNS.get((goi, rank), []))),
            # a genome whose copies on one contig were decided as separate tracks
            # (_ST492_SEPARATE): the figure marks each of those tracks with the other one
            same_genome_decided=(goi, s, contig) in _ST492_SEPARATE,
            country=ctry.get(s, "NA")))

def _st492_write_copies(out):
    p = out.replace("st492_synteny_variants", "st492_synteny_contexts")
    pd.DataFrame(_ST492_COPYROWS).to_csv(p, index=False)
    print(f"wrote {p}  ({len(_ST492_COPYROWS)} copies)")
    p = out.replace("st492_synteny_variants", "st492_synteny_dropped")
    pd.DataFrame(_ST492_DROPPED, columns=["goi", "sample", "contig", "cstart", "contig_len", "in_contig",
                                          "in_contig_len", "containment", "probe_share", "kind"]
                 ).to_csv(p, index=False)
    print(f"wrote {p}  ({len(_ST492_DROPPED)} copies dropped as duplicate contigs)")

# fragment_of may only name a parent on the same molecule and cluster, neither unresolved, with
# at least as many copies, never for a track made entirely of closed circles and never a
# doubled contig as the parent (guide 2.3), and only when exactly ONE parent qualifies (user
# decision of 2026-10-06 for every locus, guide section 8): with two, the builder picked one by
# rank alone. Molecule and cluster are the ones the track text prints (_st492_location), not the
# builder's looser mol_of. The plot script already refused the molecule and circle cases when
# drawing; this makes the table say the same. blaNDM-1: csv v5 (unresolved, GCA_029679155.1) was
# recorded as a fragment of the chromosomal csv v4, and csv v3 (parents 1 and 2) and csv v9 (1
# and 6) were grouped under csv v1. intI1: the 28-copy csv v1 was written as a fragment of a
# 1-copy track and csv v18 had five equally good parents. The qualifying parents stay listed in
# candidate_parents.
_ST492_CANDS = {}       # (goi, variant) -> every qualifying parent
_ST492_RUNS = {}        # (goi, variant) -> every track whose gene order holds this one's as a run
def _st492_parents_ok(goi, i, hits, ranked, loc_counts, runs=()):
    # every track holding the run, whatever its molecule, so the caption can say which others
    # the rule turned down and why (blaNDM-1 display v2's run also lies in v4, v6, v7 and v8)
    _ST492_RUNS[(goi, i)] = list(runs)
    ki = ranked[i - 1][1]
    def info(k): return _ST492_COPYINFO.get((goi,) + k, {})
    def mol(ks): return _st492_location(loc_counts(ks)[0], len(ks))[:2]
    ok = []
    if hits and mol(ki)[0] != "unresolved" and not all(
            info(k).get("whole_contig") and _st492_overlap(read_fna(k[0]).get(k[1], "")) for k in ki):
        ok = [j for j in hits if mol(ranked[j - 1][1]) == mol(ki) and len(ki) <= len(ranked[j - 1][1])
              and not all(info(k).get("tandem") for k in ranked[j - 1][1])]
    _ST492_CANDS[(goi, i)] = ok
    return ok if len(ok) == 1 else []

# A track holding a copy from a genome whose copies on one contig were decided as separate
# tracks (_ST492_SEPARATE) says so in its track text, inside the country parentheses:
# "(Switzerland 1; same genome as variant 22)". In this table the number is the CSV variant,
# like every variant number here; the plot prints the same words with the other track's
# DISPLAY number, so wherever display and CSV numbers agree (intI1 v21/v22) the drawn label
# and this column are the same string, and otherwise they agree through the layout map.
def _st492_same_genome_marks(allg):
    cp = pd.DataFrame(_ST492_COPYROWS)
    if not len(cp): return allg
    sg = cp[cp["same_genome_decided"]][["goi", "sample", "variant"]].drop_duplicates()
    allg = allg.copy()
    for (goi, v), d in sg.groupby(["goi", "variant"]):
        oth = sorted(set(sg[(sg["goi"] == goi) & sg["sample"].isin(d["sample"])]["variant"]) - {v})
        if not oth: continue
        note = "same genome as variant " + ", ".join(map(str, oth))
        m = (allg["goi"] == goi) & (allg["variant"] == v)
        t = allg.loc[m, "track"].iloc[0]
        t = t[:-1] + f"; {note})" if t.endswith(")") else f"{t}  ({note})"
        allg.loc[m, "track"] = t
    return allg

# Display labels. The grouping signature is fixed by then, so nothing here can merge or split
# a context; it only corrects what is printed next to an arrow, adds arrows read from the
# assembly where Bakta called none (each one named below), and records per arrow WHY a name
# carries a prime, a delta or "fam." (column note) so the caption can explain it.
#
# Gene-level rules shared by every locus (user decisions of 2026-10-06, blaNDM-1 fix round):
#  (a) the read-through ORF after blaNDM-1 is computed from the assembly (_st492_fusion)
#  (b) a grey CDS called by an AMR database at >= 80 percent subject coverage is an ARG
#  (c) an IS transposase is printed by family ("IS5 fam.") unless Bakta or ISFinder names
#      the element; a fragment of an element named only by ISFinder gets a prime (ISEc29')
#  (d) aph(3') carries its subclass from the marker call (aph(3')-VI, aph(3')-Ia)
#  (e) a grey CDS prints Bakta's gene name, or the protein symbol its product names, rather
#      than a placeholder cut from the product text ("permease", "NAD", "DNA BP")
#  (f) a CDS cut by a contig end gets a prime (_st492_contig_end_primes)
_ST492_DB = "/home/bals/Documents/dev/BACTolomew/DB/seq/markers"
_ST492_AA = dict(zip((a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"),
                     "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
def _st492_tr(s):
    return "".join(_ST492_AA.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))

_ST492_REF = {}
def _st492_ref(fasta, key):
    """First protein in a local marker database whose header contains key."""
    if (fasta, key) not in _ST492_REF:
        seq, take = [], False
        with open(f"{_ST492_DB}/{fasta}") as fh:
            for line in fh:
                if line.startswith(">"):
                    if take: break
                    take = key in line
                elif take: seq.append(line.strip())
        if not seq: _sys.exit(f"reference {key} not found in {fasta}")
        _ST492_REF[(fasta, key)] = "".join(seq)
    return _ST492_REF[(fasta, key)]

# marker_hits by (sample, contig). Every hit row is a whole Bakta CDS, so a hit belongs to a
# drawn CDS only when the two overlap reciprocally by 90 percent: the opposite-strand ORF that
# overlaps the 3' end of blaNDM-1 by 120 bp must not inherit the blaNDM-1 call.
_ST492_MH = {}
def _st492_hits(sample, contig, a, b, dbs=None):
    if not _ST492_MH:
        for (s, c), dd in mh.groupby(["Sample", "Contig"]):
            _ST492_MH[(s, f"contig_{c}")] = dd
    dd = _ST492_MH.get((sample, contig))
    if dd is None: return mh.iloc[0:0]
    lo = dd[["Start", "Stop"]].min(axis=1).astype(int); hi = dd[["Start", "Stop"]].max(axis=1).astype(int)
    ov = (hi.clip(upper=b) - lo.clip(lower=a) + 1).clip(lower=0)
    keep = (ov >= 0.9 * (b - a + 1)) & (ov >= 0.9 * (hi - lo + 1))
    if dbs is not None: keep &= dd["DB"].astype(str).isin(dbs)
    return dd[keep]

def _st492_frame(nb):
    """(focal midpoint, flipped) of a track. context() flips a track whose focal gene is on the
    minus strand; the focal arrow itself is symmetric about its midpoint, so the flip is read
    from its strand in the GFF, not from its x."""
    c = nb[nb["category"] == "center"].iloc[0]
    cm = (int(c["start"]) + int(c["end"])) / 2
    g = read_gff(c["rep_sample"])
    s0 = g[(g["contig"] == c["contig"]) & (g["start"] == int(c["start"]))]["strand"].iloc[0]
    return cm, s0 == "-"

def _st492_strand0(nb, r):
    """A row's strand on the contig (nb["strand"] is already flipped with the track)."""
    _cm, flip = _st492_frame(nb)
    return r["strand"] if not flip else {"+": "-", "-": "+"}[r["strand"]]

def _st492_newrow(nb, start, end, strand, label, category, note, product=""):
    """An arrow read from the assembly, placed with the track's own focal frame."""
    cm, flip = _st492_frame(nb)
    r = nb.iloc[0].copy()
    r["start"], r["end"], r["product"], r["gene"] = start, end, product, ""
    r["x0"], r["x1"] = (start - cm, end - cm) if not flip else (cm - end, cm - start)
    r["strand"] = strand if not flip else {"+": "-", "-": "+"}[strand]
    r["orient"] = int(r["strand"] == "+")
    r["label"], r["category"], r["note"], r["lab0"] = label, category, note, label
    return r

def _st492_cds_seq(q, a, b, strand):
    s = q[a - 1:b]
    return s if strand == "+" else _rc(s)

# ---- (a) the bleMBL::qacEdelta1 read-through ORF ----------------------------------------------
# Bakta calls a 315 bp OPPOSITE-strand "hypothetical protein" over the 3' end of blaNDM-1, and
# the base builder drew that ORF under the fusion name. The fusion itself lies on the
# blaNDM-1 strand: read from the base after the blaNDM-1 stop codon to the first in-frame
# stop, it is one codon, then BleMBL (WP_004201167.1) from its Met, running into the
# C-terminus of qacEdelta1. It is drawn from that Met to its stop and named only when both
# halves check out against sequence: an exact BleMBL N-terminal prefix of at least 30 aa
# (from the local AMRFinderPlus database) and the last 38 aa equal to the C-terminus of a
# qacEdelta1 of the same genome, or, in a genome with none, of the CARD qacEdelta1 reference
# (2026-10-07). Otherwise Bakta's ORF is drawn as a grey, unlabelled HP.
_ST492_FUSION_LAB = "ΔbleMBL::ΔqacEΔ1"
_ST492_QACE_PROTS = {}
def _st492_qace_prots(sample):
    """The genome's own qacEdelta1 proteins as [(contig, lo, hi, protein)]: Bakta's calls, or,
    in a genome where Bakta called none, every copy located by the last 20 aa of the CARD
    qacEdelta1 reference (ARO:3005010), whatever interrupts it upstream (a scaffold gap)."""
    if sample in _ST492_QACE_PROTS: return _ST492_QACE_PROTS[sample]
    g = read_gff(sample); seqs = read_fna(sample); out = []
    m = (g["gene"].astype(str).str.lower() == "qacedelta1") | \
        g["product"].astype(str).str.contains("QacE delta 1", case=False, regex=False)
    for _, r in g[m].iterrows():
        q = seqs[r["contig"]]
        out.append((r["contig"], int(r["start"]), int(r["end"]),
                    _st492_tr(_st492_cds_seq(q, int(r["start"]), int(r["end"]), r["strand"])).rstrip("*")))
    if not out:
        tail = _st492_ref("CARD_AMR.fasta", "|qacEdelta1 ")[-20:]
        for contig, q in seqs.items():
            for strand, qq in (("+", q), ("-", _rc(q))):
                for f in range(3):
                    p = _st492_tr(qq[f:])
                    for mt in _re.finditer(tail, p):
                        i1 = mt.end()                              # residues up to the C-terminus
                        if p[i1:i1 + 1] != "*": continue
                        a = f + 3 * max(0, i1 - 60); b = f + 3 * (i1 + 1)   # 0-based, half open
                        lo, hi = (a + 1, b) if strand == "+" else (len(q) - b + 1, len(q) - a)
                        out.append((contig, lo, hi, p[max(0, i1 - 60):i1]))
    _ST492_QACE_PROTS[sample] = out
    return out

def _st492_fusion(nb, sample, contig):
    q = read_fna(sample).get(contig, "")
    ble = _st492_ref("AMRF_core-AMR.fasta", "WP_004201167.1")
    drop, add = set(), []
    for i, r in nb[nb["label"] == "blaNDM-1"].iterrows():
        st, en, sd = int(r["start"]), int(r["end"]), _st492_strand0(nb, r)
        down = q[en:] if sd == "+" else _rc(q[:st - 1])
        prot = _st492_tr(down)
        # Bakta's own ORF on the other strand over the junction: replaced or demoted below
        opp = [j for j, o in nb.iterrows() if j != i and o["category"] != "center"
               and _st492_strand0(nb, o) != sd and 300 <= int(o["end"]) - int(o["start"]) + 1 <= 330
               and int(o["start"]) <= en + 400 and int(o["end"]) >= st - 400]
        ok, note = False, "fusion unconfirmed"
        if "*" not in prot:                    # the contig ends before the ORF does
            nb.at[i, "note"] = f"fusion ORF past contig end ({len(down)} bp after the stop)"
        else:
            orf = prot.split("*")[0]; k = orf.find("M")
            if 0 <= k <= 2:
                p = orf[k:]; nt = 3 * (len(p) + 1)                # Met to stop codon inclusive
                if sd == "+": a = en + 1 + 3 * k; b = a + nt - 1
                else: b = st - 1 - 3 * k; a = b - nt + 1
                lb = next((x for x in range(min(len(p), len(ble))) if p[x] != ble[x]), min(len(p), len(ble)))
                own = [qp for c, lo, hi, qp in _st492_qace_prots(sample)
                       if not (c == contig and lo <= b and hi >= a)]
                # A genome with no qacEdelta1 of its own (blaNDM-1 csv v2) is checked against the
                # CARD qacEdelta1 reference (ARO:3005010) instead, which equals every qacEdelta1
                # Bakta calls in this run (user decision of 2026-10-07)
                ref = not own
                if ref: own = [_st492_ref("CARD_AMR.fasta", "|qacEdelta1 ")]
                lq = max((next((x for x in range(min(len(p), len(qp))) if p[-1 - x] != qp[-1 - x]),
                                min(len(p), len(qp))) for qp in own), default=0)
                if lb >= 30 and lq >= 38 and lb + lq >= len(p):
                    # BleMBL codons + qacE codons; the note says when the end was checked
                    # against the reference, so the caption can say so
                    ok, note = True, f"fusion {len(p) - lq}+{lq}" + ("; CARD qacEdelta1 reference" if ref else "")
                    add.append(_st492_newrow(nb, a, b, sd, _ST492_FUSION_LAB, "ARG", note,
                                             "bleMBL::qacEdelta1 read-through ORF (from the assembly)"))
        if ok and opp:
            # the caption names the Bakta ORF this arrow replaces by its length
            L = sorted({int(nb.at[j, "end"]) - int(nb.at[j, "start"]) + 1 for j in opp})
            add[-1]["note"] = note + "; Bakta ORF " + "/".join(map(str, L)) + " bp replaced"
        for j in opp:
            if ok: drop.add(j)
            elif nb.at[j, "label"] in (_ST492_FUSION_LAB, "HP"):
                nb.at[j, "label"], nb.at[j, "category"], nb.at[j, "note"] = "HP", "HP", note
    nb = nb.drop(index=list(drop))
    if add: nb = pd.concat([nb, pd.DataFrame(add)], ignore_index=True)
    return nb.sort_values("x0").reset_index(drop=True)

# ---- (b) and (d): colour from AMR marker calls ------------------------------------------------
_ST492_AMR_DBS = ("AMRF_core-AMR", "AMRF_plus-AMR", "CARD_AMR")
_ST492_MOBILE_DBS = ("CGE-MGE-Prot_MOBILE", "ISFinder-Prot_MOBILE")
# Called by an AMR database but deliberately drawn grey, each with its caption reason (the plot
# script's KEEP_GREY carries the same text; keep the two lists in step)
_ST492_KEEP_GREY = {"cxpE": "intrinsic chloramphenicol efflux gene, AMRFinderPlus plus tier",
                    "merA": "mercury resistance, AMRFinderPlus plus tier, metal class",
                    "merR": "mercury resistance, AMRFinderPlus plus tier, metal class"}
# A resistance MUTATION in a housekeeping gene (AMRFinderPlus parC_S84F, gyrA_..., CARD
# protein-variant models "VM|") does not make the gene an acquired ARG: such calls are left
# out of rule (b), so parC on oprD keeps its grey (an open question, kept as it was).
def _st492_acquired(h):
    return h[~h["ResistanceAlleles"].astype(str).str.fullmatch(r"\w+_[A-Z]\d+[A-Z*]")
             & ~h["SubjectIdentifier"].astype(str).str.startswith("VM|")]
def _st492_amr_name(h):
    """Gene name of an AMR call: AMRFinderPlus symbol, else any database's, else the family."""
    for db in _ST492_AMR_DBS:
        g = h[(h["DB"].astype(str) == db) & h["Gene"].notna()]
        if len(g): return str(g.iloc[0]["Gene"])
    f = h[h["DB"].astype(str) == "AMRF_core-AMR"]["ResistanceGeneFamilies"].dropna()
    return str(f.iloc[0]) if len(f) else None

# ---- (c) IS elements ---------------------------------------------------------------------------
# A family-only call ("IS4 family transposase") is named from ISFinder when the hit is close
# (>= 95 percent identity over >= 80 percent of the CDS), else from the partner CDS of a
# two-gene element that names it (the istA of an "IS21-like element IS1326" istB pair), else
# printed as the family. An ISFinder name gets a prime only for a FRAGMENT: under 90 percent
# of the element's transposase matched AND the CDS shorter than 90 percent of it. intI1 display
# ISEc29 (GCA_964331735.1, 1,188 bp for a 394 aa transposase, 88 percent matched) is full
# length with a divergent C-terminus, so it keeps the plain name.
# The same prime applies to an element Bakta names itself (user decision of 2026-10-07): the
# 606 bp IS15 copies (ISFinder IS15 at 86 percent of its 234 aa transposase; a whole IS15 CDS is
# 705 bp) print IS15', like ISEc29' and ISAba1'.
def _st492_is_cut(best):
    return (float(best["Subject_coverage"]) < 90
            and float(best["Q_len"]) < 0.9 * 3 * float(best["S_len"]))

# the lengths go in the note so the caption can state the prime rule with the numbers
def _st492_is_frag_note(best):
    return (f"IS fragment {float(best['Subject_coverage']):.0f}; {float(best['Q_len']):.0f} bp CDS "
            f"for a {float(best['S_len']):.0f} aa transposase")

def _st492_is_name(sample, contig, r, lab, p, nb):
    fam = _re.fullmatch(r"(IS\w*?\d+\w*|Tn\d+) family transposase", p)
    if not fam and _re.fullmatch(r"IS[A-Za-z]*\d+[A-Za-z]*", lab):
        h = _st492_hits(sample, contig, int(r["start"]), int(r["end"]), ("ISFinder-Prot_MOBILE",))
        h = h[h["SubjectIdentifier"].astype(str).str.split("//").str[-1] == lab]
        if len(h):
            best = h.sort_values("Subject_coverage", ascending=False).iloc[0]
            if _st492_is_cut(best):
                return lab + "'", _st492_is_frag_note(best)
        return lab, ""
    if not fam or lab == "IS91/ISCR1":        # ISCR1 is named by its CGE-MGE match (label_cat)
        return lab, ""
    h0 = _st492_hits(sample, contig, int(r["start"]), int(r["end"]), ("ISFinder-Prot_MOBILE",))
    h = h0
    if len(h):
        qcov = (h["Q_end"].astype(float) - h["Q_start"].astype(float) + 1) / h["Q_len"].astype(float)
        h = h[(h["Identity"].astype(float) >= 95) & (qcov >= 0.8)]
    if len(h):
        best = h.sort_values("Subject_coverage", ascending=False).iloc[0]
        el = str(best["SubjectIdentifier"]).split("//")[-1]
        if _re.fullmatch(r"IS[A-Za-z]*\d+[A-Za-z]*", el):
            cov = float(best["Subject_coverage"])
            if _st492_is_cut(best): return el + "'", _st492_is_frag_note(best)
            return el, (f"IS element from ISFinder {el} {float(best['Identity']):.0f}/{cov:.0f}; "
                        f"{float(best['Q_len']):.0f} bp CDS for a {float(best['S_len']):.0f} aa transposase")
    f0 = fam.group(1)
    partner = [m.group(1) for _, o in nb.iterrows()
               if int(o["start"]) <= int(r["end"]) + 20 and int(o["end"]) >= int(r["start"]) - 20
               and (int(o["start"]), int(o["end"])) != (int(r["start"]), int(r["end"]))
               and (m := _re.search(rf"\b{_re.escape(f0)}-like element (IS[A-Za-z]*\d+[A-Za-z]*)",
                                    _unquote(str(o["product"]))))]
    if partner:
        return partner[0], f"IS element named by its partner CDS ({f0}-like element {partner[0]})"
    # The note keeps ISFinder's best match below the naming cut, so the caption can say "below
    # 95 percent" rather than "no match" (IS5 fam. of blaNDM-1: ISEc35 at 90 percent identity).
    # Its coverage is of the CDS, the quantity the naming cut tests (the 3,441 bp Tn3 fam. CDS of
    # intI1 carries a 528 bp IS15 match, 15 percent of it, on its other strand; the 83 percent
    # printed before was IS15's own coverage and read as nearly passing the cut).
    if len(h0):
        q0, q1 = h0["Q_start"].astype(float), h0["Q_end"].astype(float)
        h0 = h0.assign(qcov=((q1 - q0).abs() + 1) / h0["Q_len"].astype(float) * 100, rev=q1 < q0)
        b = h0.sort_values(["Identity", "qcov"], ascending=False).iloc[0]
        return (f"{f0} fam.", f"IS family; ISFinder best {str(b['SubjectIdentifier']).split('//')[-1]} "
                f"{float(b['Identity']):.0f}/{float(b['qcov']):.0f} of the CDS"
                + ("; other strand" if b["rev"] else ""))
    return f"{f0} fam.", "IS family; no ISFinder match"

# ---- (e) gene names ----------------------------------------------------------------------------
# The first protein symbol the product names (AraJ, RidA, CzcO), never one it merely contains
# ("contains ChlI-like ... domains") nor a hyphenated qualifier ("TonB-dependent").
_ST492_SYM = _re.compile(r"(?<!contains )(?<![-\w/])([A-Z][a-z]{2}[A-Z][0-9]?)(?![-\w])")
_ST492_GENE_OK = _re.compile(r"[a-z]{3}[A-Z]?[0-9]*")        # not aRO8, eAL, rRM, nCS2, rPS27A
# Products Bakta gives without a gene= although the subunit is unambiguous
# Names set on purpose upstream, which rule (e) must not replace
_ST492_KEPT_NAMES = {"DNA rep", "septicolysin"}
_ST492_PRODUCT_GENE = {"NAD(P)(+) transhydrogenase (Re/Si-specific) subunit beta": "pntB"}
def _st492_gene_name(p, gene, lab):
    """Name for a grey CDS, or None to keep the label. Bakta's gene= wins when the product names
    it ("DnaA regulatory inactivator Hda" is hda); otherwise a gene-like label the product names
    stays; otherwise the first protein symbol the product names; otherwise, for a placeholder,
    Bakta's gene=. Only the head of "X associated with Y" names X (CzcO, not czcD)."""
    if p in _ST492_PRODUCT_GENE: return _ST492_PRODUCT_GENE[p]
    ok = bool(gene) and bool(_ST492_GENE_OK.fullmatch(gene))
    if _re.search(r"\b(domain|motif)\b|DUF\d", p):
        # A domain-only product prints Bakta's gene= when that is a gene symbol (DUF882 is ycbK,
        # DUF308 hdeD, the Kef-type K+ antiporter domain rosB, all on blaOXA-72), else keeps
        # the domain name: a lower-cased domain abbreviation (rRM) is not a gene symbol
        return gene if ok else None
    head = _re.split(r"\bassociated with\b", p)[0]
    def named(sym):
        return bool(sym) and bool(_re.search(rf"(?<![-\w]){_re.escape(sym[0].upper() + sym[1:])}(?![-\w])", head))
    if ok and named(gene): return gene
    if _ST492_GENE_OK.fullmatch(lab) and named(lab): return None
    m = _ST492_SYM.search(head)
    if m: return m.group(1)[0].lower() + m.group(1)[1:]
    if _ST492_GENE_OK.fullmatch(lab) or lab in _ST492_KEPT_NAMES: return None
    return gene if ok else None

# ---- (f) CDS cut by a contig end ---------------------------------------------------------------
# Generic test, no locus or gene named. The outermost CDS of a contig is cut when the end facing
# the contig edge is its 3' end and it has no stop codon there (a partial ORF), or when that end
# is its 5' end, the reading frame stays open from its start codon all the way to the edge, and
# something says the protein is incomplete: it lacks a start codon (Prodigal ran it off the
# edge), it is at least 30 bp shorter than the commonest length of the same name drawn at this
# locus, or its best marker match covers under 95 percent of the subject. repM on the two
# closed AC237 circles of blaNDM-1 (768 of 828 bp, 93 percent of CGE-MGE repM, start 41 bp from
# where the circle was opened) is one. _st492_contained_primes, run first, catches the pieces
# this codon test cannot see.
_ST492_FROM_ASSEMBLY = ("fusion", "contig end", "scaffold gap", "ambiguity codes")
def _st492_contig_end_primes(nb, sample, contig, maxlen):
    q = read_fna(sample).get(contig, ""); L = len(q)
    g = read_gff(sample); g = g[g["contig"] == contig]
    if not L or not len(g): return nb
    first, last = int(g["start"].min()), int(g["end"].max())
    for i, r in nb.iterrows():
        lab = str(r["label"])
        if (r["category"] in ("center", "HP") or not lab or lab.endswith("'")
                or str(r["note"]).startswith(_ST492_FROM_ASSEMBLY)): continue
        st, en, sd = int(r["start"]), int(r["end"]), _st492_strand0(nb, r)
        for side in (("L",) if st == first else ()) + (("R",) if en == last else ()):
            five = (side == "L") == (sd == "+")
            if not five:
                stop = q[en - 3:en] if sd == "+" else _rc(q[st - 1:st + 2])
                cut = _ST492_AA.get(stop, "X") != "*"
            else:
                up = q[(st - 1) % 3:st - 1] if sd == "+" else _rc(q[en:])
                if sd == "-": up = up[len(up) % 3:]
                first_codon = q[st - 1:st + 2] if sd == "+" else _rc(q[en - 3:en])
                h = _st492_hits(sample, contig, st, en)
                cov = float(h["Subject_coverage"].max()) if len(h) else 100.0
                base = str(r.get("lab0", lab))
                cut = "*" not in _st492_tr(up) and (
                    first_codon not in ("ATG", "GTG", "TTG") or cov < 95
                    or (base != "HP" and en - st + 1 < maxlen.get(base, 0) - 30))
            if cut:
                nb.at[i, "label"] = lab + "'"
                nb.at[i, "note"] = "; ".join(x for x in (str(r["note"]), "contig end") if x)
                break
    return nb

# The second test for (f), by sequence (intI1 fix round): the outermost CDS on one side of a
# contig is cut when its sequence lies inside a CDS at least 10 percent longer drawn elsewhere
# at this locus (>= 90 percent of its 21-mers, either strand) and the missing part cannot fit
# between it and the contig end. It is then named after that CDS: ISEc29' for Bakta's
# "Transposase DDE domain" and "IS4 family transposase" at intI1 contig ends, DNA rep' and
# znuD2' on the blaOXA-72 circles. Not a length or position heuristic alone: an intact gene near
# a contig end is not primed. Arrows inside an inexact terminal repeat of a whole circle are
# left alone: the same sequence continues at the other contig end (blaNDM-1 csv v6 IS15).
_ST492_POOL = {}
def _st492_kmers(q, k=21):
    return {q[i:i + k] for i in range(0, len(q) - k + 1)}

def _st492_pool(goi, rep):
    """Every drawn CDS of this locus (all representatives), with its display label and its
    sequence: the references for the contained-prime and the N-gap rules."""
    if goi in _ST492_POOL: return _ST492_POOL[goi]
    out = []
    for s, nb0 in rep.values():
        nb0 = _st492_label_rules(goi, nb0.assign(rep_sample=s))
        seq = read_fna(s)
        for _, r in nb0.iterrows():
            if r["category"] in ("center", "HP") or not r["label"]: continue
            q = seq.get(r["contig"], "")[int(r["start"]) - 1:int(r["end"])]
            out.append(dict(sample=s, contig=r["contig"], start=int(r["start"]), label=r["label"],
                            category=r["category"], seq=q,
                            km=_st492_kmers(q) | _st492_kmers(_rc(q))))
    _ST492_POOL[goi] = out
    return out

def _st492_contained_primes(nb, sample, contig, pool):
    seqs = read_fna(sample); q0 = seqs.get(contig, ""); clen = len(q0)
    g = read_gff(sample); gc = g[g["contig"] == contig]
    if not clen or not len(gc): return nb
    rep_len = _st492_near_overlap(q0)[0]
    for i, r in nb.iterrows():
        lab = str(r["label"])
        if (r["category"] == "center" or lab.endswith("'")
                or str(r["note"]).startswith(_ST492_FROM_ASSEMBLY)): continue
        a, b = int(r["start"]), int(r["end"])
        if rep_len and (b <= rep_len or a > clen - rep_len): continue
        if a == int(gc["start"].min()): d = a - 1
        elif b == int(gc["end"].max()): d = clen - b
        else: continue
        q = q0[a - 1:b]; L = len(q)
        km = _st492_kmers(q)
        if not km: continue
        best = None
        for c in pool:
            if len(c["seq"]) < 1.1 * L or (c["sample"], c["contig"], c["start"]) == (sample, contig, a):
                continue
            f = len(km & c["km"]) / len(km)
            if f >= 0.9 and d < len(c["seq"]) - L and (best is None or f > best[0]):
                best = (f, c)
        if best:
            # renamed: notes on the old name no longer apply. The note keeps what the name rests
            # on (the caption states it: ISEc29' of intI1 has no marker call of its own)
            prod = _unquote(str(r.get("product", "") or "")).replace("%2C", ",").strip()
            nb.at[i, "note"] = (f"named by sequence: {best[0] * 100:.0f} percent of its 21-mers in "
                                f"{best[1]['label'].rstrip(chr(39))} of {best[1]['sample']}"
                                + (f"; Bakta: {prod}" if prod and prod != "nan" else "") + "; contig end")
            nb.at[i, "label"] = best[1]["label"].rstrip("'") + "'"
            nb.at[i, "category"] = best[1]["category"]
    return nb

# ---- named data corrections, each read from the assembly and checked ---------------------------
# blaNDM-1 display decisions of 2026-10-06. Every one is keyed by genome and contig, applies to
# whichever locus draws that contig, and stops the run if the sequence no longer says what the
# comment says.
#
# GCA_029679155.1 contig_63 (18,542 bp) ends inside msr(E): the last 1,341 bp read in frame from
# the msr(E) start codon at 17,202 without a stop, and the first 1,267 of them match a complete
# msr(E) (1,476 bp) with one mismatch. The last 74 bp (18,469-18,542) are the end of an IS15
# (IS26) copy: identical to 21,318-21,391 of the AC237 circle GCA_031157295.1 contig_90, over
# the 3' end of its IS15 transposase gene and the inverted repeat beyond it, and found again at
# the start of contig_96 (820 bp), the genome's separate IS15 contig. Whether IS15 interrupts msr(E) here or the assembly broke at
# the repeat cannot be told from the contig, so the caption only says what the sequence is.
# Bakta called no CDS there, so the track looked as if it ended after ISEc29. Drawn as msr(E)'
# over the 1,267 bp that are msr(E).
#   (genome, contig) -> (gene, expected bp, element the stretch after it must belong to)
_ST492_END_FRAGMENTS = {("GCA_029679155.1", "contig_63"): ("msr(E)", 1267, "IS15")}
def _st492_ref_cds(gene):
    """nt of the first complete CDS AMRFinderPlus calls `gene` at 100/100 in this run."""
    h = mh[(mh["Gene"].astype(str) == gene) & (mh["DB"].astype(str) == "AMRF_core-AMR")
           & (mh["Subject_coverage"] >= 100) & (mh["Identity"] >= 100)]
    h = h.sort_values(["Sample", "Contig", "Start"]).iloc[0]
    a, b = sorted((int(h["Start"]), int(h["Stop"])))
    return _st492_cds_seq(read_fna(h["Sample"])[f"contig_{h['Contig']}"], a, b, h["Strand"])

def _st492_end_fragment(nb, sample, contig):
    spec = _ST492_END_FRAGMENTS.get((sample, contig))
    if not spec: return nb
    gene, expect, elem = spec
    q = read_fna(sample)[contig]; ref = _st492_ref_cds(gene)
    gc = read_gff(sample); gc = gc[gc["contig"] == contig]
    for strand, qq in (("+", q), ("-", _rc(q))):
        i = qq.find(ref[:60])
        if i < 0: continue
        mism, n = 0, 0
        for x, y in zip(qq[i:], ref):
            if x != y:
                mism += 1
                if mism > 1: break
            n += 1
        if n != expect:
            _sys.exit(f"{sample} {contig}: {gene} fragment is {n} bp, the correction expects {expect}")
        a, b = (i + 1, i + n) if strand == "+" else (len(q) - i - n + 1, len(q) - i)
        # Drawn only as the next arrow after the window's outermost CDS: a Bakta CDS between the
        # window and the fragment means the fragment lies outside the +/- WINDOW CDS drawn.
        # intI1 reaches this contig from its other end, where ISEc29 (15,628-16,803) is CDS 15
        # and msr(E)' would be CDS 16, so there it is not drawn; blaNDM-1 draws ISEc29, then it.
        lo, hi = int(nb["start"].min()), int(nb["end"].max())
        if ((gc["start"] > hi) & (gc["end"] < a)).any() or ((gc["end"] < lo) & (gc["start"] > b)).any():
            return nb
        # what follows the match up to the contig end must lie on a Bakta-called copy of elem
        # elsewhere in the genome (within 100 bp of that CDS), or the comment above is stale
        rest = qq[i + n:]
        where = []
        if rest:
            g = read_gff(sample)
            for c2, s2 in read_fna(sample).items():
                if c2 == contig: continue
                for s2s in (s2, _rc(s2)):
                    j = s2s.find(rest)
                    if j < 0: continue
                    lo2, hi2 = (j + 1, j + len(rest)) if s2s is s2 else (len(s2) - j - len(rest) + 1, len(s2) - j)
                    gg = g[(g["contig"] == c2) & g["product"].astype(str).str.contains(elem, regex=False)
                           & (g["start"] <= hi2 + 100) & (g["end"] >= lo2 - 100)]
                    if len(gg): where.append(c2)
        if rest and not where:
            _sys.exit(f"{sample} {contig}: the {len(rest)} bp after the {gene} fragment are no longer an {elem} end")
        row = _st492_newrow(nb, a, b, strand, gene + "'", "ARG",
                            f"contig end, read from the assembly: {n:,} of {len(ref):,} bp"
                            + (f"; then {len(rest)} bp of an {elem} end" if rest else ""))
        return pd.concat([nb, pd.DataFrame([row])], ignore_index=True).sort_values("x0").reset_index(drop=True)
    _sys.exit(f"{sample} {contig}: {gene} start no longer found")

# GCA_031157325.1 contig_75 carries a 1,113 base block of IUPAC ambiguity codes at 16,033-17,508
# over the sul1 and ISCR1 next to the integron. Bakta read through it as one 1,905 bp
# "hypothetical protein" (15,998-17,902) and called only the 3' part of ISCR1, as an 840 bp
# "IS91 family transposase" (15,117-15,956). Drawn as they were, the 1,905 bp ORF put the
# sul1-ISCR1 spacer and the 5' 700 bp of ISCR1 under the name sul1'. The contig is the same
# 22,134 bp circle as GCA_031157295.1 contig_90, opened at the same base: the two are equal at
# every position where both call a base, apart from one substitution at 13,588, and every
# ambiguity code admits the sister's base. So sul1 and ISCR1 are drawn at the sister's Bakta
# coordinates (sul1 17,063-17,902, ISCR1 15,117-16,658, both minus strand), as sul1' and
# IS91/ISCR1': ORFs broken by ambiguity codes, not genes that are absent. The element keeps the
# name it has on every other track, IS91/ISCR1, with a prime (2026-10-07: "ISCR1'" read as a
# different element). Each span must hold ambiguity codes and match the sister there base for
# base, or the run stops.
#   (genome, contig) -> ((sister genome, sister contig), {(start, end): (label, category)})
_ST492_AMBIG = {("GCA_031157325.1", "contig_75"): (("GCA_031157295.1", "contig_90"),
                                                   {(17063, 17902): ("sul1'", "ARG"),
                                                    (15117, 16658): ("IS91/ISCR1'", "MGE")})}
_ST492_IUPAC = {"R": "AG", "Y": "CT", "K": "GT", "M": "AC", "S": "CG", "W": "AT", "B": "CGT",
                "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
def _st492_ambiguity(nb, sample, contig):
    spec = _ST492_AMBIG.get((sample, contig))
    if not spec: return nb
    (s2, c2), genes = spec
    q, q2 = read_fna(sample)[contig], read_fna(s2)[c2]
    g2 = read_gff(s2); g2 = g2[g2["contig"] == c2]
    add, ov = [], 0
    for (a, b), (lab, cat) in genes.items():
        cds = g2[(g2["start"] == a) & (g2["end"] == b)]
        seg, seg2 = q[a - 1:b], q2[a - 1:b]
        amb = sum(ch not in "ACGT" for ch in seg)
        same = len(q) == len(q2) and all(x == y or y in _ST492_IUPAC.get(x, "") for x, y in zip(seg, seg2))
        if len(cds) != 1 or not amb or not same:
            _sys.exit(f"{sample} {contig} {a}-{b}: ambiguity-code correction no longer matches the data")
        ov = ov + (nb["end"].clip(upper=b) - nb["start"].clip(lower=a) + 1).clip(lower=0)
        add.append(_st492_newrow(nb, a, b, cds.iloc[0]["strand"], lab, cat,
                                 f"ambiguity codes: {amb:,} in the gene; drawn at the coordinates of "
                                 f"{s2} {c2}", product=str(cds.iloc[0]["product"])))
    # Bakta's own calls lying mostly inside the redrawn genes (the read-through HP, the 3' piece
    # of ISCR1) give way; the flanking qacEdelta1, which overlaps sul1 by 7 bp, stays
    drop = nb.index[(ov >= 0.5 * (nb["end"] - nb["start"] + 1)) & (nb["category"] != "center")]
    if len(drop) != 2:
        _sys.exit(f"{sample} {contig}: expected Bakta's 2 calls over the ambiguity block, found {len(drop)}")
    nb = nb.drop(index=list(drop))
    return pd.concat([nb, pd.DataFrame(add)], ignore_index=True).sort_values("x0").reset_index(drop=True)

# GCA_031156335.1 contig_80 (22,849 bp) carries 767 bp of its AC237 circle at both ends. The copy
# at the start (1-766) equals the sister circles base for base; the copy at the end
# (22,083-22,849) has one extra A (contig position 22,626) inside the IS15 there, which shifts
# its frame, so Bakta calls a 546 bp "IS6 family transposase" (22,135-22,680) where the start
# copy has IS15 (53-658). It is the same IS15 copy drawn at both contig ends, so it is named
# IS15' (user decision of 2026-10-07), not "IS6 fam.". Named; stops the run if the contig no
# longer has a 767 bp terminal repeat one edit apart, or the CDS has moved out of it.
#   (genome, contig, start, end) -> (label, category, note)
_ST492_RENAME = {("GCA_031156335.1", "contig_80", 22135, 22680):
                 ("IS15'", "MGE", "frameshifted IS15 copy: one extra A in the terminal repeat; "
                                  "Bakta: IS6 family transposase, 546 bp")}
def _st492_renames(nb, sample, contig):
    for (s, c, a, b), (lab, cat, note) in _ST492_RENAME.items():
        if (s, c) != (sample, contig): continue
        cf = _st492_contig_circle(s, c); L = len(read_fna(s)[c])
        m = (nb["start"] == a) & (nb["end"] == b)
        if (cf["repeat"], cf["edits"]) != (767, 1) or a <= L - 767 or m.sum() != 1:
            _sys.exit(f"{s} {c} {a}-{b}: the IS15' rename no longer matches the data")
        nb.loc[m, "label"], nb.loc[m, "category"], nb.loc[m, "note"] = lab, cat, note
    return nb

# A gene split by an N scaffold gap. Bakta calls no CDS across a run of N, so the gene is
# missing from the track and the track looks like a different arrangement: qacEdelta1 in the four
# Serbian assemblies of the 21,594 bp AC237 circle (GCA_021570535.1, GCA_021570835.1,
# GCA_021571095.1, GCA_021571215.1), with a 96 to 99 bp N run and an 81 bp duplicated flank (user
# decision of 2026-10-06). For every run of >= 10 N inside the window that no CDS covers, each
# gene drawn elsewhere at this locus is tried: its first 25 bp must lie before the gap and its
# last 25 bp after it, in one orientation and within the gene's own length of the gap. A hit is
# drawn as that gene with a prime, spanning both anchors; gap_x0/gap_x1 place the gap, and the
# note says what the span holds besides the gene. Display only: the track keeps its signature.
def _st492_gap_split(nb, sample, contig, pool):
    q = read_fna(sample).get(contig, "")
    lo, hi = int(nb["start"].min()), int(nb["end"].max())
    cm, flip = _st492_frame(nb)
    X = (lambda pos: cm - pos) if flip else (lambda pos: pos - cm)
    add = []
    for m in _re.finditer(r"N{10,}", q[lo - 1:hi]):
        g0, g1 = lo + m.start(), lo + m.end() - 1
        if ((nb["start"] <= g1) & (nb["end"] >= g0)).any(): continue
        seen = set()
        for c in sorted(pool, key=lambda c: -len(c["seq"])):
            if c["label"] in seen or len(c["seq"]) < 100 or c["label"].endswith("'"): continue
            seen.add(c["label"])
            for strand, t in (("+", c["seq"]), ("-", _rc(c["seq"]))):
                reach = len(t) + 200
                left = q.rfind(t[:25], max(0, g0 - 1 - reach), g0 - 1)
                right = q.find(t[-25:], g1, g1 + reach)
                if left < 0 or right < 0: continue
                s0, s1 = left + 1, right + 25
                if ((nb["start"] < s1 - 50) & (nb["end"] > s0 + 50)).any(): continue
                # what the drawn span holds besides the gene: the N run, and any stretch the
                # assembler placed on both sides of it (the longest end of the left part that
                # is also the start of the right part), so the caption can say why it is longer
                lp, rp = q[s0 - 1:g0 - 1], q[g1:s1]
                dup = next((k for k in range(min(len(lp), len(rp)), 0, -1) if lp[-k:] == rp[:k]), 0)
                row = _st492_newrow(nb, s0, s1, strand, c["label"] + "'", c["category"],
                                    f"scaffold gap {g1 - g0 + 1} N; {dup} bp on both sides of it; "
                                    f"{s1 - s0 + 1} bp drawn, {len(t)} bp intact",
                                    product=f"{c['label']} split by an N gap")
                row["gap_x0"], row["gap_x1"] = sorted((X(g0), X(g1)))
                add.append(row)
                break
            if add and add[-1]["gap_x0"] == min(X(g0), X(g1)): break
    if add: nb = pd.concat([nb, pd.DataFrame(add)], ignore_index=True)
    return nb.sort_values("x0").reset_index(drop=True)

# ---- per-track facts for the drawing and the caption -------------------------------------------
def _st492_track_cols(nb, sample, contig):
    """gaps: N runs of 10 or more inside the drawn span, as display x ranges "x0:x1;...".
    contig_end: which drawn ends are contig ends (L, R or LR in display orientation), or
    "circle" for a whole closed circle (exact terminal overlap), whose ends are only the
    assembler's opening point, "circle, seamless" for a named closed circle whose ends join
    without overlap (_ST492_SEAMLESS), or "circle, <k> bp terminal repeat, <e> edit(s)" for a whole
    circle whose ends repeat k bp inexactly (_st492_near_overlap). An arrow lying inside such a
    repeat gets the note "terminal repeat": the same sequence is drawn at the other end too."""
    q = read_fna(sample).get(contig, ""); cm, flip = _st492_frame(nb)
    lo, hi = int(nb["start"].min()), int(nb["end"].max())
    gaps = []
    for mt in _re.finditer(r"N{10,}", q[lo - 1:hi]):
        a, b = lo + mt.start(), lo + mt.end() - 1
        x0, x1 = (a - cm, b - cm) if not flip else (cm - b, cm - a)
        gaps.append(f"{x0:g}:{x1:g}")
    g = read_gff(sample); g = g[g["contig"] == contig]
    left, right = lo <= int(g["start"].min()), hi >= int(g["end"].max())
    cf = _st492_contig_circle(sample, contig)
    rep, edits = (cf["repeat"], cf["edits"]) if left and right else (0, 0)
    doubled = any((s, c) == (sample, contig) for _g, s, c in _ST492_TANDEM)
    # An arrow inside the contig's inexact terminal repeat is noted whatever the window reaches:
    # a contig fact, so the same CDS reads the same under every locus (GCA_031156335.1 contig_80
    # IS15' is drawn whole-contig on blaNDM-1 and from one end on intI1)
    if cf["repeat"]:
        L = len(q)
        for i, r in nb.iterrows():
            if r["category"] != "center" and (int(r["end"]) <= cf["repeat"] or int(r["start"]) > L - cf["repeat"]):
                nb.at[i, "note"] = "; ".join(x for x in (str(r["note"] or ""), "terminal repeat") if x)
    if left and right and cf["overlap"]: end = "circle"
    # a named seamless circle drawn whole (blaOXA-72 csv v2): its ends are the opening point
    elif left and right and cf["seamless"]: end = "circle, seamless"
    elif rep:
        end = f"circle, {rep} bp terminal repeat, {edits} edit{'' if edits == 1 else 's'}"
    else:
        if flip: left, right = right, left
        end = ("L" if left else "") + ("R" if right else "")
        # A window that reaches only ONE end of a contig whose ends overlap (a closed circle,
        # exactly or as an inexact terminal repeat) stops where the assembler opened the
        # circle, not at an assembly break: intI1 csv v15, v16 and v17 are the 22 kb AC237
        # circles drawn whole on blaNDM-1, reached from intI1 near one contig end. The window
        # is not wrapped round the circle (that would change the +/- 9 CDS signatures).
        if end and cf["closed"] and not doubled: end += ", circle opening"
    nb["gaps"] = ";".join(gaps); nb["contig_end"] = end
    return nb

# CGE-MGE and ISFinder MOBILE calls on a drawn CDS, for the caption and the verifier. The
# CGE-MGE protein set holds every protein of its reference elements, cargo included: in intI1
# every sul1, qacEdelta1, armA, aadA and dfrA12 arrow has a CGE-MGE match (ICE6440 Sul1, Tn6279
# ArmA ...), as do repM (Tn6183 RepAciN), merC/merT (Tn502, Tn6130) and a GCN5 acetyltransferase
# (Tn6019 "unknown"). A match therefore says the CDS belongs to a known element, not that it
# moves DNA, and it does not by itself colour a CDS MGE: the colour stays with the CDS's own
# function. Recorded per arrow (column mobile) so the caption can list the grey ones. Matched
# on CDS ends within 60 bp, as the verifier matches them.
def _st492_mobile(sample, contig, start, end):
    _st492_hits(sample, contig, 0, 0)                      # fills the per-contig cache
    c = _ST492_MH.get((sample, contig))
    if c is None: return ""
    lo = c[["Start", "Stop"]].min(axis=1).astype(int); hi = c[["Start", "Stop"]].max(axis=1).astype(int)
    c = c[c["DB"].astype(str).isin(_ST492_MOBILE_DBS) & ((lo - int(start)).abs() < 60)
          & ((hi - int(end)).abs() < 60)]
    if not len(c): return ""
    best = c.sort_values(["Identity", "Subject_coverage"], ascending=False).iloc[0]
    sid = str(best["SubjectIdentifier"])
    if str(best["DB"]) == "ISFinder-Prot_MOBILE": return "ISFinder " + sid.split("//")[-1]
    f = sid.split("|")
    return f"{f[0]} {f[3].rstrip(',;')}" if len(f) >= 4 else sid

def _st492_label_rules(goi, nb):
    """Names, colours and notes of one track's Bakta CDS: every rule except those that add or
    re-read arrows from the assembly and the contig-end primes (_st492_display_labels)."""
    nb = nb.copy(); out, cats, notes = [], [], []
    sample, contig = nb["rep_sample"].iloc[0], nb["contig"].iloc[0]
    for _, r in nb.iterrows():
        lab, cat = r["label"], r["category"]
        p = _unquote(str(r["product"])).replace("%2C", ",").strip(); pl = p.lower()
        gene = str(r.get("gene", "") or "")
        note = ""
        if cat == "MGE":
            if lab == "intI1" and not _re.search(r"intI1|class 1 integron", p, _re.I):
                # the builder's "integron integrase" key also catches chromosomal tyrosine
                # recombinases; only a class 1 integrase is printed as intI1
                lab = gene if _ST492_GENE_OK.fullmatch(gene) else "integrase"
            elif "resolvase" in pl and gene == "tnpR":
                lab = "tnpR"                    # "Transposon Tn21 resolvase" is tnpR, not "Tn21"
            else:
                lab, note = _st492_is_name(sample, contig, r, lab, p, nb)
            out.append(lab); cats.append(cat); notes.append(note); continue
        if cat != "other": out.append(lab); cats.append(cat); notes.append(note); continue
        length = int(r["end"]) - int(r["start"]) + 1
        new = lab
        h = _st492_hits(sample, contig, int(r["start"]), int(r["end"]), ("CGE-MGE-Prot_MOBILE",))
        tn = [m.group(1) for s in h["SubjectIdentifier"].astype(str)
              for m in [_re.search(r"\|(Tn[pi][A-Z])\|", s)] if m]
        # A regulator Bakta names only by its family ("DNA-binding transcriptional regulator,
        # Lrp family", "DNA-binding protein, AcrR family, ...") prints as the family ("Lrp fam.",
        # "AcrR fam.") on every locus, even where Bakta's gene= gives the family's type member
        # (lrp, acrR): the gene symbol would overstate the annotation (user decision of
        # 2026-10-07). Real gene names from gene= (ycbK, hdeD, ...) are untouched.
        rfam = _re.match(r"DNA-binding (?:transcriptional regulator|protein), (\w+) family\b", p)
        if rfam:
            new = f"{rfam.group(1)} fam."; note = "regulator family"
        elif tn and lab in ("eAL", "EAL", "GGDEF", "GGDEF-EAL"):
            # a Tn21-type tnpM between intI1 and the Tn21 resolvase, which Bakta names by its EAL
            # domain: the CGE-MGE call names the gene, and a transposition-module gene is an MGE
            new = tn[0][0].lower() + tn[0][1:]; cat = "MGE"
        elif lab == "tonB" and _re.search(r"tonb[- ]dependent", pl) and "znud2" not in pl:
            # a TonB-dependent receptor piece, not the TonB protein the "tonb" key meant
            new = "TBDR frag." if length < 1500 else "TBDR"
        elif lab == "eAL":
            e = bool(_re.search(r"\beal\b", pl)); gg = "ggdef" in pl or "diguanylate" in pl
            new = "GGDEF-EAL" if (e and gg) else ("EAL" if e else "GGDEF")
        elif _re.search(r"\b(domain|motif)\b|DUF\d", p) and _ST492_GENE_OK.fullmatch(gene):
            new = gene                                 # rule (e): ycbK, not "DUF882"
        elif lab.endswith(" DCP"):
            # the bare domain name; an undefined " DCP" suffix helped nobody
            m = _re.search(r"DUF\d+", p)
            new = (m.group(0) if m else "DUF") if "duf" in pl else p.split(" domain")[0].strip()
        elif _re.fullmatch(r"NF\d+ family protein", p):
            new = ""                                   # an NCBIfam accession, not a name
        elif goi == "blaOXA-72" and gene and _re.search(r"\b(domain|motif)\b", pl):
            # Bakta's gene= here is a lower-cased DOMAIN abbreviation ("RNA recognition motif
            # (RRM) domain" carries gene=rRM), not a gene symbol. Print the domain name as the
            # " DCP" rule above does (HEPN), never a pseudo-gene. "DNA-binding ..., X family"
            # regulators are handled for every locus above ("Lrp fam.", "AcrR fam.").
            m = _re.search(rf"\b({_re.escape(gene)})\b", p, _re.I)
            if m and m.group(1).isupper(): new = m.group(1)
        else:
            g = _st492_gene_name(p, gene, lab)             # rule (e)
            if g: new = g
            # a bare cofactor word with no gene symbol is left blank rather than printed
            elif lab in ("NAD", "NADH", "NADP", "FAD", "FMN"): new = ""
        out.append(new); cats.append(cat); notes.append(note)
    nb["lab0"] = nb["label"]; nb["label"] = out; nb["category"] = cats; nb["note"] = notes

    # (b) and (d): an AMR call at >= 80 percent subject coverage makes a grey CDS an ARG, named by
    # the call; an aph(3') ARG takes its subclass from the call
    for i, r in nb.iterrows():
        if r["category"] not in ("other", "HP") and not (r["category"] == "ARG" and r["label"] == "aph(3')"):
            continue
        h0 = _st492_hits(sample, contig, int(r["start"]), int(r["end"]), _ST492_AMR_DBS)
        h0 = h0[h0["Subject_coverage"].astype(float) >= 80]
        h = _st492_acquired(h0)
        if len(h0) and not len(h):
            # only a resistance MUTATION call (parC on oprD): stays grey, with the reason the
            # caption gives (user decision of 2026-10-07: only acquired-gene calls recolour)
            nb.at[i, "note"] = "kept grey: resistance point-mutation call"; continue
        if not len(h): continue
        name = _st492_amr_name(h)
        if name in _ST492_KEEP_GREY:
            nb.at[i, "note"] = "kept grey: " + _ST492_KEEP_GREY[name]; continue
        # a call on a paralogue that the focal-copy search rejected is not a second focal copy
        if r["category"] != "ARG" and name == goi: continue
        nb.at[i, "category"] = "ARG"
        if name: nb.at[i, "label"] = name
        if r["category"] != "ARG": nb.at[i, "note"] = "AMR call"

    # A grey CDS with a mobile-element call (CGE-MGE or ISFinder, >= 80 percent subject coverage)
    # stays grey here, but the note records the call so the caption must say why (guide 2.6):
    # the MGE key reads "tnp / IS / integrase", and a plasmid replication initiator such as
    # repM (CGE-MGE RepAciN, on every AC237 track of blaNDM-1) is none of those. Whether such
    # calls should be coloured MGE instead is open (guide 2.7).
    for i, r in nb.iterrows():
        if r["category"] not in ("other", "HP") or str(nb.at[i, "note"]).startswith("kept grey"):
            continue
        h = _st492_hits(sample, contig, int(r["start"]), int(r["end"]), _ST492_MOBILE_DBS)
        h = h[h["Subject_coverage"].astype(float) >= 80]
        if len(h):
            nb.at[i, "note"] = "; ".join(x for x in (str(nb.at[i, "note"] or ""), "MOBILE call kept grey") if x)
    return nb

def _st492_display_labels(goi, nb, rep):
    nb = _st492_label_rules(goi, nb)
    nb["gap_x0"] = float("nan"); nb["gap_x1"] = float("nan")
    sample, contig = nb["rep_sample"].iloc[0], nb["contig"].iloc[0]
    nb = _st492_renames(nb, sample, contig)
    nb = _st492_fusion(nb, sample, contig)                          # (a)
    nb = _st492_ambiguity(nb, sample, contig)
    nb = _st492_end_fragment(nb, sample, contig)
    pool = _st492_pool(goi, rep)
    nb = _st492_gap_split(nb, sample, contig, pool)
    nb["lab0"] = nb["lab0"].fillna(nb["label"])
    nb = _st492_contained_primes(nb, sample, contig, pool)          # (f), by sequence
    # commonest drawn length per name (ties to the longer): sul1 has a real 885 bp variant
    lens = defaultdict(Counter)
    for _s, rnb in rep.values():
        for _, rr in rnb.iterrows():
            lens[rr["label"]][int(rr["end"]) - int(rr["start"]) + 1] += 1
    maxlen = {k: max(c.items(), key=lambda t: (t[1], t[0]))[0] for k, c in lens.items()}
    nb = _st492_contig_end_primes(nb, sample, contig, maxlen)       # (f), by codons
    nb = _st492_track_cols(nb, sample, contig)
    # the focal copy too: the caption's list of names with a CGE-MGE match must be complete
    nb["mobile"] = [_st492_mobile(sample, contig, r["start"], r["end"]) for _, r in nb.iterrows()]
    # A second copy of the gene of interest inside a window (blaOXA-72 csv v4) is drawn in the
    # focal red as the gene of interest (user decision of 2026-10-06), so the table says so too
    # instead of the "ARG" the builder's label_cat gave it.
    nb.loc[(nb["label"] == goi) & (nb["category"] != "center"), "category"] = "goi_copy"
    return nb
# ==== end of allctx corrections ============================================================

'''


def patch(src, pattern, repl, what, count=1):
    """re.subn that must hit exactly `count` times, so a reshaped builder stops the run."""
    new, k = re.subn(pattern, repl, src)
    if k != count:
        sys.exit(f"{what}: expected {count} match(es) in the base builder, found {k}; "
                 "the base builder has changed shape")
    return new


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None,
                    help="directory for the two _allctx8 CSVs (default: the builder's own)")
    args = ap.parse_args()
    if not SRC.exists():
        sys.exit(f"cannot find the base builder at {SRC}")
    src = SRC.read_text()

    # 1. lift the per-locus cap, anchored to the assignment
    src = patch(src, r"(?m)^MAX_VARIANTS\s*=\s*\d+", "MAX_VARIANTS = 10**6", "MAX_VARIANTS")

    # 2. redirect the output. Anchored on the OUT assignment itself: the basename also occurs
    #    in the docstring, which an unanchored count used to mistake for a second output.
    m = re.search(r'(?m)^OUT\s*=\s*"([^"]*/)st492_synteny_variants\.csv"', src)
    if not m:
        sys.exit("OUT assignment not found: the base builder has changed shape")
    outdir = (args.outdir.rstrip("/") + "/") if args.outdir else m.group(1)
    out = f"{outdir}st492_synteny_variants{SUFFIX}.csv"
    src = patch(src, r'(?m)^OUT\s*=\s*"[^"]*st492_synteny_variants\.csv"', f'OUT = "{out}"', "OUT")

    # 3. the corrections, each anchored and counted
    src = patch(src, r"(?m)^(GOIS = \[)", PRELUDE.replace("\\", "\\\\") + r"\1", "prelude anchor")
    src = patch(src, r"classify\(mob_call\(k\[0\], contig_of\[k\]\), k\[0\]\)",
                "_st492_classify(goi, k[0], contig_of[k], classify)", "classify call sites", 3)
    src = patch(src, r'(?m)^(    # mob_recon molecule type\. A "chromosome" call)',
                r"    groups, rep = _st492_regroup(goi, groups, rep)\n\n\1", "regroup anchor")
    src = patch(src, r'(?m)^        ctxt = ", "\.join\(f"\{c\} \{n\}" for c, n in cc\.most_common\(3\)\)$',
                "        ctxt = _st492_country_text(cc)\n"
                "        _st492_dump(goi, rank, kset, contig_of, trunc_of, parent.get(rank, []), classify)",
                "country line")
    src = patch(src, r'(?m)^        n_lab = f"n = \{len\(kset\)\} copies" \+ \(f" in \{len\(genomes\)\} genomes" '
                     r'if len\(genomes\) != len\(kset\) else ""\)$',
                "        n_lab = _st492_n_lab(len(kset), len(genomes))", "n_lab")
    # the track location text: strict-majority cluster rule (_st492_location), replacing the
    # builder's block from the location counts to the ", k chromosomal" suffix
    src = patch(src, r'(?ms)^(        locs, n_pl, n_ch, n_un = loc_counts\(kset\)\n)'
                     r'        resolved = n_pl \+ n_ch\n.*?'
                     r'^        if n_ch and molecule == "plasmid": loc \+= f", \{n_ch\} chromosomal"\n',
                r"\1        molecule, cluster, loc = _st492_location(locs, len(kset))\n", "location block")
    src = patch(src, r"(?m)^(        frames\.append\(nb\[\[)",
                r"        nb = _st492_display_labels(goi, nb, rep)\n\1", "frames.append")
    # the per-arrow note, mobile match and representative CDS coordinates (so the verifier can
    # check a drawn colour against marker_hits), the N-gap position of a gene split by a
    # scaffold gap, and the per-track gap and contig-end columns go into the table
    src = patch(src, r'"x0", "x1", "orient", "label", "category"\]\]\)',
                '"x0", "x1", "orient", "label", "category", "contig", "start", "end", "note", '
                '"mobile", "gap_x0", "gap_x1", "gaps", "contig_end"]])',
                "frames.append columns")
    src = patch(src, r"(?m)^        if hits: parent\[i\] = hits$",
                "        hits = _st492_parents_ok(goi, i, hits, ranked, loc_counts, [j for j, (s_j, _k) "
                "in enumerate(ranked, 1) if j != i and sublist(li, lab_of[s_j])])\n"
                "        if hits: parent[i] = hits", "fragment parents")
    src = patch(src, r"(?m)^(allg\.to_csv\(OUT, index=False\))$",
                r"allg = _st492_same_genome_marks(allg)\n\1\n_st492_write_copies(OUT)", "variants to_csv")

    # runpy needs a real file so the builder's own __file__-relative paths still resolve
    with tempfile.NamedTemporaryFile("w", suffix=".py", dir=str(SRC.parent),
                                     delete=False) as fh:
        fh.write(src)
        tmp = pathlib.Path(fh.name)
    try:
        print(f"[allctx] running the patched builder (no per-locus cap, output {out})",
              file=sys.stderr)
        runpy.run_path(str(tmp), run_name="__main__")
    finally:
        tmp.unlink(missing_ok=True)

    # 4. both tables must exist and agree: copies per locus in the per-copy table equal the
    #    summed n of the variant tracks
    import pandas as pd
    copies_p = pathlib.Path(out.replace("st492_synteny_variants", "st492_synteny_contexts"))
    if not pathlib.Path(out).exists() or not copies_p.exists():
        sys.exit("an output table was not written")
    v = pd.read_csv(out); c = pd.read_csv(copies_p)
    nv = v.drop_duplicates(["goi", "variant"]).groupby("goi")["n"].sum()
    nc = c.groupby("goi").size()
    if not nv.sort_index().equals(nc.reindex(nv.index).sort_index()):
        sys.exit(f"per-copy table disagrees with the variant table:\n{nv}\n{nc}")
    print(f"[allctx] per-copy table agrees with the variant table "
          f"({int(nc.sum())} copies over {len(nc)} loci)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
