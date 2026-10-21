#!/usr/bin/env python3
"""Verify the per-locus all-contexts synteny figures against their tables AND against upstream.

Every render check is made on the RENDERED files, not on the plotting parameters, because the
parameters are what is being tested. For each locus:

  tracks     the number of "variant" rows in the PDF text layer equals the number of distinct
             variants in the CSV, so no context was silently dropped or duplicated
  labels     the figure has no subtitle (since 2026-10-06), so every count is read off the
             DRAWN track labels: one "n = K cop... (countries)" label per track, display
             numbers 1..N, each label identical to its row of the layout CSV and of the
             variants CSV; the copies sum to the CSV's summed n, the genomes ("in G genomes",
             else K) sum to the per-track genome count of the per-copy table, and each
             track's printed country tallies (plus "+K other") sum to that track's genomes
  countries  a drawn "+K other" appears only on a track with more than seven countries, and
             every country never hidden in a "+K other" sums, over the drawn labels, to its
             metadata count on the locus
  chrome     no subtitle and no notes block drawn (unless ST492_SUBTITLE / ST492_NOTES = 1,
             the plot script's own switches), and no separate "additional copy" key or pink
             fill: a second focal copy is drawn in the focal red
  grammar    no "n = 1 copies" is printed
  scale      the drawn scale bar, measured from its vector path, matches its printed label,
             and matches across the whole set
  bleed      no ink touches the page edge
  page       the PDF page size equals the PNG pixel size * 72 / --dpi, within 1 pt per side

  added 2026-10-06 (rule k): per-country tallies compared per track (a genome with copies on
             two tracks is named on both); letters and cluster text against the strict-majority
             rule; no text between title and legend (any subtitle wording); legend keys only
             with a referent; DejaVu fonts only; arrow fills per category and the focal red on
             every copy of the gene; AMR-called CDS drawn as ARG (marker_hits); no glyph ink
             within OVERPRINT_PT of a gene arrow (vector outlines); no drawn arrow across an
             N run unless flagged; the caption file present, fresh and matching the figure;
             the named intI1 two-copy genome on two tracks, each marked with the other; and a
             missing layout CSV fails instead of being skipped
  added in the intI1 repair round: each drawn count label equals its variants-csv track text
             exactly (a same-genome mark compared through the layout map); a bleMBL::qacEdelta1
             read-through ORF found after a blaNDM-1 arrow (independently, from the FNA) must be
             drawn, translate to the fusion and abut blaNDM-1, with no Bakta ORF left over it;
             an ARG arrow with no Bakta CDS behind it must carry a note and be explained in the
             caption; a grey arrow with a MOBILE match must not be a mobility protein and must
             be listed in the caption; every contig the recount drops is named in the caption;
             a track drawn as a fragment is named as one in the caption
  added in the blaNDM-1 fix round (merged 2026-10-06): read from the vector PDF
             (pdftocairo -svg), no two glyphs overlapping, the dashed focal line through no
             glyph, no arrow drawn over a gene-of-interest arrow, arrow fills from the palette
             only, every legend key square with an arrow in its colour and every drawn colour
             keyed, the gene-category keys on one row; each track's location text equal to the
             one its copies give and each copy's MOB-suite call equal to contig_report.txt; the
             caption explains every primed, "fam." or read-from-the-assembly name and every
             printed abbreviation (DNA BP, CHAP, ...); no en or em dash in the figure or the
             caption (--caption-dir says where the caption files are, default --dir)
  added in the repair round after the merge: outlines of two different gene names at least
             GLYPH_CLEAR_PT apart (the hull test let "tolQ" touch "exbD"); every Bakta CDS inside
             a drawn span drawn; a copy on any closed-circle contig never truncated, and the
             caption's truncated count equal to the table's; a page of whole points, its size
             stated exactly in the caption; no other-ARG colour sentence without an other-ARG
             arrow, an MGE sentence when MGE arrows are drawn; no coloured name among the grey
             names left out; a primed name printed nowhere is said to be so; and across the set,
             equal title, header and item text heights, scale-bar ticks, bar-to-label gap and
             outer margins
  added in the stacking round (2026-10-07): the drawn order and stacks equal the explicit
             STACKS table (a copy of the plot script's), stacked rows at the tight pitch unless
             the upper row prints names of its own, every CDS of a stacked track has a
             counterpart in its root and its differing names are exactly the table's and are
             printed; blaOXA-72 keeps the display numbers of the 2026-10-06 render; the decided
             totals (blaOXA-72 84 copies in 82 genomes, blaNDM-1 17 in 17, intI1 63 in 62),
             drawn, in the per-copy table and upstream, the upstream DISTINCT genome count now
             checked also when one genome has copies on two tracks (intI1); circle facts equal
             for one contig under every locus, the named seamless circles closed and never
             truncated; the caption spells the fusion one way, its MGE sentence names only what
             is drawn and IS91/ISCR1 is glossed. The glyph-hull overlap test skips two glyphs
             drawn one after the other (a kerned pair such as "Ly"); names are tested against
             each other on their outlines.
  added in its repair round (2026-10-07): the stack root is the row next to the stack's
             names (first in the top stack, last in a lower one);
             a name difference counts only where one of the two names would be printed (grey
             names over 10 characters never are); the legend's MGE key names integrases only
             where an integrase is drawn as MGE.

Table invariants (variants CSV, per-copy CSV, layout CSV written beside each figure):

  n          n is constant within each (locus, variant) and equals that variant's copy rows
  countries  the printed country counts, plus any "+K other", sum to the track's genomes, and
             each printed count equals an independent metadata join on those genomes; a
             "+K other" never hides a country with as many genomes as one that is printed
             (the cut must not split a tie), and is used exactly when the track has more
             than seven countries (every country is named up to seven)
  grammar    no track row says "n = 1 copies" or "in 1 genomes"
  letters    each drawn C / P / U letter, and its cluster, follow the strict-majority rule
  fragments  no honoured fragment link joins two different molecules or an unresolved track,
             and a link is drawn only when exactly one parent qualifies (guide section 8)
  circles    a copy on a whole circle (exact end overlap or inexact terminal repeat) is not
             marked truncated

--upstream additionally recounts copies and genomes straight from marker_hits.parquet and the
Bakta GFF3/FNA (independently of the builder: own CDS matching, own 31-mer containment for
contig fragments) and fails if the copies summed over the drawn track labels, or the distinct
genomes of the per-copy table, disagree. carO has no marker hit and is skipped
there. It also maps every drawn arrow back to its Bakta CDS and fails on a CDS with an AMR
call drawn grey, aph(3') without its subclass, a family-only IS transposase printed as an
element without recorded evidence, a fusion name on a Bakta CDS, a mobile-element call drawn
grey with no recorded reason, "no ISFinder match" said of a transposase ISFinder does call, or
an arrow redrawn over ambiguity codes that spans none or still has a Bakta CDS drawn inside it,
and recounts the hits dropped as duplicate contigs against st492_synteny_dropped_allctx8.csv.
Copies on one contig are all counted, as the builder counts them: by the user's
decision of 2026-10-06 the two doubled AC163 contigs carry two counted blaOXA-72 copies each,
while a copy on a shorter contig that duplicates a longer one (the three Zurich genomes) is
still dropped. For blaOXA-72 it also fails if either doubled contig loses a copy, upstream,
in the per-copy table or in the drawn track, or if a Zurich genome counts more than one,
and it fails unless the drawn Swiss tallies sum to the seven Swiss blaOXA-72 genomes.

Exit status is non-zero if any check fails, so the driver can stop on a bad render.

Usage:  python3 src/scripts/verify_st492_synteny_allctx.py [--dpi 1200] [--dir DIR]
                                                            [--loci blaOXA-72,...] [--upstream]
                                                            [--csv --copies --caption-dir]
"""
import argparse
import collections
import csv
import math
import os
import pathlib
import re
import subprocess
import shutil
import sys
import tempfile

import numpy as np
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

LOCI = ["blaOXA-72", "blaNDM-1", "intI1", "blaOXA-66",
        "blaOXA-148", "blaOXA-23", "carO", "oprD"]
BAR_TOL_PCT = 2.0          # a bar is accepted within this much of its printed label
SPREAD_TOL_PCT = 2.0       # and the set's bars must agree within this much of each other
PAGE_TOL_PT = 1.0
ALL_COUNTRIES_MAX = 7      # the builder names every country up to this many per track
# The plot script's own switches; the verifier expects what the script was told to draw
SHOW_SUBTITLE = os.environ.get("ST492_SUBTITLE", "0") == "1"
SHOW_NOTES = os.environ.get("ST492_NOTES", "0") == "1"
PINK_RE = re.compile(r"rgb\(94\.11\d*%,\s*54\.90\d*%,\s*54\.90\d*%\)")  # #f08c8c

BK = "/home/bals/Documents/dev/BACTolomew/runs/sample_78eaa9126850/run_78eaa9126850/bakta"
MH = "/home/bals/Documents/dev/BACTolomew/cache/run_78eaa9126850/marker_hits.parquet"
MD = "/Data/ST492/metadata/metadata_ALL_new_merged_with_clusters_with_complex.csv"


def read_rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def expected_location(copies):
    """(letter, cluster) a track must print, from its copies' location calls (rule g): P with
    the cluster only when one cluster holds a strict majority of all its copies; C when
    chromosome calls outnumber plasmid ones and resolved calls are not outnumbered; else U."""
    locs = collections.Counter(c["location"] for c in copies)
    n_pl, n_ch = locs.get("plasmid", 0), locs.get("chromosome", 0)
    n = sum(locs.values()); n_un = n - n_pl - n_ch
    if n_pl + n_ch == 0 or n_un > n_pl + n_ch: return "U", ""
    if n_pl < n_ch: return "C", ""
    cl = collections.Counter(c["location_cluster"] for c in copies
                             if c["location"] == "plasmid" and c["location_cluster"] not in ("", "-"))
    top = sorted(cl.items(), key=lambda t: (-t[1], t[0]))[:1]
    if not top or 2 * top[0][1] <= n: return "U", ""
    return "P", top[0][0]


def table_checks(g, vrows, crows, lay, country_of):
    """Invariants of the tables the figure was drawn from. Returns (expected, problems)."""
    bad = []
    vr = [r for r in vrows if r["goi"] == g]
    cr = [r for r in crows if r["goi"] == g]
    n_of, track_of = {}, {}
    for r in vr:
        v = int(r["variant"]); n = int(r["n"])
        if n_of.setdefault(v, n) != n:
            bad.append(f"n differs between rows of csv v{v}")
        track_of[v] = r["track"]
    by_var = collections.defaultdict(list)
    for r in cr:
        by_var[int(r["variant"])].append(r)
    for v, n in n_of.items():
        if not cr: break                        # no per-copy table: nothing to join against
        if len(by_var.get(v, [])) != n:
            bad.append(f"csv v{v}: n = {n} but {len(by_var.get(v, []))} copy rows")
        genomes = sorted({r["sample"] for r in by_var.get(v, [])})
        m = re.search(r"\(([^()]*)\)\s*$", track_of[v])
        if not m:
            bad.append(f"csv v{v}: no country breakdown"); continue
        printed, rest = {}, 0
        for part in [p.strip() for p in m.group(1).split(";")[0].split(",")]:
            mo = re.fullmatch(r"\+(\d+) other", part)
            if mo: rest = int(mo.group(1)); continue
            name, k = part.rsplit(" ", 1); printed[name] = int(k)
        if sum(printed.values()) + rest != len(genomes):
            bad.append(f"csv v{v}: countries sum to {sum(printed.values()) + rest}, "
                       f"track has {len(genomes)} genomes")
        truth = collections.Counter(country_of.get(s, "unknown") for s in genomes)
        for name, k in printed.items():
            if truth.get(name, 0) != k:
                bad.append(f"csv v{v}: '{name} {k}' but metadata gives {truth.get(name, 0)}")
        if rest and len(truth) <= ALL_COUNTRIES_MAX:
            bad.append(f"csv v{v}: '+{rest} other' on a track with only {len(truth)} countries "
                       f"(every country is named up to {ALL_COUNTRIES_MAX})")
        if not rest and len(truth) > ALL_COUNTRIES_MAX:
            bad.append(f"csv v{v}: {len(truth)} countries all named, expected the top three "
                       f"and '+K other' above {ALL_COUNTRIES_MAX}")
        hidden = {c: k for c, k in truth.items() if c not in printed}
        if hidden and printed and max(hidden.values()) >= min(printed.values()):
            bad.append(f"csv v{v}: '+{rest} other' hides {hidden} although "
                       f"{min(printed.values())} is printed (a tie split at the cut)")
        if re.search(r"\bn = 1 copies\b|\bin 1 genomes\b", track_of[v]):
            bad.append(f"csv v{v}: singular count printed as plural in the track column")
    if lay:
        mol = {int(r["csv_variant"]): r["molecule"] for r in lay}
        for r in lay:
            v = int(r["csv_variant"])
            want, want_cl = expected_location(by_var.get(v, []))
            if r["letter"] != want:
                bad.append(f"display v{r['display']} / csv v{v}: letter {r['letter']}, "
                           f"copies say {want}")
            # rule (g): the cluster printed is the strict-majority cluster of ALL the track's
            # copies, with no fraction, and nothing else
            got_cl = (r["label"].rsplit("  ·  ", 1)[1] if "  ·  " in r["label"] else "")
            got_cl = got_cl.split(",")[0].strip()
            if got_cl != want_cl or re.search(r"\(\d+/\d+\)", r["label"]):
                bad.append(f"display v{r['display']} / csv v{v}: cluster text '{got_cl}', "
                           f"strict-majority rule gives '{want_cl}'")
            p = int(float(r["fragment_parent_csv"] or 0))
            if p and (mol.get(p) != r["molecule"] or r["molecule"] == "unresolved"):
                bad.append(f"display v{r['display']} / csv v{v}: fragment link to csv v{p} "
                           f"crosses molecules ({r['molecule']} -> {mol.get(p)})")
            # guide section 8: a fragment link only when exactly one parent qualifies
            cand = {str(c.get("candidate_parents") or "").strip() for c in by_var.get(v, [])}
            if p and cand != {str(p)}:
                bad.append(f"display v{r['display']} / csv v{v}: drawn as a fragment of csv v{p}, "
                           f"but the qualifying parents are {sorted(cand)}")
    for r in cr:
        # a circle's ends are where the assembler opened it, not a truncation: whether the
        # circle is drawn whole or the window only reaches its opening point (intI1 csv v15-17).
        # Tested on the sequence itself, not on the table's own flags.
        circ = int(float(r.get("end_overlap_bp") or 0)) or int(float(r.get("end_repeat_bp") or 0))
        if r.get("truncated") == "True" and (circ or (r["sample"], r["contig"]) in SEAMLESS
                                             or contig_is_circle(r["sample"], r["contig"])):
            bad.append(f"csv v{r['variant']} {r['sample']} {r['contig']}: a closed circle marked truncated")
    exp = dict(tracks=len(n_of), copies=sum(n_of.values()),
               genomes=len({r["sample"] for r in cr}) if cr else None,
               # genomes counted per track, as the drawn labels count them; equal to the
               # distinct total unless one genome has copies on two tracks
               track_genomes=(sum(len({r["sample"] for r in by_var.get(v, [])}) for v in n_of)
                              if cr else None),
               truth=collections.Counter(country_of.get(s, "unknown")
                                         for s in {r["sample"] for r in cr}),
               # countries counted per track, as the drawn labels count them: a genome with
               # copies on two tracks is tallied on both (F20)
               track_truth=sum((collections.Counter(country_of.get(s, "unknown")
                                                    for s in {r["sample"] for r in by_var.get(v, [])})
                                for v in n_of), collections.Counter()),
               csv_track={v: t for v, t in track_of.items()},
               truncated=sum(1 for r in cr if r.get("truncated") == "True") if cr else None,
               # drawn colours: a second copy of the gene (ARG in the table) is drawn in the focal red
               cats={r["category"] for r in vr if r["label"] != g},
               coloured_names={r["label"] for r in vr if r["category"] not in ("other", "HP")})
    if cr and len(cr) != exp["copies"]:
        bad.append(f"per-copy table has {len(cr)} rows, CSV n sums to {exp['copies']}")
    return exp, bad


# ---- independent upstream recount ---------------------------------------------------------
_fna = {}
def contig_is_circle(sample, contig):
    """The contig's two ends overlap: exactly (>= 20 bp, the assembler's terminal overlap) or
    as a repeat of >= 100 bp with at most 1 percent mismatches (no indels: a seed-and-extend
    scan, independent of the builder's banded edit distance; an indel there is caught by the
    builder's own test and its contig_end column, which the caption checks read)."""
    q = fna(sample).get(contig, "")
    if len(q) < 200: return False
    for k in range(min(2000, len(q) // 2), 19, -1):
        if q[:k] == q[-k:]: return True
    head = q[:30]
    for i in range(max(len(q) - 2000, len(q) // 2), len(q) - 100):
        if q[i:i + 30] == head:
            tail = q[i:]
            if sum(a != b for a, b in zip(tail, q)) <= 0.01 * len(tail): return True
    return False

def fna(sample):
    if sample not in _fna:
        d, name = {}, None
        with open(f"{BK}/{sample}/{sample}.fna") as fh:
            for line in fh:
                if line.startswith(">"): name = line[1:].split()[0]; d[name] = []
                else: d[name].append(line.strip())
        _fna[sample] = {k: "".join(v).upper() for k, v in d.items()}
    return _fna[sample]


def cds(sample):
    out = collections.defaultdict(list)
    with open(f"{BK}/{sample}/{sample}.gff3") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 8 and f[2] == "CDS":
                out[f[0]].append((int(f[3]), int(f[4]), f[8]))
    return out


_RC = str.maketrans("ACGTN", "TGCAN")
def kmers(q, k=31, step=1):
    s = set(q[i:i + k] for i in range(0, len(q) - k + 1, step))
    return s


def contained(a, b):
    """Fraction of contig a's 31-mers (sampled) present in contig b, either strand."""
    ka = kmers(a, step=7)
    if not ka: return 0.0
    kb = kmers(b) | kmers(b.translate(_RC)[::-1])
    return len(ka & kb) / len(ka)


UP_DROPPED = {}       # goi -> {(sample, contig)} the recount dropped as contained duplicates
def upstream_counts(goi, st492):
    import pandas as pd
    mh = pd.read_parquet(MH, columns=["Sample", "Contig", "Start", "Stop", "Gene", "DB"])
    h = mh[(mh["Gene"].astype(str) == goi) & mh["Sample"].isin(st492)
           & (mh["DB"].astype(str) != "COGorg24_COG")]
    allele = goi[3:] if goi.startswith("bla") else None
    copies, raw = {}, {}
    for s, hs in h.groupby("Sample"):
        if not os.path.exists(f"{BK}/{s}/{s}.gff3"): continue
        cd = cds(s); found = set()
        for _, r in hs.iterrows():
            contig = f"contig_{r['Contig']}"; lo = min(int(r["Start"]), int(r["Stop"]))
            for a, b, attr in cd.get(contig, []):
                if abs(a - lo) < 60 and (allele is None or re.search(rf"\b{re.escape(allele)}\b", attr)):
                    found.add((contig, a)); break
        # Every copy on one contig counts, including the two copies 9,830 bp apart on each
        # doubled AC163 contig: the user decided both are counted, and an assembly cannot
        # tell a tandem dimer from an artefact, so the recount does not try to.
        keep = sorted(found)
        seqs = fna(s)
        # a copy on a shorter contig that is >= 90 percent contained in a longer contig
        # carrying a kept copy is the same locus caught twice
        keep.sort(key=lambda c: -len(seqs[c[0]]))
        final = []
        for c in keep:
            if any(f[0] != c[0] and len(seqs[f[0]]) >= len(seqs[c[0]])
                   and contained(seqs[c[0]], seqs[f[0]]) >= 0.9 for f in final):
                UP_DROPPED.setdefault(goi, set()).add((s, c[0]))
                continue
            final.append(c)
        if final: copies[s] = final
        if found: raw[s] = sorted(found)
    return sum(len(v) for v in copies.values()), len(copies), copies, raw


# The user's decisions on the blaOXA-72 multi-hit genomes, checked by name so that a change
# anywhere upstream or in the builder that undoes either one fails loudly.
DOUBLED = {("GCA_964331335.1", "contig_5"): (5885, 15715),     # two counted copies each,
           ("GCA_964331905.1", "contig_5"): (6010, 15840)}     # 9,830 bp apart, own track
ZURICH = ("GCA_031875465.1", "GCA_031875745.1", "GCA_031876255.1")  # one copy each


SWISS_GENOMES = 7          # blaOXA-72 genomes from Switzerland, all to be visible on the figure


def doubled_checks(g, vrows, crows, up=None, labels=None):
    """blaOXA-72 only: both copies of each doubled contig counted and drawn, Zurich once,
    and the drawn Swiss tallies account for all seven Swiss genomes."""
    if g != "blaOXA-72": return []
    bad = []
    if labels is not None:
        ch = sum(v[2].get("Switzerland", 0) for v in labels.values())
        if ch != SWISS_GENOMES:
            bad.append(f"drawn Swiss tallies sum to {ch}, expected {SWISS_GENOMES}")
    cr = [r for r in crows if r["goi"] == g]
    vars_ = set()
    for (s, contig), starts in DOUBLED.items():
        rows = sorted(int(r["cstart"]) for r in cr if r["sample"] == s and r["contig"] == contig)
        if tuple(rows) != starts:
            bad.append(f"{s} {contig}: per-copy table holds copies at {rows}, expected {list(starts)}")
        vars_ |= {r["variant"] for r in cr if r["sample"] == s and r["contig"] == contig}
        if up is not None:
            got = sorted(a for c, a in up.get(s, []) if c == contig)
            if tuple(got) != starts:
                bad.append(f"{s} {contig}: upstream recount finds copies at {got}, "
                           f"expected {list(starts)}")
    if len(vars_) != 1:
        bad.append(f"doubled contigs fall in csv variants {sorted(vars_)}, expected one track")
    else:
        v = vars_.pop()
        others = {r["sample"] for r in cr if r["variant"] == v} - {s for s, _ in DOUBLED}
        if others: bad.append(f"doubled-contig track csv v{v} also holds {sorted(others)}")
        arrows = [r for r in vrows if r["goi"] == g and r["variant"] == v and r["label"] == g]
        xs = sorted((float(r["x0"]) + float(r["x1"])) / 2 for r in arrows)
        if (len(arrows) != 2 or {r["category"] for r in arrows} != {"center", "goi_copy"}
                or abs(abs(xs[1] - xs[0]) - 9830) > 1):
            bad.append(f"doubled-contig track csv v{v} draws {len(arrows)} {g} arrow(s) "
                       f"at {xs}, expected the focal copy and one more 9,830 bp away")
    for s in ZURICH:
        k = sum(1 for r in cr if r["sample"] == s)
        if k != 1: bad.append(f"{s}: {k} copies in the per-copy table, expected 1")
        if up is not None and len(up.get(s, [])) != 1:
            bad.append(f"{s}: upstream recount keeps {len(up.get(s, []))} copies, expected 1")
    return bad


# ---- checks added 2026-10-06 (rule k) -----------------------------------------------------
# The plot script's palette and legend wording; a drawn fill outside it, or a key without an
# arrow behind it, fails.
PAL = {"center": "#c1121f", "goi_copy": "#c1121f", "ARG": "#7048e8", "MGE": "#ff8f00",
       "TA": "#0c8599", "porin": "#1971c2", "HP": "#adb5bd", "other": "#d7dde3"}
KEYS = {"center": "gene of interest", "ARG": "other ARG", "MGE": "MGE (tnp / IS",
        "TA": "toxin-antitoxin", "porin": "porin", "HP": "hypothetical", "other": "other CDS"}
LOC_KEYS = {"C": "chromosome", "P": "plasmid", "U": "location unresolved"}
AMR_DB = {"AMRF_core-AMR", "AMRF_plus-AMR", "CARD_AMR"}
KEEP_GREY = {"cxpE", "merA", "merR"}          # the wrapper's _ST492_KEEP_GREY
ARROW_STROKE = "stroke:rgb(35.29"              # #5a6572, the arrow outline colour
OVERPRINT_PT = 0.5                             # ink closer than this to an arrow fill is overprint
WORD_RE = re.compile(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">([^<]*)</word>')


def svg_rgb(style):
    m = re.search(r"fill:rgb\(([\d.]+)%,\s*([\d.]+)%,\s*([\d.]+)%\)", style)
    return "#%02x%02x%02x" % tuple(round(float(x) * 2.55) for x in m.groups()) if m else None


def svg_path_points(d, step=4):
    """Points along an SVG path (M/L/C/Z, absolute, as cairo writes them); cubic curves are
    sampled, never their control points, which can lie outside the outline."""
    toks = re.findall(r"[MLCZ]|-?[\d.]+(?:e-?\d+)?", d)
    pts, cur, i, cmd = [], (0.0, 0.0), 0, None
    while i < len(toks):
        t = toks[i]
        if t in "MLCZ": cmd = t; i += 1; continue
        if cmd in ("M", "L"):
            cur = (float(toks[i]), float(toks[i + 1])); pts.append(cur); i += 2
        elif cmd == "C":
            c = [float(x) for x in toks[i:i + 6]]; p0 = cur
            for k in range(1, step + 1):
                u = k / step; a, b, cc, e = (1 - u) ** 3, 3 * (1 - u) ** 2 * u, 3 * (1 - u) * u ** 2, u ** 3
                pts.append((a * p0[0] + b * c[0] + cc * c[2] + e * c[4],
                            a * p0[1] + b * c[1] + cc * c[3] + e * c[5]))
            cur = (c[4], c[5]); i += 6
        else:
            i += 1
    return pts


def svg_parts(svg_text):
    """(arrows, glyph ink points): arrows as (fill hex, polygon) for every path with the arrow
    outline colour, glyphs as point lists placed where each <use> puts them."""
    sym = {m.group(1): svg_path_points(m.group(2)) for m in re.finditer(
        r'<symbol[^>]*id="([^"]+)"[^>]*>\s*<path[^>]*d="([^"]*)"', svg_text)}
    glyphs = [[(x + float(u.group(2)), y + float(u.group(3))) for x, y in sym.get(u.group(1), [])]
              for u in re.finditer(r'<use xlink:href="#([^"]+)" x="([-\d.]+)" y="([-\d.]+)"', svg_text)]
    arrows = [(svg_rgb(st), svg_path_points(d)) for st, d in SEG_RE.findall(svg_text)
              if ARROW_STROKE in st and "fill:rgb" in st]
    return arrows, [g for g in glyphs if g]


def _inside(pt, poly):
    x, y = pt; c = False
    for (x1, y1), (x2, y2) in zip(poly, poly[1:] + poly[:1]):
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1 + 1e-12) + x1: c = not c
    return c


def _seg_dist(pt, a, b):
    (px, py), (ax, ay), (bx, by) = pt, a, b
    dx, dy = bx - ax, by - ay
    t = max(0.0, min(1.0, ((px - ax) * dx + (py - ay) * dy) / (dx * dx + dy * dy + 1e-12)))
    return ((px - ax - t * dx) ** 2 + (py - ay - t * dy) ** 2) ** 0.5


def overprint_checks(arrows, glyphs):
    """Rule (i): no glyph's ink within OVERPRINT_PT of any gene arrow (its fill polygon; the
    outline stroke is 0.34 pt, so this is about the stroke's outer edge)."""
    bad, polys = [], []
    for _f, pts in arrows:
        if len(pts) < 3: continue
        xs, ys = [p[0] for p in pts], [p[1] for p in pts]
        polys.append((min(xs) - 1, max(xs) + 1, min(ys) - 1, max(ys) + 1, pts))
    for g in glyphs:
        xs, ys = [p[0] for p in g], [p[1] for p in g]
        bx = (min(xs), max(xs), min(ys), max(ys))
        for x0, x1, y0, y1, poly in polys:
            if bx[1] < x0 or bx[0] > x1 or bx[3] < y0 or bx[2] > y1: continue
            if any(_inside(q, poly) or min(_seg_dist(q, a, b) for a, b in zip(poly, poly[1:] + poly[:1]))
                   < OVERPRINT_PT for q in g):
                bad.append(f"text overprints a gene arrow at x {bx[0]:.1f}, y {bx[2]:.1f} pt")
                break
    return bad[:10] + ([f"... {len(bad) - 10} more overprints"] if len(bad) > 10 else [])


def fill_checks(g, vrows, arrows):
    """Arrow fills against the palette: every drawn arrow has its category's fill, so the count
    per fill equals the table's count per category, and every copy of the gene of interest
    (a row named after it) is the focal red."""
    bad = []
    rows = [r for r in vrows if r["goi"] == g]
    want = collections.Counter(PAL["center"] if (r["category"] == "center" or r["label"] == g)
                               else PAL.get(r["category"], "?") for r in rows)
    got = collections.Counter(f for f, _ in arrows)
    for f in sorted(set(want) | set(got)):
        if want.get(f, 0) != got.get(f, 0):
            bad.append(f"{got.get(f, 0)} arrows filled {f}, the table's categories give {want.get(f, 0)}")
    if any(r["label"] == g and r["category"] not in ("center", "goi_copy", "ARG", "MGE") for r in rows):
        bad.append(f"a copy of {g} is drawn in a non-focal category")
    return bad


def legend_checks(g, words_xy, vrows, lay):
    """Every legend key has a referent on this locus, and every drawn category / letter has a
    key. Read from the PDF text between the 'Gene category' header and the first track."""
    bad = []
    txt = " ".join(w for *_xy, w in words_xy)
    rows = [r for r in vrows if r["goi"] == g]
    cats = {("center" if r["label"] == g else r["category"]) for r in rows}
    letters = {r["letter"] for r in lay}
    hdr = [y0 for x0, y0, x1, y1, w in words_xy if w == "Location"]
    top = min((y0 for x0, y0, x1, y1, w in words_xy if w == "variant"), default=1e9)
    leg = " ".join(w for x0, y0, x1, y1, w in words_xy if y0 < top and (not hdr or y0 <= max(hdr) + 2))
    for c, k in KEYS.items():
        if (k in leg) != (c in cats):
            bad.append(f"legend key '{k}' {'drawn without a referent' if k in leg else 'missing'}")
    # the MGE key names integrases only where one is drawn as MGE (blaNDM-1's intI1)
    has_int = any(r["category"] == "MGE" and re.match(r"intI|integrase", r["label"]) for r in rows)
    if "MGE" in cats and ("integrase" in leg) != has_int:
        bad.append(f"legend MGE key {'names' if not has_int else 'misses'} integrases"
                   f"{' with no MGE integrase drawn' if not has_int else ''}")
    for l, k in (LOC_KEYS.items() if lay else []):        # letters need the layout map
        if (re.search(rf"(?<!location ){re.escape(k)}" if k != "location unresolved" else k, leg) is not None) \
                != (l in letters):
            bad.append(f"location key '{k}' {'drawn without a referent' if l not in letters else 'missing'}")
    return bad


def font_checks(pdf):
    out = subprocess.run(["pdffonts", str(pdf)], capture_output=True, text=True).stdout.splitlines()[2:]
    names = [l.split()[0] for l in out if l.strip()]
    return [f"non-DejaVu font {n}" for n in names if "DejaVu" not in n]


def subtitle_checks(words_xy):
    """No subtitle: no text line between the title and the legend's 'Gene category' header.
    Position based, so any old subtitle wording is caught, not one phrasing."""
    bad = []
    hdr = [y0 for x0, y0, x1, y1, w in words_xy if w == "Gene"]
    if not hdr: bad.append("no 'Gene category' legend header found")
    # the legend starts at its header, or (old figures) at its first key
    keys = [y0 for x0, y0, x1, y1, w in words_xy if w in ("Gene", "interest", "ARG", "variant")]
    title_bottom = min((y1 for x0, y0, x1, y1, w in words_xy if w in ("ST492:", "baumannii")), default=0)
    mid = [w for x0, y0, x1, y1, w in words_xy if title_bottom + 1 < y0 < min(keys, default=0) - 1]
    return bad + ([f"subtitle drawn: {' '.join(mid[:12])}"] if mid else [])


_mh = None
def marker_checks(g, vrows):
    """ARG colour against marker_hits (rule b), independently of the wrapper: a drawn CDS with
    an AMRFinderPlus core/plus AMR or CARD AMR call at >= 80 percent coverage must be drawn as
    an ARG (or as the focal gene), unless it is a resistance MUTATION in a housekeeping gene or
    on the kept-grey list."""
    global _mh
    import pandas as pd
    rows = [r for r in vrows if r["goi"] == g]
    if not rows or "start" not in rows[0]:
        return ["variants table has no CDS coordinates (contig/start/end): rebuild with the wrapper"]
    if _mh is None:
        _mh = pd.read_parquet(MH, columns=["Sample", "Contig", "Start", "Stop", "DB", "Gene",
                                           "Subject_coverage", "ResistanceAlleles", "SubjectIdentifier"])
        _mh = _mh[_mh["DB"].astype(str).isin(AMR_DB) & (_mh["Subject_coverage"].astype(float) >= 80)
                  & ~_mh["ResistanceAlleles"].astype(str).str.fullmatch(r"\w+_[A-Z]\d+[A-Z*]")
                  & ~_mh["SubjectIdentifier"].astype(str).str.startswith("VM|")]
        _mh = _mh.assign(contig="contig_" + _mh["Contig"].astype(str),
                         lo=_mh[["Start", "Stop"]].min(axis=1).astype(int))
    bad = []
    for r in rows:
        if not r["start"] or r["category"] in ("center", "ARG") or r["label"] == g: continue
        if r["label"] in KEEP_GREY: continue
        h = _mh[(_mh["Sample"] == r["rep_sample"]) & (_mh["contig"] == r["contig"])
                & ((_mh["lo"] - int(float(r["start"]))).abs() < 60)]
        if len(h) and not set(h["Gene"].astype(str)) & KEEP_GREY:
            bad.append(f"csv v{r['variant']} {r['label'] or '(unlabelled)'} at {r['contig']} "
                       f"{int(float(r['start']))}: AMR call drawn as {r['category']}")
    return bad


def nrun_checks(g, vrows):
    """A drawn arrow that overlaps a run of >= 10 N in its representative contig must be one
    the wrapper flagged as split by that gap (it then carries gap_x0/gap_x1)."""
    bad = []
    for r in vrows:
        if r["goi"] != g or not r.get("start"): continue
        q = fna(r["rep_sample"]).get(r["contig"], "")[int(float(r["start"])) - 1:int(float(r["end"]))]
        if re.search(r"N{10,}", q) and not r.get("gap_x0"):
            bad.append(f"csv v{r['variant']} {r['label']}: arrow spans an N run but is not flagged")
    return bad


def caption_checks(g, cap, pdf, labels, exp, pw):
    """The caption file exists, was written by this render, states the drawn totals, the
    context count and the native width, cites only drawn display numbers, and has no dashes."""
    if not cap.exists(): return [f"caption file missing: {cap.name}"]
    bad, t = [], cap.read_text()
    if cap.stat().st_mtime < pathlib.Path(pdf).stat().st_mtime - 300:
        bad.append("caption file is older than the figure")
    n_k = sum(v[0] for v in labels.values())
    m = re.search(rf"(\d+) {re.escape(g)} copies in (\d+) genomes, (\d+) contexts", t)
    if not m:
        bad.append("caption does not state 'K copies in G genomes, N contexts'")
    else:
        if int(m.group(1)) != n_k: bad.append(f"caption says {m.group(1)} copies, figure {n_k}")
        if exp["genomes"] is not None and int(m.group(2)) != exp["genomes"]:
            bad.append(f"caption says {m.group(2)} genomes, per-copy table {exp['genomes']}")
        if int(m.group(3)) != exp["tracks"]:
            bad.append(f"caption says {m.group(3)} contexts, figure has {exp['tracks']}")
    # the page is a whole number of points, and the caption gives exactly that size
    pw, ph = page_pt(pdf)
    w = re.search(r"native size: ([\d.]+) mm wide by ([\d.]+) mm tall \((\d+) x (\d+) pt\)", t)
    if not w or (int(w.group(3)), int(w.group(4))) != (round(pw or 0), round(ph or 0)) \
            or abs(float(w.group(1)) - pw * PT2MM) > 0.06 or abs(float(w.group(2)) - ph * PT2MM) > 0.06:
        bad.append(f"caption size {w.group(0) if w else 'missing'}, page is {pw} x {ph} pt")
    if pw and ph and (abs(pw - round(pw)) > 0.01 or abs(ph - round(ph)) > 0.01):
        bad.append(f"page {pw} x {ph} pt is not whole points (cairo cuts the fraction from the PDF)")
    if exp.get("truncated") is not None:
        m = re.search(r"(\d+) of (\d+) cop(?:y|ies) l(?:ies|ie) in windows cut by a contig end", t)
        if not m or int(m.group(1)) != exp["truncated"]:
            bad.append(f"caption says {m.group(1) if m else 'nothing'} copies in windows cut by a contig "
                       f"end, the per-copy table {exp['truncated']}")
    cats = exp.get("cats", set())
    if "other ARG marks" in t and "ARG" not in cats:
        bad.append("caption explains the other-ARG colour, but no other-ARG arrow is drawn")
    if "MGE" in cats and "MGE marks" not in t:
        bad.append("MGE arrows drawn but the caption does not say what MGE marks")
    m = re.search(r"Names of grey CDS left out[^:]*: ([^\n]*)", t)
    for nm in (re.findall(r"(?:^|; )(.+?) \(variants? ", m.group(1)) if m else []):
        if nm in exp.get("coloured_names", set()):
            bad.append(f"caption lists the coloured name {nm} among the grey names left out")
    for v in re.findall(r"[Vv]ariants? ((?:\d+(?:, | and )?)+)", t):
        for x in re.findall(r"\d+", v):
            if not 1 <= int(x) <= exp["tracks"]: bad.append(f"caption cites variant {x}, not drawn")
    if re.search("[\u2013\u2014]", t): bad.append("caption contains an en- or em-dash")
    if t.count("Caption text for ") != 1:
        bad.append("caption file holds more than one render (append to a stale file?)")
    return bad


# intI1, user decision of 2026-10-06: the two complete copies on GCA_964331735.1 contig_1 are
# counted, on two tracks, each marked with the other's display number.
SEPARATE = {("intI1", "GCA_964331735.1", "contig_1"): (2160506, 2840056)}


def separate_checks(g, crows, lay, labels, up=None):
    bad = []
    for (gg, s, contig), starts in SEPARATE.items():
        if gg != g: continue
        rows = [r for r in crows if r["goi"] == g and r["sample"] == s and r["contig"] == contig]
        if tuple(sorted(int(r["cstart"]) for r in rows)) != starts:
            bad.append(f"{s} {contig}: per-copy table holds {sorted(int(r['cstart']) for r in rows)}, "
                       f"expected {list(starts)}")
        if up is not None and tuple(sorted(a for c, a in up.get(s, []) if c == contig)) != starts:
            bad.append(f"{s} {contig}: upstream recount disagrees with {list(starts)}")
        disp = {int(r["csv_variant"]): int(r["display"]) for r in lay}
        dv = sorted(disp.get(int(r["variant"]), -1) for r in rows)
        if len(set(dv)) != 2: bad.append(f"{s} {contig}: copies on displays {dv}, expected two tracks")
        for a, b in ((dv[0], dv[-1]), (dv[-1], dv[0])):
            if a in labels and labels[a][6] != [b]:
                bad.append(f"display v{a}: must say 'same genome as variant {b}', says {labels[a][6]}")
    return bad


# ---- added in the stacking round (2026-10-07) -----------------------------------------------
# The plot script's STACKS table, copied here on purpose: a change to one without the other
# fails. locus -> (order of csv variants top to bottom, {root: rows}, {csv variant: own names})
STACKS = {
    "blaOXA-72": (list(range(1, 12)), {4: [2, 3, 4], 11: [10, 11]},
                  {10: {"znuD2'"}}),
    "blaNDM-1": ([4, 6, 7, 8, 1, 5, 2, 9, 3], {4: [4, 6, 7, 8, 1, 5, 2, 9, 3]},
                 {7: {"IS15", "repM'"}, 8: {"IS15", "sul1'", "IS91/ISCR1'", "repM'"},
                  1: {"qacEΔ1'"}, 5: {"msr(E)'"}}),
    # the root is the row next to the stack's shared names: first in the top stack, last in a
    # lower one (repair round of 2026-10-07)
    "intI1": ([11, 1, 6, 20, 5, 10, 18, 14, 15, 17, 16, 8, 2, 3, 23, 4, 7, 9, 12, 13, 19, 21, 22],
              {11: [11, 1, 6, 20, 5, 10], 8: [18, 14, 15, 17, 16, 8], 4: [23, 4]},
              {16: {"IS15"}, 17: {"IS15", "IS91/ISCR1'", "sul1'"}}),
}
# blaOXA-72 keeps the display numbers of the 2026-10-06 render (user decision): display N is
# csv vN, with these representatives (backup layout csv of that render)
OXA72_REPS = ["20038ENE2", "GCA_016474335.1", "GCA_001400975.1", "GCA_964331335.1", "GCA_018928145.1",
              "GCA_958320075.1", "GCA_021570155.1", "GCA_021571095.1", "GCA_021571585.1",
              "GCA_958316555.1", "SAMN09202764"]
# copies and genomes per locus, fixed by the user decisions (guide sections 7 and 8)
EXPECTED_COUNTS = {"blaOXA-72": (84, 82), "blaNDM-1": (17, 17), "intI1": (63, 62)}
# the wrapper's _ST492_SEAMLESS: closed circles whose ends join without overlap, each a rotation
# of the reference circle with its terminal overlap trimmed (tested here on the sequence)
SEAMLESS_REF = ("GCA_029679805.1", "contig_57")
SEAMLESS = {("GCA_016474335.1", "contig_35"), ("GCA_016474395.1", "contig_35"),
            ("GCA_022017315.1", "contig_40"), ("GCA_022017345.1", "contig_39"),
            ("GCA_028577815.3", "contig_47"), ("GCA_048215925.1", "contig_42"),
            ("GCA_048245695.1", "contig_39"), ("GCA_048245935.1", "contig_44")}


MAXLAB = 10                                    # the plot script's length cap on grey names
ALWAYS = {"ARG", "goi_copy", "MGE", "TA", "porin"}


def stack_checks(g, vrows, lay, drawn_names, svg_path=None):
    """The drawn stacks equal the STACKS table: rows in its order, each stacked row tagged with
    its root in the layout csv and drawn at the tight pitch, rows of different stacks a full pitch
    apart; from the table, every CDS of a stacked track has a counterpart in its root (centre
    within 300 bp or inside it, same strand) and the names that differ are exactly its own names,
    which must be printed. blaOXA-72 keeps its earlier display numbers."""
    if g not in STACKS: return []
    order, stacks, extra = STACKS[g]
    bad = []
    by_disp = {int(r["display"]): r for r in lay}
    got = [int(by_disp[d]["csv_variant"]) for d in sorted(by_disp)]
    if got != order: bad.append(f"display order (csv) {got}, the stack table says {order}")
    root_of = {v: r for r, rows in stacks.items() for v in rows}
    # the root is the row next to the shared names: first in the top stack, last in a lower one
    for root, rows in stacks.items():
        top = order.index(rows[0]) == 0
        if root != (rows[0] if top else rows[-1]):
            bad.append(f"stack of csv v{root}: the root is not the row next to its names")
    for d, r in by_disp.items():
        want = root_of.get(int(r["csv_variant"]), 0)
        if int(float(r.get("stack_root_csv") or 0)) != want:
            bad.append(f"display v{d}: stack root csv {r.get('stack_root_csv')}, the table says {want}")
    # pitch, in y-units (one unit = the full row pitch): rows of one stack at the tight pitch
    # (0.36 at the full row pitch, more when IPU shrinks), unless the upper row prints names of
    # its own, which then widen the gap under it just enough; separate tracks a full pitch apart
    ys = [float(by_disp[d]["y"]) for d in sorted(by_disp)]
    tight = min((ys[i - 1] - ys[i] for i in range(1, len(ys))), default=1)
    for i in range(1, len(ys)):
        a, b = got[i - 1], got[i]
        same = a in root_of and root_of.get(a) == root_of.get(b)
        dy = ys[i - 1] - ys[i]
        if same and not (dy < 0.6 or (a in extra and dy < 1.5)):
            bad.append(f"displays {i} and {i + 1} share a stack but sit {dy:.2f} units apart")
        if same and a not in extra and dy > tight + 0.01:
            bad.append(f"displays {i} and {i + 1}: stacked rows {dy:.2f} units apart, the tight pitch is {tight:.2f}")
        if not same and dy < 0.999: bad.append(f"displays {i} and {i + 1} are not stacked but sit {dy:.2f} units apart")
    if g == "blaOXA-72":
        reps = [by_disp[d]["rep_sample"] for d in sorted(by_disp)]
        if got != list(range(1, 12)) or reps != OXA72_REPS:
            bad.append("blaOXA-72 display numbers differ from the 2026-10-06 render")
    rows = collections.defaultdict(list)
    for r in vrows:
        if r["goi"] == g and r["category"] != "center": rows[int(r["variant"])].append(r)
    def named(lb): return bool(lb) and lb != "HP"
    # a name a track would print: grey names over MAXLAB characters are printed on no track
    def visible(lb, cat): return named(lb) and (len(lb) <= MAXLAB or cat in ALWAYS)
    for root, members in stacks.items():
        rr = [(float(x["x0"]), float(x["x1"]), x["orient"], x["label"], x["category"]) for x in rows[root]]
        for v in members:
            if v == root: continue
            dif = set()
            for x in rows[v]:
                x0, x1 = float(x["x0"]), float(x["x1"])
                c = [t for t in rr if t[2] == x["orient"] and (abs((t[0] + t[1]) / 2 - (x0 + x1) / 2) <= 300
                                                               or (x0 >= t[0] - 300 and x1 <= t[1] + 300))]
                if not c:
                    bad.append(f"csv v{v} {x['label']} at {x0:g}: no counterpart in stack root csv v{root}"); continue
                t = min(c, key=lambda t: abs((t[0] + t[1]) / 2 - (x0 + x1) / 2))
                if named(t[3]) and not named(x["label"]):
                    bad.append(f"csv v{v}: unnamed where its root csv v{root} has {t[3]}")
                if named(x["label"]) and x["label"] != t[3] and (visible(x["label"], x["category"])
                                                                 or visible(t[3], t[4])): dif.add(x["label"])
            if dif != extra.get(v, set()):
                bad.append(f"csv v{v}: differs from root csv v{root} by {sorted(dif)}, "
                           f"the stack table lists {sorted(extra.get(v, set()))}")
    for v, names in extra.items():
        for nm in names:
            stem = nm.rstrip("'").replace(" ", "")
            if not any(stem in w for w in drawn_names):
                bad.append(f"csv v{v}: its own name {nm} is not printed")
    return bad


def count_checks(g, n_copies, crows, up=None):
    """The decided totals: drawn copies, the per-copy table's distinct genomes and, with
    --upstream, the independent recount (copies and DISTINCT genomes, also when one genome has
    copies on two tracks, as on intI1)."""
    if g not in EXPECTED_COUNTS: return []
    k, gn = EXPECTED_COUNTS[g]
    bad = []
    tg = len({r["sample"] for r in crows if r["goi"] == g})
    if n_copies != k: bad.append(f"drawn copies {n_copies}, decided {k}")
    if tg != gn: bad.append(f"per-copy table has {tg} genomes, decided {gn}")
    if up is not None and (up[0], up[1]) != (k, gn):
        bad.append(f"upstream recount {up[0]} copies in {up[1]} genomes, decided {k} in {gn}")
    return bad


def circle_fact_checks(crows):
    """Circle facts belong to the contig: the same under every locus (rows of all loci), and the
    named seamless circles are closed circles, never truncated, and rotations of the reference."""
    bad = []
    cols = ("end_overlap_bp", "end_repeat_bp", "seamless_circle", "closed_circle", "circle_kind")
    if crows and "closed_circle" not in crows[0]:
        return ["per-copy table has no closed_circle column: rebuild with the wrapper"]
    seen = {}
    for r in crows:
        key = (r["sample"], r["contig"])
        val = tuple(r.get(c, "") for c in cols)
        if key in seen and seen[key][0] != val:
            bad.append(f"{key[0]} {key[1]}: circle facts differ between {seen[key][1]} and {r['goi']}")
        seen.setdefault(key, (val, r["goi"]))
    ref = fna(SEAMLESS_REF[0])[SEAMLESS_REF[1]]
    k = next((k for k in range(min(1000, len(ref) // 2), 19, -1) if ref[:k] == ref[-k:]), 0)
    C = ref[:-k] if k else ref
    for s, c in sorted(SEAMLESS):
        q = fna(s).get(c, "")
        if not (len(q) == len(C) and (q in C + C or q.translate(_RC)[::-1] in C + C)):
            bad.append(f"{s} {c}: not a rotation of the reference circle {' '.join(SEAMLESS_REF)}")
        for r in crows:
            if (r["sample"], r["contig"]) == (s, c) and (r.get("truncated") == "True" or r.get("closed_circle") != "True"):
                bad.append(f"{r['goi']} csv v{r['variant']} {s} {c}: seamless circle marked truncated or not closed")
    return bad


def wording_checks(g, cap, vrows):
    """One spelling of the fusion everywhere, and an MGE sentence that names only what is drawn."""
    t = cap.read_text() if cap.exists() else ""
    bad = []
    if re.search(r"ble-MBL|Ble-MBL", t): bad.append("caption spells the fusion's gene 'ble-MBL', the figure 'bleMBL'")
    labs = {r["label"] for r in vrows if r["goi"] == g and r["category"] == "MGE"}
    m = re.search(r"MGE marks ([^.]*)\.", t)
    if m:
        s = m.group(1)
        for word, have in (("tnpM", "tnpM" in labs), ("resolvase", "tnpR" in labs),
                           ("integrase", any(re.match(r"intI|integrase", x) for x in labs)),
                           ("transposase", any(re.match(r"IS|Tn\d+ fam", x) for x in labs))):
            if (word in s) != have:
                bad.append(f"MGE sentence {'names' if word in s else 'misses'} {word}s"
                           f"{'' if have else ' with no referent'}")
    if any(r["goi"] == g and r["label"].rstrip("'") == "IS91/ISCR1" for r in vrows) and "IS91/ISCR1: " not in t:
        bad.append("IS91/ISCR1 drawn but not glossed")
    return bad


# ---- added in the intI1 repair round (2026-10-06) -------------------------------------------
DBSEQ = "/home/bals/Documents/dev/BACTolomew/DB/seq/markers"
_AA = dict(zip((a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"),
               "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
_ref = {}
def ref_prot(fasta, key):
    if (fasta, key) not in _ref:
        seq, take = [], False
        with open(f"{DBSEQ}/{fasta}") as fh:
            for line in fh:
                if line.startswith(">"):
                    if take: break
                    take = key in line
                elif take: seq.append(line.strip())
        _ref[(fasta, key)] = "".join(seq)
    return _ref[(fasta, key)]


def translate(s):
    return "".join(_AA.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def assembly_checks(g, vrows, cap):
    """Arrows that are not Bakta CDS, checked from the sequence on their own (rule a, N gap):
    a ΔbleMBL::ΔqacEΔ1 arrow must translate, on one strand, to Met ... stop with at least 30 aa
    of the BleMBL N-terminus and the last 38 aa of the qacEdelta1 C-terminus, and must abut a
    blaNDM-1 arrow; no Bakta 300 to 330 bp ORF may still sit over the 3' end of a blaNDM-1
    arrow whose fusion was confirmed; every ARG arrow with no Bakta CDS behind it must say why
    (note) and the caption must explain it."""
    bad = []
    rows = [r for r in vrows if r["goi"] == g and r.get("start")]
    if rows and "note" not in rows[0]:
        return ["variants table has no note column: rebuild with the wrapper"]
    ble = ref_prot("AMRF_core-AMR.fasta", "WP_004201167.1")
    qtail = ref_prot("CARD_AMR.fasta", "|qacEdelta1 ")[-38:]
    t = cap.read_text() if cap.exists() else ""
    by_v = collections.defaultdict(list)
    for r in rows: by_v[r["variant"]].append(r)
    for v, rr in by_v.items():
        s = rr[0]["rep_sample"]
        bk = {(a, b) for a, b, _ in cds(s).get(rr[0]["contig"], [])}
        for r in rr:
            a, b = int(float(r["start"])), int(float(r["end"]))
            if r["label"].startswith("ΔbleMBL"):
                q = fna(s)[r["contig"]][a - 1:b]
                ok = any(p.startswith("M") and p.endswith("*") and "*" not in p[:-1]
                         and p[:30] == ble[:30] and p[:-1].endswith(qtail)
                         for p in (translate(q), translate(q.translate(_RC)[::-1])))
                if not ok: bad.append(f"csv v{v}: {r['label']} at {a}-{b} does not translate to the fusion")
                nd = [x for x in rr if x["label"] == "blaNDM-1"]
                if not any(min(abs(a - int(float(x["end"]))), abs(int(float(x["start"])) - b)) <= 12 for x in nd):
                    bad.append(f"csv v{v}: {r['label']} does not abut a blaNDM-1 arrow")
                for x in rr:
                    xa, xb = int(float(x["start"])), int(float(x["end"]))
                    if (x is not r and x["category"] in ("HP", "other", "ARG") and 300 <= xb - xa + 1 <= 330
                            and (xa, xb) in bk and xa <= b and xb >= a and not x["label"].startswith("ΔbleMBL")):
                        bad.append(f"csv v{v}: Bakta ORF {xa}-{xb} still drawn over the confirmed fusion")
            if r["label"] == "blaNDM-1":
                # the read-through ORF after this blaNDM-1, found here independently: one codon
                # after the stop, then Met..stop with the BleMBL start and the qacEdelta1 end
                q = fna(s)[r["contig"]]
                for down in (q[b:b + 600], q[max(0, a - 601):a - 1].translate(_RC)[::-1]):
                    orf = translate(down).split("*")[0]; k = orf.find("M")
                    if 0 <= k <= 2 and orf[k:k + 30] == ble[:30] and orf.endswith(qtail):
                        # rule (a) names the ORF only against the genome's OWN qacEdelta1; a
                        # genome with no qacEdelta1 CDS keeps Bakta's ORF, with the reason noted
                        own = any(re.search(r"QacE delta 1|gene=qacEdelta1", at, re.I)
                                  for cc in cds(s).values() for _a, _b, at in cc)
                        if not own and any((x.get("note") or "").startswith("fusion unconfirmed, no own qacEdelta1")
                                           for x in rr):
                            continue
                        if not any(x["label"].startswith("ΔbleMBL") and min(
                                abs(int(float(x["start"])) - b), abs(a - int(float(x["end"])))) <= 12 for x in rr):
                            bad.append(f"csv v{v}: the bleMBL::qacEdelta1 read-through ORF after blaNDM-1 "
                                       f"at {a}-{b} is not drawn")
            if r["category"] == "ARG" and (a, b) not in bk and not r.get("note"):
                bad.append(f"csv v{v}: ARG arrow {r['label']} at {a}-{b} is no Bakta CDS and carries no note")
    if any(r["label"].startswith("ΔbleMBL") for r in rows) and "read from the assembly on the blaNDM-1" not in t:
        bad.append("caption does not explain the ΔbleMBL::ΔqacEΔ1 arrow")
    for r in rows:
        if r.get("gap_x0") and "is not a Bakta CDS" not in t:
            bad.append(f"caption does not say that {r['label']} is not a Bakta CDS"); break
    return bad


_mob = None
MOBILITY_RE = re.compile(r"tnp|transposase|integrase|resolvase|recombinase|^Int|^IS\w*\d", re.I)
def mobile_checks(g, vrows, cap, lay):
    """Guide 2.6/2.7: a grey arrow with a CGE-MGE or ISFinder MOBILE call is either an MGE
    (a mobility protein matched: fail if grey) or cargo, and then the caption must list it in its
    'Drawn grey although matched there' sentence, on every display track where it is drawn."""
    global _mob
    import pandas as pd
    if _mob is None:
        _mob = pd.read_parquet(MH, columns=["Sample", "Contig", "Start", "Stop", "DB", "SubjectIdentifier"])
        _mob = _mob[_mob["DB"].astype(str).isin(("CGE-MGE-Prot_MOBILE", "ISFinder-Prot_MOBILE"))]
        _mob = _mob.assign(contig="contig_" + _mob["Contig"].astype(str),
                           lo=_mob[["Start", "Stop"]].min(axis=1).astype(int),
                           hi=_mob[["Start", "Stop"]].max(axis=1).astype(int))
    bad, need = [], collections.defaultdict(set)
    disp = {int(r["csv_variant"]): int(r["display"]) for r in lay}
    for r in vrows:
        if r["goi"] != g or not r.get("start") or r["category"] not in ("other", "HP"): continue
        a, b = int(float(r["start"])), int(float(r["end"]))
        h = _mob[(_mob["Sample"] == r["rep_sample"]) & (_mob["contig"] == r["contig"])
                 & ((_mob["lo"] - a).abs() < 60) & ((_mob["hi"] - b).abs() < 60)]
        if not len(h): continue
        prot = [str(x).split("|")[3] if str(x).count("|") >= 3 else str(x).split("//")[-1]
                for x in h["SubjectIdentifier"]]
        if any(MOBILITY_RE.search(p) for p in prot):
            bad.append(f"csv v{r['variant']} {r['label'] or '(unlabelled)'} at {a}: mobility protein "
                       f"({', '.join(sorted(set(prot)))}) drawn {r['category']}")
        need[r["label"] if r["category"] != "HP" and r["label"] else "hypothetical"].add(disp.get(int(r["variant"])))
    if need:
        t = cap.read_text() if cap.exists() else ""
        m = re.search(r"Drawn grey although matched there: ([^\n]*)", t)
        if not m: return bad + ["caption does not list the grey CDS with a MOBILE match"]
        # items read "name (variant 5; Tn502)" or "name (variants 1 and 11; Tn6183; ...)"
        items = {nm.strip(): set(map(int, re.findall(r"\d+", vv)))
                 for nm, vv in re.findall(r"(?:^|\); )([^()]+?) \(variants? ([\d, and]+);", m.group(1))}
        for lab, dv in need.items():
            key = next((k for k in items if k == lab or (lab == "hypothetical" and "hypothetical" in k)), None)
            miss = sorted(x for x in dv if key is None or x not in items[key])
            if miss: bad.append(f"caption's MOBILE list misses {lab} on display {miss}")
    return bad


def dropped_caption_checks(g, cap, up_dropped):
    """Guide 2.2: every contig dropped as a duplicate of a longer one is named in the caption."""
    t = cap.read_text() if cap.exists() else ""
    return [f"caption does not name dropped contig {s} {c}" for s, c in sorted(up_dropped)
            if not re.search(rf"{re.escape(s)} {re.escape(c)}\b", t)]


def fragment_caption_checks(g, cap, lay):
    """A track drawn as a fragment of another (tight pitch, possibly no names of its own) must
    be named as such in the caption, by display numbers."""
    t = cap.read_text() if cap.exists() else ""
    disp = {int(r["csv_variant"]): int(r["display"]) for r in lay}
    bad = []
    for r in lay:
        p = int(float(r.get("fragment_parent_csv") or 0))
        if p and f"Variant {r['display']} (n = {r['n']}) is drawn as a fragment of variant {disp.get(p)}" not in t:
            bad.append(f"caption does not state that display v{r['display']} is a fragment of v{disp.get(p)}")
    return bad


def pdf_words(pdf):
    out = subprocess.run(["pdftotext", "-bbox-layout", str(pdf), "-"],
                         capture_output=True, text=True).stdout
    return out, re.findall(r">([^<]*)</word>", out)


INK_RE = re.compile(r"stroke:rgb\(7\.4\d+%,12\.5\d+%,19\.9\d+%\)")   # #132033, the rule
SEG_RE = re.compile(r'<path\s+style="([^"]*)"\s+d="([^"]+)"')
LINE_RE = re.compile(r"M\s*([-\d.]+)\s+([-\d.]+)\s*L\s*([-\d.]+)\s+([-\d.]+)\s*Z?")
PT2MM = 25.4 / 72


def bar_mm(pdf, tmpdir):
    """Scale bar length in mm, read from the PDF's vector geometry.

    Measured off the raster instead, this is unreliable: the dark run includes the end ticks'
    stroke width and an antialiased fringe, which read about 1.3 percent long at 300 dpi and,
    because the fringe does not scale with resolution, about 1 percent SHORT at 1200. The
    same figure measured 4 percent apart between the two. The vector path is exact and
    resolution independent, so it is the only sound way to check a drawn scale.

    The bar is picked out by stroke colour: it is the only horizontal rule drawn in the ink
    colour, while the track baselines are grey and the arrow outlines are slate.
    """
    svg = pathlib.Path(tmpdir) / (pathlib.Path(pdf).stem + ".svg")
    subprocess.run(["pdftocairo", "-svg", str(pdf), str(svg)],
                   check=True, capture_output=True)
    body = svg.read_text()
    best = 0.0
    for style, d in SEG_RE.findall(body):
        if not INK_RE.search(style):
            continue
        m = LINE_RE.fullmatch(d.strip())
        if not m:
            continue
        x0, y0, x1, y1 = map(float, m.groups())
        if abs(y1 - y0) < 0.01:
            best = max(best, abs(x1 - x0))
    return best * PT2MM if best else None


LABEL_RE = re.compile(r"n = (\d+) cop(y|ies)(?: in (\d+) genomes?)?(?:\s+\(([^()]*)\))?"
                      # the cluster may carry a "(k/n)" share, and pdftotext can join the
                      # bold "variant N" to its molecule letter ("variant 1P")
                      r"(?:\s+·\s+(\S+(?: \(\d+/\d+\))?))?\s+variant\s*(\d+)")


def drawn_labels(pdf):
    """Every track label as drawn:
    {display: (copies, genomes, {country: k}, rest, text, has_countries)}.

    Read from the PDF's laid-out text, where each label sits on one line followed by its bold
    "variant N" column. Duplicate display numbers are reported, not merged."""
    out = subprocess.run(["pdftotext", "-layout", str(pdf), "-"],
                         capture_output=True, text=True).stdout
    labels, dup = {}, []
    for m in LABEL_RE.finditer(out):
        k, plural, gen, ctry, cl, disp = m.groups()
        k, disp = int(k), int(disp)
        printed, rest = {}, 0
        same = re.findall(r"same genome as variant (\d+)", ctry or "")
        for part in [p.strip() for p in (ctry or "").split(";")[0].split(",") if p.strip()]:
            mo = re.fullmatch(r"\+(\d+) other", part)
            if mo: rest = int(mo.group(1)); continue
            name, n = part.rsplit(" ", 1); printed[name] = int(n)
        text = (f"n = {k} cop{plural}" + (f" in {gen} genome" + ("s" if gen != "1" else "")
                                          if gen else "") + (f"  ({ctry})" if ctry else "")
                + (f"  ·  {cl}" if cl else ""))
        if disp in labels: dup.append(disp)
        labels[disp] = (k, int(gen) if gen else k, printed, rest, text, ctry is not None,
                        [int(x) for x in same])
    return labels, dup


def label_checks(g, labels, dup, exp, lay):
    """Counts the subtitle used to state, now read off the drawn track labels."""
    bad = []
    if dup: bad.append(f"display numbers drawn twice: {sorted(set(dup))}")
    if sorted(labels) != list(range(1, exp["tracks"] + 1)):
        bad.append(f"drawn track labels for displays {sorted(labels)}, expected 1..{exp['tracks']}")
    lay_of = {int(r["display"]): r for r in lay}
    for disp, (k, gen, printed, rest, text, has_ctry, _same) in sorted(labels.items()):
        if has_ctry and sum(printed.values()) + rest != gen:
            bad.append(f"display v{disp}: drawn countries sum to {sum(printed.values()) + rest}, "
                       f"label says {gen} genomes")
        r = lay_of.get(disp)
        if lay and r is None:
            bad.append(f"display v{disp}: not in the layout CSV"); continue
        if r is not None:
            if r["label"] != text:
                bad.append(f"display v{disp}: drawn '{text}' but layout says '{r['label']}'")
            csv_t = exp["csv_track"].get(int(r["csv_variant"]), "")
            # The count part of the csv track text ("n = K copies  (countries)") must equal the
            # drawn one exactly. The table writes a same-genome mark with the CSV number ("same
            # genome as variant 22"); it is compared through the layout map, as drawn.
            disp_of = {int(x["csv_variant"]): int(x["display"]) for x in lay}
            want = csv_t.split("  ·  ")[-1].strip()
            if "same genome as csv variant" in want:
                bad.append(f"csv v{r['csv_variant']}: track text says 'csv variant', the figure 'variant'")
            want = re.sub(r"same genome as (?:csv )?variant ([\d, ]+\d)",
                          lambda m: "same genome as variant " + ", ".join(
                              str(disp_of.get(int(x), f"csv{x}")) for x in m.group(1).split(", ")), want)
            got = text.split("  ·  ")[0].strip()
            if want != got:
                bad.append(f"display v{disp}: drawn '{got}' but csv v{r['csv_variant']} track text "
                           f"reads '{want}'")
    n_k = sum(v[0] for v in labels.values()); n_g = sum(v[1] for v in labels.values())
    if n_k != exp["copies"]:
        bad.append(f"drawn copies sum to {n_k}, CSV n sums to {exp['copies']}")
    if exp["track_genomes"] is not None and n_g != exp["track_genomes"]:
        bad.append(f"drawn genomes sum to {n_g}, per-copy table gives {exp['track_genomes']}")
    # A country that no drawn "+K other" can hide must add up, over the drawn labels, to the
    # same tally made per track from the metadata. Compared per track, not against distinct
    # genomes: a genome whose copies sit on two tracks is correctly named on both (F20).
    if exp["track_truth"] and all(v[5] for v in labels.values()):
        hiding = any(v[3] for v in labels.values())
        drawn = collections.Counter()
        for v in labels.values(): drawn.update(v[2])
        for c, t in exp["track_truth"].items():
            if drawn.get(c, 0) > t or (not hiding and drawn.get(c, 0) != t):
                bad.append(f"'{c}' drawn {drawn.get(c, 0)} times over the tracks, metadata "
                           f"gives {t} (per track)")
    return n_k, n_g, bad


def chrome_checks(pdf, raw_text, words, bar_y, boxes, tmpdir):
    """No subtitle, no notes, no pink second-copy fill or key (unless switched back on)."""
    bad = []
    if not SHOW_SUBTITLE and re.search(r"\bcontexts\b.*\ball shown\b", raw_text):
        bad.append("subtitle drawn")
    if not SHOW_NOTES:
        below = [w for y, w in boxes if y > bar_y + 1 and w.strip()]
        if below: bad.append(f"text below the scale bar (notes block?): {' '.join(below[:8])}")
    if "additional copy" in " ".join(words):
        bad.append("'additional copy of gene of interest' key drawn")
    svg = pathlib.Path(tmpdir) / (pathlib.Path(pdf).stem + ".svg")
    if svg.exists() and PINK_RE.search(svg.read_text()):
        bad.append("pink #f08c8c fill drawn (a second focal copy must use the focal red)")
    return bad


def style_metrics(raw, svg_path, page):
    """Sizes and margins that must be the same on every figure of the set, in pt: title and
    legend-header text height, scale-bar end-tick height, bar to the top of its label, bar to the
    bottom and right page edges, top of the title to the page top. Compared across the loci
    checked in one run (set_style_checks)."""
    words = [(float(a), float(b), float(c), float(e), w) for a, b, c, e, w in WORD_RE.findall(raw)]
    def h(word):
        x = next((w for w in words if w[4] == word), None)
        return round(x[3] - x[1], 2) if x else None
    out = {"title height": h("baumannii"), "legend header height": h("Gene"),
           "legend item height": h("chromosome")}
    t = next((w for w in words if w[4] == "baumannii"), None)
    if t: out["title to page top"] = round(t[1], 2)
    body = pathlib.Path(svg_path).read_text()
    bar, ticks = None, []
    for style, d in SEG_RE.findall(body):
        if not INK_RE.search(style): continue
        m = LINE_RE.fullmatch(d.strip())
        if not m: continue
        x0, y0, x1, y1 = map(float, m.groups())
        if abs(y1 - y0) < 0.01 and (bar is None or abs(x1 - x0) > abs(bar[2] - bar[0])): bar = (x0, y0, x1, y1)
        elif abs(x1 - x0) < 0.01: ticks.append((x0, min(y0, y1), max(y0, y1)))
    if bar:
        tk = [t2 - t1 for x, t1, t2 in ticks if min(abs(x - bar[0]), abs(x - bar[2])) < 0.5]
        if tk: out["bar tick height"] = round(max(tk), 2)
        kb = [w for w in words if w[4] in ("kb", "bp") and w[1] > bar[1]]
        if kb: out["bar to label top"] = round(min(w[1] for w in kb) - bar[1], 2)
        out["bar to page bottom"] = round(page[1] - bar[1], 2)
        out["bar to page right"] = round(page[0] - max(bar[0], bar[2]), 2)
    return out


def set_style_checks(styles, tol=0.25):
    bad = []
    keys = sorted({k for st in styles.values() for k in st})
    for k in keys:
        vals = {g: st.get(k) for g, st in styles.items() if st.get(k) is not None}
        # the bottom and right margins also hold the fraction of a point that makes each page
        # a whole number of points
        lim = tol + 1.0 if k in ("bar to page bottom", "bar to page right") else tol
        if len(vals) > 1 and max(vals.values()) - min(vals.values()) > lim:
            bad.append(f"{k} differs across the set: " + ", ".join(f"{g} {v:g}" for g, v in vals.items()))
    return bad


def page_pt(pdf):
    out = subprocess.run(["pdfinfo", str(pdf)], capture_output=True, text=True).stdout
    m = re.search(r"Page size:\s+([\d.]+) x ([\d.]+) pts", out)
    return (float(m.group(1)), float(m.group(2))) if m else (None, None)


# ---- drawn geometry from the vector PDF (pdftocairo -svg) -----------------------------------
# Every glyph is a <use> of an outline symbol at page coordinates (rotated names come as
# pre-rotated outlines), every gene arrow a filled polygon with the slate outline #5a6572 and
# every legend key a small filled square with no stroke, so overprints and fills can be
# checked exactly, without the raster.
PAL_HEX = {"#c1121f": "center", "#7048e8": "ARG", "#ff8f00": "MGE", "#0c8599": "TA",
       "#1971c2": "porin", "#adb5bd": "HP", "#d7dde3": "other"}
KEY_TEXT = {"center": "gene of interest", "ARG": "other ARG", "MGE": "MGE (tnp", "TA": "toxin-antitoxin",
            "porin": "porin", "HP": "hypothetical", "other": "other CDS"}
NUM_RE = re.compile(r"-?\d+(?:\.\d+)?(?:e-?\d+)?")


def rgb_hex(txt):
    p = [float(x) for x in re.findall(r"([\d.]+)%", txt)][:3]
    return "#" + "".join(f"{round(v * 2.55):02x}" for v in p)


def pairs(d):
    n = [float(x) for x in NUM_RE.findall(d)]
    return list(zip(n[0::2], n[1::2]))


def hull(pts):
    pts = sorted(set(pts))
    if len(pts) < 3: return pts
    def cross(o, a, b): return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lo, up = [], []
    for p in pts:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], p) <= 0: lo.pop()
        lo.append(p)
    for p in reversed(pts):
        while len(up) >= 2 and cross(up[-2], up[-1], p) <= 0: up.pop()
        up.append(p)
    return lo[:-1] + up[:-1]


def shrink(poly, by):
    """Move every vertex toward the centroid by `by` pt: touching outlines are not overprints."""
    if not poly: return poly
    cx = sum(p[0] for p in poly) / len(poly); cy = sum(p[1] for p in poly) / len(poly)
    out = []
    for x, y in poly:
        dx, dy = x - cx, y - cy; r = (dx * dx + dy * dy) ** 0.5
        k = max(0.0, (r - by) / r) if r else 0.0
        out.append((cx + dx * k, cy + dy * k))
    return out


def convex_overlap(a, b):
    """Separating-axis test for two convex polygons."""
    if len(a) < 3 or len(b) < 3: return False
    for poly in (a, b):
        for i in range(len(poly)):
            x1, y1 = poly[i]; x2, y2 = poly[(i + 1) % len(poly)]
            nx, ny = y2 - y1, x1 - x2
            pa = [nx * x + ny * y for x, y in a]; pb = [nx * x + ny * y for x, y in b]
            if max(pa) <= min(pb) or max(pb) <= min(pa): return False
    return True


def inside(pt, poly):
    x, y = pt; c = False
    for i in range(len(poly)):
        x1, y1 = poly[i]; x2, y2 = poly[i - 1]
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1: c = not c
    return c


def svg_geometry(svg_path):
    body = pathlib.Path(svg_path).read_text()
    sym = {}
    for sid, block in re.findall(r'<symbol[^>]*id="([^"]+)"[^>]*>(.*?)</symbol>', body, re.S):
        sym[sid] = [p for d in re.findall(r'd="([^"]*)"', block) for p in pairs(d)]
    defs_end = body.find("</defs>")
    main = body[defs_end:]
    glyphs = []                            # (index, colour, convex hull) in drawing order
    for style, uses in re.findall(r'<g style="fill:(rgb\([^)]*\))[^"]*">\s*((?:<use[^>]*/>\s*)+)</g>', main):
        for sid, x, y in re.findall(r'<use xlink:href="#([^"]+)" x="([-\d.]+)" y="([-\d.]+)"', uses):
            pts = [(px + float(x), py + float(y)) for px, py in sym.get(sid, [])]
            if len(pts) >= 3: glyphs.append((len(glyphs), rgb_hex(style), hull(pts)))
    arrows, keys, dashes = [], [], []
    for style, d in re.findall(r'<path\s+style="([^"]*)"\s+d="([^"]*)"', main):
        fm = re.search(r"fill:(rgb\([^)]*\))", style)
        if "stroke:rgb(35.29" in style and fm:
            arrows.append((rgb_hex(fm.group(1)), pairs(d)))
        elif fm and "stroke:none" in style:
            pts = pairs(d); xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
            w, h = max(xs) - min(xs), max(ys) - min(ys)
            col = rgb_hex(fm.group(1))
            if 3 < w < 10 and abs(w - h) < 0.3 and col != "#ffffff":
                keys.append((col, (min(xs), min(ys), max(xs), max(ys))))
        elif "stroke-dasharray" in style and "stroke:rgb(75.68" in style:
            dashes.append(pairs(d))
    return glyphs, arrows, keys, dashes


# Name against name, on the glyph OUTLINES (not their hulls, which shrunk by 0.25 pt let two
# names of one band touch: on blaOXA-72 display v11 the tail of the Q of "tolQ" met the "b" of
# "exbD"). Glyphs are grouped into the strings they belong to (each next glyph one advance away
# along the text direction, 0 or 40 degrees); outlines of glyphs from two different strings must
# stay GLYPH_CLEAR_PT apart.
GLYPH_CLEAR_PT = 0.3


def _subpaths(d, step=6):
    """Closed outline polylines of one glyph symbol (M starts a new one)."""
    out, cur = [], []
    toks = re.findall(r"[MLCZ]|-?[\d.]+(?:e-?\d+)?", d)
    i, cmd, pos = 0, None, (0.0, 0.0)
    while i < len(toks):
        t = toks[i]
        if t in "MLCZ":
            cmd = t; i += 1
            if t == "M" and cur: out.append(cur); cur = []
            continue
        if cmd in ("M", "L"):
            pos = (float(toks[i]), float(toks[i + 1])); cur.append(pos); i += 2
        elif cmd == "C":
            c = [float(x) for x in toks[i:i + 6]]; p0 = pos
            for k in range(1, step + 1):
                u = k / step; a, b, cc, e = (1 - u) ** 3, 3 * (1 - u) ** 2 * u, 3 * (1 - u) * u ** 2, u ** 3
                cur.append((a * p0[0] + b * c[0] + cc * c[2] + e * c[4], a * p0[1] + b * c[1] + cc * c[3] + e * c[5]))
            pos = (c[4], c[5]); i += 6
        else: i += 1
    if cur: out.append(cur)
    return out


def name_contacts(svg_path):
    """Pairs of glyphs from different strings whose outlines come within GLYPH_CLEAR_PT."""
    body = pathlib.Path(svg_path).read_text()
    sym = {}
    for sid, block in re.findall(r'<symbol[^>]*id="([^"]+)"[^>]*>(.*?)</symbol>', body, re.S):
        sym[sid] = [sp for d in re.findall(r'd="([^"]*)"', block) for sp in _subpaths(d)]
    main = body[body.find("</defs>"):]
    gl = []                                  # (string id, polylines, bbox)
    sid_n, prev, rot = -1, None, set()        # rot: strings set at 40 degrees (the gene names)
    for sid, x, y in re.findall(r'<use xlink:href="#([^"]+)" x="([-\d.]+)" y="([-\d.]+)"', main):
        x, y = float(x), float(y)
        polys = [[(px + x, py + y) for px, py in sp] for sp in sym.get(sid, []) if len(sp) >= 2]
        same = False
        if prev is not None:
            dx, dy = x - prev[0], y - prev[1]; r = (dx * dx + dy * dy) ** 0.5
            ang = math.degrees(math.atan2(-dy, dx)) if r else 0.0
            same = 0 < r <= 6.5 and (abs(ang) < 15 or abs(ang - 40) < 15)
            if same and abs(ang - 40) < 15: rot.add(sid_n + (0 if same else 1))
        if not same: sid_n += 1
        prev = (x, y)
        if not polys: continue
        xs = [q[0] for sp in polys for q in sp]; ys = [q[1] for sp in polys for q in sp]
        gl.append((sid_n, polys, (min(xs), min(ys), max(xs), max(ys))))
    def dist(a, b):
        best = 1e9
        for sa in a:
            for sb in b:
                for pt in sa:
                    for q0, q1 in zip(sb, sb[1:] + sb[:1]):
                        best = min(best, _seg_dist(pt, q0, q1))
                        if best < 1e-3: return 0.0
        return best
    def inside_any(a, b):
        return any(_inside(pt, sb) for sa in a for pt in sa[:1] for sb in b if len(sb) >= 3)
    order = sorted(range(len(gl)), key=lambda k: gl[k][2][0])
    hits = []
    c = GLYPH_CLEAR_PT
    for ai, a in enumerate(order):
        A = gl[a]
        for b in order[ai + 1:]:
            B = gl[b]
            if B[2][0] > A[2][2] + c: break
            if A[0] == B[0] or B[2][1] > A[2][3] + c or B[2][3] < A[2][1] - c: continue
            # gene names only: the title is set as three abutting runs ("A. baumannii", ...)
            if A[0] not in rot and B[0] not in rot: continue
            if dist(A[1], B[1]) < c or inside_any(A[1], B[1]) or inside_any(B[1], A[1]):
                hits.append(((A[2][0] + A[2][2]) / 2, (A[2][1] + A[2][3]) / 2))
    return hits


def geometry_checks(geo, svg_path=None):
    """No two glyphs overlapping, the focal line through no glyph, no arrow drawn over a gene of
    interest. (Glyph ink against arrows is overprint_checks, with OVERPRINT_PT to spare.)"""
    glyphs, arrows, _keys, dashes = geo
    bad = []
    small = [(i, c, shrink(h, 0.25)) for i, c, h in glyphs]
    boxes = [(min(p[0] for p in h), min(p[1] for p in h), max(p[0] for p in h), max(p[1] for p in h))
             for _i, _c, h in small]
    order = sorted(range(len(small)), key=lambda k: boxes[k][0])
    clash = 0
    for a_i, a in enumerate(order):
        for b in order[a_i + 1:]:
            if boxes[b][0] > boxes[a][2]: break
            if boxes[b][1] > boxes[a][3] or boxes[b][3] < boxes[a][1]: continue
            # two glyphs drawn one after the other are neighbours in one string: a kerned pair
            # ("Ly" of "LysR fam.") overlaps as convex hulls although the ink does not touch.
            # Contacts between two different names are tested on the outlines (name_contacts).
            if abs(small[a][0] - small[b][0]) == 1: continue
            if convex_overlap(small[a][2], small[b][2]): clash += 1
    if clash: bad.append(f"{clash} pair(s) of glyphs overlap (labels overprinting each other)")
    if svg_path is not None:
        nc = name_contacts(svg_path)
        if nc: bad.append(f"{len(nc)} glyph pair(s) of two different names within {GLYPH_CLEAR_PT} pt "
                          f"(at {', '.join(f'{x:.0f},{y:.0f}' for x, y in nc[:4])} pt)")
    cross = 0
    for seg in dashes:
        for (x1, y1), (x2, y2) in zip(seg, seg[1:]):
            n = max(2, int(abs(y2 - y1) / 0.3) + 1)
            pts = [(x1 + (x2 - x1) * t / n, y1 + (y2 - y1) * t / n) for t in range(n + 1)]
            for (i, _c, h), (bx0, by0, bx1, by1) in zip(small, boxes):
                if bx1 < min(x1, x2) or bx0 > max(x1, x2): continue
                if any(by0 <= py <= by1 and inside((px, py), h) for px, py in pts): cross += 1
    if cross: bad.append(f"dashed focal line runs through {cross} glyph(s)")
    # every gene-of-interest arrow is drawn after any arrow that overlaps it, so none hides it
    red = [k for k, (c, _p) in enumerate(arrows) if c == "#c1121f"]
    hidden = sum(1 for k in red for j in range(k + 1, len(arrows))
                 if arrows[j][0] != "#c1121f" and convex_overlap(shrink(hull(arrows[k][1]), 0.25),
                                                                 shrink(hull(arrows[j][1]), 0.25)))
    if hidden: bad.append(f"{hidden} arrow(s) drawn over a gene-of-interest arrow")
    return bad


def key_checks(geo, words):
    """Arrow fills are palette colours only, the gene-category keys sit on one row, every legend
    key square has an arrow in its colour and every drawn colour a key. (The fill counts per
    category are fill_checks.)"""
    _glyphs, arrows, keys, _d = geo
    bad = []
    drawn = collections.Counter(c for c, _p in arrows)
    off = {c: k for c, k in drawn.items() if c not in PAL_HEX}
    if off: bad.append(f"arrow fills outside the palette: {off}")
    key_cols = collections.Counter(c for c, _b in keys)
    # guide section 3: one row of gene-category keys under its header
    rows_y = sorted({round(b[1]) for _c, b in keys})
    if len(rows_y) > 1: bad.append(f"gene-category keys on {len(rows_y)} rows, the guide asks for one")
    for c in key_cols:
        if c in PAL_HEX and drawn.get(c, 0) == 0:
            bad.append(f"legend key '{KEY_TEXT[PAL_HEX[c]]}' has no arrow drawn in its colour")
    for c in drawn:
        if c in PAL_HEX and c not in key_cols:
            bad.append(f"arrows drawn in {c} ({PAL_HEX[c]}) but no legend key for it")
    text = " ".join(words)
    for cat, t in KEY_TEXT.items():
        col = next(c for c, k in PAL_HEX.items() if k == cat)
        if t in text and drawn.get(col, 0) == 0:
            bad.append(f"legend text '{t}' printed with nothing drawn in {col}")
    return bad


def dash_checks(raw_text, cap_text):
    bad = []
    if re.search("[\u2013\u2014]", raw_text): bad.append("en or em dash in the figure text")
    if cap_text and re.search("[\u2013\u2014]", cap_text): bad.append("en or em dash in the caption")
    return bad


# ---- track text against the per-copy table and MOB-suite -----------------------------------
MOB = "/home/bals/Documents/dev/BACTolomew/runs/sample_78eaa9126850/run_78eaa9126850/mobsuite/per_sample"
_mobsuite = {}
def mob_cluster(sample, contig):
    if sample not in _mobsuite:
        rows = read_rows_tsv(f"{MOB}/{sample}/contig_report.txt")
        _mobsuite[sample] = {r["contig_id"]: (r["molecule_type"], r["primary_cluster_id"]) for r in rows}
    return _mobsuite[sample].get(contig, (None, None))


def read_rows_tsv(path):
    if not os.path.exists(path): return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def location_text(copies):
    """The wrapper's location text for a track (strict-majority rule, _st492_location),
    recomputed from its copies."""
    letter, cl = expected_location(copies)
    if letter == "U": return {"location unresolved"}
    if letter == "C": return {"chromosome"}
    n_ch = sum(1 for c in copies if c["location"] == "chromosome")
    return {f"plasmid {cl}" + (f", {n_ch} chromosomal" if n_ch else "")}


def cluster_checks(g, vrows, crows, labels, lay):
    bad = []
    by_var = collections.defaultdict(list)
    for r in crows:
        if r["goi"] == g: by_var[int(r["variant"])].append(r)
    track = {int(r["variant"]): r["track"] for r in vrows if r["goi"] == g}
    for v, t in track.items():
        parts = [p.strip() for p in t.split("·")]
        if len(parts) < 3: bad.append(f"csv v{v}: track text has no location part"); continue
        want = location_text(by_var.get(v, []))
        if parts[1] not in want:
            bad.append(f"csv v{v}: location '{parts[1]}' but its copies give {sorted(want)}")
        for c in by_var.get(v, []):
            mt, cl = mob_cluster(c["sample"], c["contig"].replace("contig_", ""))
            if mt is None: mt, cl = mob_cluster(c["sample"], c["contig"])
            if mt is not None and (mt != c["mob_molecule_type"] or str(cl) != str(c["mob_primary_cluster_id"])):
                bad.append(f"csv v{v} {c['sample']} {c['contig']}: per-copy MOB call "
                           f"{c['mob_molecule_type']}/{c['mob_primary_cluster_id']}, contig_report says {mt}/{cl}")
    lay_of = {int(r["display"]): int(r["csv_variant"]) for r in lay}
    for disp, lab in labels.items():
        v = lay_of.get(disp)
        if v is None: continue
        loc = [p.strip() for p in track[v].split("·")][1]
        shown = loc.replace("location unresolved", "").replace("chromosome", "").replace("plasmid ", "").strip()
        drawn = lab[4].split("·")[-1].strip() if "·" in lab[4] else ""
        if shown != drawn:
            bad.append(f"display v{disp}: drawn cluster '{drawn}' but csv v{v} says '{loc}'")
    return bad


# ---- arrows against Bakta and marker_hits (--upstream) -------------------------------------
AMR_DBS = ("AMRF_core-AMR", "AMRF_plus-AMR", "CARD_AMR")
FROM_ASSEMBLY = ("fusion", "contig end, read from the assembly", "scaffold gap", "ambiguity codes")
MOBILE_DBS = ("CGE-MGE-Prot_MOBILE", "ISFinder-Prot_MOBILE")


def arrow_upstream_checks(g, vrows, crows):
    """Each drawn arrow is a Bakta CDS of the representative genome (or is marked as read from
    the assembly); a CDS with an AMR call at >= 80 percent coverage is never grey unless kept
    grey on purpose; aph(3') carries its subclass; a family-only IS transposase is not printed
    as a specific element; the fusion ORF lies on the blaNDM-1 strand after its stop codon."""
    import pandas as pd
    mh = pd.read_parquet(MH, columns=["Sample", "Contig", "Start", "Stop", "Gene", "DB",
                                      "Subject_coverage", "ResistanceGeneFamilies",
                                      "ResistanceAlleles", "SubjectIdentifier"])
    allmob = mh[mh["DB"].astype(str).isin(MOBILE_DBS)]
    mob = allmob[allmob["Subject_coverage"].astype(float) >= 80]
    # resistance mutations in housekeeping genes (parC_S84F, CARD "VM|" models) are not
    # acquired ARGs, as in the wrapper and marker_checks
    mh = mh[mh["DB"].astype(str).isin(AMR_DBS) & (mh["Subject_coverage"].astype(float) >= 80)
            & ~mh["ResistanceAlleles"].astype(str).str.fullmatch(r"\w+_[A-Z]\d+[A-Z*]")
            & ~mh["SubjectIdentifier"].astype(str).str.startswith("VM|")]
    bad = []
    reps = {}
    for r in vrows:
        if r["goi"] == g: reps.setdefault(int(r["variant"]), []).append(r)
    for v, rows in reps.items():
        rs = rows[0]["rep_sample"]
        cp = next((c for c in crows if c["goi"] == g and int(c["variant"]) == v and c["sample"] == rs), None)
        if cp is None: bad.append(f"csv v{v}: representative {rs} not in the per-copy table"); continue
        contig, cstart = cp["contig"], int(cp["cstart"])
        allc = cds(rs).get(contig, [])
        cen = next(((a, b, at) for a, b, at in allc if a == cstart), None)
        if cen is None: bad.append(f"csv v{v}: focal CDS not found in {rs}"); continue
        strand = {}
        with open(f"{BK}/{rs}/{rs}.gff3") as fh:
            for line in fh:
                f = line.split("\t")
                if len(f) > 8 and f[2] == "CDS" and f[0] == contig: strand[(int(f[3]), int(f[4]))] = f[6]
        cm = (cen[0] + cen[1]) / 2; flip = strand[(cen[0], cen[1])] == "-"
        hs = mh[(mh["Sample"] == rs) & (mh["Contig"].astype(str) == contig.replace("contig_", ""))]
        cen_row = next(r for r in rows if r["category"] == "center")
        # Every Bakta CDS inside the drawn span is drawn (or replaced by an arrow read from the
        # assembly over at least half of it). intI1 csv v14 drew msr(E)' beyond an undrawn ISEc29,
        # an empty backbone that read as "no gene here".
        span = []
        for r in rows:
            x0, x1 = float(r["x0"]), float(r["x1"])
            span.append((round(x0 + cm), round(x1 + cm)) if not flip else (round(cm - x1), round(cm - x0)))
        lo_s, hi_s = min(a for a, _ in span), max(b for _, b in span)
        for a, b, _at in allc:
            if a < lo_s or b > hi_s: continue
            cov = sum(max(0, min(b, d1) - max(a, d0) + 1) for d0, d1 in span)
            if cov < 0.5 * (b - a + 1):
                bad.append(f"csv v{v}: Bakta CDS {a}-{b} lies inside the drawn span of {rs} {contig} but is not drawn")
        for r in rows:
            x0, x1 = float(r["x0"]), float(r["x1"])
            a, b = (round(x0 + cm), round(x1 + cm)) if not flip else (round(cm - x1), round(cm - x0))
            note = r.get("note", "") or ""
            if note.startswith(FROM_ASSEMBLY):
                if r["label"].startswith("ΔbleMBL"):
                    # against the blaNDM-1 arrow it follows (the focal gene on blaNDM-1 tracks,
                    # any blaNDM-1 arrow on another locus's track), in drawn orientation
                    nd = [o for o in rows if o["label"] == "blaNDM-1"] or [cen_row]
                    o = min(nd, key=lambda o: abs(float(o["x0"]) + float(o["x1"]) - x0 - x1))
                    fwd = o["orient"] in ("1", "1.0", "True")
                    if r["orient"] != o["orient"] or (x0 <= float(o["x1"]) if fwd else x1 >= float(o["x0"])):
                        bad.append(f"csv v{v}: fusion ORF not on the blaNDM-1 strand after its stop")
                # an arrow redrawn over ambiguity codes must span some, and no Bakta CDS may
                # still be drawn mostly inside it
                if note.startswith("ambiguity codes"):
                    if not any(ch not in "ACGT" for ch in fna(rs)[contig][a - 1:b]):
                        bad.append(f"csv v{v}: '{r['label']}' drawn for ambiguity codes, none in {a}..{b}")
                    for o in rows:
                        if o is r or (o.get("note") or "").startswith(FROM_ASSEMBLY): continue
                        ox0, ox1 = float(o["x0"]), float(o["x1"])
                        if min(ox1, x1) - max(ox0, x0) > 0.5 * (ox1 - ox0):
                            bad.append(f"csv v{v}: '{o['label']}' still drawn inside '{r['label']}'")
                continue
            hit = next(((aa, bb, at) for aa, bb, at in allc if abs(aa - a) <= 1 and abs(bb - b) <= 1), None)
            if hit is None:
                bad.append(f"csv v{v}: arrow '{r['label']}' at {x0:g}..{x1:g} is no Bakta CDS of {rs}"); continue
            prod = re.search(r"product=([^;\n]*)", hit[2]); prod = prod.group(1) if prod else ""
            lo = hs[["Start", "Stop"]].min(axis=1); hi = hs[["Start", "Stop"]].max(axis=1)
            ov = (hi.clip(upper=b) - lo.clip(lower=a) + 1).clip(lower=0)
            amr = hs[(ov >= 0.9 * (b - a + 1)) & (ov >= 0.9 * (hi - lo + 1))]
            genes = set(amr["Gene"].dropna().astype(str)) | set(amr["ResistanceGeneFamilies"].dropna().astype(str))
            if len(amr) and r["category"] not in ("ARG", "center", "goi_copy") and not genes & KEEP_GREY:
                bad.append(f"csv v{v}: '{r['label'] or prod}' has an AMR call {sorted(genes)} but is drawn {r['category']}")
            if r["label"] == "aph(3')" and any(x.startswith("aph(3')-") for x in genes):
                bad.append(f"csv v{v}: aph(3') printed without the subclass its call gives {sorted(genes)}")
            # an element name needs the evidence the wrapper recorded (ISFinder or the partner CDS)
            if re.fullmatch(r"IS\w*?\d+\w* family transposase", prod) and re.fullmatch(r"IS\w*?\d+", r["label"]) \
                    and not note.startswith("IS element "):
                bad.append(f"csv v{v}: '{r['label']}' names an element, Bakta gives only '{prod}'")
            if g == "blaNDM-1" and r["label"].startswith("ΔbleMBL"):
                bad.append(f"csv v{v}: fusion name on a Bakta CDS ({prod}), not on the ORF read from the assembly")
            # a mobile-element call drawn grey must carry the note the caption explains it by
            ms = mob[(mob["Sample"] == rs) & (mob["Contig"].astype(str) == contig.replace("contig_", ""))]
            lo2 = ms[["Start", "Stop"]].min(axis=1); hi2 = ms[["Start", "Stop"]].max(axis=1)
            ov2 = (hi2.clip(upper=b) - lo2.clip(lower=a) + 1).clip(lower=0)
            if r["category"] in ("other", "HP") and ((ov2 >= 0.9 * (b - a + 1)) & (ov2 >= 0.9 * (hi2 - lo2 + 1))).any() \
                    and "MOBILE call" not in note:
                bad.append(f"csv v{v}: '{r['label'] or prod}' has a mobile-element call but is drawn "
                           f"{r['category']} with no reason recorded")
            # "no ISFinder match" may only be said of a transposase that has none
            if r["label"].endswith(" fam.") and "no ISFinder match" in note:
                am = allmob[(allmob["Sample"] == rs) & (allmob["Contig"].astype(str) == contig.replace("contig_", ""))
                            & (allmob["DB"].astype(str) == "ISFinder-Prot_MOBILE")]
                ov3 = (am[["Start", "Stop"]].max(axis=1).clip(upper=b) - am[["Start", "Stop"]].min(axis=1).clip(lower=a) + 1).clip(lower=0)
                isf = am[ov3 >= 0.9 * (b - a + 1)]
                if len(isf):
                    bad.append(f"csv v{v}: '{r['label']}' noted as without ISFinder match, but ISFinder calls it")
    return bad


# Gene labels that are abbreviations, which the caption must expand when printed
ABBREV = ("DNA BP", "DNA rep", "CHAP", "TBDR", "EAL", "GGDEF", "RRM", "HEPN", "CSLREA")


def rotated_names(raw, ctol=1.8, gap=2.5):
    """The 40 degree gene names as strings, rebuilt from pdftotext's word boxes (it splits a
    diagonal name into pieces of one to three letters): pieces on one text line, across which
    their box centres agree within ctol pt, and close along it are joined. Spaces are lost
    ("DNA rep" reads "DNArep"). Heuristic, used only to ask whether a name is printed at all."""
    cs, sn = math.cos(math.radians(40)), math.sin(math.radians(40))
    fr = []
    for x0, y0, x1, y1, w in re.findall(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" '
                                        r'yMax="([\d.]+)">([^<]*)</word>', raw):
        cx, cy = (float(x0) + float(x1)) / 2, (float(y0) + float(y1)) / 2
        fr.append((cx * cs - cy * sn, cx * sn + cy * cs, w.replace("&apos;", "'").replace("&amp;", "&")))
    lines = []
    for f in sorted(fr, key=lambda f: f[1]):
        for ln in lines:
            if abs(ln[0] - f[1]) <= ctol: ln[1].append(f); break
        else: lines.append([f[1], [f]])
    out = set()
    for _c, fs in lines:
        fs.sort(key=lambda f: f[0]); cur = [fs[0]]
        for f in fs[1:]:
            if f[0] - cur[-1][0] <= (len(cur[-1][2]) + len(f[2])) * 2.2 + gap: cur.append(f)
            else: out.add("".join(x[2] for x in cur)); cur = [f]
        out.add("".join(x[2] for x in cur))
    return out


def glossary_checks(g, cap, vrows, words=None):
    """The caption explains every primed, family-only or read-from-the-assembly name the table
    draws, and every printed abbreviation (DNA BP, CHAP, ...). A primed name printed nowhere
    on the figure must be found for the reader: its caption sentence says it is not printed."""
    if not cap: return []                            # a missing caption fails in caption_checks
    bad = []
    drawn = set(w.strip().replace(" ", "") for w in (words or []))
    for r in vrows:
        if r["goi"] != g: continue
        lab = r["label"]
        if lab in ("", "HP"): continue                  # never printed
        if (lab.endswith("'") or lab.endswith(" fam.") or r.get("note", "").startswith(FROM_ASSEMBLY)) \
                and lab not in cap:
            bad.append(f"caption does not explain '{lab}'")
        # rebuilt diagonal names lose spaces and can split at a prime or a gap mark, so a name
        # counts as printed when its unprimed stem is found in any rebuilt piece
        stem = lab.rstrip("'").replace(" ", "")
        if words is not None and lab.endswith("'") and not any(stem in x for x in drawn):
            sent = next((x for x in cap.split("\n") if x.startswith(lab + " (")), "")
            if "not printed" not in sent:
                bad.append(f"'{lab}' is printed nowhere, and its caption sentence does not say so")
    for ab in ABBREV:
        if ab.replace(" ", "") in drawn and f"{ab}, " not in cap:
            bad.append(f"abbreviation '{ab}' printed but not explained in the caption")
    return sorted(set(bad))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dpi", type=float, default=1200)
    ap.add_argument("--dir", default="/Data/ST492/figures2026")
    ap.add_argument("--csv", default=None)
    ap.add_argument("--copies", default=None, help="per-copy table (default: beside --csv)")
    ap.add_argument("--loci", default=",".join(LOCI), help="comma-separated loci to check")
    ap.add_argument("--upstream", action="store_true",
                    help="also recount copies and genomes from marker_hits + Bakta")
    ap.add_argument("--bp-per-mm", type=float, default=227.2727,
                    dest="bp_per_mm", help="the scale the plot script claims")
    ap.add_argument("--suffix", default="",
                    help='"" for the set with country breakdowns, "_nocountries" for the other')
    ap.add_argument("--caption-dir", default=None, dest="caption_dir",
                    help="where st492_synteny_<locus>_allctx_caption.txt is (default: --dir)")
    args = ap.parse_args()
    loci = [g.strip() for g in args.loci.split(",") if g.strip()]
    unknown = set(loci) - set(LOCI)
    if unknown: sys.exit(f"unknown loci: {sorted(unknown)}")
    d = pathlib.Path(args.dir)
    tmp = tempfile.mkdtemp(prefix="st492_verify_")
    csv_path = pathlib.Path(args.csv or d / "st492_synteny_variants_allctx8.csv")
    cp_path = pathlib.Path(args.copies or str(csv_path).replace("st492_synteny_variants",
                                                                "st492_synteny_contexts"))
    vrows = read_rows(csv_path)
    crows = read_rows(cp_path) if cp_path.exists() else []
    dr_path = pathlib.Path(str(csv_path).replace("st492_synteny_variants", "st492_synteny_dropped"))
    dr_rows = read_rows(dr_path) if dr_path.exists() else None
    if not crows: print(f"WARNING: no per-copy table at {cp_path}; table checks reduced")
    md = read_rows(MD)
    country_of = {r["Sample"]: (r["country"] or "unknown") for r in md}
    st492 = {r["Sample"] for r in md
             if re.sub(r"\.0$", "", str(r["subtype"]).strip()) == "492"}

    print(f"set: st492_synteny_<locus>_allctx{args.suffix}   loci: {', '.join(loci)}")
    print(f"{'locus':<12}{'tracks':>7}{'exp':>5}{'copies':>8}{'exp':>5}{'gen':>5}{'exp':>5}"
          f"{'bar mm':>9}{'err %':>8}{'page':>6}{'bleed':>7}  result")
    bars, failures, styles = [], [], {}
    for g in loci:
        png = d / f"st492_synteny_{g}_allctx{args.suffix}.png"
        pdf = d / f"st492_synteny_{g}_allctx{args.suffix}.pdf"
        lay_p = d / f"st492_synteny_{g}_allctx_layout.csv"
        if not png.exists() or not pdf.exists():
            failures.append(f"{g}: missing output"); print(f"{g:<12}  MISSING"); continue
        lay = read_rows(lay_p) if (lay_p.exists() and args.suffix == "") else []
        exp, bad = table_checks(g, vrows, crows, lay, country_of)
        # without the display map nothing drawn can be tied to its table row: fail, never skip
        if args.suffix == "" and not lay_p.exists():
            bad.append(f"layout CSV missing ({lay_p.name}): display numbers cannot be checked")

        raw, words = pdf_words(pdf)
        text = " ".join(w.strip() for w in words)
        # track rows only: the notes block under the scale bar also says "variant N"
        boxes = [(float(y), w) for y, w in
                 re.findall(r'<word xMin="[\d.]+" yMin="([\d.]+)"[^>]*>([^<]*)</word>', raw)]
        # the scale-bar label is a bare number followed by "kb" or "bp"; the notes may also
        # say "9,830 bp", so the deepest "bp" on the page is not necessarily the bar
        bar_y = min((y for i, (y, w) in enumerate(boxes)
                     if w.strip() in ("kb", "bp") and i
                     and re.fullmatch(r"[\d.]+", boxes[i - 1][1].strip())), default=float("inf"))
        # the bold "variant N" column only; "same genome as variant N" sits inside a label
        n_tracks = sum(1 for i, (y, w) in enumerate(boxes) if w.strip() == "variant" and y < bar_y
                       and not (i and boxes[i - 1][1].strip() == "as"))
        words_xy = [(float(a), float(b), float(c), float(e), w.strip())
                    for a, b, c, e, w in WORD_RE.findall(raw)]
        labels, dup = drawn_labels(pdf)
        n_copies, n_gen, lbad = label_checks(g, labels, dup, exp, lay)
        bad += lbad
        label_bp = next((float(words[i - 1]) * (1000 if w.strip() == "kb" else 1)
                         for i, w in enumerate(words)
                         if w.strip() in ("kb", "bp") and i
                         and re.fullmatch(r"[\d.]+", words[i - 1].strip())), None)

        im = Image.open(png)
        grey = np.array(im.convert("L"))
        ink = grey < 200
        bleed = bool(ink[0, :].any() or ink[-1, :].any()
                     or ink[:, 0].any() or ink[:, -1].any())
        pw, ph = page_pt(pdf)
        page_ok = (pw is not None and abs(im.width * 72 / args.dpi - pw) <= PAGE_TOL_PT
                   and abs(im.height * 72 / args.dpi - ph) <= PAGE_TOL_PT)

        mm = bar_mm(pdf, tmp)
        bars.append(mm)
        drawn_bp = mm * args.bp_per_mm if mm else float("nan")
        err = (drawn_bp - label_bp) / label_bp * 100 if (mm and label_bp) else float("nan")

        if n_tracks != exp["tracks"]: bad.append("track count")
        bad += chrome_checks(pdf, text, words, bar_y, boxes, tmp)
        if not SHOW_SUBTITLE: bad += subtitle_checks(words_xy)
        bad += font_checks(pdf)
        arrows, glyphs = svg_parts((pathlib.Path(tmp) / (pdf.stem + ".svg")).read_text())
        bad += fill_checks(g, vrows, arrows)
        bad += overprint_checks(arrows, glyphs)
        bad += legend_checks(g, words_xy, vrows, lay)
        bad += marker_checks(g, vrows)
        bad += nrun_checks(g, vrows)
        geo = svg_geometry(pathlib.Path(tmp) / (pdf.stem + ".svg"))
        bad += geometry_checks(geo, pathlib.Path(tmp) / (pdf.stem + ".svg")) + key_checks(geo, words)
        styles[g] = style_metrics(raw, pathlib.Path(tmp) / (pdf.stem + ".svg"), page_pt(pdf))
        if crows: bad += cluster_checks(g, vrows, crows, labels, lay)
        cap_p = pathlib.Path(args.caption_dir or d) / f"st492_synteny_{g}_allctx_caption.txt"
        cap_txt = cap_p.read_text() if cap_p.exists() else ""
        bad += dash_checks(text, cap_txt)
        if args.suffix == "":
            bad += caption_checks(g, cap_p, pdf, labels, exp, page_pt(pdf)[0] or 0)
            bad += glossary_checks(g, cap_txt, vrows, list(words) + sorted(rotated_names(raw)))
        if "n = 1 copies" in text: bad.append("'n = 1 copies' printed")
        if bleed: bad.append("ink at page edge")
        if not page_ok: bad.append(f"page {pw}x{ph} pt vs PNG {im.width}x{im.height} px at {args.dpi:g} dpi")
        if not (abs(err) < BAR_TOL_PCT): bad.append("scale bar")
        up = None
        if args.upstream and g != "carO":
            uc, ug, up, up_raw = upstream_counts(g, st492)
            multi = sorted(f"{s} {c} x{k}" for s, cs in up.items()
                           for c, k in collections.Counter(c for c, _ in cs).items() if k > 1)
            print(f"  {g}: upstream recount {uc} copies in {ug} genomes "
                  f"(drawn track labels {n_copies} copies, {n_gen} genomes)"
                  + (f"; contigs with more than one copy: {', '.join(multi)}" if multi else ""))
            if uc != n_copies: bad.append("upstream copy recount differs")
            # distinct genomes against the per-copy table, also when a genome has copies on two
            # tracks (intI1); the drawn labels, which count such a genome on both, are checked
            # against the per-track count in label_checks
            if exp["genomes"] is not None and ug != exp["genomes"]:
                bad.append(f"upstream genome recount {ug}, per-copy table {exp['genomes']}")
            bad += arrow_upstream_checks(g, vrows, crows)
            if dr_rows is not None:
                n_raw = sum(len(v) for v in up_raw.values())
                n_dr = sum(1 for r in dr_rows if r["goi"] == g)
                if n_raw - uc != n_dr:
                    bad.append(f"{n_raw - uc} hits dropped as duplicate contigs upstream, "
                               f"the dropped table lists {n_dr}")
        if args.suffix == "":
            cap = cap_p
            bad += assembly_checks(g, vrows, cap)
            bad += mobile_checks(g, vrows, cap, lay)
            bad += fragment_caption_checks(g, cap, lay)
            if up is not None: bad += dropped_caption_checks(g, cap, UP_DROPPED.get(g, set()))
        bad += doubled_checks(g, vrows, crows, up, labels)
        if args.suffix == "":
            bad += stack_checks(g, vrows, lay, {w.strip().replace(" ", "") for w in list(words) + sorted(rotated_names(raw))},
                                pathlib.Path(tmp) / (pdf.stem + ".svg"))
            bad += wording_checks(g, cap_p, vrows)
        bad += count_checks(g, n_copies, crows, (uc, ug) if up is not None else None)
        bad += separate_checks(g, crows, lay, labels, up)
        if bad: failures.append(f"{g}: " + "; ".join(bad))
        print(f"{g:<12}{n_tracks:7d}{exp['tracks']:5d}{n_copies:8d}{exp['copies']:5d}"
              f"{n_gen:5d}{(exp['track_genomes'] or -1):5d}{(mm or 0):9.2f}{err:8.2f}"
              f"{'ok' if page_ok else 'BAD':>6}{'YES' if bleed else 'no':>7}"
              f"  {'ok' if not bad else 'FAIL'}")

    shutil.rmtree(tmp, ignore_errors=True)
    # per contig, over every locus in the table
    cf = circle_fact_checks(crows)
    if cf: failures.append("circle facts: " + "; ".join(cf))
    good = [b for b in bars if b]
    if len(good) > 1:
        spread = (max(good) - min(good)) / np.mean(good) * 100
        print(f"\nbar length across the set: {min(good):.2f} to {max(good):.2f} mm, "
              f"spread {spread:.2f} percent")
        if spread > SPREAD_TOL_PCT:
            failures.append(f"scale differs by {spread:.2f} percent across the set")
    if len(styles) > 1:
        for k in sorted({k for st in styles.values() for k in st}):
            print(f"{k:<24}" + "  ".join(f"{g} {st.get(k)}" for g, st in styles.items()))
        failures += [f"set: {x}" for x in set_style_checks(styles)]

    if failures:
        print("\nFAILED:"); [print("  " + f) for f in failures]
        return 1
    print(f"\nall {len(loci)} figure(s) pass")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
