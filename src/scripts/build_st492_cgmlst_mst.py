#!/usr/bin/env python3
"""Build an ECDC/GrapeTree-style cgMLST allele minimum spanning tree for ST492.

Reads chewBBACA allele calls (results_alleles.tsv), restricts to the 93 ST492
genomes, keeps the cgMLST core (loci called in >=95% of them), and computes
pairwise ALLELIC distances ignoring loci missing in either isolate (the
Ridom/ECDC convention). Isolates with 0 allele differences are collapsed into a
single profile node (sized by isolate count). A minimum spanning tree is built
on the profile distance matrix (scipy), laid out with branch lengths set by the
square root of the allele distance, and written as node + edge tables for
plot_st492_cgmlst_mst.R (nodes are coloured by fastGEAR cluster).
"""
import os
import re
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree, connected_components

RUN = "/home/bals/Documents/dev/BACTolomew/runs/sample_78eaa9126850/run_78eaa9126850"
ALLELES = f"{RUN}/chewbbaca/out/results_alleles.tsv"
FASTGEAR = f"{RUN}/fastgear/sample_clades.tsv"
METADATA = "/Data/ST492/metadata/metadata_ALL_new_merged_with_clusters_with_complex.csv"
OUTDIR = os.environ.get("ST492_MST_OUTDIR", "/Data/ST492/figures2026")
NODES_CSV = f"{OUTDIR}/st492_cgmlst_mst_nodes.csv"
EDGES_CSV = f"{OUTDIR}/st492_cgmlst_mst_edges.csv"
META_CSV = f"{OUTDIR}/st492_cgmlst_mst_meta.csv"
CORE_CALLRATE = 0.95     # locus kept if called in >= this fraction of ST492 genomes
TREE_W_IN = 6.9          # tree box (with an inside legend) at PRINT size: figure 7.2 in (18 cm) wide
TREE_H_IN = 6.2
LAB_PT = 5.6             # edge-label font size, pt at print size

# ---- ST492 samples + metadata ----
md = pd.read_csv(METADATA)
st492 = sorted(md.loc[md["subtype"].astype(str).str.strip().isin(["492", "492.0"]), "Sample"])
fg = pd.read_csv(FASTGEAR, sep="\t").set_index("Sample")["FastGearCluster"].reindex(st492)
color_of = {s: ("Cluster NA" if pd.isna(c) else f"Cluster {int(c)}") for s, c in fg.items()}

# ---- clean allele calls (INF-x -> x; LNF/NIPH*/PLOT/ALM/ASM/LOTSC -> missing) ----
al = pd.read_csv(ALLELES, sep="\t", dtype=str).rename(columns=lambda c: "Sample" if c.lower() == "file" else c)
sub = al[al["Sample"].isin(st492)].set_index("Sample").loc[st492]


def clean(v):
    if isinstance(v, str):
        if v.startswith("INF-"):
            v = v[4:]
        if re.fullmatch(r"\d+", v):
            return int(v)
    return np.nan


M = sub.map(clean)
callrate = M.notna().mean(axis=0)
core = callrate[callrate >= CORE_CALLRATE].index
Mc = M[core].to_numpy()
n = len(st492)
n_core = len(core)

# ---- pairwise allelic distance (ignore loci missing in either isolate) ----
D = np.zeros((n, n), dtype=float)
for i in range(n):
    a = Mc[i]
    ai = ~np.isnan(a)
    for j in range(i + 1, n):
        b = Mc[j]
        both = ai & ~np.isnan(b)
        d = float(np.sum((a != b) & both))
        D[i, j] = D[j, i] = d

# ---- collapse isolates with 0 allele differences into profiles ----
zero = (D == 0).astype(int)
np.fill_diagonal(zero, 0)
n_prof, lab = connected_components(csr_matrix(zero), directed=False)
profiles = []
for p in range(n_prof):
    members = [st492[i] for i in range(n) if lab[i] == p]
    rep = min(range(n), key=lambda i: (lab[i] != p, i))  # first member index
    rep = next(i for i in range(n) if lab[i] == p)
    profiles.append({"profile": f"P{p+1:02d}", "rep": rep, "members": members})
reps = [pr["rep"] for pr in profiles]
Dp = D[np.ix_(reps, reps)]

# ---- minimum spanning tree on the profile distance matrix ----
mst = minimum_spanning_tree(csr_matrix(Dp))
mst = mst + mst.T
mst = mst.tocoo()
edges = [(int(i), int(j), float(w)) for i, j, w in zip(mst.row, mst.col, mst.data) if i < j]

# ---- unrooted equal-angle layout with SQRT-scaled branch lengths --------------
# Every branch is drawn at BASE0 + sqrt(allele distance): close to proportional,
# so splits in the tree are visible, while 1-allele branches stay long enough for
# their label. Branches above CAP_DIST are drawn at the CAP_DIST length and marked
# as shortened by the plot (a break symbol), since 1453 alleles would dwarf the rest.
np_prof = len(profiles)
BASE0 = 0.4
CAP_DIST = 200
CAP_DRAWN = 130          # a shortened branch is drawn as long as a 130-allele one


def blen_of(dist):
    return BASE0 + np.sqrt(CAP_DRAWN if dist > CAP_DIST else dist)


from collections import defaultdict
tree = defaultdict(list)
blen = {}
for a, b, w in edges:
    bl = blen_of(w)
    tree[a].append(b); tree[b].append(a)
    blen[(a, b)] = bl; blen[(b, a)] = bl
root = max(range(np_prof), key=lambda i: len(tree[i]))
parent = {root: None}
order = []
stack = [root]; seen = {root}
while stack:                                           # pre-order (parent before children)
    v = stack.pop(); order.append(v)
    for u in tree[v]:
        if u not in seen:
            seen.add(u); parent[u] = v; stack.append(u)
ntips = {}
for v in reversed(order):                              # leaves per subtree (post-order)
    kids = [u for u in tree[v] if parent.get(u) == v]
    ntips[v] = 1 if not kids else sum(ntips[u] for u in kids)
pos = np.zeros((np_prof, 2))
wedge = {root: (0.0, 2 * np.pi)}
for v in order:                                        # each child gets a tip-proportional wedge
    a0, a1 = wedge[v]
    kids = [u for u in tree[v] if parent.get(u) == v]
    tot = sum(ntips[u] for u in kids) or 1
    cur = a0
    for u in kids:
        b0, b1 = cur, cur + (a1 - a0) * ntips[u] / tot
        mid = 0.5 * (b0 + b1)
        pos[u] = pos[v] + blen[(v, u)] * np.array([np.cos(mid), np.sin(mid)])
        wedge[u] = (b0, b1)
        cur = b1
pos -= pos.mean(axis=0)

# ---- stress refinement: siblings spread apart, branch lengths kept exact -------
# Equal-angle wedges crowd the leaves of a busy hub into a narrow fan. Stress
# majorization (SMACOF) on the tree path lengths spreads them; the result only
# supplies DIRECTIONS: every branch is then rebuilt from the root at its exact
# drawn length, so lengths stay set by log(1 + allele distance).
TD = np.full((np_prof, np_prof), np.inf)
for src in range(np_prof):
    TD[src, src] = 0.0
    st = [src]
    while st:
        v = st.pop()
        for u in tree[v]:
            if np.isinf(TD[src, u]):
                TD[src, u] = TD[src, v] + blen[(v, u)]; st.append(u)
Wt = np.where(TD > 0, 1.0 / np.where(TD > 0, TD, 1) ** 2, 0.0)
X = pos.copy()
Vm = -Wt.copy(); np.fill_diagonal(Vm, Wt.sum(axis=1))
Vp = np.linalg.pinv(Vm)
for _ in range(400):
    dist = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)
    Bm = -Wt * TD / np.where(dist > 1e-9, dist, 1e-9)
    np.fill_diagonal(Bm, 0); np.fill_diagonal(Bm, -Bm.sum(axis=1))
    Xn = Vp @ Bm @ X
    if np.abs(Xn - X).max() < 1e-5:
        X = Xn; break
    X = Xn
for v in order:
    pv = parent[v]
    if pv is None:
        pos[v] = X[v]
        continue
    d = X[v] - X[pv]
    pos[v] = pos[pv] + blen[(pv, v)] * d / (np.linalg.norm(d) or 1.0)
pos -= pos.mean(axis=0)

# uniform scale (preserves exact length ratios) + node sizes
pos *= 1.5 / blen_of(1)                  # a 1-allele branch is 1.5 units long
counts = np.array([len(pr["members"]) for pr in profiles])
R_BASE = 0.3
radii = R_BASE * (1 + 0.55 * np.sqrt(counts - 1))

# ---- clearance: rotate whole subtrees about their parent until no node touches
# another node or an unrelated branch and no branches cross. A subtree rotation
# keeps every branch length exact. -----------------------------------------------
from math import cos, sin
NN_GAP, NE_GAP = 0.45, 0.35            # node-node gap and node-branch clearance, units
kids_of = defaultdict(list)
for v in order:
    if parent[v] is not None:
        kids_of[parent[v]].append(v)
def _subtree(v):
    out, st = [], [v]
    while st:
        x = st.pop(); out.append(x); st.extend(kids_of[x])
    return out
EA = np.array([a for a, b, _ in edges]); EB = np.array([b for a, b, _ in edges])


def _conflicts(p):
    """List of conflicts as (kind, node or edge, other)."""
    out = []
    d = np.linalg.norm(p[:, None, :] - p[None, :, :], axis=2)
    need = radii[:, None] + radii[None, :] + NN_GAP
    ii, jj = np.where(np.triu(d < need, 1))
    out += [("nn", int(i), int(j)) for i, j in zip(ii, jj)]
    A, B = p[EA], p[EB]; AB = B - A; den = (AB ** 2).sum(1)
    t = np.clip(((p[:, None, :] - A[None]) * AB[None]).sum(2) / den[None], 0, 1)
    dist = np.linalg.norm(p[:, None, :] - (A[None] + t[..., None] * AB[None]), axis=2)
    own = (np.arange(np_prof)[:, None] == EA[None]) | (np.arange(np_prof)[:, None] == EB[None])
    ii, ee = np.where((dist < radii[:, None] + NE_GAP) & ~own)
    out += [("ne", int(i), int(e)) for i, e in zip(ii, ee)]
    for e in range(len(edges)):                              # branch crossings
        for f in range(e + 1, len(edges)):
            if {EA[e], EB[e]} & {EA[f], EB[f]}:
                continue
            p1, p2, p3, p4 = p[EA[e]], p[EB[e]], p[EA[f]], p[EB[f]]
            c = lambda a_, b_, c_: (c_[1] - a_[1]) * (b_[0] - a_[0]) - (b_[1] - a_[1]) * (c_[0] - a_[0])
            if ((c(p3, p4, p1) > 0) != (c(p3, p4, p2) > 0)) and ((c(p1, p2, p3) > 0) != (c(p1, p2, p4) > 0)):
                out.append(("x", e, f))
    return out


def _deeper(e):
    a, b = int(EA[e]), int(EB[e])
    return b if parent.get(b) == a else a


def _rotate(p, v, ang):
    q = p.copy(); c0 = p[parent[v]]
    R = np.array([[cos(ang), -sin(ang)], [sin(ang), cos(ang)]])
    sub = _subtree(v)
    q[sub] = (p[sub] - c0) @ R.T + c0
    return q


for _ in range(300):
    conf = _conflicts(pos)
    if not conf:
        break
    kind, a, b = conf[0]
    if kind == "nn":
        movers = [a, b]
    elif kind == "ne":
        movers = [a, _deeper(b)]
    else:
        movers = [_deeper(a), _deeper(b)]
    movers = [v for v in movers if parent.get(v) is not None]
    best = (len(conf), None)
    for v in movers:
        for k in range(1, 25):
            for sgn in (1, -1):
                q = _rotate(pos, v, sgn * np.deg2rad(4 * k))
                nc = len(_conflicts(q))
                if nc < best[0]:
                    best = (nc, q)
            if best[1] is not None:
                break
    if best[1] is None:
        break
    pos = best[1]
layout_conflicts = _conflicts(pos)
pos -= pos.mean(axis=0)

# ---- rotation + legend spot: try every angle; for each, look for a free spot for
# the legend inside the tree box and keep the angle that draws the tree largest
# (rotations keep every branch length and the planarity; width >= height) -------
LEG_W_IN, LEG_H_IN = 1.6, 1.5     # matches the legend drawn by the plot script
def _bbox(p):
    lo = (p - radii[:, None]).min(axis=0); hi = (p + radii[:, None]).max(axis=0)
    return hi - lo


def _free_spot(p, lw, lh):
    """Lower-left corner of a legend box inside the tree box that holds no node or
    branch, nearest to a corner of the tree box; None if there is none."""
    lo = (p - radii[:, None]).min(axis=0); hi = (p + radii[:, None]).max(axis=0)
    if hi[0] - lo[0] < lw or hi[1] - lo[1] < lh:
        return None
    ts = np.linspace(0, 1, 30)[:, None, None]
    pts = (p[EA][None] + ts * (p[EB] - p[EA])[None]).reshape(-1, 2)
    pad = radii.max()
    xs = np.linspace(lo[0], hi[0] - lw, 24); ys = np.linspace(lo[1], hi[1] - lh, 24)
    best_ = None
    for x0 in xs:
        for y0 in ys:
            inside = ((pts[:, 0] > x0 - pad) & (pts[:, 0] < x0 + lw + pad) &
                      (pts[:, 1] > y0 - pad) & (pts[:, 1] < y0 + lh + pad)).any()
            if inside:
                continue
            dc = min(abs(x0 - lo[0]), abs(hi[0] - lw - x0)) + min(abs(y0 - lo[1]), abs(hi[1] - lh - y0))
            if best_ is None or dc < best_[0]:
                best_ = (dc, x0, y0)
    return None if best_ is None else best_[1:]


best = None
for deg in np.arange(0, 180, 2.0):
    t = np.deg2rad(deg)
    rot = pos @ np.array([[cos(t), sin(t)], [-sin(t), cos(t)]])
    rot = rot - rot.mean(axis=0)
    w_, h_ = _bbox(rot)
    if h_ > w_:
        continue
    sc = min(TREE_W_IN / w_, TREE_H_IN / h_)
    spot = _free_spot(rot, LEG_W_IN / sc, LEG_H_IN / sc)
    if spot is None:                     # legend column right of the tree
        sc = min((TREE_W_IN - LEG_W_IN - 0.15) / w_, TREE_H_IN / h_)
    if best is None or sc > best[0]:
        best = (sc, rot, spot)
unit_in, pos, spot = best
lo = (pos - radii[:, None]).min(axis=0); hi = (pos + radii[:, None]).max(axis=0)
lw, lh = LEG_W_IN / unit_in, LEG_H_IN / unit_in
leg_outside = spot is None
if leg_outside:
    leg_x0, leg_y0 = hi[0] + 0.15 / unit_in, (lo[1] + hi[1]) / 2 - lh / 2
else:
    leg_x0, leg_y0 = spot

# ---- figure scale and edge-label placement -----------------------------------
# The plot draws the tree TREE_W_IN inches wide; label boxes are sized in data
# units from that scale. Each label sits on its own branch, at the free spot
# closest to the midpoint (not on a node, not on another label).
w_, h_ = _bbox(pos)
def _lab_box(txt):
    wpt = 0.62 * LAB_PT * len(txt) + 1.6
    hpt = 1.05 * LAB_PT
    return wpt / 72 / unit_in, hpt / 72 / unit_in


def _box_hits_seg(c, bw, bh, A, B):
    # does the axis-aligned box (centre c) intersect segment AB? (sampled)
    for t in np.linspace(0, 1, 25):
        q = A + t * (B - A)
        if abs(q[0] - c[0]) < bw / 2 and abs(q[1] - c[1]) < bh / 2:
            return True
    return False
placed = []
lab_xy = {}
lab_hits = 0
order_e = sorted(range(len(edges)), key=lambda e: np.linalg.norm(pos[edges[e][0]] - pos[edges[e][1]]))
for e in order_e:
    a, b, w = edges[e]
    bw, bh = _lab_box(str(int(w)))
    best = None
    ab = pos[b] - pos[a]; nrm = np.array([-ab[1], ab[0]]) / (np.linalg.norm(ab) or 1.0)
    cands = [(t, 0.0) for t in (0.5, 0.42, 0.58, 0.34, 0.66, 0.27, 0.73)]
    cands += [(t, o) for t in (0.5, 0.4, 0.6) for o in (1, -1)]      # beside the branch
    if np.linalg.norm(ab) - radii[a] - radii[b] < bw + 0.1:     # no room on the branch
        cands = [(t, o) for t, o in cands if o != 0]
    if w > CAP_DIST:                     # the break symbol sits at the midpoint
        cands = [(t, o) for t, o in cands if abs(t - 0.5) > 0.12]
    for t, off in cands:
        c = pos[a] + t * ab + off * nrm * (0.5 * (bw if abs(nrm[0]) > abs(nrm[1]) else bh) + 0.03)
        cost = 0.02 * abs(off)
        for m in range(np_prof):                       # overlap with nodes
            dx = max(abs(c[0] - pos[m][0]) - bw / 2, 0); dy = max(abs(c[1] - pos[m][1]) - bh / 2, 0)
            if dx * dx + dy * dy < (radii[m] + 0.06) ** 2:
                cost += 1
        for (px_, py_, pw, ph) in placed:              # overlap with labels
            if abs(c[0] - px_) < (bw + pw) / 2 and abs(c[1] - py_) < (bh + ph) / 2:
                cost += 1
        for f, (c2, d2, _) in enumerate(edges):              # overlap with other branches
            if f != e and _box_hits_seg(c, bw, bh, pos[c2], pos[d2]):
                cost += 1
        cost += 0.01 * abs(t - 0.5)
        if best is None or cost < best[0]:
            best = (cost, c)
    lab_hits += int(best[0] >= 1)
    lab_xy[e] = best[1]
    placed.append((best[1][0], best[1][1], bw, bh))

# ---- node table (wide: one count column per colour value) ----
color_vals = sorted(set(color_of.values()))
rows = []
for idx, pr in enumerate(profiles):
    comp = pd.Series([color_of[m] for m in pr["members"]]).value_counts()
    row = {"profile": pr["profile"], "x": pos[idx, 0], "y": pos[idx, 1],
           "radius": float(radii[idx]),
           "n_isolates": len(pr["members"]),
           "dominant": comp.index[0],
           "members": ";".join(pr["members"])}
    for cv in color_vals:
        row[cv] = int(comp.get(cv, 0))
    rows.append(row)
nodes = pd.DataFrame(rows)
nodes.to_csv(NODES_CSV, index=False)

# ---- edge table with endpoint coords + label midpoint ----
erows = []
for k, (a, b, w) in enumerate(edges):
    erows.append({"source": profiles[a]["profile"], "target": profiles[b]["profile"],
                  "dist": int(w), "x": pos[a, 0], "y": pos[a, 1],
                  "xend": pos[b, 0], "yend": pos[b, 1],
                  "mx": lab_xy[k][0], "my": lab_xy[k][1], "capped": int(w > CAP_DIST)})
pd.DataFrame(erows).to_csv(EDGES_CSV, index=False)
pd.DataFrame([{"n_genomes": n, "core_loci": n_core, "core_callrate": CORE_CALLRATE, "n_profiles": np_prof,
               "unit_in": unit_in, "lab_pt": LAB_PT, "leg_x0": leg_x0, "leg_y0": leg_y0,
               "leg_w": lw, "leg_h": lh, "leg_outside": int(leg_outside)}]).to_csv(META_CSV, index=False)

# ---- report ----
iu = np.triu_indices(n, 1)
print(f"ST492 genomes: {n}  |  cgMLST core loci (>= {int(CORE_CALLRATE*100)}% called): {n_core}")
print(f"distinct profiles (0-diff collapsed): {np_prof}  |  MST edges: {len(edges)}")
print(f"allele distances: min={int(D[iu].min())} median={int(np.median(D[iu]))} max={int(D[iu].max())}")
print(f"MST edge distances: min={min(e[2] for e in edges):.0f} median={np.median([e[2] for e in edges]):.0f} max={max(e[2] for e in edges):.0f}")
print(f"regions: {len(color_vals)} values")
print(f"layout conflicts left (node gaps, node on branch, crossings): {len(layout_conflicts)} {layout_conflicts[:5]}")
print(f"edge labels that still touch a node or another label: {lab_hits}")
print(f"tree box: {w_:.2f} x {h_:.2f} units (aspect {w_ / h_:.2f}), {w_ * unit_in:.2f} x {h_ * unit_in:.2f} in")
print(f"legend: {'outside right' if leg_outside else 'inside the tree box'}")
print(f"wrote {NODES_CSV}\nwrote {EDGES_CSV}\nwrote {META_CSV}")
