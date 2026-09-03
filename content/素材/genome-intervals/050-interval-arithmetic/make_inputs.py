#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
050 interval-arithmetic synthetic BED dataset (pure stdlib, deterministic seed=42).

Genome: chr1 100000 bp + chr2 80000 bp (genome.txt). All intervals 0-based
half-open, same convention bedtools uses.

Sets written next to this script:
  peaks.bed      150 ChIP-like peaks, 200-800 bp, BED5; FILE ORDER SHUFFLED
                 (deliberately unsorted -> merge / -sorted footgun demos)
  genes.bed      40 non-overlapping gene bodies, 1000-3000 bp, BED6
  blacklist.bed  25 regions, 500-2000 bp, BED3
  exons.bed      80 exons, 300 bp, BED4
  spliced.bed12  25 BED12 transcripts, span 2000 bp, 3 x 150 bp blocks
                 (block starts 0/925/1850, introns ~775 bp)
  rep1/2/3.bed   3 replicate peak sets sampled from 60 shared true loci
                 (jitter +/-40 bp) + 8 private peaks each, BED5
  scores.bedgraph     dense signal grid (250 bp bin every 500 bp) for map
  ubg1/2/3.bedgraph   same grid, per-sample presence p=0.7, for unionbedg
  genome.txt          chromosome sizes

truth.json: expected values computed by an INDEPENDENT pure-python interval
algebra implemented here (sort/sweep), used later to reconcile bedtools output.
"""
import json
import os
import random

random.seed(42)
BASE = os.path.dirname(os.path.abspath(__file__))

CHROMS = [("chr1", 100000), ("chr2", 80000)]
GENOME_BP = sum(l for _, l in CHROMS)

# ---------- independent interval algebra (ground truth) ----------

def merge(ivs, d=0):
    """genomic-order merge; features merge when gap <= d (book-ended at d=0)."""
    out = []
    for c, s, e in sorted(ivs):
        if out and out[-1][0] == c and s - out[-1][2] <= d:
            if e > out[-1][2]:
                out[-1][2] = e
        else:
            out.append([c, s, e])
    return [(c, s, e) for c, s, e in out]


# NOTE: the SKILL's "merge on unsorted input silently under-merges" describes
# older bedtools; v2.31.1 refuses unsorted input with exit 1 (verified in _run.sh).


def cov(ivs):
    return sum(e - s for _, s, e in merge(ivs))


def per_chrom(ivs):
    d = {}
    for c, s, e in ivs:
        d.setdefault(c, []).append((s, e))
    return d


def overlap_stats(A, B):
    """(n A-features with >=1 B overlap, n overlapping pairs, total overlap bp)"""
    Bd = per_chrom(B)
    nA = npairs = bp = 0
    for c, s, e in A:
        hit = False
        for s2, e2 in Bd.get(c, []):
            o = min(e, e2) - max(s, s2)
            if o > 0:
                hit = True
                npairs += 1
                bp += o
        nA += 1 if hit else 0
    return nA, npairs, bp


def frac_hits(A, B, f):
    """A features having >=1 single B with per-pair overlap >= f*len(A)"""
    Bd = per_chrom(B)
    n = 0
    for c, s, e in A:
        for s2, e2 in Bd.get(c, []):
            if min(e, e2) - max(s, s2) >= f * (e - s):
                n += 1
                break
    return n


def recip_hits(A, B, f):
    """A features having >=1 single B with overlap >= f*len(A) AND >= f*len(B)"""
    Bd = per_chrom(B)
    n = 0
    for c, s, e in A:
        for s2, e2 in Bd.get(c, []):
            o = min(e, e2) - max(s, s2)
            if o >= f * (e - s) and o >= f * (e2 - s2):
                n += 1
                break
    return n


def complement(A, chroms):
    M = merge(A)
    Md = {}
    for c, s, e in M:
        Md.setdefault(c, []).append((s, e))
    gaps = []
    for c, ln in chroms:
        pos = 0
        for s, e in Md.get(c, []):
            if s > pos:
                gaps.append((c, pos, s))
            pos = max(pos, e)
        if pos < ln:
            gaps.append((c, pos, ln))
    return gaps


def multiinter_all3(sets, chroms):
    """(maximal interval count, bp) covered by ALL sets simultaneously
    (multiinter reports maximal runs, so micro-segments are re-merged)"""
    micro = []
    for c, _ln in chroms:
        lists = []
        bounds = set()
        for ivs in sets:
            L = sorted((s, e) for cc, s, e in ivs if cc == c)
            lists.append(L)
            for s, e in L:
                bounds.add(s)
                bounds.add(e)
        pts = sorted(bounds)
        for s, e in zip(pts, pts[1:]):
            mid = (s + e) / 2.0
            if all(any(a <= mid < b for a, b in L) for L in lists):
                micro.append((c, s, e))
    merged = merge(micro, 0)
    return len(merged), sum(e - s for _c, s, e in merged)


def jaccard(A, B):
    """bedtools jaccard on MERGED inputs: inter/union + n_intersections"""
    mA, mB = merge(A), merge(B)
    _nA, n_inter, inter_bp = overlap_stats(mA, mB)
    covA, covB = cov(mA), cov(mB)
    union_bp = covA + covB - inter_bp
    return inter_bp, union_bp, inter_bp / union_bp, n_inter


def map_means(A_named, B_vals):
    """per A (name given): unweighted mean of B scores of overlapping B feats."""
    B = [(c, s, e, v) for c, s, e, v in B_vals]
    Bd = {}
    for c, s, e, v in B:
        Bd.setdefault(c, []).append((s, e, v))
    out = {}
    for c, s, e, name in A_named:
        vals = [v for s2, e2, v in Bd.get(c, []) if min(e, e2) - max(s, s2) > 0]
        out[name] = (sum(vals) / len(vals)) if vals else None
    return out


# ---------- generate features ----------

def rplace(n, minlen, maxlen, allow_overlap=True, avoid=None):
    ivs = []
    tries = 0
    while len(ivs) < n:
        tries += 1
        assert tries < 200000, "placement failed"
        c, ln = random.choice(CHROMS)
        L = random.randint(minlen, maxlen)
        s = random.randint(0, ln - L)
        e = s + L
        if avoid is not None:
            bad = False
            for c2, s2, e2 in avoid:
                if c2 == c and s < e2 and e > s2:
                    bad = True
                    break
            if bad:
                continue
        ivs.append((c, s, e))
    return ivs


# peaks: BED5, shuffled file order
peaks = rplace(150, 300, 1200)
peaks_named = [(c, s, e, "peak%03d" % i, random.randint(1, 1000))
               for i, (c, s, e) in enumerate(peaks)]
peaks_rows = peaks_named[:]
random.shuffle(peaks_rows)          # footgun: file order != genomic order

# genes: BED6, non-overlapping
genes_iv = rplace(40, 1000, 3000, allow_overlap=False)
genes_named = [(c, s, e, "gene%03d" % i, random.randint(1, 100),
                random.choice("+-"))
               for i, (c, s, e) in enumerate(genes_iv)]
genes_rows = sorted(genes_named, key=lambda r: (r[0], r[1]))

blacklist = rplace(25, 500, 2000)
exons = rplace(80, 300, 300)
exons_named = [(c, s, e, "exon%03d" % i) for i, (c, s, e) in enumerate(exons)]

# BED12 spliced transcripts
BLOCK_STARTS = [0, 925, 1850]
BLOCK_LEN = 150
tx_rows, tx_blocks = [], []
for i in range(25):
    c, ln = random.choice(CHROMS)
    s = random.randint(0, ln - 2001)
    e = s + 2000
    tx_rows.append((c, s, e, "tx%03d" % i, 0, random.choice("+-"),
                    s, e, "0", 3,
                    ",".join(str(BLOCK_LEN) for _ in BLOCK_STARTS),
                    ",".join(str(BLOCK_STARTS[k]) for k in range(3)),
                    ))
    for bs in BLOCK_STARTS:
        tx_blocks.append((c, s + bs, s + bs + BLOCK_LEN, "tx%03d" % i))

# replicate peak sets from 60 shared true loci
true_loci = rplace(60, 400, 600)
rep_rows = {}
for r in range(1, 4):
    rows = []
    for c, s, e in random.sample(true_loci, 55):
        s2 = max(0, s + random.randint(-40, 40))
        rows.append((c, s2, s2 + (e - s), "rep%d_locus" % r, random.randint(100, 999)))
    for k in range(8):
        c, s, e = rplace(1, 300, 500)[0]
        rows.append((c, s, e, "rep%d_private%d" % (r, k), random.randint(100, 999)))
    rows.sort(key=lambda x: (x[0], x[1]))
    rep_rows[r] = rows

# bedGraph signal grid: 250 bp bin every 500 bp
grid = []
for c, ln in CHROMS:
    for k in range(ln // 500):
        s = k * 500
        grid.append((c, s, s + 250))
scores_rows = [(c, s, e, random.randint(0, 100)) for c, s, e in grid]
ubg_rows = {}
for r in (1, 2, 3):
    ubg_rows[r] = [(c, s, e, random.randint(0, 100))
                   for c, s, e in grid if random.random() < 0.7]

# ---------- write files ----------

def w(name, lines):
    with open(os.path.join(BASE, name), "w") as f:
        f.write("\n".join(lines) + "\n")

w("genome.txt", ["%s\t%d" % (c, l) for c, l in CHROMS])
w("peaks.bed", ["%s\t%d\t%d\t%s\t%d" % r for r in peaks_rows])
w("genes.bed", ["%s\t%d\t%d\t%s\t%d\t%s" % r for r in genes_rows])
w("blacklist.bed", ["%s\t%d\t%d" % r for r in blacklist])
w("exons.bed", ["%s\t%d\t%d\t%s" % r for r in exons_named])
w("spliced.bed12", ["%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s" % r
                    for r in tx_rows])
for r in (1, 2, 3):
    w("rep%d.bed" % r, ["%s\t%d\t%d\t%s\t%d" % x for x in rep_rows[r]])
w("scores.bedgraph", ["%s\t%d\t%d\t%d" % x for x in scores_rows])
for r in (1, 2, 3):
    w("ubg%d.bedgraph" % r, ["%s\t%d\t%d\t%d" % x for x in ubg_rows[r]])

# ---------- expected values (independent implementation) ----------

peaks_iv = [(c, s, e) for c, s, e, _n, _sc in peaks_named]
genes_iv_named = [(c, s, e, n) for c, s, e, n, _sc, _st in genes_named]
genes_iv = [(c, s, e) for c, s, e, _n in genes_iv_named]

u, npairs, wo_bp = overlap_stats(peaks_iv, genes_iv)
peaks_cov = cov(peaks_iv)
genes_cov = cov(genes_iv)

merged_peaks_d0 = merge(peaks_iv, 0)
comp = complement(peaks_iv, CHROMS)
# per-feature clip: bedtools subtract keeps A's fragments separate (no
# cross-feature merge), so the expected interval count is summed per peak
_bl = per_chrom(blacklist)
res_bp = res_ivs = 0
for c, s, e in peaks_iv:
    cur = [(s, e)]
    for s2, e2 in _bl.get(c, []):
        nxt = []
        for cs, ce in cur:
            if e2 <= cs or s2 >= ce:
                nxt.append((cs, ce))
                continue
            if s2 > cs:
                nxt.append((cs, s2))
            if e2 < ce:
                nxt.append((e2, ce))
        cur = nxt
    res_bp += sum(b - a for a, b in cur)
    res_ivs += len(cur)
subA_dropped = sum(1 for c, s, e in peaks_iv
                   if any(cc == c and s < e2 and e > s2 for cc, s2, e2 in blacklist))

inter_bp, union_bp, jac, n_inter = jaccard(peaks_iv, genes_iv)
n3, bp3 = multiinter_all3(
    [[(c, s, e) for c, s, e, _n, _sc in rep_rows[r]] for r in (1, 2, 3)], CHROMS)

tx_env = [(c, s, e) for c, s, e, *_ in tx_rows]
env_n, _p, env_bp = overlap_stats(tx_env, exons)
# -split -u counts TRANSCRIPTS with >=1 overlapping block, not blocks
blk_names = set()
blk_bp = 0
Ed = per_chrom(exons)
for c, s, e, name in tx_blocks:
    hit = False
    for s2, e2 in Ed.get(c, []):
        o = min(e, e2) - max(s, s2)
        if o > 0:
            blk_bp += o
            hit = True
    if hit:
        blk_names.add(name)
blk_n = len(blk_names)

truth = dict(
    seed=42,
    genome=dict(chroms=dict(CHROMS), total_bp=GENOME_BP),
    counts=dict(peaks=150, genes=40, blacklist=25, exons=80, transcripts=25,
                true_loci=60, grid_bins=len(grid)),
    coverage=dict(peaks_bp=peaks_cov, genes_bp=genes_cov,
                  peaks_frac=round(peaks_cov / GENOME_BP, 4)),
    intersect=dict(u=u, v=150 - u, pairs=npairs, wo_bp=wo_bp,
                   loj_lines=npairs + (150 - u),
                   wao_lines=npairs + (150 - u),
                   f05_u=frac_hits(peaks_iv, genes_iv, 0.5),
                   f05r_u=recip_hits(peaks_iv, genes_iv, 0.5),
                   f05_u_swapped=frac_hits(genes_iv, peaks_iv, 0.5)),
    subtract=dict(residual_bp=res_bp,
                  residual_ivs=res_ivs,
                  minusA_dropped=subA_dropped,
                  minusA_remaining=150 - subA_dropped),
    complement=dict(bp=GENOME_BP - peaks_cov, n=len(comp)),
    merge=dict(sorted_d0=len(merged_peaks_d0),
               sorted_d1=len(merge(peaks_iv, 1)),
               sorted_d100=len(merge(peaks_iv, 100)),
               sorted_d0_bp=peaks_cov),
    cluster=dict(n=len(merged_peaks_d0)),
    map_mean=map_means(genes_iv_named, scores_rows),
    multiinter=dict(n3=n3, bp3=bp3),
    unionbedg=dict(rows=len(set((c, s, e) for r in (1, 2, 3)
                                for c, s, e, _v in ubg_rows[r])),
                   present1=len(ubg_rows[1]), present2=len(ubg_rows[2]),
                   present3=len(ubg_rows[3])),
    jaccard=dict(inter_bp=inter_bp, union_bp=union_bp, jaccard=jac,
                 n_intersections=n_inter),
    split=dict(env_u=env_n, env_bp=env_bp, block_u=blk_n, block_bp=blk_bp),
    chrom_mismatch=dict(u=0),
)

# per-gene overlap bp sum (groupby expectation)
Bd = per_chrom(peaks_iv)
gsum = {}
for c, s, e, n in genes_iv_named:
    gsum[n] = sum(min(e, e2) - max(s, s2) for s2, e2 in Bd.get(c, [])
                  if min(e, e2) - max(s, s2) > 0)
truth["groupby_sum_bp"] = gsum

with open(os.path.join(BASE, "truth.json"), "w") as f:
    json.dump(truth, f, indent=2)

print("inputs written: peaks=%d genes=%d blacklist=%d exons=%d tx=%d "
      "grid_bins=%d" % (150, 40, 25, 80, 25, len(grid)))
print("expected: intersect -u=%d -v=%d pairs=%d wo_bp=%d "
      "f05=%d f05r=%d swapped=%d" % (u, 150 - u, npairs, wo_bp,
                                     truth["intersect"]["f05_u"],
                                     truth["intersect"]["f05r_u"],
                                     truth["intersect"]["f05_u_swapped"]))
print("expected: merge sorted d0=%d d1=%d d100=%d "
      "complement_bp=%d subtract_residual_bp=%d minusA_dropped=%d"
      % (truth["merge"]["sorted_d0"],
         truth["merge"]["sorted_d1"], truth["merge"]["sorted_d100"],
         truth["complement"]["bp"], truth["subtract"]["residual_bp"],
         subA_dropped))
print("expected: multiinter all3 n=%d bp=%d unionbedg_rows=%d "
      "jaccard=%.6f n_inter=%d split env_bp=%d block_bp=%d"
      % (n3, bp3, truth["unionbedg"]["rows"], jac, n_inter, env_bp, blk_bp))
