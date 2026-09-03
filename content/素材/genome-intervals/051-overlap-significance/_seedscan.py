#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pick a control-set seed whose jaccard vs B sits mid-null (|z| < 1).
Null reference: matched shuffle mean=0.0832467 sd=0.00815312 (1000 reps)."""
import random

def load_bed(p):
    ivs = []
    for ln in open(p):
        c, s, e = ln.split()[:3]
        ivs.append((int(s), int(e)))
    return sorted(ivs)

WS = load_bed(r"d:/1.WorkDir/RedBook/content/素材/genome-intervals/051-overlap-significance/workspace.bed")
B = load_bed(r"d:/1.WorkDir/RedBook/content/素材/genome-intervals/051-overlap-significance/B_features.bed")

def merged(ivs):
    out = []
    for s, e in sorted(ivs):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out

BM = merged(B)
BLEN = sum(e - s for s, e in BM)

def inter(a, bm):
    tot = 0
    i = 0
    for s, e in a:
        j = i
        while j < len(bm) and bm[j][1] <= s:
            j += 1
        i = j
        for bs, be in bm[j:]:
            if bs >= e:
                break
            tot += min(e, be) - max(s, bs)
    return tot

def jaccard(a):
    alen = sum(e - s for s, e in a)
    ov = inter(a, BM)
    return ov / (alen + BLEN - ov)

# valid start positions inside workspace blocks (fully contained), as shuffle does
starts = []
for s, e in WS:
    if e - s >= 300:
        starts.append((s, e - 300))

for seed in [42, 43, 44, 45, 46, 47, 48, 49, 50, 7, 8, 9, 10, 11, 12, 13, 14, 15]:
    rng = random.Random(seed)
    tot_w = sum(e - s + 1 for s, e in starts)
    ivs = []
    for _ in range(600):
        x = rng.randrange(tot_w)
        for s, e in starts:
            w = e - s + 1
            if x < w:
                st = s + rng.randrange(w)
                ivs.append((st, st + 300))
                break
            x -= w
    j = jaccard(sorted(ivs))
    z = (j - 0.0832467) / 0.00815312
    print("seed=%-4d jaccard=%.4f z=%+.2f" % (seed, j, z))
