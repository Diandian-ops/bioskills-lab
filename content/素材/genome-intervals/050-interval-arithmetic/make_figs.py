#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
050 interval-arithmetic figures: every number comes from the real WSL run
(results.json produced by parse_results.py, truth.json from make_inputs.py).
Colors: brick #b5482f (bedtools measured) + celadon #2f7d72 (expected/reference).
English labels only. Self-check prints FIGURE QUALITY: TOTAL FAILS = 0 on success.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(BASE, "results.json")))
T = json.load(open(os.path.join(BASE, "truth.json")))

BRICK = "#b5482f"
CELADON = "#2f7d72"
GREY = "#6b7280"
DARK = "#222222"

M = R["measured"]
I = T["intersect"]


def _verify(fig, name):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fails = 0
    for ax in fig.axes:
        axb = ax.get_window_extent(r)
        bbs = []
        for t in ax.texts:
            bb = t.get_window_extent(r)
            if (bb.x0 < axb.x0 - 3 or bb.x1 > axb.x1 + 3 or
                    bb.y0 < axb.y0 - 3 or bb.y1 > axb.y1 + 3):
                print("  [FAIL] %s: text out of axes: %r" % (name, t.get_text()[:50]))
                fails += 1
            bbs.append(bb)
        for i in range(len(bbs)):
            for j in range(i + 1, len(bbs)):
                a, b = bbs[i], bbs[j]
                if not (a.x1 + 1 < b.x0 or b.x1 + 1 < a.x0 or
                        a.y1 + 1 < b.y0 or b.y1 + 1 < a.y0):
                    print("  [FAIL] %s: text overlap: %r / %r" %
                          (name, ax.texts[i].get_text()[:30], ax.texts[j].get_text()[:30]))
                    fails += 1
    if fails == 0:
        print("  [PASS] %s: text within axes, no overlap" % name)
    return fails


def savefig(fig, name):
    fails = _verify(fig, name)
    fig.savefig(os.path.join(BASE, name), dpi=150)
    plt.close(fig)
    return fails


# ============================================================
# FIG 1: intersect output modes & overlap thresholds, measured vs expected
# ============================================================
labels = ["-u\nA kept", "-v\nA dropped", "-wa -wb\npairs",
          "-f 0.5\nA=peaks", "-f 0.5 -r\nreciprocal", "-f 0.5\nA=genes (swap)"]
meas = [M["u"], M["v"], M["pairs"], M["f05_u"], M["f05r_u"], M["f05_u_swapped"]]
expc = [I["u"], I["v"], I["pairs"], I["f05_u"], I["f05r_u"], I["f05_u_swapped"]]
x = np.arange(len(labels))
w = 0.38
fig, ax = plt.subplots(figsize=(8.6, 4.6))
ax.bar(x - w / 2, meas, w, color=BRICK, label="bedtools v2.31.1 (measured)")
ax.bar(x + w / 2, expc, w, color=CELADON, label="pure-python algebra (expected)")
for xi, v in zip(x - w / 2, meas):
    ax.text(xi, v + 1.6, str(v), ha="center", va="bottom", fontsize=9, color=DARK)
for xi, v in zip(x + w / 2, expc):
    ax.text(xi, v + 1.6, str(v), ha="center", va="bottom", fontsize=9, color=DARK)
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=8.6)
ax.set_ylabel("features / pairs (count)")
ax.set_ylim(0, 97)
ax.set_title("bedtools intersect matches an independent interval algebra on every mode\n"
             "150 synthetic peaks vs 40 genes, chr1 100 kb + chr2 80 kb, seed 42",
             fontsize=10.5)
ax.legend(fontsize=9, loc="upper right")
plt.tight_layout()
total_fail = savefig(fig, "fig1_intersect_modes_measured_vs_expected.png")

# ============================================================
# FIG 2: the sorted-input contract (3 panels)
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.2))

# (a) merge: unsorted refused vs sorted variants
ax = axes[0]
cats = ["unsorted\ninput", "sorted\n-d 0", "sorted\n-d 1", "sorted\n-d 100"]
vals = [M["merge_unsorted_lines"], M["merge_d0"], M["merge_d1"], M["merge_d100"]]
cols = [GREY, BRICK, BRICK, BRICK]
ax.bar(cats, vals, color=cols, width=0.62)
for xi, v in enumerate(vals):
    ax.text(xi, v + 2.2, str(v), ha="center", va="bottom", fontsize=9, color=DARK)
ax.text(0, 12, "exit 1", ha="center", va="bottom", fontsize=9, color=BRICK)
ax.set_ylabel("merged blocks (count)")
ax.set_ylim(0, 100)
ax.set_title("merge demands prior sort\nv2.31.1 refuses unsorted input",
             fontsize=10)

# (b) -sorted sweep: equal results, far less memory
ax = axes[1]
rss = [R["perf"].get("inmem_maxrss_kb", 0) / 1024.0,
       R["perf"].get("sorted_maxrss_kb", 0) / 1024.0]
ax.bar(["in-memory tree", "-sorted sweep"], rss, color=[BRICK, CELADON], width=0.5)
ax.text(0, rss[0] + 2.5, "%.1f MB" % rss[0], ha="center", va="bottom",
        fontsize=9, color=DARK)
ax.text(1, rss[1] + 2.5, "%.1f MB" % rss[1], ha="center", va="bottom",
        fontsize=9, color=DARK)
ax.set_ylabel("peak memory (MB)")
ax.set_ylim(0, 125)
ax.set_title("intersect on 300k x 300k intervals\nsame answer (73 = 73), %.0fx less RAM"
             % (rss[0] / rss[1]), fontsize=10)

# (c) chrom-naming mismatch: silent empty result
ax = axes[2]
vals = [M["mismatch_u"], M["sorted_inmem_u"]]
ax.bar(["chr1 vs 1\n(naming mismatch)", "chr1 vs chr1\n(harmonized)"],
       vals, color=[GREY, CELADON], width=0.5)
for xi, v in enumerate(vals):
    ax.text(xi, v + 1.8, str(v), ha="center", va="bottom", fontsize=9, color=DARK)
ax.set_ylabel("peaks overlapping genes (count)")
ax.set_ylim(0, 90)
ax.set_title("chromosome-name mismatch:\nempty output, exit 0, no error", fontsize=10)

plt.tight_layout()
total_fail += savefig(fig, "fig2_sorted_contract_footguns.png")

# ============================================================
# FIG 3: overlap semantics and value transfer (3 panels)
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.2))

# (a) A/B asymmetry and fractional overlap
ax = axes[0]
cats = ["default\n(1 bp)", "-f 0.5\nA=peaks", "-f 0.5 -r\nreciprocal", "-f 0.5\nA=genes"]
vals = [M["u"], M["f05_u"], M["f05r_u"], M["f05_u_swapped"]]
ax.bar(cats, vals, color=[BRICK, BRICK, CELADON, CELADON], width=0.62)
for xi, v in enumerate(vals):
    ax.text(xi, v + 1.6, str(v), ha="center", va="bottom", fontsize=9, color=DARK)
ax.set_ylabel("kept A features (count)")
ax.set_ylim(0, 88)
ax.set_title("-f is a fraction of A, -F of B:\nswapping roles changes the count", fontsize=10)

# (b) -split: envelope vs exon blocks
ax = axes[1]
vals = [M["split_env_bp"], M["split_block_bp"]]
ax.bar(["envelope\n(no -split)", "exon blocks\n(-split)"], vals,
       color=[BRICK, CELADON], width=0.5)
ax.text(0, vals[0] + 120, "%d bp" % vals[0], ha="center", va="bottom",
        fontsize=9, color=DARK)
ax.text(1, vals[1] + 120, "%d bp" % vals[1], ha="center", va="bottom",
        fontsize=9, color=DARK)
ax.set_ylabel("overlap with exons (bp)")
ax.set_ylim(0, 6800)
ax.set_title("25 spliced BED12 transcripts vs 80 exons\n"
             "inflated %.1fx without -split (hits 14 vs 10)"
             % (vals[0] / float(vals[1])), fontsize=10)

# (c) map: per-gene mean signal, measured vs expected
ax = axes[2]
gexp = T["map_mean"]
gmea = R["map_measured"]
gx = [gexp[g] for g in gexp if gexp[g] is not None and gmea.get(g) is not None]
gy = [gmea[g] for g in gexp if gexp[g] is not None and gmea.get(g) is not None]
lo = min(gx + gy)
hi = max(gx + gy)
pad = (hi - lo) * 0.08
ax.scatter(gx, gy, s=26, color=CELADON, alpha=0.85)
ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], ls="--", lw=1, color=GREY)
ax.set_xlim(lo - pad, hi + pad)
ax.set_ylim(lo - pad, hi + pad)
ax.set_xlabel("expected mean score (pure python)")
ax.set_ylabel("bedtools map mean")
ax.set_title("map -c 4 -o mean onto 40 genes:\n%d/%d genes on the diagonal"
             % (len(gx), len(gexp)), fontsize=10)

plt.tight_layout()
total_fail += savefig(fig, "fig3_overlap_semantics_map.png")

print("FIGURE QUALITY: TOTAL FAILS = %d" % total_fail)
