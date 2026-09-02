#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
010 structural-alignment 出图：基于 structural_results.json（TM-align 真跑产出）。
配色：brick-red #b5482f + celadon #2f7d72。英文标签（DejaVu Sans 无 CJK）。
出图质量自检（替代 check_figs.py）：文本包围盒落 axes 内、不重叠。
"""
import os
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = json.load(open(os.path.join(BASE, "structural_results.json")))

BRICK = "#b5482f"
CELADON = "#2f7d72"
GREY = "#6b7280"

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
    out = os.path.join(BASE, name)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return fails

def footnote(fig, text):
    fig.text(0.5, 0.012, text, ha="center", va="bottom", fontsize=8.5, color=GREY)

pairs = DATA["pairs"]
labels = [p["pair"] for p in pairs]
tm = [p["tm_score_fold_sim"] for p in pairs]
rmsd = [p["rmsd"] for p in pairs]
alen = [p["aligned_length"] for p in pairs]

# ============================================================
# Fig 1: TM-score per pair (fold similarity), with 0.5 threshold
# ============================================================
fig, ax = plt.subplots(figsize=(8.4, 5.0), dpi=150)
bars = ax.bar(labels, tm, color=CELADON, width=0.55, zorder=3)
ax.axhline(0.5, color=BRICK, lw=1.4, ls="--", zorder=4)
ax.text(len(labels) - 1, 0.515, "TM=0.5 same-fold threshold",
        color=BRICK, fontsize=9, ha="right", va="bottom")
ax.set_ylabel("TM-score (max of two length-normalised)", fontsize=10.5)
ax.set_title("TM-align fold similarity: 1UBQ vs 3 structures\n(all below 0.5 -> different folds)",
             fontsize=12, pad=12)
ax.set_ylim(0, max(tm) * 1.5 if max(tm) > 0 else 1)
for x, v in zip(range(len(labels)), tm):
    ax.text(x, v + max(tm) * 0.03, "%.3f" % v, ha="center", va="bottom",
            fontsize=9.5, color="#222222")
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "SKILL.md: fold similarity = larger of the two length-normalised TM-scores (normalised by shorter chain)")
f1 = savefig(fig, "fig1_tm_score.png")

# ============================================================
# Fig 2: RMSD per pair
# ============================================================
fig, ax = plt.subplots(figsize=(8.4, 5.0), dpi=150)
bars = ax.bar(labels, rmsd, color=BRICK, width=0.55, zorder=3)
ax.set_ylabel("RMSD (Angstrom)", fontsize=10.5)
ax.set_title("Superposition RMSD: 1UBQ vs 3 structures\nTM-align optimized RMSD (not evaluation RMSD)",
             fontsize=12, pad=12)
mx = max(rmsd) * 1.25
ax.set_ylim(0, mx)
for x, v in zip(range(len(labels)), rmsd):
    ax.text(x, v + mx * 0.02, "%.2f" % v, ha="center", va="bottom",
            fontsize=9.5, color="#222222")
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "SKILL.md: RMSD alone is misleading (length- and outlier-dependent); always report with TM-score and aligned length")
f2 = savefig(fig, "fig2_rmsd.png")

# ============================================================
# Fig 3: Aligned length per pair
# ============================================================
fig, ax = plt.subplots(figsize=(8.4, 5.0), dpi=150)
bars = ax.bar(labels, alen, color=GREY, width=0.55, zorder=3)
ax.set_ylabel("aligned residues (CA)", fontsize=10.5)
ax.set_title("TM-align aligned length: 1UBQ vs 3 structures\n1UBQ=76 residues; short aligned spans -> low fold similarity",
             fontsize=12, pad=12)
mx = max(alen) * 1.25
ax.set_ylim(0, mx)
for x, v in zip(range(len(labels)), alen):
    ax.text(x, v + mx * 0.02, str(v), ha="center", va="bottom",
            fontsize=9.5, color="#222222")
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "SKILL.md: alignment length matters; short aligned regions inflate apparent RMSD and lower TM-score")
f3 = savefig(fig, "fig3_aligned_len.png")

total_fail = f1 + f2 + f3
print("FIGURE QUALITY: TOTAL FAILS = %d" % total_fail)
import sys
sys.exit(1 if total_fail else 0)
