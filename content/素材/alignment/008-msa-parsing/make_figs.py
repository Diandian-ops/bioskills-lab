#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
008 msa-parsing 出图：基于 run.py 真跑产出的 msa_parsing_data.json
（保守度剖面直接重读 alignment.aln 计算）。
配色：brick-red #b5482f + celadon #2f7d72。英文标签。_verify 自检。
"""
import os
import json
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import AlignIO

BASE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(BASE, "msa_parsing_data.json")))
aln = AlignIO.read(os.path.join(BASE, "alignment.aln"), "fasta")

BRICK = "#b5482f"
CELADON = "#2f7d72"
GREY = "#6b7280"

def _verify(fig, name):
    fig.canvas.draw(); r = fig.canvas.get_renderer(); fails = 0
    for ax in fig.axes:
        axb = ax.get_window_extent(r); bbs = []
        for t in ax.texts:
            bb = t.get_window_extent(r)
            if (bb.x0 < axb.x0 - 3 or bb.x1 > axb.x1 + 3 or bb.y0 < axb.y0 - 3 or bb.y1 > axb.y1 + 3):
                print("  [FAIL] %s: text out of axes: %r" % (name, t.get_text()[:50])); fails += 1
            bbs.append(bb)
        for i in range(len(bbs)):
            for j in range(i + 1, len(bbs)):
                a, b = bbs[i], bbs[j]
                if not (a.x1 + 1 < b.x0 or b.x1 + 1 < a.x0 or a.y1 + 1 < b.y0 or b.y1 + 1 < a.y0):
                    print("  [FAIL] %s: text overlap: %r / %r" % (name, ax.texts[i].get_text()[:30], ax.texts[j].get_text()[:30])); fails += 1
    if fails == 0: print("  [PASS] %s: text within axes, no overlap" % name)
    return fails

def savefig(fig, name):
    fails = _verify(fig, name)
    fig.savefig(os.path.join(BASE, name), dpi=150); plt.close(fig); return fails

def footnote(fig, text):
    fig.text(0.5, 0.012, text, ha="center", va="bottom", fontsize=8.5, color=GREY)

L = D["n_cols"]; cols = list(range(L))
# 重算每列保守度（忽略空位）
cons_prof = []
for i in range(L):
    col = aln[:, i].replace('-', '')
    cons_prof.append(Counter(col).most_common(1)[0][1] / len(col) if col else 0.0)
gpc = D["gaps_per_column"]

# ============================================================
# 图 1：每列空位比例 + 0.5 阈值
# ============================================================
fig, ax = plt.subplots(figsize=(9.0, 4.8), dpi=150)
bars = ax.bar(cols, [g * 100 for g in gpc], color=CELADON, width=0.8, zorder=3)
ax.axhline(50, color=BRICK, ls="--", lw=1.3, zorder=4)
ax.text(0.4, 52, "gap fraction threshold = 50%", color=BRICK, fontsize=9, ha="left")
ax.set_xlabel("alignment column", fontsize=10.5)
ax.set_ylabel("gap fraction (%)", fontsize=10.5)
ax.set_ylim(0, max(gpc) * 100 * 1.25 + 5)
ax.set_title("Gap fraction per column (msa-parsing: gaps_per_column)\ncols 20-24 exceed 50%% -> removed by remove_gappy_columns",
             fontsize=12, pad=12)
gc = D["gappy_cols"]
if gc:
    mid = (gc[0] + gc[-1]) / 2.0
    ax.text(mid, max(gpc) * 100 + max(gpc) * 100 * 0.06, "gappy cols 20-24",
            ha="center", va="bottom", fontsize=9, color=BRICK)
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "remove_gappy_columns(0.5): 30 columns -> %d columns" % D["cleaned_cols"])
f1 = savefig(fig, "fig1_gap_profile.png")

# ============================================================
# 图 2：Henikoff 权重（近重复序列被降权）
# ============================================================
fig, ax = plt.subplots(figsize=(8.2, 5.0), dpi=150)
ids = D["seq_ids"]; hw = D["henikoff_weights"]
vals = [hw[i] for i in ids]
colors = [BRICK if i in ("s1", "s2") else CELADON for i in ids]
bars = ax.bar(ids, [v * 100 for v in vals], color=colors, width=0.6, zorder=3)
ax.axhline(20, color=GREY, ls=":", lw=1.0, zorder=2)
ax.text(4.3, 20.4, "equal weight = 20%", color=GREY, fontsize=8.5, ha="right")
ax.set_ylabel("Henikoff weight (%)", fontsize=10.5)
ax.set_ylim(0, max(vals) * 100 * 1.3)
ax.set_title("Henikoff sequence weights (msa-parsing: henikoff_weights)\nnear-duplicate s1/s2 down-weighted vs unique s3/s4/s5",
             fontsize=12, pad=12)
for i, v in zip(ids, vals):
    ax.text(i, v * 100 + max(vals) * 100 * 0.02, f"{v*100:.1f}", ha="center", va="bottom", fontsize=9, color="#222222")
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "SKILL.md: weight = sum over columns of 1/(k_c * n_{s,c}); weights sum to 1.0")
f2 = savefig(fig, "fig2_henikoff_weights.png")

# ============================================================
# 图 3：保守度剖面（标记完全保守列 + 空位区）
# ============================================================
fig, ax = plt.subplots(figsize=(9.0, 4.8), dpi=150)
ax.plot(cols, cons_prof, color=CELADON, lw=2.0, zorder=3)
ax.set_xlabel("alignment column", fontsize=10.5)
ax.set_ylabel("conservation (most-common fraction)", fontsize=10.5)
ax.set_ylim(0, 1.08)
ax.set_title("Per-column conservation profile (msa-parsing: find_conserved_positions)\nfully-conserved = %d columns; gappy region 20-24 shaded" % D["n_fully_conserved"],
             fontsize=12, pad=12)
ax.axvspan(20, 24, color=GREY, alpha=0.10, zorder=0)
# 完全保守列标 1.0
fully = [i for i in range(L) if cons_prof[i] >= 0.999]
for i in fully:
    ax.plot(i, 1.0, 'o', color=BRICK, ms=3.5, zorder=4)
ax.grid(axis="y", alpha=0.2, ls=":", lw=0.7, zorder=0)
footnote(fig, "Conserved core (cols 0-9, 25-29) = 100%%; variable region (10-19) ~40%%; gappy cols counted on non-gap chars")
f3 = savefig(fig, "fig3_conservation_profile.png")

total_fail = f1 + f2 + f3
print("FIGURE QUALITY: TOTAL FAILS = %d" % total_fail)
import sys
sys.exit(1 if total_fail else 0)
