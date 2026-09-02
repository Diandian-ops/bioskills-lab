#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
007 multiple-alignment 出图：基于 run.py 真跑产出的 msa_data.json。
配色：brick-red #b5482f + celadon #2f7d72。英文标签（DejaVu Sans 无 CJK）。
出图质量自检：文本包围盒落 axes 内、不重叠；末尾打印 FIGURE QUALITY: TOTAL FAILS = N。
"""
import os
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = json.load(open(os.path.join(BASE, "msa_data.json")))

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

# ===== 图 1：各工具比对长度 =====
tools = DATA["tools"]
labels = [t["name"] for t in tools]
lengths = [t["len"] for t in tools]
fig, ax = plt.subplots(figsize=(8.4, 5.0), dpi=150)
cols = [CELADON, CELADON, CELADON, BRICK, GREY]
bars = ax.bar(labels, lengths, color=cols, zorder=3)
ax.set_ylabel("alignment length (columns)", fontsize=10.5)
ax.set_title("Alignment length per tool\n6 human S1 serine proteases (247-304 aa)", fontsize=12, pad=12)
mx = max(lengths) * 1.18
ax.set_ylim(0, mx)
ax.set_xticks(range(len(labels)))
for b, v in zip(bars, lengths):
    ax.text(b.get_x() + b.get_width() / 2, v + mx * 0.015, str(v),
            ha="center", va="bottom", fontsize=9.5, color="#222222")
ax.set_xticklabels(labels, rotation=15, ha="right", fontsize=9)
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "Brick = MUSCLE5; grey = ClustalOmega. Real outputs from mafft/muscle/clustalo on the same 6 sequences.")
f1 = savefig(fig, "fig1_alignment_length.png")

# ===== 图 2：成对一致性热图 (mafft_linsi) =====
im = DATA["identity_matrix"]
lab = im["labels"]
M = np.array(im["matrix"], dtype=float)
fig, ax = plt.subplots(figsize=(6.6, 5.6), dpi=150)
cax = ax.imshow(M, cmap="YlOrRd", vmin=70, vmax=100)
ax.set_xticks(range(len(lab))); ax.set_yticks(range(len(lab)))
ax.set_xticklabels(lab, rotation=45, ha="right", fontsize=8.5)
ax.set_yticklabels(lab, fontsize=8.5)
ax.set_title("Pairwise identity matrix (PID2 %)\nmafft L-INS-i", fontsize=12, pad=10)
for i in range(len(lab)):
    for j in range(len(lab)):
        ax.text(j, i, "%.1f" % M[i, j], ha="center", va="center",
                fontsize=7.5, color="#222222")
cb = fig.colorbar(cax, ax=ax, fraction=0.046, pad=0.04)
cb.set_label("PID2 %", fontsize=9)
footnote(fig, "Diagonal = 100. High conservation across these cytochrome c orthologs (gene CYCS).")
f2 = savefig(fig, "fig2_identity_heatmap.png")

# ===== 图 3：逐列 gap 比例 (mafft_linsi) =====
gp = DATA["gap_profile"]["mafft_linsi"]
x = list(range(1, len(gp) + 1))
fig, ax = plt.subplots(figsize=(8.4, 4.6), dpi=150)
ax.plot(x, [g * 100 for g in gp], color=BRICK, lw=1.4, zorder=3)
ax.fill_between(x, [g * 100 for g in gp], color=BRICK, alpha=0.12)
ax.axhline(50, color=GREY, ls="--", lw=1.0, zorder=2)
ax.text(len(x) * 0.02, 52, "gap>=50% threshold", fontsize=8.5, color=GREY)
ax.set_xlabel("alignment column", fontsize=10.5)
ax.set_ylabel("gap fraction (%)", fontsize=10.5)
ax.set_title("Per-column gap fraction (mafft L-INS-i)\n6 cytochrome c", fontsize=12, pad=10)
ax.set_ylim(0, max(max(gp) * 100 * 1.2, 5))
ax.grid(alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "Few gappy columns: conserved orthologs align with minimal indels.")
f3 = savefig(fig, "fig3_gap_profile.png")

total_fail = f1 + f2 + f3
print("FIGURE QUALITY: TOTAL FAILS = %d" % total_fail)
import sys
sys.exit(1 if total_fail else 0)
