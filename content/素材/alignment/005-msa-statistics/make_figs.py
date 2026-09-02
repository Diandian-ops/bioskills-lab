#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
005 msa-statistics 出图：基于 run.py 真跑产出的 msa_statistics_data.json。
配色：brick-red #b5482f + celadon #2f7d72。英文标签。_verify 自检文本越界/重叠。
"""
import os
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

BASE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(BASE, "msa_statistics_data.json")))

BRICK = "#b5482f"
CELADON = "#2f7d72"
GREY = "#6b7280"
CMAP = LinearSegmentedColormap.from_list("brickceladon", ["#e9f3f1", CELADON, "#caa89f", BRICK])

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
                print("  [FAIL] %s: text out of axes: %r" % (name, t.get_text()[:50])); fails += 1
            bbs.append(bb)
        for i in range(len(bbs)):
            for j in range(i + 1, len(bbs)):
                a, b = bbs[i], bbs[j]
                if not (a.x1 + 1 < b.x0 or b.x1 + 1 < a.x0 or a.y1 + 1 < b.y0 or b.y1 + 1 < a.y0):
                    print("  [FAIL] %s: text overlap: %r / %r" % (name, ax.texts[i].get_text()[:30], ax.texts[j].get_text()[:30])); fails += 1
    if fails == 0:
        print("  [PASS] %s: text within axes, no overlap" % name)
    return fails

def savefig(fig, name):
    fails = _verify(fig, name)
    fig.savefig(os.path.join(BASE, name), dpi=150); plt.close(fig)
    return fails

def footnote(fig, text):
    fig.text(0.5, 0.012, text, ha="center", va="bottom", fontsize=8.5, color=GREY)

L = D["n_cols"]
cols = list(range(L))
cons = D["conservation_profile"]
ent = D["entropy"]

# ============================================================
# 图 1：保守度曲线（左轴）+ 香农熵曲线（右轴）
# ============================================================
fig, ax = plt.subplots(figsize=(9.2, 5.0), dpi=150)
ax.plot(cols, cons, color=CELADON, lw=2.0, zorder=3, label="conservation")
ax.set_xlabel("alignment column", fontsize=10.5)
ax.set_ylabel("conservation (most-common fraction)", fontsize=10, color=CELADON)
ax.set_ylim(0, 1.05)
ax.tick_params(axis="y", labelcolor=CELADON)
ax2 = ax.twinx()
ax2.plot(cols, ent, color=BRICK, lw=2.0, zorder=3, label="entropy")
ax2.set_ylabel("Shannon entropy (bits)", fontsize=10, color=BRICK)
ax2.set_ylim(0, max(ent) * 1.25)
ax2.tick_params(axis="y", labelcolor=BRICK)
ax.set_title("Conservation vs Shannon entropy across 40 columns\n(msa-statistics: column_conservation + shannon_entropy)",
             fontsize=12, pad=12)
# 空位区标出
for gc in D["gappy_cols"]:
    ax.axvline(gc, color=GREY, ls=":", lw=0.8, zorder=1)
ax.axvspan(28, 31, color=GREY, alpha=0.08, zorder=0)
ax.grid(axis="x", alpha=0.2, ls=":", lw=0.7, zorder=0)
footnote(fig, "Shaded columns 28-31 are the >50%% gap region. Conservation ignores gaps; entropy uses non-gap chars")
f1 = savefig(fig, "fig1_conservation_entropy.png")

# ============================================================
# 图 2：成对一致性矩阵（PID2）热力图
# ============================================================
fig, ax = plt.subplots(figsize=(6.6, 5.8), dpi=150)
ids = D["ids"]
M = [[D["identity_matrix"][i][j] for j in ids] for i in ids]
im = ax.imshow(M, cmap=CMAP, vmin=0, vmax=100)
ax.set_xticks(range(len(ids))); ax.set_yticks(range(len(ids)))
ax.set_xticklabels(ids, fontsize=9); ax.set_yticklabels(ids, fontsize=9)
ax.set_title("Pairwise identity matrix (PID2)\n6-sequence protein MSA", fontsize=12, pad=12)
for i in range(len(ids)):
    for j in range(len(ids)):
        v = M[i][j]
        ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=8.5,
                color="#ffffff" if (v < 35 or v > 80) else "#222222")
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cb.set_label("percent identity", fontsize=9)
ax.set_xlabel("sequence", fontsize=9.5); ax.set_ylabel("sequence", fontsize=9.5)
footnote(fig, "SKILL.md: PID2 = matches over aligned residue pairs only (highest of the four denominators)")
f2 = savefig(fig, "fig2_identity_matrix.png")

# ============================================================
# 图 3：每列空位比例 + 0.5 阈值
# ============================================================
fig, ax = plt.subplots(figsize=(9.2, 4.8), dpi=150)
gp = D["gap_profile"]
bars = ax.bar(cols, [g * 100 for g in gp], color=CELADON, width=0.8, zorder=3)
ax.axhline(50, color=BRICK, ls="--", lw=1.3, zorder=4)
ax.text(0.5, 52, "gap fraction threshold = 50%", color=BRICK, fontsize=9, ha="left")
ax.set_xlabel("alignment column", fontsize=10.5)
ax.set_ylabel("gap fraction (%)", fontsize=10.5)
ax.set_ylim(0, max(gp) * 100 * 1.2 + 5)
ax.set_title("Gap fraction per column (msa-statistics: gap_profile)\ncolumns 28-31 exceed 50%% -> flagged as gappy",
             fontsize=12, pad=12)
# 标出 >50% 的区域（一次性标注，避免相邻列文字重叠）
gc = D["gappy_cols"]
if gc:
    mid = (gc[0] + gc[-1]) / 2.0
    ax.text(mid, max(gp) * 100 + max(gp) * 100 * 0.06, "gappy cols 28-31",
            ha="center", va="bottom", fontsize=9, color=BRICK)
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "SKILL.md: columns with gap fraction >= 0.5 are alignment-artifact candidates -> remove or mask")
f3 = savefig(fig, "fig3_gap_profile.png")

total_fail = f1 + f2 + f3
print("FIGURE QUALITY: TOTAL FAILS = %d" % total_fail)
import sys
sys.exit(1 if total_fail else 0)
