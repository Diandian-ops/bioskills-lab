#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
014 star-alignment 出图：基于真实 STAR 运行产物（real_trial_stats.json）。
配色：brick-red #b5482f + celadon #2f7d72。英文标签（DejaVu Sans 无 CJK）。
出图质量自检（替代 check_figs.py）：文本包围盒落 axes 内、不重叠。
"""
import os, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(BASE, "real_trial_stats.json")))

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


total = 0

# ============================================================
# Fig 1: MAPQ of UNIQUE reads — 255 (default) vs 60 (--outSAMmapqUnique 60)
# proves the GATK-breaking default and its fix (SKILL.md claim 2)
# ============================================================
mh = D["mapq_hist"]
fig, ax = plt.subplots(figsize=(8.4, 5.0), dpi=150)
x = [0, 1]
vals = [mh["default"]["unique"], mh["mapq60"]["unique"]]
bars = ax.bar(x, vals, color=[BRICK, CELADON], width=0.55, zorder=3)
ax.set_xticks(x)
ax.set_xticklabels(["default", "mapq60"], fontsize=10.5)
ax.set_ylabel("MAPQ assigned to unique reads", fontsize=10.5)
ax.set_title("STAR MAPQ for unique reads:\n255 (default) breaks GATK  ->  60 (--outSAMmapqUnique 60)",
             fontsize=12, pad=12)
ax.set_ylim(0, max(vals) * 1.25)
for xi, v in zip(x, vals):
    ax.text(xi, v + max(vals) * 0.02, str(v), ha="center", va="bottom",
            fontsize=11, color="#222222", fontweight="bold")
ax.text(0.5, 0.90, "multimappers: MAPQ %d (%d reads) in BOTH configs"
        % (mh["default"]["multimapper"], mh["default"]["multimapper_count"]),
        transform=ax.transAxes, ha="center", fontsize=9.5, color=GREY,
        bbox=dict(boxstyle="round,pad=0.3", fc="#f4f5f7", ec=GREY, lw=0.7))
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "SKILL.md: STAR assigns MAPQ 255 (= 'unavailable') to unique reads; GATK drops it. Fix with --outSAMmapqUnique 60.")
f1 = savefig(fig, "fig1_mapq_unique_255_vs_60.png")
total += f1

# ============================================================
# Fig 2: Uniquely mapped reads % across alignment configs
# (real numbers from Log.final.out; configs are reproducible)
# ============================================================
cfg = D["configs"]
order = ["default", "mapq60", "final"]
labels = ["default", "mapq60", "final"]
uniq = [cfg[k]["uniq_pct"] for k in order]
fig, ax = plt.subplots(figsize=(8.4, 5.0), dpi=150)
bars = ax.bar(labels, uniq, color=CELADON, width=0.55, zorder=3)
ax.set_ylabel("uniquely mapped reads %", fontsize=10.5)
ax.set_ylim(0, 110)
ax.set_title("Uniquely mapped reads %% across alignment configs\n(STAR %s, %d input read pairs)"
             % (D["star_version"], cfg["default"]["input"]), fontsize=12, pad=12)
for k, b, v in zip(order, bars, uniq):
    ax.text(b.get_x() + b.get_width() / 2, v + 2, "%.2f%%" % v,
            ha="center", va="bottom", fontsize=10, color="#222222")
    ax.text(b.get_x() + b.get_width() / 2, 8,
            "avg len %.1f bp\nmis %.2f%%" % (cfg[k]["avg_len"], cfg[k]["mismatch_pct"]),
            ha="center", va="bottom", fontsize=8.5, color=GREY)
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "SKILL.md: ~0.13%% of read pairs are multimappers (MAPQ 3); they are excluded from 'uniquely mapped'.")
f2 = savefig(fig, "fig2_uniquely_mapped_pct.png")
total += f2

# ============================================================
# Fig 3: GeneCounts strandedness inference (SKILL.md claim 3)
# ReadsPerGene.out.tab columns 2/3/4 (unstranded / forward / reverse)
# ============================================================
gc = D["genecounts"]
fig, ax = plt.subplots(figsize=(8.4, 5.0), dpi=150)
keys = ["unstranded", "forward", "reverse"]
vals = [gc["unstranded"], gc["forward"], gc["reverse"]]
cols = [GREY, BRICK, CELADON]
bars = ax.bar(keys, vals, color=cols, width=0.55, zorder=3)
ax.set_ylabel("reads assigned to genes", fontsize=10.5)
ax.set_title("STAR GeneCounts strandedness inference\n(ReadsPerGene.out.tab cols 2/3/4)",
             fontsize=12, pad=12)
mx = max(vals) * 1.22
ax.set_ylim(0, mx)
for k, v in zip(keys, vals):
    ax.text(k, v + mx * 0.02, str(v), ha="center", va="bottom", fontsize=10.5, color="#222222")
ax.text(0.5, 0.90, "forward ~ reverse  ->  UNSTRANDED library; use column 2",
        transform=ax.transAxes, ha="center", fontsize=9.5, color=BRICK,
        bbox=dict(boxstyle="round,pad=0.3", fc="#fff3f0", ec=BRICK, lw=0.8))
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "SKILL.md: sum col3 vs col4; near-equal means unstranded. Wrong column roughly halves counts.")
f3 = savefig(fig, "fig3_genecounts_strandedness.png")
total += f3

print("FIGURE QUALITY: TOTAL FAILS = %d" % total)
import sys
sys.exit(1 if total else 0)
