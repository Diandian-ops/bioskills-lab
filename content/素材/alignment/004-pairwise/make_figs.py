#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
004 pairwise-alignment 出图：基于 run.py 真跑产出的 pairwise_data.json。
配色：brick-red #b5482f + celadon #2f7d72。英文标签（DejaVu Sans 无 CJK）。
出图质量自检（替代 check_figs.py）：文本包围盒落 axes 内、不重叠。
"""
import os
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = json.load(open(os.path.join(BASE, "pairwise_data.json")))

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

# ============================================================
# 图 1：全局 DNA 比对的身份构成（identities / mismatches / gaps）
# ============================================================
dg = DATA["dna_global"]
fig, ax = plt.subplots(figsize=(8.4, 5.0), dpi=150)
rows = [("Identities", dg["identities"], CELADON),
        ("Mismatches", dg["mismatches"], BRICK),
        ("Gaps", dg["gaps"], GREY)]
vals = [v for _, v, _ in rows]
labels = [l for l, _, _ in rows]
colors = [c for _, _, c in rows]
ypos = range(len(rows))[::-1]
bars = ax.barh(list(ypos), vals, color=colors, zorder=3)
ax.set_yticks(list(ypos))
ax.set_yticklabels(labels, fontsize=10)
ax.set_xlabel("residue count", fontsize=10.5)
ax.set_title("Global DNA alignment identity composition\nseq1 vs seq2 (PairwiseAligner, %d columns)" % dg["alignment_length"],
             fontsize=12, pad=12)
mx = max(vals) if max(vals) > 0 else 1
ax.set_xlim(0, mx * 1.25)
for y, v in zip(ypos, vals):
    ax.text(v + mx * 0.02, y, str(v), va="center", ha="left", fontsize=9.5, color="#222222")
ax.text(0.5, 0.92, "Percent identity (counts) = %.1f%%" % dg["pid_counts_pct"],
        transform=ax.transAxes, ha="center", fontsize=10, color=BRICK,
        bbox=dict(boxstyle="round,pad=0.3", fc="#fff3f0", ec=BRICK, lw=0.8))
ax.grid(axis="x", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "SKILL.md: counts() uses aligned non-gap positions (PID2-like). len(alignment) returns the number of sequences, not the alignment length")
f1 = savefig(fig, "fig1_identity_counts.png")

# ============================================================
# 图 2：比对模式对 score 的影响（global / local / semiglobal）
# ============================================================
ms = DATA["mode_scores"]
fig, ax = plt.subplots(figsize=(8.4, 5.0), dpi=150)
modes = ["global", "local", "semiglobal"]
mvals = [ms[m] for m in modes]
mbars = ax.bar(modes, mvals, color=[CELADON, BRICK, GREY], width=0.55, zorder=3)
ax.set_ylabel("alignment score", fontsize=10.5)
ax.set_title("Alignment mode changes the score\nlocalA vs localB (divergent flanks, conserved core)",
             fontsize=12, pad=12)
mx2 = max(mvals) * 1.25
ax.set_ylim(0, mx2)
for m, v in zip(modes, mvals):
    ax.text(m, v + mx2 * 0.02, str(v), ha="center", va="bottom", fontsize=9.5, color="#222222")
ax.annotate("local ignores divergent\nflanks -> highest score",
            xy=(1, ms["local"]), xytext=(1, mx2 * 0.72), textcoords="data",
            ha="center", va="center", fontsize=9, color=BRICK,
            bbox=dict(boxstyle="round,pad=0.3", fc="#fff3f0", ec=BRICK, lw=0.8),
            arrowprops=dict(arrowstyle="->", color=BRICK, lw=1.1))
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "SKILL.md: global forces end-to-end; local finds best region; semiglobal frees end gaps for fragment/overlap")
f2 = savefig(fig, "fig2_mode_scores.png")

total_fail = f1 + f2
print("FIGURE QUALITY: TOTAL FAILS = %d" % total_fail)
import sys
sys.exit(1 if total_fail else 0)
