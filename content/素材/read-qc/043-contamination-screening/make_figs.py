#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
043 contamination-screening: 3 figures from the real WSL run (results.json).
All numbers come from the actual kraken2 / FastQ Screen outputs on designed
mixtures of three REAL genomes (E. coli K-12 MG1655, phiX174, lambda).
Palette: brick-red #b5482f + celadon #2f7d72. English labels only.
Self-check: text bboxes must stay inside axes and not overlap; on success the
script prints exactly:  FIGURE QUALITY: TOTAL FAILS = 0
"""
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

BASE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(BASE, "results.json")))

BRICK = "#b5482f"
CELADON = "#2f7d72"
GREY = "#6b7280"
DARK = "#222222"

NAME2SHORT = {
    "Escherichia coli": "E. coli",
    "Sinsheimervirus phiX174": "PhiX",
    "Lambdavirus lambda": "Lambda",
}
SAMPLES = ["S1", "S2", "S3"]
TAXA = ["Ecoli", "PhiX", "Lambda"]


def kdet(s, taxon):
    """kraken2 species-level detected % (confidence 0.1 run)."""
    for name, pct in R["kraken"][s]["species"].items():
        if NAME2SHORT[name] == taxon:
            return pct
    return 0.0


def screen_ohog(s, genome):
    """FastQ Screen one-hit-one-genome % for one panel genome."""
    return R["screen"][s][genome]["%One_hit_one_genome"]


# ---------- figure self-check ----------
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
# FIG 1: designed contamination vs kraken2 detection (conf 0.1)
# ============================================================
groups, designed_v, detected_v = [], [], []
for s in ["S1", "S2"]:
    for taxon in TAXA:
        groups.append("%s\n%s" % (s, taxon))
        designed_v.append(R["design"][s][taxon])
        detected_v.append(kdet(s, taxon))

x = range(len(groups))
w = 0.38
fig, ax = plt.subplots(figsize=(7.6, 4.6))
ax.bar([i - w / 2 for i in x], designed_v, w, color=GREY, label="designed fraction")
ax.bar([i + w / 2 for i in x], detected_v, w, color=BRICK,
       label="kraken2 detected (conf 0.1)")
for i in x:
    ax.text(i, max(designed_v[i], detected_v[i]) + 2.2,
            "%.1f / %.2f" % (designed_v[i], detected_v[i]),
            ha="center", va="bottom", fontsize=8.5, color=DARK)
ax.set_xticks(list(x))
ax.set_xticklabels(groups, fontsize=9)
ax.set_ylabel("percent of read pairs (%)")
ax.set_ylim(0, 108)
ax.set_xlim(-0.6, len(groups) - 0.4)
ax.set_title("Designed contamination fractions are recovered at species level\n"
             "kraken2 mini real-genome DB, conf 0.1, min-hit-groups 2; "
             "labels = designed % / detected %", fontsize=10.5)
ax.legend(fontsize=9, loc="upper right")
plt.tight_layout()
total_fail = savefig(fig, "fig1_designed_vs_detected.png")

# ============================================================
# FIG 2: FastQ Screen one-hit-one-genome per panel genome
# ============================================================
ys, vals, cols, yticks, ylabels = [], [], [], [], []
y = 0
for g in ["Lambda", "PhiX", "Ecoli"]:          # bottom-to-top blocks
    block_top = y + 2
    for s in SAMPLES:
        ys.append(y)
        vals.append(screen_ohog(s, g))
        cols.append({"S1": BRICK, "S2": CELADON, "S3": GREY}[s])
        y += 1
    yticks.append(block_top - 1)
    ylabels.append(g)
    y += 1  # gap between blocks

fig, ax = plt.subplots(figsize=(7.6, 4.8))
ax.barh(ys, vals, height=0.72,
        color=cols)
for yi, v in zip(ys, vals):
    ax.text(v + 1.5, yi, "%.2f" % v, va="center", ha="left",
            fontsize=8.5, color=DARK)
ax.set_yticks(yticks)
ax.set_yticklabels(ylabels)
ax.set_xlabel("reads mapped to this genome, One hit one genome (%)")
ax.set_xlim(0, 112)
ax.set_title("FastQ Screen assigns contaminant reads to their true genome\n"
             "% one_hit_one_genome per panel genome; S1 = 5% PhiX + 5% Lambda, "
             "S2 = 15% + 10%, S3 = clean E. coli", fontsize=10.5)
ax.legend(handles=[Patch(color=BRICK, label="S1"),
                   Patch(color=CELADON, label="S2"),
                   Patch(color=GREY, label="S3 (clean)")],
          fontsize=9, loc="lower right")
plt.tight_layout()
total_fail += savefig(fig, "fig2_fastq_screen_categories.png")

# ============================================================
# FIG 3: confidence sweep - Lambda trimmed, PhiX stable
# ============================================================
confs = ["0.0", "0.05", "0.1", "0.2"]
xc = [0.0, 0.05, 0.1, 0.2]
s1_lam = [R["conf_sweep"]["S1"][c]["species"]["Lambdavirus lambda"] for c in confs]
s2_lam = [R["conf_sweep"]["S2"][c]["species"]["Lambdavirus lambda"] for c in confs]

fig, ax = plt.subplots(figsize=(7.6, 4.6))
ax.axhline(5.0, color=GREY, linestyle="--", linewidth=1)
ax.axhline(10.0, color=GREY, linestyle="--", linewidth=1)
ax.plot(xc, s1_lam, "-o", color=CELADON, markersize=5)
ax.plot(xc, s2_lam, "-o", color=BRICK, markersize=5)
for xi, v in zip(xc, s1_lam):
    ax.text(xi, v + 0.28, "%.2f" % v, ha="center", va="bottom",
            fontsize=8.5, color=CELADON)
for xi, v in zip(xc, s2_lam):
    ax.text(xi, v + 0.28, "%.2f" % v, ha="center", va="bottom",
            fontsize=8.5, color=BRICK)
ax.set_xlabel("kraken2 --confidence")
ax.set_xticks(xc)
ax.set_ylabel("Lambda detected (%)")
ax.set_ylim(4.0, 12.5)
ax.set_xlim(-0.03, 0.23)
ax.set_title("Higher confidence slightly trims the weaker Lambda signal\n"
             "S1 designed 5% Lambda, S2 10%; PhiX unaffected (5.0 / 15.0 at "
             "every confidence)", fontsize=10.5)
ax.legend(handles=[Line2D([], [], color=CELADON, marker="o", label="S1 Lambda"),
                   Line2D([], [], color=BRICK, marker="o", label="S2 Lambda"),
                   Line2D([], [], color=GREY, linestyle="--", linewidth=1,
                          label="designed level (5% / 10%)")],
          fontsize=9, loc="upper right")
plt.tight_layout()
total_fail += savefig(fig, "fig3_confidence_sweep.png")

print("FIGURE QUALITY: TOTAL FAILS = %d" % total_fail)
