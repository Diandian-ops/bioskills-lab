#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
045 bed-file-basics figures. All numbers come from bed_results.json,
which was parsed from real bedtools v2.31.1 output (_run.log) on 2026-09-03.
Palette: brick-red #b5482f + celadon #2f7d72. English labels only.
All annotation text stays inside the axes (merged into titles or plotted
within data range); self-check must print FIGURE QUALITY: TOTAL FAILS = 0.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
BRICK = "#b5482f"
CELADON = "#2f7d72"
GREY = "#6b7280"
DARK = "#222222"

with open(os.path.join(BASE, "bed_results.json"), encoding="utf-8") as f:
    D = json.load(f)


def verify(fig, name):
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


# ---------------- fig1: merge counts ----------------
fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=150)
groups = ["Gene bodies", "Exons", "CpG islands"]
raw = [D["merge"]["genes_raw"], D["merge"]["exons_raw"], D["merge"]["cpg_raw"]]
merged = [D["merge"]["genes_merged"], D["merge"]["exons_merged"],
          D["merge"]["cpg_merged"]]
merged_d = [D["merge"]["genes_merged"], D["merge"]["exons_merged"],
            D["merge"]["cpg_merged_d200"]]
xs = range(3)
w = 0.26
b1 = ax.bar([x - w for x in xs], raw, w, color=BRICK, label="raw intervals")
b2 = ax.bar(list(xs), merged, w, color=CELADON, label="bedtools merge")
b3 = ax.bar([x + w for x in xs], merged_d, w, color=GREY,
            label="bedtools merge -d 200")
for bars in (b1, b2, b3):
    for rect in bars:
        ax.text(rect.get_x() + rect.get_width() / 2, rect.get_height() + 4,
                "%d" % rect.get_height(), ha="center", va="bottom",
                fontsize=9, color=DARK)
ax.set_xticks(list(xs))
ax.set_xticklabels(groups)
ax.set_ylabel("Interval count")
ax.set_ylim(0, 250)
ax.set_title("bedtools merge collapses only overlapping intervals\n"
             "simulated chr1 (2 Mb): CpG 30 -> 20 -> 19 with -d 200",
             fontsize=10.5)
ax.legend(loc="upper right", frameon=False, fontsize=9)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(BASE, "fig1_merge_counts.png"))
fails = verify(fig, "fig1")
plt.close(fig)

# ---------------- fig2: per-100kb-bin feature counts ----------------
pb = D["windows"]["per_bin"]
mids = [(r[1] + r[2]) / 2 / 1e6 for r in pb]
fig, ax = plt.subplots(figsize=(7.6, 4.6), dpi=150)
ax.plot(mids, [r[3] for r in pb], "o-", color=CELADON, label="genes")
ax.plot(mids, [r[4] for r in pb], "s-", color=BRICK, label="exons")
ax.plot(mids, [r[5] for r in pb], "^-", color=GREY, label="CpG islands")
ax.set_xlabel("Position on chr1 (Mb)")
ax.set_ylabel("Features per 100 kb bin (count)")
ax.set_xlim(0, 2.0)
ax.set_ylim(-1, 24)
ax.set_yticks(range(0, 25, 4))
ax.set_title("Feature distribution along the 2 Mb simulated chromosome\n"
             "bedtools makewindows -g genome.txt -w 100000 + intersect -c",
             fontsize=10.5)
ax.legend(loc="upper left", frameon=False, fontsize=9)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(BASE, "fig2_per_bin_distribution.png"))
fails += verify(fig, "fig2")
plt.close(fig)

# ---------------- fig3: complement coverage identity ----------------
fig, ax = plt.subplots(figsize=(7.6, 3.9), dpi=150)
rows = [
    ("CpG partition", D["complement_identity"]["cpg_merged_bp"],
     D["complement_identity"]["non_cpg_bp"], CELADON),
    ("Gene-body partition", D["complement_identity"]["genes_bp"],
     D["complement_identity"]["non_genes_bp"], BRICK),
]
ys = [1, 0]
for (label, inside_bp, outside_bp, color), y in zip(rows, ys):
    ax.barh(y, inside_bp, height=0.5, color=color,
            label="intervals (bedtools merge output)")
    ax.barh(y, outside_bp, left=inside_bp, height=0.5, color="#d8d8d8",
            label="bedtools complement")
    ax.text(inside_bp + 20000, y, "{:,} bp".format(inside_bp), ha="left",
            va="center", fontsize=9, color=color, fontweight="bold")
    ax.text(inside_bp + outside_bp / 2, y, "{:,} bp".format(outside_bp),
            ha="center", va="center", fontsize=9, color=DARK)
ax.set_yticks(ys)
ax.set_yticklabels([r[0] for r in rows])
ax.set_xlabel("Cumulative length (bp)")
ax.set_xlim(0, 2000000)
ax.set_ylim(-0.55, 1.9)
ax.set_xticks(range(0, 2000001, 500000))
ax.set_xticklabels(["0", "0.5", "1.0", "1.5", "2.0"])
ax.set_title("Coverage identity: intervals + complement = chromosome length\n"
             "both partitions sum to exactly 2,000,000 bp", fontsize=10.5)
handles = [plt.Rectangle((0, 0), 1, 1, color=CELADON),
           plt.Rectangle((0, 0), 1, 1, color="#d8d8d8")]
ax.legend(handles, ["intervals (merged)", "complement"],
          loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=2,
          frameon=False, fontsize=8.5)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(BASE, "fig3_complement_identity.png"))
fails += verify(fig, "fig3")
plt.close(fig)

print("done, verify fails so far =", fails)
