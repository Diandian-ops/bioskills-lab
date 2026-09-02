#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""031 filtering 出图：基于 031 filtering_data.json 真实计数。"""
import os, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(BASE, "filtering_data.json")))
BRICK = "#b5482f"; CELADON = "#2f7d72"; GREY = "#6b7280"

def _verify(fig, name):
    fig.canvas.draw(); r = fig.canvas.get_renderer(); fails = 0
    for ax in fig.axes:
        axb = ax.get_window_extent(r); bbs = []
        for t in ax.texts:
            bb = t.get_window_extent(r)
            if (bb.x0 < axb.x0 - 3 or bb.x1 > axb.x1 + 3 or bb.y0 < axb.y0 - 3 or bb.y1 > axb.y1 + 3):
                print("  [FAIL] %s: %r" % (name, t.get_text()[:50])); fails += 1
            bbs.append(bb)
        for i in range(len(bbs)):
            for j in range(i + 1, len(bbs)):
                a, b = bbs[i], bbs[j]
                if not (a.x1 + 1 < b.x0 or b.x1 + 1 < a.x0 or a.y1 + 1 < b.y0 or b.y1 + 1 < a.y0):
                    print("  [FAIL] %s: overlap" % name); fails += 1
    if fails == 0: print("  [PASS] %s" % name)
    return fails

def savefig(fig, name):
    fails = _verify(fig, name)
    fig.savefig(os.path.join(BASE, name), dpi=150); plt.close(fig); return fails

def footnote(fig, text):
    fig.text(0.5, 0.012, text, ha="center", va="bottom", fontsize=8.5, color=GREY)

rows = [
    ("mapped (-F 4)", D["mapped"], CELADON),
    ("unmapped (-f 4)", D["unmapped"], BRICK),
    ("proper pair (-f 2)", D["proper_pair"], CELADON),
    ("MAPQ >= 20 (-q 20)", D["mapq_ge20"], CELADON),
    ("mapped & MAPQ>=20", D["mapped_and_q20"], BRICK),
    ("unmapped, mate mapped (-f 4 -F 8)", D["unmapped_not_mate_unmapped"], GREY),
]
labels = [r[0] for r in rows]; vals = [r[1] for r in rows]; cols = [r[2] for r in rows]
ypos = list(range(len(rows)))[::-1]
fig, ax = plt.subplots(figsize=(8.8, 5.2), dpi=150)
ax.barh(ypos, vals, color=cols, zorder=3)
ax.set_yticks(ypos); ax.set_yticklabels(labels, fontsize=9.5)
ax.set_xlabel("alignment count (of 6000 total)", fontsize=10.5)
ax.set_title("samtools view filter flags & MAPQ subsets", fontsize=12, pad=12)
ax.set_xlim(0, max(vals) * 1.18)
for y, v in zip(ypos, vals):
    ax.text(v + max(vals) * 0.012, y, str(v), va="center", ha="left", fontsize=9, color="#222222")
ax.grid(axis="x", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "Flags are bit-tests, not mutually exclusive: -f requires ALL bits, -F excludes ANY; -q is MAPQ threshold")
f1 = savefig(fig, "fig1_filter_subsets.png")

total = f1
print("FIGURE QUALITY: TOTAL FAILS = %d" % total)
import sys; sys.exit(1 if total else 0)
