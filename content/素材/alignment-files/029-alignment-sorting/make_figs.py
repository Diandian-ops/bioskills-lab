#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""029 alignment-sorting 出图：基于 029 fig_data.json 真实数组。"""
import os, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(BASE, "fig_data.json")))
BRICK = "#b5482f"; CELADON = "#2f7d72"; GREY = "#6b7280"

def _verify(fig, name):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fails = 0
    for ax in fig.axes:
        axb = ax.get_window_extent(r)
        bbs = []
        for t in ax.texts:
            bb = t.get_window_extent(r)
            if (bb.x0 < axb.x0 - 3 or bb.x1 > axb.x1 + 3 or bb.y0 < axb.y0 - 3 or bb.y1 > axb.y1 + 3):
                print("  [FAIL] %s: text out of axes: %r" % (name, t.get_text()[:50])); fails += 1
            bbs.append(bb)
        for i in range(len(bbs)):
            for j in range(i + 1, len(bbs)):
                a, b = bbs[i], bbs[j]
                if not (a.x1 + 1 < b.x0 or b.x1 + 1 < a.x0 or a.y1 + 1 < b.y0 or b.y1 + 1 < a.y0):
                    print("  [FAIL] %s: text overlap" % name); fails += 1
    if fails == 0: print("  [PASS] %s" % name)
    return fails

def savefig(fig, name):
    fails = _verify(fig, name)
    fig.savefig(os.path.join(BASE, name), dpi=150); plt.close(fig); return fails

def footnote(fig, text):
    fig.text(0.5, 0.012, text, ha="center", va="bottom", fontsize=8.5, color=GREY)

np_, sp_ = D["name_order_pos"], D["coord_order_pos"]
idx = list(range(len(np_)))
fig, ax = plt.subplots(figsize=(8.6, 5.2), dpi=150)
ax.plot(idx, np_, color=BRICK, lw=1.3, marker="o", ms=3, label="name-sorted file order")
ax.plot(idx, sp_, color=CELADON, lw=1.3, marker="s", ms=3, label="coordinate-sorted file order")
ax.set_xlabel("read index (first 60 reads)", fontsize=10.5)
ax.set_ylabel("reference position", fontsize=10.5)
ax.set_title("Sorting reshuffles read order: POS scattered -> monotonic", fontsize=12, pad=12)
ax.legend(loc="upper left", fontsize=9, framealpha=0.9)
ax.grid(alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "samtools sort -n (name) vs samtools sort (coordinate). Input was already SO:coordinate, so counts are preserved 6000 -> 6000")
f1 = savefig(fig, "fig1_sort_order.png")

fig, ax = plt.subplots(figsize=(8.6, 4.6), dpi=150)
ax.bar(["input (name/any order)", "output (coordinate)"], [6000, 6000], color=[GREY, CELADON], width=0.5, zorder=3)
ax.set_ylabel("alignment count", fontsize=10.5)
ax.set_title("Alignment count is preserved by sorting", fontsize=12, pad=12)
ax.set_ylim(0, 7200)
for i, v in enumerate([6000, 6000]):
    ax.text(i, v + 120, str(v), ha="center", va="bottom", fontsize=9.5, color="#222222")
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
f2 = savefig(fig, "fig2_count_preserved.png")

total = f1 + f2
print("FIGURE QUALITY: TOTAL FAILS = %d" % total)
import sys; sys.exit(1 if total else 0)
