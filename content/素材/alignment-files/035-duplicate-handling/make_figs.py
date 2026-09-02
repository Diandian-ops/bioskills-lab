#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""035 duplicate-handling 出图：基于 035 duplicate_data.json 真实 markdup 结果。"""
import os, json, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(BASE, "duplicate_data.json")))
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

def num(s):
    m = re.search(r"(\d+) \+ \d+ (?:primary )?duplicates", s or "")
    return int(m.group(1)) if m else 0
before = num(D.get("before_duplicate"))
after = num(D.get("after_duplicate"))

fig, ax = plt.subplots(figsize=(8.6, 4.8), dpi=150)
ax.bar(["before markdup", "after markdup"], [before, after], color=[GREY, BRICK], width=0.5, zorder=3)
ax.set_ylabel("duplicate alignments flagged", fontsize=10.5)
ax.set_title("samtools markdup: duplicates 0 -> 10", fontsize=12, pad=12)
ax.set_ylim(0, max(before, after) * 1.5 + 2)
for i, v in enumerate([before, after]):
    ax.text(i, v + 0.3, str(v), ha="center", va="bottom", fontsize=10, color="#222222")
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "Pipeline: sort -n -> fixmate -m -> sort -> markdup. rmdup is the legacy single-step alternative")
f1 = savefig(fig, "fig1_markdup_counts.png")

total = f1
print("FIGURE QUALITY: TOTAL FAILS = %d" % total)
import sys; sys.exit(1 if total else 0)
