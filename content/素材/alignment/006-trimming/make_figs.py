#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
006 alignment-trimming 出图：基于 run_tools.py 真跑产出的 trim_data.json。
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
DATA = json.load(open(os.path.join(BASE, "trim_data.json")))

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

IN = DATA["input"]["in_cols"]
gap_per_col = DATA["input"]["gap_per_col"]

# ===== 图 1：比对长度 before / after =====
tools = DATA["tools"]
names = ["input"] + [t["name"] for t in tools]
vals = [IN] + [t["out_cols"] for t in tools]
cols = [GREY] + [BRICK if "clipkit" in t["name"] else CELADON for t in tools]
fig, ax = plt.subplots(figsize=(9.0, 5.2), dpi=150)
bars = ax.bar(names, vals, color=cols, zorder=3)
ax.set_ylabel("alignment length (columns)", fontsize=10.5)
ax.set_title("Alignment length before vs after trimming\ninput MSA = 348 columns (6 human S1 serine proteases)", fontsize=12, pad=12)
mx = max(vals) * 1.15
ax.set_ylim(0, mx)
ax.set_xticks(range(len(names)))
for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width() / 2, v + mx * 0.012, str(v),
            ha="center", va="bottom", fontsize=8.5, color="#222222")
ax.set_xticklabels(names, rotation=20, ha="right", fontsize=8.5)
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "Grey = input MSA. Brick = ClipKIT, celadon = trimAl / BMGE. Real outputs from SKILL.md commands.")
f1 = savefig(fig, "fig1_alignment_length.png")

# ===== 图 2：每工具移除列数 =====
rm_names = [t["name"] for t in tools]
rm_vals = [t["removed"] for t in tools]
rm_pct = [round(100 * (1 - t["retained_frac"]), 1) for t in tools]
fig, ax = plt.subplots(figsize=(9.0, 5.2), dpi=150)
bcols = [BRICK if "clipkit" in n else CELADON for n in rm_names]
bars = ax.bar(rm_names, rm_vals, color=bcols, zorder=3)
ax.axhline(IN * 0.4, color=GREY, ls="--", lw=1.0, zorder=2)
ax.text(len(rm_names) * 0.5, IN * 0.4 + IN * 0.02, "40% rule (SKILL.md): >40% removal = too aggressive",
        fontsize=8.5, color=GREY, ha="center")
ax.set_ylabel("columns removed", fontsize=10.5)
ax.set_title("Columns removed per trimmer\ninput = 348 columns", fontsize=12, pad=12)
mx = max(rm_vals) * 1.18
ax.set_ylim(0, mx)
ax.set_xticks(range(len(rm_names)))
for b, v, p in zip(bars, rm_vals, rm_pct):
    ax.text(b.get_x() + b.get_width() / 2, v + mx * 0.012, "%d\n(%.1f%%)" % (v, p),
            ha="center", va="bottom", fontsize=8.5, color="#222222")
ax.set_xticklabels(rm_names, rotation=20, ha="right", fontsize=8.5)
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "ClipKIT kpic-smart-gap removed 64.1% here, exceeding the 40% guidance -> dataset too divergent for that mode.")
f2 = savefig(fig, "fig2_removed_per_tool.png")

# ===== 图 3：输入逐列空位 + 各工具移除列标注 =====
x = list(range(IN))
gp = [g * 100 for g in gap_per_col]
fig, ax = plt.subplots(figsize=(9.0, 4.8), dpi=150)
ax.plot(x, gp, color=GREY, lw=1.2, zorder=1, label="input per-column gap %")
rbt = DATA["removed_by_tool"].get("clipkit_kpic", [])
rtr = DATA["removed_by_tool"].get("trimal_auto", [])
if rbt:
    ax.scatter(rbt, [gp[c] for c in rbt], color=BRICK, s=14, zorder=3, label="removed by ClipKIT kpic-smart-gap")
if rtr:
    ax.scatter(rtr, [gp[c] for c in rtr], color=CELADON, s=10, marker="x", zorder=3, label="removed by trimAl -automated1")
ax.set_xlabel("input column index (0-based)", fontsize=10.5)
ax.set_ylabel("gap fraction (%)", fontsize=10.5)
ax.set_title("Input gap profile with trimmed columns flagged\nred = ClipKIT kpic-smart-gap, blue x = trimAl -automated1", fontsize=11.5, pad=10)
ax.set_ylim(0, max(max(gp) * 1.2, 5))
ax.legend(fontsize=8, loc="upper right")
ax.grid(alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "Both trimmers concentrate removals on high-gap columns, matching the skill's gap-fraction logic.")
f3 = savefig(fig, "fig3_gap_distribution.png")

total_fail = f1 + f2 + f3
print("FIGURE QUALITY: TOTAL FAILS = %d" % total_fail)
import sys
sys.exit(1 if total_fail else 0)
