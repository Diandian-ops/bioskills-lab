#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
048 coverage-analysis 出图：全部基于本次 WSL 真跑产物（results.json，无虚构数字）。
数据来源（同目录真实文件）：
  - results.json（parse_results.py 从 bedtools genomecov/coverage、samtools
    coverage/depth 输出与 design.json 设计期望解析对账所得）
配色：brick-red #b5482f + celadon #2f7d72。英文标签（DejaVu Sans 无 CJK）。
自检：文本包围盒必须落在 axes 内、互不重叠；成功时打印恰好
  FIGURE QUALITY: TOTAL FAILS = 0
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

with open(os.path.join(BASE, "results.json"), encoding="utf-8") as f:
    R = json.load(f)


# ---------- 出图自检 ----------
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


T = R["targets"]
NAMES = [t["name"] for t in T]

# ============================================================
# FIG 1: designed vs measured mean read depth per target (reconciliation)
# ============================================================
fig, ax = plt.subplots(figsize=(8.6, 4.9))
xs = range(len(T))
w = 0.36
exp_v = [t["expected_mean_read_depth"] for t in T]
mea_v = [t["measured_mean_read_depth"] for t in T]
ax.bar([x - w / 2 for x in xs], exp_v, width=w, color=CELADON,
       label="designed mean depth (from sampling plan)")
ax.bar([x + w / 2 for x in xs], mea_v, width=w, color=BRICK,
       label="measured mean depth (bedtools coverage -mean)")
for x, t in zip(xs, T):
    lab = "%.2fx / %.2fx\nbreadth %.2f%%" % (
        t["expected_mean_read_depth"], t["measured_mean_read_depth"], t["breadth_pct"])
    ax.text(x, max(t["expected_mean_read_depth"], t["measured_mean_read_depth"]) + 2.6,
            lab, ha="center", va="bottom", fontsize=8.6, color=DARK)
ax.set_xticks(list(xs))
ax.set_xticklabels(["BG\n100 kb @ ~10x", "HIGH\n50 kb @ ~60x",
                    "LOW\n30 kb @ ~1.5x", "ZERO\n20 kb @ 0x"], fontsize=9)
ax.set_ylabel("mean read depth (x)")
ax.set_ylim(0, 104)
ax.set_xlim(-0.55, len(T) - 0.45)
ax.legend(fontsize=8.5, loc="upper right", framealpha=0.92)
ax.set_title("Designed vs measured coverage over 4 engineered regions (chrS 2 Mb)\n"
             "label = designed / measured mean depth, breadth >=1x; all ratios 1.000, ZERO exactly 0",
             fontsize=10.5)
plt.tight_layout()
total_fail = savefig(fig, "fig1_design_vs_measured.png")

# ============================================================
# FIG 2: 2 kb-binned depth track along chrS with designed regions
# ============================================================
track = R["_track_2kb"]
BIN_BP = 2000
xs_t = [(i + 0.5) * BIN_BP / 1e6 for i in range(len(track))]
fig, ax = plt.subplots(figsize=(9.4, 4.7))
ax.fill_between(xs_t, track, color=CELADON, alpha=0.35, lw=0)
ax.plot(xs_t, track, color=CELADON, lw=1.0)
for (a, b), lab, col in [
    ((100000, 200000), "BG slice", GREY),
    ((500000, 550000), "HIGH", BRICK),
    ((1000000, 1030000), "LOW", GREY),
    ((1500000, 1520000), "ZERO", BRICK),
]:
    ax.axvspan(a / 1e6, b / 1e6, color=col, alpha=0.12, lw=0)
    ax.text((a + b) / 2e6, 97, lab, ha="center", va="top", fontsize=9,
            color=DARK if col == GREY else BRICK, fontweight="bold")
ax.set_xlim(0, 2.0)
ax.set_ylim(0, 106)
ax.set_xlabel("position on chrS (Mb)")
ax.set_ylabel("mean depth per 2 kb bin (x)")
ax.set_title("Coverage track along chrS (bedtools genomecov -bga / samtools depth -a):\n"
             "HIGH spike ~86x, LOW at ~1x with holes, ZERO gap at 1.50-1.52 Mb, background ~6.7x",
             fontsize=10.5)
plt.tight_layout()
total_fail += savefig(fig, "fig2_track.png")

# ============================================================
# FIG 3: depth histogram + cumulative breadth curve (genome-wide)
# ============================================================
hist = R["_hist"]
curve = R["_breadth_curve"]
fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.8, 4.7))

hd = [h[0] for h in hist if h[0] <= 45]
hb = [h[1] for h in hist if h[0] <= 45]
axL.bar(hd, hb, color=BRICK, width=0.85)
axL.set_yscale("log")
axL.set_ylim(1, 2e6)
mode_d = max(hist, key=lambda h: h[1])
axL.text(5, 6e5, "mode %dx" % mode_d[0], ha="center", va="center",
         fontsize=8.6, color=DARK)
axL.text(15, 6e5, "depth 0: %s bases (%.2f%%)" %
         (format(hist[0][1], ","), 100.0 * hist[0][1] / R["genome_dist"]["total_bases"]),
         ha="left", va="center", fontsize=8.6, color=DARK)
axL.set_xlim(-0.8, 45)
axL.set_xlabel("depth (x), capped at 45 of max 131")
axL.set_ylabel("bases at this depth (log scale)")
axL.set_title("Genome-wide depth histogram\n(genomecov default output)", fontsize=10)

cd = [c[0] for c in curve]
cb = [c[1] * 100.0 for c in curve]
axR.plot(cd, cb, color=CELADON, lw=1.8)
axR.axhline(50, color=GREY, lw=0.8, ls=":")
marks = [(1, 98.35, "1x: 98.35%", (3.5, 90.5)),
         (5, 77.79, "5x: 77.79%", (9.0, 76.0)),
         (10, 15.72, "10x: 15.72%", (14.5, 22.0)),
         (20, 2.5, "20x: 2.50%", (28.0, 10.0))]
for x, y, lab, (tx, ty) in marks:
    axR.scatter([x], [y], s=22, color=BRICK, zorder=3)
    axR.text(tx, ty, lab, ha="left", va="center", fontsize=8.6, color=DARK)
axR.text(38, 53.5, "median = 6x", ha="left", va="bottom", fontsize=8.6, color=GREY)
axR.set_xlim(0, 132)
axR.set_ylim(0, 105)
axR.set_xlabel("depth threshold (x)")
axR.set_ylabel("% of bases with depth >= threshold")
axR.set_title("Cumulative breadth curve: mean 8.51x but median 6x,\n"
              "mean/median = 1.419 (skewed by the HIGH tail)", fontsize=10)
fig.suptitle("Breadth, not mean, tells the story (2,000,000 bp chrS)", fontsize=11, y=1.0)
plt.tight_layout(rect=(0, 0, 1, 0.97))
total_fail += savefig(fig, "fig3_hist_breadth.png")

print("FIGURE QUALITY: TOTAL FAILS = %d" % total_fail)
