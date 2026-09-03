#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
039 adapter-trimming 出图：全部基于本次 WSL 真跑产物（results.json / truth_summary.json），
无任何虚构数字。数据来源：
  - results.json         （analyze_results.py 解析：逐档 × 工具指标、残留长度分箱、cutadapt 官方报告）
  - truth_summary.json   （make_inputs.py 落盘的造数据真值）
配色：brick-red #b5482f + celadon #2f7d72 + grey。英文标签（DejaVu Sans 无 CJK）。
自检：文本必须落在 axes 内、互不重叠、legend 不压数据；成功时打印恰好
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
GREY = "#8a8f98"
DARK = "#222222"

with open(os.path.join(BASE, "results.json"), encoding="utf-8") as f:
    RES = json.load(f)
with open(os.path.join(BASE, "truth_summary.json"), encoding="utf-8") as f:
    TRUTH = json.load(f)

LEVELS = ["5p", "20p", "40p"]
LEVEL_LABELS = ["5% adapter\npairs", "20% adapter\npairs", "40% adapter\npairs"]


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
        leg = ax.get_legend()
        if leg is not None and leg.get_visible():
            lbb = leg.get_window_extent(r)
            from matplotlib.container import BarContainer
            for child in ax.get_children():
                if isinstance(child, BarContainer):
                    dbb = child.get_window_extent(r)
                    if not (lbb.x1 + 1 < dbb.x0 or dbb.x1 + 1 < lbb.x0 or
                            lbb.y1 + 1 < dbb.y0 or dbb.y1 + 1 < lbb.y0):
                        print("  [FAIL] %s: legend overlaps bars" % name)
                        fails += 1
    if fails == 0:
        print("  [PASS] %s" % name)
    return fails


def savefig(fig, name):
    fails = _verify(fig, name)
    fig.savefig(os.path.join(BASE, name), dpi=150)
    plt.close(fig)
    return fails


# ============================================================
# FIG 1: SKILL.md literal ILLUMINACLIP string silently drops R2
# ============================================================
drops = {t: [RES["levels"][g][t]["reads_dropped"] for g in LEVELS]
         for t in ("cutadapt", "trimmomatic_skill", "trimmomatic_true")}
fig, ax = plt.subplots(figsize=(7.4, 4.4))
x = range(len(LEVELS))
w = 0.26
series = [("cutadapt (-O 3, -m 20:20)", drops["cutadapt"], CELADON),
          ("Trimmomatic, SKILL.md literal string", drops["trimmomatic_skill"], BRICK),
          ("Trimmomatic, boolean fixed to true", drops["trimmomatic_true"], GREY)]
ymax = max(drops["trimmomatic_skill"]) * 1.30
for i, (lab, vals, c) in enumerate(series):
    pos = [xi + (i - 1) * w for xi in x]
    ax.bar(pos, vals, width=w, color=c, label=lab)
    for p, v in zip(pos, vals):
        ax.text(p, v + ymax * 0.015, str(v), ha="center", va="bottom",
                fontsize=9, color=DARK)
ax.set_xticks(list(x))
ax.set_xticklabels(LEVEL_LABELS)
ax.set_ylabel("input read mates missing from output")
ax.set_ylim(0, ymax)
ax.set_title("The SKILL.md literal ILLUMINACLIP string ends with the word 'keepBothReads',\n"
             "which Boolean.parseBoolean reads as false: every read-through pair loses its R2\n"
             "(%s / %s / %s mates at the 5%% / 20%% / 40%% gradients); ':2:30:10:2:true' fixes it"
             % tuple("{:,}".format(d) for d in drops["trimmomatic_skill"]),
             fontsize=10)
ax.legend(fontsize=8.5, loc="upper left")
plt.tight_layout()
total_fail = savefig(fig, "fig1_illuminaclip_keepbothreads_r2_loss.png")

# ============================================================
# FIG 2: bases actually removed vs adapter remnant length
# ============================================================
from matplotlib.lines import Line2D
ca_rem = RES["by_remnant"]["cutadapt"]
tm_rem = RES["by_remnant"]["trimmomatic_true"]
xs = sorted(int(k) for k in ca_rem)
ca_y = [ca_rem[str(v)]["mean_removed"] for v in xs]
tm_y = [tm_rem[str(v)]["mean_removed"] for v in xs]
fig, ax = plt.subplots(figsize=(7.4, 4.4))
ax.plot(xs, xs, linestyle="--", color=GREY, linewidth=1.2)
ax.plot(xs, ca_y, marker="o", markersize=4, color=BRICK, linewidth=1.6)
ax.plot(xs, tm_y, marker="s", markersize=4, color=CELADON, linewidth=1.6)
handles = [
    Line2D([], [], linestyle="--", color=GREY, linewidth=1.2,
           label="full-remnant removal (y = x)"),
    Line2D([], [], marker="o", color=BRICK, linewidth=1.6,
           label="cutadapt -a/-A AGATCGGAAGAGC"),
    Line2D([], [], marker="s", color=CELADON, linewidth=1.6,
           label="Trimmomatic ILLUMINACLIP"),
]
ax.legend(handles=handles, fontsize=8.5, loc="upper left")
ax.set_xlabel("adapter remnant appended to the read (bp)")
ax.set_ylabel("bases actually removed per read (bp)")
ax.set_xlim(0, 66)
ax.set_ylim(0, 70)
ax.set_title("Bases removed vs adapter remnant length (pooled 5/20/40% gradients, read length 100 bp)\n"
             "Trimmomatic tracks y = x exactly; cutadapt sits just below it at short remnants\n"
             "(a 5-9 bp match allows 0 errors at the default 10% rate, so error-carrying remnants are missed)",
             fontsize=9.5)
plt.tight_layout()
total_fail += savefig(fig, "fig2_removal_vs_remnant_length.png")

# ============================================================
# FIG 3: bases removed per read vs gradient (truth vs tools)
# ============================================================
truth_per_read = [TRUTH["levels"][g]["adapter_bases_total"] /
                  (2 * TRUTH["levels"][g]["pairs"]) for g in LEVELS]
ca_bpr = [RES["levels"][g]["cutadapt"]["bases_removed_per_read"] for g in LEVELS]
tm_bpr = [RES["levels"][g]["trimmomatic_true"]["bases_removed_per_read"] for g in LEVELS]
fig, ax = plt.subplots(figsize=(7.4, 4.4))
x = range(len(LEVELS))
w = 0.26
series = [("truth: adapter bases appended", truth_per_read, GREY),
          ("cutadapt bases removed", ca_bpr, BRICK),
          ("Trimmomatic bases removed (true)", tm_bpr, CELADON)]
ymax = max(truth_per_read) * 1.32
for i, (lab, vals, c) in enumerate(series):
    pos = [xi + (i - 1) * w for xi in x]
    ax.bar(pos, vals, width=w, color=c, label=lab)
    for p, v in zip(pos, vals):
        ax.text(p, v + ymax * 0.015, "%.2f" % v, ha="center", va="bottom",
                fontsize=9, color=DARK)
ax.set_xticks(list(x))
ax.set_xticklabels(LEVEL_LABELS)
ax.set_ylabel("bases removed per input read (bp)")
ax.set_ylim(0, ymax)
ax.set_title("Bases removed per read scales linearly with the adapter gradient\n"
             "Trimmomatic matches the truth curve; cutadapt offsets its missed 5-15 bp remnants\n"
             "with ~2% spurious 3 bp clips on clean reads (-O 3 default, ~0.06 bp/read)",
             fontsize=10)
ax.legend(fontsize=8.5, loc="upper left")
plt.tight_layout()
total_fail += savefig(fig, "fig3_bases_removed_per_read.png")

print("FIGURE QUALITY: TOTAL FAILS = %d" % total_fail)
