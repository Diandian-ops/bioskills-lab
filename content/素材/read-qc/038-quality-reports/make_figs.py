#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
038 quality-reports 出图：全部基于本次 WSL 真跑解析产物 qc_summary.json，无任何虚构数字。
数据来源（同目录真实文件）：
  - qc_summary.json  （parse_qc.py 解析 fastqc_data.txt 与 multiqc_data/ 的真实数值）
  - 真跑样本：S1_good（高质量）/ S2_degraded（3' 质量衰减）/ S3_adapter（adapter 污染）
配色：brick-red #b5482f + celadon #2f7d72。英文标签（DejaVu Sans 无 CJK）。
说明性文字全部并入双行 set_title 或图内空白区，不放轴外。
自检：check_figs.py 通过时打印恰好 FIGURE QUALITY: TOTAL FAILS = 0
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(BASE, "qc_summary.json")

BRICK = "#b5482f"
CELADON = "#2f7d72"
GREY = "#6b7280"
DARK = "#222222"

with open(DATA, encoding="utf-8") as f:
    D = json.load(f)

SAMPLES = ["S1_good", "S2_degraded", "S3_adapter"]
S = D["samples"]
COLORS = {"S1_good": CELADON, "S2_degraded": BRICK, "S3_adapter": GREY}


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
# FIG 1: per-base mean Phred quality, 3 samples (real curves)
# ============================================================
fig, ax = plt.subplots(figsize=(7.6, 4.5))
for s in SAMPLES:
    v = S[s]
    x = v["per_base_quality_positions"]
    y = v["per_base_quality_mean"]
    ax.plot(x, y, color=COLORS[s], lw=2, label=s)
ax.axhline(20, color=BRICK, ls="--", lw=0.8, alpha=0.5)
ax.text(2, 21.2, "Q20 (fail threshold on median)", fontsize=8, color=BRICK)
ax.set_xlabel("read position (bp)")
ax.set_ylabel("mean Phred quality (fastqc_data.txt)")
ax.set_ylim(0, 48)
ax.set_xlim(1, 102)
ax.set_title("Per-base quality: only S2_degraded collapses toward the 3' end\n"
             "S1 %.1f -> %.1f, S2 %.1f -> %.1f, S3 %.1f -> %.1f; "
             "S2 mean Q at last cycle = %.1f -> module FAIL"
             % (S["S1_good"]["per_base_quality_mean"][0], S["S1_good"]["per_base_quality_mean"][-1],
                S["S2_degraded"]["per_base_quality_mean"][0], S["S2_degraded"]["per_base_quality_mean"][-1],
                S["S3_adapter"]["per_base_quality_mean"][0], S["S3_adapter"]["per_base_quality_mean"][-1],
                S["S2_degraded"]["per_base_quality_mean"][-1]),
             fontsize=10.5)
ax.legend(fontsize=9, loc="lower left", framealpha=0.9)
plt.tight_layout()
total_fail = savefig(fig, "fig1_per_base_quality.png")

# ============================================================
# FIG 2: module status matrix, 10 modules x 3 samples (real)
# ============================================================
modules = list(S["S1_good"]["module_statuses"].keys())   # 10 modules, same order
status_colors = {"pass": CELADON, "warn": "#d8a15a", "fail": BRICK}
fig, ax = plt.subplots(figsize=(7.6, 4.9))
for j, s in enumerate(SAMPLES):
    for i, m in enumerate(modules):
        st = S[s]["module_statuses"][m]
        ax.add_patch(Rectangle((j, i), 0.82, 0.82,
                               facecolor=status_colors[st], edgecolor="white"))
        ax.text(j + 0.41, i + 0.41, st[0].upper(), ha="center", va="center",
                fontsize=9, color="white", fontweight="bold")
ax.set_xlim(-0.15, 4.1)
ax.set_ylim(len(modules), -0.6)
ax.set_xticks([j + 0.41 for j in range(3)])
ax.set_xticklabels(SAMPLES, fontsize=10)
ax.set_yticks([i + 0.41 for i in range(len(modules))])
ax.set_yticklabels(modules, fontsize=8.5)
ax.tick_params(length=0)
for sp in ax.spines.values():
    sp.set_visible(False)
ax.set_title("FastQC module verdicts: 10 modules x 3 samples\n"
             "S1 all pass; S2 per-base quality FAIL; S3 overrep WARN + adapter FAIL",
             fontsize=10.5)
handles = [Patch(facecolor=status_colors["pass"], label="pass"),
           Patch(facecolor=status_colors["warn"], label="warn"),
           Patch(facecolor=status_colors["fail"], label="fail")]
ax.legend(handles=handles, fontsize=9, loc="upper right", bbox_to_anchor=(0.99, 1.0),
          framealpha=0.95)
plt.tight_layout()
total_fail += savefig(fig, "fig2_module_status_matrix.png")

# ============================================================
# FIG 3: adapter content curves + overrepresented evidence (real)
# ============================================================
fig, ax = plt.subplots(figsize=(7.6, 4.5))
for s in SAMPLES:
    v = S[s]
    x = v["adapter_positions"]
    ax.plot(x, v["adapter_total_by_position"], color=COLORS[s], lw=2, label=s)
ax.axhline(10, color=BRICK, ls="--", lw=0.9, alpha=0.6)
ax.text(2, 10.9, "10% fail threshold", fontsize=8, color=BRICK)
ax.axhline(5, color=GREY, ls=":", lw=0.9, alpha=0.7)
ax.text(2, 5.9, "5% warn threshold", fontsize=8, color=GREY)
ax.set_xlabel("read position (bp)")
ax.set_ylabel("total adapter content (%)")
ax.set_xlim(1, 102)
ax.set_ylim(0, 26)
ax.set_title("S3 read-through: adapter content climbs to %.1f%% at the 3' end (FAIL)\n"
             "top overrepresented seq %.2f%% (%d reads) = TruSeq Adapter, Index 9"
             % (S["S3_adapter"]["max_adapter_percent"],
                S["S3_adapter"]["overrepresented_max_percent"],
                S["S3_adapter"]["top_overrepresented"]["count"]),
             fontsize=10.5)
ax.legend(fontsize=9, loc="upper left", framealpha=0.95)
plt.tight_layout()
total_fail += savefig(fig, "fig3_adapter_content.png")

print("FIGURE QUALITY: TOTAL FAILS = %d" % total_fail)
