#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
040 quality-filtering 出图：全部基于本次 WSL 真跑产物（results.json，无虚构数字）。
数据来源（同目录真实文件）：
  - results.json（parse_results.py 从各工具输出 FASTQ / fastp JSON /
    cutadapt 报告 / trimmomatic 日志解析所得）
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

with open(os.path.join(BASE, "results.json"), encoding="utf-8") as f:
    R = json.load(f)

NAMES = ["fastp_filter", "fastp_cutright", "cutadapt", "trimm_sw", "trimm_maxinfo"]
LABELS = {
    "input": "input (no filtering)",
    "fastp_filter": "fastp -q 20 -u 40 -n 5 -l 36\n(filter only)",
    "fastp_cutright": "fastp --cut_right 4:20 -l 36\n(window trim)",
    "cutadapt": "cutadapt -q 20 -m 36\n(BWA running-sum)",
    "trimm_sw": "trimmomatic SLIDINGWINDOW:4:20\n+LEADING/TRAILING:3 MINLEN:36",
    "trimm_maxinfo": "trimmomatic MAXINFO:40:0.5\nMINLEN:36",
}
N_IN = R["input"]["reads"]  # 20000


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


# ============================================================
# FIG 1: reads retention per configuration (real run)
# ============================================================
rows = [("input", N_IN)] + [(n, R[n]["reads"]) for n in NAMES]
fig, ax = plt.subplots(figsize=(8.4, 4.6))
ypos = range(len(rows) - 1, -1, -1)
colors = [CELADON] + [BRICK] * (len(rows) - 1)
ax.barh(list(ypos), [v for _, v in rows], color=colors, height=0.62)
ax.set_yticks(list(ypos))
ax.set_yticklabels([LABELS[n] for n, _ in rows], fontsize=8.5)
for y, (n, v) in zip(ypos, rows):
    pct = 100.0 * v / N_IN
    ax.text(v + N_IN * 0.012, y, "%d (%.2f%%)" % (v, pct),
            va="center", ha="left", fontsize=9, color="#222222")
ax.set_xlim(0, N_IN * 1.33)
ax.set_xlabel("reads kept (of %d input)" % N_IN)
ax.set_title("Reads kept after quality filtering / trimming (same 20,000-read input)\n"
             "window trims keep ~92% of reads but shorten them; "
             "MAXINFO keeps all 20,000 by trimming harder", fontsize=10.5)
plt.tight_layout()
total_fail = savefig(fig, "fig1_reads_kept.png")

# ============================================================
# FIG 2: per-cycle mean quality profile, input vs outputs (real parse)
# ============================================================
fig, ax = plt.subplots(figsize=(8.4, 4.8))
styles = {
    "input": dict(color="#222222", lw=2.2, ls="--"),
    "fastp_filter": dict(color=GREY, lw=1.6, ls=":"),
    "fastp_cutright": dict(color=BRICK, lw=1.7),
    "cutadapt": dict(color=CELADON, lw=1.7),
    "trimm_sw": dict(color="#6fb3a9", lw=1.7),
    "trimm_maxinfo": dict(color="#d98e6b", lw=1.7, ls="-."),
}
for key in ["input", "fastp_filter", "fastp_cutright", "cutadapt", "trimm_sw", "trimm_maxinfo"]:
    prof = R[key]["per_cycle_q"]
    lab = "input" if key == "input" else LABELS[key].replace("\n", " ")
    ax.plot(range(1, len(prof) + 1), prof, label=lab, **styles[key])
ax.axhline(20, color=GREY, lw=0.8, ls=":", zorder=0)
ax.set_xlabel("cycle (5' -> 3')")
ax.set_ylabel("mean Phred quality")
ax.set_ylim(2, 42)
ax.set_xlim(1, 150)
ax.legend(fontsize=7.8, loc="lower left", framealpha=0.92)
ax.set_title("Per-cycle mean quality: input decays Q~37 -> Q~12 (plus ~8% globally low-Q reads);\n"
             "window trims cut the decayed tail, keeping quality above the Q20 line (dotted)",
             fontsize=10)
plt.tight_layout()
total_fail += savefig(fig, "fig2_per_cycle_quality.png")

# ============================================================
# FIG 3: post-trim mean Q and mean length per tool (real parse)
# ============================================================
fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.6, 4.6))
xs = range(len(NAMES))
vals_q = [R[n]["mean_q"] for n in NAMES]
vals_len = [R[n]["mean_len"] for n in NAMES]
axL.bar(xs, vals_q, color=[BRICK, CELADON, CELADON, CELADON, CELADON])
for x, v in zip(xs, vals_q):
    axL.text(x, v + 0.55, "%.2f" % v, ha="center", va="bottom", fontsize=9, color="#222222")
axL.axhline(R["input"]["mean_q"], color="#222222", ls="--", lw=1.2)
axL.text(len(NAMES) - 0.45, R["input"]["mean_q"] + 0.5, "input %.2f" % R["input"]["mean_q"],
         ha="right", va="bottom", fontsize=8.5, color="#222222")
axL.set_xticks(list(xs))
axL.set_xticklabels(["fastp\nfilter", "fastp\ncut_right", "cutadapt\n-q 20",
                     "trimm\nSW:4:20", "trimm\nMAXINFO:40:0.5"], fontsize=8.5)
axL.set_ylabel("mean Phred Q (kept bases)")
axL.set_ylim(0, 36)
axL.set_title("Post-trim mean quality", fontsize=10.5)

axR.bar(xs, vals_len, color=[BRICK, CELADON, CELADON, CELADON, CELADON])
for x, v in zip(xs, vals_len):
    axR.text(x, v + 2.2, "%.1f" % v, ha="center", va="bottom", fontsize=9, color="#222222")
axR.axhline(R["input"]["mean_len"], color="#222222", ls="--", lw=1.2)
axR.text(len(NAMES) - 0.45, R["input"]["mean_len"] - 9.5, "input 150 bp",
         ha="right", va="bottom", fontsize=8.5, color="#222222")
axR.set_xticks(list(xs))
axR.set_xticklabels(["fastp\nfilter", "fastp\ncut_right", "cutadapt\n-q 20",
                     "trimm\nSW:4:20", "trimm\nMAXINFO:40:0.5"], fontsize=8.5)
axR.set_ylabel("mean read length after filtering (bp)")
axR.set_ylim(0, 168)
axR.set_title("Post-trim mean length", fontsize=10.5)
fig.suptitle("Filter drops whole reads (mean Q 25.97, length intact); window trims raise Q to ~30 "
             "at the cost of ~1/3 of bases", fontsize=10.5, y=0.995)
plt.tight_layout(rect=(0, 0, 1, 0.94))
total_fail += savefig(fig, "fig3_q_vs_length.png")

print("FIGURE QUALITY: TOTAL FAILS = %d" % total_fail)
