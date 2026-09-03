#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
041 fastp-workflow 出图：全部基于本次 WSL 真跑产物（4 组 fastp JSON 报告），无虚构数字。
数据来源（同目录真实文件）：
  report_full.json        全流程（SKILL.md 标准 PE 工作流 + --correction）
  report_noadapter.json   消融 A：关闭接头去除
  report_nolenfilter.json 消融 B：关闭长度过滤
  report_nocut.json       消融 C：去掉 --cut_right 滑窗剪裁
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
DARK = "#222222"


def load(tag):
    with open(os.path.join(BASE, "report_%s.json" % tag), "r", encoding="utf-8") as f:
        return json.load(f)


R = {t: load(t) for t in ("full", "noadapter", "nolenfilter", "nocut")}


def m(tag, *path):
    cur = R[tag]
    for k in path:
        cur = cur[k]
    return cur


# ---------- 真实指标（full run） ----------
B_READS = m("full", "summary", "before_filtering", "total_reads")          # 60000
B_Q20 = m("full", "summary", "before_filtering", "q20_rate")
B_Q30 = m("full", "summary", "before_filtering", "q30_rate")
A_READS = m("full", "summary", "after_filtering", "total_reads")           # 50612
A_Q20 = m("full", "summary", "after_filtering", "q20_rate")
A_Q30 = m("full", "summary", "after_filtering", "q30_rate")
A_MEANLEN = m("full", "summary", "after_filtering", "read1_mean_length")
AD_TR = m("full", "adapter_cutting", "adapter_trimmed_reads")              # 12739
AD_TB = m("full", "adapter_cutting", "adapter_trimmed_bases")
LOWS = m("full", "filtering_result", "low_quality_reads")                  # 0
NS = m("full", "filtering_result", "too_many_N_reads")                     # 0
SHORTS = m("full", "filtering_result", "too_short_reads")                  # 9388
CORR_R = m("full", "filtering_result", "corrected_reads")                  # 5891
CORR_B = m("full", "filtering_result", "corrected_bases")
# 消融组
NC_LOWS = m("nocut", "filtering_result", "low_quality_reads")              # 6592
NC_NS = m("nocut", "filtering_result", "too_many_N_reads")                 # 1506
NC_SHORTS = m("nocut", "filtering_result", "too_short_reads")
NC_READS = m("nocut", "summary", "after_filtering", "total_reads")
NA_READS = m("noadapter", "summary", "after_filtering", "total_reads")
NA_Q30 = m("noadapter", "summary", "after_filtering", "q30_rate")
NL_READS = m("nolenfilter", "summary", "after_filtering", "total_reads")
NL_Q30 = m("nolenfilter", "summary", "after_filtering", "q30_rate")
NC_Q30 = m("nocut", "summary", "after_filtering", "q30_rate")
DUP = m("full", "duplication", "rate")


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
# FIG 1: before vs after (full run) - Q20/Q30 rates and reads
# ============================================================
fig, ax = plt.subplots(figsize=(7.4, 4.5))
labels = ["Q20 rate", "Q30 rate", "mean read\nlength (R1, bp)"]
before_vals = [B_Q20 * 100, B_Q30 * 100, 100.0]
after_vals = [A_Q20 * 100, A_Q30 * 100, float(A_MEANLEN)]
xs = [0, 1, 2]
w = 0.36
b1 = ax.bar([x - w / 2 for x in xs], before_vals, width=w, color=CELADON,
            label="before filtering")
b2 = ax.bar([x + w / 2 for x in xs], after_vals, width=w, color=BRICK,
            label="after filtering (full workflow)")
for x, v in zip(xs, before_vals):
    ax.text(x - w / 2, v + 2.5, "%.2f" % v if v < 100 else "100.0",
            ha="center", va="bottom", fontsize=9.5, color=DARK)
for x, v in zip(xs, after_vals):
    ax.text(x + w / 2, v + 2.5, "%.2f" % v, ha="center", va="bottom",
            fontsize=9.5, color=DARK)
ax.set_xticks(xs)
ax.set_xticklabels(labels)
ax.set_ylabel("percent (%) / bp")
ax.set_ylim(0, 125)
ax.set_yticks([0, 25, 50, 75, 100])
ax.set_title("fastp full workflow, %d reads in -> %d reads out (%.1f%% kept)\n"
             "Q30 rate: %.2f%% -> %.2f%% (+%.2f pp); corrected %d reads / %d bases"
             % (B_READS, A_READS, 100.0 * A_READS / B_READS,
                B_Q30 * 100, A_Q30 * 100, (A_Q30 - B_Q30) * 100, CORR_R, CORR_B),
             fontsize=10.5)
ax.legend(fontsize=9, loc="upper right", framealpha=0.95)
plt.tight_layout()
total_fail = savefig(fig, "fig1_before_after.png")

# ============================================================
# FIG 2: what fastp did to the reads (full run), per工序 outcome
# ============================================================
cats = ["too_short_reads\n(dropped)", "low_quality_reads\n(dropped)",
        "too_many_N_reads\n(dropped)", "adapter_trimmed_reads\n(modified, kept)",
        "corrected_reads\n(modified, kept)"]
vals = [SHORTS, LOWS, NS, AD_TR, CORR_R]
colors = [BRICK, BRICK, BRICK, CELADON, CELADON]
fig, ax = plt.subplots(figsize=(7.6, 4.6))
ypos = list(range(len(cats)))[::-1]
ax.barh(ypos, vals, color=colors, height=0.62)
for y, v in zip(ypos, vals):
    ax.text(v + max(vals) * 0.015, y, "{:,}".format(v), va="center",
            ha="left", fontsize=10, color=DARK)
ax.set_yticks(ypos)
ax.set_yticklabels(cats, fontsize=9)
ax.set_xlabel("reads (of 60,000 input reads)")
ax.set_xlim(0, max(vals) * 1.28)
ax.set_title("Full-run outcome by step (per-read counts)\n"
             "brick = dropped, celadon = trimmed/corrected but kept;\n"
             "the two zero rows: --cut_right sliding-window trim removes those tails first",
             fontsize=10)
plt.tight_layout()
total_fail += savefig(fig, "fig2_step_outcomes.png")

# ============================================================
# FIG 3: ablation - contribution of each step (4 real runs)
# ============================================================
run_names = ["full\n(SKILL workflow)", "no adapter\ntrimming",
             "no length\nfiltering", "no --cut_right\ntrim"]
kept = [A_READS, NA_READS, NL_READS, NC_READS]
q30s = [A_Q30 * 100, NA_Q30 * 100, NL_Q30 * 100, NC_Q30 * 100]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.4, 4.5))
bars1 = ax1.bar(run_names, kept, color=[CELADON, BRICK, BRICK, BRICK])
for x, v in enumerate(kept):
    ax1.text(x, v + max(kept) * 0.015, "{:,}".format(v), ha="center",
             va="bottom", fontsize=9.5, color=DARK)
ax1.set_ylabel("reads kept (after filtering)")
ax1.set_ylim(0, max(kept) * 1.25)
ax1.tick_params(axis="x", labelsize=8.5)
ax1.set_title("Reads kept: %d -> %d / %d / %d / %d" % (B_READS, *kept), fontsize=10)
bars2 = ax2.bar(run_names, q30s, color=[CELADON, BRICK, BRICK, BRICK])
for x, v in enumerate(q30s):
    ax2.text(x, v + 1.2, "%.2f" % v, ha="center", va="bottom",
             fontsize=9.5, color=DARK)
ax2.set_ylabel("Q30 rate after filtering (%)")
ax2.set_ylim(0, 78)
ax2.tick_params(axis="x", labelsize=8.5)
ax2.set_title("Q30 rate: dropping adapter trim costs %.2f pp" % (A_Q30 * 100 - NA_Q30 * 100),
              fontsize=10)
fig.suptitle("Ablation on the same input: without --cut_right the per-read filters fire "
             "(%d low-quality + %d N reads dropped)" % (NC_LOWS, NC_NS), fontsize=10.5)
plt.tight_layout(rect=(0, 0, 1, 0.94))
total_fail += savefig(fig, "fig3_ablation_steps.png")

print("duplication rate (diagnostic only): %s" % DUP)
print("FIGURE QUALITY: TOTAL FAILS = %d" % total_fail)
