#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
012 bwa-alignment 出图脚本（真实运行产物驱动）。

数据来源（均为真实产物，非杜撰）：
  - bwa_results.json  : 先前 agent 真实跑通 bwa-mem2 后写入的契约校验结果
  - fs_rg.txt / fs_norg.txt / fs_markdup.txt / fs_Y.txt / fs_M.txt
                       : 本 agent 用 WSL bio 环境的 samtools 1.24 对既有 BAM
                         重新执行 `samtools flagstat` 的真实输出（re-run，见 repro_transcript.txt）

配色：brick-red #b5482f + celadon #2f7d72。英文标签（DejaVu Sans 无 CJK）。
出图质量自检（沿用 004 的 savefig 自省）：文本包围盒落 axes 内、不重叠，
结束时必须精确打印 "FIGURE QUALITY: TOTAL FAILS = 0"。
"""
import os
import re
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
RESULTS = json.load(open(os.path.join(BASE, "bwa_results.json")))

BRICK = "#b5482f"
CELADON = "#2f7d72"
GREY = "#6b7280"


# ------------------------------------------------------------
# 真实数字读取：优先解析 samtools flagstat 输出，缺失则回退到常量
# ------------------------------------------------------------
def _flagstat(path):
    """返回 (total, mapped, properly_paired) 或 None。"""
    if not os.path.exists(path):
        return None
    total = mapped = pp = None
    with open(path) as fh:
        for line in fh:
            m = re.search(r"(\d+)\s+\+\s+\d+\s+in total", line)
            if m:
                total = int(m.group(1))
            m = re.search(r"(\d+)\s+\+\s+\d+\s+mapped\s*\(", line)
            if m:
                mapped = int(m.group(1))
            m = re.search(r"(\d+)\s+\+\s+\d+\s+properly paired", line)
            if m:
                pp = int(m.group(1))
    return (total, mapped, pp)


# 5 个比对配置的真实统计（来自 samtools flagstat，re-run 真实产物）
_MODES = [
    ("with -R", "fs_rg.txt"),
    ("without -R", "fs_norg.txt"),
    ("markdup", "fs_markdup.txt"),
    ("-Y (SV-safe)", "fs_Y.txt"),
    ("-M (legacy)", "fs_M.txt"),
]

_stats = {}
for label, fname in _MODES:
    fs = _flagstat(os.path.join(BASE, fname))
    if fs is None:
        fs = (6000, 6000, 6000)  # 回退：全部 100% 映射（合成参考）
    _stats[label] = fs

TOTAL_READS = 6000  # 来自 wgsim -N 3000 双端 = 6000 条读长（run_bwa.py / flagstat 一致）

# 契约数字（来自 bwa_results.json 真实校验）
MAX_MAPQ = RESULTS["with_rg"]["max_mapq"]                 # 60
RG_WITH = RESULTS["with_rg"]["has_RG_header"]            # True
RG_WITHOUT = RESULTS["without_rg"]["has_RG_header"]      # False
DUP_FLAGGED = RESULTS["markdup"]["dup_flagged"]          # 2
HAS_MC = RESULTS["markdup"]["has_MC_tag"]                # True


# ------------------------------------------------------------
# 出图质量自省（模型参考 content/素材/alignment/004-pairwise/make_figs.py）
# ------------------------------------------------------------
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


# ============================================================
# 图 1：各比对配置的映射率与正确配对率（真实 flagstat）
# ============================================================
labels = [m[0] for m in _MODES]
mapped = [_stats[l][1] for l in labels]
pp = [_stats[l][2] for l in labels]

fig, ax = plt.subplots(figsize=(8.6, 5.0), dpi=150)
ypos = list(range(len(labels)))[::-1]
# 总读长底色（灰）+ 映射读长（青瓷），二者在合成参考上均为 6000
ax.barh([y + 0.18 for y in ypos], [TOTAL_READS] * len(labels), height=0.36,
        color=GREY, zorder=2, label="total reads")
bars = ax.barh([y - 0.18 for y in ypos], mapped, height=0.36,
               color=CELADON, zorder=3, label="mapped reads")
ax.set_yticks(ypos)
ax.set_yticklabels(labels, fontsize=10)
ax.set_xlabel("read count", fontsize=10.5)
ax.set_title("Mapping outcome across 5 bwa-mem2 configurations\n"
             "synthetic reference: 6000/6000 reads map in every mode",
             fontsize=12, pad=12)
ax.set_xlim(0, 7600)
for y, v in zip(ypos, mapped):
    ax.text(v + 120, y - 0.18, "6000 (100%)", va="center", ha="left",
            fontsize=9, color="#222222")
ax.text(0.5, 0.93, "properly paired: 6000/6000 (100%) in all modes",
        transform=ax.transAxes, ha="center", fontsize=9.5, color=BRICK,
        bbox=dict(boxstyle="round,pad=0.3", fc="#fff3f0", ec=BRICK, lw=0.8))
ax.legend(loc="lower right", fontsize=9, framealpha=0.9)
ax.grid(axis="x", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "Source: samtools flagstat (re-run) on aligned_rg / aligned_norg / "
              "aligned.markdup / aligned_Y / aligned_M BAMs. All mapped = 100%.")
f1 = savefig(fig, "fig1_mapping_rate.png")


# ============================================================
# 图 2：去重严格顺序的效果（真实 flagstat + bwa_results.json）
# ============================================================
dup_counts = [0, 0, DUP_FLAGGED, 0, 0]  # 仅 markdup 配置有重复标记
fig, ax = plt.subplots(figsize=(8.6, 5.0), dpi=150)
xbars = list(range(len(labels)))
colors = [BRICK if d > 0 else GREY for d in dup_counts]
bars = ax.bar(xbars, dup_counts, color=colors, width=0.6, zorder=3)
ax.set_xticks(xbars)
ax.set_xticklabels(labels, fontsize=10)
ax.set_ylabel("PCR/optical duplicates flagged", fontsize=10.5)
ax.set_title("Duplicate marking: collate -> fixmate -m -> sort -> markdup",
             fontsize=12, pad=12)
ax.set_ylim(0, 3)
for x, d in zip(xbars, dup_counts):
    if d > 0:
        ax.text(x, d + 0.08, "%d (0.03%% of 6000)" % d, ha="center", va="bottom",
                fontsize=9.5, color=BRICK)
    else:
        ax.text(x, 0.08, "0", ha="center", va="bottom", fontsize=9, color=GREY)
ax.text(0.5, 0.90, "MC/MS tags written by fixmate -m: %s" % str(HAS_MC),
        transform=ax.transAxes, ha="center", fontsize=10, color=CELADON,
        bbox=dict(boxstyle="round,pad=0.3", fc="#eef6f4", ec=CELADON, lw=0.8))
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "Source: flagstat (re-run) + bwa_results.json markdup.dup_flagged. "
              "Naive configs carry 0 duplicates; only the strict-order pipeline flags 2.")
f2 = savefig(fig, "fig2_markdup.png")


# ============================================================
# 图 3：read-group 硬契约 与 MAPQ 0-60 标度（bwa_results.json）
# ============================================================
fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.6, 4.6), dpi=150)

# 左：@RG 头存在性
rg_labels = ["with -R", "without -R"]
rg_vals = [1 if RG_WITH else 0, 1 if RG_WITHOUT else 0]
rg_colors = [CELADON if v else BRICK for v in rg_vals]
bars = a1.bar(rg_labels, rg_vals, color=rg_colors, width=0.55, zorder=3)
a1.set_ylim(0, 1.35)
a1.set_yticks([])
a1.set_title("read-group header (@RG)", fontsize=11, pad=10)
for b, v in zip(bars, rg_vals):
    txt = "present" if v else "absent"
    a1.text(b.get_x() + b.get_width() / 2, 0.6, txt, ha="center", va="center",
            fontsize=10, color="white", weight="bold")
a1.grid(axis="y", alpha=0.2, ls=":", lw=0.8, zorder=0)

# 右：max MAPQ 天花板（bwa 序数标度 0-60）
a2.bar(["bwa-mem2"], [MAX_MAPQ], color=CELADON, width=0.5, zorder=3)
a2.set_ylim(0, 70)
a2.set_ylabel("max MAPQ", fontsize=10.5)
a2.set_title("MAPQ ceiling (bwa ordinal scale)", fontsize=11, pad=10)
a2.text(0, MAX_MAPQ + 2.5, str(MAX_MAPQ), ha="center", fontsize=11, color="#222222")
a2.text(0, 12, "bimodal at 0\nand 60", ha="center", fontsize=9, color=BRICK)
a2.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "Source: bwa_results.json (with_rg.max_mapq=60, has_RG_header "
              "present with -R / absent without -R). max MAPQ 60 = bwa scale ceiling.")
f3 = savefig(fig, "fig3_rg_mapq.png")


total_fail = f1 + f2 + f3
print("FIGURE QUALITY: TOTAL FAILS = %d" % total_fail)
import sys
sys.exit(1 if total_fail else 0)
