#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
013 hisat2-alignment 出图：基于真实 flagstat / hisat2 日志的真实数字。
配色：brick-red #b5482f + celadon #2f7d72。英文标签（DejaVu Sans 无 CJK）。
出图质量自检（替代 check_figs.py）：文本包围盒落 axes 内、不重叠。

数据来源（均已在盘、未重新跑 hisat2）：
  - flagstat_pe.txt / flagstat_rf.txt / flagstat_dta.txt / flagstat_junction*.txt
    （本次用 samtools 1.24 flagstat 实跑已有 BAM 得到）
  - pe.log / rf.log / dta.log / j.log / nov.log / neg.log / hisat2_results.json
    （前 agent 真跑 hisat2-build / hisat2 的留存产物）
"""
import os
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))


def read_flagstat(name):
    """返回 (total, mapped, mapped_pct, properly_paired, pp_pct)。"""
    txt = open(os.path.join(BASE, name)).read()
    total = int(re.search(r"(\d+) \+ \d+ in total", txt).group(1))
    m = re.search(r"(\d+) \+ \d+ mapped \(([\d.]+)%", txt)
    mapped, mapped_pct = int(m.group(1)), float(m.group(2))
    pp = re.search(r"(\d+) \+ \d+ properly paired", txt)
    pp_v = int(pp.group(1)) if pp else 0
    pp_pct = round(pp_v / total * 100, 2) if total else 0.0
    return total, mapped, mapped_pct, pp_v, pp_pct


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


# ============================================================
# 图 1：各模式比对率与正确配对率（PE / RF / DTA）
# 真实数（samtools flagstat 1.24 实跑已有 BAM）：
#   PE  : 5950/6000 mapped 99.17% ; properly paired 5906/6000 98.43%
#   RF  : 5950/6000 mapped 99.17% ; properly paired 5906/6000 98.43%
#   DTA : 5883/6000 mapped 98.05% ; properly paired 5782/6000 96.37%
# ============================================================
pe = read_flagstat("flagstat_pe.txt")
rf = read_flagstat("flagstat_rf.txt")
dta = read_flagstat("flagstat_dta.txt")

modes = ["PE", "RF", "DTA"]
map_pct = [pe[2], rf[2], dta[2]]
pp_pct = [pe[4], rf[4], dta[4]]

fig, ax = plt.subplots(figsize=(8.4, 5.0), dpi=150)
import numpy as np
x = np.arange(len(modes))
w = 0.36
b1 = ax.bar(x - w / 2, map_pct, w, label="Mapping rate", color=CELADON, zorder=3)
b2 = ax.bar(x + w / 2, pp_pct, w, label="Properly-paired rate", color=BRICK, zorder=3)
ax.set_xticks(x)
ax.set_xticklabels(modes, fontsize=11)
ax.set_ylabel("percent of 6000 input reads", fontsize=10.5)
ax.set_title("HISAT2 alignment: mapping & properly-paired rate by mode\n(synthetic 3000-pair PE; --dta lowers both)",
             fontsize=12, pad=12)
ax.set_ylim(94, 100.5)
for bars in (b1, b2):
    for bar in bars:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, h + 0.15, "%.2f" % h,
                ha="center", va="bottom", fontsize=9, color="#222222")
ax.legend(loc="lower left", fontsize=9.5, framealpha=0.9)
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "Sources: flagstat_pe.txt / flagstat_rf.txt / flagstat_dta.txt (samtools 1.24, real BAMs on disk)")
f1 = savefig(fig, "fig1_mapping_paired_rate.png")


# ============================================================
# 图 2：--dta 对协调比对读长占比的影响（核心论断之一）
# 真实数（pe.log / dta.log overall alignment rate 段落）：
#   PE : concordantly exactly 1 = 2953 (98.43%) ; 0 times = 47 (1.57%)
#   DTA: concordantly exactly 1 = 2891 (96.37%) ; 0 times = 109 (3.63%)
# --dta 抬高最短锚定阈值，抑制短锚定 junction 读长 -> 未比对比例翻倍。
# ============================================================
pe_ex1 = 2953
pe_0 = 47
dta_ex1 = 2891
dta_0 = 109

fig, ax = plt.subplots(figsize=(8.4, 5.0), dpi=150)
labels = ["PE", "DTA"]
ex1 = [pe_ex1, dta_ex1]
zero = [pe_0, dta_0]
b1 = ax.bar(labels, ex1, 0.5, label="Aligned concordantly exactly 1", color=CELADON, zorder=3)
b2 = ax.bar(labels, zero, 0.5, bottom=ex1, label="Aligned 0 times", color=BRICK, zorder=3)
ax.set_ylabel("read pairs (of 3000)", fontsize=10.5)
ax.set_title("--dta suppresses short-anchor junction reads\nunmapped pairs rise from 47 to 109 (x2.3)",
             fontsize=12, pad=12)
ax.set_ylim(0, 3200)
for i, (e, z) in enumerate(zip(ex1, zero)):
    ax.text(i, e / 2, "%d\n(%.2f%%)" % (e, e / 3000 * 100),
            ha="center", va="center", fontsize=9, color="white")
    ax.text(i, e + z / 2, "%d\n(%.2f%%)" % (z, z / 3000 * 100),
            ha="center", va="center", fontsize=9, color="#222222")
ax.legend(loc="upper right", fontsize=9.5, framealpha=0.9)
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "Sources: pe.log / dta.log (3000 paired reads; --dta lowers junction-read recovery)")
f2 = savefig(fig, "fig2_dta_effect.png")


# ============================================================
# 图 3：剪接 CIGAR 示意 + 两趟法证据（核心论断之二）
# 真实数：
#   junction read: 100 bp spanning intron [1000,1800] -> CIGAR 50M800N50M
#   pass1 discovered 1 novel junction (chr1 999-1800); pass2 CIGAR 50M800N50M (reused)
#   max MAPQ = 60 (pe.log via hisat2_results.json)
# ============================================================
exon_len, intron_len = 50, 800
fig, ax = plt.subplots(figsize=(8.4, 5.0), dpi=150)
ax.set_xlim(0, exon_len + intron_len + exon_len + 40)
ax.set_ylim(0, 10)
ax.axis("off")

# exon1 segment
ax.add_patch(plt.Rectangle((0, 4), exon_len, 2, color=CELADON, zorder=3))
# intron gap (hatched)
ax.add_patch(plt.Rectangle((exon_len, 4), intron_len, 2, fill=False,
                           hatch="////", edgecolor=BRICK, linewidth=1.3, zorder=3))
# exon2 segment
ax.add_patch(plt.Rectangle((exon_len + intron_len, 4), exon_len, 2, color=CELADON, zorder=3))

ax.text(exon_len / 2, 6.3, "50M", ha="center", va="center", fontsize=11, color=CELADON)
ax.text(exon_len + intron_len / 2, 6.3, "800N (intron)", ha="center", va="center",
        fontsize=11, color=BRICK)
ax.text(exon_len + intron_len + exon_len / 2, 6.3, "50M", ha="center", va="center",
        fontsize=11, color=CELADON)
ax.text((exon_len + intron_len + exon_len) / 2, 1.5,
        "HISAT2 spliced CIGAR: 50M800N50M  (one read across an 800 bp intron)",
        ha="center", va="center", fontsize=10, color="#222222")

# evidence annotations (left-aligned text, kept inside axes)
ev = [
    "Novel junctions discovered in pass1: 1  (chr1:999-1800)",
    "Pass2 CIGAR after --novel-splicesite-infile reuse: 50M800N50M (spliced)",
    "Max MAPQ = 60  (GATK-friendly; not STAR's 255)",
]
y0 = 0.90
for i, t in enumerate(ev):
    ax.text(0.02, y0 - i * 0.08, t, transform=ax.transAxes,
            ha="left", va="center", fontsize=9.5, color="#222222")
ax.set_title("Splice-aware alignment: the N in CIGAR proves intron skipping",
             fontsize=12, pad=10)
footnote(fig, "Sources: junction.bam / junction_pass2.bam (samtools view), hisat2_results.json")
f3 = savefig(fig, "fig3_splice_junction.png")


total_fail = f1 + f2 + f3
print("FIGURE QUALITY: TOTAL FAILS = %d" % total_fail)
import sys
sys.exit(1 if total_fail else 0)
