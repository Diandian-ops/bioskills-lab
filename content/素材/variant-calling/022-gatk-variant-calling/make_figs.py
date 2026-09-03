#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
022 gatk-variant-calling 出图：全部基于本次 WSL 真跑产物，无任何虚构数字。
数据来源（同目录真实文件）：
  - raw.vcf.gz / raw.g.vcf.gz     （HaplotypeCaller 标准模式 / -ERC GVCF 真跑产出）
  - _diag.log                     （诊断：修正法 mpileup 列统计、bcftools call -mv 对照、
                                   samtools stats error rate、wgsim 参数核查）
配色：brick-red #b5482f + celadon #2f7d72。英文标签（DejaVu Sans 无 CJK）。
自检：文本包围盒必须落在 axes 内、互不重叠；成功时打印恰好
  FIGURE QUALITY: TOTAL FAILS = 0
"""
import os
import gzip
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
GVCF = os.path.join(BASE, "raw.g.vcf.gz")

BRICK = "#b5482f"
CELADON = "#2f7d72"
GREY = "#6b7280"


# ---------- 真实数据解析 ----------
def parse_gvcf(path):
    """解析 gVCF：区块数、区块长度（END-POS+1）、DP>0 / DP=0 区块数。"""
    lens, dp0, dpp = [], 0, 0
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            pos, end = int(p[1]), int(p[1])
            for kv in p[7].split(";"):
                if kv.startswith("END="):
                    end = int(kv[4:])
            dp = 0
            fmt = p[8].split(":")
            smp = p[9].split(":")
            for k, v in zip(fmt, smp):
                if k == "DP" and v != ".":
                    dp = int(v)
            lens.append(end - pos + 1)
            if dp == 0:
                dp0 += 1
            else:
                dpp += 1
    return lens, dpp, dp0


lens, N_DPP, N_DP0 = parse_gvcf(GVCF)
N_BLOCKS = len(lens)                       # 2946
N_LEN1 = sum(1 for L in lens if L == 1)    # 单碱基区块
import statistics
MED = statistics.median(lens)              # 1
MEAN = statistics.mean(lens)               # ~4.07
MAX = max(lens)                            # 246
lens_sorted = sorted(lens)
P90 = lens_sorted[int(0.90 * len(lens_sorted))]   # 6
P99 = lens_sorted[int(0.99 * len(lens_sorted))]   # 61

# 诊断实测（_diag.log / samtools stats，见 repro_transcript.txt）
N_NONREF_COLS = 4659     # 修正统计法（同时排除 . 与 ,）后的非参考列
N_BCF_CALLS = 0          # bcftools mpileup | bcftools call -mv 变异记录数
N_HC_CALLS = 0           # HaplotypeCaller 标准模式变异记录数
N_REAL_LIKE = 0          # alt_reads>=10 且 alt_frac>=0.5 的列（真突变形态）
ERR_RATE = 2.192308e-02  # samtools stats error rate（wgsim -e 0.02）
MAX_ALT_FRAC = 0.185     # 全基因组最高 ALT 支持：5 reads / 18.5%（contig2:815）


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
# FIG 1: standard VCF vs gVCF record composition (real run)
# ============================================================
fig, ax = plt.subplots(figsize=(7.4, 4.4))
cats = ["standard HaplotypeCaller\n(raw.vcf.gz)", "gVCF -ERC GVCF\n(raw.g.vcf.gz)"]
var_seg = [0, 0]                    # 变异记录：两模式均为 0 条
ref_seg = [0, N_BLOCKS]             # hom-ref 区块：标准 0 / GVCF 2946
ax.bar(cats, var_seg, color=BRICK, label="variant records")
ax.bar(cats, ref_seg, bottom=var_seg, color=CELADON, label="hom-ref blocks (no-ALT, <NON_REF>)")
ax.text(0, N_BLOCKS * 0.03, "0", ha="center", va="bottom", fontsize=11, color="#222222")
ax.text(1, N_BLOCKS * 1.02, "%d" % N_BLOCKS, ha="center", va="bottom", fontsize=11, color="#222222")
ax.set_ylabel("records")
ax.set_ylim(0, N_BLOCKS * 1.22)
ax.set_title("Same BAM, two modes: standard = 0 variant records, gVCF = %d hom-ref blocks\n"
             "%d blocks with DP>0 (mean DP 35.9), %d with DP=0; all no-ALT, with <NON_REF>"
             % (N_BLOCKS, N_DPP, N_DP0),
             fontsize=10.5)
ax.legend(fontsize=9, loc="upper left")
plt.tight_layout()
total_fail = savefig(fig, "fig1_standard_vs_gvcf.png")

# ============================================================
# FIG 2: GVCF hom-ref block length distribution (real parse)
# ============================================================
bins = [("1", lambda L: L == 1), ("2-5", lambda L: 2 <= L <= 5), ("6-10", lambda L: 6 <= L <= 10),
        ("11-20", lambda L: 11 <= L <= 20), ("21-50", lambda L: 21 <= L <= 50),
        ("51-100", lambda L: 51 <= L <= 100), (">100", lambda L: L > 100)]
counts = [sum(1 for L in lens if cond(L)) for _, cond in bins]
labels = [b[0] for b in bins]
fig, ax = plt.subplots(figsize=(7.4, 4.4))
bars = ax.bar(labels, counts, color=CELADON)
for x, c in enumerate(counts):
    ax.text(x, c + max(counts) * 0.02, str(c), ha="center", va="bottom",
            fontsize=9.5, color="#222222")
ax.set_ylabel("blocks")
ax.set_xlabel("block length (bp)")
ax.set_ylim(0, max(counts) * 1.22)
ax.set_title("GVCF hom-ref block length distribution (%d blocks, 12 kb reference)\n"
             "median %g bp, mean %.1f bp, p99 %d bp, max %d bp (GQ banding)"
             % (N_BLOCKS, MED, MEAN, P99, MAX), fontsize=10.5)
plt.tight_layout()
total_fail += savefig(fig, "fig2_gvcf_block_lengths.png")

# ============================================================
# FIG 3: zero-variant diagnosis (real cross-tool comparison)
# ============================================================
metrics = ["mpileup non-ref\ncols (corrected)",
           "bcftools call -mv\nvariant records",
           "HaplotypeCaller\nvariant records",
           "real-mutation-like\ncolumns"]
vals = [N_NONREF_COLS, N_BCF_CALLS, N_HC_CALLS, N_REAL_LIKE]
colors = [CELADON, BRICK, BRICK, BRICK]
fig, ax = plt.subplots(figsize=(7.6, 4.4))
bars = ax.bar(metrics, vals, color=colors)
for x, v in enumerate(vals):
    ax.text(x, v + max(vals) * 0.02, str(v), ha="center", va="bottom",
            fontsize=10.5, color="#222222")
ax.set_ylabel("count (columns / variant records)")
ax.set_ylim(0, max(vals) * 1.25)
ax.set_title("Zero-variant diagnosis: the input carries no true variant signal\n"
             "wgsim ran with -r 0 -R 0 (no mutations); measured error rate 2.19% matches -e 0.02;\n"
             "max ALT support = 5 reads / 18.5% (contig2:815) — errors, not alleles",
             fontsize=10)
plt.tight_layout()
total_fail += savefig(fig, "fig3_zero_variant_diagnosis.png")

print("FIGURE QUALITY: TOTAL FAILS = %d" % total_fail)
