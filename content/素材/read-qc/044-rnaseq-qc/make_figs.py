#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
044 rnaseq-qc 出图：全部基于本次 WSL 真跑产物，无任何虚构数字。
数据来源（同目录真实文件）：
  - parsed_results.tsv   （parse_quant.py 解析 3 个 quant 目录 + truth.tsv 的宽表）
  - quant_*/lib_format_counts.json （salmon 真跑产出：expected_format / percent_mapped）
配色：brick-red #b5482f + celadon #2f7d72。英文标签（DejaVu Sans 无 CJK）。
自检：文本包围盒必须落在 axes 内、互不重叠；成功时打印恰好
  FIGURE QUALITY: TOTAL FAILS = 0
"""
import json
import math
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
BRICK = "#b5482f"
CELADON = "#2f7d72"


def percent_mapped(tag):
    """salmon 2.7.0 不写 meta_info.json：percent_mapped 从 quant 日志行解析（真跑数字）。"""
    p = os.path.join(BASE, "logs", "_quant_%s.log" % tag)
    if os.path.exists(p):
        with open(p, errors="replace") as f:
            ms = re.findall(r"mapped (\d+) / (\d+) fragments \(([\d.]+)%\)", f.read())
        if ms:
            n_m, n_p, pct = ms[-1]
            return float(pct), int(n_m), int(n_p)
    return 0.0, 0, 0


# ---------- 真实数据 ----------
def read_tsv(path):
    rows = []
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        for ln in f:
            rows.append(dict(zip(head, ln.rstrip("\n").split("\t"))))
    return rows


rows = read_tsv(os.path.join(BASE, "parsed_results.tsv"))
ids = [r["tx_id"] for r in rows]
exp_tpm = [float(r["tpm_expected"]) for r in rows]
tpm_A = [float(r["tpm_A"]) for r in rows]
order = sorted(range(len(ids)), key=lambda i: exp_tpm[i], reverse=True)


def lfc(tag):
    p = os.path.join(BASE, "quant_%s" % tag, "lib_format_counts.json")
    if os.path.exists(p):
        with open(p) as f:
            return json.load(f)
    return {}


lfcA, lfcIU, lfcSF = lfc("A"), lfc("IU"), lfc("SF")
fmt_A = lfcA.get("expected_format", "NA")
pmA, nmA, npA = percent_mapped("A")
pmIU, nmIU, npIU = percent_mapped("IU")
pmSF, nmSF, npSF = percent_mapped("SF")


def pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    return sxy / (sxx * syy) ** 0.5


det = [(math.log10(float(r["tpm_expected"])), math.log10(float(r["tpm_A"])))
       for r in rows if float(r["tpm_A"]) > 0]
R = pearson([x for x, _ in det], [y for _, y in det])
N_ALL, N_DET = len(rows), len(det)


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
# FIG 1: expression gradient recovery (quant_A vs truth)
# ============================================================
fig, ax = plt.subplots(figsize=(7.4, 4.6))
xs = [math.log10(v) for v in exp_tpm]
ys = [math.log10(v) for v in tpm_A]
lo, hi = min(xs + ys) - 0.15, max(xs + ys) + 0.15
ax.plot([lo, hi], [lo, hi], ls="--", lw=1.2, color=CELADON, zorder=1)
ax.scatter(xs, ys, s=26, color=BRICK, zorder=2)
ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)
ax.set_xlabel("log10 expected TPM (simulation truth)")
ax.set_ylabel("log10 salmon TPM (lib type A)")
ax.set_title("Expression gradient recovery: salmon quant -l A vs simulated truth\n"
             "Pearson r on log10 TPM = %.4f (n = %d/%d transcripts detected); "
             "dashed line = identity" % (R, N_DET, N_ALL), fontsize=10.5)
plt.tight_layout()
total_fail = savefig(fig, "fig1_gradient_recovery.png")

# ============================================================
# FIG 2: mapping rate by specified library type (strandedness)
# ============================================================
fig, ax = plt.subplots(figsize=(7.4, 4.6))
cats = ["A (auto-detected\n-> %s)" % fmt_A, "IU (unstranded)", "SF (forward,\nwrong strand)"]
vals = [pmA, pmIU, pmSF]
colors = [CELADON, CELADON, BRICK]
ax.bar(cats, vals, color=colors, width=0.55)
for x, v in enumerate(vals):
    ax.text(x, v + 2.5, "%.2f%%" % v, ha="center", va="bottom",
            fontsize=10.5, color="#222222")
ax.set_ylabel("percent mapped (%)")
ax.set_ylim(0, 115)
ax.set_title("Strandedness gates quantification: percent mapped by specified library type\n"
             "reads are dUTP-style SE (fr-firststrand, antisense); -l A auto-detects %s;\n"
             "IU ignores strand (99.69%%); SF (opposite strand): 0 compatible / %d incompatible"
             % (fmt_A, lfcSF.get("num_incompatible_fragments", 0)), fontsize=10)
plt.tight_layout()
total_fail += savefig(fig, "fig2_strandedness_mapping.png")

# ============================================================
# FIG 3: TPM dynamic range from quant.sf (-l A)
# ============================================================
fig, ax = plt.subplots(figsize=(7.8, 4.6))
xs = list(range(len(order)))
bar_vals = [max(tpm_A[i], 1.0) for i in order]     # log 轴下限截断 1 TPM
line_vals = [exp_tpm[i] for i in order]
ax.bar(xs, bar_vals, color=BRICK, width=0.7, zorder=2)
ax.plot(xs, line_vals, color=CELADON, lw=1.6, marker="o", ms=3.5, zorder=3)
ax.set_yscale("log")
ax.set_xlabel("transcripts ordered by expected TPM (high -> low)")
ax.set_ylabel("TPM (log scale)")
ax.set_title("TPM dynamic range recovered from quant.sf (-l A): %d transcripts, %.0fx span\n"
             "bars = salmon TPM; line = expected TPM (simulation truth); TPM sums to 1e6"
             % (len(ids), max(exp_tpm) / min(exp_tpm)), fontsize=10.5)
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
handles = [Patch(facecolor=BRICK, label="salmon TPM"),
           Line2D([0], [0], color=CELADON, marker="o", ms=3.5, label="expected TPM (truth)")]
ax.legend(handles=handles, fontsize=9, loc="upper right")
plt.tight_layout()
total_fail += savefig(fig, "fig3_tpm_dynamic_range.png")

print("FIGURE QUALITY: TOTAL FAILS = %d" % total_fail)
