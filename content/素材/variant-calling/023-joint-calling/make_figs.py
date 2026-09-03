"""023 joint-calling — real-data figures (1000 Genomes Phase3 chr22 slice).

Data: chr22_slice.vcf.gz (chr22:17.0-17.2 Mb, GRCh37/hs37d5, 2504 samples, 5431 sites).
Demonstrates joint calling with bcftools: split per-sample callsets, naive bcftools
merge, then quantify backfill vs a proper joint call (squared matrix).
Reads joint_calling_stats.json (produced by run.sh) + merged_naive.vcf.gz (same dir).
Usage: python make_figs.py   (run from this directory or any dir)
"""
import os
import json
import gzip

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
VCF = os.path.join(HERE, "merged_naive.vcf.gz")
JSON = os.path.join(HERE, "joint_calling_stats.json")

C_RED = "#b5482f"
C_TEAL = "#2f7d72"
C_GREY = "#555555"
C_LIGHT = "#f5f5f0"

# ------------------------------------------------------------
# savefig 自检包装（参考 004/012 模式）：文本须落在 axes 内且互不重叠，
# 结束时精确打印 "FIGURE QUALITY: TOTAL FAILS = 0" 才算通过
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
    fig.savefig(os.path.join(HERE, name), dpi=130)
    plt.close(fig)
    return fails

with open(JSON, encoding="utf-8") as f:
    D = json.load(f)

samples = D["samples"]
per = D["per_sample"]
M = D["merged_records"]
L = D["joint_records"]

variant = [per[s]["variant_sites"] for s in samples]
backfill = [per[s]["backfilled"] for s in samples]

# ================= FIG 1 — per-sample variant-site count =================
fig1, ax = plt.subplots(figsize=(8.2, 4.2))
x = range(len(samples))
bars = ax.bar(x, variant, color=C_TEAL, width=0.6)
for b, v in zip(bars, variant):
    ax.text(b.get_x() + b.get_width() / 2, v + max(variant) * 0.01, str(v),
            ha="center", fontsize=9, color="#222")
ax.set_xticks(list(x))
ax.set_xticklabels(samples, rotation=20, ha="right", fontsize=9)
ax.set_ylabel("variant sites in single-sample callset")
ax.set_title("Per-sample variant sites (simulated single-sample callsets)\n"
             "each callset reports only loci where that sample carries a variant", fontsize=11)
ax.set_ylim(0, max(variant) * 1.18)
ax.spines[["top", "right"]].set_visible(False)
ax.set_facecolor(C_LIGHT)
fig1.tight_layout()
TOTAL_FAILS = savefig(fig1, "fig1_per_sample_sites.png")

# ================= FIG 2 — backfilled genotypes per sample =================
fig2, ax = plt.subplots(figsize=(8.2, 4.2))
bars = ax.bar(x, backfill, color=C_RED, width=0.6)
for b, v in zip(bars, backfill):
    ax.text(b.get_x() + b.get_width() / 2, v + max(backfill) * 0.01, str(v),
            ha="center", fontsize=9, color="#222")
ax.set_xticks(list(x))
ax.set_xticklabels(samples, rotation=20, ha="right", fontsize=9)
ax.set_ylabel("genotypes (./*) recovered by joint calling")
ax.set_title("Per-sample missing genotypes in naive merge\n"
             "filled to 0/0 by joint calling (backfill = ./. count per sample)", fontsize=11)
ax.set_ylim(0, max(backfill) * 1.18)
ax.spines[["top", "right"]].set_visible(False)
ax.set_facecolor(C_LIGHT)
fig2.tight_layout()
TOTAL_FAILS += savefig(fig2, "fig2_backfill.png")

# ================= FIG 3 — genotype coverage matrix of naive merge =================
samples_vcf, matrix = [], []
with gzip.open(VCF, "rt") as fh:
    for line in fh:
        if line.startswith("##"):
            continue
        if line.startswith("#"):
            f = line.rstrip("\n").split("\t")
            samples_vcf = f[9:]
            continue
        f = line.rstrip("\n").split("\t")
        fmt = f[8].split(":")
        gi = fmt.index("GT")
        row = []
        for col in f[9:]:
            gt = col.split(":")[gi] if gi < len(col.split(":")) else "."
            row.append(0 if gt in ("./.", ".") else 1)
        matrix.append(row)

matrix = list(map(list, zip(*matrix)))  # now [sample][site]
nsite = len(matrix[0])
# cap columns for readability
if nsite > 110:
    step = nsite / 100
    idx = [int(i * step) for i in range(100)]
    matrix = [[r[i] for i in idx] for r in matrix]
    nsite = 100

cmap = ListedColormap([C_RED, C_TEAL])
fig3, ax = plt.subplots(figsize=(9.0, 3.6))
ax.imshow(matrix, cmap=cmap, aspect="auto", vmin=0, vmax=1)
ax.set_yticks(range(len(samples_vcf)))
ax.set_yticklabels(samples_vcf, fontsize=8)
ax.set_xlabel("loci (union of variant sites across the 6 samples)")
ax.set_title("Genotype matrix of naive bcftools merge\n"
             "teal = called genotype, red = ./. (filled to 0/0 by proper joint calling)", fontsize=11)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=C_TEAL, label="called"),
                   Patch(color=C_RED, label="./. (missing)")],
          loc="upper right", fontsize=8, frameon=False)
fig3.tight_layout()
TOTAL_FAILS += savefig(fig3, "fig3_matrix.png")
print("FIGURE QUALITY: TOTAL FAILS = %d" % TOTAL_FAILS)
import sys
sys.exit(1 if TOTAL_FAILS else 0)
