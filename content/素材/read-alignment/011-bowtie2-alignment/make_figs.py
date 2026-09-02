#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
011 bowtie2-alignment 出图：基于真实复现产物（bowtie2_results.json + samtools flagstat *.txt）。
配色：brick-red #b5482f + celadon #2f7d72。英文标签（DejaVu Sans 无 CJK）。
出图质量自检（替代 check_figs.py）：文本包围盒落 axes 内、不重叠。
真实数据来源：
  - bowtie2_results.json            （run_bowtie2.py 真跑产出）
  - flagstat_aligned_e2e.txt 等 5 个 samtools flagstat 文件（本次真跑重新抽取）
"""
import os
import re
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(BASE, "bowtie2_results.json")))

BRICK = "#b5482f"
CELADON = "#2f7d72"
GREY = "#6b7280"


def _flagstat(path):
    """Parse samtools flagstat text for real mapped % and properly-paired %."""
    txt = open(os.path.join(BASE, path)).read()
    m = re.search(r"^\s*(\d+)\s*\+\s*0 in total", txt, re.M)
    total = int(m.group(1)) if m else None
    m = re.search(r"^\s*(\d+)\s*\+\s*0 mapped \(([\d.]+)%", txt, re.M)
    mapped = (int(m.group(1)), float(m.group(2))) if m else (None, None)
    m = re.search(r"^\s*(\d+)\s*\+\s*0 properly paired \(([\d.]+)%", txt, re.M)
    pp = (int(m.group(1)), float(m.group(2))) if m else (None, None)
    return total, mapped, pp


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
# FIG 1: Alignment rate under adapter contamination (real run)
# ============================================================
clean_e2e = R["e2e"]["overall_rate"]                 # 99.23
contam_e2e = R["adapter_contam"]["e2e_rate"]         # 49.28
contam_local = R["adapter_contam"]["local_rate"]     # 99.77
recov = R["adapter_contam"]["recovered_by_local"]    # 50.49
fig, ax = plt.subplots(figsize=(8.4, 5.0), dpi=150)
groups = ["Clean reads\nend-to-end", "Adapter-contam\nend-to-end", "Adapter-contam\n--local (soft-clip)"]
vals = [clean_e2e, contam_e2e, contam_local]
colors = [CELADON, BRICK, CELADON]
bars = ax.bar(groups, vals, color=colors, zorder=3, width=0.6)
ax.set_ylabel("alignment rate (%)", fontsize=10.5)
ax.set_title("Adapter contamination drops end-to-end rate;\n--local soft-clip recovers it",
             fontsize=12, pad=12)
ax.set_ylim(0, 115)
for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width() / 2, v + 2, "%.2f%%" % v,
            ha="center", va="bottom", fontsize=10, color="#222222")
ax.annotate("+%.2f pp recovered\nby soft-clipping read ends" % recov,
            xy=(2, contam_local), xytext=(1.45, 88), textcoords="data",
            ha="center", va="center", fontsize=9.5, color=BRICK,
            bbox=dict(boxstyle="round,pad=0.3", fc="#fff3f0", ec=BRICK, lw=0.8),
            arrowprops=dict(arrowstyle="->", color=BRICK, lw=1.2))
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "Source: bowtie2_results.json e2e / adapter_contam (run_bowtie2.py, real execution)")
f1 = savefig(fig, "fig1_adapter_contam_rate.png")

# ============================================================
# FIG 2: MAPQ ceiling — Bowtie2 scale never reaches BWA 60 (real run)
# ============================================================
e2e_max = R["mapq_cap"]["e2e_max"]      # 42
local_max = R["mapq_cap"]["local_max"]  # 44
bwa = R["mapq_cap"]["bwa_equivalent"]   # 60
fig, ax = plt.subplots(figsize=(8.4, 5.0), dpi=150)
labels = ["e2e (max)", "local (max)", "BWA (equiv.)"]
caps = [e2e_max, local_max, bwa]
bars = ax.bar(labels, caps, color=[CELADON, CELADON, GREY], zorder=3, width=0.55)
ax.set_ylabel("max MAPQ", fontsize=10.5)
ax.set_title("Bowtie2 MAPQ ceiling: 42 (e2e) / 44 (local)\nnever reaches BWA's 60",
             fontsize=12, pad=12)
ax.set_ylim(0, 75)
for b, v in zip(bars, caps):
    ax.text(b.get_x() + b.get_width() / 2, v + 1.5, str(v),
            ha="center", va="bottom", fontsize=10, color="#222222")
ax.annotate("-q 60 filter would\nempty the BAM",
            xy=(2, bwa), xytext=(1.3, 30), textcoords="data",
            ha="center", va="center", fontsize=9.5, color=BRICK,
            bbox=dict(boxstyle="round,pad=0.3", fc="#fff3f0", ec=BRICK, lw=0.8),
            arrowprops=dict(arrowstyle="->", color=BRICK, lw=1.2))
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "Source: bowtie2_results.json mapq_cap (samtools view max MAPQ, real BAMs)")
f2 = savefig(fig, "fig2_mapq_ceiling.png")

# ============================================================
# FIG 3: Sensitivity presets (real run) + properly-paired rate (flagstat)
# ============================================================
pr = R["presets"]
pk = ["--very-fast", "--sensitive", "--very-sensitive"]
pv = [pr[k] for k in pk]
# real properly-paired % from fresh samtools flagstat
_, _, pp_e2e = _flagstat("flagstat_aligned_e2e.txt")
_, _, pp_local = _flagstat("flagstat_contam_local.txt")
fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.2, 4.8), dpi=150)
# left: preset alignment rate
bl = axL.bar([k.replace("--", "") for k in pk], pv, color=CELADON, zorder=3, width=0.55)
axL.set_ylabel("alignment rate (%)", fontsize=10.5)
axL.set_title("Sensitivity presets\nend-to-end, real run", fontsize=11.5, pad=10)
axL.set_ylim(90, 102)
for b, v in zip(bl, pv):
    axL.text(b.get_x() + b.get_width() / 2, v + 0.25, "%.2f%%" % v,
             ha="center", va="bottom", fontsize=9.5, color="#222222")
axL.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
# right: properly-paired rate (fragment geometry)
pp_vals = [pp_e2e[1], pp_local[1]]
br = axR.bar(["aligned_e2e", "contam_local"], pp_vals, color=BRICK, zorder=3, width=0.5)
axR.set_ylabel("properly paired (%)", fontsize=10.5)
axR.set_title("Properly-paired rate\n(fragment geometry, real flagstat)", fontsize=11.5, pad=10)
axR.set_ylim(0, 65)
for b, v, n in zip(br, pp_vals, [pp_e2e[0], pp_local[0]]):
    axR.text(b.get_x() + b.get_width() / 2, v + 1.5,
             "%.2f%%\n(%d reads)" % (v, n),
             ha="center", va="bottom", fontsize=9.5, color="#222222")
axR.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
fig.suptitle("bowtie2-alignment — faithful reproduction (real artifacts)", fontsize=12)
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
f3 = savefig(fig, "fig3_presets_properlypaired.png")

total_fail = f1 + f2 + f3
print("FIGURE QUALITY: TOTAL FAILS = %d" % total_fail)
import sys
sys.exit(1 if total_fail else 0)
