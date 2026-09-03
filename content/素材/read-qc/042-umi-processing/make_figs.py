#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
042 umi-processing 出图：全部基于本次 WSL 真跑产物，无任何虚构数字。
数据来源（同目录真实文件）：
  - results.json        （parse_results.py 产出：truth / extract / 各方法分子数）
  - stats_dir_edit_distance.tsv        （umi_tools dedup --output-stats 真跑产出）
  - stats_dir_per_umi_per_position.tsv （同上）
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
R = json.load(open(os.path.join(BASE, "results.json")))

BRICK = "#b5482f"
CELADON = "#2f7d72"
GREY = "#6b7280"
DARK = "#222222"

TRUTH = R["truth"]["n_molecules"]                       # 6000
M = R["methods"]                                        # 各方法分子估计
ERR = {k: R["methods_err_pct_" + k] for k in M}
ED = R["edit_distance"]                                 # [dist, u_obs, u_null, d_obs, d_null]
FAM = R["family_size"]                                  # [size, pre, post]
RET = R["retention_pct_directional"]                    # 36.56


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
# FIG 1: molecule estimates by strategy vs simulated truth
# ============================================================
labels = ["Simulated\ntruth", "coordinate-only\n(samtools markdup)",
          "unique\n(exact UMI)", "directional\n(umi_tools default)",
          "cluster\n(single-linkage)"]
keys = ["truth", "coordinate_only", "unique", "directional", "cluster"]
vals = [TRUTH, M["coordinate_only"], M["unique"], M["directional"], M["cluster"]]
cols = [GREY, BRICK, BRICK, CELADON, BRICK]
fig, ax = plt.subplots(figsize=(7.6, 4.6))
ax.bar(labels, vals, color=cols, width=0.62)
for x, v in enumerate(vals):
    tag = "reference" if x == 0 else "%d (%+.2f%%)" % (v, ERR[keys[x]])
    ax.text(x, v + max(vals) * 0.03, tag, ha="center", va="bottom",
            fontsize=9.5, color=DARK)
ax.set_ylabel("estimated molecules")
ax.set_ylim(0, max(vals) * 1.24)
ax.set_title("Molecule counting: (coordinate + UMI) beats coordinate alone\n"
             "6000 simulated molecules, 16445 read pairs, seed 42",
             fontsize=10.5)
plt.tight_layout()
total_fail = savefig(fig, "fig1_molecule_estimates.png")

# ============================================================
# FIG 2: UMI edit-distance distribution, unique vs directional
# ============================================================
ds = [row[0] for row in ED if row[0].isdigit() and int(row[0]) <= 8]
ds = sorted(int(d) for d in ds)
u_obs = [next(r[1] for r in ED if r[0] == str(d)) for d in ds]
d_obs = [next(r[3] for r in ED if r[0] == str(d)) for d in ds]
import numpy as np
x = np.arange(len(ds))
w = 0.38
fig, ax = plt.subplots(figsize=(7.6, 4.6))
ax.bar(x - w / 2, u_obs, w, color=BRICK, label="unique (no error model)")
ax.bar(x + w / 2, d_obs, w, color=CELADON, label="directional (default)")
ax.set_xticks(x)
ax.set_xticklabels([str(d) for d in ds])
ax.set_xlabel("edit distance between UMI pairs at the same position")
ax.set_ylabel("UMI pairs")
for xi, v in zip(x - w / 2, u_obs):
    if v > 0:
        ax.text(xi, v + max(u_obs) * 0.02, str(v), ha="center", va="bottom",
                fontsize=9, color=DARK)
for xi, v in zip(x + w / 2, d_obs):
    if v > 0:
        ax.text(xi, v + max(u_obs) * 0.02, str(v), ha="center", va="bottom",
                fontsize=9, color=DARK)
ax.set_ylim(0, max(u_obs) * 1.22)
ax.set_title("Directional folds 1-off UMI errors back into their parent\n"
             "at d=1: unique counts 674 pairs (null = 0), directional 102",
             fontsize=10.5)
ax.legend(fontsize=9, loc="upper right")
plt.tight_layout()
total_fail += savefig(fig, "fig2_edit_distance.png")

# ============================================================
# FIG 3: UMI family size distribution before vs after dedup
# ============================================================
sizes = [row[0] for row in FAM]
pre = [row[1] for row in FAM]
post = [row[2] for row in FAM]
y = np.arange(len(sizes))
fig, ax = plt.subplots(figsize=(7.6, 4.8))
ax.barh(y - w / 2, pre, w, color=BRICK, label="before dedup (exact UMI: 7036)")
ax.barh(y + w / 2, post, w, color=CELADON, label="after directional (6013)")
ax.set_yticks(y)
ax.set_yticklabels([str(s) for s in sizes])
ax.set_ylabel("reads per (position, UMI) family")
ax.set_xlabel("number of UMI families")
for yi, v in zip(y - w / 2, pre):
    ax.text(v + max(pre) * 0.012, yi, str(v), va="center", ha="left",
            fontsize=8.5, color=DARK)
for yi, v in zip(y + w / 2, post):
    ax.text(v + max(pre) * 0.012, yi, str(v), va="center", ha="left",
            fontsize=8.5, color=DARK)
ax.set_xlim(0, max(pre) * 1.16)
ax.set_title("PCR duplicates collapse: family sizes before vs after dedup\n"
             "16445 pairs -> 6013 molecules; 869 UMI-error copies folded in",
             fontsize=10.5)
ax.legend(fontsize=9, loc="upper right")
plt.tight_layout()
total_fail += savefig(fig, "fig3_family_size.png")

print("FIGURE QUALITY: TOTAL FAILS = %d" % total_fail)
