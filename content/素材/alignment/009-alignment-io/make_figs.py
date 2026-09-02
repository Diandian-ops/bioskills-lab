#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
009 alignment-io 出图：基于 run.py 真跑产出的 alignment_io_data.json。
配色：brick-red #b5482f + celadon #2f7d72。英文标签。_verify 自检。
"""
import os
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(BASE, "alignment_io_data.json")))

BRICK = "#b5482f"
CELADON = "#2f7d72"
GREY = "#6b7280"

def _verify(fig, name):
    fig.canvas.draw(); r = fig.canvas.get_renderer(); fails = 0
    for ax in fig.axes:
        axb = ax.get_window_extent(r); bbs = []
        for t in ax.texts:
            bb = t.get_window_extent(r)
            if (bb.x0 < axb.x0 - 3 or bb.x1 > axb.x1 + 3 or bb.y0 < axb.y0 - 3 or bb.y1 > axb.y1 + 3):
                print("  [FAIL] %s: text out of axes: %r" % (name, t.get_text()[:50])); fails += 1
            bbs.append(bb)
        for i in range(len(bbs)):
            for j in range(i + 1, len(bbs)):
                a, b = bbs[i], bbs[j]
                if not (a.x1 + 1 < b.x0 or b.x1 + 1 < a.x0 or a.y1 + 1 < b.y0 or b.y1 + 1 < a.y0):
                    print("  [FAIL] %s: text overlap: %r / %r" % (name, ax.texts[i].get_text()[:30], ax.texts[j].get_text()[:30])); fails += 1
    if fails == 0: print("  [PASS] %s: text within axes, no overlap" % name)
    return fails

def savefig(fig, name):
    fails = _verify(fig, name)
    fig.savefig(os.path.join(BASE, name), dpi=150); plt.close(fig); return fails

def footnote(fig, text):
    fig.text(0.5, 0.012, text, ha="center", va="bottom", fontsize=8.5, color=GREY)

# ============================================================
# 图 1：注解保留（Stockholm 主 vs FASTA 导出）
# ============================================================
ann = D["annotation"]
fig, ax = plt.subplots(figsize=(8.0, 5.2), dpi=150)
paths = ["Stockholm\n(master)", "FASTA\n(export)"]
retained = [1.0 if ann["stockholm_keeps_ss"] else 0.0,
            0.0 if ann["fasta_loses_ss"] else 1.0]
colors = [CELADON if v > 0.5 else BRICK for v in retained]
bars = ax.bar(paths, retained, color=colors, width=0.5, zorder=3)
ax.set_ylim(0, 1.25)
ax.set_ylabel("SS_cons (GC) annotation retained", fontsize=10.5)
ax.set_title("Annotation preservation: Stockholm keeps, FASTA discards\n(msa-io: Stockholm is the only format preserving GS/GR/GC)",
             fontsize=12, pad=12)
for p, v in zip(paths, retained):
    ax.text(p, v + 0.05, "KEPT" if v > 0.5 else "LOST",
            ha="center", va="bottom", fontsize=11,
            color=CELADON if v > 0.5 else BRICK, fontweight="bold")
ax.grid(axis="y", alpha=0.2, ls=":", lw=0.7, zorder=0)
footnote(fig, "Stockholm #=GC SS_cons round-trips; AlignIO.write(...,'fasta') silently drops every annotation")
f1 = savefig(fig, "fig1_annotation_survival.png")

# ============================================================
# 图 2：格式读/写支持矩阵（Bio.AlignIO 实测）
# ============================================================
sup = D["format_support"]
fmts = list(sup.keys())
read_v = [1 if sup[f]["read"] else 0 for f in fmts]
write_v = [1 if sup[f]["write"] else 0 for f in fmts]
y = range(len(fmts))[::-1]
fig, ax = plt.subplots(figsize=(8.6, 5.6), dpi=150)
import numpy as np
ax.barh([yi - 0.2 for yi in y], read_v, height=0.38, color=CELADON, label="read", zorder=3)
ax.barh([yi + 0.2 for yi in y], write_v, height=0.38, color=BRICK, label="write", zorder=3)
ax.set_yticks(list(y)); ax.set_yticklabels(fmts, fontsize=10)
ax.set_xlim(0, 1.35); ax.set_xticks([0, 1]); ax.set_xticklabels(["no", "yes"])
ax.set_xlabel("supported", fontsize=10.5)
ax.set_title("Bio.AlignIO read/write support (7 formats tested)\nall tested formats support both read and write (DNA MSA, molecule_type set)",
             fontsize=12, pad=12)
ax.legend(loc="lower right", fontsize=9)
for yi, rv, wv in zip(y, read_v, write_v):
    ax.text(rv + 0.03, yi - 0.2, "R" if rv else "-", va="center", fontsize=9, color=CELADON)
    ax.text(wv + 0.03, yi + 0.2, "W" if wv else "-", va="center", fontsize=9, color=BRICK)
ax.grid(axis="x", alpha=0.2, ls=":", lw=0.7, zorder=0)
footnote(fig, "nexus write/convert REQUIRES molecule_type; without it ValueError. A2M/A3M/MAF read via Bio.AlignIO only")
f2 = savefig(fig, "fig2_format_support.png")

total_fail = f1 + f2
print("FIGURE QUALITY: TOTAL FAILS = %d" % total_fail)
import sys
sys.exit(1 if total_fail else 0)
