#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""036 reference-operations 出图：基于 036 fig_data.json 真实 contig 长度/GC。"""
import os, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(BASE, "fig_data.json")))
BRICK = "#b5482f"; CELADON = "#2f7d72"; GREY = "#6b7280"

def _verify(fig, name):
    fig.canvas.draw(); r = fig.canvas.get_renderer(); fails = 0
    for ax in fig.axes:
        axb = ax.get_window_extent(r); bbs = []
        for t in ax.texts:
            bb = t.get_window_extent(r)
            if (bb.x0 < axb.x0 - 3 or bb.x1 > axb.x1 + 3 or bb.y0 < axb.y0 - 3 or bb.y1 > axb.y1 + 3):
                print("  [FAIL] %s: %r" % (name, t.get_text()[:50])); fails += 1
            bbs.append(bb)
        for i in range(len(bbs)):
            for j in range(i + 1, len(bbs)):
                a, b = bbs[i], bbs[j]
                if not (a.x1 + 1 < b.x0 or b.x1 + 1 < a.x0 or a.y1 + 1 < b.y0 or b.y1 + 1 < a.y0):
                    print("  [FAIL] %s: overlap" % name); fails += 1
    if fails == 0: print("  [PASS] %s" % name)
    return fails

def savefig(fig, name):
    fails = _verify(fig, name)
    fig.savefig(os.path.join(BASE, name), dpi=150); plt.close(fig); return fails

def footnote(fig, text):
    fig.text(0.5, 0.012, text, ha="center", va="bottom", fontsize=8.5, color=GREY)

contigs = D["contigs"]; lengths = D["lengths"]; gc = D["gc"]
fig, ax = plt.subplots(figsize=(8.6, 4.8), dpi=150)
ax.bar(contigs, lengths, color=CELADON, width=0.55, zorder=3)
ax.set_ylabel("contig length (bp)", fontsize=10.5)
ax.set_title("Reference contig lengths (samtools faidx)", fontsize=12, pad=12)
ax.set_ylim(0, max(lengths) * 1.20)
for c, v in zip(contigs, lengths):
    ax.text(c, v + max(lengths) * 0.015, str(v), ha="center", va="bottom", fontsize=9.5, color="#222222")
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "faidx builds .fai; region extraction (contig1:1-60) uses it. samtools dict is absent (use picard/gatk)")
f1 = savefig(fig, "fig1_contig_lengths.png")

fig, ax = plt.subplots(figsize=(8.6, 4.8), dpi=150)
ax.bar(contigs, gc, color=BRICK, width=0.55, zorder=3)
ax.set_ylabel("GC content (%)", fontsize=10.5)
ax.set_title("GC content per contig", fontsize=12, pad=12)
ax.set_ylim(0, 100)
for c, v in zip(contigs, gc):
    ax.text(c, v + 1.5, "%.1f" % v, ha="center", va="bottom", fontsize=9.5, color="#222222")
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "GC skew affects k-mer / mapping bias; reference GC is a first-pass QC signal")
f2 = savefig(fig, "fig2_gc_content.png")

total = f1 + f2
print("FIGURE QUALITY: TOTAL FAILS = %d" % total)
import sys; sys.exit(1 if total else 0)
