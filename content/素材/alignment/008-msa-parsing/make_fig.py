#!/usr/bin/env python3
"""Generate 008 figure: per-column gap fraction and conservation fraction
for the repo sample alignment (4 seqs x 21 cols). English labels only
(DejaVu Sans lacks CJK glyphs)."""

import os
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from Bio import AlignIO

HERE = os.path.dirname(os.path.abspath(__file__))
ALN = os.path.join(HERE, "sample_alignment.aln")

alignment = AlignIO.read(ALN, "clustal")
n_seq = len(alignment)
L = alignment.get_alignment_length()

gap_frac = np.array([alignment[:, i].count("-") / n_seq for i in range(L)])
cons = []
for i in range(L):
    col = alignment[:, i]
    c = Counter(col)
    mc, mn = c.most_common(1)[0]
    if mc == "-":
        c.pop("-", None)
        if c:
            mc, mn = c.most_common(1)[0]
        else:
            mn = 0
    cons.append(mn / n_seq)
cons = np.array(cons)

x = np.arange(L)
fig, ax = plt.subplots(figsize=(9, 4.2), dpi=140)
ax.plot(x, gap_frac * 100, "o-", color="#d1495b", lw=1.8, ms=5,
        label="Gap fraction (%)")
ax.plot(x, cons * 100, "s-", color="#2e7d9a", lw=1.8, ms=5,
        label="Conservation (%)")
# mark gappy columns (>50% gap)
for i in range(L):
    if gap_frac[i] >= 0.5:
        ax.axvline(i, color="#d1495b", ls=":", alpha=0.5, lw=1.2)
ax.axhline(50, color="gray", ls="--", lw=0.9, alpha=0.6)
ax.set_xlabel("Alignment column (0-based)")
ax.set_ylabel("Percent")
ax.set_title("MSA per-column profile: 4 seqs x 21 cols (repo sample)\n"
             "18/21 fully conserved; col 12 gappy (3/4 gaps)")
ax.set_ylim(-5, 105)
ax.set_xticks(x)
ax.legend(loc="lower left", frameon=False, fontsize=9)
ax.grid(axis="y", alpha=0.25)
fig.tight_layout()
out = os.path.join(HERE, "008-fig.png")
fig.savefig(out)
print(f"wrote {out}  ({L} cols: max gap%={gap_frac.max()*100:.0f}, "
      f"min cons%={(cons.min()*100):.0f})")
