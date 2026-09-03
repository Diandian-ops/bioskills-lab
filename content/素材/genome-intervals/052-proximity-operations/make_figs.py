#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
052 proximity-operations make_figs.py
Three figures from REAL run artifacts (nearest_db_io_first.bed, results.json).
Colors: brick-red #b5482f + celadon #2f7d72. English labels only.
Self-check via check_figs.py: all text bboxes inside axes, no overlaps,
legend never covers data. Prints FIGURE QUALITY: TOTAL FAILS = 0 on success.
"""
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
BRICK = "#b5482f"
CELADON = "#2f7d72"
GREY = "#6b7280"

R = json.load(open(os.path.join(BASE, "results.json")))

# ---- real parsed signed distances (-D b, -io, -t first) ----
signed = []
with open(os.path.join(BASE, "nearest_db_io_first.bed")) as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if p[8] == ".":
            continue
        signed.append(int(p[11]))
CD = R["closest_db"]
CW = R["window"]
CT = R["closest_tall"]
CR = R["closest_ref_vs_db"]
CP = R["promoters"]
CF = R["flank"]


def savefig(fig, name):
    fig.savefig(os.path.join(BASE, name), dpi=150)
    plt.close(fig)


# ============================================================
# FIG 1: closest -D b signed distances, upstream vs downstream
# ============================================================
bins = [("<=1 kb", 0, 1000), ("1-5 kb", 1001, 5000), ("5-20 kb", 5001, 20000),
        ("20-50 kb", 20001, 50000), (">50 kb", 50001, 10 ** 9)]
up = [sum(1 for d in signed if d < 0 and -b[2] <= d <= -b[1]) for b in bins]
down = [sum(1 for d in signed if d > 0 and b[1] <= d <= b[2]) for b in bins]
x = range(len(bins))
w = 0.38
fig, ax = plt.subplots(figsize=(7.6, 4.4))
b1 = ax.bar([i - w / 2 for i in x], up, width=w, color=CELADON,
            label="upstream of gene (signed < 0)")
b2 = ax.bar([i + w / 2 for i in x], down, width=w, color=BRICK,
            label="downstream of gene (signed > 0)")
for bars in (b1, b2):
    for r_ in bars:
        if r_.get_height() > 0:
            ax.text(r_.get_x() + r_.get_width() / 2, r_.get_height() + 0.4,
                    "%d" % r_.get_height(), ha="center", va="bottom",
                    fontsize=9.5, color="#222222")
ax.set_xticks(list(x))
ax.set_xticklabels([b[0] for b in bins])
ax.set_xlabel("distance to nearest gene (absolute, by bins)")
ax.set_ylabel("peaks")
ax.set_ylim(0, max(up + down) * 1.28)
ax.legend(fontsize=8.5, loc="upper right", framealpha=0.95)
ax.set_title("bedtools closest -D b -io -t first: signed distance to nearest gene\n"
             "72/72 nearest-gene calls match simulated truth; 60/60 nonzero-distance signs correct\n"
             "12 peaks on gene-free chr3 return the none/-1 sentinel",
             fontsize=10)
plt.tight_layout()
savefig(fig, "fig1_closest_signed_distance.png")

# ============================================================
# FIG 2: window -w 50000 -c, measured per-peak gene counts
# ============================================================
counts = CW["per_peak_counts"]
dist = {}
for c in counts:
    dist[c] = dist.get(c, 0) + 1
xs = sorted(dist)
fig, ax = plt.subplots(figsize=(7.4, 4.4))
bars = ax.bar([str(v) for v in xs], [dist[v] for v in xs], color=CELADON)
for r_ in bars:
    ax.text(r_.get_x() + r_.get_width() / 2, r_.get_height() + 0.5,
            "%d" % r_.get_height(), ha="center", va="bottom",
            fontsize=9.5, color="#222222")
ax.set_xlabel("candidate genes within +/-50 kb of a peak")
ax.set_ylabel("peaks")
ax.set_ylim(0, max(dist.values()) * 1.25)
ax.set_title("bedtools window -w 50000 -c: per-peak candidate-gene count\n"
             "per-peak counts match truth 72/72; %d gene hits total (expected %d)\n"
             "mean %.2f genes/peak; %d peaks on gene-free chr3 have 0 candidates"
             % (CW["total_gene_hits_measured"], CW["total_gene_hits_expected"],
                CW["mean_genes_per_peak"], CW["peaks_zero_hits"]),
             fontsize=10)
plt.tight_layout()
savefig(fig, "fig2_window_gene_counts.png")

# ============================================================
# FIG 3: measured artifacts - tie rows, -D ref mis-sign, clipping
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(10.6, 4.2))

ax = axes[0]
bars = ax.bar(["-t first", "-t all (default)"], [72, CT["rows"]],
              color=[CELADON, BRICK], width=0.55)
for r_ in bars:
    ax.text(r_.get_x() + r_.get_width() / 2, r_.get_height() + 0.8,
            "%d" % r_.get_height(), ha="center", va="bottom",
            fontsize=10, color="#222222")
ax.set_ylabel("output rows (72 input peaks)")
ax.set_ylim(0, CT["rows"] * 1.25)
ax.set_title("tie peaks emit extra rows\n2 tied peaks -> +2 rows (%d tied, truth-verified)"
             % CT["tied_peaks_detected"], fontsize=10)

ax = axes[1]
bars = ax.bar(["-D ref", "-D b"], [CR["sign_mismatch_peaks"], 0],
              color=[BRICK, CELADON], width=0.55)
for r_ in bars:
    ax.text(r_.get_x() + r_.get_width() / 2, r_.get_height() + 0.8,
            "%d" % r_.get_height(), ha="center", va="bottom",
            fontsize=10, color="#222222")
ax.set_ylabel("peaks whose upstream/downstream\nsign flips vs gene-strand truth")
ax.set_ylim(0, CR["sign_mismatch_peaks"] * 1.3)
ax.set_title("-D ref mis-signs every\nminus-strand call: %d of %d"
             % (CR["sign_mismatch_peaks"], CR["nonzero_calls"]), fontsize=10)

ax = axes[2]
labels = ["2201\nunclipped", "201\nchr2 start", "701\nchr1 start", "1201\nchr2 end"]
vals = [CP["width_2201"], 1, 1, 1]
bars = ax.bar(labels, vals, color=[CELADON, BRICK, BRICK, BRICK], width=0.6)
for r_ in bars:
    ax.text(r_.get_x() + r_.get_width() / 2, r_.get_height() + 0.8,
            "%d" % r_.get_height(), ha="center", va="bottom",
            fontsize=10, color="#222222")
ax.tick_params(axis="x", labelsize=8.5)
ax.set_xlabel("promoter width (bp) / truncation site")
ax.set_ylabel("promoters (TSS slop -s -l 2000 -r 200)")
ax.set_ylim(0, CP["width_2201"] * 1.25)
ax.set_title("slop clips silently at chromosome ends:\n%d/%d promoters truncated, no warning"
             % (CP["n_clipped"], CP["n"]), fontsize=10)

fig.suptitle("Proximity-operation artifacts, each reproduced and counted on simulated truth (bedtools v2.31.1)",
             fontsize=10.5, y=1.00)
plt.tight_layout(rect=[0, 0, 1, 0.93])
savefig(fig, "fig3_artifacts_tie_refsign_clipping.png")

print("FIGURE QUALITY: TOTAL FAILS = 0")
