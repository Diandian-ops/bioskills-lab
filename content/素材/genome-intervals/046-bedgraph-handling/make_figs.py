#!/usr/bin/env python3
"""046-bedgraph-handling: three figures from the real bedGraph/bigWig run.

fig1  raw coverage tracks of sample A vs sample B over the peak2 window
      (identical biology, ~4x deeper library) - the library-size artifact
fig2  the same window after manual RPM scaling (-scale 1e6/nreads)
fig3  bigWig round-trip fidelity: bedGraph value vs value read back from A.bw
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RED = "#b5482f"
GREEN = "#2f7d72"
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

W0, W1 = 780000, 820000  # window around peak2 (center 800000)


def window_rows(path):
    rows = []
    with open(path) as f:
        for ln in f:
            c, s, e, v = ln.split()
            s, e = int(s), int(e)
            if e > W0 and s < W1:
                rows.append((s, e, float(v)))
    return rows


def step_xy(rows):
    xs, ys = [], []
    for s, e, v in rows:
        xs.extend([s, e])
        ys.extend([v, v])
    return xs, ys


def fmt_kb(x, pos):
    return "%d" % (x // 1000)


# ---------------------------------------------------------------- fig 1
fig, ax = plt.subplots(figsize=(8.0, 4.8), dpi=150)
rows_a = window_rows("A.bedgraph")
rows_b = window_rows("B.bedgraph")
xs, ys = step_xy(rows_a)
ax.plot([x / 1000.0 for x in xs], ys, color=RED, lw=1.2,
        label="sample A raw (%d reads)" % json.load(open("meta.json"))["n_reads_A"])
xs, ys = step_xy(rows_b)
ax.plot([x / 1000.0 for x in xs], ys, color=GREEN, lw=1.2,
        label="sample B raw (%d reads, same biology)" % json.load(open("meta.json"))["n_reads_B"])
ax.set_xlim(W0 / 1000.0, W1 / 1000.0)
ax.set_ylim(0, max(v for _s, _e, v in rows_a) * 1.15)
ax.set_xlabel("position on chr1 (kb)")
ax.set_ylabel("raw bedGraph coverage (x)")
ax.set_title("Raw bedGraph heights track library size, not biology")
ax.legend(loc="upper right", frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig("fig1_library_size_artifact.png")
plt.close(fig)

# ---------------------------------------------------------------- fig 2
fig, ax = plt.subplots(figsize=(8.0, 4.8), dpi=150)
rows_as = window_rows("A_scaled.bedgraph")
rows_bs = window_rows("B_scaled.bedgraph")
xs, ys = step_xy(rows_as)
ax.plot([x / 1000.0 for x in xs], ys, color=RED, lw=1.2,
        label="sample A, RPM-scaled")
xs, ys = step_xy(rows_bs)
ax.plot([x / 1000.0 for x in xs], ys, color=GREEN, lw=1.2, alpha=0.85,
        label="sample B, RPM-scaled")
ymax = max(max(v for _s, _e, v in rows_as), max(v for _s, _e, v in rows_bs))
ax.set_xlim(W0 / 1000.0, W1 / 1000.0)
ax.set_ylim(0, ymax * 1.18)
ax.set_xlabel("position on chr1 (kb)")
ax.set_ylabel("RPM-scaled coverage (per 1e6 reads)")
ax.set_title("RPM scaling (1e6 / mapped reads) puts both samples on one scale")


def center_mean(rows, c=800000, half=250):
    tot = 0.0
    ln = 0
    for s, e, v in rows:
        lo, hi = max(s, c - half), min(e, c + half)
        if hi > lo:
            tot += v * (hi - lo)
            ln += hi - lo
    return tot / ln if ln else 0.0


ma = center_mean(rows_as)
mb = center_mean(rows_bs)
ax.text(0.03, 0.93, "peak-center mean (799.75-800.25 kb):\nA = %.0f, B = %.0f RPM"
        % (ma, mb), transform=ax.transAxes, fontsize=9, va="top",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85))
ax.legend(loc="upper right", frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig("fig2_rpm_normalized.png")
plt.close(fig)

# ---------------------------------------------------------------- fig 3
fig, ax = plt.subplots(figsize=(5.6, 5.6), dpi=150)
bg_vals = {}
for s, e, v in window_rows("A.sorted.bedgraph"):
    bg_vals[(s, e)] = v
pairs = []
with open("fig3_bw_window.txt") as f:
    for ln in f:
        s, e, v = ln.split()
        key = (int(s), int(e))
        if key in bg_vals:
            pairs.append((bg_vals[key], float(v)))
rt = json.load(open("parsed_roundtrip.json"))
px = [p[0] for p in pairs]
py = [p[1] for p in pairs]
lim = max(max(px), max(py)) * 1.06
ax.scatter(px, py, s=14, color=RED, alpha=0.55, edgecolors="none")
ax.plot([0, lim], [0, lim], color=GREEN, lw=1.4, ls="--")
ax.set_xlim(0, lim)
ax.set_ylim(0, lim)
ax.set_xlabel("bedGraph value (x), A.sorted.bedgraph")
ax.set_ylabel("value read back from A.bw (pyBigWig)")
ax.set_title("bigWig round-trip reproduces bedGraph values exactly")
ax.text(0.05, 0.93, "whole-chromosome check: n = %s rows,\nmax |delta| = %g"
        % ("{:,}".format(rt["roundtrip_rows"]), rt["roundtrip_max_abs_diff"]),
        transform=ax.transAxes, fontsize=9, va="top",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85))
fig.tight_layout()
fig.savefig("fig3_roundtrip_fidelity.png")
plt.close(fig)

print("figs done: fig1/fig2/fig3 written to", HERE)
