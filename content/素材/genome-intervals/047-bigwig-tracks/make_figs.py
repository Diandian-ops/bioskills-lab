#!/usr/bin/env python3
"""Make the three figures for 047 bigwig-tracks from the real run outputs.

Reads (same directory as this script):
  sim.sorted.bedGraph   -- simulated input track (158,000 bins)
  regions_avg.tab       -- bigWigAverageOverBed output (real tool output)
  truth.json            -- design ground truth
  _bigwiginfo.txt       -- bigWigInfo header output (real tool output)
Writes: fig1_signal_track.png, fig2_region_reconciliation.png, fig3_dilution.png
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
RED = "#b5482f"
GREEN = "#2f7d72"
GRAY = "#8a8a8a"

fig1_path = os.path.join(BASE, "fig1_signal_track.png")
fig2_path = os.path.join(BASE, "fig2_region_reconciliation.png")
fig3_path = os.path.join(BASE, "fig3_dilution.png")


def load_track():
    xs, vs = [], []
    with open(os.path.join(BASE, "sim.sorted.bedGraph")) as f:
        for line in f:
            c, s, e, v = line.split("\t")
            xs.append((int(s), int(e)))
            vs.append(float(v))
    return xs, vs


def load_avg():
    rows = []
    with open(os.path.join(BASE, "regions_avg.tab")) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) == 6:
                rows.append({
                    "name": cols[0], "size": int(cols[1]),
                    "covered": int(cols[2]), "sum": float(cols[3]),
                    "mean0": float(cols[4]), "mean": float(cols[5]),
                })
    return rows


def load_max_from_info():
    mx = None
    with open(os.path.join(BASE, "_bigwiginfo.txt")) as f:
        for line in f:
            if line.startswith("max:"):
                mx = float(line.split(":")[1].strip().replace(",", ""))
    return mx


def make_fig1():
    """Chromosome-wide signal profile (200 bp max-binned for display)."""
    xs, vs = load_track()
    # 200 bp display bins: max value per 20 consecutive 10 bp bins
    disp_x, disp_v = [], []
    cur_start, cur_max, n = None, 0.0, 0
    for (s, e), v in zip(xs, vs):
        b = s // 200
        if cur_start is None or b != cur_start:
            if cur_start is not None:
                disp_x.append(cur_start * 200 + 100)
                disp_v.append(cur_max)
            cur_start, cur_max, n = b, v, 1
        else:
            cur_max = max(cur_max, v)
    disp_x.append(cur_start * 200 + 100)
    disp_v.append(cur_max)

    fig, ax = plt.subplots(figsize=(10.0, 3.6))
    ax.axvspan(1.5e6, 1.6e6, color=GRAY, alpha=0.25, lw=0, label="gap (no data)")
    ax.axvspan(1.6e6, 2.0e6, color=GREEN, alpha=0.15, lw=0,
               label="sparse zone (20% of bins)")
    ax.plot(disp_x, disp_v, color=RED, lw=0.8,
            label="signal (200 bp max-binned)")
    ax.set_yscale("log")
    ax.set_ylim(0.05, 3000)
    ax.set_xlim(0, 2.0e6)
    ax.set_xlabel("Position on chrSim (bp)")
    ax.set_ylabel("Signal value (log scale)")
    ax.set_title("Simulated signal track chrSim: background + 6 peaks + gap/sparse zones")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False,
              fontsize=8)
    fig.tight_layout()
    fig.savefig(fig1_path, dpi=150)
    plt.close(fig)


def make_fig2():
    """bigWigAverageOverBed mean vs mean0 per region, with design truth."""
    rows = load_avg()
    truth = json.load(open(os.path.join(BASE, "truth.json")))
    names = [r["name"] for r in rows]
    mean0 = [r["mean0"] for r in rows]
    mean = [r["mean"] for r in rows]
    truth_mean = [truth[n]["mean"] if truth[n]["mean"] is not None else 0.0
                  for n in names]

    fig, ax = plt.subplots(figsize=(9.6, 4.4))
    x = range(len(names))
    w = 0.38
    ax.bar([i - w / 2 for i in x], mean0, width=w, color=GREEN,
           label="mean0 (gaps as zero)")
    ax.bar([i + w / 2 for i in x], mean, width=w, color=RED,
           label="mean (covered bases only)")
    ax.scatter(list(x), truth_mean, marker="D", s=28, color="black", zorder=5,
               label="design truth (mean)")
    ax.set_yscale("symlog", linthresh=0.1)
    ax.set_ylim(0, 900)
    ax.set_xticks(list(x))
    ax.set_xticklabels(names, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Per-region mean signal (symlog scale)")
    ax.set_title("bigWigAverageOverBed output vs design truth (all regions PASS)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False,
              fontsize=8)
    fig.tight_layout()
    fig.savefig(fig2_path, dpi=150)
    plt.close(fig)


def make_fig3():
    """Wide mean dilutes a narrow tall peak: same file, different windows."""
    rows = load_avg()
    avg = {r["name"]: r for r in rows}
    peak_max = load_max_from_info()
    vals = [
        ("Peak max\n(bigWigInfo)", peak_max, GREEN),
        ("mean over\n400 bp window", avg["peak5_400bp"]["mean"], RED),
        ("mean over\n1 Mb window", avg["wide_1Mb_peak5"]["mean"], RED),
        ("mean over\n50 kb background", avg["bg_50kb"]["mean"], GREEN),
    ]
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    xs = range(len(vals))
    ax.bar(xs, [v for _, v, _ in vals], color=[c for _, _, c in vals], width=0.6)
    ax.set_yscale("log")
    ax.set_ylim(0.1, 30000)
    ax.set_xticks(list(xs))
    ax.set_xticklabels([n for n, _, _ in vals], fontsize=8.5)
    ax.set_ylabel("Signal (log scale)")
    ax.set_title("Peak 5 (amplitude 800): window width decides the headline number")
    for i, (_, v, _) in enumerate(vals):
        ax.text(i, v * 1.6, "%.5g" % v, ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    fig.savefig(fig3_path, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    make_fig1()
    make_fig2()
    make_fig3()
    print("figures written to", BASE)
