#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_figs.py - 3 figures from the real bedtools runs (skill 051 trial).
Colors: brick #b5482f (naive/enriched emphasis), celadon #2f7d72 (matched/honest emphasis).
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
BRICK = "#b5482f"
CELADON = "#2f7d72"
N_PERM = 1000

def read_rows(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        return [dict(zip(header, ln.rstrip("\n").split("\t"))) for ln in f if ln.strip()]

summary = read_rows(os.path.join(BASE, "summary.tsv"))
fisher = {r["set"]: r for r in read_rows(os.path.join(BASE, "fisher_summary.tsv"))}
S = {(r["set"], r["null_mode"]): r for r in summary}

def nulls(mode, sset):
    with open(os.path.join(BASE, "nulls_%s_%s.tsv" % (mode, sset))) as f:
        return [float(x) for x in f if x.strip()]

# ---------------- fig 1: matched-null jaccard distributions ----------------
fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), dpi=160)
panels = [("enriched", BRICK, "Enriched query set"),
          ("random", CELADON, "Random control set")]
for ax, (sset, color, label) in zip(axes, panels):
    vals = nulls("matched", sset)
    obs = float(S[(sset, "matched")]["obs_jaccard"])
    emp_p = float(S[(sset, "matched")]["emp_p"])
    ax.hist(vals, bins=30, color=color, alpha=0.55, edgecolor="white")
    ax.axvline(obs, color=color, lw=2, ls="--")
    ax.set_xlim(0.0, 0.42)
    ax.set_xlabel("Jaccard vs B (1000 workspace-matched shuffles)")
    ax.set_ylabel("Null replicates")
    ax.set_title(label)
    hits = sum(1 for v in vals if v >= obs)
    ax.text(0.97, 0.90, "observed = %.3f\nnull mean = %.3f\nempirical p = %.3f (%d/%d)"
            % (obs, sum(vals) / len(vals), emp_p, hits, N_PERM),
            transform=ax.transAxes, ha="right", va="top", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(BASE, "fig1_null_vs_observed.png"))
plt.close(fig)

# ---------------- fig 2: analytic fisher vs matched permutation ----------------
fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), dpi=160)
sets_order = ["enriched", "random"]
labels = ["Enriched\nquery", "Random\ncontrol"]
colors = [BRICK, CELADON]

ax = axes[0]
neglog = []
for sset in sets_order:
    lf = float(fisher[sset]["py_log10p"])
    neglog.append(-lf)
bars = ax.bar([0, 1], neglog, width=0.55, color=colors)
ax.set_xticks([0, 1]); ax.set_xticklabels(labels)
ax.set_ylabel("-log10 two-tailed p")
ax.set_ylim(0, max(neglog) * 1.32)
ax.set_title("bedtools fisher (whole-genome analytic null):\nboth sets flagged")
ax.text(0, neglog[0] + max(neglog) * 0.05, "p underflows to 0\n(log10 p = %.1f,\npython cross-check)" % neglog[0],
        ha="center", va="bottom", fontsize=8.5)
ax.text(1, neglog[1] + max(neglog) * 0.05, "p = 4.5e-26\nOR = 3.43",
        ha="center", va="bottom", fontsize=8.5)

ax = axes[1]
emp = [float(S[(s, "matched")]["emp_p"]) for s in sets_order]
ax.bar([0, 1], emp, width=0.55, color=colors)
ax.set_xticks([0, 1]); ax.set_xticklabels(labels)
ax.set_ylabel("Empirical permutation p (N = 1000)")
ax.set_ylim(0, 0.78)
ax.set_title("Workspace-matched permutation null:\nonly the enriched set is significant")
ax.text(0, emp[0] + 0.04, "p = 0.001\nfold = 4.6x", ha="center", va="bottom", fontsize=8.5)
ax.text(1, emp[1] + 0.04, "p = 0.591\nfold = 0.98x", ha="center", va="bottom", fontsize=8.5)
fig.tight_layout()
fig.savefig(os.path.join(BASE, "fig2_fisher_vs_permutation.png"))
plt.close(fig)

# ---------------- fig 3: the universe choice moves the answer ----------------
fig, ax = plt.subplots(figsize=(8.6, 4.6), dpi=160)
x = [0, 1]
matched_fold = [float(S[(s, "matched")]["fold_vs_null"]) for s in sets_order]
uniform_fold = [float(S[(s, "uniform")]["fold_vs_null"]) for s in sets_order]
matched_p = [float(S[(s, "matched")]["emp_p"]) for s in sets_order]
uniform_p = [float(S[(s, "uniform")]["emp_p"]) for s in sets_order]
bw = 0.32
ax.bar([i - bw / 2 for i in x], matched_fold, width=bw, color=CELADON,
       label="matched null (shuffle inside workspace)")
ax.bar([i + bw / 2 for i in x], uniform_fold, width=bw, color=BRICK,
       label="uniform null (shuffle over 2 Mb chr)")
ax.axhline(1.0, color="#555555", lw=1, ls=":")
ax.set_yscale("log")
ax.set_ylim(0.5, 60)
ax.set_xticks(x); ax.set_xticklabels(["Enriched query set", "Random control set"])
ax.set_ylabel("Fold enrichment over null (observed / null mean jaccard)")
ax.set_title("The universe choice moves the answer: uniform null inflates every set")
for i, (fm, fu) in enumerate(zip(matched_fold, uniform_fold)):
    ax.text(i - bw / 2, fm * 1.3, "matched\n%.2fx, p=%.3f" % (fm, matched_p[i]),
            ha="center", va="bottom", fontsize=8.5)
    ax.text(i + bw / 2, fu * 1.3, "uniform\n%.2fx, p=%.3f" % (fu, uniform_p[i]),
            ha="center", va="bottom", fontsize=8.5)
ax.legend(loc="upper center", framealpha=0.95, fontsize=8.5)
fig.tight_layout()
fig.savefig(os.path.join(BASE, "fig3_universe_matters.png"))
plt.close(fig)

print("figures written")
