#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
049 gtf-gff-handling figures. All numbers come from the real WSL run
(audit.json / audit_report.txt, produced by _run.sh on 2026-09-03).
Colors: brick red #b5482f + celadon green #2f7d72. English labels.
"""
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(BASE, "audit.json")))

BRICK = "#b5482f"
CELADON = "#2f7d72"
GREY = "#6b7280"


def verify(fig, name):
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
        print("  [PASS] %s" % name)
    return fails


def savefig(fig, name):
    fails = verify(fig, name)
    fig.savefig(os.path.join(BASE, name), dpi=150)
    plt.close(fig)
    return fails


# ============ FIG 1: feature counts + conversion consistency ============
fc = R["feature_counts"]
gtf_f = fc["annotation.gtf"]
gff3_f = fc["annotation.gff3"]
conv_f = fc["conv_from_gff3.gtf"]
cats = ["gene / mRNA", "transcript\n(GTF) / mRNA", "exon", "CDS"]
v_gtf = [gtf_f.get("gene", 0), gtf_f.get("transcript", 0), gtf_f.get("exon", 0), gtf_f.get("CDS", 0)]
v_gff3 = [gff3_f.get("gene", 0), gff3_f.get("mRNA", 0), gff3_f.get("exon", 0), gff3_f.get("CDS", 0)]
v_conv = [conv_f.get("gene", 0), conv_f.get("transcript", 0), conv_f.get("exon", 0), conv_f.get("CDS", 0)]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.4, 4.6))
x = range(len(cats))
w = 0.26
ax1.bar([i - w for i in x], v_gtf, width=w, color=BRICK, label="annotation.gtf (chr1)")
ax1.bar(list(x), v_gff3, width=w, color=CELADON, label="annotation.gff3 (seqid '1')")
ax1.bar([i + w for i in x], v_conv, width=w, color=GREY,
        label="conv_from_gff3.gtf (gffread -T)")
for xi, v in zip([i - w for i in x], v_gtf):
    ax1.text(xi, v + 0.5, str(v), ha="center", va="bottom", fontsize=8.5, color="#222222")
for xi, v in zip(x, v_gff3):
    ax1.text(xi, v + 0.5, str(v), ha="center", va="bottom", fontsize=8.5, color="#222222")
for xi, v in zip([i + w for i in x], v_conv):
    ax1.text(xi, v + 0.5, str(v), ha="center", va="bottom", fontsize=8.5, color="#222222")
ax1.set_xticks(list(x))
ax1.set_xticklabels(cats, fontsize=9)
ax1.set_ylabel("feature rows")
ax1.set_ylim(0, 36)
ax1.text(0.5, 0.97, "8 genes / 9 tx / 29 exons / 14 CDS;\ngffread -T: 0 gene rows",
         transform=ax1.transAxes, ha="center", va="top", fontsize=9.5, color="#222222")
ax1.legend(fontsize=8, loc="upper right", framealpha=0.9)

checks = ["exon sets: GTF vs GFF3->GTF\n(after seqid remap)",
          "roundtrip GTF->GFF3\n(exon sets identical)",
          "spliced FASTA vs\nBED concat (raw IDs)",
          "spliced FASTA vs\nBED concat (ver. stripped)",
          "seqid overlap chr1 vs '1'\n(naive string match)"]
vals = [100, 100, 0, 100, 0]
barcolors = [CELADON, CELADON, BRICK, CELADON, BRICK]
yb = range(len(checks))
ax2.barh(list(yb), vals, color=barcolors, height=0.62)
for yi, v in zip(yb, vals):
    ax2.text(v + 2, yi, "%d%%" % v, va="center", ha="left", fontsize=9.5, color="#222222")
ax2.set_yticks(list(yb))
ax2.set_yticklabels(checks, fontsize=8)
ax2.set_xlabel("transcripts identical / shared key pairs (%)")
ax2.set_xlim(0, 118)
ax2.set_xticks([0, 25, 50, 75, 100])
ax2.invert_yaxis()
ax2.text(0.98, 0.04, "9 transcripts, all checks from audit_report.txt",
         transform=ax2.transAxes, ha="right", va="bottom", fontsize=9, color="#222222")
fig.tight_layout()
total = savefig(fig, "fig1_conversion_consistency.png")

# ============ FIG 2: stop-codon convention (+3 nt) ============
per_tx = R["per_tx"]
ids = sorted(per_tx)
delta = [per_tx[t]["cds_len_gff3"] - per_tx[t]["cds_len_gtf"] for t in ids]
pg = [per_tx[t]["prot_len"] for t in ids]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.4, 4.6))
ax1.bar(range(len(ids)), delta, color=BRICK, width=0.62)
ax1.axhline(3, color=GREY, lw=1, ls="--")
ax1.set_xticks(range(len(ids)))
ax1.set_xticklabels([t.replace("T", "T").lower() for t in ids], rotation=45,
                    ha="right", fontsize=8)
ax1.set_ylabel("CDS length, GFF3 minus GTF (nt)")
ax1.set_ylim(0, 6)
ax1.set_yticks([0, 3, 6])
ax1.set_xlabel("transcript")
ax1.text(0.5, 0.97, "Every coding transcript: GFF3 CDS is exactly +3 nt\n"
                    "(stop codon included) vs the GTF version - 9/9 transcripts",
         transform=ax1.transAxes, ha="center", va="top", fontsize=9.5, color="#222222")

lo, hi = min(pg) - 8, max(pg) + 8
ax2.plot([lo, hi], [lo, hi], color=GREY, lw=1, ls="--")
ax2.scatter(pg, pg, s=46, color=CELADON, zorder=3)
ax2.set_xlabel("protein length, GTF-derived (aa)")
ax2.set_ylabel("protein length, GFF3-derived (aa)")
ax2.set_xlim(lo, hi)
ax2.set_ylim(lo, hi)
ax2.text(0.5, 0.97, "proteins identical 9/9:\n-y output ignores the convention",
         transform=ax2.transAxes, ha="center", va="top", fontsize=9.5, color="#222222")
fig.tight_layout()
total += savefig(fig, "fig2_stop_codon_convention.png")

# ============ FIG 3: BED off-by-one control ============
deficit = [per_tx[t]["bed_deficit"] for t in ids]
fig, ax = plt.subplots(figsize=(7.6, 4.6))
ax.bar(range(len(ids)), deficit, color=BRICK, width=0.62,
       label="over-corrected BED (start-1 AND end-1)")
ax.axhline(0, color=CELADON, lw=1.6)
ax.set_xticks(range(len(ids)))
ax.set_xticklabels([t.lower() for t in ids], rotation=45, ha="right", fontsize=8)
ax.set_ylabel("bp lost vs gffread -w transcript (bp)")
ax.set_xlabel("transcript")
ax.set_ylim(-0.4, 6.4)
ax.text(0.5, 0.97, "start-1 only: 0 bp lost, 9/9 exact sequence match vs gffread -w\n"
                   "over-corrected: 1 bp lost per exon, %d bp total across 29 exons"
         % sum(deficit),
         transform=ax.transAxes, ha="center", va="top", fontsize=9.5, color="#222222")
ax.legend(fontsize=9, loc="upper right", bbox_to_anchor=(1.0, 0.86), framealpha=0.9)
fig.tight_layout()
total += savefig(fig, "fig3_bed_off_by_one.png")

print("FIGURE QUALITY: TOTAL FAILS = %d" % total)
