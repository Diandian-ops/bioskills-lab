#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
052 proximity-operations make_inputs.py
Simulated dataset for bedtools closest / window / slop / flank truth-checking.

Design (seed 20260903):
  chr1 2,000,000 bp ; chr2 1,200,000 bp ; chr3 600,000 bp (peaks only, NO genes
  -> exercises the `none` / -1 sentinel of closest).
  Genes: BED6 (chrom start end name score strand). Hand-placed edge cases:
    - chr1 + strand gene with TSS 500 bp from contig start (slop clipping at 0)
    - chr1 head-to-head bidirectional pair with a peak exactly midway (closest tie)
    - chr1 peak fully inside a gene body (-io exclusion case)
    - chr2 - strand gene with TSS near contig end (slop clipping at chrom length)
    - chr2 + strand gene starting at 0 (flank drops the left flank; promoter clipped)
  Truth values (nearest non-overlapping gene with -D b sign, tie count, window
  hit counts at w=50 kb, expected promoter intervals, expected flank counts)
  are computed here and stored in truth.json.
"""
import json
import os
import random

random.seed(20260903)

BASE = os.path.dirname(os.path.abspath(__file__))
CHROMS = {"chr1": 2000000, "chr2": 1200000, "chr3": 600000}
W_WINDOW = 50000          # window -w
PROM_UP, PROM_DOWN = 2000, 200   # slop -s -l 2000 -r 200
FLANK_B = 1000            # flank -b

genes = []  # dicts: chrom start end name score strand
gid = [0]


def add_gene(chrom, start, end, strand):
    gid[0] += 1
    genes.append({"chrom": chrom, "start": start, "end": end,
                  "name": "g%03d" % gid[0], "score": 500, "strand": strand})


def overlaps_any(chrom, start, end):
    for g in genes:
        if g["chrom"] == chrom and start < g["end"] and g["start"] < end:
            return True
    return False


# ---- hand-placed edge-case genes ----
add_gene("chr1", 500, 4500, "+")          # g001 TSS 500 bp from contig start
add_gene("chr1", 100000, 104000, "+")     # g002 bidirectional pair left (+)
add_gene("chr1", 106000, 110000, "-")     # g003 bidirectional pair right (-)
add_gene("chr1", 195000, 205000, "+")     # g004 will host the -io demo peak
add_gene("chr1", 300000, 342000, "-")     # g005 ordinary minus-strand gene
add_gene("chr2", 0, 2500, "+")            # g006 starts at 0 (flank/promoter clipping)
add_gene("chr2", 498800, 556000, "-")     # g007 bidirectional pair right (-)
add_gene("chr2", 492000, 496000, "+")     # g008 bidirectional pair left (+)
add_gene("chr2", 1198500, 1199000, "-")   # g009 TSS near contig end (clipping)

# ---- random genes, non-overlapping (min gap 3000 bp) ----
for chrom, n in (("chr1", 30), ("chr2", 20)):
    placed, tries = 0, 0
    while placed < n and tries < 20000:
        tries += 1
        L = random.randint(2000, 60000)
        s = random.randint(1000, CHROMS[chrom] - L - 1000)
        e = s + L
        if overlaps_any(chrom, s - 3000, e + 3000):
            continue
        add_gene(chrom, s, e, random.choice("+-"))
        placed += 1
    assert placed == n, "gene placement failed on %s" % chrom

# ---- peaks ----
peaks = []  # chrom start end name score


def add_peak(chrom, start, end, name):
    peaks.append({"chrom": chrom, "start": start, "end": end,
                  "name": name, "score": 500})


# designed peaks
add_peak("chr1", 104800, 105200, "peak_tie1")    # equidistant to g002/g003 (800 bp)
add_peak("chr2", 497200, 497600, "peak_tie2")    # midpoint of g008/g007
add_peak("chr1", 200000, 200400, "peak_io1")     # fully inside g004 -> -io demo
add_peak("chr1", 300, 800, "peak_edge1")         # near contig start, overlaps g001
add_peak("chr2", 0, 900, "peak_edge2")           # overlaps g006 at contig start

# promoter-proximal peaks: 10 peaks within +/-800 bp of a random gene TSS
prox_pool = [g for g in genes if g["name"] not in ("g001", "g006")]
random.shuffle(prox_pool)
for g in prox_pool[:10]:
    tss = g["start"] if g["strand"] == "+" else g["end"] - 1
    s = tss + random.randint(-800, 800)
    add_peak(g["chrom"], s, s + random.randint(150, 500),
             "peak_prox%02d" % (len([p for p in peaks if p["name"].startswith("peak_prox")]) + 1))

# distal peaks: 45 random, chr1/chr2 only, no fixed relation to genes
for i in range(45):
    chrom = random.choice(["chr1", "chr2"])
    s = random.randint(1000, CHROMS[chrom] - 1000)
    add_peak(chrom, s, s + random.randint(200, 800), "peak_dist%02d" % (i + 1))

# chr3 peaks: 12, on a chromosome with no genes -> closest sentinel `none` / -1
for i in range(12):
    s = random.randint(1000, CHROMS["chr3"] - 1000)
    add_peak("chr3", s, s + random.randint(200, 600), "peak_chr3_%02d" % (i + 1))

peaks.sort(key=lambda p: (p["chrom"], p["start"]))
genes.sort(key=lambda g: (g["chrom"], g["start"]))

# ================= truth computation =================
def gap(a_s, a_e, b_s, b_e):
    """signed-free bp distance between intervals; 0 if overlapping"""
    if a_s < b_e and b_s < a_e:
        return 0
    return a_s - b_e if a_s >= b_e else b_s - a_e


def signed_db(peak, g):
    """closest -D b convention: sign by gene strand.
    negative = peak upstream of gene (gene orientation), positive = downstream,
    0 = overlap (excluded by -io anyway).
    bedtools -d distance for non-overlapping features = bp gap + 1 (half-open
    BED ends + inclusive-base convention); verified against real run output."""
    gp = gap(peak["start"], peak["end"], g["start"], g["end"])
    if gp == 0:
        return 0
    left = peak["end"] <= g["start"]           # peak lies to lower coordinates
    upstream = (g["strand"] == "+" and left) or (g["strand"] == "-" and not left)
    return (-1 if upstream else 1) * (gp + 1)


truth_nearest = {}     # peak name -> dict(gene, dist_unsigned, dist_signed_db, strand)
truth_ties = {}        # peak name -> number of tied nearest non-overlapping genes
truth_window = {}      # peak name -> number of genes within +/-50 kb (or overlapping)
for p in peaks:
    cands = [g for g in genes if g["chrom"] == p["chrom"] and not (
        p["start"] < g["end"] and g["start"] < p["end"])]        # -io semantics
    if cands:
        dists = [abs(gap(p["start"], p["end"], g["start"], g["end"])) for g in cands]
        dmin = min(dists)
        tied = [g for g, d in zip(cands, dists) if d == dmin]
        gb = tied[0]
        truth_nearest[p["name"]] = {
            "gene": gb["name"], "strand": gb["strand"],
            "dist_unsigned": dmin, "dist_signed_db": signed_db(p, gb),
        }
        truth_ties[p["name"]] = len(tied)
    else:
        truth_nearest[p["name"]] = None
        truth_ties[p["name"]] = 0
    # window: overlap or within W_WINDOW (all genes, no -io concept in window)
    hits = [g for g in genes if g["chrom"] == p["chrom"]
            and p["start"] < g["end"] + W_WINDOW and g["start"] < p["end"] + W_WINDOW]
    truth_window[p["name"]] = len(hits)

# TSS + expected promoters (slop -s -l 2000 -r 200), incl. clipping.
# bedtools slop semantics: start -= l, end += r, then -s swaps sides for minus
# strand; TSS is a 1-bp interval so an unclipped promoter is 2201 bp wide.
tss, prom_expected = [], {}
for g in genes:
    if g["strand"] == "+":
        t0, t1 = g["start"], g["start"] + 1
        exp = [max(0, t0 - PROM_UP), min(CHROMS[g["chrom"]], t1 + PROM_DOWN)]
    else:
        t0, t1 = g["end"] - 1, g["end"]
        exp = [max(0, t0 - PROM_DOWN), min(CHROMS[g["chrom"]], t1 + PROM_UP)]
    tss.append({"chrom": g["chrom"], "start": t0, "end": t1,
                "name": g["name"], "score": g["score"], "strand": g["strand"]})
    prom_expected[g["name"]] = {"exp": exp, "width": exp[1] - exp[0], "clipped": None}
    if exp[0] == 0 and t0 - (PROM_UP if g["strand"] == "+" else PROM_DOWN) < 0:
        prom_expected[g["name"]]["clipped"] = "start"
    if exp[1] == CHROMS[g["chrom"]] and t1 + (PROM_DOWN if g["strand"] == "+" else PROM_UP) > CHROMS[g["chrom"]]:
        prom_expected[g["name"]]["clipped"] = "end"

# expected flank counts: flank drops only a ZERO-LENGTH side (a clipped-but-
# non-empty flank, e.g. start 500 -> [0,500], is still emitted); a side dies
# only when start == 0 (left) or end == chrom length (right)
expected_flanks = 0
for g in genes:
    if g["start"] > 0:
        expected_flanks += 1
    if g["end"] < CHROMS[g["chrom"]]:
        expected_flanks += 1

# -D ref vs -D b: for every nonzero-distance nearest call, ref sign = coordinate
# sign (peak left of gene -> negative), so it differs from -D b exactly when the
# gene is on the minus strand
ref_vs_db_mismatch = sum(
    1 for v in truth_nearest.values()
    if v is not None and v["dist_unsigned"] > 0 and v["strand"] == "-")

summary = {
    "n_chroms": len(CHROMS), "chrom_sizes": CHROMS,
    "n_genes": len(genes),
    "n_genes_chr1": sum(1 for g in genes if g["chrom"] == "chr1"),
    "n_genes_chr2": sum(1 for g in genes if g["chrom"] == "chr2"),
    "n_peaks": len(peaks),
    "n_peaks_chr1": sum(1 for p in peaks if p["chrom"] == "chr1"),
    "n_peaks_chr2": sum(1 for p in peaks if p["chrom"] == "chr2"),
    "n_peaks_chr3": sum(1 for p in peaks if p["chrom"] == "chr3"),
    "n_peaks_with_gene": sum(1 for v in truth_nearest.values() if v is not None),
    "n_peaks_no_gene": sum(1 for v in truth_nearest.values() if v is None),
    "n_tied_peaks": sum(1 for v in truth_ties.values() if v > 1),
    "n_minus_strand_genes": sum(1 for g in genes if g["strand"] == "-"),
    "expected_flank_regions": expected_flanks,
    "expected_flanks_if_no_clip": 2 * len(genes),
    "ref_vs_db_mismatch_expected": ref_vs_db_mismatch,
    "n_promoters_clipped": sum(1 for v in prom_expected.values() if v["clipped"]),
    "w_window": W_WINDOW, "prom_up": PROM_UP, "prom_down": PROM_DOWN,
    "flank_b": FLANK_B, "seed": 20260903,
}

# ================= write files =================
with open(os.path.join(BASE, "genes.bed"), "w") as f:
    for g in genes:
        f.write("%s\t%d\t%d\t%s\t%d\t%s\n" % (g["chrom"], g["start"], g["end"],
                                              g["name"], g["score"], g["strand"]))
with open(os.path.join(BASE, "peaks.bed"), "w") as f:
    for p in peaks:
        f.write("%s\t%d\t%d\t%s\t%d\n" % (p["chrom"], p["start"], p["end"],
                                          p["name"], p["score"]))
with open(os.path.join(BASE, "genome.txt"), "w") as f:
    for c, L in CHROMS.items():
        f.write("%s\t%d\n" % (c, L))
with open(os.path.join(BASE, "truth.json"), "w") as f:
    json.dump({"nearest": truth_nearest, "ties": truth_ties,
               "window": truth_window, "promoters": prom_expected,
               "summary": summary}, f, indent=1)

print(json.dumps(summary, indent=2))
print("make_inputs.py: wrote genes.bed (%d), peaks.bed (%d), genome.txt, truth.json"
      % (len(genes), len(peaks)))
