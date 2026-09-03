#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
052 proximity-operations parse_results.py
Compare real bedtools outputs against truth.json; write results.json.
All verification numbers reported here are real measured values.
"""
import json
import os

BASE = os.path.dirname(os.path.abspath(__file__))
W = 50000


def read_bed(path, ncols=None):
    rows = []
    with open(os.path.join(BASE, path)) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if ncols and len(p) != ncols:
                raise SystemExit("%s: expected %d cols got %d: %r" % (path, ncols, len(p), p))
            rows.append(p)
    return rows


truth = json.load(open(os.path.join(BASE, "truth.json")))
tn, tt, tw = truth["nearest"], truth["ties"], truth["window"]
S = dict(truth["summary"])

res = {"summary": S}

# ---------- closest -D b -io -t first (12 cols: A5 + B6 + dist) ----------
rows = read_bed("nearest_db_io_first.bed", 12)
n_ok_gene = n_ok_dist = n_sentinel = n_sign_ok = 0
signed_err = []
for r in rows:
    pk = r[3]
    gene_col, dist_col = r[8], r[11]
    t = tn[pk]
    if gene_col == ".":                       # sentinel row: no B on chrom
        n_sentinel += 1
        if t is None and dist_col == "-1":
            n_ok_gene += 1                    # sentinel itself is the correct call
        continue
    if t is None:
        continue
    if gene_col == t["gene"]:
        n_ok_gene += 1
    d = int(dist_col)
    if d == t["dist_signed_db"]:
        n_ok_dist += 1
    else:
        signed_err.append({"peak": pk, "bedtools": d, "truth": t["dist_signed_db"]})
    if (d < 0) == (t["dist_signed_db"] < 0):
        n_sign_ok += 1

res["closest_db"] = {
    "rows": len(rows),
    "gene_assign_correct": n_ok_gene,
    "assign_accuracy_pct": round(100.0 * n_ok_gene / S["n_peaks"], 2),
    "signed_dist_correct": n_ok_dist,
    "signed_dist_errors": signed_err[:10],
    "n_signed_dist_checked": n_ok_dist + len(signed_err),
    "sign_correct": n_sign_ok,
    "none_sentinel_rows": n_sentinel,
    "none_sentinel_expected": S["n_peaks_no_gene"],
}

# ---------- closest -d -t all (tie rows) ----------
rows_all = read_bed("nearest_all_d.bed", 12)
per_peak_rows = {}
for r in rows_all:
    per_peak_rows[r[3]] = per_peak_rows.get(r[3], 0) + 1
exp_rows = S["n_peaks"] + sum(v - 1 for v in tt.values() if v > 1)
# per-peak expected row count: ties emit `ties` rows; peaks whose chrom has no
# B still emit ONE sentinel row (`none` / -1), so expected = max(ties, 1)
tie_ok = all(per_peak_rows.get(k, 0) == max(v, 1) for k, v in tt.items())
res["closest_tall"] = {
    "rows": len(rows_all), "expected_rows": exp_rows,
    "rows_match_truth": len(rows_all) == exp_rows,
    "tied_peaks_detected": sum(1 for v in per_peak_rows.values() if v > 1),
    "tied_peaks_expected": S["n_tied_peaks"],
    "per_peak_tie_counts_match_truth": tie_ok,
    "inflated_rows": len(rows_all) - S["n_peaks"],
}

# ---------- closest -D ref (mis-sign minus strand) ----------
rows_ref = read_bed("nearest_ref_io_first.bed", 12)
ref_map = {r[3]: int(r[11]) for r in rows_ref if r[8] != "." and int(r[11]) != -1}
db_map = {r[3]: int(r[11]) for r in rows if r[8] != "." and int(r[11]) != -1}
mismatch = [pk for pk in db_map if (ref_map[pk] < 0) != (db_map[pk] < 0)]
res["closest_ref_vs_db"] = {
    "nonzero_calls": len(db_map),
    "sign_mismatch_peaks": len(mismatch),
    "expected_mismatch_minus_strand": S["ref_vs_db_mismatch_expected"],
    "match_expected": len(mismatch) == S["ref_vs_db_mismatch_expected"],
}

# ---------- closest -k 3 -d ----------
# NOTE: this run has no -io, so a B overlapping A is reported first with
# distance 0. Truth nearest excludes overlapping B (-io semantics), so the
# first row is expected to differ from truth exactly for peaks that overlap
# a gene body.
genes_raw = []
with open(os.path.join(BASE, "genes.bed")) as f:
    for line in f:
        p = line.split("\t")
        genes_raw.append((p[0], int(p[1]), int(p[2])))
rows_k3 = read_bed("top3_k3_d.bed", 12)
first_rows = {}
for r in rows_k3:
    if r[8] != ".":
        first_rows.setdefault(r[3], r)
peak_chrom = {}
for r in rows_k3:
    peak_chrom[r[3]] = (int(r[1]), int(r[2]), r[0])
k3_nearest_ok = n_overlap_peaks = 0
for pk, r in first_rows.items():
    t = tn[pk]
    if t is None:
        continue
    ps, pe, pc = peak_chrom[pk]
    overlaps = any(pc == gc and ps < ge and gs < pe for gc, gs, ge in genes_raw)
    if overlaps:
        n_overlap_peaks += 1
    elif r[8] == t["gene"]:
        k3_nearest_ok += 1
res["closest_k3"] = {
    "rows": len(rows_k3), "peaks": len(first_rows),
    "peaks_overlapping_a_gene": n_overlap_peaks,
    "first_row_gene_matches_truth_nonoverlap": k3_nearest_ok,
    "first_row_expected_match": len(first_rows) - n_overlap_peaks,
}

# ---------- window -w 50000 -c (6 cols: A5 + count) ----------
rows_w = read_bed("window_counts.bed", 6)
n_w_ok = 0
measured_counts = []
for r in rows_w:
    c = int(r[5])
    measured_counts.append(c)
    if c == tw[r[3]]:
        n_w_ok += 1
res["window"] = {
    "rows": len(rows_w),
    "count_match_peaks": n_w_ok,
    "count_match_pct": round(100.0 * n_w_ok / len(rows_w), 2),
    "total_gene_hits_measured": sum(measured_counts),
    "total_gene_hits_expected": sum(tw.values()),
    "mean_genes_per_peak": round(sum(measured_counts) / len(measured_counts), 2),
    "max_genes_per_peak": max(measured_counts),
    "peaks_zero_hits": sum(1 for c in measured_counts if c == 0),
    "per_peak_counts": measured_counts,
}

# ---------- promoters (slop -s -l 2000 -r 200 on TSS) ----------
rows_p = read_bed("promoters.bed", 6)
widths = [int(r[2]) - int(r[1]) for r in rows_p]
prom_bad = []
for r, wmeas in zip(rows_p, widths):
    exp = truth["promoters"][r[3]]["exp"]
    if [int(r[1]), int(r[2])] != exp:
        prom_bad.append({"gene": r[3], "bedtools": [int(r[1]), int(r[2])], "truth": exp})
res["promoters"] = {
    "n": len(rows_p), "n_clipped": sum(1 for w in widths if w != 2201),
    "n_clipped_expected": S["n_promoters_clipped"],
    "width_2201": sum(1 for w in widths if w == 2201),
    "intervals_match_truth": len(prom_bad) == 0,
    "mismatches": prom_bad[:5],
    "clipped_widths": sorted(w for w in widths if w != 2200),
}

# ---------- gene-body slop contrast ----------
rows_gb = read_bed("genebody_slop_b2000.bed", 6)
gb_widths = [int(r[2]) - int(r[1]) for r in rows_gb]
res["genebody_slop"] = {
    "n": len(rows_gb),
    "mean_width": round(sum(gb_widths) / len(gb_widths), 1),
    "max_width": max(gb_widths),
    "min_width": min(gb_widths),
    "promoter_mean_width": round(sum(widths) / len(widths), 1),
}

# ---------- flank ----------
rows_f = read_bed("gene_flanks.bed", 6)
res["flank"] = {
    "regions": len(rows_f),
    "expected_if_no_clip": S["expected_flanks_if_no_clip"],
    "expected_with_clip": S["expected_flank_regions"],
    "dropped": S["expected_flanks_if_no_clip"] - len(rows_f),
}

with open(os.path.join(BASE, "results.json"), "w") as f:
    json.dump(res, f, indent=1)
print(json.dumps(res, indent=1)[:6000])
print("parse_results.py done")
