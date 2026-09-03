#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""050 reconcile: parse every bedtools output produced by _run.sh, compare with
truth.json (independent pure-python interval algebra), write results.json/txt."""
import json
import os
import re

BASE = os.path.dirname(os.path.abspath(__file__))
os.chdir(BASE)

T = json.load(open("truth.json"))


def wc(fn):
    with open(fn) as f:
        return sum(1 for _ in f)


def rd(fn):
    with open(fn) as f:
        return f.read()


def cols(fn):
    with open(fn) as f:
        for ln in f:
            yield ln.rstrip("\n").split("\t")


checks = []
M = {}


def chk(name, measured, expected, tol=None):
    if tol is None:
        ok = (measured == expected)
    else:
        ok = (measured is not None and expected is not None
              and abs(measured - expected) <= tol)
    checks.append(dict(name=name, measured=measured, expected=expected, ok=bool(ok)))
    return ok


# ---------- intersect output modes ----------
M["u"] = wc("i_u.bed")
M["v"] = wc("i_v.bed")
M["pairs"] = wc("i_wawb.bed")
M["loj_lines"] = wc("i_loj.bed")
M["wao_lines"] = wc("i_wao.bed")
M["wo_bp"] = sum(int(r[-1]) for r in cols("i_wo.bed"))
c_counts = [int(r[-1]) for r in cols("i_c.bed")]
M["c_lines"] = len(c_counts)
M["c_pos"] = sum(1 for x in c_counts if x > 0)
M["c_total_hits"] = sum(c_counts)

chk("intersect -u (A kept)", M["u"], T["intersect"]["u"])
chk("intersect -v (A dropped)", M["v"], T["intersect"]["v"])
chk("intersect -wa -wb pairs", M["pairs"], T["intersect"]["pairs"])
chk("intersect -wo overlap bp", M["wo_bp"], T["intersect"]["wo_bp"])
chk("intersect -c lines", M["c_lines"], T["counts"]["peaks"])
chk("intersect -c A with count>0", M["c_pos"], T["intersect"]["u"])
chk("intersect -c total hits = pairs", M["c_total_hits"], T["intersect"]["pairs"])
chk("intersect -loj lines = pairs + no-hit A", M["loj_lines"],
    T["intersect"]["loj_lines"])
chk("intersect -wao lines = pairs + no-hit A", M["wao_lines"],
    T["intersect"]["wao_lines"])

# ---------- fractional overlap ----------
M["f05_u"] = wc("i_f05.bed")
M["f05r_u"] = wc("i_f05r.bed")
M["f05_u_swapped"] = wc("i_f05_swap.bed")
chk("-f 0.5 -u (A=peaks)", M["f05_u"], T["intersect"]["f05_u"])
chk("-f 0.5 -r -u (reciprocal)", M["f05r_u"], T["intersect"]["f05r_u"])
chk("-f 0.5 -u (A=genes, roles swapped)", M["f05_u_swapped"],
    T["intersect"]["f05_u_swapped"])

# ---------- subtract ----------
M["sub_residual_bp"] = sum(int(r[2]) - int(r[1]) for r in cols("sub_plain.bed"))
M["sub_residual_ivs"] = wc("sub_plain.bed")
M["subA_remaining"] = wc("sub_A.bed")
chk("subtract residual bp", M["sub_residual_bp"], T["subtract"]["residual_bp"])
chk("subtract residual interval count", M["sub_residual_ivs"],
    T["subtract"]["residual_ivs"])
chk("subtract -A remaining", M["subA_remaining"], T["subtract"]["minusA_remaining"])

# ---------- merge ----------
M["merge_unsorted_exit"] = int(rd("merge_unsorted.exit").strip())
M["merge_unsorted_lines"] = wc("merge_unsorted.bed")
_mue = rd("merge_unsorted.err").strip().splitlines()
M["merge_unsorted_msg"] = _mue[0] if _mue else ""
chk("merge on UNSORTED file refused (exit 1)", M["merge_unsorted_exit"], 1)
chk("merge on UNSORTED file output empty", M["merge_unsorted_lines"], 0)
M["merge_d0"] = wc("merge_d0.bed")
M["merge_d1"] = wc("merge_d1.bed")
M["merge_d100"] = wc("merge_d100.bed")
M["merge_d0_bp"] = sum(int(r[2]) - int(r[1]) for r in cols("merge_d0.bed"))
chk("merge -d 0 sorted", M["merge_d0"], T["merge"]["sorted_d0"])
chk("merge -d 1", M["merge_d1"], T["merge"]["sorted_d1"])
chk("merge -d 100", M["merge_d100"], T["merge"]["sorted_d100"])
chk("merged coverage bp (d0)", M["merge_d0_bp"], T["merge"]["sorted_d0_bp"])
score_sum_merged = sum(int(r[-1]) for r in cols("merge_co.bed"))
score_sum_raw = sum(int(r[4]) for r in cols("peaks.bed"))
chk("merge -c 5 -o sum preserves score total", score_sum_merged, score_sum_raw)

# ---------- complement / cluster ----------
M["comp_n"] = wc("comp.bed")
M["comp_bp"] = sum(int(r[2]) - int(r[1]) for r in cols("comp.bed"))
chk("complement interval count", M["comp_n"], T["complement"]["n"])
chk("complement bp", M["comp_bp"], T["complement"]["bp"])
M["cluster_n"] = max(int(r[-1]) for r in cols("cluster.bed"))
chk("cluster count = merge d0 blocks", M["cluster_n"], T["cluster"]["n"])

# ---------- map / groupby ----------
map_meas = {}
for r in cols("map_mean.bed"):
    map_meas[r[3]] = None if r[6] == "." else float(r[6])
map_exp = T["map_mean"]
# bedtools map prints means with 8 decimal places, so compare at 1e-6
bad = [g for g in map_exp
       if (map_meas.get(g) is None) != (map_exp[g] is None)
       or (map_meas.get(g) is not None
           and abs(map_meas[g] - map_exp[g]) > 1e-6)]
M["map_genes"] = len(map_meas)
M["map_mismatch"] = len(bad)
chk("map mean per gene mismatches", len(bad), 0)

gb_meas = {r[3]: int(r[4]) for r in cols("groupby.tsv")}
gb_exp = T["groupby_sum_bp"]
bad_gb = [g for g in gb_exp if gb_meas.get(g, 0) != gb_exp[g]]
M["groupby_genes"] = len(gb_meas)
M["groupby_mismatch"] = len(bad_gb)
chk("groupby per-gene overlap bp mismatches", len(bad_gb), 0)

# ---------- multiinter / unionbedg ----------
M["multi_n3"] = sum(1 for r in cols("multiinter.tsv") if r[3] == "3"
                    and not r[0].startswith("chrom"))
M["multi_bp3"] = sum(int(r[2]) - int(r[1]) for r in cols("multiinter.tsv")
                     if r[3] == "3" and not r[0].startswith("chrom"))
chk("multiinter intervals covered by all 3", M["multi_n3"], T["multiinter"]["n3"])
chk("multiinter all-3 bp", M["multi_bp3"], T["multiinter"]["bp3"])
M["unionbedg_rows"] = sum(1 for r in cols("unionbedg.tsv") if r[0] != "chrom")
chk("unionbedg rows", M["unionbedg_rows"], T["unionbedg"]["rows"])

# ---------- jaccard / fisher ----------
jl = [l for l in rd("jac.txt").strip().splitlines() if l.strip()]
hdr = jl[0].split("\t")
vals = jl[1].split("\t")
jac = dict(zip(hdr, vals))
M["jaccard"] = float(jac["jaccard"])
M["jaccard_inter_bp"] = int(jac["intersection"])
M["jaccard_union_bp"] = int(jac["union"])
M["jaccard_n_inter"] = int(jac["n_intersections"])
chk("jaccard value (6 dp)", M["jaccard"], round(T["jaccard"]["jaccard"], 6),
    tol=1e-5)
chk("jaccard intersection bp", M["jaccard_inter_bp"], T["jaccard"]["inter_bp"])
chk("jaccard union bp", M["jaccard_union_bp"], T["jaccard"]["union_bp"])
chk("jaccard n_intersections", M["jaccard_n_inter"],
    T["jaccard"]["n_intersections"])
fisher_text = rd("fisher.txt")

# ---------- -split ----------
M["split_env_u"] = wc("spl_env_u.bed")
M["split_block_u"] = wc("spl_blk_u.bed")
M["split_env_bp"] = sum(int(r[-1]) for r in cols("spl_env.wo"))
M["split_block_bp"] = sum(int(r[-1]) for r in cols("spl_blk.wo"))
chk("BED12 envelope -u", M["split_env_u"], T["split"]["env_u"])
chk("BED12 -split -u", M["split_block_u"], T["split"]["block_u"])
chk("BED12 envelope overlap bp", M["split_env_bp"], T["split"]["env_bp"])
chk("BED12 -split overlap bp", M["split_block_bp"], T["split"]["block_bp"])

# ---------- sorted contract / chrom naming ----------
M["sorted_inmem_u"] = wc("i_u_sorted.bed")
M["sorted_sweep_u"] = wc("i_u_sorted_sweep.bed")
chk("intersect -sorted -g equals in-memory result", M["sorted_sweep_u"],
    M["sorted_inmem_u"])
M["mismatch_u"] = wc("mism_u.bed")
chk("chr1-vs-1 silent empty result", M["mismatch_u"], 0)
err_unsorted = rd("err_unsorted.log").strip().splitlines()
M["err_unsorted_exit"] = int(rd("err_unsorted.exit").strip())
M["err_unsorted_msg"] = err_unsorted[-1] if err_unsorted else ""
err_rev = rd("err_revorder.log").strip().splitlines()
M["err_revorder_exit"] = int(rd("err_revorder.exit").strip())
M["err_revorder_msg"] = err_rev[-1] if err_rev else ""

# ---------- perf / pybedtools ----------
perf_lines = [l for l in rd("perf_summary.txt").strip().splitlines() if l.strip()]
perf = {}
for l in perf_lines:
    m = re.match(r"(\w+)\s+(\w+):\s*([\d.]+)", l)
    if m:
        perf["%s_%s" % (m.group(1), m.group(2))] = float(m.group(3))
pb = {}
if os.path.exists("pybedtools_numbers.txt"):
    for l in rd("pybedtools_numbers.txt").strip().splitlines():
        k, _, v = l.partition("=")
        pb[k] = v

results = dict(
    measured=M, checks=checks,
    all_ok=all(c["ok"] for c in checks),
    n_checks=len(checks), n_ok=sum(1 for c in checks if c["ok"]),
    err_unsorted_msg=M["err_unsorted_msg"], err_revorder_msg=M["err_revorder_msg"],
    fisher_text=fisher_text, perf=perf, pybedtools=pb,
    map_measured=map_meas,
)
with open("results.json", "w") as f:
    json.dump(results, f, indent=2)

with open("results.txt", "w") as f:
    f.write("== 050 interval-arithmetic: bedtools measured vs expected ==\n")
    for c in checks:
        f.write("%-46s measured=%-12s expected=%-12s %s\n"
                % (c["name"], c["measured"], c["expected"],
                   "OK" if c["ok"] else "MISMATCH"))
    f.write("\nchecks: %d/%d OK, all_ok=%s\n"
            % (results["n_ok"], results["n_checks"], results["all_ok"]))
    if bad:
        f.write("map mismatched genes: %s\n" % bad[:5])
    if bad_gb:
        f.write("groupby mismatched genes: %s\n" % bad_gb[:5])

print(rd("results.txt"))
