#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
_analyze.py - parse bedtools fisher outputs + permutation nulls, cross-validate
the analytic Fisher p-value with a pure-python exact test, write summary.tsv.

Pure-python two-sided Fisher's exact test (no scipy):
  pmf(k) = C(a+b, k) * C(c+d, d-k) / C(n, d)  over the 2x2 table (a b / c d)
  two-sided p = sum of pmf(k) <= pmf(k_obs) * (1 + 1e-9)
  computed in log space with math.lgamma; floats underflowing to 0 are reported
  alongside their log10.
"""
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))

def read_tsv(path, col_types=(str, float)):
    rows = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for ln in f:
            if not ln.strip():
                continue
            vals = ln.rstrip("\n").split("\t")
            rows.append([t(v) for t, v in zip(col_types, vals)])
    return header, rows

# ---- pure-python Fisher exact ---------------------------------------------
def logchoose(n, k):
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

def fisher_exact_2x2(a, b, c, d):
    """two-sided Fisher exact p for table [[a, b], [c, d]]; returns (p, log10p)."""
    n = a + b + c + d
    r1 = a + b
    c1 = a + c
    p_obs = logchoose(r1, a) + logchoose(n - r1, c1 - a) - logchoose(n, c1)
    total = 0.0
    for k in range(max(0, c1 - (n - r1)), min(r1, c1) + 1):
        lp = logchoose(r1, k) + logchoose(n - r1, c1 - k) - logchoose(n, c1)
        if lp <= p_obs + 1e-9:
            total += math.exp(lp - p_obs)  # scale by obs to dodge underflow
    p_rel = total  # = p / exp(p_obs)
    log10p = p_obs / math.log(10) + math.log10(total) if total > 0 else float("-inf")
    p = p_rel * math.exp(p_obs)  # underflows to 0.0 when p < ~1e-308
    return p, log10p

def parse_fisher(path):
    a = b = c = d = None
    pvals = None
    with open(path) as f:
        for ln in f:
            ln = ln.rstrip("\n")
            if ln.startswith("#     in -a"):
                cells = [x.strip() for x in ln.replace("#", "").split("|") if x.strip()]
                a, b = int(cells[1]), int(cells[2])
            elif ln.startswith("# not in -a"):
                cells = [x.strip() for x in ln.replace("#", "").split("|") if x.strip()]
                c, d = int(cells[1]), int(cells[2])
            elif ln and not ln.startswith("#") and "\t" in ln:
                parts = ln.split("\t")
                try:
                    vals = [float(x) for x in parts]
                except ValueError:
                    continue
                if len(vals) == 4:
                    pvals = vals  # left, right, two-tail, ratio
    return (a, b, c, d), pvals

def main():
    lines = []
    header = ("set", "n_int", "obs_jaccard",
              "null_mode", "null_mean", "null_median", "null_sd",
              "fold_vs_null", "zscore", "emp_p",
              "fisher_a", "fisher_b", "fisher_c", "fisher_d",
              "fisher_two_tail_p", "fisher_OR",
              "py_fisher_p", "py_fisher_log10p")
    lines.append("\t".join(header))

    fisher = {}
    frows = ["set\ta\tb\tc\td\tbedtools_two_tail_p\tbedtools_OR\tpy_fisher_p\tpy_log10p"]
    for s in ("enriched", "random"):
        table, pv = parse_fisher(os.path.join(HERE, "fisher_%s.txt" % s))
        fisher[s] = (table, pv)
        a, b, c, d = table
        p_py, log10p = fisher_exact_2x2(a, b, c, d)
        print("[%s] fisher table a=%d b=%d c=%d d=%d" % (s, a, b, c, d))
        print("[%s] bedtools: two-tail=%s ratio=%s" % (s, pv[2], pv[3]))
        print("[%s] python  : two-sided p=%.6g log10p=%.2f" % (s, p_py, log10p))
        frows.append("%s\t%d\t%d\t%d\t%d\t%.6g\t%s\t%.6g\t%.2f"
                     % (s, a, b, c, d, pv[2], pv[3], p_py, log10p))
    with open(os.path.join(HERE, "fisher_summary.tsv"), "w") as f:
        f.write("\n".join(frows) + "\n")

    _, obs_rows = read_tsv(os.path.join(HERE, "observed_jaccard.tsv"))
    obs = {name: val for name, val in obs_rows}

    for s in ("enriched", "random"):
        (a, b, c, d), pv = fisher[s]
        o = obs[s]
        for mode in ("matched", "uniform"):
            vals = []
            with open(os.path.join(HERE, "nulls_%s_%s.tsv" % (mode, s))) as f:
                for ln in f:
                    if ln.strip():
                        vals.append(float(ln))
            n = len(vals)
            mean = sum(vals) / n
            srt = sorted(vals)
            median = srt[n // 2]
            sd = (sum((x - mean) ** 2 for x in vals) / (n - 1)) ** 0.5
            hits = sum(1 for x in vals if x >= o)
            emp_p = (hits + 1) / (n + 1)
            fold = o / mean
            z = (o - mean) / sd
            lines.append("\t".join([
                s, "600", "%.6g" % o, mode,
                "%.6g" % mean, "%.6g" % median, "%.6g" % sd,
                "%.4g" % fold, "%.3f" % z, "%.6g" % emp_p,
                str(a), str(b), str(c), str(d),
                ("%.6g" % pv[2]) if pv else "NA",
                ("%.6g" % pv[3]) if pv and pv[3] == pv[3] and pv[3] != float("inf") else "inf",
                "", "",
            ]))
            print("[%s/%s] obs=%.6g null mean=%.6g median=%.6g sd=%.6g fold=%.3f z=%.2f emp_p=%.4f (hits=%d/N=%d)"
                  % (s, mode, o, mean, median, sd, fold, z, emp_p, hits, n))

    with open(os.path.join(HERE, "summary.tsv"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("wrote summary.tsv")

if __name__ == "__main__":
    main()
