#!/usr/bin/env python3
"""Reconcile bigWigAverageOverBed output against the design truth.

Reads: regions_avg.tab (name size covered sum mean0 mean), truth.json
Writes: parsed_results.tsv
Prints PASS/FAIL per region and an overall verdict.
Columns: name size covered sum mean0 mean + truth covered/sum/mean/mean0 + max rel err
"""
import json

TOL = 1e-3  # relative tolerance (bigWig stores float32)


def rel_diff(a, b):
    if b == 0:
        return 0.0 if a == 0 else float("inf")
    return abs(a - b) / abs(b)


def main():
    truth = json.load(open("truth.json"))
    rows = []
    all_ok = True
    worst = 0.0
    with open("regions_avg.tab") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            cols = line.split("\t")
            rows.append(cols)

    # bigWigAverageOverBed writes no header; columns: name size covered sum mean0 mean

    out = open("parsed_results.tsv", "w")
    out.write("name\tsize\tcovered\tcovered_truth\tsum\tsum_truth\t"
              "mean0\tmean0_truth\tmean\tmean_truth\t"
              "rel_err_mean0\trel_err_mean\tverdict\n")
    for cols in rows:
        name = cols[0]
        size, covered = int(cols[1]), int(cols[2])
        tool_sum = float(cols[3])
        tool_mean0 = float(cols[4])
        tool_mean = float(cols[5])
        t = truth[name]
        errs = [
            rel_diff(covered, t["covered"]),
            rel_diff(tool_sum, t["sum"]),
            rel_diff(tool_mean0, t["mean0"]),
        ]
        if t["mean"] is None:
            err_mean = 0.0  # no covered base in region: tool writes 0
            mean_ok = covered == 0
        else:
            err_mean = rel_diff(tool_mean, t["mean"])
            mean_ok = err_mean <= TOL
        ok = (errs[0] == 0 and errs[1] <= TOL and errs[2] <= TOL and mean_ok)
        worst = max(worst, errs[0], errs[1], errs[2], err_mean)
        all_ok = all_ok and ok
        out.write("%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.6f\t%.6f\t%.6f\t%s\t"
                  "%.2e\t%.2e\t%s\n"
                  % (name, size, covered, t["covered"], tool_sum, t["sum"],
                     tool_mean0, t["mean0"], tool_mean,
                     "NA" if t["mean"] is None else "%.6f" % t["mean"],
                     errs[2], err_mean, "PASS" if ok else "FAIL"))
        print("%-22s covered=%d/%d sum=%.2f/%.2f mean0=%.4f/%.4f mean=%.4f/%s -> %s"
              % (name, covered, t["covered"], tool_sum, t["sum"],
                 tool_mean0, t["mean0"], tool_mean,
                 "NA" if t["mean"] is None else "%.4f" % t["mean"],
                 "PASS" if ok else "FAIL"))
    out.close()
    print("worst relative error: %.2e (tolerance %.0e)" % (worst, TOL))
    print("RECONCILIATION: %s" % ("ALL PASS" if all_ok else "FAILED"))
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
