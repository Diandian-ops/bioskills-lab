#!/usr/bin/env python3
"""Generate a simulated bedGraph signal track with known ground truth.

Design (seed=42, deterministic):
  - one chromosome chrSim, 2,000,000 bp
  - 10 bp bins -> one bedGraph interval per bin
  - background: Gaussian noise around 1.0 (clipped at 0.05)
  - 6 peaks: Gaussian shape, sigma 40 bp, amplitudes 50/100/200/400/800/400
  - a 100 kb gap with no data at all (1,500,000-1,600,000)
  - a sparse zone (1,600,000-2,000,000): only every 5th bin kept (20% coverage)
  - chrom.sizes + regions.bed + truth.json written alongside
Outputs: sim.bedGraph, chrom.sizes, regions.bed, truth.json
"""
import json
import math
import random

SEED = 42
CHROM = "chrSim"
CHROM_LEN = 2_000_000
BIN = 10
SIGMA = 40.0
SPARSE_KEEP_EVERY = 5

PEAKS = [
    (200_000, 50.0),
    (400_000, 100.0),
    (700_000, 200.0),
    (1_000_000, 400.0),
    (1_250_000, 800.0),
    (1_700_000, 400.0),
]
GAP = (1_500_000, 1_600_000)
SPARSE = (1_600_000, 2_000_000)

REGIONS = [
    # (start, end, name)
    (199_800, 200_200, "peak1_400bp"),
    (999_800, 1_000_200, "peak4_400bp"),
    (1_249_800, 1_250_200, "peak5_400bp"),
    (10_000, 60_000, "bg_50kb"),
    (500_000, 1_500_000, "wide_1Mb_peak5"),
    (1_690_000, 1_710_000, "sparse_20kb_peak6"),
    (1_500_000, 1_550_000, "gap_50kb"),
]


def main():
    random.seed(SEED)
    bins = []  # (start, end, value), coordinate-sorted by construction
    n_bins = CHROM_LEN // BIN
    for i in range(n_bins):
        start = i * BIN
        if GAP[0] <= start < GAP[1]:
            continue
        if SPARSE[0] <= start < SPARSE[1] and (i % SPARSE_KEEP_EVERY != 0):
            continue
        center = start + BIN / 2.0
        v = random.gauss(1.0, 0.2)
        if v < 0.05:
            v = 0.05
        for c, amp in PEAKS:
            v += amp * math.exp(-((center - c) ** 2) / (2.0 * SIGMA ** 2))
        bins.append((start, start + BIN, round(v, 3)))

    with open("sim.bedGraph", "w") as f:
        for s, e, v in bins:
            f.write("%s\t%d\t%d\t%.3f\n" % (CHROM, s, e, v))

    with open("chrom.sizes", "w") as f:
        f.write("%s\t%d\n" % (CHROM, CHROM_LEN))

    with open("regions.bed", "w") as f:
        for s, e, name in REGIONS:
            f.write("%s\t%d\t%d\t%s\n" % (CHROM, s, e, name))

    truth = {}
    for s, e, name in REGIONS:
        tot = 0.0
        cov = 0
        for bs, be, v in bins:
            os_, oe = max(bs, s), min(be, e)
            if os_ < oe:
                tot += v * (oe - os_)
                cov += oe - os_
        truth[name] = {
            "chrom": CHROM,
            "start": s,
            "end": e,
            "size": e - s,
            "covered": cov,
            "sum": round(tot, 3),
            "mean": None if cov == 0 else round(tot / cov, 6),
            "mean0": round(tot / (e - s), 6),
        }

    with open("truth.json", "w") as f:
        json.dump(truth, f, indent=1)

    covered_bases = sum(e - s for s, e, _ in bins)
    print("bins written: %d" % len(bins))
    print("covered bases: %d / %d (%.1f%%)"
          % (covered_bases, CHROM_LEN, 100.0 * covered_bases / CHROM_LEN))
    print("max value: %.3f" % max(v for _, _, v in bins))


if __name__ == "__main__":
    main()
