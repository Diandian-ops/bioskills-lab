#!/usr/bin/env python3
"""046-bedgraph-handling: self-contained synthetic dataset.

One simulated chromosome (chr1, 2,000,000 bp) with 4 embedded enrichment
peaks. Sample A is generated from the per-base depth profile; sample B is a
random 25% subsample of the SAME read set -> identical biology, ~4x shallower
sequencing depth. All random draws use a fixed seed so the dataset is
reproducible.

Outputs (in the current directory):
  ref.fa          2 Mb reference (random ACGT, seed fixed)
  sampleA.sam     coordinate-sorted SAM, single-end 100 bp reads, MAPQ 60
  sampleB.sam     25% subsample of sampleA reads
  peaks.bed       the 4 embedded peaks (chrom/start/end/name/height/sigma)
  meta.json       run parameters (seed, read counts, peak table)
"""
import json
import math
import random

SEED = 20260903
CHROM = "chr1"
CHROM_LEN = 2000000
READ_LEN = 100
MAPQ = 60
BASE_DEPTH_A = 0.5          # background mean depth of sample A (x)
FRAC_B = 0.25               # sample B keeps this fraction of A's reads
# (center, sigma, peak height in x over background, name)
PEAKS = [
    (300000, 2000, 40.0, "peak1"),
    (800000, 2500, 60.0, "peak2"),
    (1200000, 1500, 30.0, "peak3"),
    (1650000, 3000, 50.0, "peak4"),
]

random.seed(SEED)


def poisson(lam):
    """Knuth sampler, adequate for small lam."""
    if lam <= 0:
        return 0
    L = math.exp(-lam)
    k = 0
    p = 1.0
    while True:
        p *= random.random()
        if p <= L:
            return k
        k += 1


# ---- reference sequence ----------------------------------------------------
seq_parts = []
full_60 = CHROM_LEN // 60
rem = CHROM_LEN % 60
for _ in range(full_60):
    seq_parts.append("".join(random.choice("ACGT") for _ in range(60)))
if rem:
    seq_parts.append("".join(random.choice("ACGT") for _ in range(rem)))
with open("ref.fa", "w") as f:
    f.write(">" + CHROM + "\n")
    for part in seq_parts:
        f.write(part + "\n")
refseq = "".join(seq_parts)

# ---- per-base depth profile -> read start counts ---------------------------
starts = []
for x in range(0, CHROM_LEN - READ_LEN):
    d = BASE_DEPTH_A
    for c, s, h, _n in PEAKS:
        t = (x - c) / s
        if -8.0 < t < 8.0:
            d += h * math.exp(-0.5 * t * t)
    n = poisson(d / READ_LEN)
    for _ in range(n):
        starts.append(x)

n_reads_a = len(starts)
b_mask = [random.random() < FRAC_B for _ in range(n_reads_a)]
starts_b = [p for p, keep in zip(starts, b_mask) if keep]
n_reads_b = len(starts_b)

# ---- SAM writing ------------------------------------------------------------
HDR = ("@HD\tVN:1.6\tSO:coordinate\n"
       "@SQ\tSN:%s\tLN:%d\n" % (CHROM, CHROM_LEN))
QUAL = "I" * READ_LEN


def write_sam(path, positions):
    with open(path, "w") as f:
        f.write(HDR)
        for i, p in enumerate(positions):
            seq = refseq[p:p + READ_LEN]
            f.write("read%d\t0\t%s\t%d\t%d\t%dM\t*\t0\t0\t%s\t%s\n"
                    % (i, CHROM, p + 1, MAPQ, READ_LEN, seq, QUAL))


write_sam("sampleA.sam", starts)
write_sam("sampleB.sam", starts_b)

# ---- side files --------------------------------------------------------------
with open("peaks.bed", "w") as f:
    for c, s, h, name in PEAKS:
        f.write("%s\t%d\t%d\t%s\theight=%.1fx\tsigma=%d\n"
                % (CHROM, c - 3 * s, c + 3 * s, name, h, s))

meta = {
    "seed": SEED,
    "chrom": CHROM,
    "chrom_length": CHROM_LEN,
    "read_length": READ_LEN,
    "base_depth_A": BASE_DEPTH_A,
    "subsample_fraction_B": FRAC_B,
    "n_reads_A": n_reads_a,
    "n_reads_B": n_reads_b,
    "depth_ratio_A_over_B": round(n_reads_a / n_reads_b, 3),
    "peaks": [
        {"name": name, "center": c, "sigma": s, "height_x": h}
        for c, s, h, name in PEAKS
    ],
}
with open("meta.json", "w") as f:
    json.dump(meta, f, indent=1)

print("reads A = %d, reads B = %d (ratio %.2fx)"
      % (n_reads_a, n_reads_b, n_reads_a / n_reads_b))
print("done.")
