#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""050 perf fixture: two large sorted BED files on a 250 Mb synthetic chromosome.
Independent seeds (7/8); positions strictly increasing so files are born sorted."""
import random

random.seed(7)
CHROM, LEN = "chrB1", 250000000
N = 300000

with open("big_genome.txt", "w") as f:
    f.write("%s\t%d\n" % (CHROM, LEN))

def gen(path, seed):
    rng = random.Random(seed)
    pos = 1
    with open(path, "w") as f:
        for _ in range(N):
            L = rng.randint(100, 1000)
            s = pos + rng.randint(0, 5)
            e = s + L
            if e >= LEN:
                break
            f.write("%s\t%d\t%d\n" % (CHROM, s, e))
            pos = e

gen("big_a.bed", 7)
gen("big_b.bed", 8)
print("perf fixtures: big_a.bed / big_b.bed, %d intervals each, %s %d bp" % (N, CHROM, LEN))
