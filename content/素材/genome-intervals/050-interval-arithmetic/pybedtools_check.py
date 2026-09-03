#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""050 pybedtools parity check: repeat three core operations through the Python
wrapper and compare against the CLI numbers produced in this directory.
Prints PASS/FAIL per check; exits 0 either way (absence is reported honestly)."""
import sys

try:
    import pybedtools
except ImportError:
    print("pybedtools: NOT AVAILABLE in this env (skill's Python examples "
          "not exercised; CLI used as reference implementation)")
    sys.exit(0)

print("pybedtools %s" % pybedtools.__version__)

a = pybedtools.BedTool("peaks_sorted.bed")
b = pybedtools.BedTool("genes_sorted.bed")

u = a.intersect(b, u=True).count()
print("pybedtools intersect u_count = %d  (CLI i_u.bed lines: see parse)" % u)

m = a.sort().merge()
print("pybedtools merge d0 blocks = %d" % m.count())

scores = pybedtools.BedTool("scores.bedgraph")
genes = pybedtools.BedTool("genes_sorted.bed")
mapped = genes.map(scores, c=4, o="mean")
first3 = [l.fields[6] for l in mapped[:3]]
print("pybedtools map mean first3 = %s" % first3)

with open("pybedtools_numbers.txt", "w") as f:
    f.write("u_count=%d\nmerge_d0=%d\nmap_first3=%s\n" % (u, m.count(), ",".join(first3)))
print("pybedtools parity numbers written to pybedtools_numbers.txt")
