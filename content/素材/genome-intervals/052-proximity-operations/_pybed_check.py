#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""052 pybedtools cross-check: mirror the CLI calls via the pybedtools API
(the exact python signatures quoted in SKILL.md) and diff against CLI output."""
import os
import pybedtools

BASE = os.path.dirname(os.path.abspath(__file__))
os.chdir(BASE)

peaks = pybedtools.BedTool("peaks.sorted.bed")
genes = pybedtools.BedTool("genes.sorted.bed")

# closest: python signature from SKILL.md
near = peaks.closest(genes, D="b", io=True, t="first")
cli = open("nearest_db_io_first.bed").read().strip().splitlines()
py_lines = [str(f).rstrip("\n") for f in near]
same_closest = py_lines == cli
print("closest pybedtools == CLI: %s (%d rows)" % (same_closest, len(py_lines)))

# window
win = peaks.window(genes, w=50000, c=True)
cli_w = open("window_counts.bed").read().strip().splitlines()
same_win = [str(f).rstrip("\n") for f in win] == cli_w
print("window   pybedtools == CLI: %s (%d rows)" % (same_win, len(cli_w)))

# slop promoters
tss = pybedtools.BedTool("tss.bed")
prom = tss.slop(g="genome.txt", s=True, l=2000, r=200)
cli_p = open("promoters.bed").read().strip().splitlines()
same_slop = [str(f).rstrip("\n") for f in prom] == cli_p
print("slop     pybedtools == CLI: %s (%d rows)" % (same_slop, len(cli_p)))
