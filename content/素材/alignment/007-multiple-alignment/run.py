#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
007 multiple-alignment — 真跑产出统计。
读 5 个真实比对产物 (mafft L-INS-i / --auto / FFT-NS-2, muscle -align, clustalo)，
用 BioPython AlignIO 算：比对长度、整体 gap 比例、无 gap 列数、平均成对一致性(PID2)，
并基于 mafft_linsi 产出成对一致性矩阵 + 逐列 gap 剖面。输出 msa_data.json 供 make_figs.py 出图。
"""
import os
import json
from Bio import AlignIO

BASE = os.path.dirname(os.path.abspath(__file__))

TOOLS = [
    ("mafft_linsi",  "MAFFT L-INS-i (--localpair --maxiterate 1000)"),
    ("mafft_auto",   "MAFFT --auto"),
    ("mafft_fftns2", "MAFFT FFT-NS-2 (--retree 2)"),
    ("muscle_align", "MUSCLE5 -align"),
    ("clustalo",     "ClustalOmega --force"),
]

def pid2_pair(aln, i, j):
    a = str(aln[i].seq); b = str(aln[j].seq)
    m = 0; tot = 0
    for x, y in zip(a, b):
        if x == '-' or y == '-':
            continue
        tot += 1
        if x == y:
            m += 1
    return (m / tot * 100.0) if tot else 0.0

def analyze(name, nice):
    path = os.path.join(BASE, name + ".fasta")
    aln = AlignIO.read(path, "fasta")
    n = len(aln)
    L = aln.get_alignment_length()
    chars = sum(str(s.seq).count('-') for s in aln)
    gap_frac = chars / (n * L)
    gapfree = sum(1 for c in range(L)
                  if all(str(aln[s].seq)[c] != '-' for s in range(n)))
    # avg pairwise identity over all pairs
    ids = [pid2_pair(aln, i, j) for i in range(n) for j in range(i + 1, n)]
    pid_avg = sum(ids) / len(ids) if ids else 0.0
    return {"name": name, "nice": nice, "len": L, "gap_frac": round(gap_frac, 4),
            "gapfree_cols": gapfree, "pid_avg": round(pid_avg, 2), "n_seq": n}

tools = [analyze(n, nice) for n, nice in TOOLS]

# runtime (from runtimes.txt written by the repro block)
rt = {}
rtfile = os.path.join(BASE, "runtimes.txt")
if os.path.exists(rtfile):
    for line in open(rtfile):
        if line.startswith("RUNTIME"):
            parts = line.split()
            if len(parts) == 3:
                try:
                    rt[parts[1]] = float(parts[2])
                except ValueError:
                    pass
for t in tools:
    t["runtime"] = rt.get(t["name"], None)

# identity matrix from mafft_linsi (most accurate)
aln = AlignIO.read(os.path.join(BASE, "mafft_linsi.fasta"), "fasta")
n = len(aln)
labels = [s.id for s in aln]
matrix = [[round(pid2_pair(aln, i, j), 1) for j in range(n)] for i in range(n)]

# gap profile (per column gap fraction) for mafft_linsi
L = aln.get_alignment_length()
gap_profile = [round(sum(1 for s in range(n) if str(aln[s].seq)[c] == '-') / n, 3)
               for c in range(L)]

# input stats
seqs = list(AlignIO.read(os.path.join(BASE, "sequences.fasta"), "fasta")) if False else None
with open(os.path.join(BASE, "sequences.fasta")) as f:
    raw = f.read().split(">")[1:]
raw_lens = []
for r in raw:
    lines = r.split("\n")
    raw_lens.append(sum(len(x.strip()) for x in lines[1:]))

data = {
    "input": {"n_seq": len(raw), "avg_len": round(sum(raw_lens) / len(raw))},
    "tools": tools,
    "identity_matrix": {"labels": labels, "matrix": matrix},
    "gap_profile": {"mafft_linsi": gap_profile},
}
with open(os.path.join(BASE, "msa_data.json"), "w") as f:
    json.dump(data, f, indent=2)

# ---- human-readable summary to stdout (goes into repro_transcript) ----
print("INPUT: %d sequences, avg length %d aa" % (len(raw), data["input"]["avg_len"]))
print("%-14s %-6s %-8s %-12s %-8s %s" % ("tool", "len", "gap%", "gapfree", "PID", "runtime_s"))
for t in tools:
    print("%-14s %-6d %-8.2f %-12d %-8.2f %s" % (
        t["name"], t["len"], t["gap_frac"] * 100, t["gapfree_cols"],
        t["pid_avg"], t["runtime"]))
print("pairwise identity matrix (mafft_linsi, PID2 %%) labels=%s" % labels)
for i in range(n):
    print("  ", labels[i], matrix[i])
print("msa_data.json written")
