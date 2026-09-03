#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
006 alignment-trimming — 真跑驱动 (WSL bio env)。
先用 MAFFT L-INS-i 把真实序列 (6 条人类 S1 丝氨酸蛋白酶) 比对成 MSA，
再按 SKILL.md 原命运行 ClipKIT (kpic-smart-gap / gappyout)、trimAl (-automated1 / -gappyout / -strictplus)、BMGE，
解析保留/移除的列，记录真实运行时间，写出修剪产物 + repro_transcript.txt + runtimes.txt + trim_data.json。
"""
import os
import time
import subprocess
from Bio import AlignIO

BASE = os.path.dirname(os.path.abspath(__file__))
SEQ = os.path.join(BASE, "sequences.fasta")
BMGE = "/opt/miniconda3/envs/bio/share/bmge-1.12-1/BMGE.jar"
TR = open(os.path.join(BASE, "repro_transcript.txt"), "w")
RT = open(os.path.join(BASE, "runtimes.txt"), "w")

def ver(cmd):
    try:
        return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, timeout=30).stdout.strip().splitlines()[0]
    except Exception as e:
        return "ver-fail: %s" % e

def run(name, cmd, capture_stdout=True, outfile=None):
    TR.write("########## CMD: %s ##########\n" % " ".join(cmd))
    tout = open(os.path.join(BASE, outfile), "w") if (capture_stdout and outfile) else open(os.devnull, "w")
    t0 = time.time()
    p = subprocess.run(cmd, stdout=tout, stderr=subprocess.PIPE, text=True)
    dt = time.time() - t0
    if capture_stdout and outfile:
        tout.close()
    TR.write("exit=%d  runtime=%.3fs\n" % (p.returncode, dt))
    TR.write("STDERR:\n%s\n" % (p.stderr.strip() or "(none)"))
    TR.write("RUNTIME %s %.3f\n\n" % (name, dt))
    RT.write("RUNTIME %s %.3f\n" % (name, dt))
    return p

# ---------- 0. build input MSA with MAFFT L-INS-i ----------
TR.write("########## ENVIRONMENT ##########\n")
TR.write("mafft:    %s\n" % ver(["mafft", "--version"]))
TR.write("clipkit:  %s\n" % ver(["clipkit", "--version"]))
TR.write("trimal:   %s\n" % ver(["trimal", "--version"]))
TR.write("bmge:     %s\n" % BMGE)
TR.write("biopython:%s\n" % ver(["python", "-c", "import Bio; print(Bio.__version__)"]))
TR.write("input sequences: %d\n\n" % sum(1 for l in open(SEQ) if l.startswith(">")))

run("input_msa", ["mafft", "--localpair", "--maxiterate", "1000", SEQ],
    capture_stdout=True, outfile="input_msa.fasta")
aln_in = AlignIO.read(os.path.join(BASE, "input_msa.fasta"), "fasta")
n = len(aln_in)
L_in = aln_in.get_alignment_length()
gap_per_col = [sum(1 for s in range(n) if str(aln_in[s].seq)[c] == "-") / n for c in range(L_in)]

def cols_of(fname):
    from io import StringIO
    raw = open(os.path.join(BASE, fname)).read().splitlines()
    kept = [ln for ln in raw if ln.startswith(">") or (ln and ln[0] not in ";#")]
    return AlignIO.read(StringIO("\n".join(kept)), "fasta").get_alignment_length()

results = []
# ---------- ClipKIT kpic-smart-gap (+ --log) ----------
run("clipkit_kpic", ["clipkit", "input_msa.fasta", "-m", "kpic-smart-gap",
                      "-o", "clipkit_kpic.fasta", "--log"], capture_stdout=False)
kept_kpic = []
logp = os.path.join(BASE, "clipkit_kpic.fasta.log")
if os.path.exists(logp):
    for line in open(logp):
        parts = line.rstrip("\n").split()
        if len(parts) >= 2 and parts[1].strip() == "keep":
            kept_kpic.append(int(parts[0]) - 1)  # to 0-based
results.append(("clipkit_kpic", "ClipKIT kpic-smart-gap", cols_of("clipkit_kpic.fasta"), kept_kpic))
# ---------- ClipKIT gappyout ----------
run("clipkit_gappy", ["clipkit", "input_msa.fasta", "-m", "gappyout", "-g", "0.9",
                       "-o", "clipkit_gappy.fasta"], capture_stdout=False)
results.append(("clipkit_gappy", "ClipKIT gappyout (g=0.9)", cols_of("clipkit_gappy.fasta"), None))
# ---------- trimAl -automated1 (+ -colnumbering) ----------
# -colnumbering prints the kept-column map to stdout; the trimmed alignment goes to -out.
# Capture stdout (map) to a separate file so it does not collide with the alignment file.
p = run("trimal_auto", ["trimal", "-in", "input_msa.fasta", "-out", "trimal_auto.fasta",
                         "-automated1", "-colnumbering"], capture_stdout=True, outfile="trimal_cols.txt")
kept_tr_auto = []
colmap = os.path.join(BASE, "trimal_cols.txt")
if os.path.exists(colmap):
    txt = open(colmap).read().replace("#ColumnsMap", "").replace("\n", "")
    kept_tr_auto = [int(x) for x in txt.replace(" ", "").split(",") if x.strip() != ""]
results.append(("trimal_auto", "trimAl -automated1", cols_of("trimal_auto.fasta"), kept_tr_auto))
# ---------- trimAl -gappyout ----------
run("trimal_gappy", ["trimal", "-in", "input_msa.fasta", "-out", "trimal_gappy.fasta",
                     "-gappyout"], capture_stdout=False)
results.append(("trimal_gappy", "trimAl -gappyout", cols_of("trimal_gappy.fasta"), None))
# ---------- trimAl -strictplus ----------
run("trimal_strictplus", ["trimal", "-in", "input_msa.fasta", "-out", "trimal_strictplus.fasta",
                           "-strictplus"], capture_stdout=False)
results.append(("trimal_strictplus", "trimAl -strictplus", cols_of("trimal_strictplus.fasta"), None))
# ---------- BMGE ----------
run("bmge", ["java", "-jar", BMGE, "-i", "input_msa.fasta", "-t", "AA", "-of", "bmge.fasta",
            "-h", "0.5", "-g", "0.2"], capture_stdout=False)
results.append(("bmge", "BMGE -h 0.5 -g 0.2", cols_of("bmge.fasta"), None))

# ---------- assemble stats ----------
tools = []
for name, nice, L_out, kept in results:
    removed = L_in - L_out
    frac = L_out / L_in
    rec = {"name": name, "nice": nice, "in_cols": L_in, "out_cols": L_out,
           "removed": removed, "retained_frac": round(frac, 4)}
    if kept is not None:
        rec["kept_indices"] = kept
        rec["n_kept"] = len(kept)
    tools.append(rec)

# gap profile of input + removed-column indices (clipkit/trimal_auto have maps)
removed_by_tool = {}
for name, nice, L_out, kept in results:
    if kept is not None:
        removed_by_tool[name] = sorted(set(range(L_in)) - set(kept))

data = {
    "input": {"n_seq": n, "in_cols": L_in,
              "gap_cols": sum(1 for g in gap_per_col if g > 0),
              "gap_per_col": [round(g, 3) for g in gap_per_col]},
    "tools": tools,
    "removed_by_tool": {k: v for k, v in removed_by_tool.items()},
}
import json
json.dump(data, open(os.path.join(BASE, "trim_data.json"), "w"), indent=2)

# ---------- debug summary to stdout (verbatim transcript tail) ----------
TR.write("########## TRIM SUMMARY ##########\n")
TR.write("input MSA: %d seqs, %d cols, %d gappy cols\n" % (n, L_in, data["input"]["gap_cols"]))
for t in tools:
    TR.write("%-16s %-28s in=%d out=%d removed=%d (%.1f%%)\n" % (
        t["name"], t["nice"], t["in_cols"], t["out_cols"], t["removed"],
        100 * (1 - t["retained_frac"])))
TR.close(); RT.close()
print("input MSA cols=%d  gappy_cols=%d" % (L_in, data["input"]["gap_cols"]))
for t in tools:
    print("%-16s out=%d removed=%d (%.1f%%) kept=%s" % (
        t["name"], t["out_cols"], t["removed"], 100 * (1 - t["retained_frac"]),
        t.get("n_kept", "-")))
print("repro_transcript.txt + runtimes.txt + trim_data.json written")
