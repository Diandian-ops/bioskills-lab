#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
007 multiple-alignment — 真跑驱动 (WSL bio env)。
逐条按 SKILL.md 原命运行 mafft / muscle / clustalo，记录真实 wall-clock 时间、
命令、stderr，写出比对文件 + repro_transcript.txt + runtimes.txt。
随后 run.py 读取比对文件产出统计与 msa_data.json。
"""
import os
import time
import subprocess

BASE = os.path.dirname(os.path.abspath(__file__))
SEQ = os.path.join(BASE, "sequences.fasta")
TR = open(os.path.join(BASE, "repro_transcript.txt"), "w")
RT = open(os.path.join(BASE, "runtimes.txt"), "w")

def ver(cmd):
    try:
        return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, timeout=30).stdout.strip().splitlines()[0]
    except Exception as e:
        return "ver-fail: %s" % e

def run(name, outfile, cmd, capture_stdout=True):
    TR.write("########## CMD: %s ##########\n" % " ".join(cmd))
    tout = open(os.path.join(BASE, outfile), "w") if capture_stdout else open(os.devnull, "w")
    t0 = time.time()
    p = subprocess.run(cmd, stdout=tout, stderr=subprocess.PIPE, text=True)
    dt = time.time() - t0
    if capture_stdout:
        tout.close()
    TR.write("exit=%d  runtime=%.3fs\n" % (p.returncode, dt))
    TR.write("STDERR:\n%s\n" % (p.stderr.strip() or "(none)"))
    TR.write("RUNTIME %s %.3f\n\n" % (name, dt))
    RT.write("RUNTIME %s %.3f\n" % (name, dt))
    print("%-14s exit=%d runtime=%.3fs -> %s" % (name, p.returncode, dt, outfile))

TR.write("########## ENVIRONMENT ##########\n")
TR.write("mafft:    %s\n" % ver(["mafft", "--version"]))
TR.write("muscle:   %s\n" % ver(["muscle", "-version"]))
TR.write("clustalo: %s\n" % ver(["clustalo", "--version"]))
TR.write("biopython:%s\n" % ver(["python", "-c", "import Bio; print(Bio.__version__)"]))
nseq = sum(1 for l in open(SEQ) if l.startswith(">"))
TR.write("input sequences: %d\n\n" % nseq)

# mafft writes alignment to stdout -> capture_stdout=True
run("mafft_linsi",  "mafft_linsi.fasta",  ["mafft", "--localpair", "--maxiterate", "1000", SEQ])
run("mafft_auto",   "mafft_auto.fasta",   ["mafft", "--auto", SEQ])
run("mafft_fftns2", "mafft_fftns2.fasta", ["mafft", "--retree", "2", SEQ])
# muscle / clustalo write to their own -output / -o flag -> capture_stdout=False
run("muscle_align", "muscle_align.fasta", ["muscle", "-align", SEQ, "-output", "muscle_align.fasta", "-threads", "4"], capture_stdout=False)
run("clustalo",     "clustalo.fasta",     ["clustalo", "-i", SEQ, "-o", "clustalo.fasta", "--force"], capture_stdout=False)

TR.write("########## alignment file sizes (bytes) ##########\n")
for f in ["mafft_linsi.fasta", "mafft_auto.fasta", "mafft_fftns2.fasta", "muscle_align.fasta", "clustalo.fasta"]:
    TR.write("%s  %d\n" % (f, os.path.getsize(os.path.join(BASE, f))))
TR.close(); RT.close()
print("repro_transcript.txt + runtimes.txt written")
