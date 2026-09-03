#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
045 bed-file-basics: self-contained simulated dataset (fixed seed=45).
Outputs (all in this directory):
  chr1.fa            1 chromosome "chr1", 2,000,000 bp random ACGT
  genes.bed          BED6  gene bodies (chr1, start, end, name, score, strand)
  exons.bed          BED6  exons inside genes
  cpg.bed            BED6  CpG islands (12 isolated + 6 clusters x 3 overlapping)
  peaks.bed          narrowPeak-style BED6+4 (signal, p, q, peak offset; one peak=-1)
  transcripts.bed12  BED12 transcript models built from first 8 genes' exons
  variants.vcf       30 SNVs, 1-based POS (for VCF -> BED start-1 conversion)
  landmark.gff       1-bp 1-based-closed feature at chr1:1000-1000 (round-trip test)
  reads.sam          6 minimal alignments incl. one spliced (N CIGAR) for bamtobed -split
  cpg_bare.bed       copy of cpg.bed with chrom renamed chr1 -> 1 (mismatch demo)
  cpg_crlf.bed       first 5 lines of cpg.bed with CRLF endings (CRLF demo)
Run inside WSL bio env: python3 make_inputs.py
"""
import os
import random

random.seed(45)

HERE = os.path.dirname(os.path.abspath(__file__))
CHR = "chr1"
CHR_LEN = 2_000_000


def w(name, lines):
    with open(os.path.join(HERE, name), "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")


# ---------- reference FASTA ----------
seq = "".join(random.choice("ACGT") for _ in range(CHR_LEN))
with open(os.path.join(HERE, "chr1.fa"), "w", newline="\n") as f:
    f.write(">" + CHR + "\n")
    for i in range(0, CHR_LEN, 60):
        f.write(seq[i:i + 60] + "\n")

# ---------- genes (BED6), non-overlapping walk ----------
genes = []  # (start, end, name, score, strand)
pos = 30_000
gi = 0
while pos < CHR_LEN - 30_000 and gi < 40:
    L = random.randint(3000, 15000)
    end = pos + L
    if end > CHR_LEN - 1000:
        break
    gi += 1
    genes.append((pos, end, "gene%03d" % gi, random.randint(100, 1000),
                  random.choice("+-")))
    pos = end + random.randint(20000, 45000)

w("genes.bed", ["%s\t%d\t%d\t%s\t%d\t%s" % (CHR, s, e, n, sc, st)
                for s, e, n, sc, st in genes])

# ---------- exons (BED6), 2-8 per gene, ascending inside gene body ----------
exons = []  # (start, end, name, score, strand)
gene_exons = {}  # gene name -> list of (start,end)
for s, e, n, sc, st in genes:
    glen = e - s
    n_ex = random.randint(2, 8)
    slot = max(glen // n_ex, 400)
    cur = []
    for i in range(n_ex):
        off = i * slot + (random.randint(0, slot // 3) if slot > 300 else 0)
        ex_len = random.randint(50, min(300, slot - 50))
        ex_s = s + off
        ex_e = min(ex_s + ex_len, e)
        if ex_s >= ex_e or ex_e > e:
            continue
        cur.append((ex_s, ex_e))
    if len(cur) < 2:
        cur = [(s, s + min(300, glen // 2)), (e - min(300, glen // 2), e)]
    cur.sort()
    gene_exons[n] = cur
    for j, (a, b) in enumerate(cur, 1):
        exons.append((a, b, "%s.E%d" % (n, j), sc, st))

w("exons.bed", ["%s\t%d\t%d\t%s\t%d\t%s" % (CHR, a, b, nm, sc, st)
                for a, b, nm, sc, st in exons])

# ---------- CpG islands (BED6): 12 isolated + 6 clusters of 3 overlapping ----------
cpg = []
for i in range(12):
    s = random.randint(10_000, CHR_LEN - 10_000)
    L = random.randint(500, 2000)
    cpg.append((s, s + L))
for c in range(6):
    s = random.randint(10_000, CHR_LEN - 10_000)
    for k in range(3):
        off = k * random.randint(400, 900)
        L = random.randint(600, 1500)
        cpg.append((s + off, s + off + L))
cpg.sort()

w("cpg.bed", ["%s\t%d\t%d\tcpg%03d\t%d\t." % (CHR, s, e, i + 1,
                                              random.randint(100, 1000))
              for i, (s, e) in enumerate(cpg)])

# ---------- narrowPeak-style peaks (BED6+4) ----------
np_lines = []
for i in range(10):
    s = random.randint(50_000, CHR_LEN - 50_000)
    L = random.randint(200, 800)
    peak = random.randint(50, L - 50) if i < 9 else -1
    np_lines.append("\t".join(str(x) for x in [
        CHR, s, s + L, "peak%02d" % (i + 1), random.randint(100, 1000), ".",
        round(random.uniform(1.0, 30.0), 2),
        round(random.uniform(3.0, 20.0), 2),
        round(random.uniform(0.0, 15.0), 2), peak]))
w("peaks.bed", np_lines)

# ---------- BED12 transcripts from first 8 genes ----------
# chromStart/chromEnd = first/last exon bounds so blockStarts[0] == 0
b12 = []
for s, e, n, sc, st in genes[:8]:
    blocks = gene_exons[n]
    first_s = blocks[0][0]
    last_e = blocks[-1][1]
    if n in ("gene001", "gene005"):        # non-coding: thickStart = thickEnd
        ts, te = first_s, first_s
    else:                                   # CDS = inner span
        ts, te = first_s + 100, last_e - 100
        if ts >= te:
            ts, te = first_s, last_e
    sizes = [b - a for a, b in blocks]
    starts = [a - first_s for a, b in blocks]
    b12.append("\t".join(str(x) for x in [
        CHR, first_s, last_e, "mRNA_" + n.replace("gene", ""), sc, st, ts, te, 0,
        len(blocks), ",".join(map(str, sizes)) + ",",
        ",".join(map(str, starts)) + ","]))
w("transcripts.bed12", b12)

# ---------- VCF: 30 SNVs at 1-based POS (one at POS 1000 as landmark) ----------
vcf_pos = set([1000])
while len(vcf_pos) < 30:
    vcf_pos.add(random.randint(2000, CHR_LEN - 2000))
vcf_pos = sorted(vcf_pos)
comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
head = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=2000000>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
rows = []
for i, p in enumerate(sorted(vcf_pos), 1):
    ref = seq[p - 1]
    alt = comp[ref]
    rows.append("%s\t%d\tsnp%02d\t%s\t%s\t%d\tPASS\t." %
                (CHR, p, i, ref, alt, random.randint(50, 500)))
w("variants.vcf", head + rows)

# ---------- 1-bp GFF landmark (1-based closed): chr1 1000-1000 ----------
w("landmark.gff", ["%s\tsim\tlandmark\t1000\t1000\t.\t.\t.\tID=lm1" % CHR])

# ---------- minimal SAM for bamtobed (incl. one spliced read) ----------
sq = seq[50000:50100]
sam = ["@HD\tVN:1.6\tSO:unsorted", "@SQ\tSN:chr1\tLN:%d" % CHR_LEN]
alns = [
    ("read1", 50001, "100M", sq),
    ("read2_spliced", 150001, "50M100N50M", seq[150000:150050] + seq[150150:150200]),
    ("read3", 300001, "100M", seq[300000:300100]),
    ("read4", 999951, "100M", seq[999950:1000050]),
    ("read5", 1234567, "100M", seq[1234566:1234666]),
    ("read6", 1999901, "100M", seq[1999900:2000000]),
]
for nm, p, cig, s in alns:
    sam.append("\t".join([nm, "0", CHR, str(p), "60", cig, "*", "0", "0",
                          s, "I" * len(s)]))
w("reads.sam", sam)

# ---------- failure-demo inputs ----------
w("cpg_bare.bed", [l.replace(CHR + "\t", "1\t", 1)
                   for l in ["%s\t%d\t%d\tcpg%03d\t%d\t." % (CHR, s, e, i + 1, 500)
                             for i, (s, e) in enumerate(cpg)]])

with open(os.path.join(HERE, "cpg_crlf.bed"), "w", newline="") as f:
    for i, (s, e) in enumerate(cpg[:5]):
        f.write("%s\t%d\t%d\tcpg%03d\t500\t.\r\n" % (CHR, s, e, i + 1))

# ---------- summary ----------
print("chr1_len=%d" % CHR_LEN)
print("genes=%d" % len(genes))
print("exons=%d" % len(exons))
print("cpg=%d" % len(cpg))
print("peaks=10")
print("transcripts_bed12=%d" % len(b12))
print("variants=%d" % len(rows))
print("sam_reads=%d" % len(alns))
