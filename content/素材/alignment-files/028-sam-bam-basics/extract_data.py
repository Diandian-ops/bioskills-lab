#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
028 sam-bam-basics 数据抽取（真实跑）。
环境：WSL2 Ubuntu + samtools 1.22.1（apt 安装）。
输入：content/素材/read-alignment/011-bowtie2-alignment/aligned_e2e.bam
      （bowtie2 end-to-end 真跑产物，6000 reads，单端，坐标排序）
输出：
  - 本目录自带输入：aligned_e2e.bam / .bai / reference.fa / reference.fa.fai
  - sam_basics_data.json  -> make_figs.py 用（图数据）
  - raw_for_transcript.txt -> 深拆版 第4节「严格复现」原文照搬
该脚本只在 WSL 里跑一次，产物随仓库入库；make_figs.py 之后靠 JSON 即可独立出图。
"""
import subprocess, json, re, os, shutil, collections

SRC_BAM = "/mnt/d/1.WorkDir/RedBook/content/素材/read-alignment/011-bowtie2-alignment/aligned_e2e.bam"
SRC_REF = "/mnt/d/1.WorkDir/RedBook/content/素材/read-alignment/011-bowtie2-alignment/reference.fa"
OUTDIR = "/mnt/d/1.WorkDir/RedBook/content/素材/alignment-files/028-sam-bam-basics"


def run(*args):
    return subprocess.run(list(args), capture_output=True, text=True, errors="replace")


os.makedirs(OUTDIR, exist_ok=True)

# ---- 拷贝自带输入，保证素材目录独立可复现 ----
for f in [SRC_BAM, SRC_BAM + ".bai", SRC_REF]:
    shutil.copy(f, os.path.join(OUTDIR, os.path.basename(f)))
run("samtools", "faidx", os.path.join(OUTDIR, "reference.fa"))

BAM = os.path.join(OUTDIR, "aligned_e2e.bam")
REF = os.path.join(OUTDIR, "reference.fa")

FLAG_NAMES = [(1, "PAIRED"), (2, "PROPER_PAIR"), (4, "UNMAPPED"), (8, "MATE_UNMAPPED"),
              (16, "REVERSE"), (32, "MATE_REVERSE"), (64, "READ1"), (128, "READ2"),
              (256, "SECONDARY"), (512, "QC_FAIL"), (1024, "DUP"), (2048, "SUPPLEMENTARY")]


def decode_flags(flag):
    return [n for b, n in FLAG_NAMES if flag & b]


raw = []  # 给 transcript 用的原始命令+输出


def cap(cmd_list, label):
    r = run(*cmd_list)
    raw.append("## %s\n$ %s" % (label, " ".join(cmd_list)))
    raw.append(r.stdout.rstrip("\n") if r.stdout else "(no stdout)")
    if r.stderr.strip():
        raw.append("[stderr] " + r.stderr.strip())
    raw.append("")
    return r


hdr = cap(["samtools", "view", "-H", BAM], "view -H (header)")
view = cap(["samtools", "view", BAM], "view (alignments, 前2行预览见 transcript)")
lines = view.stdout.splitlines()
total = len(lines)

flag_counts = collections.Counter()
mapq_dist = collections.Counter()
cigar_ops = collections.Counter()
readlen_dist = collections.Counter()
refs = collections.Counter()
sample_flags = []
for ln in lines:
    cols = ln.split("\t")
    if len(cols) < 11:
        continue
    flag = int(cols[1])
    mapq = int(cols[4])
    cigar = cols[5]
    flag_counts[flag] += 1
    mapq_dist[mapq] += 1
    readlen_dist[len(cols[9])] += 1
    refs[cols[2]] += 1
    for num, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        cigar_ops[op] += int(num)
    if len(sample_flags) < 10:
        sample_flags.append({
            "qname": cols[0], "flag": flag, "names": decode_flags(flag),
            "mapq": mapq, "cigar": cigar, "rname": cols[2], "pos": int(cols[3]),
        })

cigar_feat = {"contains_insertion": 0, "contains_deletion": 0,
               "contains_clip": 0, "contains_skip": 0, "clean_100M": 0}
for ln in lines:
    cols = ln.split("\t")
    if len(cols) < 11:
        continue
    cigar = cols[5]
    has_I = "I" in cigar
    has_D = "D" in cigar
    has_clip = ("S" in cigar) or ("H" in cigar)
    has_N = "N" in cigar
    if has_I:
        cigar_feat["contains_insertion"] += 1
    if has_D:
        cigar_feat["contains_deletion"] += 1
    if has_clip:
        cigar_feat["contains_clip"] += 1
    if has_N:
        cigar_feat["contains_skip"] += 1
    if cigar == "100M":
        cigar_feat["clean_100M"] += 1

cat = {
    "total": total,
    "mapped": sum(1 for ln in lines if not (int(ln.split("\t")[1]) & 4)),
    "unmapped": sum(1 for ln in lines if int(ln.split("\t")[1]) & 4),
    "secondary": sum(1 for ln in lines if int(ln.split("\t")[1]) & 256),
    "supplementary": sum(1 for ln in lines if int(ln.split("\t")[1]) & 2048),
    "duplicates": sum(1 for ln in lines if int(ln.split("\t")[1]) & 1024),
    "qcfail": sum(1 for ln in lines if int(ln.split("\t")[1]) & 512),
    "paired": sum(1 for ln in lines if int(ln.split("\t")[1]) & 1),
    "proper_pair": sum(1 for ln in lines if int(ln.split("\t")[1]) & 2),
    "reverse": sum(1 for ln in lines if int(ln.split("\t")[1]) & 16),
    "read1": sum(1 for ln in lines if int(ln.split("\t")[1]) & 64),
    "read2": sum(1 for ln in lines if int(ln.split("\t")[1]) & 128),
}

cap(["samtools", "view", "-c", BAM], "view -c (count)")
cap(["samtools", "flagstat", BAM], "flagstat")
cap(["samtools", "view", "-H", BAM, "|", "grep", "^@PG"], "view -H | grep @PG")  # placeholder, redone below
# grep 通过管道在 subprocess 里不易做，单独跑
pg = run("bash", "-c", "samtools view -H %s | grep '^@PG'" % BAM)
flagstat_raw = run("samtools", "flagstat", BAM).stdout

# flags decode 示例
for f in [0, 4, 16, 99, 147, 256, 2048]:
    r = run("samtools", "flags", str(f))
    raw.append("## samtools flags %d" % f)
    raw.append("$ samtools flags %d" % f)
    raw.append(r.stdout.strip())
    raw.append("")

# CRAM 往返演示（需参考基因组）
cram = os.path.join(OUTDIR, "demo.cram")
bam_back = os.path.join(OUTDIR, "demo_back.bam")
r1 = run("samtools", "view", "-C", "-T", REF, "-o", cram, BAM)
r2 = run("samtools", "view", "-b", "-T", REF, "-o", bam_back, cram)
cram_count = run("samtools", "view", "-c", cram).stdout.strip()
bam_back_count = run("samtools", "view", "-c", bam_back).stdout.strip()
cram_demo = {"cram_count": int(cram_count), "bam_back_count": int(bam_back_count),
             "match": cram_count == bam_back_count,
             "exit_bam2cram": r1.returncode, "exit_cram2bam": r2.returncode}
# 清理演示产物，保持目录只留自带输入 + 图数据 + 脚本
os.remove(cram)
os.remove(bam_back)

hd = [l for l in hdr.stdout.splitlines() if l.startswith("@HD")]
sq = [l for l in hdr.stdout.splitlines() if l.startswith("@SQ")]
pg_lines = pg.stdout.splitlines() if pg.stdout else []

data = {
    "bam": "aligned_e2e.bam",
    "total_reads": total,
    "header": {"HD": hd, "SQ": sq, "PG": pg_lines},
    "categories": cat,
    "mapq_dist": {str(k): v for k, v in sorted(mapq_dist.items())},
    "cigar_ops": dict(cigar_ops),
    "cigar_feat": cigar_feat,
    "readlen_dist": {str(k): v for k, v in sorted(readlen_dist.items())},
    "refs": dict(refs),
    "sample_flags": sample_flags,
    "flagstat_raw": flagstat_raw,
    "cram_demo": cram_demo,
    "samtools_version": run("samtools", "--version").stdout.splitlines()[0],
    "reference": "reference.fa",
}

with open(os.path.join(OUTDIR, "sam_basics_data.json"), "w") as f:
    json.dump(data, f, indent=2, ensure_ascii=False)

with open(os.path.join(OUTDIR, "raw_for_transcript.txt"), "w") as f:
    f.write("\n".join(raw) + "\n")

print("WROTE sam_basics_data.json  total_reads=%d" % total)
print(json.dumps(cat, indent=0))
print("mapq_dist", dict(sorted(mapq_dist.items())))
print("cigar_ops", dict(cigar_ops))
print("cram_demo", cram_demo)
print("refs", dict(refs))
