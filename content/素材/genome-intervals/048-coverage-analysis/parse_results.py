#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
048 coverage-analysis 指标解析：全部从本次 WSL 真跑产物读取，无虚构数字。
输入：design.json / cov_hist.txt / cov.bedGraph / frag.bedGraph / per_target.bed /
     per_target_mean.bed / samtools_coverage.txt / depth_all.txt.gz /
     depth_high_naive.txt / depth_high_s.txt / depth_zero.txt / _run.log(flagstat)
输出：results.json（对账表、全基因组分布、breadth、evenness、mate-overlap 对照、
     2 kb 分箱轨迹、直方图数组）；并将 cov/frag bedGraph 压缩为 .gz
"""
import gzip
import json
import math
import os

BASE = os.path.dirname(os.path.abspath(__file__))


def load_json(name):
    with open(os.path.join(BASE, name), encoding="utf-8") as f:
        return json.load(f)


def open_maybe_gz(name, mode="rt"):
    """优先读 .gz 变体（此前解析轮已压缩），否则读原文件。"""
    gz = os.path.join(BASE, name + ".gz")
    if os.path.exists(gz):
        return gzip.open(gz, mode)
    return open(os.path.join(BASE, name), encoding="utf-8")


design = load_json("design.json")
CHR = design["chrom"]
LEN = design["chrom_length_bp"]

R = {"date": "2026-09-03", "seed": design["seed"]}

# ---------- [0] 版本与比对率（来自 _run.log / _bwa.log） ----------
with open(os.path.join(BASE, "_run.log"), encoding="utf-8", errors="replace") as f:
    runlog = f.read()
ver = {}
for line in runlog.split("\n"):
    if line.startswith("bedtools v"):
        ver["bedtools"] = line.split()[0].lstrip("bedtools ").strip() if False else line.strip()
    if line.startswith("samtools "):
        ver["samtools"] = line.strip()
ver_lines = [l.strip() for l in runlog.split("\n") if "Version: 0.7" in l]
if ver_lines:
    ver["bwa"] = ver_lines[0].strip()
R["versions"] = ver

mapped = total = None
for line in runlog.split("\n"):
    if "in total" in line and line.split() and line.split()[0].isdigit():
        total = int(line.split()[0])
    if " mapped (" in line and line.split() and line.split()[0].isdigit():
        mapped = int(line.split()[0])
        mapped_pct = line.split("(")[1].split("%")[0]
if total:
    R["mapping"] = {"total_reads": total,
                    "mapped_reads": mapped,
                    "mapped_pct": float(mapped_pct)}

# ---------- [1] 逐目标对账：设计（design.json） vs 实测（per_target） ----------
mean_bed = {}
with open(os.path.join(BASE, "per_target_mean.bed"), encoding="utf-8") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        mean_bed[p[3]] = float(p[4])

targets = []
with open(os.path.join(BASE, "per_target.bed"), encoding="utf-8") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        name = p[3]
        e = design["expected_depth"][name]
        meas = mean_bed[name]
        exp_read = e["expected_mean_read_depth"]
        targets.append({
            "name": name, "start1": int(p[1]) + 1, "end1": int(p[2]),
            "length_bp": int(p[2]) - int(p[1]),
            "overlap_reads": int(p[4]), "covered_bases": int(p[5]),
            "breadth_pct": round(100.0 * float(p[7]), 2),
            "measured_mean_read_depth": round(meas, 2),
            "expected_mean_read_depth": exp_read,
            "expected_mean_frag_depth": e["expected_mean_frag_depth"],
            "measured_over_expected": round(meas / exp_read, 3) if exp_read > 0 else None,
        })
R["targets"] = targets

# ---------- [2] 全基因组直方图（cov_hist.txt，chrS 行） ----------
hist = []  # (depth, bases, fraction)
with open(os.path.join(BASE, "cov_hist.txt"), encoding="utf-8") as f:
    for line in f:
        p = line.split()
        if p[0] == CHR:
            hist.append((int(p[1]), int(p[2]), float(p[4])))
hist.sort()
tot_bases = sum(h[1] for h in hist)
gmean = sum(d * b for d, b, _ in hist) / tot_bases
cum = 0.0
median = None
d80 = None
breadth = {}
for d, b, frac in hist:
    cum += frac
    if median is None and cum >= 0.5:
        median = d
    if d80 is None and cum < 0.2:  # breadth(<d80) drops below 20% => >=d80 holds 80%+
        pass
breadth = {}
cum = 0.0
fracs = {d: f for d, b, f in hist}
# breadth(>=t) = 1 - cumfrac(<t)
cum = 0.0
b_list = []
for d, b, frac in hist:
    b_list.append((d, 1.0 - cum))
    cum += frac
def breadth_at(t):
    for d, br in b_list:
        if d >= t:
            return round(br * 100.0, 2)
    return 0.0
R["genome_dist"] = {
    "mean_depth": round(gmean, 2),
    "median_depth": median,
    "mean_over_median": round(gmean / median, 3) if median else None,
    "breadth_ge_1x_pct": breadth_at(1),
    "breadth_ge_5x_pct": breadth_at(5),
    "breadth_ge_10x_pct": breadth_at(10),
    "breadth_ge_20x_pct": breadth_at(20),
    "breadth_ge_30x_pct": breadth_at(30),
    "max_depth": hist[-1][0],
    "total_bases": tot_bases,
}
R["_breadth_curve"] = [[d, round(br, 6)] for d, br in b_list]
R["_hist"] = [[d, b] for d, b, _ in hist]

# ---------- [3] evenness：逐碱基（depth_all.txt.gz） ----------
n = 0
s1 = 0
s2 = 0
track = []
BIN = LEN // 1000
cur_sum = 0
with gzip.open(os.path.join(BASE, "depth_all.txt.gz"), "rt") as f:
    for line in f:
        p = line.split()
        d = int(p[2])
        n += 1
        s1 += d
        s2 += d * d
        cur_sum += d
        if n % BIN == 0:
            track.append(round(cur_sum / BIN, 2))
            cur_sum = 0
mean_d = s1 / n
var_d = s2 / n - mean_d * mean_d
R["evenness"] = {
    "bases": n,
    "mean_depth": round(mean_d, 3),
    "std_depth": round(math.sqrt(var_d), 3),
    "cv": round(math.sqrt(var_d) / mean_d, 3),
    "fano_factor": round(var_d / mean_d, 2),
    "poisson_ideal": 1.0,
}
R["_track_2kb"] = track

# ---------- [4] mate-overlap 对照（HIGH 区，短插入 130 bp 文库） ----------
def mean_depth_file(path, exp_lo, exp_hi):
    s = 0
    cnt = 0
    with open(os.path.join(BASE, path), encoding="utf-8") as f:
        for line in f:
            p = line.split()
            pos = int(p[1])
            if exp_lo <= pos <= exp_hi:
                s += int(p[2])
                cnt += 1
    return s / cnt

hi = design["expected_depth"]["HIGH"]
lo1, hi1 = hi["start1"], hi["end1"]
naive = mean_depth_file("depth_high_naive.txt", lo1, hi1)
dedup = mean_depth_file("depth_high_s.txt", lo1, hi1)
# 片段口径：frag.bedGraph（-pc）加权均值
s = 0.0
with open_maybe_gz("frag.bedGraph") as f:
    for line in f:
        p = line.split()
        if p[0] != CHR:
            continue
        a, b, d = int(p[1]), int(p[2]), float(p[3])
        ov = max(0, min(b, hi1) - max(a, lo1 - 1))
        if ov:
            s += d * ov
frag_mean = s / (hi1 - lo1 + 1)
R["mate_overlap_HIGH"] = {
    "insert_bp": design["insert_bp"]["high_extra"],
    "mate_overlap_bp": 2 * design["read_length_bp"] - design["insert_bp"]["high_extra"],
    "naive_read_mean_depth": round(naive, 2),
    "samtools_depth_s_mean": round(dedup, 2),
    "genomecov_pc_frag_mean": round(frag_mean, 2),
    "expected_frag_depth": hi["expected_mean_frag_depth"],
    "naive_over_frag": round(naive / frag_mean, 3),
}

# ---------- [5] 零覆盖区检出 ----------
z = design["expected_depth"]["ZERO"]
n_pos = 0
n_lines = 0
with open(os.path.join(BASE, "depth_zero.txt"), encoding="utf-8") as f:
    for line in f:
        n_lines += 1
        if int(line.split()[2]) > 0:
            n_pos += 1
# genomecov -bga：ZERO 区间内是否有非零 span
zero_spans = 0
nonzero_span_bases = 0
with open_maybe_gz("cov.bedGraph") as f:
    for line in f:
        p = line.split()
        if p[0] != CHR:
            continue
        a, b, d = int(p[1]), int(p[2]), int(p[3])
        if d > 0:
            ov = max(0, min(b, z["end1"]) - max(a, z["start1"] - 1))
            if ov:
                zero_spans += 1
                nonzero_span_bases += ov
R["zero_region"] = {
    "start1": z["start1"], "end1": z["end1"], "length_bp": z["length_bp"],
    "samtools_depth_positions": n_lines,
    "positions_with_reads": n_pos,
    "bedtools_bga_nonzero_spans": zero_spans,
    "nonzero_span_bases": nonzero_span_bases,
    "bedtools_coverage_breadth_pct": targets[3]["breadth_pct"],
}

# ---------- [6] samtools coverage 表 ----------
with open(os.path.join(BASE, "samtools_coverage.txt"), encoding="utf-8") as f:
    R["samtools_coverage"] = [l.rstrip("\n") for l in f if l.strip() and not l.startswith("#")]

with open(os.path.join(BASE, "results.json"), "w", encoding="utf-8") as f:
    json.dump(R, f, indent=2, ensure_ascii=False)

# 压缩大 bedGraph，删除逐碱基大文件（若本轮已读 .gz 则跳过）
for src in ("cov.bedGraph", "frag.bedGraph"):
        raw = os.path.join(BASE, src)
        if os.path.exists(raw):
            with open(raw, "rb") as fi, gzip.open(raw + ".gz", "wb") as fo:
                fo.write(fi.read())
            os.remove(raw)

print("== reconciliation ==")
for t in targets:
    print("%-5s read exp %.2fx / meas %.2fx (ratio %s)  frag exp %.2fx  breadth %.2f%%"
          % (t["name"], t["expected_mean_read_depth"], t["measured_mean_read_depth"],
             ("%.3f" % t["measured_over_expected"]) if t["measured_over_expected"] is not None else "n/a",
             t["expected_mean_frag_depth"], t["breadth_pct"]))
print("== genome ==")
print(json.dumps(R["genome_dist"], indent=1))
print(json.dumps(R["evenness"], indent=1))
print("== mate_overlap_HIGH ==")
print(json.dumps(R["mate_overlap_HIGH"], indent=1))
print("== zero_region ==")
print(json.dumps(R["zero_region"], indent=1))
print("== mapping ==")
print(json.dumps(R.get("mapping", {}), indent=1))
print("== versions ==")
print(json.dumps(R["versions"], indent=1))
