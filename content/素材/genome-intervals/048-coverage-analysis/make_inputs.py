#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
048 coverage-analysis 输入数据生成（seed 固定，可复现）。
设计：1 条 2 Mb 染色体 chrS，4 种覆盖形态的目标区间（对应 coverage-analysis
skill 的核心命题：mean 会掩盖分布差异，必须看 breadth / evenness）：
  - BG    chrS:100001-200000   (100 kb)  背景 ~10x，插片段 300 bp（read 深度==片段深度）
  - HIGH  chrS:500001-550000   (50 kb)   背景 10x + 额外短插入片段文库（130 bp，
                                            双端 100 bp 有 70 bp 重叠）追加，片段深度设计 ~50x，
                                            朴素 read 深度被 mate-overlap 双计推高
  - LOW   chrS:1000001-1030000 (30 kb)   设计 ~1.5x
  - ZERO  chrS:1500001-1520000 (20 kb)   设计 0x（不采样任何片段）
背景片段在全染色体均匀采样，但丢弃与 LOW / ZERO 相交的片段。
reads：双端 2x100 bp，1% 替换错误，全碱基 Q30；比对器 bwa-mem。
输出：ref.fa / reads_1.fq.gz / reads_2.fq.gz / design.json
（design.json 记录逐目标区间的设计期望深度，供对账）
"""
import gzip
import json
import os
import random

random.seed(20260903)

CHR = "chrS"
LEN = 2_000_000
READ_LEN = 100
INS_BG = 300          # 背景/LOW 插片段长度（read 不重叠）
INS_HIGH = 130        # HIGH 额外文库插片段长度（read 重叠 70 bp）
ERR = 0.01            # 替换错误率
QCHAR = chr(30 + 33)  # Q30

# 目标区间（1-based inclusive，BED 写出转 0-based half-open）
TARGETS = {
    "BG":   (100000, 200000),
    "HIGH": (500000, 550000),
    "LOW":  (1000000, 1030000),
    "ZERO": (1500000, 1520000),
}
EXCLUDE = [TARGETS["LOW"], TARGETS["ZERO"]]  # 背景采样排除区（0-based half-open）

N_BG = 65000          # 背景片段数（约 10x 于可采样区）
N_HIGH = 20000        # HIGH 额外片段数（片段深度 20000*130/50000 = 52x）
N_LOW = 150           # LOW 片段数（片段深度 150*300/30000 = 1.5x）

BASES = "ACGT"
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def overlaps_any(s, e, ivs):
    for a, b in ivs:
        if s < b and a < e:
            return True
    return False


def mutate(seq):
    out = []
    for b in seq:
        if random.random() < ERR:
            c = random.choice(BASES)
            out.append(c if c != b else COMP[c])
        else:
            out.append(b)
    return "".join(out)


def rc(seq):
    return "".join(COMP[b] for b in reversed(seq))


def main():
    # 参考序列
    ref = "".join(random.choice(BASES) for _ in range(LEN))
    with open("ref.fa", "w") as f:
        f.write(">%s\n" % CHR)
        for i in range(0, LEN, 60):
            f.write(ref[i:i + 60] + "\n")

    frag_plans = []  # (start0, insert_len)
    # 背景：均匀采样，丢弃与 LOW/ZERO 相交的片段
    n_rej = 0
    while len(frag_plans) < N_BG:
        s = random.randint(0, LEN - INS_BG)
        if overlaps_any(s, s + INS_BG, EXCLUDE):
            n_rej += 1
            continue
        frag_plans.append((s, INS_BG))
    # HIGH：短插入片段（read 重叠）
    for _ in range(N_HIGH):
        frag_plans.append((random.randint(TARGETS["HIGH"][0],
                                          TARGETS["HIGH"][1] - INS_HIGH), INS_HIGH))
    # LOW
    for _ in range(N_LOW):
        frag_plans.append((random.randint(TARGETS["LOW"][0],
                                          TARGETS["LOW"][1] - INS_BG), INS_BG))
    random.shuffle(frag_plans)

    r1 = gzip.open("reads_1.fq.gz", "wt", compresslevel=6)
    r2 = gzip.open("reads_2.fq.gz", "wt", compresslevel=6)
    # 设计期望深度数组（read 口径：重叠处双计，即朴素口径）
    depth = [0] * LEN
    for i, (s, ins) in enumerate(frag_plans):
        seq_f = ref[s:s + READ_LEN]
        seq_r = rc(ref[s + ins - READ_LEN:s + ins])
        name = "frag%06d" % i
        r1.write("@%s 1:N:0:1\n%s\n+\n%s\n" % (name, mutate(seq_f), QCHAR * READ_LEN))
        r2.write("@%s 2:N:0:1\n%s\n+\n%s\n" % (name, mutate(seq_r), QCHAR * READ_LEN))
        for p in range(s, min(s + READ_LEN, LEN)):
            depth[p] += 1
        for p in range(s + ins - READ_LEN, min(s + ins, LEN)):
            depth[p] += 1
    r1.close()
    r2.close()

    # 设计期望：逐目标区间朴素 read 深度与片段深度
    exp = {}
    for name, (a, b) in TARGETS.items():
        seg = depth[a:b]
        exp[name] = {
            "start1": a + 1, "end1": b, "length_bp": b - a,
            "expected_mean_read_depth": round(sum(seg) / (b - a), 2),
        }
    # 片段口径深度（HIGH 用 -pc / depth -s 口径）：背景与 LOW 片段 read 不重叠，
    # 片段深度==read 深度；HIGH 额外片段按 130 bp 计一次
    fdepth = [0] * LEN
    for s, ins in frag_plans:
        for p in range(s, min(s + ins, LEN)):
            fdepth[p] += 1
    for name, (a, b) in TARGETS.items():
        exp[name]["expected_mean_frag_depth"] = round(sum(fdepth[a:b]) / (b - a), 2)

    design = {
        "date": "2026-09-03", "seed": 20260903,
        "chrom": CHR, "chrom_length_bp": LEN, "read_length_bp": READ_LEN,
        "insert_bp": {"background": INS_BG, "high_extra": INS_HIGH},
        "error_rate": ERR, "fragments": {
            "background": N_BG, "high_extra": N_HIGH, "low": N_LOW,
            "background_rejected_overlapping": n_rej, "total": len(frag_plans)},
        "total_reads": len(frag_plans) * 2,
        "expected_depth": exp,
    }
    with open("design.json", "w") as f:
        json.dump(design, f, indent=2)

    print("ref.fa: %s %d bp" % (CHR, LEN))
    print("fragments: bg=%d (rejected %d) high=%d low=%d total=%d"
          % (N_BG, n_rej, N_HIGH, N_LOW, len(frag_plans)))
    print("reads: %d pairs x %d bp, err=%.3f" % (len(frag_plans), READ_LEN, ERR))
    for name in ["BG", "HIGH", "LOW", "ZERO"]:
        e = exp[name]
        print("  %s chrS:%d-%d (%d bp): read %.2fx / frag %.2fx"
              % (name, e["start1"], e["end1"], e["length_bp"],
                 e["expected_mean_read_depth"], e["expected_mean_frag_depth"]))


if __name__ == "__main__":
    main()
