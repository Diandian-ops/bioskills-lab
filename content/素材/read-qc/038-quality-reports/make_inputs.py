#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_inputs.py — 038 quality-reports 造数据（可复现，seed=20260903）

生成 3 个不同质量形态的模拟 Illumina FASTQ（各 20000 条 reads × 100 bp，gzip）：
  S1_good      高质量：Q37 上下平稳，随机序列 ~25% 每碱基，无 adapter，读段几乎无重复
  S2_degraded  质量衰减：首端 ~Q38 线性衰减到 3' 端 ~Q7（Per base sequence quality 应 FAIL）
  S3_adapter   adapter 污染：~20% reads 插入片段短于读长产生 3' read-through（Illumina
               Universal Adapter），外加 80 条完全相同的 adapter-dimer 读段（0.4%，
               Overrepresented sequences 应 WARN，duplication 抬升）

用法：python3 make_inputs.py
输出：raw_fastq/S{1,2,3}_{good,degraded,adapter}.fastq.gz
"""
import gzip
import os
import random
import sys

SEED = 20260903
N_READS = 20000
READ_LEN = 100
ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"  # Illumina Universal Adapter (TruSeq Read1)
BASES = "ACGT"

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "raw_fastq")


def rand_seq(rng, n):
    return "".join(rng.choice(BASES) for _ in range(n))


def qual_str(qs):
    return "".join(chr(33 + max(2, min(41, q))) for q in qs)


def make_good(rng):
    """高质量：Q37 上下平稳（首端略高、偶发 Q43），等概率随机碱基。
    注意：必须含 ASCII>=75 (Q>=42) 的字符，FastQC 才会判定为 Phred+33，
    否则 Q35-39 的字符落在 ASCII 68-72，会被误判为 Illumina 1.5 (Phred+64)。"""
    seq = rand_seq(rng, READ_LEN)
    qs = []
    for i in range(READ_LEN):
        base = 40 if i < 3 else 37
        # 抖动到 Q30（ASCII 63 < 64）：FastQC 只有见到 ASCII<64 的字符
        # 才会把文件判为 Phred+33，否则即使含 Q43 也判 Illumina 1.5 (+64)
        qs.append(min(43, base + rng.randint(-7, 3)))
    return seq, qual_str(qs)


def make_degraded(rng):
    """质量衰减：首端 Q42 线性降到 3' 端 Q7（带抖动）。首端 Q42 保证
    FastQC 按 Phred+33 解码（ASCII>=75），否则整文件被误判为 +64。"""
    seq = rand_seq(rng, READ_LEN)
    qs = []
    for i in range(READ_LEN):
        mean_q = 42 + (7 - 42) * i / (READ_LEN - 1)   # 42 -> 7 线性
        qs.append(int(round(mean_q)) + rng.randint(-2, 2))
    return seq, qual_str(qs)


def make_adapter(rng):
    """adapter 污染：20% 短插入 read-through；0.4% 同一 adapter-dimer 读段。
    质量首两位 Q42（保证 Phred+33 判定），其余 Q35 上下。"""
    def _qs():
        # 首两位 Q42 保证 +33 判定的辅助；其余 Q29-36（含 ASCII<64 字符，强制 Phred+33）
        return [42 if i < 2 else min(36, 35 + rng.randint(-6, 1))
                for i in range(READ_LEN)]
    roll = rng.random()
    if roll < 0.004:
        # 固定 adapter-dimer：adapter 重复填满 100 bp（所有 dimer 读段完全一致）
        seq = (ADAPTER + ADAPTER + ADAPTER + ADAPTER)[:READ_LEN]
        return seq, qual_str(_qs())
    if roll < 0.204:
        # 短插入片段 -> 3' read-through adapter
        ins = rng.randint(20, 90)
        seq = rand_seq(rng, ins) + (ADAPTER * 4)[:READ_LEN - ins]
        return seq, qual_str(_qs())
    seq = rand_seq(rng, READ_LEN)
    return seq, qual_str(_qs())


def write_fastq(path, rng, maker):
    with gzip.open(path, "wt") as f:
        for i in range(N_READS):
            seq, qual = maker(rng)
            f.write("@S3R:{:06d} 1:N:0:ATCACG\n{}\n+\n{}\n".format(i, seq, qual))


def main():
    os.makedirs(OUT, exist_ok=True)
    samples = [
        ("S1_good", make_good),
        ("S2_degraded", make_degraded),
        ("S3_adapter", make_adapter),
    ]
    # 固定 per-sample 偏移，保证跨运行可复现（str hash 受 PYTHONHASHSEED 影响不可用）
    seeds = {"S1_good": 20260903, "S2_degraded": 20260904, "S3_adapter": 20260905}
    for name, maker in samples:
        rng = random.Random(seeds[name])
        p = os.path.join(OUT, name + ".fastq.gz")
        write_fastq(p, rng, maker)
        print("wrote %s (%d reads x %d bp)" % (p, N_READS, READ_LEN))
    return 0


if __name__ == "__main__":
    sys.exit(main())
