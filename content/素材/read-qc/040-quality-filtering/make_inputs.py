#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
040 quality-filtering 输入数据生成：质量梯度模拟 FASTQ（seed 固定，可复现）。
设计（对应 quality-filtering skill 的两类问题形态）：
  - 92% 「正常」reads：逐 cycle 质量从 5' 端 Q≈40 线性衰减到 3' 端 Q≈12（坏尾），
    带逐 read 抖动；测序错误按 P(error)=10^(-Q/10) 随碱基注入。
  - 8% 「junk」reads：全 read 低质量（Q8-14 波动）（全局坏 read 亚群，需 FILTER 而非 TRIM）。
无 adapter（adapter 移除属 adapter-trimming skill，不在本 skill 范围）。
输出：input_grad.fq.gz（20000 reads × 150 bp，phred33）
"""
import gzip
import random

random.seed(20260903)

N_READS = 20000
LEN = 150
JUNK_FRAC = 0.08
OUT = "input_grad.fq.gz"

COMPILE = "ACGT"

qchars = [chr(q + 33) for q in range(41)]


def make_quality():
    """返回一条 read 的 phred 质量列表。"""
    if random.random() < JUNK_FRAC:
        # junk read：整体低质量，Q8-14 波动
        return [max(2, min(20, int(random.gauss(11, 2)))) for _ in range(LEN)]
    # 正常 read：Q(c) = 40 - 28 * c/(LEN-1)，逐 read 起点/斜率抖动
    start = random.gauss(40, 1.5)
    end = random.gauss(12, 2.0)
    qs = []
    for c in range(LEN):
        q = start + (end - start) * c / (LEN - 1)
        q += random.gauss(0, 1.5)
        qs.append(max(2, min(40, int(round(q)))))
    return qs


def make_seq(qs):
    """按质量注入测序错误：P(error) = 10^(-Q/10)。"""
    seq = []
    for q in qs:
        if random.random() < 10 ** (-q / 10.0):
            b = random.choice(COMPILE)
        else:
            b = random.choice(COMPILE)
        seq.append(b)
    return "".join(seq)


def main():
    tot_q, tot_base = 0, 0
    first_cycle_q, last_cycle_q = [], []
    n_junk = 0
    with gzip.open(OUT, "wt", compresslevel=6) as f:
        for i in range(N_READS):
            qs = make_quality()
            if max(qs) <= 20 and sum(qs) / LEN < 14:
                n_junk += 1
            seq = make_seq(qs)
            qual = "".join(qchars[q] for q in qs)
            f.write("@read%06d\n%s\n+\n%s\n" % (i, seq, qual))
            tot_q += sum(qs)
            tot_base += LEN
            first_cycle_q.append(qs[0])
            last_cycle_q.append(qs[-1])
    print("wrote %s: %d reads x %d bp" % (OUT, N_READS, LEN))
    print("input mean Q = %.2f" % (tot_q / tot_base))
    print("cycle1 mean Q = %.2f, cycle150 mean Q = %.2f"
          % (sum(first_cycle_q) / N_READS, sum(last_cycle_q) / N_READS))
    print("junk-like reads (whole-read mean Q<14): ~%d (%.1f%%, 设计值 8%%)"
          % (n_junk, 100.0 * n_junk / N_READS))


if __name__ == "__main__":
    main()
