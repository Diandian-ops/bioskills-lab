#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
041 fastp-workflow 造数据：模拟 Illumina PE FASTQ（TruSeq 接头 + 质量缺陷 + N 污染）。
设计目标：让 fastp 的每道工序都有活干
  - 接头去除（PE overlap 分析）：20% 读对 insert <= read_len，产生 read-through 接头
  - 碱基校正 --correction：25% 读对部分重叠（110-199 bp），重叠区 1% 非对称质量错误
  - 低质量过滤（-q 20 -u 40）：12% 读对整体质量沿 3' 衰减至 Q20 以下占比 >40%
  - N 过滤（-n 5）：3% 读对含 6-12 连续 N
  - 长度过滤（-l 36）：5% 读对 insert < 36，裁剪后必然过短
固定种子 20260412，可复现。
"""
import gzip
import os
import random

random.seed(20260412)

BASE = os.path.dirname(os.path.abspath(__file__))
N_PAIRS = 30000
READ_LEN = 100
ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
BASES = "ACGT"
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def revcomp(s):
    return "".join(COMPLEMENT[c] for c in reversed(s))


def rand_seq(n):
    return "".join(random.choice(BASES) for _ in range(n))


def mut_seq(seq, err_rate=0.01):
    """1% 随机测序错误（替换为其他碱基），保持可重叠不一致 -> 供 --correction。"""
    out = []
    for b in seq:
        if b != "N" and random.random() < err_rate:
            out.append(random.choice([c for c in BASES if c != b]))
        else:
            out.append(b)
    return "".join(out)


def qual_profile(n, low=False):
    """质量沿 3' 衰减；low=True 时 Q20 以下占比超 40%（触发 -q 20 -u 40 过滤）。"""
    qs = []
    for i in range(n):
        if low:
            q = 32 - 0.35 * i + random.gauss(0, 3)
        else:
            q = 36 - 0.10 * i + random.gauss(0, 3)
        q = int(max(5, min(40, q)))
        qs.append(q)
    return "".join(chr(33 + q) for q in qs)


def draw_insert():
    r = random.random()
    if r < 0.50:
        return random.randint(180, 320)   # 无重叠、无接头
    if r < 0.75:
        return random.randint(110, 199)   # 部分重叠（供 --correction）
    if r < 0.95:
        return random.randint(45, 100)    # read-through 接头
    return random.randint(20, 44)         # 极短（裁剪后触发 -l 36）


def make_pair(idx):
    """返回 (name1, seq1, qual1, name2, seq2, qual2)。"""
    insert = draw_insert()
    frag = mut_seq(rand_seq(insert))
    is_low = random.random() < 0.12          # 低质量读对
    has_ns = random.random() < 0.03          # N 污染读对

    # R1 正向读 insert 前段，不足 100 bp 补 R1 接头（read-through，接头循环填充保 100bp）
    if insert >= READ_LEN:
        s1 = frag[:READ_LEN]
    else:
        s1 = frag + (ADAPTER_R1 * 4)[:READ_LEN - insert]
    # R2 读 insert 后段的反向互补，不足补 R2 接头
    if insert >= READ_LEN:
        s2 = revcomp(frag[insert - READ_LEN:])
    else:
        s2 = revcomp(frag) + (ADAPTER_R2 * 4)[:READ_LEN - insert]

    q1 = qual_profile(len(s1), low=is_low)
    q2 = qual_profile(len(s2), low=is_low)

    # 纠错素材：重叠区（insert <= 199 时存在）注入错配——R1 高质量、R2 低质量，
    # 供 fastp --correction 用高质量碱基覆盖低质量碱基。错配位置间隔 >= 8，
    # 避免 --cut_right 的滑动窗口把低质量点切掉。
    n_fodder = 0
    if insert <= 199 and not is_low and random.random() < 0.6:
        if insert <= 100:
            f_lo, f_hi = 10, insert - 10        # 片段坐标下的重叠区
        else:
            f_lo, f_hi = insert - 100 + 10, 90
        if f_hi > f_lo:
            spots = sorted(random.sample(range(f_lo, f_hi),
                                         min(2, f_hi - f_lo)))
            last = -10
            for f in spots:
                if f - last < 8:
                    continue
                last = f
                i2 = insert - 1 - f             # R2 上对应重叠位置
                cur = s2[i2]
                s2 = s2[:i2] + random.choice([c for c in BASES if c != cur]) + s2[i2 + 1:]
                q2 = q2[:i2] + chr(33 + 7) + q2[i2 + 1:]
                if ord(q1[f]) - 33 < 35:
                    q1 = q1[:f] + chr(33 + 36) + q1[f + 1:]
                n_fodder += 1

    if has_ns:
        run = random.randint(6, 12)
        pos = random.randint(5, READ_LEN - run - 5)
        lst1, lst2 = list(s1), list(s2)
        for j in range(pos, min(pos + run, len(lst1))):
            lst1[j] = "N"
        for j in range(pos, min(pos + run, len(lst2))):
            lst2[j] = "N"
        s1, s2 = "".join(lst1), "".join(lst2)
        q1 = q1[:pos] + "2" * min(run, len(q1) - pos) + q1[pos + run:]
        q2 = q2[:pos] + "2" * min(run, len(q2) - pos) + q2[pos + run:]

    return ("@sim_pair:%06d 1:N:0:SIM" % idx, s1, q1,
            "@sim_pair:%06d 2:N:0:SIM" % idx, s2, q2)


def write_fastq(path, records):
    with gzip.open(path, "wt") as f:
        for name, seq, qual in records:
            f.write("%s\n%s\n+\n%s\n" % (name, seq, qual))


def main():
    r1_records, r2_records = [], []
    for i in range(N_PAIRS):
        n1, s1, q1, n2, s2, q2 = make_pair(i)
        r1_records.append((n1, s1, q1))
        r2_records.append((n2, s2, q2))

    write_fastq(os.path.join(BASE, "raw_R1.fq.gz"), r1_records)
    write_fastq(os.path.join(BASE, "raw_R2.fq.gz"), r2_records)

    n_adapter = sum(1 for _, s, _ in r1_records if "AGATCG" in s)
    n_ns = sum(1 for _, s, _ in r1_records if "NNNNNN" in s)
    n_badlen = sum(1 for _, s, _ in r1_records if len(s) != READ_LEN)
    print("pairs written: %d (read length %d bp, off-length reads: %d)"
          % (N_PAIRS, READ_LEN, n_badlen))
    print("R1 with adapter read-through (AGATCG...): %d (%.1f%%)"
          % (n_adapter, 100.0 * n_adapter / N_PAIRS))
    print("R1 with N-run (>=6 Ns): %d (%.1f%%)" % (n_ns, 100.0 * n_ns / N_PAIRS))
    print("low-quality pairs (design rate 12 pct): about %d" % int(N_PAIRS * 0.12))


if __name__ == "__main__":
    main()
