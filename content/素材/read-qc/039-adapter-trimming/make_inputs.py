#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
039 adapter-trimming 造数据：模拟 Illumina TruSeq PE read-through 接头梯度数据。

设计（固定随机种子 39，可复现）：
  - 3 个梯度档：5% / 20% / 40% 的 read 对为短插入片段（read-through，
    R1/R2 两条读段都带 3' 接头），其余为长插入片段（不带接头）。
  - 读长 100 bp；read-through 对的插入片段长度均匀分布于 [67, 95] bp
    （保证剪完接头后读长 >= 36，不与 MINLEN 混淆），
    接头残留长度 = 100 - 插入片段长，即 [5, 33] bp——不超过接头序列
    本身的长度（TruSeq RT 接头 33 bp），与真实 read-through 一致。
  - 非read-through 对的插入片段长度均匀分布于 [120, 280] bp。
  - 接头序列用 SKILL.md「Verified Adapter Sequences」表的原版序列：
      R1 3': AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
      R2 3': AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
  - 全读段 1% 替换型测序错误（碱基质量相应降到 Q20），其余 Q25-Q40。

产出（同目录）：
  grad{5p,20p,40p}_R1.fq.gz / grad{5p,20p,40p}_R2.fq.gz
  truth.tsv.gz           （逐读段真值：read_id, has_adapter, adapter_bases, read_len）
  truth_summary.json     （汇总真值）
"""
import gzip
import json
import os
import random

BASE = os.path.dirname(os.path.abspath(__file__))
SEED = 39
N_PAIRS = 20000          # 每档 read 对数
READ_LEN = 100
ERROR_RATE = 0.01

ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"   # TruSeq R1 3'（SKILL.md 原版）
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"   # TruSeq R2 3'（SKILL.md 原版）

LEVELS = [("5p", 0.05), ("20p", 0.20), ("40p", 0.40)]
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def mutate(seq, rng):
    """1% 替换错误，返回 (序列, 质量字符串)。错误位 Q20，正确位 Q25-Q40。"""
    out, quals = [], []
    for b in seq:
        if rng.random() < ERROR_RATE:
            out.append(rng.choice([c for c in "ACGT" if c != b]))
            quals.append(chr(33 + 20))
        else:
            out.append(b)
            quals.append(chr(33 + rng.randint(25, 40)))
    return "".join(out), "".join(quals)


def make_pair(idx, level_frac, rng):
    """返回 (r1_seq, r1_qual, r2_seq, r2_qual, has_adapter, adapter_bases)。"""
    if rng.random() < level_frac:
        ins = rng.randint(36, 95)                     # read-through 短插入
        ad1 = ADAPTER_R1[: READ_LEN - ins]
        ad2 = ADAPTER_R2[: READ_LEN - ins]
        s1 = rng_dna(ins, rng) + ad1
        s2 = revcomp(rng_dna_cache[idx % len(rng_dna_cache)])[:ins] + ad2
        has_ad, ad_bases = 1, READ_LEN - ins
    else:
        ins = rng.randint(120, 280)                   # 长插入，无接头
        s1 = rng_dna(ins, rng)[:READ_LEN]
        s2 = revcomp(rng_dna_cache[idx % len(rng_dna_cache)])[:READ_LEN]
        has_ad, ad_bases = 0, 0
    # 注意：为保证 R2 是 R1 插入片段的反补（同一分子），重新生成——见下
    return s1, s2, has_ad, ad_bases


def rng_dna(n, rng):
    return "".join(rng.choices("ACGT", k=n))


def build_pair(idx, level_frac, rng):
    """一条物理片段：insert 序列同时产生 R1（正向）与 R2（反向互补）。"""
    if rng.random() < level_frac:
        ins_len = rng.randint(67, 95)
        insert = rng_dna(ins_len, rng)
        ad1 = ADAPTER_R1[: READ_LEN - ins_len]
        ad2 = ADAPTER_R2[: READ_LEN - ins_len]
        r1 = insert + ad1
        r2 = revcomp(insert) + ad2
        has_ad, ad_bases = 1, READ_LEN - ins_len
    else:
        ins_len = rng.randint(120, 280)
        insert = rng_dna(ins_len, rng)
        r1 = insert[:READ_LEN]
        r2 = revcomp(insert)[:READ_LEN]
        has_ad, ad_bases = 0, 0
    return r1, r2, has_ad, ad_bases


def main():
    rng = random.Random(SEED)
    summary = {"seed": SEED, "n_pairs": N_PAIRS, "read_len": READ_LEN,
               "error_rate": ERROR_RATE,
               "adapters": {"R1_3p": ADAPTER_R1, "R2_3p": ADAPTER_R2},
               "levels": {}}
    truth_rows = []

    for tag, frac in LEVELS:
        r1_path = os.path.join(BASE, "grad%s_R1.fq.gz" % tag)
        r2_path = os.path.join(BASE, "grad%s_R2.fq.gz" % tag)
        n_ad_pairs = 0
        ad_bases_total = 0
        ad_reads = 0
        remnant_hist = {}
        with gzip.open(r1_path, "wt", compresslevel=6) as f1, \
             gzip.open(r2_path, "wt", compresslevel=6) as f2:
            for i in range(N_PAIRS):
                r1, r2, has_ad, ad_bases = build_pair(i, frac, rng)
                q1 = mutate(r1, rng)
                q2 = mutate(r2, rng)
                # PE 惯例：R1/R2 共享同一 read name（cutadapt 校验配对时要求一致）
                name = "grad%s:%06d" % (tag, i)
                f1.write("@%s\n%s\n+\n%s\n" % (name, q1[0], q1[1]))
                f2.write("@%s\n%s\n+\n%s\n" % (name, q2[0], q2[1]))
                truth_rows.append((name, "R1", has_ad, ad_bases, len(r1)))
                truth_rows.append((name, "R2", has_ad, ad_bases, len(r2)))
                if has_ad:
                    n_ad_pairs += 1
                    ad_reads += 2
                    ad_bases_total += 2 * ad_bases
                    remnant_hist[ad_bases] = remnant_hist.get(ad_bases, 0) + 2

        summary["levels"][tag] = {
            "fraction_pairs": frac,
            "pairs": N_PAIRS,
            "read_through_pairs": n_ad_pairs,
            "adapter_reads": ad_reads,
            "adapter_reads_pct": round(100.0 * ad_reads / (2 * N_PAIRS), 2),
            "adapter_bases_total": ad_bases_total,
            "mean_adapter_bases_per_adapter_read":
                round(ad_bases_total / ad_reads, 2) if ad_reads else 0.0,
            "bases_total": 2 * N_PAIRS * READ_LEN,
        }
        print("level %s: read-through pairs %d/%d, adapter reads %.2f%%, "
              "adapter bases %d" % (tag, n_ad_pairs, N_PAIRS,
                                    summary["levels"][tag]["adapter_reads_pct"],
                                    ad_bases_total))

    with gzip.open(os.path.join(BASE, "truth.tsv.gz"), "wt") as f:
        f.write("read_id\tmate\thas_adapter\tadapter_bases\tread_len\n")
        for rid, mate, has_ad, ad_bases, rlen in truth_rows:
            f.write("%s\t%s\t%d\t%d\t%d\n" % (rid, mate, has_ad, ad_bases, rlen))

    with open(os.path.join(BASE, "truth_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print("wrote truth.tsv.gz (%d rows) and truth_summary.json" % len(truth_rows))


if __name__ == "__main__":
    main()
