#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
044 rnaseq-qc 素材生成：合成小型转录组 + dUTP 风格单端 reads（固定随机种子，可复现）。
产出（同目录）：
  transcripts.fa   24 条模拟转录本（长度 500-3000 nt，ACGT 随机序列，无同源序列）
  truth.tsv        每条转录本：长度、分子占比真值 frac、期望 TPM(=frac*1e6)、期望 reads 数
  reads_se.fq.gz   单端 100 bp reads 共 80000 条，1% 替换错误
方向设定：dUTP / fr-firststrand 的 SE 口径 —— 读段取转录本正义链窗口的反向互补
（读段位于反义链），对应 SKILL.md 文库表中的 salmon SR。对照库型由 salmon 侧
（-l A / -l IU / -l SF）切换，数据只生成一份。
采样口径：reads 按片段占比 frac*length 分配（长转录本产生更多片段，位置均匀），
因此期望 TPM = frac * 1e6，与 salmon TPM 直接可比。
"""
import gzip
import os
import random

BASE = os.path.dirname(os.path.abspath(__file__))
SEED_TX = 44
SEED_READS = 44044
N_TX = 24
READ_LEN = 100
N_READS = 80000
ERR_RATE = 0.01          # 每碱基替换错误率
FRAC_LO, FRAC_HI = -3.5, -0.5   # 分子占比真值 log10 均匀取样区间（跨 3 个数量级）

rng = random.Random(SEED_TX)
COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


# ---------- 1) 转录本序列与表达梯度真值 ----------
lengths = [rng.randrange(500, 3001) for _ in range(N_TX)]
raw = sorted((10 ** rng.uniform(FRAC_LO, FRAC_HI) for _ in range(N_TX)), reverse=True)
tot = sum(raw)
fracs = [v / tot for v in raw]

txs = []
for i in range(N_TX):
    tid = "tx%04d" % (i + 1)
    seq = "".join(rng.choice("ACGT") for _ in range(lengths[i]))
    txs.append((tid, seq))

# 期望 reads 按片段占比 frac*L 归一分配（长转录本产生更多片段）
wsum = sum(fracs[i] * lengths[i] for i in range(N_TX))
exp_reads = [N_READS * fracs[i] * lengths[i] / wsum for i in range(N_TX)]
with open(os.path.join(BASE, "transcripts.fa"), "w") as fa, \
     open(os.path.join(BASE, "truth.tsv"), "w") as tsv:
    tsv.write("tx_id\tlength\tfrac\ttpm_expected\treads_expected\n")
    for i in range(N_TX):
        fa.write(">%s len=%d frac=%.6g\n%s\n"
                 % (txs[i][0], lengths[i], fracs[i], txs[i][1]))
        tsv.write("%s\t%d\t%.8g\t%.2f\t%.2f\n"
                  % (txs[i][0], lengths[i], fracs[i],
                     fracs[i] * 1e6, exp_reads[i]))

# ---------- 2) dUTP 风格 SE reads ----------
rng2 = random.Random(SEED_READS)
alloc = [int(x) for x in exp_reads]
left = N_READS - sum(alloc)
order = sorted(range(N_TX), key=lambda i: exp_reads[i] - int(exp_reads[i]), reverse=True)
for k in range(left):
    alloc[order[k % N_TX]] += 1

recs = []
rid = 0
for i in range(N_TX):
    tid, seq = txs[i]
    L = len(seq)
    for _ in range(alloc[i]):
        start = rng2.randrange(0, L - READ_LEN + 1)
        read = rc(seq[start:start + READ_LEN])   # dUTP SE：读段在反义链
        if ERR_RATE > 0:
            bases = []
            for b in read:
                if rng2.random() < ERR_RATE:
                    alt = rng2.choice([c for c in "ACGT" if c != b])
                    bases.append(alt)
                else:
                    bases.append(b)
            read = "".join(bases)
        rid += 1
        recs.append("@read%06d\n%s\n+\n%s\n" % (rid, read, "I" * READ_LEN))
rng2.shuffle(recs)

with gzip.open(os.path.join(BASE, "reads_se.fq.gz"), "wt") as fq:
    fq.write("".join(recs))

print("transcripts: %d (lengths %d-%d nt, total %d bp)"
      % (N_TX, min(lengths), max(lengths), sum(lengths)))
print("frac range: %.3g - %.3g (span %.0fx), expected TPM %.0f - %.0f"
      % (min(fracs), max(fracs), max(fracs) / min(fracs),
         min(fracs) * 1e6, max(fracs) * 1e6))
print("reads: %d SE %d bp, err %.2f, min/max per-tx expected reads %.1f / %.1f"
      % (N_READS, READ_LEN, ERR_RATE, min(exp_reads), max(exp_reads)))
print("orientation: dUTP-style SE (read = reverse complement of transcript)")
print("wrote transcripts.fa, truth.tsv, reads_se.fq.gz")
