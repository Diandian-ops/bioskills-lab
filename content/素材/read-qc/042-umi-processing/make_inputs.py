#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
042 umi-processing 模拟数据生成（纯 stdlib，确定性 seed=42）。

生成结构：
  - 参考：2 条 contig（chrS1/chrS2），各 20000 bp 随机 ACGT。
  - 原始分子：6000 个，各带 12 nt 随机 UMI（4^12 = 16,777,216 空间）。
  - PCR 重复族：每分子扩增 1~8 份拷贝（加权抽样），模拟读取深度。
  - UMI 错误：第 2 份起的拷贝以 8% 概率发生 1 碱基 UMI 突变
    （directional 法的 count-gradient 规则应把这些 1-off 子代折回母本）。
  - 同坐标不同 UMI 干扰对：300 对分子共享完全相同的 (contig, start)——
    其中 150 对 UMI 汉明距离=1 且家族规模 3~4（count 相近，directional 不合并、
    cluster 按单连通会合并、unique 精确匹配不合并）；另外 150 对 UMI 汉明距离>=4。
  - N 碱基干扰：1% 读对的 UMI 含 1 个 N 碱基（测 extract 对含 N UMI 的处置）。
  - 测序错误：基因组部分每碱基 1% 替换（UMI 部分除设计的 N 外保持干净）。
  - 所有分子固定 + 链（简化，不影响 (坐标+UMI) 去重逻辑）。

输出（同目录）：ref.fa / R1.fq.gz / R2.fq.gz / molecules.tsv / truth.json
读段命名：mol{i:05d}c{j}（extract 后 umi_tools 会把 UMI 以 _UMI 追加到读名）。
"""
import gzip
import json
import os
import random

random.seed(42)
BASE = os.path.dirname(os.path.abspath(__file__))

CONTIGS = [("chrS1", 20000), ("chrS2", 20000)]
FRAG_LEN = 150          # 每分子片段长
R1_GENOMIC = 88         # R1 = UMI(12) + 88 nt 基因组 = 100 nt
R2_LEN = 100            # R2 = frag[50:150] 的反向互补
UMI_LEN = 12
N_MOL = 6000
N_PAIRS_SAME_COORD = 300
N_PAIRS_HAMMING1 = 150  # 其中汉明距离=1 的干扰对（cluster 合并诱饵）
UMI_ERR_PROB = 0.08     # 第 2 份起每拷贝发生 1 碱基 UMI 突变的概率
SEQ_ERR_PROB = 0.01     # 基因组部分每碱基替换率
N_FRAC = 0.01           # UMI 含 1 个 N 碱基的读对比例
FAMILY_SIZES = [1, 2, 3, 4, 5, 6, 8]
FAMILY_WEIGHTS = [30, 25, 20, 10, 8, 4, 3]

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def rand_seq(n):
    return "".join(random.choices("ACGT", k=n))


def hamming1(s):
    """返回与 s 汉明距离恰为 1 的序列。"""
    p = random.randrange(len(s))
    alt = random.choice([b for b in "ACGT" if b != s[p]])
    return s[:p] + alt + s[p + 1:]


def hamming(s, t):
    return sum(1 for a, b in zip(s, t) if a != b)


# ---------- 参考序列 ----------
genome = {}
with open(os.path.join(BASE, "ref.fa"), "w") as f:
    for name, ln in CONTIGS:
        seq = rand_seq(ln)
        genome[name] = seq
        f.write(">%s\n" % name)
        for i in range(0, ln, 60):
            f.write(seq[i:i + 60] + "\n")

# ---------- 分子 ----------
molecules = []          # dict: id, contig, start, umi, family_size, pair_group, pair_role, pair_type
pair_group = 0
n_h1_done = 0

# 300 对同坐标干扰分子（成对占 600 个分子名额）
special_coords = []
for _ in range(N_PAIRS_SAME_COORD):
    contig = random.choice(CONTIGS)[0]
    start = random.randrange(0, genome[contig].__len__() - FRAG_LEN - 10)
    special_coords.append((contig, start))

for ci, (contig, start) in enumerate(special_coords):
    pair_group += 1
    umi_a = rand_seq(UMI_LEN)
    if n_h1_done < N_PAIRS_HAMMING1:
        umi_b = hamming1(umi_a)
        ptype = "hamming1"
        n_h1_done += 1
    else:
        while True:
            umi_b = rand_seq(UMI_LEN)
            if hamming(umi_a, umi_b) >= 4:
                break
        ptype = "diverse"
    for role, umi in (("a", umi_a), ("b", umi_b)):
        molecules.append(dict(
            mol_id="mol%05d" % len(molecules), contig=contig, start=start,
            umi=umi, family_size=random.choice([3, 4]),
            pair_group="P%03d" % pair_group, pair_role=role, pair_type=ptype))

# 其余普通分子
for _ in range(N_MOL - len(molecules)):
    contig = random.choice(CONTIGS)[0]
    start = random.randrange(0, genome[contig].__len__() - FRAG_LEN - 10)
    molecules.append(dict(
        mol_id="mol%05d" % len(molecules), contig=contig, start=start,
        umi=rand_seq(UMI_LEN),
        family_size=random.choices(FAMILY_SIZES, weights=FAMILY_WEIGHTS)[0],
        pair_group="", pair_role="", pair_type=""))

random.shuffle(molecules)

# ---------- 生成读对 ----------
r1_out = gzip.open(os.path.join(BASE, "R1.fq.gz"), "wt", newline="\n")
r2_out = gzip.open(os.path.join(BASE, "R2.fq.gz"), "wt", newline="\n")

n_pairs = 0
n_n_umi_pairs = 0
n_umi_err_copies = 0
family_hist = {}

for m in molecules:
    frag = genome[m["contig"]][m["start"]:m["start"] + FRAG_LEN]
    fam = m["family_size"]
    family_hist[fam] = family_hist.get(fam, 0) + 1
    for j in range(fam):
        umi = m["umi"]
        if j > 0 and random.random() < UMI_ERR_PROB:
            umi = hamming1(umi)
            n_umi_err_copies += 1
        # 1% 读对：UMI 内引入 1 个 N 碱基
        has_n = random.random() < N_FRAC
        if has_n:
            p = random.randrange(UMI_LEN)
            umi = umi[:p] + "N" + umi[p + 1:]
            n_n_umi_pairs += 1
        # R1：UMI + 基因组前 88 nt（带测序错误）
        g1 = "".join(b if random.random() >= SEQ_ERR_PROB
                     else random.choice([c for c in "ACGT" if c != b])
                     for b in frag[:R1_GENOMIC])
        r1 = umi + g1
        # R2：frag[50:150] 反向互补（带测序错误）
        g2 = "".join(b if random.random() >= SEQ_ERR_PROB
                     else random.choice([c for c in "ACGT" if c != b])
                     for b in frag[50:FRAG_LEN])
        r2 = rc(g2)
        name = "%sc%d" % (m["mol_id"], j)
        qual = "I" * len(r1)
        r1_out.write("@%s\n%s\n+\n%s\n" % (name, r1, qual))
        r2_out.write("@%s\n%s\n+\n%s\n" % (name, r2, "I" * len(r2)))
        n_pairs += 1

r1_out.close()
r2_out.close()

# ---------- 元信息 ----------
with open(os.path.join(BASE, "molecules.tsv"), "w") as f:
    f.write("mol_id\tcontig\tstart\tumi\tfamily_size\tpair_group\tpair_role\tpair_type\n")
    for m in molecules:
        f.write("%s\t%s\t%d\t%s\t%d\t%s\t%s\t%s\n" % (
            m["mol_id"], m["contig"], m["start"], m["umi"],
            m["family_size"], m["pair_group"], m["pair_role"], m["pair_type"]))

distinct_coords = len(set((m["contig"], m["start"]) for m in molecules))
n_h1_pairs = sum(1 for m in molecules if m["pair_type"] == "hamming1") // 2
truth = dict(
    seed=42, n_molecules=len(molecules), n_pairs=n_pairs,
    n_umi_err_copies=n_umi_err_copies, n_n_umi_pairs=n_n_umi_pairs,
    umi_len=UMI_LEN, family_hist={str(k): v for k, v in sorted(family_hist.items())},
    distinct_coords=distinct_coords, n_pairs_same_coord=N_PAIRS_SAME_COORD,
    n_pairs_hamming1=n_h1_pairs, umi_space=4 ** UMI_LEN,
    umi_err_prob=UMI_ERR_PROB, seq_err_prob=SEQ_ERR_PROB, n_frac=N_FRAC)
with open(os.path.join(BASE, "truth.json"), "w") as f:
    json.dump(truth, f, indent=2)

print("molecules=%d pairs=%d umi_err_copies=%d n_umi_pairs=%d distinct_coords=%d"
      % (len(molecules), n_pairs, n_umi_err_copies, n_n_umi_pairs, distinct_coords))
print("family_hist=%s" % truth["family_hist"])
