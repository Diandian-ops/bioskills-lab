#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
004 pairwise-alignment 真跑脚本：严格复现 bioSkills pairwise-alignment SKILL.md
给出的 Bio.Align.PairwiseAligner 用法。输入为同目录 sequences.fasta（含 SKILL.md
示例 DNA 序列 ACCGGTAACGTAG / ACCGTTAACGAAG 及构造的蛋白/局部比对样例）。

产出：
  - repro_transcript.txt   执行 SKILL.md 命令的真实 stdout
  - pairwise_data.json     供 make_figs.py 出图的结构化数据

全部命令源于 SKILL.md：Creating an Aligner / Performing Alignments /
Alignment Counts / Percent Identity: Definitions Matter / Substitution Matrix from
Alignment / Export Alignment to Different Formats / Local & Semiglobal 配置。
"""
import os, json
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq
from Bio import SeqIO

BASE = os.path.dirname(os.path.abspath(__file__))
records = list(SeqIO.parse(os.path.join(BASE, "sequences.fasta"), "fasta"))
seq_by_id = {r.id: str(r.seq) for r in records}

out = []
def log(s=""):
    out.append(str(s))
    print(s)

# ============================================================
# 1. 全局 DNA 比对（SKILL.md "Creating an Aligner" DNA 配置 + "Performing Alignments"）
# ============================================================
log("="*70)
log("1. Global DNA alignment  (PairwiseAligner, mode='global')")
log("   match=2 mismatch=-1 open_gap=-10 extend_gap=-0.5")
log("="*70)
aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1,
                          open_gap_score=-10, extend_gap_score=-0.5)
seq1 = Seq(seq_by_id['seq1'])
seq2 = Seq(seq_by_id['seq2'])
alignments = aligner.align(seq1, seq2)
alignment = alignments[0]
log(f"Found {len(alignments)} optimal alignment(s)")
log(alignment)
log(f"alignment.score      = {alignment.score}")
log(f"alignment.shape      = {alignment.shape}")
log(f"len(alignment)       = {len(alignment)}")
counts = alignment.counts()
log(f"counts.identities    = {counts.identities}")
log(f"counts.mismatches    = {counts.mismatches}")
log(f"counts.gaps          = {counts.gaps}")
total_aligned = counts.identities + counts.mismatches
pid_counts = counts.identities / total_aligned * 100
log(f"Percent identity (counts, PID2-like) = {pid_counts:.1f}%")

# 四个 PID 定义（SKILL.md "Percent Identity: Definitions Matter"，分母不同）
tgt = str(alignment[0, :]); qry = str(alignment[1, :])
def pid_defs(a, b):
    matches = sum(x == y and x != '-' for x, y in zip(a, b))
    p1 = matches / sum(x != '-' or y != '-' for x, y in zip(a, b))
    p2 = matches / sum(x != '-' and y != '-' for x, y in zip(a, b))
    p3 = matches / min(len(a.replace('-', '')), len(b.replace('-', '')))
    p4 = matches / ((len(a.replace('-', '')) + len(b.replace('-', ''))) / 2)
    return p1, p2, p3, p4
pid1, pid2, pid3, pid4 = pid_defs(tgt, qry)
log(f"PID1 (aligned+internal gaps) = {pid1*100:.1f}%")
log(f"PID2 (residue pairs only)    = {pid2*100:.1f}%")
log(f"PID3 (shorter length)        = {pid3*100:.1f}%")
log(f"PID4 (mean length)           = {pid4*100:.1f}%")
spread = (max(pid1, pid2, pid3, pid4) - min(pid1, pid2, pid3, pid4)) * 100
log(f"PID spread on this alignment = {spread:.1f}%")

# 替代矩阵（SKILL.md "Substitution Matrix from Alignment"）
substitutions = alignment.substitutions
log("alignment.substitutions (rows=target, cols=query):")
log(str(substitutions))

# 导出格式（SKILL.md "Export Alignment to Different Formats"）
for fmt in ('fasta', 'clustal', 'psl', 'sam'):
    log(f"format(alignment, '{fmt}'):\n{format(alignment, fmt)}")

# ============================================================
# 2. 蛋白比对（BLOSUM62 + 仿射空位，SKILL.md "Protein Alignment"）
# ============================================================
log("="*70)
log("2. Protein alignment  (BLOSUM62, open=-11 extend=-1)")
log("="*70)
blosum62 = substitution_matrices.load('BLOSUM62')
paligner = PairwiseAligner(mode='global', substitution_matrix=blosum62,
                           open_gap_score=-11, extend_gap_score=-1)
pseq1 = Seq(seq_by_id['protA']); pseq2 = Seq(seq_by_id['protB'])
pal = paligner.align(pseq1, pseq2)[0]
log(pal)
pcounts = pal.counts()
log(f"protein score={pal.score} identities={pcounts.identities} "
    f"mismatches={pcounts.mismatches} gaps={pcounts.gaps}")

# ============================================================
# 3. 局部 vs 全局 vs 半全局（SKILL.md "Local" / "Semiglobal" 配置）
#    同一对序列：两端发散、中间保守 -> 模式显著影响 score
# ============================================================
log("="*70)
log("3. Mode comparison on localA/localB (divergent flanks, conserved core)")
log("="*70)
la = Seq(seq_by_id['localA']); lb = Seq(seq_by_id['localB'])

g_mode = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1,
                         open_gap_score=-10, extend_gap_score=-0.5)
g_score = g_mode.score(la, lb)

l_mode = PairwiseAligner(mode='local', match_score=2, mismatch_score=-1,
                         open_gap_score=-10, extend_gap_score=-0.5)
l_score = l_mode.score(la, lb)

s_mode = PairwiseAligner(mode='global')
s_mode.end_gap_score = 0.0   # 双端自由空位 = 半全局
s_score = s_mode.score(la, lb)
log(f"global    score = {g_score}")
log(f"local     score = {l_score}")
log(f"semiglobal score = {s_score}")

# ============================================================
# 写出
# ============================================================
with open(os.path.join(BASE, "repro_transcript.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(out) + "\n")

data = {
    "dna_global": {
        "score": float(alignment.score),
        "identities": int(counts.identities),
        "mismatches": int(counts.mismatches),
        "gaps": int(counts.gaps),
        "alignment_length": int(alignment.shape[1]),
        "pid_counts_pct": round(pid_counts, 2),
    },
    "pid_defs_pct": {
        "pid1": round(pid1 * 100, 2),
        "pid2": round(pid2 * 100, 2),
        "pid3": round(pid3 * 100, 2),
        "pid4": round(pid4 * 100, 2),
    },
    "pid_spread_pct": round(spread, 2),
    "protein": {
        "score": float(pal.score),
        "identities": int(pcounts.identities),
        "mismatches": int(pcounts.mismatches),
        "gaps": int(pcounts.gaps),
        "matrix": "BLOSUM62",
    },
    "mode_scores": {
        "global": float(g_score),
        "local": float(l_score),
        "semiglobal": float(s_score),
    },
}
with open(os.path.join(BASE, "pairwise_data.json"), "w", encoding="utf-8") as f:
    json.dump(data, f, indent=2)
log("\nWROTE pairwise_data.json + repro_transcript.txt")
