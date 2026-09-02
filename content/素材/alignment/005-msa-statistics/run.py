#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
005 msa-statistics 真跑脚本：严格复现 bioSkills msa-statistics SKILL.md 给出的统计函数。
输入为同目录 alignment.fasta（构造的小蛋白 MSA，6 条 × 40 列，含保守核心/可变区/空位区）。

产出：
  - repro_transcript.txt   执行 SKILL.md 函数的真实 stdout
  - msa_statistics_data.json  供 make_figs.py 出图
"""
import os, json, math
from Bio import AlignIO
from Bio.Align import substitution_matrices
from collections import Counter

BASE = os.path.dirname(os.path.abspath(__file__))
alignment = AlignIO.read(os.path.join(BASE, "alignment.fasta"), "fasta")

out = []
def log(s=""):
    out.append(str(s)); print(s)

ROBINSON_BACKGROUND = {
    'A': 0.0780, 'R': 0.0512, 'N': 0.0427, 'D': 0.0530, 'C': 0.0193,
    'Q': 0.0419, 'E': 0.0629, 'G': 0.0738, 'H': 0.0224, 'I': 0.0526,
    'L': 0.0922, 'K': 0.0596, 'M': 0.0224, 'F': 0.0399, 'P': 0.0508,
    'S': 0.0712, 'T': 0.0584, 'W': 0.0133, 'Y': 0.0327, 'V': 0.0653,
}

# ============================================================
# 1. 成对一致性（SKILL.md "Calculate Identity Between Two Sequences" + PID 四定义）
# ============================================================
log("="*70); log("1. Pairwise identity (PID1-4) on sp1 vs sp2"); log("="*70)
def pairwise_identity(seq1, seq2, method='pid1'):
    matches = sum(a == b and a != '-' for a, b in zip(seq1, seq2))
    if method == 'pid1':
        denom = sum(a != '-' or b != '-' for a, b in zip(seq1, seq2))
    elif method == 'pid2':
        denom = sum(a != '-' and b != '-' for a, b in zip(seq1, seq2))
    elif method == 'pid3':
        denom = min(len(seq1.replace('-', '')), len(seq2.replace('-', '')))
    elif method == 'pid4':
        denom = (len(seq1.replace('-', '')) + len(seq2.replace('-', ''))) / 2
    return matches / denom if denom > 0 else 0

seq1, seq2 = str(alignment[0].seq), str(alignment[1].seq)
pid = {m: round(pairwise_identity(seq1, seq2, m) * 100, 1) for m in ['pid1','pid2','pid3','pid4']}
for m in ['pid1','pid2','pid3','pid4']:
    log(f"{m}: {pid[m]}%")

# 全 N×N 身份矩阵（pid2，仅残基对）
log("="*70); log("2. Identity matrix (PID2) for all 6 sequences"); log("="*70)
ids = [r.id for r in alignment]
seqs = [str(r.seq) for r in alignment]
idmat = {}
for i, a in enumerate(seqs):
    for j, b in enumerate(seqs):
        idmat.setdefault(ids[i], {})[ids[j]] = round(pairwise_identity(a, b, 'pid2') * 100, 1)
for i in ids:
    log(" ".join(f"{idmat[i][j]:5.1f}" for j in ids))

# ============================================================
# 3. 保守度（SKILL.md "Conservation Score"）
# ============================================================
log("="*70); log("3. Per-column conservation + average + profile(window=10)"); log("="*70)
def column_conservation(alignment, col_idx, ignore_gaps=True):
    column = alignment[:, col_idx]
    if ignore_gaps:
        column = column.replace('-', '')
    if not column:
        return 0.0
    counts = Counter(column)
    return counts.most_common(1)[0][1] / len(column)
def average_conservation(alignment, ignore_gaps=True):
    scores = [column_conservation(alignment, i, ignore_gaps) for i in range(alignment.get_alignment_length())]
    return sum(scores) / len(scores)
def conservation_profile(alignment, window=10):
    prof = []
    L = alignment.get_alignment_length()
    for i in range(L):
        s = max(0, i - window // 2); e = min(L, i + window // 2)
        prof.append(sum(column_conservation(alignment, j) for j in range(s, e)) / (e - s))
    return prof
avg_cons = average_conservation(alignment)
log(f"Average conservation: {avg_cons*100:.1f}%")
profile = conservation_profile(alignment, window=10)

# ============================================================
# 4. 信息量（SKILL.md "Information Content" + Shannon entropy）
# ============================================================
log("="*70); log("4. Shannon entropy + KL information content per column"); log("="*70)
def shannon_entropy(column, ignore_gaps=True):
    if ignore_gaps: column = column.replace('-', '')
    if not column: return 0.0
    counts = Counter(column); total = len(column); ent = 0.0
    for c in counts.values():
        p = c / total
        if p > 0: ent -= p * math.log2(p)
    return ent
def information_content(column, background, ignore_gaps=True):
    if ignore_gaps: column = column.replace('-', '')
    if not column: return 0.0
    counts = Counter(column); total = len(column)
    return sum((c / total) * math.log2((c / total) / background.get(r, 1e-9))
               for r, c in counts.items() if c > 0)
L = alignment.get_alignment_length()
entropy = [round(shannon_entropy(alignment[:, i]), 3) for i in range(L)]
infoc = [round(information_content(alignment[:, i], ROBINSON_BACKGROUND), 3) for i in range(L)]
log(f"mean entropy = {sum(entropy)/L:.3f} bits; mean IC = {sum(infoc)/L:.3f} bits")

# ============================================================
# 5. 空位统计（SKILL.md "Gap Statistics"）
# ============================================================
log("="*70); log("5. Gap fraction per column"); log("="*70)
gap_profile = [round(alignment[:, i].count('-') / len(alignment), 3) for i in range(L)]
avg_gaps = sum(gap_profile) / L
log(f"Average gap fraction: {avg_gaps*100:.1f}%")
gappy_cols = [i for i in range(L) if gap_profile[i] >= 0.5]
log(f"Columns with gap fraction >= 0.5: {gappy_cols}")

# ============================================================
# 6. 比对质量（SKILL.md "Alignment Quality Metrics"）
# ============================================================
log("="*70); log("6. alignment_score + sum_of_pairs (BLOSUM62)"); log("="*70)
def alignment_score(alignment, match=1, mismatch=-1, gap=-2):
    total = 0
    for ci in range(alignment.get_alignment_length()):
        col = alignment[:, ci]
        for i, c1 in enumerate(col):
            for c2 in col[i+1:]:
                if c1 == '-' or c2 == '-': total += gap
                elif c1 == c2: total += match
                else: total += mismatch
    return total
def sum_of_pairs(alignment, substitution_matrix=None):
    if substitution_matrix is None:
        substitution_matrix = substitution_matrices.load('BLOSUM62')
    total = 0.0
    for ci in range(alignment.get_alignment_length()):
        col = alignment[:, ci]
        for i, c1 in enumerate(col):
            for c2 in col[i+1:]:
                if c1 == '-' or c2 == '-': continue
                try: total += substitution_matrix[c1, c2]
                except (KeyError, IndexError): continue
    return total
ascore = alignment_score(alignment)
sp_score = sum_of_pairs(alignment)
log(f"alignment_score = {ascore}")
log(f"sum_of_pairs (BLOSUM62) = {sp_score:.1f}")

# ============================================================
# 7. 距离矩阵（SKILL.md "Distance Correction Models" -> DistanceCalculator blosum62）
# ============================================================
log("="*70); log("7. DistanceCalculator(blosum62)"); log("="*70)
from Bio.Phylo.TreeConstruction import DistanceCalculator
calc = DistanceCalculator('blosum62')
dm = calc.get_distance(alignment)
dmat = {dm.names[i]: {dm.names[j]: round(float(dm[i][j]), 4) for j in range(len(dm.names))}
        for i in range(len(dm.names))}
log(str(dm))

# ============================================================
# 写出
# ============================================================
with open(os.path.join(BASE, "repro_transcript.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(out) + "\n")

data = {
    "n_seqs": len(alignment),
    "n_cols": L,
    "pid_sp1_sp2": pid,
    "identity_matrix": idmat,
    "ids": ids,
    "average_conservation_pct": round(avg_cons * 100, 1),
    "conservation_profile": [round(x, 3) for x in profile],
    "entropy": entropy,
    "info_content": infoc,
    "gap_profile": gap_profile,
    "average_gap_fraction_pct": round(avg_gaps * 100, 1),
    "gappy_cols": gappy_cols,
    "alignment_score": ascore,
    "sum_of_pairs": round(sp_score, 1),
    "distance_matrix": dmat,
}
with open(os.path.join(BASE, "msa_statistics_data.json"), "w", encoding="utf-8") as f:
    json.dump(data, f, indent=2)
log("\nWROTE msa_statistics_data.json + repro_transcript.txt")
