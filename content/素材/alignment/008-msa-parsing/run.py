#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
008 msa-parsing 真跑脚本：严格复现 bioSkills msa-parsing SKILL.md 的解析/过滤/统计函数。
输入为同目录 alignment.aln（构造的小 DNA MSA，5 条 × 30 列，含保守核心/可变区/
空位区，且 s1/s2 近重复以演示 Henikoff 权重）。

产出：
  - repro_transcript.txt   执行 SKILL.md 函数的真实 stdout
  - msa_parsing_data.json  供 make_figs.py 出图
"""
import os, json, re
import numpy as np
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter

BASE = os.path.dirname(os.path.abspath(__file__))
alignment = AlignIO.read(os.path.join(BASE, "alignment.aln"), "fasta")

out = []
def log(s=""):
    out.append(str(s)); print(s)

# ============================================================
# 1. 加载与基本信息（SKILL.md "Loading Alignments" / "Extracting Sequence Information"）
# ============================================================
log("="*70); log("1. Load + sequence ids / get by id"); log("="*70)
log(f"{len(alignment)} sequences, {alignment.get_alignment_length()} columns")
seq_ids = [r.id for r in alignment]
log("seq_ids = " + str(seq_ids))
def get_sequence_by_id(alignment, seq_id):
    for r in alignment:
        if r.id == seq_id:
            return r
    return None
t = get_sequence_by_id(alignment, 's1')
log(f"get_sequence_by_id('s1') -> {str(t.seq)}")

# ============================================================
# 2. 保守位点（SKILL.md "Find Conserved Positions"）
# ============================================================
log("="*70); log("2. find_conserved_positions (threshold 1.0 / 0.8)"); log("="*70)
def find_conserved_positions(alignment, threshold=1.0):
    conserved = []
    for ci in range(alignment.get_alignment_length()):
        col = alignment[:, ci]
        counts = Counter(col)
        mc, mc_count = counts.most_common(1)[0]
        if mc != '-' and mc_count / len(alignment) >= threshold:
            conserved.append((ci, mc))
    return conserved
fully = find_conserved_positions(alignment, 1.0)
mostly = find_conserved_positions(alignment, 0.8)
log(f"fully conserved (>=1.0): {len(fully)} columns")
log(f"mostly conserved (>=0.8): {len(mostly)} columns")

# ============================================================
# 3. 空位分析（SKILL.md "Gap Analysis"）
# ============================================================
log("="*70); log("3. per-sequence gaps / per-column gaps / gappy columns"); log("="*70)
gap_counts = [(r.id, str(r.seq).count('-')) for r in alignment]
for sid, g in gap_counts:
    log(f"{sid}: {g} gaps")
def gaps_per_column(alignment):
    return [alignment[:, i].count('-') for i in range(alignment.get_alignment_length())]
gpc = gaps_per_column(alignment)
def find_gappy_columns(alignment, threshold=0.5):
    num = len(alignment); gappy = []
    for ci in range(alignment.get_alignment_length()):
        if alignment[:, ci].count('-') / num >= threshold:
            gappy.append(ci)
    return gappy
gappy_cols = find_gappy_columns(alignment, 0.5)
log(f"gappy columns (>=0.5): {gappy_cols}")
def remove_gappy_columns(alignment, threshold=0.5):
    num = len(alignment); keep = []
    for ci in range(alignment.get_alignment_length()):
        if alignment[:, ci].count('-') / num < threshold:
            keep.append(ci)
    recs = []
    for r in alignment:
        new = ''.join(str(r.seq)[i] for i in keep)
        recs.append(SeqRecord(Seq(new), id=r.id, description=r.description))
    return MultipleSeqAlignment(recs)
cleaned = remove_gappy_columns(alignment, 0.5)
log(f"after remove_gappy_columns: {len(cleaned)} seqs, {cleaned.get_alignment_length()} columns")

# ============================================================
# 4. 一致序列（SKILL.md "Consensus Sequence"）
# ============================================================
log("="*70); log("4. consensus_sequence (threshold 0.5)"); log("="*70)
def consensus_sequence(alignment, threshold=0.5, gap_char='-', ambiguous='N'):
    cons = []
    for ci in range(alignment.get_alignment_length()):
        col = alignment[:, ci]
        counts = Counter(col); mc, mcc = counts.most_common(1)[0]
        if mc == gap_char:
            counts.pop(gap_char, None)
            if counts: mc, mcc = counts.most_common(1)[0]
            else: mc = gap_char
        cons.append(mc if mcc / len(alignment) >= threshold else ambiguous)
    return ''.join(cons)
cons = consensus_sequence(alignment, 0.5)
log(f"consensus = {cons}")

# ============================================================
# 5. 提取无空位区域（SKILL.md "Extract Ungapped Regions"）
# ============================================================
log("="*70); log("5. extract_ungapped_regions(ref_idx=0)"); log("="*70)
def extract_ungapped_regions(alignment, ref_idx=0):
    ref = str(alignment[ref_idx].seq)
    ungapped = [i for i, ch in enumerate(ref) if ch != '-']
    recs = []
    for r in alignment:
        new = ''.join(str(r.seq)[i] for i in ungapped)
        recs.append(SeqRecord(Seq(new), id=r.id, description=r.description))
    return MultipleSeqAlignment(recs)
ungapped = extract_ungapped_regions(alignment, 0)
log(f"ungapped regions: {len(ungapped)} seqs, {ungapped.get_alignment_length()} columns")

# ============================================================
# 6. 过滤（SKILL.md "Sequence Filtering"）
# ============================================================
log("="*70); log("6. filter_by_gap_content / remove_duplicates"); log("="*70)
def filter_by_gap_content(alignment, max_gap_fraction=0.1):
    return MultipleSeqAlignment([r for r in alignment
                                 if str(r.seq).count('-') / len(r.seq) <= max_gap_fraction])
fg = filter_by_gap_content(alignment, 0.1)
log(f"filter_by_gap_content(0.1): {len(fg)} seqs kept")
def remove_duplicates(alignment):
    seen = set()
    return MultipleSeqAlignment([r for r in alignment
                                 if not (str(r.seq) in seen or seen.add(str(r.seq)))])
rd = remove_duplicates(alignment)
log(f"remove_duplicates: {len(rd)} unique seqs")

# ============================================================
# 7. 坐标映射（SKILL.md "Vectorized Coordinate Mapping"）
# ============================================================
log("="*70); log("7. coordinate_map(alignment[0])"); log("="*70)
def coordinate_map(record):
    chars = np.frombuffer(str(record.seq).encode('ascii'), dtype=np.uint8)
    is_res = chars != ord('-')
    seq_to_aln = np.flatnonzero(is_res)
    aln_to_seq = np.where(is_res, np.cumsum(is_res) - 1, -1)
    return seq_to_aln, aln_to_seq
s2a, a2s = coordinate_map(alignment[0])
log(f"residue 10 -> alignment column {int(s2a[10])}")
log(f"alignment column 15 -> residue index {int(a2s[15])}")

# ============================================================
# 8. Henikoff 权重（SKILL.md "Henikoff Sequence Weights"）
# ============================================================
log("="*70); log("8. henikoff_weights"); log("="*70)
def henikoff_weights(alignment):
    arr = np.array([list(str(r.seq)) for r in alignment])
    w = np.zeros(len(alignment))
    for ci in range(arr.shape[1]):
        residues, inverse, counts = np.unique(arr[:, ci], return_inverse=True, return_counts=True)
        if '-' in residues: continue
        w += 1.0 / (len(residues) * counts[inverse])
    return w / w.sum()
hw = henikoff_weights(alignment)
hw_map = {alignment[i].id: round(float(hw[i]), 4) for i in range(len(alignment))}
for k, v in hw_map.items(): log(f"{k}: weight = {v}")
log(f"weight sum = {sum(hw_map.values()):.4f}")
log(f"note: s1/s2 are near-duplicates -> both down-weighted vs unique s3/s4/s5")

# ============================================================
# 写出
# ============================================================
with open(os.path.join(BASE, "repro_transcript.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(out) + "\n")

data = {
    "n_seqs": len(alignment),
    "n_cols": alignment.get_alignment_length(),
    "seq_ids": seq_ids,
    "gaps_per_sequence": {k: g for k, g in gap_counts},
    "gaps_per_column": gpc,
    "gappy_cols": gappy_cols,
    "n_fully_conserved": len(fully),
    "n_mostly_conserved": len(mostly),
    "cleaned_cols": cleaned.get_alignment_length(),
    "consensus": cons,
    "ungapped_cols": ungapped.get_alignment_length(),
    "filtered_kept": len(fg),
    "unique_seqs": len(rd),
    "henikoff_weights": hw_map,
    "coord_sample": {"residue10_to_col": int(s2a[10]), "col15_to_residue": int(a2s[15])},
}
with open(os.path.join(BASE, "msa_parsing_data.json"), "w", encoding="utf-8") as f:
    json.dump(data, f, indent=2)
log("\nWROTE msa_parsing_data.json + repro_transcript.txt")
