'''pairwise-alignment skill 严格复现
严格使用 skill 自带的方法与代码模式（protein_alignment.py / alignment_from_file.py / SKILL.md）：
  - 蛋白比对标准配置：BLOSUM62 + global + open=-11/extend=-1
  - 读 FASTA（skill 的“Working with SeqRecord Objects”模式）
  - 演示 SKILL.md 点名的两个坑：默认 gap=0 陷阱、PID 四口径
不引入 skill 之外的任何分析。
'''
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

# ---- skill: Working with SeqRecord Objects / alignment_from_file.py ----
records = list(SeqIO.parse('sequences.fasta', 'fasta'))
seq1, seq2 = records[0].seq, records[1].seq
print(f'序列1: {records[0].id}  (长度 {len(seq1)})')
print(f'序列2: {records[1].id}  (长度 {len(seq2)})')

# ===== 坑一演示：PairwiseAligner 默认 gap penalty = 0 =====
# SKILL.md: “PairwiseAligner() with no arguments uses match_score=1, mismatch_score=0,
#            open_gap_score=0, extend_gap_score=0 ... Always specify gap penalties explicitly”
print('\n========== 坑一：默认 gap penalty = 0 的陷阱 ==========')
blosum62 = substitution_matrices.load('BLOSUM62')

aligner_default = PairwiseAligner(mode='global', substitution_matrix=blosum62)  # 不显式设 gap
aligner_recommended = PairwiseAligner(mode='global', substitution_matrix=blosum62,
                                       open_gap_score=-11, extend_gap_score=-1)

a_def = aligner_default.align(seq1, seq2)[0]
a_rec = aligner_recommended.align(seq1, seq2)[0]

c_def = a_def.counts()
c_rec = a_rec.counts()
L_def = a_def.shape[1]   # 对齐长度 = 列数（len(alignment) 返回的是序列数=2，不能用）
L_rec = a_rec.shape[1]
print(f'默认配置(未设gap)   : align_len={L_def:3d}  gaps={c_def.gaps:3d}  score={a_def.score:.1f}')
print(f'推荐配置(-11/-1)    : align_len={L_rec:3d}  gaps={c_rec.gaps:3d}  score={a_rec.score:.1f}')
print(f'默认比推荐多出的 gap : {c_def.gaps - c_rec.gaps}')

# ===== 标准复现：用 skill 的 protein_alignment.py 配置输出比对 =====
print('\n========== 标准蛋白比对 (skill protein_alignment.py 配置) ==========')
print(f'Score: {a_rec.score}\n')
print(a_rec)

# ===== 坑二演示：Percent Identity 四口径 =====
# SKILL.md “Percent Identity: Definitions Matter” 表的精确定义
print('\n========== 坑二：Percent Identity 四口径 ==========')
ident = c_rec.identities
mm = c_rec.mismatches
gaps = c_rec.gaps
L = len(a_rec)                       # alignment length
n_pairs = ident + mm                 # aligned residue pairs (excl gaps)
L_min = min(len(seq1), len(seq2))    # shorter sequence length
L_mean = (len(seq1) + len(seq2)) / 2 # mean sequence length

pid1 = ident / L * 100               # 对齐长度(含gap) 作分母 — gap-aware, 保守
pid2 = ident / n_pairs * 100         # 非gap配对位置 作分母 — 总是最高
pid3 = ident / L_min * 100           # 较短序列长度 作分母 — 长度归一
pid4 = ident / L_mean * 100          # 平均序列长度 作分母 — 与结构相似度最相关

print(f'identities={ident}  mismatches={mm}  gaps={gaps}  align_len={L_rec}')
print(f'PID1 (含gap对齐长) : {pid1:.1f}%')
print(f'PID2 (非gap配对)   : {pid2:.1f}%   <- counts() 近似此口径')
print(f'PID3 (较短序列长)   : {pid3:.1f}%')
print(f'PID4 (平均序列长)   : {pid4:.1f}%')
print(f'四口径极差         : {max(pid1,pid2,pid3,pid4)-min(pid1,pid2,pid3,pid4):.1f} pct')

# 写一份对齐文本供笔记引用
with open('alignment_rec.txt', 'w') as f:
    f.write(str(a_rec))
print('\n[已写出 alignment_rec.txt]')
