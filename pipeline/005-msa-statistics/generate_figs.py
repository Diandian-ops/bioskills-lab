'''Generate real-data figures for the msa-statistics skill reproduction.
Uses ONLY the skill's own functions (shannon_entropy, information_content)
plus the WRONG uniform-background variant to demonstrate the pitfall.'''
from Bio import AlignIO
from collections import Counter
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.expanduser('~/.workbuddy/skills/bioSkills-figure-quality'))
from fig_quality import C_WRONG, C_CORRECT, C_DATA, C_MUTED

# ---- skill's own definitions (verbatim from entropy_analysis.py) ----
ROBINSON_BACKGROUND = {
    'A': 0.0780, 'R': 0.0512, 'N': 0.0427, 'D': 0.0530, 'C': 0.0193,
    'Q': 0.0419, 'E': 0.0629, 'G': 0.0738, 'H': 0.0224, 'I': 0.0526,
    'L': 0.0922, 'K': 0.0596, 'M': 0.0224, 'F': 0.0399, 'P': 0.0508,
    'S': 0.0712, 'T': 0.0584, 'W': 0.0133, 'Y': 0.0327, 'V': 0.0653,
}

def shannon_entropy(column, ignore_gaps=True):
    if ignore_gaps:
        column = column.replace('-', '')
    if not column:
        return 0.0
    counts = Counter(column)
    total = len(column)
    entropy = 0.0
    for count in counts.values():
        p = count / total
        if p > 0:
            entropy -= p * math.log2(p)
    return entropy

def information_content(column, background, ignore_gaps=True):
    if ignore_gaps:
        column = column.replace('-', '')
    if not column:
        return 0.0
    counts = Counter(column)
    total = len(column)
    ic = 0.0
    for residue, count in counts.items():
        observed = count / total
        expected = background.get(residue, 1e-9)
        if observed > 0:
            ic += observed * math.log2(observed / expected)
    return ic

aln = AlignIO.read('alignment.fasta', 'fasta')
L = aln.get_alignment_length()

H = [shannon_entropy(str(aln[:, i])) for i in range(L)]
IC_rob = [information_content(str(aln[:, i]), ROBINSON_BACKGROUND) for i in range(L)]
UNIF = {a: 1/20 for a in ROBINSON_BACKGROUND}
IC_uni = [information_content(str(aln[:, i]), UNIF) for i in range(L)]

print(f"MSA: {len(aln)} seqs x {L} cols")
print(f"avg Shannon H = {sum(H)/L:.3f} bits (max possible = log2(20) = {math.log2(20):.3f})")
print(f"avg IC (Robinson) = {sum(IC_rob)/L:.3f} bits")
print(f"avg IC (uniform)  = {sum(IC_uni)/L:.3f} bits")
print(f"fully conserved cols (H==0): {sum(1 for h in H if h==0)}")

# ================= FIGURE =================
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 7.5))
fig.patch.set_facecolor('white')
for ax in (ax1, ax2):
    ax.set_facecolor('white')

x = np.arange(L)
ax1.plot(x, H, color=C_DATA, lw=1.3, label='Shannon entropy H')
ax1.plot(x, IC_rob, color=C_WRONG, lw=1.3, label='Info content (KL vs Robinson)')
ax1.axhline(math.log2(20), color=C_MUTED, ls='--', lw=0.8, label='max H = log2(20) = 4.32')
ax1.set_title('A. Real globin MSA — per-column entropy & information content', fontsize=10)
ax1.set_xlabel('alignment column', fontsize=9)
ax1.set_ylabel('bits', fontsize=9)
ax1.legend(fontsize=8, loc='upper right')
ax1.set_ylim(-0.2, 7)

# Bottom: the pitfall — fully-conserved column of different residues
residues = ['Leu(L)', 'Ala(A)', 'Gly(G)', 'Trp(W)']
bg = [ROBINSON_BACKGROUND[r[0]] for r in residues]
ic_uni_val = [math.log2(20)] * 4          # uniform: all identical
ic_rob_val = [math.log2(1/b) for b in bg]  # Robinson KL: differs
xpos = np.arange(len(residues))
w = 0.38
# 语义配色：uniform=砖红(WRONG) / Robinson=青瓷绿(correct)，与 004 一致
ax2.bar(xpos - w/2, ic_uni_val, w, color=C_WRONG, label='Uniform background (WRONG for protein)')
ax2.bar(xpos + w/2, ic_rob_val, w, color=C_CORRECT, label='Robinson 1991 background (correct)')
for i, v in enumerate(ic_rob_val):
    ax2.text(xpos[i]+w/2, v+0.08, f'{v:.2f}', ha='center', fontsize=8, color=C_CORRECT)
ax2.set_xticks(xpos)
ax2.set_xticklabels(residues, fontsize=9)
ax2.set_ylabel('IC for a fully-conserved column (bits)', fontsize=9)
ax2.set_title('B. Pitfall: protein IC must use empirical (Robinson) background, not uniform', fontsize=10)
ax2.legend(fontsize=8, loc='upper left')
ax2.set_ylim(0, 7)

plt.tight_layout()
plt.savefig('005-fig.png', dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
print('saved 005-fig.png')
