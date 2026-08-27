#!/usr/bin/env python3
"""Real reproduction of bioSkills alignment/msa-parsing.

Reads the repo-bundled sample alignment (CLUSTAL), then exercises the
skill's documented code patterns: loading, column/gap analysis, consensus,
filtering, coordinate mapping, Henikoff weights, Neff, and MI/APC.

Input : sample_alignment.aln (CLUSTAL, 4 seqs x 21 cols) from
        content/库/bioSkills/alignment/alignment-io/examples/
Run   : python run_msa_parsing.py
"""
import os
from collections import Counter

import numpy as np

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

HERE = os.path.dirname(os.path.abspath(__file__))
ALN = os.path.join(HERE, "sample_alignment.aln")

print("=" * 70)
print("MSA PARSING — real reproduction (bioSkills alignment/msa-parsing)")
print("=" * 70)

# --- Loading -------------------------------------------------------------
# The repo sample is CLUSTAL format, so read with format='clustal'
# (the skill's examples hard-code format='fasta' for a .fasta input).
alignment = AlignIO.read(ALN, "clustal")
n_seq = len(alignment)
n_col = alignment.get_alignment_length()
print(f"\n[Loading] {n_seq} sequences, {n_col} columns  (format=clustal)")

# --- Extract sequence information ---------------------------------------
seq_ids = [record.id for record in alignment]
sequences = [str(record.seq) for record in alignment]
print(f"[IDs] {seq_ids}")


def get_sequence_by_id(aln, seq_id):
    for record in aln:
        if record.id == seq_id:
            return record
    return None


target = get_sequence_by_id(alignment, "seq3")
print(f"[get_sequence_by_id seq3] {str(target.seq)}")

# --- Column-wise analysis ------------------------------------------------
def find_conserved_positions(aln, threshold=1.0):
    conserved = []
    for col_idx in range(aln.get_alignment_length()):
        column = aln[:, col_idx]
        counts = Counter(column)
        most_common_char, most_common_count = counts.most_common(1)[0]
        if most_common_char != "-":
            conservation = most_common_count / len(aln)
            if conservation >= threshold:
                conserved.append((col_idx, most_common_char, round(conservation, 3)))
    return conserved


fully = find_conserved_positions(alignment, threshold=1.0)
mostly = find_conserved_positions(alignment, threshold=0.8)
print(f"\n[Conserved] fully (>=100%): {len(fully)} positions")
for c in fully:
    print(f"    col {c[0]:>2}  {c[1]}  (cons={c[2]})")
print(f"[Conserved] mostly (>=80%): {len(mostly)} positions -> "
      f"{[c[0] for c in mostly]}")

# --- Gap analysis --------------------------------------------------------
gap_counts = [(record.id, str(record.seq).count("-")) for record in alignment]
print("\n[Gap per sequence]")
for sid, g in gap_counts:
    print(f"    {sid}: {g} gaps ({g / n_col * 100:.1f}%)")


def gaps_per_column(aln):
    return [aln[:, i].count("-") for i in range(aln.get_alignment_length())]


gap_profile = gaps_per_column(alignment)
print(f"[Gap per column] {gap_profile}")


def find_gappy_columns(aln, threshold=0.5):
    num = len(aln)
    return [i for i in range(aln.get_alignment_length())
            if aln[:, i].count("-") / num >= threshold]


gappy = find_gappy_columns(alignment, threshold=0.5)
print(f"[Gappy columns] (>50% gaps): {gappy}")


def remove_gappy_columns(aln, threshold=0.5):
    num = len(aln)
    keep = [i for i in range(aln.get_alignment_length())
            if aln[:, i].count("-") / num < threshold]
    new_records = []
    for record in aln:
        new_seq = "".join(str(record.seq)[i] for i in keep)
        new_records.append(SeqRecord(Seq(new_seq), id=record.id,
                                     description=record.description))
    return MultipleSeqAlignment(new_records), keep


cleaned, kept_cols = remove_gappy_columns(alignment, threshold=0.5)
print(f"[Remove gappy] kept {len(kept_cols)}/{n_col} columns -> "
      f"{cleaned.get_alignment_length()} cols, {len(cleaned)} seqs")

# --- Consensus sequence --------------------------------------------------
def consensus_sequence(aln, threshold=0.5, gap_char="-", ambiguous="N"):
    consensus = []
    for col_idx in range(aln.get_alignment_length()):
        column = aln[:, col_idx]
        counts = Counter(column)
        most_common_char, most_common_count = counts.most_common(1)[0]
        if most_common_char == gap_char:
            counts.pop(gap_char, None)
            if counts:
                most_common_char, most_common_count = counts.most_common(1)[0]
            else:
                most_common_char = gap_char
        if most_common_count / len(aln) >= threshold:
            consensus.append(most_common_char)
        else:
            consensus.append(ambiguous)
    return "".join(consensus)


cons = consensus_sequence(alignment, threshold=0.5)
print(f"\n[Consensus @0.5] {cons}")
cons70 = consensus_sequence(alignment, threshold=0.7)
print(f"[Consensus @0.7] {cons70}")

# --- Extract ungapped regions from reference ----------------------------
def extract_ungapped_regions(aln, ref_idx=0):
    ref_seq = str(aln[ref_idx].seq)
    ungapped_cols = [i for i, ch in enumerate(ref_seq) if ch != "-"]
    new_records = []
    for record in aln:
        new_seq = "".join(str(record.seq)[i] for i in ungapped_cols)
        new_records.append(SeqRecord(Seq(new_seq), id=record.id,
                                     description=record.description))
    return MultipleSeqAlignment(new_records)


ungapped = extract_ungapped_regions(alignment, ref_idx=0)
print(f"\n[Ungapped from ref=0] {ungapped.get_alignment_length()} cols "
      f"(ref seq1 had 1 gap removed)")

# --- Sequence filtering -------------------------------------------------
import re


def filter_by_id(aln, pattern):
    regex = re.compile(pattern)
    return MultipleSeqAlignment([r for r in aln if regex.search(r.id)])


def filter_by_gap_content(aln, max_gap_fraction=0.1):
    return MultipleSeqAlignment(
        [r for r in aln if str(r.seq).count("-") / len(r.seq) <= max_gap_fraction])


def remove_duplicates(aln):
    seen = set()
    return MultipleSeqAlignment(
        [r for r in aln if not (str(r.seq) in seen or seen.add(str(r.seq)))])


f_seq12 = filter_by_id(alignment, r"seq[12]")
f_lowgap = filter_by_gap_content(alignment, max_gap_fraction=0.05)
print(f"\n[filter_by_id seq[12]] {len(f_seq12)} seqs -> {[r.id for r in f_seq12]}")
print(f"[filter_by_gap <=5%] {len(f_lowgap)} seqs -> {[r.id for r in f_lowgap]}")
print(f"[remove_duplicates] {len(remove_duplicates(alignment))} seqs "
      f"(no dups in sample)")

# --- Vectorized coordinate mapping --------------------------------------
def coordinate_map(record):
    chars = np.frombuffer(str(record.seq).encode("ascii"), dtype=np.uint8)
    is_residue = chars != ord("-")
    seq_to_aln = np.flatnonzero(is_residue)
    aln_to_seq = np.where(is_residue, np.cumsum(is_residue) - 1, -1)
    return seq_to_aln, aln_to_seq


s2a, a2s = coordinate_map(alignment[0])
print(f"\n[Coordinate map seq1] residue positions in aln: {s2a.tolist()}")
print(f"    aln col -> seq pos (gap=-1): {a2s.tolist()}")
print(f"    residue at aln col 5 -> seq pos {a2s[5]} "
      f"('{str(alignment[0].seq)[5]}')")

# --- Henikoff sequence weights ------------------------------------------
def henikoff_weights(aln):
    seq_array = np.array([list(str(r.seq)) for r in aln])
    weights = np.zeros(len(aln))
    for col_idx in range(seq_array.shape[1]):
        residues, inverse, counts = np.unique(
            seq_array[:, col_idx], return_inverse=True, return_counts=True)
        if "-" in residues:
            continue
        weights += 1.0 / (len(residues) * counts[inverse])
    return weights / weights.sum()


hw = henikoff_weights(alignment)
print("\n[Henikoff weights]")
for sid, w in zip(seq_ids, hw):
    print(f"    {sid}: {w:.4f}  (sum={hw.sum():.4f})")

# --- Effective sequence number (Neff) -----------------------------------
def _identity(a, b):
    # identity over columns where both are residues
    match = 0
    total = 0
    for ca, cb in zip(a, b):
        if ca == "-" or cb == "-":
            continue
        total += 1
        if ca == cb:
            match += 1
    return match / total if total else 0.0


def neff(aln, threshold=0.62):
    seqs = [str(r.seq) for r in aln]
    clusters = []  # list of (representative_idx, members)
    for i, s in enumerate(seqs):
        placed = False
        for cl in clusters:
            rep = seqs[cl[0]]
            if _identity(rep, s) >= threshold:
                cl[1].append(i)
                placed = True
                break
        if not placed:
            clusters.append([i, [i]])
    sizes = [len(cl[1]) for cl in clusters]
    return float(sum(1.0 / sz for sz in sizes)), sizes


ne, sizes = neff(alignment, threshold=0.62)
print(f"\n[Neff @id=0.62] Neff={ne:.3f}  (L={n_col}, Neff/L={ne / n_col:.3f})")
print(f"    clusters: {sizes}  (4 seqs -> {len(sizes)} clusters)")

# --- Mutual information with APC ----------------------------------------
def mi_matrix_apc(aln):
    seq_array = np.array([list(str(r.seq)) for r in aln])
    L = seq_array.shape[1]
    # frequency of each residue per column
    cols = []
    for j in range(L):
        uniq, cnt = np.unique(seq_array[:, j], return_counts=True)
        p = cnt / cnt.sum()
        cols.append(dict(zip(uniq.tolist(), p.tolist())))
    mi = np.zeros((L, L))
    for i in range(L):
        for j in range(i + 1, L):
            # joint over positions excluding any gap in either column
            mask = (seq_array[:, i] != "-") & (seq_array[:, j] != "-")
            sub = seq_array[mask]
            if len(sub) < 2:
                continue
            uniq_i, cnt_i = np.unique(sub[:, i], return_counts=True)
            uniq_j, cnt_j = np.unique(sub[:, j], return_counts=True)
            pi = dict(zip(uniq_i.tolist(), cnt_i / cnt_i.sum()))
            pj = dict(zip(uniq_j.tolist(), cnt_j / cnt_j.sum()))
            # joint
            pairs, pcnt = np.unique(list(zip(sub[:, i], sub[:, j])),
                                    axis=0, return_counts=True)
            pij = dict(zip([(a, b) for a, b in pairs], pcnt / pcnt.sum()))
            val = 0.0
            for (a, b), pab in pij.items():
                pa = pi[a]
                pb = pj[b]
                if pa > 0 and pb > 0:
                    val += pab * np.log2(pab / (pa * pb))
            mi[i, j] = mi[j, i] = val
    col_mean = mi.mean(axis=0)
    overall = col_mean.mean()
    apc = np.outer(col_mean, col_mean) / overall if overall > 0 else np.zeros_like(mi)
    return mi, mi - apc


mi_raw, mi_apc_mat = mi_matrix_apc(alignment)
top_raw = np.unravel_index(np.argmax(mi_raw), mi_raw.shape)
top_apc = np.unravel_index(np.argmax(mi_apc_mat), mi_apc_mat.shape)
print(f"\n[MI / APC] L={n_col} cols")
print(f"    NOTE: L={n_col} < 100 -> APC over-correction caveat applies; "
      f"raw MI shown, APC reported for completeness")
print(f"    top raw-MI pair  (col {top_raw[0]},{top_raw[1]}) = {mi_raw[top_raw]:.4f} bits")
print(f"    top APC-MI pair  (col {top_apc[0]},{top_apc[1]}) = {mi_apc_mat[top_apc]:.4f} bits")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
