#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
009 alignment-io 真跑脚本：严格复现 bioSkills alignment-io SKILL.md 的读写/转换/
注解保留函数。输入为同目录 alignment.fasta（构造的小 DNA MSA，4 条 × 20 列）。
不依赖 pyhmmer（SKILL.md 流式部分未实测）。

产出：
  - repro_transcript.txt   执行 SKILL.md 函数的真实 stdout
  - alignment_io_data.json  供 make_figs.py 出图
"""
import os, json, tempfile
from Bio import AlignIO
from Bio import Align
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

BASE = os.path.dirname(os.path.abspath(__file__))
alignment = AlignIO.read(os.path.join(BASE, "alignment.fasta"), "fasta")
for _r in alignment:
    _r.annotations['molecule_type'] = 'DNA'

out = []
def log(s=""):
    out.append(str(s)); print(s)

# ============================================================
# 1. 读 / 写 / 转换（SKILL.md "Reading/Writing/Format Conversion"）
# ============================================================
log("="*70); log("1. read / write / convert between formats"); log("="*70)
for fmt in ('fasta', 'clustal', 'phylip-relaxed', 'stockholm', 'nexus'):
    p = os.path.join(BASE, f"demo.{fmt.split('-')[0]}")
    AlignIO.write(alignment, p, fmt)
    back = AlignIO.read(p, fmt)
    log(f"write+read {fmt:14s}: {len(back)} seqs, {back.get_alignment_length()} cols")
# 直接一步转换
c1 = AlignIO.convert(os.path.join(BASE, "demo.clustal"), 'clustal',
                     os.path.join(BASE, "demo_from_clu.phy"), 'phylip-relaxed')
log(f"AlignIO.convert clustal->phylip-relaxed: {c1} alignment(s)")
# nexus 转换必须带 molecule_type，否则抛 ValueError（真实踩坑）
try:
    AlignIO.convert(os.path.join(BASE, "demo.stockholm"), 'stockholm',
                    os.path.join(BASE, "demo_from_sto.nex"), 'nexus')
    log("AlignIO.convert stockholm->nexus (no molecule_type): OK")
except ValueError as e:
    log(f"AlignIO.convert stockholm->nexus (no molecule_type): FAILED -> {e}")
# 带字母表指定
c3 = AlignIO.convert(os.path.join(BASE, "demo.stockholm"), 'stockholm',
                     os.path.join(BASE, "demo_from_sto2.nex"), 'nexus', molecule_type='DNA')
log(f"AlignIO.convert stockholm->nexus (molecule_type=DNA): {c3} alignment(s)")

# ============================================================
# 2. 格式支持实测（SKILL.md "Format Coverage Map"）
# ============================================================
log("="*70); log("2. format read/write support test (Bio.AlignIO)"); log("="*70)
fmts = ['fasta', 'clustal', 'phylip', 'phylip-relaxed', 'phylip-sequential', 'stockholm', 'nexus']
support = {}
for fmt in fmts:
    rec = {'read': False, 'write': False}
    try:
        tmp = os.path.join(BASE, f"_t.{fmt.replace('-', '')}")
        AlignIO.write(alignment, tmp, fmt); rec['write'] = True
        AlignIO.read(tmp, fmt); rec['read'] = True
        os.remove(tmp)
    except Exception as e:
        log(f"  {fmt}: {type(e).__name__}: {e}")
    support[fmt] = rec
    log(f"{fmt:16s} read={rec['read']} write={rec['write']}")

# ============================================================
# 3. Stockholm 注解保留（SKILL.md "Stockholm Format Annotations" + round-trip 坑）
# ============================================================
log("="*70); log("3. Stockholm annotations: kept vs lost after FASTA export"); log("="*70)
# 构造带注解的 Stockholm MSA
aln2 = MultipleSeqAlignment([SeqRecord(Seq(str(r.seq)), id=r.id, description=r.description) for r in alignment])
L = aln2.get_alignment_length()
aln2.column_annotations['secondary_structure'] = '(' * 10 + ')' * 10   # SS_cons, 长度=L
for r in aln2:
    r.annotations['source'] = 'demo'
sto_path = os.path.join(BASE, "annotated.sto")
AlignIO.write(aln2, sto_path, 'stockholm')
# 读回（Stockholm 主）
sto_back = AlignIO.read(sto_path, 'stockholm')
ss_cons = sto_back.column_annotations.get('secondary_structure')
gs_in_sto = any('source' in r.annotations for r in sto_back)
log(f"Stockholm read-back: SS_cons = {ss_cons}; per-seq GS(source) present = {gs_in_sto}")
kept_in_sto = ss_cons is not None
# 导出 FASTA（注解应丢失）
fa_path = os.path.join(BASE, "annotated.fasta")
AlignIO.write(aln2, fa_path, 'fasta')
fa_back = AlignIO.read(fa_path, 'fasta')
lost_ss = 'secondary_structure' not in fa_back.column_annotations
lost_gs = not any('source' in r.annotations for r in fa_back)
log(f"FASTA read-back: SS_cons present = {not lost_ss} (expected False); GS present = {not lost_gs} (expected False)")
# 现代 Bio.Align API
clu_path = os.path.join(BASE, "demo.clustal")
mod = Align.read(clu_path, 'clustal')
log(f"Bio.Align.Align.read(clustal): {len(mod)} seqs, {mod.length} cols")
Align.write(mod, os.path.join(BASE, "demo_modern.fasta"), 'fasta')
log("Bio.Align.Align.write(fasta): ok")

# ============================================================
# 4. 程序化构造比对（SKILL.md "Creating Alignments Programmatically"）
# ============================================================
log("="*70); log("4. create MultipleSeqAlignment from code"); log("="*70)
recs = [SeqRecord(Seq('ACTGACTGACTG'), id='seq1'),
        SeqRecord(Seq('ACTGACT-ACTG'), id='seq2'),
        SeqRecord(Seq('ACTG-CTGACTG'), id='seq3')]
built = MultipleSeqAlignment(recs)
AlignIO.write(built, os.path.join(BASE, "new_alignment.fasta"), 'fasta')
log(f"built alignment: {len(built)} seqs, {built.get_alignment_length()} cols")

# ============================================================
# 写出
# ============================================================
with open(os.path.join(BASE, "repro_transcript.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(out) + "\n")

data = {
    "n_seqs": len(alignment),
    "n_cols": L,
    "format_support": support,
    "annotation": {
        "stockholm_keeps_ss": bool(kept_in_sto),
        "fasta_loses_ss": bool(lost_ss),
        "fasta_loses_gs": bool(lost_gs),
        "ss_cons": ss_cons,
    },
}
with open(os.path.join(BASE, "alignment_io_data.json"), "w", encoding="utf-8") as f:
    json.dump(data, f, indent=2)
log("\nWROTE alignment_io_data.json + repro_transcript.txt")
