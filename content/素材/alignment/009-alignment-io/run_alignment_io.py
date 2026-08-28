"""009 alignment-io 真实试用复现脚本。

忠实复现 bioSkills alignment-io SKILL.md 文档化用法。
样本：仓库 content/库/bioSkills/alignment/alignment-io/examples/sample_alignment.aln
      （CLUSTAL, 4 条序列 × 21 列，含 * 保守标记行）
环境：BioPython 1.88 / managed venv。

覆盖点：
  1. 读取（AlignIO.read / parse / list）
  2. 写出（single / multiple / handle）
  3. 格式转换（AlignIO.convert + molecule_type）
  4. 访问与切片（迭代 / 索引 / 列切片 / 组合切片）
  5. 程序化构建（MultipleSeqAlignment + SeqRecord）
  6. PHYLIP 严格 10 字符截断静默合并 footgun
  7. phylip-relaxed 保留长名
  8. Stockholm 注释保留（GS/GR/GC）与转 FASTA 静默丢弃（round-trip caveat）
  9. A2M/A3M 大小写编码 → 提取 match 列
 10. 现代 Bio.Align API（Align.read → .counts() / .substitutions）
 11. MAF 块坐标转换（maf_to_plus_strand_coords）
 12. pyhmmer.easel.MSAFile 流式（未安装 → 诚实声明）
"""
from __future__ import annotations
import io
import json
from pathlib import Path

from Bio import AlignIO, Align
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

HERE = Path(__file__).parent
SAMPLE = HERE / "sample_alignment.aln"


def section(t: str) -> None:
    print("\n" + "=" * 72)
    print(t)
    print("=" * 72)


def ensure_dna(rec: SeqRecord) -> SeqRecord:
    """Biopython 1.78+ 写出部分格式要求 molecule_type。"""
    if "molecule_type" not in rec.annotations:
        rec.annotations["molecule_type"] = "DNA"
    return rec


def safe_write(alignment, path: Path, fmt: str) -> bool:
    try:
        AlignIO.write(alignment, str(path), fmt)
        return True
    except (ValueError, TypeError) as e:
        print(f"  [write {fmt}] 跳过：{e}")
        return False


# ----------------------------------------------------------------------------
section("1. 读取 single / parse / list")
aln = AlignIO.read(str(SAMPLE), "clustal")
print(f"AlignIO.read(clustal) -> {len(aln)} 条序列, {aln.get_alignment_length()} 列")
print("序列 ID / 长度：")
for rec in aln:
    print(f"  {rec.id}: {len(rec.seq)} bp")

# 多比对文件 parse（用同一文件演示 list(parse)）
alignments = list(AlignIO.parse(str(SAMPLE), "clustal"))
print(f"list(AlignIO.parse) -> {len(alignments)} 个比对")

# ----------------------------------------------------------------------------
section("2. 写出 single / multiple / handle")
safe_write(aln, HERE / "out.fasta", "fasta")
safe_write(aln, HERE / "out.phy", "phylip-relaxed")
# handle 形式
with open(HERE / "out_handle.clustal", "w") as fh:
    AlignIO.write(aln, fh, "clustal")
print("写出 out.fasta / out.phy / out_handle.clustal")

# ----------------------------------------------------------------------------
section("3. 格式转换 AlignIO.convert + molecule_type")
AlignIO.convert(str(SAMPLE), "clustal", str(HERE / "conv.sto"), "stockholm")
AlignIO.convert(str(SAMPLE), "clustal", str(HERE / "conv.nex"), "nexus", molecule_type="DNA")
print("convert: clustal -> stockholm / nexus(带 molecule_type)")

# ----------------------------------------------------------------------------
section("4. 访问与切片")
print("迭代序列：")
for rec in aln:
    print(f"  {rec.id}: {rec.seq}")
print(f"索引 aln[0].id = {aln[0].id}; aln[-1].id = {aln[-1].id}")
col_slice = aln[:, 5:15]
print(f"列切片 aln[:,5:15] -> {col_slice.get_alignment_length()} 列")
col5 = aln[:, 5]
print(f"单列 aln[:,5] -> {col5}")

# ----------------------------------------------------------------------------
section("5. 程序化构建比对")
built = MultipleSeqAlignment([
    SeqRecord(Seq("ACTGACTGACTG"), id="seq1"),
    SeqRecord(Seq("ACTGACT-ACTG"), id="seq2"),
    SeqRecord(Seq("ACTG-CTGACTG"), id="seq3"),
])
print(f"手工构建 -> {len(built)} 条, {built.get_alignment_length()} 列")
safe_write(built, HERE / "built.fasta", "fasta")

# ----------------------------------------------------------------------------
section("6. PHYLIP 严格 10 字符截断行为（Bio 1.88 实测）")
# 6a. 截断但不冲突：两条前缀不同的长名 -> 截断到 10 字符，写成功
distinct = MultipleSeqAlignment([
    ensure_dna(SeqRecord(Seq("ACGTACGTACGTACGT"), id="Homo_sapiens_chr1")),
    ensure_dna(SeqRecord(Seq("ACGTACGTACGTACGT"), id="Mus_musculus_chr1")),
])
AlignIO.write(distinct, str(HERE / "strict_trunc.phy"), "phylip")
trunc_back = AlignIO.read(str(HERE / "strict_trunc.phy"), "phylip")
print(f"截断(不冲突)读回 ID: {[r.id for r in trunc_back]}  (长名被截到 10 字符)")

# 6b. 冲突：两条仅 10 字符前缀相同的长名
#     Bio 1.88 实测：直接抛 ValueError 而非静默合并（防御性，比 skill 描述的"静默"更安全）
collide = MultipleSeqAlignment([
    ensure_dna(SeqRecord(Seq("ACGTACGTACGTACGT"), id="Homo_sapiens_chr1")),
    ensure_dna(SeqRecord(Seq("ACGTACGTACGTACGT"), id="Homo_sapiens_chr2")),
])
try:
    AlignIO.write(collide, str(HERE / "strict_collide.phy"), "phylip")
    cb = AlignIO.read(str(HERE / "strict_collide.phy"), "phylip")
    print(f"冲突情况读回序列数: {len(cb)}  (未报错)")
except ValueError as e:
    print(f"冲突情况 Bio 1.88 直接抛错（非静默合并）: {e}")

# ----------------------------------------------------------------------------
section("7. phylip-relaxed 保留长名")
AlignIO.write(collide, str(HERE / "relaxed.phy"), "phylip-relaxed")
relaxed_back = AlignIO.read(str(HERE / "relaxed.phy"), "phylip-relaxed")
relaxed_ids = [r.id for r in relaxed_back]
print(f"phylip-relaxed 读回 ID: {relaxed_ids}  (长名保留)")

# ----------------------------------------------------------------------------
section("8. Stockholm 注释保留 vs 转 FASTA 静默丢弃")
sto_text = """# STOCKHOLM 1.0
#=GF ID    TEST_FAMILY
seq1    ACGTACG
seq2    ACGTACG
seq3    ACGTACG
#=GS seq1 OS Homo sapiens
#=GR seq1 SS .HHH...
#=GC SS_cons ..HHH..
//
"""
sto_path = HERE / "test.sto"
sto_path.write_text(sto_text)
sto_aln = AlignIO.read(str(sto_path), "stockholm")
gs = {r.id: r.annotations for r in sto_aln}
gr = {r.id: r.letter_annotations.get("secondary_structure") for r in sto_aln}
gc = sto_aln.column_annotations.get("secondary_structure")
print("Stockholm 读回注释：")
print(f"  GS (record.annotations): {gs}")
print(f"  GR (letter_annotations SS): {gr}")
print(f"  GC (column_annotations SS_cons): {gc}")

# 转 FASTA 再读回，检查注释是否保留
AlignIO.write(sto_aln, str(HERE / "sto2fasta.fasta"), "fasta")
fasta_back = AlignIO.read(str(HERE / "sto2fasta.fasta"), "fasta")
fas_gs = any(bool(r.annotations) for r in fasta_back)
fas_gr = any("secondary_structure" in r.letter_annotations for r in fasta_back)
fas_gc = "secondary_structure" in fasta_back.column_annotations
print("Stockholm -> FASTA -> 读回：")
print(f"  GS 保留? {fas_gs}; GR 保留? {fas_gr}; GC 保留? {fas_gc}")
print("  (三者皆 False = 注释被静默丢弃，符合 round-trip caveat)")

# ----------------------------------------------------------------------------
section("9. A2M/A3M 大小写编码 → 提取 match 列")
a2m = MultipleSeqAlignment([
    ensure_dna(SeqRecord(Seq("AC-EFGH"), id="q")),     # 插入列处为 '-'
    ensure_dna(SeqRecord(Seq("ACdEFGH"), id="hit1")),  # 小写 d = insert 列
])
match_only = [
    "".join(c for c in str(r.seq) if c.isupper() or c == "-")
    for r in a2m
]
print("A2M 原始序列：")
for r in a2m:
    print(f"  {r.id}: {r.seq}")
print("提取 match 列（仅大写 + '-'）：")
for r, m in zip(a2m, match_only):
    print(f"  {r.id}: {m}")

# ----------------------------------------------------------------------------
section("10. 现代 Bio.Align API：Align.read / .counts / .substitutions")
mod = Align.read(str(SAMPLE), "clustal")
print(f"Align.read -> type {type(mod).__name__}, {len(mod)} 条, {mod.length} 列")
counts = mod.counts()
print(f".counts() -> type {type(counts).__name__}")
subs = mod.substitutions  # Bio 1.88：属性(Array)，非方法
print(f".substitutions -> type {type(subs).__name__}; 维度 {getattr(subs, 'shape', 'N/A')}")

# ----------------------------------------------------------------------------
section("11. MAF 块坐标转换（maf_to_plus_strand_coords）")
maf_text = """##maf version=1
a score=0
s ref.chr1  0 7 + 100 ACGTACG
s qry.chr2 10 7 -  50 ACGTACG
"""
def maf_to_plus_strand_coords(row_anno):
    if row_anno["strand"] == "-":
        return row_anno["srcSize"] - row_anno["start"] - row_anno["size"]
    return row_anno["start"]

maf_blocks = list(AlignIO.parse(io.StringIO(maf_text), "maf"))
print(f"MAF 解析出 {len(maf_blocks)} 个 block")

# Biopython 1.88 实测：MAF 行的 strand 解析为 int (-1 / 1)，源名在 rec.id
print("逐行注释（注意 strand 类型与源名字段）：")
for block in maf_blocks:
    for rec in block:
        a = rec.annotations
        print(f"  rec.id={rec.id!r} annotations={a}")

# SKILL.md 原函数（按字符串 '-' 比较）——在 Bio 1.88 下会判错
print("SKILL.md 原函数 maf_to_plus_strand_coords（strand == '-'）：")
for block in maf_blocks:
    for rec in block:
        a = rec.annotations
        plus = maf_to_plus_strand_coords(a)
        print(f"  {rec.id}: strand={a['strand']!r} -> plus_start={plus}  "
              f"(期望 minus 行 = {a['srcSize'] - a['start'] - a['size']})")

# 修正版：兼容 int 与 str 的 strand
def maf_to_plus_strand_coords_fixed(row_anno):
    strand = row_anno["strand"]
    if strand in ("-", -1):
        return row_anno["srcSize"] - row_anno["start"] - row_anno["size"]
    return row_anno["start"]

print("修正版（兼容 int/-1）：")
for block in maf_blocks:
    for rec in block:
        a = rec.annotations
        plus = maf_to_plus_strand_coords_fixed(a)
        print(f"  {rec.id}: strand={a['strand']!r} -> plus_start={plus}")

# ----------------------------------------------------------------------------
section("12. pyhmmer.easel.MSAFile 流式（大 Stockholm 数据库）")
try:
    import pyhmmer  # noqa: F401
    print("pyhmmer 已安装，可流式读取 Pfam-A.full / BFD。")
except ImportError:
    print("pyhmmer 未安装 -> 流式读取演示无法真跑（诚实声明）。")
    print("  pip install pyhmmer 后可按 SKILL.md 示例流式处理 Pfam-A.full。")

# ----------------------------------------------------------------------------
section("注释存活矩阵（供出图）")
# 以 test.sto（带 GS/GR/GC）为源，转换到各格式后回读，检测三类注释是否存活
survival = {}
targets = {
    "stockholm": "stockholm",
    "fasta": "fasta",
    "clustal": "clustal",
    "phylip-relaxed": "phylip-relaxed",
    "nexus": "nexus",
}
src = AlignIO.read(str(sto_path), "stockholm")
for r in src:
    ensure_dna(r)
for name, fmt in targets.items():
    outp = HERE / f"_surv_{name}.txt"
    ok = safe_write(src, outp, fmt)
    if not ok:
        survival[name] = {"GS": False, "GR": False, "GC": False, "written": False}
        continue
    try:
        back = AlignIO.read(str(outp), fmt)
        gs = any(bool(r.annotations) for r in back)
        gr = any("secondary_structure" in r.letter_annotations for r in back)
        gc = "secondary_structure" in back.column_annotations
        survival[name] = {"GS": gs, "GR": gr, "GC": gc, "written": True}
    except Exception as e:  # noqa: BLE001
        survival[name] = {"GS": False, "GR": False, "GC": False, "written": False, "err": str(e)}
print(json.dumps(survival, indent=2))
(HERE / "annotation_survival.json").write_text(json.dumps(survival, indent=2))
print("\n009 复现完成。")
