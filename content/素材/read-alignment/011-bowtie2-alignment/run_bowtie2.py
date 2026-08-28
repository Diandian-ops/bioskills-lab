#!/usr/bin/env python3
"""011 bowtie2-alignment — 忠实复现 bioSkills SKILL.md 的方法与结论。

覆盖点：
  - bowtie2-build 索引构建（-x 期望 basename 而非文件名，含负向测试）
  - end-to-end(默认) vs --local 的核心生物学决策
  - 适配器污染 reads 下 end-to-end 掉率 / --local 软截回血（SKILL 关键洞察1）
  - ChIP(-q30 -F1804 --no-mixed --no-discordant) / ATAC(--local --dovetail -X2000) 旗标
  - 灵敏度预设 --very-fast | --sensitive | --very-sensitive
  - MAPQ 上限实测（e2e 42 / local 44，永远到不了 BWA 的 60）
  - 多映射 -k 5
所有结果写入 bowtie2_results.json，供笔记与出图引用。
"""
import os
import re
import subprocess
import json

BASE = os.path.dirname(os.path.abspath(__file__))
BT2 = "/Applications/anaconda3/envs/bioaligners/bin/bowtie2"
BUILD = "/Applications/anaconda3/envs/bioaligners/bin/bowtie2-build"
SAM = "/Applications/anaconda3/bin/samtools"

REF = f"{BASE}/reference.fa"
IDX = f"{BASE}/reference_index"          # basename（正确用法）
IDX_FILE = f"{BASE}/reference_index.1.bt2"  # 文件名（错误用法）
R1 = f"{BASE}/reads_1.fq"
R2 = f"{BASE}/reads_2.fq"

ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"  # Illumina TruSeq adapter


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True, **kw)


def overall_rate(stderr):
    m = re.search(r"([\d.]+)% overall alignment rate", stderr)
    return float(m.group(1)) if m else None


def max_mapq(bam):
    out = sh(f"{SAM} view -F 4 {bam}")
    qs = []
    for line in out.stdout.splitlines():
        if line:
            qs.append(int(line.split("\t")[4]))
    return max(qs) if qs else None


def mapped_count(bam):
    out = sh(f"{SAM} view -F 4 {bam}")
    return sum(1 for _ in out.stdout.splitlines() if _)


results = {}

# ---------- 1. 建索引 ----------
print("[1] bowtie2-build index ...")
r = sh(f"{BUILD} --threads 8 {REF} {IDX}")
print("    build stderr tail:", (r.stderr.strip().splitlines()[-1:] or [""])[0])
print("    rc:", r.returncode)
artifacts = [f for f in os.listdir(BASE) if f.startswith("reference_index") and f.endswith(".bt2")]
print("    index files:", len(artifacts))
results["index"] = {"rc": r.returncode, "bt2_files": len(artifacts)}

# ---------- 2. 基础 PE 比对 (end-to-end 默认) ----------
print("[2] basic end-to-end PE ...")
bam_e2e = f"{BASE}/aligned_e2e.bam"
r = sh(f"{BT2} -p 8 -x {IDX} -1 {R1} -2 {R2} 2>{BASE}/e2e.log | {SAM} sort -@4 -o {bam_e2e} -")
rate_e2e = overall_rate(open(f"{BASE}/e2e.log").read())
sh(f"{SAM} index {bam_e2e}")
results["e2e"] = {"overall_rate": rate_e2e, "max_mapq": max_mapq(bam_e2e),
                  "mapped_reads": mapped_count(bam_e2e)}
print(f"    rate={rate_e2e}%  maxMAPQ={results['e2e']['max_mapq']}")

# ---------- 3. 负向测试：-x 传文件名 ----------
print("[3] negative: -x with a .bt2 filename ...")
r = sh(f"{BT2} -x {IDX_FILE} -1 {R1} -2 {R2} -S {BASE}/neg.sam 2>{BASE}/neg.log")
err = open(f"{BASE}/neg.log").read().strip().splitlines()[-1:] or [""]
results["neg_x_filename"] = {"rc": r.returncode, "last_err": err[0]}
print("    rc:", r.returncode, "|", err[0])

# ---------- 4. 适配器污染 reads 生成 (端读接 20bp 接头) ----------
print("[4] simulate adapter contamination on 3' end of all R1 reads ...")
contam1 = f"{BASE}/reads_contam_1.fq"
with open(R1) as fh, open(contam1, "w") as out:
    while True:
        h = fh.readline()
        if not h:
            break
        seq = fh.readline().rstrip("\n")
        plus = fh.readline()
        qual = fh.readline().rstrip("\n")
        seq2 = seq + ADAPTER
        qual2 = qual + "I" * len(ADAPTER)
        out.write(h + seq2 + "\n" + plus + qual2 + "\n")
print("    wrote", contam1)

# ---------- 5. 污染 reads: end-to-end vs --local ----------
print("[5] adapter-contam: end-to-end vs --local ...")
bam_ce = f"{BASE}/contam_e2e.bam"
r1 = sh(f"{BT2} -p 8 -x {IDX} -1 {contam1} -2 {R2} 2>{BASE}/contam_e2e.log | {SAM} sort -@4 -o {bam_ce} -")
rate_ce = overall_rate(open(f"{BASE}/contam_e2e.log").read())
sh(f"{SAM} index {bam_ce}")

bam_cl = f"{BASE}/contam_local.bam"
r2 = sh(f"{BT2} -p 8 --local -x {IDX} -1 {contam1} -2 {R2} 2>{BASE}/contam_local.log | {SAM} sort -@4 -o {bam_cl} -")
rate_cl = overall_rate(open(f"{BASE}/contam_local.log").read())
sh(f"{SAM} index {bam_cl}")

results["adapter_contam"] = {
    "e2e_rate": rate_ce, "local_rate": rate_cl,
    "recovered_by_local": round((rate_cl or 0) - (rate_ce or 0), 2),
    "local_max_mapq": max_mapq(bam_cl),
}
print(f"    e2e={rate_ce}%  local={rate_cl}%  recovered={results['adapter_contam']['recovered_by_local']}pp")

# ---------- 6. ChIP-seq 旗标 ----------
print("[6] ChIP-seq flags ...")
bam_chip = f"{BASE}/chip.bam"
r = sh(f"{BT2} -p 8 --very-sensitive --no-mixed --no-discordant "
       f"--rg-id s1 --rg SM:s1 --rg PL:ILLUMINA --rg LB:lib1 "
       f"-x {IDX} -1 {R1} -2 {R2} 2>{BASE}/chip.log | "
       f"{SAM} view -bS -q 30 -F 1804 - | {SAM} sort -@4 -o {bam_chip} -")
rate_chip = overall_rate(open(f"{BASE}/chip.log").read())
sh(f"{SAM} index {bam_chip}")
results["chip"] = {"overall_rate": rate_chip, "max_mapq": max_mapq(bam_chip),
                   "mapped_after_q30_F1804": mapped_count(bam_chip)}
print(f"    rate={rate_chip}%  maxMAPQ={results['chip']['max_mapq']}")

# ---------- 7. ATAC-seq 旗标 (local + dovetail + X2000) ----------
print("[7] ATAC-seq flags ...")
bam_atac = f"{BASE}/atac.bam"
r = sh(f"{BT2} -p 8 --very-sensitive --local --dovetail -X 2000 --no-mixed --no-discordant "
       f"-x {IDX} -1 {R1} -2 {R2} 2>{BASE}/atac.log | "
       f"{SAM} view -bS -q 30 -F 1804 - | {SAM} sort -@4 -o {bam_atac} -")
rate_atac = overall_rate(open(f"{BASE}/atac.log").read())
sh(f"{SAM} index {bam_atac}")
results["atac"] = {"overall_rate": rate_atac, "max_mapq": max_mapq(bam_atac),
                   "mapped_after_q30_F1804": mapped_count(bam_atac)}
print(f"    rate={rate_atac}%  maxMAPQ={results['atac']['max_mapq']}")

# ---------- 8. 灵敏度预设 ----------
print("[8] sensitivity presets ...")
presets = {}
for p in ["--very-fast", "--sensitive", "--very-sensitive"]:
    bam_p = f"{BASE}/preset_{p.lstrip('-')}.bam"
    r = sh(f"{BT2} {p} -p 8 -x {IDX} -1 {R1} -2 {R2} 2>{BASE}/preset_{p.lstrip('-')}.log | {SAM} sort -@4 -o {bam_p} -")
    presets[p] = overall_rate(open(f"{BASE}/preset_{p.lstrip('-')}.log").read())
results["presets"] = presets
print("   ", presets)

# ---------- 9. 多映射 -k 5 ----------
print("[9] multi-mapping -k 5 ...")
out_sam = f"{BASE}/multimap.sam"
r = sh(f"{BT2} -k 5 -p 8 -x {IDX} -1 {R1} -2 {R2} -S {out_sam} 2>{BASE}/multimap.log")
# 统计被报告多次的 read（secondary 数）
sec = sh(f"{SAM} view -f 256 {out_sam} 2>/dev/null | wc -l").stdout.strip()
results["multimap_k5"] = {"rc": r.returncode, "secondary_alignments": int(sec or 0)}
print("    secondary(二次) alignments:", results["multimap_k5"]["secondary_alignments"])

# ---------- 10. MAPQ 上限：e2e vs local 实测 ----------
results["mapq_cap"] = {
    "e2e_max": results["e2e"]["max_mapq"],
    "local_max": results["adapter_contam"]["local_max_mapq"],
    "bwa_equivalent": 60,
    "note": "Bowtie2 MAPQ 上限 42(e2e)/44(local)，永远到不了 BWA 的 60；-q 60 会把 BAM 清空",
}

with open(f"{BASE}/bowtie2_results.json", "w") as f:
    json.dump(results, f, indent=2, ensure_ascii=False)
print("\n[done] wrote bowtie2_results.json")
print(json.dumps(results, indent=2, ensure_ascii=False))
