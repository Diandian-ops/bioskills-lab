#!/usr/bin/env python3
"""012 bwa-alignment — 忠实复现 bioSkills SKILL.md（bwa-mem2 主工具）。

覆盖点：
  - bwa-mem2 index 建索引（产物 .0123/.amb/.ann/.bwt.2bit.64/.pac，与 bwa index 不互通）
  - read group 硬契约：带 -R 与不带 -R 的 @RG 头差异（GATK 契约）
  - 去重严格顺序 collate -> fixmate -m -> sort -> markdup（验证 MC/MS 标签落地）
  - -Y（软截 supplementary 保全程）vs -M（降级为 secondary，SV 不可用）
  - -K 100000000 跨线程可复现（两次跑指纹一致）
  - 索引错配负向测试：bwa-mem2 用 bwa 原版索引 -> 报错
结果写入 bwa_results.json。
"""
import os, re, subprocess, json, hashlib

BASE = os.path.dirname(os.path.abspath(__file__))
MEM2 = "/Applications/anaconda3/envs/bioaligners/bin/bwa-mem2"
BWA = "/Applications/anaconda3/envs/bioaligners/bin/bwa"
SAM = "/Applications/anaconda3/bin/samtools"
REF = f"{BASE}/reference.fa"
PREFIX = f"{BASE}/reference.fa"   # bwa-mem2 index 默认以参考名为索引前缀
R1 = f"{BASE}/reads_1.fq"
R2 = f"{BASE}/reads_2.fq"
RG = r"'@RG\tID:s1\tSM:s1\tPL:ILLUMINA\tLB:lib1'"


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True, **kw)


def head_has_rg(bam):
    return "@RG" in sh(f"{SAM} view -H {bam}").stdout


def md5_aln(bam):
    out = sh(f"{SAM} view -F 4 {bam} | cut -f1-9").stdout
    return hashlib.md5(out.encode()).hexdigest()


def flag_count(bam, flagbits):
    n = 0
    for line in sh(f"{SAM} view {bam}").stdout.splitlines():
        if line:
            f = int(line.split("\t")[1])
            if f & flagbits:
                n += 1
    return n


def max_mapq(bam):
    qs = [int(l.split("\t")[4]) for l in sh(f"{SAM} view -F 4 {bam}").stdout.splitlines() if l]
    return max(qs) if qs else None


def idx_files():
    exts = (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
    return [f for f in os.listdir(BASE) if f.startswith("reference.fa") and f.endswith(exts)]


results = {}

# 1. 建索引
print("[1] bwa-mem2 index ...")
r = sh(f"{MEM2} index {REF} 2>{BASE}/idx.log")
results["index"] = {"rc": r.returncode, "files": len(idx_files())}
print("   rc", r.returncode, "index files:", len(idx_files()))

# 2. 带 read group 比对
print("[2] align WITH read group -> sorted BAM ...")
bam_rg = f"{BASE}/aligned_rg.bam"
r = sh(f"{MEM2} mem -t 8 -R {RG} {PREFIX} {R1} {R2} 2>{BASE}/align_rg.log | {SAM} sort -@4 -o {bam_rg} -")
sh(f"{SAM} index {bam_rg}")
results["with_rg"] = {"rc": r.returncode, "has_RG_header": head_has_rg(bam_rg),
                      "max_mapq": max_mapq(bam_rg),
                      "mapped_reads": int(sh(f"{SAM} view -F 4 -c {bam_rg}").stdout.strip())}
print("   has @RG:", head_has_rg(bam_rg), "maxMAPQ:", results["with_rg"]["max_mapq"])

# 3. 不带 read group 比对（契约负向）
print("[3] align WITHOUT read group (contract negative) ...")
bam_norg = f"{BASE}/aligned_norg.bam"
r = sh(f"{MEM2} mem -t 8 {PREFIX} {R1} {R2} 2>{BASE}/align_norg.log | {SAM} sort -@4 -o {bam_norg} -")
sh(f"{SAM} index {bam_norg}")
results["without_rg"] = {"has_RG_header": head_has_rg(bam_norg)}
print("   has @RG (should be False):", head_has_rg(bam_norg))

# 4. 去重严格顺序
print("[4] dedup strict ordering: collate -> fixmate -m -> sort -> markdup ...")
bam_dup = f"{BASE}/aligned.markdup.bam"
r = sh(f"{MEM2} mem -t 8 -R {RG} {PREFIX} {R1} {R2} 2>{BASE}/md.log | "
       f"{SAM} collate -@4 -O -u - | {SAM} fixmate -m -@4 -u - - | "
       f"{SAM} sort -@4 -u - | {SAM} markdup -@4 - {bam_dup}")
sh(f"{SAM} index {bam_dup}")
mc = sh(f"{SAM} view {bam_dup} | grep -m1 'MC:Z' | wc -l").stdout.strip()
results["markdup"] = {"rc": r.returncode, "has_MC_tag": int(mc or 0) > 0,
                      "dup_flagged": flag_count(bam_dup, 1024)}
print("   has MC tag:", results["markdup"]["has_MC_tag"], "dup flagged:", results["markdup"]["dup_flagged"])

# 5. -Y vs -M（SV 信号语义）
print("[5] -Y (supplementary kept) vs -M (secondary) ...")
bam_y = f"{BASE}/aligned_Y.bam"
sh(f"{MEM2} mem -t 8 -Y -R {RG} {PREFIX} {R1} {R2} 2>{BASE}/y.log | {SAM} sort -@4 -o {bam_y} -")
bam_m = f"{BASE}/aligned_M.bam"
sh(f"{MEM2} mem -t 8 -M -R {RG} {PREFIX} {R1} {R2} 2>{BASE}/m.log | {SAM} sort -@4 -o {bam_m} -")
results["Y_vs_M"] = {
    "Y_supplementary_0x800": flag_count(bam_y, 2048),
    "M_secondary_0x100": flag_count(bam_m, 256),
    "note": "合成参考无 SV 断点，二者均为 0；差异在语义：SV 时 -M 把分裂片段降为 secondary(0x100) 会被 SV caller 跳过，-Y 保留 supplementary(0x800) 全程序列",
}
print("   -Y supp:", results["Y_vs_M"]["Y_supplementary_0x800"],
      "| -M sec:", results["Y_vs_M"]["M_secondary_0x100"])

# 6. -K 可复现
print("[6] -K 100000000 reproducibility ...")
b1, b2 = f"{BASE}/repro1.bam", f"{BASE}/repro2.bam"
sh(f"{MEM2} mem -t 8 -K 100000000 -R {RG} {PREFIX} {R1} {R2} 2>{BASE}/repro1.log | {SAM} sort -@4 -o {b1} -")
sh(f"{MEM2} mem -t 8 -K 100000000 -R {RG} {PREFIX} {R1} {R2} 2>{BASE}/repro2.log | {SAM} sort -@4 -o {b2} -")
m1, m2 = md5_aln(b1), md5_aln(b2)
results["reproducible_K"] = {"md5_1": m1, "md5_2": m2, "identical": m1 == m2}
print("   identical:", m1 == m2)

# 7. 索引错配负向测试：bwa 原版索引给 bwa-mem2 用
print("[7] negative: bwa-mem2 on a bwa-original index ...")
sh(f"mkdir -p {BASE}/bwaonly && {BWA} index -p {BASE}/bwaonly/ref {REF} 2>{BASE}/bwaidx.log")
r = sh(f"{MEM2} mem -t 8 {BASE}/bwaonly/ref {R1} {R2} 2>{BASE}/mismatch.log")
err = open(f"{BASE}/mismatch.log").read().strip().splitlines()[-1:] or [""]
results["idx_mismatch"] = {"rc": r.returncode, "last_err": err[0]}
print("   rc:", r.returncode, "|", err[0])

with open(f"{BASE}/bwa_results.json", "w") as f:
    json.dump(results, f, indent=2, ensure_ascii=False)
print("\n[done] wrote bwa_results.json")
print(json.dumps(results, indent=2, ensure_ascii=False))
