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
import os, re, shutil, subprocess, json, hashlib

BASE = os.path.dirname(os.path.abspath(__file__))


def tool(env_var, name, legacy=()):
    """定位外部二进制：环境变量 > PATH > 旧版绝对路径（历史 mac 安装）。

    跨平台：不在代码里写死 /Applications/anaconda3/... 这类 macOS 路径，
    Windows 原生（conda）同样能靠 PATH 或环境变量定位。
    """
    p = os.environ.get(env_var)
    if p and os.path.exists(p):
        return p
    p = shutil.which(name)
    if p:
        return p
    for lp in legacy:
        if os.path.exists(lp):
            return lp
    raise SystemExit(
        "[missing-tool] 未找到 %s。请安装后用环境变量指定路径：\n"
        "    set %s=C:\\path\\to\\%s          (Windows)\n"
        "    export %s=/path/to/%s           (macOS/Linux)"
        % (name, env_var, name, env_var, name))


MEM2 = tool("BWA_MEM2", "bwa-mem2",
            ("/Applications/anaconda3/envs/bioaligners/bin/bwa-mem2",))
BWA = tool("BWA", "bwa",
           ("/Applications/anaconda3/envs/bioaligners/bin/bwa",))
SAM = tool("SAMTOOLS", "samtools",
           ("/Applications/anaconda3/bin/samtools",))
REF = f"{BASE}/reference.fa"
PREFIX = f"{BASE}/reference.fa"   # bwa-mem2 index 默认以参考名为索引前缀
R1 = f"{BASE}/reads_1.fq"
R2 = f"{BASE}/reads_2.fq"
# 注意：不再用外层单引号包裹。POSIX shell 靠单引号保留 \t，
# 但 Windows cmd.exe 不识别单引号，会把引号一起传给 bwa-mem2。
# 该串本身不含空格，两种 shell 下裸写都安全。
RG = r"@RG\tID:s1\tSM:s1\tPL:ILLUMINA\tLB:lib1"


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True, **kw)


def head_has_rg(bam):
    return "@RG" in sh(f"{SAM} view -H {bam}").stdout


def md5_aln(bam):
    # 原为 `samtools view -F 4 <bam> | cut -f1-9`。
    # cut 属 Unix coreutils，Windows 原生环境没有，故改为 Python 内截取（语义等价）：
    # 取每行前 9 个字段；不足 9 段时 cut 也是原样输出，join 行为一致。
    raw = sh(f"{SAM} view -F 4 {bam}").stdout
    lines = ["\t".join(l.split("\t")[:9]) for l in raw.splitlines() if l]
    out = "".join(l + "\n" for l in lines)
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
# 原为 `samtools view <bam> | grep -m1 'MC:Z' | wc -l`：
# grep/wc 同样不是 Windows 原生命令。语义是「是否存在至少一行含 MC:Z」。
raw_dup = sh(f"{SAM} view {bam_dup}").stdout
has_mc = any("MC:Z" in l for l in raw_dup.splitlines())
results["markdup"] = {"rc": r.returncode, "has_MC_tag": has_mc,
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
