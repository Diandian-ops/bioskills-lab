#!/usr/bin/env python3
"""013 hisat2-alignment — 忠实复现 bioSkills SKILL.md（HISAT2 剪接比对）。

覆盖点：
  - hisat2-build 建索引（.ht2 文件，-x 期望 basename）
  - 基础 PE 比对 + --rna-strandness RF（XS 链标签）
  - 剪接能力：手工跨 intron 的 reads -> HISAT2 产出含 N 的剪接 CIGAR（核心演示）
  - --dta（转录本组装模式，文档转述 + 跑通）
  - 两趟法：--novel-splicesite-outfile 发现 junction，--novel-splicesite-infile 复用
  - -x 传 .ht2 文件名负向测试
结果写入 hisat2_results.json。
"""
import os, re, shutil, subprocess, json

BASE = os.path.dirname(os.path.abspath(__file__))


def tool(env_var, name, legacy=()):
    """定位外部二进制：环境变量 > PATH > 旧版绝对路径（历史 mac 安装）。

    跨平台：不写死 /Applications/anaconda3/... 这类 macOS 路径，
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


H2 = tool("HISAT2", "hisat2",
          ("/Applications/anaconda3/envs/bioaligners/bin/hisat2",))
BUILD = tool("HISAT2_BUILD", "hisat2-build",
             ("/Applications/anaconda3/envs/bioaligners/bin/hisat2-build",))
SAM = tool("SAMTOOLS", "samtools",
           ("/Applications/anaconda3/bin/samtools",))
REF = f"{BASE}/reference.fa"
IDX = f"{BASE}/hisat2_index"   # basename
R1 = f"{BASE}/reads_1.fq"
R2 = f"{BASE}/reads_2.fq"
JR = f"{BASE}/junction_read.fa"


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True, **kw)


def cigar_of(bam, readname=None):
    out = sh(f"{SAM} view {bam}").stdout
    for line in out.splitlines():
        if not line:
            continue
        if readname and not line.startswith(readname):
            continue
        return line.split("\t")[5]
    return None


def has_xs(bam):
    return "XS:A:" in sh(f"{SAM} view {bam}").stdout


def overall_rate(log):
    m = re.search(r"([\d.]+)% overall alignment rate", open(log).read())
    return float(m.group(1)) if m else None


def max_mapq(bam):
    qs = [int(l.split("\t")[4]) for l in sh(f"{SAM} view -F 4 {bam}").stdout.splitlines() if l]
    return max(qs) if qs else None


results = {}

# 1. 建索引
print("[1] hisat2-build index ...")
r = sh(f"{BUILD} -p 8 {REF} {IDX} 2>{BASE}/idx.log")
arts = [f for f in os.listdir(BASE) if f.startswith("hisat2_index") and f.endswith(".ht2")]
results["index"] = {"rc": r.returncode, "ht2_files": len(arts)}
print("   rc", r.returncode, "ht2:", len(arts))

# 2. 基础 PE 比对
print("[2] basic PE alignment ...")
bam_pe = f"{BASE}/aligned_pe.bam"
r = sh(f"{H2} -p 8 -x {IDX} -1 {R1} -2 {R2} 2>{BASE}/pe.log | {SAM} sort -@4 -o {bam_pe} -")
sh(f"{SAM} index {bam_pe}")
results["pe"] = {"overall_rate": overall_rate(f"{BASE}/pe.log"), "max_mapq": max_mapq(bam_pe)}
print("   rate", results["pe"]["overall_rate"], "maxMAPQ", results["pe"]["max_mapq"])

# 3. 剪接能力：跨 intron 的 reads -> 含 N 的 CIGAR
print("[3] spliced read -> N in CIGAR ...")
bam_j = f"{BASE}/junction.bam"
r = sh(f"{H2} -p 8 -x {IDX} -U {JR} -f 2>{BASE}/j.log | {SAM} sort -@4 -o {bam_j} -")
sh(f"{SAM} index {bam_j}")
cig = cigar_of(bam_j, "junction1")
spliced = bool(cig) and "N" in cig
results["splice"] = {"cigar": cig, "is_spliced": spliced}
print("   cigar:", cig, "| spliced:", spliced)

# 4. --rna-strandness RF
print("[4] --rna-strandness RF (XS tag) ...")
bam_rf = f"{BASE}/aligned_rf.bam"
r = sh(f"{H2} -p 8 --rna-strandness RF -x {IDX} -1 {R1} -2 {R2} 2>{BASE}/rf.log | {SAM} sort -@4 -o {bam_rf} -")
sh(f"{SAM} index {bam_rf}")
results["strand_rf"] = {"has_XS_tag": has_xs(bam_rf)}
print("   has XS:", has_xs(bam_rf))

# 5. --dta（转录本组装模式）
print("[5] --dta transcript-assembly mode ...")
bam_dta = f"{BASE}/aligned_dta.bam"
r = sh(f"{H2} -p 8 --dta -x {IDX} -1 {R1} -2 {R2} 2>{BASE}/dta.log | {SAM} sort -@4 -o {bam_dta} -")
results["dta"] = {"rc": r.returncode, "overall_rate": overall_rate(f"{BASE}/dta.log")}
print("   rc", r.returncode, "rate", results["dta"]["overall_rate"])

# 6. 两趟法：发现并复用 novel junction
print("[6] two-pass: discover then reuse novel splice site ...")
nov = f"{BASE}/junction.novel.txt"
r1 = sh(f"{H2} -p 8 -x {IDX} -U {JR} -f --novel-splicesite-outfile {nov} 2>{BASE}/nov.log -S /dev/null")
nov_lines = open(nov).read().strip().splitlines() if os.path.exists(nov) else []
bam_j2 = f"{BASE}/junction_pass2.bam"
r2 = sh(f"{H2} -p 8 -x {IDX} -U {JR} -f --novel-splicesite-infile {nov} 2>{BASE}/nov2.log | {SAM} sort -@4 -o {bam_j2} -")
sh(f"{SAM} index {bam_j2}")
cig2 = cigar_of(bam_j2, "junction1")
results["two_pass"] = {"discovered_junctions": len(nov_lines),
                       "pass2_cigar": cig2,
                       "pass2_spliced": bool(cig2) and "N" in cig2}
print("   discovered:", len(nov_lines), "| pass2 cigar:", cig2)

# 7. -x 传 .ht2 文件名负向测试
print("[7] negative: -x with a .ht2 filename ...")
ht2file = sorted([f for f in os.listdir(BASE) if f.startswith('hisat2_index') and f.endswith('.ht2')])[0]
r = sh(f"{H2} -p 8 -x {BASE}/{ht2file} -1 {R1} -2 {R2} 2>{BASE}/neg.log -S {BASE}/neg.sam")
err = open(f"{BASE}/neg.log").read().strip().splitlines()[-1:] or [""]
results["neg_x_file"] = {"rc": r.returncode, "last_err": err[0]}
print("   rc", r.returncode, "|", err[0])

with open(f"{BASE}/hisat2_results.json", "w") as f:
    json.dump(results, f, indent=2, ensure_ascii=False)
print("\n[done] wrote hisat2_results.json")
print(json.dumps(results, indent=2, ensure_ascii=False))
