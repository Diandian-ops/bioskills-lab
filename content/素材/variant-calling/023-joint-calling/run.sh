#!/bin/bash
# 023 joint-calling 真实跑脚本（WSL / bcftools）
# 数据源：同目录 chr22_slice.vcf.gz（1000G Phase3 chr22 切片，2504 样本 5431 位点）
# 链路：拆逐样本 callset（-e 'GT="ref"' 模拟单样本调用只报变异位点）
#       -> bcftools merge 朴素合并（反模式）
#       -> 与真·联合调用（6 样本子集，平方化矩阵）对比，量化回填
set -e
HERE=/mnt/d/1.WorkDir/RedBook/content/素材/variant-calling/023-joint-calling
IN=$HERE/chr22_slice.vcf.gz
LOG=$HERE/run_log.txt
JSON=$HERE/joint_calling_stats.json
: > "$LOG"

echo "=== 0. 输入概览 ===" | tee -a "$LOG"
TOT=$(bcftools view -H "$IN" | wc -l)
echo "input_records=$TOT" | tee -a "$LOG"
NSAMP=$(bcftools query -l "$IN" | wc -l)
echo "input_samples=$NSAMP" | tee -a "$LOG"

# 取 6 个分布于族群名单的样本（避免首 6 个同族群、变异位点高度重叠导致退化）
bcftools query -l "$IN" > "$HERE/all_samples.txt"
SAMPLES=$(sed -n '1p;420p;840p;1260p;1680p;2100p' "$HERE/all_samples.txt" | paste -sd, -)
echo "SAMPLES=$SAMPLES" | tee -a "$LOG"

# 1) 逐样本 callset：仅保留该样本为非参考（携带变异）的位点，模拟单样本调用产物
for s in $(echo "$SAMPLES" | tr ',' ' '); do
  # 两步法：先 -s 选出样本，再用 -e 'GT="ref"' 排除该样本为纯合参考的位点
  # （一步法 -s 与 -e 同时写会让 GT 表达式在全部样本上求值，得到错误结果）
  bcftools view -s "$s" "$IN" -Ou | bcftools view -e 'GT="ref"' -Oz -o "$HERE/per_${s}.vcf.gz"
  bcftools index -t "$HERE/per_${s}.vcf.gz"
  n=$(bcftools view -H "$HERE/per_${s}.vcf.gz" | wc -l)
  echo "per_sample_variant_sites $s = $n" | tee -a "$LOG"
done

# 2) 朴素合并（bcftools merge 单样本 callset 的反模式）
bcftools merge -Oz -o "$HERE/merged_naive.vcf.gz" "$HERE"/per_*.vcf.gz
bcftools index -t "$HERE/merged_naive.vcf.gz"
M=$(bcftools view -H "$HERE/merged_naive.vcf.gz" | wc -l)
echo "merged_naive_records=$M" | tee -a "$LOG"

# 3) 真·联合调用（6 样本子集，全部位点全样本定型，平方化矩阵）
bcftools view -s "$SAMPLES" "$IN" -Oz -o "$HERE/joint6.vcf.gz"
bcftools index -t "$HERE/joint6.vcf.gz"
L=$(bcftools view -H "$HERE/joint6.vcf.gz" | wc -l)
echo "joint6_records=$L" | tee -a "$LOG"

# 4) 逐样本 ./. 统计（直接数 GT 字段，不依赖 PSC 列语义；awk 保证 0 计数也退出 0）
for s in $(echo "$SAMPLES" | tr ',' ' '); do
  # 朴素合并里该样本的 ./. 数 = 联合调用回填到 0/0 的位点数
  miss=$(bcftools view -s "$s" "$HERE/merged_naive.vcf.gz" -Ou | bcftools query -f '[%GT]\n' | awk '$0=="./."{c++} END{print c+0}')
  echo "naive_missing $s = $miss" | tee -a "$LOG"
  # 真·联合里该样本缺失数（应为 0，验证平方化）
  jmiss=$(bcftools view -s "$s" "$HERE/joint6.vcf.gz" -Ou | bcftools query -f '[%GT]\n' | awk '$0=="./."{c++} END{print c+0}')
  echo "joint_missing $s = $jmiss" | tee -a "$LOG"
done

# 5) 合并前后 GT 一致性：仅在该样本被调用的位点上逐位点比对
#    raw = 字符串直接比；norm = 等位无序比（消除 1/2 与 2/1 相位/排序差，反映真实生物学一致性）
echo "=== GT consistency (called sites only) ===" | tee -a "$LOG"
normGT() { awk -F'\t' '{n=split($2,a,"[|/]"); if(a[1]>a[2]){t=a[1];a[1]=a[2];a[2]=t}; print $1"\t"a[1]"/"a[2]}'; }
for s in $(echo "$SAMPLES" | tr ',' ' '); do
  bcftools view -s "$s" "$HERE/merged_naive.vcf.gz" -Ou | bcftools query -f '%POS\t[%GT]' > /tmp/m_$s.txt
  bcftools view -s "$s" "$HERE/joint6.vcf.gz" -Ou | bcftools query -f '%POS\t[%GT]' > /tmp/j_$s.txt
  called=$(wc -l < /tmp/m_$s.txt)
  raw=$(join -t $'\t' -j 1 <(sort /tmp/m_$s.txt) <(sort /tmp/j_$s.txt) | awk -F'\t' '($2!="./.") && ($2!=$3){c++} END{print c+0}')
  normc=$(join -t $'\t' -j 1 <(sort /tmp/m_$s.txt | normGT) <(sort /tmp/j_$s.txt | normGT) | awk -F'\t' '($2!="./.") && ($2!=$3){c++} END{print c+0}')
  echo "gtcheck $s called=$called raw_mismatch=$raw norm_mismatch=$normc" | tee -a "$LOG"
done

# 6) 汇总 JSON（供 make_figs.py，自包含、无需再调用 bcftools）
python3 - "$HERE" "$LOG" "$JSON" <<'PY'
import sys, json, re
here, logp, jsonp = sys.argv[1], sys.argv[2], sys.argv[3]
lines = open(logp, encoding="utf-8").read().splitlines()

def grab(pat):
    for ln in lines:
        m = re.search(pat, ln)
        if m: return m.group(1)
    return None

samples = grab(r"^SAMPLES=(.+)").split(",")
M = int(grab(r"^merged_naive_records=(\d+)"))
L = int(grab(r"^joint6_records=(\d+)"))
per = {}
for ln in lines:
    m = re.match(r"^per_sample_variant_sites (\S+) = (\d+)$", ln)
    if m: per.setdefault(m.group(1), {})["variant_sites"] = int(m.group(2))
    m = re.match(r"^naive_missing (\S+) = (\d+)$", ln)
    if m: per.setdefault(m.group(1), {})["naive_missing"] = int(m.group(2))
    m = re.match(r"^joint_missing (\S+) = (\d+)$", ln)
    if m: per.setdefault(m.group(1), {})["joint_missing"] = int(m.group(2))
    m = re.match(r"^gtcheck (\S+) called=(\d+) raw_mismatch=(\d+) norm_mismatch=(\d+)$", ln)
    if m:
        per.setdefault(m.group(1), {})["called"] = int(m.group(2))
        per.setdefault(m.group(1), {})["gt_raw_mismatch"] = int(m.group(3))
        per.setdefault(m.group(1), {})["gt_norm_mismatch"] = int(m.group(4))
# 回填数 = 朴素合并里该样本的 ./. 数（真·联合调用把这些格子填成 0/0 等真实基因型）
for s in per:
    per[s]["backfilled"] = per[s].get("naive_missing", 0)
    # GT 一致性（等位无序）= 1 - norm_mismatch / called
    c = per[s].get("called", 0)
    nm = per[s].get("gt_norm_mismatch", 0)
    per[s]["gt_consistency"] = round(100.0 * (c - nm) / c, 4) if c else 100.0

data = {
    "source": "chr22_slice.vcf.gz (1000G Phase3 chr22:17.0-17.2Mb, GRCh37/hs37d5)",
    "samples": samples,
    "input_records": int(grab(r"^input_records=(\d+)")),
    "input_samples": int(grab(r"^input_samples=(\d+)")),
    "merged_records": M,
    "joint_records": L,
    "per_sample": per,
    "total_missing_naive": sum(v.get("naive_missing",0) for v in per.values()),
    "total_genotypes_joint": L * len(samples),
    "total_genotypes_naive": M * len(samples),
}
json.dump(data, open(jsonp,"w"), indent=2)
print("WROTE", jsonp)
print(json.dumps(data, indent=2))
PY

echo "DONE"
