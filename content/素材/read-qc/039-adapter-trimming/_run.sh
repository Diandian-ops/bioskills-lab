#!/usr/bin/env bash
# 039 adapter-trimming 真跑脚本：在 WSL Ubuntu bio-qc 环境执行。
# 用法: wsl.exe -d Ubuntu -u root -- bash /mnt/d/1.WorkDir/RedBook/content/素材/read-qc/039-adapter-trimming/_run.sh
set -euo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"
: > "$DIR/_run.log"
exec >> "$DIR/_run.log" 2>&1

echo "=== 039 adapter-trimming real run ==="
date
echo

source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio-qc

echo "--- tool versions ---"
cutadapt --version
java -version 2>&1 | head -1
ls -d "$CONDA_PREFIX"/share/trimmomatic-*/
echo "--- built-in adapter files (SKILL.md: ls \$CONDA_PREFIX/share/trimmomatic-*/adapters/) ---"
ls "$CONDA_PREFIX"/share/trimmomatic-*/adapters/
ADAPTERS_FA="$(ls "$CONDA_PREFIX"/share/trimmomatic-*/adapters/TruSeq3-PE-2.fa | head -1)"
echo "--- adapter fasta used: $ADAPTERS_FA ---"
cat "$ADAPTERS_FA"
echo

echo "--- step 1: generate gradient inputs (make_inputs.py, seed 39) ---"
python3 make_inputs.py
echo

for G in 5p 20p 40p; do
  echo "=== cutadapt grad $G (SKILL.md PE 口径: -a/-A AGATCGGAAGAGC -m 20:20) ==="
  # cutadapt 5.2 的运行报告输出到 stdout（实测 stderr 为空），两者一并捕获
  cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 20:20 \
    -o "ca_${G}_R1.fq.gz" -p "ca_${G}_R2.fq.gz" \
    "grad${G}_R1.fq.gz" "grad${G}_R2.fq.gz" \
    > "cutadapt_report_${G}.txt" 2>&1
  tail -32 "cutadapt_report_${G}.txt"
  echo

  echo "=== trimmomatic grad $G · A: SKILL.md 字面口径 (ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads) ==="
  trimmomatic PE -phred33 -threads 4 \
    "grad${G}_R1.fq.gz" "grad${G}_R2.fq.gz" \
    "tm_${G}_R1_p.fq.gz" "tm_${G}_R1_u.fq.gz" \
    "tm_${G}_R2_p.fq.gz" "tm_${G}_R2_u.fq.gz" \
    "ILLUMINACLIP:${ADAPTERS_FA}:2:30:10:2:keepBothReads" MINLEN:36
  echo

  echo "=== trimmomatic grad $G · B: 修正布尔位 (ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:true) ==="
  trimmomatic PE -phred33 -threads 4 \
    "grad${G}_R1.fq.gz" "grad${G}_R2.fq.gz" \
    "tmk_${G}_R1_p.fq.gz" "tmk_${G}_R1_u.fq.gz" \
    "tmk_${G}_R2_p.fq.gz" "tmk_${G}_R2_u.fq.gz" \
    "ILLUMINACLIP:${ADAPTERS_FA}:2:30:10:2:true" MINLEN:36
  echo
done

echo "--- step 3: parse results ---"
python3 analyze_results.py
echo
echo "=== done ==="
date
