#!/usr/bin/env bash
# 040 quality-filtering 真跑脚本（WSL Ubuntu / bio-qc 环境）
# 用法：wsl.exe -d Ubuntu -u root -- bash /mnt/d/.../040-quality-filtering/_run.sh
set -u
BASE="$(cd "$(dirname "$0")" && pwd)"
cd "$BASE"
exec > >(tee -a _run.log) 2>&1

source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio-qc

echo "==== [0] versions ===="
date
fastp --version
cutadapt --version
trimmomatic -version 2>&1 | head -1 || java -jar "$(command -v trimmomatic)" -version 2>&1 | head -1
python3 --version

echo "==== [1] generate quality-gradient FASTQ (seed 20260903) ===="
python3 make_inputs.py

echo "==== [2] fastp per-read filter (-q 20 -u 40 -n 5 -l 36) ===="
fastp -i input_grad.fq.gz -o out_fastp_filter.fq.gz \
      -q 20 -u 40 -n 5 -l 36 \
      -j _fastp_filter.json -h _fastp_filter.html
echo "fastp_filter exit: $?"

echo "==== [3] fastp window trim --cut_right (SKILL 口径: --cut_right --cut_window_size 4 --cut_mean_quality 20 -l 36) ===="
fastp -i input_grad.fq.gz -o out_fastp_cutright.fq.gz \
      --cut_right --cut_window_size 4 --cut_mean_quality 20 -l 36 \
      -j _fastp_cutright.json -h _fastp_cutright.html
echo "fastp_cutright exit: $?"

echo "==== [4] cutadapt -q 20 -m 36 (BWA running-sum 3' trim) ===="
cutadapt -q 20 -m 36 -o out_cutadapt.fq.gz input_grad.fq.gz > _cutadapt_report.txt
echo "cutadapt exit: $?"

echo "==== [5] trimmomatic LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 ===="
trimmomatic SE -phred33 input_grad.fq.gz out_trimm_sw.fq.gz \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 > _trimm_sw.log 2>&1
echo "trimmomatic SW exit: $?"; cat _trimm_sw.log

echo "==== [6] trimmomatic MAXINFO:40:0.5 MINLEN:36 ===="
trimmomatic SE -phred33 input_grad.fq.gz out_trimm_maxinfo.fq.gz \
      MAXINFO:40:0.5 MINLEN:36 > _trimm_maxinfo.log 2>&1
echo "trimmomatic MAXINFO exit: $?"; cat _trimm_maxinfo.log

echo "==== [7] parse metrics ===="
python3 parse_results.py
echo "==== done ===="
date
