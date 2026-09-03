#!/usr/bin/env bash
# 041 fastp-workflow 真实复现主脚本（在 WSL Ubuntu bio-qc 环境内执行）
set -x
cd "$(dirname "$0")"

source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio-qc

echo "=== [1] make simulated PE FASTQ (adapter + low-quality + N pollution) ==="
python3 make_inputs.py

echo "=== [2] fastp FULL run (SKILL.md standard PE workflow + --correction) ==="
fastp -i raw_R1.fq.gz -I raw_R2.fq.gz -o clean_R1.fq.gz -O clean_R2.fq.gz \
      --detect_adapter_for_pe --correction \
      --cut_right --cut_window_size 4 --cut_mean_quality 20 \
      -q 20 -l 36 -w 8 \
      -h report_full.html -j report_full.json

echo "=== [3] ablation A: --disable_adapter_trimming ==="
fastp -i raw_R1.fq.gz -I raw_R2.fq.gz -o clean_noadapter_R1.fq.gz -O clean_noadapter_R2.fq.gz \
      --detect_adapter_for_pe --correction \
      --cut_right --cut_window_size 4 --cut_mean_quality 20 \
      -q 20 -l 36 -w 8 \
      --disable_adapter_trimming \
      -h report_noadapter.html -j report_noadapter.json

echo "=== [4] ablation B: --disable_length_filtering ==="
fastp -i raw_R1.fq.gz -I raw_R2.fq.gz -o clean_nolen_R1.fq.gz -O clean_nolen_R2.fq.gz \
      --detect_adapter_for_pe --correction \
      --cut_right --cut_window_size 4 --cut_mean_quality 20 \
      -q 20 -l 36 -w 8 \
      --disable_length_filtering \
      -h report_nolenfilter.html -j report_nolenfilter.json

echo "=== [5] ablation C: drop --cut_right (per-read quality/N filters fire) ==="
fastp -i raw_R1.fq.gz -I raw_R2.fq.gz -o clean_nocut_R1.fq.gz -O clean_nocut_R2.fq.gz \
      --detect_adapter_for_pe --correction \
      -q 20 -l 36 -w 8 \
      -h report_nocut.html -j report_nocut.json

echo "=== [6] parse fastp JSON reports -> _metrics.tsv ==="
python3 parse_reports.py

echo "=== [7] output files ==="
ls -la
echo "ALL DONE"
