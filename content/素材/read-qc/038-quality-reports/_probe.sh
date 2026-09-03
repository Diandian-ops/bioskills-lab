#!/usr/bin/env bash
# 038 环境探针：确认 bio-qc / bio 环境可用工具与版本（输出全部重定向日志）
set -u
source /opt/miniconda3/etc/profile.d/conda.sh
LOG=/mnt/d/1.WorkDir/RedBook/content/素材/read-qc/038-quality-reports/_env_probe.log
{
  echo "== probe date: $(date +%F) =="
  conda activate bio-qc
  echo "-- bio-qc env --"
  fastqc --version
  multiqc --version
  fastp --version 2>&1
  samtools --version | head -1
  echo "-- bio env (dwgsim) --"
  conda activate bio
  dwgsim 2>&1 | head -2
} > "$LOG" 2>&1
cat "$LOG"
