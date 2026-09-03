#!/usr/bin/env bash
# 只重跑解析步骤（fastqc/multiqc 产物已落盘且确定）
set -u
source /opt/miniconda3/etc/profile.d/conda.sh
BASE=/mnt/d/1.WorkDir/RedBook/content/素材/read-qc/038-quality-reports
cd "$BASE"
conda activate bio-qc
python3 parse_qc.py
