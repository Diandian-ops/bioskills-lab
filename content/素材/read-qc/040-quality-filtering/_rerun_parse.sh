#!/usr/bin/env bash
# 040 仅重跑解析步骤（工具输出已存在）
set -u
BASE="$(cd "$(dirname "$0")" && pwd)"
cd "$BASE"
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio-qc
python3 parse_results.py
