#!/usr/bin/env bash
# 1) 安装 bcftools（支撑 023 joint-calling）  2) 探测 gatk4 可安装性（022）
set -u
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio

echo "===== A. 安装 bcftools ====="
timeout 900 conda install -y -c bioconda bcftools 2>&1 | tail -12
echo "--- 安装后校验 ---"
command -v bcftools >/dev/null 2>&1 && bcftools --version 2>&1 | head -1 || echo "bcftools STILL MISSING"

echo ""
echo "===== B. gatk4 可安装性探测 (dry-run, 不安装) ====="
timeout 900 conda install -y --dry-run -c bioconda gatk4 2>&1 | tail -20
echo "GATK_DRYRUN_EXIT=$?"
