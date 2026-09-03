#!/usr/bin/env bash
# 创建独立的 bio-gatk 环境，避免把 bio 环境的 openjdk 25 降级到 17。
set -u
source /opt/miniconda3/etc/profile.d/conda.sh

LOG=/mnt/d/1.WorkDir/RedBook/content/素材/variant-calling/022-gatk-variant-calling/_setup_env.log
echo "=== start $(date -u) ===" | tee -a "$LOG"

conda create -y -n bio-gatk -c bioconda -c conda-forge gatk4 openjdk=17 \
  2>&1 | tee -a "$LOG"

echo "=== create exit: $? ===" | tee -a "$LOG"

# 验证
conda run -n bio-gatk bash -lc 'echo JAVA:; java -version 2>&1 | head -2; echo GATK:; gatk --version 2>&1 | tail -3' \
  2>&1 | tee -a "$LOG"

echo "=== done $(date -u) ===" | tee -a "$LOG"
