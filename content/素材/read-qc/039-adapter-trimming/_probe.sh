#!/usr/bin/env bash
# 探针2：验证 ILLUMINACLIP 第6段布尔位的真实语义
set -uo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"
exec >> "$DIR/_probe.log" 2>&1

source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio-qc
ADAPTERS_FA="$(ls "$CONDA_PREFIX"/share/trimmomatic-*/adapters/TruSeq3-PE-2.fa | head -1)"

echo
echo "=== probe B: 5-param form ILLUMINACLIP:fa:2:30:10:false ==="
trimmomatic PE -phred33 -threads 4 \
  grad5p_R1.fq.gz grad5p_R2.fq.gz \
  tm4b_5p_R1_p.fq.gz tm4b_5p_R1_u.fq.gz tm4b_5p_R2_p.fq.gz tm4b_5p_R2_u.fq.gz \
  "ILLUMINACLIP:${ADAPTERS_FA}:2:30:10:false" MINLEN:36 || echo "probe B exit=$?"

echo
echo "=== probe D: 6-param form ILLUMINACLIP:fa:2:30:10:2:true ==="
trimmomatic PE -phred33 -threads 4 \
  grad5p_R1.fq.gz grad5p_R2.fq.gz \
  tm4d_5p_R1_p.fq.gz tm4d_5p_R1_u.fq.gz tm4d_5p_R2_p.fq.gz tm4d_5p_R2_u.fq.gz \
  "ILLUMINACLIP:${ADAPTERS_FA}:2:30:10:2:true" MINLEN:36 || echo "probe D exit=$?"
echo "=== probe2 done ==="
