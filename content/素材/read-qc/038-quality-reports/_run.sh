#!/usr/bin/env bash
# 038 quality-reports 主流程：造数据 -> fastqc -> multiqc -> 解析真实数值
# 全部工具输出重定向到 _run.log；在 WSL Ubuntu 的 bio-qc 环境中执行
set -u
source /opt/miniconda3/etc/profile.d/conda.sh
BASE=/mnt/d/1.WorkDir/RedBook/content/素材/read-qc/038-quality-reports
cd "$BASE"
conda activate bio-qc

{
  echo "== 038 quality-reports run: $(date '+%F_%T') =="
  echo "python: $(python3 --version 2>&1)"
  echo "fastqc: $(fastqc --version 2>&1)"
  echo "multiqc: $(multiqc --version 2>&1)"

  echo
  echo "[1/4] make inputs (make_inputs.py, seed=20260903)"
  python3 make_inputs.py || exit 1
  ls -la raw_fastq/

  echo
  echo "[2/4] fastqc per-file QC"
  mkdir -p qc/raw && rm -rf qc/raw && mkdir -p qc/raw
  fastqc -t 3 -o qc/raw raw_fastq/*.fastq.gz || exit 1

  echo
  echo "[3/4] multiqc aggregation"
  rm -rf qc/multiqc
  multiqc qc/raw -o qc/multiqc -f || exit 1
  ls -la qc/multiqc/multiqc_data/ | head -20

  echo
  echo "[4/4] parse real values (parse_qc.py)"
  python3 parse_qc.py || exit 1

  echo
  echo "== done: $(date '+%F_%T') =="
} > _run.log 2>&1
echo "exit=$?"
tail -60 _run.log
