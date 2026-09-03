#!/usr/bin/env bash
# 042 重跑解析与汇总（前置步骤产物已就绪，不重复跑比对/dedup）
set -uo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio
echo "coordonly R1 (excl dup) = $(samtools view -c -f 0x40 -F 0xD00 coordonly.bam)"
python3 parse_results.py 2>&1 | tee 09_parse.log
{
  echo "# 042 umi-processing repro transcript"
  echo "# generated: $(date -Is)"
  for f in 00_env_versions.log 01_make_inputs.log 02_extract.log 03_align.log 04_flagstat.log 05_markdup_stats.log 06_dedup_directional.log 07_dedup_unique.log 08_dedup_cluster.log 09_parse.log results.txt; do
    [ -f "$f" ] && { echo ""; echo "===== $f ====="; cat "$f"; }
  done
} > repro_transcript.txt
echo "DONE"
