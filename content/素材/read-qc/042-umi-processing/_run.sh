#!/usr/bin/env bash
# 042 umi-processing 真跑主流程：extract -> bwa 比对 -> dedup x3 + markdup 对照 -> 解析
set -uo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"
source /opt/miniconda3/etc/profile.d/conda.sh

echo "== [0] 环境版本 =="
conda list -n bio-umi 2>/dev/null | grep -E "^(umi_tools|python) " | tee 00_env_versions.log
conda list -n bio 2>/dev/null | grep -E "^(bwa|samtools) " | tee -a 00_env_versions.log

echo "== [1] 造数据 =="
python3 make_inputs.py 2>&1 | tee 01_make_inputs.log

echo "== [2] umi_tools extract（bio-umi）=="
conda activate bio-umi
umi_tools extract --stdin=R1.fq.gz --read2-in=R2.fq.gz \
    --stdout=R1_umi.fq.gz --read2-out=R2_umi.fq.gz \
    --bc-pattern=NNNNNNNNNNNN --log=extract.log 2>&1 | tee 02_extract.log
echo "---- extract.log ----"
cat extract.log

echo "== [3] bwa 比对（bio）=="
conda deactivate
conda activate bio
bwa index ref.fa 2>&1 | tail -2 | tee 03_align.log
bwa mem -t 4 -R "@RG\tID:rg1\tSM:sample1\tLB:lib1\tPL:ILLUMINA" \
    ref.fa R1_umi.fq.gz R2_umi.fq.gz 2>>03_align.log | samtools sort -o sorted.bam
samtools index sorted.bam
samtools flagstat sorted.bam | tee 04_flagstat.log

echo "== [4] 坐标法对照：samtools markdup（bio）=="
samtools sort -n -o namesort.bam sorted.bam
samtools fixmate -m namesort.bam fixmate.bam
samtools sort -o possort.bam fixmate.bam
samtools markdup -s possort.bam coordonly.bam 2>05_markdup_stats.log
echo "markdup exit=$?"
cat 05_markdup_stats.log
echo "coordonly.bam primary R1 = $(samtools view -c -f 0x40 -F 0x900 coordonly.bam)"

echo "== [5] umi_tools dedup x3（bio-umi）=="
conda deactivate
conda activate bio-umi
umi_tools dedup -I sorted.bam -S dedup_directional.bam --paired \
    --output-stats=stats_dir --log=dedup_directional.log 2>&1 | tee 06_dedup_directional.log
umi_tools dedup -I sorted.bam -S dedup_unique.bam --paired --method=unique \
    --log=dedup_unique.log 2>&1 | tee 07_dedup_unique.log
umi_tools dedup -I sorted.bam -S dedup_cluster.bam --paired --method=cluster \
    --log=dedup_cluster.log 2>&1 | tee 08_dedup_cluster.log
for m in directional unique cluster; do
  echo "dedup_$m primary R1 = $(samtools view -c -f 0x40 -F 0x900 dedup_$m.bam 2>/dev/null || echo NA)"
done
echo "---- dedup_directional.log ----"
cat dedup_directional.log

echo "== [6] 解析 =="
conda deactivate
conda activate bio
python3 parse_results.py 2>&1 | tee 09_parse.log

echo "== [7] 汇总 repro_transcript.txt =="
{
  echo "# 042 umi-processing repro transcript"
  echo "# generated: $(date -Is)"
  for f in 00_env_versions.log 01_make_inputs.log 02_extract.log extract.log \
           03_align.log 04_flagstat.log 05_markdup_stats.log \
           06_dedup_directional.log dedup_directional.log \
           07_dedup_unique.log 08_dedup_cluster.log 09_parse.log results.txt; do
    [ -f "$f" ] && { echo ""; echo "===== $f ====="; cat "$f"; }
  done
} > repro_transcript.txt
echo "DONE"
