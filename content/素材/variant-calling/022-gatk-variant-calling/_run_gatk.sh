#!/usr/bin/env bash
# 真跑 GATK HaplotypeCaller（标准模式 + -ERC GVCF 模式），并收集可量化统计。
# 全部输出落盘到素材目录，命令本身不向终端喷大段文本。
set -u
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio-gatk

MAT=/mnt/d/1.WorkDir/RedBook/content/素材/variant-calling/022-gatk-variant-calling
REF=$MAT/reference.fa
BAM=$MAT/aligned_e2e.bam
BAMRG=$MAT/aligned_rg.bam
LOG=$MAT/_run_gatk.log
: > "$LOG"

echo "=== gatk version ===" | tee -a "$LOG"
gatk --version 2>&1 | tail -3 | tee -a "$LOG"
echo "java:" | tee -a "$LOG"
java -version 2>&1 | head -2 | tee -a "$LOG"

echo "=== CreateSequenceDictionary ===" | tee -a "$LOG"
gatk CreateSequenceDictionary -R "$REF" -O "$MAT/reference.dict" 2>&1 | tail -5 | tee -a "$LOG"
ls -la "$MAT/reference.dict" | tee -a "$LOG"

echo "=== add read group (GATK 需要 @RG) ===" | tee -a "$LOG"
samtools addreplacerg -r "@RG\tID:rg1\tSM:sample1\tPL:ILLUMINA\tLB:lib1" -o "$BAMRG" "$BAM" 2>&1 | tee -a "$LOG"
samtools index "$BAMRG" 2>&1 | tee -a "$LOG"
samtools view -H "$BAMRG" | grep -m1 "@RG" | tee -a "$LOG"

echo "=== HaplotypeCaller 标准模式 ===" | tee -a "$LOG"
gatk HaplotypeCaller -R "$REF" -I "$BAMRG" -O "$MAT/raw.vcf.gz" \
  --native-pair-hmm-threads 4 2>&1 | tail -8 | tee -a "$LOG"
ls -la "$MAT/raw.vcf.gz" | tee -a "$LOG"
bcftools index -t "$MAT/raw.vcf.gz" 2>&1 | tee -a "$LOG"

echo "=== HaplotypeCaller -ERC GVCF 模式 ===" | tee -a "$LOG"
gatk HaplotypeCaller -R "$REF" -I "$BAMRG" -ERC GVCF -O "$MAT/raw.g.vcf.gz" \
  --native-pair-hmm-threads 4 2>&1 | tail -8 | tee -a "$LOG"
ls -la "$MAT/raw.g.vcf.gz" | tee -a "$LOG"
bcftools index -t "$MAT/raw.g.vcf.gz" 2>&1 | tee -a "$LOG"

echo "=== bcftools stats: 标准模式 ===" | tee -a "$LOG"
bcftools stats "$MAT/raw.vcf.gz" > "$MAT/_stats_standard.txt" 2>&1
grep -E "^SN|^TSTV|^SiS" "$MAT/_stats_standard.txt" | tee -a "$LOG"

echo "=== bcftools stats: GVCF 模式 ===" | tee -a "$LOG"
bcftools stats "$MAT/raw.g.vcf.gz" > "$MAT/_stats_gvcf.txt" 2>&1
grep -E "^SN|^TSTV" "$MAT/_stats_gvcf.txt" | tee -a "$LOG"

echo "=== NON_REF 检查（GVCF 符号等位）===" | tee -a "$LOG"
grep -c "<NON_REF>" "$MAT/raw.g.vcf.gz" 2>/dev/null || echo "0"
zgrep -m1 "<NON_REF>" "$MAT/raw.g.vcf.gz" 2>/dev/null | head -1 | tee -a "$LOG"

echo "=== done $(date -u) ===" | tee -a "$LOG"
