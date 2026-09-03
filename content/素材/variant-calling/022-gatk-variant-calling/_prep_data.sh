#!/usr/bin/env bash
# 数据准备：把用到的 reference + BAM 复制进素材目录，生成 .fai；并探一下 BAM 里有无真实变异信号。
set -u
source /opt/miniconda3/etc/profile.d/conda.sh

MAT=/mnt/d/1.WorkDir/RedBook/content/素材/variant-calling/022-gatk-variant-calling
SRC=/mnt/d/1.WorkDir/RedBook/content/素材/read-alignment/011-bowtie2-alignment
LOG=$MAT/_prep_data.log
: > "$LOG"

echo "=== copy data ===" | tee -a "$LOG"
cp -v "$SRC/reference.fa"        "$MAT/reference.fa"        | tee -a "$LOG"
cp -v "$SRC/aligned_e2e.bam"     "$MAT/aligned_e2e.bam"     | tee -a "$LOG"
cp -v "$SRC/aligned_e2e.bam.bai" "$MAT/aligned_e2e.bam.bai" | tee -a "$LOG"

echo "=== faidx (bio env samtools) ===" | tee -a "$LOG"
conda run -n bio bash -lc "samtools faidx $MAT/reference.fa" 2>&1 | tee -a "$LOG"
ls -la "$MAT/reference.fa.fai" | tee -a "$LOG"

echo "=== BAM header ===" | tee -a "$LOG"
conda run -n bio bash -lc "samtools view -H $MAT/aligned_e2e.bam" 2>&1 | tee -a "$LOG"

echo "=== flagstat ===" | tee -a "$LOG"
conda run -n bio bash -lc "samtools flagstat $MAT/aligned_e2e.bam" 2>&1 | tee -a "$LOG"

echo "=== mpileup: 统计有覆盖的碱基列数 & 非纯合参考的候选变异列 ===" | tee -a "$LOG"
conda run -n bio bash -lc "samtools mpileup -f $MAT/reference.fa $MAT/aligned_e2e.bam 2>/dev/null" > "$MAT/_mpileup.txt" 2>/dev/null
TOTAL_COLS=$(wc -l < "$MAT/_mpileup.txt")
# 含 . , 之外字符（即非纯参考）的列视为候选变异
VAR_COLS=$(awk 'length($5)>0 && $5 !~ /^[.]*$/ {c++} END{print c+0}' "$MAT/_mpileup.txt")
echo "total pileup cols (covered): $TOTAL_COLS" | tee -a "$LOG"
echo "cols with non-ref base(s) (candidate variants): $VAR_COLS" | tee -a "$LOG"

echo "=== done ===" | tee -a "$LOG"
