#!/usr/bin/env bash
# 048 coverage-analysis: 真跑脚本（bwa + samtools + bedtools，全部按 SKILL.md 口径）
set -e
cd "$(dirname "$0")"
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio
exec > _run.log 2>&1

echo "== versions =="
bedtools --version
samtools --version | head -2
bwa 2>&1 | grep -i "^Version" || true
echo "mosdepth: $(which mosdepth 2>/dev/null || echo NOT INSTALLED)"

echo "== [1] input data ="
python3 make_inputs.py

echo "== [2] index + align ="
bwa index ref.fa > /dev/null
samtools faidx ref.fa
bwa mem -t 4 -R "@RG\tID:rg1\tSM:sample1\tPL:ILLUMINA\tLB:lib1" \
    ref.fa reads_1.fq.gz reads_2.fq.gz 2> _bwa.log \
    | samtools sort -o aln.bam -
samtools index aln.bam
samtools flagstat aln.bam

echo "== [3] bedtools genomecov (SKILL.md 口径) ="
# -bga: bedGraph 含零覆盖区段
bedtools genomecov -ibam aln.bam -bga > cov.bedGraph
wc -l cov.bedGraph
# 裸默认 = 5 列直方图（不是 track）
bedtools genomecov -ibam aln.bam > cov_hist.txt
head -5 cov_hist.txt
tail -3 cov_hist.txt
# -pc: 片段口径（mate 重叠只计一次）
bedtools genomecov -ibam aln.bam -pc -bg > frag.bedGraph
wc -l frag.bedGraph

echo "== [4] bedtools coverage: -a targets(BED) -b reads(BAM) ="
# targets.bed：0-based half-open
printf "chrS\t100000\t200000\tBG\nchrS\t500000\t550000\tHIGH\nchrS\t1000000\t1030000\tLOW\nchrS\t1500000\t1520000\tZERO\n" > targets.bed
bedtools coverage -a targets.bed -b aln.bam > per_target.bed
cat per_target.bed
wc -l per_target.bed
bedtools coverage -a targets.bed -b aln.bam -mean > per_target_mean.bed
cat per_target_mean.bed

echo "== [5] samtools coverage / depth ="
samtools coverage aln.bam | tee samtools_coverage.txt
# 逐碱基（全染色体，含 -a 零深度位）
samtools depth -a aln.bam | gzip > depth_all.txt.gz
# HIGH 区：朴素 vs -s（重叠 mate 只计一次）
samtools depth -a -r chrS:500001-550000 aln.bam > depth_high_naive.txt
samtools depth -a -s -r chrS:500001-550000 aln.bam > depth_high_s.txt
# ZERO 区：-a 应输出全零
samtools depth -a -r chrS:1500001-1520000 aln.bam > depth_zero.txt
wc -l depth_high_naive.txt depth_high_s.txt depth_zero.txt
awk '$3>0' depth_zero.txt | wc -l

echo "== done ="
