#!/usr/bin/env bash
# 051 overlap-significance: real bedtools runs per SKILL.md
# Stage 1: make inputs, versions, bedtools fisher (analytic screen)
# Stage 2: permutation null (bedtools shuffle + jaccard, N=1000 per arm)
set -euo pipefail

DIR="/mnt/d/1.WorkDir/RedBook/content/素材/genome-intervals/051-overlap-significance"
cd "$DIR"

source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio

echo "===== versions ====="
date -u
bedtools --version
python3 --version

echo "===== stage 1: inputs ====="
python3 make_inputs.py

echo "===== head of inputs ====="
head -3 workspace.bed B_features.bed A_enriched.bed A_random.bed genome.txt

echo "===== stage 1: bedtools fisher (analytic 2x2 screen, whole-genome null) ====="
bedtools fisher -a A_enriched.bed -b B_features.bed -g genome.txt > fisher_enriched.txt
bedtools fisher -a A_random.bed  -b B_features.bed -g genome.txt > fisher_random.txt
echo "--- fisher A_enriched vs B ---"; cat fisher_enriched.txt
echo "--- fisher A_random  vs B ---"; cat fisher_random.txt

echo "===== observed jaccard ====="
echo -e "set\tjaccard" | tee observed_jaccard.tsv
for S in enriched random; do
  J=$(bedtools jaccard -a A_${S}.bed -b B_features.bed | tail -1 | cut -f3)
  echo -e "${S}\t${J}" | tee -a observed_jaccard.tsv
done

echo "===== stage 2: permutation null (shuffle + jaccard, N=1000 per arm) ====="
N=1000
for S in enriched random; do
  for MODE in matched uniform; do
    OUT="nulls_${MODE}_${S}.tsv"
    : > "$OUT"
    echo "arm: mode=${MODE} set=${S}"
    for i in $(seq 1 "$N"); do
      if [ "$MODE" = "matched" ]; then
        bedtools shuffle -i A_${S}.bed -g genome.txt -incl workspace.bed -chrom -seed "$i" \
          | LC_ALL=C sort -k1,1 -k2,2n \
          | bedtools jaccard -a - -b B_features.bed | tail -1 | cut -f3 >> "$OUT"
      else
        bedtools shuffle -i A_${S}.bed -g genome.txt -chrom -seed "$i" \
          | LC_ALL=C sort -k1,1 -k2,2n \
          | bedtools jaccard -a - -b B_features.bed | tail -1 | cut -f3 >> "$OUT"
      fi
    done
    echo "  done: $OUT ($(wc -l < "$OUT") values)"
  done
done

echo "===== stage 2 done ====="
date -u
wc -l nulls_*.tsv
