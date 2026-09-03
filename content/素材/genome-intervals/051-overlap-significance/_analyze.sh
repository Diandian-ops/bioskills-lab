#!/usr/bin/env bash
set -euo pipefail
DIR="/mnt/d/1.WorkDir/RedBook/content/素材/genome-intervals/051-overlap-significance"
cd "$DIR"
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio
echo -e "set\tjaccard" | tee observed_jaccard.tsv
for S in enriched random; do
  J=$(bedtools jaccard -a A_${S}.bed -b B_features.bed | tail -1 | cut -f3)
  echo -e "${S}\t${J}" | tee -a observed_jaccard.tsv
done
echo "===== analyze ====="
python3 _analyze.py
