#!/usr/bin/env bash
set -euo pipefail
DIR="/mnt/d/1.WorkDir/RedBook/content/素材/genome-intervals/051-overlap-significance"
cd "$DIR"
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio
echo "== B union bases =="
bedtools merge -i B_features.bed | awk '{s+=$3-$2} END {print s}'
echo "== workspace bases =="
awk '{s+=$3-$2} END {print s}' workspace.bed
echo "== A_random intersect B bases =="
bedtools intersect -a A_random.bed -b B_features.bed | awk '{s+=$3-$2} END {print s}'
echo "== A_enriched intersect B bases =="
bedtools intersect -a A_enriched.bed -b B_features.bed | awk '{s+=$3-$2} END {print s}'
echo "== how much of shuffled A lands outside workspace (300 reps, matched) =="
: > diag_outside.tsv
for i in $(seq 1 300); do
  bedtools shuffle -i A_random.bed -g genome.txt -incl workspace.bed -chrom -seed "$i" \
    | bedtools intersect -a - -b workspace.bed -v \
    | awk '{s+=$3-$2} END {print s+0}' >> diag_outside.tsv
done
python3 - <<'EOF'
vals = [int(x) for x in open('diag_outside.tsv')]
n = len(vals)
mean = sum(vals)/n
print("shuffled-A bases outside workspace: mean %.1f bp of 180000 (%.2f%%), nonzero %d/%d"
      % (mean, 100*mean/180000, sum(1 for v in vals if v>0), n))
EOF
echo "== shuffle placement vs tile placement: 300 shuffles of A_random, mean intersect with B =="
: > diag_inter.tsv
for i in $(seq 1 300); do
  bedtools shuffle -i A_random.bed -g genome.txt -incl workspace.bed -chrom -seed "$i" \
    | LC_ALL=C sort -k1,1 -k2,2n \
    | bedtools intersect -a - -b B_features.bed \
    | awk '{s+=$3-$2} END {print s+0}' >> diag_inter.tsv
done
python3 - <<'EOF'
vals = [int(x) for x in open('diag_inter.tsv')]
mean = sum(vals)/len(vals)
print("shuffle->B inter: mean %.0f bp (obs A_random->B = 36038-equivalent below)" % mean)
EOF
bedtools intersect -a A_random.bed -b B_features.bed | awk '{s+=$3-$2} END {print "obs A_random->B inter:", s}'
