#!/usr/bin/env bash
# 052 proximity-operations: real bedtools run (closest/window/slop/flank)
set -x
cd /mnt/d/1.WorkDir/RedBook/content/素材/genome-intervals/052-proximity-operations
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio

bedtools --version > _versions.txt 2>&1
python3 -c "import pybedtools" >/dev/null 2>_pybedtools_check.txt \
  && echo "pybedtools AVAILABLE" >> _pybedtools_check.txt \
  || echo "pybedtools MISSING -> CLI only (noted honestly in notes)" >> _pybedtools_check.txt

python3 make_inputs.py > _make_inputs.log 2>&1

bedtools sort -i genes.bed > genes.sorted.bed
bedtools sort -i peaks.bed > peaks.sorted.bed

# --- closest: skill main recipe (signed by gene strand, ignore overlaps, one row/peak)
bedtools closest -a peaks.sorted.bed -b genes.sorted.bed -D b -io -t first > nearest_db_io_first.bed
# --- closest: default -t all with -d (tie double-counting demo)
bedtools closest -a peaks.sorted.bed -b genes.sorted.bed -d -t all > nearest_all_d.bed
# --- closest: -D ref (coordinate sign only; mis-signs minus-strand genes)
bedtools closest -a peaks.sorted.bed -b genes.sorted.bed -D ref -io -t first > nearest_ref_io_first.bed
# --- closest: k=3 nearest with unsigned distance
bedtools closest -a peaks.sorted.bed -b genes.sorted.bed -k 3 -d > top3_k3_d.bed

# --- window: all genes within +/-50 kb, counted per peak
bedtools window -a peaks.sorted.bed -b genes.sorted.bed -w 50000 -c > window_counts.bed

# --- promoters: collapse genes to TSS (skill awk), then strand-aware slop
awk -v OFS='\t' '{ if ($6=="+") print $1,$2,$2+1,$4,$5,$6; else print $1,$3-1,$3,$4,$5,$6 }' genes.bed > tss.bed
bedtools slop -i tss.bed -g genome.txt -s -l 2000 -r 200 > promoters.bed
# --- wrong-promoter contrast: slop -b 2000 on gene BODIES
bedtools slop -i genes.bed -g genome.txt -b 2000 > genebody_slop_b2000.bed
# --- flank: regions beside each gene, body dropped
bedtools flank -i genes.bed -g genome.txt -s -b 1000 > gene_flanks.bed

python3 _pybed_check.py > _pybed_check.log 2>&1
python3 parse_results.py > _parse.log 2>&1
echo "RUN DONE exit=$?"
