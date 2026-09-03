#!/usr/bin/env bash
# _run.sh -- 049 gtf-gff-handling real run (gffread 0.12.9 / bedtools 2.31.1)
set -u
D="/mnt/d/1.WorkDir/RedBook/content/素材/genome-intervals/049-gtf-gff-handling"
cd "$D" || exit 1
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio

: > _run.log
exec > >(tee -a _run.log) 2>&1

echo "### date: $(date '+%Y-%m-%d %H:%M:%S')"
echo "### versions"
gffread --version
bedtools --version
python3 --version

echo "### [1] generate inputs"
python3 make_inputs.py

echo "### [2] conversions"
echo "-- GFF3 -> GTF (-T)"
gffread annotation.gff3 -T -o conv_from_gff3.gtf
echo "exit=$?"
echo "-- GTF -> GFF3 (default output)"
gffread annotation.gtf -o conv_from_gtf.gff3
echo "exit=$?"
head -3 conv_from_gff3.gtf
echo "--"
head -4 conv_from_gtf.gff3

echo "### [3] extract from GTF (seqid chr1 == genome)"
gffread -w tx_from_gtf.fa -g genome.fa annotation.gtf
gffread -x cds_from_gtf.fa -g genome.fa annotation.gtf
gffread -y prot_from_gtf.fa -g genome.fa annotation.gtf
grep -c '>' tx_from_gtf.fa cds_from_gtf.fa prot_from_gtf.fa

echo "### [4] extract from GFF3 as-is (seqid '1' vs genome 'chr1')"
gffread -w tx_from_gff3_asis.fa -g genome.fa annotation.gff3
echo "exit=$?"
if [ -f tx_from_gff3_asis.fa ]; then grep -c '>' tx_from_gff3_asis.fa; wc -c tx_from_gff3_asis.fa; fi

echo "### [5] BED per skill: start-1 only (+ over-corrected control), seqid remap"
python3 make_bed.py
bedtools getfasta -fi genome.fa -bed exons.bed -s -name -fo exons_bed_sense.fa
bedtools getfasta -fi genome.fa -bed exons_overcorrected.bed -s -name -fo exons_bed_overcor.fa

echo "### [6] extract from renamed GFF3 (seqid 1 -> chr1)"
gffread -w tx_from_gff3.fa -g genome.fa annotation_chr1.gff3
gffread -x cds_from_gff3.fa -g genome.fa annotation_chr1.gff3
gffread -y prot_from_gff3.fa -g genome.fa annotation_chr1.gff3
grep -c '>' tx_from_gff3.fa cds_from_gff3.fa prot_from_gff3.fa

echo "### [7] audit: 9-column walk, phase, namespaces, consistency"
python3 audit.py

echo "### [8] file inventory"
ls -la
echo "### done"
