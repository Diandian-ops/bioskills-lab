#!/usr/bin/env bash
set -e
cd /mnt/d/1.WorkDir/RedBook/content/素材/alignment-files/028-sam-bam-basics
BAM=aligned_e2e.bam
REF=reference.fa

echo "===== 028 sam-bam-basics repro transcript ====="
echo ""
echo "## samtools view -H (header)"
echo "$ samtools view -H $BAM"
samtools view -H $BAM | head -10
echo ""
echo "## samtools view (first 2 alignments)"
echo "$ samtools view $BAM | head -2"
samtools view $BAM | head -2
echo ""
echo "## samtools view -c (count)"
echo "$ samtools view -c $BAM"
samtools view -c $BAM
echo ""
echo "## samtools flagstat"
echo "$ samtools flagstat $BAM"
samtools flagstat $BAM
echo ""
echo "## samtools flags decode examples"
for f in 0 4 16 99 147 256 2048; do
  echo "$ samtools flags $f"
  samtools flags $f
  echo ""
done
echo "## region query (contig1:1-100)"
echo "$ samtools view $BAM contig1:1-100 | head -3"
samtools view $BAM contig1:1-100 | head -3
echo ""
echo "## BAM -> SAM -> BAM roundtrip"
echo "$ samtools view -h -o demo.sam $BAM && samtools view -c demo.sam"
samtools view -h -o demo.sam $BAM
samtools view -c demo.sam
echo "$ samtools view -b -o demo_back.bam demo.sam && samtools view -c demo_back.bam"
samtools view -b -o demo_back.bam demo.sam
samtools view -c demo_back.bam
rm -f demo.sam demo_back.bam
echo ""
echo "## BAM -> CRAM -> BAM roundtrip"
echo "$ samtools view -C -T $REF -o demo.cram $BAM && samtools view -c demo.cram"
samtools view -C -T $REF -o demo.cram $BAM
samtools view -c demo.cram
echo "$ samtools view -b -T $REF -o demo_back2.bam demo.cram && samtools view -c demo_back2.bam"
samtools view -b -T $REF -o demo_back2.bam demo.cram
samtools view -c demo_back2.bam
rm -f demo.cram demo_back2.bam
