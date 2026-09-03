#!/usr/bin/env bash
# 045 bed-file-basics: real-run driver (WSL Ubuntu, conda env bio)
# All tool output goes to files in this directory; full transcript -> _run.log
set -u
DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"
exec > _run.log 2>&1

source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio

echo "=== 0) versions ==="
bedtools --version
samtools --version | head -2
python3 --version
python3 -c "import pybedtools; print('pybedtools', pybedtools.__version__)" 2>&1 | tail -1
python3 -c "import pyranges; print('pyranges', pyranges.__version__)" 2>&1 | tail -1

echo "=== 1) inputs (make_inputs.py, seed=45) ==="
python3 make_inputs.py

echo "=== 2) genome.txt from the SAME reference FASTA ==="
samtools faidx chr1.fa
cut -f1,2 chr1.fa.fai > genome.txt
cat genome.txt

echo "=== 3) validation: field count / inverted / chrom names ==="
for f in genes.bed exons.bed cpg.bed peaks.bed transcripts.bed12; do
  echo "-- $f field counts (expect one unique value):"
  awk -F'\t' '{print NF}' "$f" | sort -u
done
echo "-- negative or inverted intervals (expect none):"
awk -F'\t' '$2 < 0 || $2 >= $3' genes.bed exons.bed cpg.bed
echo "(rc=$?)"
echo "-- chrom names per file:"
for f in genes.bed exons.bed cpg.bed; do
  echo -n "$f: "; cut -f1 "$f" | sort -u | tr '\n' ' '; echo
done

echo "=== 4) sort: bedtools sort vs coreutils sort -k1,1 -k2,2n ==="
bedtools sort -i cpg.bed > cpg.sorted.bed
bedtools sort -i exons.bed > exons.sorted.bed
bedtools sort -i genes.bed > genes.sorted.bed
sort -k1,1 -k2,2n cpg.bed > cpg.coreutils.bed
if diff -q cpg.sorted.bed cpg.coreutils.bed > /dev/null; then
  echo "bedtools sort == sort -k1,1 -k2,2n : IDENTICAL"
else
  echo "bedtools sort == sort -k1,1 -k2,2n : DIFFER"
  diff cpg.sorted.bed cpg.coreutils.bed | head -5
fi
tac cpg.sorted.bed > cpg.unsorted.bed
wc -l cpg.bed cpg.sorted.bed cpg.unsorted.bed

echo "=== 5) slop needs genome.txt; boundary clipping at chromosome ends ==="
printf 'chr1\t0\t100\tedge_left\t500\t.\nchr1\t1999500\t2000000\tedge_right\t500\t.\n' > edge.bed
bedtools slop -i edge.bed -g genome.txt -b 500 > edge.slop.bed
cat edge.bed edge.slop.bed
bedtools slop -i genes.sorted.bed -g genome.txt -b 500 > genes.slop500.bed
wc -l genes.slop500.bed
awk -F'\t' '$2 == 0 || $3 == 2000000' edge.slop.bed genes.slop500.bed > slop_clipped.txt
wc -l slop_clipped.txt

echo "=== 6) flank: 1 kb upstream promoters, overlap CpG ==="
bedtools flank -i genes.sorted.bed -g genome.txt -l 1000 -r 0 > promoters.bed
wc -l promoters.bed
bedtools intersect -a promoters.bed -b cpg.sorted.bed -c > promoters_cpg.txt
awk -F'\t' '$7 >= 1' promoters_cpg.txt | wc -l
bedtools intersect -a promoters.bed -b cpg.sorted.bed -wao > promoters_cpg_pairs.txt
wc -l promoters_cpg_pairs.txt

echo "=== 7) merge: only overlapping/book-ended intervals collapse ==="
bedtools merge -i cpg.sorted.bed > cpg.merged.bed
bedtools merge -i cpg.sorted.bed -d 200 > cpg.merged.d200.bed
bedtools merge -i genes.sorted.bed > genes.merged.bed
bedtools merge -i exons.sorted.bed > exons.merged.bed
wc -l cpg.bed cpg.merged.bed cpg.merged.d200.bed genes.bed genes.merged.bed exons.bed exons.merged.bed
awk -F'\t' '{s += $3 - $2} END {print "cpg_merged_bp=" s}' cpg.merged.bed
awk -F'\t' '{s += $3 - $2} END {print "cpg_merged_d200_bp=" s}' cpg.merged.d200.bed

echo "=== 8) intersect: exons vs CpG islands ==="
bedtools intersect -a exons.sorted.bed -b cpg.sorted.bed -c > exons_cpg.txt
awk -F'\t' '$7 >= 1' exons_cpg.txt | wc -l
bedtools intersect -a exons.sorted.bed -b cpg.sorted.bed > exons_cpg_pairs.bed
wc -l exons_cpg_pairs.bed
awk -F'\t' '{s += $3 - $2} END {print s}' exons_cpg_pairs.bed

echo "=== 9) complement: coverage identity vs genome ==="
bedtools complement -i cpg.sorted.bed -g genome.txt > non_cpg.bed
wc -l non_cpg.bed
echo "-- identity check: cpg_merged_bp + non_cpg_bp must equal 2000000"
awk -F'\t' '{s += $3 - $2} END {print "cpg_merged_bp=" s}' cpg.merged.bed
awk -F'\t' '{s += $3 - $2} END {print "non_cpg_bp=" s}' non_cpg.bed
bedtools complement -i genes.sorted.bed -g genome.txt > non_genes.bed
awk -F'\t' '{s += $3 - $2} END {print "genes_bp=" s}' genes.sorted.bed
awk -F'\t' '{s += $3 - $2} END {print "non_genes_bp=" s}' non_genes.bed

echo "=== 9b) shuffle: random relocation, also needs genome.txt ==="
bedtools shuffle -i cpg.sorted.bed -g genome.txt -seed 45 > cpg.shuffled.bed
wc -l cpg.shuffled.bed
head -3 cpg.shuffled.bed
bedtools intersect -a cpg.sorted.bed -b cpg.shuffled.bed -u | wc -l

echo "=== 10) makewindows: fixed vs sliding ==="
bedtools makewindows -g genome.txt -w 10000 > windows_fixed.bed
bedtools makewindows -g genome.txt -w 10000 -s 5000 -i winnum > windows_slide.bed
wc -l windows_fixed.bed windows_slide.bed
tail -2 windows_slide.bed
bedtools makewindows -g genome.txt -w 100000 > bins100k.bed
echo "-- per-100kb-bin feature counts -> per_bin.tsv"
bedtools intersect -a bins100k.bed -b genes.sorted.bed -c | cut -f4 > bin_genes.txt
bedtools intersect -a bins100k.bed -b exons.sorted.bed -c | cut -f4 > bin_exons.txt
bedtools intersect -a bins100k.bed -b cpg.sorted.bed -c | cut -f4 > bin_cpg.txt
echo -e "bin\tbin_start\tbin_end\tgenes\texons\tcpg" > per_bin.tsv
paste bin_genes.txt bin_exons.txt bin_cpg.txt \
  | awk -F'\t' '{printf "%d\t%d\t%d\t%d\t%d\t%d\n", NR, (NR-1)*100000, NR*100000, $1, $2, $3}' >> per_bin.tsv
cat per_bin.tsv

echo "=== 11) getfasta: extract sequences, verify len == end - start ==="
bedtools getfasta -fi chr1.fa -bed exons.sorted.bed -name -fo exons.fa
grep -c '^>' exons.fa
python3 - <<'PY'
recs = {}
name, seq = None, []
with open('exons.fa') as f:
    for line in f:
        if line.startswith('>'):
            if name: recs[name] = ''.join(seq)
            name, seq = line[1:].strip(), []
        else:
            seq.append(line.strip())
    if name: recs[name] = ''.join(seq)
vals = list(recs.values())
n = lens_ok = i = 0
for line in open('exons.sorted.bed'):
    c = line.rstrip('\n').split('\t')
    s, e = int(c[1]), int(c[2])
    n += 1
    if len(vals[i]) == e - s:
        lens_ok += 1
    else:
        print('LEN MISMATCH', i, s, e, len(vals[i]))
    i += 1
print('getfasta_records=%d len_eq_end_minus_start=%d' % (n, lens_ok))
first_name = list(recs)[0]
print('first=%s len=%d' % (first_name, len(recs[first_name])))
PY

echo "=== 12) VCF(1-based) -> BED(0-based): start-1, end unchanged ==="
grep -v '^#' variants.vcf | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2}' > variants.bed
wc -l variants.bed
head -3 variants.bed
echo "-- 1-bp round trip: GFF 1000-1000 (1-based closed) -> BED 999-1000"
awk 'BEGIN{OFS="\t"} {print $1, $4-1, $4}' landmark.gff > landmark.bed
cat landmark.bed
echo "-- base identity: samtools faidx chr1:1000-1000 vs VCF POS=1000 REF"
samtools faidx chr1.fa chr1:1000-1000 > landmark_base.txt
cat landmark_base.txt
grep -P '^chr1\t1000\t' variants.vcf
bedtools getfasta -fi chr1.fa -bed landmark.bed -fo landmark_base_bedtools.txt
cat landmark_base_bedtools.txt

echo "=== 13) bamtobed: SAM 1-based POS -> BED 0-based start (POS-1) ==="
samtools view -bS reads.sam > reads.bam
bedtools bamtobed -i reads.bam > alignments.bed
cat alignments.bed
echo "-- spliced read (50M100N50M) with -split:"
bedtools bamtobed -i reads.bam -split > spliced.bed
grep 'read2_spliced' spliced.bed

echo "=== 14) BED12 block invariants + bed12tobed6 ==="
bedtools bed12tobed6 -i transcripts.bed12 > transcripts_exons.bed6
wc -l transcripts_exons.bed6
python3 - <<'PY'
import sys
ok_all = True
n_models = 0
n_blocks = 0
for line in open('transcripts.bed12'):
    c = line.rstrip('\n').split('\t')
    s, e = int(c[1]), int(c[2])
    bs = [int(x) for x in c[11].rstrip(',').split(',')]
    sz = [int(x) for x in c[10].rstrip(',').split(',')]
    n_models += 1
    n_blocks += len(bs)
    checks = [
        bs[0] == 0,
        bs[-1] + sz[-1] == e - s,
        all(bs[i] < bs[i+1] for i in range(len(bs)-1)),
        all(bs[i] + sz[i] <= e - s for i in range(len(bs))),
        int(c[6]) >= s and int(c[7]) <= e,
    ]
    if not all(checks):
        ok_all = False
        print('INVARIANT FAIL', c[3], checks)
print('bed12_models=%d total_blocks=%d invariants_all_pass=%s'
      % (n_models, n_blocks, ok_all))
PY

echo "=== 15) narrowPeak: summit = chromStart + peak (0-based offset, col 10) ==="
awk 'BEGIN{OFS="\t"} {printf "%s\t%d\t%d\t%s\t%d\n", $1, $2, $3, $4, $2 + $10}' \
  peaks.bed > summits.tsv
cat summits.tsv
awk -F'\t' '$10 == -1 {print "peak_offset_-1 (not assigned): " $4}' peaks.bed

echo "=== 16) failure demo A: chrom-name mismatch (chr1 vs 1) is SILENT ==="
bedtools intersect -a cpg.sorted.bed -b cpg_bare.bed > mismatch_result.bed
echo "intersect lines: $(wc -l < mismatch_result.bed)  rc=$?"
cut -f1 cpg_bare.bed | sort -u

echo "=== 17) failure demo B: CRLF endings ==="
cat -A cpg_crlf.bed | head -2 > crlf_catA.txt
cat crlf_catA.txt
awk -F'\t' '{if ($NF ~ /\r$/) n++} END {print "crlf_lines=" n + 0}' cpg_crlf.bed
awk -F'\t' '{print "last_field=[" $NF "] on line " NR; exit}' cpg_crlf.bed
sed 's/\r$//' cpg_crlf.bed > cpg_crlf_fixed.bed
if diff -q cpg_crlf_fixed.bed <(head -5 cpg.sorted.bed) > /dev/null; then
  echo "CRLF-fixed file == first 5 lines of cpg.sorted.bed : IDENTICAL"
fi

echo "=== 18) failure demo C: -sorted with unsorted input errors out ==="
bedtools intersect -a cpg.unsorted.bed -b exons.sorted.bed -sorted \
  > sorted_demo.out 2> sorted_demo.err
echo "rc=$?"
cat sorted_demo.err
echo "-- same pair without -sorted (in-memory path tolerates any order):"
bedtools intersect -a cpg.unsorted.bed -b exons.sorted.bed | wc -l

echo "=== 19) parse results into bed_results.json ==="
python3 _parse.py

echo "=== 20) done ==="
ls -la
