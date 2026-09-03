#!/usr/bin/env bash
# 046-bedgraph-handling: real-run pipeline.
# bedtools genomecov -> bedGraph -> LC_COLLATE=C sort -> bedGraphToBigWig,
# plus bigWig round-trip verification, unionbedg matrix, and four
# error-contract tests against bedGraphToBigWig.
set -u
DIR=/mnt/d/1.WorkDir/RedBook/content/素材/genome-intervals/046-bedgraph-handling
cd "$DIR"
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio

echo "=== 046 bedgraph-handling real run $(date +%F) ==="

echo "--- versions ---"
bedtools --version
samtools --version | head -1
bedGraphToBigWig 2>&1 | head -1
if command -v bigWigToBedGraph >/dev/null 2>&1; then bigWigToBedGraph 2>&1 | head -1; else echo "bigWigToBedGraph: NOT INSTALLED in this env"; fi
bigWigInfo 2>&1 | head -1
python3 -c "import pyBigWig; print('pyBigWig', pyBigWig.__version__)"
if command -v bamCoverage >/dev/null 2>&1; then bamCoverage --version; else echo "bamCoverage: NOT INSTALLED in this env (deeptools path documented, not run)"; fi

echo "--- step 1: inputs ---"
python3 make_inputs.py
samtools faidx ref.fa
cut -f1,2 ref.fa.fai > chrom.sizes
echo "chrom.sizes:"; cat chrom.sizes
samtools view -bS sampleA.sam | samtools sort -o sampleA.bam -
samtools index sampleA.bam
samtools view -bS sampleB.sam | samtools sort -o sampleB.bam -
samtools index sampleB.bam
NA=$(samtools view -c sampleA.bam)
NB=$(samtools view -c sampleB.bam)
echo "mapped reads: A=$NA  B=$NB"
samtools depth -a -r chr1 sampleA.bam | awk '{s+=$3} END{printf "sampleA mean depth = %.4fx over 2,000,000 bp\n", s/NR}'
samtools depth -a -r chr1 sampleB.bam | awk '{s+=$3} END{printf "sampleB mean depth = %.4fx over 2,000,000 bp\n", s/NR}'

echo "--- step 2: genomecov -bg vs -bga ---"
bedtools genomecov -ibam sampleA.bam -bg > A.bedgraph
bedtools genomecov -ibam sampleA.bam -bga > A_bga.bedgraph
bedtools genomecov -ibam sampleB.bam -bg > B.bedgraph
bedtools genomecov -ibam sampleB.bam -bga > B_bga.bedgraph
wc -l A.bedgraph A_bga.bedgraph B.bedgraph B_bga.bedgraph
awk '{if($4>m)m=$4} END{print "max raw coverage in A_bga.bedgraph = " m "x"}' A_bga.bedgraph
echo "bp tiled by -bga (should equal 2000000): $(awk '{s+=$3-$2} END{print s}' A_bga.bedgraph)"
echo "bp covered >0 by -bg (sample A): $(awk '{s+=$3-$2} END{print s}' A.bedgraph)"

echo "--- step 3: manual RPM scaling (-scale 1e6/nreads) ---"
SA=$(awk "BEGIN{printf \"%.10g\", 1000000/$NA}")
SB=$(awk "BEGIN{printf \"%.10g\", 1000000/$NB}")
echo "RPM scale factors: A=$SA  B=$SB"
bedtools genomecov -ibam sampleA.bam -bg -scale $SA > A_scaled.bedgraph
bedtools genomecov -ibam sampleB.bam -bg -scale $SB > B_scaled.bedgraph
head -2 A.bedgraph; head -2 A_scaled.bedgraph

echo "--- step 4: sort + bedGraphToBigWig ---"
LC_COLLATE=C sort -k1,1 -k2,2n A.bedgraph > A.sorted.bedgraph
LC_COLLATE=C sort -k1,1 -k2,2n B.bedgraph > B.sorted.bedgraph
bedGraphToBigWig A.sorted.bedgraph chrom.sizes A.bw && echo "A.bw OK"
bedGraphToBigWig B.sorted.bedgraph chrom.sizes B.bw && echo "B.bw OK"
ls -l A.bw B.bw A.sorted.bedgraph

echo "--- step 5: bigWig round-trip (pyBigWig full-chromosome readback + bigWigInfo) ---"
echo "note: bigWigToBedGraph is not installed in this env; the round-trip is verified"
echo "      with pyBigWig, which returns the underlying bedGraph rows via bw.intervals()."
bigWigInfo A.bw | head -12
python3 - <<'PYEOF'
import json
import pyBigWig

bw = pyBigWig.open("A.bw")
print("bigWig chroms:", bw.chroms())
print("bw.stats chr1:780000-820000 mean =", bw.stats("chr1", 780000, 820000, type="mean")[0])

bw_rows = bw.intervals("chr1", 0, 2000000)
bg_rows = []
with open("A.sorted.bedgraph") as f:
    for ln in f:
        c, s, e, v = ln.split()
        bg_rows.append((int(s), int(e), float(v)))

n = min(len(bw_rows), len(bg_rows))
maxd = 0.0
for (s1, e1, v1), (s2, e2, v2) in zip(bw_rows, bg_rows):
    assert s1 == s2 and e1 == e2, "interval mismatch at %d" % s1
    d = abs(v1 - v2)
    if d > maxd:
        maxd = d
print("round-trip rows compared = %d (bigWig %d vs bedGraph %d), max |delta| = %g"
      % (n, len(bw_rows), len(bg_rows), maxd))

rows = bw.intervals("chr1", 780000, 820000)
with open("fig3_bw_window.txt", "w") as out:
    for s, e, v in rows:
        out.write("%d\t%d\t%.6f\n" % (s, e, v))
print("bigWig intervals in 780000-820000 =", len(rows))
bw.close()

json.dump({"roundtrip_rows": n, "roundtrip_max_abs_diff": maxd,
           "bigwig_intervals_total": len(bw_rows)},
          open("parsed_roundtrip.json", "w"), indent=1)
PYEOF

echo "--- step 6: unionbedg multi-sample matrix ---"
bedtools unionbedg -i A.bedgraph B.bedgraph -header -names A B > matrix_union.txt
wc -l matrix_union.txt
head -3 matrix_union.txt
tail -1 matrix_union.txt

echo "--- step 7: error-contract tests against bedGraphToBigWig ---"
echo "[test 7a] unsorted input (reverse order):"
tac A.sorted.bedgraph > A.unsorted.bdg
if bedGraphToBigWig A.unsorted.bdg chrom.sizes err_unsorted.bw 2> err_unsorted.txt; then
  echo "UNEXPECTED: no error raised on unsorted input"
else
  echo "exit=$? ; message: $(head -1 err_unsorted.txt)"
fi

echo "[test 7b] chrom.sizes shorter than the data (wrong assembly proxy):"
printf 'chr1\t1000000\n' > chrom_wrong.sizes
if bedGraphToBigWig A.sorted.bedgraph chrom_wrong.sizes err_wrongsize.bw 2> err_wrongsize.txt; then
  echo "UNEXPECTED: no error raised on short chrom.sizes"
else
  echo "exit=$? ; message: $(head -1 err_wrongsize.txt)"
fi

echo "[test 7c] overlapping intervals (two tracks concatenated):"
cat A.bedgraph B.bedgraph | LC_COLLATE=C sort -k1,1 -k2,2n > AB_overlap.bdg
if bedGraphToBigWig AB_overlap.bdg chrom.sizes err_overlap.bw 2> err_overlap.txt; then
  echo "UNEXPECTED: no error raised on overlapping input"
else
  echo "exit=$? ; message: $(head -1 err_overlap.txt)"
fi
bedtools merge -i AB_overlap.bdg -d 0 -c 4 -o max > AB_merged.bdg
wc -l AB_overlap.bdg AB_merged.bdg
bedGraphToBigWig AB_merged.bdg chrom.sizes AB.bw && echo "AB.bw OK after merge (-c 4 -o max)"

echo "[test 7d] chrom naming mismatch (bedGraph says '1', chrom.sizes says 'chr1'):"
sed 's/^chr1/1/' A.sorted.bedgraph > A_named1.bdg
if bedGraphToBigWig A_named1.bdg chrom.sizes err_nomatch.bw 2> err_nomatch.txt; then
  echo "no error raised; checking what landed in the bigWig:"
  python3 - <<'PYEOF'
import pyBigWig
bw = pyBigWig.open("err_nomatch.bw")
print("bigWig chroms after name mismatch:", bw.chroms())
bw.close()
PYEOF
else
  echo "exit=$? ; message: $(head -1 err_nomatch.txt)"
fi

echo "--- step 8: file inventory ---"
ls -l *.bedgraph *.bw *.bdg matrix_union.txt 2>/dev/null
echo "=== END OF RUN ==="
