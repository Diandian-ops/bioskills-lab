#!/usr/bin/env bash
# 047 bigwig-tracks real run: bedGraph -> bigWig -> bigWigInfo -> bigWigAverageOverBed
set -u
cd "$(dirname "$0")"
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio

echo "== tool versions =="
bedGraphToBigWig 2>&1 | head -1
bigWigInfo 2>&1 | head -1
bigWigAverageOverBed 2>&1 | head -1
bedtools --version
python3 --version

echo "== step 1: make simulated inputs (seed=42) =="
python3 make_inputs.py

echo "== step 2: sort bedGraph (skill requirement) =="
export LC_ALL=C
sort -k1,1 -k2,2n sim.bedGraph > sim.sorted.bedGraph
sort -c -k1,1 -k2,2n sim.sorted.bedGraph && echo "sorted OK"
wc -l sim.bedGraph sim.sorted.bedGraph chrom.sizes regions.bed

echo "== step 3: build bigWig =="
bedGraphToBigWig sim.sorted.bedGraph chrom.sizes track.bw
echo "bedGraphToBigWig exit: $?"
ls -la sim.bedGraph sim.sorted.bedGraph track.bw

echo "== step 4: bigWigInfo (header sanity check) =="
bigWigInfo track.bw | tee _bigwiginfo.txt
echo "-- bigWigInfo -chroms --"
bigWigInfo -chroms track.bw | tee _bigwiginfo_chroms.txt || echo "bigWigInfo -chroms failed"

echo "== step 5: bigWigAverageOverBed (per-region means) =="
bigWigAverageOverBed track.bw regions.bed regions_avg.tab
echo "bigWigAverageOverBed exit: $?"
cat regions_avg.tab

echo "== step 6: reconcile against design truth =="
python3 reconcile.py
echo "reconcile exit: $?"

echo "== done =="
