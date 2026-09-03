#!/usr/bin/env bash
# 043 contamination-screening: fetch three small real genomes from NCBI/ENA
set -uo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$DIR/refs"
cd "$DIR/refs"

ok() { grep -q '^>' "$1" 2>/dev/null && grep -q '[ACGT]' "$1"; }

fetch() {
  local acc="$1" out="$2"
  if ok "$out"; then echo "-- $out already valid, skip"; return 0; fi
  for attempt in 1 2; do
    echo "--- fetching $acc -> $out (attempt $attempt)"
    curl -s --max-time 600 \
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$acc&rettype=fasta&retmode=text" \
      -o "$out"
    if ok "$out"; then return 0; fi
    sleep 3
  done
  return 1
}

fetch NC_000913.3 ecoli.fna || exit 1    # Escherichia coli K-12 MG1655, 4.64 Mb (sample of origin)
fetch NC_001422.1 phix.fna  || exit 1    # Escherichia phage phiX174, 5.4 kb (contaminant A; .2 flaky via eutils, .1 verified)
fetch NC_001416.1 lambda.fna || exit 1   # Enterobacteria phage lambda, 48.5 kb (contaminant B)

for f in ecoli.fna phix.fna lambda.fna; do
  echo "== $f: $(grep -c '>' $f) record(s), $(wc -c < $f) bytes, $(grep -v '>' $f | tr -d '\n' | wc -c) bp"
  grep '>' $f
done
