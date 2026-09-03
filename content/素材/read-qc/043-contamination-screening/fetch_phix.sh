#!/usr/bin/env bash
# fetch phiX174 via multiple fallback sources
set -uo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR/refs"

try_eutils() {
  curl -s --max-time 120 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$1&rettype=fasta&retmode=text" -o phix.fna
}

try_ena() {
  curl -s --max-time 120 "https://www.ebi.ac.uk/ena/browser/api/fasta/$1" -o phix.fna
}

ok() { grep -q '^>' phix.fna 2>/dev/null && grep -q 'ACGT' phix.fna; }

for acc in NC_001422.2 NC_001422.1 NC_001422; do
  for i in 1 2; do
    echo "-- eutils try acc=$acc attempt=$i"
    try_eutils "$acc"
    if ok; then echo "OK via eutils $acc"; head -1 phix.fna; wc -c phix.fna; exit 0; fi
    sleep 3
  done
  echo "-- ena try acc=$acc"
  try_ena "$acc"
  if ok; then echo "OK via ENA $acc"; head -1 phix.fna; wc -c phix.fna; exit 0; fi
done

echo "FAILED all sources"
exit 1
