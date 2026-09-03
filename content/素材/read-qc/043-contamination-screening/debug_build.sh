#!/usr/bin/env bash
set -uo pipefail
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio-qc
DIR="/mnt/d/1.WorkDir/RedBook/content/素材/read-qc/043-contamination-screening"
LIBEXEC=/opt/miniconda3/envs/bio-qc/share/kraken2-2.17.1/libexec
TT="$DIR/debug_tt"
cd "$DIR/db/k2_mini"
SEED=11111111111111111111111111111111110011001100110011001100110011

echo "== T1: fixed real taxonomy + real map + tiny fasta (taxid 511145)"
printf "NC_TEST.1\t511145\n" > "$TT/m5.map"
printf ">NC_TEST.1 dummy\n%s\n" "$(python3 -c 'print("ACGT"*40)')" > "$TT/test5.fna"
cat "$TT/test5.fna" | "$LIBEXEC/build_db" -k 35 -l 31 -S $SEED -H "$TT/t1.h" -t "$TT/t1.t" -o "$TT/t1.o" -n taxonomy/ -m "$TT/m5.map" -c 1000 -p 1 2>&1

echo "== T2: fixed real taxonomy + real map + REAL ecoli.fna"
cat library/added/sUiaAcWP4b.fna | "$LIBEXEC/build_db" -k 35 -l 31 -S $SEED -H "$TT/t2.h" -t "$TT/t2.t" -o "$TT/t2.o" -n taxonomy/ -m seqid2taxid.map -c 2159542 -p 4 2>&1
ls -la "$TT/t2.h"
