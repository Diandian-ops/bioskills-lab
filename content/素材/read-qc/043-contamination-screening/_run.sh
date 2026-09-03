#!/usr/bin/env bash
# 043 contamination-screening: full real-run pipeline
# kraken2 + FastQ Screen on three designed mixture samples
set -euo pipefail
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio-qc
DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"
mkdir -p kraken2 screen

echo "########## STAGE 0: environment"
kraken2 --version | head -1
fastq_screen --version 2>&1 | head -1
multiqc --version
bowtie2 --version | head -1
python3 --version

echo "########## STAGE 1: fetch real genomes"
bash "$DIR/fetch_genomes.sh"

echo "########## STAGE 2: generate designed mixture reads"
python3 "$DIR/make_inputs.py"

echo "########## STAGE 3: build mini kraken2 DB + bowtie2 panel"
bash "$DIR/build_db.sh"

echo "########## STAGE 4: kraken2 main runs (confidence 0.1, min-hit-groups 2)"
for s in S1 S2 S3; do
  echo "--- kraken2 $s"
  kraken2 --db "$DIR/db/k2_mini" --threads 4 --confidence 0.1 --minimum-hit-groups 2 \
    --paired --use-names --gzip-compressed \
    --report "kraken2/${s}.kreport" --output "kraken2/${s}.kraken" \
    "reads/${s}_1.fastq.gz" "reads/${s}_2.fastq.gz" | tail -4
done

echo "########## STAGE 5: kraken2 confidence sweep (S1, S2)"
for s in S1 S2; do
  for c in 0.0 0.05 0.1 0.2; do
    echo "--- kraken2 $s confidence $c"
    kraken2 --db "$DIR/db/k2_mini" --threads 4 --confidence "$c" --minimum-hit-groups 2 \
      --paired --use-names --gzip-compressed \
      --report "kraken2/${s}.conf${c}.kreport" --output /dev/null \
      "reads/${s}_1.fastq.gz" "reads/${s}_2.fastq.gz" | tail -3
  done
done

echo "########## STAGE 6: FastQ Screen"
for s in S1 S2 S3; do
  echo "--- fastq_screen $s"
  fastq_screen --conf "$DIR/fastq_screen.conf" --threads 4 --outdir screen \
    "reads/${s}_1.fastq.gz" "reads/${s}_2.fastq.gz" | tail -20
done

echo "########## STAGE 7: multiqc"
multiqc screen -o multiqc --force 2>&1 | tail -3

echo "########## STAGE 8: parse results"
python3 "$DIR/parse_results.py"

echo "########## DONE"
