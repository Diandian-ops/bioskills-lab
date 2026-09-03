#!/usr/bin/env bash
# 043 contamination-screening: build mini real-genome kraken2 DB + bowtie2 panel
# NOTE: kraken2-build --download-taxonomy uses rsync (blocked here); we instead
# fetch taxdump over HTTPS and hand-build seqid2taxid.map for our 3 accessions.
set -euo pipefail
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio-qc
DIR="$(cd "$(dirname "$0")" && pwd)"
DB="$DIR/db/k2_mini"

# NOTE: keep db/taxdump.tar.gz if already downloaded (18 MB); only wipe build dirs
rm -rf "$DB" "$DIR/db/bt2"
mkdir -p "$DB/taxonomy" "$DB/library" "$DIR/db/bt2"

echo "=== taxids from GenBank headers (eutils)"
taxid_for() {
  curl -s --max-time 300 \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$1&rettype=gb&retmode=text" \
    | head -c 200000 | grep -oE 'taxon:[0-9]+' | head -1 | cut -d: -f2
}
TID_ECOLI=$(taxid_for NC_000913.3)
TID_PHIX=$(taxid_for NC_001422.1)
TID_LAMBDA=$(taxid_for NC_001416.1)
echo "ecoli=$TID_ECOLI phix=$TID_PHIX lambda=$TID_LAMBDA"
[ -n "$TID_ECOLI" ] && [ -n "$TID_PHIX" ] && [ -n "$TID_LAMBDA" ]

echo "=== seqid2taxid.map"
printf 'NC_000913.3\t%s\n' "$TID_ECOLI" >  "$DB/seqid2taxid.map"
printf 'NC_001422.1\t%s\n' "$TID_PHIX" >> "$DB/seqid2taxid.map"
printf 'NC_001416.1\t%s\n' "$TID_LAMBDA" >> "$DB/seqid2taxid.map"
cat "$DB/seqid2taxid.map"

echo "=== taxonomy subset via eutils (real NCBI lineage)"
# full taxdump.tar.gz (75 MB) downloads at ~18 KB/s on this network (>1 h);
# eutils Taxonomy API is fast and provides the same real lineage data, so we
# build a complete-ancestry subset covering exactly our 3 taxids.
python3 "$DIR/build_taxsubset.py"
ls -la "$DB/taxonomy/nodes.dmp" "$DB/taxonomy/names.dmp"

echo "=== kraken2-build: add libraries"
kraken2-build --db "$DB" --add-to-library "$DIR/refs/ecoli.fna"
kraken2-build --db "$DB" --add-to-library "$DIR/refs/phix.fna"
kraken2-build --db "$DB" --add-to-library "$DIR/refs/lambda.fna"

echo "=== kraken2-build: build"
kraken2-build --db "$DB" --build --threads 4 --kmer-len 35 --minimizer-len 31

echo "=== bowtie2 indexes for FastQ Screen"
for g in ecoli phix lambda; do
  bowtie2-build --threads 4 "$DIR/refs/$g.fna" "$DIR/db/bt2/$g" > /dev/null
done

BOWTIE2_BIN="$(which bowtie2)"
printf 'BOWTIE2\t%s\nTHREADS\t4\nDATABASE\tEcoli\t%s\nDATABASE\tPhiX\t%s\nDATABASE\tLambda\t%s\n' \
  "$BOWTIE2_BIN" "$DIR/db/bt2/ecoli" "$DIR/db/bt2/phix" "$DIR/db/bt2/lambda" > "$DIR/fastq_screen.conf"
cat "$DIR/fastq_screen.conf"

echo "=== kraken2-inspect (top)"
kraken2-inspect --db "$DB" | head -12
echo "=== DB sizes"
du -sh "$DB" "$DIR/db/bt2"
