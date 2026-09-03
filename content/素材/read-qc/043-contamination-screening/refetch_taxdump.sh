#!/usr/bin/env bash
# re-download corrupted taxdump.tar.gz (previous partial download was truncated)
set -euo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
TAR="$DIR/db/taxdump.tar.gz"
rm -f "$TAR"
echo "expected size via HEAD:"
curl -sI --max-time 60 "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz" | tr -d "\r" | grep -i "content-length" || true
curl -s --max-time 1800 "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz" -o "$TAR"
echo "downloaded bytes: $(stat -c %s "$TAR")"
gzip -t "$TAR" && echo "GZIP OK"
tar -tzf "$TAR" > /dev/null && echo "TAR OK"
