#!/usr/bin/env bash
# download taxdump.tar.gz with resume + retries
set -uo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
OUT="$DIR/db/taxdump.tar.gz"
mkdir -p "$DIR/db"
URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
for i in 1 2 3 4 5 6 7 8; do
  echo "-- attempt $i (resume)"
  curl -sS -C - --retry 3 --retry-delay 2 --max-time 900 "$URL" -o "$OUT"
  sz=$(wc -c < "$OUT")
  echo "   size now: $sz"
  if tar -tzf "$OUT" > /dev/null 2>&1; then
    echo "OK: taxdump.tar.gz complete ($sz bytes)"
    exit 0
  fi
  sleep 3
done
echo "FAILED"
exit 1
