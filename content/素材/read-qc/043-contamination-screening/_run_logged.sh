#!/usr/bin/env bash
# wrapper: run _run.sh and tee transcript to repro_transcript.txt
DIR="$(cd "$(dirname "$0")" && pwd)"
bash "$DIR/_run.sh" > "$DIR/repro_transcript.txt" 2>&1
echo "EXIT=$?"
