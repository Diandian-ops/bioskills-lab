#!/usr/bin/env bash
set -uo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio
python3 - <<'EOF'
import json, subprocess
from collections import defaultdict
groups = defaultdict(list)
for ln in open("molecules.tsv").readlines()[1:]:
    p = ln.rstrip("\n").split("\t")
    if p[5]:
        groups[p[5]].append((p[1], int(p[2]), p[3], p[7]))

def keys(bam):
    out = subprocess.run("samtools view %s" % bam, shell=True, capture_output=True, text=True).stdout
    s = set()
    for ln in out.splitlines():
        p = ln.split("\t")
        flag = int(p[1])
        if not (flag & 0x40) or (flag & 0xD00):
            continue
        s.add((p[2], int(p[3]), p[0].split("_")[-1]))
    return s

for name, bam in (("directional","dedup_directional.bam"),("unique","dedup_unique.bam"),("cluster","dedup_cluster.bam"),("coordinate_only","coordonly.bam")):
    ks = keys(bam)
    h1b = dvb = h1t = dvt = 0
    for g in groups.values():
        kept = [((c, st+1, u) in ks) for c, st, u, t in g]
        if g[0][3] == "hamming1":
            h1t += 1; h1b += all(kept)
        else:
            dvt += 1; dvb += all(kept)
    print("%-16s hamming1 both kept %d/%d   diverse both kept %d/%d" % (name, h1b, h1t, dvb, dvt))
EOF
