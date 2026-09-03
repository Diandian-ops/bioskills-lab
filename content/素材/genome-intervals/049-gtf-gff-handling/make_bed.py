#!/usr/bin/env python3
# make_bed.py -- GTF -> BED per skill: subtract 1 from start ONLY, end unchanged.
# Also emits the over-corrected variant (start-1 AND end-1) as a negative control,
# and the seqid-remapped GFF3 (1 -> chr1) used to fix the namespace mismatch.
import json
import os

OUT = os.path.dirname(os.path.abspath(__file__))
truth = json.load(open(os.path.join(OUT, "truth.json")))

exons_rows, over_rows, genes_rows = [], [], []
for txid, t in sorted(truth["tx"].items()):
    exons = t["exons"]
    ordered = exons if t["strand"] == "+" else sorted(exons, reverse=True)
    for s, e in ordered:  # transcriptional order so concatenation = spliced sense seq
        exons_rows.append((truth["chr"], s - 1, e, txid, ".", t["strand"]))
        over_rows.append((truth["chr"], s - 1, e - 1, txid, ".", t["strand"]))
    gs = min(x[0] for x in exons)
    ge = max(x[1] for x in exons)
    genes_rows.append((truth["chr"], gs - 1, ge, txid, ".", t["strand"]))

for name, rows in [("exons.bed", exons_rows), ("exons_overcorrected.bed", over_rows),
                   ("genes.bed", genes_rows)]:
    with open(os.path.join(OUT, name), "w") as f:
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")

# seqid remap: the only fix a namespace mismatch needs
with open(os.path.join(OUT, "annotation.gff3")) as f:
    lines = f.readlines()
with open(os.path.join(OUT, "annotation_chr1.gff3"), "w") as f:
    for ln in lines:
        if ln.startswith("#"):
            f.write(ln)
        else:
            p = ln.rstrip("\n").split("\t")
            p[0] = truth["chr"]
            f.write("\t".join(p) + "\n")

print("bed rows: exons=%d overcorrected=%d genes=%d; remapped gff3 written"
      % (len(exons_rows), len(over_rows), len(genes_rows)))
