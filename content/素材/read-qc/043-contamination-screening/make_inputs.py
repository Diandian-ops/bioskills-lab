#!/usr/bin/env python3
"""043 contamination-screening: generate paired-end reads with exactly designed
contamination fractions from three real reference genomes."""
import gzip
import os
import random

BASE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(BASE, "refs")
OUT = os.path.join(BASE, "reads")
os.makedirs(OUT, exist_ok=True)

random.seed(43)

READ_LEN = 150
INSERT_MEAN, INSERT_SD = 300, 40
ERR_RATE = 0.01          # per-base substitution error
Q_GOOD = (28, 38)        # quality range for correct bases
Q_ERR = (8, 18)          # quality range for erroneous bases
TOTAL = 50000            # read pairs per sample

GENOMES = {"Ecoli": "ecoli.fna", "PhiX": "phix.fna", "Lambda": "lambda.fna"}

MIXES = {
    "S1": {"Ecoli": 0.90, "PhiX": 0.05, "Lambda": 0.05},
    "S2": {"Ecoli": 0.75, "PhiX": 0.15, "Lambda": 0.10},
    "S3": {"Ecoli": 1.00, "PhiX": 0.00, "Lambda": 0.00},
}

COMP = str.maketrans("ACGTN", "TGCAN")


def load_fasta(path):
    parts = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                parts.append(line.strip())
    return "".join(parts).upper()


def revcomp(s):
    return s.translate(COMP)[::-1]


def make_pair(seq):
    frag_len = min(len(seq), max(READ_LEN * 2, int(random.gauss(INSERT_MEAN, INSERT_SD))))
    start = random.randrange(0, len(seq) - frag_len + 1)
    frag = seq[start:start + frag_len]

    def build(r):
        bases, quals = [], []
        for b in r:
            if random.random() < ERR_RATE:
                bases.append(random.choice([c for c in "ACGT" if c != b]))
                quals.append(random.randint(*Q_ERR))
            else:
                bases.append(b)
                quals.append(random.randint(*Q_GOOD))
        return "".join(bases), "".join(chr(33 + q) for q in quals)

    r1, q1 = build(frag[:READ_LEN])
    r2, q2 = build(revcomp(frag[-READ_LEN:]))
    return r1, q1, r2, q2


def write_fastq(path, records):
    with gzip.open(path, "wt", compresslevel=6) as fh:
        for i, (r, q) in enumerate(records, 1):
            fh.write(f"@{os.path.basename(path).rsplit('.', 2)[0]}:{i}\n{r}\n+\n{q}\n")


def main():
    genomes = {k: load_fasta(os.path.join(REF, v)) for k, v in GENOMES.items()}
    for k, v in genomes.items():
        print(f"loaded {k}: {len(v)} bp")

    design_rows = ["sample\ttaxon\tpairs\tdesigned_pct"]
    for sample, mix in MIXES.items():
        buckets = []
        for taxon, frac in mix.items():
            n = round(frac * TOTAL)
            if n == 0:
                continue
            reads1, reads2 = [], []
            for _ in range(n):
                r1, q1, r2, q2 = make_pair(genomes[taxon])
                reads1.append((r1, q1))
                reads2.append((r2, q2))
            buckets.append((taxon, reads1, reads2))
            design_rows.append(f"{sample}\t{taxon}\t{n}\t{frac * 100:.1f}")
            print(f"{sample} {taxon}: {n} pairs ({frac * 100:.1f}%)")

        # interleave contaminants randomly among host reads
        order = []
        for taxon, reads1, reads2 in buckets:
            order.extend((taxon, i) for i in range(len(reads1)))
        random.shuffle(order)

        idx = {t: 0 for t in genomes}
        bucket_map = {b[0]: b for b in buckets}
        w1, w2 = [], []
        for taxon, i in order:
            w1.append(bucket_map[taxon][1][i])
            w2.append(bucket_map[taxon][2][i])

        write_fastq(os.path.join(OUT, f"{sample}_1.fastq.gz"), w1)
        write_fastq(os.path.join(OUT, f"{sample}_2.fastq.gz"), w2)
        print(f"{sample}: wrote {len(w1)} pairs")

    with open(os.path.join(BASE, "design.tsv"), "w") as fh:
        fh.write("\n".join(design_rows) + "\n")
    print("design.tsv written")


if __name__ == "__main__":
    main()
