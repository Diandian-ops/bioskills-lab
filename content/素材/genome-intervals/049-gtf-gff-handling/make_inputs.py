#!/usr/bin/env python3
# make_inputs.py -- self-contained mini annotation set for 049 gtf-gff-handling
# Fixed seed 49. One 2 Mb chromosome (chr1), 8 genes / 12 transcripts.
# Deliberate convention splits between the GTF and GFF3 versions:
#   - GTF: seqid "chr1", GENCODE-style attributes (gene_type), versioned IDs
#     (G0001.1 / T0001.1), CDS EXCLUDES the stop codon.
#   - GFF3: seqid "1" (accession-style namespace), Ensembl/RefSeq-style
#     attributes (gene_biotype), unversioned IDs, CDS INCLUDES the stop codon.
# Both carry correctly computed phase on every CDS segment.
import json
import os
import random

SEED = 49
CHR = "chr1"
GFF3_SEQID = "1"
CHR_LEN = 2_000_000
SRC = "test"
OUT = os.path.dirname(os.path.abspath(__file__))

rng = random.Random(SEED)


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


# ---------- layout ----------
genes = []
for i in range(8):
    gstart = 50_000 + i * 220_000 + rng.randint(0, 20_000)
    strand = "+" if i % 3 == 0 else "-"
    genes.append({"id": "G%04d" % (i + 1), "start": gstart, "strand": strand,
                  "n_tx": 2 if i % 2 == 0 else 1})

STOP = {"TAA", "TAG", "TGA"}


def make_exons(strand, region):
    n = rng.randint(2, 5)
    exons = []
    pos = region[0]
    for k in range(n):
        elen = rng.randint(150, 400)
        exons.append([pos, pos + elen - 1])
        pos += elen + rng.randint(500, 2000)
    return exons


def tx_len(exons):
    return sum(e[1] - e[0] + 1 for e in exons)


def tx_coords(exons):
    """map transcript coordinate (1-based) -> genomic position"""
    m = []
    for s, e in exons:
        m.extend(range(s, e + 1))
    return m


coding_design = {}  # genomic pos (1-based) -> base, for the primary CDS window
truth_tx = {}

for g in genes:
    region = (g["start"], g["start"] + 14_000)
    exonsA = make_exons(g["strand"], region)
    if exonsA[-1][1] > CHR_LEN:
        raise SystemExit("layout overflow")
    txs = [{"id": "T%04d" % (g["start"] % 10 + 1), "exons": exonsA}]
    # a fixed transcript numbering: per-gene index
    txs[0]["id"] = "T%s01" % g["id"][1:]
    if g["n_tx"] == 2 and len(exonsA) >= 3:
        skip = rng.randint(1, len(exonsA) - 2)  # cassette exon
        exonsB = [e for k, e in enumerate(exonsA) if k != skip]
        txs.append({"id": "T%s02" % g["id"][1:], "exons": exonsB})

    # coding window on primary transcript A (transcript coordinates, 1-based)
    u = rng.randint(50, 120)                      # 5' UTR
    L = rng.choice([150, 198, 249, 303, 360, 447])  # incl. stop, multiple of 3
    assert u + L + 30 <= tx_len(exonsA)
    cdsA_tx = (u + 1, u + L)                      # window INCLUDING stop

    for t in txs:
        coords = tx_coords(t["exons"])
        idx = {p: k for k, p in enumerate(coords)}
        # CDS genomic positions = window positions present in this transcript
        win = [coords[p - 1] for p in range(cdsA_tx[0], cdsA_tx[1] + 1)
               if p - 1 < len(coords)]
        segs = []
        for p in win:
            if segs and segs[-1][1] == p - 1:
                segs[-1][1] = p
            else:
                segs.append([p, p])
        t["cds_full"] = segs        # includes stop codon (GFF3 convention)
        # GTF convention: CDS excludes the transcriptionally-last 3 nt (stop).
        # transcriptional end = rightmost genomic end on '+', leftmost on '-'.
        t["cds_nostop"] = []
        if g["strand"] == "+":
            first, last = segs[0], segs[-1]
            rem = last[1] - last[0] + 1
            for s in segs[:-1]:
                t["cds_nostop"].append(list(s))
            if rem > 3:
                t["cds_nostop"].append([last[0], last[1] - 3])
        else:
            first, last = segs[0], segs[-1]
            rem = first[1] - first[0] + 1
            if rem > 3:
                t["cds_nostop"].append([first[0] + 3, first[1]])
            for s in segs[1:]:
                t["cds_nostop"].append(list(s))
        t["strand"] = g["strand"]

        # design coding sequence once, on transcript A's full window
        if t is txs[0]:
            cdsL_full = sum(e - s + 1 for s, e in segs)
            codons = ["ATG"]
            for _ in range(cdsL_full // 3 - 2):
                while True:
                    c = "".join(rng.choice("ACGT") for _ in range(3))
                    if c not in STOP:
                        break
                codons.append(c)
            codons.append("TAA")
            seq = "".join(codons)
            if g["strand"] == "-":
                seq = rc(seq)
            for p, b in zip(win, seq):
                coding_design[p] = b

        truth_tx[t["id"]] = {
            "gene": g["id"], "strand": g["strand"],
            "exons": [list(e) for e in t["exons"]],
            "cds_full": t["cds_full"], "cds_nostop": t["cds_nostop"],
        }
    g["txs"] = txs

# ---------- genome ----------
seq = [rng.choice("ACGT") for _ in range(CHR_LEN)]
for p, b in coding_design.items():
    seq[p - 1] = b
seq = "".join(seq)

with open(os.path.join(OUT, "genome.fa"), "w") as f:
    f.write(">%s\n" % CHR)
    for i in range(0, CHR_LEN, 60):
        f.write(seq[i:i + 60] + "\n")


# ---------- phases ----------
def assign_phase(segs, strand):
    """phase = bases to trim from the segment's 5' (transcriptional) end."""
    ph = []
    cum = 0
    ordered = segs if strand == "+" else sorted(segs, reverse=True)
    for s, e in ordered:
        ph.append((3 - (cum % 3)) % 3)
        cum += e - s + 1
    return ph


def fmt_gtf_attrs(d):
    return " ".join('%s "%s";' % (k, v) for k, v in d)


def fmt_gff3_attrs(d):
    return ";".join("%s=%s" % (k, v) for k, v in d)


# ---------- GTF ----------
gtf_rows = []
for g in genes:
    gstart = min(e[0] for t in g["txs"] for e in t["exons"])
    gend = max(e[1] for t in g["txs"] for e in t["exons"])
    gtf_rows.append((gstart, [CHR, SRC, "gene", str(gstart), str(gend), ".", g["strand"],
                              ".", fmt_gtf_attrs([("gene_id", g["id"] + ".1"),
                                                  ("gene_type", "protein_coding"),
                                                  ("gene_name", "gene_" + g["id"])])]))
    for t in g["txs"]:
        ts = min(e[0] for e in t["exons"])
        te = max(e[1] for e in t["exons"])
        base = [("gene_id", g["id"] + ".1"), ("transcript_id", t["id"] + ".1"),
                ("gene_type", "protein_coding")]
        gtf_rows.append((ts, [CHR, SRC, "transcript", str(ts), str(te), ".", g["strand"],
                              ".", fmt_gtf_attrs(base)]))
        for k, e in enumerate(sorted(t["exons"])):
            gtf_rows.append((e[0], [CHR, SRC, "exon", str(e[0]), str(e[1]), ".", g["strand"],
                                    ".", fmt_gtf_attrs(base + [("exon_number", str(k + 1))])]))
        phs = assign_phase(t["cds_nostop"], g["strand"])
        ordered = t["cds_nostop"] if g["strand"] == "+" else sorted(t["cds_nostop"], reverse=True)
        for (s, e), ph in zip(ordered, phs):
            gtf_rows.append((s, [CHR, SRC, "CDS", str(s), str(e), ".", g["strand"], str(ph),
                                 fmt_gtf_attrs(base + [("exon_number", "1")])]))

gtf_rows.sort(key=lambda r: r[0])
with open(os.path.join(OUT, "annotation.gtf"), "w") as f:
    for _, row in gtf_rows:
        f.write("\t".join(row) + "\n")

# ---------- GFF3 ----------
with open(os.path.join(OUT, "annotation.gff3"), "w") as f:
    f.write("##gff-version 3\n")
    f.write("##sequence-region %s 1 %d\n" % (GFF3_SEQID, CHR_LEN))
    for g in genes:
        gstart = min(e[0] for t in g["txs"] for e in t["exons"])
        gend = max(e[1] for t in g["txs"] for e in t["exons"])
        f.write("\t".join([GFF3_SEQID, SRC, "gene", str(gstart), str(gend), ".",
                           g["strand"], ".",
                           fmt_gff3_attrs([("ID", "gene:" + g["id"]),
                                           ("Name", "gene_" + g["id"]),
                                           ("gene_biotype", "protein_coding")])]) + "\n")
        for t in g["txs"]:
            ts = min(e[0] for e in t["exons"])
            te = max(e[1] for e in t["exons"])
            f.write("\t".join([GFF3_SEQID, SRC, "mRNA", str(ts), str(te), ".",
                               g["strand"], ".",
                               fmt_gff3_attrs([("ID", "transcript:" + t["id"]),
                                               ("Parent", "gene:" + g["id"]),
                                               ("gene_biotype", "protein_coding")])]) + "\n")
            for k, e in enumerate(sorted(t["exons"])):
                f.write("\t".join([GFF3_SEQID, SRC, "exon", str(e[0]), str(e[1]), ".",
                                   g["strand"], ".",
                                   fmt_gff3_attrs([("ID", "exon:%s:%d" % (t["id"], k + 1)),
                                                   ("Parent", "transcript:" + t["id"]),
                                                   ("rank", str(k + 1))])]) + "\n")
            phs = assign_phase(t["cds_full"], g["strand"])
            ordered = t["cds_full"] if g["strand"] == "+" else sorted(t["cds_full"], reverse=True)
            for (s, e), ph in zip(ordered, phs):
                f.write("\t".join([GFF3_SEQID, SRC, "CDS", str(s), str(e), ".",
                                   g["strand"], str(ph),
                                   fmt_gff3_attrs([("ID", "cds:%s" % t["id"]),
                                                   ("Parent", "transcript:" + t["id"])])]) + "\n")

# ---------- truth ----------
truth = {"seed": SEED, "chr": CHR, "chr_len": CHR_LEN, "gff3_seqid": GFF3_SEQID,
         "n_genes": len(genes), "n_tx": len(truth_tx),
         "tx": truth_tx}
with open(os.path.join(OUT, "truth.json"), "w") as f:
    json.dump(truth, f, indent=1)

n_exons = sum(len(t["exons"]) for t in truth_tx.values())
n_cds_gtf = sum(len(t["cds_nostop"]) for t in truth_tx.values())
n_cds_gff = sum(len(t["cds_full"]) for t in truth_tx.values())
print("seed=%d chr_len=%d genes=%d tx=%d exons=%d CDS_segs_gtf=%d CDS_segs_gff3=%d designed_bases=%d"
      % (SEED, CHR_LEN, len(genes), len(truth_tx), n_exons, n_cds_gtf, n_cds_gff,
         len(coding_design)))
