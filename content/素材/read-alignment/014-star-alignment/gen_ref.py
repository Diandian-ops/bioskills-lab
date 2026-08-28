"""Generate a small synthetic reference + GTF for STAR real-trial reproduction.

chr1 carries a single intron (exon1 = 1-1000, exon2 = 1801-4000, intron = 1001-1800)
so a hand-built junction-spanning read can be placed across the gap with an N CIGAR.
chr2/3/4 are flat contigs each carrying one gene, to give GeneCounts multiple rows.
"""
import os, random

BASE = os.path.dirname(os.path.abspath(__file__))

# EXON / INTRON layout on chr1 (1-based, inclusive)
CHR1_LEN = 4000
EXON1_END = 1000          # exon1 = [1, 1000]
EXON2_START = 1801        # exon2 = [1801, 4000]; intron = [1001, 1800] (800 bp)
CONTIG_LEN = 3000         # chr2/3/4


def rseq(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def main(seed=20260827):
    chr1 = rseq(CHR1_LEN, seed)
    chr2 = rseq(CONTIG_LEN, seed + 1)
    chr3 = rseq(CONTIG_LEN, seed + 2)
    chr4 = rseq(CONTIG_LEN, seed + 3)

    fa = (f">chr1 gene_with_intron\n{chr1}\n"
          f">chr2 gene_flat2\n{chr2}\n"
          f">chr3 gene_flat3\n{chr3}\n"
          f">chr4 gene_flat4\n{chr4}\n")
    with open(f"{BASE}/reference.fa", "w") as f:
        f.write(fa)

    gtf = (
        'chr1\tstar_sim\tgene\t1\t4000\t.\t+\t.\tgene_id "GENE1"; transcript_id "TR1";\n'
        'chr1\tstar_sim\texon\t1\t1000\t.\t+\t.\tgene_id "GENE1"; transcript_id "TR1";\n'
        'chr1\tstar_sim\texon\t1801\t4000\t.\t+\t.\tgene_id "GENE1"; transcript_id "TR1";\n'
        'chr2\tstar_sim\tgene\t1\t3000\t.\t+\t.\tgene_id "GENE2"; transcript_id "TR2";\n'
        'chr2\tstar_sim\texon\t1\t3000\t.\t+\t.\tgene_id "GENE2"; transcript_id "TR2";\n'
        'chr3\tstar_sim\tgene\t1\t3000\t.\t+\t.\tgene_id "GENE3"; transcript_id "TR3";\n'
        'chr3\tstar_sim\texon\t1\t3000\t.\t+\t.\tgene_id "GENE3"; transcript_id "TR3";\n'
        'chr4\tstar_sim\tgene\t1\t3000\t.\t+\t.\tgene_id "GENE4"; transcript_id "TR4";\n'
        'chr4\tstar_sim\texon\t1\t3000\t.\t+\t.\tgene_id "GENE4"; transcript_id "TR4";\n'
    )
    with open(f"{BASE}/annotation.gtf", "w") as f:
        f.write(gtf)

    return {"chr1": chr1, "chr2": chr2, "chr3": chr3, "chr4": chr4}


if __name__ == "__main__":
    main()
