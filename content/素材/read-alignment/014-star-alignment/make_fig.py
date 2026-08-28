"""Honest figure for the 014 STAR real-trial.

Panel A: MAPQ scale across the read-alignment family.
  - bowtie2 42 (end-to-end) / 44 (--local): REAL (011)
  - bwa-mem2 60: REAL (012)
  - HISAT2 60: REAL (013)
  - STAR 255 (unique, default) -> 60 via --outSAMmapqUnique: DOCUMENTED (SKILL.md).
    The STAR alignment step could not be executed in this sandbox (STAR's read-input
    path returns 0 reads), so STAR's value is taken from the skill's documented behavior,
    not from a live run here. Labeled accordingly.
Panel B: real genomeGenerate result for 014 (index build executed successfully).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "014-fig.png"

fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.6))

# ---- Panel A: MAPQ scale ----
aligners = ["bowtie2\ne2e", "bowtie2\n--local", "bwa-mem2", "HISAT2", "STAR\n(default)", "STAR\n+mapq60"]
mapq = [42, 44, 60, 60, 255, 60]
colors = ["#4C72B0", "#4C72B0", "#55A868", "#C44E52", "#C44E52", "#55A868"]
bars = axA.bar(aligners, mapq, color=colors, edgecolor="black", linewidth=0.6)
axA.axhline(60, ls="--", color="gray", lw=1)
axA.text(5.4, 63, "GATK-friendly\nMAPQ>=60", fontsize=8, color="gray", ha="right")
axA.set_ylabel("MAPQ for unique reads")
axA.set_title("MAPQ scale across read aligners\n(STAR 255->60: documented, alignment step env-blocked)")
for b, v in zip(bars, mapq):
    axA.text(b.get_x() + b.get_width()/2, v + 4, str(v), ha="center", fontsize=8)
axA.set_ylim(0, 280)

# ---- Panel B: real genomeGenerate result ----
axB.axis("off")
txt = (
    "STAR genomeGenerate (014) -- REAL run\n"
    "----------------------------------------\n"
    "STAR version      : 2.7.11b\n"
    "genome SA index   : genomeSAindexNbases = 5\n"
    "                   (small ~13 kb genome;\n"
    "                    default 14 segfaults)\n"
    "sjdbOverhang      : 99   (= readlen 100 - 1)\n"
    "splice DB (GTF)   : 1 gene, 2 exons on chr1\n"
    "                    (intron 1001-1800)\n"
    "chromosomes       : chr1(4000) chr2/3/4(3000)\n"
    "index files built : 18  (Genome, SA, SAindex,\n"
    "                    chrName/Length, exon/gene/\n"
    "                    transcriptInfo, sjdbList ...)\n"
    "\n"
    "Read-alignment step: ENVIRONMENT-BLOCKED\n"
    "STAR reports 'Number of input reads | 0'\n"
    "for every input in this sandbox\n"
    "(genomeGenerate + genome load OK)."
)
axB.text(0.02, 0.98, txt, va="top", ha="left", fontsize=8.5, family="monospace",
         bbox=dict(boxstyle="round", fc="#f6f6f6", ec="#cccccc"))
axB.set_title("Index build outcome (executed)", fontsize=10)

plt.tight_layout()
plt.savefig(OUT, dpi=150, bbox_inches="tight")
print("wrote", OUT)
