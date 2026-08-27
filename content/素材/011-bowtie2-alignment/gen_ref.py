import random

random.seed(20260411)
bases = "ACGT"
contigs = []
for i in range(4):
    seq = "".join(random.choice(bases) for _ in range(3000))
    wrapped = "\n".join(seq[j:j + 70] for j in range(0, len(seq), 70))
    contigs.append(f">contig{i + 1} synthetic_reference_segment_{i + 1}\n" + wrapped)
with open("reference.fa", "w") as f:
    f.write("\n".join(contigs) + "\n")
total = sum(len(c) - c.count("\n") - c.split("\n", 1)[0].count(">") for c in contigs)
print("reference.fa written; contigs:", len(contigs), "total bp:", 4 * 3000)
