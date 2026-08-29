import os, collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CSQ = os.path.join(HERE, "csq_blocks.tsv")

def parse(path):
    cons = collections.Counter()
    per_pos = collections.Counter()
    genes = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            key = f[3] + ":" + f[4]
            blocks = [b for b in f[5].split(",") if b]
            per_pos[key] += len(blocks)
            for b in blocks:
                parts = b.split("|")
                cons[parts[0]] += 1
                if len(parts) > 1 and parts[1]:
                    genes.add(parts[1])
    return cons, per_pos, genes

cons, per_pos, genes = parse(CSQ)

# ---- fig1: consequence type distribution ----
order = sorted(cons.items(), key=lambda kv: kv[1])
names = [k for k, _ in order]
vals = [v for _, v in order]
fig, ax = plt.subplots(figsize=(7.6, 4.6))
bars = ax.barh(names, vals, color="#2f7d72")
for b, v in zip(bars, vals):
    ax.text(v + 1.0, b.get_y() + b.get_height()/2, str(v),
            va="center", ha="left", fontsize=9.5, color="#222222")
ax.set_xlabel("consequence blocks (bcftools csq, per sample-haplotype)")
ax.set_xlim(0, max(vals)*1.20)
ax.set_title("Variant consequence distribution, chr22:23.0-23.2 Mb\n"
              "5 samples x 2 haplotypes | Ensembl GRCh37.87 transcripts | %d genes hit" % len(genes))
plt.tight_layout()
fig.savefig(os.path.join(HERE, "fig1_consequence_types.png"), dpi=130)
plt.close(fig)

# ---- fig2: transcripts hit per variant (why transcript selection matters) ----
hist = collections.Counter(per_pos.values())
ks = sorted(hist)
counts = [hist[k] for k in ks]
total = sum(counts)
multi = sum(c for k, c in hist.items() if k >= 2)
fig, ax = plt.subplots(figsize=(7.6, 4.4))
bars = ax.bar([str(k) for k in ks], counts, color="#b5482f")
for b, v in zip(bars, counts):
    ax.text(b.get_x()+b.get_width()/2, v+0.4, str(v), ha="center", fontsize=9.5, color="#222222")
ax.set_xlabel("transcript blocks hit per variant position")
ax.set_ylabel("variant positions")
ax.set_ylim(0, max(counts)*1.22)
ax.set_title("One variant, many transcripts: %d of %d annotated positions (%.1f%%) hit >=2 blocks"
             % (multi, total, 100.0*multi/total))
plt.tight_layout()
fig.savefig(os.path.join(HERE, "fig2_multi_transcript.png"), dpi=130)
plt.close(fig)
print("fig1_consequence_types.png + fig2_multi_transcript.png written")
print("consequences:", dict(cons))
print("genes:", len(genes), "| positions:", total, "| multi-transcript:", multi)
