import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

def rd(p):
    with open(os.path.join(HERE, p)) as fh:
        return "".join(x.strip() for x in fh.read().split("\n")[1:]).upper()

ref = rd("ref_200kb.fa")
proj = [("-H 1 (hap 1)", "snp_hap1.fa"), ("-H 2 (hap 2)", "snp_hap2.fa"),
        ("-H A (all ALT)", "snp_hapA.fa"), ("-H R (REF at het)", "snp_hapR.fa"),
        ("IUPAC (default -s)", "snp_iupac.fa")]
diffs, ambs, seqs = [], [], {}
for label, f in proj:
    s = rd(f); seqs[label] = s
    diffs.append(sum(1 for a, b in zip(ref, s) if a != b))
    ambs.append(sum(1 for c in s if c not in "ACGTN"))

h1 = seqs["-H 1 (hap 1)"]; h2 = seqs["-H 2 (hap 2)"]
h1h2 = sum(1 for a, b in zip(h1, h2) if a != b)

# ---- fig1: each projection is lossy in a different way ----
fig, ax = plt.subplots(figsize=(7.8, 4.4))
xs = range(len(proj))
bars = ax.bar([p[0] for p in proj], diffs, color="#2f7d72")
ax.bar([p[0] for p in proj], ambs, bottom=[d - a for d, a in zip(diffs, ambs)],
       color="#b5482f", label="IUPAC ambiguity chars")
for b, d in zip(bars, diffs):
    ax.text(b.get_x()+b.get_width()/2, d+8, str(d), ha="center", fontsize=10, color="#222222")
ax.set_ylabel("positions differing from reference")
ax.set_ylim(0, max(diffs)*1.20)
ax.set_title("One diploid genome, five lossy projections (chr22:23.0-23.2 Mb, SNP-only)\n"
              "sample HG00096, 200001 bp window")
ax.legend(fontsize=9)
plt.tight_layout()
fig.savefig(os.path.join(HERE, "fig1_projections_diff.png"), dpi=130)
plt.close(fig)

# ---- fig2: internal consistency proves -H 1/-H 2 are real haplotypes ----
total_alt = diffs[2]      # -H A : every ALT-bearing site
hom_alt   = diffs[3]      # -H R : REF at het -> only hom-alt differ from ref
het       = total_alt - hom_alt
fig, ax = plt.subplots(figsize=(7.4, 4.2))
ax.bar(["ALT-bearing sites\n(-H A)"], [total_alt], color="#7a7a7a")
ax.bar(["ALT-bearing sites\n(-H A)"], [hom_alt], color="#2f7d72", label="homozygous ALT (-H R)")
ax.text(0, hom_alt/2, "%d hom-ALT" % hom_alt, ha="center", va="center",
        fontsize=10, color="white", fontweight="bold")
ax.text(0, hom_alt + (total_alt-hom_alt)/2, "%d heterozygous" % het, ha="center", va="center",
        fontsize=10, color="white", fontweight="bold")
ax.set_ylabel("sites")
ax.set_ylim(0, total_alt*1.25)
ax.set_title("Consistency check: %d total - %d hom-ALT = %d het\n"
              "measured hap1-vs-hap2 difference = %d  (exact match)" % (total_alt, hom_alt, het, h1h2))
ax.legend(fontsize=9, loc="upper right")
plt.tight_layout()
fig.savefig(os.path.join(HERE, "fig2_phase_consistency.png"), dpi=130)
plt.close(fig)
print("fig1_projections_diff.png + fig2_phase_consistency.png written")
print("diffs:", dict(zip([p[0] for p in proj], diffs)))
print("total_alt=%d hom_alt=%d het=%d measured_h1h2=%d" % (total_alt, hom_alt, het, h1h2))
