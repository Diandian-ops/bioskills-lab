import os, subprocess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
VCF = os.path.join(HERE, "chr22_slice.vcf.gz")

def count(expr):
    p = subprocess.run(["bcftools", "view", "-i", expr, VCF], capture_output=True)
    p2 = subprocess.run(["bcftools", "view", "-H"], input=p.stdout, capture_output=True)
    return p2.stdout.count(b"\n")

raw     = count("SVLEN >= 50")
abssvl  = count("ABS(SVLEN) >= 50")
endlen  = count("INFO/END-POS+1 >= 50")
endlen0 = count("INFO/END-POS >= 50")

# ---- fig1: the SVLEN representation trap ----
labels = ["SVLEN >= 50", "ABS(SVLEN) >= 50",
          "END - POS + 1 >= 50", "END - POS >= 50"]
vals = [raw, abssvl, endlen, endlen0]
colors = ["#b5482f", "#b5482f", "#2f7d72", "#2f7d72"]
fig, ax = plt.subplots(figsize=(7.8, 4.4))
bars = ax.bar(labels, vals, color=colors)
for b, v in zip(bars, vals):
    ax.text(b.get_x()+b.get_width()/2, v+0.12, str(v), ha="center", fontsize=11, color="#222222")
ax.set_ylabel("SV records retained")
ax.set_ylim(0, max(vals)*1.35 if max(vals) else 1)
ax.set_title("SV length filtering: SVLEN is absent in this callset\n"
              "chr22:17.0-17.2 Mb, 5 SV records (3 CNV + 2 DEL)")
plt.tight_layout()
fig.savefig(os.path.join(HERE, "fig1_svlen_trap.png"), dpi=130)
plt.close(fig)

# ---- fig2: SV span with breakpoint confidence (CIPOS) ----
# POS is the last unaffected anchor base; span length = END - POS
recs = [
    (17006397, 17022734, "CNV", None),
    (17026339, 17027748, "DEL", (-76, 50)),
    (17034911, 17076240, "CNV", None),
    (17071709, 17306387, "CNV", None),
    (17107196, 17108416, "DEL", (-366, 0)),
]
fig, ax = plt.subplots(figsize=(8.0, 4.2))
for i, (pos, end, st, cipos) in enumerate(recs):
    y = len(recs) - 1 - i
    ax.plot([pos, end], [y, y], linewidth=6,
            color="#b5482f" if st == "DEL" else "#2f7d72", solid_capstyle="butt")
    if cipos:
        ax.plot([pos+cipos[0], pos+cipos[1]], [y, y], linewidth=12,
                color="#b5482f", alpha=0.25, solid_capstyle="butt")
        ax.text(end, y + 0.22, "  imprecise  CIPOS=%d,%d" % cipos, fontsize=8.5, color="#b5482f")
    else:
        ax.text(end, y + 0.22, "  no CIPOS", fontsize=8.5, color="#666666")
ax.set_yticks(range(len(recs)))
ax.set_yticklabels(["%d\n%s" % (r[0], r[2]) for r in reversed(recs)], fontsize=9)
ax.set_ylim(-0.75, len(recs) + 0.35)
ax.set_xlabel("chr22 position (GRCh37)")
ax.set_title("SV spans and breakpoint confidence: shaded bands = CIPOS uncertainty")
ax.set_xlim(16990000, 17340000)
plt.tight_layout()
fig.savefig(os.path.join(HERE, "fig2_sv_spans_cipos.png"), dpi=130)
plt.close(fig)
print("fig1_svlen_trap.png + fig2_sv_spans_cipos.png written")
print("counts:", dict(zip(labels, vals)))
