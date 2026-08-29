import os, subprocess, collections, statistics
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
VCF = os.path.join(HERE, "annotated.vcf.gz")

def q(fmt):
    return subprocess.run(["bcftools", "query", "-f", fmt, VCF],
                          capture_output=True, text=True).stdout.splitlines()

sig = q('%INFO/CLNSIG\n'); rev = q('%INFO/CLNREVSTAT\n')
af  = q('%INFO/AF\n')
pops = q('%INFO/EAS_AF\t%INFO/EUR_AF\t%INFO/AFR_AF\t%INFO/AMR_AF\t%INFO/SAS_AF\n')

def stars(r):
    if "practice_guideline" in r: return 4
    if "reviewed_by_expert_panel" in r: return 3
    if "multiple_submitters" in r and "no_conflicts" in r: return 2
    if "criteria_provided" in r: return 1
    return 0

def f(x):
    try: return float(x)
    except Exception: return None

# ---- fig1: ClinVar review status = star rating ----
rows = [(s, r) for s, r in zip(sig, rev) if s != "." and r != "."]
cnt = collections.Counter(stars(r) for _, r in rows)
labels = ["4-star\npractice_guideline", "3-star\nexpert panel",
          "2-star\nmulti-submitter", "1-star\n(lead only)"]
vals = [cnt.get(4, 0), cnt.get(3, 0), cnt.get(2, 0), cnt.get(1, 0) + cnt.get(0, 0)]
colors = ["#2f7d72", "#2f7d72", "#7a7a7a", "#b5482f"]
fig, ax = plt.subplots(figsize=(7.6, 4.4))
bars = ax.bar(labels, vals, color=colors)
for b, v in zip(bars, vals):
    ax.text(b.get_x()+b.get_width()/2, v+3, str(v), ha="center", fontsize=10, color="#222222")
ax.set_ylabel("annotated variant sites")
ax.set_ylim(0, max(vals)*1.22)
ax.set_title("ClinVar assertions are leads, not evidence\n"
              "chr22:21.2-21.4 Mb | %d of %d sites carried a ClinVar record" % (len(rows), len(sig)))
plt.tight_layout()
fig.savefig(os.path.join(HERE, "fig1_clinvar_stars.png"), dpi=130)
plt.close(fig)

# ---- fig2: population-max AF vs global AF ----
ratios = []
for a, p in zip(af, pops):
    g = f(a)
    vals = [f(v) for v in p.split("\t")]
    vals = [v for v in vals if v is not None]
    if g is None or not vals or g <= 0:
        continue
    ratios.append(max(vals) / g)
nb = 30
lo, hi = min(ratios), min(max(ratios), 8.0)
w = (hi - lo) / nb
hist = [0]*nb
for r in ratios:
    b = min(nb-1, int((min(r, hi) - lo)/w))
    hist[b] += 1
centers = [lo + (i+0.5)*w for i in range(nb)]
fig, ax = plt.subplots(figsize=(7.6, 4.4))
ax.bar(centers, hist, width=w*0.9, color="#b5482f")
med = statistics.median(ratios)
ax.axvline(med, color="#222222", linestyle="--", linewidth=1.2)
ax.text(med + 0.12, max(hist)*0.82, "median = %.2f x" % med, fontsize=10, color="#222222")
ax.axvline(2.0, color="#2f7d72", linestyle=":", linewidth=1.2)
ax.text(2.12, max(hist)*0.60, "2x line", fontsize=9, color="#2f7d72")
ax.set_xlabel("population-max AF / global AF")
ax.set_ylabel("variant sites")
ax.set_title("Global AF dilutes ancestry-specific frequency\n"
              "%.1f%% of %d sites have popmax >= 2x the global AF"
              % (100.0*sum(1 for r in ratios if r >= 2)/len(ratios), len(ratios)))
plt.tight_layout()
fig.savefig(os.path.join(HERE, "fig2_popmax_vs_global_af.png"), dpi=130)
plt.close(fig)
print("fig1_clinvar_stars.png + fig2_popmax_vs_global_af.png written")
print("annotated sites:", len(rows), "| ratios n:", len(ratios), "| median ratio: %.2f" % med)
