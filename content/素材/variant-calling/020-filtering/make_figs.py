import os, subprocess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "chr22_slice.vcf.gz")

def bcft(*args):
    return subprocess.run(["bcftools"] + list(args), capture_output=True, text=True).stdout

def titv(path):
    out = bcft("stats", path)
    for line in out.splitlines():
        if line.startswith("TSTV"):
            f = line.split("\t")
            # TSTV [1]id [2]ts [3]tv [4]ts/tv [5]ts/tv (1st ALT)
            return float(f[4]) if len(f) > 4 else None
    return None

def nrec(path):
    return int(bcft("view", "-H", path).count("\n"))

THRESH = [0.0, 0.0001, 0.001, 0.01, 0.05]
counts, ratios = [], []
for t in THRESH:
    if t == 0.0:
        src = RAW
        counts.append(nrec(src)); ratios.append(titv(src))
    else:
        p = os.path.join(HERE, "_tmp_%s.vcf.gz" % str(t).replace(".", ""))
        bcft("view", "-i", "INFO/AF >= %s" % t, RAW, "-Oz", "-o", p)
        bcft("index", "-f", p)
        counts.append(nrec(p)); ratios.append(titv(p))
        os.remove(p); os.remove(p + ".csi")

labels = ["all", "AF>=1e-4", "AF>=1e-3", "AF>=1e-2", "AF>=5e-2"]

# ---- fig1: Ti/Tv vs AF threshold ----
fig, ax = plt.subplots(figsize=(7.4, 4.3))
bars = ax.bar(labels, ratios, color="#2f7d72")
ax.axhspan(2.0, 2.1, color="#b5482f", alpha=0.15, zorder=0)
ax.axhspan(3.0, 3.3, color="#7a7a7a", alpha=0.12, zorder=0)
for b, v in zip(bars, ratios):
    ax.text(b.get_x()+b.get_width()/2, v+0.03, "%.2f" % v, ha="center", fontsize=10, color="#222222")
ax.set_ylabel("Ti/Tv")
ax.set_ylim(0, 3.6)
ax.set_title("Filter validation: Ti/Tv rises as rare (singleton-rich) variants are removed\n"
              "chr22:17.0-17.2 Mb, 2504 samples, 5431 sites | shaded: WGS 2.0-2.1 / WES 3.0-3.3")
plt.tight_layout()
fig.savefig(os.path.join(HERE, "fig1_titv_after_filter.png"), dpi=130)
plt.close(fig)

# ---- fig2: records retained vs dropped ----
dropped = [counts[0]-c for c in counts]
fig, ax = plt.subplots(figsize=(7.4, 4.3))
ax.bar(labels, counts, color="#2f7d72", label="retained")
ax.bar(labels, dropped, bottom=counts, color="#b5482f", alpha=0.75, label="removed")
for i, c in enumerate(counts):
    ax.text(i, c+60, str(c), ha="center", fontsize=9.5, color="#222222")
ax.set_ylabel("variant records")
ax.set_ylim(0, counts[0]*1.18)
ax.set_title("bcftools filter on AF: retained vs removed")
ax.legend(fontsize=9, loc="upper right")
plt.tight_layout()
fig.savefig(os.path.join(HERE, "fig2_retained_vs_removed.png"), dpi=130)
plt.close(fig)
print("fig1_titv_after_filter.png + fig2_retained_vs_removed.png written")
print("Ti/Tv:", list(zip(labels, ratios)))
print("counts:", list(zip(labels, counts)))
