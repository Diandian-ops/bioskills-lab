import os, gzip
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = HERE  # 输入数据随本素材目录自带，不依赖 pipeline 工作区

# ---- fig1: isec partition sizes (raw vs norm) ----
def count(path):
    n = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            if not line.startswith("#"):
                n += 1
    return n

p0 = os.path.join(WORK, "isec", "0000.vcf.gz")   # private to raw
p1 = os.path.join(WORK, "isec", "0001.vcf.gz")   # private to norm
p2 = os.path.join(WORK, "isec", "0002.vcf.gz")   # shared
n0, n1, n2 = count(p0), count(p1), count(p2)

fig, ax = plt.subplots(figsize=(7, 4.2))
labels = ["raw-only\n(private to raw)", "norm-only\n(private to norm)", "shared\n(both)"]
vals = [n0, n1, n2]
colors = ["#2f7d72", "#b5482f", "#7a7a7a"]
bars = ax.bar(labels, vals, color=colors)
ax.set_ylabel("variant records")
ax.set_title("bcftools isec: raw vs normalized chr22 slice (chr22:17.0-17.2 Mb)\n"
              "multi-allelic split adds 80 norm-only records (same sites, different ALT tuples)")
for b, v in zip(bars, vals):
    ax.text(b.get_x()+b.get_width()/2, v+40, str(v), ha="center", fontsize=10, color="#222222")
ax.set_ylim(0, max(vals)*1.15)
plt.tight_layout()
fig.savefig(os.path.join(HERE, "fig1_isec_partitions.png"), dpi=130)
plt.close(fig)

# ---- fig2: stale AC/AN before/after fill-tags (first 12 records subset to 5 samples) ----
def an_values(path, n=12):
    out = []
    with gzip.open(path, "rt") as f:
        cnt = 0
        for line in f:
            if line.startswith("#"):
                continue
            info = line.split("\t")[7]
            d = dict(kv.split("=", 1) for kv in info.split(";") if "=" in kv)
            out.append(int(d.get("AN", 0)))
            cnt += 1
            if cnt >= n:
                break
    return out

stale = os.path.join(WORK, "subset5_stale.vcf.gz")
fixed = os.path.join(WORK, "subset5_fixed2.vcf.gz")
sa = an_values(stale); fa = an_values(fixed)
xs = list(range(1, len(sa)+1))

fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(xs, sa, marker="o", color="#b5482f", label="before (stale AN from 2504 samples)", linewidth=1.5)
ax.plot(xs, fa, marker="s", color="#2f7d72", label="after fill-tags (5 samples)", linewidth=1.5)
ax.set_xlabel("record index (subset of 5 samples)")
ax.set_ylabel("AN (allele count denominator)")
ax.set_title("Stale AC/AN after subsetting: fix with bcftools +fill-tags")
ax.legend(fontsize=9)
ax.set_ylim(0, max(sa)*1.1)
plt.tight_layout()
fig.savefig(os.path.join(HERE, "fig2_stale_af_an.png"), dpi=130)
plt.close(fig)
print("fig1_isec_partitions.png + fig2_stale_af_an.png written")
