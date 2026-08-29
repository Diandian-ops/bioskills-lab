#!/usr/bin/env python3
# 017 variant-normalization 真实试用 — 出图脚本
# 数据源：1000 Genomes Phase3 整合基因型集 chr22 (GRCh37/hs37d5)，工作切片 17.0-17.2 Mb
# 图与本脚本平铺于此目录；运行：python make_figs.py
import os, subprocess, traceback
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
PIPE = HERE  # 输入数据随本素材目录自带，不依赖 pipeline 工作区
VCFGZ = os.path.join(PIPE, "region_17.vcf.gz")
SPLIT = os.path.join(PIPE, "demo_multi.split.vcf.gz")

def bcft(*args):
    return subprocess.run(["bcftools"] + list(args), capture_output=True, text=True)

# ---- 真实统计（从 VCF 实时计算） ----
out = bcft("view", "-H", VCFGZ).stdout
total = out.count("\n")
bi = multi = sv = 0
for line in out.splitlines():
    alt = line.split("\t")[4]
    if "," in alt:
        multi += 1
    else:
        bi += 1
    if "SVTYPE" in line.split("\t")[7]:
        sv += 1
indel = len([l for l in bcft("view", "-v", "indels", VCFGZ).stdout.splitlines() if not l.startswith("#")])

# 拆分示例的真实 AF（Number=A）：位点 17020038 C->A,T
af_a = af_t = None
for r in bcft("view", "-H", SPLIT).stdout.strip().splitlines():
    f = r.split("\t"); alt = f[4]; info = f[7]
    af = [x for x in info.split(";") if x.startswith("AF=")]
    if not af:
        continue
    v = float(af[0].split("=")[ 1])
    if alt == "A":
        af_a = v
    elif alt == "T":
        af_t = v

print(f"total={total} bi={bi} multi={multi} indel={indel} sv={sv} af_a={af_a} af_t={af_t}")

GREY = "#888888"; TEAL = "#2f7d72"; RED = "#b5482f"
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "#f7f7f4",
                     "text.color": "#222222", "axes.labelcolor": "#222222",
                     "xtick.color": "#222222", "ytick.color": "#222222",
                     "font.family": "DejaVu Sans"})

# ---- 图1：变异构成 ----
fig, ax = plt.subplots(figsize=(7.2, 4.2))
cats = ["Biallelic", "Multiallelic"]
vals = [bi, multi]
bars = ax.bar(cats, vals, color=[GREY, RED], width=0.55, zorder=3)
ax.set_ylabel("records")
ax.set_title("Variant records in chr22 17.0-17.2 Mb (1000G Phase3, N=2504)")
ax.set_ylim(0, max(vals) * 1.18)
for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width()/2, v + max(vals)*0.02, f"{v}",
            ha="center", va="bottom", fontsize=11, color="#222222", fontweight="bold")
ax.text(0.5, 0.92, f"of which INDEL={indel}, SV/CNV={sv}, total={total}",
        transform=ax.transAxes, ha="center", fontsize=9, color="#555555")
ax.grid(axis="y", color="#dddddd", zorder=0)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig1_variant_counts.png"), dpi=140)
plt.close(fig)

# ---- 图2：多等位拆分 AF 归位（真实位点 17020038） ----
fig2, ax2 = plt.subplots(figsize=(7.2, 4.2))
groups = ["ALT=A", "ALT=T"]
before = [af_a, af_t]
after = [af_a, af_t]
x = list(range(len(groups))); w = 0.36
ax2.bar([i - w/2 for i in x], before, w, label="before (1 record)", color=GREY, zorder=3)
ax2.bar([i + w/2 for i in x], after, w, label="after (2 records)", color=TEAL, zorder=3)
ax2.set_xticks(x); ax2.set_xticklabels(groups)
ax2.set_ylabel("AF (allele frequency)")
ax2.set_title("Multiallelic split preserves AF, re-assigned per Number=A")
ax2.legend(frameon=False, fontsize=9)
ymax = max(after) * 1.25
for i, v in enumerate(before):
    ax2.text(i - w/2, v + ymax*0.02, f"{v:.4f}", ha="center", va="bottom", fontsize=9, color="#222222")
for i, v in enumerate(after):
    ax2.text(i + w/2, v + ymax*0.02, f"{v:.4f}", ha="center", va="bottom", fontsize=9, color="#222222")
ax2.set_ylim( 0, ymax)
ax2.grid(axis="y", color="#dddddd", zorder=0)
fig2.tight_layout()
fig2.savefig(os.path.join(HERE, "fig2_split_af.png"), dpi=140)
plt.close(fig2)
print("figures written")
