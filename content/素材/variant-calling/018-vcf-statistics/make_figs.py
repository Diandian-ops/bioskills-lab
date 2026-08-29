"""018 vcf-statistics — 真实数据图（1000 Genomes Phase3 chr22 切片，2504 样本）。

数据：chr22_slice.vcf.gz（chr22:17.0-17.2 Mb, GRCh37/hs37d5, 2504 样本, 5431 位点）
工具：bcftools 1.24（stats / stats -s -）
用法：python make_figs.py  （与脚本同目录放 chr22_slice.vcf.gz）
"""
import os
import subprocess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
VCF = os.path.join(HERE, "chr22_slice.vcf.gz")

C_RED = "#b5482f"      # 砖红
C_TEAL = "#2f7d72"     # 青瓷绿
C_GREY = "#555555"
C_LIGHT = "#f5f5f0"

def run(args):
    r = subprocess.run(["bcftools", "stats"] + args + [VCF],
                      capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr)
    return r.stdout

# ---------- 1. 总体统计 ----------
stats = run([])
trans = transv = None
for ln in stats.splitlines():
    if ln.startswith("TSTV"):
        f = ln.split("\t")
        # TSTV [1]id [2]nTransitions [3]nTransversions [4]TsTv ...
        trans = int(f[2]); transv = int(f[3]); tstv = float(f[4])
snp_n = indel_n = multi_n = other_n = None
for ln in stats.splitlines():
    if ln.startswith("SN"):
        f = ln.split("\t")
        if "number of SNPs:" in f[2]: snp_n = int(f[3])
        elif "number of indels:" in f[2]: indel_n = int(f[3])
        elif "number of multiallelic sites:" in f[2]: multi_n = int(f[3])
        elif "number of others:" in f[2]: other_n = int(f[3])

sing_ts = sing_tv = None
for ln in stats.splitlines():
    if ln.startswith("SiS"):
        f = ln.split("\t")
        # SiS [1]id [2]nSingletons [3]nTransitions [4]nTransversions ...
        sing_ts = int(f[3]); sing_tv = int(f[4])

# ---------- 2. 每样本 het/hom-alt 比值 ----------
psc = run(["-s", "-"])
ratios = []
for ln in psc.splitlines():
    if not ln.startswith("PSC"):
        continue
    f = ln.split("\t")
    nalt = int(f[4]); het = int(f[5])
    if nalt > 0:
        ratios.append(het / nalt)

# ================= FIG 1 — Ti/Tv 与变异构成 =================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 4.0))

# 左：Ti/Tv 柱状
ax1.bar(["Transitions", "Transversions"], [trans, transv], color=[C_RED, C_TEAL], width=0.55)
ax1.set_title("Ti/Tv — 2Mb chr22 slice (2504 samples)", fontsize=12, color="#222")
ax1.set_ylabel("variant count", color="#222")
for i, v in enumerate([trans, transv]):
    ax1.text(i, v + 40, str(v), ha="center", fontsize=10, color="#222")
ax1.text(0.5, max(trans, transv) * 0.62, f"overall Ti/Tv = {tstv:.2f}",
         ha="center", fontsize=11, color=C_GREY)
ax1.text(0.5, 0.14, "WGS expected Ti/Tv ~2.0-2.1; WES ~3.0-3.3",
         transform=ax1.transAxes, ha="center", fontsize=8.5, color=C_GREY)
ax1.spines[["top", "right"]].set_visible(False)
ax1.set_facecolor(C_LIGHT)

# 右：SNV/indel/多等位/其它 堆叠信息
cats = ["SNPs", "Indels", "Multiallelic", "SV/CNV"]
vals = [snp_n, indel_n, multi_n, other_n]
ax2.bar(cats, vals, color=[C_RED, C_TEAL, "#caa45a", C_GREY], width=0.6)
ax2.set_title("Variant composition (5431 records)", fontsize=12, color="#222")
ax2.set_ylabel("count", color="#222")
for i, v in enumerate(vals):
    ax2.text(i, v + 60, str(v), ha="center", fontsize=10, color="#222")
ax2.spines[["top", "right"]].set_visible(False)
ax2.set_facecolor(C_LIGHT)

fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig1_titv_composition.png"), dpi=130)
print("wrote fig1_titv_composition.png")

# ================= FIG 2 — 每样本 het/hom-alt 比值分布 =================
fig2, ax = plt.subplots(figsize=(8.5, 4.2))
nb = 30
lo = min(ratios); hi = max(ratios); w = (hi - lo) / nb
hist = [0] * nb
for r in ratios:
    b = min(nb - 1, int((r - lo) / w))
    hist[b] += 1
centers = [lo + (i + 0.5) * w for i in range(nb)]
ax.bar(centers, hist, width=w * 0.9, color=C_TEAL)
ax.set_title("Per-sample het / hom-alt ratio (n=2504)", fontsize=12, color="#222")
ax.set_xlabel("het / non-ref-hom ratio", color="#222")
ax.set_ylabel("number of samples", color="#222")
ax.axvline(1.67, color=C_RED, ls="--", lw=1.2)
ax.text(1.67, max(hist) * 0.9, f" median 1.67", color=C_RED, fontsize=9,
        va="top", ha="left")
ax.text(0.99, 0.95, f"min {lo:.2f}  max {hi:.2f}  mean {sum(ratios)/len(ratios):.2f}",
        transform=ax.transAxes, ha="right", fontsize=9, color=C_GREY)
ax.spines[["top", "right"]].set_visible(False)
ax.set_facecolor(C_LIGHT)
fig2.tight_layout()
fig2.savefig(os.path.join(HERE, "fig2_hethom_dist.png"), dpi=130)
print("wrote fig2_hethom_dist.png")
print("DONE")
