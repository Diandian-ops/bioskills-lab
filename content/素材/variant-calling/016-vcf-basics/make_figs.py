#!/usr/bin/env python3
"""016 vcf-basics 出图：全部基于 015 真跑产出的 variants.vcf.gz（9 变异）。
配色沿用 bioSkills 真实试用约定：brick-red #b5482f + celadon #2f7d72（色盲友好）。
标签用英文（DejaVu Sans 无 CJK 字形，中文叙述放在 md 笔记里）。

出图质量约束（由 check_figs.py 程序化校验）：
  - 所有文本包围盒必须落在 axes 内（不裁切）
  - 文本之间不得互相重叠
  - legend 不得压住数据点
数据点位置一律保持真实值，不做 jitter。
"""
import os, subprocess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BRICK = "#b5482f"
CELADON = "#2f7d72"
GREY = "#6b7280"
OUT = os.path.dirname(os.path.abspath(__file__))  # 图与脚本平铺，与 008/011 等素材目录一致
VCF_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "variants.vcf.gz")

# ---- 取数：GQ 由 PL 推导（SKILL.md: GQ = 两个最小 PL 之差，上限 99）----
# 用 bcftools query 取字段（与其余素材脚本一致，避免 cyvcf2 依赖）
def _bcft(*args):
    return subprocess.run(["bcftools"] + list(args),
                          capture_output=True, text=True).stdout

_FMT = '%POS\t%QUAL\t%INFO/DP\t[%AD]\t[%PL]\t[%DP]\t%TYPE\n'
rows = []
for line in _bcft("query", "-f", _FMT, VCF_PATH).strip().splitlines():
    pos_s, qual_s, info_dp_s, ad_s, pl_s, fmt_dp_s, vtype = line.split("\t")
    pl = [int(x) for x in pl_s.split(",")]
    s = sorted(pl)
    gq = min(s[1] - s[0], 99)
    ad = [int(x) for x in ad_s.split(",")]
    rows.append(dict(
        pos=int(pos_s), gt="0/1", pl=pl, idx0=pl.index(min(pl)), gq=gq,
        qual=float(qual_s), ad=ad, sum_ad=sum(ad),
        fmt_dp=int(fmt_dp_s), info_dp=int(info_dp_s),
        ab=(ad[1] / (ad[0] + ad[1])) if (ad[0] + ad[1]) else 0.0,
        vtype=vtype.lower(),
    ))


def _hit(a, b, pad=1.0):
    return not (a.x1 + pad < b.x0 or b.x1 + pad < a.x0
                or a.y1 + pad < b.y0 or b.y1 + pad < a.y0)


def place_labels(ax, pts, candidates, fontsize=9, color="#333333"):
    """贪心避让式标注：为每个点依次尝试候选偏移，取第一个既不越界也不与
    已放置标签相交的位置。全部失败则退回第一个候选（check_figs.py 会报 FAIL）。"""
    fig = ax.figure
    r = fig.canvas.get_renderer()
    fig.canvas.draw()
    axb = ax.get_window_extent(r)
    placed = []
    for (x, y, lab) in pts:
        chosen = None
        for (dx, dy) in candidates:
            t = ax.annotate(lab, (x, y), textcoords="offset points",
                            xytext=(dx, dy), fontsize=fontsize, color=color)
            fig.canvas.draw()
            bb = t.get_window_extent(r)
            ok = (bb.x0 >= axb.x0 + 3 and bb.x1 <= axb.x1 - 3
                  and bb.y0 >= axb.y0 + 3 and bb.y1 <= axb.y1 - 3
                  and not any(_hit(bb, p) for p in placed))
            if ok:
                placed.append(bb)
                chosen = t
                break
            t.remove()
        if chosen is None:
            ax.annotate(lab, (x, y), textcoords="offset points",
                        xytext=candidates[0], fontsize=fontsize, color=color)
    return placed


CANDIDATES = [(9, 6), (9, -14), (-32, 6), (9, 19), (-32, -15),
              (22, -25), (-50, 7), (9, -27), (30, 12)]

# ============================================================
# 图 1：QUAL vs GQ —— 两者回答不同问题，不可互换
# ============================================================
fig, ax = plt.subplots(figsize=(8.2, 5.4), dpi=150)
pos = [r["pos"] for r in rows]
qual = [r["qual"] for r in rows]
gq = [r["gq"] for r in rows]

ax.scatter(qual, gq, s=95, color=BRICK, edgecolor="white", linewidth=1.2,
           zorder=3, label="variant site")

ax.axhline(99, ls="--", lw=1.1, color=GREY, zorder=1)
ax.set_xlabel("QUAL   site-level: is there ANY variant at this position?", fontsize=10.5)
ax.set_ylabel("GQ   genotype-level: is THIS genotype call right?", fontsize=10.5)
ax.set_title("QUAL and GQ answer different questions  (9 sites, sample HG00100)",
             fontsize=12.5, pad=14)
ax.grid(alpha=0.25, ls=":", lw=0.8)
ax.set_xlim(min(qual) - 22, max(qual) + 30)
ax.set_ylim(20, 118)
ax.legend(loc="lower left", frameon=True, fontsize=9)
# 关键：先定死坐标范围 + tight_layout，再做标签避让，
# 否则 place_labels 用的是自动缩放坐标系，最终渲染位置会整体偏移。
fig.tight_layout()

# GQ cap 99：y 用 data 坐标（get_yaxis_transform 中 x 为 axes 分数、y 为 data 值）
ax.text(0.006, 99.6, "GQ cap 99  (6 of 9 sites saturate here)",
        transform=ax.get_yaxis_transform(), va="bottom", ha="left",
        fontsize=8.5, color=GREY)

place_labels(ax, [(r["qual"], r["gq"], str(r["pos"])) for r in rows], CANDIDATES)

# 两个方向相反的反例：放在图中两处空白区（中部空白 / 左下空白）
by_gq = min(rows, key=lambda r: r["gq"])      # 828: QUAL 180.83 高 / GQ 35 最低
by_qual = min(rows, key=lambda r: r["qual"])  # 2564: QUAL 55.40 最低 / GQ 84
ax.annotate("QUAL %.1f (high)\nbut GQ %d (lowest)" % (by_gq["qual"], by_gq["gq"]),
            xy=(by_gq["qual"], by_gq["gq"]), xytext=(142, 34), textcoords="data",
            ha="left", va="center", fontsize=9, color=BRICK,
            arrowprops=dict(arrowstyle="->", color=BRICK, lw=1.2))
ax.annotate("QUAL %.1f (lowest)\nbut GQ %d" % (by_qual["qual"], by_qual["gq"]),
            xy=(by_qual["qual"], by_qual["gq"]), xytext=(47, 68), textcoords="data",
            ha="left", va="center", fontsize=9, color=CELADON,
            arrowprops=dict(arrowstyle="->", color=CELADON, lw=1.2))
fig.savefig(os.path.join(OUT, "fig1_qual_vs_gq.png"))
plt.close(fig)

# ============================================================
# 图 2：sum(AD) vs INFO/DP —— sum(AD) 常小于 DP，属预期
# ============================================================
fig, ax = plt.subplots(figsize=(8.6, 5.2), dpi=150)
x = range(len(rows))
sum_ad = [r["sum_ad"] for r in rows]
info_dp = [r["info_dp"] for r in rows]

w = 0.38
ax.bar([i - w / 2 for i in x], sum_ad, w, color=CELADON,
       label="sum(AD)   allelic depths", zorder=3)
ax.bar([i + w / 2 for i in x], info_dp, w, color=BRICK,
       label="INFO/DP   raw read depth", zorder=3)

# 标出差值 > 0 的位点
place_labels(ax,
             [(i, max(r["info_dp"], r["sum_ad"]), "+%d" % (r["info_dp"] - r["sum_ad"]))
              for i, r in enumerate(rows) if r["info_dp"] - r["sum_ad"] > 0],
             [(0, 6), (0, 14), (14, 6), (-14, 6)],
             fontsize=9, color=BRICK)

ax.set_xticks(list(x))
ax.set_xticklabels([str(r["pos"]) for r in rows])
ax.set_xlabel("variant position (chr17)", fontsize=10.5)
ax.set_ylabel("read count", fontsize=10.5)
ax.set_title("sum(AD) is often lower than INFO/DP - expected, not an error",
             fontsize=12.5, pad=14)
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8)
ax.legend(loc="upper left", frameon=True, fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig2_ad_vs_dp.png"))
plt.close(fig)

# ============================================================
# 图 3：等位基因平衡 AB —— 偏离 0.5 的 het 需警惕
# ============================================================
fig, ax = plt.subplots(figsize=(8.6, 5.0), dpi=150)
ab = [r["ab"] for r in rows]
colors = [BRICK if (a < 0.2 or a > 0.8) else CELADON for a in ab]
ax.bar([str(r["pos"]) for r in rows], ab, color=colors, width=0.62, zorder=3)

ax.axhline(0.5, ls="-", lw=1.3, color="#333333", zorder=4)
ax.axhspan(0.8, 1.0, color=BRICK, alpha=0.08, zorder=0)
ax.axhspan(0.0, 0.2, color=BRICK, alpha=0.08, zorder=0)

# 数值标签放在柱底内侧（白字），避开 0.5 基线与可疑区标注
for i, r in enumerate(rows):
    ax.text(i, 0.02, "%.2f" % r["ab"], ha="center", va="bottom",
            fontsize=9.5, color="white", fontweight="bold", zorder=5)

# 右侧空白区放说明文字（柱体最高 0.83，1.12 上限之上留白）
ax.text(len(rows) - 0.45, 0.505, "0.50 = expected balance for a het call",
        va="bottom", ha="right", fontsize=9, color="#333333")
ax.text(len(rows) - 0.45, 0.815, "suspicious zone (AB > 0.8)",
        va="bottom", ha="right", fontsize=8.5, color=BRICK)

ax.set_ylim(0, 1.12)
ax.set_ylabel("allele balance  =  alt_AD / (ref_AD + alt_AD)", fontsize=10.5)
ax.set_xlabel("variant position (chr17)", fontsize=10.5)
ax.set_title("Allele balance of het calls  (red = outside 0.2-0.8)",
             fontsize=12.5, pad=14)
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig3_allele_balance.png"))
plt.close(fig)

print("figures written to", OUT)
for f in sorted(os.listdir(OUT)):
    print("  ", f)
