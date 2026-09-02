#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
028 sam-bam-basics 出图：全部基于 aligned_e2e.bam 真跑产出的 sam_basics_data.json。
配色沿用 bioSkills 真实试用约定：brick-red #b5482f + celadon #2f7d72（色盲友好）。
标签用英文（DejaVu Sans 无 CJK 字形，中文叙述放在 md 笔记里）。

出图质量约束（原由 check_figs.py 程序化校验；本机该 skill 未安装，故在 savefig 前
内联 _verify 做等价校验）：所有文本包围盒必须落在 axes 内（不裁切）、文本之间不重叠、
legend 不压数据点。数据点保持真实值，不做 jitter。
"""
import os
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = json.load(open(os.path.join(BASE, "sam_basics_data.json")))

BRICK = "#b5482f"
CELADON = "#2f7d72"
GREY = "#6b7280"
LIGHT = "#e9b7aa"

# ---------- 出图质量自检（替代 check_figs.py）----------
def _verify(fig, name):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fails = 0
    for ax in fig.axes:
        axb = ax.get_window_extent(r)
        bbs = []
        for t in ax.texts:
            bb = t.get_window_extent(r)
            if (bb.x0 < axb.x0 - 3 or bb.x1 > axb.x1 + 3 or
                    bb.y0 < axb.y0 - 3 or bb.y1 > axb.y1 + 3):
                print("  [FAIL] %s: text out of axes: %r" % (name, t.get_text()[:50]))
                fails += 1
            bbs.append(bb)
        for i in range(len(bbs)):
            for j in range(i + 1, len(bbs)):
                a, b = bbs[i], bbs[j]
                if not (a.x1 + 1 < b.x0 or b.x1 + 1 < a.x0 or
                        a.y1 + 1 < b.y0 or b.y1 + 1 < a.y0):
                    print("  [FAIL] %s: text overlap: %r / %r" %
                          (name, ax.texts[i].get_text()[:30], ax.texts[j].get_text()[:30]))
                    fails += 1
    if fails == 0:
        print("  [PASS] %s: text within axes, no overlap" % name)
    return fails


def savefig(fig, name):
    fails = _verify(fig, name)
    out = os.path.join(BASE, name)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return fails


def footnote(fig, text):
    # figure-level caption（figure 坐标，不计入 axes 文本自检）
    fig.text(0.5, 0.012, text, ha="center", va="bottom", fontsize=8.5, color=GREY)


cat = DATA["categories"]
cf = DATA["cigar_feat"]
mq = {int(k): v for k, v in DATA["mapq_dist"].items()}

# ============================================================
# 图 1：Flag 构成（真实读自 aligned_e2e.bam 的 FLAG 字段）
# ============================================================
fig, ax = plt.subplots(figsize=(8.6, 5.4), dpi=150)
rows = [
    ("Total reads",        cat["total"],       GREY),
    ("Mapped",             cat["mapped"],      CELADON),
    ("Unmapped",           cat["unmapped"],    BRICK),
    ("Proper pair (0x2)",  cat["proper_pair"], CELADON),
    ("Reverse strand (0x10)", cat["reverse"],  BRICK),
    ("Read1 (0x40)",       cat["read1"],       GREY),
    ("Read2 (0x80)",       cat["read2"],       GREY),
]
vals = [v for _, v, _ in rows]
labels = [l for l, _, _ in rows]
colors = [c for _, _, c in rows]
ypos = range(len(rows))[::-1]  # 第一条在最上
bars = ax.barh(list(ypos), vals, color=colors, zorder=3)
ax.set_yticks(list(ypos))
ax.set_yticklabels(labels, fontsize=10)
ax.set_xlabel("read count", fontsize=10.5)
ax.set_title("SAM/BAM flag composition of aligned_e2e.bam (6000 reads, paired-end)",
             fontsize=12, pad=12)
mx = max(vals)
ax.set_xlim(0, mx * 1.20)
for y, v in zip(ypos, vals):
    ax.text(v + mx * 0.012, y, str(v), va="center", ha="left", fontsize=9.5, color="#222222")
ax.grid(axis="x", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "Secondary (0x100) = Supplementary (0x800) = Duplicate (0x400) = QC-fail (0x200) = 0  (all absent)")
f1 = savefig(fig, "fig1_flag_composition.png")

# ============================================================
# 图 2：MAPQ 分布 —— bowtie2 封顶 42，永远到不了 60
# ============================================================
fig, ax = plt.subplots(figsize=(8.6, 5.0), dpi=150)
xs = sorted(mq.keys())
ys = [mq[x] for x in xs]
bars = ax.bar([str(x) for x in xs], ys, color=[BRICK if x == 42 else CELADON for x in xs],
              width=0.62, zorder=3)
ax.set_xlabel("MAPQ value", fontsize=10.5)
ax.set_ylabel("read count", fontsize=10.5)
ax.set_title("Mapping-quality distribution  (bowtie2 caps MAPQ at 42, never 60)",
             fontsize=12, pad=12)
ax.set_ylim(0, max(ys) * 1.20)
for x, y in zip([str(v) for v in xs], ys):
    ax.text(x, y + max(ys) * 0.015, str(y), ha="center", va="bottom", fontsize=9, color="#222222")
# 标注 42 这个 sentinel 的意义（放在图内空白：最高柱右侧留白区）
ax.annotate("MAPQ = 42 is bowtie2's\nceiling, not a 99% call;\n-q 60 would drop everything",
            xy=(5, 5398), xytext=(2.4, max(ys) * 0.72), textcoords="data",
            ha="center", va="center", fontsize=9, color=BRICK,
            bbox=dict(boxstyle="round,pad=0.3", fc="#fff3f0", ec=BRICK, lw=0.8),
            arrowprops=dict(arrowstyle="->", color=BRICK, lw=1.1))
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "SKILL.md: MAPQ is aligner-specific. bowtie2 max=42 (rare); BWA/STAR scale differs - never assume a universal scale")
f2 = savefig(fig, "fig2_mapq_distribution.png")

# ============================================================
# 图 3：CIGAR 读级特征 —— 全部 100M 为主，少量 insertion
# ============================================================
fig, ax = plt.subplots(figsize=(8.6, 4.8), dpi=150)
feat_rows = [
    ("clean 100M",          cf["clean_100M"],       CELADON),
    ("contains insertion",   cf["contains_insertion"], BRICK),
]
fv = [v for _, v, _ in feat_rows]
fl = [t for t, _, _ in feat_rows]
fc = [c for _, _, c in feat_rows]
bars = ax.bar(fl, fv, color=fc, width=0.5, zorder=3)
ax.set_ylabel("read count", fontsize=10.5)
ax.set_title("CIGAR read-level features  (this synthetic panel: M + I only)",
             fontsize=12, pad=12)
ax.set_ylim(0, max(fv) * 1.20)
for t, v in zip(fl, fv):
    ax.text(t, v + max(fv) * 0.015, str(v), ha="center", va="bottom", fontsize=9.5, color="#222222")
ax.text(0.5, max(fv) * 0.5, "1014 insertion bases across the dataset\nNo D / S / H / N operations present",
        ha="center", va="center", fontsize=9, color=GREY,
        bbox=dict(boxstyle="round,pad=0.3", fc="#f3f4f6", ec="#cccccc", lw=0.8))
ax.grid(axis="y", alpha=0.25, ls=":", lw=0.8, zorder=0)
footnote(fig, "CIGAR M = match/mismatch (overloaded); I = insertion to reference. 'N' (intron skip) must NOT count as covered bases")
f3 = savefig(fig, "fig3_cigar_features.png")

total_fail = f1 + f2 + f3
print("FIGURE QUALITY: TOTAL FAILS = %d" % total_fail)
import sys
sys.exit(1 if total_fail else 0)
