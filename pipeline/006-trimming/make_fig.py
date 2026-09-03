import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pathlib import Path
import sys, os

# 引入 figure-quality helper（legend 安全区 / 无文字重叠 / 语义配色）
sys.path.insert(0, os.path.expanduser('~/.workbuddy/skills/bioSkills-figure-quality'))
from fig_quality import save_clean, place_legend, caption, C_CORRECT, C_WRONG

BASE = "/Users/zhangdiandian/RedBook"
LOG = Path(BASE + "/pipeline/006-trimming/trimmed.fasta.log")
OUT = Path(BASE + "/content/素材/alignment/006-alignment-trimming/006-fig.png")

# 解析 clipkit --log：每行 "pos keep|trim classification gap_proportion"
cols = []
for line in LOG.read_text().splitlines():
    line = line.strip()
    if not line:
        continue
    parts = line.split()
    if len(parts) < 2:
        continue
    cols.append((int(parts[0]), parts[1]))

n = len(cols)
kept = sum(1 for _, s in cols if s == "keep")
trimmed = n - kept
pct = kept / n * 100

fig, ax = plt.subplots(figsize=(10, 2.8))
for i, (pos, status) in enumerate(cols):
    color = C_CORRECT if status == "keep" else C_WRONG
    ax.add_patch(plt.Rectangle((i, 0), 1, 1, color=color, ec="white", lw=0.6))
    ax.text(i + 0.5, 0.5, str(pos), ha="center", va="center",
            fontsize=9, color="white" if status == "trim" else "black")

ax.set_xlim(0, n)
ax.set_ylim(0, 1)
ax.set_yticks([])
ax.set_xlabel("Alignment column", fontsize=9)
ax.set_title(
    f"ClipKIT kpic-smart-gap \u2014 {n} sites \u2192 {kept} kept ({pct:.1f}%), {trimmed} trimmed",
    fontsize=10, pad=8,
)

# legend 用 fig_quality helper 放到 plot 区域外右侧，不压任何东西
place_legend(ax, loc='outside right',
             handles=[Patch(color=C_CORRECT, label="kept"),
                      Patch(color=C_WRONG, label="trimmed")])

# 底部说明用 fig 级 caption，在 x 轴标签之下
caption(fig, f"Retained {pct:.1f}% \u2014 within safe-trim regime (<20% removed, per SKILL.md 20/40% rule)")

save_clean(fig, OUT, dpi=160)
print(f"SAVED {OUT} | kept={kept} trimmed={trimmed} retained={pct:.1f}%")
