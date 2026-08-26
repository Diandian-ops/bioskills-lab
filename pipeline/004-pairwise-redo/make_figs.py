'''按 pairwise-alignment skill 点名的两个坑出图（不引入 skill 外分析）'''
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.expanduser('~/.workbuddy/skills/bioSkills-figure-quality'))
from fig_quality import C_WRONG, C_CORRECT

OUT = '/Users/zhangdiandian/RedBook/content/素材/004-pairwise'

# ---- 图1：坑一 默认 gap=0 陷阱 ----
fig, ax = plt.subplots(figsize=(6.2, 3.4), dpi=150)
cats = ['Alignment\nlength', 'Gap count']
default = [172, 55]
recommended = [149, 9]
x = range(len(cats))
w = 0.36
b1 = ax.bar([i - w/2 for i in x], default, w, label='PairwiseAligner() default (gap=0)', color=C_WRONG)
b2 = ax.bar([i + w/2 for i in x], recommended, w, label='Recommended open=-11 / extend=-1', color=C_CORRECT)
ax.set_xticks(list(x)); ax.set_xticklabels(cats)
ax.set_title('Pitfall 1: default gap penalty = 0', fontsize=12, fontweight='bold')
# legend 放到 plot 区域外（右外），避免压住任何数据或标注
ax.legend(fontsize=8, loc='upper left', bbox_to_anchor=(1.01, 1.0), frameon=False)
for b in list(b1) + list(b2):
    ax.text(b.get_x()+b.get_width()/2, b.get_height()+1, str(int(b.get_height())), ha='center', va='bottom', fontsize=9)
ax.spines[['top','right']].set_visible(False)
# 关键结论放图底部 caption（不进 plot 区域，不和任何元素撞）
fig.text(0.5, 0.02, '→ default adds 46 spurious gaps (+15.4%)',
         ha='center', fontsize=9, fontweight='bold', color=C_WRONG)
fig.tight_layout()
fig.subplots_adjust(bottom=0.16)
fig.savefig(f'{OUT}/004-fig-gappitrap.png', bbox_inches='tight')
plt.close(fig)

# ---- 图2：坑二 PID 四口径 ----
fig, ax = plt.subplots(figsize=(6.2, 3.6), dpi=150)
pids = [43.6, 46.4, 45.8, 45.0]
labels = ['PID1\n(gap-aware)', 'PID2\n(non-gap pairs)', 'PID3\n(shorter len)', 'PID4\n(mean len)']
bars = ax.bar(labels, pids, color=[C_CORRECT, C_WRONG, C_CORRECT, C_CORRECT])
ax.set_ylim(40, 50)
ax.set_ylabel('Percent identity (%)')
ax.set_title('Pitfall 2: Percent identity has 4 definitions', fontsize=12, fontweight='bold')
for b, v in zip(bars, pids):
    ax.text(b.get_x()+b.get_width()/2, v+0.2, f'{v:.1f}%', ha='center', fontsize=9)
ax.spines[['top','right']].set_visible(False)
# 底部说明放到 figure 级、x 轴标签之下（绝不压 PID1-PID4 轴标签）
fig.text(0.5, 0.02, 'Same alignment → 43.6% ~ 46.4% (range 2.8 pct; up to 11.5% on other pairs)',
         ha='center', fontsize=7.5, color='#555')
fig.tight_layout()
fig.subplots_adjust(bottom=0.20)
fig.savefig(f'{OUT}/004-fig-pid.png', bbox_inches='tight')
plt.close(fig)
print('OK figures written')
