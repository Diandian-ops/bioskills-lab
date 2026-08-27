"""Pairwise TM-score matrix heatmap for the 010 structural-alignment trial.
Real data from run_structural_alignment.py (TMalign v20220412).
Convention: TM > 0.5 = same fold -> red (中国股票/相似度上行=红); low -> green.
White background, dark text (light theme).
"""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/Users/zhangdiandian/RedBook/content/素材/010-structural-alignment"
res = json.loads(open(f"{BASE}/structural_results.json").read())

labels = ["1ubq", "1ubi", "1fmb"]
# Build symmetric matrix from TMalign results (larger normalised TM-score)
tm = res["tmalign"]
def get(a, b):
    if a == b:
        return 1.0
    key = f"{a}|{b}"
    key2 = f"{b}|{a}"
    if key in tm:
        return tm[key]["tm_score"]
    if key2 in tm:
        return tm[key2]["tm_score"]
    return None

M = [[get(a, b) for b in labels] for a in labels]

fig, ax = plt.subplots(figsize=(5.6, 4.8))
# Red = high similarity (TM large); Green = low. Reverse a red-green ramp.
cmap = plt.cm.RdYlGn_r  # high TM -> red, low -> green
im = ax.imshow(M, cmap=cmap, vmin=0.3, vmax=1.0)

ax.set_xticks(range(len(labels)))
ax.set_yticks(range(len(labels)))
ax.set_xticklabels(labels, color="#222")
ax.set_yticklabels(labels, color="#222")
ax.set_title("Pairwise TM-score (structure similarity)", color="#222", fontsize=13)
ax.set_xlabel("Structure B", color="#222")
ax.set_ylabel("Structure A", color="#222")

for i in range(len(labels)):
    for j in range(len(labels)):
        v = M[i][j]
        txt = "1.000" if i == j else f"{v:.3f}"
        ax.text(j, i, txt, ha="center", va="center", color="#111", fontsize=12)

# colorbar with TM thresholds annotated
cbar = fig.colorbar(im, ax=ax)
cbar.set_label("TM-score (>0.5 = same fold)", color="#222")
cbar.ax.yaxis.set_tick_params(color="#222")
for t in cbar.ax.get_yticklabels():
    t.set_color("#222")

# annotate threshold line at 0.5
ax.text(2.42, 0.0, "TM=0.5 threshold", color="#444", fontsize=9)

plt.tight_layout()
out = f"{BASE}/010-fig.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
print("wrote", out)
