"""009 出图：注释存活矩阵热图（来自 run_alignment_io.py 的真实回读结果）。

数据来源：annotation_survival.json
  以 test.sto（Stockholm, 7 列, 含 GS/GR/GC 注释）= 源，
  转换到各格式后回读，检测三类注释是否存活。
标签统一用英文，避免 matplotlib 默认字体缺 CJK 字形。
"""
from __future__ import annotations
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = Path(__file__).parent
data = json.loads((HERE / "annotation_survival.json").read_text())

formats = ["stockholm", "nexus", "fasta", "clustal", "phylip-relaxed"]
ann_types = ["GS", "GR", "GC"]
ann_labels = {
    "GS": "GS\n(seq-level)",
    "GR": "GR\n(residue)",
    "GC": "GC\n(column)",
}
fmt_labels = {
    "stockholm": "Stockholm",
    "nexus": "NEXUS",
    "fasta": "FASTA",
    "clustal": "Clustal",
    "phylip-relaxed": "PHYLIP-relaxed",
}

# 构建矩阵：1=存活, 0=丢失
mat = []
for fmt in formats:
    row = [1 if data[fmt].get(at) else 0 for at in ann_types]
    mat.append(row)

fig, ax = plt.subplots(figsize=(6.4, 4.4), dpi=150)
cmap = matplotlib.colors.ListedColormap(["#e63946", "#2a9d8f"])
ax.imshow(mat, cmap=cmap, vmin=0, vmax=1, aspect="auto")

for i, fmt in enumerate(formats):
    for j, at in enumerate(ann_types):
        txt = "KEPT" if mat[i][j] == 1 else "DROPPED"
        ax.text(j, i, txt, ha="center", va="center",
                color="white", fontsize=11, fontweight="bold")

ax.set_xticks(range(len(ann_types)))
ax.set_xticklabels([ann_labels[at] for at in ann_types], fontsize=10)
ax.set_yticks(range(len(formats)))
ax.set_yticklabels([fmt_labels[f] for f in formats], fontsize=11)

ax.set_title("009 alignment-io: annotation survival matrix\n"
             "(Stockholm source -> write to format -> read back)",
             fontsize=12, fontweight="bold", pad=12)
ax.set_xlabel("Annotation type", fontsize=10)
ax.set_ylabel("Target format", fontsize=10)

legend = [Patch(facecolor="#2a9d8f", label="preserved"),
          Patch(facecolor="#e63946", label="silently dropped")]
ax.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, -0.16),
          ncol=2, frameon=False, fontsize=9)

ax.set_xticks([x - 0.5 for x in range(1, len(ann_types))], minor=True)
ax.set_yticks([y - 0.5 for y in range(1, len(formats))], minor=True)
ax.grid(which="minor", color="white", linewidth=2)
ax.tick_params(which="minor", length=0)

plt.tight_layout()
out = HERE / "009-fig.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"wrote {out} ({out.stat().st_size} bytes)")
