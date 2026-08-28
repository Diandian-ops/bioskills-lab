#!/usr/bin/env python3
"""013 fig: HISAT2 spliced alignment (real-run data). English labels to avoid CJK tofu."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import json, os

BASE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(f"{BASE}/hisat2_results.json"))
fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
teal, red, gray = "#1f9e89", "#d1495b", "#bbbbbb"

# Panel A: spliced CIGAR schematic
a = axes[0]
a.set_xlim(0, 100); a.set_ylim(0, 10); a.set_yticks([])
a.add_patch(plt.Rectangle((0, 3), 50, 4, color=teal, alpha=0.8))
a.add_patch(plt.Rectangle((58, 3), 42, 4, color=teal, alpha=0.8))
a.text(25, 8, "exon1 (50M)", ha="center", fontsize=9)
a.text(79, 8, "exon2 (50M)", ha="center", fontsize=9)
a.text(29, 1.2, "intron 800N", ha="center", fontsize=9, color=red)
a.set_title("Spliced alignment\nCIGAR: 50M800N50M")
a.axhline(5, color=red, lw=1, ls="--", xmin=0.5/100, xmax=0.58)

# Panel B: MAPQ ceiling
b = axes[1]
bars = b.bar(["HISAT2", "STAR"], [R["pe"]["max_mapq"], 255], color=[teal, gray])
b.set_ylim(0, 280); b.set_ylabel("unique-read MAPQ")
b.set_title("MAPQ scale (GATK-friendly)")
for bar, v in zip(bars, [R["pe"]["max_mapq"], 255]):
    b.text(bar.get_x()+bar.get_width()/2, v+5, str(v), ha="center", fontsize=9)
b.text(0.5, -0.28, "HISAT2=60 (no 255 reassignment)", transform=b.transAxes,
       ha="center", fontsize=8, color=teal)

# Panel C: two-pass novel junction
c = axes[2]
c.axis("off")
c.text(0.5, 0.75, "Two-pass novel junction", ha="center", fontsize=10, weight="bold")
c.text(0.5, 0.50, f"discovered: {R['two_pass']['discovered_junctions']} junction",
       ha="center", fontsize=9)
c.text(0.5, 0.30, f"pass2 CIGAR: {R['two_pass']['pass2_cigar']}", ha="center", fontsize=9, color=teal)
c.text(0.5, 0.10, f"alignment rate (PE): {R['pe']['overall_rate']}%", ha="center", fontsize=9)

fig.suptitle("hisat2-alignment — faithful reproduction", fontsize=11)
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
out = f"{BASE}/013-fig.png"
fig.savefig(out, dpi=150)
print("wrote", out)
