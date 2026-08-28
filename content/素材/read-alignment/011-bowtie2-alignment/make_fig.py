#!/usr/bin/env python3
"""011 fig: bowtie2 key behaviors (real-run data). English labels to avoid CJK tofu."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import json, os

BASE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(f"{BASE}/bowtie2_results.json"))

fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
teal, red, gray = "#1f9e89", "#d1495b", "#bbbbbb"

# Panel A: adapter contamination, e2e vs local
a = axes[0]
vals = [R["adapter_contam"]["e2e_rate"], R["adapter_contam"]["local_rate"]]
bars = a.bar(["end-to-end", "--local"], vals, color=[red, teal])
a.set_ylim(0, 105); a.set_ylabel("alignment rate (%)")
a.set_title("Adapter-contaminated reads")
for b, v in zip(bars, vals):
    a.text(b.get_x()+b.get_width()/2, v+1.5, f"{v:.1f}%", ha="center", fontsize=9)
a.text(0.5, -0.28, f"+{R['adapter_contam']['recovered_by_local']} pp by soft-clip",
       transform=a.transAxes, ha="center", fontsize=8, color=teal)

# Panel B: MAPQ cap
b = axes[1]
caps = [R["mapq_cap"]["e2e_max"], R["mapq_cap"]["local_max"], R["mapq_cap"]["bwa_equivalent"]]
bars = b.bar(["e2e", "local", "BWA"], caps, color=[teal, teal, gray])
b.set_ylim(0, 70); b.set_ylabel("max MAPQ")
b.set_title("MAPQ ceiling")
for bar, v in zip(bars, caps):
    b.text(bar.get_x()+bar.get_width()/2, v+1.5, str(v), ha="center", fontsize=9)
b.text(0.5, -0.28, "Bowtie2 never reaches 60", transform=b.transAxes,
       ha="center", fontsize=8, color=red)

# Panel C: sensitivity presets
c = axes[2]
pr = R["presets"]
keys = ["--very-fast", "--sensitive", "--very-sensitive"]
vals = [pr[k] for k in keys]
bars = c.bar([k.replace("--", "") for k in keys], vals, color=teal)
c.set_ylim(90, 102); c.set_ylabel("alignment rate (%)")
c.set_title("Sensitivity presets")
for bar, v in zip(bars, vals):
    c.text(bar.get_x()+bar.get_width()/2, v+0.3, f"{v:.1f}%", ha="center", fontsize=9)

fig.suptitle("bowtie2-alignment — faithful reproduction", fontsize=11)
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
out = f"{BASE}/011-fig.png"
fig.savefig(out, dpi=150)
print("wrote", out)
