#!/usr/bin/env python3
"""012 fig: bwa-mem2 key behaviors (real-run data). English labels to avoid CJK tofu."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import json, os

BASE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(f"{BASE}/bwa_results.json"))
fig, axes = plt.subplots(2, 2, figsize=(11, 7))
teal, red, gray = "#1f9e89", "#d1495b", "#bbbbbb"

# (0,0) MAPQ ceiling: bwa 60 vs bowtie2 42/44
a = axes[0, 0]
vals = [R["with_rg"]["max_mapq"], 42, 44]
bars = a.bar(["bwa-mem2", "bowtie2 e2e", "bowtie2 local"], vals, color=[teal, red, red])
a.set_ylim(0, 70); a.set_ylabel("max MAPQ")
a.set_title("MAPQ ceiling (012 vs 011)")
for bar, v in zip(bars, vals):
    a.text(bar.get_x()+bar.get_width()/2, v+1.5, str(v), ha="center", fontsize=9)

# (0,1) read-group contract
a = axes[0, 1]
has = R["with_rg"]["has_RG_header"]; nolang = R["without_rg"]["has_RG_header"]
bars = a.bar(["with -R", "without -R"], [1 if has else 0, 1 if nolang else 0],
            color=[teal, red])
a.set_ylim(0, 1.3); a.set_yticks([]); a.set_title("read-group header (@RG)")
a.text(0, 0.15, "present" if has else "absent", ha="center", fontsize=9, color="white")
a.text(1, 0.15, "present" if nolang else "absent", ha="center", fontsize=9, color="white")

# (1,0) dedup strict ordering
a = axes[1, 0]
a.axis("off")
steps = "collate -> fixmate -m -> sort -> markdup"
a.text(0.5, 0.8, "duplicate-marking order", ha="center", fontsize=10, weight="bold")
a.text(0.5, 0.55, steps, ha="center", fontsize=9, color=teal)
a.text(0.5, 0.32, f"MC tag written: {R['markdup']['has_MC_tag']}", ha="center", fontsize=9)
a.text(0.5, 0.12, f"duplicates flagged: {R['markdup']['dup_flagged']}", ha="center", fontsize=9, color=red)

# (1,1) -K reproducibility
a = axes[1, 1]
a.axis("off")
a.text(0.5, 0.7, "-K 100000000 reproducibility", ha="center", fontsize=10, weight="bold")
a.text(0.5, 0.45, f"run1 md5 = {R['reproducible_K']['md5_1'][:12]}...", ha="center", fontsize=8)
a.text(0.5, 0.30, f"run2 md5 = {R['reproducible_K']['md5_2'][:12]}...", ha="center", fontsize=8)
a.text(0.5, 0.10, f"identical: {R['reproducible_K']['identical']}", ha="center", fontsize=10, color=teal)

fig.suptitle("bwa-alignment (bwa-mem2) — faithful reproduction", fontsize=11)
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
out = f"{BASE}/012-fig.png"
fig.savefig(out, dpi=150)
print("wrote", out)
