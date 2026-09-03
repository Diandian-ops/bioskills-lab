#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
042 umi-processing 结果解析（在 bio 环境跑，借助 samtools；stdlib only）。
输入：extract.log / R1_umi.fq.gz / sorted.bam / dedup_*.bam / coordonly.bam /
      molecules.tsv / truth.json / stats_dir_*.tsv
输出：results.json / results.txt（供笔记与出图引用的真实数字）
"""
import gzip
import json
import os
import subprocess
from collections import Counter

BASE = os.path.dirname(os.path.abspath(__file__))
os.chdir(BASE)


def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip()


R = {}

# ---------- truth ----------
R["truth"] = json.load(open("truth.json"))

# ---------- extract 结果 ----------
umi_len_hist = Counter()
n_extracted = 0
n_umi_with_n = 0
distinct_umi = set()
with gzip.open("R1_umi.fq.gz", "rt") as f:
    for i, line in enumerate(f):
        if i % 4 != 0:
            continue
        name = line[1:].strip()
        n_extracted += 1
        umi = name.split("_")[-1]
        umi_len_hist[len(umi)] += 1
        distinct_umi.add(umi)
        if "N" in umi:
            n_umi_with_n += 1
R["extract"] = dict(
    reads_extracted=n_extracted, umi_len_hist={str(k): v for k, v in sorted(umi_len_hist.items())},
    n_umi_with_n=n_umi_with_n, distinct_umi_seqs=len(distinct_umi))
R["extract_log_lines"] = [l.rstrip() for l in open("extract.log", encoding="utf-8", errors="replace")
                          if l.strip() and not l.startswith("#")]

# ---------- 比对 ----------
flagstat = sh("samtools flagstat sorted.bam")
mapped = int(sh("samtools view -c -F 0x904 sorted.bam"))
total = int(sh("samtools view -c -F 0x900 sorted.bam"))
R["align"] = dict(reads_total=total, reads_mapped=mapped,
                  mapped_pct=round(100.0 * mapped / total, 2) if total else None,
                  flagstat_head=flagstat.splitlines()[:8])

# markdup 坐标法统计
mk = open("05_markdup_stats.log").read()
R["markdup"] = {}
for key in ("DUPLICATE PAIR", "ESTIMATED_LIBRARY_SIZE"):
    for ln in mk.splitlines():
        if ln.startswith(key + ":"):
            R["markdup"][key] = ln.split(":", 1)[1].strip()

# ---------- 各方法分子估计（primary R1 计数，排除 dup/secondary/supplementary）----------
def count_r1(bam):
    return int(sh("samtools view -c -f 0x40 -F 0xD00 " + bam))

R["methods"] = {
    "coordinate_only": count_r1("coordonly.bam"),
    "unique": count_r1("dedup_unique.bam"),
    "directional": count_r1("dedup_directional.bam"),
    "cluster": count_r1("dedup_cluster.bam"),
}
n_truth = R["truth"]["n_molecules"]
for k, v in R["methods"].items():
    R["methods_err_pct_" + k] = round(100.0 * (v - n_truth) / n_truth, 2)
R["dedup_directional_log"] = [l.rstrip() for l in open("dedup_directional.log", encoding="utf-8", errors="replace")
                              if l.strip() and not l.startswith("#")]
pairs_in = n_extracted
R["retention_pct_directional"] = round(100.0 * R["methods"]["directional"] / pairs_in, 2)

# ---------- 同坐标不同 UMI 干扰对的实证检查 ----------
mol_rows = [l.rstrip("\n").split("\t") for l in open("molecules.tsv")][1:]
mol = [dict(mol_id=m[0], contig=m[1], start=int(m[2]), umi=m[3],
            group=m[5], role=m[6], ptype=m[7]) for m in mol_rows]
groups = {}
for m in mol:
    if m["group"]:
        groups.setdefault(m["group"], []).append(m)
R["n_true_keys"] = len(set((m["contig"], m["start"], m["umi"]) for m in mol))


def bam_coord_umi_keys(bam):
    """从 BAM R1 primary 提取 (contig, POS, umi-from-name)。"""
    out = subprocess.Popen("samtools view %s" % bam, shell=True,
                           stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    keys = set()
    for ln in out.stdout:
        p = ln.split("\t")
        flag = int(p[1])
        if not (flag & 0x40) or (flag & 0xD00):
            continue
        keys.add((p[2], int(p[3]), p[0].split("_")[-1]))
    out.stdout.close()
    out.wait()
    return keys


def pair_keep_stats(keys):
    """keys: BAM 里的 (contig, POS1based, observed_umi)；分子保留判定 =
    该分子的母本 UMI 是否有代表读段（exact match，错误子代不计）。"""
    h1 = h1_both = dv = 0
    for g, ms in groups.items():
        kept = [((m["contig"], m["start"] + 1, m["umi"]) in keys) for m in ms]
        if ms[0]["ptype"] == "hamming1":
            h1 += 1
            h1_both += (all(kept))
        else:
            dv += 1
    return h1, h1_both, dv


for name, bam in (("directional", "dedup_directional.bam"),
                  ("unique", "dedup_unique.bam"),
                  ("cluster", "dedup_cluster.bam")):
    h1, both, dv = pair_keep_stats(bam_coord_umi_keys(bam))
    R["pairs_" + name] = dict(hamming1_pairs=h1, hamming1_both_kept=both, diverse_pairs=dv)

# ---------- stats 文件（directional）----------
# stats_dir_edit_distance.tsv 表头: unique unique_null directional directional_null edit_distance
ed_rows = []
for ln in open("stats_dir_edit_distance.tsv"):
    if ln.startswith("#") or ln.startswith("unique\t"):
        continue
    p = ln.split()
    if len(p) >= 5:
        ed_rows.append([p[4], int(p[0]), float(p[1]), int(p[2]), float(p[3])])
R["edit_distance"] = ed_rows
if ed_rows:
    R["edit_distance_single_umi"] = dict(unique=ed_rows[0][1], directional=ed_rows[0][3])

# stats_dir_per_umi_per_position.tsv 表头: counts instances_pre instances_post
fam_rows = []
for ln in open("stats_dir_per_umi_per_position.tsv"):
    if ln.startswith("#") or ln.startswith("counts"):
        continue
    p = ln.split()
    if len(p) >= 3:
        try:
            fam_rows.append([int(p[0]), int(p[1]), int(p[2])])
        except ValueError:
            pass
R["family_size"] = fam_rows
R["family_size_totals"] = dict(pre=sum(r[1] for r in fam_rows), post=sum(r[2] for r in fam_rows))

# ---------- 落盘 ----------
json.dump(R, open("results.json", "w"), indent=2)
with open("results.txt", "w") as f:
    f.write("== truth ==\n" + json.dumps(R["truth"], indent=2) + "\n")
    f.write("== extract ==\n" + json.dumps(R["extract"], indent=2) + "\n")
    f.write("extract.log 结果行:\n")
    for l in R["extract_log_lines"]:
        f.write("  %s\n" % l)
    f.write("== align ==\n" + json.dumps(R["align"], indent=2) + "\n")
    f.write("== markdup ==\n" + json.dumps(R["markdup"], indent=2) + "\n")
    f.write("== methods (molecule estimates, truth=%d) ==\n" % n_truth)
    for k, v in R["methods"].items():
        f.write("  %-16s %6d  (err %+7.2f%%)\n" % (k, v, R["methods_err_pct_" + k]))
    f.write("retention directional = %.2f%% of read pairs\n" % R["retention_pct_directional"])
    f.write("== hamming1 干扰对实证（both kept / total pairs）==\n")
    for name in ("directional", "unique", "cluster"):
        d = R["pairs_" + name]
        f.write("  %-12s hamming1 both kept: %d/%d  diverse pairs: %d\n"
                % (name, d["hamming1_both_kept"], d["hamming1_pairs"], d["diverse_pairs"]))
    f.write("== edit_distance: distance | unique_obs unique_null | dir_obs dir_null ==\n")
    for row in ed_rows:
        f.write("  %-11s unique=%-6d null=%-8.2f dir=%-6d null=%.2f\n" % tuple(row))
    f.write("== family size: counts instances_pre instances_post ==\n")
    for row in fam_rows:
        f.write("  size=%d pre=%d post=%d\n" % tuple(row))
    f.write("  totals pre=%d post=%d\n" % (R["family_size_totals"]["pre"],
                                           R["family_size_totals"]["post"]))

print(open("results.txt").read())
