#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
044 rnaseq-qc 解析：三个 salmon quant 目录 + truth.tsv -> parsed_results.tsv + _summary.txt。
纯标准库；只读本次真跑落盘文件，不做任何插值。
"""
import json
import os
import re

BASE = os.path.dirname(os.path.abspath(__file__))
TAGS = ["A", "IU", "SF"]


def read_tsv(path):
    rows = []
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        for ln in f:
            rows.append(dict(zip(head, ln.rstrip("\n").split("\t"))))
    return rows


def load_tag(tag):
    """读取 quant_<tag>/：quant.sf + lib_format_counts.json；返回 (rows, info)。"""
    d = os.path.join(BASE, "quant_%s" % tag)
    info = {"exists": os.path.isdir(d)}
    rows = []
    qsf = os.path.join(d, "quant.sf")
    if os.path.exists(qsf):
        with open(qsf) as f:
            head = f.readline().rstrip("\n").split("\t")
            for ln in f:
                p = dict(zip(head, ln.rstrip("\n").split("\t")))
                rows.append((p["Name"], float(p["TPM"]), float(p["NumReads"])))
    lfc = os.path.join(d, "lib_format_counts.json")
    if os.path.exists(lfc):
        with open(lfc) as f:
            j = json.load(f)
        info["expected_format"] = j.get("expected_format", "NA")
        info["compatible_fragment_ratio"] = j.get("compatible_fragment_ratio", None)
        info["num_compatible_fragments"] = j.get("num_compatible_fragments", None)
        info["num_incompatible_fragments"] = j.get("num_incompatible_fragments", None)
        info["num_assigned_fragments"] = j.get("num_assigned_fragments", None)
        info["observed_SR"] = j.get("SR", None)
    # salmon 2.7.0 不写 meta_info.json：percent_mapped 从 quant 日志行解析
    logp = os.path.join(BASE, "logs", "_quant_%s.log" % tag)
    if os.path.exists(logp):
        with open(logp, errors="replace") as f:
            ms = re.findall(r"mapped (\d+) / (\d+) fragments \(([\d.]+)%\)", f.read())
        if ms:
            n_m, n_p, pct = ms[-1]
            info["num_mapped"] = int(n_m)
            info["num_processed"] = int(n_p)
            info["percent_mapped"] = float(pct)
    return rows, info


def pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    return sxy / (sxx * syy) ** 0.5 if sxx > 0 and syy > 0 else float("nan")


def rankdata(vs):
    """平均秩（处理并列）。"""
    order = sorted(range(len(vs)), key=lambda i: vs[i])
    ranks = [0.0] * len(vs)
    i = 0
    while i < len(vs):
        j = i
        while j + 1 < len(vs) and vs[order[j + 1]] == vs[order[i]]:
            j += 1
        r = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[order[k]] = r
        i = j + 1
    return ranks


truth = {r["tx_id"]: r for r in read_tsv(os.path.join(BASE, "truth.tsv"))}
data, infos = {}, {}
for tag in TAGS:
    rows, info = load_tag(tag)
    data[tag] = {name: (tpm, nr) for name, tpm, nr in rows}
    infos[tag] = info

# ---------- 逐转录本宽表 ----------
out = ["tx_id\tlength\tfrac\ttpm_expected\treads_expected"
       "\ttpm_A\tnumreads_A\ttpm_IU\tnumreads_IU\ttpm_SF\tnumreads_SF"]
for tid in sorted(truth):
    t = truth[tid]
    row = [tid, t["length"], "%.8g" % float(t["frac"]), t["tpm_expected"], t["reads_expected"]]
    for tag in TAGS:
        tpm, nr = data[tag].get(tid, (0.0, 0.0))
        row += ["%.4f" % tpm, "%d" % int(nr)]
    out.append("\t".join(row))
with open(os.path.join(BASE, "parsed_results.tsv"), "w") as f:
    f.write("\n".join(out) + "\n")

# ---------- 汇总统计 ----------
lines = []
for tag in TAGS:
    info = infos[tag]
    lines.append("== quant_%s ==" % tag)
    if not info["exists"]:
        # 运行失败（预期：-l SF 对反义 reads 不兼容）。从日志取错误信息。
        logp = os.path.join(BASE, "logs", "_quant_%s.log" % tag)
        tail = ""
        if os.path.exists(logp):
            with open(logp, errors="replace") as f:
                content = f.read().strip()
            err_lines = [ln for ln in content.splitlines()
                         if "error" in ln.lower() or "abort" in ln.lower() or "fail" in ln.lower()]
            tail = (err_lines[-1] if err_lines else (content.splitlines()[-1] if content else ""))[:200]
        lines.append("  run did not produce output; log error line: %s" % tail)
        continue
    lines.append("  expected_format = %s ; mapped = %s/%s (%s%%)"
                 % (info.get("expected_format"), info.get("num_mapped"),
                    info.get("num_processed"), info.get("percent_mapped")))
    lines.append("  lib_format_counts: compatible_ratio = %s ; compatible = %s ; "
                 "incompatible = %s ; assigned = %s ; SR-count = %s"
                 % (info.get("compatible_fragment_ratio"),
                    info.get("num_compatible_fragments"),
                    info.get("num_incompatible_fragments"),
                    info.get("num_assigned_fragments"),
                    info.get("observed_SR")))
    if data[tag]:
        tpm_sum = sum(v[0] for v in data[tag].values())
        nr_sum = sum(v[1] for v in data[tag].values())
        lines.append("  quant.sf: TPM sum = %.1f ; NumReads total = %d ; records = %d"
                     % (tpm_sum, nr_sum, len(data[tag])))

# 相关性（用 -l A 结果 vs 真值）
det = [(float(t["tpm_expected"]), data["A"][tid][0])
       for tid, t in truth.items() if data["A"].get(tid, (0, 0))[0] > 0]
n_all, n_det = len(truth), len(det)
log_exp = [__import__("math").log10(x) for x, _ in det]
log_obs = [__import__("math").log10(y) for _, y in det]
r = pearson(log_exp, log_obs)
exp_all = [float(t["tpm_expected"]) for t in truth.values()]
obs_all = [data["A"].get(tid, (0, 0))[0] for tid in truth]
rho = pearson(rankdata(exp_all), rankdata(obs_all))
ratios = sorted(y / x for x, y in det)
med_ratio = ratios[len(ratios) // 2]
lines.append("== gradient recovery (quant_A vs truth) ==")
lines.append("  detected %d/%d transcripts (TPM>0); Pearson r on log10 TPM = %.4f ; "
             "Spearman rho (all, tie-averaged) = %.4f" % (n_det, n_all, r, rho))
lines.append("  median TPM_observed/TPM_expected = %.4f ; ratio range %.4f - %.4f"
             % (med_ratio, ratios[0], ratios[-1]))

txt = "\n".join(lines)
with open(os.path.join(BASE, "_summary.txt"), "w") as f:
    f.write(txt + "\n")
print(txt)
