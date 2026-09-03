#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
040 quality-filtering 指标解析：从真跑产物（各工具输出 FASTQ / fastp JSON /
cutadapt 报告 / trimmomatic 日志）提取真实指标，落盘 results.json。
指标：保留 reads 数/率、保留 bases 率、修剪后平均 Q、修剪后平均读长、
逐 cycle 平均质量曲线、读长分布（分箱）。
"""
import gzip
import json
import os

BASE = os.path.dirname(os.path.abspath(__file__))

CONFIGS = [
    ("fastp_filter", "out_fastp_filter.fq.gz"),
    ("fastp_cutright", "out_fastp_cutright.fq.gz"),
    ("cutadapt", "out_cutadapt.fq.gz"),
    ("trimm_sw", "out_trimm_sw.fq.gz"),
    ("trimm_maxinfo", "out_trimm_maxinfo.fq.gz"),
]

LEN_BINS = [(36, 60), (61, 85), (86, 110), (111, 135), (136, 150)]


def fastq_stats(path):
    """统计一个 FASTQ(gz)：reads、平均 Q、平均长度、逐 cycle 平均 Q、长度分箱。"""
    n = 0
    tot_q = 0
    tot_len = 0
    cyc_q = []
    cyc_n = []
    bins = [0] * len(LEN_BINS)
    with gzip.open(path, "rt") as f:
        while True:
            name = f.readline()
            if not name:
                break
            seq = f.readline().rstrip("\n")
            f.readline()
            qual = f.readline().rstrip("\n")
            n += 1
            L = len(seq)
            tot_len += L
            tot_q += sum(ord(c) - 33 for c in qual)
            for c, ch in enumerate(qual):
                if c >= len(cyc_q):
                    cyc_q.append(0.0)
                    cyc_n.append(0)
                cyc_q[c] += ord(ch) - 33
                cyc_n[c] += 1
            for i, (lo, hi) in enumerate(LEN_BINS):
                if lo <= L <= hi:
                    bins[i] += 1
                    break
    return {
        "reads": n,
        "mean_q": round(tot_q / max(tot_len, 1), 2),
        "mean_len": round(tot_len / max(n, 1), 2),
        "per_cycle_q": [round(q / c, 2) for q, c in zip(cyc_q, cyc_n)],
        "len_bins": bins,
    }


def main():
    res = {"input": fastq_stats(os.path.join(BASE, "input_grad.fq.gz"))}
    res["input"]["len_bin_labels"] = ["%d-%d" % b for b in LEN_BINS]

    # fastp JSON 报告
    for tag in ("filter", "cutright"):
        with open(os.path.join(BASE, "_fastp_%s.json" % tag)) as f:
            j = json.load(f)
        s = j["summary"]
        res["fastp_%s_report" % tag] = {
            "before_reads": s["before_filtering"]["total_reads"],
            "before_bases": s["before_filtering"]["total_bases"],
            "after_reads": s["after_filtering"]["total_reads"],
            "after_bases": s["after_filtering"]["total_bases"],
            "after_q30_bases": s["after_filtering"]["q30_bases"],
            "after_q30_rate": round(s["after_filtering"]["q30_rate"], 4),
            "filtering_result": j.get("filtering_result", {}),
        }

    # cutadapt 报告
    with open(os.path.join(BASE, "_cutadapt_report.txt")) as f:
        ca = f.read()
    res["cutadapt_report"] = {}
    for line in ca.splitlines():
        if line.startswith("Total reads processed"):
            res["cutadapt_report"]["reads_processed"] = int(line.split(":")[1].strip().replace(",", ""))
        if line.startswith("Reads written"):
            res["cutadapt_report"]["reads_written"] = int(line.split(":")[1].split("(")[0].strip().replace(",", ""))
        if line.startswith("Total basepairs processed"):
            res["cutadapt_report"]["bp_processed"] = int(line.split(":")[1].strip().replace(",", "").split()[0])
        if line.startswith("Quality-trimmed"):
            res["cutadapt_report"]["bp_quality_trimmed"] = int(line.split(":")[1].split("(")[0].strip().replace(",", "").split()[0])
        if line.startswith("Total written"):
            res["cutadapt_report"]["bp_written"] = int(line.split(":")[1].split("(")[0].strip().replace(",", "").split()[0])

    # trimmomatic 日志
    for tag in ("sw", "maxinfo"):
        log = os.path.join(BASE, "_trimm_%s.log" % tag)
        if os.path.exists(log):
            with open(log) as f:
                txt = f.read()
            for line in txt.splitlines():
                if "Input Reads" in line:
                    res["trimm_%s_log" % tag] = line.strip()

    # 各工具输出 FASTQ 真实统计
    for name, fn in CONFIGS:
        res[name] = fastq_stats(os.path.join(BASE, fn))

    n_in = res["input"]["reads"]
    for name, _ in CONFIGS:
        res[name]["retention_pct"] = round(100.0 * res[name]["reads"] / n_in, 2)
        res[name]["bases_retention_pct"] = round(
            100.0 * res[name]["mean_len"] * res[name]["reads"] /
            (res["input"]["mean_len"] * n_in), 2)

    with open(os.path.join(BASE, "results.json"), "w") as f:
        json.dump(res, f, indent=2)

    # 摘要打印
    print("input: %d reads, mean Q %.2f, mean len %.1f"
          % (n_in, res["input"]["mean_q"], res["input"]["mean_len"]))
    for name, _ in CONFIGS:
        r = res[name]
        print("%-15s kept %6d (%5.2f%%)  mean Q %.2f  mean len %.1f  bases kept %.2f%%"
              % (name, r["reads"], r["retention_pct"], r["mean_q"],
                 r["mean_len"], r["bases_retention_pct"]))
    print("fastp_filter filtering_result:", res["fastp_filter_report"]["filtering_result"])
    print("fastp_cutright filtering_result:", res["fastp_cutright_report"]["filtering_result"])
    print("cutadapt:", res["cutadapt_report"])


if __name__ == "__main__":
    main()
