#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
041 fastp-workflow 报告解析：从三组 fastp JSON 报告提取 before/after 真实指标，
落盘 _metrics.tsv（制表符分隔，供笔记与出图引用）。
"""
import json
import os

BASE = os.path.dirname(os.path.abspath(__file__))
RUNS = [
    ("full", "report_full.json"),
    ("noadapter", "report_noadapter.json"),
    ("nolenfilter", "report_nolenfilter.json"),
    ("nocut", "report_nocut.json"),
]

# (JSON 路径, 展示名)
FIELDS = [
    (("summary", "before_filtering", "total_reads"), "before_total_reads"),
    (("summary", "before_filtering", "total_bases"), "before_total_bases"),
    (("summary", "before_filtering", "q20_rate"), "before_q20_rate"),
    (("summary", "before_filtering", "q30_rate"), "before_q30_rate"),
    (("summary", "before_filtering", "read1_mean_length"), "before_read1_mean_len"),
    (("summary", "after_filtering", "total_reads"), "after_total_reads"),
    (("summary", "after_filtering", "total_bases"), "after_total_bases"),
    (("summary", "after_filtering", "q20_rate"), "after_q20_rate"),
    (("summary", "after_filtering", "q30_rate"), "after_q30_rate"),
    (("summary", "after_filtering", "read1_mean_length"), "after_read1_mean_len"),
    (("adapter_cutting", "adapter_trimmed_reads"), "adapter_trimmed_reads"),
    (("adapter_cutting", "adapter_trimmed_bases"), "adapter_trimmed_bases"),
    (("filtering_result", "passed_filter_reads"), "passed_filter_reads"),
    (("filtering_result", "low_quality_reads"), "low_quality_reads"),
    (("filtering_result", "too_many_N_reads"), "too_many_N_reads"),
    (("filtering_result", "too_short_reads"), "too_short_reads"),
    (("filtering_result", "corrected_reads"), "corrected_reads"),
    (("filtering_result", "corrected_bases"), "corrected_bases"),
    (("duplication", "rate"), "duplication_rate"),
]


def get_nested(d, path):
    cur = d
    for k in path:
        if not isinstance(cur, dict) or k not in cur:
            return None
        cur = cur[k]
    return cur


def main():
    rows = {}
    for tag, fn in RUNS:
        with open(os.path.join(BASE, fn), "r", encoding="utf-8") as f:
            rep = json.load(f)
        vals = {}
        for path, name in FIELDS:
            vals[name] = get_nested(rep, path)
        rows[tag] = vals

    with open(os.path.join(BASE, "_metrics.tsv"), "w", encoding="utf-8") as f:
        header = ["metric"] + [t for t, _ in RUNS]
        f.write("\t".join(header) + "\n")
        for _, name in FIELDS:
            line = [name]
            for t, _ in RUNS:
                v = rows[t][name]
                line.append("" if v is None else ("%.6g" % v if isinstance(v, float) else str(v)))
            f.write("\t".join(line) + "\n")

    # stdout 摘要
    for t, _ in RUNS:
        v = rows[t]
        print("== %s ==" % t)
        for name in [n for _, n in FIELDS]:
            print("  %-24s %s" % (name, v[name]))
    print("written: _metrics.tsv")


if __name__ == "__main__":
    main()
