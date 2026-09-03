#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
039 adapter-trimming 结果解析：以逐读段真值（truth.tsv.gz）为基准，
统一口径评估 cutadapt 与 Trimmomatic（两种 ILLUMINACLIP 口径）在 5%/20%/40% 接头梯度上的表现。

工具变体：
  cutadapt            SKILL.md PE 口径: -a/-A AGATCGGAAGAGC -m 20:20
  trimmomatic_skill   SKILL.md 字面参数串 ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads
  trimmomatic_true    修正布尔位   ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:true

判定口径（三工具完全一致，按 read_id 对齐，R1/R2 分别处理）：
  - 命中（hit）     : 输出中该 read_id 存在且序列变短
  - 漏检（untouched）: 输出中该 read_id 存在且长度不变
  - 丢弃（dropped） : 输出中该 read_id 不存在（MINLEN/keepBothReads 行为等）

产出: results.json / results.tsv / bins.tsv
"""
import gzip
import json
import os
from collections import defaultdict

BASE = os.path.dirname(os.path.abspath(__file__))
LEVELS = ["5p", "20p", "40p"]
TOOLS = ["cutadapt", "trimmomatic_skill", "trimmomatic_true"]
BINS = [(5, 9), (10, 15), (16, 24), (25, 33)]


def read_fq_lens(path, store):
    """把 fastq.gz 中 read_id -> seq_len 存入 store。"""
    with gzip.open(path, "rt") as f:
        while True:
            head = f.readline()
            if not head:
                break
            seq = f.readline().rstrip("\n")
            f.readline()
            f.readline()
            store[head[1:].rstrip("\n")] = len(seq)


def main():
    truth = {}          # (read_id, mate) -> (has_adapter, adapter_bases, read_len)
    with gzip.open(os.path.join(BASE, "truth.tsv.gz"), "rt") as f:
        f.readline()
        for line in f:
            rid, mate, has_ad, ad_bases, rlen = line.rstrip("\n").split("\t")
            truth[(rid, mate)] = (int(has_ad), int(ad_bases), int(rlen))

    results = {"levels": {}, "bins": {}, "cutadapt_report": {}}
    tsv_rows = ["level\ttool\treads_in\treads_out\treads_dropped\t"
                "truth_adapter_reads\tadapter_hit_reads\trecall_pct\t"
                "false_clip_reads\ttruth_adapter_bases\tbases_removed\t"
                "bases_removed_pct\tbases_removed_per_read"]
    bin_tot = defaultdict(dict)   # bin -> {tool_hits: x, tool_truth: y, ...}
    by_remnant = defaultdict(lambda: defaultdict(list))   # tool -> adapter_bases -> [removed bp]
    bin_rows = ["bin\tn_truth\tn_cutadapt_hit\tcutadapt_recall_pct\t"
                "n_trimmomatic_hit\ttrimmomatic_recall_pct"]

    for lev in LEVELS:
        ins = {"R1": {}, "R2": {}}
        read_fq_lens(os.path.join(BASE, "grad%s_R1.fq.gz" % lev), ins["R1"])
        read_fq_lens(os.path.join(BASE, "grad%s_R2.fq.gz" % lev), ins["R2"])

        ca = {"R1": {}, "R2": {}}
        read_fq_lens(os.path.join(BASE, "ca_%s_R1.fq.gz" % lev), ca["R1"])
        read_fq_lens(os.path.join(BASE, "ca_%s_R2.fq.gz" % lev), ca["R2"])
        tms = {"R1": {}, "R2": {}}
        tmk = {"R1": {}, "R2": {}}
        for mate, tag in (("R1", "R1"), ("R2", "R2")):
            for part in ("p", "u"):
                read_fq_lens(os.path.join(BASE, "tm_%s_%s_%s.fq.gz" % (lev, tag, part)),
                             tms[mate])
                read_fq_lens(os.path.join(BASE, "tmk_%s_%s_%s.fq.gz" % (lev, tag, part)),
                             tmk[mate])
        outs = {"cutadapt": ca, "trimmomatic_skill": tms, "trimmomatic_true": tmk}

        n_in = len(ins["R1"]) + len(ins["R2"])
        lev_res = {"reads_in": n_in}
        for tool in TOOLS:
            adapter_hits = untouched = dropped = false_clip = 0
            truth_ad_reads = truth_ad_bases = 0
            bases_in = bases_out = 0
            for mate in ("R1", "R2"):
                out = outs[tool][mate]
                for rid, ilen in ins[mate].items():
                    has_ad, ad_bases, _ = truth[(rid, mate)]
                    bases_in += ilen
                    if rid not in out:
                        dropped += 1
                        continue
                    olen = out[rid]
                    bases_out += olen
                    hit = olen < ilen
                    if hit:
                        if not has_ad:
                            false_clip += 1
                    else:
                        untouched += 1
                    if has_ad:
                        truth_ad_reads += 1
                        truth_ad_bases += ad_bases
                        if hit:
                            adapter_hits += 1
                        key = next(b for b in BINS if b[0] <= ad_bases <= b[1])
                        bt = bin_tot[key]
                        if tool in ("cutadapt", "trimmomatic_true"):
                            bt[tool] = bt.get(tool, 0) + (1 if hit else 0)
                            bt[tool + "_truth"] = bt.get(tool + "_truth", 0) + 1
                            by_remnant[tool][ad_bases].append(ilen - olen)
            bases_removed = bases_in - bases_out
            recall = 100.0 * adapter_hits / truth_ad_reads if truth_ad_reads else 0.0
            lev_res[tool] = {
                "reads_out": len(outs[tool]["R1"]) + len(outs[tool]["R2"]),
                "reads_dropped": dropped,
                "truth_adapter_reads": truth_ad_reads,
                "adapter_hit_reads": adapter_hits,
                "recall_pct": round(recall, 2),
                "untouched_reads": untouched,
                "false_clip_reads": false_clip,
                "false_clip_pct_of_clean_reads":
                    round(100.0 * false_clip / (n_in - truth_ad_reads), 2),
                "truth_adapter_bases": truth_ad_bases,
                "bases_removed": bases_removed,
                "bases_removed_pct": round(100.0 * bases_removed / bases_in, 4),
                "bases_removed_per_read": round(bases_removed / n_in, 3),
                "bases_in": bases_in,
            }
            tsv_rows.append(
                "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.2f\t%d\t%d\t%d\t%.4f\t%.3f"
                % (lev, tool, n_in, lev_res[tool]["reads_out"], dropped,
                   truth_ad_reads, adapter_hits, recall, false_clip, truth_ad_bases,
                   bases_removed, lev_res[tool]["bases_removed_pct"],
                   lev_res[tool]["bases_removed_per_read"]))
        results["levels"][lev] = lev_res

    for b in BINS:
        key = "%d-%d" % b
        bt = bin_tot[b]
        n_truth = bt["cutadapt_truth"]          # 与 trimmomatic_true_truth 相同
        ca_hit, tm_hit = bt["cutadapt"], bt["trimmomatic_true"]
        bin_rows.append("%s\t%d\t%d\t%.2f\t%d\t%.2f"
                        % (key, n_truth, ca_hit, 100.0 * ca_hit / n_truth,
                           tm_hit, 100.0 * tm_hit / n_truth))
        results["bins"][key] = {
            "n_truth": n_truth,
            "cutadapt_hit": ca_hit,
            "cutadapt_recall_pct": round(100.0 * ca_hit / n_truth, 2),
            "trimmomatic_hit": tm_hit,
            "trimmomatic_recall_pct": round(100.0 * tm_hit / n_truth, 2),
        }

    # cutadapt 自带报告的官方口径数字（5.2 报告输出在 stdout）
    for lev in LEVELS:
        rep = {}
        with open(os.path.join(BASE, "cutadapt_report_%s.txt" % lev),
                  encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if s.startswith("Total read pairs"):
                    rep["total_read_pairs"] = s.split()[4]
                elif s.startswith("Read 1 with adapter"):
                    rep["read1_with_adapter"] = {"pct": s.split("(")[1].split("%")[0].strip(),
                                                 "count": s.split()[4]}
                elif s.startswith("Read 2 with adapter"):
                    rep["read2_with_adapter"] = {"pct": s.split("(")[1].split("%")[0].strip(),
                                                 "count": s.split()[4]}
                elif s.startswith("Pairs written"):
                    rep["pairs_written"] = {"count": s.split()[4],
                                            "pct": s.split()[5].strip("()%")}
                elif s.startswith("Bases removed from R1"):
                    rep["bases_removed_R1"] = s.split()[4]
                elif s.startswith("Bases removed from R2"):
                    rep["bases_removed_R2"] = s.split()[4]
        results["cutadapt_report"][lev] = rep

    # 逐残留长度的平均实际剪除量（合并三档；cutadapt vs trimmomatic_true）
    results["by_remnant"] = {}
    for tool in ("cutadapt", "trimmomatic_true"):
        results["by_remnant"][tool] = {
            str(ab): {"n": len(v), "mean_removed": round(sum(v) / len(v), 2)}
            for ab, v in sorted(by_remnant[tool].items())
        }

    with open(os.path.join(BASE, "results.json"), "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)
    for name, rows in (("results.tsv", tsv_rows), ("bins.tsv", bin_rows)):
        with open(os.path.join(BASE, name), "w", encoding="utf-8") as f:
            f.write("\n".join(rows) + "\n")
    print("== per-level ==")
    print(json.dumps(results["levels"], indent=2))
    print("== bins ==")
    print("\n".join(bin_rows))
    print("== cutadapt report ==")
    print(json.dumps(results["cutadapt_report"], indent=2))
    print("wrote results.json, results.tsv, bins.tsv")


if __name__ == "__main__":
    main()
