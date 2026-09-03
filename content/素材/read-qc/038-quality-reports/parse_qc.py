#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parse_qc.py — 解析本次真跑的 FastQC zip（fastqc_data.txt）与 MultiQC 数据目录，
抽取真实数值，落盘 qc_summary.json（逐位点曲线 + 模块判定）与 qc_summary.tsv（每样本一行）。

只解析真实文件中的真实数字，不做任何插值或虚构。
"""
import glob
import json
import os
import re
import zipfile

BASE = os.path.dirname(os.path.abspath(__file__))
QC_RAW = os.path.join(BASE, "qc", "raw")
MQC_DIR = os.path.join(BASE, "qc", "multiqc", "multiqc_data")

TABLE_MODULES = {
    "Basic Statistics": "basic",
    "Per base sequence quality": "per_base_quality",
    "Per base sequence content": "per_base_content",
    "Sequence duplication levels": "duplication_levels",
    "Overrepresented sequences": "overrepresented",
    "Adapter Content": "adapter_content",
}


def parse_fastqc_data(text):
    out = {"modules": {}, "statuses": {}}
    m = re.search(r"#Total Deduplicated Percentage\t([\d.]+)", text)
    if m:
        out["dedup_pct"] = float(m.group(1))
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith(">>") and not ln.startswith(">>END_MODULE"):
            name, status = ln[2:].split("\t")
            out["statuses"][name] = status
            i += 1
            rows = []
            header = None
            while i < len(lines) and not lines[i].startswith(">>END_MODULE"):
                row = lines[i]
                if header is None:
                    header = row.split("\t")
                else:
                    rows.append(row.split("\t"))
                i += 1
            key = TABLE_MODULES.get(name)
            if key and header:
                cols = header[1:]
                if name == "Basic Statistics":
                    # 两列 key/value 表（Variable | Value），按键值对存
                    kv = {r[0]: (r[1] if len(r) > 1 else "") for r in rows}
                    out["modules"][key] = {"name": name, "kv": kv}
                else:
                    cols = [c.lstrip("%").lstrip("#") for c in header]
                    table = {c: [] for c in cols}
                    for r in rows:
                        for c, v in zip(cols, r):
                            if c in table:
                                table[c].append(v)
                    out["modules"][key] = {"name": name, "columns": cols, "rows": table}
            if name == "Overrepresented sequences":
                # 表头: Sequence Count Percentage Possible Source
                recs = []
                for r in rows:
                    if len(r) >= 4:
                        recs.append({"sequence": r[0], "count": int(r[1]),
                                     "percent": float(r[2]), "source": r[3]})
                out["modules"]["overrepresented"] = {"name": name, "records": recs}
        i += 1
    return out


def to_float_list(xs):
    out = []
    for x in xs:
        try:
            out.append(float(x))
        except ValueError:
            out.append(None)
    return out


def main():
    summary = {"samples": {}, "multiqc": {}}
    for zp in sorted(glob.glob(os.path.join(QC_RAW, "*_fastqc.zip"))):
        with zipfile.ZipFile(zp) as z:
            data_name = [n for n in z.namelist() if n.endswith("fastqc_data.txt")][0]
            text = z.read(data_name).decode("utf-8")
        d = parse_fastqc_data(text)
        sname = os.path.basename(zp).replace("_fastqc.zip", "")

        basic = d["modules"]["basic"]["kv"]
        total = int(basic["Total Sequences"])
        gc = float(basic["%GC"])
        # Total Deduplicated Percentage 在 duplication 模块头部（Basic Statistics 无此项）
        dedup = d.get("dedup_pct")

        pbq = d["modules"]["per_base_quality"]["rows"]
        # FastQC 将位置分组（1..9 单碱基，之后成对如 10-11）；取组末端作为真实 x
        pbq_pos = [int(p.split("-")[-1]) for p in pbq["Base"]]
        mean_q = to_float_list(pbq["Mean"])
        med_q = to_float_list(pbq["Median"])
        lq = to_float_list(pbq["Lower Quartile"])

        pbc = d["modules"]["per_base_content"]["rows"]
        pbc_pos = [int(p.split("-")[-1]) for p in pbc["Base"]]
        content = {b: to_float_list(pbc[b]) for b in ("A", "C", "G", "T")}

        adp = d["modules"]["adapter_content"]["rows"]
        # 取各 adapter 列随位置的百分比，汇总为该位置全部 adapter 合计
        adp_pos = [int(p.split("-")[-1]) for p in adp["Position"]]
        adp_cols = [c for c in adp.keys() if c != "Position"]
        adp_total = [0.0] * len(adp["Position"])
        for c in adp_cols:
            vals = to_float_list(adp[c])
            for j, v in enumerate(vals):
                if v is not None:
                    adp_total[j] += v
        max_adapter = max(adp_total)

        over = d["modules"].get("overrepresented", {}).get("records", [])
        over_max_pct = max([r["percent"] for r in over], default=0.0)
        over_total_reads = sum(r["count"] for r in over)

        summary["samples"][sname] = {
            "total_sequences": total,
            "percent_gc": gc,
            "dedup_percentage": dedup,
            "mean_quality_first10": round(sum(mean_q[:10]) / 10, 2),
            "mean_quality_last10": round(sum(mean_q[-10:]) / 10, 2),
            "min_mean_quality": min(q for q in mean_q if q is not None),
            "per_base_quality_positions": pbq_pos,
            "per_base_quality_mean": mean_q,
            "per_base_quality_median": med_q,
            "per_base_quality_lower_quartile": lq,
            "per_base_content_positions": pbc_pos,
            "per_base_content": content,
            "adapter_total_by_position": [round(v, 2) for v in adp_total],
            "adapter_positions": adp_pos,
            "max_adapter_percent": round(max_adapter, 2),
            "overrepresented_count": len(over),
            "overrepresented_max_percent": over_max_pct,
            "overrepresented_total_reads": over_total_reads,
            "top_overrepresented": over[0] if over else None,
            "module_statuses": d["statuses"],
        }

    # MultiQC general stats（真实抓取值，与 fastqc_data.txt 交叉核对）
    gs_path = os.path.join(MQC_DIR, "multiqc_general_stats.txt")
    if os.path.exists(gs_path):
        with open(gs_path, encoding="utf-8") as f:
            header = f.readline().rstrip("\n").split("\t")
            for ln in f:
                parts = ln.rstrip("\n").split("\t")
                summary["multiqc"][parts[0]] = dict(zip(header[1:], parts[1:]))

    with open(os.path.join(BASE, "qc_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    # TSV：每样本一行
    tsv_cols = ["sample", "total_sequences", "percent_gc", "dedup_percentage",
                "mean_quality_first10", "mean_quality_last10", "min_mean_quality",
                "max_adapter_percent", "overrepresented_max_percent",
                "n_status_pass", "n_status_warn", "n_status_fail"]
    with open(os.path.join(BASE, "qc_summary.tsv"), "w", encoding="utf-8") as f:
        f.write("\t".join(tsv_cols) + "\n")
        for sname, s in summary["samples"].items():
            st = s["module_statuses"]
            row = [sname, str(s["total_sequences"]), str(s["percent_gc"]),
                   str(s["dedup_percentage"]), str(s["mean_quality_first10"]),
                   str(s["mean_quality_last10"]), str(s["min_mean_quality"]),
                   str(s["max_adapter_percent"]), str(s["overrepresented_max_percent"]),
                   str(sum(1 for v in st.values() if v == "PASS")),
                   str(sum(1 for v in st.values() if v == "WARN")),
                   str(sum(1 for v in st.values() if v == "FAIL"))]
            f.write("\t".join(row) + "\n")

    print("parsed %d samples; multiqc rows: %d" % (len(summary["samples"]), len(summary["multiqc"])))
    for sname, s in summary["samples"].items():
        print("%s: %d reads, GC %.1f%%, dedup %s%%, Q(first10) %.2f, Q(last10) %.2f, "
              "adapter max %.2f%%, overrep max %.2f%%" % (
                  sname, s["total_sequences"], s["percent_gc"], s["dedup_percentage"],
                  s["mean_quality_first10"], s["mean_quality_last10"],
                  s["max_adapter_percent"], s["overrepresented_max_percent"]))
        for m, v in s["module_statuses"].items():
            print("    %-32s %s" % (m, v))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
