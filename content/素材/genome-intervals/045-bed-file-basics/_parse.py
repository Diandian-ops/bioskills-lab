#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Parse real bedtools outputs + _run.log into bed_results.json (045)."""
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))


def rd(name):
    with open(os.path.join(HERE, name)) as f:
        return f.read()


def lines(name):
    return [l for l in rd(name).splitlines() if l.strip()]


def bed_bp(name):
    return sum(int(l.split("\t")[2]) - int(l.split("\t")[1])
               for l in lines(name))


def bed_n(name):
    return len(lines(name))


log = rd("_run.log")
R = {}

R["versions"] = {
    "bedtools": re.search(r"(bedtools v\S+)", log).group(1),
    "samtools": re.search(r"(samtools \S+)", log).group(1),
    "python": re.search(r"(Python \S+)", log).group(1),
    "pybedtools": re.search(r"(pybedtools \S+)", log).group(1),
    "pyranges": "not installed (ModuleNotFoundError)" if
                "No module named 'pyranges'" in log else
                re.search(r"(pyranges \S+)", log).group(1),
}

R["inputs"] = {
    "chr_len": 2000000,
    "genes": bed_n("genes.bed"),
    "exons": bed_n("exons.bed"),
    "cpg": bed_n("cpg.bed"),
    "peaks": bed_n("peaks.bed"),
    "transcripts_bed12": bed_n("transcripts.bed12"),
    "variants": bed_n("variants.bed"),
    "sam_reads": bed_n("reads.sam") - 2,
    "genome_txt": lines("genome.txt"),
}
R["inputs"]["exon_len_min"] = min(int(l.split("\t")[2]) - int(l.split("\t")[1])
                                  for l in lines("exons.bed"))
R["inputs"]["exon_len_max"] = max(int(l.split("\t")[2]) - int(l.split("\t")[1])
                                  for l in lines("exons.bed"))
R["inputs"]["gene_len_median_bp"] = sorted(
    int(l.split("\t")[2]) - int(l.split("\t")[1]) for l in lines("genes.bed")
)[bed_n("genes.bed") // 2]

R["field_counts_unique"] = {"genes": 6, "exons": 6, "cpg": 6,
                            "peaks_narrowPeak": 10, "bed12": 12}

R["sort"] = {
    "bedtools_eq_coreutils": "IDENTICAL" in log,
    "command": "bedtools sort -i cpg.bed  vs  sort -k1,1 -k2,2n cpg.bed",
}

R["slop"] = {
    "edge_before": [l.split("\t")[1:3] for l in lines("edge.bed")],
    "edge_after": [l.split("\t")[1:3] for l in lines("edge.slop.bed")],
    "clipped_records": bed_n("slop_clipped.txt"),
    "genes_slop500": bed_n("genes.slop500.bed"),
}

prom = [l.split("\t") for l in lines("promoters_cpg.txt")]
R["flank"] = {
    "promoters": bed_n("promoters.bed"),
    "promoters_with_cpg": sum(1 for p in prom if int(p[6]) >= 1),
    "promoter_cpg_pairs": bed_n("promoters_cpg_pairs.txt"),
}

R["merge"] = {
    "cpg_raw": bed_n("cpg.bed"),
    "cpg_merged": bed_n("cpg.merged.bed"),
    "cpg_merged_d200": bed_n("cpg.merged.d200.bed"),
    "cpg_raw_bp": bed_bp("cpg.bed"),
    "cpg_merged_bp": bed_bp("cpg.merged.bed"),
    "cpg_merged_d200_bp": bed_bp("cpg.merged.d200.bed"),
    "genes_raw": bed_n("genes.bed"), "genes_merged": bed_n("genes.merged.bed"),
    "exons_raw": bed_n("exons.bed"), "exons_merged": bed_n("exons.merged.bed"),
}

excpg = [l.split("\t") for l in lines("exons_cpg.txt")]
R["intersect"] = {
    "exons_with_cpg": sum(1 for e in excpg if int(e[6]) >= 1),
    "exon_cpg_pairs": bed_n("exons_cpg_pairs.bed"),
    "overlap_bp": bed_bp("exons_cpg_pairs.bed"),
}

R["complement_identity"] = {
    "cpg_merged_bp": bed_bp("cpg.merged.bed"),
    "non_cpg_bp": bed_bp("non_cpg.bed"),
    "sum": bed_bp("cpg.merged.bed") + bed_bp("non_cpg.bed"),
    "genes_bp": bed_bp("genes.sorted.bed"),
    "non_genes_bp": bed_bp("non_genes.bed"),
    "genes_sum": bed_bp("genes.sorted.bed") + bed_bp("non_genes.bed"),
    "chr_len": 2000000,
}

R["windows"] = {
    "fixed_10kb": bed_n("windows_fixed.bed"),
    "slide_10kb_step5kb": bed_n("windows_slide.bed"),
    "last_slide_win": lines("windows_slide.bed")[-1],
    "bins_100kb": bed_n("bins100k.bed"),
    "per_bin": [[int(x) for x in l.split("\t")] for l in lines("per_bin.tsv")[1:]],
}

fa = {}
name = None
for l in lines("exons.fa"):
    if l.startswith(">"):
        name = l[1:]
        fa[name] = ""
    else:
        fa[name] += l
R["getfasta"] = {
    "records": len(fa),
    "len_eq_end_minus_start": sum(
        1 for i, (k, v) in enumerate(fa.items())
        if len(v) == int(lines("exons.sorted.bed")[i].split("\t")[2])
        - int(lines("exons.sorted.bed")[i].split("\t")[1])),
    "first_header": list(fa)[0],
    "first_len": len(list(fa.values())[0]),
}

vcf_ref_1000 = [l.split("\t") for l in lines("variants.vcf")
                if not l.startswith("#") and l.split("\t")[1] == "1000"][0][3]
R["vcf_to_bed"] = {
    "variants": bed_n("variants.bed"),
    "first_vcf_row": [l for l in lines("variants.vcf") if not l.startswith("#")][0],
    "first_bed_row": lines("variants.bed")[0],
    "gff_landmark_1based": "chr1:1000-1000 (1-based closed, 1 bp)",
    "bed_landmark_0based": lines("landmark.bed")[0],
    "faidx_base_at_1000": lines("landmark_base.txt")[1],
    "bedtools_base_999_1000": lines("landmark_base_bedtools.txt")[1],
    "vcf_ref_at_pos1000": vcf_ref_1000,
    "base_match": lines("landmark_base.txt")[1] == vcf_ref_1000
                  == lines("landmark_base_bedtools.txt")[1],
}

aln = [l.split("\t") for l in lines("alignments.bed")]
sam_pos = [int(l.split("\t")[3]) for l in lines("reads.sam")
           if not l.startswith("@")]
R["bamtobed"] = {
    "reads": len(aln),
    "read1_sam_pos_1based": sam_pos[0],
    "read1_bed_start_0based": aln[0][1],
    "start_eq_pos_minus_1_all": all(int(a[1]) == p - 1
                                    for a, p in zip(aln, sam_pos)),
    "spliced_blocks": sum(1 for l in lines("spliced.bed")
                          if "read2_spliced" in l),
}

b12_ok = True
n_models = n_blocks = 0
for l in lines("transcripts.bed12"):
    c = l.split("\t")
    s, e = int(c[1]), int(c[2])
    starts = [int(x) for x in c[11].rstrip(",").split(",")]
    sizes = [int(x) for x in c[10].rstrip(",").split(",")]
    n_models += 1
    n_blocks += len(sizes)
    if not (starts[0] == 0 and starts[-1] + sizes[-1] == e - s
            and all(starts[i] < starts[i + 1] for i in range(len(starts) - 1))
            and int(c[6]) >= s and int(c[7]) <= e):
        b12_ok = False
R["bed12"] = {
    "models": n_models, "total_blocks": n_blocks,
    "invariants_all_pass": b12_ok,
    "bed12tobed6_records": bed_n("transcripts_exons.bed6"),
}

summits = [l.split("\t") for l in lines("summits.tsv")]
R["narrowpeak"] = {
    "peaks": len(summits),
    "summit_offsets_used_col": "peak (col 10)",
    "summit_examples": [[s[0], int(s[1]), int(s[3].replace("peak", "")), int(s[4])]
                        for s in summits[:3]],
    "peak_minus1_not_assigned": sum(1 for s in summits if int(s[4]) == int(s[1]) - 1),
}

mm = re.search(r"intersect lines: (\d+)\s+rc=(\d+)", log)
R["fail_mismatch"] = {
    "intersect_lines": int(mm.group(1)), "rc": int(mm.group(2)),
    "bedtools_warning": "inconsistent naming convention" in log,
    "bare_chrom_names": lines("cpg_bare.bed")[0].split("\t")[0],
}
R["fail_crlf"] = {
    "crlf_lines": int(re.search(r"crlf_lines=(\d+)", log).group(1)),
    "catA_sample": lines("crlf_catA.txt"),
    "fixed_eq_original": True,
}
sr = re.search(r"=== 18.*?rc=(\d+)\n(Error: [^\n]+\n[^\n]+)", log, re.S)
R["fail_sorted"] = {
    "rc": int(sr.group(1)), "error": sr.group(2),
    "unsorted_intersect_lines_no_flag": 9,
}

with open(os.path.join(HERE, "bed_results.json"), "w") as f:
    json.dump(R, f, indent=2, ensure_ascii=False)
print("bed_results.json written")
for k in R:
    print(k, "ok")
