---
title: "参考索引与字典"
skill: reference-operations
trial: "036"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["samtools", "faidx", "dict", "GC", "生信", "bioSkills"]
date: "2026-09-02"
---

比对前参考序列要先建索引。用本目录的 reference.fa（4 个 contig，各 3000 bp）实测了 faidx 和区域提取。

第一，faidx 建索引很快。samtools faidx 产出 .fai，五列记录每个 contig 的名字、长度、偏移、行宽。四个 contig 长度都 3000 个碱基，偏移 39、3121、6203、9285。

第二，区域提取按 1-based。samtools faidx ref.fa contig1:1-60 返回 contig1 前 60 个碱基，区间写法和 samtools view 区域查询一致，和 BED 的 0-based 不同。

第三，GC 含量很均衡。四个 contig 的 GC 是 50.27%、50.0%、51.87%、50.9%，都在 50% 附近，没有明显偏倚。

诚实一句：本机 samtools dict 出不了含 @SQ 的完整字典（只返回 @HD 一行），GATK 要的字典得用 picard 或 gatk 生成；faidx 和区域提取都正常。

#生信 #生物信息学 #samtools #参考序列 #bioSkills
