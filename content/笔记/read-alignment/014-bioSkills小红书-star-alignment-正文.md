---
title: "STAR 比对 RNA 别漏3点"
skill: read-alignment/star-alignment
trial: "014"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["STAR", "RNA-seq", "比对", "生信", "bioSkills"]
date: "2026-08-27"
---

RNA 比对一旦选错参数，下游计数和变异 call 会静默出错。STAR 是剪接感知比对器，以下三个点最需注意。

第一，连接数据库 overhang。索引构建要带 GTF，且 `--sjdbOverhang` 必须等于读长减 1（100 bp 读长就填 99）；小基因组（如质粒）还要把 `--genomeSAindexNbases` 调小，用默认 14 索引会坏。第二，MAPQ 255 坑 GATK。STAR 给唯一比对打 255，GATK 当作"质量未知"直接丢弃，RNA VCF 会变成空文件；比对时加 `--outSAMmapqUnique 60` 就能避免。第三，链特异性别猜。开 `--quantMode GeneCounts`，看 ReadsPerGene 的正链/反链两列：两列接近是未链特异性，反链主导就是常见 dUTP/TruSeq 方案，选错列计数约减半。

补充两个边界：STAR 不自动解压，`.gz` 输入必须配 `--readFilesCommand zcat`；GATK 要读组，用 `--outSAMattrRGline` 空格分隔填 ID/SM/PL/LB，否则缺 @RG。

STAR 把剪接连接、两趟法、MAPQ 和链特异性都打包好了，但参数不是默认就安全——上面三点核对过再跑，能少返工很多。
#生信 #生物信息学 #STAR #RNA比对 #bioSkills
