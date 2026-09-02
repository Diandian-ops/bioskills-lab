---
title: "BAM 三要素"
skill: sam-bam-basics
trial: "028"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["samtools", "BAM", "SAM", "CRAM", "生信", "bioSkills"]
date: "2026-09-01"
---

BAM 文件不是黑盒，看懂三个字段就够了：FLAG、MAPQ、CIGAR。用 bowtie2 真跑出来的 paired-end BAM（6000 reads）实测了一遍，把这三个字段的实测结果记下来。

第一，FLAG 是位掩码，不要硬背，用 samtools flags 解码。比如 99 = PAIRED + PROPER_PAIR + MREVERSE + READ1，147 = PAIRED + PROPER_PAIR + REVERSE + READ2。实测 6000 条 reads 里，5954 条 mapped，46 条 unmapped，3000 条 read1 和 3000 条 read2，2976 条在反向链。Secondary（0x100）和 Supplementary（0x800）在本数据里都是 0，但含义完全不同：secondary 是同一 read 的备选位置，SNV/indel 常用 -F 256 去掉；supplementary 是拆分比对的一部分，做 SV 或融合检测时必须保留。

第二，MAPQ 不是通用质量分，不同 aligner 刻度不一样。bowtie2 的 MAPQ 最高只到 42，这 6000 条 reads 里有 5398 条集中在 42。所以用 -q 60 过滤 bowtie2 的 BAM 会把所有结果都删掉。STAR 的“唯一比对” sentinel 是 255，BWA 才是 60，不能直接套同一个阈值。

第三，CIGAR 描述比对动作。实测 5815 条 read 是 clean 100M，139 条含 insertion（共 1014 bp），没有 deletion、soft-clip、hard-clip、skip。M 是 match/mismatch 的并集，N 是 RNA-seq 的 intron skip，算覆盖时必须排除，不能当成真实覆盖区域。

打开 BAM 先看 header：samtools view -H 能拿到 VN/SO、contig 名、@PG 程序链。格式转换用 -b（BAM）、-C（CRAM）、-h（带 header）；CRAM 必须配合同一参考基因组，用 -T reference.fa。区域查询是 1-based closed，和 pysam 的 0-based half-open 不一样，混用会在边界处出错。

#生信 #生物信息学 #BAM #samtools #bioSkills
