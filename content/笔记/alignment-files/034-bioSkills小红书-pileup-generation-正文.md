---
title: "Pileup 深度统计"
skill: pileup-generation
trial: "034"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["samtools", "mpileup", "pileup", "depth", "生信", "bioSkills"]
date: "2026-09-02"
---

想看每个位置被多少条 read 盖住，用 samtools mpileup。用 bowtie2 真跑产物（6000 条 reads）实测。

第一，文本 pileup 能出。samtools mpileup 输出逐位点六列：contig、位置、参考碱基、深度、碱基堆叠、质量。前几行深度从 1 倍涨到 3 倍，说明 read 在该区域逐步叠加。

第二，深度有均值和峰值。文本 pileup 共 11962 个位点，平均深度 25.414 倍，最大深度 54 倍，覆盖相对均匀、没有极端堆积。

第三，BCF 这条路在本机断了。较新版本 samtools 移除了 mpileup -g 生成 BCF 的功能，本机也没装 bcftools，所以只产出文本 pileup，BCF 文件是 0 个字节。要出 BCF/VCF 得另装 bcftools mpileup。

#生信 #生物信息学 #samtools #BAM #bioSkills
