---
title: "索引与区域查询"
skill: alignment-indexing
trial: "030"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["samtools", "BAM", "index", "bai", "csi", "生信", "bioSkills"]
date: "2026-09-02"
---

BAM 不建索引就没法随机查区域。用 bowtie2 真跑产物（6000 条 reads）实测了 samtools index 的两种索引。

第一，BAI 和 CSI 都建。samtools index 默认出 .bai（占 336 个字节），加 -c 出 .csi（占 155 个字节）。本输入 4 个 contig 都只有 3000 bp，两种都够用；超长 contig 必须用 CSI。

第二，每 contig 比对数看 idxstats。四个 contig 比对上的 reads 是 1485、1491、1490、1488 条，末行星号汇总了 34 条无法定位的 unmapped。

第三，区域查询计数含 unmapped。samtools view -c contig1 返回 1486 条，比 idxstats 的 1485 多 1 条，多出的就是落在 contig1 上但本身没比上的那条 read。别把这点差异当成文件坏了。

#生信 #生物信息学 #samtools #BAM #bioSkills
