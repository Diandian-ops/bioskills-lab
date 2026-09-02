---
title: "引物剪切处理"
skill: amplicon-clipping
trial: "037"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["samtools", "ampliconclip", "soft-clip", "primer", "生信", "bioSkills"]
date: "2026-09-02"
---

扩增子测序的引物区段要剪掉，免得干扰变异检出。用 bowtie2 真跑产物（6000 条 reads）跑了 samtools ampliconclip。

第一，primer BED 是最小真实示例。本例用 contig1 上的 2 个 primer 区间（坐标 100 到 200、250 到 350），坐标来自 BAM 实际覆盖区域，不是某用户的完整 panel。

第二，剪切后 soft-clip 从 0 条变 53 条。ampliconclip 输出 FILTERED 0、FAILED 0、WRITTEN 6000 条，全部写出；比对上 read 从 5954 条变 5948 条，被剪的碱基不再计入比对。

第三，剪完要重排序。clip 改了 read 的比对起点，输出不再坐标排序，直接 samtools index 会报 unsorted 失败，得先 sort 再 index。

诚实一句：primer BED 是构造的最小示例，真实 panel 规模更大，clip 比例会更高。

#生信 #生物信息学 #samtools #BAM #bioSkills
