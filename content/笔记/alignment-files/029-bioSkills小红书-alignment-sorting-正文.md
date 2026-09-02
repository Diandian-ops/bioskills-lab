---
title: "排序顺序控制"
skill: alignment-sorting
trial: "029"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["samtools", "BAM", "sort", "alignment", "生信", "bioSkills"]
date: "2026-09-02"
---

比对文件不是随便排的，顺序写在 header 里。用 bowtie2 真跑产物（6000 条 reads）实测了一遍 samtools sort 的两种模式。

第一，坐标排序是默认。samtools sort 不带头参时按参考坐标排，@HD 写 SO:coordinate。本输入本就坐标排序，再排一次 @HD 还是 SO:coordinate，记录数保持 6000 条进、6000 条出，排序不改内容只改顺序。

第二，名称排序靠 -n。samtools sort -n 按读段名排，@HD 变成 SO:queryname、SS:queryname:natural。排完前 5 个位置从 3、4、6、17、21 变成 327、3、4、388、6，说明文件内顺序确实按名字重排了。

第三，下游决定用哪种。建索引、pileup 要坐标排序；fixmate 和 markdup 要名称排序（成对 read 相邻）。排序前先看 @HD 确认 SO，避免拿错顺序的 BAM 去建索引。

#生信 #生物信息学 #samtools #BAM #bioSkills
