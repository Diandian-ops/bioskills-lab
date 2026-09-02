---
title: "重复标记处理"
skill: duplicate-handling
trial: "035"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["samtools", "markdup", "rmdup", "duplicate", "生信", "bioSkills"]
date: "2026-09-02"
---

扩增子或 PCR 会有重复 read，要标记掉。用 bowtie2 真跑产物（6000 条 reads）走了一遍 markdup 标准流程。

第一，流程四步不能省。sort -n 把成对 read 排一起，fixmate -m 写 MC 和 ms tag，再 sort 回坐标序，最后 markdup 标记重复；fixmate 跳过的话 markdup 不会判重。

第二，标记前 0 条、标记后 10 条。输入本身 FLAG 里没有 duplicate 位（0 条），markdup 后标记出 10 条疑似 PCR 或光学重复，增量从 0 条到 10 条。

第三，标记和删除两回事。markdup 默认只设 1 个重复标记位（0x400），保留记录，下游用 -F 1024 过滤；要直接删加 -r。旧版 rmdup 不依赖 fixmate tag，判重不如 markdup 准，新流程别用。

#生信 #生物信息学 #samtools #BAM #bioSkills
