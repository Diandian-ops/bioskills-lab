---
title: "过滤与子集抽取"
skill: alignment-filtering
trial: "031"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["samtools", "BAM", "filter", "MAPQ", "FLAG", "生信", "bioSkills"]
date: "2026-09-02"
---

BAM 里不是每条都要留。用 bowtie2 真跑产物（6000 条 reads）实测了 samtools view 的过滤组合。

第一，-F 4 取比对上、-f 4 取没比对上。实测比对上 5954 条、没比对上 46 条；-f 2 取 proper-pair 得到 3042 条。

第二，MAPQ 阈值按比对器定。-q 20 取 MAPQ 不低于 20 的得到 5924 条。bowtie2 的 MAPQ 上限是 42，所以 -q 60 会把全部结果删光（60 是 BWA 上限，约为本输入上限的 1.4 倍）；STAR 唯一比对的标记是 255 不是 60，阈值不能混用。

第三，组合过滤看交集。-F 4 -q 20（比对上且 MAPQ≥20）也是 5924 条，和 -q 20 相同，因为那 46 条没比对上的 MAPQ 都是 0，已经被 -q 20 排掉了。-f 4 -F 8 取「没比对上但 mate 比对上了」的 read 有 12 条。

#生信 #生物信息学 #samtools #BAM #bioSkills
