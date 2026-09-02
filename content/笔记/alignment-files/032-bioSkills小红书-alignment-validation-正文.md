---
title: "比对完整性校验"
skill: alignment-validation
trial: "032"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["samtools", "BAM", "quickcheck", "flagstat", "stats", "生信", "bioSkills"]
date: "2026-09-02"
---

BAM 交出去前得先确认没坏。本机 samtools 没有 validate 子命令，改用 quickcheck + flagstat + stats 实测了一遍（输入 6000 条 reads）。

第一，quickcheck 最快。samtools quickcheck 退出码 0，没有发现 1 个表层损坏，说明文件能打开、magic 正确、索引在。

第二，flagstat 看自洽。6000 条里比对上 5954 条（99.23%）、没比对上 46 条，两者相加正好 6000；proper-pair 3042 条（50.70%），singletons 12 条（0.20%）。

第三，stats 给细节。error rate 约 2.19%，插入片段平均 498.0、标准差 48.5，inward pairs 2971 条，和 flagstat 数值互证。

诚实一句：validate 在本机报 unrecognized command，逐记录的坐标与 CIGAR 校验没法做，上面这套是轻量替代。

#生信 #生物信息学 #samtools #BAM #bioSkills
