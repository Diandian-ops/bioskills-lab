---
title: "统计与覆盖汇总"
skill: bam-statistics
trial: "033"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["samtools", "BAM", "flagstat", "coverage", "bedcov", "生信", "bioSkills"]
date: "2026-09-02"
---

比对完要看数字画像。用 bowtie2 真跑产物（6000 条 reads）跑了 flagstat + stats + coverage + bedcov。

第一，比对率看 flagstat。5954 条比对上（99.23%），proper-pair 3042 条（50.70%）；插入片段平均 498.0、标准差 48.5，符合 paired-end 预期。

第二，覆盖看 coverage。四个 contig 覆盖碱基数 2997、2995、2995、2993 个，覆盖率 99.77% 到 99.9%，平均深度约 49.4 到 49.6 倍，覆盖很均匀。

第三，区间覆盖看 bedcov。四个 contig 全区间累计覆盖量 148260、148837、148712、148577 个碱基，合计 594386 个，和总比对碱基数一致。

covbases 是覆盖碱基数、meandepth 是平均深度，两个不是一个意思，别混。

#生信 #生物信息学 #samtools #BAM #bioSkills
