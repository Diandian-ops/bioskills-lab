---
title: "人为混进 10% 污染，筛查找得到吗"
skill: contamination-screening
trial: "043"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["污染筛查", "FastQScreen", "kraken2", "read-qc", "bioSkills", "生信"]
date: "2026-09-03"
---

<!-- META
标题建议: 人为混进 10% 污染，筛查找得到吗
封面卡：3 份各 50000 对 reads、按 90/5/5 等比例混合三种真实基因组，kraken2 与 FastQ Screen 各查一遍，设计比例 vs 检出比例逐项对账
/META -->

封面卡：3 份各 50000 对 reads、按 90/5/5 与 75/15/10 的比例混合三种真实基因组，kraken2 与 FastQ Screen 各查一遍，设计比例 vs 检出比例逐项对账。

发现一：设计比例被原样检出。三个基因组各按已知比例混合（S1 = 90% E. coli + 5% PhiX + 5% lambda，S2 = 75% + 15% + 10%），kraken2（confidence 0.1）的物种级检出与设计最大偏差 0.15 个百分点——S2 的 lambda 设计 10%、检出 9.85%；阴性对照（100% E. coli）检出 99.98%，零假阳性物种。

发现二：PhiX 精确、lambda 差一截，差值有去处。FastQ Screen 的唯一命中列（One_hit_one_genome）对 PhiX 给出与设计完全相等的 5.00% 和 15.00%；lambda 只给出 4.48% 和 8.99%，少的 0.5–1 个百分点全部落进「同时命中多个基因组」类别——同源区 reads 被分到多重命中而非丢失，这正是要看唯一命中列、不能只看总比对率的原因。

发现三：提高 confidence，代价由弱信号承担。把 confidence 从 0.0 提到 0.2，PhiX 检出纹丝不动（5.00% / 15.00%），lambda 弱污染从 4.98% 修剪到 4.82%（S1）。压掉痕迹物种长尾的同时，最先牺牲的是弱污染信号——剂量要按样本类型权衡。

结论：物种筛查能可靠回答「有哪些物种」，且比例可定量复核；但它回答不了「是谁的 DNA」——同种样本交换与 index hopping 在物种筛查里呈现完全干净的单物种轮廓，人类或单物种队列必须补 SNP 指纹检查（somalier / NGSCheckMate）。本次用 3 个 NCBI 真实基因组自建 13 MB 迷你库即可完成整套定量基准；自建库时分类树 dmp 文件必须用标准「TAB+竖线+TAB」分隔，否则库建得完、样本全未分类。

#生信 #生物信息学 #kraken2 #FastQScreen #污染筛查 #bioSkills（工具版本=kraken2 2.17.1 / FastQ Screen 0.16.0）
