---
title: "MSA 先清两样再分析"
skill: msa-parsing
trial: "008"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["Biopython", "AlignIO", "msa-parsing", "生信", "bioSkills"]
date: "2026-09-03"
---

封面卡：比对完别急着统计，先清空位列和重复序列

用 BioPython 的 AlignIO 真跑 msa-parsing SKILL.md 的解析与过滤函数，输入是 5 条 × 30 列的小 DNA MSA（多序列比对），把预处理要点记下来。

发现一：空位比例超 50% 的列是伪影，先删。本例 cols 20-24 有 3/5 序列带空位（gap=0.6），被标成 gappy 列，remove_gappy_columns 一下把 30 列裁到 25 列。这种列不做下游统计，否则拉低整列指标。

发现二：近重复序列会污染每列指标，用 Henikoff 权重降权。本例 s1 和 s2 只在 1 个位点不同，是近重复。henikoff_weights 给它们 0.175 和 0.170，唯一序列 s3/s4/s5 是 0.215-0.220，权重和正好 1.0。系统发育结构数据上每列统计都要乘这个权重，否则被多数 clade 主导。

发现三：一致序列低于阈值的列会出模糊字符。本例可变区没有残基达到 0.5，一致序列在中间一段全是 N（DNA）或 X（蛋白），不能当真实残基用。坐标映射用 numpy 向量化版（cumsum）是 O(1) 查找，空位列返回 -1，比逐字符遍历快得多。

结论：比对解析不是读进来就完事。gappy 列先移除、冗余序列用 Henikoff 权重校正、一致序列的 N/X 别当真。核心解析函数本机全部跑通，MI-APC、Neff、A2M/A3M 在 examples 里没逐行复现。

#生信 #生物信息学 #Biopython #MSA #bioSkills
