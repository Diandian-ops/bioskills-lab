---
title: "MSA 统计先看这三样"
skill: msa-statistics
trial: "005"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["Biopython", "AlignIO", "msa-statistics", "生信", "bioSkills"]
date: "2026-09-03"
---

封面卡：比对完了怎么判断好坏？保守度、空位、一致性三张图

用 BioPython 的 AlignIO 真跑了一遍 msa-statistics SKILL.md 的统计函数，输入是一个 6 条 × 40 列的小蛋白 MSA，把关键指标记下来。

发现一：保守度曲线和熵曲线是镜像的。保守核心列保守度 100%、熵接近 0 bits；可变区保守度掉到约 17%、熵升到 2.5 bits 左右。蛋白的信息量别用均匀熵，要用 KL 散度对 Robinson & Robinson 1991 经验背景算——本例平均熵 1.03 bits，但 KL 信息量平均 3.34 bits，差出三倍，因为氨基酸频率本来就不均匀。

发现二：空位比例超过 50% 的列是伪影信号。本例 cols 28-31 有 4/6 序列带空位（gap=0.67），被直接标记为 gappy。这种列别拿去做下游统计，先移除或掩码。平均空位比例只有 6.7%，但集中的那几列才是问题。

发现三：一致性有四个分母（PID1-4），能差 11.5%。本例无空位等长对四个都是 60% 重合，但只要出现内部空位或长度不对称就会分开，写结论一定标清用哪一种。sum-of-pairs 用 BLOSUM62 矩阵按 matrix[c1,c2] 访问，遇到 U、J 这种字母表外字符会抛 IndexError，要 try/except 跳过。

结论：比对质量不是看一张总表，而是逐列看保守度、空位、一致性。gappy 列先处理，PID 定义要标注，蛋白信息量用 KL 背景。核心统计函数本机全部跑通，Capra-Singh JSD、PSSM、Neff、MI-APC 在 examples 里没逐行复现。

#生信 #生物信息学 #Biopython #MSA #bioSkills
