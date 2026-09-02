---
title: "双序列比对怎么选模式"
skill: pairwise-alignment
trial: "004"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["Biopython", "PairwiseAligner", "双序列比对", "生信", "bioSkills"]
date: "2026-09-03"
---

封面卡：两条序列怎么比对？全局/局部/半全局，分数差一倍

用 Biopython 的 PairwiseAligner 真跑了一遍 SKILL.md 的三种模式，把实测结果记下来。

发现一：模式选错，分数天差地别。两条「两端发散、中间保守」的序列，全局比对 score 只有 12.0，局部比对是 20.0，半全局是 11.0。局部比对会跳过发散的两端、只给保守核心打分，所以最高。长度差很多或者两端不保守的序列，别用全局，改局部或半全局（自由端空位）。

发现二：用替换矩阵必须显式设空位罚分。PairwiseAligner 默认的 open/extend 都是 0，配上 BLOSUM62 这种正分矩阵，会冒出一堆无意义的短空位。蛋白比对要用 open=-11、extend=-1（BLASTP 默认）。实测 protA/protB 用 BLOSUM62 比对，score=41.0，8 个一致、1 个错配。

发现三：百分比一致性有四个分母（PID1-4），同一个比对能差到 11.5%。本报告的无空位等长例子四个都是 84.6% 重合，但只要出现内部空位或长度不对称就会分开，下游比较一定写清用哪一种。另外 len(alignment) 返回的是序列条数不是比对长度，比对长度用 alignment.shape[1]。

结论：PairwiseAligner 默认适合小于 10 kb 的交互和脚本。空位罚分记得显式给；全球/局部/半全局按序列特征选；PID 定义要标注。比对器核心流程在本机全部跑通，库性能基准和 Karlin-Altschul 显著性未实测。

#生信 #生物信息学 #Biopython #双序列比对 #bioSkills
