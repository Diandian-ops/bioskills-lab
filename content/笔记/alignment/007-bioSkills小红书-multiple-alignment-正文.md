---
title: "三款比对工具实测对比"
skill: multiple-alignment
trial: "007"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["MAFFT", "MUSCLE5", "ClustalOmega", "多序列比对", "bioSkills", "生信"]
date: "2026-09-03"
---

<!-- META
标题建议: 三款比对工具实测对比
封面卡：比对工具那么多，MAFFT / MUSCLE5 / ClustalOmega 到底差多少？拿真实序列跑一遍看数字
/META -->

封面卡：同一个蛋白集，MAFFT、MUSCLE5、ClustalOmega 三款比对工具结果差多少？拿 6 条真实人类蛋白酶序列真跑了一遍看数字。

发现一：三款工具的比对长度几乎一样。6 条人类 S1 丝氨酸蛋白酶（长度 247–304 残基）比对后，MAFFT L-INS-i 给 348 列、MUSCLE5 `-align` 给 349 列、ClustalOmega 给 342 列，整体空位都在 22%–24% 区间。差异只是个位数列，同一数据集下主流工具结论一致，不用纠结哪款绝对更准。

发现二：`mafft --auto` 在小数据上等价于 L-INS-i。SKILL.md 说 `--auto` 会在 200 条处把算法从迭代优化翻成单次渐进，但本例只有 6 条，没触发降级，产物和显式 L-INS-i 逐字节一致（都 348 列、平均一致性 46.95%）。所以别盲信 `--auto` 省事——大集合一定要显式写算法并记录版本。

发现三：运行时间差很多。MUSCLE5 最快（0.051 秒），ClustalOmega 最慢（1.245 秒），差约 24 倍。序列少看不出，上到几千条这条就决定你能不能跑完。另外 trypsin 三个旁系之间一致性约 85%，和 chymotrypsin/elastase 之间跌到 33%–44%——分歧越大，Loop 区比对越不可靠，下游结论要保守。

结论：日常小数据 MAFFT L-INS-i 或 MUSCLE5 都稳；大集合用 FFT-NS-2 / `-super5` 并显式指定算法；发表前记得做 GUIDANCE2 或 MUSCLE5 集成评估逐列置信度。本机三款工具全部跑通，T-Coffee、密码子感知、置信度评估等扩展路径未逐行复现，结论按 SKILL.md 陈述。

#生信 #生物信息学 #多序列比对 #MAFFT #MUSCLE5 #bioSkills
