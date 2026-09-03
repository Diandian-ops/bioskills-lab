---
title: "同坐标的两条 reads 是重复吗"
skill: umi-processing
trial: "042"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["UMI", "umi_tools", "去重", "PCR重复", "bioSkills", "生信"]
date: "2026-09-03"
---

<!-- META
标题建议: 同坐标的两条 reads 是重复吗
封面卡：PCR 重复数错会怎样？6000 个模拟分子、16445 对 reads，四种去重方法对同一份 BAM 各数一遍，再和 truth 对账
/META -->

封面卡：PCR 重复数错会怎样？6000 个模拟分子、16445 对 reads（每个分子带 12 nt UMI），四种去重方法对同一份 BAM 各数一遍，再和 truth 对账。

发现一：只看坐标会错杀。模拟数据里故意放了 300 对共享同一坐标、但 UMI 不同的分子对，samtools markdup 按坐标去重只数出 5313 个分子，比 truth 6000 少 11.45%——坐标相同不等于同一分子，UMI 才是分子的身份证。

发现二：精确匹配 UMI 会多计。869 份 PCR 拷贝的 UMI 带着 1 个碱基的测序错误，unique 方法（无错误模型）数出 7036 个「分子」，比 truth 多 17.27%。编辑距离分布给出铁证：同位置 UMI 对在距离 1 处观测到 674 对，而随机 UMI 的零假设是 0 对——这 674 对全是错误抄本，不是新分子。

发现三：directional（umi_tools 默认）最接近 truth：6013 个，偏差 +0.22%，读对保留率 36.56%。它是唯一带错误模型的方法，把 1-off 邻居按计数梯度（父本 ≥ 子代约 2 倍）折回母本。cluster 方法单连通全并，数出 5891 个（少 1.82%），150 对同坐标干扰对里只保住 1 对，合并得激进。

结论：有 UMI 的文库，去重键是（坐标 + UMI）、方法用默认 directional；坐标去重留给非 UMI 的 DNA 文库。注意两个路由禁区：非 UMI bulk RNA-seq 的重复坐标是生物学事实，不做 dedup；CellRanger/STARsolo 的输出已自带 UMI 折叠，不要再处理一遍。安装走 conda（bioconda，python 3.11 环境）；pip 装在 python 3.13 下构建失败，且 `umi_tools version` 子命令在 1.1.6 是损坏的，版本用 conda list 核对。

#生信 #生物信息学 #UMI #umi_tools #PCR #bioSkills
