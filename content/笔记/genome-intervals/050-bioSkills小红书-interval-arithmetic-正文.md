---
title: "bedtools：算术没错，错在前提"
skill: interval-arithmetic
trial: "050"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["bedtools", "BED", "基因组区间", "intersect", "bioSkills", "生信"]
date: "2026-09-03"
---

<!-- META
标题建议: bedtools：算术没错，错在前提
封面卡：150 条合成峰 vs 40 条基因，intersect 七种模式 + merge/subtract/map 逐个试，再跟一份独立 python 区间代数对账 40 项
/META -->

封面卡：150 条合成峰 vs 40 条基因，intersect 七种模式 + merge/subtract/map 逐个试，再跟一份独立 python 区间代数对账 40 项。工具是 bedtools v2.31.1，数据是 2 条模拟染色体共 180000 bp。

发现一：输出模式只改打印，不改重叠。同一对文件，-u 保留 73 条、-v 剔除 77 条、-wa -wb 打 84 行配对、重叠合计 49087 bp，七种模式与独立实现的期望值逐项相等（对账 40 项 vs 40 项全对）。一个细节：-loj 不是每条 A 一行，而是每个重叠对一行、无命中的 A 补 NULL，所以行数 = 84 行 + 77 行 = 161 行，跟 -wao 一样。

发现二：算术精确，前提会静默出错。三个实测：① 未排序文件喂给 merge，v2.31.1 直接报错退出（exit 1），输出 0 块——先 sort 再 merge 是硬纪律；② 染色体名 chr1 → 1 后，intersect 返回 0 条、exit 0，只有 stderr 一条命名 WARNING，计数为 0 全靠对账才能发现；③ -sorted 扫掠省内存 10.6 倍（30 万 × 30 万条区间，103.1 MB → 9.8 MB），结果与内存版完全一致（73 条 = 73 条）。

发现三：两个口径差异最大的地方。-f 是 A 的比例、-F 是 B 的比例：默认 1 bp 口径 73 条，-f 0.5（A=peaks）剩 59 条，-f 0.5 -r 双向各半只剩 7 条，交换 -a/-b 角色后 6 条——阈值加在谁身上直接改变答案。另一个是 -split：25 条带内含子的 BED12 转录本 vs 80 条外显子，按整条包络算重叠 5590 bp，按外显子块算只剩 1519 bp，包络口径虚高 3.7 倍，命中条数也多 4 条（14 条 vs 10 条）。

结论：区间运算本身是精确的集合代数，bedtools、pybedtools 算的是同一套几何；出错的是四个前提——输入排序、染色体命名一致、-split、-g genome 文件。跑交集前先统一命名、先 sort、RNA 类数据加 -split，然后用独立实现对账一遍，四步做完，区间运算基本不会错。

#生信 #生物信息学 #bedtools #BED #基因组区间 #bioSkills
