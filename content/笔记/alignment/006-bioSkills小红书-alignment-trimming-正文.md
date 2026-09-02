---
title: "三款修剪工具实测对比"
skill: alignment-trimming
trial: "006"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["ClipKIT", "trimAl", "BMGE", "多序列比对", "bioSkills", "生信"]
date: "2026-09-03"
---

<!-- META
标题建议: 三款修剪工具实测对比
封面卡：比对建完要不要修剪？ClipKIT / trimAl / BMGE 拿同一个真实比对跑一遍，看各砍掉多少列
/META -->

封面卡：比对建完要不要修剪？ClipKIT、trimAl、BMGE 拿同一个真实比对（6 条人类蛋白酶，348 列）跑了一遍，看各砍掉多少列。

发现一：激进度差很多。同一份 348 列的比对，ClipKIT `kpic-smart-gap` 砍掉 223 列（剩 125，移除 64.1%），trimAl `-strictplus` 砍 158 列（剩 190，45.4%），而 ClipKIT `gappyout`、trimAl `-automated1`、BMGE 只砍 99–108 列（剩 240–249，约 29%–31%）。同一个输入，模式决定砍多少。

发现二：SKILL.md 的 20%/40% 规则是真的硬指标。移除超过 40% 列就偏激、会连信号一起删。本例 `kpic-smart-gap`（64.1%）和 `-strictplus`（45.4%）都超了——因为这 6 条是 trypsin 和 chymotrypsin/elastase 混在一起（跨旁系一致性只有 33%–44%），比对本身分歧大，推荐默认模式在这类分歧集上偏激进、过砍。日常近缘集用默认没问题，分歧集要换轻模式或干脆不修。

发现三：修剪落点都集中在大空位列。ClipKIT 的 `--log` 和 trimAl 的 `-colnumbering` 能给出保留列索引，被删的列几乎全在空位比例高的区域——这正说明修剪在按空位分数干活，不是乱删。下游要做逐位点分析记得留着这个列映射。

结论：HMM（隐马尔可夫模型）建库用 trimAl `-gappyout`，系统发育超级矩阵用 ClipKIT `kpic-smart-gap`，深原核用 BMGE `-h 0.4`；发表级复现别用 trimAl `-automated1`（跨版本行为漂移），显式写模式并记版本。三款工具本机全跑通，Divvier、HMMcleaner、Gblocks、TCS/MACSE 等扩展路径未逐行复现。

#生信 #生物信息学 #多序列比对 #ClipKIT #trimAl #bioSkills
