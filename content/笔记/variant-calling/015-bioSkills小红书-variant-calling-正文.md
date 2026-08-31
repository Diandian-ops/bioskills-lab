---
title: "bcftools 变异检出跑通，三个关键参数"
skill: variant-calling
trial: "015"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["bcftools", "variant-calling", "SNP", "生信", "bioSkills"]
date: "2026-08-28"
---

比对做完下一步就是变异检出。bcftools 的 `mpileup | call` 两步法是最快的胚系 SNP/indel 检出方案之一，但参数选不对要么漏检要么假阳性。用 samtools 官方测试数据真跑了一遍，三个点最值得注意。

第一，质量过滤 `-q` 和 `-Q`。`-q 20` 过掉 MAPQ 低于 20 的读长（比对不可靠），`-Q 20` 过掉碱基质量低于 20 的碱基（测序噪声）。这次用的测试 BAM 读长 MAPQ 大部分在 60，两个过滤基本没影响——但如果你的数据覆盖度低或比对质量差，这两个参数就是控制假阳性的主要旋钮。第二，必须带注释 `-a FORMAT/DP,FORMAT/AD`。DP 是总深度，AD 是 REF 和 ALT 各自的等位基因深度，下游算 VAF、判断杂合纯合全靠它，不加的话 VCF 基本没法做质控过滤。第三，引擎别乱选。bcftools 是逐位点独立判定的位置模型，速度快不需要训练数据，但 indel 在同源重复和低复杂度区域精度明显弱于 GATK HaplotypeCaller 和 DeepVariant——后者会丢弃原始比对重新构建单倍型。简单场景用 bcftools 足够，人类 WGS 要高精度 indel 就换引擎。

补充两个边界：BAQ 默认开启（抑制 indel 附近假阳性 SNP，`-B` 才关闭）；不同来源的 VCF 在合并或比较前必须跑 `bcftools norm -f ref.fa` 做左对齐规范化，否则同一变异因表示不一致被当成不同记录。

一行命令就能出结果，但理解每个参数在做什么，才能在 0 变异的时候快速定位是数据问题还是参数问题。
#生信 #生物信息学 #变异检出 #bcftools #bioSkills
