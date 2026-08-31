---
title: "VCF 里 QUAL 和 GQ 的区别"
skill: vcf-basics
trial: "016"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["bcftools", "cyvcf2", "VCF", "BCF", "生信", "bioSkills"]
date: "2026-08-28"
---

VCF 是最常见、也最容易读错的文件格式。同一个文件里 QUAL 和 GQ 都标着「质量」，管的却是两件事。拿 bcftools 真跑出来的单样本 VCF（9 个变异位点）把这套字段解读流程走了一遍，三个点值得记下来。

第一，QUAL 和 GQ 的层级不同。QUAL 是位点级的，回答「这个位点有没有变异」，数值会随总深度放大，所以高覆盖的假象位点也能拿到很高的 QUAL；GQ 是每样本基因型级的，回答「这个基因型判对了没」，上限压在 99。这 9 个位点里有 6 个的 GQ 顶到上限、失去了区分度，剩下 3 个才真正拉开差距，而且排序是冲突的：有一个位点 QUAL 排第 4 高、GQ 却是全组最低，另一个位点 QUAL 全组最低、GQ 却比它高出一大截。过滤位点级噪声看 QUAL（更准确的是按深度归一化的 QD），判定基因型能不能信看 GQ，两者换不了。

第二，sum(AD) 小于 DP 是正常的。DP 计入跨越这个位点的全部读长，包括碱基质量低、比对模糊的那些；AD 只算能明确支持某个等位的读长，REF 排在前。9 个位点里 4 个存在差值，最大的一处差 4 条读长，其余 5 个相等。另外杂合位点的等位基因平衡 AB = alt 深度 / (ref 深度 + alt 深度) 得自己从 AD 推，真杂合应在 0.5 附近，这次有两个位点偏离到 0.83，落在需要回头看比对的区间。

第三，header 没声明的字段，bcftools 会直接报错而不是静默返回空。查 %INFO/AF 和 %GQ 都遇到 `Error: no such tag defined in the VCF header`，退出码 255。AF 用 `bcftools +fill-tags -- -t AF` 补齐（单样本杂合位点补出来是 0.5，等于 AC/AN），GQ 则从 PL 推导：GQ 就是两个最小 PL 值之差。还有一处细节：PL 里最小值 0 所在的下标就是 caller 判定的基因型，多等位拆分后要按 `k*(k+1)/2+j` 重新取值，不能按下标位置切。

读 VCF 字段前先确认两件事——这个字段是位点级还是每样本级，它在 header 里的 Number 是 A / R / G / . 中的哪一个。层级和 Number 定下来，字段的含义才定下来。
#生信 #生物信息学 #VCF #基因型 #bioSkills
