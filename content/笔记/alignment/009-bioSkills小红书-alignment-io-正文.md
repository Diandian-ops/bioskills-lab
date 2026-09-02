---
title: "比对格式：注释存活看格式"
skill: alignment-io
trial: "009"
type: "xhs-body"
category: "bioSkills 真实试用"
tags: ["Biopython", "AlignIO", "alignment-io", "生信", "bioSkills"]
date: "2026-09-02"
---

封面卡：比对格式不是随便转，带注解的 MSA 转错格式会静默丢信息

用 BioPython 的 AlignIO 真跑 alignment-io SKILL.md 的读写与转换函数，输入是多序列比对（MSA）示例：4 条 × 20 列的小 DNA 比对，把格式互操作要点记下来。

发现一：注释存活看格式。Stockholm 能保留 SS_cons 这类列注解（本例 20 列二级结构串 round-trip 后还在），FASTA 导出后 SS_cons 和 GS 序列注解一起静默丢失。带注解的 MSA 一定要留 Stockholm 主文件，别拿 FASTA 当权威源。

发现二：NEXUS 转换必须带字母表。把 Stockholm 转 NEXUS 时，不声明 molecule_type 会直接抛 ValueError；带上 molecule_type='DNA' 才成功（1 个比对写出）。DNA/RNA/蛋白序列写 NEXUS 前记得先设字母表。

发现三：格式支持很全但现代 API 取列数用法不同。7 种格式（fasta/clustal/phylip/phylip-relaxed/phylip-sequential/stockholm/nexus）实测 read=write=True。Bio.Align 现代 API 的 Align.read 返回 Alignment 对象，取列数要用 .length，老 SKILL.md 的 get_alignment_length() 在新对象上不存在；程序化用 MultipleSeqAlignment 也能拼出 3 条 × 12 列 MSA。

结论：格式互操作先想清楚下游要什么——RAxML/IQ-TREE 用 phylip-relaxed、MrBayes 用 nexus、HMMER/Pfam 用 stockholm。带注解留 Stockholm，NEXUS 设 molecule_type，现代 API 用 .length。MAF 负链坐标、A2M/A3M、pyhmmer 流式在 examples 里没逐行复现。

#生信 #生物信息学 #Biopython #MSA #bioSkills
