# 009｜bioSkills alignment-io：比对格式转换的注释掉失对比

<!--
META
用途: 009 alignment-io 小红书帖子「正文文本框」文案，与出图源稿配套，不进站点、不参与 md2card。
标题建议: 比对格式转换的注释处理
/ META
-->

做 MSA（多序列比对）下游工具要吃特定格式：FASTA、PHYLIP、Stockholm、NEXUS 各有不同。bioSkills 的 alignment-io 就是处理这些格式读写的，但它的实测结果表明，注释在不同格式之间不是天然保留。

拿示例数据试了下，做了两件事：把同一个 Stockholm 源文件写成多种格式再读回，看 GS/GR/GC 注释还在不在；再用长序列名测试 PHYLIP 严格模式，并检查了 MAF 负链坐标转换。

实测结果对比如下：

- Stockholm 自己：三类注释全部保留。FASTA、Clustal、PHYLIP-relaxed：三类全部丢失。NEXUS：只保留 GS，GR/GC 丢失。
- PHYLIP strict 会把 ID 硬截成 10 字符；两条序列名若在前 10 字符相同，BioPython 1.88 直接抛 ValueError，不是静默合并；phylip-relaxed 才能保留长名。
- MAF 负链坐标转换要小心 strand 字段：SKILL.md 示例按字符串 '-' 判断，可 BioPython 1.88 实际解析成整数 -1，原函数对负链返回的 plus-strand 起始坐标错误（10 bp vs 正确值 33 bp）。

所以做下游分析时，带注释的 MSA 务必留一份 Stockholm 主文件，转格式前先确认注释要不要；写 PHYLIP 前检查 ID 长度；MAF 坐标转换时把 strand 同时按字符串和整数处理。

#生信 #生物信息学 #多序列比对 #MSA #bioSkills
