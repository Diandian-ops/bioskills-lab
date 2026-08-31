<!--
META
用途: 004 pairwise-alignment 小红书帖子「正文文本框」文案，与出图源稿配套，不进站点、不参与 md2card。
标题建议: pairwise 比对缺口罚分的默认行为
/ META
-->

拿到两条同源序列，想知道差多少、哪些位置保守，第一步是做两两比对（pairwise alignment）。bioSkills 的 pairwise-alignment 就是教怎么把这一步做对的。

拿示例数据（人血红蛋白 α 链 vs β 链，NCBI P69905 / P68871）试了下，用 Biopython 的 PairwiseAligner，踩到两个点：

一是默认不设罚分，gap 全 0，比对会疯狂插无意义短 gap 凑分。同样这两条序列，默认对齐出 172 列、塞了 55 个 gap；按 BLOSUM62 推荐 open=-11 / extend=-1，只有 149 列、9 个 gap——多出来 46 个假 gap（+15.4%）。

二是相似度（PID）有四种定义，同一份比对能差出 2.8 个百分点，极端情况达 11.5%。报 PID 必须说清口径。

所以 pairwise 比对别用默认值，BLOSUM62 配 open=-11 / extend=-1，并明确 PID 口径。

#生信 #生物信息学 #序列比对 #PairwiseAlignment #bioSkills
