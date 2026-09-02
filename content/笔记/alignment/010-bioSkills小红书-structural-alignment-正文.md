# 结构比对工具链实测

<!--
META
用途: 010 structural-alignment 小红书帖子「正文文本框」文案，与出图源稿配套，不进站点、不参与 md2card。
标题建议: 结构比对 TM-align 实测
/META
-->

序列一致性太低（低于 25%）时，序列比对已经不可靠，这时候要看结构。bioSkills 的 structural-alignment 就是解决「结构怎么比、怎么判同源」的一套工具链，覆盖 TM-align、Foldseek、US-align、DALI、Foldmason 这些常用方案。

拿 RCSB 上真实的四个小球蛋白结构做了实测：1UBQ（泛素，76 个残基）、1CRN（46 个残基）、1ENH（54 个残基）、1R69（63 个残基），互为跨折叠对照。

实测结论如下：

- TM-align 两两打分（1UBQ 分别对比其余三个）：TM-score 依次为 0.361、0.389、0.337，全部低于 0.5 折叠阈值，判定为不同折叠；对应 RMSD 约 3.1–3.4 Å，比对长度仅 32–40 个残基。
- TM-score 比 RMSD 更能反映折叠相似度：本例 RMSD 都不算极端高，但比对长度短（远低于链长 46–76）， pulled down 了 TM-score，所以判折叠必须同时看 TM-score 和比对长度。
- Foldseek easy-search 对这四个本地结构做自比，默认 3Di+AA 检索只返回自身命中（相似度 1.000），没有任何跨结构命中——和 TM-align「都不同折叠」的结论一致。低相似度跨折叠探测要开 `--alignment-type 1`（TMalign 全局）或放宽阈值。

做结构比对时记住三点：TM-score 取两个长度归一化分数里的较大值（按较短链归一化）；TM>0.5 才算同折叠，别只看 RMSD；Foldseek 默认检索口径偏保守，跨折叠候选要显式全局对齐复核。

#生信 #生物信息学 #结构比对 #TMalign #Foldseek #bioSkills
