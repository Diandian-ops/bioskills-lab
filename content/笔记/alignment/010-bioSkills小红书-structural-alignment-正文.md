# 010｜bioSkills structural-alignment：结构比对工具实测

<!--
META
用途: 010 structural-alignment 小红书帖子「正文文本框」文案，与出图源稿配套，不进站点、不参与 md2card。
标题建议: 结构比对 TM-align 用法概览
/META
-->

序列一致性太低（低于 25%）时，序列比对已经不可靠，这时候要看结构。bioSkills 的 structural-alignment 就是解决「结构怎么比、怎么判同源」的一套工具链，覆盖了 Superimposer、TM-align、Foldseek 这些常用方案。

拿 RCSB 上真实的结构做了实测：1ubq 和 1ubi 都是泛素（ubiquitin），1fmb 是另一个折叠，作为对照。

实测结论如下：

- Bio.PDB.Superimposer 已知残基对应时直接算 RMSD：1ubq↔1ubi 仅 0.09 Å（76 个 CA 几乎完全重合），说明二者几乎同一个结构。
- TM-align 两两打分：1ubq↔1ubi 的 TM-score 0.999（>0.5，同一折叠）；1ubq↔1fmb 的 TM-score 0.408（<0.5，不同折叠）。TM-score 比 RMSD 更能反映折叠相似度。
- Foldseek easy-cluster 把 1ubq 和 1ubi 聚到一起、1fmb 单独成簇，和 TM-align 结论一致；但 Foldseek 默认检索会过滤掉 TM-score 低于 0.5 阈值 的跨折叠候选，低相似度探测要额外开全局对齐复核。

踩到一个文档坑：SKILL.md 里给的 `--format-output` 用了 `tmscore`，实际二进制会报「Format code tmscore does not exist」，正确应写 `alntmscore`/`qtmscore`/`ttmscore`。

做结构比对时记住：TM-score>0.5 才算同折叠，别只看 RMSD；Foldseek 表格列名别照抄文档里的 tmscore。

#生信 #生物信息学 #结构比对 #TMalign #Foldseek #bioSkills
