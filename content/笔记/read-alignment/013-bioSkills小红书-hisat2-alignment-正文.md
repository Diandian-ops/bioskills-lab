# HISAT2 对比 STAR：低内存 RNA 剪接比对首选

<!--
META
用途: 013 hisat2-alignment 小红书帖子「正文文本框」文案，与出图源稿配套，不进站点、不参与 md2card。
标题建议: HISAT2：低内存 RNA 比对
/META
-->

做 RNA-seq 比对，STAR 功能最全但吃内存（内存对比：STAR 约 30 GB、HISAT2 约 7 GB，差约 4 倍）；机器内存紧张、或只需要剪接比对，HISAT2 是更省的选择——它的层级图索引只用约 7 GB（内存对比差约 4 倍）就能做近 STAR 的剪接比对。

用合成参考 + 手工跨 junction 读长复现了全部命令，三个实测发现：

① 剪接能力是真的。手工拼一条跨 800 bp 内含子的 100 bp 读长（前 50 bp 在 exon1 末端、后 50 bp 在 exon2 起始，共一条），HISAT2 把它对齐成 `50M800N50M`——中间的 N 就是跳过的 800 bp 内含子（一条）。一句话：它能把跨外显子 junction 的读长正确拼回去。

② MAPQ 给唯一定位读长打 60（=MAPQ 上限），对 GATK 友好。这一点和 bwa 一致，但和 STAR 不同：STAR 用 255 表示唯一映射（=STAR 标度），进 GATK 前还得重新赋值。HISAT2 直接给 60（=MAPQ 上限），下游变异调用少一步麻烦。

③ `--dta` 只给转录本组装用，普通计数别开。这个模式刻意抑制短锚定的 junction 读长，是给 StringTie / Cufflinks 组装用的；直接拿去做 featureCounts / htseq 计数会丢掉可用 junction 读长（实测比对率从 99.17% 掉到 98.05%，未比对读长翻倍到 109 对）。链特异性也要按建库化学设（dUTP/TruSeq 用 RF），设错计数约减半。

结论：内存受限的 RNA 比对优先 HISAT2（~7 GB 图索引；内存对比 STAR 约 4 倍）；需要原生基因计数、融合检测或最高新 junction 灵敏度才上 STAR。allele-robust 映射用 SNP 图索引；`--dta` 仅用于组装流程；链特异性务必按建库设对。
#生信 #生物信息学 #HISAT2对比STAR #RNA比对 #bioSkills
