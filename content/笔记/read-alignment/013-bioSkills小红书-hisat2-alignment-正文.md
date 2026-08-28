# 013｜bioSkills hisat2-alignment：低内存 RNA 剪接比对首选

<!--
META
用途: 013 hisat2-alignment 小红书帖子「正文文本框」文案，与出图源稿配套，不进站点、不参与 md2card。
标题建议: RNA 比对选 HISAT2 省内存
/META
-->

做 RNA-seq 比对，STAR 功能最全但吃内存（人类参考约 30 GB）。机器内存紧张、或只需要剪接比对，HISAT2 是更省的选择——它的层级图索引只用约 7 GB 就能做近 STAR 的剪接比对。

实测三个发现：

① 剪接能力是真的。手工拼一条跨 800 bp 内含子的 100 bp 读长（前 50 bp 在 exon1 末端、后 50 bp 在 exon2 起始），HISAT2 把它对齐成 `50M800N50M`——中间的 N 就是跳过的内含子。一句话：它能把跨外显子 junction 的读长正确拼回去。

② MAPQ 给唯一定位读长打 60，对 GATK 友好。这一点和 bwa 一致，但和 STAR 不同：STAR 用 255 表示唯一映射，进 GATK 前还得重新赋值。HISAT2 直接给 60，下游变异调用少一步麻烦。

③ `--dta` 只给转录本组装用，普通计数别开。这个模式刻意抑制短锚定的 junction 读长，是给 StringTie / Cufflinks 组装用的；直接拿去做 featureCounts / htseq 计数会丢掉可用 junction 读长。链特异性也要按建库化学设（dUTP/TruSeq 用 RF），设错计数约减半。

结论：内存受限的 RNA 比对优先 HISAT2（~7 GB 图索引）；需要原生基因计数、融合检测或最高新 junction 灵敏度才上 STAR。allele-robust 映射用 SNP 图索引；`--dta` 仅用于组装流程；链特异性务必按建库设对。
