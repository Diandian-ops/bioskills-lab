# 011｜bioSkills bowtie2-alignment：短读长比对模式怎么选

<!--
META
用途: 011 bowtie2-alignment 小红书帖子「正文文本框」文案，与出图源稿配套，不进站点、不参与 md2card。
标题建议: bowtie2 比对模式与适用场景
/META
-->

做生信短读长比对，bowtie2 是 ChIP / ATAC 这些峰实验的默认比对器。它最关键的不是「跑通」，而是选对模式：整条读长必须对齐（end-to-end 默认），还是读长末端允许软截（--local）。

实测三个发现：

① 接头污染的 reads 用默认 end-to-end，比对率从 99% 掉到 49%；切到 --local 软截后回血到 99.8%。差出来的 50 个百分点，就是被误罚的好核心。带接头读穿的样本优先 --local，或先剪切再比对。

② Bowtie2 的 MAPQ 上限是 42（end-to-end）/ 44（local），永远到不了 BWA 的 60。直接照搬变异流程里的「MAPQ ≥ 60 唯一映射」过滤，会把 BAM 清空。ChIP / ATAC 按惯例用 -q 30 丢多映射即可。

③ 峰实验要打开碎片几何旗标：--no-mixed --no-discordant 只留一致配对，ATAC 再加 --dovetail -X 2000 接纳短碎片跨核小体配对。这些旗标决定下游峰调用器看到的坐标，比核心比对更影响结果。

结论：bowtie2 比对 = 为峰实验选对模式与几何旗标。不确定时 ChIP 用 --very-sensitive 配 --no-mixed --no-discordant 再 -q 30；ATAC 切 --local --dovetail -X 2000。MAPQ 阈值跟着比对器走，别套 BWA 的 60。
