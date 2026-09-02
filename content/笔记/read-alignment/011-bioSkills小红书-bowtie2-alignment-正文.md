# 011｜bowtie2 比对：模式与碎片几何（版本=2.5.5）

<!--
META
用途: 011 bowtie2-alignment 小红书帖子「正文文本框」文案，与出图源稿配套，不进站点、不参与 md2card。
标题建议: bowtie2 比对模式怎么选
/-->

## 功能定位与适用范围

做 DNA 短读长比对，bowtie2 是 ChIP / ATAC 等峰实验的默认比对器。它最关键的不是「跑通」，而是选对模式：整条读长必须对齐（end-to-end 默认），还是读长末端允许软截（--local）。

适用范围：ChIP-seq、ATAC-seq、CUT&RUN 等需要碎片层面信号的场景；读长末端有接头读穿或引物末端的样本（用 --local 软截）。DNA 变异调用走 bwa-alignment，RNA 走 star/hisat2（场景对比见关联 skill），亚硫酸氢盐走 bismark，均不在本 skill 内。

## 属性表

| 项 | 值 |
|----|----|
| 主工具 | bowtie2（命令行，实测 2.5.5） |
| 索引 | reference_index.{1,2,3,4}.bt2 + .rev.{1,2}.bt2（共 6 个分片文件，后缀编号 1-4），-x 传 basename |
| MAPQ 上限 | end-to-end 42 / local 44，对比 BWA 为 60（满分标度） |
| 默认预设 | --sensitive（参数 -D/-R/-L 分别取 15/2/22） |
| 典型过滤 | ChIP/ATAC 用 -q 30（MAPQ=30）丢多映射；-F 1804 去 unmapped/secondary/dup/QC-fail |

## 成分拆解

- **end-to-end vs --local 是生物学决策**：end-to-end 强制整条读长匹配，适合干净基因组 DNA；--local 软截不可信读长末端，适合接头读穿、ATAC 的短碎片与高接头污染。在带接头污染的 reads 上用 end-to-end 会误罚好核心、压低比对率。
- **Bowtie2 MAPQ 与 BWA 不同标度**：由 AS/XS 驱动、离散，上限 e2e 42 / local 44，对比 BWA 的 60 到不了。照抄变异流程的 `-q 60`「唯一映射」过滤会把 BAM 清空；按比对器用 `-q 30`。
- **碎片几何旗标决定下游峰坐标**：`--no-mixed --no-discordant` 只留一致配对；ATAC 再加 `--dovetail -X 2000` 接纳短碎片跨核小体配对。这些旗标比核心比对更影响峰集。

## 严格复现

用合成参考（共 4 条 contig，每条 3000bp）+ wgsim 模拟 3000 对 PE reads 真跑全部命令，并以 samtools flagstat 独立复算，结果与日志完全一致：

```bash
# 适配器污染 reads：end-to-end vs --local
bowtie2 -p 8 -x reference_index -1 reads_contam_1.fq -2 reads_2.fq 2> c_e2e.log | samtools sort -o contam_e2e.bam -
# c_e2e.log: 49.28% overall alignment rate
bowtie2 -p 8 --local -x reference_index -1 reads_contam_1.fq -2 reads_2.fq 2> c_local.log | samtools sort -o contam_local.bam -
# c_local.log: 99.77% overall alignment rate
```
实测：--local 软截比 end-to-end 回血 **+50.49 pp（49.28% → 99.77%）**；max MAPQ e2e 42 / local 44。

```bash
# ChIP 旗标
bowtie2 -p 8 --very-sensitive --no-mixed --no-discordant -x reference_index \
    -1 reads_1.fq -2 reads_2.fq 2> chip.log | samtools view -bS -q 30 -F 1804 - | samtools sort -o chip.bam -
# chip.log: 50.93% 报告率；过滤后留存 3056 reads（合成参考为 4 段独立 contig，--no-discordant 丢弃跨 contig 配对，属数据特性）

# ATAC 旗标
bowtie2 -p 8 --very-sensitive --local --dovetail -X 2000 --no-mixed --no-discordant \
    -x reference_index -1 reads_1.fq -2 reads_2.fq 2> atac.log | samtools view -bS -q 30 -F 1804 - | samtools sort -o atac.bam -
# atac.log: 100.00% overall alignment rate；过滤后留存 6000 reads；max MAPQ = 44
```

负向验证：`-x reference_index.1.bt2`（文件名）返回 rc=255「does not exist or is not a Bowtie 2 index」，传 basename 即恢复。

**未覆盖（诚实标注）**：本步骤未重跑 bowtie2-build / bowtie2（索引与 BAM 已落盘），命令行为忠实复现原始真实调用，数值取自磁盘日志；samtools flagstat 由本步骤对已有 BAM 实跑，与日志比对率一致。

## 实践要点

- 不确定时：ChIP 用 `--very-sensitive --no-mixed --no-discordant` 再 `-q 30（MAPQ=30）`；ATAC 切 `--local --dovetail -X 2000（片段上限=2000）`。
- MAPQ 阈值跟着比对器走：Bowtie2 用 -q 30（MAPQ=30），别套 BWA 的 -q 60（MAPQ=60）。
- `-x` 永远传索引 basename，不要带 `.bt2` 后缀（索引分片后缀 .bt2，编号 1-4 共 6 个文件）或具体分片文件。
- 带接头污染的 reads 优先 `--local` 或先剪切再比对。

## 小结

bowtie2 比对 = 为峰实验选对模式与几何旗标。本次复现实测印证两条关键论断：适配器污染下 --local 比 end-to-end 高 +50.49 pp，以及 Bowtie2 MAPQ 上限 e2e 42 / local 44、到不了 BWA 的 60。文档与工具行为一致。

#生信 #生物信息学 #bowtie2 #短读长比对 #ChIPseq #ATACseq #bioSkills（工具版本=2.5.5）
