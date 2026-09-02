# 014｜STAR RNA 比对三处易错（版本=2.7.10b）

<!--
META
用途: 014 star-alignment 小红书帖子「正文文本框」文案，与出图源稿配套，不进站点、不参与 md2card。
标题建议: STAR RNA 比对别漏三点
/-->

## 功能定位与适用范围

做 RNA-seq 比对，STAR 是剪接感知比对器的代表：读长跨外显子连接处时，它用 CIGAR 里的 `N` 跳过内含子。它最关键的不是「跑通」，而是三个直接影响下游计数与变异 call 的点。

适用范围：需要基因组坐标的 RNA 分析（异构体、融合、RNA 变异、剪接 QC、单细胞）。DNA 走 bwa/bowtie2；内存受限走 hisat2；只做已知转录本差异表达可走 salmon/kallisto。

## 属性表

| 项 | 值 |
|----|----|
| 主工具 | STAR（命令行，实测 2.7.10b） |
| 索引 | star_index/ 共 16 个文件；`--sjdbOverhang` = 读长-1 |
| 默认 MAPQ | 唯一比对 = 255（语义为"质量未知"），多比对 = 3 / 1 / 0 |
| 链特异性 | `--quantMode GeneCounts` 免费回报未链/正链/反链三列 |
| 压缩输入 | 不自动解压，`.gz` 需 `--readFilesCommand zcat` |

## 成分拆解

- **连接数据库 overhang 要对**：索引构建要带 GTF，且 `--sjdbOverhang` 必须等于读长减 1（100 bp 读长就填 99）；小基因组还要把 `--genomeSAindexNbases` 调小，用默认 14 索引会坏。
- **MAPQ 255 坑 GATK**：STAR 给唯一比对打 255，GATK 当作"质量未知"直接丢弃，RNA VCF 会变成空文件；比对时加 `--outSAMmapqUnique 60` 就能避免。多比对读长 MAPQ 只有 3/1/0，是纯位点计数——照抄 DNA 的 `-q 10` 过滤会删掉全部多比对。
- **链特异性别猜**：开 `--quantMode GeneCounts`，看 ReadsPerGene 的正链列/反链列：两列接近是未链特异性，反链主导就是常见 dUTP/TruSeq 方案，选错列计数约减半。

## 严格复现

用合成参考（约 13 kb，chr1 含 800 bp 内含子）+ wgsim 模拟 2999 对（reads 数=2999）PE reads，在 WSL conda `bio` 环境真跑全部命令，并以 samtools flagstat 独立复算，结果与日志一致：

```bash
# 唯一比对 MAPQ：默认 255 vs 加 --outSAMmapqUnique 60
STAR ... --readFilesIn reads_1.fq reads_2.fq --outSAMtype BAM SortedByCoordinate
# Log.final.out: 输入 2999 对，唯一比对率 = 99.87%
STAR ... --outSAMmapqUnique 60
# samtools 看 MAPQ：默认 5990 条 @255，改后 5990 条 @60；多比对恒为 8 条 @3
```

实测印证两条关键论断：唯一比对默认 MAPQ=255、加参数改写为 60（多比对 8 条恒为 MAPQ=3）；跨外显子读长得到 `50M800N50M` 的剪接 CIGAR（N=跳过区域，跨度=800 bp 内含子）。

```bash
# 链特异性：ReadsPerGene.out.tab 跳过 4 行头，比正链列(col3) vs 反链列(col4)
# 实测：未链列=2917，正链列=1451，反链列=1466 -> 正≈反，未链特异性，取 col2
```

负向验证：`.gz` 输入不带 `--readFilesCommand zcat`，STAR 报 `wrong read ID line format` 直接退出，证明它不自动解压。

**未覆盖（诚实标注）**：素材目录一组陈旧 `dbg*` 残桩曾是沙箱空读长（reads 数=0）输出，本次未采用；全部结论来自 WSL 环境真实重跑的 run_* 产物，flagstat 二次复算与日志吻合。

## 实践要点

- 索引：`--sjdbOverhang` 设为读长减 1；小基因组降 `--genomeSAindexNbases`，否则索引损坏。
- GATK 路径：比对即加 `--outSAMmapqUnique 60`（MAPQ=60），避免 255（MAPQ=255）被丢弃导致空 VCF。
- 计数：别对 RNA 做 MAPQ 过滤（会删多比对）；链特异性以 GeneCounts 两列判断，错列约减半。
- 输入压缩：STAR 不自动解压，`.gz` 必须 `--readFilesCommand zcat`。
- 读组：用 `--outSAMattrRGline`（空格分隔）填 ID/SM/PL/LB，否则 GATK 缺 @RG。

## 小结

STAR 把剪接连接、MAPQ、链特异性都打包了，但参数不是默认就安全。本次复现实测印证三点：唯一比对 MAPQ 默认 255、加 `--outSAMmapqUnique 60` 改 60；跨外显子读长得到 `50M800N50M` 剪接 CIGAR；GeneCounts 正链列（1451）≈ 反链列（1466）→ 未链特异性文库取 col2。文档与工具行为一致。

#生信 #生物信息学 #STAR #RNA比对 #剪接比对 #bioSkills（工具版本=2.7.10b）
