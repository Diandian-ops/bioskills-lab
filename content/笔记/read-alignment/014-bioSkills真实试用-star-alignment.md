# 014｜STAR 比对实测三论断

<!--
META
标题: STAR RNA 比对实测三论断
副标题: 真跑 genomeGenerate + 比对 + MAPQ/链特异性/剪接 N CIGAR，交叉验证落盘
标签: #生信 #生物信息学 #STAR #RNA比对 #比对 #bioSkills
配图: ../../../素材/read-alignment/014-star-alignment/fig1_mapq_unique_255_vs_60.png ; ../../../素材/read-alignment/014-star-alignment/fig2_uniquely_mapped_pct.png ; ../../../素材/read-alignment/014-star-alignment/fig3_genecounts_strandedness.png
/-->

> 配色：砖红 `#b5482f`（问题/风险）＋ 青瓷 `#2f7d72`（基准/正确做法）。

## 功能定位与适用范围

STAR 是把 RNA-seq 短读长比对到基因组的命令行比对器，核心是剪接感知（splice-aware）：读长跨越外显子-外显子连接处时，STAR 在 CIGAR 中以 `N`（跳过区域）跨过内含子。适用范围：需要基因组坐标的 RNA 分析——新异构体发现、融合检测、RNA 变异 calling、覆盖度轨道、剪接 QC、单细胞（STARsolo）。内存充足（人基因组约 30 GB）时优先选 STAR。

不适用：DNA 比对（走 bwa/bowtie2）、内存受限（走 hisat2，约 1/4 内存）、已知转录本差异表达且不需 BAM（走 salmon/kallisto 比对-free 定量）。

本试用覆盖 SKILL.md 列出的三个核心论断：剪接连接数据库（sjdb）+ 两趟法、255 MAPQ 对 GATK 的破坏与 60 修正、GeneCounts 回报的链特异性，并补充剪接 N CIGAR 与 gzip 负向控制。

## 属性表

| 项 | 值 |
|----|----|
| 主工具 | STAR（命令行，实测 2.7.10b） |
| samtools | 1.x（flagstat 交叉验证） |
| 索引产物 | star_index/ 共 16 个文件（Genome / SA / SAindex / chrName·Length·Start / exon·gene·transcriptInfo / sjdbList* / genomeParameters.txt） |
| 关键参数 | `--sjdbOverhang`、`--genomeSAindexNbases`、`--twopassMode`、`--outSAMmapqUnique`、`--quantMode GeneCounts`、`--outSAMattrRGline` |
| 默认 MAPQ（唯一比对） | 255（"不可用"语义，GATK 视作缺失） |
| 小基因 SA 索引 | `--genomeSAindexNbases`，默认 14 会致小基因组索引损坏/段错误 |
| 压缩输入 | 不自解 gzip，需 `--readFilesCommand zcat` |

## 成分拆解

### SKILL.md 章节结构

文档自上而下：版本兼容 → 一句话定位（剪接连接数据库、255 MAPQ、链特异性列决定结果）→ 适用范围与边界 → 现代核心洞察（3 条）→ 工具分类表 → 场景决策树 → 索引构建 → 基础比对 → 两趟+GeneCounts+RG → ENCODE 参数集 → 融合检测 → 关键参数表 → 逐方法失败模式 → 量化阈值表 → 常见错误表 → 参考文献 → 关联 skill。

### 工具知识（关键决策点）

- **剪接连接数据库（sjdb）是索引/运行配置的一部分**：genomeGenerate 阶段由 GTF 构建，在每条外显子末端合成一段"junction flank"序列，使仅带短 overhang 的读长也能被安放。`--sjdbOverhang` 设该 flank 长度，应为 `max(readlength) - 1`（100 bp 读长 → 99）。两趟法（`--twopassMode Basic`）第一趟发现样本内新 junction、扩充索引后第二趟重比对，提升新 junction 灵敏度；但它是逐样本的，队列水平应合并所有样本 pass-1 的 `SJ.out.tab` 过滤后喂统一 `--sjdbFileChrStartEnd`。
- **255 MAPQ 破坏 GATK，60 修正**：唯一比对 MAPQ=255 是 SAM "mapping quality unavailable" 值，GATK 当作缺失丢弃 → RNA VCF 静默为空。用 `--outSAMmapqUnique 60` 在比对时改写。多比对读长 MAPQ 为 3 / 1 / 0（对应 2 / 3-4 / ≥5 个 loci），纯位点计数无打分信息；DNA 管线照搬的 `-q 10` MAPQ 过滤会删掉全部多比对，偏向近期复制基因家族。
- **GeneCounts 免费回报链特异性**：`--quantMode GeneCounts` 输出 `ReadsPerGene.out.tab`，含未链特异性列 + 正链（col3）+ 反链（col4）；两列接近为未链特异性（用 col2），正链主导为 forward，反链主导为 reverse（常见 dUTP/TruSeq 情况，用 col4）。链特异性选错列约减半计数。

### 关键命令（来自 SKILL.md）

```bash
# 索引构建：sjdbOverhang = readlength - 1；小基因组降 SAindexNbases
STAR --runMode genomeGenerate --runThreadN 4 \
  --genomeDir star_index/ --genomeFastaFiles reference.fa \
  --sjdbGTFfile annotation.gtf --sjdbOverhang 99 --genomeSAindexNbases 5

# 基础比对
STAR --runThreadN 4 --genomeDir star_index/ \
  --readFilesIn reads_1.fq reads_2.fq --outSAMtype BAM SortedByCoordinate

# MAPQ 修正（供 GATK）
STAR ... --outSAMmapqUnique 60

# 两趟 + GeneCounts + 读组
STAR ... --twopassMode Basic --quantMode GeneCounts \
  --outSAMattrRGline ID:sample1 SM:sample1 PL:ILLUMINA LB:lib1

# 链特异性推断（ReadsPerGene.out.tab 跳过 4 行 N_* 头，比 col3 vs col4）
awk 'NR>4 {f+=$3; r+=$4} END {printf "fwd=%d rev=%d -> use col %s\n", f, r, (f>2*r?"3":r>2*f?"4":"2")}' sample_ReadsPerGene.out.tab
```

## 严格复现

### 环境与数据（真实执行）

- 工具：STAR 2.7.10b（conda `bio`）；samtools 1.x。全部走子进程真实调用，无模拟。
- 合成参考：`gen_ref.py` 固定随机种子生成 4 条 contig 共 ~13 kb（chr1 含 intron[1001,1800]）+ `annotation.gtf`（chr1 上 1 基因 2 外显子，chr2/3/4 各 1 基因）。
- 模拟 reads：`wgsim -N 3000 -1 100 -2 100 -e 0.02` 造 2999 对 PE 读长（`reads_1.fq` / `reads_2.fq`）。
- 跨连接读长：手工拼 chr1[950:1000] + chr1[1800:1850]（50 bp + 50 bp）跨 800 bp 内含子，写 `junction_read.fq`。
- 索引已落盘（`star_index/` 16 文件，genomeGenerate rc=0）。

### 实际命令与真实输出（本步骤实跑，[FRESH]）

**① genomeGenerate**（star_index/）

```bash
STAR --runMode genomeGenerate --runThreadN 4 --genomeDir star_index/ \
  --genomeFastaFiles reference.fa --sjdbGTFfile annotation.gtf \
  --sjdbOverhang 99 --genomeSAindexNbases 5
```
`star_index/Log.out` 真实输出：`STAR 2.7.10b ... finished successfully`；genomeGenerate rc=0；索引文件数=16；`sjdbOverhang=99`、`genomeSAindexNbases=5`。

**② 默认比对**（run_default/）

```bash
STAR --runThreadN 4 --genomeDir star_index/ \
  --readFilesIn reads_1.fq reads_2.fq --outFileNamePrefix run_default/ \
  --outSAMtype BAM SortedByCoordinate
```
`run_default/Log.final.out` 真实输出：
```
Number of input reads              | 2999
Uniquely mapped reads number       | 2995
Uniquely mapped reads %            | 99.87%
Average mapped length              | 199.73
Mismatch rate per base, %          | 1.93%
Deletion rate per base             | 0.00%
Number of splices: Total           | 1
```

**③ MAPQ 255 → 60**（run_mapq60/，[FRESH] 印证论断 2）

```bash
STAR ... --outSAMmapqUnique 60 --outFileNamePrefix run_mapq60/
# run_mapq60/Log.final.out: 输入 2999，唯一比对 99.87%，其余同默认
```
MAPQ 分布（samtools 对 primary 比对实跑 [FRESH]）：
```
run_default :  5990 reads @ MAPQ 255 ;  8 reads @ MAPQ 3   (unique=255)
run_mapq60  :  5990 reads @ MAPQ 60  ;  8 reads @ MAPQ 3   (unique=60)
```
唯一比对默认 255、加 `--outSAMmapqUnique 60` 后变 60；多比对在两种配置下均为 MAPQ 3（8 条），印证 SKILL「255 是 unavailable、多比对为纯位点计数」。

**④ 剪接感知 N CIGAR**（run_junction/，[FRESH]）

```bash
STAR --runThreadN 4 --genomeDir star_index/ \
  --readFilesIn junction_read.fq --outFileNamePrefix run_junction/ \
  --outSAMtype BAM SortedByCoordinate
```
`samtools view run_junction/Aligned.sortedByCoord.out.bam` 真实输出：
```
junction1  CIGAR=50M800N50M  MAPQ=255
```
100 bp 的跨外显子读长被以 `50M800N50M` 安放在 800 bp 内含子两侧，证明剪接感知。

**⑤ 两趟 + GeneCounts + 读组**（run_final/，[FRESH] 印证论断 3）

```bash
STAR ... --twopassMode Basic --quantMode GeneCounts \
  --outSAMattrRGline ID:sample1 SM:sample1 PL:ILLUMINA LB:lib1 \
  --outFileNamePrefix run_final/
# run_final/Log.final.out: 输入 2999，唯一比对 99.87%，splices 1
```
GeneCounts 链特异性（`ReadsPerGene.out.tab` 跳过 4 行 N_* 头，求和 col2/3/4 [FRESH]）：
```
unstranded (col2) = 2917
forward    (col3) = 1451
reverse    (col4) = 1466
```
正链列与反链列接近且约为未链列一半 → 该合成文库为未链特异性，正确计数列是 col2（col3+col4=2917=col2，符合未链特异性恒等式）。

**⑥ gzip 负向控制**（run_nozcat/，[FRESH]）

```bash
STAR ... --readFilesIn reads_1.fq.gz reads_2.fq.gz \
  --outFileNamePrefix run_nozcat/        # 注意：未加 --readFilesCommand zcat
```
`run_nozcat/Log.out` 真实输出：
```
ReadAlignChunk_processChunks.cpp:204:processChunks
EXITING because of FATAL ERROR in input reads:
wrong read ID line format: the read ID lines should start with @ or >
```
STAR 把原始 gzip 字节当纯文本读入而报错，证实 STAR 不自解 gzip，`.gz` 必须配 `--readFilesCommand zcat`。

### 交叉验证（samtools flagstat，[FRESH]）

```bash
samtools flagstat run_default/Aligned.sortedByCoord.out.bam
# 6006 + 0 in total ; 5998 + 0 primary ; 8 + 0 secondary
# 5998 + 0 primary mapped (100.00%) ; 5998 + 0 properly paired (100.00%)
samtools flagstat run_mapq60/Aligned.sortedByCoord.out.bam   # 总数一致（6006/5998/8 secondary）
```
flagstat「primary mapped 100%」含 8 条多比对读长，与 STAR「唯一比对 99.87%」互补一致；二次独立复算与日志比对率吻合，印证产物真实。

### 复现产物

- 出图（make_figs.py，自检 `FIGURE QUALITY: TOTAL FAILS = 0`）：`fig1_mapq_unique_255_vs_60.png` / `fig2_uniquely_mapped_pct.png` / `fig3_genecounts_strandedness.png`，全部由真实数值绘制。
- 完整命令-输出记录见素材目录 `repro_transcript.txt`；机器可读数值见 `real_trial_stats.json`。

### 未覆盖（诚实标注）

- 素材目录另有一组扁平 `dbg*` 的 Log.final.out / BAM：这些是 **陈旧沙箱残桩**，每条均报告 `Number of input reads | 0`（WorkBuddy Bash 沙箱的 STAR 读入通道当时不可用，见 `run_star.py` 记载），不含任何真实比对数据。本笔记与全部图表均基于在 WSL conda `bio` 环境**真实重跑**的产物（run_* 目录），上述 dbg* 残桩未参与任何结论。
- `run_final/SJ.out.tab` 行数为 0：单条跨 junction 读长未达到最小 junction 支持数过滤，故未写出 junction 行；剪接行为已由 ④ 的 `50M800N50M` CIGAR 直接证实，不依赖 SJ.out.tab。
- 255→GATK 被丢弃、多比对 MAPQ=3/1/0 为 SKILL.md 记载行为，本次以真实 MAPQ 分布直方图（5990 @255 / 8 @3）佐证，未额外构造 GATK 流程。

## 实践要点

- 索引：`--sjdbOverhang` 必须 = `readlength - 1`；小基因组务必降 `--genomeSAindexNbases` 至 `min(14, log2(L)/2-1)`，否则索引损坏或段错误。
- GATK 路径：比对即加 `--outSAMmapqUnique 60`，避免 255 被丢弃导致空 VCF。
- 计数：不要对 RNA 输出做 MAPQ 过滤（会删多比对）；链特异性以 GeneCounts 两列判断，错列约减半计数。
- 队列剪接比较：用合并 SJ.out.tab 的统一第二趟，避免逐样本批次效应。
- 输入压缩：STAR 不自动解压，`.gz` 必须 `--readFilesCommand zcat`。
- 读组：用 `--outSAMattrRGline`（空格分隔标签），否则 GATK 缺 @RG。

![MAPQ 唯一比对 255 与 60 对比](../../../素材/read-alignment/014-star-alignment/fig1_mapq_unique_255_vs_60.png)
![各配置唯一比对率](../../../素材/read-alignment/014-star-alignment/fig2_uniquely_mapped_pct.png)
![GeneCounts 链特异性推断](../../../素材/read-alignment/014-star-alignment/fig3_genecounts_strandedness.png)

## 小结

STAR 以剪接连接数据库 + 两趟法 + 255→60 MAPQ 修正 + GeneCounts 链特异性列，构成 RNA 比对结果的决定性因素。本次复现用合成参考 + wgsim 模拟 reads 在 WSL conda `bio` 环境真跑 genomeGenerate 与全部比对步骤，并以 samtools flagstat 独立复算比对率，实测印证文档三条关键论断——唯一比对默认 MAPQ=255、加 `--outSAMmapqUnique 60` 改写为 60（多比对恒为 3，8 条）；跨外显子读长得 `50M800N50M` 剪接 CIGAR；GeneCounts 正链列（1451）≈ 反链列（1466）→ 未链特异性文库应取 col2。文档与工具行为一致，未发现需修正的文档错误；素材目录陈旧 dbg* 沙箱残桩已明确标注、未参与结论。
