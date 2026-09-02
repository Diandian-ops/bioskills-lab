# 012｜DNA 短读长比对（含 3 条核心契约）

<!--
META
用途: 012 bwa-alignment 小红书帖子「正文文本框」文案，与出图源稿配套，不进站点、不参与 md2card。
标题建议: bwa 比对必写 read group
/META
-->

做 DNA 短读长比对，bwa-mem2（常见工具）是变异调用流程的默认比对器。它最关键的不是「跑通」，而是比对前就把下游要用的契约写好——否则后续 GATK 直接报错或样本被静默混在一起。

本次复现用合成参考（24000 bp）+ wgsim 生成 3000 对 PE 读长（共 6000 条）执行了全部命令，结果如下。

## 功能定位与适用范围

bwa-mem2 把 DNA 短读长（单端 / 双端）比对到参考基因组，是 bwa-mem 的加速版，输出近一致；服务于变异调用、覆盖度、ChIP/ATAC、SV 检测。RNA 走 star/hisat2，长读长与亚硫酸氢盐不在此范围。

## 属性表

| 项 | 值 |
|----|----|
| 主工具 | bwa-mem2（常见比对器） |
| 索引产物 | .0123 / .amb / .ann / .bwt.2bit.64 / .pac（共 5 个文件，与 bwa index 的 .bwt/.sa 不互通） |
| MAPQ 标度 | 0–60，双峰集中在 0（多映射）与 60（无竞争位点） |
| 去重顺序 | collate → fixmate -m → sort → markdup |

## 成分拆解

三个会静默决定下游成败的决策：

① read group 是硬契约，不是装饰。比对时用 `-R` 注入 `@RG`（含样本名 SM、误差单元 ID、平台 PL、文库 LB），生成的 BAM 头里就有 @RG；不写，BAM 头里就完全没有 @RG。GATK 靠 SM 分组样本、靠 ID/LB 建模误差和去重，缺了要么报错要么行为错误，事后用 Picard 补是一次完整 BAM 重写。

② MAPQ 标度到 60，和 bowtie2 不一样。同一批数据实测 bwa-mem2 的 max MAPQ=60（bwa 序数标度到顶）；bowtie2 上限只有 42（end-to-end）/ 44（local）。所以「MAPQ 高值唯一映射」过滤在 bwa 上合理，照搬去压 bowtie2 的 BAM 会清空——阈值要跟着比对器走。

③ 去重顺序错了会静默标错。正确顺序是 collate（按名）→ fixmate -m（写 MC 标签）→ sort（坐标）→ markdup。fixmate 不带 `-m` 或拿坐标序 BAM 直接 markdup，重复标记会悄悄出错；SV 流程还要用 `-Y` 保留分裂读长全程序列，不靠 `-M`（会把 SV 证据降级藏起来）。

## 严格复现

真实命令（节选）：

```bash
# 注入 read group，流式排为坐标序 BAM
bwa-mem2 mem -t 8 -R '@RG\tID:s1\tSM:s1\tPL:ILLUMINA\tLB:lib1' \
    reference.fa r1.fq r2.fq | samtools sort -@4 -o aligned_rg.bam -

# 去重严格顺序
bwa-mem2 mem -t 8 -R '@RG...' r1 r2 | \
    samtools collate -@4 -O -u - | samtools fixmate -m -@4 -u - - | \
    samtools sort -@4 -u - | samtools markdup -@4 - aligned.markdup.bam
```

真实输出（节选）：

```text
# 带 -R 比对：read 6000 sequences (600000 bp)，Processed 6000 reads
# 不带 -R：BAM 头无 @RG（契约负向印证）
# 索引错配：ERROR! Unable to open the file: .../ref.bwt.2bit.64
```

复现得到的映射统计（重跑 samtools flagstat，五个配置一致）：6000/6000 映射率 100%、6000/6000 正确配对率 100%；仅 markdup 配置经严格顺序标出 2 个重复标记（占 6000 的 0.03%）。`-K 100000000` 两次跑对齐指纹 md5 完全一致。

## 实践要点

- 比对即注入 `@RG`（SM/ID/PL/LB），别等 GATK 报错再补。
- 人类数据用 GRCh38 + decoy（hs38d1，常见分析集），HLA 敏感场景加 ALT 与 `.alt` 做 ALT-aware 映射。
- 去重严格走 `collate → fixmate -m → sort → markdup`；amplicon/PCR 别做坐标去重（用 UMI）。
- SV 用 `-Y` 不用 `-M`；要可复现加 `-K 100000000`（常见做法）。
- bwa 与 bwa-mem2 索引工具配对，别混用（常见）。

## 小结

bwa 比对 = 为下游变异调用写对契约。比对即注入 @RG；人类数据用 GRCh38 + decoy；去重严格走 collate→fixmate -m→sort→markdup；SV 用 -Y 不用 -M；要可复现加 -K 100000000。本次复现实测印证：带 -R 的 BAM 含 @RG、缺 -R 则无；max MAPQ=60；-K 两次跑指纹一致；bwa 原版索引喂给 bwa-mem2 直接报错。

#生信 #生物信息学 #bwa #短读长比对 #bioSkills
