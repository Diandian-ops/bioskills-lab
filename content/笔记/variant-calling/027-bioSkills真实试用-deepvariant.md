# 027 · bioSkills 真实试用：deepvariant（深度学习变异检测）

## 功能定位与适用范围

`deepvariant` 讲解**用 Google DeepVariant 检测胚系 SNP 与 indel**——它把变异检测重新表述为对多通道 pileup 张量做的 CNN 图像分类。

- **适用**：选择平台特异的模型（WGS、WES、PACBIO、ONT_R104、HYBRID_PACBIO_ILLUMINA）；一次性 `run_deepvariant` 与三阶段（make_examples / call_variants / postprocess_variants）流程的选择；用 GPU 加速 call_variants；用 DeepTrio 做家系/三人组与新发变异检测；用 GLnexus（**不是** GenotypeGVCFs）联合基因分型 DeepVariant 的 gVCF。
- **不适用**：体细胞调用（用 DeepSomatic）。

## 属性表（本次环境）

| 项 | 值 |
|---|---|
| DeepVariant | **未安装** |
| Docker | 可执行文件存在（`/usr/local/bin/docker`），但**守护进程未运行**（`docker info` 失败） |
| 比对数据 | **无 BAM/CRAM** |
| 参考基因组 | 未本地部署（017 起改用远程 range 提取） |
| 可真跑部分 | **无** |
| 结论性质 | 全篇为文档口径与命令模板，未做任何真跑 |

**未真跑的直接原因**：DeepVariant 以容器分发（`docker pull google/deepvariant:1.6.1`，镜像数 GB），需要（a）运行中的容器守护进程、（b）去重标记且排序索引的 BAM/CRAM、（c）本地参考 FASTA 及其索引。本环境三项**均不满足**，且模型推理对 WGS 需约 4–6 h CPU / 64 GB 内存。因此本篇为严谨的概念稿，不做任何伪实测。

## 成分拆解

### 1. 治理原则

DeepVariant 用训练好的卷积神经网络替代参数化的 HMM（隐马尔可夫模型，即按状态转移对序列做概率建模的经典统计模型）/贝叶斯基因型判定器，把 pileup 图像分类为 hom-ref / het / hom-alt。两个后果驱动所有下游决策：

1. **事后没有可手工调参的统计过滤器。** CNN 已经输出一个校准过的 FILTER 列（可信变异为 `PASS`，被判为纯合参考的位点为 `RefCall`）。在 DeepVariant 输出之上叠加 GATK 硬过滤（QD/FS/MQ/SOR 阈值）或 VQSR，删除的是真阳性而非假阳性——那些注释在 VCF 里根本不存在。调用后的处理只限于 QUAL/GQ 阈值化、规范化与区域限制。
2. **网络是从原始碱基质量中学到错误模型的**，因此在上游跑 BQSR 会消耗时间并**略微降低** DeepVariant 的准确性。DeepVariant 自身的指导是跳过 BQSR。输入要求只是一份排序、索引、去重标记的 BAM/CRAM——没有更多。

DeepVariant 只做胚系变异。体细胞调用用 DeepSomatic；二倍体基因型类别无法表示亚克隆等位比例。

### 2. 工作流程（三阶段）

1. **`make_examples`**（CPU 密集，运行时瓶颈）扫描 BAM，在「非参考支持通过宽松的、为召回率调优的筛选」的候选位点处，把每个候选渲染成多通道 pileup 图像，写入分片的 TFRecords。行是读段、列是参考位置；通道编码读段碱基身份、碱基质量、比对质量、链方向、该读段是否支持候选等位、以及该碱基是否与参考不同。Illumina 模型额外加插入片段大小通道；长读段模型加单倍型通道。具体张量维度随版本变化——任何公开的数字都应视为示意。用 `--num_shards` 并行化。
2. **`call_variants`** 用训练好的 Inception 系列 CNN 对每个样本做推理，输出三类基因型似然。这是**唯一**可 GPU 加速的阶段。
3. **`postprocess_variants`** 排序 CNN 输出、解析多等位、把似然转换为 VCF/gVCF。

这种基于图像的设计正是 DeepVariant 在 indel 与困难上下文（同聚物、串联重复、低复杂度区域）上胜过参数化调用器的原因：CNN 学到了 pileup 几何中启发式过滤器看不到的视觉模式。模型是**平台特异的**，因为测序仪的错误模式（Illumina 的替换、ONT 的同聚物 indel）在图像上看起来不同，而每个模型学到的是其训练平台的伪影分布。

### 3. 模型选择

`--model_type` 是承重参数：用错模型会**静默**降低准确性，因为 CNN 期望在 pileup 中看到平台特异的错误模式，而它**不会报错**。

| `--model_type` | 适用 | 训练于 | 失效/降级场景 |
|---|---|---|---|
| `WGS` | Illumina 短读长 WGS | 30–50× PCR-free Illumina | 不加 `--regions` 用于外显子、用于长读长、用于 PCR 扩增子数据 |
| `WES` | Illumina 外显子/靶向 | 捕获外显子 | 不提供 `--regions` BED 就运行（会浪费数小时扫描靶外基因组） |
| `PACBIO` | PacBio HiFi（CCS） | HiFi，每读段 Q30+ | 用于 CLR 读段（Q10–15 的错误谱模型从未见过） |
| `ONT_R104` | ONT R10.4+ 化学 | R10.4 simplex/duplex | 用于 R9.4 数据（改用 Clair3 的 R9.4 模型）；精度仍低于 HiFi |
| `HYBRID_PACBIO_ILLUMINA` | 同时有 HiFi 与 Illumina 的样本 | 混合 HiFi + Illumina | 只有单一平台时 |

### 4. DeepVariant vs GATK vs DRAGEN

- **DeepVariant**——开源工具中 indel 准确性最好、困难区域/长读长表现最好；换模型即可跨平台泛化；无需调节过滤器。indel、困难区域与长读长的默认选择。
- **GATK HaplotypeCaller**——每个参数可审计、成熟的联合调用与参考置信度平方化、有监管先例。需要 GenomicsDB 扩展的超大队列、或已在 GATK 上验证过的临床流程宜选它。
- **DRAGEN**——FPGA 加速，30× 基因组约 20–25 分钟，在难以比对的基准中胜出；有硬件或云时按吞吐量优选。

选择取决于队列规模、平台、可审计性与吞吐量，而非单纯准确性。

### 5. 安装与运行（命令模板）

```bash
docker pull google/deepvariant:1.6.1              # GPU 版：google/deepvariant:1.6.1-gpu
singularity pull docker://google/deepvariant:1.6.1  # Singularity 替代

# 一次性运行（务必同时产出 gVCF，便于后续用 GLnexus 做队列联合调用而不重跑 DeepVariant）
docker run -v "${PWD}:/input" -v "${PWD}/output:/output" \
    google/deepvariant:1.6.1 \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=/input/reference.fa \
    --reads=/input/sample.bam \
    --output_vcf=/output/sample.vcf.gz \
    --output_gvcf=/output/sample.g.vcf.gz \
    --num_shards=16
```

外显子/靶向调用必须加 `--regions`（否则 WES 模型会扫描靶外基因组数小时）：

```bash
docker run -v "${PWD}:/data" google/deepvariant:1.6.1 \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WES --ref=/data/reference.fa --reads=/data/exome.bam \
    --regions=/data/targets.bed --output_vcf=/data/exome.vcf.gz --num_shards=8
```

三阶段分离（用于自定义分片、断点续跑、混用 CPU/GPU 节点）：

```bash
docker run -v "${PWD}:/data" google/deepvariant:1.6.1 \
    /opt/deepvariant/bin/make_examples --mode calling \
    --ref /data/reference.fa --reads /data/sample.bam \
    --examples /data/examples.tfrecord.gz --gvcf /data/gvcf.tfrecord.gz

docker run -v "${PWD}:/data" google/deepvariant:1.6.1 \
    /opt/deepvariant/bin/call_variants \
    --outfile /data/call_variants.tfrecord.gz \
    --examples /data/examples.tfrecord.gz --checkpoint /opt/models/wgs

docker run -v "${PWD}:/data" google/deepvariant:1.6.1 \
    /opt/deepvariant/bin/postprocess_variants \
    --ref /data/reference.fa --infile /data/call_variants.tfrecord.gz \
    --outfile /data/output.vcf.gz --gvcf_outfile /data/output.g.vcf.gz \
    --nonvariant_site_tfrecord_path /data/gvcf.tfrecord.gz
```

### 6. GPU 加速

GPU 加速只惠及 `call_variants`（CNN 推理）；`make_examples` 与 `postprocess_variants` 是 CPU 密集的，随 `--num_shards` 扩展。对大型队列，在 CPU 节点上跨样本并行往往比排队等 GPU 更划算。

```bash
docker run --gpus all -v "${PWD}:/data" google/deepvariant:1.6.1-gpu \
    /opt/deepvariant/bin/run_deepvariant --model_type=WGS \
    --ref=/data/reference.fa --reads=/data/sample.bam \
    --output_vcf=/data/output.vcf.gz --num_shards=16
```

### 7. DeepTrio（家系 / 三人组调用）

DeepTrio 把 pileup 图像扩展为同时覆盖先证者与双亲，使 CNN 能学到遗传上下文并直接调用新发变异。这优于朴素的三人组相减——后者的「新发集合」被逐样本独立错误产生的假阳性主导。适用于家系研究、孟德尔一致性工作与新发变异发现。

```bash
docker run -v "${PWD}:/data" google/deepvariant:deeptrio-1.6.1 \
    /opt/deepvariant/bin/run_deeptrio --model_type=WGS \
    --ref=/data/reference.fa \
    --reads_child=/data/child.bam --reads_parent1=/data/father.bam --reads_parent2=/data/mother.bam \
    --sample_name_child=CHILD --sample_name_parent1=FATHER --sample_name_parent2=MOTHER \
    --output_vcf_child=/data/child.vcf.gz --output_vcf_parent1=/data/father.vcf.gz \
    --output_vcf_parent2=/data/mother.vcf.gz \
    --output_gvcf_child=/data/child.g.vcf.gz --output_gvcf_parent1=/data/father.g.vcf.gz \
    --output_gvcf_parent2=/data/mother.g.vcf.gz --num_shards=16
```

三个逐样本 gVCF 再用 GLnexus 合并为一个三人组 VCF；联合上下文正是支撑孟德尔违规与新发率分析的东西。

### 8. 用 GLnexus 做联合调用

DeepVariant 的 gVCF 要用 **GLnexus** 联合基因分型，**不是** GATK GenotypeGVCFs——GLnexus 跨逐样本 gVCF 做等位置一（allele unification），并随样本加入增量增长，避免全队列重处理。

```bash
docker run -v "${PWD}:/data" quay.io/mlin/glnexus:v1.4.1 \
    /usr/local/bin/glnexus_cli --config DeepVariantWGS /data/*.g.vcf.gz \
    | bcftools view - -Oz -o cohort.vcf.gz
```

| GLnexus `--config` | 用途 | 说明 |
|---|---|---|
| `DeepVariantWGS` | Illumina WGS gVCF | 多数 WGS 队列的默认 |
| `DeepVariantWES` | Illumina 外显子 gVCF | 针对更高深度、更窄区域调优 |
| `DeepVariant_unfiltered` | 保留全部变异位点 | 研究探索；假阳性更多，对三人组/新发分析有用（`RefCall` 位点在此有意义） |

### 9. 输出与质控

DeepVariant 输出已被 CNN 过滤（FILTER 列为 `PASS` / `RefCall`）。**不要**应用 GATK 硬过滤或 VQSR。合法的后处理只有 QUAL/GQ 阈值化、规范化与区域限制：

```bash
bcftools stats output.vcf.gz | grep TSTV          # Ti/Tv：WGS 期望 2.0-2.1，WES 3.0-3.3
bcftools view -i 'QUAL>20 && FMT/GQ>20' output.vcf.gz -Oz -o filtered.vcf.gz
```

### 10. 基准测试与 GIAB 循环性告诫

对 GIAB 真值集做基准时，应使用单倍型感知的比较器（hap.py + vcfeval），限制在置信区域 BED 内，并按区域难度分层：

```bash
docker run -v "${PWD}:/data" jmcdani20/hap.py:latest \
    /opt/hap.py/bin/hap.py /data/HG002_GRCh38_truth.vcf.gz /data/deepvariant_output.vcf.gz \
    -f /data/HG002_confident.bed -r /data/reference.fa -o /data/benchmark \
    --engine=vcfeval --threads 16
```

**承重的告诫**：DeepVariant 在 GIAB 真值集上训练（主要是 HG001），随后又习惯性地在 GIAB 样本上做基准。当训练集与测试集都来自 HG001–HG007 时，一个「F1 = 0.999」的宣传数字，有相当一部分是在度量对真值集特异性的记忆，而非泛化能力。诚实的做法是：看重留出样本的表现（在 HG001/3/4/5/6/7 上训练、在 HG002 上测试——precisionFDA V2 正是通过对半盲的父母 HG003/HG004 打分来做到这一点）；报告困难区域与 CMRG 分层，而不是一个全基因组数字；在临床部署前，在与人群匹配的、已充分表征的材料上验证，而不是相信一个已发表的 GIAB F1。**一个只报告全局 F1、既无分层也无留出或非 GIAB 样本的基准，不具备决策价值。**

### 11. 与其他调用器的近似准确性

GIAB HG002/HG003/HG004（GRCh38）上的近似 F1；具体数值随样本、覆盖度与版本变化。在简单 SNP 上每个现代调用器都超过 F1 0.999，因此有决策价值的差距在 indel 与困难区域。

| 调用器 | SNP F1 | Indel F1 | 速度（30× WGS） | 说明 |
|---|---|---|---|---|
| DeepVariant | ~0.999 | ~0.993 | ~4–6 h CPU，~1–2 h GPU | 开源工具中 indel 准确性最高；无 GPU 时较慢 |
| GATK HaplotypeCaller | ~0.999 | ~0.989 | ~4–8 h CPU | 可审计；联合调用生态成熟 |
| Strelka2 | ~0.998 | ~0.960 | ~1–2 h CPU | 快；已不再积极维护 |
| Clair3 | ~0.998 | ~0.980 | ~8 h（50× ONT） | 长读段表现强；开发中 |

### 12. 资源需求

| 数据 | 内存 | CPU 时间 | GPU 时间 | 说明 |
|---|---|---|---|---|
| WGS 30× | 64 GB | ~4–6 h | ~1–2 h | `--num_shards` 使 make_examples 线性扩展 |
| WES | 32 GB | ~30 min | ~10 min | 靶向区域更小 |
| PacBio HiFi 30× | 64 GB | ~3–5 h | ~1–2 h | 读段更少但更长 |
| ONT 50× | 64 GB | ~6–8 h | ~2–3 h | 错误率更高 → 候选位点更多 |

## 未覆盖（诚实标注）

本篇**未做任何真跑**。本环境缺少三项必要条件：运行中的容器守护进程（`docker info` 失败）、去重标记的 BAM/CRAM、本地参考 FASTA 与索引。因此以下内容仅为文档口径与命令模板，未经验证：

- 容器拉取与 `--model_type` 令牌的实际枚举（新版本会新增模型并重命名镜像标签，须以所用容器的 `--helpfull` 为准）。
- `run_deepvariant` 与三阶段流程的实际产出、耗时与张量维度。
- GPU 加速比、DeepTrio 的新发调用、GLnexus 合并。
- Ti/Tv 与 F1 的实际测量。

## 实践要点

- **模型必须匹配测序仪**：`--model_type` 用错不会报错，只会静默降级。
- **不要在 DeepVariant 输出上叠加 GATK 硬过滤或 VQSR**：它没有 QD/FS/MQ 注释，CNN 已过滤；只做 QUAL/GQ 阈值化。
- **上游跳过 BQSR**：网络是从原始质量中学的，BQSR 会略微降低准确性。
- **WES 必须给 `--regions`**，否则模型会扫描靶外基因组数小时。
- **GPU 只加速 call_variants**；CPU 阶段靠 `--num_shards`。
- **队列合并用 GLnexus，不用 GenotypeGVCFs**：DeepVariant 的 gVCF 不是 GATK 参考置信度 gVCF。
- **三人组新发变异用 DeepTrio**，不要做朴素的三人组相减。
- **读 GIAB F1 要带留出于与分层**：全局 F1 部分度量的是对真值集的记忆。

## 小结

deepvariant 的机制核心是「把 pileup 渲染成图像、用 CNN 分类基因型」，由此带来两条刚性后果：输出已由 CNN 校准因而**禁止叠加 GATK 过滤器**，以及上游**应跳过 BQSR**。本篇因环境缺容器守护进程与 BAM，未能做任何真跑；为避免伪实测，全篇按文档口径完整记录了机制、模型选择、命令模板、资源需求与 GIAB 循环性告诫，并明确标注未验证范围。

（本篇无真跑产物，故未建素材目录。）
