<!--
META
标题建议: DeepVariant 过滤与 BQSR 处理
/ META
-->

# DeepVariant 两条铁律：别叠 GATK 过滤，别做 BQSR

（编号 027 · 本篇为概念稿：本机 docker 守护进程未运行、无 BAM 与本地参考，未做真跑）

**要点一｜它的输出已经过滤过了**
DeepVariant 把每个候选位点的 pileup 渲染成多通道图像（行是读段、列是参考位置，通道编码碱基、质量、比对质量、链方向、是否支持候选等位、是否与参考不同），再用训练好的 CNN 分类成 hom-ref / het / hom-alt。

所以它的 FILTER 列已经是 CNN 的校准结果（`PASS` 或 `RefCall`）。**在它上面叠 GATK 硬过滤（QD/FS/MQ/SOR）或 VQSR，删掉的是真阳性**——因为那些注释在 DeepVariant 的 VCF 里根本不存在。合法的后处理只有三样：QUAL/GQ 阈值化、规范化、区域限制。

**要点二｜上游要跳过 BQSR**
网络是从**原始**碱基质量里学到错误模型的。跑 BQSR 不仅浪费时间，还会**略微降低** DeepVariant 的准确性。输入要求只有一份排序、索引、去重标记的 BAM/CRAM——没有更多。

另外模型必须匹配测序仪（`WGS` / `WES` / `PACBIO` / `ONT_R104` / `HYBRID_PACBIO_ILLUMINA`）。用错模型**不会报错**，只会静默降级，因为 CNN 期望看到的是该平台特有的错误模式。外显子还必须给 `--regions`，否则模型会花数小时扫描靶外基因组。

**要点三｜队列合并用 GLnexus，别用 GenotypeGVCFs**
DeepVariant 的 gVCF 不是 GATK 的参考置信度 gVCF，走 GenotypeGVCFs 会失败。正确路径是 GLnexus（预设 `DeepVariantWGS` / `DeepVariantWES`），它随样本加入增量增长，不必全队列重跑。参考吞吐：2504 样本的 chr22 合并，GLnexus 0.84 h vs GATK 6.83 h。

GPU 只加速 `call_variants`（CNN 推理）这一阶段；`make_examples` 与 `postprocess_variants` 是 CPU 密集的，靠 `--num_shards` 扩展。

顺带一条读文献的提醒：DeepVariant 在 GIAB 上训练，又习惯在 GIAB 上做基准，所以那个「F1 = 0.999」有一部分是在度量对真值集的记忆。要看留出样本（在 HG002 上测）和困难区域/CMRG 分层，而不是一个全局数字。

**结论**：DeepVariant 的价值在 indel 与困难区域；代价是必须改掉从 GATK 带来的两个肌肉记忆——不要后置过滤，不要前置 BQSR。
#生信 #生物信息学 #DeepVariant #变异检出 #bioSkills
