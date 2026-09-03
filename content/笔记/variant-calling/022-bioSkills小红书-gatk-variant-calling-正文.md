# GATK 真跑：0 条变异，恰是最对的一次输出

（编号 022 · 真跑日期 2026-09-03：WSL Ubuntu + GATK 4.6.2.0；合成参考 4 条 contig × 3000 bp；模拟比对 BAM 覆盖约 50×）

**要点一｜同一个 BAM，标准模式 0 条，GVCF 却有 2946 条**
本次完整跑了一遍 GATK 链路：bioconda 装 GATK 4.6.2.0 → `samtools addreplacerg` 补读组 → HaplotypeCaller 两种模式各跑一次。结果：

| 模式 | 产出 | bcftools stats |
|---|---|---|
| 标准 HaplotypeCaller | raw.vcf.gz | 变异记录 0 条 |
| -ERC GVCF | raw.g.vcf.gz | 2946 条 hom-ref 区块，每条带 `<NON_REF>` 符号等位 |

这不是 bug，正好实测了 GVCF 参考置信度模型的含义：无变异输入下，标准模式沉默；GVCF 对**每个位置**照样输出「此处为纯合参考」的置信度（2942 个区块 DP>0，均值 35.9），这就是队列联合基因型能区分「确证参考」与「无数据」的物理基础。

**要点二｜0 条变异的根因：输入本来就没有突变（三重验证）**
查 011 的模拟命令发现 reads 是 `wgsim -e 0.02 -r 0 -R 0` 生成的——**-r 0 = 零突变**，全基因组只有 2% 测序错误。三个独立证据闭环：

- 跨工具对照：`bcftools call -mv` 同样输出 0 条；
- 逐列扫描：全基因组最高 ALT 支持度仅 5 reads / 18.5%，没有任何位点具备真突变形态；
- samtools stats 实测 error rate = 2.19%，与 -e 0.02 对上。

结论：散在的随机测序错误不构成等位证据，两个调用器都给出了「无变异输入」的正确答案。顺带修了一个统计坑：数 mpileup 非参考列必须同时排除 `.` 和 `,`（反向匹配），否则 11152 虚高，修正后是 4659 / 11962 列。

**要点三｜环境不再卡 GitHub**
2026-09-01 评估时卡在「GitHub 拿不到 GATK jar」，本次经 bioconda 通道 `conda create -n bio-gatk -c bioconda -c conda-forge gatk4 openjdk=17` 一次装成，原阻塞结论作废。独立环境顺便把 openjdk 钉在 17，不污染主环境。

**结论**：拿到调用器的空输出，先怀疑数据、再怀疑参数——用第二个工具交叉验证、看位点支持度形态、核 error rate，三步定位根因。本次 0 条变异恰好完整展示了「调用器在无变异输入下的两种正确行为」。
#生信 #生物信息学 #GATK #HaplotypeCaller #GVCF #bioSkills
