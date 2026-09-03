<!-- META
用途: 小红书帖子正文（bioSkills 038 quality-reports）
标题建议: FastQC 质控红绿灯的正确读法
封面卡文案: 3 份模拟 FASTQ × 10 个 FastQC 模块，30 次红绿灯判定全记录
-->

# FastQC 质控报告：红绿灯之外还要读图

（编号 038 · 真跑日期 2026-09-03：WSL Ubuntu + FastQC v0.12.1 + MultiQC v1.35；3 份模拟 FASTQ，各 20000 条 reads × 100 bp）

**要点一｜同一套流程，三种形态，30 次判定 = 30 次灯亮对**
S1_good 十个模块全 pass；S2_degraded 的 3' 端均值 Q 从 40.8 衰减到 7.0，per-base quality 直接 fail；S3_adapter 末端 adapter 含量 19.51%（fail 阈值 10%），overrepresented 同时 warn。10 个模块 × 3 份样本共 30 次判定，与三种数据形态一一对应——受控输入下灯的语义完全可预测。

**要点二｜0.365% 的 duplication，全部来自 74 条同源序列**
S3 里埋了 80 条完全相同的 adapter-dimer reads，FastQC 抓到其中 74 条（占总 reads 0.37%），样本去重率 100.0% → 99.635%。这正好演示「duplication 是读段级指标」：同样形态在真实数据里可能是 PCR 偏差，也可能是高表达分子的真实丰度，去不去重要看复杂度证据，不是见重就去。

**要点三｜FastQC 会猜错质量编码，读图前先查 Encoding 字段**
第一版 S1 的质量平直在 Q35–39 且无低值字符，FastQC 判成 Illumina 1.5（Phred+64），均值 Q 整体平移 31，实测 Q(first10) 35.81 → 7.26，模块直接 fail——数据本身没变。加入 Q30 低值抖动复测，三样本 Encoding 统一为 Sanger / Illumina 1.9，Q(first10) 7.26 → 35.81。

**结论**：红绿灯按随机 WGS 文库校准，受控样本能全亮对，换实验类型就要先读图再读灯。定位污染最省力的证据链是两模块交叉印证：adapter 曲线 + overrepresented 序列（实测自动标注来源 = TruSeq Adapter, Index 9）。
#生信 #生物信息学 #FastQC #MultiQC #质控 #bioSkills
