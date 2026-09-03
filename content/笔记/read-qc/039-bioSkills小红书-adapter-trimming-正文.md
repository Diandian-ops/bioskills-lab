<!--
META
编号: 039
模块: read-qc / adapter-trimming
真跑日期: 2026-09-03
环境: WSL Ubuntu · cutadapt 5.2 vs Trimmomatic 0.41
标题建议: 接头修剪差一个词，R2 成对消失
/META
-->

# 差一个词，每个 read-through 对丢一条 R2

（编号 039 · 真跑 2026-09-03 · WSL Ubuntu · cutadapt 5.2 vs Trimmomatic 0.41 · 合成梯度数据每档 20,000 对，接头读段占比 5.12% / 20.14% / 40.48% 三档）

**封面卡｜一个词的差别：keepBothReads vs true —— 40% 档 R2 差 8,095 条**

**发现一｜SKILL.md 字面参数串：丢掉的 R2 为 1,025 / 4,029 / 8,095 条，精确等于 read-through 对数**

Trimmomatic 字面口径 ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads 跑三档梯度，成对保留率 94.88% / 79.86% / 59.52%，单 mate 输出 1,025 / 4,029 / 8,095 条 —— 与造数据真值里的 read-through 对数逐一相等，即每个短插入片段读对恰好丢掉 R2。根因：末位参数由 Boolean.parseBoolean 解析，单词 keepBothReads 被解析为 false，palindrome 模式检出 read-through 后就把 R2 当冗余丢弃。末位改为 true 重跑：三档丢失均为 0 条，去除碱基量 39,402 / 153,854 / 309,810 bp，与真值完全一致。

**发现二｜cutadapt 默认参数召回 98.3%-98.4%，缺口集中在 5–9 bp 短残留**

三档召回率 98.34% / 98.32% / 98.44%，逐箱对比：5–9 bp 残留箱召回 93.40%（4,205/4,502 条），10–15 bp 箱起回升到 99.59% 以上；对照 Trimmomatic palindrome 模式四个残留长度箱均 100.00%（它从 R1/R2 重叠反推，不依赖逐条匹配）。另一面是默认 -O 3 在干净读段上的误剪：825 / 646 / 456 条，占干净读段 2.17% / 2.02% / 1.92%，约合 0.06 bp/read 的碱基代价。

**发现三｜Trimmomatic 0.41 参数个数对比：只认 6 参数形态，4/5 参数直接报错**

参数个数对比实测：4 参数形态（…:10:true）与 5 参数形态（…:10:false）均抛 NumberFormatException（第 5 位按整数解析）；完整 6 参数 …:2:30:10:2:true 才能运行。换工具或换版本时，ILLUMINACLIP 参数串需要整体核对，不是只改布尔位。

**结论**：正确配置下两种工具都能把接头剪干净（三档丢失 0 条 vs 0 条），实测差异集中在两处——参数写法（字面词 vs 显式布尔值，直接决定 R2 去留）与 5–9 bp 短残留（93.40% vs 100.00%）。把 R2 丢失定位到「每个 read-through 对丢一条」，靠的是用 12 万行逐读段真值逐 read_id 对齐，而不是看工具自带报告。

#生信 #生物信息学 #cutadapt #Trimmomatic #接头去除 #bioSkills
