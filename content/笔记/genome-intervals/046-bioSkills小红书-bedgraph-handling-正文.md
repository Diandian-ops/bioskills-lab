<!-- META
用途: 小红书帖子正文（bioSkills 046 bedgraph-handling）
标题建议: reads砍到1/4，峰高也砍到1/4
封面卡文案: 2 Mb 模拟染色体 + 4 个峰，bedGraph 转 bigWig，27,427 行对账零误差
-->

# bedGraph 轨道：峰高差 4 倍，可能只是测序深度差 4 倍

（编号 046 · 真跑日期 2026-09-03：WSL Ubuntu + bedtools v2.31.1 + bedGraphToBigWig v2.10 + pyBigWig 0.3.25；模拟 chr1 长 2,000,000 bp，样本 A 20,591 条 reads、样本 B 为其 25% 子样本）

**要点一｜同一生物学，raw 轨道峰高差约 4 倍**
把同批 reads 随机砍掉 75% 当样本 B（20,591 → 5,155 条，深度比 3.99x），两条 bedGraph 在浏览器里峰形一模一样，高度却差约 4 倍（A 峰区最高 78x vs B 最高 24x）。bedGraph 第 4 列本质是测序深度的函数：不归一化就跨样本比高低，比的全是库深。

**要点二｜一行 -scale 把两条轨道拉回同一标尺**
bedtools genomecov 的 `-scale 1000000/reads数` 就是最朴素的 RPM：A 乘 48.564907、B 乘 193.986421（比值恰好等于深度比 3.994），缩放后两轨道重合。但注意 CPM/RPKM/BPM/RPGC 这类库深归一化全都预设「总信号守恒」——遇到全局水平变化（比如组蛋白修饰敲低）会把真变化摊平成假无变化，那条路只有上机前的 spike-in 能救。

**要点三｜转 bigWig 的四条红线，v2.10 全部显式报错**
拿违规输入逐个试了 bedGraphToBigWig：未排序 →「is not sorted at line 2」；chrom.sizes 偏短 →「End coordinate 1000122 bigger than chr1 size of 1000000」；两条轨道拼接出重叠 →「overlapping regions in bedGraph line 10」，用 bedtools merge -c 4 -o max 折叠（34,862 → 5,944 行）后一次通过；染色体改名 chr1→1 →「1 is not found in chromosome sizes file」。转换成功的轨道做无损对账：27,427 行逐行回读，max delta = 0，文件从 613 KB 压到 234 KB 还附带索引。

**结论**：轨道比较前先裁决归一化合法性，chrom.sizes 永远从比对用的那个 FASTA 生成（samtools faidx + cut -f1,2），排序用 LC_COLLATE=C，转换后用 pyBigWig 做一次逐行对账再分发。
#生信 #生物信息学 #bedGraph #bigWig #ChIPseq #bioSkills
