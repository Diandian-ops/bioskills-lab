<!--
META
编号: 045
模块: genome-intervals / bed-file-basics
真跑日期: 2026-09-03
环境: WSL Ubuntu · bedtools v2.31.1
标题建议: BED文件为什么会差1个碱基？真跑实测
/META
-->

# start 是裸整数：BED 坐标约定真跑校验

（编号 045 · 真跑 2026-09-03 · WSL Ubuntu · bedtools v2.31.1 · 模拟 chr1 全长 2,000,000 bp：40 基因 / 215 外显子 / 30 CpG 岛 / 10 narrowPeak / 8 条 BED12）

**封面卡｜start 是裸整数，约定不存在文件里——BED 比 GTF/VCF 少的那个 1，全靠人记**

**发现一｜坐标换算只动 start：BED start = 1-based 起点减 1（单位 bp），end 数值不变**

BED 是 0-based 半开区间，GTF/VCF/SAM 是 1-based 闭区间，转换规则只有一句：start 减 1，end 不变。实测三层验证——215 个外显子全部经 `bedtools getfasta` 提取，序列长度与 end−start 逐条相等（215/215，首条 30,739−30,450=289 bp）；构造 1 bp 地标 GFF `1000-1000` 转 BED 得 `999 1000`，faidx、getfasta、VCF REF 三条途径取到的碱基一致（T）；30 条 VCF 变异转 BED、6 条 SAM 比对经 `bamtobed`，start = POS−1 全部成立。经典错误是「两端都减 1」，会把每个基因体做短 1 bp 且不报任何错。

**发现二｜三类静默失败真跑出两种，一种会报错**

- 染色体名不匹配（`chr1` vs `1`）：intersect 输出 0 行、退出码 0——完全合法的空结果，看起来像「真的没有重叠」。bedtools 2.31.1 只在 stderr 打了 naming convention 警告，不能依赖它兜底。
- CRLF 行尾：Windows 经手的文件 `cat -A` 实测行尾 `^M$`，`\r` 粘在最后一列上（实测 5/5 行），`sed 's/\r$//'` 一行修复。
- `-sorted` 遇乱序输入：直接报错退出（rc=1，Error: Sorted input specified...），去掉 `-sorted` 同样输入正常出 9 行。有报错的路径容易定位；名字不匹配这种全程无报错，更难察觉。

**发现三｜merge 只折叠重叠区间；划分 + 补集恒等于染色体长**

CpG 岛 30 → merge 20 → merge -d 200 得 19；无重叠的基因区 40 → 40、外显子 215 → 215（merge 后计数不变本身就是「无重叠」的快速检验）。覆盖恒等式实测精确成立：merge 后 CpG 30,862 bp + 补集 1,969,138 bp = 2,000,000 bp；基因区 365,304 bp + 补集 1,646,696 bp = 2,000,000 bp——genome.txt 现场从同一 FASTA 生成（samtools faidx + cut -f1,2），slop 越界自动裁回染色体端点。另附一个实测错误案例：narrowPeak 的 summit 取第 10 列（peak 偏移，575,414+383=575,797 bp），首跑误用第 8 列 pValue，错误配置算出合法但错误的数值——BED6+4 无表头，列语义纯靠位置记忆。

**结论**：BED 的一切疑难都回到「约定不在数据里」——转换靠 start−1、验证靠 1 bp 地标往返、防静默失败靠比对染色体命名和现场生成 genome 文件。本次真跑 bedtools 11 个子命令，坐标守恒性质（215/215 长度、36/36 次 POS−1、两条 2,000,000 bp 恒等式）全部成立，脚本与数据已落盘可复现。

#生信 #生物信息学 #bedtools #BED格式 #基因组注释 #bioSkills
