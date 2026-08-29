# VCF 拼接合并：bcftools 组合操作的三个真坑

（编号 019 · 真实数据：1000G Phase3 chr22 / 2504 样本 / 5431 位点）

**要点一｜先规范，再组合**
`bcftools merge / isec / concat` 都以「染色体+位置+REF+ALT」四元组为键。一个未左对齐的 indel、未拆的多等位、未分解的 MNP，会被当成「不同变异」静默错叠。`isec` 默认要求 ALT 完全匹配。实测：原始 5431 条、norm 后 5471 条，isec 显示 norm 多出 80 条「私有」记录——多等位被拆开后，同一位点变成了不同的 ALT 元组。组合前务必 `bcftools norm -m-any -f ref`。

**要点二｜merge 不是联合基因型**
单样本 VCF 合并时，某样本「没这个位点」和「纯合参考」分不清。`bcftools merge` 默认填 `./.`，加 `-0` 会填 `0/0`——都是猜测，因为 merge 对该位点没有证据。真要区分「确证 hom-ref」与「无数据」，得用 gVCF 联合基因型（GenotypeGVCFs）。不要拿 merge 当 joint calling。

**要点三｜子集后 AC/AN 会失真**
`bcftools view -s` 默认会重算 AC/AN，但一旦用 `-I` 或别的方式绕过，统计量就停在原始样本数。实测：抽 5 个样本后，首记录仍显示 `AN=5008`（原 2504 样本），实际应为 `10`（5×2）。一行修复：`bcftools +fill-tags -- -t AC,AN,AF`。子集后下游 allele frequency 会被旧基数污染，这一步不能省。

**结论**：vcf-manipulation 的核心是「统一表示 + 统一样本名/ contig」后再组合；规范化与 fill-tags 是避免静默错误的前置条件。
