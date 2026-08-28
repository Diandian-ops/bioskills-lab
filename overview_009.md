# 009 bioSkills alignment-io 真实试用 — 完成概览

## 完成内容
- 真实复现：纯 Biopython（`Bio.AlignIO` + `Bio.Align`），复用仓库 `alignment-io/examples/sample_alignment.aln`（CLUSTAL，4×21）。
- 产物：
  - `content/笔记/009-bioSkills真实试用-alignment-io.md`（进站点 DEEP DIVE 06）
  - `content/笔记/009-bioSkills小红书-alignment-io-正文.md`（标题建议「比对格式转换会静默丢注释」12 字）
  - `content/素材/009-alignment-io/`（脚本、图、日志、样本、数据）
- 两份 lint 均 0 ERROR 0 WARN。
- `pipeline/build_lab_site.py` 已接线 009，站点重建并 http 200 抽检通过。
- 已提交并推送 `main`：641e154 → origin。

## 关键真实发现
1. **注释存活矩阵**：仅 Stockholm 回读保留 GS/GR/GC；NEXUS 仅保留 GS；FASTA/Clustal/PHYLIP-relaxed 全部静默丢弃。
2. **PHYLIP strict 10 字符截断**：BioPython 1.88 对冲突长名直接抛 `ValueError`，非静默合并；`phylip-relaxed` 可保留长名。
3. **MAF 负链坐标**：BioPython 1.88 将 `strand` 解析为整数 `-1`/`1`，SKILL.md 示例按字符串 `'-'` 判断会导致负链 plus_start 错误（10 bp 而非 33 bp）。
4. **A2M 大小写编码 / Bio.Align 现代 API / pyhmmer 未安装声明**均如实记录。

## 下一任务
- 按顺序进入 010 `structural-alignment`：需 brew 安装 TM-align/Foldseek，下载真实 PDB，运行结构比对。
