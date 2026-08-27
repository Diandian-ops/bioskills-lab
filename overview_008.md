# 008｜bioSkills msa-parsing 真实试用 — 概览

按既定模式（`004/005/006/007`）完成 alignment 第 5 个 skill **msa-parsing** 的真实试用，已推送 `origin/main`（`b1e9786`）。

## 做了什么
- **读 SKILL.md + 探示例**：msa-parsing 为纯 Biopython（AlignIO）薄封装，含保守性/gap/共识/过滤/坐标映射/Henikoff/Neff/MI-APC。
- **真实复现**：用仓库自带示例 `alignment-io/examples/sample_alignment.aln`（CLUSTAL，4×21）跑通全部函数；环境为 managed venv（Bio 1.88 / numpy 2.1.3）。
- **两份笔记**（gate_lint 均 0/0）：
  - 深拆版 `008-bioSkills真实试用-msa-parsing.md` → 进站点 DEEP DIVE 05
  - 小红书正文 `008-bioSkills小红书-msa-parsing-正文.md`（标题「小样本上多序列比对加权会退化」14 字）
- **出图**：`008-fig.png`（逐列 gap%/保守% 曲线）+ 复现脚本/日志/样本，落 `content/素材/008-msa-parsing/`。
- **站点重建**：`build_lab_site.py` 接线 008（ALIGN_SUBS / 侧栏 active+skip / 总览 DONE 008 / 已完成 5 个），重建后 http 200 抽检全过。
- **提交推送**：`b1e9786`（8 文件，668 行）。`output/` 仍 gitignore，站点本地生成。

## 关键发现（如实写入笔记）
- 小样本下加权指标退化：Henikoff 全 0.25、Neff=0.25（4 条聚 1 簇，Neff/L=0.012）、MI 全 0 —— 恰好印证 skill 标注的适用边界。
- 用 007 真实 102×3 MAFFT 对齐二次验证同一套函数：Henikoff 0.31/0.31/0.38（非均匀）、50 gappy 列 → 在真实 gap 富集数据上正常分化。

## 出图（已完成，附说明）
- **小红书卡片已出图**：用 `008-bioSkills小红书-msa-parsing-正文.md`（正文版，不引用图片）跑通 trial2xhs，导出 2 张 PNG → `output/xhs-cards/008-bioSkills小红书-msa-parsing-正文/notebook-standalone-hd/`。
- **真实试用版出图受阻（md2card 应用层 bug，非本机/网络问题）**：`008-bioSkills真实试用-msa-parsing.md` 引用 `008-fig.png`，md2card 在绑定本地图片后，导出 preflight 等待图片 ready 永远不置位，提交按钮 60s 始终 `disabled`（无 console 报错、非单纯慢）。正文版不引用图片故顺利。小红书卡片用正文版出图更契合发布场景；若需真实试用版含图卡片，须修 md2card 图片 preflight 或绕过绑图。

## 剩余 alignment 候选
- `alignment-io`（读写字）、`structural-alignment`（需 Foldseek/TM-align，成本更高）。
