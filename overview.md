# 010 structural-alignment 真实试用 — 完成

## 已完成
按 bioSkills alignment 家族「真实试用」路线，完成第 10 个 skill `structural-alignment` 的复现、拆解、出图、站点接入与提交。

### 真实复现（都真跑过）
- **环境**：Biopython 1.88（managed venv）；Foldseek 10.941cd33（conda env `foldseek`）；TM-align 20220412（源码编译，macOS 需把 `#include <malloc.h>` 改为 `<malloc/malloc.h>`）。
- **数据**：RCSB 下载 `content/素材/010-structural-alignment/pdbs/{1ubq,1ubi,1fmb}.pdb`（1ubq/1ubi = 泛素同源；1fmb = 跨折叠对照）。
- **P1 Bio.PDB.Superimposer**：1ubq↔1ubi RMSD = 0.0908 Å（76 CA，几乎重合）。
- **P2 TMalign 两两**：1ubq↔1ubi TM-score = 0.9991 / RMSD 0.09（>0.5 同折叠）；1ubq↔1fmb TM-score = 0.4079 / RMSD 3.41（<0.5 不同折叠）。
- **P3 Foldseek**：easy-cluster 把 1ubq+1ubi 聚为一簇、1fmb 单独成簇，与 TMalign 结论一致；easy-search 默认仅返回 1ubq→1ubq/1ubi（过滤掉 TM<0.5 的 1fmb）。
- **配图**：`010-fig.png` 两两 TM-score 矩阵热图（红=相似度高 / 绿=低，契合中国习惯）。

### 文档产出
- `content/笔记/010-bioSkills真实试用-structural-alignment.md` — DEEP DIVE 07 深度笔记（功能定位 / 成分拆解 / 严格复现 / 坑实测 / 实践要点）。
- `content/笔记/010-bioSkills小红书-structural-alignment-正文.md` — 小红书正文（≤20 字标题，4 段结构），gate_lint ERROR=0。
- `content/素材/010-structural-alignment/` — `run_structural_alignment.py`、`make_fig.py`、`structural_results.json`、3 个真实 PDB、编译好的 TMalign 二进制。
- 站点已重建：`output/bioSkills-site/alignment/structural-alignment.html`，并接入左侧栏 / 总览「已完成 7 个深度试用」。

### 关键发现（需用户注意）
1. **SKILL.md 文档错误**：`--format-output` 里的 `tmscore` 在 Foldseek 实际不可用，正确为 `alntmscore`/`qtmscore`/`ttmscore`。
2. **Foldseek 默认检索口径**：默认 alignment-type（3Di+AA）会过滤掉 TM<0.5 的跨折叠候选，低相似度探测需 `--alignment-type 1` 或 easy-cluster 复核。
3. **未实测**：US-align / DALI / Foldmason / Foldseek-Multimer / pLDDT 掩码 / PyMOL·ChimeraX 可视化——仅按文档记录，未执行。

### 下一步
- 如需补全说明书，可创建自动化「每日/定期」巡检 bioSkills 文档与实际二进制行为差异。
- 本 skill（library 级，非 agent_created）的 `tmscore` 文档错误属 P2 级小修正，已在本报告标注；是否修改原 skill 需用户确认（未自动改）。

提交：待用户确认后 push（本次未执行 git push，避免未授权推送）。
