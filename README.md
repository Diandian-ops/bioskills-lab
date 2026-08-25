# 小红书 · 生信工程师工作笔记 · 内容仓库

> 定位：**一个生信工程师的工作笔记**。随手记录日常分析、踩坑、用的工具、对方法的吐槽。
> 作者：男。口吻保持中性、真实、第一人称「我」，不堆术语、只讲能跑通的。
> 账号名 / Bio 统一「工作笔记」调性，不漂移。

---

## 目录结构（源 / 产物分层）

设计原则：**源可再生不了，产物一条命令重建**。源（笔记 / 素材 / 库 / 构建脚本）永不进产物目录。

```
RedBook/
├── README.md                  # 本文件：总索引 + 命名规范（每加一篇在此追加一行）
├── content/                   # 【源】所有不可再生内容，单一事实源
│   ├── 笔记/                  # 正式成稿的笔记（按编号命名，md 驱动站点 + 小红书）
│   │   ├── 001-生信AI工作流笔记.md
│   │   ├── 002-bioSkills使用体验.md
│   │   ├── 003-bioSkills拆解01-alignment专题.md
│   │   ├── 003-bioSkills拆解01-alignment专题_dingcard.md   # 历史遗留，DingCard 已停用
│   │   ├── 004-bioSkills真实试用-pairwise-alignment.md     # 完整试用版
│   │   ├── 004-bioSkills小红书-pairwise-alignment.md       # 引流精简版
│   │   ├── 005-bioSkills真实试用-msa-statistics.md         # 完整试用版
│   │   └── 005-bioSkills小红书-msa-statistics.md           # 引流精简版
│   ├── 素材/                  # 图片按所属笔记分文件夹，互不打扰
│   │   ├── 002-bioSkills/
│   │   ├── 003-alignment/
│   │   ├── 004-pairwise/
│   │   └── 005-msa-statistics/
│   └── 库/                    # 参考仓库 / 外部资源（本地 clone）
│       ├── bioSkills/         # GPTomics/bioSkills (562 skills, 65 categories, 本地统计)
│       │   # 上游: https://github.com/GPTomics/bioSkills
│       │   # 用途: skills 拆解系列的源码对照与真跑验证
│       └── obsidian-wechat-converter/
├── pipeline/                  # 【源】构建脚本 + 真实复现实验（不可再生，但属"生产工具"）
│   ├── build_lab_site.py      # 站点生成器：读 content/笔记/*.md → output/bioSkills-site/
│   ├── make_cover.py          # 封面卡生成器（品牌标题卡）
│   ├── 004-pairwise-redo/     # pairwise-alignment 严格复现脚本
│   └── 005-msa-statistics/    # msa-statistics 严格复现脚本
└── output/                      # 【纯产物】一条命令可重建，删了不心疼
    └── bioSkills-site/        # 多文件静态站点（index/alignment/004/005 + assets 三图）
```

> ⚠️ **源 / 产物不可混居**：`build_lab_site.py` 重生成时会先 `rmtree(output/bioSkills-site)` 再重建。所有源都在 `content/` 与 `pipeline/`，绝不在 `output/` 内，否则会被误删。

---

## 命名规范

| 类型 | 规则 | 示例 |
|------|------|------|
| 笔记 | `编号-主题.md`，编号从 001 递增，跨系列连续 | `002-bioSkills使用体验.md` |
| 素材 | 按所属笔记建文件夹 `编号-主题/`，图片保留原名或加前缀 `cover-` | `004-pairwise/cover-xxx.png` |
| 真实试用 | `编号-bioSkills真实试用-<skill>.md`（完整深挖版） | `004-bioSkills真实试用-pairwise-alignment.md` |
| 引流精简 | `编号-bioSkills小红书-<skill>.md`（封面卡 + 实测，3–6 卡） | `004-bioSkills小红书-pairwise-alignment.md` |

> 笔记内引用图片用相对路径，例如 `../素材/004-pairwise/xxx.png`，移动文件后记得同步改路径。

---

## 笔记索引

| 编号 | 主题 | 系列 | 状态 | 文件 | 封面/图 |
|------|------|------|------|------|------|
| 001 | 生信 × AI 工作流 | 生信×AI | 可发 | [content/笔记/001-生信AI工作流笔记.md](content/笔记/001-生信AI工作流笔记.md) | content/素材/001-生信AI/（早期 AI 插画，已弃用） |
| 002 | bioSkills 使用体验 | skills 体验与架构 | 可发 | [content/笔记/002-bioSkills使用体验.md](content/笔记/002-bioSkills使用体验.md) | content/素材/002-bioSkills/ |
| 003 | bioSkills 拆解 · alignment（图文版） | skills 专项拆解 | 可发(配4张真实数据图) | [content/笔记/003-bioSkills拆解01-alignment专题.md](content/笔记/003-bioSkills拆解01-alignment专题.md) | content/素材/003-alignment/ (4张) |
| 004 | bioSkills 真实试用 · pairwise-alignment（完整） | skills 专项拆解 | 已入站点 | [content/笔记/004-bioSkills真实试用-pairwise-alignment.md](content/笔记/004-bioSkills真实试用-pairwise-alignment.md) | content/素材/004-pairwise/ |
| 004 | bioSkills 小红书 · pairwise-alignment（引流） | skills 专项拆解 | 引流精简 | [content/笔记/004-bioSkills小红书-pairwise-alignment.md](content/笔记/004-bioSkills小红书-pairwise-alignment.md) | content/素材/004-pairwise/ |
| 005 | bioSkills 真实试用 · msa-statistics（完整） | skills 专项拆解 | 已入站点 | [content/笔记/005-bioSkills真实试用-msa-statistics.md](content/笔记/005-bioSkills真实试用-msa-statistics.md) | content/素材/005-msa-statistics/ |
| 005 | bioSkills 小红书 · msa-statistics（引流） | skills 专项拆解 | 引流精简 | [content/笔记/005-bioSkills小红书-msa-statistics.md](content/笔记/005-bioSkills小红书-msa-statistics.md) | content/素材/005-msa-statistics/ |

**系列规划**
- `skills 使用体验`：从 bioSkills 仓库出发，逐个 skill 真跑、记录 AI 实际产出。**alignment 已发(003, 图文版)**；004/005 已深究成完整试用笔记并入站点。下期预告：variant-calling（体细胞变异）实测。
- `工具拆解`：fastp / samtools / cellsnp …… 一个一个拆（早期试写未采用，已清理）。
- `报错日记` / `项目复盘`：随缘补。

---

## 交付链路（受控闭环）

选题(S1) → 真跑复现 Gate(S2a) → 写笔记(S2b, `redbook-bio-note-writer`) → 转站点(S3, `bioSkills-lab-site`) → 小红书引流(S4, 原 DingCard **已停用**) → 互链闭环(S5)。

- **主交付物** = 笔记（md）+ 站点（`output/bioSkills-site/`）。
- **DingCard 已停用**：排版效果不好，不再产出小红书 ZIP。站点为当前唯一对外载体，本地预览 `python -m http.server 8123` 于 `output/bioSkills-site`。

---

## 如何新增一篇

1. 在 `content/笔记/` 新建 `00X-主题.md`（编号取当前最大 +1）。
2. 在 `content/素材/` 新建同名文件夹 `00X-主题/`，放该篇图片。
3. 笔记内用相对路径引用图片：`../素材/00X-主题/封面.png`。
4. 文案固定结构：**痛点切入 → 我的解法/体验 → 可抄命令或截图 → 避坑提醒**。不堆互动话术（关注引导 / 下期预告 / 评论区抛问题），只讲内容本身。
5. 发布前把正文里 `@XXX` 换成你的账号；封面在小红书编辑页加一行大字标题。
6. 回来在上方「笔记索引」表追加一行。
7. 重生成站点：`python pipeline/build_lab_site.py`（managed venv，已装 `markdown`）。

---

## 更新记录

- 2026-08-25（目录治理收尾）：① 源/产物分层固化——顶层收敛为 `content/`(源: 库/笔记/素材) · `pipeline/`(构建脚本+复现实验) · `output/`(纯产物 bioSkills-site)；② 删除早前 DingCard 全部遗留 ZIP/解压目录/预览图，环境纯净；③ 清理零依赖死物（早期 AI 插画、未发布草稿、规划存档、死代码）；④ 原 `实验/` 改名为 `pipeline/`（装的是正式构建脚本，非"实验"）；⑤ 原 `笔记/素材/库` 收拢进 `content/`，全文引用同步。端到端重生成验证通过（站点重建 + 三图正确拷入 `assets/`）。`bioSkills-trial-pipeline/SKILL.md` 的 S4 仍残留 DingCard 引用，按"先提不自动改"约定待拍板。
- 2026-08-25（DingCard 尾巴收口）：① `bioSkills-trial-pipeline/SKILL.md` 的 4 处 S4 引用改写为"已停用"，不再声称调 `dingcard-layout-export` 出 ZIP；② 删除 `content/库/DingCard/`（114MB 外部 git clone，可重建）；README 目录树移除该目录行。删后重跑站点构建成功，证明零硬依赖。注：`dingcard-layout-export` skill 自身（已标⚠停用）仍残留指向已删目录的悬空路径，仅供历史参考、不影响活链路。
- 2026-08-25：生成 bioSkills 真实数据封面（562 skills / 65 categories 分布图）；新增《bioSkills 合集简介》；修正 README 库统计数字（562/65）与文案结构（去掉互动话术要求）。
- 2026-08-24：初始化仓库，归档现有物料（2 篇正式笔记 + 1 篇草稿 + 7 张图），建立目录结构与本索引。
- 2026-08-23(晚)：003 从纯文字深度稿重排为**图文发布版**（配 4 张 matplotlib 真实数据图：MSA 着色 / 黄昏区曲线 / 保守度熵剖面 / 算法交叉图）；clone `库/bioSkills/` 作本地源码对照；README 补 `库/` 说明与系列进度。
