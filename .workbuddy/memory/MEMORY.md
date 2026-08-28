# RedBook 项目长期约定

## bioSkills 真实试用纪律（006 起硬性，2026-08-25 用户明确）
- 做 bioSkills 真实试用（004/005/006...）时，严格按对应 skill 的 SKILL.md 规定的方法/命令/参数执行，**完全复现**，不自行添加 skill 未提及的参数、步骤或内容。
- skill 文档是唯一事实源：真跑只验证 skill 写的方法能否跑通、出什么图，不发挥、不夹带私货。
- 用户原话："skills要规定好了 完全复现 不要 自己加内容"。
- 数据输入优先用 bioSkills 仓库自带示例（如 `content/库/bioSkills/alignment/alignment-io/examples/sample_alignment.aln`），不自造玩具数据。
- **两份文件纪律（006 后硬性，用户原话"真实试用也要像005一样 skill成分拆解，只有小红书才精简"）**：每个主题必须同时产出
  - `content/笔记/00X-bioSkills真实试用-<topic>.md` = **深拆版**（对齐 005 风格：功能定位与适用范围+属性表+成分拆解[文件结构/各脚本/工具知识/API/经验坑]+严格复现[环境/数据/标准输出/坑实测]+实践要点+小结），进 `output/bioSkills-site/` 站点。
  - `content/笔记/00X-bioSkills小红书-<topic>.md` = **引流精简版**（封面卡+2-3发现+结论，不展开成分拆解），仅作发布稿、不进站点。
  - 不可只产其一。规范已写入 `redbook-bio-note-writer/SKILL.md` 第八节（8.1 深拆模板 / 8.2 精简模板 + 两文件纪律表）。
- **写作风格：客观描述、不拟人、不主观评判（2026-08-25 用户明确，原话"语句带情感带不和谐，要用描述和真实实践视角"）**：笔记用描述性语言陈述 skill 的功能/输入/输出/实践结论，避免把 skill 当主体做主观评判。
  - 禁用："这个 skill 教对了什么""这个 skill 不负责做 X""它教的是""生信老手踩过…后来人不用重复踩"等拟人/主观表述。
  - 中性替代：标题用「功能定位与适用范围」「实践要点」；边界用「适用范围：本 skill 的输入为…，比对构建由 multiple-alignment 覆盖，不在本 skill 范围内」；内容用「内容覆盖：…」。
  - 规范已固化进 `redbook-bio-note-writer/SKILL.md` 第四节（客观描述纪律）+ 8.1/8.2 模板标题。004/005/006 三篇笔记已全量改写为中性写法。
- **分类目录层纪律（2026-08-28 用户确认「复用 bioSkills 分类名 + 现在就重排」）**：004–015 已从扁平移入按 bioSkills 顶层分类名分层的子目录，从根消除后期堆砌。
  - `content/笔记/<cat>/` = 两文件（真实试用 + 小红书）按 `alignment` / `read-alignment` / `variant-calling` 分层；001–003 早期杂项保持扁平不动。
  - `content/素材/<cat>/` = 同构分层，与笔记一一对应。
  - `output/xhs-cards/<cat>/` = 同构分层（output 不进库）。
  - 图片引用：笔记在 `<cat>/` 下时正文用 `../../素材/<cat>/00X-主题/图.png`（退两层到仓库根再进 素材/<cat>）；文字提及素材用 `content/素材/<cat>/00X-主题/`。
  - 站点构建器 `pipeline/build_lab_site.py` 已**配置化**：`TRIALS` 列表 + `cat_page()` 自动扫描 `content/笔记/<cat>/`，**新增一篇只需往 TRIALS 加一行**，不再手写三处硬编码映射。
  - 发布脚本 `md2card_automation.py` / `trial2xhs.py` 从输入 md 路径推导 cat 层，产物落到 `output/xhs-cards/<cat>/<slug>/`，素材根退两层到 `content/素材/<cat>/`。
  - 新增分类时：建 `content/笔记/<新cat>/` + `content/素材/<新cat>/`，在 `TRIALS` 加条目即可；README 索引表同步加行。
  - 提交 `5658729`（716 files）已完成该重排 + 配置化。
- **小红书文件名带 `-正文` 后缀（易漏）**：实际文件是 `00X-bioSkills小红书-<topic>-正文.md`，README/索引引用时必须带后缀（004–015 曾有 12 条链接因此全挂，2026-08-28 修复）。
- **素材目录须自带数据可独立复现（016 起）**：`content/素材/<cat>/00X-/` 内除图外，还要放输入数据 + `make_figs.py` + `repro_transcript.txt`，脚本用 `os.path.dirname(os.path.abspath(__file__))` 取路径、图与脚本**平铺**（对齐 008/011）。`pipeline/00X-*/` 工作区保留中间产物、**不入库**（015/016 皆如此，用户「不希望延伸」）。
- **站点侧栏需同步 `FAMILY_SUBS`（易漏）**：`build_lab_site.py` 里 `TRIALS` 加行后，若该 cat 还没有 `FAMILY_SUBS` 条目，侧栏不会高亮也不会出现子页链接（015 曾因此从未被高亮）。新增 cat 的第二篇时补一条 `VARCALL_SUBS` 式列表即可。

## 出图质量：程序化自检 Gate（016 起）
- 当前模型读不了图，**出图后一律跑** `python ~/.workbuddy/skills/bioSkills-figure-quality/check_figs.py make_figs.py`，要求 `TOTAL FAILS = 0`（退出码可接 CI）。它 hook `Figure.savefig`，落盘前量所有 Text 的 window bbox，判定「文本越界 / 文本互相重叠 / legend 压数据点」。
- `place_labels()`（fig_quality.py）必须在 `set_xlim/set_ylim` + `tight_layout()` **之后**调用：顺序反了避让会算在自动缩放坐标系下，渲染时坐标一变标签整体偏移，且自检与渲染同处错误坐标系会报**假 PASS**。
- `ax.get_yaxis_transform()` 中 **x 是 axes 分数、y 是 data 值**，混用会把文字甩到轴外被裁。
- 数据点扎堆时（016 有 6/9 位点 GQ 顶到 99 上限、两个 QUAL 仅差 0.16）固定 `xytext` 必然重叠，只能贪心避让；柱状图数值标签与顶部参考线说明冲突时，把数值标签放**柱底内侧**（白字加粗）。

## 笔记机检 gate_lint（2026-08-28 改版）
- `pipeline/gate_lint.py` 已跳过 ``` 围栏代码块、YAML frontmatter、markdown 表格行——此前命令行/报错原文被当叙述句判 ADJ_REPEAT、frontmatter 被判 BARE_NUM，噪声淹没真信号。
- 现状基线：004–016 全部 `ERROR=0 WARN=0`；002 的 METAPHOR 是早期遗留（001–003 保持扁平不动，不改）。改完工具或新增笔记后跑一次全量回归确认无新增 ERROR。
