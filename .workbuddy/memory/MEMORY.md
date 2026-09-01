# RedBook 项目长期约定

## 1. bioSkills 真实试用纪律（006 起硬性）
- **workbuddy-skills 库**：remote=Diandian-ops/workbuddy-skills.git，本地 `D:/1.WorkDir/workbuddy-skills`（2026-09-01 新 clone），单 main 直接提交，skill 平铺目录。**git 身份坑**：全机所有仓库均无 user.name/email config，历史作者统一 `Diandian-ops <Diandian-ops@users.noreply.github.com>`，新 clone 库提交前须手工设 local identity。bio-note-fidelity-audit 已入库（165c6c4）。
- **完全复现**：严格按对应 skill 的 SKILL.md 规定的方法/命令/参数执行，不自行添加 skill 未提及的参数、步骤或内容。用户原话："skills要规定好了 完全复现 不要 自己加内容"。数据优先用 bioSkills 仓库自带示例，不自造玩具数据。
- **两份文件纪律**：每个主题必须同时产出
  - `content/笔记/<cat>/00X-bioSkills真实试用-<topic>.md` = **深拆版**（功能定位与适用范围 + 属性表 + 成分拆解 + 严格复现 + 实践要点 + 小结），进站点。
  - `content/笔记/<cat>/00X-bioSkills小红书-<topic>-正文.md` = **引流精简版**（封面卡 + 2-3 发现 + 结论），仅作发布稿、不进站点。**注意文件名带 `-正文` 后缀**（004–015 曾有 12 条链接因此全挂）。
- **客观描述、不拟人、不主观评判**：禁用"这个 skill 教对了什么""它教的是""生信老手踩过…"等表述；标题用「功能定位与适用范围」「实践要点」，边界用「适用范围：…不在本 skill 范围内」。
- **分类目录层**：`content/笔记/<cat>/` 与 `content/素材/<cat>/` 按 bioSkills 顶层分类名分层（alignment / read-alignment / variant-calling）；001–003 保持扁平不动。图片引用用 `../../素材/<cat>/00X-主题/图.png`。
- **素材目录须自带数据可独立复现**：`content/素材/<cat>/00X-/` 内含输入数据 + `make_figs.py` + `repro_transcript.txt`，脚本用 `os.path.dirname(os.path.abspath(__file__))` 自引用、图与脚本平铺。**禁止引用 `pipeline/`**（不入库，clone 后必跑不起来）。验收：从 `/tmp` 跑 `python <素材>/make_figs.py`。
- **站点侧栏需同步 `FAMILY_SUBS`**：`build_lab_site.py` 的 `TRIALS` 加行后，若该 cat 无 `FAMILY_SUBS` 条目则侧栏不高亮（015 曾踩）。
- **未执行必须标注**：未真跑/未安装的部分写进「未覆盖（诚实标注）」章节，或结尾段落等价披露「相关结论按 SKILL.md 原文陈述，未以构造数据演示」。**无素材目录 + 无标注 = 高风险**。
- **复现困难时的止损纪律（2026-09-01 用户拍板）**：工具链可下但输入数据结构性困难时，**不折腾装链**，转为「详细拆解 SKILL.md 内容 + 专设『复现难度评估』节如实说明障碍」。已执行于 022/023/027。关键环境事实：1000G 比对已全迁 **GRCh38DH 全基因组 CRAM**（每样本 8–15GB、无 chr 子集），旧 GRCh37 BAM 路径实测 404 → 需要 BAM 输入的 skill 在本机均无法真跑；DeepVariant 无 Windows 原生形态；本机无 WSL/Docker/samtools（2026-09-01 实测）。`pipeline/022-tools/downloads/` 留有 gatk zip(667MB)+jdk17(190MB)，gitignored，去留待用户定。

## 2. 出图与机检 Gate
- **出图自检**：当前模型读不了图，出图后一律跑 `python ~/.workbuddy/skills/bioSkills-figure-quality/check_figs.py make_figs.py`，要求 `TOTAL FAILS = 0`。
- `place_labels()` 必须在 `set_xlim/set_ylim` + `tight_layout()` **之后**调用，否则报假 PASS。`ax.get_yaxis_transform()` 中 x 是 axes 分数、y 是 data 值，混用会甩出轴外。
- **高频 FAIL：轴外底部注释**（017–020 连续 4 次）：`ax.text(0.5, -0.22, ...)` 必判越界 → 改为并入双行 `set_title()` 或放图内空白区（y≈0.05）。顶部同理须留余量。
- **gate_lint**：`pipeline/gate_lint.py` 已跳过围栏代码块、frontmatter、表格行。基线 004–027 全 `ERROR=0 WARN=0`。表头不算上下文——单位要写进单元格（如 `1409 bp`）。
- **`.gitignore` 只忽略 015 及之后工作区**：`pipeline/01[5-9]-*/` 与 `pipeline/02[0-9]-*/`。**不可 blanket 写成 `pipeline/0*/`**（004–006 早期工作区已入库）。

## 3. 小红书发布台 xhs-site
- 入口：双击 `output/xhs-site/index.html`（file://）或 `launch_console.bat`（8899；**8899 常被别的服务占用**，可改起 `serve_console.py 8900` 避开）。
- `scripts/` = 工具集：`build_xhs_console.py`（生成器=唯一真相源，改 UI 须改它再重跑）、`serve_console.py`、`md2card_patch.py`（Playwright 驱动出图）。
- 卡片落点 `output/xhs-site/xhs-cards/<cat>/<slug>/<theme>-standalone-hd/*.png`，四主题 = notebook/minimal/research/editorial。`output/` 整段 gitignored。
- **改发布台 UI 后必须**：改 `build_xhs_console.py` → 重跑生成器 → 重启服务进程（Windows 锁 .py，先 `taskkill` 旧 PID）。

## 4. md2card 出图运维（首选路径：前台 `python -c` 编排器）
- 起服务：`cd content/库/md2card && npm run dev`（5173，VITE ready）。出图：`<venv>/python.exe output/xhs-site/scripts/md2card_patch.py <真实试用md绝对路径> --theme <t> --cover standalone --out output/xhs-site/xhs-cards --url http://localhost:5173/app/ --hd`。
- **批量首选**：前台 `python -c` 内联编排器，`subprocess.run([VENV, PATCH, md, '--theme', th, ...])` 逐主题串行，用 `os.listdir` 判「有 PNG 且无 .nometa.md」跳过已完成。不受 120s 后台上限约束，37 主题一次跑完（`gen=37 ok=37 fail=0`）。
  13. **脚本不清理旧图**（2026-09-01 踩）：改稿重出时旧 PNG 残留，`ls | wc -l` 虚高（看着仍 22 张，实际新稿 11 张）。→ **改稿重出前必须 `rm -rf <slug>/*-standalone-hd/`**；判断真实张数用清空后的结果，或用 `os.path.getmtime` 按时间戳区分新旧。
  14. **改稿出图（降张数）**：摘录稿放 `output/xhs-site/card-ms/<cat>/<原slug>.md`（父目录名即 cat，文件名=原 slug 以便覆盖原卡片），出图加 `--assets <素材目录>` 显式绑图（摘录稿不在 `content/笔记/<cat>/` 下，自动推导会失败）。**正式笔记不动**。

- **坑速查**：
  1. **Git Bash 路径转换**：`/d/...` 被 MSYS 转成 `D:\d\...` → 用 PowerShell 或反斜杠 Windows 路径。
  2. **/api/regen 中文路径乱码** → 直接 venv python argv 传中文路径，不走接口。
  3. **后台任务 ~120s 硬上限**：必被杀并标 failed，长批量须分窗幂等续跑。
  4. **前台/后台 FS 隔离**（现象，机制未明）：前台 Write/heredoc 写的文件后台看不到；后台 PowerShell `Set-Content` 写真实 FS。编排脚本用前台 `python -c` 内联规避。
  5. **安全扫描拦截**（现象）：后台 `Start-Job`/`Start-Process` 被静默拦截；嵌 `%VAR%` 报 cmd 语法错。
  6. **md2card 单实例**：5173 共享，必须顺序串行，不可并行。
  7. **orphan-temp 重做**：cap 杀进程留 `.nometa.md` + `.zip` 无 PNG → 每窗前 `find -delete` 清理，且完成判定须含「无 .nometa.md」。
  8. **pre-check hang**：导出前 `wait_for_function(timeout=20000)` 恒 hang，已降到 `3000`（output/ 脚本，gitignored 可改），单主题 25s→8s。
  9. **枚举必须只取「真实试用」md**：`Where-Object { $_.Name -like '*真实试用*' }`，否则给小红书版出图并污染目录（用 `find <xhs-cards> -type d -name '*小红书*' -exec rm -rf {} +` 清）。
  10. **PowerShell 5.1 `Join-Path` 传 4 个位置参数抛异常**（非"静默丢弃"）：`ParameterBindingException` 被 `SilentlyContinue` 吞 → 返回 `$null` → 校验恒 false → 已完成主题被反复重做。用嵌套 `Join-Path`；**批量脚本不要用 SilentlyContinue 掩盖参数绑定错误**。
  11. **长编排脚本"1s 空转 no-op"**（根因未查明）：打印了枚举行却零 GEN，前后台 PowerShell 均偶发 → 改用前台 `python -c`。
  12. **`--density` 对张数无效**（2026-09-01 实测）：016 四主题 `balanced`→`compact` 张数不变（22→22）。md2card 分页由**内容长度/代码块数量**决定，降张数只能精简源稿。

## 5. 小红书张数硬约束（2026-09-01 用户补充）
- **平台硬上限 18 张**；用户心理预期 **10–16 张**。
- 当前 24 篇×4 主题共 96 组，张数 8–22（均值 12.4）。**016 = 22 张（超上限，发不出去）**、021 = 17 张（超预期）、024/026/027 = 16 张（边界）、012/013/018 = 8–9 张（略低但可发）。
- 超标处理只能精简源稿（或做卡片专用摘录稿，不动正式笔记），调 density 无效。

## 6. ⚠️ 环境现状：生信工具链已丢失（2026-09-01 探测）
- `bcftools` / `samtools` / `bgzip` / `tabix` **全盘搜索不到**；无 conda/mamba/miniforge 目录；WSL 未安装（无分发版）；PATH 无生信目录。
- 但素材 transcript 显示 **2026-08-28 尚在 `bcftools 1.24` 环境跑通** → 工具链在 8/28 之后被清理。
- **影响面**：① 022/023/027 补跑无基础链；② 已跑过的 9 篇素材 `make_figs.py`（依赖 bcftools）当前不可复现。
- 其他探测：java 仅 IDE 内置 jre21（`C:\Users\Admin\.antigravity\extensions\redhat.java-1.52.0-win32-x64\jre\21.0.9`）；docker/singularity/apptainer 均无；网络 **github 返回 000（不通）**、broad 返回 200（通）。→ DeepVariant（仅 Docker 分发）基本不可行；GATK 需 jar（GitHub 不通）+ Java + 参考基因组 + BAM。
- VCF 解析统一用 `bcftools query` 而非 cyvcf2（cyvcf2 在受管 venv 重建后已丢失）。

## 7. 表述纪律：实测 vs 推断（2026-09-01 用户质疑后确立）
- 凡汇报结论必须能指出**证据来源**（命令 stdout、文件内容、`os.listdir` 统计、curl 状态码 = 硬证据）；推断须显式标注"我推测/可能"。
- **禁止把推断写成机制**。无实测支撑时只能写**现象**并注明「根因未查明」。已发生一次自我修正：坑10 原把"异常被吞返回 null"写成"静默丢弃参数"。坑4 "overlay"、坑5 "安全扫描"同为推断标签，按现象对待。

## 8. 笔记保真度核查（可复用流程）
- 已固化为用户级 skill **`bio-note-fidelity-audit`**（`~/.workbuddy/skills/`）：`python scripts/audit_fidelity.py --root <repo> --cat <cat> --from 016 --to 027`。
- 五步：命令/参数溯源 → 人工复核标记项 → 实测证据链 → 未执行诚实标注 → 表述客观性与偏差标注。
- **机器初筛只做筛选不做判定**：016–027 核查中标记 5 项，人工核实后 **4 项误报**（环境版本行被当命令、`gatk VariantFiltration` 因换行未连续匹配）。
- 另一坑：Bash `cat` 读中文路径文件显示为空 → 须用 python `io.open` 复核（020 transcript 因此一度误判无内容）。
- 2026-09-01 已核查 016–027：**12 篇均未发现无出处编造**。报告 `output/audit/笔记可信度核查-016-027.md`（不入库），脚本 `pipeline/audit/audit_note_vs_skill.py`（入库）。
