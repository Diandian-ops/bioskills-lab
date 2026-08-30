# 换机复现任务（交给新 agent 执行）

你要在**一台新的 Windows 11 电脑**上，把 `bioskills-lab` 这个项目仓库完整、可复现地搭起来并验证能正常复现。

仓库本身**已经做过跨平台改造**，你不需要修复任何代码，只做「拉取 + 装依赖 + 验证」：
- 路径不再硬编码 macOS（`/Applications/anaconda3/...` 已改为三级查找）
- 入库的 152 个二进制文件已有 `.gitattributes` 声明 `binary` 兜底
- 有 `pipeline/check_env.py` 一键自检（纯标准库，装依赖前就能跑）

## 前置条件（先确认）
- Windows 11，已装 `git`（建议 2.40+）
- 已装 `conda`（Miniconda / Anaconda 任一；没有先去装 Miniconda）
- 能联网（要 clone 两个 GitHub 仓库 + conda 装包）
- 用 **PowerShell** 执行下方命令（不是 cmd）

## 一次性执行（复制整段到 PowerShell 粘贴，中途不用停）

```powershell
# 0. 关 autocrlf —— 必须在 clone 之前执行，否则二进制（bam/gz/png）会被改写损坏
git config --global core.autocrlf false

# 1. 拉主仓（已含全部笔记/素材/出图，约 18MB；不要打包拷贝，直接 clone）
git clone https://github.com/Diandian-ops/bioskills-lab.git RedBook
cd RedBook

# 2. 拉外部参考仓（建站必需，已被 gitignore 排除，必须单独 clone）
git clone https://github.com/GPTomics/bioSkills.git content/库/bioSkills

# 3. 装依赖（conda 一条命令搞定；无 conda 见下方故障兜底）
conda install -c conda-forge python=3.12 bcftools samtools matplotlib markdown -y

# 4. 验证：FAIL=0 即环境就绪
python pipeline/check_env.py

# 5. 重生成站点（产物在 output/bioSkills-site/）
python pipeline/build_lab_site.py
```

## 成功标准（必须全部满足才算完成）
1. 第 4 步 `python pipeline/check_env.py` 输出 **`FAIL=0`**（`WARN` 可忽略，不影响复现）。
2. 第 5 步 `build_lab_site.py` 无报错退出，生成 `output/bioSkills-site/`。
3. 站点页数与图数符合预期（参考基线：28 页 / 34 图）。

## 可选深度验证（建议做一次，确认素材可独立复现）
```powershell
cd /tmp
python <RedBook仓库路径>/content/素材/variant-calling/018-vcf-statistics/make_figs.py
```
应无报错并生成 PNG。所有素材脚本均零网络依赖、从任意 cwd 可跑。

## 故障兜底
- **`010–014` 六个脚本在真机 Windows 报管道相关错**：不要改代码，直接改用 **WSL2**（Linux 环境，与开发机 mac 行为一致，最稳）。在 WSL2 里重新走上面流程即可。
- **conda 装包慢/失败**：换 `winget install bcftools samtools` + `pip install matplotlib markdown`（Python 用系统或微软商店版）。
- **`content/库/bioSkills` clone 失败/太慢**：它只是建站必需，站点浏览和笔记阅读不依赖它。可先跳过第 2 步；但 `build_lab_site.py` 会提示缺它，补回即可。
- **`core.autocrlf` 忘了提前关导致二进制损坏**：删除整个仓库目录重新 clone（这次先关 autocrlf）。

## 不要做的事
- 不要用 U 盘 / 压缩包拷贝旧机器文件来「转移」，一律用 `git clone`。
- 不要修改仓库目录结构、脚本或笔记内容——你只负责把环境搭起来并验证。
- 不要 `git push` 或改动远程——这是只读搭建任务。

## 完成后汇报
- `check_env.py` 的完整输出（尤其 FAIL / WARN 计数）
- `build_lab_site.py` 是否成功、生成了几页几图
- 是在原生 PowerShell 还是 WSL2 下完成
- 任何与上面「成功标准」不符的偏差
