#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""RedBook / bioSkills 真实试用实验室 — 跨平台环境自检。

用途：换机（macOS <-> Windows）后，先跑这个脚本确认「这台机器能不能完整复现」，
再去做 clone / 装依赖 / 出图 / 建站。把「忘了 clone 外部仓库」这类隐性阻塞
在动手前就暴露出来，而不是等到 build_lab_site.py 抛 FileNotFoundError。

设计约束（针对 Windows 原生环境）：
  1. 纯标准库 —— 必须在 `pip install -r requirements.txt` **之前**就能跑。
  2. 不使用 shell=True —— Windows 是 cmd.exe，管道/引号语义与 POSIX 不同。
  3. 中文输出做编码兜底 —— Windows 控制台默认非 UTF-8，直接 print 中文会抛异常。
  4. 全部路径用 os.path / pathlib —— 不假设 "/" 分隔符。

用法：
    python pipeline/check_env.py            # 离线自检（默认，含小红书链路 WARN）
    python pipeline/check_env.py --network  # 额外探测 htslib 远程数据读取能力
    python pipeline/check_env.py --xhs      # 严格门禁：小红书出图链路缺失即 FAIL

退出码：
    0  无 FAIL（可复现；WARN 不影响）
    1  存在 FAIL（当前状态下某些环节跑不通）
"""

import argparse
import importlib
import os
import subprocess
import sys

# --- Windows 控制台编码兜底：中文路径/提示在 cp936 下会 UnicodeEncodeError ---
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

PASS, WARN, FAIL = "PASS", "WARN", "FAIL"
_ORDER = {PASS: 0, WARN: 1, FAIL: 2}


class Result:
    def __init__(self, status, name, detail, hint=""):
        self.status = status
        self.name = name
        self.detail = detail
        self.hint = hint


def _run(cmd, timeout=30, cwd=None):
    """跨平台执行外部命令。列表形式，不经 shell。

    注意 cwd：git 类命令必须显式指定仓库根目录，
    否则读到的是「当前工作目录」的配置而非本仓库的（曾因此漏检 core.autocrlf）。
    """
    try:
        p = subprocess.run(cmd, capture_output=True, timeout=timeout, cwd=cwd)
    except FileNotFoundError:
        return None, "", "命令不存在"
    except subprocess.TimeoutExpired:
        return None, "", "执行超时"
    except OSError as e:
        return None, "", "执行失败: %s" % e
    out = p.stdout.decode("utf-8", "replace") if p.stdout else ""
    err = p.stderr.decode("utf-8", "replace") if p.stderr else ""
    return p.returncode, out, err


def _first_line(text):
    for line in text.splitlines():
        if line.strip():
            return line.strip()
    return ""


# ---------------------------------------------------------------- 检查项

def check_python():
    v = sys.version_info
    ver = "%d.%d.%d" % (v.major, v.minor, v.micro)
    if (v.major, v.minor) < (3, 9):
        return Result(FAIL, "Python 版本", ver, "需要 Python >= 3.9")
    return Result(PASS, "Python 版本", ver, "")


def check_packages():
    out = []
    status = PASS
    for mod, pipname in (("matplotlib", "matplotlib"), ("markdown", "markdown")):
        try:
            m = __import__(mod)
            ver = getattr(m, "__version__", "unknown")
            out.append("%s %s" % (mod, ver))
        except ImportError:
            status = FAIL
            out.append("%s 未安装" % mod)
    hint = "pip install -r requirements.txt" if status == FAIL else ""
    return Result(status, "Python 依赖", " / ".join(out), hint)


def _tool_version(tool, min_recommended=None):
    """查找外部二进制并取版本。返回 (Result)。"""
    from shutil import which
    path = which(tool)
    if not path:
        return Result(FAIL, "外部工具 %s" % tool, "未找到",
                      "conda install -c conda-forge %s" % tool)
    rc, out, err = _run([tool, "--version"])
    ver = _first_line(out) or _first_line(err) or "已找到但无法取版本"
    hint = ""
    if rc is None:
        return Result(WARN, "外部工具 %s" % tool, "%s（%s）" % (path, ver), "")
    return Result(PASS, "外部工具 %s" % tool, "%s  @ %s" % (ver, path), hint)


def check_tools():
    res = []
    res.append(_tool_version("bcftools"))
    # 以下是可选工具：只有部分素材脚本用到
    for optional in ("samtools", "bgzip", "tabix"):
        r = _tool_version(optional)
        if r.status == FAIL:
            r.status = WARN
            r.detail = "未找到（仅部分脚本需要）"
        res.append(r)
    return res


def check_external_clones():
    """content/库/ 下的外部仓库：被 .gitignore 忽略，换机必须重新 clone。"""
    target = os.path.join(ROOT, "content", "库", "bioSkills")
    res = []
    if os.path.isdir(target) and os.listdir(target):
        n_cat = len([d for d in os.listdir(target)
                     if os.path.isdir(os.path.join(target, d)) and not d.startswith(".")])
        res.append(Result(PASS, "外部 clone bioSkills", "已就绪，%d 个分类" % n_cat, ""))
    else:
        res.append(Result(
            FAIL, "外部 clone bioSkills", "缺失",
            "build_lab_site.py 会在此崩溃。执行：\n"
            "        git clone https://github.com/GPTomics/bioSkills.git content/库/bioSkills"))
    return res


def _pw_browsers_path():
    """Playwright 浏览器缓存目录（按平台；可被环境变量覆盖）。"""
    env = os.environ.get("PLAYWRIGHT_BROWSERS_PATH")
    if env:
        return env
    if os.name == "nt":
        return os.path.join(os.path.expandvars("%LOCALAPPDATA%"), "ms-playwright")
    if sys.platform == "darwin":
        return os.path.expanduser("~/Library/Caches/ms-playwright")
    return os.path.expanduser("~/.cache/ms-playwright")


def _has_chromium():
    base = _pw_browsers_path()
    if not os.path.isdir(base):
        return False
    return any(d.startswith("chromium") for d in os.listdir(base))


def check_xhs_deps(strict=False):
    """小红书出图链路（md→图文）依赖：md2card dev server + playwright + chromium。

    默认非阻塞（缺失判 WARN，不影响建站）；--xhs 时升级为 FAIL，
    作为「能否跑 trial2xhs.py」的严格门禁。
    """
    res = []
    # 外部 clone
    for name, url in (
        ("md2card", "https://github.com/haodongcui/md2card.git"),
        ("obsidian-wechat-converter",
         "https://github.com/DavidLam-oss/obsidian-wechat-converter.git"),
    ):
        p = os.path.join(ROOT, "content", "库", name)
        if os.path.isdir(p) and os.listdir(p):
            res.append(Result(PASS, "xhs 外部 clone %s" % name, "已就绪", ""))
        else:
            res.append(Result(
                WARN, "xhs 外部 clone %s" % name, "缺失",
                "git clone %s content/库/%s" % (url, name)))

    # node / npm（md2card dev server 需要）
    for tool in ("node", "npm"):
        r = _tool_version(tool)
        if r.status == FAIL:
            r.status = WARN
            r.detail = "未找到（md2card dev server 需要）"
            r.hint = ("conda install -c conda-forge nodejs   或   "
                      "winget install OpenJS.NodeJS")
        res.append(r)

    # playwright 模块
    try:
        importlib.import_module("playwright")
        res.append(Result(PASS, "xhs playwright", "已安装", ""))
    except ImportError:
        res.append(Result(
            WARN, "xhs playwright", "未安装",
            "pip install playwright pillow && playwright install chromium"))

    # chromium 浏览器
    if _has_chromium():
        res.append(Result(PASS, "xhs chromium", "已安装", ""))
    else:
        res.append(Result(WARN, "xhs chromium", "未安装", "playwright install chromium"))

    if strict:
        for r in res:
            if r.status == WARN:
                r.status = FAIL
    return res


def check_gitattributes():
    p = os.path.join(ROOT, ".gitattributes")
    if not os.path.isfile(p):
        return Result(WARN, ".gitattributes", "缺失",
                      "二进制文件可能被 core.autocrlf 损坏，应随仓库提供")
    return Result(PASS, ".gitattributes", "已就位", "")


def check_autocrlf():
    rc, out, err = _run(["git", "config", "--get", "core.autocrlf"], cwd=ROOT)
    val = _first_line(out).lower() if rc == 0 else ""
    if rc != 0 or val in ("", "false", "input"):
        return Result(PASS, "git core.autocrlf", val or "未设置(等同 false)", "")
    if val == "true":
        return Result(WARN, "git core.autocrlf", "true",
                      "Windows 上建议 git config --global core.autocrlf false，"
                      "配合 .gitattributes 使用")
    return Result(WARN, "git core.autocrlf", val, "非常规取值，请确认")


def check_material_scripts():
    """素材脚本自包含性：每个脚本与其数据应同目录，且能用 HERE 定位。"""
    mat_root = os.path.join(ROOT, "content", "素材")
    if not os.path.isdir(mat_root):
        return Result(FAIL, "素材目录", "缺失", "")
    found, bad = 0, []
    for dirpath, dirnames, filenames in os.walk(mat_root):
        dirnames[:] = [d for d in dirnames if not d.startswith(".")]
        for fn in filenames:
            if fn in ("make_figs.py", "make_fig.py"):
                found += 1
                p = os.path.join(dirpath, fn)
                try:
                    with open(p, "r", encoding="utf-8") as f:
                        src = f.read()
                except OSError:
                    bad.append(os.path.relpath(p, ROOT))
                    continue
                # 素材脚本不得回溯引用 pipeline/ 工作区（该目录不入库，clone 后必然跑不起来）。
                # 真实写法是 os.path.join(HERE, "..", "..", ..., "pipeline", ...) 的分散参数，
                # 源码里并不存在字面量 "../../"，所以按「同一行同时出现 .. 与 pipeline」判定。
                for line in src.splitlines():
                    if "pipeline" in line and ('".."' in line or "'..'" in line):
                        bad.append(os.path.relpath(p, ROOT))
                        break
    if bad:
        return Result(FAIL, "素材脚本自包含", "%d 个脚本引用了 pipeline/" % len(bad),
                      "；".join(bad[:3]))
    return Result(PASS, "素材脚本自包含", "已扫描 %d 个出图脚本，无 pipeline/ 硬引用" % found, "")


def check_network_htslib():
    """探测 htslib 是否支持远程 URL（libcurl）。只在 --network 时执行。

    注意：仓库内的素材脚本零网络依赖（数据已随素材目录入库），
    这里只是为「需要重新拉取公开数据时」提供能力确认。
    """
    url = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
           "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")
    rc, out, err = _run(["bcftools", "view", "-h", url], timeout=90)
    if rc == 0 and out.startswith("##fileformat"):
        return Result(PASS, "htslib 远程读取", "支持（可直接读 https VCF）", "")
    msg = (err or out or "").strip().splitlines()
    return Result(WARN, "htslib 远程读取", "不可用" if rc else "未知",
                  "仅影响重新下载公开数据；仓库内素材复现不需要。" + (msg[0][:80] if msg else ""))


# ---------------------------------------------------------------- 主流程

def main():
    ap = argparse.ArgumentParser(description="RedBook 跨平台环境自检")
    ap.add_argument("--network", action="store_true",
                    help="额外探测 htslib 远程数据读取能力（需要联网，较慢）")
    ap.add_argument("--xhs", action="store_true",
                    help="严格门禁：把小红书出图链路(node/npm/playwright/chromium/md2card)"
                         "的依赖从 WARN 升级为 FAIL，用于确认能跑 trial2xhs.py")
    args = ap.parse_args()

    results = []
    results.append(check_python())
    results.append(check_packages())
    results.extend(check_tools())
    results.extend(check_external_clones())
    results.append(check_gitattributes())
    results.append(check_autocrlf())
    results.append(check_material_scripts())
    results.extend(check_xhs_deps(strict=args.xhs))
    if args.network:
        results.append(check_network_htslib())

    name_w = max(len(r.name) for r in results) + 2
    print("=" * 78)
    print("RedBook 环境自检   (root: %s)" % ROOT)
    print("平台: %s" % sys.platform)
    print("=" * 78)
    for r in results:
        print("[%-4s] %-*s %s" % (r.status, name_w, r.name, r.detail))
        if r.hint:
            for line in r.hint.split("\n"):
                print("         -> %s" % line)
    print("-" * 78)

    n_fail = sum(1 for r in results if r.status == FAIL)
    n_warn = sum(1 for r in results if r.status == WARN)
    n_pass = sum(1 for r in results if r.status == PASS)
    print("PASS=%d  WARN=%d  FAIL=%d" % (n_pass, n_warn, n_fail))

    if n_fail:
        print("结论：当前环境存在阻塞项，按上面的 -> 提示处理后重跑。")
        return 1
    if n_warn:
        print("结论：可复现（有 %d 项非阻塞告警，不影响出图/机检/建站）。" % n_warn)
    else:
        print("结论：环境就绪，可完整复现。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
