#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
trial2xhs.py — bioSkills 笔记 → 小红书出图 的编排器（human-in-the-loop 的自动化段）。

把「可自动化」的两步串成一条命令：
  1) Gate 机检：先跑 gate_lint 拦下机械违规（ERROR 级不让进出图）
  2) md2card 出图：lint 通过后用 Playwright 驱动本地 md2card 工作台，
     按标准化默认（notebook 主题 / standalone 独立封面 / --hd 2160x2880）出图，
     产物落 output/xhs-site/xhs-cards/<分类>/<slug>/notebook-standalone-hd/
     （分类从笔记路径 content/笔记/<分类>/<md> 自动推导；001-003 无分类则落到 <out>/<slug>/）

设计原则（对齐用户判定）：
  - 只自动化「格式/质量门槛 + 出图」，不自动化「选题 / 写内容 / 发帖」，
    那三步必须人来做（S1/S3 内容需人写、S5 无 API 需手动发）。
  - lint 是质量闸门：ERROR 直接拦下，不让半成品进出图。

用法：
  python pipeline/trial2xhs.py <00X-bioSkills真实试用-<topic>.md> [--force] [--no-server-check]
                              [--theme notebook] [--cover standalone] [--density balanced] [--hd]

前置：
  - md2card dev server: cd content/库/md2card && npm i && npm run dev  (http://localhost:5173/app/)
  - venv(playwright): ~/.workbuddy/binaries/python/envs/default
"""
import argparse
import os
import subprocess
import sys
import urllib.request

VENV_PY = os.path.expanduser("~/.workbuddy/binaries/python/envs/default/bin/python")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GATE = os.path.join(ROOT, "pipeline", "gate_lint.py")
MD2CARD = os.path.join(ROOT, "pipeline", "md2card_automation.py")
URL = "http://localhost:5173/app/"
OUT_ROOT = os.path.join(ROOT, "output", "xhs-site", "xhs-cards")


def run_gate(md):
    print("─── [Gate 1/2] 机检笔记可读性 + 雷区 ───")
    r = subprocess.run([VENV_PY, GATE, md])
    # 退出码 1 = 有 ERROR
    return r.returncode


def server_ok(url):
    try:
        urllib.request.urlopen(url, timeout=3)
        return True
    except Exception:
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("md", help="真实试用笔记 md（如 00X-bioSkills真实试用-<topic>.md）")
    ap.add_argument("--force", action="store_true", help="Gate 有 ERROR 也强行继续（不推荐）")
    ap.add_argument("--no-server-check", action="store_true", help="跳过 md2card dev server 探活")
    ap.add_argument("--theme", default="notebook",
                    choices=["minimal", "research", "editorial", "notebook"])
    ap.add_argument("--cover", default="standalone",
                    choices=["integrated", "standalone", "none"])
    ap.add_argument("--density", default="balanced",
                    choices=["relaxed", "balanced", "compact"])
    ap.add_argument("--hd", action="store_true", default=True, help="高清原图(2160x2880)，默认开")
    ap.add_argument("--no-hd", dest="hd", action="store_false", help="关闭高清，用标准发布(1080x1440)")
    args = ap.parse_args()

    md = os.path.abspath(args.md)
    if not os.path.isfile(md):
        print(f"[!] 找不到 md: {md}", file=sys.stderr)
        sys.exit(2)

    # ---- Gate 1: lint ----
    rc = run_gate(md)
    if rc == 1 and not args.force:
        print("\n[✗] Gate 未通过（存在 ERROR 级机械违规）。先改完再出图。")
        print("    可加 --force 强行继续，但不推荐。")
        sys.exit(2)
    if rc == 1 and args.force:
        print("\n[!] --force：忽略 ERROR 继续出图（产物可能含违规表述）")
    elif rc == 0:
        print("[✓] Gate 通过，无机械违规")

    # ---- Gate 2: md2card server 探活 ----
    if not args.no_server_check:
        print("\n─── [Gate 2/2] md2card dev server 探活 ───")
        if not server_ok(URL):
            print(f"[✗] 连不上 {URL}")
            print("    请先启动：cd content/库/md2card && npm i && npm run dev")
            print("    （或加 --no-server-check 跳过本检查）")
            sys.exit(3)
        print(f"[✓] {URL} 可达")

    # ---- 出图 ----
    print("\n─── 出图（md2card：%s / %s / %s）───" % (
        args.theme, args.cover, "hd" if args.hd else "std"))
    cmd = [VENV_PY, MD2CARD, md,
           "--theme", args.theme, "--cover", args.cover,
           "--density", args.density, "--out", OUT_ROOT]
    if args.hd:
        cmd.append("--hd")
    r = subprocess.run(cmd)
    if r.returncode != 0:
        print("[✗] 出图失败，见上方 md2card 日志")
        sys.exit(r.returncode)

    # 产物路径提示
    import re
    base = os.path.basename(md)
    base = re.sub(r"\.(md|markdown|txt)$", "", base, flags=re.I)
    base = re.sub(r"[^\w\u4e00-\u9fff-]+", "_", base).strip("_") or "note"
    variant = f"{args.theme}-{args.cover}" + ("-hd" if args.hd else "")
    cat = os.path.basename(os.path.dirname(os.path.abspath(md)))
    out_parts = [OUT_ROOT]
    if cat and cat != "笔记":
        out_parts.append(cat)
    out_parts += [base, variant]
    out_dir = os.path.join(*out_parts)
    print(f"\n[✓] 完成。产物目录：\n    {out_dir}")
    sys.exit(0)


if __name__ == "__main__":
    main()
