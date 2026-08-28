#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
md2card_automation.py — 用 Playwright 驱动本地 md2card 工作台，把 bioSkills 笔记 md
自动转成小红书 3:4 多卡 PNG。出图源稿建议用「真实试用」版 md（内容全、与站点一致）。

前提：
  - content/库/md2card 已 `npm install` 并 `npm run dev` 跑在 http://localhost:5173/app/
  - venv 含 playwright + Pillow：~/.workbuddy/binaries/python/envs/default

用法：
  python pipeline/md2card_automation.py <note.md> [--theme notebook] [--density balanced]
                                         [--cover standalone] [--canvas 3:4] [--hd]
                                         [--assets DIR] [--out output/xhs-cards] [--url URL]

产物（ZIP 下载后已清理，只留 PNG）：
  <out>/<分类>/<slug>/<theme>-<cover>[-hd]/01.png, 02.png, ...
  例：output/xhs-cards/alignment/005-bioSkills真实试用-msa-statistics/notebook-standalone-hd/01.png
  （分类从笔记路径 content/笔记/<分类>/<md> 自动推导；001-003 无分类则落到 <out>/<slug>/）
"""
import argparse
import os
import re
import sys
import zipfile
import shutil
import glob

from playwright.sync_api import sync_playwright

VENV = os.path.expanduser("~/.workbuddy/binaries/python/envs/default")
if os.path.exists(os.path.join(VENV, "bin", "python")):
    os.environ.setdefault("PLAYWRIGHT_BROWSERS_PATH",
                          os.path.expanduser("~/Library/Caches/ms-playwright"))

DENSITY_TEXT = {"relaxed": "舒展", "balanced": "技术平衡", "compact": "紧凑"}
COVER_TEXT = {"integrated": "融合首页", "standalone": "独立封面", "none": "无封面"}


def slugify(path: str) -> str:
    base = os.path.basename(path)
    base = re.sub(r"\.(md|markdown|txt)$", "", base, flags=re.I)
    base = re.sub(r"[^\w\u4e00-\u9fff-]+", "_", base)
    return base.strip("_") or "note"


def find_assets_dir(md_path: str, explicit: str | None) -> str | None:
    if explicit:
        return explicit if os.path.isdir(explicit) else None
    note_dir = os.path.dirname(os.path.abspath(md_path))
    slug = slugify(md_path)
    # 分类 = 笔记所在子目录名（content/笔记/<cat>/<md> -> cat）；
    # 若笔记直接放在 content/笔记/ 下（早期 001-003），则无分类层。
    cat = os.path.basename(note_dir)
    # 笔记在 content/笔记/<cat>/<md>，退两层到 content/，再进 素材
    assets_root = os.path.abspath(os.path.join(note_dir, "..", "..", "素材"))
    if cat and cat != "笔记" and os.path.isdir(os.path.join(assets_root, cat)):
        assets_root = os.path.join(assets_root, cat)
    if not os.path.isdir(assets_root):
        return None
    m = re.match(r"^(\d+)", slug)
    num = m.group(1) if m else ""
    candidates = [os.path.join(assets_root, slug)]
    if num:
        last_seg = slug.rsplit("-", 1)[-1]
        candidates.append(os.path.join(assets_root, f"{num}-{last_seg}"))
        for d in sorted(os.listdir(assets_root)):
            if d.startswith(num + "-") and os.path.isdir(os.path.join(assets_root, d)):
                candidates.append(os.path.join(assets_root, d))
    candidates.append(assets_root)
    for c in candidates:
        c = os.path.abspath(c)
        if os.path.isdir(c) and glob.glob(os.path.join(c, "*.png")):
            return c
    return None


def run(md_path, theme, density, cover, canvas, hd, assets, out, url):
    md_path = os.path.abspath(md_path)
    if not os.path.isfile(md_path):
        print(f"[!] 找不到 md: {md_path}", file=sys.stderr)
        return 1

    assets_dir = find_assets_dir(md_path, assets)
    print(f"[*] 笔记: {md_path}")
    print(f"[*] 素材目录: {assets_dir or '(无，跳过绑图)'}")

    slug = slugify(md_path)
    variant = f"{theme}-{cover}" + ("-hd" if hd else "")
    # 分类层：从笔记所在子目录推导，使产物落到 output/xhs-cards/<cat>/<slug>/<variant>
    cat = os.path.basename(os.path.dirname(md_path))
    out_parts = [out]
    if cat and cat != "笔记":
        out_parts.append(cat)
    out_parts += [slug, variant]
    out_dir = os.path.abspath(os.path.join(*out_parts))
    os.makedirs(out_dir, exist_ok=True)

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        # 删除 showDirectoryPicker，强制走隐藏 input 回退（避免原生目录选择器拦截）
        context = browser.new_context(
            accept_downloads=True,
            viewport={"width": 1440, "height": 1000},
        )
        context.add_init_script("delete window.showDirectoryPicker;")
        page = context.new_page()
        page.goto(url, wait_until="networkidle")
        page.wait_for_selector("text=Markdown 编辑器", timeout=30000)
        print("[*] 工作台已加载")

        # 1) 导入 Markdown（触发 loadFile -> 自动弹图片绑定对话框）
        md_input = page.locator('label.recommended-file-button input[type="file"]')
        md_input.set_input_files(md_path)
        print("[*] 已灌入 markdown，等待分页/图片引用解析…")
        page.wait_for_timeout(1500)

        # 2) 绑图（若对话框弹出，直接 set 第2个隐藏 input = markdown asset folder）
        dlg = page.locator(".markdown-asset-dialog")
        if dlg.count() and dlg.is_visible():
            print("[*] 检测到本地图片引用对话框，自动绑图…")
            hidden = page.locator("input.visually-hidden")
            if assets_dir:
                # markdownAssetFolderInputRef 带 webkitdirectory，必须传目录路径
                hidden.nth(1).set_input_files(assets_dir)
                page.wait_for_timeout(1500)
                try:
                    dlg.wait_for(state="hidden", timeout=8000)
                    print(f"[*] 绑图完成，素材目录: {assets_dir}")
                except Exception:
                    print("[!] 对话框未自动关闭，检查素材名是否与 md 引用匹配")
            else:
                print("[!] 无素材目录，关闭对话框")
                page.locator(".markdown-asset-dialog .dialog-cancel").click()
        else:
            print("[*] 无未绑定图片引用，直接进入排版")

        # 3) 预设排版（切到"排版设置"面板，按分类逐项正确点击）
        page.locator('.workbench-tabs button:has-text("排版设置")').click()
        page.wait_for_timeout(300)

        # 主题（默认 theme 分类已激活）
        page.locator(f".card-theme-option.theme-{theme}").click()
        page.wait_for_timeout(150)

        # 密度 -> 切到 排版(layout) 分类
        page.locator('.settings-category-nav button:has-text("排版")').click()
        page.wait_for_timeout(200)
        page.locator(f'.density-preset-choices button:has-text("{DENSITY_TEXT[density]}")').click()
        page.wait_for_timeout(150)

        # 封面 + 画布比例 -> 切到 画布(canvas) 分类
        page.locator('.settings-category-nav button:has-text("画布")').click()
        page.wait_for_timeout(200)
        page.locator(f'.cover-mode-control button:has-text("{COVER_TEXT[cover]}")').click()
        page.wait_for_timeout(150)
        page.locator(f'button:has-text("{canvas}")').first.click()
        page.wait_for_timeout(300)

        # 4) 导出
        page.locator("button.top-export-button").click()
        page.wait_for_selector(".export-dialog", timeout=10000)
        print("[*] 导出对话框已开，选择清晰度…")
        if hd:
            page.locator('.export-quality-option:has-text("高清原图")').click()
        else:
            page.locator('.export-quality-option:has-text("标准发布")').click()
        page.wait_for_timeout(500)

        # 等待 preflight ready -> 提交按钮可用
        submit = page.locator(".export-dialog .dialog-submit")
        page.wait_for_function(
            "document.querySelector('.export-dialog .dialog-submit') && "
            "!document.querySelector('.export-dialog .dialog-submit').disabled",
            timeout=30000,
        )
        print("[*] 导出前检查通过，开始下载 ZIP…")
        with page.expect_download() as dl_info:
            submit.click()
        download = dl_info.value
        zip_path = os.path.join(out_dir, f"{slug}.zip")
        download.save_as(zip_path)
        print(f"[*] 已下载 ZIP: {zip_path}")

        # 解压
        with zipfile.ZipFile(zip_path) as z:
            z.extractall(out_dir)
        pngs = sorted(glob.glob(os.path.join(out_dir, "*.png")),
                      key=lambda x: os.path.basename(x))
        print(f"[✓] 共导出 {len(pngs)} 张 PNG -> {out_dir}")
        for f in pngs:
            print("   ", os.path.relpath(f, out_dir))
        # 清理 ZIP 残留（产物只留 PNG）
        try:
            os.remove(zip_path)
            print(f"[*] 已清理 ZIP: {zip_path}")
        except OSError:
            pass
        browser.close()
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("md", help="笔记 md 路径（建议传 00X-bioSkills真实试用-<topic>.md）")
    ap.add_argument("--theme", default="notebook",
                    choices=["minimal", "research", "editorial", "notebook"])
    ap.add_argument("--density", default="balanced",
                    choices=["relaxed", "balanced", "compact"])
    ap.add_argument("--cover", default="integrated",
                    choices=["integrated", "standalone", "none"])
    ap.add_argument("--canvas", default="3:4", choices=["3:4", "2:3"])
    ap.add_argument("--hd", action="store_true", help="高清原图(2160x2880) 而非标准发布(1080x1440)")
    ap.add_argument("--assets", default=None, help="素材目录（默认自动探测 content/素材/<stem>）")
    ap.add_argument("--out", default="output/xhs-cards",
                   help="输出根目录，实际产物落在 <out>/<slug>/<theme>-<cover>[-hd]/")
    ap.add_argument("--url", default="http://localhost:5173/app/")
    args = ap.parse_args()
    sys.exit(run(args.md, args.theme, args.density, args.cover,
                 args.canvas, args.hd, args.assets, args.out, args.url))


if __name__ == "__main__":
    main()
