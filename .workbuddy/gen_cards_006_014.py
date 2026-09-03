#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""前台串行编排器：对 alignment 006/007/010 + read-alignment 011-014 的卡片摘录稿
跑 md2card_patch.py（4 主题）。md2card 单实例，必须串行；已完成（有 PNG 且无 .nometa.md）自动跳过。

用法：先确保 md2card dev server 已起（content/库/md2card 下 `npm run dev`，占用 5173），
再 `python gen_cards_006_014.py`。

摘录稿（card-ms）需事先由 lead 从深拆笔记生成并落到：
  output/xhs-site/card-ms/alignment/<slug>.md
  output/xhs-site/card-ms/read-alignment/<slug>.md
素材目录：
  content/素材/alignment/<slug>/
  content/素材/read-alignment/<slug>/
"""
import os, subprocess, glob, sys

REPO = r"D:/1.WorkDir/RedBook"
VENV = r"C:/Users/Admin/.workbuddy/binaries/python/envs/default/Scripts/python.exe"
PATCH = os.path.join(REPO, "output/xhs-site/scripts/md2card_patch.py")
OUT = os.path.join(REPO, "output/xhs-site/xhs-cards")
URL = "http://localhost:5173/app/"

# (cat, slug) —— slug 与 build_lab_site.py TRIALS 中一致
JOBS = [
    ("alignment",       "006-alignment-trimming"),
    ("alignment",       "007-multiple-alignment"),
    ("alignment",       "010-structural-alignment"),
    ("read-alignment",  "011-bowtie2-alignment"),
    ("read-alignment",  "012-bwa-alignment"),
    ("read-alignment",  "013-hisat2-alignment"),
    ("read-alignment",  "014-star-alignment"),
]
THEMES = ["notebook", "minimal", "research", "editorial"]

CARD_MS_ROOT = os.path.join(REPO, "output/xhs-site/card-ms")
ASSETS_ROOT = os.path.join(REPO, "content/素材")


def is_done(out_dir):
    pngs = glob.glob(os.path.join(out_dir, "*.png"))
    nometa = glob.glob(os.path.join(out_dir, ".*.nometa.md"))
    return bool(pngs) and not nometa


gen = ok = skip = fail = missing = 0
for cat, sk in JOBS:
    md = os.path.join(CARD_MS_ROOT, cat, sk + ".md")
    assets = os.path.join(ASSETS_ROOT, cat, sk)
    if not os.path.isfile(md):
        print(f"MISSING card-ms {cat}/{sk} -> {md}")
        missing += 1
        continue
    if not os.path.isdir(assets):
        print(f"WARN no assets dir {assets} (cards may lack figures)")
    for th in THEMES:
        variant = f"{th}-standalone-hd"
        out_dir = os.path.join(OUT, cat, sk, variant)
        if is_done(out_dir):
            print(f"SKIP {cat}/{sk}/{th} (already {len(glob.glob(os.path.join(out_dir,'*.png')))} png)")
            skip += 1
            continue
        cmd = [VENV, PATCH, md, "--theme", th, "--cover", "standalone",
               "--out", OUT, "--url", URL, "--hd", "--assets", assets]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, errors="replace", timeout=180)
            pngs = glob.glob(os.path.join(out_dir, "*.png"))
            if r.returncode == 0 and pngs and not glob.glob(os.path.join(out_dir, ".*.nometa.md")):
                print(f"OK   {cat}/{sk}/{th} -> {len(pngs)} png")
                ok += 1
            else:
                print(f"FAIL {cat}/{sk}/{th} rc={r.returncode} png={len(pngs)}")
                print((r.stdout + r.stderr)[-400:])
                fail += 1
        except subprocess.TimeoutExpired:
            print(f"FAIL {cat}/{sk}/{th} TIMEOUT")
            fail += 1
        gen += 1

print(f"\nSUMMARY gen={gen} ok={ok} skip={skip} fail={fail} missing={missing}")
sys.exit(1 if (fail or missing) else 0)
