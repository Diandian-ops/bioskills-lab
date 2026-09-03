#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""前台串行编排器：对 alignment-files 029-037 的卡片摘录稿跑 md2card_patch.py（4 主题）。
md2card 单实例，必须串行；已完成（有 PNG 且无 .nometa.md）自动跳过。"""
import os, subprocess, glob, sys

REPO = r"D:/1.WorkDir/RedBook"
VENV = r"C:/Users/Admin/.workbuddy/binaries/python/envs/default/Scripts/python.exe"
PATCH = os.path.join(REPO, "output/xhs-site/scripts/md2card_patch.py")
CARD_MS = os.path.join(REPO, "output/xhs-site/card-ms/alignment-files")
ASSETS = os.path.join(REPO, "content/素材/alignment-files")
OUT = os.path.join(REPO, "output/xhs-site/xhs-cards")
URL = "http://localhost:5173/app/"

SKILLS = [
    "029-alignment-sorting", "030-indexing", "031-filtering", "032-validation",
    "033-bam-statistics", "034-pileup-generation", "035-duplicate-handling",
    "036-reference-operations", "037-amplicon-clipping",
]
THEMES = ["notebook", "minimal", "research", "editorial"]

def is_done(out_dir):
    pngs = glob.glob(os.path.join(out_dir, "*.png"))
    nometa = glob.glob(os.path.join(out_dir, ".*.nometa.md"))
    return bool(pngs) and not nometa

gen = ok = skip = fail = 0
for sk in SKILLS:
    md = os.path.join(CARD_MS, sk + ".md")
    assets = os.path.join(ASSETS, sk)
    for th in THEMES:
        variant = f"{th}-standalone-hd"
        out_dir = os.path.join(OUT, "alignment-files", sk, variant)
        if is_done(out_dir):
            print(f"SKIP {sk}/{th} (already {len(glob.glob(os.path.join(out_dir,'*.png')))} png)")
            skip += 1
            continue
        cmd = [VENV, PATCH, md, "--theme", th, "--cover", "standalone",
               "--out", OUT, "--url", URL, "--hd", "--assets", assets]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, errors="replace", timeout=180)
            pngs = glob.glob(os.path.join(out_dir, "*.png"))
            if r.returncode == 0 and pngs:
                print(f"OK   {sk}/{th} -> {len(pngs)} png")
                ok += 1
            else:
                print(f"FAIL {sk}/{th} rc={r.returncode} png={len(pngs)}")
                print((r.stdout + r.stderr)[-400:])
                fail += 1
        except subprocess.TimeoutExpired:
            print(f"FAIL {sk}/{th} TIMEOUT")
            fail += 1
        gen += 1

print(f"\nSUMMARY gen={gen} ok={ok} skip={skip} fail={fail}")
sys.exit(1 if fail else 0)
