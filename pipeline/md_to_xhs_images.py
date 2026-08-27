#!/usr/bin/env python3
"""把 bioSkills 站内深拆版 markdown 渲染成竖版 3:4 小红书多卡图。

链路：读 md → 去除 META 注释 → markdown 转 HTML（套竖版 3:4 模板）
     → Playwright 全页截图 → 按 1080×1440(@2x=2160×2880) 切成多张 3:4 图。

用法：
    python pipeline/md_to_xhs_images.py <md_path> [--out DIR]
    python pipeline/md_to_xhs_images.py content/笔记/006-bioSkills真实试用-alignment-trimming.md

依赖（隔离 venv）：markdown / playwright / Pillow，且已 `playwright install chromium`。
"""
import os
import re
import argparse
from PIL import Image
import markdown
from playwright.sync_api import sync_playwright

# BASE 推导：优先环境变量，否则脚本位置向上两级（本仓库根）
BASE = os.environ.get("REDBOOK_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)

# 小红书 3:4，@2x 高清输出（物理像素）
W, H = 1080, 1440
SCALE = 2
PW, PH = W * SCALE, H * SCALE  # 2160 × 2880

TPL = """<!doctype html><html lang="zh"><head><meta charset="utf-8">
<style>
*{{box-sizing:border-box}}
body{{margin:0;background:#fff;color:#1a1a1a;
 font-family:-apple-system,BlinkMacSystemFont,"PingFang SC","Heiti SC","Microsoft YaHei",sans-serif;
 line-height:1.85;font-size:34px;font-weight:400;}}
.wrap{{width:__W__px;margin:0 auto;padding:60px 50px 90px;}}
h1{{font-size:60px;line-height:1.3;font-weight:800;color:#b5341f;margin:0 0 30px;
 border-bottom:6px solid #b5341f;padding-bottom:20px;}}
h2{{font-size:46px;font-weight:800;color:#b5341f;margin:50px 0 20px;
 padding-left:20px;border-left:10px solid #b5341f;}}
h3{{font-size:38px;font-weight:700;color:#222;margin:32px 0 14px;}}
p{{margin:18px 0;}}
ul,ol{{margin:16px 0;padding-left:1.4em;}}
li{{margin:12px 0;}}
strong{{color:#b5341f;font-weight:700;}}
a{{color:#b5341f;}}
code{{background:#f1f1f1;padding:2px 8px;border-radius:6px;font-size:26px;
 font-family:"SF Mono",Menlo,Consolas,monospace;}}
pre{{background:#f6f6f6;border:1px solid #e6e6e6;border-radius:12px;
 padding:22px 24px;font-size:26px;line-height:1.6;
 font-family:"SF Mono",Menlo,Consolas,monospace;white-space:pre-wrap;word-break:break-word;}}
pre code{{background:none;padding:0;font-size:26px;}}
table{{border-collapse:collapse;width:100%;margin:22px 0;font-size:26px;}}
th,td{{border:1px solid #ddd;padding:12px 14px;text-align:left;}}
th{{background:#b5341f;color:#fff;font-weight:700;}}
tr:nth-child(even) td{{background:#faf6f4;}}
img{{display:block;max-width:100%;height:auto;margin:26px auto;border-radius:10px;
 box-shadow:0 2px 10px rgba(0,0,0,.08);}}
blockquote{{border-left:6px solid #ccc;margin:16px 0;padding:8px 18px;color:#666;}}
hr{{border:none;border-top:2px dashed #e0d8d4;margin:34px 0;}}
</style></head><body><div class="wrap">{body}</div></body></html>"""


def strip_meta(text: str) -> str:
    """去除 <!-- ... --> 注释块（META 段等，不进发布稿）。"""
    return re.sub(r"<!--.*?-->", "", text, flags=re.DOTALL)


def md_to_html(md_path: str) -> str:
    text = open(md_path, encoding="utf-8").read()
    text = strip_meta(text)
    # 图片相对路径 ../素材/ -> 绝对 file://，确保 Playwright 能加载本地图
    text = text.replace("../素材/", f"file://{BASE}/content/素材/")
    body = markdown.markdown(text, extensions=["tables", "fenced_code"])
    return TPL.replace("__W__", str(PW)).replace("{body}", body)


def render_and_slice(html_path: str, out_dir: str):
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page(viewport={"width": PW, "height": PH})
        page.goto("file://" + html_path)
        page.wait_for_timeout(900)  # 等图片/字体渲染
        full = os.path.join(out_dir, "_full.png")
        page.screenshot(path=full, full_page=True)
        browser.close()

    im = Image.open(full)
    width, height = im.size
    card_h = PH
    n = max(1, (height + card_h - 1) // card_h)
    paths = []
    for i in range(n):
        # 最后一张贴底，保证内容连续、不出现顶部空白
        top = max(0, height - card_h) if i == n - 1 else i * card_h
        bottom = min(top + card_h, height)
        crop = im.crop((0, top, width, bottom))
        if crop.height < card_h:  # 末张不足一屏，白底补齐
            canvas = Image.new("RGB", (width, card_h), (255, 255, 255))
            canvas.paste(crop, (0, 0))
            crop = canvas
        out = os.path.join(out_dir, f"card_{i + 1:02d}.png")
        crop.save(out)
        paths.append(out)
    return paths


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("md", help="站内深拆版 markdown 路径")
    ap.add_argument("--out", help="输出目录（默认 output/xhs-images/<文件名>/）")
    args = ap.parse_args()

    md_path = os.path.abspath(args.md)
    if not os.path.isfile(md_path):
        raise SystemExit(f"md 不存在: {md_path}")
    basename = os.path.splitext(os.path.basename(md_path))[0]
    out_dir = args.out or os.path.join(BASE, "output", "xhs-images", basename)
    os.makedirs(out_dir, exist_ok=True)

    html_path = os.path.join(out_dir, basename + ".html")
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(md_to_html(md_path))

    paths = render_and_slice(html_path, out_dir)
    print(f"产出 {len(paths)} 张 3:4 图 -> {out_dir}")
    for p in paths:
        print("  ", os.path.basename(p))


if __name__ == "__main__":
    main()
