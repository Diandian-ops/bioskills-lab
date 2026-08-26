#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成 bioSkills 试用小红书封面卡：skill 名 + 功能要点。

portrait 3:4 (1080x1440)，白底清爽风格（适配小红书信息流）。
中文用 PingFang SC。

用法：
  python make_cover.py --out content/素材/006-trimming/006-cover.png \\
    --name "alignment-trimming" \\
    --line "比对修剪别乱砍" --line "默认模式只动 gap 列" \\
    --cap "kpic-smart-gap 只动 gap 列" --cap "--log 让修剪可被审计"
"""
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['PingFang SC', 'Songti SC', 'Heiti SC']
plt.rcParams['axes.unicode_minus'] = False

# 小红书友好配色：白底 + 暖红强调（不刺眼、不沉闷）
ACCENT = '#d94a38'      # 暖红（比旧 #b5341f 更亮，适合手机屏）
BG     = '#ffffff'      # 纯白底
TEXT_D = '#1a1a1a'      # 主文字
TEXT_M = '#555555'      # 次要文字
TEXT_L = '#999999'      # 弱化文字


def make_cover(out, name, lines, caps, accent=ACCENT):
    fig, ax = plt.subplots(figsize=(5.4, 7.2), dpi=200)
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)
    ax.axis('off')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # 左侧色条（细，品牌标识）
    ax.add_patch(plt.Rectangle((0.04, 0.88), 0.03, 0.08,
                               color=accent, transform=ax.transAxes))

    # skill 名（等宽）
    name_fs = 30 if len(name) > 15 else 36
    ax.text(0.10, 0.82, name, fontsize=name_fs, color=TEXT_D,
            fontweight='bold', fontfamily='monospace', transform=ax.transAxes)

    # 功能描述行
    y = 0.70
    for ln in lines:
        ax.text(0.10, y, ln, fontsize=16, color=TEXT_M, transform=ax.transAxes)
        y -= 0.065

    # 分隔线（轻）
    y -= 0.02
    ax.plot([0.10, 0.90], [y, y], color='#e8e8e8', lw=1, transform=ax.transAxes)
    y -= 0.05

    # 要点列表（自然语句，轻量分隔符）
    if caps:
        for cap in caps:
            ax.text(0.12, y, '\u00b7', fontsize=16, color=accent,
                    transform=ax.transAxes, va='center')
            ax.text(0.17, y, cap, fontsize=14.5, color=TEXT_M,
                    transform=ax.transAxes, va='center')
            y -= 0.058

    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close(fig)
    print('cover ->', out)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--out', required=True)
    p.add_argument('--name', required=True)
    p.add_argument('--line', action='append', default=[])
    p.add_argument('--cap', action='append', default=[])
    p.add_argument('--accent', default=ACCENT)
    args = p.parse_args()
    make_cover(args.out, args.name, args.line, args.cap, accent=args.accent)
