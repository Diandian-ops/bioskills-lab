#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成 bioSkills 试用小红书封面卡（品牌标题卡）：skill 名 + 一句话功能 + 核心能力点。

portrait 3:4 (1080x1440)，编辑档案风暖纸 #eee8db + 红强调 #b5341f，中文用 PingFang SC。
这是「小红书精简模板」的封面生成器，以后每个 trial 直接套。

用法：
  python make_cover.py --out content/素材/004-pairwise/004-cover.png \
    --name "pairwise-alignment" --tag "bioSkills 真实试用" \
    --line "两条序列的最优比对" --line "教你怎么用 Biopython 跑通" \
    --cap "选矩阵：BLOSUM62" --cap "设 gap：open=-11 / extend=-1" --cap "算 PID / 经验 p 值"
"""
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['PingFang SC', 'Songti SC', 'Heiti SC']
plt.rcParams['axes.unicode_minus'] = False


def make_cover(out, name, tag, lines, caps, accent='#b5341f', paper='#eee8db'):
    fig, ax = plt.subplots(figsize=(5.4, 7.2), dpi=200)
    fig.patch.set_facecolor(paper)
    ax.set_facecolor(paper)
    ax.axis('off')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # 顶部红色短条 + tag
    ax.add_patch(plt.Rectangle((0.08, 0.93), 0.12, 0.012,
                               color=accent, transform=ax.transAxes))
    ax.text(0.08, 0.875, tag, fontsize=15, color=accent,
            fontweight='bold', transform=ax.transAxes)

    # skill 名（等宽，长名自动缩小）
    name_fs = 30 if len(name) > 15 else 38
    ax.text(0.08, 0.73, name, fontsize=name_fs, color='#1a1a1a',
            fontweight='bold', fontfamily='monospace', transform=ax.transAxes)

    # 一句话功能（多行）
    y = 0.62
    for ln in lines:
        ax.text(0.08, y, ln, fontsize=17, color='#3a3a3a', transform=ax.transAxes)
        y -= 0.07

    # 分隔线
    y -= 0.01
    ax.plot([0.08, 0.92], [y, y], color=accent, lw=2, transform=ax.transAxes)
    y -= 0.05

    # 核心能力点
    ax.text(0.08, y, '核心能力', fontsize=13, color=accent,
            fontweight='bold', transform=ax.transAxes)
    y -= 0.055
    for cap in caps:
        ax.text(0.10, y, '●', fontsize=11, color=accent,
                transform=ax.transAxes, va='center')
        ax.text(0.15, y, cap, fontsize=15, color='#222',
                transform=ax.transAxes, va='center')
        y -= 0.055

    # 底部脚注
    ax.text(0.08, 0.05, '真实运行 · 非演示稿', fontsize=11, color='#888',
            style='italic', transform=ax.transAxes)

    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close(fig)
    print('cover ->', out)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--out', required=True)
    p.add_argument('--name', required=True)
    p.add_argument('--tag', default='bioSkills 真实试用')
    p.add_argument('--line', action='append', default=[])
    p.add_argument('--cap', action='append', default=[])
    p.add_argument('--accent', default='#b5341f')
    args = p.parse_args()
    make_cover(args.out, args.name, args.tag, args.line, args.cap, accent=args.accent)
