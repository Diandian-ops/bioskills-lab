#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gate_lint.py — bioSkills 笔记 / 小红书正文 的「可读性 + 雷区」机检器。

把 redbook-bio-note-writer（第十节可读性 Gate）与 redbook-xhs-body
（雷区清单 + 自审 Gate）里**机械化可检**的条目转成规则，扫完输出带行号报告。
ERROR = 硬违规，必须改（编排器会据此拦下，不让半成品进出图/发布）；
WARN = 需人眼复核的语义项（通顺/跳跃/反而误用/裸数字/相邻重复/缩写未释）。

用法：
  python pipeline/gate_lint.py <file.md> [--json] [--quiet]
退出码：有 ERROR → 1；无 ERROR → 0。
"""
import argparse
import json
import re
import sys

# ---- 规则定义 ----------------------------------------------------------------

# 口语隐喻 / 拟人（note-writer 10.1 + 历史踩坑，ERROR）
RE_METAPHOR = re.compile(r"踩坑|这个坑|那些坑|陷阱|撞墙|坑确实存在|坑点|实战踩坑|核心坑|避坑|排坑|会产出多少假")

# 生硬搭配「套…(假设/公式/背景/模型)」但排除「套用」（ERROR）
RE_STIFF = re.compile(r"套(?!用).{0,5}(假设|公式|背景|模型)")

# 小标题说教词（note-writer 10.1 / 第四节，ERROR，仅查标题行）
RE_HEADING_PREACHY = re.compile(r"^(#{1,6})\s.*?(不能|千万别|别用|别套|会产出多少假|不是一刀切)")

# 正文证明感（xhs-body 雷区1，ERROR）
RE_PROOF = re.compile(r"原样跑了一遍|我拿.{0,12}跑了|真实跑过|亲自跑|我跑了一遍|我们跑了一遍|我跑通了")

# 「反而」误用（两 skill 都点名，WARN，需人确认是否真反转）
RE_FANER = re.compile(r"反而")

# 未解释缩写（note-writer 10.1 术语首现必解释，WARN）
ABBR = ["IC", "PID", "JSD", "KL", "PSSM", "HMM", "MSA", "BLOSUM", "PSI-BLAST", "SP-score"]
EXPL_WORDS = ["信息量", "一致性", "保守", "距离", "概率", "矩阵", "打分", "比对"]

# 裸数字上下文标记（xhs-body 雷区2，WARN，仅正文模式）
NUM_CTX = re.compile(r"(bits|列|个|%|％|倍|行|条|张|基因|残基|位点|vs|对比|→|高于|低于|差|相对|更|分别|（常见）|（稀有）|（保守）|常见|稀有|保守|=)")
# 元数据/说明行豁免（这些行里的数字是序号/编号，不算裸数据）
META_LINE = re.compile(r"(用途|文案|标题建议|META|建议|序号|小红书|帖子|正文文本)")


def is_body(path: str) -> bool:
    return path.endswith("-正文.md") or "-小红书-" in path


def find_line(lines, idx):
    return idx + 1


def lint(text, path):
    lines = text.split("\n")
    findings = []  # (severity, line, rule, snippet, suggestion)

    body = is_body(path)

    in_meta = False  # META 块（作者自用建议/元信息）内的行不卡发布正文规则
    for i, ln in enumerate(lines):
        if re.search(r"<!--|META", ln):
            in_meta = True
        if not in_meta:
            # --- ERROR: 口语隐喻 ---
            m = RE_METAPHOR.search(ln)
            if m:
                findings.append(("ERROR", i + 1, "METAPHOR", ln.strip(),
                                 "改用客观陈述，如「该配置存在偏差」「错误配置导致…」"))
            # --- ERROR: 生硬搭配 套… ---
            m = RE_STIFF.search(ln)
            if m:
                findings.append(("ERROR", i + 1, "STIFF_COLLOC", ln.strip(),
                                 "「套…假设/公式」→「套用」或换陈述"))
            # --- ERROR: 小标题说教词 ---
            if RE_HEADING_PREACHY.search(ln):
                findings.append(("ERROR", i + 1, "HEADING_PREACHY", ln.strip(),
                                 "客观描述事实，去掉「不能/别/千万别」等说教词"))
            # --- ERROR(xhs-body): 证明感 ---
            if body and RE_PROOF.search(ln):
                findings.append(("ERROR", i + 1, "PROOF_BOAST", ln.strip(),
                                 "改为「拿示例数据试了下」等轻描淡写表述"))
            # --- WARN: 反而 ---
            if RE_FANER.search(ln):
                findings.append(("WARN", i + 1, "FANER", ln.strip(),
                                 "确认是否真有反转；仅陈述事实（如稀有残基信息量更高）不能用反而"))
            # --- WARN(xhs-body): 裸数字 ---
            if body:
                has_num = bool(re.search(r"\d", ln))
                if has_num and not NUM_CTX.search(ln) and not META_LINE.search(ln) and len(ln.strip()) > 4:
                    findings.append(("WARN", i + 1, "BARE_NUM", ln.strip(),
                                     "数字需带 标签+单位+对照，零背景读者能懂"))
        if re.search(r"/\s*META|-->\s*$", ln):
            in_meta = False

    # --- WARN: 未解释缩写（首现） ---
    for ab in ABBR:
        pat = re.compile(rf"\b{ab}\b")
        for i, ln in enumerate(lines):
            if pat.search(ln):
                pos = ln.find(ab)
                window = ln[max(0, pos - 15): pos + 30]
                explained = ("（" in window) or ("(" in window) or any(w in window for w in EXPL_WORDS)
                if not explained:
                    findings.append(("WARN", i + 1, "ABBR_UNEXPL", ln.strip(),
                                     f"「{ab}」首现需括号注或一句大白话解释"))
                break  # 只看首现

    # --- WARN: 相邻句重复同一观点 ---
    sents = [s.strip() for s in re.split(r"[。！？\n]", text) if len(s.strip()) > 8]
    for a, b in zip(sents, sents[1:]):
        # 跳过代码/公式片段（含 = 或括号密集），避免误报
        if re.search(r"[=\[\]()]", a) and re.search(r"[=\[\]()]", b):
            continue
        def bigrams(s):
            s = re.sub(r"\s+", "", s)
            return set(s[j:j + 2] for j in range(len(s) - 1))
        ga, gb = bigrams(a), bigrams(b)
        if ga and gb:
            jac = len(ga & gb) / len(ga | gb)
            if jac > 0.5:
                findings.append(("WARN", -1, "ADJ_REPEAT", f"{a[:20]}… / {b[:20]}…",
                                 "相邻两句信息重叠过高，每段只推进一个逻辑环节"))

    return findings


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("md")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--quiet", action="store_true", help="只打印 ERROR")
    args = ap.parse_args()

    try:
        with open(args.md, encoding="utf-8") as f:
            text = f.read()
    except OSError as e:
        print(f"[!] 读不了: {e}", file=sys.stderr)
        sys.exit(2)

    findings = lint(text, args.md)
    errors = [f for f in findings if f[0] == "ERROR"]
    warns = [f for f in findings if f[0] == "WARN"]

    if args.json:
        print(json.dumps([
            {"severity": s, "line": l, "rule": r, "snippet": sn, "suggestion": su}
            for (s, l, r, sn, su) in findings
        ], ensure_ascii=False, indent=2))
    else:
        mode = "小红书正文" if is_body(args.md) else "笔记(真实试用/通用)"
        print(f"=== gate_lint · {mode} · {args.md} ===")
        if not findings:
            print("✓ 无机械违规")
        else:
            show = errors if args.quiet else findings
            for sev, line, rule, snip, sug in show:
                loc = f"L{line}" if line > 0 else "—"
                print(f"[{sev}] {loc} | {rule}\n    > {snip}\n    → {sug}")
        print(f"\n统计: ERROR={len(errors)}  WARN={len(warns)}")

    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
