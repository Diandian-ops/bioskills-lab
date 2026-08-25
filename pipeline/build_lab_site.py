#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成 bioSkills 真实试用实验室 - 多文件静态站点（多层级侧栏 + 内容居中 + 移动抽屉）。"""
import html, os, shutil, re, markdown as _md

BASE = "/Users/zhangdiandian/RedBook"
BIO = BASE + "/content/库/bioSkills"
IMG004 = BASE + "/content/素材/004-pairwise"
IMG005 = BASE + "/content/素材/005-msa-statistics"
EXP004 = BASE + "/pipeline/004-pairwise-redo"
EXP005 = BASE + "/pipeline/005-msa-statistics"
SITE = BASE + "/output/bioSkills-site"
ASSETS = SITE + "/assets"

def code(path):
    with open(path, "r", encoding="utf-8") as f:
        return html.escape(f.read().rstrip("\n"))

def load_note(md_path):
    """读取 笔记（单一事实源）-> (html, meta, needed_imgs)。
    剥离 META 注释块（并解析 标题/副标题），YAML frontmatter 作兜底；
    图片相对路径改写并收集源文件，表格/fenced 代码扩展，代码块包 .codewrap。"""
    raw = open(md_path, "r", encoding="utf-8").read()
    meta = {}
    m = re.search(r"<!--\s*META(.*?)/META\s*-->", raw, flags=re.S)
    if m:
        for line in m.group(1).splitlines():
            if ":" in line:
                k, v = line.split(":", 1)
                meta[k.strip()] = v.strip()
        raw = raw[:m.start()] + raw[m.end():]
    elif raw.startswith("---"):
        raw = re.sub(r"^---.*?\n---\n", "", raw, count=1, flags=re.S)
    needed = []
    def repl_img(mm):
        rel = mm.group(2)
        src = os.path.normpath(os.path.join(os.path.dirname(md_path), rel))
        if os.path.exists(src):
            needed.append(src)
            return "![%s](../assets/%s)" % (mm.group(1), os.path.basename(src))
        return mm.group(0)
    raw = re.sub(r"!\[([^\]]*)\]\(([^)]+)\)", repl_img, raw)
    html = _md.markdown(raw, extensions=["tables", "fenced_code"])
    html = html.replace("<pre>", '<div class="codewrap"><pre>').replace("</pre>", "</pre></div>")
    # 笔记里的图是裸 <img>（在 <p> 内），包进 <figure> 以复用 figure 边框样式 + 放大交互
    html = re.sub(r"<p>\s*(<img[^>]*>)\s*</p>", r"<figure>\1</figure>", html)
    # META 存在时 masthead 已渲染 标题(h1)+副标题，正文自带的 H1+引导段是给 DingCard 独立渲染用的
    # → 站点上去重（DingCard 读原始 md 不受影响）；无 标题 时保留正文 H1 作兜底
    # 引导段可能含内嵌标签(<strong>/<code>)，用 .*? 而非 [^<]* 以兼容
    if meta.get("标题"):
        html = re.sub(r"^<h1>.*?</h1>\s*", "", html, count=1, flags=re.S)
        html = re.sub(r"^<p>.*?</p>\s*(<hr>\s*)?", "", html, count=1, flags=re.S)
    return html, meta, needed


cats = []
for name in sorted(os.listdir(BIO)):
    p = os.path.join(BIO, name)
    if os.path.isdir(p) and not name.startswith('.'):
        n = len([d for d in os.listdir(p) if os.path.isdir(os.path.join(p, d)) and not d.startswith('.')])
        if n > 0:
            cats.append((name, n))
DONE = {"alignment"}
NCAT = len(cats)
NTOTAL = sum(n for _, n in cats)

# ---------- 领域分类法（domain taxonomy）：64 category -> 10 领域 ----------
DOMAIN_META = [
    ("序列基础与读取", "#6366f1"),
    ("转录组与差异",   "#0ea5e9"),
    ("单细胞与空间",   "#14b8a6"),
    ("表观与调控",     "#f59e0b"),
    ("基因组与变异",   "#ef4444"),
    ("蛋白与结构",     "#8b5cf6"),
    ("微生物与生态",   "#22c55e"),
    ("多组学与系统",   "#ec4899"),
    ("免疫与临床",     "#f97316"),
    ("数据与工具",     "#64748b"),
]
DOMAIN_OF = {
    "alignment":"序列基础与读取","alignment-files":"序列基础与读取","read-alignment":"序列基础与读取","read-qc":"序列基础与读取","primer-design":"序列基础与读取","restriction-analysis":"序列基础与读取","sequence-io":"序列基础与读取","sequence-manipulation":"序列基础与读取",
    "alternative-splicing":"转录组与差异","differential-expression":"转录组与差异","expression-matrix":"转录组与差异","ribo-seq":"转录组与差异","rna-quantification":"转录组与差异","small-rna-seq":"转录组与差异","epitranscriptomics":"转录组与差异","rna-structure":"转录组与差异",
    "flow-cytometry":"单细胞与空间","imaging-mass-cytometry":"单细胞与空间","single-cell":"单细胞与空间","spatial-transcriptomics":"单细胞与空间","tcr-bcr-analysis":"单细胞与空间",
    "atac-seq":"表观与调控","chip-seq":"表观与调控","clip-seq":"表观与调控","methylation-analysis":"表观与调控","hi-c-analysis":"表观与调控","gene-regulatory-networks":"表观与调控",
    "comparative-genomics":"基因组与变异","copy-number":"基因组与变异","crispr-screens":"基因组与变异","genome-annotation":"基因组与变异","genome-assembly":"基因组与变异","genome-engineering":"基因组与变异","genome-intervals":"基因组与变异","long-read-sequencing":"基因组与变异","phasing-imputation":"基因组与变异","phylogenetics":"基因组与变异","population-genetics":"基因组与变异","variant-calling":"基因组与变异",
    "chemoinformatics":"蛋白与结构","proteomics":"蛋白与结构","structural-biology":"蛋白与结构",
    "ecological-genomics":"微生物与生态","metagenomics":"微生物与生态","microbiome":"微生物与生态",
    "causal-genomics":"多组学与系统","metabolomics":"多组学与系统","multi-omics-integration":"多组学与系统","pathway-analysis":"多组学与系统","systems-biology":"多组学与系统","temporal-genomics":"多组学与系统",
    "clinical-biostatistics":"免疫与临床","clinical-databases":"免疫与临床","epidemiological-genomics":"免疫与临床","immunoinformatics":"免疫与临床","liquid-biopsy":"免疫与临床",
    "clawhub-installer":"数据与工具","data-visualization":"数据与工具","database-access":"数据与工具","machine-learning":"数据与工具","reporting":"数据与工具","workflow-management":"数据与工具","workflows":"数据与工具","experimental-design":"数据与工具",
}
dmap = {}
for name, n in cats:
    dmap.setdefault(DOMAIN_OF.get(name, "数据与工具"), []).append((name, n))
TREE = {"name":"bioSkills","children":[]}
for dname, dcolor in DOMAIN_META:
    items = sorted(dmap.get(dname, []), key=lambda x: -x[1])
    TREE["children"].append({
        "name": dname, "color": dcolor,
        "children": [{"name":c,"value":n,"done": c in DONE,
                      "link": ("alignment/index.html" if c=="alignment" else None),
                      "skills": sorted(d for d in os.listdir(os.path.join(BIO,c)) if os.path.isdir(os.path.join(BIO,c, d)) and not d.startswith('.'))} for c,n in items]
    })
TREE_JSON = __import__("json").dumps(TREE, ensure_ascii=False)

# 复现代码由笔记自身承载（笔记内已含复现结果与图示），不再单独注入 pipeline 下的 .py

chips = []
for name, n in cats:
    cls = "chip chip-done" if name in DONE else "chip"
    chips.append('<span class="%s" title="%s">%s <b>x%d</b></span>' % (cls, name, name, n))
CAT_MAP = "".join(chips)

ALIGN_SUBS = [("pairwise", "pairwise-alignment", "alignment/pairwise-alignment.html"),
              ("msa", "msa-statistics", "alignment/msa-statistics.html")]

def sidebar(active, prefix=""):
    s = '<aside class="side" id="side"><div class="brand">BIO / LAB</div>'
    s += '<nav class="sidenav" aria-label="站点导航">'
    s += '<a class="%s" href="%sindex.html">实验室首页</a>' % ('cur' if active == 'index' else '', prefix)
    # 按领域全景分组（中文领域头 + 配色点）
    for dname, dcolor in DOMAIN_META:
        items = dmap.get(dname, [])
        dom_active = any(c == active or (active in ('pairwise', 'msa') and c == 'alignment') for c, _ in items)
        dtot = sum(n for _, n in items)
        s += '<details class="navgrp%s" data-key="D:%s"><summary><span class="dot" style="background:%s"></span>%s<span class="cnt">×%d</span></summary>' % (
            ' open' if dom_active else '', dname, dcolor, dname, dtot)
        for cname, cn in items:
            is_active = (cname == active) or (active in ('pairwise', 'msa') and cname == 'alignment')
            if cname in DONE:
                s += '<details class="navsub%s" data-key="C:%s"><summary class="%s">%s<span class="cnt">×%d</span></summary>' % (
                    ' open' if is_active else '', cname, 'cur' if is_active else '', cname, cn)
                for key, label, href in ALIGN_SUBS:
                    s += '<a class="sub %s" href="%s%s">%s</a>' % ('cur' if active == key else '', prefix, href, label)
                sk = sorted(d for d in os.listdir(os.path.join(BIO, cname))
                            if os.path.isdir(os.path.join(BIO, cname, d)) and not d.startswith('.'))
                for sn in sk:
                    if sn in ('pairwise-alignment', 'msa-statistics'):
                        continue
                    s += '<span class="sub nav-dis">%s</span>' % sn
                s += '</details>'
            elif is_active:
                s += '<details class="navsub open" data-key="C:%s"><summary class="cur">%s<span class="cnt">×%d</span></summary>' % (cname, cname, cn)
                sk = sorted(d for d in os.listdir(os.path.join(BIO, cname))
                            if os.path.isdir(os.path.join(BIO, cname, d)) and not d.startswith('.'))
                for sn in sk[:24]:
                    s += '<span class="sub nav-dis">%s</span>' % sn
                if len(sk) > 24:
                    s += '<span class="sub nav-dis">… 其余 %d 个</span>' % (len(sk) - 24)
                s += '</details>'
            else:
                s += '<div class="navcat"><span class="dot" style="background:%s"></span>%s<span class="cnt">×%d</span></div>' % (dcolor, cname, cn)
        s += '</details>'
    s += '</nav></aside>'
    return s

STYLE = """/* 文档站：深色多层级侧栏 + 浅色内容(居中) + 靛蓝强调 */
:root{ --bg:#f6f7fb; --side:#0f1623; --side2:#182338; --side-text:#c3ccdb; --side-dim:#7c879b;
  --ink:#15202e; --ink-soft:#475569; --muted:#7b8794; --line:#e6e9f0; --accent:#6d5efc; --accent-d:#5338e0;
  --accent-soft:#efecfe; --card:#ffffff; --green:#16a34a; }
*{box-sizing:border-box;}
html{scroll-behavior:smooth;}
body{margin:0;background:linear-gradient(180deg,#fcfcfd 0%,#f4f5f8 100%);color:var(--ink);
  font-family:-apple-system,"PingFang SC","Noto Sans CJK SC",sans-serif;line-height:1.75;}
.side{position:fixed;top:0;left:0;width:264px;height:100vh;overflow-y:auto;background:linear-gradient(168deg,#111b2c 0%,#0c1320 60%,#0a0f1a 100%);border-right:1px solid rgba(109,94,252,.22);padding:22px 12px 40px;z-index:70;}
.brand{font:700 16px/1 "Courier New",monospace;letter-spacing:.18em;color:#c7d2fe;margin:4px 8px 18px;}
.sidenav{display:flex;flex-direction:column;gap:2px;font-size:14px;}
.sidenav>a{color:var(--side-text);text-decoration:none;padding:8px 12px;border-radius:6px;}
.sidenav>a:hover{background:var(--side2);color:#fff;}
.sidenav>a.cur{background:rgba(99,102,241,.20);color:#c7d2fe;font-weight:600;}
.navgrp{margin:4px 0;}
.navgrp>summary{list-style:none;cursor:pointer;padding:8px 12px;border-radius:6px;color:#cbd5e1;font-weight:600;font-size:14px;user-select:none;}
.navgrp>summary::-webkit-details-marker{display:none;}
.navgrp>summary::before{content:"▸";display:inline-block;margin-right:8px;font-size:10px;transition:transform .15s;color:var(--accent);}
.navgrp[open]>summary::before{transform:rotate(90deg);}
.navgrp>summary:hover{background:var(--side2);}
.navgrp>summary.cur{color:#c7d2fe;}
.sub{display:block;padding:7px 12px 7px 32px;font-size:13px;color:var(--side-dim);text-decoration:none;border-radius:6px;}
.sub:hover{background:var(--side2);color:#fff;}
.sub.cur{background:rgba(99,102,241,.20);color:#c7d2fe;font-weight:600;}
.nav-dis{display:block;padding:6px 12px 6px 32px;font-size:13px;color:#566072;opacity:.55;}
.cnt{color:#4b5563;font-weight:400;}
.navgrp>summary{display:flex;align-items:center;}
.navgrp>summary .dot{flex:none;width:9px;height:9px;border-radius:50%;margin-right:9px;}
.navgrp>summary .cnt{margin-left:auto;margin-right:4px;}
.navsub{margin:2px 0 3px 8px;border-left:2px solid var(--side2);}
.navsub>summary{list-style:none;cursor:pointer;display:flex;align-items:center;padding:6px 12px 6px 14px;font-size:13px;color:#cbd5e1;border-radius:6px;user-select:none;}
.navsub>summary::-webkit-details-marker{display:none;}
.navsub>summary::before{content:"▸";display:inline-block;margin-right:6px;font-size:9px;transition:transform .15s;color:var(--accent);}
.navsub[open]>summary::before{transform:rotate(90deg);}
.navsub>summary.cur{color:#c7d2fe;font-weight:600;}
.navsub>summary:hover{background:var(--side2);}
.navsub>summary .cnt{margin-left:auto;margin-right:2px;}
.navsub .sub{padding-left:20px;font-size:12.5px;}
.navcat{display:flex;align-items:center;padding:6px 12px 6px 14px;font-size:13px;color:var(--side-dim);}
.navcat .dot{flex:none;width:7px;height:7px;border-radius:50%;margin-right:9px;opacity:.6;}
.navcat .cnt{margin-left:auto;margin-right:2px;}
.content{margin-left:264px;}
.wrap{max-width:760px;margin:0 auto;padding:38px 28px 90px;}
.menubtn{display:none;position:fixed;top:12px;left:12px;z-index:90;width:42px;height:42px;border:none;border-radius:8px;
  background:var(--accent);color:#fff;font-size:18px;cursor:pointer;box-shadow:0 2px 10px rgba(0,0,0,.2);}
.crumb{font-size:13px;color:var(--muted);margin:0 0 18px;}
.crumb a{color:var(--accent);text-decoration:none;}
.masthead{border-bottom:2px solid var(--accent);padding:0 0 22px;margin-bottom:10px;}
.kicker{font:600 12px/1 "Courier New",monospace;letter-spacing:.3em;color:var(--accent);text-transform:uppercase;}
h1{font-family:"Songti SC","SimSun",serif;font-size:33px;font-weight:700;margin:10px 0 6px;color:#0f172a;line-height:1.25;}
h2{font-family:"Songti SC","SimSun",serif;font-size:25px;font-weight:700;margin:8px 0 4px;color:#0f172a;scroll-margin-top:24px;}
.sub{color:var(--muted);font-size:15px;margin:0;}
.lead{margin:18px 0 0;font-size:16px;color:var(--ink-soft);}
.stats{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:24px 0;}
.stat{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:16px;text-align:center;box-shadow:0 1px 2px rgba(0,0,0,.03);}
.stat .v{font:700 26px/1 "Courier New",monospace;color:var(--accent-d);}
.stat .l{color:var(--muted);font-size:12px;margin-top:6px;}
.map{display:flex;flex-wrap:wrap;gap:7px;margin:16px 0;}
.chip{background:var(--card);border:1px solid var(--line);border-radius:14px;padding:5px 11px;font-size:12.5px;color:var(--muted);}
.chip b{color:#475569;font-weight:700;}
.chip-done{background:var(--accent);border-color:var(--accent);color:#fff;}
.chip-done b{color:#e0e7ff;font-weight:700;}
.tm-wrap{position:relative;margin:18px 0;}
#tm svg{display:block;}
.tm-tip{position:absolute;display:none;pointer-events:none;z-index:6;background:#0f172a;color:#e2e8f0;font-size:12.5px;line-height:1.55;padding:8px 11px;border-radius:8px;box-shadow:0 8px 24px rgba(0,0,0,.28);max-width:220px;}
.tm-tip b{color:#fff;}
.tm-legend{display:flex;flex-wrap:wrap;gap:8px 16px;margin:14px 0 2px;}
.tm-lg{display:inline-flex;align-items:center;gap:6px;font-size:12.5px;color:var(--ink-soft);}
.tm-lg i{width:12px;height:12px;border-radius:3px;display:inline-block;}
.tm-lg b{color:#0f172a;font-family:"Courier New",monospace;}
.sklist{border:1px solid var(--line);background:var(--card);border-radius:12px;overflow:hidden;margin:18px 0;}
.sklist a,.sklist div{padding:12px 16px;border-bottom:1px solid var(--line);display:flex;justify-content:space-between;align-items:center;text-decoration:none;color:var(--ink);}
.sklist a:last-child,.sklist div:last-child{border-bottom:none;}
.sklist a:hover{background:var(--accent-soft);}
.tag{display:inline-block;font:600 11px/1 "Courier New",monospace;padding:3px 8px;border-radius:4px;letter-spacing:.05em;}
.tag-done{background:var(--accent);color:#fff;}
.tag-use{background:var(--green);color:#fff;}
.tag-todo{background:#e2e8f0;color:#64748b;}
section{border-top:1px solid var(--line);margin-top:42px;padding-top:28px;scroll-margin-top:24px;}
.sec-kicker{font:600 12px/1 "Courier New",monospace;letter-spacing:.28em;color:var(--accent);text-transform:uppercase;}
.sec-one{color:var(--muted);font-size:15px;margin:0 0 18px;}
h3{font-size:17px;margin:26px 0 8px;color:#0f172a;scroll-margin-top:24px;}
p{margin:10px 0;color:var(--ink-soft);}
ul{margin:8px 0;padding-left:22px;} li{margin:5px 0;color:var(--ink-soft);}
.tree{background:#1e2430;color:#d6deeb;font:13px/1.6 "Courier New",monospace;padding:14px 16px;border-radius:8px;overflow-x:auto;white-space:pre;}
.cards{display:grid;grid-template-columns:1fr 1fr;gap:14px;margin:16px 0;}
.card{border:1px solid var(--line);border-radius:12px;padding:14px 16px;background:var(--card);box-shadow:0 1px 2px rgba(0,0,0,.03);}
.card .n{font:700 22px/1 "Courier New",monospace;color:var(--accent-d);}
.card .t{font-weight:600;margin:6px 0 4px;color:#0f172a;}
.card .d{color:var(--muted);font-size:13px;margin:0;}
figure{margin:22px 0;text-align:center;}
figure img, .content img{max-width:100%;border:1px solid var(--line);border-radius:8px;background:#fff;box-shadow:0 1px 3px rgba(0,0,0,.05);cursor:zoom-in;transition:filter .2s ease, transform .2s ease;}
figure:hover img{filter:brightness(1.03);}
figcaption{color:var(--muted);font-size:12.5px;margin-top:8px;}
.codewrap{position:relative;margin:12px 0;}
pre{margin:0;border-radius:8px;}
pre code{font:13px/1.55 "Courier New",monospace;padding:14px 16px;display:block;border-radius:8px;}
.copy{position:absolute;top:8px;right:8px;z-index:3;width:30px;height:30px;display:inline-flex;align-items:center;justify-content:center;background:rgba(15,23,42,.78);color:#fff;border:none;border-radius:7px;cursor:pointer;opacity:0;transition:opacity .15s ease, background .15s ease;backdrop-filter:blur(2px);}
.codewrap:hover .copy,.tbl-wrap:hover .copy{opacity:1;}
.copy:hover{background:var(--accent-d);}
.copy.ok{background:var(--green);}
@media (hover:none){ .codewrap .copy,.tbl-wrap .copy{opacity:1;} }
.tbl-wrap{position:relative;margin:16px 0;overflow-x:auto;border:1px solid var(--line);border-radius:12px;background:var(--card);}
.tbl-wrap table{width:100%;border-collapse:collapse;font-size:14px;margin:0;}
.tbl-wrap th,.tbl-wrap td{padding:9px 12px;border-bottom:1px solid var(--line);text-align:left;color:var(--ink-soft);}
.tbl-wrap th{background:var(--accent-soft);color:var(--accent-d);font-weight:600;}
.tbl-wrap tr:last-child td{border-bottom:none;}
.res{border:1px solid var(--line);background:var(--accent-soft);border-radius:12px;padding:13px 16px;font-size:14px;color:var(--ink);}
.res b{color:var(--accent-d);}
.note{background:var(--card);border:1px solid var(--line);border-left:4px solid var(--accent);border-radius:12px;padding:13px 16px;font-size:14px;color:var(--ink-soft);}
.todo{border:1px dashed var(--accent);border-radius:12px;padding:16px 18px;background:#fafaff;}
footer{margin-top:60px;border-top:1px solid var(--line);padding-top:22px;color:var(--muted);font-size:13px;}
.totop{position:fixed;right:22px;bottom:22px;z-index:60;width:42px;height:42px;border-radius:50%;border:none;background:var(--accent);color:#fff;font-size:18px;cursor:pointer;display:none;box-shadow:0 2px 10px rgba(0,0,0,.18);}
/* (modal 弹窗样式已随回退移除) */
.lb{position:fixed;inset:0;z-index:200;background:rgba(15,23,42,.82);display:none;align-items:center;justify-content:center;padding:24px;cursor:zoom-out;backdrop-filter:blur(3px);}
.lb.open{display:flex;}
.lb img{max-width:94vw;max-height:92vh;border-radius:10px;box-shadow:0 30px 80px rgba(0,0,0,.5);background:#fff;}
@media(max-width:900px){
  .side{transform:translateX(-100%);transition:transform .22s ease;box-shadow:4px 0 24px rgba(0,0,0,.3);}
  .side.open{transform:none;}
  .content{margin-left:0;}
  .wrap{padding:64px 18px 80px;}
  .stats{grid-template-columns:repeat(2,1fr);}
  .cards{grid-template-columns:1fr;}
  .menubtn{display:block;}
}
"""
COPY_JS = """<script src="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/highlight.min.js"></script>
<script>
  hljs.highlightAll();
  var COPY_ICON='<svg viewBox="0 0 24 24" width="15" height="15" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><rect x="9" y="9" width="11" height="11" rx="2"/><path d="M5 15V5a2 2 0 0 1 2-2h10"/></svg>';
  var OK_ICON='<svg viewBox="0 0 24 24" width="15" height="15" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"><polyline points="20 6 9 17 4 12"/></svg>';
  function mkCopy(){
    var b=document.createElement('button');
    b.type='button'; b.className='copy'; b.setAttribute('aria-label','复制'); b.innerHTML=COPY_ICON;
    return b;
  }
  function bindCopy(btn, getText){
    btn.addEventListener('click', function(e){
      e.stopPropagation();
      navigator.clipboard.writeText(getText()).then(function(){
        btn.innerHTML=OK_ICON; btn.classList.add('ok');
        setTimeout(function(){ btn.innerHTML=COPY_ICON; btn.classList.remove('ok'); }, 1400);
      });
    });
  }
  // 代码块：图标复制按钮
  document.querySelectorAll('pre code').forEach(function(block){
    var wrap=block.parentNode, btn=mkCopy();
    bindCopy(btn, function(){ return block.textContent; });
    wrap.appendChild(btn);
  });
  // 表格：图标复制按钮（复制为 TSV，可直接贴入 Excel）
  document.querySelectorAll('table').forEach(function(tbl){
    var wrap=document.createElement('div'); wrap.className='tbl-wrap';
    tbl.parentNode.insertBefore(wrap, tbl); wrap.appendChild(tbl);
    var btn=mkCopy();
    bindCopy(btn, function(){
      return [].slice.call(tbl.querySelectorAll('tr')).map(function(tr){
        return [].slice.call(tr.children).map(function(td){ return td.innerText.replace(/\\t/g,' ').replace(/\\n/g,' '); }).join('\\t');
      }).join('\\n');
    });
    wrap.appendChild(btn);
  });
  // 图片点击放大（lightbox）
  var lb=document.createElement('div'); lb.className='lb'; lb.innerHTML='<img alt="放大图">';
  document.body.appendChild(lb);
  var lbImg=lb.querySelector('img');
  document.querySelectorAll('.content img').forEach(function(im){
    im.addEventListener('click', function(){ lbImg.src=im.src; lb.classList.add('open'); });
  });
  lb.addEventListener('click', function(){ lb.classList.remove('open'); });
  document.addEventListener('keydown', function(e){ if(e.key==='Escape') lb.classList.remove('open'); });
  // 回到顶部 / 移动抽屉 / 侧栏展开状态持久化
  var tb=document.getElementById('totop');
  window.addEventListener('scroll',function(){tb.style.display=window.scrollY>600?'block':'none';});
  tb.addEventListener('click',function(){window.scrollTo({top:0,behavior:'smooth'});});
  var mb=document.getElementById('menubtn'), sd=document.getElementById('side');
  if(mb){ mb.addEventListener('click',function(){ sd.classList.toggle('open'); }); }
  document.querySelectorAll('.sidenav a').forEach(function(a){ a.addEventListener('click',function(){ sd.classList.remove('open'); }); });
  (function(){
    var KEY='bioLabNavState';
    function load(){ try{ return JSON.parse(localStorage.getItem(KEY)||'{}'); }catch(e){ return {}; } }
    function save(st){ try{ localStorage.setItem(KEY, JSON.stringify(st)); }catch(e){} }
    var st = load();
    var defs = document.querySelectorAll('.sidenav details[data-key]');
    defs.forEach(function(d){
      if (d.querySelector('.cur')) { d.open = true; return; }
      var k = d.getAttribute('data-key');
      if (k in st) { d.open = st[k]; }
    });
    defs.forEach(function(d){
      d.addEventListener('toggle', function(){
        var k = d.getAttribute('data-key');
        st[k] = d.open;
        save(st);
      });
    });
  })();
</script>"""

def page(prefix, active, content_html):
    return ('<!DOCTYPE html>\n<html lang="zh-CN">\n<head>\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
            '<title>bioSkills 真实试用实验室</title>\n'
            '<link rel="stylesheet" href="%sstyle.css">\n'
            '<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/styles/github.min.css">\n'
            '</head>\n<body>\n' % prefix
            + sidebar(active, prefix)
            + '<button class="menubtn" id="menubtn" aria-label="菜单">☰</button>'
            + '<div class="content"><div class="wrap">' + content_html + '</div></div>'
            + '<button class="totop" id="totop">↑</button>' + COPY_JS + '</body></html>')

# ---------- 首页 ----------
idx = '<div class="masthead"><div class="kicker">BIO / LAB</div><h1>bioSkills 真实试用实验室</h1>'
idx += '<p class="sub">基于 GPTomics/bioSkills 的可复现实践 - 严格按原方案复现 + 成分拆解 + 真实数据图</p>'
idx += '<p class="lead">这不是 AI 生成的演示稿，是<strong>真跑通</strong>的结果。bioSkills 是个教 AI coding agents 做生信任务的 skill 仓库（已 archived，结构仍有效）。这个实验室挑真实 skill，严格按它自己的方案复现，再拆开"零件"告诉你它到底教什么、坑在哪。</p></div>'
idx += '<section id="overview"><div class="sec-kicker">OVERVIEW</div><h2>bioSkills 是什么</h2>'
idx += '<div class="stats"><div class="stat"><div class="v">%d</div><div class="l">个 category</div></div><div class="stat"><div class="v">%d</div><div class="l">个 skill 总量</div></div><div class="stat"><div class="v">archived</div><div class="l">结构有效，鼓励 fork</div></div><div class="stat"><div class="v">1</div><div class="l">category 已试用</div></div></div>' % (NCAT, NTOTAL)
idx += '<div class="note"><b>定位</b>：A collection of skills that guide AI coding agents (Claude Code / Codex / Gemini / OpenCode / OpenClaw) through common bioinformatics tasks. 每个 skill 含代码范式、最佳实践、示例。</div>'
idx += '<p style="margin-top:14px;">方法论：<strong>真实数据 → 严格按 skill 原方案复现 → 拆零件（文件结构/参考脚本/核心 API/经验封装）→ 出真实数据图 → 写进站点</strong>。从左侧目录进入任意分支。</p></section>'
# ---------- 类别地图：领域着色的 Treemap（面积=skill数 / 颜色=领域）----------
TREEMAP_JS = """(function(){
  const TREE = %s;
  const tip = document.getElementById('tm-tip');
  const WRAP = document.querySelector('.tm-wrap');
  const W = 920, H = 600, PADTOP = 20, PIN = 3;
  const root = d3.hierarchy(TREE)
    .sum(function(d){ return d.value || 0; })
    .sort(function(a,b){ return b.value - a.value; });
  d3.treemap().size([W, H]).paddingOuter(2).paddingTop(PADTOP).paddingInner(PIN)(root);
  const svg = d3.select('#tm').append('svg')
    .attr('viewBox','0 0 '+W+' '+H).attr('preserveAspectRatio','xMidYMid meet')
    .style('width','100%%').style('height','auto').style('display','block');
  // 领域头（写在每域顶部留白带）
  root.children.forEach(function(d){
    svg.append('text').attr('x', d.x0 + 6).attr('y', d.y0 + 14)
      .attr('fill', d.data.color).attr('font-size', 13).attr('font-weight', 700)
      .text(d.data.name + '  ·  ' + d.value);
  });
  // 叶子
  root.leaves().forEach(function(leaf){
    const c = leaf.data, d = leaf.parent.data;
    const w = leaf.x1 - leaf.x0, h = leaf.y1 - leaf.y0;
    if (w < 1 || h < 1) return;
    const base = d3.color(d.color);
    const sibs = leaf.parent.children;
    const idx = sibs.indexOf(leaf);
    const frac = sibs.length > 1 ? idx / (sibs.length - 1) : 0;
    const fill = base.brighter(0.2 + frac * 0.9).formatHex();
    const rect = svg.append('rect')
      .attr('x', leaf.x0).attr('y', leaf.y0).attr('width', w).attr('height', h)
      .attr('fill', fill).attr('stroke', '#fff').attr('stroke-width', 1).attr('rx', 2)
      .style('cursor', c.link ? 'pointer' : 'default');
    if (c.done){ rect.attr('stroke', '#ef4444').attr('stroke-width', 2.5); }
    if (w > 46 && h > 24){
      const fs = (w > 150 && h > 40) ? 12.5 : 10.5;
      const maxC = Math.max(4, Math.floor((w - 10) / (fs * 0.6)));
      const label = c.name.length > maxC ? c.name.slice(0, maxC - 1) + '…' : c.name;
      svg.append('text').attr('x', leaf.x0 + 5).attr('y', leaf.y0 + (h > 34 ? 16 : h / 2 + 4))
        .attr('fill', '#0f172a').attr('font-size', fs).attr('font-weight', 600).text(label);
    } else if (h > 1.5 * w && h >= 22){
      const fs = 10.5, cx = leaf.x0 + w / 2, cy = leaf.y0 + h / 2;
      const maxR = Math.max(4, Math.floor((h - 10) / (fs * 1.05)));
      const label = c.name.length > maxR ? c.name.slice(0, maxR - 1) + '…' : c.name;
      svg.append('text').attr('x', cx).attr('y', cy).attr('fill', '#0f172a')
        .attr('font-size', fs).attr('font-weight', 600).attr('text-anchor', 'middle')
        .attr('transform', 'rotate(-90 ' + cx + ' ' + cy + ')').text(label);
    }
    rect.on('mousemove', function(e){
      const wr = WRAP.getBoundingClientRect();
      tip.style.display = 'block';
      let mx = e.clientX - wr.left + 14, my = e.clientY - wr.top + 14;
      if (mx + 230 > wr.width) mx = wr.width - 230;
      tip.style.left = mx + 'px'; tip.style.top = my + 'px';
      tip.innerHTML = '<b>' + c.name + '</b><br>领域：' + d.name + '<br>skill 数：×' + c.value + (c.done ? '<br><span style="color:#fca5a5">已深度试用 ✓</span>' : '');
    }).on('mouseleave', function(){ tip.style.display = 'none'; })
      .on('click', function(){ if (c.link) window.location.href = c.link; });
  });
  // 悬停 tooltip 看明细；红框(已深度试用) 点击进入 004/005 深度篇（见上方 rect click 处理）
})();""" % TREE_JSON

legend = "".join('<span class="tm-lg"><i style="background:%s"></i>%s <b>%d</b></span>' % (color, name, sum(n for c,n in dmap.get(name, [])))
                for name, color in DOMAIN_META)
map_sec = ('<section id="map"><div class="sec-kicker">CATEGORY MAP</div><h2>类别地图 · 领域全景</h2>'
 '<p class="sec-one">%d 个 category、共 %d 个 skill，归入 <b>10 个领域</b>。<b>色块面积 = 该 category 的 skill 数</b>（体量感），<b>颜色 = 所属领域</b>；<b>悬停</b>查看领域 / skill 数 / 是否深度试用。<span style="color:#ef4444;font-weight:700">红框</span> = 已深度试用，点击进入 004/005 深度篇。</p>'
 '<div class="tm-wrap"><div id="tm"></div><div class="tm-tip" id="tm-tip"></div></div>'
 '<div class="tm-legend">%s</div>'
 '<script src="https://cdnjs.cloudflare.com/ajax/libs/d3/7.9.0/d3.min.js"></script>'
 '<script>%s</script>'
 '</section>') % (NCAT, NTOTAL, legend, TREEMAP_JS)
idx += map_sec
idx += '<footer><p>所有数据均来自真实运行，代码可在对应 pipeline/ 目录复现。单篇图文版已发小红书。</p><p>bioSkills 真实试用实验室 - 持续更新中</p></footer>'
idx = page("", "index", idx)

# ---------- alignment 总览 ----------
al = '<div class="crumb"><a href="../index.html">实验室首页</a> › alignment</div>'
al += '<div class="masthead"><div class="sec-kicker">CATEGORY</div><h1>alignment</h1><p class="sub">序列比对家族 - 7 个 skill，已完成 2 个深度试用</p></div>'
al += '<div class="sklist">'
al += '<a href="pairwise-alignment.html"><span>pairwise-alignment</span><span class="tag tag-done">DONE 004</span></a>'
al += '<a href="msa-statistics.html"><span>msa-statistics</span><span class="tag tag-done">DONE 005</span></a>'
for name in ["alignment-io", "alignment-trimming", "msa-parsing", "multiple-alignment", "structural-alignment"]:
    al += '<div><span>%s</span><span class="tag tag-todo">待做</span></div>' % name
al += '</div><div class="todo" style="margin-top:18px;"><p style="margin:0;">同一套方法论持续补完：真实数据 → 严格按 skill 复现 → 成分拆解 → 出图。每完成一个 skill，上方列表自动点亮。</p></div>'
al += '<footer><p><a href="../index.html">← 返回实验室首页</a></p></footer>'
al = page("../", "alignment", al)

def skill_page(prefix, active, crumb, sec_kicker, h1, sec_one, body_html, code_blocks):
    c = '<div class="crumb"><a href="%sindex.html">实验室首页</a> › <a href="alignment/index.html">alignment</a> › %s</div>' % (prefix, crumb)
    c += '<div class="masthead"><div class="sec-kicker">%s</div><h1>%s</h1><p class="sub">%s</p></div>' % (sec_kicker, h1, sec_one)
    c += body_html
    if code_blocks:
        c += '<h2>本实验室真实复现代码</h2><p class="sec-one">以下代码来自 pipeline/ 目录，与 SKILL.md 配置严格一致，可直接复现。</p>'
    for cb in code_blocks:
        c += '<div class="codewrap"><pre><code class="language-python">%s</code></pre></div>' % cb
    c += '<footer><p><a href="alignment/index.html">← 返回 alignment 分支</a></p></footer>'
    return page(prefix, active, c)

NOTE004 = BASE + "/content/笔记/004-bioSkills真实试用-pairwise-alignment.md"
NOTE005 = BASE + "/content/笔记/005-bioSkills真实试用-msa-statistics.md"
body004, meta004, need004 = load_note(NOTE004)
p004 = skill_page("../", "pairwise", "pairwise-alignment", "DEEP DIVE 01",
                 meta004.get("标题", "pairwise-alignment"),
                 meta004.get("副标题", ""), body004, [])
body005, meta005, need005 = load_note(NOTE005)
p005 = skill_page("../", "msa", "msa-statistics", "DEEP DIVE 02",
                 meta005.get("标题", "msa-statistics"),
                 meta005.get("副标题", ""), body005, [])

# ---------- 写出 ----------
shutil.rmtree(SITE, ignore_errors=True)
os.makedirs(SITE); os.makedirs(ASSETS); os.makedirs(SITE + "/alignment")
for _src in sorted(set(need004 + need005)):
    shutil.copy(_src, ASSETS + "/" + os.path.basename(_src))
with open(SITE + "/style.css", "w", encoding="utf-8") as f: f.write(STYLE)
with open(SITE + "/index.html", "w", encoding="utf-8") as f: f.write(idx)
with open(SITE + "/alignment/index.html", "w", encoding="utf-8") as f: f.write(al)
with open(SITE + "/alignment/pairwise-alignment.html", "w", encoding="utf-8") as f: f.write(p004)
with open(SITE + "/alignment/msa-statistics.html", "w", encoding="utf-8") as f: f.write(p005)
print("站点已生成 -> " + SITE)
for root, _, files in os.walk(SITE):
    for fn in sorted(files):
        p = os.path.join(root, fn)
        print("  %s  (%d KB)" % (os.path.relpath(p, SITE), os.path.getsize(p)//1024))
