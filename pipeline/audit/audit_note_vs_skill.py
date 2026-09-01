# -*- coding: utf-8 -*-
"""核查 v2：全文扫描命令行，按 tool(+subcommand) 与 flag 双向溯源到 SKILL.md。"""
import os, re, glob, io

ROOT = 'D:/1.WorkDir/RedBook'
NOTES = os.path.join(ROOT, 'content', '笔记', 'variant-calling')
SKILLS = os.path.join(ROOT, 'content', '库', 'bioSkills', 'variant-calling')

PAIRS = [
    ('016', 'vcf-basics'), ('017', 'variant-normalization'),
    ('018', 'vcf-statistics'), ('019', 'vcf-manipulation'),
    ('020', 'filtering-best-practices'), ('021', 'variant-annotation'),
    ('022', 'gatk-variant-calling'), ('023', 'joint-calling'),
    ('024', 'structural-variant-calling'), ('025', 'consensus-sequences'),
    ('026', 'clinical-interpretation'), ('027', 'deepvariant'),
]

TOOLS = set('''bcftools samtools gatk bgzip tabix vcftools bedtools vcfanno snpeff snpsift
python python3 Rscript vt delly manta lumpy svaba whatshap deepvariant docker singularity
apptainer vep annovar freebayes varscan mutect2 strelka lofreq conda pip curl wget awk sed
grep cut sort uniq head tail wc cat less echo cd mkdir cp mv rm ls find zcat gunzip make
minimap2 bwa picard jupyter quarto R'''.split())

FLAG = re.compile(r'(?<![\w-])(--[A-Za-z][\w-]*)(?![\w-])')

def scan_cmds(text):
    """全文扫描：以工具名开头的行（允许前导 $ / 反引号 / 空格 / 列表符）。"""
    out = []
    for raw in text.split('\n'):
        s = raw.strip()
        if not s or s.startswith('#') or s.startswith('```'):
            continue
        s = s.lstrip('-*>').strip()
        s = s.strip('`').strip()
        if s.startswith('$ '):
            s = s[2:].strip()
        parts = s.split()
        if not parts:
            continue
        base = os.path.basename(parts[0]).lower()
        if base in TOOLS:
            out.append(s)
    return out

def sub_of(cmd):
    parts = cmd.split()
    base = os.path.basename(parts[0]).lower()
    sub = ''
    if base == 'gatk' and len(parts) > 1 and not parts[1].startswith('-'):
        sub = parts[1]
    elif base in ('bcftools', 'samtools', 'bedtools', 'vcftools') and len(parts) > 1 and not parts[1].startswith('-'):
        sub = parts[1]
    return (base + ' ' + sub).strip()

def main():
    lines = []
    for num, topic in PAIRS:
        ncand = glob.glob(os.path.join(NOTES, num + '-*' + topic + '.md'))
        skillp = os.path.join(SKILLS, topic, 'SKILL.md')
        if not ncand:
            lines.append('## %s %s NOTE_NOT_FOUND' % (num, topic)); continue
        if not os.path.exists(skillp):
            lines.append('## %s %s SKILL_NOT_FOUND' % (num, topic)); continue
        note = io.open(ncand[0], encoding='utf-8').read()
        skill = io.open(skillp, encoding='utf-8').read()
        skl = skill.lower()
        cmds = scan_cmds(note)
        sub_miss, flag_miss = [], []
        for c in cmds:
            sc = sub_of(c)
            if sc and sc not in skl:
                sub_miss.append((sc, c))
            fl = FLAG.findall(c)
            miss = [f for f in fl if f.lower() not in skl]
            if miss:
                flag_miss.append((c, sorted(set(miss))))
        lines.append('## %s %s' % (num, topic))
        lines.append('  cmds=%d  subcmd_untraceable=%d  flag_untraceable=%d' % (len(cmds), len(sub_miss), len(flag_miss)))
        for sc, c in sub_miss[:15]:
            lines.append('  [SUB?] %-28s | %s' % (sc, c[:110]))
        for c, miss in flag_miss[:15]:
            lines.append('  [FLAG?] %-58s MISSING=%s' % (c[:58], ','.join(miss)))
        lines.append('')
    out = os.path.join(ROOT, 'pipeline', 'audit', 'audit_raw.txt')
    io.open(out, 'w', encoding='utf-8').write('\n'.join(lines))
    print('\n'.join(lines))

main()
