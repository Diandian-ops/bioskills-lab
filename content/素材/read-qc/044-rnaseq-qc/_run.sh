#!/bin/bash
# 044 rnaseq-qc 真跑主流程（WSL Ubuntu，bio-qc 环境）
# 口径：严格按 SKILL.md —— salmon index -> salmon quant -l A 自动判定链向，
#       另跑 -l IU（无链）与 -l SF（错链）对照验证 SKILL.md 文库表；
#       RSeQC/Picard/Qualimap/RNA-SeQC 本机缺失，对应环节不跑（诚实标注）。
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio-qc

echo "== [1/6] make_inputs.py =="
python3 make_inputs.py 2>&1 | tee logs/_make_inputs.log

echo "== [2/6] salmon index =="
salmon index -t transcripts.fa -i salmon_index 2>&1 | tee logs/_index.log

echo "== [3/6] salmon quant -l A (auto-detect, per SKILL.md CLI) =="
salmon quant -i salmon_index -l A -r reads_se.fq.gz -o quant_A 2>&1 | tee logs/_quant_A.log

echo "== [4/6] salmon quant -l IU (unstranded control) =="
salmon quant -i salmon_index -l IU -r reads_se.fq.gz -o quant_IU 2>&1 | tee logs/_quant_IU.log

echo "== [5/6] salmon quant -l SF (wrong-strand control; non-zero exit is an expected outcome) =="
salmon quant -i salmon_index -l SF -r reads_se.fq.gz -o quant_SF 2>&1 | tee logs/_quant_SF.log \
    || echo "[expected] -l SF exited non-zero on antisense reads (see log)"

echo "== [6/6] parse =="
python3 parse_quant.py 2>&1 | tee logs/_parse.log
# 出图（make_figs.py）与本汇总（_make_repro.sh）在受管 venv / 汇总阶段执行：
# bio-qc 环境无 matplotlib，图由 Windows 侧受管 venv python 生成（logs/_figs.log）。

# ---------- 汇总 repro_transcript.txt ----------
bash _make_repro.sh

# ---------- 汇总 repro_transcript.txt ----------
{
    echo "### 044 rnaseq-qc repro transcript（WSL Ubuntu + bio-qc，salmon $(salmon --version 2>&1)）"
    echo "生成时间：$(date -Is)"
    echo
    echo "== 环境核查（_env_check.log）=="
    cat _env_check.log
    echo
    echo "== step1 make_inputs =="
    cat logs/_make_inputs.log
    echo
    echo "== step2 salmon index（日志尾部）=="
    tail -n 12 logs/_index.log
    echo
    for tag in A IU SF; do
        echo "== step salmon quant -l $tag（完整日志）=="
        cat logs/_quant_$tag.log
        echo
    done
    for d in quant_A quant_IU quant_SF; do
        if [ -f "$d/lib_format_counts.json" ]; then
            echo "== $d/lib_format_counts.json =="
            cat "$d/lib_format_counts.json"
            echo
        fi
    done
    echo "== step6 parse_quant =="
    cat logs/_parse.log
    echo
    echo "== parsed_results.tsv =="
    cat parsed_results.tsv
    echo
    echo "== make_figs（Windows 侧受管 venv 生成，logs/_figs.log）=="
    if [ -f logs/_figs.log ]; then cat logs/_figs.log; else echo "(待生成)"; fi
} > repro_transcript.txt
echo "repro_transcript.txt written"
