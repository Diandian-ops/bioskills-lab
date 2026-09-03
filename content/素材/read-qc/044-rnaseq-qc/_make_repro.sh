#!/bin/bash
# 044 rnaseq-qc：仅重建 repro_transcript.txt（在 make_figs 落盘 logs/_figs.log 后调用）
set -euo pipefail
cd "$(dirname "$0")"
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio-qc
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
