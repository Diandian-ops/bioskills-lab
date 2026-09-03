#!/bin/bash
# 044 rnaseq-qc 环境核查：确认 salmon/dwgsim 及 RSeQC 系工具是否可用
source /opt/miniconda3/etc/profile.d/conda.sh
echo "=== bio-qc env ==="
conda activate bio-qc
echo "salmon: $(command -v salmon) $(salmon --version 2>&1)"
echo "fastp:  $(command -v fastp) $(fastp --version 2>&1)"
echo "samtools: $(command -v samtools) $(samtools --version 2>&1 | head -1)"
for t in infer_experiment.py geneBody_coverage.py tin.py read_distribution.py read_duplication.py junction_saturation.py picard qualimap rnaseqc how_are_we_stranded_here.py; do
  printf "%-24s %s\n" "$t:" "$(command -v $t || echo MISSING)"
done
echo ""
echo "=== bio env ==="
conda deactivate
conda activate bio
echo "dwgsim: $(command -v dwgsim) $(dwgsim --version 2>&1 | head -1)"
for t in infer_experiment.py salmon picard qualimap; do
  printf "%-24s %s\n" "$t:" "$(command -v $t || echo MISSING)"
done
