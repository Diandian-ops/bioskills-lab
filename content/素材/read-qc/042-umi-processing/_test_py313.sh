#!/usr/bin/env bash
# 042 附属实验：实测 umi_tools 1.1.6 在 python 3.13 下的打包/安装表现（属性表依据）
set -uo pipefail
OUT="/mnt/d/1.WorkDir/RedBook/content/素材/read-qc/042-umi-processing"
source /opt/miniconda3/etc/profile.d/conda.sh
conda create -y -n umi-py313 python=3.13 > "$OUT/py313_create.log" 2>&1
conda activate umi-py313
pip install umi-tools > "$OUT/py313_pip.log" 2>&1
RC=$?
echo "pip install exit=$RC"
tail -30 "$OUT/py313_pip.log"
# 清理临时环境
conda deactivate
conda env remove -y -n umi-py313 > /dev/null 2>&1
echo "cleanup done"
