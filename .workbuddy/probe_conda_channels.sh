#!/usr/bin/env bash
# 探测 conda 通道可达性（决定 022/023 能否装包）
set -u
echo "===== CONDA CHANNEL REACHABILITY ====="
for u in \
  "https://conda.anaconda.org/bioconda/linux-64/repodata.json" \
  "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/linux-64/repodata.json" \
  "https://repo.anaconda.com/pkgs/main/linux-64/repodata.json" \
  "https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/linux-64/repodata.json" ; do
  code=$(timeout 20 curl -s -o /dev/null -w "%{http_code}" "$u" 2>/dev/null || echo FAIL)
  printf "%-72s %s\n" "$(echo "$u" | sed 's|https://||')" "$code"
done

echo ""
echo "===== CONDA CONFIG (已有通道) ====="
conda config --show channels 2>/dev/null | head -10
echo "--- .condarc ---"
cat ~/.condarc 2>/dev/null | head -20 || echo "(无 .condarc)"

echo ""
echo "===== 能否解析 bcftools / gatk4 (dry-run, 不安装) ====="
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio
timeout 180 conda install -y --dry-run -c bioconda bcftools 2>&1 | tail -6
