#!/usr/bin/env bash
# 022 / 023 / 027 补跑可行性【只读】探测 — 在 WSL bio 环境执行
# 用法: wsl.exe -d Ubuntu -u root -- bash /mnt/d/1.WorkDir/RedBook/.workbuddy/probe_022_027.sh
set -u
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio

echo "===== 1. TOOL AVAILABILITY ====="
for t in gatk bcftools samtools bgzip tabix java python3; do
  printf "%-12s" "$t"
  if command -v "$t" >/dev/null 2>&1; then
    printf "FOUND   "
    case "$t" in
      bcftools) bcftools --version 2>&1 | head -1 ;;
      samtools) samtools --version 2>&1 | head -1 ;;
      bgzip|tabix) echo "(htslib bundled)" ;;
      gatk) timeout 60 gatk --version 2>&1 | grep -iE "version|Genome Analysis" | head -2 ;;
      java) java -version 2>&1 | head -1 ;;
      python3) python3 -c "import sys; print('python', sys.version.split()[0])" ;;
    esac
  else
    echo "MISSING"
  fi
done

echo ""
echo "===== 2. DEEPVARIANT ====="
if command -v dv_make_examples >/dev/null 2>&1; then
  echo "dv_make_examples FOUND: $(command -v dv_make_examples)"
else
  echo "dv_make_examples MISSING"
fi
python3 -c "import deepvariant; print('deepvariant python module FOUND')" 2>&1 | head -1
command -v docker >/dev/null 2>&1 && echo "docker FOUND" || echo "docker MISSING (deepvariant 官方唯一分发途径)"

echo ""
echo "===== 3. NETWORK (能否下载 gatk jar / deepvariant) ====="
timeout 15 curl -s -o /dev/null -w "github      = %{http_code}\n" https://github.com || echo "github      = FAIL"
timeout 15 curl -s -o /dev/null -w "broadinstitute = %{http_code}\n" https://software.broadinstitute.org || echo "broadinstitute = FAIL"
timeout 15 curl -s -o /dev/null -w "pypi        = %{http_code}\n" https://pypi.org || echo "pypi        = FAIL"

echo ""
echo "===== 4. 现有 1000G chr22 数据 (015 素材目录, 022/023/027 可复用) ====="
D=/mnt/d/1.WorkDir/RedBook/content/素材/variant-calling
ls -la "$D/015-variant-calling/" 2>/dev/null | head -15
echo "--- 024/025/026 供参考 ---"
ls -la "$D/024-sv-calling/" 2>/dev/null | head -8
