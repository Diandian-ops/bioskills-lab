#!/bin/bash
source /opt/miniconda3/etc/profile.d/conda.sh
conda init bash >/dev/null 2>&1 || true
export PIP_OVERRIDE=1
# 用 TUNA 镜像加速（已验证可达），并强制绕过 defaults 的 Anaconda ToS 卡死
cat > /opt/miniconda3/.condarc <<'RC'
channels:
  - conda-forge
  - bioconda
show_channel_urls: true
default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
custom_channels:
  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  pytorch: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
RC
echo "[cfg] condarc written"
# 已存在则跳过创建
if conda env list | grep -q "/opt/miniconda3/envs/bio$"; then
  echo "[skip] bio env already exists"
else
  echo "[1/2] creating bio env (mirror + override-channels + retries)..."
  for i in 1 2 3; do
    echo "=== attempt $i ==="
    conda create -y -n bio --override-channels -c conda-forge -c bioconda \
      python=3.11 biopython numpy matplotlib \
      mafft muscle clustalo clipkit trimal bmge \
      foldseek tmalign foldmason \
      bowtie2 bwa hisat2 star samtools dwgsim bedtools 2>&1 | tail -8
    if conda env list | grep -q "/opt/miniconda3/envs/bio$"; then echo "ENV_OK"; break; fi
    sleep 5
  done
fi
echo "[2/2] verify tools in bio env ..."
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio
for t in mafft muscle clustalo clipkit trimal bmge foldseek tmalign foldmason bowtie2 bwa hisat2 STAR samtools dwgsim bedtools python; do
  printf "%-12s %s\n" "$t" "$(command -v $t || echo MISSING)"
done
echo "BIO_ENV_DONE"
