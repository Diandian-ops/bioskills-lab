#!/bin/bash
source /opt/miniconda3/etc/profile.d/conda.sh
conda init bash >/dev/null 2>&1 || true
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
channel_alias: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
proxy_servers:
  http: http://127.0.0.1:7897
  https: http://127.0.0.1:7897
remote_max_retries: 5
remote_connect_timeout_secs: 30
RC
echo "[cfg] condarc written (TUNA mirror + proxy)"
if conda env list | grep -qE "^bio\s+/opt/miniconda3/envs/bio"; then
  echo "[skip] bio env already exists"
else
  for i in 1 2 3; do
    echo "=== attempt $i ==="
    conda create -y -n bio -c conda-forge -c bioconda \
      python=3.11 biopython numpy matplotlib \
      mafft muscle clustalo clipkit trimal bmge \
      foldseek tmalign foldmason \
      bowtie2 bwa hisat2 star samtools dwgsim bedtools 2>&1 | tail -6
    if conda env list | grep -qE "^bio\s+/opt/miniconda3/envs/bio"; then echo "ENV_OK"; break; fi
    sleep 5
  done
fi
echo "[verify]"
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio
for t in mafft muscle clustalo clipkit trimal bmge foldseek tmalign foldmason bowtie2 bwa hisat2 STAR samtools dwgsim bedtools python; do
  printf "%-12s %s\n" "$t" "$(command -v $t || echo MISSING)"
done
echo "BIO_ENV_DONE"
