echo "=== proxy env ==="
echo "https_proxy=$https_proxy"; echo "http_proxy=$http_proxy"
echo "=== connectivity (curl -sI, 12s timeout each) ==="
for u in \
  "https://conda.anaconda.org/bioconda/" \
  "https://conda.anaconda.org/conda-forge/" \
  "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/" \
  "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/" \
  "https://repo.anaconda.com/pkgs/main/" ; do
  code=$(curl -sI --max-time 12 -o /dev/null -w "%{http_code}" "$u" 2>/dev/null || echo "ERR")
  echo "$code  $u"
done
echo "=== conda present? ==="
source /opt/miniconda3/etc/profile.d/conda.sh 2>/dev/null && conda --version && conda info --envs
