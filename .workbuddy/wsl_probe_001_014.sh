set -e
echo "=== WSL distro ==="; wsl.exe -l -v 2>/dev/null || echo "no wsl"
echo "=== existing tools in WSL ==="
for t in conda mamba mafft muscle clustalo bowtie2 bwa hisat2 STAR clipkit trimal bmge foldseek tmalign foldmason; do
  p=$(command -v $t 2>/dev/null || echo "-")
  echo "$t: $p"
done
echo "=== conda channels reachable? ==="
if command -v conda >/dev/null 2>&1; then
  conda search -c bioconda mafft --override-channels 2>&1 | tail -3 || echo "conda search failed"
else
  echo "conda not installed in WSL"
fi
