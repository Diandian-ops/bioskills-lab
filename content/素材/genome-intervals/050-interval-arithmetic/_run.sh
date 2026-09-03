#!/usr/bin/env bash
# 050 interval-arithmetic: real bedtools run following SKILL.md scope.
# Everything the skill mentions that is testable here gets exercised:
#   intersect output modes / -f -F -r / subtract(-A) / merge(-d,-c,-o) /
#   complement / cluster / map / groupby / multiinter / unionbedg /
#   jaccard / fisher / -sorted contract / chrom-naming footgun / -split.
set -uo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio

echo "== [0] environment versions =="
{
  echo "bedtools: $(bedtools --version)"
  echo "python3: $(python3 --version 2>&1)"
  conda list -n bio 2>/dev/null | grep -E "^(bedtools|pybedtools) "
} | tee 00_env.log

echo "== [1] make inputs =="
python3 make_inputs.py 2>&1 | tee 01_make_inputs.log

echo "== [2] sorted working copies =="
bedtools sort -g genome.txt -i peaks.bed     > peaks_sorted.bed
bedtools sort -g genome.txt -i genes.bed     > genes_sorted.bed
bedtools sort -g genome.txt -i exons.bed     > exons_sorted.bed
bedtools sort -g genome.txt -i scores.bedgraph > scores_sorted.bedgraph
wc -l peaks.bed peaks_sorted.bed genes_sorted.bed exons_sorted.bed | tee 02_sort.log

echo "== [3] intersect output modes =="
{
  bedtools intersect -a peaks_sorted.bed -b genes_sorted.bed -u     > i_u.bed
  bedtools intersect -a peaks_sorted.bed -b genes_sorted.bed -v     > i_v.bed
  bedtools intersect -a peaks_sorted.bed -b genes_sorted.bed -c     > i_c.bed
  bedtools intersect -a peaks_sorted.bed -b genes_sorted.bed -wa -wb > i_wawb.bed
  bedtools intersect -a peaks_sorted.bed -b genes_sorted.bed -loj    > i_loj.bed
  bedtools intersect -a peaks_sorted.bed -b genes_sorted.bed -wo     > i_wo.bed
  bedtools intersect -a peaks_sorted.bed -b genes_sorted.bed -wao    > i_wao.bed
  for f in i_u i_v i_wawb i_loj i_wo i_wao i_c; do
    echo "$f lines: $(wc -l < $f.bed)"
  done
} 2>&1 | tee 03_intersect.log

echo "== [4] fractional overlap -f / -F / -r, role asymmetry =="
{
  bedtools intersect -a peaks_sorted.bed -b genes_sorted.bed -f 0.5 -u > i_f05.bed
  bedtools intersect -a peaks_sorted.bed -b genes_sorted.bed -f 0.5 -r -u > i_f05r.bed
  bedtools intersect -a genes_sorted.bed -b peaks_sorted.bed -f 0.5 -u > i_f05_swap.bed
  echo "f05(A=peaks)=$(wc -l < i_f05.bed)  f05r=$(wc -l < i_f05r.bed)  f05_swapped(A=genes)=$(wc -l < i_f05_swap.bed)"
} 2>&1 | tee 04_frac.log

echo "== [5] subtract / merge / complement / cluster =="
{
  bedtools subtract -a peaks_sorted.bed -b blacklist.bed > sub_plain.bed
  bedtools subtract -a peaks_sorted.bed -b blacklist.bed -A > sub_A.bed
  # footgun probe: merge on the UNSORTED file (v2.31.1 refuses, exit 1)
  set +e
  bedtools merge -i peaks.bed > merge_unsorted.bed 2> merge_unsorted.err
  rc=$?; echo "$rc" > merge_unsorted.exit
  set -u
  echo "merge_unsorted exit=$rc  err: $(head -1 merge_unsorted.err)"
  bedtools merge -i peaks_sorted.bed           > merge_d0.bed
  bedtools merge -d 1 -i peaks_sorted.bed      > merge_d1.bed
  bedtools merge -d 100 -i peaks_sorted.bed    > merge_d100.bed
  bedtools merge -c 4,5 -o distinct,sum -i peaks_sorted.bed > merge_co.bed
  bedtools complement -i peaks_sorted.bed -g genome.txt > comp.bed
  bedtools sort -i peaks.bed | bedtools cluster -d 0 > cluster.bed
  echo "sub_plain ivs=$(wc -l < sub_plain.bed)  sub_A remaining=$(wc -l < sub_A.bed)"
  echo "merge blocks: unsorted=$(wc -l < merge_unsorted.bed) d0=$(wc -l < merge_d0.bed) d1=$(wc -l < merge_d1.bed) d100=$(wc -l < merge_d100.bed)"
  echo "complement ivs=$(wc -l < comp.bed)"
  echo "cluster ids max: $(awk '{if($6>m)m=$6}END{print m}' cluster.bed)"
  echo "--- merge_co head ---"; head -3 merge_co.bed
} 2>&1 | tee 05_setops.log

echo "== [6] map / groupby =="
{
  bedtools map -a genes_sorted.bed -b scores_sorted.bedgraph -c 4 -o mean > map_mean.bed
  bedtools intersect -a genes_sorted.bed -b peaks_sorted.bed -wo \
    | bedtools groupby -g 1,2,3,4 -c 12 -o sum > groupby.tsv
  echo "--- map head ---";  head -3 map_mean.bed
  echo "--- groupby head ---"; head -3 groupby.tsv
} 2>&1 | tee 06_map_groupby.log

echo "== [7] multiinter / unionbedg =="
{
  bedtools multiinter -header -names rep1 rep2 rep3 -i rep1.bed rep2.bed rep3.bed > multiinter.tsv
  bedtools unionbedg -header -names rep1 rep2 rep3 -i ubg1.bedgraph ubg2.bedgraph ubg3.bedgraph > unionbedg.tsv
  echo "multiinter all-3 rows: $(awk 'NR>1 && $4==3' multiinter.tsv | wc -l)"
  echo "multiinter all-3 bp:   $(awk 'NR>1 && $4==3{s+=$3-$2}END{print s}' multiinter.tsv)"
  echo "unionbedg rows: $(wc -l < unionbedg.tsv)"
  head -3 unionbedg.tsv
} 2>&1 | tee 07_multi_union.log

echo "== [8] jaccard / fisher (on merged inputs, mechanics only) =="
{
  bedtools merge -i genes_sorted.bed > genes_d0.bed
  bedtools jaccard -a merge_d0.bed -b genes_d0.bed -g genome.txt > jac.txt
  bedtools fisher  -a merge_d0.bed -b genes_d0.bed -g genome.txt > fisher.txt
  cat jac.txt
  echo "--- fisher ---"; cat fisher.txt
} 2>&1 | tee 08_jaccard_fisher.log

echo "== [9] -split: BED12 envelope vs blocks =="
{
  bedtools intersect -a spliced.bed12 -b exons_sorted.bed -u        > spl_env_u.bed
  bedtools intersect -a spliced.bed12 -b exons_sorted.bed -u -split > spl_blk_u.bed
  bedtools intersect -a spliced.bed12 -b exons_sorted.bed -wo       > spl_env.wo
  bedtools intersect -a spliced.bed12 -b exons_sorted.bed -wo -split > spl_blk.wo
  echo "transcripts hitting exons: envelope=$(wc -l < spl_env_u.bed) split=$(wc -l < spl_blk_u.bed)"
  echo "overlap bp: envelope=$(awk '{s+=$NF}END{print s}' spl_env.wo) split=$(awk '{s+=$NF}END{print s}' spl_blk.wo)"
} 2>&1 | tee 09_split.log

echo "== [10] -sorted contract: error cases + correct sweep =="
{
  bedtools intersect -a peaks.bed -b genes_sorted.bed -sorted -u > /dev/null 2> err_unsorted.log
  rc=$?; echo "$rc" > err_unsorted.exit; echo "exit_unsorted=$rc"
  cat err_unsorted.log
  # reversed chromosome order (chr2 block before chr1) vs genome.txt order
  awk '$1=="chr2"' peaks_sorted.bed > po1.tmp
  awk '$1=="chr1"' peaks_sorted.bed > po2.tmp
  cat po1.tmp po2.tmp > peaks_revorder.bed
  awk '$1=="chr2"' genes_sorted.bed > go1.tmp
  awk '$1=="chr1"' genes_sorted.bed > go2.tmp
  cat go1.tmp go2.tmp > genes_revorder.bed
  bedtools intersect -a peaks_revorder.bed -b genes_revorder.bed -sorted -g genome.txt -u > /dev/null 2> err_revorder.log
  rc=$?; echo "$rc" > err_revorder.exit; echo "exit_revorder=$rc"
  cat err_revorder.log
  # correct low-memory sweep, pinned by -g
  bedtools intersect -a peaks_sorted.bed -b genes_sorted.bed -sorted -g genome.txt -u > i_u_sorted_sweep.bed
  bedtools intersect -a peaks_sorted.bed -b genes_sorted.bed -u > i_u_sorted.bed
  echo "sweep_u=$(wc -l < i_u_sorted_sweep.bed) inmem_u=$(wc -l < i_u_sorted.bed)"
  rm -f po1.tmp po2.tmp go1.tmp go2.tmp
} 2>&1 | tee 10_sorted_contract.log

echo "== [11] chrom-naming footgun (chr1 vs 1) =="
{
  sed 's/^chr//' peaks_sorted.bed > peaks_noprefix.bed
  set +e
  bedtools intersect -a peaks_noprefix.bed -b genes_sorted.bed -u > mism_u.bed 2> mism_err.log
  echo "exit=$?"
  echo "mismatch -u lines: $(wc -l < mism_u.bed)  stderr bytes: $(wc -c < mism_err.log)"
  set -u
} 2>&1 | tee 11_mismatch.log

echo "== [12] pybedtools parity =="
python3 pybedtools_check.py 2>&1 | tee 11_pybedtools.log

echo "== [13] perf: in-memory tree vs -sorted sweep (300k x 300k) =="
{
  python3 perf_gen.py
  t0=$(date +%s%N)
  bedtools intersect -a big_a.bed -b big_b.bed -c > /dev/null
  t1=$(date +%s%N)
  awk -v a="$t0" -v b="$t1" 'BEGIN{printf "inmem elapsed_s: %.2f\n", (b-a)/1e9}'
  t0=$(date +%s%N)
  bedtools intersect -a big_a.bed -b big_b.bed -sorted -g big_genome.txt -c > /dev/null
  t1=$(date +%s%N)
  awk -v a="$t0" -v b="$t1" 'BEGIN{printf "sorted elapsed_s: %.2f\n", (b-a)/1e9}'
  if [ -x /usr/bin/time ]; then
    /usr/bin/time -v bedtools intersect -a big_a.bed -b big_b.bed -c > /dev/null 2> perf_inmem.time
    /usr/bin/time -v bedtools intersect -a big_a.bed -b big_b.bed -sorted -g big_genome.txt -c > /dev/null 2> perf_sorted.time
    echo "inmem maxrss_kb: $(grep 'Maximum resident' perf_inmem.time | grep -o '[0-9]*')"
    echo "sorted maxrss_kb: $(grep 'Maximum resident' perf_sorted.time | grep -o '[0-9]*')"
  else
    echo "/usr/bin/time not available: max RSS not measured"
  fi
} 2>&1 | tee perf_summary.txt | tee 12_perf.log

echo "== [14] reconcile measured vs expected =="
python3 parse_results.py 2>&1 | tee 13_parse.log

echo "== [15] assemble repro_transcript.txt =="
{
  echo "# 050 interval-arithmetic repro transcript"
  echo "# generated: $(date -Is)"
  for f in 00_env.log 01_make_inputs.log 02_sort.log 03_intersect.log \
           04_frac.log 05_setops.log 06_map_groupby.log 07_multi_union.log \
           08_jaccard_fisher.log 09_split.log 10_sorted_contract.log \
           11_mismatch.log 11_pybedtools.log 12_perf.log perf_summary.txt \
           13_parse.log results.txt jac.txt fisher.txt err_unsorted.log \
           err_revorder.log merge_unsorted.err pybedtools_numbers.txt; do
    [ -f "$f" ] && { echo ""; echo "===== $f ====="; cat "$f"; }
  done
} > repro_transcript.txt
echo "DONE"
