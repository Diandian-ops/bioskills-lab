#!/usr/bin/env bash
# 022 诊断：HaplotypeCaller 标准模式 0 变异的根因验证（全部输出落盘，不喷终端）。
set -u
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bio

MAT=/mnt/d/1.WorkDir/RedBook/content/素材/variant-calling/022-gatk-variant-calling
LOG=$MAT/_diag.log
: > "$LOG"

echo "=== [1] 修正统计法重算 mpileup 非参考列（同时排除 . 与 ,）===" | tee -a "$LOG"
awk '{
  b=$5; gsub(/\^./,"",b); gsub(/\$/,"",b);
  while (match(b,/[+-][0-9]+/)) {
    n=substr(b,RSTART+1,RLENGTH-1)+0;
    b=substr(b,1,RSTART-1) substr(b,RSTART+RLENGTH+1+n);
  }
  if (b ~ /[ACGTNacgtn]/) c++;
} END{print "total pileup cols (covered): " NR; print "cols with non-ref base(s), corrected: " c+0}' \
  "$MAT/_mpileup.txt" | tee -a "$LOG"

echo "=== [2] bcftools mpileup | bcftools call -mv（跨工具对照）===" | tee -a "$LOG"
bcftools mpileup -f "$MAT/reference.fa" -a FORMAT/DP,FORMAT/AD "$MAT/aligned_rg.bam" 2>/dev/null \
  | bcftools call -mv -Oz -o "$MAT/_bcftools_calls.vcf.gz" 2>> "$LOG"
bcftools index -tf "$MAT/_bcftools_calls.vcf.gz" 2>> "$LOG"
bcftools stats "$MAT/_bcftools_calls.vcf.gz" > "$MAT/_diag_bcftools_stats.txt" 2>&1
grep "^SN" "$MAT/_diag_bcftools_stats.txt" | head -12 | tee -a "$LOG"
echo "--- called sites: CHROM POS REF ALT QUAL DP AD ---" | tee -a "$LOG"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%DP\t%AD]\n' "$MAT/_bcftools_calls.vcf.gz" | tee -a "$LOG"

echo "=== [3] 候选真实突变扫描：每列最高 ALT 支持度 TOP10 ===" | tee -a "$LOG"
MAT=$MAT python3 - <<'PY' | tee -a "$LOG"
import os, re, collections
indel = re.compile(r'[+-](\d+)')
cols = []
with open(os.path.join(os.environ['MAT'], '_mpileup.txt')) as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        chrom, pos, ref, dp, bases = p[0], int(p[1]), p[2], int(p[3]), p[4]
        b = re.sub(r'\^.', '', bases).replace('$', '')
        while True:
            m = indel.search(b)
            if not m: break
            n = int(m.group(1)); s = m.start()
            b = b[:s] + b[m.end():][n:]
        cnt = collections.Counter(ch.lower() for ch in b if ch in 'ACGTacgtnN')
        alt_n = sum(cnt.values()); alt_frac = alt_n / dp if dp else 0.0
        alt_top = cnt.most_common(1)[0][0] if cnt else '-'
        cols.append((alt_n, alt_frac, chrom, pos, ref, dp, alt_top))
cols.sort(reverse=True)
print('chrom\tpos\tref\tDP\ttop_alt\talt_reads\talt_frac')
for alt_n, alt_frac, chrom, pos, ref, dp, alt_top in cols[:10]:
    print(f'{chrom}\t{pos}\t{ref}\t{dp}\t{alt_top}\t{alt_n}\t{alt_frac:.3f}')
n_strong = sum(1 for alt_n, alt_frac, *_ in cols if alt_n >= 10 and alt_frac >= 0.5)
print(f'cols with alt_reads>=10 AND alt_frac>=0.5 (real-mutation-like): {n_strong}')
PY

echo "=== [4] BAM error rate / 平均碱基质量（samtools stats 独立复核）===" | tee -a "$LOG"
samtools stats "$MAT/aligned_rg.bam" > "$MAT/_diag_bamstats.txt" 2>&1
grep -E "^SN\s+(error rate|average quality|bases mapped \(cigar\)|reads mapped|average length|reads unmapped)" \
  "$MAT/_diag_bamstats.txt" | tee -a "$LOG"

echo "=== [5] GVCF 覆盖度分布（HC 是否真的看到了 reads）===" | tee -a "$LOG"
bcftools query -f '%POS\t%INFO/END\t[%DP]\n' "$MAT/raw.g.vcf.gz" | \
  awk '{d=$3; if(d==0) z++; else {s+=d; c++}} END{print "GVCF blocks total: " NR; print "blocks with DP=0: " z+0; print "blocks with DP>0: " c+0; if(c>0) printf "mean DP of covered blocks: %.1f\n", s/c}' | tee -a "$LOG"

echo "=== [6] wgsim 读名内嵌错配数（含 e=0.02 测序错误）===" | tee -a "$LOG"
samtools view "$MAT/aligned_rg.bam" | awk -F'\t' '{split($1,a,"_"); n=(a[4]+0)+(a[5]+0); s+=n; c++; if(n==0)z0++} END{print "reads: " c; printf "mean header mismatch count per read: %.2f\n", s/c; print "reads with 0 header mismatch: " z0+0}' | tee -a "$LOG"

echo "=== done $(date -u) ===" | tee -a "$LOG"
