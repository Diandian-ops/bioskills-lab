#!/bin/bash
IN=/mnt/d/1.WorkDir/RedBook/content/素材/variant-calling/023-joint-calling/chr22_slice.vcf.gz
echo "all sites for HG00096:"; bcftools view -s HG00096 "$IN" | wc -l
echo "two-step non-ref count HG00096:"; bcftools view -s HG00096 "$IN" -Ou | bcftools view -e 'GT="ref"' -Ou | bcftools view -H | wc -l
echo -n "positions HG00096: "; bcftools view -s HG00096 "$IN" -Ou | bcftools view -e 'GT="ref"' -Ou | bcftools query -f '%POS\n' | tr '\n' ' '; echo
echo -n "positions HG01280: "; bcftools view -s HG01280 "$IN" -Ou | bcftools view -e 'GT="ref"' -Ou | bcftools query -f '%POS\n' | tr '\n' ' '; echo
