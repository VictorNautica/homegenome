#!/bin/bash

# start with a one-person dir
# of per-chr imputation results

dir=./imputed_vcfs_HRC_0730
id=1
bcftools concat $dir/{1..22}.vcf.gz -Oz > $dir/all.vcf.gz
tabix $dir/all.vcf.gz
bcftools view $dir/all.vcf.gz -R all_effects_snps.txt -Ov | \
	awk '$0!~/^#/{split($10, gta, ":"); gt=gta[1]; ds=gta[3];
	print $1, $2, $3, $4, $5, gt, ds}' | gzip -c > allchr_${id}.txt.gz
