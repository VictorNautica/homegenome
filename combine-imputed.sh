#!/bin/bash

# USAGE: bash combine-imputed.sh 
# Put this script in a directory with Sanger-imputed VCFs, i.e. 1.vcf.gz-22.vcf.gz.
# Only autosomes are used currently.
# NOTE: each person's VCFs must be in a separate directory!
#
# Run this script (once for each person).
#
# It will output a merged file in the same directory, marked with the ID argument.
# You can delete the merged file after you run all the conversion scripts.

set -e

dir=./

if [ ! -f "$dir/1.vcf.gz" ];
then
	echo "File $dir/1.vcf.gz not found - check if you are in the correct folder."
	exit 1
fi

echo -n "Enter desired user ID (a number or several letters) + ENTER: "
read id
if [ -f "/srv/shiny-server/homegenome/grs/grs_${id}.tsv" ];
then
	echo "File grs_${id}.tsv already exists, choose a different user ID"
	exit 1
fi	
if [ -f "/srv/shiny-server/homegenome/grs/cv_${id}.tsv" ];
then
	echo "File cv_${id}.tsv already exists, choose a different user ID"
	exit 1
fi

## All checks passed, continue:

bcftools concat $dir/{1..22}.vcf.gz -Oz > $dir/all.vcf.gz
tabix $dir/all.vcf.gz
bcftools view $dir/all.vcf.gz \
	-R /srv/shiny-server/homegenome/all_effects_snps.txt -Ov | \
	awk '$0!~/^#/{split($10, gta, ":"); gt=gta[1]; ds=gta[3];
	print $1, $2, $3, $4, $5, gt, ds}' | gzip -c > allchr_${id}.txt.gz

## Proceed to merge with GRS, clinvar, and encrypt
