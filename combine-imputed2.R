options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(sodium)

## This script processes individual full-genome files into GRS table.
## USAGE: Rscript combine-imputed2.R
## To be run separately for each person, in a directory with that person's VCFs.

# read CV and UKBB data
dirsrv = "/srv/shiny-server/homegenome/"
cv = read.table(paste0(dirsrv, "clinvar_filtered.tsv"),
				h=T, sep="\t", quote="", comment.char = "")
ukbbSum = read.table(paste0(dirsrv, "ukbb_summarizedgrs.tsv"),
				h=T, sep="\t", quote="")
ukbb = read.table(paste0(dirsrv, "all_effects.tsv"), h=F, sep="\t", quote="")
colnames(ukbb) = c("CHR", "POS", "REF", "ALT", "rs",
				   "n", "ac", "acy", "beta", "se", "tstat", "p", "pheno")
afs = read.table(paste0(dirsrv, "afs.tsv"), h=T, sep="\t")

# find any genome file (as produced by combine-imputed.sh)
# and get individual ID from this file
files = list.files(pattern="^allchr")
if(!length(files)){
	print("ERROR: this folder must contain a genotype file allchr_*.txt.gz")
	quit(1)
}
id = strsplit(files, "_")[[1]]
id = strsplit(id[2], "\\.")[[1]][1]

print(sprintf("working on individual with ID %s", i))

# read VCF
vcf = read.table(gzfile(files[i]))
colnames(vcf) = c("CHR", "POS", "rsHRC", "REF", "ALT", "gt", "ds")

# merge with clinvar
vcfc = inner_join(vcf, cv, by=c("CHR"="Chromosome", "POS"="Start",
								"REF"="ReferenceAllele", "ALT"="AlternateAllele"))
vcfc = select(vcfc, -one_of(c("Type", "LastEvaluated", "nsv.esv..dbVar.", "Assembly",
							  "Cytogenetic", "NumberSubmitters", "Guidelines", "SubmitterCategories",
							  "RCVaccession", "ChromosomeAccession", "Stop")))
vcfc = mutate(vcfc, Name = sapply(strsplit(Name, ":"), "[[", 2)) %>%
	mutate(Name = substring(Name, 3))

# drop frequent genotypes
vcfc = left_join(vcfc, afs, by=c("CHR", "POS", "REF", "ALT")) %>%
	mutate(pU = dbinom(round(ds), 2, af), pH = dbinom(round(ds), 2, AF_GLOBAL)) %>%
	filter(pU < 0.9 | pH < 0.9)
if(!nrow(vcfc)){
	print("ERROR: something went wrong, CV merge output is empty")
	quit(1)
}

fileCV = paste0("cv_", i, ".tsv")
write.table(vcfc, fileCV, quote=F, row.names=F, col.names=T, sep="\t")

# merge with summaries and get GRSs
vcf = inner_join(vcf, ukbb, by=c("CHR", "POS", "REF", "ALT")) %>%
	group_by(pheno) %>%
	summarize(GRS = sum(beta*ds))
if(!nrow(vcf)){
	print("ERROR: something went wrong, UKBB merge output is empty")
	quit(1)
}

fileGRS = paste0("grs_", i, ".tsv")
write.table(vcf, fileGRS, quote=F, row.names=F, col.names=T, sep="\t")

# copy both output files to server directory
system(paste("cp", fileCV, dirsrv))
system(paste("cp", fileGRS, dirsrv))
system(paste0("chgrp shiny ", dirsrv, fileCV))
system(paste0("chgrp shiny ", dirsrv, fileGRS))
system(paste0("chmod g+r ", dirsrv, fileCV))
system(paste0("chmod g+r ", dirsrv, fileGRS))

