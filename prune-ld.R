#!/usr/bin/Rscript
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)

setwd("~/Documents/23andme/assoc/")
origs = list.files(pattern = "assoc.tsv$")

allphenos = NULL
for(orig in origs){
	print(sprintf("working on file %s", orig))
	phenonum = strsplit(orig, "\\.")[[1]][1]
	
	df = read.table(orig, h=F)
	df$pheno = phenonum
	df = separate(df, V1, c("chr", "pos", "ref", "alt"), convert = T) %>%
		arrange(chr, pos)
	df = group_by(df, chr) %>%
		mutate(d=pos-lag(pos, default=-Inf), clumpid=cumsum(d>1e5)) %>%
		group_by(chr, clumpid) %>%
		filter(rank(V9, ties.method = "random")==1)
	print(sprintf("%d lines remaining", nrow(df)))
	
	allphenos = bind_rows(allphenos, df[,1:13])
}
colnames(allphenos)[5:12] = c("rsid", "nsamples", "ac", "acy", "beta", "se", "tstat", "p")
write.table(allphenos, "../all_effects.tsv", quote=F, col.names=F, row.names=F, sep="\t")

