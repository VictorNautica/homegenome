options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
setwd("~/Documents/23andme/")

# Various snippets for processing for UKBB and clinvar files,
# before merging them with individual VCFs.

cats = read.table("UKBB_pheno_categories.tsv", h=T, sep="\t", quote="", fill = T)

## category decoding:
catDict = data.frame(codes=c("c", "nc", "d", "cd", "b", "f", "a", "sd", "ps", "e", "ey", "p", "vt",
							 "ci", "h", "ia", "ii", "ik", "io", "is",
							 "d0", "nc0", "0"),
				 categories=c("Cancer code, self-reported", "Non-cancer illness code, self-reported",
				 			 "Diagnoses - main ICD10", "Underlying (primary) cause of death",
				 			 "Lifestyle, behavior", "Fitness, activity", "Anthropometry",
				 			 "Smoking, drinking", "Psychology, psychiatry", "Intelligence, education",
				 			 "Eyesight", "Pain", "Vitamin and mineral supplements",
				 			 "Cancer, ICD10", "Hearing", "Infections and parasites, ICD10",
				 			 "Blood and cardiovascular diseases", "Teeth and digestive disorders",
				 			 "Reproductive health and childbirth", "Fractures, external injuries",
				 			 "Diagnoses, ICD10, not classified further",
				 			 "Illnesses, self-reported, not classified further", "Various"),
				 definition=c(rep("(from biobank)", 4),
				 			 "Various questions on life history and lifestyle choices",
				 			 "Questions on activity and exericse habits, some fitness measurements",
				 			 "Body size, weight, obesity measurements, some parameters obtained during exercise",
				 			 "", "Symptoms of depression, suicide, anxiety, behavior disorders",
				 			 "Cognitive performance measures, qualifications, work type",
				 			 "Eye diseases, glass usage", "", "",
				 			 "As reported in causes of death or ICD10 diagnoses",
				 			 "Reported hearing strength and ear disorders", "", "ICD10 or self-reported",
				 			 "ICD10 or self-reported",
				 			 "Reproductive diseases, fertility, general questions on pregnancy",
				 			 "", "", "", ""))
write.table(catDict, file="UKBB_category_dict.tsv", row.names=F, col.names=T, quote=F, sep="\t")


## get allele frequencies in UKBB and HRC
ukbb = read.table("all_effects.tsv", h=F, sep="\t")
colnames(ukbb) = c("CHR", "POS", "REF", "ALT", "rs",
				   "n", "ac", "acy", "beta", "se", "tstat", "p", "pheno")
# store all SNPs
write.table(unique(ukbb[,1:2]), "all_effects_snps.txt", sep="\t", quote=F, row.names=F, col.names=F)

## (extract these SNPs from 1000G in CEU/FIN/GBR and get allele freqs)


ref = read.table(pipe("cut -f2 all_effects.tsv | sort | uniq |
			  awk 'FNR==NR{a[$1]; next} $2 in a' - ~/data/geno/references/HRC.r1-1.GRCh37.wgs.mac5.sites.tab"))
colnames(ref) = c("CHR", "POS", "rsHRC", "REF", "ALT", "acHRC", "anHRC",
				  "AF_GLOBAL", "acNO1000G", "anNO1000G", "afNO1000G", "AA")
ref = filter(ref, CHR!="X")
ref$CHR = as.numeric(ref$CHR)
ukbb2 = inner_join(ukbb, ref, by=c("CHR", "POS", "REF", "ALT"))

ukbb2$af = ukbb2$ac/ukbb2$n/2
ukbb2$af[is.na(ukbb2$af)] = 0

afs = group_by(ukbb2, CHR, POS, REF, ALT) %>%
	summarize(af=min(af), AF_GLOBAL=min(AF_GLOBAL))
write.table(afs, "afs.tsv", col.names=T, row.names=F, quote=F, sep="\t")

# get expect distr for GRSs in UKBB
ukbbSum = group_by(ukbb2, pheno) %>%
	summarize(GRSmeanU = sum(2 * beta * af),
			  GRSsdU = sqrt(sum(2 * beta^2 * af * (1-af))),
			  GRSmeanH = sum(2 * beta * AF_GLOBAL),
			  GRSsdH = sqrt(sum(2 * beta^2 * AF_GLOBAL * (1-AF_GLOBAL))))

# fix phenotype category tags and attach (together with pheno description)
cats$Categories[cats$Categories=="d"] = "d,d0"
cats$Categories[cats$Categories=="nc"] = "nc,nc0"
cats$Categories[cats$Categories==""] = "0"
ukbbSum = left_join(ukbbSum, cats, by=c("pheno"="Phenotype.code"))

write.table(ukbbSum, "ukbb_summarizedgrs.tsv", col.names=T, row.names=F, quote=F, sep="\t")

### LOAD CLINVAR
cv = read.table("clinvarsummary.txt", h=T, sep="\t", quote="", comment.char = "")
cv = filter(cv, Chromosome %in% 1:22)
cv$Chromosome = as.numeric(cv$Chromosome)
cv$PhenotypeList = gsub("not specified[;]*", "", cv$PhenotypeList)
cv$PhenotypeList = gsub("not provided[;]*", "", cv$PhenotypeList)
cv = filter(cv, ReferenceAllele %in% c("A", "C", "G", "T"),
			AlternateAllele %in% c("A", "C", "G", "T")) %>%
	filter(Assembly=="GRCh37") %>%
	filter(PhenotypeList!="")

# remove variants of no interesting consequences
cv = rename(cv, Clinical = ClinicalSignificance, Gene = GeneSymbol) %>%
	filter(Clinical != "Uncertain significance",
		   !Clinical %in% c("-", "not provided", "other"),
		   !grepl("^Benign", Clinical),
		   !grepl("^Likely benign", Clinical))

# rename categories
cv$Clinical[grepl("^Conflicting", cv$Clinical)] = "unclear"
cv$Clinical[grepl("^Uncertain", cv$Clinical)] = "unclear"
cv$Clinical[grepl("protective", cv$Clinical)] = "protective"
cv$Clinical[grepl("^Pathogenic", cv$Clinical)] = "pathogenic"
cv$Clinical[grepl("^Likely pathogenic", cv$Clinical)] = "pathogenic"
cv$Clinical[grepl("^Affects", cv$Clinical)] = "association"
cv$Clinical[grepl("^drug", cv$Clinical)] = "drug response"
cv$Clinical[grepl("^risk", cv$Clinical)] = "association"
cv$Clinical[grepl("^association", cv$Clinical)] = "association"

write.table(cv, "clinvar_filtered.tsv", col.names=T, row.names=F, quote=F, sep="\t")
