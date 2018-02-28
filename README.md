# homegenome
Shiny app for viewing and analyzing 23andMe results

## Files
* UKBB manifest, incl. wget commands for autodownload
* UKBB heritability estimates
* UKBB phenotype collection notes
* own genome, Sanger-imputed to 22 separate autosomal VCFs
* external allele frequencies (HRC)

## Files created by this soft
* dictionary of categories
* login tables (w/ hashed passwords)
* UKBB summary statistics, LD-pruned
* encrypted personal GRS and ClinVar summaries
* various temporary intermediates

## Scripts
1. `filter-phenofiles.sh` to drop boring phenotypes (medications...)
2. `download-summaries.sh` to wget all the rest, and filter p<1e-5
3. `prune-ld.R` to clump SNPs based on distance and UKBB association results
4. `process-ukbb-cv.R` to summarize UKBB and ClinVar reference tables
5. `combine-imputed.sh` to extract dosages from your Sanger-imputed VCFs
6. `combine-imputed2.R` to find your GRSs and pathogenic markers
7. `app.R` - the actual UI and server for Shiny viewer

## User process
1. HRC-impute your genotypes in Sanger
2. Have R and packages `dplyr`, `tidyr`, `sodium`
3. Download and run `bash combine-imputed.sh`
4. Download and run `Rscript combine-imputed2.sh`
5. Browse to the app, login and enjoy
