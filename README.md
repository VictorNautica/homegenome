# homegenome
Shiny app for viewing and analyzing 23andMe results

## Files
* UKBB manifest, incl. wget commands for autodownload
* UKBB heritability estimates
* UKBB phenotype collection notes
* own genome, Sanger-imputed to 22 separate autosomal VCFs
* tables to decode your own IDs
* dictionary of categories
* External allele frequencies (HRC)

## Scripts
1. `filter-phenofiles.sh` to drop boring phenotypes (medications...)
2. `download-summaries.sh` to wget all the rest, and filter p<1e-5
3. `prune-ld.R` to clump SNPs based on distance and UKBB association results
4. `combine-imputed.sh` to extract dosages from Sanger-imputed VCFs
5. `combine-imputed2.R` to summarize your dosages into GRSs based on UKBB results
6. `app.R` - the actual UI and server for Shiny viewer
