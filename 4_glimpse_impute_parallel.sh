# Imputes low-coverage data using 1000 Genomes as the reference panel.
# Specify which chromosomes to be imputed on input 2.
# Chromosome X imputation requires text file with sample ploidy.

parallel --jobs 16 bash glimpse_impute.sh ::: vcf/1x/* ::: chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX:60001-2699520 chrX:154931044-155260560 && parallel --jobs 16 bash glimpse_impute_X.sh ::: vcf/1x/* ::: chrX:2699521-154931043
