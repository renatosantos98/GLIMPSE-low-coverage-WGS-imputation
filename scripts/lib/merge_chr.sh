#!/usr/bin/env bash
# Create list of imputed chunks for each chromosome and ligate them
# Change to use a different number of CPU threads.
threads=1

SAMPLE=${1#"GLIMPSE_ligated/"};
OUT=imputed_vcf/${SAMPLE}.vcf.gz;
bcftools concat --threads ${threads} -Ou --file-list ${1}/list.txt | bcftools sort -m 1G -Oz -o ${OUT};
bcftools index -f ${OUT};
