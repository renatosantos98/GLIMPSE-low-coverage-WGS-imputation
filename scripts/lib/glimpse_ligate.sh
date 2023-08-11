#!/usr/bin/env bash
# Create list of imputed chunks for each chromosome and ligate them
# Change to use a different number of CPU threads.
threads=1

SAMPLE=${1}
REGION=${2}
ls ${SAMPLE}/${REGION}/*.bcf > ${SAMPLE}/${REGION}/list.txt;
mkdir -p GLIMPSE_ligated/${SAMPLE#"GLIMPSE_imputed/"};
GLIMPSE_ligate --thread ${threads} --input ${SAMPLE}/${REGION}/list.txt --output GLIMPSE_ligated/${SAMPLE#"GLIMPSE_imputed/"}/${REGION}.merged.bcf;
bcftools index --threads ${threads} -f GLIMPSE_ligated/${SAMPLE#"GLIMPSE_imputed/"}/${REGION}.merged.bcf;
