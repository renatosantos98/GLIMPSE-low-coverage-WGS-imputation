#!/usr/bin/env bash
# Imputes low-coverage data using 1000 Genomes as the reference panel.
# Requires a /maps directory with genetic maps for each chromosome.
# Will raise a Segmentation fault error if the sample is haploid. This is known to not affect the imputation process. See https://github.com/odelaneau/GLIMPSE/issues/92.

# Change to use a different number of CPU threads.
threads=1

# GLIMPSE_phase Imputation
SAMPLE=${1}
SAMPLENAME=${SAMPLE#"vcf/1x/"}
VCF=${SAMPLE}/X.vcf.gz;
REF=1000G/1000G.chrX.bcf;
MAP=maps/X:2699521-154931043.b37.gmap;
mkdir -p GLIMPSE_imputed/${SAMPLENAME}/X:2699521-154931043;
while IFS="" read -r LINE || [ -n "$LINE" ];
    do printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
    IRG=$(echo $LINE | cut -d" " -f3)
    ORG=$(echo $LINE | cut -d" " -f4)
    OUT=GLIMPSE_imputed/${SAMPLENAME}/X:2699521-154931043/${ID}.bcf
    GLIMPSE_phase --thread ${threads} --input ${VCF} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT} --samples-file scripts/lib/samples_X_ploidy.txt;
    bcftools index --threads ${threads} -f ${OUT};
done < chunks/chunks.X:2699521-154931043.txt;
