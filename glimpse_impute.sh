# Imputes low-coverage data using 1000 Genomes as the reference panel.
# Requires a /maps directory with genetic maps for each chromosome.

# Change to use a different number of CPU threads.
threads=1

# GLIMPSE_phase Imputation
SAMPLE=${1}
SAMPLENAME=${SAMPLE#"vcf/1x/"}
REGION=${2}
VCF=${SAMPLE}/${REGION}.vcf.gz;
REF=1000G/1000G.${REGION%":"*}.bcf;
MAP=maps/${REGION%":"*}.b37.gmap;
mkdir -p GLIMPSE_imputed/${SAMPLENAME}/${REGION};
while IFS="" read -r LINE || [ -n "$LINE" ];
    do printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
    IRG=$(echo $LINE | cut -d" " -f3)
    ORG=$(echo $LINE | cut -d" " -f4)
    OUT=GLIMPSE_imputed/${SAMPLENAME}/${REGION}/${ID}.bcf
    GLIMPSE_phase --thread ${threads} --input ${VCF} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT};
    bcftools index --threads ${threads} -f ${OUT};
done < chunks/chunks.${REGION}.txt;
