# Imputes low-coverage data using 1000 Genomes as the reference panel.
# Requires a /maps directory with genetic maps for each chromosome.

# Change to use a different number of CPU threads.
threads=16

# Whole Genome
chromosomes=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 )

# Chromosome 22 only
# chromosomes=( chr22 )

GLIMPSE_phase Imputation
for SAMPLE in vcf/1x/*;
    do echo sample: ${SAMPLE};
    for REGION in "${chromosomes[@]}";
        do VCF=${SAMPLE}/${REGION}.vcf.gz;
        echo vcf: ${VCF};
        REF=1000G/1000G.${REGION}.bcf;
        MAP=maps/${REGION}.b37.gmap;
        mkdir -p GLIMPSE_imputed/${SAMPLE#"vcf/1x/"}/${REGION};
        while IFS="" read -r LINE || [ -n "$LINE" ];
            do printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
            IRG=$(echo $LINE | cut -d" " -f3)
            ORG=$(echo $LINE | cut -d" " -f4)
            OUT=GLIMPSE_imputed/${SAMPLE#"vcf/1x/"}/${REGION}/${ID}.bcf
            GLIMPSE_phase --input ${VCF} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT};
            bcftools index $threads -f ${OUT};
        done < chunks/chunks.${REGION}.txt;
    done;
done

# Create list of imputed chunks for each chromosome and ligate them
for SAMPLE in GLIMPSE_imputed/*;
    do
    for REGION in "${chromosomes[@]}";
        do ls ${SAMPLE}/${REGION}/*.bcf > ${SAMPLE}/${REGION}/list.txt;
        mkdir -p GLIMPSE_ligated/${SAMPLE#"GLIMPSE_imputed/"};
        GLIMPSE_ligate --thread ${threads} --input ${SAMPLE}/${REGION}/list.txt --output GLIMPSE_ligated/${SAMPLE#"GLIMPSE_imputed/"}/${REGION}.merged.bcf;
        bcftools index --threads ${threads} -f GLIMPSE_ligated/${SAMPLE#"GLIMPSE_imputed/"}/${REGION}.merged.bcf;
    done;
done
