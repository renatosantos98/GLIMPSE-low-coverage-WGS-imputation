# Create list of imputed chunks for each chromosome and ligate them

chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX:60001-2699520 chrX:2699521-154931043 chrX:154931044-155260560)

# Change to use a different number of CPU threads.
threads=1

for SAMPLE in GLIMPSE_imputed/*;
    do
    for REGION in "${chromosomes[@]}";
        do ls ${SAMPLE}/${REGION}/*.bcf > ${SAMPLE}/${REGION}/list.txt;
        mkdir -p GLIMPSE_ligated/${SAMPLE#"GLIMPSE_imputed/"};
        GLIMPSE_ligate --thread ${threads} --input ${SAMPLE}/${REGION}/list.txt --output GLIMPSE_ligated/${SAMPLE#"GLIMPSE_imputed/"}/${REGION}.merged.bcf;
        bcftools index --threads ${threads} -f GLIMPSE_ligated/${SAMPLE#"GLIMPSE_imputed/"}/${REGION}.merged.bcf;
    done;
done
