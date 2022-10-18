# Calculates genotype likelihoods of low-coverage BAM files and stored them in VCF files.
# To change sample names it requires a 'samples.txt' file with the format:
# Old_sample_name New_sample_name
# The two fields are separated by a space.
# Especially important to ensure 1x and 30x vcf files have the same sample name on the header to run GLIMPSE_concordance.

# Change to use a different number of CPU threads.
threads=16

# Change 1x to 30x for high-coverage files.
coverage=1x

# Whole Genome
chromosomes=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 )

# Chromosome 22 only
# chromosomes=( chr22 )

# Calculation of genotype likelihoods
for BAM in seq/${coverage}/*.bam;
    do echo bam: ${BAM};
    for REGION in "${chromosomes[@]}";
        do VCF=1000G/1000G.${REGION}.sites.vcf.gz;
        TSV=1000G/1000G.${REGION}.sites.tsv.gz;
        REFGEN=reference/hg19.fa;
        SAMPLE=${BAM#"seq/"${coverage}"/"};
        SAMPLE=${SAMPLE%.bam};
        echo sample: ${SAMPLE};
        mkdir -p vcf/${coverage}/${SAMPLE};
        OUT=vcf/${coverage}/${SAMPLE}/${REGION}.vcf.gz;
        echo Calling ${REGION} of ${SAMPLE}...;
        bcftools mpileup --threads $threads -f ${REFGEN} -I -E -a 'FORMAT/DP' -T ${VCF} -r ${REGION} ${BAM} -Ou | bcftools call --threads $threads --ploidy GRCh37 -Aim -C alleles -T ${TSV} -Oz -o ${OUT};
        echo Indexing ${SAMPLE}...;
        bcftools index --threads $threads -f ${OUT};
    done;
done

chromosomeX=( chrX:60001-2699520 chrX:2699521-154931043 chrX:154931044-155260560 )

for BAM in seq/${coverage}/*.bam;
    do echo bam: ${BAM};
    for REGION in "${chromosomeX[@]}";
        do VCF=1000G/1000G.chrX.sites.vcf.gz;
        TSV=1000G/1000G.chrX.sites.tsv.gz;
        REFGEN=reference/hg19.fa;
        SAMPLE=${BAM#"seq/"${coverage}"/"};
        SAMPLE=${SAMPLE%.bam};
        echo sample: ${SAMPLE};
        mkdir -p vcf/${coverage}/${SAMPLE};
        OUT=vcf/${coverage}/${SAMPLE}/${REGION}.vcf.gz;
        echo Calling ${REGION} of ${SAMPLE}...;
        bcftools mpileup --threads $threads -f ${REFGEN} -I -E -a 'FORMAT/DP' -T ${VCF} -r ${REGION} ${BAM} -Ou | bcftools call --threads $threads --ploidy GRCh37 -Aim -C alleles -T ${TSV} -Oz -o ${OUT};
        echo Indexing ${SAMPLE}...;
        bcftools index --threads $threads -f ${OUT};
    done;
done
