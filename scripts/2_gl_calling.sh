#!/usr/bin/env bash
# Calculates genotype likelihoods of low-coverage CRAM files and stores them in VCF files.
# To change sample names it requires a 'samples.txt' file in the scripts folder with the format:
# "Old_sample_name" "New_sample_name"
# The two fields are separated by a space.
# Especially important to ensure 1x and 30x vcf files have the same sample name on the header to run GLIMPSE_concordance.

# Change to use a different number of CPU threads.
threads=16

# Used to organise validation files into different coverage levels.
coverage=1x

# Calculation of genotype likelihoods

# Severe Covid Samples
for BAM in bam/22D*/22D*.recal.cram;
    do echo bam: ${BAM};
    SAMPLE=${BAM#"bam/22D"*"/"};
    SAMPLE=${SAMPLE%.recal.cram};
    echo sample: ${SAMPLE};
    mkdir -p vcf/${coverage}/${SAMPLE};
    grep "_${SAMPLE}" scripts/lib/samples_sex.txt > ploidy.txt
    for REGION in {1..22} X;
        do VCF=1000G/1000G.chr${REGION%":"*}.sites.vcf.gz;
        TSV=1000G/1000G.chr${REGION%":"*}.sites.tsv.gz;
        REFGEN=/mnt/e/Sarek/references/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta;
        OUT=vcf/${coverage}/${SAMPLE}/${REGION}.vcf.gz;
        echo Calling chromosome ${REGION} of ${SAMPLE}...;
        bcftools mpileup --threads ${threads} -f ${REFGEN} -I -E -a 'FORMAT/DP' -T ${VCF} -r ${REGION} ${BAM} -Ou | bcftools call --threads ${threads} --ploidy GRCh37 -S ploidy.txt -Aim -C alleles -T ${TSV} -Oz -o ${OUT};
        echo Indexing ${OUT}...;
        bcftools index --threads 4 -f ${OUT};
    done;
    rm -f ploidy.txt;
done;

coverage=("1x" "30x")
# Validation genome only
for coverage in "${coverage[@]}";
do for PLATFORM in bgi illumina;
        do for BAM in bam/${PLATFORM}/${coverage}/*.cram;
            do echo bam: ${BAM};
            SAMPLE=${BAM#bam/${PLATFORM}/${coverage}/};
            SAMPLE=${SAMPLE%_${PLATFORM}.cram};
            echo sample: ${SAMPLE};
            mkdir -p vcf/${coverage}/${SAMPLE}/${PLATFORM}/;
            grep "${SAMPLE}" scripts/lib/samples_sex.txt > ploidy.txt
            for REGION in X;
                do VCF=1000G/1000G.chr${REGION%":"*}.sites.vcf.gz;
                TSV=1000G/1000G.chr${REGION%":"*}.sites.tsv.gz;
                REFGEN=/mnt/e/Sarek/references/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta;
                OUT=vcf/${coverage}/${SAMPLE}/${PLATFORM}/${REGION}.vcf.gz;
                echo Calling chromosome ${REGION} of ${SAMPLE}...;
                bcftools mpileup --threads ${threads} -f ${REFGEN} -I -E -a 'FORMAT/DP' -T ${VCF} -r ${REGION} ${BAM} -Ou | bcftools reheader --threads ${threads} -s scripts/lib/samples.txt | bcftools call --threads ${threads} --ploidy GRCh37 -S ploidy.txt -Aim -C alleles -T ${TSV} -Oz -o ${OUT};
                echo Indexing ${OUT}...;
                bcftools index --threads 4 -f ${OUT};
            done;
            rm -f ploidy.txt;
        done;
    done;
done;
