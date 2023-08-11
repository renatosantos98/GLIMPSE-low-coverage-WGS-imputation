#!/usr/bin/env bash
# Imputes low-coverage data using 1000 Genomes as the reference panel.
# Specify which chromosomes to be imputed on input 2.
# Chromosome X imputation requires text file with sample ploidy.

# Change to use a different number of CPU threads.
threads=16

# Imputation of Autosomes
parallel --jobs ${threads} bash scripts/lib/glimpse_impute.sh ::: vcf/1x/22D* ::: {1..22} X:60001-2699520 X:154931044-155260560;

# Imputation of X non-PAR region
parallel --jobs ${threads} bash scripts/lib/glimpse_impute_X.sh ::: vcf/1x/22D*;

# Ligation of Chunks
parallel --jobs ${threads} bash scripts/lib/glimpse_ligate.sh ::: GLIMPSE_imputed/22D* ::: {1..22} X:60001-2699520 X:2699521-154931043 X:154931044-155260560;

# Merging of separate chromosome files into one single file
for SAMPLE in GLIMPSE_ligated/22D*;
    do ls ${SAMPLE}/*.merged.bcf > ${SAMPLE}/list.txt;
done

parallel --jobs ${threads} bash scripts/lib/merge_chr.sh ::: GLIMPSE_ligated/22D*


# Validation genome only
# Imputation of Autosomes
parallel --jobs ${threads} bash scripts/lib/glimpse_impute.sh ::: vcf/1x/IBS001/{bgi,illumina} ::: X:60001-2699520 X:154931044-155260560;

# Imputation of X non-PAR region
parallel --jobs ${threads} bash scripts/lib/glimpse_impute_X.sh ::: vcf/1x/IBS001/{bgi,illumina};

# Ligation of Chunks
parallel --jobs ${threads} bash scripts/lib/glimpse_ligate.sh ::: GLIMPSE_imputed/IBS001/{bgi,illumina} ::: X:60001-2699520 X:2699521-154931043 X:154931044-155260560;

# Merging of separate chromosome files into one single file
for SAMPLE in GLIMPSE_ligated/IBS001/{bgi,illumina};
    do ls ${SAMPLE}/*.merged.bcf > ${SAMPLE}/list.txt;
done

mkdir -p imputed_vcf;

parallel --jobs ${threads} bash scripts/lib/merge_chr.sh ::: GLIMPSE_ligated/IBS001/{bgi,illumina}
