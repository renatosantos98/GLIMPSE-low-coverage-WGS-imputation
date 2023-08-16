#!/usr/bin/env bash
# Prepares vcf files for ancestry analysis and performs principal component analysis using a GTM model.

# Change to use a different number of CPU threads.
threads=16

mkdir -p pca/{raw,pruned/{ALL,IBS},merged,extracted,severe_covid,1000G_pca,IBS_pca}

cd pca/raw

rsync -avzP rsync://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{1..22}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz{,.tbi} .
rsync -avzP rsync://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped .

reference="/mnt/e/Sarek/references/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta"

for chr in {1..22}; do
    input_vcf="ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    output_bcf="ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf"
    IBS_output_bcf="IBS.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf"
    # Normalize VCF file so multi-allelic calls are split and indels are left-aligned compared to reference genome
    bcftools norm --threads ${threads} -m-any --check-ref w -f "${reference}" "${input_vcf}" -Ou | \
    # Unset IDs will be set to CHROM:POS:REF:ALT*
    bcftools annotate --threads ${threads} -x ID -I +'%CHROM:%POS:%REF:%ALT' -Ou | \
    # Remove duplicate records
    bcftools norm --threads ${threads} -Ob --rm-dup both -o "${output_bcf}"
    # Index the resulting BCF file
    bcftools index -f --threads ${threads} "${output_bcf}"
    # Filter output file to include only IBS individuals
    bcftools view -S ../../scripts/lib/IBS_1000G.txt --force-samples -Ob -o "${IBS_output_bcf}" "${output_bcf}"
    # Index the resulting BCF file
    bcftools index -f --threads ${threads} "${IBS_output_bcf}"
done

# Convert the BCF files to PLINK format
for chr in {1..22}; do
    plink --noweb \
      --bcf ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf \
      --keep-allele-order \
      --vcf-idspace-to _ \
      --const-fid \
      --allow-extra-chr 0 \
      --split-x b37 no-fail \
      --make-bed \
      --out ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes

    plink --noweb \
      --bcf IBS.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf \
      --keep-allele-order \
      --vcf-idspace-to _ \
      --const-fid \
      --allow-extra-chr 0 \
      --split-x b37 no-fail \
      --make-bed \
      --out IBS.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes
done

find . -name "*.bim" | grep -e "ALL" > ForMergeALL.list
sed -i 's/.bim//g' ForMergeALL.list
find . -name "*.bim" | grep -e "IBS" > ForMergeIBS.list
sed -i 's/.bim//g' ForMergeIBS.list

# Merge projects into single PLINK files
plink --merge-list ForMergeALL.list --out ../merged/ALL --memory 122000
plink --merge-list ForMergeIBS.list --out ../merged/IBS --memory 122000

cd ../..
# Merge all filtered severe COVID patient VCF files
bcftools merge --threads ${threads} -r $(seq -s, 1 22) --missing-to-ref -Ou imputed_vcf/22D*.vcf.gz | \
# Normalize VCF file so multi-allelic calls are split and indels are left-aligned compared to reference genome
bcftools norm --threads ${threads} -m-any --check-ref w -f "${reference}" -Ou | \
# Unset IDs will be set to CHROM:POS:REF:ALT*
bcftools annotate --threads ${threads} -x ID -I +'%CHROM:%POS:%REF:%ALT' -Ou | \
# Remove duplicate records
bcftools norm --threads ${threads} -Ob --rm-dup both -o pca/severe_covid/severe_covid_merged.bcf
bcftools index --threads ${threads} pca/severe_covid/severe_covid_merged.bcf

# Convert the BCF files to PLINK format
cd pca/severe_covid
# plink --bcf severe_covid_merged.bcf --make-bed --out severe_covid_merged

# Determine common variants between the severe COVID and 1000 Genomes and IBS datasets
plink --bfile severe_covid_merged --bmerge ../merged/ALL --out common_variants_severe_covid_ALL --memory 122000
plink --bfile severe_covid_merged --bmerge ../merged/IBS --out common_variants_severe_covid_IBS --memory 122000

# Extract the common variants
plink --bfile severe_covid_merged --extract common_variants_severe_covid_ALL.bim --make-bed --out ../extracted/severe_covid_extracted_ALL --memory 122000
plink --bfile severe_covid_merged --extract common_variants_severe_covid_IBS.bim --make-bed --out ../extracted/severe_covid_extracted_IBS --memory 122000
plink --bfile ../merged/ALL --extract common_variants_severe_covid_ALL.bim --make-bed --out ../extracted/ALL_extracted --memory 122000
plink --bfile ../merged/IBS --extract common_variants_severe_covid_IBS.bim --make-bed --out ../extracted/IBS_extracted --memory 122000

# Create common variants datasets for PCA
cd ../extracted
plink --bfile severe_covid_extracted_ALL --bmerge ALL_extracted --out severe_covid_ALL_merged --memory 122000
plink --bfile severe_covid_extracted_IBS --bmerge IBS_extracted --out severe_covid_IBS_merged --memory 122000

# Prune variants by MAF and VIF
cd ..
plink --bfile extracted/severe_covid_ALL_merged --maf 0.10 --indep 50 5 1.5 --out pruned/severe_covid_ALL_pruned --memory 122000
plink --bfile extracted/severe_covid_ALL_merged --extract pruned/severe_covid_ALL_pruned.prune.in --make-bed --out pruned/severe_covid_ALL_pruned --memory 122000
plink --bfile extracted/severe_covid_IBS_merged --maf 0.10 --indep 50 5 1.5 --out pruned/severe_covid_IBS_pruned --memory 122000
plink --bfile extracted/severe_covid_IBS_merged --extract pruned/severe_covid_IBS_pruned.prune.in --make-bed --out pruned/severe_covid_IBS_pruned --memory 122000

# Perform PCA
cd 1000G_pca
plink --bfile ../pruned/severe_covid_ALL_pruned --pca --memory 122000

cd ../IBS_pca
plink --bfile ../pruned/severe_covid_IBS_pruned --pca --memory 122000
cd ../..

python scripts/lib/pca.py
