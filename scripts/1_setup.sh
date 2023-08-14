#!/usr/bin/env bash
# All scripts are expected to run from the base level directory, above "scripts/" (i.e., "../scripts").
# Prepares the imputation reference panel and imputation chunks.

# Change to use a different number of CPU threads.
threads=16

# Create /1000G directory
mkdir -p 1000G

# Download 1000 Genomes reference panel files into /1000G
echo Downloading 1000 Genomes reference panel...
rsync -avzP rsync://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{1..22}.phase3*.genotypes.vcf.gz{,.tbi} 1000G/
rsync -avzP rsync://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3*.genotypes.vcf.gz{,.tbi} 1000G/

# Reference Panel QC
echo Performing 1000G reference panel QC...
for finput in 1000G/ALL.chr*.vcf.gz;
    do chromosome=`echo $finput | grep -oP chr[[:alnum:]]+`;
    foutput=1000G/1000G.$chromosome.bcf;
    echo \(1000 Genomes\) Performing QC on $chromosome...;
    bcftools view --threads $threads -m 2 -M 2 -v snps $finput -Ob -o $foutput;
    bcftools index --threads 4 -f $foutput;
done

# Remove useless files
rm -f 1000G/ALL.chr*.vcf.gz{,.tbi}

# Extract variable positions on the reference panel
echo Extracting variable positions on the reference panel...
for finput in 1000G/1000G.*.bcf;
    do foutput=${finput%.bcf}.sites.vcf.gz;
    region=`echo $finput | grep -oP chr[[:alnum:]]+`;
    echo \(1000 Genomes GL Calculation VCF\) Processing $finput...;
    bcftools view --threads $threads -G -m 2 -M 2 -v snps -r ${region#chr} $finput -Oz -o $foutput;
    bcftools index --threads 4 -f $foutput;
done

for finput in 1000G/*sites.vcf.gz;
    do foutput=${finput%.vcf.gz}.tsv.gz;
    echo \(1000 Genomes GL Calculation TSV\) Processing $finput...;
    bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' $finput | bgzip -@ 16 -c > $foutput;
    tabix -s1 -b2 -e2 $foutput;
done

# Create chunks directory
mkdir -p chunks

# Prepare GLIMPSE Imputation Chunks
echo Creating GLIMPSE imputation chunks...
for region in {1..22};
    do echo Creating chr$region imputation chunks...;
    GLIMPSE_chunk --thread 16 --input 1000G/1000G.chr${region}.sites.vcf.gz --region ${region} --window-size 2000000 --buffer-size 200000 --output chunks/chunks.${region}.txt;
done

for region in X:60001-2699520 X:2699521-154931043 X:154931044-155260560;
    do echo Creating chr$region imputation chunks...;
    GLIMPSE_chunk --thread 16 --input 1000G/1000G.chrX.sites.vcf.gz --region ${region} --window-size 2000000 --buffer-size 200000 --output chunks/chunks.${region}.txt;
done
