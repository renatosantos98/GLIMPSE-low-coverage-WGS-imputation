# Creates the conda environment, installs the required tools, and prepares a genome reference and

# Change to use a different number of CPU threads.
threads=16

# Create environment - requires conda to be installed
echo Creating environment...
conda create -n glimpse
conda activate glimpse
mamba install glimpse-bio bwa samtools bcftools matplotlib numpy

# Download hg19 files into /reference
echo Downloading hg19 reference...
mkdir -p reference
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ reference/

# Unzip and concatenate hg19 into a single file
echo Creating hg19 reference file...
gunzip -c reference/chr*.fa.gz > reference/hg19.fa

# Remove useless files
rm -f reference/chr*.fa.gz

# Index hg19 reference
echo Indexing hg19 reference file...
bwa index -a bwtsw reference/hg19.fa
samtools faidx reference/hg19.fa

# Create /1000G directory
mkdir -p 1000G

# Download 1000 Genomes reference panel files into /1000G
echo Downloading 1000 Genomes reference panel...
rsync -avzP rsync://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr*.vcf.gz{,.tbi} 1000G/

# Reference Panel QC
for finput in 1000G/ALL.chr*.vcf.gz;
    do chromosome=`echo $finput | grep -oP chr[[:alnum:]]+`;
    foutput=1000G/1000G.$chromosome.temp.bcf;
    echo \(1000 Genomes\) Performing QC on $chromosome...;
    bcftools view --threads $threads -m 2 -M 2 -v snps $finput -Ob -o $foutput;
done

# Create file used to convert chromosome names to GRCh37 standard
echo -e "1   chr1\n2   chr2\n3   chr3\n4   chr4\n5   chr5\n6   chr6\n7   chr7\n8   chr8\n9   chr9\n10   chr10\n11   chr11\n12   chr12\n13   chr13\n14   chr14\n15   chr15\n16   chr16\n17   chr17\n18   chr18\n19   chr19\n20   chr20\n21   chr21\n22   chr22\nX   chrX\nY   chrY\nMT  chrM" > annotate.chr.txt

# Convert chromosome names on 1000G bcf files
for finput in 1000G/1000G.*.temp.bcf;
    do foutput=${finput%.temp.bcf}.bcf;
    echo \(1000 Genomes\) Changing chromosome names on $finput...;
    bcftools annotate --threads $threads --rename-chrs annotate.chr.txt -Ob -o $foutput $finput;
    bcftools index --threads $threads -f $foutput;
    rm -f $finput;
done

rm -f 1000G/ALL.chr*.vcf.gz{,.tbi}

# Extract variable positions in the reference panel
for finput in 1000G/1000G.*.bcf;
    do foutput=${finput%.bcf}.sites.vcf.gz;
    region=`echo $finput | grep -oP chr[[:alnum:]]+`;
    echo \(1000 Genomes GL Calculation VCF\) Processing $finput...;
    bcftools view --threads $threads -G -m 2 -M 2 -v snps -r $region $finput -Oz -o $foutput;
    bcftools index --threads $threads -f $foutput;
done

for finput in 1000G/*sites.vcf.gz;
    do foutput=${finput%.vcf.gz}.tsv.gz;
    echo \(1000 Genomes GL Calculation TSV\) Processing $finput...;
    bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' $finput | bgzip -@ 16 -c > $foutput;
    tabix -s1 -b2 -e2 $foutput;
done

# Prepare GLIMPSE Imputation Chunks
chromosomes=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 )

for region in "${chromosomes[@]}";
    do GLIMPSE_chunk --thread 16 --input 1000G/1000G.${region}.sites.vcf.gz --region ${region} --window-size 2000000 --buffer-size 200000 --output chunks/chunks.${region}.txt;
done

chromosomeX=( chrX:60001-2699520 chrX:2699521-154931043 chrX:154931044-155260560 )

for region in "${chromosomeX[@]}";
    do GLIMPSE_chunk --thread 16 --input 1000G/1000G.chrX.sites.vcf.gz --region ${region} --window-size 2000000 --buffer-size 200000 --output chunks/chunks.${region}.txt;
done
