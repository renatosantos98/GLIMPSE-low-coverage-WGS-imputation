# Aligns fastq files to the hg19 reference genome.
# Requires fastq files to be in the directory /fastq.
# Saves aligned files to the /seq/1x directory.

# Change to use a different number of CPU threads.
threads=16
sortthreads=`expr $threads / 2`

# Change 1x to 30x for high-coverage files.
coverage=1x

mkdir -p seq/${coverage}

for fastq1 in fastq/22*1PFW/clean_data/22*_1.fq.gz;
    do fastq2=${fastq1%1.fq.gz}2.fq.gz;
    echo fastq1: ${fastq1};
    echo fastq2: ${fastq2};
    samoutput=${fastq1#"fastq/"};
    samoutput=${samoutput%"/clean_data/"22*_1.fq.gz}.sam;
    echo sam: ${samoutput};
    bamoutput=${samoutput#"fastq/"};
    bamoutput=seq/${coverage}/${bamoutput%.sam}.bam;
    echo bam: ${bamoutput};
    echo Aligning ${samoutput%.sam}...;
    bwa mem -t $threads reference/hg19.fa $fastq1 $fastq2 > $samoutput;
    echo Converting and sorting $samoutput to bam format...;
    samtools view --threads $sortthreads -b -h --fai-reference reference/hg19.fa.fai $samoutput | samtools sort --threads $sortthreads -o $bamoutput;
    rm -f $samoutput;
    echo Indexing $bamoutput...;
    samtools index -@ $threads $bamoutput;
done
