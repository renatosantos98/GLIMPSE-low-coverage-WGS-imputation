# Aligns fastq files to the hg19 reference genome.
# Requires fastq files to be in the directory /fastq.
# Saves aligned files to the /seq/1x directory.

# Change to use a different number of CPU threads.
threads=16
sortthreads=`expr $threads / 2`

# Change 1x to 30x for high-coverage files.
coverage=1x

mkdir -p seq/1x

for finput in fastq/*.fastq.gz;
    do samoutput=${finput%.fastq.gz}.sam;
    bamoutput=${samoutput#fastq/};
    bamoutput=seq/${coverage}/${bamoutput%.sam}.bam;
    echo Aligning $finput...;
    bwa mem -t $threads reference/hg19.fa $finput > $samoutput;
    echo Sorting $samoutput...;
    samtools view --threads $sortthreads -b -h --fai-reference reference/hg19.fa.fai $samoutput | samtools sort --threads $sortthreads -o $bamoutput;
    rm -f $samoutput;
    echo Indexing $bamoutput...;
    samtools index --threads $threads $bamoutput;
done
