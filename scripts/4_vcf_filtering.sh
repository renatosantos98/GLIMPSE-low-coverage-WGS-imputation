#!/usr/bin/env bash
# Filters variants in VCF files for MAF and FORMAT/GP scores.

# Change to use a different number of CPU threads.
threads=16

mkdir -p filtered_vcf

# Keeps only imputed variants with Minor Allele Frequency higher than 2% and Genotype Probability higher than 80%.
for VCF in imputed_vcf/22D*.vcf.gz;
    do OUT=filtered_vcf/${VCF#"imputed_vcf/"};
    bcftools view --threads ${threads} -i "INFO/RAF[0]>0.02 && INFO/RAF[0]<0.95 && FORMAT/GP[0:*]>0.80" -Oz9 -o ${OUT} ${VCF};
    bcftools index --threads 2 -f ${OUT};
    bcftools stats --threads ${threads} -v -r "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X" ${OUT} > ${OUT%.vcf.gz}.stats.txt;
done

# Validation
for PLATFORM in bgi illumina;
    do bcftools view --threads ${threads} -i "INFO/RAF[0]>0.02 && INFO/RAF[0]<0.95 && FORMAT/GP[0:*]>0.80" -Oz9 -o filtered_vcf/IBS001/${PLATFORM}.vcf.gz imputed_vcf/IBS001/${PLATFORM}.vcf.gz;
    bcftools index --threads 2 -f filtered_vcf/IBS001/${PLATFORM}.vcf.gz
    bcftools stats --threads ${threads} -v -r "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X" filtered_vcf/IBS001/${PLATFORM}.vcf.gz > filtered_vcf/IBS001/${PLATFORM}.stats.txt;
done;

# Array to store the extracted numbers
numbers=()

rm -f filtered_vcf/snps.txt

# Loop through all files matching the glob pattern
for file in filtered_vcf/22D*.stats.txt; do
    # Read line 26 and extract the rightmost number
    number=$(sed -n '26s/.*\t\([0-9]\+\)$/\1/p' "$file")

    # Add the number to the array
    numbers+=("$number")

    # Output the file name and the extracted number to the output file
    echo "File: $file, Number of SNPs: $number" >> filtered_vcf/snps.txt
done

# Calculate the average number
sum=0
for number in "${numbers[@]}"; do
    sum=$(($sum + $number))
done
average=$(awk "BEGIN {printf \"%.2f\", $sum / ${#numbers[@]}}")

# Calculate the 95% confidence interval
n=${#numbers[@]}
mean=$(awk "BEGIN {printf \"%.2f\", $sum / $n}")
variance=0
for number in "${numbers[@]}"; do
    diff=$(echo "$number - $mean" | bc)
    variance=$(echo "$variance + ($diff * $diff)" | bc)
done
stddev=$(awk "BEGIN {printf \"%.2f\", sqrt($variance / $n)}")
margin=$(awk "BEGIN {printf \"%.2f\", 1.96 * ($stddev / sqrt($n))}")
lower_limit=$(awk "BEGIN {printf \"%.2f\", $mean - $margin}")
upper_limit=$(awk "BEGIN {printf \"%.2f\", $mean + $margin}")

# Output the average number and the confidence interval
echo "Average Number: $average" >> filtered_vcf/snps.txt
echo "95% Confidence Interval: [$lower_limit, $upper_limit]" >> filtered_vcf/snps.txt
