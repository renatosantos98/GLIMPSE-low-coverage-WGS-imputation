# Scripts

This directory contains the source code used for producing the results and figures for the publication. It consists of a pipeline of 6 shell scripts, which make use of supporting files and scripts inside the [`lib`](lib) directory. A citable permanent version of these scripts is archived in a [Figshare repository](https://doi.org/10.25452/figshare.plus.21679799).

## Imputation Pipeline

The following are the functions of each script in the imputation pipeline developed.

[`1_setup.sh`](1_setup.sh) - Sets up the imputation reference panel and chunks. Downloads the 1000 Genomes genotypes and extracts variable positions from it as files that can read by `bcftools` and `GLIMPSE1`. Then calculates imputation chunks out of the variable positions extracted before.

[`2_gl_calling.sh`](2_gl_calling.sh) - Calculates genotype likelihoods with `bcftools mpileup` for given CRAM files and saves them to VCF files. Uses the file [`lib/samples_sex.txt`](lib/samples_sex.txt) to determine chromosome X ploidy for each sample. Also uses the file [`lib/samples.txt`](lib/samples.txt) to normalise the sample name VCF header field in the validation files with `bcftools reheader`, required to run `GLIMPSE_concordance`.

[`3_glimpse_impute_parallel.sh`](3_glimpse_impute_parallel.sh) - Imputes calculated genotypes for each chunk with `GLIMPSE1`, using the 1000 Genomes data as the reference panel, and ligates chunks together to create separate files for each chromosome. Then, merges chromosomes together to create whole genome (chromosomes 1 to 22 and X) files.

[`4_vcf_filtering.sh`](4_vcf_filtering.sh) - Filters variants in the imputed VCF files by allele frequency (INFO/RAF higher than 2% and lower than 95%) and genotype probabilities (FORMAT/GP higher than 80%) [^1]. After all the imputed VCF files are filtered, calculates the average number of SNVs and 95% confidence interval for the entire sample cohort.

[`5_glimpse_concordance.sh`](5_glimpse_concordance.sh) - Calculates squared Pearson correlation between high-coverage and imputed dosages across chromosomes 1 to 22 and X, and writes them to several text files. Then runs [`lib/concordance_plot.py`](lib/concordance_plot.py) to create the concordance plots used in the publication.

[`6_pca.sh`](6_pca.sh) - Using `plink`, determines population-specific genetic markers from the 1000 Genomes reference panel, then extracts those variants from the imputed files, and finally runs [`lib/pca.py`](lib/pca.py) to plot the samples against the 1000 Genomes superpopulations.

[^1]: Sousa da Mota, B., Rubinacci, S., Cruz Dávalos, D.I. et al. Imputation of ancient human genomes. Nat Commun 14, 3660 (2023). https://doi.org/10.1038/s41467-023-39202-0
