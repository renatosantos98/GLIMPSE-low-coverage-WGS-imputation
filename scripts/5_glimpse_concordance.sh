#!/usr/bin/env bash
# Change to use a different number of CPU threads.
threads=16

# Results validation
# Example
# GLIMPSE_concordance --thread $threads --input GLIMPSE_concordance/concordance.lst --minDP 8 --output GLIMPSE_concordance/output --minPROB 0.9999 --bins 0.00000 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.10 0.15 0.20 0.3  0.4 0.50

# Requires file concordance.lst with the following parameters
# 'Region' 'Allele frequencies' 'The validation 30x vcf dataset' 'The imputed bcf data'
# chr22 1000G/1000G.chr22.sites.vcf.gz vcf/30x/IBS001/chr22.vcf.gz GLIMPSE_ligated/IBS001/chr22.merged.bcf

# Post-imputation validation
# BGI 1x vs BGI 30x
GLIMPSE_concordance --thread ${threads} --input GLIMPSE_concordance/IBS001/bgi1_bgi30/concordance.txt --minDP 8 --output GLIMPSE_concordance/IBS001/bgi1_bgi30/output --minPROB 0.9999 --bins 0.00000 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.10 0.15 0.20 0.3 0.4 0.50

# Illumina 1x vs Illumina 30x
GLIMPSE_concordance --thread ${threads} --input GLIMPSE_concordance/IBS001/illumina1_illumina30/concordance.txt --minDP 8 --output GLIMPSE_concordance/IBS001/illumina1_illumina30/output --minPROB 0.9999 --bins 0.00000 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.10 0.15 0.20 0.3 0.4 0.50

# BGI 1x vs Illumina 30x
GLIMPSE_concordance --thread ${threads} --input GLIMPSE_concordance/IBS001/bgi1_illumina30/concordance.txt --minDP 8 --output GLIMPSE_concordance/IBS001/bgi1_illumina30/output --minPROB 0.9999 --bins 0.00000 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.10 0.15 0.20 0.3 0.4 0.50

# Illumina 1x vs BGI 30x
GLIMPSE_concordance --thread ${threads} --input GLIMPSE_concordance/IBS001/illumina1_bgi30/concordance.txt --minDP 8 --output GLIMPSE_concordance/IBS001/illumina1_bgi30/output --minPROB 0.9999 --bins 0.00000 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.10 0.15 0.20 0.3 0.4 0.50

cd GLIMPSE_concordance/IBS001
python ../../scripts/7_concordance_plot.py
cd ../..

# Post-filtering validation
# BGI 1x vs BGI 30x
GLIMPSE_concordance --thread ${threads} --input GLIMPSE_concordance/IBS001_filtered/bgi1_bgi30/concordance.txt --minDP 8 --output GLIMPSE_concordance/IBS001_filtered/bgi1_bgi30/output --minPROB 0.9999 --bins 0.00000 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.10 0.15 0.20 0.3 0.4 0.50

# Illumina 1x vs Illumina 30x
GLIMPSE_concordance --thread ${threads} --input GLIMPSE_concordance/IBS001_filtered/illumina1_illumina30/concordance.txt --minDP 8 --output GLIMPSE_concordance/IBS001_filtered/illumina1_illumina30/output --minPROB 0.9999 --bins 0.00000 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.10 0.15 0.20 0.3 0.4 0.50

# BGI 1x vs Illumina 30x
GLIMPSE_concordance --thread ${threads} --input GLIMPSE_concordance/IBS001_filtered/bgi1_illumina30/concordance.txt --minDP 8 --output GLIMPSE_concordance/IBS001_filtered/bgi1_illumina30/output --minPROB 0.9999 --bins 0.00000 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.10 0.15 0.20 0.3 0.4 0.50

# Illumina 1x vs BGI 30x
GLIMPSE_concordance --thread ${threads} --input GLIMPSE_concordance/IBS001_filtered/illumina1_bgi30/concordance.txt --minDP 8 --output GLIMPSE_concordance/IBS001_filtered/illumina1_bgi30/output --minPROB 0.9999 --bins 0.00000 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.10 0.15 0.20 0.3 0.4 0.50

cd GLIMPSE_concordance/IBS001_filtered
python ../../scripts/lib/concordance_plot.py
cd ../..
