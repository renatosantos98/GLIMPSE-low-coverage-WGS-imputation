# Change to use a different number of CPU threads.
threads=16

# Results validation
GLIMPSE_concordance --thread $threads --input GLIMPSE_concordance/concordance.lst --minDP 8 --output GLIMPSE_concordance/output --minPROB 0.9999 --bins 0.00000 0.00100 0.00200 0.00500 0.01000 0.05000 0.10000 0.20000 0.50000

# Requires file concordance.lst with the following parameters
# 'Region' 'Allele frequencies at each site' 'The validation 30x dataset' 'The imputed data'
# chr22 1000G/1000G.chr22.sites.vcf.gz vcf/30x/manuel_corpas/chr22.vcf.gz GLIMPSE_ligated/manuel_corpas/chr22.merged.bcf
