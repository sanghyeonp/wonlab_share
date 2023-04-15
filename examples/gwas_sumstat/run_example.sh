#!/bin/bash

script=../../src/gwas_sumstat/gwas_sumstat_ukb_neale.R

# Options for --retain-col: [RSID, CHR, POS, ALT, REF, EAF, MAF, BETA, SE, PVAL, N]
# Can specify other columns in the summary statistics.

Rscript ${script} --gwas ./data/sumstat_ukbneale_100.tsv \
                --retain-col RSID CHR POS ALT REF EAF BETA SE PVAL N\
                --rename-col SNP chr pos A1 A2 Freq b se p n\
                --maf-filter 0.001
