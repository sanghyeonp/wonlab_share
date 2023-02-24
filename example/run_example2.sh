#!/bin/bash

Rscript ../cohen_kappa.R --ref_gwas ./data/GWAS1.subset.csv \
                    --ref_gwas_delim comma \
                    --ref_gwas_snp SNP \
                    --ref_gwas_beta est \
                    --ref_name Trait1 \
                    --alt_gwas ./data/GWAS2.subset.tsv \
                    --alt_gwas_delim tab \
                    --alt_gwas_snp SNP \
                    --alt_gwas_beta Effect \
                    --alt_name Trait2 \
                    --snplist_FUMA ./data/leadSNPs.txt \
                    --outf example2 \
                    --verbose \
                    --rds \
                    --table
