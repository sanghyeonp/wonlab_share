#!/bin/bash

python3 ../liftover_auto.py --file data/gwas_catalog_v1.0_hg38.tsv \
                            --delim tab \
                            --snp_col SNP \
                            --chr_col Chr \
                            --pos_col Pos \
                            --build_from 38 \
                            --build_to 37 \
                            --keep_intermediate \
                            --unlifted_snplist \
                            --outf gwas_catalog_v1.0_hg37.tsv \
                            --verbose
