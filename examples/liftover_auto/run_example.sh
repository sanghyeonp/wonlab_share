#!/bin/bash
# script=/data1/sanghyeon/wonlab_contribute/combined/src/liftover_auto/liftover_auto.py
python3 ../../src/liftover_auto/liftover_auto.py --file data/gwas_catalog_v1.0_hg38.tsv \
                            --delim tab \
                            --snp-col SNP \
                            --chr-col Chr \
                            --pos-col Pos \
                            --build-from 38 \
                            --build-to 37 \
                            --keep-intermediate \
                            --unlifted-snplist \
                            --outf gwas_catalog_v1.0_hg37.tsv \
                            --verbose