#!/bin/bash
script=/data1/sanghyeon/wonlab_contribute/combined/src/annovar/annovar_rsid_map.R

Rscript ${script} --gwas ./data/example.1k \
                --delim-in tab \
                --chr-col CHR \
                --pos-col POS \
                --ref-col Allele1 \
                --alt-col Allele2 \
                --genome-build 19 \
                --dbgap-build 150 \
                --nthread 5 \
                --save-mapping-file \
                --outf example.1k
