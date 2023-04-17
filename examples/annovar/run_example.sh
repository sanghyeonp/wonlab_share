#!/bin/bash

script=../../src/annovar/run_annovar.py

# python3 ${script} --file ./data/example.1k \
#                 --delim-in tab \
#                 --infer-chr-pos-ref-alt SNPID CHR:POS:REF:ALT : \
#                 --save-unannotated-snp \
#                 --save-flipped-snp \
#                 --save-mapping-file

python3 ${script} --file ./data/example.1k \
                --delim-in tab \
                --chr-col CHR \
                --pos-col POS \
                --ref-col Allele1 \
                --alt-col Allele2 \
                --save-unannotated-snp \
                --save-flipped-snp \
                --save-mapping-file \
                --delete-intermediate-files