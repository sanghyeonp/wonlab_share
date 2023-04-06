#!/bin/bash

script=../../src/annovar/run_annovar.py

# python3 ${script} --file ./data/example.1k \
#                 --delim_in tab \
#                 --infer_chr_pos_ref_alt SNPID CHR:POS:REF:ALT : \
#                 --save_unannotated_snp \
#                 --save_flipped_snp \
#                 --save_mapping_file

python3 ${script} --file ./data/example.1k \
                --delim_in tab \
                --chr_col CHR \
                --pos_col POS \
                --ref_col Allele1 \
                --alt_col Allele2 \
                --save_unannotated_snp \
                --save_flipped_snp \
                --save_mapping_file \
                --delete_intermediate_files