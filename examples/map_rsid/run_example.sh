#!/bin/bash

script=../../src/map_rsid/map_rsid.py 

python3 ${script} --file data/example_sumstat_1k.tsv.gz \
                --build 37 \
                --file_compression gzip \
                --delim_in tab \
                --infer_chr_pos variant CHR:POS:REF:ALT : \
                --keep_refmap_position \
                --keep_refmap_allele \
                --outf example1.out.tsv \
                --delim_out tab

# Generate mapping file that can be used for other GWAS summary statistics of the same cohort
python3 ${script} --file data/example_sumstat_1k.tsv.gz \
                --build 37 \
                --file_compression gzip \
                --delim_in tab \
                --infer_chr_pos variant CHR:POS:REF:ALT : \
                --keep_refmap_position \
                --keep_refmap_allele \
                --outf example1.out.mapping.tsv \
                --delim_out tab \
                --keep_col variant SNP_ref