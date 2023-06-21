#!/bin/bash

# script=/data1/sanghyeon/wonlab_contribute/combined/src/compute_ld/compute_ld.py


# python3 ../../src/compute_ld/compute_ld.py \
#     --file-in ./data/example_data_50.txt \
#     --delim-in whitespace \
#     --threads 10 \
#     --reference 1kg \
#     --variant-identifier rsid \
#     --ancestry European \
#     --outf example_data_multi_pairwise.csv \
#     --delim-out comma


python3 ../../src/compute_ld/compute_ld.py \
    --file-in ./data/example_data_20.txt \
    --delim-in whitespace \
    --all-possible-pairs \
    --threads 15 \
    --reference 1kg \
    --variant-identifier rsid \
    --ancestry European \
    --outf example_data_multi_allpairwise.csv \
    --delim-out comma
