#!/bin/bash

# script=/data1/sanghyeon/wonlab_contribute/combined/src/map_snp_info/map_snp_info.py

script=../../src/map_snp_info/map_snp_info.py

python3 ${script} --file ./data/example_snps_1k.txt \
                --rsid-col variant_id \
                --build 37 \
                --outf abc.example_snps_1k.txt
