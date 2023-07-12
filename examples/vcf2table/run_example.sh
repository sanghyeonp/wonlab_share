#!/bin/bash

script=/data1/sanghyeon/wonlab_contribute/combined/src/vcf2table/vcf2table.py


## VCF to table where input file is not compressed
python3 ../../src/vcf2table/vcf2table.py --file ./data/example.vcf \
                                        --delim-in tab \
                                        --outf vcf2table.example.1.txt


## VCF to table where input file is compressed by gzip
python3 ../../src/vcf2table/vcf2table.py --file ./data/example.vcf.gz \
                                        --delim-in tab \
                                        --compression-in gzip \
                                        --outf vcf2table.example.2.txt
