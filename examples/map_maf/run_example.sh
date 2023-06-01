#!/bin/bash

# script=/data1/sanghyeon/wonlab_contribute/combined/src/map_maf/map_maf.py

python3 ../../src/map_maf/map_maf.py --file ./data/random_1k.txt \
                                    --delim-in whitespace \
                                    --chr-col CHR \
                                    --pos-col POS \
                                    --ancestry GLOBAL EAS EUR AMR AFR SAS \
                                    --outf maf.random_1k.csv \
                                    --delim-out comma
                                    