#!/bin/bash

# script=/data1/sanghyeon/wonlab_contribute/combined/src/map_af/map_af.py

python3 ../../src/map_af/map_af.py --file ./data/random_1k.txt \
                                    --delim-in whitespace \
                                    --chr-col CHR \
                                    --pos-col POS \
                                    --a1-col A1 \
                                    --ancestry GLOBAL EAS EUR AMR AFR SAS \
                                    --outf af.random_1k.csv \
                                    --delim-out comma
                                    