#!/bin/bash
# script=/data1/sanghyeon/wonlab_contribute/combined/src/mcc/mcc.py

python3  ../../src/mcc/mcc.py --file data/inputf.tsv \
                    --p-col P P2 \
                    --delim-in tab \
                    --outf output.csv \
                    --delim-out comma \
                    --log \
                    --verbose