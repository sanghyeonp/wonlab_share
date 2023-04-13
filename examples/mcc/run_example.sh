#!/bin/bash

python3  ../../src/mcc/mcc.py --file data/inputf.tsv \
                    --p_col P P2 \
                    --delim_in tab \
                    --outf output.csv \
                    --delim_out comma \
                    --log \
                    --verbose