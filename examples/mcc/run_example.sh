#!/bin/bash

python3  ../../src/mcc/mcc.py --file data/inputf.tsv \
                    --p-col P P2 \
                    --delim-in tab \
                    --outf output.csv \
                    --delim-out comma \
                    --log \
                    --verbose