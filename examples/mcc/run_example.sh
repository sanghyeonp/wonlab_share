#!/bin/bash

python3  ../../src/mcc/mcc.py --file data/inputf.tsv \
                    --p_col P \
                    --delim tab \
                    --outf output.tsv \
                    --log \
                    --verbose