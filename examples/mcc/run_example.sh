#!/bin/bash

python3  ../mcc.py --file data/inputf.tsv \
                    --p_col P \
                    --delim tab \
                    --outf output.tsv \
                    --log \
                    --verbose