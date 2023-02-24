#!/bin/bash

python3  ../mcc.py --file data/inputf.tsv \
                    --pval P \
                    --delim tab \
                    --outf output.tsv \
                    --log \
                    --verbose